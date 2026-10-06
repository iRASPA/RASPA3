module;

module spatial_decomposition_device_kernels;

// Device source of the cluster pair kernel: the list build, the pruning and the pair evaluation. Written in the
// kernel dialect of kernel_sources.ixx (the GLOBAL / LOCAL / KERNEL_GROUP_SIZE / FLOAT3 ... macros, defined by
// the backend that compiles the source; OpenCL C 1.2 today), so that the physics has one source per backend
// family. The pair kernel is the device transcription of the lane loop of ClusterPairKernel<float>::compute
// (cluster_kernel.ixx), so that the two stay comparable line by line.
//
// Slots: the atoms in coarse-cell order (cells at least the list cutoff wide), Morton-sorted within a cell and
// padded per cell to a multiple of 8, so that i-clusters (8 slots) and j-clusters (4 slots) are compact and never
// span two cells. A dummy slot carries the molecule NO_ATOM and a far-away build position.
//
// Lists: per i-cluster a strided row of at most `blocksPerCluster` blocks (j-cluster, 32-bit pair mask with bit
// a * 4 + b for the pair (8 I + a, 4 J + b)) and a count: the outer list (cutoff + Verlet skin), built from the
// 27-cell stencil of the i-cluster's cell with the positions of the binning. The pair kernel works on the lane
// lists derived from it at the current positions every few steps: per work-item (a, b) of the i-cluster the
// j-clusters J of its listed pairs within cutoff + prune skin, so that the kernel evaluates listed pairs only
// (the pairs of a block are about half real pairs: dummies, same-molecule pairs and pairs beyond the cutoff
// make up the rest). Both are *full* lists: every pair appears in the lists of both its clusters, so the pair
// kernel accumulates only the force on the i atom (no atomics, no reduction across work-groups) at twice the
// pair evaluations; the host halves the energies and the strain.
//
// Periodicity: the positions are the wrapped positions of the atoms (no ghost images); the pair vector is taken
// with the minimum image of the current cell (valid while cutoff + skin <= half the smallest perpendicular width,
// which the cell list checks). The list build uses the translation of the stencil offset, and the minimum image
// where the stencil wraps around a small axis.
//
// Same-molecule pairs: the list build leaves out the excluded and the scaled (1-4) pairs of a molecule
// (IntraMolecularExclusions::isExcludedFromPairList; the bonded kernel evaluates the scaled pairs), given as a
// CSR table over the atoms in system order (`exclusionStart` per atom, bit 31 set when every pair of the
// molecule is excluded, i.e. a rigid molecule; `exclusionPartner` the partner atoms) with the system atom of
// every slot (`slotAtom`). Every other same-molecule pair is listed and evaluated like a pair of two molecules.
const char* const deviceKernelPairSource = R"CLC(
#define CLUSTER_I 8
#define CLUSTER_J 4
#define GROUP_SIZE 32
#define PARTIALS 11
#define STENCIL 27
#define NO_ATOM 0xFFFFFFFFu
#define ALL_EXCLUDED_BIT 0x80000000u

typedef struct
{
  float cell[9];         // ax ay az bx by bz cx cy cz (columns of the cell matrix)
  float inverseCell[9];
  float cutOffVDWSquared;
  float cutOffChargeSquared;
  float alpha;
  float alphaSquared;
  float alphaOverSqrtPi;
  float coulombFactor;
  float innerCutoffSquared;  // (cutoff + prune skin)^2 of the inner list
  uint useCharge;
  uint orthorhombic;
  uint numberOfTypes;
  uint blocksPerCluster;     // stride of the outer list rows
  uint pairsPerLane;         // capacity of a lane list
  uint switchMode;           // Lennard-Jones switching on [r_s, rc]: 0 none, 1 potential switch, 2 force switch
  float switchDistanceSquared;
  float switchDistance;
  float switchInverseWidth;  // 1 / (rc - r_s)
  float switchInverseCutOff3;  // rc^-3
  float switchA12;           // rc^6 / (rc^6 - r_s^6)
  float switchA6;            // rc^3 / (rc^3 - r_s^3)
  uint padding[3];
} Parameters;

typedef struct
{
  float cell[9];
  float inverseCell[9];
  float listCutoffSquared;   // (cutoff + Verlet skin)^2 of the outer list
  uint gridX;
  uint gridY;
  uint gridZ;
  uint blocksPerCluster;
  uint orthorhombic;
} BuildParameters;

DEVICE_FUNCTION float3 minimumImage(float3 dr, CONSTANT const float* cell, CONSTANT const float* inverseCell,
                           uint orthorhombic)
{
  if (orthorhombic)
  {
    dr.x -= cell[0] * rint(dr.x * inverseCell[0]);
    dr.y -= cell[4] * rint(dr.y * inverseCell[4]);
    dr.z -= cell[8] * rint(dr.z * inverseCell[8]);
    return dr;
  }
  float3 s;
  s.x = inverseCell[0] * dr.x + inverseCell[3] * dr.y + inverseCell[6] * dr.z;
  s.y = inverseCell[1] * dr.x + inverseCell[4] * dr.y + inverseCell[7] * dr.z;
  s.z = inverseCell[2] * dr.x + inverseCell[5] * dr.y + inverseCell[8] * dr.z;
  s -= rint(s);
  float3 r;
  r.x = cell[0] * s.x + cell[3] * s.y + cell[6] * s.z;
  r.y = cell[1] * s.x + cell[4] * s.y + cell[7] * s.z;
  r.z = cell[2] * s.x + cell[5] * s.y + cell[8] * s.z;
  return r;
}

// Bounding box of every j-cluster over its real atoms; w of the minimum is 1 when the cluster has real atoms.
KERNEL void clusterBounds(GLOBAL const float4* RESTRICT buildPosition,  // x, y, z, molecule bits per slot
                            VALUE_ARG(uint, numberOfJClusters),
                            GLOBAL float4* RESTRICT clusterMin,
                            GLOBAL float4* RESTRICT clusterMax KERNEL_INDEX_ARGS)
{
  const uint J = GLOBAL_ID();
  if (J >= numberOfJClusters) return;
  float3 lo = FLOAT3(FLT_MAX, FLT_MAX, FLT_MAX);
  float3 hi = FLOAT3(-FLT_MAX, -FLT_MAX, -FLT_MAX);
  float real = 0.0f;
  for (uint b = 0; b < CLUSTER_J; ++b)
  {
    const float4 p = buildPosition[J * CLUSTER_J + b];
    if (AS_UINT(p.w) != NO_ATOM)
    {
      lo = fmin(lo, p.xyz);
      hi = fmax(hi, p.xyz);
      real = 1.0f;
    }
  }
  clusterMin[J] = FLOAT4(lo, real);
  clusterMax[J] = FLOAT4(hi, 0.0f);
}

// Writes the blocks found by the work-items of a round in order into the row of the work-group (the row keeps
// its true count; blocks beyond the capacity are dropped, the host grows the rows and rebuilds).
DEVICE_FUNCTION uint appendBlocks(LOCAL uint* flags, uint id, uint J, uint mask, uint count, uint capacity,
                         GLOBAL uint* RESTRICT rowCluster, GLOBAL uint* RESTRICT rowMask)
{
  flags[id] = (mask != 0u) ? 1u : 0u;
  LOCAL_BARRIER();
  uint before = 0u;
  uint all = 0u;
  for (uint k = 0; k < GROUP_SIZE; ++k)
  {
    const uint f = flags[k];
    all += f;
    before += (k < id) ? f : 0u;
  }
  if (mask != 0u && count + before < capacity)
  {
    rowCluster[count + before] = J;
    rowMask[count + before] = mask;
  }
  LOCAL_BARRIER();
  return count + all;
}

// Whether the same-molecule pair (atomI, atomJ) (system atom indices) is left out of the lists: excluded or
// scaled. A linear scan of the short partner list of atomI (1-2, 1-3 and 1-4 partners; a rigid molecule is
// flagged instead of listed).
DEVICE_FUNCTION bool excludedFromList(uint atomI, uint atomJ, GLOBAL const uint* RESTRICT exclusionStart,
                                      GLOBAL const uint* RESTRICT exclusionPartner)
{
  const uint start = exclusionStart[atomI];
  if (start & ALL_EXCLUDED_BIT) return true;
  const uint end = exclusionStart[atomI + 1] & ~ALL_EXCLUDED_BIT;
  for (uint k = start; k < end; ++k)
  {
    if (exclusionPartner[k] == atomJ) return true;
  }
  return false;
}

// The outer list of one i-cluster per work-group: the candidates are the j-clusters of the 27 cells around the
// cluster's cell (each neighbour cell once, also when the stencil wraps around a small axis); a work-item takes
// one candidate per round, tests its bounding box and builds the 32-bit mask of the pairs within the list cutoff
// (all pairs of two molecules, and the same-molecule pairs that are neither excluded nor scaled).
KERNEL_GROUP_SIZE(GROUP_SIZE)
void buildList(GLOBAL const float4* RESTRICT buildPosition,
               GLOBAL const uint* RESTRICT slotAtom,          // system atom per slot (NO_ATOM for a dummy slot)
               GLOBAL const uint* RESTRICT exclusionStart,    // per system atom: CSR start of its left-out partners
               GLOBAL const uint* RESTRICT exclusionPartner,  // the left-out partners (system atoms)
               GLOBAL const uint* RESTRICT cellSlotStart,     // first slot per coarse cell (cells + 1)
               GLOBAL const uint* RESTRICT cellOfCluster,     // coarse cell per i-cluster
               GLOBAL const float4* RESTRICT clusterMin,
               GLOBAL const float4* RESTRICT clusterMax,
               CONSTANT const BuildParameters* bp,
               GLOBAL uint* RESTRICT outerCluster,            // rows of blocksPerCluster
               GLOBAL uint* RESTRICT outerMask,
               GLOBAL uint* RESTRICT outerCount KERNEL_INDEX_ARGS)
{
  const uint I = GROUP_ID();
  const uint id = LOCAL_ID();
  const uint capacity = bp->blocksPerCluster;
  GLOBAL uint* RESTRICT rowCluster = outerCluster + I * capacity;
  GLOBAL uint* RESTRICT rowMask = outerMask + I * capacity;

  LOCAL_DECL(float4 pi[CLUSTER_I]);
  LOCAL_DECL(uint atomI[CLUSTER_I]);
  LOCAL_DECL(uint stencilCell[STENCIL]);
  LOCAL_DECL(float4 stencilShift[STENCIL]);  // translation of the neighbour cell; w != 0: image not fixed (minimum image)
  LOCAL_DECL(uint stencilFirst[STENCIL]);    // first candidate index of the neighbour cell
  LOCAL_DECL(uint totalCandidates);
  LOCAL_DECL(uint flags[GROUP_SIZE]);

  const float4 min0 = clusterMin[2 * I];
  const float4 min1 = clusterMin[2 * I + 1];
  const float4 max0 = clusterMax[2 * I];
  const float4 max1 = clusterMax[2 * I + 1];
  if (min0.w == 0.0f && min1.w == 0.0f)
  {
    if (id == 0) outerCount[I] = 0u;
    return;
  }
  const float3 minI = fmin(min0.xyz, min1.xyz);
  const float3 maxI = fmax(max0.xyz, max1.xyz);
  const float listCutoffSquared = bp->listCutoffSquared;

  if (id < CLUSTER_I)
  {
    pi[id] = buildPosition[I * CLUSTER_I + id];
    atomI[id] = slotAtom[I * CLUSTER_I + id];
  }
  if (id < STENCIL)
  {
    const uint c = cellOfCluster[I];
    const int gx = (int)bp->gridX;
    const int gy = (int)bp->gridY;
    const int gz = (int)bp->gridZ;
    const int cx = (int)(c % bp->gridX);
    const int cy = (int)((c / bp->gridX) % bp->gridY);
    const int cz = (int)(c / (bp->gridX * bp->gridY));
    const int dx = (int)(id % 3u) - 1;
    const int dy = (int)((id / 3u) % 3u) - 1;
    const int dz = (int)(id / 9u) - 1;
    // on an axis with fewer than 3 cells several offsets reach the same cell: keep the first, and take the
    // minimum image for the pairs of such a cell
    const bool first = (dx - gx < -1) && (dy - gy < -1) && (dz - gz < -1);
    const bool ambiguous = (dx - gx >= -1) || (dx + gx <= 1) || (dy - gy >= -1) || (dy + gy <= 1) ||
                           (dz - gz >= -1) || (dz + gz <= 1);
    int nx = cx + dx, ny = cy + dy, nz = cz + dz;
    int wx = 0, wy = 0, wz = 0;
    if (nx < 0) { nx += gx; wx = -1; } else if (nx >= gx) { nx -= gx; wx = 1; }
    if (ny < 0) { ny += gy; wy = -1; } else if (ny >= gy) { ny -= gy; wy = 1; }
    if (nz < 0) { nz += gz; wz = -1; } else if (nz >= gz) { nz -= gz; wz = 1; }
    const uint n = (uint)((nz * gy + ny) * gx + nx);
    const float fx = (float)wx, fy = (float)wy, fz = (float)wz;
    float4 shift;
    shift.x = bp->cell[0] * fx + bp->cell[3] * fy + bp->cell[6] * fz;
    shift.y = bp->cell[1] * fx + bp->cell[4] * fy + bp->cell[7] * fz;
    shift.z = bp->cell[2] * fx + bp->cell[5] * fy + bp->cell[8] * fz;
    shift.w = ambiguous ? 1.0f : 0.0f;
    stencilCell[id] = n;
    stencilShift[id] = shift;
    stencilFirst[id] = first ? (cellSlotStart[n + 1] - cellSlotStart[n]) / CLUSTER_J : 0u;
  }
  LOCAL_BARRIER();
  if (id == 0)
  {
    uint sum = 0u;
    for (uint s = 0; s < STENCIL; ++s)
    {
      const uint n = stencilFirst[s];
      stencilFirst[s] = sum;
      sum += n;
    }
    totalCandidates = sum;
  }
  LOCAL_BARRIER();
  const uint total = totalCandidates;

  uint count = 0u;
  for (uint q0 = 0; q0 < total; q0 += GROUP_SIZE)
  {
    const uint q = q0 + id;
    uint J = 0u;
    uint mask = 0u;
    if (q < total)
    {
      uint s = 0u;  // the stencil entry owning candidate q: the last one starting at or before q
      for (uint k = 1; k < STENCIL; ++k) s = (stencilFirst[k] <= q) ? k : s;
      J = cellSlotStart[stencilCell[s]] / CLUSTER_J + (q - stencilFirst[s]);
      const float4 shift = stencilShift[s];
      const bool ambiguous = shift.w != 0.0f;
      const float4 lo = clusterMin[J];
      bool candidate = lo.w != 0.0f;
      if (candidate && !ambiguous)
      {
        // distance between the bounding boxes, with the translation applied to the j box
        const float3 minJ = lo.xyz + shift.xyz;
        const float3 maxJ = clusterMax[J].xyz + shift.xyz;
        const float3 gap = fmax(fmax(FLOAT3(0.0f, 0.0f, 0.0f), minI - maxJ), minJ - maxI);
        candidate = dot(gap, gap) < listCutoffSquared;
      }
      if (candidate)
      {
        for (uint b = 0; b < CLUSTER_J; ++b)
        {
          const float4 pj = buildPosition[J * CLUSTER_J + b];
          const uint moleculeJ = AS_UINT(pj.w);
          if (moleculeJ == NO_ATOM) continue;
          const uint atomJ = slotAtom[J * CLUSTER_J + b];
          for (uint a = 0; a < CLUSTER_I; ++a)
          {
            const float4 pa = pi[a];
            const uint moleculeI = AS_UINT(pa.w);
            if (moleculeI == NO_ATOM) continue;
            if (moleculeI == moleculeJ)
            {
              // the same atom (the clusters of one cell overlap) or a left-out same-molecule pair
              const uint ai = atomI[a];
              if (ai == atomJ || excludedFromList(ai, atomJ, exclusionStart, exclusionPartner)) continue;
            }
            float3 dr = pa.xyz - pj.xyz - shift.xyz;
            if (ambiguous) dr = minimumImage(dr, bp->cell, bp->inverseCell, bp->orthorhombic);
            if (dot(dr, dr) < listCutoffSquared) mask |= 1u << (a * CLUSTER_J + b);
          }
        }
      }
    }
    count = appendBlocks(flags, id, J, mask, count, capacity, rowCluster, rowMask);
  }
  if (id == 0) outerCount[I] = count;
}

// The pair list of one i-cluster per work-group from its outer row at the current positions: work-item (a, b),
// a = id / 4, b = id % 4, collects the j-clusters J of the listed pairs (8 I + a, 4 J + b) within the inner cutoff
// into its own lane list (layout [I][k][lane], so that the 32 lanes read consecutive words). The count is exact;
// entries beyond the lane capacity are dropped and the host grows the lists and compacts again.
KERNEL_GROUP_SIZE(GROUP_SIZE)
void compactList(GLOBAL const float4* RESTRICT position,
                 GLOBAL const uint* RESTRICT outerCluster,
                 GLOBAL const uint* RESTRICT outerMask,
                 GLOBAL const uint* RESTRICT outerCount,
                 CONSTANT const Parameters* p,
                 GLOBAL uint* RESTRICT pairList,
                 GLOBAL uint* RESTRICT laneCount KERNEL_INDEX_ARGS)
{
  const uint I = GROUP_ID();
  const uint id = LOCAL_ID();
  const uint a = id >> 2;
  const uint b = id & 3u;
  const uint capacity = p->blocksPerCluster;
  const uint laneCapacity = p->pairsPerLane;
  GLOBAL const uint* RESTRICT rowCluster = outerCluster + I * capacity;
  GLOBAL const uint* RESTRICT rowMask = outerMask + I * capacity;
  GLOBAL uint* RESTRICT lane = pairList + (size_t)I * laneCapacity * GROUP_SIZE + id;
  const float3 pi = position[I * CLUSTER_I + a].xyz;
  const uint laneBit = 1u << id;
  const float innerCutoffSquared = p->innerCutoffSquared;
  const uint orthorhombic = p->orthorhombic;
  const float ax = p->cell[0], ay = p->cell[4], az = p->cell[8];
  const float iax = p->inverseCell[0], iay = p->inverseCell[4], iaz = p->inverseCell[8];

  const uint total = outerCount[I];
  uint count = 0u;
  for (uint q = 0; q < total; ++q)
  {
    const uint mask = rowMask[q];
    if ((mask & laneBit) == 0u) continue;
    const uint J = rowCluster[q];
    float3 dr = pi - position[J * CLUSTER_J + b].xyz;
    if (orthorhombic)
    {
      dr.x -= ax * rint(dr.x * iax);
      dr.y -= ay * rint(dr.y * iay);
      dr.z -= az * rint(dr.z * iaz);
    }
    else
    {
      dr = minimumImage(dr, p->cell, p->inverseCell, 0u);
    }
    if (dot(dr, dr) < innerCutoffSquared)
    {
      if (count < laneCapacity) lane[count * GROUP_SIZE] = J;
      ++count;
    }
  }
  laneCount[I * GROUP_SIZE + id] = count;
}

// One work-group of 32 work-items per i-cluster; work-item (a, b) = (id / 4, id % 4) evaluates the pairs
// (8 I + a, 4 J + b) of its lane list. Every trip of the loop evaluates a listed pair (the lanes idle only beyond
// their own count, up to the longest lane of the group). Energies and the strain derivative sum
// g_ij (x) (x_i - x_j) are accumulated per work-item and reduced per work-group.
KERNEL_GROUP_SIZE(GROUP_SIZE)
void clusterPairs(GLOBAL const float4* RESTRICT position,      // x, y, z, q per slot (wrapped positions)
                  GLOBAL const uint* RESTRICT typeOf,          // pseudo-atom type per slot
                  GLOBAL const uint* RESTRICT pairList,        // lane lists [I][k][lane] of j-clusters
                  GLOBAL const uint* RESTRICT laneCount,       // entries per lane
                  GLOBAL const float4* RESTRICT lennardJones,  // per type pair: 4 epsilon, sigma^6, shift, 0
                  CONSTANT const Parameters* p,
                  GLOBAL float4* RESTRICT force,               // gradient per slot
                  GLOBAL float* RESTRICT partials              // per i-cluster: eVDW, eCharge, 9 strain terms
                  KERNEL_INDEX_ARGS)
{
  const uint I = GROUP_ID();
  const uint id = LOCAL_ID();
  const uint a = id >> 2;
  const uint b = id & 3u;
  const uint i = I * CLUSTER_I + a;
  const float4 pi = position[i];
  const uint rowOffset = typeOf[i] * p->numberOfTypes;
  const uint laneCapacity = p->pairsPerLane;
  GLOBAL const uint* RESTRICT lane = pairList + (size_t)I * laneCapacity * GROUP_SIZE + id;

  LOCAL_DECL(uint counts[GROUP_SIZE]);
  const uint n = min(laneCount[I * GROUP_SIZE + id], laneCapacity);
  counts[id] = n;
  LOCAL_BARRIER();
  uint end = 0u;
  for (uint k = 0; k < GROUP_SIZE; ++k) end = max(end, counts[k]);

  const float cutOffVDWSquared = p->cutOffVDWSquared;
  const float cutOffChargeSquared = p->cutOffChargeSquared;
  const float alpha = p->alpha;
  const float alphaSquared = p->alphaSquared;
  const float alphaOverSqrtPi = p->alphaOverSqrtPi;
  const float coulombFactor = p->coulombFactor;
  const uint useCharge = p->useCharge;
  const uint orthorhombic = p->orthorhombic;
  const uint switchMode = p->switchMode;
  const float switchDistanceSquared = p->switchDistanceSquared;
  const float switchDistance = p->switchDistance;
  const float switchInverseWidth = p->switchInverseWidth;
  const float switchInverseCutOff3 = p->switchInverseCutOff3;
  const float switchA12 = p->switchA12;
  const float switchA6 = p->switchA6;
  const float ax = p->cell[0], ay = p->cell[4], az = p->cell[8];
  const float iax = p->inverseCell[0], iay = p->inverseCell[4], iaz = p->inverseCell[8];

  // Abramowitz-Stegun 7.1.26: erfc(x) = t (a1 + t (a2 + t (a3 + t (a4 + t a5)))) exp(-x^2), t = 1 / (1 + p x)
  const float asP = 0.3275911f;
  const float asA1 = 0.254829592f;
  const float asA2 = -0.284496736f;
  const float asA3 = 1.421413741f;
  const float asA4 = -1.453152027f;
  const float asA5 = 1.061405429f;
  const float minimumRR = 1e-4f;

  float3 fi = FLOAT3(0.0f, 0.0f, 0.0f);
  float eVDW = 0.0f;
  float eCharge = 0.0f;
  float sxx = 0.0f, syx = 0.0f, szx = 0.0f, sxy = 0.0f, syy = 0.0f, szy = 0.0f, sxz = 0.0f, syz = 0.0f, szz = 0.0f;

  // one listed pair: the Lennard-Jones and real-space Ewald terms with the i-side gradient and the strain
  // derivative accumulated; `inList` is 0 when this lane is past its own count (padded to the longest lane)
#define PAIR_TERM(pj, lj, inList)                                                                             \
  {                                                                                                           \
    float3 dr = pi.xyz - (pj).xyz;                                                                            \
    if (orthorhombic)                                                                                         \
    {                                                                                                         \
      dr.x -= ax * rint(dr.x * iax);                                                                          \
      dr.y -= ay * rint(dr.y * iay);                                                                          \
      dr.z -= az * rint(dr.z * iaz);                                                                          \
    }                                                                                                         \
    else                                                                                                      \
    {                                                                                                         \
      dr = minimumImage(dr, p->cell, p->inverseCell, 0u);                                                     \
    }                                                                                                         \
    const float rr = dot(dr, dr);                                                                             \
    const float rrSafe = fmax(rr, minimumRR);                                                                 \
    const float invRR = 1.0f / rrSafe;                                                                        \
    const float maskVDW = (rr < cutOffVDWSquared) ? (inList) : 0.0f;                                          \
    const float invRR3 = invRR * invRR * invRR;                                                               \
    const float rri3 = (lj).y * invRR3;                                                                       \
    const float rri6 = rri3 * rri3;                                                                           \
    float uVDW = (lj).x * (rri6 - rri3) - (lj).z;                                                             \
    float fVDW = 12.0f * (lj).x * rri3 * (0.5f - rri3) * invRR;                                               \
    if (switchMode != 0u && rr > switchDistanceSquared)                                                       \
    {                                                                                                         \
      const float rSw = sqrt(rrSafe);                                                                         \
      const float invRSw = rSw * invRR;                                                                       \
      if (switchMode == 1u)                                                                                   \
      {                                                                                                       \
        const float x = (rSw - switchDistance) * switchInverseWidth;                                          \
        const float x2 = x * x;                                                                               \
        const float oneMinusX = 1.0f - x;                                                                     \
        const float sw = 1.0f + x2 * x * (-10.0f + x * (15.0f - 6.0f * x));                                   \
        const float dsw = -30.0f * switchInverseWidth * x2 * oneMinusX * oneMinusX;                           \
        const float uPlain = (lj).x * (rri6 - rri3);                                                          \
        fVDW = fVDW * sw + uPlain * dsw * invRSw;                                                             \
        uVDW = uPlain * sw;                                                                                   \
      }                                                                                                       \
      else                                                                                                    \
      {                                                                                                       \
        const float c6 = (lj).x * (lj).y;                                                                     \
        const float c12 = c6 * (lj).y;                                                                        \
        const float invR3 = invRSw * invRR;                                                                   \
        const float invR4 = invR3 * invRSw;                                                                   \
        const float d3 = invR3 - switchInverseCutOff3;                                                        \
        const float d6 = invR3 * invR3 - switchInverseCutOff3 * switchInverseCutOff3;                         \
        uVDW = switchA12 * c12 * d6 * d6 - switchA6 * c6 * d3 * d3;                                           \
        fVDW = (-12.0f * switchA12 * c12 * d6 * invR3 * invR4 + 6.0f * switchA6 * c6 * d3 * invR4) * invRSw;  \
      }                                                                                                       \
    }                                                                                                         \
    eVDW += maskVDW * uVDW;                                                                                   \
    float factor = maskVDW * fVDW;                                                                            \
    if (useCharge)                                                                                            \
    {                                                                                                         \
      const float qq = pi.w * (pj).w;                                                                         \
      const float maskCharge = (rr < cutOffChargeSquared && qq != 0.0f) ? (inList) : 0.0f;                    \
      const float r = sqrt(rrSafe);                                                                           \
      const float invR = r * invRR;                                                                           \
      const float tt = 1.0f / (1.0f + asP * alpha * r);                                                       \
      const float gauss = exp(-(alphaSquared * rrSafe));                                                      \
      const float erfcValue = tt * (asA1 + tt * (asA2 + tt * (asA3 + tt * (asA4 + tt * asA5)))) * gauss;      \
      const float u = erfcValue * invR;                                                                       \
      const float dudrr = -((alphaOverSqrtPi * gauss + 0.5f * erfcValue * invR) * invRR);                     \
      const float prefactor = coulombFactor * (maskCharge * qq);                                              \
      eCharge += prefactor * u;                                                                               \
      factor += 2.0f * (prefactor * dudrr);                                                                   \
    }                                                                                                         \
    const float3 g = factor * dr;                                                                             \
    fi += g;                                                                                                  \
    sxx += g.x * dr.x;                                                                                        \
    syx += g.y * dr.x;                                                                                        \
    szx += g.z * dr.x;                                                                                        \
    sxy += g.x * dr.y;                                                                                        \
    syy += g.y * dr.y;                                                                                        \
    szy += g.z * dr.y;                                                                                        \
    sxz += g.x * dr.z;                                                                                        \
    syz += g.y * dr.z;                                                                                        \
    szz += g.z * dr.z;                                                                                        \
  }

  // two listed pairs per iteration: the j-slot loads of the second are in flight while the first is evaluated
  const uint endPairs = end & ~1u;
  for (uint k = 0; k < endPairs; k += 2)
  {
    const bool active0 = k < n;
    const bool active1 = (k + 1u) < n;
    const uint J0 = active0 ? lane[k * GROUP_SIZE] : 0u;
    const uint J1 = active1 ? lane[(k + 1u) * GROUP_SIZE] : 0u;
    const uint j0 = J0 * CLUSTER_J + b;
    const uint j1 = J1 * CLUSTER_J + b;
    const float4 pj0 = position[j0];
    const float4 pj1 = position[j1];
    const float4 lj0 = lennardJones[rowOffset + typeOf[j0]];
    const float4 lj1 = lennardJones[rowOffset + typeOf[j1]];
    PAIR_TERM(pj0, lj0, active0 ? 1.0f : 0.0f)
    PAIR_TERM(pj1, lj1, active1 ? 1.0f : 0.0f)
  }
  if (endPairs < end)
  {
    const bool active = endPairs < n;
    const uint J = active ? lane[endPairs * GROUP_SIZE] : 0u;
    const uint j = J * CLUSTER_J + b;
    const float4 pj = position[j];
    const float4 lj = lennardJones[rowOffset + typeOf[j]];
    PAIR_TERM(pj, lj, active ? 1.0f : 0.0f)
  }
#undef PAIR_TERM

  // the force on atom i: the four j-lane partial sums of its row
  LOCAL_DECL(float3 rowForce[GROUP_SIZE]);
  LOCAL_DECL(float acc[PARTIALS][GROUP_SIZE]);
  rowForce[id] = fi;
  acc[0][id] = eVDW;
  acc[1][id] = eCharge;
  acc[2][id] = sxx;
  acc[3][id] = syx;
  acc[4][id] = szx;
  acc[5][id] = sxy;
  acc[6][id] = syy;
  acc[7][id] = szy;
  acc[8][id] = sxz;
  acc[9][id] = syz;
  acc[10][id] = szz;
  LOCAL_BARRIER();
  if (b == 0)
  {
    const float3 f = rowForce[id] + rowForce[id + 1] + rowForce[id + 2] + rowForce[id + 3];
    force[i] = FLOAT4(f.x, f.y, f.z, 0.0f);
  }
  if (id < PARTIALS)
  {
    float sum = 0.0f;
    for (uint n = 0; n < GROUP_SIZE; ++n) sum += acc[id][n];
    partials[I * PARTIALS + id] = sum;
  }
}
)CLC";
