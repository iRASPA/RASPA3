module;

module spatial_decomposition_opencl_pair_kernel;

// OpenCL C (1.2) source of the device side of the cluster pair kernel: the list build, the pruning and the pair
// evaluation. The pair kernel is the device transcription of the lane loop of ClusterPairKernel<float>::compute
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
const char* const openclPairKernelSource = R"CLC(
#define CLUSTER_I 8
#define CLUSTER_J 4
#define GROUP_SIZE 32
#define PARTIALS 11
#define STENCIL 27
#define NO_ATOM 0xFFFFFFFFu

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
  uint padding[2];
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

inline float3 minimumImage(float3 dr, __constant const float* cell, __constant const float* inverseCell,
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
__kernel void clusterBounds(__global const float4* restrict buildPosition,  // x, y, z, molecule bits per slot
                            const uint numberOfJClusters,
                            __global float4* restrict clusterMin,
                            __global float4* restrict clusterMax)
{
  const uint J = get_global_id(0);
  if (J >= numberOfJClusters) return;
  float3 lo = (float3)(FLT_MAX, FLT_MAX, FLT_MAX);
  float3 hi = (float3)(-FLT_MAX, -FLT_MAX, -FLT_MAX);
  float real = 0.0f;
  for (uint b = 0; b < CLUSTER_J; ++b)
  {
    const float4 p = buildPosition[J * CLUSTER_J + b];
    if (as_uint(p.w) != NO_ATOM)
    {
      lo = fmin(lo, p.xyz);
      hi = fmax(hi, p.xyz);
      real = 1.0f;
    }
  }
  clusterMin[J] = (float4)(lo, real);
  clusterMax[J] = (float4)(hi, 0.0f);
}

// Writes the blocks found by the work-items of a round in order into the row of the work-group (the row keeps
// its true count; blocks beyond the capacity are dropped, the host grows the rows and rebuilds).
inline uint appendBlocks(__local uint* flags, uint id, uint J, uint mask, uint count, uint capacity,
                         __global uint* restrict rowCluster, __global uint* restrict rowMask)
{
  flags[id] = (mask != 0u) ? 1u : 0u;
  barrier(CLK_LOCAL_MEM_FENCE);
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
  barrier(CLK_LOCAL_MEM_FENCE);
  return count + all;
}

// The outer list of one i-cluster per work-group: the candidates are the j-clusters of the 27 cells around the
// cluster's cell (each neighbour cell once, also when the stencil wraps around a small axis); a work-item takes
// one candidate per round, tests its bounding box and builds the 32-bit mask of the pairs within the list cutoff
// between different molecules.
__kernel __attribute__((reqd_work_group_size(GROUP_SIZE, 1, 1)))
void buildList(__global const float4* restrict buildPosition,
               __global const uint* restrict cellSlotStart,     // first slot per coarse cell (cells + 1)
               __global const uint* restrict cellOfCluster,     // coarse cell per i-cluster
               __global const float4* restrict clusterMin,
               __global const float4* restrict clusterMax,
               __constant const BuildParameters* bp,
               __global uint* restrict outerCluster,            // rows of blocksPerCluster
               __global uint* restrict outerMask,
               __global uint* restrict outerCount)
{
  const uint I = get_group_id(0);
  const uint id = get_local_id(0);
  const uint capacity = bp->blocksPerCluster;
  __global uint* restrict rowCluster = outerCluster + I * capacity;
  __global uint* restrict rowMask = outerMask + I * capacity;

  __local float4 pi[CLUSTER_I];
  __local uint stencilCell[STENCIL];
  __local float4 stencilShift[STENCIL];  // translation of the neighbour cell; w != 0: image not fixed (minimum image)
  __local uint stencilFirst[STENCIL];    // first candidate index of the neighbour cell
  __local uint totalCandidates;
  __local uint flags[GROUP_SIZE];

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

  if (id < CLUSTER_I) pi[id] = buildPosition[I * CLUSTER_I + id];
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
  barrier(CLK_LOCAL_MEM_FENCE);
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
  barrier(CLK_LOCAL_MEM_FENCE);
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
        const float3 gap = fmax(fmax((float3)(0.0f, 0.0f, 0.0f), minI - maxJ), minJ - maxI);
        candidate = dot(gap, gap) < listCutoffSquared;
      }
      if (candidate)
      {
        for (uint b = 0; b < CLUSTER_J; ++b)
        {
          const float4 pj = buildPosition[J * CLUSTER_J + b];
          const uint moleculeJ = as_uint(pj.w);
          if (moleculeJ == NO_ATOM) continue;
          for (uint a = 0; a < CLUSTER_I; ++a)
          {
            const float4 pa = pi[a];
            const uint moleculeI = as_uint(pa.w);
            if (moleculeI == NO_ATOM || moleculeI == moleculeJ) continue;
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
__kernel __attribute__((reqd_work_group_size(GROUP_SIZE, 1, 1)))
void compactList(__global const float4* restrict position,
                 __global const uint* restrict outerCluster,
                 __global const uint* restrict outerMask,
                 __global const uint* restrict outerCount,
                 __constant const Parameters* p,
                 __global uint* restrict pairList,
                 __global uint* restrict laneCount)
{
  const uint I = get_group_id(0);
  const uint id = get_local_id(0);
  const uint a = id >> 2;
  const uint b = id & 3u;
  const uint capacity = p->blocksPerCluster;
  const uint laneCapacity = p->pairsPerLane;
  __global const uint* restrict rowCluster = outerCluster + I * capacity;
  __global const uint* restrict rowMask = outerMask + I * capacity;
  __global uint* restrict lane = pairList + (size_t)I * laneCapacity * GROUP_SIZE + id;
  const float3 pi = position[I * CLUSTER_I + a].xyz;
  const uint laneBit = 1u << id;
  const float innerCutoffSquared = p->innerCutoffSquared;

  const uint total = outerCount[I];
  uint count = 0u;
  for (uint q = 0; q < total; ++q)
  {
    const uint mask = rowMask[q];
    if ((mask & laneBit) == 0u) continue;
    const uint J = rowCluster[q];
    const float3 dr = minimumImage(pi - position[J * CLUSTER_J + b].xyz, p->cell, p->inverseCell, p->orthorhombic);
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
__kernel __attribute__((reqd_work_group_size(GROUP_SIZE, 1, 1)))
void clusterPairs(__global const float4* restrict position,      // x, y, z, q per slot (wrapped positions)
                  __global const uint* restrict typeOf,          // pseudo-atom type per slot
                  __global const uint* restrict pairList,        // lane lists [I][k][lane] of j-clusters
                  __global const uint* restrict laneCount,       // entries per lane
                  __global const float* restrict lennardJones,   // per type pair: 4 epsilon, sigma^6, shift
                  __constant const Parameters* p,
                  __global float4* restrict force,               // gradient per slot
                  __global float* restrict partials)             // per i-cluster: eVDW, eCharge, 9 strain terms
{
  const uint I = get_group_id(0);
  const uint id = get_local_id(0);
  const uint a = id >> 2;
  const uint b = id & 3u;
  const uint i = I * CLUSTER_I + a;
  const float4 pi = position[i];
  const uint rowOffset = typeOf[i] * p->numberOfTypes;
  const uint laneCapacity = p->pairsPerLane;
  __global const uint* restrict lane = pairList + (size_t)I * laneCapacity * GROUP_SIZE + id;

  __local uint counts[GROUP_SIZE];
  const uint n = min(laneCount[I * GROUP_SIZE + id], laneCapacity);
  counts[id] = n;
  barrier(CLK_LOCAL_MEM_FENCE);
  uint end = 0u;
  for (uint k = 0; k < GROUP_SIZE; ++k) end = max(end, counts[k]);

  const float cutOffVDWSquared = p->cutOffVDWSquared;
  const float cutOffChargeSquared = p->cutOffChargeSquared;
  const float alpha = p->alpha;
  const float alphaSquared = p->alphaSquared;
  const float alphaOverSqrtPi = p->alphaOverSqrtPi;
  const float coulombFactor = p->coulombFactor;
  const uint useCharge = p->useCharge;

  // Abramowitz-Stegun 7.1.26: erfc(x) = t (a1 + t (a2 + t (a3 + t (a4 + t a5)))) exp(-x^2), t = 1 / (1 + p x)
  const float asP = 0.3275911f;
  const float asA1 = 0.254829592f;
  const float asA2 = -0.284496736f;
  const float asA3 = 1.421413741f;
  const float asA4 = -1.453152027f;
  const float asA5 = 1.061405429f;
  const float minimumRR = 1e-4f;

  float3 fi = (float3)(0.0f, 0.0f, 0.0f);
  float eVDW = 0.0f;
  float eCharge = 0.0f;
  float sxx = 0.0f, syx = 0.0f, szx = 0.0f, sxy = 0.0f, syy = 0.0f, szy = 0.0f, sxz = 0.0f, syz = 0.0f, szz = 0.0f;

  for (uint k = 0; k < end; ++k)
  {
    const bool active = k < n;
    const uint J = active ? lane[k * GROUP_SIZE] : 0u;
    const uint j = J * CLUSTER_J + b;
    const float4 pj = position[j];
    const float inList = active ? 1.0f : 0.0f;

    float3 dr = pi.xyz - pj.xyz;
    dr = minimumImage(dr, p->cell, p->inverseCell, p->orthorhombic);
    const float rr = dot(dr, dr);
    const float rrSafe = fmax(rr, minimumRR);
    const float invRR = 1.0f / rrSafe;

    const uint t = 3u * (rowOffset + typeOf[j]);
    const float e4 = lennardJones[t];
    const float s6 = lennardJones[t + 1];
    const float sh = lennardJones[t + 2];
    const float maskVDW = (rr < cutOffVDWSquared) ? inList : 0.0f;
    const float invRR3 = invRR * invRR * invRR;
    const float rri3 = s6 * invRR3;
    const float rri6 = rri3 * rri3;
    eVDW += maskVDW * (e4 * (rri6 - rri3) - sh);
    float factor = maskVDW * (12.0f * e4 * rri3 * (0.5f - rri3) * invRR);

    if (useCharge)
    {
      const float qq = pi.w * pj.w;
      const float maskCharge = (rr < cutOffChargeSquared && qq != 0.0f) ? inList : 0.0f;
      const float r = sqrt(rrSafe);
      const float invR = r * invRR;
      const float tt = 1.0f / (1.0f + asP * alpha * r);
      const float gauss = exp(-(alphaSquared * rrSafe));
      const float erfcValue = tt * (asA1 + tt * (asA2 + tt * (asA3 + tt * (asA4 + tt * asA5)))) * gauss;
      const float u = erfcValue * invR;
      const float dudrr = -((alphaOverSqrtPi * gauss + 0.5f * erfcValue * invR) * invRR);
      const float prefactor = coulombFactor * (maskCharge * qq);
      eCharge += prefactor * u;
      factor += 2.0f * (prefactor * dudrr);
    }

    const float3 g = factor * dr;
    fi += g;
    sxx += g.x * dr.x;
    syx += g.y * dr.x;
    szx += g.z * dr.x;
    sxy += g.x * dr.y;
    syy += g.y * dr.y;
    szy += g.z * dr.y;
    sxz += g.x * dr.z;
    syz += g.y * dr.z;
    szz += g.z * dr.z;
  }

  // the force on atom i: the four j-lane partial sums of its row
  __local float3 rowForce[GROUP_SIZE];
  __local float acc[PARTIALS][GROUP_SIZE];
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
  barrier(CLK_LOCAL_MEM_FENCE);
  if (b == 0)
  {
    const float3 f = rowForce[id] + rowForce[id + 1] + rowForce[id + 2] + rowForce[id + 3];
    force[i] = (float4)(f.x, f.y, f.z, 0.0f);
  }
  if (id < PARTIALS)
  {
    float sum = 0.0f;
    for (uint n = 0; n < GROUP_SIZE; ++n) sum += acc[id][n];
    partials[I * PARTIALS + id] = sum;
  }
}
)CLC";
