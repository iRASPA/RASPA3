module;

module cbmc_device_kernels;

// Device source of the multiple-try CBMC grower: the external (framework + inter-molecular) energy of every trial
// molecule of a batch against the static background of the move, in single precision.
//
// One work-group per trial set (a complete rigid molecule of `atomsPerSet` atoms): the lanes stride over the
// background atoms, each lane accumulating the Lennard-Jones and real-space Ewald terms of its background atoms
// with all atoms of the set, and the group reduces the lane sums in a fixed tree. The result is therefore a
// deterministic function of the inputs (no atomics, no launch-order dependence), which the proposal scheme
// relies on: the energy of a configuration must be the same whether it is evaluated as a fresh trial or as the
// retraced existing molecule.
//
// Exclusions follow CBMC::computeInterMolecularEnergy: a background atom whose molecule id equals the trial
// molecule's id, or the move's skipped molecule, is left out. Framework atoms carry the framework bit in their
// type word and use the framework cutoff; the other atoms the molecule cutoff. An overlap (a Lennard-Jones pair
// energy above the overlap criterion, as on the host) is counted in a separate partial, so the host sees it
// without comparing against a sentinel energy.
const char* const cbmcDeviceGrowKernelSource = R"CLC(
#define GROW_GROUP 64
#define GROW_MAX_ATOMS 32
#define GROW_PARTIALS 3
#define FRAMEWORK_BIT 0x80000000u
#define TYPE_MASK 0x7FFFFFFFu

typedef struct
{
  float cell[9];         // ax ay az bx by bz cx cy cz (columns of the cell matrix)
  float inverseCell[9];
  float cutOffFrameworkVDWSquared;
  float cutOffMoleculeVDWSquared;
  float cutOffChargeSquared;
  float alpha;
  float coulombFactor;
  float overlapCriteria;
  uint useCharge;
  uint orthorhombic;
  uint numberOfTypes;
  uint numberOfBackground;
  uint atomsPerSet;
  uint trialMoleculeId;
  uint skipMoleculeId;
  uint padding[5];
} GrowParameters;

DEVICE_FUNCTION float3 growMinimumImage(float3 dr, CONSTANT const float* cell, CONSTANT const float* inverseCell,
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

// Partials per trial set: Lennard-Jones energy, real-space Ewald energy, number of overlapping pairs.
KERNEL_GROUP_SIZE(GROW_GROUP)
void growEnergies(GLOBAL const float4* RESTRICT background,          // per background atom: position, charge
                  GLOBAL const uint* RESTRICT backgroundType,         // per background atom: type | framework bit
                  GLOBAL const uint* RESTRICT backgroundMolecule,     // per background atom: molecule id
                  GLOBAL const float4* RESTRICT trial,                // per trial atom (set-major): position, charge
                  GLOBAL const uint* RESTRICT trialType,              // per trial atom: type
                  GLOBAL const float4* RESTRICT lennardJones,         // per type pair: 4 eps, sigma^2, shift, 0
                  CONSTANT const GrowParameters* p,
                  GLOBAL float* RESTRICT partials KERNEL_INDEX_ARGS)
{
  LOCAL_DECL(float4 setPosition[GROW_MAX_ATOMS]);
  LOCAL_DECL(uint setType[GROW_MAX_ATOMS]);
  LOCAL_DECL(float scratch[GROW_PARTIALS * GROW_GROUP]);

  const uint lid = LOCAL_ID();
  const uint set = GROUP_ID();
  const uint n = p->atomsPerSet;
  if (lid < n)
  {
    setPosition[lid] = trial[set * n + lid];
    setType[lid] = trialType[set * n + lid];
  }
  LOCAL_BARRIER();

  const uint useCharge = p->useCharge;
  const uint orthorhombic = p->orthorhombic;
  const uint numberOfTypes = p->numberOfTypes;
  const float cutOffChargeSquared = p->cutOffChargeSquared;
  const float alpha = p->alpha;
  const float coulombFactor = p->coulombFactor;
  const float overlapCriteria = p->overlapCriteria;
  const float minimumRR = 1e-8f;

  float eVDW = 0.0f;
  float eCharge = 0.0f;
  float overlaps = 0.0f;

  for (uint j = lid; j < p->numberOfBackground; j += GROW_GROUP)
  {
    const uint molecule = backgroundMolecule[j];
    if (molecule == p->trialMoleculeId || molecule == p->skipMoleculeId) continue;
    const float4 pj = background[j];
    const uint typeWord = backgroundType[j];
    const uint typeJ = typeWord & TYPE_MASK;
    const float cutOffVDWSquared =
        (typeWord & FRAMEWORK_BIT) ? p->cutOffFrameworkVDWSquared : p->cutOffMoleculeVDWSquared;

    for (uint i = 0; i < n; ++i)
    {
      const float4 pi = setPosition[i];
      float3 dr = pi.xyz - pj.xyz;
      dr = growMinimumImage(dr, p->cell, p->inverseCell, orthorhombic);
      const float rr = dot(dr, dr);
      const float rrSafe = fmax(rr, minimumRR);
      if (rr < cutOffVDWSquared)
      {
        const float4 lj = lennardJones[setType[i] * numberOfTypes + typeJ];
        const float sigmaOverRR = lj.y / rrSafe;  // (sigma/r)^2
        const float rri3 = sigmaOverRR * sigmaOverRR * sigmaOverRR;
        const float rri6 = rri3 * rri3;
        const float u = lj.x * (rri6 - rri3) - lj.z;
        if (u > overlapCriteria)
        {
          overlaps += 1.0f;
        }
        else
        {
          eVDW += u;
        }
      }
      if (useCharge && rr < cutOffChargeSquared)
      {
        const float qq = pi.w * pj.w;
        if (qq != 0.0f)
        {
          const float r = sqrt(rrSafe);
          eCharge += coulombFactor * qq * erfc(alpha * r) / r;
        }
      }
    }
  }

  scratch[lid] = eVDW;
  scratch[GROW_GROUP + lid] = eCharge;
  scratch[2 * GROW_GROUP + lid] = overlaps;
  for (uint stride = GROW_GROUP / 2; stride > 0; stride >>= 1)
  {
    LOCAL_BARRIER();
    if (lid < stride)
    {
      for (uint q = 0; q < GROW_PARTIALS; ++q)
      {
        scratch[q * GROW_GROUP + lid] += scratch[q * GROW_GROUP + lid + stride];
      }
    }
  }
  LOCAL_BARRIER();
  if (lid < GROW_PARTIALS) partials[set * GROW_PARTIALS + lid] = scratch[lid * GROW_GROUP];
}
)CLC";
