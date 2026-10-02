module;

module spatial_decomposition_opencl_bonded;

// OpenCL C (1.2) source of the per-molecule terms on the device: the Ewald self and intramolecular exclusion
// corrections (Interactions::addChargeSelfEnergy / addIntraMolecularChargeExclusionGradient), the bonds, bends,
// torsions and improper torsions (the potentialEnergyGradientStrain functions of the intramolecular potentials,
// transcribed case by case), the intramolecular Lennard-Jones and Coulomb pairs (VanDerWaalsPotential,
// CoulombPotential), and the atomic-to-molecular virial correction of the non-bonded gradients.
//
// One work-item per slot (atom). It gathers everything that acts on its atom: the exclusion pairs with the other
// atoms of its molecule, and the bonded terms it takes part in, listed per atom of the component (a term is thus
// evaluated once per atom it involves; the energy and the exclusion strain are counted by one of them). No
// atomics and no reduction across work-items; the gradient is added to the force of the slot, which at this point
// holds the pair + mesh gradient (the kernel runs after the pair kernel and the mesh interpolation). The positions
// are the float positions relative to the first atom of the molecule, packed by the host in double, so that the
// stiff bonded terms do not see the rounding of the absolute positions.
const char* const openclBondedKernelSource = R"CLC(
#define BONDED_GROUP 64
#define BONDED_PARTIALS 26
#define NO_ATOM 0xFFFFFFFFu
#define TWO_PI 6.283185307179586f
#define TWO_OVER_SQRT_PI 1.1283791670955126f
#define RADIANS_TO_DEGREES 57.29577951308232f

typedef struct
{
  uint kind;  // 0 bond, 1 bend, 2 torsion, 3 improper torsion, 4 intramolecular Lennard-Jones, 5 intramolecular Coulomb
  uint type;  // the BondType / BendType / TorsionType value
  uint atoms[4];
  float parameters[6];
} Term;

typedef struct
{
  float alpha;
  float twoAlphaOverSqrtPi;
  float selfPrefactor;  // C alpha / sqrt(pi)
  float coulombFactor;
  uint useCharge;
  uint numberOfSlots;
  uint padding[2];
} BondedParameters;

// erf(x)/x - 2/sqrt(pi) without cancellation: the Taylor series for small x (x <= 0.5: the omitted terms are
// below 1e-8 relative), the difference otherwise (at most a factor 12 of cancellation at x = 0.5)
inline float erfOverXMinusLimit(float x)
{
  const float x2 = x * x;
  if (x2 <= 0.25f)
  {
    const float series = -1.0f / 3.0f +
                         x2 * (0.1f + x2 * (-1.0f / 42.0f +
                                            x2 * (1.0f / 216.0f + x2 * (-1.0f / 1320.0f + x2 * (1.0f / 9360.0f)))));
    return TWO_OVER_SQRT_PI * x2 * series;
  }
  return (erf(x) - TWO_OVER_SQRT_PI * x) / x;
}

// BondPotential (distancePotentialEnergyGradientStrain): energy and DF with gradient_A = DF dr, dr = posA - posB
inline float bondTerm(uint type, const float* P, float3 posA, float3 posB, float3* gA, float3* gB)
{
  const float3 dr = posA - posB;
  const float rr = dot(dr, dr);
  const float r = sqrt(rr);
  float U = 0.0f;
  float DF = 0.0f;
  float temp, temp2, r1;
  switch (type)
  {
    case 2:  // Harmonic
      U = 0.5f * P[0] * (r - P[1]) * (r - P[1]);
      DF = P[0] * (r - P[1]) / r;
      break;
    case 3:  // CoreShellSpring
      U = 0.5f * P[0] * rr;
      DF = P[0];
      break;
    case 4:  // Morse
      temp = exp(P[1] * (P[2] - r));
      U = P[0] * ((1.0f - temp) * (1.0f - temp) - 1.0f);
      DF = 2.0f * P[0] * P[1] * (1.0f - temp) * temp / r;
      break;
    case 5:  // LJ_12_6
      temp = 1.0f / (rr * rr * rr);
      U = P[0] * temp * temp - P[1] * temp;
      DF = 6.0f * (P[1] * temp - 2.0f * P[0] * temp * temp) / rr;
      break;
    case 6:  // LennardJones
      temp = (P[1] / rr) * (P[1] / rr) * (P[1] / rr);
      U = 4.0f * P[0] * (temp * (temp - 1.0f));
      DF = 24.0f * P[0] * (temp * (1.0f - 2.0f * temp)) / rr;
      break;
    case 7:  // Buckingham
      temp = P[2] / (rr * rr * rr);
      temp2 = P[0] * exp(-P[1] * r);
      U = -temp + temp2;
      DF = (6.0f / rr) * temp - P[1] * temp2 / r;
      break;
    case 8:  // RestrainedHarmonic
      r1 = r - P[1];
      temp = min(fabs(r1), P[2]);
      U = 0.5f * P[0] * temp * temp + P[0] * P[2] * max(fabs(r1) - P[2], 0.0f);
      DF = -P[0] * copysign(temp, r1) / r;
      break;
    case 9:  // Quartic
      temp = r - P[1];
      temp2 = temp * temp;
      U = 0.5f * P[0] * temp2 + (1.0f / 3.0f) * P[2] * temp * temp2 + 0.25f * P[3] * temp2 * temp2;
      DF = temp * (P[0] + P[2] * temp + P[3] * temp2) / r;
      break;
    case 10:  // CFF_Quartic
      temp = r - P[1];
      temp2 = temp * temp;
      U = P[0] * temp2 + P[2] * temp * temp2 + P[3] * temp2 * temp2;
      DF = temp * (2.0f * P[0] + 3.0f * P[2] * temp + 4.0f * P[3] * temp2) / r;
      break;
    case 11:  // MM3
      temp = r - P[1];
      temp2 = temp * temp;
      U = P[0] * temp2 * (1.0f - 2.55f * temp + (7.0f / 12.0f) * 2.55f * 2.55f * temp2);
      DF = P[0] * (2.0f + 2.55f * (4.0f * 2.55f * (7.0f / 12.0f) * temp - 3.0f) * temp) * temp / r;
      break;
    default:  // None, Fixed
      break;
  }
  if (!(r > 0.0f)) DF = 0.0f;
  *gA = DF * dr;
  *gB = -DF * dr;
  return U;
}

// BendPotential::potentialEnergyGradientStrain
inline float bendTerm(uint type, const float* P, float3 posA, float3 posB, float3 posC, float3* gA, float3* gB,
                      float3* gC)
{
  float3 dr_ab = posA - posB;
  const float r_ab = length(dr_ab);
  dr_ab /= r_ab;
  float3 dr_cb = posC - posB;
  const float r_cb = length(dr_cb);
  dr_cb /= r_cb;
  const float cos_theta = clamp(dot(dr_ab, dr_cb), -1.0f, 1.0f);
  const float theta = acos(cos_theta);
  const float DTDX = -1.0f / sqrt(1.0f - cos_theta * cos_theta);
  float U = 0.0f;
  float DF = 0.0f;
  float temp, temp2;
  switch (type)
  {
    case 2:  // Harmonic
    case 3:  // CoreShell
      U = 0.5f * P[0] * (theta - P[1]) * (theta - P[1]);
      DF = P[0] * (theta - P[1]) * DTDX;
      break;
    case 4:  // Quartic
      temp = theta - P[1];
      temp2 = temp * temp;
      U = 0.5f * P[0] * temp2 + (1.0f / 3.0f) * P[2] * temp * temp2 + 0.25f * P[3] * temp2 * temp2;
      DF = (P[0] * temp + P[2] * temp2 + P[3] * temp * temp2) * DTDX;
      break;
    case 5:  // CFF_Quartic
      temp = theta - P[1];
      temp2 = temp * temp;
      U = P[0] * temp2 + P[2] * temp * temp2 + P[3] * temp2 * temp2;
      DF = (2.0f * P[0] * temp + 3.0f * P[2] * temp2 + 4.0f * P[3] * temp * temp2) * DTDX;
      break;
    case 6:  // HarmonicCosine
      temp = cos_theta - P[1];
      U = 0.5f * P[0] * temp * temp;
      DF = P[0] * temp;
      break;
    case 7:  // Cosine
      temp = P[1] * theta - P[2];
      U = P[0] * (1.0f + cos(temp));
      DF = -P[0] * P[1] * sin(temp) * DTDX;
      break;
    case 8:  // Tafipolsky
      U = 0.5f * P[0] * (1.0f + cos(theta)) * (1.0f + cos(2.0f * theta));
      DF = P[0] * cos_theta * (2.0f + 3.0f * cos_theta);
      break;
    case 9:   // MM3
    case 10:  // MM3_inplane
      temp = (theta - P[1]) * RADIANS_TO_DEGREES;
      temp2 = temp * temp;
      U = P[0] * temp2 * (1.0f - 0.014f * temp + 5.6e-5f * temp2 - 7.0e-7f * temp * temp2 + 2.2e-8f * temp2 * temp2);
      DF = P[0] * RADIANS_TO_DEGREES *
           (2.0f - (3.0f * 0.014f - (4.0f * 5.6e-5f - (5.0f * 7.0e-7f - 6.0f * 2.2e-8f * temp) * temp) * temp) * temp) *
           temp * DTDX;
      break;
    default:  // Rigid, Fixed
      break;
  }
  *gA = DF * (dr_cb - cos_theta * dr_ab) / r_ab;
  *gC = DF * (dr_ab - cos_theta * dr_cb) / r_cb;
  *gB = -(*gA + *gC);
  return U;
}

// TorsionPotential::potentialEnergyGradientStrain
inline float torsionTerm(uint type, const float* P, float3 posA, float3 posB, float3 posC, float3 posD, float3* gA,
                         float3* gB, float3* gC, float3* gD)
{
  const float3 Dab = posA - posB;
  const float3 Dcb = posC - posB;
  const float rbc = length(Dcb);
  const float3 Dcb_unit = Dcb / rbc;
  const float3 Ddc = posD - posC;
  const float dot_ab = dot(Dab, Dcb_unit);
  const float dot_cd = dot(Ddc, Dcb_unit);
  float3 dr = Dab - dot_ab * Dcb_unit;
  const float r = length(dr);
  dr /= r;
  float3 ds = Ddc - dot_cd * Dcb_unit;
  const float s = length(ds);
  ds /= s;
  const float cos_phi = clamp(dot(dr, ds), -1.0f, 1.0f);
  const float cos_phi2 = cos_phi * cos_phi;
  const float3 Pb = cross(Dab, Dcb_unit);
  const float3 Pc = cross(Dcb_unit, Ddc);
  const float orientation = dot(Dcb_unit, cross(Pb, Pc));
  float phi = copysign(acos(cos_phi), orientation);
  const float sin_phi = copysign(max(1.0e-8f, fabs(sin(phi))), sin(phi));
  float U = 0.0f;
  float DF = 0.0f;
  float temp, shifted_cos_phi, shifted_cos_phi2, shifted_sin_phi;
  switch (type)
  {
    case 1:  // Harmonic
      phi -= P[1];
      phi -= rint(phi / TWO_PI) * TWO_PI;
      U = 0.5f * P[0] * phi * phi;
      DF = -P[0] * phi / sin_phi;
      break;
    case 2:  // HarmonicCosine
      U = 0.5f * P[0] * (cos_phi - P[1]) * (cos_phi - P[1]);
      DF = P[0] * (cos_phi - P[1]);
      break;
    case 3:   // ThreeCosine
    case 12:  // MM3
      U = 0.5f * P[0] * (1.0f + cos_phi) + P[1] * (1.0f - cos_phi2) +
          0.5f * P[2] * (1.0f - 3.0f * cos_phi + 4.0f * cos_phi * cos_phi2);
      DF = 0.5f * P[0] - 2.0f * P[1] * cos_phi + 1.5f * P[2] * (4.0f * cos_phi2 - 1.0f);
      break;
    case 4:  // RyckaertBellemans
      U = P[0] - P[1] * cos_phi + P[2] * cos_phi2 - P[3] * cos_phi * cos_phi2 + P[4] * cos_phi2 * cos_phi2 -
          P[5] * cos_phi2 * cos_phi2 * cos_phi;
      DF = -P[1] + 2.0f * P[2] * cos_phi - 3.0f * P[3] * cos_phi2 + 4.0f * P[4] * cos_phi2 * cos_phi -
           5.0f * P[5] * cos_phi2 * cos_phi2;
      break;
    case 5:  // TraPPE
      U = P[0] + (1.0f + cos_phi) * (P[1] + P[3] - 2.0f * (cos_phi - 1.0f) * (P[2] - 2.0f * P[3] * cos_phi));
      DF = P[1] - 4.0f * P[2] * cos_phi + 3.0f * P[3] * (4.0f * cos_phi2 - 1.0f);
      break;
    case 6:  // TraPPE_Extended
      U = P[0] - P[2] + P[4] + (P[1] - 3.0f * P[3]) * cos_phi + (2.0f * P[2] - 8.0f * P[4]) * cos_phi2 +
          4.0f * P[3] * cos_phi2 * cos_phi + 8.0f * P[4] * cos_phi2 * cos_phi2;
      DF = P[1] - 3.0f * P[3] + 4.0f * (P[2] - 4.0f * P[4]) * cos_phi + 12.0f * P[3] * cos_phi2 +
           32.0f * P[4] * cos_phi2 * cos_phi;
      break;
    case 7:  // ModifiedTraPPE
      phi -= P[4];
      phi -= rint(phi / TWO_PI) * TWO_PI;
      shifted_cos_phi = cos(phi);
      shifted_sin_phi = sin(phi);
      shifted_cos_phi2 = shifted_cos_phi * shifted_cos_phi;
      U = P[0] + P[1] + P[3] + (P[1] - 3.0f * P[3]) * shifted_cos_phi - 2.0f * P[2] * shifted_cos_phi2 +
          4.0f * P[3] * shifted_cos_phi * shifted_cos_phi2;
      DF = ((P[1] - 3.0f * P[3]) * shifted_sin_phi - 4.0f * P[2] * shifted_cos_phi * shifted_sin_phi +
            12.0f * P[3] * shifted_cos_phi2 * shifted_sin_phi) /
           sin_phi;
      break;
    case 8:  // CVFF
      temp = P[1] * phi - P[2];
      U = P[0] * (1.0f + cos(temp));
      DF = P[0] * P[1] * sin(temp) / sin_phi;
      break;
    case 9:  // CFF
      U = P[0] * (1.0f - cos_phi) + 2.0f * P[1] * (1.0f - cos_phi2) +
          P[2] * (1.0f + 3.0f * cos_phi - 4.0f * cos_phi * cos_phi2);
      DF = -P[0] - 4.0f * P[1] * cos_phi + 3.0f * P[2] * (1.0f - 4.0f * cos_phi2);
      break;
    case 10:  // CFF2
      U = P[0] * (1.0f + cos_phi) + P[2] + cos_phi * (-3.0f * P[2] + 2.0f * cos_phi * (P[1] + 2.0f * P[2] * cos_phi));
      DF = P[0] - 3.0f * P[2] + 4.0f * cos_phi * (P[1] + 3.0f * P[2] * cos_phi);
      break;
    case 11:  // OPLS
      U = 0.5f * (P[0] + (1.0f + cos_phi) * (P[1] + P[3] - 2.0f * (cos_phi - 1.0f) * (P[2] - 2.0f * P[3] * cos_phi)));
      DF = 0.5f * P[1] - 2.0f * P[2] * cos_phi + 1.5f * P[3] * (4.0f * cos_phi2 - 1.0f);
      break;
    case 13:  // FourierSeries
      U = 0.5f * (P[0] + 2.0f * P[1] + P[2] + P[4] + 2.0f * P[5] + (P[0] - 3.0f * P[2] + 5.0f * P[4]) * cos_phi -
                  2.0f * (P[1] - 4.0f * P[3] + 9.0f * P[5]) * cos_phi2 +
                  4.0f * (P[2] - 5.0f * P[4]) * cos_phi2 * cos_phi - 8.0f * (P[3] - 6.0f * P[5]) * cos_phi2 * cos_phi2 +
                  16.0f * P[4] * cos_phi2 * cos_phi2 * cos_phi - 32.0f * P[5] * cos_phi2 * cos_phi2 * cos_phi2);
      DF = 0.5f * (P[0] - 3.0f * P[2] + 5.0f * P[4]) - 2.0f * (P[1] - 4.0f * P[3] + 9.0f * P[5]) * cos_phi +
           6.0f * (P[2] - 5.0f * P[4]) * cos_phi2 - 16.0f * (P[3] - 6.0f * P[5]) * cos_phi2 * cos_phi +
           40.0f * P[4] * cos_phi2 * cos_phi2 - 96.0f * P[5] * cos_phi2 * cos_phi * cos_phi;
      break;
    case 14:  // FourierSeries2
      U = 0.5f * (P[2] + 2.0f * P[3] + P[4] - 3.0f * P[2] * cos_phi + 5.0f * P[4] * cos_phi + P[0] * (1.0f + cos_phi) +
                  2.0f * (P[1] - P[1] * cos_phi2 +
                          cos_phi2 * (P[5] * (3.0f - 4.0f * cos_phi2) * (3.0f - 4.0f * cos_phi2) +
                                      4.0f * P[3] * (cos_phi2 - 1.0f) +
                                      2.0f * cos_phi * (P[2] + P[4] * (4.0f * cos_phi2 - 5.0f)))));
      DF = 0.5f * P[0] + P[2] * (6.0f * cos_phi2 - 1.5f) +
           P[4] * (2.5f - 30.0f * cos_phi2 + 40.0f * cos_phi2 * cos_phi2) +
           cos_phi * (-2.0f * P[1] + P[3] * (16.0f * cos_phi2 - 8.0f) +
                      P[5] * (18.0f - 96.0f * cos_phi2 + 96.0f * cos_phi2 * cos_phi2));
      break;
    case 16:  // Polynomial: U = sum_i p_i c^i, DF = sum_i i p_i c^(i-1)
      U = P[5];
      DF = 0.0f;
      for (int i = 4; i >= 0; --i)
      {
        DF = DF * cos_phi + U;
        U = U * cos_phi + P[i];
      }
      break;
    default:  // Fixed, CVFFBlocked
      break;
  }
  const float d = dot_ab / rbc;
  const float e = dot_cd / rbc;
  const float3 dtA = (ds - cos_phi * dr) / r;
  const float3 dtD = (dr - cos_phi * ds) / s;
  const float3 dtB = dtA * (d - 1.0f) + e * dtD;
  const float3 dtC = -dtD * (e + 1.0f) - d * dtA;
  *gA = DF * dtA;
  *gB = DF * dtB;
  *gC = DF * dtC;
  *gD = DF * dtD;
  return U;
}

inline void reducePartials(__local float* scratch, uint count, __global float* restrict partials)
{
  const uint lid = get_local_id(0);
  const uint groupSize = get_local_size(0);
  for (uint stride = groupSize / 2; stride > 0; stride >>= 1)
  {
    barrier(CLK_LOCAL_MEM_FENCE);
    if (lid < stride)
    {
      for (uint q = 0; q < count; ++q) scratch[q * groupSize + lid] += scratch[q * groupSize + lid + stride];
    }
  }
  barrier(CLK_LOCAL_MEM_FENCE);
  if (lid < count) partials[get_group_id(0) * count + lid] = scratch[lid * groupSize];
}

// Intramolecular pairs: kind 4 Lennard-Jones with P[0] = scaling 4 epsilon, P[1] = sigma^2; kind 5 Coulomb with
// P[0] = scaling C qA qB (VanDerWaalsPotential / CoulombPotential::potentialEnergyGradientStrain)
inline float pairTerm(uint kind, const float* P, float3 posA, float3 posB, float3* gA, float3* gB)
{
  const float3 dr = posA - posB;
  const float rr = dot(dr, dr);
  float U, DF;
  if (kind == 4)
  {
    const float s = P[1] / rr;
    const float t = s * s * s;
    U = P[0] * (t * (t - 1.0f));
    DF = 6.0f * P[0] * (t * (1.0f - 2.0f * t)) / rr;
  }
  else
  {
    const float r = sqrt(rr);
    U = P[0] / r;
    DF = -P[0] / (rr * r);
  }
  *gA = DF * dr;
  *gB = -DF * dr;
  return U;
}

// Partials per work-group: net-charge self term, reduced exclusion term (see the kernel), bond, bend, torsion,
// improper torsion energies; the exclusion strain
// derivative and the virial correction (ax ay az bx by bz cx cy cz each); intramolecular Lennard-Jones and Coulomb
// energies.
__kernel __attribute__((reqd_work_group_size(BONDED_GROUP, 1, 1)))
void bondedAtoms(__global const float4* restrict position,         // x, y, z, charge per slot
                 __global const float4* restrict relative,         // position relative to the molecule's first atom
                 __global const uint* restrict typeOf,             // pseudo-atom type per slot
                 __global const uint* restrict slotMolecule,       // (molecule << 8) | index in the molecule, or NO_ATOM
                 __global const uint4* restrict moleculeInfo,      // first original atom, atoms, atom offset of the component
                 __global const uint* restrict slotOfOriginal,     // slot of every atom in the original order
                 __global const uint* restrict atomTermStart,      // CSR offsets of the terms per component atom
                 __global const uint* restrict atomTerms,          // (term << 2) | role
                 __global const Term* restrict terms,
                 __global const float* restrict massOfType,
                 __constant const BondedParameters* p,
                 __global float4* restrict force,
                 __global float* restrict partials)
{
  __local float scratch[BONDED_PARTIALS * BONDED_GROUP];
  const uint lid = get_local_id(0);
  const uint slot = get_global_id(0);
  float acc[BONDED_PARTIALS];
  for (uint q = 0; q < BONDED_PARTIALS; ++q) acc[q] = 0.0f;

  const uint packed = (slot < p->numberOfSlots) ? slotMolecule[slot] : NO_ATOM;
  if (packed != NO_ATOM)
  {
    const uint m = packed >> 8;
    const uint a = packed & 0xFFu;
    const uint4 info = moleculeInfo[m];
    const uint first = info.x;
    const uint n = info.y;
    const uint atomOffset = info.z;
    const float3 ra = relative[slot].xyz;
    const float qa = position[slot].w;
    const bool charged = p->useCharge != 0;

    // Self energy and exclusion pairs with the other atoms of the molecule, center of mass. The self energy
    // (-C alpha/sqrt(pi) sum q^2) and the exclusion energy (-C sum_{a<b} qa qb erf(alpha r)/r) cancel almost
    // completely for a neutral molecule; their sum is evaluated without the cancellation as
    //   -C sum_{a<b} qa qb [erf(alpha r)/r - 2 alpha/sqrt(pi)] - C alpha/sqrt(pi) (sum_a qa)^2
    // (exact identity): partial 0 holds the net-charge term, partial 1 the reduced pair term. The host
    // separates the self energy (evaluated in double) from the sum.
    float3 exclusionGradient = (float3)(0.0f, 0.0f, 0.0f);
    float3 com = (float3)(0.0f, 0.0f, 0.0f);
    float totalMass = 0.0f;
    float moleculeCharge = 0.0f;
    for (uint b = 0; b < n; ++b)
    {
      const uint sb = slotOfOriginal[first + b];
      const float3 rb = relative[sb].xyz;
      const float mass = massOfType[typeOf[sb]];
      com += mass * rb;
      totalMass += mass;
      if (!charged) continue;
      const float qb = position[sb].w;
      moleculeCharge += qb;
      if (b == a) continue;
      const float3 dr = ra - rb;
      const float rr = dot(dr, dr);
      const float r = sqrt(rr);
      const float x = p->alpha * r;
      const float potential = erf(x) / r;
      const float gaussian = p->twoAlphaOverSqrtPi * exp(-x * x);
      const float firstDerivativeFactor = (gaussian - potential) / rr;
      const float prefactor = p->coulombFactor * qa * qb;
      const float gradientFactor = prefactor * firstDerivativeFactor;
      exclusionGradient -= gradientFactor * dr;
      if (b > a)
      {
        acc[1] -= prefactor * p->alpha * erfOverXMinusLimit(x);
        const float3 f = -gradientFactor * dr;
        acc[6] += f.x * dr.x;
        acc[7] += f.x * dr.y;
        acc[8] += f.x * dr.z;
        acc[9] += f.y * dr.x;
        acc[10] += f.y * dr.y;
        acc[11] += f.y * dr.z;
        acc[12] += f.z * dr.x;
        acc[13] += f.z * dr.y;
        acc[14] += f.z * dr.z;
      }
    }
    com /= totalMass;
    if (charged) acc[0] = -p->selfPrefactor * qa * moleculeCharge;

    // virial correction of the non-bonded gradient (pairs + mesh from the force buffer, plus the exclusions)
    float4 f = force[slot];
    const float3 nonbonded = f.xyz + exclusionGradient;
    const float3 arm = ra - com;
    acc[15] = arm.x * nonbonded.x;
    acc[16] = arm.x * nonbonded.y;
    acc[17] = arm.x * nonbonded.z;
    acc[18] = arm.y * nonbonded.x;
    acc[19] = arm.y * nonbonded.y;
    acc[20] = arm.y * nonbonded.z;
    acc[21] = arm.z * nonbonded.x;
    acc[22] = arm.z * nonbonded.y;
    acc[23] = arm.z * nonbonded.z;

    // bonded terms of this atom
    float3 bonded = (float3)(0.0f, 0.0f, 0.0f);
    const uint begin = atomTermStart[atomOffset + a];
    const uint end = atomTermStart[atomOffset + a + 1];
    for (uint t = begin; t < end; ++t)
    {
      const uint reference = atomTerms[t];
      const uint role = reference & 3u;
      const Term term = terms[reference >> 2];
      const float3 posA = relative[slotOfOriginal[first + term.atoms[0]]].xyz;
      const float3 posB = relative[slotOfOriginal[first + term.atoms[1]]].xyz;
      float3 gA = (float3)(0.0f, 0.0f, 0.0f), gB = gA, gC = gA, gD = gA;
      float U;
      if (term.kind == 0)
      {
        U = bondTerm(term.type, term.parameters, posA, posB, &gA, &gB);
      }
      else if (term.kind == 1)
      {
        const float3 posC = relative[slotOfOriginal[first + term.atoms[2]]].xyz;
        U = bendTerm(term.type, term.parameters, posA, posB, posC, &gA, &gB, &gC);
      }
      else if (term.kind <= 3)
      {
        const float3 posC = relative[slotOfOriginal[first + term.atoms[2]]].xyz;
        const float3 posD = relative[slotOfOriginal[first + term.atoms[3]]].xyz;
        U = torsionTerm(term.type, term.parameters, posA, posB, posC, posD, &gA, &gB, &gC, &gD);
      }
      else
      {
        U = pairTerm(term.kind, term.parameters, posA, posB, &gA, &gB);
      }
      bonded += (role == 0) ? gA : (role == 1) ? gB : (role == 2) ? gC : gD;
      if (role == 0) acc[(term.kind <= 3) ? 2 + term.kind : 20 + term.kind] += U;
    }

    f.x = nonbonded.x + bonded.x;
    f.y = nonbonded.y + bonded.y;
    f.z = nonbonded.z + bonded.z;
    force[slot] = f;
  }

  for (uint q = 0; q < BONDED_PARTIALS; ++q) scratch[q * BONDED_GROUP + lid] = acc[q];
  reducePartials(scratch, BONDED_PARTIALS, partials);
}
)CLC";
