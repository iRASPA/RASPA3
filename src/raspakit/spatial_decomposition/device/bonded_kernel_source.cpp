module;

module spatial_decomposition_device_kernels;

// Device source (kernel dialect of kernel_sources.ixx) of the per-molecule terms: the Ewald self and intramolecular exclusion
// corrections (Interactions::addChargeSelfEnergy / addIntraMolecularChargeExclusionGradient), the bonds, bends,
// torsions and improper torsions (the potentialEnergyGradientStrain functions of the intramolecular potentials,
// transcribed case by case), the intramolecular Lennard-Jones and Coulomb pairs (VanDerWaalsPotential,
// CoulombPotential), and the atomic-to-molecular virial correction of the non-bonded gradients.
//
// Two kernels. bondedTerms: one work-item per term instance (a term of a component in one of its molecules; the
// instances of a molecule are consecutive and ordered by kind, so that the work-items of a SIMD group mostly run
// the same code), which evaluates the term once and writes the gradient on each of its atoms to a per-instance
// slot of `termGradient`, with the energies reduced per work-group. bondedAtoms: one work-item per slot (atom),
// which evaluates the exclusion pairs with the other atoms of its molecule, the virial correction, and gathers
// the gradients of the terms its atom takes part in (listed per atom of the component). No atomics; the gradient
// is added to the force of the slot, which at this point holds the pair + mesh gradient (the kernels run after the
// pair kernel and the mesh interpolation). The positions are the float positions relative to the first atom of
// the molecule, indexed by the atom's index in the system (so the atoms of a molecule are consecutive), packed by
// the host in double so that the stiff bonded terms do not see the rounding of the absolute positions; the w
// component carries the charge.
const char* const deviceKernelBondedSource = R"CLC(
#define BONDED_GROUP 64
#define ATOM_PARTIALS 20
#define TERM_GROUP 64
#define TERM_PARTIALS 8
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
  uint numberOfInstances;  // term instances over all molecules
  uint atomPartialOffset;  // first float of the per-atom partials in the partial buffer (after the term partials)
} BondedParameters;

typedef struct
{
  uint firstAtom;      // index in the system of the first atom of the molecule
  uint numberOfAtoms;
  uint atomOffset;     // offset of the component's atoms in the per-atom tables
  uint termOffset;     // offset of the component's terms in the term table
  uint instanceBase;   // first term instance of the molecule
  uint gradientBase;   // first gradient slot of the molecule
  uint padding[2];
} MoleculeInfo;

// erf(x)/x - 2/sqrt(pi) without cancellation: the Taylor series for small x (x <= 0.5: the omitted terms are
// below 1e-8 relative), the difference otherwise (at most a factor 12 of cancellation at x = 0.5)
DEVICE_FUNCTION float erfOverXMinusLimit(float x)
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
DEVICE_FUNCTION float bondTerm(uint type, const float* P, float3 posA, float3 posB, float3* gA, float3* gB)
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
DEVICE_FUNCTION float bendTerm(uint type, const float* P, float3 posA, float3 posB, float3 posC, float3* gA, float3* gB,
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
DEVICE_FUNCTION float torsionTerm(uint type, const float* P, float3 posA, float3 posB, float3 posC, float3 posD, float3* gA,
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

DEVICE_FUNCTION void reducePartials(LOCAL float* scratch, uint count, GLOBAL float* RESTRICT partials)
{
  const uint lid = LOCAL_ID();
  const uint groupSize = LOCAL_SIZE();
  for (uint stride = groupSize / 2; stride > 0; stride >>= 1)
  {
    LOCAL_BARRIER();
    if (lid < stride)
    {
      for (uint q = 0; q < count; ++q) scratch[q * groupSize + lid] += scratch[q * groupSize + lid + stride];
    }
  }
  LOCAL_BARRIER();
  if (lid < count) partials[GROUP_ID() * count + lid] = scratch[lid * groupSize];
}

// Intramolecular pairs: kind 4 Lennard-Jones with P[0] = scaling 4 epsilon, P[1] = sigma^2; kind 5 Coulomb with
// P[0] = scaling C qA qB (VanDerWaalsPotential / CoulombPotential::potentialEnergyGradientStrain)
DEVICE_FUNCTION float pairTerm(uint kind, const float* P, float3 posA, float3 posB, float3* gA, float3* gB)
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

// One work-item per term instance: evaluates the term once and writes the gradient on each of its atoms to the
// gradient slots of the instance (gradientBase of the molecule + the term's offset within the component's block,
// one float4 per participating atom). Partials per work-group: bond, bend, torsion, improper torsion,
// intramolecular Lennard-Jones and Coulomb energies (indexed by kind).
KERNEL_GROUP_SIZE(TERM_GROUP)
void bondedTerms(GLOBAL const float4* RESTRICT relative,            // per atom (system order): x, y, z relative to
                                                                      // the first atom of the molecule, charge
                 GLOBAL const uint* RESTRICT instanceMolecule,      // molecule of every term instance
                 GLOBAL const MoleculeInfo* RESTRICT moleculeInfo,
                 GLOBAL const Term* RESTRICT terms,                 // the terms of the components
                 GLOBAL const uint* RESTRICT gradientOffset,        // per term: offset of its gradient slots
                 CONSTANT const BondedParameters* p,
                 GLOBAL float4* RESTRICT termGradient,
                 GLOBAL float* RESTRICT partials)
{
  LOCAL float scratch[TERM_PARTIALS * TERM_GROUP];
  const uint lid = LOCAL_ID();
  const uint g = GLOBAL_ID();
  float acc[TERM_PARTIALS];
  for (uint q = 0; q < TERM_PARTIALS; ++q) acc[q] = 0.0f;

  if (g < p->numberOfInstances)
  {
    const uint m = instanceMolecule[g];
    const MoleculeInfo info = moleculeInfo[m];
    const uint t = info.termOffset + (g - info.instanceBase);
    const Term term = terms[t];
    GLOBAL const float4* RESTRICT pos = relative + info.firstAtom;
    GLOBAL float4* RESTRICT out = termGradient + info.gradientBase + gradientOffset[t];
    const float3 posA = pos[term.atoms[0]].xyz;
    const float3 posB = pos[term.atoms[1]].xyz;
    float3 gA = FLOAT3(0.0f, 0.0f, 0.0f), gB = gA, gC = gA, gD = gA;
    float U;
    if (term.kind == 0)
    {
      U = bondTerm(term.type, term.parameters, posA, posB, &gA, &gB);
    }
    else if (term.kind == 1)
    {
      const float3 posC = pos[term.atoms[2]].xyz;
      U = bendTerm(term.type, term.parameters, posA, posB, posC, &gA, &gB, &gC);
      out[2] = FLOAT4(gC, 0.0f);
    }
    else if (term.kind <= 3)
    {
      const float3 posC = pos[term.atoms[2]].xyz;
      const float3 posD = pos[term.atoms[3]].xyz;
      U = torsionTerm(term.type, term.parameters, posA, posB, posC, posD, &gA, &gB, &gC, &gD);
      out[2] = FLOAT4(gC, 0.0f);
      out[3] = FLOAT4(gD, 0.0f);
    }
    else
    {
      U = pairTerm(term.kind, term.parameters, posA, posB, &gA, &gB);
    }
    out[0] = FLOAT4(gA, 0.0f);
    out[1] = FLOAT4(gB, 0.0f);
    acc[min(term.kind, (uint)(TERM_PARTIALS - 1))] = U;
  }

  for (uint q = 0; q < TERM_PARTIALS; ++q) scratch[q * TERM_GROUP + lid] = acc[q];
  reducePartials(scratch, TERM_PARTIALS, partials);
}

// One work-item per slot (atom): the Ewald self and exclusion corrections with the other atoms of its molecule,
// the virial correction of the non-bonded gradient, and the gradients of the term instances it takes part in,
// gathered from the slots the term kernel wrote. Partials per work-group: net-charge self term, reduced exclusion
// term (see below), the exclusion strain derivative and the virial correction (ax ay az bx by bz cx cy cz each).
KERNEL_GROUP_SIZE(BONDED_GROUP)
void bondedAtoms(GLOBAL const float4* RESTRICT relative,          // per atom (system order), see bondedTerms
                 GLOBAL const uint* RESTRICT slotMolecule,        // (molecule << 8) | index in the molecule, or NO_ATOM
                 GLOBAL const MoleculeInfo* RESTRICT moleculeInfo,
                 GLOBAL const float* RESTRICT massOfAtom,         // per atom (system order)
                 GLOBAL const uint* RESTRICT atomGradientStart,   // CSR offsets of the gradient slots per component atom
                 GLOBAL const uint* RESTRICT atomGradients,       // gradient slot (within the molecule's block)
                 GLOBAL const float4* RESTRICT termGradient,
                 CONSTANT const BondedParameters* p,
                 GLOBAL float4* RESTRICT force,
                 GLOBAL float* RESTRICT partials)
{
  LOCAL float scratch[ATOM_PARTIALS * BONDED_GROUP];
  const uint lid = LOCAL_ID();
  const uint slot = GLOBAL_ID();
  float acc[ATOM_PARTIALS];
  for (uint q = 0; q < ATOM_PARTIALS; ++q) acc[q] = 0.0f;

  const uint packed = (slot < p->numberOfSlots) ? slotMolecule[slot] : NO_ATOM;
  if (packed != NO_ATOM)
  {
    const uint m = packed >> 8;
    const uint a = packed & 0xFFu;
    const MoleculeInfo info = moleculeInfo[m];
    const uint first = info.firstAtom;
    const uint n = info.numberOfAtoms;
    GLOBAL const float4* RESTRICT pos = relative + first;
    GLOBAL const float* RESTRICT mass = massOfAtom + first;
    const float4 pa = pos[a];
    const float3 ra = pa.xyz;
    const float qa = pa.w;
    const bool charged = p->useCharge != 0;

    // Self energy and exclusion pairs with the other atoms of the molecule, center of mass. The self energy
    // (-C alpha/sqrt(pi) sum q^2) and the exclusion energy (-C sum_{a<b} qa qb erf(alpha r)/r) cancel almost
    // completely for a neutral molecule; their sum is evaluated without the cancellation as
    //   -C sum_{a<b} qa qb [erf(alpha r)/r - 2 alpha/sqrt(pi)] - C alpha/sqrt(pi) (sum_a qa)^2
    // (exact identity): partial 0 holds the net-charge term, partial 1 the reduced pair term. The host
    // separates the self energy (evaluated in double) from the sum.
    float3 exclusionGradient = FLOAT3(0.0f, 0.0f, 0.0f);
    float3 com = FLOAT3(0.0f, 0.0f, 0.0f);
    float totalMass = 0.0f;
    float moleculeCharge = 0.0f;
    for (uint b = 0; b < n; ++b)
    {
      const float4 pb = pos[b];
      const float mb = mass[b];
      com += mb * pb.xyz;
      totalMass += mb;
      if (!charged) continue;
      moleculeCharge += pb.w;
      if (b == a) continue;
      const float3 dr = ra - pb.xyz;
      const float rr = dot(dr, dr);
      const float r = sqrt(rr);
      const float x = p->alpha * r;
      const float potential = erf(x) / r;
      const float gaussian = p->twoAlphaOverSqrtPi * exp(-x * x);
      const float firstDerivativeFactor = (gaussian - potential) / rr;
      const float prefactor = p->coulombFactor * qa * pb.w;
      const float gradientFactor = prefactor * firstDerivativeFactor;
      exclusionGradient -= gradientFactor * dr;
      if (b > a)
      {
        acc[1] -= prefactor * p->alpha * erfOverXMinusLimit(x);
        const float3 f = -gradientFactor * dr;
        acc[2] += f.x * dr.x;
        acc[3] += f.x * dr.y;
        acc[4] += f.x * dr.z;
        acc[5] += f.y * dr.x;
        acc[6] += f.y * dr.y;
        acc[7] += f.y * dr.z;
        acc[8] += f.z * dr.x;
        acc[9] += f.z * dr.y;
        acc[10] += f.z * dr.z;
      }
    }
    com /= totalMass;
    if (charged) acc[0] = -p->selfPrefactor * qa * moleculeCharge;

    // virial correction of the non-bonded gradient (pairs + mesh from the force buffer, plus the exclusions)
    float4 f = force[slot];
    const float3 nonbonded = f.xyz + exclusionGradient;
    const float3 arm = ra - com;
    acc[11] = arm.x * nonbonded.x;
    acc[12] = arm.x * nonbonded.y;
    acc[13] = arm.x * nonbonded.z;
    acc[14] = arm.y * nonbonded.x;
    acc[15] = arm.y * nonbonded.y;
    acc[16] = arm.y * nonbonded.z;
    acc[17] = arm.z * nonbonded.x;
    acc[18] = arm.z * nonbonded.y;
    acc[19] = arm.z * nonbonded.z;

    // the gradients of the terms of this atom
    float3 bonded = FLOAT3(0.0f, 0.0f, 0.0f);
    GLOBAL const float4* RESTRICT gradients = termGradient + info.gradientBase;
    const uint begin = atomGradientStart[info.atomOffset + a];
    const uint end = atomGradientStart[info.atomOffset + a + 1];
    for (uint t = begin; t < end; ++t) bonded += gradients[atomGradients[t]].xyz;

    f.x = nonbonded.x + bonded.x;
    f.y = nonbonded.y + bonded.y;
    f.z = nonbonded.z + bonded.z;
    force[slot] = f;
  }

  for (uint q = 0; q < ATOM_PARTIALS; ++q) scratch[q * BONDED_GROUP + lid] = acc[q];
  reducePartials(scratch, ATOM_PARTIALS, partials + p->atomPartialOffset);
}
)CLC";
