module;

module spatial_decomposition_device_kernels;

// Device source (kernel dialect of kernel_sources.ixx) of the polymer-bond overlay of the Polymerization driver:
// the device transcription of Interactions::computePolymerLinkGradient (interactions/polymer_link.cpp). The
// source is compiled appended to the bonded source (bonded_kernel_source.cpp), whose term functions (bondTerm,
// bendTerm, torsionTerm) and reducePartials it reuses; its own names carry the 'polymer' prefix.
//
// The overlay is an addition on top of the complete pair + mesh + bonded gradients, evaluated after the bonded
// kernels: the bonds formed during the run with their bends and torsions (ramped in with a factor lambda per
// bond), the restraints of the stereocentres, the pair corrections of the pairs whose non-bonded treatment the
// new bonds change (value(wanted) - value(counted) for a pair the pair lists still hold, the whole ramped value
// for a pair the lists leave out; the pair kernel and the mesh are left untouched), and the
// 'replacement - original' overrides of the monomer terms that contain a retyped atom. The terms refer to the
// atoms by their index in the system; the geometry of a term is read from the wrapped slot positions of the pair
// lists (through the slot of every atom) and imaged consecutively along the chain of the term, as the host does
// with the stored positions.
//
// Two kernels. polymerTerms: one work-item per term evaluates it and writes the gradient on each of its atoms to
// the four gradient slots of the term; the energies (by slot of RunningEnergy) and the strain derivative of the
// inter-molecular terms are reduced per work-group. polymerAtoms: one work-item per participating atom gathers
// the gradients of its terms (listed per atom, the inter-molecular ones first), adds the sum to the force of the
// atom's slot and reduces the atomic-to-molecular virial correction (position relative to the molecule's center
// of mass, from the bonded kernels' relative positions and molecule centers, times the inter-molecular gradient).
// No atomics.
const char* const deviceKernelPolymerSource = R"CLC(
#define POLYMER_GROUP 64
#define POLYMER_TERM_PARTIALS 17
#define POLYMER_ATOM_PARTIALS 9

typedef struct
{
  uint kind;      // 0 link bond, 1 link bend, 2 link torsion, 3 override bond, 4 override bend, 5 override torsion,
                  // 6 pair correction, 7 stereocentre restraint
  uint type;      // the potential type (link; override: replacement); pair: the wanted treatment (0 full, 1 scaled,
                  // 2 excluded)
  uint type2;     // override: the original type; pair: the counted treatment
  uint ramp;      // index of the ramp factor of the term
  uint atoms[4];  // system atom indices
  uint flags;     // bit 0: inter-molecular (strain derivative and virial correction); bit 1: the pair energies go to
                  // the intramolecular slot; bit 2: the pair lists leave the pair out (the term supplies
                  // lambda wanted + (1 - lambda) counted instead of lambda (wanted - counted))
  uint padding[3];
  float parameters[6];   // link / replacement parameters; pair: wanted (f 4 eps, sigma^6, f shift, f_coulomb);
                         // stereocentre: k, chi V_0
  float parameters2[6];  // original parameters; pair: counted
} PolymerTerm;

typedef struct
{
  float cell[9];  // ax ay az bx by bz cx cy cz (columns of the cell matrix)
  float inverseCell[9];
  float cutOffVDWSquared;
  float cutOffChargeSquared;
  float alpha;
  float alphaSquared;
  float alphaOverSqrtPi;
  float twoAlphaOverSqrtPi;
  float coulombFactor;
  uint useCharge;
  uint orthorhombic;
  uint numberOfTerms;
  uint numberOfAtoms;  // participating atoms
  uint switchMode;     // Lennard-Jones switching of the full pairs (as the pair kernel): 0 none, 1 potential, 2 force
  float switchDistanceSquared;
  float switchDistance;
  float switchInverseWidth;
  float switchInverseCutOff3;
  float switchA12;
  float switchA6;
} PolymerParameters;

DEVICE_FUNCTION float3 polymerMinimumImage(float3 dr, CONSTANT const PolymerParameters* p)
{
  if (p->orthorhombic)
  {
    dr.x -= p->cell[0] * rint(dr.x * p->inverseCell[0]);
    dr.y -= p->cell[4] * rint(dr.y * p->inverseCell[4]);
    dr.z -= p->cell[8] * rint(dr.z * p->inverseCell[8]);
    return dr;
  }
  float3 s;
  s.x = p->inverseCell[0] * dr.x + p->inverseCell[3] * dr.y + p->inverseCell[6] * dr.z;
  s.y = p->inverseCell[1] * dr.x + p->inverseCell[4] * dr.y + p->inverseCell[7] * dr.z;
  s.z = p->inverseCell[2] * dr.x + p->inverseCell[5] * dr.y + p->inverseCell[8] * dr.z;
  s -= rint(s);
  float3 r;
  r.x = p->cell[0] * s.x + p->cell[3] * s.y + p->cell[6] * s.z;
  r.y = p->cell[1] * s.x + p->cell[4] * s.y + p->cell[7] * s.z;
  r.z = p->cell[2] * s.x + p->cell[5] * s.y + p->cell[8] * s.z;
  return r;
}

// One pair in one of its three treatments (polymer_link.cpp: pairValue): the energies by slot (van der Waals,
// real-space Coulomb, Ewald exclusion) and the gradient factor (gradient on the first atom = factor dr).
//   0 full:     the pair kernel's Lennard-Jones (switched when the force field switches) inside the VDW cutoff
//               and the Ewald real-space term inside the Coulomb cutoff, P = (4 eps, sigma^6, shift)
//   1 scaled:   f times the plain Lennard-Jones inside the VDW cutoff (P = f 4 eps, sigma^6, f shift) and inside
//               the Coulomb cutoff the bare f_c C qq / r minus the long-range part C qq erf(alpha r) / r the
//               Fourier sum counted (P[3] = f_c)
//   2 excluded: the long-range part the Fourier sum counted is removed, -C qq erf(alpha r) / r (no cutoff)
DEVICE_FUNCTION void polymerPairState(uint state, PRIVATE const float* P, float rr, float qq,
                                      CONSTANT const PolymerParameters* p, PRIVATE float* uVDW,
                                      PRIVATE float* uCoulomb, PRIVATE float* uExclusion, PRIVATE float* factor)
{
  *uVDW = 0.0f;
  *uCoulomb = 0.0f;
  *uExclusion = 0.0f;
  *factor = 0.0f;
  if (!(rr > 0.0f)) return;
  const float r = sqrt(rr);
  const float invR = 1.0f / r;
  const float invRR = invR * invR;
  const bool charged = p->useCharge != 0u && qq != 0.0f;
  if (state == 0u)
  {
    if (rr < p->cutOffVDWSquared && P[0] != 0.0f)
    {
      const float invRR3 = invRR * invRR * invRR;
      const float rri3 = P[1] * invRR3;
      const float rri6 = rri3 * rri3;
      float u = P[0] * (rri6 - rri3) - P[2];
      float f = 12.0f * P[0] * rri3 * (0.5f - rri3) * invRR;
      if (p->switchMode != 0u && rr > p->switchDistanceSquared)
      {
        if (p->switchMode == 1u)
        {
          const float x = (r - p->switchDistance) * p->switchInverseWidth;
          const float x2 = x * x;
          const float oneMinusX = 1.0f - x;
          const float sw = 1.0f + x2 * x * (-10.0f + x * (15.0f - 6.0f * x));
          const float dsw = -30.0f * p->switchInverseWidth * x2 * oneMinusX * oneMinusX;
          const float uPlain = P[0] * (rri6 - rri3);
          f = f * sw + uPlain * dsw * invR;
          u = uPlain * sw;
        }
        else
        {
          const float c6 = P[0] * P[1];
          const float c12 = c6 * P[1];
          const float invR3 = invR * invRR;
          const float invR4 = invR3 * invR;
          const float d3 = invR3 - p->switchInverseCutOff3;
          const float d6 = invR3 * invR3 - p->switchInverseCutOff3 * p->switchInverseCutOff3;
          u = p->switchA12 * c12 * d6 * d6 - p->switchA6 * c6 * d3 * d3;
          f = (-12.0f * p->switchA12 * c12 * d6 * invR3 * invR4 + 6.0f * p->switchA6 * c6 * d3 * invR4) * invR;
        }
      }
      *uVDW = u;
      *factor += f;
    }
    if (charged && rr < p->cutOffChargeSquared)
    {
      const float x = p->alpha * r;
      const float erfcValue = erfc(x);
      const float gauss = exp(-x * x);
      const float prefactor = p->coulombFactor * qq;
      *uCoulomb = prefactor * erfcValue * invR;
      *factor -= prefactor * (p->twoAlphaOverSqrtPi * gauss * invR + erfcValue * invRR) * invR;
    }
  }
  else if (state == 1u)
  {
    if (rr < p->cutOffVDWSquared && P[0] != 0.0f)
    {
      const float t = P[1] * invRR * invRR * invRR;
      *uVDW = P[0] * (t * (t - 1.0f)) - P[2];
      *factor += 6.0f * P[0] * (t * (1.0f - 2.0f * t)) * invRR;
    }
    if (charged && rr < p->cutOffChargeSquared)
    {
      const float x = p->alpha * r;
      const float potential = erf(x) * invR;
      const float gaussian = p->twoAlphaOverSqrtPi * exp(-x * x);
      const float prefactor = p->coulombFactor * qq;
      *uCoulomb = prefactor * (P[3] * invR - potential);
      *factor -= prefactor * (P[3] * invR * invRR + (gaussian - potential) * invRR);
    }
  }
  else if (charged)
  {
    const float x = p->alpha * r;
    const float potential = erf(x) * invR;
    const float gaussian = p->twoAlphaOverSqrtPi * exp(-x * x);
    const float prefactor = p->coulombFactor * qq;
    *uExclusion = -prefactor * potential;
    *factor -= prefactor * (gaussian - potential) * invRR;
  }
}

// ChiralVolume::compute: the normalised chiral volume V = (b_p x b_n) . b_z / (|b_p| |b_n| |b_z|) of the three
// substituent vectors from the centre and its gradients on the centre (g0) and the substituents (g1..g3).
DEVICE_FUNCTION float polymerChiralVolume(float3 bp, float3 bn, float3 bz, PRIVATE float3* g0, PRIVATE float3* g1,
                                          PRIVATE float3* g2, PRIVATE float3* g3)
{
  const float lp = length(bp);
  const float ln = length(bn);
  const float lz = length(bz);
  const float norm = lp * ln * lz;
  *g0 = FLOAT3(0.0f, 0.0f, 0.0f);
  *g1 = *g0;
  *g2 = *g0;
  *g3 = *g0;
  if (!(norm > 0.0f)) return 0.0f;
  const float3 pn = cross(bp, bn);
  const float3 nz = cross(bn, bz);
  const float3 zp = cross(bz, bp);
  const float value = dot(pn, bz) / norm;
  *g1 = nz / norm - (value / (lp * lp)) * bp;
  *g2 = zp / norm - (value / (ln * ln)) * bn;
  *g3 = pn / norm - (value / (lz * lz)) * bz;
  *g0 = -(*g1 + *g2 + *g3);
  return value;
}

// One work-item per term. Partials per work-group: the energies crossLink, bond, bend, torsion,
// moleculeMoleculeVDW, moleculeMoleculeCharge, ewald_exclusion, intraVDW (0..7) and the strain derivative of the
// inter-molecular terms (8..16: ax ay az bx by bz cx cy cz as sum_i (r_i - r_0) (x) g_i).
KERNEL_GROUP_SIZE(POLYMER_GROUP)
void polymerTerms(GLOBAL const float4* RESTRICT position,     // per slot: wrapped position, charge
                  GLOBAL const uint* RESTRICT slotOfAtom,     // per atom (system order): its slot
                  GLOBAL const PolymerTerm* RESTRICT terms,
                  GLOBAL const float* RESTRICT lambda,        // per ramp: the ramp factor
                  CONSTANT const PolymerParameters* p,
                  GLOBAL float4* RESTRICT termGradient,       // four slots per term
                  GLOBAL float* RESTRICT partials KERNEL_INDEX_ARGS)
{
  LOCAL_DECL(float scratch[POLYMER_TERM_PARTIALS * POLYMER_GROUP]);
  const uint lid = LOCAL_ID();
  const uint g = GLOBAL_ID();
  float acc[POLYMER_TERM_PARTIALS];
  for (uint q = 0; q < POLYMER_TERM_PARTIALS; ++q) acc[q] = 0.0f;

  if (g < p->numberOfTerms)
  {
    const PolymerTerm term = terms[g];
    const float lam = lambda[term.ramp];
    GLOBAL float4* RESTRICT out = termGradient + 4 * g;
    const uint kind = term.kind;
    const uint count = (kind == 0u || kind == 3u || kind == 6u) ? 2u : (kind == 1u || kind == 4u) ? 3u : 4u;
    const bool inter = (term.flags & 1u) != 0u;
    const bool unlisted = (term.flags & 4u) != 0u;

    // the geometry of the term: the first atom at its wrapped position, the others imaged consecutively along the
    // chain (the stereocentre: its substituents imaged from the centre)
    float4 pos4[4];
    float3 pos[4];
    pos4[0] = position[slotOfAtom[term.atoms[0]]];
    pos[0] = pos4[0].xyz;
    pos[1] = pos[0];
    pos[2] = pos[0];
    pos[3] = pos[0];
    for (uint i = 1; i < count; ++i)
    {
      pos4[i] = position[slotOfAtom[term.atoms[i]]];
      const float3 reference = (kind == 7u) ? pos[0] : pos[i - 1];
      pos[i] = reference + polymerMinimumImage(pos4[i].xyz - reference, p);
    }

    // named gradients (not an array: the device compilers mishandle the address of an indexed vector element)
    float3 gA = FLOAT3(0.0f, 0.0f, 0.0f), gB = gA, gC = gA, gD = gA;
    if (lam != 0.0f || unlisted)
    {
      float scale = lam;  // applied to the gradients below (the pair correction applies it in its factor)
      float P[6], Q[6];
      for (uint k = 0; k < 6; ++k)
      {
        P[k] = term.parameters[k];
        Q[k] = term.parameters2[k];
      }
      if (kind == 0u)
      {
        const float U = bondTerm(term.type, P, pos[0], pos[1], &gA, &gB);
        acc[0] = lam * U;
      }
      else if (kind == 1u)
      {
        const float U = bendTerm(term.type, P, pos[0], pos[1], pos[2], &gA, &gB, &gC);
        acc[0] = lam * U;
      }
      else if (kind == 2u)
      {
        const float U = torsionTerm(term.type, P, pos[0], pos[1], pos[2], pos[3], &gA, &gB, &gC, &gD);
        acc[0] = lam * U;
      }
      else if (kind == 3u)
      {
        float3 oA, oB;
        const float uNew = bondTerm(term.type, P, pos[0], pos[1], &gA, &gB);
        const float uOld = bondTerm(term.type2, Q, pos[0], pos[1], &oA, &oB);
        gA -= oA;
        gB -= oB;
        acc[1] = lam * (uNew - uOld);
      }
      else if (kind == 4u)
      {
        float3 oA, oB, oC;
        const float uNew = bendTerm(term.type, P, pos[0], pos[1], pos[2], &gA, &gB, &gC);
        const float uOld = bendTerm(term.type2, Q, pos[0], pos[1], pos[2], &oA, &oB, &oC);
        gA -= oA;
        gB -= oB;
        gC -= oC;
        acc[2] = lam * (uNew - uOld);
      }
      else if (kind == 5u)
      {
        float3 oA, oB, oC, oD;
        const float uNew = torsionTerm(term.type, P, pos[0], pos[1], pos[2], pos[3], &gA, &gB, &gC, &gD);
        const float uOld = torsionTerm(term.type2, Q, pos[0], pos[1], pos[2], pos[3], &oA, &oB, &oC, &oD);
        gA -= oA;
        gB -= oB;
        gC -= oC;
        gD -= oD;
        acc[3] = lam * (uNew - uOld);
      }
      else if (kind == 6u)
      {
        const float3 dr = pos[0] - pos[1];
        const float rr = dot(dr, dr);
        const float qq = pos4[0].w * pos4[1].w;
        float wVDW, wCoulomb, wExclusion, wFactor, cVDW, cCoulomb, cExclusion, cFactor;
        polymerPairState(term.type, P, rr, qq, p, &wVDW, &wCoulomb, &wExclusion, &wFactor);
        polymerPairState(term.type2, Q, rr, qq, p, &cVDW, &cCoulomb, &cExclusion, &cFactor);
        // the pair lists supply the counted value of a listed pair, nothing of an unlisted one
        const float countedWeight = unlisted ? (1.0f - lam) : -lam;
        acc[((term.flags & 2u) != 0u) ? 7u : 4u] = lam * wVDW + countedWeight * cVDW;
        acc[5] = lam * wCoulomb + countedWeight * cCoulomb;
        acc[6] = lam * wExclusion + countedWeight * cExclusion;
        const float factor = lam * wFactor + countedWeight * cFactor;
        gA = factor * dr;
        gB = -gA;
        scale = 1.0f;
      }
      else
      {
        const float3 bp = pos[1] - pos[0];
        const float3 bn = pos[2] - pos[0];
        const float3 bz = pos[3] - pos[0];
        const float volume = polymerChiralVolume(bp, bn, bz, &gA, &gB, &gC, &gD);
        const float deviation = volume - P[1];
        acc[0] = lam * 0.5f * P[0] * deviation * deviation;
        const float factor = P[0] * deviation;
        gA *= factor;
        gB *= factor;
        gC *= factor;
        gD *= factor;
      }
      gA *= scale;
      gB *= scale;
      gC *= scale;
      gD *= scale;
      if (inter)
      {
        // strain derivative sum_i (r_i - r_0) (x) g_i (the first atom contributes nothing; the gradients of the
        // atoms a term does not have are zero)
        const float3 armB = pos[1] - pos[0];
        const float3 armC = pos[2] - pos[0];
        const float3 armD = pos[3] - pos[0];
        acc[8] += armB.x * gB.x + armC.x * gC.x + armD.x * gD.x;
        acc[9] += armB.x * gB.y + armC.x * gC.y + armD.x * gD.y;
        acc[10] += armB.x * gB.z + armC.x * gC.z + armD.x * gD.z;
        acc[11] += armB.y * gB.x + armC.y * gC.x + armD.y * gD.x;
        acc[12] += armB.y * gB.y + armC.y * gC.y + armD.y * gD.y;
        acc[13] += armB.y * gB.z + armC.y * gC.z + armD.y * gD.z;
        acc[14] += armB.z * gB.x + armC.z * gC.x + armD.z * gD.x;
        acc[15] += armB.z * gB.y + armC.z * gC.y + armD.z * gD.y;
        acc[16] += armB.z * gB.z + armC.z * gC.z + armD.z * gD.z;
      }
    }
    out[0] = FLOAT4(gA, 0.0f);
    out[1] = FLOAT4(gB, 0.0f);
    out[2] = FLOAT4(gC, 0.0f);
    out[3] = FLOAT4(gD, 0.0f);
  }

  for (uint q = 0; q < POLYMER_TERM_PARTIALS; ++q) scratch[q * POLYMER_GROUP + lid] = acc[q];
  reducePartials(scratch, POLYMER_TERM_PARTIALS, partials, lid, POLYMER_GROUP, GROUP_ID());
}

// One work-item per participating atom: gathers the gradients of its terms (the inter-molecular ones first in its
// list), adds the sum to the force of its slot and reduces the virial correction (r - r_com) (x) g of the
// inter-molecular part (ax ay az bx by bz cx cy cz).
KERNEL_GROUP_SIZE(POLYMER_GROUP)
void polymerAtoms(GLOBAL const uint* RESTRICT atoms,           // the participating atoms (system order indices)
                  GLOBAL const uint* RESTRICT refStart,        // CSR offsets of the gradient slots per atom
                  GLOBAL const uint* RESTRICT refLinkEnd,      // per atom: end of its inter-molecular slots
                  GLOBAL const uint* RESTRICT refs,            // gradient slots (4 * term + role)
                  GLOBAL const float4* RESTRICT termGradient,
                  GLOBAL const uint* RESTRICT slotOfAtom,      // per atom (system order): its slot
                  GLOBAL const float4* RESTRICT relative,      // per atom: position relative to the molecule's first atom
                  GLOBAL const uint* RESTRICT moleculeOfAtom,  // per atom: its molecule
                  GLOBAL const float4* RESTRICT moleculeCenter,  // per molecule: mass-weighted center (relative frame)
                  CONSTANT const PolymerParameters* p,
                  GLOBAL float4* RESTRICT force,
                  GLOBAL float* RESTRICT partials KERNEL_INDEX_ARGS)
{
  LOCAL_DECL(float scratch[POLYMER_ATOM_PARTIALS * POLYMER_GROUP]);
  const uint lid = LOCAL_ID();
  const uint g = GLOBAL_ID();
  float acc[POLYMER_ATOM_PARTIALS];
  for (uint q = 0; q < POLYMER_ATOM_PARTIALS; ++q) acc[q] = 0.0f;

  if (g < p->numberOfAtoms)
  {
    const uint atom = atoms[g];
    const uint begin = refStart[g];
    const uint linkEnd = refLinkEnd[g];
    const uint end = refStart[g + 1];
    float3 link = FLOAT3(0.0f, 0.0f, 0.0f);
    float3 other = FLOAT3(0.0f, 0.0f, 0.0f);
    for (uint k = begin; k < linkEnd; ++k) link += termGradient[refs[k]].xyz;
    for (uint k = linkEnd; k < end; ++k) other += termGradient[refs[k]].xyz;

    const float3 arm = relative[atom].xyz - moleculeCenter[moleculeOfAtom[atom]].xyz;
    acc[0] = arm.x * link.x;
    acc[1] = arm.x * link.y;
    acc[2] = arm.x * link.z;
    acc[3] = arm.y * link.x;
    acc[4] = arm.y * link.y;
    acc[5] = arm.y * link.z;
    acc[6] = arm.z * link.x;
    acc[7] = arm.z * link.y;
    acc[8] = arm.z * link.z;

    const uint slot = slotOfAtom[atom];
    float4 f = force[slot];
    f.x += link.x + other.x;
    f.y += link.y + other.y;
    f.z += link.z + other.z;
    force[slot] = f;
  }

  for (uint q = 0; q < POLYMER_ATOM_PARTIALS; ++q) scratch[q * POLYMER_GROUP + lid] = acc[q];
  reducePartials(scratch, POLYMER_ATOM_PARTIALS, partials, lid, POLYMER_GROUP, GROUP_ID());
}
)CLC";
