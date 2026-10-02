module;

module spatial_decomposition_device_kernels;

// The resident MD step (device/resident.ixx): the velocity-Verlet integrator of rigid and flexible molecules on
// the device, in the kernel dialect of kernel_sources.ixx. The integration variables (positions, velocities,
// centers of mass, quaternions and their momenta) are double-float numbers: a (hi, lo) pair of floats with ~48
// bits of mantissa ("emulated double"), stored as two float4 arrays. The forces are the single-precision device
// forces of the step. The program must be compiled with DeviceMath::Strict: the error-free transformations below
// need every operation as written (explicit fma, no reassociation).
//
// Per step the host launches, in order:
//   residentAtomsA, residentMoleculesA, residentPack          first half: scale + kick, drift, free rotor,
//                                                              cartesian positions, slot positions + displacements
//   (the device pair / mesh / bonded chain of DeviceStep)
//   residentTorques, residentAtomsB, residentMoleculesB       second half: molecular gradients and torques, kick
//                                                              into the other velocity buffers, kinetic energies
// The thermostat scaling of the host (Nose-Hoover chain) enters as the factors (scaleT, scaleR) of the first
// half; the second half writes the kicked velocities to separate buffers so that the host can repeat it with the
// corrected forces when a list overflowed, without undoing a kick.
const char* const deviceKernelResidentSource = R"KERNEL(
#ifdef __OPENCL_VERSION__
#pragma OPENCL FP_CONTRACT OFF
#endif

#define ATOM_GROUP 256
#define MOLECULE_GROUP 64

// ---------------------------------------------------------------------------------------------------------------
// double-float arithmetic: x = hi + lo with |lo| <= ulp(hi) / 2
// ---------------------------------------------------------------------------------------------------------------
DEVICE_FUNCTION float2 dfTwoSum(float a, float b)
{
  const float s = a + b;
  const float bb = s - a;
  const float e = (a - (s - bb)) + (b - bb);
  return FLOAT2(s, e);
}

DEVICE_FUNCTION float2 dfQuickTwoSum(float a, float b)
{
  const float s = a + b;
  const float e = b - (s - a);
  return FLOAT2(s, e);
}

DEVICE_FUNCTION float2 dfAdd(float2 a, float2 b)
{
  float2 s = dfTwoSum(a.x, b.x);
  const float2 t = dfTwoSum(a.y, b.y);
  s.y += t.x;
  s = dfQuickTwoSum(s.x, s.y);
  s.y += t.y;
  return dfQuickTwoSum(s.x, s.y);
}

DEVICE_FUNCTION float2 dfNeg(float2 a) { return FLOAT2(-a.x, -a.y); }

DEVICE_FUNCTION float2 dfSub(float2 a, float2 b) { return dfAdd(a, dfNeg(b)); }

DEVICE_FUNCTION float2 dfMul(float2 a, float2 b)
{
  const float p = a.x * b.x;
  float e = fma(a.x, b.x, -p);
  e = fma(a.x, b.y, e);
  e = fma(a.y, b.x, e);
  return dfQuickTwoSum(p, e);
}

DEVICE_FUNCTION float2 dfMulF(float2 a, float b)
{
  const float p = a.x * b;
  float e = fma(a.x, b, -p);
  e = fma(a.y, b, e);
  return dfQuickTwoSum(p, e);
}

DEVICE_FUNCTION float2 dfHalf(float2 a) { return FLOAT2(0.5f * a.x, 0.5f * a.y); }

DEVICE_FUNCTION float2 dfOf(float a) { return FLOAT2(a, 0.0f); }

DEVICE_FUNCTION float dfValue(float2 a) { return a.x + a.y; }

// sin and cos of a double-float: reduction by multiples of pi / 2, Taylor series on |r| <= pi / 4 (the terms
// beyond r^19 / 19! are below 1e-19)
DEVICE_FUNCTION void dfSinCos(float2 x, PRIVATE float2* s, PRIVATE float2* c)
{
  const float2 halfPi = FLOAT2(1.5707963705062866f, -4.371138828673793e-08f);
  const float k = rint(x.x * 0.6366197466850281f);
  const float2 r = dfSub(x, dfMulF(halfPi, k));
  const float2 r2 = dfMul(r, r);

  float2 p = FLOAT2(-8.220635078476521e-18f, -1.681478096855898e-25f);
  p = dfAdd(FLOAT2(2.8114573589663704e-15f, -1.0462084739763658e-22f), dfMul(r2, p));
  p = dfAdd(FLOAT2(-7.647163609812713e-13f, -1.2200710471178288e-20f), dfMul(r2, p));
  p = dfAdd(FLOAT2(1.6059044372074283e-10f, -5.352526511562726e-18f), dfMul(r2, p));
  p = dfAdd(FLOAT2(-2.5052107943679403e-08f, -4.4176230446483665e-16f), dfMul(r2, p));
  p = dfAdd(FLOAT2(2.7557318844628753e-06f, 3.793571224297229e-14f), dfMul(r2, p));
  p = dfAdd(FLOAT2(-0.00019841270113829523f, 2.725596874933456e-12f), dfMul(r2, p));
  p = dfAdd(FLOAT2(0.008333333767950535f, -4.34617203337595e-10f), dfMul(r2, p));
  p = dfAdd(FLOAT2(-0.1666666716337204f, 4.967053879312289e-09f), dfMul(r2, p));
  const float2 sinR = dfAdd(r, dfMul(r, dfMul(r2, p)));

  float2 q = FLOAT2(-1.5619206814541513e-16f, -1.5404471465941993e-24f);
  q = dfAdd(FLOAT2(4.7794772561329454e-14f, 7.62544404448643e-22f), dfMul(r2, q));
  q = dfAdd(FLOAT2(-1.147074536050896e-11f, -2.372207689231238e-19f), dfMul(r2, q));
  q = dfAdd(FLOAT2(2.0876755879584152e-09f, 1.1082839809204342e-16f), dfMul(r2, q));
  q = dfAdd(FLOAT2(-2.755731998149713e-07f, 7.575112209051195e-15f), dfMul(r2, q));
  q = dfAdd(FLOAT2(2.4801587642286904e-05f, -3.40699609366682e-13f), dfMul(r2, q));
  q = dfAdd(FLOAT2(-0.0013888889225199819f, 3.3631094437103215e-11f), dfMul(r2, q));
  q = dfAdd(FLOAT2(0.0416666679084301f, -1.2417634698280722e-09f), dfMul(r2, q));
  q = dfAdd(FLOAT2(-0.5f, 0.0f), dfMul(r2, q));
  const float2 cosR = dfAdd(dfOf(1.0f), dfMul(r2, q));

  const int quadrant = ((int)k) & 3;
  if (quadrant == 0)
  {
    *s = sinR;
    *c = cosR;
  }
  else if (quadrant == 1)
  {
    *s = cosR;
    *c = dfNeg(sinR);
  }
  else if (quadrant == 2)
  {
    *s = dfNeg(sinR);
    *c = dfNeg(cosR);
  }
  else
  {
    *s = dfNeg(cosR);
    *c = sinR;
  }
}

typedef struct
{
  float2 x, y, z;
} DF3;

// quaternion (r; ix, iy, iz); stored as float4 (ix, iy, iz, r)
typedef struct
{
  float2 r, x, y, z;
} DF4;

DEVICE_FUNCTION DF3 df3Load(GLOBAL const float4* hi, GLOBAL const float4* lo, uint i)
{
  const float4 h = hi[i];
  const float4 l = lo[i];
  DF3 v;
  v.x = FLOAT2(h.x, l.x);
  v.y = FLOAT2(h.y, l.y);
  v.z = FLOAT2(h.z, l.z);
  return v;
}

DEVICE_FUNCTION void df3Store(GLOBAL float4* hi, GLOBAL float4* lo, uint i, DF3 v, float w)
{
  hi[i] = FLOAT4(v.x.x, v.y.x, v.z.x, w);
  lo[i] = FLOAT4(v.x.y, v.y.y, v.z.y, 0.0f);
}

DEVICE_FUNCTION DF4 df4Load(GLOBAL const float4* hi, GLOBAL const float4* lo, uint i)
{
  const float4 h = hi[i];
  const float4 l = lo[i];
  DF4 q;
  q.x = FLOAT2(h.x, l.x);
  q.y = FLOAT2(h.y, l.y);
  q.z = FLOAT2(h.z, l.z);
  q.r = FLOAT2(h.w, l.w);
  return q;
}

DEVICE_FUNCTION void df4Store(GLOBAL float4* hi, GLOBAL float4* lo, uint i, DF4 q)
{
  hi[i] = FLOAT4(q.x.x, q.y.x, q.z.x, q.r.x);
  lo[i] = FLOAT4(q.x.y, q.y.y, q.z.y, q.r.y);
}

DEVICE_FUNCTION float3 df3Value(DF3 v) { return FLOAT3(dfValue(v.x), dfValue(v.y), dfValue(v.z)); }

DEVICE_FUNCTION DF3 df3Scale(DF3 v, float2 s)
{
  DF3 r;
  r.x = dfMul(v.x, s);
  r.y = dfMul(v.y, s);
  r.z = dfMul(v.z, s);
  return r;
}

DEVICE_FUNCTION DF3 df3Add(DF3 a, DF3 b)
{
  DF3 r;
  r.x = dfAdd(a.x, b.x);
  r.y = dfAdd(a.y, b.y);
  r.z = dfAdd(a.z, b.z);
  return r;
}

DEVICE_FUNCTION DF3 df3Sub(DF3 a, DF3 b)
{
  DF3 r;
  r.x = dfSub(a.x, b.x);
  r.y = dfSub(a.y, b.y);
  r.z = dfSub(a.z, b.z);
  return r;
}

DEVICE_FUNCTION float2 df3Dot(DF3 a, DF3 b)
{
  return dfAdd(dfAdd(dfMul(a.x, b.x), dfMul(a.y, b.y)), dfMul(a.z, b.z));
}

DEVICE_FUNCTION DF3 df3Cross(DF3 a, DF3 b)
{
  DF3 r;
  r.x = dfSub(dfMul(a.y, b.z), dfMul(a.z, b.y));
  r.y = dfSub(dfMul(a.z, b.x), dfMul(a.x, b.z));
  r.z = dfSub(dfMul(a.x, b.y), dfMul(a.y, b.x));
  return r;
}

// v -= factor * g (g single precision, factor and v double-float)
DEVICE_FUNCTION DF3 df3KickScaled(DF3 v, float2 scale, float2 factor, float3 g)
{
  DF3 r;
  r.x = dfSub(dfMul(v.x, scale), dfMulF(factor, g.x));
  r.y = dfSub(dfMul(v.y, scale), dfMulF(factor, g.y));
  r.z = dfSub(dfMul(v.z, scale), dfMulF(factor, g.z));
  return r;
}

DEVICE_FUNCTION DF4 df4KickScaled(DF4 p, float2 scale, float2 factor, float4 t)
{
  DF4 r;
  r.x = dfSub(dfMul(p.x, scale), dfMulF(factor, t.x));
  r.y = dfSub(dfMul(p.y, scale), dfMulF(factor, t.y));
  r.z = dfSub(dfMul(p.z, scale), dfMulF(factor, t.z));
  r.r = dfSub(dfMul(p.r, scale), dfMulF(factor, t.w));
  return r;
}

// rotation of the body-fixed vector v into the laboratory frame by the unit quaternion q (simd_quatd operator*):
// 2 (u . v) u + (r^2 - u . u) v + 2 r (u x v)
DEVICE_FUNCTION DF3 df3Rotate(DF4 q, DF3 v)
{
  DF3 u;
  u.x = q.x;
  u.y = q.y;
  u.z = q.z;
  const float2 uv2 = dfMulF(df3Dot(u, v), 2.0f);
  const float2 w = dfSub(dfMul(q.r, q.r), df3Dot(u, u));
  const float2 r2 = dfMulF(q.r, 2.0f);
  const DF3 cross = df3Cross(u, v);
  DF3 r;
  r.x = dfAdd(dfAdd(dfMul(uv2, u.x), dfMul(w, v.x)), dfMul(r2, cross.x));
  r.y = dfAdd(dfAdd(dfMul(uv2, u.y), dfMul(w, v.y)), dfMul(r2, cross.y));
  r.z = dfAdd(dfAdd(dfMul(uv2, u.z), dfMul(w, v.z)), dfMul(r2, cross.z));
  return r;
}

// one NO_SQUISH rotation about body axis k (Rigid::NoSquishRotate), zeta = dt (...) / I_k / 4
DEVICE_FUNCTION void noSquishRotate(int k, float2 dt, float2 inverseInertia, PRIVATE DF4* p, PRIVATE DF4* q)
{
  float2 zeta;
  if (k == 1)
  {
    zeta = dfSub(dfAdd(dfSub(dfMul(p->x, q->r), dfMul(p->r, q->x)), dfMul(p->y, q->z)), dfMul(p->z, q->y));
  }
  else if (k == 2)
  {
    zeta = dfAdd(dfAdd(dfSub(dfNeg(dfMul(p->r, q->y)), dfMul(p->x, q->z)), dfMul(p->y, q->r)), dfMul(p->z, q->x));
  }
  else
  {
    zeta = dfAdd(dfSub(dfAdd(dfNeg(dfMul(p->r, q->z)), dfMul(p->x, q->y)), dfMul(p->y, q->x)), dfMul(p->z, q->r));
  }
  zeta = dfMulF(dfMul(dfMul(dt, zeta), inverseInertia), 0.25f);
  float2 s, c;
  dfSinCos(zeta, &s, &c);

  DF4 pn, qn;
  if (k == 1)
  {
    pn.r = dfSub(dfMul(c, p->r), dfMul(s, p->x));
    pn.x = dfAdd(dfMul(c, p->x), dfMul(s, p->r));
    pn.y = dfAdd(dfMul(c, p->y), dfMul(s, p->z));
    pn.z = dfSub(dfMul(c, p->z), dfMul(s, p->y));
    qn.r = dfSub(dfMul(c, q->r), dfMul(s, q->x));
    qn.x = dfAdd(dfMul(c, q->x), dfMul(s, q->r));
    qn.y = dfAdd(dfMul(c, q->y), dfMul(s, q->z));
    qn.z = dfSub(dfMul(c, q->z), dfMul(s, q->y));
  }
  else if (k == 2)
  {
    pn.r = dfSub(dfMul(c, p->r), dfMul(s, p->y));
    pn.x = dfSub(dfMul(c, p->x), dfMul(s, p->z));
    pn.y = dfAdd(dfMul(c, p->y), dfMul(s, p->r));
    pn.z = dfAdd(dfMul(c, p->z), dfMul(s, p->x));
    qn.r = dfSub(dfMul(c, q->r), dfMul(s, q->y));
    qn.x = dfSub(dfMul(c, q->x), dfMul(s, q->z));
    qn.y = dfAdd(dfMul(c, q->y), dfMul(s, q->r));
    qn.z = dfAdd(dfMul(c, q->z), dfMul(s, q->x));
  }
  else
  {
    pn.r = dfSub(dfMul(c, p->r), dfMul(s, p->z));
    pn.x = dfAdd(dfMul(c, p->x), dfMul(s, p->y));
    pn.y = dfSub(dfMul(c, p->y), dfMul(s, p->x));
    pn.z = dfAdd(dfMul(c, p->z), dfMul(s, p->r));
    qn.r = dfSub(dfMul(c, q->r), dfMul(s, q->z));
    qn.x = dfAdd(dfMul(c, q->x), dfMul(s, q->y));
    qn.y = dfSub(dfMul(c, q->y), dfMul(s, q->x));
    qn.z = dfAdd(dfMul(c, q->z), dfMul(s, q->r));
  }
  *p = pn;
  *q = qn;
}

// Rigid::NoSquishFreeRotorOrderTwo: five symmetric sweeps 3, 2, 1, 2, 3 (dtTenth = dt / 10, dtFifth = dt / 5)
DEVICE_FUNCTION void noSquishFreeRotor(float2 dtTenth, float2 dtFifth, float2 invIx, float2 invIy, float2 invIz,
                                       PRIVATE DF4* p, PRIVATE DF4* q)
{
  for (int i = 0; i < 5; ++i)
  {
    noSquishRotate(3, dtTenth, invIz, p, q);
    noSquishRotate(2, dtTenth, invIy, p, q);
    noSquishRotate(1, dtFifth, invIx, p, q);
    noSquishRotate(2, dtTenth, invIy, p, q);
    noSquishRotate(3, dtTenth, invIz, p, q);
  }
}

// rotational kinetic energy of (p, q): omega = (-p.r, p.u) q / 2 / I, E = sum I omega^2 / 2
DEVICE_FUNCTION float2 rotationalEnergy(DF4 p, DF4 q, float4 inertiaHi, float4 inertiaLo, float4 inverseHi,
                                        float4 inverseLo)
{
  const float2 pr = dfNeg(p.r);
  const float2 wx = dfSub(dfAdd(dfAdd(dfMul(pr, q.x), dfMul(p.x, q.r)), dfMul(p.y, q.z)), dfMul(p.z, q.y));
  const float2 wy = dfAdd(dfAdd(dfSub(dfMul(pr, q.y), dfMul(p.x, q.z)), dfMul(p.y, q.r)), dfMul(p.z, q.x));
  const float2 wz = dfAdd(dfSub(dfAdd(dfMul(pr, q.z), dfMul(p.x, q.y)), dfMul(p.y, q.x)), dfMul(p.z, q.r));
  const float2 ox = dfMul(dfHalf(wx), FLOAT2(inverseHi.x, inverseLo.x));
  const float2 oy = dfMul(dfHalf(wy), FLOAT2(inverseHi.y, inverseLo.y));
  const float2 oz = dfMul(dfHalf(wz), FLOAT2(inverseHi.z, inverseLo.z));
  float2 e = dfMul(FLOAT2(inertiaHi.x, inertiaLo.x), dfMul(ox, ox));
  e = dfAdd(e, dfMul(FLOAT2(inertiaHi.y, inertiaLo.y), dfMul(oy, oy)));
  e = dfAdd(e, dfMul(FLOAT2(inertiaHi.z, inertiaLo.z), dfMul(oz, oz)));
  return dfHalf(e);
}

// ---------------------------------------------------------------------------------------------------------------
// work-group reductions
// ---------------------------------------------------------------------------------------------------------------
DEVICE_FUNCTION void reduceSum(LOCAL float2* scratch, float2 value, uint lid, uint size, GLOBAL float2* out,
                               uint group)
{
  scratch[lid] = value;
  LOCAL_BARRIER();
  for (uint s = size / 2; s > 0; s >>= 1)
  {
    if (lid < s) scratch[lid] = dfAdd(scratch[lid], scratch[lid + s]);
    LOCAL_BARRIER();
  }
  if (lid == 0) out[group] = scratch[0];
}

DEVICE_FUNCTION void reduceMax(LOCAL float2* scratch, float2 value, uint lid, uint size, GLOBAL float2* out,
                               uint group)
{
  scratch[lid] = value;
  LOCAL_BARRIER();
  for (uint s = size / 2; s > 0; s >>= 1)
  {
    if (lid < s) scratch[lid] = fmax(scratch[lid], scratch[lid + s]);
    LOCAL_BARRIER();
  }
  if (lid == 0) out[group] = scratch[0];
}

// ---------------------------------------------------------------------------------------------------------------
// first half of the step
// ---------------------------------------------------------------------------------------------------------------

// flexible atoms: v = scaleT v - dt/2 g / m, r += dt v (the atoms of rigid molecules are set by the molecules)
KERNEL_GROUP_SIZE(ATOM_GROUP)
void residentAtomsA(GLOBAL float4* RESTRICT positionHi, GLOBAL float4* RESTRICT positionLo,
                    GLOBAL float4* RESTRICT velocityHi, GLOBAL float4* RESTRICT velocityLo,
                    GLOBAL const float4* RESTRICT atomMass, GLOBAL const uint4* RESTRICT atomInfo,
                    GLOBAL const uint* RESTRICT slotOfAtom, GLOBAL const float4* RESTRICT force,
                    VALUE_ARG(uint, numberOfAtoms), VALUE_ARG(float2, scaleT), VALUE_ARG(float2, halfDt),
                    VALUE_ARG(float2, dt) KERNEL_INDEX_ARGS)
{
  const uint i = GLOBAL_ID();
  if (i >= numberOfAtoms) return;
  const uint4 info = atomInfo[i];
  if (info.w == 0u) return;
  const float4 g = force[slotOfAtom[i]];
  const float4 m = atomMass[i];
  const float2 invMass = FLOAT2(m.z, m.w);
  const float2 factor = dfMul(halfDt, invMass);
  DF3 v = df3Load(velocityHi, velocityLo, i);
  v = df3KickScaled(v, scaleT, factor, g.xyz);
  const float charge = positionHi[i].w;
  DF3 r = df3Load(positionHi, positionLo, i);
  r = df3Add(r, df3Scale(v, dt));
  df3Store(velocityHi, velocityLo, i, v, 0.0f);
  df3Store(positionHi, positionLo, i, r, charge);
}

// molecules: V = scaleT V - dt/2 G / M, COM += dt V; rigid: P = scaleR P - dt/2 T, free rotor, atom positions
// COM + q ref; flexible: COM = mass-weighted mean of the (drifted) atoms
KERNEL_GROUP_SIZE(MOLECULE_GROUP)
void residentMoleculesA(GLOBAL float4* RESTRICT positionHi, GLOBAL float4* RESTRICT positionLo,
                        GLOBAL const float4* RESTRICT atomMass, GLOBAL float4* RESTRICT comHi,
                        GLOBAL float4* RESTRICT comLo, GLOBAL float4* RESTRICT velocityHi,
                        GLOBAL float4* RESTRICT velocityLo, GLOBAL float4* RESTRICT orientationHi,
                        GLOBAL float4* RESTRICT orientationLo, GLOBAL float4* RESTRICT momentumHi,
                        GLOBAL float4* RESTRICT momentumLo, GLOBAL const float4* RESTRICT moleculeGradient,
                        GLOBAL const float4* RESTRICT moleculeTorque, GLOBAL const float4* RESTRICT moleculeMass,
                        GLOBAL const uint4* RESTRICT moleculeInfo, GLOBAL const float4* RESTRICT componentInertia,
                        GLOBAL const uint* RESTRICT componentReferenceOffset,
                        GLOBAL const float4* RESTRICT referenceHi, GLOBAL const float4* RESTRICT referenceLo,
                        VALUE_ARG(uint, numberOfMolecules), VALUE_ARG(float2, scaleT), VALUE_ARG(float2, scaleR),
                        VALUE_ARG(float2, halfDt), VALUE_ARG(float2, dt), VALUE_ARG(float2, dtTenth),
                        VALUE_ARG(float2, dtFifth) KERNEL_INDEX_ARGS)
{
  const uint m = GLOBAL_ID();
  if (m >= numberOfMolecules) return;
  const uint4 info = moleculeInfo[m];
  const float4 mass = moleculeMass[m];
  const float2 invMass = FLOAT2(mass.z, mass.w);
  const float4 g = moleculeGradient[m];

  DF3 velocity = df3Load(velocityHi, velocityLo, m);
  velocity = df3KickScaled(velocity, scaleT, dfMul(halfDt, invMass), g.xyz);
  DF3 com = df3Load(comHi, comLo, m);
  com = df3Add(com, df3Scale(velocity, dt));
  df3Store(velocityHi, velocityLo, m, velocity, 0.0f);

  if (info.w != 0u)
  {
    DF4 p = df4Load(momentumHi, momentumLo, m);
    DF4 q = df4Load(orientationHi, orientationLo, m);
    p = df4KickScaled(p, scaleR, halfDt, moleculeTorque[m]);
    const uint c = info.z;
    const float4 inverseHi = componentInertia[4u * c + 2u];
    const float4 inverseLo = componentInertia[4u * c + 3u];
    noSquishFreeRotor(dtTenth, dtFifth, FLOAT2(inverseHi.x, inverseLo.x), FLOAT2(inverseHi.y, inverseLo.y),
                      FLOAT2(inverseHi.z, inverseLo.z), &p, &q);
    df4Store(momentumHi, momentumLo, m, p);
    df4Store(orientationHi, orientationLo, m, q);

    const uint offset = componentReferenceOffset[c];
    for (uint b = 0; b < info.y; ++b)
    {
      const uint atom = info.x + b;
      const DF3 reference = df3Load(referenceHi, referenceLo, offset + b);
      const DF3 position = df3Add(com, df3Rotate(q, reference));
      df3Store(positionHi, positionLo, atom, position, positionHi[atom].w);
    }
  }
  else
  {
    DF3 weighted;
    weighted.x = dfOf(0.0f);
    weighted.y = dfOf(0.0f);
    weighted.z = dfOf(0.0f);
    for (uint b = 0; b < info.y; ++b)
    {
      const uint atom = info.x + b;
      const float4 am = atomMass[atom];
      weighted = df3Add(weighted, df3Scale(df3Load(positionHi, positionLo, atom), FLOAT2(am.x, am.y)));
    }
    com = df3Scale(weighted, invMass);
  }
  df3Store(comHi, comLo, m, com, 0.0f);
}

// slot positions (wrapped: position minus the box translation of the binning) and charges for the pair, mesh and
// bonded kernels, the positions relative to the first atom of the molecule for the bonded kernel, and per group
// the largest squared displacement since the list build and since the list compaction
KERNEL_GROUP_SIZE(ATOM_GROUP)
void residentPack(GLOBAL const float4* RESTRICT positionHi, GLOBAL const float4* RESTRICT positionLo,
                  GLOBAL const float4* RESTRICT translationHi, GLOBAL const float4* RESTRICT translationLo,
                  GLOBAL const uint* RESTRICT slotOfAtom, GLOBAL const uint4* RESTRICT atomInfo,
                  GLOBAL const float4* RESTRICT buildPosition, GLOBAL const float4* RESTRICT compactReference,
                  GLOBAL float4* RESTRICT slotPosition, GLOBAL float4* RESTRICT relative,
                  GLOBAL float2* RESTRICT partials, VALUE_ARG(uint, numberOfAtoms),
                  VALUE_ARG(uint, writeRelative) KERNEL_INDEX_ARGS)
{
  LOCAL float2 scratch[ATOM_GROUP];
  const uint i = GLOBAL_ID();
  float2 displacement = FLOAT2(0.0f, 0.0f);
  if (i < numberOfAtoms)
  {
    const DF3 r = df3Load(positionHi, positionLo, i);
    const DF3 t = df3Load(translationHi, translationLo, i);
    const float3 p = df3Value(df3Sub(r, t));
    const uint slot = slotOfAtom[i];
    const float charge = positionHi[i].w;
    slotPosition[slot] = FLOAT4(p.x, p.y, p.z, charge);
    const float3 db = p - buildPosition[slot].xyz;
    const float3 dc = p - compactReference[slot].xyz;
    displacement = FLOAT2(dot(db, db), dot(dc, dc));
    if (writeRelative != 0u)
    {
      const DF3 first = df3Load(positionHi, positionLo, atomInfo[i].z);
      const float3 d = df3Value(df3Sub(r, first));
      relative[i] = FLOAT4(d.x, d.y, d.z, charge);
    }
  }
  reduceMax(scratch, displacement, LOCAL_ID(), ATOM_GROUP, partials, GROUP_ID());
}

// ---------------------------------------------------------------------------------------------------------------
// second half of the step
// ---------------------------------------------------------------------------------------------------------------

// per molecule the gradient G = sum g and, for rigid molecules, the orientation gradient -2 q (0, torque) with
// torque = sum (R(q) (g_b - G m_b / M)) x ref_b (single precision, like the forces)
KERNEL_GROUP_SIZE(MOLECULE_GROUP)
void residentTorques(GLOBAL const float4* RESTRICT force, GLOBAL const uint* RESTRICT slotOfAtom,
                     GLOBAL const float4* RESTRICT atomMass, GLOBAL const float4* RESTRICT orientationHi,
                     GLOBAL const float4* RESTRICT moleculeMass, GLOBAL const uint4* RESTRICT moleculeInfo,
                     GLOBAL const uint* RESTRICT componentReferenceOffset,
                     GLOBAL const float4* RESTRICT referenceHi, GLOBAL float4* RESTRICT moleculeGradient,
                     GLOBAL float4* RESTRICT moleculeTorque, VALUE_ARG(uint, numberOfMolecules) KERNEL_INDEX_ARGS)
{
  const uint m = GLOBAL_ID();
  if (m >= numberOfMolecules) return;
  const uint4 info = moleculeInfo[m];
  float3 total = FLOAT3(0.0f, 0.0f, 0.0f);
  for (uint b = 0; b < info.y; ++b) total += force[slotOfAtom[info.x + b]].xyz;
  moleculeGradient[m] = FLOAT4(total.x, total.y, total.z, 0.0f);
  if (info.w == 0u)
  {
    moleculeTorque[m] = FLOAT4(0.0f, 0.0f, 0.0f, 0.0f);
    return;
  }
  const float4 q = orientationHi[m];
  const float r = q.w, ix = q.x, iy = q.y, iz = q.z;
  // rows of R(q) (double3x3::buildRotationMatrix, laboratory -> body frame)
  const float3 row0 = FLOAT3(2.0f * (r * r + ix * ix) - 1.0f, 2.0f * (ix * iy + r * iz), 2.0f * (ix * iz - r * iy));
  const float3 row1 = FLOAT3(2.0f * (ix * iy - r * iz), 2.0f * (iy * iy + r * r) - 1.0f, 2.0f * (iy * iz + r * ix));
  const float3 row2 = FLOAT3(2.0f * (ix * iz + r * iy), 2.0f * (iy * iz - r * ix), 2.0f * (r * r + iz * iz) - 1.0f);
  const float invMass = moleculeMass[m].z;
  const uint offset = componentReferenceOffset[info.z];
  float3 torque = FLOAT3(0.0f, 0.0f, 0.0f);
  for (uint b = 0; b < info.y; ++b)
  {
    const uint atom = info.x + b;
    const float3 f = force[slotOfAtom[atom]].xyz - total * (atomMass[atom].x * invMass);
    const float3 body = FLOAT3(dot(row0, f), dot(row1, f), dot(row2, f));
    torque += cross(body, referenceHi[offset + b].xyz);
  }
  // -2 q (0, torque)
  const float3 u = FLOAT3(ix, iy, iz);
  const float3 uxt = cross(u, torque);
  moleculeTorque[m] = -2.0f * FLOAT4(r * torque.x + uxt.x, r * torque.y + uxt.y, r * torque.z + uxt.z, -dot(u, torque));
}

// flexible atoms: v_out = v_in - dt/2 g / m; per group the translational kinetic energy sum m v^2 / 2
KERNEL_GROUP_SIZE(ATOM_GROUP)
void residentAtomsB(GLOBAL const float4* RESTRICT velocityInHi, GLOBAL const float4* RESTRICT velocityInLo,
                    GLOBAL float4* RESTRICT velocityOutHi, GLOBAL float4* RESTRICT velocityOutLo,
                    GLOBAL const float4* RESTRICT atomMass, GLOBAL const uint4* RESTRICT atomInfo,
                    GLOBAL const uint* RESTRICT slotOfAtom, GLOBAL const float4* RESTRICT force,
                    GLOBAL float2* RESTRICT partials, VALUE_ARG(uint, numberOfAtoms),
                    VALUE_ARG(float2, halfDt) KERNEL_INDEX_ARGS)
{
  LOCAL float2 scratch[ATOM_GROUP];
  const uint i = GLOBAL_ID();
  float2 kinetic = FLOAT2(0.0f, 0.0f);
  if (i < numberOfAtoms && atomInfo[i].w != 0u)
  {
    const float4 g = force[slotOfAtom[i]];
    const float4 m = atomMass[i];
    DF3 v = df3Load(velocityInHi, velocityInLo, i);
    v = df3KickScaled(v, dfOf(1.0f), dfMul(halfDt, FLOAT2(m.z, m.w)), g.xyz);
    df3Store(velocityOutHi, velocityOutLo, i, v, 0.0f);
    kinetic = dfHalf(dfMul(FLOAT2(m.x, m.y), df3Dot(v, v)));
  }
  reduceSum(scratch, kinetic, LOCAL_ID(), ATOM_GROUP, partials, GROUP_ID());
}

// molecules: V_out = V_in - dt/2 G / M; rigid: P_out = P_in - dt/2 T; per group the translational kinetic energy
// of the rigid molecules (M V^2 / 2) and the rotational kinetic energy (partials[2 group], partials[2 group + 1])
KERNEL_GROUP_SIZE(MOLECULE_GROUP)
void residentMoleculesB(GLOBAL const float4* RESTRICT velocityInHi, GLOBAL const float4* RESTRICT velocityInLo,
                        GLOBAL float4* RESTRICT velocityOutHi, GLOBAL float4* RESTRICT velocityOutLo,
                        GLOBAL const float4* RESTRICT momentumInHi, GLOBAL const float4* RESTRICT momentumInLo,
                        GLOBAL float4* RESTRICT momentumOutHi, GLOBAL float4* RESTRICT momentumOutLo,
                        GLOBAL const float4* RESTRICT orientationHi, GLOBAL const float4* RESTRICT orientationLo,
                        GLOBAL const float4* RESTRICT moleculeGradient,
                        GLOBAL const float4* RESTRICT moleculeTorque, GLOBAL const float4* RESTRICT moleculeMass,
                        GLOBAL const uint4* RESTRICT moleculeInfo, GLOBAL const float4* RESTRICT componentInertia,
                        GLOBAL float2* RESTRICT partials, VALUE_ARG(uint, numberOfMolecules),
                        VALUE_ARG(float2, halfDt) KERNEL_INDEX_ARGS)
{
  LOCAL float2 scratchT[MOLECULE_GROUP];
  LOCAL float2 scratchR[MOLECULE_GROUP];
  const uint m = GLOBAL_ID();
  float2 kineticT = FLOAT2(0.0f, 0.0f);
  float2 kineticR = FLOAT2(0.0f, 0.0f);
  if (m < numberOfMolecules)
  {
    const uint4 info = moleculeInfo[m];
    const float4 mass = moleculeMass[m];
    DF3 velocity = df3Load(velocityInHi, velocityInLo, m);
    velocity = df3KickScaled(velocity, dfOf(1.0f), dfMul(halfDt, FLOAT2(mass.z, mass.w)), moleculeGradient[m].xyz);
    df3Store(velocityOutHi, velocityOutLo, m, velocity, 0.0f);
    DF4 p = df4Load(momentumInHi, momentumInLo, m);
    if (info.w != 0u)
    {
      p = df4KickScaled(p, dfOf(1.0f), halfDt, moleculeTorque[m]);
      kineticT = dfHalf(dfMul(FLOAT2(mass.x, mass.y), df3Dot(velocity, velocity)));
      const DF4 q = df4Load(orientationHi, orientationLo, m);
      const uint c = info.z;
      kineticR = rotationalEnergy(p, q, componentInertia[4u * c], componentInertia[4u * c + 1u],
                                  componentInertia[4u * c + 2u], componentInertia[4u * c + 3u]);
    }
    df4Store(momentumOutHi, momentumOutLo, m, p);
  }
  reduceSum(scratchT, kineticT, LOCAL_ID(), MOLECULE_GROUP, partials, 2u * GROUP_ID());
  reduceSum(scratchR, kineticR, LOCAL_ID(), MOLECULE_GROUP, partials, 2u * GROUP_ID() + 1u);
}

// ---------------------------------------------------------------------------------------------------------------
// velocity scaling (the deferred thermostat factor before a download of the state)
// ---------------------------------------------------------------------------------------------------------------
KERNEL_GROUP_SIZE(ATOM_GROUP)
void residentScaleAtoms(GLOBAL float4* RESTRICT velocityHi, GLOBAL float4* RESTRICT velocityLo,
                        GLOBAL const uint4* RESTRICT atomInfo, VALUE_ARG(uint, numberOfAtoms),
                        VALUE_ARG(float2, scaleT) KERNEL_INDEX_ARGS)
{
  const uint i = GLOBAL_ID();
  if (i >= numberOfAtoms || atomInfo[i].w == 0u) return;
  const DF3 v = df3Scale(df3Load(velocityHi, velocityLo, i), scaleT);
  df3Store(velocityHi, velocityLo, i, v, 0.0f);
}

KERNEL_GROUP_SIZE(MOLECULE_GROUP)
void residentScaleMolecules(GLOBAL float4* RESTRICT velocityHi, GLOBAL float4* RESTRICT velocityLo,
                            GLOBAL float4* RESTRICT momentumHi, GLOBAL float4* RESTRICT momentumLo,
                            GLOBAL const uint4* RESTRICT moleculeInfo, VALUE_ARG(uint, numberOfMolecules),
                            VALUE_ARG(float2, scaleT), VALUE_ARG(float2, scaleR) KERNEL_INDEX_ARGS)
{
  const uint m = GLOBAL_ID();
  if (m >= numberOfMolecules) return;
  const DF3 v = df3Scale(df3Load(velocityHi, velocityLo, m), scaleT);
  df3Store(velocityHi, velocityLo, m, v, 0.0f);
  if (moleculeInfo[m].w != 0u)
  {
    DF4 p = df4Load(momentumHi, momentumLo, m);
    p.r = dfMul(p.r, scaleR);
    p.x = dfMul(p.x, scaleR);
    p.y = dfMul(p.y, scaleR);
    p.z = dfMul(p.z, scaleR);
    df4Store(momentumHi, momentumLo, m, p);
  }
}
)KERNEL";
