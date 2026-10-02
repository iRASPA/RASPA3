module;

module spatial_decomposition_opencl_mesh;

// OpenCL C (1.2) source of the particle-mesh Ewald sum on the device: the device transcription of PPPM
// (pppm.cpp) over the slots of the pair kernel.
//
//  - spreadCharges: cardinal B-spline charge assignment of every charged slot onto the mesh. The device has no
//    floating-point atomics (OpenCL 1.2), so the mesh is accumulated in 32-bit fixed point with integer atomics:
//    order independent, hence bit-reproducible.
//  - The 3D transform is real-to-complex: the charge mesh is real, so only the half spectrum kz = 0..Kz/2 is
//    computed and stored ((ix Ky + iy) Hz + kz with Hz = Kz/2 + 1). All transforms are unnormalized (like FFTW) and
//    use one Stockham autosort FFT (radices 2, 3, 4, 5, twiddles from a table) per line in local memory; a
//    work-group handles a tile of adjacent lines so that the global loads coalesce on the strided axes.
//  - fftRealForward: the z lines. Reads the fixed-point charge mesh (converting it to single precision and zeroing
//    it for the next step), packs the even and odd samples into a complex line of length M = Kz/2, transforms it
//    and unpacks the half spectrum X[k] = E[k] + W^k O[k] (E, O: the spectra of the even and odd samples).
//  - fftLines: the complex y and x lines of the half spectrum, forward and backward by a direction flag.
//  - fftRealBackward: the inverse of fftRealForward, writing the real potential mesh.
//  - applyInfluence: multiplies the half spectrum with G(m) = C 2 pi / V exp(-k^2 / 4 alpha^2) / k^2 |b|^2 computed
//    on the fly (so a cell change costs nothing extra) and reduces per work-group the reciprocal energy, its
//    strain derivative and the single-ion sums of the net-charge correction (the interior kz planes count twice
//    for their conjugate partners).
//  - interpolateForces: the gradient dE/dr of every charged slot from the potential mesh with the analytic
//    B-spline derivatives, added to the force of the slot.
//
// The B-spline order is a compile-time constant (MESH_ORDER, set when the program is built) so that the weight
// arrays live in registers.
const char* const openclMeshKernelSource = R"CLC(
#define NO_ATOM 0xFFFFFFFFu
#define TWO_PI 6.283185307179586f
#define INFLUENCE_GROUP 128
#define INFLUENCE_PARTIALS 14
#define INFLUENCE_POINTS 8

typedef struct
{
  float inverseCell[9];  // ax ay az bx by bz cx cy cz (columns) of the inverse cell matrix
  float scale;           // fixed-point scale of the charge mesh
  float inverseScale;
  float prefactor;       // C 2 pi / V
  float alphaFactor;     // -1 / (4 alpha^2)
  float inverseFourAlphaSquared;
  uint meshX, meshY, meshZ;
  uint numberOfSlots;
  uint padding[6];
} MeshParameters;

inline float3 fractional(float3 r, __constant const float* ic)
{
  float3 s;
  s.x = ic[0] * r.x + ic[3] * r.y + ic[6] * r.z;
  s.y = ic[1] * r.x + ic[4] * r.y + ic[7] * r.z;
  s.z = ic[2] * r.x + ic[5] * r.y + ic[8] * r.z;
  return s - floor(s);
}

// M_p(w + j), j = 0..p-1, and the derivatives M_{p-1}(w + j) - M_{p-1}(w + j - 1) (PPPM::bsplineWeights)
inline void bsplineWeights(float w, float* weights, float* derivatives)
{
#pragma unroll
  for (uint j = 0; j < MESH_ORDER; ++j)
  {
    weights[j] = 0.0f;
    derivatives[j] = 0.0f;
  }
  weights[0] = w;
  weights[1] = 1.0f - w;
#pragma unroll
  for (uint k = 3; k <= MESH_ORDER; ++k)
  {
    if (k == MESH_ORDER)
    {
      derivatives[0] = weights[0];
#pragma unroll
      for (uint j = 1; j < MESH_ORDER; ++j) derivatives[j] = weights[j] - weights[j - 1];
    }
    const float inverse = 1.0f / (float)(k - 1);
#pragma unroll
    for (uint j = k - 1; j > 0; --j)
    {
      const float u = w + (float)j;
      weights[j] = (u * weights[j] + ((float)k - u) * weights[j - 1]) * inverse;
    }
    weights[0] = w * weights[0] * inverse;
  }
}

// anchor (mesh point at or below the atom) and the fractional offset along one axis
inline int anchor(float s, uint K, float* w)
{
  const float u = min(s * (float)K, nextafter((float)K, 0.0f));
  const int k0 = min((int)u, (int)K - 1);
  *w = u - (float)k0;
  return k0;
}

__kernel void spreadCharges(__global const float4* restrict position,  // x, y, z, charge per slot
                            __constant const MeshParameters* p,
                            __global int* restrict mesh)
{
  const uint slot = get_global_id(0);
  if (slot >= p->numberOfSlots) return;
  const float4 r = position[slot];
  const float q = r.w;
  if (q == 0.0f) return;
  const float3 s = fractional(r.xyz, p->inverseCell);
  float wx, wy, wz;
  const int kx0 = anchor(s.x, p->meshX, &wx);
  const int ky0 = anchor(s.y, p->meshY, &wy);
  const int kz0 = anchor(s.z, p->meshZ, &wz);
  float mx[MESH_ORDER], my[MESH_ORDER], mz[MESH_ORDER], dummy[MESH_ORDER];
  bsplineWeights(wx, mx, dummy);
  bsplineWeights(wy, my, dummy);
  bsplineWeights(wz, mz, dummy);
  const int Kx = (int)p->meshX, Ky = (int)p->meshY, Kz = (int)p->meshZ;
  const float scale = q * p->scale;
#pragma unroll
  for (uint a = 0; a < MESH_ORDER; ++a)
  {
    int ix = kx0 - (int)a;
    if (ix < 0) ix += Kx;
    const float qx = scale * mx[a];
#pragma unroll
    for (uint b = 0; b < MESH_ORDER; ++b)
    {
      int iy = ky0 - (int)b;
      if (iy < 0) iy += Ky;
      const float qxy = qx * my[b];
      __global int* row = mesh + ((uint)ix * (uint)Ky + (uint)iy) * (uint)Kz;
#pragma unroll
      for (uint c = 0; c < MESH_ORDER; ++c)
      {
        int iz = kz0 - (int)c;
        if (iz < 0) iz += Kz;
        atomic_add(row + iz, (int)rint(qxy * mz[c]));
      }
    }
  }
}

inline float2 cmul(float2 a, float2 b) { return (float2)(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x); }
inline float2 rotate90(float2 a, float sign) { return (float2)(-sign * a.y, sign * a.x); }  // sign * i * a
inline float2 conjugate(float2 a) { return (float2)(a.x, -a.y); }

// One Stockham stage of radix R over the tile: j-th butterfly of the t-th line (Govindaraju et al., 2008).
// `sign` is -1 for the forward and +1 for the backward transform.
#define STAGE_HEAD(R)                                                                 \
  const uint NR = N / R;                                                              \
  const uint work = NR * tile;                                                        \
  const uint twiddleStride = NR / Ns;                                                 \
  for (uint idx = get_local_id(0); idx < work; idx += get_local_size(0))              \
  {                                                                                   \
    const uint t = idx % tile;                                                        \
    const uint j = idx / tile;                                                        \
    const uint k = j % Ns;                                                            \
    const uint idxD = (j / Ns) * Ns * R + k;                                          \
    float2 v[R];                                                                      \
    for (uint r = 0; r < R; ++r) v[r] = in[(j + r * NR) * tile + t];                  \
    for (uint r = 1; r < R; ++r)                                                      \
    {                                                                                 \
      float2 tw = twiddle[r * k * twiddleStride];                                     \
      tw.y *= -sign;                                                                  \
      v[r] = cmul(v[r], tw);                                                          \
    }

#define STAGE_TAIL(R)                                                       \
    for (uint r = 0; r < R; ++r) out[(idxD + r * Ns) * tile + t] = v[r]; \
  }

inline void stageRadix2(__local const float2* in, __local float2* out, __global const float2* twiddle, uint N,
                        uint Ns, uint tile, float sign)
{
  STAGE_HEAD(2)
  const float2 a = v[0], b = v[1];
  v[0] = a + b;
  v[1] = a - b;
  STAGE_TAIL(2)
}

inline void stageRadix3(__local const float2* in, __local float2* out, __global const float2* twiddle, uint N,
                        uint Ns, uint tile, float sign)
{
  STAGE_HEAD(3)
  const float2 t1 = v[1] + v[2];
  const float2 t2 = v[0] - 0.5f * t1;
  const float2 t3 = rotate90(v[1] - v[2], sign * 0.8660254037844386f);
  v[0] = v[0] + t1;
  v[1] = t2 + t3;
  v[2] = t2 - t3;
  STAGE_TAIL(3)
}

inline void stageRadix4(__local const float2* in, __local float2* out, __global const float2* twiddle, uint N,
                        uint Ns, uint tile, float sign)
{
  STAGE_HEAD(4)
  const float2 a0 = v[0] + v[2], a1 = v[0] - v[2];
  const float2 b0 = v[1] + v[3], b1 = rotate90(v[1] - v[3], sign);
  v[0] = a0 + b0;
  v[1] = a1 + b1;
  v[2] = a0 - b0;
  v[3] = a1 - b1;
  STAGE_TAIL(4)
}

inline void stageRadix5(__local const float2* in, __local float2* out, __global const float2* twiddle, uint N,
                        uint Ns, uint tile, float sign)
{
  STAGE_HEAD(5)
  const float c1 = 0.30901699437494745f, c2 = -0.8090169943749473f;
  const float s1 = sign * 0.9510565162951535f, s2 = sign * 0.5877852522924732f;
  const float2 t1 = v[1] + v[4], t2 = v[2] + v[3], t3 = v[1] - v[4], t4 = v[2] - v[3];
  const float2 a1 = v[0] + c1 * t1 + c2 * t2;
  const float2 a2 = v[0] + c2 * t1 + c1 * t2;
  const float2 b1 = s1 * t3 + s2 * t4;
  const float2 b2 = s2 * t3 - s1 * t4;
  v[0] = v[0] + t1 + t2;
  v[1] = a1 + rotate90(b1, 1.0f);
  v[4] = a1 - rotate90(b1, 1.0f);
  v[2] = a2 + rotate90(b2, 1.0f);
  v[3] = a2 - rotate90(b2, 1.0f);
  STAGE_TAIL(5)
}

// The Stockham stages over the tile held in bufferA (bufferB: scratch); `radixCode` holds the radix of stage s in
// bits 2s..2s+1 (0: 2, 1: 3, 2: 4, 3: 5). Returns the buffer holding the result (bufferA or bufferB).
inline __local float2* runStages(__local float2* bufferA, __local float2* bufferB,
                                 __global const float2* restrict twiddle, uint N, uint radixCode, uint stages,
                                 uint tile, float sign)
{
  __local float2* in = bufferA;
  __local float2* out = bufferB;
  uint Ns = 1;
  for (uint s = 0; s < stages; ++s)
  {
    const uint R = 2 + ((radixCode >> (2 * s)) & 3u);
    switch (R)
    {
      case 2: stageRadix2(in, out, twiddle, N, Ns, tile, sign); break;
      case 3: stageRadix3(in, out, twiddle, N, Ns, tile, sign); break;
      case 4: stageRadix4(in, out, twiddle, N, Ns, tile, sign); break;
      default: stageRadix5(in, out, twiddle, N, Ns, tile, sign); break;
    }
    barrier(CLK_LOCAL_MEM_FENCE);
    __local float2* swap = in;
    in = out;
    out = swap;
    Ns *= R;
  }
  return in;
}

// Forward real transform of the z lines: `lineCount` lines of 2M fixed-point values (line l at l * 2M of the
// charge mesh, zeroed on reading) to the half spectrum of M + 1 complex values at l * (M + 1). The even and odd
// samples are packed into one complex line Z of length M: with E = (Z[k] + conj Z[M-k]) / 2 and
// O = -i (Z[k] - conj Z[M-k]) / 2 the spectrum is X[k] = E + W^k O, W = exp(-2 pi i / 2M) (`halfTwiddle`, k = 0..M).
__kernel void fftRealForward(__global int2* restrict fixedPoint, __global float2* restrict data,
                             __global const float2* restrict twiddle, __global const float2* restrict halfTwiddle,
                             const uint M, const uint radixCode, const uint stages, const uint lineCount,
                             const uint tile, const float inverseScale, __local float2* bufferA,
                             __local float2* bufferB)
{
  const uint lid = get_local_id(0);
  const uint groupSize = get_local_size(0);
  const uint line0 = get_group_id(0) * tile;
  const uint valid = min(tile, lineCount - line0);
  for (uint e = lid; e < M * tile; e += groupSize)
  {
    const uint j = e % M;
    const uint t = e / M;
    float2 value = (float2)(0.0f, 0.0f);
    if (t < valid)
    {
      const uint index = (line0 + t) * M + j;
      const int2 pair = fixedPoint[index];
      value = (float2)((float)pair.x, (float)pair.y) * inverseScale;
      fixedPoint[index] = (int2)(0, 0);
    }
    bufferA[j * tile + t] = value;
  }
  barrier(CLK_LOCAL_MEM_FENCE);
  __local const float2* Z = runStages(bufferA, bufferB, twiddle, M, radixCode, stages, tile, -1.0f);
  const uint H = M + 1;
  for (uint e = lid; e < H * tile; e += groupSize)
  {
    const uint k = e % H;
    const uint t = e / H;
    if (t >= valid) continue;
    const float2 zk = Z[(k == M ? 0u : k) * tile + t];
    const float2 zm = conjugate(Z[(k == 0 ? 0u : M - k) * tile + t]);
    const float2 even = 0.5f * (zk + zm);
    const float2 d = zk - zm;
    const float2 odd = (float2)(0.5f * d.y, -0.5f * d.x);
    data[(line0 + t) * H + k] = even + cmul(halfTwiddle[k], odd);
  }
}

// Inverse of fftRealForward: the half spectrum (M + 1 complex values per line) to 2M real values per line of the
// potential mesh (written as M pairs), unnormalized like the complex transforms (the result is 2M times the
// inverse DFT): E = (X[k] + conj X[M-k]) / 2, O = conj(W^k) (X[k] - conj X[M-k]) / 2, Z = E + i O, and after the
// M-point inverse transform z[2j] + i z[2j+1] = 2 Z[j].
__kernel void fftRealBackward(__global const float2* restrict data, __global float2* restrict potential,
                              __global const float2* restrict twiddle, __global const float2* restrict halfTwiddle,
                              const uint M, const uint radixCode, const uint stages, const uint lineCount,
                              const uint tile, __local float2* bufferA, __local float2* bufferB)
{
  const uint lid = get_local_id(0);
  const uint groupSize = get_local_size(0);
  const uint line0 = get_group_id(0) * tile;
  const uint valid = min(tile, lineCount - line0);
  const uint H = M + 1;
  for (uint e = lid; e < H * tile; e += groupSize)
  {
    const uint k = e % H;
    const uint t = e / H;
    bufferB[k * tile + t] = (t < valid) ? data[(line0 + t) * H + k] : (float2)(0.0f, 0.0f);
  }
  barrier(CLK_LOCAL_MEM_FENCE);
  for (uint e = lid; e < M * tile; e += groupSize)
  {
    const uint k = e % M;
    const uint t = e / M;
    const float2 xk = bufferB[k * tile + t];
    const float2 xm = conjugate(bufferB[(M - k) * tile + t]);
    const float2 even = 0.5f * (xk + xm);
    const float2 odd = cmul(conjugate(halfTwiddle[k]), 0.5f * (xk - xm));
    bufferA[k * tile + t] = (float2)(even.x - odd.y, even.y + odd.x);
  }
  barrier(CLK_LOCAL_MEM_FENCE);
  __local const float2* z = runStages(bufferA, bufferB, twiddle, M, radixCode, stages, tile, 1.0f);
  for (uint e = lid; e < M * tile; e += groupSize)
  {
    const uint j = e % M;
    const uint t = e / M;
    if (t < valid) potential[(line0 + t) * M + j] = 2.0f * z[j * tile + t];
  }
}

// Complex lines of length N along one axis of the half spectrum: element j of line (outer, inner) is at
// outer * outerStride + inner * lineStride + j * axisStride. A work-group transforms `tile` consecutive inner
// lines.
__kernel void fftLines(__global float2* restrict data, __global const float2* restrict twiddle, const uint N,
                       const uint radixCode, const uint stages, const uint axisStride, const uint lineStride,
                       const uint innerCount, const uint outerStride, const uint tile, const uint tilesPerOuter,
                       const float sign, __local float2* bufferA, __local float2* bufferB)
{
  const uint lid = get_local_id(0);
  const uint groupSize = get_local_size(0);
  const uint g = get_group_id(0);
  const uint outer = g / tilesPerOuter;
  const uint inner0 = (g - outer * tilesPerOuter) * tile;
  const uint valid = min(tile, innerCount - inner0);
  const uint base = outer * outerStride + inner0 * lineStride;
  const uint total = N * tile;
  // the loads run along the contiguous direction: along the line when the axis is contiguous, else across the tile
  if (axisStride == 1)
  {
    for (uint e = lid; e < total; e += groupSize)
    {
      const uint j = e % N;
      const uint t = e / N;
      bufferA[j * tile + t] = (t < valid) ? data[base + t * lineStride + j] : (float2)(0.0f, 0.0f);
    }
  }
  else
  {
    for (uint e = lid; e < total; e += groupSize)
    {
      const uint t = e % tile;
      const uint j = e / tile;
      bufferA[e] = (t < valid) ? data[base + t * lineStride + j * axisStride] : (float2)(0.0f, 0.0f);
    }
  }
  barrier(CLK_LOCAL_MEM_FENCE);
  __local const float2* in = runStages(bufferA, bufferB, twiddle, N, radixCode, stages, tile, sign);

  if (axisStride == 1)
  {
    for (uint e = lid; e < total; e += groupSize)
    {
      const uint j = e % N;
      const uint t = e / N;
      if (t < valid) data[base + t * lineStride + j] = in[j * tile + t];
    }
  }
  else
  {
    for (uint e = lid; e < total; e += groupSize)
    {
      const uint t = e % tile;
      const uint j = e / tile;
      if (t < valid) data[base + t * lineStride + j * axisStride] = in[e];
    }
  }
}

// Tree reduction of `count` per-item values held in scratch[q * groupSize + lid] into partials[group * count + q].
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

// G(m) F(m) over the half spectrum kz = 0..Kz/2, with the per-group sums of: the reciprocal energy sum_m G |F|^2,
// its strain derivative (symmetric: xx xy xz yy yz zz), the single-ion sum sum_m bare(m) and the single-ion strain
// tensor (symmetric); the planes 0 < kz < Kz/2 stand for their conjugate partners too and count twice. Every
// work-item handles INFLUENCE_POINTS consecutive points (one index decomposition, then carries) so that the
// per-group reduction is amortized.
__kernel __attribute__((reqd_work_group_size(INFLUENCE_GROUP, 1, 1)))
void applyInfluence(__global float2* restrict data, __global const float* restrict moduliX,
                    __global const float* restrict moduliY, __global const float* restrict moduliZ,
                    __constant const MeshParameters* p, __global float* restrict partials)
{
  __local float scratch[INFLUENCE_PARTIALS * INFLUENCE_GROUP];
  const uint lid = get_local_id(0);
  const uint Kx = p->meshX, Ky = p->meshY, Kz = p->meshZ;
  const uint Hz = Kz / 2 + 1;
  const uint count = Kx * Ky * Hz;
  __constant const float* ic = p->inverseCell;
  float acc[INFLUENCE_PARTIALS];
  for (uint q = 0; q < INFLUENCE_PARTIALS; ++q) acc[q] = 0.0f;
  const uint first = get_global_id(0) * INFLUENCE_POINTS;
  uint mz = first % Hz;
  uint my = (first / Hz) % Ky;
  uint mx = first / (Hz * Ky);
  for (uint i = first; i < min(first + INFLUENCE_POINTS, count); ++i)
  {
    const uint pz = mz, py = my, px = mx;
    if (++mz == Hz)
    {
      mz = 0;
      if (++my == Ky)
      {
        my = 0;
        ++mx;
      }
    }
    if (i == 0)
    {
      data[0] = (float2)(0.0f, 0.0f);
      continue;
    }
    const float sx = (float)((px > Kx / 2) ? (int)px - (int)Kx : (int)px);
    const float sy = (float)((py > Ky / 2) ? (int)py - (int)Ky : (int)py);
    const float sz = (float)pz;  // pz <= Kz / 2
    // k = 2 pi (sx rowX + sy rowY + sz rowZ) with the rows (ax bx cx), (ay by cy), (az bz cz) of the inverse cell
    float3 k;
    k.x = TWO_PI * (sx * ic[0] + sy * ic[1] + sz * ic[2]);
    k.y = TWO_PI * (sx * ic[3] + sy * ic[4] + sz * ic[5]);
    k.z = TWO_PI * (sx * ic[6] + sy * ic[7] + sz * ic[8]);
    const float ksq = dot(k, k);
    const float weight = (pz == 0 || 2 * pz == Kz) ? 1.0f : 2.0f;
    const float g0 = p->prefactor * exp(p->alphaFactor * ksq) / ksq;
    const float g = g0 * moduliX[px] * moduliY[py] * moduliZ[pz];
    const float2 F = data[i];
    data[i] = g * F;
    const float e = weight * g * (F.x * F.x + F.y * F.y);
    const float bare = weight * g0;
    const float shape = 2.0f * (1.0f / ksq + p->inverseFourAlphaSquared);
    const float fac = shape * e;
    const float ionFac = shape * bare;
    acc[0] += e;
    acc[1] += e - fac * k.x * k.x;
    acc[2] -= fac * k.x * k.y;
    acc[3] -= fac * k.x * k.z;
    acc[4] += e - fac * k.y * k.y;
    acc[5] -= fac * k.y * k.z;
    acc[6] += e - fac * k.z * k.z;
    acc[7] += bare;
    acc[8] += bare - ionFac * k.x * k.x;
    acc[9] -= ionFac * k.x * k.y;
    acc[10] -= ionFac * k.x * k.z;
    acc[11] += bare - ionFac * k.y * k.y;
    acc[12] -= ionFac * k.y * k.z;
    acc[13] += bare - ionFac * k.z * k.z;
  }
  for (uint q = 0; q < INFLUENCE_PARTIALS; ++q) scratch[q * INFLUENCE_GROUP + lid] = acc[q];
  reducePartials(scratch, INFLUENCE_PARTIALS, partials);
}

// dE/dr_i = 2 q_i sum_nodes phi(node) d(Mx My Mz)/dr_i, added to the force of the slot
__kernel void interpolateForces(__global const float4* restrict position, __global const float* restrict potential,
                                __constant const MeshParameters* p, __global float4* restrict force)
{
  const uint slot = get_global_id(0);
  if (slot >= p->numberOfSlots) return;
  const float4 r = position[slot];
  const float q = r.w;
  if (q == 0.0f) return;
  const float3 s = fractional(r.xyz, p->inverseCell);
  float wx, wy, wz;
  const int kx0 = anchor(s.x, p->meshX, &wx);
  const int ky0 = anchor(s.y, p->meshY, &wy);
  const int kz0 = anchor(s.z, p->meshZ, &wz);
  float mx[MESH_ORDER], my[MESH_ORDER], mz[MESH_ORDER], dx[MESH_ORDER], dy[MESH_ORDER], dz[MESH_ORDER];
  bsplineWeights(wx, mx, dx);
  bsplineWeights(wy, my, dy);
  bsplineWeights(wz, mz, dz);
  const int Kx = (int)p->meshX, Ky = (int)p->meshY, Kz = (int)p->meshZ;
  float gx = 0.0f, gy = 0.0f, gz = 0.0f;
#pragma unroll
  for (uint a = 0; a < MESH_ORDER; ++a)
  {
    int ix = kx0 - (int)a;
    if (ix < 0) ix += Kx;
#pragma unroll
    for (uint b = 0; b < MESH_ORDER; ++b)
    {
      int iy = ky0 - (int)b;
      if (iy < 0) iy += Ky;
      __global const float* row = potential + ((uint)ix * (uint)Ky + (uint)iy) * (uint)Kz;
      float sumW = 0.0f, sumD = 0.0f;
#pragma unroll
      for (uint c = 0; c < MESH_ORDER; ++c)
      {
        int iz = kz0 - (int)c;
        if (iz < 0) iz += Kz;
        const float phi = row[iz];
        sumW += phi * mz[c];
        sumD += phi * dz[c];
      }
      gx += dx[a] * my[b] * sumW;
      gy += mx[a] * dy[b] * sumW;
      gz += mx[a] * my[b] * sumD;
    }
  }
  __constant const float* ic = p->inverseCell;
  const float fx = gx * (float)Kx, fy = gy * (float)Ky, fz = gz * (float)Kz;
  const float factor = 2.0f * q;
  float4 f = force[slot];
  f.x += factor * (fx * ic[0] + fy * ic[1] + fz * ic[2]);
  f.y += factor * (fx * ic[3] + fy * ic[4] + fz * ic[5]);
  f.z += factor * (fx * ic[6] + fy * ic[7] + fz * ic[8]);
  force[slot] = f;
}
)CLC";
