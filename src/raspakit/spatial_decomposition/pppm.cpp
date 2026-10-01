module;

#include <fftw3.h>

module spatial_decomposition_pppm;

import std;

import int3;
import double3;
import double3x3;
import simulationbox;

namespace
{
std::once_flag fftwThreadsOnce;
std::mutex fftwPlannerMutex;

inline bool sameMatrix(const double3x3& a, const double3x3& b)
{
  return a.ax == b.ax && a.ay == b.ay && a.az == b.az && a.bx == b.bx && a.by == b.by && a.bz == b.bz && a.cx == b.cx &&
         a.cy == b.cy && a.cz == b.cz;
}
}  // namespace

PPPM::~PPPM() { release(); }

void PPPM::release()
{
  if (forwardPlan)
  {
    fftw_destroy_plan(static_cast<fftw_plan>(forwardPlan));
    forwardPlan = nullptr;
  }
  if (backwardPlan)
  {
    fftw_destroy_plan(static_cast<fftw_plan>(backwardPlan));
    backwardPlan = nullptr;
  }
  if (chargeMesh)
  {
    fftw_free(chargeMesh);
    chargeMesh = nullptr;
  }
  threadBuffers.clear();
  if (potential)
  {
    fftw_free(potential);
    potential = nullptr;
  }
  if (spectrum)
  {
    fftw_free(spectrum);
    spectrum = nullptr;
  }
}

std::size_t PPPM::realSize() const
{
  return static_cast<std::size_t>(mesh.x) * static_cast<std::size_t>(mesh.y) * static_cast<std::size_t>(mesh.z);
}

std::size_t PPPM::complexSize() const
{
  return static_cast<std::size_t>(mesh.x) * static_cast<std::size_t>(mesh.y) *
         (static_cast<std::size_t>(mesh.z) / 2 + 1);
}

std::size_t PPPM::nextFFTFriendly(std::size_t n)
{
  std::size_t candidate = std::max<std::size_t>(1, n);
  while (true)
  {
    std::size_t remainder = candidate;
    for (const std::size_t prime : {2uz, 3uz, 5uz})
    {
      while (remainder % prime == 0) remainder /= prime;
    }
    if (remainder == 1) return candidate;
    ++candidate;
  }
}

void PPPM::bsplineWeights(std::size_t p, double w, std::span<double> weights, std::span<double> derivatives)
{
  // M_2(w + j): j = 0 -> w, j = 1 -> 1 - w
  for (std::size_t j = 0; j < p; ++j)
  {
    weights[j] = 0.0;
    derivatives[j] = 0.0;
  }
  weights[0] = w;
  weights[1] = 1.0 - w;
  for (std::size_t k = 3; k <= p; ++k)
  {
    if (k == p)
    {
      // M_p'(u) = M_{p-1}(u) - M_{p-1}(u - 1), evaluated before the last recursion step overwrites M_{p-1}
      derivatives[0] = weights[0];
      for (std::size_t j = 1; j < p; ++j) derivatives[j] = weights[j] - weights[j - 1];
    }
    const double inverse = 1.0 / static_cast<double>(k - 1);
    for (std::size_t j = k - 1; j > 0; --j)
    {
      const double u = w + static_cast<double>(j);
      weights[j] = (u * weights[j] + (static_cast<double>(k) - u) * weights[j - 1]) * inverse;
    }
    weights[0] = w * weights[0] * inverse;
  }
}

void PPPM::initialize(const SimulationBox& box, double alphaValue, double meshSpacing, std::size_t orderValue,
                      std::size_t numberOfThreads, double coulombConversionFactor)
{
  release();

  alpha = alphaValue;
  order = std::clamp<std::size_t>(orderValue, 3, 7);
  conversionFactor = coulombConversionFactor;

  const double lengthA = box.cell[0].length();
  const double lengthB = box.cell[1].length();
  const double lengthC = box.cell[2].length();
  auto size = [&](double length) -> std::int32_t
  {
    const std::size_t minimum =
        std::max<std::size_t>(2 * order, static_cast<std::size_t>(std::ceil(length / meshSpacing)));
    return static_cast<std::int32_t>(nextFFTFriendly(minimum));
  };
  mesh = int3(size(lengthA), size(lengthB), size(lengthC));

  const std::size_t threads = std::max<std::size_t>(1, numberOfThreads);
  const std::size_t n = realSize();
  chargeMesh = fftw_alloc_real(n);
  std::fill_n(chargeMesh, n, 0.0);
  potential = fftw_alloc_real(n);
  std::fill_n(potential, n, 0.0);
  spectrum = fftw_alloc_complex(complexSize());

  threadBuffers.assign(threads, ThreadBuffer{});
  for (ThreadBuffer& buffer : threadBuffers)
  {
    buffer.histogramX.assign(static_cast<std::size_t>(mesh.x), 0u);
    buffer.histogramY.assign(static_cast<std::size_t>(mesh.y), 0u);
    buffer.histogramZ.assign(static_cast<std::size_t>(mesh.z), 0u);
  }

  {
    std::scoped_lock lock(fftwPlannerMutex);
    std::call_once(fftwThreadsOnce, [] { fftw_init_threads(); });
    fftw_plan_with_nthreads(static_cast<int>(threads));
    forwardPlan =
        fftw_plan_dft_r2c_3d(mesh.x, mesh.y, mesh.z, chargeMesh, static_cast<fftw_complex*>(spectrum), FFTW_MEASURE);
    backwardPlan =
        fftw_plan_dft_c2r_3d(mesh.x, mesh.y, mesh.z, static_cast<fftw_complex*>(spectrum), potential, FFTW_MEASURE);
  }
  std::fill_n(chargeMesh, n, 0.0);

  computeBsplineModuli();
  computeInfluenceFunction(box);
}

void PPPM::computeBsplineModuli()
{
  std::array<double, 8> weights{};
  std::array<double, 8> derivatives{};
  bsplineWeights(order, 0.0, std::span<double>(weights.data(), order), std::span<double>(derivatives.data(), order));
  // weights[j] = M_p(j); the denominator uses M_p(k + 1), k = 0..p-2

  auto moduli = [&](std::int32_t K) -> std::vector<double>
  {
    std::vector<double> result(static_cast<std::size_t>(K));
    for (std::int32_t m = 0; m < K; ++m)
    {
      std::complex<double> denominator{0.0, 0.0};
      for (std::size_t k = 0; k + 2 <= order; ++k)
      {
        const double angle =
            2.0 * std::numbers::pi * static_cast<double>(m) * static_cast<double>(k) / static_cast<double>(K);
        denominator += weights[k + 1] * std::complex<double>(std::cos(angle), std::sin(angle));
      }
      const double normSquared = std::norm(denominator);
      result[static_cast<std::size_t>(m)] = normSquared > 1e-14 ? 1.0 / normSquared : 0.0;
    }
    // odd orders have a zero at the Nyquist frequency; interpolate the modulus from the neighbours there
    for (std::int32_t m = 0; m < K; ++m)
    {
      if (result[static_cast<std::size_t>(m)] == 0.0)
      {
        const std::size_t previous = static_cast<std::size_t>((m - 1 + K) % K);
        const std::size_t next = static_cast<std::size_t>((m + 1) % K);
        result[static_cast<std::size_t>(m)] = 0.5 * (result[previous] + result[next]);
      }
    }
    return result;
  };
  bsplineModulusX = moduli(mesh.x);
  bsplineModulusY = moduli(mesh.y);
  bsplineModulusZ = moduli(mesh.z);
}

bool PPPM::influenceFunctionOutdated(const SimulationBox& box, double alphaValue) const
{
  return !sameMatrix(box.inverseCell, inverseCellAtBuild) || box.volume != volume || alphaValue != alpha;
}

void PPPM::beginInfluenceFunction(const SimulationBox& box, double alphaValue, std::size_t numberOfThreads)
{
  alpha = alphaValue;
  inverseCell = box.inverseCell;
  volume = box.volume;
  if (influence.size() != complexSize()) influence.assign(complexSize(), 0.0);
  const std::size_t threads = std::max<std::size_t>(1, numberOfThreads);
  partialIonSum.assign(threads, 0.0);
  partialIonStrain.assign(threads, double3x3{});
}

void PPPM::computeInfluenceSlab(std::size_t thread, std::size_t numberOfThreads)
{
  const std::size_t threads = std::max<std::size_t>(1, numberOfThreads);
  if (thread >= threads) return;
  const std::int32_t planes = mesh.x;
  const std::int32_t first = static_cast<std::int32_t>((static_cast<std::size_t>(planes) * thread) / threads);
  const std::int32_t last = static_cast<std::int32_t>((static_cast<std::size_t>(planes) * (thread + 1)) / threads);

  const double3 rowX(inverseCell.ax, inverseCell.bx, inverseCell.cx);
  const double3 rowY(inverseCell.ay, inverseCell.by, inverseCell.cy);
  const double3 rowZ(inverseCell.az, inverseCell.bz, inverseCell.cz);
  const double prefactor = conversionFactor * 2.0 * std::numbers::pi / volume;
  const double alphaFactor = -0.25 / (alpha * alpha);
  const std::size_t halfZ = static_cast<std::size_t>(mesh.z) / 2 + 1;

  double ionSum = 0.0;
  double3x3 ionTensor{};
  for (std::int32_t mx = first; mx < last; ++mx)
  {
    const std::int32_t sx = (mx > mesh.x / 2) ? mx - mesh.x : mx;
    const double3 kx = 2.0 * std::numbers::pi * static_cast<double>(sx) * rowX;
    for (std::int32_t my = 0; my < mesh.y; ++my)
    {
      const std::int32_t sy = (my > mesh.y / 2) ? my - mesh.y : my;
      const double3 kxy = kx + 2.0 * std::numbers::pi * static_cast<double>(sy) * rowY;
      for (std::int32_t mz = 0; mz < static_cast<std::int32_t>(halfZ); ++mz)
      {
        const std::size_t index =
            (static_cast<std::size_t>(mx) * static_cast<std::size_t>(mesh.y) + static_cast<std::size_t>(my)) * halfZ +
            static_cast<std::size_t>(mz);
        if (mx == 0 && my == 0 && mz == 0)
        {
          influence[index] = 0.0;
          continue;
        }
        const double3 k = kxy + 2.0 * std::numbers::pi * static_cast<double>(mz) * rowZ;
        const double ksq = double3::dot(k, k);
        const double bare = prefactor * std::exp(alphaFactor * ksq) / ksq;
        const double weight = (mz == 0 || (mesh.z % 2 == 0 && mz == mesh.z / 2)) ? 1.0 : 2.0;
        ionSum += weight * bare;
        influence[index] = bare * bsplineModulusX[static_cast<std::size_t>(mx)] *
                           bsplineModulusY[static_cast<std::size_t>(my)] *
                           bsplineModulusZ[static_cast<std::size_t>(mz)];

        // strain response of one wave vector: bare [I - 2 (1/k^2 + 1/(4 alpha^2)) k k^T]
        const double fac = 2.0 * (1.0 / ksq + 0.25 / (alpha * alpha)) * weight * bare;
        ionTensor.ax += weight * bare - fac * k.x * k.x;
        ionTensor.bx += -fac * k.x * k.y;
        ionTensor.cx += -fac * k.x * k.z;
        ionTensor.ay += -fac * k.y * k.x;
        ionTensor.by += weight * bare - fac * k.y * k.y;
        ionTensor.cy += -fac * k.y * k.z;
        ionTensor.az += -fac * k.z * k.x;
        ionTensor.bz += -fac * k.z * k.y;
        ionTensor.cz += weight * bare - fac * k.z * k.z;
      }
    }
  }
  partialIonSum[thread] = ionSum;
  partialIonStrain[thread] = ionTensor;
}

void PPPM::finishInfluenceFunction()
{
  singleIonSum = 0.0;
  ionStrain = double3x3{};
  for (std::size_t t = 0; t < partialIonSum.size(); ++t)
  {
    singleIonSum += partialIonSum[t];
    ionStrain += partialIonStrain[t];
  }
  inverseCellAtBuild = inverseCell;
}

void PPPM::computeInfluenceFunction(const SimulationBox& box)
{
  beginInfluenceFunction(box, alpha, 1);
  computeInfluenceSlab(0, 1);
  finishInfluenceFunction();
}

void PPPM::updateBox(const SimulationBox& box, double alphaValue)
{
  if (influenceFunctionOutdated(box, alphaValue))
  {
    alpha = alphaValue;
    computeInfluenceFunction(box);
  }
}

namespace
{
// The smallest periodic arc of mesh planes (start, length) that holds every plane with a non-zero count, extended
// by `halo` planes below the first one (the B-spline support of an anchor reaches the planes anchor - order + 1
// .. anchor). The arc is the complement of the largest circular run of empty planes; it degenerates to the whole
// axis when the extended arc would cover it. Returns length 0 when no plane is occupied.
std::pair<std::int32_t, std::int32_t> coveringArc(std::span<const std::uint32_t> histogram, std::int32_t halo)
{
  const std::int32_t K = static_cast<std::int32_t>(histogram.size());
  std::int32_t first = -1;
  for (std::int32_t k = 0; k < K; ++k)
  {
    if (histogram[static_cast<std::size_t>(k)] != 0)
    {
      first = k;
      break;
    }
  }
  if (first < 0) return {0, 0};

  // longest run of empty planes, scanned circularly from an occupied plane so that no run is split
  std::int32_t longestGap = 0;
  std::int32_t longestGapEnd = first;  // the occupied plane that follows the longest gap
  std::int32_t gap = 0;
  for (std::int32_t step = 1; step <= K; ++step)
  {
    const std::int32_t k = (first + step) % K;
    if (histogram[static_cast<std::size_t>(k)] == 0)
    {
      ++gap;
    }
    else
    {
      if (gap > longestGap)
      {
        longestGap = gap;
        longestGapEnd = k;
      }
      gap = 0;
    }
  }
  std::int32_t start = longestGapEnd;
  std::int32_t length = K - longestGap;

  length += halo;
  if (length >= K) return {0, K};
  start = ((start - halo) % K + K) % K;
  return {start, length};
}
}  // namespace

void PPPM::spread(std::size_t thread, std::span<const std::uint32_t> atoms, const double* x, const double* y,
                  const double* z, const double* charge, const double* scalingCoulomb)
{
  ThreadBuffer& buffer = threadBuffers[thread];
  const std::size_t p = order;
  const std::int32_t Kx = mesh.x;
  const std::int32_t Ky = mesh.y;
  const std::int32_t Kz = mesh.z;
  const bool direct = (threadBuffers.size() == 1);

  // pass 1: anchors (the mesh point at or below the atom along every axis) and their counts per plane
  buffer.anchors.clear();
  std::fill(buffer.histogramX.begin(), buffer.histogramX.end(), 0u);
  std::fill(buffer.histogramY.begin(), buffer.histogramY.end(), 0u);
  std::fill(buffer.histogramZ.begin(), buffer.histogramZ.end(), 0u);
  for (const std::uint32_t i : atoms)
  {
    const double q = scalingCoulomb[i] * charge[i];
    if (q == 0.0) continue;
    const double3 r(x[i], y[i], z[i]);
    double3 s = inverseCell * r;
    s.x -= std::floor(s.x);
    s.y -= std::floor(s.y);
    s.z -= std::floor(s.z);
    const double ux = std::min(s.x * static_cast<double>(Kx), std::nextafter(static_cast<double>(Kx), 0.0));
    const double uy = std::min(s.y * static_cast<double>(Ky), std::nextafter(static_cast<double>(Ky), 0.0));
    const double uz = std::min(s.z * static_cast<double>(Kz), std::nextafter(static_cast<double>(Kz), 0.0));
    buffer.anchors.insert(buffer.anchors.end(), {q, ux, uy, uz});
    ++buffer.histogramX[static_cast<std::size_t>(ux)];
    ++buffer.histogramY[static_cast<std::size_t>(uy)];
    ++buffer.histogramZ[static_cast<std::size_t>(uz)];
  }

  // the sub-box of the mesh touched by these atoms: the anchors' arcs plus the order - 1 planes of the spline
  // support below them; with one thread the FFT input itself is the target
  double* grid = nullptr;
  if (direct)
  {
    buffer.start = int3(0, 0, 0);
    buffer.length = mesh;
    grid = chargeMesh;
    std::fill_n(chargeMesh, realSize(), 0.0);
  }
  else
  {
    const std::int32_t halo = static_cast<std::int32_t>(p) - 1;
    const auto [sx, lx] = coveringArc(buffer.histogramX, halo);
    const auto [sy, ly] = coveringArc(buffer.histogramY, halo);
    const auto [sz, lz] = coveringArc(buffer.histogramZ, halo);
    buffer.start = int3(sx, sy, sz);
    buffer.length = int3(lx, ly, lz);
    if (buffer.anchors.empty())
    {
      buffer.length = int3(0, 0, 0);
      return;
    }
    const std::size_t size = static_cast<std::size_t>(lx) * static_cast<std::size_t>(ly) * static_cast<std::size_t>(lz);
    buffer.values.assign(size, 0.0);
    grid = buffer.values.data();
  }

  const std::int32_t Ly = buffer.length.y;
  const std::int32_t Lz = buffer.length.z;
  const std::int32_t startX = buffer.start.x;
  const std::int32_t startY = buffer.start.y;
  const std::int32_t startZ = buffer.start.z;

  // pass 2: spread into the buffer with local (sub-box) indices
  std::array<double, 8> wx{}, wy{}, wz{}, dummy{};
  std::array<std::size_t, 8> ix{}, iy{}, iz{};
  const double* anchor = buffer.anchors.data();
  const std::size_t count = buffer.anchors.size() / 4;
  for (std::size_t n = 0; n < count; ++n, anchor += 4)
  {
    const double q = anchor[0];
    const double ux = anchor[1];
    const double uy = anchor[2];
    const double uz = anchor[3];
    const std::int32_t kx0 = static_cast<std::int32_t>(ux);
    const std::int32_t ky0 = static_cast<std::int32_t>(uy);
    const std::int32_t kz0 = static_cast<std::int32_t>(uz);
    bsplineWeights(p, ux - static_cast<double>(kx0), std::span<double>(wx.data(), p),
                   std::span<double>(dummy.data(), p));
    bsplineWeights(p, uy - static_cast<double>(ky0), std::span<double>(wy.data(), p),
                   std::span<double>(dummy.data(), p));
    bsplineWeights(p, uz - static_cast<double>(kz0), std::span<double>(wz.data(), p),
                   std::span<double>(dummy.data(), p));
    for (std::size_t j = 0; j < p; ++j)
    {
      // mesh index (k0 - j) mod K, relative to the start of the arc; inside the arc by construction
      const std::int32_t dj = static_cast<std::int32_t>(j);
      ix[j] = static_cast<std::size_t>(((kx0 - dj - startX) % Kx + Kx) % Kx);
      iy[j] = static_cast<std::size_t>(((ky0 - dj - startY) % Ky + Ky) % Ky);
      iz[j] = static_cast<std::size_t>(((kz0 - dj - startZ) % Kz + Kz) % Kz);
    }
    for (std::size_t a = 0; a < p; ++a)
    {
      const double qx = q * wx[a];
      const std::size_t offsetX = ix[a] * static_cast<std::size_t>(Ly);
      for (std::size_t b = 0; b < p; ++b)
      {
        const double qxy = qx * wy[b];
        double* row = grid + (offsetX + iy[b]) * static_cast<std::size_t>(Lz);
        for (std::size_t c = 0; c < p; ++c)
        {
          row[iz[c]] += qxy * wz[c];
        }
      }
    }
  }
}

void PPPM::reduceMeshes(std::size_t thread, std::size_t numberOfThreads)
{
  if (threadBuffers.size() <= 1) return;
  const std::int32_t Kx = mesh.x;
  const std::int32_t Ky = mesh.y;
  const std::int32_t Kz = mesh.z;
  const std::size_t planeSize = static_cast<std::size_t>(Ky) * static_cast<std::size_t>(Kz);
  const std::int32_t firstPlane = static_cast<std::int32_t>((static_cast<std::size_t>(Kx) * thread) / numberOfThreads);
  const std::int32_t lastPlane =
      static_cast<std::int32_t>((static_cast<std::size_t>(Kx) * (thread + 1)) / numberOfThreads);

  for (std::int32_t kx = firstPlane; kx < lastPlane; ++kx)
  {
    double* plane = chargeMesh + static_cast<std::size_t>(kx) * planeSize;
    std::fill_n(plane, planeSize, 0.0);

    for (const ThreadBuffer& buffer : threadBuffers)
    {
      if (buffer.length.x == 0) continue;
      const std::int32_t lx = ((kx - buffer.start.x) % Kx + Kx) % Kx;
      if (lx >= buffer.length.x) continue;
      const std::int32_t Ly = buffer.length.y;
      const std::int32_t Lz = buffer.length.z;
      // the z-arc of the buffer as at most two contiguous runs of the mesh axis
      const std::int32_t firstRun = std::min(Lz, Kz - buffer.start.z);
      const std::int32_t secondRun = Lz - firstRun;
      const double* slab = buffer.values.data() +
                           static_cast<std::size_t>(lx) * static_cast<std::size_t>(Ly) * static_cast<std::size_t>(Lz);
      for (std::int32_t ly = 0; ly < Ly; ++ly)
      {
        const std::int32_t ky = (buffer.start.y + ly) % Ky;
        double* destination = plane + static_cast<std::size_t>(ky) * static_cast<std::size_t>(Kz);
        const double* source = slab + static_cast<std::size_t>(ly) * static_cast<std::size_t>(Lz);
        double* firstDestination = destination + buffer.start.z;
        for (std::int32_t lz = 0; lz < firstRun; ++lz) firstDestination[lz] += source[lz];
        const double* wrapped = source + firstRun;
        for (std::int32_t lz = 0; lz < secondRun; ++lz) destination[lz] += wrapped[lz];
      }
    }
  }
}

double PPPM::solve(bool withVirial)
{
  fftw_execute(static_cast<fftw_plan>(forwardPlan));

  fftw_complex* F = static_cast<fftw_complex*>(spectrum);
  const std::size_t halfZ = static_cast<std::size_t>(mesh.z) / 2 + 1;
  const std::size_t n = complexSize();
  double energy = 0.0;
  reciprocalStrain = double3x3{};

  if (!withVirial)
  {
    for (std::size_t index = 0; index < n; ++index)
    {
      const std::size_t mz = index % halfZ;
      const double g = influence[index];
      const double re = F[index][0];
      const double im = F[index][1];
      const double weight = (mz == 0 || (mesh.z % 2 == 0 && mz == static_cast<std::size_t>(mesh.z) / 2)) ? 1.0 : 2.0;
      energy += weight * g * (re * re + im * im);
      F[index][0] = g * re;
      F[index][1] = g * im;
    }
  }
  else
  {
    const double3 rowX(inverseCell.ax, inverseCell.bx, inverseCell.cx);
    const double3 rowY(inverseCell.ay, inverseCell.by, inverseCell.cy);
    const double3 rowZ(inverseCell.az, inverseCell.bz, inverseCell.cz);
    const double inverseFourAlphaSquared = 0.25 / (alpha * alpha);
    double3x3 tensor{};
    std::size_t index = 0;
    for (std::int32_t mx = 0; mx < mesh.x; ++mx)
    {
      const std::int32_t sx = (mx > mesh.x / 2) ? mx - mesh.x : mx;
      const double3 kx = 2.0 * std::numbers::pi * static_cast<double>(sx) * rowX;
      for (std::int32_t my = 0; my < mesh.y; ++my)
      {
        const std::int32_t sy = (my > mesh.y / 2) ? my - mesh.y : my;
        const double3 kxy = kx + 2.0 * std::numbers::pi * static_cast<double>(sy) * rowY;
        for (std::int32_t mz = 0; mz < static_cast<std::int32_t>(halfZ); ++mz, ++index)
        {
          const double g = influence[index];
          const double re = F[index][0];
          const double im = F[index][1];
          F[index][0] = g * re;
          F[index][1] = g * im;
          if (g == 0.0) continue;
          const double weight = (mz == 0 || (mesh.z % 2 == 0 && mz == mesh.z / 2)) ? 1.0 : 2.0;
          const double e = weight * g * (re * re + im * im);
          energy += e;
          const double3 k = kxy + 2.0 * std::numbers::pi * static_cast<double>(mz) * rowZ;
          const double ksq = double3::dot(k, k);
          const double fac = 2.0 * (1.0 / ksq + inverseFourAlphaSquared) * e;
          tensor.ax += e - fac * k.x * k.x;
          tensor.bx += -fac * k.x * k.y;
          tensor.cx += -fac * k.x * k.z;
          tensor.ay += -fac * k.y * k.x;
          tensor.by += e - fac * k.y * k.y;
          tensor.cy += -fac * k.y * k.z;
          tensor.az += -fac * k.z * k.x;
          tensor.bz += -fac * k.z * k.y;
          tensor.cz += e - fac * k.z * k.z;
        }
      }
    }
    reciprocalStrain = tensor;
  }

  fftw_execute(static_cast<fftw_plan>(backwardPlan));
  return energy;
}

void PPPM::interpolate(std::span<const std::uint32_t> atoms, const double* x, const double* y, const double* z,
                       const double* charge, const double* scalingCoulomb, double* fx, double* fy, double* fz) const
{
  const std::size_t p = order;
  const std::size_t Kx = static_cast<std::size_t>(mesh.x);
  const std::size_t Ky = static_cast<std::size_t>(mesh.y);
  const std::size_t Kz = static_cast<std::size_t>(mesh.z);
  const double3 rowX(inverseCell.ax, inverseCell.bx, inverseCell.cx);
  const double3 rowY(inverseCell.ay, inverseCell.by, inverseCell.cy);
  const double3 rowZ(inverseCell.az, inverseCell.bz, inverseCell.cz);
  const double3 scaledRowX = static_cast<double>(Kx) * rowX;
  const double3 scaledRowY = static_cast<double>(Ky) * rowY;
  const double3 scaledRowZ = static_cast<double>(Kz) * rowZ;
  std::array<double, 8> wx{}, wy{}, wz{}, dx{}, dy{}, dz{};
  std::array<std::size_t, 8> ix{}, iy{}, iz{};

  for (const std::uint32_t i : atoms)
  {
    const double q = scalingCoulomb[i] * charge[i];
    if (q == 0.0) continue;
    const double3 r(x[i], y[i], z[i]);
    double3 s = inverseCell * r;
    s.x -= std::floor(s.x);
    s.y -= std::floor(s.y);
    s.z -= std::floor(s.z);
    const double ux = s.x * static_cast<double>(Kx);
    const double uy = s.y * static_cast<double>(Ky);
    const double uz = s.z * static_cast<double>(Kz);
    const std::size_t kx0 = std::min(static_cast<std::size_t>(ux), Kx - 1);
    const std::size_t ky0 = std::min(static_cast<std::size_t>(uy), Ky - 1);
    const std::size_t kz0 = std::min(static_cast<std::size_t>(uz), Kz - 1);
    bsplineWeights(p, ux - static_cast<double>(kx0), std::span<double>(wx.data(), p), std::span<double>(dx.data(), p));
    bsplineWeights(p, uy - static_cast<double>(ky0), std::span<double>(wy.data(), p), std::span<double>(dy.data(), p));
    bsplineWeights(p, uz - static_cast<double>(kz0), std::span<double>(wz.data(), p), std::span<double>(dz.data(), p));
    for (std::size_t j = 0; j < p; ++j)
    {
      ix[j] = (kx0 + Kx - j) % Kx;
      iy[j] = (ky0 + Ky - j) % Ky;
      iz[j] = (kz0 + Kz - j) % Kz;
    }

    // dE/du_alpha = 2 q sum_nodes phi * d(Mx My Mz)/du_alpha
    double gx = 0.0, gy = 0.0, gz = 0.0;
    for (std::size_t a = 0; a < p; ++a)
    {
      const std::size_t offsetX = ix[a] * Ky;
      for (std::size_t b = 0; b < p; ++b)
      {
        const double* row = potential + (offsetX + iy[b]) * Kz;
        double sumW = 0.0, sumD = 0.0;
        for (std::size_t c = 0; c < p; ++c)
        {
          const double phi = row[iz[c]];
          sumW += phi * wz[c];
          sumD += phi * dz[c];
        }
        gx += dx[a] * wy[b] * sumW;
        gy += wx[a] * dy[b] * sumW;
        gz += wx[a] * wy[b] * sumD;
      }
    }
    const double factor = 2.0 * q;
    const double3 gradient = factor * (gx * scaledRowX + gy * scaledRowY + gz * scaledRowZ);
    fx[i] += gradient.x;
    fy[i] += gradient.y;
    fz[i] += gradient.z;
  }
}

std::string PPPM::status() const
{
  std::string result = std::format("    particle-mesh Ewald: mesh {} x {} x {}, B-spline order {}, alpha {:.6f} A^-1\n",
                                   mesh.x, mesh.y, mesh.z, order, alpha);
  if (threadBuffers.size() > 1)
  {
    std::size_t bufferPoints = 0;
    for (const ThreadBuffer& buffer : threadBuffers)
    {
      bufferPoints += static_cast<std::size_t>(buffer.length.x) * static_cast<std::size_t>(buffer.length.y) *
                      static_cast<std::size_t>(buffer.length.z);
    }
    result += std::format("    charge assignment: {} per-thread sub-box buffers, together {:.2f} x the mesh\n",
                          threadBuffers.size(),
                          realSize() > 0 ? static_cast<double>(bufferPoints) / static_cast<double>(realSize()) : 0.0);
  }
  return result;
}
