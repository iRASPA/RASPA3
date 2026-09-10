module;

#include <fftw3.h>

module energy_shared_nldft;

import std;

import uint3;
import double3;
import double3x3;
import unit_cell;
import energy_shared_linear_probe;
import energy_shared_bet_surface_area;

// The uniform limit of the functional and the FFT machinery of the inhomogeneous solve. Everything here
// works in Kelvin energy units (k_B = 1): densities in molecules / Å³, pressures in K / Å³, converted to
// pascal only at the interface.

namespace
{
constexpr double pressureKA3ToPascal = 1.380649e7;  // (K / Å³) × k_B / Å³ in Pa

unsigned nldftWorkerCount()
{
  unsigned n = std::thread::hardware_concurrency();
  return n == 0 ? 1u : n;
}

template <typename Body>
void nldftParallelFor(std::size_t n, unsigned workers, Body &&body)
{
  if (workers <= 1 || n < 2048)
  {
    body(0, n);
    return;
  }
  std::vector<std::thread> pool;
  pool.reserve(workers);
  const std::size_t chunk = (n + workers - 1) / workers;
  for (unsigned t = 0; t < workers; ++t)
  {
    const std::size_t begin = static_cast<std::size_t>(t) * chunk;
    const std::size_t end = std::min(n, begin + chunk);
    if (begin < end) pool.emplace_back([&body, begin, end] { body(begin, end); });
  }
  for (std::thread &thread : pool) thread.join();
}

int nldftWrap(int i, int n)
{
  int r = i % n;
  return r < 0 ? r + n : r;
}

double nldftWrap01(double x)
{
  x -= std::floor(x);
  return x < 0.0 ? x + 1.0 : x;
}

void nldftCicAdd(double *grid, std::size_t nx, std::size_t ny, std::size_t nz, double fx, double fy, double fz,
                 double weight)
{
  if (!(weight > 0.0) || !(weight < 1.0e300)) return;
  fx = nldftWrap01(fx);
  fy = nldftWrap01(fy);
  fz = nldftWrap01(fz);
  const double x = fx * static_cast<double>(nx);
  const double y = fy * static_cast<double>(ny);
  const double z = fz * static_cast<double>(nz);
  const int i0 = static_cast<int>(std::floor(x));
  const int j0 = static_cast<int>(std::floor(y));
  const int k0 = static_cast<int>(std::floor(z));
  const int i1 = i0 + 1;
  const int j1 = j0 + 1;
  const int k1 = k0 + 1;
  const double tx = x - static_cast<double>(i0);
  const double ty = y - static_cast<double>(j0);
  const double tz = z - static_cast<double>(k0);
  const int ii0 = nldftWrap(i0, static_cast<int>(nx));
  const int ii1 = nldftWrap(i1, static_cast<int>(nx));
  const int jj0 = nldftWrap(j0, static_cast<int>(ny));
  const int jj1 = nldftWrap(j1, static_cast<int>(ny));
  const int kk0 = nldftWrap(k0, static_cast<int>(nz));
  const int kk1 = nldftWrap(k1, static_cast<int>(nz));
  const double w000 = (1.0 - tx) * (1.0 - ty) * (1.0 - tz);
  const double w100 = tx * (1.0 - ty) * (1.0 - tz);
  const double w010 = (1.0 - tx) * ty * (1.0 - tz);
  const double w110 = tx * ty * (1.0 - tz);
  const double w001 = (1.0 - tx) * (1.0 - ty) * tz;
  const double w101 = tx * (1.0 - ty) * tz;
  const double w011 = (1.0 - tx) * ty * tz;
  const double w111 = tx * ty * tz;
  auto at = [&](int i, int j, int k) -> double &
  { return grid[(static_cast<std::size_t>(k) * ny + static_cast<std::size_t>(j)) * nx + static_cast<std::size_t>(i)]; };
  at(ii0, jj0, kk0) += weight * w000;
  at(ii1, jj0, kk0) += weight * w100;
  at(ii0, jj1, kk0) += weight * w010;
  at(ii1, jj1, kk0) += weight * w110;
  at(ii0, jj0, kk1) += weight * w001;
  at(ii1, jj0, kk1) += weight * w101;
  at(ii0, jj1, kk1) += weight * w011;
  at(ii1, jj1, kk1) += weight * w111;
}

double nldftCicSample(const double *grid, std::size_t nx, std::size_t ny, std::size_t nz, double fx, double fy,
                      double fz)
{
  fx = nldftWrap01(fx);
  fy = nldftWrap01(fy);
  fz = nldftWrap01(fz);
  const double x = fx * static_cast<double>(nx);
  const double y = fy * static_cast<double>(ny);
  const double z = fz * static_cast<double>(nz);
  const int i0 = static_cast<int>(std::floor(x));
  const int j0 = static_cast<int>(std::floor(y));
  const int k0 = static_cast<int>(std::floor(z));
  const double tx = x - static_cast<double>(i0);
  const double ty = y - static_cast<double>(j0);
  const double tz = z - static_cast<double>(k0);
  const int ii0 = nldftWrap(i0, static_cast<int>(nx));
  const int ii1 = nldftWrap(i0 + 1, static_cast<int>(nx));
  const int jj0 = nldftWrap(j0, static_cast<int>(ny));
  const int jj1 = nldftWrap(j0 + 1, static_cast<int>(ny));
  const int kk0 = nldftWrap(k0, static_cast<int>(nz));
  const int kk1 = nldftWrap(k0 + 1, static_cast<int>(nz));
  auto at = [&](int i, int j, int k)
  { return grid[(static_cast<std::size_t>(k) * ny + static_cast<std::size_t>(j)) * nx + static_cast<std::size_t>(i)]; };
  const double c00 = at(ii0, jj0, kk0) * (1.0 - tx) + at(ii1, jj0, kk0) * tx;
  const double c10 = at(ii0, jj1, kk0) * (1.0 - tx) + at(ii1, jj1, kk0) * tx;
  const double c01 = at(ii0, jj0, kk1) * (1.0 - tx) + at(ii1, jj0, kk1) * tx;
  const double c11 = at(ii0, jj1, kk1) * (1.0 - tx) + at(ii1, jj1, kk1) * tx;
  const double c0 = c00 * (1.0 - ty) + c10 * ty;
  const double c1 = c01 * (1.0 - ty) + c11 * ty;
  return c0 * (1.0 - tz) + c1 * tz;
}

// The uniform fluid: Carnahan-Starling hard spheres (the White Bear functional's own uniform limit) plus
// the mean-field attraction (1/2) a rho².
struct BulkModel
{
  double T{0.0};
  double packingFactor{0.0};  // π d³ / 6 per hard-sphere site
  double sites{1.0};          // 1 for spherical FMT; 2 for a TraPPE N2 dumbbell
  double a{0.0};              // int u_att d³r, K Å³, negative

  double packing(double rho) const { return sites * packingFactor * rho; }

  double muOf(double rho) const
  {
    double e = packing(rho);
    double o = 1.0 - e;
    return T * std::log(rho) + sites * T * e * (8.0 - 9.0 * e + 3.0 * e * e) / (o * o * o) + a * rho;
  }
  double pressureOf(double rho) const
  {
    double e = packing(rho);
    double o = 1.0 - e;
    double zCS = (1.0 + e + e * e - e * e * e) / (o * o * o);
    return T * rho * (1.0 - sites + sites * zCS) + 0.5 * a * rho * rho;
  }
  double muDerivative(double rho) const
  {
    double e = packing(rho);
    double o = 1.0 - e;
    double deta = sites * packingFactor;
    return T / rho + sites * T * deta * (8.0 - 2.0 * e) / (o * o * o * o) + a;
  }
};

// The WCA split of the Lennard-Jones pair: constant -epsilon inside the minimum, the bare tail outside,
// truncated at the cutoff. Used for spherical Ravikovitch N2.
double attractiveTail(double r, double sigma, double epsilon, double rMin, double cutoff)
{
  if (r >= cutoff) return 0.0;
  if (r < rMin) return -epsilon;
  double s3 = (sigma / r) * (sigma / r) * (sigma / r);
  double s6 = s3 * s3;
  return 4.0 * epsilon * (s6 * s6 - s6);
}

// Full 12-6 from contact (u=0 at r=σ) to the cutoff. The r < σ core is left to FMT; integrating it
// makes a diverge. This is TraPPE N–N as written, not the WCA flattening to −ε inside rmin.
double fullLennardJones(double r, double sigma, double epsilon, double cutoff)
{
  if (!(r >= sigma) || r >= cutoff) return 0.0;
  double s3 = (sigma / r) * (sigma / r) * (sigma / r);
  double s6 = s3 * s3;
  return 4.0 * epsilon * (s6 * s6 - s6);
}

double ljTailPrimitive(double r, double sigma, double epsilon)
{
  double s3 = sigma * sigma * sigma;
  double s6 = s3 * s3;
  double s12 = s6 * s6;
  double r3 = r * r * r;
  double r9 = r3 * r3 * r3;
  return 4.0 * epsilon * (-s12 / (9.0 * r9) + s6 / (3.0 * r3));
}

// a = int u_att d³r, analytic. WCA includes the −ε ball inside rMin.
double meanFieldIntegral(double sigma, double epsilon, double rMin, double cutoff)
{
  double core = -epsilon * rMin * rMin * rMin / 3.0;
  return 4.0 * std::numbers::pi * (core + ljTailPrimitive(cutoff, sigma, epsilon) - ljTailPrimitive(rMin, sigma, epsilon));
}

// a = int u_LJ d³r from r = σ to the cutoff (full 12-6, no WCA core).
double meanFieldIntegralFullLJ(double sigma, double epsilon, double cutoff)
{
  return 4.0 * std::numbers::pi * (ljTailPrimitive(cutoff, sigma, epsilon) - ljTailPrimitive(sigma, sigma, epsilon));
}

bool siteAttraction(const NLDFTOptions &options)
{
  return !options.attractiveSiteOffsets.empty() && options.attractiveSigma > 0.0 && options.attractiveEpsilon > 0.0;
}

double hardSphereVolume(double d) { return std::numbers::pi * d * d * d / 6.0; }

// Lens overlap of two spheres of diameter d whose centres are `separation` apart.
double spherePairOverlap(double d, double separation)
{
  const double R = 0.5 * d;
  if (!(separation > 0.0) || separation >= 2.0 * R) return 0.0;
  return std::numbers::pi * (2.0 * R - separation) * (2.0 * R - separation) * (separation + 4.0 * R) /
         (12.0 * separation);
}

// Union volume of the dispersing-site spheres along the molecular axis. Pairwise inclusion-exclusion
// is exact for a dumbbell; for three or more sites it ignores triple overlap (N2 has two).
double fusedSiteVolume(double d, const std::vector<double> &offsets)
{
  const double V = hardSphereVolume(d);
  if (offsets.size() <= 1) return V;
  double unionV = static_cast<double>(offsets.size()) * V;
  for (std::size_t i = 0; i < offsets.size(); ++i)
    for (std::size_t j = i + 1; j < offsets.size(); ++j)
      unionV -= spherePairOverlap(d, std::fabs(offsets[j] - offsets[i]));
  return std::max(unionV, V);
}

BulkModel bulkModelOf(const NLDFTOptions &options)
{
  BulkModel model;
  model.T = options.temperature;
  const bool dumbbell = !options.dumbbellSiteOffsets.empty() && options.dumbbellSiteDiameter > 0.0;
  const double d = dumbbell ? options.dumbbellSiteDiameter : options.hardSphereDiameter;
  model.packingFactor = hardSphereVolume(d);
  if (dumbbell)
  {
    // Independent-site η = n_sites ρ π d³/6 overcounts the fused dumbbell (two TraPPE N spheres
    // overlap by ~9 Å³). That fluid is supercritical at 77 K and μ_hs oscillates in 8-rings.
    // Scale the CS segment count so η = ρ V_union, matching the uniform White Bear n3 of
    // overlap-weighted site deposition.
    const double unionV = fusedSiteVolume(d, options.dumbbellSiteOffsets);
    model.sites = unionV / model.packingFactor;
  }
  else
  {
    model.sites = 1.0;
  }
  if (model.sites < 1.0) model.sites = 1.0;
  if (siteAttraction(options))
  {
    const double nAtt = static_cast<double>(options.attractiveSiteOffsets.size());
    const double aSite =
        meanFieldIntegralFullLJ(options.attractiveSigma, options.attractiveEpsilon, options.cutoff);
    model.a = nAtt * nAtt * aSite;
  }
  else
  {
    double rMin = std::pow(2.0, 1.0 / 6.0) * options.sigma;
    model.a = meanFieldIntegral(options.sigma, options.epsilon, rMin, options.cutoff);
  }
  return model;
}

// The two spinodals of the subcritical loop, by scanning mu'(rho) for sign changes and bisecting each.
bool findSpinodals(const BulkModel &model, double &lowSpinodal, double &highSpinodal)
{
  const double rhoMax = 0.60 / (model.sites * model.packingFactor);
  const std::size_t points = 20000;
  double previousRho = rhoMax * 1.0e-8;
  double previousSign = model.muDerivative(previousRho);
  std::vector<double> crossings;
  for (std::size_t i = 1; i <= points; ++i)
  {
    double rho = rhoMax * static_cast<double>(i) / static_cast<double>(points);
    double sign = model.muDerivative(rho);
    if ((previousSign > 0.0) != (sign > 0.0))
    {
      double lo = previousRho;
      double hi = rho;
      for (int iteration = 0; iteration < 80; ++iteration)
      {
        double mid = 0.5 * (lo + hi);
        if ((model.muDerivative(mid) > 0.0) == (previousSign > 0.0))
          lo = mid;
        else
          hi = mid;
      }
      crossings.push_back(0.5 * (lo + hi));
    }
    previousRho = rho;
    previousSign = sign;
  }
  if (crossings.size() < 2) return false;
  lowSpinodal = crossings.front();
  highSpinodal = crossings.back();
  return true;
}

// The density at which the given branch reaches the chemical potential mu, by bisection on a monotone
// stretch.
double densityAtMu(const BulkModel &model, double mu, double rhoLow, double rhoHigh)
{
  double lo = rhoLow;
  double hi = rhoHigh;
  for (int iteration = 0; iteration < 100; ++iteration)
  {
    double mid = 0.5 * (lo + hi);
    if (model.muOf(mid) < mu)
      lo = mid;
    else
      hi = mid;
  }
  return 0.5 * (lo + hi);
}

double gasDensityAtPressure(const BulkModel &model, double pressure, double rhoHigh)
{
  double lo = 1.0e-18;
  double hi = rhoHigh;
  for (int iteration = 0; iteration < 100; ++iteration)
  {
    double mid = 0.5 * (lo + hi);
    if (model.pressureOf(mid) < pressure)
      lo = mid;
    else
      hi = mid;
  }
  return 0.5 * (lo + hi);
}

// White Bear phi3(n3) and its derivative, with the small-n3 series where the exact form loses digits.
void phi3Of(double n3, double &value, double &derivative)
{
  constexpr double c = 1.0 / (36.0 * std::numbers::pi);
  double om = 1.0 - n3;
  if (n3 < 1.0e-3)
  {
    double g = 1.5 - n3 / 3.0 - n3 * n3 / 12.0;
    double gp = -1.0 / 3.0 - n3 / 6.0;
    value = c * g / (om * om);
    derivative = c * (gp / (om * om) + 2.0 * g / (om * om * om));
  }
  else
  {
    double lnOm = std::log(om);
    double numerator = n3 + om * om * lnOm;
    double denominator = n3 * n3 * om * om;
    double numeratorPrime = 1.0 - 2.0 * om * lnOm - om;
    double denominatorPrime = 2.0 * n3 * om * om - 2.0 * n3 * n3 * om;
    value = c * numerator / denominator;
    derivative = c * (numeratorPrime * denominator - numerator * denominatorPrime) / (denominator * denominator);
  }
}

std::once_flag fftwThreadsFlag;

}  // namespace


NLDFTBulk nldftBulkCoexistence(const NLDFTOptions &options)
{
  NLDFTBulk bulk;
  BulkModel model = bulkModelOf(options);
  bulk.meanFieldIntegral = model.a;

  double lowSpinodal = 0.0;
  double highSpinodal = 0.0;
  if (!findSpinodals(model, lowSpinodal, highSpinodal))
  {
    // No van der Waals loop at this T (the HS repulsion outruns the WCA tail). Keep the same μ(ρ)
    // and P(ρ) the pore FMT uses, and measure x against experimental nitrogen P0.
    bulk.experimentalSaturation = true;
    bulk.saturationPressure = nitrogenSaturationPressure;
    bulk.liquidDensity = 0.45 / (model.sites * model.packingFactor);
    bulk.gasDensity = gasDensityAtPressure(model, bulk.saturationPressure / pressureKA3ToPascal, bulk.liquidDensity);
    return bulk;
  }

  const double rhoMax = 0.60 / (model.sites * model.packingFactor);
  double muLow = model.muOf(highSpinodal);   // the loop's minimum
  double muHigh = model.muOf(lowSpinodal);   // the loop's maximum
  for (int iteration = 0; iteration < 100; ++iteration)
  {
    double mu = 0.5 * (muLow + muHigh);
    double rhoGas = densityAtMu(model, mu, 1.0e-18, lowSpinodal);
    double rhoLiquid = densityAtMu(model, mu, highSpinodal, rhoMax);
    if (model.pressureOf(rhoLiquid) > model.pressureOf(rhoGas))
      muHigh = mu;
    else
      muLow = mu;
  }
  double mu = 0.5 * (muLow + muHigh);
  bulk.gasDensity = densityAtMu(model, mu, 1.0e-18, lowSpinodal);
  bulk.liquidDensity = densityAtMu(model, mu, highSpinodal, rhoMax);
  bulk.saturationPressure = model.pressureOf(bulk.gasDensity) * pressureKA3ToPascal;

  // Fused-dumbbell CS with the Ravikovitch WCA tail is only barely subcritical (P0 ~ 10 atm).
  // GCMC and the report Henry use experimental nitrogen P0; so does x when the model P0 is not
  // atmosphere-like. μ(ρ) and the liquid density stay those of the same functional.
  const bool dumbbell = !options.dumbbellSiteOffsets.empty() && options.dumbbellSiteDiameter > 0.0;
  if ((dumbbell || siteAttraction(options)) &&
      (bulk.saturationPressure < 5.0e4 || bulk.saturationPressure > 2.0e5))
  {
    bulk.experimentalSaturation = true;
    bulk.saturationPressure = nitrogenSaturationPressure;
    bulk.gasDensity = gasDensityAtPressure(model, bulk.saturationPressure / pressureKA3ToPascal, lowSpinodal);
  }
  return bulk;
}


double nldftGasDensity(const NLDFTOptions &options, const NLDFTBulk &bulk, double x)
{
  BulkModel model = bulkModelOf(options);
  double lowSpinodal = 0.0;
  double highSpinodal = 0.0;
  const bool hasLoop = findSpinodals(model, lowSpinodal, highSpinodal);
  if (!(bulk.saturationPressure > 0.0)) return 0.0;
  const double rhoHigh = hasLoop ? lowSpinodal : std::max(bulk.liquidDensity, 1.0e-6);
  double pressure = x * bulk.saturationPressure / pressureKA3ToPascal;
  return gasDensityAtPressure(model, pressure, rhoHigh);
}


void nldftApplyHenryPressureWindow(NLDFTOptions &options, double henryMoleculesPerCellAtExperimentalP0,
                                   double cellVolumeAngstrom3, double isothermSaturationPressurePa)
{
  // K_H [molecules / cell / Pa] from the reported occupancy at experimental P0.
  const double henryCoefficient =
      (henryMoleculesPerCellAtExperimentalP0 > 0.0 && nitrogenSaturationPressure > 0.0)
          ? henryMoleculesPerCellAtExperimentalP0 / nitrogenSaturationPressure
          : 0.0;
  const double packingCapacity =
      (cellVolumeAngstrom3 > 0.0) ? cellVolumeAngstrom3 / nitrogenLiquidVolume : 0.0;

  double lowestPressure = nldftPressureWindowFloorPa;
  if (henryCoefficient > 0.0 && packingCapacity > 0.0)
  {
    lowestPressure =
        std::min(nldftPressureWindowFloorPa, nldftPressureWindowFillFraction * packingCapacity / henryCoefficient);
  }

  const double p0 =
      (isothermSaturationPressurePa > 0.0) ? isothermSaturationPressurePa : nitrogenSaturationPressure;
  double highestPressure = std::max(p0, lowestPressure);
  // Keep at least a minimal log span so a single-point request still has distinct endpoints when
  // the caller later densifies the grid; WHAM uses a factor-of-five floor between ladder ends.
  if (highestPressure < lowestPressure * 5.0) highestPressure = lowestPressure * 5.0;

  options.xLow = lowestPressure / p0;
  options.xHigh = highestPressure / p0;
  if (options.xHigh < options.xLow) std::swap(options.xLow, options.xHigh);
}


NLDFTGridIsotherm nldftIsothermOnGrid(std::span<const float> energyKelvin, uint3 gridSize, const UnitCell &cell,
                                      const NLDFTOptions &options, std::span<const float> orientationKelvin,
                                      std::size_t nOrientations)
{
  NLDFTGridIsotherm result;
  result.bulk = nldftBulkCoexistence(options);
  if (!(result.bulk.saturationPressure > 0.0) || energyKelvin.empty()) return result;

  std::size_t nx = gridSize.x;
  std::size_t ny = gridSize.y;
  std::size_t nz = gridSize.z;
  if (energyKelvin.size() != nx * ny * nz) return result;
  const bool molecular = nOrientations > 1 && orientationKelvin.size() == energyKelvin.size() * nOrientations;
  result.molecular = molecular;
  result.numberOfOrientations = molecular ? nOrientations : 1;
  const std::size_t M = result.numberOfOrientations;

  // Decimate to the spacing the density profile actually needs; every FFT in every iteration pays for
  // each voxel kept. Each coarse voxel takes the free energy -kT ln <exp(-U/kT)> of its block rather than
  // a point sample: a channel minimum narrower than the coarse spacing then keeps its exact Boltzmann
  // weight (the Henry integral is conserved), it is only drawn a little wider.
  auto strideFor = [&](std::size_t n, double length)
  {
    if (options.targetSpacing <= 0.0 || n == 0) return std::size_t{1};
    auto limit = static_cast<std::size_t>(options.targetSpacing * static_cast<double>(n) / length);
    for (std::size_t s = std::max<std::size_t>(limit, 1); s > 1; --s)
    {
      if (n % s == 0) return s;
    }
    return std::size_t{1};
  };
  const std::size_t strideX = strideFor(nx, cell.lengthA);
  const std::size_t strideY = strideFor(ny, cell.lengthB);
  const std::size_t strideZ = strideFor(nz, cell.lengthC);
  std::vector<float> decimated;
  std::span<const float> energyView = energyKelvin;
  if (strideX * strideY * strideZ > 1)
  {
    std::size_t fullX = nx;
    std::size_t fullY = ny;
    nx /= strideX;
    ny /= strideY;
    nz /= strideZ;
    decimated.resize(nx * ny * nz);
    const double blockCount = static_cast<double>(strideX * strideY * strideZ);
    for (std::size_t k = 0; k < nz; ++k)
    {
      for (std::size_t j = 0; j < ny; ++j)
      {
        for (std::size_t i = 0; i < nx; ++i)
        {
          double weight = 0.0;
          for (std::size_t dk = 0; dk < strideZ; ++dk)
          {
            for (std::size_t dj = 0; dj < strideY; ++dj)
            {
              for (std::size_t di = 0; di < strideX; ++di)
              {
                double u = static_cast<double>(
                    energyKelvin[((k * strideZ + dk) * fullY + j * strideY + dj) * fullX + i * strideX + di]);
                weight += std::exp(std::max(-u / options.temperature, -700.0));
              }
            }
          }
          decimated[(k * ny + j) * nx + i] =
              static_cast<float>(-options.temperature * std::log(std::max(weight / blockCount, 1.0e-300)));
        }
      }
    }
    energyView = decimated;
  }

  std::vector<float> decimatedOrient;
  std::span<const float> orientView = orientationKelvin;
  if (molecular && strideX * strideY * strideZ > 1)
  {
    const std::size_t fullX = gridSize.x;
    const std::size_t fullY = gridSize.y;
    decimatedOrient.resize(nx * ny * nz * M);
    const double blockCount = static_cast<double>(strideX * strideY * strideZ);
    for (std::size_t k = 0; k < nz; ++k)
    {
      for (std::size_t j = 0; j < ny; ++j)
      {
        for (std::size_t i = 0; i < nx; ++i)
        {
          for (std::size_t o = 0; o < M; ++o)
          {
            double weight = 0.0;
            for (std::size_t dk = 0; dk < strideZ; ++dk)
            {
              for (std::size_t dj = 0; dj < strideY; ++dj)
              {
                for (std::size_t di = 0; di < strideX; ++di)
                {
                  const std::size_t vox =
                      ((k * strideZ + dk) * fullY + j * strideY + dj) * fullX + i * strideX + di;
                  double u = static_cast<double>(orientationKelvin[vox * nOrientations + o]);
                  weight += std::exp(std::max(-u / options.temperature, -700.0));
                }
              }
            }
            decimatedOrient[((k * ny + j) * nx + i) * M + o] =
                static_cast<float>(-options.temperature * std::log(std::max(weight / blockCount, 1.0e-300)));
          }
        }
      }
    }
    orientView = decimatedOrient;
  }

  const std::size_t N = nx * ny * nz;
  const std::size_t nComplexX = nx / 2 + 1;
  const std::size_t nC = nz * ny * nComplexX;

  const double T = options.temperature;
  const double beta = 1.0 / T;
  const bool dumbbell = molecular && !options.dumbbellSiteOffsets.empty() && options.dumbbellSiteDiameter > 0.0;
  const bool siteAtt = molecular && siteAttraction(options);
  const double d = dumbbell ? options.dumbbellSiteDiameter : options.hardSphereDiameter;
  const double R = 0.5 * d;
  const double dV = cell.volume / static_cast<double>(N);
  const double s0 = 1.0 / (std::numbers::pi * d * d);   // n0 = s0 n2
  const double s1 = 1.0 / (2.0 * std::numbers::pi * d);  // n1 = s1 n2

  BulkModel model = bulkModelOf(options);
  double lowSpinodal = 0.0;
  double highSpinodal = 0.0;
  const bool hasLoop = findSpinodals(model, lowSpinodal, highSpinodal);
  const double rhoGasCap = hasLoop ? lowSpinodal : std::max(result.bulk.liquidDensity, 1.0e-6);
  // Contact peaks on a coarse voxel can be several times the bulk liquid (a 1-D file is a few voxels
  // across). They cannot be hundreds of times: rhoPeakCap = 4 /Å³ was 200× liquid and let a single
  // well voxel hold a molecule, which is how TraPPE N2 reached 280 / cell. The cell-mean density
  // cannot: eta = rho * pi d³ / 6 lives below 1, and 0.85 is already past any physical pore fill.
  const double rhoPeakCap = 8.0 * std::max(result.bulk.liquidDensity, 0.017);
  const double rhoMeanCap = 0.85 / (model.sites * model.packingFactor);

  std::call_once(fftwThreadsFlag,
                 []
                 {
                   fftw_init_threads();
                   unsigned threads = std::max(1u, std::min(8u, std::thread::hardware_concurrency()));
                   fftw_plan_with_nthreads(static_cast<int>(threads));
                 });

  // Everything the transforms touch comes from fftw_malloc so one plan serves every array.
  auto realFree = [](double *p) { fftw_free(p); };
  auto complexFree = [](fftw_complex *p) { fftw_free(p); };
  using RealArray = std::unique_ptr<double[], decltype(realFree)>;
  using ComplexArray = std::unique_ptr<fftw_complex[], decltype(complexFree)>;
  auto makeReal = [&]() { return RealArray(fftw_alloc_real(N), realFree); };
  auto makeComplex = [&]() { return ComplexArray(fftw_alloc_complex(nC), complexFree); };

  RealArray rho = makeReal();
  RealArray rhoSite = makeReal();
  RealArray n3B = makeReal();     // n3 in, dPhi/dn3 out
  RealArray n2A = makeReal();     // n2 in, combined scalar derivative out
  RealArray vx = makeReal();      // n2v in, combined vector derivative out
  RealArray vy = makeReal();
  RealArray vz = makeReal();
  RealArray phiAtt = makeReal();
  RealArray muHS = makeReal();
  ComplexArray rhoHat = makeComplex();
  ComplexArray scratch = makeComplex();
  ComplexArray acc = makeComplex();
  if (siteAtt)
  {
    for (std::size_t i = 0; i < N; ++i) phiAtt[i] = 0.0;
  }

  fftw_plan forward = fftw_plan_dft_r2c_3d(static_cast<int>(nz), static_cast<int>(ny), static_cast<int>(nx),
                                           rho.get(), rhoHat.get(), FFTW_MEASURE);
  fftw_plan inverse = fftw_plan_dft_c2r_3d(static_cast<int>(nz), static_cast<int>(ny), static_cast<int>(nx),
                                           scratch.get(), muHS.get(), FFTW_MEASURE);

  // The weight spectra: radial transforms evaluated on the reciprocal lattice of whatever cell this is.
  std::vector<double> w2Hat(nC), w3Hat(nC), uHat(nC), kX(nC), kY(nC), kZ(nC), kAbs(nC);
  double kMax = 0.0;
  for (std::size_t iz = 0; iz < nz; ++iz)
  {
    double mz = (iz <= nz / 2) ? static_cast<double>(iz) : static_cast<double>(iz) - static_cast<double>(nz);
    for (std::size_t iy = 0; iy < ny; ++iy)
    {
      double my = (iy <= ny / 2) ? static_cast<double>(iy) : static_cast<double>(iy) - static_cast<double>(ny);
      for (std::size_t ix = 0; ix < nComplexX; ++ix)
      {
        std::size_t index = (iz * ny + iy) * nComplexX + ix;
        double3 m(static_cast<double>(ix), my, mz);
        double3 k = 2.0 * std::numbers::pi * transposedMultiply(cell.inverseCell, m);
        kX[index] = k.x;
        kY[index] = k.y;
        kZ[index] = k.z;
        double kk = std::sqrt(k.x * k.x + k.y * k.y + k.z * k.z);
        kAbs[index] = kk;
        kMax = std::max(kMax, kk);
        if (kk < 1.0e-12)
        {
          w3Hat[index] = 4.0 * std::numbers::pi * R * R * R / 3.0;
          w2Hat[index] = 4.0 * std::numbers::pi * R * R;
        }
        else
        {
          double kR = kk * R;
          w3Hat[index] = 4.0 * std::numbers::pi * (std::sin(kR) - kR * std::cos(kR)) / (kk * kk * kk);
          w2Hat[index] = 4.0 * std::numbers::pi * R * std::sin(kR) / kk;
        }
      }
    }
  }

  // The pair's radial transform, tabulated once over |k| and read back per lattice point:
  // u(k) = (4 pi / k) int_0^rc u(r) r sin(kr) dr, Simpson. Site-site TraPPE uses the 12-6 from σ;
  // spherical N2 keeps the Ravikovitch WCA tail. u(k=0) is the pair integral (a_NN), not n² a_NN.
  {
    const bool ljSites = siteAtt;
    const double pairSigma = ljSites ? options.attractiveSigma : options.sigma;
    const double pairEpsilon = ljSites ? options.attractiveEpsilon : options.epsilon;
    const double rMin = std::pow(2.0, 1.0 / 6.0) * pairSigma;
    const double nAtt =
        ljSites ? static_cast<double>(options.attractiveSiteOffsets.size()) : 1.0;
    const std::size_t tableSize = 8192;
    const std::size_t nr = 4000;  // even
    const double dr = options.cutoff / static_cast<double>(nr);
    std::vector<double> uTable(tableSize);
    uTable[0] = (nAtt > 0.0) ? model.a / (nAtt * nAtt) : model.a;
    for (std::size_t t = 1; t < tableSize; ++t)
    {
      double k = kMax * static_cast<double>(t) / static_cast<double>(tableSize - 1);
      double sum = 0.0;
      for (std::size_t i = 0; i <= nr; ++i)
      {
        double r = dr * static_cast<double>(i);
        double weight = (i == 0 || i == nr) ? 1.0 : ((i % 2 == 1) ? 4.0 : 2.0);
        const double u = ljSites ? fullLennardJones(r, pairSigma, pairEpsilon, options.cutoff)
                                 : attractiveTail(r, pairSigma, pairEpsilon, rMin, options.cutoff);
        sum += weight * u * r * std::sin(k * r);
      }
      uTable[t] = 4.0 * std::numbers::pi * (sum * dr / 3.0) / k;
    }
    for (std::size_t index = 0; index < nC; ++index)
    {
      double position = kAbs[index] / kMax * static_cast<double>(tableSize - 1);
      std::size_t bin = std::min(static_cast<std::size_t>(position), tableSize - 2);
      double fraction = position - static_cast<double>(bin);
      uHat[index] = uTable[bin] * (1.0 - fraction) + uTable[bin + 1] * fraction;
    }
  }

  // beta U for the Euler-Lagrange map. Overlap (βU > 690) is dropped. Spherical ρ clips wells at
  // wellFloorInKT: a Helmholtz PMF can be a Coulomb hole a spherical density cannot occupy. Molecular
  // ρ(r, ω) keeps the TraPPE well (FMT packing is the bound); clipping it at 12 kT made every μ see
  // the same e^12 contact layer and stalled FER at ~4 / cell against GCMC ~10. Henry is always unclipped.
  const double wellFloor = options.wellFloorInKT > 0.0 ? options.wellFloorInKT : nldftWellFloorInKT;
  std::vector<double> betaU(N * M);
  std::vector<char> accessible(N, 0);
  std::vector<std::size_t> liveStart(N + 1, 0);
  std::vector<double> liveLogWeight;
  std::vector<std::uint16_t> liveOrient;
  if (molecular)
  {
    liveLogWeight.reserve(N * M / 8);
    liveOrient.reserve(N * M / 8);
    for (std::size_t i = 0; i < N; ++i)
    {
      liveStart[i] = liveLogWeight.size();
      for (std::size_t o = 0; o < M; ++o)
      {
        const double bu = std::min(beta * static_cast<double>(orientView[i * M + o]), 700.0);
        betaU[i * M + o] = bu;
        if (bu < 699.0) accessible[i] = 1;
        if (bu > 690.0) continue;
        liveLogWeight.push_back(-bu);
        liveOrient.push_back(static_cast<std::uint16_t>(o));
      }
    }
    liveStart[N] = liveLogWeight.size();
  }
  else
  {
    for (std::size_t i = 0; i < N; ++i)
    {
      const double bu = std::min(beta * static_cast<double>(energyView[i]), 700.0);
      betaU[i] = std::max(bu, -wellFloor);
      accessible[i] = bu < 699.0 ? 1 : 0;
    }
  }

  const bool needDeposit = dumbbell || siteAtt;
  const std::vector<double> &geomOffsets =
      siteAtt ? options.attractiveSiteOffsets : options.dumbbellSiteOffsets;
  const std::size_t nSites = geomOffsets.size();
  const double siteScale =
      (dumbbell && nSites > 0) ? model.sites / static_cast<double>(nSites) : 1.0;
  const double depositScale = dumbbell ? siteScale : 1.0;
  std::vector<double3> siteDelta;
  if (needDeposit && nSites > 0)
  {
    siteDelta.resize(M * nSites);
    const std::vector<double3> axes = orientationSet(M, options.headTailSymmetric);
    for (std::size_t o = 0; o < M; ++o)
    {
      for (std::size_t s = 0; s < nSites; ++s)
        siteDelta[o * nSites + s] = cell.inverseCell * (geomOffsets[s] * axes[o]);
    }
  }
  for (std::size_t i = 0; i < N; ++i) muHS[i] = 0.0;
  RealArray muHSPrev = makeReal();
  for (std::size_t i = 0; i < N; ++i) muHSPrev[i] = 0.0;
  bool haveMuHS = false;

  auto sampleSites = [&](const double *grid, std::size_t o, double fx, double fy, double fz) -> double
  {
    double sum = 0.0;
    for (std::size_t s = 0; s < nSites; ++s)
    {
      const double3 df = siteDelta[o * nSites + s];
      sum += nldftCicSample(grid, nx, ny, nz, fx + df.x, fy + df.y, fz + df.z);
    }
    return sum;
  };
  auto dumbbellMu = [&](std::size_t o, double fx, double fy, double fz) -> double
  { return siteScale * sampleSites(muHS.get(), o, fx, fy, fz); };

  const unsigned workers = nldftWorkerCount();
  std::vector<std::vector<double>> siteLocal;
  if (needDeposit)
  {
    siteLocal.resize(workers);
    for (auto &local : siteLocal) local.assign(N, 0.0);
  }

  const double inverseN = 1.0 / static_cast<double>(N);
  auto convolveTo = [&](const std::vector<double> &spectrum, double *out)
  {
    for (std::size_t index = 0; index < nC; ++index)
    {
      scratch[index][0] = rhoHat[index][0] * spectrum[index] * inverseN;
      scratch[index][1] = rhoHat[index][1] * spectrum[index] * inverseN;
    }
    fftw_execute_dft_c2r(inverse, scratch.get(), out);
  };
  // n2v_j = IFFT( -i k_j w3 rho ): (a + b i)(-i c) = b c - a c i.
  auto vectorWeightTo = [&](const std::vector<double> &kComponent, double *out)
  {
    for (std::size_t index = 0; index < nC; ++index)
    {
      double c = w3Hat[index] * kComponent[index] * inverseN;
      scratch[index][0] = rhoHat[index][1] * c;
      scratch[index][1] = -rhoHat[index][0] * c;
    }
    fftw_execute_dft_c2r(inverse, scratch.get(), out);
  };

  auto eulerLagrange = [&](double muOverT, std::vector<double> &updated) -> double
  {
    fftw_execute_dft_r2c(forward, rho.get(), rhoHat.get());
    if (!siteAtt) convolveTo(uHat, phiAtt.get());
    if (needDeposit)
    {
      std::fill(rhoSite.get(), rhoSite.get() + static_cast<std::ptrdiff_t>(N), 0.0);
      for (auto &local : siteLocal) std::fill(local.begin(), local.end(), 0.0);
      const std::size_t depositChunk = (N + workers - 1) / std::max(workers, 1u);
      nldftParallelFor(N, workers, [&](std::size_t begin, std::size_t end)
      {
        const unsigned thread = static_cast<unsigned>(
            std::min(begin / std::max(depositChunk, std::size_t{1}), static_cast<std::size_t>(workers - 1)));
        double *local = siteLocal[thread].data();
        for (std::size_t i = begin; i < end; ++i)
        {
          if (rho[i] < 1.0e-18) continue;
          const std::size_t beginO = liveStart[i];
          const std::size_t endO = liveStart[i + 1];
          if (beginO >= endO) continue;
          const std::size_t iz = i / (nx * ny);
          const std::size_t iy = (i / nx) % ny;
          const std::size_t ix = i % nx;
          const double fx = (static_cast<double>(ix) + 0.5) / static_cast<double>(nx);
          const double fy = (static_cast<double>(iy) + 0.5) / static_cast<double>(ny);
          const double fz = (static_cast<double>(iz) + 0.5) / static_cast<double>(nz);
          const std::size_t nLive = endO - beginO;
          std::vector<double> logw(nLive);
          double maxLog = -1.0e300;
          for (std::size_t k = beginO; k < endO; ++k)
          {
            const std::size_t o = liveOrient[k];
            double extra = 0.0;
            if (dumbbell) extra += dumbbellMu(o, fx, fy, fz);
            if (siteAtt) extra += beta * sampleSites(phiAtt.get(), o, fx, fy, fz);
            logw[k - beginO] = liveLogWeight[k] - extra;
            maxLog = std::max(maxLog, logw[k - beginO]);
          }
          double Z = 0.0;
          for (double lw : logw) Z += std::exp(std::min(lw - maxLog, 80.0));
          const double invZ = (Z > 0.0) ? 1.0 / Z : 0.0;
          for (std::size_t k = beginO; k < endO; ++k)
          {
            const double occ = rho[i] * std::exp(std::min(logw[k - beginO] - maxLog, 80.0)) * invZ;
            if (!(occ > 0.0)) continue;
            const std::size_t o = liveOrient[k];
            for (std::size_t s = 0; s < nSites; ++s)
            {
              const double3 df = siteDelta[o * nSites + s];
              nldftCicAdd(local, nx, ny, nz, fx + df.x, fy + df.y, fz + df.z, depositScale * occ);
            }
          }
        }
      });
      for (unsigned t = 0; t < workers; ++t)
        for (std::size_t i = 0; i < N; ++i) rhoSite[i] += siteLocal[t][i];
      if (dumbbell) fftw_execute_dft_r2c(forward, rhoSite.get(), rhoHat.get());
      if (siteAtt && dumbbell) convolveTo(uHat, phiAtt.get());
    }
    convolveTo(w3Hat, n3B.get());
    convolveTo(w2Hat, n2A.get());
    vectorWeightTo(kX, vx.get());
    vectorWeightTo(kY, vy.get());
    vectorWeightTo(kZ, vz.get());

    nldftParallelFor(N, workers, [&](std::size_t begin, std::size_t end)
    {
      for (std::size_t i = begin; i < end; ++i)
      {
        double n3 = std::clamp(n3B[i], 0.0, 0.999);
        double n2 = std::max(n2A[i], 0.0);
        double wx = vx[i];
        double wy = vy[i];
        double wz = vz[i];
        double vSquared = wx * wx + wy * wy + wz * wz;
        double n2Squared = n2 * n2;
        if (vSquared > 0.99 * n2Squared)
        {
          double scale = (n2Squared > 0.0) ? std::sqrt(0.99 * n2Squared / vSquared) : 0.0;
          wx *= scale;
          wy *= scale;
          wz *= scale;
          vSquared = 0.99 * n2Squared;
        }
        double o = 1.0 / (1.0 - n3);
        double phi3 = 0.0;
        double phi3Prime = 0.0;
        phi3Of(n3, phi3, phi3Prime);

        double n0 = s0 * n2;
        double n1 = s1 * n2;
        double dPhi0 = -std::log(1.0 - n3);
        double dPhi1 = n2 * o;
        double dPhi2 = n1 * o + 3.0 * (n2Squared - vSquared) * phi3;
        double dPhi3 = n0 * o + (n1 * n2 - s1 * vSquared) * o * o +
                       (n2 * n2Squared - 3.0 * n2 * vSquared) * phi3Prime;

        n2A[i] = s0 * dPhi0 + s1 * dPhi1 + dPhi2;
        n3B[i] = dPhi3;
        double vectorFactor = -2.0 * s1 * o - 6.0 * n2 * phi3;
        vx[i] = vectorFactor * wx;
        vy[i] = vectorFactor * wy;
        vz[i] = vectorFactor * wz;
      }
    });

    fftw_execute_dft_r2c(forward, n2A.get(), scratch.get());
    for (std::size_t index = 0; index < nC; ++index)
    {
      acc[index][0] = scratch[index][0] * w2Hat[index];
      acc[index][1] = scratch[index][1] * w2Hat[index];
    }
    fftw_execute_dft_r2c(forward, n3B.get(), scratch.get());
    for (std::size_t index = 0; index < nC; ++index)
    {
      acc[index][0] += scratch[index][0] * w3Hat[index];
      acc[index][1] += scratch[index][1] * w3Hat[index];
    }
    const std::vector<double> *kComponents[3] = {&kX, &kY, &kZ};
    double *vectorFields[3] = {vx.get(), vy.get(), vz.get()};
    for (int component = 0; component < 3; ++component)
    {
      fftw_execute_dft_r2c(forward, vectorFields[component], scratch.get());
      const std::vector<double> &kc = *kComponents[component];
      for (std::size_t index = 0; index < nC; ++index)
      {
        double c = w3Hat[index] * kc[index];
        acc[index][0] += -scratch[index][1] * c;
        acc[index][1] += scratch[index][0] * c;
      }
    }
    for (std::size_t index = 0; index < nC; ++index)
    {
      scratch[index][0] = acc[index][0] * inverseN;
      scratch[index][1] = acc[index][1] * inverseN;
    }
    fftw_execute_dft_c2r(inverse, scratch.get(), muHS.get());
    if (dumbbell)
    {
      if (haveMuHS)
      {
        for (std::size_t i = 0; i < N; ++i) muHS[i] = 0.5 * muHS[i] + 0.5 * muHSPrev[i];
      }
      for (std::size_t i = 0; i < N; ++i) muHSPrev[i] = muHS[i];
      haveMuHS = true;
    }
    if (siteAtt && !dumbbell)
    {
      fftw_execute_dft_r2c(forward, rhoSite.get(), rhoHat.get());
      convolveTo(uHat, phiAtt.get());
    }

    constexpr double tiny = 1.0e-300;
    constexpr double logTiny = -690.0;
    std::vector<double> threadSumSq(workers, 0.0);
    const std::size_t chunk = (N + workers - 1) / std::max(workers, 1u);
    nldftParallelFor(N, workers, [&](std::size_t begin, std::size_t end)
    {
      double sumSq = 0.0;
      for (std::size_t i = begin; i < end; ++i)
      {
        const double shiftAtt = (siteAtt && molecular) ? muOverT : muOverT - beta * phiAtt[i];
        double value = tiny;
        if (molecular)
        {
          const std::size_t beginO = liveStart[i];
          const std::size_t endO = liveStart[i + 1];
          if (beginO < endO)
          {
            const std::size_t iz = i / (nx * ny);
            const std::size_t iy = (i / nx) % ny;
            const std::size_t ix = i % nx;
            const double fx = (static_cast<double>(ix) + 0.5) / static_cast<double>(nx);
            const double fy = (static_cast<double>(iy) + 0.5) / static_cast<double>(ny);
            const double fz = (static_cast<double>(iz) + 0.5) / static_cast<double>(nz);
            auto extraOf = [&](std::size_t o) -> double
            {
              double extra = dumbbell ? dumbbellMu(o, fx, fy, fz) : muHS[i];
              if (siteAtt) extra += beta * sampleSites(phiAtt.get(), o, fx, fy, fz);
              return extra;
            };
            double maxLog = -1.0e300;
            for (std::size_t k = beginO; k < endO; ++k)
              maxLog = std::max(maxLog, shiftAtt + liveLogWeight[k] - extraOf(liveOrient[k]));
            double sum = 0.0;
            for (std::size_t k = beginO; k < endO; ++k)
            {
              const double extra = extraOf(liveOrient[k]);
              sum += std::exp(std::min(shiftAtt + liveLogWeight[k] - extra - maxLog, 80.0));
            }
            const double logVal = maxLog + std::log(std::max(sum, tiny)) - std::log(static_cast<double>(M));
            value = (logVal >= 80.0) ? rhoPeakCap : std::max(std::exp(logVal), tiny);
          }
        }
        else
        {
          double exponent = shiftAtt - muHS[i] - betaU[i];
          value = (exponent < logTiny) ? tiny : std::exp(std::min(exponent, 80.0));
        }
        value = std::min(value, rhoPeakCap);
        const double diff = value - rho[i];
        sumSq += diff * diff;
        updated[i] = value;
      }
      const unsigned thread = static_cast<unsigned>(std::min(begin / std::max(chunk, std::size_t{1}),
                                                            static_cast<std::size_t>(workers - 1)));
      threadSumSq[thread] = sumSq;
    });
    double sumSq = 0.0;
    for (unsigned t = 0; t < workers; ++t) sumSq += threadSumSq[t];
    return std::sqrt(sumSq / static_cast<double>(N));
  };

  auto meanDensity = [&]()
  {
    double sum = 0.0;
    for (std::size_t i = 0; i < N; ++i) sum += rho[i];
    return sum / static_cast<double>(N);
  };
  auto loadingOf = [&]()
  {
    double sum = 0.0;
    for (std::size_t i = 0; i < N; ++i) sum += rho[i];
    return sum * dV;
  };

  std::deque<double> targets;
  {
    const std::size_t points = static_cast<std::size_t>(std::max(options.pressurePoints, 1));
    double logMin = std::log(options.xLow);
    double logMax = std::log(options.xHigh);
    for (std::size_t i = 0; i < points; ++i)
    {
      double fraction = (points > 1) ? static_cast<double>(i) / static_cast<double>(points - 1) : 0.0;
      targets.push_back(std::exp(logMin + (logMax - logMin) * fraction));
    }
  }
  const double minLogStep = (std::log(options.xHigh) - std::log(options.xLow)) /
                            static_cast<double>(4 * std::max(options.pressurePoints, 1));

  std::vector<double> updated(N);
  std::vector<double> lastConverged(N);
  std::vector<double> trial(N);

  bool haveProfile = false;
  bool hasPrev = false;
  double xPrev = 0.0;
  double nPrev = 0.0;
  std::size_t extraPoints = 0;
  const std::size_t maxExtra = static_cast<std::size_t>(std::max(options.pressurePoints, 1));
  result.isotherm.reserve(targets.size() + maxExtra);

  auto orientationBoltzmann = [&](std::size_t i, double cap) -> double
  {
    if (!molecular) return std::exp(std::min(-betaU[i], cap));
    double sum = 0.0;
    for (std::size_t k = liveStart[i]; k < liveStart[i + 1]; ++k)
      sum += std::exp(std::min(liveLogWeight[k], cap));
    return sum / static_cast<double>(M);
  };

  auto seedDensity = [&](double rhoBulk)
  {
    for (std::size_t i = 0; i < N; ++i)
    {
      if (options.seedLiquid)
      {
        rho[i] = std::min(result.bulk.liquidDensity * orientationBoltzmann(i, 0.0), rhoPeakCap);
      }
      else
      {
        rho[i] = std::min(rhoBulk * orientationBoltzmann(i, 8.0), rhoPeakCap);
      }
    }
  };

  while (!targets.empty())
  {
    double x = targets.front();
    targets.pop_front();
    double pressure = x * result.bulk.saturationPressure / pressureKA3ToPascal;
    double rhoBulk = gasDensityAtPressure(model, pressure, rhoGasCap);
    double muOverT = model.muOf(rhoBulk) / T;

    if (!haveProfile)
      seedDensity(rhoBulk);
    else
      for (std::size_t i = 0; i < N; ++i) rho[i] = lastConverged[i];

    double mixing = options.mixing;
    double acceptedResidual = std::numeric_limits<double>::max();
    bool converged = false;
    for (std::size_t iteration = 0; iteration < options.maxIterations; ++iteration)
    {
      double residual = eulerLagrange(muOverT, updated);
      double densityScale = std::max(result.bulk.liquidDensity, meanDensity());
      double nUpdated = 0.0;
      for (std::size_t i = 0; i < N; ++i) nUpdated += updated[i];
      nUpdated *= dV;
      const double nNow = loadingOf();
      const double nScale = std::max(std::max(nNow, nUpdated), 1.0e-12);
      if (residual < options.tolerance * densityScale || std::fabs(nUpdated - nNow) < 1.0e-3 * nScale)
      {
        converged = true;
        break;
      }

      double currentMean = meanDensity();
      if (residual < acceptedResidual)
      {
        acceptedResidual = residual;
        mixing = std::min(mixing * 1.05, 0.5);
      }
      else
      {
        mixing = std::max(mixing * 0.5, 0.02);
      }

      double updatedMean = 0.0;
      for (std::size_t i = 0; i < N; ++i) updatedMean += updated[i];
      updatedMean /= static_cast<double>(N);
      if (updatedMean > 1.5 * currentMean)
        mixing = std::min(std::max(mixing, 0.2), 0.5);

      double trialMean = 0.0;
      double accessibleMean = 0.0;
      std::size_t accessibleCount = 0;
      for (std::size_t i = 0; i < N; ++i)
      {
        trial[i] = std::min((1.0 - mixing) * rho[i] + mixing * updated[i], rhoPeakCap);
        if (trial[i] < 0.0) trial[i] = 0.0;
        trialMean += trial[i];
        if (accessible[i])
        {
          accessibleMean += trial[i];
          ++accessibleCount;
        }
      }
      trialMean /= static_cast<double>(N);
      if (accessibleCount > 0)
      {
        accessibleMean /= static_cast<double>(accessibleCount);
        if (accessibleMean > rhoMeanCap)
        {
          double scale = rhoMeanCap / accessibleMean;
          for (std::size_t i = 0; i < N; ++i)
          {
            if (accessible[i]) trial[i] *= scale;
          }
        }
      }
      else if (trialMean > rhoMeanCap && trialMean > 0.0)
      {
        double scale = rhoMeanCap / trialMean;
        for (std::size_t i = 0; i < N; ++i) trial[i] *= scale;
      }
      for (std::size_t i = 0; i < N; ++i) rho[i] = trial[i];
    }

    double loading = loadingOf();
    if (!converged)
    {
      ++result.unconverged;
      if (haveProfile)
      {
        for (std::size_t i = 0; i < N; ++i) rho[i] = lastConverged[i];
        loading = loadingOf();
      }
    }

    double fillStep = 0.2 * result.bulk.liquidDensity * cell.volume;
    if (hasPrev && loading > nPrev + fillStep && extraPoints < maxExtra)
    {
      double logStep = std::log(x) - std::log(xPrev);
      if (logStep > minLogStep && minLogStep > 0.0)
      {
        targets.push_front(x);
        targets.push_front(std::exp(0.5 * (std::log(x) + std::log(xPrev))));
        ++extraPoints;
        if (haveProfile)
          for (std::size_t i = 0; i < N; ++i) rho[i] = lastConverged[i];
        continue;
      }
    }

    // Adsorption sweep: a capillary fill that Picard then desorbs is kept. FER reached n ≈ 10 (the
    // GCMC plateau) and the next x fell back to the four-molecule single-file branch.
    bool heldAdsorption = false;
    if (hasPrev && haveProfile && nPrev - loading > std::max(0.2 * nPrev, 0.05))
    {
      for (std::size_t i = 0; i < N; ++i) rho[i] = lastConverged[i];
      loading = nPrev;
      heldAdsorption = true;
    }

    if (!heldAdsorption && converged && meanDensity() <= rhoMeanCap * 1.01)
    {
      for (std::size_t i = 0; i < N; ++i) lastConverged[i] = rho[i];
      haveProfile = true;
    }

    result.isotherm.push_back(IsothermPoint{x, loading});
    xPrev = x;
    nPrev = loading;
    hasPrev = true;
  }

  // Desorption / equilibrium branch. GCMC samples the equilibrium state; the adsorption sweep of
  // mean-field DFT delays 10-ring fill until the pore spinodal (FER: n ≈ 4 until x ≈ 0.27, then 11).
  // Walking down from the filled profile keeps that loading into the BET window, which is the
  // comparison GCMC BET 250–350 m²/g is making.
  if (haveProfile && result.isotherm.size() > 4 &&
      result.isotherm.back().moleculesPerCell >
          1.5 * result.isotherm.front().moleculesPerCell + 1.0)
  {
    std::vector<double> xDown;
    xDown.reserve(result.isotherm.size());
    for (auto it = result.isotherm.rbegin(); it != result.isotherm.rend(); ++it)
      xDown.push_back(it->relativePressure);
    std::vector<IsothermPoint> desorption;
    desorption.reserve(xDown.size());
    desorption.push_back(result.isotherm.back());
    for (std::size_t k = 1; k < xDown.size(); ++k)
    {
      const double x = xDown[k];
      const double pressure = x * result.bulk.saturationPressure / pressureKA3ToPascal;
      const double rhoBulk = gasDensityAtPressure(model, pressure, rhoGasCap);
      const double muOverT = model.muOf(rhoBulk) / T;
      for (std::size_t i = 0; i < N; ++i) rho[i] = lastConverged[i];

      double mixing = options.mixing;
      double acceptedResidual = std::numeric_limits<double>::max();
      bool converged = false;
      for (std::size_t iteration = 0; iteration < options.maxIterations; ++iteration)
      {
        double residual = eulerLagrange(muOverT, updated);
        double densityScale = std::max(result.bulk.liquidDensity, meanDensity());
        double nUpdated = 0.0;
        for (std::size_t i = 0; i < N; ++i) nUpdated += updated[i];
        nUpdated *= dV;
        const double nNow = loadingOf();
        const double nScale = std::max(std::max(nNow, nUpdated), 1.0e-12);
        if (residual < options.tolerance * densityScale || std::fabs(nUpdated - nNow) < 1.0e-3 * nScale)
        {
          converged = true;
          break;
        }
        if (residual < acceptedResidual)
        {
          acceptedResidual = residual;
          mixing = std::min(mixing * 1.05, 0.5);
        }
        else
        {
          mixing = std::max(mixing * 0.5, 0.02);
        }
        for (std::size_t i = 0; i < N; ++i)
        {
          trial[i] = std::min((1.0 - mixing) * rho[i] + mixing * updated[i], rhoPeakCap);
          if (trial[i] < 0.0) trial[i] = 0.0;
        }
        for (std::size_t i = 0; i < N; ++i) rho[i] = trial[i];
      }
      double loading = loadingOf();
      if (!converged)
      {
        ++result.unconverged;
        for (std::size_t i = 0; i < N; ++i) rho[i] = lastConverged[i];
        loading = loadingOf();
      }
      else if (meanDensity() <= rhoMeanCap * 1.01)
      {
        for (std::size_t i = 0; i < N; ++i) lastConverged[i] = rho[i];
      }
      desorption.push_back(IsothermPoint{x, loading});
    }
    std::reverse(desorption.begin(), desorption.end());
    result.isotherm = std::move(desorption);
  }

  fftw_destroy_plan(forward);
  fftw_destroy_plan(inverse);

  return result;
}


NLDFTExternalFieldStats nldftMaskAndClipExternalPotential(std::span<float> energyKelvin, double temperature,
                                                         float ceilingKelvin, double wellFloorInKT)
{
  NLDFTExternalFieldStats stats;
  if (energyKelvin.empty()) return stats;
  if (temperature <= 0.0) temperature = 77.355;
  if (!(wellFloorInKT > 0.0)) wellFloorInKT = nldftWellFloorInKT;
  stats.wellFloorKelvin = wellFloorInKT * temperature;
  const float floorKelvin = static_cast<float>(-stats.wellFloorKelvin);
  for (std::size_t i = 0; i < energyKelvin.size(); ++i)
  {
    if (energyKelvin[i] < floorKelvin)
    {
      energyKelvin[i] = floorKelvin;
      ++stats.clippedVoxels;
    }
    else if (energyKelvin[i] > ceilingKelvin)
    {
      energyKelvin[i] = ceilingKelvin;
    }
  }
  return stats;
}
