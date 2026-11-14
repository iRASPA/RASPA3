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

// Cloud-in-cell deposit and sample on the density grid. Grid point i stands at fractional i / n, which is
// where the energy field put its voxel too, so a site at the voxel's own position lands on that one point
// and is read back from it alone; a site a fraction of a voxel off is shared with the next point in each
// axis by that fraction. Asking at (i + 1/2) / n instead would put every site half a voxel off the field it
// is to feel, and split a site that is not off at all over the eight points around it.
void nldftCicDeposit(double *grid, std::size_t nx, std::size_t ny, std::size_t nz, double fx, double fy, double fz,
                     double weight)
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

void nldftCicAdd(double *grid, std::size_t nx, std::size_t ny, std::size_t nz, double fx, double fy, double fz,
                 double weight)
{
  if (!(weight > 0.0) || !(weight < 1.0e300)) return;
  nldftCicDeposit(grid, nx, ny, nz, fx, fy, fz, weight);
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
    // never zero: the window is spaced logarithmically, and an enormous (or infinite) Henry coefficient of an
    // unshielded well took the lower end to nought and the whole isotherm with its logarithm
    lowestPressure = std::max(
        std::min(nldftPressureWindowFloorPa, nldftPressureWindowFillFraction * packingCapacity / henryCoefficient),
        1.0e-30);
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
  // The block's free energy is worked out from its lowest energy down: u_min - kT ln <exp(-(u - u_min)/kT)>.
  // Every term of the average is then between nought and one and the lowest is one itself, so nothing
  // overflows and the average is never nought, however deep a well is and however high a wall. Summing
  // exp(-u/kT) as it stands overflows where a state lies more than 709 kT down, which is no depth a real
  // well has but is one an orientation energy with a charge in it is free to have, and the block came out at
  // minus infinity. It also made a block that is wall all through 691 kT high, where the wall itself is at
  // the ceiling.
  const std::size_t voxelsInBlock = strideX * strideY * strideZ;
  auto blockFreeEnergies = [&](std::span<const float> fine, std::size_t states)
  {
    const std::size_t fullX = gridSize.x;
    const std::size_t fullY = gridSize.y;
    std::vector<float> coarse(nx * ny * nz * states);
    const double blockCount = static_cast<double>(voxelsInBlock);
    auto eachOfBlock = [&](std::size_t i, std::size_t j, std::size_t k, std::size_t state, auto &&take)
    {
      for (std::size_t dk = 0; dk < strideZ; ++dk)
      {
        for (std::size_t dj = 0; dj < strideY; ++dj)
        {
          for (std::size_t di = 0; di < strideX; ++di)
          {
            const std::size_t voxel = ((k * strideZ + dk) * fullY + j * strideY + dj) * fullX + i * strideX + di;
            take(static_cast<double>(fine[voxel * states + state]));
          }
        }
      }
    };
    for (std::size_t k = 0; k < nz; ++k)
    {
      for (std::size_t j = 0; j < ny; ++j)
      {
        for (std::size_t i = 0; i < nx; ++i)
        {
          for (std::size_t state = 0; state < states; ++state)
          {
            double lowest = std::numeric_limits<double>::max();
            eachOfBlock(i, j, k, state, [&](double u) { lowest = std::min(lowest, u); });
            double weight = 0.0;
            eachOfBlock(i, j, k, state, [&](double u) { weight += std::exp(-(u - lowest) / options.temperature); });
            coarse[((k * ny + j) * nx + i) * states + state] =
                static_cast<float>(lowest - options.temperature * std::log(weight / blockCount));
          }
        }
      }
    }
    return coarse;
  };

  std::vector<float> decimated;
  std::span<const float> energyView = energyKelvin;
  std::vector<float> decimatedOrient;
  std::span<const float> orientView = orientationKelvin;
  if (voxelsInBlock > 1)
  {
    nx /= strideX;
    ny /= strideY;
    nz /= strideZ;
    decimated = blockFreeEnergies(energyKelvin, 1);
    energyView = decimated;
    if (molecular)
    {
      decimatedOrient = blockFreeEnergies(orientationKelvin, nOrientations);
      orientView = decimatedOrient;
    }
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
  // How close the packing may come to one before the functional gives up on it. A site holding one molecule has n3
  // exponentially close to one in the chemical potential, so this is well past any occupancy that shows in a loading.
  constexpr double largestPacking = 1.0 - 1.0e-6;

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

  std::vector<double> w2Hat(nC), w3Hat(nC), uHat(nC), vHatX(nC), vHatY(nC), vHatZ(nC), kAbs(nC);
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
        double kk = std::sqrt(k.x * k.x + k.y * k.y + k.z * k.z);
        kAbs[index] = kk;
        kMax = std::max(kMax, kk);
      }
    }
  }

  // The weights of the hard-sphere functional, laid onto the grid in real space and transformed from there.
  //
  // The obvious alternative is their analytic radial transforms on the reciprocal lattice, which is what a grid cannot
  // take: the transform of the shell, 4 pi R sin(kR) / k, has hardly decayed by the grid's last wavevector, and cutting
  // it off there rings. On a zeolite that left n2 negative over a quarter of the cell and |n2v| at six times n2. The
  // functional's repulsion is carried by n2, so where n2 had gone negative the fluid felt none, whatever the packing,
  // and the density went there. No solution was to be had past a few molecules per cage.
  //
  // A kernel of non-negative weights cannot do that. With rho >= 0, each of n3 and n2 is a sum of non-negative terms,
  // and because the vector weight is the same shell samples times their unit vectors, |n2v| <= n2 by the triangle
  // inequality, exactly, on any grid. Samples come in antipodal pairs, so the scalar kernels transform to a real
  // spectrum and the vector kernel to an imaginary one, and each is normalised at k = 0, so the uniform fluid and with
  // it the whole of the bulk are untouched. Sharing a sample among its eight grid points smooths the kernel by about a
  // voxel, which is second order in the spacing.
  {
    const double spacing = std::min({cell.lengthA / static_cast<double>(nx), cell.lengthB / static_cast<double>(ny),
                                     cell.lengthC / static_cast<double>(nz)});
    const double sampleSpacing = 0.25 * spacing;
    auto forEachOnShell = [&](double radius, auto &&visit)
    {
      const std::size_t pairs = std::max<std::size_t>(
          64, static_cast<std::size_t>(2.0 * std::numbers::pi * radius * radius / (sampleSpacing * sampleSpacing)));
      const double goldenAngle = std::numbers::pi * (3.0 - std::sqrt(5.0));
      for (std::size_t q = 0; q < pairs; ++q)
      {
        const double z = (static_cast<double>(q) + 0.5) / static_cast<double>(pairs);
        const double r = std::sqrt(std::max(0.0, 1.0 - z * z));
        const double angle = goldenAngle * static_cast<double>(q);
        const double3 direction(r * std::cos(angle), r * std::sin(angle), z);
        const double share = 0.5 / static_cast<double>(pairs);
        visit(direction, share);
        visit(-1.0 * direction, share);
      }
    };

    RealArray shell = makeReal();
    RealArray ball = makeReal();
    RealArray shellX = makeReal();
    RealArray shellY = makeReal();
    RealArray shellZ = makeReal();
    for (std::size_t i = 0; i < N; ++i)
    {
      shell[i] = 0.0;
      ball[i] = 0.0;
      shellX[i] = 0.0;
      shellY[i] = 0.0;
      shellZ[i] = 0.0;
    }
    auto depositAt = [&](RealArray &grid, const double3 &position, double weight)
    {
      const double3 f = cell.inverseCell * position;
      nldftCicDeposit(grid.get(), nx, ny, nz, f.x, f.y, f.z, weight);
    };
    const double shellArea = 4.0 * std::numbers::pi * R * R;
    forEachOnShell(R, [&](const double3 &direction, double share)
    {
      const double3 position = R * direction;
      depositAt(shell, position, shellArea * share);
      depositAt(shellX, position, shellArea * share * direction.x);
      depositAt(shellY, position, shellArea * share * direction.y);
      depositAt(shellZ, position, shellArea * share * direction.z);
    });
    const std::size_t radii = std::max<std::size_t>(16, static_cast<std::size_t>(2.0 * R / spacing));
    double ballVolume = 0.0;
    for (std::size_t q = 0; q < radii; ++q)
    {
      const double r = R * (static_cast<double>(q) + 0.5) / static_cast<double>(radii);
      const double volume = 4.0 * std::numbers::pi * r * r * R / static_cast<double>(radii);
      ballVolume += volume;
      forEachOnShell(r, [&](const double3 &direction, double share)
                     { depositAt(ball, r * direction, volume * share); });
    }
    const double ballScale = (4.0 * std::numbers::pi * R * R * R / 3.0) / ballVolume;
    for (std::size_t i = 0; i < N; ++i) ball[i] *= ballScale;

    // A kernel transforms as a sum over e^{-ikr}, so an odd one is purely imaginary; vectorWeightTo takes c in -i c.
    auto transform = [&](RealArray &grid, std::vector<double> &spectrum, int part, double sign)
    {
      fftw_execute_dft_r2c(forward, grid.get(), rhoHat.get());
      for (std::size_t index = 0; index < nC; ++index) spectrum[index] = sign * rhoHat[index][part];
    };
    transform(shell, w2Hat, 0, 1.0);
    transform(ball, w3Hat, 0, 1.0);
    transform(shellX, vHatX, 1, -1.0);
    transform(shellY, vHatY, 1, -1.0);
    transform(shellZ, vHatZ, 1, -1.0);
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
    // a 1x1x1 grid has only k = 0 (kMax = 0): the table then spans an arbitrary finite range and every lattice
    // point reads its first entry, instead of 0 / 0
    const double kTable = kMax > 0.0 ? kMax : 1.0;
    std::vector<double> uTable(tableSize);
    uTable[0] = (nAtt > 0.0) ? model.a / (nAtt * nAtt) : model.a;
    for (std::size_t t = 1; t < tableSize; ++t)
    {
      double k = kTable * static_cast<double>(t) / static_cast<double>(tableSize - 1);
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
      double position = kAbs[index] / kTable * static_cast<double>(tableSize - 1);
      std::size_t bin = std::min(static_cast<std::size_t>(position), tableSize - 2);
      double fraction = position - static_cast<double>(bin);
      uHat[index] = uTable[bin] * (1.0 - fraction) + uTable[bin + 1] * fraction;
    }
  }

  // beta U for the Euler-Lagrange map. Overlap (βU > 690) is dropped. Spherical ρ clips wells at
  // wellFloorInKT: a Helmholtz PMF can be a Coulomb hole a spherical density cannot occupy. Molecular
  // ρ(r, ω) keeps the TraPPE well (FMT packing is the bound); clipping it at 12 kT made every μ see
  // the same e^12 contact layer and stalled FER at ~4 / cell against GCMC ~10. Henry is always unclipped.
  //
  // A molecular density is a density of states, one for each orientation at each voxel, and only those a molecule can
  // be found in are kept. One that lies more than largestLiveBetaU above zero is not among them: the field of the
  // other molecules lowers it by ten kT at the most and the reservoir is no denser than the saturated gas, so every
  // such state in a cell together holds less than e^-40 of a molecule. On a framework that is most of the states there
  // are, the walls being most of the cell, and what is left is what the descent has to carry.
  constexpr double largestLiveBetaU = 50.0;
  const double wellFloor = options.wellFloorInKT > 0.0 ? options.wellFloorInKT : nldftWellFloorInKT;
  std::vector<double> betaU(molecular ? 0 : N);
  std::vector<std::size_t> liveStart(N + 1, 0);
  std::vector<double> liveLogWeight;
  std::vector<std::uint16_t> liveOrient;
  // the live orientations are kept as 16-bit indices: more orientations than that wrapped silently
  if (molecular && M > static_cast<std::size_t>(std::numeric_limits<std::uint16_t>::max()) + 1)
  {
    throw std::invalid_argument(std::format("NLDFT: {} orientations is more than the {} the orientation index holds\n",
                                            M,
                                            static_cast<std::size_t>(std::numeric_limits<std::uint16_t>::max()) + 1));
  }
  if (molecular)
  {
    for (std::size_t i = 0; i < N; ++i)
    {
      liveStart[i] = liveLogWeight.size();
      for (std::size_t o = 0; o < M; ++o)
      {
        const double bu = beta * static_cast<double>(orientView[i * M + o]);
        if (!(bu <= largestLiveBetaU)) continue;
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
    }
  }
  const std::size_t K = liveLogWeight.size();
  const double logM = std::log(static_cast<double>(M));

  const bool needDeposit = dumbbell || siteAtt;
  const std::vector<double> &geomOffsets =
      siteAtt ? options.attractiveSiteOffsets : options.dumbbellSiteOffsets;
  const std::size_t nSites = geomOffsets.size();
  const double siteScale =
      (dumbbell && nSites > 0) ? model.sites / static_cast<double>(nSites) : 1.0;
  const double depositScale = dumbbell ? siteScale : 1.0;

  // The sites are deposited scaled, so that the packing of the site density is the fused molecule's and the
  // hard-sphere part matches the bulk's. The attraction counts every site in full: in the uniform fluid the
  // sum over a molecule's sites of the convolved field has to come to a rho, the mean-field term of the bulk
  // chemical potential, which is n_att^2 a_site rho and so wants the n_att rho of sites unscaled. Convolving
  // the scaled deposit as it stood left the pore fluid short of its attraction by the scale, about a quarter
  // for the nitrogen dumbbell, and the uniform liquid was no longer stationary at the bulk's own chemical
  // potential. The scaled deposit is reused, with the scale taken back out of the pair spectrum.
  std::vector<double> uHatSites = uHat;
  if (siteAtt && dumbbell && depositScale > 0.0)
  {
    for (double &u : uHatSites) u /= depositScale;
  }
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

  const unsigned workers = nldftWorkerCount();
  std::vector<std::vector<double>> siteLocal;
  if (needDeposit)
  {
    siteLocal.resize(workers);
    for (auto &local : siteLocal) local.assign(N, 0.0);
  }

  // Sums and loops over the voxels, or over the states, shared among the workers.
  auto workerOf = [&](std::size_t begin, std::size_t count) -> std::size_t
  {
    const std::size_t chunk = (count + workers - 1) / std::max(workers, 1u);
    return std::min(begin / std::max(chunk, std::size_t{1}), static_cast<std::size_t>(workers - 1));
  };
  auto sumOf = [&](std::size_t count, auto &&term) -> double
  {
    std::vector<double> sums(workers, 0.0);
    nldftParallelFor(count, workers, [&](std::size_t begin, std::size_t end)
    {
      double sum = 0.0;
      for (std::size_t i = begin; i < end; ++i) sum += term(i);
      sums[workerOf(begin, count)] = sum;
    });
    return std::accumulate(sums.begin(), sums.end(), 0.0);
  };
  auto eachOf = [&](std::size_t count, auto &&apply)
  {
    nldftParallelFor(count, workers, [&](std::size_t begin, std::size_t end)
    {
      for (std::size_t i = begin; i < end; ++i) apply(i);
    });
  };
  auto sumOver = [&](auto &&term) -> double { return sumOf(N, term); };
  auto eachVoxel = [&](auto &&apply) { eachOf(N, apply); };

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
  // n2v_j = IFFT( -i c_j rho ): (a + b i)(-i c) = b c - a c i.
  auto vectorWeightTo = [&](const std::vector<double> &spectrum, double *out)
  {
    for (std::size_t index = 0; index < nC; ++index)
    {
      double c = spectrum[index] * inverseN;
      scratch[index][0] = rhoHat[index][1] * c;
      scratch[index][1] = -rhoHat[index][0] * c;
    }
    fftw_execute_dft_c2r(inverse, scratch.get(), out);
  };

  // The molecular density itself is the molecules in each live state. The density at a voxel, rho, is their sum there,
  // and the density of sites, rhoSite, is each state laid down at its sites. Placing the states makes the three one,
  // and nothing else changes any of them, so the orientations a field is worked out from are those of the density it
  // is the field of. Sharing the molecules at a voxel among the orientations by the fields of the step before, as an
  // iteration in rho alone has to, left them a function of the path until it had converged.
  std::vector<double> rhoState(K, 0.0);
  RealArray siteField = needDeposit ? makeReal() : RealArray(nullptr, realFree);
  constexpr double smallestDeposit = 1.0e-40;  // molecules / Å³ in one state, below which it is laid down nowhere
  auto fractionOf = [&](std::size_t i) -> double3
  {
    return double3(static_cast<double>(i % nx) / static_cast<double>(nx),
                   static_cast<double>((i / nx) % ny) / static_cast<double>(ny),
                   static_cast<double>(i / (nx * ny)) / static_cast<double>(nz));
  };
  auto placeStates = [&]()
  {
    eachVoxel([&](std::size_t i)
    {
      double sum = 0.0;
      for (std::size_t k = liveStart[i]; k < liveStart[i + 1]; ++k) sum += rhoState[k];
      rho[i] = sum;
    });
    if (!needDeposit) return;

    for (auto &local : siteLocal) std::fill(local.begin(), local.end(), 0.0);
    nldftParallelFor(N, workers, [&](std::size_t begin, std::size_t end)
    {
      double *local = siteLocal[workerOf(begin, N)].data();
      for (std::size_t i = begin; i < end; ++i)
      {
        if (liveStart[i] == liveStart[i + 1]) continue;
        const double3 f = fractionOf(i);
        for (std::size_t k = liveStart[i]; k < liveStart[i + 1]; ++k)
        {
          const double held = rhoState[k];
          if (!(held > smallestDeposit)) continue;
          const std::size_t o = liveOrient[k];
          for (std::size_t s = 0; s < nSites; ++s)
          {
            const double3 df = siteDelta[o * nSites + s];
            nldftCicAdd(local, nx, ny, nz, f.x + df.x, f.y + df.y, f.z + df.z, depositScale * held);
          }
        }
      }
    });
    eachVoxel([&](std::size_t i)
    {
      double sum = 0.0;
      for (unsigned t = 0; t < workers; ++t) sum += siteLocal[t][i];
      rhoSite[i] = sum;
    });
  };

  // The third term of the free energy density is phi3(n3) n2^3 (1 - xi^2)^3, with xi = |n2v| / n2, the form of
  // Rosenfeld, Schmidt, Lowen and Tarazona (Phys. Rev. E 55, 4245 (1997)), rather than phi3(n3) (n2^3 - 3 n2 n2v^2).
  // The two agree where xi is small, and in the uniform fluid exactly. They part where a density peak is narrow, and
  // that is what a zeolite makes: phi3 diverges as n3 goes to one, and past xi^2 = 1/3 the other form goes with it to
  // minus infinity, so the more tightly the fluid packs the lower the free energy, with no bottom to it. Here every
  // term of the density is non-negative, the functional is bounded below, and a pocket that fits one molecule holds
  // one.
  //
  // The weights keep n2^2 - n2v^2 from going negative but for rounding, which is the one thing the guard here is for.
  struct ThirdTerm
  {
    double value{0.0};
    double byN2{0.0};
    double byVector{0.0};  // the derivative by n2v_j is this times n2v_j
  };
  auto thirdTermOf = [](double n2, double vSquared) -> ThirdTerm
  {
    if (!(n2 > 1.0e-100)) return {};
    const double gap = std::max(n2 * n2 - vSquared, 0.0);
    const double inverse = 1.0 / n2;
    const double gapOverN2Squared = gap * inverse * inverse;
    const double gapSquaredOverN2Cubed = gap * gapOverN2Squared * inverse;
    return {.value = gap * gapSquaredOverN2Cubed,
            .byN2 = 3.0 * gapOverN2Squared * gapOverN2Squared * (n2 * n2 + vSquared),
            .byVector = -6.0 * gapSquaredOverN2Cubed};
  };

  // The fields of the density in hand: the attractive one, and the hard-sphere one, which is the derivative of the
  // White Bear free energy by the density it is a functional of. For a dumbbell that is the density of sites.
  auto fieldsOf = [&]()
  {
    fftw_execute_dft_r2c(forward, rho.get(), rhoHat.get());
    if (!siteAtt) convolveTo(uHat, phiAtt.get());
    if (dumbbell) fftw_execute_dft_r2c(forward, rhoSite.get(), rhoHat.get());
    if (siteAtt && dumbbell) convolveTo(uHatSites, phiAtt.get());
    convolveTo(w3Hat, n3B.get());
    convolveTo(w2Hat, n2A.get());
    vectorWeightTo(vHatX, vx.get());
    vectorWeightTo(vHatY, vy.get());
    vectorWeightTo(vHatZ, vz.get());

    nldftParallelFor(N, workers, [&](std::size_t begin, std::size_t end)
    {
      for (std::size_t i = begin; i < end; ++i)
      {
        double n3 = std::clamp(n3B[i], 0.0, largestPacking);
        double n2 = std::max(n2A[i], 0.0);
        double wx = vx[i];
        double wy = vy[i];
        double wz = vz[i];
        double vSquared = wx * wx + wy * wy + wz * wz;
        double o = 1.0 / (1.0 - n3);
        double phi3 = 0.0;
        double phi3Prime = 0.0;
        phi3Of(n3, phi3, phi3Prime);

        double n0 = s0 * n2;
        double n1 = s1 * n2;
        double dPhi0 = -std::log(1.0 - n3);
        double dPhi1 = n2 * o;
        const ThirdTerm third = thirdTermOf(n2, vSquared);
        double dPhi2 = n1 * o + third.byN2 * phi3;
        double dPhi3 = n0 * o + (n1 * n2 - s1 * vSquared) * o * o + third.value * phi3Prime;

        n2A[i] = s0 * dPhi0 + s1 * dPhi1 + dPhi2;
        n3B[i] = dPhi3;
        double vectorFactor = -2.0 * s1 * o + third.byVector * phi3;
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
    const std::vector<double> *vectorSpectra[3] = {&vHatX, &vHatY, &vHatZ};
    double *vectorFields[3] = {vx.get(), vy.get(), vz.get()};
    for (int component = 0; component < 3; ++component)
    {
      fftw_execute_dft_r2c(forward, vectorFields[component], scratch.get());
      const std::vector<double> &spectrum = *vectorSpectra[component];
      for (std::size_t index = 0; index < nC; ++index)
      {
        double c = spectrum[index];
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
    if (siteAtt && !dumbbell)
    {
      fftw_execute_dft_r2c(forward, rhoSite.get(), rhoHat.get());
      convolveTo(uHat, phiAtt.get());
    }
  };

  constexpr double logDensityFloor = -690.0;

  // What the Euler-Lagrange equation makes of a spherical density, and how far from it the density is.
  auto eulerLagrange = [&](double muOverT, std::vector<double> &updated) -> double
  {
    fieldsOf();
    constexpr double tiny = 1.0e-300;
    return std::sqrt(sumOver([&](std::size_t i)
    {
      const double exponent = muOverT - beta * phiAtt[i] - muHS[i] - betaU[i];
      const double value = (exponent < logDensityFloor) ? tiny : std::exp(std::min(exponent, 80.0));
      const double difference = value - rho[i];
      updated[i] = value;
      return difference * difference;
    }) / static_cast<double>(N));
  };

  // The same for a molecular density, state by state: the logarithm of what the equation makes of each, the gradient
  // of the grand potential by the logarithm of what it holds, and how far the two are apart over all of them, against
  // the density itself. A state feels the fields at its sites where the functional is one of sites and at its centre
  // where it is one of centres.
  auto stateGradientAt = [&](double muOverT, const std::vector<double> &logValues, std::vector<double> &out,
                             std::vector<double> &logUpdate) -> double
  {
    fieldsOf();
    if (needDeposit)
    {
      eachVoxel([&](std::size_t i)
      {
        siteField[i] = (dumbbell ? siteScale * muHS[i] : 0.0) + (siteAtt ? beta * phiAtt[i] : 0.0);
      });
    }
    std::vector<double> apart(workers, 0.0);
    std::vector<double> held(workers, 0.0);
    nldftParallelFor(N, workers, [&](std::size_t begin, std::size_t end)
    {
      double sumApart = 0.0;
      double sumHeld = 0.0;
      for (std::size_t i = begin; i < end; ++i)
      {
        if (liveStart[i] == liveStart[i + 1]) continue;
        const double3 f = fractionOf(i);
        const double atTheCentre = (dumbbell ? 0.0 : muHS[i]) + (siteAtt ? 0.0 : beta * phiAtt[i]);
        for (std::size_t k = liveStart[i]; k < liveStart[i + 1]; ++k)
        {
          const double field =
              atTheCentre + (needDeposit ? sampleSites(siteField.get(), liveOrient[k], f.x, f.y, f.z) : 0.0);
          const double logValue = std::clamp(muOverT - logM + liveLogWeight[k] - field, logDensityFloor, 80.0);
          const double difference = std::exp(logValue) - rhoState[k];
          sumApart += difference * difference;
          sumHeld += rhoState[k] * rhoState[k];
          logUpdate[k] = logValue;
          out[k] = rhoState[k] * dV * (logValues[k] - logValue);
        }
      }
      apart[workerOf(begin, N)] = sumApart;
      held[workerOf(begin, N)] = sumHeld;
    });
    return std::sqrt(std::accumulate(apart.begin(), apart.end(), 0.0) /
                     std::max(std::accumulate(held.begin(), held.end(), 0.0), 1.0e-300));
  };

  auto loadingOf = [&]()
  {
    double sum = 0.0;
    for (std::size_t i = 0; i < N; ++i) sum += rho[i];
    return sum * dV;
  };

  // The grand potential of the profile in hand, in kT per cell, the functional the updates above are the stationary
  // point of, term by term:
  //
  //   Omega = sum_w int rho_w (ln(M rho_w) - 1 + beta U_w - mu) + int Phi(n[rho_site]) + 1/2 beta int rho_att phi
  //
  // rho_w being the molecules in orientation w, Phi the White Bear free energy density of the same weighted densities
  // the update differentiates, under the same clamps, and phi the attractive field. Differentiating this with respect
  // to rho_w gives back that update, which is what lets a descent be held to it. The sites of the attraction count in
  // full where those of the hard-sphere deposit are scaled, so the scale comes back out here as the update takes it out
  // of the pair spectrum.
  RealArray omegaPhi = makeReal();
  auto grandPotentialOf = [&](double muOverT) -> double
  {
    fftw_execute_dft_r2c(forward, dumbbell ? rhoSite.get() : rho.get(), rhoHat.get());
    convolveTo(w3Hat, n3B.get());
    convolveTo(w2Hat, n2A.get());
    vectorWeightTo(vHatX, vx.get());
    vectorWeightTo(vHatY, vy.get());
    vectorWeightTo(vHatZ, vz.get());

    const double hardSphere = sumOver([&](std::size_t i)
    {
      double n3 = std::clamp(n3B[i], 0.0, largestPacking);
      double n2 = std::max(n2A[i], 0.0);
      double vSquared = vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i];
      double phi3 = 0.0;
      double phi3Prime = 0.0;
      phi3Of(n3, phi3, phi3Prime);
      double n0 = s0 * n2;
      double n1 = s1 * n2;
      return -n0 * std::log(1.0 - n3) + (n1 * n2 - s1 * vSquared) / (1.0 - n3) +
             thirdTermOf(n2, vSquared).value * phi3;
    });

    double attraction = 0.0;
    if (siteAtt)
    {
      if (!dumbbell) fftw_execute_dft_r2c(forward, rhoSite.get(), rhoHat.get());
      convolveTo(dumbbell ? uHatSites : uHat, omegaPhi.get());
      attraction = sumOver([&](std::size_t i) { return rhoSite[i] * omegaPhi[i]; }) / depositScale;
    }
    else
    {
      if (dumbbell) fftw_execute_dft_r2c(forward, rho.get(), rhoHat.get());
      convolveTo(uHat, omegaPhi.get());
      attraction = sumOver([&](std::size_t i) { return rho[i] * omegaPhi[i]; });
    }

    const double idealAndField =
        molecular ? sumOf(K, [&](std::size_t k)
                    {
                      const double held = rhoState[k];
                      return held > 0.0 ? held * (std::log(held) + logM - 1.0 - liveLogWeight[k] - muOverT) : 0.0;
                    })
                  : sumOver([&](std::size_t i)
                    {
                      return rho[i] > 0.0 ? rho[i] * (std::log(rho[i]) - 1.0 + betaU[i] - muOverT) : 0.0;
                    });
    return dV * (idealAndField + hardSphere + 0.5 * beta * attraction);
  };

  // How far the profile last solved for was from the equation, against the density it is measured by.
  double relativeResidual = 0.0;

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

  // What the descent goes down in: the logarithm of the density at each voxel, or for a molecular density of what each
  // live state holds.
  const std::size_t unknowns = molecular ? K : N;
  double *const density = molecular ? rhoState.data() : rho.get();

  // The last converged profile, which is all a warm start starts from. The fields are worked out from the density each
  // time they are wanted, so there is nothing else to carry from one pressure to the next, and nothing a point that
  // failed can leave behind for the one after it.
  std::vector<double> updated(molecular ? 0 : N);
  std::vector<double> lastConverged(unknowns);
  auto keepConverged = [&]() { eachOf(unknowns, [&](std::size_t v) { lastConverged[v] = density[v]; }); };
  auto restoreConverged = [&]()
  {
    eachOf(unknowns, [&](std::size_t v) { density[v] = lastConverged[v]; });
    if (molecular) placeStates();
  };

  bool haveProfile = false;
  bool hasPrev = false;
  double xPrev = 0.0;
  double nPrev = 0.0;
  std::size_t extraPoints = 0;
  const std::size_t maxExtra = static_cast<std::size_t>(std::max(options.pressurePoints, 1));
  std::vector<NLDFTBranchPoint> adsorption;
  adsorption.reserve(targets.size() + maxExtra);

  // The Henry profile of the gas, or the liquid's with the wells taken out, which is what a sweep starts from. A deep
  // well is held to e^8 of the gas: the descent backs a start off until it does not overlap itself, and there is no
  // call to make it back off far.
  auto seedDensity = [&](double rhoBulk)
  {
    const double base = options.seedLiquid ? result.bulk.liquidDensity : rhoBulk;
    const double deepest = options.seedLiquid ? 0.0 : 8.0;
    if (molecular)
    {
      eachOf(K, [&](std::size_t k)
      {
        rhoState[k] = base / static_cast<double>(M) * std::exp(std::min(liveLogWeight[k], deepest));
      });
      placeStates();
      return;
    }
    eachVoxel([&](std::size_t i) { rho[i] = base * std::exp(std::min(-betaU[i], deepest)); });
  };

  // The solve: going down the grand potential itself with L-BFGS over the logarithm of the density. Omega is a true
  // functional of the density, with the gradient rho dV (ln rho - ln rho_update), so every step can be held to it: a
  // step that does not lower Omega is cut, and that keeps the profile on the right side of the barrier at n3 = 1, where
  // the functional gets steeper than any fixed-point iteration can follow. The first direction, scaled by the ideal
  // part of the Hessian, is simply a step to the update in ln rho. Log variables let a pore fill in a few steps, which
  // damping in rho takes many to do.
  //
  // A molecular density goes down the same way, in what each of its states holds. An iteration in the density at a
  // voxel has to share that density among the orientations by fields that are not yet those of the density it will
  // make, which is no functional of anything until it has converged, and past a few molecules a cell it stopped moving
  // well short of that. With the states as the unknowns there is nothing to share out.
  constexpr std::size_t deepestHistory = 10;
  constexpr double historyBytes = 4.0e9;  // ten pairs over the states of a large cell would be most of a machine
  const std::size_t historyDepth = std::clamp(
      static_cast<std::size_t>(historyBytes / (16.0 * static_cast<double>(std::max(unknowns, 1uz)))), 3uz,
      deepestHistory);
  constexpr double largestRise = 2.0;
  // molecules / Å³, a part in a million of the liquid, and for a state its share of that
  const double thatCounts = 1.0e-8 / static_cast<double>(M);
  const double logThatCounts = std::log(thatCounts);
  std::vector<std::vector<double>> stepHistory(historyDepth, std::vector<double>(unknowns));
  std::vector<std::vector<double>> gradientHistory(historyDepth, std::vector<double>(unknowns));
  std::vector<double> logRho(unknowns), gradient(unknowns), direction(unknowns), logTrial(unknowns),
      gradientTrial(unknowns), logUpdate(unknowns);
  std::array<double, deepestHistory> curvature{};
  std::array<double, deepestHistory> projection{};

  // Putting a set of logarithms in as the density, with all that is kept in step with it.
  auto place = [&](const std::vector<double> &logValues)
  {
    eachOf(unknowns, [&](std::size_t v) { density[v] = std::exp(logValues[v]); });
    if (molecular) placeStates();
  };

  auto packingIsPossible = [&]() -> bool
  {
    fftw_execute_dft_r2c(forward, dumbbell ? rhoSite.get() : rho.get(), rhoHat.get());
    convolveTo(w3Hat, omegaPhi.get());
    return *std::max_element(omegaPhi.get(), omegaPhi.get() + static_cast<std::ptrdiff_t>(N)) < largestPacking;
  };

  // The gradient of the grand potential by the logarithms, with the logarithm of the update of each alongside it, and
  // what comes back is the residual: how far the density is from its update, against the density.
  auto gradientAt = [&](double muOverT, const std::vector<double> &logValues, std::vector<double> &out) -> double
  {
    if (molecular) return stateGradientAt(muOverT, logValues, out, logUpdate);
    const double residual = eulerLagrange(muOverT, updated);
    eachVoxel([&](std::size_t i)
    {
      logUpdate[i] = std::max(std::log(std::max(updated[i], 1.0e-300)), logDensityFloor);
      out[i] = rho[i] * dV * (logValues[i] - logUpdate[i]);
    });
    const double scale = std::sqrt(sumOver([&](std::size_t i) { return rho[i] * rho[i]; }) / static_cast<double>(N));
    return residual / std::max(scale, 1.0e-300);
  };

  constexpr std::size_t steadyWindow = 100;
  const double steadyTolerance = options.steadyTolerance > 0.0 ? options.steadyTolerance : 0.1 * options.tolerance;
  auto descend = [&](double muOverT) -> bool
  {
    eachOf(unknowns, [&](std::size_t v)
    {
      logRho[v] = std::max(std::log(std::max(density[v], 1.0e-300)), logDensityFloor);
    });
    place(logRho);

    // A starting guess can overlap itself, which a seed at a deep well does: back the whole profile off until it does
    // not.
    while (!packingIsPossible())
    {
      eachOf(unknowns, [&](std::size_t v) { logRho[v] -= 0.1; });
      place(logRho);
    }

    // What holds next to nothing feels the rest and is not felt by it, so its minimum, the rest as it is, is the update
    // itself, and that is where it is sent. It cannot be left to the gradient: that is its density times how far it is
    // from the update, and a density that has gone to nothing has none, however far off it is. A pore that fills
    // pushes the density between its molecules down by tens of orders of magnitude on the way, and on the way back
    // down in pressure those voxels stayed where they were, the residual with them, for as long as the iteration was
    // let run. The update is that of the profile in hand, since the last gradient taken is always the one kept.
    auto toTheUpdateWhereNothingIsHeld = [&]()
    {
      eachOf(unknowns, [&](std::size_t v)
      {
        if (logRho[v] > logThatCounts) return;
        direction[v] = logUpdate[v] - logRho[v];
      });
    };

    std::size_t stored = 0;
    std::size_t head = 0;
    double hessianScale = 1.0;
    double residual = gradientAt(muOverT, logRho, gradient);
    double omega = grandPotentialOf(muOverT);
    double omegaBefore = omega;
    double loadingBefore = loadingOf();
    for (std::size_t iteration = 0; iteration < options.maxIterations; ++iteration)
    {
      relativeResidual = residual;
      if (relativeResidual < options.tolerance) return true;
      // A molecule packed close to solid has a mode stiff enough that the residual, the exponential of what is left,
      // is still coming down hundreds of steps after Omega and the loading have stopped showing it, and there is one
      // such mode per molecule, more than the history holds. So this does not wait for the residual. What it ends on is
      // short of where the residual would have ended by more than the change it watches: that is over a hundred steps,
      // and the tail can be several hundred long, which on a cage packed to solid left the loading a part and a half in
      // ten thousand out. It is for a profile that is still moving, slowly. One that has stopped because it cannot move
      // is not converged, and what did that is dealt with where the direction is set.
      if (iteration > 0 && iteration % steadyWindow == 0)
      {
        const double loadingNow = loadingOf();
        const double allowed = steadyTolerance * std::max(loadingNow, 1.0);
        if (std::fabs(omega - omegaBefore) < allowed && std::fabs(loadingNow - loadingBefore) < allowed) return true;
        omegaBefore = omega;
        loadingBefore = loadingNow;
      }

      eachOf(unknowns, [&](std::size_t v) { direction[v] = gradient[v]; });
      for (std::size_t k = 0; k < stored; ++k)
      {
        const std::size_t j = (head + historyDepth - 1 - k) % historyDepth;
        const std::vector<double> &step = stepHistory[j];
        const std::vector<double> &change = gradientHistory[j];
        projection[j] = curvature[j] * sumOf(unknowns, [&](std::size_t v) { return step[v] * direction[v]; });
        eachOf(unknowns, [&](std::size_t v) { direction[v] -= projection[j] * change[v]; });
      }
      const double scale = stored > 0 ? hessianScale : 1.0;
      eachOf(unknowns, [&](std::size_t v) { direction[v] *= scale / (std::max(density[v], 1.0e-300) * dV); });
      for (std::size_t k = stored; k-- > 0;)
      {
        const std::size_t j = (head + historyDepth - 1 - k) % historyDepth;
        const std::vector<double> &step = stepHistory[j];
        const std::vector<double> &change = gradientHistory[j];
        const double back = curvature[j] * sumOf(unknowns, [&](std::size_t v) { return change[v] * direction[v]; });
        eachOf(unknowns, [&](std::size_t v) { direction[v] += (projection[j] - back) * step[v]; });
      }
      eachOf(unknowns, [&](std::size_t v) { direction[v] = -direction[v]; });
      toTheUpdateWhereNothingIsHeld();

      double slope = sumOf(unknowns, [&](std::size_t v) { return gradient[v] * direction[v]; });
      if (!(slope < 0.0))
      {
        stored = 0;
        eachOf(unknowns, [&](std::size_t v) { direction[v] = logUpdate[v] - logRho[v]; });
        slope = sumOf(unknowns, [&](std::size_t v) { return gradient[v] * direction[v]; });
      }

      // The step is measured by what holds anything. One state of no account, with a direction as large as its density
      // is small, could otherwise have every other taking steps too small to see.
      double largestChange = 0.0;
      for (std::size_t v = 0; v < unknowns; ++v)
      {
        if (logRho[v] > logThatCounts) largestChange = std::max(largestChange, std::fabs(direction[v]));
      }

      double stepLength = std::min(1.0, largestRise / std::max(largestChange, 1.0e-12));
      bool accepted = false;
      double omegaTrial = omega;
      double residualTrial = residual;
      for (int trial = 0; trial < 40; ++trial)
      {
        eachOf(unknowns, [&](std::size_t v)
        {
          // What holds next to nothing may come up as far as where it starts to count in a single step, being of no
          // account until then, and no further than any other from there on. It must not overflow: the release build
          // has no NaNs to catch, and a comparison of Omega against one says nothing.
          const double rise = std::max(largestRise, logThatCounts - logRho[v]);
          const double change = std::clamp(stepLength * direction[v], -largestRise, rise);
          logTrial[v] = std::max(logRho[v] + change, logDensityFloor);
        });
        place(logTrial);
        if (packingIsPossible())
        {
          omegaTrial = grandPotentialOf(muOverT);
          if (omegaTrial <= omega + 1.0e-4 * stepLength * slope + 1.0e-11 * std::fabs(omega))
          {
            residualTrial = gradientAt(muOverT, logTrial, gradientTrial);
            accepted = true;
            break;
          }
        }
        stepLength *= 0.5;
      }

      if (!accepted)
      {
        place(logRho);
        if (stored == 0) return false;
        stored = 0;
        residual = gradientAt(muOverT, logRho, gradient);
        continue;
      }

      // The history is of what counts. The rest is sent to the update whatever the history says, and the ideal term in
      // ln rho curves the wrong way that far below its minimum, which is no curvature to learn from.
      std::vector<double> &step = stepHistory[head];
      std::vector<double> &change = gradientHistory[head];
      eachOf(unknowns, [&](std::size_t v)
      {
        const bool counts = logRho[v] > logThatCounts;
        step[v] = counts ? logTrial[v] - logRho[v] : 0.0;
        change[v] = counts ? gradientTrial[v] - gradient[v] : 0.0;
      });
      const double stepAlongChange = sumOf(unknowns, [&](std::size_t v) { return step[v] * change[v]; });
      if (stepAlongChange > 1.0e-14)
      {
        // The ideal curvature alone is too soft by orders of magnitude at high density, so the first step of every
        // iteration would otherwise be cut back many times over, and every cut costs a grand potential.
        const double changeSquared = sumOf(unknowns, [&](std::size_t v)
        {
          return change[v] * change[v] / (std::max(density[v], 1.0e-300) * dV);
        });
        hessianScale = stepAlongChange / std::max(changeSquared, 1.0e-300);
        curvature[head] = 1.0 / stepAlongChange;
        head = (head + 1) % historyDepth;
        stored = std::min(stored + 1, historyDepth);
      }
      logRho.swap(logTrial);
      gradient.swap(gradientTrial);
      omega = omegaTrial;
      residual = residualTrial;
    }
    return false;
  };

  while (!targets.empty())
  {
    double x = targets.front();
    targets.pop_front();
    double pressure = x * result.bulk.saturationPressure / pressureKA3ToPascal;
    double rhoBulk = gasDensityAtPressure(model, pressure, rhoGasCap);
    double muOverT = model.muOf(rhoBulk) / T;

    if (!haveProfile)
    {
      seedDensity(rhoBulk);
    }
    else
    {
      restoreConverged();
    }

    const bool converged = descend(muOverT);

    double loading = loadingOf();
    if (!converged)
    {
      ++result.unconverged;
      if (haveProfile)
      {
        restoreConverged();
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
        if (haveProfile) restoreConverged();
        continue;
      }
    }

    // Adsorption sweep: a fill that the next pressure then lets go of is kept, and that pressure is given no state of
    // its own. FER reached n ≈ 10 (the GCMC plateau) and the next x fell back to the four-molecule single-file branch.
    bool heldAdsorption = false;
    if (hasPrev && haveProfile && nPrev - loading > std::max(0.2 * nPrev, 0.05))
    {
      restoreConverged();
      loading = nPrev;
      heldAdsorption = true;
    }

    NLDFTBranchPoint point{.relativePressure = x, .moleculesPerCell = loading, .chemicalPotential = muOverT};
    if (!heldAdsorption && converged)
    {
      point.grandPotential = grandPotentialOf(muOverT);
      point.residual = relativeResidual;
      point.stationary = true;
      keepConverged();
      haveProfile = true;
    }

    adsorption.push_back(point);
    xPrev = x;
    nPrev = loading;
    hasPrev = true;
  }

  result.adsorptionIsotherm = adsorption;
  result.desorptionIsotherm = adsorption;

  // Desorption branch, walked down from the filled profile. The sweep up of mean-field DFT delays a fill until the
  // pore's spinodal (FER: n = 4 until x = 0.27, then 11), and the sweep down delays the emptying just as far the other
  // way, so neither is the state GCMC samples. Between them they hold the two states the comparison below is among.
  if (haveProfile && adsorption.size() > 4 &&
      adsorption.back().moleculesPerCell > 1.5 * adsorption.front().moleculesPerCell + 1.0)
  {
    std::vector<double> xDown;
    xDown.reserve(adsorption.size());
    for (auto it = adsorption.rbegin(); it != adsorption.rend(); ++it) xDown.push_back(it->relativePressure);
    std::vector<NLDFTBranchPoint> desorption;
    desorption.reserve(xDown.size());
    desorption.push_back(adsorption.back());
    for (std::size_t k = 1; k < xDown.size(); ++k)
    {
      const double x = xDown[k];
      const double pressure = x * result.bulk.saturationPressure / pressureKA3ToPascal;
      const double rhoBulk = gasDensityAtPressure(model, pressure, rhoGasCap);
      const double muOverT = model.muOf(rhoBulk) / T;
      restoreConverged();

      const bool converged = descend(muOverT);
      double loading = loadingOf();
      NLDFTBranchPoint point{.relativePressure = x, .chemicalPotential = muOverT};
      if (!converged)
      {
        ++result.unconverged;
        restoreConverged();
        loading = loadingOf();
      }
      else
      {
        point.grandPotential = grandPotentialOf(muOverT);
        point.residual = relativeResidual;
        point.stationary = true;
        keepConverged();
      }
      point.moleculesPerCell = loading;
      desorption.push_back(point);
    }
    std::reverse(desorption.begin(), desorption.end());
    result.desorptionIsotherm = std::move(desorption);
  }

  // The equilibrium isotherm: at each pressure, the state of lower grand potential. A point that is not a state at its
  // own pressure has no say, and where a branch has none the other answers.
  //
  // Only two different states are there to choose between. Away from a loop the sweeps are one curve, apart by no more
  // than the tolerance they converged to, and then the rounding in Omega has nothing to say and the sweep up is kept,
  // so that the isotherm does not hop between two readings of the same state.
  constexpr double sameStateTolerance = 1.0e-3;
  result.isotherm.reserve(adsorption.size());
  for (std::size_t p = 0; p < adsorption.size(); ++p)
  {
    const NLDFTBranchPoint &up = result.adsorptionIsotherm[p];
    const NLDFTBranchPoint &down = result.desorptionIsotherm[p];
    const bool apart = std::fabs(down.moleculesPerCell - up.moleculesPerCell) >
                       sameStateTolerance * std::max(up.moleculesPerCell, 1.0e-12);

    bool takeDown = false;
    if (up.stationary && down.stationary)
      takeDown = apart && down.grandPotential < up.grandPotential;
    else
      takeDown = down.stationary && !up.stationary;

    if (!up.stationary && !down.stationary) ++result.pointsWithoutASolution;

    const NLDFTBranchPoint &chosen = takeDown ? down : up;
    result.isotherm.push_back(IsothermPoint{chosen.relativePressure, chosen.moleculesPerCell});
    if (takeDown && apart && result.equilibriumTransition == 0.0) result.equilibriumTransition = chosen.relativePressure;
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
