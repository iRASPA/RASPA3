module;

module energy_shared_bet_surface_area;

import std;

import double3;
import double3x3;
import unit_cell;
import crystal;
import energy_shared_well_field;
import energy_shared_well_surface;

BETSurfaceArea::BETSurfaceArea() {}

BETSurfaceArea::~BETSurfaceArea() {}


namespace
{
// How many sites the filament is split into, in order of depth and with equal capacity each. Enough that
// the depth profile of a channel is resolved into a rise the BET line can be read off, few enough that
// every bin still carries a real share of the capacity rather than a single deep voxel.
constexpr std::size_t numberOfFilamentEnergyBins = 16;

// A volume per unit cell in Å³ over a molar mass in g/mol, read as a volume per gram in mL/g: Avogadro's
// number times the cubic Angstrom in cm³.
constexpr double cellVolumeToMillilitrePerGram = 0.60221419947;

// The n-layer BET occupancy per site (Brunauer, Emmett and Teller's finite form). At layers -> infinity it
// is the classical Cx/((1-x)(1-x+Cx)); at layers = 1 it collapses to Langmuir. Between pore walls the stack
// terminates when the pore is full, and it is that termination that turns the isotherm type I and stops
// n(1-x) rising at the filling step. `layers` need not be an integer: it interpolates the two forms.
double betOccupancy(double x, double C, double layers)
{
  if (!(x > 0.0) || x >= 1.0) return 0.0;
  if (!(layers > 0.0))
  {
    double denom = (1.0 - x) * (1.0 - x + C * x);
    if (!(denom > 0.0)) return 0.0;
    return C * x / denom;
  }
  double xN = std::pow(x, layers);
  double numerator = C * x * (1.0 - (layers + 1.0) * xN + layers * xN * x);
  double denominator = (1.0 - x) * (1.0 + (C - 1.0) * x - C * xN * x);
  if (!(std::fabs(denominator) > 0.0)) return 0.0;
  return numerator / denominator;
}

double langmuirOccupancy(double x, double C)
{
  if (!(x > 0.0)) return 0.0;
  return C * x / (1.0 + C * x);
}

// The relative pressures the model isotherm is tabulated on. The grid has to start below the knee or the
// fit is handed nothing but the plateau: in the Henry limit n = K x, so the pores fill around
// x = n_sat / K, and a strongly binding micropore puts that far below the 1e-5 that a free surface wants.
// Four decades of headroom covers the rise and the lowest candidate window, which reaches a thousandth of
// its ceiling. The conventional start stands wherever the knee is above it, so an open surface is
// tabulated exactly as before.
std::vector<double> relativePressureGrid(double saturation, double henryMoleculesPerCell)
{
  constexpr double conventionalStart = 1.0e-5;
  constexpr double end = 0.35;

  double start = conventionalStart;
  if (henryMoleculesPerCell > 0.0 && saturation > 0.0)
  {
    const double knee = saturation / henryMoleculesPerCell;
    start = std::clamp(std::min(conventionalStart, 1.0e-4 * knee), 1.0e-30, conventionalStart);
  }

  // Enough points to keep the resolution per decade roughly what 200 points over 4.5 decades used to give.
  const double decades = std::log10(end / start);
  const int numberOfPoints = std::clamp(static_cast<int>(std::lround(45.0 * decades)), 200, 4000);

  const double logMin = std::log(start);
  const double logMax = std::log(end);
  std::vector<double> grid;
  grid.reserve(static_cast<std::size_t>(numberOfPoints));
  for (int i = 0; i < numberOfPoints; ++i)
  {
    grid.push_back(
        std::exp(logMin + (logMax - logMin) * static_cast<double>(i) / static_cast<double>(numberOfPoints - 1)));
  }
  return grid;
}

struct LinearFit
{
  double slope{0.0};
  double intercept{0.0};
  double rSquared{0.0};
  bool ok{false};
};

LinearFit fitLine(std::span<const IsothermPoint> points, double xLow, double xHigh)
{
  LinearFit fit;
  double sumX = 0.0;
  double sumY = 0.0;
  double sumXX = 0.0;
  double sumXY = 0.0;
  double sumYY = 0.0;
  std::size_t count = 0;

  for (const IsothermPoint &point : points)
  {
    if (point.relativePressure < xLow || point.relativePressure > xHigh) continue;
    double n = point.moleculesPerCell;
    double oneMinus = 1.0 - point.relativePressure;
    if (!(n > 0.0) || !(oneMinus > 0.0)) continue;
    double x = point.relativePressure;
    double y = x / (n * oneMinus);
    if (!std::isfinite(y)) continue;
    sumX += x;
    sumY += y;
    sumXX += x * x;
    sumXY += x * y;
    sumYY += y * y;
    ++count;
  }

  if (count < 3) return fit;
  double inv = 1.0 / static_cast<double>(count);
  double meanX = sumX * inv;
  double meanY = sumY * inv;
  double varX = sumXX - static_cast<double>(count) * meanX * meanX;
  double covXY = sumXY - static_cast<double>(count) * meanX * meanY;
  if (!(std::fabs(varX) > 0.0)) return fit;
  fit.slope = covXY / varX;
  fit.intercept = meanY - fit.slope * meanX;
  double varY = sumYY - static_cast<double>(count) * meanY * meanY;
  if (varY > 0.0)
  {
    double r = covXY / std::sqrt(varX * varY);
    fit.rSquared = r * r;
  }
  fit.ok = std::isfinite(fit.slope) && std::isfinite(fit.intercept);
  return fit;
}

constexpr double angstromSquaredToSquareMetrePerMol = 6.0221419947e3;

// The Rouquerol window and the BET line, fitted to whatever isotherm sits in `result.isotherm`. Fills the
// fit fields (n_m, C, window, r², areas) and nothing else.
//
// The ceiling of any admissible window is measured off the isotherm rather than assumed: the largest x to
// which n(1-x) still rises, Rouquerol's first criterion. On an open surface that rise outlives the
// textbook 0.30 and the conventional window survives the cap; in a micropore it stops where the pores
// fill --- below 0.05 in the small IRMOFs of Walton and Snurr (JACS 129, 8552 (2007)) --- and the window
// has to follow it down, because a fit over the free-surface convention there straddles the filling step
// and returns a negative C. The remaining criteria are checked per candidate window: positive intercept
// and a positive capacity.
//
// Everything below the cap satisfies the first criterion already; the candidates scale with the cap so a
// micropore window keeps the shape of the conventional one (0.05-0.30 is a high of 0.30 with lows at 1/6,
// 1/10, 1/30 of it) at whatever height the isotherm allows. The deepest lows are for the windows of
// ultramicroporous materials, which reach several decades below their (already low) ceiling: Bae's MFI
// window spans 5e-5 to 1e-2. Among the admissible windows the straightest wins; among windows straight to
// within a factor of two of the best, the widest and highest wins, so a near-ideal isotherm reports the
// conventional window.
// The statistical thickness of the nitrogen film that a relative pressure grows on non-porous silica, in
// Angstrom: the reference a t-plot reads a porous sample against. The four branches are the fit Galarneau
// et al. give for a measured silica reference (Catalysis Research 2 (2022) 29, eqs 2-5). Outside the range
// they cover the reference says nothing, and neither does the t-plot.
std::optional<double> statisticalFilmThickness(double x)
{
  auto polynomial = [x](std::initializer_list<double> coefficients)
  {
    double value = 0.0;
    double power = 1.0;
    for (const double coefficient : coefficients)
    {
      value += coefficient * power;
      power *= x;
    }
    return value;
  };

  if (x < 0.009 || x > 0.90) return std::nullopt;
  if (x <= 0.125)
  {
    return polynomial({1.62973, 76.4748, -2171.7914, 41734.77357, -465290.41181, 2.72432e6, -6.43708e6});
  }
  if (x <= 0.60) return 3.07721 + 5.64019 * x;
  if (x <= 0.75)
  {
    return polynomial({4592.05803, -38117.31548, 131602.19741, -241680.40239, 249079.8569, -136632.44762, 31182.4149});
  }
  return polynomial({2098.4, -10711.0, 18954.0, -9197.5, -10624.0, 14046.0, -4553.0});
}

// The t-plot: adsorbed volume against the film thickness the same relative pressure would have grown on a
// non-porous surface. A micropore fills long before the reference grows its first layer, so the plot rises
// steeply and then flattens, and the standard reading takes the intercept of the flat stretch as the
// micropore volume and its slope as the surface a film is still free to grow on --- the mesopore plus
// external area. Experimental work reads it over 3.2 to 4.2 Å of thickness.
//
// On a periodic crystal there is no external surface at all, so the slope is not a measurement of one: it
// says how far the reading is disturbed by the filling itself. That is the bias experimental work has to
// calibrate away with reference materials, and here it can simply be looked at.
void applyTPlot(BETSurfaceArea &result, double mass, double liquidVolume)
{
  constexpr double lowestThickness = 3.2;
  constexpr double highestThickness = 4.2;

  if (!(mass > 0.0) || !(liquidVolume > 0.0)) return;

  std::vector<std::pair<double, double>> plot;
  for (const IsothermPoint &point : result.isotherm)
  {
    const std::optional<double> thickness = statisticalFilmThickness(point.relativePressure);
    if (!thickness.has_value()) continue;
    if (thickness.value() < lowestThickness || thickness.value() > highestThickness) continue;
    plot.emplace_back(thickness.value(), point.moleculesPerCell * liquidVolume * cellVolumeToMillilitrePerGram / mass);
  }
  if (plot.size() < 3) return;

  const double count = static_cast<double>(plot.size());
  double sumT = 0.0;
  double sumV = 0.0;
  double sumTT = 0.0;
  double sumTV = 0.0;
  for (const auto &[thickness, volume] : plot)
  {
    sumT += thickness;
    sumV += volume;
    sumTT += thickness * thickness;
    sumTV += thickness * volume;
  }
  const double determinant = count * sumTT - sumT * sumT;
  if (!(std::fabs(determinant) > 0.0)) return;

  const double slope = (count * sumTV - sumT * sumV) / determinant;
  const double intercept = (sumTT * sumV - sumT * sumTV) / determinant;

  double totalSquares = 0.0;
  double residualSquares = 0.0;
  const double mean = sumV / count;
  for (const auto &[thickness, volume] : plot)
  {
    totalSquares += (volume - mean) * (volume - mean);
    residualSquares += (volume - slope * thickness - intercept) * (volume - slope * thickness - intercept);
  }

  result.tPlotMicroporeVolume = intercept;
  // A slope in mL/g per Angstrom is an area: a film of thickness t on S m²/g holds S t of volume.
  result.tPlotExternalArea = 1.0e4 * slope;
  result.tPlotRSquared = (totalSquares > 0.0) ? 1.0 - residualSquares / totalSquares : 1.0;
  result.tPlotNumberOfPoints = plot.size();
}

void applyWindowFit(BETSurfaceArea &result, double mass, double cellVolume, double crossSection, double liquidVolume)
{
  if (result.isotherm.empty()) return;

  // The Gurvich reading: whatever the isotherm holds where it has stopped rising, taken as liquid.
  // It is the number reported next to every BET area, and unlike the area it rests on no fit --- the pore
  // either holds that much liquid or it does not.
  result.saturationLoading = 0.0;
  for (const IsothermPoint &point : result.isotherm)
  {
    result.saturationLoading = std::max(result.saturationLoading, point.moleculesPerCell);
  }
  if (mass > 0.0 && liquidVolume > 0.0)
  {
    result.microporeVolume =
        result.saturationLoading * liquidVolume * cellVolumeToMillilitrePerGram / mass;
  }
  applyTPlot(result, mass, liquidVolume);

  // Rouquerol's first criterion is that n(1-x) rises, so the ceiling of any admissible window is where it
  // stops rising: the maximum. Walking up until the rise breaks by a fixed per-step fraction instead makes
  // the ceiling an artifact of how finely the isotherm happens to be sampled. On a saturated plateau
  // n(1-x) falls by only x(r-1) between neighbouring points of a log grid of ratio r, so a per-step slack
  // of 1e-3 puts the ceiling near x = 1e-3/(r-1) whatever the material is --- around 0.02 here, decades
  // above where a micropore fills. Every candidate window then lies on the plateau, where x/(n(1-x)) is a
  // straight line through the origin, and the fit turns on the sign of an intercept that is numerically
  // zero: the same isotherm reads either the Gurvich capacity or nothing at all.
  double xCap = 0.0;
  {
    // n(1-x) is nonnegative, so -1 is below every value the scan can see (an infinity sentinel is not safe
    // under this build's fast-math flags).
    double best = -1.0;
    for (const IsothermPoint &point : result.isotherm)
    {
      best = std::max(best, point.moleculesPerCell * (1.0 - point.relativePressure));
    }
    // The first point to reach the maximum, so that noise out on the plateau of a simulated isotherm
    // cannot drag the ceiling past the knee.
    xCap = result.isotherm.front().relativePressure;
    for (const IsothermPoint &point : result.isotherm)
    {
      if (point.moleculesPerCell * (1.0 - point.relativePressure) >= best * 0.999)
      {
        xCap = point.relativePressure;
        break;
      }
    }
  }
  xCap = std::min(xCap, 0.30);

  double xLow = result.isotherm.front().relativePressure;
  double xHigh = xCap;
  LinearFit chosen{};
  bool found = false;

  struct Candidate
  {
    LinearFit fit;
    double low;
    double high;
  };
  std::vector<Candidate> admissible;

  for (double highFraction : {1.0, 0.75, 0.5, 0.35, 0.25})
  {
    double high = xCap * highFraction;
    for (double lowFraction : {1.0 / 6.0, 0.1, 1.0 / 30.0, 0.01, 1.0e-3})
    {
      double low = high * lowFraction;
      LinearFit fit = fitLine(result.isotherm, low, high);
      if (!fit.ok || !(fit.intercept > 0.0)) continue;
      if (!(1.0 / (fit.slope + fit.intercept) > 0.0)) continue;
      if (!(1.0 + fit.slope / fit.intercept > 0.0)) continue;
      // A monolayer cannot hold more than the pore holds when it is full. The condition costs nothing on
      // an open surface, where the stack keeps growing and the isotherm has no saturation to exceed, and
      // it is what keeps a window from being fitted to the filling step of a microporous isotherm: the
      // step is where n rises fastest, a line through it extrapolates to more than the pore contains, and
      // the resulting area is arbitrary rather than wrong-by-a-little. Faujasite's model isotherm was read
      // as 181 molecules per cell against the 112 it holds, and 1529 m²/g against a 975 m²/g reference.
      if (result.saturationLoading > 0.0 && 1.0 / (fit.slope + fit.intercept) > result.saturationLoading)
      {
        continue;
      }
      admissible.push_back(Candidate{fit, low, high});
    }
  }

  if (!admissible.empty())
  {
    double bestMisfit = 1.0;
    for (const Candidate &candidate : admissible)
    {
      bestMisfit = std::min(bestMisfit, 1.0 - candidate.fit.rSquared);
    }
    const double tolerated = 2.0 * bestMisfit + 1.0e-12;
    const Candidate *pick = nullptr;
    for (const Candidate &candidate : admissible)
    {
      if (1.0 - candidate.fit.rSquared > tolerated) continue;
      if (pick == nullptr || candidate.high > pick->high ||
          (candidate.high == pick->high && candidate.low < pick->low))
      {
        pick = &candidate;
      }
    }
    chosen = pick->fit;
    xLow = pick->low;
    xHigh = pick->high;
    result.monolayerCapacity = 1.0 / (chosen.slope + chosen.intercept);
    result.cConstant = 1.0 + chosen.slope / chosen.intercept;
    found = true;
  }

  if (!found)
  {
    // No admissible window at all: fit the whole rising region, which is Walton and Snurr's reading of the
    // criteria when the step and the monolayer cannot be told apart. A negative C is not a BET number.
    LinearFit fit = fitLine(result.isotherm, xLow, xCap);
    chosen = fit;
    xHigh = xCap;
    if (fit.ok && fit.intercept > 0.0 && std::fabs(fit.slope + fit.intercept) > 0.0)
    {
      const double c = 1.0 + fit.slope / fit.intercept;
      if (c > 0.0 && 1.0 / (fit.slope + fit.intercept) > 0.0)
      {
        result.monolayerCapacity = 1.0 / (fit.slope + fit.intercept);
        result.cConstant = c;
      }
    }
  }

  // Still nothing: no line through this isotherm is a BET line, and the reason is always the same. The
  // pore fills within about a decade of x, so the ceiling Rouquerol's first criterion puts on the window
  // sits at the top of the step, everything below it is the step, and across a step x/(n(1-x)) falls
  // rather than rises. Mordenite's model isotherm goes from one molecule per cell to its full nine
  // between x = 1.6e-8 and 2e-7 and is then flat to x = 1, and every candidate window in it fits a
  // negative monolayer.
  //
  // Such an isotherm is not unreadable, though: it is Type I, and Type I has an exact BET reading in the
  // limit of large C. The monolayer is the plateau, because the pore fills once and stops, and C follows
  // from the other end of the isotherm, the Henry slope n/x as x -> 0, which for BET is n_m C. Both ends
  // are properties of the isotherm rather than of a fitted window, so nothing here depends on where the
  // window would have been put.
  if (result.monolayerCapacity <= 0.0 && result.saturationLoading > 0.0)
  {
    const IsothermPoint &lowest = result.isotherm.front();
    result.monolayerCapacity = result.saturationLoading;
    result.cConstant = (lowest.relativePressure > 0.0)
                           ? lowest.moleculesPerCell / (lowest.relativePressure * result.saturationLoading)
                           : 0.0;
    result.plateauReading = true;
    xLow = lowest.relativePressure;
    xHigh = xCap;
    chosen.rSquared = 0.0;
  }

  result.windowLow = xLow;
  result.windowHigh = xHigh;
  result.rSquared = chosen.rSquared;

  double betArea = result.monolayerCapacity * crossSection;
  if (mass > 0.0)
  {
    result.gravimetricArea = betArea * angstromSquaredToSquareMetrePerMol / mass;
  }
  if (cellVolume > 0.0)
  {
    result.volumetricArea = 1.0e4 * betArea / cellVolume;
  }
}

// The grand potential of one lattice site at filling theta, in a field t = mu - U and a mean-field
// nitrogen-nitrogen coupling W (negative). The entropy term is the ideal lattice-gas mixing entropy.
double sitePotential(double theta, double t, double W, double kT)
{
  double mixing = theta * std::log(theta) + (1.0 - theta) * std::log(1.0 - theta);
  return -t * theta + 0.5 * W * theta * theta + kT * mixing;
}

// Bisection for the stationarity condition kT ln(theta/(1-theta)) + W theta = t on an interval where the
// left side is increasing.
double solveOccupancyBranch(double lo, double hi, double t, double W, double kT)
{
  for (int iteration = 0; iteration < 80; ++iteration)
  {
    double mid = 0.5 * (lo + hi);
    double g = kT * std::log(mid / (1.0 - mid)) + W * mid - t;
    if (g < 0.0)
      lo = mid;
    else
      hi = mid;
  }
  return 0.5 * (lo + hi);
}

// The equilibrium occupancy of one site: the global minimiser of the grand potential. Below the mean-field
// critical temperature (|W| > 4 kT) the stationarity condition has a dilute and a filled branch; the
// Maxwell construction keeps whichever has the lower potential, which is the condensation approximation.
double latticeOccupancy(double t, double W, double kT)
{
  constexpr double edge = 1.0e-18;
  if (!(std::fabs(W) > 4.0 * kT))
  {
    return solveOccupancyBranch(edge, 1.0 - edge, t, W, kT);
  }
  double disc = std::sqrt(1.0 - 4.0 * kT / std::fabs(W));
  double thetaLow = 0.5 * (1.0 - disc);   // upper spinodal of the dilute branch
  double thetaHigh = 0.5 * (1.0 + disc);  // lower spinodal of the filled branch
  auto g = [&](double theta) { return kT * std::log(theta / (1.0 - theta)) + W * theta - t; };
  bool hasDilute = g(thetaLow) >= 0.0;
  bool hasFilled = g(thetaHigh) <= 0.0;
  if (hasDilute && hasFilled)
  {
    double dilute = solveOccupancyBranch(edge, thetaLow, t, W, kT);
    double filled = solveOccupancyBranch(thetaHigh, 1.0 - edge, t, W, kT);
    return (sitePotential(dilute, t, W, kT) <= sitePotential(filled, t, W, kT)) ? dilute : filled;
  }
  if (hasFilled) return solveOccupancyBranch(thetaHigh, 1.0 - edge, t, W, kT);
  return solveOccupancyBranch(edge, thetaLow, t, W, kT);
}

// What share of the filament is a pore the contact sheet has not already drawn. The two meshes bound the
// same merged well: the sheet is the part of its boundary where the field still crosses the level, the
// filament overlay is the whole of it. Where the channel is wide enough for a sheet the two areas coincide
// and the medial curve inside threads the very molecules the sheet has counted; where the well has closed
// over into a tube there is no sheet left and the curve is the only reading of the pore. MFI is the first
// case (335 Å² of sheet against 347 Å² of filament boundary, so 3% is new) and ferrierite the second
// (4 Å² against 115 Å², so 96% is). Areas rather than capacities, because this is a question about which
// surface was drawn, not about how molecules pack on it.
double unsheetedShareOfFilament(double sheetArea, double filamentBoundaryArea)
{
  if (!(filamentBoundaryArea > 0.0)) return 1.0;

  return std::max(0.0, 1.0 - std::min(1.0, sheetArea / filamentBoundaryArea));
}

}  // namespace


GeometricSites geometricAdsorptionSites(std::span<const SheetPatch> patches,
                                        std::span<const FilamentVoxel> filament, double filamentBoundaryArea)
{
  GeometricSites result;

  double sheetArea = 0.0;
  for (const SheetPatch &patch : patches)
  {
    if (!(patch.area > 0.0)) continue;
    LatticeSite site;
    site.capacity = patch.area / nitrogenCrossSection;
    site.energy = (patch.energy[0] + patch.energy[1] + patch.energy[2]) / 3.0;
    result.sites.push_back(site);
    result.sheetCapacity += site.capacity;
    sheetArea += patch.area;
  }

  bool geometricPacking = false;
  for (const FilamentVoxel &voxel : filament)
  {
    if (voxel.area > 0.0 || voxel.length > 0.0) geometricPacking = true;
  }

  struct Contribution
  {
    double capacity{0.0};
    double energy{0.0};
  };
  std::vector<Contribution> contributions;
  contributions.reserve(filament.size());

  for (const FilamentVoxel &voxel : filament)
  {
    if (!(voxel.volume > 0.0) && !(voxel.length > 0.0) && !(voxel.area > 0.0)) continue;
    // A 2-D merged well packs on its midplane (area / 16.2), a 1-D file along its length (L / v_L^{1/3}),
    // and volume / v_L is what remains when the caller knows only a blob of centres. The tube of centres
    // is not a bound on any of this: it is the locus a molecule's centre may occupy, which for a single
    // file is thinner than the molecule by whatever the channel leaves it, so bounding a file by
    // volume / v_L undercounts by that factor --- ferrierite's channels hold four molecules per cell
    // through a tube of 25 Å³, which reads as 0.44.
    double cap = 0.0;
    if (voxel.area > 0.0)
    {
      cap = voxel.area / nitrogenCrossSection;
      result.filamentMidplaneCapacity += cap;
    }
    else if (voxel.length > 0.0)
    {
      cap = voxel.length / nitrogenLinearSpacing;
      result.filamentFileCapacity += cap;
    }
    else if (!geometricPacking)
    {
      cap = voxel.volume / nitrogenLiquidVolume;
      result.filamentBlobCapacity += cap;
    }
    if (!(cap > 0.0)) continue;
    result.filamentCapacity += cap;
    contributions.push_back(Contribution{cap, voxel.energy});
  }

  // Only the part of the pore the sheet has not already drawn is the filament's to count.
  const double unsheeted = unsheetedShareOfFilament(sheetArea, filamentBoundaryArea);
  result.filamentCapacity *= unsheeted;
  result.filamentFileCapacity *= unsheeted;
  result.filamentMidplaneCapacity *= unsheeted;
  result.filamentBlobCapacity *= unsheeted;
  for (Contribution &contribution : contributions) contribution.capacity *= unsheeted;

  // The filament is a stretch of channel, not one adsorption site. Its voxels run from the level the well
  // surface is drawn at down to the deepest point of the channel, which in ferrierite is seventeen kT of
  // spread, and the contact sheet is a chain of triangles with the same kind of spread. Put the whole
  // filament capacity at one energy and it fills within a kT of that energy, and a vertical step has no
  // BET line through it: x/(n(1-x)) drops two decades across the step, every window containing it fits a
  // negative monolayer, and the fit retreats to whatever plateau lies below --- 0.07 of ferrierite's 4.09
  // molecules per cell. Binning by equal capacity keeps the spread that broadens the step into something
  // an isotherm can be read off, while leaving every bin enough capacity that no single deep voxel carries
  // the Henry match by itself, which is what per-voxel sites used to do.
  if (result.filamentCapacity > 0.0 && !contributions.empty())
  {
    std::ranges::sort(contributions, {}, &Contribution::energy);

    const double capacityPerBin = result.filamentCapacity / static_cast<double>(numberOfFilamentEnergyBins);
    double binCapacity = 0.0;
    double binEnergySum = 0.0;
    for (std::size_t index = 0; index < contributions.size(); ++index)
    {
      double remaining = contributions[index].capacity;
      const double energy = contributions[index].energy;
      const bool last = (index + 1 == contributions.size());
      while (remaining > 0.0)
      {
        // A voxel that straddles a boundary is split between the two bins, so the bins hold equal capacity
        // however the voxel capacities happen to fall.
        const double room = capacityPerBin - binCapacity;
        const double taken = (!last && remaining > room) ? room : remaining;
        binCapacity += taken;
        binEnergySum += taken * energy;
        remaining -= taken;
        if (binCapacity >= capacityPerBin || (last && remaining <= 0.0))
        {
          LatticeSite site;
          site.capacity = binCapacity;
          site.energy = binEnergySum / binCapacity;
          result.sites.push_back(site);
          binCapacity = 0.0;
          binEnergySum = 0.0;
        }
      }
    }
  }

  return result;
}


GeometricSites fieldAdsorptionSites(const WellField &field, double packingDistance)
{
  GeometricSites result;
  if (field.energy.empty()) return result;

  const std::size_t nx = static_cast<std::size_t>(field.gridSize.x);
  const std::size_t ny = static_cast<std::size_t>(field.gridSize.y);
  const std::size_t nz = static_cast<std::size_t>(field.gridSize.z);
  if (nx == 0uz || ny == 0uz || nz == 0uz) return result;

  const double3x3 cell = field.unitCell.cell;
  const double packingSquared = packingDistance * packingDistance;

  // Deepest first. The molecule that enters at infinite dilution sits at the bottom of the field, and each
  // one after it takes the deepest point still open to it, so the order molecules are placed in is the
  // order they adsorb in and the energy each carries is the energy of the point it took. Only the
  // attractive region is offered: a voxel the field does not bind at cannot hold a molecule at 77 K, and
  // the blocked pockets are already sitting at the blocked energy, well above zero.
  std::vector<std::uint32_t> order;
  for (std::uint32_t voxel = 0u; voxel < static_cast<std::uint32_t>(field.energy.size()); ++voxel)
  {
    if (field.energy[voxel] < 0.0f) order.push_back(voxel);
  }
  if (order.empty()) return result;
  std::ranges::sort(order, {}, [&](std::uint32_t voxel) { return field.energy[voxel]; });

  // Bins at least a packing distance across in every direction, so a molecule too close to the point being
  // tried has to be in one of the twenty-seven bins around it. In a cell narrow enough to hold fewer than
  // three bins across an axis the same bin is reached more than once; the indices are made unique first.
  const double3 widths = field.unitCell.perpendicularWidths();
  const std::size_t bx = std::max(1uz, static_cast<std::size_t>(widths.x / packingDistance));
  const std::size_t by = std::max(1uz, static_cast<std::size_t>(widths.y / packingDistance));
  const std::size_t bz = std::max(1uz, static_cast<std::size_t>(widths.z / packingDistance));

  std::vector<double3> placed;
  std::vector<std::vector<std::uint32_t>> bins(bx * by * bz);
  std::vector<std::size_t> neighbourhood;
  neighbourhood.reserve(27uz);

  for (std::uint32_t voxel : order)
  {
    const std::size_t i = voxel % nx;
    const std::size_t j = (voxel / nx) % ny;
    const std::size_t k = voxel / (nx * ny);
    const double3 point(static_cast<double>(i) / static_cast<double>(nx),
                        static_cast<double>(j) / static_cast<double>(ny),
                        static_cast<double>(k) / static_cast<double>(nz));

    const std::size_t ib = std::min(bx - 1uz, static_cast<std::size_t>(point.x * static_cast<double>(bx)));
    const std::size_t jb = std::min(by - 1uz, static_cast<std::size_t>(point.y * static_cast<double>(by)));
    const std::size_t kb = std::min(bz - 1uz, static_cast<std::size_t>(point.z * static_cast<double>(bz)));

    neighbourhood.clear();
    for (std::size_t dk = 0uz; dk < 3uz; ++dk)
    {
      const std::size_t kk = (kb + bz + dk - 1uz) % bz;
      for (std::size_t dj = 0uz; dj < 3uz; ++dj)
      {
        const std::size_t jj = (jb + by + dj - 1uz) % by;
        for (std::size_t di = 0uz; di < 3uz; ++di)
        {
          const std::size_t ii = (ib + bx + di - 1uz) % bx;
          neighbourhood.push_back((kk * by + jj) * bx + ii);
        }
      }
    }
    std::ranges::sort(neighbourhood);
    const auto duplicates = std::ranges::unique(neighbourhood);
    neighbourhood.erase(duplicates.begin(), duplicates.end());

    bool clash = false;
    for (std::size_t bin : neighbourhood)
    {
      for (std::uint32_t other : bins[bin])
      {
        double3 ds = point - placed[other];
        ds.x -= std::round(ds.x);
        ds.y -= std::round(ds.y);
        ds.z -= std::round(ds.z);
        if ((cell * ds).length_squared() < packingSquared)
        {
          clash = true;
          break;
        }
      }
      if (clash) break;
    }
    if (clash) continue;

    const std::size_t bin = (kb * by + jb) * bx + ib;
    bins[bin].push_back(static_cast<std::uint32_t>(placed.size()));
    placed.push_back(point);
    result.sites.push_back(LatticeSite{1.0, static_cast<double>(field.energy[voxel])});
  }

  // Who attracts whom. Every image within the coupling shell counts, not the nearest one only: in a cell
  // narrower than twice the shell a molecule really does have two neighbours in the same direction, and
  // the crystal is what is being modelled. The pairs are few --- a few hundred molecules to a cell --- so
  // they are found by looking at all of them.
  const std::ptrdiff_t reachX = static_cast<std::ptrdiff_t>(std::ceil(nitrogenCouplingDistance / widths.x));
  const std::ptrdiff_t reachY = static_cast<std::ptrdiff_t>(std::ceil(nitrogenCouplingDistance / widths.y));
  const std::ptrdiff_t reachZ = static_cast<std::ptrdiff_t>(std::ceil(nitrogenCouplingDistance / widths.z));
  const double couplingSquared = nitrogenCouplingDistance * nitrogenCouplingDistance;

  result.neighbours.assign(placed.size(), {});
  for (std::size_t i = 0uz; i < placed.size(); ++i)
  {
    for (std::size_t j = 0uz; j < placed.size(); ++j)
    {
      const double3 ds = placed[i] - placed[j];
      for (std::ptrdiff_t sz = -reachZ; sz <= reachZ; ++sz)
      {
        for (std::ptrdiff_t sy = -reachY; sy <= reachY; ++sy)
        {
          for (std::ptrdiff_t sx = -reachX; sx <= reachX; ++sx)
          {
            const double3 shifted(ds.x + static_cast<double>(sx), ds.y + static_cast<double>(sy),
                                  ds.z + static_cast<double>(sz));
            const double rr = (cell * shifted).length_squared();
            if (rr > 1e-12 && rr < couplingSquared) result.neighbours[i].push_back(static_cast<std::uint32_t>(j));
          }
        }
      }
    }
  }

  // One molecule to a site, so the capacity is the count. The mesh split into sheet and filament has no
  // meaning here --- nothing was drawn --- and the whole of it is reported as the blob term, which is what
  // the volume reading of the pore has always been called.
  result.filamentBlobCapacity = static_cast<double>(result.sites.size());
  result.filamentCapacity = result.filamentBlobCapacity;

  return result;
}


LatticeIsotherm latticeGasIsotherm(std::span<const LatticeSite> sites, double henryMoleculesPerCell,
                                   double thermalEnergy, double heatOfLiquefaction,
                                   std::span<const std::vector<std::uint32_t>> neighbours)
{
  LatticeIsotherm model;
  if (sites.empty() || !(thermalEnergy > 0.0)) return model;

  const double kT = thermalEnergy;
  const double W = -2.0 * heatOfLiquefaction;

  // The bare lattice Henry slope, sum over sites of capacity exp((W/2 - U)/kT); the anchor scales the
  // activity so the model's x -> 0 slope is the exact grid limit. The shift absorbs the vibrational
  // prefactor the lattice entropy misses.
  double henryRaw = 0.0;
  double saturation = 0.0;
  for (const LatticeSite &site : sites)
  {
    if (!(site.capacity > 0.0)) continue;
    saturation += site.capacity;
    henryRaw += site.capacity * std::exp(std::min((0.5 * W - site.energy) / kT, 700.0));
  }
  const double f = (henryRaw > 0.0 && henryMoleculesPerCell > 0.0) ? henryMoleculesPerCell / henryRaw : 1.0;
  model.henryPrefactor = f;
  model.saturationCapacity = saturation;

  // mu(x) = W/2 + kT ln(f x): coexistence of bulk liquid (U = 0, theta -> 1) with its vapour at x = 1,
  // shifted by the Henry anchor.
  auto chemicalPotential = [&](double x) { return 0.5 * W + kT * std::log(f * x); };

  const bool coupled = (neighbours.size() == sites.size());

  auto totalAt = [&](double x)
  {
    const double mu = chemicalPotential(x);
    double n = 0.0;
    for (const LatticeSite &site : sites)
    {
      if (!(site.capacity > 0.0)) continue;
      n += site.capacity * latticeOccupancy(mu - site.energy, W, kT);
    }
    return n;
  };

  // The pair energy that makes a molecule with a liquid's worth of neighbours feel the liquid's binding:
  // the mean-field energy per molecule at coordination z_L is z_L epsilon / 2 = -q_L, which is also what
  // puts the coexistence of the bulk at mu = W/2, the same chemical potential the anchor above uses. A
  // site with fewer neighbours is bound by proportionately less and fills at a higher pressure.
  const double pairEnergy = W / nitrogenLiquidCoordination;

  auto relax = [&](std::vector<double> &theta, double mu)
  {
    constexpr std::size_t maximumIterations = 1000uz;
    constexpr double mixing = 0.4;
    constexpr double tolerance = 1e-10;
    for (std::size_t iteration = 0uz; iteration < maximumIterations; ++iteration)
    {
      double largestChange = 0.0;
      for (std::size_t i = 0uz; i < sites.size(); ++i)
      {
        double field = 0.0;
        for (std::uint32_t j : neighbours[i]) field += theta[j];
        const double argument = std::clamp((mu - sites[i].energy - pairEnergy * field) / kT, -700.0, 700.0);
        const double target = 1.0 / (1.0 + std::exp(-argument));
        largestChange = std::max(largestChange, std::fabs(target - theta[i]));
        theta[i] += mixing * (target - theta[i]);
      }
      if (largestChange < tolerance) break;
    }
  };

  // Below the transition the equations have both a dilute and a filled solution and the iteration keeps
  // whichever it started nearest; the grand potential says which is the state and which the metastable
  // one, so both starts are relaxed and compared. This is the Maxwell construction of the per-site model,
  // done once for the whole graph instead of site by site.
  auto grandPotential = [&](const std::vector<double> &theta, double mu)
  {
    double omega = 0.0;
    for (std::size_t i = 0uz; i < sites.size(); ++i)
    {
      const double t = std::clamp(theta[i], 1e-300, 1.0 - 1e-16);
      omega += (sites[i].energy - mu) * t;
      omega += kT * (t * std::log(t) + (1.0 - t) * std::log(1.0 - t));
      double field = 0.0;
      for (std::uint32_t j : neighbours[i]) field += theta[j];
      omega += 0.5 * pairEnergy * t * field;
    }
    return omega;
  };

  const std::vector<double> grid = relativePressureGrid(saturation, henryMoleculesPerCell);
  model.isotherm.reserve(grid.size());

  if (!coupled)
  {
    for (double x : grid) model.isotherm.push_back(IsothermPoint{x, totalAt(x)});
    return model;
  }

  // Each branch is followed from the end where it is the only solution --- the dilute one from the bottom
  // of the pressure range upwards, the filled one from the top downwards --- so that neither is lost to
  // the other on the way through the transition. Every point then has both, and the one with the lower
  // grand potential is the state.
  const std::size_t points = grid.size();
  std::vector<double> ascendingLoading(points, 0.0);
  std::vector<double> ascendingPotential(points, 0.0);
  std::vector<double> descendingLoading(points, 0.0);
  std::vector<double> descendingPotential(points, 0.0);

  auto loadingOf = [&](const std::vector<double> &theta)
  {
    double n = 0.0;
    for (std::size_t i = 0uz; i < sites.size(); ++i) n += sites[i].capacity * theta[i];
    return n;
  };

  std::vector<double> theta(sites.size(), 0.0);
  for (std::size_t k = 0uz; k < points; ++k)
  {
    const double mu = chemicalPotential(grid[k]);
    relax(theta, mu);
    ascendingLoading[k] = loadingOf(theta);
    ascendingPotential[k] = grandPotential(theta, mu);
  }

  theta.assign(sites.size(), 1.0);
  for (std::size_t k = points; k-- > 0uz;)
  {
    const double mu = chemicalPotential(grid[k]);
    relax(theta, mu);
    descendingLoading[k] = loadingOf(theta);
    descendingPotential[k] = grandPotential(theta, mu);
  }

  for (std::size_t k = 0uz; k < points; ++k)
  {
    const bool dilute = ascendingPotential[k] <= descendingPotential[k];
    model.isotherm.push_back(IsothermPoint{grid[k], dilute ? ascendingLoading[k] : descendingLoading[k]});
  }

  return model;
}


BETSurfaceArea BETSurfaceArea::fromIsotherm(std::vector<IsothermPoint> isotherm, double mass, double cellVolume,
                                            double crossSection, double liquidVolume)
{
  BETSurfaceArea result;
  result.isotherm = std::move(isotherm);
  applyWindowFit(result, mass, cellVolume, crossSection, liquidVolume);
  return result;
}

BETSurfaceArea BETSurfaceArea::fromIsothermFixedWindow(std::vector<IsothermPoint> isotherm, double mass,
                                                       double cellVolume, double windowLow, double windowHigh,
                                                       double crossSection, double liquidVolume)
{
  BETSurfaceArea result;
  result.isotherm = std::move(isotherm);
  if (result.isotherm.empty() || !(windowHigh > windowLow)) return result;

  result.saturationLoading = 0.0;
  for (const IsothermPoint &point : result.isotherm)
  {
    result.saturationLoading = std::max(result.saturationLoading, point.moleculesPerCell);
  }
  if (mass > 0.0 && liquidVolume > 0.0)
  {
    result.microporeVolume =
        result.saturationLoading * liquidVolume * cellVolumeToMillilitrePerGram / mass;
  }

  result.windowLow = windowLow;
  result.windowHigh = windowHigh;
  result.plateauReading = false;

  const LinearFit fit = fitLine(result.isotherm, windowLow, windowHigh);
  result.rSquared = fit.rSquared;
  if (!fit.ok || !(fit.intercept > 0.0) || !(std::fabs(fit.slope + fit.intercept) > 0.0))
  {
    return result;
  }
  const double monolayer = 1.0 / (fit.slope + fit.intercept);
  const double c = 1.0 + fit.slope / fit.intercept;
  if (!(monolayer > 0.0) || !(c > 0.0)) return result;
  if (result.saturationLoading > 0.0 && monolayer > result.saturationLoading) return result;

  result.monolayerCapacity = monolayer;
  result.cConstant = c;
  const double betArea = result.monolayerCapacity * crossSection;
  if (mass > 0.0)
  {
    result.gravimetricArea = betArea * angstromSquaredToSquareMetrePerMol / mass;
  }
  if (cellVolume > 0.0)
  {
    result.volumetricArea = 1.0e4 * betArea / cellVolume;
  }
  return result;
}


BETSurfaceArea BETSurfaceArea::fromSamples(std::span<const SheetPatch> patches,
                                           std::span<const FilamentVoxel> filament, double henryMoleculesPerCell,
                                           double thermalEnergy, double heatOfLiquefaction, double mass,
                                           double cellVolume, double multilayerRoom, double filamentBoundaryArea)
{
  BETSurfaceArea result;
  result.multilayerRoom = multilayerRoom;

  struct Site
  {
    double capacity{0.0};
    double rawC{0.0};
    bool multilayer{true};
  };

  std::vector<Site> sites;
  sites.reserve(patches.size() + filament.size());

  const double beta = (thermalEnergy > 0.0) ? 1.0 / thermalEnergy : 0.0;
  double sheetArea = 0.0;
  double sheetCapacity = 0.0;
  double sheetHenryRaw = 0.0;

  for (const SheetPatch &patch : patches)
  {
    if (!(patch.area > 0.0)) continue;
    Site site;
    site.capacity = patch.area / nitrogenCrossSection;
    site.multilayer = true;
    double energy = (patch.energy[0] + patch.energy[1] + patch.energy[2]) / 3.0;
    double exponent = std::min((-energy - heatOfLiquefaction) * beta, 700.0);
    site.rawC = (beta > 0.0) ? std::exp(exponent) : 1.0;
    sites.push_back(site);
    sheetArea += patch.area;
    sheetCapacity += site.capacity;
    sheetHenryRaw += site.capacity * site.rawC;
  }

  bool geometricPacking = false;
  for (const FilamentVoxel &voxel : filament)
  {
    if (voxel.area > 0.0 || voxel.length > 0.0) geometricPacking = true;
  }

  double filamentHenryRaw = 0.0;
  bool filamentIsPlanar = false;
  for (const FilamentVoxel &voxel : filament)
  {
    if (!(voxel.volume > 0.0) && !(voxel.length > 0.0) && !(voxel.area > 0.0)) continue;
    // A 2-D merged well packs on the midplane (area / 16.2). A 1-D file packs along its length
    // (L / v_L^{1/3}). Volume / v_L is what remains when the caller only knows a blob of centres
    // (unit tests, or a filament that never grew a mesh). The tube of centres bounds none of them: it is
    // where a centre may sit, which along a file is thinner than the molecule itself.
    double cap = 0.0;
    if (voxel.area > 0.0)
    {
      cap = voxel.area / nitrogenCrossSection;
      filamentIsPlanar = true;
      result.filamentMidplaneCapacity += cap;
    }
    else if (voxel.length > 0.0)
    {
      cap = voxel.length / nitrogenLinearSpacing;
      result.filamentFileCapacity += cap;
    }
    else if (!geometricPacking)
    {
      cap = voxel.volume / nitrogenLiquidVolume;
      result.filamentBlobCapacity += cap;
    }
    if (!(cap > 0.0)) continue;
    double exponent = std::min((-voxel.energy - heatOfLiquefaction) * beta, 700.0);
    double rawC = (beta > 0.0) ? std::exp(exponent) : 1.0;
    result.filamentCapacity += cap;
    filamentHenryRaw += cap * rawC;
  }

  // The sheet and the filament read one pore two ways, so the filament keeps only the share of it the
  // sheet has not drawn. Without that, MFI's medial curve adds 30 molecules per cell of channel centre
  // line to the 21 its channel walls already carry and the model isotherm saturates at more than twice
  // the loading an explicit simulation reaches. The Henry weight shrinks with the capacity, so C is
  // unchanged and only the room behind it moves.
  const double unsheeted = unsheetedShareOfFilament(sheetArea, filamentBoundaryArea);
  result.filamentCapacity *= unsheeted;
  result.filamentFileCapacity *= unsheeted;
  result.filamentMidplaneCapacity *= unsheeted;
  result.filamentBlobCapacity *= unsheeted;
  filamentHenryRaw *= unsheeted;

  if (result.filamentCapacity > 0.0)
  {
    // One merged-well site, one C. Per-voxel Langmuir lets the deepest voxel dominate henryRaw; f then
    // crushes every other voxel and the isotherm never saturates. A single site has C = Henry / n_fil
    // when there is no sheet, which is the Type I filling GCMC BET reads as n_m. Next to a contact
    // sheet the filament is more wall, not a deeper trap: using the sheet's C stops a core from
    // stealing the Henry match.
    Site site;
    site.capacity = result.filamentCapacity;
    site.multilayer = filamentIsPlanar && sheetArea > 0.0;
    site.rawC = (sheetCapacity > 0.0) ? sheetHenryRaw / sheetCapacity
                                      : filamentHenryRaw / result.filamentCapacity;
    sites.push_back(site);
  }

  double henryRaw = 0.0;
  for (const Site &site : sites) henryRaw += site.capacity * site.rawC;
  double f = (henryRaw > 0.0 && henryMoleculesPerCell > 0.0) ? henryMoleculesPerCell / henryRaw : 1.0;
  result.henryAnchor = f;

  auto nOfX = [&](double x)
  {
    double n = 0.0;
    for (const Site &site : sites)
    {
      double C = f * site.rawC;
      n += site.capacity * (site.multilayer ? betOccupancy(x, C, multilayerRoom) : langmuirOccupancy(x, C));
    }
    return n;
  };

  double saturation = 0.0;
  for (const Site &site : sites) saturation += site.capacity;

  const std::vector<double> grid = relativePressureGrid(saturation, henryMoleculesPerCell);
  result.isotherm.reserve(grid.size());
  for (double x : grid)
  {
    result.isotherm.push_back(IsothermPoint{x, nOfX(x)});
  }

  applyWindowFit(result, mass, cellVolume, nitrogenCrossSection, nitrogenLiquidVolume);

  auto gravimetricOf = [&](double molecules)
  {
    return (mass > 0.0) ? molecules * nitrogenCrossSection * angstromSquaredToSquareMetrePerMol / mass : 0.0;
  };
  result.sheetCapacity = sheetCapacity;
  result.sheetGravimetricArea = gravimetricOf(sheetCapacity);
  result.filamentFileGravimetricArea = gravimetricOf(result.filamentFileCapacity);
  result.filamentMidplaneGravimetricArea = gravimetricOf(result.filamentMidplaneCapacity);

  return result;
}
