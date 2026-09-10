module;

module energy_well_field;

import std;

import int3;
import uint3;
import double3;
import double3x3;
import unit_cell;
import crystal;
import pair_interactions;
import blocking_spheres;
import energy_shared_linear_probe;
import energy_shared_electrostatic_potential_grid;
import energy_shared_well_field;
import energy_shared_blocking_mask;
import structure_parallel;

namespace
{
double smoothPotentialAt(std::span<const float> smoothPotential, double3 s, uint3 gridSize)
{
  s.x -= std::floor(s.x);
  s.y -= std::floor(s.y);
  s.z -= std::floor(s.z);

  double3 scaled(s.x * static_cast<double>(gridSize.x), s.y * static_cast<double>(gridSize.y),
                 s.z * static_cast<double>(gridSize.z));

  std::size_t lowX = static_cast<std::size_t>(std::floor(scaled.x)) % gridSize.x;
  std::size_t lowY = static_cast<std::size_t>(std::floor(scaled.y)) % gridSize.y;
  std::size_t lowZ = static_cast<std::size_t>(std::floor(scaled.z)) % gridSize.z;

  double fracX = scaled.x - std::floor(scaled.x);
  double fracY = scaled.y - std::floor(scaled.y);
  double fracZ = scaled.z - std::floor(scaled.z);

  std::size_t highX = (lowX + 1) % gridSize.x;
  std::size_t highY = (lowY + 1) % gridSize.y;
  std::size_t highZ = (lowZ + 1) % gridSize.z;

  auto at = [&](std::size_t i, std::size_t j, std::size_t k)
  { return static_cast<double>(smoothPotential[(k * gridSize.y + j) * gridSize.x + i]); };
  auto mix = [](double a, double b, double t) { return a + t * (b - a); };

  double c00 = mix(at(lowX, lowY, lowZ), at(highX, lowY, lowZ), fracX);
  double c10 = mix(at(lowX, highY, lowZ), at(highX, highY, lowZ), fracX);
  double c01 = mix(at(lowX, lowY, highZ), at(highX, lowY, highZ), fracX);
  double c11 = mix(at(lowX, highY, highZ), at(highX, highY, highZ), fracX);

  return mix(mix(c00, c10, fracY), mix(c01, c11, fracY), fracZ);
}

struct Neighbourhood
{
  UnitCell unitCell{};
  std::vector<double3> positions{};
  std::vector<ProbeSite> sites{};
  std::vector<SitePair> pairs{};
  std::vector<double3> axes{};
  std::vector<double3> fractionalOffsets{};
  std::size_t numberOfOrientations{1};
  double thermalEnergy{0.0};
  double halfSpan{0.0};
  double extraReach{0.0};
  std::vector<BlockingSphere> spheres{};
  int3 shells{0, 0, 0};
  double cutOff{12.0};
  double cutOffSquared{144.0};
  double longestCutOff{12.0};
  double ceiling{0.0};
  double blockedEnergyPerAngstrom{0.0};
  double3 widths{};

  bool useCharges{false};
  double coulombCutOffSquared{0.0};
  double dummyCore2{wellDummyCoreRadius * wellDummyCoreRadius};
  ScreenedCoulomb screened{};
  std::span<const float> smoothPotential{};
  uint3 potentialGridSize{0, 0, 0};

  void widenFor(double distance)
  {
    this->extraReach = distance;
    double reach = this->longestCutOff + this->halfSpan + distance;
    this->shells = int3(static_cast<std::int32_t>(std::floor(reach / this->widths.x + 0.5)),
                        static_cast<std::int32_t>(std::floor(reach / this->widths.y + 0.5)),
                        static_cast<std::int32_t>(std::floor(reach / this->widths.z + 0.5)));
  }

  std::pair<double *, double *> scratch() const
  {
    static thread_local std::vector<double> buffer;
    if (buffer.size() < 2 * this->numberOfOrientations) buffer.assign(2 * this->numberOfOrientations, 0.0);
    return {buffer.data(), buffer.data() + this->numberOfOrientations};
  }

  const std::vector<NearbyImage> &gather(double3 centre) const
  {
    static thread_local std::vector<NearbyImage> nearby;
    nearby.clear();

    const double reach = this->longestCutOff + this->halfSpan + this->extraReach;
    const double reachSquared = reach * reach;
    for (std::size_t iatom = 0; iatom < this->positions.size(); ++iatom)
    {
      double3 ds = centre - this->positions[iatom];
      ds.x -= std::rint(ds.x);
      ds.y -= std::rint(ds.y);
      ds.z -= std::rint(ds.z);

      for (std::int32_t a = -this->shells.x; a <= this->shells.x; ++a)
      {
        for (std::int32_t b = -this->shells.y; b <= this->shells.y; ++b)
        {
          for (std::int32_t c = -this->shells.z; c <= this->shells.z; ++c)
          {
            double3 t = ds + double3(static_cast<double>(a), static_cast<double>(b), static_cast<double>(c));
            double far = std::max({std::abs(t.x) * this->widths.x, std::abs(t.y) * this->widths.y,
                                   std::abs(t.z) * this->widths.z});
            if (far > reach) continue;

            double3 dr = this->unitCell.cell * t;
            double rr = double3::dot(dr, dr);
            if (rr >= reachSquared) continue;
            nearby.push_back(NearbyImage{dr, iatom});
          }
        }
      }
    }
    return nearby;
  }

  template <bool WithCharges, bool WantClearance, bool WithLJ>
  void addPair(double rr, const SitePair &pair, bool dispersion, double &total, double &nearest,
               char *skipSmooth, char *dummyHit, char *overlapHit) const
  {
    if constexpr (WithLJ)
    {
      if (dispersion && rr < this->cutOffSquared)
      {
        double clamped = std::max(rr, 1.0e-8);
        double ratio = pair.sigma2 / clamped;
        double ratio3 = ratio * ratio * ratio;
        total += std::min(pair.epsilon4 * ratio3 * (ratio3 - 1.0) - pair.shift, this->ceiling);

        // r < σ is the TraPPE overlap (U_LJ > 0). rmin = 2^{1/6}σ is the well, not a wall.
        if (overlapHit != nullptr && rr < pair.sigma2) *overlapHit = 1;

        if constexpr (WantClearance) nearest = std::min(nearest, std::sqrt(clamped) - pair.contact);
      }
    }

    // Massless charge sites (N_com) have no VDW. A 1 Å floor is only to keep 1/r off the nucleus;
    // it is not a contact wall. Dispersing N atoms gate Coulomb at r < σ, not at rmin.
    const bool dummyCore = !dispersion && rr < this->dummyCore2;
    if (dummyCore)
    {
      if (dummyHit != nullptr) *dummyHit = 1;
      if (skipSmooth != nullptr) *skipSmooth = 1;
    }

    if constexpr (WithCharges)
    {
      if (!dummyCore && rr < this->coulombCutOffSquared && pair.chargeProduct != 0.0)
      {
        total += std::min(pair.chargeProduct * this->screened.at(rr), this->ceiling);
      }
    }
  }

  template <bool WithCharges, bool WantClearance, bool WithLJ>
  void walkAll(const std::vector<NearbyImage> &nearby, double3 shift, double *energies, double *clearances,
               char *skipSmooth, char *dummyHit, char *overlapHit, const char *skipElectro) const
  {
    const std::size_t numberOfSites = this->sites.size();
    const std::size_t orientations = this->numberOfOrientations;
    const double3 *directions = this->axes.data();

    std::fill(energies, energies + orientations, 0.0);
    if constexpr (WantClearance) std::fill(clearances, clearances + orientations, 1.0e10);
    if (dummyHit != nullptr) std::fill(dummyHit, dummyHit + orientations, 0);
    if (overlapHit != nullptr) std::fill(overlapHit, overlapHit + orientations, 0);

    double commonEnergy = 0.0;
    double commonClearance = 1.0e10;
    static thread_local std::vector<char> commonSkip;
    static thread_local std::vector<char> commonDummy;
    static thread_local std::vector<char> commonOverlap;
    commonSkip.assign(numberOfSites, 0);
    commonDummy.assign(numberOfSites, 0);
    commonOverlap.assign(numberOfSites, 0);

    const double outerReach = this->longestCutOff + this->halfSpan;
    const double outerReachSquared = outerReach * outerReach;

    for (const NearbyImage &image : nearby)
    {
      const double3 base = image.dr + shift;
      if (double3::dot(base, base) > outerReachSquared) continue;

      const SitePair *here = this->pairs.data() + image.atom * numberOfSites;

      for (std::size_t s = 0; s < numberOfSites; ++s)
      {
        const SitePair pair = here[s];
        const double offset = this->sites[s].offset;
        const bool dispersion = this->sites[s].dispersion;

        if (offset == 0.0)
        {
          this->addPair<WithCharges, WantClearance, WithLJ>(double3::dot(base, base), pair, dispersion, commonEnergy,
                                                            commonClearance, commonSkip.data() + s,
                                                            commonDummy.data() + s, commonOverlap.data() + s);
          continue;
        }

        for (std::size_t o = 0; o < orientations; ++o)
        {
          if (skipElectro != nullptr && skipElectro[o]) continue;
          double3 d = base + offset * directions[o];
          char *skip = skipSmooth == nullptr ? nullptr : skipSmooth + o * numberOfSites + s;
          char *dummy = dummyHit == nullptr ? nullptr : dummyHit + o;
          char *overlap = overlapHit == nullptr ? nullptr : overlapHit + o;
          this->addPair<WithCharges, WantClearance, WithLJ>(double3::dot(d, d), pair, dispersion, energies[o],
                                                            clearances[o], skip, dummy, overlap);
        }
      }
    }

    for (std::size_t o = 0; o < orientations; ++o)
    {
      if (skipElectro != nullptr && skipElectro[o]) continue;
      energies[o] += commonEnergy;
    }
    if constexpr (WantClearance)
    {
      for (std::size_t o = 0; o < orientations; ++o) clearances[o] = std::min(clearances[o], commonClearance);
    }
    if (skipSmooth != nullptr)
    {
      for (std::size_t o = 0; o < orientations; ++o)
      {
        for (std::size_t s = 0; s < numberOfSites; ++s)
        {
          if (commonSkip[s]) skipSmooth[o * numberOfSites + s] = 1;
        }
      }
    }
    if (dummyHit != nullptr)
    {
      for (std::size_t s = 0; s < numberOfSites; ++s)
      {
        if (!commonDummy[s]) continue;
        for (std::size_t o = 0; o < orientations; ++o) dummyHit[o] = 1;
      }
    }
    if (overlapHit != nullptr)
    {
      for (std::size_t s = 0; s < numberOfSites; ++s)
      {
        if (!commonOverlap[s]) continue;
        for (std::size_t o = 0; o < orientations; ++o) overlapHit[o] = 1;
      }
    }
  }

  void addSmoothPotential(double3 point, double *energies, const char *skipSmooth, const char *skipElectro) const
  {
    const std::size_t numberOfSites = this->sites.size();
    for (std::size_t o = 0; o < this->numberOfOrientations; ++o)
    {
      if (skipElectro != nullptr && skipElectro[o]) continue;
      const double3 *fractional = this->fractionalOffsets.data() + o * numberOfSites;
      double sum = 0.0;
      for (std::size_t s = 0; s < numberOfSites; ++s)
      {
        if (this->sites[s].charge == 0.0) continue;
        if (skipSmooth != nullptr && skipSmooth[o * numberOfSites + s]) continue;
        sum += this->sites[s].charge *
               smoothPotentialAt(this->smoothPotential, point + fractional[s], this->potentialGridSize);
      }
      energies[o] += sum;
    }
  }

  template <bool WantClearance>
  double reduce(const std::vector<NearbyImage> &nearby, double3 centre, double3 shiftFractional,
                double3 shiftCartesian, double &clearance, std::size_t &bestOrientation,
                float *orientationOut = nullptr) const
  {
    auto [energies, clearances] = this->scratch();

    // Soft TraPPE: rmin is the LJ well, not a wall. Keep U_LJ for r ≥ σ (including the inner
    // shoulder σ < r < rmin that GCMC samples). Gate Coulomb only at true overlap — either N
    // inside σ, or the dummy inside 1 Å — so a quadrupole cannot be evaluated with one N buried
    // and N_com still on. Overlapping poses keep the repulsive LJ; they are not ceilinged.
    static thread_local std::vector<char> dummyHit;
    static thread_local std::vector<char> overlapHit;
    static thread_local std::vector<char> skipElectro;
    dummyHit.assign(this->numberOfOrientations, 0);
    overlapHit.assign(this->numberOfOrientations, 0);
    skipElectro.assign(this->numberOfOrientations, 0);
    this->walkAll<false, true, true>(nearby, shiftCartesian, energies, clearances, nullptr, dummyHit.data(),
                                    overlapHit.data(), nullptr);
    for (std::size_t o = 0; o < this->numberOfOrientations; ++o)
    {
      skipElectro[o] = (dummyHit[o] || overlapHit[o]) ? 1 : 0;
    }

    if (this->useCharges)
    {
      static thread_local std::vector<double> lj;
      static thread_local std::vector<char> skipSmooth;
      lj.assign(energies, energies + this->numberOfOrientations);
      skipSmooth.assign(this->numberOfOrientations * this->sites.size(), 0);
      this->walkAll<true, false, false>(nearby, shiftCartesian, energies, clearances, skipSmooth.data(), nullptr,
                                       nullptr, skipElectro.data());
      this->addSmoothPotential(centre + shiftFractional, energies, skipSmooth.data(), skipElectro.data());
      for (std::size_t o = 0; o < this->numberOfOrientations; ++o)
      {
        if (dummyHit[o])
          energies[o] = this->ceiling;
        else
          energies[o] = std::min(lj[o] + energies[o], this->ceiling);
      }
    }
    else
    {
      for (std::size_t o = 0; o < this->numberOfOrientations; ++o)
      {
        if (dummyHit[o]) energies[o] = this->ceiling;
      }
    }

    double least = std::numeric_limits<double>::max();
    double sum = 0.0;
    double best = 1.0e10;
    bestOrientation = 0;
    const double beta = this->thermalEnergy > 0.0 ? 1.0 / this->thermalEnergy : 0.0;
    for (std::size_t o = 0; o < this->numberOfOrientations; ++o)
    {
      double energy = std::min(energies[o], this->ceiling);
      if (orientationOut != nullptr) orientationOut[o] = static_cast<float>(energy);

      if constexpr (WantClearance)
      {
        if (o == 0 || clearances[o] > best)
        {
          best = clearances[o];
          bestOrientation = o;
        }
      }

      if (beta > 0.0)
      {
        if (energy < least)
        {
          sum = sum * std::exp(-beta * (least - energy)) + 1.0;
          least = energy;
        }
        else
        {
          sum += std::exp(-beta * (energy - least));
        }
      }
      else
      {
        least = std::min(least, energy);
      }
    }

    clearance = best;
    if (beta > 0.0)
    {
      return std::min(
          least - this->thermalEnergy * std::log(sum / static_cast<double>(this->numberOfOrientations)),
          this->ceiling);
    }
    return least;
  }

  double energyAlong(const std::vector<NearbyImage> &nearby, double3 centre, double3 shiftFractional,
                     double3 shiftCartesian) const
  {
    double clearance = 0.0;
    std::size_t ignored = 0;
    return this->reduce<false>(nearby, centre, shiftFractional, shiftCartesian, clearance, ignored);
  }

  double3 softminDirection(const std::vector<NearbyImage> &nearby, std::size_t orientation, double nearest,
                           double &reliability) const
  {
    double3 directionSum{};
    double weightSum = 0.0;
    const std::size_t numberOfSites = this->sites.size();
    const double3 axis = this->axes[orientation];
    for (const NearbyImage &image : nearby)
    {
      const SitePair *here = this->pairs.data() + image.atom * numberOfSites;
      for (std::size_t s = 0; s < numberOfSites; ++s)
      {
        if (!this->sites[s].dispersion) continue;
        double3 d = image.dr + this->sites[s].offset * axis;
        double r = d.length();
        if (!(r > 1.0e-6)) continue;
        double weighted = r - here[s].contact;
        if (weighted - nearest > 6.0 * wellSoftminTau) continue;
        double w = std::exp(-(weighted - nearest) / wellSoftminTau);
        directionSum += (-d / r) * w;
        weightSum += w;
      }
    }
    if (!(weightSum > 0.0))
    {
      reliability = 0.0;
      return double3{};
    }
    reliability = directionSum.length() / weightSum;
    return double3::normalize(directionSum);
  }
};

Neighbourhood makeNeighbourhood(const PairInteractions &interactions, const Crystal &framework,
                                const LinearProbe &probe, std::size_t numberOfOrientations, double thermalEnergy,
                                std::span<const BlockingSphere> blockingSpheres, double blockedEnergyPerAngstrom,
                                double ceiling, const ElectrostaticPotentialGrid *potential, double coulombFactor)
{
  Neighbourhood neighbourhood;
  neighbourhood.unitCell = framework.unitCell;
  neighbourhood.spheres.assign(blockingSpheres.begin(), blockingSpheres.end());
  neighbourhood.blockedEnergyPerAngstrom = blockedEnergyPerAngstrom;
  neighbourhood.ceiling = ceiling;
  neighbourhood.thermalEnergy = thermalEnergy;
  neighbourhood.cutOff = interactions.cutOffVDW;
  neighbourhood.cutOffSquared = neighbourhood.cutOff * neighbourhood.cutOff;
  neighbourhood.widths = framework.unitCell.perpendicularWidths();

  neighbourhood.useCharges = probe.isCharged() && potential != nullptr && potential->numberOfVoxels() > 0;
  double coulombCutOff = 0.0;
  if (neighbourhood.useCharges)
  {
    if (coulombFactor == 0.0)
    {
      throw std::runtime_error(
          "Well field: the near half of the electrostatic sum needs the same charge-to-energy conversion "
          "the far half was built with\n");
    }
    coulombCutOff = potential->cutOff;
    neighbourhood.coulombCutOffSquared = coulombCutOff * coulombCutOff;
    neighbourhood.screened.build(potential->alpha, neighbourhood.coulombCutOffSquared);
    neighbourhood.smoothPotential = std::span<const float>(potential->smoothPotential);
    neighbourhood.potentialGridSize = potential->gridSize;
  }
  neighbourhood.longestCutOff = std::max(neighbourhood.cutOff, coulombCutOff);

  const std::size_t numberOfAtoms = framework.atoms.size();
  neighbourhood.positions.reserve(numberOfAtoms);
  for (std::size_t i = 0; i < numberOfAtoms; ++i)
  {
    neighbourhood.positions.push_back(framework.fractionalPositions.empty()
                                          ? framework.unitCell.inverseCell * framework.atoms[i].position
                                          : framework.fractionalPositions[i]);
  }

  std::vector<std::size_t> kept;
  for (std::size_t s = 0; s < probe.sites.size(); ++s)
  {
    const LinearProbe::Site &site = probe.sites[s];
    double charge = neighbourhood.useCharges ? site.charge : 0.0;
    bool carriesDispersion = false;
    for (std::size_t i = 0; i < numberOfAtoms; ++i)
    {
      carriesDispersion =
          carriesDispersion || interactions(site.type, framework.atoms[i].type).strengthParameter != 0.0;
    }
    if (!carriesDispersion && charge == 0.0) continue;

    kept.push_back(s);
    neighbourhood.sites.push_back(ProbeSite{site.offset, charge, carriesDispersion});
  }
  if (!std::ranges::any_of(neighbourhood.sites, &ProbeSite::dispersion))
  {
    throw std::runtime_error(
        std::format("Well field: probe '{}' has no site with any dispersion against this framework\n", probe.name));
  }

  const std::size_t numberOfSites = neighbourhood.sites.size();
  neighbourhood.pairs.resize(numberOfAtoms * numberOfSites);
  for (std::size_t i = 0; i < numberOfAtoms; ++i)
  {
    for (std::size_t s = 0; s < numberOfSites; ++s)
    {
      const PairParameters &pair = interactions(probe.sites[kept[s]].type, framework.atoms[i].type);
      neighbourhood.pairs[i * numberOfSites + s] =
          SitePair{.epsilon4 = 4.0 * pair.strengthParameter,
                   .sigma2 = pair.sizeParameter * pair.sizeParameter,
                   .contact = wellContactPrefactor * pair.sizeParameter,
                   .shift = pair.shift,
                   .chargeProduct = coulombFactor * neighbourhood.sites[s].charge * framework.atoms[i].charge};
    }
  }
  // Massless charge sites have no LJ. Their only electrostatic core is 1 Å: enough to keep 1/r off the
  // nucleus, not a contact wall. VDW overlap of the dispersing sites is decided before Coulomb is added.
  const double dummyCore2 = wellDummyCoreRadius * wellDummyCoreRadius;
  for (std::size_t i = 0; i < numberOfAtoms; ++i)
  {
    for (std::size_t s = 0; s < numberOfSites; ++s)
    {
      if (!neighbourhood.sites[s].dispersion) neighbourhood.pairs[i * numberOfSites + s].sigma2 = dummyCore2;
    }
  }

  double span = 0.0;
  for (const ProbeSite &site : neighbourhood.sites) span = std::max(span, std::abs(site.offset));
  neighbourhood.halfSpan = span;
  neighbourhood.axes = (span > 0.0 && numberOfOrientations > 1)
                           ? orientationSet(numberOfOrientations, probe.headTailSymmetric)
                           : std::vector<double3>{double3(0.0, 0.0, 1.0)};
  neighbourhood.numberOfOrientations = neighbourhood.axes.size();
  neighbourhood.fractionalOffsets.resize(neighbourhood.axes.size() * numberOfSites);
  for (std::size_t o = 0; o < neighbourhood.axes.size(); ++o)
  {
    for (std::size_t s = 0; s < numberOfSites; ++s)
    {
      neighbourhood.fractionalOffsets[o * numberOfSites + s] =
          framework.unitCell.inverseCell * (neighbourhood.sites[s].offset * neighbourhood.axes[o]);
    }
  }

  neighbourhood.widenFor(0.0);
  return neighbourhood;
}

double3 wrapUnit(double3 fractional)
{
  fractional.x -= std::floor(fractional.x);
  fractional.y -= std::floor(fractional.y);
  fractional.z -= std::floor(fractional.z);
  return fractional;
}

double3 quantizeFractional(double3 fractional)
{
  auto quantize = [](double x)
  {
    double r = std::rint(x * 1048576.0);
    if (r == 1048576.0) r = 0.0;
    return r / 1048576.0;
  };
  return double3(quantize(fractional.x), quantize(fractional.y), quantize(fractional.z));
}

double smoothstep(double edge0, double edge1, double x)
{
  if (edge1 <= edge0) return x >= edge1 ? 1.0 : 0.0;
  double t = std::clamp((x - edge0) / (edge1 - edge0), 0.0, 1.0);
  return t * t * (3.0 - 2.0 * t);
}


void refineVertices(std::vector<double3> &corners, std::vector<double> &energies, const Neighbourhood &neighbourhood,
                    float isovalue)
{
  energies.assign(corners.size(), 0.0);
  const double3x3 inverseCell = neighbourhood.unitCell.inverseCell;
  const float iso = isovalue;

  forEachIndex(corners.size(), workersAvailable(),
               [&](std::size_t, std::size_t vertex)
               {
                 double3 point = quantizeFractional(wrapUnit(corners[vertex]));

                 const std::vector<NearbyImage> &nearby = neighbourhood.gather(point);

                 double nearest = 0.0;
                 std::size_t bestOrientation = 0;
                 double energyHere =
                     neighbourhood.reduce<true>(nearby, point, double3{}, double3{}, nearest, bestOrientation);
                 energies[vertex] = energyHere;

                 if (energyHere > static_cast<double>(iso) - 0.02 * std::fabs(static_cast<double>(iso)))
                 {
                   return;
                 }

                 double pocket = blockingSphereDistance(point, neighbourhood.unitCell, neighbourhood.spheres);
                 if (pocket < nearest) return;

                 double reliability = 0.0;
                 double3 direction = neighbourhood.softminDirection(nearby, bestOrientation, nearest, reliability);
                 double span = 0.7 * smoothstep(0.25, 0.6, reliability);
                 if (span < 0.05 || direction.length_squared() < 1.0e-16) return;

                 double3 rayFractional = inverseCell * direction;

                 auto energyAlong = [&](double s)
                 { return neighbourhood.energyAlong(nearby, point, rayFractional * s, direction * s); };

                 constexpr int coarse = 14;
                 double sBest = 0.0;
                 double uBest = energyHere;
                 for (int i = -coarse; i <= coarse; ++i)
                 {
                   double s = span * static_cast<double>(i) / static_cast<double>(coarse);
                   double u = energyAlong(s);
                   if (u < uBest)
                   {
                     uBest = u;
                     sBest = s;
                   }
                 }
                 if (std::fabs(sBest) >= span - 0.5 * span / static_cast<double>(coarse)) return;

                 constexpr double invphi = 0.6180339887;
                 double a = sBest - span / static_cast<double>(coarse);
                 double b = sBest + span / static_cast<double>(coarse);
                 double x1 = b - invphi * (b - a);
                 double x2 = a + invphi * (b - a);
                 double f1 = energyAlong(x1);
                 double f2 = energyAlong(x2);
                 for (int iteration = 0; iteration < 20; ++iteration)
                 {
                   if (f1 < f2)
                   {
                     b = x2;
                     x2 = x1;
                     f2 = f1;
                     x1 = b - invphi * (b - a);
                     f1 = energyAlong(x1);
                   }
                   else
                   {
                     a = x1;
                     x1 = x2;
                     f1 = f2;
                     x2 = a + invphi * (b - a);
                     f2 = energyAlong(x2);
                   }
                 }
                 double s = 0.5 * (a + b);
                 corners[vertex] = wrapUnit(point + rayFractional * s);
                 energies[vertex] = energyAlong(s);
               });
}
}  // namespace


WellField WellFieldCPU::compute(const PairInteractions &interactions, const Crystal &framework,
                                const LinearProbe &probe, uint3 gridSize, std::size_t numberOfOrientations,
                                double thermalEnergy, std::span<const BlockingSphere> blockingSpheres,
                                double blockedEnergyPerAngstrom, double ceiling,
                                const ElectrostaticPotentialGrid *potential, double coulombFactor)
{
  if (gridSize.x == 0 || gridSize.y == 0 || gridSize.z == 0)
  {
    throw std::runtime_error("Well field: the grid must have at least one point along each axis\n");
  }
  if (potential != nullptr && potential->numberOfVoxels() > 0 &&
      (potential->gridSize.x != gridSize.x || potential->gridSize.y != gridSize.y ||
       potential->gridSize.z != gridSize.z))
  {
    throw std::runtime_error("Well field: the electrostatic potential is on a different grid than this one\n");
  }

  std::chrono::steady_clock::time_point time_begin = std::chrono::steady_clock::now();

  Neighbourhood neighbourhood = makeNeighbourhood(interactions, framework, probe, numberOfOrientations, thermalEnergy,
                                                  blockingSpheres, blockedEnergyPerAngstrom, ceiling, potential,
                                                  coulombFactor);

  WellField field;
  field.gridSize = gridSize;
  field.unitCell = framework.unitCell;
  field.probe = probe;
  field.probeName = probe.name;
  field.numberOfOrientations = neighbourhood.numberOfOrientations;
  field.thermalEnergy = thermalEnergy;
  field.cutOff = neighbourhood.cutOff;
  field.ceiling = ceiling;
  field.chargesIncluded = neighbourhood.useCharges;
  field.chargesIgnored = probe.isCharged() && !neighbourhood.useCharges;
  if (neighbourhood.useCharges)
  {
    field.ewaldAlpha = potential->alpha;
    field.numberOfWaveVectors = potential->numberOfWaveVectors;
  }

  const std::size_t numberOfVoxels = gridSize.x * gridSize.y * gridSize.z;
  field.energy.assign(numberOfVoxels, 0.0f);
  field.distance.assign(numberOfVoxels, 1.0e10f);
  field.reliability.assign(numberOfVoxels, 1.0f);
  const std::size_t nOrient = field.numberOfOrientations;
  // Cap is a soft memory guard for U(r, ω). 128³ × 128 orientations is ~1.07 GB of floats and
  // is the default --nldft --molecule N2 grid; the old 250M cutoff silently dropped that case
  // onto spherical ρ(Helmholtz). Allow up to ~4 GB (1e9 floats) on typical workstation RAM.
  if (nOrient > 1 && numberOfVoxels * nOrient <= 1'000'000'000)
    field.orientationEnergy.assign(numberOfVoxels * nOrient, 0.0f);

  const double invX = 1.0 / static_cast<double>(gridSize.x);
  const double invY = 1.0 / static_cast<double>(gridSize.y);
  const double invZ = 1.0 / static_cast<double>(gridSize.z);

  forEachBlock(gridSize.z, workersAvailable(),
               [&](std::size_t, std::size_t begin, std::size_t end)
               {
                 for (std::size_t iz = begin; iz < end; ++iz)
                 {
                   for (std::size_t iy = 0; iy < gridSize.y; ++iy)
                   {
                     for (std::size_t ix = 0; ix < gridSize.x; ++ix)
                     {
                       double3 fractional(static_cast<double>(ix) * invX, static_cast<double>(iy) * invY,
                                          static_cast<double>(iz) * invZ);

                       double pocket = blockingSphereDistance(fractional, neighbourhood.unitCell, neighbourhood.spheres);

                       const std::vector<NearbyImage> &nearby = neighbourhood.gather(fractional);

                       double bestClearance = 1.0e10;
                       std::size_t bestOrientation = 0;
                       const std::size_t voxel = (iz * gridSize.y + iy) * gridSize.x + ix;
                       float *orientationOut =
                           field.orientationEnergy.empty() ? nullptr : field.orientationEnergy.data() + voxel * nOrient;
                       double value = neighbourhood.reduce<true>(nearby, fractional, double3{}, double3{},
                                                                 bestClearance, bestOrientation, orientationOut);

                       double reliability = 1.0;
                       if (bestClearance < 1.0e9)
                       {
                         neighbourhood.softminDirection(nearby, bestOrientation, bestClearance, reliability);
                       }

                       double energy = pocket < 0.0
                                           ? std::min(-pocket * neighbourhood.blockedEnergyPerAngstrom, neighbourhood.ceiling)
                                           : value;
                       if (pocket < 0.0 && orientationOut != nullptr)
                       {
                         for (std::size_t o = 0; o < nOrient; ++o) orientationOut[o] = static_cast<float>(energy);
                       }
                       field.energy[voxel] = static_cast<float>(energy);
                       field.distance[voxel] = static_cast<float>(std::min(bestClearance, pocket));
                       field.reliability[voxel] = static_cast<float>(reliability);
                     }
                   }
                 }
               });

  std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - time_begin;
  field.seconds = elapsed.count();
  return field;
}


void WellFieldCPU::refineVertices(std::vector<double3> &corners, std::vector<double> &energies,
                                 const PairInteractions &interactions, const Crystal &framework,
                                 const NeighbourhoodParameters &parameters,
                                 std::span<const BlockingSphere> blockingSpheres, double iso)
{
  Neighbourhood neighbourhood =
      makeNeighbourhood(interactions, framework, parameters.probe, parameters.numberOfOrientations,
                        parameters.thermalEnergy, blockingSpheres, parameters.blockedEnergyPerAngstrom,
                        parameters.ceiling, parameters.potential, parameters.coulombFactor);
  neighbourhood.widenFor(parameters.extraReach > 0.0 ? parameters.extraReach : wellRefinementReach);
  ::refineVertices(corners, energies, neighbourhood, static_cast<float>(iso));
}
