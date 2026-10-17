module;

module interactions_neighbourhood;

import std;

import double3;
import double3x3;
import atom;
import running_energy;
import simulationbox;
import forcefield;
import vdwparameters;
import units;
import mc_cell_list;
import potential_pair_derivatives;
import potential_pair_vdw;
import potential_pair_coulomb;
import potential_ewald_real_space_table;
import simd_lanes;
import threadpool;

void Interactions::Neighbourhood::append(const Atom& atom, const double3& image, std::uint32_t recordedMoleculeId)
{
  x.push_back(image.x);
  y.push_back(image.y);
  z.push_back(image.z);
  charge.push_back(atom.charge);
  scalingVDW.push_back(atom.scalingVDW);
  scalingCoulomb.push_back(atom.scalingCoulomb);
  type.push_back(static_cast<std::uint32_t>(atom.type));
  groupId.push_back(static_cast<std::uint8_t>(atom.groupId));
  moleculeId.push_back(recordedMoleculeId);
  moleculeIdValue.push_back(static_cast<double>(recordedMoleculeId));
  if (atom.scalingVDW != 1.0 || atom.scalingCoulomb != 1.0 || atom.groupId != 0) fullyCoupled = false;
  ++count;
}

void Interactions::Neighbourhood::clear()
{
  count = 0;
  x.clear();
  y.clear();
  z.clear();
  charge.clear();
  scalingVDW.clear();
  scalingCoulomb.clear();
  type.clear();
  groupId.clear();
  moleculeId.clear();
  moleculeIdValue.clear();
  invalidateRows();
}

void Interactions::Neighbourhood::invalidateRows()
{
  // the row buffers keep their storage (no reallocation and no zero fill per gather): invalidated by type
  for (TypeRow& row : rows) row.type = paddingMoleculeId;
}

const Interactions::Neighbourhood::TypeRow& Interactions::Neighbourhood::lennardJonesRow(const ForceField& forceField,
                                                                                         std::uint32_t typeB) const
{
  TypeRow* free = nullptr;
  for (TypeRow& candidate : rows)
  {
    if (candidate.type == typeB) return candidate;
    if (candidate.type == paddingMoleculeId && !free) free = &candidate;
  }

  TypeRow& row = free ? *free : rows.emplace_back();
  row.type = typeB;
  const std::size_t padded = paddedSize();
  if (row.epsilon4.size() < padded)
  {
    row.epsilon4.resize(padded);
    row.sigma6.resize(padded);
    row.shift.resize(padded);
  }
  const std::size_t types = forceField.numberOfPseudoAtoms;
  const VDWParameters* parameters = forceField.data.data() + static_cast<std::size_t>(typeB) * types;
  for (std::size_t j = 0; j < padded; ++j)
  {
    const VDWParameters& pair = parameters[type[j]];
    if (pair.type == VDWParameters::Type::LennardJones)
    {
      const double sigma2 = pair.parameters.y * pair.parameters.y;
      row.epsilon4[j] = 4.0 * pair.parameters.x;
      row.sigma6[j] = sigma2 * sigma2 * sigma2;
      row.shift[j] = pair.shift;
    }
    else
    {
      row.epsilon4[j] = 0.0;
      row.sigma6[j] = 1.0;
      row.shift[j] = 0.0;
    }
  }
  return row;
}

void Interactions::Neighbourhood::pad()
{
  // records far from the centre (outside any cut-off) with the reserved molecule id; a lane kernel runs over whole
  // blocks and masks them by that id (and, with resolved images, by distance as well)
  constexpr double far = 1.0e6;
  const std::size_t padded = (count + padding - 1) / padding * padding;
  x.resize(padded, centre.x + far);
  y.resize(padded, centre.y + far);
  z.resize(padded, centre.z + far);
  charge.resize(padded, 0.0);
  scalingVDW.resize(padded, 1.0);
  scalingCoulomb.resize(padded, 1.0);
  type.resize(padded, 0u);
  groupId.resize(padded, static_cast<std::uint8_t>(0));
  moleculeId.resize(padded, paddingMoleculeId);
  moleculeIdValue.resize(padded, static_cast<double>(paddingMoleculeId));
}

void Interactions::Neighbourhood::gatherAll(std::span<const Atom> moleculeAtoms, double pairCutOff,
                                            std::optional<std::uint32_t> skipMolecule)
{
  gatherAllImpl(moleculeAtoms, pairCutOff, skipMolecule, Background::Molecules);
}

void Interactions::Neighbourhood::gatherAllFramework(std::span<const Atom> frameworkAtoms, double pairCutOff)
{
  gatherAllImpl(frameworkAtoms, pairCutOff, std::nullopt, Background::Framework);
}

void Interactions::Neighbourhood::gather(const SimulationBox& simulationBox, const MCCellList* cellList,
                                         std::span<const Atom> moleculeAtoms, const double3& batchCentre,
                                         double radiusOfBatch, double pairCutOff,
                                         std::optional<std::uint32_t> skipMolecule)
{
  gatherImpl(simulationBox, cellList, moleculeAtoms, batchCentre, radiusOfBatch, pairCutOff, skipMolecule,
             Background::Molecules);
}

void Interactions::Neighbourhood::gatherFramework(const SimulationBox& simulationBox, const MCCellList* cellList,
                                                  std::span<const Atom> frameworkAtoms, const double3& batchCentre,
                                                  double radiusOfBatch, double pairCutOff)
{
  gatherImpl(simulationBox, cellList, frameworkAtoms, batchCentre, radiusOfBatch, pairCutOff, std::nullopt,
             Background::Framework);
}

void Interactions::Neighbourhood::gatherAllImpl(std::span<const Atom> atoms, double pairCutOff,
                                                std::optional<std::uint32_t> skipMolecule, Background kind)
{
  background = kind;
  centre = double3(0.0, 0.0, 0.0);
  batchRadius = 0.0;
  cutOff = pairCutOff;
  nearestImageUnique = false;
  imagesResolved = false;
  fullyCoupled = true;
  allAtoms = true;

  clear();

  const bool framework = kind == Background::Framework;
  const std::uint32_t skipped = skipMolecule.value_or(paddingMoleculeId);
  for (const Atom& atom : atoms)
  {
    if (atom.moleculeId == skipped) continue;
    append(atom, atom.position, framework ? frameworkMoleculeId : atom.moleculeId);
  }
  numberOfCandidates = atoms.size();
  pad();
}

void Interactions::Neighbourhood::gatherFromCells(const SimulationBox& simulationBox, const MCCellList& cellList,
                                                  double radius, bool framework, std::uint32_t skipped)
{
  const std::array<bool, 3> whole = cellList.wholeAxes(radius);
  const bool anyWhole = whole[0] || whole[1] || whole[2];
  const double radiusSquared = radius * radius;
  const double3& widths = cellList.perpendicularWidths;
  const double3 inverseWidths(1.0 / widths.x, 1.0 / widths.y, 1.0 / widths.z);

  // Scratch across the cells: the slot and the position (relative to the centre) of every record kept. Sized
  // once to the atom count (an upper bound on the records of any gather), so the cells write by index.
  thread_local std::vector<std::uint32_t> keptSlot{};
  thread_local std::vector<double> keptX{}, keptY{}, keptZ{};
  if (keptSlot.size() < cellList.numberOfAtoms)
  {
    keptSlot.resize(cellList.numberOfAtoms);
    keptX.resize(cellList.numberOfAtoms);
    keptY.resize(cellList.numberOfAtoms);
    keptZ.resize(cellList.numberOfAtoms);
  }
  std::uint32_t* __restrict slots = keptSlot.data();
  double* __restrict kx = keptX.data();
  double* __restrict ky = keptY.data();
  double* __restrict kz = keptZ.data();

  const double* __restrict X = cellList.slotX.data();
  const double* __restrict Y = cellList.slotY.data();
  const double* __restrict Z = cellList.slotZ.data();
  const std::uint32_t* __restrict ids = cellList.slotMoleculeId.data();

  // Per-cell scratch for the displacements and squared distances of a bucket (the first pass streams through the
  // bucket without a dependency between the records; the second compacts the records kept).
  thread_local std::vector<double> bucketDX{}, bucketDY{}, bucketDZ{}, bucketRR{};
  if (bucketDX.size() < cellList.capacity)
  {
    bucketDX.resize(cellList.capacity);
    bucketDY.resize(cellList.capacity);
    bucketDZ.resize(cellList.capacity);
    bucketRR.resize(cellList.capacity);
  }
  double* __restrict dxs = bucketDX.data();
  double* __restrict dys = bucketDY.data();
  double* __restrict dzs = bucketDZ.data();
  double* __restrict rrs = bucketRR.data();

  // the minimum image along the whole axes only (a zero width leaves the axis alone: the shift is exact there)
  const bool perAxisImage = anyWhole && cellList.orthorhombic;
  const double wrapX = whole[0] ? widths.x : 0.0, wrapY = whole[1] ? widths.y : 0.0, wrapZ = whole[2] ? widths.z : 0.0;
  const double unwrapX = whole[0] ? inverseWidths.x : 0.0, unwrapY = whole[1] ? inverseWidths.y : 0.0,
               unwrapZ = whole[2] ? inverseWidths.z : 0.0;

  std::size_t candidates = 0;
  std::size_t m = 0;  // the records kept so far
  cellList.forEachCellWithin(centre, radius,
                             [&](std::size_t cell, const double3& shift)
                             {
                               const std::size_t base = cell * cellList.capacity;
                               const std::size_t n = cellList.counts[cell];
                               candidates += n;
                               const double ox = shift.x - centre.x, oy = shift.y - centre.y, oz = shift.z - centre.z;
                               const double* __restrict cx = X + base;
                               const double* __restrict cy = Y + base;
                               const double* __restrict cz = Z + base;

                               // first pass: the displacements to the centre
                               if (!anyWhole)
                               {
                                 for (std::size_t k = 0; k < n; ++k)
                                 {
                                   const double dx = cx[k] + ox, dy = cy[k] + oy, dz = cz[k] + oz;
                                   dxs[k] = dx;
                                   dys[k] = dy;
                                   dzs[k] = dz;
                                   rrs[k] = dx * dx + dy * dy + dz * dz;
                                 }
                               }
                               else if (perAxisImage)
                               {
                                 for (std::size_t k = 0; k < n; ++k)
                                 {
                                   double dx = cx[k] + ox, dy = cy[k] + oy, dz = cz[k] + oz;
                                   dx -= wrapX * std::round(dx * unwrapX);
                                   dy -= wrapY * std::round(dy * unwrapY);
                                   dz -= wrapZ * std::round(dz * unwrapZ);
                                   dxs[k] = dx;
                                   dys[k] = dy;
                                   dzs[k] = dz;
                                   rrs[k] = dx * dx + dy * dy + dz * dz;
                                 }
                               }
                               else
                               {
                                 for (std::size_t k = 0; k < n; ++k)
                                 {
                                   const double3 d = simulationBox.applyPeriodicBoundaryConditions(
                                       double3(cx[k] + ox, cy[k] + oy, cz[k] + oz));
                                   dxs[k] = d.x;
                                   dys[k] = d.y;
                                   dzs[k] = d.z;
                                   rrs[k] = double3::dot(d, d);
                                 }
                               }

                               // second pass: the records within the radius (the skipped molecule left out), written by
                               // index
                               for (std::size_t k = 0; k < n; ++k)
                               {
                                 const bool keep = (rrs[k] <= radiusSquared) & (ids[base + k] != skipped);
                                 slots[m] = static_cast<std::uint32_t>(base + k);
                                 kx[m] = dxs[k];
                                 ky[m] = dys[k];
                                 kz[m] = dzs[k];
                                 m += keep;
                               }
                             });
  numberOfCandidates = candidates;
  if (m == 0) return;

  // one block copy of the records kept out of the list's arrays
  count = m;
  x.resize(count);
  y.resize(count);
  z.resize(count);
  charge.resize(count);
  scalingVDW.resize(count);
  scalingCoulomb.resize(count);
  type.resize(count);
  groupId.resize(count);
  moleculeId.resize(count);
  moleculeIdValue.resize(count);
  for (std::size_t i = 0; i < m; ++i)
  {
    x[i] = centre.x + kx[i];
    y[i] = centre.y + ky[i];
    z[i] = centre.z + kz[i];
  }
  for (std::size_t i = 0; i < m; ++i)
  {
    const std::size_t slot = slots[i];
    charge[i] = cellList.slotCharge[slot];
    type[i] = cellList.slotType[slot];
  }
  if (cellList.numberOfPartiallyCoupledAtoms == 0)
  {
    // every atom of the list is fully coupled: the constants, no per-slot copy
    std::fill_n(scalingVDW.data(), m, 1.0);
    std::fill_n(scalingCoulomb.data(), m, 1.0);
    std::fill_n(groupId.data(), m, static_cast<std::uint8_t>(0));
  }
  else
  {
    bool coupled = true;
    for (std::size_t i = 0; i < m; ++i)
    {
      const std::size_t slot = slots[i];
      const double sVDW = cellList.slotScalingVDW[slot];
      const double sCoulomb = cellList.slotScalingCoulomb[slot];
      const std::uint8_t group = cellList.slotGroupId[slot];
      scalingVDW[i] = sVDW;
      scalingCoulomb[i] = sCoulomb;
      groupId[i] = group;
      coupled &= (sVDW == 1.0) & (sCoulomb == 1.0) & (group == 0);
    }
    if (!coupled) fullyCoupled = false;
  }
  if (framework)
  {
    std::fill_n(moleculeId.data(), m, frameworkMoleculeId);
    std::fill_n(moleculeIdValue.data(), m, static_cast<double>(frameworkMoleculeId));
  }
  else
  {
    for (std::size_t i = 0; i < m; ++i)
    {
      const std::size_t slot = slots[i];
      moleculeId[i] = ids[slot];
      moleculeIdValue[i] = cellList.slotMoleculeIdValue[slot];
    }
  }
}

void Interactions::Neighbourhood::gatherImpl(const SimulationBox& simulationBox, const MCCellList* cellList,
                                             std::span<const Atom> atoms, const double3& batchCentre,
                                             double radiusOfBatch, double pairCutOff,
                                             std::optional<std::uint32_t> skipMolecule, Background kind)
{
  background = kind;
  centre = batchCentre;
  batchRadius = radiusOfBatch;
  cutOff = pairCutOff;
  const double radius = batchRadius + cutOff;
  const double3 widths = simulationBox.perpendicularWidths();
  const double minimumWidth = std::min({widths.x, widths.y, widths.z});
  const double halfWidth = 0.5 * minimumWidth;
  nearestImageUnique = radius < halfWidth;
  // beyond half the box several images of an atom can lie within the gather radius; all of them are gathered
  // when the cut-off itself keeps at most one image of an atom within reach of a trial position, and the
  // radius stays below the box width (the images 'forEachImageWithin' enumerates -- the nearest one and its
  // neighbours across the box -- are then all of them; a wider batch, e.g. a displacement over the whole box,
  // falls back to the nearest image per atom with the minimum image taken per pair in the kernels)
  const bool allImages = !nearestImageUnique && cutOff <= halfWidth && radius < minimumWidth;
  imagesResolved = nearestImageUnique || allImages;
  fullyCoupled = true;
  allAtoms = false;

  const bool framework = kind == Background::Framework;
  const std::uint32_t skipped = skipMolecule.value_or(paddingMoleculeId);

  if (cellList && cellList->enabled && cellList->valid && cellList->numberOfAtoms == atoms.size())
  {
    // Block gather: every stencil cell within reach is filtered and copied as a block out of the list's
    // structure-of-arrays slots. The records of a cell are placed at 'position + shift' (exact along the axes
    // the stencil does not wrap) and at the minimum image to the centre along the whole axes; so the records
    // are the images the pairs use when the nearest image is unique, and otherwise the kernels take the
    // minimum image per pair (the cut-off then exceeds half the box along some axis). The record arrays are
    // resized in place (not cleared first: a resize from empty would zero-fill them before the copy).
    imagesResolved = nearestImageUnique;
    count = 0;
    invalidateRows();
    gatherFromCells(simulationBox, *cellList, radius, framework, skipped);
    pad();
    return;
  }

  clear();

  // no usable cell list: every atom is a candidate, placed with the minimum image to the centre (and its further
  // images when the radius reaches them)
  for (const Atom& atom : atoms)
  {
    if (atom.moleculeId == skipped) continue;
    const std::uint32_t recorded = framework ? frameworkMoleculeId : atom.moleculeId;
    const double3 d = simulationBox.applyPeriodicBoundaryConditions(atom.position - centre);
    MCCellList::forEachImageWithin(simulationBox, widths, d, radius, allImages,
                                   [&](const double3& image) { append(atom, centre + image, recorded); });
  }
  numberOfCandidates = atoms.size();
  pad();
}

bool Interactions::Neighbourhood::covers(std::span<const Atom> atoms, double cutOffVDW, double cutOffCoulomb,
                                         bool useCharge) const
{
  const double queryCutOff = std::max(cutOffVDW, useCharge ? cutOffCoulomb : 0.0);
  if (queryCutOff > cutOff) return false;
  if (allAtoms) return true;
  // a hair of tolerance for the rounding of the centroid
  const double limit = batchRadius + 1e-9 * std::max(1.0, batchRadius);
  const double limitSquared = limit * limit;
  for (const Atom& atom : atoms)
  {
    const double3 d = atom.position - centre;
    if (double3::dot(d, d) > limitSquared) return false;
  }
  return true;
}

std::pair<double3, double> Interactions::Neighbourhood::boundingSphere(std::span<const std::span<const Atom>> sets)
{
  double3 sum(0.0, 0.0, 0.0);
  std::size_t count = 0;
  for (std::span<const Atom> set : sets)
  {
    for (const Atom& atom : set)
    {
      sum += atom.position;
      ++count;
    }
  }
  if (count == 0) return {double3(0.0, 0.0, 0.0), 0.0};
  const double3 centre = sum / static_cast<double>(count);
  double radiusSquared = 0.0;
  for (std::span<const Atom> set : sets)
  {
    for (const Atom& atom : set)
    {
      const double3 d = atom.position - centre;
      radiusSquared = std::max(radiusSquared, double3::dot(d, d));
    }
  }
  return {centre, std::sqrt(radiusSquared)};
}

// ---------------------------------------------------------------------------------------------------
// The generic kernel: any pair potential, soft-core scaling, dU/dlambda; the summation order of the brute
// force (one background atom after the other per trial atom).
// ---------------------------------------------------------------------------------------------------

namespace
{
/// Books the summed pair energies as molecule-molecule or framework-molecule terms, after the background.
inline void bookPairEnergies(const Interactions::Neighbourhood& neighbourhood, RunningEnergy& energySum,
                             double energyVDW, double energyCharge)
{
  if (neighbourhood.background == Interactions::Neighbourhood::Background::Framework)
  {
    energySum.frameworkMoleculeVDW = energyVDW;
    energySum.frameworkMoleculeCharge = energyCharge;
    return;
  }
  energySum.moleculeMoleculeVDW = energyVDW;
  energySum.moleculeMoleculeCharge = energyCharge;
}

template <bool ResolvedImages>
[[nodiscard]] std::optional<RunningEnergy> genericNeighbourhoodEnergy(
    const ForceField& forceField, const SimulationBox& simulationBox, const Interactions::Neighbourhood& neighbourhood,
    std::span<const Atom> atoms, double cutOffVDW, double cutOffCoulomb, double sign, bool checkOverlap) noexcept
{
  const bool useCharge = forceField.useCharge;
  const double overlapCriteria = forceField.energyOverlapCriteria;
  const double cutOffVDWSquared = cutOffVDW * cutOffVDW;
  const double cutOffChargeSquared = cutOffCoulomb * cutOffCoulomb;

  const std::size_t n = neighbourhood.size();
  const double* __restrict xs = neighbourhood.x.data();
  const double* __restrict ys = neighbourhood.y.data();
  const double* __restrict zs = neighbourhood.z.data();
  const double* __restrict charges = neighbourhood.charge.data();
  const double* __restrict scalingsVDW = neighbourhood.scalingVDW.data();
  const double* __restrict scalingsCoulomb = neighbourhood.scalingCoulomb.data();
  const std::uint32_t* __restrict types = neighbourhood.type.data();
  const std::uint8_t* __restrict groupIds = neighbourhood.groupId.data();
  const std::uint32_t* __restrict moleculeIds = neighbourhood.moleculeId.data();

  RunningEnergy energySum{};
  double energyVDW = 0.0;
  double energyCharge = 0.0;

  for (const Atom& atom : atoms)
  {
    const std::uint32_t molB = atom.moleculeId;
    const double3 posB = atom.position;
    const std::size_t typeB = static_cast<std::size_t>(atom.type);
    const std::uint8_t groupIdB = atom.groupId;
    const double scalingVDWB = atom.scalingVDW;
    const double scalingCoulombB = atom.scalingCoulomb;
    const double chargeB = atom.charge;

    for (std::size_t j = 0; j < n; ++j)
    {
      if (moleculeIds[j] == molB) continue;

      double3 dr(xs[j] - posB.x, ys[j] - posB.y, zs[j] - posB.z);
      if constexpr (!ResolvedImages)
      {
        dr = simulationBox.applyPeriodicBoundaryConditions(dr);
      }
      const double rr = double3::dot(dr, dr);

      if (rr < cutOffVDWSquared)
      {
        const Potentials::PairDerivatives<0> factors = Potentials::potentialVDW<0>(
            forceField, scalingsVDW[j], scalingVDWB, rr, static_cast<std::size_t>(types[j]), typeB);
        if (checkOverlap && factors.energy > overlapCriteria) return std::nullopt;
        energyVDW += factors.energy;
        energySum.addDudlambdaVDW(groupIds[j], groupIdB, scalingsVDW[j], scalingVDWB, sign * factors.dUdlambda);
      }
      if (useCharge && rr < cutOffChargeSquared)
      {
        if (charges[j] * chargeB != 0.0 && scalingsCoulomb[j] * scalingCoulombB != 0.0)
        {
          const double r = std::sqrt(rr);
          const Potentials::PairDerivatives<0> factors =
              Potentials::potentialCoulomb<0>(forceField, scalingsCoulomb[j], scalingCoulombB, r, charges[j], chargeB);
          energyCharge += factors.energy;
          energySum.addDudlambdaCharge(groupIds[j], groupIdB, scalingsCoulomb[j], scalingCoulombB,
                                       sign * factors.dUdlambda);
        }
      }
    }
  }

  bookPairEnergies(neighbourhood, energySum, sign * energyVDW, sign * energyCharge);
  return energySum;
}

// ---------------------------------------------------------------------------------------------------
// The lane kernel: plain 12-6 Lennard-Jones (optionally shifted) plus the tabulated Ewald real-space term
// between fully coupled atoms, 'lanes' background records at a time. The same arithmetic as the cluster kernel
// of the spatial-decomposition MD engine (energies only, no forces).
// ---------------------------------------------------------------------------------------------------

/// Two SIMD registers of doubles per operation (4 lanes on NEON / SSE, 8 on AVX2): enough independent work per
/// iteration to hide the latency of the divisions and gathers.
constexpr std::size_t lanes = 2 * SIMD::simdBytes / sizeof(double);
constexpr std::uint32_t paddingMoleculeId = Interactions::Neighbourhood::paddingMoleculeId;
static_assert(Interactions::Neighbourhood::padding % lanes == 0, "the neighbourhood padding must cover a block");
using Vec = SIMD::Lanes<double, lanes>;

/// How the kernel treats the periodic images of the records.
enum class Images
{
  Resolved,     ///< The records are the images to use: no minimum image in the kernel.
  Rectangular,  ///< The minimum image per lane, per axis (a rectangular box).
  General,      ///< The minimum image per lane in the general cell-matrix form.
};

/// The minimum image per lane: for a rectangular box the per-axis form of
/// SimulationBox::applyPeriodicBoundaryConditions (to rounding), otherwise the general cell-matrix form.
template <Images Mode>
[[nodiscard]] std::optional<RunningEnergy> laneNeighbourhoodEnergy(
    const ForceField& forceField, const SimulationBox& simulationBox, const Interactions::Neighbourhood& neighbourhood,
    std::span<const Atom> atoms, double cutOffVDW, double cutOffCoulomb, double sign, bool checkOverlap) noexcept
{
  const bool useCharge = forceField.useCharge;
  const double3x3& cell = simulationBox.cell;
  const double3x3& inverseCell = simulationBox.inverseCell;
  const double lengthX = cell.ax, lengthY = cell.by, lengthZ = cell.cz;
  const double inverseLengthX = inverseCell.ax, inverseLengthY = inverseCell.by, inverseLengthZ = inverseCell.cz;
  const double overlapCriteria = forceField.energyOverlapCriteria;
  const double cutOffVDWSquared = cutOffVDW * cutOffVDW;
  const double cutOffChargeSquared = cutOffCoulomb * cutOffCoulomb;
  const double coulombFactor = Units::CoulombicConversionFactor;

  const std::size_t padded = neighbourhood.paddedSize();
  const double* __restrict xs = neighbourhood.x.data();
  const double* __restrict ys = neighbourhood.y.data();
  const double* __restrict zs = neighbourhood.z.data();
  const double* __restrict charges = neighbourhood.charge.data();
  const double* __restrict moleculeIds = neighbourhood.moleculeIdValue.data();
  constexpr double paddingId = static_cast<double>(paddingMoleculeId);

  const EwaldRealSpaceTable& table = forceField.ewaldRealSpaceTable;
  const double* __restrict tableCoefficients = table.coefficients.data();
  const double tableInverseSpacing = table.inverseSpacing;
  const double tableRRMin = table.rrMin;
  const double tableMaximumPosition = static_cast<double>(table.value.size() >= 2 ? table.value.size() - 2 : 0);

  // pairs are never this close between molecules; the clamp keeps the masked lanes (same molecule, padding)
  // finite
  constexpr double minimumRR = 1.0e-4;

  double energyVDW = 0.0;
  double energyCharge = 0.0;

  for (const Atom& atom : atoms)
  {
    const Interactions::Neighbourhood::TypeRow& row =
        neighbourhood.lennardJonesRow(forceField, static_cast<std::uint32_t>(atom.type));
    const double* __restrict rowEpsilon4 = row.epsilon4.data();
    const double* __restrict rowSigma6 = row.sigma6.data();
    const double* __restrict rowShift = row.shift.data();

    const double molB = static_cast<double>(atom.moleculeId);
    const double xB = atom.position.x, yB = atom.position.y, zB = atom.position.z;
    const double chargeB = atom.charge;
    const bool chargedB = useCharge && chargeB != 0.0;

    Vec eVDW = Vec::zero();
    Vec eCharge = Vec::zero();
    Vec maximumVDW = Vec::zero();
    Vec anyBelow = Vec::zero();  // an interacting lane closer than the table range was met

    for (std::size_t j = 0; j < padded; j += lanes)
    {
      Vec dx = Vec::load(xs + j) - xB;
      Vec dy = Vec::load(ys + j) - yB;
      Vec dz = Vec::load(zs + j) - zB;
      if constexpr (Mode == Images::Rectangular)
      {
        dx -= lengthX * round(dx * inverseLengthX);
        dy -= lengthY * round(dy * inverseLengthY);
        dz -= lengthZ * round(dz * inverseLengthZ);
      }
      else if constexpr (Mode == Images::General)
      {
        Vec sx = inverseCell.m11 * dx + inverseCell.m12 * dy + inverseCell.m13 * dz;
        Vec sy = inverseCell.m21 * dx + inverseCell.m22 * dy + inverseCell.m23 * dz;
        Vec sz = inverseCell.m31 * dx + inverseCell.m32 * dy + inverseCell.m33 * dz;
        sx -= round(sx);
        sy -= round(sy);
        sz -= round(sz);
        dx = cell.m11 * sx + cell.m12 * sy + cell.m13 * sz;
        dy = cell.m21 * sx + cell.m22 * sy + cell.m23 * sz;
        dz = cell.m31 * sx + cell.m32 * sy + cell.m33 * sz;
      }
      const Vec rr = dx * dx + dy * dy + dz * dz;

      const Vec ids = Vec::load(moleculeIds + j);
      const Vec otherMolecule = notEqual(ids, molB) * notEqual(ids, paddingId);

      const Vec rrSafe = max(rr, minimumRR);
      const Vec invRR = reciprocal(rrSafe);

      // Lennard-Jones: 4 eps [ (sigma^2/r^2)^6 - (sigma^2/r^2)^3 ] - shift, the parameters of this trial type
      // against every record laid out along the records
      const Vec maskVDW = otherMolecule * less(rr, cutOffVDWSquared);
      const Vec e4 = Vec::load(rowEpsilon4 + j);
      const Vec s6 = Vec::load(rowSigma6 + j);
      const Vec sh = Vec::load(rowShift + j);
      const Vec rri3 = s6 * (invRR * invRR * invRR);
      const Vec uVDW = maskVDW * (e4 * (rri3 * (rri3 - 1.0)) - sh);
      eVDW += uVDW;
      maximumVDW = max(maximumVDW, uVDW);

      if (chargedB)
      {
        const Vec qq = Vec::load(charges + j) * chargeB;
        const Vec maskCharge = otherMolecule * less(rr, cutOffChargeSquared) * nonZero(qq);

        // tabulated erfc(alpha r) / r in r^2: the cubic of the interval, one coefficient record per lane
        // (a lane below the table range reads the first node and is corrected after the loop)
        Vec position = (rrSafe - tableRRMin) * tableInverseSpacing;
        position = min(max(position, 0.0), tableMaximumPosition);
        const Vec node = floor(position);
        const Vec t = position - node;
        std::uint32_t index[lanes];
        node.truncate(index);
        Vec a, b, c, d;
        Vec::gather4(tableCoefficients, index, a, b, c, d);
        const Vec phi = ((a * t + b) * t + c) * t + d;
        anyBelow = max(anyBelow, maskCharge * less(rr, tableRRMin));
        eCharge += maskCharge * (coulombFactor * qq) * phi;
      }
    }

    if (checkOverlap && maximumVDW.maximum() > overlapCriteria) return std::nullopt;
    energyVDW += eVDW.sum();
    energyCharge += eCharge.sum();

    if (anyBelow.maximum() != 0.0) [[unlikely]]
    {
      // interacting pairs closer than the table range (within 1 Angstrom): replace the value of the first node
      // the lanes read by the exact term
      for (std::size_t j = 0; j < padded; ++j)
      {
        if (moleculeIds[j] == molB || moleculeIds[j] == paddingId) continue;
        double3 d(xs[j] - xB, ys[j] - yB, zs[j] - zB);
        if constexpr (Mode != Images::Resolved) d = simulationBox.applyPeriodicBoundaryConditions(d);
        const double rr = double3::dot(d, d);
        const double qq = charges[j] * chargeB;
        if (!(rr < tableRRMin) || !(rr < cutOffChargeSquared) || qq == 0.0) continue;
        double exactPhi, exactDerivative;
        table.exact(std::max(rr, minimumRR), exactPhi, exactDerivative);
        energyCharge += coulombFactor * qq * (exactPhi - tableCoefficients[3]);
      }
    }
  }

  RunningEnergy energySum{};
  bookPairEnergies(neighbourhood, energySum, sign * energyVDW, sign * energyCharge);
  return energySum;
}

/// The conditions of 'Interactions::neighbourhoodFastKernelApplies' for a background whose atoms are all fully
/// coupled when 'backgroundFullyCoupled'.
[[nodiscard]] bool fastKernelApplies(const ForceField& forceField, bool backgroundFullyCoupled,
                                     std::span<const Atom> atoms) noexcept
{
  if (!backgroundFullyCoupled) return false;
  const std::size_t types = forceField.numberOfPseudoAtoms;
  if (types == 0 || forceField.data.size() < types * types) return false;

  if (forceField.useCharge)
  {
    // Ewald through the real-space table of the force field (fully coupled pairs read it in the scalar kernel
    // too); the other charge methods and a table that does not span the cut-off go through the generic kernel
    const EwaldRealSpaceTable& table = forceField.ewaldRealSpaceTable;
    if (forceField.chargeMethod != ForceField::ChargeMethod::Ewald || !forceField.useEwaldRealSpaceTable ||
        table.alpha != forceField.EwaldAlpha || !table.spans(forceField.cutOffCoulomb))
    {
      return false;
    }
  }

  for (const Atom& atom : atoms)
  {
    if (atom.scalingVDW != 1.0 || atom.scalingCoulomb != 1.0 || atom.groupId != 0) return false;
    const std::size_t typeB = static_cast<std::size_t>(atom.type);
    if (typeB >= types) return false;
    // plain (possibly shifted) Lennard-Jones or no interaction between this type and every type
    const VDWParameters* row = forceField.data.data() + typeB * types;
    for (std::size_t a = 0; a < types; ++a)
    {
      if (row[a].type != VDWParameters::Type::LennardJones && row[a].type != VDWParameters::Type::None) return false;
    }
  }
  return true;
}
}  // namespace

bool Interactions::neighbourhoodFastKernelApplies(const ForceField& forceField, const Neighbourhood& neighbourhood,
                                                  std::span<const Atom> atoms) noexcept
{
  return fastKernelApplies(forceField, neighbourhood.fullyCoupled, atoms);
}

std::optional<RunningEnergy> Interactions::computeNeighbourhoodEnergy(
    const ForceField& forceField, const SimulationBox& simulationBox, const Neighbourhood& neighbourhood,
    std::span<const Atom> atoms, double cutOffVDW, double cutOffCoulomb, double sign, bool checkOverlap) noexcept
{
  if (neighbourhoodFastKernelApplies(forceField, neighbourhood, atoms))
  {
    if (neighbourhood.imagesResolved)
    {
      return laneNeighbourhoodEnergy<Images::Resolved>(forceField, simulationBox, neighbourhood, atoms, cutOffVDW,
                                                       cutOffCoulomb, sign, checkOverlap);
    }
    if (simulationBox.type == SimulationBox::Type::Rectangular)
    {
      return laneNeighbourhoodEnergy<Images::Rectangular>(forceField, simulationBox, neighbourhood, atoms, cutOffVDW,
                                                          cutOffCoulomb, sign, checkOverlap);
    }
    return laneNeighbourhoodEnergy<Images::General>(forceField, simulationBox, neighbourhood, atoms, cutOffVDW,
                                                    cutOffCoulomb, sign, checkOverlap);
  }
  if (neighbourhood.imagesResolved)
  {
    return genericNeighbourhoodEnergy<true>(forceField, simulationBox, neighbourhood, atoms, cutOffVDW, cutOffCoulomb,
                                            sign, checkOverlap);
  }
  return genericNeighbourhoodEnergy<false>(forceField, simulationBox, neighbourhood, atoms, cutOffVDW, cutOffCoulomb,
                                           sign, checkOverlap);
}

// ---------------------------------------------------------------------------------------------------
// The single-atom kernel straight out of the cells.
// ---------------------------------------------------------------------------------------------------

namespace
{
/// Per-thread scratch of the single-atom kernel: the records within the cut-off (their squared distance and
/// slot), and -- for the lane pass -- their charge and the Lennard-Jones parameters of the trial type against
/// them, padded to whole lane blocks.
struct AtomKernelScratch
{
  std::vector<double> bucketRR{};  ///< The squared distances of one bucket's records.
  std::vector<double> rr{};
  std::vector<std::uint32_t> slot{};
  std::vector<double> charge{}, epsilon4{}, sigma6{}, shift{};
  std::vector<double> typeEpsilon4{}, typeSigma6{}, typeShift{};  ///< By type, for the trial type.

  /// Room for every atom of the list (an upper bound on the records kept) plus a block of padding.
  void reserve(std::size_t numberOfAtoms, std::size_t capacity)
  {
    if (bucketRR.size() < capacity) bucketRR.resize(capacity);
    const std::size_t size = numberOfAtoms + lanes;
    if (rr.size() >= size) return;
    rr.resize(size);
    slot.resize(size);
    charge.resize(size);
    epsilon4.resize(size);
    sigma6.resize(size);
    shift.resize(size);
  }

  void typeTable(const ForceField& forceField, std::size_t typeB)
  {
    const std::size_t types = forceField.numberOfPseudoAtoms;
    typeEpsilon4.resize(types);
    typeSigma6.resize(types);
    typeShift.resize(types);
    const VDWParameters* parameters = forceField.data.data() + typeB * types;
    for (std::size_t a = 0; a < types; ++a)
    {
      const VDWParameters& pair = parameters[a];
      if (pair.type == VDWParameters::Type::LennardJones)
      {
        const double sigma2 = pair.parameters.y * pair.parameters.y;
        typeEpsilon4[a] = 4.0 * pair.parameters.x;
        typeSigma6[a] = sigma2 * sigma2 * sigma2;
        typeShift[a] = pair.shift;
      }
      else
      {
        typeEpsilon4[a] = 0.0;
        typeSigma6[a] = 1.0;
        typeShift[a] = 0.0;
      }
    }
  }
};

AtomKernelScratch& atomKernelScratch()
{
  thread_local AtomKernelScratch scratch{};
  return scratch;
}

/// The parallel arrays the kept records are read from (a cell list or an all-atom neighbourhood).
struct RecordArrays
{
  const std::uint32_t* type;
  const double* charge;
  const double* scalingVDW;
  const double* scalingCoulomb;
  const std::uint8_t* groupId;
};

/// The second pass of the single-atom kernels: the energy of 'atom' against the 'm' records the first pass kept
/// in the scratch (their squared distances 'rr' and slots), in SIMD lanes when 'fastKernel' (the conditions of
/// 'fastKernelApplies') and through the generic pair arithmetic otherwise; booked as framework or molecule terms.
std::optional<RunningEnergy> evaluateKeptRecords(const ForceField& forceField, const RecordArrays& records,
                                                 const Atom& atom, std::size_t m, double cutOffVDW,
                                                 double cutOffCoulomb, double sign, bool checkOverlap,
                                                 bool fastKernel, bool framework, AtomKernelScratch& scratch) noexcept
{
  RunningEnergy energySum{};
  double energyVDW = 0.0;
  double energyCharge = 0.0;
  const std::size_t typeB = static_cast<std::size_t>(atom.type);
  const double chargeB = atom.charge;
  const std::uint32_t* __restrict slotType = records.type;
  const double* __restrict slotCharge = records.charge;
  const bool useCharge = forceField.useCharge;
  const double cutOffVDWSquared = cutOffVDW * cutOffVDW;
  const double cutOffChargeSquared = cutOffCoulomb * cutOffCoulomb;
  const double overlapCriteria = forceField.energyOverlapCriteria;
  double* __restrict rrs = scratch.rr.data();
  const std::uint32_t* __restrict slots = scratch.slot.data();

  if (!fastKernel)
  {
    // the generic pair arithmetic (soft-core scaling, dU/dlambda), record by record
    const double scalingVDWB = atom.scalingVDW;
    const double scalingCoulombB = atom.scalingCoulomb;
    const std::uint8_t groupIdB = atom.groupId;
    for (std::size_t i = 0; i < m; ++i)
    {
      const std::size_t slot = slots[i];
      const double rr = rrs[i];
      if (rr < cutOffVDWSquared)
      {
        const double scalingA = records.scalingVDW[slot];
        const Potentials::PairDerivatives<0> factors = Potentials::potentialVDW<0>(
            forceField, scalingA, scalingVDWB, rr, static_cast<std::size_t>(slotType[slot]), typeB);
        if (checkOverlap && factors.energy > overlapCriteria) return std::nullopt;
        energyVDW += factors.energy;
        energySum.addDudlambdaVDW(records.groupId[slot], groupIdB, scalingA, scalingVDWB,
                                  sign * factors.dUdlambda);
      }
      if (useCharge && rr < cutOffChargeSquared)
      {
        const double chargeA = slotCharge[slot];
        const double scalingA = records.scalingCoulomb[slot];
        if (chargeA * chargeB != 0.0 && scalingA * scalingCoulombB != 0.0)
        {
          const double r = std::sqrt(rr);
          const Potentials::PairDerivatives<0> factors =
              Potentials::potentialCoulomb<0>(forceField, scalingA, scalingCoulombB, r, chargeA, chargeB);
          energyCharge += factors.energy;
          energySum.addDudlambdaCharge(records.groupId[slot], groupIdB, scalingA, scalingCoulombB,
                                       sign * factors.dUdlambda);
        }
      }
    }
  }
  else
  {
    // Second pass: the charges and the Lennard-Jones parameters (by type) of the records kept, the block padded
    // with far, inert records; then the lane kernel over them.
    scratch.typeTable(forceField, typeB);
    const double* __restrict typeEpsilon4 = scratch.typeEpsilon4.data();
    const double* __restrict typeSigma6 = scratch.typeSigma6.data();
    const double* __restrict typeShift = scratch.typeShift.data();
    double* __restrict qs = scratch.charge.data();
    double* __restrict e4s = scratch.epsilon4.data();
    double* __restrict s6s = scratch.sigma6.data();
    double* __restrict shs = scratch.shift.data();
    const std::size_t padded = (m + lanes - 1) / lanes * lanes;
    for (std::size_t i = 0; i < m; ++i)
    {
      const std::size_t slot = slots[i];
      const std::uint32_t type = slotType[slot];
      qs[i] = slotCharge[slot];
      e4s[i] = typeEpsilon4[type];
      s6s[i] = typeSigma6[type];
      shs[i] = typeShift[type];
    }
    for (std::size_t i = m; i < padded; ++i)
    {
      rrs[i] = 1.0e12;
      qs[i] = 0.0;
      e4s[i] = 0.0;
      s6s[i] = 1.0;
      shs[i] = 0.0;
    }

    const EwaldRealSpaceTable& table = forceField.ewaldRealSpaceTable;
    const double* __restrict tableCoefficients = table.coefficients.data();
    const double tableInverseSpacing = table.inverseSpacing;
    const double tableRRMin = table.rrMin;
    const double tableMaximumPosition = static_cast<double>(table.value.size() >= 2 ? table.value.size() - 2 : 0);
    const double coulombFactor = Units::CoulombicConversionFactor;
    constexpr double minimumRR = 1.0e-4;
    const bool chargedB = useCharge && chargeB != 0.0;

    Vec eVDW = Vec::zero();
    Vec eCharge = Vec::zero();
    Vec maximumVDW = Vec::zero();
    Vec anyBelow = Vec::zero();
    for (std::size_t j = 0; j < padded; j += lanes)
    {
      const Vec rr = Vec::load(rrs + j);
      const Vec rrSafe = max(rr, minimumRR);
      const Vec invRR = reciprocal(rrSafe);

      const Vec maskVDW = less(rr, cutOffVDWSquared);
      const Vec e4 = Vec::load(e4s + j);
      const Vec s6 = Vec::load(s6s + j);
      const Vec sh = Vec::load(shs + j);
      const Vec rri3 = s6 * (invRR * invRR * invRR);
      const Vec uVDW = maskVDW * (e4 * (rri3 * (rri3 - 1.0)) - sh);
      eVDW += uVDW;
      maximumVDW = max(maximumVDW, uVDW);

      if (chargedB)
      {
        const Vec qq = Vec::load(qs + j) * chargeB;
        const Vec maskCharge = less(rr, cutOffChargeSquared) * nonZero(qq);
        Vec position = (rrSafe - tableRRMin) * tableInverseSpacing;
        position = min(max(position, 0.0), tableMaximumPosition);
        const Vec node = floor(position);
        const Vec t = position - node;
        std::uint32_t index[lanes];
        node.truncate(index);
        Vec a, b, c, d;
        Vec::gather4(tableCoefficients, index, a, b, c, d);
        const Vec phi = ((a * t + b) * t + c) * t + d;
        anyBelow = max(anyBelow, maskCharge * less(rr, tableRRMin));
        eCharge += maskCharge * (coulombFactor * qq) * phi;
      }
    }
    if (checkOverlap && maximumVDW.maximum() > overlapCriteria) return std::nullopt;
    energyVDW = eVDW.sum();
    energyCharge = eCharge.sum();

    if (anyBelow.maximum() != 0.0) [[unlikely]]
    {
      // interacting pairs closer than the table range: the first node the lanes read replaced by the exact term
      for (std::size_t i = 0; i < m; ++i)
      {
        const double rr = rrs[i];
        const double qq = qs[i] * chargeB;
        if (!(rr < tableRRMin) || !(rr < cutOffChargeSquared) || qq == 0.0) continue;
        double exactPhi, exactDerivative;
        table.exact(std::max(rr, minimumRR), exactPhi, exactDerivative);
        energyCharge += coulombFactor * qq * (exactPhi - tableCoefficients[3]);
      }
    }
  }

  if (framework)
  {
    energySum.frameworkMoleculeVDW = sign * energyVDW;
    energySum.frameworkMoleculeCharge = sign * energyCharge;
  }
  else
  {
    energySum.moleculeMoleculeVDW = sign * energyVDW;
    energySum.moleculeMoleculeCharge = sign * energyCharge;
  }
  return energySum;
}

/// The single-atom kernel; with 'PairOnce' only the records with an atom index above 'ownIndex' are kept, so a
/// sum over every atom of the list visits each pair once (the full-energy sum).
template <bool PairOnce>
std::optional<RunningEnergy> cellListAtomEnergy(const ForceField& forceField, const SimulationBox& simulationBox,
                                                const MCCellList& cellList, const Atom& atom, double cutOffVDW,
                                                double cutOffCoulomb, double sign, bool checkOverlap,
                                                std::optional<std::uint32_t> skipMolecule,
                                                Interactions::Neighbourhood::Background background,
                                                [[maybe_unused]] std::uint32_t ownIndex) noexcept
{
  using Interactions::Neighbourhood;
  const bool useCharge = forceField.useCharge;
  const double cutOff = std::max(cutOffVDW, useCharge ? cutOffCoulomb : 0.0);
  const double cutOffSquared = cutOff * cutOff;
  const bool framework = background == Neighbourhood::Background::Framework;
  const double3 centre = atom.position;

  // the exclusions of the gather and the kernel: the atoms of the skipped molecule, and (for a molecule
  // background) those of the trial atom's own molecule
  const std::uint32_t skipped = skipMolecule.value_or(paddingMoleculeId);
  const std::uint32_t molB = framework ? paddingMoleculeId : atom.moleculeId;

  // the image placement of the gather ('Neighbourhood::gatherFromCells'): the cell's shift along the axes the
  // stencil does not wrap, the minimum image to the atom along the whole axes; a cut-off reaching half the box
  // leaves the minimum image to the pair (as the kernels do for a gather of that radius)
  const std::array<bool, 3> whole = cellList.wholeAxes(cutOff);
  const bool anyWhole = whole[0] || whole[1] || whole[2];
  const bool perAxisImage = anyWhole && cellList.orthorhombic;
  const double3& widths = cellList.perpendicularWidths;
  const double wrapX = whole[0] ? widths.x : 0.0, wrapY = whole[1] ? widths.y : 0.0, wrapZ = whole[2] ? widths.z : 0.0;
  const double unwrapX = whole[0] ? 1.0 / widths.x : 0.0, unwrapY = whole[1] ? 1.0 / widths.y : 0.0,
               unwrapZ = whole[2] ? 1.0 / widths.z : 0.0;
  const double halfWidth = 0.5 * std::min({widths.x, widths.y, widths.z});
  const bool minimumImagePerPair = !(cutOff < halfWidth);

  AtomKernelScratch& scratch = atomKernelScratch();
  scratch.reserve(cellList.numberOfAtoms, cellList.capacity);
  double* __restrict bucketRR = scratch.bucketRR.data();
  double* __restrict rrs = scratch.rr.data();
  std::uint32_t* __restrict slots = scratch.slot.data();

  const double* __restrict X = cellList.slotX.data();
  const double* __restrict Y = cellList.slotY.data();
  const double* __restrict Z = cellList.slotZ.data();
  const std::uint32_t* __restrict slotId = cellList.slotMoleculeId.data();
  [[maybe_unused]] const std::uint32_t* __restrict slotIndex = cellList.slotAtomIndex.data();

  // First pass, over the cells within the cut-off: the squared distance of every record of the stencil, the
  // records within the cut-off (the excluded molecules left out) compacted by index.
  std::size_t m = 0;
  cellList.forEachCellWithin(centre, cutOff,
                             [&](std::size_t cell, const double3& shift)
                             {
                               const std::size_t base = cell * cellList.capacity;
                               const std::size_t n = cellList.counts[cell];
                               const double ox = shift.x - centre.x, oy = shift.y - centre.y, oz = shift.z - centre.z;
                               const double* __restrict cx = X + base;
                               const double* __restrict cy = Y + base;
                               const double* __restrict cz = Z + base;
                               const std::uint32_t* __restrict ids = slotId + base;
                               [[maybe_unused]] const std::uint32_t* __restrict indices = slotIndex + base;

                               // the squared distances of the bucket, streamed without a dependency between
                               // the records (vectorizable), then the compaction of the records kept
                               if (!anyWhole && !minimumImagePerPair)
                               {
                                 for (std::size_t k = 0; k < n; ++k)
                                 {
                                   const double dx = cx[k] + ox, dy = cy[k] + oy, dz = cz[k] + oz;
                                   bucketRR[k] = dx * dx + dy * dy + dz * dz;
                                 }
                               }
                               else if (perAxisImage && !minimumImagePerPair)
                               {
                                 for (std::size_t k = 0; k < n; ++k)
                                 {
                                   double dx = cx[k] + ox, dy = cy[k] + oy, dz = cz[k] + oz;
                                   dx -= wrapX * std::round(dx * unwrapX);
                                   dy -= wrapY * std::round(dy * unwrapY);
                                   dz -= wrapZ * std::round(dz * unwrapZ);
                                   bucketRR[k] = dx * dx + dy * dy + dz * dz;
                                 }
                               }
                               else
                               {
                                 for (std::size_t k = 0; k < n; ++k)
                                 {
                                   const double3 d = simulationBox.applyPeriodicBoundaryConditions(
                                       double3(cx[k] + ox, cy[k] + oy, cz[k] + oz));
                                   bucketRR[k] = double3::dot(d, d);
                                 }
                               }
                               for (std::size_t k = 0; k < n; ++k)
                               {
                                 const double rr = bucketRR[k];
                                 bool kept = (rr <= cutOffSquared) & (ids[k] != molB) & (ids[k] != skipped);
                                 if constexpr (PairOnce) kept = kept & (indices[k] > ownIndex);
                                 rrs[m] = rr;
                                 slots[m] = static_cast<std::uint32_t>(base + k);
                                 m += kept;
                               }
                             });

  const RecordArrays records{cellList.slotType.data(), cellList.slotCharge.data(), cellList.slotScalingVDW.data(),
                             cellList.slotScalingCoulomb.data(), cellList.slotGroupId.data()};
  const bool fastKernel = fastKernelApplies(forceField, cellList.numberOfPartiallyCoupledAtoms == 0, {&atom, 1});
  return evaluateKeptRecords(forceField, records, atom, m, cutOffVDW, cutOffCoulomb, sign, checkOverlap, fastKernel,
                             framework, scratch);
}
}  // namespace

std::optional<RunningEnergy> Interactions::computeCellListAtomEnergy(
    const ForceField& forceField, const SimulationBox& simulationBox, const MCCellList& cellList, const Atom& atom,
    double cutOffVDW, double cutOffCoulomb, double sign, bool checkOverlap, std::optional<std::uint32_t> skipMolecule,
    Neighbourhood::Background background) noexcept
{
  return cellListAtomEnergy<false>(forceField, simulationBox, cellList, atom, cutOffVDW, cutOffCoulomb, sign,
                                   checkOverlap, skipMolecule, background, 0);
}

namespace
{
/**
 * \brief The sum over the atoms of a per-atom energy, in chunks of consecutive atoms: on the thread pool (the
 * helper threads plus the calling thread pull the next chunk from a shared counter) when the pool is active and
 * the system is large enough for the dispatch to pay off, on the calling thread otherwise. The chunk sums are
 * added in chunk order, so the result does not depend on the scheduling.
 *
 * 'perAtom(first, last, sum)' adds the energies of the atoms [first, last) to 'sum'; it runs on several threads
 * at once and must only use per-thread scratch ('atomKernelScratch') and read-only shared data.
 */
template <typename PerAtom>
RunningEnergy chunkedAtomSum(std::size_t numberOfAtoms, PerAtom&& perAtom) noexcept
{
  constexpr std::size_t chunk = 32;
  constexpr std::size_t minimumAtomsForThreads = 512;
  const std::size_t numberOfChunks = (numberOfAtoms + chunk - 1) / chunk;

  auto& pool = ThreadPool::ThreadPool<ThreadPool::details::default_function_type, std::jthread>::instance();
  const std::size_t helpers =
      pool.getThreadingType() == ThreadPool::ThreadingType::ThreadPool ? pool.getThreadCount() : 0;
  if (helpers == 0 || numberOfAtoms < minimumAtomsForThreads)
  {
    RunningEnergy sum{};
    perAtom(0, numberOfAtoms, sum);
    return sum;
  }

  std::vector<RunningEnergy> partial(numberOfChunks);
  std::atomic<std::size_t> next{0};
  auto drain = [&]()
  {
    for (;;)
    {
      const std::size_t c = next.fetch_add(1, std::memory_order_relaxed);
      if (c >= numberOfChunks) return;
      perAtom(c * chunk, std::min(numberOfAtoms, (c + 1) * chunk), partial[c]);
    }
  };
  std::vector<std::future<void>> pending(helpers);
  for (std::size_t i = 0; i != helpers; ++i) pending[i] = pool.enqueue(drain);
  drain();
  for (std::future<void>& task : pending) task.get();

  RunningEnergy sum{};
  for (const RunningEnergy& part : partial) sum += part;
  return sum;
}
}  // namespace

RunningEnergy Interactions::computeCellListInterMolecularEnergy(const ForceField& forceField,
                                                                const SimulationBox& simulationBox,
                                                                const MCCellList& cellList,
                                                                std::span<const Atom> moleculeAtoms) noexcept
{
  const double cutOffVDW = forceField.cutOffMoleculeVDW;
  const double cutOffCoulomb = forceField.cutOffCoulomb;
  return chunkedAtomSum(moleculeAtoms.size(),
                        [&](std::size_t first, std::size_t last, RunningEnergy& sum)
                        {
                          for (std::size_t i = first; i < last; ++i)
                          {
                            // no overlap check: the kernel then always returns an energy
                            sum += *cellListAtomEnergy<true>(forceField, simulationBox, cellList, moleculeAtoms[i],
                                                             cutOffVDW, cutOffCoulomb, 1.0, false, std::nullopt,
                                                             Neighbourhood::Background::Molecules,
                                                             static_cast<std::uint32_t>(i));
                          }
                        });
}

RunningEnergy Interactions::computeAllPairsInterMolecularEnergy(const ForceField& forceField,
                                                                const SimulationBox& simulationBox,
                                                                std::span<const Atom> moleculeAtoms) noexcept
{
  const bool useCharge = forceField.useCharge;
  const double cutOffVDW = forceField.cutOffMoleculeVDW;
  const double cutOffCoulomb = forceField.cutOffCoulomb;
  const double cutOff = std::max(cutOffVDW, useCharge ? cutOffCoulomb : 0.0);
  const double cutOffSquared = cutOff * cutOff;
  const std::size_t numberOfAtoms = moleculeAtoms.size();

  // the structure-of-arrays copy of the atoms, in atom order (the padding records are never visited); gathered
  // once here and read by every thread of the sum
  thread_local Neighbourhood neighbourhood{};
  neighbourhood.gatherAll(moleculeAtoms, cutOff);
  const double* __restrict xs = neighbourhood.x.data();
  const double* __restrict ys = neighbourhood.y.data();
  const double* __restrict zs = neighbourhood.z.data();
  const RecordArrays records{neighbourhood.type.data(), neighbourhood.charge.data(), neighbourhood.scalingVDW.data(),
                             neighbourhood.scalingCoulomb.data(), neighbourhood.groupId.data()};
  const bool fastKernel = fastKernelApplies(forceField, neighbourhood.fullyCoupled, moleculeAtoms);

  const bool rectangular = simulationBox.type == SimulationBox::Type::Rectangular;
  const double3x3& cell = simulationBox.cell;
  const double3x3& inverseCell = simulationBox.inverseCell;
  const double lengthX = cell.ax, lengthY = cell.by, lengthZ = cell.cz;
  const double inverseLengthX = inverseCell.ax, inverseLengthY = inverseCell.by, inverseLengthZ = inverseCell.cz;

  // every atom against the records behind its molecule (the molecule's atoms are contiguous): each pair once,
  // the pairs within the molecule never visited
  auto perAtom = [&](std::size_t first, std::size_t last, RunningEnergy& sum)
  {
    AtomKernelScratch& scratch = atomKernelScratch();
    scratch.reserve(numberOfAtoms, numberOfAtoms);
    double* __restrict allRR = scratch.bucketRR.data();
    double* __restrict rrs = scratch.rr.data();
    std::uint32_t* __restrict slots = scratch.slot.data();

    for (std::size_t i = first; i < last; ++i)
    {
      const Atom& atom = moleculeAtoms[i];
      const double xB = atom.position.x, yB = atom.position.y, zB = atom.position.z;
      // the first record behind the atom's molecule
      std::size_t end = i + 1;
      while (end < numberOfAtoms && moleculeAtoms[end].moleculeId == atom.moleculeId) ++end;

      // First pass: the minimum-image squared distances to the records behind the molecule, streamed without a
      // dependency between the records, then the compaction of those within the cut-off.
      if (rectangular)
      {
        for (std::size_t j = end; j < numberOfAtoms; ++j)
        {
          double dx = xs[j] - xB, dy = ys[j] - yB, dz = zs[j] - zB;
          dx -= lengthX * std::round(dx * inverseLengthX);
          dy -= lengthY * std::round(dy * inverseLengthY);
          dz -= lengthZ * std::round(dz * inverseLengthZ);
          allRR[j] = dx * dx + dy * dy + dz * dz;
        }
      }
      else
      {
        for (std::size_t j = end; j < numberOfAtoms; ++j)
        {
          const double dx = xs[j] - xB, dy = ys[j] - yB, dz = zs[j] - zB;
          double sx = inverseCell.m11 * dx + inverseCell.m12 * dy + inverseCell.m13 * dz;
          double sy = inverseCell.m21 * dx + inverseCell.m22 * dy + inverseCell.m23 * dz;
          double sz = inverseCell.m31 * dx + inverseCell.m32 * dy + inverseCell.m33 * dz;
          sx -= std::round(sx);
          sy -= std::round(sy);
          sz -= std::round(sz);
          const double ex = cell.m11 * sx + cell.m12 * sy + cell.m13 * sz;
          const double ey = cell.m21 * sx + cell.m22 * sy + cell.m23 * sz;
          const double ez = cell.m31 * sx + cell.m32 * sy + cell.m33 * sz;
          allRR[j] = ex * ex + ey * ey + ez * ez;
        }
      }
      std::size_t m = 0;
      for (std::size_t j = end; j < numberOfAtoms; ++j)
      {
        const double rr = allRR[j];
        rrs[m] = rr;
        slots[m] = static_cast<std::uint32_t>(j);
        m += (rr <= cutOffSquared);
      }

      sum += *evaluateKeptRecords(forceField, records, atom, m, cutOffVDW, cutOffCoulomb, 1.0, false, fastKernel,
                                  false, scratch);
    }
  };
  return chunkedAtomSum(numberOfAtoms, perAtom);
}
