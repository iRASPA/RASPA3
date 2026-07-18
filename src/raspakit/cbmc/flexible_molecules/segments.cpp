module;

module cbmc_segments;

import std;

import atom;
import component;
import connectivity_table;
import intra_molecular_potentials;

// Determine the next group-aware growth step given the set of already-placed beads.
//
// Mirrors 'ConnectivityTable::nextBeads' for the flexible case, but when the frontier bond leads
// into a rigid group the whole (still unplaced) group is emitted as a single rigid-unit segment.
static CBMC::GrowSegment nextGroupAwareStep(const Component &component, const std::vector<std::size_t> &placedBeads)
{
  const ConnectivityTable &connectivity = component.connectivityTable;
  std::size_t numberOfBeads = connectivity.numberOfBeads;

  std::vector<bool> placed(numberOfBeads, false);
  for (std::size_t bead : placedBeads) placed[bead] = true;

  // 'filteredInteractions' takes mutable spans, so keep a mutable copy of the placed beads.
  std::vector<std::size_t> placedVec(placedBeads.begin(), placedBeads.end());

  // Search for the first frontier bond (a placed bead 'k' connected to an unplaced bead 'j'),
  // scanning placed beads in order and neighbors by ascending index (identical ordering to
  // 'ConnectivityTable::nextBeads' so the flexible case is unchanged).
  std::optional<std::size_t> currentBeadOpt{};
  std::optional<std::size_t> firstNextBead{};
  for (std::size_t k : placedBeads)
  {
    for (std::size_t j = 0; j != numberOfBeads; ++j)
    {
      if (connectivity[j, k] && !placed[j])
      {
        currentBeadOpt = k;
        firstNextBead = j;
        break;
      }
    }
    if (currentBeadOpt.has_value()) break;
  }

  if (!currentBeadOpt.has_value())
  {
    throw std::runtime_error(std::format("Error in CBMC: No bead can be grown\n"));
  }

  std::size_t currentBead = currentBeadOpt.value();
  std::size_t nextBead = firstNextBead.value();

  std::optional<std::size_t> rigidGroup = component.rigidGroupContaining(nextBead);

  if (rigidGroup.has_value())
  {
    bool anchorInGroup = component.rigidGroupContaining(currentBead) == rigidGroup;

    if (!anchorInGroup)
    {
      // The frontier bond leads from a placed bead outside the group into the group. Do not place
      // the whole group at once: first grow the single connecting atom as an ordinary flexible bond
      // (so its junction bond, and its bend/torsion when a previous bead exists, are sampled). Once
      // this atom is fixed the remainder of the rigid group is grown as a rigid body hinged on it.
      std::vector<std::size_t> singleNextBead{nextBead};

      std::optional<std::size_t> previousBead{};
      for (std::size_t i = 0; i != numberOfBeads; ++i)
      {
        if (connectivity[i, currentBead] && placed[i])
        {
          previousBead = i;
          break;
        }
      }

      Potentials::IntraMolecularPotentials intra =
          component.intraMolecularPotentials.filteredInteractions(numberOfBeads, placedVec, singleNextBead);

      return {previousBead, currentBead, singleNextBead, false, std::nullopt, intra};
    }

    // The anchor ('currentBead') is a placed atom of the group: grow the remaining unplaced atoms of
    // the group as one rigid body hinged on the anchor.
    const MoleculeGroup &group = component.groups[rigidGroup.value()];
    std::vector<std::size_t> nextBeads{};
    nextBeads.reserve(group.atoms.size());
    for (std::size_t atom : group.atoms)
    {
      if (!placed[atom]) nextBeads.push_back(atom);
    }

    // Reference bead for the bend/torsion bias: a placed neighbor of 'currentBead' outside the group
    // (deterministically the first such neighbor). Absent only when the group is the growth seed.
    std::optional<std::size_t> previousBead{};
    for (std::size_t i = 0; i != numberOfBeads; ++i)
    {
      if (connectivity[i, currentBead] && placed[i] && component.atomGroupIds[i] != rigidGroup.value())
      {
        previousBead = i;
        break;
      }
    }

    Potentials::IntraMolecularPotentials intra =
        component.intraMolecularPotentials.filteredInteractions(numberOfBeads, placedVec, nextBeads);

    return {previousBead, currentBead, nextBeads, true, rigidGroup, intra};
  }

  // Flexible segment. Reproduce 'ConnectivityTable::nextBeads': determine the previous bead and,
  // when there is a previous bead, grow all unplaced (flexible) neighbors of 'currentBead' together
  // (branch point); when there is no previous bead grow only the single frontier bead.
  std::vector<std::size_t> nextBeads{};
  std::size_t numberOfPreviousBeads{};
  std::optional<std::size_t> previousBead{};
  for (std::size_t i = 0; i != numberOfBeads; ++i)
  {
    if (!connectivity[i, currentBead]) continue;

    if (!placed[i])
    {
      // Rigid-group neighbors are grown as their own segment, never mixed into a flexible branch.
      if (!component.rigidGroupContaining(i).has_value())
      {
        nextBeads.push_back(i);
      }
    }
    else
    {
      ++numberOfPreviousBeads;
      // Keep the first (lowest-index) placed neighbor as the deterministic bend/torsion reference.
      if (!previousBead.has_value()) previousBead = i;
    }
  }

  if (numberOfPreviousBeads == 0)
  {
    std::vector<std::size_t> singleNextBead{nextBead};
    Potentials::IntraMolecularPotentials intra =
        component.intraMolecularPotentials.filteredInteractions(numberOfBeads, placedVec, singleNextBead);
    return {std::nullopt, currentBead, singleNextBead, false, std::nullopt, intra};
  }

  // More than one placed neighbor is only allowed when 'currentBead' is part of an already-placed
  // rigid group (e.g. a flexible tail growing off a ring atom): pick the first placed neighbor as
  // the reference. For a purely flexible anchor this signals an unsupported flexible ring.
  if (numberOfPreviousBeads > 1 && !component.rigidGroupContaining(currentBead).has_value())
  {
    throw std::runtime_error(std::format("Error in CBMC: Multiple previous beads\n"));
  }

  Potentials::IntraMolecularPotentials intra =
      component.intraMolecularPotentials.filteredInteractions(numberOfBeads, placedVec, nextBeads);
  return {previousBead, currentBead, nextBeads, false, std::nullopt, intra};
}

std::vector<CBMC::GrowSegment> CBMC::buildGrowSegments(const Component &component,
                                                       const std::vector<std::size_t> &beadsAlreadyPlaced)
{
  std::size_t numberOfBeads = component.connectivityTable.numberOfBeads;

  std::vector<GrowSegment> segments{};
  std::vector<std::size_t> placed(beadsAlreadyPlaced.begin(), beadsAlreadyPlaced.end());

  while (placed.size() < numberOfBeads)
  {
    GrowSegment segment = nextGroupAwareStep(component, placed);
    placed.insert(placed.end(), segment.nextBeads.begin(), segment.nextBeads.end());
    segments.push_back(std::move(segment));
  }

  return segments;
}
