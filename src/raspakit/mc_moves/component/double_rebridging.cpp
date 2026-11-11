module;

module mc_moves_double_rebridging;

import std;

import component;
import atom;
import molecule;
import double3;
import double3x3;
import simulationbox;
import randomnumbers;
import system;
import running_energy;
import interactions_framework_molecule;
import interactions_intermolecular;
import interactions_ewald;
import interactions_external_field;
import mc_moves_move_types;
import mc_moves_cputime;
import mc_moves_concerted_rotation_geometry;
import mc_moves_bridging_common;

std::optional<RunningEnergy> MC_Moves::intramolecularDoubleRebridgingMove(RandomNumber& random, System& system,
                                                                          std::size_t selectedComponent,
                                                                          std::size_t selectedMolecule)
{
  Move::Types move = Move::Types::IntramolecularDoubleRebridging;
  Component& component = system.components[selectedComponent];

  // Restrictions: the polarization energy update is not implemented; molecules without a valid site
  // pair cannot rebridge. Neither changes during a simulation, so the attempts are not counted as
  // trials.
  const Component::BridgingTopology& topology = component.bridgingTopology();
  if (topology.sitePairs.empty() || system.forceField.computePolarization) return std::nullopt;

  // A cross-linked molecule can not rebridge: the move relabels the repeat units between the two bridged sites,
  // which would silently move a linked reactive site to another atom (the cross-link table refers to atom indices).
  if (system.crossLinks.moleculeIsLinked(selectedComponent, selectedMolecule)) return std::nullopt;

  component.mc_moves_statistics.addTrial(move);

  const std::vector<Component::BridgingTopology::Unit>& units = topology.units;
  const std::size_t pairIndex =
      static_cast<std::size_t>(random.uniform() * static_cast<double>(topology.sitePairs.size()));
  const std::array<std::size_t, 2> pair = topology.sitePairs[pairIndex];
  const bool reprojectSideGroups = topology.sitePairReprojectsSideGroups[pairIndex];
  const std::size_t a = pair[0], b = pair[1];
  const double maximumDistance = topology.maximumBridgeDistance;

  std::span<Atom> molecule_atoms = system.spanOfMolecule(selectedComponent, selectedMolecule);
  Molecule& molecule = system.moleculeData[system.moleculeIndexOfComponent(selectedComponent, selectedMolecule)];
  const std::vector<Atom> oldAtoms(molecule_atoms.begin(), molecule_atoms.end());

  auto position = [&](std::size_t unit) { return oldAtoms[units[unit].backboneAtom].position; };

  // Reach criterion, symmetric between the states: the old trimers span (a, a+4) and (b, b+4), the new
  // ones span (a, b) and (a+4, b+4). The reverse move checks the same four distances.
  if ((position(a) - position(a + 4)).length() > maximumDistance ||
      (position(b) - position(b + 4)).length() > maximumDistance ||
      (position(a) - position(b)).length() > maximumDistance ||
      (position(a + 4) - position(b + 4)).length() > maximumDistance)
  {
    return std::nullopt;
  }

  // Closure 1: the head dimer (a-1, a) bridges onto the far end of the segment, (b, b-1), which
  // becomes the new (a+4, a+5); the third chart atom b-2 exists when the segment has three or more
  // units. Closure 2: the near end of the segment, now (a+5, a+4) in reversed order, bridges onto
  // the tail dimer (b+4, b+5) that stays in place.
  const bool segmentHasThreeUnits = (b >= a + 6);
  const std::optional<double3> tailA8 = (b + 6 < units.size()) ? std::optional<double3>(position(b + 6)) : std::nullopt;

  const Bridging::Anchors oldAnchors1{
      .a1 = position(a - 1),
      .a2 = position(a),
      .a6 = position(a + 4),
      .a7 = position(a + 5),
      .a8 = segmentHasThreeUnits ? std::optional<double3>(position(a + 6)) : std::nullopt};
  const Bridging::Anchors newAnchors1{
      .a1 = position(a - 1),
      .a2 = position(a),
      .a6 = position(b),
      .a7 = position(b - 1),
      .a8 = segmentHasThreeUnits ? std::optional<double3>(position(b - 2)) : std::nullopt};
  const Bridging::Anchors oldAnchors2{
      .a1 = position(b - 1), .a2 = position(b), .a6 = position(b + 4), .a7 = position(b + 5), .a8 = tailA8};
  const Bridging::Anchors newAnchors2{
      .a1 = position(a + 5), .a2 = position(a + 4), .a6 = position(b + 4), .a7 = position(b + 5), .a8 = tailA8};
  const ConcertedRotation::Trimer oldTrimer1{position(a + 1), position(a + 2), position(a + 3)};
  const ConcertedRotation::Trimer oldTrimer2{position(b + 1), position(b + 2), position(b + 3)};

  const std::optional<Bridging::Closure> closure1 = Bridging::rebridge(random, oldAnchors1, oldTrimer1, newAnchors1);
  if (!closure1.has_value()) return std::nullopt;
  const std::optional<Bridging::Closure> closure2 = Bridging::rebridge(random, oldAnchors2, oldTrimer2, newAnchors2);
  if (!closure2.has_value()) return std::nullopt;

  // The two states chart the (rigidly fixed) segment from opposite ends; the densities are relative
  // to the same rigid-motion measure of the segment only after dividing by the chart volume elements.
  const double logChart = std::log(Bridging::chartVolumeElement(newAnchors1.a6, newAnchors1.a7, newAnchors1.a8)) -
                          std::log(Bridging::chartVolumeElement(oldAnchors1.a6, oldAnchors1.a7, oldAnchors1.a8));

  // Trial chain: the backbone of the segment a+4 ... b reversed slot by slot (the units are congruent,
  // so the side atoms map in order), the two trimers re-bridged with their side groups carried along.
  std::vector<Atom> trialAtoms = oldAtoms;
  for (std::size_t k = a + 4; k <= b; ++k)
  {
    const Component::BridgingTopology::Unit& target = units[k];
    const Component::BridgingTopology::Unit& source = units[a + 4 + b - k];
    trialAtoms[target.backboneAtom].position = oldAtoms[source.backboneAtom].position;
  }
  if (!Bridging::placeTrimer(trialAtoms, oldAtoms, units, a, oldAnchors1, oldTrimer1, newAnchors1, closure1->trimer))
  {
    return std::nullopt;
  }
  if (!Bridging::placeTrimer(trialAtoms, oldAtoms, units, b, oldAnchors2, oldTrimer2, newAnchors2, closure2->trimer))
  {
    return std::nullopt;
  }

  // The side groups of the reversed segment. The chain now runs through every unit in the opposite direction
  // (its previous and next backbone neighbours are interchanged). Carrying the side atoms rigidly from the local
  // frame (atom, previous, next) of the source unit to the frame of the target slot, as the trimer units are,
  // keeps the handedness of every unit with respect to the chain direction; copying them in place inverts it
  // (the two in-plane frame vectors are interchanged, which reverses the frame normal). The topology selects
  // the map that takes the declared stereo pattern onto itself ('sitePairReprojectsSideGroups'). Either map is
  // a proper rigid motion per unit whose inverse the reverse move applies (the frames are functions of the
  // backbone positions, which it restores exactly), so no Jacobian enters.
  for (std::size_t k = a + 4; k <= b; ++k)
  {
    const Component::BridgingTopology::Unit& target = units[k];
    if (target.sideAtoms.empty()) continue;
    const std::size_t s = a + 4 + b - k;
    const Component::BridgingTopology::Unit& source = units[s];
    if (!reprojectSideGroups)
    {
      for (std::size_t p = 0; p != target.sideAtoms.size(); ++p)
      {
        trialAtoms[target.sideAtoms[p]].position = oldAtoms[source.sideAtoms[p]].position;
      }
      continue;
    }
    const double3 oldCenter = position(s);
    const std::optional<double3x3> oldFrame =
        ConcertedRotation::localFrame(oldCenter, position(s - 1), position(s + 1));
    const double3 newCenter = trialAtoms[target.backboneAtom].position;
    const std::optional<double3x3> newFrame =
        ConcertedRotation::localFrame(newCenter, trialAtoms[units[k - 1].backboneAtom].position,
                                      trialAtoms[units[k + 1].backboneAtom].position);
    // a collinear unit of the segment has no frame to carry its side group with: refuse the move
    if (!oldFrame.has_value() || !newFrame.has_value()) return std::nullopt;
    for (std::size_t p = 0; p != target.sideAtoms.size(); ++p)
    {
      trialAtoms[target.sideAtoms[p]].position = ConcertedRotation::transformRigidly(
          oldCenter, *oldFrame, newCenter, *newFrame, oldAtoms[source.sideAtoms[p]].position);
    }
  }

  if (system.insideBlockedPockets(component, trialAtoms)) return std::nullopt;

  // a declared chiral centre is a hard constraint: the units of the segment and of the trimers keep their
  // handedness, but a centre at an anchor, whose neighbour a3 or a5 moves while its side group stays, can invert
  if (!component.intraMolecularPotentials.preservesChiralCenters(trialAtoms)) return std::nullopt;

  // The atoms whose positions change are the two trimers with their side groups and, when re-projected, the
  // side groups of the reversed segment; the backbone of the segment (and its side groups when copied in place)
  // is the same set of (type, position) pairs before and after, so it drops out of every non-bonded difference
  // with the surroundings.
  std::vector<Atom> movedOld{}, movedNew{};
  for (std::size_t site : {a, b})
  {
    for (std::size_t m = 1; m <= 3; ++m)
    {
      const Component::BridgingTopology::Unit& unit = units[site + m];
      movedOld.push_back(oldAtoms[unit.backboneAtom]);
      movedNew.push_back(trialAtoms[unit.backboneAtom]);
      for (std::size_t atom : unit.sideAtoms)
      {
        movedOld.push_back(oldAtoms[atom]);
        movedNew.push_back(trialAtoms[atom]);
      }
    }
  }
  for (std::size_t k = a + 4; reprojectSideGroups && k <= b; ++k)
  {
    for (std::size_t atom : units[k].sideAtoms)
    {
      movedOld.push_back(oldAtoms[atom]);
      movedNew.push_back(trialAtoms[atom]);
    }
  }

  std::optional<RunningEnergy> externalFieldMolecule =
      timed(system, component, move, Move::Timing::ExternalFieldMolecule,
            [&]
            {
              return Interactions::computeExternalFieldEnergyDifference(
                  system.hasExternalField, system.forceField, system.simulationBox,
                  system.externalFieldInterpolationGrid, movedNew, movedOld);
            });
  if (!externalFieldMolecule.has_value()) return std::nullopt;

  std::optional<RunningEnergy> frameworkMolecule =
      timed(system, component, move, Move::Timing::FrameworkMolecule,
            [&]
            {
              return Interactions::computeFrameworkMoleculeEnergyDifference(
                  system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
                  system.frameworkCellList(), system.spanOfFrameworkAtoms(), movedNew, movedOld);
            });
  if (!frameworkMolecule.has_value()) return std::nullopt;

  std::optional<RunningEnergy> interMolecule = timed(system, component, move, Move::Timing::MoleculeMolecule,
                                                     [&]
                                                     {
                                                       return Interactions::computeInterMolecularEnergyDifference(
                                                           system.forceField, system.simulationBox, system.cellList(),
                                                           system.spanOfMoleculeAtoms(), movedNew, movedOld);
                                                     });
  if (!interMolecule.has_value()) return std::nullopt;

  // Ewald: the full chain is passed so that the intramolecular exclusion terms of the trimers with the
  // rest of the chain are updated (the reversed segment contributes identical terms on both sides).
  RunningEnergy ewaldFourierEnergy =
      timed(system, component, move, Move::Timing::Ewald,
            [&]
            {
              return Interactions::energyDifferenceEwaldFourier(system.eik_x, system.eik_y, system.eik_z, system.eik_xy,
                                                                system.storedEik, system.trialEik, system.forceField,
                                                                system.simulationBox, system.components, trialAtoms, molecule_atoms);
            });

  // Intramolecular: the bonded terms across the four seams and the intramolecular non-bonded terms
  // change. Recomputing all internal terms keeps the bookkeeping exact.
  RunningEnergy internalDifference = component.intraMolecularPotentials.computeInternalEnergies(system.forceField, system.simulationBox, trialAtoms) -
                                     component.intraMolecularPotentials.computeInternalEnergies(system.forceField, system.simulationBox, molecule_atoms);

  RunningEnergy energyDifference = externalFieldMolecule.value() + frameworkMolecule.value() + interMolecule.value() +
                                   ewaldFourierEnergy + internalDifference;

  component.mc_moves_statistics.addConstructed(move);

  const double logAcceptance =
      -system.beta * energyDifference.potentialEnergy() + closure1->logWeight + closure2->logWeight + logChart;
  if (random.uniform() < std::exp(logAcceptance))
  {
    component.mc_moves_statistics.addAccepted(move);

    Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.trialEik);

    std::copy(trialAtoms.cbegin(), trialAtoms.cend(), molecule_atoms.begin());
    // The atoms of the reversed segment did not move: the chain runs through them the other way, so each unit's
    // atoms now sit in the mirrored slots. Their persistent tags (and dynamics) go with them, as the exchanged
    // tails of the double bridging do; the re-bridged trimers were placed anew and keep theirs.
    const std::size_t first = static_cast<std::size_t>(molecule_atoms.data() - system.atomData.data());
    for (std::size_t k = a + 4; k < a + 4 + b - k; ++k)
    {
      const Component::BridgingTopology::Unit& unit = units[k];
      const Component::BridgingTopology::Unit& mirror = units[a + 4 + b - k];
      system.exchangeAtomTags(first + unit.backboneAtom, first + mirror.backboneAtom);
      for (std::size_t p = 0; p != unit.sideAtoms.size(); ++p)
      {
        system.exchangeAtomTags(first + unit.sideAtoms[p], first + mirror.sideAtoms[p]);
      }
    }
    // in-place change of one molecule: keep the molecule cell list consistent, like the other in-place moves
    system.cellListAtomsMoved(molecule_atoms);
    molecule.centerOfMassPosition = component.computeCenterOfMass(molecule_atoms);

    return energyDifference;
  }
  return std::nullopt;
}
