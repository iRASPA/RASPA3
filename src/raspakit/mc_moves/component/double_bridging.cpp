module;

module mc_moves_double_bridging;

import std;

import component;
import atom;
import molecule;
import double3;
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

std::optional<RunningEnergy> MC_Moves::doubleBridgingMove(RandomNumber& random, System& system,
                                                          std::size_t selectedComponent, std::size_t selectedMolecule)
{
  Move::Types move = Move::Types::DoubleBridging;
  Component& component = system.components[selectedComponent];

  // Restrictions: the polarization energy update is not implemented; molecules without a valid site
  // cannot bridge; only whole molecules take part (the exchanged tails must carry the same scaling),
  // so at least two integer molecules are needed. None of these change with the configuration, so
  // the attempts are not counted as trials.
  const Component::BridgingTopology& topology = component.bridgingTopology();
  const std::size_t numberOfFractional = system.numberOfFractionalMoleculesPerComponent[selectedComponent];
  const std::size_t numberOfMolecules = system.numberOfMoleculesPerComponent[selectedComponent];
  if (topology.sites.empty() || system.forceField.computePolarization || numberOfMolecules < numberOfFractional + 2 ||
      selectedMolecule < numberOfFractional)
  {
    return std::nullopt;
  }

  component.mc_moves_statistics.addTrial(move);

  const std::vector<Component::BridgingTopology::Unit>& units = topology.units;
  const std::size_t site =
      topology.sites[static_cast<std::size_t>(random.uniform() * static_cast<double>(topology.sites.size()))];
  const double maximumDistance = topology.maximumBridgeDistance;
  const SimulationBox& box = system.simulationBox;

  auto backbone = [&](std::span<const Atom> atoms, std::size_t unit)
  { return atoms[units[unit].backboneAtom].position; };

  // The anchor atoms a2 (unit 'site') and a6 (unit 'site + 4') of every integer molecule. They never
  // move; the partner criterion is a function of these positions only, so the criterion in the new
  // state follows from relabelling the two chains.
  struct AnchorPair
  {
    double3 a2{}, a6{};
  };
  std::vector<AnchorPair> anchors(numberOfMolecules);
  for (std::size_t m = numberOfFractional; m != numberOfMolecules; ++m)
  {
    const std::span<const Atom> atoms = system.spanOfMolecule(selectedComponent, m);
    anchors[m] = AnchorPair{backbone(atoms, site), backbone(atoms, site + 4)};
  }
  auto withinReach = [&](const double3& a, const double3& b)
  { return box.applyPeriodicBoundaryConditions(a - b).length() <= maximumDistance; };
  // Symmetric in the two chains: chain p can bridge onto the tail of q and q onto the tail of p.
  auto isCandidate = [&](const AnchorPair& p, const AnchorPair& q)
  { return withinReach(p.a2, q.a6) && withinReach(q.a2, p.a6); };
  auto countCandidates = [&](std::size_t self, std::span<const AnchorPair> all)
  {
    std::size_t count = 0;
    for (std::size_t m = numberOfFractional; m != all.size(); ++m)
    {
      if (m != self && isCandidate(all[self], all[m])) ++count;
    }
    return count;
  };

  // The chain's own trimer must satisfy the reach criterion, so that the reverse move (which bridges
  // the partner's head onto this trimer's tail) finds the partner among its candidates.
  const AnchorPair& anchorsI = anchors[selectedMolecule];
  if ((anchorsI.a2 - anchorsI.a6).length() > maximumDistance) return std::nullopt;

  std::vector<std::size_t> candidates{};
  for (std::size_t m = numberOfFractional; m != numberOfMolecules; ++m)
  {
    if (m != selectedMolecule && isCandidate(anchorsI, anchors[m])) candidates.push_back(m);
  }
  if (candidates.empty()) return std::nullopt;
  const std::size_t partnerMolecule =
      candidates[static_cast<std::size_t>(random.uniform() * static_cast<double>(candidates.size()))];
  const AnchorPair& anchorsJ = anchors[partnerMolecule];
  if ((anchorsJ.a2 - anchorsJ.a6).length() > maximumDistance) return std::nullopt;
  const std::size_t candidatesI = candidates.size();
  const std::size_t candidatesJ = countCandidates(partnerMolecule, anchors);

  // Molecules are stored unwrapped: shift the partner by the lattice vector that brings its a6 to the
  // minimum image of the selected chain's a2, so that both new chains are contiguous in space.
  const double3 separation = anchorsJ.a6 - anchorsI.a2;
  const double3 shift = box.applyPeriodicBoundaryConditions(separation) - separation;
  if ((anchorsJ.a2 + shift - anchorsI.a6).length() > maximumDistance) return std::nullopt;

  std::span<Atom> atomsI = system.spanOfMolecule(selectedComponent, selectedMolecule);
  std::span<Atom> atomsJ = system.spanOfMolecule(selectedComponent, partnerMolecule);
  const std::vector<Atom> oldI(atomsI.begin(), atomsI.end());
  std::vector<Atom> oldJ(atomsJ.begin(), atomsJ.end());
  for (Atom& atom : oldJ) atom.position += shift;

  // The two closures: each chain's head dimer (a1, a2) bridges onto the other chain's tail (a6, a7,
  // a8) with the internal geometry of its own excised trimer.
  const std::optional<std::size_t> a8Unit =
      (site + 6 < units.size()) ? std::optional<std::size_t>(site + 6) : std::nullopt;
  auto anchorsOf = [&](std::span<const Atom> atoms)
  {
    return Bridging::Anchors{
        .a1 = backbone(atoms, site - 1),
        .a2 = backbone(atoms, site),
        .a6 = backbone(atoms, site + 4),
        .a7 = backbone(atoms, site + 5),
        .a8 = a8Unit.has_value() ? std::optional<double3>(backbone(atoms, a8Unit.value())) : std::nullopt};
  };
  auto trimerOf = [&](std::span<const Atom> atoms)
  {
    return ConcertedRotation::Trimer{backbone(atoms, site + 1), backbone(atoms, site + 2), backbone(atoms, site + 3)};
  };

  const Bridging::Anchors oldAnchorsI = anchorsOf(oldI);
  const Bridging::Anchors oldAnchorsJ = anchorsOf(oldJ);
  const ConcertedRotation::Trimer oldTrimerI = trimerOf(oldI);
  const ConcertedRotation::Trimer oldTrimerJ = trimerOf(oldJ);
  const Bridging::Anchors newAnchorsI{oldAnchorsI.a1, oldAnchorsI.a2, oldAnchorsJ.a6, oldAnchorsJ.a7, oldAnchorsJ.a8};
  const Bridging::Anchors newAnchorsJ{oldAnchorsJ.a1, oldAnchorsJ.a2, oldAnchorsI.a6, oldAnchorsI.a7, oldAnchorsI.a8};

  const std::optional<Bridging::Closure> closureI = Bridging::rebridge(random, oldAnchorsI, oldTrimerI, newAnchorsI);
  if (!closureI.has_value()) return std::nullopt;
  const std::optional<Bridging::Closure> closureJ = Bridging::rebridge(random, oldAnchorsJ, oldTrimerJ, newAnchorsJ);
  if (!closureJ.has_value()) return std::nullopt;

  // Trial chains: tails (units site + 4 and beyond) exchanged slot by slot, trimers re-bridged.
  std::vector<Atom> trialI = oldI;
  std::vector<Atom> trialJ = oldJ;
  for (std::size_t k = site + 4; k != units.size(); ++k)
  {
    trialI[units[k].backboneAtom].position = oldJ[units[k].backboneAtom].position;
    trialJ[units[k].backboneAtom].position = oldI[units[k].backboneAtom].position;
    for (std::size_t atom : units[k].sideAtoms)
    {
      trialI[atom].position = oldJ[atom].position;
      trialJ[atom].position = oldI[atom].position;
    }
  }
  Bridging::placeTrimer(trialI, oldI, units, site, oldAnchorsI, oldTrimerI, newAnchorsI, closureI->trimer);
  Bridging::placeTrimer(trialJ, oldJ, units, site, oldAnchorsJ, oldTrimerJ, newAnchorsJ, closureJ->trimer);

  // Partner-selection probabilities in the new state: only the anchor pairs of the two chains change
  // (each keeps its a2 and takes the other's a6). Both remain candidates of each other.
  std::vector<AnchorPair> newAnchors = anchors;
  newAnchors[selectedMolecule] = AnchorPair{anchorsI.a2, anchorsJ.a6 + shift};
  newAnchors[partnerMolecule] = AnchorPair{anchorsJ.a2 + shift, anchorsI.a6};
  const std::size_t newCandidatesI = countCandidates(selectedMolecule, newAnchors);
  const std::size_t newCandidatesJ = countCandidates(partnerMolecule, newAnchors);
  const double logSelection =
      std::log(1.0 / static_cast<double>(newCandidatesI) + 1.0 / static_cast<double>(newCandidatesJ)) -
      std::log(1.0 / static_cast<double>(candidatesI) + 1.0 / static_cast<double>(candidatesJ));

  if (system.insideBlockedPockets(component, trialI) || system.insideBlockedPockets(component, trialJ))
  {
    return std::nullopt;
  }

  // The atoms whose positions change: the two trimers with their side groups. Every other atom keeps
  // its position (and, for the framework, the external field and the other molecules, its identity
  // does not matter), so the non-bonded differences with the surroundings involve these atoms only.
  std::vector<Atom> movedOld{}, movedNew{};
  auto collectTrimerAtoms = [&](std::vector<Atom>& destination, const std::vector<Atom>& source)
  {
    for (std::size_t m = 1; m <= 3; ++m)
    {
      const Component::BridgingTopology::Unit& unit = units[site + m];
      destination.push_back(source[unit.backboneAtom]);
      for (std::size_t atom : unit.sideAtoms) destination.push_back(source[atom]);
    }
  };
  collectTrimerAtoms(movedOld, oldI);
  collectTrimerAtoms(movedOld, oldJ);
  collectTrimerAtoms(movedNew, trialI);
  collectTrimerAtoms(movedNew, trialJ);

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
                  system.spanOfFrameworkAtoms(), movedNew, movedOld);
            });
  if (!frameworkMolecule.has_value()) return std::nullopt;

  // Molecule-molecule: the trimers against every other molecule, plus the interaction between the
  // two chains themselves (pairs change between intra- and intermolecular with the tail exchange, so
  // the full pair interaction of the two chains is evaluated before and after).
  std::optional<RunningEnergy> interMolecule =
      timed(system, component, move, Move::Timing::MoleculeMolecule,
            [&]() -> std::optional<RunningEnergy>
            {
              RunningEnergy sum{};
              for (std::size_t c = 0; c != system.components.size(); ++c)
              {
                for (std::size_t m = 0; m != system.numberOfMoleculesPerComponent[c]; ++m)
                {
                  if (c == selectedComponent && (m == selectedMolecule || m == partnerMolecule)) continue;
                  std::optional<RunningEnergy> difference = Interactions::computeInterMolecularEnergyDifference(
                      system.forceField, system.simulationBox, system.spanOfMolecule(c, m), movedNew, movedOld);
                  if (!difference.has_value()) return std::nullopt;
                  sum += difference.value();
                }
              }
              std::optional<RunningEnergy> pairNew = Interactions::computeInterMolecularEnergyDifference(
                  system.forceField, system.simulationBox, trialI, trialJ, std::span<const Atom>{});
              if (!pairNew.has_value()) return std::nullopt;
              std::optional<RunningEnergy> pairOld = Interactions::computeInterMolecularEnergyDifference(
                  system.forceField, system.simulationBox, oldI, std::span<const Atom>{}, oldJ);
              if (!pairOld.has_value()) return std::nullopt;
              return sum + pairNew.value() + pairOld.value();
            });
  if (!interMolecule.has_value()) return std::nullopt;

  // Ewald: the Fourier sum only sees the trimer displacements (the exchanged tails contribute
  // identical terms on both sides), but the intramolecular exclusion terms change with the chain
  // membership, so both chains are passed in full.
  std::vector<Atom> allNew = trialI;
  allNew.insert(allNew.end(), trialJ.begin(), trialJ.end());
  std::vector<Atom> allOld = oldI;
  allOld.insert(allOld.end(), oldJ.begin(), oldJ.end());
  RunningEnergy ewaldFourierEnergy =
      timed(system, component, move, Move::Timing::Ewald,
            [&]
            {
              return Interactions::energyDifferenceEwaldFourier(system.eik_x, system.eik_y, system.eik_z, system.eik_xy,
                                                                system.storedEik, system.trialEik, system.forceField,
                                                                system.simulationBox, allNew, allOld);
            });

  // Intramolecular: bonded terms across the seams and the intramolecular non-bonded terms of the
  // relabelled chains. Recomputing all internal terms of both chains keeps the bookkeeping exact.
  RunningEnergy internalDifference = component.intraMolecularPotentials.computeInternalEnergies(trialI) +
                                     component.intraMolecularPotentials.computeInternalEnergies(trialJ) -
                                     component.intraMolecularPotentials.computeInternalEnergies(oldI) -
                                     component.intraMolecularPotentials.computeInternalEnergies(oldJ);

  // Tail corrections cancel exactly: the multiset of pseudo-atom types in the box is unchanged.
  RunningEnergy energyDifference = externalFieldMolecule.value() + frameworkMolecule.value() + interMolecule.value() +
                                   ewaldFourierEnergy + internalDifference;

  component.mc_moves_statistics.addConstructed(move);

  const double logAcceptance =
      -system.beta * energyDifference.potentialEnergy() + closureI->logWeight + closureJ->logWeight + logSelection;
  if (random.uniform() < std::exp(logAcceptance))
  {
    component.mc_moves_statistics.addAccepted(move);

    Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.trialEik);

    std::copy(trialI.cbegin(), trialI.cend(), atomsI.begin());
    std::copy(trialJ.cbegin(), trialJ.cend(), atomsJ.begin());
    system.moleculeData[system.moleculeIndexOfComponent(selectedComponent, selectedMolecule)].centerOfMassPosition =
        component.computeCenterOfMass(atomsI);
    system.moleculeData[system.moleculeIndexOfComponent(selectedComponent, partnerMolecule)].centerOfMassPosition =
        component.computeCenterOfMass(atomsJ);

    return energyDifference;
  }
  return std::nullopt;
}
