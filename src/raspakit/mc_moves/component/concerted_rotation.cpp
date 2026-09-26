module;

module mc_moves_concerted_rotation;

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

namespace
{
// Two trimers are the same closure when every atom agrees to this distance [Angstrom]; the solver
// converges to ~1e-12, the tolerance only has to separate distinct solutions.
constexpr double sameTrimerTolerance = 1e-6;

bool sameTrimer(const ConcertedRotation::Trimer &trimer, const double3 &a3, const double3 &a4, const double3 &a5)
{
  return (trimer.a3 - a3).length() < sameTrimerTolerance && (trimer.a4 - a4).length() < sameTrimerTolerance &&
         (trimer.a5 - a5).length() < sameTrimerTolerance;
}
}  // namespace

std::optional<RunningEnergy> MC_Moves::concertedRotationMove(RandomNumber &random, System &system,
                                                             std::size_t selectedComponent,
                                                             std::size_t selectedMolecule)
{
  Move::Types move = Move::Types::ConcertedRotation;
  Component &component = system.components[selectedComponent];

  // Restriction: the polarization energy update is not implemented. Molecules without a valid
  // window cannot use this move. Neither condition changes during a simulation, so the attempts
  // are not counted as trials.
  const std::vector<Component::ConcertedRotationWindow> &windows = component.concertedRotationWindows();
  if (windows.empty() || system.forceField.computePolarization)
  {
    return std::nullopt;
  }

  // Mixed step sizes for the driver angle (see the pivot and crankshaft moves): channel 0 is the
  // adaptive small-step channel, channel 1 the full randomization pinned at pi.
  std::size_t channel = (random.uniform() < component.concertedRotationRandomizationFraction) ? 1uz : 0uz;

  component.mc_moves_statistics.addTrial(move, channel);

  std::span<Atom> molecule_atoms = system.spanOfMolecule(selectedComponent, selectedMolecule);
  Molecule &molecule = system.moleculeData[system.moleculeIndexOfComponent(selectedComponent, selectedMolecule)];

  // The window is chosen uniformly; the reverse move uses the same window with the negated driver
  // angle, so the window and angle proposals are symmetric.
  const Component::ConcertedRotationWindow &window =
      windows[static_cast<std::size_t>(random.uniform() * static_cast<double>(windows.size()))];
  const std::array<std::size_t, 8> &backbone = window.backbone;

  std::array<double3, 8> oldPositions{};
  for (std::size_t i = 0; i != 8; ++i) oldPositions[i] = molecule_atoms[backbone[i]].position;
  const std::optional<double3> a8 =
      window.a8.has_value() ? std::optional<double3>(molecule_atoms[window.a8.value()].position) : std::nullopt;

  // The driver: rotate a2 about the a0-a1 axis. Bond a1a2 and bend a0a1a2 are invariant.
  double maxAngle = component.mc_moves_statistics.getMaxChange(move, channel);
  double driverAngle = maxAngle * 2.0 * (random.uniform() - 0.5);
  const double3 driverAxis = (oldPositions[1] - oldPositions[0]).normalized();
  const double3 newA2 = ConcertedRotation::rotateAboutAxis(oldPositions[1], driverAxis, driverAngle, oldPositions[2]);

  // Re-bridge the trimer with the window's invariant internal geometry, measured from the old state.
  const std::span<const double3, 7> oldBackbone(oldPositions.data() + 1, 7);
  const ConcertedRotation::BackboneGeometry geometry = ConcertedRotation::BackboneGeometry::fromPositions(oldBackbone);

  const std::vector<ConcertedRotation::Trimer> forwardSolutions =
      ConcertedRotation::rebridge(oldPositions[1], newA2, oldPositions[6], oldPositions[7], geometry);
  if (forwardSolutions.empty())
  {
    return std::nullopt;
  }
  const ConcertedRotation::Trimer &chosen =
      forwardSolutions[static_cast<std::size_t>(random.uniform() * static_cast<double>(forwardSolutions.size()))];

  std::array<double3, 8> newPositions = oldPositions;
  newPositions[2] = newA2;
  newPositions[3] = chosen.a3;
  newPositions[4] = chosen.a4;
  newPositions[5] = chosen.a5;

  // The reverse move rotates a2 back and solves the same closure; it must find the old trimer among
  // its solutions (their number enters the acceptance rule). A miss is a numerical failure of the
  // root search: reject rather than break reversibility.
  const std::vector<ConcertedRotation::Trimer> backwardSolutions =
      ConcertedRotation::rebridge(oldPositions[1], oldPositions[2], oldPositions[6], oldPositions[7], geometry);
  const bool oldTrimerFound =
      std::any_of(backwardSolutions.begin(), backwardSolutions.end(), [&](const ConcertedRotation::Trimer &trimer)
                  { return sameTrimer(trimer, oldPositions[3], oldPositions[4], oldPositions[5]); });
  if (!oldTrimerFound)
  {
    return std::nullopt;
  }

  const double jacobianOld = ConcertedRotation::closureJacobian(oldBackbone, a8);
  const double jacobianNew =
      ConcertedRotation::closureJacobian(std::span<const double3, 7>(newPositions.data() + 1, 7), a8);
  if (!(jacobianOld > 0.0) || !(jacobianNew > 0.0))
  {
    return std::nullopt;
  }

  // Trial positions: the four backbone atoms, and their side groups carried rigidly with the local
  // frame (previous backbone atom, atom, next backbone atom), which the move maps congruently.
  std::vector<Atom> trialAtoms(molecule_atoms.begin(), molecule_atoms.end());
  for (std::size_t i = 2; i <= 5; ++i)
  {
    trialAtoms[backbone[i]].position = newPositions[i];
    const std::vector<std::size_t> &substituent = window.substituents[i - 2];
    if (substituent.empty()) continue;
    const double3x3 oldFrame =
        ConcertedRotation::localFrame(oldPositions[i], oldPositions[i - 1], oldPositions[i + 1]);
    const double3x3 newFrame =
        ConcertedRotation::localFrame(newPositions[i], newPositions[i - 1], newPositions[i + 1]);
    for (std::size_t atom : substituent)
    {
      trialAtoms[atom].position = ConcertedRotation::transformRigidly(oldPositions[i], oldFrame, newPositions[i],
                                                                      newFrame, molecule_atoms[atom].position);
    }
  }

  if (system.insideBlockedPockets(component, trialAtoms))
  {
    return std::nullopt;
  }

  // Compute external field energy contribution
  std::optional<RunningEnergy> externalFieldMolecule =
      timed(system, component, move, Move::Timing::ExternalFieldMolecule,
            [&]
            {
              return Interactions::computeExternalFieldEnergyDifference(
                  system.hasExternalField, system.forceField, system.simulationBox,
                  system.externalFieldInterpolationGrid, trialAtoms, molecule_atoms);
            });
  if (!externalFieldMolecule.has_value()) return std::nullopt;

  // Compute framework-molecule energy contribution
  std::optional<RunningEnergy> frameworkMolecule =
      timed(system, component, move, Move::Timing::FrameworkMolecule,
            [&]
            {
              return Interactions::computeFrameworkMoleculeEnergyDifference(
                  system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
                  system.spanOfFrameworkAtoms(), trialAtoms, molecule_atoms);
            });
  if (!frameworkMolecule.has_value()) return std::nullopt;

  // Compute molecule-molecule energy contribution
  std::optional<RunningEnergy> interMolecule =
      timed(system, component, move, Move::Timing::MoleculeMolecule,
            [&]
            {
              return Interactions::computeInterMolecularEnergyDifference(
                  system.forceField, system.simulationBox, system.spanOfMoleculeAtoms(), trialAtoms, molecule_atoms);
            });
  if (!interMolecule.has_value()) return std::nullopt;

  // Compute Ewald energy contribution
  RunningEnergy ewaldFourierEnergy =
      timed(system, component, move, Move::Timing::Ewald,
            [&]
            {
              return Interactions::energyDifferenceEwaldFourier(system.eik_x, system.eik_y, system.eik_z, system.eik_xy,
                                                                system.storedEik, system.trialEik, system.forceField,
                                                                system.simulationBox, trialAtoms, molecule_atoms);
            });

  // Intramolecular energy contribution: the window's bond lengths and backbone bends are invariant,
  // but seven torsions, the bends and torsions to side groups at a1 and a6, and the intramolecular
  // non-bonded energy change. Recomputing all internal terms keeps the bookkeeping exact.
  RunningEnergy internalDifference = component.intraMolecularPotentials.computeInternalEnergies(trialAtoms) -
                                     component.intraMolecularPotentials.computeInternalEnergies(molecule_atoms);

  RunningEnergy energyDifference = externalFieldMolecule.value() + frameworkMolecule.value() + interMolecule.value() +
                                   ewaldFourierEnergy + internalDifference;

  component.mc_moves_statistics.addConstructed(move, channel);

  // Acceptance: Metropolis on the energy, times the solution-count ratio and the closure-Jacobian
  // ratio of the deterministic map (see the module documentation).
  const double logAcceptance = -system.beta * energyDifference.potentialEnergy() +
                               std::log(static_cast<double>(forwardSolutions.size())) -
                               std::log(static_cast<double>(backwardSolutions.size())) + std::log(jacobianOld) -
                               std::log(jacobianNew);
  if (random.uniform() < std::exp(logAcceptance))
  {
    component.mc_moves_statistics.addAccepted(move, channel);

    Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.trialEik);

    std::copy(trialAtoms.cbegin(), trialAtoms.cend(), molecule_atoms.begin());
    molecule.centerOfMassPosition = component.computeCenterOfMass(molecule_atoms);

    return energyDifference;
  }
  return std::nullopt;
}

