module;

module mc_moves_pivot_cbmc;

import std;

import component;
import atom;
import molecule;
import double3;
import double3x3;
import simd_quatd;
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

std::optional<RunningEnergy> MC_Moves::pivotCBMCMove(RandomNumber &random, System &system,
                                                     std::size_t selectedComponent, std::size_t selectedMolecule)
{
  Move::Types move = Move::Types::PivotCBMC;
  Component &component = system.components[selectedComponent];

  // Restriction: the polarization energy update is not implemented. Molecules without a valid
  // pivot axis (fully rigid molecules, pure rings, diatomics) cannot use this move. Neither
  // condition changes during a simulation, so the attempts are not counted as trials.
  const std::vector<Component::PivotBond> &pivotBonds = component.pivotBonds();
  if (pivotBonds.empty() || system.forceField.computePolarization)
  {
    return std::nullopt;
  }

  // Mixed step sizes, as in the plain pivot: channel 0 draws the trial angles within the adaptive
  // window, channel 1 (a fraction pivotCBMCRandomizationFraction of the attempts, window pinned at
  // pi) randomizes them completely. Separate statistics keep the randomizations from biasing the
  // adaptive window of the small-step channel.
  std::size_t channel = (random.uniform() < component.pivotCBMCRandomizationFraction) ? 1uz : 0uz;

  component.mc_moves_statistics.addTrial(move, channel);

  std::span<Atom> molecule_atoms = system.spanOfMolecule(selectedComponent, selectedMolecule);
  Molecule &molecule = system.moleculeData[system.moleculeIndexOfComponent(selectedComponent, selectedMolecule)];

  // Select a pivot axis uniformly from the precomputed valid bonds. The rotated part is the
  // smaller side of the molecule, chosen deterministically, so the reverse move selects the same
  // bond and the same part.
  const std::size_t bondIndex = static_cast<std::size_t>(random.uniform() * static_cast<double>(pivotBonds.size()));
  component.mc_moves_statistics.addSubTrial(move, bondIndex, pivotBonds.size());
  const Component::PivotBond &pivotBond = pivotBonds[bondIndex];
  const std::array<std::size_t, 2> &bond = pivotBond.bond;
  const std::vector<std::size_t> &rotatedAtoms = pivotBond.rotatedAtoms;
  const std::span<const std::size_t> movedIndices(rotatedAtoms);

  const std::size_t numberOfTrials = std::max(1uz, component.pivotCBMCNumberOfTrialAngles);
  const double maxAngle = component.mc_moves_statistics.getMaxChange(move, channel);
  const double3 axisOrigin = molecule_atoms[bond[0]].position;
  const double3 rotationAxis = (molecule_atoms[bond[1]].position - axisOrigin).normalized();

  // Rigid rotation of the selected part about the bond axis by 'angle' (relative to the current
  // configuration), written into 'trialAtoms' (which holds a copy of the molecule).
  auto rotateInto = [&](std::vector<Atom> &trialAtoms, double angle)
  {
    simd_quatd q = simd_quatd::fromAxisAngle(angle, rotationAxis);
    double3x3 rotationMatrix = double3x3::buildRotationMatrixInverse(q);
    for (std::size_t index : rotatedAtoms)
    {
      trialAtoms[index].position = axisOrigin + rotationMatrix * (molecule_atoms[index].position - axisOrigin);
    }
  };

  std::vector<Atom> movedOld;
  movedOld.reserve(rotatedAtoms.size());
  for (std::size_t index : rotatedAtoms) movedOld.push_back(molecule_atoms[index]);
  const RunningEnergy internalOld = component.intraMolecularPotentials.computeInternalEnergies(molecule_atoms);

  // The bias energy of a trial relative to the current configuration: everything except the Ewald
  // Fourier part. Only the rotated part interacts differently with its surroundings, so the
  // external-field, framework and intermolecular differences are evaluated for those atoms alone;
  // the intramolecular terms use the whole molecule (the torsions through the pivot bond and the
  // intramolecular non-bonded terms change). Nullopt signals an overlap or a blocked pocket, i.e.
  // a trial of zero weight.
  auto biasEnergyDifference = [&](const std::vector<Atom> &trialAtoms) -> std::optional<RunningEnergy>
  {
    if (system.insideBlockedPockets(component, trialAtoms)) return std::nullopt;

    std::vector<Atom> movedNew;
    movedNew.reserve(rotatedAtoms.size());
    for (std::size_t index : rotatedAtoms) movedNew.push_back(trialAtoms[index]);

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

    std::optional<RunningEnergy> interMolecule =
        timed(system, component, move, Move::Timing::MoleculeMolecule,
              [&]
              {
                return Interactions::computeInterMolecularEnergyDifference(
                    system.forceField, system.simulationBox, system.cellList(), system.spanOfMoleculeAtoms(),
                    movedNew, movedOld);
              });
    if (!interMolecule.has_value()) return std::nullopt;

    RunningEnergy internalDifference =
        component.intraMolecularPotentials.computeInternalEnergies(trialAtoms) - internalOld;

    return externalFieldMolecule.value() + frameworkMolecule.value() + interMolecule.value() + internalDifference +
           system.crossLinkEnergyDifference(selectedComponent, selectedMolecule, trialAtoms, molecule_atoms);
  };

  // Log-sum-exp of the Boltzmann factors of a set of trials; nullopt when every trial has zero weight.
  // (Written without infinities: the code is compiled with -ffast-math.)
  auto logRosenbluthWeight = [](std::span<const std::optional<double>> logWeights) -> std::optional<double>
  {
    std::optional<double> maximum{};
    for (const std::optional<double> &w : logWeights)
    {
      if (w.has_value()) maximum = maximum.has_value() ? std::max(maximum.value(), w.value()) : w.value();
    }
    if (!maximum.has_value()) return std::nullopt;
    double sum = 0.0;
    for (const std::optional<double> &w : logWeights)
    {
      if (w.has_value()) sum += std::exp(w.value() - maximum.value());
    }
    return maximum.value() + std::log(sum);
  };

  // Forward: k trial angles, their bias energies and log-weights -beta dU_i.
  std::vector<Atom> trialAtoms(molecule_atoms.begin(), molecule_atoms.end());
  std::vector<double> trialAngles(numberOfTrials);
  std::vector<std::optional<RunningEnergy>> trialEnergies(numberOfTrials);
  std::vector<std::optional<double>> trialLogWeights(numberOfTrials);
  for (std::size_t i = 0; i != numberOfTrials; ++i)
  {
    trialAngles[i] = maxAngle * 2.0 * (random.uniform() - 0.5);
    rotateInto(trialAtoms, trialAngles[i]);
    trialEnergies[i] = biasEnergyDifference(trialAtoms);
    if (trialEnergies[i].has_value())
    {
      trialLogWeights[i] = -system.beta * trialEnergies[i]->potentialEnergy();
    }
  }
  const std::optional<double> logWeightNew = logRosenbluthWeight(trialLogWeights);
  if (!logWeightNew.has_value()) return std::nullopt;

  // Select one trial with probability w_i / W_new.
  std::size_t selected = numberOfTrials;
  {
    const double target = random.uniform();
    double cumulative = 0.0;
    for (std::size_t i = 0; i != numberOfTrials; ++i)
    {
      if (!trialLogWeights[i].has_value()) continue;
      cumulative += std::exp(trialLogWeights[i].value() - logWeightNew.value());
      selected = i;
      if (target < cumulative) break;
    }
  }
  const double selectedAngle = trialAngles[selected];
  const RunningEnergy biasDifference = trialEnergies[selected].value();

  // Reverse: the old configuration (zero energy difference, unit weight) plus k - 1 trial angles
  // drawn about the NEW configuration, i.e. selectedAngle + psi_j relative to the old one.
  std::vector<std::optional<double>> reverseLogWeights(numberOfTrials);
  reverseLogWeights[0] = 0.0;
  for (std::size_t j = 1; j < numberOfTrials; ++j)
  {
    const double psi = maxAngle * 2.0 * (random.uniform() - 0.5);
    rotateInto(trialAtoms, selectedAngle + psi);
    std::optional<RunningEnergy> energy = biasEnergyDifference(trialAtoms);
    if (energy.has_value())
    {
      reverseLogWeights[j] = -system.beta * energy->potentialEnergy();
    }
  }
  const double logWeightOld = logRosenbluthWeight(reverseLogWeights).value();

  // The selected trial configuration and its Ewald Fourier correction.
  rotateInto(trialAtoms, selectedAngle);
  RunningEnergy ewaldFourierEnergy =
      timed(system, component, move, Move::Timing::Ewald,
            [&]
            {
              return Interactions::energyDifferenceEwaldFourierMovedAtoms(
                  system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.trialEik,
                  system.forceField, system.simulationBox, trialAtoms, molecule_atoms, movedIndices);
            });

  RunningEnergy energyDifference = biasDifference + ewaldFourierEnergy;

  component.mc_moves_statistics.addConstructed(move, channel);
  component.mc_moves_statistics.addSubConstructed(move, bondIndex);

  // Acceptance: the Rosenbluth ratio times the Metropolis factor of the part not in the bias.
  const double logAcceptance =
      logWeightNew.value() - logWeightOld - system.beta * ewaldFourierEnergy.potentialEnergy();
  if (logAcceptance >= 0.0 || random.uniform() < std::exp(logAcceptance))
  {
    component.mc_moves_statistics.addAccepted(move, channel);
    component.mc_moves_statistics.addSubAccepted(move, bondIndex);

    Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.trialEik);

    std::copy(trialAtoms.cbegin(), trialAtoms.cend(), molecule_atoms.begin());
    system.cellListAtomsMoved(molecule_atoms);
    molecule.centerOfMassPosition = component.computeCenterOfMass(molecule_atoms);

    return energyDifference;
  }
  return std::nullopt;
}
