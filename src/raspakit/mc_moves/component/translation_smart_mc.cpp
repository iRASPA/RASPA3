module;

module mc_moves_translation_smart_mc;

import std;

import component;
import molecule;
import atom;
import atom_dynamics;
import double3;
import randomnumbers;
import system;
import running_energy;
import interactions_framework_molecule;
import interactions_intermolecular;
import interactions_ewald;
import interactions_external_field;
import interactions_polarization;
import mc_moves_move_types;

// Computes the total (center-of-mass) force on the selected molecule at the configuration given by
// 'selectedAtoms'. Contributions: framework (real space + interpolation grid), inter-molecular real space, and
// the reciprocal-space Ewald force evaluated from the supplied total structure factor 'structureFactor'
// (S_old = system.storedEik for the current configuration, S_new = system.trialEik for the trial configuration).
//
// Intramolecular self/exclusion Ewald corrections are omitted on purpose: they generate internal forces that
// cancel in the net force of a rigid translation. Polarization forces are not included in the drift; the bias
// only needs to be a deterministic function of the configuration for the Metropolis-Hastings correction to keep
// the move exact, and the exact polarization contribution is still part of the accepted energy difference.
static double3 computeMoleculeForce(
    System &system, std::span<const Atom> selectedAtoms,
    const std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>> &structureFactor)
{
  std::vector<AtomDynamics> dynamics(selectedAtoms.size());

  Interactions::computeFrameworkMoleculeGradient(system.forceField, system.simulationBox,
                                                 system.spanOfFrameworkAtoms(), selectedAtoms, dynamics,
                                                 system.interpolationGrids);

  Interactions::computeInterMolecularGradientMolecule(system.forceField, system.simulationBox,
                                                      system.spanOfMoleculeAtoms(), selectedAtoms, dynamics);

  Interactions::computeEwaldFourierGradientSingleMolecule(system.eik_x, system.eik_y, system.eik_z, system.eik_xy,
                                                          structureFactor, system.forceField, system.simulationBox,
                                                          selectedAtoms, dynamics);

  double3 gradient{};
  for (const AtomDynamics &d : dynamics) gradient += d.gradient;

  // The force is the negative gradient.
  return -gradient;
}

std::optional<RunningEnergy> MC_Moves::translationSmartMCMove(RandomNumber &random, System &system,
                                                               std::size_t selectedComponent,
                                                               std::size_t selectedMolecule)
{
  std::span<Atom> molecule_atoms = system.spanOfMolecule(selectedComponent, selectedMolecule);
  Molecule &molecule = system.moleculeData[system.moleculeIndexOfComponent(selectedComponent, selectedMolecule)];

  std::chrono::steady_clock::time_point time_begin, time_end;
  Move::Types move = Move::Types::TranslationSmartMC;
  Component &component = system.components[selectedComponent];

  component.mc_moves_statistics.addTrial(move);

  // 'sigma' is the standard deviation of the Gaussian part of the trial displacement; the drift coefficient
  // follows the smart-MC relation b = beta * sigma^2 / 2.
  double sigma = component.mc_moves_statistics.getMaxChange(move);
  double b = 0.5 * system.beta * sigma * sigma;

  // Force on the selected molecule in the current configuration (uses the maintained structure factor S_old).
  time_begin = std::chrono::steady_clock::now();
  double3 forceOld = computeMoleculeForce(system, molecule_atoms, system.storedEik);
  time_end = std::chrono::steady_clock::now();
  component.mc_moves_cputime[move][Move::Timing::Integration] += (time_end - time_begin);
  system.mc_moves_cputime[move][Move::Timing::Integration] += (time_end - time_begin);

  // Biased trial displacement of the center of mass: Delta r = b * F_old + sigma * xi.
  double3 xi(random.Gaussian(), random.Gaussian(), random.Gaussian());
  double3 displacement = b * forceOld + sigma * xi;

  std::pair<Molecule, std::vector<Atom>> trialMolecule =
      component.translate(molecule, molecule_atoms, displacement);

  // Reject when the trial position lies inside a blocked pocket. The configuration was not modified.
  if (system.insideBlockedPockets(component, trialMolecule.second))
  {
    return std::nullopt;
  }

  std::vector<double3> electricFieldMoleculeNew(molecule_atoms.size());
  std::vector<double3> electricFieldMoleculeOld(molecule_atoms.size());

  // Exact energy difference of the move, using the same incremental machinery as the regular translation move.
  time_begin = std::chrono::steady_clock::now();
  std::optional<RunningEnergy> externalFieldMolecule = Interactions::computeExternalFieldEnergyDifference(
      system.hasExternalField, system.forceField, system.simulationBox, system.externalFieldInterpolationGrid,
      trialMolecule.second, molecule_atoms);
  time_end = std::chrono::steady_clock::now();
  component.mc_moves_cputime[move][Move::Timing::ExternalFieldMolecule] += (time_end - time_begin);
  system.mc_moves_cputime[move][Move::Timing::ExternalFieldMolecule] += (time_end - time_begin);
  if (!externalFieldMolecule.has_value()) return std::nullopt;

  time_begin = std::chrono::steady_clock::now();
  std::optional<RunningEnergy> frameworkMolecule;
  if (system.forceField.computePolarization)
  {
    frameworkMolecule = Interactions::computeFrameworkMoleculeEnergyDifference(
        system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
        system.spanOfFrameworkAtoms(), electricFieldMoleculeNew, electricFieldMoleculeOld, trialMolecule.second,
        molecule_atoms);
  }
  else
  {
    frameworkMolecule = Interactions::computeFrameworkMoleculeEnergyDifference(
        system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
        system.spanOfFrameworkAtoms(), trialMolecule.second, molecule_atoms);
  }
  time_end = std::chrono::steady_clock::now();
  component.mc_moves_cputime[move][Move::Timing::FrameworkMolecule] += (time_end - time_begin);
  system.mc_moves_cputime[move][Move::Timing::FrameworkMolecule] += (time_end - time_begin);
  if (!frameworkMolecule.has_value()) return std::nullopt;

  time_begin = std::chrono::steady_clock::now();
  std::vector<double3> electricFieldNeighborDelta;
  std::optional<RunningEnergy> interMolecule;
  if (system.forceField.computePolarization && !system.forceField.omitInterPolarization)
  {
    electricFieldNeighborDelta.assign(system.spanOfMoleculeAtoms().size(), double3(0.0, 0.0, 0.0));
    interMolecule = Interactions::computeInterMolecularPolarizationElectricFieldDifference(
        system.forceField, system.simulationBox, electricFieldNeighborDelta, electricFieldMoleculeNew,
        electricFieldMoleculeOld, system.spanOfMoleculeAtoms(), trialMolecule.second, molecule_atoms);
  }
  else
  {
    interMolecule = Interactions::computeInterMolecularEnergyDifference(
        system.forceField, system.simulationBox, system.spanOfMoleculeAtoms(), trialMolecule.second, molecule_atoms);
  }
  time_end = std::chrono::steady_clock::now();
  component.mc_moves_cputime[move][Move::Timing::MoleculeMolecule] += (time_end - time_begin);
  system.mc_moves_cputime[move][Move::Timing::MoleculeMolecule] += (time_end - time_begin);
  if (!interMolecule.has_value()) return std::nullopt;

  // Ewald energy difference. This also fills 'system.trialEik' with the updated total structure factor S_new
  // (S_new = S_old - S_molecule_old + S_molecule_new); most terms cancel, so only the moved molecule contributes.
  time_begin = std::chrono::steady_clock::now();
  RunningEnergy ewaldFourierEnergy;
  if (system.forceField.computePolarization)
  {
    ewaldFourierEnergy = Interactions::energyDifferenceEwaldFourier(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
        system.trialEik, system.forceField, system.simulationBox, electricFieldMoleculeNew, electricFieldMoleculeOld,
        trialMolecule.second, molecule_atoms);
  }
  else
  {
    ewaldFourierEnergy = Interactions::energyDifferenceEwaldFourier(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.trialEik, system.forceField,
        system.simulationBox, trialMolecule.second, molecule_atoms);
  }
  time_end = std::chrono::steady_clock::now();
  component.mc_moves_cputime[move][Move::Timing::Ewald] += (time_end - time_begin);
  system.mc_moves_cputime[move][Move::Timing::Ewald] += (time_end - time_begin);

  RunningEnergy polarizationDifference;
  if (system.forceField.computePolarization)
  {
    polarizationDifference = Interactions::computePolarizationEnergyDifference(
        system.forceField, electricFieldMoleculeNew, electricFieldMoleculeOld, trialMolecule.second, molecule_atoms);

    if (!system.forceField.omitInterPolarization)
    {
      polarizationDifference += Interactions::computePolarizationEnergyNeighborDifference(
          system.forceField, system.spanOfMoleculeElectricField(), electricFieldNeighborDelta,
          system.spanOfMoleculeAtoms());
    }
  }

  RunningEnergy energyDifference = externalFieldMolecule.value() + frameworkMolecule.value() + interMolecule.value() +
                                   ewaldFourierEnergy + polarizationDifference;

  // Force on the selected molecule in the trial configuration (uses the updated structure factor S_new). The
  // other molecules are unchanged, so the current 'spanOfMoleculeAtoms()' (still holding the old positions of
  // the moved molecule) is used; the self-interaction is excluded by molecule id.
  time_begin = std::chrono::steady_clock::now();
  double3 forceNew = computeMoleculeForce(system, trialMolecule.second, system.trialEik);
  time_end = std::chrono::steady_clock::now();
  component.mc_moves_cputime[move][Move::Timing::Integration] += (time_end - time_begin);
  system.mc_moves_cputime[move][Move::Timing::Integration] += (time_end - time_begin);

  component.mc_moves_statistics.addConstructed(move);

  // Metropolis-Hastings acceptance with the asymmetric-proposal (Hastings) correction. The forward proposal
  // (o -> n) draws Delta r ~ N(b * F_old, sigma^2); the reverse proposal (n -> o) draws -Delta r ~ N(b * F_new,
  // sigma^2). log[q(n->o)/q(o->n)] = (|Delta r - b F_old|^2 - |Delta r + b F_new|^2) / (2 sigma^2).
  double3 forwardDeviation = displacement - b * forceOld;
  double3 reverseDeviation = displacement + b * forceNew;
  double logBias = (forwardDeviation.length_squared() - reverseDeviation.length_squared()) / (2.0 * sigma * sigma);

  if (random.uniform() < std::exp(-system.beta * energyDifference.potentialEnergy() + logBias))
  {
    component.mc_moves_statistics.addAccepted(move);

    Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.trialEik);

    std::copy(trialMolecule.second.cbegin(), trialMolecule.second.cend(), molecule_atoms.begin());
    molecule = trialMolecule.first;

    if (system.forceField.computePolarization)
    {
      if (!system.forceField.omitInterPolarization)
      {
        std::span<double3> storedElectricField = system.spanOfMoleculeElectricField();
        for (std::size_t i = 0; i < storedElectricField.size(); ++i)
        {
          storedElectricField[i] += electricFieldNeighborDelta[i];
        }
      }
      std::span<double3> electricFieldMolecule = system.spanElectricFieldOld(selectedComponent, selectedMolecule);
      std::copy(electricFieldMoleculeNew.begin(), electricFieldMoleculeNew.end(), electricFieldMolecule.begin());
    }

    return energyDifference;
  }

  return std::nullopt;
}
