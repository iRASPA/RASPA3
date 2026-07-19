module;

module mc_moves_translation_rotation_smart_mc;

import std;

import component;
import molecule;
import atom;
import atom_dynamics;
import double3;
import simd_quatd;
import randomnumbers;
import system;
import running_energy;
import interactions_framework_molecule;
import interactions_intermolecular;
import interactions_ewald;
import interactions_external_field;
import interactions_polarization;
import mc_moves_move_types;

// Fills per-atom energy gradients on the selected molecule (framework + inter + Ewald). A single
// gradient evaluation provides both the net force (for the translational drift) and the torque (for
// the rotational drift). Polarization forces are omitted from the drift; the polarization energy
// change still enters Delta U. The bias only needs to be a deterministic function of the
// configuration for the Metropolis-Hastings correction to keep the move exact.
static void computeMoleculeGradients(
    System &system, std::span<const Atom> selectedAtoms,
    const std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>> &structureFactor,
    std::vector<AtomDynamics> &dynamics)
{
  dynamics.assign(selectedAtoms.size(), AtomDynamics{});

  Interactions::computeFrameworkMoleculeGradient(system.forceField, system.simulationBox,
                                                 system.spanOfFrameworkAtoms(), selectedAtoms, dynamics,
                                                 system.interpolationGrids);

  Interactions::computeInterMolecularGradientMolecule(system.forceField, system.simulationBox,
                                                      system.spanOfMoleculeAtoms(), selectedAtoms, dynamics);

  Interactions::computeEwaldFourierGradientSingleMolecule(system.eik_x, system.eik_y, system.eik_z, system.eik_xy,
                                                          structureFactor, system.forceField, system.simulationBox,
                                                          selectedAtoms, dynamics);
}

// Net force on the molecule: F = sum_i -grad_i.
static double3 computeMoleculeForce(std::span<const AtomDynamics> dynamics)
{
  double3 gradient{};
  for (const AtomDynamics &d : dynamics) gradient += d.gradient;
  return -gradient;
}

// Lab-frame torque about the rotation pivot: tau = sum_i (r_i - pivot) x F_i with F = -grad.
// Rigid molecules rotate about the COM; flexible molecules about the starting bead (as Component::rotate).
static double3 computeMoleculeTorque(const Component &component, const Molecule &molecule,
                                     std::span<const Atom> selectedAtoms, std::span<const AtomDynamics> dynamics)
{
  const double3 pivot =
      component.rigid ? molecule.centerOfMassPosition : selectedAtoms[component.startingBead].position;

  double3 torque{};
  for (std::size_t i = 0; i != selectedAtoms.size(); ++i)
  {
    torque += double3::cross(selectedAtoms[i].position - pivot, -dynamics[i].gradient);
  }
  return torque;
}

static simd_quatd quaternionFromAngularDisplacement(double3 deltaPhi)
{
  const double angle = deltaPhi.length();
  if (angle < 1.0e-14)
  {
    return simd_quatd(0.0, 0.0, 0.0, 1.0);
  }
  return simd_quatd::fromAxisAngle(angle, deltaPhi * (1.0 / angle));
}

std::optional<RunningEnergy> MC_Moves::translationRotationSmartMCMove(RandomNumber &random, System &system,
                                                                      std::size_t selectedComponent,
                                                                      std::size_t selectedMolecule)
{
  std::span<Atom> molecule_atoms = system.spanOfMolecule(selectedComponent, selectedMolecule);
  Molecule &molecule = system.moleculeData[system.moleculeIndexOfComponent(selectedComponent, selectedMolecule)];

  std::chrono::steady_clock::time_point time_begin, time_end;
  Move::Types move = Move::Types::TranslationRotationSmartMC;
  Component &component = system.components[selectedComponent];

  component.mc_moves_statistics.addTrial(move, 0);
  component.mc_moves_statistics.addTrial(move, 1);

  // Monoatomic species have no rotational degrees of freedom; use translation smart MC instead.
  if (molecule.numberOfAtoms < 2)
  {
    return std::nullopt;
  }

  // Two independent step sizes: 'sigmaTranslation' (Angstrom) for the Gaussian part of the trial
  // displacement, 'sigmaRotation' (radians) for the Gaussian part of the trial rotation vector.
  // Each drift coefficient follows the smart-MC relation b = beta * sigma^2 / 2.
  double sigmaTranslation = component.mc_moves_statistics.getMaxChange(move, 0);
  double sigmaRotation = component.mc_moves_statistics.getMaxChange(move, 1);
  double bTranslation = 0.5 * system.beta * sigmaTranslation * sigmaTranslation;
  double bRotation = 0.5 * system.beta * sigmaRotation * sigmaRotation;

  // Force and torque on the selected molecule in the current configuration, from a single
  // gradient evaluation (uses the maintained structure factor S_old).
  time_begin = std::chrono::steady_clock::now();
  std::vector<AtomDynamics> dynamicsOld;
  computeMoleculeGradients(system, molecule_atoms, system.storedEik, dynamicsOld);
  double3 forceOld = computeMoleculeForce(dynamicsOld);
  double3 torqueOld = computeMoleculeTorque(component, molecule, molecule_atoms, dynamicsOld);
  time_end = std::chrono::steady_clock::now();
  component.mc_moves_cputime[move][Move::Timing::Integration] += (time_end - time_begin);
  system.mc_moves_cputime[move][Move::Timing::Integration] += (time_end - time_begin);

  // Biased trial displacement and rotation vector: Delta r = b_t * F_old + sigma_t * xi_t and
  // Delta phi = b_r * tau_old + sigma_r * xi_r, with independent Gaussian noise vectors.
  double3 xiTranslation(random.Gaussian(), random.Gaussian(), random.Gaussian());
  double3 displacement = bTranslation * forceOld + sigmaTranslation * xiTranslation;

  double3 xiRotation(random.Gaussian(), random.Gaussian(), random.Gaussian());
  double3 angularDisplacement = bRotation * torqueOld + sigmaRotation * xiRotation;
  simd_quatd rotationQuaternion = quaternionFromAngularDisplacement(angularDisplacement);

  // Translate first, then rotate about the moved pivot. The pivot is a material point of the
  // molecule, so the two operations commute and the reverse proposal (-Delta r, -Delta phi)
  // maps the trial state exactly back onto the old state.
  std::pair<Molecule, std::vector<Atom>> translatedMolecule =
      component.translate(molecule, molecule_atoms, displacement);
  std::pair<Molecule, std::vector<Atom>> trialMolecule =
      component.rotate(translatedMolecule.first, translatedMolecule.second, rotationQuaternion);

  // Reject when the trial position lies inside a blocked pocket. The configuration was not modified.
  if (system.insideBlockedPockets(component, trialMolecule.second))
  {
    return std::nullopt;
  }

  std::vector<double3> electricFieldMoleculeNew(molecule_atoms.size());
  std::vector<double3> electricFieldMoleculeOld(molecule_atoms.size());

  // Exact energy difference of the move, using the same incremental machinery as the regular moves.
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

  // Force and torque on the selected molecule in the trial configuration (uses the updated structure
  // factor S_new). The other molecules are unchanged, so the current 'spanOfMoleculeAtoms()' (still
  // holding the old positions of the moved molecule) is used; the self-interaction is excluded by molecule id.
  time_begin = std::chrono::steady_clock::now();
  std::vector<AtomDynamics> dynamicsNew;
  computeMoleculeGradients(system, trialMolecule.second, system.trialEik, dynamicsNew);
  double3 forceNew = computeMoleculeForce(dynamicsNew);
  double3 torqueNew = computeMoleculeTorque(component, trialMolecule.first, trialMolecule.second, dynamicsNew);
  time_end = std::chrono::steady_clock::now();
  component.mc_moves_cputime[move][Move::Timing::Integration] += (time_end - time_begin);
  system.mc_moves_cputime[move][Move::Timing::Integration] += (time_end - time_begin);

  component.mc_moves_statistics.addConstructed(move, 0);
  component.mc_moves_statistics.addConstructed(move, 1);

  // Metropolis-Hastings acceptance with the asymmetric-proposal (Hastings) correction. The translational
  // and rotational noises are independent, so the forward/reverse proposal ratio factorizes and the
  // log-bias is the sum of the two per-channel terms:
  // log[q(n->o)/q(o->n)] = (|Delta r - b_t F_old|^2 - |Delta r + b_t F_new|^2) / (2 sigma_t^2)
  //                      + (|Delta phi - b_r tau_old|^2 - |Delta phi + b_r tau_new|^2) / (2 sigma_r^2).
  double3 forwardDeviationTranslation = displacement - bTranslation * forceOld;
  double3 reverseDeviationTranslation = displacement + bTranslation * forceNew;
  double logBiasTranslation =
      (forwardDeviationTranslation.length_squared() - reverseDeviationTranslation.length_squared()) /
      (2.0 * sigmaTranslation * sigmaTranslation);

  double3 forwardDeviationRotation = angularDisplacement - bRotation * torqueOld;
  double3 reverseDeviationRotation = angularDisplacement + bRotation * torqueNew;
  double logBiasRotation = (forwardDeviationRotation.length_squared() - reverseDeviationRotation.length_squared()) /
                           (2.0 * sigmaRotation * sigmaRotation);

  double logBias = logBiasTranslation + logBiasRotation;

  if (random.uniform() < std::exp(-system.beta * energyDifference.potentialEnergy() + logBias))
  {
    component.mc_moves_statistics.addAccepted(move, 0);
    component.mc_moves_statistics.addAccepted(move, 1);

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
