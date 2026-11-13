module;

module mc_moves_rotation_smart_mc;

import std;

import ewald_structure_factor;

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
import mc_moves_cputime;

// Fills per-atom energy gradients on the selected molecule (framework + inter + Ewald), matching the
// cheap single-molecule path used by translation smart MC. Polarization forces are omitted from the
// drift; the polarization energy change still enters Delta U.
static void computeMoleculeGradients(
    System &system, std::span<const Atom> selectedAtoms,
    const EwaldStructureFactor &structureFactor,
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

double MC_Moves::logRotationProposalRatio(double3 angularDisplacement, double3 driftOld, double3 driftNew, double sigma)
{
  const double twoSigmaSquared = 2.0 * sigma * sigma;
  const double angle = angularDisplacement.length();
  if (angle < 1.0e-14)
  {
    // the identity rotation: the principal term alone (its aliases have no axis)
    return (driftOld.length_squared() - driftNew.length_squared()) / twoSigmaSquared;
  }
  const double3 axis = angularDisplacement * (1.0 / angle);

  // The rotation is reached by every alias sign * axis * (theta + 2 pi k) of the drawn angle, and the density sums
  // them in log space. The alias nearest the drift dominates, k ~ (sign axis . drift - theta) / 2 pi, and that is
  // where each sum is centred, four aliases either side (at a step size of 1.5 radians the fifth lies more than
  // twenty standard deviations out). Centring both on the drawn angle missed the dominant alias of the reverse move
  // whenever |theta + axis . driftNew| exceeded nine pi, as a large (repulsive) torque gives. No infinite sentinel:
  // the fast-math flags of the build fold comparisons with infinity.
  auto logTerm = [&](double theta, const double3& drift, double sign)
  {
    const double3 vector = axis * (sign * theta);
    return 2.0 * std::log(std::abs(theta)) - (vector - drift).length_squared() / twoSigmaSquared;
  };
  auto logAdd = [](double a, double b)
  {
    const double high = std::max(a, b);
    return high + std::log(std::exp(a - high) + std::exp(b - high));
  };
  auto logSum = [&](const double3& drift, double sign)
  {
    const double twoPi = 2.0 * std::numbers::pi;
    const double centre = std::round((sign * double3::dot(axis, drift) - angle) / twoPi);
    double sum = 0.0;
    bool first = true;
    for (int j = -4; j <= 4; ++j)
    {
      const double theta = angle + twoPi * (centre + static_cast<double>(j));
      if (std::abs(theta) < 1.0e-12) continue;  // a vanishing weight theta^2
      const double term = logTerm(theta, drift, sign);
      sum = first ? term : logAdd(sum, term);
      first = false;
    }
    return sum;
  };
  return logSum(driftNew, -1.0) - logSum(driftOld, 1.0);
}

std::optional<RunningEnergy> MC_Moves::rotationSmartMCMove(RandomNumber &random, System &system,
                                                           std::size_t selectedComponent,
                                                           std::size_t selectedMolecule)
{
  std::span<Atom> molecule_atoms = system.spanOfMolecule(selectedComponent, selectedMolecule);
  Molecule &molecule = system.moleculeData[system.moleculeIndexOfComponent(selectedComponent, selectedMolecule)];

  std::chrono::steady_clock::time_point time_begin, time_end;
  Move::Types move = Move::Types::RotationSmartMC;
  Component &component = system.components[selectedComponent];

  // Restriction: monoatomic species have no rotational degrees of freedom. This does not change during a
  // simulation, so the attempts are not counted as trials (counted, they would report a vanishing acceptance
  // and shrink the step for nothing).
  if (molecule.numberOfAtoms < 2)
  {
    return std::nullopt;
  }

  component.mc_moves_statistics.addTrial(move);

  // 'sigma' is the angular standard deviation (radians) of the Gaussian part; the drift follows
  // the smart-MC relation b = beta * sigma^2 / 2.
  double sigma = component.mc_moves_statistics.getMaxChange(move);
  double b = 0.5 * system.beta * sigma * sigma;

  time_begin = std::chrono::steady_clock::now();
  std::vector<AtomDynamics> dynamicsOld;
  computeMoleculeGradients(system, molecule_atoms, system.storedEik, dynamicsOld);
  double3 torqueOld = computeMoleculeTorque(component, molecule, molecule_atoms, dynamicsOld);
  time_end = std::chrono::steady_clock::now();
  component.mc_moves_cputime[move][Move::Timing::Integration] += (time_end - time_begin);
  system.mc_moves_cputime[move][Move::Timing::Integration] += (time_end - time_begin);

  // Biased trial rotation vector: Delta phi = b * tau_old + sigma * xi, applied as a quaternion.
  double3 xi(random.Gaussian(), random.Gaussian(), random.Gaussian());
  double3 angularDisplacement = b * torqueOld + sigma * xi;
  simd_quatd rotationQuaternion = quaternionFromAngularDisplacement(angularDisplacement);

  std::pair<Molecule, std::vector<Atom>> trialMolecule =
      component.rotate(molecule, molecule_atoms, rotationQuaternion);

  if (system.insideBlockedPockets(component, trialMolecule.second))
  {
    return std::nullopt;
  }

  std::vector<double3> electricFieldMoleculeNew(molecule_atoms.size());
  std::vector<double3> electricFieldMoleculeOld(molecule_atoms.size());

  std::optional<RunningEnergy> externalFieldMolecule =
      timed(system, component, move, Move::Timing::ExternalFieldMolecule,
            [&]
            {
              return Interactions::computeExternalFieldEnergyDifference(
                  system.hasExternalField, system.forceField, system.simulationBox,
                  system.externalFieldInterpolationGrid, trialMolecule.second, molecule_atoms);
            });
  if (!externalFieldMolecule.has_value()) return std::nullopt;

  std::optional<RunningEnergy> frameworkMolecule =
      timed(system, component, move, Move::Timing::FrameworkMolecule,
            [&]() -> std::optional<RunningEnergy>
            {
              if (system.forceField.computePolarization)
              {
                return Interactions::computeFrameworkMoleculeEnergyDifference(
                    system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
                    system.spanOfFrameworkAtoms(), electricFieldMoleculeNew, electricFieldMoleculeOld,
                    trialMolecule.second, molecule_atoms);
              }
              return Interactions::computeFrameworkMoleculeEnergyDifference(
                  system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
                  system.frameworkCellList(), system.spanOfFrameworkAtoms(), trialMolecule.second, molecule_atoms);
            });
  if (!frameworkMolecule.has_value()) return std::nullopt;

  time_begin = std::chrono::steady_clock::now();
  std::vector<double3> electricFieldNeighborDelta;
  std::optional<RunningEnergy> interMolecule;
  if (system.forceField.computePolarization && !system.forceField.omitInterPolarization)
  {
    electricFieldNeighborDelta.assign(system.spanOfMoleculeAtoms().size(), double3(0.0, 0.0, 0.0));
    interMolecule = Interactions::computeInterMolecularPolarizationElectricFieldDifference(
        system.forceField, system.simulationBox, electricFieldNeighborDelta, electricFieldMoleculeNew,
        electricFieldMoleculeOld, system.spanOfMoleculeAtoms(), trialMolecule.second, molecule_atoms,
        !system.crossLinks.moleculeIsLinked(selectedComponent, selectedMolecule));
  }
  else
  {
    interMolecule = Interactions::computeInterMolecularEnergyDifference(
        system.forceField, system.simulationBox, system.cellList(), system.spanOfMoleculeAtoms(), trialMolecule.second,
        molecule_atoms, !system.crossLinks.moleculeIsLinked(selectedComponent, selectedMolecule));
  }
  time_end = std::chrono::steady_clock::now();
  component.mc_moves_cputime[move][Move::Timing::MoleculeMolecule] += (time_end - time_begin);
  system.mc_moves_cputime[move][Move::Timing::MoleculeMolecule] += (time_end - time_begin);
  if (!interMolecule.has_value()) return std::nullopt;

  RunningEnergy ewaldFourierEnergy =
      timed(system, component, move, Move::Timing::Ewald,
            [&]() -> RunningEnergy
            {
              if (system.forceField.computePolarization)
              {
                return Interactions::energyDifferenceEwaldFourier(
                    system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik,
                    system.storedEik, system.trialEik, system.forceField, system.simulationBox, system.components,
                    electricFieldMoleculeNew, electricFieldMoleculeOld, trialMolecule.second, molecule_atoms);
              }
              return Interactions::energyDifferenceEwaldFourier(
                  system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.trialEik,
                  system.forceField, system.simulationBox, system.components, trialMolecule.second, molecule_atoms);
            });

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

  RunningEnergy energyDifference =
      externalFieldMolecule.value() + frameworkMolecule.value() + interMolecule.value() + ewaldFourierEnergy +
      polarizationDifference +
      system.crossLinkEnergyDifference(selectedComponent, selectedMolecule, trialMolecule.second, molecule_atoms);

  time_begin = std::chrono::steady_clock::now();
  std::vector<AtomDynamics> dynamicsNew;
  computeMoleculeGradients(system, trialMolecule.second, system.trialEik, dynamicsNew);
  double3 torqueNew = computeMoleculeTorque(component, trialMolecule.first, trialMolecule.second, dynamicsNew);
  time_end = std::chrono::steady_clock::now();
  component.mc_moves_cputime[move][Move::Timing::Integration] += (time_end - time_begin);
  system.mc_moves_cputime[move][Move::Timing::Integration] += (time_end - time_begin);

  component.mc_moves_statistics.addConstructed(move);

  // Metropolis-Hastings correction for the asymmetric rotation proposal (the aliased rotation vectors included)
  const double logBias = logRotationProposalRatio(angularDisplacement, b * torqueOld, b * torqueNew, sigma);

  if (random.uniform() < std::exp(-system.beta * energyDifference.potentialEnergy() + logBias))
  {
    component.mc_moves_statistics.addAccepted(move);

    Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.trialEik);

    std::copy(trialMolecule.second.cbegin(), trialMolecule.second.cend(), molecule_atoms.begin());
    system.cellListAtomsMoved(molecule_atoms);
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
