module;

module mc_moves_parallel_tempering_swap;

import std;

import component;
import atom;
import framework;
import double3;
import double3x3;
import simd_quatd;
import simulationbox;
import cbmc;
import randomnumbers;
import system;
import energy_factor;
import energy_status;
import energy_status_inter;
import running_energy;
import property_lambda_probability_histogram;
import property_widom;
import property_loading;
import averages;
import forcefield;
import interactions_framework_molecule;
import interactions_intermolecular;
import interactions_ewald;
import interactions_external_field;
import mc_moves_move_types;

namespace
{

bool sameAtomDefinition(const Atom& atomA, const Atom& atomB)
{
  return atomA.position == atomB.position && atomA.charge == atomB.charge && atomA.type == atomB.type;
}

bool sameHamiltonian(const ForceField& forceFieldA, const ForceField& forceFieldB)
{
  // Automatic Ewald wave-vector bounds are box-derived caches, not Hamiltonian parameters.
  // ForceField::temperature is likewise the state point used to derive the pair coefficients;
  // temperature-dependent Hamiltonians still compare unequal through those derived coefficients.
  ForceField normalizedB = forceFieldB;
  normalizedB.temperature = forceFieldA.temperature;
  if (forceFieldA.automaticEwald && forceFieldB.automaticEwald)
  {
    normalizedB.EwaldAlpha = forceFieldA.EwaldAlpha;
    normalizedB.numberOfWaveVectors = forceFieldA.numberOfWaveVectors;
    normalizedB.reciprocalIntegerCutOffSquared = forceFieldA.reciprocalIntegerCutOffSquared;
    normalizedB.reciprocalCutOffSquared = forceFieldA.reciprocalCutOffSquared;
  }
  return forceFieldA == normalizedB;
}

bool compatibleMobileTopology(const System& systemA, const System& systemB)
{
  auto hasFractionalSlots = [](const System& system)
  {
    return std::ranges::any_of(system.numberOfFractionalMoleculesPerComponent,
                               [](std::size_t count) { return count != 0; });
  };

  // This implementation deliberately supports only rigid, whole-molecule replicas.
  // Flexible intramolecular definitions, reaction coordinates, and fractional lambda
  // state do not yet have a complete cross-replica compatibility comparator. Rejecting
  // them before the random draw is conservative and preserves detailed balance.
  if (!sameHamiltonian(systemA.forceField, systemB.forceField) || systemA.hasExternalField || systemB.hasExternalField ||
      systemA.components.size() != systemB.components.size() ||
      systemA.numberOfFrameworkAtoms != systemB.numberOfFrameworkAtoms ||
      !systemA.reactions.list.empty() || !systemB.reactions.list.empty() ||
      hasFractionalSlots(systemA) || hasFractionalSlots(systemB))
  {
    return false;
  }

  if (systemA.framework.has_value() != systemB.framework.has_value())
  {
    return false;
  }
  if (systemA.framework.has_value())
  {
    const Framework& frameworkA = systemA.framework.value();
    const Framework& frameworkB = systemB.framework.value();
    if (frameworkA.name != frameworkB.name || frameworkA.simulationBox != frameworkB.simulationBox ||
        frameworkA.numberOfUnitCells != frameworkB.numberOfUnitCells ||
        frameworkA.atoms.size() != frameworkB.atoms.size() || systemA.simulationBox != systemB.simulationBox)
    {
      return false;
    }
    for (std::size_t i = 0; i < frameworkA.atoms.size(); ++i)
    {
      if (!sameAtomDefinition(frameworkA.atoms[i], frameworkB.atoms[i]))
      {
        return false;
      }
    }
  }

  for (std::size_t componentId = 0; componentId < systemA.components.size(); ++componentId)
  {
    const Component& componentA = systemA.components[componentId];
    const Component& componentB = systemB.components[componentId];
    if (!componentA.rigid || !componentB.rigid || componentA.name != componentB.name ||
        componentA.atoms.size() != componentB.atoms.size())
    {
      return false;
    }
    for (std::size_t atomId = 0; atomId < componentA.atoms.size(); ++atomId)
    {
      if (!sameAtomDefinition(componentA.atoms[atomId], componentB.atoms[atomId]))
      {
        return false;
      }
    }
  }
  return true;
}

template <typename T>
void swapMobileTail(std::vector<T>& dataA, std::size_t fixedSizeA, std::vector<T>& dataB, std::size_t fixedSizeB)
{
  std::vector<T> mobileA(std::make_move_iterator(dataA.begin() + static_cast<std::ptrdiff_t>(fixedSizeA)),
                         std::make_move_iterator(dataA.end()));
  std::vector<T> mobileB(std::make_move_iterator(dataB.begin() + static_cast<std::ptrdiff_t>(fixedSizeB)),
                         std::make_move_iterator(dataB.end()));
  dataA.erase(dataA.begin() + static_cast<std::ptrdiff_t>(fixedSizeA), dataA.end());
  dataB.erase(dataB.begin() + static_cast<std::ptrdiff_t>(fixedSizeB), dataB.end());
  dataA.insert(dataA.end(), std::make_move_iterator(mobileB.begin()), std::make_move_iterator(mobileB.end()));
  dataB.insert(dataB.end(), std::make_move_iterator(mobileA.begin()), std::make_move_iterator(mobileA.end()));
}

void rebuildConfigurationDerivedState(System& system)
{
  system.forceField.initializeEwaldParameters(system.simulationBox);
  system.eik_x.clear();
  system.eik_y.clear();
  system.eik_z.clear();
  system.eik_xy.clear();
  system.storedEik.clear();
  system.fixedFrameworkStoredEik.clear();
  system.totalEik.clear();
  system.precomputeTotalRigidEnergy();
  system.runningEnergies = system.computeTotalEnergies();
  system.totalEik = system.storedEik;
  system.CoulombicFourierEnergySingleIon = Interactions::computeEwaldFourierEnergySingleIon(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.forceField, system.simulationBox,
      double3(0.0, 0.0, 0.0), 1.0);
  system.loadings =
      LoadingData(system.components.size(), system.numberOfIntegerMoleculesPerComponent, system.simulationBox);
  system.updateMoleculeAtomInformation();
  system.computeNumberOfPseudoAtoms();
  system.netCharge = system.netChargeFramework + system.netChargeAdsorbates;
  system.checkMoleculeIds();
}

}  // namespace

std::optional<std::pair<RunningEnergy, RunningEnergy>> MC_Moves::ParallelTemperingSwap(RandomNumber &random,
                                                                                       System &systemA, System &systemB)
{
  std::chrono::steady_clock::time_point time_begin, time_end;
  Move::Types move = Move::Types::ParallelTempering;

  // Update swap move counts for both systems
  systemA.mc_moves_statistics.addTrial(move);

  // A complete cross-Hamiltonian evaluator is not available here. Reject incompatible
  // Hamiltonians and topologies before drawing an acceptance variate rather than mixing
  // partial cross energies with full running energies.
  if (!compatibleMobileTopology(systemA, systemB))
  {
    return std::nullopt;
  }

  double acc = std::exp((systemB.beta - systemA.beta) *
                        (systemB.runningEnergies.potentialEnergy() - systemA.runningEnergies.potentialEnergy()));

  if (systemA.pressure != systemB.pressure)
  {
    /// Ref: "Hyper-parallel tempering Monte Carlo: Application to the Lennard-Jones fluid and the
    /// restricted primitive model",  G. Yan and J.J. de Pablo, JCP, 111(21): 9509-9516, 1999

    // Adjust acceptance probability for pressure differences
    time_begin = std::chrono::steady_clock::now();
    const std::ptrdiff_t moleculeDifference =
        static_cast<std::ptrdiff_t>(systemB.numberOfIntegerMolecules()) -
        static_cast<std::ptrdiff_t>(systemA.numberOfIntegerMolecules());
    acc *= std::pow(systemB.pressure / systemA.pressure, static_cast<double>(moleculeDifference));
    time_end = std::chrono::steady_clock::now();

    systemA.mc_moves_cputime[move][Move::Timing::Fugacity] += (time_end - time_begin);
  }

  // Update constructed move counts for both systems
  systemA.mc_moves_statistics.addConstructed(move);

  // Apply acceptance/rejection rule
  if (random.uniform() < acc)
  {
    // Update accepted move counts for both systems
    systemA.mc_moves_statistics.addAccepted(move);

    // Swap configuration-owned state. Thermodynamic state, force fields, move controls,
    // accumulated statistics, and property samplers remain attached to their replicas.
    swapMobileTail(systemA.atomData, systemA.numberOfFrameworkAtoms, systemB.atomData, systemB.numberOfFrameworkAtoms);
    swapMobileTail(systemA.atomDynamics, systemA.numberOfFrameworkAtoms, systemB.atomDynamics,
                   systemB.numberOfFrameworkAtoms);
    std::swap(systemA.moleculeData, systemB.moleculeData);
    if (!systemA.framework.has_value())
    {
      std::swap(systemA.simulationBox, systemB.simulationBox);
    }
    std::swap(systemA.numberOfMoleculesPerComponent, systemB.numberOfMoleculesPerComponent);
    std::swap(systemA.numberOfIntegerMoleculesPerComponent, systemB.numberOfIntegerMoleculesPerComponent);
    swapMobileTail(systemA.electricPotential, systemA.numberOfFrameworkAtoms, systemB.electricPotential,
                   systemB.numberOfFrameworkAtoms);
    swapMobileTail(systemA.electricField, systemA.numberOfFrameworkAtoms, systemB.electricField,
                   systemB.numberOfFrameworkAtoms);
    swapMobileTail(systemA.electricFieldNew, systemA.numberOfFrameworkAtoms, systemB.electricFieldNew,
                   systemB.numberOfFrameworkAtoms);
    std::swap(systemA.netChargeAdsorbates, systemB.netChargeAdsorbates);
    std::swap(systemA.netChargePerComponent, systemB.netChargePerComponent);
    std::swap(systemA.translationalCenterOfMassConstraint, systemB.translationalCenterOfMassConstraint);
    std::swap(systemA.translationalDegreesOfFreedom, systemB.translationalDegreesOfFreedom);
    std::swap(systemA.rotationalDegreesOfFreedom, systemB.rotationalDegreesOfFreedom);

    rebuildConfigurationDerivedState(systemA);
    rebuildConfigurationDerivedState(systemB);

    return std::make_pair(systemA.runningEnergies, systemB.runningEnergies);
  }

  return std::nullopt;
}
