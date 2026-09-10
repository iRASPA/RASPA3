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

bool anyNonzero(const std::vector<std::size_t>& counts)
{
  return std::ranges::any_of(counts, [](std::size_t count) { return count != 0; });
}

// Pair-swap, group-swap, Gibbs and reaction fractional slots need extra discrete state
// and standard-state factors that this move does not evaluate. GC and pair-GC
// (SwapCFCMC / SwapCBCFCMC) slots are supported when both replicas match.
bool hasUnsupportedFractionalSlots(const System& system)
{
  return anyNonzero(system.numberOfPairSwapFractionalMoleculesPerComponent_CFCMC) ||
         anyNonzero(system.numberOfPairSwapCBFractionalMoleculesPerComponent_CFCMC) ||
         anyNonzero(system.numberOfGroupSwapFractionalMoleculesPerComponent_CFCMC) ||
         anyNonzero(system.numberOfGroupSwapCBFractionalMoleculesPerComponent_CFCMC) ||
         anyNonzero(system.numberOfGibbsSwapFractionalMoleculesPerComponent_CFCMC) ||
         anyNonzero(system.numberOfGibbsFractionalMoleculesPerComponent_CFCMC) ||
         anyNonzero(system.numberOfParallelReactionFractionalMoleculesPerComponent_CFCMC) ||
         anyNonzero(system.numberOfSerialReactionFractionalMoleculesPerComponent_CFCMC);
}

bool matchingGCFractionalLayout(const System& systemA, const System& systemB)
{
  return systemA.numberOfFractionalMoleculesPerComponent == systemB.numberOfFractionalMoleculesPerComponent &&
         systemA.numberOfGCFractionalMoleculesPerComponent_CFCMC ==
             systemB.numberOfGCFractionalMoleculesPerComponent_CFCMC &&
         systemA.numberOfPairGCFractionalMoleculesPerComponent_CFCMC ==
             systemB.numberOfPairGCFractionalMoleculesPerComponent_CFCMC;
}

bool matchingLambdaGrids(const System& systemA, const System& systemB)
{
  for (std::size_t componentId = 0; componentId < systemA.components.size(); ++componentId)
  {
    const Component& componentA = systemA.components[componentId];
    const Component& componentB = systemB.components[componentId];
    if (!componentA.hasFractionalMolecule && !componentB.hasFractionalMolecule) continue;
    if (componentA.hasFractionalMolecule != componentB.hasFractionalMolecule) return false;
    if (componentA.lambdaGC.numberOfSamplePoints != componentB.lambdaGC.numberOfSamplePoints) return false;
    if (componentA.lambdaGC.biasFactor.size() != componentB.lambdaGC.biasFactor.size()) return false;
    if (componentA.lambdaGC.currentBin >= componentA.lambdaGC.biasFactor.size()) return false;
    if (componentB.lambdaGC.currentBin >= componentB.lambdaGC.biasFactor.size()) return false;
  }
  return true;
}

bool compatibleMobileTopology(const System& systemA, const System& systemB)
{
  if (!sameHamiltonian(systemA.forceField, systemB.forceField) || systemA.hasExternalField || systemB.hasExternalField ||
      systemA.components.size() != systemB.components.size() ||
      systemA.numberOfFrameworkAtoms != systemB.numberOfFrameworkAtoms ||
      !systemA.reactions.list.empty() || !systemB.reactions.list.empty() ||
      hasUnsupportedFractionalSlots(systemA) || hasUnsupportedFractionalSlots(systemB) ||
      !matchingGCFractionalLayout(systemA, systemB) || !matchingLambdaGrids(systemA, systemB))
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

std::optional<double> tmmcLogBias(const System& replica, const System& configuration)
{
  if (!replica.tmmc.doTMMC || !replica.tmmc.useBias) return 0.0;
  if (!replica.tmmc.useTMBias && !replica.tmmc.useWangLandau) return 0.0;
  if (replica.components.empty() || configuration.components.empty()) return 0.0;

  const std::size_t moleculeCount = configuration.numberOfIntegerMoleculesPerComponent.front();
  const std::size_t lambdaBin =
      replica.tmmc.lambdaChain() ? configuration.components.front().lambdaGC.currentBin : 0uz;
  const std::size_t index = replica.tmmc.chainIndex(moleculeCount, lambdaBin);
  if (index >= replica.tmmc.bias.size()) return std::nullopt;
  return replica.tmmc.bias[index];
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
  system.trialEik.clear();
  system.precomputeTotalRigidEnergy();
  system.runningEnergies = system.computeTotalEnergies();
  system.trialEik = system.storedEik;
  system.CoulombicFourierEnergySingleIon = Interactions::computeEwaldFourierEnergySingleIon(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.forceField, system.simulationBox,
      double3(0.0, 0.0, 0.0), 1.0);
  system.loadings =
      LoadingData(system.components.size(), system.numberOfIntegerMoleculesPerComponent, system.simulationBox);
  system.updateMoleculeAtomInformation();
  system.computeNumberOfPseudoAtoms();
  system.computeTailCorrectionCounts();
  system.netCharge = system.netChargeFramework + system.netChargeAdsorbates;
  system.checkMoleculeIds();
  if (system.tmmc.doTMMC && !system.components.empty())
  {
    system.tmmc.currentLambdaBin = system.components.front().lambdaGC.currentBin;
  }
}

}  // namespace

std::optional<double> MC_Moves::ParallelTemperingLogAcceptance(const System& systemA, const System& systemB)
{
  if (!compatibleMobileTopology(systemA, systemB))
  {
    return std::nullopt;
  }

  // Symmetric exchange of configurations X_A ↔ X_B between ensembles (β_A, f_A) and (β_B, f_B).
  // For a shared temperature-independent Hamiltonian the energy term is
  // (β_B − β_A)(U(X_B) − U(X_A)). The activity uses integer molecule counts only: the
  // fractional molecule is already in U and in the replica-local λ-bias.
  //
  //     log R = (β_B − β_A)(U_B − U_A)
  //           + Σ_i (N_B,i − N_A,i) log(a_A,i / a_B,i)
  //           + (β_B P_B − β_A P_A)(V_B − V_A)          [variable-cell / no-framework only]
  //           + Σ_q [B_A(λ_q(X_B)) − B_A(λ_q(X_A)) + B_B(λ_q(X_A)) − B_B(λ_q(X_B))]
  //           + B^{TM}_A(X_B) − B^{TM}_A(X_A) + B^{TM}_B(X_A) − B^{TM}_B(X_B)
  double logR = (systemB.beta - systemA.beta) *
                (systemB.runningEnergies.potentialEnergy() - systemA.runningEnergies.potentialEnergy());

  for (std::size_t componentId = 0; componentId < systemA.components.size(); ++componentId)
  {
    const std::ptrdiff_t moleculeDifference =
        static_cast<std::ptrdiff_t>(systemB.numberOfIntegerMoleculesPerComponent[componentId]) -
        static_cast<std::ptrdiff_t>(systemA.numberOfIntegerMoleculesPerComponent[componentId]);
    if (moleculeDifference == 0) continue;

    const Component& componentA = systemA.components[componentId];
    const Component& componentB = systemB.components[componentId];
    const double fugacityA = componentA.molFraction * componentA.fugacityCoefficient.value_or(1.0) * systemA.pressure;
    const double fugacityB = componentB.molFraction * componentB.fugacityCoefficient.value_or(1.0) * systemB.pressure;
    const double activityA = systemA.beta * fugacityA;
    const double activityB = systemB.beta * fugacityB;
    if (!(activityA > 0.0) || !(activityB > 0.0))
    {
      return std::nullopt;
    }
    logR += static_cast<double>(moleculeDifference) * (std::log(activityA) - std::log(activityB));
  }

  // Volume travels with the configuration only when there is no framework. The PV term
  // belongs to an isobaric ensemble; a fixed-framework μVT box keeps its cell.
  if (!systemA.framework.has_value())
  {
    logR += (systemB.beta * systemB.pressure - systemA.beta * systemA.pressure) *
            (systemB.simulationBox.volume - systemA.simulationBox.volume);
  }

  for (std::size_t componentId = 0; componentId < systemA.components.size(); ++componentId)
  {
    const Component& componentA = systemA.components[componentId];
    const Component& componentB = systemB.components[componentId];
    if (!componentA.hasFractionalMolecule) continue;
    const std::size_t binA = componentA.lambdaGC.currentBin;
    const std::size_t binB = componentB.lambdaGC.currentBin;
    logR += componentA.lambdaGC.biasFactor[binB] - componentA.lambdaGC.biasFactor[binA] +
            componentB.lambdaGC.biasFactor[binA] - componentB.lambdaGC.biasFactor[binB];
  }

  const std::optional<double> tmmcAOnA = tmmcLogBias(systemA, systemA);
  const std::optional<double> tmmcAOnB = tmmcLogBias(systemA, systemB);
  const std::optional<double> tmmcBOnB = tmmcLogBias(systemB, systemB);
  const std::optional<double> tmmcBOnA = tmmcLogBias(systemB, systemA);
  if (!tmmcAOnA.has_value() || !tmmcAOnB.has_value() || !tmmcBOnB.has_value() || !tmmcBOnA.has_value())
  {
    return std::nullopt;
  }
  logR += *tmmcAOnB - *tmmcAOnA + *tmmcBOnA - *tmmcBOnB;

  return logR;
}

std::optional<std::pair<RunningEnergy, RunningEnergy>> MC_Moves::ParallelTemperingSwap(RandomNumber &random,
                                                                                       System &systemA, System &systemB)
{
  std::chrono::steady_clock::time_point time_begin, time_end;
  Move::Types move = Move::Types::ParallelTempering;

  systemA.mc_moves_statistics.addTrial(move);

  time_begin = std::chrono::steady_clock::now();
  const std::optional<double> logAcceptance = ParallelTemperingLogAcceptance(systemA, systemB);
  time_end = std::chrono::steady_clock::now();
  systemA.mc_moves_cputime[move][Move::Timing::Fugacity] += (time_end - time_begin);

  if (!logAcceptance.has_value())
  {
    return std::nullopt;
  }

  systemA.mc_moves_statistics.addConstructed(move);

  const double logUniform = std::log(std::max(random.uniform(), std::numeric_limits<double>::min()));
  if (logUniform < *logAcceptance)
  {
    systemA.mc_moves_statistics.addAccepted(move);

    // Swap configuration-owned state. Thermodynamic state, force fields, learned λ/TMMC
    // biases, move controls, accumulated statistics, and property samplers stay put.
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
    std::swap(systemA.containsTheFractionalMolecule, systemB.containsTheFractionalMolecule);
    for (std::size_t componentId = 0; componentId < systemA.components.size(); ++componentId)
    {
      std::swap(systemA.components[componentId].lambdaGC.currentBin,
                systemB.components[componentId].lambdaGC.currentBin);
    }

    rebuildConfigurationDerivedState(systemA);
    rebuildConfigurationDerivedState(systemB);

    return std::make_pair(systemA.runningEnergies, systemB.runningEnergies);
  }

  return std::nullopt;
}
