module;

module mc_moves_parallel_tempering_swap;

import std;

import ewald_structure_factor;

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
import vdwparameters;
import interactions_framework_molecule;
import interactions_intermolecular;
import interactions_internal;
import molecule;
import interactions_ewald;
import interactions_cross_link;
import interactions_external_field;
import mc_moves_move_types;
import mc_moves_cputime;
import integrators_compute;
import integrators_update;
import thermobarostat;
import elastic_constants;

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
  // An automatic Coulomb cut-off (half the box, the default) is a box-derived cache as well: replicas at different
  // volumes (NPT replica exchange) would otherwise never compare equal and every swap would be refused.
  if (forceFieldA.cutOffCoulombAutomatic && forceFieldB.cutOffCoulombAutomatic)
  {
    normalizedB.cutOffCoulomb = forceFieldA.cutOffCoulomb;
  }
  // The same for the automatic Van der Waals cut-offs, whose shifts, tail corrections and switching constants follow
  // them: after the first volume move every swap of an NPT run with 'CutOffVDW': 'Auto' was refused otherwise.
  bool vdwCutOffNormalized = false;
  if (forceFieldA.cutOffMoleculeVDWAutomatic && forceFieldB.cutOffMoleculeVDWAutomatic)
  {
    vdwCutOffNormalized = vdwCutOffNormalized || normalizedB.cutOffMoleculeVDW != forceFieldA.cutOffMoleculeVDW;
    normalizedB.cutOffMoleculeVDW = forceFieldA.cutOffMoleculeVDW;
  }
  if (forceFieldA.cutOffFrameworkVDWAutomatic && forceFieldB.cutOffFrameworkVDWAutomatic)
  {
    vdwCutOffNormalized = vdwCutOffNormalized || normalizedB.cutOffFrameworkVDW != forceFieldA.cutOffFrameworkVDW;
    normalizedB.cutOffFrameworkVDW = forceFieldA.cutOffFrameworkVDW;
  }
  if (vdwCutOffNormalized) normalizedB.preComputeCutOffDependentParameters();
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
    // A pinned lambda (fixed-lambda thermodynamic integration) is a Hamiltonian parameter of the replica,
    // not configuration state: the fractional atoms carry their scaling factors with the configuration,
    // so an exchange between replicas pinned at different lambdas would leave replica A holding atoms
    // scaled at lambda_B. Replicas at different pinned lambdas exchange through 'LambdaExchange'.
    if (componentA.fixedLambdaBin != componentB.fixedLambdaBin) return false;
    if (componentA.fixedLambdaBin.has_value() && componentA.fixedLambdaCoordinate != componentB.fixedLambdaCoordinate)
    {
      return false;
    }
    if (componentA.lambdaGC.numberOfSamplePoints != componentB.lambdaGC.numberOfSamplePoints) return false;
    if (componentA.lambdaGC.biasFactor.size() != componentB.lambdaGC.biasFactor.size()) return false;
    if (componentA.lambdaGC.currentBin >= componentA.lambdaGC.biasFactor.size()) return false;
    if (componentB.lambdaGC.currentBin >= componentB.lambdaGC.biasFactor.size()) return false;
  }
  return true;
}

// Solute tempering: the replicas share the Hamiltonian except for the scaled pair parameters and
// charges of the tempered component, which the exchange accounts for explicitly.
// A replica that was never scaled (lambda = 1, no component recorded) is compatible with a scaled one.
std::optional<std::size_t> temperedComponent(const System& systemA, const System& systemB)
{
  return systemA.soluteTemperingComponent.has_value() ? systemA.soluteTemperingComponent
                                                      : systemB.soluteTemperingComponent;
}

bool sameSoluteTempering(const System& systemA, const System& systemB)
{
  if (systemA.soluteTemperingComponent.has_value() && systemB.soluteTemperingComponent.has_value())
  {
    return systemA.soluteTemperingComponent == systemB.soluteTemperingComponent;
  }
  return true;
}

bool sameHamiltonianUpToSoluteTempering(const System& systemA, const System& systemB)
{
  const std::optional<std::size_t> componentId = temperedComponent(systemA, systemB);
  if (!componentId.has_value())
  {
    return sameHamiltonian(systemA.forceField, systemB.forceField);
  }
  // compare everything but the pair tables ('data' and the 1-4 table 'data14') and the charges of the
  // solute's pseudo-atom types: these are exactly the fields ForceField::scaleSoluteInteractions scales
  ForceField normalizedB = systemB.forceField;
  if (normalizedB.data.size() != systemA.forceField.data.size() ||
      normalizedB.data14.size() != systemA.forceField.data14.size() ||
      normalizedB.pseudoAtoms.size() != systemA.forceField.pseudoAtoms.size())
  {
    return false;
  }
  const std::vector<bool> soluteTypes = systemA.pseudoAtomTypesOfComponent(componentId.value());
  const std::size_t n = systemA.forceField.numberOfPseudoAtoms;
  for (std::size_t i = 0; i < n; ++i)
  {
    for (std::size_t j = 0; j < n; ++j)
    {
      if (soluteTypes[i] || soluteTypes[j])
      {
        normalizedB.data[i * n + j] = systemA.forceField.data[i * n + j];
        if (!normalizedB.data14.empty()) normalizedB.data14[i * n + j] = systemA.forceField.data14[i * n + j];
      }
    }
    if (soluteTypes[i]) normalizedB.pseudoAtoms[i].charge = systemA.forceField.pseudoAtoms[i].charge;
  }
  return sameHamiltonian(systemA.forceField, normalizedB);
}

bool sameAtomDefinitionUpToCharge(const Atom& atomA, const Atom& atomB)
{
  return atomA.position == atomB.position && atomA.type == atomB.type;
}

bool compatibleMobileTopology(const System& systemA, const System& systemB)
{
  if (!sameSoluteTempering(systemA, systemB)) return false;
  if (!sameHamiltonianUpToSoluteTempering(systemA, systemB) || systemA.hasExternalField || systemB.hasExternalField ||
      systemA.components.size() != systemB.components.size() ||
      systemA.numberOfFrameworkAtoms != systemB.numberOfFrameworkAtoms ||
      !systemA.reactions.list.empty() || !systemB.reactions.list.empty() ||
      hasUnsupportedFractionalSlots(systemA) || hasUnsupportedFractionalSlots(systemB) ||
      !matchingGCFractionalLayout(systemA, systemB) || !matchingLambdaGrids(systemA, systemB) ||
      systemA.crossLinks.bondTypes.size() != systemB.crossLinks.bondTypes.size())
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

  // Flexible components travel with the configuration like rigid ones: their conformation is in the
  // atom positions, their topology (connectivity, bonded terms, rigid fragments) is in the component,
  // which is the same definition in every replica, and their cross-links are swapped alongside. What
  // stays with the replica is only what the replica derived for its own temperature (growth plans,
  // recoil references, ideal-gas reservoirs), which is bias, not state.
  for (std::size_t componentId = 0; componentId < systemA.components.size(); ++componentId)
  {
    const Component& componentA = systemA.components[componentId];
    const Component& componentB = systemB.components[componentId];
    if (componentA.rigid != componentB.rigid || componentA.name != componentB.name ||
        componentA.atoms.size() != componentB.atoms.size() ||
        componentA.numberOfRigidFragments() != componentB.numberOfRigidFragments())
    {
      return false;
    }
    const std::optional<std::size_t> temperedId = temperedComponent(systemA, systemB);
    const bool tempered = temperedId.has_value() && temperedId.value() == componentId;
    for (std::size_t atomId = 0; atomId < componentA.atoms.size(); ++atomId)
    {
      const bool same = tempered ? sameAtomDefinitionUpToCharge(componentA.atoms[atomId], componentB.atoms[atomId])
                                 : sameAtomDefinition(componentA.atoms[atomId], componentB.atoms[atomId]);
      if (!same)
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

  // the chain is of the tracked component (the adsorbate), not of whatever component comes first
  const std::size_t componentId = replica.tmmcComponent();
  if (componentId >= configuration.components.size()) return 0.0;
  const std::size_t moleculeCount = configuration.numberOfIntegerMoleculesPerComponent[componentId];
  const std::size_t lambdaBin =
      replica.tmmc.lambdaChain() ? configuration.components[componentId].lambdaGC.currentBin : 0uz;
  const std::size_t index = replica.tmmc.chainIndex(moleculeCount, lambdaBin);
  if (index >= replica.tmmc.bias.size()) return std::nullopt;
  return replica.tmmc.bias[index];
}

// Solute tempering: E_target(X) - E_holder(X) for the configuration X held by 'holder', i.e. the energy
// change of switching the solute's pair parameters, charges and intramolecular potentials from the
// Hamiltonian of 'holder' to that of 'target' at fixed positions. Only the solute-involving terms
// differ: the solute's intramolecular energy and its pair interactions with everything (real space,
// tail corrections, framework, Ewald Fourier/self/exclusion). The cost is that of a few
// single-molecule energy evaluations. Empty when an energy routine reports an overlap.
std::optional<double> soluteHamiltonianChange(const System& holder, const System& target);
}  // namespace

std::optional<double> MC_Moves::ParallelTemperingSoluteHamiltonianChange(const System& holder, const System& target)
{
  return soluteHamiltonianChange(holder, target);
}

namespace
{
std::optional<double> soluteHamiltonianChange(const System& holder, const System& target)
{
  const std::optional<std::size_t> temperedId = temperedComponent(holder, target);
  if (!temperedId.has_value() || holder.soluteTemperingLambda == target.soluteTemperingLambda)
  {
    return 0.0;
  }
  const std::size_t componentId = temperedId.value();
  const std::size_t numberOfSoluteMolecules = holder.numberOfMoleculesPerComponent[componentId];
  if (numberOfSoluteMolecules == 0uz) return 0.0;
  const double chargeFactor = std::sqrt(target.soluteTemperingLambda / holder.soluteTemperingLambda);

  // The target Hamiltonian evaluated in the holder's box: its automatic cut-offs and Ewald parameters are box-derived
  // and were derived for the target's own box (NPT replicas differ in volume), so they are derived again for this one
  // (a cut-off beyond half the holder's box broke the minimum image, and a real-space alpha differing from the
  // holder's reciprocal one mixed two Ewald splittings).
  ForceField targetForceField = target.forceField;
  targetForceField.initializeAutomaticCutOff(holder.simulationBox);
  targetForceField.initializeEwaldParameters(holder.simulationBox);

  // the configuration with the solute charges of the target Hamiltonian
  std::span<const Atom> atoms = holder.spanOfMoleculeAtoms();
  std::vector<Atom> atomsTarget(atoms.begin(), atoms.end());
  std::vector<Atom> soluteHolder{};
  std::vector<Atom> soluteTarget{};
  for (Atom& atom : atomsTarget)
  {
    if (static_cast<std::size_t>(atom.componentId) != componentId) continue;
    soluteHolder.push_back(atom);
    atom.charge *= chargeFactor;
    soluteTarget.push_back(atom);
  }

  // pair interactions of the solute molecules with all other molecules; the solute-solute pairs are
  // counted twice in the molecule-by-molecule sum and corrected with the solute-only pair sum
  auto soluteInterEnergy = [&](const ForceField& forceField, std::span<const Atom> all,
                               std::span<const Atom> solute) -> std::optional<double>
  {
    std::optional<RunningEnergy> withAll =
        Interactions::computeInterMolecularEnergyDifference(forceField, holder.simulationBox, all, solute, {});
    if (!withAll.has_value()) return std::nullopt;
    RunningEnergy soluteSolute = Interactions::computeInterMolecularEnergy(forceField, holder.simulationBox, solute);
    return withAll->potentialEnergy() - soluteSolute.potentialEnergy();
  };
  const std::optional<double> interTarget = soluteInterEnergy(targetForceField, atomsTarget, soluteTarget);
  const std::optional<double> interHolder = soluteInterEnergy(holder.forceField, atoms, soluteHolder);
  if (!interTarget.has_value() || !interHolder.has_value()) return std::nullopt;
  double change = interTarget.value() - interHolder.value();

  // tail corrections: the solvent-solvent entries are equal in both force fields and cancel
  change += Interactions::computeInterMolecularTailEnergy(targetForceField, holder.simulationBox, atoms).tail -
            Interactions::computeInterMolecularTailEnergy(holder.forceField, holder.simulationBox, atoms).tail;

  if (holder.framework.has_value())
  {
    std::span<const Atom> frameworkAtoms = holder.spanOfFrameworkAtoms();
    // the grids of a replica tabulate the framework interaction of its own Hamiltonian (System::copyInterpolationGridsFrom)
    std::optional<RunningEnergy> frameworkTarget = Interactions::computeFrameworkMoleculeEnergyDifference(
        targetForceField, holder.simulationBox, target.interpolationGrids, holder.framework, frameworkAtoms,
        soluteTarget, {});
    std::optional<RunningEnergy> frameworkHolder = Interactions::computeFrameworkMoleculeEnergyDifference(
        holder.forceField, holder.simulationBox, holder.interpolationGrids, holder.framework, frameworkAtoms,
        soluteHolder, {});
    if (!frameworkTarget.has_value() || !frameworkHolder.has_value()) return std::nullopt;
    change += frameworkTarget->potentialEnergy() - frameworkHolder->potentialEnergy();
    change +=
        Interactions::computeFrameworkMoleculeTailEnergy(targetForceField, holder.simulationBox, frameworkAtoms, atoms)
            .tail -
        Interactions::computeFrameworkMoleculeTailEnergy(holder.forceField, holder.simulationBox, frameworkAtoms, atoms)
            .tail;
  }

  // Ewald: Fourier, self, net-charge and intramolecular-exclusion terms of the rescaled solute charges
  // (the Ewald parameters are the same in both replicas: same box, same cut-off)
  {
    std::vector<std::complex<double>> eik_x{};
    std::vector<std::complex<double>> eik_y{};
    std::vector<std::complex<double>> eik_z{};
    std::vector<std::complex<double>> eik_xy{};
    EwaldStructureFactor storedEik = holder.storedEik;
    EwaldStructureFactor trialEik{};
    change += Interactions::energyDifferenceEwaldFourier(eik_x, eik_y, eik_z, eik_xy, storedEik, trialEik,
                                                         holder.forceField, holder.simulationBox, holder.components,
                                                         soluteTarget, soluteHolder, holder.ewaldMoveNetCharge(),
                                                         holder.ewaldMoveChargeDerivative(soluteHolder))
                  .potentialEnergy();
  }

  // intramolecular energy of the solute molecules with the scaled potentials of either replica
  std::size_t firstMolecule = 0uz;
  for (std::size_t i = 0; i < componentId; ++i) firstMolecule += holder.numberOfMoleculesPerComponent[i];
  std::span<const Molecule> soluteMolecules{&holder.moleculeData[firstMolecule], numberOfSoluteMolecules};
  change += Interactions::computeIntraMolecularEnergy(targetForceField, holder.simulationBox,
                                                      target.components[componentId].intraMolecularPotentials,
                                                      soluteMolecules, atomsTarget)
                .potentialEnergy() -
            Interactions::computeIntraMolecularEnergy(holder.forceField, holder.simulationBox,
                                                      holder.components[componentId].intraMolecularPotentials,
                                                      soluteMolecules, atoms)
                .potentialEnergy();

  // The cross-links of the configuration (their table travels with it): their bonded terms are not tempered, but the
  // pairs they exclude carry the scaled solute pair parameters and charges, so their correction differs between the
  // two Hamiltonians.
  if (!holder.crossLinks.empty())
  {
    change += Interactions::computeCrossLinkEnergy(targetForceField, holder.simulationBox, target.components,
                                                   holder.numberOfMoleculesPerComponent, atomsTarget, holder.crossLinks)
                  .potentialEnergy() -
              Interactions::computeCrossLinkEnergy(holder.forceField, holder.simulationBox, holder.components,
                                                   holder.numberOfMoleculesPerComponent, atoms, holder.crossLinks)
                  .potentialEnergy();
  }

  return change;
}

// After a swap the replica holds momenta sampled at the temperature of the partner replica:
// rescale them to the replica's own temperature and recompute everything the integrator derives
// from the configuration (gradients, kinetic energies, extended-system energies).
void rebuildMolecularDynamicsState(System& system, double velocityScaling)
{
  Integrators::scaleVelocities(system.moleculeData, system.spanOfMoleculeAtoms(), system.spanOfMoleculeDynamics(),
                               system.components, {velocityScaling, velocityScaling}, system.framework,
                               system.spanOfFrameworkDynamics(), system.spanOfGroupData(),
                               system.spanOfFrameworkGroupData());

  // the degrees of freedom travelled with the configuration (equal for equal molecule counts, but
  // the swap also supports differing counts); the chain state of the heat bath is kept
  if (system.thermostat.has_value())
  {
    system.thermostat->refreshDegreesOfFreedom(system.translationalDegreesOfFreedom,
                                               system.rotationalDegreesOfFreedom,
                                               system.translationalCenterOfMassConstraint);
  }

  system.precomputeTotalGradients();
  Integrators::updateCenterOfMassAndQuaternionGradients(system.moleculeData, system.spanOfMoleculeAtoms(),
                                                        system.spanOfMoleculeDynamics(), system.components,
                                                        system.spanOfGroupData(), system.framework,
                                                        system.spanOfFrameworkDynamics(),
                                                        system.spanOfFrameworkGroupData());
  system.runningEnergies.translationalKineticEnergy = Integrators::computeTranslationalKineticEnergy(
      system.moleculeData, system.spanOfMoleculeAtoms(), system.spanOfMoleculeDynamics(), system.components,
      system.framework, system.spanOfFrameworkAtoms(), system.spanOfFrameworkDynamics(), &system.forceField,
      system.spanOfGroupData(), system.spanOfFrameworkGroupData());
  system.runningEnergies.rotationalKineticEnergy =
      Integrators::computeRotationalKineticEnergy(system.moleculeData, system.components, system.spanOfGroupData(),
                                                  system.framework, system.spanOfFrameworkGroupData());
  if (system.thermostat.has_value())
  {
    system.runningEnergies.NoseHooverEnergy = system.thermostat->getEnergy();
  }
  if (system.thermobarostat.has_value())
  {
    // the box travelled with the configuration (no framework): the barostat's volume coordinate and the
    // molecule-based degrees of freedom it couples to are re-derived from the post-swap state, as the
    // molecular-dynamics driver does when it (re)starts from a configuration
    system.thermobarostat->refreshDegreesOfFreedom(
        barostatTranslationalDegreesOfFreedom(system, system.thermobarostat->coupling),
        system.translationalCenterOfMassConstraint, system.simulationBox.volume);
    system.runningEnergies.thermobarostatEnergy = system.thermobarostat->energy(system.simulationBox.volume);
  }

  // the conserved (extended-system) energy is discontinuous across a swap: restart the drift
  // bookkeeping from the post-swap state
  system.conservedEnergy = system.runningEnergies.conservedEnergy();
  system.referenceEnergy = system.conservedEnergy;
}

}  // namespace

std::optional<std::string> MC_Moves::parallelTemperingUnsupportedReason(const System& system, bool soluteTempering)
{
  if (system.framework.has_value() && system.framework->hasMobileAtoms())
  {
    return std::string(
        "a framework with mobile (flexible) atoms is not supported: the replica-exchange swap "
        "exchanges the adsorbates only, and a replica would keep its own framework configuration "
        "with the energies of its partner");
  }
  for (const VDWParameters& parameter : system.forceField.data)
  {
    if (parameter.type == VDWParameters::Type::FeynmannHibbs)
    {
      return std::string(
          "temperature-dependent pair potentials (Feynman-Hibbs) are not supported: the replicas "
          "then have different Hamiltonians, and the acceptance rule of the swap would need the "
          "energy of each configuration under the potential of the other temperature");
    }
  }
  if (std::ranges::any_of(system.components,
                          [](const Component& component)
                          {
                            return component.idealGasRosenbluthWeight.has_value() &&
                                   component.idealGasRosenbluthWeight.value() != 1.0;
                          }))
  {
    return std::string(
        "an 'IdealGasRosenbluthWeight' is the intramolecular partition function at one temperature and one "
        "Hamiltonian: copied to every replica of a temperature or solute-tempering ladder it is wrong at the "
        "others, and the insertions and the swaps then sample the wrong chemical potentials");
  }
  // the insertions of a flexible molecule draw (or grow) its conformation from the intramolecular Boltzmann
  // distribution, so a replica samples the activity beta f / Z_intra(beta); the swap knows beta f only and the
  // ratio Z_intra(beta_B) / Z_intra(beta_A) per molecule exchanged is missing from its acceptance
  if (std::ranges::any_of(system.components, [](const Component& component)
                          { return component.swappable && !component.rigid && component.atoms.size() > 1uz; }))
  {
    return std::string(
        "an open ensemble with a flexible molecule is not supported: a replica samples the activity beta f over "
        "the intramolecular partition function at its own temperature, which the swap of configurations with "
        "different molecule counts does not know");
  }
  if (soluteTempering && system.forceField.computePolarization)
  {
    return std::string(
        "solute tempering is not supported together with polarization: the change of the solute "
        "Hamiltonian at a swap has no polarization term");
  }
  return std::nullopt;
}

std::optional<double> MC_Moves::ParallelTemperingLogAcceptance(const System& systemA, const System& systemB)
{
  // the cases the drivers refuse up front (see parallelTemperingUnsupportedReason)
  const bool soluteTempering = temperedComponent(systemA, systemB).has_value();
  if (MC_Moves::parallelTemperingUnsupportedReason(systemA, soluteTempering).has_value() ||
      MC_Moves::parallelTemperingUnsupportedReason(systemB, soluteTempering).has_value())
  {
    return std::nullopt;
  }
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

  // Solute tempering: the replicas hold different Hamiltonians H_A, H_B (lambda ladder). With
  // Δ_A(X_B) = E_A(X_B) − E_B(X_B) the energy change of configuration X_B under the Hamiltonian of A,
  //
  //     log R = β_A [E_A(X_A) − E_A(X_B)] + β_B [E_B(X_B) − E_B(X_A)]
  //           = (β_B − β_A)(U_B − U_A) − β_A Δ_A(X_B) − β_B Δ_B(X_A)
  //
  // which reduces to the shared-Hamiltonian term above when the lambdas are equal.
  if (systemA.soluteTemperingLambda != systemB.soluteTemperingLambda)
  {
    const std::optional<double> changeAOfB = soluteHamiltonianChange(systemB, systemA);
    const std::optional<double> changeBOfA = soluteHamiltonianChange(systemA, systemB);
    if (!changeAOfB.has_value() || !changeBOfA.has_value()) return std::nullopt;
    logR -= systemA.beta * changeAOfB.value() + systemB.beta * changeBOfA.value();
  }

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
    // the activity the insertions sample: beta f over the ideal-gas Rosenbluth weight <W_IG>(beta) of the
    // component (its intramolecular partition function), with which the insertions and deletions normalise
    const double activityA = systemA.beta * fugacityA / componentA.idealGasRosenbluthWeight.value_or(1.0);
    const double activityB = systemB.beta * fugacityB / componentB.idealGasRosenbluthWeight.value_or(1.0);
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
  Move::Types move = Move::Types::ParallelTempering;

  systemA.mc_moves_statistics.addTrial(move);

  const std::optional<double> logAcceptance =
      timed(systemA, move, Move::Timing::Fugacity, [&] { return ParallelTemperingLogAcceptance(systemA, systemB); });

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
    systemA.exchangeConfiguration(systemB);

    // the atoms travelled with their charges; with solute tempering the charges belong to the Hamiltonian
    // of the replica (lambda), so the solute charges are rescaled to the receiving replica
    if (systemA.soluteTemperingLambda != systemB.soluteTemperingLambda)
    {
      // a replica that was never scaled is the lambda = 1 member of the ladder: record the component
      if (!systemA.soluteTemperingComponent.has_value()) systemA.soluteTemperingComponent = systemB.soluteTemperingComponent;
      if (!systemB.soluteTemperingComponent.has_value()) systemB.soluteTemperingComponent = systemA.soluteTemperingComponent;
      systemA.rescaleSoluteCharges(std::sqrt(systemA.soluteTemperingLambda / systemB.soluteTemperingLambda));
      systemB.rescaleSoluteCharges(std::sqrt(systemB.soluteTemperingLambda / systemA.soluteTemperingLambda));
    }
    systemA.rebuildConfigurationDerivedState();
    systemB.rebuildConfigurationDerivedState();

    return std::make_pair(systemA.runningEnergies, systemB.runningEnergies);
  }

  return std::nullopt;
}

std::optional<std::pair<RunningEnergy, RunningEnergy>> MC_Moves::ParallelTemperingSwapMolecularDynamics(
    RandomNumber &random, System &systemA, System &systemB)
{
  const double temperatureA = systemA.temperature;
  const double temperatureB = systemB.temperature;

  // the configurational acceptance rule is the same as for Monte Carlo: with the momenta rescaled
  // below, the kinetic parts of the Boltzmann factors cancel exactly
  if (!ParallelTemperingSwap(random, systemA, systemB).has_value())
  {
    return std::nullopt;
  }

  // replica A now holds the momenta generated at T_B (and vice versa)
  rebuildMolecularDynamicsState(systemA, std::sqrt(temperatureA / temperatureB));
  rebuildMolecularDynamicsState(systemB, std::sqrt(temperatureB / temperatureA));

  return std::make_pair(systemA.runningEnergies, systemB.runningEnergies);
}
