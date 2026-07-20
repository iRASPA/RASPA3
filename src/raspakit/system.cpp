module;

module system;

import std;

import archive;
import randomnumbers;
import stringutils;
import int3;
import uint3;
import double3;
import double3x3;
import double3x3x3;
import simd_quatd;
import cubic;
import atom;
import framework;
import component;
import cbmc_move_statistics;
import simulationbox;
import forcefield;
import units;
import property_loading;
import averages;
import skparser;
import skposcarparser;
import skstructure;
import skatom;
import skcell;
import sample_movies;
import property_enthalpy;
import property_pressure;
import energy_dudlambda;
import energy_status;
import energy_status_inter;
import energy_status_intra;
import property_simulationbox;
import average_energy_type;
import property_energy;
import property_partial_molar_properties;
import property_lambda_probability_histogram;
import property_widom;
import property_temperature;
import property_msd;
import running_energy;
import threadpool;
// import isotherm;
// import multi_site_isotherm;
// import pressure_range;
import bond_potential;
import move_statistics;
import mc_moves_probabilities;
import mc_moves_move_types;
import mc_moves_cputime;
import reaction;
import reactions;
import cbmc;
import cbmc_chain_data;
import interactions_framework_molecule;
import interactions_framework_molecule_grid;
import interactions_intermolecular;
import interactions_pair_kernel;
import interactions_ewald;
import interactions_internal;
import interactions_external_field;
import interactions_external_field_grid;
import interactions_polarization_derivatives;
import equation_of_states;
import thermostat;
import thermobarostat;
import json;
import integrators;
import integrators_compute;
import integrators_update;
import interpolation_energy_grid;
import property_number_of_molecules_evolution;
import property_volume_evolution;
import property_conserved_energy_evolution;
import minimization_cell_layout;
#if !(defined(__has_include) && __has_include(<mdspan>))
// import mdspan;
#endif

// construct System programmatically
/*! \brief Brief description.
 *         Brief description continued.
 *
 *  Detailed description starts here.
 */
System::System(ForceField forcefield, std::optional<SimulationBox> box, bool hasExternalField, double T,
               std::optional<double> P, double heliumVoidFraction, std::optional<Framework> f, std::vector<Component> c,
               std::vector<std::vector<double3>> initialpositions, std::vector<std::size_t> initialNumberOfMolecules,
               std::size_t numberOfBlocks, const MCMoveProbabilities& systemProbabilities)
    : temperature(T),
      pressure(P.value_or(0.0) / Units::PressureConversionFactor),
      input_pressure(P.value_or(0.0)),
      beta(1.0 / (Units::KB * T)),
      heliumVoidFraction(heliumVoidFraction),
      framework(f),
      components(c),
      loadings(c.size()),
      swappableComponents(),
      initialNumberOfMolecules(initialNumberOfMolecules),
      numberOfMoleculesPerComponent(c.size()),
      numberOfIntegerMoleculesPerComponent(c.size()),
      numberOfFractionalMoleculesPerComponent(c.size()),
      numberOfGCFractionalMoleculesPerComponent_CFCMC(c.size()),
      numberOfPairGCFractionalMoleculesPerComponent_CFCMC(c.size()),
      numberOfPairSwapFractionalMoleculesPerComponent_CFCMC(c.size()),
      numberOfPairSwapCBFractionalMoleculesPerComponent_CFCMC(c.size()),
      numberOfGibbsSwapFractionalMoleculesPerComponent_CFCMC(c.size()),
      numberOfGibbsFractionalMoleculesPerComponent_CFCMC(c.size()),
      numberOfParallelReactionFractionalMoleculesPerComponent_CFCMC(c.size()),
      numberOfSerialReactionFractionalMoleculesPerComponent_CFCMC(c.size()),
      numberOfReactionFractionalMoleculesPerComponent_CFCMC(),
      idealGasEnergiesPerComponent(c.size()),
      forceField(forcefield),
      hasExternalField(hasExternalField),
      numberOfPseudoAtoms(c.size(), std::vector<std::size_t>(forceField.pseudoAtoms.size())),
      totalNumberOfPseudoAtoms(forceField.pseudoAtoms.size()),
      atomData({}),
      moleculeData({}),
      runningEnergies(),
      currentEnergyStatus(1, f.has_value() ? 1 : 0, c.size()),
      netChargePerComponent(c.size()),
      mc_moves_probabilities(systemProbabilities),
      mc_moves_statistics(),
      reactions(),
      tmmc(),
      averageEnergies(numberOfBlocks, 1, f.has_value() ? 1 : 0, c.size()),
      averageLoadings(numberOfBlocks, c.size()),
      averageEnthalpiesOfAdsorption(numberOfBlocks, c.size()),
      averagePartialMolarProperties(numberOfBlocks, c.size()),
      averageTemperature(numberOfBlocks),
      averageTranslationalTemperature(numberOfBlocks),
      averageRotationalTemperature(numberOfBlocks),
      averagePressure(numberOfBlocks),
      averageSimulationBox(numberOfBlocks),
      interpolationGrids(forceField.pseudoAtoms.size() + 1, std::nullopt)
{
  input_pressureTensorDiagonal = double3(input_pressure, input_pressure, input_pressure);
  pressureTensorDiagonal = double3(pressure, pressure, pressure);

  // Temperature-dependent potentials (Feynman-Hibbs) require the external temperature;
  // recompute the derived constants, shifts, and tail-corrections with the system temperature.
  if (forceField.temperature != T)
  {
    forceField.temperature = T;
    forceField.preComputeDerivedParameters();
    forceField.preComputePotentialShift();
    forceField.preComputeTailCorrection();
  }

  if (box.has_value())
  {
    simulationBox = box.value();
  }

  removeRedundantMoves();
  determineSwappableComponents();
  determineFractionalComponents();
  assignDUdlambdaGroups();
  rescaleMoveProbabilities();
  rescaleMolarFractions();
  computeNumberOfPseudoAtoms();
  computeTailCorrectionCounts();

  createFrameworks();
  if (framework.has_value())
  {
    simulationBox = framework->simulationBox.scaled(framework->numberOfUnitCells);
  }

  forceField.initializeEwaldParameters(simulationBox);

  CoulombicFourierEnergySingleIon = Interactions::computeEwaldFourierEnergySingleIon(
      eik_x, eik_y, eik_z, eik_xy, forceField, simulationBox, double3(0.0, 0.0, 0.0), 1.0);

  precomputeTotalRigidEnergy();

  translationalCenterOfMassConstraint = 0;
  translationalDegreesOfFreedom = 0;
  rotationalDegreesOfFreedom = 0;
  if (framework && framework->hasMobileAtoms())
  {
    if (framework->isMixed())
    {
      translationalDegreesOfFreedom += 3 * framework->flexibleAtomCount + 3 * framework->numberOfRigidGroups();
      for (const FrameworkGroup& group : framework->groups)
      {
        if (group.isRigidBody()) rotationalDegreesOfFreedom += group.rotationalDegreesOfFreedom;
      }
    }
    else
    {
      translationalDegreesOfFreedom += 3 * numberOfFrameworkAtoms;
    }
  }

  createInitialMolecules(initialpositions);
  computeTailCorrectionCounts();

  // Build the per-component ideal-gas conformation reservoirs used to seed CBMC growth. Done after the
  // initial molecules are placed so their placement keeps using the cold-start seed (unchanged initial
  // geometry); the reservoir is only consulted by the production Monte-Carlo moves.
  buildConformationReservoirs();

  equationOfState = EquationOfState(EquationOfState::Type::PengRobinson, EquationOfState::MixingRules::VanDerWaals, T,
                                    P.value_or(0.0), simulationBox, heliumVoidFraction, components);

  averageEnthalpiesOfAdsorption.resize(swappableComponents.size());
  averagePartialMolarProperties.resize(swappableComponents.size());
}

void System::createFrameworks()
{
  netChargeFramework = 0.0;
  if (framework.has_value())
  {
    const std::vector<Atom>& atoms = framework->atoms;
    for (const Atom& atom : atoms)
    {
      atomData.push_back(atom);
      atomDynamics.push_back(AtomDynamics{});
      electricPotential.push_back(0.0);
      electricField.push_back(double3(0.0, 0.0, 0.0));
      electricFieldNew.push_back(double3(0.0, 0.0, 0.0));
    }
    numberOfFrameworkAtoms += atoms.size();
    // Lab-fixed prefix for Ewald: all atoms when fully rigid, Fixed groups when mixed, else none.
    numberOfRigidFrameworkAtoms += framework->numberOfFixedAtoms();
    netChargeFramework += framework->netCharge;
    netCharge += framework->netCharge;
  }
}

void System::rebuildForFramework(const Framework& newFramework, const SimulationBox& newSimulationBox)
{
  const std::size_t previousTranslationalFrameworkDof =
      (framework.has_value() && framework->hasMobileAtoms())
          ? (framework->isMixed() ? 3 * framework->flexibleAtomCount + 3 * framework->numberOfRigidGroups()
                                  : 3 * numberOfFrameworkAtoms)
          : 0;
  std::size_t previousRotationalFrameworkDof = 0;
  if (framework.has_value() && framework->isMixed())
  {
    for (const FrameworkGroup& group : framework->groups)
    {
      if (group.isRigidBody()) previousRotationalFrameworkDof += group.rotationalDegreesOfFreedom;
    }
  }
  const std::size_t previousFrameworkAtoms = numberOfFrameworkAtoms;

  // Detach any guest-molecule storage (the suffix after the framework prefix) so it can be re-appended once
  // the new framework atoms are in place. For a framework-only system these are empty.
  const auto moleculeOffset = static_cast<std::vector<Atom>::difference_type>(previousFrameworkAtoms);
  std::vector<Atom> moleculeAtoms(atomData.begin() + moleculeOffset, atomData.end());
  std::vector<AtomDynamics> moleculeDynamics(atomDynamics.begin() + moleculeOffset, atomDynamics.end());

  framework = newFramework;
  simulationBox = newSimulationBox;

  const std::vector<Atom>& frameworkAtoms = framework->atoms;
  numberOfFrameworkAtoms = frameworkAtoms.size();
  numberOfRigidFrameworkAtoms = framework->numberOfFixedAtoms();

  // Rebuild the per-atom storage as [new framework atoms] ++ [preserved molecule atoms]; the field/potential
  // buffers are reset to the new size (they are recomputed on the next energy/field evaluation).
  atomData.assign(frameworkAtoms.begin(), frameworkAtoms.end());
  atomData.insert(atomData.end(), moleculeAtoms.begin(), moleculeAtoms.end());
  atomDynamics.assign(frameworkAtoms.size(), AtomDynamics{});
  atomDynamics.insert(atomDynamics.end(), moleculeDynamics.begin(), moleculeDynamics.end());
  electricPotential.assign(atomData.size(), 0.0);
  electricField.assign(atomData.size(), double3(0.0, 0.0, 0.0));
  electricFieldNew.assign(atomData.size(), double3(0.0, 0.0, 0.0));

  // Swap the old framework net-charge contribution for the new one (adsorbate contribution is unchanged).
  netCharge -= netChargeFramework;
  netChargeFramework = framework->netCharge;
  netCharge += netChargeFramework;

  // Adjust framework degrees of freedom for the replacement host.
  translationalDegreesOfFreedom -= previousTranslationalFrameworkDof;
  rotationalDegreesOfFreedom -= previousRotationalFrameworkDof;
  if (framework->hasMobileAtoms())
  {
    if (framework->isMixed())
    {
      translationalDegreesOfFreedom += 3 * framework->flexibleAtomCount + 3 * framework->numberOfRigidGroups();
      for (const FrameworkGroup& group : framework->groups)
      {
        if (group.isRigidBody()) rotationalDegreesOfFreedom += group.rotationalDegreesOfFreedom;
      }
    }
    else
    {
      translationalDegreesOfFreedom += 3 * numberOfFrameworkAtoms;
    }
  }

  // Note: pseudo-atom counts intentionally track only component (guest) atoms, matching the constructor which
  // computes them before the framework atoms are appended; the preserved guest counts stay valid, so framework
  // atoms are not (re)counted here.

  forceField.initializeEwaldParameters(simulationBox);
  CoulombicFourierEnergySingleIon = Interactions::computeEwaldFourierEnergySingleIon(
      eik_x, eik_y, eik_z, eik_xy, forceField, simulationBox, double3(0.0, 0.0, 0.0), 1.0);

  precomputeTotalRigidEnergy();
}

void System::determineSwappableComponents()
{
  for (std::size_t componentId{0}; Component& component : components)
  {
    if (component.mc_moves_probabilities.getProbability(Move::Types::Swap) > 0.0 ||
        component.mc_moves_probabilities.getProbability(Move::Types::SwapCBMC) > 0.0 ||
        component.mc_moves_probabilities.getProbability(Move::Types::PairSwapCBMC) > 0.0 ||
        component.mc_moves_probabilities.getProbability(Move::Types::PairSwap) > 0.0 ||
        component.mc_moves_probabilities.getProbability(Move::Types::PairSwapCFCMC) > 0.0 ||
        component.mc_moves_probabilities.getProbability(Move::Types::PairSwapCBCFCMC) > 0.0 ||
        component.mc_moves_probabilities.getProbability(Move::Types::SwapCFCMC) > 0.0 ||
        component.mc_moves_probabilities.getProbability(Move::Types::SwapCBCFCMC) > 0.0)
    {
      component.swappable = true;
    }

    if (component.mc_moves_probabilities.getProbability(Move::Types::GibbsSwapCBMC) > 0.0 ||
        component.mc_moves_probabilities.getProbability(Move::Types::GibbsSwapCFCMC) > 0.0 ||
        component.mc_moves_probabilities.getProbability(Move::Types::GibbsSwapCBCFCMC) > 0.0 ||
        component.mc_moves_probabilities.getProbability(Move::Types::GibbsConventionalCFCMC) > 0.0 ||
        component.mc_moves_probabilities.getProbability(Move::Types::GibbsConventionalCBCFCMC) > 0.0)
    {
      component.swappable = true;
    }

    if (component.swappable)
    {
      swappableComponents.push_back(componentId);
    }

    ++componentId;
  }
}

// determine the required number of fractional molecules

void System::rescaleMoveProbabilities()
{
  for (Component& component : components)
  {
    component.mc_moves_probabilities.join(mc_moves_probabilities);
  }
}

void System::removeRedundantMoves()
{
  for (Component& component : components)
  {
    component.mc_moves_probabilities.removeRedundantMoves();
  }
}

void System::optimizeMCMoves()
{
  mc_moves_statistics.optimizeMCMoves();
  for (Component& component : components)
  {
    component.mc_moves_statistics.optimizeMCMoves();

    // Adapt the internal CBMC / ring-closure Monte-Carlo step sizes (per bead) towards their target
    // acceptance ratios.
    for (CBMCMoveStatistics& cbmcStatistics : component.cbmc_moves_statistics)
    {
      cbmcStatistics.optimize();
    }
  }
}

void System::rescaleMolarFractions()
{
  double totalMolfraction = 0.0;
  double numberOfSwappableComponents = 0.0;
  for (const Component& component : components)
  {
    if (component.swappable)
    {
      totalMolfraction += component.molFraction;
      numberOfSwappableComponents += 1.0;
    }
  }

  if (totalMolfraction > 0.0)
  {
    for (Component& component : components)
    {
      if (component.swappable)
      {
        component.molFraction /= totalMolfraction;
      }
    }
  }
  else
  {
    for (Component& component : components)
    {
      if (component.swappable)
      {
        component.molFraction /= numberOfSwappableComponents;
      }
    }
  }
}

void System::computeNumberOfPseudoAtoms()
{
  for (std::size_t i = 0; i != components.size(); ++i)
  {
    std::fill(numberOfPseudoAtoms[i].begin(), numberOfPseudoAtoms[i].end(), 0);
  }
  std::fill(totalNumberOfPseudoAtoms.begin(), totalNumberOfPseudoAtoms.end(), 0);

  for (const Atom& atom : atomData)
  {
    std::size_t componentId = static_cast<std::size_t>(atom.componentId);
    std::size_t type = static_cast<std::size_t>(atom.type);
    numberOfPseudoAtoms[componentId][type] += 1;
    totalNumberOfPseudoAtoms[type] += 1;
  }
}

void System::computeTailCorrectionCounts()
{
  std::size_t numberOfPseudoAtomTypes = forceField.pseudoAtoms.size();

  effectiveNumberOfPseudoAtomsVDW.assign(numberOfPseudoAtomTypes, 0.0);
  for (std::size_t group = 0; group < maximumNumberOfDUDlambdaGroups; ++group)
  {
    fractionalPseudoAtomCountsPerGroup[group].assign(numberOfPseudoAtomTypes, 0.0);
  }

  for (const Atom& atom : spanOfMoleculeAtoms())
  {
    std::size_t type = static_cast<std::size_t>(atom.type);
    effectiveNumberOfPseudoAtomsVDW[type] += atom.scalingVDW;
    if (atom.groupId != 0)
    {
      fractionalPseudoAtomCountsPerGroup[static_cast<std::size_t>(atom.groupId) - 1][type] += 1.0;
    }
  }
}

void System::addAtomToTailCorrectionCounts(const Atom& atom)
{
  std::size_t type = static_cast<std::size_t>(atom.type);
  effectiveNumberOfPseudoAtomsVDW[type] += atom.scalingVDW;
  if (atom.groupId != 0)
  {
    fractionalPseudoAtomCountsPerGroup[static_cast<std::size_t>(atom.groupId) - 1][type] += 1.0;
  }
}

void System::removeAtomFromTailCorrectionCounts(const Atom& atom)
{
  std::size_t type = static_cast<std::size_t>(atom.type);
  effectiveNumberOfPseudoAtomsVDW[type] -= atom.scalingVDW;
  if (atom.groupId != 0)
  {
    fractionalPseudoAtomCountsPerGroup[static_cast<std::size_t>(atom.groupId) - 1][type] -= 1.0;
  }
}

void System::sampleProperties(std::size_t systemId, std::size_t currentBlock, std::size_t currentCycle)
{
  std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
  double w = weight();

  averageSimulationBox.addSample(currentBlock, simulationBox, w);

  double translationalKineticEnergy = Integrators::computeTranslationalKineticEnergy(
      moleculeData, spanOfMoleculeAtoms(), spanOfMoleculeDynamics(), components, framework, spanOfFrameworkAtoms(),
      spanOfFrameworkDynamics(), &forceField, spanOfGroupData());
  double translationalTemperature =
      2.0 * translationalKineticEnergy /
      (Units::KB * static_cast<double>(translationalDegreesOfFreedom - translationalCenterOfMassConstraint));
  averageTranslationalTemperature.addSample(currentBlock, translationalTemperature, w);

  double rotationalKineticEnergy =
      Integrators::computeRotationalKineticEnergy(moleculeData, components, spanOfGroupData());
  double rotationalTemperature =
      rotationalDegreesOfFreedom > 0
          ? 2.0 * rotationalKineticEnergy / (Units::KB * static_cast<double>(rotationalDegreesOfFreedom))
          : 0.0;
  averageRotationalTemperature.addSample(currentBlock, rotationalTemperature, w);

  double overallTemperature =
      2.0 * (translationalKineticEnergy + rotationalKineticEnergy) /
      (Units::KB * static_cast<double>(translationalDegreesOfFreedom - translationalCenterOfMassConstraint +
                                       rotationalDegreesOfFreedom));
  averageTemperature.addSample(currentBlock, overallTemperature, w);

  loadings = LoadingData(components.size(), numberOfIntegerMoleculesPerComponent, simulationBox);
  averageLoadings.addSample(currentBlock, loadings, w);

  EnthalpyOfAdsorptionTerms enthalpyTerms = EnthalpyOfAdsorptionTerms(
      swappableComponents, numberOfIntegerMoleculesPerComponent, runningEnergies.potentialEnergy(), temperature);
  averageEnthalpiesOfAdsorption.addSample(currentBlock, enthalpyTerms, w);

  PartialMolarPropertiesTerms partialMolarTerms =
      PartialMolarPropertiesTerms(swappableComponents, numberOfIntegerMoleculesPerComponent,
                                  runningEnergies.potentialEnergy(), simulationBox.volume);
  averagePartialMolarProperties.addSample(currentBlock, partialMolarTerms, w);

  std::size_t numberOfMolecules =
      std::accumulate(numberOfIntegerMoleculesPerComponent.begin(), numberOfIntegerMoleculesPerComponent.end(), 0uz);
  double currentIdealPressure = static_cast<double>(numberOfMolecules) / (beta * simulationBox.volume);

  averagePressure.addSample(currentBlock, currentIdealPressure, currentExcessPressureTensor, w);

  for (std::size_t componentId{0}; Component& component : components)
  {
    double componentDensity =
        static_cast<double>(numberOfIntegerMoleculesPerComponent[componentId]) / simulationBox.volume;

    double lambda = component.lambdaGC.lambdaValue();
    double dudlambda = currentDUdlambda(lambda, component.lambdaGC.dUdlambdaGroupId);
    component.lambdaGC.sampleHistogram(currentBlock, componentDensity, dudlambda, containsTheFractionalMolecule, w);

    if (usesGibbsConventionalCFCMC())
    {
      const double gibbsLambda = component.lambdaGibbs.lambdaValue();
      const double gibbsDudlambda = currentDUdlambda(gibbsLambda, component.lambdaGibbs.dUdlambdaGroupId);
      component.lambdaGibbs.sampleHistogram(currentBlock, componentDensity, gibbsDudlambda, true, w);
    }

    if (componentDrivesPairSwapLambda(componentId, Move::Types::PairSwapCFCMC))
    {
      const double pairLambda = component.lambdaPairSwap.lambdaValue();
      component.lambdaPairSwap.sampleHistogram(
          currentBlock, componentDensity, currentDUdlambda(pairLambda, component.lambdaPairSwap.dUdlambdaGroupId),
          containsTheFractionalMolecule, w);
    }

    if (componentDrivesPairSwapLambda(componentId, Move::Types::PairSwapCBCFCMC))
    {
      const double pairLambda = component.lambdaPairSwapCB.lambdaValue();
      component.lambdaPairSwapCB.sampleHistogram(
          currentBlock, componentDensity, currentDUdlambda(pairLambda, component.lambdaPairSwapCB.dUdlambdaGroupId),
          containsTheFractionalMolecule, w);
    }

    ++componentId;
  }

  reactionLambdaSampleProductionHistograms(currentBlock, w);

  updateSamplePDBMovie(systemId, currentCycle);

  if (writeLammpsData.has_value())
  {
    writeLammpsData->update(currentCycle, components, atomData, atomDynamics, moleculeData, simulationBox, forceField,
                            numberOfIntegerMoleculesPerComponent, framework);
  }

  if (propertyConventionalRadialDistributionFunction.has_value())
  {
    propertyConventionalRadialDistributionFunction->sample(simulationBox, spanOfFrameworkAtoms(), spanOfMoleculeAtoms(),
                                                           currentCycle, currentBlock);
  }

  // Force-based RDF is not sampled here. Monte Carlo calls sampleForceBasedRDFWithFullGradients()
  // (full U, including intramolecular). Molecular dynamics reuses integrator forces via
  // sampleForceBasedRDFFromCurrentGradients().

  if (propertyMoleculeProperties.has_value())
  {
    propertyMoleculeProperties->sample(components, numberOfMoleculesPerComponent, spanOfMoleculeAtoms(), currentCycle,
                                       currentBlock);
  }

  if (averageEnergyHistogram.has_value())
  {
    averageEnergyHistogram->addSample(
        currentBlock, currentCycle,
        {runningEnergies.potentialEnergy(), runningEnergies.frameworkMoleculeVDW + runningEnergies.moleculeMoleculeVDW,
         runningEnergies.frameworkMoleculeCharge + runningEnergies.moleculeMoleculeCharge +
             runningEnergies.ewald_fourier + runningEnergies.ewald_self + runningEnergies.ewald_exclusion,
         runningEnergies.polarization},
        w);
  }

  if (averageNumberOfMoleculesHistogram.has_value())
  {
    averageNumberOfMoleculesHistogram->addSample(currentBlock, currentCycle, numberOfIntegerMoleculesPerComponent, w);
  }

  if (propertyMSD.has_value())
  {
    propertyMSD->addSample(currentCycle, moleculeData);
  }

  if (propertyVACF.has_value())
  {
    propertyVACF->addSample(currentCycle, moleculeData);
  }

  if (propertyDensityGrid.has_value())
  {
    propertyDensityGrid->sample(framework, simulationBox, spanOfMoleculeAtoms(), currentCycle);
  }

  std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();

  mc_moves_cputime.propertySampling += (t2 - t1);
}

void System::samplePropertiesEvolution(std::size_t absoluteCurrentCycle)
{
  if (propertyNumberOfMoleculesEvolution.has_value())
  {
    propertyNumberOfMoleculesEvolution->addSample(absoluteCurrentCycle, numberOfIntegerMoleculesPerComponent);
  }
  if (propertyVolumeEvolution.has_value())
  {
    propertyVolumeEvolution->addSample(absoluteCurrentCycle, simulationBox.volume);
  }
  if (propertyConservedEnergyEvolution.has_value())
  {
    propertyConservedEnergyEvolution->addSample(absoluteCurrentCycle, runningEnergies);
  }
}

void System::precomputeTotalRigidEnergy() noexcept
{
  Interactions::precomputeEwaldFourierRigid(eik_x, eik_y, eik_z, eik_xy, fixedFrameworkStoredEik, forceField,
                                            simulationBox, spanOfRigidFrameworkAtoms());
}

void System::precomputeTotalGradients() noexcept
{
  runningEnergies = Integrators::updateGradients(
      moleculeData, spanOfMoleculeAtoms(), spanOfMoleculeDynamics(), spanOfFrameworkAtoms(), forceField, simulationBox,
      components, eik_x, eik_y, eik_z, eik_xy, trialEik, fixedFrameworkStoredEik, interpolationGrids,
      numberOfMoleculesPerComponent, framework, spanOfFrameworkDynamics());
}

RunningEnergy System::computeTotalEnergies() noexcept
{
  if (fixedFrameworkStoredEik.empty())
  {
    precomputeTotalRigidEnergy();
  }

  std::span<const Atom> frameworkAtomPositions = spanOfFrameworkAtoms();
  std::span<Atom> moleculeAtomPositions = spanOfMoleculeAtoms();

  RunningEnergy runningIntraEnergy{};
  std::size_t index = 0;
  for (std::size_t i = 0; i < components.size(); ++i)
  {
    if (numberOfMoleculesPerComponent[i] > 0)
    {
      std::span<const Molecule> span_molecules = {&moleculeData[index], numberOfMoleculesPerComponent[i]};
      runningIntraEnergy += Interactions::computeIntraMolecularEnergy(components[i].intraMolecularPotentials,
                                                                      span_molecules, spanOfMoleculeAtoms());
    }

    index += numberOfMoleculesPerComponent[i];
  }
  if (framework && framework->hasMobileAtoms())
  {
    runningIntraEnergy += Interactions::computeFrameworkIntraMolecularEnergy(forceField, *framework, simulationBox,
                                                                             frameworkAtomPositions);
  }

  if (forceField.computePolarization)
  {
    std::span<double3> moleculeElectricField = spanOfMoleculeElectricField();

    std::fill(moleculeElectricField.begin(), moleculeElectricField.end(), double3(0.0, 0.0, 0.0));

    RunningEnergy frameworkMoleculeEnergy = Interactions::computeFrameworkMoleculeElectricField(
        forceField, simulationBox, moleculeElectricField, frameworkAtomPositions, moleculeAtomPositions);

    // When molecule-molecule polarization is requested, the inter-molecular electric field has to be added
    // to the stored field so that the polarization energy is consistent with the incremental Monte-Carlo moves
    // (which maintain the same field). Otherwise only the framework field contributes and the plain
    // inter-molecular energy is sufficient.
    RunningEnergy intermolecularEnergy;
    if (forceField.omitInterPolarization)
    {
      intermolecularEnergy =
          Interactions::computeInterMolecularEnergy(forceField, simulationBox, moleculeAtomPositions);
    }
    else
    {
      intermolecularEnergy = Interactions::computeInterMolecularElectricField(
          forceField, simulationBox, moleculeElectricField, moleculeAtomPositions);
    }

    RunningEnergy frameworkMoleculeTailEnergy = Interactions::computeFrameworkMoleculeTailEnergy(
        forceField, simulationBox, frameworkAtomPositions, moleculeAtomPositions);
    RunningEnergy intermolecularTailEnergy =
        Interactions::computeInterMolecularTailEnergy(forceField, simulationBox, moleculeAtomPositions);

    RunningEnergy ewaldEnergy = Interactions::computeEwaldFourierElectricField(
        eik_x, eik_y, eik_z, eik_xy, fixedFrameworkStoredEik, storedEik, forceField, simulationBox,
        moleculeElectricField, components, numberOfMoleculesPerComponent, moleculeAtomPositions);

    RunningEnergy polarizationEnergy = computePolarizationEnergy();

    RunningEnergy externalFieldEnergy;
    Interactions::computeExternalFieldEnergy(hasExternalField, forceField, simulationBox, moleculeAtomPositions,
                                             externalFieldEnergy, externalFieldInterpolationGrid);

    return frameworkMoleculeEnergy + intermolecularEnergy + frameworkMoleculeTailEnergy + intermolecularTailEnergy +
           ewaldEnergy + polarizationEnergy + runningIntraEnergy + externalFieldEnergy;
  }
  else
  {
    RunningEnergy frameworkMoleculeEnergy = Interactions::computeFrameworkMoleculeEnergy(
        forceField, simulationBox, interpolationGrids, framework, frameworkAtomPositions, moleculeAtomPositions);
    RunningEnergy intermolecularEnergy =
        Interactions::computeInterMolecularEnergy(forceField, simulationBox, moleculeAtomPositions);

    RunningEnergy frameworkMoleculeTailEnergy = Interactions::computeFrameworkMoleculeTailEnergy(
        forceField, simulationBox, frameworkAtomPositions, moleculeAtomPositions);
    RunningEnergy intermolecularTailEnergy =
        Interactions::computeInterMolecularTailEnergy(forceField, simulationBox, moleculeAtomPositions);

    RunningEnergy ewaldEnergy = Interactions::computeEwaldFourierEnergy(
        eik_x, eik_y, eik_z, eik_xy, fixedFrameworkStoredEik, storedEik, forceField, simulationBox, components,
        numberOfMoleculesPerComponent, moleculeAtomPositions, netChargeFramework);

    RunningEnergy externalFieldEnergy;
    Interactions::computeExternalFieldEnergy(hasExternalField, forceField, simulationBox, moleculeAtomPositions,
                                             externalFieldEnergy, externalFieldInterpolationGrid);

    return frameworkMoleculeEnergy + intermolecularEnergy + frameworkMoleculeTailEnergy + intermolecularTailEnergy +
           ewaldEnergy + runningIntraEnergy + externalFieldEnergy;
  }
}

RunningEnergy System::computePolarizationEnergy() noexcept
{
  RunningEnergy energy{};

  std::span<const Atom> moleculeAtomPositions = spanOfMoleculeAtoms();
  std::span<double3> moleculeElectricField = spanOfMoleculeElectricField();

  for (std::size_t i = 0; i < moleculeAtomPositions.size(); ++i)
  {
    std::size_t type = moleculeAtomPositions[i].type;
    // The polarization coupling is scaled by the atom's Coulomb scaling so that a fractional (CFCMC)
    // molecule decouples from the field as lambda decreases (matching the incremental moves).
    double polarizability = moleculeAtomPositions[i].scalingCoulomb *
                            forceField.pseudoAtoms[type].polarizability / Units::CoulombicConversionFactor;
    energy.polarization -= 0.5 * polarizability * double3::dot(moleculeElectricField[i], moleculeElectricField[i]);
  }

  return energy;
}

void System::computeTotalElectrostaticPotential() noexcept
{
  if (fixedFrameworkStoredEik.empty())
  {
    precomputeTotalRigidEnergy();
  }

  std::span<Atom> frameworkAtomPositions = spanOfFrameworkAtoms();
  std::span<Atom> moleculeAtomPositions = spanOfMoleculeAtoms();
  std::span<double> moleculeElectrostaticPotential = spanOfMoleculeElectrostaticPotential();

  std::fill(moleculeElectrostaticPotential.begin(), moleculeElectrostaticPotential.end(), 0.0);

  Interactions::computeInterMolecularElectrostaticPotential(forceField, simulationBox, moleculeElectrostaticPotential,
                                                            moleculeAtomPositions);

  Interactions::computeFrameworkMoleculeElectrostaticPotential(
      forceField, simulationBox, moleculeElectrostaticPotential, frameworkAtomPositions, moleculeAtomPositions);

  Interactions::computeEwaldFourierElectrostaticPotential(
      eik_x, eik_y, eik_z, eik_xy, fixedFrameworkStoredEik, storedEik, moleculeElectrostaticPotential, forceField,
      simulationBox, components, numberOfMoleculesPerComponent, moleculeAtomPositions);
}

void System::computeTotalElectricField() noexcept
{
  if (fixedFrameworkStoredEik.empty())
  {
    precomputeTotalRigidEnergy();
  }

  std::span<Atom> frameworkAtomPositions = spanOfFrameworkAtoms();
  std::span<Atom> moleculeAtomPositions = spanOfMoleculeAtoms();
  std::span<double3> moleculeElectricField = spanOfMoleculeElectricField();

  std::fill(moleculeElectricField.begin(), moleculeElectricField.end(), double3(0.0, 0.0, 0.0));

  Interactions::computeInterMolecularElectricField(forceField, simulationBox, moleculeElectricField,
                                                   moleculeAtomPositions);

  Interactions::computeFrameworkMoleculeElectricField(forceField, simulationBox, moleculeElectricField,
                                                      frameworkAtomPositions, moleculeAtomPositions);

  Interactions::computeEwaldFourierElectricField(eik_x, eik_y, eik_z, eik_xy, fixedFrameworkStoredEik, storedEik,
                                                 forceField, simulationBox, moleculeElectricField, components,
                                                 numberOfMoleculesPerComponent, moleculeAtomPositions);
}

std::pair<EnergyStatus, double3x3> System::computeMolecularPressure() noexcept
{
  // Scratch buffer so molecular-pressure sampling does not mutate live MD site gradients.
  // Strain-derivative routines write intermolecular/framework forces here for the atomic-to-molecular
  // virial correction. Intramolecular bonded forces are intentionally omitted (they cancel in the
  // molecular virial) and must not overwrite AtomDynamics used by Velocity-Verlet.
  std::vector<AtomDynamics> pressureMoleculeDynamics(spanOfMoleculeAtoms().size());
  std::span<AtomDynamics> pressureDynamics(pressureMoleculeDynamics);

  const std::span<const Atom> moleculeAtoms = spanOfMoleculeAtoms();

  // Polarization rides along with the real-space Coulomb strain loops: they gather the per-atom electric
  // field and its cell-strain response (molecular COM scaling) in one pass, so the pair walk is not
  // repeated. The Ewald reciprocal framework field and the final contraction into the polarization energy
  // and strain tensor happen afterwards in computePolarizationMolecularPressureStrain. The polarization
  // forces are intentionally kept out of 'pressureDynamics': the polarization strain below already uses
  // COM arms, so it must not participate in the atomic-to-molecular virial correction.
  const bool gatherPolarization = forceField.computePolarization && forceField.useCharge && !moleculeAtoms.empty();
  std::vector<double3> polarizationField;
  std::vector<std::array<double3, 9>> polarizationFieldStrain;
  std::vector<double3> polarizationComOffset;
  std::vector<double> polarizationPolarizability;
  Interactions::PolarizationFieldStrain polarizationGather{};
  const Interactions::PolarizationFieldStrain* frameworkGather = nullptr;
  const Interactions::PolarizationFieldStrain* interGather = nullptr;
  if (gatherPolarization)
  {
    polarizationField.assign(moleculeAtoms.size(), double3(0.0, 0.0, 0.0));
    polarizationFieldStrain.assign(moleculeAtoms.size(), {});
    polarizationComOffset.assign(moleculeAtoms.size(), double3(0.0, 0.0, 0.0));
    polarizationPolarizability.resize(moleculeAtoms.size());
    for (std::size_t i = 0; i < moleculeAtoms.size(); ++i)
    {
      // Scaled by the atom's Coulomb scaling: fractional (CFCMC) molecules decouple from the field.
      polarizationPolarizability[i] =
          moleculeAtoms[i].scalingCoulomb *
          forceField.pseudoAtoms[static_cast<std::size_t>(moleculeAtoms[i].type)].polarizability /
          Units::CoulombicConversionFactor;
    }

    // Mass-weighted COM offsets from the current atom positions for every molecule (rigid and flexible),
    // matching the COM-scaling volume move.
    for (const Molecule& molecule : moleculeData)
    {
      double totalMass = 0.0;
      double3 com(0.0, 0.0, 0.0);
      for (std::size_t k = 0; k < molecule.numberOfAtoms; ++k)
      {
        const Atom& atom = moleculeAtoms[molecule.atomIndex + k];
        const double mass = forceField.pseudoAtoms[static_cast<std::size_t>(atom.type)].mass;
        com += mass * atom.position;
        totalMass += mass;
      }
      com = com / totalMass;
      for (std::size_t k = 0; k < molecule.numberOfAtoms; ++k)
      {
        polarizationComOffset[molecule.atomIndex + k] = moleculeAtoms[molecule.atomIndex + k].position - com;
      }
    }

    polarizationGather = {std::span<double3>(polarizationField),
                          std::span<std::array<double3, 9>>(polarizationFieldStrain),
                          std::span<const double3>(polarizationComOffset),
                          std::span<const double>(polarizationPolarizability)};
    frameworkGather = &polarizationGather;
    // The inter-molecular field obeys the omit flags of the polarization model; the routine's own
    // omitInterInteractions early-out covers the remaining case.
    if (!forceField.omitInterPolarization) interGather = &polarizationGather;
  }

  std::pair<EnergyStatus, double3x3> pressureInfo = Interactions::computeFrameworkMoleculeEnergyStrainDerivative(
      forceField, framework, interpolationGrids, components, simulationBox, spanOfFrameworkAtoms(),
      spanOfMoleculeAtoms(), pressureDynamics, frameworkGather);

  pressureInfo.first.translationalKineticEnergy = runningEnergies.translationalKineticEnergy;
  pressureInfo.first.rotationalKineticEnergy = runningEnergies.rotationalKineticEnergy;
  pressureInfo.first.noseHooverEnergy = runningEnergies.NoseHooverEnergy;

  pressureInfo = pairSum(pressureInfo,
                         Interactions::computeInterMolecularEnergyStrainDerivative(
                             forceField, components, simulationBox, spanOfMoleculeAtoms(), pressureDynamics,
                             interGather));

  pressureInfo = pairSum(pressureInfo,
                         Interactions::computeEwaldFourierEnergyStrainDerivative(
                             eik_x, eik_y, eik_z, eik_xy, fixedFrameworkStoredEik, storedEik, forceField, simulationBox,
                             framework, components, numberOfMoleculesPerComponent, spanOfMoleculeAtoms(),
                             pressureDynamics, netChargeFramework, netChargePerComponent));

  std::size_t molecule_index = 0;
  for (std::size_t i = 0; i < components.size(); ++i)
  {
    if (numberOfMoleculesPerComponent[i] > 0)
    {
      std::span<const Molecule> span_molecules = {&moleculeData[molecule_index], numberOfMoleculesPerComponent[i]};
      RunningEnergy runningIntraEnergy = Interactions::computeIntraMolecularEnergy(
          components[i].intraMolecularPotentials, span_molecules, spanOfMoleculeAtoms());

      pressureInfo.first.intraComponentEnergies[i].bond += runningIntraEnergy.bond;
      pressureInfo.first.intraComponentEnergies[i].ureyBradley += runningIntraEnergy.ureyBradley;
      pressureInfo.first.intraComponentEnergies[i].bend += runningIntraEnergy.bend;
      pressureInfo.first.intraComponentEnergies[i].inversionBend += runningIntraEnergy.inversionBend;
      pressureInfo.first.intraComponentEnergies[i].outOfPlaneBend += runningIntraEnergy.outOfPlaneBend;
      pressureInfo.first.intraComponentEnergies[i].torsion += runningIntraEnergy.torsion;
      pressureInfo.first.intraComponentEnergies[i].improperTorsion += runningIntraEnergy.improperTorsion;
      pressureInfo.first.intraComponentEnergies[i].bondBond += runningIntraEnergy.bondBond;
      pressureInfo.first.intraComponentEnergies[i].bondBend += runningIntraEnergy.bondBend;
      pressureInfo.first.intraComponentEnergies[i].bondTorsion += runningIntraEnergy.bondTorsion;
      pressureInfo.first.intraComponentEnergies[i].bendBend += runningIntraEnergy.bendBend;
      pressureInfo.first.intraComponentEnergies[i].bendTorsion += runningIntraEnergy.bendTorsion;
      pressureInfo.first.intraComponentEnergies[i].vanDerWaals += runningIntraEnergy.intraVDW;
      pressureInfo.first.intraComponentEnergies[i].coulomb += runningIntraEnergy.intraCoul;

      // Intramolecular potentials (bonds, bends, torsions, ...) must NOT contribute to the molecular
      // (center-of-mass based) pressure: internal forces sum to zero over each molecule and cancel in the
      // molecular virial. Their strain derivative is therefore intentionally not accumulated here, and their
      // gradients are deliberately kept out of 'pressureDynamics' so they cannot pollute the
      // atomic-to-molecular correction term computed below. Only the energy is recorded (above) for reporting.
    }

    molecule_index += numberOfMoleculesPerComponent[i];
  }

  if (gatherPolarization)
  {
    // Complete the gathered real-space field with the Ewald reciprocal framework contribution and contract
    // into the polarization energy and its (unsymmetrized) strain-derivative tensor; the symmetrization at
    // the end of this function averages the off-diagonal pairs, matching the symmetric strain generators.
    const auto [polarizationEnergy, polarizationStrain] = Interactions::computePolarizationMolecularPressureStrain(
        *this, polarizationField, polarizationFieldStrain, polarizationComOffset, polarizationPolarizability);

    pressureInfo.first.polarizationEnergy = EnergyDuDlambda(polarizationEnergy, 0.0);
    pressureInfo.second += polarizationStrain;
  }

  pressureInfo.first.sumTotal();

  double pressureTailCorrection = 0.0;
  // Tail correction to the (excess) pressure virial. The per-pair integral 'tailCorrectionPressure' equals
  // Integrate[U'(r) r^3, {r, rc, Inf}]; the isotropic virial tail correction is -(2 pi / 3 V) Sum_ij n_i n_j <integral>
  // (see RASPA2 CalculateTailCorrection). The overall minus sign and the 1/3 are essential: the van der Waals tail is
  // attractive and must LOWER the pressure. The extra global negation of 'pressureInfo.second' below turns the
  // '-= pressureTailCorrection' into the correct negative diagonal contribution.
  double preFactor = -2.0 * std::numbers::pi / (3.0 * simulationBox.volume);
  for (std::vector<Atom>::iterator it1 = atomData.begin(); it1 != atomData.end(); ++it1)
  {
    std::size_t typeA = static_cast<std::size_t>(it1->type);
    double scalingVDWA = it1->scalingVDW;

    pressureTailCorrection += scalingVDWA * scalingVDWA * preFactor * forceField(typeA, typeA).tailCorrectionPressure;

    for (std::vector<Atom>::iterator it2 = it1 + 1; it2 != atomData.end(); ++it2)
    {
      std::size_t typeB = static_cast<std::size_t>(it2->type);
      double scalingVDWB = it2->scalingVDW;

      pressureTailCorrection +=
          scalingVDWA * scalingVDWB * 2.0 * preFactor * forceField(typeA, typeB).tailCorrectionPressure;
    }
  }

  pressureInfo.second.ax -= pressureTailCorrection;
  pressureInfo.second.by -= pressureTailCorrection;
  pressureInfo.second.cz -= pressureTailCorrection;

  // Correct rigid molecule contribution using the constraints forces.
  // Molecule::atomIndex indexes the molecule-atom span (framework atoms excluded).
  double3x3 correctionTerm{};
  for (Molecule& molecule : moleculeData)
  {
    const std::span<const Atom> span = moleculeAtoms.subspan(molecule.atomIndex, molecule.numberOfAtoms);
    const std::span<AtomDynamics> spanDynamics = {&pressureMoleculeDynamics[molecule.atomIndex],
                                                  molecule.numberOfAtoms};

    double totalMass = 0.0;
    double3 com(0.0, 0.0, 0.0);
    for (const Atom& atom : span)
    {
      double mass = forceField.pseudoAtoms[static_cast<std::size_t>(atom.type)].mass;
      com += mass * atom.position;
      totalMass += mass;
    }
    com = com / totalMass;

    for (std::size_t k = 0; k < span.size(); ++k)
    {
      const double3 position = span[k].position;
      const double3 gradient = spanDynamics[k].gradient;

      correctionTerm.ax += (position.x - com.x) * gradient.x;
      correctionTerm.ay += (position.x - com.x) * gradient.y;
      correctionTerm.az += (position.x - com.x) * gradient.z;

      correctionTerm.bx += (position.y - com.y) * gradient.x;
      correctionTerm.by += (position.y - com.y) * gradient.y;
      correctionTerm.bz += (position.y - com.y) * gradient.z;

      correctionTerm.cx += (position.z - com.z) * gradient.x;
      correctionTerm.cy += (position.z - com.z) * gradient.y;
      correctionTerm.cz += (position.z - com.z) * gradient.z;
    }
  }

  pressureInfo.second = -(pressureInfo.second - correctionTerm);

  double temp = 0.5 * (pressureInfo.second.ay + pressureInfo.second.bx);
  pressureInfo.second.ay = pressureInfo.second.bx = temp;
  temp = 0.5 * (pressureInfo.second.az + pressureInfo.second.cx);
  pressureInfo.second.az = pressureInfo.second.cx = temp;
  temp = 0.5 * (pressureInfo.second.bz + pressureInfo.second.cy);
  pressureInfo.second.bz = pressureInfo.second.cy = temp;

  return pressureInfo;
}

void System::checkCartesianPositions()
{
  std::span<Atom> moleculeAtomPositions = spanOfMoleculeAtoms();

  std::size_t index{};
  for (Molecule& molecule : moleculeData)
  {
    std::span<Atom> span = std::span(&moleculeAtomPositions[index], molecule.numberOfAtoms);
    if (components[molecule.componentId].rigid)
    {
      simd_quatd q = molecule.orientation;
      double3x3 M = double3x3::buildRotationMatrixInverse(q);

      for (std::size_t i = 0; i != span.size(); i++)
      {
        double3 expandedPosition =
            molecule.centerOfMassPosition + M * components[molecule.componentId].atoms[i].position;
        if ((std::abs(span[i].position.x - expandedPosition.x) > 1e-5) ||
            (std::abs(span[i].position.y - expandedPosition.y) > 1e-5) ||
            (std::abs(span[i].position.z - expandedPosition.z) > 1e-5))
        {
          throw std::runtime_error(
              std::format("Difference detected between atom position ({} {} {}) and position generated from "
                          "quaternion ({} {} {})\n",
                          span[i].position.x, span[i].position.y, span[i].position.z, expandedPosition.x,
                          expandedPosition.y, expandedPosition.z));
        }
      }
    }
    index += molecule.numberOfAtoms;
  }
}

void System::setThermostat(const std::optional<Thermostat>& thermo)
{
  if (thermo.has_value())
  {
    thermostat = Thermostat(temperature, timeStep, translationalDegreesOfFreedom, rotationalDegreesOfFreedom,
                            thermo->thermostatChainLength, thermo->numberOfYoshidaSuzukiSteps,
                            thermo->timeScaleParameterThermostat);
  }
}

void System::setThermobarostat(const std::optional<Thermobarostat>& barostat)
{
  if (!barostat.has_value())
  {
    thermobarostat.reset();
    return;
  }
  molecularDynamicsEnsemble = barostat->ensemble;
  thermobarostat = Thermobarostat(barostat->ensemble, barostat->cellType, barostat->monoclinicAngle, temperature,
                                  pressure, timeStep, translationalDegreesOfFreedom, barostat->chainLength,
                                  barostat->numberOfYoshidaSuzukiSteps, barostat->timeScaleParameterBarostat);
  thermobarostat->numberOfRespaSteps = barostat->numberOfRespaSteps;
}

void System::setSamplePDBMovie(const std::optional<SampleMovie>& movie)
{
  if (movie.has_value())
  {
    samplePDBMovie = movie;
  }
}

void System::updateSamplePDBMovie(std::size_t systemId, std::size_t currentCycle)
{
  if (samplePDBMovie.has_value())
  {
    samplePDBMovie->update(forceField, systemId, simulationBox, spanOfMoleculeAtoms(), components,
                           numberOfMoleculesPerComponent, currentCycle, spanOfFrameworkAtoms());
  }
}

void System::setNumberOfMoleculesHistogram(const std::optional<PropertyNumberOfMoleculesHistogram>& hist)
{
  if (hist.has_value())
  {
    averageNumberOfMoleculesHistogram = PropertyNumberOfMoleculesHistogram(
        hist->numberOfBlocks, components.size(), hist->range, hist->sampleEvery, hist->writeEvery);
  }
}

void System::setAverageEnergyHistogram(const std::optional<PropertyEnergyHistogram>& hist)
{
  if (hist.has_value())
  {
    averageEnergyHistogram = hist;
  }
}

void System::setPropertyDensityGrid(const std::optional<PropertyDensityGrid>& grid)
{
  if (grid.has_value())
  {
    propertyDensityGrid = grid;
  }
}

void System::setPropertyNumberOfMoleculesEvolution(std::optional<PropertyNumberOfMoleculesEvolution> property)
{
  if (property.has_value())
  {
    propertyNumberOfMoleculesEvolution = property;
  }
}

void System::setPropertyVolumeEvolution(std::optional<PropertyVolumeEvolution> property)
{
  if (property.has_value())
  {
    propertyVolumeEvolution = property;
  }
}

void System::setPropertyConservedEnergyEvolution(std::optional<PropertyConservedEnergyEvolution> property)
{
  if (property.has_value())
  {
    propertyConservedEnergyEvolution = property;
  }
}

void System::setPropertyConventionalRDF(const std::optional<PropertyConventionalRadialDistributionFunction>& rdf)
{
  if (rdf.has_value())
  {
    propertyConventionalRadialDistributionFunction = PropertyConventionalRadialDistributionFunction(
        5, forceField.pseudoAtoms.size(), rdf->numberOfBins, 12.0, rdf->sampleEvery, rdf->writeEvery);
  }
}

void System::setPropertyRDF(const std::optional<PropertyRadialDistributionFunction>& rdf)
{
  if (rdf.has_value())
  {
    propertyRadialDistributionFunction = PropertyRadialDistributionFunction(
        5, forceField.pseudoAtoms.size(), rdf->numberOfBins, rdf->range, rdf->sampleEvery, rdf->writeEvery);
  }
}

bool System::forceBasedRDFSampleDue(std::size_t currentCycle) const
{
  return propertyRadialDistributionFunction.has_value() && propertyRadialDistributionFunction->sampleEvery > 0uz &&
         (currentCycle % propertyRadialDistributionFunction->sampleEvery == 0uz);
}

void System::sampleForceBasedRDFFromCurrentGradients(std::size_t currentCycle, std::size_t currentBlock)
{
  if (!propertyRadialDistributionFunction.has_value()) return;
  propertyRadialDistributionFunction->sample(simulationBox, spanOfFrameworkAtoms(), spanOfFrameworkDynamics(),
                                             moleculeData, spanOfMoleculeAtoms(), spanOfMoleculeDynamics(),
                                             currentCycle, currentBlock);
}

void System::sampleForceBasedRDFWithFullGradients(std::size_t currentCycle, std::size_t currentBlock)
{
  if (!forceBasedRDFSampleDue(currentCycle)) return;
  precomputeTotalGradients();
  sampleForceBasedRDFFromCurrentGradients(currentCycle, currentBlock);
}

void System::setPropertyMSD(const std::optional<PropertyMeanSquaredDisplacement>& msd)
{
  if (msd.has_value())
  {
    propertyMSD = PropertyMeanSquaredDisplacement(numberOfMoleculesPerComponent, moleculeData.size(), timeStep,
                                                  msd->numberOfBlockElementsMSD, msd->sampleEvery, msd->writeEvery);
  }
}

void System::setPropertyVACF(const std::optional<PropertyVelocityAutoCorrelationFunction>& vacf)
{
  if (vacf.has_value())
  {
    propertyVACF = PropertyVelocityAutoCorrelationFunction(numberOfMoleculesPerComponent, moleculeData.size(), timeStep,
                                                           vacf->numberOfBuffersVACF, vacf->bufferLengthVACF,
                                                           vacf->sampleEvery, vacf->writeEvery);
  }
}

Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const System& s)
{
  archive << s.versionNumber;

  archive << s.temperature;
  archive << s.pressure;
  archive << s.input_pressure;
  archive << s.pressureTensorDiagonal;
  archive << s.input_pressureTensorDiagonal;
  archive << s.beta;
  archive << static_cast<std::uint8_t>(s.cellMinimizationType);
  archive << static_cast<std::uint8_t>(s.monoclinicAngleType);

  archive << s.heliumVoidFraction;

  archive << s.numberOfFrameworks;
  archive << s.numberOfFrameworkAtoms;
  archive << s.numberOfRigidFrameworkAtoms;

  archive << s.framework;
  archive << s.components;

  archive << s.equationOfState;
  archive << static_cast<std::uint8_t>(s.molecularDynamicsEnsemble);
  archive << s.thermostat;
  archive << s.thermobarostat;

  archive << s.loadings;

  archive << s.swappableComponents;
  archive << s.initialNumberOfMolecules;

  archive << s.numberOfMoleculesPerComponent;
  archive << s.numberOfIntegerMoleculesPerComponent;
  archive << s.numberOfFractionalMoleculesPerComponent;
  archive << s.numberOfGCFractionalMoleculesPerComponent_CFCMC;
  archive << s.numberOfPairGCFractionalMoleculesPerComponent_CFCMC;
  archive << s.numberOfPairSwapFractionalMoleculesPerComponent_CFCMC;
  archive << s.numberOfPairSwapCBFractionalMoleculesPerComponent_CFCMC;
  archive << s.numberOfGibbsSwapFractionalMoleculesPerComponent_CFCMC;
  archive << s.numberOfGibbsFractionalMoleculesPerComponent_CFCMC;
  archive << s.numberOfParallelReactionFractionalMoleculesPerComponent_CFCMC;
  archive << s.numberOfSerialReactionFractionalMoleculesPerComponent_CFCMC;
  archive << s.numberOfReactionFractionalMoleculesPerComponent_CFCMC;

  archive << s.idealGasEnergiesPerComponent;

  archive << s.forceField;
  archive << s.hasExternalField;

  archive << s.numberOfPseudoAtoms;
  archive << s.totalNumberOfPseudoAtoms;

  archive << s.translationalCenterOfMassConstraint;
  archive << s.translationalDegreesOfFreedom;
  archive << s.rotationalDegreesOfFreedom;

  archive << s.timeStep;

  archive << s.simulationBox;
  archive << s.containsTheFractionalMolecule;

  archive << s.atomData;
  archive << s.atomDynamics;
  archive << s.moleculeData;
  archive << s.electricPotential;
  archive << s.electricField;
  archive << s.electricFieldNew;

  archive << s.conservedEnergy;
  archive << s.referenceEnergy;
  archive << s.accumulatedDrift;

  archive << s.rigidEnergies;
  archive << s.runningEnergies;

  archive << s.currentExcessPressureTensor;
  archive << s.currentEnergyStatus;

  archive << s.numberOfHybridMCSteps;

  archive << s.eik_xy;
  archive << s.eik_x;
  archive << s.eik_y;
  archive << s.eik_z;
  archive << s.storedEik;
  archive << s.fixedFrameworkStoredEik;
  archive << s.trialEik;
  archive << s.CoulombicFourierEnergySingleIon;
  archive << s.netCharge;
  archive << s.netChargeFramework;
  archive << s.netChargeAdsorbates;
  archive << s.netChargePerComponent;

  archive << s.mc_moves_probabilities;
  archive << s.mc_moves_statistics;
  archive << s.mc_moves_cputime;

  archive << s.reactions;
  archive << s.tmmc;

  archive << s.averageEnergies;
  archive << s.averageLoadings;
  archive << s.averageEnthalpiesOfAdsorption;
  archive << s.averagePartialMolarProperties;
  archive << s.averageTemperature;
  archive << s.averageTranslationalTemperature;
  archive << s.averageRotationalTemperature;
  archive << s.averagePressure;
  archive << s.averageSimulationBox;
  archive << s.propertyElasticConstantsFluctuation;
  archive << s.elasticConstantsSampleEvery;

  archive << s.samplePDBMovie;

  archive << s.propertyConventionalRadialDistributionFunction;
  archive << s.propertyRadialDistributionFunction;
  archive << s.propertyDensityGrid;
  archive << s.averageEnergyHistogram;
  archive << s.averageNumberOfMoleculesHistogram;
  archive << s.propertyMoleculeProperties;
  archive << s.propertyMSD;
  archive << s.propertyVACF;
  archive << s.writeLammpsData;

  archive << s.propertyNumberOfMoleculesEvolution;
  archive << s.propertyVolumeEvolution;
  archive << s.propertyConservedEnergyEvolution;

  archive << s.interpolationGrids;
  archive << s.externalFieldInterpolationGrid;

#if DEBUG_ARCHIVE
  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, System& s)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > s.versionNumber)
  {
    const std::source_location& location = std::source_location::current();
    throw std::runtime_error(
        std::format("Invalid version reading 'System' at line {} in file {}\n", location.line(), location.file_name()));
  }

  archive >> s.temperature;
  archive >> s.pressure;
  archive >> s.input_pressure;
  archive >> s.pressureTensorDiagonal;
  archive >> s.input_pressureTensorDiagonal;
  archive >> s.beta;
  if (versionNumber >= 3)
  {
    std::uint8_t cellMinimizationType;
    std::uint8_t monoclinicAngleType;
    archive >> cellMinimizationType;
    archive >> monoclinicAngleType;
    s.cellMinimizationType = static_cast<CellMinimizationType>(cellMinimizationType);
    s.monoclinicAngleType = static_cast<MonoclinicAngleType>(monoclinicAngleType);
  }

  archive >> s.heliumVoidFraction;

  archive >> s.numberOfFrameworks;
  archive >> s.numberOfFrameworkAtoms;
  archive >> s.numberOfRigidFrameworkAtoms;

  archive >> s.framework;
  s.numberOfRigidFrameworkAtoms = s.framework && s.framework->rigid ? s.numberOfFrameworkAtoms : 0;
  archive >> s.components;

  archive >> s.equationOfState;
  if (versionNumber >= 4)
  {
    std::uint8_t molecularDynamicsEnsemble;
    archive >> molecularDynamicsEnsemble;
    s.molecularDynamicsEnsemble = static_cast<MolecularDynamicsEnsemble>(molecularDynamicsEnsemble);
  }
  archive >> s.thermostat;
  if (versionNumber >= 4) archive >> s.thermobarostat;

  archive >> s.loadings;

  archive >> s.swappableComponents;
  archive >> s.initialNumberOfMolecules;

  archive >> s.numberOfMoleculesPerComponent;
  archive >> s.numberOfIntegerMoleculesPerComponent;
  archive >> s.numberOfFractionalMoleculesPerComponent;
  archive >> s.numberOfGCFractionalMoleculesPerComponent_CFCMC;
  archive >> s.numberOfPairGCFractionalMoleculesPerComponent_CFCMC;
  archive >> s.numberOfPairSwapFractionalMoleculesPerComponent_CFCMC;
  archive >> s.numberOfPairSwapCBFractionalMoleculesPerComponent_CFCMC;
  archive >> s.numberOfGibbsSwapFractionalMoleculesPerComponent_CFCMC;
  archive >> s.numberOfGibbsFractionalMoleculesPerComponent_CFCMC;
  archive >> s.numberOfParallelReactionFractionalMoleculesPerComponent_CFCMC;
  archive >> s.numberOfSerialReactionFractionalMoleculesPerComponent_CFCMC;
  archive >> s.numberOfReactionFractionalMoleculesPerComponent_CFCMC;

  archive >> s.idealGasEnergiesPerComponent;

  archive >> s.forceField;
  archive >> s.hasExternalField;

  archive >> s.numberOfPseudoAtoms;
  archive >> s.totalNumberOfPseudoAtoms;

  archive >> s.translationalCenterOfMassConstraint;
  archive >> s.translationalDegreesOfFreedom;
  archive >> s.rotationalDegreesOfFreedom;

  archive >> s.timeStep;

  archive >> s.simulationBox;
  archive >> s.containsTheFractionalMolecule;

  archive >> s.atomData;
  archive >> s.atomDynamics;
  archive >> s.moleculeData;
  archive >> s.electricPotential;
  archive >> s.electricField;
  archive >> s.electricFieldNew;

  archive >> s.conservedEnergy;
  archive >> s.referenceEnergy;
  archive >> s.accumulatedDrift;

  archive >> s.rigidEnergies;
  archive >> s.runningEnergies;

  archive >> s.currentExcessPressureTensor;
  archive >> s.currentEnergyStatus;

  archive >> s.numberOfHybridMCSteps;

  archive >> s.eik_xy;
  archive >> s.eik_x;
  archive >> s.eik_y;
  archive >> s.eik_z;
  archive >> s.storedEik;
  archive >> s.fixedFrameworkStoredEik;
  archive >> s.trialEik;
  archive >> s.CoulombicFourierEnergySingleIon;
  archive >> s.netCharge;
  archive >> s.netChargeFramework;
  archive >> s.netChargeAdsorbates;
  archive >> s.netChargePerComponent;

  archive >> s.mc_moves_probabilities;
  archive >> s.mc_moves_statistics;
  archive >> s.mc_moves_cputime;

  archive >> s.reactions;
  archive >> s.tmmc;

  archive >> s.averageEnergies;
  archive >> s.averageLoadings;
  archive >> s.averageEnthalpiesOfAdsorption;
  archive >> s.averagePartialMolarProperties;
  archive >> s.averageTemperature;
  archive >> s.averageTranslationalTemperature;
  archive >> s.averageRotationalTemperature;
  archive >> s.averagePressure;
  archive >> s.averageSimulationBox;
  if (versionNumber >= 5)
  {
    archive >> s.propertyElasticConstantsFluctuation;
    archive >> s.elasticConstantsSampleEvery;
  }

  archive >> s.samplePDBMovie;

  archive >> s.propertyConventionalRadialDistributionFunction;
  archive >> s.propertyRadialDistributionFunction;
  archive >> s.propertyDensityGrid;
  archive >> s.averageEnergyHistogram;
  archive >> s.averageNumberOfMoleculesHistogram;
  archive >> s.propertyMoleculeProperties;
  archive >> s.propertyMSD;
  archive >> s.propertyVACF;
  archive >> s.writeLammpsData;

  archive >> s.propertyNumberOfMoleculesEvolution;
  archive >> s.propertyVolumeEvolution;
  archive >> s.propertyConservedEnergyEvolution;

  // archive >> s.columnNumberOfGridPoints;
  // archive >> s.columnTotalPressure;
  // archive >> s.columnPressureGradient;
  // archive >> s.columnVoidFraction;
  // archive >> s.columnParticleDensity;
  // archive >> s.columnEntranceVelocity;
  // archive >> s.columnLength;
  // archive >> s.columnTimeStep;
  // archive >> s.columnNumberOfTimeSteps;
  // archive >> s.columnAutoNumberOfTimeSteps;
  // archive >> s.mixturePredictionMethod;
  // archive >> s.pressure_range;
  // archive >> s.numberOfCarrierGases;
  // archive >> s.carrierGasComponent;
  // archive >> s.maxIsothermTerms;

  archive >> s.interpolationGrids;
  archive >> s.externalFieldInterpolationGrid;

  // Derived quantities: rebuild the aggregated tail-correction counts from the restored atoms.
  s.computeTailCorrectionCounts();

#if DEBUG_ARCHIVE
  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("System: Error in binary restart\n"));
  }
#endif

  return archive;
}

void System::writeRestartFile(std::size_t systemId)
{
  nlohmann::json json;

  if (!framework.has_value())
  {
    json["SimulationBox"] = simulationBox;
  }

  for (std::size_t component_id = 0; component_id < components.size(); ++component_id)
  {
    std::vector<double3> positions{};

    if (numberOfMoleculesPerComponent[component_id] > 0)
    {
      positions = spanOfIntegerAtomsOfComponent(component_id) | std::views::transform(&Atom::position) |
                  std::ranges::to<std::vector>();
    }
    json[components[component_id].name] = positions;
  }

  std::string fileNameString = std::format("output/restart_{}_{}.s{}.json", temperature, input_pressure, systemId);
  std::ofstream file(fileNameString);
  file << json.dump(2);
}

std::string System::repr() const { return std::string("system test"); }
