module;

module system;

import std;

import atom;
import component;
import double3;
import forcefield;
import framework;
import property_loading;
import interactions_ewald;

// Solute tempering (REST2): the Hamiltonian of one component is scaled per replica. The scaling is
// stored in three places that have to stay consistent: the pair parameters of the solute's
// pseudo-atom types in the force field, the component definition (atom charges, intramolecular
// potentials) used when molecules are grown, and the charges of the molecules present in the
// system. The framework and the other components are untouched.

std::vector<bool> System::pseudoAtomTypesOfComponent(std::size_t componentId) const
{
  std::vector<bool> types(forceField.numberOfPseudoAtoms, false);
  for (const Atom& atom : components[componentId].atoms)
  {
    types[static_cast<std::size_t>(atom.type)] = true;
  }
  return types;
}

void System::scaleSoluteHamiltonian(std::size_t componentId, double lambda)
{
  if (componentId >= components.size())
  {
    throw std::runtime_error("[System]: solute tempering: component index out of range\n");
  }
  if (!(lambda > 0.0))
  {
    throw std::runtime_error("[System]: solute tempering: lambda must be positive\n");
  }
  if (soluteTemperingComponent.has_value() && soluteTemperingComponent.value() != componentId)
  {
    throw std::runtime_error("[System]: solute tempering: only one tempered component is supported\n");
  }

  // The pair parameters are scaled per pseudo-atom type, so the types of the solute must be exclusive
  // to the solute.
  const std::vector<bool> soluteTypes = pseudoAtomTypesOfComponent(componentId);
  auto sharedTypes = [&](std::span<const Atom> atoms) -> std::vector<std::string>
  {
    std::set<std::string> names{};
    for (const Atom& atom : atoms)
    {
      if (soluteTypes[static_cast<std::size_t>(atom.type)])
      {
        names.insert(forceField.pseudoAtoms[static_cast<std::size_t>(atom.type)].name);
      }
    }
    return {names.begin(), names.end()};
  };
  for (std::size_t otherId = 0; otherId < components.size(); ++otherId)
  {
    if (otherId == componentId) continue;
    const std::vector<std::string> shared = sharedTypes(components[otherId].atoms);
    if (!shared.empty())
    {
      throw std::runtime_error(std::format(
          "[System]: solute tempering: component '{}' shares the pseudo-atom type(s) {} with component '{}'; "
          "give the tempered component its own pseudo-atom types (duplicate the definitions in the force field "
          "under new names)\n",
          components[componentId].name, shared, components[otherId].name));
    }
  }
  if (framework.has_value())
  {
    const std::vector<std::string> shared = sharedTypes(framework->atoms);
    if (!shared.empty())
    {
      throw std::runtime_error(std::format(
          "[System]: solute tempering: component '{}' shares the pseudo-atom type(s) {} with the framework\n",
          components[componentId].name, shared));
    }
  }

  forceField.scaleSoluteInteractions(soluteTypes, lambda);
  components[componentId].scaleSoluteHamiltonian(lambda);

  soluteTemperingComponent = componentId;
  soluteTemperingLambda *= lambda;

  rescaleSoluteCharges(std::sqrt(lambda));
  rebuildConfigurationDerivedState();
}

void System::rescaleSoluteCharges(double factor)
{
  if (!soluteTemperingComponent.has_value()) return;
  const std::size_t componentId = soluteTemperingComponent.value();

  for (Atom& atom : spanOfMoleculeAtoms())
  {
    if (static_cast<std::size_t>(atom.componentId) == componentId)
    {
      atom.charge *= factor;
    }
  }

  netChargePerComponent[componentId] *= factor;
  netChargeAdsorbates = std::accumulate(netChargePerComponent.begin(), netChargePerComponent.end(), 0.0);
  netCharge = netChargeFramework + netChargeAdsorbates;
}

void System::rebuildConfigurationDerivedState()
{
  forceField.initializeEwaldParameters(simulationBox);
  eik_x.clear();
  eik_y.clear();
  eik_z.clear();
  eik_xy.clear();
  storedEik.clear();
  fixedFrameworkStoredEik.clear();
  trialEik.clear();
  precomputeTotalRigidEnergy();
  runningEnergies = computeTotalEnergies();
  trialEik = storedEik;
  CoulombicFourierEnergySingleIon = Interactions::computeEwaldFourierEnergySingleIon(
      eik_x, eik_y, eik_z, eik_xy, forceField, simulationBox, double3(0.0, 0.0, 0.0), 1.0);
  loadings = LoadingData(components.size(), numberOfIntegerMoleculesPerComponent, simulationBox);
  updateMoleculeAtomInformation();
  computeNumberOfPseudoAtoms();
  computeTailCorrectionCounts();
  netCharge = netChargeFramework + netChargeAdsorbates;
  checkMoleculeIds();
  if (tmmc.doTMMC && !components.empty())
  {
    tmmc.currentLambdaBin = components.front().lambdaGC.currentBin;
  }
}
