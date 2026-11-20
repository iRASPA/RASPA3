module;

module system_preparation;

import std;

import system;
import ideal_gas_reference;

void Preparation::prepareConformations(System& system)
{
  // the CBMC ideal-gas conformation reservoirs are Boltzmann samples at the system temperature
  system.buildConformationReservoirs();
  // so is the recoil-growth openness reference; after the reservoirs, which its grows start from
  system.buildRecoilReferenceConformations();
}

void Preparation::prepareSystems(std::span<System> systems, const IdealGasSettings& idealGas)
{
  prepareIdealGasRosenbluthWeights(systems, idealGas);
}

void Preparation::prepareStatePoint(System& system, const IdealGasSettings& idealGas)
{
  prepareConformations(system);
  prepareIdealGasRosenbluthWeights(std::span<System>(&system, 1uz), idealGas);
}

std::vector<std::vector<double>> Preparation::prepareReplicas(std::span<System> systems,
                                                             const ReplicaSettings& settings)
{
  for (System& system : systems) prepareConformations(system);

  if (settings.acrossTemperatures)
  {
    return prepareIdealGasReferences(systems, settings.ladderTemperatures, settings.otherTemperatures,
                                     settings.idealGas);
  }

  // each replica has a Hamiltonian of its own: its weight is that of its own molecules
  prepareIdealGasRosenbluthWeights(systems, settings.idealGas);
  const std::size_t numberOfComponents = systems.empty() ? 0uz : systems.front().components.size();
  return std::vector<std::vector<double>>(numberOfComponents,
                                          std::vector<double>(settings.otherTemperatures.size(), 0.0));
}
