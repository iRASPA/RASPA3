module;

export module system_preparation;

import std;

import system;
export import ideal_gas_reference;

/// What a run needs worked out before its first cycle: quantities that follow from the input and are either given
/// or computed. A driver calls one of the functions below once its systems are at their state points; a quantity
/// that was given is used as it is, one that was not is computed, and the status of the component tells which.
/// A run that is resumed from a binary restart file does not come here: what was computed is in the file.
export namespace Preparation
{
/// The settings of the preparation of the replicas of one system.
struct ReplicaSettings
{
  /// The ladder of temperatures the driver replicated the system over (empty for systems declared one by one).
  std::span<const double> ladderTemperatures{};
  /// Temperatures that are reweighted to but not sampled.
  std::span<const double> otherTemperatures{};
  /// Whether the replicas are compared across temperatures through the configuration integral of a flexible
  /// molecule. Not for the scaled Hamiltonians of solute tempering, where such a molecule is refused.
  bool acrossTemperatures{true};
  IdealGasSettings idealGas{};
};

/// The conformations of the isolated molecules of 'system' at its temperature: the ideal-gas conformation
/// reservoirs that seed the growth, and the recoil-growth openness reference (a no-op when recoil growth is off).
void prepareConformations(System& system);

/// Prepares systems that stand on their own, each at the temperature it was made with (whose conformations the
/// construction of the system has built): the ideal-gas Rosenbluth weight of the flexible molecules that need one
/// and were given none.
void prepareSystems(std::span<System> systems, const IdealGasSettings& idealGas = {});

/// Prepares one system whose temperature was changed after it was made: its conformations, and the ideal-gas
/// Rosenbluth weights that were computed at the temperature before.
void prepareStatePoint(System& system, const IdealGasSettings& idealGas = {});

/// Prepares the replicas of one system, each brought to its own state point by the driver: the conformations of
/// each, and the ideal-gas references of the flexible molecules (see 'prepareIdealGasReferences', whose result is
/// returned: per component, the ln of the configuration integral at 'otherTemperatures').
std::vector<std::vector<double>> prepareReplicas(std::span<System> systems, const ReplicaSettings& settings = {});
}  // namespace Preparation
