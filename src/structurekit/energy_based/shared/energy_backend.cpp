module;

module energy_shared_energy_backend;

import std;

import uint3;
import crystal;
import pair_interactions;

import energy_shared_linear_probe;
import energy_shared_molecular_energy_grid;
import energy_shared_electrostatic_potential_grid;

namespace
{
// Everything a landscape depends on. The framework is taken by name, by how many atoms are in it and by the
// volume it occupies rather than by its contents: two frameworks agreeing on all three within one run are
// the same framework read twice, and a run that reads a second structure changes the name.
struct FieldRequest
{
  std::string backend;
  std::string frameworkName;
  std::string probeName;
  std::size_t numberOfAtoms{0};
  double volume{0.0};
  uint3 gridSize{0, 0, 0};
  std::size_t numberOfOrientations{0};
  double temperature{0.0};
  bool electrostatics{false};
  double relativePrecision{0.0};

  bool sameAs(const FieldRequest &other) const
  {
    return this->backend == other.backend && this->frameworkName == other.frameworkName &&
           this->probeName == other.probeName && this->numberOfAtoms == other.numberOfAtoms &&
           this->volume == other.volume && this->gridSize.x == other.gridSize.x &&
           this->gridSize.y == other.gridSize.y && this->gridSize.z == other.gridSize.z &&
           this->numberOfOrientations == other.numberOfOrientations && this->temperature == other.temperature &&
           this->electrostatics == other.electrostatics && this->relativePrecision == other.relativePrecision;
  }
};

std::mutex heldMutex;
bool haveHeld = false;
FieldRequest heldRequest;
MolecularField heldField;
}  // namespace


void forgetMolecularField()
{
  std::scoped_lock lock(heldMutex);
  haveHeld = false;
  heldRequest = FieldRequest{};
  heldField = MolecularField{};
}


MolecularField buildMolecularField(const EnergyBackend &backend, const PairInteractions &interactions,
                                   const Crystal &framework, const LinearProbe &probe, uint3 gridSize,
                                   std::size_t numberOfOrientations, double temperature, bool useElectrostatics,
                                   double relativePrecision)
{
  MolecularField field;

  bool wantsElectrostatics = useElectrostatics && probe.isCharged();

  FieldRequest request;
  request.backend = backend.name;
  request.frameworkName = framework.name;
  request.probeName = probe.name;
  request.numberOfAtoms = framework.atoms.size();
  request.volume = framework.unitCell.volume;
  request.gridSize = gridSize;
  request.numberOfOrientations = numberOfOrientations;
  request.temperature = temperature;
  request.electrostatics = wantsElectrostatics;
  request.relativePrecision = relativePrecision;

  std::scoped_lock lock(heldMutex);
  if (haveHeld && heldRequest.sameAs(request))
  {
    field = heldField;
    field.reused = true;
    return field;
  }

  // Whatever is held goes before the next is built, so that two landscapes are never in memory at once.
  haveHeld = false;
  heldField = MolecularField{};

  if (wantsElectrostatics)
  {
    field.potential = backend.electrostaticPotentialGrid(interactions, framework, gridSize, relativePrecision);

    if (std::getenv("RASPA_VERIFY_ELECTROSTATICS") != nullptr)
    {
      std::ofstream check;
      check.open(framework.name + ".electrostatics.check.txt");
      check << ElectrostaticPotentialGrid::splitIndependenceCheck(interactions, framework);
      check.close();
    }
  }

  field.grid = backend.molecularEnergyGrid(interactions, framework, probe, gridSize, numberOfOrientations, temperature,
                                           wantsElectrostatics ? &field.potential : nullptr);

  heldRequest = request;
  heldField = field;
  haveHeld = true;

  return field;
}
