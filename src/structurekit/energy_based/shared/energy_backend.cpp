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
// A digest of a run of numbers, so that a request can say what its force field, probe and framework were
// without holding copies of them. A collision is as unlikely as a hash collision ever is and costs a reused
// field, so the digest is widened with the counts as well as the values.
void stir(std::size_t &digest, double value)
{
  digest ^= std::hash<double>{}(value) + 0x9e3779b97f4a7c15ULL + (digest << 6) + (digest >> 2);
}

void stir(std::size_t &digest, std::size_t value)
{
  digest ^= std::hash<std::size_t>{}(value) + 0x9e3779b97f4a7c15ULL + (digest << 6) + (digest >> 2);
}

std::size_t digestOf(const PairInteractions &interactions)
{
  std::size_t digest = 0;
  stir(digest, interactions.numberOfTypes);
  stir(digest, interactions.cutOffVDW);
  stir(digest, interactions.cutOffCoulomb);
  for (double charge : interactions.charges) stir(digest, charge);
  for (const PairParameters &pair : interactions.parameters)
  {
    stir(digest, pair.sizeParameter);
    stir(digest, pair.strengthParameter);
    stir(digest, pair.shift);
  }
  return digest;
}

std::size_t digestOf(const LinearProbe &probe)
{
  std::size_t digest = 0;
  stir(digest, probe.sites.size());
  stir(digest, static_cast<std::size_t>(probe.headTailSymmetric ? 1 : 0));
  for (const LinearProbe::Site &site : probe.sites)
  {
    stir(digest, site.type);
    stir(digest, site.offset);
    stir(digest, site.charge);
  }
  return digest;
}

std::size_t digestOf(const Crystal &framework)
{
  std::size_t digest = 0;
  stir(digest, framework.atoms.size());
  for (std::size_t row = 0; row < 3; ++row)
  {
    stir(digest, framework.unitCell.cell[row].x);
    stir(digest, framework.unitCell.cell[row].y);
    stir(digest, framework.unitCell.cell[row].z);
  }
  for (const CrystalAtom &atom : framework.atoms)
  {
    stir(digest, atom.position.x);
    stir(digest, atom.position.y);
    stir(digest, atom.position.z);
    stir(digest, atom.charge);
    stir(digest, atom.type);
  }
  // Where the atoms are is held twice, in Ångström and in fractions of the cell, and it is the fractions
  // the builders of the fields read. The two agree in a framework that was read from a file, but nothing
  // makes them, and a framework changed in the one and not in the other is a different landscape.
  stir(digest, framework.fractionalPositions.size());
  for (const auto &fractional : framework.fractionalPositions)
  {
    stir(digest, fractional.x);
    stir(digest, fractional.y);
    stir(digest, fractional.z);
  }
  return digest;
}

// Everything a landscape depends on. The framework, the force field and the probe are each taken by name
// and by a digest of their contents, so that a second structure of the same name, a changed cutoff or a
// probe given other parameters is seen as a different request and not handed the field held from the last.
struct FieldRequest
{
  std::string backend;
  std::string blocking;
  std::string frameworkName;
  std::string probeName;
  std::size_t numberOfAtoms{0};
  double volume{0.0};
  uint3 gridSize{0, 0, 0};
  std::size_t numberOfOrientations{0};
  double temperature{0.0};
  bool electrostatics{false};
  double relativePrecision{0.0};
  std::size_t interactionsDigest{0};
  std::size_t probeDigest{0};
  std::size_t frameworkDigest{0};

  bool sameAs(const FieldRequest &other) const
  {
    return this->backend == other.backend && this->blocking == other.blocking &&
           this->frameworkName == other.frameworkName && this->probeName == other.probeName &&
           this->numberOfAtoms == other.numberOfAtoms && this->volume == other.volume &&
           this->gridSize.x == other.gridSize.x && this->gridSize.y == other.gridSize.y &&
           this->gridSize.z == other.gridSize.z && this->numberOfOrientations == other.numberOfOrientations &&
           this->temperature == other.temperature && this->electrostatics == other.electrostatics &&
           this->relativePrecision == other.relativePrecision && this->interactionsDigest == other.interactionsDigest &&
           this->probeDigest == other.probeDigest && this->frameworkDigest == other.frameworkDigest;
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
  request.blocking = backend.blocking;
  request.frameworkName = framework.name;
  request.probeName = probe.name;
  request.numberOfAtoms = framework.atoms.size();
  request.volume = framework.unitCell.volume;
  request.gridSize = gridSize;
  request.numberOfOrientations = numberOfOrientations;
  request.temperature = temperature;
  request.electrostatics = wantsElectrostatics;
  request.relativePrecision = relativePrecision;
  request.interactionsDigest = digestOf(interactions);
  request.probeDigest = digestOf(probe);
  request.frameworkDigest = digestOf(framework);

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
