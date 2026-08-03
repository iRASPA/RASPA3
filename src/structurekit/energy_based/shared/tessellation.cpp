module;

module energy_shared_tessellation;

import std;

import uint3;
import double3;
import double3x3;
import skspacegroupdatabase;
import crystal;
import pair_interactions;
import units;
import energy_shared_isosurface;
import energy_shared_energy_backend;
import energy_shared_linear_probe;
import energy_shared_probe_energy_grid;
import energy_shared_molecular_energy_grid;
import energy_shared_electrostatic_potential_grid;


EnergyTessellation::EnergyTessellation() {}


EnergyTessellation::~EnergyTessellation() {}


namespace
{
// What the two routes differ in, which is nothing but the field, the labels, and what the report calls them.
// The division itself is one piece of arithmetic and is written once.
struct Division
{
  std::span<const float> energy;
  std::span<const std::int32_t> label;
  uint3 gridSize;

  std::string fileName;
  std::string backendName;

  // The lines describing the probe, which are the only part of the header that is not common to both.
  std::vector<std::string> probeLines;

  // What the field and its labels cost, which is one figure for a single atom and two for a molecule.
  std::vector<std::string> timingLines;

  // What the label is decided on, which for a molecule is not quite what it is for a single atom.
  std::vector<std::string> ruleLines;
};
}  // namespace


// The division itself: every point of the cell handed to an atom, and the volume, the void, the surface and
// the binding that fall in each atom's share added up. Neither the arithmetic nor the report knows whether
// the label came from one atom's pull or from a molecule's, which is what keeps the two comparable.
static void divideAmongAtoms(EnergyTessellation &tessellation, const Crystal &framework, double isoValue,
                             const EnergyBackend &backend, const Division &division)
{
  const uint3 gridSize = division.gridSize;
  const std::size_t numberOfVoxels = division.energy.size();
  if (numberOfVoxels == 0) return;

  std::chrono::steady_clock::time_point time_begin = std::chrono::steady_clock::now();

  const std::size_t numberOfAtoms = framework.atoms.size();
  tessellation.atoms.assign(numberOfAtoms, EnergyAtomShare{});
  for (EnergyAtomShare &share : tessellation.atoms) share.deepestEnergy = 0.0;

  const double volume = framework.unitCell.volume;
  const double voxelVolume = volume / static_cast<double>(numberOfVoxels);

  for (std::size_t v = 0; v < numberOfVoxels; ++v)
  {
    const std::int32_t atom = division.label[v];
    if (atom < 0 || static_cast<std::size_t>(atom) >= numberOfAtoms) continue;
    EnergyAtomShare &share = tessellation.atoms[static_cast<std::size_t>(atom)];

    ++share.numberOfVoxels;

    const double here = static_cast<double>(division.energy[v]);
    if (here < isoValue) share.voidVolume += voxelVolume;
    if (here < 0.0) share.bindingEnergy += here;
    share.deepestEnergy = std::min(share.deepestEnergy, here);
  }

  // The surface of the void, cut up by the same labels. A triangle is given to the atom holding the grid
  // point nearest its centre, which is the same rule the geometric route uses and is unambiguous as long as
  // a triangle stays inside one voxel, which marching cubes guarantees.
  std::vector<double3> corners = backend.isosurfaceTriangles(division.energy, gridSize, isoValue);

  const double3x3 cell = framework.unitCell.cell;
  const double largestPlausible = largestPlausibleTriangleArea(cell, gridSize);

  auto voxelIndex = [&](std::size_t i, std::size_t j, std::size_t k)
  { return (k * gridSize.y + j) * gridSize.x + i; };

  for (std::size_t i = 0; i + 2 < corners.size(); i += 3)
  {
    const double3 a = corners[i];
    const double3 b = corners[i + 1];
    const double3 c = corners[i + 2];

    const double3 pa = cell * a;
    const double3 pb = cell * b;
    const double3 pc = cell * c;
    const double area = 0.5 * double3::cross(pb - pa, pc - pa).length();
    if (!std::isfinite(area)) continue;

    // A triangle is confined to a single voxel, so anything larger is numerical debris rather than surface.
    if (area > largestPlausible)
    {
      ++tessellation.numberOfRejectedTriangles;
      continue;
    }

    ++tessellation.numberOfTriangles;
    tessellation.totalArea += area;

    const double3 centre((a.x + b.x + c.x) / 3.0, (a.y + b.y + c.y) / 3.0, (a.z + b.z + c.z) / 3.0);
    const std::size_t i0 = static_cast<std::size_t>(centre.x * static_cast<double>(gridSize.x)) % gridSize.x;
    const std::size_t j0 = static_cast<std::size_t>(centre.y * static_cast<double>(gridSize.y)) % gridSize.y;
    const std::size_t k0 = static_cast<std::size_t>(centre.z * static_cast<double>(gridSize.z)) % gridSize.z;

    const std::int32_t atom = division.label[voxelIndex(i0, j0, k0)];
    if (atom >= 0 && static_cast<std::size_t>(atom) < numberOfAtoms)
    {
      tessellation.atoms[static_cast<std::size_t>(atom)].area += area;
    }
    else
    {
      tessellation.undecidedArea += area;
    }
  }

  const double mass = framework.mass;
  const double toGravimetric = (mass > 0.0) ? Units::Angstrom * Units::Angstrom * Units::AvogadroConstant / mass : 0.0;

  for (EnergyAtomShare &share : tessellation.atoms)
  {
    share.volume = static_cast<double>(share.numberOfVoxels) * voxelVolume;
    share.volumeFraction = share.volume / volume;
    share.gravimetricArea = share.area * toGravimetric;

    tessellation.totalVoidVolume += share.voidVolume;
    tessellation.totalBindingEnergy += share.bindingEnergy;
  }
  tessellation.totalGravimetricArea = tessellation.totalArea * toGravimetric;

  std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - time_begin;
  tessellation.seconds = elapsed.count();

  std::ofstream myfile;
  myfile.open(division.fileName);
  std::print(myfile, "# Volume and area divided among the atoms by strongest attraction\n");
  std::print(myfile, "# Crystal: {}\n", framework.name);
  std::print(myfile, "# Space-group HM-symbol: {}\n",
             SKSpaceGroupDataBase::spaceGroupData[framework.spaceGroupHallNumber].HMString());
  std::print(myfile, "# Number of framework atoms: {}\n", numberOfAtoms);
  std::print(myfile, "# Crystal volume: {} [Å³]\n", volume);
  std::print(myfile, "# Crystal mass: {} [g/mol]\n", mass);
  for (const std::string &line : division.probeLines) std::print(myfile, "{}\n", line);
  std::print(myfile, "# Grid: {} x {} x {} points\n", gridSize.x, gridSize.y, gridSize.z);
  std::print(myfile, "# Iso-value: {} [internal], {:.4f} [K]\n", isoValue, isoValue * Units::EnergyToKelvin);
  std::print(myfile, "# Triangles: {} kept, {} discarded as too large for a voxel\n", tessellation.numberOfTriangles,
             tessellation.numberOfRejectedTriangles);
  for (const std::string &line : division.timingLines) std::print(myfile, "{}\n", line);
  std::print(myfile, "# Timing ({}): {} [s] for the division\n", division.backendName, tessellation.seconds);
  std::print(myfile, "#\n");
  for (const std::string &line : division.ruleLines) std::print(myfile, "{}\n", line);
  std::print(myfile, "#\n");
  std::print(myfile, "# Unlike the geometric division, this one depends on the probe. A probe that binds\n");
  std::print(myfile, "# strongly to one element and weakly to another will hand the first element ground that a\n");
  std::print(myfile, "# hard sphere of the same size would have given to the second. An element the guest\n");
  std::print(myfile, "# cannot reach, silicon in a pure-silica framework being the standing example, is left\n");
  std::print(myfile, "# with almost nothing here while the geometric division still gives it a share of the\n");
  std::print(myfile, "# cell. That holds whether the model leaves its dispersion out, as the oxygen-only ones\n");
  std::print(myfile, "# do, or gives it the weak short-ranged one of TraPPE-zeo: either way the oxygens that\n");
  std::print(myfile, "# coordinate it are nearer to the guest and pull harder.\n");
  std::print(myfile, "#\n");
  std::print(myfile, "# The boundaries between cells are where two atoms tie, and in a symmetric framework the\n");
  std::print(myfile, "# ties run all along them. Which atom wins is then settled by rounding, so the two\n");
  std::print(myfile, "# backends hand about one boundary point in a hundred to different atoms at this grid\n");
  std::print(myfile, "# size, falling off as the reciprocal of it. Totals are unaffected; single atoms' shares\n");
  std::print(myfile, "# are, by a few percent, and should be read with that in mind.\n");
  std::print(myfile, "#\n");
  std::print(myfile, "# column 1: atom\n");
  std::print(myfile, "# column 2: its share of the cell [Å³]\n");
  std::print(myfile, "# column 3: as a fraction of the cell\n");
  std::print(myfile, "# column 4: the void within its share [Å³]\n");
  std::print(myfile, "# column 5: the iso-surface within its share [Å²]\n");
  std::print(myfile, "# column 6: the same per unit mass [m²/g]\n");
  std::print(myfile, "# column 7: the binding it accounts for [K]\n");
  std::print(myfile, "# column 8: the deepest point of its share [K]\n");
  std::print(myfile, "#   atom   volume [Å³]     fraction     void [Å³]     area [Å²]      [m²/g]"
                     "     binding [K]     deepest [K]\n");
  for (std::size_t atom = 0; atom < numberOfAtoms; ++atom)
  {
    const EnergyAtomShare &share = tessellation.atoms[atom];
    std::print(myfile, "CrystalAtom: {:5} {:13.5f} {:12.8f} {:13.5f} {:13.5f} {:11.4f} {:15.4f} {:15.4f}\n", atom,
               share.volume, share.volumeFraction, share.voidVolume, share.area, share.gravimetricArea,
               share.bindingEnergy * Units::EnergyToKelvin, share.deepestEnergy * Units::EnergyToKelvin);
  }

  std::print(myfile, "\n");
  std::print(myfile, "{:36} {:16.5f} [Å²] {:12.4f} [m²/g]\n", "Total iso-surface:", tessellation.totalArea,
             tessellation.totalGravimetricArea);
  std::print(myfile, "{:36} {:16.5f} [Å²]\n", "Surface in no atom's cell:", tessellation.undecidedArea);
  std::print(myfile, "{:36} {:16.5f} [Å³] {:12.8f}\n", "Total void at this level:", tessellation.totalVoidVolume,
             tessellation.totalVoidVolume / volume);
  std::print(myfile, "{:36} {:16.4f} [K]\n", "Total binding accounted for:",
             tessellation.totalBindingEnergy * Units::EnergyToKelvin);
  myfile.close();
}


void EnergyTessellation::run(const EnergyBackend &backend, const PairInteractions &interactions, const Crystal &framework,
                             std::string probePseudoAtom, double isoValue, uint3 gridSize)
{
  this->grid = backend.probeEnergyGrid(interactions, framework, probePseudoAtom, gridSize);
  this->isoValue = isoValue;
  this->isMolecular = false;

  if (this->grid.numberOfVoxels() == 0) return;

  Division division;
  division.energy = this->grid.energy;
  division.label = this->grid.strongestAtom;
  division.gridSize = gridSize;
  division.fileName =
      framework.name + "." + this->grid.probeName + ".energy.shares." + this->grid.backend + ".txt";
  division.backendName = this->grid.backend;
  division.timingLines = {
      std::format("# Timing ({}): {} [s] for the field and its labels", this->grid.backend, this->grid.seconds)};
  division.probeLines = {
      std::format("# Probe: {}, epsilon {} [K], sigma {} [Å]", this->grid.probeName,
                  this->grid.probeEpsilon * Units::EnergyToKelvin, this->grid.probeSigma),
      std::format("# Cutoff: {} [Å]", this->grid.cutOff)};
  division.ruleLines = {
      "# A point of the cell belongs to the atom whose own contribution to the energy there is",
      "# the most negative: the atom pulling hardest on the probe, not the one whose surface is",
      "# nearest. The two questions have the same answer through most of a pocket and different",
      "# answers in a window, where a nearby atom near the bottom of its well pulls weakly and",
      "# one further out on the steep flank of the potential does not."};

  divideAmongAtoms(*this, framework, isoValue, backend, division);

  // The whole-cell division, in the same form the clearance route writes it, so the two can be laid side by
  // side.
  this->grid.writeTessellation(framework);
}


void EnergyTessellation::run(const EnergyBackend &backend, const PairInteractions &interactions, const Crystal &framework,
                             const LinearProbe &probe, double isoValue, uint3 gridSize,
                             std::size_t numberOfOrientations, double temperature, bool useElectrostatics,
                             double relativePrecision)
{
  MolecularField field = buildMolecularField(backend, interactions, framework, probe, gridSize, numberOfOrientations,
                                             temperature, useElectrostatics, relativePrecision);
  this->molecularGrid = field.grid;
  this->isoValue = isoValue;
  this->isMolecular = true;

  if (this->molecularGrid.numberOfVoxels() == 0) return;

  std::chrono::steady_clock::time_point time_begin = std::chrono::steady_clock::now();
  const ElectrostaticPotentialGrid *potential =
      (field.potential.numberOfVoxels() > 0) ? &field.potential : nullptr;
  std::vector<std::int32_t> labels = backend.molecularStrongestAtom(
      interactions, framework, probe, gridSize, this->molecularGrid.numberOfOrientations, temperature, potential);
  std::chrono::duration<double> labelling = std::chrono::steady_clock::now() - time_begin;

  Division division;
  division.energy = this->molecularGrid.freeEnergy;
  division.label = labels;
  division.gridSize = gridSize;
  division.fileName =
      framework.name + "." + probe.name + ".energy.shares." + this->molecularGrid.backend + ".txt";
  division.backendName = this->molecularGrid.backend;
  division.timingLines = {
      std::format("# Timing ({}): {} [s] for the landscape, {} [s] for the labels", this->molecularGrid.backend,
                  this->molecularGrid.seconds, labelling.count())};
  division.probeLines = {
      std::format("# Molecule: {}, {} sites, {:.4f} [Å] end to end", probe.name, probe.sites.size(), probe.length()),
      std::format("# Orientations sampled: {} over the {}", this->molecularGrid.numberOfOrientations,
                  this->molecularGrid.overHemisphere ? "hemisphere, the molecule being the same end for end"
                                                     : "sphere"),
      std::format("# Temperature: {} [K]", temperature),
      std::format("# Cutoff: {} [Å]", this->molecularGrid.cutOff),
      this->molecularGrid.chargesIncluded
          ? std::format("# Electrostatics: Ewald, split at alpha = {:.5f} [1/Å], {} wave vectors",
                        this->molecularGrid.ewaldAlpha, this->molecularGrid.numberOfWaveVectors)
          : std::string("# Electrostatics: none")};
  division.ruleLines = {
      "# A point of the cell belongs to the atom whose own contribution to the energy there is",
      "# the most negative, averaged over the ways the molecule can be turned at that point. The",
      "# average is weighted by exp(-U/kT) and not taken flat, which is the whole of the",
      "# definition worth arguing over: a molecule in a window can point almost nowhere, and an",
      "# average over the orientations it cannot take up is an average over things that never",
      "# happen. Weighted, the question asked is which atom pulls hardest on the molecule as it",
      "# actually sits.",
      "#",
      "# Only the part of the energy that belongs to one atom decides it: the dispersion, and the",
      "# near half of the electrostatic sum, which is a sum over pairs. The far half is a sum over",
      "# waves and belongs to the framework as a whole, so it takes no part in the division. It is",
      "# in the weights, which are about how the molecule sits rather than about who holds it.",
      "#",
      "# The field the volumes and the surface are read from is the free energy, so the void here",
      "# is the void of the molecule and not of its best-oriented self."};

  if (this->molecularGrid.chargesIncluded)
  {
    division.ruleLines.insert(
        division.ruleLines.end(),
        {"#",
         "# Read the division below with the charges in mind, because for a charged molecule they",
         "# decide it. The near half of the electrostatic sum is damped by erfc, which at four",
         "# Ångström and the splitting above is a few percent of the bare Coulomb term, so what an",
         "# atom is credited with here is the screened part of its pull and not the whole of it.",
         "# Screening falls off with distance faster than the Coulomb term does, so the division is",
         "# tipped towards near atoms by more than the physics alone would tip it.",
         "#",
         "# The effect of this is not small and is not subtle. Silicon in a pure-silica framework",
         "# takes almost none of the cell with the charges off, not for want of dispersion, which",
         "# the usual models do give it, but because the oxygens coordinating it stand between it",
         "# and the guest and a dispersion term falls off as the sixth power. A Coulomb term does",
         "# not, and silicon carries twice the charge of oxygen, so with the charges on it takes",
         "# most of the cell. Both are statements about the model rather than about the framework,",
         "# and neither should be quoted without it."});
  }

  divideAmongAtoms(*this, framework, isoValue, backend, division);

  // The whole-cell division, in the same form the other routes write it.
  ProbeEnergyGrid asField;
  asField.gridSize = gridSize;
  asField.unitCell = framework.unitCell;
  asField.energy = this->molecularGrid.freeEnergy;
  asField.strongestAtom = labels;
  asField.probeName = probe.name;
  asField.cutOff = this->molecularGrid.cutOff;
  asField.ceiling = this->molecularGrid.ceiling;
  asField.backend = this->molecularGrid.backend;
  asField.writeTessellation(framework);
}
