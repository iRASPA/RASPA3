module;

module energy_shared_probe_energy_grid;

import std;

import double3;
import unit_cell;
import crystal;
import units;

ProbeEnergyGrid::ProbeEnergyGrid() {}

ProbeEnergyGrid::~ProbeEnergyGrid() {}

double3 ProbeEnergyGrid::fractionalPosition(std::size_t i, std::size_t j, std::size_t k) const
{
  return double3(static_cast<double>(i) / static_cast<double>(this->gridSize.x),
                 static_cast<double>(j) / static_cast<double>(this->gridSize.y),
                 static_cast<double>(k) / static_cast<double>(this->gridSize.z));
}

double3 ProbeEnergyGrid::cartesianPosition(std::size_t i, std::size_t j, std::size_t k) const
{
  return this->unitCell.cell * this->fractionalPosition(i, j, k);
}

double ProbeEnergyGrid::minimumEnergy() const
{
  if (this->energy.empty()) return 0.0;
  return static_cast<double>(*std::ranges::min_element(this->energy));
}

double3 ProbeEnergyGrid::spacing() const
{
  double3 a = this->unitCell.cell[0];
  double3 b = this->unitCell.cell[1];
  double3 c = this->unitCell.cell[2];

  return double3(a.length() / static_cast<double>(this->gridSize.x),
                 b.length() / static_cast<double>(this->gridSize.y),
                 c.length() / static_cast<double>(this->gridSize.z));
}


void ProbeEnergyGrid::writeTessellation(const Crystal &framework) const
{
  if (this->strongestAtom.empty()) return;

  double3 spacing = this->spacing();
  const double voxel = this->unitCell.volume / static_cast<double>(this->numberOfVoxels());

  // How much of the cell falls to each atom, and how much of the energy. The first is what a tessellation is
  // for; the second says how much of what the probe feels anywhere in the cell each atom is answerable for,
  // which the geometric cells cannot report because they are not built from anything that adds up.
  std::vector<std::size_t> voxelsPerAtom(framework.atoms.size(), 0);
  std::vector<double> energyPerAtom(framework.atoms.size(), 0.0);
  for (std::size_t v = 0; v < this->strongestAtom.size(); ++v)
  {
    const std::int32_t atom = this->strongestAtom[v];
    if (atom < 0 || static_cast<std::size_t>(atom) >= voxelsPerAtom.size()) continue;
    ++voxelsPerAtom[static_cast<std::size_t>(atom)];

    // Only where the probe is actually held, since a point buried in a wall sits at the ceiling and would
    // otherwise swamp every real number in the column.
    const double here = static_cast<double>(this->energy[v]);
    if (here < 0.0) energyPerAtom[static_cast<std::size_t>(atom)] += here;
  }

  std::ofstream myfile;
  myfile.open(framework.name + "." + this->probeName + ".energy.tessellation." + this->backend + ".txt");
  std::print(myfile, "# Tessellation of the cell by strongest attraction (probe energy grid)\n");
  std::print(myfile, "# Crystal: {}\n", framework.name);
  std::print(myfile, "# Number of framework atoms: {}\n", framework.atoms.size());
  std::print(myfile, "# Crystal volume: {} [Å³]\n", this->unitCell.volume);
  std::print(myfile, "# Probe: {}, epsilon {} [K], sigma {} [Å]\n", this->probeName,
             this->probeEpsilon * Units::EnergyToKelvin, this->probeSigma);
  std::print(myfile, "# Grid: {} x {} x {} points, spacing {:.5f} x {:.5f} x {:.5f} [Å]\n", this->gridSize.x,
             this->gridSize.y, this->gridSize.z, spacing.x, spacing.y, spacing.z);
  std::print(myfile, "# Volume one grid point stands for: {} [Å³]\n", voxel);
  std::print(myfile, "# Cutoff: {} [Å]\n", this->cutOff);
  std::print(myfile, "# Timing ({}): {} [s]\n", this->backend, this->seconds);
  std::print(myfile, "#\n");
  std::print(myfile, "# A point belongs to the atom whose own contribution to the energy there, summed over\n");
  std::print(myfile, "# that atom's periodic images within the cutoff, is the most negative: the atom pulling\n");
  std::print(myfile, "# hardest on the probe rather than the one whose surface is nearest. Where no atom is\n");
  std::print(myfile, "# attractive at all, which is the inside of a wall, the largest contribution is taken\n");
  std::print(myfile, "# instead, and that is the nearest atom.\n");
  std::print(myfile, "#\n");
  std::print(myfile, "# In a pocket these are the nearest-surface cells, since the nearest atom is also the\n");
  std::print(myfile, "# strongest. In a window between two cavities they need not be: a nearby atom at the\n");
  std::print(myfile, "# bottom of its well pulls weakly while one further out on the steep flank dominates, and\n");
  std::print(myfile, "# it is those places that decide what a molecule can do.\n");
  std::print(myfile, "#\n");
  std::print(myfile, "# column 1: atom\n");
  std::print(myfile, "# column 2: grid points it pulls hardest on\n");
  std::print(myfile, "# column 3: the volume that comes to [Å³]\n");
  std::print(myfile, "# column 4: as a fraction of the cell\n");
  std::print(myfile, "# column 5: the binding energy it accounts for, summed over its own points [K]\n");
  std::print(myfile, "#   atom    points   volume [Å³]     fraction     binding [K]\n");
  for (std::size_t atom = 0; atom < voxelsPerAtom.size(); ++atom)
  {
    const double share = static_cast<double>(voxelsPerAtom[atom]) * voxel;
    std::print(myfile, "CrystalAtom: {:6} {:9} {:13.5f} {:12.8f} {:15.4f}\n", atom, voxelsPerAtom[atom], share,
               share / this->unitCell.volume, energyPerAtom[atom] * Units::EnergyToKelvin);
  }

  std::print(myfile, "\n");
  std::print(myfile, "# The grid itself: the atom holding each point and what the probe feels there, in K.\n");
  std::print(myfile, "# Ordered with x varying fastest.\n");
  std::print(myfile, "#   i    j    k     atom     energy [K]\n");
  for (std::size_t k = 0; k < this->gridSize.z; ++k)
  {
    for (std::size_t j = 0; j < this->gridSize.y; ++j)
    {
      for (std::size_t i = 0; i < this->gridSize.x; ++i)
      {
        const std::size_t voxelIndex = this->voxelIndex(i, j, k);
        std::print(myfile, "{:5} {:4} {:4} {:8} {:15.5f}\n", i, j, k, this->strongestAtom[voxelIndex],
                   static_cast<double>(this->energy[voxelIndex]) * Units::EnergyToKelvin);
      }
    }
  }

  myfile.close();
}
