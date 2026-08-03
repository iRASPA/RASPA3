module;

module energy_shared_molecular_energy_barrier;

import std;

import int3;
import uint3;
import double3;
import skspacegroupdatabase;
import crystal;
import pair_interactions;
import units;
import energy_shared_linear_probe;
import energy_shared_energy_backend;
import energy_shared_molecular_energy_grid;
import energy_shared_energy_barrier;
import energy_shared_electrostatic_potential_grid;


MolecularEnergyBarrier::MolecularEnergyBarrier() {}


MolecularEnergyBarrier::~MolecularEnergyBarrier() {}


void MolecularEnergyBarrier::run(const EnergyBackend &backend, const PairInteractions &interactions, const Crystal &framework, const LinearProbe &probe,
                                 uint3 gridSize, std::size_t numberOfOrientations, double temperature,
                                 bool useElectrostatics, double relativePrecision)
{
  MolecularField field = buildMolecularField(backend, interactions, framework, probe, gridSize, numberOfOrientations,
                                             temperature, useElectrostatics, relativePrecision);
  this->potential = field.potential;
  this->grid = field.grid;

  if (this->grid.numberOfVoxels() == 0) return;

  this->fromFreeEnergy = EnergyBarrier::fromField(this->grid.gridSize, this->grid.freeEnergy, temperature);
  this->fromMinimumEnergy = EnergyBarrier::fromField(this->grid.gridSize, this->grid.minimumEnergy, temperature);
  this->orientationalPenalty = this->fromFreeEnergy.percolationBarrier - this->fromMinimumEnergy.percolationBarrier;

  const double toKelvin = Units::EnergyToKelvin;
  const double kT = Units::KB * temperature;
  const double never = std::numeric_limits<double>::max();
  double3 spacing = this->grid.spacing();

  auto asKelvin = [&](double energy) { return energy * toKelvin; };
  auto asKJPerMol = [&](double energy) { return energy * Units::EnergyToKJPerMol; };

  std::ofstream myfile;
  myfile.open(framework.name + "." + probe.name + ".grid.barrier." + this->grid.backend + ".txt");
  std::print(myfile, "# Percolation barrier for a rigid linear molecule, from its energy landscape\n");
  std::print(myfile, "# Crystal: {}\n", framework.name);
  std::print(myfile, "# Space-group HM-symbol: {}\n",
             SKSpaceGroupDataBase::spaceGroupData[framework.spaceGroupHallNumber].HMString());
  std::print(myfile, "# Number of framework atoms: {}\n", framework.atoms.size());
  std::print(myfile, "# Crystal volume: {} [Å³]\n", framework.unitCell.volume);
  std::print(myfile, "# Molecule: {}, {} sites, {:.4f} [Å] end to end\n", probe.name, probe.sites.size(),
             probe.length());
  for (const LinearProbe::Site &site : probe.sites)
  {
    std::print(myfile, "#   site {:8} at {:+8.4f} [Å]  strength/kʙ: {:10.4f} [K]  size: {:8.4f} [Å]  charge: {:+8.5f}\n",
               site.name, site.offset, interactions[site.type].strengthParameter * toKelvin,
               interactions[site.type].sizeParameter, site.charge);
  }
  std::print(myfile, "# Orientations sampled: {} over the {}\n", this->grid.numberOfOrientations,
             this->grid.overHemisphere ? "hemisphere, the molecule being the same end for end"
                                       : "whole sphere");
  std::print(myfile, "# Cutoff: {} [Å], periodic images searched: {} x {} x {} either way\n", this->grid.cutOff,
             this->grid.numberOfImageShells.x, this->grid.numberOfImageShells.y, this->grid.numberOfImageShells.z);
  std::print(myfile, "# Grid: {} x {} x {} points, spacing {:.5f} x {:.5f} x {:.5f} [Å]\n", this->grid.gridSize.x,
             this->grid.gridSize.y, this->grid.gridSize.z, spacing.x, spacing.y, spacing.z);
  std::print(myfile, "# Temperature: {} [K]\n", temperature);
  std::print(myfile, "# Timing ({}): {} [s] for the landscape\n", this->grid.backend, this->grid.seconds);
  std::print(myfile, "# CPU Timing: {} [s] and {} [s] for the two sweeps\n", this->fromFreeEnergy.sweepSeconds,
             this->fromMinimumEnergy.sweepSeconds);

  if (this->grid.chargesIncluded)
  {
    std::print(myfile, "# Electrostatics: Ewald, split at alpha = {:.5f} [1/Å], {} wave vectors, {} relative\n",
               this->potential.alpha, this->potential.numberOfWaveVectors, this->potential.relativePrecision);
    std::print(myfile, "# Coulomb cutoff: {} [Å], net charge in the cell: {:+.3e}\n", this->potential.cutOff,
               this->potential.netCharge);
    std::print(myfile, "# Timing ({}): {} [s] for the framework potential\n", this->potential.backend,
               this->potential.seconds);

    if (this->potential.largestFrameworkCharge == 0.0)
    {
      std::print(myfile, "#\n");
      std::print(myfile, "# WARNING: every framework atom carries zero charge, so the potential is zero and the\n");
      std::print(myfile, "# electrostatics below amount to nothing. A CIF with no charge column leaves the atoms\n");
      std::print(myfile, "# neutral; pass --charges-from pseudo-atoms to take them from the force field instead.\n");
    }
  }

  if (this->grid.chargesIgnored)
  {
    std::print(myfile, "#\n");
    std::print(myfile, "# WARNING: this molecule carries partial charges and they have not been acted on. Only\n");
    std::print(myfile, "# the dispersion and repulsion are in the landscape below. For a molecule whose behaviour\n");
    std::print(myfile, "# in a framework is largely electrostatic, carbon dioxide above all, the numbers here are\n");
    std::print(myfile, "# not the ones to quote: they are the shape of the molecule without its quadrupole.\n");
  }

  std::print(myfile, "#\n");
  std::print(myfile, "# A molecule with a shape has an energy for every way it can be turned, and the barrier\n");
  std::print(myfile, "# depends on which of those is taken as the energy at a point. Two readings are given.\n");
  std::print(myfile, "#\n");
  std::print(myfile, "# The free energy, -kT ln <exp(-U/kT)> over orientations, is the landscape a molecule that\n");
  std::print(myfile, "# turns freely on the way through actually climbs. It is the one to quote. The minimum\n");
  std::print(myfile, "# over orientations is the landscape it would climb if it were always turned the best way,\n");
  std::print(myfile, "# which is the zero-temperature limit and a lower bound: it grants the best orientation\n");
  std::print(myfile, "# everywhere and charges nothing for taking it up.\n");
  std::print(myfile, "#\n");
  std::print(myfile, "# The gap between the two is what the molecule pays in orientation rather than in energy.\n");
  std::print(myfile, "# A wide channel costs nothing, since the molecule can point almost anywhere in it. A\n");
  std::print(myfile, "# window it must thread costs a great deal, and none of that cost appears in the geometry\n");
  std::print(myfile, "# or in a single-site energy.\n");
  std::print(myfile, "\n");

  auto line = [&](const std::string &label, double value)
  {
    if (value == never)
    {
      std::print(myfile, "{:34} never\n", label);
      return;
    }
    std::print(myfile, "{:34} {:14.4f} [K] {:12.4f} [kJ/mol] {:10.2f} kT\n", label, asKelvin(value),
               asKJPerMol(value), value / kT);
  };

  // The dimensionality is read off at exactly the barrier, and a framework whose channels are alike in more
  // than one direction opens the second of them a hair above the first, which makes a two-dimensional
  // network look one-dimensional. Say what is open within a kT of the barrier beside it.
  auto dimensionality = [&](const EnergyBarrier &barrier)
  {
    std::print(myfile, "Dimensionality at the barrier:     {}\n", barrier.dimensionalityAtBarrier);

    const auto [within, cost] = barrier.dimensionalityWithin(kT);
    if (within > barrier.dimensionalityAtBarrier)
    {
      std::print(myfile, "Dimensionality within a kT of it:  {}, the last of them {:.4f} [K] above\n", within,
                 asKelvin(cost));
    }
  };

  std::print(myfile, "# From the orientational free energy\n");
  if (this->fromFreeEnergy.percolates)
  {
    line("Percolation barrier:", this->fromFreeEnergy.percolationBarrier);
    dimensionality(this->fromFreeEnergy);
    line("Deepest site anywhere:", this->fromFreeEnergy.deepestWell);
    line("Deepest site on the path:", this->fromFreeEnergy.deepestWellOnPath);
    line("Activation energy to leave:",
         this->fromFreeEnergy.percolationBarrier - this->fromFreeEnergy.deepestWellOnPath);
    std::print(myfile, "Boltzmann factor exp(-Ea/kT):      {:14.6e}\n",
               std::exp(-(this->fromFreeEnergy.percolationBarrier - this->fromFreeEnergy.deepestWellOnPath) / kT));
  }
  else
  {
    std::print(myfile, "No path runs through the crystal at any energy, which cannot happen while every\n");
    std::print(myfile, "grid point takes part and so points at a fault rather than at the structure.\n");
  }

  std::print(myfile, "\n# From the minimum over orientations, the zero-temperature bound\n");
  if (this->fromMinimumEnergy.percolates)
  {
    line("Percolation barrier:", this->fromMinimumEnergy.percolationBarrier);
    dimensionality(this->fromMinimumEnergy);
    line("Deepest site anywhere:", this->fromMinimumEnergy.deepestWell);
    line("Deepest site on the path:", this->fromMinimumEnergy.deepestWellOnPath);
    line("Activation energy to leave:",
         this->fromMinimumEnergy.percolationBarrier - this->fromMinimumEnergy.deepestWellOnPath);
  }

  std::print(myfile, "\n# What the orientation costs\n");
  line("Free energy above best-oriented:", this->orientationalPenalty);
  std::print(myfile, "# As a factor on the hopping rate:  {:14.6e}\n", std::exp(-this->orientationalPenalty / kT));

  std::print(myfile, "\n# The barrier to running in at least one, two, and three directions, from the free\n");
  std::print(myfile, "# energy. The sweep passes each on its way up, so they cost nothing beyond the first.\n");
  for (std::size_t dimension = 0; dimension < 3; ++dimension)
  {
    line(std::format("Barrier in {} direction(s) or more:", dimension + 1),
         this->fromFreeEnergy.barrierByDimension[dimension]);
  }

  myfile.close();
}
