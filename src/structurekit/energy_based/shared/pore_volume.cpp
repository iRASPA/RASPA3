module;

module energy_shared_pore_volume;

import std;

import uint3;
import double3;
import skspacegroupdatabase;
import crystal;
import pair_interactions;
import units;
import grid_percolation;
import grid_connected_components;
import energy_shared_linear_probe;
import energy_shared_energy_backend;
import energy_shared_molecular_energy_grid;
import energy_shared_electrostatic_potential_grid;


EnergyPoreVolume::EnergyPoreVolume() {}


EnergyPoreVolume::~EnergyPoreVolume() {}


void EnergyPoreVolume::run(const EnergyBackend &backend, const PairInteractions &interactions, const Crystal &framework,
                           const LinearProbe &probe, double isoValue, uint3 gridSize,
                           std::size_t numberOfOrientations, double temperature, double thresholdInKT,
                           bool useElectrostatics, double relativePrecision)
{
  MolecularField field = buildMolecularField(backend, interactions, framework, probe, gridSize, numberOfOrientations,
                                             temperature, useElectrostatics, relativePrecision);
  this->grid = field.grid;
  this->potential = field.potential;
  this->isoValue = isoValue;
  this->temperature = temperature;
  this->thresholdInKT = thresholdInKT;

  const std::size_t numberOfVoxels = this->grid.numberOfVoxels();
  if (numberOfVoxels == 0) return;

  std::chrono::steady_clock::time_point time_begin = std::chrono::steady_clock::now();

  const double total = static_cast<double>(numberOfVoxels);
  const double volume = framework.unitCell.volume;

  std::vector<float> openness(numberOfVoxels);
  for (std::size_t v = 0; v < numberOfVoxels; ++v) openness[v] = -this->grid.freeEnergy[v];
  const double level = -isoValue;

  GridComponents components = GridComponents::compute(gridSize, openness, level);

  std::size_t channelVoxels = 0;
  std::size_t pocketVoxels = 0;
  for (const GridPore &pore : components.pores)
  {
    if (pore.isChannel)
      channelVoxels += pore.numberOfVoxels;
    else
      pocketVoxels += pore.numberOfVoxels;
  }

  this->numberOfChannels = components.numberOfChannels;
  this->numberOfPockets = components.numberOfPockets;
  this->dimensionality = components.dimensionality;

  this->voidFraction = static_cast<double>(components.numberOfOpenVoxels) / total;
  this->accessibleFraction = static_cast<double>(channelVoxels) / total;
  this->inaccessibleFraction = static_cast<double>(pocketVoxels) / total;

  this->voidVolume = this->voidFraction * volume;
  this->accessibleVolume = this->accessibleFraction * volume;
  this->inaccessibleVolume = this->inaccessibleFraction * volume;

  // The same region, divided by what it costs to leave. The sweep takes in the whole of the field and not
  // only the region, since leaving is a matter of climbing out of the region for a while.
  GridPercolation swept = sweepPercolation(gridSize, openness, std::numeric_limits<float>::lowest(), true);

  const float noWayOut = std::numeric_limits<float>::lowest();
  const double kT = Units::KB * temperature;

  std::vector<float> escapeOfPore(components.pores.size(), noWayOut);
  for (std::size_t v = 0; v < numberOfVoxels; ++v)
  {
    const std::int32_t pore = components.voxelPore[v];
    if (pore < 0) continue;
    escapeOfPore[static_cast<std::size_t>(pore)] =
        std::max(escapeOfPore[static_cast<std::size_t>(pore)], swept.escapeOpenness[v]);
  }

  std::size_t reachableVoxels = 0;
  std::size_t sealedVoxels = 0;
  for (std::size_t poreIndex = 0; poreIndex < components.pores.size(); ++poreIndex)
  {
    const GridPore &pore = components.pores[poreIndex];
    if (pore.isChannel)
    {
      reachableVoxels += pore.numberOfVoxels;
      continue;
    }

    const float escape = escapeOfPore[poreIndex];
    const double barrier =
        (escape == noWayOut || kT <= 0.0)
            ? std::numeric_limits<double>::max()
            : (-static_cast<double>(escape) - (-pore.largestOpenness)) / kT;

    if (barrier > thresholdInKT)
    {
      sealedVoxels += pore.numberOfVoxels;
    }
    else
    {
      reachableVoxels += pore.numberOfVoxels;
      this->largestReachableBarrierInKT = std::max(this->largestReachableBarrierInKT, barrier);
    }
  }

  this->reachableFraction = static_cast<double>(reachableVoxels) / total;
  this->sealedFraction = static_cast<double>(sealedVoxels) / total;
  this->reachableVolume = this->reachableFraction * volume;
  this->sealedVolume = this->sealedFraction * volume;

  std::size_t leastVoxels = 0;
  for (std::size_t v = 0; v < numberOfVoxels; ++v)
  {
    if (static_cast<double>(this->grid.minimumEnergy[v]) <= isoValue) ++leastVoxels;
  }
  this->minimumEnergyVoidFraction = static_cast<double>(leastVoxels) / total;
  this->minimumEnergyVoidVolume = this->minimumEnergyVoidFraction * volume;

  // The thermodynamic reading. Everything takes part, the wall included, where it contributes nothing.
  double sum = 0.0;
  if (temperature > 0.0)
  {
    const double beta = 1.0 / (Units::KB * temperature);
    for (std::size_t v = 0; v < numberOfVoxels; ++v)
    {
      sum += std::exp(-beta * static_cast<double>(this->grid.freeEnergy[v]));
    }
  }
  this->boltzmannFraction = sum / total;
  this->boltzmannVolume = this->boltzmannFraction * volume;
  this->readsAsFraction = this->boltzmannFraction <= 1.0;

  double toGravimetric = (Units::Angstrom * Units::Angstrom * Units::Angstrom * 1.0e6) * Units::AvogadroConstant /
                         framework.mass;
  this->frameworkDensity = framework.mass / (volume * Units::Angstrom * Units::Angstrom * Units::Angstrom *
                                                     1.0e6 * Units::AvogadroConstant);

  this->gravimetricVoidVolume = this->voidVolume * toGravimetric;
  this->gravimetricAccessibleVolume = this->accessibleVolume * toGravimetric;
  this->gravimetricBoltzmannVolume = this->boltzmannVolume * toGravimetric;
  this->gravimetricReachableVolume = this->reachableVolume * toGravimetric;

  std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - time_begin;
  this->seconds = elapsed.count();

  double3 spacing = this->grid.spacing();
  const double voxelVolume = volume / total;

  std::ofstream myfile;
  myfile.open(framework.name + "." + probe.name + ".energy.av." + this->grid.backend + ".txt");
  std::print(myfile, "# Pore volume from the energy landscape\n");
  std::print(myfile, "# Crystal: {}\n", framework.name);
  std::print(myfile, "# Space-group HM-symbol: {}\n",
             SKSpaceGroupDataBase::spaceGroupData[framework.spaceGroupHallNumber].HMString());
  std::print(myfile, "# Number of framework atoms: {}\n", framework.atoms.size());
  std::print(myfile, "# Crystal volume: {} [Å³]\n", volume);
  std::print(myfile, "# Crystal mass: {} [g/mol]\n", framework.mass);
  std::print(myfile, "# Crystal density: {} [g/cm³]\n", this->frameworkDensity);
  std::print(myfile, "# Molecule: {}, {} sites, {:.4f} [Å] end to end\n", probe.name, probe.sites.size(),
             probe.length());
  std::print(myfile, "# Orientations sampled: {} over the {}\n", this->grid.numberOfOrientations,
             this->grid.overHemisphere ? "hemisphere, the molecule being the same end for end" : "whole sphere");
  std::print(myfile, "# Grid: {} x {} x {} points, spacing {:.5f} x {:.5f} x {:.5f} [Å], far face left out\n",
             gridSize.x, gridSize.y, gridSize.z, spacing.x, spacing.y, spacing.z);
  std::print(myfile, "# Temperature: {} [K]\n", temperature);
  std::print(myfile, "# Cutoff: {} [Å]\n", this->grid.cutOff);
  if (this->grid.chargesIncluded)
  {
    std::print(myfile, "# Electrostatics: Ewald, split at alpha = {:.5f} [1/Å], {} wave vectors\n",
               this->potential.alpha, this->potential.numberOfWaveVectors);
  }
  if (this->grid.chargesIgnored)
  {
    std::print(myfile, "#\n");
    std::print(myfile, "# WARNING: this molecule carries partial charges and they have not been acted on.\n");
  }
  std::print(myfile, "# Iso-value: {} [internal], {:.4f} [K]\n", this->isoValue,
             this->isoValue * Units::EnergyToKelvin);
  std::print(myfile, "# Timing ({}): {} [s] for the landscape\n", this->grid.backend, this->grid.seconds);
  std::print(myfile, "# Timing ({}): {} [s] for the counting\n", backend.name, this->seconds);
  std::print(myfile, "#\n");

  std::print(myfile, "# The geometric reading: the points at which the molecule is at least as well off as in\n");
  std::print(myfile, "# the gas, counted, and divided by whether the piece holding them runs through the\n");
  std::print(myfile, "# crystal. This is the number to set beside a clearance-grid or Voronoi table, and it\n");
  std::print(myfile, "# moves when the iso-value moves.\n");
  std::print(myfile, "Void volume:         fraction {}  {} [Å³]  {} [cm³/g]\n", this->voidFraction, this->voidVolume,
             this->gravimetricVoidVolume);
  std::print(myfile, "Accessible volume:   fraction {}  {} [Å³]  {} [cm³/g]\n", this->accessibleFraction,
             this->accessibleVolume, this->gravimetricAccessibleVolume);
  std::print(myfile, "Inaccessible volume: fraction {}  {} [Å³]  {} [cm³/g]\n", this->inaccessibleFraction,
             this->inaccessibleVolume, this->inaccessibleVolume * toGravimetric);
  std::print(myfile, "Best orientation:    fraction {}  {} [Å³]  {} [cm³/g]\n", this->minimumEnergyVoidFraction,
             this->minimumEnergyVoidVolume, this->minimumEnergyVoidVolume * toGravimetric);
  std::print(myfile, "# Grid points below the level: {} of {}, one point standing for {:.6f} [Å³]\n",
             components.numberOfOpenVoxels, numberOfVoxels, voxelVolume);
  std::print(myfile, "# Channels: {}, pockets: {}, of those {} holding a single grid point\n",
             this->numberOfChannels, this->numberOfPockets, components.numberOfSinglePointPockets);
  std::print(myfile, "# Directions the void runs in: {}\n", this->dimensionality);
  std::print(myfile, "\n");

  std::print(myfile, "# The same region divided again, by what it costs the molecule to leave a piece rather\n");
  std::print(myfile, "# than by whether a path out exists at the level. On a landscape this is the division\n");
  std::print(myfile, "# that means something, and on a framework with tight windows it is not the one above.\n");
  std::print(myfile, "# A cage the molecule enters by paying a few kT is counted inaccessible by connectivity\n");
  std::print(myfile, "# and is filled in any experiment; only a barrier it will not pay in the time it is\n");
  std::print(myfile, "# given really shuts it out.\n");
  std::print(myfile, "Reachable volume:    fraction {}  {} [Å³]  {} [cm³/g]\n", this->reachableFraction,
             this->reachableVolume, this->gravimetricReachableVolume);
  std::print(myfile, "Sealed volume:       fraction {}  {} [Å³]  {} [cm³/g]\n", this->sealedFraction,
             this->sealedVolume, this->sealedVolume * toGravimetric);
  std::print(myfile, "# Sealed above {} kT, which is {:.2f} [K] at this temperature.\n", this->thresholdInKT,
             this->thresholdInKT * kT * Units::EnergyToKelvin);
  std::print(myfile, "# The steepest climb anything counted reachable had to make was {:.3f} kT, so the\n",
             this->largestReachableBarrierInKT);
  std::print(myfile, "# threshold is {} the answer here.\n",
             (this->largestReachableBarrierInKT < 0.5 * this->thresholdInKT) ? "nowhere near deciding"
                                                                            : "close enough to be worth quoting with");
  std::print(myfile, "# The blocking-spheres report goes through the same division piece by piece.\n");
  std::print(myfile, "\n");

  std::print(myfile, "# The thermodynamic reading: the integral of exp(-A/kT) over the whole cell, the wall\n");
  std::print(myfile, "# included, where it contributes nothing. No level enters it and no point is in or out.\n");
  std::print(myfile, "# It is the volume the molecule would fill if it were spread over the framework in\n");
  std::print(myfile, "# proportion to how long it spends in each place, and it is what a Henry coefficient is\n");
  std::print(myfile, "# built from.\n");
  std::print(myfile, "Boltzmann volume:    fraction {}  {} [Å³]  {} [cm³/g]\n", this->boltzmannFraction,
             this->boltzmannVolume, this->gravimetricBoltzmannVolume);
  if (!this->readsAsFraction)
  {
    std::print(myfile, "# Above one, so this is not a fraction of anything. The molecule is bound, and a place\n");
    std::print(myfile, "# it sits at constantly is worth more than the room it takes up. Read it as an excess\n");
    std::print(myfile, "# Boltzmann factor; a void fraction is a thing to measure with helium.\n");
  }
  std::print(myfile, "\n");

  const double ratio = (this->voidFraction > 0.0) ? this->boltzmannFraction / this->voidFraction : 0.0;
  std::print(myfile, "# The two against each other: {:.4f} of the geometric reading.\n", ratio);
  std::print(myfile, "#\n");
  std::print(myfile, "# Above one, the framework's room is worth more to this molecule than its size: the\n");
  std::print(myfile, "# places it can be, it is held in. Below one, the room is there and the molecule has\n");
  std::print(myfile, "# little reason to use it. The two come together for a molecule barely held anywhere,\n");
  std::print(myfile, "# which is why a void volume is quoted for helium and why quoting one for anything that\n");
  std::print(myfile, "# adsorbs is a question about which of these two was meant.\n");

  myfile.close();
}
