module;

module energy_shared_pore_analysis;

import std;

import uint3;
import double3;
import skspacegroupdatabase;
import crystal;
import pair_interactions;
import units;
import grid_pore_size;
import grid_percolation;
import grid_connected_components;
import energy_shared_linear_probe;
import energy_shared_energy_backend;
import energy_shared_energy_barrier;
import energy_shared_molecular_energy_grid;
import energy_shared_electrostatic_potential_grid;


EnergyPoreAnalysis::EnergyPoreAnalysis() {}


EnergyPoreAnalysis::~EnergyPoreAnalysis() {}


void EnergyPoreAnalysis::run(const EnergyBackend &backend, const PairInteractions &interactions, const Crystal &framework,
                             const LinearProbe &probe, double isoValue, uint3 gridSize,
                             std::size_t numberOfOrientations, double temperature, bool useElectrostatics,
                             double relativePrecision)
{
  MolecularField field = buildMolecularField(backend, interactions, framework, probe, gridSize, numberOfOrientations,
                                             temperature, useElectrostatics, relativePrecision);
  this->grid = field.grid;
  this->potential = field.potential;
  this->isoValue = isoValue;
  this->temperature = temperature;

  const std::size_t numberOfVoxels = this->grid.numberOfVoxels();
  if (numberOfVoxels == 0) return;

  std::chrono::steady_clock::time_point time_begin = std::chrono::steady_clock::now();

  // The levels, which need nothing named.
  this->levels = EnergyBarrier::fromField(gridSize, this->grid.freeEnergy, temperature);
  this->minimumEnergyLevels = EnergyBarrier::fromField(gridSize, this->grid.minimumEnergy, temperature);
  this->orientationalPenalty = this->levels.percolationBarrier - this->minimumEnergyLevels.percolationBarrier;

  // The lengths, which need one. The landscape is turned over into an openness and the region below the
  // iso-value becomes the region above the negated level, as everywhere else here.
  std::vector<float> openness(numberOfVoxels);
  for (std::size_t v = 0; v < numberOfVoxels; ++v) openness[v] = -this->grid.freeEnergy[v];
  const double level = -isoValue;

  this->channels = GridComponents::compute(gridSize, openness, level);

  // How much room there is at each point of that region, which is what a diameter is read off. Outside the
  // region it is negative, so a floor of zero is exactly the region taking part and nothing else.
  std::vector<float> distance = distanceToIsosurface(gridSize, framework.unitCell, openness, level);
  GridPercolation swept = sweepPercolation(gridSize, distance, 0.0f);

  double widest = 0.0;
  for (std::size_t v = 0; v < numberOfVoxels; ++v) widest = std::max(widest, static_cast<double>(distance[v]));

  this->diameters.includedSphereDiameter = 2.0 * std::max(0.0, widest);
  this->diameters.numberOfVoidVoxels = swept.numberOfVoxels;
  this->diameters.percolates = swept.percolates;
  this->diameters.dimensionalityAtThreshold = swept.dimensionalityAtThreshold;
  this->diameters.seconds = swept.seconds;

  if (swept.percolates)
  {
    this->diameters.freeSphereDiameter = 2.0 * swept.percolationOpenness;
    this->diameters.includedAlongFreePathDiameter = 2.0 * swept.highestOpennessOnPath;
  }
  for (std::size_t dimension = 0; dimension < 3; ++dimension)
  {
    const double reach = swept.opennessByDimension[dimension];
    this->diameters.freeSphereDiameterByDimension[dimension] =
        (reach == std::numeric_limits<double>::lowest()) ? 0.0 : 2.0 * reach;
  }

  std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - time_begin;
  this->seconds = elapsed.count();

  const double kT = Units::KB * temperature;
  double3 spacing = this->grid.spacing();

  auto inKelvin = [](double energy) { return energy * Units::EnergyToKelvin; };

  std::ofstream myfile;
  myfile.open(framework.name + "." + probe.name + ".energy.res." + this->grid.backend + ".txt");
  std::print(myfile, "# Pore diameters from the energy landscape (Di, Df, Dif)\n");
  std::print(myfile, "# Crystal: {}\n", framework.name);
  std::print(myfile, "# Space-group HM-symbol: {}\n",
             SKSpaceGroupDataBase::spaceGroupData[framework.spaceGroupHallNumber].HMString());
  std::print(myfile, "# Number of framework atoms: {}\n", framework.atoms.size());
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
  std::print(myfile, "# Timing ({}): {} [s] for the landscape\n", this->grid.backend, this->grid.seconds);
  std::print(myfile, "# Timing ({}): {} [s] for the sweeps\n", backend.name, this->seconds);
  std::print(myfile, "#\n");

  std::print(myfile, "# The three diameters are levels of a field rather than lengths, and on an energy\n");
  std::print(myfile, "# landscape they stay levels: nothing has to be named to ask for them, and they come\n");
  std::print(myfile, "# out in Kelvin. These are the answers to quote.\n");
  std::print(myfile, "#\n");
  std::print(myfile, "# Every energy is given in K and in kT at {} K, a barrier being worth knowing chiefly\n",
             temperature);
  std::print(myfile, "# as a multiple of kT.\n");
  std::print(myfile, "#\n");
  std::print(myfile, "{:44} {:14.4f} K {:10.3f} kT\n", "Di, the deepest well:", inKelvin(this->levels.deepestWell),
             (kT > 0.0) ? this->levels.deepestWell / kT : 0.0);
  if (this->levels.percolates)
  {
    std::print(myfile, "{:44} {:14.4f} K {:10.3f} kT\n",
               "Df, the percolation barrier:", inKelvin(this->levels.percolationBarrier),
               (kT > 0.0) ? this->levels.percolationBarrier / kT : 0.0);
    std::print(myfile, "{:44} {:14.4f} K {:10.3f} kT\n",
               "Dif, the deepest well on that path:", inKelvin(this->levels.deepestWellOnPath),
               (kT > 0.0) ? this->levels.deepestWellOnPath / kT : 0.0);
    std::print(myfile, "{:44} {:14}\n", "Directions it runs in there:", this->levels.dimensionalityAtBarrier);

    const auto [within, cost] = this->levels.dimensionalityWithin(kT);
    if (within > this->levels.dimensionalityAtBarrier)
    {
      std::print(myfile,
                 "# Another {} open{} {:.4f} K above that, which is {:.4f} kT and nothing at this temperature,\n"
                 "# so read the network as {}-dimensional.\n",
                 (within - this->levels.dimensionalityAtBarrier == 1) ? "direction" : "two directions",
                 (within - this->levels.dimensionalityAtBarrier == 1) ? "s" : "", cost * Units::EnergyToKelvin,
                 cost / kT, within);
    }
  }
  else
  {
    std::print(myfile, "The landscape never runs through the crystal, at any energy.\n");
  }
  std::print(myfile, "{:44} {:14.4f} K {:10.3f} kT\n",
             "The same barrier, best orientation:", inKelvin(this->minimumEnergyLevels.percolationBarrier),
             (kT > 0.0) ? this->minimumEnergyLevels.percolationBarrier / kT : 0.0);
  std::print(myfile, "{:44} {:14.4f} K {:10.3f} kT\n",
             "What turning costs at the window:", inKelvin(this->orientationalPenalty),
             (kT > 0.0) ? this->orientationalPenalty / kT : 0.0);
  std::print(myfile, "\n");

  std::print(myfile, "# The barrier to running in at least one, two and three independent directions, which\n");
  std::print(myfile, "# the same pass gives on its way up.\n");
  for (int dimension = 1; dimension <= 3; ++dimension)
  {
    const double barrier = this->levels.barrierByDimension[static_cast<std::size_t>(dimension - 1)];
    if (barrier < std::numeric_limits<double>::max())
    {
      std::print(myfile, "Df in {} direction(s) or more:                {:14.4f} K {:10.3f} kT\n", dimension,
                 inKelvin(barrier), (kT > 0.0) ? barrier / kT : 0.0);
    }
    else
    {
      std::print(myfile, "Df in {} direction(s) or more:                          none\n", dimension);
    }
  }
  std::print(myfile, "\n");

  std::print(myfile, "# The same three as lengths, for setting beside a geometric table. A length needs a\n");
  std::print(myfile, "# boundary and the landscape has none, so these are measured against the contour at the\n");
  std::print(myfile, "# iso-value below, and they are statements about that level as much as about the\n");
  std::print(myfile, "# framework. Moving the level moves all three. The levels above do not move.\n");
  std::print(myfile, "#\n");
  std::print(myfile, "# The contour is found by walking out to where the field crosses the level and\n");
  std::print(myfile, "# interpolating between the two grid points that straddle it, so a distance is right to\n");
  std::print(myfile, "# about half a voxel and only a finer grid improves it.\n");
  std::print(myfile, "#\n");
  std::print(myfile, "# Iso-value: {} [internal], {:.4f} [K]\n", this->isoValue, inKelvin(this->isoValue));
  std::print(myfile, "# Grid points inside the contour: {} of {}\n", this->diameters.numberOfVoidVoxels,
             numberOfVoxels);
  std::print(myfile, "Di (largest included sphere):            {} [Å]\n", this->diameters.includedSphereDiameter);
  std::print(myfile, "Df (largest free sphere):                {} [Å]\n", this->diameters.freeSphereDiameter);
  std::print(myfile, "Dif (included sphere along free path):   {} [Å]\n",
             this->diameters.includedAlongFreePathDiameter);
  if (!this->diameters.percolates)
  {
    std::print(myfile, "At this level the region does not percolate: every piece closes inside one cell.\n");
    std::print(myfile, "Note that the barrier above says the molecule gets through all the same, by going\n");
    std::print(myfile, "above the level. That is the difference between the two halves of this file.\n");
  }
  for (int dimension = 1; dimension <= 3; ++dimension)
  {
    const double diameter = this->diameters.freeSphereDiameterByDimension[static_cast<std::size_t>(dimension - 1)];
    if (diameter > 0.0)
      std::print(myfile, "Df in {} direction(s) or more:            {} [Å]\n", dimension, diameter);
    else
      std::print(myfile, "Df in {} direction(s) or more:            none\n", dimension);
  }
  myfile.close();

  std::ofstream channelFile;
  channelFile.open(framework.name + "." + probe.name + ".energy.chan." + this->grid.backend + ".txt");
  std::print(channelFile, "# Channels and pockets of the energy landscape\n");
  std::print(channelFile, "# Crystal: {}\n", framework.name);
  std::print(channelFile, "# Molecule: {}, {} sites\n", probe.name, probe.sites.size());
  std::print(channelFile, "# Grid: {} x {} x {} points, spacing {:.5f} x {:.5f} x {:.5f} [Å]\n", gridSize.x,
             gridSize.y, gridSize.z, spacing.x, spacing.y, spacing.z);
  std::print(channelFile, "# Temperature: {} [K]\n", temperature);
  std::print(channelFile, "# Iso-value: {} [internal], {:.4f} [K]\n", this->isoValue, inKelvin(this->isoValue));
  std::print(channelFile, "# Grid points below it: {} of {}\n", this->channels.numberOfOpenVoxels, numberOfVoxels);
  std::print(channelFile, "# Timing: {} [s]\n", this->channels.seconds);
  std::print(channelFile, "# A pore's dimensionality is the number of independent lattice translations that\n");
  std::print(channelFile, "# bring it back to itself: none for a pocket, and one, two or three for a channel,\n");
  std::print(channelFile, "# a layer, and a network.\n");
  std::print(channelFile, "#\n");
  std::print(channelFile, "# Which pieces there are is a property of the level, and unlike the geometric route\n");
  std::print(channelFile, "# the level was chosen. Raising it takes in the windows the molecule crosses by\n");
  std::print(channelFile, "# paying a few kT, and pockets join up into channels as it goes. Whether a piece\n");
  std::print(channelFile, "# here is really sealed is a question about barriers, which the blocking-spheres\n");
  std::print(channelFile, "# report answers; this one only says what the level does.\n");
  std::print(channelFile, "Number of channels: {}\n", this->channels.numberOfChannels);
  std::print(channelFile, "Number of pockets:  {}\n", this->channels.numberOfPockets);
  std::print(channelFile, "  of those, pockets holding a single grid point: {}\n",
             this->channels.numberOfSinglePointPockets);
  std::print(channelFile, "Pore system dimensionality: {}\n", this->channels.dimensionality);
  for (std::size_t i = 0; i < this->channels.pores.size(); ++i)
  {
    const GridPore &pore = this->channels.pores[i];
    std::print(channelFile, "  pore {}: {} dimensionality={} points={} deepest={:.4f} [K]\n", i,
               pore.isChannel ? "channel" : "pocket", pore.dimensionality, pore.numberOfVoxels,
               inKelvin(-pore.largestOpenness));
  }
  channelFile.close();
}
