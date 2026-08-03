module;

module energy_shared_blocking_spheres;

import std;

import int3;
import uint3;
import double3;
import unit_cell;
import skspacegroupdatabase;
import crystal;
import pair_interactions;
import units;
import grid_pore_size;
import grid_percolation;
import grid_connected_components;
import grid_blocking_cover;
import energy_shared_linear_probe;
import energy_shared_energy_backend;
import energy_shared_molecular_energy_grid;
import energy_shared_electrostatic_potential_grid;


EnergyBlockingSpheres::EnergyBlockingSpheres() {}


EnergyBlockingSpheres::~EnergyBlockingSpheres() {}


void EnergyBlockingSpheres::run(const EnergyBackend &backend, const PairInteractions &interactions,
                                const Crystal &framework, const LinearProbe &probe, double isoValue,
                                uint3 gridSize, std::size_t numberOfOrientations, double temperature,
                                double thresholdInKT, bool useElectrostatics, double relativePrecision)
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

  // Everything downstream reads a field as an openness, larger meaning more open, so the landscape is turned
  // over once here and the level with it.
  std::vector<float> openness(numberOfVoxels);
  for (std::size_t v = 0; v < numberOfVoxels; ++v) openness[v] = -this->grid.freeEnergy[v];
  const double level = -isoValue;

  // The pieces of the region the molecule can occupy.
  GridComponents components = GridComponents::compute(gridSize, openness, level);

  // The way out of each of them. The sweep has to take in the whole of the field and not only the region,
  // since leaving is exactly a matter of climbing out of the region for a while, and a sweep floored at the
  // level would have no route to offer anything.
  GridPercolation sweep =
      sweepPercolation(gridSize, openness, std::numeric_limits<float>::lowest(), true);

  const float noWayOut = std::numeric_limits<float>::lowest();
  const double never = std::numeric_limits<double>::max();
  const double kT = Units::KB * temperature;

  // Every point of one piece leaves by the same route, the pieces only growing as the level falls, so this
  // is a reading of one number per piece rather than an average of many.
  std::vector<float> escapeOfPore(components.pores.size(), noWayOut);
  for (std::size_t v = 0; v < numberOfVoxels; ++v)
  {
    const std::int32_t pore = components.voxelPore[v];
    if (pore < 0) continue;
    escapeOfPore[static_cast<std::size_t>(pore)] =
        std::max(escapeOfPore[static_cast<std::size_t>(pore)], sweep.escapeOpenness[v]);
  }

  const double voxelVolume = framework.unitCell.volume / static_cast<double>(numberOfVoxels);

  std::vector<std::uint8_t> needsCover(components.pores.size(), 0);
  std::size_t blockedVoxels = 0;
  std::size_t leakyVoxels = 0;
  std::size_t voidVoxels = 0;

  for (std::size_t poreIndex = 0; poreIndex < components.pores.size(); ++poreIndex)
  {
    const GridPore &pore = components.pores[poreIndex];
    voidVoxels += pore.numberOfVoxels;
    if (pore.isChannel) continue;

    EnergyCavity cavity;
    cavity.pore = poreIndex;
    cavity.numberOfVoxels = pore.numberOfVoxels;
    cavity.volume = static_cast<double>(pore.numberOfVoxels) * voxelVolume;
    cavity.deepestEnergy = -pore.largestOpenness;
    cavity.centreFractional = pore.centroidFractional;

    const float escape = escapeOfPore[poreIndex];
    if (escape == noWayOut)
    {
      cavity.escapeEnergy = never;
      cavity.escapeBarrier = never;
      cavity.barrierInKT = never;
      cavity.blocked = true;
    }
    else
    {
      cavity.escapeEnergy = -static_cast<double>(escape);
      cavity.escapeBarrier = cavity.escapeEnergy - cavity.deepestEnergy;
      cavity.barrierInKT = (kT > 0.0) ? cavity.escapeBarrier / kT : never;
      cavity.blocked = cavity.barrierInKT > thresholdInKT;
    }

    if (cavity.blocked)
    {
      needsCover[poreIndex] = 1;
      ++this->numberOfBlockedCavities;
      blockedVoxels += pore.numberOfVoxels;
    }
    else
    {
      ++this->numberOfLeakyCavities;
      leakyVoxels += pore.numberOfVoxels;
    }

    this->cavities.push_back(cavity);
  }

  this->numberOfChannels = components.numberOfChannels;
  this->voidFraction = static_cast<double>(voidVoxels) / static_cast<double>(numberOfVoxels);
  this->blockedFraction = static_cast<double>(blockedVoxels) / static_cast<double>(numberOfVoxels);
  this->leakyFraction = static_cast<double>(leakyVoxels) / static_cast<double>(numberOfVoxels);

  // Roomiest first, so the report reads from the cavities that hold something down to the ones that barely
  // exist.
  std::ranges::sort(this->cavities, [](const EnergyCavity &a, const EnergyCavity &b)
                    { return a.numberOfVoxels > b.numberOfVoxels; });

  // Where the molecule is entitled to be: the channels, and the cavities it gets out of on its own. A sphere
  // that reached into either would block room the framework really offers.
  std::vector<std::uint8_t> reachable(numberOfVoxels, 0);
  for (std::size_t v = 0; v < numberOfVoxels; ++v)
  {
    const std::int32_t pore = components.voxelPore[v];
    if (pore < 0) continue;
    if (needsCover[static_cast<std::size_t>(pore)] == 0) reachable[v] = 1;
  }

  std::chrono::duration<double> sweepElapsed = std::chrono::steady_clock::now() - time_begin;
  this->sweepSeconds = sweepElapsed.count();

  // A sphere is grown from the middle of a cavity outwards, and how deep inside it a point lies is what
  // orders them. The distance to the contour is the energy landscape's answer to that question.
  std::vector<float> distance = distanceToIsosurface(gridSize, framework.unitCell, openness, level);

  GridBlockingCover cover = coverPockets(gridSize, framework.unitCell, components.voxelPore, components.pores,
                                     needsCover, distance, reachable);

  this->spheres = cover.spheres;
  this->numberOfClippedSpheres = cover.numberOfClippedSpheres;
  this->numberOfRefusedPoints = cover.numberOfRefusedPoints;
  this->coverSeconds = cover.seconds;

  // What a simulation reads: nothing but the numbers, a comment being more than the format allows. The
  // molecule is in the name because unlike the geometric answer this one is a property of the pair.
  std::ofstream blockFile;
  blockFile.open(framework.name + "." + probe.name + ".block");
  std::print(blockFile, "{}\n", this->spheres.size());
  for (const GridBlockingSphere &sphere : this->spheres)
  {
    std::print(blockFile, "{} {} {} {}\n", sphere.centreFractional.x, sphere.centreFractional.y,
               sphere.centreFractional.z, sphere.radius);
  }
  blockFile.close();

  double3 spacing = this->grid.spacing();

  std::ofstream report;
  report.open(framework.name + "." + probe.name + ".energy.block." + this->grid.backend + ".txt");
  std::print(report, "# Blocking spheres (energy landscape)\n");
  std::print(report, "# Crystal: {}\n", framework.name);
  std::print(report, "# Space-group HM-symbol: {}\n",
             SKSpaceGroupDataBase::spaceGroupData[framework.spaceGroupHallNumber].HMString());
  std::print(report, "# Number of framework atoms: {}\n", framework.atoms.size());
  std::print(report, "# Crystal volume: {} [Å³]\n", framework.unitCell.volume);
  std::print(report, "# Molecule: {}, {} sites, {:.4f} [Å] end to end\n", probe.name, probe.sites.size(),
             probe.length());
  std::print(report, "# Orientations sampled: {} over the {}\n", this->grid.numberOfOrientations,
             this->grid.overHemisphere ? "hemisphere, the molecule being the same end for end" : "whole sphere");
  std::print(report, "# Grid: {} x {} x {} points, spacing {:.5f} x {:.5f} x {:.5f} [Å], far face left out\n",
             gridSize.x, gridSize.y, gridSize.z, spacing.x, spacing.y, spacing.z);
  std::print(report, "# Temperature: {} [K]\n", temperature);
  std::print(report, "# Cutoff: {} [Å]\n", this->grid.cutOff);
  if (this->grid.chargesIncluded)
  {
    std::print(report, "# Electrostatics: Ewald, split at alpha = {:.5f} [1/Å], {} wave vectors\n",
               this->potential.alpha, this->potential.numberOfWaveVectors);
  }
  if (this->grid.chargesIgnored)
  {
    std::print(report, "#\n");
    std::print(report, "# WARNING: this molecule carries partial charges and they have not been acted on.\n");
  }
  std::print(report, "# Iso-value: {} [internal], {:.4f} [K]\n", this->isoValue,
             this->isoValue * Units::EnergyToKelvin);
  std::print(report, "# Blocked above: {} kT, which is {:.2f} [K] at this temperature\n", this->thresholdInKT,
             this->thresholdInKT * kT * Units::EnergyToKelvin);
  std::print(report, "# Timing ({}): {} [s] for the landscape\n", this->grid.backend, this->grid.seconds);
  std::print(report, "# Timing ({}): {} [s] for the pieces and their barriers\n", backend.name, this->sweepSeconds);
  std::print(report, "# Timing ({}): {} [s] for the spheres\n", backend.name, this->coverSeconds);
  std::print(report, "# Spheres written to {}.{}.block: {}\n", framework.name, probe.name, this->spheres.size());
  std::print(report, "#\n");

  std::print(report, "{:44} {:14}\n", "Pieces that run through the crystal:", this->numberOfChannels);
  std::print(report, "{:44} {:14}\n", "Pieces that do not, and are blocked:", this->numberOfBlockedCavities);
  std::print(report, "{:44} {:14}\n", "Pieces that do not, and are not:", this->numberOfLeakyCavities);
  std::print(report, "{:44} {:14.8f}\n", "Of the cell, void at this level:", this->voidFraction);
  std::print(report, "{:44} {:14.8f}\n", "Of the cell, blocked:", this->blockedFraction);
  std::print(report, "{:44} {:14.8f}\n", "Of the cell, sealed but reachable anyway:", this->leakyFraction);
  std::print(report, "{:44} {:14}\n", "Spheres stopped by reachable room:", this->numberOfClippedSpheres);
  std::print(report, "{:44} {:14}\n", "Points too near reachable room to cover:", this->numberOfRefusedPoints);
  std::print(report, "\n");

  std::print(report, "# The last of the three counts is what this route has to say that the geometric one has\n");
  std::print(report, "# not. Those pieces do not run anywhere, so a hard sphere put into one could never get\n");
  std::print(report, "# out and the geometric route blocks them. A molecule gets out of them, in this case in\n");
  std::print(report, "# fewer than exp({}) attempts, so blocking them would throw away uptake the framework\n",
             this->thresholdInKT);
  std::print(report, "# really has. They are left open and listed below with the rest.\n");
  std::print(report, "#\n");
  std::print(report, "# How much of that is a real difference and how much is the threshold can be read off\n");
  std::print(report, "# the barriers themselves, which is why they are all here rather than only the verdict.\n");
  std::print(report, "# A framework whose barriers are hundreds of kT and tens of kT and nothing between is\n");
  std::print(report, "# one where the threshold decided nothing; a framework with a barrier sitting on it is\n");
  std::print(report, "# one where the answer should be quoted with the threshold beside it.\n");
  std::print(report, "#\n");

  if (!this->cavities.empty())
  {
    std::print(report, "# Every piece that does not run through the crystal. The barrier is the climb from the\n");
    std::print(report, "# bottom of the piece to the lowest energy at which a path runs from it to the outside.\n");
    std::print(report, "# column 1-3: the centre of the piece, fractional\n");
    std::print(report, "# column 4: the volume it holds [Å³]\n");
    std::print(report, "# column 5: grid points in it\n");
    std::print(report, "# column 6: the deepest energy in it [K]\n");
    std::print(report, "# column 7: the energy it has to reach to leave [K]\n");
    std::print(report, "# column 8: the climb between the two [K]\n");
    std::print(report, "# column 9: the same, in kT\n");
    std::print(report, "#       s_a         s_b         s_c   volume [Å³]  points   deepest [K]     escape [K]"
                       "    barrier [K]  barrier [kT]  verdict\n");
    for (const EnergyCavity &cavity : this->cavities)
    {
      const bool sealed = cavity.escapeBarrier >= never;
      std::print(report, "Piece: {:11.7f} {:11.7f} {:11.7f} {:12.5f} {:7} {:14.4f} ", cavity.centreFractional.x,
                 cavity.centreFractional.y, cavity.centreFractional.z, cavity.volume, cavity.numberOfVoxels,
                 cavity.deepestEnergy * Units::EnergyToKelvin);
      if (sealed)
      {
        std::print(report, "{:>14} {:>14} {:>13}  ", "no way out", "-", "-");
      }
      else
      {
        std::print(report, "{:14.4f} {:14.4f} {:13.4f}  ", cavity.escapeEnergy * Units::EnergyToKelvin,
                   cavity.escapeBarrier * Units::EnergyToKelvin, cavity.barrierInKT);
      }
      std::print(report, "{}\n", cavity.blocked ? "blocked" : "gets out");
    }
    std::print(report, "\n");

    // What the threshold is worth: the same division made again at a spread of thresholds, so that a reader
    // can see at once whether the answer turned on it.
    std::print(report, "# The same division at other thresholds. A barrier holds a molecule for about\n");
    std::print(report, "# exp(barrier/kT) attempts at a lattice frequency, so the row to read is the one whose\n");
    std::print(report, "# time is the time the framework was given to fill.\n");
    std::print(report, "#  kT      time held      pieces blocked   of the cell\n");
    const std::array<double, 6> thresholds{10.0, 20.0, 25.0, 30.0, 35.0, 40.0};
    for (double threshold : thresholds)
    {
      std::size_t blocked = 0;
      std::size_t voxels = 0;
      for (const EnergyCavity &cavity : this->cavities)
      {
        if (!(cavity.barrierInKT > threshold)) continue;
        ++blocked;
        voxels += cavity.numberOfVoxels;
      }

      // The time a barrier of this many kT holds for, at a lattice frequency of 1e12 per second.
      const double seconds = std::exp(threshold) / 1.0e12;
      std::print(report, "{:5.1f} {:14.3e} s {:14} {:14.8f}\n", threshold, seconds, blocked,
                 static_cast<double>(voxels) / static_cast<double>(numberOfVoxels));
    }
    std::print(report, "\n");
  }

  if (!this->spheres.empty())
  {
    std::print(report, "# The spheres, which are a covering of whatever the blocked pieces turned out to be.\n");
    std::print(report, "# Each is grown from the middle of its piece until it reaches the furthest point of it,\n");
    std::print(report, "# and cut back at the nearest point the molecule can reach on its own.\n");
    std::print(report, "# column 1-3: centre, fractional\n");
    std::print(report, "# column 4: radius [Å]\n");
    std::print(report, "# column 5: which piece\n");
    std::print(report, "# column 6: grid points of that piece this sphere covered\n");
    std::print(report, "# column 7: what stopped it growing\n");
    std::print(report, "#       s_a         s_b         s_c   radius [Å]   piece    points  stopped by\n");
    for (const GridBlockingSphere &sphere : this->spheres)
    {
      std::print(report, "Sphere: {:11.7f} {:11.7f} {:11.7f} {:11.5f} {:7} {:9}  {}\n", sphere.centreFractional.x,
                 sphere.centreFractional.y, sphere.centreFractional.z, sphere.radius, sphere.pocket,
                 sphere.voxelsCovered, sphere.clipped ? "reachable room" : "the piece");
    }
  }

  report.close();
}
