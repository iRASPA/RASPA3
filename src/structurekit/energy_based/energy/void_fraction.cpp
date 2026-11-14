module;

module energy_void_fraction;

import std;

import opencl;
import float4;
import double2;
import double4;
import randomnumbers;
import skspacegroupdatabase;
import pair_interactions;
import crystal;
import unit_cell;
import units;
import blocking_spheres;
import energy_shared_blocking_mask;
import energy_shared_probe_energy_grid;
#if !(defined(__has_include) && __has_include(<mdspan>))
//import mdspan;
#endif

void EnergyVoidFraction::run(const PairInteractions &interactions, const Crystal &framework, std::string probePseudoAtom, 
                             std::optional<std::size_t> numberOfIterations, std::optional<std::size_t> numberOfInnersteps,
                             double temperature, std::span<const BlockingSphere> blockingSpheres,
                             double blockedEnergyPerAngstrom, double ceiling)
{
  if (temperature <= 0.0)
  {
    throw std::runtime_error("EnergyVoidFraction: a Boltzmann average needs a temperature above zero\n");
  }
  // Nought for either is the usual figure, as on the well-surface and BET routes: a ramp of nought and a
  // ceiling of nought make a pocket worth nothing at all, which is open void, and nobody who hands over
  // spheres is asking for that.
  if (blockedEnergyPerAngstrom == 0.0)
  {
    blockedEnergyPerAngstrom = blockedEnergyPerAngstromInKelvin * Units::KelvinToEnergy;
  }
  if (ceiling == 0.0) ceiling = probeEnergyCeilingInKelvin * Units::KelvinToEnergy;
  const double beta = 1.0 / (Units::KB * temperature);
  RandomNumber random{std::nullopt};
  std::chrono::steady_clock::time_point time_begin, time_end;

  std::optional<std::size_t> probe_type = interactions.findType(probePseudoAtom);
  if (!probe_type.has_value())
  {
    throw std::runtime_error(std::format("MC_SurfaceArea: Unknown probe-atom type\n"));
  }

  std::size_t number_of_iterations = numberOfIterations.value_or(1000);
  std::size_t number_of_inner_steps = numberOfInnersteps.value_or(1000);

  time_begin = std::chrono::steady_clock::now();

  // the fractional positions of the atoms, from the Cartesian ones when the crystal does not carry them (one built
  // by hand: indexing its empty 'fractionalPositions' by the atoms read past it)
  std::vector<double3> fractionalPositions = framework.fractionalPositionsOfAtoms();
  if (fractionalPositions.size() != framework.size())
  {
    fractionalPositions.resize(framework.size());
    for (std::size_t i = 0; i < framework.size(); ++i)
    {
      fractionalPositions[i] = double3::fract(framework.unitCell.inverseCell * framework.atoms[i].position);
    }
  }

  double cutoff = interactions.cutOffVDW;
  int3 numberOfReplicas = framework.unitCell.smallestNumberOfUnitCellsForMinimumImagesConvention(cutoff);
  UnitCell simulation_box = framework.unitCell.scaled(numberOfReplicas);
  double3x3 simulation_cell = simulation_box.cell;

  double sumBoltzmannWeight{};
  std::size_t count{};
  for (std::size_t l = 0; l < number_of_iterations; ++l)
  {
    for (std::size_t m = 0; m < number_of_inner_steps; ++m)
    {
      double3 random_position = {random.uniform(), random.uniform(), random.uniform()};

      // A point in a pocket the probe cannot reach has the pocket's energy, as it has on the grids, and what the
      // atoms round it would have made of it is not asked.
      const double pocket = blockingSpheres.empty()
                                ? 1.0
                                : blockingSphereDistance(random_position, framework.unitCell, blockingSpheres);
      if (pocket < 0.0)
      {
        sumBoltzmannWeight += std::exp(-beta * std::min(-pocket * blockedEnergyPerAngstrom, ceiling));
        count++;
        continue;
      }

      double energy = 0.0;
      for (std::size_t atomIndex = 0; atomIndex < framework.size(); ++atomIndex)
      {
        double3 atomPosition = fractionalPositions[atomIndex];
        std::size_t typeB = framework.atoms[atomIndex].type;
        const PairParameters &pair = interactions(probe_type.value(), typeB);
        double epsilonTimesFour = 4.0 * pair.strengthParameter;
        double sigmaSquared = pair.sizeParameter * pair.sizeParameter;
        for (std::size_t i = 0; i < static_cast<std::size_t>(numberOfReplicas.x); i++)
        {
          for (std::size_t j = 0; j < static_cast<std::size_t>(numberOfReplicas.y); j++)
          {
            for (std::size_t k = 0; k < static_cast<std::size_t>(numberOfReplicas.z); k++)
            {
              double3 translation(static_cast<double>(i), static_cast<double>(j), static_cast<double>(k));
              double3 ds = (random_position - atomPosition + translation) / numberOfReplicas;

              ds.x -= std::rint(ds.x);
              ds.y -= std::rint(ds.y);
              ds.z -= std::rint(ds.z);

              double3 dr = simulation_cell * ds;

              // A point on an atom's centre has no pair energy; held off zero it is buried far above any
              // weight, as in the grid routes.
              double rr = std::max(double3::dot(dr, dr), 1.0e-6);
              if (rr < cutoff * cutoff)
              {
                // Shifted Lennard-Jones, as everywhere else here: the pair's own mixed size and strength,
                // less whatever the truncation was shifted by.
                double temp = rr / sigmaSquared;
                double ratio3 = 1.0 / (temp * temp * temp);

                energy += epsilonTimesFour * (ratio3 * (ratio3 - 1.0)) - pair.shift;
              }
            }
          }
        }
      }

      sumBoltzmannWeight += std::exp(-beta * energy);
      count++;
    }
  }

  time_end = std::chrono::steady_clock::now();

  std::chrono::duration<double> timing = time_end - time_begin;

  this->numberOfInsertions = count;
  this->boltzmannAverage = count > 0 ? sumBoltzmannWeight / static_cast<double>(count) : 0.0;

  std::ofstream myfile;
  myfile.open(framework.name + ".energy.vf.cpu.txt");
  std::print(myfile, "# Void-fraction using energy-based method\n");
  std::print(myfile, "# Crystal: {}\n", framework.name);
  std::print(myfile, "# Space-group Hall-number: {}\n", framework.spaceGroupHallNumber);
  std::print(myfile, "# Space-group Hall-symbol: {}\n", SKSpaceGroupDataBase::spaceGroupData[framework.spaceGroupHallNumber].HallString());
  std::print(myfile, "# Space-group HM-symbol: {}\n", SKSpaceGroupDataBase::spaceGroupData[framework.spaceGroupHallNumber].HMString());
  std::print(myfile, "# Space-group IT number: {}\n", SKSpaceGroupDataBase::spaceGroupData[framework.spaceGroupHallNumber].number());
  std::print(myfile, "# Number of framework atoms: {}\n", framework.atoms.size());
  std::print(myfile, "# Crystal volume: {} [Å³]\n", framework.unitCell.volume);
  std::print(myfile, "# Crystal mass: {} [g/mol]\n", framework.mass);
  std::print(myfile, "# Crystal density: {} [kg/m³]\n", 1e-3 * framework.mass /
      (framework.unitCell.volume * Units::Angstrom * Units::Angstrom * Units::Angstrom * Units::AvogadroConstant));
  std::print(myfile, "# Probe atom: {} strength-value/kʙ: {} [K] size-value: {} [Å]\n", probePseudoAtom, 
      interactions[probe_type.value()].strengthParameter * Units::EnergyToKelvin, 
      interactions[probe_type.value()].sizeParameter);
  std::print(myfile, "# Cutoff: {} Å\n", cutoff);
  std::print(myfile, "# Number of unit cells: {}x{}x{}\n", numberOfReplicas.x, numberOfReplicas.y, numberOfReplicas.z);
  std::print(myfile, "# Temperature: {} [K]\n", temperature);
  std::print(myfile, "# Blocking spheres: {}\n", blockingSpheres.size());
  std::print(myfile, "# CPU Timing: {} [s]\n", timing.count());
  myfile << this->boltzmannAverage << std::endl;
  myfile.close();
}
