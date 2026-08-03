module;

module mc_surface_area;

import std;

import double3;
import randomnumbers;
import sampled_structure;

void MC_SurfaceArea::run(const SampledStructure &structure, const SampledProbe &probe,
                         std::optional<std::size_t> numberOfIterations,
                         std::optional<std::size_t> numberOfInnerSteps)
{
  RandomNumber random{std::nullopt};

  std::chrono::steady_clock::time_point time_begin = std::chrono::steady_clock::now();

  std::size_t number_of_iterations = numberOfIterations.value_or(100);
  std::size_t number_of_inner_steps = numberOfInnerSteps.value_or(1000);

  double accumulated_surface_area{};
  for (std::size_t i = 0; i < number_of_iterations; ++i)
  {
    for (std::size_t atom_index = 0; atom_index < structure.size(); ++atom_index)
    {
      double equilibrium_distance = structure.radii[atom_index];

      double total_trials{};
      double counted{};
      for (std::size_t j = 0; j < number_of_inner_steps; ++j)
      {
        double3 vec = random.randomVectorOnUnitSphere();

        double3 position = structure.positions[atom_index] + equilibrium_distance * vec;
        if (!structure.overlaps(position, atom_index))
        {
          counted += 1.0;
        }

        total_trials += 1.0;
      }

      accumulated_surface_area +=
          (counted / total_trials) * 4.0 * std::numbers::pi * equilibrium_distance * equilibrium_distance;
    }
  }

  std::chrono::duration<double> timing = std::chrono::steady_clock::now() - time_begin;

  this->seconds = timing.count();
  this->surfaceArea = accumulated_surface_area / static_cast<double>(number_of_iterations);

  std::ofstream myfile;
  myfile.open(structure.name + ".mc.sa.cpu.txt");
  std::print(myfile, "# Surface area using Monte Carlo-based method\n");
  structure.writeHeader(myfile);
  probe.writeHeader(myfile);
  std::print(myfile, "# Number of iterations: {}\n", number_of_iterations);
  std::print(myfile, "# Number of inner-steps (sample points per atom): {}\n", number_of_inner_steps);
  std::print(myfile, "# CPU Timing: {} [s]\n", this->seconds);
  std::print(myfile, "{} [Å²]\n", this->surfaceArea);
  std::print(myfile, "{} [m²/cm³]\n", 1.0e4 * this->surfaceArea / structure.unitCell.volume);
  std::print(myfile, "{} [m²/g]\n", this->surfaceArea * structure.gravimetricFactor());
  myfile.close();
}
