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

  double sum{};
  double sum_of_squares{};
  for (std::size_t i = 0; i < number_of_iterations; ++i)
  {
    double surface_area{};
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

      surface_area +=
          (counted / total_trials) * 4.0 * std::numbers::pi * equilibrium_distance * equilibrium_distance;
    }

    // One independent reading of the total area: keep sum and sum of squares for the mean and its error.
    sum += surface_area;
    sum_of_squares += surface_area * surface_area;
  }

  std::chrono::duration<double> timing = std::chrono::steady_clock::now() - time_begin;

  this->seconds = timing.count();
  double n = static_cast<double>(number_of_iterations);
  this->surfaceArea = n > 0.0 ? sum / n : 0.0;

  // Student-t critical values for a two-sided 95% interval (df = 1..20); beyond that the normal limit.
  constexpr std::array<double, 21> student_t_95{
      0.0,   12.71, 4.303, 3.182, 2.776, 2.571, 2.447, 2.365, 2.306, 2.262, 2.228,
      2.201, 2.179, 2.160, 2.145, 2.131, 2.120, 2.110, 2.101, 2.093, 2.086};
  this->surfaceAreaError = 0.0;
  if (number_of_iterations >= 3)
  {
    double standard_error = std::sqrt((sum_of_squares - sum * sum / n) / (n * (n - 1.0)));
    std::size_t degrees_of_freedom = number_of_iterations - 1;
    double t = degrees_of_freedom < student_t_95.size() ? student_t_95[degrees_of_freedom] : 1.959963984540054;
    this->surfaceAreaError = t * standard_error;
  }

  std::ofstream myfile;
  myfile.open(structure.name + ".mc.sa.cpu.txt");
  std::print(myfile, "# Surface area using Monte Carlo-based method\n");
  structure.writeHeader(myfile);
  probe.writeHeader(myfile);
  std::print(myfile, "# Number of iterations: {}\n", number_of_iterations);
  std::print(myfile, "# Number of inner-steps (sample points per atom): {}\n", number_of_inner_steps);
  std::print(myfile, "# CPU Timing: {} [s]\n", this->seconds);
  std::print(myfile, "# The area is the mean over independent passes; the error beside it is the\n");
  std::print(myfile, "# half-width of the 95% confidence interval (Student-t times the standard\n");
  std::print(myfile, "# error of the mean from the sum and sum of squares of the passes).\n");
  std::print(myfile, "{} +/- {} [Å²]\n", this->surfaceArea, this->surfaceAreaError);
  std::print(myfile, "{} +/- {} [m²/cm³]\n", 1.0e4 * this->surfaceArea / structure.unitCell.volume,
             1.0e4 * this->surfaceAreaError / structure.unitCell.volume);
  std::print(myfile, "{} +/- {} [m²/g]\n", this->surfaceArea * structure.gravimetricFactor(),
             this->surfaceAreaError * structure.gravimetricFactor());
  myfile.close();
}
