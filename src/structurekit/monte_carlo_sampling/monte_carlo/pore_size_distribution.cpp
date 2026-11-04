module;

module mc_pore_size_distribution;

import std;

import double3;
import double3x3;
import randomnumbers;
import sampled_structure;

void MC_PoreSizeDistribution::run(const SampledStructure &structure, std::optional<std::size_t> numberOfIterations,
                                  std::optional<std::size_t> numberOfInnerSteps, std::optional<double> maximumRange)
{
  RandomNumber random{std::nullopt};

  std::size_t number_of_iterations = numberOfIterations.value_or(10000);
  std::size_t number_of_inner_steps = numberOfInnerSteps.value_or(10000);

  std::chrono::steady_clock::time_point time_begin = std::chrono::steady_clock::now();

  double delta_r = maximumRange.value_or(10.0) / static_cast<double>(numberOfBins);

  // A second run is a run of its own and does not add to the first.
  std::ranges::fill(histogram, 0.0);
  std::ranges::fill(histogram_cummulative, 0.0);

  for (std::size_t i = 0; i < number_of_iterations; ++i)
  {
    double3 sA = double3(random.uniform(), random.uniform(), random.uniform());
    double3 posA = structure.unitCell.cell * sA;

    // A point of the void is covered by the sphere about itself if by no other, so that is where its size starts.
    std::optional<double> ownRadius = structure.freeRadius(posA);
    if (!ownRadius.has_value()) continue;

    double largest_radius = ownRadius.value();
    for (std::size_t j = 0; j < number_of_inner_steps; ++j)
    {
      double3 sB = double3(random.uniform(), random.uniform(), random.uniform());
      double3 posB = structure.unitCell.cell * sB;

      std::optional<double> radius = structure.freeRadius(posB);

      if (!radius.has_value()) continue;

      double3 dr = structure.unitCell.applyPeriodicBoundaryConditions(posA - posB);
      double rr = double3::dot(dr, dr);

      if (rr > radius.value() * radius.value()) continue;

      largest_radius = std::max(largest_radius, radius.value());
    }

    if (!(delta_r > 0.0)) continue;

    // A size past the end of the range is drawn in no bin, and is at or above every one of them.
    const double bin = largest_radius / delta_r;
    const std::size_t index = bin < static_cast<double>(numberOfBins) ? static_cast<std::size_t>(bin) : numberOfBins;
    if (index < numberOfBins) histogram[index]++;
    for (std::size_t k = 0; k <= index && k < numberOfBins; ++k)
    {
      histogram_cummulative[k]++;
    }
  }

  std::chrono::duration<double> timing = std::chrono::steady_clock::now() - time_begin;

  this->seconds = timing.count();

  std::ofstream myfile;
  myfile.open(structure.name + ".mc.psd.cpu.txt");
  std::print(myfile, "# Pore-size distribution using Mont Carlo-based method\n");
  structure.writeHeader(myfile);
  std::print(myfile, "# Number of iterations: {}\n", number_of_iterations);
  std::print(myfile, "# Number of inner-steps (sample points per atom): {}\n", number_of_inner_steps);
  std::print(myfile, "# CPU Timing: {} [s]\n", this->seconds);
  std::print(myfile, "# column 1: diameter d [A]\n");
  std::print(myfile, "# column 2: PSD, the share of the cell per unit of diameter [1/A]\n");
  std::print(myfile, "# column 3: cumulative pore volume\n");
  std::print(myfile, "# value at d=0 is related to the void-fraction\n");

  // The bins are of radius and the curve is against the diameter, so a bin is twice as wide as `delta_r` here.
  double normalization = 1.0 / static_cast<double>(number_of_iterations);
  double width = 2.0 * delta_r;
  for (std::size_t index = 0; index < numberOfBins; ++index)
  {
    std::print(myfile, "{} {} {}\n", width * (static_cast<double>(index) + 0.5),
               width > 0.0 ? histogram[index] * normalization / width : 0.0,
               histogram_cummulative[index] * normalization);
  }

  myfile.close();
}
