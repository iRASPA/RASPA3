module;

module integration_surface_area;

import std;

import double3;
import sampled_structure;

void Integration_SurfaceArea::run(const SampledStructure &structure, const SampledProbe &probe,
                                  std::optional<std::size_t> numberOfSlices)
{
  std::chrono::steady_clock::time_point time_begin = std::chrono::steady_clock::now();

  std::size_t number_of_slices = numberOfSlices.value_or(1024);

  double accumulated_surface_area{};
  for (std::size_t atom_index = 0; atom_index < structure.size(); ++atom_index)
  {
    double equilibrium_distance = structure.radii[atom_index];

    double total_trials{};
    double counted{};
    for (std::size_t stackNumber = 0; stackNumber <= number_of_slices; ++stackNumber)
    {
      for (std::size_t sliceNumber = 0; sliceNumber < number_of_slices; ++sliceNumber)
      {
        double theta =
            static_cast<double>(sliceNumber) * 2.0 * std::numbers::pi / static_cast<double>(number_of_slices);

        double u = static_cast<double>(stackNumber) / static_cast<double>(number_of_slices);
        double phi = std::acos(2.0 * u - 1.0);

        double sinTheta = std::sin(theta);
        double sinPhi = std::sin(phi);
        double cosTheta = std::cos(theta);
        double cosPhi = std::cos(phi);
        double3 unit_vector{sinPhi * cosTheta, sinPhi * sinTheta, cosPhi};

        double3 position = structure.positions[atom_index] + equilibrium_distance * unit_vector;
        if (!structure.overlaps(position, atom_index))
        {
          counted += 1.0;
        }

        total_trials += 1.0;
      }
    }

    accumulated_surface_area +=
        (counted / total_trials) * 4.0 * std::numbers::pi * equilibrium_distance * equilibrium_distance;
  }

  std::chrono::duration<double> timing = std::chrono::steady_clock::now() - time_begin;

  this->seconds = timing.count();
  this->surfaceArea = accumulated_surface_area;

  std::ofstream myfile;
  myfile.open(structure.name + ".integration.sa.cpu.txt");
  std::print(myfile, "# Surface area using integration method\n");
  structure.writeHeader(myfile);
  probe.writeHeader(myfile);
  std::print(myfile, "# Number of integration points per atom: {}\n", (number_of_slices + 1) * number_of_slices);
  std::print(myfile, "# CPU Timing: {} [s]\n", this->seconds);
  std::print(myfile, "{} [Å²]\n", this->surfaceArea);
  std::print(myfile, "{} [m²/cm³]\n", 1.0e4 * this->surfaceArea / structure.unitCell.volume);
  std::print(myfile, "{} [m²/g]\n", this->surfaceArea * structure.gravimetricFactor());
  myfile.close();
}
