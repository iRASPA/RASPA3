module;

module integration_surface_area;

import std;

import double3;
import sampled_structure;
import structure_parallel;

void Integration_SurfaceArea::run(const SampledStructure &structure, const SampledProbe &probe,
                                  std::optional<std::size_t> numberOfSlices)
{
  std::chrono::steady_clock::time_point time_begin = std::chrono::steady_clock::now();

  std::size_t number_of_slices = numberOfSlices.value_or(1024);
  // the estimate divides by it (no slices gave 0 / 0, a not-a-number reported as the area)
  if (number_of_slices == 0)
  {
    throw std::invalid_argument("Integration_SurfaceArea: the number of slices must be at least one\n");
  }

  // The directions are the same for every atom: their sines and cosines are tabulated once per stack and per
  // slice (as they were computed, so the directions are the same numbers), instead of four trigonometric calls
  // and an arc cosine per direction per atom.
  // The stacks are equal shares of the height along the axis, which are equal shares of the area, and each is
  // taken at its middle. Taken at its edges there is one stack too many, and a whole stack at each pole.
  std::vector<double> sinPhi(number_of_slices), cosPhi(number_of_slices), sinTheta(number_of_slices),
      cosTheta(number_of_slices);
  for (std::size_t stackNumber = 0; stackNumber < number_of_slices; ++stackNumber)
  {
    double u = (static_cast<double>(stackNumber) + 0.5) / static_cast<double>(number_of_slices);
    double phi = std::acos(2.0 * u - 1.0);
    sinPhi[stackNumber] = std::sin(phi);
    cosPhi[stackNumber] = std::cos(phi);
  }
  for (std::size_t sliceNumber = 0; sliceNumber < number_of_slices; ++sliceNumber)
  {
    double theta = static_cast<double>(sliceNumber) * 2.0 * std::numbers::pi / static_cast<double>(number_of_slices);
    sinTheta[sliceNumber] = std::sin(theta);
    cosTheta[sliceNumber] = std::cos(theta);
  }

  // The atoms are independent: dealt out over the pool in blocks, each atom's share kept in its own slot, and the
  // total taken in atom order afterwards, so the result is the one of the single loop.
  std::vector<double> areaOfAtom(structure.size(), 0.0);
  forEachBlock(structure.size(), workersAvailable(),
               [&](std::size_t, std::size_t begin, std::size_t end)
               {
                 for (std::size_t atom_index = begin; atom_index < end; ++atom_index)
                 {
                   double equilibrium_distance = structure.radii[atom_index];

                   double total_trials{};
                   double counted{};
                   for (std::size_t stackNumber = 0; stackNumber < number_of_slices; ++stackNumber)
                   {
                     for (std::size_t sliceNumber = 0; sliceNumber < number_of_slices; ++sliceNumber)
                     {
                       double3 unit_vector{sinPhi[stackNumber] * cosTheta[sliceNumber],
                                           sinPhi[stackNumber] * sinTheta[sliceNumber], cosPhi[stackNumber]};

                       double3 position = structure.positions[atom_index] + equilibrium_distance * unit_vector;
                       if (!structure.overlaps(position, atom_index))
                       {
                         counted += 1.0;
                       }

                       total_trials += 1.0;
                     }
                   }

                   areaOfAtom[atom_index] =
                       (counted / total_trials) * 4.0 * std::numbers::pi * equilibrium_distance * equilibrium_distance;
                 }
               });

  double accumulated_surface_area{};
  for (double area : areaOfAtom) accumulated_surface_area += area;

  std::chrono::duration<double> timing = std::chrono::steady_clock::now() - time_begin;

  this->seconds = timing.count();
  this->surfaceArea = accumulated_surface_area;

  std::ofstream myfile;
  myfile.open(structure.name + ".integration.sa.cpu.txt");
  std::print(myfile, "# Surface area using integration method\n");
  structure.writeHeader(myfile);
  probe.writeHeader(myfile);
  std::print(myfile, "# Number of integration points per atom: {}\n", number_of_slices * number_of_slices);
  std::print(myfile, "# CPU Timing: {} [s]\n", this->seconds);
  std::print(myfile, "{} [Å²]\n", this->surfaceArea);
  std::print(myfile, "{} [m²/cm³]\n", 1.0e4 * this->surfaceArea / structure.unitCell.volume);
  std::print(myfile, "{} [m²/g]\n", this->surfaceArea * structure.gravimetricFactor());
  myfile.close();
}
