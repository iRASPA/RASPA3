module;

module energy_shared_energy_barrier;

import std;

import int3;
import uint3;
import double3;
import double3x3;
import skspacegroupdatabase;
import crystal;
import pair_interactions;
import units;
import energy_shared_probe_energy_grid;
import energy_shared_energy_backend;
import grid_percolation;


namespace
{
// An opt-in check of the field a backend built, against the same sum taken on the processor in double
// precision and over a wider range of images than the builder allows itself. It is the builder's shortcuts
// that are on trial here: the float arithmetic, the image range the perpendicular widths settle, and the
// test that throws out an image before building it.
void writeFieldCheck(const PairInteractions &interactions, const Crystal &framework, const ProbeEnergyGrid &energyGrid)
{
  std::optional<std::size_t> probeType = interactions.findType(energyGrid.probeName);
  std::vector<double3> fractionalPositions = framework.fractionalPositions;
  double3x3 cell = framework.unitCell.cell;
  const double cutOffSquared = energyGrid.cutOff * energyGrid.cutOff;
  const std::int32_t extra = 1;

  std::ofstream dump;
  dump.open(framework.name + ".energyfield.check.txt");
  std::print(dump, "# i j k  kernel[internal]  reference[internal]  relative difference\n");
  for (std::size_t k = 0; k < energyGrid.gridSize.z; k += 13)
  {
    for (std::size_t j = 0; j < energyGrid.gridSize.y; j += 11)
    {
      for (std::size_t i = 0; i < energyGrid.gridSize.x; i += 7)
      {
        double3 s = energyGrid.fractionalPosition(i, j, k);

        double reference = 0.0;
        for (std::size_t atom = 0; atom < fractionalPositions.size(); ++atom)
        {
          std::size_t atomType = framework.atoms[atom].type;
          double sigma = interactions(probeType.value(), atomType).sizeParameter;
          double epsilon = interactions(probeType.value(), atomType).strengthParameter;

          double3 ds = s - fractionalPositions[atom];
          ds.x -= std::rint(ds.x);
          ds.y -= std::rint(ds.y);
          ds.z -= std::rint(ds.z);

          for (std::int32_t a = -energyGrid.numberOfImageShells.x - extra;
               a <= energyGrid.numberOfImageShells.x + extra; ++a)
          {
            for (std::int32_t b = -energyGrid.numberOfImageShells.y - extra;
                 b <= energyGrid.numberOfImageShells.y + extra; ++b)
            {
              for (std::int32_t c = -energyGrid.numberOfImageShells.z - extra;
                   c <= energyGrid.numberOfImageShells.z + extra; ++c)
              {
                double3 t = ds + double3(static_cast<double>(a), static_cast<double>(b), static_cast<double>(c));
                double3 dr = cell * t;
                double rr = double3::dot(dr, dr);
                if (rr >= cutOffSquared) continue;
                rr = std::max(rr, 1.0e-6);

                double ratio3 = (sigma * sigma / rr) * (sigma * sigma / rr) * (sigma * sigma / rr);
                reference += std::min(4.0 * epsilon * ratio3 * (ratio3 - 1.0), energyGrid.ceiling);
              }
            }
          }
        }
        reference = std::min(reference, energyGrid.ceiling);

        double kernel = static_cast<double>(energyGrid.energy[energyGrid.voxelIndex(i, j, k)]);
        double scale = std::max({std::fabs(kernel), std::fabs(reference), 1.0});
        std::print(dump, "{} {} {} {:.10e} {:.10e} {:.3e}\n", i, j, k, kernel, reference,
                   std::fabs(kernel - reference) / scale);
      }
    }
  }
  dump.close();
}
}  // namespace


EnergyBarrier::EnergyBarrier() {}


EnergyBarrier::~EnergyBarrier() {}

std::pair<int, double> EnergyBarrier::dimensionalityWithin(double margin) const
{
  if (!this->percolates) return {0, 0.0};

  int reached = this->dimensionalityAtBarrier;
  double cost = 0.0;
  for (int dimension = this->dimensionalityAtBarrier + 1; dimension <= 3; ++dimension)
  {
    const double next = this->barrierByDimension[static_cast<std::size_t>(dimension - 1)];
    if (next >= std::numeric_limits<double>::max()) break;

    const double above = next - this->percolationBarrier;
    if (above > margin) break;

    reached = dimension;
    cost = above;
  }
  return {reached, cost};
}


void EnergyBarrier::run(const EnergyBackend &backend, const PairInteractions &interactions, const Crystal &framework,
                        std::string probePseudoAtom, uint3 gridSize, double atTemperature)
{
  ProbeEnergyGrid energyGrid = backend.probeEnergyGrid(interactions, framework, probePseudoAtom, gridSize);
  this->run(interactions, framework, energyGrid, atTemperature);
}


EnergyBarrier EnergyBarrier::fromField(uint3 gridSize, std::span<const float> energy, double atTemperature)
{
  EnergyBarrier barrier;
  barrier.temperature = atTemperature;
  if (energy.empty()) return barrier;

  // The sweep reads its field as an openness, larger meaning more open, and an energy is the other way
  // about: the probe passes where the energy is low. Negating is the whole of the translation between them.
  std::vector<float> openness(energy.size());
  std::ranges::transform(energy, openness.begin(), [](float value) { return -value; });

  // Every point takes part, including those buried in the framework. There is no threshold to set and none
  // to get wrong: a path through an atom is a path whose highest point is enormous, and the sweep passes it
  // over in favour of the window long before it reaches that height. Letting every point in is what makes
  // the answer the global saddle rather than the best route within some region chosen beforehand.
  GridPercolation swept = sweepPercolation(gridSize, openness, std::numeric_limits<float>::lowest());

  barrier.numberOfVoxels = swept.numberOfVoxels;
  barrier.sweepSeconds = swept.seconds;
  barrier.percolates = swept.percolates;
  barrier.dimensionalityAtBarrier = swept.dimensionalityAtThreshold;
  barrier.deepestWell = -swept.highestOpenness;

  if (swept.percolates)
  {
    barrier.percolationBarrier = -swept.percolationOpenness;
    barrier.deepestWellOnPath = -swept.highestOpennessOnPath;
  }

  for (std::size_t dimension = 0; dimension < 3; ++dimension)
  {
    double value = swept.opennessByDimension[dimension];
    barrier.barrierByDimension[dimension] =
        (value == std::numeric_limits<double>::lowest()) ? std::numeric_limits<double>::max() : -value;
  }

  return barrier;
}


void EnergyBarrier::run(const PairInteractions &interactions, const Crystal &framework,
                        const ProbeEnergyGrid &energyGrid, double atTemperature)
{
  if (energyGrid.numberOfVoxels() == 0)
  {
    this->grid = energyGrid;
    this->temperature = atTemperature;
    return;
  }

  if (std::getenv("RASPA_VERIFY_ENERGY_FIELD") != nullptr) writeFieldCheck(interactions, framework, energyGrid);

  *this = EnergyBarrier::fromField(energyGrid.gridSize, energyGrid.energy, atTemperature);
  this->grid = energyGrid;

  const double never = std::numeric_limits<double>::max();

  const double toKelvin = Units::EnergyToKelvin;
  const double kT = Units::KB * atTemperature;
  double3 spacing = energyGrid.spacing();

  auto asKelvin = [&](double energy) { return energy * toKelvin; };
  auto asKJPerMol = [&](double energy) { return energy * Units::EnergyToKJPerMol; };

  std::ofstream myfile;
  myfile.open(framework.name + ".grid.barrier." + energyGrid.backend + ".txt");
  std::print(myfile, "# Percolation barrier from the probe energy grid\n");
  std::print(myfile, "# Crystal: {}\n", framework.name);
  std::print(myfile, "# Space-group Hall-number: {}\n", framework.spaceGroupHallNumber);
  std::print(myfile, "# Space-group HM-symbol: {}\n",
             SKSpaceGroupDataBase::spaceGroupData[framework.spaceGroupHallNumber].HMString());
  std::print(myfile, "# Number of framework atoms: {}\n", framework.atoms.size());
  std::print(myfile, "# Crystal volume: {} [Å³]\n", framework.unitCell.volume);
  std::print(myfile, "# Probe atom: {} strength-value/kʙ: {} [K] size-value: {} [Å]\n", energyGrid.probeName,
             energyGrid.probeEpsilon * toKelvin, energyGrid.probeSigma);
  std::print(myfile, "# Cutoff: {} [Å], periodic images searched: {} x {} x {} either way\n", energyGrid.cutOff,
             energyGrid.numberOfImageShells.x, energyGrid.numberOfImageShells.y, energyGrid.numberOfImageShells.z);
  std::print(myfile, "# Grid: {} x {} x {} points, spacing {:.5f} x {:.5f} x {:.5f} [Å]\n", energyGrid.gridSize.x,
             energyGrid.gridSize.y, energyGrid.gridSize.z, spacing.x, spacing.y, spacing.z);
  std::print(myfile, "# Grid points swept: {}\n", this->numberOfVoxels);
  std::print(myfile, "# Temperature the barrier is read against: {} [K]\n", atTemperature);
  std::print(myfile, "# Timing ({}): {} [s] for the energy field\n", energyGrid.backend, energyGrid.seconds);
  std::print(myfile, "# CPU Timing: {} [s] for the sweep\n", this->sweepSeconds);
  std::print(myfile, "#\n");
  std::print(myfile, "# The barrier is the lowest energy at which a connected path runs from one cell to the\n");
  std::print(myfile, "# next: the probe cannot get through the crystal without being raised to it somewhere, and\n");
  std::print(myfile, "# it can get through if it is. It is the saddle of the landscape, and it is the global one,\n");
  std::print(myfile, "# since the sweep passes every point of the field and takes the least height at which the\n");
  std::print(myfile, "# passage first exists rather than the height of some path chosen in advance.\n");
  std::print(myfile, "#\n");
  std::print(myfile, "# It is read off a grid, so it is reached from above: the true saddle falls between grid\n");
  std::print(myfile, "# points and a coarse grid samples the window off its centre, which can only overstate the\n");
  std::print(myfile, "# height. The energy is also steep near a wall, so a barrier through a narrow window moves\n");
  std::print(myfile, "# more with the grid than a barrier through a wide one, and is worth refining twice to see.\n");
  std::print(myfile, "\n");

  if (this->percolates)
  {
    std::print(myfile, "Percolation barrier:            {:14.4f} [K] {:12.4f} [kJ/mol] {:10.2f} kT\n",
               asKelvin(this->percolationBarrier), asKJPerMol(this->percolationBarrier),
               this->percolationBarrier / kT);
    std::print(myfile, "Dimensionality at the barrier:  {}\n", this->dimensionalityAtBarrier);

    // Read off at exactly the barrier, a framework whose channels are alike in more than one direction shows
    // one of them, the rest opening a hair above. Say what is open within a kT beside it.
    const auto [within, cost] = this->dimensionalityWithin(kT);
    if (within > this->dimensionalityAtBarrier)
    {
      std::print(myfile, "Dimensionality within a kT:     {}, the last of them {:.4f} [K] above\n", within,
                 asKelvin(cost));
    }
  }
  else
  {
    std::print(myfile, "No path runs through the crystal at any energy, which cannot happen while every\n");
    std::print(myfile, "grid point takes part and so points at a fault rather than at the structure.\n");
  }

  std::print(myfile, "Deepest site anywhere:          {:14.4f} [K] {:12.4f} [kJ/mol] {:10.2f} kT\n",
             asKelvin(this->deepestWell), asKJPerMol(this->deepestWell), this->deepestWell / kT);

  if (this->percolates)
  {
    std::print(myfile, "Deepest site on the path:       {:14.4f} [K] {:12.4f} [kJ/mol] {:10.2f} kT\n",
               asKelvin(this->deepestWellOnPath), asKJPerMol(this->deepestWellOnPath),
               this->deepestWellOnPath / kT);

    // What a molecule sitting in the deepest site of the percolating network has to find in order to leave
    // it, which is the quantity a hopping rate is written in terms of.
    double activation = this->percolationBarrier - this->deepestWellOnPath;
    std::print(myfile, "\n# The height a probe resting in the deepest site of that network has to climb to leave\n");
    std::print(myfile, "# it, which is the activation energy a hopping rate is written in terms of. The factor\n");
    std::print(myfile, "# beside it is exp(-activation/kT), the part of an Arrhenius rate this route can supply;\n");
    std::print(myfile, "# the attempt frequency in front of it is a dynamical quantity and is not read off a\n");
    std::print(myfile, "# static grid.\n");
    std::print(myfile, "Activation energy to leave:     {:14.4f} [K] {:12.4f} [kJ/mol] {:10.2f} kT\n",
               asKelvin(activation), asKJPerMol(activation), activation / kT);
    std::print(myfile, "Boltzmann factor exp(-Ea/kT):   {:14.6e}\n", std::exp(-activation / kT));
  }

  std::print(myfile, "\n# The barrier to running in at least one, two, and three directions. The sweep passes\n");
  std::print(myfile, "# each of these on its way up, so they come from the same pass as the barrier itself. A\n");
  std::print(myfile, "# structure whose channels run one way has a low barrier along them and a high one across.\n");
  for (std::size_t dimension = 0; dimension < 3; ++dimension)
  {
    double barrier = this->barrierByDimension[dimension];
    if (barrier == never)
    {
      std::print(myfile, "Barrier in {} direction(s) or more: never\n", dimension + 1);
    }
    else
    {
      std::print(myfile, "Barrier in {} direction(s) or more: {:14.4f} [K] {:12.4f} [kJ/mol] {:10.2f} kT\n",
                 dimension + 1, asKelvin(barrier), asKJPerMol(barrier), barrier / kT);
    }
  }

  myfile.close();
}
