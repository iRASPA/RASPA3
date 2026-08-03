module;

module voronoi_pore_diameters;

import std;

import double3;
import crystal;
import pair_interactions;
import voronoi_network;
import pore_diameters;

void VoronoiPoreDiameters::run(const PairInteractions& interactions, const Crystal& framework)
{
  std::chrono::steady_clock::time_point time_begin = std::chrono::steady_clock::now();

  std::vector<double3> fractionalPositions;
  std::vector<double> radii;
  fractionalPositions.reserve(framework.atoms.size());
  radii.reserve(framework.atoms.size());
  for (const CrystalAtom& atom : framework.atoms)
  {
    fractionalPositions.push_back(framework.unitCell.inverseCell * atom.position);
    std::size_t type = atom.type;
    radii.push_back(0.5 * interactions(type, type).sizeParameter);
  }

  network = VoronoiNetwork::create(framework.unitCell, fractionalPositions, radii);
  result = PoreDiameters::compute(network);

  std::chrono::duration<double> timing = std::chrono::steady_clock::now() - time_begin;

  std::ofstream myfile;
  myfile.open(framework.name + ".voronoi.res.txt");
  std::print(myfile, "# Pore diameters from the Voronoi network (Di, Df, Dif)\n");
  std::print(myfile, "# Crystal: {}\n", framework.name);
  std::print(myfile, "# Space-group Hall-number: {}\n", framework.spaceGroupHallNumber);
  std::print(myfile, "# Number of framework atoms: {}\n", framework.atoms.size());
  std::print(myfile, "# Number of Voronoi nodes: {}\n", network.nodes.size());
  std::print(myfile, "# Number of Voronoi edges (directed): {}\n", network.edges.size());
  std::print(myfile, "# CPU Timing: {} [s]\n", timing.count());
  std::print(myfile, "Di (largest included sphere):            {} [Å]\n", result.includedSphereDiameter);
  std::print(myfile, "Df (largest free sphere):                {} [Å]\n", result.freeSphereDiameter);
  std::print(myfile, "Dif (included sphere along free path):   {} [Å]\n", result.includedAlongFreePathDiameter);
  myfile.close();
}
