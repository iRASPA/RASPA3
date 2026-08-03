module;

module mc_pore_diameters;

import std;

import pore_diameters;
import sampled_structure;
import sampled_roadmap;
import mc_backend;

void MC_PoreDiameters::run(const SampledStructure &structure, std::optional<std::size_t> numberOfIterations,
                           std::optional<std::size_t> numberOfInnerSteps)
{
  run(structure, SampledRoadmap::build(structure, samplingBackendCPU(), numberOfIterations, numberOfInnerSteps));
}

void MC_PoreDiameters::run(const SampledStructure &structure, const SampledRoadmap &roadmap)
{
  std::chrono::steady_clock::time_point time_begin = std::chrono::steady_clock::now();

  this->roadmap = roadmap;
  this->result = PoreDiameters::compute(roadmap.network);
  this->percolates = this->result.freeSphereDiameter > 0.0;

  std::chrono::duration<double> timing = std::chrono::steady_clock::now() - time_begin;
  this->seconds = timing.count();

  std::ofstream myfile;
  myfile.open(std::format("{}.mc.res.{}.txt", structure.name, roadmap.backend));
  std::print(myfile, "# Pore diameters (Di, Df, Dif) by sampling the void\n");
  structure.writeHeader(myfile);
  roadmap.writeHeader(myfile);
  std::print(myfile, "# Timing (diameters, on the processor either way): {} [s]\n", this->seconds);
  std::print(myfile, "# Every diameter here is a lower bound: it is what a sphere was shown to do, and a\n");
  std::print(myfile, "# larger sample can only find more. Run two sizes to see how much is left.\n");
  std::print(myfile, "Di (largest included sphere):            {} [Å]\n", this->result.includedSphereDiameter);
  if (this->percolates)
  {
    std::print(myfile, "Df (largest free sphere):                {} [Å]\n", this->result.freeSphereDiameter);
    std::print(myfile, "Dif (included sphere along free path):   {} [Å]\n",
               this->result.includedAlongFreePathDiameter);
  }
  else
  {
    std::print(myfile, "Df (largest free sphere):                no path across the crystal was found\n");
    std::print(myfile, "Dif (included sphere along free path):   no path across the crystal was found\n");
  }
  myfile.close();
}
