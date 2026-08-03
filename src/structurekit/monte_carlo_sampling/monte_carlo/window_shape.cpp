module;

module mc_window_shape;

import std;

import pore_diameters;
import pore_window;
import sampled_structure;
import sampled_roadmap;
import mc_backend;

void MC_WindowShape::run(const SampledStructure &structure, std::optional<std::size_t> numberOfIterations,
                         std::optional<std::size_t> numberOfInnerSteps)
{
  run(structure, SampledRoadmap::build(structure, samplingBackendCPU(), numberOfIterations, numberOfInnerSteps));
}

void MC_WindowShape::run(const SampledStructure &structure, const SampledRoadmap &roadmap)
{
  std::chrono::steady_clock::time_point time_begin = std::chrono::steady_clock::now();

  this->roadmap = roadmap;

  PoreDiameters diameters = PoreDiameters::compute(roadmap.network);
  this->percolates = diameters.freeSphereDiameter > 0.0;
  this->freeSphereDiameter = diameters.freeSphereDiameter;

  this->freeSphere = freeSphereWindow(roadmap.network);
  this->channels = channelWindows(roadmap.network, roadmap.components);

  std::chrono::duration<double> timing = std::chrono::steady_clock::now() - time_begin;
  this->seconds = timing.count();

  std::ofstream myfile;
  myfile.open(std::format("{}.mc.window.{}.txt", structure.name, roadmap.backend));
  std::print(myfile, "# The shape of the pore windows, by sampling the void\n");
  structure.writeHeader(myfile);
  roadmap.writeHeader(myfile);
  std::print(myfile, "# Timing (windows, on the processor either way): {} [s]\n", this->seconds);
  std::print(myfile, "# Each window is the plane across the passage where the bottleneck sits. Df is the\n");
  std::print(myfile, "# inscribed circle of it and says nothing about how far it reaches the other ways,\n");
  std::print(myfile, "# which these two pairs of numbers do: the free width is the shortest and the longest\n");
  std::print(myfile, "# free chord through the bottleneck, and the ellipse is the largest-area one that fits\n");
  std::print(myfile, "# there, which may be narrower than Df one way in return for being longer the other.\n");
  std::print(myfile, "# Both are lower bounds, each direction being followed only as far as the first atom.\n");
  std::print(myfile, "#\n");
  std::print(myfile, "# Where the plane is depends on the sample, the bottleneck being the narrowest point of\n");
  std::print(myfile, "# the narrowest hop the sample found; what is measured in it does not. Neither is a\n");
  std::print(myfile, "# criterion for a non-spherical molecule, whose passage turns on its own shape and on\n");
  std::print(myfile, "# how it may turn along the way. They describe the window.\n");

  auto writeWindow = [&](const PoreWindow &measured)
  {
    std::print(myfile, " free width {:.5f} - {:.5f} [Å] ellipse {:.5f} x {:.5f} [Å] ring {} atoms{}\n",
               measured.smallestFreeWidth, measured.largestFreeWidth, measured.minorAxis, measured.majorAxis,
               measured.boundingAtoms,
               measured.clipped ? " (open in some direction, cut off at the reach of the cell)" : "");
    std::print(myfile,
               "  window at ({:.5f}, {:.5f}, {:.5f}) [Å], across ({:.5f}, {:.5f}, {:.5f}), major axis along "
               "({:.5f}, {:.5f}, {:.5f})\n",
               measured.position.x, measured.position.y, measured.position.z, measured.normal.x,
               measured.normal.y, measured.normal.z, measured.majorAxisDirection.x, measured.majorAxisDirection.y,
               measured.majorAxisDirection.z);
  };

  std::print(myfile, "Df ({:.5f} [Å]) window:", this->freeSphereDiameter);
  if (this->freeSphere.measured)
    writeWindow(this->freeSphere);
  else if (!this->percolates)
    std::print(myfile, " none, no path across the crystal was found\n");
  else
    std::print(myfile, " the bottleneck is shut, leaving no window to measure\n");

  std::print(myfile, "\n# Per channel, at the bottleneck of that channel's own widest percolating path.\n");
  std::print(myfile, "# Channels are numbered as in the .mc.chan file.\n");
  for (const ChannelWindow &channel : this->channels)
  {
    std::print(myfile, "channel {} (dimensionality {}): Df {:.5f} [Å]\n", channel.poreIndex,
               channel.dimensionality, channel.freeSphereDiameter);
    if (channel.window.measured)
      writeWindow(channel.window);
    else
      std::print(myfile, " window not measured (bottleneck clearance {:.5f} [Å])\n", channel.window.freeRadius);
  }
  myfile.close();
}
