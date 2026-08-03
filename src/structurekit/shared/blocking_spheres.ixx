module;

export module blocking_spheres;

import std;

import double3;

// A sphere covering part of a pocket a probe cannot reach, for a simulation to reject insertions in.
//
// A blocking sphere has two jobs and they pull against each other: between them the spheres of a pocket
// have to hold all of it, or a molecule still gets into the part they miss, and none of them may reach a
// channel, or the simulation loses part of its own pore. Every route here answers to both; what they differ
// in is what they measure the two against.
export struct BlockingSphere
{
  double3 centerFractional;
  double radius;  // Å
};

// The `.block` file a simulation reads back, which is a count and one line of fractional centre and radius
// per sphere and nothing else, no comment being allowed in it. What the spheres are and where they came
// from goes in a report of the route's own beside it.
export void writeBlockFile(const std::string& frameworkName, const std::vector<BlockingSphere>& spheres);
