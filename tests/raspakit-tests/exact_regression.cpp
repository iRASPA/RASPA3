#include <gtest/gtest.h>

import std;

import int3;
import double3;
import simulationbox;
import atom;
import forcefield;
import framework;
import pore_accessibility;
import exact_boundary_components;
import exact_surface_patches;
import exact_solvent_excluded;
import exact_pore_size_distribution;

// The exact analyses run end to end on a real framework, against what they answered before.
//
// Every other test of these analyses checks a property: a volume against a closed form, a derivative against
// the volume it is the derivative of, a predicate against a case worked out by hand. Those are the tests worth
// having and they are elsewhere. What none of them does is notice that a restructuring has moved an answer,
// because each of them holds one piece of the machinery against something outside it and a change that moves
// the whole pipeline together moves the check with it. So this file does the one thing they cannot: it runs the
// pipeline the command line runs, on a framework with the awkward features in it, and holds the answers against
// the answers of the day they were recorded.
//
// ITQ-29 is all-silica LTA and is the smallest structure that has everything these analyses divide by. Its
// large cages form a three-dimensional channel system, and its sodalite cages are sealed to anything larger
// than a six-ring, so the framework has an accessible pore volume, an inaccessible one, closed surfaces and
// running-away ones all at the same probe. The inaccessible volume is also the most delicate number the
// analyses produce, being what the divergence theorem leaves after two larger terms have cancelled.
//
// The tolerance is relative and is 1e-12. It is not tighter than that on purpose. The last two or three digits
// of these numbers are not determined: they move by a part in 1e14 between one thread and eight, because a sum
// reduced over partials is added up in a different order, and they have been seen to move by a part in 1e15
// between two builds of the same source, because the compiler is free to contract a multiply and an add into
// one instruction in one of them and not the other. A test that demanded equality would fail for those reasons
// and teach everyone to ignore it. A part in 1e12 is some three orders of magnitude above that noise and many
// orders below any change of substance, quadrature and root-finding tolerances here being nearer 1e-9.
//
// When a change is meant to move these numbers, run `--gtest_also_run_disabled_tests` with
// `--gtest_filter=exact_regression.DISABLED_record` and paste what it prints over the block below.

namespace
{

struct Recorded
{
  const char* name;
  double expected;
  double measured;
};

// A relative tolerance with a floor, so that a quantity which is legitimately zero is held to an absolute
// tolerance rather than to no tolerance at all.
void expectRecorded(const std::vector<Recorded>& values)
{
  for (const Recorded& value : values)
  {
    const double tolerance = 1.0e-12 * std::max(1.0, std::abs(value.expected));
    EXPECT_NEAR(value.measured, value.expected, tolerance)
        << value.name << " moved by " << (value.measured - value.expected) << ", which is "
        << (value.measured - value.expected) / std::max(1.0e-300, std::abs(value.expected)) << " of it";
  }
}

struct Structure
{
  SimulationBox box;
  std::vector<double3> fractionalPositions;
  std::vector<double> radii;
};

Structure itq29()
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  Framework framework = Framework::makeITQ29(forceField, int3(1, 1, 1));

  Structure structure;
  structure.box = framework.simulationBox;
  for (const Atom& atom : framework.unitCellAtoms)
  {
    structure.fractionalPositions.push_back(framework.simulationBox.inverseCell * atom.position);
    const std::size_t type = static_cast<std::size_t>(atom.type);
    structure.radii.push_back(0.5 * forceField(type, type).sizeParameter());
  }
  return structure;
}

// The probe and the sweep resolution of the recorded run. Both are arguments of the answers, so both are
// written down here rather than left to a default that may move.
constexpr double probeRadius = 1.2;
constexpr std::size_t subdivisions = 1;

// The same call sequence the surface-area and void-volume commands make.
struct Pipeline
{
  BoundaryComponents components;
  MeasuredPatches patches;
  SolventExcludedGeometry geometry;
};

Pipeline runPipeline(const Structure& structure)
{
  const PoreAccessibility accessibility =
      PoreAccessibility::create(structure.box, structure.fractionalPositions, structure.radii, probeRadius);

  Pipeline run;
  run.components = boundaryComponents(accessibility);
  const std::vector<ComponentVerdict> verdicts = boundaryComponentVerdicts(accessibility, run.components);
  run.patches = exactAccessibleSurfaceAreaByComponent(accessibility, run.components, verdicts, subdivisions);
  run.geometry =
      solventExcludedGeometry(accessibility, probeRadius, run.components, verdicts, run.patches, subdivisions);
  return run;
}

// A short sweep: enough diameters to cross both of the framework's pore sizes, and few enough, at few enough
// refinements of the intervals between them, that the whole of it costs a second or two. The refinements are
// what corner the spikes, so the recorded `largestDiameter` is only ever as sharp as this number makes it.
constexpr double maximumDiameter = 12.0;
constexpr std::size_t numberOfBins = 6;
constexpr std::size_t refinements = 4;

PoreSizeDistributionCurve runCurve(const Structure& structure)
{
  auto build = [&](double inflation)
  { return PoreAccessibility::create(structure.box, structure.fractionalPositions, structure.radii, inflation); };

  return exactPoreSizeDistribution(build, structure.box.volume, maximumDiameter, numberOfBins, subdivisions,
                                   probeRadius, refinements);
}

}  // namespace

TEST(exact_regression, itq29_surface_and_volume)
{
  const Pipeline run = runPipeline(itq29());
  const SolventExcludedGeometry& geometry = run.geometry;

  // Counts, which are not floating point and are held exactly. A change in any of them is a change in what the
  // boundary was found to be made of, whatever the volumes then come to. Two connected surfaces: the wall of
  // the channel system and the wall of the sodalite cage. The six degenerate vertices are the framework's own
  // symmetry, four spheres meeting at a point where a less symmetric structure would have three; the clipped
  // ones are concave patches a neighbouring probe position reaches into, which at this probe is most of them.
  EXPECT_EQ(run.components.numberOfComponents, 2uz);
  EXPECT_EQ(run.patches.diagnostics.unplacedArcs, 0uz);
  EXPECT_EQ(geometry.diagnostics.cuspedArcs, 0uz);
  EXPECT_EQ(geometry.diagnostics.clippedVertices, 166uz);
  EXPECT_EQ(geometry.diagnostics.degenerateVertices, 6uz);
  EXPECT_EQ(geometry.diagnostics.vanishedVertices, 0uz);
  EXPECT_EQ(geometry.diagnostics.discardedCorners, 0uz);

  expectRecorded({
      {"accessible area", 265.70743809340621, run.patches.accessible},
      {"inaccessible area", 45.628037154540024, run.patches.inaccessible},
      {"undecided area", 0.0, run.patches.undecided},

      {"accessible volume", 1331.0028094825439, geometry.accessibleVolume},
      {"shell volume", 485.54748713646171, geometry.shellVolume},
      {"excluded volume", 845.45532234608208, geometry.excludedVolume},
      {"pore volume", 825.76537707962939, geometry.poreVolume},

      {"accessible pore volume", 703.88411171993005, geometry.accessiblePoreVolume},
      {"inaccessible pore volume", 121.88126535969934, geometry.inaccessiblePoreVolume},
      {"undecided pore volume", 0.0, geometry.undecidedPoreVolume},

      {"distribution", 59.828140947402701, geometry.distribution},
      {"toroidal distribution", 26.129503970415435, geometry.toroidalDistribution},
      {"concave distribution", 33.69863697698726, geometry.concaveDistribution},

      {"convex area", 104.35344184210972, geometry.area.convex},
      {"saddle area", 244.94747946764224, geometry.area.saddle},
      {"concave area", 142.57755748834552, geometry.area.concave},
  });

  // The three kinds of patch are the whole of the excluded surface and nothing else, so their fractions add to
  // one. This is not a recorded number but a statement about the decomposition, and it holds to round-off.
  EXPECT_NEAR(geometry.area.convexFraction() + geometry.area.saddleFraction() + geometry.area.concaveFraction(),
              1.0, 1.0e-12);

  // The pore volume is divided and not sampled, so the three parts add back to it exactly.
  EXPECT_NEAR(geometry.accessiblePoreVolume + geometry.inaccessiblePoreVolume + geometry.undecidedPoreVolume,
              geometry.poreVolume, 1.0e-9);
}

TEST(exact_regression, itq29_pore_size_distribution)
{
  const PoreSizeDistributionCurve curve = runCurve(itq29());

  ASSERT_EQ(curve.points.size(), numberOfBins);

  expectRecorded({
      {"void volume", 926.80710761972296, curve.voidVolume},
      {"probe accessible volume", 703.88411171993005, curve.probeAccessibleVolume},

      {"integral", 0.26938068862486148, curve.integral},
      {"singular weight", 0.73451983718933156, curve.singularWeight},
      {"largest diameter", 9.536321810724246, curve.largestDiameter},

      {"probe accessible integral", 0.18349738844711069, curve.probeAccessibleIntegral},
      {"probe accessible singular weight", 0.81836703708498715, curve.probeAccessibleSingularWeight},
      {"probe accessible largest diameter", 10.312302665106525, curve.probeAccessibleLargestDiameter},
  });

  // The whole of the void has some pore size, so the continuous part and the spikes come to one between them.
  // How near they come is a property of the grid and not of the analysis: the continuous part is integrated by
  // the trapezium rule over the rows, and six rows over twelve angstroms is a coarse grid deliberately, so four
  // parts in a thousand is what it is worth here. It is the recorded value above that pins this down tightly;
  // this only says that the two parts still account for the void.
  EXPECT_NEAR(curve.integral + curve.singularWeight, 1.0, 5.0e-3);

  // Nothing is left at the end of the range: the largest sphere in the framework is inside it.
  EXPECT_EQ(curve.truncatedWeight, 0.0);

  // Every row of the curve, so that a change anywhere along it is caught and not only at the ends.
  const std::vector<double> cumulative = {0.96809426869770654,  0.8722189426117215,     0.80730914790974595,
                                          0.66204771743467605,  0.63689538105696941,    -2.2079708519334414e-15};
  const std::vector<double> distribution = {0.050566473333284907, 0.031778333346549538, 0.032548376722029336,
                                            0.015318164306166983, 0.010383488910987832, 0.0};

  std::vector<Recorded> rows;
  std::vector<std::string> names;
  names.reserve(2 * curve.points.size());
  for (std::size_t i = 0; i < curve.points.size(); ++i)
  {
    names.push_back(std::format("cumulative at row {}", i));
    rows.push_back({names.back().c_str(), cumulative[i], curve.points[i].cumulative});
    names.push_back(std::format("distribution at row {}", i));
    rows.push_back({names.back().c_str(), distribution[i], curve.points[i].distribution});
  }
  expectRecorded(rows);
}

// Not a test. Prints the block above as it would be written, for the day one of these answers is meant to
// change; see the note at the top of this file.
TEST(exact_regression, DISABLED_record)
{
  const Structure structure = itq29();
  const Pipeline run = runPipeline(structure);
  const SolventExcludedGeometry& geometry = run.geometry;

  auto line = [](const char* name, double value) { std::print("      {{\"{}\", {:.17g}, }},\n", name, value); };

  auto count = [](const char* name, std::size_t value) { std::print("  {} = {}uz\n", name, value); };
  count("numberOfComponents", run.components.numberOfComponents);
  count("unplacedArcs", run.patches.diagnostics.unplacedArcs);
  count("cuspedArcs", geometry.diagnostics.cuspedArcs);
  count("clippedVertices", geometry.diagnostics.clippedVertices);
  count("degenerateVertices", geometry.diagnostics.degenerateVertices);
  count("vanishedVertices", geometry.diagnostics.vanishedVertices);
  count("discardedCorners", geometry.diagnostics.discardedCorners);

  line("accessible area", run.patches.accessible);
  line("inaccessible area", run.patches.inaccessible);
  line("undecided area", run.patches.undecided);
  line("accessible volume", geometry.accessibleVolume);
  line("shell volume", geometry.shellVolume);
  line("excluded volume", geometry.excludedVolume);
  line("pore volume", geometry.poreVolume);
  line("accessible pore volume", geometry.accessiblePoreVolume);
  line("inaccessible pore volume", geometry.inaccessiblePoreVolume);
  line("undecided pore volume", geometry.undecidedPoreVolume);
  line("distribution", geometry.distribution);
  line("toroidal distribution", geometry.toroidalDistribution);
  line("concave distribution", geometry.concaveDistribution);
  line("convex area", geometry.area.convex);
  line("saddle area", geometry.area.saddle);
  line("concave area", geometry.area.concave);

  const PoreSizeDistributionCurve curve = runCurve(structure);
  line("void volume", curve.voidVolume);
  line("probe accessible volume", curve.probeAccessibleVolume);
  line("integral", curve.integral);
  line("singular weight", curve.singularWeight);
  line("largest diameter", curve.largestDiameter);
  line("probe accessible integral", curve.probeAccessibleIntegral);
  line("probe accessible singular weight", curve.probeAccessibleSingularWeight);
  line("probe accessible largest diameter", curve.probeAccessibleLargestDiameter);
  line("truncated weight", curve.truncatedWeight);
  line("integral plus singular weight", curve.integral + curve.singularWeight);

  std::print("  const std::vector<double> cumulative = {{");
  for (const PoreSizeDistributionPoint& point : curve.points) std::print("{:.17g}, ", point.cumulative);
  std::print("}};\n  const std::vector<double> distribution = {{");
  for (const PoreSizeDistributionPoint& point : curve.points) std::print("{:.17g}, ", point.distribution);
  std::print("}};\n");
}
