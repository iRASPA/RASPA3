#include <gtest/gtest.h>

import std;

import int3;
import double3;
import simulationbox;
import randomnumbers;
import forcefield;
import skapolloniusdiagram;
import apollonius_network;
import apollonius_accessibility;
import voronoi_network;
import voronoi_pore_diameters;
import voronoi_channels;
import pore_accessibility;
import voronoi_accessible_volume;
import voronoi_surface_area;
import pore_window;

// The pore properties zeo++ reports, read off the Apollonius diagram instead of the radical Voronoi
// network. The two networks are the same kind of object, so what is checked here is that the
// translation from the diagram keeps every node and arc, and that the properties that come out are
// the ones the diagram says they should be.

// Simple-cubic lattice of one atom per cell. With all radii equal the Apollonius diagram is the
// ordinary Voronoi diagram, so the analytic answers for the cube apply:
//   Di = sqrt(3) a - 2r (the corner of the cube), Df = sqrt(2) a - 2r (the middle of a cube edge).
// The corner is where eight cells meet, which makes this the most degenerate input there is: one
// vertex touching eight sites at once.
TEST(apollonius_pore_analysis, simple_cubic_pore_diameters)
{
  double a = 5.0;
  double r = 1.0;
  SimulationBox box(a, a, a);

  ApolloniusPoreNetwork pores = ApolloniusPoreNetwork::create(box, {double3(0.0, 0.0, 0.0)}, {r});

  EXPECT_TRUE(pores.networkIsComplete());
  EXPECT_EQ(pores.network.nodes.size(), 1);
  EXPECT_NEAR(pores.network.nodes[0].radius, 0.5 * std::sqrt(3.0) * a - r, 1.0e-6);

  PoreDiameters diameters = PoreDiameters::compute(pores.network);
  EXPECT_NEAR(diameters.includedSphereDiameter, std::sqrt(3.0) * a - 2.0 * r, 1.0e-6);
  EXPECT_NEAR(diameters.freeSphereDiameter, std::sqrt(2.0) * a - 2.0 * r, 1.0e-6);
  EXPECT_NEAR(diameters.includedAlongFreePathDiameter, std::sqrt(3.0) * a - 2.0 * r, 1.0e-6);

  ChannelAnalysis channels = ChannelAnalysis::compute(pores.network, 0.5);
  EXPECT_EQ(channels.numberOfChannels, 1);
  EXPECT_EQ(channels.numberOfPockets, 0);
  ASSERT_EQ(channels.pores.size(), 1);
  EXPECT_EQ(channels.pores[0].dimensionality, 3);
}

// Every vertex of the diagram becomes one node and every arc becomes one edge each way, carrying the
// diagram's own clearance and bottleneck. Nothing is recomputed on the way across, so a discrepancy
// here is a translation that lost or altered part of the diagram.
TEST(apollonius_pore_analysis, network_holds_the_whole_diagram)
{
  double a = 12.0;
  SimulationBox box(a, a, a);

  RandomNumber random{std::optional<std::size_t>(42)};
  std::vector<double3> fractionalPositions;
  std::vector<double> radii;
  for (std::size_t i = 0; i < 25; ++i)
  {
    fractionalPositions.push_back(double3(random.uniform(), random.uniform(), random.uniform()));
    radii.push_back(0.8 + 0.9 * random.uniform());
  }

  SKApolloniusDiagram diagram =
      SKApolloniusDiagram::create(box.cell, fractionalPositions, radii, 1, SKApolloniusRegion::FreeSpace);
  ApolloniusPoreNetwork pores = ApolloniusPoreNetwork::create(box, fractionalPositions, radii);

  ASSERT_EQ(pores.network.nodes.size(), diagram.vertices.size());
  for (std::size_t v = 0; v < diagram.vertices.size(); ++v)
  {
    EXPECT_DOUBLE_EQ(pores.network.nodes[v].radius, diagram.vertices[v].radius);
    // A vertex is the maximum of the clearance in its own right, so the node needs no ascent from it.
    EXPECT_DOUBLE_EQ(pores.network.nodes[v].maximalRadius, diagram.vertices[v].radius);
  }

  std::size_t numberOfRings = 0;
  for (const SKApolloniusEdge& edge : diagram.edges)
    if (edge.isLoop) ++numberOfRings;
  EXPECT_EQ(pores.numberOfRings, numberOfRings);
  EXPECT_EQ(pores.numberOfArcs, diagram.edges.size() - numberOfRings);
  ASSERT_EQ(pores.network.edges.size(), 2 * pores.numberOfArcs);

  // Each arc appears both ways round, with the lattice shift reversed along with it.
  std::multiset<std::tuple<std::size_t, std::size_t, int, int, int>> directed;
  for (const VoronoiEdge& edge : pores.network.edges)
  {
    directed.insert({edge.from, edge.to, edge.delta.x, edge.delta.y, edge.delta.z});
  }
  for (const VoronoiEdge& edge : pores.network.edges)
  {
    EXPECT_EQ(directed.count({edge.to, edge.from, -edge.delta.x, -edge.delta.y, -edge.delta.z}),
              directed.count({edge.from, edge.to, edge.delta.x, edge.delta.y, edge.delta.z}));
  }

  // Every atom that a vertex touches is reachable from that atom, at the image of the vertex the
  // tangency belongs to, which is what the accessibility test walks.
  for (std::size_t v = 0; v < diagram.vertices.size(); ++v)
  {
    const SKApolloniusVertex& vertex = diagram.vertices[v];
    for (std::size_t s = 0; s < vertex.siteIndices.size(); ++s)
    {
      std::size_t atomIndex = vertex.siteIndices[s];
      bool found = false;
      for (const auto& [nodeIndex, relative] : pores.network.atomNodeVectors[atomIndex])
      {
        if (nodeIndex != v) continue;
        // The tangency puts the vertex at distance radius + clearance from the atom it touches.
        if (std::abs(relative.length() - (radii[atomIndex] + vertex.radius)) < 1.0e-8) found = true;
      }
      EXPECT_TRUE(found) << "vertex " << v << " is not reachable from site " << atomIndex;
    }
  }
}

// What the diagram buys over the radical network. The largest sphere that fits anywhere is the
// diagram's deepest vertex, exactly; and because the widest path through the structure runs along the
// arcs of the Apollonius diagram, no path found on the radical network can be wider. Both diameters
// must therefore come out at least as large as the radical network's, and Di must agree with what the
// diagram itself reports.
TEST(apollonius_pore_analysis, diameters_are_at_least_the_radical_ones)
{
  double a = 12.0;
  SimulationBox box(a, a, a);

  RandomNumber random{std::optional<std::size_t>(7)};
  std::vector<double3> fractionalPositions;
  std::vector<double> radii;
  for (std::size_t i = 0; i < 25; ++i)
  {
    fractionalPositions.push_back(double3(random.uniform(), random.uniform(), random.uniform()));
    radii.push_back(0.8 + 0.9 * random.uniform());
  }

  ApolloniusPoreNetwork pores = ApolloniusPoreNetwork::create(box, fractionalPositions, radii);
  ASSERT_TRUE(pores.networkIsComplete());
  PoreDiameters apollonius = PoreDiameters::compute(pores.network);

  SKApolloniusDiagram diagram =
      SKApolloniusDiagram::create(box.cell, fractionalPositions, radii, 1, SKApolloniusRegion::FreeSpace);
  EXPECT_NEAR(apollonius.includedSphereDiameter, 2.0 * diagram.largestEmptySphereRadius(), 1.0e-9);

  VoronoiNetwork radical = VoronoiNetwork::create(box, fractionalPositions, radii);
  PoreDiameters voronoi = PoreDiameters::compute(radical);

  EXPECT_GE(apollonius.includedSphereDiameter, voronoi.includedSphereDiameter - 1.0e-6);
  EXPECT_GE(apollonius.freeSphereDiameter, voronoi.freeSphereDiameter - 1.0e-6);

  // Orderings that hold of any pore network.
  EXPECT_GE(apollonius.includedSphereDiameter, apollonius.freeSphereDiameter);
  EXPECT_GE(apollonius.includedAlongFreePathDiameter, apollonius.freeSphereDiameter);
}

// The clearance of a point: how far a probe centred there is from the nearest surface, negative if it
// is inside an atom. The radii handed in are already inflated by the probe.
static double clearanceAt(const SimulationBox& box, const std::vector<double3>& fractionalPositions,
                          const std::vector<double>& radii, const double3& point)
{
  double best = std::numeric_limits<double>::max();
  for (std::size_t i = 0; i < fractionalPositions.size(); ++i)
  {
    double3 delta = box.applyPeriodicBoundaryConditions(box.cell * fractionalPositions[i] - point);
    best = std::min(best, delta.length() - radii[i]);
  }
  return best;
}

// The surface-area and volume samplers ask the classifier one question per point, so the classifier is
// what has to be right. The cells of the Apollonius diagram are cut by the clearance, so the atom a
// point is attributed to is the one really nearest to it, and whether the point is solid is decided by
// that same atom. Both are checked here against a direct search over every atom: no point may be left
// undecidable, and the solid/void split must be exact.
//
// The configuration is the one that empties cells of the radical diagram: a small atom squeezed
// between two large ones, close enough that the radical planes against both neighbours pass it on the
// wrong side. The Apollonius diagram of the same atoms gives the small atom the cell it deserves.
TEST(apollonius_pore_analysis, classifier_is_exact_about_solid_with_heterogeneous_radii)
{
  double a = 12.0;
  SimulationBox box(a, a, a);
  const double largeRadius = 2.2;
  const double smallRadius = 0.35;
  const double probe = 0.4;
  const double separation = 2.0;

  std::vector<double3> fractionalPositions;
  std::vector<double> radii;
  for (double offset : {0.25, 0.75})
  {
    double3 centre(0.5, offset, offset);
    fractionalPositions.push_back(centre - double3(separation / a, 0.0, 0.0));
    radii.push_back(largeRadius);
    fractionalPositions.push_back(centre + double3(separation / a, 0.0, 0.0));
    radii.push_back(largeRadius);
    fractionalPositions.push_back(centre);
    radii.push_back(smallRadius);
  }

  ASSERT_LT(separation, std::sqrt(std::pow(largeRadius + probe, 2.0) - std::pow(smallRadius + probe, 2.0)));

  std::vector<double> inflatedRadii;
  for (double r : radii) inflatedRadii.push_back(r + probe);

  ApolloniusAccessibility classifier = ApolloniusAccessibility::create(box, fractionalPositions, radii, probe);
  ASSERT_TRUE(classifier.diagram.networkIsComplete());

  // A grid rather than random points, so that a failure is reproducible and the whole cell is covered.
  const int steps = 30;
  std::size_t solidPoints = 0;
  std::size_t voidPoints = 0;
  for (int ix = 0; ix < steps; ++ix)
  {
    for (int iy = 0; iy < steps; ++iy)
    {
      for (int iz = 0; iz < steps; ++iz)
      {
        double3 fractional((ix + 0.5) / steps, (iy + 0.5) / steps, (iz + 0.5) / steps);
        double3 point = box.cell * fractional;

        double clearance = clearanceAt(box, fractionalPositions, inflatedRadii, point);
        PointClassification classification = classifier.accessibility.classify(point);

        ASSERT_FALSE(classification.resample) << "point " << ix << "," << iy << "," << iz << " was undecidable";

        // Points sitting on a surface may fall either way; everything else is decided.
        if (std::abs(clearance) < 1.0e-6) continue;
        EXPECT_EQ(classification.inside, clearance < 0.0)
            << "point " << ix << "," << iy << "," << iz << " at clearance " << clearance;

        if (clearance < 0.0)
          ++solidPoints;
        else
          ++voidPoints;
      }
    }
  }

  // Guard against a vacuous test: the cell must hold both solid and void.
  EXPECT_GT(solidPoints, 0u);
  EXPECT_GT(voidPoints, 0u);
}

// With every radius the same the Apollonius diagram is the ordinary Voronoi diagram, which is what the
// radical diagram also becomes, so the two classifiers are of the same network and the same sampled
// area and volume must come out of both, to within the noise of the sampling.
TEST(apollonius_pore_analysis, sampled_area_and_volume_match_the_radical_ones_for_equal_radii)
{
  double a = 12.0;
  SimulationBox box(a, a, a);
  const double radius = 1.6;
  const double probe = 0.35;

  RandomNumber random{std::optional<std::size_t>(11)};
  std::vector<double3> fractionalPositions;
  std::vector<double> radii;
  for (std::size_t i = 0; i < 14; ++i)
  {
    fractionalPositions.push_back(double3(random.uniform(), random.uniform(), random.uniform()));
    radii.push_back(radius);
  }

  ApolloniusAccessibility apollonius = ApolloniusAccessibility::create(box, fractionalPositions, radii, probe);
  ASSERT_TRUE(apollonius.diagram.networkIsComplete());
  PoreAccessibility radical = PoreAccessibility::create(box, fractionalPositions, radii, probe);

  VolumeSample apolloniusVolume = sampleAccessibleVolume(apollonius.accessibility, 200000);
  VolumeSample radicalVolume = sampleAccessibleVolume(radical, 200000);
  EXPECT_NEAR(apolloniusVolume.accessibleFraction, radicalVolume.accessibleFraction, 0.01);
  EXPECT_NEAR(apolloniusVolume.inaccessibleFraction, radicalVolume.inaccessibleFraction, 0.01);

  // Guard against a vacuous test: there has to be void to find.
  EXPECT_GT(apolloniusVolume.accessibleFraction + apolloniusVolume.inaccessibleFraction, 0.05);

  SurfaceAreaSample apolloniusArea = sampleAccessibleSurfaceArea(apollonius.accessibility, 50);
  SurfaceAreaSample radicalArea = sampleAccessibleSurfaceArea(radical, 50);
  double totalArea = radicalArea.accessible + radicalArea.inaccessible;
  EXPECT_GT(totalArea, 0.0);
  EXPECT_NEAR(apolloniusArea.accessible, radicalArea.accessible, 0.02 * totalArea);
  EXPECT_NEAR(apolloniusArea.inaccessible, radicalArea.inaccessible, 0.02 * totalArea);
}

// A window whose shape is known in advance: a ring of atoms on an ellipse of semi-axes 5 and 9, all in
// one plane, with the point to measure from at its centre. The free cross-section is then the inside of
// the ring, so the smallest free chord is across the narrow way and the largest across the wide way, and
// both are the ring's own semi-axes less an atom radius, doubled. The ring is closed, its atoms being
// closer together than their own diameter, so no direction escapes and nothing is cut off.
TEST(apollonius_pore_analysis, window_of_a_known_ring_comes_out_as_the_ring)
{
  double a = 60.0;
  SimulationBox box(a, a, a);
  const double atomRadius = 1.6;
  const std::size_t ringSize = 24;
  const double3 centre = box.cell * double3(0.5, 0.5, 0.5);

  auto ringNetwork = [&](double firstSemiAxis, double secondSemiAxis)
  {
    VoronoiNetwork network;
    network.simulationBox = box;
    for (std::size_t k = 0; k < ringSize; ++k)
    {
      double angle = 2.0 * std::numbers::pi * static_cast<double>(k) / static_cast<double>(ringSize);
      double3 position =
          centre + double3(firstSemiAxis * std::cos(angle), secondSemiAxis * std::sin(angle), 0.0);
      network.atomPositionsFractional.push_back(double3::fract(box.inverseCell * position));
      network.atomRadii.push_back(atomRadius);
    }
    return network;
  };

  VoronoiNetwork elliptical = ringNetwork(5.0, 9.0);
  PoreWindow window = PoreWindow::measure(elliptical, centre, double3(0.0, 0.0, 1.0));

  ASSERT_TRUE(window.measured);
  EXPECT_FALSE(window.clipped);

  // The ring has an atom on each of its own axes, so the narrow way across is exact: the tightest chord
  // runs straight at the two nearest atoms. The wide way is a little wider than the ring's own axis,
  // since a chord aimed between two neighbouring atoms passes where their surfaces cross rather than
  // where either surface is nearest, and that is farther out.
  EXPECT_NEAR(window.freeRadius, 5.0 - atomRadius, 1.0e-9);
  EXPECT_NEAR(window.smallestFreeWidth, 2.0 * (5.0 - atomRadius), 1.0e-6);
  EXPECT_GE(window.largestFreeWidth, 2.0 * (9.0 - atomRadius) - 1.0e-6);
  EXPECT_LE(window.largestFreeWidth, 2.0 * (9.0 - atomRadius) + 0.3);

  // The inscribed ellipse is the ring itself, up to the scalloping a ring of discrete atoms leaves.
  EXPECT_NEAR(window.minorAxis, 2.0 * (5.0 - atomRadius), 0.5);
  EXPECT_NEAR(window.majorAxis, 2.0 * (9.0 - atomRadius), 0.7);
  EXPECT_GE(window.majorAxis, window.minorAxis);
  EXPECT_NEAR(std::abs(double3::dot(window.majorAxisDirection, double3(0.0, 1.0, 0.0))), 1.0, 1.0e-2);

  // Every atom of the ring bounds the window; a few may be missed between sampled directions.
  EXPECT_GE(window.boundingAtoms, 20u);
  EXPECT_LE(window.boundingAtoms, ringSize);

  // A round ring is the control: the two numbers must collapse onto the one Df already gives.
  VoronoiNetwork circular = ringNetwork(6.0, 6.0);
  PoreWindow round = PoreWindow::measure(circular, centre, double3(0.0, 0.0, 1.0));

  ASSERT_TRUE(round.measured);
  EXPECT_NEAR(round.freeRadius, 6.0 - atomRadius, 1.0e-9);
  EXPECT_NEAR(round.smallestFreeWidth, 2.0 * round.freeRadius, 1.0e-6);
  EXPECT_GE(round.largestFreeWidth, 2.0 * round.freeRadius - 1.0e-6);
  EXPECT_LE(round.largestFreeWidth, 2.0 * round.freeRadius + 0.4);
  EXPECT_NEAR(round.minorAxis, 2.0 * round.freeRadius, 0.2);
  EXPECT_NEAR(round.majorAxis, round.minorAxis, 0.3);
}

// The window of a channel, taken from the diagram rather than handed in. The simple-cubic lattice puts
// its bottleneck at the middle of a cube edge, where four atoms are equally near, so the tightest chord
// across the passage runs between two of them and is the free sphere Df itself. The other directions in
// that plane run down the corridors between the atoms and never meet one, which is what `clipped` is
// for: the plane across this bottleneck is not a window ringed by atoms but a crossing of channels.
TEST(apollonius_pore_analysis, channel_window_of_the_simple_cubic_bottleneck)
{
  double a = 5.0;
  double r = 1.0;
  SimulationBox box(a, a, a);

  ApolloniusPoreNetwork pores = ApolloniusPoreNetwork::create(box, {double3(0.0, 0.0, 0.0)}, {r});
  ChannelAnalysis channels = ChannelAnalysis::compute(pores.network, 0.5);
  std::vector<ChannelWindow> windows = channelWindows(pores.network, channels);

  ASSERT_EQ(windows.size(), 1u);
  const ChannelWindow& channel = windows[0];
  EXPECT_EQ(channel.dimensionality, 3);
  EXPECT_NEAR(channel.freeSphereDiameter, std::sqrt(2.0) * a - 2.0 * r, 1.0e-6);

  ASSERT_TRUE(channel.window.measured);
  // The clearance is recomputed here from every atom, so it agreeing with the arc's bottleneck says the
  // diagram put the bottleneck where it says the bottleneck is.
  EXPECT_NEAR(2.0 * channel.window.freeRadius, channel.freeSphereDiameter, 1.0e-6);
  EXPECT_NEAR(channel.window.smallestFreeWidth, channel.freeSphereDiameter, 1.0e-6);
  EXPECT_TRUE(channel.window.clipped);
  EXPECT_GE(channel.window.boundingAtoms, 4u);
}

// The orderings that hold of any window: Df is one of the chords through the bottleneck, so it is the
// shortest of them at most, and the ellipse cannot be wider one way than the other way round. Checked
// over a random structure, where the channels are of no particular shape.
TEST(apollonius_pore_analysis, channel_windows_bracket_the_free_sphere)
{
  double a = 12.0;
  SimulationBox box(a, a, a);

  RandomNumber random{std::optional<std::size_t>(7)};
  std::vector<double3> fractionalPositions;
  std::vector<double> radii;
  for (std::size_t i = 0; i < 25; ++i)
  {
    fractionalPositions.push_back(double3(random.uniform(), random.uniform(), random.uniform()));
    radii.push_back(0.8 + 0.9 * random.uniform());
  }

  ApolloniusPoreNetwork pores = ApolloniusPoreNetwork::create(box, fractionalPositions, radii);
  ASSERT_TRUE(pores.networkIsComplete());
  ChannelAnalysis channels = ChannelAnalysis::compute(pores.network, 1.0);
  std::vector<ChannelWindow> windows = channelWindows(pores.network, channels);
  ASSERT_FALSE(windows.empty());

  PoreDiameters diameters = PoreDiameters::compute(pores.network);
  for (const ChannelWindow& channel : windows)
  {
    // A channel is a part of the structure, so its own free sphere is at most the structure's.
    EXPECT_LE(channel.freeSphereDiameter, diameters.freeSphereDiameter + 1.0e-6);
    ASSERT_TRUE(channel.window.measured);
    EXPECT_NEAR(2.0 * channel.window.freeRadius, channel.freeSphereDiameter, 1.0e-6);
    EXPECT_GE(channel.window.smallestFreeWidth, channel.freeSphereDiameter - 1.0e-6);
    EXPECT_GE(channel.window.largestFreeWidth, channel.window.smallestFreeWidth - 1.0e-9);
    EXPECT_GE(channel.window.majorAxis, channel.window.minorAxis - 1.0e-9);
    EXPECT_LE(channel.window.majorAxis, channel.window.largestFreeWidth + 1.0e-6);
    EXPECT_NEAR(channel.window.normal.length(), 1.0, 1.0e-9);
  }
}

// Accuracy and cost of the Apollonius diagram against the radical network it replaces, on a directory
// of P1 CIFs, with the workload zeo++ is given for the same structures:
//   network -res -chan 1.2 -sa 1.2 1.2 2000 -vol 1.2 1.2 50000 <cif>
// Both backends are sampled with the same seed and the same points, so a difference in area or volume
// is a difference in what the two diagrams call accessible, not sampling noise. Disabled by default.
TEST(apollonius_pore_analysis, DISABLED_compare_p1_cifs_apollonius_vs_radical)
{
  const std::string directory = "/Users/dubbelda/Research/Codes/new-zeo";
  const double probeRadius = 1.2;
  const std::size_t surfaceSamplesPerAtom = 2000;
  const std::size_t volumeSamples = 50000;

  ForceField forceField = ForceField::makeZeoPlusPlusForceField(12.0, true, false, false);

  std::vector<std::string> files;
  for (const auto& entry : std::filesystem::directory_iterator(directory))
    if (entry.path().extension() == ".cif" && entry.path().stem().string().ends_with("-P1"))
      files.push_back(entry.path().string());
  std::sort(files.begin(), files.end());
  ASSERT_FALSE(files.empty());

  auto clock = []() { return std::chrono::steady_clock::now(); };
  auto ms = [](auto d) { return std::chrono::duration<double, std::milli>(d).count(); };

  // The area and volume estimates, over whichever classifier is handed in, from the same points.
  auto sample = [&](const PoreAccessibility& accessibility, const SimulationBox& box)
  {
    RandomNumber random{std::optional<std::size_t>(42)};
    double accessibleArea = 0.0;
    double inaccessibleArea = 0.0;
    for (std::size_t i = 0; i < accessibility.atomPositions.size(); ++i)
    {
      double inflatedRadius = accessibility.atomRadii[i];
      double sphereArea = 4.0 * std::numbers::pi * inflatedRadius * inflatedRadius;
      std::size_t accessibleCount = 0;
      std::size_t inaccessibleCount = 0;
      for (std::size_t s = 0; s < surfaceSamplesPerAtom; ++s)
      {
        double3 point = accessibility.atomPositions[i] + inflatedRadius * random.randomVectorOnUnitSphere();
        if (accessibility.overlapsAtom(point, i)) continue;
        PointClassification classification = accessibility.classify(point);
        if (classification.resample || classification.inside) continue;
        if (classification.accessible)
          ++accessibleCount;
        else
          ++inaccessibleCount;
      }
      double perSample = sphereArea / static_cast<double>(surfaceSamplesPerAtom);
      accessibleArea += static_cast<double>(accessibleCount) * perSample;
      inaccessibleArea += static_cast<double>(inaccessibleCount) * perSample;
    }

    std::size_t accessibleVolumeCount = 0;
    std::size_t inaccessibleVolumeCount = 0;
    for (std::size_t s = 0; s < volumeSamples; ++s)
    {
      double3 point = box.cell * double3(random.uniform(), random.uniform(), random.uniform());
      PointClassification classification = accessibility.classify(point);
      if (classification.inside || classification.resample) continue;
      if (classification.accessible)
        ++accessibleVolumeCount;
      else
        ++inaccessibleVolumeCount;
    }
    double volume = box.volume;
    return std::tuple<double, double, double, double>{
        accessibleArea, inaccessibleArea,
        static_cast<double>(accessibleVolumeCount) / static_cast<double>(volumeSamples) * volume,
        static_cast<double>(inaccessibleVolumeCount) / static_cast<double>(volumeSamples) * volume};
  };

  for (const std::string& file : files)
  {
    std::string name = std::filesystem::path(file).stem().string();

    // The (uniform) P1 CIF: cell parameters by keyword, atoms as 7-column rows
    // (label, type_symbol, forcefield_label, fract_x, fract_y, fract_z, charge).
    std::ifstream stream(file);
    ASSERT_TRUE(stream.good()) << "cannot open " << file;
    double a{}, b{}, c{}, alpha{}, beta{}, gamma{};
    std::vector<double3> fractionalPositions;
    std::vector<double> radii;
    std::string line;
    while (std::getline(stream, line))
    {
      std::istringstream iss(line);
      std::string first;
      if (!(iss >> first)) continue;
      if (first == "_cell_length_a") iss >> a;
      else if (first == "_cell_length_b") iss >> b;
      else if (first == "_cell_length_c") iss >> c;
      else if (first == "_cell_angle_alpha") iss >> alpha;
      else if (first == "_cell_angle_beta") iss >> beta;
      else if (first == "_cell_angle_gamma") iss >> gamma;
      else if (first[0] != '_' && first != "loop_" && !first.starts_with("data_"))
      {
        std::string element, label;
        double x, y, z, charge;
        if (iss >> element >> label >> x >> y >> z >> charge)
        {
          std::optional<std::size_t> pseudoType = forceField.findPseudoAtom(element);
          ASSERT_TRUE(pseudoType.has_value()) << "unknown element " << element << " in " << name;
          fractionalPositions.push_back(double3(x, y, z));
          radii.push_back(0.5 * forceField(pseudoType.value(), pseudoType.value()).sizeParameter());
        }
      }
    }
    ASSERT_FALSE(fractionalPositions.empty()) << name;

    SimulationBox::Type type =
        (std::abs(alpha - 90.0) > 1.0e-3 || std::abs(beta - 90.0) > 1.0e-3 || std::abs(gamma - 90.0) > 1.0e-3)
            ? SimulationBox::Type::Triclinic
            : SimulationBox::Type::Rectangular;
    SimulationBox box(a, b, c, alpha * std::numbers::pi / 180.0, beta * std::numbers::pi / 180.0,
                      gamma * std::numbers::pi / 180.0, type);

    for (const std::string& backend : {std::string("radical"), std::string("apollonius")})
    {
      bool useApollonius = (backend == "apollonius");

      auto t0 = clock();
      VoronoiNetwork network = useApollonius ? ApolloniusPoreNetwork::create(box, fractionalPositions, radii).network
                                             : VoronoiNetwork::create(box, fractionalPositions, radii);
      auto t1 = clock();
      PoreDiameters diameters = PoreDiameters::compute(network);
      ChannelAnalysis channels = ChannelAnalysis::compute(network, probeRadius);
      auto t2 = clock();

      // The area and volume are read off a second diagram, of the atoms inflated by the probe.
      PoreAccessibility accessibility =
          useApollonius
              ? ApolloniusAccessibility::create(box, fractionalPositions, radii, probeRadius).accessibility
              : PoreAccessibility::create(box, fractionalPositions, radii, probeRadius);
      auto t3 = clock();
      auto [asa, nasa, av, nav] = sample(accessibility, box);
      auto t4 = clock();

      std::string dims;
      for (const VoronoiPore& pore : channels.pores)
        if (pore.isChannel) dims += std::format("{}", pore.dimensionality);

      std::cout << std::format(
          "RESULT {} {} atoms= {} Di= {:.5f} Df= {:.5f} Dif= {:.5f} chan= {} dims= [{}] "
          "ASA= {:.2f} NASA= {:.2f} AV= {:.2f} NAV= {:.2f} "
          "t_net= {:.1f} t_chan= {:.1f} t_acc= {:.1f} t_sample= {:.1f} t_total= {:.1f}\n",
          name, backend, fractionalPositions.size(), diameters.includedSphereDiameter, diameters.freeSphereDiameter,
          diameters.includedAlongFreePathDiameter, channels.numberOfChannels, dims, asa, nasa, av, nav, ms(t1 - t0),
          ms(t2 - t1), ms(t3 - t2), ms(t4 - t3), ms(t4 - t0));
      std::cout.flush();
    }
  }
}
