#include <gtest/gtest.h>

import std;

import int3;
import double3;
import simulationbox;
import atom;
import randomnumbers;
import forcefield;
import framework;
import voronoi_network;
import voronoi_pore_diameters;
import voronoi_channels;
import pore_accessibility;

// Simple-cubic lattice of one atom per cell (lattice a, atom radius r). The Voronoi cell
// is a cube; all corners are periodic images of a single node at the cube centre.
//   Di  = 2*(sqrt(3)/2 a - r) = sqrt(3) a - 2r   (vertex is farthest point)
//   Df  = 2*(a/sqrt(2)   - r) = sqrt(2) a - 2r    (edge midpoint is the bottleneck)
//   Dif = Di                                       (only one node lies on the path)
TEST(voronoi_analysis, simple_cubic_pore_diameters)
{
  double a = 5.0;
  double r = 1.0;
  SimulationBox box(a, a, a);
  VoronoiNetwork network = VoronoiNetwork::create(box, {double3(0.0, 0.0, 0.0)}, {r});

  EXPECT_EQ(network.nodes.size(), 1);

  PoreDiameters diameters = PoreDiameters::compute(network);
  EXPECT_NEAR(diameters.includedSphereDiameter, std::sqrt(3.0) * a - 2.0 * r, 1.0e-6);
  EXPECT_NEAR(diameters.freeSphereDiameter, std::sqrt(2.0) * a - 2.0 * r, 1.0e-6);
  EXPECT_NEAR(diameters.includedAlongFreePathDiameter, std::sqrt(3.0) * a - 2.0 * r, 1.0e-6);

  // Ordering that must always hold.
  EXPECT_GE(diameters.includedSphereDiameter, diameters.freeSphereDiameter);
  EXPECT_GE(diameters.includedAlongFreePathDiameter, diameters.freeSphereDiameter);
}

// The connected void network of a simple-cubic lattice percolates in all three directions.
TEST(voronoi_analysis, simple_cubic_is_three_dimensional_channel)
{
  double a = 5.0;
  SimulationBox box(a, a, a);
  VoronoiNetwork network = VoronoiNetwork::create(box, {double3(0.0, 0.0, 0.0)}, {1.0});

  ChannelAnalysis channels = ChannelAnalysis::compute(network, 0.5);
  EXPECT_EQ(channels.numberOfChannels, 1);
  EXPECT_EQ(channels.numberOfPockets, 0);
  ASSERT_EQ(channels.pores.size(), 1);
  EXPECT_TRUE(channels.pores[0].isChannel);
  EXPECT_EQ(channels.pores[0].dimensionality, 3);
}

TEST(voronoi_analysis, lattice_vector_rank)
{
  EXPECT_EQ(latticeVectorRank({}), 0);
  EXPECT_EQ(latticeVectorRank({int3(1, 0, 0)}), 1);
  EXPECT_EQ(latticeVectorRank({int3(1, 0, 0), int3(2, 0, 0)}), 1);
  EXPECT_EQ(latticeVectorRank({int3(1, 0, 0), int3(0, 1, 0)}), 2);
  EXPECT_EQ(latticeVectorRank({int3(1, 0, 0), int3(0, 1, 0), int3(1, 1, 0)}), 2);
  EXPECT_EQ(latticeVectorRank({int3(1, 0, 0), int3(0, 1, 0), int3(0, 0, 1)}), 3);
}

// Accessibility classifier on the simple-cubic system: the atom interior is solid, the
// cage centre is accessible void.
TEST(voronoi_analysis, accessibility_classification)
{
  double a = 5.0;
  double r = 1.0;
  double probe = 0.5;
  SimulationBox box(a, a, a);
  PoreAccessibility accessibility = PoreAccessibility::create(box, {double3(0.0, 0.0, 0.0)}, {r}, probe);

  // Point well inside the inflated atom.
  PointClassification insidePoint = accessibility.classify(double3(0.2, 0.0, 0.0));
  EXPECT_TRUE(insidePoint.inside);

  // Cage centre: void and accessible (network percolates in 3D).
  PointClassification centre = accessibility.classify(double3(0.5 * a, 0.5 * a, 0.5 * a));
  EXPECT_FALSE(centre.inside);
  if (!centre.resample)
  {
    EXPECT_TRUE(centre.accessible);
  }
}

// The Apollonius vertex primitive: the sphere tangent to and outside four given spheres.
TEST(voronoi_analysis, apollonius_tangent_sphere)
{
  // A regular tetrahedron of equal spheres. The tangent sphere sits at the circumcentre, here
  // the origin, with radius equal to the circumradius less the common atom radius.
  std::array<double3, 4> centres{double3(1.0, 1.0, 1.0), double3(1.0, -1.0, -1.0), double3(-1.0, 1.0, -1.0),
                                 double3(-1.0, -1.0, 1.0)};
  std::array<double, 4> equalRadii{0.5, 0.5, 0.5, 0.5};

  std::vector<ApolloniusSphere> spheres = apolloniusTangentSpheres(centres, equalRadii);
  ASSERT_EQ(spheres.size(), 1u);
  EXPECT_NEAR(spheres[0].radius, std::sqrt(3.0) - 0.5, 1.0e-12);
  EXPECT_NEAR(spheres[0].centre.length(), 0.0, 1.0e-12);

  // Unequal radii have no simple closed-form reference, so check the defining property instead:
  // the solution must touch every one of the four spheres.
  std::array<double, 4> unequalRadii{0.9, 0.4, 0.6, 0.25};
  std::vector<ApolloniusSphere> weighted = apolloniusTangentSpheres(centres, unequalRadii);
  EXPECT_FALSE(weighted.empty());
  for (const ApolloniusSphere& sphere : weighted)
  {
    for (std::size_t i = 0; i < 4; ++i)
    {
      EXPECT_NEAR((sphere.centre - centres[i]).length(), unequalRadii[i] + sphere.radius, 1.0e-9);
    }
  }

  // Coplanar centres leave the tangent sphere undetermined along the normal direction.
  std::array<double3, 4> coplanar{double3(0.0, 0.0, 0.0), double3(1.0, 0.0, 0.0), double3(0.0, 1.0, 0.0),
                                  double3(1.0, 1.0, 0.0)};
  EXPECT_TRUE(apolloniusTangentSpheres(coplanar, equalRadii).empty());
}

// The largest included sphere must be the true peak of the clearance field. A radical Voronoi
// vertex is only the peak when the surrounding radii are equal; with heterogeneous radii the peak
// moves off the vertex to an Apollonius vertex. Validated against a brute-force scan of the
// clearance field, which the network must match rather than merely bound.
TEST(voronoi_analysis, largest_included_sphere_matches_brute_force_clearance_peak)
{
  double a = 14.0;
  SimulationBox box(a, a, a);

  std::mt19937 generator(11);
  std::uniform_real_distribution<double> uniform(0.0, 1.0);
  std::uniform_real_distribution<double> radiusDistribution(0.4, 2.2);

  std::vector<double3> fractionalPositions;
  std::vector<double> radii;
  for (std::size_t i = 0; i < 40; ++i)
  {
    fractionalPositions.push_back(double3(uniform(generator), uniform(generator), uniform(generator)));
    radii.push_back(radiusDistribution(generator));
  }

  auto clearanceAt = [&](const double3& fractional)
  {
    double clearance = std::numeric_limits<double>::max();
    for (std::size_t j = 0; j < fractionalPositions.size(); ++j)
    {
      double3 delta = fractionalPositions[j] - fractional;
      delta -= double3(std::round(delta.x), std::round(delta.y), std::round(delta.z));
      clearance = std::min(clearance, (box.cell * delta).length() - radii[j]);
    }
    return clearance;
  };

  // Reference peak by exhaustive scan followed by a shrinking local search from each of the best
  // grid points. Several starts are needed because the clearance field has many basins.
  std::size_t resolution = 60;
  std::vector<std::pair<double, double3>> scan;
  for (std::size_t ix = 0; ix < resolution; ++ix)
    for (std::size_t iy = 0; iy < resolution; ++iy)
      for (std::size_t iz = 0; iz < resolution; ++iz)
      {
        double3 fractional(static_cast<double>(ix) / resolution, static_cast<double>(iy) / resolution,
                           static_cast<double>(iz) / resolution);
        scan.emplace_back(clearanceAt(fractional), fractional);
      }
  std::size_t starts = 64;
  std::partial_sort(scan.begin(), scan.begin() + static_cast<std::ptrdiff_t>(starts), scan.end(),
                    [](const auto& lhs, const auto& rhs) { return lhs.first > rhs.first; });
  double gridClearance = scan.front().first;

  double bestClearance = -std::numeric_limits<double>::max();
  for (std::size_t start = 0; start < starts; ++start)
  {
    double3 current = scan[start].second;
    double currentClearance = scan[start].first;
    for (double step = 1.0 / static_cast<double>(resolution); step > 1.0e-9; step *= 0.5)
    {
      bool improved = true;
      while (improved)
      {
        improved = false;
        for (int ix = -1; ix <= 1; ++ix)
          for (int iy = -1; iy <= 1; ++iy)
            for (int iz = -1; iz <= 1; ++iz)
            {
              double3 trial = current + step * double3(ix, iy, iz);
              double clearance = clearanceAt(trial);
              if (clearance > currentClearance)
              {
                currentClearance = clearance;
                current = trial;
                improved = true;
              }
            }
      }
    }
    bestClearance = std::max(bestClearance, currentClearance);
  }

  VoronoiNetwork network = VoronoiNetwork::create(box, fractionalPositions, radii);
  ASSERT_FALSE(network.nodes.empty());

  double networkRadius = 0.0;
  double clearanceRadius = 0.0;
  for (const VoronoiNode& node : network.nodes)
  {
    networkRadius = std::max(networkRadius, node.maximalRadius);
    clearanceRadius = std::max(clearanceRadius, node.radius);
    EXPECT_GE(node.maximalRadius, node.radius - 1.0e-9);

    // The decisive check: the reported sphere must actually exist. Measuring the clearance at the
    // reported centre independently proves the radius is realised and not merely asserted, so no
    // node can claim more room than the structure has.
    double3 fractional = box.inverseCell * node.maximalPosition;
    EXPECT_NEAR(clearanceAt(fractional), node.maximalRadius, 1.0e-9);
  }

  // The refinement must earn its place: the radical vertices fall short of the peak, and the
  // Apollonius ascent must close the gap rather than merely narrow it.
  EXPECT_LT(clearanceRadius, gridClearance) << "structure does not exercise the refinement";
  EXPECT_GE(networkRadius, bestClearance - 1.0e-6);
  EXPECT_NEAR(network.largestIncludedSphereDiameter(), 2.0 * networkRadius, 1.0e-12);
}

// Node and edge radii must be true clearances, min over every atom of |x - x_j| - r_j, not just
// over the atoms bounding the feature. Radical adjacency is decided by the power distance, which
// orders atoms differently once the radii differ, so a non-bounding atom can be the nearest in
// clearance. Checked against brute force over all atoms and their 27 nearest images.
TEST(voronoi_analysis, radii_are_true_clearances_for_heterogeneous_radii)
{
  double a = 14.0;
  SimulationBox box(a, a, a);

  std::mt19937 generator(7);
  std::uniform_real_distribution<double> uniform(0.0, 1.0);
  std::uniform_real_distribution<double> radiusDistribution(0.4, 2.2);

  std::vector<double3> fractionalPositions;
  std::vector<double> radii;
  for (std::size_t i = 0; i < 40; ++i)
  {
    fractionalPositions.push_back(double3(uniform(generator), uniform(generator), uniform(generator)));
    radii.push_back(radiusDistribution(generator));
  }
  ASSERT_GT(*std::max_element(radii.begin(), radii.end()) - *std::min_element(radii.begin(), radii.end()), 1.0)
      << "radii must be spread out enough for power and clearance ordering to differ";

  VoronoiNetwork network = VoronoiNetwork::create(box, fractionalPositions, radii);
  ASSERT_FALSE(network.nodes.empty());
  ASSERT_FALSE(network.edges.empty());

  auto clearanceAt = [&](const double3& point)
  {
    double clearance = std::numeric_limits<double>::max();
    for (std::size_t j = 0; j < fractionalPositions.size(); ++j)
    {
      for (int ix = -1; ix <= 1; ++ix)
        for (int iy = -1; iy <= 1; ++iy)
          for (int iz = -1; iz <= 1; ++iz)
          {
            double3 image = box.cell * (fractionalPositions[j] + double3(static_cast<double>(ix),
                                                                         static_cast<double>(iy),
                                                                         static_cast<double>(iz)));
            clearance = std::min(clearance, (point - image).length() - radii[j]);
          }
    }
    return clearance;
  };

  // Overestimating is the error that matters: it claims more room than exists. Undershooting by
  // up to the node merge tolerance is expected and safe, because a node keeps the smallest
  // contribution of the vertices merged into it while storing only one of their positions.
  const double mergeTolerance = 0.02;
  for (const VoronoiNode& node : network.nodes)
  {
    double clearance = clearanceAt(node.position);
    EXPECT_LE(node.radius, clearance + 1.0e-9) << "node radius claims more room than the true clearance";
    EXPECT_GE(node.radius, clearance - mergeTolerance) << "node radius is below the merge tolerance band";
  }

  // An edge radius is the tightest clearance anywhere on the segment, so sampling the segment
  // may only ever find room, never less than the stored value.
  for (const VoronoiEdge& edge : network.edges)
  {
    double3 from = network.nodes[edge.from].position;
    double3 to = network.nodes[edge.to].position + box.cell * double3(static_cast<double>(edge.delta.x),
                                                                      static_cast<double>(edge.delta.y),
                                                                      static_cast<double>(edge.delta.z));
    double tightest = std::numeric_limits<double>::max();
    const std::size_t samples = 64;
    for (std::size_t s = 0; s <= samples; ++s)
    {
      double t = static_cast<double>(s) / static_cast<double>(samples);
      tightest = std::min(tightest, clearanceAt(from + t * (to - from)));
    }
    EXPECT_GE(tightest, edge.radius - 1.0e-6) << "a sphere of the stored radius does not fit along the edge";
  }
}

// The network is a radical (power) diagram, in which a small atom engulfed by a larger one has
// a completely empty cell. Selecting the cell to test by Euclidean nearest atom can then land
// on a site with no vertices, leaving the point undecidable; callers discard such points, which
// biases sampled areas and volumes low. Selecting by power distance always lands in a non-empty
// cell, so no point may come back flagged for resampling.
TEST(voronoi_analysis, classify_decides_every_point_with_heterogeneous_radii)
{
  double a = 12.0;
  SimulationBox box(a, a, a);
  const double largeRadius = 2.2;
  const double smallRadius = 0.35;
  const double probe = 0.4;

  // A small atom midway between two large ones, each at separation s. The radical plane
  // against a large neighbour sits at (s^2 + r^2 - R^2) / 2s from the small atom, measured
  // towards that neighbour, so once s < sqrt(R^2 - r^2) that offset is negative and the cell
  // is pushed to the far side. Two such neighbours on opposite sides impose contradictory
  // constraints and the power cell is empty. Radii are inflated by the probe, so the bound is
  // sqrt(2.6^2 - 0.75^2) = 2.49 A and a separation of 2 A qualifies.
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

  PoreAccessibility accessibility = PoreAccessibility::create(box, fractionalPositions, radii, probe);

  // Guard against the test going vacuous: the configuration must really produce empty cells.
  std::size_t atomsWithoutNodes = 0;
  for (const auto& nodeVectors : accessibility.network.atomNodeVectors)
    if (nodeVectors.empty()) ++atomsWithoutNodes;
  ASSERT_GT(atomsWithoutNodes, 0u) << "configuration does not produce an empty power cell";

  RandomNumber random{std::optional<std::size_t>(42)};
  for (std::size_t s = 0; s < 20000; ++s)
  {
    double3 point = box.cell * double3(random.uniform(), random.uniform(), random.uniform());
    PointClassification classification = accessibility.classify(point);
    ASSERT_FALSE(classification.resample) << "sample " << s << " was left undecidable";
  }
}

// End-to-end on a real zeolite (ITQ-29 / LTA): the full pipeline should build a network,
// find a percolating channel system, and produce a consistent Di >= Df ordering.
TEST(voronoi_analysis, itq29_pipeline)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  Framework framework = Framework::makeITQ29(forceField, int3(1, 1, 1));

  std::vector<double3> fractionalPositions;
  std::vector<double> radii;
  for (const Atom& atom : framework.unitCellAtoms)
  {
    fractionalPositions.push_back(framework.simulationBox.inverseCell * atom.position);
    std::size_t type = static_cast<std::size_t>(atom.type);
    radii.push_back(0.5 * forceField(type, type).sizeParameter());
  }

  VoronoiNetwork network = VoronoiNetwork::create(framework.simulationBox, fractionalPositions, radii);
  EXPECT_GT(network.nodes.size(), 0);
  EXPECT_GT(network.edges.size(), 0);

  PoreDiameters diameters = PoreDiameters::compute(network);
  EXPECT_GT(diameters.includedSphereDiameter, 0.0);
  EXPECT_GT(diameters.freeSphereDiameter, 0.0);
  EXPECT_GE(diameters.includedSphereDiameter, diameters.freeSphereDiameter - 1.0e-9);
  EXPECT_GE(diameters.includedAlongFreePathDiameter, diameters.freeSphereDiameter - 1.0e-9);

  // A modest probe still percolates through the LTA pore system in 3D.
  ChannelAnalysis channels = ChannelAnalysis::compute(network, 1.0);
  std::size_t maxDimensionality = 0;
  for (const VoronoiPore& pore : channels.pores)
    if (pore.isChannel) maxDimensionality = std::max<std::size_t>(maxDimensionality, static_cast<std::size_t>(pore.dimensionality));
  EXPECT_GE(channels.numberOfChannels, 1);
  EXPECT_EQ(maxDimensionality, 3);
}

// Timing of the full pipeline on FAU (all-silica Y), matching the workload of the zeo++
// runs on FAU_SI.cif (probe 1.2 Å, 2000 surface samples/atom, 50000 volume samples).
// Disabled by default; run with --gtest_also_run_disabled_tests.
TEST(voronoi_analysis, DISABLED_fau_timing_vs_zeopp)
{
  // Same radii table as zeo++ (CLI option --zeo++): radius = sigma / 2.
  ForceField forceField = ForceField::makeZeoPlusPlusForceField(12.0, true, false, false);
  Framework framework = Framework::makeFAU(forceField, int3(1, 1, 1));

  std::vector<double3> fractionalPositions;
  std::vector<double> radii;
  for (const Atom& atom : framework.unitCellAtoms)
  {
    fractionalPositions.push_back(framework.simulationBox.inverseCell * atom.position);
    std::size_t type = static_cast<std::size_t>(atom.type);
    radii.push_back(0.5 * forceField(type, type).sizeParameter());
  }
  double probeRadius = 1.2;

  auto clock = []() { return std::chrono::steady_clock::now(); };
  auto ms = [](auto d) { return std::chrono::duration<double, std::milli>(d).count(); };

  auto t0 = clock();
  VoronoiNetwork network = VoronoiNetwork::create(framework.simulationBox, fractionalPositions, radii);
  auto t1 = clock();
  PoreDiameters diameters = PoreDiameters::compute(network);
  auto t2 = clock();
  ChannelAnalysis channels = ChannelAnalysis::compute(network, probeRadius);
  auto t3 = clock();
  PoreAccessibility accessibility =
      PoreAccessibility::create(framework.simulationBox, fractionalPositions, radii, probeRadius);
  auto t4 = clock();

  RandomNumber random{std::nullopt};
  std::size_t samplesPerAtom = 2000;
  for (std::size_t i = 0; i < accessibility.atomPositions.size(); ++i)
  {
    double inflated = accessibility.atomRadii[i];
    for (std::size_t s = 0; s < samplesPerAtom; ++s)
    {
      double3 point = accessibility.atomPositions[i] + inflated * random.randomVectorOnUnitSphere();
      volatile auto c = accessibility.classify(point);
      (void)c;
    }
  }
  auto t5 = clock();
  std::size_t volumeSamples = 50000;
  for (std::size_t s = 0; s < volumeSamples; ++s)
  {
    double3 point = framework.simulationBox.cell * double3(random.uniform(), random.uniform(), random.uniform());
    volatile auto c = accessibility.classify(point);
    (void)c;
  }
  auto t6 = clock();

  std::cout << "FAU atoms: " << framework.unitCellAtoms.size() << ", Voronoi nodes: " << network.nodes.size() << "\n";
  std::cout << "Di/Df/Dif [Å]: " << diameters.includedSphereDiameter << " " << diameters.freeSphereDiameter << " "
            << diameters.includedAlongFreePathDiameter << "\n";
  std::cout << std::format("network build (-res core): {:.1f} ms\n", ms(t1 - t0));
  std::cout << std::format("pore diameters:            {:.1f} ms\n", ms(t2 - t1));
  std::cout << std::format("channel analysis (-chan):  {:.1f} ms\n", ms(t3 - t2));
  std::cout << std::format("accessibility setup:       {:.1f} ms\n", ms(t4 - t3));
  std::cout << std::format("surface area (-sa 2000):   {:.1f} ms\n", ms(t5 - t4));
  std::cout << std::format("accessible volume (-vol):  {:.1f} ms\n", ms(t6 - t5));
  std::cout << std::format("TOTAL:                     {:.1f} ms\n", ms(t6 - t0));
}

// Dump the RASPA3 Voronoi network for one P1 CIF (path in RASPA_DUMP_CIF) in the same
// text layout as zeo++'s `-nt2` output, for direct node/edge comparison.
TEST(voronoi_analysis, DISABLED_dump_network_nt2)
{
  const char* cifPath = std::getenv("RASPA_DUMP_CIF");
  const char* outPath = std::getenv("RASPA_DUMP_OUT");
  ASSERT_NE(cifPath, nullptr);
  ASSERT_NE(outPath, nullptr);

  ForceField forceField = ForceField::makeZeoPlusPlusForceField(12.0, true, false, false);

  std::ifstream stream(cifPath);
  ASSERT_TRUE(stream.good());
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
        ASSERT_TRUE(pseudoType.has_value());
        fractionalPositions.push_back(double3(x, y, z));
        radii.push_back(0.5 * forceField(pseudoType.value(), pseudoType.value()).sizeParameter());
      }
    }
  }
  SimulationBox::Type type = (std::abs(alpha - 90.0) > 1.0e-3 || std::abs(beta - 90.0) > 1.0e-3 ||
                              std::abs(gamma - 90.0) > 1.0e-3)
                                 ? SimulationBox::Type::Triclinic
                                 : SimulationBox::Type::Rectangular;
  SimulationBox simulationBox(a, b, c, alpha * std::numbers::pi / 180.0, beta * std::numbers::pi / 180.0,
                              gamma * std::numbers::pi / 180.0, type);

  VoronoiNetwork network = VoronoiNetwork::create(simulationBox, fractionalPositions, radii);

  std::ofstream out(outPath);
  out << "Vertex table:\n";
  for (std::size_t i = 0; i < network.nodes.size(); ++i)
  {
    const VoronoiNode& node = network.nodes[i];
    out << i << " " << node.position.x << " " << node.position.y << " " << node.position.z << " " << node.radius
        << "\n";
  }
  out << "Edge table:\n";
  for (const VoronoiEdge& edge : network.edges)
  {
    out << edge.from << " -> " << edge.to << " " << edge.radius << " " << edge.delta.x << " " << edge.delta.y << " "
        << edge.delta.z << " " << edge.length << "\n";
  }
}

// Accuracy + performance comparison against zeo++ 0.3 on a directory of P1 CIFs.
// Workload matches the zeo++ invocation:
//   network -res -chan 1.2 -sa 1.2 1.2 2000 -vol 1.2 1.2 50000 <cif>
// Prints one machine-parseable line per structure. Disabled by default.
TEST(voronoi_analysis, DISABLED_compare_p1_cifs_vs_zeopp)
{
  const std::string directory = "/Users/dubbelda/test_cifs/P1";
  const double probeRadius = 1.2;
  const std::size_t surfaceSamplesPerAtom = 2000;
  const std::size_t volumeSamples = 50000;

  ForceField forceField = ForceField::makeZeoPlusPlusForceField(12.0, true, false, false);

  std::vector<std::string> files;
  for (const auto& entry : std::filesystem::directory_iterator(directory))
    if (entry.path().extension() == ".cif") files.push_back(entry.path().string());
  std::sort(files.begin(), files.end());
  ASSERT_FALSE(files.empty());

  auto clock = []() { return std::chrono::steady_clock::now(); };
  auto ms = [](auto d) { return std::chrono::duration<double, std::milli>(d).count(); };

  for (const std::string& file : files)
  {
    std::string name = std::filesystem::path(file).stem().string();

    // Parse the (uniform) P1 CIF: cell parameters by keyword, atoms as 7-column rows
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

    SimulationBox::Type type = (std::abs(alpha - 90.0) > 1.0e-3 || std::abs(beta - 90.0) > 1.0e-3 ||
                                std::abs(gamma - 90.0) > 1.0e-3)
                                   ? SimulationBox::Type::Triclinic
                                   : SimulationBox::Type::Rectangular;
    SimulationBox simulationBox(a, b, c, alpha * std::numbers::pi / 180.0, beta * std::numbers::pi / 180.0,
                                gamma * std::numbers::pi / 180.0, type);

    auto t0 = clock();
    VoronoiNetwork network = VoronoiNetwork::create(simulationBox, fractionalPositions, radii);
    auto t1 = clock();
    PoreDiameters diameters = PoreDiameters::compute(network);
    ChannelAnalysis channels = ChannelAnalysis::compute(network, probeRadius);
    auto t2 = clock();
    PoreAccessibility accessibility =
        PoreAccessibility::create(simulationBox, fractionalPositions, radii, probeRadius);
    auto t3 = clock();

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
    auto t4 = clock();

    std::size_t accessibleVolumeCount = 0;
    std::size_t inaccessibleVolumeCount = 0;
    for (std::size_t s = 0; s < volumeSamples; ++s)
    {
      double3 point = simulationBox.cell * double3(random.uniform(), random.uniform(), random.uniform());
      PointClassification classification = accessibility.classify(point);
      if (classification.inside || classification.resample) continue;
      if (classification.accessible)
        ++accessibleVolumeCount;
      else
        ++inaccessibleVolumeCount;
    }
    auto t5 = clock();

    double volume = simulationBox.volume;
    double accessibleVolume = static_cast<double>(accessibleVolumeCount) / static_cast<double>(volumeSamples) * volume;
    double inaccessibleVolume =
        static_cast<double>(inaccessibleVolumeCount) / static_cast<double>(volumeSamples) * volume;

    std::string dims;
    for (const VoronoiPore& pore : channels.pores)
      if (pore.isChannel) dims += std::format("{}", pore.dimensionality);

    std::cout << std::format(
        "RESULT {} atoms= {} Di= {:.5f} Df= {:.5f} Dif= {:.5f} chan= {} dims= [{}] "
        "ASA= {:.2f} NASA= {:.2f} AV= {:.2f} NAV= {:.2f} "
        "t_net= {:.1f} t_chan= {:.1f} t_acc= {:.1f} t_sa= {:.1f} t_vol= {:.1f} t_total= {:.1f}\n",
        name, fractionalPositions.size(), diameters.includedSphereDiameter, diameters.freeSphereDiameter,
        diameters.includedAlongFreePathDiameter, channels.numberOfChannels, dims, accessibleArea, inaccessibleArea,
        accessibleVolume, inaccessibleVolume, ms(t1 - t0), ms(t2 - t1), ms(t3 - t2), ms(t4 - t3), ms(t5 - t4),
        ms(t5 - t0));
  }
}

// Full cross-check against zeo++ 0.3 on every CIF shipped with it, using the same radial
// (power) Voronoi diagram and the same radii table (--zeo++). Reports Di/Df/Dif, channel
// count and dimensionality, accessible/inaccessible surface area and accessible volume for
// a probe radius of 1.2 Å. Disabled by default; run with --gtest_also_run_disabled_tests.
TEST(voronoi_analysis, DISABLED_compare_all_cifs_vs_zeopp)
{
  // Read zeo++'s own CSSR expansion of each structure (produced with `network -cssr`) so
  // both codes operate on byte-for-byte identical atom sets and cells; this isolates the
  // Voronoi/power algorithm from any CIF symmetry-expansion differences. Radii are taken
  // from the same table (element -> sigma/2) as zeo++.
  const std::string directory = "/Users/dubbelda/Research/raspa3/.tmp_bench/";
  const std::vector<std::string> names{"AFX_SI", "DDR", "EMT", "ERI_SI", "FAU_SI", "LTA4A", "LTN"};
  const double probeRadius = 1.2;
  const std::size_t surfaceSamplesPerAtom = 5000;
  const std::size_t volumeSamples = 2000000;

  ForceField forceField = ForceField::makeZeoPlusPlusForceField(12.0, true, false, false);

  for (const std::string& name : names)
  {
    std::ifstream stream(directory + name + ".cssr");
    ASSERT_TRUE(stream.good()) << "cannot open " << name << ".cssr";

    // CSSR: line 1 = a b c; line 2 = alpha beta gamma ...; line 3 = numAtoms; line 4 = name;
    // then one atom per line: index element x y z <8 connectivity ints> charge.
    std::string line;
    std::getline(stream, line);
    double a, b, c;
    { std::istringstream iss(line); iss >> a >> b >> c; }
    std::getline(stream, line);
    double alpha, beta, gamma;
    { std::istringstream iss(line); iss >> alpha >> beta >> gamma; }
    std::getline(stream, line);
    std::size_t numberOfAtoms = 0;
    { std::istringstream iss(line); iss >> numberOfAtoms; }
    std::getline(stream, line);  // structure name

    SimulationBox::Type type = (std::abs(alpha - 90.0) > 1.0e-3 || std::abs(beta - 90.0) > 1.0e-3 ||
                                std::abs(gamma - 90.0) > 1.0e-3)
                                   ? SimulationBox::Type::Triclinic
                                   : SimulationBox::Type::Rectangular;
    SimulationBox simulationBox(a, b, c, alpha * std::numbers::pi / 180.0, beta * std::numbers::pi / 180.0,
                                gamma * std::numbers::pi / 180.0, type);

    std::vector<double3> fractionalPositions;
    std::vector<double> radii;
    for (std::size_t i = 0; i < numberOfAtoms; ++i)
    {
      std::getline(stream, line);
      std::istringstream iss(line);
      std::size_t index;
      std::string element;
      double x, y, z;
      iss >> index >> element >> x >> y >> z;
      std::optional<std::size_t> pseudoType = forceField.findPseudoAtom(element);
      ASSERT_TRUE(pseudoType.has_value()) << "unknown element " << element << " in " << name;
      fractionalPositions.push_back(double3(x, y, z));
      radii.push_back(0.5 * forceField(pseudoType.value(), pseudoType.value()).sizeParameter());
    }
    VoronoiNetwork network = VoronoiNetwork::create(simulationBox, fractionalPositions, radii);
    PoreDiameters diameters = PoreDiameters::compute(network);
    ChannelAnalysis channels = ChannelAnalysis::compute(network, probeRadius);

    std::string dims;
    for (const VoronoiPore& pore : channels.pores)
      if (pore.isChannel) dims += std::format("{} ", pore.dimensionality);

    PoreAccessibility accessibility =
        PoreAccessibility::create(simulationBox, fractionalPositions, radii, probeRadius);

    RandomNumber random{std::optional<std::size_t>(42)};

    // Accessible / inaccessible surface area (same scheme as VoronoiSurfaceArea::run).
    double accessibleArea = 0.0;
    double inaccessibleArea = 0.0;
    std::size_t surfaceKept = 0;
    std::size_t surfaceResampled = 0;
    std::size_t surfaceInside = 0;
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
        if (classification.resample) ++surfaceResampled;
        if (classification.inside) ++surfaceInside;
        if (!classification.resample && !classification.inside) ++surfaceKept;
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

    // Accessible / inaccessible volume (same scheme as VoronoiAccessibleVolume::run).
    std::size_t accessibleVolumeCount = 0;
    std::size_t inaccessibleVolumeCount = 0;
    std::size_t volumeResampled = 0;
    std::size_t volumeInside = 0;
    for (std::size_t s = 0; s < volumeSamples; ++s)
    {
      double3 point = simulationBox.cell * double3(random.uniform(), random.uniform(), random.uniform());
      PointClassification classification = accessibility.classify(point);
      if (classification.resample) ++volumeResampled;
      if (classification.inside) ++volumeInside;
      if (classification.inside || classification.resample) continue;
      if (classification.accessible)
        ++accessibleVolumeCount;
      else
        ++inaccessibleVolumeCount;
    }
    double volume = simulationBox.volume;
    double accessibleVolume = static_cast<double>(accessibleVolumeCount) / static_cast<double>(volumeSamples) * volume;
    double inaccessibleVolume =
        static_cast<double>(inaccessibleVolumeCount) / static_cast<double>(volumeSamples) * volume;

    // Blocking spheres: reproduce VoronoiBlockingSpheres::run. Sample the void, split into
    // accessible points and inaccessible points grouped by pocket, then greedily cover each
    // pocket with spheres. Report the number of distinct pockets that received samples and
    // the number of blocking spheres, to compare with zeo++ `-block`.
    auto periodicDistance = [&](const double3& p, const double3& q)
    { return simulationBox.applyPeriodicBoundaryConditions(p - q).length(); };
    auto mostDenseIndex = [&](const std::vector<double3>& pts) -> std::size_t
    {
      std::size_t count = std::min<std::size_t>(pts.size(), 1000);
      if (count <= 1) return 0;
      double meanDistance = 0.0;
      std::size_t pairs = 0;
      for (std::size_t p = 0; p < count; ++p)
        for (std::size_t q = p + 1; q < count; ++q)
        {
          meanDistance += periodicDistance(pts[p], pts[q]);
          ++pairs;
        }
      meanDistance /= static_cast<double>(std::max<std::size_t>(1, pairs));
      double sigmaSquared = std::max(1.0e-6, meanDistance * meanDistance);
      std::size_t best = 0;
      double bestDensity = -1.0;
      for (std::size_t p = 0; p < count; ++p)
      {
        double density = 0.0;
        for (std::size_t q = 0; q < count; ++q)
          if (p != q)
          {
            double d = periodicDistance(pts[p], pts[q]);
            density += std::exp(-d * d / sigmaSquared);
          }
        if (density > bestDensity) { bestDensity = density; best = p; }
      }
      return best;
    };

    const double overshoot = 0.1;
    std::size_t blockSamples = std::min<std::size_t>(1200000, 2000 * numberOfAtoms);  // ~zeo++ -block density
    std::vector<double3> accessiblePoints;
    std::map<std::int32_t, std::vector<double3>> pocketPoints;
    for (std::size_t s = 0; s < blockSamples; ++s)
    {
      double3 point = simulationBox.cell * double3(random.uniform(), random.uniform(), random.uniform());
      PointClassification classification = accessibility.classify(point);
      if (classification.inside || classification.resample) continue;
      if (classification.accessible)
        accessiblePoints.push_back(point);
      else if (classification.poreId >= 0)
        pocketPoints[classification.poreId].push_back(point);
    }
    std::size_t pocketsSampled = pocketPoints.size();
    std::size_t blockingSpheres = 0;
    for (auto& [poreId, pts] : pocketPoints)
    {
      std::vector<double3> remaining = pts;
      while (!remaining.empty())
      {
        double3 center = remaining[mostDenseIndex(remaining)];
        double furthestPocket = 0.0;
        for (const double3& p : remaining) furthestPocket = std::max(furthestPocket, periodicDistance(center, p));
        double closestChannel = std::numeric_limits<double>::max();
        for (const double3& p : accessiblePoints)
          closestChannel = std::min(closestChannel, periodicDistance(center, p));
        double radius;
        if (accessiblePoints.empty())
          radius = furthestPocket + probeRadius + overshoot;
        else if (furthestPocket < closestChannel)
          radius = std::min(furthestPocket + probeRadius + overshoot, closestChannel - (probeRadius + overshoot));
        else
          radius = std::max(overshoot, closestChannel - (probeRadius + overshoot));
        ++blockingSpheres;
        std::vector<double3> survivors;
        for (const double3& p : remaining)
          if (periodicDistance(center, p) >= radius) survivors.push_back(p);
        if (survivors.size() == remaining.size()) break;
        remaining = std::move(survivors);
      }
    }

    std::cout << std::format(
        "{:<7} atoms={:4d} nodes={:5d}  Di={:8.4f} Df={:8.4f} Dif={:8.4f}  chan={} dims=[{}] "
        "ASA={:9.2f} NASA={:8.2f}  AV={:9.2f} NAV={:8.2f}  pockets(net)={} pockets(sampled)={} blockSpheres={}\n",
        name, numberOfAtoms, network.nodes.size(), diameters.includedSphereDiameter,
        diameters.freeSphereDiameter, diameters.includedAlongFreePathDiameter, channels.numberOfChannels, dims,
        accessibleArea, inaccessibleArea, accessibleVolume, inaccessibleVolume, channels.numberOfPockets,
        pocketsSampled, blockingSpheres);
    std::cout << std::format("        surface samples: kept={} resampled={} inside={}  -> {:.2f}% of non-overlapping\n",
                             surfaceKept, surfaceResampled, surfaceInside,
                             100.0 * static_cast<double>(surfaceResampled + surfaceInside) /
                                 static_cast<double>(std::max<std::size_t>(1, surfaceKept + surfaceResampled + surfaceInside)));
    std::cout << std::format("        volume  samples: void={} resampled={} inside={}  -> void fraction {:.4f}\n",
                             accessibleVolumeCount + inaccessibleVolumeCount, volumeResampled, volumeInside,
                             static_cast<double>(accessibleVolumeCount + inaccessibleVolumeCount) /
                                 static_cast<double>(volumeSamples));

    // A node radius is the smallest clearance to the atoms whose cells meet at that node. In a
    // radical diagram an atom that does not touch the node can still be nearer in clearance
    // terms, which would make the stored radius an overestimate. Compare against a brute-force
    // minimum over every atom and its 27 nearest images to size that effect.
    double worstOverestimate = 0.0;
    double correctedMaximumRadius = 0.0;
    for (const VoronoiNode& node : network.nodes)
    {
      double trueClearance = std::numeric_limits<double>::max();
      for (std::size_t j = 0; j < fractionalPositions.size(); ++j)
      {
        for (int ix = -1; ix <= 1; ++ix)
        {
          for (int iy = -1; iy <= 1; ++iy)
          {
            for (int iz = -1; iz <= 1; ++iz)
            {
              double3 image = simulationBox.cell *
                              (fractionalPositions[j] + double3(static_cast<double>(ix), static_cast<double>(iy),
                                                                static_cast<double>(iz)));
              trueClearance = std::min(trueClearance, (node.position - image).length() - radii[j]);
            }
          }
        }
      }
      worstOverestimate = std::max(worstOverestimate, node.radius - trueClearance);
      correctedMaximumRadius = std::max(correctedMaximumRadius, trueClearance);
    }
    std::cout << std::format("        node clearance : Di={:.5f} Di_corrected={:.5f}  worst overestimate={:.3e} A\n",
                             diameters.includedSphereDiameter, 2.0 * correctedMaximumRadius, worstOverestimate);
  }
}
