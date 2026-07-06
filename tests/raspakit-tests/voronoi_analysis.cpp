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
import voronoi_accessibility;

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
  VoronoiAccessibility accessibility = VoronoiAccessibility::create(box, {double3(0.0, 0.0, 0.0)}, {r}, probe);

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
  VoronoiAccessibility accessibility =
      VoronoiAccessibility::create(framework.simulationBox, fractionalPositions, radii, probeRadius);
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
    VoronoiAccessibility accessibility =
        VoronoiAccessibility::create(simulationBox, fractionalPositions, radii, probeRadius);
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

    VoronoiAccessibility accessibility =
        VoronoiAccessibility::create(simulationBox, fractionalPositions, radii, probeRadius);

    RandomNumber random{std::optional<std::size_t>(42)};

    // Accessible / inaccessible surface area (same scheme as VoronoiSurfaceArea::run).
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

    // Accessible / inaccessible volume (same scheme as VoronoiAccessibleVolume::run).
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
        "NASA={:8.2f} NAV={:8.2f}  pockets(net)={} pockets(sampled)={} blockSpheres={}\n",
        name, numberOfAtoms, network.nodes.size(), diameters.includedSphereDiameter,
        diameters.freeSphereDiameter, diameters.includedAlongFreePathDiameter, channels.numberOfChannels, dims,
        inaccessibleArea, inaccessibleVolume, channels.numberOfPockets, pocketsSampled, blockingSpheres);
  }
}
