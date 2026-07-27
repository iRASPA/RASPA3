module;

module apollonius_surface_area;

import std;

import double3;
import atom;
import framework;
import forcefield;
import units;
import apollonius_accessibility;
import exact_surface_patches;
import voronoi_surface_area;

void ApolloniusSurfaceArea::run(const ForceField& forceField, const Framework& framework,
                                std::string probePseudoAtom, Method method,
                                std::optional<std::size_t> samplesPerAtom,
                                std::optional<std::size_t> subdivisions)
{
  std::chrono::steady_clock::time_point time_begin = std::chrono::steady_clock::now();

  std::optional<std::size_t> probeType = forceField.findPseudoAtom(probePseudoAtom);
  if (!probeType.has_value())
  {
    throw std::runtime_error("ApolloniusSurfaceArea: Unknown probe-atom type\n");
  }
  double probeRadius = 0.5 * forceField[probeType.value()].sizeParameter();

  std::vector<double3> fractionalPositions;
  std::vector<double> radii;
  fractionalPositions.reserve(framework.unitCellAtoms.size());
  radii.reserve(framework.unitCellAtoms.size());
  for (const Atom& atom : framework.unitCellAtoms)
  {
    fractionalPositions.push_back(framework.simulationBox.inverseCell * atom.position);
    std::size_t type = static_cast<std::size_t>(atom.type);
    radii.push_back(0.5 * forceField(type, type).sizeParameter());
  }

  ApolloniusAccessibility classifier =
      ApolloniusAccessibility::create(framework.simulationBox, fractionalPositions, radii, probeRadius);

  const std::size_t density = samplesPerAtom.value_or(50);  // per Å² (zeo++ default)
  const std::size_t panels = std::max<std::size_t>(1, subdivisions.value_or(1));
  ExactSurfaceAreaSample measured;

  if (method == Method::Exact)
  {
    measured = exactAccessibleSurfaceAreaByComponent(classifier.accessibility, panels);
    accessibleSurfaceArea = measured.accessible;
    inaccessibleSurfaceArea = measured.inaccessible;
    undecidedSurfaceArea = measured.undecided;
  }
  else
  {
    SurfaceAreaSample sample = sampleAccessibleSurfaceArea(classifier.accessibility, density);
    accessibleSurfaceArea = sample.accessible;
    inaccessibleSurfaceArea = sample.inaccessible;
    undecidedSurfaceArea = 0.0;
  }

  std::chrono::duration<double> timing = std::chrono::steady_clock::now() - time_begin;

  double volume = framework.simulationBox.volume;
  double toGravimetric = Units::Angstrom * Units::Angstrom * Units::AvogadroConstant / framework.unitCellMass;

  std::ofstream myfile;
  myfile.open(framework.name + ".apollonius.sa.txt");
  if (method == Method::Exact)
  {
    std::print(myfile, "# Accessible / inaccessible surface area (Apollonius, exact)\n");
  }
  else
  {
    std::print(myfile, "# Accessible / inaccessible surface area (Apollonius + Monte Carlo)\n");
  }
  std::print(myfile, "# Framework: {}\n", framework.name);
  std::print(myfile, "# Probe atom: {} radius: {} [Å]\n", probePseudoAtom, probeRadius);
  if (method == Method::Exact)
  {
    std::print(myfile, "# Quadrature: {}-point Gauss-Legendre per half panel, {} panel(s) per smooth piece\n",
               exactQuadratureOrder, panels);
    std::print(myfile, "# Surface patches measured: {} arcs\n", measured.numberOfArcs);
    std::print(myfile,
               "# Connected surfaces: {}, of which {} run away through the crystal, {} seal off void and {} "
               "are clusters of atoms the network was asked about\n",
               measured.numberOfSurfaces, measured.runawaySurfaces, measured.sealedSurfaces,
               measured.clusterSurfaces);
  }
  else
  {
    std::print(myfile, "# Sample density: {} [points/Å²]\n", density);
  }
  std::print(myfile, "# Framework volume: {} [Å³]\n", volume);
  classifier.diagram.writeHeader(myfile);
  std::print(myfile, "# CPU Timing: {} [s]\n", timing.count());
  std::print(myfile, "Accessible surface area:   {} [Å²]  {} [m²/cm³]  {} [m²/g]\n", accessibleSurfaceArea,
             1.0e4 * accessibleSurfaceArea / volume, accessibleSurfaceArea * toGravimetric);
  std::print(myfile, "Inaccessible surface area: {} [Å²]  {} [m²/cm³]  {} [m²/g]\n", inaccessibleSurfaceArea,
             1.0e4 * inaccessibleSurfaceArea / volume, inaccessibleSurfaceArea * toGravimetric);
  if (undecidedSurfaceArea > 0.0)
  {
    std::print(myfile, "Undecided surface area:    {} [Å²]  {} [m²/cm³]  {} [m²/g]\n", undecidedSurfaceArea,
               1.0e4 * undecidedSurfaceArea / volume, undecidedSurfaceArea * toGravimetric);
  }
  std::print(myfile, "Total surface area:        {} [Å²]\n",
             accessibleSurfaceArea + inaccessibleSurfaceArea + undecidedSurfaceArea);
  myfile.close();
}
