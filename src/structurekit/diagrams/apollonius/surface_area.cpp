module;

module apollonius_surface_area;

import std;

import double3;
import crystal;
import pair_interactions;
import units;
import apollonius_accessibility;
import pore_accessibility;
import exact_surface_patches;
import exact_boundary_components;
import exact_solvent_excluded;
import voronoi_surface_area;

void ApolloniusSurfaceArea::run(const PairInteractions& interactions, const Crystal& framework,
                                std::string probePseudoAtom, Method method,
                                std::optional<std::size_t> samplesPerAtom,
                                std::optional<std::size_t> subdivisions, std::string reachabilityProbe)
{
  std::chrono::steady_clock::time_point time_begin = std::chrono::steady_clock::now();

  std::optional<std::size_t> probeType = interactions.findType(probePseudoAtom);
  if (!probeType.has_value())
  {
    throw std::runtime_error("ApolloniusSurfaceArea: Unknown probe-atom type\n");
  }
  double probeRadius = 0.5 * interactions[probeType.value()].sizeParameter;

  const std::string reachName = reachabilityProbe.empty() ? probePseudoAtom : reachabilityProbe;
  const bool hybrid = reachName != probePseudoAtom;
  double reachRadius = probeRadius;
  if (hybrid)
  {
    std::optional<std::size_t> reachType = interactions.findType(reachName);
    if (!reachType.has_value())
    {
      throw std::runtime_error("ApolloniusSurfaceArea: Unknown reachability probe-atom type\n");
    }
    reachRadius = 0.5 * interactions[reachType.value()].sizeParameter;
    if (reachRadius > probeRadius + 1.0e-12)
    {
      throw std::runtime_error(
          "ApolloniusSurfaceArea: reachability probe must be no larger than the surface probe\n");
    }
  }

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

  ApolloniusAccessibility classifier =
      ApolloniusAccessibility::create(framework.unitCell, fractionalPositions, radii, probeRadius);
  std::optional<ApolloniusAccessibility> reachClassifier;
  if (hybrid)
  {
    reachClassifier = ApolloniusAccessibility::create(framework.unitCell, fractionalPositions, radii, reachRadius);
  }
  const PoreAccessibility* reachPtr = reachClassifier ? &reachClassifier->accessibility : nullptr;

  const std::size_t density = samplesPerAtom.value_or(50);  // per Å² (zeo++ default)
  const std::size_t panels = std::max<std::size_t>(1, subdivisions.value_or(1));
  MeasuredPatches measured;

  if (method == Method::Exact)
  {
    // Decomposed once, used twice: for the accessible area surface by surface, and for the excluded surface
    // behind it, whose convex, saddle and concave pieces hang off the same patches, creases and wedges.
    BoundaryComponents components = boundaryComponents(classifier.accessibility);
    std::vector<ComponentVerdict> verdicts =
        boundaryComponentVerdicts(classifier.accessibility, components, reachPtr);
    const SurfaceSidePolicy policy =
        hybrid ? SurfaceSidePolicy::network : SurfaceSidePolicy::geometric;

    measured = exactAccessibleSurfaceAreaByComponent(classifier.accessibility, components, verdicts, panels,
                                                     SurfaceMoments::volume, policy);
    accessibleSurfaceArea = measured.accessible;
    inaccessibleSurfaceArea = measured.inaccessible;
    undecidedSurfaceArea = measured.undecided;

    excludedSurface =
        solventExcludedGeometry(classifier.accessibility, probeRadius, components, verdicts, measured, panels);
  }
  else
  {
    SurfaceAreaSample sample = sampleAccessibleSurfaceArea(classifier.accessibility, density, reachPtr);
    accessibleSurfaceArea = sample.accessible;
    inaccessibleSurfaceArea = sample.inaccessible;
    undecidedSurfaceArea = 0.0;
  }

  std::chrono::duration<double> timing = std::chrono::steady_clock::now() - time_begin;

  double volume = framework.unitCell.volume;
  double toGravimetric = Units::Angstrom * Units::Angstrom * Units::AvogadroConstant / framework.mass;

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
  std::print(myfile, "# Crystal: {}\n", framework.name);
  std::print(myfile, "# Probe atom (surface geometry): {} radius: {} [Å]\n", probePseudoAtom, probeRadius);
  if (hybrid)
  {
    std::print(myfile, "# Probe that labels reachable vs sealed: {} radius: {} [Å]\n", reachName, reachRadius);
  }
  if (method == Method::Exact)
  {
    std::print(myfile, "# Quadrature: {}-point Gauss-Legendre per half panel, {} panel(s) per smooth piece\n",
               exactQuadratureOrder, panels);
    std::print(myfile, "# Surface patches measured: {} arcs\n", measured.diagnostics.numberOfArcs);
    std::print(myfile,
               "# Connected surfaces: {}, of which {} run away through the crystal, {} seal off void and {} "
               "are clusters of atoms the network was asked about\n",
               measured.numberOfSurfaces, measured.runawaySurfaces, measured.sealedSurfaces,
               measured.clusterSurfaces);

    // A surface is periodic under a subgroup of the lattice, and how many directions that subgroup spans is
    // how many directions the pore behind the surface runs away in. Integer arithmetic on the translations the
    // decomposition has already accumulated, so it is decided rather than resolved, and no pore network is
    // consulted for it.
    std::print(myfile, "# Pore system: {}{}-dimensional, from the periodicity of the surfaces themselves\n",
               measured.clusterSurfaces > 0 ? "at least " : "", measured.dimensionality());
    for (std::size_t rank = 1; rank < 4; ++rank)
    {
      if (measured.surfacesOfDimension[rank] == 0) continue;
      std::print(myfile, "#   {} surface(s) running away in {} direction(s), {} Å² of wall\n",
                 measured.surfacesOfDimension[rank], rank, measured.areaOfDimension[rank]);
    }
    if (measured.surfacesOfDimension[0] > 0)
    {
      std::print(myfile, "#   {} bounded surface(s), {} Å² of wall\n", measured.surfacesOfDimension[0],
                 measured.areaOfDimension[0]);
    }
  }
  else
  {
    std::print(myfile, "# Sample density: {} [points/Å²]\n", density);
  }
  std::print(myfile, "# Crystal volume: {} [Å³]\n", volume);
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
  if (method == Method::Exact)
  {
    writeExcludedSurfaceAreas(myfile, excludedSurface);
  }
  myfile.close();
}
