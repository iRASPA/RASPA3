module;

module voronoi_surface_area;

import std;

import double3;
import randomnumbers;
import skspacegroupdatabase;
import crystal;
import pair_interactions;
import units;
import voronoi_network;
import pore_accessibility;
import exact_surface_patches;
import exact_boundary_components;
import exact_solvent_excluded;

SurfaceAreaSample sampleAccessibleSurfaceArea(const PoreAccessibility& accessibility, std::size_t density)
{
  RandomNumber random{samplingSeed};
  SurfaceAreaSample sample;

  for (std::size_t i = 0; i < accessibility.atomPositions.size(); ++i)
  {
    double inflatedRadius = accessibility.atomRadii[i];
    double sphereArea = 4.0 * std::numbers::pi * inflatedRadius * inflatedRadius;
    std::size_t numberOfSamples =
        std::max<std::size_t>(1, static_cast<std::size_t>(sphereArea * static_cast<double>(density)));

    std::size_t accessibleCount = 0;
    std::size_t inaccessibleCount = 0;
    for (std::size_t s = 0; s < numberOfSamples; ++s)
    {
      PointClassification classification;
      bool buried = false;
      for (std::size_t attempt = 0; attempt < resampleLimit; ++attempt)
      {
        double3 point = accessibility.atomPositions[i] + inflatedRadius * random.randomVectorOnUnitSphere();

        // A point under another inflated atom carries no area, and that is an answer rather than a
        // failure to decide, so it is kept as one instead of being drawn again.
        buried = accessibility.overlapsAtom(point, i);
        if (buried) break;

        classification = accessibility.classify(point);
        if (!classification.resample) break;
      }
      if (buried || classification.resample || classification.inside) continue;
      if (classification.accessible)
        ++accessibleCount;
      else
        ++inaccessibleCount;
    }

    double perSample = sphereArea / static_cast<double>(numberOfSamples);
    sample.accessible += static_cast<double>(accessibleCount) * perSample;
    sample.inaccessible += static_cast<double>(inaccessibleCount) * perSample;
  }

  return sample;
}

void writeExcludedSurfaceAreas(std::ostream& stream, const SolventExcludedGeometry& geometry)
{
  std::print(stream, "#\n");
  std::print(stream,
             "# The wall itself, at the same probe. The area above is the sheet the probe's centre traces over the\n"
             "# framework; the surface the probe touches is the excluded one, and it is made of three kinds of piece.\n"
             "# Where the probe rests against one atom the surface is that atom seen bare, convex. Where it rolls\n"
             "# along the crease between two it is a piece of a torus, saddle-shaped. Where it is wedged against\n"
             "# three or more it stands still and the surface is a piece of the probe's own sphere, concave. The\n"
             "# three account for the wall exactly, so the fractions below add to one along each row.\n");
  std::print(stream,
             "#\n"
             "# What the split says is how much of the wall is atom the probe can reach and how much is\n"
             "# corrugation it has to bridge. It is also what the pore-size distribution is carried by: only the\n"
             "# reentrant part moves as the probe grows, so a wall with none has no distribution at all.\n");
  std::print(stream, "#\n");
  std::print(stream, "# The four columns after the area are fractions of it, and the first three of them add to one.\n");
  std::print(stream, "#                              area [Å²]     convex     saddle    concave  reentrant\n");

  auto row = [&](const char* name, const ExcludedSurfaceAreas& areas)
  {
    if (areas.total() <= 0.0) return;
    std::print(stream, "{:<26} {:13.5f} {:10.6f} {:10.6f} {:10.6f} {:10.6f}\n", name, areas.total(),
               areas.convexFraction(), areas.saddleFraction(), areas.concaveFraction(), areas.reentrantFraction());
  };

  row("Excluded surface, all:", geometry.area);
  row("  facing channels:", geometry.accessibleArea);
  row("  facing pockets:", geometry.inaccessibleArea);
  row("  undecided:", geometry.undecidedArea);

  const SolventExcludedGeometry::Diagnostics& counts = geometry.diagnostics;
  std::print(stream, "# Pieces: {} creases, {} of them folded back on themselves at a cusp, and {} wedges, {} of\n",
             counts.numberOfArcs, counts.cuspedArcs, counts.numberOfVertices, counts.clippedVertices);
  std::print(stream, "# them cut into by a neighbouring probe. Excluded volume: {} [Å³]\n", geometry.excludedVolume);
}


void VoronoiSurfaceArea::run(const PairInteractions& interactions, const Crystal& framework, std::string probePseudoAtom,
                             Method method, std::optional<std::size_t> samplesPerAtom,
                             std::optional<std::size_t> subdivisions)
{
  std::chrono::steady_clock::time_point time_begin = std::chrono::steady_clock::now();

  std::optional<std::size_t> probeType = interactions.findType(probePseudoAtom);
  if (!probeType.has_value())
  {
    throw std::runtime_error("VoronoiSurfaceArea: Unknown probe-atom type\n");
  }
  double probeRadius = 0.5 * interactions[probeType.value()].sizeParameter;

  std::vector<double3> fractionalPositions;
  std::vector<double> radii;
  for (const CrystalAtom& atom : framework.atoms)
  {
    fractionalPositions.push_back(framework.unitCell.inverseCell * atom.position);
    std::size_t type = atom.type;
    radii.push_back(0.5 * interactions(type, type).sizeParameter);
  }

  PoreAccessibility accessibility =
      PoreAccessibility::create(framework.unitCell, fractionalPositions, radii, probeRadius);

  const std::size_t density = samplesPerAtom.value_or(50);  // per Å² (zeo++ default)
  const std::size_t panels = std::max<std::size_t>(1, subdivisions.value_or(1));
  MeasuredPatches measured;

  if (method == Method::Exact)
  {
    // The boundary is decomposed once and used twice: for the accessible area, surface by surface, and for the
    // excluded surface behind it, whose three kinds of patch hang off the same patches, creases and wedges.
    BoundaryComponents components = boundaryComponents(accessibility);
    std::vector<ComponentVerdict> verdicts = boundaryComponentVerdicts(accessibility, components);

    measured = exactAccessibleSurfaceAreaByComponent(accessibility, components, verdicts, panels);
    accessibleSurfaceArea = measured.accessible;
    inaccessibleSurfaceArea = measured.inaccessible;
    undecidedSurfaceArea = measured.undecided;

    excludedSurface = solventExcludedGeometry(accessibility, probeRadius, components, verdicts, measured, panels);
  }
  else
  {
    SurfaceAreaSample sample = sampleAccessibleSurfaceArea(accessibility, density);
    accessibleSurfaceArea = sample.accessible;
    inaccessibleSurfaceArea = sample.inaccessible;
    undecidedSurfaceArea = 0.0;
  }

  std::chrono::duration<double> timing = std::chrono::steady_clock::now() - time_begin;

  double volume = framework.unitCell.volume;
  double toGravimetric = Units::Angstrom * Units::Angstrom * Units::AvogadroConstant / framework.mass;

  std::ofstream myfile;
  myfile.open(framework.name + ".voronoi.sa.txt");
  std::print(myfile, "# Accessible / inaccessible surface area (Voronoi, {})\n",
             method == Method::Exact ? "exact" : "Monte Carlo");
  std::print(myfile, "# Crystal: {}\n", framework.name);
  std::print(myfile, "# Probe atom: {} radius: {} [Å]\n", probePseudoAtom, probeRadius);
  if (method == Method::Exact)
  {
    std::print(myfile, "# Quadrature: {} panel(s) per smooth piece\n", panels);
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
  std::print(myfile, "# CPU Timing: {} [s]\n", timing.count());
  std::print(myfile, "Accessible surface area:   {} [Å²]  {} [m²/cm³]  {} [m²/g]\n", accessibleSurfaceArea,
             1.0e4 * accessibleSurfaceArea / volume, accessibleSurfaceArea * toGravimetric);
  std::print(myfile, "Inaccessible surface area: {} [Å²]  {} [m²/cm³]  {} [m²/g]\n", inaccessibleSurfaceArea,
             1.0e4 * inaccessibleSurfaceArea / volume, inaccessibleSurfaceArea * toGravimetric);
  if (method == Method::Exact)
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
