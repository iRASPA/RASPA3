module;

module voronoi_surface_area;

import std;

import double3;
import randomnumbers;
import skspacegroupdatabase;
import atom;
import framework;
import forcefield;
import units;
import voronoi_network;
import voronoi_accessibility;

void VoronoiSurfaceArea::run(const ForceField& forceField, const Framework& framework, std::string probePseudoAtom,
                             std::optional<std::size_t> samplesPerAtom)
{
  RandomNumber random{std::nullopt};
  std::chrono::system_clock::time_point time_begin = std::chrono::system_clock::now();

  std::optional<std::size_t> probeType = forceField.findPseudoAtom(probePseudoAtom);
  if (!probeType.has_value())
  {
    throw std::runtime_error("VoronoiSurfaceArea: Unknown probe-atom type\n");
  }
  double probeRadius = 0.5 * forceField[probeType.value()].sizeParameter();

  std::vector<double3> fractionalPositions;
  std::vector<double> radii;
  for (const Atom& atom : framework.unitCellAtoms)
  {
    fractionalPositions.push_back(framework.simulationBox.inverseCell * atom.position);
    std::size_t type = static_cast<std::size_t>(atom.type);
    radii.push_back(0.5 * forceField(type, type).sizeParameter());
  }

  VoronoiAccessibility accessibility =
      VoronoiAccessibility::create(framework.simulationBox, fractionalPositions, radii, probeRadius);

  const std::size_t density = samplesPerAtom.value_or(50);  // per Å² (zeo++ default)

  double accessibleArea = 0.0;
  double inaccessibleArea = 0.0;

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
      double3 point = accessibility.atomPositions[i] + inflatedRadius * random.randomVectorOnUnitSphere();

      // Reject the point if it lies inside any other inflated atom.
      if (accessibility.overlapsAtom(point, i)) continue;

      PointClassification classification = accessibility.classify(point);
      if (classification.resample || classification.inside) continue;
      if (classification.accessible)
        ++accessibleCount;
      else
        ++inaccessibleCount;
    }

    double perSample = sphereArea / static_cast<double>(numberOfSamples);
    accessibleArea += static_cast<double>(accessibleCount) * perSample;
    inaccessibleArea += static_cast<double>(inaccessibleCount) * perSample;
  }

  accessibleSurfaceArea = accessibleArea;
  inaccessibleSurfaceArea = inaccessibleArea;

  std::chrono::duration<double> timing = std::chrono::system_clock::now() - time_begin;

  double volume = framework.simulationBox.volume;
  double toGravimetric = Units::Angstrom * Units::Angstrom * Units::AvogadroConstant / framework.unitCellMass;

  std::ofstream myfile;
  myfile.open(framework.name + ".voronoi.sa.txt");
  std::print(myfile, "# Accessible / inaccessible surface area (Voronoi + Monte Carlo)\n");
  std::print(myfile, "# Framework: {}\n", framework.name);
  std::print(myfile, "# Probe atom: {} radius: {} [Å]\n", probePseudoAtom, probeRadius);
  std::print(myfile, "# Sample density: {} [points/Å²]\n", density);
  std::print(myfile, "# Framework volume: {} [Å³]\n", volume);
  std::print(myfile, "# CPU Timing: {} [s]\n", timing.count());
  std::print(myfile, "Accessible surface area:   {} [Å²]  {} [m²/cm³]  {} [m²/g]\n", accessibleSurfaceArea,
             1.0e4 * accessibleSurfaceArea / volume, accessibleSurfaceArea * toGravimetric);
  std::print(myfile, "Inaccessible surface area: {} [Å²]  {} [m²/cm³]  {} [m²/g]\n", inaccessibleSurfaceArea,
             1.0e4 * inaccessibleSurfaceArea / volume, inaccessibleSurfaceArea * toGravimetric);
  std::print(myfile, "Total surface area:        {} [Å²]\n", accessibleSurfaceArea + inaccessibleSurfaceArea);
  myfile.close();
}
