module;

module apollonius_accessible_volume;

import std;

import double3;
import atom;
import framework;
import forcefield;
import units;
import apollonius_accessibility;
import voronoi_accessible_volume;

void ApolloniusAccessibleVolume::run(const ForceField& forceField, const Framework& framework,
                                     std::string probePseudoAtom, std::optional<std::size_t> numberOfSamples)
{
  std::chrono::steady_clock::time_point time_begin = std::chrono::steady_clock::now();

  std::optional<std::size_t> probeType = forceField.findPseudoAtom(probePseudoAtom);
  if (!probeType.has_value())
  {
    throw std::runtime_error("ApolloniusAccessibleVolume: Unknown probe-atom type\n");
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

  double volume = framework.simulationBox.volume;
  std::size_t samples = numberOfSamples.value_or(static_cast<std::size_t>(200.0 * volume));  // 200 per Å³
  samples = std::max<std::size_t>(1, samples);

  VolumeSample sample = sampleAccessibleVolume(classifier.accessibility, samples);
  accessibleVolumeFraction = sample.accessibleFraction;
  inaccessibleVolumeFraction = sample.inaccessibleFraction;
  accessibleVolume = accessibleVolumeFraction * volume;
  inaccessibleVolume = inaccessibleVolumeFraction * volume;

  std::chrono::duration<double> timing = std::chrono::steady_clock::now() - time_begin;

  double densityCrystal = framework.unitCellMass / (volume * Units::Angstrom * Units::Angstrom * Units::Angstrom *
                                                    1.0e6 * Units::AvogadroConstant);
  double toGravimetric = (Units::Angstrom * Units::Angstrom * Units::Angstrom * 1.0e6) * Units::AvogadroConstant /
                         framework.unitCellMass;

  std::ofstream myfile;
  myfile.open(framework.name + ".apollonius.av.txt");
  std::print(myfile, "# Accessible / inaccessible volume (Apollonius + Monte Carlo)\n");
  std::print(myfile, "# Framework: {}\n", framework.name);
  std::print(myfile, "# Probe atom: {} radius: {} [Å]\n", probePseudoAtom, probeRadius);
  std::print(myfile, "# Number of samples: {}\n", samples);
  std::print(myfile, "# Framework volume: {} [Å³]\n", volume);
  std::print(myfile, "# Framework density: {} [g/cm³]\n", densityCrystal);
  classifier.diagram.writeHeader(myfile);
  std::print(myfile, "# CPU Timing: {} [s]\n", timing.count());
  std::print(myfile, "Accessible volume:   fraction {}  {} [Å³]  {} [cm³/g]\n", accessibleVolumeFraction,
             accessibleVolume, accessibleVolume * toGravimetric);
  std::print(myfile, "Inaccessible volume: fraction {}  {} [Å³]  {} [cm³/g]\n", inaccessibleVolumeFraction,
             inaccessibleVolume, inaccessibleVolume * toGravimetric);
  myfile.close();
}
