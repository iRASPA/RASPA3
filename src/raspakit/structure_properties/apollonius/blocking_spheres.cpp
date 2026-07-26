module;

module apollonius_blocking_spheres;

import std;

import double3;
import atom;
import framework;
import forcefield;
import apollonius_accessibility;
import voronoi_blocking_spheres;

void ApolloniusBlockingSpheres::run(const ForceField& forceField, const Framework& framework,
                                    std::string probePseudoAtom, std::optional<std::size_t> numberOfSamples)
{
  std::optional<std::size_t> probeType = forceField.findPseudoAtom(probePseudoAtom);
  if (!probeType.has_value())
  {
    throw std::runtime_error("ApolloniusBlockingSpheres: Unknown probe-atom type\n");
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
  std::size_t samples = numberOfSamples.value_or(static_cast<std::size_t>(200.0 * volume));

  spheres = computeBlockingSpheres(classifier.accessibility, samples);

  std::ofstream myfile;
  myfile.open(framework.name + ".block");
  std::print(myfile, "{}\n", spheres.size());
  for (const BlockingSphere& sphere : spheres)
  {
    std::print(myfile, "{} {} {} {}\n", sphere.centerFractional.x, sphere.centerFractional.y,
               sphere.centerFractional.z, sphere.radius);
  }
  myfile.close();
}
