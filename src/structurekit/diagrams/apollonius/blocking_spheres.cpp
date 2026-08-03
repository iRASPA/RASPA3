module;

module apollonius_blocking_spheres;

import std;

import double3;
import crystal;
import pair_interactions;
import apollonius_accessibility;
import exact_void_split;
import voronoi_blocking_spheres;

void ApolloniusBlockingSpheres::run(const PairInteractions& interactions, const Crystal& framework,
                                    std::string probePseudoAtom, std::optional<std::size_t> numberOfSamples)
{
  std::optional<std::size_t> probeType = interactions.findType(probePseudoAtom);
  if (!probeType.has_value())
  {
    throw std::runtime_error("ApolloniusBlockingSpheres: Unknown probe-atom type\n");
  }
  double probeRadius = 0.5 * interactions[probeType.value()].sizeParameter;

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

  double volume = framework.unitCell.volume;

  ExactVoidSplit split = exactVoidSplitByComponents(classifier.accessibility, volume);
  fallbackReason = measuredSpheresRefused(split);
  if (fallbackReason.empty())
  {
    measured = true;
    pockets = split.pockets;
    spheres = exactBlockingSpheres(split);
  }
  else
  {
    std::size_t samples = numberOfSamples.value_or(static_cast<std::size_t>(200.0 * volume));
    spheres = computeBlockingSpheres(classifier.accessibility, samples);
  }

  writeBlockingSpheres(framework.name, "apollonius", probePseudoAtom, probeRadius, spheres, pockets, fallbackReason);
}
