module;

module apollonius_accessible_volume;

import std;

import double3;
import atom;
import framework;
import forcefield;
import units;
import apollonius_accessibility;
import exact_union_volume;
import exact_void_split;
import voronoi_accessible_volume;

void ApolloniusAccessibleVolume::run(const ForceField& forceField, const Framework& framework,
                                     std::string probePseudoAtom, Method method,
                                     std::optional<std::size_t> numberOfSamples,
                                     std::optional<std::size_t> subdivisions)
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
  const std::size_t panels = std::max<std::size_t>(1, subdivisions.value_or(1));

  // One classified sweep of the patches carries both halves of the answer: the volume of the union, hence
  // the whole of the void, and the volume of every sealed pore, from the surface that bounds it.
  bool sampled = method != Method::Exact;
  if (method == Method::Exact)
  {
    ExactVoidSplit split = exactVoidSplitByComponents(classifier.accessibility, volume, panels);
    voidVolume = split.voidVolume;
    voidFraction = voidVolume / volume;
    closureDefect = split.closureDefect;
    numberOfPockets = split.numberOfPockets;
    numberOfSurfaces = split.numberOfSurfaces;
    provedSurfaces = split.provedSurfaces;
    provedPockets = split.provedPockets;
    numberOfEnclosedSolids = split.numberOfEnclosedSolids;
    signDisagreements = split.signDisagreements;
    signDisagreementVolume = split.signDisagreementVolume;
    splitRejection = split.rejection;

    if (split.reliable)
    {
      splitMeasured = true;
      accessibleVolume = split.accessibleVolume;
      inaccessibleVolume = split.inaccessibleVolume;
      accessibleVolumeFraction = accessibleVolume / volume;
      inaccessibleVolumeFraction = inaccessibleVolume / volume;
    }
    else
    {
      // Surface that faces no pore, or a pocket whose boundary did not close: the arcs cannot be trusted
      // to divide the void, so the division falls back on the sampling it was to have replaced. The total
      // stands either way, belonging to the union of the atoms rather than to the diagram.
      sampled = true;
    }
  }

  if (sampled)
  {
    VolumeSample sample = sampleAccessibleVolume(classifier.accessibility, samples);
    accessibleVolumeFraction = sample.accessibleFraction;
    inaccessibleVolumeFraction = sample.inaccessibleFraction;
    accessibleVolume = accessibleVolumeFraction * volume;
    inaccessibleVolume = inaccessibleVolumeFraction * volume;

    if (method != Method::Exact)
    {
      voidFraction = accessibleVolumeFraction + inaccessibleVolumeFraction;
      voidVolume = voidFraction * volume;
    }
  }

  std::chrono::duration<double> timing = std::chrono::steady_clock::now() - time_begin;

  double densityCrystal = framework.unitCellMass / (volume * Units::Angstrom * Units::Angstrom * Units::Angstrom *
                                                    1.0e6 * Units::AvogadroConstant);
  double toGravimetric = (Units::Angstrom * Units::Angstrom * Units::Angstrom * 1.0e6) * Units::AvogadroConstant /
                         framework.unitCellMass;

  std::ofstream myfile;
  myfile.open(framework.name + ".apollonius.av.txt");
  std::print(myfile, "# Void volume (Apollonius, {})\n", method == Method::Exact ? "exact" : "Monte Carlo");
  std::print(myfile, "# Framework: {}\n", framework.name);
  std::print(myfile, "# Probe atom: {} radius: {} [Å]\n", probePseudoAtom, probeRadius);
  if (!splitMeasured) std::print(myfile, "# Number of samples: {}\n", samples);
  if (method == Method::Exact)
  {
    std::print(myfile, "# Quadrature behind the patch areas: {} panel(s) per smooth piece\n", panels);
    std::print(myfile, "# Pockets measured: {}\n", numberOfPockets);
    std::print(myfile, "# Connected surfaces: {}, of which {} settled by a free ball; pockets so settled: {}\n",
               numberOfSurfaces, provedSurfaces, provedPockets);
    if (numberOfEnclosedSolids > 0)
    {
      std::print(myfile, "# Clusters standing inside a pocket, whose room was subtracted: {}\n",
                 numberOfEnclosedSolids);
    }
    if (signDisagreements > 0)
    {
      std::print(myfile, "# Pockets the network takes for channels: {} ({} Å³)\n", signDisagreements,
                 signDisagreementVolume);
    }
    if (splitMeasured)
    {
      std::print(myfile, "# Closure defect over the pockets: {} (zero for a boundary that closes)\n", closureDefect);
    }
    else
    {
      std::print(myfile, "# The measured division was rejected ({}); the division below is sampled\n",
                 splitRejection);
    }
  }
  std::print(myfile, "# Framework volume: {} [Å³]\n", volume);
  std::print(myfile, "# Framework density: {} [g/cm³]\n", densityCrystal);
  classifier.diagram.writeHeader(myfile);
  std::print(myfile, "# CPU Timing: {} [s]\n", timing.count());
  std::print(myfile, "Void volume:         fraction {}  {} [Å³]  {} [cm³/g]{}\n", voidFraction, voidVolume,
             voidVolume * toGravimetric, method == Method::Exact ? "  (measured)" : "");
  std::print(myfile, "Accessible volume:   fraction {}  {} [Å³]  {} [cm³/g]{}\n", accessibleVolumeFraction,
             accessibleVolume, accessibleVolume * toGravimetric, splitMeasured ? "  (measured)" : "");
  std::print(myfile, "Inaccessible volume: fraction {}  {} [Å³]  {} [cm³/g]{}\n", inaccessibleVolumeFraction,
             inaccessibleVolume, inaccessibleVolume * toGravimetric, splitMeasured ? "  (measured)" : "");
  myfile.close();
}
