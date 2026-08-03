module;

module brute_force_validation;

import std;

import double3;
import crystal;
import pair_interactions;
import units;
import unit_cell;

import apollonius_accessibility;
import apollonius_network;
import pore_diameters;
import blocking_spheres;
import voronoi_blocking_spheres;
import exact_boundary_components;
import exact_surface_patches;
import exact_solvent_excluded;
import exact_void_split;

import brute_force_structure;
import brute_force_voxels;
import brute_force_diameters;
import brute_force_surface_area;
import brute_force_solvent_excluded;
import brute_force_pore_volume;
import brute_force_blocking_pockets;

namespace
{
// The atoms of one cell, as the brute force wants them: where they are and how far out of bounds they are.
BruteForceStructure structureWith(const Crystal &framework, std::vector<double> radii)
{
  return BruteForceStructure::create(framework.name, framework.unitCell,
                                     framework.cartesianPositions(), std::move(radii));
}

// The atoms' own radii, which is what the void is measured against when the probe is a point.
std::vector<double> ownRadii(const PairInteractions &interactions, const Crystal &framework)
{
  std::vector<double> radii;
  radii.reserve(framework.atoms.size());
  for (const CrystalAtom &atom : framework.atoms)
  {
    std::size_t type = atom.type;
    radii.push_back(0.5 * interactions(type, type).sizeParameter);
  }
  return radii;
}

double probeRadiusOf(const PairInteractions &interactions, const std::string &name)
{
  std::optional<std::size_t> type = interactions.findType(name);
  if (!type.has_value())
  {
    throw std::runtime_error(std::format("BruteForceValidation: unknown probe atom '{}'\n", name));
  }
  return 0.5 * interactions[type.value()].sizeParameter;
}

std::vector<double> inflatedBy(const std::vector<double> &radii, double probeRadius)
{
  std::vector<double> out;
  out.reserve(radii.size());
  for (double radius : radii) out.push_back(radius + probeRadius);
  return out;
}

void writeCheck(std::ostream &stream, const BruteForceCheck &check)
{
  if (!check.applicable)
  {
    std::print(stream, "  {:<44} {}\n", check.property, check.basis);
    return;
  }

  std::print(stream, "  {:<44} exact {:>16.6f}   brute force {:>16.6f}   difference {:>+13.6f} [{}]   {}\n",
             check.property, check.exact, check.bruteForce, check.difference(), check.units,
             check.agrees() ? "agrees" : "DISAGREES");
  std::print(stream, "  {:<44} allowed {:>14.6f}   {}\n", "", check.tolerance, check.basis);
}
}  // namespace

std::size_t BruteForceValidation::numberOfDisagreements() const
{
  return static_cast<std::size_t>(std::ranges::count_if(this->checks, [](const BruteForceCheck &check)
                                                        { return !check.agrees(); }));
}

void BruteForceValidation::run(const PairInteractions &interactions, const Crystal &framework,
                               const std::string &surfaceProbe, const std::string &voidProbe,
                               const BruteForceSettings &settings)
{
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

  const double volume = framework.unitCell.volume;
  const std::size_t panels = std::max<std::size_t>(1, settings.subdivisions);

  std::vector<double> bareRadii = ownRadii(interactions, framework);

  std::vector<double3> fractionalPositions;
  fractionalPositions.reserve(framework.atoms.size());
  for (const CrystalAtom &atom : framework.atoms)
  {
    fractionalPositions.push_back(framework.unitCell.inverseCell * atom.position);
  }

  double surfaceProbeRadius = probeRadiusOf(interactions, surfaceProbe);
  double voidProbeRadius = probeRadiusOf(interactions, voidProbe);

  // ---- the pore diameters, against the Apollonius network -------------------------------------------
  {
    ApolloniusPoreNetwork pores =
        ApolloniusPoreNetwork::create(framework.unitCell, fractionalPositions, bareRadii);
    PoreDiameters exact = PoreDiameters::compute(pores.network);

    BruteForceStructure bare = structureWith(framework, bareRadii);
    BruteForceVoxels voxels = BruteForceVoxels::build(bare, settings.spacing);
    this->diameters = BruteForceDiameters::compute(bare, voxels);

    double grain = std::max({voxels.spacing.x, voxels.spacing.y, voxels.spacing.z});

    // Di is walked to a local maximum from the roomiest points of the grid, so it is not limited by the
    // spacing and the two should stand on top of each other.
    this->checks.push_back(BruteForceCheck{
        .property = "Di, largest included sphere",
        .units = "Å",
        .exact = exact.includedSphereDiameter,
        .bruteForce = this->diameters.includedSphereDiameter,
        .tolerance = 1.0e-4,
        .basis = "the walk uphill converges to the same peak, so only rounding should separate them"});

    // Df and Dif come off grid lines that pass near the bottleneck rather than through it, so they come out
    // small and creep up as the grid is refined. Two grid steps is what that is worth.
    this->checks.push_back(
        BruteForceCheck{.property = "Df, largest free sphere",
                        .units = "Å",
                        .exact = exact.freeSphereDiameter,
                        .bruteForce = this->diameters.freeSphereDiameter,
                        .tolerance = 4.0 * grain,
                        .basis = std::format("two grid steps either way, the grid being {:.3f} Å", grain),
                        .applicable = this->diameters.percolates || exact.freeSphereDiameter > 0.0});

    // The one part of it that is not an estimate: every hop was credited only with the width a sphere is
    // shown to manage along the whole of it, so a sphere of this size demonstrably gets through, and an
    // exact Df below it would be an exact Df that is wrong.
    this->checks.push_back(BruteForceCheck{
        .property = "Df is at least what was demonstrated",
        .units = "Å",
        .exact = exact.freeSphereDiameter,
        .bruteForce = this->diameters.guaranteedFreeSphereDiameter,
        .tolerance = std::max(exact.freeSphereDiameter - this->diameters.guaranteedFreeSphereDiameter, 0.0),
        .basis = "a lower bound with nothing estimated in it: the exact Df may exceed it and not fall short",
        .applicable = this->diameters.percolates});

    this->checks.push_back(
        BruteForceCheck{.property = "Dif, included sphere along the free path",
                        .units = "Å",
                        .exact = exact.includedAlongFreePathDiameter,
                        .bruteForce = this->diameters.includedAlongFreePathDiameter,
                        .tolerance = 4.0 * grain,
                        .basis = std::format("two grid steps either way, the grid being {:.3f} Å", grain),
                        .applicable = this->diameters.percolates || exact.freeSphereDiameter > 0.0});
  }

  // ---- the surface area and its decomposition, against the exact sweep ------------------------------
  {
    ApolloniusAccessibility classifier = ApolloniusAccessibility::create(
        framework.unitCell, fractionalPositions, bareRadii, surfaceProbeRadius);

    BoundaryComponents components = boundaryComponents(classifier.accessibility);
    std::vector<ComponentVerdict> verdicts = boundaryComponentVerdicts(classifier.accessibility, components);
    MeasuredPatches measured =
        exactAccessibleSurfaceAreaByComponent(classifier.accessibility, components, verdicts, panels);

    BruteForceStructure bare = structureWith(framework, bareRadii);
    BruteForceStructure inflated = structureWith(framework, inflatedBy(bareRadii, surfaceProbeRadius));
    BruteForceVoxels voxels = BruteForceVoxels::build(inflated, settings.spacing);

    this->surfaceNecksProved = voxels.necksProved;
    this->surfaceNecksTried = voxels.necksTried;

    this->surfaceArea = BruteForceSurfaceArea::compute(inflated, voxels, settings.samplesPerAtom,
                                                       !settings.skipSolventExcluded);

    // Four standard errors is a one-in-sixteen-thousand coincidence, which is the point at which a
    // disagreement is worth looking at rather than worth repeating with more points.
    double areaTolerance = 4.0 * this->surfaceArea.totalAreaError;

    this->checks.push_back(
        BruteForceCheck{.property = "Surface area, total",
                        .units = "Å²",
                        .exact = measured.total(),
                        .bruteForce = this->surfaceArea.totalArea,
                        .tolerance = areaTolerance,
                        .basis = std::format("four standard errors of {} directions per atom",
                                             settings.samplesPerAtom)});

    // Splitting it needs the flood as well as the sample. Where stepping off the wall found no void the
    // grid had labelled, the area behind that point went to neither side, so what it could have been worth
    // is exactly the allowance the split needs.
    double splitSlack = this->surfaceArea.undecidedArea + 0.005 * measured.total();

    this->checks.push_back(BruteForceCheck{
        .property = "Surface area, reachable from outside",
        .units = "Å²",
        .exact = measured.accessible,
        .bruteForce = this->surfaceArea.accessibleArea,
        .tolerance = areaTolerance + splitSlack,
        .basis = std::format("the same, and {:.4f} Å² of wall that stepping off it found no labelled void "
                             "behind, so that it went to neither side",
                             this->surfaceArea.undecidedArea)});

    this->checks.push_back(
        BruteForceCheck{.property = "Surface area, sealed in pockets",
                        .units = "Å²",
                        .exact = measured.inaccessible,
                        .bruteForce = this->surfaceArea.inaccessibleArea,
                        .tolerance = areaTolerance + splitSlack,
                        .basis = "the same"});

    if (!settings.skipSolventExcluded)
    {
      SolventExcludedGeometry excluded = solventExcludedGeometry(classifier.accessibility, surfaceProbeRadius,
                                                                 components, verdicts, measured, panels);

      this->exactVertices = excluded.diagnostics.numberOfVertices;
      this->exactClippedVertices = excluded.diagnostics.clippedVertices;
      this->exactDegenerateVertices = excluded.diagnostics.degenerateVertices;
      this->exactVanishedVertices = excluded.diagnostics.vanishedVertices;
      this->exactDiscardedCorners = excluded.diagnostics.discardedCorners;

      this->solventExcluded =
          BruteForceSolventExcluded::compute(bare, inflated, this->surfaceArea, surfaceProbeRadius,
                                             settings.samplesPerAtom, settings.creaseSteps,
                                             settings.cornerSamples);

      this->checks.push_back(
          BruteForceCheck{.property = "Excluded surface, convex",
                          .units = "Å²",
                          .exact = excluded.area.convex,
                          .bruteForce = this->solventExcluded.convexArea,
                          .tolerance = 4.0 * this->solventExcluded.convexAreaError,
                          .basis = "four standard errors; the same directions as the surface area, weighted "
                                   "by the atom's own radius instead of the inflated one"});

      // The saddle sweep is a quadrature rather than a sample, so what limits it is the step in the two
      // angles rather than a count of points, and its error falls off as the square of the step.
      this->checks.push_back(
          BruteForceCheck{.property = "Excluded surface, saddle",
                          .units = "Å²",
                          .exact = excluded.area.saddle,
                          .bruteForce = this->solventExcluded.saddleArea,
                          .tolerance = 0.01 * std::max(excluded.area.saddle, 1.0),
                          .basis = "one per cent, the sweep round each crease being a midpoint rule in two "
                                   "angles rather than a sample"});

      this->checks.push_back(
          BruteForceCheck{.property = "Excluded surface, concave",
                          .units = "Å²",
                          .exact = excluded.area.concave,
                          .bruteForce = this->solventExcluded.concaveArea,
                          .tolerance = 4.0 * this->solventExcluded.concaveAreaError +
                                       0.01 * std::max(excluded.area.concave, 1.0),
                          .basis = "four standard errors of the draw over each corner, and a per cent for "
                                   "the corners themselves"});

      this->checks.push_back(
          BruteForceCheck{.property = "Excluded surface, total",
                          .units = "Å²",
                          .exact = excluded.area.total(),
                          .bruteForce = this->solventExcluded.totalArea(),
                          .tolerance = 4.0 * this->solventExcluded.convexAreaError +
                                       4.0 * this->solventExcluded.concaveAreaError +
                                       0.01 * std::max(excluded.area.saddle + excluded.area.concave, 1.0),
                          .basis = "the three above added, and their allowances with them"});
    }

    // ---- the blocking spheres, held against the flood ------------------------------------------------
    {
      ExactVoidSplit split = exactVoidSplitByComponents(classifier.accessibility, volume, panels);

      this->blockingPockets = BruteForceBlockingPockets::compute(inflated, voxels);

      // Reported and not judged. How many pieces the void falls into is the one thing a grid is entitled to
      // disagree about: a neck it cannot resolve splits one pocket into two, and a neck it can resolve but
      // the classifier cannot joins two into one. What the spheres have to get right is not the count but
      // the cover, and that is the check below.
      this->checks.push_back(BruteForceCheck{
          .property = "Sealed pockets, how many",
          .units = "",
          .exact = static_cast<double>(split.numberOfPockets),
          .bruteForce = static_cast<double>(this->blockingPockets.pockets.size()),
          .applicable = false,
          .basis = std::format("exact {}, brute force {} (reported, not judged: a grid may split a pocket "
                               "at a neck it cannot resolve)",
                               split.numberOfPockets, this->blockingPockets.pockets.size())});

      std::string refused = measuredSpheresRefused(split);
      if (refused.empty())
      {
        std::vector<BlockingSphere> spheres = exactBlockingSpheres(split);

        std::vector<std::pair<double3, double>> asPairs;
        asPairs.reserve(spheres.size());
        for (const BlockingSphere &sphere : spheres)
        {
          asPairs.emplace_back(sphere.centerFractional, sphere.radius);
        }

        this->blockingAudit = auditBlockingSpheres(inflated, voxels, asPairs);

        // Not a comparison of two numbers but a property the spheres either have or do not: they must cover
        // every scrap of void a molecule cannot reach and none of the void it can. A sphere that fails the
        // first lets a molecule be inserted where it could never have got to; one that fails the second
        // blocks off part of a pore that was open.
        this->checks.push_back(BruteForceCheck{
            .property = "Blocking spheres leave no sealed void open",
            .units = "Å³",
            .exact = 0.0,
            .bruteForce = this->blockingAudit.volumeLeftOpen,
            .tolerance = 4.0 * voxels.voxelVolume,
            .basis = "sealed void outside every sphere, which a simulation would insert into"});

        this->checks.push_back(BruteForceCheck{
            .property = "Blocking spheres cover no open void",
            .units = "Å³",
            .exact = 0.0,
            .bruteForce = this->blockingAudit.volumeCoveredUp,
            .tolerance = 4.0 * voxels.voxelVolume,
            .basis = "void in a channel that the spheres refuse, which a simulation would then never fill"});
      }
    }
  }

  // ---- the void fraction and its split, against the exact one ---------------------------------------
  {
    ApolloniusAccessibility classifier =
        ApolloniusAccessibility::create(framework.unitCell, fractionalPositions, bareRadii, voidProbeRadius);

    ExactVoidSplit split = exactVoidSplitByComponents(classifier.accessibility, volume, panels);

    BruteForceStructure inflated = structureWith(framework, inflatedBy(bareRadii, voidProbeRadius));
    BruteForceVoxels voxels = BruteForceVoxels::build(inflated, settings.spacing);

    this->poreVolume = BruteForcePoreVolume::compute(inflated, voxels, settings.volumePoints);

    double volumeTolerance = 4.0 * this->poreVolume.voidFractionError * volume;

    this->checks.push_back(
        BruteForceCheck{.property = "Pore volume, total",
                        .units = "Å³",
                        .exact = split.voidVolume,
                        .bruteForce = this->poreVolume.voidVolume,
                        .tolerance = volumeTolerance,
                        .basis = std::format("four standard errors of {} points thrown at the cell",
                                             this->poreVolume.numberOfPoints)});

    // The split is the flood's, not the sample's, so it carries the grid's reach as well as the sample's
    // error. Void too thin for the grid to have labelled at all went to neither side, and that volume is
    // measured rather than guessed at, so it is what the allowance is made of.
    double splitSlack = this->poreVolume.unresolvedVolume + 0.005 * std::max(split.voidVolume, 1.0);

    this->checks.push_back(
        BruteForceCheck{.property = "Pore volume, reachable from outside",
                        .units = "Å³",
                        .exact = split.accessibleVolume,
                        .bruteForce = this->poreVolume.accessibleVolume,
                        .tolerance = volumeTolerance + splitSlack,
                        .basis = std::format("the same, and {:.4f} Å³ of void too thin for the grid to "
                                             "have labelled, which went to neither side",
                                             this->poreVolume.unresolvedVolume)});

    this->checks.push_back(BruteForceCheck{.property = "Pore volume, sealed in pockets",
                                           .units = "Å³",
                                           .exact = split.inaccessibleVolume,
                                           .bruteForce = this->poreVolume.inaccessibleVolume,
                                           .tolerance = volumeTolerance + splitSlack,
                                           .basis = "the same"});
  }

  std::chrono::duration<double> timing = std::chrono::steady_clock::now() - begin;
  this->seconds = timing.count();

  // ---- the report ------------------------------------------------------------------------------------
  std::ofstream report(framework.name + ".brute-force.txt");

  std::print(report, "# The exact routes checked against brute force\n");
  std::print(report, "#\n");
  std::print(report, "# Every number on the left is what the exact geometry says. Every number on the right\n");
  std::print(report, "# was worked out again from the atom positions and radii alone, on a grid and with a\n");
  std::print(report, "# great many random points, by code that shares no diagram, no network and no pore\n");
  std::print(report, "# classifier with the left-hand side. The allowance under each pair is what the\n");
  std::print(report, "# right-hand side can resolve, not what the left-hand side is hoped to be worth.\n");
  std::print(report, "#\n");
  std::print(report, "# Crystal: {}\n", framework.name);
  std::print(report, "# Atoms in the unit cell: {}\n", framework.atoms.size());
  std::print(report, "# Cell volume: {} [Å³]\n", volume);
  std::print(report, "# Probe for the surface, the excluded surface and the blocking spheres: {}, radius {} [Å]\n",
             surfaceProbe, surfaceProbeRadius);
  std::print(report, "# Probe for the void: {}, radius {} [Å]\n", voidProbe, voidProbeRadius);
  std::print(report, "# Grid the void is flooded on: about {} Å between voxel centres\n", settings.spacing);
  std::print(report, "# Directions per atom: {}\n", settings.samplesPerAtom);
  std::print(report, "# Points thrown at the cell: {}\n", this->poreVolume.numberOfPoints);
  std::print(report, "# Quadrature of the exact route: {} panel(s) per smooth piece\n", panels);
  std::print(report, "# CPU Timing: {} [s]\n", this->seconds);
  std::print(report, "\n");

  std::size_t disagreements = this->numberOfDisagreements();
  std::print(report, "{} of {} checks agree.\n\n",
             this->checks.size() - disagreements, this->checks.size());

  for (const BruteForceCheck &check : this->checks) writeCheck(report, check);

  std::print(report, "\n# What the brute force found on its own\n");
  std::print(report, "Void voxels flooded into {} channel(s) and {} pocket(s)\n", this->poreVolume.numberOfChannels,
             this->poreVolume.numberOfPockets);
  std::print(report, "Void fraction from the points: {} ± {}\n", this->poreVolume.voidFraction,
             this->poreVolume.voidFractionError);
  std::print(report, "Void fraction from the voxels: {}\n", this->poreVolume.voidFractionFromVoxels);
  std::print(report,
             "Necks too narrow for the grid to see, proved passable by a straight line: {} of {} tried on "
             "the void grid, {} of {} on the probe's\n",
             this->poreVolume.necksProved, this->poreVolume.necksTried, this->surfaceNecksProved,
             this->surfaceNecksTried);
  std::print(report, "Void too thin for the grid to have labelled: {} Å³ of {} Å³\n",
             this->poreVolume.unresolvedVolume, this->poreVolume.voidVolume);
  std::print(report, "Wall with no labelled void behind it: {} Å² of {} Å²\n", this->surfaceArea.undecidedArea,
             this->surfaceArea.totalArea);
  std::print(report, "Widest path runs away in {} direction(s)\n", this->diameters.dimensionality);
  std::print(report, "The walk uphill improved on the best voxel by {} Å\n", this->diameters.walkGainedForDi);

  if (!settings.skipSolventExcluded)
  {
    std::print(report, "Creases swept: {}, corners found: {} (the exact route had {} vertices, {} clipped, "
                       "{} degenerate, {} vanished, {} corners discarded)\n",
               this->solventExcluded.numberOfPairs, this->solventExcluded.numberOfCorners,
               this->exactVertices, this->exactClippedVertices, this->exactDegenerateVertices,
               this->exactVanishedVertices, this->exactDiscardedCorners);
    std::print(report, "Corners by how many atoms the probe rests on: ");
    for (std::size_t touching = 3; touching < this->solventExcluded.cornersByContacts.size(); ++touching)
    {
      if (this->solventExcluded.cornersByContacts[touching] == 0) continue;
      std::print(report, "{} on {} ({:.3f} Å²)  ", this->solventExcluded.cornersByContacts[touching], touching,
                 this->solventExcluded.concaveByContacts[touching]);
    }
    std::print(report, "\n");
    std::print(report, "Corners with another within a thousandth of an Å: {}, closest pair {} Å\n",
               this->solventExcluded.crowdedCorners, this->solventExcluded.closestCorners);
    std::print(report,
               "Generated and then discarded as buried under a probe that could reach it: "
               "{} Å² convex, {} Å² saddle, {} Å² concave\n",
               this->solventExcluded.buriedConvexArea, this->solventExcluded.buriedSaddleArea,
               this->solventExcluded.buriedConcaveArea);
  }

  if (this->blockingAudit.numberOfSpheres > 0)
  {
    std::print(report, "Blocking spheres audited: {}, worst reach into a channel {} Å\n",
               this->blockingAudit.numberOfSpheres, this->blockingAudit.worstOverreach);
  }

  report.close();
}
