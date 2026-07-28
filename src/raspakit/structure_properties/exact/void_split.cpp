module;

module exact_void_split;

import std;

import int3;
import double3;
import voronoi_channels;
import voronoi_accessibility;
import exact_surface_patches;
import exact_union_volume;
import exact_boundary_components;

// How large a closure defect is still quadrature rather than a boundary that failed to close. The
// integral of the normal over a closed surface is zero, and what is left of it here is the error of the
// latitude quadrature, which on a structure whose pockets are properly labelled comes out at parts in
// 1e15 of the area. A defect orders above that means arcs bounding one pore were labelled with another,
// so that neither pore's boundary is whole and neither volume is that pore's.
constexpr double closureTolerance = 1.0e-6;

// How far a point is from the surfaces around it: the farthest one named surface gets from it, and the nearest
// any of the surfaces facing a channel comes to it.
//
// A patch lies on the sphere of one atom, so on it the distance to the point depends only on the direction, and
// it grows with the component of the direction along the line from the point to the atom's centre. So an
// extremum over a patch is at the direction straight toward or away from the point, where the patch holds it,
// and otherwise on the patch's boundary --- and along a bounding circle that component is one sinusoid in the
// angle, taking its extreme value at a turning point or at the end of the arc. Testing those candidates is
// exact, and nothing is stepped over.
//
// The two are taken in different frames, which is not an inconsistency but the difference between the two
// questions. The reach is about one surface, which closes only in its own frame, so its patches are carried
// there by the translations the decomposition accumulated. The distance to a channel is about real space, where
// what matters is the copy of an atom that happens to lie nearest, so those patches are taken at their nearest
// image.
struct SurfaceDistances
{
  double reach{0.0};
  double channel{std::numeric_limits<double>::max()};  // no infinities: the build turns them off
};

SurfaceDistances surfaceDistances(const VoronoiAccessibility& accessibility, const BoundaryComponents& components,
                                  const std::vector<std::uint8_t>& facesChannel, std::int32_t label,
                                  const double3& centre)
{
  SurfaceDistances result;
  for (std::size_t atomIndex = 0; atomIndex < components.atoms.size(); ++atomIndex)
  {
    const SphereBoundary& boundary = components.atoms[atomIndex];
    if (boundary.buried) continue;

    const double radius = accessibility.atomRadii[atomIndex];

    // The nearest copy of this atom, for the channel distance, and the least the whole of its sphere could come
    // to the point. Where that already exceeds the nearest patch found so far, no patch of this atom can beat
    // it, which is what keeps the search over the channels local.
    const double3 nearestImage =
        centre + accessibility.simulationBox.applyPeriodicBoundaryConditions(accessibility.atomPositions[atomIndex] -
                                                                            centre);
    const double sphereFloor = std::abs((nearestImage - centre).length() - radius);

    for (std::size_t patch = 0; patch < components.componentOfPatch[atomIndex].size(); ++patch)
    {
      const std::int32_t of = components.componentOfPatch[atomIndex][patch];
      const bool onSurface = (of == label);
      const bool onChannel = of >= 0 && facesChannel[static_cast<std::size_t>(of)] != 0 && sphereFloor < result.channel;
      if (!onSurface && !onChannel) continue;

      const int3 offset = components.offsetOfPatch[atomIndex][patch];
      const double3 inFrame =
          accessibility.atomPositions[atomIndex] +
          accessibility.simulationBox.cell *
              double3(static_cast<double>(offset.x), static_cast<double>(offset.y), static_cast<double>(offset.z));
      const double3 sphereCentre = onSurface ? inFrame : nearestImage;
      const double3 outward = sphereCentre - centre;

      auto consider = [&](const double3& unitDirection)
      {
        double distance = (sphereCentre + unitDirection * radius - centre).length();
        if (onSurface) result.reach = std::max(result.reach, distance);
        if (onChannel) result.channel = std::min(result.channel, distance);
      };

      // The two directions the sphere gets farthest from and nearest to the point in. Each has to be shown to be
      // exposed before it is asked which patch it is on: a direction is buried where it lies inside the cap some
      // neighbour cuts off, and a sphere left with a single patch answers that patch for every direction,
      // buried or not. Where a direction is not on the patch the extremum is on the patch's edges instead, and
      // taking it anyway would put the answer out at the whole sphere's, most of which is buried.
      double distance = outward.length();
      if (distance > 1.0e-12)
      {
        for (double sign : {1.0, -1.0})
        {
          double3 along = outward * (sign / distance);
          bool exposed = std::ranges::none_of(boundary.circles, [&](const SphereCircle& circle)
                                             { return double3::dot(along, circle.axis) > circle.cosineHalfAngle; });
          if (exposed && components.patchOfDirection(atomIndex, along) == static_cast<std::int32_t>(patch))
          {
            consider(along);
          }
        }
      }

      for (const SphereCircle& circle : boundary.circles)
      {
        for (std::size_t arc = 0; arc < circle.arcPatch.size(); ++arc)
        {
          if (circle.arcPatch[arc] != static_cast<std::int32_t>(patch)) continue;

          // A circle no crossing cuts is one arc closing on itself, where both turning points are in range.
          double begin = circle.cornerAngles.empty() ? 0.0 : circle.cornerAngles[arc];
          double end = circle.cornerAngles.empty()
                           ? 2.0 * std::numbers::pi
                           : circle.cornerAngles[(arc + 1) % circle.cornerAngles.size()];
          if (end <= begin) end += 2.0 * std::numbers::pi;

          consider(circle.direction(begin));
          consider(circle.direction(end));

          double turning = std::atan2(double3::dot(outward, circle.second), double3::dot(outward, circle.first));
          for (int half = 0; half < 2; ++half)
          {
            for (int lift = 0; lift < 3; ++lift)
            {
              double angle = turning + static_cast<double>(half) * std::numbers::pi +
                             static_cast<double>(lift - 1) * 2.0 * std::numbers::pi;
              if (angle >= begin && angle <= end) consider(circle.direction(angle));
            }
          }
        }
      }
    }
  }
  return result;
}

ExactVoidSplit exactVoidSplit(const VoronoiAccessibility& accessibility, const ExactSurfaceAreaSample& patches,
                              double cellVolume)
{
  ExactVoidSplit split;
  split.voidVolume = cellVolume - unionOfBallsVolume(accessibility, patches);
  split.undecidedArea = patches.undecided;

  std::size_t openPockets = 0;
  double worstArea = 0.0;
  for (std::size_t poreId = 0; poreId < patches.pores.size(); ++poreId)
  {
    const VoronoiPore& pore = accessibility.channels.pores[poreId];
    const PoreBoundaryMoments& moments = patches.pores[poreId];
    if (pore.isChannel || moments.area <= 0.0) continue;

    split.inaccessibleVolume -= (moments.radiusWeightedArea + moments.originWeighted) / 3.0;
    ++split.numberOfPockets;

    double defect = moments.vectorArea.length() / moments.area;
    if (defect > closureTolerance)
    {
      ++openPockets;
      if (defect > split.closureDefect) worstArea = moments.area;
    }
    split.closureDefect = std::max(split.closureDefect, defect);
  }

  // The channels take the rest. Nothing is measured twice, the void being what the union leaves and the
  // pockets being what the classifier sealed off within it.
  split.accessibleVolume = split.voidVolume - split.inaccessibleVolume;

  // A pocket volume can only come out negative, or larger than the void it is part of, through a boundary
  // that did not close, which the defect should already have caught; both are checked because either makes
  // the answer useless and neither should pass silently. The comparison against the void is given the slack
  // the two carry themselves, a structure whose void is all pockets having them equal and not deserving to
  // be rejected over the last digits of either.
  const double slack = 1.0e-6 * std::max(split.voidVolume, 1.0);
  if (split.undecidedArea > 0.0)
  {
    split.rejection = std::format("{} Å² of surface faces no pore", split.undecidedArea);
  }
  else if (openPockets > 0)
  {
    // The defect localises the disagreement, which is its use: it names the pocket whose boundary is
    // incomplete, and the area of that boundary says how much of the surface went to the wrong pore.
    split.rejection =
        std::format("{} of {} pocket boundaries do not close, the worst by {} of its {} Å²", openPockets,
                    split.numberOfPockets, split.closureDefect, worstArea);
  }
  else if (split.inaccessibleVolume < -slack)
  {
    split.rejection = std::format("a pocket volume came out negative, {} Å³", split.inaccessibleVolume);
  }
  else if (split.inaccessibleVolume > split.voidVolume + slack)
  {
    split.rejection = std::format("the pockets hold more than the void, {} against {} Å³", split.inaccessibleVolume,
                                  split.voidVolume);
  }
  split.reliable = split.rejection.empty();

  // Within the slack above, which is round-off rather than disagreement.
  if (split.reliable)
  {
    split.inaccessibleVolume = std::clamp(split.inaccessibleVolume, 0.0, split.voidVolume);
    split.accessibleVolume = split.voidVolume - split.inaccessibleVolume;
  }

  return split;
}

ExactVoidSplit exactVoidSplit(const VoronoiAccessibility& accessibility, double cellVolume,
                              std::size_t subdivisions)
{
  return exactVoidSplit(accessibility, exactAccessibleSurfaceArea(accessibility, subdivisions, true), cellVolume);
}


ExactVoidSplit exactVoidSplitByComponents(const VoronoiAccessibility& accessibility,
                                          const BoundaryComponents& components,
                                          const std::vector<ComponentLabel>& labels,
                                          const ExactSurfaceAreaSample& patches, double cellVolume)
{
  ExactVoidSplit split;
  split.voidVolume = cellVolume - unionOfBallsVolume(accessibility, patches);
  split.undecidedArea = patches.undecided;
  split.numberOfSurfaces = components.numberOfComponents;

  // The same points the sweep took its moments about, since the centroid below is read back against them.
  const std::vector<double3> origins = surfaceMomentOrigins(accessibility, components);

  // Which surfaces have the accessible void on their far side, by the same argument the area is divided by: a
  // surface that runs away through the crystal has the void running away with it, one that closes around void
  // seals it in, and one that closes around solid is a cluster whose surroundings only the network can place.
  // Taken first because a pocket needs it about every surface but its own.
  std::vector<std::uint8_t> facesChannel(components.numberOfComponents, 0);
  for (std::size_t label = 0; label < patches.components.size(); ++label)
  {
    const PoreBoundaryMoments& moments = patches.components[label];
    if (moments.area <= 0.0) continue;
    if (components.componentPercolates[label] != 0)
    {
      facesChannel[label] = 1;
    }
    else if (-(moments.radiusWeightedArea + moments.originWeighted) <= 0.0)
    {
      facesChannel[label] = (labels[label].decided && labels[label].accessible) ? 1 : 0;
    }
  }

  std::size_t openPockets = 0;
  double worstArea = 0.0;
  for (std::size_t label = 0; label < patches.components.size(); ++label)
  {
    const ComponentLabel& answer = labels[label];
    if (answer.proved) ++split.provedSurfaces;

    const PoreBoundaryMoments& moments = patches.components[label];
    if (moments.area <= 0.0) continue;

    // A surface that closes on a translate of itself bounds nothing at all, so it is neither a pocket nor a
    // correction to one, and there is nothing about it to check.
    if (components.componentPercolates[label] != 0) continue;

    // Whether it closes comes first, since the volume it encloses and even the sign of that volume mean
    // nothing until it does.
    double defect = moments.vectorArea.length() / moments.area;
    if (defect > closureTolerance)
    {
      ++openPockets;
      if (defect > split.closureDefect) worstArea = moments.area;
    }
    split.closureDefect = std::max(split.closureDefect, defect);

    double enclosed = -(moments.radiusWeightedArea + moments.originWeighted) / 3.0;
    if ((enclosed > 0.0) != (answer.decided && !answer.accessible))
    {
      ++split.signDisagreements;
      split.signDisagreementVolume += std::abs(enclosed);
    }

    if (enclosed > 0.0)
    {
      // Void, and enclosed, so sealed: a pocket, whatever the network makes of it.
      ++split.numberOfPockets;
      if (answer.proved) ++split.provedPockets;

      // Where the pocket is, from the same arcs. The moments were taken about a point of the surface itself,
      // carried into the surface's own frame, so the centroid comes back in that frame too and is brought
      // home before anything is asked about it.
      PocketGeometry pocket;
      pocket.volume = enclosed;
      pocket.area = moments.area;
      pocket.equivalentRadius = std::cbrt(3.0 * enclosed / (4.0 * std::numbers::pi));

      double3 centre = origins[label] + moments.enclosedFirstMoment * (1.0 / enclosed);
      pocket.centreFractional = double3::fract(accessibility.simulationBox.inverseCell * centre);
      pocket.centre = accessibility.simulationBox.cell * pocket.centreFractional;

      // The clearance is measured against every atom, and where the centre is inside the pocket the nearest of
      // them is one of the pocket's own walls: so it is the exact distance from the centre to the boundary and
      // not a bound on it, and the ball of it lies within the pocket. A pocket bent round a corner can have its
      // centroid outside itself, in the framework, where the clearance is negative, or in a channel beyond,
      // where a ball about it would block a pore. Neither leaves a ball to write, and both say so with a zero.
      pocket.centreInChannel = accessibility.provablyAccessible(pocket.centre);
      pocket.freeRadius =
          pocket.centreInChannel ? 0.0 : std::max(0.0, accessibility.clearance(pocket.centre));

      SurfaceDistances distances =
          surfaceDistances(accessibility, components, facesChannel, static_cast<std::int32_t>(label), centre);
      pocket.coveringRadius = distances.reach;
      pocket.channelRadius = distances.channel;
      split.pockets.push_back(pocket);
    }
    else
    {
      // Solid. It is room taken out of a pocket if it stands inside one, and nothing at all if it stands in
      // open void, which is the one thing here the geometry of this surface alone cannot say.
      if (!answer.decided || answer.accessible) continue;
      ++split.numberOfEnclosedSolids;
    }

    split.inaccessibleVolume += enclosed;
  }

  split.accessibleVolume = split.voidVolume - split.inaccessibleVolume;

  const double slack = 1.0e-6 * std::max(split.voidVolume, 1.0);
  if (split.undecidedArea > 0.0)
  {
    // Only a cluster of atoms standing in the void leaves the question to the network, so this is a cluster
    // the network could not place rather than an arc it could not place.
    split.rejection = std::format("{} Å² of surface stands round a cluster the network cannot place",
                                  split.undecidedArea);
  }
  else if (openPockets > 0)
  {
    // Where the surfaces are the pieces, a defect is no longer a disagreement between labels: it is either
    // quadrature or a surface the decomposition failed to join up, and the area names which one.
    split.rejection = std::format("{} bounded surfaces do not close, the worst by {} of its {} Å²", openPockets,
                                  split.closureDefect, worstArea);
  }
  else if (split.inaccessibleVolume < -slack)
  {
    split.rejection = std::format("a pocket volume came out negative, {} Å³", split.inaccessibleVolume);
  }
  else if (split.inaccessibleVolume > split.voidVolume + slack)
  {
    split.rejection = std::format("the pockets hold more than the void, {} against {} Å³", split.inaccessibleVolume,
                                  split.voidVolume);
  }
  split.reliable = split.rejection.empty();

  if (split.reliable)
  {
    split.inaccessibleVolume = std::clamp(split.inaccessibleVolume, 0.0, split.voidVolume);
    split.accessibleVolume = split.voidVolume - split.inaccessibleVolume;
  }

  return split;
}


ExactVoidSplit exactVoidSplitByComponents(const VoronoiAccessibility& accessibility, double cellVolume,
                                          std::size_t subdivisions)
{
  BoundaryComponents components = boundaryComponents(accessibility);
  std::vector<ComponentLabel> labels = labelBoundaryComponents(accessibility, components);
  ExactSurfaceAreaSample patches =
      exactAccessibleSurfaceAreaByComponent(accessibility, components, labels, subdivisions);
  return exactVoidSplitByComponents(accessibility, components, labels, patches, cellVolume);
}
