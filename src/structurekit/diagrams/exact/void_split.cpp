module;

module exact_void_split;

import std;

import int3;
import double3;
import voronoi_channels;
import pore_accessibility;
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
// what matters is the copy of an atom that happens to lie nearest, so those patches are taken over every copy
// that can come near enough.
struct SurfaceDistances
{
  double reach{0.0};
  double channel{std::numeric_limits<double>::max()};  // no infinities: the build turns them off
};

// What the two routes below settle in the same way, once each has said the piece only it can say.
//
// `reason` is that piece: the route's own objection to its own division, and empty where it has none. What is
// left is the same question either way. A pocket volume can only come out negative, or larger than the void
// it is part of, through a boundary that did not close, which the defect should already have caught; both are
// checked because either makes the answer useless and neither should pass silently. The comparison against
// the void is given the slack the two carry themselves, a structure whose void is all pockets having them
// equal and not deserving to be rejected over the last digits of either.
void settleSplit(ExactVoidSplit& split, std::string reason)
{
  split.accessibleVolume = split.voidVolume - split.inaccessibleVolume;

  const double slack = 1.0e-6 * std::max(split.voidVolume, 1.0);
  if (!reason.empty())
  {
    split.rejection = std::move(reason);
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
}

SurfaceDistances surfaceDistances(const PoreAccessibility& accessibility, const BoundaryComponents& components,
                                  const std::vector<std::uint8_t>& facesChannel, std::int32_t component,
                                  const double3& centre)
{
  SurfaceDistances result;

  // The extrema of the distance from the centre to one patch of the sphere of radius `radius` about
  // `sphereCentre`, each handed to `consider` as the unit direction from the sphere's centre. Every patch is
  // the same on every copy of its atom, so this serves whichever copy is wanted.
  auto forEachExtremum = [&](std::size_t atomIndex, std::size_t patch, const double3& sphereCentre, auto&& consider)
  {
    const SphereBoundary& boundary = components.atoms[atomIndex];
    const double3 outward = sphereCentre - centre;

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
  };

  auto facesAChannel = [&](std::int32_t of) { return of >= 0 && facesChannel[static_cast<std::size_t>(of)] != 0; };

  for (std::size_t atomIndex = 0; atomIndex < components.atoms.size(); ++atomIndex)
  {
    const SphereBoundary& boundary = components.atoms[atomIndex];
    if (boundary.buried) continue;

    const double radius = accessibility.atomRadii[atomIndex];

    // The pocket's own surface, in the frame it was measured in.
    for (std::size_t patch = 0; patch < components.componentOfPatch[atomIndex].size(); ++patch)
    {
      if (components.componentOfPatch[atomIndex][patch] != component) continue;
      const int3 offset = components.offsetOfPatch[atomIndex][patch];
      const double3 inFrame =
          accessibility.atomPositions[atomIndex] +
          accessibility.unitCell.cell *
              double3(static_cast<double>(offset.x), static_cast<double>(offset.y), static_cast<double>(offset.z));
      forEachExtremum(atomIndex, patch, inFrame, [&](const double3& unitDirection)
                      { result.reach = std::max(result.reach, (inFrame + unitDirection * radius - centre).length()); });
    }

    // A first channel distance, from the copy of this atom that the fractional rounding puts nearest. The least
    // the whole of that sphere could come to the point says whether it is worth asking at all.
    const double3 nearestImage =
        centre + accessibility.unitCell.applyPeriodicBoundaryConditions(accessibility.atomPositions[atomIndex] -
                                                                            centre);
    const double sphereFloor = std::abs((nearestImage - centre).length() - radius);
    if (sphereFloor >= result.channel) continue;
    for (std::size_t patch = 0; patch < components.componentOfPatch[atomIndex].size(); ++patch)
    {
      if (!facesAChannel(components.componentOfPatch[atomIndex][patch])) continue;
      forEachExtremum(atomIndex, patch, nearestImage, [&](const double3& unitDirection)
                      { result.channel = std::min(result.channel, (nearestImage + unitDirection * radius - centre).length()); });
    }
  }

  // Every copy of every atom whose sphere can come nearer than that, which the first pass does not see: a patch
  // that faces away from the point on the nearest copy faces it on the next, and in a skewed cell the copy the
  // rounding picks is not always the nearest. The copies within the distance found so far (plus a radius) are
  // all that can improve on it, and the walk over the cells finds exactly those.
  if (result.channel < 0.5 * std::numeric_limits<double>::max())
  {
    const double reach = result.channel + accessibility.maximumAtomRadius;
    for (const NeighbourImage& image : accessibility.neighbourAtomImages(centre, reach))
    {
      if (components.atoms[image.index].buried) continue;
      const double radius = accessibility.atomRadii[image.index];
      if (image.delta.length() - radius >= result.channel) continue;
      const double3 sphereCentre = centre + image.delta;
      for (std::size_t patch = 0; patch < components.componentOfPatch[image.index].size(); ++patch)
      {
        if (!facesAChannel(components.componentOfPatch[image.index][patch])) continue;
        forEachExtremum(image.index, patch, sphereCentre, [&](const double3& unitDirection)
                        { result.channel = std::min(result.channel, (sphereCentre + unitDirection * radius - centre).length()); });
      }
    }
  }
  return result;
}

ExactVoidSplit exactVoidSplit(const PoreAccessibility& accessibility, const MeasuredPatches& patches,
                              double cellVolume)
{
  ExactVoidSplit split;
  split.voidVolume = cellVolume - unionOfBallsVolume(accessibility, patches);
  split.undecidedArea = patches.undecided;

  std::size_t openBoundaries = 0;
  double worstArea = 0.0;
  for (std::size_t poreId = 0; poreId < patches.pores.size(); ++poreId)
  {
    const VoronoiPore& pore = accessibility.channels.pores[poreId];
    const BoundaryMoments& moments = patches.pores[poreId];
    if (pore.isChannel || moments.area <= 0.0) continue;

    split.inaccessibleVolume -= (moments.radiusWeightedArea + moments.originWeighted) / 3.0;
    ++split.numberOfPockets;

    double defect = moments.vectorArea.length() / moments.area;
    if (defect > closureTolerance)
    {
      ++openBoundaries;
      if (defect > split.closureDefect) worstArea = moments.area;
    }
    split.closureDefect = std::max(split.closureDefect, defect);
  }

  // The channels take the rest. Nothing is measured twice, the void being what the union leaves and the
  // pockets being what the classifier sealed off within it.
  std::string reason;
  if (split.undecidedArea > 0.0)
  {
    reason = std::format("{} Å² of surface faces no pore", split.undecidedArea);
  }
  else if (openBoundaries > 0)
  {
    // The defect localises the disagreement, which is its use: it names the pocket whose boundary is
    // incomplete, and the area of that boundary says how much of the surface went to the wrong pore.
    reason = std::format("{} of {} pocket boundaries do not close, the worst by {} of its {} Å²", openBoundaries,
                         split.numberOfPockets, split.closureDefect, worstArea);
  }
  settleSplit(split, std::move(reason));
  return split;
}

ExactVoidSplit exactVoidSplit(const PoreAccessibility& accessibility, double cellVolume,
                              std::size_t subdivisions)
{
  return exactVoidSplit(accessibility, exactAccessibleSurfaceAreaByPore(accessibility, subdivisions), cellVolume);
}


ExactVoidSplit exactVoidSplitByComponents(const PoreAccessibility& accessibility,
                                          const BoundaryComponents& components,
                                          const std::vector<ComponentVerdict>& verdicts,
                                          const MeasuredPatches& patches, double cellVolume)
{
  // Every pocket here is reported with a centre, and the centre is the first moment of the region divided by
  // its volume, so patches measured without one have nothing to give and would report every pocket at the
  // point its moments were taken about. Said here rather than discovered there.
  if (patches.moments != SurfaceMoments::andCentre)
  {
    throw std::runtime_error(
        "exactVoidSplitByComponents: the patches were measured without the first moment, so the pockets have no "
        "centre; sweep them with SurfaceMoments::andCentre\n");
  }

  ExactVoidSplit split;
  split.voidVolume = cellVolume - unionOfBallsVolume(accessibility, patches, components);
  split.undecidedArea = patches.undecided;
  split.numberOfSurfaces = components.numberOfComponents;

  // The same points the sweep took its moments about, since the centroid below is read back against them.
  const std::vector<double3> origins = surfaceMomentOrigins(accessibility, components);

  // Which surfaces have the accessible void on their far side, by the same division the area and the excluded
  // volume use. Taken first because a pocket needs it about every surface but its own.
  const std::vector<SurfaceSide> sides = surfaceSides(components, patches, verdicts, patches.sidePolicy);
  std::vector<std::uint8_t> facesChannel(components.numberOfComponents, 0);
  for (std::size_t component = 0; component < components.numberOfComponents; ++component)
  {
    facesChannel[component] = (sides[component].side > 0) ? 1 : 0;
  }

  // The shortest vector from `from` to any copy of `point`, which is the cell's to find: in a skewed cell it is
  // not always the copy nearest in fractions.
  auto nearestDelta = [&](const double3& point, const double3& from)
  { return accessibility.unitCell.applyPeriodicBoundaryConditions(point - from); };

  // Everything about a pocket that is read off its centre, from a centre given in the frame its surface was
  // measured in. Done once with the centroid of the wall alone, and again for a pocket with a cluster inside
  // it, whose centroid is elsewhere.
  auto placePocket = [&](PocketGeometry& pocket, std::size_t component, const double3& centre)
  {
    pocket.equivalentRadius = std::cbrt(3.0 * pocket.volume / (4.0 * std::numbers::pi));
    pocket.centreFractional = double3::fract(accessibility.unitCell.inverseCell * centre);
    pocket.centre = accessibility.unitCell.cell * pocket.centreFractional;

    // The clearance is measured against every atom, and where the centre is inside the pocket the nearest of
    // them is one of the pocket's own walls: so it is the exact distance from the centre to the boundary and
    // not a bound on it, and the ball of it lies within the pocket. A pocket bent round a corner can have its
    // centroid outside itself, in the framework, where the clearance is negative, or in a channel beyond,
    // where a ball about it would block a pore. Neither leaves a ball to write, and both say so with a zero.
    pocket.centreInChannel = accessibility.provablyAccessible(pocket.centre);
    pocket.freeRadius = pocket.centreInChannel ? 0.0 : std::max(0.0, accessibility.clearance(pocket.centre));

    SurfaceDistances distances =
        surfaceDistances(accessibility, components, facesChannel, static_cast<std::int32_t>(component), centre);
    pocket.coveringRadius = distances.reach;
    pocket.channelRadius = distances.channel;
  };

  // A cluster of atoms standing inside a pocket. Its surface encloses solid, and what it encloses is room the
  // pocket does not have; but which pocket that is can only be asked once the pockets are known, so the
  // clusters are kept until then.
  struct EnclosedSolid
  {
    double volume;         // Å³, of the solid, positive
    double area;           // Å², of its surface
    double3 centre;        // its centroid, in the frame its surface was measured in
    double3 surfacePoint;  // a point of its surface, in the same frame
    std::int32_t poreId;   // the pore the network says the void round it belongs to
  };
  std::vector<EnclosedSolid> solids;
  std::vector<std::size_t> pocketComponent;  // per pocket, in the order of `split.pockets`
  std::vector<double3> pocketCentre;         // per pocket, the wall's centroid in the surface's own frame
  std::vector<std::int32_t> pocketPoreId;

  std::size_t openBoundaries = 0;
  double worstArea = 0.0;
  for (std::size_t component = 0; component < patches.components.size(); ++component)
  {
    const ComponentVerdict& answer = verdicts[component];
    if (answer.proved) ++split.provedSurfaces;

    const BoundaryMoments& moments = patches.components[component];
    if (moments.area <= 0.0) continue;

    // A surface that closes on a translate of itself bounds nothing at all, so it is neither a pocket nor a
    // correction to one, and there is nothing about it to check.
    if (components.componentPercolates[component] != 0) continue;

    // Whether it closes comes first, since the volume it encloses and even the sign of that volume mean
    // nothing until it does.
    double defect = moments.vectorArea.length() / moments.area;
    if (defect > closureTolerance)
    {
      ++openBoundaries;
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

      double3 centre = origins[component] + moments.enclosedFirstMoment * (1.0 / enclosed);
      placePocket(pocket, component, centre);
      split.pockets.push_back(pocket);
      pocketComponent.push_back(component);
      pocketCentre.push_back(centre);
      pocketPoreId.push_back(answer.decided ? answer.poreId : -1);
    }
    else
    {
      // Solid. It is room taken out of a pocket if it stands inside one, and nothing at all if it stands in
      // open void, which is the one thing here the geometry of this surface alone cannot say.
      if (!answer.decided || answer.accessible) continue;
      ++split.numberOfEnclosedSolids;

      const auto [atomIndex, patch] = components.componentRepresentative[component];
      const SphereBoundary& boundary = components.atoms[atomIndex];
      double3 surfacePoint = origins[component];
      if (patch < boundary.patchRepresentative.size())
      {
        surfacePoint += boundary.patchRepresentative[patch] * accessibility.atomRadii[atomIndex];
      }
      solids.push_back(EnclosedSolid{-enclosed, moments.area,
                                     origins[component] + moments.enclosedFirstMoment * (1.0 / enclosed),
                                     surfacePoint, answer.poreId});
    }

    split.inaccessibleVolume += enclosed;
  }

  // Each cluster charged to the pocket it stands in. The cluster lies inside the region the pocket's wall
  // bounds, and the farthest point of a bounded region from anywhere is on its boundary, so a point of the
  // cluster's surface is within the wall's own reach of the wall's centroid: an exact test, which most clusters
  // pass for one pocket only. Where the reaches of two pockets overlap, the network says which, the cluster's
  // surface and the wall round it facing the same pore. A cluster neither settles is counted and left where it
  // is, in the total but in no pocket, and the blocking spheres refuse the pockets on that count.
  std::vector<std::vector<std::size_t>> solidsOfPocket(split.pockets.size());
  for (std::size_t solid = 0; solid < solids.size(); ++solid)
  {
    std::vector<std::size_t> candidates;
    for (std::size_t p = 0; p < split.pockets.size(); ++p)
    {
      const double within = nearestDelta(solids[solid].surfacePoint, pocketCentre[p]).length();
      const double reach = split.pockets[p].coveringRadius;
      if (within <= reach + 1.0e-9 * std::max(reach, 1.0)) candidates.push_back(p);
    }
    if (candidates.size() > 1 && solids[solid].poreId >= 0)
    {
      std::erase_if(candidates, [&](std::size_t p) { return pocketPoreId[p] != solids[solid].poreId; });
    }
    if (candidates.size() == 1)
    {
      solidsOfPocket[candidates.front()].push_back(solid);
    }
    else
    {
      ++split.unplacedSolids;
    }
  }

  for (std::size_t p = 0; p < split.pockets.size(); ++p)
  {
    if (solidsOfPocket[p].empty()) continue;

    // The region is the wall's less the clusters', so its first moment is the wall's less theirs, each cluster
    // taken at the copy that stands inside this wall. What is read off the centre is then read off again.
    PocketGeometry& pocket = split.pockets[p];
    double3 moment = pocketCentre[p] * pocket.volume;
    for (std::size_t solid : solidsOfPocket[p])
    {
      const EnclosedSolid& cluster = solids[solid];
      const double3 inside = pocketCentre[p] + nearestDelta(cluster.centre, pocketCentre[p]);
      moment -= inside * cluster.volume;
      pocket.volume -= cluster.volume;
      pocket.area += cluster.area;
    }
    if (pocket.volume > 0.0)
    {
      placePocket(pocket, pocketComponent[p], moment * (1.0 / pocket.volume));
    }
  }

  std::string reason;
  if (split.undecidedArea > 0.0)
  {
    // Only a cluster of atoms standing in the void leaves the question to the network, so this is a cluster
    // the network could not place rather than an arc it could not place.
    reason = std::format("{} Å² of surface stands round a cluster the network cannot place", split.undecidedArea);
  }
  else if (openBoundaries > 0)
  {
    // Where the surfaces are the pieces, a defect is no longer a disagreement between verdicts: it is either
    // quadrature or a surface the decomposition failed to join up, and the area names which one.
    reason = std::format("{} bounded surfaces do not close, the worst by {} of its {} Å²", openBoundaries,
                         split.closureDefect, worstArea);
  }
  settleSplit(split, std::move(reason));
  return split;
}


ExactVoidSplit exactVoidSplitByComponents(const PoreAccessibility& accessibility, double cellVolume,
                                          std::size_t subdivisions)
{
  BoundaryComponents components = boundaryComponents(accessibility);
  std::vector<ComponentVerdict> verdicts = boundaryComponentVerdicts(accessibility, components);
  MeasuredPatches patches = exactAccessibleSurfaceAreaByComponent(accessibility, components, verdicts,
                                                                  subdivisions, SurfaceMoments::andCentre);
  return exactVoidSplitByComponents(accessibility, components, verdicts, patches, cellVolume);
}
