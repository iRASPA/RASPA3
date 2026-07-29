module;

module exact_surface_patches;

import std;

import double3;
import int3;
import voronoi_channels;
import voronoi_accessibility;
import exact_boundary_components;
import exact_sphere_sweep;


// Where an arc's side is read off, when it is read off the connected surface it belongs to rather than
// asked about arc by arc. `origins` is the point each surface's moments are taken about, in that surface's
// own frame.
struct ComponentRoute
{
  const BoundaryComponents* components{nullptr};
  const std::vector<double3>* origins{nullptr};
};


// The caps covering one atom's sphere, as the sweep wants them. Returns false where a neighbour swallows the
// sphere whole, in which case the atom carries no exposed surface at all.
//
// Where the boundary has already been decomposed, both the circles and the crossings among them are taken
// from it rather than found again. They are the same circles: the decomposition asks the same neighbour query
// over the same reach and prunes the contained discs by the same rule, so this is the same list in the same
// order, already built. What that saves is the neighbour query and the prune, and then the search for the
// crossings, which between them are the whole of the setup cost of a sphere. `crossings` is left null where
// there is no decomposition to take them from, and the sweep looks for them itself.
bool sweepCirclesOfAtom(const VoronoiAccessibility& accessibility, std::size_t atomIndex,
                        const BoundaryComponents* components, std::vector<SweepCircle>& circles,
                        std::vector<double3>& axes, const std::vector<double3>*& crossings)
{
  circles.clear();
  axes.clear();
  crossings = nullptr;

  if (components != nullptr)
  {
    const SphereBoundary& boundary = components->atoms[atomIndex];
    if (boundary.buried) return false;

    circles.reserve(boundary.circles.size());
    for (const SphereCircle& known : boundary.circles)
    {
      SweepCircle circle;
      circle.axis = known.axis;
      circle.cosineHalfAngle = known.cosineHalfAngle;
      circle.sineHalfAngle = known.sineHalfAngle;
      circle.halfAngle = known.halfAngle;
      circles.push_back(circle);
      axes.push_back(circle.axis);
    }
    crossings = &boundary.crossings;
    return true;
  }

  const double radius = accessibility.atomRadii[atomIndex];
  const double3 centre = accessibility.atomPositions[atomIndex];

  // Every sphere that can reach this one, its own periodic images included.
  for (const auto& [delta, neighbourRadius] :
       accessibility.neighbourAtoms(centre, radius + accessibility.maximumAtomRadius))
  {
    double distance = delta.length();
    if (distance < 1.0e-12)
    {
      // The sphere itself, which the query returns along with the rest, or another sitting on top of it.
      // Only a strictly larger one covers anything, and then it covers everything.
      if (neighbourRadius > radius) return false;
      continue;
    }

    // The circle of intersection sits at this polar angle from the neighbour's direction. At least one means
    // the neighbour reaches nothing of this sphere; at most minus one means it swallows it whole.
    double cosineHalfAngle =
        (radius * radius + distance * distance - neighbourRadius * neighbourRadius) / (2.0 * radius * distance);
    std::optional<SweepCircle> circle = makeSweepCircle(delta * (1.0 / distance), cosineHalfAngle);
    if (!circle.has_value()) continue;
    if (circle->cosineHalfAngle <= -1.0) return false;

    circles.push_back(circle.value());
    axes.push_back(circle->axis);
  }

  // A circle whose disc lies inside another's carries no boundary and no coverage of its own, and leaving it
  // in would only add latitudes at which nothing happens.
  if (circles.size() > 1)
  {
    pruneContainedDiscs(circles);
    axes.clear();
    for (const SweepCircle& circle : circles) axes.push_back(circle.axis);
  }
  return true;
}


// The exposed area of one probe-inflated sphere, added to `sample` split by the classifier.
void measureSphere(const VoronoiAccessibility& accessibility, std::size_t atomIndex, std::size_t subdivisions,
                   bool classifyArcs, const ComponentRoute* route, std::vector<SweepCircle>& circles,
                   SweepWorkspace& work, ExactSurfaceAreaSample& sample)
{
  const double radius = accessibility.atomRadii[atomIndex];
  const double3 centre = accessibility.atomPositions[atomIndex];

  // What an area on this sphere is worth on the bare sphere inside it, the two patches being the same set of
  // directions seen from the same centre.
  const double bareRatio = (radius > 0.0) ? accessibility.bareRadius(atomIndex) / radius : 0.0;
  const double bareWeight = bareRatio * bareRatio;

  // What an area on this sphere is worth as volume between it and the bare sphere: the cone from the centre
  // out to the inflated radius, less the cone to the bare one, over the same solid angle. It is accumulated
  // here rather than assembled from the totals because the weight is the atom's own and a sum over patches
  // cannot be reweighted after the fact.
  const double shellWeight = radius / 3.0 * (1.0 - bareRatio * bareRatio * bareRatio);

  const BoundaryComponents* known = (route != nullptr) ? route->components : nullptr;
  const std::vector<double3>* crossings = nullptr;
  if (!sweepCirclesOfAtom(accessibility, atomIndex, known, circles, work.axes, crossings)) return;

  const std::array<double3, 3> frame = sweepFrame(work.axes);
  const double3 firstAxis = frame[0];
  const double3 secondAxis = frame[1];
  const double3 polarAxis = frame[2];

  sweepExposedLatitudes(
      circles, frame, crossings, subdivisions, work,
      [&](const LatitudeGap& gap)
      {
        const double area = radius * radius * gap.sineLatitude * gap.span * gap.weight;
        ++sample.numberOfArcs;
        sample.area += area;
        sample.radiusWeightedArea += radius * area;
        sample.convexArea += bareWeight * area;
        sample.shellVolume += shellWeight * area;
        if (atomIndex < sample.atomArea.size()) sample.atomArea[atomIndex] += area;

        if (!classifyArcs) return;

        // The integral of the outward normal over this arc. Its azimuthal part is elementary, being the
        // difference of a sine and of a cosine between the ends of the same gap whose length gave the area,
        // so the moments cost three sums and no new geometry.
        const double3 normalIntegral = firstAxis * (gap.sineLatitude * (gap.sineEnd - gap.sineBegin)) +
                                       secondAxis * (gap.sineLatitude * (gap.cosineBegin - gap.cosineEnd)) +
                                       polarAxis * (gap.cosineLatitude * gap.span);
        const double3 arcVectorArea = normalIntegral * (radius * radius * gap.sineLatitude * gap.weight);

        // The azimuthal integrals of the products of two of the normal's components, which are what the first
        // moment of the enclosed region needs beyond the normal itself. They are the same two endpoints once
        // more, at twice the angle, which the double angle formulas give from what the gap already carries.
        const double sineTwiceEnd = 2.0 * gap.sineEnd * gap.cosineEnd;
        const double sineTwiceBegin = 2.0 * gap.sineBegin * gap.cosineBegin;
        const double cosineTwiceEnd = 2.0 * gap.cosineEnd * gap.cosineEnd - 1.0;
        const double cosineTwiceBegin = 2.0 * gap.cosineBegin * gap.cosineBegin - 1.0;

        const double spread = 0.25 * (sineTwiceEnd - sineTwiceBegin);
        const double integralCosine = gap.sineEnd - gap.sineBegin;
        const double integralSine = gap.cosineBegin - gap.cosineEnd;
        const double integralCosineCosine = 0.5 * gap.span + spread;
        const double integralSineSine = 0.5 * gap.span - spread;
        const double integralSineCosine = 0.25 * (cosineTwiceBegin - cosineTwiceEnd);

        // Everything one arc adds to the surface it belongs to, given where that surface's moments are taken
        // from. Both routes below add the same thing and differ only in which surface it goes to and which
        // copy of the cell the arc is carried into.
        auto addArcTo = [&](PoreBoundaryMoments& moments, const double3& shifted, const double3& origin)
        {
          const double3 delta = shifted - origin;

          // The tensor of the arc applied to `delta`, in the frame the sweep is done in.
          const double alongFirst = double3::dot(delta, firstAxis);
          const double alongSecond = double3::dot(delta, secondAxis);
          const double alongPolar = double3::dot(delta, polarAxis);
          double3 tensorTimesDelta =
              firstAxis * (gap.sineLatitude * gap.sineLatitude * alongFirst * integralCosineCosine +
                           gap.sineLatitude * gap.sineLatitude * alongSecond * integralSineCosine +
                           gap.sineLatitude * gap.cosineLatitude * alongPolar * integralCosine) +
              secondAxis * (gap.sineLatitude * gap.sineLatitude * alongFirst * integralSineCosine +
                            gap.sineLatitude * gap.sineLatitude * alongSecond * integralSineSine +
                            gap.sineLatitude * gap.cosineLatitude * alongPolar * integralSine) +
              polarAxis * (gap.cosineLatitude * gap.sineLatitude * alongFirst * integralCosine +
                           gap.cosineLatitude * gap.sineLatitude * alongSecond * integralSine +
                           gap.cosineLatitude * gap.cosineLatitude * alongPolar * gap.span);
          tensorTimesDelta = tensorTimesDelta * (radius * radius * gap.sineLatitude * gap.weight);

          moments.area += area;
          moments.radiusWeightedArea += radius * area;
          moments.convexArea += bareWeight * area;
          moments.shellVolume += shellWeight * area;
          moments.originWeighted += double3::dot(delta, arcVectorArea);
          moments.vectorArea += arcVectorArea;
          moments.enclosedFirstMoment +=
              arcVectorArea * (-0.5 * (double3::dot(delta, delta) + radius * radius)) - tensorTimesDelta * radius;
        };

        if (route != nullptr)
        {
          // Which patch the arc lies on. Either end of the gap is a point of that patch's own edge, where a
          // bounding circle crosses this latitude, and the arcs of the circles carry the patch they bound: so
          // the lookup is exact wherever there is such an end. A latitude no circle cuts has none, and only
          // there is a path on the sphere looked for instead.
          std::int32_t patch = -1;
          if (gap.begin > 0.0)
          {
            patch = route->components->patchOfCirclePoint(atomIndex, gap.atBegin(frame));
          }
          if (patch < 0 && gap.end < 2.0 * std::numbers::pi)
          {
            patch = route->components->patchOfCirclePoint(atomIndex, gap.atEnd(frame));
          }
          if (patch < 0)
          {
            patch = route->components->patchOfDirection(atomIndex, gap.atMiddle(frame));
          }

          const std::int32_t label =
              (patch < 0) ? -1 : route->components->componentOfPatch[atomIndex][static_cast<std::size_t>(patch)];
          if (label < 0)
          {
            sample.undecided += area;
            ++sample.unplacedArcs;
            sample.unplacedArea += area;
          }
          else
          {
            // Which side the arc is on is not settled here. It is settled for the whole surface once the
            // surface has been measured, since what decides it is the volume the surface encloses and that is
            // not known until then.

            // The arc is on an atom of the home cell, but the surface it belongs to runs through several
            // copies of the cell, and it closes only once the pieces are brought into one frame. Which
            // translation does that for this patch was settled when the surface was assembled, so nothing
            // here has to be inferred from where the arc happens to be.
            const int3 offset = route->components->offsetOfPatch[atomIndex][static_cast<std::size_t>(patch)];
            const double3 shifted =
                centre + accessibility.simulationBox.cell * double3(static_cast<double>(offset.x),
                                                                    static_cast<double>(offset.y),
                                                                    static_cast<double>(offset.z));

            addArcTo(sample.components[static_cast<std::size_t>(label)], shifted,
                     (*route->origins)[static_cast<std::size_t>(label)]);
          }
        }
        else
        {
          const double3 outward = gap.atMiddle(frame);
          const PointClassification classification = accessibility.classify(centre + outward * radius);

          if (classification.resample || classification.inside)
            sample.undecided += area;
          else if (classification.accessible)
            sample.accessible += area;
          else
            sample.inaccessible += area;

          if (classification.poreId >= 0)
          {
            // The arc faces one lift of its pore, and the pore's boundary closes only in the frame its own
            // nodes are in, so the arc is carried there before it is added. Taking the origin at the pore's
            // own first node rather than at some fixed point keeps |x - c| small, which matters: what
            // survives of the choice of origin is the closure defect times the distance to it.
            const VoronoiPore& pore = accessibility.channels.pores[static_cast<std::size_t>(classification.poreId)];
            const double3 origin = accessibility.network.nodes[pore.nodeIndices.front()].position;
            const double3 shifted =
                centre + accessibility.simulationBox.cell *
                             double3(static_cast<double>(classification.latticeOffset.x),
                                     static_cast<double>(classification.latticeOffset.y),
                                     static_cast<double>(classification.latticeOffset.z));

            addArcTo(sample.pores[static_cast<std::size_t>(classification.poreId)], shifted, origin);
          }
        }
      });
}


ExactSurfaceAreaSample exactAccessibleSurfaceArea(const VoronoiAccessibility& accessibility,
                                                  std::size_t subdivisions, bool classifyArcs)
{
  ExactSurfaceAreaSample sample;
  if (classifyArcs) sample.pores.assign(accessibility.channels.pores.size(), PoreBoundaryMoments{});
  sample.atomArea.assign(accessibility.atomPositions.size(), 0.0);

  std::vector<SweepCircle> circles;
  SweepWorkspace work;
  for (std::size_t i = 0; i < accessibility.atomPositions.size(); ++i)
  {
    measureSphere(accessibility, i, subdivisions, classifyArcs, nullptr, circles, work, sample);
  }
  return sample;
}


std::vector<double3> surfaceMomentOrigins(const VoronoiAccessibility& accessibility,
                                          const BoundaryComponents& components)
{
  // Each surface's moments are taken about a point of the surface itself, carried into the surface's own
  // frame. The choice cannot change a closed surface's volume; what it changes is how much of the closure
  // defect shows up in it, and a near point is what keeps that small.
  std::vector<double3> origins(components.numberOfComponents, double3(0.0, 0.0, 0.0));
  for (std::size_t label = 0; label < components.numberOfComponents; ++label)
  {
    const auto [atomIndex, patch] = components.componentRepresentative[label];
    if (patch >= components.offsetOfPatch[atomIndex].size()) continue;
    int3 offset = components.offsetOfPatch[atomIndex][patch];
    origins[label] = accessibility.atomPositions[atomIndex] +
                     accessibility.simulationBox.cell * double3(static_cast<double>(offset.x),
                                                                static_cast<double>(offset.y),
                                                                static_cast<double>(offset.z));
  }
  return origins;
}


ExactSurfaceAreaSample exactAccessibleSurfaceAreaByComponent(const VoronoiAccessibility& accessibility,
                                                             const BoundaryComponents& components,
                                                             const std::vector<ComponentLabel>& labels,
                                                             std::size_t subdivisions)
{
  ExactSurfaceAreaSample sample;
  sample.components.assign(components.numberOfComponents, PoreBoundaryMoments{});
  sample.atomArea.assign(accessibility.atomPositions.size(), 0.0);

  std::vector<double3> origins = surfaceMomentOrigins(accessibility, components);

  ComponentRoute route;
  route.components = &components;
  route.origins = &origins;

  std::vector<SweepCircle> circles;
  SweepWorkspace work;
  for (std::size_t i = 0; i < accessibility.atomPositions.size(); ++i)
  {
    measureSphere(accessibility, i, subdivisions, true, &route, circles, work, sample);
  }

  // Which side each surface is on, now that each has been measured whole. Two of the three cases are settled
  // by the surface itself. A surface that closes on a translate of itself runs away through the crystal, and
  // so does the void along it, so its area is reachable. A surface that closes on itself has the void on one
  // side and the solid on the other, and the sign of the volume it encloses says which: positive is void, and
  // void enclosed by a closed surface can reach nothing outside it, so that area is not reachable. Only a
  // surface enclosing solid --- a cluster of atoms with the void outside it --- leaves a question for the
  // network, which is whether the void around the cluster is itself sealed.
  for (std::size_t label = 0; label < components.numberOfComponents; ++label)
  {
    const PoreBoundaryMoments& moments = sample.components[label];
    if (moments.area <= 0.0) continue;

    ++sample.numberOfSurfaces;

    // How many directions this surface runs away in, which is the dimensionality of the pore behind it.
    std::size_t rank = static_cast<std::size_t>(std::clamp(components.componentDimensionality[label], 0, 3));
    ++sample.surfacesOfDimension[rank];
    sample.areaOfDimension[rank] += moments.area;

    bool reachable;
    if (components.componentPercolates[label] != 0)
    {
      ++sample.runawaySurfaces;
      reachable = true;
    }
    else if (-(moments.radiusWeightedArea + moments.originWeighted) > 0.0)
    {
      ++sample.sealedSurfaces;
      reachable = false;
    }
    else
    {
      ++sample.clusterSurfaces;
      const ComponentLabel& answer = labels[label];
      if (!answer.decided)
      {
        sample.undecided += moments.area;
        continue;
      }
      reachable = answer.accessible;
    }

    if (reachable)
      sample.accessible += moments.area;
    else
      sample.inaccessible += moments.area;
  }
  return sample;
}


ExactSurfaceAreaSample exactAccessibleSurfaceAreaByComponent(const VoronoiAccessibility& accessibility,
                                                             std::size_t subdivisions)
{
  BoundaryComponents components = boundaryComponents(accessibility);
  std::vector<ComponentLabel> labels = labelBoundaryComponents(accessibility, components);
  return exactAccessibleSurfaceAreaByComponent(accessibility, components, labels, subdivisions);
}
