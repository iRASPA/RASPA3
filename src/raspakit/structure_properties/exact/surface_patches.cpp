module;

module exact_surface_patches;

import std;

import double3;
import int3;
import voronoi_channels;
import pore_accessibility;
import exact_boundary_components;
import exact_parallel;
import exact_sphere_sweep;


// What the sweep does with each arc beyond measuring it, which is the one thing its three callers differ in.
// The geometry below is the same in all three, and so is every total it accumulates.
enum class ArcAttribution : std::uint8_t
{
  // Nothing. The area and the moments of the whole boundary are wanted and no division of them, which is
  // what a volume asks for and the whole of what can be had without consulting a network at all.
  none,

  // Ask the network about each arc on its own, and collect the arcs pore by pore. See the note on
  // `exactAccessibleSurfaceAreaByPore`: this is the reference route and the expensive one.
  perArc,

  // Read the side off the connected surface the arc lies on, and collect the arcs surface by surface. This
  // is the route everything outside the tests takes, and the only one that needs a `ComponentRoute`.
  perSurface,
};

// Where an arc's side is read off under `perSurface`, and nothing under the other two. `origins` is the point
// each surface's moments are taken about, in that surface's own frame.
struct ComponentRoute
{
  const BoundaryComponents* components{nullptr};
  const std::vector<double3>* origins{nullptr};
};


// The caps covering one atom's sphere, as the sweep wants them. Returns false where a neighbour swallows the
// sphere whole, in which case the atom carries no exposed surface at all.
//
// Where the boundary has already been decomposed, both the circles and the crossings among them are taken
// from it rather than found again. They are the same circles: the decomposition enumerates the caps through
// `eachCapCoveringAtom`, as the branch below does, and prunes the contained discs by the same rule, so this
// is the same list in the same order, already built. What that saves is the neighbour query and the prune,
// and then the search for the crossings, which between them are the whole of the setup cost of a sphere.
// `crossings` is left null where there is no decomposition to take them from, and the sweep looks for them
// itself.
bool sweepCirclesOfAtom(const PoreAccessibility& accessibility, std::size_t atomIndex,
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

  const bool exposed = eachCapCoveringAtom(accessibility, atomIndex,
                                           [&](const double3& axis, double cosineHalfAngle, std::size_t, int3)
                                           {
                                             SweepCircle circle;
                                             circle.axis = axis;
                                             circle.cosineHalfAngle = cosineHalfAngle;
                                             circle.halfAngle = std::acos(cosineHalfAngle);
                                             circle.sineHalfAngle = std::sin(circle.halfAngle);
                                             circles.push_back(circle);
                                           });
  if (!exposed) return false;

  // A circle whose disc lies inside another's carries no boundary and no coverage of its own, and leaving it
  // in would only add latitudes at which nothing happens.
  pruneContainedDiscs(circles);
  for (const SweepCircle& circle : circles) axes.push_back(circle.axis);
  return true;
}


// One sphere, as the sweep of it wants it: the atom, the frame the sweep runs in, and the two weights that
// carry an area on this inflated sphere back to the bare sphere inside it.
struct SweptSphere
{
  const PoreAccessibility& accessibility;
  std::size_t atomIndex{0};
  double3 centre{0.0, 0.0, 0.0};
  double radius{0.0};
  double bareWeight{0.0};   // an area here, as area on the bare sphere
  double shellWeight{0.0};  // and as the volume between the two
  std::array<double3, 3> frame{double3{1.0, 0.0, 0.0}, double3{0.0, 1.0, 0.0}, double3{0.0, 0.0, 1.0}};
  SurfaceMoments moments{SurfaceMoments::volume};
};

// One exposed arc, reduced to everything a surface can want of it that does not depend on which surface that
// is or where the surface's moments are taken from. What remains is `addTo`, and it is the same computation on
// both of the routes below: they differ only in which surface the arc goes to and which copy of the cell it is
// carried into.
struct SweptArc
{
  double area{0.0};
  double3 vectorArea{0.0, 0.0, 0.0};

  double sineLatitude{0.0};
  double cosineLatitude{0.0};
  double span{0.0};
  double weight{0.0};

  // The azimuthal integrals of the normal's components and of the products of two of them, which are what the
  // first moment of the enclosed region needs beyond the normal itself. Left at zero, and not paid for, unless
  // the sweep was asked for the first moment.
  double integralCosine{0.0};
  double integralSine{0.0};
  double integralCosineCosine{0.0};
  double integralSineSine{0.0};
  double integralSineCosine{0.0};

  void addTo(const SweptSphere& sphere, BoundaryMoments& moments, const double3& shifted,
             const double3& origin) const;
};

// The area of a gap is the one thing every route wants of it, and the one thing the cheapest route wants at
// all, so it is taken separately and handed back in.
SweptArc sweptArc(const SweptSphere& sphere, const LatitudeGap& gap, double area)
{
  const double3& firstAxis = sphere.frame[0];
  const double3& secondAxis = sphere.frame[1];
  const double3& polarAxis = sphere.frame[2];

  SweptArc arc;
  arc.area = area;
  arc.sineLatitude = gap.sineLatitude;
  arc.cosineLatitude = gap.cosineLatitude;
  arc.span = gap.span;
  arc.weight = gap.weight;

  // The integral of the outward normal over this arc. Its azimuthal part is elementary, being the difference
  // of a sine and of a cosine between the ends of the same gap whose length gave the area, so the moments cost
  // three sums and no new geometry.
  const double3 normalIntegral = firstAxis * (gap.sineLatitude * (gap.sineEnd - gap.sineBegin)) +
                                 secondAxis * (gap.sineLatitude * (gap.cosineBegin - gap.cosineEnd)) +
                                 polarAxis * (gap.cosineLatitude * gap.span);
  arc.vectorArea = normalIntegral * (sphere.radius * sphere.radius * gap.sineLatitude * gap.weight);

  if (sphere.moments == SurfaceMoments::volume) return arc;

  // The products are the same two endpoints once more, at twice the angle, which the double angle formulas
  // give from what the gap already carries.
  const double sineTwiceEnd = 2.0 * gap.sineEnd * gap.cosineEnd;
  const double sineTwiceBegin = 2.0 * gap.sineBegin * gap.cosineBegin;
  const double cosineTwiceEnd = 2.0 * gap.cosineEnd * gap.cosineEnd - 1.0;
  const double cosineTwiceBegin = 2.0 * gap.cosineBegin * gap.cosineBegin - 1.0;

  const double spread = 0.25 * (sineTwiceEnd - sineTwiceBegin);
  arc.integralCosine = gap.sineEnd - gap.sineBegin;
  arc.integralSine = gap.cosineBegin - gap.cosineEnd;
  arc.integralCosineCosine = 0.5 * gap.span + spread;
  arc.integralSineSine = 0.5 * gap.span - spread;
  arc.integralSineCosine = 0.25 * (cosineTwiceBegin - cosineTwiceEnd);
  return arc;
}

// Everything this arc adds to the surface it belongs to, given where that surface's moments are taken from.
void SweptArc::addTo(const SweptSphere& sphere, BoundaryMoments& moments, const double3& shifted,
                     const double3& origin) const
{
  const double radius = sphere.radius;
  const double3 delta = shifted - origin;

  moments.area += area;
  moments.radiusWeightedArea += radius * area;
  moments.convexArea += sphere.bareWeight * area;
  moments.shellVolume += sphere.shellWeight * area;
  moments.originWeighted += double3::dot(delta, vectorArea);
  moments.vectorArea += vectorArea;

  if (sphere.moments == SurfaceMoments::volume) return;

  const double3& firstAxis = sphere.frame[0];
  const double3& secondAxis = sphere.frame[1];
  const double3& polarAxis = sphere.frame[2];

  // The tensor of the arc applied to `delta`, in the frame the sweep is done in.
  const double alongFirst = double3::dot(delta, firstAxis);
  const double alongSecond = double3::dot(delta, secondAxis);
  const double alongPolar = double3::dot(delta, polarAxis);
  double3 tensorTimesDelta = firstAxis * (sineLatitude * sineLatitude * alongFirst * integralCosineCosine +
                                          sineLatitude * sineLatitude * alongSecond * integralSineCosine +
                                          sineLatitude * cosineLatitude * alongPolar * integralCosine) +
                             secondAxis * (sineLatitude * sineLatitude * alongFirst * integralSineCosine +
                                           sineLatitude * sineLatitude * alongSecond * integralSineSine +
                                           sineLatitude * cosineLatitude * alongPolar * integralSine) +
                             polarAxis * (cosineLatitude * sineLatitude * alongFirst * integralCosine +
                                          cosineLatitude * sineLatitude * alongSecond * integralSine +
                                          cosineLatitude * cosineLatitude * alongPolar * span);
  tensorTimesDelta = tensorTimesDelta * (radius * radius * sineLatitude * weight);

  moments.enclosedFirstMoment +=
      vectorArea * (-0.5 * (double3::dot(delta, delta) + radius * radius)) - tensorTimesDelta * radius;
}

// The arc put to the connected surface it lies on, which the decomposition already knows.
void attributeToSurface(const SweptSphere& sphere, const SweptArc& arc, const LatitudeGap& gap,
                        const ComponentRoute& route, MeasuredPatches& patches)
{
  const std::size_t atomIndex = sphere.atomIndex;

  // Which patch the arc lies on. Either end of the gap is a point of that patch's own edge, where a bounding
  // circle crosses this latitude, and the arcs of the circles carry the patch they bound: so the lookup is
  // exact wherever there is such an end. A latitude no circle cuts has none, and only there is a path on the
  // sphere looked for instead.
  std::int32_t patch = -1;
  if (gap.begin > 0.0)
  {
    patch = route.components->patchOfCirclePoint(atomIndex, gap.atBegin(sphere.frame));
  }
  if (patch < 0 && gap.end < 2.0 * std::numbers::pi)
  {
    patch = route.components->patchOfCirclePoint(atomIndex, gap.atEnd(sphere.frame));
  }
  if (patch < 0)
  {
    patch = route.components->patchOfDirection(atomIndex, gap.atMiddle(sphere.frame));
  }

  const std::int32_t component =
      (patch < 0) ? -1 : route.components->componentOfPatch[atomIndex][static_cast<std::size_t>(patch)];
  if (component < 0)
  {
    patches.undecided += arc.area;
    ++patches.diagnostics.unplacedArcs;
    patches.diagnostics.unplacedArea += arc.area;
    return;
  }

  // Which side the arc is on is not settled here. It is settled for the whole surface once the surface has
  // been measured, since what decides it is the volume the surface encloses and that is not known until then.

  // The arc is on an atom of the home cell, but the surface it belongs to runs through several copies of the
  // cell, and it closes only once the pieces are brought into one frame. Which translation does that for this
  // patch was settled when the surface was assembled, so nothing here has to be inferred from where the arc
  // happens to be.
  const int3 offset = route.components->offsetOfPatch[atomIndex][static_cast<std::size_t>(patch)];
  const double3 shifted =
      sphere.centre + sphere.accessibility.simulationBox.cell * double3(static_cast<double>(offset.x),
                                                                        static_cast<double>(offset.y),
                                                                        static_cast<double>(offset.z));

  arc.addTo(sphere, patches.components[static_cast<std::size_t>(component)], shifted,
            (*route.origins)[static_cast<std::size_t>(component)]);
}

// The arc put to the classifier on its own, and to whichever pore of the network that lands it in.
void attributeToPore(const SweptSphere& sphere, const SweptArc& arc, const LatitudeGap& gap,
                     MeasuredPatches& patches)
{
  const PoreAccessibility& accessibility = sphere.accessibility;
  const double3 outward = gap.atMiddle(sphere.frame);
  const PointClassification classification = accessibility.classify(sphere.centre + outward * sphere.radius);

  if (classification.resample || classification.inside)
    patches.undecided += arc.area;
  else if (classification.accessible)
    patches.accessible += arc.area;
  else
    patches.inaccessible += arc.area;

  if (classification.poreId < 0) return;

  // The arc faces one lift of its pore, and the pore's boundary closes only in the frame its own nodes are in,
  // so the arc is carried there before it is added. Taking the origin at the pore's own first node rather than
  // at some fixed point keeps |x - c| small, which matters: what survives of the choice of origin is the
  // closure defect times the distance to it.
  const VoronoiPore& pore = accessibility.channels.pores[static_cast<std::size_t>(classification.poreId)];
  const double3 origin = accessibility.network.nodes[pore.nodeIndices.front()].position;
  const double3 shifted = sphere.centre + accessibility.simulationBox.cell *
                                              double3(static_cast<double>(classification.latticeOffset.x),
                                                      static_cast<double>(classification.latticeOffset.y),
                                                      static_cast<double>(classification.latticeOffset.z));

  arc.addTo(sphere, patches.pores[static_cast<std::size_t>(classification.poreId)], shifted, origin);
}

// The exposed area of one probe-inflated sphere, added to `patches` and attributed as `attribution` says.
// `route` is given exactly under `ArcAttribution::perSurface` and is null otherwise.
void measureSphere(const PoreAccessibility& accessibility, std::size_t atomIndex, std::size_t subdivisions,
                   ArcAttribution attribution, const ComponentRoute* route, SurfaceMoments wanted,
                   std::vector<SweepCircle>& circles, SweepWorkspace& work, MeasuredPatches& patches)
{
  const double radius = accessibility.atomRadii[atomIndex];

  const BoundaryComponents* known = (route != nullptr) ? route->components : nullptr;
  const std::vector<double3>* crossings = nullptr;
  if (!sweepCirclesOfAtom(accessibility, atomIndex, known, circles, work.axes, crossings)) return;

  // What an area on this sphere is worth on the bare sphere inside it, the two patches being the same set of
  // directions seen from the same centre; and what it is worth as volume between the two, which is the cone
  // from the centre out to the inflated radius less the cone to the bare one, over the same solid angle. The
  // volume is accumulated arc by arc rather than assembled from the totals because the weight is the atom's
  // own and a sum over patches cannot be reweighted after the fact.
  const double bareRatio = (radius > 0.0) ? accessibility.bareRadius(atomIndex) / radius : 0.0;

  const SweptSphere sphere{.accessibility = accessibility,
                           .atomIndex = atomIndex,
                           .centre = accessibility.atomPositions[atomIndex],
                           .radius = radius,
                           .bareWeight = bareRatio * bareRatio,
                           .shellWeight = radius / 3.0 * (1.0 - bareRatio * bareRatio * bareRatio),
                           .frame = sweepFrame(work.axes),
                           .moments = wanted};

  sweepExposedLatitudes(
      circles, sphere.frame, crossings, subdivisions, work,
      [&](const LatitudeGap& gap)
      {
        const double area = radius * radius * gap.sineLatitude * gap.span * gap.weight;
        ++patches.diagnostics.numberOfArcs;
        patches.area += area;
        patches.radiusWeightedArea += radius * area;
        if (atomIndex < patches.atomArea.size()) patches.atomArea[atomIndex] += area;

        if (attribution == ArcAttribution::none) return;

        const SweptArc arc = sweptArc(sphere, gap, area);
        if (attribution == ArcAttribution::perSurface)
        {
          attributeToSurface(sphere, arc, gap, *route, patches);
        }
        else
        {
          attributeToPore(sphere, arc, gap, patches);
        }
      });
}


void addMomentsTo(BoundaryMoments& into, const BoundaryMoments& part)
{
  into.area += part.area;
  into.radiusWeightedArea += part.radiusWeightedArea;
  into.originWeighted += part.originWeighted;
  into.vectorArea += part.vectorArea;
  into.convexArea += part.convexArea;
  into.shellVolume += part.shellVolume;
  into.enclosedFirstMoment += part.enclosedFirstMoment;
}

// One worker's share of a sweep, added to another's. Everything the sweep accumulates and nothing else: the
// fields the caller fills before or after the sweep are none of this function's business, and the two
// partials were sized alike from the same blank.
void addSweptTo(MeasuredPatches& into, const MeasuredPatches& part)
{
  into.accessible += part.accessible;
  into.inaccessible += part.inaccessible;
  into.undecided += part.undecided;
  into.area += part.area;
  into.radiusWeightedArea += part.radiusWeightedArea;

  into.diagnostics.numberOfArcs += part.diagnostics.numberOfArcs;
  into.diagnostics.unplacedArcs += part.diagnostics.unplacedArcs;
  into.diagnostics.unplacedArea += part.diagnostics.unplacedArea;

  for (std::size_t i = 0; i < into.atomArea.size(); ++i) into.atomArea[i] += part.atomArea[i];
  for (std::size_t i = 0; i < into.pores.size(); ++i) addMomentsTo(into.pores[i], part.pores[i]);
  for (std::size_t i = 0; i < into.components.size(); ++i) addMomentsTo(into.components[i], part.components[i]);
}

// Every sphere swept, into a copy of `blank` for each worker, and the copies added together in worker order.
// On one worker that is exactly the loop it has always been: arc by arc into a single accumulator, and no
// reduction to reassociate anything. On more than one the answer can move in its last digits, which is what
// a reduction costs and what `exact_parallel` says about it.
MeasuredPatches sweepAtoms(const PoreAccessibility& accessibility, std::size_t subdivisions,
                           ArcAttribution attribution, const ComponentRoute* route, SurfaceMoments wanted,
                           const MeasuredPatches& blank)
{
  const std::size_t atoms = accessibility.atomPositions.size();
  const std::size_t workers = workersAvailable();

  if (workers == 1)
  {
    MeasuredPatches patches = blank;
    std::vector<SweepCircle> circles;
    SweepWorkspace work;
    for (std::size_t i = 0; i < atoms; ++i)
    {
      measureSphere(accessibility, i, subdivisions, attribution, route, wanted, circles, work, patches);
    }
    return patches;
  }

  std::vector<MeasuredPatches> partials(workers, blank);
  std::vector<std::vector<SweepCircle>> circles(workers);
  std::vector<SweepWorkspace> work(workers);
  forEachIndex(atoms, workers,
               [&](std::size_t worker, std::size_t atom)
               {
                 measureSphere(accessibility, atom, subdivisions, attribution, route, wanted, circles[worker],
                               work[worker], partials[worker]);
               });

  MeasuredPatches patches = std::move(partials.front());
  for (std::size_t worker = 1; worker < workers; ++worker) addSweptTo(patches, partials[worker]);
  return patches;
}

// The sweep behind both of the routes that consult no decomposition, which differ only in whether each arc
// is put to the classifier on its own.
MeasuredPatches sweepEveryAtom(const PoreAccessibility& accessibility, std::size_t subdivisions,
                               ArcAttribution attribution)
{
  MeasuredPatches blank;
  if (attribution == ArcAttribution::perArc)
  {
    blank.pores.assign(accessibility.channels.pores.size(), BoundaryMoments{});
  }
  blank.atomArea.assign(accessibility.atomPositions.size(), 0.0);

  // The area route collects no moments at all and the arc-by-arc route collects them pore by pore, where the
  // volume is read and the centroid is not: nothing asks either of these for a first moment.
  return sweepAtoms(accessibility, subdivisions, attribution, nullptr, SurfaceMoments::volume, blank);
}


MeasuredPatches exactSurfaceArea(const PoreAccessibility& accessibility, std::size_t subdivisions)
{
  return sweepEveryAtom(accessibility, subdivisions, ArcAttribution::none);
}


MeasuredPatches exactAccessibleSurfaceAreaByPore(const PoreAccessibility& accessibility,
                                                 std::size_t subdivisions)
{
  return sweepEveryAtom(accessibility, subdivisions, ArcAttribution::perArc);
}


std::vector<SurfaceSide> surfaceSides(const BoundaryComponents& components, const MeasuredPatches& patches,
                                      const std::vector<ComponentVerdict>& verdicts)
{
  std::vector<SurfaceSide> sides(components.numberOfComponents, SurfaceSide{});
  for (std::size_t component = 0; component < components.numberOfComponents; ++component)
  {
    if (component >= patches.components.size()) continue;
    const BoundaryMoments& moments = patches.components[component];
    if (moments.area <= 0.0) continue;

    if (components.componentPercolates[component] != 0)
    {
      sides[component] = {SurfaceClosure::runaway, 1};
    }
    else if (-(moments.radiusWeightedArea + moments.originWeighted) > 0.0)
    {
      sides[component] = {SurfaceClosure::sealed, -1};
    }
    else
    {
      const ComponentVerdict& verdict = verdicts[component];
      sides[component] = {SurfaceClosure::cluster, !verdict.decided ? 0 : (verdict.accessible ? 1 : -1)};
    }
  }
  return sides;
}


std::vector<double3> surfaceMomentOrigins(const PoreAccessibility& accessibility,
                                          const BoundaryComponents& components)
{
  // Each surface's moments are taken about a point of the surface itself, carried into the surface's own
  // frame. The choice cannot change a closed surface's volume; what it changes is how much of the closure
  // defect shows up in it, and a near point is what keeps that small.
  std::vector<double3> origins(components.numberOfComponents, double3(0.0, 0.0, 0.0));
  for (std::size_t component = 0; component < components.numberOfComponents; ++component)
  {
    const auto [atomIndex, patch] = components.componentRepresentative[component];
    if (patch >= components.offsetOfPatch[atomIndex].size()) continue;
    int3 offset = components.offsetOfPatch[atomIndex][patch];
    origins[component] = accessibility.atomPositions[atomIndex] +
                     accessibility.simulationBox.cell * double3(static_cast<double>(offset.x),
                                                                static_cast<double>(offset.y),
                                                                static_cast<double>(offset.z));
  }
  return origins;
}


MeasuredPatches exactAccessibleSurfaceAreaByComponent(const PoreAccessibility& accessibility,
                                                      const BoundaryComponents& components,
                                                      const std::vector<ComponentVerdict>& verdicts,
                                                      std::size_t subdivisions, SurfaceMoments wanted)
{
  MeasuredPatches blank;
  blank.components.assign(components.numberOfComponents, BoundaryMoments{});
  blank.atomArea.assign(accessibility.atomPositions.size(), 0.0);
  blank.moments = wanted;

  std::vector<double3> origins = surfaceMomentOrigins(accessibility, components);

  ComponentRoute route;
  route.components = &components;
  route.origins = &origins;

  MeasuredPatches patches =
      sweepAtoms(accessibility, subdivisions, ArcAttribution::perSurface, &route, wanted, blank);

  // Which side each surface is on, now that each has been measured whole. The counts alongside say what
  // settled each of them, which is how much of this division is geometry's and how much is the network's.
  const std::vector<SurfaceSide> sides = surfaceSides(components, patches, verdicts);
  for (std::size_t component = 0; component < components.numberOfComponents; ++component)
  {
    const SurfaceSide& side = sides[component];
    if (side.closure == SurfaceClosure::empty) continue;

    const BoundaryMoments& moments = patches.components[component];
    ++patches.numberOfSurfaces;

    // How many directions this surface runs away in, which is the dimensionality of the pore behind it.
    std::size_t rank = static_cast<std::size_t>(std::clamp(components.componentDimensionality[component], 0, 3));
    ++patches.surfacesOfDimension[rank];
    patches.areaOfDimension[rank] += moments.area;

    switch (side.closure)
    {
      case SurfaceClosure::runaway:
        ++patches.runawaySurfaces;
        break;
      case SurfaceClosure::sealed:
        ++patches.sealedSurfaces;
        break;
      case SurfaceClosure::cluster:
        ++patches.clusterSurfaces;
        break;
      case SurfaceClosure::empty:
        break;
    }

    if (side.side > 0)
      patches.accessible += moments.area;
    else if (side.side < 0)
      patches.inaccessible += moments.area;
    else
      patches.undecided += moments.area;
  }
  return patches;
}


MeasuredPatches exactAccessibleSurfaceAreaByComponent(const PoreAccessibility& accessibility,
                                                      std::size_t subdivisions, SurfaceMoments wanted)
{
  BoundaryComponents components = boundaryComponents(accessibility);
  std::vector<ComponentVerdict> verdicts = boundaryComponentVerdicts(accessibility, components);
  return exactAccessibleSurfaceAreaByComponent(accessibility, components, verdicts, subdivisions, wanted);
}
