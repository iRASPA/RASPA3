module;

module exact_void_split;

import std;

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
