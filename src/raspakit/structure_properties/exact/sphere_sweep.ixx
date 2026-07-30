module;

export module exact_sphere_sweep;

import std;

import double3;

// Sweeping the exposed part of one sphere, latitude by latitude.
//
// Three of the exact analyses need the part of a sphere that no other sphere covers. The surface area sweeps
// the boundary of the union of the inflated atoms; the solvent excluded geometry sweeps the probe's own
// sphere at a vertex, clipped by the neighbouring probes; and the decomposition into patches cuts the same
// sphere along the same circles to find out what is connected to what. What each measures differs, and where
// each measures it does not: that is what lives here.
//
// The region is a sphere less a set of spherical caps. Cut it into circles of latitude in a frame chosen so
// that no cap sits on the pole, and on each circle the covered part is a union of arcs, one per cap that
// reaches that latitude, each centred on that cap's own azimuth. What is left between them is exposed, in
// closed form, and the integral over latitude is smooth between the latitudes at which the arcs appear,
// vanish, or run into one another. Those latitudes are the breakpoints, and between them a Gauss-Legendre
// rule is exact to round-off.

// How many nodes the latitude rule uses on each half of each smooth piece. Ten integrate a polynomial of
// degree nineteen exactly, and what the pieces and the endpoint substitution leave to integrate is analytic
// on each of them, so ten already puts the error at round-off. Raising it buys nothing and costs whatever the
// caller does per node. Exported to be reported with a result rather than to be chosen by a caller.
export inline constexpr std::size_t exactQuadratureOrder = 10uz;

// Gaps in a circle of latitude shorter than this are dropped. They are the seams where two caps meet almost
// tangentially: they carry no area worth measuring and are the one place where the closed form for an arc
// loses digits.
export inline constexpr double sweepGapTolerance = 1.0e-12;

// Covered has to mean covered with room to spare. A framework is symmetric, so three spheres meeting in one
// point is the ordinary case rather than a coincidence, and a crossing that a third cap happens to pass
// exactly through is still a crossing. Deciding that on the sign of a rounding error loses it, and with it an
// edge of a patch or a latitude at which an integrand turns. A crossing wrongly kept costs one more panel of
// quadrature; a crossing wrongly dropped costs digits.
export inline constexpr double capCoverTolerance = 1.0e-9;

// How long a panel of the latitude rule may be where it sits next to a pole, as a multiple of the room
// between the piece it belongs to and that pole. See `panelBoundaries`, where the number is used and argued.
export inline constexpr double poleClearance = 4.0;

// A unit vector perpendicular to `axis`, chosen so that the cross product behind it is well conditioned.
export double3 perpendicularTo(const double3& axis);

// A polar angle folded back into [0, pi], which is where the extreme latitudes of a cap live: a cap whose
// axis lies at polar angle beta turns back at |beta - theta| and at beta + theta, and the second of those
// reaches past the far pole and comes back when beta + theta exceeds pi.
export double foldedPolarAngle(double angle);

// The frame a sphere is swept in, given the axes of the caps covering it. Latitude slicing degenerates for a
// cap whose axis sits on the polar axis, so the polar axis is chosen to be as far as possible from every one
// of them, out of a fixed set of candidates: the answer must not depend on which run computed it. Returned as
// `{first, second, polar}`, the first two being where the azimuth is measured from.
export std::array<double3, 3> sweepFrame(const std::vector<double3>& axes);

// Whether the disc of the inner cap lies within the disc of the outer, so that the inner one bounds nothing
// of its own and covers nothing the outer does not.
//
// The condition is that the angle between the axes plus the inner half angle is at most the outer one, taken
// on the cosines rather than on the angles: `acos` is the most expensive thing in a sweep, and this is asked
// of every pair of caps on every sphere.
export inline bool discWithinDisc(double cosineBetweenAxes, double cosineInner, double sineInner,
                                  double cosineOuter, double sineOuter)
{
  if (cosineOuter > cosineInner) return false;  // the outer half angle has to be the larger
  return cosineBetweenAxes >= cosineOuter * cosineInner + sineOuter * sineInner;
}

// Drops the caps whose discs lie inside another's. Anything carrying `axis`, `cosineHalfAngle` and
// `sineHalfAngle` will do, which is what lets the decomposition and the two sweeps share this although each
// keeps a different amount of other detail on its own circles.
export template <typename Circle>
void pruneContainedDiscs(std::vector<Circle>& circles)
{
  if (circles.size() < 2) return;

  std::vector<bool> redundant(circles.size(), false);
  for (std::size_t i = 0; i < circles.size(); ++i)
  {
    for (std::size_t j = 0; j < circles.size(); ++j)
    {
      if (i == j || redundant[j]) continue;
      if (discWithinDisc(double3::dot(circles[i].axis, circles[j].axis), circles[i].cosineHalfAngle,
                         circles[i].sineHalfAngle, circles[j].cosineHalfAngle, circles[j].sineHalfAngle))
      {
        redundant[i] = true;
        break;
      }
    }
  }

  std::size_t kept = 0;
  for (std::size_t i = 0; i < circles.size(); ++i)
  {
    if (!redundant[i]) circles[kept++] = circles[i];
  }
  circles.resize(kept);
}

// A point where two of the caps cross that no third one covers, with the two caps it belongs to.
//
// These are the only places where the exposed part of a circle can begin or end: a crossing hidden inside a
// third cap has that cap's surface on both sides of it and interrupts nothing. They are the corners of the
// exposed region, and so both the corners of its patches and the latitudes at which the length of an exposed
// circle of latitude stops being analytic.
export struct CapCrossing
{
  std::size_t firstCircle{0};
  std::size_t secondCircle{0};
  double3 direction{0.0, 0.0, 0.0};
};

// Every crossing of two of the caps that no third one covers, collected into `crossings`. Generic over the
// circle type for the same reason as the prune above; only `axis` and `cosineHalfAngle` are read.
export template <typename Circle>
void uncoveredCrossings(const std::vector<Circle>& circles, std::vector<CapCrossing>& crossings)
{
  crossings.clear();
  for (std::size_t j = 0; j + 1 < circles.size(); ++j)
  {
    for (std::size_t k = j + 1; k < circles.size(); ++k)
    {
      const Circle& first = circles[j];
      const Circle& second = circles[k];
      const double alignment = double3::dot(first.axis, second.axis);
      const double denominator = 1.0 - alignment * alignment;
      if (denominator < 1.0e-14) continue;  // parallel axes: concentric circles never cross

      const double alongFirst = (first.cosineHalfAngle - alignment * second.cosineHalfAngle) / denominator;
      const double alongSecond = (second.cosineHalfAngle - alignment * first.cosineHalfAngle) / denominator;
      const double outOfPlaneSquared =
          (1.0 - alongFirst * first.cosineHalfAngle - alongSecond * second.cosineHalfAngle) / denominator;
      if (outOfPlaneSquared <= 0.0) continue;  // the circles miss one another

      const double3 inPlane = first.axis * alongFirst + second.axis * alongSecond;
      const double3 outOfPlane = double3::cross(first.axis, second.axis) * std::sqrt(outOfPlaneSquared);
      for (std::size_t side = 0; side < 2; ++side)
      {
        double3 direction = (side == 0) ? inPlane + outOfPlane : inPlane - outOfPlane;
        direction = direction * (1.0 / direction.length());

        bool covered = false;
        for (std::size_t l = 0; l < circles.size() && !covered; ++l)
        {
          if (l == j || l == k) continue;
          covered = double3::dot(direction, circles[l].axis) > circles[l].cosineHalfAngle + capCoverTolerance;
        }
        if (!covered) crossings.push_back(CapCrossing{j, k, direction});
      }
    }
  }
}

// The same, for a caller that asks for them once rather than for every sphere of a structure.
export template <typename Circle>
std::vector<CapCrossing> uncoveredCrossings(const std::vector<Circle>& circles)
{
  std::vector<CapCrossing> crossings;
  uncoveredCrossings(circles, crossings);
  return crossings;
}

// One cap of the sphere being swept, placed in the frame the sweep is done in.
//
// The two latitudes are where the cap's boundary turns back on itself: between them it crosses each circle of
// latitude in a genuine arc, outside them the latitude is either wholly inside the cap or wholly outside it.
// Which of the two depends only on whether the cap reaches over the pole on that side, so the pair of flags
// settles it. That is what lets a whole piece of the integral be skipped where some cap buries it, and lets
// the rest be integrated against a handful of caps rather than against all of them.
export struct SweepCircle
{
  double3 axis{0.0, 0.0, 0.0};
  double cosineHalfAngle{0.0};
  double sineHalfAngle{0.0};
  double halfAngle{0.0};

  double polarAngle{0.0};  // of the axis, in the sweep frame
  double cosinePolar{0.0};
  double sinePolar{0.0};
  double azimuth{0.0};  // of the axis, in the same frame
  double cosineAzimuth{1.0};
  double sineAzimuth{0.0};

  double lowestLatitude{0.0};
  double highestLatitude{0.0};
  bool reachesOverPole{false};
  bool reachesOverAntipole{false};
};

// A cap from an axis and the cosine of its half angle, or nothing where the cap covers none of the sphere.
// A cap covering all of it comes back too and is the caller's to notice, by `cosineHalfAngle <= -1`: what an
// entirely covered sphere means is the one thing the callers do not agree on.
export std::optional<SweepCircle> makeSweepCircle(const double3& axis, double cosineHalfAngle);

// One exposed stretch of one circle of latitude, as the sweep hands it to whatever is being measured.
//
// The endpoints come with their own cosine and sine, and those are not recovered from the angles. The arcs a
// cap covers are built by turning the cap's own azimuth through the half width, which is an angle addition on
// quantities the sweep already holds, so the direction at either end of a gap costs a few multiplications and
// no trigonometry at all. Everything the moments of an arc need is a product of these.
export struct LatitudeGap
{
  double sineLatitude{0.0};
  double cosineLatitude{0.0};
  double weight{0.0};  // of the latitude rule, so that an integral is a sum of `value * weight`

  double begin{0.0};  // azimuth of the near end
  double end{0.0};    // and of the far one
  double span{0.0};   // `end - begin`, the angular length of the gap

  double cosineBegin{1.0};
  double sineBegin{0.0};
  double cosineEnd{1.0};
  double sineEnd{0.0};

  // The direction at a given azimuth of this latitude, in the frame the sweep is done in.
  double3 at(const std::array<double3, 3>& frame, double cosineAzimuth, double sineAzimuth) const
  {
    return frame[0] * (sineLatitude * cosineAzimuth) + frame[1] * (sineLatitude * sineAzimuth) +
           frame[2] * cosineLatitude;
  }
  double3 atBegin(const std::array<double3, 3>& frame) const { return at(frame, cosineBegin, sineBegin); }
  double3 atEnd(const std::array<double3, 3>& frame) const { return at(frame, cosineEnd, sineEnd); }

  // The direction at the middle of the gap. The only thing here that costs trigonometry, so it is left to
  // the callers that want it rather than taken for every gap.
  double3 atMiddle(const std::array<double3, 3>& frame) const
  {
    const double middle = begin + 0.5 * span;
    return at(frame, std::cos(middle), std::sin(middle));
  }
};

// One arc of a circle of latitude that a cap covers, with the direction at either end carried alongside the
// angle so that the gaps between the arcs need no trigonometry of their own.
export struct CoveredArc
{
  double begin{0.0};
  double end{0.0};
  double cosineBegin{1.0};
  double sineBegin{0.0};
  double cosineEnd{1.0};
  double sineEnd{0.0};
};

// Scratch the sweep needs, kept by the caller so that a structure's worth of spheres costs one allocation
// rather than one per sphere.
export struct SweepWorkspace
{
  std::vector<double3> axes;            // of the caps, for choosing the frame
  std::vector<CapCrossing> crossings;   // where two of them cross uncovered
  std::vector<double> breakpoints;      // the latitudes between which the integrand is analytic
  std::vector<std::size_t> cutting;     // the caps that cut the piece being integrated
  std::vector<double> panels;           // the latitudes that piece is cut at for the quadrature
  std::vector<CoveredArc> covered;      // of one circle of latitude
};

// Places every cap in the frame, puts them in order of their own azimuth, and collects the latitudes at which
// the exposed length of a circle of latitude stops being analytic: where a cap's boundary turns back, and
// where two of them cross uncovered. The breakpoints are left sorted and always bracketed by 0 and pi.
//
// `knownCrossings`, where it is given, is that second set of latitudes already found, which spares the search
// for them: it is cubic in the number of caps, and a caller that has decomposed the same sphere into patches
// has been over every one of those pairs on its own account. It has to be the crossings of these very caps.
//
// It finds the same crossings, and finds each of them to within a couple of the last bits: which cap of a
// pair is reached first decides which way round the products in a cross product are written, and the compiler
// contracts those into fused multiply-adds, which round one of the two and not the other. What that moves is
// where a panel of the quadrature begins and ends, never how much region there is to integrate, and the
// panels are cut finer than they need to be already. Two neighbouring panels in place of one is the whole of
// the consequence.
export void prepareSweep(std::vector<SweepCircle>& circles, const std::array<double3, 3>& frame,
                         const std::vector<double3>* knownCrossings, SweepWorkspace& work);

// Gauss-Legendre nodes and weights on the unit interval, found once by Newton's method on the Legendre
// polynomial with the usual Chebyshev starting guess.
export struct GaussRule
{
  std::array<double, exactQuadratureOrder> nodes{};
  std::array<double, exactQuadratureOrder> weights{};
};

export const GaussRule& unitIntervalGaussRule();

// The latitudes one smooth piece is cut at, into `subdivisions` panels of equal length and then graded
// towards an end that stands close to a pole. Left in `panels`, first entry `begin` and last `end`.
//
// Each half of a panel is anchored at an end and integrated against a squared variable, which is what makes
// the square root a cap's width leaves at that end smooth, and it is why one panel per piece is enough
// nearly everywhere. What the substitution cannot reach is a singularity lying just outside the panel, and
// there is one at each pole: the half width of an arc is a ratio with the sine of the latitude beneath it,
// so the integrand carries on past the end of a piece towards a pole the piece never reaches and turns
// singular where that sine vanishes.
//
// A piece can end very close to a pole. The frame keeps the pole away from the cap axes, but a cap turns
// back at its axis give or take its half angle, and nothing in the choice of frame looks at the half angle:
// a cap whose axis stands a radian and a half off the pole can still turn back a fortieth of a radian from
// it. A panel reaching from there across the whole piece is then integrating a function whose nearest
// singularity is thirty times nearer than the panel is long, and that costs three or four digits of the
// thirteen the rule otherwise gives.
//
// Halving the end panel until it is no longer than `poleClearance` times the room left to the pole puts the
// singularity back out of the rule's reach, and it does so for the panels behind it too: each is about as
// long as its own distance from the pole, so the ratio that matters is the same on all of them. It is
// geometric in the room left, so a few halvings settle it however little there is, and it is done only
// where the room is small --- a piece clear of both poles keeps exactly the panels the caller asked for and
// the arithmetic it had before. `cut` says whether any cap cuts this piece; where none does there is
// nothing but the sine of the latitude to integrate and that is analytic at the poles as anywhere else.
//
// Four is the multiple because the loss sets in below about eight and is gone by four, measured against the
// closed form for a single cap turned to put a turning latitude at every distance from a pole.
export void panelBoundaries(double begin, double end, std::size_t subdivisions, bool cut,
                            std::vector<double>& panels);

// The covered arcs of one latitude, ordered by where they begin.
//
// They arrive very nearly in that order already: the caps are taken in order of their own azimuth and each
// covers an arc centred on it, so only the one or two arcs that run over zero are out of place. An insertion
// sort passes over such a list once, where a general sort pays its setup every time, and this is called once
// per node of the quadrature, which is the innermost thing there is here.
export inline void sortByBeginning(std::vector<CoveredArc>& arcs)
{
  for (std::size_t i = 1; i < arcs.size(); ++i)
  {
    const CoveredArc held = arcs[i];
    std::size_t j = i;
    while (j > 0 && arcs[j - 1].begin > held.begin)
    {
      arcs[j] = arcs[j - 1];
      --j;
    }
    arcs[j] = held;
  }
}

// Walks the exposed part of the sphere and calls `measure(gap)` for every exposed stretch of every circle of
// latitude the quadrature visits. `subdivisions` cuts each smooth piece into that many panels of equal
// length; one settles an area to some thirteen digits and raising it is a check on that rather than an
// improvement to it. Panels beyond those are added where a piece ends near a pole, whatever the caller asks
// for, since there one long panel is not enough on its own: see `panelBoundaries`.
// `knownCrossings` is as in `prepareSweep`, and null for a caller with nothing already found.
//
// The whole of the geometry is here and none of what is being measured: the area of a gap is
// `radius * radius * gap.sineLatitude * gap.span * gap.weight`, and its moments are products of the cosines
// and sines the gap carries.
export template <typename Measure>
void sweepExposedLatitudes(std::vector<SweepCircle>& circles, const std::array<double3, 3>& frame,
                           const std::vector<double3>* knownCrossings, std::size_t subdivisions,
                           SweepWorkspace& work, Measure&& measure)
{
  prepareSweep(circles, frame, knownCrossings, work);

  const GaussRule& rule = unitIntervalGaussRule();
  const double twoPi = 2.0 * std::numbers::pi;
  const std::size_t parts = std::max<std::size_t>(1, subdivisions);

  for (std::size_t piece = 0; piece + 1 < work.breakpoints.size(); ++piece)
  {
    const double pieceBegin = work.breakpoints[piece];
    const double pieceEnd = work.breakpoints[piece + 1];
    if (pieceEnd - pieceBegin < 1.0e-14) continue;

    // Which caps cut the latitudes of this piece, and whether any of them covers them entirely. Both are
    // settled anywhere inside the piece, its ends being exactly the latitudes at which either can change, so
    // they are settled once here rather than at every node of the quadrature.
    const double interior = 0.5 * (pieceBegin + pieceEnd);
    work.cutting.clear();
    bool buried = false;
    for (std::size_t i = 0; i < circles.size(); ++i)
    {
      const SweepCircle& circle = circles[i];
      if (interior > circle.lowestLatitude && interior < circle.highestLatitude)
      {
        work.cutting.push_back(i);
      }
      else if ((interior <= circle.lowestLatitude && circle.reachesOverPole) ||
               (interior >= circle.highestLatitude && circle.reachesOverAntipole))
      {
        buried = true;
        break;
      }
    }
    if (buried) continue;

    panelBoundaries(pieceBegin, pieceEnd, parts, !work.cutting.empty(), work.panels);

    for (std::size_t panel = 0; panel + 1 < work.panels.size(); ++panel)
    {
      const double begin = work.panels[panel];
      const double end = work.panels[panel + 1];
      const double middle = 0.5 * (begin + end);

      // The uncovered length of a latitude leaves a piece end like the square root of the distance to it, a
      // cap appearing or vanishing there with a width of that shape. Anchoring the two halves at the ends and
      // substituting a square makes the integrand smooth in the new variable, which is what leaves the
      // quadrature error at round-off rather than at the square root of it.
      for (std::size_t half = 0; half < 2; ++half)
      {
        const double anchor = (half == 0) ? begin : end;
        const double span = (half == 0) ? middle - begin : end - middle;
        const double direction = (half == 0) ? 1.0 : -1.0;

        for (std::size_t node = 0; node < exactQuadratureOrder; ++node)
        {
          const double parameter = rule.nodes[node];
          const double latitude = anchor + direction * span * parameter * parameter;
          const double sineLatitude = std::sin(latitude);
          if (sineLatitude <= 0.0) continue;

          LatitudeGap gap;
          gap.sineLatitude = sineLatitude;
          gap.cosineLatitude = std::cos(latitude);
          gap.weight = 2.0 * span * parameter * rule.weights[node];

          // Each cutting cap covers one arc of this latitude, centred on its own azimuth and reaching a half
          // width to either side. The direction at an end of that arc is the cap's own direction turned
          // through the half width, which is an angle addition and no trigonometry of its own: the cap
          // carries the cosine and sine of its azimuth, and the half width brings its own. The caps are
          // visited in order of their azimuth, so the arcs come out sorted but for the ones running over
          // zero, which are entered as their two pieces to keep the sweep below on [0, 2pi).
          work.covered.clear();
          for (std::size_t i : work.cutting)
          {
            const SweepCircle& circle = circles[i];
            const double cosineHalfWidth = (circle.cosineHalfAngle - gap.cosineLatitude * circle.cosinePolar) /
                                           (sineLatitude * circle.sinePolar);
            if (cosineHalfWidth >= 1.0) continue;  // the cap does not reach this latitude after all

            // At or below minus one the cap covers the whole of the latitude, which is the half width being
            // a straight angle and the turn through it a reflection.
            const bool whole = cosineHalfWidth <= -1.0;
            const double halfWidth = whole ? std::numbers::pi : std::acos(cosineHalfWidth);
            const double cosineHalfWidthClamped = whole ? -1.0 : cosineHalfWidth;
            const double sineHalfWidth =
                whole ? 0.0 : std::sqrt(std::max(0.0, 1.0 - cosineHalfWidth * cosineHalfWidth));

            CoveredArc arc;
            arc.begin = circle.azimuth - halfWidth;
            arc.end = circle.azimuth + halfWidth;
            arc.cosineBegin = circle.cosineAzimuth * cosineHalfWidthClamped + circle.sineAzimuth * sineHalfWidth;
            arc.sineBegin = circle.sineAzimuth * cosineHalfWidthClamped - circle.cosineAzimuth * sineHalfWidth;
            arc.cosineEnd = circle.cosineAzimuth * cosineHalfWidthClamped - circle.sineAzimuth * sineHalfWidth;
            arc.sineEnd = circle.sineAzimuth * cosineHalfWidthClamped + circle.cosineAzimuth * sineHalfWidth;

            // Split at zero where it runs over, the direction there being along the first axis of the frame.
            if (arc.begin < 0.0)
            {
              work.covered.push_back(CoveredArc{arc.begin + twoPi, twoPi, arc.cosineBegin, arc.sineBegin, 1.0, 0.0});
              work.covered.push_back(CoveredArc{0.0, arc.end, 1.0, 0.0, arc.cosineEnd, arc.sineEnd});
            }
            else if (arc.end > twoPi)
            {
              work.covered.push_back(CoveredArc{arc.begin, twoPi, arc.cosineBegin, arc.sineBegin, 1.0, 0.0});
              work.covered.push_back(CoveredArc{0.0, arc.end - twoPi, 1.0, 0.0, arc.cosineEnd, arc.sineEnd});
            }
            else
            {
              work.covered.push_back(arc);
            }
          }
          sortByBeginning(work.covered);

          // What is left of the latitude, arc by arc: each gap is one connected stretch of exposed surface.
          double cursor = 0.0;
          double cosineCursor = 1.0;
          double sineCursor = 0.0;
          for (std::size_t arc = 0; arc <= work.covered.size(); ++arc)
          {
            const bool last = (arc == work.covered.size());
            const double gapEnd = last ? twoPi : work.covered[arc].begin;

            if (gapEnd - cursor > sweepGapTolerance)
            {
              gap.begin = cursor;
              gap.end = gapEnd;
              gap.span = gapEnd - cursor;
              gap.cosineBegin = cosineCursor;
              gap.sineBegin = sineCursor;
              gap.cosineEnd = last ? 1.0 : work.covered[arc].cosineBegin;
              gap.sineEnd = last ? 0.0 : work.covered[arc].sineBegin;
              measure(gap);
            }

            if (!last && work.covered[arc].end > cursor)
            {
              cursor = work.covered[arc].end;
              cosineCursor = work.covered[arc].cosineEnd;
              sineCursor = work.covered[arc].sineEnd;
            }
          }
        }
      }
    }
  }
}
