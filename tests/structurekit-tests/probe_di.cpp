#include <gtest/gtest.h>
import std;
import double3;
import unit_cell;
import randomnumbers;
import pore_accessibility;
import apollonius_accessibility;
import voronoi_network;
import pore_diameters;

// A re-run of the ascent from one node with the rules varied: search reach, candidate count, jump bound.
struct Variant { double reachExtra; std::size_t considered; int jumpRule; std::size_t iterations; };
double ascend(const PoreAccessibility& access, const UnitCell& box, const std::vector<double3>& f, const std::vector<double>& r,
              double3 start, double startRadius, const Variant& v)
{
  double rmax = 0.0; for (double x : r) rmax = std::max(rmax, x);
  double maximal = startRadius; double3 point = start;
  struct C { double3 centre; double radius; double clearance; };
  for (std::size_t it = 0; it < v.iterations; ++it)
  {
    std::vector<C> cand;
    for (std::size_t j = 0; j < f.size(); ++j)
      for (int ix = -1; ix <= 1; ++ix) for (int iy = -1; iy <= 1; ++iy) for (int iz = -1; iz <= 1; ++iz)
      {
        double3 image = box.cell * (f[j] + double3(ix, iy, iz));
        double dist = (point - image).length();
        if (dist <= maximal + rmax + v.reachExtra) cand.push_back({image, r[j], dist - r[j]});
      }
    if (cand.size() < 4) return maximal;
    std::size_t considered = std::min(v.considered, cand.size());
    std::partial_sort(cand.begin(), cand.begin() + static_cast<std::ptrdiff_t>(considered), cand.end(), [](const C& a, const C& b) { return a.clearance < b.clearance; });
    double best = maximal; double3 bestPoint = point;
    for (std::size_t i0 = 0; i0 < considered; ++i0) for (std::size_t i1 = i0 + 1; i1 < considered; ++i1) for (std::size_t i2 = i1 + 1; i2 < considered; ++i2) for (std::size_t i3 = i2 + 1; i3 < considered; ++i3)
    {
      for (const ApolloniusSphere& sphere : apolloniusTangentSpheres({cand[i0].centre, cand[i1].centre, cand[i2].centre, cand[i3].centre}, {cand[i0].radius, cand[i1].radius, cand[i2].radius, cand[i3].radius}))
      {
        if (sphere.radius <= best) continue;
        double jump = (sphere.centre - point).length();
        if (v.jumpRule == 0 && jump > startRadius) continue;
        if (v.jumpRule == 1)
        {
          // the segment from point to the centre must be free
          bool free = true;
          for (int k = 1; k <= 20 && free; ++k) if (access.clearance(point + (sphere.centre - point) * (k / 20.0)) <= 0.0) free = false;
          if (!free) continue;
        }
        double3 wrapped = box.cell * double3::fract(box.inverseCell * sphere.centre);
        double clearance = access.clearance(wrapped);
        if (clearance > best) { best = clearance; bestPoint = wrapped; }
      }
    }
    if (best <= maximal + 1.0e-12) return maximal;
    maximal = best; point = bestPoint;
  }
  return maximal;
}

TEST(probe_di, ascent_failures)
{
  std::vector<Variant> variants = {{0.0, 8, 1, 16}, {2.5, 8, 1, 16}, {1.0e9, 8, 1, 16}, {0.0, 12, 1, 16}, {2.5, 12, 1, 16}, {1.0e9, 12, 1, 16},
                                   {0.0, 8, 1, 64}, {2.5, 8, 1, 64}, {1.0e9, 8, 1, 64}, {0.0, 12, 1, 64}, {2.5, 12, 1, 64}, {1.0e9, 12, 1, 64}, {5.0, 10, 1, 32}};
  std::size_t noNodeReachedByPositive = 0, noNodeReachedByBuried = 0, buriedWithFreeSeed = 0, buriedTotal = 0;
  std::vector<std::size_t> fixedBy(variants.size(), 0);
  std::vector<double3> F; std::vector<double> R;
  std::size_t total = 0, peaks = 0, failures = 0, unreachable = 0, matched = 0;
  double worstShortfall = 0.0, worstUnreachable = 0.0;
  RandomNumber directions{std::optional<std::size_t>(7)};
  for (int family = 0; family < 3; ++family)
  for (std::size_t seed = 1; seed <= 40; ++seed)
  {
    RandomNumber random{std::optional<std::size_t>(seed + 100 * family)};
    const double edge = (family == 2) ? 9.0 : 12.0;
    const std::size_t count = (family == 0) ? 24 : (family == 1 ? 40 : 14);
    std::vector<double3>& f = F; std::vector<double>& r = R; f.clear(); r.clear();
    for (std::size_t i = 0; i < count; ++i)
    {
      f.push_back(double3(random.uniform(), random.uniform(), random.uniform()));
      double u = random.uniform();
      r.push_back(family == 1 ? (u < 0.5 ? 0.5 : 2.5) : 0.5 + 2.0 * u);
    }
    UnitCell box(edge, edge, edge);
    PoreAccessibility radical, apollonius;
    try { radical = PoreAccessibility::create(box, f, r, 0.0); apollonius = ApolloniusAccessibility::create(box, f, r, 0.0).accessibility; } catch (std::exception&) { continue; }
    ++total;
    PoreDiameters dr = PoreDiameters::compute(radical.network);
    PoreDiameters da = PoreDiameters::compute(apollonius.network);
    bool difDiffers = da.includedAlongFreePathDiameter - dr.includedAlongFreePathDiameter > 1.0e-6;
    if (difDiffers)
    {
      // The Apollonius peak that sets Dif, and what the radical network makes of it.
      std::vector<std::size_t> all(apollonius.network.nodes.size()); std::iota(all.begin(), all.end(), 0uz);
      PercolatingPath pa = widestPercolatingPath(apollonius.network, all);
      std::size_t best = pa.componentNodes.front();
      for (std::size_t n : pa.componentNodes) if (apollonius.network.nodes[n].maximalRadius > apollonius.network.nodes[best].maximalRadius) best = n;
      const VoronoiNode& peak = apollonius.network.nodes[best];
      std::vector<std::size_t> allR(radical.network.nodes.size()); std::iota(allR.begin(), allR.end(), 0uz);
      PercolatingPath pr = widestPercolatingPath(radical.network, allR);
      bool matchedAnywhere = false, matchedInComponent = false; std::size_t insideCount = 0; double bestInsideR = -1.0;
      for (std::size_t n = 0; n < radical.network.nodes.size(); ++n)
      {
        const VoronoiNode& node = radical.network.nodes[n];
        double3 dm = box.applyPeriodicBoundaryConditions(node.maximalPosition - peak.position);
        bool m = dm.length() < 1.0e-4 && std::abs(node.maximalRadius - peak.radius) < 1.0e-6;
        if (m) { matchedAnywhere = true; if (std::ranges::find(pr.componentNodes, n) != pr.componentNodes.end()) matchedInComponent = true; }
        double3 dp = box.applyPeriodicBoundaryConditions(node.position - peak.position);
        if (dp.length() < peak.radius) { ++insideCount; bestInsideR = std::max(bestInsideR, node.maximalRadius); if (node.radius <= 0.0) std::print("   buried node inside Dif peak ball, radius {:.4f}\n", node.radius); }
      }
      std::print("DIF family {} seed {} apollonius Dif peak r {:.4f}; radical Dif {:.4f} Df {:.4f} vs apollonius Df {:.4f}; matched anywhere {} in radical component {}; radical nodes inside its ball {} best maximal {:.4f}; component sizes radical {} apollonius {}\n",
                 family, seed, peak.radius, 0.5 * dr.includedAlongFreePathDiameter, 0.5 * dr.freeSphereDiameter, 0.5 * da.freeSphereDiameter, matchedAnywhere, matchedInComponent, insideCount, bestInsideR, pr.componentNodes.size(), pa.componentNodes.size());
    }

    for (const VoronoiNode& node : radical.network.nodes)
    {
      if (node.radius > 0.0 || node.atomIndices.size() < 4) continue;
      ++buriedTotal;
      std::array<double3, 4> c; std::array<double, 4> rr;
      for (std::size_t a = 0; a < 4; ++a) { std::size_t j = node.atomIndices[a]; double3 image = box.cell * f[j]; c[a] = node.position + box.applyPeriodicBoundaryConditions(image - node.position); rr[a] = r[j]; }
      for (const ApolloniusSphere& sphere : apolloniusTangentSpheres(c, rr))
      {
        double3 wrapped = box.cell * double3::fract(box.inverseCell * sphere.centre);
        if (radical.clearance(wrapped) > 0.0) { ++buriedWithFreeSeed; break; }
      }
    }
    for (const VoronoiNode& peak : apollonius.network.nodes)
    {
      if (peak.radius <= 1.0e-6) continue;
      // Local maximum of the clearance? Every small step must lower it.
      bool isMax = true;
      const double here = radical.clearance(peak.position);
      for (int k = 0; k < 60 && isMax; ++k)
      {
        double3 d(directions.uniform() - 0.5, directions.uniform() - 0.5, directions.uniform() - 0.5);
        d = d * (1.0e-3 / d.length());
        if (radical.clearance(peak.position + d) > here + 1.0e-12) isMax = false;
      }
      if (!isMax) continue;
      ++peaks;
      bool anyMatch = false; bool anyInside = false; double bestInside = -1.0; std::vector<const VoronoiNode*> inside;
      for (const VoronoiNode& node : radical.network.nodes)
      {
        double3 dm = box.applyPeriodicBoundaryConditions(node.maximalPosition - peak.position);
        if (dm.length() < 1.0e-4 && std::abs(node.maximalRadius - peak.radius) < 1.0e-6) anyMatch = true;
        if (node.radius <= 0.0) continue;
        double3 dp = box.applyPeriodicBoundaryConditions(node.position - peak.position);
        if (dp.length() < peak.radius) { anyInside = true; bestInside = std::max(bestInside, node.maximalRadius); inside.push_back(&node); }
      }
      if (anyMatch) { ++matched; continue; }
      if (anyInside && bestInside >= peak.radius - 1.0e-6) { ++matched; continue; }
      if (anyInside)
      {
        ++failures; worstShortfall = std::max(worstShortfall, peak.radius - bestInside);
        for (std::size_t vi = 0; vi < variants.size(); ++vi)
        {
          double reached = -1.0;
          for (const VoronoiNode* node : inside) reached = std::max(reached, ascend(radical, box, f, r, node->position, node->radius, variants[vi]));
          if (reached >= peak.radius - 1.0e-6) ++fixedBy[vi];
        }
        if (false) std::print("ASCENT SHORT family {} seed {} peak r {:.4f} at ({:.3f},{:.3f},{:.3f}) best inside {:.4f} difDiffers {}\n", family, seed, peak.radius, peak.position.x, peak.position.y, peak.position.z, bestInside, difDiffers);
      }
      else
      {
        ++unreachable; worstUnreachable = std::max(worstUnreachable, peak.radius);
        const Variant full{1.0e9, 12, 1, 64};
        bool byPositive = false, byBuried = false;
        for (const VoronoiNode& node : radical.network.nodes)
        {
          if (node.radius > 0.0)
          {
            double reached = ascend(radical, box, f, r, node.position, node.radius, full);
            // reached the same peak? compare radius only (a cheap proxy) and position via a final check
            if (std::abs(reached - peak.radius) < 1.0e-6) byPositive = true;
          }
          else
          {
            // seed from the tangent sphere of the vertex's own atoms, at their images nearest the vertex
            if (node.atomIndices.size() < 4) continue;
            std::array<double3, 4> c; std::array<double, 4> rr;
            for (std::size_t a = 0; a < 4; ++a)
            {
              std::size_t j = node.atomIndices[a];
              double3 image = box.cell * f[j];
              image = node.position + box.applyPeriodicBoundaryConditions(image - node.position);
              c[a] = image; rr[a] = r[j];
            }
            for (const ApolloniusSphere& sphere : apolloniusTangentSpheres(c, rr))
            {
              double3 wrapped = box.cell * double3::fract(box.inverseCell * sphere.centre);
              double clearance = radical.clearance(wrapped);
              if (clearance <= 0.0) continue;
              double reached = ascend(radical, box, f, r, wrapped, clearance, full);
              if (std::abs(reached - peak.radius) < 1.0e-6) byBuried = true;
            }
          }
          if (byPositive) break;
        }
        if (byPositive) ++noNodeReachedByPositive; else if (byBuried) ++noNodeReachedByBuried;
        if (false) std::print("NO NODE family {} seed {} peak r {:.4f} difDiffers {}\n", family, seed, peak.radius, difDiffers);
      }
    }
  }
  for (std::size_t vi = 0; vi < variants.size(); ++vi) std::print("variant {} (reachExtra {:.0e} considered {} jumpRule {} iterations {}): fixed {} of {}\n", vi, variants[vi].reachExtra, variants[vi].considered, variants[vi].jumpRule, variants[vi].iterations, fixedBy[vi], failures);
  std::print("no-node peaks reached from a positive node under the full variant {}, otherwise from a buried vertex's tangent sphere {}; buried nodes {} of which with a free tangent seed {}\n", noNodeReachedByPositive, noNodeReachedByBuried, buriedTotal, buriedWithFreeSeed);
  std::print("structures {} true peaks {} matched {} ascentShort {} worstShortfall {:.4f} noNodeInside {} worstUnreachable {:.4f}\n", total, peaks, matched, failures, worstShortfall, unreachable, worstUnreachable);
}
