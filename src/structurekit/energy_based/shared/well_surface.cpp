module;

module energy_shared_well_surface;

import std;

import uint3;
import double3;
import unit_cell;
import crystal;
import pair_interactions;
import blocking_spheres;
import energy_shared_linear_probe;
import energy_shared_electrostatic_potential_grid;
import energy_shared_energy_backend;
import energy_shared_well_field;
import energy_shared_isosurface;
import energy_isosurface;
import energy_well_field;

// The measurement itself, in the units the field is held in. Nothing here knows what a Kelvin is: the report
// and the driver next door do that, and the reason they are a separate translation unit is written up there.

namespace
{
std::vector<double3> trianglesOf(const EnergyBackend *backend, std::span<const float> field, uint3 gridSize,
                                double isoValue)
{
  if (backend != nullptr) return backend->isosurfaceTriangles(field, gridSize, isoValue);
  return EnergyIsosurface::trianglesOfIsosurface(field, gridSize, isoValue);
}

// Connected components of the filament mesh by periodically welded vertices; components below
// wellFilamentMinimumArea are dropped.
std::vector<double3> removeFilamentSpecks(const UnitCell &unitCell, std::vector<double3> corners)
{
  const std::size_t triangles = corners.size() / 3;
  if (triangles == 0) return {};

  std::vector<std::size_t> parent(triangles);
  std::iota(parent.begin(), parent.end(), 0);

  auto findRoot = [&](std::size_t i)
  {
    std::size_t r = i;
    while (parent[r] != r) r = parent[r];
    std::size_t w = i;
    while (parent[w] != r)
    {
      std::size_t nx = parent[w];
      parent[w] = r;
      w = nx;
    }
    return r;
  };

  struct Key
  {
    std::int32_t x;
    std::int32_t y;
    std::int32_t z;
    bool operator==(const Key &) const = default;
  };
  struct KeyHash
  {
    std::size_t operator()(const Key &key) const
    {
      std::size_t h = static_cast<std::size_t>(key.x);
      h = h * 6364136223846793005ull + static_cast<std::size_t>(key.y);
      h = h * 6364136223846793005ull + static_cast<std::size_t>(key.z);
      return h;
    }
  };

  auto quantize = [](double x) -> std::int32_t
  {
    double wrapped = x - std::floor(x);
    std::int32_t r = static_cast<std::int32_t>(std::rint(wrapped * 1048576.0));
    return r == 1048576 ? 0 : r;
  };

  std::unordered_map<Key, std::size_t, KeyHash> seen;
  seen.reserve(3 * triangles);
  for (std::size_t t = 0; t < triangles; ++t)
  {
    for (std::size_t v = 0; v < 3; ++v)
    {
      const double3 &p = corners[3 * t + v];
      Key key{quantize(p.x), quantize(p.y), quantize(p.z)};
      auto [it, inserted] = seen.emplace(key, t);
      if (!inserted)
      {
        std::size_t a = findRoot(t);
        std::size_t b = findRoot(it->second);
        if (a != b) parent[std::max(a, b)] = std::min(a, b);
      }
    }
  }

  std::unordered_map<std::size_t, double> area;
  for (std::size_t t = 0; t < triangles; ++t)
  {
    double3 p0 = unitCell.cell * corners[3 * t];
    double3 p1 = unitCell.cell * corners[3 * t + 1];
    double3 p2 = unitCell.cell * corners[3 * t + 2];
    area[findRoot(t)] += 0.5 * double3::cross(p1 - p0, p2 - p0).length();
  }

  std::vector<double3> kept;
  kept.reserve(corners.size());
  for (std::size_t t = 0; t < triangles; ++t)
  {
    if (area[findRoot(t)] >= wellFilamentMinimumArea)
    {
      kept.push_back(corners[3 * t]);
      kept.push_back(corners[3 * t + 1]);
      kept.push_back(corners[3 * t + 2]);
    }
  }
  return kept;
}

void refineThrough(std::vector<double3> &corners, std::vector<double> &energies, const EnergyBackend *backend,
                   const PairInteractions &interactions, const Crystal &framework, const WellField &field,
                   std::span<const BlockingSphere> blockingSpheres, const ElectrostaticPotentialGrid *potential,
                   double coulombFactor, double iso)
{
  NeighbourhoodParameters parameters;
  parameters.probe = field.probe;
  parameters.numberOfOrientations = field.numberOfOrientations;
  parameters.thermalEnergy = field.thermalEnergy;
  parameters.extraReach = wellRefinementReach;
  parameters.ceiling = field.ceiling;
  parameters.potential = potential;
  parameters.coulombFactor = coulombFactor;

  if (backend != nullptr && backend->refineWellVertices)
  {
    backend->refineWellVertices(corners, energies, interactions, framework, parameters, blockingSpheres, iso);
    return;
  }
  WellFieldCPU::refineVertices(corners, energies, interactions, framework, parameters, blockingSpheres, iso);
}

double3 principalAxis(double xx, double yy, double zz, double xy, double xz, double yz, double lambda)
{
  double3 r1(xx - lambda, xy, xz);
  double3 r2(xy, yy - lambda, yz);
  double3 r3(xz, yz, zz - lambda);
  double3 v = double3::cross(r1, r2);
  if (v.length_squared() < 1.0e-24) v = double3::cross(r1, r3);
  if (v.length_squared() < 1.0e-24) v = double3::cross(r2, r3);
  if (v.length_squared() < 1.0e-24) return double3(1.0, 0.0, 0.0);
  return double3::normalize(v);
}

// Local dimensionality of the merged-well set: a 1-D file along a channel axis, or a 2-D midplane where two
// walls have closed the sheet. The ridge is the Gelb-Gubbins inscribed radius of {U < iso} --- the first
// half of the energy PSD, the same distance-to-iso the pore-size curve is built from --- not the distance
// transform of the reliability overlay. That overlay is a noisy subset of the energy void: holes in rel
// put every voxel next to a boundary, the hop-count ridge is the whole blob, and a 1-D graph of it is
// tens of parallel files. A 1-D core that still has a contact sheet next door is the interior of a lined
// pore and is dropped.
//
// Packing is measured on that classified set, not on a cylinder of the whole overlay. A wrapping tube is a
// 1-D file even when the overlay is a flattened ribbon; a slab that wraps two axes with one component per
// slice is a 2-D midplane. File length is whole channels times whole sites along the period.
void assignFilamentPacking(WellSurface &result, const WellField &field, const UnitCell &unitCell,
                           std::span<const std::size_t> filamentIndex, std::span<const double3> filamentCorners)
{
  const std::size_t nFilament = result.filamentVoxels.size();
  if (nFilament == 0 || filamentIndex.size() != nFilament) return;

  const int nx = static_cast<int>(field.gridSize.x);
  const int ny = static_cast<int>(field.gridSize.y);
  const int nz = static_cast<int>(field.gridSize.z);
  const std::size_t nVoxels = field.numberOfVoxels();
  auto wrap = [](int v, int n) -> int
  {
    int r = v % n;
    return r < 0 ? r + n : r;
  };
  auto linear = [&](int x, int y, int z) -> std::size_t
  { return (static_cast<std::size_t>(z) * static_cast<std::size_t>(ny) + static_cast<std::size_t>(y)) *
           static_cast<std::size_t>(nx) + static_cast<std::size_t>(x); };
  auto decode = [&](std::size_t i, int &x, int &y, int &z)
  {
    x = static_cast<int>(i % static_cast<std::size_t>(nx));
    y = static_cast<int>((i / static_cast<std::size_t>(nx)) % static_cast<std::size_t>(ny));
    z = static_cast<int>(i / (static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny)));
  };

  std::vector<std::uint32_t> slotOf(nVoxels, std::numeric_limits<std::uint32_t>::max());
  std::vector<double3> cartesian(nFilament);
  for (std::size_t s = 0; s < nFilament; ++s)
  {
    slotOf[filamentIndex[s]] = static_cast<std::uint32_t>(s);
    int ix, iy, iz;
    decode(filamentIndex[s], ix, iy, iz);
    double3 frac((static_cast<double>(ix) + 0.5) / static_cast<double>(nx),
                 (static_cast<double>(iy) + 0.5) / static_cast<double>(ny),
                 (static_cast<double>(iz) + 0.5) / static_cast<double>(nz));
    cartesian[s] = unitCell.cell * frac;
  }

  const int d6[6][3] = {{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}};
  const float step[6] = {
      static_cast<float>((unitCell.cell * double3(1.0 / static_cast<double>(nx), 0.0, 0.0)).length()),
      static_cast<float>((unitCell.cell * double3(1.0 / static_cast<double>(nx), 0.0, 0.0)).length()),
      static_cast<float>((unitCell.cell * double3(0.0, 1.0 / static_cast<double>(ny), 0.0)).length()),
      static_cast<float>((unitCell.cell * double3(0.0, 1.0 / static_cast<double>(ny), 0.0)).length()),
      static_cast<float>((unitCell.cell * double3(0.0, 0.0, 1.0 / static_cast<double>(nz))).length()),
      static_cast<float>((unitCell.cell * double3(0.0, 0.0, 1.0 / static_cast<double>(nz))).length())};

  const float iso = static_cast<float>(result.isoValue);
  std::vector<std::uint8_t> inVoid(nVoxels, 0);
  for (std::size_t i = 0; i < nVoxels; ++i)
  {
    if (field.energy[i] < iso) inVoid[i] = 1;
  }

  // Anisotropic inscribed radius of {U < iso}, in Å: Gelb-Gubbins distance to the energy iso.
  constexpr float unreached = std::numeric_limits<float>::max();
  std::vector<float> inscribed(nVoxels, unreached);
  using Node = std::pair<float, std::size_t>;
  std::priority_queue<Node, std::vector<Node>, std::greater<Node>> closer;
  for (std::size_t i = 0; i < nVoxels; ++i)
  {
    if (!inVoid[i]) continue;
    int ix, iy, iz;
    decode(i, ix, iy, iz);
    float seed = unreached;
    for (int n = 0; n < 6; ++n)
    {
      std::size_t nb = linear(wrap(ix + d6[n][0], nx), wrap(iy + d6[n][1], ny), wrap(iz + d6[n][2], nz));
      if (inVoid[nb]) continue;
      seed = std::min(seed, 0.5f * step[n]);
    }
    if (seed < unreached)
    {
      inscribed[i] = seed;
      closer.push({seed, i});
    }
  }
  while (!closer.empty())
  {
    auto [d, i] = closer.top();
    closer.pop();
    if (d > inscribed[i]) continue;
    int ix, iy, iz;
    decode(i, ix, iy, iz);
    for (int n = 0; n < 6; ++n)
    {
      std::size_t nb = linear(wrap(ix + d6[n][0], nx), wrap(iy + d6[n][1], ny), wrap(iz + d6[n][2], nz));
      if (!inVoid[nb]) continue;
      float next = d + step[n];
      if (next < inscribed[nb])
      {
        inscribed[nb] = next;
        closer.push({next, nb});
      }
    }
  }

  std::vector<std::uint8_t> ridge(nFilament, 0);
  for (std::size_t s = 0; s < nFilament; ++s)
  {
    const float r = inscribed[filamentIndex[s]];
    bool localMax = r > 0.0f && r < unreached;
    if (!localMax)
    {
      ridge[s] = 0;
      continue;
    }
    int ix, iy, iz;
    decode(filamentIndex[s], ix, iy, iz);
    for (int n = 0; n < 6; ++n)
    {
      std::size_t nb = linear(wrap(ix + d6[n][0], nx), wrap(iy + d6[n][1], ny), wrap(iz + d6[n][2], nz));
      if (inscribed[nb] >= unreached) continue;
      if (inscribed[nb] > r)
      {
        localMax = false;
        break;
      }
    }
    ridge[s] = localMax ? 1 : 0;
  }

  auto gyrationEigenvalues = [](double xx, double yy, double zz, double xy, double xz, double yz, double &l1,
                                double &l2, double &l3)
  {
    double q = (xx + yy + zz) / 3.0;
    double xxq = xx - q, yyq = yy - q, zzq = zz - q;
    double p2 = xxq * xxq + yyq * yyq + zzq * zzq + 2.0 * (xy * xy + xz * xz + yz * yz);
    if (!(p2 > 1.0e-24))
    {
      l1 = l2 = l3 = q;
      return;
    }
    double p = std::sqrt(p2 / 6.0);
    double invp = 1.0 / p;
    double bxx = xxq * invp, byy = yyq * invp, bzz = zzq * invp;
    double bxy = xy * invp, bxz = xz * invp, byz = yz * invp;
    double det = bxx * (byy * bzz - byz * byz) - bxy * (bxy * bzz - byz * bxz) + bxz * (bxy * byz - byy * bxz);
    double r = std::min(1.0, std::max(-1.0, det * 0.5));
    double phi = std::acos(r) / 3.0;
    l1 = q + 2.0 * p * std::cos(phi);
    l3 = q + 2.0 * p * std::cos(phi + 2.0 * std::numbers::pi / 3.0);
    l2 = 3.0 * q - l1 - l3;
    if (l2 < l3) std::swap(l2, l3);
    if (l1 < l2) std::swap(l1, l2);
    if (l2 < l3) std::swap(l2, l3);
  };

  const double minSpacing = std::min({unitCell.lengthA / static_cast<double>(nx),
                                      unitCell.lengthB / static_cast<double>(ny),
                                      unitCell.lengthC / static_cast<double>(nz)});
  constexpr double neighbourhood = 4.0;
  constexpr double neighbourhoodSquared = neighbourhood * neighbourhood;
  const int lim = std::min(8, std::max(1, static_cast<int>(std::ceil(neighbourhood / std::max(minSpacing, 1.0e-6)))));
  std::vector<std::uint8_t> ridgePlanar(nFilament, 0);
  std::vector<double3> tangent(nFilament);
  std::vector<double3> planeNormal(nFilament);
  for (std::size_t s = 0; s < nFilament; ++s)
  {
    if (!ridge[s]) continue;
    int ix, iy, iz;
    decode(filamentIndex[s], ix, iy, iz);
    std::vector<double3> local;
    local.reserve(64);
    for (int dz = -lim; dz <= lim; ++dz)
    {
      for (int dy = -lim; dy <= lim; ++dy)
      {
        for (int dx = -lim; dx <= lim; ++dx)
        {
          std::size_t nb = linear(wrap(ix + dx, nx), wrap(iy + dy, ny), wrap(iz + dz, nz));
          std::uint32_t other = slotOf[nb];
          if (other == std::numeric_limits<std::uint32_t>::max() || !ridge[other]) continue;
          double3 dr = unitCell.applyPeriodicBoundaryConditions(cartesian[other] - cartesian[s]);
          if (double3::dot(dr, dr) > neighbourhoodSquared) continue;
          local.push_back(dr);
        }
      }
    }
    if (local.size() < 6) continue;
    double inv = 1.0 / static_cast<double>(local.size());
    double meanx = 0.0, meany = 0.0, meanz = 0.0;
    for (const double3 &p : local)
    {
      meanx += p.x;
      meany += p.y;
      meanz += p.z;
    }
    meanx *= inv;
    meany *= inv;
    meanz *= inv;
    double xx = 0.0, yy = 0.0, zz = 0.0, xy = 0.0, xz = 0.0, yz = 0.0;
    for (const double3 &p : local)
    {
      double cx = p.x - meanx, cy = p.y - meany, cz = p.z - meanz;
      xx += cx * cx;
      yy += cy * cy;
      zz += cz * cz;
      xy += cx * cy;
      xz += cx * cz;
      yz += cy * cz;
    }
    xx *= inv;
    yy *= inv;
    zz *= inv;
    xy *= inv;
    xz *= inv;
    yz *= inv;
    double l1, l2, l3;
    gyrationEigenvalues(xx, yy, zz, xy, xz, yz, l1, l2, l3);
    tangent[s] = principalAxis(xx, yy, zz, xy, xz, yz, l1);
    planeNormal[s] = principalAxis(xx, yy, zz, xy, xz, yz, l3);
    if (l1 > 0.0 && l2 > 0.45 * l1) ridgePlanar[s] = 1;
  }

  std::vector<std::uint8_t> sheetAdjacent(nFilament, 0);
  std::vector<std::uint8_t> nearPlanar(nFilament, 0);
  for (std::size_t s = 0; s < nFilament; ++s)
  {
    int ix, iy, iz;
    decode(filamentIndex[s], ix, iy, iz);
    for (const int *d : d6)
    {
      std::size_t nb = linear(wrap(ix + d[0], nx), wrap(iy + d[1], ny), wrap(iz + d[2], nz));
      if (field.distance[nb] > 0.0f) sheetAdjacent[s] = 1;
    }
    for (int dz = -1; dz <= 1 && !nearPlanar[s]; ++dz)
    {
      for (int dy = -1; dy <= 1 && !nearPlanar[s]; ++dy)
      {
        for (int dx = -1; dx <= 1; ++dx)
        {
          std::size_t nb = linear(wrap(ix + dx, nx), wrap(iy + dy, ny), wrap(iz + dz, nz));
          std::uint32_t other = slotOf[nb];
          if (other == std::numeric_limits<std::uint32_t>::max()) continue;
          if (ridgePlanar[other])
          {
            nearPlanar[s] = 1;
            break;
          }
        }
      }
    }
  }

  // A 1-D file is a wrapping channel, not a graph of overlay voxels. Wrapping is read from the
  // medial tight-channel set (d ≤ 0 and low reliability), not from {U < iso}: a barely-fitting
  // probe pinches the energy well into pockets, but the 8-ring still wraps and GCMC still sits
  // one molecule per file. Components are then classified by how that channel wraps.
  std::vector<std::uint8_t> inChannel(nVoxels, 0);
  for (std::size_t i = 0; i < nVoxels; ++i)
  {
    if (field.distance[i] <= 0.0f &&
        field.reliability[i] <= static_cast<float>(wellFilamentReliabilityThreshold))
      inChannel[i] = 1;
  }

  std::vector<int> channelComp(nVoxels, -1);
  struct ComponentWrap
  {
    bool a{false};
    bool b{false};
    bool c{false};
  };
  std::vector<ComponentWrap> componentWraps;
  std::vector<int> imageX(nVoxels, 0), imageY(nVoxels, 0), imageZ(nVoxels, 0);
  for (std::size_t seed = 0; seed < nVoxels; ++seed)
  {
    if (!inChannel[seed] || channelComp[seed] >= 0) continue;
    const int id = static_cast<int>(componentWraps.size());
    ComponentWrap wrapFlags;
    std::queue<std::size_t> walk;
    channelComp[seed] = id;
    imageX[seed] = imageY[seed] = imageZ[seed] = 0;
    walk.push(seed);
    while (!walk.empty())
    {
      const std::size_t i = walk.front();
      walk.pop();
      int ix, iy, iz;
      decode(i, ix, iy, iz);
      for (int dz = -1; dz <= 1; ++dz)
      {
        for (int dy = -1; dy <= 1; ++dy)
        {
          for (int dx = -1; dx <= 1; ++dx)
          {
            if (dx == 0 && dy == 0 && dz == 0) continue;
            const int jx = ix + dx, jy = iy + dy, jz = iz + dz;
            const std::size_t j = linear(wrap(jx, nx), wrap(jy, ny), wrap(jz, nz));
            if (!inChannel[j]) continue;
            const int nxImage = imageX[i] + (jx >= nx ? 1 : jx < 0 ? -1 : 0);
            const int nyImage = imageY[i] + (jy >= ny ? 1 : jy < 0 ? -1 : 0);
            const int nzImage = imageZ[i] + (jz >= nz ? 1 : jz < 0 ? -1 : 0);
            if (channelComp[j] >= 0)
            {
              if (channelComp[j] != id) continue;
              if (nxImage != imageX[j]) wrapFlags.a = true;
              if (nyImage != imageY[j]) wrapFlags.b = true;
              if (nzImage != imageZ[j]) wrapFlags.c = true;
              continue;
            }
            channelComp[j] = id;
            imageX[j] = nxImage;
            imageY[j] = nyImage;
            imageZ[j] = nzImage;
            walk.push(j);
          }
        }
      }
    }
    componentWraps.push_back(wrapFlags);
  }

  std::vector<int> componentOf(nFilament, -1);
  for (std::size_t s = 0; s < nFilament; ++s) componentOf[s] = channelComp[filamentIndex[s]];

  std::vector<std::vector<std::size_t>> channelMembers(componentWraps.size());
  for (std::size_t i = 0; i < nVoxels; ++i)
  {
    if (channelComp[i] >= 0) channelMembers[static_cast<std::size_t>(channelComp[i])].push_back(i);
  }

  auto filesAlong = [&](int axis, int nPlane, int nU, int nV, int component) -> double
  {
    if (nPlane <= 0 || component < 0) return 0.0;
    if (static_cast<std::size_t>(component) >= channelMembers.size()) return 0.0;
    double sum = 0.0;
    int nUsed = 0;
    std::vector<std::uint8_t> plane(static_cast<std::size_t>(nU) * static_cast<std::size_t>(nV), 0);
    std::vector<int> parent(static_cast<std::size_t>(nU) * static_cast<std::size_t>(nV), -1);
    auto at = [&](int u, int v) -> std::size_t
    { return static_cast<std::size_t>(v) * static_cast<std::size_t>(nU) + static_cast<std::size_t>(u); };
    auto find = [&](int i) -> int
    {
      int r = i;
      while (parent[static_cast<std::size_t>(r)] != r) r = parent[static_cast<std::size_t>(r)];
      while (parent[static_cast<std::size_t>(i)] != r)
      {
        int nxt = parent[static_cast<std::size_t>(i)];
        parent[static_cast<std::size_t>(i)] = r;
        i = nxt;
      }
      return r;
    };
    for (int p = 0; p < nPlane; ++p)
    {
      std::fill(plane.begin(), plane.end(), 0);
      std::fill(parent.begin(), parent.end(), -1);
      for (std::size_t idx : channelMembers[static_cast<std::size_t>(component)])
      {
        int ix, iy, iz;
        decode(idx, ix, iy, iz);
        int u = 0, v = 0, planeIndex = 0;
        if (axis == 0)
        {
          planeIndex = ix;
          u = iy;
          v = iz;
        }
        else if (axis == 1)
        {
          planeIndex = iy;
          u = ix;
          v = iz;
        }
        else
        {
          planeIndex = iz;
          u = ix;
          v = iy;
        }
        if (planeIndex != p) continue;
        plane[at(u, v)] = 1;
        parent[at(u, v)] = static_cast<int>(at(u, v));
      }
      bool any = false;
      for (int v = 0; v < nV; ++v)
      {
        for (int u = 0; u < nU; ++u)
        {
          if (!plane[at(u, v)]) continue;
          any = true;
          for (int dv = -1; dv <= 1; ++dv)
          {
            for (int du = -1; du <= 1; ++du)
            {
              if (du == 0 && dv == 0) continue;
              const int uu = wrap(u + du, nU), vv = wrap(v + dv, nV);
              if (!plane[at(uu, vv)]) continue;
              int a = find(static_cast<int>(at(u, v)));
              int b = find(static_cast<int>(at(uu, vv)));
              if (a != b) parent[static_cast<std::size_t>(std::max(a, b))] = std::min(a, b);
            }
          }
        }
      }
      if (!any) continue;
      int ncc = 0;
      for (int v = 0; v < nV; ++v)
      {
        for (int u = 0; u < nU; ++u)
        {
          if (!plane[at(u, v)]) continue;
          if (find(static_cast<int>(at(u, v))) == static_cast<int>(at(u, v))) ++ncc;
        }
      }
      sum += static_cast<double>(ncc);
      ++nUsed;
    }
    if (nUsed == 0) return 0.0;
    return sum / static_cast<double>(nUsed);
  };

  constexpr double fileSpacing = 3.864107;
  auto packedLength = [&](double nFilesMean, double period, double transverseA, double transverseB) -> double
  {
    if (!(period > 0.0) || !(nFilesMean > 0.0)) return 0.0;
    double nFiles = std::round(nFilesMean);
    if (nFiles < 1.0) nFiles = 1.0;
    const double maxFiles =
        std::max(1.0, std::round(std::min(transverseA, transverseB) / fileSpacing));
    nFiles = std::min(nFiles, maxFiles);
    const double sites = std::max(1.0, std::round(period / fileSpacing));
    return nFiles * sites * fileSpacing;
  };

  std::vector<std::uint8_t> kind(nFilament, 1);
  std::vector<std::uint8_t> slitComponent(componentWraps.size(), 0);
  double nFilesA = 0.0, nFilesB = 0.0, nFilesC = 0.0;
  const double a = unitCell.lengthA, b = unitCell.lengthB, cLen = unitCell.lengthC;
  for (std::size_t c = 0; c < componentWraps.size(); ++c)
  {
    const ComponentWrap w = componentWraps[c];
    const int nWrap = (w.a ? 1 : 0) + (w.b ? 1 : 0) + (w.c ? 1 : 0);
    const double nA = w.a ? filesAlong(0, nx, ny, nz, static_cast<int>(c)) : 0.0;
    const double nB = w.b ? filesAlong(1, ny, nx, nz, static_cast<int>(c)) : 0.0;
    const double nC = w.c ? filesAlong(2, nz, nx, ny, static_cast<int>(c)) : 0.0;
    const bool slit = nWrap >= 2 && (!w.a || nA < 1.5) && (!w.b || nB < 1.5) && (!w.c || nC < 1.5) &&
                      ((!w.a && a <= b && a <= cLen) || (!w.b && b <= a && b <= cLen) ||
                       (!w.c && cLen <= a && cLen <= b));
    slitComponent[c] = slit ? 1 : 0;
    if (slit) continue;
    bool hasWell = false;
    for (std::size_t s = 0; s < nFilament && !hasWell; ++s)
    {
      if (componentOf[s] == static_cast<int>(c)) hasWell = true;
    }
    if (!hasWell) continue;
    auto addFiles = [&](double n) -> double
    {
      double f = std::round(n);
      return f < 1.0 ? 1.0 : f;
    };
    if (nWrap == 3)
    {
      const double pa = w.a ? a : 1.0e300, pb = w.b ? b : 1.0e300, pc = w.c ? cLen : 1.0e300;
      if (pa <= pb && pa <= pc) nFilesA += addFiles(nA);
      else if (pb <= pc) nFilesB += addFiles(nB);
      else nFilesC += addFiles(nC);
    }
    else if (nWrap == 1)
    {
      if (w.a) nFilesA += addFiles(nA);
      if (w.b) nFilesB += addFiles(nB);
      if (w.c) nFilesC += addFiles(nC);
    }
    else
    {
      bool packed = false;
      if (w.a && nA >= 1.5)
      {
        nFilesA += addFiles(nA);
        packed = true;
      }
      if (w.b && nB >= 1.5)
      {
        nFilesB += addFiles(nB);
        packed = true;
      }
      if (w.c && nC >= 1.5)
      {
        nFilesC += addFiles(nC);
        packed = true;
      }
      if (!packed)
      {
        const double pa = w.a ? a : 1.0e300, pb = w.b ? b : 1.0e300, pc = w.c ? cLen : 1.0e300;
        if (pa <= pb && pa <= pc) nFilesA += 1.0;
        else if (pb <= pc) nFilesB += 1.0;
        else nFilesC += 1.0;
      }
    }
  }
  auto capFiles = [&](double n, double t0, double t1) -> double
  {
    if (!(n > 0.0)) return 0.0;
    const double maxFiles = std::max(1.0, std::round(std::min(t0, t1) / fileSpacing));
    return std::min(n, maxFiles);
  };
  nFilesA = capFiles(nFilesA, b, cLen);
  nFilesB = capFiles(nFilesB, a, cLen);
  nFilesC = capFiles(nFilesC, a, b);
  const double length1D = packedLength(nFilesA, a, b, cLen) + packedLength(nFilesB, b, a, cLen) +
                          packedLength(nFilesC, cLen, a, b);

  for (std::size_t s = 0; s < nFilament; ++s)
  {
    const int c = componentOf[s];
    if (c >= 0 && slitComponent[static_cast<std::size_t>(c)])
    {
      kind[s] = 2;
      continue;
    }
    const ComponentWrap w = (c >= 0) ? componentWraps[static_cast<std::size_t>(c)] : ComponentWrap{};
    const int nWrap = (w.a ? 1 : 0) + (w.b ? 1 : 0) + (w.c ? 1 : 0);
    if (nWrap >= 1)
    {
      kind[s] = 1;
      continue;
    }
    if (nearPlanar[s])
      kind[s] = 2;
    else if (sheetAdjacent[s])
      kind[s] = 0;
    else
      kind[s] = 1;
  }

  double volume1D = 0.0, volume2D = 0.0;
  for (std::size_t s = 0; s < nFilament; ++s)
  {
    if (kind[s] == 1) volume1D += result.filamentVoxels[s].volume;
    if (kind[s] == 2) volume2D += result.filamentVoxels[s].volume;
  }

  // The midplane is the planar ridge: each ridge voxel contributes the grid face perpendicular to the
  // plane normal. Half the overlay is the same surface measured as a mesh, used when the ridge is empty.
  const double voxelVolume = result.filamentVoxels.front().volume;
  const double3 ax = unitCell.cell * double3(1.0 / static_cast<double>(nx), 0.0, 0.0);
  const double3 ay = unitCell.cell * double3(0.0, 1.0 / static_cast<double>(ny), 0.0);
  const double3 az = unitCell.cell * double3(0.0, 0.0, 1.0 / static_cast<double>(nz));
  double ridgeFace = 0.0;
  for (std::size_t s = 0; s < nFilament; ++s)
  {
    if (!ridge[s] || kind[s] != 2) continue;
    double3 nrm = planeNormal[s];
    if (nrm.length_squared() < 0.25) nrm = double3(1.0, 0.0, 0.0);
    nrm = double3::normalize(nrm);
    double thickness = std::max({std::abs(double3::dot(nrm, ax)), std::abs(double3::dot(nrm, ay)),
                                 std::abs(double3::dot(nrm, az))});
    if (thickness > 0.0) ridgeFace += voxelVolume / thickness;
  }

  double overlayPlanar = 0.0;
  auto wrap01 = [](double x) { return x - std::floor(x); };
  for (std::size_t i = 0; i + 2 < filamentCorners.size(); i += 3)
  {
    double3 fc = (filamentCorners[i] + filamentCorners[i + 1] + filamentCorners[i + 2]) * (1.0 / 3.0);
    int cx = static_cast<int>(std::floor(wrap01(fc.x) * static_cast<double>(nx))) % nx;
    int cy = static_cast<int>(std::floor(wrap01(fc.y) * static_cast<double>(ny))) % ny;
    int cz = static_cast<int>(std::floor(wrap01(fc.z) * static_cast<double>(nz))) % nz;
    if (cx < 0) cx += nx;
    if (cy < 0) cy += ny;
    if (cz < 0) cz += nz;
    int votesPlanar = 0;
    int votesFile = 0;
    for (int dz = -2; dz <= 2; ++dz)
    {
      for (int dy = -2; dy <= 2; ++dy)
      {
        for (int dx = -2; dx <= 2; ++dx)
        {
          std::uint32_t slot =
              slotOf[linear(wrap(cx + dx, nx), wrap(cy + dy, ny), wrap(cz + dz, nz))];
          if (slot == std::numeric_limits<std::uint32_t>::max()) continue;
          if (kind[slot] == 2)
            ++votesPlanar;
          else if (kind[slot] == 1)
            ++votesFile;
        }
      }
    }
    if (votesPlanar > votesFile && votesPlanar > 0)
    {
      double3 p0 = unitCell.cell * filamentCorners[i];
      double3 p1 = unitCell.cell * filamentCorners[i + 1];
      double3 p2 = unitCell.cell * filamentCorners[i + 2];
      overlayPlanar += 0.5 * double3::cross(p1 - p0, p2 - p0).length();
    }
  }
  double area2D = ridgeFace;
  if (!(area2D > 0.0)) area2D = 0.5 * overlayPlanar;
  if (!(area2D > 0.0) && volume2D > 0.0) area2D = volume2D / std::cbrt(voxelVolume);

  result.filamentRidgeLength = length1D;
  result.filamentMedialArea = area2D;

  for (std::size_t s = 0; s < nFilament; ++s)
  {
    FilamentVoxel &voxel = result.filamentVoxels[s];
    voxel.length = 0.0;
    voxel.area = 0.0;
    if (kind[s] == 1 && volume1D > 0.0) voxel.length = voxel.volume / volume1D * length1D;
    if (kind[s] == 2 && volume2D > 0.0) voxel.area = voxel.volume / volume2D * area2D;
  }
}
}  // namespace


WellField computeWellField(const PairInteractions &interactions, const Crystal &framework, const LinearProbe &probe,
                           uint3 gridSize, std::size_t numberOfOrientations, double thermalEnergy,
                           std::span<const BlockingSphere> blockingSpheres, double blockedEnergyPerAngstrom,
                           double ceiling, const ElectrostaticPotentialGrid *potential, double coulombFactor)
{
  return WellFieldCPU::compute(interactions, framework, probe, gridSize, numberOfOrientations, thermalEnergy,
                               blockingSpheres, blockedEnergyPerAngstrom, ceiling, potential, coulombFactor);
}


float effectiveTrimIsovalue(const WellField &field, double isovalue)
{
  float iso = static_cast<float>(isovalue);
  double deepest = field.deepestEnergy();
  if (!(iso > static_cast<float>(deepest)))
  {
    iso = 0.25f * static_cast<float>(deepest);
  }
  return iso;
}


WellSurface wellSurfaceOfField(const Crystal &framework, const PairInteractions &interactions,
                               const WellField &field, double isoValue, double energyScale, double thermalEnergy,
                               std::span<const BlockingSphere> blockingSpheres, const EnergyBackend *backend,
                               const ElectrostaticPotentialGrid *potential, double coulombFactor)
{
  WellSurface result;
  result.isoValue = isoValue;
  result.energyScale = energyScale;
  result.thermalEnergy = thermalEnergy;

  if (field.numberOfVoxels() == 0) return result;

  std::chrono::steady_clock::time_point time_begin = std::chrono::steady_clock::now();

  const float iso = effectiveTrimIsovalue(field, isoValue);
  result.isoValue = static_cast<double>(iso);

  const std::size_t numberOfVoxels = field.numberOfVoxels();
  std::vector<float> combined(numberOfVoxels, 0.0f);
  for (std::size_t i = 0; i < numberOfVoxels; ++i)
  {
    combined[i] = std::max(-field.distance[i], static_cast<float>(energyScale) * (field.energy[i] - iso));
  }

  std::vector<double3> corners = trianglesOf(backend, combined, field.gridSize, 0.0);

  std::vector<double> energies;
  refineThrough(corners, energies, backend, interactions, framework, field, blockingSpheres, potential,
                coulombFactor, static_cast<double>(iso));
  for (double energy : energies)
  {
    if (energy > static_cast<double>(iso) - 0.02 * std::fabs(static_cast<double>(iso)))
    {
      ++result.numberOfTrimmedVertices;
    }
  }

  const double beta = (thermalEnergy > 0.0) ? 1.0 / thermalEnergy : 0.0;
  const double largestPlausible = largestPlausibleTriangleArea(framework.unitCell.cell, field.gridSize);
  const double deepest = field.deepestEnergy();
  const double capThreshold =
      static_cast<double>(iso) - std::max({0.02 * std::abs(static_cast<double>(iso)), 0.01 * std::abs(deepest), 1.0e-6});

  result.patches.reserve(corners.size() / 3);
  result.deepestWell = std::numeric_limits<double>::infinity();
  double depthMoment = 0.0;
  result.area = 0.0;
  result.numberOfTriangles = 0;
  result.numberOfRejectedTriangles = 0;

  for (std::size_t i = 0; i + 2 < corners.size(); i += 3)
  {
    double3 p1 = framework.unitCell.cell * corners[i];
    double3 p2 = framework.unitCell.cell * corners[i + 1];
    double3 p3 = framework.unitCell.cell * corners[i + 2];
    double area = 0.5 * double3::cross(p2 - p1, p3 - p1).length();
    if (!(std::isfinite(area) && area < largestPlausible))
    {
      ++result.numberOfRejectedTriangles;
      continue;
    }

    SheetPatch patch;
    patch.area = area;
    patch.energy[0] = energies[i];
    patch.energy[1] = energies[i + 1];
    patch.energy[2] = energies[i + 2];
    double meanEnergy = (patch.energy[0] + patch.energy[1] + patch.energy[2]) / 3.0;
    if (!(meanEnergy < capThreshold)) continue;

    result.patches.push_back(patch);
    result.area += area;
    ++result.numberOfTriangles;

    double weight = 1.0;
    if (beta > 0.0)
    {
      double w0 = std::exp(std::min(-patch.energy[0] * beta, 700.0));
      double w1 = std::exp(std::min(-patch.energy[1] * beta, 700.0));
      double w2 = std::exp(std::min(-patch.energy[2] * beta, 700.0));
      weight = (w0 + w1 + w2) / 3.0;
    }
    result.weightedArea += area * weight;
    depthMoment += area * meanEnergy;
    result.deepestWell = std::min({result.deepestWell, patch.energy[0], patch.energy[1], patch.energy[2]});
  }

  if (result.area > 0.0) result.meanDepth = depthMoment / result.area;
  if (!std::isfinite(result.deepestWell)) result.deepestWell = 0.0;

  if (framework.mass > 0.0)
  {
    constexpr double angstromSquaredToSquareMetrePerMol = 6.0221419947e3;
    result.gravimetricArea = result.area * angstromSquaredToSquareMetrePerMol / framework.mass;
    result.gravimetricWeightedArea = result.weightedArea * angstromSquaredToSquareMetrePerMol / framework.mass;
  }
  if (framework.unitCell.volume > 0.0)
  {
    result.volumetricArea = 1.0e4 * result.area / framework.unitCell.volume;
  }

  std::vector<float> filamentField(numberOfVoxels, 0.0f);
  std::size_t interiorPoints = 0;
  const double voxelVolume = framework.unitCell.volume / static_cast<double>(numberOfVoxels);
  std::vector<std::size_t> filamentIndex;
  for (std::size_t i = 0; i < numberOfVoxels; ++i)
  {
    float value = std::max(field.reliability[i] - static_cast<float>(wellFilamentReliabilityThreshold),
                           std::max(field.distance[i], static_cast<float>(energyScale) * (field.energy[i] - iso)));
    filamentField[i] = value;
    if (value < 0.0f)
    {
      ++interiorPoints;
      FilamentVoxel voxel;
      voxel.volume = voxelVolume;
      voxel.energy = static_cast<double>(field.energy[i]);
      result.filamentVoxels.push_back(voxel);
      filamentIndex.push_back(i);
      result.filamentVolume += voxelVolume;
    }
  }

  std::vector<double3> filamentCorners;
  if (interiorPoints > 0)
  {
    filamentCorners = trianglesOf(backend, filamentField, field.gridSize, 0.0);
    filamentCorners = removeFilamentSpecks(framework.unitCell, std::move(filamentCorners));
    IsosurfaceArea filament = accumulateTriangleAreas(framework.unitCell.cell, field.gridSize, filamentCorners);
    result.filamentArea = filament.area;
    result.numberOfFilamentTriangles = filament.numberOfTriangles;
  }

  // Diagnostic cylinder of the whole overlay. Packing uses the ridge and the split mesh instead: a
  // pancake is not a tube, and mixing the two into one A²/(4πV) contaminates the 1-D file.
  if (result.filamentVolume > 0.0 && result.filamentArea > 0.0)
  {
    result.filamentLength =
        result.filamentArea * result.filamentArea / (4.0 * std::numbers::pi * result.filamentVolume);
  }

  assignFilamentPacking(result, field, framework.unitCell, filamentIndex, filamentCorners);

  std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - time_begin;
  result.seconds = elapsed.count();
  return result;
}
