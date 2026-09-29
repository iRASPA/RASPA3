module;

module spatial_decomposition_domain_decomposition;

import std;

import int3;
import double3;

int3 DomainDecomposition::chooseGrid(std::size_t threads, int3 cells, double3 widths, std::optional<int3> requested)
{
  const std::size_t numberOfThreads = std::max<std::size_t>(1, threads);

  if (requested.has_value())
  {
    const int3 g = requested.value();
    if (g.x < 1 || g.y < 1 || g.z < 1 ||
        static_cast<std::size_t>(g.x) * static_cast<std::size_t>(g.y) * static_cast<std::size_t>(g.z) !=
            numberOfThreads)
    {
      throw std::runtime_error(
          std::format("[Spatial decomposition]: the product of 'DomainGrid' [{}, {}, {}] must equal the number of "
                      "threads ({})\n",
                      g.x, g.y, g.z, numberOfThreads));
    }
    if (g.x > cells.x || g.y > cells.y || g.z > cells.z)
    {
      throw std::runtime_error(
          std::format("[Spatial decomposition]: 'DomainGrid' [{}, {}, {}] exceeds the cell grid [{}, {}, {}] (every "
                      "sub-domain must be at least one cell of the neighbour search wide); use fewer threads along "
                      "that axis or a smaller 'VerletSkin'\n",
                      g.x, g.y, g.z, cells.x, cells.y, cells.z));
    }
    return g;
  }

  std::optional<int3> best{};
  double bestArea = std::numeric_limits<double>::max();
  const std::size_t t = numberOfThreads;
  for (std::size_t nx = 1; nx <= t; ++nx)
  {
    if (t % nx != 0) continue;
    const std::size_t remainder = t / nx;
    for (std::size_t ny = 1; ny <= remainder; ++ny)
    {
      if (remainder % ny != 0) continue;
      const std::size_t nz = remainder / ny;
      if (nx > static_cast<std::size_t>(cells.x) || ny > static_cast<std::size_t>(cells.y) ||
          nz > static_cast<std::size_t>(cells.z))
      {
        continue;
      }
      // Surface of one sub-domain: 2 (ly lz + lx lz + lx ly) with l_i = width_i / n_i
      const double lx = widths.x / static_cast<double>(nx);
      const double ly = widths.y / static_cast<double>(ny);
      const double lz = widths.z / static_cast<double>(nz);
      const double area = ly * lz + lx * lz + lx * ly;
      if (area < bestArea - 1e-12 * std::max(1.0, bestArea))
      {
        bestArea = area;
        best = int3(static_cast<std::int32_t>(nx), static_cast<std::int32_t>(ny), static_cast<std::int32_t>(nz));
      }
    }
  }

  if (!best.has_value())
  {
    throw std::runtime_error(std::format(
        "[Spatial decomposition]: {} threads cannot be laid out on a cell grid of [{}, {}, {}] cells "
        "(at most {} sub-domains fit in this box, every sub-domain must be at least one cell wide along every axis); "
        "use fewer threads or a smaller 'VerletSkin'\n",
        numberOfThreads, cells.x, cells.y, cells.z,
        static_cast<std::size_t>(cells.x) * static_cast<std::size_t>(cells.y) * static_cast<std::size_t>(cells.z)));
  }
  return best.value();
}

namespace
{
// Cuts a set of coordinates (given by index into `values`) into `parts` groups of equal size; the cut positions
// are the midpoints between the neighbouring coordinates of consecutive groups. Returns parts + 1 positions
// starting with 0 and ending with 1. Empty or tiny sets fall back to equidistant cuts.
std::vector<double> balancedCuts(std::vector<std::uint32_t>& indices, const std::vector<double>& values, int parts)
{
  std::vector<double> cuts(static_cast<std::size_t>(parts) + 1);
  cuts.front() = 0.0;
  cuts.back() = 1.0;
  if (parts == 1) return cuts;
  const std::size_t count = indices.size();
  if (count < static_cast<std::size_t>(parts))
  {
    for (int p = 1; p < parts; ++p) cuts[static_cast<std::size_t>(p)] = static_cast<double>(p) / parts;
    return cuts;
  }
  std::sort(indices.begin(), indices.end(), [&](std::uint32_t a, std::uint32_t b) { return values[a] < values[b]; });
  for (int p = 1; p < parts; ++p)
  {
    const std::size_t boundary = (static_cast<std::size_t>(p) * count) / static_cast<std::size_t>(parts);
    double cut = 0.5 * (values[indices[boundary - 1]] + values[indices[boundary]]);
    // keep the cuts strictly increasing
    cut = std::max(cut, std::nextafter(cuts[static_cast<std::size_t>(p) - 1], 2.0));
    cuts[static_cast<std::size_t>(p)] = std::min(cut, 1.0);
  }
  return cuts;
}

int slabOf(const std::vector<double>& cuts, double value)
{
  // the slab p with cuts[p] <= value < cuts[p + 1]
  const auto it = std::upper_bound(cuts.begin(), cuts.end(), value);
  const int p = static_cast<int>(std::distance(cuts.begin(), it)) - 1;
  return std::clamp(p, 0, static_cast<int>(cuts.size()) - 2);
}
}  // namespace

void DomainDecomposition::build(int3 g, std::span<const double3> fractionalPositions)
{
  grid = g;
  const std::size_t count = fractionalPositions.size();
  std::vector<double> sx(count), sy(count), sz(count);
  for (std::size_t i = 0; i < count; ++i)
  {
    sx[i] = fractionalPositions[i].x;
    sy[i] = fractionalPositions[i].y;
    sz[i] = fractionalPositions[i].z;
  }

  std::vector<std::uint32_t> all(count);
  std::iota(all.begin(), all.end(), 0u);
  planesX = balancedCuts(all, sx, grid.x);

  planesY.assign(static_cast<std::size_t>(grid.x), {});
  planesZ.assign(static_cast<std::size_t>(grid.x), std::vector<std::vector<double>>(static_cast<std::size_t>(grid.y)));
  std::vector<std::vector<std::uint32_t>> slabs(static_cast<std::size_t>(grid.x));
  for (const std::uint32_t i : all) slabs[static_cast<std::size_t>(slabOf(planesX, sx[i]))].push_back(i);
  for (int ix = 0; ix < grid.x; ++ix)
  {
    std::vector<std::uint32_t>& slab = slabs[static_cast<std::size_t>(ix)];
    planesY[static_cast<std::size_t>(ix)] = balancedCuts(slab, sy, grid.y);
    std::vector<std::vector<std::uint32_t>> columns(static_cast<std::size_t>(grid.y));
    for (const std::uint32_t i : slab)
    {
      columns[static_cast<std::size_t>(slabOf(planesY[static_cast<std::size_t>(ix)], sy[i]))].push_back(i);
    }
    for (int iy = 0; iy < grid.y; ++iy)
    {
      planesZ[static_cast<std::size_t>(ix)][static_cast<std::size_t>(iy)] =
          balancedCuts(columns[static_cast<std::size_t>(iy)], sz, grid.z);
    }
  }
}

std::uint32_t DomainDecomposition::ownerOf(const double3& fractional) const
{
  const int ix = slabOf(planesX, fractional.x);
  const int iy = slabOf(planesY[static_cast<std::size_t>(ix)], fractional.y);
  const int iz = slabOf(planesZ[static_cast<std::size_t>(ix)][static_cast<std::size_t>(iy)], fractional.z);
  return static_cast<std::uint32_t>((iz * grid.y + iy) * grid.x + ix);
}

std::string DomainDecomposition::status() const
{
  return std::format("sub-domain grid {} x {} x {} (cuts balanced on the atom count)", grid.x, grid.y, grid.z);
}
