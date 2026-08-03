module;

module sampled_roadmap;

import std;

import int3;
import double3;
import double3x3;
import unit_cell;
import randomnumbers;
import voronoi_network;
import channel_analysis;
import sampled_structure;
import sampling_backend;

namespace
{
// The kept samples sorted into a periodic grid of bins, so that the ones near a sample can be found without
// looking at all of them. The bins are at least as wide as the link distance in every direction, which is
// what makes the twenty-seven bins around a sample's own the only ones that can hold a link of it.
struct SampleGrid
{
  int3 counts{1, 1, 1};
  std::vector<std::size_t> start;     // where each bin's samples begin in `contents`
  std::vector<std::size_t> contents;  // sample indices, bin by bin

  SampleGrid(const UnitCell &unitCell, const std::vector<double3> &fractional, double spacing)
  {
    double3 widths = unitCell.perpendicularWidths();

    // Never fewer than one bin along an axis, and never so many that the grid outgrows what it indexes.
    auto along = [&](double width)
    {
      double wanted = spacing > 0.0 ? std::floor(width / spacing) : 1.0;
      return static_cast<std::int32_t>(std::clamp(wanted, 1.0, 512.0));
    };
    counts = int3(along(widths.x), along(widths.y), along(widths.z));

    std::size_t numberOfBins = static_cast<std::size_t>(counts.x) * static_cast<std::size_t>(counts.y) *
                               static_cast<std::size_t>(counts.z);

    std::vector<std::size_t> tally(numberOfBins, 0);
    std::vector<std::size_t> bin(fractional.size(), 0);
    for (std::size_t index = 0; index < fractional.size(); ++index)
    {
      bin[index] = binOf(fractional[index]);
      ++tally[bin[index]];
    }

    start.assign(numberOfBins + 1, 0);
    std::partial_sum(tally.begin(), tally.end(), start.begin() + 1);

    std::vector<std::size_t> cursor(start.begin(), start.end() - 1);
    contents.resize(fractional.size());
    for (std::size_t index = 0; index < fractional.size(); ++index)
    {
      contents[cursor[bin[index]]++] = index;
    }
  }

  std::size_t binOf(const double3 &fractional) const
  {
    auto axis = [](double coordinate, std::int32_t count)
    {
      double wrapped = coordinate - std::floor(coordinate);
      std::int32_t index = static_cast<std::int32_t>(wrapped * static_cast<double>(count));
      return std::clamp(index, std::int32_t{0}, count - 1);
    };

    std::int32_t i = axis(fractional.x, counts.x);
    std::int32_t j = axis(fractional.y, counts.y);
    std::int32_t k = axis(fractional.z, counts.z);

    return (static_cast<std::size_t>(k) * static_cast<std::size_t>(counts.y) + static_cast<std::size_t>(j)) *
               static_cast<std::size_t>(counts.x) +
           static_cast<std::size_t>(i);
  }

  // The bin indices along one axis that neighbour `index`, without repeating any: the three around it where
  // there are three to have, and otherwise every bin there is, the cell being too thin to divide.
  static void neighboursAlong(std::int32_t index, std::int32_t count, std::vector<std::int32_t> &into)
  {
    into.clear();
    if (count >= 3)
    {
      into.push_back((index + count - 1) % count);
      into.push_back(index);
      into.push_back((index + 1) % count);
    }
    else
    {
      for (std::int32_t i = 0; i < count; ++i) into.push_back(i);
    }
  }
};

// How many distinct pockets are carried into the stage that tries every pair of them as a hop. That stage
// is quadratic in this number and each pair costs a pass over the atoms, so it is bounded; a crystal with
// more pockets than this in a unit cell keeps the roomiest of them, which are the ones the diameters are
// decided by.
constexpr std::size_t maximumNumberOfPockets = 1024;

// How many of a pocket's hops are tried again with the straight line bent out of the way. Every hop of
// every pair is a great many, and the ones worth the trouble are the ones the straight line already nearly
// manages, so they are ranked by that and the best of them taken.
constexpr std::size_t refinedHopsPerPocket = 16;

// A way that has been asked for but not yet measured: which two nodes it runs between, which image of the
// second it reaches, and the vector from the first to it.
struct Candidate
{
  std::size_t from;
  std::size_t to;
  int3 delta;
  double3 displacement;
  std::uint32_t seed;
};

// Ways are measured in batches, so that whatever is measuring them has enough to do to be worth asking.
// Cases are handed in one at a time and come back out one at a time, through `consume`; what happens in
// between is that they are held until there are enough of them, sent off together, and handed back in the
// order they were given. Nothing here knows whether that happened on this processor or on a device.
class WayBatch
{
 public:
  WayBatch(const SamplingBackend &backend, const SampledStructure &structure, bool bend)
      : backend(backend), structure(structure), bend(bend)
  {
  }

  template <typename Consume>
  void add(const Candidate &candidate, const double3 &origin, Consume &&consume)
  {
    this->pending.push_back(candidate);
    this->origins.push_back(origin);

    if (this->pending.size() >= samplingBatchSize) this->flush(consume);
  }

  template <typename Consume>
  void flush(Consume &&consume)
  {
    if (this->pending.empty()) return;

    this->displacements.clear();
    this->displacements.reserve(this->pending.size());
    for (const Candidate &candidate : this->pending) this->displacements.push_back(candidate.displacement);

    this->ways.assign(this->pending.size(), SampledWay{});

    if (this->bend)
    {
      this->seeds.clear();
      this->seeds.reserve(this->pending.size());
      for (const Candidate &candidate : this->pending) this->seeds.push_back(candidate.seed);

      this->backend.widestWays(this->structure, this->origins, this->displacements, this->seeds,
                               samplingRefinementDepth, this->ways);
    }
    else
    {
      this->backend.straightWays(this->structure, this->origins, this->displacements, this->ways);
    }

    for (std::size_t i = 0; i < this->pending.size(); ++i) consume(this->pending[i], this->ways[i]);

    this->pending.clear();
    this->origins.clear();
  }

 private:
  const SamplingBackend &backend;
  const SampledStructure &structure;
  bool bend;

  std::vector<Candidate> pending;
  std::vector<double3> origins;
  std::vector<double3> displacements;
  std::vector<std::uint32_t> seeds;
  std::vector<SampledWay> ways;
};
}  // namespace

SampledRoadmap SampledRoadmap::build(const SampledStructure &structure, const SamplingBackend &backend,
                                     std::optional<std::size_t> numberOfIterations,
                                     std::optional<std::size_t> numberOfInnerSteps)
{
  SampledRoadmap self;
  VoronoiNetwork &network = self.network;

  RandomNumber random{std::nullopt};

  std::size_t number_of_iterations = numberOfIterations.value_or(100000);
  std::size_t number_of_inner_steps = numberOfInnerSteps.value_or(200);

  self.backend = backend.name;

  std::chrono::steady_clock::time_point time_begin = std::chrono::steady_clock::now();

  const UnitCell &unitCell = structure.unitCell;
  const double3x3 &cell = unitCell.cell;

  // Throw the points. They are drawn in one stream, so that which point is which does not depend on how the
  // work was shared out, and the clearance at each is worked out afterwards, that being the part worth
  // handing over in bulk.
  std::vector<double3> drawn(number_of_iterations);
  std::vector<double3> drawnCartesian(number_of_iterations);
  for (std::size_t i = 0; i < number_of_iterations; ++i)
  {
    drawn[i] = double3(random.uniform(), random.uniform(), random.uniform());
    drawnCartesian[i] = cell * drawn[i];
  }

  std::vector<double> drawnClearance(number_of_iterations, 0.0);
  for (std::size_t begin = 0; begin < number_of_iterations; begin += samplingBatchSize)
  {
    std::size_t count = std::min(samplingBatchSize, number_of_iterations - begin);
    backend.clearances(structure, std::span(drawnCartesian).subspan(begin, count),
                       std::span(drawnClearance).subspan(begin, count));
  }

  // Keep the ones that landed in the void. What share of them that is, is the void fraction, and it comes
  // out of this route whether it was asked for or not.
  std::vector<double3> fractional;
  std::vector<double3> position;
  std::vector<double> radius;
  for (std::size_t i = 0; i < number_of_iterations; ++i)
  {
    if (drawnClearance[i] <= 0.0) continue;

    fractional.push_back(drawn[i]);
    position.push_back(drawnCartesian[i]);
    radius.push_back(drawnClearance[i]);
  }

  drawn.clear();
  drawn.shrink_to_fit();
  drawnCartesian.clear();
  drawnCartesian.shrink_to_fit();
  drawnClearance.clear();
  drawnClearance.shrink_to_fit();

  std::size_t numberOfNodes = position.size();

  self.numberOfSamples = number_of_iterations;
  self.numberOfVoidSamples = numberOfNodes;
  self.voidFraction =
      number_of_iterations > 0 ? static_cast<double>(numberOfNodes) / static_cast<double>(number_of_iterations) : 0.0;

  // How far apart two samples may be and still be joined. The points fall at a density of one per cell
  // volume over sample count wherever they can fall at all, so a multiple of that spacing sets the number of
  // hops a sample has: about seventeen at this multiple. More hops are never worse -- a way that is offered
  // more ways round can only find a wider one -- but they cost, and beyond this the sample roadmap is not
  // where the bottlenecks are coming from anyway.
  const double linkFactor = 1.6;
  self.linkDistance =
      number_of_iterations > 0
          ? linkFactor * std::cbrt(unitCell.volume / static_cast<double>(number_of_iterations))
          : 0.0;

  SampleGrid grid(unitCell, fractional, self.linkDistance);

  // The roadmap, in the shape every pore network here is kept in, so that whatever reads a diagram reads it
  // by the same arithmetic.
  network.unitCell = unitCell;
  network.atomRadii = structure.radii;
  network.atomPositionsFractional.reserve(structure.size());
  for (const double3 &atom : structure.positions)
  {
    double3 s = unitCell.inverseCell * atom;
    network.atomPositionsFractional.push_back(
        double3(s.x - std::floor(s.x), s.y - std::floor(s.y), s.z - std::floor(s.z)));
  }

  network.nodes.resize(numberOfNodes);
  for (std::size_t i = 0; i < numberOfNodes; ++i)
  {
    VoronoiNode &node = network.nodes[i];
    node.position = position[i];
    node.fractional = fractional[i];
    node.radius = radius[i];
    node.maximalRadius = radius[i];
    node.maximalPosition = position[i];
  }

  // Every hop goes in both ways round, as the network's own builders store theirs, and carries where its
  // narrowest point is and which way it runs there. The diagram routes can only supply that for an
  // Apollonius edge; here every hop has it, because every hop was measured along rather than solved for.
  auto addLink = [&](std::size_t from, std::size_t to, const int3 &delta, const SampledWay &way, double length)
  {
    VoronoiEdge forward;
    forward.from = from;
    forward.to = to;
    forward.delta = delta;
    forward.radius = way.radius;
    forward.length = length;
    forward.bottleneckPosition = way.position;
    forward.bottleneckDirection = way.direction;
    forward.hasBottleneckGeometry = true;
    network.edges.push_back(forward);

    VoronoiEdge backward = forward;
    backward.from = to;
    backward.to = from;
    backward.delta = int3(-delta.x, -delta.y, -delta.z);
    network.edges.push_back(backward);
  };

  auto keepOpenHop = [&](const Candidate &candidate, const SampledWay &way)
  {
    if (way.radius <= 0.0) return;
    addLink(candidate.from, candidate.to, candidate.delta, way, candidate.displacement.length());
  };

  // The image of one sample that lies nearest another, and which image that is.
  auto nearestImage = [&](const double3 &from, const double3 &to)
  {
    double3 dr = to - from;
    double3 s = unitCell.inverseCell * dr;
    int3 image(static_cast<std::int32_t>(std::rint(s.x)), static_cast<std::int32_t>(std::rint(s.y)),
               static_cast<std::int32_t>(std::rint(s.z)));
    double3 displacement =
        dr - cell * double3(static_cast<double>(image.x), static_cast<double>(image.y),
                            static_cast<double>(image.z));
    return std::pair<int3, double3>{int3(-image.x, -image.y, -image.z), displacement};
  };

  // The bins of the grid that can hold a neighbour of a point, given where in the grid it sits.
  auto neighbourhoodOf = [&](const double3 &fractionalPoint, std::vector<std::size_t> &bins)
  {
    auto axisIndex = [](double coordinate, std::int32_t count)
    {
      double fraction = coordinate - std::floor(coordinate);
      return std::clamp(static_cast<std::int32_t>(fraction * static_cast<double>(count)), std::int32_t{0},
                        count - 1);
    };

    std::vector<std::int32_t> alongX, alongY, alongZ;
    SampleGrid::neighboursAlong(axisIndex(fractionalPoint.x, grid.counts.x), grid.counts.x, alongX);
    SampleGrid::neighboursAlong(axisIndex(fractionalPoint.y, grid.counts.y), grid.counts.y, alongY);
    SampleGrid::neighboursAlong(axisIndex(fractionalPoint.z, grid.counts.z), grid.counts.z, alongZ);

    bins.clear();
    for (std::int32_t k : alongZ)
    {
      for (std::int32_t j : alongY)
      {
        for (std::int32_t i : alongX)
        {
          bins.push_back((static_cast<std::size_t>(k) * static_cast<std::size_t>(grid.counts.y) +
                          static_cast<std::size_t>(j)) *
                             static_cast<std::size_t>(grid.counts.x) +
                         static_cast<std::size_t>(i));
        }
      }
    }
  };

  // Join the samples that are near enough to each other, wherever a sphere can be slid from one to the
  // other. Both ends being in the void says nothing about what lies between them, which is the whole of
  // what is being asked here.
  {
    WayBatch hops(backend, structure, false);
    std::vector<std::size_t> bins;

    for (std::size_t i = 0; i < numberOfNodes; ++i)
    {
      neighbourhoodOf(fractional[i], bins);

      for (std::size_t bin : bins)
      {
        for (std::size_t slot = grid.start[bin]; slot < grid.start[bin + 1]; ++slot)
        {
          std::size_t other = grid.contents[slot];

          // Each pair is looked at once, by the lower of the two.
          if (other <= i) continue;

          auto [delta, displacement] = nearestImage(position[i], position[other]);
          if (double3::dot(displacement, displacement) >= self.linkDistance * self.linkDistance) continue;

          hops.add(Candidate{.from = i, .to = other, .delta = delta, .displacement = displacement, .seed = 0},
                   position[i], keepOpenHop);
        }
      }
    }

    hops.flush(keepOpenHop);
  }

  // A sample no hop of its leads anywhere roomier is the top of its pocket as far as the sample can see.
  // Those are the ones worth walking further uphill; the rest are already beaten by one that is.
  std::vector<char> isPeak(numberOfNodes, 1);
  for (const VoronoiEdge &edge : network.edges)
  {
    bool beaten = radius[edge.to] > radius[edge.from] ||
                  (radius[edge.to] == radius[edge.from] && edge.to < edge.from);
    if (beaten) isPeak[edge.from] = 0;
  }

  std::vector<std::size_t> peaks;
  for (std::size_t i = 0; i < numberOfNodes; ++i)
  {
    if (isPeak[i] != 0) peaks.push_back(i);
  }
  self.numberOfPeaks = peaks.size();

  // Walk each peak uphill, so that what is reported is the true deepest point of the pocket rather than the
  // best point that happened to be drawn.
  {
    std::vector<SampledPeak> starts(peaks.size());
    std::vector<std::uint32_t> seeds(peaks.size());
    for (std::size_t i = 0; i < peaks.size(); ++i)
    {
      starts[i] = SampledPeak{.radius = radius[peaks[i]], .position = position[peaks[i]]};
      seeds[i] = static_cast<std::uint32_t>(peaks[i]);
    }

    std::vector<SampledPeak> ends(peaks.size());
    for (std::size_t begin = 0; begin < peaks.size(); begin += samplingBatchSize)
    {
      std::size_t count = std::min(samplingBatchSize, peaks.size() - begin);
      backend.walksUphill(structure, std::span(starts).subspan(begin, count),
                          std::span(seeds).subspan(begin, count), number_of_inner_steps,
                          std::span(ends).subspan(begin, count));
    }

    for (std::size_t i = 0; i < peaks.size(); ++i)
    {
      network.nodes[peaks[i]].maximalRadius = ends[i].radius;
      network.nodes[peaks[i]].maximalPosition = ends[i].position;
    }
  }

  // The peaks that walked to the same place are the same pocket found twice, and what is left after
  // dropping every one that ended up inside a roomier one is a list of the distinct pockets of the crystal.
  // Keeping the roomiest of them is what bounds the work of the stage below, which is quadratic in how many
  // there are.
  std::vector<std::size_t> pocketOf;
  {
    std::vector<std::size_t> ordered = peaks;
    std::sort(ordered.begin(), ordered.end(), [&](std::size_t a, std::size_t b)
              { return network.nodes[a].maximalRadius > network.nodes[b].maximalRadius; });

    for (std::size_t node : ordered)
    {
      if (pocketOf.size() >= maximumNumberOfPockets) break;

      const VoronoiNode &candidate = network.nodes[node];

      bool insideARoomierOne = false;
      for (std::size_t kept : pocketOf)
      {
        const VoronoiNode &pocket = network.nodes[kept];
        double3 dr =
            unitCell.applyPeriodicBoundaryConditions(candidate.maximalPosition - pocket.maximalPosition);
        if (double3::dot(dr, dr) < pocket.maximalRadius * pocket.maximalRadius)
        {
          insideARoomierOne = true;
          break;
        }
      }

      if (!insideARoomierOne) pocketOf.push_back(node);
    }
  }
  self.numberOfPocketCentres = pocketOf.size();

  // The deepest point of each pocket joins the roadmap as a node of its own, and every pair of them is
  // tried as a way from one to the other.
  //
  // This is what the sampled roadmap cannot reach by having more points in it. A hop between two samples is
  // as good as the line between them happens to be, and a line drawn between two points that merely fell
  // near the middles of neighbouring pockets misses the middle of the window between them by about the
  // spacing of the sample, however fine that spacing is made. The line between the two pocket centres does
  // not: it is the line the passage was built around, and its narrowest point is the window itself. So the
  // bottlenecks are found by the same handful of long ways the diagram routes read them from, while the
  // sample underneath goes on covering everything those ways leave out.
  std::size_t firstPocketNode = network.nodes.size();
  for (std::size_t node : pocketOf)
  {
    double3 s = unitCell.inverseCell * network.nodes[node].maximalPosition;
    double3 wrapped(s.x - std::floor(s.x), s.y - std::floor(s.y), s.z - std::floor(s.z));

    VoronoiNode centre;
    centre.fractional = wrapped;
    centre.position = cell * wrapped;
    centre.radius = network.nodes[node].maximalRadius;
    centre.maximalRadius = centre.radius;
    centre.maximalPosition = centre.position;
    network.nodes.push_back(centre);
  }

  // Long enough for a pocket to reach its own next repeat, which is the whole of what a straight channel
  // through the crystal is, and no longer: past one repeat there is always a nearer image of the same end.
  const double pocketLinkDistance =
      std::max({unitCell.lengthA, unitCell.lengthB, unitCell.lengthC});

  // Every pocket against every other, straight first. The ones the straight line comes nearest to managing
  // are the ones that stand for a real passage, open or merely bent, so those are the ones tried again with
  // the line allowed to bend.
  {
    // Per pocket, its straight hops ranked by how nearly they succeeded, so that only the best of them are
    // paid for again.
    std::vector<std::vector<std::pair<double, Candidate>>> straightOf(pocketOf.size());

    WayBatch straight(backend, structure, false);
    auto rank = [&](const Candidate &candidate, const SampledWay &way)
    { straightOf[candidate.from - firstPocketNode].push_back({way.radius, candidate}); };

    for (std::size_t a = 0; a < pocketOf.size(); ++a)
    {
      std::size_t from = firstPocketNode + a;
      const double3 &origin = network.nodes[from].position;

      for (std::size_t b = a; b < pocketOf.size(); ++b)
      {
        std::size_t to = firstPocketNode + b;

        for (std::int32_t nz = -1; nz <= 1; ++nz)
        {
          for (std::int32_t ny = -1; ny <= 1; ++ny)
          {
            for (std::int32_t nx = -1; nx <= 1; ++nx)
            {
              // A pocket paired with itself in this cell is no hop, and a pocket paired with itself the
              // other way round is the hop already made.
              if (a == b && (nx < 0 || (nx == 0 && (ny < 0 || (ny == 0 && nz <= 0))))) continue;

              double3 displacement = network.nodes[to].position +
                                     cell * double3(static_cast<double>(nx), static_cast<double>(ny),
                                                    static_cast<double>(nz)) -
                                     origin;

              if (double3::dot(displacement, displacement) > pocketLinkDistance * pocketLinkDistance) continue;

              straight.add(Candidate{.from = from,
                                     .to = to,
                                     .delta = int3(nx, ny, nz),
                                     .displacement = displacement,
                                     .seed = 0},
                           origin, rank);
            }
          }
        }
      }
    }
    straight.flush(rank);

    // The best few of each pocket's hops, bent; the rest kept as the straight line left them.
    WayBatch bent(backend, structure, true);
    std::uint32_t seed = 0;

    for (std::size_t a = 0; a < pocketOf.size(); ++a)
    {
      std::vector<std::pair<double, Candidate>> &hops = straightOf[a];

      std::size_t refined = std::min(refinedHopsPerPocket, hops.size());
      std::partial_sort(hops.begin(), hops.begin() + static_cast<std::ptrdiff_t>(refined), hops.end(),
                        [](const auto &p, const auto &q) { return p.first > q.first; });

      for (std::size_t rankIndex = 0; rankIndex < hops.size(); ++rankIndex)
      {
        if (rankIndex < refined)
        {
          Candidate candidate = hops[rankIndex].second;
          candidate.seed = seed++;
          bent.add(candidate, network.nodes[candidate.from].position, keepOpenHop);
        }
        else if (hops[rankIndex].first > 0.0)
        {
          // Nothing more will be learnt about these, so they go in as they are. The direction is the line
          // itself and the narrowest point is where the sphere was smallest along it, which is what the
          // straight pass already found; it is recovered from the geometry rather than kept, the ranking
          // above having no use for it.
          const Candidate &candidate = hops[rankIndex].second;
          double length = candidate.displacement.length();
          SampledWay way{.radius = hops[rankIndex].first,
                         .position = network.nodes[candidate.from].position + 0.5 * candidate.displacement,
                         .direction = length > 0.0 ? (1.0 / length) * candidate.displacement
                                                   : double3(0.0, 0.0, 1.0)};
          addLink(candidate.from, candidate.to, candidate.delta, way, length);
        }
      }

      hops.clear();
      hops.shrink_to_fit();
    }
    bent.flush(keepOpenHop);
  }

  // And each pocket into the sample around it, so that the pockets are part of the same roadmap rather than
  // a graph beside it.
  {
    WayBatch hops(backend, structure, false);
    std::vector<std::size_t> bins;

    for (std::size_t a = 0; a < pocketOf.size(); ++a)
    {
      std::size_t from = firstPocketNode + a;
      const double3 &origin = network.nodes[from].position;

      neighbourhoodOf(network.nodes[from].fractional, bins);

      for (std::size_t bin : bins)
      {
        for (std::size_t slot = grid.start[bin]; slot < grid.start[bin + 1]; ++slot)
        {
          std::size_t other = grid.contents[slot];

          auto [delta, displacement] = nearestImage(origin, position[other]);
          if (double3::dot(displacement, displacement) >= self.linkDistance * self.linkDistance) continue;

          hops.add(
              Candidate{.from = from, .to = other, .delta = delta, .displacement = displacement, .seed = 0},
              origin, keepOpenHop);
        }
      }
    }

    hops.flush(keepOpenHop);
  }

  // What is left over after all that may be a pocket or it may be a piece of a channel the sample failed to
  // join up, and the difference matters more than anywhere else it comes up: a piece mistaken for a pocket
  // is counted as room a molecule cannot reach and, if the structure is being blocked, has a sphere written
  // over part of a working pore.
  //
  // The two are told apart by trying. From the roomiest node of a component that does not run anywhere, a
  // way is sought to the roomiest node of every other component within reach, bent where the straight line
  // will not do, exactly as between the pockets above. A way that opens joins them and the piece was never
  // a pocket; a way that stays shut leaves the pocket a pocket, and now for a reason rather than for want
  // of a hop. Joining changes what the components are, so it is done twice over.
  //
  // No way is ever invented, every one of them being measured, so this can only ever join and never split.
  for (std::size_t round = 0; round < 2; ++round)
  {
    ChannelAnalysis components = ChannelAnalysis::compute(network, 0.0);
    if (components.numberOfPockets == 0 || components.pores.size() < 2) break;

    // One node stands for each component: the roomiest of them, which is the likeliest end of a way and the
    // furthest from the walls at either end of it.
    std::vector<std::size_t> representative;
    representative.reserve(components.pores.size());
    for (const VoronoiPore &pore : components.pores)
    {
      std::size_t best = pore.nodeIndices.front();
      for (std::size_t node : pore.nodeIndices)
      {
        if (network.nodes[node].radius > network.nodes[best].radius) best = node;
      }
      representative.push_back(best);
    }

    std::size_t joined = 0;
    auto keepJoin = [&](const Candidate &candidate, const SampledWay &way)
    {
      if (way.radius <= 0.0) return;
      addLink(candidate.from, candidate.to, candidate.delta, way, candidate.displacement.length());
      ++joined;
    };

    WayBatch ways(backend, structure, true);
    std::uint32_t seed = static_cast<std::uint32_t>(numberOfNodes + round * components.pores.size());

    for (std::size_t pore = 0; pore < components.pores.size(); ++pore)
    {
      if (components.pores[pore].isChannel) continue;

      std::size_t from = representative[pore];
      const double3 &origin = network.nodes[from].position;

      for (std::size_t other = 0; other < components.pores.size(); ++other)
      {
        if (other == pore) continue;

        std::size_t to = representative[other];

        auto [delta, displacement] = nearestImage(origin, network.nodes[to].position);
        if (double3::dot(displacement, displacement) > pocketLinkDistance * pocketLinkDistance) continue;

        ways.add(
            Candidate{.from = from, .to = to, .delta = delta, .displacement = displacement, .seed = seed++},
            origin, keepJoin);
      }
    }
    ways.flush(keepJoin);

    self.numberOfJoins += joined;
    if (joined == 0) break;
  }

  self.numberOfLinks = network.edges.size() / 2;
  self.firstPocketNode = firstPocketNode;
  self.numberOfInnerSteps = number_of_inner_steps;

  // How the roadmap finally falls apart, and which of the pieces the sample was fine enough to speak for.
  self.components = ChannelAnalysis::compute(network, 0.0);
  self.poreIsResolved.assign(self.components.pores.size(), 0);

  for (std::size_t pore = 0; pore < self.components.pores.size(); ++pore)
  {
    const VoronoiPore &piece = self.components.pores[pore];
    if (piece.isChannel)
    {
      self.poreIsResolved[pore] = 1;
      ++self.numberOfChannels;
      continue;
    }

    // A pocket has a walk uphill ending in it. Those ends are the nodes from `firstPocketNode` on, so a
    // piece holding one of them holds a genuine local maximum of the room available and is a cavity; a
    // piece holding none is a stray the sample could not resolve.
    bool resolved =
        std::ranges::any_of(piece.nodeIndices, [&](std::size_t node) { return node >= firstPocketNode; });

    self.poreIsResolved[pore] = resolved ? 1 : 0;
    if (resolved)
      ++self.numberOfPockets;
    else
      ++self.numberOfUnresolvedPieces;
  }

  std::chrono::duration<double> timing = std::chrono::steady_clock::now() - time_begin;
  self.seconds = timing.count();

  return self;
}

bool SampledRoadmap::isReachable(std::size_t node) const
{
  std::int32_t pore = this->components.nodePoreId[node];
  return pore >= 0 && this->components.pores[static_cast<std::size_t>(pore)].isChannel;
}

bool SampledRoadmap::isShutIn(std::size_t node) const
{
  std::int32_t pore = this->components.nodePoreId[node];
  if (pore < 0) return false;

  std::size_t index = static_cast<std::size_t>(pore);
  return !this->components.pores[index].isChannel && this->poreIsResolved[index] != 0;
}

void SampledRoadmap::writeHeader(std::ostream &stream) const
{
  std::print(stream, "# Number of iterations (points thrown): {}\n", this->numberOfSamples);
  std::print(stream, "# Number of inner-steps (uphill steps per peak): {}\n", this->numberOfInnerSteps);
  std::print(stream, "# Number of samples in the void: {}\n", this->numberOfVoidSamples);
  std::print(stream, "# Void fraction of the sample: {}\n", this->voidFraction);
  std::print(stream, "# Link distance: {} [Å]\n", this->linkDistance);
  std::print(stream, "# Number of links: {}\n", this->numberOfLinks);
  std::print(stream, "# Number of peaks walked uphill: {}\n", this->numberOfPeaks);
  std::print(stream, "# Number of distinct pocket centres found: {}\n", this->numberOfPocketCentres);
  std::print(stream, "# Number of ways found between components left apart: {}\n", this->numberOfJoins);
  std::print(stream, "# Pieces of the void: {} channel(s), {} pocket(s), {} too thinly sampled to say\n",
             this->numberOfChannels, this->numberOfPockets, this->numberOfUnresolvedPieces);
  std::print(stream, "# Roadmap built on: {}\n", this->backend);
  std::print(stream, "# Timing (roadmap): {} [s]\n", this->seconds);
}
