module;

module property_polymer_chains;

import std;

import archive;
import double3;
import atom;
import forcefield;
import component;
import simulationbox;
import json;
import polymer_links;
import property_molecule_shape;

// ---- Accumulator ----------------------------------------------------------------------------------------------------

void PropertyPolymerChains::Accumulator::add(const Moments &values)
{
  ++count;
  for (std::size_t m = 0; m < NumberOfMoments; ++m)
  {
    sum[m] += values[m];
    sumOfSquares[m] += values[m] * values[m];
  }
}

PropertyPolymerChains::Accumulator &PropertyPolymerChains::Accumulator::operator+=(const Accumulator &other)
{
  count += other.count;
  for (std::size_t m = 0; m < NumberOfMoments; ++m)
  {
    sum[m] += other.sum[m];
    sumOfSquares[m] += other.sumOfSquares[m];
  }
  return *this;
}

double PropertyPolymerChains::Accumulator::standardDeviation(Moment moment) const
{
  if (count < 2) return 0.0;
  const double n = static_cast<double>(count);
  const double average = sum[moment] / n;
  return std::sqrt(std::max(0.0, (sumOfSquares[moment] - n * average * average) / (n - 1.0)));
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PropertyPolymerChains::Key &k)
{
  archive << k.length;
  archive << k.topology;
  archive << k.living;
  archive << k.conversionBin;
  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyPolymerChains::Key &k)
{
  archive >> k.length;
  archive >> k.topology;
  archive >> k.living;
  archive >> k.conversionBin;
  return archive;
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PropertyPolymerChains::Accumulator &a)
{
  archive << a.count;
  archive << a.sum;
  archive << a.sumOfSquares;
  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyPolymerChains::Accumulator &a)
{
  archive >> a.count;
  archive >> a.sum;
  archive >> a.sumOfSquares;
  return archive;
}

// ---- property -------------------------------------------------------------------------------------------------------

PropertyPolymerChains::PropertyPolymerChains(std::size_t sampleEvery, std::size_t writeEvery, Sampling sampling,
                                             std::size_t lagSteps, std::size_t numberOfConversionBins,
                                             std::size_t minimumLength, bool massWeighted)
    : sampleEvery(sampleEvery),
      writeEvery(writeEvery),
      sampling(sampling),
      lagSteps(lagSteps),
      numberOfConversionBins(numberOfConversionBins),
      minimumLength(minimumLength),
      massWeighted(massWeighted)
{
  if (sampleEvery == 0 || writeEvery == 0 || numberOfConversionBins == 0)
  {
    throw std::invalid_argument(
        "PropertyPolymerChains: the sampling and writing intervals and the number of conversion bins must be "
        "positive");
  }
}

std::string_view PropertyPolymerChains::topologyName(Topology topology)
{
  switch (topology)
  {
    case Topology::Linear:
      return "linear";
    case Topology::Branched:
      return "branched";
    case Topology::Cyclic:
      return "cyclic";
  }
  return "unknown";
}

std::string_view PropertyPolymerChains::samplingName(Sampling sampling)
{
  return sampling == Sampling::Periodic ? "Periodic" : "BeforeReaction";
}

std::vector<PropertyPolymerChains::Chain> PropertyPolymerChains::findChains(const PolymerLinkTable &table,
                                                                           const std::vector<Component> &components,
                                                                           const SimulationBox &simulationBox,
                                                                           std::span<const Atom> atoms,
                                                                           std::span<const std::uint32_t> seeds)
{
  std::vector<Chain> chains{};
  if (!table.activated || table.bonds.empty()) return chains;

  // the walks start at the seeds, or at the atoms of every polymer bond (a molecule without one is not a chain)
  std::vector<std::uint32_t> starts(seeds.begin(), seeds.end());
  if (seeds.empty())
  {
    for (const PolymerBond &bond : table.bonds)
    {
      starts.push_back(static_cast<std::uint32_t>(table.layout.index(bond.a)));
      starts.push_back(static_cast<std::uint32_t>(table.layout.index(bond.b)));
    }
  }

  // Per atom the chain it was assigned to and its index in that chain. Each walk unwraps its chain by adding the
  // minimum-image bond vectors; a bond between two visited atoms whose unwrapped positions are not one bond apart
  // closes a cycle through a periodic image, a chain that spans the box. Per chain the local adjacency is kept for
  // the end-to-end sweep.
  constexpr std::size_t unassigned = std::numeric_limits<std::size_t>::max();
  std::vector<std::size_t> chainOf(table.layout.numberOfAtoms, unassigned);
  std::vector<std::size_t> localIndex(table.layout.numberOfAtoms, unassigned);
  std::vector<std::vector<std::vector<std::size_t>>> adjacency{};
  for (std::uint32_t start : starts)
  {
    if (start >= chainOf.size() || chainOf[start] != unassigned) continue;
    const std::size_t id = chains.size();
    Chain &chain = chains.emplace_back();
    std::vector<std::vector<std::size_t>> &neighbours = adjacency.emplace_back();
    chainOf[start] = id;
    localIndex[start] = 0;
    chain.atoms.push_back(start);
    chain.positions.push_back(atoms[start].position);
    neighbours.emplace_back();
    for (std::size_t head = 0; head < chain.atoms.size(); ++head)
    {
      const std::uint32_t current = chain.atoms[head];
      for (std::uint32_t partner : table.neighboursOf(current, components))
      {
        const double3 bond = simulationBox.applyPeriodicBoundaryConditions(atoms[partner].position -
                                                                           atoms[current].position);
        if (chainOf[partner] == unassigned)
        {
          chainOf[partner] = id;
          localIndex[partner] = chain.atoms.size();
          chain.atoms.push_back(partner);
          chain.positions.push_back(chain.positions[head] + bond);
          neighbours.emplace_back();
        }
        else if ((chain.positions[localIndex[partner]] - chain.positions[head] - bond).length() > 1.0)
        {
          chain.spanning = true;
        }
        neighbours[head].push_back(localIndex[partner]);
      }
    }
  }

  // The units of a chain are the pieces of molecules it holds: a molecule whose bond broke (an initiator split into
  // two radicals, a scission) is two pieces, which can end up in two chains or at the two ends of one. A piece is
  // the part of a molecule that its own bonds connect within the chain (the root of 'unitOf' per atom). Counted by
  // molecule, a chain that two pieces of one initiator end counted one unit too few and was taken for a cycle, and the
  // state of the groups of one piece was that of the chain of the other piece as well.
  std::vector<std::size_t> unitOf(table.layout.numberOfAtoms, unassigned);
  for (std::size_t id = 0; Chain &chain : chains)
  {
    std::vector<std::size_t> parent(chain.atoms.size());
    std::iota(parent.begin(), parent.end(), std::size_t{0});
    const auto root = [&](std::size_t k)
    {
      while (parent[k] != k) k = parent[k] = parent[parent[k]];
      return k;
    };
    for (std::size_t k = 0; k < chain.atoms.size(); ++k)
    {
      for (std::size_t next : adjacency[id][k])
      {
        if (table.layout.moleculeOf(chain.atoms[k]) == table.layout.moleculeOf(chain.atoms[next]))
        {
          parent[root(k)] = root(next);
        }
      }
    }
    for (std::size_t k = 0; k < chain.atoms.size(); ++k) unitOf[chain.atoms[k]] = chain.atoms[root(k)];
    ++id;
  }

  // per chain its units (an atom of the unit -> polymer bonds) and its last change
  std::vector<std::map<std::size_t, std::size_t>> units(chains.size());
  for (std::size_t id = 0; Chain &chain : chains)
  {
    for (std::uint32_t atom : chain.atoms) units[id].try_emplace(unitOf[atom], 0uz);
    chain.length = units[id].size();
    ++id;
  }
  for (const PolymerBond &bond : table.bonds)
  {
    const std::size_t a = table.layout.index(bond.a);
    const std::size_t id = chainOf[a];
    if (id == unassigned) continue;
    Chain &chain = chains[id];
    ++chain.polymerBonds;
    ++units[id][unitOf[a]];
    ++units[id][unitOf[table.layout.index(bond.b)]];
    chain.lastChange = std::max(chain.lastChange, bond.formedAtStep);
  }
  // a broken bond changed the chains of both its pieces
  for (const BrokenBond &bond : table.brokenBonds)
  {
    for (const std::size_t atom : {table.layout.index(bond.a), table.layout.index(bond.b)})
    {
      const std::size_t id = chainOf[atom];
      if (id != unassigned) chains[id].lastChange = std::max(chains[id].lastChange, bond.brokenAtStep);
    }
  }
  for (const Retyping &retyping : table.retypings)
  {
    const std::size_t id = chainOf[table.layout.index(retyping.atom)];
    if (id != unassigned) chains[id].lastChange = std::max(chains[id].lastChange, retyping.step);
  }

  for (std::size_t id = 0; Chain &chain : chains)
  {
    bool branched = false;
    for (const auto &[unit, bonds] : units[id])
    {
      branched = branched || bonds > 2;
      // the groups of the piece (those whose first atom it holds): the chain is living when one is active, and their
      // state changes are changes
      const PolymerAtomRef ref = table.layout.ref(static_cast<std::uint32_t>(unit));
      for (std::size_t g = 0; g < table.groupsPerMolecule[ref.componentId]; ++g)
      {
        const std::size_t group = table.groupIndex(ref.componentId, ref.moleculeIndex, g);
        if (unitOf[table.atomOf(group, 0, components)] != unit) continue;
        chain.living = chain.living || table.active(group);
        if (table.reacted(group)) chain.lastChange = std::max(chain.lastChange, table.groupStep[group]);
      }
    }
    // a tree of n pieces has n - 1 polymer bonds; one more (or one within a piece) closes a cycle
    chain.topology = chain.polymerBonds + 1 > chain.length ? Topology::Cyclic
                     : branched                            ? Topology::Branched
                                                           : Topology::Linear;

    // the end-to-end path: the farthest atom (in bonds) from the first, and the farthest atom from that one
    auto farthest = [&](std::size_t source) -> std::pair<std::size_t, std::size_t>
    {
      std::vector<std::size_t> distance(chain.atoms.size(), unassigned);
      std::vector<std::size_t> queue{source};
      distance[source] = 0;
      for (std::size_t head = 0; head < queue.size(); ++head)
      {
        for (std::size_t next : adjacency[id][queue[head]])
        {
          if (distance[next] != unassigned) continue;
          distance[next] = distance[queue[head]] + 1;
          queue.push_back(next);
        }
      }
      return {queue.back(), distance[queue.back()]};
    };
    const std::size_t end = farthest(0).first;
    const auto [other, bonds] = farthest(end);
    chain.contourBonds = bonds;
    chain.endToEnd = chain.positions[other] - chain.positions[end];
    ++id;
  }

  // only chains: a seed in a molecule without polymer bonds gives none
  std::erase_if(chains, [](const Chain &chain) { return chain.polymerBonds == 0; });
  return chains;
}

PropertyPolymerChains::Moments PropertyPolymerChains::measure(const Chain &chain, const ForceField &forceField,
                                                              std::span<const Atom> atoms, bool massWeighted)
{
  std::vector<double> weights(chain.atoms.size(), 1.0);
  double mass = 0.0;
  for (std::size_t i = 0; i < chain.atoms.size(); ++i)
  {
    const double m = forceField.pseudoAtoms[atoms[chain.atoms[i]].type].mass;
    mass += m;
    if (massWeighted) weights[i] = m;
  }
  double total = std::accumulate(weights.begin(), weights.end(), 0.0);
  if (total <= 0.0)
  {
    std::ranges::fill(weights, 1.0);
    total = static_cast<double>(weights.size());
  }
  for (double &w : weights) w /= total;

  const PropertyMoleculeShape::Descriptors d = PropertyMoleculeShape::computeDescriptors(chain.positions, weights);
  Moments values{};
  values[RadiusOfGyrationSquared] = d.radiusOfGyrationSquared;
  values[RadiusOfGyration] = std::sqrt(d.radiusOfGyrationSquared);
  values[Lambda1] = d.eigenvalues.x;
  values[Lambda2] = d.eigenvalues.y;
  values[Lambda3] = d.eigenvalues.z;
  values[Asphericity] = d.asphericity;
  values[Acylindricity] = d.acylindricity;
  values[ShapeAnisotropy] = d.shapeAnisotropy;
  values[Prolateness] = d.prolateness;
  values[EndToEndSquared] = double3::dot(chain.endToEnd, chain.endToEnd);
  values[EndToEnd] = chain.endToEnd.length();
  values[ContourBonds] = static_cast<double>(chain.contourBonds);
  values[NumberOfAtoms] = static_cast<double>(chain.atoms.size());
  values[Mass] = mass;
  return values;
}

std::size_t PropertyPolymerChains::conversionBin(double conversion) const
{
  const double x = std::clamp(conversion, 0.0, 1.0) * static_cast<double>(numberOfConversionBins);
  return std::min(static_cast<std::size_t>(x), numberOfConversionBins - 1);
}

void PropertyPolymerChains::sample(const PolymerLinkTable &table, const ForceField &forceField,
                                   const std::vector<Component> &components, const SimulationBox &simulationBox,
                                   std::span<const Atom> atoms, std::span<const std::uint32_t> seeds)
{
  if (!table.activated) return;
  ++numberOfSamplingEvents;
  const std::size_t bin = conversionBin(table.conversion());
  const std::uint64_t settled = table.relaxationSteps + lagSteps;
  for (const Chain &chain : findChains(table, components, simulationBox, atoms, seeds))
  {
    if (chain.length < minimumLength) continue;
    if (chain.spanning)
    {
      ++skippedSpanning;
      continue;
    }
    if (table.currentStep < chain.lastChange + settled)
    {
      ++skippedRamping;
      continue;
    }
    accumulators[Key{chain.length, chain.topology, chain.living, bin}].add(measure(chain, forceField, atoms,
                                                                                    massWeighted));
  }
}

std::map<std::pair<std::size_t, PropertyPolymerChains::Topology>, PropertyPolymerChains::Accumulator>
PropertyPolymerChains::pooled() const
{
  std::map<std::pair<std::size_t, Topology>, Accumulator> result{};
  for (const auto &[key, accumulator] : accumulators) result[{key.length, key.topology}] += accumulator;
  return result;
}

void PropertyPolymerChains::writeOutput(std::size_t systemId, std::size_t currentCycle, const PolymerLinkTable &table,
                                        bool force) const
{
  if (!force && currentCycle % writeEvery != 0uz) return;
  if (!table.activated) return;

  std::filesystem::create_directory("polymer_chains");

  // the columns shared by both text files, after the key columns
  auto header = [&](std::ostream &stream, std::size_t column)
  {
    std::print(stream, "# column {}: samples [-]\n", column++);
    for (std::size_t m = 0; m < NumberOfMoments; ++m)
    {
      std::print(stream, "# column {}: <{}> [{}]\n", column++, momentNames[m], momentUnits[m]);
      std::print(stream, "# column {}: standard deviation of {} [{}]\n", column++, momentNames[m], momentUnits[m]);
    }
    std::print(stream, "# column {}: <R^2>/<Rg^2> [-]\n", column);
  };
  auto row = [&](std::ostream &stream, const Accumulator &accumulator)
  {
    std::print(stream, " {}", accumulator.count);
    for (std::size_t m = 0; m < NumberOfMoments; ++m)
    {
      std::print(stream, " {:.6g} {:.6g}", accumulator.mean(static_cast<Moment>(m)),
                 accumulator.standardDeviation(static_cast<Moment>(m)));
    }
    const double rg2 = accumulator.mean(RadiusOfGyrationSquared);
    std::print(stream, " {:.6g}\n", rg2 > 0.0 ? accumulator.mean(EndToEndSquared) / rg2 : 0.0);
  };
  auto preamble = [&](std::ostream &stream)
  {
    std::print(stream, "# size and shape of the polymer chains per chain length (molecules per chain), measured during the\n");
    std::print(stream, "# polymerization once a chain has kept its topology for the force ramp ({} steps) plus {} steps\n",
               table.relaxationSteps, lagSteps);
    std::print(stream, "# sampling: {}{}, {} weights; step {}, conversion {:.6f}\n", samplingName(sampling),
               sampling == Sampling::Periodic ? std::format(" every {} MD steps", sampleEvery) : std::string{},
               massWeighted ? "mass" : "uniform bead", table.currentStep, table.conversion());
    std::print(stream, "# {} sampling events, {} chains skipped while ramping, {} spanning the periodic boundaries\n",
               numberOfSamplingEvents, skippedRamping, skippedSpanning);
    std::print(stream, "# the standard deviations are the widths of the distributions, not the errors of the means; combine\n");
    std::print(stream, "# independent runs from the JSON file for the errors (raspalib.polymer_chains)\n");
  };

  {
    std::ofstream stream(std::format("polymer_chains/polymer_chains.s{}.txt", systemId));
    preamble(stream);
    std::print(stream, "# pooled over living and dead chains and over the conversion\n");
    std::print(stream, "# column 1: chain length [molecules]\n");
    std::print(stream, "# column 2: topology (linear, branched, cyclic)\n");
    header(stream, 3);
    for (const auto &[key, accumulator] : pooled())
    {
      std::print(stream, "{} {}", key.first, topologyName(key.second));
      row(stream, accumulator);
    }
  }

  {
    std::ofstream stream(std::format("polymer_chains/polymer_chains_resolved.s{}.txt", systemId));
    preamble(stream);
    std::print(stream, "# column 1: chain length [molecules]\n");
    std::print(stream, "# column 2: topology (linear, branched, cyclic)\n");
    std::print(stream, "# column 3: living (1: a group of the chain is active) [-]\n");
    std::print(stream, "# column 4: lower edge of the conversion bin [-]\n");
    std::print(stream, "# column 5: upper edge of the conversion bin [-]\n");
    header(stream, 6);
    const double width = 1.0 / static_cast<double>(numberOfConversionBins);
    for (const auto &[key, accumulator] : accumulators)
    {
      std::print(stream, "{} {} {} {:.4f} {:.4f}", key.length, topologyName(key.topology), key.living ? 1 : 0,
                 static_cast<double>(key.conversionBin) * width, static_cast<double>(key.conversionBin + 1) * width);
      row(stream, accumulator);
    }
  }

  // the raw sums: runs are combined by adding them
  nlohmann::json json{};
  json["version"] = 1;
  json["system"] = systemId;
  json["step"] = table.currentStep;
  json["conversion"] = table.conversion();
  json["settings"] = {{"sampling", samplingName(sampling)},
                      {"sampleEvery", sampleEvery},
                      {"relaxationSteps", table.relaxationSteps},
                      {"lagSteps", lagSteps},
                      {"conversionBins", numberOfConversionBins},
                      {"minimumLength", minimumLength},
                      {"massWeighted", massWeighted}};
  json["samplingEvents"] = numberOfSamplingEvents;
  json["skippedRamping"] = skippedRamping;
  json["skippedSpanning"] = skippedSpanning;
  json["moments"] = std::vector<std::string>(momentNames.begin(), momentNames.end());
  json["units"] = std::vector<std::string>(momentUnits.begin(), momentUnits.end());
  nlohmann::json entries = nlohmann::json::array();
  for (const auto &[key, accumulator] : accumulators)
  {
    entries.push_back({{"length", key.length},
                       {"topology", topologyName(key.topology)},
                       {"living", key.living},
                       {"conversionBin", key.conversionBin},
                       {"count", accumulator.count},
                       {"sum", accumulator.sum},
                       {"sumOfSquares", accumulator.sumOfSquares}});
  }
  json["entries"] = entries;
  std::ofstream stream(std::format("polymer_chains/polymer_chains.s{}.json", systemId));
  stream << json.dump(1) << "\n";
}

std::string PropertyPolymerChains::printStatus() const
{
  std::uint64_t samples = 0;
  std::size_t longest = 0;
  std::set<std::size_t> lengths{};
  for (const auto &[key, accumulator] : accumulators)
  {
    samples += accumulator.count;
    longest = std::max(longest, key.length);
    lengths.insert(key.length);
  }
  return std::format(
      "    polymer chains:                  {} chain samples of {} lengths (longest {}), {} skipped while ramping, {} "
      "spanning\n",
      samples, lengths.size(), longest, skippedRamping, skippedSpanning);
}

std::string PropertyPolymerChains::printSettings() const
{
  std::ostringstream stream;
  std::print(stream, "Polymer chains:\n");
  std::print(stream, "    size and shape per chain length, measured {} after the force ramp plus {} steps\n",
             sampling == Sampling::Periodic ? std::format("every {} MD steps", sampleEvery)
                                            : std::string("before every reaction that changes a chain (and at the end)"),
             lagSteps);
  std::print(stream, "    weights: {}, minimum length: {} molecules, conversion bins: {}\n",
             massWeighted ? "mass" : "uniform bead", minimumLength, numberOfConversionBins);
  std::print(stream, "    write every: {}\n", writeEvery);
  std::print(stream, "\n");
  return stream.str();
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PropertyPolymerChains &p)
{
  archive << p.versionNumber;
  archive << p.sampleEvery;
  archive << p.writeEvery;
  archive << p.sampling;
  archive << p.lagSteps;
  archive << p.numberOfConversionBins;
  archive << p.minimumLength;
  archive << p.massWeighted;
  archive << p.accumulators;
  archive << p.numberOfSamplingEvents;
  archive << p.skippedRamping;
  archive << p.skippedSpanning;
#if DEBUG_ARCHIVE
  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif
  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyPolymerChains &p)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > p.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'PropertyPolymerChains' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }
  archive >> p.sampleEvery;
  archive >> p.writeEvery;
  archive >> p.sampling;
  archive >> p.lagSteps;
  archive >> p.numberOfConversionBins;
  archive >> p.minimumLength;
  archive >> p.massWeighted;
  archive >> p.accumulators;
  archive >> p.numberOfSamplingEvents;
  archive >> p.skippedRamping;
  archive >> p.skippedSpanning;
#if DEBUG_ARCHIVE
  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("PropertyPolymerChains: Error in binary restart\n"));
  }
#endif
  return archive;
}
