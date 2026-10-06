module;

module intra_molecular_exclusions;

import std;

import archive;
import connectivity_table;

IntraMolecularExclusions IntraMolecularExclusions::allPairs(std::size_t numberOfAtoms)
{
  std::vector<std::array<std::size_t, 2>> excluded{};
  excluded.reserve(numberOfAtoms * (numberOfAtoms - 1) / 2);
  for (std::size_t i = 0; i + 1 < numberOfAtoms; ++i)
  {
    for (std::size_t j = i + 1; j < numberOfAtoms; ++j)
    {
      excluded.push_back({i, j});
    }
  }
  IntraMolecularExclusions exclusions = fromPairs(numberOfAtoms, std::move(excluded));
  exclusions.allExcluded = true;
  return exclusions;
}

IntraMolecularExclusions IntraMolecularExclusions::fromPairs(std::size_t numberOfAtoms,
                                                             std::vector<std::array<std::size_t, 2>> excluded)
{
  IntraMolecularExclusions exclusions{};
  exclusions.numberOfAtoms = numberOfAtoms;

  for (std::array<std::size_t, 2> &pair : excluded)
  {
    if (pair[0] >= numberOfAtoms || pair[1] >= numberOfAtoms || pair[0] == pair[1])
    {
      throw std::runtime_error(
          std::format("[IntraMolecularExclusions]: invalid excluded pair [{}, {}] for {} atoms\n", pair[0], pair[1],
                      numberOfAtoms));
    }
    if (pair[1] < pair[0]) std::swap(pair[0], pair[1]);
  }
  std::ranges::sort(excluded);
  const auto duplicates = std::ranges::unique(excluded);
  excluded.erase(duplicates.begin(), duplicates.end());

  exclusions.pairs.reserve(excluded.size());
  std::vector<std::vector<std::uint32_t>> partnersPerAtom(numberOfAtoms);
  for (const std::array<std::size_t, 2> &pair : excluded)
  {
    const std::uint32_t a = static_cast<std::uint32_t>(pair[0]);
    const std::uint32_t b = static_cast<std::uint32_t>(pair[1]);
    exclusions.pairs.push_back({a, b});
    partnersPerAtom[a].push_back(b);
    partnersPerAtom[b].push_back(a);
  }

  exclusions.offsets.assign(numberOfAtoms + 1, 0);
  exclusions.partners.reserve(2 * excluded.size());
  for (std::size_t i = 0; i < numberOfAtoms; ++i)
  {
    std::ranges::sort(partnersPerAtom[i]);
    exclusions.partners.insert(exclusions.partners.end(), partnersPerAtom[i].begin(), partnersPerAtom[i].end());
    exclusions.offsets[i + 1] = static_cast<std::uint32_t>(exclusions.partners.size());
  }
  exclusions.pairListOffsets = exclusions.offsets;
  exclusions.pairListPartners = exclusions.partners;

  return exclusions;
}

void IntraMolecularExclusions::setScaledPairs(std::vector<ScaledPair> scaled)
{
  for (ScaledPair &pair : scaled)
  {
    if (pair.atomA >= numberOfAtoms || pair.atomB >= numberOfAtoms || pair.atomA == pair.atomB)
    {
      throw std::runtime_error(std::format("[IntraMolecularExclusions]: invalid scaled pair [{}, {}] for {} atoms\n",
                                           pair.atomA, pair.atomB, numberOfAtoms));
    }
    if (pair.atomB < pair.atomA) std::swap(pair.atomA, pair.atomB);
    if (isExcluded(pair.atomA, pair.atomB))
    {
      throw std::runtime_error(
          std::format("[IntraMolecularExclusions]: scaled pair [{}, {}] is an excluded pair\n", pair.atomA, pair.atomB));
    }
  }
  std::ranges::sort(scaled, {}, [](const ScaledPair &p) { return std::array<std::uint32_t, 2>{p.atomA, p.atomB}; });
  const auto duplicates = std::ranges::unique(scaled, {}, [](const ScaledPair &p) { return std::array<std::uint32_t, 2>{p.atomA, p.atomB}; });
  scaled.erase(duplicates.begin(), duplicates.end());
  scaledPairs = std::move(scaled);

  std::vector<std::vector<std::uint32_t>> partnersPerAtom(numberOfAtoms);
  for (std::size_t i = 0; i < numberOfAtoms; ++i)
  {
    const std::span<const std::uint32_t> excluded = partnersOf(i);
    partnersPerAtom[i].assign(excluded.begin(), excluded.end());
  }
  for (const ScaledPair &pair : scaledPairs)
  {
    partnersPerAtom[pair.atomA].push_back(pair.atomB);
    partnersPerAtom[pair.atomB].push_back(pair.atomA);
  }

  pairListOffsets.assign(numberOfAtoms + 1, 0);
  pairListPartners.clear();
  for (std::size_t i = 0; i < numberOfAtoms; ++i)
  {
    std::ranges::sort(partnersPerAtom[i]);
    const auto duplicates = std::ranges::unique(partnersPerAtom[i]);
    partnersPerAtom[i].erase(duplicates.begin(), duplicates.end());
    pairListPartners.insert(pairListPartners.end(), partnersPerAtom[i].begin(), partnersPerAtom[i].end());
    pairListOffsets[i + 1] = static_cast<std::uint32_t>(pairListPartners.size());
  }
}

std::vector<std::array<std::size_t, 2>> IntraMolecularExclusions::topologicalPairs(
    const ConnectivityTable &connectivity)
{
  std::vector<std::array<std::size_t, 2>> pairs{};
  for (const std::array<std::size_t, 2> &bond : connectivity.findAllBonds())
  {
    pairs.push_back({std::min(bond[0], bond[1]), std::max(bond[0], bond[1])});
  }
  for (const std::array<std::size_t, 3> &bend : connectivity.findAllBends())
  {
    pairs.push_back({std::min(bend[0], bend[2]), std::max(bend[0], bend[2])});
  }
  std::ranges::sort(pairs);
  const auto duplicates = std::ranges::unique(pairs);
  pairs.erase(duplicates.begin(), duplicates.end());
  return pairs;
}

std::vector<std::array<std::size_t, 2>> IntraMolecularExclusions::pairs14(const ConnectivityTable &connectivity)
{
  std::vector<std::array<std::size_t, 2>> topological = topologicalPairs(connectivity);
  std::vector<std::array<std::size_t, 2>> pairs{};
  for (const std::array<std::size_t, 4> &torsion : connectivity.findAllTorsions())
  {
    const std::array<std::size_t, 2> pair{std::min(torsion[0], torsion[3]), std::max(torsion[0], torsion[3])};
    if (pair[0] == pair[1]) continue;
    if (std::ranges::binary_search(topological, pair)) continue;
    pairs.push_back(pair);
  }
  std::ranges::sort(pairs);
  const auto duplicates = std::ranges::unique(pairs);
  pairs.erase(duplicates.begin(), duplicates.end());
  return pairs;
}

std::string IntraMolecularExclusions::printStatus() const
{
  std::ostringstream stream;
  std::print(stream, "    Non-bonded exclusions: {} pairs excluded, {} pairs interacting, {} pairs scaled\n",
             pairs.size(), numberOfNonExcludedPairs(), scaledPairs.size());
  for (const ScaledPair &pair : scaledPairs)
  {
    std::print(stream, "        scaled pair ({}, {}): vdW {:g}, Coulomb {:g}{}\n", pair.atomA, pair.atomB,
               pair.scalingVDW, pair.scalingCoulomb, pair.pair14 ? ", 1-4 parameters" : "");
  }
  return stream.str();
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const IntraMolecularExclusions &e)
{
  archive << e.versionNumber;

  archive << e.numberOfAtoms;
  archive << e.allExcluded;
  archive << e.pairs;
  archive << e.offsets;
  archive << e.partners;
  archive << e.pairListOffsets;
  archive << e.pairListPartners;

  archive << e.scaledPairs.size();
  for (const IntraMolecularExclusions::ScaledPair &pair : e.scaledPairs)
  {
    archive << pair.atomA;
    archive << pair.atomB;
    archive << pair.scalingVDW;
    archive << pair.scalingCoulomb;
    archive << pair.pair14;
  }

#if DEBUG_ARCHIVE
  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, IntraMolecularExclusions &e)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > e.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'IntraMolecularExclusions' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> e.numberOfAtoms;
  archive >> e.allExcluded;
  archive >> e.pairs;
  archive >> e.offsets;
  archive >> e.partners;
  archive >> e.pairListOffsets;
  archive >> e.pairListPartners;

  std::size_t numberOfScaledPairs{};
  archive >> numberOfScaledPairs;
  e.scaledPairs.resize(numberOfScaledPairs);
  for (IntraMolecularExclusions::ScaledPair &pair : e.scaledPairs)
  {
    archive >> pair.atomA;
    archive >> pair.atomB;
    archive >> pair.scalingVDW;
    archive >> pair.scalingCoulomb;
    if (versionNumber >= 3) archive >> pair.pair14;
  }

#if DEBUG_ARCHIVE
  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("IntraMolecularExclusions: Error in binary restart\n"));
  }
#endif

  return archive;
}
