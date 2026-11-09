// Generates the CIF test structures in testdata/spacegroups: for every Hall setting the hardest few of a hundred
// candidate structures whose symmetry is exactly that setting, then a few with atoms on special positions and a few
// with noisy coordinates.
// Usage: generate_symmetrykit_testdata [output directory] [structures kept per setting] [candidates per setting]
//
// A candidate puts two to four orbits of random general positions, each of a different element, in a cell whose
// metric is a random one symmetrised over the rotations of the setting. A single orbit can have more symmetry than
// its group (one atom per cell is centrosymmetric, a lone triangle of P 3 has a mirror plane), but two generic
// orbits of different elements only share the group's own operations. Half of the candidates are drawn to be hard:
// their metric is a cubic one nudged by a little, so that the lattice is almost but not quite more symmetric than
// the group, and their atoms sit just off special positions, so that symmetry images almost but not quite coincide.
// Every candidate is checked with the space-group search and dropped when its symmetry is not the setting's (the
// lattice symmetry is checked separately to be that of the crystal family). The candidates are ranked by how close
// their metric is to a higher symmetry, how close their atoms come to each other and how many atoms they have, and
// the hardest are written.
//
// Two further kinds are written after those, half as many of each: 'special' candidates add one or two orbits of
// atoms on special positions (fixed points of one or two of the setting's operations, so with a smaller orbit)
// before the general ones, and 'noisy' candidates displace every atom by up to a tenth of a symmetry precision of
// 0.05 Angstrom and jolt the cell by a little, and are to be identified at that precision, which is written to the
// file as _raspa_symmetry_precision.

import std;

import int3;
import int3x3;
import double3;
import double3x3;
import skdefinitions;
import skrotationmatrix;
import skseitzintegermatrix;
import skintegersymmetryoperationset;
import skpointgroup;
import skspacegroup;
import skspacegroupsetting;
import skspacegroupdatabase;
import sksymmetrycell;

namespace
{
constexpr std::array<const char*, 52> elements{
    "H",  "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Na", "Mg", "Al", "Si", "P",  "S",  "Cl", "K",  "Ca", "Ti",
    "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Rb", "Sr", "Zr", "Nb", "Mo",
    "Ag", "Cd", "In", "Sn", "Sb", "I",  "Cs", "Ba", "La", "Ce", "W",  "Pt", "Au", "Pb", "Bi", "U"};

// the operations of a setting in a fixed order, so that the output does not depend on hashing
std::vector<SKSeitzIntegerMatrix> sortedOperations(const SKSpaceGroupSetting& setting)
{
  const auto& set = setting.fullSeitzMatrices().operations;
  std::vector<SKSeitzIntegerMatrix> operations(set.begin(), set.end());
  std::sort(operations.begin(), operations.end(),
            [](const SKSeitzIntegerMatrix& a, const SKSeitzIntegerMatrix& b)
            {
              const int3x3& p = a.rotation.int3x3_m;
              const int3x3& q = b.rotation.int3x3_m;
              const int3 ta = a.normalizedTranslation();
              const int3 tb = b.normalizedTranslation();
              return std::tie(p.m11, p.m12, p.m13, p.m21, p.m22, p.m23, p.m31, p.m32, p.m33, ta.x, ta.y, ta.z) <
                     std::tie(q.m11, q.m12, q.m13, q.m21, q.m22, q.m23, q.m31, q.m32, q.m33, tb.x, tb.y, tb.z);
            });
  return operations;
}

// the number of rotations of the lattice of the setting's crystal family (a trigonal group in hexagonal axes has
// the hexagonal lattice, in rhombohedral axes the rhombohedral one)
std::size_t latticeOrder(const SKSpaceGroupSetting& setting)
{
  switch (SKPointGroup::pointGroupData[setting.pointGroupNumber()].holohedry())
  {
    case Holohedry::triclinic:
      return 2;
    case Holohedry::monoclinic:
      return 4;
    case Holohedry::orthorhombic:
      return 8;
    case Holohedry::tetragonal:
      return 16;
    case Holohedry::trigonal:
      return setting.qualifier() == "R" ? 12 : 24;
    case Holohedry::hexagonal:
      return 24;
    case Holohedry::cubic:
      return 48;
    default:
      return 0;
  }
}

enum class Kind
{
  plain,
  hard,
  special,
  noisy
};
constexpr const char* kindNames[] = {"plain", "hard", "special", "noisy"};

struct Structure
{
  double3x3 cell;
  std::vector<std::tuple<double3, std::size_t, double>> atoms;  // position, element index, occupancy
  double metricGap{};         // how far the metric is from a higher symmetry (relative lengths, angles over 90)
  double shortestDistance{};  // the closest two atoms come [Angstrom]
  Kind kind{};
  double precision{};  // the symmetry precision the structure is to be identified at [Angstrom]
};

// the exact image of a point under an operation, without wrapping into the cell
double3 exactImage(const SKSeitzIntegerMatrix& operation, const double3& point)
{
  return double3x3(operation.rotation.int3x3_m) * point + double3(operation.translation) / 24.0;
}

// the operations whose cyclic group has a fixed point: those whose power is the identity, not a lattice
// translation as for a screw or glide (a setting with none, like P 21 21 21, has no special positions)
std::vector<const SKSeitzIntegerMatrix*> fixedPointOperations(const std::vector<SKSeitzIntegerMatrix>& operations)
{
  std::vector<const SKSeitzIntegerMatrix*> closed;
  for (const SKSeitzIntegerMatrix& operation : operations)
  {
    if (operation.rotation == SKRotationMatrix::identity) continue;
    const double3 probe(0.1234, 0.2345, 0.3456);
    double3 image = probe;
    for (int power = 0; power < 6; ++power)
    {
      image = exactImage(operation, image);
      if ((image - probe).length_squared() < 1.0e-18)
      {
        closed.push_back(&operation);
        break;
      }
    }
  }
  return closed;
}

// a point fixed by one or two random operations of the setting, so on a special position: the average of a random
// point over the cyclic group of an operation is fixed by it, and alternately averaging over two operations
// converges to a point fixed by both when their fixed sets intersect
std::optional<double3> specialSite(const std::vector<SKSeitzIntegerMatrix>& operations, std::mt19937_64& rng)
{
  std::uniform_real_distribution<double> fraction(0.0, 1.0);
  std::uniform_int_distribution<int> coin(0, 1);
  const std::vector<const SKSeitzIntegerMatrix*> closed = fixedPointOperations(operations);
  if (closed.empty()) return std::nullopt;

  std::uniform_int_distribution<std::size_t> pickClosed(0, closed.size() - 1);
  std::vector<const SKSeitzIntegerMatrix*> chosen{closed[pickClosed(rng)]};
  if (coin(rng)) chosen.push_back(closed[pickClosed(rng)]);

  double3 site(fraction(rng), fraction(rng), fraction(rng));
  for (int iteration = 0; iteration < 60; ++iteration)
  {
    for (const SKSeitzIntegerMatrix* operation : chosen)
    {
      double3 sum = site;
      double3 image = site;
      std::size_t count = 1;
      for (int power = 1; power < 6; ++power)
      {
        image = exactImage(*operation, image);
        if ((image - site).length_squared() < 1.0e-18) break;
        sum += image;
        ++count;
      }
      site = sum / double(count);
    }
  }
  for (const SKSeitzIntegerMatrix* operation : chosen)
  {
    const double3 d = exactImage(*operation, site) - site;
    if ((d - double3(std::rint(d.x), std::rint(d.y), std::rint(d.z))).length_squared() > 1.0e-16) return std::nullopt;
  }
  return double3::fract(site);
}

// how close the cell is to a higher symmetry than the group forces: the smallest relative difference between two
// lengths that are not equal and the smallest distance of an angle from 60, 90 or 120 degrees that is not exactly
// there; the smaller, the harder the structure is to assign
double metricGap(const SKSymmetryCell& cell)
{
  const double lengths[3] = {cell.a(), cell.b(), cell.c()};
  const double angles[3] = {cell.alpha() * 180.0 / std::numbers::pi, cell.beta() * 180.0 / std::numbers::pi,
                            cell.gamma() * 180.0 / std::numbers::pi};
  double gap = 1.0;
  for (std::size_t i = 0; i < 3; ++i)
    for (std::size_t j = i + 1; j < 3; ++j)
    {
      const double relative = std::fabs(lengths[i] - lengths[j]) / std::max(lengths[i], lengths[j]);
      if (relative > 1.0e-6) gap = std::min(gap, relative);
    }
  for (const double angle : angles)
    for (const double special : {60.0, 90.0, 120.0})
    {
      const double difference = std::fabs(angle - special) / 90.0;
      if (difference > 1.0e-8) gap = std::min(gap, difference);
    }
  return gap;
}

std::optional<Structure> generate(const SKSpaceGroupSetting& setting,
                                  const std::vector<SKSeitzIntegerMatrix>& operations, std::mt19937_64& rng, Kind kind)
{
  // the cell has to hold up to four orbits of 'order' atoms at the spacing: its edge grows with the order
  const bool hard = kind == Kind::hard;
  const double precision = kind == Kind::noisy ? 0.05 : 1.0e-5;
  const double spacing = hard ? 0.5 : 0.8;
  std::uniform_real_distribution<double> noise(-precision / 10.0, precision / 10.0);
  const double shortestEdge = std::max(5.0, std::cbrt(3.0 * 4.0 * double(operations.size())));
  std::uniform_real_distribution<double> length(shortestEdge, std::max(25.0, shortestEdge + 10.0));
  std::uniform_real_distribution<double> angle(60.0, 120.0);
  std::uniform_real_distribution<double> fraction(0.0, 1.0);
  std::uniform_int_distribution<std::size_t> orbitCount(2, 4);
  std::uniform_int_distribution<std::size_t> element(0, elements.size() - 1);
  std::uniform_real_distribution<double> lengthNudge(0.02, 0.5);
  std::uniform_real_distribution<double> angleNudge(0.05, 2.0);
  std::uniform_real_distribution<double> positionNudge(0.005, 0.05);
  std::uniform_int_distribution<int> coin(0, 1);
  std::uniform_int_distribution<std::size_t> specialIndex(0, 6);
  constexpr double specialValues[7] = {0.0, 0.25, 0.5, 0.75, 1.0 / 3.0, 2.0 / 3.0, 0.125};
  const auto sign = [&] { return coin(rng) ? 1.0 : -1.0; };

  for (int attempt = 0; attempt < 1000; ++attempt)
  {
    double a, b, c, alpha, beta, gamma;
    if (hard)
    {
      // a cubic metric nudged by a little, so that the lattice is almost of a higher symmetry than the group's
      a = length(rng);
      b = a + sign() * lengthNudge(rng);
      c = a + sign() * lengthNudge(rng);
      alpha = 90.0 + sign() * angleNudge(rng);
      beta = 90.0 + sign() * angleNudge(rng);
      gamma = 90.0 + sign() * angleNudge(rng);
    }
    else
    {
      // a generic cell: angles away from 90, lengths apart
      a = length(rng), b = length(rng), c = length(rng);
      alpha = angle(rng), beta = angle(rng), gamma = angle(rng);
      if (std::fabs(a - b) < 0.5 || std::fabs(b - c) < 0.5 || std::fabs(a - c) < 0.5) continue;
      if (std::fabs(alpha - 90.0) < 3.0 || std::fabs(beta - 90.0) < 3.0 || std::fabs(gamma - 90.0) < 3.0) continue;
    }

    // symmetrised over the setting's rotations: the metric then has the lattice symmetry of the crystal family
    SKSymmetryCell seed(a, b, c, alpha, beta, gamma);
    const double3x3 metric = seed.metricTensor();
    double3 columns[3] = {double3(0.0, 0.0, 0.0), double3(0.0, 0.0, 0.0), double3(0.0, 0.0, 0.0)};
    for (const SKSeitzIntegerMatrix& operation : operations)
    {
      const double3x3 rotation(operation.rotation.int3x3_m);
      const double3x3 transformed = rotation.transpose() * metric * rotation;
      for (std::size_t j = 0; j < 3; ++j) columns[j] += transformed[j];
    }
    const double count = double(operations.size());
    const double3x3 symmetrisedMetric(columns[0] / count, columns[1] / count, columns[2] / count);
    double3x3 cell = SKSymmetryCell::createFromMetricTensor(symmetrisedMetric).unitCell();
    if (kind == Kind::noisy)
    {
      // the cell jolted by a little (a tenth of the precision), so that its symmetry is only approximate
      const SKSymmetryCell exact = SKSymmetryCell::createFromUnitCell(cell);
      cell = SKSymmetryCell(exact.a() + noise(rng), exact.b() + noise(rng), exact.c() + noise(rng),
                            (exact.alpha() + noise(rng) / exact.a()) * 180.0 / std::numbers::pi,
                            (exact.beta() + noise(rng) / exact.b()) * 180.0 / std::numbers::pi,
                            (exact.gamma() + noise(rng) / exact.c()) * 180.0 / std::numbers::pi)
                 .unitCell();
    }
    if (!(cell.determinant() > 50.0)) continue;
    if (SKSymmetryCell::findLatticeSymmetry(cell, std::max(precision, 1.0e-4)).rotations().size() !=
        latticeOrder(setting))
      continue;

    // orbits of general positions, each of a different element, apart from each other
    Structure structure{cell, {}, 0.0, 0.0, kind, precision};
    const std::size_t orbits = orbitCount(rng);
    // special candidates (and half of the noisy ones) start with one or two orbits on special positions
    const std::size_t specialOrbits =
        kind == Kind::special || (kind == Kind::noisy && coin(rng)) ? std::size_t(1 + coin(rng)) : 0;
    bool valid = true;
    std::vector<std::size_t> usedTypes;
    for (std::size_t orbit = 0; orbit < orbits + specialOrbits && valid; ++orbit)
    {
      const bool special = orbit < specialOrbits;
      double3 site(fraction(rng), fraction(rng), fraction(rng));
      if (special)
      {
        const std::optional<double3> fixed = specialSite(operations, rng);
        if (!fixed)
        {
          valid = false;
          break;
        }
        site = *fixed;
      }
      if (hard)
      {
        // just off a special position, so that symmetry images almost coincide
        site = double3(specialValues[specialIndex(rng)] + sign() * positionNudge(rng),
                       specialValues[specialIndex(rng)] + sign() * positionNudge(rng),
                       specialValues[specialIndex(rng)] + sign() * positionNudge(rng));
      }
      std::size_t type = element(rng);
      while (std::find(usedTypes.begin(), usedTypes.end(), type) != usedTypes.end()) type = element(rng);
      usedTypes.push_back(type);

      std::vector<double3> images;
      for (const SKSeitzIntegerMatrix& operation : operations)
      {
        const double3 image = double3::fract(operation * site);
        const bool seen = std::any_of(
            images.begin(), images.end(),
            [&](const double3& other)
            {
              const double3 d = image - other;
              return (d - double3(std::rint(d.x), std::rint(d.y), std::rint(d.z))).length_squared() < 1.0e-12;
            });
        if (!seen) images.push_back(image);
      }
      if (special ? images.size() >= operations.size() : images.size() != operations.size())
      {
        valid = false;  // a general position where a special one was meant, or the other way round
        break;
      }
      for (const double3& image : images) structure.atoms.emplace_back(image, type, 1.0);
    }
    if (!valid) continue;

    if (kind == Kind::noisy)
    {
      // every atom displaced by a little (a tenth of the precision), so that the symmetry is only approximate
      for (auto& atom : structure.atoms)
      {
        std::get<0>(atom) =
            double3::fract(std::get<0>(atom) + cell.inverse() * double3(noise(rng), noise(rng), noise(rng)));
      }
    }

    // every pair of atoms apart; the closest approach is kept for the ranking
    structure.shortestDistance = std::numeric_limits<double>::max();
    for (std::size_t i = 0; i < structure.atoms.size() && valid; ++i)
      for (std::size_t j = i + 1; j < structure.atoms.size(); ++j)
      {
        double3 d = std::get<0>(structure.atoms[i]) - std::get<0>(structure.atoms[j]);
        d = d - double3(std::rint(d.x), std::rint(d.y), std::rint(d.z));
        const double distance = (cell * d).length();
        structure.shortestDistance = std::min(structure.shortestDistance, distance);
        if (distance < spacing)
        {
          valid = false;
          break;
        }
      }
    if (!valid) continue;
    structure.metricGap = metricGap(SKSymmetryCell::createFromUnitCell(cell));

    // the symmetry of the structure is the setting's group, no more and no less, both by element and with the
    // elements disregarded (as the search with partial occupancies allowed does: two atoms of P 1 are then
    // centrosymmetric)
    for (const bool allowPartialOccupancies : {false, true})
    {
      const std::optional<SKSpaceGroup::FoundSpaceGroupInfo> found =
          SKSpaceGroup::findSpaceGroup(cell, structure.atoms, allowPartialOccupancies, precision);
      if (!found || SKSpaceGroup(found->HallNumber).spaceGroupSetting().number() != setting.number()) valid = false;
    }
    if (!valid) continue;
    return structure;
  }
  return std::nullopt;
}

// the operations are not listed: the Hall symbol determines them, and listing them made the files five times larger
void writeCIF(const std::filesystem::path& path, const SKSpaceGroupSetting& setting, const Structure& structure,
              std::size_t hallNumber, std::size_t index)
{
  const SKSymmetryCell cell = SKSymmetryCell::createFromUnitCell(structure.cell);
  std::ofstream file(path);
  file << std::format("data_hall_{:03}_{:03}\n", hallNumber, index);
  file << std::format("# generated test structure: space group {} ({}), Hall setting {}\n", setting.number(),
                      setting.HMString(), hallNumber);
  file << std::format("# {} candidate: metric gap {:.3e}, closest atoms {:.3f} Angstrom, {} atoms\n",
                      kindNames[std::size_t(structure.kind)], structure.metricGap, structure.shortestDistance,
                      structure.atoms.size());
  file << std::format("_raspa_symmetry_precision {:.6g}\n", structure.precision);
  file << std::format("_raspa_hall_number {}\n", hallNumber);
  file << std::format("_raspa_point_group_number {}\n", setting.pointGroupNumber());
  file << std::format("_space_group_IT_number {}\n", setting.number());
  file << std::format("_symmetry_Int_Tables_number {}\n", setting.number());
  file << std::format("_symmetry_space_group_name_H-M '{}'\n", setting.HMString());
  file << std::format("_symmetry_space_group_name_Hall '{}'\n", setting.HallString());
  file << std::format("_cell_length_a {:.8f}\n_cell_length_b {:.8f}\n_cell_length_c {:.8f}\n", cell.a(), cell.b(),
                      cell.c());
  file << std::format("_cell_angle_alpha {:.8f}\n_cell_angle_beta {:.8f}\n_cell_angle_gamma {:.8f}\n",
                      cell.alpha() * 180.0 / std::numbers::pi, cell.beta() * 180.0 / std::numbers::pi,
                      cell.gamma() * 180.0 / std::numbers::pi);
  file << "loop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n"
          "_atom_site_fract_z\n_atom_site_occupancy\n";
  std::size_t label = 1;
  for (const auto& [position, type, occupancy] : structure.atoms)
  {
    file << std::format("{}{} {} {:.8f} {:.8f} {:.8f} {:.4f}\n", elements[type], label++, elements[type], position.x,
                        position.y, position.z, occupancy);
  }
}
}  // namespace

int main(int argc, char* argv[])
{
  const std::filesystem::path outputDirectory = argc > 1 ? argv[1] : "testdata/spacegroups";
  const std::size_t keptPerSetting = argc > 2 ? std::stoul(argv[2]) : 10;
  const std::size_t candidatesPerSetting = argc > 3 ? std::stoul(argv[3]) : 100;

  std::size_t written = 0;
  for (std::size_t hallNumber = 1; hallNumber <= 530; ++hallNumber)
  {
    const SKSpaceGroupSetting& setting = SKSpaceGroupDataBase::spaceGroupData[hallNumber];
    const std::vector<SKSeitzIntegerMatrix> operations = sortedOperations(setting);
    const std::filesystem::path directory = outputDirectory / std::format("{:03}", setting.number());
    std::filesystem::create_directories(directory);

    // the plain and hard candidates in one ranking, then the special and the noisy ones in their own, half as many
    const std::array<std::tuple<std::vector<Kind>, std::size_t, std::size_t>, 3> rounds{
        std::tuple{std::vector<Kind>{Kind::plain, Kind::hard}, candidatesPerSetting, keptPerSetting},
        std::tuple{std::vector<Kind>{Kind::special}, candidatesPerSetting / 2, keptPerSetting / 2},
        std::tuple{std::vector<Kind>{Kind::noisy}, candidatesPerSetting / 2, keptPerSetting / 2}};
    std::size_t seed = 0;
    std::size_t fileIndex = 0;
    for (const auto& [kinds, candidateCount, keptCount] : rounds)
    {
      if (kinds.front() == Kind::special && fixedPointOperations(operations).empty()) continue;
      std::vector<Structure> candidates;
      for (std::size_t index = 0; index < candidateCount; ++index, ++seed)
      {
        std::mt19937_64 rng(hallNumber * 1000 + seed);
        if (std::optional<Structure> structure = generate(setting, operations, rng, kinds[index % kinds.size()]))
        {
          candidates.push_back(std::move(*structure));
        }
      }
      if (candidates.size() < keptCount)
      {
        std::cerr << std::format("only {} of {} {} candidates valid for Hall setting {} ({})\n", candidates.size(),
                                 candidateCount, kindNames[std::size_t(kinds.front())], hallNumber, setting.HMString());
      }

      // rank: the sum of the ranks by metric gap (small first), closest approach (small first) and atom count (large
      // first); the lowest sums are the hardest
      std::vector<double> score(candidates.size(), 0.0);
      const auto addRank = [&](auto key)
      {
        std::vector<std::size_t> order(candidates.size());
        std::iota(order.begin(), order.end(), 0uz);
        std::sort(order.begin(), order.end(),
                  [&](std::size_t x, std::size_t y) { return key(candidates[x]) < key(candidates[y]); });
        for (std::size_t rank = 0; rank < order.size(); ++rank) score[order[rank]] += double(rank);
      };
      addRank([](const Structure& structure) { return structure.metricGap; });
      addRank([](const Structure& structure) { return structure.shortestDistance; });
      addRank([](const Structure& structure) { return -double(structure.atoms.size()); });
      std::vector<std::size_t> order(candidates.size());
      std::iota(order.begin(), order.end(), 0uz);
      std::sort(order.begin(), order.end(), [&](std::size_t x, std::size_t y) { return score[x] < score[y]; });

      for (std::size_t index = 0; index < std::min(keptCount, order.size()); ++index, ++fileIndex)
      {
        writeCIF(directory / std::format("hall_{:03}_{:03}.cif", hallNumber, fileIndex), setting,
                 candidates[order[index]], hallNumber, fileIndex);
        ++written;
      }
    }
  }
  std::cout << std::format("wrote {} structures to {}\n", written, outputDirectory.string());
  return 0;
}
