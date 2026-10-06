#include <gtest/gtest.h>

import std;

import amber_prmtop_reader;
import atom;
import component;
import double3;
import forcefield;
import json;
import running_energy;
import simulationbox;
import system;
import units;

// A hand-built AMBER topology: N-methylacetamide with ff14SB-like parameters (plus a periodicity-4 term, an
// off-phase term and two impropers), three TIP3P waters with the SETTLE-style O-H/O-H/H-H bond topology, Na+
// and Cl- (with a non-Lorentz-Berthelot Na+/Cl- pair). The reference energies are evaluated here straight from
// the prmtop arrays with the AMBER functional forms; RASPA must reproduce them with the water as rigid bodies.
namespace
{
const char *prmtopText = R"PRMTOP(%VERSION  VERSION_STAMP = V0001.000  DATE = 01/01/26  00:00:00
%FLAG TITLE
%FORMAT(20a4)
NMA + 3 TIP3P + NaCl (RASPA test)
%FLAG POINTERS
%FORMAT(10I8)
      23      11      16       4      14       4      20       2       0       0
      60       7       4       4       2       9      10       9      11       0
       0       0       0       0       0       0       0       0       3       0
       0
%FLAG ATOM_NAME
%FORMAT(20a4)
CH3 HH31HH32HH33C   O   N   H   CH3 HH31HH32HH33O   H1  H2  O   H1  H2  O   H1  
H2  Na+ Cl- 
%FLAG CHARGE
%FORMAT(5E16.8)
 -6.67300626E+00  2.04636429E+00  2.04636429E+00  2.04636429E+00  1.08823576E+01
 -1.03484442E+01 -7.57501011E+00  4.95464337E+00 -2.71512270E+00  1.77849648E+00
  1.77849648E+00  1.77849648E+00 -1.51973982E+01  7.59869910E+00  7.59869910E+00
 -1.51973982E+01  7.59869910E+00  7.59869910E+00 -1.51973982E+01  7.59869910E+00
  7.59869910E+00  1.82223000E+01 -1.82223000E+01
%FLAG ATOMIC_NUMBER
%FORMAT(10I8)
       6       1       1       1       6       8       7       1       6       1
       1       1       8       1       1       8       1       1       8       1
       1      11      17
%FLAG MASS
%FORMAT(5E16.8)
  1.20100000E+01  1.00800000E+00  1.00800000E+00  1.00800000E+00  1.20100000E+01
  1.60000000E+01  1.40100000E+01  1.00800000E+00  1.20100000E+01  1.00800000E+00
  1.00800000E+00  1.00800000E+00  1.60000000E+01  1.00800000E+00  1.00800000E+00
  1.60000000E+01  1.00800000E+00  1.00800000E+00  1.60000000E+01  1.00800000E+00
  1.00800000E+00  2.29900000E+01  3.54500000E+01
%FLAG ATOM_TYPE_INDEX
%FORMAT(10I8)
       1       2       2       2       3       4       5       6       1       7
       7       7       8       9       9       8       9       9       8       9
       9      10      11
%FLAG NUMBER_EXCLUDED_ATOMS
%FORMAT(10I8)
       8       5       4       3       7       3       5       4       3       2
       1       1       2       1       1       2       1       1       2       1
       1       1       1
%FLAG NONBONDED_PARM_INDEX
%FORMAT(10I8)
       1       2       4       7      11      16      22      29      37      46
      56       2       3       5       8      12      17      23      30      38
      47      57       4       5       6       9      13      18      24      31
      39      48      58       7       8       9      10      14      19      25
      32      40      49      59      11      12      13      14      15      20
      26      33      41      50      60      16      17      18      19      20
      21      27      34      42      51      61      22      23      24      25
      26      27      28      35      43      52      62      29      30      31
      32      33      34      35      36      44      53      63      37      38
      39      40      41      42      43      44      45      54      64      46
      47      48      49      50      51      52      53      54      55      65
      56      57      58      59      60      61      62      63      64      65
      66
%FLAG RESIDUE_LABEL
%FORMAT(20a4)
ACE NME WAT WAT WAT Na+ Cl- 
%FLAG RESIDUE_POINTER
%FORMAT(10I8)
       1       7      13      16      19      22      23
%FLAG BOND_FORCE_CONSTANT
%FORMAT(5E16.8)
  3.40000000E+02  3.17000000E+02  5.70000000E+02  4.90000000E+02  4.34000000E+02
  3.37000000E+02  3.40000000E+02  5.53000000E+02  5.53000000E+02
%FLAG BOND_EQUIL_VALUE
%FORMAT(5E16.8)
  1.09000000E+00  1.52200000E+00  1.22900000E+00  1.33500000E+00  1.01000000E+00
  1.44900000E+00  1.09000000E+00  9.57200000E-01  1.51360000E+00
%FLAG ANGLE_FORCE_CONSTANT
%FORMAT(5E16.8)
  3.50000000E+01  5.00000000E+01  8.00000000E+01  7.00000000E+01  8.00000000E+01
  5.00000000E+01  5.00000000E+01  5.00000000E+01  5.00000000E+01  3.50000000E+01
%FLAG ANGLE_EQUIL_VALUE
%FORMAT(5E16.8)
  1.91113553E+00  1.91113553E+00  2.10137642E+00  2.03505391E+00  2.14500965E+00
  2.09439510E+00  2.12755636E+00  2.06018665E+00  1.91113553E+00  1.91113553E+00
%FLAG DIHEDRAL_FORCE_CONSTANT
%FORMAT(5E16.8)
  8.00000000E-01  8.00000000E-02  0.00000000E+00  2.50000000E+00  2.00000000E+00
  3.00000000E-01  1.50000000E-01  1.05000000E+01  1.00000000E+00
%FLAG DIHEDRAL_PERIODICITY
%FORMAT(5E16.8)
  1.00000000E+00  3.00000000E+00  2.00000000E+00  2.00000000E+00  1.00000000E+00
  4.00000000E+00  3.00000000E+00  2.00000000E+00  2.00000000E+00
%FLAG DIHEDRAL_PHASE
%FORMAT(5E16.8)
  0.00000000E+00  3.14159265E+00  0.00000000E+00  3.14159265E+00  0.00000000E+00
  3.14159265E+00  1.04719755E+00  3.14159265E+00  3.14159265E+00
%FLAG SCEE_SCALE_FACTOR
%FORMAT(5E16.8)
  1.20000000E+00  1.20000000E+00  1.20000000E+00  1.20000000E+00  1.20000000E+00
  1.20000000E+00  1.20000000E+00  1.20000000E+00  1.20000000E+00
%FLAG SCNB_SCALE_FACTOR
%FORMAT(5E16.8)
  2.00000000E+00  2.00000000E+00  2.00000000E+00  2.00000000E+00  2.00000000E+00
  2.00000000E+00  2.00000000E+00  2.00000000E+00  2.00000000E+00
%FLAG SOLTY
%FORMAT(5E16.8)
  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00
  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00
  0.00000000E+00
%FLAG LENNARD_JONES_ACOEF
%FORMAT(5E16.8)
  1.04308023E+06  9.71708117E+04  7.51607703E+03  9.24822270E+05  8.61541883E+04
  8.19971662E+05  6.47841731E+05  5.44261042E+04  5.74393458E+05  3.79876399E+05
  9.95480466E+05  8.96776989E+04  8.82619071E+05  6.06829342E+05  9.44293233E+05
  2.56678134E+03  1.07193646E+02  2.27577561E+03  1.02595236E+03  2.12601181E+03
  1.39982777E-01  6.78771368E+04  4.98586848E+03  6.01816484E+04  3.69471530E+04
  6.20665997E+04  5.94667300E+01  3.25969625E+03  7.85890042E+05  6.91773368E+04
  6.96790708E+05  4.72934643E+05  7.42364908E+05  1.52094041E+03  4.75728872E+04
  5.81935564E+05  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00
  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00
  1.49995867E+05  1.09119307E+04  1.32990267E+05  8.12116828E+04  1.36919347E+05
  1.25820601E+02  7.11465296E+03  1.04820861E+05  0.00000000E+00  1.55205818E+04
  3.47876374E+06  3.96588227E+05  3.08436311E+06  2.41912766E+06  3.44482863E+06
  1.95780471E+04  2.92681455E+05  2.78932950E+06  0.00000000E+00  8.49440877E+05
  9.24719470E+06
%FLAG LENNARD_JONES_BCOEF
%FORMAT(5E16.8)
  6.75612247E+02  1.26919150E+02  2.17257828E+01  5.99015525E+02  1.12529845E+02
  5.31102864E+02  6.26720080E+02  1.11805549E+02  5.55666448E+02  5.64885984E+02
  7.36907417E+02  1.36131731E+02  6.53361429E+02  6.77220874E+02  8.01323529E+02
  2.06278363E+01  2.59456373E+00  1.82891803E+01  1.53505284E+01  2.09604198E+01
  9.37598976E-02  1.06076943E+02  1.76949863E+01  9.40505980E+01  9.21192136E+01
  1.13252061E+02  1.93248820E+00  1.43076527E+01  6.36687196E+02  1.16264660E+02
  5.64503554E+02  5.81361517E+02  6.90894667E+02  1.72393904E+01  9.64152120E+01
  5.94825035E+02  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00
  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00
  2.42242668E+02  4.02144725E+01  2.14778698E+02  2.09807371E+02  2.58405233E+02
  4.31824675E+00  3.24719551E+01  2.19857570E+02  0.00000000E+00  7.36779156E+01
  9.31819600E+02  1.93646595E+02  8.26175679E+02  9.14638042E+02  1.03528777E+03
  4.30253727E+01  1.66355651E+02  9.05890802E+02  0.00000000E+00  4.96397934E+02
  1.14737423E+03
%FLAG BONDS_INC_HYDROGEN
%FORMAT(10I8)
       0       3       1       0       6       1       0       9       1      18
      21       5      24      27       7      24      30       7      24      33
       7      36      39       8      36      42       8      39      42       9
      45      48       8      45      51       8      48      51       9      54
      57       8      54      60       8      57      60       9
%FLAG BONDS_WITHOUT_HYDROGEN
%FORMAT(10I8)
       0      12       2      12      15       3      12      18       4      18
      24       6
%FLAG ANGLES_INC_HYDROGEN
%FORMAT(10I8)
       3       0       6       1       3       0       9       1       6       0
       9       1       3       0      12       2       6       0      12       2
       9       0      12       2      12      18      21       6      21      18
      24       8      18      24      27       9      18      24      30       9
      18      24      33       9      27      24      30      10      27      24
      33      10      30      24      33      10
%FLAG ANGLES_WITHOUT_HYDROGEN
%FORMAT(10I8)
       0      12      15       3       0      12      18       4      15      12
      18       5      12      18      24       7
%FLAG DIHEDRALS_INC_HYDROGEN
%FORMAT(10I8)
       3       0      12      15       1       3       0     -12      15       2
       3       0      12      18       3       6       0      12      15       1
       6       0     -12      15       2       6       0      12      18       3
       9       0      12      15       1       9       0     -12      15       2
       9       0      12      18       3       0      12      18      21       4
      15      12      18      21       4      15      12     -18      21       5
      12      18      24      27       6      21      18      24      27       7
      12      18      24      30       6      21      18      24      30       7
      12      18      24      33       6      21      18      24      33       7
       0      18     -12     -15       8      12      24     -18     -21       9
%FLAG DIHEDRALS_WITHOUT_HYDROGEN
%FORMAT(10I8)
       0      12      18      24       4      15      12      18      24       4
%FLAG EXCLUDED_ATOMS_LIST
%FORMAT(10I8)
       2       3       4       5       6       7       8       9       3       4
       5       6       7       4       5       6       7       5       6       7
       6       7       8       9      10      11      12       7       8       9
       8       9      10      11      12       9      10      11      12      10
      11      12      11      12      12       0      14      15      15       0
      17      18      18       0      20      21      21       0       0       0
%FLAG HBOND_ACOEF
%FORMAT(5E16.8)

%FLAG HBOND_BCOEF
%FORMAT(5E16.8)

%FLAG HBCUT
%FORMAT(5E16.8)

%FLAG AMBER_ATOM_TYPE
%FORMAT(20a4)
CT  HC  HC  HC  C   O   N   H   CT  H1  H1  H1  OW  HW  HW  OW  HW  HW  OW  HW  
HW  Na+ Cl- 
%FLAG TREE_CHAIN_CLASSIFICATION
%FORMAT(20a4)
BLA BLA BLA BLA BLA BLA BLA BLA BLA BLA BLA BLA BLA BLA BLA BLA BLA BLA BLA BLA 
BLA BLA BLA 
%FLAG JOIN_ARRAY
%FORMAT(10I8)
       0       0       0       0       0       0       0       0       0       0
       0       0       0       0       0       0       0       0       0       0
       0       0       0
%FLAG IROTAT
%FORMAT(10I8)
       0       0       0       0       0       0       0       0       0       0
       0       0       0       0       0       0       0       0       0       0
       0       0       0
%FLAG RADIUS_SET
%FORMAT(20a4)
modified Bondi radii (mbondi)
%FLAG RADII
%FORMAT(5E16.8)
  1.50000000E+00  1.50000000E+00  1.50000000E+00  1.50000000E+00  1.50000000E+00
  1.50000000E+00  1.50000000E+00  1.50000000E+00  1.50000000E+00  1.50000000E+00
  1.50000000E+00  1.50000000E+00  1.50000000E+00  1.50000000E+00  1.50000000E+00
  1.50000000E+00  1.50000000E+00  1.50000000E+00  1.50000000E+00  1.50000000E+00
  1.50000000E+00  1.50000000E+00  1.50000000E+00
%FLAG SCREEN
%FORMAT(5E16.8)
  8.00000000E-01  8.00000000E-01  8.00000000E-01  8.00000000E-01  8.00000000E-01
  8.00000000E-01  8.00000000E-01  8.00000000E-01  8.00000000E-01  8.00000000E-01
  8.00000000E-01  8.00000000E-01  8.00000000E-01  8.00000000E-01  8.00000000E-01
  8.00000000E-01  8.00000000E-01  8.00000000E-01  8.00000000E-01  8.00000000E-01
  8.00000000E-01  8.00000000E-01  8.00000000E-01
%FLAG IPOL
%FORMAT(10I8)
       0
)PRMTOP";

const char *inpcrdText = R"INPCRD(NMA + 3 TIP3P + NaCl (RASPA test)
    23
   0.0000000   0.0000000   0.0000000  -0.3728020   0.7846325  -0.6583848
  -0.3399187  -0.9949593  -0.3232819  -0.3870374   0.2609587   0.9739110
   1.5300000   0.0000000   0.0000000   2.1686472   1.0467398   0.1845685
   2.1383473  -1.1939487  -0.0000000   3.1547012  -1.2310874  -0.0777527
   1.4022048  -2.4314084   0.1711370   0.8014257  -2.6130800  -0.7200196
   2.0949877  -3.2602100   0.3168754   0.7556903  -2.3229495   1.0544578
   5.5000000   1.5000000   3.0000000   5.2981886   1.0891778   2.1593282
   4.7381494   1.3111838   3.5478722  -4.0000000   3.2000000  -2.5000000
  -4.0852900   3.2616200  -1.5486008  -3.7562293   4.0834285  -2.7763365
   1.5000000  -4.5000000  -3.2000000   1.9569810  -3.6677743  -3.3216580
   0.8771621  -4.3301625  -2.4932752  -4.0000000  -3.5000000   3.0000000
   5.0000000  -3.0000000  -1.5000000
)INPCRD";

class TemporaryDirectory
{
 public:
  TemporaryDirectory()
  {
    static std::size_t counter{};
    path = std::filesystem::temp_directory_path() /
           std::format("raspa3-prmtop-test-{}-{}", std::random_device{}(), counter++);
    std::filesystem::create_directories(path);
  }
  ~TemporaryDirectory()
  {
    std::error_code ignored;
    std::filesystem::remove_all(path, ignored);
  }
  std::filesystem::path path;
};

void writeText(const std::filesystem::path &file, const char *text)
{
  std::ofstream stream(file);
  stream << text;
}

struct TestFiles
{
  TemporaryDirectory directory{};
  std::filesystem::path prmtop{};
  std::filesystem::path inpcrd{};
  TestFiles(const char *coordinates = inpcrdText)
  {
    prmtop = directory.path / "nma.prmtop";
    inpcrd = directory.path / "nma.inpcrd";
    writeText(prmtop, prmtopText);
    writeText(inpcrd, coordinates);
  }
};

double3 toDouble3(const std::array<double, 3> &p) { return double3(p[0], p[1], p[2]); }

double dihedralAngle(const double3 &a, const double3 &b, const double3 &c, const double3 &d)
{
  // IUPAC convention (AMBER): trans = 180 degrees, sign by the right-hand rule about b -> c
  double3 b1 = b - a, b2 = c - b, b3 = d - c;
  double3 n1 = double3::cross(b1, b2);
  double3 n2 = double3::cross(b2, b3);
  double3 m = double3::cross(n1, b2.normalized());
  return std::atan2(double3::dot(m, n2), double3::dot(n1, n2));
}

double bendAngle(const double3 &a, const double3 &b, const double3 &c)
{
  double3 u = (a - b).normalized(), v = (c - b).normalized();
  return std::acos(std::clamp(double3::dot(u, v), -1.0, 1.0));
}

struct Reference
{
  double bond{}, angle{}, dihedral{};
  double intraLJ{}, intraCoulomb{};  // within a molecule: scaled 1-4 plus every non-excluded pair
  double interLJ{}, interCoulomb{};
};

// Energies in RASPA's internal units from the AMBER functional forms.
Reference referenceEnergies(const AMBER::Prmtop &prmtop, const std::vector<std::array<double, 3>> &positions)
{
  const double kcal = Units::KCalPerMolToEnergy;
  const std::size_t n = prmtop.numberOfAtoms();
  const std::size_t ntypes = prmtop.numberOfTypes();
  Reference reference{};

  // molecule of every atom (connected components of the bond graph)
  std::vector<std::size_t> molecule(n);
  std::iota(molecule.begin(), molecule.end(), std::size_t{0});
  auto root = [&](std::size_t a)
  {
    while (molecule[a] != a) a = molecule[a];
    return a;
  };
  for (const std::vector<std::int64_t> *list : {&prmtop.bondsIncHydrogen, &prmtop.bondsWithoutHydrogen})
    for (std::size_t i = 0; i + 2 < list->size(); i += 3)
    {
      std::size_t a = root(static_cast<std::size_t>((*list)[i] / 3)), b = root(static_cast<std::size_t>((*list)[i + 1] / 3));
      if (a != b) molecule[std::max(a, b)] = std::min(a, b);
    }

  for (const std::vector<std::int64_t> *list : {&prmtop.bondsIncHydrogen, &prmtop.bondsWithoutHydrogen})
    for (std::size_t i = 0; i + 2 < list->size(); i += 3)
    {
      std::size_t a = static_cast<std::size_t>((*list)[i] / 3), b = static_cast<std::size_t>((*list)[i + 1] / 3);
      std::size_t type = static_cast<std::size_t>((*list)[i + 2] - 1);
      double r = (toDouble3(positions[a]) - toDouble3(positions[b])).length();
      reference.bond += prmtop.bondForceConstant[type] * std::pow(r - prmtop.bondEquilValue[type], 2) * kcal;
    }
  for (const std::vector<std::int64_t> *list : {&prmtop.anglesIncHydrogen, &prmtop.anglesWithoutHydrogen})
    for (std::size_t i = 0; i + 3 < list->size(); i += 4)
    {
      std::size_t a = static_cast<std::size_t>((*list)[i] / 3), b = static_cast<std::size_t>((*list)[i + 1] / 3),
                  c = static_cast<std::size_t>((*list)[i + 2] / 3);
      std::size_t type = static_cast<std::size_t>((*list)[i + 3] - 1);
      double theta = bendAngle(toDouble3(positions[a]), toDouble3(positions[b]), toDouble3(positions[c]));
      reference.angle += prmtop.angleForceConstant[type] * std::pow(theta - prmtop.angleEquilValue[type], 2) * kcal;
    }

  auto lennardJones = [&](std::size_t a, std::size_t b, double r)
  {
    std::size_t ta = static_cast<std::size_t>(prmtop.atomTypeIndex[a] - 1), tb = static_cast<std::size_t>(prmtop.atomTypeIndex[b] - 1);
    std::int64_t index = prmtop.nonbondedParmIndex[ta * ntypes + tb];
    if (index <= 0) return 0.0;
    double A = prmtop.lennardJonesACoef[static_cast<std::size_t>(index - 1)];
    double B = prmtop.lennardJonesBCoef[static_cast<std::size_t>(index - 1)];
    return (A / std::pow(r, 12) - B / std::pow(r, 6)) * kcal;
  };
  auto coulomb = [&](std::size_t a, std::size_t b, double r)
  { return Units::CoulombicConversionFactor * prmtop.charges[a] * prmtop.charges[b] / r; };

  for (const std::vector<std::int64_t> *list : {&prmtop.dihedralsIncHydrogen, &prmtop.dihedralsWithoutHydrogen})
    for (std::size_t i = 0; i + 4 < list->size(); i += 5)
    {
      std::size_t a = static_cast<std::size_t>(std::abs((*list)[i]) / 3), b = static_cast<std::size_t>(std::abs((*list)[i + 1]) / 3),
                  c = static_cast<std::size_t>(std::abs((*list)[i + 2]) / 3), d = static_cast<std::size_t>(std::abs((*list)[i + 3]) / 3);
      std::size_t type = static_cast<std::size_t>((*list)[i + 4] - 1);
      double phi = dihedralAngle(toDouble3(positions[a]), toDouble3(positions[b]), toDouble3(positions[c]), toDouble3(positions[d]));
      reference.dihedral += prmtop.dihedralForceConstant[type] *
                            (1.0 + std::cos(prmtop.dihedralPeriodicity[type] * phi - prmtop.dihedralPhase[type])) * kcal;
      bool computes14 = (*list)[i + 2] >= 0 && (*list)[i + 3] >= 0;
      if (computes14)
      {
        double r = (toDouble3(positions[a]) - toDouble3(positions[d])).length();
        reference.intraLJ += lennardJones(a, d, r) / prmtop.scnbScaleFactor[type];
        reference.intraCoulomb += coulomb(a, d, r) / prmtop.sceeScaleFactor[type];
      }
    }

  // non-bonded pairs outside the exclusion list
  std::vector<std::set<std::size_t>> excluded(n);
  {
    std::size_t offset = 0;
    for (std::size_t a = 0; a < n; ++a)
    {
      std::size_t count = static_cast<std::size_t>(prmtop.numberExcludedAtoms[a]);
      for (std::size_t k = 0; k < count; ++k)
      {
        std::int64_t b = prmtop.excludedAtomsList[offset + k];
        if (b > 0) excluded[a].insert(static_cast<std::size_t>(b - 1));
      }
      offset += count;
    }
  }
  for (std::size_t a = 0; a < n; ++a)
    for (std::size_t b = a + 1; b < n; ++b)
    {
      if (excluded[a].contains(b)) continue;
      double r = (toDouble3(positions[a]) - toDouble3(positions[b])).length();
      if (root(a) == root(b))
      {
        reference.intraLJ += lennardJones(a, b, r);
        reference.intraCoulomb += coulomb(a, b, r);
      }
      else
      {
        reference.interLJ += lennardJones(a, b, r);
        reference.interCoulomb += coulomb(a, b, r);
      }
    }
  return reference;
}

Component readComponent(const ForceField &forceField, const std::filesystem::path &directory, const std::string &name,
                        std::size_t componentId)
{
  std::filesystem::path stem = directory / name;
  return Component(Component::Type::Adsorbate, componentId, forceField, name, stem.string(), 5, 21);
}
}  // namespace

TEST(amber_prmtop, parses_sections)
{
  TestFiles files;
  AMBER::Prmtop prmtop = AMBER::parsePrmtop(files.prmtop);

  EXPECT_EQ(prmtop.title, "NMA + 3 TIP3P + NaCl (RASPA test)");
  EXPECT_EQ(prmtop.numberOfAtoms(), 23uz);
  EXPECT_EQ(prmtop.numberOfTypes(), 11uz);
  EXPECT_FALSE(prmtop.hasBox());
  ASSERT_EQ(prmtop.charges.size(), 23uz);
  EXPECT_NEAR(prmtop.charges[0], -0.3662, 1e-6);
  EXPECT_NEAR(prmtop.charges[21], 1.0, 1e-6);
  EXPECT_EQ(prmtop.atomNames[1], "HH31");
  EXPECT_EQ(prmtop.amberAtomTypes[21], "Na+");
  EXPECT_EQ(prmtop.amberAtomTypes[22], "Cl-");
  ASSERT_EQ(prmtop.residueLabels.size(), 7uz);
  EXPECT_EQ(prmtop.residueLabels[2], "WAT");
  EXPECT_EQ(prmtop.bondsIncHydrogen.size(), 3uz * 16);
  EXPECT_EQ(prmtop.bondsWithoutHydrogen.size(), 3uz * 4);
  EXPECT_EQ(prmtop.dihedralsIncHydrogen.size(), 5uz * 20);
  EXPECT_EQ(prmtop.excludedAtomsList.size(), 60uz);
  EXPECT_NEAR(prmtop.angleEquilValue[0], 109.5 * std::numbers::pi / 180.0, 1e-7);
  EXPECT_NEAR(prmtop.dihedralPhase[6], 60.0 * std::numbers::pi / 180.0, 1e-7);

  AMBER::Coordinates coordinates = AMBER::readCoordinates(files.inpcrd, 23);
  ASSERT_EQ(coordinates.positions.size(), 23uz);
  EXPECT_TRUE(coordinates.velocities.empty());
  EXPECT_FALSE(coordinates.box.has_value());
  EXPECT_NEAR(coordinates.positions[4][0], 1.53, 1e-9);
}

TEST(amber_prmtop, chamber_charges_use_the_charmm_coulomb_constant)
{
  // a CHAMBER prmtop (CTITLE) stores q * sqrt(332.0716) instead of q * 18.2223
  TestFiles files;
  std::string text = prmtopText;
  text.replace(text.find("%FLAG TITLE"), std::string("%FLAG TITLE").size(), "%FLAG CTITLE");
  writeText(files.prmtop, text.c_str());
  AMBER::Prmtop prmtop = AMBER::parsePrmtop(files.prmtop);
  EXPECT_EQ(prmtop.title, "NMA + 3 TIP3P + NaCl (RASPA test)");
  EXPECT_NEAR(prmtop.charges[21], 18.2223 / std::sqrt(332.0716), 1e-12);
  EXPECT_NEAR(prmtop.charges[0], -6.67300626 / std::sqrt(332.0716), 1e-12);
}

TEST(amber_prmtop, conversion_layout)
{
  TestFiles files;
  AMBER::ReadResult result = AMBER::readPrmtop(files.prmtop, files.inpcrd);

  ASSERT_EQ(result.components.size(), 4uz);
  const AMBER::ReadComponent &protein = result.components[0];
  const AMBER::ReadComponent &water = result.components[1];
  EXPECT_EQ(protein.name, "protein");
  EXPECT_EQ(protein.atomsPerMolecule, 12uz);
  EXPECT_EQ(protein.count, 1uz);
  EXPECT_FALSE(protein.rigid);
  EXPECT_EQ(water.name, "WAT");
  EXPECT_EQ(water.atomsPerMolecule, 3uz);
  EXPECT_EQ(water.count, 3uz);
  EXPECT_TRUE(water.rigid);
  EXPECT_FALSE(water.definition.contains("Connectivity"));
  EXPECT_EQ(result.components[2].name, "Na_plus");
  EXPECT_EQ(result.components[3].name, "Cl-");
  EXPECT_TRUE(result.components[2].rigid);

  // rigid water carries the equilibrium geometry: r(OH) = 0.9572, r(HH) = 1.5136
  {
    const nlohmann::json &atoms = water.definition["PseudoAtoms"];
    double3 o(atoms[0][1][0], atoms[0][1][1], atoms[0][1][2]);
    double3 h1(atoms[1][1][0], atoms[1][1][1], atoms[1][1][2]);
    double3 h2(atoms[2][1][0], atoms[2][1][1], atoms[2][1][2]);
    EXPECT_NEAR((o - h1).length(), 0.9572, 1e-9);
    EXPECT_NEAR((o - h2).length(), 0.9572, 1e-9);
    EXPECT_NEAR((h1 - h2).length(), 1.5136, 1e-9);
    EXPECT_EQ(atoms[0][0], "OW");
  }

  // bonded terms of the protein
  const nlohmann::json &definition = protein.definition;
  EXPECT_EQ(definition["Bonds"].size(), 11uz);
  EXPECT_EQ(definition["Bends"].size(), 18uz);
  EXPECT_EQ(definition["Torsions"].size(), 16uz);
  EXPECT_EQ(definition["ImproperTorsions"].size(), 5uz);  // 3 off-phase terms + 2 impropers
  EXPECT_NEAR(definition["Intra14VanDerWaalsScalingValue"].get<double>(), 0.5, 1e-12);
  EXPECT_NEAR(definition["Intra14ChargeChargeScalingValue"].get<double>(), 1.0 / 1.2, 1e-12);

  const double K = Units::KCalPerMolToEnergy * Units::EnergyToKelvin;
  // CT-HC bond: k = 340 kcal/mol/A^2 -> p_0 = 2k
  EXPECT_EQ(definition["Bonds"][0][1], "HARMONIC");
  EXPECT_NEAR(definition["Bonds"][0][2][0].get<double>(), 2.0 * 340.0 * K, 1e-6);
  EXPECT_NEAR(definition["Bonds"][0][2][1].get<double>(), 1.090, 1e-12);
  // HC-CT-C-O: 0.8 (n=1, 0) + 0.08 (n=3, 180) -> polynomial in cos(phi)
  bool foundTorsion = false;
  for (const nlohmann::json &torsion : definition["Torsions"])
  {
    if (torsion[0] != nlohmann::json({1, 0, 4, 5})) continue;
    foundTorsion = true;
    EXPECT_EQ(torsion[1], "POLYNOMIAL");
    EXPECT_NEAR(torsion[2][0].get<double>(), (0.8 + 0.08) * K, 1e-6);
    EXPECT_NEAR(torsion[2][1].get<double>(), (0.8 + 3.0 * 0.08) * K, 1e-6);
    EXPECT_NEAR(torsion[2][2].get<double>(), 0.0, 1e-9);
    EXPECT_NEAR(torsion[2][3].get<double>(), -4.0 * 0.08 * K, 1e-6);
  }
  EXPECT_TRUE(foundTorsion);

  // force field: one pseudo-atom per AMBER type, the Na+/Cl- pair deviates from Lorentz-Berthelot
  EXPECT_EQ(result.forceField["PseudoAtoms"].size(), 11uz);
  ASSERT_TRUE(result.forceField.contains("BinaryInteractions"));
  ASSERT_EQ(result.forceField["BinaryInteractions"].size(), 1uz);
  EXPECT_NEAR(result.forceField["BinaryInteractions"][0]["parameters"][0].get<double>(),
              1.3 * std::sqrt(0.0874393 * 0.0355910) * K, 1e-6);
  EXPECT_EQ(result.forceField["TruncationMethod"], "truncated");
  EXPECT_FALSE(result.forceField["TailCorrections"].get<bool>());

  // vacuum: direct Coulomb, box around the molecules, coordinates seeded through restart.json
  EXPECT_EQ(result.simulation["Systems"][0]["ChargeMethod"], "Coulomb");
  EXPECT_EQ(result.simulation["Systems"][0]["RestartFileName"], "restart.json");
  EXPECT_EQ(result.simulation["Components"][1]["CreateNumberOfMolecules"].get<std::size_t>(), 0uz);
  ASSERT_TRUE(result.restart.contains("SimulationBox"));
  EXPECT_EQ(result.restart["WAT"].size(), 9uz);
  EXPECT_EQ(result.restart["protein"].size(), 12uz);
  EXPECT_GT(result.boxLengths[0], 2.0 * result.forceField["CutOffVDW"].get<double>());
}

TEST(amber_prmtop, single_point_energy_matches_amber_functional_forms)
{
  TestFiles files;
  AMBER::ReadResult result = AMBER::readPrmtop(files.prmtop, files.inpcrd);
  TemporaryDirectory output;
  AMBER::writeRaspaInput(result, output.path);

  ForceField forceField((output.path / "force_field.json").string());
  forceField.chargeMethod = ForceField::ChargeMethod::Coulomb;  // as simulation.json's "ChargeMethod"
  forceField.useCharge = true;

  std::vector<Component> components{};
  std::vector<std::vector<double3>> initialPositions{};
  for (std::size_t i = 0; i < result.components.size(); ++i)
  {
    components.push_back(readComponent(forceField, output.path, result.components[i].name, i));
    initialPositions.push_back(result.restart[result.components[i].name].get<std::vector<double3>>());
  }
  EXPECT_FALSE(components[0].rigid);
  EXPECT_TRUE(components[1].rigid);
  EXPECT_EQ(components[0].intraMolecularPotentials.improperTorsions.size(), 5uz);

  SimulationBox box(result.boxLengths[0], result.boxLengths[1], result.boxLengths[2]);
  System system(forceField, box, false, 300.0, 1e5, 1.0, {}, components, initialPositions, {0, 0, 0, 0}, 5);
  ASSERT_EQ(system.spanOfMoleculeAtoms().size(), 23uz);

  // the rigid waters were fitted onto the AMBER coordinates
  {
    std::span<const Atom> atoms = system.spanOfMoleculeAtoms();
    for (std::size_t i = 12; i < 21; ++i)
    {
      EXPECT_NEAR(atoms[i].position.x, result.positions[i][0], 1e-6);
      EXPECT_NEAR(atoms[i].position.y, result.positions[i][1], 1e-6);
      EXPECT_NEAR(atoms[i].position.z, result.positions[i][2], 1e-6);
    }
  }

  RunningEnergy energy = system.computeTotalEnergies();

  AMBER::Prmtop prmtop = AMBER::parsePrmtop(files.prmtop);
  Reference reference = referenceEnergies(prmtop, result.positions);

  auto tolerance = [](double value) { return 1e-6 * std::max(1.0, std::abs(value)); };
  EXPECT_NEAR(energy.bond, reference.bond, tolerance(reference.bond));
  EXPECT_NEAR(energy.bend, reference.angle, tolerance(reference.angle));
  EXPECT_NEAR(energy.torsion + energy.improperTorsion, reference.dihedral, tolerance(reference.dihedral));
  EXPECT_NEAR(energy.intraVDW, reference.intraLJ, tolerance(reference.intraLJ));
  EXPECT_NEAR(energy.intraCoul, reference.intraCoulomb, tolerance(reference.intraCoulomb));
  EXPECT_NEAR(energy.moleculeMoleculeVDW, reference.interLJ, tolerance(reference.interLJ));
  EXPECT_NEAR(energy.moleculeMoleculeCharge, reference.interCoulomb, tolerance(reference.interCoulomb));

  // the terms are not trivially zero
  EXPECT_GT(reference.bond * Units::EnergyToKelvin, 1.0);
  EXPECT_GT(reference.angle * Units::EnergyToKelvin, 1.0);
  EXPECT_GT(std::abs(reference.dihedral) * Units::EnergyToKelvin, 1.0);
  EXPECT_GT(std::abs(reference.interCoulomb) * Units::EnergyToKelvin, 1.0);
}

TEST(amber_prmtop, periodic_coordinates_select_ewald)
{
  std::string withBox = std::string(inpcrdText) +
                        "  30.0000000  30.0000000  30.0000000  90.0000000  90.0000000  90.0000000\n";
  TestFiles files(withBox.c_str());
  AMBER::ReadResult result = AMBER::readPrmtop(files.prmtop, files.inpcrd);

  EXPECT_EQ(result.simulation["Systems"][0]["ChargeMethod"], "Ewald");
  EXPECT_NEAR(result.boxLengths[0], 30.0, 1e-12);
  EXPECT_NEAR(result.restart["SimulationBox"]["length-c"].get<double>(), 30.0, 1e-12);
  EXPECT_NEAR(result.forceField["CutOffVDW"].get<double>(), 10.0, 1e-12);
  EXPECT_TRUE(result.forceField.contains("EwaldPrecision"));
  // coordinates are not shifted in a periodic system
  EXPECT_NEAR(result.positions[4][0], 1.53, 1e-9);

  TemporaryDirectory output;
  AMBER::writeRaspaInput(result, output.path);
  EXPECT_TRUE(std::filesystem::exists(output.path / "restart.json"));
  EXPECT_TRUE(std::filesystem::exists(output.path / "WAT.json"));
  EXPECT_TRUE(std::filesystem::exists(output.path / "conversion_warnings.txt"));
}
