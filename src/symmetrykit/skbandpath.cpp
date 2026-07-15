module;

module skbandpath;

import std;

import double3;
import double3x3;
import int3x3;

import skdefinitions;
import sksymmetrycell;
import sktransformationmatrix;
import skspacegroup;
import skpointgroup;

namespace
{

// ---------------------------------------------------------------------------------------------------------------------
// Small geometry helpers (lattice vectors are stored as columns of a double3x3).
// ---------------------------------------------------------------------------------------------------------------------

struct CellParameters
{
  double a, b, c;
  double cosAlpha, cosBeta, cosGamma;
};

CellParameters cellParameters(const double3x3& cell)
{
  double3 v1 = cell[0];
  double3 v2 = cell[1];
  double3 v3 = cell[2];
  double a = v1.length();
  double b = v2.length();
  double c = v3.length();
  return CellParameters{a,
                        b,
                        c,
                        double3::dot(v2, v3) / (b * c),
                        double3::dot(v1, v3) / (a * c),
                        double3::dot(v1, v2) / (a * b)};
}

// Build a cell (columns are lattice vectors) from cell parameters (angles in radians).
double3x3 cellFromParameters(double a, double b, double c, double alpha, double beta, double gamma)
{
  double ca = std::cos(alpha);
  double cb = std::cos(beta);
  double cg = std::cos(gamma);
  double sg = std::sin(gamma);
  double volumeFactor = std::sqrt(std::max(0.0, 1.0 - ca * ca - cb * cb - cg * cg + 2.0 * ca * cb * cg));
  double3 va(a, 0.0, 0.0);
  double3 vb(b * cg, b * sg, 0.0);
  double3 vc(c * cb, c * (ca - cb * cg) / sg, c * volumeFactor / sg);
  return double3x3(va, vb, vc);
}

// Reciprocal lattice (columns are reciprocal vectors), without the 2*pi factor: (cell^{-1})^T.
double3x3 reciprocalLattice(const double3x3& cell) { return cell.inverse().transpose(); }

// Cosines of the angles between the columns of a matrix.
struct AngleCosines
{
  double alpha, beta, gamma;
};

AngleCosines columnAngleCosines(const double3x3& m)
{
  double3 v1 = m[0];
  double3 v2 = m[1];
  double3 v3 = m[2];
  return AngleCosines{double3::dot(v2, v3) / (v2.length() * v3.length()),
                      double3::dot(v1, v3) / (v1.length() * v3.length()),
                      double3::dot(v1, v2) / (v1.length() * v2.length())};
}

// ---------------------------------------------------------------------------------------------------------------------
// Bravais-lattice classification from the space-group number.
// ---------------------------------------------------------------------------------------------------------------------

// Crystal-family letter: 'a' triclinic, 'm' monoclinic, 'o' orthorhombic, 't' tetragonal,
// 'h' trigonal/hexagonal, 'c' cubic.
char crystalFamilyLetter(std::size_t number)
{
  if (number <= 2) return 'a';
  if (number <= 15) return 'm';
  if (number <= 74) return 'o';
  if (number <= 142) return 't';
  if (number <= 194) return 'h';
  return 'c';
}

// Conventional-setting centring letter (P, C, A, F, I, R) for each space group (1..230).
char centringLetter(std::size_t number)
{
  auto in = [number](std::initializer_list<std::size_t> list)
  { return std::find(list.begin(), list.end(), number) != list.end(); };

  // monoclinic
  if (in({5, 8, 9, 12, 15})) return 'C';
  // orthorhombic
  if (in({20, 21, 35, 36, 37, 63, 64, 65, 66, 67, 68})) return 'C';
  if (in({38, 39, 40, 41})) return 'A';
  if (in({22, 42, 43, 69, 70})) return 'F';
  if (in({23, 24, 44, 45, 46, 71, 72, 73, 74})) return 'I';
  // tetragonal
  if (in({79, 80, 82, 87, 88, 97, 98, 107, 108, 109, 110, 119, 120, 121, 122, 139, 140, 141, 142})) return 'I';
  // trigonal (rhombohedral)
  if (in({146, 148, 155, 160, 161, 166, 167})) return 'R';
  // cubic
  if (in({196, 202, 203, 209, 210, 216, 219, 225, 226, 227, 228})) return 'F';
  if (in({197, 199, 204, 206, 211, 214, 217, 220, 229, 230})) return 'I';
  return 'P';
}

// ---------------------------------------------------------------------------------------------------------------------
// Change-of-basis matrix P converting conventional to primitive cell: (a_p b_p c_p) = (a b c) P.
// (HPKOT Table 3; the mC and oA conventions differ from the default spglib primitive cell.)
// ---------------------------------------------------------------------------------------------------------------------
double3x3 primitiveTransformation(const std::string& bravais)
{
  if (bravais == "cP" || bravais == "tP" || bravais == "hP" || bravais == "oP" || bravais == "mP" || bravais == "aP")
  {
    return double3x3::identity();
  }
  if (bravais == "cF" || bravais == "oF")
  {
    return double3x3(double3(0.0, 1.0, 1.0), double3(1.0, 0.0, 1.0), double3(1.0, 1.0, 0.0)) * 0.5;
  }
  if (bravais == "cI" || bravais == "tI" || bravais == "oI")
  {
    return double3x3(double3(-1.0, 1.0, 1.0), double3(1.0, -1.0, 1.0), double3(1.0, 1.0, -1.0)) * 0.5;
  }
  if (bravais == "hR")
  {
    return double3x3(double3(2.0, 1.0, 1.0), double3(-1.0, 1.0, 1.0), double3(-1.0, -2.0, 1.0)) * (1.0 / 3.0);
  }
  if (bravais == "oC")
  {
    return double3x3(double3(1.0, -1.0, 0.0), double3(1.0, 1.0, 0.0), double3(0.0, 0.0, 2.0)) * 0.5;
  }
  if (bravais == "oA")
  {
    return double3x3(double3(0.0, 1.0, -1.0), double3(0.0, 1.0, 1.0), double3(2.0, 0.0, 0.0)) * 0.5;
  }
  if (bravais == "mC")
  {
    return double3x3(double3(1.0, 1.0, 0.0), double3(-1.0, 1.0, 0.0), double3(0.0, 0.0, 2.0)) * 0.5;
  }
  return double3x3::identity();
}

// ---------------------------------------------------------------------------------------------------------------------
// Triclinic (aP) standardization: Niggli-reduce the reciprocal cell, then bring the reciprocal angles to an
// all-acute (aP3) or all-obtuse (aP2) form (HPKOT section on the aP lattice).
// Returns (extended-Bravais symbol, standardized real conventional cell in columns).
// ---------------------------------------------------------------------------------------------------------------------
std::pair<std::string, double3x3> standardizeTriclinic(const double3x3& conventionalCell)
{
  // Niggli-reduce the reciprocal lattice.
  double3x3 reciprocal0 = reciprocalLattice(conventionalCell);
  SKSymmetryCell reciprocalCell = SKSymmetryCell::createFromUnitCell(reciprocal0);
  auto reduced = reciprocalCell.computeReducedNiggliCellAndChangeOfBasisMatrix();

  double3x3 reciprocal2 = reciprocal0;
  if (reduced)
  {
    // Reduced (oriented) reciprocal vectors: B_reduced = B_0 * M, with M the integer change-of-basis matrix.
    reciprocal2 = reciprocal0 * double3x3(reduced->second.transformation);
  }
  double3x3 realCell2 = reciprocalLattice(reciprocal2);

  // Choose the axis permutation M2 that makes the |k_i k_j cos(k_angle)| product smallest.
  AngleCosines recip2 = columnAngleCosines(reciprocal2);
  double ka = reciprocal2[0].length();
  double kb = reciprocal2[1].length();
  double kc = reciprocal2[2].length();
  std::array<double, 3> conditions = {std::fabs(kb * kc * recip2.alpha), std::fabs(kc * ka * recip2.beta),
                                       std::fabs(ka * kb * recip2.gamma)};
  std::size_t smallest =
      static_cast<std::size_t>(std::distance(conditions.begin(), std::min_element(conditions.begin(), conditions.end())));

  double3x3 M2 = double3x3::identity();
  if (smallest == 0)
  {
    M2 = double3x3(double3(0.0, 1.0, 0.0), double3(0.0, 0.0, 1.0), double3(1.0, 0.0, 0.0));
  }
  else if (smallest == 1)
  {
    M2 = double3x3(double3(0.0, 0.0, 1.0), double3(1.0, 0.0, 0.0), double3(0.0, 1.0, 0.0));
  }
  double3x3 realCell3 = realCell2 * M2;

  // Sign-flip matrix M3 to make the reciprocal angles all-acute or all-obtuse.
  AngleCosines recip3 = columnAngleCosines(reciprocalLattice(realCell3));
  bool a = recip3.alpha > 0.0;
  bool b = recip3.beta > 0.0;
  bool g = recip3.gamma > 0.0;
  double3x3 M3 = double3x3::identity();
  if ((a && b && g) || (!a && !b && !g))
  {
    M3 = double3x3(double3(1.0, 0.0, 0.0), double3(0.0, 1.0, 0.0), double3(0.0, 0.0, 1.0));
  }
  else if ((a && !b && !g) || (!a && b && g))
  {
    M3 = double3x3(double3(1.0, 0.0, 0.0), double3(0.0, -1.0, 0.0), double3(0.0, 0.0, -1.0));
  }
  else if ((!a && b && !g) || (a && !b && g))
  {
    M3 = double3x3(double3(-1.0, 0.0, 0.0), double3(0.0, 1.0, 0.0), double3(0.0, 0.0, -1.0));
  }
  else
  {
    M3 = double3x3(double3(-1.0, 0.0, 0.0), double3(0.0, -1.0, 0.0), double3(0.0, 0.0, 1.0));
  }
  double3x3 realCellFinal = realCell3 * M3;

  AngleCosines recipFinal = columnAngleCosines(reciprocalLattice(realCellFinal));
  std::string ext = (recipFinal.alpha <= 0.0 && recipFinal.beta <= 0.0 && recipFinal.gamma <= 0.0) ? "aP2" : "aP3";
  return {ext, realCellFinal};
}

// ---------------------------------------------------------------------------------------------------------------------
// Determine the extended (case-specific) Bravais-lattice symbol.
// ---------------------------------------------------------------------------------------------------------------------
std::string extendedBravaisSymbol(const std::string& bravais, std::size_t number, const CellParameters& p,
                                  double threshold)
{
  (void)threshold;
  double a = p.a, b = p.b, c = p.c;
  double cosBeta = p.cosBeta;
  double sin2Beta = 1.0 - cosBeta * cosBeta;

  if (bravais == "cP") return (number <= 206) ? "cP1" : "cP2";
  if (bravais == "cF") return (number <= 206) ? "cF1" : "cF2";
  if (bravais == "cI") return "cI1";
  if (bravais == "tP") return "tP1";
  if (bravais == "tI") return (c <= a) ? "tI1" : "tI2";
  if (bravais == "oP") return "oP1";
  if (bravais == "oF")
  {
    double ia = 1.0 / (a * a);
    double ib = 1.0 / (b * b);
    double ic = 1.0 / (c * c);
    if (ia > ib + ic) return "oF1";
    if (ic > ia + ib) return "oF2";
    return "oF3";
  }
  if (bravais == "oI")
  {
    // Index of the longest of the three axes: c -> 1, a -> 2, b -> 3.
    if (c >= a && c >= b) return "oI1";
    if (a >= b && a >= c) return "oI2";
    return "oI3";
  }
  if (bravais == "oC") return (a <= b) ? "oC1" : "oC2";
  if (bravais == "oA") return (b <= c) ? "oA1" : "oA2";
  if (bravais == "hP")
  {
    static const std::array<std::size_t, 15> hp1 = {143, 144, 145, 146, 147, 148, 149, 151, 153, 157, 159, 160, 161, 162, 163};
    return (std::find(hp1.begin(), hp1.end(), number) != hp1.end()) ? "hP1" : "hP2";
  }
  if (bravais == "hR") return (std::sqrt(3.0) * a <= std::sqrt(2.0) * c) ? "hR1" : "hR2";
  if (bravais == "mP") return "mP1";
  if (bravais == "mC")
  {
    if (b < a * std::sqrt(sin2Beta)) return "mC1";
    double test = -a * cosBeta / c + a * a * sin2Beta / (b * b);
    return (test <= 1.0) ? "mC2" : "mC3";
  }
  return "";  // aP handled separately
}

// ---------------------------------------------------------------------------------------------------------------------
// Article-derived tables of labelled k-points and recommended paths.
// ---------------------------------------------------------------------------------------------------------------------
struct Tables
{
  std::vector<SKKVector> points;
  std::vector<std::pair<std::string, std::string>> segments;

  void p(std::string label, double x, double y, double z) { points.push_back({std::move(label), double3(x, y, z)}); }
  void s(std::string from, std::string to) { segments.emplace_back(std::move(from), std::move(to)); }
};

Tables buildTables(const std::string& ext, const CellParameters& cp)
{
  const double a = cp.a, b = cp.b, c = cp.c;
  const double cosBeta = cp.cosBeta;
  const double sin2Beta = 1.0 - cosBeta * cosBeta;

  Tables t;

  if (ext == "cP1" || ext == "cP2")
  {
    t.p("GAMMA", 0, 0, 0);
    t.p("R", 0.5, 0.5, 0.5);
    t.p("M", 0.5, 0.5, 0);
    t.p("X", 0, 0.5, 0);
    t.p("X_1", 0.5, 0, 0);
    t.s("GAMMA", "X"); t.s("X", "M"); t.s("M", "GAMMA"); t.s("GAMMA", "R"); t.s("R", "X"); t.s("R", "M");
    if (ext == "cP1") t.s("M", "X_1");
  }
  else if (ext == "cF1" || ext == "cF2")
  {
    t.p("GAMMA", 0, 0, 0);
    t.p("X", 0.5, 0, 0.5);
    t.p("L", 0.5, 0.5, 0.5);
    t.p("W", 0.5, 0.25, 0.75);
    t.p("W_2", 0.75, 0.25, 0.5);
    t.p("K", 0.375, 0.375, 0.75);
    t.p("U", 0.625, 0.25, 0.625);
    t.s("GAMMA", "X"); t.s("X", "U"); t.s("K", "GAMMA"); t.s("GAMMA", "L"); t.s("L", "W"); t.s("W", "X");
    if (ext == "cF1") t.s("X", "W_2");
  }
  else if (ext == "cI1")
  {
    t.p("GAMMA", 0, 0, 0);
    t.p("H", 0.5, -0.5, 0.5);
    t.p("P", 0.25, 0.25, 0.25);
    t.p("N", 0, 0, 0.5);
    t.s("GAMMA", "H"); t.s("H", "N"); t.s("N", "GAMMA"); t.s("GAMMA", "P"); t.s("P", "H"); t.s("P", "N");
  }
  else if (ext == "tP1")
  {
    t.p("GAMMA", 0, 0, 0);
    t.p("Z", 0, 0, 0.5);
    t.p("M", 0.5, 0.5, 0);
    t.p("A", 0.5, 0.5, 0.5);
    t.p("R", 0, 0.5, 0.5);
    t.p("X", 0, 0.5, 0);
    t.s("GAMMA", "X"); t.s("X", "M"); t.s("M", "GAMMA"); t.s("GAMMA", "Z"); t.s("Z", "R"); t.s("R", "A");
    t.s("A", "Z"); t.s("X", "R"); t.s("M", "A");
  }
  else if (ext == "tI1")
  {
    double H = (1.0 + c * c / (a * a)) / 4.0;
    t.p("GAMMA", 0, 0, 0);
    t.p("M", -0.5, 0.5, 0.5);
    t.p("X", 0, 0, 0.5);
    t.p("P", 0.25, 0.25, 0.25);
    t.p("Z", H, H, -H);
    t.p("Z_0", -H, 1.0 - H, H);
    t.p("N", 0, 0.5, 0);
    t.s("GAMMA", "X"); t.s("X", "M"); t.s("M", "GAMMA"); t.s("GAMMA", "Z"); t.s("Z_0", "M"); t.s("X", "P");
    t.s("P", "N"); t.s("N", "GAMMA");
  }
  else if (ext == "tI2")
  {
    double H = (1.0 + a * a / (c * c)) / 4.0;
    double Z = a * a / (2.0 * c * c);
    t.p("GAMMA", 0, 0, 0);
    t.p("M", 0.5, 0.5, -0.5);
    t.p("X", 0, 0, 0.5);
    t.p("P", 0.25, 0.25, 0.25);
    t.p("N", 0, 0.5, 0);
    t.p("S_0", -H, H, H);
    t.p("S", H, 1.0 - H, -H);
    t.p("R", -Z, Z, 0.5);
    t.p("G", 0.5, 0.5, -Z);
    t.s("GAMMA", "X"); t.s("X", "P"); t.s("P", "N"); t.s("N", "GAMMA"); t.s("GAMMA", "M"); t.s("M", "S");
    t.s("S_0", "GAMMA"); t.s("X", "R"); t.s("G", "M");
  }
  else if (ext == "oP1")
  {
    t.p("GAMMA", 0, 0, 0);
    t.p("X", 0.5, 0, 0);
    t.p("Z", 0, 0, 0.5);
    t.p("U", 0.5, 0, 0.5);
    t.p("Y", 0, 0.5, 0);
    t.p("S", 0.5, 0.5, 0);
    t.p("T", 0, 0.5, 0.5);
    t.p("R", 0.5, 0.5, 0.5);
    t.s("GAMMA", "X"); t.s("X", "S"); t.s("S", "Y"); t.s("Y", "GAMMA"); t.s("GAMMA", "Z"); t.s("Z", "U");
    t.s("U", "R"); t.s("R", "T"); t.s("T", "Z"); t.s("X", "U"); t.s("Y", "T"); t.s("S", "R");
  }
  else if (ext == "oF1")
  {
    double J = (1.0 + a * a / (b * b) - a * a / (c * c)) / 4.0;
    double H = (1.0 + a * a / (b * b) + a * a / (c * c)) / 4.0;
    t.p("GAMMA", 0, 0, 0);
    t.p("T", 1, 0.5, 0.5);
    t.p("Z", 0.5, 0.5, 0);
    t.p("Y", 0.5, 0, 0.5);
    t.p("SIGMA_0", 0, H, H);
    t.p("U_0", 1, 1.0 - H, 1.0 - H);
    t.p("A_0", 0.5, 0.5 + J, J);
    t.p("C_0", 0.5, 0.5 - J, 1.0 - J);
    t.p("L", 0.5, 0.5, 0.5);
    t.s("GAMMA", "Y"); t.s("Y", "T"); t.s("T", "Z"); t.s("Z", "GAMMA"); t.s("GAMMA", "SIGMA_0"); t.s("U_0", "T");
    t.s("Y", "C_0"); t.s("A_0", "Z"); t.s("GAMMA", "L");
  }
  else if (ext == "oF2")
  {
    double J = (1.0 + c * c / (a * a) - c * c / (b * b)) / 4.0;
    double K = (1.0 + c * c / (a * a) + c * c / (b * b)) / 4.0;
    t.p("GAMMA", 0, 0, 0);
    t.p("T", 0, 0.5, 0.5);
    t.p("Z", 0.5, 0.5, 1);
    t.p("Y", 0.5, 0, 0.5);
    t.p("LAMBDA_0", K, K, 0);
    t.p("Q_0", 1.0 - K, 1.0 - K, 1);
    t.p("G_0", 0.5 - J, 1.0 - J, 0.5);
    t.p("H_0", 0.5 + J, J, 0.5);
    t.p("L", 0.5, 0.5, 0.5);
    t.s("GAMMA", "T"); t.s("T", "Z"); t.s("Z", "Y"); t.s("Y", "GAMMA"); t.s("GAMMA", "LAMBDA_0"); t.s("Q_0", "Z");
    t.s("T", "G_0"); t.s("H_0", "Y"); t.s("GAMMA", "L");
  }
  else if (ext == "oF3")
  {
    double H = (1.0 + a * a / (b * b) - a * a / (c * c)) / 4.0;
    double K = (1.0 + b * b / (a * a) - b * b / (c * c)) / 4.0;
    double P = (1.0 + c * c / (b * b) - c * c / (a * a)) / 4.0;
    t.p("GAMMA", 0, 0, 0);
    t.p("T", 0, 0.5, 0.5);
    t.p("Z", 0.5, 0.5, 0);
    t.p("Y", 0.5, 0, 0.5);
    t.p("A_0", 0.5, 0.5 + H, H);
    t.p("C_0", 0.5, 0.5 - H, 1.0 - H);
    t.p("B_0", 0.5 + K, 0.5, K);
    t.p("D_0", 0.5 - K, 0.5, 1.0 - K);
    t.p("G_0", P, 0.5 + P, 0.5);
    t.p("H_0", 1.0 - P, 0.5 - P, 0.5);
    t.p("L", 0.5, 0.5, 0.5);
    t.s("GAMMA", "Y"); t.s("Y", "C_0"); t.s("A_0", "Z"); t.s("Z", "B_0"); t.s("D_0", "T"); t.s("T", "G_0");
    t.s("H_0", "Y"); t.s("T", "GAMMA"); t.s("GAMMA", "Z"); t.s("GAMMA", "L");
  }
  else if (ext == "oI1")
  {
    double Z = (1.0 + a * a / (c * c)) / 4.0;
    double H = (1.0 + b * b / (c * c)) / 4.0;
    double D = (b * b - a * a) / (4.0 * c * c);
    double N = (a * a + b * b) / (4.0 * c * c);
    t.p("GAMMA", 0, 0, 0);
    t.p("X", 0.5, 0.5, -0.5);
    t.p("S", 0.5, 0, 0);
    t.p("R", 0, 0.5, 0);
    t.p("T", 0, 0, 0.5);
    t.p("W", 0.25, 0.25, 0.25);
    t.p("SIGMA_0", -Z, Z, Z);
    t.p("F_2", Z, 1.0 - Z, -Z);
    t.p("Y_0", H, -H, H);
    t.p("U_0", 1.0 - H, H, -H);
    t.p("L_0", -N, N, 0.5 - D);
    t.p("M_0", N, -N, 0.5 + D);
    t.p("J_0", 0.5 - D, 0.5 + D, -N);
    t.s("GAMMA", "X"); t.s("X", "F_2"); t.s("SIGMA_0", "GAMMA"); t.s("GAMMA", "Y_0"); t.s("U_0", "X");
    t.s("GAMMA", "R"); t.s("R", "W"); t.s("W", "S"); t.s("S", "GAMMA"); t.s("GAMMA", "T"); t.s("T", "W");
  }
  else if (ext == "oI2")
  {
    double Z = (1.0 + b * b / (a * a)) / 4.0;
    double H = (1.0 + c * c / (a * a)) / 4.0;
    double D = (c * c - b * b) / (4.0 * a * a);
    double N = (b * b + c * c) / (4.0 * a * a);
    t.p("GAMMA", 0, 0, 0);
    t.p("X", -0.5, 0.5, 0.5);
    t.p("S", 0.5, 0, 0);
    t.p("R", 0, 0.5, 0);
    t.p("T", 0, 0, 0.5);
    t.p("W", 0.25, 0.25, 0.25);
    t.p("Y_0", Z, -Z, Z);
    t.p("U_2", -Z, Z, 1.0 - Z);
    t.p("LAMBDA_0", H, H, -H);
    t.p("G_2", -H, 1.0 - H, H);
    t.p("K", 0.5 - D, -N, N);
    t.p("K_2", 0.5 + D, N, -N);
    t.p("K_4", -N, 0.5 - D, 0.5 + D);
    t.s("GAMMA", "X"); t.s("X", "U_2"); t.s("Y_0", "GAMMA"); t.s("GAMMA", "LAMBDA_0"); t.s("G_2", "X");
    t.s("GAMMA", "R"); t.s("R", "W"); t.s("W", "S"); t.s("S", "GAMMA"); t.s("GAMMA", "T"); t.s("T", "W");
  }
  else if (ext == "oI3")
  {
    double Z = (1.0 + c * c / (b * b)) / 4.0;
    double Y = (1.0 + a * a / (b * b)) / 4.0;
    double D = (a * a - c * c) / (4.0 * b * b);
    double M = (c * c + a * a) / (4.0 * b * b);
    t.p("GAMMA", 0, 0, 0);
    t.p("X", 0.5, -0.5, 0.5);
    t.p("S", 0.5, 0, 0);
    t.p("R", 0, 0.5, 0);
    t.p("T", 0, 0, 0.5);
    t.p("W", 0.25, 0.25, 0.25);
    t.p("SIGMA_0", -Y, Y, Y);
    t.p("F_0", Y, -Y, 1.0 - Y);
    t.p("LAMBDA_0", Z, Z, -Z);
    t.p("G_0", 1.0 - Z, -Z, Z);
    t.p("V_0", M, 0.5 - D, -M);
    t.p("H_0", -M, 0.5 + D, M);
    t.p("H_2", 0.5 + D, -M, 0.5 - D);
    t.s("GAMMA", "X"); t.s("X", "F_0"); t.s("SIGMA_0", "GAMMA"); t.s("GAMMA", "LAMBDA_0"); t.s("G_0", "X");
    t.s("GAMMA", "R"); t.s("R", "W"); t.s("W", "S"); t.s("S", "GAMMA"); t.s("GAMMA", "T"); t.s("T", "W");
  }
  else if (ext == "oC1" || ext == "oA1")
  {
    double X = (ext == "oC1") ? (1.0 + a * a / (b * b)) / 4.0 : (1.0 + b * b / (c * c)) / 4.0;
    t.p("GAMMA", 0, 0, 0);
    t.p("Y", -0.5, 0.5, 0);
    t.p("T", -0.5, 0.5, 0.5);
    t.p("Z", 0, 0, 0.5);
    t.p("S", 0, 0.5, 0);
    t.p("R", 0, 0.5, 0.5);
    t.p("SIGMA_0", X, X, 0);
    t.p("C_0", -X, 1.0 - X, 0);
    t.p("A_0", X, X, 0.5);
    t.p("E_0", -X, 1.0 - X, 0.5);
    t.s("GAMMA", "Y"); t.s("Y", "C_0"); t.s("SIGMA_0", "GAMMA"); t.s("GAMMA", "Z"); t.s("Z", "A_0"); t.s("E_0", "T");
    t.s("T", "Y"); t.s("GAMMA", "S"); t.s("S", "R"); t.s("R", "Z"); t.s("Z", "T");
  }
  else if (ext == "oC2" || ext == "oA2")
  {
    double X = (ext == "oC2") ? (1.0 + b * b / (a * a)) / 4.0 : (1.0 + c * c / (b * b)) / 4.0;
    t.p("GAMMA", 0, 0, 0);
    t.p("Y", 0.5, 0.5, 0);
    t.p("T", 0.5, 0.5, 0.5);
    t.p("T_2", 0.5, 0.5, -0.5);
    t.p("Z", 0, 0, 0.5);
    t.p("Z_2", 0, 0, -0.5);
    t.p("S", 0, 0.5, 0);
    t.p("R", 0, 0.5, 0.5);
    t.p("R_2", 0, 0.5, -0.5);
    t.p("DELTA_0", -X, X, 0);
    t.p("F_0", X, 1.0 - X, 0);
    t.p("B_0", -X, X, 0.5);
    t.p("B_2", -X, X, -0.5);
    t.p("G_0", X, 1.0 - X, 0.5);
    t.p("G_2", X, 1.0 - X, -0.5);
    t.s("GAMMA", "Y"); t.s("Y", "F_0"); t.s("DELTA_0", "GAMMA"); t.s("GAMMA", "Z"); t.s("Z", "B_0"); t.s("G_0", "T");
    t.s("T", "Y"); t.s("GAMMA", "S"); t.s("S", "R"); t.s("R", "Z"); t.s("Z", "T");
  }
  else if (ext == "hP1" || ext == "hP2")
  {
    const double third = 1.0 / 3.0;
    t.p("GAMMA", 0, 0, 0);
    t.p("A", 0, 0, 0.5);
    t.p("K", third, third, 0);
    t.p("H", third, third, 0.5);
    t.p("H_2", third, third, -0.5);
    t.p("M", 0.5, 0, 0);
    t.p("L", 0.5, 0, 0.5);
    t.s("GAMMA", "M"); t.s("M", "K"); t.s("K", "GAMMA"); t.s("GAMMA", "A"); t.s("A", "L"); t.s("L", "H");
    t.s("H", "A"); t.s("L", "M"); t.s("H", "K");
    if (ext == "hP1") t.s("K", "H_2");
  }
  else if (ext == "hR1")
  {
    double D = a * a / (4.0 * c * c);
    double Y = 5.0 / 6.0 - 2.0 * D;
    double N = 1.0 / 3.0 + D;
    t.p("GAMMA", 0, 0, 0);
    t.p("T", 0.5, 0.5, 0.5);
    t.p("L", 0.5, 0, 0);
    t.p("L_2", 0, -0.5, 0);
    t.p("L_4", 0, 0, -0.5);
    t.p("F", 0.5, 0, 0.5);
    t.p("F_2", 0.5, 0.5, 0);
    t.p("S_0", N, -N, 0);
    t.p("S_2", 1.0 - N, 0, N);
    t.p("S_4", N, 0, -N);
    t.p("S_6", 1.0 - N, N, 0);
    t.p("H_0", 0.5, -1.0 + Y, 1.0 - Y);
    t.p("H_2", Y, 1.0 - Y, 0.5);
    t.p("H_4", Y, 0.5, 1.0 - Y);
    t.p("H_6", 0.5, 1.0 - Y, -1.0 + Y);
    t.p("M_0", N, -1.0 + Y, N);
    t.p("M_2", 1.0 - N, 1.0 - Y, 1.0 - N);
    t.p("M_4", Y, N, N);
    t.p("M_6", 1.0 - N, 1.0 - N, 1.0 - Y);
    t.p("M_8", N, N, -1.0 + Y);
    t.s("GAMMA", "T"); t.s("T", "H_2"); t.s("H_0", "L"); t.s("L", "GAMMA"); t.s("GAMMA", "S_0"); t.s("S_2", "F");
    t.s("F", "GAMMA");
  }
  else if (ext == "hR2")
  {
    double Z = 1.0 / 6.0 - c * c / (9.0 * a * a);
    double H = 0.5 - 2.0 * Z;
    double N = 0.5 + Z;
    t.p("GAMMA", 0, 0, 0);
    t.p("T", 0.5, -0.5, 0.5);
    t.p("P_0", H, -1.0 + H, H);
    t.p("P_2", H, H, H);
    t.p("R_0", 1.0 - H, -H, -H);
    t.p("M", 1.0 - N, -N, 1.0 - N);
    t.p("M_2", N, -1.0 + N, -1.0 + N);
    t.p("L", 0.5, 0, 0);
    t.p("F", 0.5, -0.5, 0);
    t.s("GAMMA", "L"); t.s("L", "T"); t.s("T", "P_0"); t.s("P_2", "GAMMA"); t.s("GAMMA", "F");
  }
  else if (ext == "mP1")
  {
    double Y = (1.0 + a / c * cosBeta) / (2.0 * sin2Beta);
    double N = 0.5 + Y * c * cosBeta / a;
    t.p("GAMMA", 0, 0, 0);
    t.p("Z", 0, 0.5, 0);
    t.p("B", 0, 0, 0.5);
    t.p("B_2", 0, 0, -0.5);
    t.p("Y", 0.5, 0, 0);
    t.p("Y_2", -0.5, 0, 0);
    t.p("C", 0.5, 0.5, 0);
    t.p("C_2", -0.5, 0.5, 0);
    t.p("D", 0, 0.5, 0.5);
    t.p("D_2", 0, 0.5, -0.5);
    t.p("A", -0.5, 0, 0.5);
    t.p("E", -0.5, 0.5, 0.5);
    t.p("H", -Y, 0, 1.0 - N);
    t.p("H_2", -1.0 + Y, 0, N);
    t.p("H_4", -Y, 0, -N);
    t.p("M", -Y, 0.5, 1.0 - N);
    t.p("M_2", -1.0 + Y, 0.5, N);
    t.p("M_4", -Y, 0.5, -N);
    t.s("GAMMA", "Z"); t.s("Z", "D"); t.s("D", "B"); t.s("B", "GAMMA"); t.s("GAMMA", "A"); t.s("A", "E");
    t.s("E", "Z"); t.s("Z", "C_2"); t.s("C_2", "Y_2"); t.s("Y_2", "GAMMA");
  }
  else if (ext == "mC1")
  {
    double kZ = (2.0 + a / c * cosBeta) / (4.0 * sin2Beta);
    double kH = 0.5 - 2.0 * kZ * c * cosBeta / a;
    double kS = 0.75 - b * b / (4.0 * a * a * sin2Beta);
    double kP = kS - (0.75 - kS) * a * cosBeta / c;
    t.p("GAMMA", 0, 0, 0);
    t.p("Y_2", -0.5, 0.5, 0);
    t.p("Y_4", 0.5, -0.5, 0);
    t.p("A", 0, 0, 0.5);
    t.p("M_2", -0.5, 0.5, 0.5);
    t.p("V", 0.5, 0, 0);
    t.p("V_2", 0, 0.5, 0);
    t.p("L_2", 0, 0.5, 0.5);
    t.p("C", 1.0 - kS, 1.0 - kS, 0);
    t.p("C_2", -1.0 + kS, kS, 0);
    t.p("C_4", kS, -1.0 + kS, 0);
    t.p("D", -1.0 + kP, kP, 0.5);
    t.p("D_2", 1.0 - kP, 1.0 - kP, 0.5);
    t.p("E", -1.0 + kZ, 1.0 - kZ, 1.0 - kH);
    t.p("E_2", -kZ, kZ, kH);
    t.p("E_4", kZ, -kZ, 1.0 - kH);
    t.s("GAMMA", "C"); t.s("C_2", "Y_2"); t.s("Y_2", "GAMMA"); t.s("GAMMA", "M_2"); t.s("M_2", "D"); t.s("D_2", "A");
    t.s("A", "GAMMA"); t.s("L_2", "GAMMA"); t.s("GAMMA", "V_2");
  }
  else if (ext == "mC2")
  {
    double kZ = (a * a / (b * b) + (1.0 + a / c * cosBeta) / sin2Beta) / 4.0;
    double kM = (1.0 + a * a / (b * b)) / 4.0;
    double kD = -a * c * cosBeta / (2.0 * b * b);
    double kX = 0.5 - 2.0 * kZ * c * cosBeta / a;
    double kP = 1.0 + kZ - 2.0 * kM;
    double kS = kX - 2.0 * kD;
    t.p("GAMMA", 0, 0, 0);
    t.p("Y", 0.5, 0.5, 0);
    t.p("A", 0, 0, 0.5);
    t.p("M", 0.5, 0.5, 0.5);
    t.p("V_2", 0, 0.5, 0);
    t.p("L_2", 0, 0.5, 0.5);
    t.p("F", -1.0 + kP, 1.0 - kP, 1.0 - kS);
    t.p("F_2", 1.0 - kP, kP, kS);
    t.p("F_4", kP, 1.0 - kP, 1.0 - kS);
    t.p("H", -kZ, kZ, kX);
    t.p("H_2", kZ, 1.0 - kZ, 1.0 - kX);
    t.p("H_4", kZ, -kZ, 1.0 - kX);
    t.p("G", -kM, kM, kD);
    t.p("G_2", kM, 1.0 - kM, -kD);
    t.p("G_4", kM, -kM, -kD);
    t.p("G_6", 1.0 - kM, kM, kD);
    t.s("GAMMA", "Y"); t.s("Y", "M"); t.s("M", "A"); t.s("A", "GAMMA"); t.s("L_2", "GAMMA"); t.s("GAMMA", "V_2");
  }
  else if (ext == "mC3")
  {
    double kZ = (a * a / (b * b) + (1.0 + a / c * cosBeta) / sin2Beta) / 4.0;
    double kR = 1.0 - kZ * b * b / (a * a);
    double kE = 0.5 - 2.0 * kZ * c * cosBeta / a;
    double kF = kE / 2.0 + a * a / (4.0 * b * b) + a * c * cosBeta / (2.0 * b * b);
    double kU = 2.0 * kF - kZ;
    double kW = c / (2.0 * a * cosBeta) * (1.0 - 4.0 * kU + a * a * sin2Beta / (b * b));
    double kD = -0.25 + kW / 2.0 - kZ * c * cosBeta / a;
    t.p("GAMMA", 0, 0, 0);
    t.p("Y", 0.5, 0.5, 0);
    t.p("A", 0, 0, 0.5);
    t.p("M_2", -0.5, 0.5, 0.5);
    t.p("V", 0.5, 0, 0);
    t.p("V_2", 0, 0.5, 0);
    t.p("L_2", 0, 0.5, 0.5);
    t.p("I", -1.0 + kR, kR, 0.5);
    t.p("I_2", 1.0 - kR, 1.0 - kR, 0.5);
    t.p("K", -kU, kU, kW);
    t.p("K_2", -1.0 + kU, 1.0 - kU, 1.0 - kW);
    t.p("K_4", 1.0 - kU, kU, kW);
    t.p("H", -kZ, kZ, kE);
    t.p("H_2", kZ, 1.0 - kZ, 1.0 - kE);
    t.p("H_4", kZ, -kZ, 1.0 - kE);
    t.p("N", -kF, kF, kD);
    t.p("N_2", kF, 1.0 - kF, -kD);
    t.p("N_4", kF, -kF, -kD);
    t.p("N_6", 1.0 - kF, kF, kD);
    t.s("GAMMA", "A"); t.s("A", "I_2"); t.s("I", "M_2"); t.s("M_2", "GAMMA"); t.s("GAMMA", "Y"); t.s("L_2", "GAMMA");
    t.s("GAMMA", "V_2");
  }
  else if (ext == "aP2")
  {
    t.p("GAMMA", 0, 0, 0);
    t.p("Z", 0, 0, 0.5);
    t.p("Y", 0, 0.5, 0);
    t.p("X", 0.5, 0, 0);
    t.p("V", 0.5, 0.5, 0);
    t.p("U", 0.5, 0, 0.5);
    t.p("T", 0, 0.5, 0.5);
    t.p("R", 0.5, 0.5, 0.5);
    t.s("GAMMA", "X"); t.s("Y", "GAMMA"); t.s("GAMMA", "Z"); t.s("R", "GAMMA"); t.s("GAMMA", "T"); t.s("U", "GAMMA");
    t.s("GAMMA", "V");
  }
  else if (ext == "aP3")
  {
    t.p("GAMMA", 0, 0, 0);
    t.p("Z", 0, 0, 0.5);
    t.p("Y", 0, 0.5, 0);
    t.p("Y_2", 0, -0.5, 0);
    t.p("X", 0.5, 0, 0);
    t.p("V_2", 0.5, -0.5, 0);
    t.p("U_2", -0.5, 0, 0.5);
    t.p("T_2", 0, -0.5, 0.5);
    t.p("R_2", -0.5, -0.5, 0.5);
    t.s("GAMMA", "X"); t.s("Y", "GAMMA"); t.s("GAMMA", "Z"); t.s("R_2", "GAMMA"); t.s("GAMMA", "T_2");
    t.s("U_2", "GAMMA"); t.s("GAMMA", "V_2");
  }

  return t;
}

// Augment the path with the corresponding -k segments (used when the structure has neither inversion
// nor assumed time-reversal symmetry).
void augmentWithTimeReversal(Tables& t)
{
  std::vector<SKKVector> primed;
  for (const SKKVector& kv : t.points)
  {
    if (kv.label == "GAMMA") continue;
    primed.push_back(SKKVector{kv.label + "'", double3(-kv.coordinates.x, -kv.coordinates.y, -kv.coordinates.z)});
  }
  std::vector<std::pair<std::string, std::string>> original = t.segments;
  for (const auto& [from, to] : original)
  {
    std::string newFrom = (from == "GAMMA") ? from : from + "'";
    std::string newTo = (to == "GAMMA") ? to : to + "'";
    t.segments.emplace_back(newFrom, newTo);
  }
  t.points.insert(t.points.end(), primed.begin(), primed.end());
}

}  // namespace

std::optional<double3> SKBandPath::coordinatesForLabel(std::string_view label) const
{
  for (const SKKVector& kv : kPoints)
  {
    if (kv.label == label) return kv.coordinates;
  }
  return std::nullopt;
}

double3 SKBandPath::reciprocalFractionalInAnalyzedCell(const double3& primitiveReciprocalFractional) const
{
  // Reciprocal fractional coordinates transform contravariantly with respect to the real-space change of
  // basis. With (a_p b_p c_p) = (a b c) Q and (analyzed) = (a b c) M, one has (a_p) = (analyzed) M^{-1} Q, so
  //   f_analyzed = ( (M^{-1} Q)^T )^{-1} f_primitive = M^T (Q^T)^{-1} f_primitive.
  double3 fConventional = conventionalToPrimitiveTransformation.transpose().inverse() * primitiveReciprocalFractional;
  return conventionalToAnalyzedTransformation.transpose() * fConventional;
}

std::optional<SKBandPath> SKBrillouinZonePath(double3x3 conventionalCell, std::size_t spaceGroupNumber,
                                              bool hasInversionSymmetry, bool withTimeReversal, double threshold)
{
  if (spaceGroupNumber < 1 || spaceGroupNumber > 230) return std::nullopt;

  char family = crystalFamilyLetter(spaceGroupNumber);
  char centring = centringLetter(spaceGroupNumber);
  std::string bravais = std::string{family} + std::string{centring};

  CellParameters params = cellParameters(conventionalCell);

  std::string ext;
  double3x3 usedConventionalCell = conventionalCell;
  if (bravais == "aP")
  {
    auto [extAP, standardizedCell] = standardizeTriclinic(conventionalCell);
    ext = extAP;
    usedConventionalCell = standardizedCell;
    params = cellParameters(usedConventionalCell);
  }
  else
  {
    ext = extendedBravaisSymbol(bravais, spaceGroupNumber, params, threshold);
  }
  if (ext.empty()) return std::nullopt;

  Tables tables = buildTables(ext, params);
  if (tables.points.empty()) return std::nullopt;

  bool augmented = !hasInversionSymmetry && !withTimeReversal;
  if (augmented) augmentWithTimeReversal(tables);

  SKBandPath result;
  result.bravaisLattice = bravais;
  result.extendedBravaisLattice = ext;
  result.hasInversionSymmetry = hasInversionSymmetry;
  result.timeReversalAugmented = augmented;
  result.kPoints = std::move(tables.points);
  result.segments = std::move(tables.segments);

  double3x3 P = primitiveTransformation(bravais);
  result.primitiveTransformationMatrix = P;
  result.conventionalCell = usedConventionalCell;
  result.primitiveCell = usedConventionalCell * P;
  result.primitiveReciprocalCell = (2.0 * std::numbers::pi) * reciprocalLattice(result.primitiveCell);

  // Q maps the (possibly re-standardized) primitive cell back onto the conventional cell that was passed in;
  // it equals P except for the triclinic case where the cell is re-standardized.
  result.conventionalToPrimitiveTransformation = conventionalCell.inverse() * result.primitiveCell;
  result.conventionalToAnalyzedTransformation = double3x3::identity();

  return result;
}

std::optional<SKBandPath> SKBrillouinZonePath(double3x3 unitCell,
                                              std::vector<std::tuple<double3, std::size_t, double>> atoms,
                                              bool allowPartialOccupancies, double symmetryPrecision,
                                              bool withTimeReversal, double threshold)
{
  std::optional<SKSpaceGroup::FoundSpaceGroupInfo> found =
      SKSpaceGroup::findSpaceGroup(unitCell, atoms, allowPartialOccupancies, symmetryPrecision);
  if (!found) return std::nullopt;

  SKSpaceGroup spaceGroup(found->HallNumber);
  std::size_t spaceGroupNumber = spaceGroup.spaceGroupSetting().number();
  std::size_t pointGroupNumber = spaceGroup.spaceGroupSetting().pointGroupNumber();

  bool hasInversionSymmetry = false;
  if (pointGroupNumber < SKPointGroup::pointGroupData.size())
  {
    hasInversionSymmetry = SKPointGroup::pointGroupData[pointGroupNumber].centrosymmetric();
  }

  // Build a canonical conventional cell (columns are lattice vectors) from the standardized cell parameters.
  const SKSymmetryCell& c = found->cell;
  double3x3 conventionalCell = cellFromParameters(c.a(), c.b(), c.c(), c.alpha(), c.beta(), c.gamma());

  std::optional<SKBandPath> path =
      SKBrillouinZonePath(conventionalCell, spaceGroupNumber, hasInversionSymmetry, withTimeReversal, threshold);
  if (path)
  {
    // transformationMatrix maps the conventional cell onto the analyzed unit cell: unitCell = conventional * M.
    path->conventionalToAnalyzedTransformation = found->transformationMatrix;
  }
  return path;
}
