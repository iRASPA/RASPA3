module;

module cmap_potential;

import std;

import archive;
import units;
import double3;
import double3x3;

namespace
{
// The derivatives at the nodes of the C^2 periodic cubic spline through equidistant periodic data
// y_0, ..., y_{n-1} (spacing h): m_{i-1} + 4 m_i + m_{i+1} = 3 (y_{i+1} - y_{i-1}) / h. The cyclic tridiagonal
// system is solved with the Sherman-Morrison correction of the plain tridiagonal solve.
std::vector<double> periodicSplineDerivatives(std::span<const double> y, double h)
{
  const std::size_t n = y.size();
  std::vector<double> rhs(n);
  for (std::size_t i = 0; i < n; ++i)
  {
    rhs[i] = 3.0 * (y[(i + 1) % n] - y[(i + n - 1) % n]) / h;
  }
  if (n < 3)
  {
    // degenerate (n = 1, 2): the spline is flat
    return std::vector<double>(n, 0.0);
  }

  const auto tridiagonal = [n](std::span<const double> diagonal, std::span<const double> r)
  {
    std::vector<double> x(n), gamma(n);
    double beta = diagonal[0];
    x[0] = r[0] / beta;
    for (std::size_t j = 1; j < n; ++j)
    {
      gamma[j] = 1.0 / beta;
      beta = diagonal[j] - gamma[j];
      x[j] = (r[j] - x[j - 1]) / beta;
    }
    for (std::size_t j = n - 1; j-- > 0;)
    {
      x[j] -= gamma[j + 1] * x[j + 1];
    }
    return x;
  };

  const double corner = -4.0;
  std::vector<double> diagonal(n, 4.0);
  diagonal[0] = 4.0 - corner;
  diagonal[n - 1] = 4.0 - 1.0 / corner;
  std::vector<double> x = tridiagonal(diagonal, rhs);
  std::vector<double> u(n, 0.0);
  u[0] = corner;
  u[n - 1] = 1.0;
  const std::vector<double> z = tridiagonal(diagonal, u);
  const double factor = (x[0] + x[n - 1] / corner) / (1.0 + z[0] + z[n - 1] / corner);
  for (std::size_t i = 0; i < n; ++i) x[i] -= factor * z[i];
  return x;
}

// alpha = A F A^T of the bicubic patch on the unit square (Wikipedia 'Bicubic interpolation')
constexpr std::array<std::array<double, 4>, 4> patchMatrix{
    {{1.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 1.0, 0.0}, {-3.0, 3.0, -2.0, -1.0}, {2.0, -2.0, 1.0, 1.0}}};

void addStrain(double3x3 &strain, const double3 &d, const double3 &g)
{
  strain.ax += d.x * g.x;
  strain.bx += d.y * g.x;
  strain.cx += d.z * g.x;
  strain.ay += d.x * g.y;
  strain.by += d.y * g.y;
  strain.cy += d.z * g.y;
  strain.az += d.x * g.z;
  strain.bz += d.y * g.z;
  strain.cz += d.z * g.z;
}
}  // namespace

CMAPMap::CMAPMap(std::string name, std::size_t resolution, std::vector<double> energies)
    : name(std::move(name)), resolution(resolution), energies(std::move(energies))
{
  if (this->resolution < 3)
  {
    throw std::runtime_error(
        std::format("[CMAP]: map '{}' has a resolution of {}, at least 3 grid points per angle are required\n",
                    this->name, this->resolution));
  }
  if (this->energies.size() != this->resolution * this->resolution)
  {
    throw std::runtime_error(std::format("[CMAP]: map '{}' has {} energies, expected {} ({} x {})\n", this->name,
                                         this->energies.size(), this->resolution * this->resolution,
                                         this->resolution, this->resolution));
  }
  buildSpline();
}

void CMAPMap::buildSpline()
{
  const std::size_t n = resolution;
  coefficients.assign(n * n, std::array<double, 16>{});
  if (n == 0) return;
  const double h = 2.0 * std::numbers::pi / static_cast<double>(n);

  const auto at = [n](std::size_t i, std::size_t j) { return (i % n) * n + (j % n); };

  // first derivatives along the two angles, the cross derivative from a spline through dE/dphi along psi
  std::vector<double> dPhi(n * n), dPsi(n * n), dPhiPsi(n * n), line(n);
  for (std::size_t j = 0; j < n; ++j)
  {
    for (std::size_t i = 0; i < n; ++i) line[i] = energies[at(i, j)];
    const std::vector<double> m = periodicSplineDerivatives(line, h);
    for (std::size_t i = 0; i < n; ++i) dPhi[at(i, j)] = m[i];
  }
  for (std::size_t i = 0; i < n; ++i)
  {
    for (std::size_t j = 0; j < n; ++j) line[j] = energies[at(i, j)];
    const std::vector<double> m = periodicSplineDerivatives(line, h);
    for (std::size_t j = 0; j < n; ++j) dPsi[at(i, j)] = m[j];
  }
  for (std::size_t i = 0; i < n; ++i)
  {
    for (std::size_t j = 0; j < n; ++j) line[j] = dPhi[at(i, j)];
    const std::vector<double> m = periodicSplineDerivatives(line, h);
    for (std::size_t j = 0; j < n; ++j) dPhiPsi[at(i, j)] = m[j];
  }

  // the bicubic patch of every cell from the values, gradients and cross derivatives at its four corners
  for (std::size_t i = 0; i < n; ++i)
  {
    for (std::size_t j = 0; j < n; ++j)
    {
      std::array<std::array<double, 4>, 4> F{};
      for (std::size_t a = 0; a < 2; ++a)
      {
        for (std::size_t b = 0; b < 2; ++b)
        {
          const std::size_t corner = at(i + a, j + b);
          F[a][b] = energies[corner];
          F[a][b + 2] = h * dPsi[corner];
          F[a + 2][b] = h * dPhi[corner];
          F[a + 2][b + 2] = h * h * dPhiPsi[corner];
        }
      }
      // alpha = A F A^T
      std::array<std::array<double, 4>, 4> AF{};
      for (std::size_t r = 0; r < 4; ++r)
        for (std::size_t c = 0; c < 4; ++c)
          for (std::size_t k = 0; k < 4; ++k) AF[r][c] += patchMatrix[r][k] * F[k][c];
      std::array<double, 16> &cell = coefficients[at(i, j)];
      for (std::size_t k = 0; k < 4; ++k)
        for (std::size_t l = 0; l < 4; ++l)
        {
          double value = 0.0;
          for (std::size_t c = 0; c < 4; ++c) value += AF[k][c] * patchMatrix[l][c];
          cell[k * 4 + l] = value;
        }
    }
  }
}

CMAPMap::Evaluation CMAPMap::evaluate(double phi, double psi) const
{
  Evaluation result{};
  const std::size_t n = resolution;
  if (n == 0 || coefficients.size() != n * n) return result;
  const double h = 2.0 * std::numbers::pi / static_cast<double>(n);
  const double size = static_cast<double>(n);

  const auto locate = [&](double angle, std::size_t &index, double &fraction)
  {
    double x = (angle + std::numbers::pi) / h;
    x -= std::floor(x / size) * size;
    if (x < 0.0) x += size;
    if (x >= size) x -= size;
    double whole = std::floor(x);
    if (whole >= size) whole = size - 1.0;
    if (whole < 0.0) whole = 0.0;
    index = static_cast<std::size_t>(whole);
    fraction = x - whole;
  };
  std::size_t i, j;
  double t, u;
  locate(phi, i, t);
  locate(psi, j, u);

  const std::array<double, 16> &c = coefficients[i * n + j];
  const std::array<double, 4> tPowers{1.0, t, t * t, t * t * t};
  const std::array<double, 4> uPowers{1.0, u, u * u, u * u * u};
  double e{}, dt{}, du{}, dtt{}, dtu{}, duu{};
  for (std::size_t k = 0; k < 4; ++k)
  {
    for (std::size_t l = 0; l < 4; ++l)
    {
      const double coefficient = c[k * 4 + l];
      e += coefficient * tPowers[k] * uPowers[l];
      if (k > 0) dt += static_cast<double>(k) * coefficient * tPowers[k - 1] * uPowers[l];
      if (l > 0) du += static_cast<double>(l) * coefficient * tPowers[k] * uPowers[l - 1];
      if (k > 1) dtt += static_cast<double>(k * (k - 1)) * coefficient * tPowers[k - 2] * uPowers[l];
      if (k > 0 && l > 0) dtu += static_cast<double>(k * l) * coefficient * tPowers[k - 1] * uPowers[l - 1];
      if (l > 1) duu += static_cast<double>(l * (l - 1)) * coefficient * tPowers[k] * uPowers[l - 2];
    }
  }
  result.energy = e;
  result.dPhi = dt / h;
  result.dPsi = du / h;
  result.dPhiPhi = dtt / (h * h);
  result.dPhiPsi = dtu / (h * h);
  result.dPsiPsi = duu / (h * h);
  return result;
}

void CMAPMap::scaleEnergy(double factor)
{
  for (double &energy : energies) energy *= factor;
  for (std::array<double, 16> &cell : coefficients)
    for (double &coefficient : cell) coefficient *= factor;
}

std::string CMAPMap::print() const
{
  const auto [minimum, maximum] = std::ranges::minmax_element(energies);
  return std::format("CMAP map '{}': {} x {} grid, energies from {:.3f} to {:.3f} K", name, resolution, resolution,
                     energies.empty() ? 0.0 : *minimum * Units::EnergyToKelvin,
                     energies.empty() ? 0.0 : *maximum * Units::EnergyToKelvin);
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const CMAPMap &m)
{
  archive << m.versionNumber;
  archive << m.name;
  archive << m.resolution;
  archive << m.energies;
#if DEBUG_ARCHIVE
  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif
  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, CMAPMap &m)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > m.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'CMAPMap' at line {} in file {}\n", location.line(),
                                         location.file_name()));
  }
  archive >> m.name;
  archive >> m.resolution;
  archive >> m.energies;
#if DEBUG_ARCHIVE
  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("CMAPMap: Error in binary restart\n"));
  }
#endif
  m.buildSpline();
  return archive;
}

std::string CMAPPotential::print() const
{
  return std::format("CMAP ({}, {}, {}, {}, {}) map {}", identifiers[0], identifiers[1], identifiers[2],
                     identifiers[3], identifiers[4], mapIndex);
}

double CMAPPotential::dihedralAngle(const double3 &posA, const double3 &posB, const double3 &posC,
                                    const double3 &posD)
{
  const double3 b1 = posB - posA;
  const double3 b2 = posC - posB;
  const double3 b3 = posD - posC;
  const double3 m = double3::cross(b1, b2);
  const double3 n = double3::cross(b2, b3);
  return std::atan2(std::sqrt(double3::dot(b2, b2)) * double3::dot(b1, n), double3::dot(m, n));
}

std::pair<double, std::array<double3, 4>> CMAPPotential::dihedralAngleAndGradient(const double3 &posA,
                                                                                  const double3 &posB,
                                                                                  const double3 &posC,
                                                                                  const double3 &posD)
{
  const double3 b1 = posB - posA;
  const double3 b2 = posC - posB;
  const double3 b3 = posD - posC;
  const double3 m = double3::cross(b1, b2);
  const double3 n = double3::cross(b2, b3);
  const double b2Squared = double3::dot(b2, b2);
  const double b2Length = std::sqrt(b2Squared);
  const double phi = std::atan2(b2Length * double3::dot(b1, n), double3::dot(m, n));

  // Blondel & Karplus (1996); the central atoms by translation and rotation invariance
  const double mSquared = std::max(double3::dot(m, m), 1.0e-20);
  const double nSquared = std::max(double3::dot(n, n), 1.0e-20);
  const double3 gradientA = (-b2Length / mSquared) * m;
  const double3 gradientD = (b2Length / nSquared) * n;
  const double p = double3::dot(b1, b2) / b2Squared;
  const double q = double3::dot(b3, b2) / b2Squared;
  const double3 gradientB = -gradientA - p * gradientA + q * gradientD;
  const double3 gradientC = -gradientD + p * gradientA - q * gradientD;
  return {phi, {gradientA, gradientB, gradientC, gradientD}};
}

double CMAPPotential::calculateEnergy(const CMAPMap &map, const double3 &posA, const double3 &posB,
                                      const double3 &posC, const double3 &posD, const double3 &posE) const
{
  const double phi = dihedralAngle(posA, posB, posC, posD);
  const double psi = dihedralAngle(posB, posC, posD, posE);
  return map.evaluate(phi, psi).energy;
}

std::tuple<double, std::array<double3, 5>, double3x3> CMAPPotential::potentialEnergyGradientStrain(
    const CMAPMap &map, const double3 &posA, const double3 &posB, const double3 &posC, const double3 &posD,
    const double3 &posE) const
{
  const auto [phi, dPhi] = dihedralAngleAndGradient(posA, posB, posC, posD);
  const auto [psi, dPsi] = dihedralAngleAndGradient(posB, posC, posD, posE);
  const CMAPMap::Evaluation value = map.evaluate(phi, psi);

  std::array<double3, 5> gradient{};
  for (std::size_t k = 0; k < 4; ++k)
  {
    gradient[k] += value.dPhi * dPhi[k];
    gradient[k + 1] += value.dPsi * dPsi[k];
  }

  double3x3 strain{};
  addStrain(strain, posA - posB, gradient[0]);
  addStrain(strain, posC - posB, gradient[2]);
  addStrain(strain, posD - posB, gradient[3]);
  addStrain(strain, posE - posB, gradient[4]);

  return {value.energy, gradient, strain};
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const CMAPPotential &p)
{
  archive << p.versionNumber;
  archive << p.identifiers;
  archive << p.mapIndex;
#if DEBUG_ARCHIVE
  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif
  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, CMAPPotential &p)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > p.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'CMAPPotential' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }
  archive >> p.identifiers;
  archive >> p.mapIndex;
#if DEBUG_ARCHIVE
  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("CMAPPotential: Error in binary restart\n"));
  }
#endif
  return archive;
}
