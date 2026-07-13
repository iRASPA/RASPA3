module;

module minimization_internal_coordinate_hessian;

import std;

import double3;

namespace
{
using Minimization::CoordinateBlock;

double component(const double3 &v, std::size_t axis) { return (&v.x)[axis]; }

CoordinateBlock scaledBlock(const CoordinateBlock &block, double scale)
{
  CoordinateBlock result{};
  for (std::size_t a = 0; a < 3; ++a)
  {
    for (std::size_t b = 0; b < 3; ++b)
    {
      result[a][b] = scale * block[a][b];
    }
  }
  return result;
}

CoordinateBlock transposed(const CoordinateBlock &block)
{
  CoordinateBlock result{};
  for (std::size_t a = 0; a < 3; ++a)
  {
    for (std::size_t b = 0; b < 3; ++b)
    {
      result[a][b] = block[b][a];
    }
  }
  return result;
}

CoordinateBlock sumBlocks(const CoordinateBlock &first, const CoordinateBlock &second, double scaleFirst = 1.0,
                          double scaleSecond = 1.0)
{
  CoordinateBlock result{};
  for (std::size_t a = 0; a < 3; ++a)
  {
    for (std::size_t b = 0; b < 3; ++b)
    {
      result[a][b] = scaleFirst * first[a][b] + scaleSecond * second[a][b];
    }
  }
  return result;
}
}  // namespace

Minimization::InternalCoordinate Minimization::distanceInternalCoordinate(const double3 &posI, const double3 &posJ)
{
  InternalCoordinate coordinate{};
  coordinate.numAtoms = 2;

  const double3 difference = posI - posJ;
  const double r = std::sqrt(double3::dot(difference, difference));
  const double3 unitVector = difference / r;

  coordinate.value = r;
  coordinate.gradient[0] = unitVector;
  coordinate.gradient[1] = -unitVector;

  // d²r/dx dy = (delta_ab - u_a u_b) / r on the same atom, negated across atoms.
  CoordinateBlock block{};
  for (std::size_t a = 0; a < 3; ++a)
  {
    for (std::size_t b = 0; b < 3; ++b)
    {
      block[a][b] = ((a == b ? 1.0 : 0.0) - component(unitVector, a) * component(unitVector, b)) / r;
    }
  }
  coordinate.hessian[0][0] = block;
  coordinate.hessian[1][1] = block;
  coordinate.hessian[0][1] = scaledBlock(block, -1.0);
  coordinate.hessian[1][0] = coordinate.hessian[0][1];

  return coordinate;
}

Minimization::InternalCoordinate Minimization::angleInternalCoordinate(const double3 &posA, const double3 &posB,
                                                                       const double3 &posC)
{
  InternalCoordinate coordinate{};
  coordinate.numAtoms = 3;

  const double3 RabVector = posA - posB;
  const double rab = std::sqrt(double3::dot(RabVector, RabVector));
  const double3 Rab = RabVector / rab;

  const double3 RbcVector = posC - posB;
  const double rbc = std::sqrt(double3::dot(RbcVector, RbcVector));
  const double3 Rbc = RbcVector / rbc;

  const double cosTheta = std::clamp(double3::dot(Rab, Rbc), -1.0, 1.0);
  const double sinTheta = std::max(1.0e-8, std::sqrt(1.0 - cosTheta * cosTheta));
  coordinate.value = std::acos(cosTheta);

  // First derivatives of cos(theta).
  const double3 dtA = (Rbc - cosTheta * Rab) / rab;
  const double3 dtC = (Rab - cosTheta * Rbc) / rbc;
  const double3 dtB = -(dtA + dtC);

  // Second derivatives of cos(theta), same expressions as the bend Hessian with DF = 1.
  CoordinateBlock AA{};
  CoordinateBlock CC{};
  CoordinateBlock AC{};
  for (std::size_t a = 0; a < 3; ++a)
  {
    for (std::size_t b = 0; b < 3; ++b)
    {
      AA[a][b] = (cosTheta * component(Rab, a) * component(Rab, b) / rab - component(dtA, a) * component(Rab, b) -
                  component(dtA, b) * component(Rab, a)) /
                     rab -
                 (a == b ? cosTheta / (rab * rab) : 0.0);
      CC[a][b] = (cosTheta * component(Rbc, a) * component(Rbc, b) / rbc - component(dtC, a) * component(Rbc, b) -
                  component(dtC, b) * component(Rbc, a)) /
                     rbc -
                 (a == b ? cosTheta / (rbc * rbc) : 0.0);
      AC[a][b] = (cosTheta * component(Rab, a) * component(Rbc, b) - component(Rab, a) * component(Rab, b) -
                  component(Rbc, a) * component(Rbc, b) + (a == b ? 1.0 : 0.0)) /
                 (rab * rbc);
    }
  }

  // Remaining blocks follow from translational invariance (block rows sum to zero).
  const CoordinateBlock CA = transposed(AC);
  const CoordinateBlock AB = sumBlocks(AA, AC, -1.0, -1.0);
  const CoordinateBlock CB = sumBlocks(CC, CA, -1.0, -1.0);
  const CoordinateBlock BB = sumBlocks(sumBlocks(AA, CC), sumBlocks(AC, CA));

  const std::array<double3, 3> cosGradients = {dtA, dtB, dtC};
  const std::array<std::array<CoordinateBlock, 3>, 3> cosHessian = {
      {{AA, AB, AC}, {transposed(AB), BB, transposed(CB)}, {CA, CB, CC}}};

  // Chain rule from cos(theta) to theta: dtheta/dc = -1/sin, d²theta/dc² = -cos/sin³.
  const double firstFactor = -1.0 / sinTheta;
  const double secondFactor = -cosTheta / (sinTheta * sinTheta * sinTheta);
  for (std::size_t i = 0; i < 3; ++i)
  {
    coordinate.gradient[i] = firstFactor * cosGradients[i];
    for (std::size_t j = 0; j < 3; ++j)
    {
      for (std::size_t a = 0; a < 3; ++a)
      {
        for (std::size_t b = 0; b < 3; ++b)
        {
          coordinate.hessian[i][j][a][b] = firstFactor * cosHessian[i][j][a][b] +
                                           secondFactor * component(cosGradients[i], a) *
                                               component(cosGradients[j], b);
        }
      }
    }
  }

  return coordinate;
}

Minimization::InternalCoordinate Minimization::dihedralInternalCoordinate(const double3 &posA, const double3 &posB,
                                                                          const double3 &posC, const double3 &posD)
{
  InternalCoordinate coordinate{};
  coordinate.numAtoms = 4;

  const double3 Dab = posA - posB;
  const double3 DcbRaw = posC - posB;
  const double rbc = std::sqrt(double3::dot(DcbRaw, DcbRaw));
  const double3 Dcb = DcbRaw / rbc;
  const double3 Ddc = posD - posC;

  const double dotAB = double3::dot(Dab, Dcb);
  const double dotCD = double3::dot(Ddc, Dcb);

  double3 dr = Dab - dotAB * Dcb;
  const double r = std::sqrt(double3::dot(dr, dr));
  dr /= r;

  double3 ds = Ddc - dotCD * Dcb;
  const double s = std::sqrt(double3::dot(ds, ds));
  ds /= s;

  const double cosPhi = std::clamp(double3::dot(dr, ds), -1.0, 1.0);
  coordinate.value = cosPhi;

  const double3 Pb = double3::cross(Dab, Dcb);
  const double3 Pc = double3::cross(Dcb, Ddc);
  const double sign = double3::dot(Dcb, double3::cross(Pb, Pc));
  coordinate.phi = std::copysign(std::acos(cosPhi), sign);
  const double sinPhiRaw = std::sin(coordinate.phi);
  coordinate.sinPhi = std::copysign(std::max(1.0e-8, std::fabs(sinPhiRaw)), sinPhiRaw);

  const double d = dotAB / rbc;
  const double e = dotCD / rbc;

  const double3 dtA = (ds - cosPhi * dr) / r;
  const double3 dtD = (dr - cosPhi * ds) / s;
  const double3 dtB = dtA * (d - 1.0) + e * dtD;
  const double3 dtC = -dtD * (e + 1.0) - d * dtA;

  coordinate.gradient[0] = dtA;
  coordinate.gradient[1] = dtB;
  coordinate.gradient[2] = dtC;
  coordinate.gradient[3] = dtD;

  // Second derivatives of cos(phi); RASPA2 torsion Hessian expressions with DF = 1.
  const double3 DIL = Dcb / rbc;
  const double3 DDJ = ((2.0 * d - 1.0) * Dcb - Dab / rbc) / rbc;
  const double3 DDK = -(2.0 * d * Dcb - Dab / rbc) / rbc;
  const double3 DEJ = (2.0 * e * Dcb - Ddc / rbc) / rbc;
  const double3 DEK = -((2.0 * e + 1.0) * Dcb - Ddc / rbc) / rbc;

  CoordinateBlock AA{};
  CoordinateBlock DD{};
  CoordinateBlock AD{};
  for (std::size_t a = 0; a < 3; ++a)
  {
    for (std::size_t b = 0; b < 3; ++b)
    {
      const double delta = (a == b ? 1.0 : 0.0);
      AA[a][b] = (cosPhi * (component(dr, a) * component(dr, b) + component(Dcb, a) * component(Dcb, b) - delta) / r -
                  component(dr, a) * component(dtA, b) - component(dr, b) * component(dtA, a)) /
                 r;
      DD[a][b] = (cosPhi * (component(ds, a) * component(ds, b) + component(Dcb, a) * component(Dcb, b) - delta) / s -
                  component(ds, a) * component(dtD, b) - component(ds, b) * component(dtD, a)) /
                 s;
      AD[a][b] = (cosPhi * component(dr, a) * component(ds, b) - component(dr, a) * component(dr, b) -
                  component(ds, a) * component(ds, b) - component(Dcb, a) * component(Dcb, b) + delta) /
                 (r * s);
    }
  }

  CoordinateBlock AB{};
  CoordinateBlock AC{};
  CoordinateBlock BD{};
  CoordinateBlock CD{};
  for (std::size_t a = 0; a < 3; ++a)
  {
    for (std::size_t b = 0; b < 3; ++b)
    {
      AB[a][b] = AA[a][b] * (d - 1.0) + AD[a][b] * e + component(DIL, a) * component(dtA, b);
      AC[a][b] = -AA[a][b] * d - AD[a][b] * (e + 1.0) - component(DIL, a) * component(dtA, b);
      BD[a][b] = AD[a][b] * (d - 1.0) + DD[a][b] * e + component(DIL, b) * component(dtD, a);
      CD[a][b] = -AD[a][b] * d - DD[a][b] * (e + 1.0) - component(DIL, b) * component(dtD, a);
    }
  }

  CoordinateBlock BB{};
  CoordinateBlock CC{};
  CoordinateBlock BC{};
  for (std::size_t a = 0; a < 3; ++a)
  {
    for (std::size_t b = 0; b < 3; ++b)
    {
      if (b >= a)
      {
        BB[a][b] = AB[a][b] * (d - 1.0) + BD[b][a] * e + component(DDJ, b) * component(dtA, a) +
                   component(DEJ, b) * component(dtD, a);
        CC[a][b] = -CD[a][b] * (e + 1.0) - AC[b][a] * d - component(DDK, a) * component(dtA, b) -
                   component(DEK, a) * component(dtD, b);
      }
      BC[a][b] = -AB[b][a] * d - BD[a][b] * (e + 1.0) - component(DDJ, a) * component(dtA, b) -
                 component(DEJ, a) * component(dtD, b);
    }
  }
  for (std::size_t a = 0; a < 3; ++a)
  {
    for (std::size_t b = 0; b < a; ++b)
    {
      BB[a][b] = BB[b][a];
      CC[a][b] = CC[b][a];
    }
  }

  coordinate.hessian[0][0] = AA;
  coordinate.hessian[1][1] = BB;
  coordinate.hessian[2][2] = CC;
  coordinate.hessian[3][3] = DD;
  coordinate.hessian[0][1] = AB;
  coordinate.hessian[1][0] = transposed(AB);
  coordinate.hessian[0][2] = AC;
  coordinate.hessian[2][0] = transposed(AC);
  coordinate.hessian[0][3] = AD;
  coordinate.hessian[3][0] = transposed(AD);
  coordinate.hessian[1][2] = BC;
  coordinate.hessian[2][1] = transposed(BC);
  coordinate.hessian[1][3] = BD;
  coordinate.hessian[3][1] = transposed(BD);
  coordinate.hessian[2][3] = CD;
  coordinate.hessian[3][2] = transposed(CD);

  return coordinate;
}

Minimization::InternalCoordinate Minimization::dihedralAngleInternalCoordinate(const double3 &posA, const double3 &posB,
                                                                               const double3 &posC, const double3 &posD)
{
  // Reuse the cos(phi) coordinate and apply the phi = f(cos(phi)) chain rule.
  InternalCoordinate cosine = dihedralInternalCoordinate(posA, posB, posC, posD);

  InternalCoordinate coordinate{};
  coordinate.numAtoms = cosine.numAtoms;
  coordinate.value = cosine.phi;
  coordinate.phi = cosine.phi;
  coordinate.sinPhi = cosine.sinPhi;

  const double firstFactor = -1.0 / cosine.sinPhi;
  const double secondFactor = -cosine.value / (cosine.sinPhi * cosine.sinPhi * cosine.sinPhi);
  for (std::size_t i = 0; i < 4; ++i)
  {
    coordinate.gradient[i] = firstFactor * cosine.gradient[i];
    for (std::size_t j = 0; j < 4; ++j)
    {
      for (std::size_t a = 0; a < 3; ++a)
      {
        for (std::size_t b = 0; b < 3; ++b)
        {
          coordinate.hessian[i][j][a][b] = firstFactor * cosine.hessian[i][j][a][b] +
                                           secondFactor * component(cosine.gradient[i], a) *
                                               component(cosine.gradient[j], b);
        }
      }
    }
  }

  return coordinate;
}
