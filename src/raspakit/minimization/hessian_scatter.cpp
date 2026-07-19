module;

module minimization_hessian_scatter;

import std;

import double3x3;
import minimization_rigid_kinematics;

namespace
{
constexpr std::array<double, 3> drComponents(const double3& dr) { return {dr.x, dr.y, dr.z}; }

std::optional<std::size_t> positionDofBase(const MinimizationDofLayout& layout, std::size_t moleculeIndex,
                                           std::size_t localAtomIndex, bool rigid)
{
  if (rigid)
  {
    return layout.atomRigidComDof(moleculeIndex, localAtomIndex);
  }
  return layout.flexibleAtomDof(moleculeIndex, localAtomIndex, MinimizationDofAxis::X);
}

std::optional<std::size_t> orientationDofBase(const MinimizationDofLayout& layout, std::size_t moleculeIndex,
                                              std::size_t localAtomIndex)
{
  return layout.atomRigidOrientationDof(moleculeIndex, localAtomIndex);
}

void addMatrixBlock(GeneralizedHessian& hessian, std::optional<std::size_t> rowBase,
                    std::optional<std::size_t> columnBase, const double3x3& matrix, double scale = 1.0)
{
  if (!rowBase || !columnBase)
  {
    return;
  }

  const std::array<std::array<double, 3>, 3> values = {
      {{matrix.ax, matrix.ay, matrix.az}, {matrix.bx, matrix.by, matrix.bz}, {matrix.cx, matrix.cy, matrix.cz}}};
  for (std::size_t row = 0; row < 3; ++row)
  {
    for (std::size_t column = 0; column < 3; ++column)
    {
      hessian.add(*rowBase + row, *columnBase + column, scale * values[row][column]);
    }
  }
}

void scatterSameEntityBlock(GeneralizedHessian& hessian, std::size_t baseDof, double f1, double f2, const double3& dr)
{
  const auto components = drComponents(dr);
  for (std::size_t alpha = 0; alpha < 3; ++alpha)
  {
    for (std::size_t beta = 0; beta < 3; ++beta)
    {
      const double value = f2 * components[alpha] * components[beta] + (alpha == beta ? f1 : 0.0);
      hessian.add(baseDof + alpha, baseDof + beta, value);
    }
  }
}

void scatterCrossEntityBlock(GeneralizedHessian& hessian, std::size_t baseDofA, std::size_t baseDofB, double f1,
                             double f2, const double3& dr)
{
  const auto components = drComponents(dr);
  for (std::size_t alpha = 0; alpha < 3; ++alpha)
  {
    for (std::size_t beta = 0; beta < 3; ++beta)
    {
      const double value = f2 * components[alpha] * components[beta] + (alpha == beta ? f1 : 0.0);
      hessian.add(baseDofA + alpha, baseDofB + beta, -value);
      hessian.add(baseDofB + beta, baseDofA + alpha, -value);
    }
  }
}

double3x3 transposeBlock(const double3x3& matrix)
{
  double3x3 transposed{};
  transposed.ax = matrix.ax;
  transposed.ay = matrix.bx;
  transposed.az = matrix.cx;
  transposed.bx = matrix.ay;
  transposed.by = matrix.by;
  transposed.bz = matrix.cy;
  transposed.cx = matrix.az;
  transposed.cy = matrix.bz;
  transposed.cz = matrix.cz;
  return transposed;
}

/*
 * Position-orientation coupling block (Eq. S56 of Dubbeldam, Krishna, Snurr 2009):
 * entry (row = position axis alpha, column = orientation axis beta) is
 *   f2 (dr . DVec_beta) dr_alpha + f1 [DVec_beta]_alpha.
 */
double3x3 positionOrientationBlock(double f1, double f2, const double3& dr,
                                   const Minimization::RigidAtomDerivatives& derivatives)
{
  const std::array<double3, 3> dVec = {derivatives.dVecX, derivatives.dVecY, derivatives.dVecZ};

  double3x3 block{};
  std::array<double*, 9> entries = {&block.ax, &block.ay, &block.az, &block.bx, &block.by,
                                    &block.bz, &block.cx, &block.cy, &block.cz};
  const std::array<double, 3> drComponentsArray = drComponents(dr);
  for (std::size_t alpha = 0; alpha < 3; ++alpha)
  {
    for (std::size_t beta = 0; beta < 3; ++beta)
    {
      const double projection = double3::dot(dr, dVec[beta]);
      const double dVecBetaAlpha = (&dVec[beta].x)[alpha];
      *entries[alpha * 3 + beta] = f2 * projection * drComponentsArray[alpha] + f1 * dVecBetaAlpha;
    }
  }
  return block;
}

void scatterCenterOfMassOrientation(GeneralizedHessian& hessian, std::optional<std::size_t> positionBaseI,
                                    std::optional<std::size_t> orientationBaseI,
                                    std::optional<std::size_t> positionBaseJ,
                                    std::optional<std::size_t> orientationBaseJ,
                                    const Minimization::RigidAtomDerivatives* derivativesI,
                                    const Minimization::RigidAtomDerivatives* derivativesJ, double f1, double f2,
                                    const double3& dr)
{
  // dr = posA - posB. Same-molecule blocks enter with +, cross-molecule blocks with -;
  // both orders are written since the generalized Hessian stores the full dense matrix.
  if (derivativesI != nullptr)
  {
    const double3x3 blockI = positionOrientationBlock(f1, f2, dr, *derivativesI);
    const double3x3 blockITransposed = transposeBlock(blockI);
    addMatrixBlock(hessian, positionBaseI, orientationBaseI, blockI);
    addMatrixBlock(hessian, orientationBaseI, positionBaseI, blockITransposed);
    addMatrixBlock(hessian, orientationBaseI, positionBaseJ, blockITransposed, -1.0);
    addMatrixBlock(hessian, positionBaseJ, orientationBaseI, blockI, -1.0);
  }

  if (derivativesJ != nullptr)
  {
    const double3x3 blockJ = positionOrientationBlock(f1, f2, dr, *derivativesJ);
    const double3x3 blockJTransposed = transposeBlock(blockJ);
    addMatrixBlock(hessian, positionBaseJ, orientationBaseJ, blockJ);
    addMatrixBlock(hessian, orientationBaseJ, positionBaseJ, blockJTransposed);
    addMatrixBlock(hessian, positionBaseI, orientationBaseJ, blockJ, -1.0);
    addMatrixBlock(hessian, orientationBaseJ, positionBaseI, blockJTransposed, -1.0);
  }
}

/*
 * Orientation-orientation blocks (Eq. S55 of Dubbeldam, Krishna, Snurr 2009).
 * Same molecule: f1 (dr . DDVec_ab) sign(dr wrt this molecule) +
 *                f2 (dr . DVec_a)(dr . DVec_b) + f1 (DVec_a . DVec_b).
 * Cross: -[f2 (dr . DVecI_a)(dr . DVecJ_b) + f1 (DVecI_a . DVecJ_b)].
 */
void scatterOrientationOrientation(GeneralizedHessian& hessian, std::optional<std::size_t> orientationBaseI,
                                   std::optional<std::size_t> orientationBaseJ,
                                   const Minimization::RigidAtomDerivatives* derivativesI,
                                   const Minimization::RigidAtomDerivatives* derivativesJ, double f1, double f2,
                                   const double3& dr)
{
  auto dVecs = [](const Minimization::RigidAtomDerivatives& derivatives)
  { return std::array<double3, 3>{derivatives.dVecX, derivatives.dVecY, derivatives.dVecZ}; };

  // ddVec[alpha][beta], symmetric in (alpha, beta)
  auto ddVecs = [](const Minimization::RigidAtomDerivatives& derivatives)
  {
    return std::array<std::array<double3, 3>, 3>{{{derivatives.ddVecAX, derivatives.ddVecAY, derivatives.ddVecAZ},
                                                  {derivatives.ddVecAY, derivatives.ddVecBY, derivatives.ddVecBZ},
                                                  {derivatives.ddVecAZ, derivatives.ddVecBZ, derivatives.ddVecCZ}}};
  };

  // dr = posA - posB, so d(dr)/domegaI = +DVecI and d(dr)/domegaJ = -DVecJ;
  // the DDVec term enters with + for molecule I and - for molecule J.
  auto sameMoleculeBlock = [&](const Minimization::RigidAtomDerivatives& derivatives, double ddSign)
  {
    const auto dVec = dVecs(derivatives);
    const auto ddVec = ddVecs(derivatives);
    double3x3 block{};
    std::array<double*, 9> entries = {&block.ax, &block.ay, &block.az, &block.bx, &block.by,
                                      &block.bz, &block.cx, &block.cy, &block.cz};
    for (std::size_t alpha = 0; alpha < 3; ++alpha)
    {
      for (std::size_t beta = 0; beta < 3; ++beta)
      {
        *entries[alpha * 3 + beta] = ddSign * f1 * double3::dot(dr, ddVec[alpha][beta]) +
                                     f2 * double3::dot(dr, dVec[alpha]) * double3::dot(dr, dVec[beta]) +
                                     f1 * double3::dot(dVec[alpha], dVec[beta]);
      }
    }
    return block;
  };

  if (derivativesI != nullptr && orientationBaseI)
  {
    addMatrixBlock(hessian, orientationBaseI, orientationBaseI, sameMoleculeBlock(*derivativesI, 1.0));
  }

  if (derivativesJ != nullptr && orientationBaseJ)
  {
    addMatrixBlock(hessian, orientationBaseJ, orientationBaseJ, sameMoleculeBlock(*derivativesJ, -1.0));
  }

  if (derivativesI != nullptr && derivativesJ != nullptr && orientationBaseI && orientationBaseJ)
  {
    const auto dVecI = dVecs(*derivativesI);
    const auto dVecJ = dVecs(*derivativesJ);
    double3x3 cross{};
    std::array<double*, 9> entries = {&cross.ax, &cross.ay, &cross.az, &cross.bx, &cross.by,
                                      &cross.bz, &cross.cx, &cross.cy, &cross.cz};
    for (std::size_t alpha = 0; alpha < 3; ++alpha)
    {
      for (std::size_t beta = 0; beta < 3; ++beta)
      {
        *entries[alpha * 3 + beta] = -f2 * double3::dot(dr, dVecI[alpha]) * double3::dot(dr, dVecJ[beta]) -
                                     f1 * double3::dot(dVecI[alpha], dVecJ[beta]);
      }
    }
    addMatrixBlock(hessian, orientationBaseI, orientationBaseJ, cross);
    addMatrixBlock(hessian, orientationBaseJ, orientationBaseI, transposeBlock(cross));
  }
}

void addBendCartesianBlock(GeneralizedHessian& hessian, std::optional<std::size_t> baseI,
                           std::optional<std::size_t> baseJ, const double3& dtI, const double3& dtJ, double ddf,
                           const double3x3& d2)
{
  if (!baseI || !baseJ)
  {
    return;
  }

  const std::array<double, 3> dti = {dtI.x, dtI.y, dtI.z};
  const std::array<double, 3> dtj = {dtJ.x, dtJ.y, dtJ.z};
  // RASPA2 letter convention: letter is the row (a=x, b=y, c=z), suffix is the column.
  const std::array<std::array<double, 3>, 3> d2Values = {
      {{d2.ax, d2.ay, d2.az}, {d2.bx, d2.by, d2.bz}, {d2.cx, d2.cy, d2.cz}}};

  for (std::size_t row = 0; row < 3; ++row)
  {
    for (std::size_t column = 0; column < 3; ++column)
    {
      hessian.add(*baseI + row, *baseJ + column, ddf * dti[row] * dtj[column] + d2Values[row][column]);
    }
  }
}

void addBendOffDiagonalBlock(GeneralizedHessian& hessian, std::optional<std::size_t> baseI,
                             std::optional<std::size_t> baseJ, const double3& dtI, const double3& dtJ, double ddf,
                             const double3x3& d2)
{
  addBendCartesianBlock(hessian, baseI, baseJ, dtI, dtJ, ddf, d2);
  addBendCartesianBlock(hessian, baseJ, baseI, dtJ, dtI, ddf, transposeBlock(d2));
}
}  // namespace

void Minimization::scatterAtomicPositionPosition(GeneralizedHessian& hessian, const MinimizationDofLayout& layout,
                                                 std::size_t moleculeA, std::size_t localAtomA, std::size_t moleculeB,
                                                 std::size_t localAtomB, double f1, double f2, const double3& dr)
{
  const auto baseA = layout.flexibleAtomDof(moleculeA, localAtomA, MinimizationDofAxis::X);
  const auto baseB = layout.flexibleAtomDof(moleculeB, localAtomB, MinimizationDofAxis::X);
  if (!baseA || !baseB)
  {
    return;
  }

  scatterAtomicPositionPositionByDof(hessian, *baseA, *baseB, f1, f2, dr);
}

void Minimization::scatterAtomicPositionPositionByDof(GeneralizedHessian& hessian, std::size_t baseA, std::size_t baseB,
                                                      double f1, double f2, const double3& dr)
{
  if (baseA == baseB)
  {
    scatterSameEntityBlock(hessian, baseA, f1, f2, dr);
    return;
  }

  scatterSameEntityBlock(hessian, baseA, f1, f2, dr);
  scatterSameEntityBlock(hessian, baseB, f1, f2, dr);
  scatterCrossEntityBlock(hessian, baseA, baseB, f1, f2, dr);
}

void Minimization::scatterInteractionHessian(GeneralizedHessian& hessian, const MinimizationDofLayout& layout,
                                             const RigidDerivativeCache& rigidCache, std::size_t moleculeA,
                                             std::size_t localAtomA, bool rigidA, const double3& posA,
                                             const double3& comA, std::size_t moleculeB, std::size_t localAtomB,
                                             bool rigidB, const double3& posB, const double3& comB, double f1,
                                             double f2, const double3& dr)
{
  const auto positionBaseA = positionDofBase(layout, moleculeA, localAtomA, rigidA);
  const auto positionBaseB = positionDofBase(layout, moleculeB, localAtomB, rigidB);
  if (!positionBaseA || !positionBaseB)
  {
    return;
  }

  if (moleculeA == moleculeB && localAtomA == localAtomB && rigidA == rigidB)
  {
    scatterSameEntityBlock(hessian, *positionBaseA, f1, f2, dr);
    return;
  }

  scatterSameEntityBlock(hessian, *positionBaseA, f1, f2, dr);
  scatterSameEntityBlock(hessian, *positionBaseB, f1, f2, dr);
  scatterCrossEntityBlock(hessian, *positionBaseA, *positionBaseB, f1, f2, dr);

  const Minimization::RigidAtomDerivatives* derivativesA = rigidA ? &rigidCache.atom(moleculeA, localAtomA) : nullptr;
  const Minimization::RigidAtomDerivatives* derivativesB = rigidB ? &rigidCache.atom(moleculeB, localAtomB) : nullptr;

  scatterCenterOfMassOrientation(hessian, positionBaseA, orientationDofBase(layout, moleculeA, localAtomA),
                                 positionBaseB, orientationDofBase(layout, moleculeB, localAtomB), derivativesA,
                                 derivativesB, f1, f2, dr);
  scatterOrientationOrientation(hessian, orientationDofBase(layout, moleculeA, localAtomA),
                                orientationDofBase(layout, moleculeB, localAtomB), derivativesA, derivativesB, f1, f2,
                                dr);

  (void)posA;
  (void)comA;
  (void)posB;
  (void)comB;
}

void Minimization::scatterFrameworkMoleculeHessian(GeneralizedHessian& hessian, const MinimizationDofLayout& layout,
                                                   const RigidDerivativeCache& rigidCache, std::size_t moleculeIndex,
                                                   std::size_t localAtom, bool rigid, double f1, double f2,
                                                   const double3& dr)
{
  const auto positionBase = positionDofBase(layout, moleculeIndex, localAtom, rigid);
  if (!positionBase)
  {
    return;
  }

  scatterSameEntityBlock(hessian, *positionBase, f1, f2, dr);

  if (rigid)
  {
    const Minimization::RigidAtomDerivatives* derivatives = &rigidCache.atom(moleculeIndex, localAtom);
    const auto orientationBase = orientationDofBase(layout, moleculeIndex, localAtom);
    scatterCenterOfMassOrientation(hessian, positionBase, orientationBase, std::nullopt, std::nullopt, derivatives,
                                   nullptr, f1, f2, dr);
    scatterOrientationOrientation(hessian, orientationBase, std::nullopt, derivatives, nullptr, f1, f2, dr);
  }
}

void Minimization::scatterFlexibleFrameworkMoleculeHessian(GeneralizedHessian& hessian,
                                                           const MinimizationDofLayout& layout,
                                                           const RigidDerivativeCache& rigidCache,
                                                           std::size_t frameworkAtom, std::size_t moleculeIndex,
                                                           std::size_t localAtom, bool moleculeRigid, double f1,
                                                           double f2, const double3& dr)
{
  const auto moleculePositionBase = positionDofBase(layout, moleculeIndex, localAtom, moleculeRigid);
  const auto frameworkPositionBase = layout.frameworkAtomDof(frameworkAtom, MinimizationDofAxis::X);
  if (!moleculePositionBase || !frameworkPositionBase) return;

  scatterSameEntityBlock(hessian, *moleculePositionBase, f1, f2, dr);
  scatterSameEntityBlock(hessian, *frameworkPositionBase, f1, f2, dr);
  scatterCrossEntityBlock(hessian, *moleculePositionBase, *frameworkPositionBase, f1, f2, dr);

  if (moleculeRigid)
  {
    const RigidAtomDerivatives* derivatives = &rigidCache.atom(moleculeIndex, localAtom);
    const auto orientationBase = orientationDofBase(layout, moleculeIndex, localAtom);
    scatterCenterOfMassOrientation(hessian, moleculePositionBase, orientationBase, frameworkPositionBase, std::nullopt,
                                   derivatives, nullptr, f1, f2, dr);
    scatterOrientationOrientation(hessian, orientationBase, std::nullopt, derivatives, nullptr, f1, f2, dr);
  }
}

void Minimization::scatterBendHessian(GeneralizedHessian& hessian, const MinimizationDofLayout& layout,
                                      std::size_t moleculeIndex, std::size_t localAtomA, std::size_t localAtomB,
                                      std::size_t localAtomC, const BendHessianGeometry& geometry)
{
  const auto baseA = layout.flexibleAtomDof(moleculeIndex, localAtomA, MinimizationDofAxis::X);
  const auto baseB = layout.flexibleAtomDof(moleculeIndex, localAtomB, MinimizationDofAxis::X);
  const auto baseC = layout.flexibleAtomDof(moleculeIndex, localAtomC, MinimizationDofAxis::X);
  if (!baseA || !baseB || !baseC)
  {
    return;
  }

  scatterBendHessianByDof(hessian, {*baseA, *baseB, *baseC}, geometry);
}

void Minimization::scatterBendHessianByDof(GeneralizedHessian& hessian, const std::array<std::size_t, 3>& bases,
                                           const BendHessianGeometry& geometry)
{
  const std::optional<std::size_t> baseA = bases[0];
  const std::optional<std::size_t> baseB = bases[1];
  const std::optional<std::size_t> baseC = bases[2];

  const double DF = geometry.DF;
  const double DDF = geometry.DDF;
  const double3& dtA = geometry.dtA;
  const double3& dtB = geometry.dtB;
  const double3& dtC = geometry.dtC;
  const double3& rab = geometry.Rab;
  const double3& rbc = geometry.Rbc;
  const double cosTheta = geometry.cosTheta;
  const double rabLength = geometry.rab;
  const double rbcLength = geometry.rbc;

  double3x3 d2I{};
  d2I.ax = -DF * (2.0 * dtA.x * rab.x + cosTheta * (1.0 - rab.x * rab.x) / rabLength) / rabLength;
  d2I.by = -DF * (2.0 * dtA.y * rab.y + cosTheta * (1.0 - rab.y * rab.y) / rabLength) / rabLength;
  d2I.cz = -DF * (2.0 * dtA.z * rab.z + cosTheta * (1.0 - rab.z * rab.z) / rabLength) / rabLength;
  d2I.ay = DF * (cosTheta * rab.x * rab.y / rabLength - dtA.x * rab.y - dtA.y * rab.x) / rabLength;
  d2I.az = DF * (cosTheta * rab.x * rab.z / rabLength - dtA.x * rab.z - dtA.z * rab.x) / rabLength;
  d2I.bz = DF * (cosTheta * rab.y * rab.z / rabLength - dtA.y * rab.z - dtA.z * rab.y) / rabLength;
  d2I.bx = d2I.ay;
  d2I.cx = d2I.az;
  d2I.cy = d2I.bz;

  double3x3 d2K{};
  d2K.ax = -DF * (2.0 * dtC.x * rbc.x + cosTheta * (1.0 - rbc.x * rbc.x) / rbcLength) / rbcLength;
  d2K.by = -DF * (2.0 * dtC.y * rbc.y + cosTheta * (1.0 - rbc.y * rbc.y) / rbcLength) / rbcLength;
  d2K.cz = -DF * (2.0 * dtC.z * rbc.z + cosTheta * (1.0 - rbc.z * rbc.z) / rbcLength) / rbcLength;
  d2K.ay = DF * (cosTheta * rbc.x * rbc.y / rbcLength - dtC.x * rbc.y - dtC.y * rbc.x) / rbcLength;
  d2K.az = DF * (cosTheta * rbc.x * rbc.z / rbcLength - dtC.x * rbc.z - dtC.z * rbc.x) / rbcLength;
  d2K.bz = DF * (cosTheta * rbc.y * rbc.z / rbcLength - dtC.y * rbc.z - dtC.z * rbc.y) / rbcLength;
  d2K.bx = d2K.ay;
  d2K.cx = d2K.az;
  d2K.cy = d2K.bz;

  double3x3 d2IK{};
  d2IK.ax = DF * (cosTheta * rab.x * rbc.x - rab.x * rab.x - rbc.x * rbc.x + 1.0) / (rabLength * rbcLength);
  d2IK.ay = DF * (cosTheta * rab.x * rbc.y - rab.x * rab.y - rbc.x * rbc.y) / (rabLength * rbcLength);
  d2IK.az = DF * (cosTheta * rab.x * rbc.z - rab.x * rab.z - rbc.x * rbc.z) / (rabLength * rbcLength);
  d2IK.bx = DF * (cosTheta * rab.y * rbc.x - rab.y * rab.x - rbc.y * rbc.x) / (rabLength * rbcLength);
  d2IK.by = DF * (cosTheta * rab.y * rbc.y - rab.y * rab.y - rbc.y * rbc.y + 1.0) / (rabLength * rbcLength);
  d2IK.bz = DF * (cosTheta * rab.y * rbc.z - rab.y * rab.z - rbc.y * rbc.z) / (rabLength * rbcLength);
  d2IK.cx = DF * (cosTheta * rab.z * rbc.x - rab.z * rab.x - rbc.z * rbc.x) / (rabLength * rbcLength);
  d2IK.cy = DF * (cosTheta * rab.z * rbc.y - rab.z * rab.y - rbc.z * rbc.y) / (rabLength * rbcLength);
  d2IK.cz = DF * (cosTheta * rab.z * rbc.z - rab.z * rab.z - rbc.z * rbc.z + 1.0) / (rabLength * rbcLength);

  addBendCartesianBlock(hessian, baseA, baseA, dtA, dtA, DDF, d2I);

  double3x3 d2B{};
  d2B.ax = d2I.ax + d2K.ax + 2.0 * d2IK.ax;
  d2B.ay = d2I.ay + d2K.ay + d2IK.ay + d2IK.bx;
  d2B.by = d2I.by + d2K.by + 2.0 * d2IK.by;
  d2B.az = d2I.az + d2K.az + d2IK.az + d2IK.cx;
  d2B.bz = d2I.bz + d2K.bz + d2IK.bz + d2IK.cy;
  d2B.cz = d2I.cz + d2K.cz + 2.0 * d2IK.cz;
  d2B.bx = d2B.ay;
  d2B.cx = d2B.az;
  d2B.cy = d2B.bz;
  addBendCartesianBlock(hessian, baseB, baseB, dtB, dtB, DDF, d2B);
  addBendCartesianBlock(hessian, baseC, baseC, dtC, dtC, DDF, d2K);

  double3x3 d2AB{};
  d2AB.ax = -d2I.ax - d2IK.ax;
  d2AB.ay = -d2I.ay - d2IK.ay;
  d2AB.az = -d2I.az - d2IK.az;
  d2AB.bx = -d2I.ay - d2IK.bx;
  d2AB.by = -d2I.by - d2IK.by;
  d2AB.bz = -d2I.bz - d2IK.bz;
  d2AB.cx = -d2I.az - d2IK.cx;
  d2AB.cy = -d2I.bz - d2IK.cy;
  d2AB.cz = -d2I.cz - d2IK.cz;
  addBendOffDiagonalBlock(hessian, baseA, baseB, dtA, dtB, DDF, d2AB);

  addBendOffDiagonalBlock(hessian, baseA, baseC, dtA, dtC, DDF, d2IK);

  // CB block: -(D2K + D2IK^T); D2IK rows refer to atom A, so it enters transposed here.
  double3x3 d2CB{};
  d2CB.ax = -d2K.ax - d2IK.ax;
  d2CB.ay = -d2K.ay - d2IK.bx;
  d2CB.az = -d2K.az - d2IK.cx;
  d2CB.bx = -d2K.bx - d2IK.ay;
  d2CB.by = -d2K.by - d2IK.by;
  d2CB.bz = -d2K.bz - d2IK.cy;
  d2CB.cx = -d2K.cx - d2IK.az;
  d2CB.cy = -d2K.cy - d2IK.bz;
  d2CB.cz = -d2K.cz - d2IK.cz;
  addBendOffDiagonalBlock(hessian, baseC, baseB, dtC, dtB, DDF, d2CB);
}

void Minimization::scatterTorsionHessian(GeneralizedHessian& hessian, const MinimizationDofLayout& layout,
                                         std::size_t moleculeIndex, std::size_t localAtomA, std::size_t localAtomB,
                                         std::size_t localAtomC, std::size_t localAtomD,
                                         const TorsionHessianGeometry& geometry)
{
  const auto baseA = layout.flexibleAtomDof(moleculeIndex, localAtomA, MinimizationDofAxis::X);
  const auto baseB = layout.flexibleAtomDof(moleculeIndex, localAtomB, MinimizationDofAxis::X);
  const auto baseC = layout.flexibleAtomDof(moleculeIndex, localAtomC, MinimizationDofAxis::X);
  const auto baseD = layout.flexibleAtomDof(moleculeIndex, localAtomD, MinimizationDofAxis::X);
  if (!baseA || !baseB || !baseC || !baseD)
  {
    return;
  }

  scatterTorsionHessianByDof(hessian, {*baseA, *baseB, *baseC, *baseD}, geometry);
}

void Minimization::scatterTorsionHessianByDof(GeneralizedHessian& hessian, const std::array<std::size_t, 4>& bases,
                                              const TorsionHessianGeometry& geometry)
{
  const std::optional<std::size_t> baseA = bases[0];
  const std::optional<std::size_t> baseB = bases[1];
  const std::optional<std::size_t> baseC = bases[2];
  const std::optional<std::size_t> baseD = bases[3];

  const double DF = geometry.DF;
  const double DDF = geometry.DDF;
  const double3& dtA = geometry.dtA;
  const double3& dtB = geometry.dtB;
  const double3& dtC = geometry.dtC;
  const double3& dtD = geometry.dtD;
  const double3& Dab = geometry.Dab;
  const double3& Dcb = geometry.Dcb;  // unit vector
  const double3& Ddc = geometry.Ddc;
  const double rbc = geometry.rbc;
  const double3& dr = geometry.dr;  // unit vector
  const double3& ds = geometry.ds;  // unit vector
  const double r = geometry.r;
  const double s = geometry.s;
  const double d = geometry.d;
  const double e = geometry.e;
  const double cosPhi = geometry.cosPhi;

  // Derivatives of d = (Dab.Dcb)/rbc and e = (Ddc.Dcb)/rbc scaled by DF (RASPA2 internal_hessian.c).
  const double3 DIL = DF * Dcb / rbc;
  const double3 DDJ = DF * ((2.0 * d - 1.0) * Dcb - Dab / rbc) / rbc;
  const double3 DDK = -DF * (2.0 * d * Dcb - Dab / rbc) / rbc;
  const double3 DEJ = DF * (2.0 * e * Dcb - Ddc / rbc) / rbc;
  const double3 DEK = -DF * ((2.0 * e + 1.0) * Dcb - Ddc / rbc) / rbc;

  // Diagonal block for atom A (symmetric; RASPA2 computes the upper triangle only).
  double3x3 d2I{};
  d2I.ax = DF * (cosPhi * (dr.x * dr.x + Dcb.x * Dcb.x - 1.0) / r - 2.0 * dr.x * dtA.x) / r;
  d2I.by = DF * (cosPhi * (dr.y * dr.y + Dcb.y * Dcb.y - 1.0) / r - 2.0 * dr.y * dtA.y) / r;
  d2I.cz = DF * (cosPhi * (dr.z * dr.z + Dcb.z * Dcb.z - 1.0) / r - 2.0 * dr.z * dtA.z) / r;
  d2I.ay = DF * (cosPhi * (dr.x * dr.y + Dcb.x * Dcb.y) / r - dr.x * dtA.y - dr.y * dtA.x) / r;
  d2I.az = DF * (cosPhi * (dr.x * dr.z + Dcb.x * Dcb.z) / r - dr.x * dtA.z - dr.z * dtA.x) / r;
  d2I.bz = DF * (cosPhi * (dr.y * dr.z + Dcb.y * Dcb.z) / r - dr.y * dtA.z - dr.z * dtA.y) / r;
  d2I.bx = d2I.ay;
  d2I.cx = d2I.az;
  d2I.cy = d2I.bz;

  // Diagonal block for atom D (symmetric).
  double3x3 d2L{};
  d2L.ax = DF * (cosPhi * (ds.x * ds.x + Dcb.x * Dcb.x - 1.0) / s - 2.0 * ds.x * dtD.x) / s;
  d2L.by = DF * (cosPhi * (ds.y * ds.y + Dcb.y * Dcb.y - 1.0) / s - 2.0 * ds.y * dtD.y) / s;
  d2L.cz = DF * (cosPhi * (ds.z * ds.z + Dcb.z * Dcb.z - 1.0) / s - 2.0 * ds.z * dtD.z) / s;
  d2L.ay = DF * (cosPhi * (ds.x * ds.y + Dcb.x * Dcb.y) / s - ds.x * dtD.y - ds.y * dtD.x) / s;
  d2L.az = DF * (cosPhi * (ds.x * ds.z + Dcb.x * Dcb.z) / s - ds.x * dtD.z - ds.z * dtD.x) / s;
  d2L.bz = DF * (cosPhi * (ds.y * ds.z + Dcb.y * Dcb.z) / s - ds.y * dtD.z - ds.z * dtD.y) / s;
  d2L.bx = d2L.ay;
  d2L.cx = d2L.az;
  d2L.cy = d2L.bz;

  // AD off-diagonal block (full 3x3; letter is the row, suffix is the column).
  double3x3 d2IL{};
  d2IL.ax = DF * (cosPhi * dr.x * ds.x - dr.x * dr.x - ds.x * ds.x - Dcb.x * Dcb.x + 1.0) / (r * s);
  d2IL.ay = DF * (cosPhi * dr.x * ds.y - dr.x * dr.y - ds.x * ds.y - Dcb.x * Dcb.y) / (r * s);
  d2IL.az = DF * (cosPhi * dr.x * ds.z - dr.x * dr.z - ds.x * ds.z - Dcb.x * Dcb.z) / (r * s);
  d2IL.bx = DF * (cosPhi * dr.y * ds.x - dr.y * dr.x - ds.y * ds.x - Dcb.y * Dcb.x) / (r * s);
  d2IL.by = DF * (cosPhi * dr.y * ds.y - dr.y * dr.y - ds.y * ds.y - Dcb.y * Dcb.y + 1.0) / (r * s);
  d2IL.bz = DF * (cosPhi * dr.y * ds.z - dr.y * dr.z - ds.y * ds.z - Dcb.y * Dcb.z) / (r * s);
  d2IL.cx = DF * (cosPhi * dr.z * ds.x - dr.z * dr.x - ds.z * ds.x - Dcb.z * Dcb.x) / (r * s);
  d2IL.cy = DF * (cosPhi * dr.z * ds.y - dr.z * dr.y - ds.z * ds.y - Dcb.z * Dcb.y) / (r * s);
  d2IL.cz = DF * (cosPhi * dr.z * ds.z - dr.z * dr.z - ds.z * ds.z - Dcb.z * Dcb.z + 1.0) / (r * s);

  // AB off-diagonal block.
  double3x3 d2IJ{};
  d2IJ.ax = d2I.ax * (d - 1.0) + d2IL.ax * e + DIL.x * dtA.x;
  d2IJ.ay = d2I.ay * (d - 1.0) + d2IL.ay * e + DIL.x * dtA.y;
  d2IJ.az = d2I.az * (d - 1.0) + d2IL.az * e + DIL.x * dtA.z;
  d2IJ.bx = d2I.ay * (d - 1.0) + d2IL.bx * e + DIL.y * dtA.x;
  d2IJ.by = d2I.by * (d - 1.0) + d2IL.by * e + DIL.y * dtA.y;
  d2IJ.bz = d2I.bz * (d - 1.0) + d2IL.bz * e + DIL.y * dtA.z;
  d2IJ.cx = d2I.az * (d - 1.0) + d2IL.cx * e + DIL.z * dtA.x;
  d2IJ.cy = d2I.bz * (d - 1.0) + d2IL.cy * e + DIL.z * dtA.y;
  d2IJ.cz = d2I.cz * (d - 1.0) + d2IL.cz * e + DIL.z * dtA.z;

  // AC off-diagonal block.
  double3x3 d2IK{};
  d2IK.ax = -d2I.ax * d - d2IL.ax * (e + 1.0) - DIL.x * dtA.x;
  d2IK.ay = -d2I.ay * d - d2IL.ay * (e + 1.0) - DIL.x * dtA.y;
  d2IK.az = -d2I.az * d - d2IL.az * (e + 1.0) - DIL.x * dtA.z;
  d2IK.bx = -d2I.ay * d - d2IL.bx * (e + 1.0) - DIL.y * dtA.x;
  d2IK.by = -d2I.by * d - d2IL.by * (e + 1.0) - DIL.y * dtA.y;
  d2IK.bz = -d2I.bz * d - d2IL.bz * (e + 1.0) - DIL.y * dtA.z;
  d2IK.cx = -d2I.az * d - d2IL.cx * (e + 1.0) - DIL.z * dtA.x;
  d2IK.cy = -d2I.bz * d - d2IL.cy * (e + 1.0) - DIL.z * dtA.y;
  d2IK.cz = -d2I.cz * d - d2IL.cz * (e + 1.0) - DIL.z * dtA.z;

  // BD off-diagonal block.
  double3x3 d2JL{};
  d2JL.ax = d2IL.ax * (d - 1.0) + d2L.ax * e + DIL.x * dtD.x;
  d2JL.ay = d2IL.ay * (d - 1.0) + d2L.ay * e + DIL.y * dtD.x;
  d2JL.az = d2IL.az * (d - 1.0) + d2L.az * e + DIL.z * dtD.x;
  d2JL.bx = d2IL.bx * (d - 1.0) + d2L.ay * e + DIL.x * dtD.y;
  d2JL.by = d2IL.by * (d - 1.0) + d2L.by * e + DIL.y * dtD.y;
  d2JL.bz = d2IL.bz * (d - 1.0) + d2L.bz * e + DIL.z * dtD.y;
  d2JL.cx = d2IL.cx * (d - 1.0) + d2L.az * e + DIL.x * dtD.z;
  d2JL.cy = d2IL.cy * (d - 1.0) + d2L.bz * e + DIL.y * dtD.z;
  d2JL.cz = d2IL.cz * (d - 1.0) + d2L.cz * e + DIL.z * dtD.z;

  // CD off-diagonal block.
  double3x3 d2KL{};
  d2KL.ax = -d2IL.ax * d - d2L.ax * (e + 1.0) - DIL.x * dtD.x;
  d2KL.ay = -d2IL.ay * d - d2L.ay * (e + 1.0) - DIL.y * dtD.x;
  d2KL.az = -d2IL.az * d - d2L.az * (e + 1.0) - DIL.z * dtD.x;
  d2KL.bx = -d2IL.bx * d - d2L.ay * (e + 1.0) - DIL.x * dtD.y;
  d2KL.by = -d2IL.by * d - d2L.by * (e + 1.0) - DIL.y * dtD.y;
  d2KL.bz = -d2IL.bz * d - d2L.bz * (e + 1.0) - DIL.z * dtD.y;
  d2KL.cx = -d2IL.cx * d - d2L.az * (e + 1.0) - DIL.x * dtD.z;
  d2KL.cy = -d2IL.cy * d - d2L.bz * (e + 1.0) - DIL.y * dtD.z;
  d2KL.cz = -d2IL.cz * d - d2L.cz * (e + 1.0) - DIL.z * dtD.z;

  // Diagonal block for atom B (symmetric).
  double3x3 d2J{};
  d2J.ax = d2IJ.ax * (d - 1.0) + d2JL.ax * e + DDJ.x * dtA.x + DEJ.x * dtD.x;
  d2J.by = d2IJ.by * (d - 1.0) + d2JL.by * e + DDJ.y * dtA.y + DEJ.y * dtD.y;
  d2J.cz = d2IJ.cz * (d - 1.0) + d2JL.cz * e + DDJ.z * dtA.z + DEJ.z * dtD.z;
  d2J.ay = d2IJ.ay * (d - 1.0) + d2JL.bx * e + DDJ.y * dtA.x + DEJ.y * dtD.x;
  d2J.az = d2IJ.az * (d - 1.0) + d2JL.cx * e + DDJ.z * dtA.x + DEJ.z * dtD.x;
  d2J.bz = d2IJ.bz * (d - 1.0) + d2JL.cy * e + DDJ.z * dtA.y + DEJ.z * dtD.y;
  d2J.bx = d2J.ay;
  d2J.cx = d2J.az;
  d2J.cy = d2J.bz;

  // Diagonal block for atom C (symmetric).
  double3x3 d2K{};
  d2K.ax = -d2KL.ax * (e + 1.0) - d2IK.ax * d - DDK.x * dtA.x - DEK.x * dtD.x;
  d2K.by = -d2KL.by * (e + 1.0) - d2IK.by * d - DDK.y * dtA.y - DEK.y * dtD.y;
  d2K.cz = -d2KL.cz * (e + 1.0) - d2IK.cz * d - DDK.z * dtA.z - DEK.z * dtD.z;
  d2K.ay = -d2KL.ay * (e + 1.0) - d2IK.bx * d - DDK.x * dtA.y - DEK.x * dtD.y;
  d2K.az = -d2KL.az * (e + 1.0) - d2IK.cx * d - DDK.x * dtA.z - DEK.x * dtD.z;
  d2K.bz = -d2KL.bz * (e + 1.0) - d2IK.cy * d - DDK.y * dtA.z - DEK.y * dtD.z;
  d2K.bx = d2K.ay;
  d2K.cx = d2K.az;
  d2K.cy = d2K.bz;

  // BC off-diagonal block.
  double3x3 d2JK{};
  d2JK.ax = -d2IJ.ax * d - d2JL.ax * (e + 1.0) - DDJ.x * dtA.x - DEJ.x * dtD.x;
  d2JK.ay = -d2IJ.bx * d - d2JL.ay * (e + 1.0) - DDJ.x * dtA.y - DEJ.x * dtD.y;
  d2JK.az = -d2IJ.cx * d - d2JL.az * (e + 1.0) - DDJ.x * dtA.z - DEJ.x * dtD.z;
  d2JK.bx = -d2IJ.ay * d - d2JL.bx * (e + 1.0) - DDJ.y * dtA.x - DEJ.y * dtD.x;
  d2JK.by = -d2IJ.by * d - d2JL.by * (e + 1.0) - DDJ.y * dtA.y - DEJ.y * dtD.y;
  d2JK.bz = -d2IJ.cy * d - d2JL.bz * (e + 1.0) - DDJ.y * dtA.z - DEJ.y * dtD.z;
  d2JK.cx = -d2IJ.az * d - d2JL.cx * (e + 1.0) - DDJ.z * dtA.x - DEJ.z * dtD.x;
  d2JK.cy = -d2IJ.bz * d - d2JL.cy * (e + 1.0) - DDJ.z * dtA.y - DEJ.z * dtD.y;
  d2JK.cz = -d2IJ.cz * d - d2JL.cz * (e + 1.0) - DDJ.z * dtA.z - DEJ.z * dtD.z;

  addBendCartesianBlock(hessian, baseA, baseA, dtA, dtA, DDF, d2I);
  addBendCartesianBlock(hessian, baseB, baseB, dtB, dtB, DDF, d2J);
  addBendCartesianBlock(hessian, baseC, baseC, dtC, dtC, DDF, d2K);
  addBendCartesianBlock(hessian, baseD, baseD, dtD, dtD, DDF, d2L);

  addBendOffDiagonalBlock(hessian, baseA, baseB, dtA, dtB, DDF, d2IJ);
  addBendOffDiagonalBlock(hessian, baseA, baseC, dtA, dtC, DDF, d2IK);
  addBendOffDiagonalBlock(hessian, baseA, baseD, dtA, dtD, DDF, d2IL);
  addBendOffDiagonalBlock(hessian, baseB, baseC, dtB, dtC, DDF, d2JK);
  addBendOffDiagonalBlock(hessian, baseB, baseD, dtB, dtD, DDF, d2JL);
  addBendOffDiagonalBlock(hessian, baseC, baseD, dtC, dtD, DDF, d2KL);
}

namespace
{
// (generalized DOF index, dp/dq column) contributions mapping a Cartesian site displacement of one
// atom onto its generalized degrees of freedom. Flexible atoms carry three Cartesian columns; atoms
// driven by a rigid body carry three center-of-mass columns plus three orientation columns (dVec).
std::vector<std::pair<std::size_t, double3>> siteEntries(const Minimization::HessianSite& site)
{
  std::vector<std::pair<std::size_t, double3>> entries;
  if (!site.positionBase)
  {
    return entries;
  }
  entries.reserve(site.orientationBase ? 6 : 3);
  entries.emplace_back(*site.positionBase + 0, double3(1.0, 0.0, 0.0));
  entries.emplace_back(*site.positionBase + 1, double3(0.0, 1.0, 0.0));
  entries.emplace_back(*site.positionBase + 2, double3(0.0, 0.0, 1.0));
  if (site.orientationBase && site.derivatives != nullptr)
  {
    entries.emplace_back(*site.orientationBase + 0, site.derivatives->dVecX);
    entries.emplace_back(*site.orientationBase + 1, site.derivatives->dVecY);
    entries.emplace_back(*site.orientationBase + 2, site.derivatives->dVecZ);
  }
  return entries;
}

// Project a Cartesian per-term Hessian (a dense 3N x 3N matrix, atom-major) onto the generalized
// DOFs of the participating sites and add the gradient-curvature term on rigid orientation blocks.
template <std::size_t N>
void projectCartesianTerm(GeneralizedHessian& hessian, const std::array<Minimization::HessianSite, N>& sites,
                          const std::array<double3, N>& gradients, const GeneralizedHessian& local)
{
  std::array<std::vector<std::pair<std::size_t, double3>>, N> entries;
  for (std::size_t i = 0; i < N; ++i)
  {
    entries[i] = siteEntries(sites[i]);
  }

  for (std::size_t i = 0; i < N; ++i)
  {
    for (std::size_t j = 0; j < N; ++j)
    {
      for (const auto& [dofI, weightI] : entries[i])
      {
        for (const auto& [dofJ, weightJ] : entries[j])
        {
          double value = 0.0;
          const std::array<double, 3> wi = {weightI.x, weightI.y, weightI.z};
          const std::array<double, 3> wj = {weightJ.x, weightJ.y, weightJ.z};
          for (std::size_t r = 0; r < 3; ++r)
          {
            for (std::size_t c = 0; c < 3; ++c)
            {
              value += wi[r] * local(3 * i + r, 3 * j + c) * wj[c];
            }
          }
          hessian.add(dofI, dofJ, value);
        }
      }
    }
  }

  // Gradient curvature: sum_i g_i . d2p_i/domega_a domega_b on each rigid orientation block.
  for (std::size_t i = 0; i < N; ++i)
  {
    if (!sites[i].orientationBase || sites[i].derivatives == nullptr)
    {
      continue;
    }
    const Minimization::RigidAtomDerivatives& d = *sites[i].derivatives;
    const std::array<std::array<double3, 3>, 3> ddVec = {{{d.ddVecAX, d.ddVecAY, d.ddVecAZ},
                                                          {d.ddVecAY, d.ddVecBY, d.ddVecBZ},
                                                          {d.ddVecAZ, d.ddVecBZ, d.ddVecCZ}}};
    const std::size_t base = *sites[i].orientationBase;
    for (std::size_t a = 0; a < 3; ++a)
    {
      for (std::size_t b = 0; b < 3; ++b)
      {
        hessian.add(base + a, base + b, double3::dot(gradients[i], ddVec[a][b]));
      }
    }
  }
}
}  // namespace

Minimization::HessianSite Minimization::makeHessianSite(const MinimizationDofLayout& layout,
                                                        const RigidDerivativeCache& rigidCache,
                                                        std::size_t moleculeIndex, std::size_t localAtom)
{
  HessianSite site{};
  if (const auto comBase = layout.atomRigidComDof(moleculeIndex, localAtom))
  {
    site.positionBase = comBase;
    site.orientationBase = layout.atomRigidOrientationDof(moleculeIndex, localAtom);
    site.derivatives = &rigidCache.atom(moleculeIndex, localAtom);
  }
  else
  {
    site.positionBase = layout.flexibleAtomDof(moleculeIndex, localAtom, MinimizationDofAxis::X);
  }
  return site;
}

Minimization::HessianSite Minimization::makeFrameworkHessianSite(const MinimizationDofLayout& layout,
                                                                 const RigidDerivativeCache& rigidCache,
                                                                 std::size_t frameworkAtom)
{
  HessianSite site{};
  if (const auto comBase = layout.frameworkAtomRigidComDof(frameworkAtom))
  {
    site.positionBase = comBase;
    site.orientationBase = *comBase + 3;
    site.derivatives = &rigidCache.frameworkAtom(frameworkAtom);
  }
  else
  {
    site.positionBase = layout.frameworkAtomDof(frameworkAtom, MinimizationDofAxis::X);
  }
  return site;
}

void Minimization::scatterRadialHessianSites(GeneralizedHessian& hessian, const std::array<HessianSite, 2>& sites,
                                             const std::array<double3, 2>& gradients, double f1, double f2,
                                             const double3& dr)
{
  GeneralizedHessian local(6, 0);
  scatterAtomicPositionPositionByDof(local, 0, 3, f1, f2, dr);
  projectCartesianTerm<2>(hessian, sites, gradients, local);
}

void Minimization::scatterBendHessianSites(GeneralizedHessian& hessian, const std::array<HessianSite, 3>& sites,
                                           const std::array<double3, 3>& gradients, const BendHessianGeometry& geometry)
{
  GeneralizedHessian local(9, 0);
  scatterBendHessianByDof(local, {0, 3, 6}, geometry);
  projectCartesianTerm<3>(hessian, sites, gradients, local);
}

void Minimization::scatterTorsionHessianSites(GeneralizedHessian& hessian, const std::array<HessianSite, 4>& sites,
                                              const std::array<double3, 4>& gradients,
                                              const TorsionHessianGeometry& geometry)
{
  GeneralizedHessian local(12, 0);
  scatterTorsionHessianByDof(local, {0, 3, 6, 9}, geometry);
  projectCartesianTerm<4>(hessian, sites, gradients, local);
}

void Minimization::scatterCartesianTermSites(GeneralizedHessian& hessian, std::span<const HessianSite> sites,
                                             std::span<const double3> gradients,
                                             const GeneralizedHessian& localCartesian)
{
  const std::size_t n = sites.size();

  std::vector<std::vector<std::pair<std::size_t, double3>>> entries(n);
  for (std::size_t i = 0; i < n; ++i)
  {
    entries[i] = siteEntries(sites[i]);
  }

  for (std::size_t i = 0; i < n; ++i)
  {
    for (std::size_t j = 0; j < n; ++j)
    {
      for (const auto& [dofI, weightI] : entries[i])
      {
        const std::array<double, 3> wi = {weightI.x, weightI.y, weightI.z};
        for (const auto& [dofJ, weightJ] : entries[j])
        {
          const std::array<double, 3> wj = {weightJ.x, weightJ.y, weightJ.z};
          double value = 0.0;
          for (std::size_t r = 0; r < 3; ++r)
          {
            for (std::size_t c = 0; c < 3; ++c)
            {
              value += wi[r] * localCartesian(3 * i + r, 3 * j + c) * wj[c];
            }
          }
          hessian.add(dofI, dofJ, value);
        }
      }
    }
  }

  // Gradient curvature: sum_i g_i . d2p_i/domega_a domega_b on each rigid orientation block.
  for (std::size_t i = 0; i < n; ++i)
  {
    if (!sites[i].orientationBase || sites[i].derivatives == nullptr)
    {
      continue;
    }
    const RigidAtomDerivatives& d = *sites[i].derivatives;
    const std::array<std::array<double3, 3>, 3> ddVec = {{{d.ddVecAX, d.ddVecAY, d.ddVecAZ},
                                                          {d.ddVecAY, d.ddVecBY, d.ddVecBZ},
                                                          {d.ddVecAZ, d.ddVecBZ, d.ddVecCZ}}};
    const std::size_t base = *sites[i].orientationBase;
    for (std::size_t a = 0; a < 3; ++a)
    {
      for (std::size_t b = 0; b < 3; ++b)
      {
        hessian.add(base + a, base + b, double3::dot(gradients[i], ddVec[a][b]));
      }
    }
  }
}

void Minimization::removeRigidOffsetStrainGradientSites(GeneralizedHessian& hessian,
                                                        const MinimizationDofLayout& layout,
                                                        const RigidDerivativeCache& rigidCache,
                                                        std::size_t moleculeIndex,
                                                        std::span<const std::size_t> localAtoms,
                                                        std::span<const double3> positions,
                                                        std::span<const double3> gradients)
{
  for (std::size_t i = 0; i < localAtoms.size(); ++i)
  {
    const std::size_t localAtom = localAtoms[i];
    if (!layout.atomRigidComDof(moleculeIndex, localAtom).has_value())
    {
      continue;
    }
    const double3 offset = positions[i] - rigidCache.bodyCenterOfMass(moleculeIndex, localAtom);
    const double3& gradient = gradients[i];
    hessian.strainGradient().ax -= offset.x * gradient.x;
    hessian.strainGradient().bx -= offset.y * gradient.x;
    hessian.strainGradient().cx -= offset.z * gradient.x;
    hessian.strainGradient().ay -= offset.x * gradient.y;
    hessian.strainGradient().by -= offset.y * gradient.y;
    hessian.strainGradient().cy -= offset.z * gradient.y;
    hessian.strainGradient().az -= offset.x * gradient.z;
    hessian.strainGradient().bz -= offset.y * gradient.z;
    hessian.strainGradient().cz -= offset.z * gradient.z;
  }
}

void Minimization::scatterAtomicPositionStrainIsotropic(GeneralizedHessian& hessian,
                                                        const MinimizationDofLayout& layout, std::size_t moleculeA,
                                                        std::size_t localAtomA, std::size_t moleculeB,
                                                        std::size_t localAtomB, double f1, double f2, const double3& dr)
{
  if (hessian.numStrain() != 1)
  {
    return;
  }

  const auto baseA = layout.flexibleAtomDof(moleculeA, localAtomA, MinimizationDofAxis::X);
  const auto baseB = layout.flexibleAtomDof(moleculeB, localAtomB, MinimizationDofAxis::X);
  if (!baseA || !baseB)
  {
    return;
  }

  const std::array<double, 3> components = drComponents(dr);
  constexpr std::array<MinimizationDofAxis, 3> axes = {MinimizationDofAxis::X, MinimizationDofAxis::Y,
                                                       MinimizationDofAxis::Z};

  const std::array<double, 3> positionStrainA = {
      f2 * components[0] * components[0] * components[0] + 2.0 * f1 * components[0] +
          f2 * components[1] * components[1] * components[0] + f2 * components[2] * components[2] * components[0],
      f2 * components[0] * components[0] * components[1] + f2 * components[1] * components[1] * components[1] +
          2.0 * f1 * components[1] + f2 * components[2] * components[2] * components[1],
      f2 * components[0] * components[0] * components[2] + f2 * components[1] * components[1] * components[2] +
          f2 * components[2] * components[2] * components[2] + 2.0 * f1 * components[2]};

  for (std::size_t axis = 0; axis < 3; ++axis)
  {
    if (auto dofA = layout.flexibleAtomDof(moleculeA, localAtomA, axes[axis]))
    {
      hessian.addPositionStrain(*dofA, 0, positionStrainA[axis]);
    }
    if (auto dofB = layout.flexibleAtomDof(moleculeB, localAtomB, axes[axis]))
    {
      hessian.addPositionStrain(*dofB, 0, -positionStrainA[axis]);
    }
  }
}

void Minimization::scatterSitePositionStrainIsotropic(GeneralizedHessian& hessian, const MinimizationDofLayout& layout,
                                                      const RigidDerivativeCache& rigidCache, std::size_t moleculeIndex,
                                                      std::size_t localAtom, bool rigid, double sign, double f1,
                                                      double f2, const double3& dr, const double3& drStrainDerivative)
{
  if (hessian.numStrain() != 1)
  {
    return;
  }

  const double drDotStrain = double3::dot(dr, drStrainDerivative);
  const double3 positionValue = f1 * (dr + drStrainDerivative) + (f2 * drDotStrain) * dr;

  if (rigid)
  {
    if (const auto positionBase = layout.atomRigidComDof(moleculeIndex, localAtom))
    {
      hessian.addPositionStrain(*positionBase + 0, 0, sign * positionValue.x);
      hessian.addPositionStrain(*positionBase + 1, 0, sign * positionValue.y);
      hessian.addPositionStrain(*positionBase + 2, 0, sign * positionValue.z);
    }
    if (const auto orientationBase = layout.atomRigidOrientationDof(moleculeIndex, localAtom))
    {
      const RigidAtomDerivatives& derivatives = rigidCache.atom(moleculeIndex, localAtom);
      const std::array<double3, 3> dVec = {derivatives.dVecX, derivatives.dVecY, derivatives.dVecZ};
      for (std::size_t axis = 0; axis < 3; ++axis)
      {
        const double value =
            f2 * drDotStrain * double3::dot(dr, dVec[axis]) + f1 * double3::dot(drStrainDerivative, dVec[axis]);
        hessian.addPositionStrain(*orientationBase + axis, 0, sign * value);
      }
    }
    return;
  }

  constexpr std::array<MinimizationDofAxis, 3> axes = {MinimizationDofAxis::X, MinimizationDofAxis::Y,
                                                       MinimizationDofAxis::Z};
  for (std::size_t axis = 0; axis < 3; ++axis)
  {
    if (const auto dof = layout.flexibleAtomDof(moleculeIndex, localAtom, axes[axis]))
    {
      hessian.addPositionStrain(*dof, 0, sign * (&positionValue.x)[axis]);
    }
  }
}

void Minimization::scatterAtomicStrainStrainIsotropic(GeneralizedHessian& hessian, double f1, double f2,
                                                      const double3& dr, const double3& posA, const double3& comA,
                                                      const double3& posB, const double3& comB, bool rigidA,
                                                      bool rigidB)
{
  if (hessian.numStrain() != 1)
  {
    return;
  }

  // Rigid internal offsets do not scale. Thus, under exp(epsilon) isotropic strain,
  //
  //   dr(epsilon) = exp(epsilon) c + delta,  c = dr - offsetA + offsetB,
  //   dr'(0) = dr''(0) = c.
  //
  // For a radial potential with Cartesian Hessian H = f2 dr dr^T + f1 I,
  //
  //   d2U/dEpsilon2 = c^T H c + grad(U).c
  //                  = f2 (dr.c)^2 + f1 (c.c + dr.c).
  //
  // The last term is the contribution specific to exp(epsilon), absent for linear 1+epsilon scaling.
  const double3 offsetA = rigidA ? posA - comA : double3{};
  const double3 offsetB = rigidB ? posB - comB : double3{};
  const double3 strainDerivative = dr - offsetA + offsetB;
  const double projection = double3::dot(dr, strainDerivative);
  const double value =
      f2 * projection * projection + f1 * (double3::dot(strainDerivative, strainDerivative) + projection);
  hessian.addStrainStrain(0, 0, value);
}
