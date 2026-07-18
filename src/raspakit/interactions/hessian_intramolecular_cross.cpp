module;

module interactions_hessian_intramolecular;

import std;

import double3;
import double3x3;
import units;
import atom;
import component;
import intra_molecular_potentials;
import bond_bond_potential;
import bond_bend_potential;
import bend_bend_potential;
import bond_torsion_potential;
import bend_torsion_potential;
import inversion_bend_potential;
import out_of_plane_bend_potential;
import generalized_hessian;
import minimization_dof_layout;
import minimization_internal_coordinate_hessian;
import minimization_hessian_scatter;
import minimization_rigid_kinematics;
import interactions_intramolecular_hyper_dual;

using namespace Interactions::HyperDual;

namespace
{
// Group-aware scatter of a dense Cartesian per-term Hessian (M-body, M = 3 or 4) for a
// semi-flexible molecule: flexible atoms keep their Cartesian degrees of freedom, atoms driven by
// a rigid group couple through the group's center of mass and orientation. The dense 4x4 block of
// 3x3 sub-blocks is projected onto the generalized degrees of freedom (J^T H J plus the
// gradient-curvature term on rigid orientations), and the strain gradient is corrected to the
// molecular (rigid-body) convention by removing the non-scaling internal-offset virial of every
// rigid-group atom.
void scatterSemiFlexibleCrossTerm(GeneralizedHessian& hessian, const MinimizationDofLayout& layout,
                                  const Minimization::RigidDerivativeCache& rigidCache, std::size_t moleculeIndex,
                                  std::span<const std::size_t> termAtoms, std::span<const double3> positions,
                                  std::span<const double3> gradients,
                                  const std::array<std::array<CoordinateBlock, 4>, 4>& full)
{
  const std::size_t M = termAtoms.size();

  GeneralizedHessian local(3 * M, 0);
  for (std::size_t ti = 0; ti < M; ++ti)
  {
    for (std::size_t tj = 0; tj < M; ++tj)
    {
      for (std::size_t a = 0; a < 3; ++a)
      {
        for (std::size_t b = 0; b < 3; ++b)
        {
          local.add(3 * ti + a, 3 * tj + b, full[ti][tj][a][b]);
        }
      }
    }
  }

  std::vector<Minimization::HessianSite> sites(M);
  for (std::size_t t = 0; t < M; ++t)
  {
    sites[t] = Minimization::makeHessianSite(layout, rigidCache, moleculeIndex, termAtoms[t]);
  }

  Minimization::scatterCartesianTermSites(hessian, sites, gradients, local);
  Minimization::removeRigidOffsetStrainGradientSites(hessian, layout, rigidCache, moleculeIndex, termAtoms, positions,
                                                     gradients);
}
}  // namespace

RunningEnergy Interactions::computeIntraMolecularBondBondHessian(
    std::span<const Molecule> moleculeData, std::span<const Atom> atoms, std::span<const Component> components,
    const MinimizationDofLayout& layout, GeneralizedHessian& hessian, std::span<AtomDynamics> dynamics)
{
  RunningEnergy energies{};

  const Minimization::RigidDerivativeCache rigidCache =
      Minimization::RigidDerivativeCache::build(moleculeData, components, atoms);

  for (std::size_t moleculeIndex = 0; moleculeIndex < moleculeData.size(); ++moleculeIndex)
  {
    const Molecule& molecule = moleculeData[moleculeIndex];
    if (components[molecule.componentId].rigid)
    {
      continue;
    }
    const bool semiFlexible = components[molecule.componentId].isSemiFlexible();

    const Potentials::IntraMolecularPotentials& potentials = components[molecule.componentId].intraMolecularPotentials;
    std::span<const Atom> atom_span = {&atoms[molecule.atomIndex], molecule.numberOfAtoms};
    std::span<AtomDynamics> dynamics_span = {&dynamics[molecule.atomIndex], molecule.numberOfAtoms};

    for (const BondBondPotential& term : potentials.bondBonds)
    {
      const std::size_t A = term.identifiers[0];
      const std::size_t B = term.identifiers[1];
      const std::size_t C = term.identifiers[2];

      auto [energy, gradient, strain] =
          term.potentialEnergyGradientStrain(atom_span[A].position, atom_span[B].position, atom_span[C].position);
      energies.bondBond += energy;
      dynamics_span[A].gradient += gradient[0];
      dynamics_span[B].gradient += gradient[1];
      dynamics_span[C].gradient += gradient[2];
      hessian.strainGradient() += strain;

      const InternalCoordinate rab =
          Minimization::distanceInternalCoordinate(atom_span[A].position, atom_span[B].position);
      const InternalCoordinate rbc =
          Minimization::distanceInternalCoordinate(atom_span[C].position, atom_span[B].position);

      const Dual<2> qab = Dual<2>::variable(rab.value, 0);
      const Dual<2> qbc = Dual<2>::variable(rbc.value, 1);
      Dual<2> U{};
      switch (term.type)
      {
        case BondBondType::CVFF:
        case BondBondType::CFF:
          U = term.parameters[0] * (qab - term.parameters[1]) * (qbc - term.parameters[2]);
          break;
      }

      std::array<double, 3> dUdq{};
      std::array<std::array<double, 3>, 3> d2{};
      extractDerivatives(U, dUdq, d2);

      if (semiFlexible)
      {
        const auto full = assembleCrossCartesianHessian({rab, rbc}, {{{0, 1}}, {{2, 1}}}, dUdq, d2);
        const std::array<std::size_t, 3> ids{A, B, C};
        const std::array<double3, 3> positions{atom_span[A].position, atom_span[B].position, atom_span[C].position};
        const std::array<double3, 3> grads{gradient[0], gradient[1], gradient[2]};
        scatterSemiFlexibleCrossTerm(hessian, layout, rigidCache, moleculeIndex, ids, positions, grads, full);
        continue;
      }

      scatterCrossHessian(hessian, layout, moleculeIndex, {A, B, C}, {rab, rbc}, {{{0, 1}}, {{2, 1}}},
                          {atom_span[A].position, atom_span[B].position, atom_span[C].position, double3{}}, dUdq, d2);
    }
  }

  return energies;
}

RunningEnergy Interactions::computeIntraMolecularBondBendHessian(
    std::span<const Molecule> moleculeData, std::span<const Atom> atoms, std::span<const Component> components,
    const MinimizationDofLayout& layout, GeneralizedHessian& hessian, std::span<AtomDynamics> dynamics)
{
  RunningEnergy energies{};

  const Minimization::RigidDerivativeCache rigidCache =
      Minimization::RigidDerivativeCache::build(moleculeData, components, atoms);

  for (std::size_t moleculeIndex = 0; moleculeIndex < moleculeData.size(); ++moleculeIndex)
  {
    const Molecule& molecule = moleculeData[moleculeIndex];
    if (components[molecule.componentId].rigid)
    {
      continue;
    }
    const bool semiFlexible = components[molecule.componentId].isSemiFlexible();

    const Potentials::IntraMolecularPotentials& potentials = components[molecule.componentId].intraMolecularPotentials;
    std::span<const Atom> atom_span = {&atoms[molecule.atomIndex], molecule.numberOfAtoms};
    std::span<AtomDynamics> dynamics_span = {&dynamics[molecule.atomIndex], molecule.numberOfAtoms};

    for (const BondBendPotential& term : potentials.bondBends)
    {
      const std::size_t A = term.identifiers[0];
      const std::size_t B = term.identifiers[1];
      const std::size_t C = term.identifiers[2];

      auto [energy, gradient, strain] =
          term.potentialEnergyGradientStrain(atom_span[A].position, atom_span[B].position, atom_span[C].position);
      energies.bondBend += energy;
      dynamics_span[A].gradient += gradient[0];
      dynamics_span[B].gradient += gradient[1];
      dynamics_span[C].gradient += gradient[2];
      hessian.strainGradient() += strain;

      const InternalCoordinate rab =
          Minimization::distanceInternalCoordinate(atom_span[A].position, atom_span[B].position);
      const InternalCoordinate rbc =
          Minimization::distanceInternalCoordinate(atom_span[C].position, atom_span[B].position);
      const InternalCoordinate theta =
          Minimization::angleInternalCoordinate(atom_span[A].position, atom_span[B].position, atom_span[C].position);

      const Dual<3> qab = Dual<3>::variable(rab.value, 0);
      const Dual<3> qbc = Dual<3>::variable(rbc.value, 1);
      const Dual<3> qth = Dual<3>::variable(theta.value, 2);
      const auto& p = term.parameters;
      Dual<3> U{};
      switch (term.type)
      {
        case BondBendType::CVFF:
        case BondBendType::CFF:
          U = (qth - p[0]) * (p[1] * (qab - p[2]) + p[3] * (qbc - p[4]));
          break;
        case BondBendType::MM3:
          U = p[0] * ((qab - p[1]) + (qbc - p[2])) * Units::RadiansToDegrees * (qth - p[3]);
          break;
        case BondBendType::TruncatedHarmonic:
        {
          const Dual<3> dth = qth - p[1];
          const Dual<3> exponent = (powDual(qab, 8.0) + powDual(qbc, 8.0)) * (-1.0 / std::pow(p[2], 8));
          U = 0.5 * p[0] * dth * dth * expDual(exponent);
          break;
        }
        case BondBendType::ScreenedHarmonic:
        {
          const Dual<3> dth = qth - p[1];
          const Dual<3> exponent = (qab * (-1.0 / p[2])) + (qbc * (-1.0 / p[3]));
          U = 0.5 * p[0] * dth * dth * expDual(exponent);
          break;
        }
        case BondBendType::ScreenedVessal:
        {
          const double pi = std::numbers::pi;
          const Dual<3> thetaMinusPi = qth - pi;
          const Dual<3> bracket = ((p[1] - pi) * (p[1] - pi)) - (thetaMinusPi * thetaMinusPi);
          const Dual<3> exponent = (qab * (-1.0 / p[2])) + (qbc * (-1.0 / p[3]));
          U = (Dual<3>::constant(p[0] / 8.0) * reciprocal(thetaMinusPi * thetaMinusPi)) * (bracket * bracket) *
              expDual(exponent);
          break;
        }
        case BondBendType::TruncatedVessal:
        {
          const double pi = std::numbers::pi;
          const Dual<3> dth = qth - p[1];
          const Dual<3> dth2 = qth + (p[1] - 2.0 * pi);
          const Dual<3> theTerm = powDual(qth, p[2]) * (dth * dth) * (dth2 * dth2);
          const double constTerm = 0.5 * p[2] * std::pow(pi, p[2] - 1.0) * std::pow(pi - p[1], 3.0);
          const Dual<3> theBracket = theTerm - (dth * dth) * constTerm;
          const Dual<3> exponent = (powDual(qab, 8.0) + powDual(qbc, 8.0)) * (-1.0 / std::pow(p[3], 8));
          U = p[0] * theBracket * expDual(exponent);
          break;
        }
      }

      std::array<double, 3> dUdq{};
      std::array<std::array<double, 3>, 3> d2{};
      extractDerivatives(U, dUdq, d2);

      if (semiFlexible)
      {
        const auto full =
            assembleCrossCartesianHessian({rab, rbc, theta}, {{{0, 1}}, {{2, 1}}, {{0, 1, 2}}}, dUdq, d2);
        const std::array<std::size_t, 3> ids{A, B, C};
        const std::array<double3, 3> positions{atom_span[A].position, atom_span[B].position, atom_span[C].position};
        const std::array<double3, 3> grads{gradient[0], gradient[1], gradient[2]};
        scatterSemiFlexibleCrossTerm(hessian, layout, rigidCache, moleculeIndex, ids, positions, grads, full);
        continue;
      }

      scatterCrossHessian(hessian, layout, moleculeIndex, {A, B, C}, {rab, rbc, theta},
                          {{{0, 1}}, {{2, 1}}, {{0, 1, 2}}},
                          {atom_span[A].position, atom_span[B].position, atom_span[C].position, double3{}}, dUdq, d2);
    }
  }

  return energies;
}

RunningEnergy Interactions::computeIntraMolecularBendBendHessian(
    std::span<const Molecule> moleculeData, std::span<const Atom> atoms, std::span<const Component> components,
    const MinimizationDofLayout& layout, GeneralizedHessian& hessian, std::span<AtomDynamics> dynamics)
{
  RunningEnergy energies{};

  const Minimization::RigidDerivativeCache rigidCache =
      Minimization::RigidDerivativeCache::build(moleculeData, components, atoms);

  for (std::size_t moleculeIndex = 0; moleculeIndex < moleculeData.size(); ++moleculeIndex)
  {
    const Molecule& molecule = moleculeData[moleculeIndex];
    if (components[molecule.componentId].rigid)
    {
      continue;
    }
    const bool semiFlexible = components[molecule.componentId].isSemiFlexible();

    const Potentials::IntraMolecularPotentials& potentials = components[molecule.componentId].intraMolecularPotentials;
    std::span<const Atom> atom_span = {&atoms[molecule.atomIndex], molecule.numberOfAtoms};
    std::span<AtomDynamics> dynamics_span = {&dynamics[molecule.atomIndex], molecule.numberOfAtoms};

    for (const BendBendPotential& term : potentials.bendBends)
    {
      const std::size_t A = term.identifiers[0];
      const std::size_t B = term.identifiers[1];
      const std::size_t C = term.identifiers[2];
      const std::size_t D = term.identifiers[3];

      auto [energy, gradient, strain] = term.potentialEnergyGradientStrain(
          atom_span[A].position, atom_span[B].position, atom_span[C].position, atom_span[D].position);
      energies.bendBend += energy;
      dynamics_span[A].gradient += gradient[0];
      dynamics_span[B].gradient += gradient[1];
      dynamics_span[C].gradient += gradient[2];
      dynamics_span[D].gradient += gradient[3];
      hessian.strainGradient() += strain;

      const InternalCoordinate theta1 =
          Minimization::angleInternalCoordinate(atom_span[A].position, atom_span[B].position, atom_span[C].position);
      const InternalCoordinate theta2 =
          Minimization::angleInternalCoordinate(atom_span[A].position, atom_span[B].position, atom_span[D].position);

      const Dual<2> q1 = Dual<2>::variable(theta1.value, 0);
      const Dual<2> q2 = Dual<2>::variable(theta2.value, 1);
      const auto& p = term.parameters;
      Dual<2> U{};
      switch (term.type)
      {
        case BendBendType::CVFF:
        case BendBendType::CFF:
          U = p[0] * (q1 - p[1]) * (q2 - p[2]);
          break;
        case BendBendType::MM3:
          U = (-p[0] * Units::RadiansToDegrees * Units::RadiansToDegrees) * (q1 - p[1]) * (q2 - p[2]);
          break;
      }

      std::array<double, 3> dUdq{};
      std::array<std::array<double, 3>, 3> d2{};
      extractDerivatives(U, dUdq, d2);

      if (semiFlexible)
      {
        const auto full =
            assembleCrossCartesianHessian({theta1, theta2}, {{{0, 1, 2}}, {{0, 1, 3}}}, dUdq, d2);
        const std::array<std::size_t, 4> ids{A, B, C, D};
        const std::array<double3, 4> positions{atom_span[A].position, atom_span[B].position, atom_span[C].position,
                                               atom_span[D].position};
        const std::array<double3, 4> grads{gradient[0], gradient[1], gradient[2], gradient[3]};
        scatterSemiFlexibleCrossTerm(hessian, layout, rigidCache, moleculeIndex, ids, positions, grads, full);
        continue;
      }

      scatterCrossHessian(hessian, layout, moleculeIndex, {A, B, C, D}, {theta1, theta2}, {{{0, 1, 2}}, {{0, 1, 3}}},
                          {atom_span[A].position, atom_span[B].position, atom_span[C].position, atom_span[D].position},
                          dUdq, d2);
    }
  }

  return energies;
}

RunningEnergy Interactions::computeIntraMolecularBondTorsionHessian(
    std::span<const Molecule> moleculeData, std::span<const Atom> atoms, std::span<const Component> components,
    const MinimizationDofLayout& layout, GeneralizedHessian& hessian, std::span<AtomDynamics> dynamics)
{
  RunningEnergy energies{};

  const Minimization::RigidDerivativeCache rigidCache =
      Minimization::RigidDerivativeCache::build(moleculeData, components, atoms);

  for (std::size_t moleculeIndex = 0; moleculeIndex < moleculeData.size(); ++moleculeIndex)
  {
    const Molecule& molecule = moleculeData[moleculeIndex];
    if (components[molecule.componentId].rigid)
    {
      continue;
    }
    const bool semiFlexible = components[molecule.componentId].isSemiFlexible();

    const Potentials::IntraMolecularPotentials& potentials = components[molecule.componentId].intraMolecularPotentials;
    std::span<const Atom> atom_span = {&atoms[molecule.atomIndex], molecule.numberOfAtoms};
    std::span<AtomDynamics> dynamics_span = {&dynamics[molecule.atomIndex], molecule.numberOfAtoms};

    for (const BondTorsionPotential& term : potentials.bondTorsions)
    {
      const std::size_t A = term.identifiers[0];
      const std::size_t B = term.identifiers[1];
      const std::size_t C = term.identifiers[2];
      const std::size_t D = term.identifiers[3];

      auto [energy, gradient, strain] = term.potentialEnergyGradientStrain(
          atom_span[A].position, atom_span[B].position, atom_span[C].position, atom_span[D].position);
      energies.bondTorsion += energy;
      dynamics_span[A].gradient += gradient[0];
      dynamics_span[B].gradient += gradient[1];
      dynamics_span[C].gradient += gradient[2];
      dynamics_span[D].gradient += gradient[3];
      hessian.strainGradient() += strain;

      const InternalCoordinate rbc =
          Minimization::distanceInternalCoordinate(atom_span[B].position, atom_span[C].position);
      const InternalCoordinate cphi = Minimization::dihedralInternalCoordinate(
          atom_span[A].position, atom_span[B].position, atom_span[C].position, atom_span[D].position);

      const Dual<2> qr = Dual<2>::variable(rbc.value, 0);
      const Dual<2> qc = Dual<2>::variable(cphi.value, 1);
      const auto& p = term.parameters;
      Dual<2> U{};
      switch (term.type)
      {
        case BondTorsionType::MM3:
        {
          const Dual<2> temp = qr - p[3];
          const Dual<2> c2 = qc * qc;
          U = p[0] * temp * qc + p[1] * temp * (2.0 * c2 - 1.0) + p[2] * temp * ((4.0 * c2 * qc) - 3.0 * qc);
          break;
        }
      }

      std::array<double, 3> dUdq{};
      std::array<std::array<double, 3>, 3> d2{};
      extractDerivatives(U, dUdq, d2);

      if (semiFlexible)
      {
        const auto full =
            assembleCrossCartesianHessian({rbc, cphi}, {{{1, 2}}, {{0, 1, 2, 3}}}, dUdq, d2);
        const std::array<std::size_t, 4> ids{A, B, C, D};
        const std::array<double3, 4> positions{atom_span[A].position, atom_span[B].position, atom_span[C].position,
                                               atom_span[D].position};
        const std::array<double3, 4> grads{gradient[0], gradient[1], gradient[2], gradient[3]};
        scatterSemiFlexibleCrossTerm(hessian, layout, rigidCache, moleculeIndex, ids, positions, grads, full);
        continue;
      }

      scatterCrossHessian(hessian, layout, moleculeIndex, {A, B, C, D}, {rbc, cphi}, {{{1, 2}}, {{0, 1, 2, 3}}},
                          {atom_span[A].position, atom_span[B].position, atom_span[C].position, atom_span[D].position},
                          dUdq, d2);
    }
  }

  return energies;
}

RunningEnergy Interactions::computeIntraMolecularBendTorsionHessian(
    std::span<const Molecule> moleculeData, std::span<const Atom> atoms, std::span<const Component> components,
    const MinimizationDofLayout& layout, GeneralizedHessian& hessian, std::span<AtomDynamics> dynamics)
{
  RunningEnergy energies{};

  const Minimization::RigidDerivativeCache rigidCache =
      Minimization::RigidDerivativeCache::build(moleculeData, components, atoms);

  for (std::size_t moleculeIndex = 0; moleculeIndex < moleculeData.size(); ++moleculeIndex)
  {
    const Molecule& molecule = moleculeData[moleculeIndex];
    if (components[molecule.componentId].rigid)
    {
      continue;
    }
    const bool semiFlexible = components[molecule.componentId].isSemiFlexible();

    const Potentials::IntraMolecularPotentials& potentials = components[molecule.componentId].intraMolecularPotentials;
    std::span<const Atom> atom_span = {&atoms[molecule.atomIndex], molecule.numberOfAtoms};
    std::span<AtomDynamics> dynamics_span = {&dynamics[molecule.atomIndex], molecule.numberOfAtoms};

    for (const BendTorsionPotential& term : potentials.bendTorsions)
    {
      const std::size_t A = term.identifiers[0];
      const std::size_t B = term.identifiers[1];
      const std::size_t C = term.identifiers[2];
      const std::size_t D = term.identifiers[3];

      auto [energy, gradient, strain] = term.potentialEnergyGradientStrain(
          atom_span[A].position, atom_span[B].position, atom_span[C].position, atom_span[D].position);
      energies.bendTorsion += energy;
      dynamics_span[A].gradient += gradient[0];
      dynamics_span[B].gradient += gradient[1];
      dynamics_span[C].gradient += gradient[2];
      dynamics_span[D].gradient += gradient[3];
      hessian.strainGradient() += strain;

      const InternalCoordinate theta1 =
          Minimization::angleInternalCoordinate(atom_span[A].position, atom_span[B].position, atom_span[C].position);
      const InternalCoordinate theta2 =
          Minimization::angleInternalCoordinate(atom_span[B].position, atom_span[C].position, atom_span[D].position);

      // Smoothed uses the signed angle phi; all other types are polynomials in cos(phi).
      const bool usesPhi = (term.type == BendTorsionType::Smoothed);
      const InternalCoordinate torsion =
          usesPhi ? Minimization::dihedralAngleInternalCoordinate(atom_span[A].position, atom_span[B].position,
                                                                  atom_span[C].position, atom_span[D].position)
                  : Minimization::dihedralInternalCoordinate(atom_span[A].position, atom_span[B].position,
                                                             atom_span[C].position, atom_span[D].position);

      const Dual<3> q1 = Dual<3>::variable(theta1.value, 0);
      const Dual<3> q2 = Dual<3>::variable(theta2.value, 1);
      const Dual<3> qt = Dual<3>::variable(torsion.value, 2);  // cos(phi) or phi
      const auto& p = term.parameters;
      Dual<3> U{};
      switch (term.type)
      {
        case BendTorsionType::CVFF:
        case BendTorsionType::CFF:
          U = p[0] * (q1 - p[1]) * (q2 - p[2]) * qt;
          break;
        case BendTorsionType::Smoothed:
          U = p[0] * (1.0 + cosDual(p[1] * qt - p[2])) * smoothingDual(q1) * smoothingDual(q2);
          break;
        case BendTorsionType::SmoothedThreeCosine:
        {
          const Dual<3> c2 = qt * qt;
          U = (0.5 * p[0] * (1.0 + qt) + p[1] * (1.0 - c2) + 0.5 * p[2] * (1.0 - 3.0 * qt + 4.0 * qt * c2)) *
              smoothingDual(q1) * smoothingDual(q2);
          break;
        }
        case BendTorsionType::Nicholas:
        {
          const Dual<3> c2 = qt * qt;
          U = (0.5 * p[0] * (1.0 + qt) + p[1] * (1.0 - c2) + 0.5 * p[2] * (1.0 - 3.0 * qt + 4.0 * qt * c2)) *
              smoothingDual(q1);
          break;
        }
        case BendTorsionType::SmoothedCFF:
        {
          const Dual<3> c2 = qt * qt;
          U = (p[0] * (1.0 - qt) + 2.0 * p[1] * (1.0 - c2) + p[2] * (1.0 + 3.0 * qt - 4.0 * qt * c2)) *
              smoothingDual(q1) * smoothingDual(q2);
          break;
        }
        case BendTorsionType::SmoothedCFF2:
          U = (p[0] * (1.0 + qt) + p[2] + qt * (-3.0 * p[2] + 2.0 * qt * (p[1] + 2.0 * p[2] * qt))) *
              smoothingDual(q1) * smoothingDual(q2);
          break;
        case BendTorsionType::SmoothedCFF3:
          U = p[0] * (q1 - p[1]) * (q2 - p[2]) * qt * smoothingDual(q1) * smoothingDual(q2);
          break;
      }

      std::array<double, 3> dUdq{};
      std::array<std::array<double, 3>, 3> d2{};
      extractDerivatives(U, dUdq, d2);

      if (semiFlexible)
      {
        const auto full = assembleCrossCartesianHessian(
            {theta1, theta2, torsion}, {{{0, 1, 2}}, {{1, 2, 3}}, {{0, 1, 2, 3}}}, dUdq, d2);
        const std::array<std::size_t, 4> ids{A, B, C, D};
        const std::array<double3, 4> positions{atom_span[A].position, atom_span[B].position, atom_span[C].position,
                                               atom_span[D].position};
        const std::array<double3, 4> grads{gradient[0], gradient[1], gradient[2], gradient[3]};
        scatterSemiFlexibleCrossTerm(hessian, layout, rigidCache, moleculeIndex, ids, positions, grads, full);
        continue;
      }

      scatterCrossHessian(hessian, layout, moleculeIndex, {A, B, C, D}, {theta1, theta2, torsion},
                          {{{0, 1, 2}}, {{1, 2, 3}}, {{0, 1, 2, 3}}},
                          {atom_span[A].position, atom_span[B].position, atom_span[C].position, atom_span[D].position},
                          dUdq, d2);
    }
  }

  return energies;
}

RunningEnergy Interactions::computeIntraMolecularInversionBendHessian(
    std::span<const Molecule> moleculeData, std::span<const Atom> atoms, std::span<const Component> components,
    const MinimizationDofLayout& layout, GeneralizedHessian& hessian, std::span<AtomDynamics> dynamics)
{
  RunningEnergy energies{};

  const Minimization::RigidDerivativeCache rigidCache =
      Minimization::RigidDerivativeCache::build(moleculeData, components, atoms);

  for (std::size_t moleculeIndex = 0; moleculeIndex < moleculeData.size(); ++moleculeIndex)
  {
    const Molecule& molecule = moleculeData[moleculeIndex];
    if (components[molecule.componentId].rigid)
    {
      continue;
    }
    const bool semiFlexible = components[molecule.componentId].isSemiFlexible();

    const Potentials::IntraMolecularPotentials& potentials = components[molecule.componentId].intraMolecularPotentials;
    std::span<const Atom> atom_span = {&atoms[molecule.atomIndex], molecule.numberOfAtoms};
    std::span<AtomDynamics> dynamics_span = {&dynamics[molecule.atomIndex], molecule.numberOfAtoms};

    for (const InversionBendPotential& term : potentials.inversionBends)
    {
      const std::size_t A = term.identifiers[0];
      const std::size_t B = term.identifiers[1];
      const std::size_t C = term.identifiers[2];
      const std::size_t D = term.identifiers[3];

      auto [energy, gradient, strain] = term.potentialEnergyGradientStrain(
          atom_span[A].position, atom_span[B].position, atom_span[C].position, atom_span[D].position);
      energies.inversionBend += energy;
      dynamics_span[A].gradient += gradient[0];
      dynamics_span[B].gradient += gradient[1];
      dynamics_span[C].gradient += gradient[2];
      dynamics_span[D].gradient += gradient[3];
      hessian.strainGradient() += strain;

      // Differentiate the full Cartesian inversion-bend energy through 12 hyper-dual variables
      // (three per atom); the inversion angle chi is not a distance/angle/dihedral coordinate.
      const std::array<double3, 4> positions = {atom_span[A].position, atom_span[B].position, atom_span[C].position,
                                                atom_span[D].position};
      auto seed = [&](std::size_t atomSlot) -> DualVec<12>
      {
        return {Dual<12>::variable((&positions[atomSlot].x)[0], atomSlot * 3 + 0),
                Dual<12>::variable((&positions[atomSlot].x)[1], atomSlot * 3 + 1),
                Dual<12>::variable((&positions[atomSlot].x)[2], atomSlot * 3 + 2)};
      };
      const DualVec<12> posA = seed(0);
      const DualVec<12> posB = seed(1);
      const DualVec<12> posC = seed(2);
      const DualVec<12> posD = seed(3);

      const DualVec<12> Rab = posA - posB;
      const DualVec<12> Rbc = posC - posB;
      const DualVec<12> Rbd = posD - posB;
      const DualVec<12> Rcd = posD - posC;
      const DualVec<12> Rad = posD - posA;

      const bool useAcdPlane = term.type == InversionBendType::Harmonic2 ||
                               term.type == InversionBendType::HarmonicCosine2 ||
                               term.type == InversionBendType::Planar2 || term.type == InversionBendType::MM3;

      Dual<12> c{};
      if (useAcdPlane)
      {
        const Dual<12> temp = dotDual(Rad, Rcd);
        c = dotDual(Rcd, Rcd) * dotDual(Rad, Rad) - temp * temp;
      }
      else
      {
        const Dual<12> temp = dotDual(Rbc, Rbd);
        c = dotDual(Rbc, Rbc) * dotDual(Rbd, Rbd) - temp * temp;
      }

      const Dual<12> rab2 = dotDual(Rab, Rab);
      const Dual<12> e = dotDual(Rab, crossDual(Rbd, Rbc));
      const Dual<12> cosChi = sqrtDual(rab2 - (e * e) / c) / sqrtDual(rab2);

      const auto& p = term.parameters;
      Dual<12> U{};
      switch (term.type)
      {
        case InversionBendType::Harmonic:
        case InversionBendType::Harmonic2:
        {
          const Dual<12> t = acosDual(cosChi) - p[1];
          U = 0.5 * p[0] * t * t;
          break;
        }
        case InversionBendType::HarmonicCosine:
        case InversionBendType::HarmonicCosine2:
        {
          const Dual<12> t = cosChi - p[1];
          U = 0.5 * p[0] * t * t;
          break;
        }
        case InversionBendType::Planar:
        case InversionBendType::Planar2:
          U = p[0] * (1.0 - cosChi);
          break;
        case InversionBendType::MM3:
        {
          const Dual<12> t = (acosDual(cosChi) - p[1]) * Units::RadiansToDegrees;
          const Dual<12> t2 = t * t;
          U = p[0] * t2 * (1.0 - 0.014 * t + 5.6e-5 * t2 - 7.0e-7 * (t * t2) + 2.2e-8 * (t2 * t2));
          break;
        }
      }

      std::array<double3, 4> analyticGradient{};
      std::array<std::array<CoordinateBlock, 4>, 4> analyticHessian{};
      for (std::size_t i = 0; i < 4; ++i)
      {
        analyticGradient[i] = {U.g[i * 3 + 0], U.g[i * 3 + 1], U.g[i * 3 + 2]};
        for (std::size_t j = 0; j < 4; ++j)
        {
          for (std::size_t a = 0; a < 3; ++a)
          {
            for (std::size_t b = 0; b < 3; ++b)
            {
              analyticHessian[i][j][a][b] = U.h[i * 3 + a][j * 3 + b];
            }
          }
        }
      }

      if (semiFlexible)
      {
        const std::array<std::size_t, 4> ids{A, B, C, D};
        scatterSemiFlexibleCrossTerm(hessian, layout, rigidCache, moleculeIndex, ids, positions, analyticGradient,
                                     analyticHessian);
        continue;
      }

      scatterCartesianFourBody(hessian, layout, moleculeIndex, {A, B, C, D}, positions, analyticGradient,
                               analyticHessian);
    }
  }

  return energies;
}

RunningEnergy Interactions::computeIntraMolecularOutOfPlaneBendHessian(
    std::span<const Molecule> moleculeData, std::span<const Atom> atoms, std::span<const Component> components,
    const MinimizationDofLayout& layout, GeneralizedHessian& hessian, std::span<AtomDynamics> dynamics)
{
  RunningEnergy energies{};

  const Minimization::RigidDerivativeCache rigidCache =
      Minimization::RigidDerivativeCache::build(moleculeData, components, atoms);

  for (std::size_t moleculeIndex = 0; moleculeIndex < moleculeData.size(); ++moleculeIndex)
  {
    const Molecule& molecule = moleculeData[moleculeIndex];
    if (components[molecule.componentId].rigid)
    {
      continue;
    }
    const bool semiFlexible = components[molecule.componentId].isSemiFlexible();

    const Potentials::IntraMolecularPotentials& potentials = components[molecule.componentId].intraMolecularPotentials;
    std::span<const Atom> atom_span = {&atoms[molecule.atomIndex], molecule.numberOfAtoms};
    std::span<AtomDynamics> dynamics_span = {&dynamics[molecule.atomIndex], molecule.numberOfAtoms};

    for (const OutOfPlaneBendPotential& term : potentials.outOfPlaneBends)
    {
      const std::size_t A = term.identifiers[0];
      const std::size_t B = term.identifiers[1];
      const std::size_t C = term.identifiers[2];
      const std::size_t D = term.identifiers[3];

      auto [energy, gradient, strain] = term.potentialEnergyGradientStrain(
          atom_span[A].position, atom_span[B].position, atom_span[C].position, atom_span[D].position);
      energies.outOfPlaneBend += energy;
      dynamics_span[A].gradient += gradient[0];
      dynamics_span[B].gradient += gradient[1];
      dynamics_span[C].gradient += gradient[2];
      dynamics_span[D].gradient += gradient[3];
      hessian.strainGradient() += strain;

      // Differentiate the full Cartesian out-of-plane-bend energy through 12 hyper-dual variables
      // (three per atom). The energy expression below must mirror OutOfPlaneBendPotential::calculateEnergy:
      // keeping the two in lockstep makes the Hessian automatically correct for every functional form.
      const std::array<double3, 4> positions = {atom_span[A].position, atom_span[B].position, atom_span[C].position,
                                                atom_span[D].position};
      auto seed = [&](std::size_t atomSlot) -> DualVec<12>
      {
        return {Dual<12>::variable((&positions[atomSlot].x)[0], atomSlot * 3 + 0),
                Dual<12>::variable((&positions[atomSlot].x)[1], atomSlot * 3 + 1),
                Dual<12>::variable((&positions[atomSlot].x)[2], atomSlot * 3 + 2)};
      };
      [[maybe_unused]] const DualVec<12> posA = seed(0);
      [[maybe_unused]] const DualVec<12> posB = seed(1);
      [[maybe_unused]] const DualVec<12> posC = seed(2);
      [[maybe_unused]] const DualVec<12> posD = seed(3);

      Dual<12> U{};
      switch (term.type)
      {
        case OutOfPlaneBendType::Harmonic:
          // OutOfPlaneBendPotential::calculateEnergy currently returns 0 for the harmonic form; the
          // hyper-dual energy mirrors that, so this term contributes a zero gradient and Hessian.
          U = Dual<12>::constant(0.0);
          break;
      }

      std::array<double3, 4> analyticGradient{};
      std::array<std::array<CoordinateBlock, 4>, 4> analyticHessian{};
      for (std::size_t i = 0; i < 4; ++i)
      {
        analyticGradient[i] = {U.g[i * 3 + 0], U.g[i * 3 + 1], U.g[i * 3 + 2]};
        for (std::size_t j = 0; j < 4; ++j)
        {
          for (std::size_t a = 0; a < 3; ++a)
          {
            for (std::size_t b = 0; b < 3; ++b)
            {
              analyticHessian[i][j][a][b] = U.h[i * 3 + a][j * 3 + b];
            }
          }
        }
      }

      if (semiFlexible)
      {
        const std::array<std::size_t, 4> ids{A, B, C, D};
        scatterSemiFlexibleCrossTerm(hessian, layout, rigidCache, moleculeIndex, ids, positions, analyticGradient,
                                     analyticHessian);
        continue;
      }

      scatterCartesianFourBody(hessian, layout, moleculeIndex, {A, B, C, D}, positions, analyticGradient,
                               analyticHessian);
    }
  }

  return energies;
}
