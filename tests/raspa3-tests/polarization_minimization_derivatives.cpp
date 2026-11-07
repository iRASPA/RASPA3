#include <gtest/gtest.h>

import std;

import atom;
import bond_potential;
import component;
import connectivity_table;
import double3;
import forcefield;
import generalized_hessian;
import interactions_polarization_derivatives;
import intra_molecular_potentials;
import minimization_cell_layout;
import minimization_dof_layout;
import minimization_evaluate_derivatives;
import minimization_generalized_coordinates;
import simulationbox;
import system;

namespace
{
std::pair<std::vector<double>, std::vector<double>> evaluate(System& system, const MinimizationDofLayout& layout,
                                                             bool withHessian)
{
  GeneralizedHessian hessian(layout.numDofs(), 0);
  std::vector<double> gradient(layout.numDofs());
  DerivativeCapabilities capabilities{.energy = true, .gradient = true, .hessianPositionPosition = withHessian};
  DerivativeResults results{.gradient = gradient, .hessian = hessian};
  evaluateDerivatives(system, layout, capabilities, results);
  return {gradient, std::vector<double>(hessian.positionPosition().begin(), hessian.positionPosition().end())};
}

// A flexible charged ion (the source of the field) and a polarizable, uncharged semi-flexible molecule: a rigid
// bent three-site core plus a tail atom on a harmonic bond. The core is placed rotated, off its reference frame.
System makeIonAndPolarizableSemiFlexibleMolecule(CellMinimizationType cellType)
{
  ForceField forceField({{"P", false, 15.0, 0.8, 0.0, 8, false}, {"Q", false, 14.0, 0.0, 1.5, 8, false}},
                        {{60.0, 3.0}, {40.0, 3.2}}, ForceField::MixingRule::Lorentz_Berthelot, 9.0, 9.0, 9.0, true,
                        false, true);
  forceField.computePolarization = true;
  forceField.omitInterPolarization = false;
  forceField.omitInterInteractions = false;
  const std::uint16_t typeP = static_cast<std::uint16_t>(*forceField.findPseudoAtom("P"));
  const std::uint16_t typeQ = static_cast<std::uint16_t>(*forceField.findPseudoAtom("Q"));

  Component ion(forceField, "ion", 100.0, 1e6, 0.2, {Atom({0.0, 0.0, 0.0}, 0.8, 1.0, 0, typeP, 0, false, false)},
                ConnectivityTable(1), Potentials::IntraMolecularPotentials{}, 5, 21);
  ion.rigid = false;

  ConnectivityTable connectivityTable(4);
  for (const auto& [i, j] : {std::pair{0uz, 1uz}, std::pair{1uz, 2uz}, std::pair{2uz, 3uz}})
  {
    connectivityTable[i, j] = true;
    connectivityTable[j, i] = true;
  }
  Potentials::IntraMolecularPotentials potentials;
  potentials.bonds.emplace_back(std::array<std::size_t, 2>{2, 3}, BondType::Harmonic, std::vector<double>{1500.0, 1.5});
  Component molecule(forceField, "semi-flexible", 100.0, 1e6, 0.2,
                     {Atom({-1.2, 0.5, 0.0}, 0.0, 1.0, 0, typeQ, 0, false, false),
                      Atom({0.0, -0.7, 0.0}, 0.0, 1.0, 0, typeQ, 0, false, false),
                      Atom({1.2, 0.5, 0.0}, 0.0, 1.0, 0, typeQ, 0, false, false),
                      Atom({2.4, -0.5, 0.7}, 0.0, 1.0, 0, typeQ, 0, false, false)},
                     connectivityTable, potentials, 5, 21);
  molecule.buildFragmentGraph({{0, 1, 2}});

  System system(forceField, SimulationBox(23.0, 21.0, 19.0), false, 300.0, 1.5e5, 1.0, {}, {ion, molecule}, {}, {1, 1},
                5);
  std::span<Atom> atoms = system.spanOfMoleculeAtoms();
  atoms[0].position = double3(7.3, 9.1, 8.4);
  // the rigid core rotated by 0.4 rad about z around its first atom, translated next to the ion
  const double c = std::cos(0.4);
  const double s = std::sin(0.4);
  const double3 origin(11.0, 11.2, 10.1);
  for (std::size_t k = 0; k < 4; ++k)
  {
    const double3 local = system.components[1].atoms[k].position;
    atoms[1 + k].position = origin + double3(c * local.x - s * local.y, s * local.x + c * local.y, local.z);
  }
  system.cellMinimizationType = cellType;
  return system;
}
}  // namespace

TEST(POLARIZATION_MINIMIZATION_DERIVATIVES, gradient_only_cell_gradient_includes_the_polarization_stress)
{
  // A gradient-only evaluation must give the same cell gradient as one that also builds the Hessian.
  System system = makeIonAndPolarizableSemiFlexibleMolecule(CellMinimizationType::Regular);
  const CellMinimizationLayout cellLayout =
      makeCellMinimizationLayout(system.cellMinimizationType, system.monoclinicAngleType);
  const MinimizationDofLayout layout =
      buildMinimizationDofLayout(system.moleculeData, system.components, 0, cellLayout.size());
  ASSERT_GT(layout.numberOfCellDofs(), 0uz);

  System gradientOnlySystem = system;
  const std::vector<double> gradientOnly = evaluate(gradientOnlySystem, layout, false).first;
  const std::vector<double> withHessian = evaluate(system, layout, true).first;
  for (std::size_t dof = 0; dof < layout.numDofs(); ++dof)
  {
    EXPECT_NEAR(gradientOnly[dof], withHessian[dof], 1e-10 * std::max(1.0, std::abs(withHessian[dof])))
        << "dof " << dof;
  }
}

TEST(POLARIZATION_MINIMIZATION_DERIVATIVES, semi_flexible_rigid_group_polarization_stress_matches_finite_difference)
{
  // Under cell strain the atoms of a rigid group move with the group's center of mass, so the polarization
  // stress must use the offsets from that center of mass.
  System system = makeIonAndPolarizableSemiFlexibleMolecule(CellMinimizationType::Regular);
  const CellMinimizationLayout cellLayout =
      makeCellMinimizationLayout(system.cellMinimizationType, system.monoclinicAngleType);
  const MinimizationDofLayout layout =
      buildMinimizationDofLayout(system.moleculeData, system.components, 0, cellLayout.size());
  const std::vector<std::uint8_t> movable(system.spanOfMoleculeAtoms().size(), 1);

  const Interactions::PolarizationDerivatives reference =
      Interactions::computePolarizationDerivatives(system, movable, cellLayout.bases);
  constexpr double delta = 2.0e-5;
  for (std::size_t a = 0; a < cellLayout.size(); ++a)
  {
    double energies[2];
    for (std::size_t side = 0; side < 2; ++side)
    {
      std::vector<double> displacement(layout.numDofs(), 0.0);
      displacement[*layout.cellDof(a)] = side == 0 ? delta : -delta;
      System deformed = system;
      applyGeneralizedDisplacement(deformed, layout, displacement);
      energies[side] = Interactions::computePolarizationDerivatives(deformed, movable).energy;
    }
    const double numerical = (energies[0] - energies[1]) / (2.0 * delta);
    EXPECT_NEAR(reference.cellGradient[a], numerical,
                3e-4 * std::max({1.0, std::abs(numerical), std::abs(reference.cellGradient[a])}))
        << "strain " << a;
  }
}

TEST(POLARIZATION_MINIMIZATION_DERIVATIVES, semi_flexible_rigid_group_orientation_hessian_matches_finite_difference)
{
  // The orientation degrees of freedom of a rigid group pick up the gradient-curvature term of the polarization
  // energy just like those of a whole rigid molecule.
  System system = makeIonAndPolarizableSemiFlexibleMolecule(CellMinimizationType::Fixed);
  const MinimizationDofLayout layout = buildMinimizationDofLayout(system.moleculeData, system.components, 0, 0);
  const std::optional<std::size_t> orientation = layout.atomRigidOrientationDof(1, 0);
  ASSERT_TRUE(orientation.has_value());

  System referenceSystem = system;
  const std::vector<double> hessian = evaluate(referenceSystem, layout, true).second;
  constexpr double delta = 1.0e-5;
  const std::size_t n = layout.numDofs();
  for (std::size_t axis = 0; axis < 3; ++axis)
  {
    const std::size_t column = *orientation + axis;
    std::vector<double> gradients[2];
    for (std::size_t side = 0; side < 2; ++side)
    {
      std::vector<double> displacement(n, 0.0);
      displacement[column] = side == 0 ? delta : -delta;
      System displaced = system;
      applyGeneralizedDisplacement(displaced, layout, displacement);
      gradients[side] = evaluate(displaced, layout, false).first;
    }
    for (std::size_t row = 0; row < n; ++row)
    {
      const double numerical = (gradients[0][row] - gradients[1][row]) / (2.0 * delta);
      const double analytic = hessian[row * n + column];
      EXPECT_NEAR(analytic, numerical, 1e-4 * std::max({1.0, std::abs(numerical), std::abs(analytic)}))
          << "row " << row << ", orientation axis " << axis;
    }
  }
}
