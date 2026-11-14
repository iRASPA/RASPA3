#include <gtest/gtest.h>

import std;

import int3;
import double3;
import double3x3;
import simd_quatd;
import atom;
import atom_dynamics;
import molecule;
import vdwparameters;
import forcefield;
import framework;
import component;
import connectivity_table;
import intra_molecular_potentials;
import bond_potential;
import bend_potential;
import torsion_potential;
import bond_bend_potential;
import bend_torsion_potential;
import system;
import simulationbox;
import running_energy;
import generalized_hessian;
import minimization_dof_layout;
import minimization_evaluate_derivatives;
import interactions_hessian_framework_molecule;
import interactions_framework_molecule;

namespace
{
void setRigidMoleculePositions(System& system, std::size_t moleculeIndex)
{
  Molecule& molecule = system.moleculeData[moleculeIndex];
  const Component& component = system.components[molecule.componentId];
  const double3x3 rotation = double3x3::buildRotationMatrixInverse(molecule.orientation);
  const double3 com = molecule.centerOfMassPosition;
  for (std::size_t localAtom = 0; localAtom < molecule.numberOfAtoms; ++localAtom)
  {
    system.spanOfMoleculeAtoms()[molecule.atomIndex + localAtom].position =
        com + rotation * component.atoms[localAtom].position;
  }
}

simd_quatd quatFromRotationVector(const double3& omega)
{
  const double angle = std::sqrt(double3::dot(omega, omega));
  if (angle < 1e-30)
  {
    return simd_quatd(0.0, 0.0, 0.0, 1.0);
  }
  return simd_quatd::fromAxisAngle(angle, omega / angle);
}

void useSecondOrderTaylorShiftedLennardJones(ForceField& forceField)
{
  for (VDWParameters& parameters : forceField.data)
  {
    if (parameters.type == VDWParameters::Type::LennardJones)
    {
      parameters.type = VDWParameters::Type::LennardJonesSecondOrderTaylorShifted;
    }
  }
  forceField.preComputeDerivedParameters();
  forceField.preComputePotentialShift();
}

// Checks the single isotropic strain blocks of 'evaluateDerivatives' against central differences of the energy of the
// state 'apply(epsilon, dof, delta)' puts the system in: every coordinate scaled with exp(epsilon) (the box as well)
// and the degree of freedom 'dof' displaced by 'delta'.
void expectIsotropicStrainBlocksMatchFiniteDifference(
    System& system, const MinimizationDofLayout& layout,
    const std::function<void(double, std::optional<std::size_t>, double)>& apply)
{
  const SimulationBox baseBox = system.simulationBox;
  const auto energyAt = [&](double epsilon, std::optional<std::size_t> dof, double delta)
  {
    apply(epsilon, dof, delta);
    const double scale = std::exp(epsilon);
    system.simulationBox = SimulationBox(scale * baseBox.cell.ax, scale * baseBox.cell.by, scale * baseBox.cell.cz);
    GeneralizedHessian scratch(layout.numDofs(), 0);
    DerivativeResults results{.gradient = {}, .hessian = scratch};
    evaluateDerivatives(system, layout, DerivativeCapabilities{.energy = true}, results);
    return results.energy;
  };

  apply(0.0, std::nullopt, 0.0);
  system.simulationBox = baseBox;
  GeneralizedHessian hessian(layout.numDofs(), 1);
  DerivativeResults analytic{.gradient = {}, .hessian = hessian};
  evaluateDerivatives(
      system, layout,
      DerivativeCapabilities{.energy = true, .hessianPositionStrain = true, .hessianStrainStrain = true}, analytic);

  constexpr double strainStep = 1.0e-4;
  constexpr double positionStep = 1.0e-4;
  constexpr double tolerance = 1.0e-4;
  const double numericalStrainStrain =
      (energyAt(strainStep, std::nullopt, 0.0) - 2.0 * energyAt(0.0, std::nullopt, 0.0) +
       energyAt(-strainStep, std::nullopt, 0.0)) /
      (strainStep * strainStep);
  EXPECT_NEAR(numericalStrainStrain, hessian.strainStrain()[0],
              tolerance * std::max({1.0, std::abs(numericalStrainStrain), std::abs(hessian.strainStrain()[0])}))
      << "strain-strain";

  for (std::size_t dof = 0; dof < layout.numDofs(); ++dof)
  {
    const double numerical = (energyAt(strainStep, dof, positionStep) - energyAt(-strainStep, dof, positionStep) -
                              energyAt(strainStep, dof, -positionStep) + energyAt(-strainStep, dof, -positionStep)) /
                             (4.0 * strainStep * positionStep);
    const double value = hessian.positionStrain()[dof];
    EXPECT_NEAR(numerical, value, tolerance * std::max({1.0, std::abs(numerical), std::abs(value)}))
        << "position-strain dof=" << dof;
  }

  apply(0.0, std::nullopt, 0.0);
  system.simulationBox = baseBox;
}
}  // namespace

// A bonded pair of flexible framework atoms and a flexible single-site molecule; charged, the bonded framework pair
// is an Ewald exclusion pair.
void expectFlexibleFrameworkIsotropicStrainBlocksMatchFiniteDifference(bool charged)
{
  ForceField forceField({{"A", false, 12.0, 0.0, 0.0, 6, false}}, {{120.0, 3.4}},
                        ForceField::MixingRule::Lorentz_Berthelot, 10.0, 10.0, 10.0, true, false, charged);
  useSecondOrderTaylorShiftedLennardJones(forceField);
  const SimulationBox box(30.0, 30.0, 30.0);
  const Atom positiveAtom({0.25, 0.25, 0.25}, charged ? 0.6 : 0.0, 1.0, 0, 0, 0, 0, true);
  const Atom negativeAtom({0.25, 0.25, 0.25}, charged ? -0.6 : 0.0, 1.0, 0, 0, 0, 0, true);
  Framework framework(forceField, "flexible-framework", box, 1, {positiveAtom, negativeAtom},
                      {positiveAtom, negativeAtom}, {1, 1, 1});
  framework.rigid = false;
  framework.intraMolecularPotentials.bonds = {BondPotential({0, 1}, BondType::Harmonic, {2000.0, 1.5})};
  framework.intraMolecularImageShifts.bonds = {{{int3(0, 0, 0), int3(0, 0, 0)}}};
  framework.connectivityTable = ConnectivityTable(2);
  framework.connectivityTable[0, 1] = true;
  framework.connectivityTable[1, 0] = true;

  ConnectivityTable connectivity(1);
  Component component(forceField, "single-site", 100.0, 1.0e6, 0.1,
                      {Atom({0.0, 0.0, 0.0}, charged ? 0.4 : 0.0, 1.0, 0, 0, 0, false, false)}, connectivity,
                      Potentials::IntraMolecularPotentials{}, 0, 0);
  component.rigid = false;
  System system(forceField, box, false, 300.0, 1.0e4, 1.0, {framework}, {component}, {}, {1}, 5);

  const std::array<double3, 3> basePositions = {double3(10.0, 10.0, 10.0), double3(11.6, 10.3, 9.8),
                                                double3(13.4, 11.1, 10.4)};
  std::array<Atom*, 3> atoms = {&system.spanOfFrameworkAtoms()[0], &system.spanOfFrameworkAtoms()[1],
                                &system.spanOfMoleculeAtoms()[0]};
  const MinimizationDofLayout layout =
      buildMinimizationDofLayout(system.moleculeData, system.components, system.spanOfFrameworkAtoms().size());
  ASSERT_EQ(layout.numDofs(), 9u);
  std::array<std::size_t, 9> dofOf{};
  for (std::size_t axis = 0; axis < 3; ++axis)
  {
    const MinimizationDofAxis dofAxis = static_cast<MinimizationDofAxis>(axis);
    dofOf[*layout.frameworkAtomDof(0, dofAxis)] = 0 * 3 + axis;
    dofOf[*layout.frameworkAtomDof(1, dofAxis)] = 1 * 3 + axis;
    dofOf[*layout.flexibleAtomDof(0, 0, dofAxis)] = 2 * 3 + axis;
  }

  const auto apply = [&](double epsilon, std::optional<std::size_t> dof, double delta)
  {
    for (std::size_t i = 0; i < atoms.size(); ++i)
    {
      double3 position = basePositions[i];
      if (dof && dofOf[*dof] / 3 == i) (&position.x)[dofOf[*dof] % 3] += delta;
      atoms[i]->position = std::exp(epsilon) * position;
    }
  };
  expectIsotropicStrainBlocksMatchFiniteDifference(system, layout, apply);
}

TEST(minimization_hessian_framework, flexible_framework_isotropic_strain_blocks_match_finite_difference)
{
  expectFlexibleFrameworkIsotropicStrainBlocksMatchFiniteDifference(false);
}

TEST(minimization_hessian_framework, charged_flexible_framework_isotropic_strain_blocks_match_finite_difference)
{
  expectFlexibleFrameworkIsotropicStrainBlocksMatchFiniteDifference(true);
}

TEST(minimization_hessian_framework, semi_flexible_molecule_isotropic_strain_blocks_match_finite_difference)
{
  ForceField forceField({{"O", false, 16.0, 0.0, 0.0, 8, false}}, {{80.0, 3.0}},
                        ForceField::MixingRule::Lorentz_Berthelot, 8.0, 8.0, 8.0, false, false, false);
  ConnectivityTable connectivityTable(4);
  for (std::size_t i = 0; i < 3; ++i)
  {
    connectivityTable[i, i + 1] = true;
    connectivityTable[i + 1, i] = true;
  }
  // every term couples the rigid group {0, 1, 2} to the flexible atom 3, whose rigid offsets do not scale
  Potentials::IntraMolecularPotentials potentials;
  potentials.bonds = {BondPotential({2, 3}, BondType::Harmonic, {1500.0, 1.1})};
  potentials.bends = {BendPotential({1, 2, 3}, BendType::Harmonic, {800.0, 110.0})};
  potentials.torsions = {TorsionPotential({0, 1, 2, 3}, TorsionType::TraPPE, {0.0, 120.0, -30.0, 200.0})};
  const std::uint16_t type = static_cast<std::uint16_t>(*forceField.findPseudoAtom("O"));
  Component component(forceField, "semi-flexible", 64.0, 1.0e6, 0.02,
                      {Atom({-0.6, 0.25, 0.0}, 0.0, 1.0, 0, type, 0, false, false),
                       Atom({0.0, -0.35, 0.0}, 0.0, 1.0, 0, type, 0, false, false),
                       Atom({0.6, 0.25, 0.0}, 0.0, 1.0, 0, type, 0, false, false),
                       Atom({1.2, -0.25, 0.35}, 0.0, 1.0, 0, type, 0, false, false)},
                      connectivityTable, potentials, 0, 0);
  component.buildFragmentGraph({{0, 1, 2}});
  const std::array<double3, 4> basePositions = {double3(5.0, 5.0, 5.0), double3(5.6, 4.4, 5.0), double3(6.2, 5.0, 5.0),
                                                double3(7.0, 4.3, 5.6)};
  System system(forceField, SimulationBox(20.0, 20.0, 20.0), false, 300.0, 1e4, 1.0, {}, {component},
                {{basePositions[0], basePositions[1], basePositions[2], basePositions[3]}}, {0}, 5);
  ASSERT_TRUE(system.components[0].isSemiFlexible());
  system.initializeGroupData();

  const MinimizationDofLayout layout = buildMinimizationDofLayout(system.moleculeData, system.components);
  ASSERT_EQ(layout.numDofs(), 9u);
  const std::size_t comBase = *layout.atomRigidComDof(0, 0);
  const std::size_t orientationBase = *layout.atomRigidOrientationDof(0, 0);
  const std::size_t flexibleBase = *layout.flexibleAtomDof(0, 3, MinimizationDofAxis::X);
  const double3 baseCom = (basePositions[0] + basePositions[1] + basePositions[2]) / 3.0;

  // The rigid group moves its center of mass with the strain and rotates its offsets about lab axes
  // (Rodrigues); the flexible atom scales.
  const auto apply = [&](double epsilon, std::optional<std::size_t> dof, double delta)
  {
    double3 comShift{};
    double3 omega{};
    double3 flexibleShift{};
    if (dof && *dof >= comBase && *dof < comBase + 3) (&comShift.x)[*dof - comBase] = delta;
    if (dof && *dof >= orientationBase && *dof < orientationBase + 3) (&omega.x)[*dof - orientationBase] = delta;
    if (dof && *dof >= flexibleBase && *dof < flexibleBase + 3) (&flexibleShift.x)[*dof - flexibleBase] = delta;
    const double scale = std::exp(epsilon);
    const double angle = std::sqrt(double3::dot(omega, omega));
    std::span<Atom> atoms = system.spanOfMoleculeAtoms();
    for (std::size_t i = 0; i < 3; ++i)
    {
      double3 offset = basePositions[i] - baseCom;
      if (angle > 0.0)
      {
        const double3 axis = omega / angle;
        offset = offset * std::cos(angle) + double3::cross(axis, offset) * std::sin(angle) +
                 axis * (double3::dot(axis, offset) * (1.0 - std::cos(angle)));
      }
      atoms[i].position = scale * (baseCom + comShift) + offset;
    }
    atoms[3].position = scale * (basePositions[3] + flexibleShift);
  };
  expectIsotropicStrainBlocksMatchFiniteDifference(system, layout, apply);
}

TEST(minimization_hessian_framework, semi_flexible_cross_terms_isotropic_strain_blocks_match_finite_difference)
{
  ForceField forceField({{"O", false, 16.0, 0.0, 0.0, 8, false}}, {{80.0, 3.0}},
                        ForceField::MixingRule::Lorentz_Berthelot, 8.0, 8.0, 8.0, false, false, false);
  ConnectivityTable connectivityTable(4);
  for (std::size_t i = 0; i < 3; ++i)
  {
    connectivityTable[i, i + 1] = true;
    connectivityTable[i + 1, i] = true;
  }
  // the cross terms couple the rigid group {0, 1, 2} to the flexible atom 3, whose rigid offsets do not scale
  Potentials::IntraMolecularPotentials potentials;
  potentials.bonds = {BondPotential({2, 3}, BondType::Harmonic, {1500.0, 1.1})};
  // the cross terms (scattered by the semi-flexible cross-term route, which had no isotropic strain blocks)
  potentials.bondBends = {BondBendPotential({1, 2, 3}, BondBendType::CVFF, {109.5, 4000.0, 1.0, 4000.0, 1.1})};
  potentials.bendTorsions = {BendTorsionPotential({0, 1, 2, 3}, BendTorsionType::CVFF, {500.0, 109.5, 109.5})};
  const std::uint16_t type = static_cast<std::uint16_t>(*forceField.findPseudoAtom("O"));
  Component component(forceField, "semi-flexible", 64.0, 1.0e6, 0.02,
                      {Atom({-0.6, 0.25, 0.0}, 0.0, 1.0, 0, type, 0, false, false),
                       Atom({0.0, -0.35, 0.0}, 0.0, 1.0, 0, type, 0, false, false),
                       Atom({0.6, 0.25, 0.0}, 0.0, 1.0, 0, type, 0, false, false),
                       Atom({1.2, -0.25, 0.35}, 0.0, 1.0, 0, type, 0, false, false)},
                      connectivityTable, potentials, 0, 0);
  component.buildFragmentGraph({{0, 1, 2}});
  const std::array<double3, 4> basePositions = {double3(5.0, 5.0, 5.0), double3(5.6, 4.4, 5.0), double3(6.2, 5.0, 5.0),
                                                double3(7.0, 4.3, 5.6)};
  System system(forceField, SimulationBox(20.0, 20.0, 20.0), false, 300.0, 1e4, 1.0, {}, {component},
                {{basePositions[0], basePositions[1], basePositions[2], basePositions[3]}}, {0}, 5);
  ASSERT_TRUE(system.components[0].isSemiFlexible());
  system.initializeGroupData();

  const MinimizationDofLayout layout = buildMinimizationDofLayout(system.moleculeData, system.components);
  ASSERT_EQ(layout.numDofs(), 9u);
  const std::size_t comBase = *layout.atomRigidComDof(0, 0);
  const std::size_t orientationBase = *layout.atomRigidOrientationDof(0, 0);
  const std::size_t flexibleBase = *layout.flexibleAtomDof(0, 3, MinimizationDofAxis::X);
  const double3 baseCom = (basePositions[0] + basePositions[1] + basePositions[2]) / 3.0;

  // The rigid group moves its center of mass with the strain and rotates its offsets about lab axes
  // (Rodrigues); the flexible atom scales.
  const auto apply = [&](double epsilon, std::optional<std::size_t> dof, double delta)
  {
    double3 comShift{};
    double3 omega{};
    double3 flexibleShift{};
    if (dof && *dof >= comBase && *dof < comBase + 3) (&comShift.x)[*dof - comBase] = delta;
    if (dof && *dof >= orientationBase && *dof < orientationBase + 3) (&omega.x)[*dof - orientationBase] = delta;
    if (dof && *dof >= flexibleBase && *dof < flexibleBase + 3) (&flexibleShift.x)[*dof - flexibleBase] = delta;
    const double scale = std::exp(epsilon);
    const double angle = std::sqrt(double3::dot(omega, omega));
    std::span<Atom> atoms = system.spanOfMoleculeAtoms();
    for (std::size_t i = 0; i < 3; ++i)
    {
      double3 offset = basePositions[i] - baseCom;
      if (angle > 0.0)
      {
        const double3 axis = omega / angle;
        offset = offset * std::cos(angle) + double3::cross(axis, offset) * std::sin(angle) +
                 axis * (double3::dot(axis, offset) * (1.0 - std::cos(angle)));
      }
      atoms[i].position = scale * (baseCom + comShift) + offset;
    }
    atoms[3].position = scale * (basePositions[3] + flexibleShift);
  };
  expectIsotropicStrainBlocksMatchFiniteDifference(system, layout, apply);
}

TEST(minimization_hessian_framework, flexible_framework_molecule_cross_block_matches_finite_difference)
{
  ForceField forceField({{"A", false, 12.0, 0.0, 0.0, 6, false}}, {{120.0, 3.4}},
                        ForceField::MixingRule::Lorentz_Berthelot, 10.0, 10.0, 10.0, true, false, false);
  useSecondOrderTaylorShiftedLennardJones(forceField);
  const SimulationBox box(20.0, 20.0, 20.0);
  const Atom frameworkAtom({0.25, 0.25, 0.25}, 0.0, 1.0, 0, 0, 0, 0, true);
  Framework framework(forceField, "flexible-framework", box, 1, {frameworkAtom}, {frameworkAtom}, {1, 1, 1});
  framework.rigid = false;

  ConnectivityTable connectivity(1);
  Component component(forceField, "single-site", 100.0, 1.0e6, 0.1,
                      {Atom({0.0, 0.0, 0.0}, 0.0, 1.0, 0, 0, 0, false, false)}, connectivity,
                      Potentials::IntraMolecularPotentials{}, 0, 0);
  component.rigid = false;
  System system(forceField, box, false, 300.0, 1.0e4, 1.0, {framework}, {component}, {}, {1}, 5);
  system.spanOfFrameworkAtoms()[0].position = {5.0, 5.0, 5.0};
  system.spanOfMoleculeAtoms()[0].position = {8.8, 5.2, 5.1};

  const MinimizationDofLayout layout =
      buildMinimizationDofLayout(system.moleculeData, system.components, system.spanOfFrameworkAtoms().size());
  ASSERT_EQ(layout.numDofs(), 6u);
  GeneralizedHessian hessian(layout.numDofs(), 0);
  Interactions::computeFrameworkMoleculeHessian(
      system.forceField, system.simulationBox, system.moleculeData, system.components, system.spanOfFrameworkAtoms(),
      system.spanOfMoleculeAtoms(), layout, hessian, system.spanOfMoleculeDynamics(), system.spanOfFrameworkDynamics());

  const auto energy = [&]()
  {
    const RunningEnergy value = Interactions::computeFrameworkMoleculeEnergy(
        system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
        system.spanOfFrameworkAtoms(), system.spanOfMoleculeAtoms());
    return value.frameworkMoleculeVDW + value.frameworkMoleculeCharge;
  };
  constexpr double delta = 2.0e-4;
  constexpr double tolerance = 5.0e-2;
  std::array<Atom*, 2> atoms{&system.spanOfFrameworkAtoms()[0], &system.spanOfMoleculeAtoms()[0]};
  const double referenceEnergy = energy();
  for (std::size_t row = 0; row < 6; ++row)
  {
    for (std::size_t column = 0; column < 6; ++column)
    {
      Atom& atomA = *atoms[row / 3];
      Atom& atomB = *atoms[column / 3];
      double numerical{};
      if (row == column)
      {
        (&atomA.position.x)[row % 3] += delta;
        const double plus = energy();
        (&atomA.position.x)[row % 3] -= 2.0 * delta;
        const double minus = energy();
        (&atomA.position.x)[row % 3] += delta;
        numerical = (plus - 2.0 * referenceEnergy + minus) / (delta * delta);
      }
      else
      {
        (&atomA.position.x)[row % 3] += delta;
        (&atomB.position.x)[column % 3] += delta;
        const double ePP = energy();
        (&atomB.position.x)[column % 3] -= 2.0 * delta;
        const double ePM = energy();
        (&atomA.position.x)[row % 3] -= 2.0 * delta;
        const double eMM = energy();
        (&atomB.position.x)[column % 3] += 2.0 * delta;
        const double eMP = energy();
        (&atomA.position.x)[row % 3] += delta;
        (&atomB.position.x)[column % 3] -= delta;
        numerical = (ePP - ePM - eMP + eMM) / (4.0 * delta * delta);
      }
      EXPECT_NEAR(hessian(row, column), numerical, tolerance) << "row=" << row << " column=" << column;
    }
  }
}

TEST(minimization_hessian_framework, flexible_md_direct_gradient_matches_finite_difference)
{
  ForceField forceField({{"P", false, 15.0, 0.8, 0.0, 8, false}, {"N", false, 14.0, -0.4, 0.0, 7, false}},
                        {{80.0, 3.0}, {55.0, 3.2}}, ForceField::MixingRule::Lorentz_Berthelot, 10.0, 10.0, 10.0, true,
                        false, true);
  useSecondOrderTaylorShiftedLennardJones(forceField);
  const SimulationBox box(24.0, 24.0, 24.0);
  const Atom frameworkAtom({0.25, 0.25, 0.25}, 0.8, 1.0, 0, 0, 0, 0, true);
  Framework framework(forceField, "charged-flexible-framework", box, 1, {frameworkAtom}, {frameworkAtom},
                      int3(1, 1, 1));
  framework.rigid = false;
  Component component(forceField, "charged-site", 100.0, 1.0e6, 0.1,
                      {Atom({0.0, 0.0, 0.0}, -0.4, 1.0, 0, 1, 0, false, false)}, ConnectivityTable(1),
                      Potentials::IntraMolecularPotentials{}, 0, 0);
  component.rigid = false;
  System system(forceField, box, false, 300.0, 1.0e4, 1.0, framework, {component}, {}, {1}, 5);
  system.spanOfFrameworkAtoms()[0].position = {7.0, 8.0, 9.0};
  system.spanOfMoleculeAtoms()[0].position = {10.8, 8.4, 8.7};

  std::span<AtomDynamics> moleculeDynamics = system.spanOfMoleculeDynamics();
  std::span<AtomDynamics> frameworkDynamics = system.spanOfFrameworkDynamics();
  const RunningEnergy analytic = Interactions::computeFrameworkMoleculeGradient(
      system.forceField, system.simulationBox, system.spanOfFrameworkAtoms(), system.spanOfMoleculeAtoms(),
      moleculeDynamics, system.interpolationGrids, system.framework, frameworkDynamics);
  EXPECT_NE(analytic.frameworkMoleculeVDW, 0.0);
  EXPECT_NE(analytic.frameworkMoleculeCharge, 0.0);

  const auto energy = [&]()
  {
    return Interactions::computeFrameworkMoleculeEnergy(system.forceField, system.simulationBox,
                                                        system.interpolationGrids, system.framework,
                                                        system.spanOfFrameworkAtoms(), system.spanOfMoleculeAtoms())
        .potentialEnergy();
  };
  std::array<Atom*, 2> atoms{&system.spanOfFrameworkAtoms()[0], &system.spanOfMoleculeAtoms()[0]};
  std::array<AtomDynamics*, 2> dynamics{&frameworkDynamics[0], &moleculeDynamics[0]};
  constexpr double delta = 1.0e-5;
  constexpr double relativeTolerance = 2.0e-5;
  for (std::size_t atom = 0; atom < atoms.size(); ++atom)
  {
    for (std::size_t axis = 0; axis < 3; ++axis)
    {
      double& coordinate = (&atoms[atom]->position.x)[axis];
      coordinate += delta;
      const double plus = energy();
      coordinate -= 2.0 * delta;
      const double minus = energy();
      coordinate += delta;
      const double numerical = (plus - minus) / (2.0 * delta);
      EXPECT_NEAR((&dynamics[atom]->gradient.x)[axis], numerical,
                  relativeTolerance * std::max(1.0, std::abs(numerical)))
          << "atom=" << atom << " axis=" << axis;
    }
  }
  EXPECT_NE(double3::dot(frameworkDynamics[0].gradient, frameworkDynamics[0].gradient), 0.0);
  EXPECT_NE(double3::dot(moleculeDynamics[0].gradient, moleculeDynamics[0].gradient), 0.0);
  EXPECT_NEAR(frameworkDynamics[0].gradient.x + moleculeDynamics[0].gradient.x, 0.0, 1.0e-12);
  EXPECT_NEAR(frameworkDynamics[0].gradient.y + moleculeDynamics[0].gradient.y, 0.0, 1.0e-12);
  EXPECT_NEAR(frameworkDynamics[0].gradient.z + moleculeDynamics[0].gradient.z, 0.0, 1.0e-12);
}

TEST(minimization_hessian_framework, rigid_co2_in_itq29_vdw_matches_finite_difference)
{
  const double delta = 4e-4;
  const double relativeTolerance = 1e-5;
  // Roundoff floor of the FD reference: the total framework energy is large, so the
  // four-point cancellation leaves noise of order eps * E / (4 delta^2).
  const double absoluteTolerance = 1e-5;

  ForceField forceField = ForceField::makeZeoliteForceField(11.8, true, false, true);
  useSecondOrderTaylorShiftedLennardJones(forceField);
  // The finite-difference reference differentiates the (order-0) energy twice; the rotational Hessian entries of
  // the rigid molecule are a heavy cancellation of pair terms, which amplifies the second-derivative error of the
  // tabulated Ewald real-space term above this test's tolerance. Compare against the closed-form energy.
  forceField.useEwaldRealSpaceTable = false;
  // 2x2x2 unit cells so the box exceeds twice the cutoff (minimum image stays consistent).
  Framework framework = Framework::makeITQ29(forceField, int3(2, 2, 2));
  Component co2 = Component::makeCO2(forceField, 0, true);

  System system = System(forceField, std::nullopt, false, 300.0, 1e4, 1.0, {framework}, {co2}, {}, {1}, 5);

  const double3 baseCom(5.93355, 7.93355, 5.93355);
  const simd_quatd baseOrientation =
      simd_quatd::fromAxisAngle(0.4, double3(1.0, -2.0, 1.5) / std::sqrt(7.25)).normalized();

  system.moleculeData[0].centerOfMassPosition = baseCom;
  system.moleculeData[0].orientation = baseOrientation;
  setRigidMoleculePositions(system, 0);

  const MinimizationDofLayout layout = buildMinimizationDofLayout(system.moleculeData, system.components);
  ASSERT_EQ(layout.numDofs(), 6u);

  GeneralizedHessian hessian(layout.numDofs(), 0);
  Interactions::computeFrameworkMoleculeHessian(system.forceField, system.simulationBox, system.moleculeData,
                                                system.components, system.spanOfFrameworkAtoms(),
                                                system.spanOfMoleculeAtoms(), layout, hessian);

  auto energyAtDisplacement = [&](std::span<const double> displacement)
  {
    Molecule& molecule = system.moleculeData[0];
    molecule.centerOfMassPosition = baseCom + double3(displacement[0], displacement[1], displacement[2]);
    const double3 omega(displacement[3], displacement[4], displacement[5]);
    molecule.orientation = (quatFromRotationVector(omega) * baseOrientation).normalized();
    setRigidMoleculePositions(system, 0);

    const RunningEnergy energy = Interactions::computeFrameworkMoleculeEnergy(
        system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
        system.spanOfFrameworkAtoms(), system.spanOfMoleculeAtoms());
    return energy.frameworkMoleculeVDW + energy.frameworkMoleculeCharge;
  };

  std::array<double, 6> displacement{};
  for (std::size_t row = 0; row < 6; ++row)
  {
    for (std::size_t col = 0; col < 6; ++col)
    {
      const auto dofRow = layout.rigidMoleculeDof(0, static_cast<RigidDof>(row));
      const auto dofCol = layout.rigidMoleculeDof(0, static_cast<RigidDof>(col));
      ASSERT_TRUE(dofRow);
      ASSERT_TRUE(dofCol);

      auto finiteDifference = [&](double step)
      {
        displacement.fill(0.0);
        displacement[row] += step;
        displacement[col] += step;
        const double ePP = energyAtDisplacement(displacement);

        displacement.fill(0.0);
        displacement[row] += step;
        displacement[col] -= step;
        const double ePM = energyAtDisplacement(displacement);

        displacement.fill(0.0);
        displacement[row] -= step;
        displacement[col] += step;
        const double eMP = energyAtDisplacement(displacement);

        displacement.fill(0.0);
        displacement[row] -= step;
        displacement[col] -= step;
        const double eMM = energyAtDisplacement(displacement);
        return (ePP - ePM - eMP + eMM) / (4.0 * step * step);
      };
      const double numerical = (4.0 * finiteDifference(delta) - finiteDifference(2.0 * delta)) / 3.0;
      const double analytic = hessian(*dofRow, *dofCol);
      const double scale = std::max({1.0, std::abs(numerical), std::abs(analytic)});
      EXPECT_NEAR(numerical, analytic, absoluteTolerance + relativeTolerance * scale)
          << "row=" << row << " col=" << col;
    }
  }
}

TEST(minimization_hessian_framework, rigid_co2_in_itq29_strain_matches_finite_difference)
{
  const double strainStep = 1.2e-3;
  const double relativeTolerance = 1e-5;
  // Roundoff floor: the total framework-molecule energy is large, so the finite-difference
  // cancellation leaves absolute noise.
  const double absoluteTolerance = 1e-5;

  ForceField forceField = ForceField::makeZeoliteForceField(11.8, true, false, true);
  useSecondOrderTaylorShiftedLennardJones(forceField);
  Framework framework = Framework::makeITQ29(forceField, int3(2, 2, 2));
  Component co2 = Component::makeCO2(forceField, 0, false);

  System system = System(forceField, std::nullopt, false, 300.0, 1e4, 1.0, {framework}, {co2}, {}, {1}, 5);

  const double3 baseCom(5.93355, 7.93355, 5.93355);
  const simd_quatd baseOrientation =
      simd_quatd::fromAxisAngle(0.4, double3(1.0, -2.0, 1.5) / std::sqrt(7.25)).normalized();

  // Base framework positions and cell; under strain both scale with exp(epsilon), matching the
  // affine convention of the analytic strain blocks.
  std::span<Atom> frameworkAtoms = system.spanOfFrameworkAtoms();
  std::vector<double3> baseFrameworkPositions(frameworkAtoms.size());
  for (std::size_t i = 0; i < frameworkAtoms.size(); ++i)
  {
    baseFrameworkPositions[i] = frameworkAtoms[i].position;
  }
  const SimulationBox baseBox = system.simulationBox;

  std::array<double, 6> displacement{};
  auto applyState = [&](double strainExponent)
  {
    const double factor = std::exp(strainExponent);
    Molecule& molecule = system.moleculeData[0];
    molecule.centerOfMassPosition = factor * (baseCom + double3(displacement[0], displacement[1], displacement[2]));
    const double3 omega(displacement[3], displacement[4], displacement[5]);
    molecule.orientation = (quatFromRotationVector(omega) * baseOrientation).normalized();
    setRigidMoleculePositions(system, 0);
    for (std::size_t i = 0; i < frameworkAtoms.size(); ++i)
    {
      frameworkAtoms[i].position = factor * baseFrameworkPositions[i];
    }
  };
  applyState(0.0);

  const MinimizationDofLayout layout = buildMinimizationDofLayout(system.moleculeData, system.components);
  ASSERT_EQ(layout.numDofs(), 6u);

  GeneralizedHessian hessian(layout.numDofs(), 1);
  const RunningEnergy analyticEnergy = Interactions::computeFrameworkMoleculeHessian(
      system.forceField, system.simulationBox, system.moleculeData, system.components, system.spanOfFrameworkAtoms(),
      system.spanOfMoleculeAtoms(), layout, hessian);

  auto energyAt = [&](double strainExponent)
  {
    applyState(strainExponent);
    const double factor = std::exp(strainExponent);
    const SimulationBox strainedBox(factor * baseBox.cell.ax, factor * baseBox.cell.by, factor * baseBox.cell.cz);
    const RunningEnergy energy = Interactions::computeFrameworkMoleculeEnergy(
        system.forceField, strainedBox, system.interpolationGrids, system.framework, system.spanOfFrameworkAtoms(),
        system.spanOfMoleculeAtoms());
    return energy.frameworkMoleculeVDW + energy.frameworkMoleculeCharge;
  };

  // Consistency: the Hessian routine and the FD evaluator must see the same base state.
  EXPECT_NEAR(analyticEnergy.frameworkMoleculeVDW + analyticEnergy.frameworkMoleculeCharge, energyAt(0.0), 1e-8)
      << "base-state energy mismatch between analytic Hessian and FD evaluator";

  // strain-strain block.
  {
    displacement.fill(0.0);
    auto finiteDifference = [&](double step)
    {
      const double eZero = energyAt(0.0);
      const double ePlus = energyAt(step);
      const double eMinus = energyAt(-step);
      return (ePlus - 2.0 * eZero + eMinus) / (step * step);
    };
    const double numerical = (4.0 * finiteDifference(strainStep) - finiteDifference(2.0 * strainStep)) / 3.0;
    const double analytic = hessian.strainStrain()[0];
    const double scale = std::max({1.0, std::abs(numerical), std::abs(analytic)});
    EXPECT_NEAR(numerical, analytic, absoluteTolerance + relativeTolerance * scale) << "strain-strain";
  }

  // position-strain blocks over the rigid COM and orientation DOFs.
  auto gradientAt = [&](double strainExponent)
  {
    displacement.fill(0.0);
    applyState(strainExponent);
    const double factor = std::exp(strainExponent);
    system.simulationBox = SimulationBox(factor * baseBox.cell.ax, factor * baseBox.cell.by, factor * baseBox.cell.cz);

    std::vector<double> gradient(layout.numDofs());
    GeneralizedHessian scratchHessian(layout.numDofs(), 0);
    DerivativeCapabilities capabilities{.energy = false, .gradient = true, .hessianPositionPosition = true};
    DerivativeResults derivatives{.gradient = gradient, .hessian = scratchHessian};
    evaluateDerivatives(system, layout, capabilities, derivatives);
    for (std::size_t slot = 0; slot < 3; ++slot)
    {
      gradient[slot] *= factor;
    }
    return gradient;
  };
  constexpr double gradientStep = 5e-5;
  const std::vector<double> gradientPlus = gradientAt(gradientStep);
  const std::vector<double> gradientMinus = gradientAt(-gradientStep);
  for (std::size_t slot = 0; slot < 6; ++slot)
  {
    const auto dof = layout.rigidMoleculeDof(0, static_cast<RigidDof>(slot));
    ASSERT_TRUE(dof);

    const double numerical = (gradientPlus[*dof] - gradientMinus[*dof]) / (2.0 * gradientStep);
    const double analytic = hessian.positionStrain()[*dof];
    const double scale = std::max({1.0, std::abs(numerical), std::abs(analytic)});
    EXPECT_NEAR(numerical, analytic, absoluteTolerance + relativeTolerance * scale) << "position-strain dof=" << *dof;
  }

  // Restore base state.
  displacement.fill(0.0);
  applyState(0.0);
  system.simulationBox = baseBox;
}

TEST(minimization_hessian_framework, single_framework_atom_rigid_co2_strain_step_convergence)
{
  ForceField forceField = ForceField::makeZeoliteForceField(11.8, true, false, true);
  useSecondOrderTaylorShiftedLennardJones(forceField);
  const Atom frameworkSite(double3(0.25, 0.25, 0.25), 0.0, 1.0, 0, 3, 0, false, false);
  Framework framework(forceField, "single-site", SimulationBox(40.0, 40.0, 40.0), 1, {frameworkSite}, {frameworkSite},
                      int3(1, 1, 1));
  Component co2 = Component::makeCO2(forceField, 0, false);
  System system = System(forceField, std::nullopt, false, 300.0, 1e4, 1.0, {framework}, {co2}, {}, {1}, 5);

  const double3 frameworkBase(10.0, 10.0, 10.0);
  const double3 comBase(14.0, 12.0, 11.0);
  const simd_quatd orientationBase =
      simd_quatd::fromAxisAngle(0.6, double3(1.0, -2.0, 0.7) / std::sqrt(5.49)).normalized();
  system.spanOfFrameworkAtoms()[0].position = frameworkBase;

  std::array<double, 6> displacement{};
  auto applyState = [&](double epsilon)
  {
    const double scale = std::exp(epsilon);
    system.spanOfFrameworkAtoms()[0].position = scale * frameworkBase;
    Molecule& molecule = system.moleculeData[0];
    molecule.centerOfMassPosition = scale * (comBase + double3(displacement[0], displacement[1], displacement[2]));
    molecule.orientation =
        (quatFromRotationVector(double3(displacement[3], displacement[4], displacement[5])) * orientationBase)
            .normalized();
    setRigidMoleculePositions(system, 0);
  };
  applyState(0.0);

  const MinimizationDofLayout layout = buildMinimizationDofLayout(system.moleculeData, system.components);
  ASSERT_EQ(layout.numDofs(), 6u);
  GeneralizedHessian hessian(layout.numDofs(), 1);
  std::span<AtomDynamics> dynamics = system.spanOfMoleculeDynamics();
  Interactions::computeFrameworkMoleculeHessian(system.forceField, system.simulationBox, system.moleculeData,
                                                system.components, system.spanOfFrameworkAtoms(),
                                                system.spanOfMoleculeAtoms(), layout, hessian, dynamics);

  constexpr double gradientStep = 1e-5;
  std::span<Atom> moleculeAtoms = system.spanOfMoleculeAtoms();
  for (std::size_t atom = 0; atom < moleculeAtoms.size(); ++atom)
  {
    for (std::size_t axis = 0; axis < 3; ++axis)
    {
      double& coordinate = (&moleculeAtoms[atom].position.x)[axis];
      coordinate += gradientStep;
      const RunningEnergy plus = Interactions::computeFrameworkMoleculeEnergy(
          system.forceField, SimulationBox(40.0, 40.0, 40.0), system.interpolationGrids, system.framework,
          system.spanOfFrameworkAtoms(), moleculeAtoms);
      coordinate -= 2.0 * gradientStep;
      const RunningEnergy minus = Interactions::computeFrameworkMoleculeEnergy(
          system.forceField, SimulationBox(40.0, 40.0, 40.0), system.interpolationGrids, system.framework,
          system.spanOfFrameworkAtoms(), moleculeAtoms);
      coordinate += gradientStep;
      const double numerical = ((plus.frameworkMoleculeVDW + plus.frameworkMoleculeCharge) -
                                (minus.frameworkMoleculeVDW + minus.frameworkMoleculeCharge)) /
                               (2.0 * gradientStep);
      EXPECT_NEAR((&dynamics[atom].gradient.x)[axis], numerical, 1e-5 * std::max(1.0, std::abs(numerical)));
    }
  }

  auto energyAt = [&](double epsilon)
  {
    applyState(epsilon);
    const double scale = std::exp(epsilon);
    const SimulationBox box(40.0 * scale, 40.0 * scale, 40.0 * scale);
    const RunningEnergy energy = Interactions::computeFrameworkMoleculeEnergy(
        system.forceField, box, system.interpolationGrids, system.framework, system.spanOfFrameworkAtoms(),
        system.spanOfMoleculeAtoms());
    return energy.frameworkMoleculeVDW + energy.frameworkMoleculeCharge;
  };

  // A step sweep distinguishes a bad analytic expression from FD truncation/roundoff.
  // Each estimate must converge to the same analytic value over the stable central-difference range.
  constexpr std::array<double, 4> strainSteps = {5e-4, 4e-4, 3e-4, 2e-4};
  constexpr double positionStep = 2e-5;
  for (double strainStep : strainSteps)
  {
    displacement.fill(0.0);
    const double numericalStrain =
        (energyAt(strainStep) - 2.0 * energyAt(0.0) + energyAt(-strainStep)) / (strainStep * strainStep);
    EXPECT_NEAR(numericalStrain, hessian.strainStrain()[0],
                1e-5 * std::max({1.0, std::abs(numericalStrain), std::abs(hessian.strainStrain()[0])}))
        << "strainStep=" << strainStep;

    for (std::size_t slot = 0; slot < 6; ++slot)
    {
      const auto dof = layout.rigidMoleculeDof(0, static_cast<RigidDof>(slot));
      ASSERT_TRUE(dof);
      displacement.fill(0.0);
      displacement[slot] = positionStep;
      const double ePP = energyAt(strainStep);
      const double ePM = energyAt(-strainStep);
      displacement[slot] = -positionStep;
      const double eMP = energyAt(strainStep);
      const double eMM = energyAt(-strainStep);
      const double numerical = (ePP - ePM - eMP + eMM) / (4.0 * positionStep * strainStep);
      const double analytic = hessian.positionStrain()[*dof];
      EXPECT_NEAR(numerical, analytic, 1e-5 * std::max({1.0, std::abs(numerical), std::abs(analytic)}))
          << "slot=" << slot << " strainStep=" << strainStep;
    }
  }
}
