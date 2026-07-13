#include <gtest/gtest.h>

import std;

import double3;
import double3x3;
import simd_quatd;
import atom;
import atom_dynamics;
import molecule;
import forcefield;
import component;
import system;
import simulationbox;
import connectivity_table;
import intra_molecular_potentials;
import running_energy;
import generalized_hessian;
import minimization_dof_layout;
import interactions_ewald;
import interactions_hessian_ewald;

namespace
{
// Scratch buffers plus a callable evaluating the Ewald energy (Fourier + exclusion) for the
// current atom positions and a given simulation box. The self energy is omitted: it is
// independent of all degrees of freedom and its large constant offset would dominate the
// cancellation noise of the finite differences.
struct EwaldEnergyEvaluator
{
  std::vector<std::complex<double>> eik_x;
  std::vector<std::complex<double>> eik_y;
  std::vector<std::complex<double>> eik_z;
  std::vector<std::complex<double>> eik_xy;
  std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>> fixedFrameworkStoredEik;
  std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>> storedEik;

  double operator()(const System &system, const SimulationBox &box)
  {
    RunningEnergy energy = Interactions::computeEwaldFourierEnergy(
        eik_x, eik_y, eik_z, eik_xy, fixedFrameworkStoredEik, storedEik, system.forceField, box, system.components,
        system.numberOfMoleculesPerComponent, system.spanOfMoleculeAtoms(), 0.0);
    return energy.ewald_fourier + energy.ewald_exclusion;
  }
};

ForceField makeChargedForceField()
{
  return ForceField({{"P", false, 15.0, 0.8, 0.0, 8, false},
                     {"N", false, 14.0, -0.4, 0.0, 8, false}},
                    {{60.0, 3.0}, {40.0, 3.2}},
                    ForceField::MixingRule::Lorentz_Berthelot, 11.0, 11.0, 11.0, true, false, true);
}

// Two flexible bent 3-site molecules with zero net charge each.
System makeFlexiblePairSystem(const ForceField &forceField)
{
  ConnectivityTable connectivityTable(3);
  connectivityTable[0, 1] = true;
  connectivityTable[1, 0] = true;
  connectivityTable[1, 2] = true;
  connectivityTable[2, 1] = true;

  Component component = Component(forceField, "flexibleIon", 100.0, 1e6, 0.2,
                                  {Atom({-0.9, 0.4, 0.1}, -0.4, 1.0, 0, 1, 0, false, false),
                                   Atom({0.0, 0.0, 0.0}, 0.8, 1.0, 0, 0, 0, false, false),
                                   Atom({1.0, 0.3, -0.2}, -0.4, 1.0, 0, 1, 0, false, false)},
                                  connectivityTable, Potentials::IntraMolecularPotentials{}, 5, 21);
  component.rigid = false;

  System system =
      System(forceField, SimulationBox(25.0, 25.0, 25.0), false, 300.0, 1e4, 1.0, {}, {component}, {}, {2}, 5);

  const std::array<double3, 6> positions = {double3(9.1, 10.4, 10.1),  double3(10.0, 10.0, 10.0),
                                            double3(11.0, 10.3, 9.8),  double3(12.6, 11.9, 10.7),
                                            double3(13.5, 11.5, 10.6), double3(14.4, 11.8, 11.3)};
  std::span<Atom> atoms = system.spanOfMoleculeAtoms();
  for (std::size_t i = 0; i < positions.size(); ++i)
  {
    atoms[i].position = positions[i];
  }
  return system;
}

System makeSingleIonSystem(const ForceField &forceField)
{
  Component component =
      Component(forceField, "singleIon", 100.0, 1e6, 0.2,
                {Atom({0.0, 0.0, 0.0}, 0.8, 1.0, 0, 0, 0, false, false)}, ConnectivityTable(1),
                Potentials::IntraMolecularPotentials{}, 5, 21);
  component.rigid = false;

  System system =
      System(forceField, SimulationBox(25.0, 25.0, 25.0), false, 300.0, 1e4, 1.0, {}, {component}, {}, {1}, 5);
  system.spanOfMoleculeAtoms()[0].position = double3(8.3, 11.1, 13.7);
  return system;
}

struct FlexibleDofLabel
{
  std::size_t atom;  // global atom index
  std::size_t axis;
};

std::vector<FlexibleDofLabel> buildFlexibleDofLabels(const MinimizationDofLayout &layout, const System &system)
{
  std::vector<FlexibleDofLabel> labels(layout.numDofs());
  for (std::size_t moleculeIndex = 0; moleculeIndex < system.moleculeData.size(); ++moleculeIndex)
  {
    const Molecule &molecule = system.moleculeData[moleculeIndex];
    for (std::size_t localAtom = 0; localAtom < molecule.numberOfAtoms; ++localAtom)
    {
      for (std::size_t axis = 0; axis < 3; ++axis)
      {
        const auto dof = layout.flexibleAtomDof(moleculeIndex, localAtom, static_cast<MinimizationDofAxis>(axis));
        EXPECT_TRUE(dof.has_value());
        labels[*dof] = {molecule.atomIndex + localAtom, axis};
      }
    }
  }
  return labels;
}

void setRigidMoleculePositions(System &system, std::size_t moleculeIndex)
{
  Molecule &molecule = system.moleculeData[moleculeIndex];
  const Component &component = system.components[molecule.componentId];
  const double3x3 rotation = double3x3::buildRotationMatrixInverse(molecule.orientation);
  const double3 com = molecule.centerOfMassPosition;
  for (std::size_t localAtom = 0; localAtom < molecule.numberOfAtoms; ++localAtom)
  {
    system.spanOfMoleculeAtoms()[molecule.atomIndex + localAtom].position =
        com + rotation * component.atoms[localAtom].position;
  }
}

simd_quatd quatFromRotationVector(const double3 &omega)
{
  const double angle = std::sqrt(double3::dot(omega, omega));
  if (angle < 1e-30)
  {
    return simd_quatd(0.0, 0.0, 0.0, 1.0);
  }
  return simd_quatd::fromAxisAngle(angle, omega / angle);
}

struct RigidState
{
  double3 com;
  simd_quatd orientation;
};

// Orientation displacements premultiply the base orientation with the exponential map of the
// rotation vector; the center of mass scales with the isotropic strain.
void applyDisplacedState(System &system, std::size_t moleculeIndex, const RigidState &base,
                         std::span<const double> displacement, double strainFactor)
{
  Molecule &molecule = system.moleculeData[moleculeIndex];
  molecule.centerOfMassPosition =
      strainFactor * (base.com + double3(displacement[0], displacement[1], displacement[2]));
  const double3 omega(displacement[3], displacement[4], displacement[5]);
  molecule.orientation = (quatFromRotationVector(omega) * base.orientation).normalized();
  setRigidMoleculePositions(system, moleculeIndex);
}
}  // namespace

TEST(minimization_hessian_ewald, flexible_charged_pair_matches_finite_difference)
{
  // delta balances roundoff (the large erf exclusion energy cancels over 4*delta^2) against
  // O(delta^2) truncation of the oscillatory Fourier terms.
  const double delta = 5e-4;
  const double relativeTolerance = 5e-5;

  ForceField forceField = makeChargedForceField();
  System system = makeFlexiblePairSystem(forceField);

  const MinimizationDofLayout layout = buildMinimizationDofLayout(system.moleculeData, system.components);
  ASSERT_EQ(layout.numDofs(), 18u);

  GeneralizedHessian hessian(layout.numDofs(), 1);
  std::vector<AtomDynamics> dynamics(system.spanOfMoleculeAtoms().size());
  Interactions::computeEwaldFourierHessian(system, layout, hessian, dynamics);

  EwaldEnergyEvaluator evaluator{};
  auto energy = [&]() { return evaluator(system, system.simulationBox); };

  std::span<Atom> atoms = system.spanOfMoleculeAtoms();
  const std::vector<FlexibleDofLabel> labels = buildFlexibleDofLabels(layout, system);

  auto perturb = [&](const FlexibleDofLabel &label, double amount)
  { (&atoms[label.atom].position.x)[label.axis] += amount; };

  // Gradient check.
  for (std::size_t dof = 0; dof < labels.size(); ++dof)
  {
    perturb(labels[dof], delta);
    const double ePlus = energy();
    perturb(labels[dof], -2.0 * delta);
    const double eMinus = energy();
    perturb(labels[dof], delta);
    const double numerical = (ePlus - eMinus) / (2.0 * delta);
    const double analytic = (&dynamics[labels[dof].atom].gradient.x)[labels[dof].axis];
    const double scale = std::max({1.0, std::abs(numerical), std::abs(analytic)});
    EXPECT_NEAR(numerical, analytic, relativeTolerance * scale) << "gradient dof=" << dof;
  }

  // Position-position block.
  for (std::size_t row = 0; row < labels.size(); ++row)
  {
    for (std::size_t col = 0; col < labels.size(); ++col)
    {
      perturb(labels[row], delta);
      perturb(labels[col], delta);
      const double ePP = energy();
      perturb(labels[col], -2.0 * delta);
      const double ePM = energy();
      perturb(labels[row], -2.0 * delta);
      const double eMM = energy();
      perturb(labels[col], 2.0 * delta);
      const double eMP = energy();
      perturb(labels[row], delta);
      perturb(labels[col], -delta);

      const double numerical = (ePP - ePM - eMP + eMM) / (4.0 * delta * delta);
      const double analytic = hessian(row, col);
      const double scale = std::max({1.0, std::abs(numerical), std::abs(analytic)});
      EXPECT_NEAR(numerical, analytic, relativeTolerance * scale) << "row=" << row << " col=" << col;
    }
  }
}

TEST(minimization_hessian_ewald, single_ion_bogusz_correction_cancels_strain_derivatives)
{
  ForceField forceField = makeChargedForceField();
  System system = makeSingleIonSystem(forceField);
  const double3 basePosition = system.spanOfMoleculeAtoms()[0].position;
  const double3 baseBoxLengths(25.0, 25.0, 25.0);

  const MinimizationDofLayout layout = buildMinimizationDofLayout(system.moleculeData, system.components);
  ASSERT_EQ(layout.numDofs(), 3u);
  GeneralizedHessian hessian(layout.numDofs(), 1);
  std::vector<AtomDynamics> dynamics(1);
  Interactions::computeEwaldFourierHessian(system, layout, hessian, dynamics);

  // For one ion, Fourier + Bogusz equals the constant self term. Consequently every position,
  // position-strain, and strain derivative vanishes. This catches the Bogusz derivatives without
  // relying on cancellation between unrelated particles.
  constexpr double tolerance = 2e-10;
  EXPECT_NEAR(dynamics[0].gradient.x, 0.0, tolerance);
  EXPECT_NEAR(dynamics[0].gradient.y, 0.0, tolerance);
  EXPECT_NEAR(dynamics[0].gradient.z, 0.0, tolerance);
  for (std::size_t i = 0; i < layout.numDofs(); ++i)
  {
    EXPECT_NEAR(hessian.positionStrain()[i], 0.0, tolerance) << "position-strain dof=" << i;
    for (std::size_t j = 0; j < layout.numDofs(); ++j)
    {
      EXPECT_NEAR(hessian(i, j), 0.0, tolerance) << "row=" << i << " col=" << j;
    }
  }
  EXPECT_NEAR(hessian.strainStrain()[0], 0.0, tolerance);
  const double3x3 &strain = hessian.strainGradient();
  EXPECT_NEAR(strain.ax, 0.0, tolerance);
  EXPECT_NEAR(strain.ay, 0.0, tolerance);
  EXPECT_NEAR(strain.az, 0.0, tolerance);
  EXPECT_NEAR(strain.bx, 0.0, tolerance);
  EXPECT_NEAR(strain.by, 0.0, tolerance);
  EXPECT_NEAR(strain.bz, 0.0, tolerance);
  EXPECT_NEAR(strain.cx, 0.0, tolerance);
  EXPECT_NEAR(strain.cy, 0.0, tolerance);
  EXPECT_NEAR(strain.cz, 0.0, tolerance);

  EwaldEnergyEvaluator evaluator{};
  auto energyAtStrain = [&](double epsilon)
  {
    const double scale = std::exp(epsilon);
    system.spanOfMoleculeAtoms()[0].position = scale * basePosition;
    return evaluator(system, SimulationBox(scale * baseBoxLengths.x, scale * baseBoxLengths.y,
                                           scale * baseBoxLengths.z));
  };

  // A comparatively large step suppresses cancellation noise from the constant self-energy
  // offset; the exact one-ion result is strain-independent for every step size.
  constexpr double strainStep = 1e-2;
  const double eMinus = energyAtStrain(-strainStep);
  const double eZero = energyAtStrain(0.0);
  const double ePlus = energyAtStrain(strainStep);
  EXPECT_NEAR((ePlus - eMinus) / (2.0 * strainStep), 0.0, 1e-8);
  EXPECT_NEAR((ePlus - 2.0 * eZero + eMinus) / (strainStep * strainStep), 0.0, 2e-6);

  system.spanOfMoleculeAtoms()[0].position = basePosition;
}

TEST(minimization_hessian_ewald, flexible_charged_pair_strain_matches_finite_difference)
{
  const double delta = 1e-4;
  const double strainStep = 1e-3;
  const double relativeTolerance = 1e-4;

  ForceField forceField = makeChargedForceField();
  System system = makeFlexiblePairSystem(forceField);

  const MinimizationDofLayout layout = buildMinimizationDofLayout(system.moleculeData, system.components);
  GeneralizedHessian hessian(layout.numDofs(), 1);
  std::vector<AtomDynamics> dynamics(system.spanOfMoleculeAtoms().size());
  Interactions::computeEwaldFourierHessian(system, layout, hessian, dynamics);

  EwaldEnergyEvaluator evaluator{};
  std::span<Atom> atoms = system.spanOfMoleculeAtoms();
  const std::vector<FlexibleDofLabel> labels = buildFlexibleDofLabels(layout, system);
  const double3 boxLengths(25.0, 25.0, 25.0);

  std::vector<double3> savedPositions(atoms.size());
  auto evaluate = [&](std::optional<FlexibleDofLabel> perturb, double positionDelta, double strainExponent) -> double
  {
    for (std::size_t atom = 0; atom < atoms.size(); ++atom) savedPositions[atom] = atoms[atom].position;
    if (perturb) (&atoms[perturb->atom].position.x)[perturb->axis] += positionDelta;
    const double factor = std::exp(strainExponent);
    for (std::size_t atom = 0; atom < atoms.size(); ++atom) atoms[atom].position = atoms[atom].position * factor;
    const SimulationBox strainedBox(factor * boxLengths.x, factor * boxLengths.y, factor * boxLengths.z);
    const double energy = evaluator(system, strainedBox);
    for (std::size_t atom = 0; atom < atoms.size(); ++atom) atoms[atom].position = savedPositions[atom];
    return energy;
  };

  // strain-strain block.
  {
    const double eZero = evaluate(std::nullopt, 0.0, 0.0);
    const double ePlus = evaluate(std::nullopt, 0.0, strainStep);
    const double eMinus = evaluate(std::nullopt, 0.0, -strainStep);
    const double numerical = (ePlus - 2.0 * eZero + eMinus) / (strainStep * strainStep);
    const double analytic = hessian.strainStrain()[0];
    const double scale = std::max({1.0, std::abs(numerical), std::abs(analytic)});
    EXPECT_NEAR(numerical, analytic, relativeTolerance * scale) << "strain-strain";
  }

  // position-strain blocks.
  for (std::size_t dof = 0; dof < labels.size(); ++dof)
  {
    const double ePP = evaluate(labels[dof], delta, strainStep);
    const double ePM = evaluate(labels[dof], delta, -strainStep);
    const double eMP = evaluate(labels[dof], -delta, strainStep);
    const double eMM = evaluate(labels[dof], -delta, -strainStep);
    const double numerical = (ePP - ePM - eMP + eMM) / (4.0 * delta * strainStep);
    const double analytic = hessian.positionStrain()[dof];
    const double scale = std::max({1.0, std::abs(numerical), std::abs(analytic)});
    EXPECT_NEAR(numerical, analytic, relativeTolerance * scale) << "position-strain dof=" << dof;
  }
}

TEST(minimization_hessian_ewald, rigid_co2_pair_matches_finite_difference)
{
  const double delta = 1e-4;
  const double relativeTolerance = 1e-5;

  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  Component co2 = Component::makeCO2(forceField, 0, false);

  System system = System(forceField, SimulationBox(25.0, 25.0, 25.0), false, 300.0, 1e4, 1.0, {}, {co2}, {}, {2}, 5);

  const std::array<RigidState, 2> baseStates = {
      RigidState{double3(9.0, 10.0, 10.0),
                 simd_quatd::fromAxisAngle(0.4, double3(1.0, 2.0, 3.0) / std::sqrt(14.0)).normalized()},
      RigidState{double3(12.5, 10.7, 10.4),
                 simd_quatd::fromAxisAngle(-0.8, double3(0.5, -1.0, 2.0) / std::sqrt(5.25)).normalized()}};

  std::array<double, 12> displacement{};
  auto applyState = [&](double strainExponent)
  {
    const double factor = std::exp(strainExponent);
    for (std::size_t moleculeIndex = 0; moleculeIndex < 2; ++moleculeIndex)
    {
      applyDisplacedState(system, moleculeIndex, baseStates[moleculeIndex],
                          std::span<const double>(displacement).subspan(moleculeIndex * 6, 6), factor);
    }
  };
  applyState(0.0);

  const MinimizationDofLayout layout = buildMinimizationDofLayout(system.moleculeData, system.components);
  ASSERT_EQ(layout.numDofs(), 12u);

  GeneralizedHessian hessian(layout.numDofs(), 1);
  std::vector<AtomDynamics> dynamics(system.spanOfMoleculeAtoms().size());
  Interactions::computeEwaldFourierHessian(system, layout, hessian, dynamics);

  // Global DOF index -> slot in the 12-vector displacement.
  std::array<std::size_t, 12> dofToSlot{};
  for (std::size_t moleculeIndex = 0; moleculeIndex < 2; ++moleculeIndex)
  {
    for (std::size_t local = 0; local < 6; ++local)
    {
      const auto dof = layout.rigidMoleculeDof(moleculeIndex, static_cast<RigidDof>(local));
      ASSERT_TRUE(dof.has_value());
      dofToSlot[*dof] = moleculeIndex * 6 + local;
    }
  }

  EwaldEnergyEvaluator evaluator{};
  const double3 boxLengths(25.0, 25.0, 25.0);
  auto energyAt = [&](double strainExponent)
  {
    applyState(strainExponent);
    const double factor = std::exp(strainExponent);
    const SimulationBox strainedBox(factor * boxLengths.x, factor * boxLengths.y, factor * boxLengths.z);
    const double energy = evaluator(system, strainedBox);
    return energy;
  };

  // Position-position block over all rigid DOFs.
  for (std::size_t row = 0; row < 12; ++row)
  {
    for (std::size_t col = 0; col < 12; ++col)
    {
      const std::size_t rowSlot = dofToSlot[row];
      const std::size_t colSlot = dofToSlot[col];

      displacement.fill(0.0);
      displacement[rowSlot] += delta;
      displacement[colSlot] += delta;
      const double ePP = energyAt(0.0);

      displacement.fill(0.0);
      displacement[rowSlot] += delta;
      displacement[colSlot] -= delta;
      const double ePM = energyAt(0.0);

      displacement.fill(0.0);
      displacement[rowSlot] -= delta;
      displacement[colSlot] += delta;
      const double eMP = energyAt(0.0);

      displacement.fill(0.0);
      displacement[rowSlot] -= delta;
      displacement[colSlot] -= delta;
      const double eMM = energyAt(0.0);

      const double numerical = (ePP - ePM - eMP + eMM) / (4.0 * delta * delta);
      const double analytic = hessian(row, col);
      const double scale = std::max({1.0, std::abs(numerical), std::abs(analytic)});
      EXPECT_NEAR(numerical, analytic, relativeTolerance * scale) << "row=" << row << " col=" << col;
    }
  }

  // Strain blocks.
  const double strainStep = 1e-3;
  const double strainTolerance = 1e-4;
  {
    displacement.fill(0.0);
    const double eZero = energyAt(0.0);
    const double ePlus = energyAt(strainStep);
    const double eMinus = energyAt(-strainStep);
    const double numerical = (ePlus - 2.0 * eZero + eMinus) / (strainStep * strainStep);
    const double analytic = hessian.strainStrain()[0];
    const double scale = std::max({1.0, std::abs(numerical), std::abs(analytic)});
    EXPECT_NEAR(numerical, analytic, strainTolerance * scale) << "strain-strain";
  }

  for (std::size_t dof = 0; dof < 12; ++dof)
  {
    const std::size_t slot = dofToSlot[dof];

    displacement.fill(0.0);
    displacement[slot] += delta;
    const double ePP = energyAt(strainStep);
    const double ePM = energyAt(-strainStep);
    displacement.fill(0.0);
    displacement[slot] -= delta;
    const double eMP = energyAt(strainStep);
    const double eMM = energyAt(-strainStep);

    const double numerical = (ePP - ePM - eMP + eMM) / (4.0 * delta * strainStep);
    const double analytic = hessian.positionStrain()[dof];
    const double scale = std::max({1.0, std::abs(numerical), std::abs(analytic)});
    EXPECT_NEAR(numerical, analytic, strainTolerance * scale) << "position-strain dof=" << dof;
  }

  // Restore base state.
  displacement.fill(0.0);
  applyState(0.0);
}
