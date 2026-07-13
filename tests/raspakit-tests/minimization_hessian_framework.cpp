#include <gtest/gtest.h>

import std;

import int3;
import double3;
import double3x3;
import simd_quatd;
import atom;
import atom_dynamics;
import molecule;
import forcefield;
import framework;
import component;
import system;
import simulationbox;
import running_energy;
import generalized_hessian;
import minimization_dof_layout;
import interactions_hessian_framework_molecule;
import interactions_framework_molecule;
import potential_energy_vdw;
import potential_energy_coulomb;

namespace
{
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
}  // namespace

TEST(minimization_hessian_framework, rigid_co2_in_itq29_vdw_matches_finite_difference)
{
  const double delta = 1e-4;
  const double relativeTolerance = 1e-5;
  // Roundoff floor of the FD reference: the total framework energy is large, so the
  // four-point cancellation leaves noise of order eps * E / (4 delta^2).
  const double absoluteTolerance = 1e-3;

  ForceField forceField = ForceField::makeZeoliteForceField(11.8, true, false, true);
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
  Interactions::computeFrameworkMoleculeHessian(system, layout, hessian);

  auto energyAtDisplacement = [&](std::span<const double> displacement)
  {
    Molecule &molecule = system.moleculeData[0];
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

      displacement.fill(0.0);
      displacement[row] += delta;
      displacement[col] += delta;
      const double ePP = energyAtDisplacement(displacement);

      displacement.fill(0.0);
      displacement[row] += delta;
      displacement[col] -= delta;
      const double ePM = energyAtDisplacement(displacement);

      displacement.fill(0.0);
      displacement[row] -= delta;
      displacement[col] += delta;
      const double eMP = energyAtDisplacement(displacement);

      displacement.fill(0.0);
      displacement[row] -= delta;
      displacement[col] -= delta;
      const double eMM = energyAtDisplacement(displacement);

      const double numerical = (ePP - ePM - eMP + eMM) / (4.0 * delta * delta);
      const double analytic = hessian(*dofRow, *dofCol);
      const double scale = std::max({1.0, std::abs(numerical), std::abs(analytic)});
      EXPECT_NEAR(numerical, analytic, absoluteTolerance + relativeTolerance * scale)
          << "row=" << row << " col=" << col;
    }
  }
}

TEST(minimization_hessian_framework, rigid_co2_in_itq29_strain_fixed_neighbor_list_matches_finite_difference)
{
  const double delta = 1e-4;
  const double strainStep = 3e-4;
  const double relativeTolerance = 1e-4;
  // Roundoff floor: the total framework-molecule energy is large, so the finite-difference
  // cancellation leaves absolute noise.
  const double absoluteTolerance = 5e-3;

  ForceField forceField = ForceField::makeZeoliteForceField(11.8, true, false, true);
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
    Molecule &molecule = system.moleculeData[0];
    molecule.centerOfMassPosition =
        factor * (baseCom + double3(displacement[0], displacement[1], displacement[2]));
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
  const RunningEnergy analyticEnergy = Interactions::computeFrameworkMoleculeHessian(system, layout, hessian);

  struct BasePair
  {
    std::size_t frameworkAtom;
    std::size_t moleculeAtom;
    bool vdw;
    bool charge;
  };
  std::vector<BasePair> basePairs;
  // Differentiate the smooth energy represented by the Hessian: the interaction set at the
  // expansion point. Rebuilding the cutoff list at +/- strain is not a valid FD oracle when a
  // pair crosses the hard cutoff; in this ITQ29 state two VDW pairs do so at epsilon=+3e-4.
  const double vdwCutoffSquared =
      system.forceField.cutOffFrameworkVDW * system.forceField.cutOffFrameworkVDW;
  const double chargeCutoffSquared = system.forceField.cutOffCoulomb * system.forceField.cutOffCoulomb;
  std::span<Atom> moleculeAtoms = system.spanOfMoleculeAtoms();
  for (std::size_t frameworkAtom = 0; frameworkAtom < frameworkAtoms.size(); ++frameworkAtom)
  {
    for (std::size_t moleculeAtom = 0; moleculeAtom < moleculeAtoms.size(); ++moleculeAtom)
    {
      const double3 dr = baseBox.applyPeriodicBoundaryConditions(moleculeAtoms[moleculeAtom].position -
                                                                 frameworkAtoms[frameworkAtom].position);
      const double rr = double3::dot(dr, dr);
      const bool vdw = rr < vdwCutoffSquared;
      const bool charge = system.forceField.useCharge && rr < chargeCutoffSquared;
      if (vdw || charge)
      {
        basePairs.push_back({frameworkAtom, moleculeAtom, vdw, charge});
      }
    }
  }

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

  auto fixedBasePairEnergyAt = [&](double strainExponent)
  {
    applyState(strainExponent);
    const double factor = std::exp(strainExponent);
    const SimulationBox strainedBox(factor * baseBox.cell.ax, factor * baseBox.cell.by, factor * baseBox.cell.cz);
    double energy = 0.0;
    for (const BasePair &pair : basePairs)
    {
      const Atom &frameworkAtom = frameworkAtoms[pair.frameworkAtom];
      const Atom &moleculeAtom = moleculeAtoms[pair.moleculeAtom];
      const double3 dr =
          strainedBox.applyPeriodicBoundaryConditions(moleculeAtom.position - frameworkAtom.position);
      const double rr = double3::dot(dr, dr);
      if (pair.vdw)
      {
        energy += Potentials::potentialVDWEnergy(
                      system.forceField, frameworkAtom.scalingVDW, moleculeAtom.scalingVDW, rr,
                      static_cast<std::size_t>(frameworkAtom.type), static_cast<std::size_t>(moleculeAtom.type))
                      .energy;
      }
      if (pair.charge)
      {
        energy += Potentials::potentialCoulombEnergy(
                      system.forceField, frameworkAtom.scalingCoulomb, moleculeAtom.scalingCoulomb, std::sqrt(rr),
                      frameworkAtom.charge, moleculeAtom.charge)
                      .energy;
      }
    }
    return energy;
  };

  // Consistency: the Hessian routine and the FD evaluator must see the same base state.
  EXPECT_NEAR(analyticEnergy.frameworkMoleculeVDW + analyticEnergy.frameworkMoleculeCharge, energyAt(0.0), 1e-8)
      << "base-state energy mismatch between analytic Hessian and FD evaluator";

  // strain-strain block.
  {
    displacement.fill(0.0);
    const double eZero = fixedBasePairEnergyAt(0.0);
    const double ePlus = fixedBasePairEnergyAt(strainStep);
    const double eMinus = fixedBasePairEnergyAt(-strainStep);
    const double numerical = (ePlus - 2.0 * eZero + eMinus) / (strainStep * strainStep);
    const double analytic = hessian.strainStrain()[0];
    const double scale = std::max({1.0, std::abs(numerical), std::abs(analytic)});
    EXPECT_NEAR(numerical, analytic, absoluteTolerance + relativeTolerance * scale) << "strain-strain";
  }

  // position-strain blocks over the rigid COM and orientation DOFs.
  for (std::size_t slot = 0; slot < 6; ++slot)
  {
    const auto dof = layout.rigidMoleculeDof(0, static_cast<RigidDof>(slot));
    ASSERT_TRUE(dof);

    displacement.fill(0.0);
    displacement[slot] += delta;
    const double ePP = fixedBasePairEnergyAt(strainStep);
    const double ePM = fixedBasePairEnergyAt(-strainStep);
    displacement.fill(0.0);
    displacement[slot] -= delta;
    const double eMP = fixedBasePairEnergyAt(strainStep);
    const double eMM = fixedBasePairEnergyAt(-strainStep);

    const double numerical = (ePP - ePM - eMP + eMM) / (4.0 * delta * strainStep);
    const double analytic = hessian.positionStrain()[*dof];
    const double scale = std::max({1.0, std::abs(numerical), std::abs(analytic)});
    EXPECT_NEAR(numerical, analytic, absoluteTolerance + relativeTolerance * scale)
        << "position-strain dof=" << *dof;
  }

  // Restore base state.
  displacement.fill(0.0);
  applyState(0.0);
}

TEST(minimization_hessian_framework, single_framework_atom_rigid_co2_strain_step_convergence)
{
  ForceField forceField = ForceField::makeZeoliteForceField(11.8, true, false, true);
  const Atom frameworkSite(double3(0.25, 0.25, 0.25), 0.0, 1.0, 0, 3, 0, false, false);
  Framework framework(forceField, "single-site", SimulationBox(40.0, 40.0, 40.0), 1, {frameworkSite},
                      {frameworkSite}, int3(1, 1, 1));
  Component co2 = Component::makeCO2(forceField, 0, false);
  System system =
      System(forceField, std::nullopt, false, 300.0, 1e4, 1.0, {framework}, {co2}, {}, {1}, 5);

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
    Molecule &molecule = system.moleculeData[0];
    molecule.centerOfMassPosition =
        scale * (comBase + double3(displacement[0], displacement[1], displacement[2]));
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
  Interactions::computeFrameworkMoleculeHessian(system, layout, hessian, dynamics);

  constexpr double gradientStep = 1e-5;
  std::span<Atom> moleculeAtoms = system.spanOfMoleculeAtoms();
  for (std::size_t atom = 0; atom < moleculeAtoms.size(); ++atom)
  {
    for (std::size_t axis = 0; axis < 3; ++axis)
    {
      double &coordinate = (&moleculeAtoms[atom].position.x)[axis];
      coordinate += gradientStep;
      const RunningEnergy plus = Interactions::computeFrameworkMoleculeEnergy(
          system.forceField, SimulationBox(40.0, 40.0, 40.0), system.interpolationGrids, system.framework,
          system.spanOfFrameworkAtoms(), moleculeAtoms);
      coordinate -= 2.0 * gradientStep;
      const RunningEnergy minus = Interactions::computeFrameworkMoleculeEnergy(
          system.forceField, SimulationBox(40.0, 40.0, 40.0), system.interpolationGrids, system.framework,
          system.spanOfFrameworkAtoms(), moleculeAtoms);
      coordinate += gradientStep;
      const double numerical =
          ((plus.frameworkMoleculeVDW + plus.frameworkMoleculeCharge) -
           (minus.frameworkMoleculeVDW + minus.frameworkMoleculeCharge)) /
          (2.0 * gradientStep);
      EXPECT_NEAR((&dynamics[atom].gradient.x)[axis], numerical,
                  1e-5 * std::max(1.0, std::abs(numerical)));
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
  constexpr std::array<double, 4> strainSteps = {2e-3, 1e-3, 5e-4, 2.5e-4};
  constexpr double positionStep = 2e-5;
  for (double strainStep : strainSteps)
  {
    displacement.fill(0.0);
    const double numericalStrain =
        (energyAt(strainStep) - 2.0 * energyAt(0.0) + energyAt(-strainStep)) / (strainStep * strainStep);
    EXPECT_NEAR(numericalStrain, hessian.strainStrain()[0],
                2e-4 * std::max({1.0, std::abs(numericalStrain), std::abs(hessian.strainStrain()[0])}))
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
      EXPECT_NEAR(numerical, analytic, 2e-4 * std::max({1.0, std::abs(numerical), std::abs(analytic)}))
          << "slot=" << slot << " strainStep=" << strainStep;
    }
  }
}
