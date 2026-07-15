#include <gtest/gtest.h>

import std;

import atom;
import bend_potential;
import bond_potential;
import component;
import connectivity_table;
import double3;
import double3x3;
import forcefield;
import framework;
import generalized_hessian;
import interactions_polarization_derivatives;
import intra_molecular_potentials;
import int3;
import minimization_cell_layout;
import minimization_dof_layout;
import minimization_evaluate_derivatives;
import minimization_generalized_coordinates;
import minimization;
import minimization_options;
import molecule;
import simulationbox;
import simd_quatd;
import system;
import torsion_potential;
import vdwparameters;

namespace
{
struct Evaluation
{
  double energy{};
  std::vector<double> gradient;
  std::vector<double> hessian;
};

Evaluation evaluate(System& system, const MinimizationDofLayout& layout)
{
  GeneralizedHessian hessian(layout.numDofs(), 0);
  std::vector<double> gradient(layout.numDofs());
  DerivativeCapabilities capabilities{.energy = true, .gradient = true, .hessianPositionPosition = true};
  DerivativeResults results{.gradient = gradient, .hessian = hessian};
  evaluateDerivatives(system, layout, capabilities, results);
  return Evaluation{results.energy, gradient,
                    std::vector<double>(hessian.positionPosition().begin(), hessian.positionPosition().end())};
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

System makeLennardJonesPair(CellMinimizationType type)
{
  ForceField forceField = ForceField::makeZeoliteForceField(10.0, true, false, true);
  useSecondOrderTaylorShiftedLennardJones(forceField);
  Component methane = Component::makeMethane(forceField, 0);
  methane.rigid = false;
  System system =
      System(forceField, SimulationBox(24.0, 22.0, 20.0), false, 300.0, 2.0e5, 1.0, {}, {methane}, {}, {2}, 5);
  system.spanOfMoleculeAtoms()[0].position = double3(8.1, 9.3, 10.2);
  system.spanOfMoleculeAtoms()[1].position = double3(11.8, 10.7, 9.6);
  system.cellMinimizationType = type;
  return system;
}

ForceField makeChargedForceField()
{
  return ForceField({{"P", false, 15.0, 0.8, 0.0, 8, false}, {"N", false, 14.0, -0.8, 0.0, 8, false}},
                    {{60.0, 3.0}, {40.0, 3.2}}, ForceField::MixingRule::Lorentz_Berthelot, 9.0, 9.0, 9.0, true, false,
                    true);
}

System makeChargedPair()
{
  ForceField forceField = makeChargedForceField();
  Component positive =
      Component(forceField, "positive", 100.0, 1e6, 0.2, {Atom({0.0, 0.0, 0.0}, 0.8, 1.0, 0, 0, 0, false, false)},
                ConnectivityTable(1), Potentials::IntraMolecularPotentials{}, 5, 21);
  Component negative =
      Component(forceField, "negative", 100.0, 1e6, 0.2, {Atom({0.0, 0.0, 0.0}, -0.8, 1.0, 0, 1, 0, false, false)},
                ConnectivityTable(1), Potentials::IntraMolecularPotentials{}, 5, 21);
  positive.rigid = false;
  negative.rigid = false;
  System system = System(forceField, SimulationBox(23.0, 21.0, 19.0), false, 300.0, 1.5e5, 1.0, {},
                         {positive, negative}, {}, {1, 1}, 5);
  system.spanOfMoleculeAtoms()[0].position = double3(7.3, 9.1, 8.4);
  system.spanOfMoleculeAtoms()[1].position = double3(12.2, 11.4, 10.7);
  system.cellMinimizationType = CellMinimizationType::Regular;
  return system;
}

void regenerateRigidMolecule(System& system, std::size_t moleculeIndex)
{
  Molecule& molecule = system.moleculeData[moleculeIndex];
  const Component& component = system.components[molecule.componentId];
  const double3x3 rotation = double3x3::buildRotationMatrixInverse(molecule.orientation);
  for (std::size_t atom = 0; atom < molecule.numberOfAtoms; ++atom)
  {
    system.spanOfMoleculeAtoms()[molecule.atomIndex + atom].position =
        molecule.centerOfMassPosition + rotation * component.atoms[atom].position;
  }
}

System makeRigidChargedPair()
{
  ForceField forceField = makeChargedForceField();
  Component molecule = Component(
      forceField, "rigid-neutral", 100.0, 1e6, 0.2,
      {Atom({-1.0, 0.0, 0.0}, -0.4, 1.0, 0, 1, 0, false, false), Atom({0.0, 0.0, 0.0}, 0.8, 1.0, 0, 0, 0, false, false),
       Atom({1.0, 0.0, 0.0}, -0.4, 1.0, 0, 1, 0, false, false)},
      ConnectivityTable(3), Potentials::IntraMolecularPotentials{}, 5, 21);
  molecule.rigid = true;
  System system =
      System(forceField, SimulationBox(23.0, 21.0, 19.0), false, 300.0, 1.5e5, 1.0, {}, {molecule}, {}, {2}, 5);
  system.moleculeData[0].centerOfMassPosition = double3(7.3, 9.1, 8.4);
  system.moleculeData[1].centerOfMassPosition = double3(12.2, 11.4, 10.7);
  system.moleculeData[0].orientation = simd_quatd::fromAxisAngle(0.35, double3(0.0, 1.0, 0.0));
  system.moleculeData[1].orientation = simd_quatd::fromAxisAngle(-0.27, double3(0.0, 0.0, 1.0));
  regenerateRigidMolecule(system, 0);
  regenerateRigidMolecule(system, 1);
  system.cellMinimizationType = CellMinimizationType::Regular;
  return system;
}

// Two flexible, polarizable single-atom ions in a box (no framework): the polarization field is entirely
// real-space, so this validates the analytic real-space cell-strain derivatives (stress, position-strain,
// strain-strain), including the many-body Term A coupling.
System makePolarizableFlexiblePair()
{
  ForceField forceField({{"P", false, 15.0, 0.8, 1.5, 8, false}, {"N", false, 14.0, -0.8, 1.2, 8, false}},
                        {{60.0, 3.0}, {40.0, 3.2}}, ForceField::MixingRule::Lorentz_Berthelot, 9.0, 9.0, 9.0, true, false,
                        true);
  forceField.computePolarization = true;
  forceField.omitInterPolarization = false;
  forceField.omitInterInteractions = false;
  Component positive =
      Component(forceField, "positive", 100.0, 1e6, 0.2, {Atom({0.0, 0.0, 0.0}, 0.8, 1.0, 0, 0, 0, false, false)},
                ConnectivityTable(1), Potentials::IntraMolecularPotentials{}, 5, 21);
  Component negative =
      Component(forceField, "negative", 100.0, 1e6, 0.2, {Atom({0.0, 0.0, 0.0}, -0.8, 1.0, 0, 1, 0, false, false)},
                ConnectivityTable(1), Potentials::IntraMolecularPotentials{}, 5, 21);
  positive.rigid = false;
  negative.rigid = false;
  System system = System(forceField, SimulationBox(23.0, 21.0, 19.0), false, 300.0, 1.5e5, 1.0, {},
                         {positive, negative}, {}, {1, 1}, 5);
  system.spanOfMoleculeAtoms()[0].position = double3(7.3, 9.1, 8.4);
  system.spanOfMoleculeAtoms()[1].position = double3(12.2, 11.4, 10.7);
  system.cellMinimizationType = CellMinimizationType::Regular;
  return system;
}

// A rigid charged framework plus two flexible, polarizable ions with Ewald: the polarization field now has
// a fixed-framework reciprocal (structure-factor) contribution in addition to the real-space part, so this
// validates the analytic reciprocal-space (k / volume) cell-strain derivatives.
System makePolarizableRigidFrameworkPair()
{
  ForceField forceField({{"P", false, 15.0, 0.4, 1.5, 8, false}, {"N", false, 14.0, -0.4, 1.2, 8, false}},
                        {{60.0, 3.0}, {40.0, 3.2}}, ForceField::MixingRule::Lorentz_Berthelot, 11.0, 11.0, 11.0, true,
                        false, true);
  forceField.computePolarization = true;
  forceField.omitInterPolarization = false;
  forceField.omitInterInteractions = false;
  forceField.automaticEwald = false;
  forceField.EwaldAlpha = 0.3;
  forceField.numberOfWaveVectors = int3(8, 8, 8);

  const SimulationBox box(24.0, 24.0, 24.0);
  const std::size_t typeP = *forceField.findPseudoAtom("P");
  const std::size_t typeN = *forceField.findPseudoAtom("N");

  std::vector<Atom> frameworkAtoms{Atom({0.25, 0.5, 0.5}, 0.4, 1.0, 0, static_cast<std::uint16_t>(typeP), 0, 0, true),
                                   Atom({0.75, 0.5, 0.5}, -0.4, 1.0, 0, static_cast<std::uint16_t>(typeN), 0, 0, true)};
  Framework framework(forceField, "rigid-charged", box, 1, frameworkAtoms, frameworkAtoms, int3(1, 1, 1));
  framework.rigid = true;

  Component cation(forceField, "cation", 30.0, 1.0e6, 0.011,
                   {Atom({0.0, 0.0, 0.0}, 0.4, 1.0, 0, static_cast<std::uint16_t>(typeP), 0, false, false)},
                   ConnectivityTable(1), Potentials::IntraMolecularPotentials{}, 0, 0);
  cation.rigid = false;
  Component anion(forceField, "anion", 30.0, 1.0e6, 0.011,
                  {Atom({0.0, 0.0, 0.0}, -0.4, 1.0, 0, static_cast<std::uint16_t>(typeN), 0, false, false)},
                  ConnectivityTable(1), Potentials::IntraMolecularPotentials{}, 0, 0);
  anion.rigid = false;

  System system(forceField, box, false, 300.0, 1e4, 1.0, {framework}, {cation, anion}, {}, {1, 1}, 5);
  std::span<Atom> atoms = system.spanOfMoleculeAtoms();
  atoms[0].position = double3(6.0, 14.0, 12.0);
  atoms[1].position = double3(9.0, 11.5, 13.5);
  system.cellMinimizationType = CellMinimizationType::Regular;
  return system;
}

// A rigid, polarizable three-site molecule pair in a box (no framework): the polarizable field points sit
// off the center of mass, so this validates the analytic real-space cell-strain derivatives for rigid
// molecules (center-of-mass separation arm; pure force change projected onto center-of-mass and orientation
// DOFs plus the B_a . grad kinematic term).
System makePolarizableRigidMoleculePair()
{
  ForceField forceField({{"P", false, 15.0, 0.8, 1.5, 8, false}, {"N", false, 14.0, -0.4, 1.2, 8, false}},
                        {{60.0, 3.0}, {40.0, 3.2}}, ForceField::MixingRule::Lorentz_Berthelot, 9.0, 9.0, 9.0, true, false,
                        true);
  forceField.computePolarization = true;
  forceField.omitInterPolarization = false;
  forceField.omitInterInteractions = false;
  const std::size_t typeP = *forceField.findPseudoAtom("P");
  const std::size_t typeN = *forceField.findPseudoAtom("N");
  Component molecule =
      Component(forceField, "rigid-polar", 100.0, 1e6, 0.2,
                {Atom({-1.0, 0.0, 0.0}, -0.4, 1.0, 0, static_cast<std::uint16_t>(typeN), 0, false, false),
                 Atom({0.0, 0.0, 0.0}, 0.8, 1.0, 0, static_cast<std::uint16_t>(typeP), 0, false, false),
                 Atom({1.0, 0.0, 0.0}, -0.4, 1.0, 0, static_cast<std::uint16_t>(typeN), 0, false, false)},
                ConnectivityTable(3), Potentials::IntraMolecularPotentials{}, 5, 21);
  molecule.rigid = true;
  System system =
      System(forceField, SimulationBox(23.0, 21.0, 19.0), false, 300.0, 1.5e5, 1.0, {}, {molecule}, {}, {2}, 5);
  system.moleculeData[0].centerOfMassPosition = double3(7.3, 9.1, 8.4);
  system.moleculeData[1].centerOfMassPosition = double3(12.2, 11.4, 10.7);
  system.moleculeData[0].orientation = simd_quatd::fromAxisAngle(0.35, double3(0.0, 1.0, 0.0));
  system.moleculeData[1].orientation = simd_quatd::fromAxisAngle(-0.27, double3(0.0, 0.0, 1.0));
  regenerateRigidMolecule(system, 0);
  regenerateRigidMolecule(system, 1);
  system.cellMinimizationType = CellMinimizationType::Regular;
  return system;
}

// A rigid, polarizable three-site molecule inside a rigid charged framework with Ewald: the polarizable
// field points belong to a rigid molecule, so this validates the reciprocal-space cell-strain derivatives
// including the body-offset phase terms of a rigid field point (the phase k.r_A deforms only via the COM).
// The molecule atoms are left uncharged (only polarizable): the framework reciprocal field is still sampled
// at the rigid molecule atoms, but the (separately validated) framework-molecule Ewald charge interaction is
// switched off so this test isolates the polarization reciprocal contribution.
System makePolarizableRigidMoleculeInFramework()
{
  ForceField forceField({{"P", false, 15.0, 0.4, 1.5, 8, false}, {"N", false, 14.0, -0.4, 1.2, 8, false}},
                        {{60.0, 3.0}, {40.0, 3.2}}, ForceField::MixingRule::Lorentz_Berthelot, 11.0, 11.0, 11.0, true,
                        false, true);
  forceField.computePolarization = true;
  forceField.omitInterPolarization = false;
  forceField.omitInterInteractions = false;
  forceField.automaticEwald = false;
  forceField.EwaldAlpha = 0.3;
  forceField.numberOfWaveVectors = int3(8, 8, 8);

  const SimulationBox box(24.0, 24.0, 24.0);
  const std::size_t typeP = *forceField.findPseudoAtom("P");
  const std::size_t typeN = *forceField.findPseudoAtom("N");

  std::vector<Atom> frameworkAtoms{Atom({0.25, 0.5, 0.5}, 0.4, 1.0, 0, static_cast<std::uint16_t>(typeP), 0, 0, true),
                                   Atom({0.75, 0.5, 0.5}, -0.4, 1.0, 0, static_cast<std::uint16_t>(typeN), 0, 0, true)};
  Framework framework(forceField, "rigid-charged", box, 1, frameworkAtoms, frameworkAtoms, int3(1, 1, 1));
  framework.rigid = true;

  Component molecule =
      Component(forceField, "rigid-polar", 30.0, 1.0e6, 0.011,
                {Atom({-1.0, 0.0, 0.0}, 0.0, 1.0, 0, static_cast<std::uint16_t>(typeN), 0, false, false),
                 Atom({0.0, 0.0, 0.0}, 0.0, 1.0, 0, static_cast<std::uint16_t>(typeP), 0, false, false),
                 Atom({1.0, 0.0, 0.0}, 0.0, 1.0, 0, static_cast<std::uint16_t>(typeN), 0, false, false)},
                ConnectivityTable(3), Potentials::IntraMolecularPotentials{}, 0, 0);
  molecule.rigid = true;

  System system(forceField, box, false, 300.0, 1e4, 1.0, {framework}, {molecule}, {}, {1}, 5);
  system.moleculeData[0].centerOfMassPosition = double3(12.0, 6.0, 12.0);
  system.moleculeData[0].orientation = simd_quatd::fromAxisAngle(0.5, double3(0.57735, 0.57735, 0.57735));
  regenerateRigidMolecule(system, 0);
  system.cellMinimizationType = CellMinimizationType::Regular;
  return system;
}

// Finite-difference check of the polarization cell-strain derivatives in isolation: it differences only the
// polarization energy of computePolarizationDerivatives under a homogeneous cell deformation, so it validates
// the analytic polarization stress (cellGradient) and strain-strain block directly, without depending on the
// other interactions' cell-strain Hessians (the framework-molecule Ewald/vdW rigid-body strain Hessian is a
// separate, pre-existing code path).
void expectPolarizationStrainMatchesFiniteDifference(System system, double gradientTolerance,
                                                     double hessianTolerance)
{
  const CellMinimizationLayout cellLayout =
      makeCellMinimizationLayout(system.cellMinimizationType, system.monoclinicAngleType);
  const std::size_t numberOfFrameworkAtoms = system.spanOfFrameworkAtoms().size();
  const std::size_t numberOfMoleculeAtoms = system.spanOfMoleculeAtoms().size();
  const bool flexibleFramework = system.framework.has_value() && !system.framework->rigid;
  const std::size_t frameworkDofs = flexibleFramework ? numberOfFrameworkAtoms : 0;
  const MinimizationDofLayout layout =
      buildMinimizationDofLayout(system.moleculeData, system.components, frameworkDofs, cellLayout.size());

  std::vector<std::uint8_t> movable(numberOfFrameworkAtoms + numberOfMoleculeAtoms, 0);
  for (std::size_t atom = 0; atom < numberOfFrameworkAtoms; ++atom) movable[atom] = flexibleFramework ? 1 : 0;
  for (std::size_t atom = 0; atom < numberOfMoleculeAtoms; ++atom) movable[numberOfFrameworkAtoms + atom] = 1;

  const Interactions::PolarizationDerivatives reference =
      Interactions::computePolarizationDerivatives(system, movable, cellLayout.bases);
  constexpr double delta = 2.0e-5;

  const auto polarizationEnergy = [&](const std::vector<double>& displacement)
  {
    System deformed = system;
    applyGeneralizedDisplacement(deformed, layout, displacement);
    return Interactions::computePolarizationDerivatives(deformed, movable).energy;
  };
  // Evaluate the reference along the same deformation path (zero displacement) so the rigid-molecule
  // regeneration in applyGeneralizedDisplacement does not introduce a constant offset in the 2-point
  // diagonal second difference.
  const double referenceEnergy = polarizationEnergy(std::vector<double>(layout.numDofs(), 0.0));

  const std::size_t n = cellLayout.size();
  for (std::size_t a = 0; a < n; ++a)
  {
    const std::size_t columnA = *layout.cellDof(a);
    std::vector<double> plus(layout.numDofs(), 0.0);
    std::vector<double> minus(layout.numDofs(), 0.0);
    plus[columnA] = delta;
    minus[columnA] = -delta;
    const double ePlus = polarizationEnergy(plus);
    const double eMinus = polarizationEnergy(minus);

    const double numericalGradient = (ePlus - eMinus) / (2.0 * delta);
    const double gradientScale = std::max({1.0, std::abs(numericalGradient), std::abs(reference.cellGradient[a])});
    EXPECT_NEAR(reference.cellGradient[a], numericalGradient, gradientTolerance * gradientScale) << "strain a=" << a;

    for (std::size_t b = 0; b < n; ++b)
    {
      double numericalHessian{};
      if (a == b)
      {
        numericalHessian = (ePlus - 2.0 * referenceEnergy + eMinus) / (delta * delta);
      }
      else
      {
        const std::size_t columnB = *layout.cellDof(b);
        std::vector<double> pp(layout.numDofs(), 0.0);
        std::vector<double> pm(layout.numDofs(), 0.0);
        std::vector<double> mp(layout.numDofs(), 0.0);
        std::vector<double> mm(layout.numDofs(), 0.0);
        pp[columnA] = pp[columnB] = delta;
        pm[columnA] = delta;
        pm[columnB] = -delta;
        mp[columnA] = -delta;
        mp[columnB] = delta;
        mm[columnA] = mm[columnB] = -delta;
        numericalHessian = (polarizationEnergy(pp) - polarizationEnergy(pm) - polarizationEnergy(mp) +
                            polarizationEnergy(mm)) /
                           (4.0 * delta * delta);
      }
      const double analyticHessian = reference.strainStrain[a * n + b];
      const double hessianScale = std::max({1.0, std::abs(numericalHessian), std::abs(analyticHessian)});
      EXPECT_NEAR(analyticHessian, numericalHessian, hessianTolerance * hessianScale)
          << "strain a=" << a << " b=" << b;
    }
  }
}

void expectCellDerivativesMatchFiniteDifference(System system, double gradientTolerance, double hessianTolerance)
{
  const CellMinimizationLayout cellLayout =
      makeCellMinimizationLayout(system.cellMinimizationType, system.monoclinicAngleType);
  const std::size_t frameworkDofs =
      system.framework && !system.framework->rigid ? system.spanOfFrameworkAtoms().size() : 0;
  const MinimizationDofLayout layout =
      buildMinimizationDofLayout(system.moleculeData, system.components, frameworkDofs, cellLayout.size());
  const Evaluation reference = evaluate(system, layout);
  constexpr double delta = 2.0e-5;

  for (std::size_t cellCoordinate = 0; cellCoordinate < cellLayout.size(); ++cellCoordinate)
  {
    const std::size_t column = *layout.cellDof(cellCoordinate);
    std::vector<double> plusDisplacement(layout.numDofs(), 0.0);
    std::vector<double> minusDisplacement(layout.numDofs(), 0.0);
    plusDisplacement[column] = delta;
    minusDisplacement[column] = -delta;
    System plusSystem = system;
    System minusSystem = system;
    applyGeneralizedDisplacement(plusSystem, layout, plusDisplacement);
    applyGeneralizedDisplacement(minusSystem, layout, minusDisplacement);
    const Evaluation plus = evaluate(plusSystem, layout);
    const Evaluation minus = evaluate(minusSystem, layout);

    const double numericalGradient = (plus.energy - minus.energy) / (2.0 * delta);
    const double gradientScale = std::max({1.0, std::abs(numericalGradient), std::abs(reference.gradient[column])});
    EXPECT_NEAR(reference.gradient[column], numericalGradient, gradientTolerance * gradientScale)
        << "cell coordinate=" << cellCoordinate;

    for (std::size_t row = 0; row < layout.numDofs(); ++row)
    {
      double numericalHessian{};
      if (row == column)
      {
        numericalHessian = (plus.energy - 2.0 * reference.energy + minus.energy) / (delta * delta);
      }
      else
      {
        std::vector<double> pp(layout.numDofs(), 0.0);
        std::vector<double> pm(layout.numDofs(), 0.0);
        std::vector<double> mp(layout.numDofs(), 0.0);
        std::vector<double> mm(layout.numDofs(), 0.0);
        pp[row] = pp[column] = delta;
        pm[row] = delta;
        pm[column] = -delta;
        mp[row] = -delta;
        mp[column] = delta;
        mm[row] = mm[column] = -delta;
        System ppSystem = system;
        System pmSystem = system;
        System mpSystem = system;
        System mmSystem = system;
        applyGeneralizedDisplacement(ppSystem, layout, pp);
        applyGeneralizedDisplacement(pmSystem, layout, pm);
        applyGeneralizedDisplacement(mpSystem, layout, mp);
        applyGeneralizedDisplacement(mmSystem, layout, mm);
        numericalHessian = (evaluate(ppSystem, layout).energy - evaluate(pmSystem, layout).energy -
                            evaluate(mpSystem, layout).energy + evaluate(mmSystem, layout).energy) /
                           (4.0 * delta * delta);
      }
      const double analyticHessian = reference.hessian[row * layout.numDofs() + column];
      const double hessianScale = std::max({1.0, std::abs(numericalHessian), std::abs(analyticHessian)});
      EXPECT_NEAR(analyticHessian, numericalHessian, hessianTolerance * hessianScale)
          << "row=" << row << " cell coordinate=" << cellCoordinate;
    }
  }
}
}  // namespace

TEST(minimization_variable_cell, symmetric_logarithmic_strain_layouts_and_aliases)
{
  EXPECT_EQ(makeCellMinimizationLayout(CellMinimizationType::Fixed, MonoclinicAngleType::Beta).size(), 0u);
  EXPECT_EQ(makeCellMinimizationLayout(CellMinimizationType::Isotropic, MonoclinicAngleType::Beta).size(), 1u);
  EXPECT_EQ(makeCellMinimizationLayout(CellMinimizationType::Anisotropic, MonoclinicAngleType::Beta).size(), 3u);
  EXPECT_EQ(makeCellMinimizationLayout(CellMinimizationType::Monoclinic, MonoclinicAngleType::Beta).size(), 4u);
  EXPECT_EQ(makeCellMinimizationLayout(CellMinimizationType::Regular, MonoclinicAngleType::Beta).size(), 6u);

  const CellMinimizationLayout regularAlias =
      makeCellMinimizationLayout(CellMinimizationType::RegularUpperTriangle, MonoclinicAngleType::Beta);
  const CellMinimizationLayout monoclinicAlias =
      makeCellMinimizationLayout(CellMinimizationType::MonoclinicUpperTriangle, MonoclinicAngleType::Gamma);
  EXPECT_TRUE(regularAlias.isCompatibilityAlias());
  EXPECT_TRUE(monoclinicAlias.isCompatibilityAlias());
  EXPECT_EQ(regularAlias.size(), 6u);
  EXPECT_EQ(monoclinicAlias.size(), 4u);
  EXPECT_DOUBLE_EQ(monoclinicAlias.bases.back().m12, 0.5);
  EXPECT_DOUBLE_EQ(monoclinicAlias.bases.back().m21, 0.5);
}

TEST(minimization_variable_cell, logarithmic_strain_preserves_positive_volume_and_rigid_geometry)
{
  ForceField forceField = ForceField::makeZeoliteForceField(10.0, true, false, true);
  Component co2 = Component::makeCO2(forceField, 0, false);
  System system = System(forceField, SimulationBox(20.0, 21.0, 22.0), false, 300.0, 0.0, 1.0, {}, {co2}, {}, {1}, 5);
  system.cellMinimizationType = CellMinimizationType::Regular;
  system.moleculeData[0].centerOfMassPosition = double3(8.0, 9.0, 10.0);
  const double initialBond =
      (system.spanOfMoleculeAtoms()[0].position - system.spanOfMoleculeAtoms()[1].position).length();
  const CellMinimizationLayout cellLayout =
      makeCellMinimizationLayout(system.cellMinimizationType, system.monoclinicAngleType);
  const MinimizationDofLayout layout =
      buildMinimizationDofLayout(system.moleculeData, system.components, 0, cellLayout.size());
  std::vector<double> displacement(layout.numDofs(), 0.0);
  displacement[*layout.cellDof(0)] = -0.2;
  displacement[*layout.cellDof(1)] = 0.08;
  displacement[*layout.cellDof(5)] = 0.1;
  applyGeneralizedDisplacement(system, layout, displacement);

  EXPECT_GT(system.simulationBox.volume, 0.0);
  const double finalBond =
      (system.spanOfMoleculeAtoms()[0].position - system.spanOfMoleculeAtoms()[1].position).length();
  EXPECT_NEAR(finalBond, initialBond, 1.0e-12);
}

TEST(minimization_variable_cell, all_cell_modes_vdw_gradient_and_hessian_match_finite_difference)
{
  for (CellMinimizationType type :
       {CellMinimizationType::Isotropic, CellMinimizationType::Anisotropic, CellMinimizationType::Monoclinic,
        CellMinimizationType::Regular, CellMinimizationType::RegularUpperTriangle,
        CellMinimizationType::MonoclinicUpperTriangle})
  {
    System system = makeLennardJonesPair(type);
    system.monoclinicAngleType = MonoclinicAngleType::Gamma;
    expectCellDerivativesMatchFiniteDifference(std::move(system), 2.0e-5, 2.0e-4);
  }
}

TEST(minimization_variable_cell, charged_regular_gradient_and_hessian_match_finite_difference)
{
  expectCellDerivativesMatchFiniteDifference(makeChargedPair(), 3.0e-4, 2.0e-3);
}

TEST(minimization_variable_cell, periodic_image_blocks_match_finite_difference)
{
  System system = makeLennardJonesPair(CellMinimizationType::Regular);
  system.spanOfMoleculeAtoms()[0].position = double3(0.7, 10.0, 10.0);
  system.spanOfMoleculeAtoms()[1].position = double3(23.0, 10.2, 9.8);
  expectCellDerivativesMatchFiniteDifference(std::move(system), 3.0e-5, 6.0e-4);
}

TEST(minimization_variable_cell, periodic_framework_molecule_blocks_match_finite_difference)
{
  ForceField forceField = ForceField::makeZeoliteForceField(10.0, true, false, true);
  useSecondOrderTaylorShiftedLennardJones(forceField);
  Component methane = Component::makeMethane(forceField, 0);
  methane.rigid = false;
  const SimulationBox box(24.0, 22.0, 20.0);
  const Atom frameworkAtom({0.02, 0.45, 0.5}, 0.0, 1.0, 0, 0, 0, 0, true);
  Framework framework(forceField, "periodic-framework", box, 1, {frameworkAtom}, {frameworkAtom}, int3(1, 1, 1));
  framework.rigid = true;
  System system(forceField, box, false, 300.0, 2.0e5, 1.0, framework, {methane}, {}, {1}, 5);
  system.spanOfFrameworkAtoms()[0].position = double3(0.5, 10.0, 10.0);
  system.spanOfMoleculeAtoms()[0].position = double3(23.0, 10.2, 9.8);
  system.cellMinimizationType = CellMinimizationType::Regular;
  expectCellDerivativesMatchFiniteDifference(std::move(system), 3.0e-5, 8.0e-4);
}

TEST(minimization_variable_cell, periodic_framework_bond_blocks_match_finite_difference)
{
  ForceField forceField = ForceField::makeZeoliteForceField(10.0, true, false, true);
  const SimulationBox box(10.0, 10.0, 10.0);
  const Atom atomA({0.05, 0.5, 0.5}, 0.0, 1.0, 0, 0, 0, 0, true);
  const Atom atomB({0.95, 0.5, 0.5}, 0.0, 1.0, 0, 0, 0, 0, true);
  Framework framework(forceField, "periodic-bond", box, 1, {atomA, atomB}, {atomA, atomB}, int3(1, 1, 1));
  framework.rigid = false;
  framework.intraMolecularPotentials.bonds = {BondPotential({0, 1}, BondType::Harmonic, {2000.0, 1.2})};
  framework.intraMolecularImageShifts.bonds = {{{int3(0, 0, 0), int3(-1, 0, 0)}}};
  System system(forceField, box, false, 300.0, 0.0, 1.0, framework, {}, {}, {}, 5);
  system.spanOfFrameworkAtoms()[0].position = double3(0.5, 5.0, 5.0);
  system.spanOfFrameworkAtoms()[1].position = double3(9.5, 5.0, 5.0);
  system.cellMinimizationType = CellMinimizationType::Regular;
  expectCellDerivativesMatchFiniteDifference(std::move(system), 1.0e-6, 1.0e-3);
}

TEST(minimization_variable_cell, periodic_framework_bend_and_torsion_blocks_match_finite_difference)
{
  ForceField forceField = ForceField::makeZeoliteForceField(10.0, true, false, true);
  const SimulationBox box(10.0, 10.0, 10.0);
  const std::vector<Atom> atoms = {
      Atom({0.95, 0.5, 0.5}, 0.0, 1.0, 0, 0, 0, 0, true), Atom({0.05, 0.5, 0.5}, 0.0, 1.0, 0, 0, 0, 0, true),
      Atom({0.05, 0.6, 0.5}, 0.0, 1.0, 0, 0, 0, 0, true), Atom({0.05, 0.6, 0.6}, 0.0, 1.0, 0, 0, 0, 0, true)};
  Framework framework(forceField, "periodic-bend-torsion", box, 1, atoms, atoms, int3(1, 1, 1));
  framework.rigid = false;
  framework.intraMolecularPotentials.bends = {BendPotential({0, 1, 2}, BendType::Harmonic, {1500.0, 100.0})};
  framework.intraMolecularPotentials.torsions = {
      TorsionPotential({0, 1, 2, 3}, TorsionType::TraPPE, {0.0, 120.0, -30.0, 200.0})};
  framework.intraMolecularImageShifts.bends = {{{int3(-1, 0, 0), int3(0, 0, 0), int3(0, 0, 0)}}};
  framework.intraMolecularImageShifts.torsions = {{{int3(-1, 0, 0), int3(0, 0, 0), int3(0, 0, 0), int3(0, 0, 0)}}};
  System system(forceField, box, false, 300.0, 0.0, 1.0, framework, {}, {}, {}, 5);
  system.spanOfFrameworkAtoms()[0].position = double3(9.5, 5.0, 5.0);
  system.spanOfFrameworkAtoms()[1].position = double3(0.5, 5.0, 5.0);
  system.spanOfFrameworkAtoms()[2].position = double3(0.5, 6.0, 5.0);
  system.spanOfFrameworkAtoms()[3].position = double3(0.5, 6.0, 6.0);
  system.cellMinimizationType = CellMinimizationType::Regular;
  expectCellDerivativesMatchFiniteDifference(std::move(system), 2.0e-5, 2.0e-3);
}

TEST(minimization_variable_cell, rigid_charged_mixed_blocks_match_finite_difference)
{
  expectCellDerivativesMatchFiniteDifference(makeRigidChargedPair(), 3.0e-4, 5.0e-3);
}

TEST(minimization_variable_cell, polarization_real_space_matches_finite_difference)
{
  expectCellDerivativesMatchFiniteDifference(makePolarizableFlexiblePair(), 3.0e-4, 3.0e-3);
}

TEST(minimization_variable_cell, polarization_reciprocal_framework_matches_finite_difference)
{
  expectCellDerivativesMatchFiniteDifference(makePolarizableRigidFrameworkPair(), 5.0e-4, 5.0e-3);
}

TEST(minimization_variable_cell, polarization_rigid_molecule_real_space_matches_finite_difference)
{
  expectCellDerivativesMatchFiniteDifference(makePolarizableRigidMoleculePair(), 3.0e-4, 5.0e-3);
}

// Isolated polarization stress + strain-strain for a rigid molecule (real-space only): a direct check of the
// analytic center-of-mass-arm derivatives independent of the other interactions' cell-strain Hessians.
TEST(minimization_variable_cell, polarization_rigid_molecule_real_space_strain_matches_finite_difference)
{
  expectPolarizationStrainMatchesFiniteDifference(makePolarizableRigidMoleculePair(), 3.0e-4, 3.0e-3);
}

// Isolated polarization stress + strain-strain for a rigid molecule sampling the fixed-framework reciprocal
// field: validates the rigid field-point body-offset phase terms of the reciprocal strain derivatives. The
// framework-molecule Ewald/vdW rigid-body strain Hessian (a separate, pre-existing path) is not exercised
// here because only the polarization energy is differenced.
TEST(minimization_variable_cell, polarization_rigid_molecule_reciprocal_strain_matches_finite_difference)
{
  expectPolarizationStrainMatchesFiniteDifference(makePolarizableRigidMoleculeInFramework(), 5.0e-4, 5.0e-3);
}

TEST(minimization_variable_cell, fixed_cell_keeps_original_layout_and_energy)
{
  System system = makeLennardJonesPair(CellMinimizationType::Fixed);
  const MinimizationDofLayout oldLayout = buildMinimizationDofLayout(system.moleculeData, system.components);
  const CellMinimizationLayout fixedLayout =
      makeCellMinimizationLayout(system.cellMinimizationType, system.monoclinicAngleType);
  const MinimizationDofLayout newLayout =
      buildMinimizationDofLayout(system.moleculeData, system.components, 0, fixedLayout.size());
  EXPECT_EQ(oldLayout.numDofs(), newLayout.numDofs());
  EXPECT_EQ(newLayout.numberOfCellDofs(), 0u);
  System oldSystem = system;
  System newSystem = system;
  EXPECT_DOUBLE_EQ(evaluate(oldSystem, oldLayout).energy, evaluate(newSystem, newLayout).energy);
}

TEST(minimization_variable_cell, driver_applies_a_clipped_cell_step)
{
  System system = makeLennardJonesPair(CellMinimizationType::Isotropic);
  const double initialVolume = system.simulationBox.volume;
  MinimizationOptions options{};
  options.maximumNumberOfSteps = 2;
  options.maximumStepLength = 0.2;
  options.maximumCellStepLength = 0.02;
  options.convergenceFactor = 0.0;
  options.rmsGradientTolerance = 1.0e-20;
  options.maxGradientTolerance = 1.0e-20;
  Minimization minimization(options, {system}, false);
  minimization.setup();
  EXPECT_THROW(minimization.runPhase(), std::runtime_error);
  EXPECT_GT(minimization.systems[0].simulationBox.volume, 0.0);
  EXPECT_NE(minimization.systems[0].simulationBox.volume, initialVolume);
}
