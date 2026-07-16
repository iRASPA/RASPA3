#include <gtest/gtest.h>

import std;

import double3;
import double3x3;
import simd_quatd;
import atom;
import atom_dynamics;
import molecule;
import forcefield;
import framework;
import int3;
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
// Evaluates the Wolf "reciprocal" energy contribution, which for the finite-cutoff shifted methods is just the
// intra-molecular exclusion / completion energy (there is no Fourier part). The strain-independent self energy
// is deliberately excluded so its large constant offset does not dominate the finite-difference cancellation.
struct WolfExclusionEnergyEvaluator
{
  std::vector<std::complex<double>> eik_x;
  std::vector<std::complex<double>> eik_y;
  std::vector<std::complex<double>> eik_z;
  std::vector<std::complex<double>> eik_xy;
  std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>> fixedFrameworkStoredEik;
  std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>> storedEik;

  double operator()(const System& system, const SimulationBox& box)
  {
    RunningEnergy energy = Interactions::computeEwaldFourierEnergy(
        eik_x, eik_y, eik_z, eik_xy, fixedFrameworkStoredEik, storedEik, system.forceField, box, system.components,
        system.numberOfMoleculesPerComponent, system.spanOfMoleculeAtoms(), 0.0);
    return energy.ewald_fourier + energy.ewald_exclusion;
  }
};

ForceField makeWolfForceField()
{
  ForceField forceField({{"P", false, 15.0, 0.8, 0.0, 8, false}, {"N", false, 14.0, -0.4, 0.0, 8, false}},
                        {{60.0, 3.0}, {40.0, 3.2}}, ForceField::MixingRule::Lorentz_Berthelot, 11.0, 11.0, 11.0, true,
                        false, true);
  forceField.chargeMethod = ForceField::ChargeMethod::Wolf;
  forceField.omitEwaldFourier = true;
  forceField.EwaldAlpha = 0.25;
  return forceField;
}

// Two flexible bent 3-site molecules with zero net charge each (so the exclusion completion is non-trivial).
System makeFlexiblePairSystem(const ForceField& forceField)
{
  ConnectivityTable connectivityTable(3);
  connectivityTable[0, 1] = true;
  connectivityTable[1, 0] = true;
  connectivityTable[1, 2] = true;
  connectivityTable[2, 1] = true;

  Component component = Component(
      forceField, "flexibleIon", 100.0, 1e6, 0.2,
      {Atom({-0.9, 0.4, 0.1}, -0.4, 1.0, 0, 1, 0, false, false), Atom({0.0, 0.0, 0.0}, 0.8, 1.0, 0, 0, 0, false, false),
       Atom({1.0, 0.3, -0.2}, -0.4, 1.0, 0, 1, 0, false, false)},
      connectivityTable, Potentials::IntraMolecularPotentials{}, 5, 21);
  component.rigid = false;

  System system =
      System(forceField, SimulationBox(25.0, 25.0, 25.0), false, 300.0, 1e4, 1.0, {}, {component}, {}, {2}, 5);

  const std::array<double3, 6> positions = {double3(9.1, 10.4, 10.1),  double3(10.0, 10.0, 10.0),
                                            double3(11.0, 10.3, 9.8),  double3(12.6, 11.9, 10.7),
                                            double3(13.5, 11.5, 10.6), double3(14.4, 11.8, 11.3)};
  std::span<Atom> atoms = system.spanOfMoleculeAtoms();
  for (std::size_t i = 0; i < positions.size(); ++i) atoms[i].position = positions[i];
  return system;
}

struct FlexibleDofLabel
{
  std::size_t atom;
  std::size_t axis;
};

std::vector<FlexibleDofLabel> buildFlexibleDofLabels(const MinimizationDofLayout& layout, const System& system)
{
  std::vector<FlexibleDofLabel> labels(layout.numDofs());
  for (std::size_t moleculeIndex = 0; moleculeIndex < system.moleculeData.size(); ++moleculeIndex)
  {
    const Molecule& molecule = system.moleculeData[moleculeIndex];
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
}  // namespace

// The Wolf intra-molecular exclusion / completion Hessian must reproduce finite differences of the exclusion
// energy for flexible multi-site molecules, closing the previous gradient/Hessian inconsistency (the Hessian
// path used to add only the self energy for Wolf). Position-position, position-strain and strain-strain blocks
// are all checked; cross-molecule blocks must vanish (there is no reciprocal coupling for Wolf).
TEST(minimization_hessian_wolf, flexible_charged_pair_matches_finite_difference)
{
  const double delta = 5e-4;
  const double strainStep = 1e-3;
  const double relativeTolerance = 5e-5;

  ForceField forceField = makeWolfForceField();
  System system = makeFlexiblePairSystem(forceField);

  const MinimizationDofLayout layout = buildMinimizationDofLayout(system.moleculeData, system.components);
  ASSERT_EQ(layout.numDofs(), 18u);

  GeneralizedHessian hessian(layout.numDofs(), 1);
  std::vector<AtomDynamics> dynamics(system.spanOfMoleculeAtoms().size());
  const RunningEnergy analytic = Interactions::computeEwaldFourierHessian(
      system.forceField, system.simulationBox, system.framework, system.fixedFrameworkStoredEik,
      system.netChargeFramework, system.moleculeData, system.components, system.spanOfFrameworkAtoms(),
      system.spanOfMoleculeAtoms(), layout, hessian, dynamics);
  EXPECT_EQ(analytic.ewald_fourier, 0.0);
  EXPECT_NE(analytic.ewald_exclusion, 0.0);

  WolfExclusionEnergyEvaluator evaluator{};
  std::span<Atom> atoms = system.spanOfMoleculeAtoms();
  const std::vector<FlexibleDofLabel> labels = buildFlexibleDofLabels(layout, system);
  const double3 boxLengths(25.0, 25.0, 25.0);

  auto energy = [&]() { return evaluator(system, system.simulationBox); };
  auto perturb = [&](const FlexibleDofLabel& label, double amount)
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
    const double value = (&dynamics[labels[dof].atom].gradient.x)[labels[dof].axis];
    const double scale = std::max({1.0, std::abs(numerical), std::abs(value)});
    EXPECT_NEAR(numerical, value, relativeTolerance * scale) << "gradient dof=" << dof;
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
      const double value = hessian(row, col);
      const double scale = std::max({1.0, std::abs(numerical), std::abs(value)});
      EXPECT_NEAR(numerical, value, relativeTolerance * scale) << "row=" << row << " col=" << col;
    }
  }

  // Strain blocks (isotropic exp(epsilon) scaling of positions and cell).
  std::vector<double3> savedPositions(atoms.size());
  auto evaluate = [&](std::optional<FlexibleDofLabel> perturbDof, double positionDelta, double strainExponent) -> double
  {
    for (std::size_t atom = 0; atom < atoms.size(); ++atom) savedPositions[atom] = atoms[atom].position;
    if (perturbDof) (&atoms[perturbDof->atom].position.x)[perturbDof->axis] += positionDelta;
    const double factor = std::exp(strainExponent);
    for (std::size_t atom = 0; atom < atoms.size(); ++atom) atoms[atom].position = atoms[atom].position * factor;
    const SimulationBox strainedBox(factor * boxLengths.x, factor * boxLengths.y, factor * boxLengths.z);
    const double value = evaluator(system, strainedBox);
    for (std::size_t atom = 0; atom < atoms.size(); ++atom) atoms[atom].position = savedPositions[atom];
    return value;
  };

  {
    const double eZero = evaluate(std::nullopt, 0.0, 0.0);
    const double ePlus = evaluate(std::nullopt, 0.0, strainStep);
    const double eMinus = evaluate(std::nullopt, 0.0, -strainStep);
    const double numerical = (ePlus - 2.0 * eZero + eMinus) / (strainStep * strainStep);
    const double value = hessian.strainStrain()[0];
    const double scale = std::max({1.0, std::abs(numerical), std::abs(value)});
    EXPECT_NEAR(numerical, value, 1e-4 * scale) << "strain-strain";
  }

  for (std::size_t dof = 0; dof < labels.size(); ++dof)
  {
    const double ePP = evaluate(labels[dof], delta, strainStep);
    const double ePM = evaluate(labels[dof], delta, -strainStep);
    const double eMP = evaluate(labels[dof], -delta, strainStep);
    const double eMM = evaluate(labels[dof], -delta, -strainStep);
    const double numerical = (ePP - ePM - eMP + eMM) / (4.0 * delta * strainStep);
    const double value = hessian.positionStrain()[dof];
    const double scale = std::max({1.0, std::abs(numerical), std::abs(value)});
    EXPECT_NEAR(numerical, value, 1e-4 * scale) << "position-strain dof=" << dof;
  }
}

// For rigid molecules the intra-molecular exclusion term is invariant under the center-of-mass and
// orientation degrees of freedom, so it must contribute energy but no Hessian curvature (matching the Ewald
// treatment). The gradient on the rigid degrees of freedom must likewise be zero.
TEST(minimization_hessian_wolf, rigid_pair_has_energy_but_no_hessian)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  forceField.chargeMethod = ForceField::ChargeMethod::Wolf;
  forceField.omitEwaldFourier = true;
  forceField.EwaldAlpha = 0.25;

  Component co2 = Component::makeCO2(forceField, 0, true);
  System system = System(forceField, SimulationBox(25.0, 25.0, 25.0), false, 300.0, 1e4, 1.0, {}, {co2}, {}, {2}, 5);

  std::span<Atom> atoms = system.spanOfMoleculeAtoms();
  const std::array<double3, 6> positions = {double3(9.0, 10.0, 8.851),  double3(9.0, 10.0, 10.0),
                                            double3(9.0, 10.0, 11.149), double3(13.0, 11.0, 9.851),
                                            double3(13.0, 11.0, 11.0),  double3(13.0, 11.0, 12.149)};
  for (std::size_t i = 0; i < positions.size(); ++i) atoms[i].position = positions[i];

  const MinimizationDofLayout layout = buildMinimizationDofLayout(system.moleculeData, system.components);
  ASSERT_EQ(layout.numDofs(), 12u);

  GeneralizedHessian hessian(layout.numDofs(), 1);
  std::vector<AtomDynamics> dynamics(system.spanOfMoleculeAtoms().size());
  const RunningEnergy analytic = Interactions::computeEwaldFourierHessian(
      system.forceField, system.simulationBox, system.framework, system.fixedFrameworkStoredEik,
      system.netChargeFramework, system.moleculeData, system.components, system.spanOfFrameworkAtoms(),
      system.spanOfMoleculeAtoms(), layout, hessian, dynamics);

  EXPECT_NE(analytic.ewald_exclusion, 0.0);
  EXPECT_NE(analytic.ewald_self, 0.0);

  for (std::size_t i = 0; i < layout.numDofs(); ++i)
  {
    for (std::size_t j = 0; j < layout.numDofs(); ++j)
    {
      EXPECT_NEAR(hessian(i, j), 0.0, 1e-12) << "row=" << i << " col=" << j;
    }
  }
  for (const AtomDynamics& dyn : dynamics)
  {
    EXPECT_NEAR(dyn.gradient.x, 0.0, 1e-12);
    EXPECT_NEAR(dyn.gradient.y, 0.0, 1e-12);
    EXPECT_NEAR(dyn.gradient.z, 0.0, 1e-12);
  }
}

// A flexible framework's bonded (1-2, 1-3, 1-4) pairs are excluded from the real-space shifted pair sum and
// need the completion q_i q_j (V(r) - 1/r). This validates the framework gradient and the framework
// position-position Hessian block against finite differences of that exclusion energy (mirroring the erf-based
// Ewald flexible-framework treatment, now provided for the Wolf variants as well).
TEST(minimization_hessian_wolf, flexible_framework_excluded_pairs_match_finite_difference)
{
  const double delta = 2e-4;
  const double relativeTolerance = 1e-4;

  ForceField forceField = makeWolfForceField();
  const SimulationBox box(25.0, 25.0, 25.0);

  const std::uint16_t typeP = static_cast<std::uint16_t>(*forceField.findPseudoAtom("P"));
  const std::uint16_t typeN = static_cast<std::uint16_t>(*forceField.findPseudoAtom("N"));

  // A charged three-atom chain (net charge zero) so bonds 0-1, 1-2 and the 0-2 bend are all excluded pairs.
  std::vector<Atom> frameworkAtoms = {Atom(double3(0.30, 0.40, 0.50), 0.8, 1.0, 0, typeP, 0, 0, 1),
                                      Atom(double3(0.34, 0.40, 0.50), -0.4, 1.0, 0, typeN, 0, 0, 1),
                                      Atom(double3(0.39, 0.41, 0.50), -0.4, 1.0, 0, typeN, 0, 0, 1)};
  Framework framework(forceField, "wolf-flexible-framework", box, 1, frameworkAtoms, frameworkAtoms, int3(1, 1, 1));
  framework.rigid = false;
  ConnectivityTable connectivity(3);
  connectivity[0, 1] = true;
  connectivity[1, 0] = true;
  connectivity[1, 2] = true;
  connectivity[2, 1] = true;
  framework.connectivityTable = connectivity;

  Component component(forceField, "ion", 100.0, 1e6, 0.2, {Atom(double3(0.0, 0.0, 0.0), -0.4, 1.0, 0, typeN, 0, 0, 0)},
                     ConnectivityTable(1), Potentials::IntraMolecularPotentials{}, 5, 21);
  component.rigid = false;
  System system(forceField, box, false, 300.0, 1e4, 1.0, {framework}, {component}, {}, {0}, 5);

  const MinimizationDofLayout layout =
      buildMinimizationDofLayout(system.moleculeData, system.components, system.spanOfFrameworkAtoms().size());
  ASSERT_EQ(layout.numDofs(), 9u);

  GeneralizedHessian hessian(layout.numDofs(), 1);
  std::vector<AtomDynamics> moleculeDynamics(system.spanOfMoleculeAtoms().size());
  std::vector<AtomDynamics> frameworkDynamics(system.spanOfFrameworkAtoms().size());
  const RunningEnergy analytic = Interactions::computeEwaldFourierHessian(
      system.forceField, system.simulationBox, system.framework, system.fixedFrameworkStoredEik,
      system.netChargeFramework, system.moleculeData, system.components, system.spanOfFrameworkAtoms(),
      system.spanOfMoleculeAtoms(), layout, hessian, moleculeDynamics, frameworkDynamics);
  EXPECT_EQ(analytic.ewald_fourier, 0.0);
  EXPECT_NE(analytic.ewald_exclusion, 0.0);
  EXPECT_NE(analytic.ewald_self, 0.0);

  // Finite-difference reference: the gradient routine returns the same framework exclusion energy and forces.
  auto exclusionEnergy = [&]() -> double
  {
    std::vector<std::complex<double>> ex, ey, ez, exy;
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>> total, fixed;
    std::vector<AtomDynamics> mScratch(system.spanOfMoleculeAtoms().size());
    std::vector<AtomDynamics> fScratch(system.spanOfFrameworkAtoms().size());
    const RunningEnergy value = Interactions::computeEwaldFourierGradient(
        ex, ey, ez, exy, total, fixed, system.forceField, system.simulationBox, system.components,
        system.numberOfMoleculesPerComponent, system.spanOfMoleculeAtoms(), mScratch, system.netChargeFramework,
        system.framework, system.spanOfFrameworkAtoms(), fScratch);
    return value.ewald_fourier + value.ewald_exclusion;
  };

  std::span<Atom> atoms = system.spanOfFrameworkAtoms();
  auto perturb = [&](std::size_t atom, std::size_t axis, double amount) { (&atoms[atom].position.x)[axis] += amount; };
  auto dofOf = [&](std::size_t atom, std::size_t axis)
  { return *layout.frameworkAtomDof(atom, static_cast<MinimizationDofAxis>(axis)); };

  // Gradient check.
  for (std::size_t atom = 0; atom < atoms.size(); ++atom)
  {
    for (std::size_t axis = 0; axis < 3; ++axis)
    {
      perturb(atom, axis, delta);
      const double ePlus = exclusionEnergy();
      perturb(atom, axis, -2.0 * delta);
      const double eMinus = exclusionEnergy();
      perturb(atom, axis, delta);
      const double numerical = (ePlus - eMinus) / (2.0 * delta);
      const double value = (&frameworkDynamics[atom].gradient.x)[axis];
      const double scale = std::max({1.0, std::abs(numerical), std::abs(value)});
      EXPECT_NEAR(numerical, value, relativeTolerance * scale) << "gradient atom=" << atom << " axis=" << axis;
    }
  }

  // Position-position block.
  for (std::size_t rowAtom = 0; rowAtom < atoms.size(); ++rowAtom)
  {
    for (std::size_t rowAxis = 0; rowAxis < 3; ++rowAxis)
    {
      for (std::size_t colAtom = 0; colAtom < atoms.size(); ++colAtom)
      {
        for (std::size_t colAxis = 0; colAxis < 3; ++colAxis)
        {
          perturb(rowAtom, rowAxis, delta);
          perturb(colAtom, colAxis, delta);
          const double ePP = exclusionEnergy();
          perturb(colAtom, colAxis, -2.0 * delta);
          const double ePM = exclusionEnergy();
          perturb(rowAtom, rowAxis, -2.0 * delta);
          const double eMM = exclusionEnergy();
          perturb(colAtom, colAxis, 2.0 * delta);
          const double eMP = exclusionEnergy();
          perturb(rowAtom, rowAxis, delta);
          perturb(colAtom, colAxis, -delta);

          const double numerical = (ePP - ePM - eMP + eMM) / (4.0 * delta * delta);
          const double value = hessian(dofOf(rowAtom, rowAxis), dofOf(colAtom, colAxis));
          const double scale = std::max({1.0, std::abs(numerical), std::abs(value)});
          EXPECT_NEAR(numerical, value, relativeTolerance * scale)
              << "row=(" << rowAtom << "," << rowAxis << ") col=(" << colAtom << "," << colAxis << ")";
        }
      }
    }
  }
}
