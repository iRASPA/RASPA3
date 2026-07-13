#include <gtest/gtest.h>

import std;

import double3;
import atom;
import atom_dynamics;
import molecule;
import forcefield;
import component;
import system;
import simulationbox;
import connectivity_table;
import intra_molecular_potentials;
import bond_bond_potential;
import bond_bend_potential;
import bend_bend_potential;
import bond_torsion_potential;
import bend_torsion_potential;
import inversion_bend_potential;
import generalized_hessian;
import minimization_dof_layout;
import interactions_hessian_intramolecular;

namespace
{
ForceField makeTestForceField()
{
  return ForceField({{"CH4", false, 16.04246, 0.0, 0.0, 8, false},
                     {"CH3", false, 15.03452, 0.0, 0.0, 8, false},
                     {"CH2", false, 14.02658, 0.0, 0.0, 8, false},
                     {"CH", false, 13.01864, 0.0, 0.0, 8, false},
                     {"C", false, 12.0, 0.0, 0.0, 8, false}},
                    {{158.5, 3.72}, {108.0, 3.76}, {56.0, 3.96}, {17.0, 4.67}, {0.8, 6.38}},
                    ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true, false, false);
}

void perturbAtom(std::span<Atom> atoms, std::size_t atom, std::size_t axis, double delta)
{
  (&atoms[atom].position.x)[axis] += delta;
}

struct DofLabel
{
  std::size_t atom;
  std::size_t axis;
};

std::vector<DofLabel> buildDofLabels(const MinimizationDofLayout &layout, std::size_t numAtoms)
{
  std::vector<DofLabel> labels(layout.numDofs());
  for (std::size_t dof = 0; dof < layout.numDofs(); ++dof)
  {
    bool found = false;
    for (std::size_t atom = 0; atom < numAtoms && !found; ++atom)
    {
      for (std::size_t axis = 0; axis < 3; ++axis)
      {
        if (layout.flexibleAtomDof(0, atom, static_cast<MinimizationDofAxis>(axis)) == dof)
        {
          labels[dof] = {atom, axis};
          found = true;
          break;
        }
      }
    }
    EXPECT_TRUE(found);
  }
  return labels;
}

void compareHessianToFiniteDifference(std::span<Atom> atoms, const std::vector<DofLabel> &labels,
                                      const GeneralizedHessian &hessian,
                                      const std::function<double()> &energyFunction, double delta,
                                      double relativeTolerance)
{
  const std::size_t numDofs = labels.size();
  for (std::size_t row = 0; row < numDofs; ++row)
  {
    for (std::size_t col = 0; col < numDofs; ++col)
    {
      const DofLabel rowLabel = labels[row];
      const DofLabel colLabel = labels[col];
      double numerical{};

      if (row == col)
      {
        const double eZero = energyFunction();
        perturbAtom(atoms, rowLabel.atom, rowLabel.axis, delta);
        const double ePlus = energyFunction();
        perturbAtom(atoms, rowLabel.atom, rowLabel.axis, -2.0 * delta);
        const double eMinus = energyFunction();
        perturbAtom(atoms, rowLabel.atom, rowLabel.axis, delta);
        numerical = (ePlus - 2.0 * eZero + eMinus) / (delta * delta);
      }
      else
      {
        perturbAtom(atoms, rowLabel.atom, rowLabel.axis, delta);
        perturbAtom(atoms, colLabel.atom, colLabel.axis, delta);
        const double ePP = energyFunction();
        perturbAtom(atoms, colLabel.atom, colLabel.axis, -2.0 * delta);
        const double ePM = energyFunction();
        perturbAtom(atoms, rowLabel.atom, rowLabel.axis, -2.0 * delta);
        const double eMM = energyFunction();
        perturbAtom(atoms, colLabel.atom, colLabel.axis, 2.0 * delta);
        const double eMP = energyFunction();
        perturbAtom(atoms, rowLabel.atom, rowLabel.axis, delta);
        perturbAtom(atoms, colLabel.atom, colLabel.axis, -delta);
        numerical = (ePP - ePM - eMP + eMM) / (4.0 * delta * delta);
      }

      const double analytic = hessian(row, col);
      const double scale = std::max({1.0, std::abs(numerical), std::abs(analytic)});
      EXPECT_NEAR(numerical, analytic, relativeTolerance * scale) << "row=" << row << " col=" << col;
    }
  }
}

// Finite-difference check of the isotropic strain blocks. The strain scalar epsilon scales all
// atom positions by exp(epsilon) about the origin, matching the convention of the pair terms.
void compareStrainToFiniteDifference(std::span<Atom> atoms, std::size_t numAtoms, const std::vector<DofLabel> &labels,
                                     const GeneralizedHessian &hessian,
                                     const std::function<double()> &energyFunction, double delta, double strainStep,
                                     double relativeTolerance)
{
  std::vector<double3> savedPositions(numAtoms);

  auto evaluate = [&](std::optional<DofLabel> perturb, double positionDelta, double strainExponent) -> double
  {
    for (std::size_t atom = 0; atom < numAtoms; ++atom) savedPositions[atom] = atoms[atom].position;
    if (perturb) (&atoms[perturb->atom].position.x)[perturb->axis] += positionDelta;
    const double factor = std::exp(strainExponent);
    for (std::size_t atom = 0; atom < numAtoms; ++atom) atoms[atom].position = atoms[atom].position * factor;
    const double energy = energyFunction();
    for (std::size_t atom = 0; atom < numAtoms; ++atom) atoms[atom].position = savedPositions[atom];
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
    const DofLabel label = labels[dof];
    const double ePP = evaluate(label, delta, strainStep);
    const double ePM = evaluate(label, delta, -strainStep);
    const double eMP = evaluate(label, -delta, strainStep);
    const double eMM = evaluate(label, -delta, -strainStep);
    const double numerical = (ePP - ePM - eMP + eMM) / (4.0 * delta * strainStep);
    const double analytic = hessian.positionStrain()[dof];
    const double scale = std::max({1.0, std::abs(numerical), std::abs(analytic)});
    EXPECT_NEAR(numerical, analytic, relativeTolerance * scale) << "position-strain dof=" << dof;
  }
}

// A bent 4-atom chain with non-degenerate bond, bend, and dihedral geometry.
System makeChainSystem(const ForceField &forceField, const Potentials::IntraMolecularPotentials &intraMolecularPotentials)
{
  ConnectivityTable connectivityTable(4);
  connectivityTable[0, 1] = true;
  connectivityTable[1, 0] = true;
  connectivityTable[1, 2] = true;
  connectivityTable[2, 1] = true;
  connectivityTable[2, 3] = true;
  connectivityTable[3, 2] = true;

  Component component = Component(forceField, "chain", 425.125, 3796000.0, 0.201,
                                  {Atom({-0.5, 1.4, 0.3}, 0.0, 1.0, 0, 1, 0, false, false),
                                   Atom({0.0, 0.0, 0.0}, 0.0, 1.0, 0, 2, 0, false, false),
                                   Atom({1.54, 0.0, 0.0}, 0.0, 1.0, 0, 2, 0, false, false),
                                   Atom({2.1, 1.2, 1.0}, 0.0, 1.0, 0, 1, 0, false, false)},
                                  connectivityTable, intraMolecularPotentials, 5, 21);
  component.rigid = false;

  return System(forceField, SimulationBox(25.0, 25.0, 25.0), false, 300.0, 1e4, 1.0, {}, {component}, {}, {1}, 5);
}

void runCrossTest(const Potentials::IntraMolecularPotentials &intraMolecularPotentials,
                  const std::function<void(std::span<const Molecule>, std::span<Atom>, std::span<const Component>,
                                           const MinimizationDofLayout &, GeneralizedHessian &,
                                           std::span<AtomDynamics>)> &computeFunction,
                  const std::function<double(std::span<const Atom>, const Potentials::IntraMolecularPotentials &)> &
                      energyEvaluator,
                  double delta = 1e-5, double relativeTolerance = 5e-3)
{
  ForceField forceField = makeTestForceField();
  System system = makeChainSystem(forceField, intraMolecularPotentials);

  const MinimizationDofLayout layout = buildMinimizationDofLayout(system.moleculeData, system.components);
  ASSERT_EQ(layout.numDofs(), 12u);

  // One isotropic strain degree of freedom so the strain blocks are populated as well.
  GeneralizedHessian hessian(layout.numDofs(), 1);
  std::span<Atom> atoms = system.spanOfMoleculeAtoms();
  std::span<AtomDynamics> dynamics = system.spanOfMoleculeDynamics();

  computeFunction(system.moleculeData, atoms, system.components, layout, hessian, dynamics);

  const Potentials::IntraMolecularPotentials &potentials = system.components[0].intraMolecularPotentials;
  auto energyFunction = [&]() -> double { return energyEvaluator(atoms, potentials); };

  const std::vector<DofLabel> labels = buildDofLabels(layout, 4);
  compareHessianToFiniteDifference(atoms, labels, hessian, energyFunction, delta, relativeTolerance);
  // The strain finite differences use a moderate step: scale-invariant angle/dihedral coordinates
  // make the true strain coupling exactly zero, so an overly small step would let floating-point
  // cancellation dominate, while the exponential bond-bend terms need a step small enough to keep
  // truncation error in check.
  compareStrainToFiniteDifference(atoms, 4, labels, hessian, energyFunction, 2e-3, 2e-3, relativeTolerance);
}
}  // namespace

TEST(minimization_hessian_cross, bond_bond_cvff)
{
  Potentials::IntraMolecularPotentials potentials{};
  potentials.bondBonds = {BondBondPotential({0, 1, 2}, BondBondType::CVFF, {50000.0, 1.5, 1.5})};

  runCrossTest(
      potentials, &Interactions::computeIntraMolecularBondBondHessian,
      [](std::span<const Atom> atoms, const Potentials::IntraMolecularPotentials &p)
      {
        double energy = 0.0;
        for (const BondBondPotential &t : p.bondBonds)
          energy += t.calculateEnergy(atoms[t.identifiers[0]].position, atoms[t.identifiers[1]].position,
                                      atoms[t.identifiers[2]].position);
        return energy;
      });
}

TEST(minimization_hessian_cross, bond_bend_cvff)
{
  Potentials::IntraMolecularPotentials potentials{};
  potentials.bondBends = {BondBendPotential({0, 1, 2, 3}, BondBendType::CVFF, {109.5, 40000.0, 1.5, 40000.0, 1.5})};

  runCrossTest(
      potentials, &Interactions::computeIntraMolecularBondBendHessian,
      [](std::span<const Atom> atoms, const Potentials::IntraMolecularPotentials &p)
      {
        double energy = 0.0;
        for (const BondBendPotential &t : p.bondBends)
          energy += t.calculateEnergy(atoms[t.identifiers[0]].position, atoms[t.identifiers[1]].position,
                                      atoms[t.identifiers[2]].position, atoms[t.identifiers[3]].position);
        return energy;
      });
}

TEST(minimization_hessian_cross, bond_bend_mm3)
{
  Potentials::IntraMolecularPotentials potentials{};
  potentials.bondBends = {BondBendPotential({0, 1, 2, 3}, BondBendType::MM3, {1.2, 1.5, 1.5, 109.5})};

  runCrossTest(
      potentials, &Interactions::computeIntraMolecularBondBendHessian,
      [](std::span<const Atom> atoms, const Potentials::IntraMolecularPotentials &p)
      {
        double energy = 0.0;
        for (const BondBendPotential &t : p.bondBends)
          energy += t.calculateEnergy(atoms[t.identifiers[0]].position, atoms[t.identifiers[1]].position,
                                      atoms[t.identifiers[2]].position, atoms[t.identifiers[3]].position);
        return energy;
      });
}

TEST(minimization_hessian_cross, bond_bend_truncated_harmonic)
{
  Potentials::IntraMolecularPotentials potentials{};
  potentials.bondBends = {BondBendPotential({0, 1, 2, 3}, BondBendType::TruncatedHarmonic, {60000.0, 109.5, 2.6})};

  runCrossTest(
      potentials, &Interactions::computeIntraMolecularBondBendHessian,
      [](std::span<const Atom> atoms, const Potentials::IntraMolecularPotentials &p)
      {
        double energy = 0.0;
        for (const BondBendPotential &t : p.bondBends)
          energy += t.calculateEnergy(atoms[t.identifiers[0]].position, atoms[t.identifiers[1]].position,
                                      atoms[t.identifiers[2]].position, atoms[t.identifiers[3]].position);
        return energy;
      });
}

TEST(minimization_hessian_cross, bond_bend_screened_harmonic)
{
  Potentials::IntraMolecularPotentials potentials{};
  potentials.bondBends = {BondBendPotential({0, 1, 2, 3}, BondBendType::ScreenedHarmonic, {60000.0, 109.5, 1.2, 1.2})};

  runCrossTest(
      potentials, &Interactions::computeIntraMolecularBondBendHessian,
      [](std::span<const Atom> atoms, const Potentials::IntraMolecularPotentials &p)
      {
        double energy = 0.0;
        for (const BondBendPotential &t : p.bondBends)
          energy += t.calculateEnergy(atoms[t.identifiers[0]].position, atoms[t.identifiers[1]].position,
                                      atoms[t.identifiers[2]].position, atoms[t.identifiers[3]].position);
        return energy;
      });
}

TEST(minimization_hessian_cross, bond_bend_screened_vessal)
{
  Potentials::IntraMolecularPotentials potentials{};
  potentials.bondBends = {BondBendPotential({0, 1, 2, 3}, BondBendType::ScreenedVessal, {60000.0, 109.5, 1.2, 1.2})};

  runCrossTest(
      potentials, &Interactions::computeIntraMolecularBondBendHessian,
      [](std::span<const Atom> atoms, const Potentials::IntraMolecularPotentials &p)
      {
        double energy = 0.0;
        for (const BondBendPotential &t : p.bondBends)
          energy += t.calculateEnergy(atoms[t.identifiers[0]].position, atoms[t.identifiers[1]].position,
                                      atoms[t.identifiers[2]].position, atoms[t.identifiers[3]].position);
        return energy;
      });
}

TEST(minimization_hessian_cross, bond_bend_truncated_vessal)
{
  Potentials::IntraMolecularPotentials potentials{};
  potentials.bondBends = {BondBendPotential({0, 1, 2, 3}, BondBendType::TruncatedVessal, {60000.0, 109.5, 4.0, 2.6})};

  runCrossTest(
      potentials, &Interactions::computeIntraMolecularBondBendHessian,
      [](std::span<const Atom> atoms, const Potentials::IntraMolecularPotentials &p)
      {
        double energy = 0.0;
        for (const BondBendPotential &t : p.bondBends)
          energy += t.calculateEnergy(atoms[t.identifiers[0]].position, atoms[t.identifiers[1]].position,
                                      atoms[t.identifiers[2]].position, atoms[t.identifiers[3]].position);
        return energy;
      });
}

TEST(minimization_hessian_cross, bend_bend_cvff)
{
  Potentials::IntraMolecularPotentials potentials{};
  potentials.bendBends = {BendBendPotential({0, 1, 2, 3}, BendBendType::CVFF, {50000.0, 109.5, 109.5})};

  runCrossTest(
      potentials, &Interactions::computeIntraMolecularBendBendHessian,
      [](std::span<const Atom> atoms, const Potentials::IntraMolecularPotentials &p)
      {
        double energy = 0.0;
        for (const BendBendPotential &t : p.bendBends)
          energy += t.calculateEnergy(atoms[t.identifiers[0]].position, atoms[t.identifiers[1]].position,
                                      atoms[t.identifiers[2]].position, atoms[t.identifiers[3]].position);
        return energy;
      });
}

TEST(minimization_hessian_cross, bend_bend_mm3)
{
  Potentials::IntraMolecularPotentials potentials{};
  potentials.bendBends = {BendBendPotential({0, 1, 2, 3}, BendBendType::MM3, {1.5, 109.5, 109.5})};

  runCrossTest(
      potentials, &Interactions::computeIntraMolecularBendBendHessian,
      [](std::span<const Atom> atoms, const Potentials::IntraMolecularPotentials &p)
      {
        double energy = 0.0;
        for (const BendBendPotential &t : p.bendBends)
          energy += t.calculateEnergy(atoms[t.identifiers[0]].position, atoms[t.identifiers[1]].position,
                                      atoms[t.identifiers[2]].position, atoms[t.identifiers[3]].position);
        return energy;
      });
}

TEST(minimization_hessian_cross, bond_torsion_mm3)
{
  Potentials::IntraMolecularPotentials potentials{};
  potentials.bondTorsions = {BondTorsionPotential({0, 1, 2, 3}, BondTorsionType::MM3, {2.0, 1.0, 0.5, 1.54})};

  runCrossTest(
      potentials, &Interactions::computeIntraMolecularBondTorsionHessian,
      [](std::span<const Atom> atoms, const Potentials::IntraMolecularPotentials &p)
      {
        double energy = 0.0;
        for (const BondTorsionPotential &t : p.bondTorsions)
          energy += t.calculateEnergy(atoms[t.identifiers[0]].position, atoms[t.identifiers[1]].position,
                                      atoms[t.identifiers[2]].position, atoms[t.identifiers[3]].position);
        return energy;
      });
}

TEST(minimization_hessian_cross, bend_torsion_cvff)
{
  Potentials::IntraMolecularPotentials potentials{};
  potentials.bendTorsions = {BendTorsionPotential({0, 1, 2, 3}, BendTorsionType::CVFF, {5000.0, 109.5, 109.5})};

  runCrossTest(
      potentials, &Interactions::computeIntraMolecularBendTorsionHessian,
      [](std::span<const Atom> atoms, const Potentials::IntraMolecularPotentials &p)
      {
        double energy = 0.0;
        for (const BendTorsionPotential &t : p.bendTorsions)
          energy += t.calculateEnergy(atoms[t.identifiers[0]].position, atoms[t.identifiers[1]].position,
                                      atoms[t.identifiers[2]].position, atoms[t.identifiers[3]].position);
        return energy;
      });
}

TEST(minimization_hessian_cross, bend_torsion_smoothed_three_cosine)
{
  Potentials::IntraMolecularPotentials potentials{};
  potentials.bendTorsions = {
      BendTorsionPotential({0, 1, 2, 3}, BendTorsionType::SmoothedThreeCosine, {1000.0, 500.0, 800.0})};

  runCrossTest(
      potentials, &Interactions::computeIntraMolecularBendTorsionHessian,
      [](std::span<const Atom> atoms, const Potentials::IntraMolecularPotentials &p)
      {
        double energy = 0.0;
        for (const BendTorsionPotential &t : p.bendTorsions)
          energy += t.calculateEnergy(atoms[t.identifiers[0]].position, atoms[t.identifiers[1]].position,
                                      atoms[t.identifiers[2]].position, atoms[t.identifiers[3]].position);
        return energy;
      });
}

TEST(minimization_hessian_cross, bend_torsion_nicholas)
{
  Potentials::IntraMolecularPotentials potentials{};
  potentials.bendTorsions = {BendTorsionPotential({0, 1, 2, 3}, BendTorsionType::Nicholas, {1000.0, 500.0, 800.0})};

  runCrossTest(
      potentials, &Interactions::computeIntraMolecularBendTorsionHessian,
      [](std::span<const Atom> atoms, const Potentials::IntraMolecularPotentials &p)
      {
        double energy = 0.0;
        for (const BendTorsionPotential &t : p.bendTorsions)
          energy += t.calculateEnergy(atoms[t.identifiers[0]].position, atoms[t.identifiers[1]].position,
                                      atoms[t.identifiers[2]].position, atoms[t.identifiers[3]].position);
        return energy;
      });
}

TEST(minimization_hessian_cross, bend_torsion_smoothed_cff)
{
  Potentials::IntraMolecularPotentials potentials{};
  potentials.bendTorsions = {BendTorsionPotential({0, 1, 2, 3}, BendTorsionType::SmoothedCFF, {1000.0, 500.0, 800.0})};

  runCrossTest(
      potentials, &Interactions::computeIntraMolecularBendTorsionHessian,
      [](std::span<const Atom> atoms, const Potentials::IntraMolecularPotentials &p)
      {
        double energy = 0.0;
        for (const BendTorsionPotential &t : p.bendTorsions)
          energy += t.calculateEnergy(atoms[t.identifiers[0]].position, atoms[t.identifiers[1]].position,
                                      atoms[t.identifiers[2]].position, atoms[t.identifiers[3]].position);
        return energy;
      });
}

TEST(minimization_hessian_cross, bend_torsion_smoothed_cff3)
{
  Potentials::IntraMolecularPotentials potentials{};
  potentials.bendTorsions = {
      BendTorsionPotential({0, 1, 2, 3}, BendTorsionType::SmoothedCFF3, {5000.0, 109.5, 109.5})};

  runCrossTest(
      potentials, &Interactions::computeIntraMolecularBendTorsionHessian,
      [](std::span<const Atom> atoms, const Potentials::IntraMolecularPotentials &p)
      {
        double energy = 0.0;
        for (const BendTorsionPotential &t : p.bendTorsions)
          energy += t.calculateEnergy(atoms[t.identifiers[0]].position, atoms[t.identifiers[1]].position,
                                      atoms[t.identifiers[2]].position, atoms[t.identifiers[3]].position);
        return energy;
      });
}

TEST(minimization_hessian_cross, bend_torsion_smoothed_cff2)
{
  Potentials::IntraMolecularPotentials potentials{};
  potentials.bendTorsions = {BendTorsionPotential({0, 1, 2, 3}, BendTorsionType::SmoothedCFF2, {1000.0, 500.0, 800.0})};

  runCrossTest(
      potentials, &Interactions::computeIntraMolecularBendTorsionHessian,
      [](std::span<const Atom> atoms, const Potentials::IntraMolecularPotentials &p)
      {
        double energy = 0.0;
        for (const BendTorsionPotential &t : p.bendTorsions)
          energy += t.calculateEnergy(atoms[t.identifiers[0]].position, atoms[t.identifiers[1]].position,
                                      atoms[t.identifiers[2]].position, atoms[t.identifiers[3]].position);
        return energy;
      });
}

TEST(minimization_hessian_cross, bend_torsion_smoothed_phi)
{
  Potentials::IntraMolecularPotentials potentials{};
  potentials.bendTorsions = {BendTorsionPotential({0, 1, 2, 3}, BendTorsionType::Smoothed, {1500.0, 3.0, 0.0})};

  runCrossTest(
      potentials, &Interactions::computeIntraMolecularBendTorsionHessian,
      [](std::span<const Atom> atoms, const Potentials::IntraMolecularPotentials &p)
      {
        double energy = 0.0;
        for (const BendTorsionPotential &t : p.bendTorsions)
          energy += t.calculateEnergy(atoms[t.identifiers[0]].position, atoms[t.identifiers[1]].position,
                                      atoms[t.identifiers[2]].position, atoms[t.identifiers[3]].position);
        return energy;
      });
}

namespace
{
void runInversionBendTest(InversionBendType type, std::vector<double> parameters)
{
  Potentials::IntraMolecularPotentials potentials{};
  potentials.inversionBends = {InversionBendPotential({0, 1, 2, 3}, type, parameters)};

  runCrossTest(
      potentials, &Interactions::computeIntraMolecularInversionBendHessian,
      [](std::span<const Atom> atoms, const Potentials::IntraMolecularPotentials &p)
      {
        double energy = 0.0;
        for (const InversionBendPotential &t : p.inversionBends)
          energy += t.calculateEnergy(atoms[t.identifiers[0]].position, atoms[t.identifiers[1]].position,
                                      atoms[t.identifiers[2]].position, atoms[t.identifiers[3]].position);
        return energy;
      });
}
}  // namespace

TEST(minimization_hessian_cross, inversion_bend_harmonic)
{
  runInversionBendTest(InversionBendType::Harmonic, {60000.0, 25.0});
}

TEST(minimization_hessian_cross, inversion_bend_harmonic_cosine)
{
  runInversionBendTest(InversionBendType::HarmonicCosine, {60000.0, 25.0});
}

TEST(minimization_hessian_cross, inversion_bend_planar)
{
  runInversionBendTest(InversionBendType::Planar, {60000.0});
}

TEST(minimization_hessian_cross, inversion_bend_harmonic2)
{
  runInversionBendTest(InversionBendType::Harmonic2, {60000.0, 25.0});
}

TEST(minimization_hessian_cross, inversion_bend_harmonic_cosine2)
{
  runInversionBendTest(InversionBendType::HarmonicCosine2, {60000.0, 25.0});
}

TEST(minimization_hessian_cross, inversion_bend_planar2)
{
  runInversionBendTest(InversionBendType::Planar2, {60000.0});
}

TEST(minimization_hessian_cross, inversion_bend_mm3)
{
  runInversionBendTest(InversionBendType::MM3, {0.5, 25.0});
}
