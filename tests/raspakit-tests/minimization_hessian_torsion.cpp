#include <gtest/gtest.h>

import std;

import double3;
import atom;
import atom_dynamics;
import forcefield;
import component;
import system;
import simulationbox;
import connectivity_table;
import intra_molecular_potentials;
import bond_potential;
import torsion_potential;
import dihedral_derivatives;
import urey_bradley_potential;
import van_der_waals_potential;
import vdwparameters;
import coulomb_potential;
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

const std::array<double3, 4> gaucheButane{double3(-0.5, 1.4, 0.3), double3(0.0, 0.0, 0.0), double3(1.54, 0.0, 0.0),
                                          double3(2.1, 1.2, 1.0)};

// all four atoms in the z = 0 plane, phi = 0
const std::array<double3, 4> cisButane{double3(-0.5, 1.4, 0.0), double3(0.0, 0.0, 0.0), double3(1.54, 0.0, 0.0),
                                       double3(2.04, 1.4, 0.0)};

System makeButaneSystem(const ForceField& forceField, TorsionType torsionType, std::vector<double> parameters,
                        const std::array<double3, 4>& positions)
{
  ConnectivityTable connectivityTable(4);
  connectivityTable[0, 1] = true;
  connectivityTable[1, 0] = true;
  connectivityTable[1, 2] = true;
  connectivityTable[2, 1] = true;
  connectivityTable[2, 3] = true;
  connectivityTable[3, 2] = true;

  Potentials::IntraMolecularPotentials intraMolecularPotentials{};
  intraMolecularPotentials.torsions = {TorsionPotential({0, 1, 2, 3}, torsionType, std::move(parameters))};

  Component component = Component(
      forceField, "butane", 425.125, 3796000.0, 0.201,
      {Atom(positions[0], 0.0, 1.0, 0, 1, 0, false, false), Atom(positions[1], 0.0, 1.0, 0, 2, 0, false, false),
       Atom(positions[2], 0.0, 1.0, 0, 2, 0, false, false), Atom(positions[3], 0.0, 1.0, 0, 1, 0, false, false)},
      connectivityTable, intraMolecularPotentials, 5, 21);
  component.rigid = false;

  return System(forceField, SimulationBox(25.0, 25.0, 25.0), false, 300.0, 1e4, 1.0, {}, {component}, {}, {1}, 5);
}

void runTorsionFiniteDifferenceTest(TorsionType torsionType, std::vector<double> parameters,
                                    const std::array<double3, 4>& positions = gaucheButane)
{
  const double delta = 1e-5;
  const double relativeTolerance = 5e-3;

  ForceField forceField = makeTestForceField();
  System system = makeButaneSystem(forceField, torsionType, std::move(parameters), positions);

  const MinimizationDofLayout layout = buildMinimizationDofLayout(system.moleculeData, system.components);
  ASSERT_EQ(layout.numDofs(), 12u);

  GeneralizedHessian hessian(layout.numDofs(), 0);
  std::span<Atom> atoms = system.spanOfMoleculeAtoms();
  std::span<AtomDynamics> dynamics = system.spanOfMoleculeDynamics();

  Interactions::computeIntraMolecularTorsionHessian(system.moleculeData, atoms, system.components, layout, hessian,
                                                    dynamics);

  const Potentials::IntraMolecularPotentials &potentials =
      system.components[0].intraMolecularPotentials;
  auto energyFunction = [&]() -> double
  {
    double energy = 0.0;
    for (const TorsionPotential &torsion : potentials.torsions)
    {
      energy += torsion.calculateEnergy(atoms[torsion.identifiers[0]].position, atoms[torsion.identifiers[1]].position,
                                        atoms[torsion.identifiers[2]].position, atoms[torsion.identifiers[3]].position);
    }
    return energy;
  };

  const std::vector<DofLabel> labels = buildDofLabels(layout, 4);
  compareHessianToFiniteDifference(atoms, labels, hessian, energyFunction, delta, relativeTolerance);
}
}  // namespace

TEST(minimization_hessian_torsion, trappe_butane_flexible_analytic_matches_finite_difference)
{
  runTorsionFiniteDifferenceTest(TorsionType::TraPPE, {0.0, 355.03, -68.19, 791.32});
}

TEST(minimization_hessian_torsion, harmonic_dihedral_flexible_analytic_matches_finite_difference)
{
  runTorsionFiniteDifferenceTest(TorsionType::Harmonic, {5000.0, 65.0});
}

TEST(minimization_hessian_torsion, harmonic_cosine_flexible_analytic_matches_finite_difference)
{
  runTorsionFiniteDifferenceTest(TorsionType::HarmonicCosine, {800.0, 60.0});
}

// A potential in the signed angle phi at a planar (cis) dihedral: the chain rule through cos(phi) divided by
// sin(phi) = 0 there, so the out-of-plane force U'(phi) grad(phi) and the curvature U'' grad(phi) grad(phi)^T were lost
// (zeroed or NaN). The exact phi derivatives keep both.
TEST(minimization_hessian_torsion, phi_potentials_at_a_cis_dihedral_match_finite_difference)
{
  runTorsionFiniteDifferenceTest(TorsionType::Harmonic, {5000.0, 65.0}, cisButane);
  runTorsionFiniteDifferenceTest(TorsionType::ModifiedTraPPE, {0.0, 355.03, -68.19, 791.32, 40.0}, cisButane);
  runTorsionFiniteDifferenceTest(TorsionType::CVFF, {1000.0, 2.0, 50.0}, cisButane);
}

TEST(minimization_hessian_torsion, phi_potential_gradient_at_a_cis_dihedral_matches_finite_difference)
{
  for (const auto& [type, parameters] : std::vector<std::pair<TorsionType, std::vector<double>>>{
           {TorsionType::Harmonic, {5000.0, 65.0}},
           {TorsionType::ModifiedTraPPE, {0.0, 355.03, -68.19, 791.32, 40.0}},
           {TorsionType::CVFF, {1000.0, 2.0, 50.0}}})
  {
    const TorsionPotential torsion({0, 1, 2, 3}, type, parameters);
    std::array<double3, 4> positions = cisButane;
    const auto [energy, gradient, strain] =
        torsion.potentialEnergyGradientStrain(positions[0], positions[1], positions[2], positions[3]);
    EXPECT_NEAR(energy, torsion.calculateEnergy(positions[0], positions[1], positions[2], positions[3]), 1e-10);

    const double delta = 1e-6;
    double outOfPlane = 0.0;
    for (std::size_t atom = 0; atom < 4; ++atom)
    {
      for (std::size_t axis = 0; axis < 3; ++axis)
      {
        (&positions[atom].x)[axis] += delta;
        const double ePlus = torsion.calculateEnergy(positions[0], positions[1], positions[2], positions[3]);
        (&positions[atom].x)[axis] -= 2.0 * delta;
        const double eMinus = torsion.calculateEnergy(positions[0], positions[1], positions[2], positions[3]);
        (&positions[atom].x)[axis] += delta;
        const double numerical = (ePlus - eMinus) / (2.0 * delta);
        const double analytic = (&gradient[atom].x)[axis];
        EXPECT_NEAR(numerical, analytic, 1e-5 * std::max(1.0, std::abs(numerical)))
            << static_cast<int>(type) << " atom=" << atom << " axis=" << axis;
        if (axis == 2) outOfPlane += std::abs(analytic);
      }
    }
    EXPECT_GT(outOfPlane, 1.0) << static_cast<int>(type);  // the force that rotates the dihedral out of the plane
  }
}

// The exact phi derivatives against finite differences at a gauche, a cis and a near-trans dihedral, and the angle
// against the IUPAC sign convention of the energy routines.
TEST(dihedral_derivatives, match_finite_differences_and_the_iupac_sign)
{
  const std::array<double3, 4> nearTrans{double3(-0.5, 1.4, 0.0), double3(0.0, 0.0, 0.0), double3(1.54, 0.0, 0.0),
                                         double3(2.04, -1.4, 0.02)};
  for (const std::array<double3, 4>& geometry : {gaucheButane, cisButane, nearTrans})
  {
    std::array<double3, 4> p = geometry;
    const Potentials::DihedralDerivatives exact = Potentials::dihedralAngleDerivatives(p[0], p[1], p[2], p[3], true);
    EXPECT_NEAR(exact.phi, Potentials::dihedralAngle(p[0], p[1], p[2], p[3]), 1e-14);

    // reference IUPAC angle: copysign(acos(cos phi), -Rcb . ((Rab x Rcb) x (Rcb x Rdc)))
    const double3 Dab = p[0] - p[1];
    const double3 Dcb = (p[2] - p[1]).normalized();
    const double3 Ddc = p[3] - p[2];
    const double3 dr = (Dab - double3::dot(Dab, Dcb) * Dcb).normalized();
    const double3 ds = (Ddc - double3::dot(Ddc, Dcb) * Dcb).normalized();
    const double reference =
        std::copysign(std::acos(std::clamp(double3::dot(dr, ds), -1.0, 1.0)),
                      -double3::dot(Dcb, double3::cross(double3::cross(Dab, Dcb), double3::cross(Dcb, Ddc))));
    EXPECT_NEAR(exact.phi, reference, 1e-7);

    const double delta = 1e-5;
    auto phiAt = [&](std::size_t i, double shift)
    {
      (&p[i / 3].x)[i % 3] += shift;
      const double phi = Potentials::dihedralAngle(p[0], p[1], p[2], p[3]);
      (&p[i / 3].x)[i % 3] -= shift;
      // unwrap across the branch cut at +-pi
      return phi - std::rint((phi - exact.phi) / (2.0 * std::numbers::pi)) * 2.0 * std::numbers::pi;
    };
    for (std::size_t i = 0; i < 12; ++i)
    {
      const double numerical = (phiAt(i, delta) - phiAt(i, -delta)) / (2.0 * delta);
      EXPECT_NEAR(numerical, (&exact.gradient[i / 3].x)[i % 3], 1e-7) << "gradient i=" << i;

      for (std::size_t j = 0; j < 12; ++j)
      {
        (&p[j / 3].x)[j % 3] += delta;
        const Potentials::DihedralDerivatives plus =
            Potentials::dihedralAngleDerivatives(p[0], p[1], p[2], p[3], false);
        (&p[j / 3].x)[j % 3] -= 2.0 * delta;
        const Potentials::DihedralDerivatives minus =
            Potentials::dihedralAngleDerivatives(p[0], p[1], p[2], p[3], false);
        (&p[j / 3].x)[j % 3] += delta;
        const double second = ((&plus.gradient[i / 3].x)[i % 3] - (&minus.gradient[i / 3].x)[i % 3]) / (2.0 * delta);
        EXPECT_NEAR(second, exact.hessian[i][j], 1e-6) << "hessian i=" << i << " j=" << j;
        EXPECT_NEAR(exact.hessian[i][j], exact.hessian[j][i], 1e-12);
      }
    }
  }
}

TEST(minimization_hessian_intramolecular_pairs, urey_bradley_vdw_coulomb_match_finite_difference)
{
  const double delta = 1e-5;
  const double relativeTolerance = 5e-3;

  ForceField forceField = makeTestForceField();

  ConnectivityTable connectivityTable(3);
  connectivityTable[0, 1] = true;
  connectivityTable[1, 0] = true;
  connectivityTable[1, 2] = true;
  connectivityTable[2, 1] = true;

  Potentials::IntraMolecularPotentials intraMolecularPotentials{};
  intraMolecularPotentials.ureyBradleys = {UreyBradleyPotential({0, 2}, UreyBradleyType::Harmonic, {40000.0, 2.6})};
  intraMolecularPotentials.vanDerWaals = {
      VanDerWaalsPotential({0, 2}, VDWParameters::Type::LennardJones, {120.0, 2.5}, 1.0)};
  intraMolecularPotentials.coulombs = {CoulombPotential({0, 2}, CoulombType::Coulomb, 0.4, -0.35, 1.0)};

  Component component = Component(forceField, "triatomic", 369.825, 4247660.0, 0.1524,
                                  {Atom({-0.3, 1.45, 0.2}, 0.4, 1.0, 0, 1, 0, false, false),
                                   Atom({0.0, 0.0, 0.0}, 0.0, 1.0, 0, 2, 0, false, false),
                                   Atom({1.5, -0.2, 0.4}, -0.35, 1.0, 0, 1, 0, false, false)},
                                  connectivityTable, intraMolecularPotentials, 5, 21);
  component.rigid = false;

  System system = System(forceField, SimulationBox(25.0, 25.0, 25.0), false, 300.0, 1e4, 1.0, {}, {component}, {}, {1},
                         5);

  const MinimizationDofLayout layout = buildMinimizationDofLayout(system.moleculeData, system.components);
  ASSERT_EQ(layout.numDofs(), 9u);

  GeneralizedHessian hessian(layout.numDofs(), 0);
  std::span<Atom> atoms = system.spanOfMoleculeAtoms();
  std::span<AtomDynamics> dynamics = system.spanOfMoleculeDynamics();

  Interactions::computeIntraMolecularUreyBradleyHessian(system.moleculeData, atoms, system.components, layout, hessian,
                                                        dynamics);
  Interactions::computeIntraMolecularVanDerWaalsHessian(system.forceField, system.simulationBox, system.moleculeData,
                                                        atoms, system.components, layout, hessian, dynamics);
  Interactions::computeIntraMolecularCoulombHessian(system.forceField, system.simulationBox, system.moleculeData, atoms,
                                                    system.components, layout, hessian, dynamics);

  const Potentials::IntraMolecularPotentials &potentials = system.components[0].intraMolecularPotentials;
  auto energyFunction = [&]() -> double
  {
    double energy = 0.0;
    for (const UreyBradleyPotential &ureyBradley : potentials.ureyBradleys)
    {
      energy += ureyBradley.calculateEnergy(atoms[ureyBradley.identifiers[0]].position,
                                            atoms[ureyBradley.identifiers[1]].position);
    }
    for (const VanDerWaalsPotential &vanDerWaals : potentials.vanDerWaals)
    {
      energy += vanDerWaals.calculateEnergy(atoms[vanDerWaals.identifiers[0]].position,
                                            atoms[vanDerWaals.identifiers[1]].position);
    }
    for (const CoulombPotential &coulomb : potentials.coulombs)
    {
      energy += coulomb.calculateEnergy(atoms[coulomb.identifiers[0]].position, atoms[coulomb.identifiers[1]].position);
    }
    return energy;
  };

  const std::vector<DofLabel> labels = buildDofLabels(layout, 3);
  compareHessianToFiniteDifference(atoms, labels, hessian, energyFunction, delta, relativeTolerance);
}
