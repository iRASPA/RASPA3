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
import running_energy;
import intra_molecular_potentials;
import bond_potential;
import bend_potential;
import torsion_potential;
import van_der_waals_potential;
import coulomb_potential;
import integrators;
import integrators_update;
import integrators_compute;

// The MD integrator obtains the intramolecular forces from
// IntraMolecularPotentials::computeInternalGradient (per-potential potentialEnergyGradientStrain),
// whereas Monte Carlo uses computeInternalEnergies (per-potential calculateEnergy). These tests
// check that both paths describe the same potential for every flexible term, including the 1-4
// scaled Coulomb and the 1-5+ Lennard-Jones/Coulomb pairs that the bond/bend/torsion tests do
// not cover, and that the resulting Velocity-Verlet dynamics conserve energy.
namespace
{
ForceField makeChargedChainForceField()
{
  return ForceField({{"A", false, 14.0, 0.3, 0.0, 6, false}, {"B", false, 16.0, -0.3, 0.0, 8, false}},
                    {{80.0, 3.5}, {50.0, 3.0}}, ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, false,
                    false, true);
}

// A six-bead A-B-A-B-A-B chain (net charge zero) with harmonic bonds and bends, a TraPPE-extended
// torsion series with a non-zero cos(4 phi) term, Lennard-Jones on the 1-5 and 1-6 pairs, and
// Coulomb on the 1-4 (scaled by 0.5), 1-5 and 1-6 pairs.
Component makeChargedChain(const ForceField& forceField)
{
  ConnectivityTable connectivityTable(6);
  for (std::size_t i = 0; i < 5; ++i)
  {
    connectivityTable[i, i + 1] = true;
    connectivityTable[i + 1, i] = true;
  }

  Potentials::IntraMolecularPotentials potentials{};
  for (std::size_t i = 0; i < 5; ++i)
  {
    potentials.bonds.push_back(BondPotential({i, i + 1}, BondType::Harmonic, {96500.0, 1.5}));
  }
  for (std::size_t i = 0; i < 4; ++i)
  {
    potentials.bends.push_back(BendPotential({i, i + 1, i + 2}, BendType::Harmonic, {62500.0, 112.0}));
  }
  for (std::size_t i = 0; i < 3; ++i)
  {
    potentials.torsions.push_back(TorsionPotential({i, i + 1, i + 2, i + 3}, TorsionType::TraPPE_Extended,
                                                   {2029.99, -751.83, -538.95, -22.10, -51.27}));
  }

  const std::array<double, 6> charges{0.3, -0.3, 0.3, -0.3, 0.3, -0.3};
  auto lennardJones = [&](std::size_t a, std::size_t b)
  {
    // Lorentz-Berthelot of the two pseudo-atom types (A even, B odd), in Kelvin and Angstrom.
    const double epsilonA = (a % 2 == 0) ? 80.0 : 50.0;
    const double epsilonB = (b % 2 == 0) ? 80.0 : 50.0;
    const double sigmaA = (a % 2 == 0) ? 3.5 : 3.0;
    const double sigmaB = (b % 2 == 0) ? 3.5 : 3.0;
    return VanDerWaalsPotential({a, b}, VanDerWaalsType::LennardJones,
                                {std::sqrt(epsilonA * epsilonB), 0.5 * (sigmaA + sigmaB)}, 1.0);
  };
  potentials.vanDerWaals = {lennardJones(0, 4), lennardJones(1, 5), lennardJones(0, 5)};
  potentials.coulombs = {CoulombPotential({0, 3}, CoulombType::Coulomb, charges[0], charges[3], 0.5),
                         CoulombPotential({1, 4}, CoulombType::Coulomb, charges[1], charges[4], 0.5),
                         CoulombPotential({2, 5}, CoulombType::Coulomb, charges[2], charges[5], 0.5),
                         CoulombPotential({0, 4}, CoulombType::Coulomb, charges[0], charges[4], 1.0),
                         CoulombPotential({1, 5}, CoulombType::Coulomb, charges[1], charges[5], 1.0),
                         CoulombPotential({0, 5}, CoulombType::Coulomb, charges[0], charges[5], 1.0)};

  std::vector<Atom> atoms{};
  for (std::size_t i = 0; i < 6; ++i)
  {
    atoms.push_back(Atom({0.0, 0.0, 0.0}, charges[i], 1.0, 0, static_cast<std::uint16_t>(i % 2), 0, false, false));
  }

  Component component = Component(forceField, "charged-chain", 500.0, 3.0e6, 0.3, atoms, connectivityTable,
                                  potentials, 5, 21);
  component.rigid = false;
  return component;
}

// A non-planar, strained conformation so that every term exerts a force.
std::vector<double3> strainedConformation()
{
  const double3 origin(10.0, 10.0, 10.0);
  return {origin + double3(0.00, 0.00, 0.00), origin + double3(1.55, 0.00, 0.00), origin + double3(2.05, 1.40, 0.10),
          origin + double3(3.50, 1.55, 0.35), origin + double3(4.05, 2.85, -0.40), origin + double3(5.40, 3.10, 0.25)};
}

double3 finiteDifferenceGradient(const Potentials::IntraMolecularPotentials& potentials, std::span<Atom> atoms,
                                 std::size_t index, double delta)
{
  const double3 saved = atoms[index].position;
  double3 gradient{};
  for (std::size_t k = 0; k < 3; ++k)
  {
    double3 step{};
    (k == 0 ? step.x : k == 1 ? step.y : step.z) = delta;
    atoms[index].position = saved + step;
    const double plus = potentials.computeInternalEnergies(atoms).potentialEnergy();
    atoms[index].position = saved - step;
    const double minus = potentials.computeInternalEnergies(atoms).potentialEnergy();
    atoms[index].position = saved;
    (k == 0 ? gradient.x : k == 1 ? gradient.y : gradient.z) = (plus - minus) / (2.0 * delta);
  }
  return gradient;
}

void expectGradientMatchesFiniteDifference(const Potentials::IntraMolecularPotentials& potentials,
                                           std::span<Atom> atoms, std::span<AtomDynamics> dynamics,
                                           std::string_view term)
{
  const double delta = 1.0e-5;

  for (AtomDynamics& dynamic : dynamics) dynamic.gradient = double3(0.0, 0.0, 0.0);
  const RunningEnergy gradientPathEnergy = potentials.computeInternalGradient(atoms, dynamics);
  const RunningEnergy energyPathEnergy = potentials.computeInternalEnergies(atoms);

  EXPECT_NEAR(gradientPathEnergy.potentialEnergy(), energyPathEnergy.potentialEnergy(),
              1.0e-9 * std::max(1.0, std::abs(energyPathEnergy.potentialEnergy())))
      << term << ": gradient and energy paths disagree on the energy";

  double maximumGradient = 0.0;
  for (std::size_t i = 0; i < atoms.size(); ++i)
  {
    maximumGradient = std::max({maximumGradient, std::abs(dynamics[i].gradient.x), std::abs(dynamics[i].gradient.y),
                                std::abs(dynamics[i].gradient.z)});
  }
  ASSERT_GT(maximumGradient, 1.0) << term << ": test conformation exerts no force";

  const double tolerance = 1.0e-5 * (1.0 + maximumGradient);
  for (std::size_t i = 0; i < atoms.size(); ++i)
  {
    const double3 numerical = finiteDifferenceGradient(potentials, atoms, i, delta);
    EXPECT_NEAR(dynamics[i].gradient.x, numerical.x, tolerance) << term << ": atom " << i << " x";
    EXPECT_NEAR(dynamics[i].gradient.y, numerical.y, tolerance) << term << ": atom " << i << " y";
    EXPECT_NEAR(dynamics[i].gradient.z, numerical.z, tolerance) << term << ": atom " << i << " z";
  }
}
}  // namespace

TEST(MC_intramolecular_gradient, flexible_charged_chain_per_term_matches_finite_difference)
{
  ForceField forceField = makeChargedChainForceField();
  Component component = makeChargedChain(forceField);
  System system = System(forceField, SimulationBox(25.0, 25.0, 25.0), false, 300.0, 1e4, 1.0, {}, {component},
                         {strainedConformation()}, {0}, 5);

  std::span<Atom> atoms = system.spanOfMoleculeAtoms();
  std::span<AtomDynamics> dynamics = system.spanOfMoleculeDynamics();
  ASSERT_EQ(atoms.size(), 6uz);
  const Potentials::IntraMolecularPotentials& full = system.components[0].intraMolecularPotentials;
  ASSERT_EQ(full.vanDerWaals.size(), 3uz);
  ASSERT_EQ(full.coulombs.size(), 6uz);

  Potentials::IntraMolecularPotentials bonds{};
  bonds.bonds = full.bonds;
  expectGradientMatchesFiniteDifference(bonds, atoms, dynamics, "bonds");

  Potentials::IntraMolecularPotentials bends{};
  bends.bends = full.bends;
  expectGradientMatchesFiniteDifference(bends, atoms, dynamics, "bends");

  Potentials::IntraMolecularPotentials torsions{};
  torsions.torsions = full.torsions;
  expectGradientMatchesFiniteDifference(torsions, atoms, dynamics, "torsions");

  Potentials::IntraMolecularPotentials vanDerWaals{};
  vanDerWaals.vanDerWaals = full.vanDerWaals;
  expectGradientMatchesFiniteDifference(vanDerWaals, atoms, dynamics, "intramolecular van der Waals");

  Potentials::IntraMolecularPotentials coulombs{};
  coulombs.coulombs = full.coulombs;
  expectGradientMatchesFiniteDifference(coulombs, atoms, dynamics, "intramolecular Coulomb");

  expectGradientMatchesFiniteDifference(full, atoms, dynamics, "all terms");
}

TEST(integrators_flexible_adsorbate, charged_flexible_chain_nve_conserves_energy)
{
  ForceField forceField = makeChargedChainForceField();
  Component component = makeChargedChain(forceField);
  System system = System(forceField, SimulationBox(25.0, 25.0, 25.0), false, 300.0, 1e4, 1.0, {}, {component},
                         {strainedConformation()}, {0}, 5);
  system.timeStep = 1.0e-4;

  std::span<AtomDynamics> dynamics = system.spanOfMoleculeDynamics();
  dynamics[0].velocity = {0.05, -0.02, 0.01};
  dynamics[1].velocity = {-0.04, 0.03, -0.01};
  dynamics[2].velocity = {0.02, 0.01, 0.03};
  dynamics[3].velocity = {-0.01, -0.04, 0.02};
  dynamics[4].velocity = {0.03, 0.02, -0.05};
  dynamics[5].velocity = {-0.05, 0.00, 0.00};

  RunningEnergy initial = Integrators::updateGradients(
      system.moleculeData, system.spanOfMoleculeAtoms(), system.spanOfMoleculeDynamics(),
      system.spanOfFrameworkAtoms(), system.forceField, system.simulationBox, system.components, system.eik_x,
      system.eik_y, system.eik_z, system.eik_xy, system.trialEik, system.fixedFrameworkStoredEik,
      system.interpolationGrids, system.numberOfMoleculesPerComponent, system.framework,
      system.spanOfFrameworkDynamics());
  initial.translationalKineticEnergy = Integrators::computeTranslationalKineticEnergy(
      system.moleculeData, system.spanOfMoleculeAtoms(), system.spanOfMoleculeDynamics(), system.components,
      system.framework, system.spanOfFrameworkAtoms(), system.spanOfFrameworkDynamics(), &system.forceField);
  const double initialEnergy = initial.conservedEnergy();
  const double3 initialPosition = system.spanOfMoleculeAtoms()[0].position;

  double maximumDrift = 0.0;
  for (std::size_t step = 0; step < 400; ++step)
  {
    const RunningEnergy current = Integrators::velocityVerlet(
        system.moleculeData, system.spanOfMoleculeAtoms(), system.spanOfMoleculeDynamics(), system.components,
        system.timeStep, system.thermostat, system.spanOfFrameworkAtoms(), system.forceField, system.simulationBox,
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.trialEik, system.fixedFrameworkStoredEik,
        system.interpolationGrids, system.numberOfMoleculesPerComponent, system.framework,
        system.spanOfFrameworkDynamics());
    maximumDrift = std::max(maximumDrift, std::abs(current.conservedEnergy() - initialEnergy));
  }

  // The chain must actually have moved (bond and pair forces are large in the strained conformation).
  EXPECT_GT((system.spanOfMoleculeAtoms()[0].position - initialPosition).length(), 1.0e-3);
  EXPECT_LT(maximumDrift / std::max(1.0, std::abs(initialEnergy)), 1.0e-5);
}
