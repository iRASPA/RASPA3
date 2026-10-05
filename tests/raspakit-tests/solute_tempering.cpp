#include <gtest/gtest.h>

import std;

import double3;
import int3;
import units;
import atom;
import forcefield;
import vdwparameters;
import pseudo_atom;
import potential_pair_vdw;
import potential_pair_derivatives;
import component;
import system;
import simulationbox;
import running_energy;
import randomnumbers;
import connectivity_table;
import intra_molecular_potentials;
import bond_potential;
import urey_bradley_potential;
import bend_potential;
import inversion_bend_potential;
import out_of_plane_bend_potential;
import torsion_potential;
import bond_bond_potential;
import bond_bend_potential;
import bond_torsion_potential;
import bend_bend_potential;
import bend_torsion_potential;
import van_der_waals_potential;
import coulomb_potential;
import mc_moves_parallel_tempering_swap;

// Solute tempering (REST2) realizes the per-replica Hamiltonian E_lambda = lambda (U_ss + U_intra) +
// sqrt(lambda) U_sw + U_ww by scaling the energy parameters of the solute pseudo-atom types, the solute
// charges and the solute's intramolecular potentials. That requires every potential to be homogeneous
// of degree one in the parameters that 'scaleEnergy' touches, and the system-level scaling to agree
// with the energy difference the exchange acceptance rule evaluates.

namespace
{
const double3 posA{0.1, 1.3, 0.2};
const double3 posB{0.0, 0.0, 0.0};
const double3 posC{1.5, 0.1, -0.1};
const double3 posD{2.0, 1.2, 0.9};

// generic positive parameter values: (force constant [K], distance [A], angle [deg], ...)
std::vector<double> genericParameters(std::size_t count)
{
  static const std::array<double, 6> values{1200.0, 1.53, 115.0, 2.1, 0.8, 1.3};
  std::vector<double> parameters(count);
  for (std::size_t i = 0; i < count; ++i) parameters[i] = values[i % values.size()];
  return parameters;
}

void expectScaledEnergy(double scaled, double original, double lambda, const std::string& label)
{
  ASSERT_TRUE(std::isfinite(original)) << label;
  ASSERT_TRUE(std::isfinite(scaled)) << label;
  const double tolerance = std::max(1e-12 * std::abs(original), 1e-9);
  EXPECT_NEAR(scaled, lambda * original, tolerance) << label;
}
}  // namespace

// every intramolecular potential type: U(lambda-scaled parameters) == lambda U

TEST(solute_tempering, bond_potentials_are_homogeneous)
{
  const double lambda = 0.37;
  for (std::size_t t = 0; t < BondPotential::numberOfBondParameters.size(); ++t)
  {
    BondPotential potential({0, 1}, static_cast<BondType>(t), genericParameters(BondPotential::numberOfBondParameters[t]));
    BondPotential scaled = potential;
    scaled.scaleEnergy(lambda);
    expectScaledEnergy(scaled.calculateEnergy(posA, posB), potential.calculateEnergy(posA, posB), lambda,
                       std::format("bond type {}", t));
  }
}

TEST(solute_tempering, urey_bradley_potentials_are_homogeneous)
{
  const double lambda = 0.37;
  for (std::size_t t = 0; t < UreyBradleyPotential::numberOfUreyBradleyParameters.size(); ++t)
  {
    UreyBradleyPotential potential({0, 1}, static_cast<UreyBradleyType>(t),
                                   genericParameters(UreyBradleyPotential::numberOfUreyBradleyParameters[t]));
    UreyBradleyPotential scaled = potential;
    scaled.scaleEnergy(lambda);
    expectScaledEnergy(scaled.calculateEnergy(posA, posC), potential.calculateEnergy(posA, posC), lambda,
                       std::format("urey-bradley type {}", t));
  }
}

TEST(solute_tempering, bend_potentials_are_homogeneous)
{
  const double lambda = 0.37;
  for (std::size_t t = 0; t < BendPotential::numberOfBendParameters.size(); ++t)
  {
    BendPotential potential({0, 1, 2}, static_cast<BendType>(t),
                            genericParameters(BendPotential::numberOfBendParameters[t]));
    BendPotential scaled = potential;
    scaled.scaleEnergy(lambda);
    expectScaledEnergy(scaled.calculateEnergy(posA, posB, posC, posD), potential.calculateEnergy(posA, posB, posC, posD),
                       lambda, std::format("bend type {}", t));
  }
}

TEST(solute_tempering, inversion_bend_potentials_are_homogeneous)
{
  const double lambda = 0.37;
  for (std::size_t t = 0; t < InversionBendPotential::numberOfInversionBendParameters.size(); ++t)
  {
    InversionBendPotential potential({0, 1, 2, 3}, static_cast<InversionBendType>(t),
                                     genericParameters(InversionBendPotential::numberOfInversionBendParameters[t]));
    InversionBendPotential scaled = potential;
    scaled.scaleEnergy(lambda);
    expectScaledEnergy(scaled.calculateEnergy(posA, posB, posC, posD), potential.calculateEnergy(posA, posB, posC, posD),
                       lambda, std::format("inversion-bend type {}", t));
  }
}

TEST(solute_tempering, out_of_plane_bend_potentials_are_homogeneous)
{
  const double lambda = 0.37;
  for (std::size_t t = 0; t < OutOfPlaneBendPotential::numberOfOutOfPlaneBendParameters.size(); ++t)
  {
    OutOfPlaneBendPotential potential({0, 1, 2, 3}, static_cast<OutOfPlaneBendType>(t),
                                      genericParameters(OutOfPlaneBendPotential::numberOfOutOfPlaneBendParameters[t]));
    OutOfPlaneBendPotential scaled = potential;
    scaled.scaleEnergy(lambda);
    expectScaledEnergy(scaled.calculateEnergy(posA, posB, posC, posD), potential.calculateEnergy(posA, posB, posC, posD),
                       lambda, std::format("out-of-plane-bend type {}", t));
  }
}

TEST(solute_tempering, torsion_potentials_are_homogeneous)
{
  const double lambda = 0.37;
  for (std::size_t t = 0; t < TorsionPotential::numberOfTorsionParameters.size(); ++t)
  {
    TorsionPotential potential({0, 1, 2, 3}, static_cast<TorsionType>(t),
                               genericParameters(TorsionPotential::numberOfTorsionParameters[t]));
    TorsionPotential scaled = potential;
    scaled.scaleEnergy(lambda);
    expectScaledEnergy(scaled.calculateEnergy(posA, posB, posC, posD), potential.calculateEnergy(posA, posB, posC, posD),
                       lambda, std::format("torsion type {}", t));
  }
}

TEST(solute_tempering, cross_term_potentials_are_homogeneous)
{
  const double lambda = 0.37;
  for (std::size_t t = 0; t < BondBondPotential::numberOfBondBondParameters.size(); ++t)
  {
    BondBondPotential potential({0, 1, 2}, static_cast<BondBondType>(t),
                                genericParameters(BondBondPotential::numberOfBondBondParameters[t]));
    BondBondPotential scaled = potential;
    scaled.scaleEnergy(lambda);
    expectScaledEnergy(scaled.calculateEnergy(posA, posB, posC), potential.calculateEnergy(posA, posB, posC), lambda,
                       std::format("bond-bond type {}", t));
  }
  for (std::size_t t = 0; t < BondBendPotential::numberOfBondBendParameters.size(); ++t)
  {
    BondBendPotential potential({0, 1, 2, 3}, static_cast<BondBendType>(t),
                                genericParameters(BondBendPotential::numberOfBondBendParameters[t]));
    BondBendPotential scaled = potential;
    scaled.scaleEnergy(lambda);
    expectScaledEnergy(scaled.calculateEnergy(posA, posB, posC, posD), potential.calculateEnergy(posA, posB, posC, posD),
                       lambda, std::format("bond-bend type {}", t));
  }
  for (std::size_t t = 0; t < BondTorsionPotential::numberOfBondTorsionParameters.size(); ++t)
  {
    BondTorsionPotential potential({0, 1, 2, 3}, static_cast<BondTorsionType>(t),
                                   genericParameters(BondTorsionPotential::numberOfBondTorsionParameters[t]));
    BondTorsionPotential scaled = potential;
    scaled.scaleEnergy(lambda);
    expectScaledEnergy(scaled.calculateEnergy(posA, posB, posC, posD), potential.calculateEnergy(posA, posB, posC, posD),
                       lambda, std::format("bond-torsion type {}", t));
  }
  for (std::size_t t = 0; t < BendBendPotential::numberOfBendBendParameters.size(); ++t)
  {
    BendBendPotential potential({0, 1, 2, 3}, static_cast<BendBendType>(t),
                                genericParameters(BendBendPotential::numberOfBendBendParameters[t]));
    BendBendPotential scaled = potential;
    scaled.scaleEnergy(lambda);
    expectScaledEnergy(scaled.calculateEnergy(posA, posB, posC, posD), potential.calculateEnergy(posA, posB, posC, posD),
                       lambda, std::format("bend-bend type {}", t));
  }
  for (std::size_t t = 0; t < BendTorsionPotential::numberOfBendTorsionParameters.size(); ++t)
  {
    BendTorsionPotential potential({0, 1, 2, 3}, static_cast<BendTorsionType>(t),
                                   genericParameters(BendTorsionPotential::numberOfBendTorsionParameters[t]));
    BendTorsionPotential scaled = potential;
    scaled.scaleEnergy(lambda);
    expectScaledEnergy(scaled.calculateEnergy(posA, posB, posC, posD), potential.calculateEnergy(posA, posB, posC, posD),
                       lambda, std::format("bend-torsion type {}", t));
  }
}

TEST(solute_tempering, intramolecular_pair_lists_scale_as_lambda)
{
  const double lambda = 0.37;
  Potentials::IntraMolecularPotentials potentials{};
  potentials.vanDerWaals = {VanDerWaalsPotential({0, 3}, VanDerWaalsType::LennardJones, {98.0, 3.75}, 0.5)};
  potentials.coulombs = {CoulombPotential({0, 3}, CoulombType::Coulomb, 0.3, -0.3, 0.5)};

  Potentials::IntraMolecularPotentials scaled = potentials;
  scaled.scaleEnergy(lambda);

  expectScaledEnergy(scaled.vanDerWaals[0].calculateEnergy(posA, posD),
                     potentials.vanDerWaals[0].calculateEnergy(posA, posD), lambda, "intramolecular van der Waals");
  expectScaledEnergy(scaled.coulombs[0].calculateEnergy(posA, posD), potentials.coulombs[0].calculateEnergy(posA, posD),
                     lambda, "intramolecular Coulomb");
}

// every van der Waals potential type: scaling the solute type's table entries multiplies the pair energy
// by lambda (solute-solute), sqrt(lambda) (solute-solvent) and leaves the solvent-solvent pair alone
TEST(solute_tempering, force_field_pair_energies_scale_per_pair_class)
{
  struct Case
  {
    VDWParameters::Type type;
    std::vector<double> values;
    double distance;
  };
  const std::vector<Case> cases{
      {VDWParameters::Type::LennardJones, {119.8, 3.405}, 3.8},
      {VDWParameters::Type::BuckingHam, {3.0e5, 3.5, 1.2e4}, 3.0},
      {VDWParameters::Type::Morse, {500.0, 1.5, 2.0}, 3.0},
      {VDWParameters::Type::MM3, {120.0, 3.5}, 3.9},
      {VDWParameters::Type::BornHugginsMeyer, {6.0e4, 3.15, 2.34, 1.0e6, 2.0e6}, 4.5},
      {VDWParameters::Type::LennardJonesShiftedForce, {119.8, 3.405}, 3.8},
      {VDWParameters::Type::LennardJonesSecondOrderTaylorShifted, {119.8, 3.405}, 3.8},
      {VDWParameters::Type::Potential12_6, {6.0e7, 2.5e4}, 4.0},
      {VDWParameters::Type::Potential12_6_2_0, {6.0e7, -2.5e4, -100.0, 0.0}, 4.0},
      {VDWParameters::Type::CFF9_6, {1.0e6, 2.0e4}, 4.0},
      {VDWParameters::Type::CFFEpsilonSigma, {120.0, 3.4}, 3.8},
      {VDWParameters::Type::MatsuokaClementiYoshimine, {2.0e5, 3.0, -1.0e3, 1.2}, 3.0},
      {VDWParameters::Type::Generic, {3.0e5, 3.2, 1.0e3, 1.0e4, 5.0e4, 1.0e5}, 3.0},
      {VDWParameters::Type::PellenqNicholson, {1.0e5, 3.28, 1.0e4, 5.0e4, 2.0e5}, 2.0},
      {VDWParameters::Type::HydratedIonWater, {1.0e7, 3.5, 1.0e3, 1.0e4, 1.0e6}, 5.0},
      {VDWParameters::Type::Mie, {6.0e7, 12.0, 2.5e4, 6.0}, 4.0},
      {VDWParameters::Type::WeeksChandlerAndersen, {119.8, 3.405}, 3.7},
      {VDWParameters::Type::RepulsiveHarmonic, {500.0, 2.5}, 2.0},
  };

  const double lambda = 0.37;
  for (const Case& c : cases)
  {
    ForceField forceField({PseudoAtom("S", false, 12.0, 0.0, 0.0, 6, false), PseudoAtom("W", false, 12.0, 0.0, 0.0, 6, false)},
                          {VDWParameters(c.type, c.values), VDWParameters(c.type, c.values)},
                          ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true, false, false);
    ForceField scaled = forceField;
    scaled.scaleSoluteInteractions({true, false}, lambda);

    const double rr = c.distance * c.distance;
    const std::string name = VDWParameters::nameOfType(c.type);
    expectScaledEnergy(Potentials::potentialVDW<0>(scaled, 1.0, 1.0, rr, 0, 0).energy,
                       Potentials::potentialVDW<0>(forceField, 1.0, 1.0, rr, 0, 0).energy, lambda, name + " (S,S)");
    expectScaledEnergy(Potentials::potentialVDW<0>(scaled, 1.0, 1.0, rr, 0, 1).energy,
                       Potentials::potentialVDW<0>(forceField, 1.0, 1.0, rr, 0, 1).energy, std::sqrt(lambda),
                       name + " (S,W)");
    expectScaledEnergy(Potentials::potentialVDW<0>(scaled, 1.0, 1.0, rr, 1, 1).energy,
                       Potentials::potentialVDW<0>(forceField, 1.0, 1.0, rr, 1, 1).energy, 1.0, name + " (W,W)");
    EXPECT_NE(Potentials::potentialVDW<0>(forceField, 1.0, 1.0, rr, 0, 0).energy, 0.0) << name;
  }
}

// system level: a charged flexible 4-site chain (the solute, with its own pseudo-atom types) in CO2,
// with Ewald summation and tail corrections

namespace
{
ForceField makeForceField()
{
  ForceField forceField({{"O_co2", false, 15.9994, -0.3256, 0.0, 8, true},
                         {"C_co2", false, 12.0, 0.6512, 0.0, 6, true},
                         {"CH3_s", false, 15.03452, 0.0, 0.0, 6, true},
                         {"CH2_s", false, 14.02658, 0.0, 0.0, 6, true}},
                        {{85.671, 3.017}, {29.933, 2.745}, {98.0, 3.75}, {46.0, 3.95}},
                        ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, false, true, true);
  forceField.automaticEwald = false;
  forceField.EwaldAlpha = 0.25;
  forceField.numberOfWaveVectors = int3(6, 6, 6);
  return forceField;
}

Component makeChain(const ForceField& forceField, std::size_t componentId, std::uint16_t endType, std::uint16_t midType,
                    std::string name)
{
  ConnectivityTable connectivityTable(4);
  connectivityTable[0, 1] = true;
  connectivityTable[1, 0] = true;
  connectivityTable[1, 2] = true;
  connectivityTable[2, 1] = true;
  connectivityTable[2, 3] = true;
  connectivityTable[3, 2] = true;

  Potentials::IntraMolecularPotentials potentials{};
  potentials.bonds = {BondPotential({0, 1}, BondType::Harmonic, {96500.0, 1.54}),
                      BondPotential({1, 2}, BondType::Harmonic, {96500.0, 1.54}),
                      BondPotential({2, 3}, BondType::Harmonic, {96500.0, 1.54})};
  potentials.bends = {BendPotential({0, 1, 2}, BendType::Harmonic, {62500.0, 114.0}),
                      BendPotential({1, 2, 3}, BendType::Harmonic, {62500.0, 114.0})};
  potentials.torsions = {TorsionPotential({0, 1, 2, 3}, TorsionType::TraPPE, {0.0, 355.03, -68.19, 791.32})};
  potentials.vanDerWaals = {VanDerWaalsPotential({0, 3}, VanDerWaalsType::LennardJones, {98.0, 3.75}, 0.5)};
  potentials.coulombs = {CoulombPotential({0, 3}, CoulombType::Coulomb, 0.3, -0.3, 0.5)};

  const std::uint8_t c = static_cast<std::uint8_t>(componentId);
  Component component =
      Component(forceField, std::move(name), 425.125, 3796000.0, 0.201,
                {Atom({-0.5, 1.4, 0.3}, 0.3, 1.0, 0, endType, c, false, false),
                 Atom({0.0, 0.0, 0.0}, -0.3, 1.0, 0, midType, c, false, false),
                 Atom({1.54, 0.0, 0.0}, 0.3, 1.0, 0, midType, c, false, false),
                 Atom({2.1, 1.2, 1.0}, -0.3, 1.0, 0, endType, c, false, false)},
                connectivityTable, potentials, 5, 21);
  component.rigid = false;
  return component;
}

// the chain near the origin, three CO2 molecules around it (within the cut-off) and one far away
void placeMolecules(System& system, double3 shift)
{
  std::span<Atom> atoms = system.spanOfMoleculeAtoms();
  ASSERT_EQ(atoms.size(), 4uz + 4uz * 3uz);

  atoms[0].position = double3(-0.5, 1.4, 0.3) + shift;
  atoms[1].position = double3(0.0, 0.0, 0.0) + shift;
  atoms[2].position = double3(1.54, 0.0, 0.0) + shift;
  atoms[3].position = double3(2.1, 1.2, 1.0) + shift;

  const std::array<double3, 4> centers{double3(4.5, 0.5, 0.0), double3(-3.5, 2.0, 1.0), double3(0.5, -4.0, 2.5),
                                       double3(12.0, 12.0, 12.0)};
  const std::array<double3, 4> axes{double3(0.0, 0.0, 1.0), double3(1.0, 0.0, 0.0), double3(0.6, 0.8, 0.0),
                                    double3(0.0, 1.0, 0.0)};
  for (std::size_t m = 0; m != centers.size(); ++m)
  {
    atoms[4 + 3 * m + 0].position = centers[m] + 1.149 * axes[m];
    atoms[4 + 3 * m + 1].position = centers[m];
    atoms[4 + 3 * m + 2].position = centers[m] - 1.149 * axes[m];
  }
  system.rebuildConfigurationDerivedState();
}

System makeSystem(double3 shift = double3(0.0, 0.0, 0.0))
{
  ForceField forceField = makeForceField();
  Component chain = makeChain(forceField, 0, 2, 3, "chain");
  Component co2 = Component::makeCO2(forceField, 1, true);
  System system = System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e4, 1.0, {}, {chain, co2}, {},
                         {1, 4}, 5);
  placeMolecules(system, shift);
  return system;
}

double potentialEnergy(const System& system) { return system.runningEnergies.potentialEnergy(); }
}  // namespace

TEST(solute_tempering, lambda_one_leaves_the_system_unchanged)
{
  System system = makeSystem();
  System scaled = system;
  scaled.scaleSoluteHamiltonian(0, 1.0);

  EXPECT_DOUBLE_EQ(potentialEnergy(scaled), potentialEnergy(system));
  EXPECT_EQ(scaled.soluteTemperingComponent, std::optional<std::size_t>{0});
  EXPECT_DOUBLE_EQ(scaled.soluteTemperingLambda, 1.0);
  std::span<const Atom> a = system.spanOfMoleculeAtoms();
  std::span<const Atom> b = scaled.spanOfMoleculeAtoms();
  for (std::size_t i = 0; i < a.size(); ++i) EXPECT_DOUBLE_EQ(a[i].charge, b[i].charge);
}

TEST(solute_tempering, scaling_a_shared_pseudo_atom_type_throws)
{
  ForceField forceField = makeForceField();
  Component chainA = makeChain(forceField, 0, 2, 3, "chainA");
  Component chainB = makeChain(forceField, 1, 2, 3, "chainB");  // shares CH3_s / CH2_s with chainA
  System system = System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e4, 1.0, {}, {chainA, chainB}, {},
                         {1, 1}, 5);
  EXPECT_THROW(system.scaleSoluteHamiltonian(0, 0.5), std::runtime_error);
}

// E(lambda) = U_ww + lambda (U_ss + U_intra) + sqrt(lambda) U_sw: the full system energy (real-space
// pair interactions, tail corrections, Ewald Fourier/self/exclusion, intramolecular) must lie on that
// two-term curve; fit the coefficients from three lambda values and predict a fourth
TEST(solute_tempering, total_energy_follows_the_rest2_scaling_law)
{
  System reference = makeSystem();
  const std::array<double, 4> lambdas{1.0, 0.64, 0.36, 0.16};
  std::array<double, 4> energies{};
  for (std::size_t k = 0; k < lambdas.size(); ++k)
  {
    System scaled = reference;
    scaled.scaleSoluteHamiltonian(0, lambdas[k]);
    energies[k] = potentialEnergy(scaled);

    // consistency of the stored running energies with a fresh total-energy evaluation
    RunningEnergy recomputed = scaled.computeTotalEnergies();
    EXPECT_NEAR(recomputed.potentialEnergy(), energies[k], 1e-8 * std::abs(energies[k]) + 1e-8);
  }

  // solve for (c0, c1, c2) in E = c0 + c1 lambda + c2 sqrt(lambda) from the first three points
  auto row = [](double lambda) { return std::array<double, 3>{1.0, lambda, std::sqrt(lambda)}; };
  std::array<std::array<double, 3>, 3> A{row(lambdas[0]), row(lambdas[1]), row(lambdas[2])};
  std::array<double, 3> b{energies[0], energies[1], energies[2]};
  // Gaussian elimination
  for (std::size_t i = 0; i < 3; ++i)
  {
    for (std::size_t j = i + 1; j < 3; ++j)
    {
      const double f = A[j][i] / A[i][i];
      for (std::size_t k = i; k < 3; ++k) A[j][k] -= f * A[i][k];
      b[j] -= f * b[i];
    }
  }
  std::array<double, 3> c{};
  for (std::size_t i = 3; i-- > 0;)
  {
    double s = b[i];
    for (std::size_t k = i + 1; k < 3; ++k) s -= A[i][k] * c[k];
    c[i] = s / A[i][i];
  }
  const auto r = row(lambdas[3]);
  const double predicted = c[0] * r[0] + c[1] * r[1] + c[2] * r[2];
  EXPECT_NEAR(energies[3], predicted, 1e-7 * std::abs(energies[0]) + 1e-7);

  // the solute-solute (+ intramolecular) and solute-solvent parts are not trivially zero in this setup
  EXPECT_GT(std::abs(c[1]), 1.0);
  EXPECT_GT(std::abs(c[2]), 1.0);
}

TEST(solute_tempering, scaling_is_multiplicative)
{
  System once = makeSystem();
  once.scaleSoluteHamiltonian(0, 0.2);

  System twice = makeSystem();
  twice.scaleSoluteHamiltonian(0, 0.5);
  twice.scaleSoluteHamiltonian(0, 0.4);

  EXPECT_NEAR(twice.soluteTemperingLambda, once.soluteTemperingLambda, 1e-15);
  EXPECT_NEAR(potentialEnergy(twice), potentialEnergy(once), 1e-9 * std::abs(potentialEnergy(once)) + 1e-9);
}

// the exchange acceptance evaluates E_target(X) - E_holder(X); for the configuration held by the
// lambda = 1 replica this is the energy change of scaling that very system
TEST(solute_tempering, solute_hamiltonian_change_matches_rescaled_system_energy)
{
  System systemA = makeSystem();
  System systemB = systemA;
  systemB.scaleSoluteHamiltonian(0, 0.3);

  const double expected = potentialEnergy(systemB) - potentialEnergy(systemA);
  const std::optional<double> changeAToB = MC_Moves::ParallelTemperingSoluteHamiltonianChange(systemA, systemB);
  const std::optional<double> changeBToA = MC_Moves::ParallelTemperingSoluteHamiltonianChange(systemB, systemA);
  ASSERT_TRUE(changeAToB.has_value());
  ASSERT_TRUE(changeBToA.has_value());
  EXPECT_NEAR(changeAToB.value(), expected, 1e-8 * std::abs(expected) + 1e-8);
  EXPECT_NEAR(changeBToA.value(), -expected, 1e-8 * std::abs(expected) + 1e-8);

  // identical configurations: the acceptance rule reduces to zero (equal temperatures)
  const std::optional<double> logR = MC_Moves::ParallelTemperingLogAcceptance(systemA, systemB);
  ASSERT_TRUE(logR.has_value());
  EXPECT_NEAR(logR.value(), 0.0, 1e-8);
}

// different configurations: log R = -beta [ (E_A(X_B) - E_B(X_B)) + (E_B(X_A) - E_A(X_A)) ] at equal T,
// with E_A(X_B) obtained by brute force (rescaling a copy of B back to A's Hamiltonian)
TEST(solute_tempering, log_acceptance_matches_brute_force_hamiltonian_switch)
{
  System systemA = makeSystem();
  System systemB = makeSystem(double3(0.7, -0.4, 1.1));
  systemB.scaleSoluteHamiltonian(0, 0.3);

  System aAtB = systemB;  // configuration X_B, Hamiltonian A
  aAtB.scaleSoluteHamiltonian(0, 1.0 / 0.3);
  System bAtA = systemA;  // configuration X_A, Hamiltonian B
  bAtA.scaleSoluteHamiltonian(0, 0.3);

  const double deltaAofB = potentialEnergy(aAtB) - potentialEnergy(systemB);
  const double deltaBofA = potentialEnergy(bAtA) - potentialEnergy(systemA);
  const double expected = -systemA.beta * deltaAofB - systemB.beta * deltaBofA;

  const std::optional<double> logR = MC_Moves::ParallelTemperingLogAcceptance(systemA, systemB);
  ASSERT_TRUE(logR.has_value());
  EXPECT_NEAR(logR.value(), expected, 1e-7 * std::abs(expected) + 1e-7);
  EXPECT_GT(std::abs(expected), 1e-3);
}

// after an accepted exchange the configurations have travelled but each replica keeps its Hamiltonian:
// the solute charges carry the receiving replica's scaling and the running energies are consistent
TEST(solute_tempering, swap_keeps_the_hamiltonian_with_the_replica)
{
  // nearly identical configurations so that the exchange has an appreciable acceptance probability
  System systemA = makeSystem();
  System systemB = makeSystem(double3(0.01, -0.005, 0.02));
  systemB.scaleSoluteHamiltonian(0, 0.3);

  const double energyAofXB = [&] {
    System s = systemB;
    s.scaleSoluteHamiltonian(0, 1.0 / 0.3);
    return potentialEnergy(s);
  }();
  const double energyBofXA = [&] {
    System s = systemA;
    s.scaleSoluteHamiltonian(0, 0.3);
    return potentialEnergy(s);
  }();
  const double3 positionB0 = systemB.spanOfMoleculeAtoms()[0].position;

  bool accepted = false;
  for (std::size_t seed = 1; seed < 500 && !accepted; ++seed)
  {
    RandomNumber random(seed);
    accepted = MC_Moves::ParallelTemperingSwap(random, systemA, systemB).has_value();
  }
  ASSERT_TRUE(accepted);

  EXPECT_DOUBLE_EQ(systemA.soluteTemperingLambda, 1.0);
  EXPECT_DOUBLE_EQ(systemB.soluteTemperingLambda, 0.3);
  EXPECT_EQ(systemA.spanOfMoleculeAtoms()[0].position.x, positionB0.x);

  std::span<const Atom> atomsA = systemA.spanOfMoleculeAtoms();
  std::span<const Atom> atomsB = systemB.spanOfMoleculeAtoms();
  EXPECT_NEAR(atomsA[0].charge, 0.3, 1e-14);
  EXPECT_NEAR(atomsB[0].charge, 0.3 * std::sqrt(0.3), 1e-14);

  EXPECT_NEAR(potentialEnergy(systemA), energyAofXB, 1e-8 * std::abs(energyAofXB) + 1e-8);
  EXPECT_NEAR(potentialEnergy(systemB), energyBofXA, 1e-8 * std::abs(energyBofXA) + 1e-8);
  EXPECT_NEAR(systemA.computeTotalEnergies().potentialEnergy(), potentialEnergy(systemA),
              1e-8 * std::abs(energyAofXB) + 1e-8);
  EXPECT_NEAR(systemB.computeTotalEnergies().potentialEnergy(), potentialEnergy(systemB),
              1e-8 * std::abs(energyBofXA) + 1e-8);
}
