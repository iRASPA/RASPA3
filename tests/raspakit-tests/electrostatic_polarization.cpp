#include <gtest/gtest.h>

import std;

import int3;
import double3;
import double3x3;
import units;
import atom;
import molecule;
import pseudo_atom;
import vdwparameters;
import forcefield;
import framework;
import component;
import system;
import simulationbox;
import running_energy;
import randomnumbers;
import interactions_intermolecular;
import interactions_framework_molecule;
import interactions_ewald;
import energy_status;
import mc_moves_translation;
import mc_moves_rotation;
import mc_moves_gibbs_identity_change;
import mc_moves_lambda_exchange;

TEST(electrostatic_polarization, Test_2_CO2_in_ITQ_29_2x2x2)
{
  ForceField forceField = ForceField::makeZeoliteForceField(11.8, true, false, true);

  forceField.computePolarization = true;
  forceField.omitInterPolarization = false;
  forceField.omitInterInteractions = false;
  forceField.automaticEwald = false;
  forceField.EwaldAlpha = 0.25;
  forceField.numberOfWaveVectors = int3(8, 8, 8);

  Framework f = Framework::makeITQ29(forceField, int3(2, 2, 2));
  Component c = Component::makeCO2(forceField, 0, true);

  System system = System(forceField, std::nullopt, false, 300.0, 1e4, 1.0, {f}, {c}, {}, {2}, 5);

  std::span<Atom> spanOfMoleculeAtoms = system.spanOfMoleculeAtoms();
  spanOfMoleculeAtoms[0].position = double3(5.93355, 7.93355, 5.93355 + 1.149);
  spanOfMoleculeAtoms[0].scalingCoulomb = 1.0;
  spanOfMoleculeAtoms[1].position = double3(5.93355, 7.93355, 5.93355 + 0.0);
  spanOfMoleculeAtoms[1].scalingCoulomb = 1.0;
  spanOfMoleculeAtoms[2].position = double3(5.93355, 7.93355, 5.93355 - 1.149);
  spanOfMoleculeAtoms[2].scalingCoulomb = 1.0;

  spanOfMoleculeAtoms[3].position = double3(5.93355, 3.93355, 5.93355 + 1.149);
  spanOfMoleculeAtoms[3].scalingCoulomb = 1.0;
  spanOfMoleculeAtoms[4].position = double3(5.93355, 3.93355, 5.93355 + 0.0);
  spanOfMoleculeAtoms[4].scalingCoulomb = 1.0;
  spanOfMoleculeAtoms[5].position = double3(5.93355, 3.93355, 5.93355 - 1.149);
  spanOfMoleculeAtoms[5].scalingCoulomb = 1.0;

  RunningEnergy energy = system.computeTotalEnergies();

  // Total polarization energy for the framework + reciprocal + molecule-molecule (real-space) field model.
  // The intramolecular reciprocal-exclusion field is intentionally not part of the field (the reciprocal field
  // is built from the fixed framework structure factor only), and the molecule-molecule contribution is included
  // because omitInterPolarization == false.
  EXPECT_NEAR(energy.polarization * Units::EnergyToKelvin, -1.3034969457316241, 1e-6);
}

static double maxFieldDifference(std::span<const double3> a, std::span<const double3> b)
{
  double m = 0.0;
  for (std::size_t i = 0; i < a.size(); ++i)
  {
    m = std::max(m, (a[i] - b[i]).length());
  }
  return m;
}

// Lambda exchange (Hamiltonian parallel tempering of pinned fractional molecules) for an ion pair, i.e. two
// fractional molecules rescaled in place one after the other. Both replicas are polarizable: the energy, the
// stored field and the dU/dlambda accumulators have to agree with a full recomputation after the exchange. The
// pinned bins straddle lambda = 1/2, where the Coulomb coupling of the fractional molecules switches on, so the
// total-energy field routines and the move differences must apply the same skip of Coulomb-decoupled pairs.
namespace
{
System makePinnedPairSystem(std::size_t fixedBin, std::uint64_t seed, bool polarization)
{
  ForceField forceField =
      ForceField({{"Na+", false, 22.99, 1.0, 1.5, 11, false}, {"Cl-", false, 35.45, -1.0, 3.0, 17, false}},
                 {{15.0966, 2.65755}, {142.562, 3.51932}}, ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0,
                 12.0, true, false, true);
  forceField.computePolarization = polarization;
  forceField.omitInterPolarization = false;

  Component na = Component::makeIon(forceField, 0, "Na+", forceField.findPseudoAtom("Na+").value(), 1.0);
  na.pairComponentId = 1;
  na.maximumPairDistance = 8.0;
  na.fixedLambdaBin = fixedBin;
  na.fixedLambdaCoordinate = Component::FixedLambdaCoordinate::PairSwap;
  na.lambdaPairSwap.computeDUdlambda = true;

  Component cl = Component::makeIon(forceField, 1, "Cl-", forceField.findPseudoAtom("Cl-").value(), -1.0);
  cl.pairComponentId = 0;
  cl.maximumPairDistance = 8.0;

  // very high temperature so that the exchange is accepted regardless of the energy change
  System system = System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 1.0e8, 1e5, 1.0, {}, {na, cl}, {},
                         {6, 6}, 5);
  RandomNumber random(seed);
  for (Atom& atom : system.spanOfMoleculeAtoms())
  {
    atom.position = double3(30.0 * random.uniform(), 30.0 * random.uniform(), 30.0 * random.uniform());
  }
  system.runningEnergies = system.computeTotalEnergies();
  return system;
}

void expectLambdaExchangeConsistent(bool polarization)
{
  System systemA = makePinnedPairSystem(3, 7, polarization);
  System systemB = makePinnedPairSystem(15, 11, polarization);

  const double lambdaA = systemA.components[0].lambdaPairSwap.lambdaValue();
  const double lambdaB = systemB.components[0].lambdaPairSwap.lambdaValue();
  ASSERT_LT(lambdaA, 0.5);
  ASSERT_GT(lambdaB, 0.5);
  ASSERT_EQ(systemA.numberOfPairSwapFractionalMoleculesPerComponent_CFCMC[0], 1uz);
  ASSERT_EQ(systemA.numberOfPairSwapFractionalMoleculesPerComponent_CFCMC[1], 1uz);

  RandomNumber random(42);
  std::optional<std::pair<RunningEnergy, RunningEnergy>> result = MC_Moves::LambdaExchange(random, systemA, systemB);
  ASSERT_TRUE(result.has_value());

  EXPECT_DOUBLE_EQ(systemA.components[0].lambdaPairSwap.lambdaValue(), lambdaB);
  EXPECT_DOUBLE_EQ(systemB.components[0].lambdaPairSwap.lambdaValue(), lambdaA);

  for (System* s : {&systemA, &systemB})
  {
    std::span<double3> maintained = s->spanOfMoleculeElectricField();
    std::vector<double3> snapshot(maintained.begin(), maintained.end());

    RunningEnergy recomputed = s->computeTotalEnergies();
    RunningEnergy drift = s->runningEnergies - recomputed;

    if (polarization)
    {
      EXPECT_LT(maxFieldDifference(snapshot, s->spanOfMoleculeElectricField()), 1e-8);
      EXPECT_GT(std::abs(recomputed.polarization), 0.0);
    }
    EXPECT_NEAR(drift.polarization, 0.0, 1e-6);
    EXPECT_NEAR(drift.moleculeMoleculeVDW, 0.0, 1e-6);
    EXPECT_NEAR(drift.moleculeMoleculeCharge, 0.0, 1e-6);
    EXPECT_NEAR(drift.ewald_fourier, 0.0, 1e-6);
    EXPECT_NEAR(drift.potentialEnergy(), 0.0, 1e-6);
    EXPECT_NEAR(drift.totalDudlambdaVDW(), 0.0, 1e-6);
    EXPECT_NEAR(drift.totalDudlambdaCharge(), 0.0, 1e-6);
    EXPECT_NEAR(drift.totalDudlambdaEwald(), 0.0, 1e-6);
  }
}
}  // namespace

TEST(electrostatic_polarization, lambda_exchange_pair_cfcmc_without_polarization)
{
  expectLambdaExchangeConsistent(false);
}

TEST(electrostatic_polarization, lambda_exchange_pair_cfcmc_with_polarization)
{
  expectLambdaExchangeConsistent(true);
}

TEST(electrostatic_polarization, stored_field_consistency_translation)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  forceField.computePolarization = true;
  forceField.omitInterPolarization = false;
  forceField.omitInterInteractions = false;
  forceField.omitEwaldFourier = false;

  Framework f = Framework::makeITQ29(forceField, int3(2, 2, 2));
  Component c = Component::makeCO2(forceField, 0, true);
  System system = System(forceField, std::nullopt, false, 300.0, 1e4, 1.0, {f}, {c}, {}, {10}, 5);

  RandomNumber random(42);

  RunningEnergy running = system.computeTotalEnergies();

  for (std::size_t step = 0; step < 400; ++step)
  {
    std::size_t selectedMolecule = system.randomMoleculeOfComponent(random, 0);

    std::optional<RunningEnergy> e =
        MC_Moves::translationMove(random, system, 0, selectedMolecule);
    if (e.has_value()) running += e.value();
  }

  // snapshot the incrementally-maintained field, then rebuild it from scratch and compare
  std::span<double3> maintained = system.spanOfMoleculeElectricField();
  std::vector<double3> snapshot(maintained.begin(), maintained.end());

  system.computeTotalElectricField();
  std::span<double3> rebuilt = system.spanOfMoleculeElectricField();

  EXPECT_LT(maxFieldDifference(snapshot, rebuilt), 1e-8);

  RunningEnergy recomputed = system.computeTotalEnergies();
  EXPECT_NEAR(running.polarization - recomputed.polarization, 0.0, 1e-6);
  EXPECT_NEAR(running.potentialEnergy() - recomputed.potentialEnergy(), 0.0, 1e-6);
}

// Thermodynamic integration: the instantaneous dU/dlambda of a fractional molecule must match the
// numerical lambda-derivative of the total potential energy, including the polarization part (the
// polarization coupling is scaled by the atom's Coulomb scaling, and the field the fractional
// molecule produces at other molecules is scaled as well).
TEST(electrostatic_polarization, dudlambda_matches_finite_difference)
{
  ForceField forceField = ForceField::makeZeoliteForceField(11.8, true, false, true);
  forceField.computePolarization = true;
  forceField.omitInterPolarization = false;
  forceField.omitInterInteractions = false;
  forceField.automaticEwald = false;
  forceField.EwaldAlpha = 0.25;
  forceField.numberOfWaveVectors = int3(8, 8, 8);

  Framework f = Framework::makeITQ29(forceField, int3(2, 2, 2));
  Component c = Component::makeCO2(forceField, 0, true);
  System system = System(forceField, std::nullopt, false, 300.0, 1e4, 1.0, {f}, {c}, {}, {6}, 5);

  const double lambda = 0.7;
  std::span<Atom> atoms = system.spanOfMoleculeAtoms();
  for (std::size_t i = 0; i < 3; ++i)
  {
    atoms[i].setScalingToFractional(lambda, 1);
  }

  system.runningEnergies = system.computeTotalEnergies();
  double analytic = system.currentDUdlambda(lambda, 1) * Units::EnergyToKelvin;

  const double delta = 1e-5;
  for (std::size_t i = 0; i < 3; ++i) atoms[i].setScaling(lambda + 0.5 * delta);
  double forward = system.computeTotalEnergies().potentialEnergy();
  for (std::size_t i = 0; i < 3; ++i) atoms[i].setScaling(lambda - 0.5 * delta);
  double backward = system.computeTotalEnergies().potentialEnergy();

  double numeric = ((forward - backward) / delta) * Units::EnergyToKelvin;

  EXPECT_NEAR(analytic, numeric, std::max(1e-4, 1e-6 * std::abs(numeric)));

  // The polarization contribution itself must be non-trivial, otherwise this test would not
  // distinguish the polarization-aware dU/dlambda from the plain VDW/charge/Ewald bookkeeping.
  for (std::size_t i = 0; i < 3; ++i) atoms[i].setScaling(lambda);
  system.runningEnergies = system.computeTotalEnergies();
  double withoutPolarization = system.runningEnergies.dudlambda(lambda, 1) * Units::EnergyToKelvin;
  EXPECT_GT(std::abs(analytic - withoutPolarization), 0.1);  // ~1.19 K for this configuration
}

// The Gibbs identity change switches a molecule in each box; with molecule-molecule polarization the real-space
// field the switched molecule feels from its neighbours and the field change it causes on them must be part of
// the energy and of the stored field, exactly as in the single-box identity change.
TEST(electrostatic_polarization, stored_field_consistency_gibbs_identity_change)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  forceField.computePolarization = true;
  forceField.omitInterPolarization = false;
  forceField.omitInterInteractions = false;
  forceField.omitEwaldFourier = false;

  Framework f = Framework::makeITQ29(forceField, int3(2, 2, 2));
  Component a = Component::makeCO2(forceField, 0, true);
  Component b = Component::makeCO2(forceField, 1, true);
  a.gibbsIdentityChanges = {1};
  b.gibbsIdentityChanges = {0};
  ASSERT_EQ(a.type, b.type);

  System systemI = System(forceField, std::nullopt, false, 300.0, 1e4, 1.0, {f}, {a, b}, {}, {6, 4}, 5);
  System systemII = System(forceField, std::nullopt, false, 300.0, 1e4, 1.0, {f}, {a, b}, {}, {4, 6}, 5);

  RandomNumber random(17);
  RunningEnergy runningI = systemI.computeTotalEnergies();
  RunningEnergy runningII = systemII.computeTotalEnergies();
  systemI.runningEnergies = runningI;
  systemII.runningEnergies = runningII;

  std::size_t accepted = 0;
  for (std::size_t step = 0; step < 300; ++step)
  {
    const std::size_t componentA = random.uniform_integer(0, 1);
    std::optional<std::pair<RunningEnergy, RunningEnergy>> e =
        MC_Moves::GibbsIdentityChangeMove_CBMC(random, systemI, systemII, componentA);
    if (e.has_value())
    {
      ++accepted;
      runningI += e->first;
      runningII += e->second;
    }
  }
  ASSERT_GE(accepted, 5uz);

  for (System* system : {&systemI, &systemII})
  {
    std::span<double3> maintained = system->spanOfMoleculeElectricField();
    std::vector<double3> snapshot(maintained.begin(), maintained.end());
    system->computeTotalElectricField();
    EXPECT_LT(maxFieldDifference(snapshot, system->spanOfMoleculeElectricField()), 1e-8);
  }

  RunningEnergy recomputedI = systemI.computeTotalEnergies();
  RunningEnergy recomputedII = systemII.computeTotalEnergies();
  EXPECT_NEAR(runningI.polarization - recomputedI.polarization, 0.0, 1e-6);
  EXPECT_NEAR(runningI.potentialEnergy() - recomputedI.potentialEnergy(), 0.0, 1e-6);
  EXPECT_NEAR(runningII.polarization - recomputedII.polarization, 0.0, 1e-6);
  EXPECT_NEAR(runningII.potentialEnergy() - recomputedII.potentialEnergy(), 0.0, 1e-6);
}
