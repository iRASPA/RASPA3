#include <gtest/gtest.h>

import std;

import int3;
import double3;
import atom;
import forcefield;
import framework;
import component;
import system;
import simulationbox;
import running_energy;
import monte_carlo;
import equation_of_states;
import isotherm_bet;
import dctst_wham;
import mc_moves_move_types;
import mc_moves_probabilities;
import mc_moves_acceptance;
import mc_moves_parallel_tempering_swap;

// Regression tests: the overlap test of the dcTST WHAM solver, the samples of a state in a bin it forbids, the
// fractional molecule that the CFCMC and CB/CFCMC swaps share, the roots of the cubic below the co-volume, the
// replica-exchange swap under an external field, the biased acceptance in log space, the tracked adsorbate and
// the occupancy scout.

namespace
{

Component makeCO2(const ForceField& forceField, const MCMoveProbabilities& probabilities)
{
  return Component(forceField, "CO2", 304.1282, 7377300.0, 0.22394,
                   {Atom({0, 0, 1.149}, -0.3256, 1.0, 0, 4, 0, false, false),
                    Atom({0, 0, 0.000}, 0.6512, 1.0, 0, 3, 0, false, false),
                    Atom({0, 0, -1.149}, -0.3256, 1.0, 0, 4, 0, false, false)},
                   {}, {}, 5, 21, probabilities, std::nullopt, false);
}

}  // namespace

// Two states with disjoint histograms and a finite bias on each other's bin. WHAM knows the samples of a bin by
// their sum only and has a fixed point here, g = {0, 0}, on which it converges at once; each state is credited with
// the share e^-10 of the samples of the other, far above the old threshold of 1e-12.
TEST(WHAM_OVERLAP_AND_SHARED_SLOT_FIXES, dctst_wham_needs_a_bin_that_both_states_sampled)
{
  const DcTSTWham::Result disjoint =
      DcTSTWham::solve({{1000.0, 0.0}, {0.0, 1000.0}}, {{0.0, 10.0}, {10.0, 0.0}}, 1e-10, 10000);
  EXPECT_FALSE(disjoint.connected);
  EXPECT_FALSE(disjoint.converged);

  // one shared bin measures the difference
  const DcTSTWham::Result overlapping =
      DcTSTWham::solve({{1000.0, 3.0}, {2.0, 1000.0}}, {{0.0, 10.0}, {10.0, 0.0}}, 1e-10, 10000);
  EXPECT_TRUE(overlapping.connected);
  EXPECT_TRUE(overlapping.converged);

  // a chain of three: the outer states are linked through the middle one
  const DcTSTWham::Result chain =
      DcTSTWham::solve({{100.0, 50.0, 0.0, 0.0}, {0.0, 60.0, 40.0, 0.0}, {0.0, 0.0, 30.0, 90.0}},
                       {{0.0, 1.0, 4.0, 9.0}, {1.0, 0.0, 1.0, 4.0}, {9.0, 4.0, 1.0, 0.0}}, 1e-10, 10000);
  EXPECT_TRUE(chain.connected);
  EXPECT_TRUE(chain.converged);
}

// The samples of a state in a bin it forbids itself (and that another state allows) were counted in its n_k and in
// the M_b of the bin, while the state is left out of the mixture there. They are no samples of the model: the
// solution is that of the data without them.
TEST(WHAM_OVERLAP_AND_SHARED_SLOT_FIXES, dctst_wham_drops_the_samples_of_a_state_in_a_bin_it_forbids)
{
  const std::vector<std::vector<double>> bias = {{0.0, 0.0, 1e30}, {0.5, 0.0, -0.5}};
  const DcTSTWham::Result with = DcTSTWham::solve({{100.0, 80.0, 500.0}, {30.0, 60.0, 90.0}}, bias, 1e-12, 10000);
  const DcTSTWham::Result without = DcTSTWham::solve({{100.0, 80.0, 0.0}, {30.0, 60.0, 90.0}}, bias, 1e-12, 10000);
  ASSERT_TRUE(with.converged);
  ASSERT_TRUE(without.converged);
  for (std::size_t k = 0; k < 2uz; ++k) EXPECT_NEAR(with.freeEnergies[k], without.freeEnergies[k], 1e-9);
  for (std::size_t b = 0; b < 3uz; ++b) EXPECT_NEAR(with.logDensity[b], without.logDensity[b], 1e-9);

  // a state whose samples all lie in bins it forbids has none: it is no reference, and it links nothing
  const DcTSTWham::Result stranded = DcTSTWham::solve({{0.0, 0.0, 500.0}, {30.0, 60.0, 90.0}}, bias, 1e-12, 10000);
  EXPECT_TRUE(stranded.converged);
  EXPECT_NEAR(stranded.freeEnergies[1], 0.0, 1e-12);
}

// A state without samples that forbids every bin with a density: its largest term is about -1e30, far above the
// 'log of zero', and its free energy came out as +1e30 instead of the 0 of a state that takes no part.
TEST(WHAM_OVERLAP_AND_SHARED_SLOT_FIXES, dctst_wham_gives_an_unsampled_forbidden_state_a_zero_free_energy)
{
  const DcTSTWham::Result result = DcTSTWham::solve({{0.0, 0.0}, {100.0, 50.0}, {40.0, 110.0}},
                                                    {{1e30, 1e30}, {0.0, 0.0}, {0.5, -0.5}}, 1e-10, 10000);
  EXPECT_TRUE(result.converged);
  EXPECT_NEAR(result.freeEnergies[0], 0.0, 1e-6);  // (zero up to the gauge shift of the last iteration)
  EXPECT_LT(std::abs(result.freeEnergies[2]), 10.0);
}

// SwapCFCMC and SwapCBCFCMC move along the one coordinate 'lambdaGC' of the component: one fractional molecule.
// Each used to get its own while both read and wrote 'lambdaGC.currentBin'.
TEST(WHAM_OVERLAP_AND_SHARED_SLOT_FIXES, cfcmc_and_cbcfcmc_swaps_share_one_fractional_molecule)
{
  const ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  const Framework framework = Framework::makeFAU(forceField, int3(1, 1, 1));

  MCMoveProbabilities probabilities;
  probabilities.setProbability(Move::Types::Translation, 1.0);
  probabilities.setProbability(Move::Types::Rotation, 1.0);
  probabilities.setProbability(Move::Types::SwapCFCMC, 1.0);
  probabilities.setProbability(Move::Types::SwapCBCFCMC, 1.0);
  // a Widom move of either kind adds nothing next to a swap of either kind
  probabilities.setProbability(Move::Types::WidomCFCMC, 1.0);
  probabilities.setProbability(Move::Types::WidomCBCFCMC, 1.0);

  System system(forceField, std::nullopt, false, 300.0, 1e5, 1.0, {framework}, {makeCO2(forceField, probabilities)}, {},
                {10}, 5);
  EXPECT_EQ(system.numberOfFractionalMoleculesPerComponent[0], 1uz);
  EXPECT_EQ(system.numberOfIntegerMoleculesPerComponent[0], 10uz);
  EXPECT_EQ(system.numberOfMoleculesPerComponent[0], 11uz);
  EXPECT_EQ(system.indexOfFractionalMoleculeForMove(Move::Types::SwapCFCMC, 0),
            system.indexOfFractionalMoleculeForMove(Move::Types::SwapCBCFCMC, 0));
  EXPECT_EQ(system.indexOfFractionalMoleculeForMove(Move::Types::WidomCFCMC, 0),
            system.indexOfFractionalMoleculeForMove(Move::Types::WidomCBCFCMC, 0));
  EXPECT_TRUE(system.usesGCLambda(0));
  EXPECT_EQ(system.components[0].mc_moves_probabilities.getProbability(Move::Types::WidomCFCMC), 0.0);
  EXPECT_EQ(system.components[0].mc_moves_probabilities.getProbability(Move::Types::WidomCBCFCMC), 0.0);

  // the CB move alone has its molecule in that slot as well
  {
    MCMoveProbabilities onlyCB;
    onlyCB.setProbability(Move::Types::SwapCBCFCMC, 1.0);
    System other(forceField, std::nullopt, false, 300.0, 1e5, 1.0, {framework}, {makeCO2(forceField, onlyCB)}, {}, {2},
                 5);
    EXPECT_EQ(other.numberOfFractionalMoleculesPerComponent[0], 1uz);
    EXPECT_EQ(other.indexOfFractionalMoleculeForMove(Move::Types::SwapCBCFCMC, 0), 0uz);
    EXPECT_TRUE(other.usesGCLambda(0));
  }

  // both moves in one run: the running energies stay those of the configuration
  std::vector<System> systems{system};
  MonteCarlo mc = MonteCarlo({30, 0, 10, 10, 1000, 10000, 5000, 5000}, systems, 42uz, 5, false);
  mc.run();
  for (System& s : mc.systems)
  {
    EXPECT_EQ(s.numberOfFractionalMoleculesPerComponent[0], 1uz);
    const RunningEnergy drift = s.runningEnergies - s.computeTotalEnergies();
    EXPECT_NEAR(drift.potentialEnergy(), 0.0, 1e-6);
    EXPECT_NEAR(drift.frameworkMoleculeVDW, 0.0, 1e-6);
    EXPECT_NEAR(drift.moleculeMoleculeVDW, 0.0, 1e-6);
    EXPECT_NEAR(drift.frameworkMoleculeCharge, 0.0, 1e-6);
    EXPECT_NEAR(drift.moleculeMoleculeCharge, 0.0, 1e-6);
    EXPECT_NEAR(drift.ewald_fourier, 0.0, 1e-6);
    EXPECT_NEAR(drift.ewald_self, 0.0, 1e-6);
    EXPECT_NEAR(drift.ewald_exclusion, 0.0, 1e-6);
    EXPECT_NEAR(drift.tail, 0.0, 1e-6);
  }
}

// Hydrogen at 298 K and 1 bar: A < B + B^2, so the Peng-Robinson cubic has three real roots, two of them below
// the co-volume B (one negative). The fugacity coefficient of those is the exponential of the logarithm of a
// negative number; the fluid is the one root above B.
TEST(WHAM_OVERLAP_AND_SHARED_SLOT_FIXES, roots_below_the_covolume_are_no_fluid)
{
  const std::vector<EquationOfState::FluidInput> hydrogen = {{33.19, 1.313e6, -0.216, 1.0, true}};
  for (const EquationOfState::Type type :
       {EquationOfState::Type::PengRobinson, EquationOfState::Type::SoaveRedlichKwong})
  {
    for (const double pressure : {1.0e5, 1.0e7, 1.0e9})
    {
      const std::vector<EquationOfState::FluidResult> result = EquationOfState::computeFluidProperties(
          298.0, pressure, hydrogen, type, EquationOfState::MixingRules::VanDerWaals);
      ASSERT_EQ(result.size(), 1uz);
      ASSERT_TRUE(result.front().fugacityCoefficient.has_value());
      const double phi = result.front().fugacityCoefficient.value();
      const double Z = result.front().compressibility;
      EXPECT_TRUE(phi > 0.5 && phi < 1.0e6) << "phi = " << phi << " at " << pressure << " Pa";
      EXPECT_TRUE(Z > 0.9 && Z < 100.0) << "Z = " << Z << " at " << pressure << " Pa";
    }
  }
  const std::vector<EquationOfState::FluidResult> ambient = EquationOfState::computeFluidProperties(
      298.0, 1.0e5, hydrogen, EquationOfState::Type::PengRobinson, EquationOfState::MixingRules::VanDerWaals);
  EXPECT_NEAR(ambient.front().compressibility, 1.0003, 5e-4);
  EXPECT_NEAR(ambient.front().fugacityCoefficient.value(), 1.0003, 5e-4);
}

// Two replicas in the same external field: its energy travels with the configuration in the running energies.
// Any external field used to refuse the exchange.
TEST(WHAM_OVERLAP_AND_SHARED_SLOT_FIXES, replicas_in_the_same_external_field_exchange)
{
  auto make = [](double temperature, bool externalField, double boxLength)
  {
    ForceField forceField({{"X", false, 16.0, 0.0, 0.0, 1, false}}, {{150.0, 3.7}},
                          ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true, false, false);
    forceField.potentialEnergySurfaceType = ForceField::PotentialEnergySurfaceType::SecondOrderPolynomialTestFunction;
    MCMoveProbabilities probabilities;
    probabilities.setProbability(Move::Types::Translation, 1.0);
    const Component particle(forceField, "particle", 190.564, 4599200.0, 0.01142,
                             {Atom({0.0, 0.0, 0.0}, 0.0, 1.0, 0, 0, 0, false, false)}, {}, {}, 5, 21, probabilities);
    System system(forceField, SimulationBox(boxLength, boxLength, boxLength), externalField, temperature, 1e5, 1.0, {},
                  {particle}, {}, {4}, 5);
    system.runningEnergies = system.computeTotalEnergies();
    return system;
  };

  const System a = make(300.0, true, 30.0);
  const System b = make(350.0, true, 30.0);
  const std::optional<double> logR = MC_Moves::ParallelTemperingLogAcceptance(a, b);
  ASSERT_TRUE(logR.has_value());
  EXPECT_NEAR(logR.value(),
              (b.beta - a.beta) * (b.runningEnergies.potentialEnergy() - a.runningEnergies.potentialEnergy()), 1e-9);

  // a field in one replica only, or the field laid out in another box, is another Hamiltonian
  EXPECT_FALSE(MC_Moves::ParallelTemperingLogAcceptance(a, make(350.0, false, 30.0)).has_value());
  EXPECT_FALSE(MC_Moves::ParallelTemperingLogAcceptance(a, make(350.0, true, 32.0)).has_value());

  // Unless the replicas sample their volume: the box is then part of the state and travels with the
  // configuration, the two boxes differ from the first accepted volume change on, and the field, evaluated in the
  // box each configuration is in, is one Hamiltonian in both. (Refused, an isobaric ladder in a field never
  // exchanged again after that change.)
  System isobaricA = make(300.0, true, 30.0);
  System isobaricB = make(350.0, true, 32.0);
  for (System* system : {&isobaricA, &isobaricB})
  {
    system->mc_moves_probabilities.setProbability(Move::Types::VolumeChange, 0.1);
  }
  EXPECT_TRUE(MC_Moves::ParallelTemperingLogAcceptance(isobaricA, isobaricB).has_value());
  // one of them alone does not make the other's box a coordinate
  isobaricB.mc_moves_probabilities.setProbability(Move::Types::VolumeChange, 0.0);
  EXPECT_FALSE(MC_Moves::ParallelTemperingLogAcceptance(isobaricA, isobaricB).has_value());
}

// ln P and the bias compensate each other: each is beyond the range of an exponential, their sum is not.
TEST(WHAM_OVERLAP_AND_SHARED_SLOT_FIXES, biased_acceptance_is_formed_in_log_space)
{
  // exp(-800) * exp(800) is 0 * inf
  EXPECT_DOUBLE_EQ(MC_Moves::biasedAcceptance(-800.0, 800.0).sampling, 1.0);
  EXPECT_DOUBLE_EQ(MC_Moves::biasedAcceptance(800.0, -800.0).sampling, 1.0);
  EXPECT_NEAR(MC_Moves::biasedAcceptance(-800.0, 799.0).sampling, std::exp(-1.0), 1e-12);
  // the physical probability stays a number the transition matrix can clamp
  EXPECT_GE(MC_Moves::biasedAcceptance(-800.0, 800.0).physical, 0.0);
  EXPECT_LT(MC_Moves::biasedAcceptance(-800.0, 800.0).physical, 1e-300);
  EXPECT_GT(MC_Moves::biasedAcceptance(800.0, -800.0).physical, 1e300);
  EXPECT_LT(MC_Moves::biasedAcceptance(800.0, -800.0).physical, 1.7e308);
  // ordinary values are unchanged
  EXPECT_DOUBLE_EQ(MC_Moves::biasedAcceptance(-1.5, 0.25).physical, std::exp(-1.5));
  EXPECT_DOUBLE_EQ(MC_Moves::biasedAcceptance(-1.5, 0.25).sampling, std::exp(-1.25));
  // with a prefactor: 2 exp(-800) exp(800), and a prefactor of zero is a move that can not be accepted
  EXPECT_NEAR(MC_Moves::biasedAcceptance(2.0, -800.0, 800.0).sampling, 2.0, 1e-12);
  EXPECT_EQ(MC_Moves::biasedAcceptance(0.0, 5.0, 800.0).sampling, 0.0);
  EXPECT_EQ(MC_Moves::biasedAcceptance(0.0, 5.0, 800.0).physical, 0.0);
}

// The tracked component is the adsorbate whose count changes: an adsorbate of a fixed number listed before it is
// part of the environment.
TEST(WHAM_OVERLAP_AND_SHARED_SLOT_FIXES, tracked_component_is_the_swappable_adsorbate)
{
  const ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  MCMoveProbabilities fixedProbabilities;
  fixedProbabilities.setProbability(Move::Types::Translation, 1.0);
  MCMoveProbabilities openProbabilities;
  openProbabilities.setProbability(Move::Types::Translation, 1.0);
  openProbabilities.setProbability(Move::Types::Swap, 1.0);

  Component solvent = makeCO2(forceField, fixedProbabilities);
  solvent.name = "solvent";
  const Component probe = makeCO2(forceField, openProbabilities);
  {
    const System system(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e5, 1.0, {}, {solvent, probe}, {},
                        {3, 0}, 5);
    ASSERT_FALSE(system.components[0].swappable);
    ASSERT_TRUE(system.components[1].swappable);
    EXPECT_EQ(system.tmmcComponent(), 1uz);
  }
  // without a swappable adsorbate it is the first one, as before
  {
    const System system(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e5, 1.0, {}, {solvent, solvent}, {},
                        {3, 0}, 5);
    EXPECT_EQ(system.tmmcComponent(), 0uz);
  }
  // and the first of two swappable ones
  {
    const System system(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e5, 1.0, {}, {probe, probe}, {},
                        {0, 0}, 5);
    EXPECT_EQ(system.tmmcComponent(), 0uz);
  }
}

// The scout fills the empty box: the molecules the input created are removed first. Started far above the
// equilibrium loading (about 13 molecules of an ideal gas here) the loading fell block after block, and with the
// one-sided test every falling block counted towards the plateau.
TEST(WHAM_OVERLAP_AND_SHARED_SLOT_FIXES, bet_scout_starts_from_the_empty_host)
{
  const ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  // no translation: a molecule leaves only by a deletion, so that the start is felt for long
  MCMoveProbabilities probabilities;
  probabilities.setProbability(Move::Types::Swap, 1.0);
  Component probe = makeCO2(forceField, probabilities);
  probe.crossSection = 16.2;
  probe.liquidVolume = 57.7;
  probe.saturationPressure = 101325.0;
  const System system(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e5, 1.0, {}, {probe}, {}, {120}, 5);
  ASSERT_EQ(system.numberOfIntegerMoleculesPerComponent[0], 120uz);

  const NitrogenBETScoutPoint point = scoutNitrogenBETOccupancy(system, 2.0e6, 1000);
  EXPECT_GT(point.meanOccupancy, 6.0);
  EXPECT_LT(point.meanOccupancy, 40.0);
  // the peak is taken over every block: it saw the initial molecules
  EXPECT_LT(point.peakOccupancy, 60uz);
}
