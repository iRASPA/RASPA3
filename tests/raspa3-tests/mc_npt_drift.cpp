#include <gtest/gtest.h>

import std;

import double3;
import units;
import atom;
import forcefield;
import framework;
import component;
import system;
import monte_carlo;
import simulationbox;
import running_energy;
import mc_moves_move_types;
import mc_moves_probabilities;
import mc_moves_statistics;
import move_statistics;
import mc_moves_volume;
import randomnumbers;
import int3;
import ewald_structure_factor;
import connectivity_table;
import intra_molecular_potentials;
import bond_potential;

TEST(MC_NPT_DRIFT, anisotropic_volume_change)
{
  const ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);

  MCMoveProbabilities probabilities = MCMoveProbabilities();
  probabilities.setProbability(Move::Types::Translation, 1.0);

  Component co2 = Component(forceField, "CO2", 304.1282, 7377300.0, 0.22394,
                            {Atom({0, 0, 1.149}, -0.3256, 1.0, 0, 4, 0, false, false),
                             Atom({0, 0, 0.000}, 0.6512, 1.0, 0, 3, 0, false, false),
                             Atom({0, 0, -1.149}, -0.3256, 1.0, 0, 4, 0, false, false)},
                            {}, {}, 5, 21, probabilities, std::nullopt, false);

  SimulationBox box = SimulationBox(30.0, 30.0, 30.0, 100.0 * (std::numbers::pi / 180.0),
                                    95.0 * (std::numbers::pi / 180.0), 75.0 * (std::numbers::pi / 180.0));

  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::AnisotropicVolumeChange, 0.01);

  System system = System(forceField, box, false, 300.0, 1e4, 1.0, {}, {co2}, {}, {10}, 5, systemProbabilities);
  system.pressureTensorDiagonal = double3(1e4, 1.2e4, 0.8e4) / Units::PressureConversionFactor;

  std::vector<System> systems{system};
  MonteCarlo mc = MonteCarlo({20, 0, 5, 5, 1000, 10000, 5000, 5000}, systems, 42uz, 5, false);
  mc.run();

  for (System& s : mc.systems)
  {
    RunningEnergy recomputedEnergies = s.computeTotalEnergies();
    RunningEnergy drift = s.runningEnergies - recomputedEnergies;

    EXPECT_NEAR(drift.potentialEnergy(), 0.0, 1e-6);
    EXPECT_NEAR(drift.moleculeMoleculeVDW, 0.0, 1e-6);
    EXPECT_NEAR(drift.moleculeMoleculeCharge, 0.0, 1e-6);
    EXPECT_NEAR(drift.ewald_fourier, 0.0, 1e-6);
    EXPECT_NEAR(drift.ewald_self, 0.0, 1e-6);
    EXPECT_NEAR(drift.ewald_exclusion, 0.0, 1e-6);
    EXPECT_NEAR(drift.tail, 0.0, 1e-6);
  }
}

TEST(MC_NPT_DRIFT, isotropic_volume_change)
{
  const ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);

  MCMoveProbabilities probabilities = MCMoveProbabilities();
  probabilities.setProbability(Move::Types::Translation, 1.0);

  Component co2 = Component(forceField, "CO2", 304.1282, 7377300.0, 0.22394,
                            {Atom({0, 0, 1.149}, -0.3256, 1.0, 0, 4, 0, false, false),
                             Atom({0, 0, 0.000}, 0.6512, 1.0, 0, 3, 0, false, false),
                             Atom({0, 0, -1.149}, -0.3256, 1.0, 0, 4, 0, false, false)},
                            {}, {}, 5, 21, probabilities, std::nullopt, false);

  Component methane =
      Component(forceField, "methane", 190.564, 4599200, 0.01142, {Atom({0, 0, 0}, 0.0, 1.0, 0, 2, 1, false, false)},
                {}, {}, 5, 21, probabilities, std::nullopt, false);

  SimulationBox box = SimulationBox(30.0, 30.0, 30.0, 100.0 * (std::numbers::pi / 180.0),
                                    95.0 * (std::numbers::pi / 180.0), 75.0 * (std::numbers::pi / 180.0));

  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::VolumeChange, 0.01);

  System system =
      System(forceField, box, false, 300.0, 1e4, 1.0, {}, {co2, methane}, {}, {10, 15}, 5, systemProbabilities);

  std::vector<System> systems{system};
  MonteCarlo mc = MonteCarlo({20, 0, 5, 5, 1000, 10000, 5000, 5000}, systems, 42uz, 5, false);
  mc.run();

  for (System& s : mc.systems)
  {
    RunningEnergy recomputedEnergies = s.computeTotalEnergies();
    RunningEnergy drift = s.runningEnergies - recomputedEnergies;

    EXPECT_NEAR(drift.potentialEnergy(), 0.0, 1e-6);
    EXPECT_NEAR(drift.moleculeMoleculeVDW, 0.0, 1e-6);
    EXPECT_NEAR(drift.moleculeMoleculeCharge, 0.0, 1e-6);
    EXPECT_NEAR(drift.ewald_fourier, 0.0, 1e-6);
    EXPECT_NEAR(drift.ewald_self, 0.0, 1e-6);
    EXPECT_NEAR(drift.ewald_exclusion, 0.0, 1e-6);
    EXPECT_NEAR(drift.tail, 0.0, 1e-6);
  }
}

// With "CutOff": "auto" every volume change moves the VDW cut-off with the box. The shifts and tail corrections
// must follow: both on the trial box (the accepted running energy is the total at the new cut-off) and on the
// restore after a rejected trial (the subsequent translation differences use the shifts of the old cut-off).
TEST(MC_NPT_DRIFT, automatic_cutoff_follows_the_volume_moves)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, true, false);
  forceField.cutOffFrameworkVDWAutomatic = true;
  forceField.cutOffMoleculeVDWAutomatic = true;

  MCMoveProbabilities probabilities = MCMoveProbabilities();
  probabilities.setProbability(Move::Types::Translation, 1.0);

  Component methane =
      Component(forceField, "methane", 190.564, 4599200, 0.01142, {Atom({0, 0, 0}, 0.0, 1.0, 0, 6, 0, false, false)},
                {}, {}, 5, 21, probabilities, std::nullopt, false);

  // dense (liquid-like) state: the box stays small so the automatic cut-off is short and the shifts are sizeable
  SimulationBox box = SimulationBox(16.0, 16.0, 16.0);

  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::VolumeChange, 1.0);

  System system = System(forceField, box, false, 120.0, 1e8, 1.0, {}, {methane}, {}, {60}, 5, systemProbabilities);
  system.mc_moves_statistics.setMaxChange(Move::Types::VolumeChange, 0.5);

  // the shifts and tail corrections of every pair are those of the current cut-off
  auto expectCutOffDependentParametersUpToDate = [](const System& s)
  {
    ForceField reference = s.forceField;
    reference.preComputeCutOffDependentParameters();
    for (std::size_t pair = 0; pair < reference.data.size(); ++pair)
    {
      EXPECT_DOUBLE_EQ(s.forceField.data[pair].shift, reference.data[pair].shift);
      EXPECT_DOUBLE_EQ(s.forceField.data[pair].tailCorrectionEnergy, reference.data[pair].tailCorrectionEnergy);
      EXPECT_DOUBLE_EQ(s.forceField.data[pair].tailCorrectionPressure, reference.data[pair].tailCorrectionPressure);
    }
  };

  // Single volume trials: after every trial, accepted or rejected, the shifts and tail corrections must belong
  // to the cut-off the system was left with. (In a full run an accepted volume move resets the running energy to
  // the exact total, so stale shifts left by a rejected trial would only show up in the energy drift when the
  // last trial happened to be rejected; here they are checked trial by trial.)
  {
    System single = system;
    single.forceField.initializeAutomaticCutOff(single.simulationBox);
    single.runningEnergies = single.computeTotalEnergies();
    expectCutOffDependentParametersUpToDate(single);
    const double shiftBefore = single.forceField.data[6 * forceField.numberOfPseudoAtoms + 6].shift;

    RandomNumber random(7);
    std::size_t accepted = 0, rejected = 0;
    for (std::size_t trial = 0; trial < 64; ++trial)
    {
      const double cutOffBefore = single.forceField.cutOffMoleculeVDW;
      if (const std::optional<RunningEnergy> energy = MC_Moves::volumeMove(random, single))
      {
        ++accepted;
        single.runningEnergies += energy.value();
        EXPECT_NE(single.forceField.cutOffMoleculeVDW, cutOffBefore);
      }
      else
      {
        ++rejected;
        EXPECT_DOUBLE_EQ(single.forceField.cutOffMoleculeVDW, cutOffBefore);
      }
      expectCutOffDependentParametersUpToDate(single);

      const RunningEnergy drift = single.runningEnergies - single.computeTotalEnergies();
      EXPECT_NEAR(drift.potentialEnergy(), 0.0, 1e-6);
    }
    EXPECT_GT(accepted, 0uz);
    EXPECT_GT(rejected, 0uz);
    // the shift is sizeable and has actually changed at this density
    EXPECT_NE(single.forceField.data[6 * forceField.numberOfPseudoAtoms + 6].shift, shiftBefore);
    EXPECT_GT(std::abs(shiftBefore), 0.01);
  }

  // Full run with translations in between the volume moves: the cut-off follows the box and nothing drifts.
  std::vector<System> systems{system};
  MonteCarlo mc = MonteCarlo({200, 0, 20, 20, 1000, 10000, 5000, 5000}, systems, 42uz, 5, false);
  mc.run();

  for (System& s : mc.systems)
  {
    // both accepted and rejected volume trials occurred
    const MoveStatistics<double>& volumeStatistics =
        std::get<MoveStatistics<double>>(s.mc_moves_statistics[Move::Types::VolumeChange]);
    EXPECT_GT(volumeStatistics.totalAccepted, 0.0);
    EXPECT_GT(volumeStatistics.totalCounts, volumeStatistics.totalAccepted);

    const double expectedCutOff =
        0.5 * std::min({s.simulationBox.perpendicularWidths().x, s.simulationBox.perpendicularWidths().y,
                        s.simulationBox.perpendicularWidths().z}) -
        std::numeric_limits<double>::epsilon();
    EXPECT_DOUBLE_EQ(s.forceField.cutOffMoleculeVDW, expectedCutOff);
    EXPECT_NE(s.forceField.cutOffMoleculeVDW, 12.0);
    expectCutOffDependentParametersUpToDate(s);

    RunningEnergy recomputedEnergies = s.computeTotalEnergies();
    RunningEnergy drift = s.runningEnergies - recomputedEnergies;
    EXPECT_NEAR(drift.potentialEnergy(), 0.0, 1e-6);
    EXPECT_NEAR(drift.moleculeMoleculeVDW, 0.0, 1e-6);
    EXPECT_NEAR(drift.tail, 0.0, 1e-6);
  }
}

TEST(MC_NPT_DRIFT, isotropic_volume_change_flexible)
{
  const ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);

  MCMoveProbabilities probabilities = MCMoveProbabilities();
  probabilities.setProbability(Move::Types::Translation, 1.0);

  // Flexible (bonded) two-site "ethane"-like molecule: a harmonic bond makes the component flexible, so it is
  // grown/moved with intramolecular potentials. Volume moves must scale it by its center of mass, leaving the
  // internal geometry (and thus the bond energy) unchanged.
  ConnectivityTable connectivityTable(2);
  connectivityTable[0, 1] = true;
  connectivityTable[1, 0] = true;

  Potentials::IntraMolecularPotentials intraMolecularPotentials{};
  intraMolecularPotentials.bonds = {BondPotential({0, 1}, BondType::Harmonic, {96500.0, 1.54})};

  Component ethane = Component(forceField, "ethane", 305.33, 4871800.0, 0.0993,
                               {Atom({0.0, 0.0, 0.77}, 0.0, 1.0, 0, 6, 0, false, false),
                                Atom({0.0, 0.0, -0.77}, 0.0, 1.0, 0, 6, 0, false, false)},
                               connectivityTable, intraMolecularPotentials, 5, 21, probabilities, std::nullopt, false);

  SimulationBox box = SimulationBox(30.0, 30.0, 30.0, 100.0 * (std::numbers::pi / 180.0),
                                    95.0 * (std::numbers::pi / 180.0), 75.0 * (std::numbers::pi / 180.0));

  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::VolumeChange, 0.01);

  System system = System(forceField, box, false, 300.0, 1e4, 1.0, {}, {ethane}, {}, {20}, 5, systemProbabilities);

  std::vector<System> systems{system};
  MonteCarlo mc = MonteCarlo({20, 0, 5, 5, 1000, 10000, 5000, 5000}, systems, 42uz, 5, false);
  mc.run();

  for (System& s : mc.systems)
  {
    RunningEnergy recomputedEnergies = s.computeTotalEnergies();
    RunningEnergy drift = s.runningEnergies - recomputedEnergies;

    EXPECT_NEAR(drift.potentialEnergy(), 0.0, 1e-6);
    EXPECT_NEAR(drift.moleculeMoleculeVDW, 0.0, 1e-6);
    EXPECT_NEAR(drift.moleculeMoleculeCharge, 0.0, 1e-6);
    EXPECT_NEAR(drift.ewald_fourier, 0.0, 1e-6);
    EXPECT_NEAR(drift.ewald_self, 0.0, 1e-6);
    EXPECT_NEAR(drift.ewald_exclusion, 0.0, 1e-6);
    EXPECT_NEAR(drift.tail, 0.0, 1e-6);
    EXPECT_NEAR(drift.bond, 0.0, 1e-6);
  }
}

TEST(MC_NPT_DRIFT, anisotropic_volume_change_flexible)
{
  const ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);

  MCMoveProbabilities probabilities = MCMoveProbabilities();
  probabilities.setProbability(Move::Types::Translation, 1.0);

  // Flexible (bonded) two-site "ethane"-like molecule (see isotropic_volume_change_flexible). Anisotropic volume
  // moves scale each molecule by its center of mass, so bond lengths (and the bond energy) must be preserved even
  // when the three box directions scale by different factors.
  ConnectivityTable connectivityTable(2);
  connectivityTable[0, 1] = true;
  connectivityTable[1, 0] = true;

  Potentials::IntraMolecularPotentials intraMolecularPotentials{};
  intraMolecularPotentials.bonds = {BondPotential({0, 1}, BondType::Harmonic, {96500.0, 1.54})};

  Component ethane = Component(forceField, "ethane", 305.33, 4871800.0, 0.0993,
                               {Atom({0.0, 0.0, 0.77}, 0.0, 1.0, 0, 6, 0, false, false),
                                Atom({0.0, 0.0, -0.77}, 0.0, 1.0, 0, 6, 0, false, false)},
                               connectivityTable, intraMolecularPotentials, 5, 21, probabilities, std::nullopt, false);

  SimulationBox box = SimulationBox(30.0, 30.0, 30.0, 100.0 * (std::numbers::pi / 180.0),
                                    95.0 * (std::numbers::pi / 180.0), 75.0 * (std::numbers::pi / 180.0));

  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::AnisotropicVolumeChange, 0.01);

  System system = System(forceField, box, false, 300.0, 1e4, 1.0, {}, {ethane}, {}, {20}, 5, systemProbabilities);
  system.pressureTensorDiagonal = double3(1e4, 1.2e4, 0.8e4) / Units::PressureConversionFactor;

  std::vector<System> systems{system};
  MonteCarlo mc = MonteCarlo({20, 0, 5, 5, 1000, 10000, 5000, 5000}, systems, 42uz, 5, false);
  mc.run();

  for (System& s : mc.systems)
  {
    RunningEnergy recomputedEnergies = s.computeTotalEnergies();
    RunningEnergy drift = s.runningEnergies - recomputedEnergies;

    EXPECT_NEAR(drift.potentialEnergy(), 0.0, 1e-6);
    EXPECT_NEAR(drift.moleculeMoleculeVDW, 0.0, 1e-6);
    EXPECT_NEAR(drift.moleculeMoleculeCharge, 0.0, 1e-6);
    EXPECT_NEAR(drift.ewald_fourier, 0.0, 1e-6);
    EXPECT_NEAR(drift.ewald_self, 0.0, 1e-6);
    EXPECT_NEAR(drift.ewald_exclusion, 0.0, 1e-6);
    EXPECT_NEAR(drift.tail, 0.0, 1e-6);
    EXPECT_NEAR(drift.bond, 0.0, 1e-6);
  }
}

namespace
{
System makeDiluteCO2System()
{
  const ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  Component co2 = Component(forceField, "CO2", 304.1282, 7377300.0, 0.22394,
                            {Atom({0, 0, 1.149}, -0.3256, 1.0, 0, 4, 0, false, false),
                             Atom({0, 0, 0.000}, 0.6512, 1.0, 0, 3, 0, false, false),
                             Atom({0, 0, -1.149}, -0.3256, 1.0, 0, 4, 0, false, false)},
                            {}, {}, 5, 21, {}, std::nullopt, false);
  System system(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e4, 1.0, {}, {co2}, {}, {10}, 5);
  system.runningEnergies = system.computeTotalEnergies();
  return system;
}
}  // namespace

TEST(MC_NPT_DRIFT, anisotropic_pressure_work_is_a_state_function)
{
  // The pressure work of the anisotropic volume move must be the difference of a potential of the box lengths:
  // the forward and the reverse work cancel (detailed balance) and the work is independent of the path. The
  // former first-order form V_o sum_a p_aa (s_a - 1) satisfies neither, not even for a hydrostatic tensor.
  System system = makeDiluteCO2System();
  system.pressureTensorDiagonal = double3(1e4, 1.2e4, 0.8e4) / Units::PressureConversionFactor;
  ASSERT_EQ(system.referenceSimulationBox.volume, system.simulationBox.volume);

  const SimulationBox a = system.simulationBox;
  const SimulationBox b = a.scaled(double3(1.07, 0.95, 1.02));
  const SimulationBox c = b.scaled(double3(0.90, 1.10, 1.05));

  const double ab = MC_Moves::anisotropicPressureWork(system, a, b);
  const double ba = MC_Moves::anisotropicPressureWork(system, b, a);
  const double bc = MC_Moves::anisotropicPressureWork(system, b, c);
  const double ac = MC_Moves::anisotropicPressureWork(system, a, c);
  ASSERT_NE(ab, 0.0);
  const double tolerance = 1e-12 * (std::abs(ab) + std::abs(bc));
  EXPECT_NEAR(ab + ba, 0.0, tolerance);
  EXPECT_NEAR(ab + bc, ac, tolerance);

  // the deviatoric part does contribute: the work differs from the hydrostatic p dV
  const double hydrostatic =
      (system.pressureTensorDiagonal.x + system.pressureTensorDiagonal.y + system.pressureTensorDiagonal.z) / 3.0;
  EXPECT_GT(std::abs(ab - hydrostatic * (b.volume - a.volume)), 1e-6 * std::abs(ab));

  // the first-order form is not reversible
  const auto linearized = [&](const SimulationBox& o, const SimulationBox& n)
  {
    const double3& p = system.pressureTensorDiagonal;
    return o.volume * (p.x * (n.lengthA / o.lengthA - 1.0) + p.y * (n.lengthB / o.lengthB - 1.0) +
                       p.z * (n.lengthC / o.lengthC - 1.0));
  };
  EXPECT_GT(std::abs(linearized(a, b) + linearized(b, a)), 1e-4 * std::abs(ab));

  // a hydrostatic tensor recovers the NpT work p (V_n - V_o) exactly, independent of the reference cell
  system.pressureTensorDiagonal = double3(system.pressure, system.pressure, system.pressure);
  system.referenceSimulationBox = c;
  const double isotropic = MC_Moves::anisotropicPressureWork(system, a, b);
  EXPECT_NEAR(isotropic, system.pressure * (b.volume - a.volume), 1e-12 * std::abs(isotropic));
  EXPECT_GT(std::abs(linearized(a, b) - isotropic), 1e-4 * std::abs(isotropic));
}

TEST(MC_NPT_DRIFT, volume_move_rederives_the_ewald_wave_vectors_for_the_new_box)
{
  // A user-fixed Coulomb cut-off with automatic Ewald parameters: alpha is box independent but the number of
  // wave vectors follows the box widths, so an accepted volume change must re-derive it (as every driver does
  // at stage entry); a rejected one must leave the force field and the system untouched.
  std::size_t accepted = 0;
  std::size_t rejected = 0;
  for (std::size_t seed = 0; seed < 24 && (accepted < 3 || rejected < 3); ++seed)
  {
    System system = makeDiluteCO2System();
    ASSERT_FALSE(system.forceField.cutOffCoulombAutomatic);
    ASSERT_TRUE(system.forceField.automaticEwald);
    system.mc_moves_statistics.setMaxChange(Move::Types::VolumeChange, 1.0);
    // the pressure term decides between acceptance and rejection of the large volume changes proposed here
    system.pressure = (seed % 2 == 0) ? 0.0 : 1.0e3;

    const SimulationBox boxBefore = system.simulationBox;
    const int3 waveVectorsBefore = system.forceField.numberOfWaveVectors;
    const double alphaBefore = system.forceField.EwaldAlpha;
    const RunningEnergy energyBefore = system.runningEnergies;
    const EwaldStructureFactor eikBefore(system.storedEik);
    const std::span<const Atom> atoms = system.spanOfMoleculeAtoms();
    const std::vector<Atom> atomsBefore(atoms.begin(), atoms.end());

    RandomNumber random(seed);
    const std::optional<RunningEnergy> result = MC_Moves::volumeMove(random, system);
    if (result.has_value())
    {
      ++accepted;
      system.runningEnergies += result.value();
      EXPECT_NE(system.simulationBox.volume, boxBefore.volume);

      ForceField reference = system.forceField;
      reference.numberOfWaveVectors = waveVectorsBefore;
      reference.initializeEwaldParameters(system.simulationBox);
      EXPECT_EQ(system.forceField.numberOfWaveVectors.x, reference.numberOfWaveVectors.x);
      EXPECT_EQ(system.forceField.numberOfWaveVectors.y, reference.numberOfWaveVectors.y);
      EXPECT_EQ(system.forceField.numberOfWaveVectors.z, reference.numberOfWaveVectors.z);
      EXPECT_EQ(system.forceField.reciprocalIntegerCutOffSquared, reference.reciprocalIntegerCutOffSquared);
      EXPECT_DOUBLE_EQ(system.forceField.EwaldAlpha, alphaBefore);
      EXPECT_DOUBLE_EQ(system.forceField.cutOffCoulomb, 12.0);
      if (std::abs(std::log(system.simulationBox.volume / boxBefore.volume)) > 0.6)
      {
        EXPECT_NE(system.forceField.numberOfWaveVectors.x, waveVectorsBefore.x);
      }

      const RunningEnergy exact = system.computeTotalEnergies();
      const double scale = std::max(1.0, std::abs(exact.potentialEnergy()));
      EXPECT_NEAR(system.runningEnergies.potentialEnergy(), exact.potentialEnergy(), 1.0e-8 * scale);
      EXPECT_NEAR(system.runningEnergies.ewald_fourier, exact.ewald_fourier, 1.0e-8 * scale);
    }
    else
    {
      ++rejected;
      EXPECT_EQ(system.simulationBox.volume, boxBefore.volume);
      EXPECT_EQ(system.forceField.numberOfWaveVectors.x, waveVectorsBefore.x);
      EXPECT_EQ(system.forceField.numberOfWaveVectors.y, waveVectorsBefore.y);
      EXPECT_EQ(system.forceField.numberOfWaveVectors.z, waveVectorsBefore.z);
      EXPECT_EQ(system.forceField.EwaldAlpha, alphaBefore);
      EXPECT_EQ(system.runningEnergies.potentialEnergy(), energyBefore.potentialEnergy());
      EXPECT_TRUE(system.storedEik == eikBefore);
      for (std::size_t k = 0; k < atoms.size(); ++k)
      {
        EXPECT_EQ(atoms[k].position.x, atomsBefore[k].position.x);
        EXPECT_EQ(atoms[k].position.y, atomsBefore[k].position.y);
        EXPECT_EQ(atoms[k].position.z, atomsBefore[k].position.z);
      }
    }
  }
  EXPECT_GE(accepted, 1uz);
  EXPECT_GE(rejected, 1uz);
}

// A rigid framework keeps its atoms where they are while the box changes, but its stored Fourier sums belong to
// the reciprocal vectors, Ewald alpha and wave-vector count of the old box. A volume trial rebuilds them for the
// trial box, so that the accepted energy is the exact energy of the new box, and a refused trial restores them.
TEST(MC_NPT_DRIFT, volume_move_rebuilds_the_rigid_framework_fourier_sums_for_the_new_box)
{
  const ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  const Framework framework = Framework::makeFAU(forceField, int3(1, 1, 1));
  Component co2 = Component(forceField, "CO2", 304.1282, 7377300.0, 0.22394,
                            {Atom({0, 0, 1.149}, -0.3256, 1.0, 0, 4, 0, false, false),
                             Atom({0, 0, 0.000}, 0.6512, 1.0, 0, 3, 0, false, false),
                             Atom({0, 0, -1.149}, -0.3256, 1.0, 0, 4, 0, false, false)},
                            {}, {}, 5, 21, {}, std::nullopt, false);

  std::size_t accepted = 0;
  std::size_t rejected = 0;
  for (std::size_t seed = 0; seed < 24 && (accepted < 2 || rejected < 2); ++seed)
  {
    System system(forceField, std::nullopt, false, 300.0, 1e4, 1.0, framework, {co2}, {}, {8}, 5);
    system.runningEnergies = system.computeTotalEnergies();
    ASSERT_FALSE(system.fixedFrameworkStoredEik.empty());
    system.mc_moves_statistics.setMaxChange(Move::Types::VolumeChange, 0.3);
    system.pressure = (seed % 2 == 0) ? 0.0 : 1.0e3;

    const EwaldStructureFactor frameworkEikBefore(system.fixedFrameworkStoredEik);
    const RunningEnergy energyBefore = system.runningEnergies;

    RandomNumber random(seed);
    const std::optional<RunningEnergy> result = MC_Moves::volumeMove(random, system);
    if (result.has_value())
    {
      ++accepted;
      system.runningEnergies += result.value();

      // the exact energy of the new box, with the framework sums computed afresh for it
      System fresh = system;
      fresh.fixedFrameworkStoredEik.clear();
      const RunningEnergy exact = fresh.computeTotalEnergies();
      const double scale = std::max(1.0, std::abs(exact.potentialEnergy()));
      EXPECT_NEAR(system.runningEnergies.ewald_fourier, exact.ewald_fourier, 1.0e-8 * scale);
      EXPECT_NEAR(system.runningEnergies.potentialEnergy(), exact.potentialEnergy(), 1.0e-8 * scale);
      EXPECT_FALSE(system.fixedFrameworkStoredEik == frameworkEikBefore);
    }
    else
    {
      ++rejected;
      EXPECT_EQ(system.runningEnergies.potentialEnergy(), energyBefore.potentialEnergy());
      EXPECT_TRUE(system.fixedFrameworkStoredEik == frameworkEikBefore);
    }
  }
  EXPECT_GE(accepted, 1uz);
  EXPECT_GE(rejected, 1uz);
}

TEST(MC_NPT_DRIFT, volume_move_keeps_the_external_field_energy)
{
  // The external-field energy changes with the rigid translation of the molecules and must be part of the
  // acceptance and of the returned energy (a volume move that leaves it out zeroes it in the running energy).
  const auto makeSystem = []()
  {
    ForceField forceField({{"O", false, 16.0, 0.0, 0.0, 8, false}}, {{80.0, 3.0}},
                          ForceField::MixingRule::Lorentz_Berthelot, 8.0, 8.0, 8.0, false, false, false);
    forceField.potentialEnergySurfaceType = ForceField::PotentialEnergySurfaceType::MullerBrown;
    const std::uint16_t type = static_cast<std::uint16_t>(*forceField.findPseudoAtom("O"));
    Component atomic(forceField, "atomic", 64.0, 1.0e6, 0.02,
                     {Atom({0.0, 0.0, 0.0}, 0.0, 1.0, 0, type, 0, false, false)}, {}, {}, 5, 21);
    System system(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e4, 1.0, {}, {atomic},
                  {{double3(0.6, 0.2, 5.0), double3(-0.3, 1.0, 12.0), double3(1.2, 0.4, 20.0)}}, {0}, 5);
    system.hasExternalField = true;
    system.pressure = 0.0;
    system.mc_moves_statistics.setMaxChange(Move::Types::VolumeChange, 0.3);
    system.runningEnergies = system.computeTotalEnergies();
    return system;
  };

  std::size_t accepted = 0;
  for (std::size_t seed = 0; seed < 16 && accepted < 3; ++seed)
  {
    System system = makeSystem();
    const RunningEnergy exactOld = system.computeTotalEnergies();
    ASSERT_NE(exactOld.externalFieldVDW, 0.0);

    RandomNumber random(seed);
    const std::optional<RunningEnergy> result = MC_Moves::volumeMove(random, system);
    if (!result.has_value()) continue;
    ++accepted;
    system.runningEnergies += result.value();

    const RunningEnergy exactNew = system.computeTotalEnergies();
    EXPECT_NE(exactNew.externalFieldVDW, exactOld.externalFieldVDW);
    const double scale = std::max(1.0, std::abs(exactNew.potentialEnergy()));
    EXPECT_NEAR(system.runningEnergies.externalFieldVDW, exactNew.externalFieldVDW, 1.0e-10 * scale);
    EXPECT_NEAR(system.runningEnergies.potentialEnergy(), exactNew.potentialEnergy(), 1.0e-10 * scale);
  }
  EXPECT_GE(accepted, 1uz);
}

TEST(MC_NPT_DRIFT, volume_jacobian_counts_scaled_fractional_molecules)
{
  // Ideal-gas NPT: <V> = (n + 1) / (βP) with n the number of centres of mass the move scales. A reserved
  // CFCMC slot is such a centre of mass, so n includes it; counting only integer molecules would sample
  // half this volume (empty-box NPT).
  const ForceField forceField({{"A", false, 16.0, 0.0, 0.0, 8, false}}, {{0.0, 3.0}},
                              ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, false, false, false);
  MCMoveProbabilities componentProbabilities{};
  componentProbabilities.setProbability(Move::Types::SwapCFCMC, 1.0);
  Component atomic(forceField, "atomic", 64.0, 1.0e6, 0.02, {Atom({0.0, 0.0, 0.0}, 0.0, 1.0, 0, 0, 0, false, false)},
                   {}, {}, 5, 21, componentProbabilities, std::nullopt, false);
  System system(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e4, 1.0, {}, {atomic}, {}, {0}, 5);
  ASSERT_EQ(system.numberOfIntegerMoleculesPerComponent[0], 0uz);
  ASSERT_EQ(system.numberOfFractionalMoleculesPerComponent[0], 1uz);
  ASSERT_EQ(system.moleculeData.size(), 1uz);
  EXPECT_DOUBLE_EQ(MC_Moves::numberOfScaledMolecules(system), 1.0);

  system.pressure = 2.0 / (system.beta * system.simulationBox.volume);
  system.runningEnergies = system.computeTotalEnergies();
  EXPECT_NEAR(system.runningEnergies.potentialEnergy(), 0.0, 1.0e-12);
  system.mc_moves_statistics.setMaxChange(Move::Types::VolumeChange, 0.4);

  RandomNumber random(7);
  double sumV = 0.0;
  constexpr std::size_t samples = 8000;
  for (std::size_t i = 0; i != samples; ++i)
  {
    if (const std::optional<RunningEnergy> energy = MC_Moves::volumeMove(random, system))
    {
      system.runningEnergies += energy.value();
    }
    sumV += system.simulationBox.volume;
  }
  const double expected = 2.0 / (system.beta * system.pressure);
  EXPECT_NEAR(sumV / static_cast<double>(samples), expected, 0.2 * expected);
}
