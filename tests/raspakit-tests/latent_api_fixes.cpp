#include <gtest/gtest.h>

import std;

import int3;
import double3;
import units;
import atom;
import atom_dynamics;
import pseudo_atom;
import vdwparameters;
import forcefield;
import framework;
import component;
import connectivity_table;
import intra_molecular_potentials;
import simulationbox;
import system;
import running_energy;
import interpolation_energy_grid;
import interactions_intermolecular;
import interactions_polarization;
import interactions_framework_molecule;
import integrators_cputime;
import json;
import bond_potential;
import distance_potential_gradient_hessian_strain;
import bend_potential;
import urey_bradley_potential;
import cbmc_chain_common;
import cbmc_grow_step;
import randomnumbers;
import property_lambda_probability_histogram;
import property_volume_evolution;
import averages;
import mc_moves_cputime;
import mc_moves_move_types;
import mc_moves_probabilities;
import molecule;
import potential_intra_pair;

namespace
{
ForceField makeChargedForceField()
{
  return ForceField({PseudoAtom("A", false, 12.0, 0.4, 1.5, 6, false)}, {VDWParameters(120.0, 3.4)},
                    ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true, false, true);
}
}  // namespace

// With the intermolecular polarization switched off, the field-difference routines still return the van der Waals
// and Coulomb energy difference (they returned zero, while their contract is the energy of the move) and leave the
// field buffers alone.
TEST(latent_api_fixes, field_difference_energies_do_not_depend_on_the_polarization_switch)
{
  ForceField forceField = makeChargedForceField();
  forceField.omitInterPolarization = true;
  const SimulationBox box(30.0, 30.0, 30.0);
  const std::vector<Atom> atoms{Atom({5.0, 5.0, 5.0}, 0.4, 1.0, 0, 0, 0, 0, false),
                                Atom({9.0, 5.0, 5.0}, -0.4, 1.0, 1, 0, 0, 0, false)};
  const std::vector<Atom> newAtoms{Atom({8.4, 5.6, 5.2}, -0.4, 1.0, 1, 0, 0, 0, false)};
  const std::span<const Atom> oldAtoms(&atoms[1], 1);

  const std::optional<RunningEnergy> reference =
      Interactions::computeInterMolecularEnergyDifference(forceField, box, atoms, newAtoms, oldAtoms);
  ASSERT_TRUE(reference.has_value());
  ASSERT_NE(reference->moleculeMoleculeVDW, 0.0);
  ASSERT_NE(reference->moleculeMoleculeCharge, 0.0);

  std::vector<double3> fieldMolecules(atoms.size());
  std::vector<double3> fieldMolecule(newAtoms.size());
  const std::optional<RunningEnergy> withField = Interactions::computeInterMolecularElectricFieldDifference(
      forceField, box, fieldMolecules, fieldMolecule, atoms, newAtoms, oldAtoms);
  ASSERT_TRUE(withField.has_value());
  EXPECT_NEAR(withField->moleculeMoleculeVDW, reference->moleculeMoleculeVDW, 1e-10);
  EXPECT_NEAR(withField->moleculeMoleculeCharge, reference->moleculeMoleculeCharge, 1e-10);

  std::vector<double3> neighbourDelta(atoms.size());
  std::vector<double3> fieldNew(newAtoms.size());
  std::vector<double3> fieldOld(newAtoms.size());
  const std::optional<RunningEnergy> polarizationVariant =
      Interactions::computeInterMolecularPolarizationElectricFieldDifference(forceField, box, neighbourDelta, fieldNew,
                                                                             fieldOld, atoms, newAtoms, oldAtoms);
  ASSERT_TRUE(polarizationVariant.has_value());
  EXPECT_NEAR(polarizationVariant->moleculeMoleculeVDW, reference->moleculeMoleculeVDW, 1e-10);
  EXPECT_NEAR(polarizationVariant->moleculeMoleculeCharge, reference->moleculeMoleculeCharge, 1e-10);

  for (const std::vector<double3>* buffer : {&fieldMolecules, &fieldMolecule, &neighbourDelta, &fieldNew, &fieldOld})
  {
    for (const double3& field : *buffer) EXPECT_EQ(double3::dot(field, field), 0.0);
  }
}

// The three-argument polarization difference takes the new field first; an insertion passes no old field, which
// gives no old term (the old field was read past its end).
TEST(latent_api_fixes, polarization_difference_of_an_insertion_has_no_old_term)
{
  const ForceField forceField = makeChargedForceField();
  std::vector<Atom> atoms{Atom({5.0, 5.0, 5.0}, 0.4, 1.0, 0, 0, 0, 0, false)};
  std::vector<double3> newField{double3(0.3, -0.2, 0.1)};
  std::vector<double3> oldField{double3(0.1, 0.1, 0.1)};

  const double polarizability = forceField.pseudoAtoms[0].polarizability / Units::CoulombicConversionFactor;
  const double newTerm = -0.5 * polarizability * double3::dot(newField[0], newField[0]);
  const double oldTerm = 0.5 * polarizability * double3::dot(oldField[0], oldField[0]);

  EXPECT_NEAR(Interactions::computePolarizationEnergyDifference(forceField, newField, {}, atoms).polarization, newTerm,
              1e-14);
  EXPECT_NEAR(Interactions::computePolarizationEnergyDifference(forceField, newField, oldField, atoms).polarization,
              newTerm + oldTerm, 1e-14);
}

// A caller without grids (an empty grid vector) gets the pair sum: the per-type grid was indexed without a bounds
// check, and the Coulomb grid read without checking it exists.
TEST(latent_api_fixes, framework_molecule_gradient_without_grids_evaluates_the_pairs)
{
  ForceField forceField = makeChargedForceField();
  const SimulationBox box(24.0, 24.0, 24.0);
  const Atom frameworkAtom({0.25, 0.25, 0.25}, 0.4, 1.0, 0, 0, 0, 0, true);
  Framework framework(forceField, "framework", box, 1, {frameworkAtom}, {frameworkAtom}, int3(1, 1, 1));
  Component component(forceField, "site", 100.0, 1.0e6, 0.1, {Atom({0.0, 0.0, 0.0}, -0.4, 1.0, 0, 0, 0, false, false)},
                      ConnectivityTable(1), Potentials::IntraMolecularPotentials{}, 0, 0);
  System system(forceField, box, false, 300.0, 1.0e4, 1.0, framework, {component}, {}, {1}, 5);
  system.spanOfFrameworkAtoms()[0].position = {7.0, 8.0, 9.0};
  system.spanOfMoleculeAtoms()[0].position = {10.6, 8.4, 8.7};

  const RunningEnergy reference = Interactions::computeFrameworkMoleculeEnergy(
      system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
      system.spanOfFrameworkAtoms(), system.spanOfMoleculeAtoms());
  ASSERT_NE(reference.frameworkMoleculeVDW, 0.0);

  const std::vector<std::optional<InterpolationEnergyGrid>> noGrids{};
  const RunningEnergy gradient = Interactions::computeFrameworkMoleculeGradient(
      system.forceField, system.simulationBox, system.spanOfFrameworkAtoms(), system.spanOfMoleculeAtoms(),
      system.spanOfMoleculeDynamics(), noGrids, system.framework);
  EXPECT_NEAR(gradient.frameworkMoleculeVDW, reference.frameworkMoleculeVDW, 1e-10);
  EXPECT_NEAR(gradient.frameworkMoleculeCharge, reference.frameworkMoleculeCharge, 1e-10);
}

// The molecular-dynamics difference table is defined (it was declared only), and the MD outputs carry the
// thermobarostat and conserved energies; the timing total leaves out the rotations nested in the free rotor.
TEST(latent_api_fixes, md_energy_outputs_carry_every_term_and_the_timing_total_counts_once)
{
  RunningEnergy energy{};
  energy.moleculeMoleculeVDW = 1.0;
  energy.translationalKineticEnergy = 2.0;
  energy.rotationalKineticEnergy = 0.5;
  energy.NoseHooverEnergy = 0.25;
  energy.thermobarostatEnergy = 0.125;
  RunningEnergy recomputed = energy;
  recomputed.thermobarostatEnergy = 0.0;

  const std::string difference = energy.printMDDiff(recomputed);
  EXPECT_NE(difference.find("Thermobarostat enthalpy"), std::string::npos);
  EXPECT_NE(difference.find("Conserved energy"), std::string::npos);
  EXPECT_NE(energy.printMD().find("Thermobarostat enthalpy"), std::string::npos);
  EXPECT_NE(energy.repr().find("Thermobarostat enthalpy"), std::string::npos);

  const auto json = energy.jsonMD();
  EXPECT_TRUE(json.contains("Conserved energy [K]"));
  EXPECT_TRUE(json.contains("Baro/thermostat energy [K]"));
  EXPECT_TRUE(json.contains("framework-molecule Real [K]"));
  EXPECT_NEAR(json["Conserved energy [K]"].get<double>(), Units::EnergyToKelvin * energy.conservedEnergy(), 1e-9);

  IntegratorsCPUTime time{};
  time.noSquishFreeRotorOrderTwo = std::chrono::duration<double>(3.0);
  time.noSquishRotate = std::chrono::duration<double>(2.0);
  EXPECT_DOUBLE_EQ(time.total().count(), 3.0);
}

// A two-parameter bond and bend left the unused parameters as whatever the stack held: two equal potentials
// compared unequal, and the normalization caches keyed on the parameters missed (or, with NaN bytes, broke the
// ordering of their maps).
TEST(latent_api_fixes, unused_bond_and_bend_parameters_are_zero)
{
  for (int repeat = 0; repeat < 4; ++repeat)
  {
    const BondPotential bond({0, 1}, BondType::Harmonic, {1000.0, 1.5});
    EXPECT_EQ(bond.parameters[2], 0.0);
    EXPECT_EQ(bond.parameters[3], 0.0);
    EXPECT_TRUE(bond == BondPotential({0, 1}, BondType::Harmonic, {1000.0, 1.5}));
    const BendPotential bend({0, 1, 2}, BendType::Harmonic, {62500.0, 114.0});
    EXPECT_EQ(bend.parameters[2], 0.0);
    EXPECT_EQ(bend.parameters[3], 0.0);
    EXPECT_TRUE(bend == BendPotential({0, 1, 2}, BendType::Harmonic, {62500.0, 114.0}));
  }
}

// A lambda-bias file is resampled onto the bins at the nearest of its points: truncating i delta (M - 1), which
// lands just below an integer for some bins in IEEE arithmetic, read the value of the preceding point into those
// bins even for a file of the same size. (The fast-math build may fold the product to i exactly, so there this is a
// guard of the identity rather than a reproduction.)
TEST(latent_api_fixes, lambda_bias_file_of_the_same_size_is_read_bin_for_bin)
{
  const std::filesystem::path path = std::filesystem::temp_directory_path() / "raspa3_lambda_bias.json";
  {
    std::ofstream file(path);
    file << "{\"bias\": [";
    for (int i = 0; i < 19; ++i) file << (i == 0 ? "" : ", ") << i;
    file << "]}";
  }
  // 19 bins: i delta (M - 1) lands below i for i = 7, 11, 14 and 15
  PropertyLambdaProbabilityHistogram histogram(5, 19);
  histogram.readBiasingFile(path);
  for (std::size_t i = 0; i < 19; ++i) EXPECT_EQ(histogram.biasFactor[i], static_cast<double>(i)) << i;
  std::filesystem::remove(path);
}

// A periodic write of an evolution comes after the sample of its own cycle: writing at cycle 4 with a sample every
// cycle has the samples 0 to 4, five of them, and the write printed only the first four.
TEST(latent_api_fixes, evolution_written_during_the_run_includes_the_sample_of_its_cycle)
{
  const std::filesystem::path original = std::filesystem::current_path();
  const std::filesystem::path directory = std::filesystem::temp_directory_path() / "raspa3_volume_evolution";
  std::filesystem::remove_all(directory);
  std::filesystem::create_directories(directory);
  std::filesystem::current_path(directory);
  {
    PropertyVolumeEvolution evolution(10, 1, 4);
    for (std::size_t cycle = 0; cycle <= 4; ++cycle) evolution.addSample(cycle, 1000.0 + static_cast<double>(cycle));
    evolution.writeOutput(0, 4);
  }
  std::ifstream file("volume_evolution/volume_evolution.s0.txt");
  std::size_t lines = 0;
  for (std::string line; std::getline(file, line);)
    if (!line.empty() && line[0] != '#') ++lines;
  file.close();
  std::filesystem::current_path(original);
  std::filesystem::remove_all(directory);
  EXPECT_EQ(lines, 5uz);
}

// The confidence interval of the defined blocks takes the Student-t value of their own degrees of freedom: three
// blocks (two degrees of freedom) took the value of the configured five blocks (t = 2.776 instead of 4.303).
TEST(latent_api_fixes, mean_confidence_uses_the_student_t_of_its_blocks)
{
  std::vector<double> blocks{1.0, 2.0, 3.0};
  const auto [mean, confidence] = meanConfidence(blocks);
  EXPECT_DOUBLE_EQ(mean, 2.0);
  const double standardError = std::sqrt(2.0 / 6.0);
  EXPECT_NEAR(confidence, studentTCriticalValue(2) * standardError, 1e-12);
  EXPECT_GT(confidence, 4.0 * standardError);
}

// The reaction moves book their times on the system, but the system report listed the system and the cross-system
// moves only (and a component never holds them): their times were reported nowhere.
TEST(latent_api_fixes, reaction_move_times_are_reported_with_the_system)
{
  MCMoveCpuTime timings;
  timings[Move::Types::ReactionCFCMC][Move::Timing::Total] = std::chrono::duration<double>(1.25);
  EXPECT_NE(
      timings.writeMCMoveCPUTimeStatistics().find(Move::moveNames[std::to_underlying(Move::Types::ReactionCFCMC)]),
      std::string::npos);
  EXPECT_TRUE(timings.jsonSystemMCMoveCPUTimeStatistics().contains(
      Move::moveNames[std::to_underlying(Move::Types::ReactionCFCMC)]));
}

// The lambda coordinate of the Gibbs CFCMC move is sampled (System::samplePropertiesMC) and now also reported, in
// the text and the JSON move statistics; only the grand-canonical lambda was.
TEST(latent_api_fixes, gibbs_lambda_statistics_are_reported)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  Component methane = Component::makeMethane(forceField, 0);
  methane.mc_moves_probabilities.setProbability(Move::Types::GibbsConventionalCFCMC, 1.0);
  System system =
      System(forceField, SimulationBox(25.0, 25.0, 25.0), false, 300.0, 1e4, 1.0, {}, {methane}, {}, {2}, 5);

  EXPECT_NE(system.writeMCMoveStatistics().find("GibbsCFCMC lambda statistics"), std::string::npos);
  const nlohmann::json json = system.jsonMCMoveStatistics();
  EXPECT_TRUE(json[system.components[0].name]["lambdaStatistics"].contains("GibbsCFCMC"));
}

// A massless probe molecule has no centre of mass: its record takes the geometric centre, as the Component
// constructor does, instead of 0/0 (NaN).
TEST(latent_api_fixes, massless_probe_molecule_record_takes_the_geometric_centre)
{
  ForceField forceField({{"X", false, 0.0, 0.0, 0.0, 0, false}}, {{10.0, 3.0}},
                        ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, false, false, false);
  const Component probe(forceField, "probe", 1.0, 1.0, 0.0,
                        {Atom(double3(0.0, 0.0, 0.0), 0.0, 1.0, 0.0, 0u, 0u, 0u, 0u, false)}, {}, {}, 5, 21);
  const std::vector<Atom> atoms{Atom(double3(1.0, 2.0, 3.0), 0.0, 1.0, 0.0, 0u, 0u, 0u, 0u, false)};
  const double3 com = probe.createMoleculeRecord(atoms).centerOfMassPosition;
  EXPECT_EQ(std::to_string(com.x), std::to_string(1.0));
  EXPECT_EQ(std::to_string(com.y), std::to_string(2.0));
  EXPECT_EQ(std::to_string(com.z), std::to_string(3.0));
}

// The scaled 1-4 Coulomb pair of a flexible framework, with the framework's own exclusion correction (Ewald:
// -C qq erf(alpha r)/r; real-space methods: C qq (V - 1/r) inside the cutoff) added on top, is the molecule's scaled
// pair: s C qq / r inside the cutoff, beyond it the long-range erf part under Ewald and nothing under Wolf, the
// scaled cut-off potential under a plain cut-off. The framework used the bare s C qq / r at every r and method.
TEST(latent_api_fixes, framework_scaled_14_coulomb_is_the_molecule_pair)
{
  ForceField forceField = ForceField::makeZeoliteForceField(6.0, true, false, true);
  forceField.EwaldAlpha = 0.3;
  const double s = 0.5;
  const double qA = 0.4;
  const double qB = -0.6;
  const double C = Units::CoulombicConversionFactor * qA * qB;
  for (const double r : {3.0, 8.0})
  {
    const bool inside = r < forceField.cutOffCoulomb;
    forceField.chargeMethod = ForceField::ChargeMethod::Ewald;
    EXPECT_NEAR(Potentials::frameworkScaledPairCoulomb<0>(forceField, s, r, qA, qB).energy,
                inside ? s * C / r : C * std::erf(0.3 * r) / r, 1e-10 * std::abs(C))
        << "Ewald at r = " << r;
    forceField.chargeMethod = ForceField::ChargeMethod::Wolf;
    EXPECT_NEAR(Potentials::frameworkScaledPairCoulomb<0>(forceField, s, r, qA, qB).energy, inside ? s * C / r : 0.0,
                1e-10 * std::abs(C))
        << "Wolf at r = " << r;

    // the gradient factor is the derivative of the energy, (dU/dr)/r
    for (const ForceField::ChargeMethod method : {ForceField::ChargeMethod::Ewald, ForceField::ChargeMethod::Wolf})
    {
      forceField.chargeMethod = method;
      const double h = 1e-6;
      const double derivative = (Potentials::frameworkScaledPairCoulomb<0>(forceField, s, r + h, qA, qB).energy -
                                 Potentials::frameworkScaledPairCoulomb<0>(forceField, s, r - h, qA, qB).energy) /
                                (2.0 * h);
      EXPECT_NEAR(r * Potentials::frameworkScaledPairCoulomb<1>(forceField, s, r, qA, qB).firstDerivativeFactor,
                  derivative, 1e-6 * std::abs(C));
    }
  }
}

// A Rigid bend takes its angle in degrees, as a Fixed one: the growth places the bead at that angle (in radians) and
// the cone normalization is log(2 pi sin theta0). The degrees were kept, the cone was at 109.47 rad and a bend at 180
// degrees had the logarithm of a negative sine.
TEST(latent_api_fixes, rigid_bend_angle_is_converted_to_radians)
{
  RandomNumber random(7);
  for (const double degrees : {109.47, 180.0 - 1e-6})
  {
    const BendPotential bend({0, 1, 2}, BendType::Rigid, {degrees});
    EXPECT_NEAR(bend.generateBendAngle(random, 1.0), degrees * Units::DegreesToRadians, 1e-12);
    const double normalization = bend.logBoltzmannConeNormalization(1.0);
    EXPECT_EQ(std::to_string(normalization),
              std::to_string(std::log(2.0 * std::numbers::pi * std::sin(degrees * Units::DegreesToRadians))));
  }
}

// A core-shell spring is (1/2) p_0 r^2: its Hessian is p_0 I at every separation, also with the shell on the core
// (r = 0, its place at zero field) and at a separation whose cube underflows. It was nought at r = 0 and not a number
// just above it.
TEST(latent_api_fixes, core_shell_spring_hessian_is_p0_at_zero_separation)
{
  std::array<double, maximumNumberOfBondParameters> parameters{};
  parameters[0] = 1500.0;
  for (const double separation : {0.0, 1e-160, 0.3})
  {
    const auto [energy, gradient, strain, f1, f2] = Potentials::Internal::distancePotentialEnergyGradientHessianStrain(
        BondType::CoreShellSpring, parameters, double3(separation, 0.0, 0.0), double3(0.0, 0.0, 0.0), false);
    EXPECT_EQ(std::to_string(f1), std::to_string(1500.0)) << separation;
    EXPECT_EQ(std::to_string(f2), std::to_string(0.0)) << separation;
  }
}

// On an overlap the polarization field difference still completes its field deltas (the removed old position
// included) before it reports the overlap: the moves that apply the deltas whatever the energy is were left with the
// new position's field added and the old one's never taken out.
TEST(latent_api_fixes, polarization_field_difference_is_complete_on_an_overlap)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  forceField.chargeMethod = ForceField::ChargeMethod::Coulomb;
  forceField.computePolarization = true;
  forceField.omitInterPolarization = false;
  const SimulationBox box(30.0, 30.0, 30.0);
  const std::vector<Atom> neighbours{Atom(double3(10.0, 10.0, 10.0), 0.5, 1.0, 1.0, 0, 6, 0, 0, false)};
  const std::vector<Atom> newAtoms{Atom(double3(10.2, 10.0, 10.0), -0.5, 1.0, 1.0, 1, 6, 0, 0, false)};
  const std::vector<Atom> oldAtoms{Atom(double3(15.0, 10.0, 10.0), -0.5, 1.0, 1.0, 1, 6, 0, 0, false)};

  auto deltas = [&](bool checkOverlap)
  {
    std::vector<double3> neighbour(1), moleculeNew(1), moleculeOld(1);
    const std::optional<RunningEnergy> energy = Interactions::computeInterMolecularPolarizationElectricFieldDifference(
        forceField, box, neighbour, moleculeNew, moleculeOld, neighbours, newAtoms, oldAtoms, checkOverlap);
    return std::pair{energy.has_value(), neighbour[0]};
  };
  const auto [checkedHasValue, checked] = deltas(true);
  const auto [uncheckedHasValue, unchecked] = deltas(false);
  EXPECT_FALSE(checkedHasValue);
  EXPECT_TRUE(uncheckedHasValue);
  EXPECT_NEAR(checked.x, unchecked.x, 1e-12 * std::abs(unchecked.x));
}

// The internal terms among the beads of the seed fragment (a rigid molecule's Urey-Bradley term here) are a constant
// of its geometry and leave the CBMC weight alone; an attach step still weighs the terms it does not sample.
TEST(latent_api_fixes, cbmc_seed_step_leaves_its_constant_internal_terms_out_of_the_weight)
{
  std::vector<Atom> atoms{Atom(double3(0.0, 0.0, 0.0), 0.0, 1.0, 1.0, 0, 0, 0, 0, false),
                          Atom(double3(2.0, 0.0, 0.0), 0.0, 1.0, 1.0, 0, 0, 0, 0, false)};
  CBMC::GrowStep step;
  step.kind = CBMC::GrowStep::Kind::PlaceSeedFragment;
  step.nextBeads = {0, 1};
  step.intra.ureyBradleys.push_back(UreyBradleyPotential({0, 1}, UreyBradleyType::Harmonic, {1000.0, 1.5}));
  EXPECT_EQ(CBMC::unsampledInternalLogFactor(1.0, step, atoms), 0.0);
  step.kind = CBMC::GrowStep::Kind::AttachFragment;
  step.rigidBody = true;
  EXPECT_LT(CBMC::unsampledInternalLogFactor(1.0, step, atoms), 0.0);
}
