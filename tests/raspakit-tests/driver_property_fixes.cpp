#include <gtest/gtest.h>

#include "../test_support.hpp"

import std;

import int3;
import double3;
import units;
import atom;
import atom_dynamics;
import forcefield;
import component;
import connectivity_table;
import intra_molecular_potentials;
import bond_potential;
import simulationbox;
import system;
import thermobarostat;
import minimization_cell_layout;
import dctst_barrier;
import sample_movies;
import property_energy_histogram;
import property_density_grid;
import property_rdf;
import mc_moves_move_types;
import mc_moves_parallel_tempering_swap;
import thermostat;
import property_loading;
import property_number_of_molecules_histogram;
import simulation_schedule;
import monte_carlo;

namespace
{
ForceField makeZeroForceField()
{
  return ForceField({{"CH2", false, 14.03, 0.0, 0.0, 6, false}}, {{0.0, 3.95}},
                    ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true, true, false);
}

Component makeDimer(const ForceField& forceField, bool flexible)
{
  ConnectivityTable connectivityTable(2);
  Potentials::IntraMolecularPotentials intraMolecularPotentials{};
  if (flexible)
  {
    connectivityTable[0, 1] = true;
    connectivityTable[1, 0] = true;
    intraMolecularPotentials.bonds = {BondPotential({0, 1}, BondType::Harmonic, {5000.0, 1.5})};
  }
  return Component(
      forceField, flexible ? "flexible-dimer" : "rigid-dimer", 305.33, 4871800.0, 0.0993,
      {Atom({0.0, 0.0, 0.0}, 0.0, 1.0, 0, 0, 0, false, false), Atom({1.5, 0.0, 0.0}, 0.0, 1.0, 0, 0, 0, false, false)},
      connectivityTable, intraMolecularPotentials, 5, 21);
}

// Runs 'body' with the working directory set to a fresh temporary directory.
template <typename F>
void inTemporaryDirectory(F&& body)
{
  const TemporaryDirectory directory;
  const std::filesystem::path original = std::filesystem::current_path();
  std::filesystem::current_path(directory.path());
  body();
  std::filesystem::current_path(original);
}
}  // namespace

// The insertions of a flexible molecule sample the activity beta f / Z_intra(beta): the swap of configurations between
// temperatures does not know Z_intra, so an open ensemble with a flexible molecule is refused (a rigid one is not).
TEST(driver_property_fixes, parallel_tempering_refuses_an_open_ensemble_of_flexible_molecules)
{
  const ForceField forceField = makeZeroForceField();
  for (const bool flexible : {false, true})
  {
    Component component = makeDimer(forceField, flexible);
    ASSERT_EQ(component.rigid, !flexible);
    component.mc_moves_probabilities.setProbability(Move::Types::Swap, 1.0);
    const System system(forceField, SimulationBox(25.0, 25.0, 25.0), false, 300.0, 1e4, 1.0, {}, {component}, {}, {0},
                        5);
    EXPECT_EQ(MC_Moves::parallelTemperingUnsupportedReason(system, false).has_value(), flexible);
  }
}

// A single flexible molecule of several atoms has its centre-of-mass momentum removed as well.
TEST(driver_property_fixes, a_single_flexible_molecule_conserves_its_center_of_mass_momentum)
{
  const ForceField forceField = makeZeroForceField();
  for (const bool flexible : {false, true})
  {
    const std::vector<std::vector<double3>> positions = {{double3(5.0, 5.0, 5.0), double3(6.5, 5.0, 5.0)}};
    const System system(forceField, SimulationBox(25.0, 25.0, 25.0), false, 300.0, 1e4, 1.0, {},
                        {makeDimer(forceField, flexible)}, positions, {0}, 5);
    ASSERT_EQ(system.numberOfMolecules(), 1uz);
    EXPECT_EQ(system.conservesCenterOfMassMomentum(), flexible);
  }
}

// Under the zero centre-of-mass momentum constraint the conserved energy of the barostat carries -k_B T ln V (the
// potential of the constraint force 3 k_B T on the strain).
TEST(driver_property_fixes, thermobarostat_energy_carries_the_center_of_mass_constraint_term)
{
  Thermobarostat state(MolecularDynamicsEnsemble::NPT, CellMinimizationType::Isotropic, MonoclinicAngleType::Beta,
                       300.0, 2.0, 0.001, 30, 3, 1, 0.5);
  state.refreshDegreesOfFreedom(30, 3, 8000.0);
  EXPECT_NEAR(state.centerOfMassConstraintForce(), 3.0 * Units::KB * 300.0, 1e-12);
  const double difference = state.energy(9000.0) - state.energy(8000.0);
  EXPECT_NEAR(difference, 2.0 * 1000.0 - Units::KB * 300.0 * std::log(9000.0 / 8000.0), 1e-9);

  state.refreshDegreesOfFreedom(30, 0, 8000.0);
  EXPECT_NEAR(state.energy(9000.0) - state.energy(8000.0), 2.0 * 1000.0, 1e-9);
}

// The tagged molecule is placed on the barrier plane itself, not on the nearest of its periodic images k lambda (the
// cages of the other molecules are counted from q = 0).
TEST(driver_property_fixes, dctst_projection_is_onto_the_barrier_plane)
{
  const SimulationBox box(40.0, 20.0, 20.0);
  const DcTSTBarrier barrier(box, double3(0.5, 0.5, 0.5), double3(1.0, 0.0, 0.0), 10.0);
  for (const double q : {-14.0, -9.0, -3.0, 4.0, 8.0, 12.5})
  {
    const double3 r = barrier.point + double3(q, 3.0, -2.0);
    const double3 projected = barrier.projectOntoPlane(box, r, 0.5);
    EXPECT_NEAR(barrier.projection(box, projected), 0.5, 1e-10) << "q = " << q;
  }
}

// The force-based RDF integrates from g(0) = 0 with its absolute normalisation: a single fixed-site pair gives
// g = -beta V / (4 pi P) F . r_hat / r^2 beyond the pair distance (it was forced to 1 there by g = 1 - S(r)/S(0)).
TEST(driver_property_fixes, force_based_rdf_has_the_absolute_normalisation)
{
  const SimulationBox box(20.0, 20.0, 20.0);
  const double beta = 1.0 / (Units::KB * 300.0);
  std::vector<Atom> framework(1);
  framework[0].position = double3(5.0, 5.0, 5.0);
  std::vector<Atom> molecules(1);
  molecules[0].position = double3(7.0, 5.0, 5.0);
  // the gradient of the molecule site, so that g = 0.37 beyond r = 2
  const double gx = -0.37 * 16.0 * std::numbers::pi / (beta * box.volume);
  const std::vector<double3> moleculeGradients{double3(gx, 0.0, 0.0)};

  PropertyRadialDistributionFunction rdf(1, 1, 40, 10.0, 1, 1);
  rdf.sample(box, beta, framework, std::vector<double3>(1), std::vector<std::uint8_t>(1), molecules, moleculeGradients,
             0, 0);
  const auto [average, error] = rdf.result(0, 0);
  for (std::size_t bin = 0; bin != average.size(); ++bin)
  {
    const double r = (static_cast<double>(bin) + 0.5) * 0.25;
    EXPECT_NEAR(average[bin], r > 2.0 ? 0.37 : 0.0, 1e-10) << "r = " << r;
  }
}

// The sites of a rigid molecule move with the molecule: equal and opposite (constraint-like) site forces of a rigid
// molecule carry no RDF signal.
TEST(driver_property_fixes, force_based_rdf_takes_the_body_force_of_a_rigid_molecule)
{
  const ForceField forceField = makeZeroForceField();
  const std::vector<std::vector<double3>> positions = {
      {double3(5.0, 5.0, 5.0), double3(6.5, 5.0, 5.0), double3(5.0, 8.0, 5.0), double3(6.5, 8.0, 5.0)}};
  System system(forceField, SimulationBox(20.0, 20.0, 20.0), false, 300.0, 1e4, 1.0, {}, {makeDimer(forceField, false)},
                positions, {0}, 5);
  ASSERT_EQ(system.numberOfMolecules(), 2uz);
  system.setPropertyRDF(PropertyRadialDistributionFunction(1, 1, 40, 10.0, 1, 1));
  std::span<AtomDynamics> dynamics = system.spanOfMoleculeDynamics();
  dynamics[0].gradient = double3(0.0, 2.5, 0.0);
  dynamics[1].gradient = double3(0.0, -2.5, 0.0);
  dynamics[2].gradient = double3(0.0, 0.0, 0.0);
  dynamics[3].gradient = double3(0.0, 0.0, 0.0);
  system.sampleForceBasedRDFFromCurrentGradients(0, 0);

  const auto [average, error] = system.propertyRadialDistributionFunction->result(0, 0);
  for (std::size_t bin = 0; bin != average.size(); ++bin) EXPECT_NEAR(average[bin], 0.0, 1e-12) << bin;
}

// The number density of the density grid is normalised by the average volume of the samples, not by the volume of
// the box at the time of writing.
TEST(driver_property_fixes, density_grid_number_density_uses_the_sampled_volume)
{
  const ForceField forceField = makeZeroForceField();
  const std::vector<Component> components{makeDimer(forceField, false)};
  inTemporaryDirectory(
      [&]
      {
        PropertyDensityGrid grid(0, 1, int3(4, 4, 4), 1, 1, {}, PropertyDensityGrid::Normalization::NumberDensity,
                                 PropertyDensityGrid::Binning::Standard);
        std::vector<Atom> atoms(1);
        atoms[0].position = double3(1.0, 1.0, 1.0);
        grid.sample(std::nullopt, SimulationBox(10.0, 10.0, 10.0), atoms, 0);
        grid.writeOutput(0, SimulationBox(20.0, 20.0, 20.0), forceField, std::nullopt, components, 0);

        std::ifstream stream(std::format("density_grids/grid_cell_component_{}.s0.cube", components[0].name));
        ASSERT_TRUE(stream.good());
        std::string line;
        for (std::size_t i = 0; i != 6; ++i) std::getline(stream, line);
        double maximum{};
        double value{};
        while (stream >> value) maximum = std::max(maximum, value);
        // one atom in a voxel of 1000 / 64 A^3
        EXPECT_NEAR(maximum, 64.0 / 1000.0, 1e-12);
      });
}

// The movie file of a run resumed from a binary restart is not truncated when the writer is constructed.
TEST(driver_property_fixes, sample_movie_construction_keeps_an_existing_movie)
{
  inTemporaryDirectory(
      [&]
      {
        std::filesystem::create_directory("movies");
        {
          std::ofstream stream("movies/movie.s0.pdb");
          stream << "MODEL 1\n";
        }
        const SampleMovie movie(0, 1, false, std::nullopt);
        std::ifstream stream("movies/movie.s0.pdb");
        const std::string contents((std::istreambuf_iterator<char>(stream)), std::istreambuf_iterator<char>());
        EXPECT_EQ(contents, "MODEL 1\n");
      });
}

// The energy-histogram range is printed in the displayed unit (Kelvin), as it was given.
TEST(driver_property_fixes, energy_histogram_prints_its_range_in_kelvin)
{
  const PropertyEnergyHistogram histogram(1, 10, {-1000.0, 500.0}, 1, std::nullopt);
  const std::string settings = histogram.printSettings();
  const std::size_t position = settings.find("valueRange: (");
  ASSERT_NE(position, std::string::npos);
  EXPECT_NEAR(std::stod(settings.substr(position + 13)), -1000.0, 1e-9);
}

// Out-of-range indices are refused (they aliased the next row of channels or read past the averages).
TEST(driver_property_fixes, property_accessors_refuse_out_of_range_indices)
{
  const PropertyRadialDistributionFunction rdf(1, 2, 10, 5.0, 1, 1);
  EXPECT_NO_THROW(static_cast<void>(rdf.result(1, 1)));
  EXPECT_THROW(static_cast<void>(rdf.result(0, 2)), std::out_of_range);
  EXPECT_THROW(static_cast<void>(rdf.result(2, 0)), std::out_of_range);

  const PropertyLoading loading(5, 1);
  EXPECT_NO_THROW(static_cast<void>(loading.averageLoadingNumberOfMolecules(0)));
  EXPECT_THROW(static_cast<void>(loading.averageLoadingNumberOfMolecules(1)), std::out_of_range);
}

// An empty or reversed molecule-count range is refused (second - first wrapped around).
TEST(driver_property_fixes, number_of_molecules_histogram_refuses_an_empty_range)
{
  EXPECT_THROW(PropertyNumberOfMoleculesHistogram(5, {10uz, 0uz}, 1, std::nullopt), std::invalid_argument);
  EXPECT_THROW(PropertyNumberOfMoleculesHistogram(5, 1, {3uz, 3uz}, 1, std::nullopt), std::invalid_argument);
  EXPECT_EQ(PropertyNumberOfMoleculesHistogram(5, {2uz, 12uz}, 1, std::nullopt).numberOfBins, 10uz);
}

// A driver with more blocks than its systems wrote past the block averages of the systems.
TEST(driver_property_fixes, monte_carlo_refuses_a_block_count_different_from_its_systems)
{
  const ForceField forceField = makeZeroForceField();
  const System system(forceField, SimulationBox(25.0, 25.0, 25.0), false, 300.0, 1e4, 1.0, {},
                      {makeDimer(forceField, false)}, {}, {0}, 5);
  EXPECT_NO_THROW(MonteCarlo(SimulationSchedule{10}, {system}, 12, 5, false));
  EXPECT_THROW(MonteCarlo(SimulationSchedule{10}, {system}, 12, 8, false), std::invalid_argument);
}

// None clears the thermostat and the movie, as it clears every other optional property; the property setters size
// what they are given for the system.
TEST(driver_property_fixes, system_setters_clear_on_none_and_size_for_the_system)
{
  const ForceField forceField = makeZeroForceField();
  System system(forceField, SimulationBox(25.0, 25.0, 25.0), false, 300.0, 1e4, 1.0, {}, {makeDimer(forceField, false)},
                {}, {0}, 5);

  system.setThermostat(Thermostat(5, 5, 0.15));
  ASSERT_TRUE(system.thermostat.has_value());
  system.setThermostat(std::nullopt);
  EXPECT_FALSE(system.thermostat.has_value());

  system.setSamplePDBMovie(SampleMovie(0, 10, true, std::nullopt));
  ASSERT_TRUE(system.samplePDBMovie.has_value());
  system.setSamplePDBMovie(std::nullopt);
  EXPECT_FALSE(system.samplePDBMovie.has_value());

  system.setPropertyDensityGrid(PropertyDensityGrid(
      0, 3, int3(4, 4, 4), 1, 1, {}, PropertyDensityGrid::Normalization::Max, PropertyDensityGrid::Binning::Standard));
  ASSERT_TRUE(system.propertyDensityGrid.has_value());
  EXPECT_EQ(system.propertyDensityGrid->numberOfComponents, 1uz);
  EXPECT_THROW(system.setPropertyDensityGrid(PropertyDensityGrid(0, 1, int3(4, 4, 4), 1, 1, {7},
                                                                 PropertyDensityGrid::Normalization::Max,
                                                                 PropertyDensityGrid::Binning::Standard)),
               std::out_of_range);

  system.setAverageEnergyHistogram(PropertyEnergyHistogram(9, 16, {-1000.0, 0.0}, 1, std::nullopt));
  ASSERT_TRUE(system.averageEnergyHistogram.has_value());
  EXPECT_EQ(system.averageEnergyHistogram->histogram.numberOfBlocks, 5uz);
}
