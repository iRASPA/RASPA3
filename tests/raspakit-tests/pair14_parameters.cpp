#include <gtest/gtest.h>

import std;

import double3;
import double3x3;
import units;
import atom;
import atom_dynamics;
import pseudo_atom;
import vdwparameters;
import forcefield;
import component;
import system;
import simulationbox;
import connectivity_table;
import running_energy;
import intra_molecular_exclusions;
import intra_molecular_potentials;
import bond_potential;
import bend_potential;
import torsion_potential;
import van_der_waals_potential;
import potential_pair_derivatives;
import potential_pair_vdw;
import integrators_update;
import spatial_decomposition_settings;
import spatial_decomposition_force_engine;
import spatial_decomposition_device_step;

// CHARMM-style force fields give a pseudo-atom separate Lennard-Jones parameters for its 1-4 pairs. The force
// field keeps them in a second pair table ('parameters14' in the force-field file, ForceField::pair14); the 1-4
// pairs of a component then use that table (also at a 1-4 scaling of one), in every evaluation path: the energy
// and gradient routines, the CBMC pair terms, and the spatial-decomposition engine.
namespace
{
// A-B-A-B chain: regular (epsilon, sigma) in Kelvin / Angstrom, and much softer 1-4 parameters
ForceField makeForceField(bool with14)
{
  ForceField forceField({{"A", false, 14.0, 0.0, 0.0, 6, false}, {"B", false, 16.0, 0.0, 0.0, 8, false}},
                        {{80.0, 3.5}, {50.0, 3.0}}, ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0,
                        false, false, false);
  if (with14)
  {
    forceField.setPair14SelfInteraction(0, VDWParameters(20.0, 3.0));
    forceField.setPair14SelfInteraction(1, VDWParameters(10.0, 2.6));
    forceField.applyMixingRule();
    forceField.preComputeDerivedParameters();
    forceField.preComputePotentialShift();
    forceField.preComputeTailCorrection();
  }
  return forceField;
}

Component makeChain(const ForceField& forceField, double scaling14)
{
  ConnectivityTable connectivityTable(4);
  for (std::size_t i = 0; i < 3; ++i)
  {
    connectivityTable[i, i + 1] = true;
    connectivityTable[i + 1, i] = true;
  }
  Potentials::IntraMolecularPotentials potentials{};
  for (std::size_t i = 0; i < 3; ++i)
  {
    potentials.bonds.push_back(BondPotential({i, i + 1}, BondType::Harmonic, {96500.0, 1.5}));
  }
  for (std::size_t i = 0; i < 2; ++i)
  {
    potentials.bends.push_back(BendPotential({i, i + 1, i + 2}, BendType::Harmonic, {62500.0, 112.0}));
  }
  potentials.torsions.push_back(
      TorsionPotential({0, 1, 2, 3}, TorsionType::TraPPE_Extended, {2029.99, -751.83, -538.95, -22.10, -51.27}));

  std::vector<Atom> atoms{};
  for (std::size_t i = 0; i < 4; ++i)
  {
    atoms.push_back(Atom({0.0, 0.0, 0.0}, 0.0, 1.0, 0, static_cast<std::uint16_t>(i % 2), 0, false, false));
  }
  Component component = Component(forceField, "chain", 500.0, 3.0e6, 0.3, atoms, connectivityTable, potentials, 5, 21);
  component.rigid = false;
  component.intra14VanDerWaalsScaling = scaling14;
  component.buildIntraMolecularNonBondedPairs(forceField);
  return component;
}

std::vector<double3> conformation()
{
  const double3 origin(10.0, 10.0, 10.0);
  return {origin + double3(0.00, 0.00, 0.00), origin + double3(1.55, 0.00, 0.00), origin + double3(2.05, 1.40, 0.10),
          origin + double3(3.30, 1.55, 0.95)};
}

double lennardJones(double epsilonKelvin, double sigma, double r)
{
  const double s6 = std::pow(sigma / r, 6.0);
  return 4.0 * epsilonKelvin * Units::KelvinToEnergy * (s6 * s6 - s6);
}
}  // namespace

TEST(pair14_parameters, force_field_mixes_and_shifts_the_14_table)
{
  const ForceField plain = makeForceField(false);
  EXPECT_FALSE(plain.hasPair14Parameters());
  EXPECT_EQ(&plain.pair14(0, 1), &plain(0, 1));

  const ForceField forceField = makeForceField(true);
  ASSERT_TRUE(forceField.hasPair14Parameters());
  // the regular table is untouched
  EXPECT_NEAR(forceField(0, 0).parameters.x * Units::EnergyToKelvin, 80.0, 1e-10);
  EXPECT_NEAR(forceField(0, 1).parameters.x * Units::EnergyToKelvin, std::sqrt(80.0 * 50.0), 1e-10);
  EXPECT_NEAR(forceField(0, 1).parameters.y, 3.25, 1e-12);
  // the 1-4 table: the self terms as given, the cross terms by Lorentz-Berthelot
  EXPECT_NEAR(forceField.pair14(0, 0).parameters.x * Units::EnergyToKelvin, 20.0, 1e-10);
  EXPECT_NEAR(forceField.pair14(1, 1).parameters.y, 2.6, 1e-12);
  EXPECT_NEAR(forceField.pair14(0, 1).parameters.x * Units::EnergyToKelvin, std::sqrt(20.0 * 10.0), 1e-10);
  EXPECT_NEAR(forceField.pair14(0, 1).parameters.y, 2.8, 1e-12);
  EXPECT_EQ(forceField.pair14(0, 1), forceField.pair14(1, 0));
  EXPECT_EQ(forceField.pair14(0, 1).type, VDWParameters::Type::LennardJones);

  // a shifted force field shifts the 1-4 table at the same cut-off
  ForceField shifted({{"A", false, 14.0, 0.0, 0.0, 6, false}}, {{80.0, 3.5}},
                     ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true, false, false);
  shifted.setPair14SelfInteraction(0, VDWParameters(20.0, 3.0));
  shifted.applyMixingRule();
  shifted.preComputeDerivedParameters();
  shifted.preComputePotentialShift();
  EXPECT_NEAR(shifted.pair14(0, 0).shift, lennardJones(20.0, 3.0, 12.0), 1e-12);
  EXPECT_NEAR(shifted(0, 0).shift, lennardJones(80.0, 3.5, 12.0), 1e-12);
}

TEST(pair14_parameters, force_field_file_reads_parameters14)
{
  const std::filesystem::path directory = std::filesystem::temp_directory_path() / "raspa_pair14_test";
  std::filesystem::create_directories(directory);
  {
    std::ofstream file(directory / "force_field.json");
    file << R"({
  "MixingRule": "Lorentz-Berthelot",
  "TruncationMethod": "truncated",
  "TailCorrections": false,
  "CutOff": 12.0,
  "ChargeMethod": "None",
  "PseudoAtoms": [
    {"name": "A", "framework": false, "element": "C", "mass": 14.0, "charge": 0.0},
    {"name": "B", "framework": false, "element": "O", "mass": 16.0, "charge": 0.0},
    {"name": "C", "framework": false, "element": "N", "mass": 14.0, "charge": 0.0}
  ],
  "SelfInteractions": [
    {"name": "A", "type": "lennard-jones", "parameters": [80.0, 3.5], "parameters14": [20.0, 3.0]},
    {"name": "B", "type": "lennard-jones", "parameters": [50.0, 3.0], "parameters14": [10.0, 2.6]},
    {"name": "C", "type": "lennard-jones", "parameters": [60.0, 3.2]}
  ],
  "BinaryInteractions": [
    {"names": ["A", "C"], "type": "lennard-jones", "parameters14": [5.0, 2.0]}
  ]
})";
  }
  const ForceField forceField((directory / "force_field.json").string());
  std::filesystem::remove_all(directory);

  ASSERT_TRUE(forceField.hasPair14Parameters());
  EXPECT_NEAR(forceField(0, 0).parameters.x * Units::EnergyToKelvin, 80.0, 1e-10);
  EXPECT_NEAR(forceField.pair14(0, 0).parameters.x * Units::EnergyToKelvin, 20.0, 1e-10);
  EXPECT_NEAR(forceField.pair14(1, 1).parameters.x * Units::EnergyToKelvin, 10.0, 1e-10);
  // a pseudo-atom without parameters14 keeps its regular parameters in the 1-4 table
  EXPECT_NEAR(forceField.pair14(2, 2).parameters.x * Units::EnergyToKelvin, 60.0, 1e-10);
  // mixed 1-4 cross terms
  EXPECT_NEAR(forceField.pair14(0, 1).parameters.x * Units::EnergyToKelvin, std::sqrt(20.0 * 10.0), 1e-10);
  EXPECT_NEAR(forceField.pair14(1, 2).parameters.x * Units::EnergyToKelvin, std::sqrt(10.0 * 60.0), 1e-10);
  // the explicit 1-4 binary interaction, both ways; the regular A-C pair stays mixed
  EXPECT_NEAR(forceField.pair14(0, 2).parameters.x * Units::EnergyToKelvin, 5.0, 1e-10);
  EXPECT_NEAR(forceField.pair14(2, 0).parameters.y, 2.0, 1e-12);
  EXPECT_NEAR(forceField(0, 2).parameters.x * Units::EnergyToKelvin, std::sqrt(80.0 * 60.0), 1e-10);
}

TEST(pair14_parameters, component_keeps_14_pairs_apart_and_evaluates_them_with_the_14_table)
{
  const ForceField plain = makeForceField(false);
  const ForceField forceField = makeForceField(true);
  const SimulationBox box(30.0, 30.0, 30.0);

  // without a 1-4 table, a 1-4 pair at scaling one is an ordinary pair (not listed)
  Component ordinary = makeChain(plain, 1.0);
  EXPECT_TRUE(ordinary.intraMolecularPotentials.exclusions.scaledPairs.empty());

  // with a 1-4 table, the 1-4 pair is listed with its flag, also at scaling one
  Component chain = makeChain(forceField, 1.0);
  const IntraMolecularExclusions& exclusions = chain.intraMolecularPotentials.exclusions;
  ASSERT_EQ(exclusions.scaledPairs.size(), 1uz);
  EXPECT_EQ(exclusions.scaledPairs[0].atomA, 0u);
  EXPECT_EQ(exclusions.scaledPairs[0].atomB, 3u);
  EXPECT_EQ(exclusions.scaledPairs[0].scalingVDW, 1.0);
  EXPECT_TRUE(exclusions.scaledPairs[0].pair14);
  EXPECT_TRUE(exclusions.isExcludedFromPairList(0, 3));
  EXPECT_FALSE(exclusions.isExcluded(0, 3));
  EXPECT_TRUE(exclusions.scalingOf(0, 3).pair14);
  EXPECT_FALSE(exclusions.scalingOf(0, 2).pair14);

  // the energy of the 1-4 pair (the only non-excluded pair) is the Lennard-Jones of the 1-4 table
  std::vector<Atom> atoms = chain.atoms;
  const std::vector<double3> positions = conformation();
  for (std::size_t i = 0; i < 4; ++i) atoms[i].position = positions[i];
  const double r = (positions[3] - positions[0]).length();
  const RunningEnergy energy = chain.intraMolecularPotentials.computeInternalIntraVanDerWaalsEnergies(forceField, box, atoms);
  EXPECT_NEAR(energy.intraVDW, lennardJones(std::sqrt(20.0 * 10.0), 2.8, r), 1e-10 * std::abs(energy.intraVDW));
  const RunningEnergy regular = ordinary.intraMolecularPotentials.computeInternalIntraVanDerWaalsEnergies(plain, box, atoms);
  EXPECT_NEAR(regular.intraVDW, lennardJones(std::sqrt(80.0 * 50.0), 3.25, r), 1e-10 * std::abs(regular.intraVDW));
  EXPECT_NE(energy.intraVDW, regular.intraVDW);

  // the half-scaled 1-4 pair: half the 1-4 table
  Component half = makeChain(forceField, 0.5);
  ASSERT_EQ(half.intraMolecularPotentials.exclusions.scaledPairs.size(), 1uz);
  EXPECT_TRUE(half.intraMolecularPotentials.exclusions.scaledPairs[0].pair14);
  const RunningEnergy halfEnergy = half.intraMolecularPotentials.computeInternalIntraVanDerWaalsEnergies(forceField, box, atoms);
  EXPECT_NEAR(halfEnergy.intraVDW, 0.5 * energy.intraVDW, 1e-10 * std::abs(energy.intraVDW));

  // the explicit pair term (CBMC) carries the 1-4 parameters and the flag
  std::vector<VanDerWaalsPotential> terms{};
  chain.intraMolecularPotentials.forEachVanDerWaalsTerm([&](const VanDerWaalsPotential& pair) { terms.push_back(pair); });
  ASSERT_EQ(terms.size(), 1uz);
  EXPECT_TRUE(terms[0].pair14);
  EXPECT_EQ(terms[0].parameters, forceField.pair14(0, 1));
  Potentials::IntraMolecularPotentials explicitTerms{};
  explicitTerms.vanDerWaals = terms;
  EXPECT_NEAR(explicitTerms.computeInternalIntraVanDerWaalsEnergies(forceField, box, atoms).intraVDW, energy.intraVDW,
              1e-10 * std::abs(energy.intraVDW));

  // the gradient path agrees with the energy, and with finite differences
  std::vector<AtomDynamics> dynamics(4);
  const RunningEnergy gradientEnergy = chain.intraMolecularPotentials.computeInternalGradient(forceField, box, atoms, dynamics);
  EXPECT_NEAR(gradientEnergy.intraVDW, energy.intraVDW, 1e-10 * std::abs(energy.intraVDW));
  const double delta = 1e-5;
  for (std::size_t k = 0; k < 3; ++k)
  {
    double3 step{};
    (k == 0 ? step.x : k == 1 ? step.y : step.z) = delta;
    atoms[3].position = positions[3] + step;
    const double plus = chain.intraMolecularPotentials.computeInternalEnergies(forceField, box, atoms).potentialEnergy();
    atoms[3].position = positions[3] - step;
    const double minus = chain.intraMolecularPotentials.computeInternalEnergies(forceField, box, atoms).potentialEnergy();
    atoms[3].position = positions[3];
    const double numerical = (plus - minus) / (2.0 * delta);
    const double analytical = k == 0 ? dynamics[3].gradient.x : k == 1 ? dynamics[3].gradient.y : dynamics[3].gradient.z;
    EXPECT_NEAR(analytical, numerical, 1e-5 * (1.0 + std::abs(numerical)));
  }
}

TEST(pair14_parameters, spatial_decomposition_engine_uses_the_14_table)
{
  const ForceField forceField = makeForceField(true);
  Component chain = makeChain(forceField, 1.0);
  const SimulationBox box(30.0, 30.0, 30.0);

  // two chains, apart from each other
  std::vector<double3> positions = conformation();
  for (const double3& position : conformation()) positions.push_back(position + double3(8.0, 3.0, -2.0));
  System system = System(forceField, box, false, 300.0, 1e4, 1.0, {}, {chain}, {positions}, {0}, 5);
  ASSERT_EQ(system.spanOfMoleculeAtoms().size(), 8uz);

  const RunningEnergy reference = Integrators::updateGradients(
      system.moleculeData, system.spanOfMoleculeAtoms(), system.spanOfMoleculeDynamics(), system.spanOfFrameworkAtoms(),
      system.forceField, system.simulationBox, system.components, system.eik_x, system.eik_y, system.eik_z,
      system.eik_xy, system.trialEik, system.fixedFrameworkStoredEik, system.interpolationGrids,
      system.numberOfMoleculesPerComponent, system.framework, system.spanOfFrameworkDynamics(), &system.crossLinks);
  std::vector<double3> referenceGradient{};
  for (const AtomDynamics& dynamics : system.spanOfMoleculeDynamics()) referenceGradient.push_back(dynamics.gradient);
  const double r = (positions[3] - positions[0]).length();
  EXPECT_NEAR(reference.intraVDW, 2.0 * lennardJones(std::sqrt(20.0 * 10.0), 2.8, r),
              1e-10 * std::abs(reference.intraVDW));

  auto check = [&](SpatialDecompositionSettings settings, double tolerance, std::string_view label)
  {
    SpatialDecompositionForceEngine engine(settings);
    engine.initialize(system);
    const RunningEnergy energy = engine.computeGradients(system, true);
    EXPECT_NEAR(energy.intraVDW, reference.intraVDW, tolerance * std::abs(reference.intraVDW)) << label;
    EXPECT_NEAR(energy.potentialEnergy(), reference.potentialEnergy(), tolerance * std::abs(reference.potentialEnergy()))
        << label;
    std::size_t i = 0;
    for (const AtomDynamics& dynamics : system.spanOfMoleculeDynamics())
    {
      EXPECT_NEAR((dynamics.gradient - referenceGradient[i]).length(), 0.0,
                  10.0 * tolerance * (1.0 + referenceGradient[i].length()))
          << label << " atom " << i;
      ++i;
    }
    return engine.usesDeviceBonded();
  };

  for (std::size_t threads : {1uz, 2uz})
  {
    SpatialDecompositionSettings settings;
    settings.numberOfThreads = threads;
    settings.verletSkin = 1.0;
    check(settings, 1e-9, std::format("host, {} threads", threads));
  }

  // the device bonded kernels take the 1-4 parameters of the pair from the same table
  for (PairDevice pairDevice : {PairDevice::OpenCL, PairDevice::Metal, PairDevice::CUDA})
  {
    if (!DeviceStep::available(pairDevice)) continue;
    SpatialDecompositionSettings settings;
    settings.numberOfThreads = 2;
    settings.verletSkin = 1.0;
    settings.pairDevice = pairDevice;
    EXPECT_TRUE(check(settings, 1e-5, pairDeviceName(pairDevice)));
  }
}
