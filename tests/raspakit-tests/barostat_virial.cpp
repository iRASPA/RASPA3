#include <gtest/gtest.h>
#include <unistd.h>

import std;

import int3;
import double3;
import double3x3;
import framework;
import units;
import atom;
import atom_dynamics;
import pseudo_atom;
import vdwparameters;
import forcefield;
import component;
import system;
import simulationbox;
import energy_status;
import running_energy;
import potential_coulomb_real_space;
import interactions_intermolecular;
import bond_potential;
import bend_potential;
import torsion_potential;
import connectivity_table;
import intra_molecular_potentials;
import molecule;
import randomnumbers;
import integrators_update;
import elastic_constants;
import thermobarostat;
import interactions_ewald;
import minimization_cell_layout;
import interactions_internal;
import molecular_dynamics;
import ewald_structure_factor;
import bond_potential;

// The molecular-dynamics barostat couples to the centers of mass of rigid molecules and to every atom
// of a flexible molecule individually. Its driving virial must therefore be the derivative of the total
// energy under the affine scaling of exactly those points, which for flexible molecules is the atomic
// virial including the bonded forces, not the molecular (center-of-mass) virial that the Monte Carlo
// volume move and the reported pressure use. Feeding the molecular virial together with the atomic
// kinetic energy over-estimates the pressure by (N_atoms - N_molecules) k T / V (about 1700 bar for a
// liquid of 2000 HDDA) and blows up an NPT run; 'computeBarostatVirial' converts one into the other.

namespace
{
void useSecondOrderTaylorShiftedLennardJones(ForceField &forceField)
{
  for (VDWParameters &parameters : forceField.data)
  {
    if (parameters.type == VDWParameters::Type::LennardJones)
    {
      parameters.type = VDWParameters::Type::LennardJonesSecondOrderTaylorShifted;
    }
  }
  forceField.preComputeDerivedParameters();
  forceField.preComputePotentialShift();
}

// A periodic box of flexible united-atom heptane (harmonic bonds and bends, TraPPE torsions) placed whole
// on a grid as strained zig-zag chains, so bonded and non-bonded forces are both non-zero. The Lennard-Jones
// potential is Taylor-shifted at the cutoff so the finite differences are smooth, and the tail corrections
// are off: the analytic tail pressure is the virial of an unshifted truncated potential (it contains the
// impulsive term (2 pi / 3) rho^2 r_c^3 u(r_c)), which -dU_tail/dV of a shifted potential does not have.
System makeStrainedHeptaneBox(double boxLength, std::size_t nPerDim, RandomNumber &random)
{
  ForceField forceField = ForceField({{"CH3", false, 15.03452, 0.0, 0.0, 8, false},
                                      {"CH2", false, 14.02658, 0.0, 0.0, 8, false}},
                                     {{98.0, 3.75}, {46.0, 3.95}}, ForceField::MixingRule::Lorentz_Berthelot, 12.0,
                                     12.0, 12.0, false, false, false);
  forceField.useCharge = false;
  forceField.omitEwaldFourier = true;
  useSecondOrderTaylorShiftedLennardJones(forceField);

  ConnectivityTable connectivityTable(7);
  Potentials::IntraMolecularPotentials potentials{};
  for (std::size_t i = 0; i < 6; ++i)
  {
    connectivityTable[i, i + 1] = true;
    connectivityTable[i + 1, i] = true;
    potentials.bonds.emplace_back(std::array<std::size_t, 2>{i, i + 1}, BondType::Harmonic,
                                  std::vector<double>{96500.0, 1.54});
  }
  for (std::size_t i = 0; i < 5; ++i)
  {
    potentials.bends.emplace_back(std::array<std::size_t, 3>{i, i + 1, i + 2}, BendType::Harmonic,
                                  std::vector<double>{62500.0, 114.0});
  }
  for (std::size_t i = 0; i < 4; ++i)
  {
    potentials.torsions.emplace_back(std::array<std::size_t, 4>{i, i + 1, i + 2, i + 3}, TorsionType::TraPPE,
                                     std::vector<double>{0.0, 355.03, -68.19, 791.32});
  }

  std::vector<Atom> atoms;
  for (std::size_t a = 0; a < 7; ++a)
  {
    const std::uint16_t type = (a == 0 || a == 6) ? 0 : 1;
    atoms.push_back(Atom({0.0, 0.0, 0.0}, 0.0, 1.0, 0, type, 0, false, false));
  }
  Component heptane(forceField, "heptane", 540.13, 2736000.0, 0.349, atoms, connectivityTable, potentials, 5, 21);

  const std::size_t nMolecules = nPerDim * nPerDim * nPerDim;
  System system = System(forceField, SimulationBox(boxLength, boxLength, boxLength), false, 300.0, 1e5, 1.0, {},
                         {heptane}, {}, {nMolecules}, 5);

  std::span<Atom> atomData = system.spanOfMoleculeAtoms();
  const double spacing = boxLength / static_cast<double>(nPerDim);
  std::size_t index = 0;
  for (std::size_t ix = 0; ix < nPerDim; ++ix)
    for (std::size_t iy = 0; iy < nPerDim; ++iy)
      for (std::size_t iz = 0; iz < nPerDim; ++iz)
      {
        const double3 origin((static_cast<double>(ix) + 0.5) * spacing, (static_cast<double>(iy) + 0.5) * spacing,
                             (static_cast<double>(iz) + 0.5) * spacing);
        for (std::size_t a = 0; a < 7; ++a)
        {
          // a zig-zag with 1.54 A bonds and ~114 degree bends, strained by a random displacement
          const double z = (static_cast<double>(a) - 3.0) * 1.29;
          const double x = (a % 2 == 0) ? 0.0 : 0.84;
          const double3 noise(random.uniform() - 0.5, random.uniform() - 0.5, random.uniform() - 0.5);
          atomData[index * 7 + a].position = origin + double3(x, 0.0, z) + 0.15 * noise;
        }
        ++index;
      }
  return system;
}

// Total potential energy (inter-molecular, tail and intra-molecular) of a configuration in a box.
double totalEnergy(const System &system, const SimulationBox &box, std::span<const Atom> atoms)
{
  double energy = (Interactions::computeInterMolecularEnergy(system.forceField, box, atoms) +
                   Interactions::computeInterMolecularTailEnergy(system.forceField, box, atoms))
                      .potentialEnergy();
  for (const Molecule &molecule : system.moleculeData)
  {
    const Component &component = system.components[molecule.componentId];
    energy += component.intraMolecularPotentials
                  .computeInternalEnergies(system.forceField, box, atoms.subspan(molecule.atomIndex, molecule.numberOfAtoms))
                  .potentialEnergy();
  }
  return energy;
}

// -dU/dV under the affine scaling of every atom with the box (the coupling of the barostat to a flexible
// molecule); the molecules are whole, so scaling the raw positions is the scaling of the cell.
double atomicScalingExcessPressure(const System &system)
{
  const double delta = 1e-5;
  const SimulationBox boxPlus = system.simulationBox.scaled(std::cbrt(1.0 + delta));
  const SimulationBox boxMinus = system.simulationBox.scaled(std::cbrt(1.0 - delta));
  const std::span<const Atom> atoms = std::as_const(system).spanOfMoleculeAtoms();
  std::vector<Atom> plus(atoms.begin(), atoms.end());
  std::vector<Atom> minus(atoms.begin(), atoms.end());
  for (Atom &atom : plus) atom.position = std::cbrt(1.0 + delta) * atom.position;
  for (Atom &atom : minus) atom.position = std::cbrt(1.0 - delta) * atom.position;
  return -(totalEnergy(system, boxPlus, plus) - totalEnergy(system, boxMinus, minus)) /
         (boxPlus.volume - boxMinus.volume);
}

// -dU/dV under the scaling of the centers of mass only (the coupling of the Monte Carlo volume move); the
// intra-molecular energy does not change.
double centerOfMassScalingExcessPressure(System &system)
{
  const double delta = 1e-5;
  const SimulationBox boxPlus = system.simulationBox.scaled(std::cbrt(1.0 + delta));
  const SimulationBox boxMinus = system.simulationBox.scaled(std::cbrt(1.0 - delta));
  const auto plus = system.scaledCenterOfMassPositions(system.simulationBox, boxPlus);
  const auto minus = system.scaledCenterOfMassPositions(system.simulationBox, boxMinus);
  return -(totalEnergy(system, boxPlus, plus.second) - totalEnergy(system, boxMinus, minus.second)) /
         (boxPlus.volume - boxMinus.volume);
}

void computeGradients(System &system)
{
  Integrators::updateGradients(
      system.moleculeData, system.spanOfMoleculeAtoms(), system.spanOfMoleculeDynamics(), system.spanOfFrameworkAtoms(),
      system.forceField, system.simulationBox, system.components, system.eik_x, system.eik_y, system.eik_z,
      system.eik_xy, system.trialEik, system.fixedFrameworkStoredEik, system.interpolationGrids,
      system.numberOfMoleculesPerComponent, system.framework, system.spanOfFrameworkDynamics(), &system.crossLinks);
}
}  // namespace

// For flexible molecules the barostat virial is the atomic virial: its trace / 3V must be the finite
// difference -dU/dV under the affine scaling of all atoms, while the molecular virial reproduces the
// center-of-mass scaling. The two differ by the bonded and the intra-molecular part of the virial, which
// is far from zero for strained chains.
TEST(BAROSTAT_VIRIAL, flexible_heptane_matches_atomic_scaling_finite_difference)
{
  RandomNumber random(17);
  System system = makeStrainedHeptaneBox(36.0, 3, random);
  computeGradients(system);

  const double volume = system.simulationBox.volume;
  const double3x3 molecularVirial = system.computeMolecularPressure().second;
  const double3x3 barostatVirial = computeBarostatVirial(system, molecularVirial, BarostatCoupling::Atomic);

  const double molecularExcessPressure = molecularVirial.trace() / (3.0 * volume);
  const double barostatExcessPressure = barostatVirial.trace() / (3.0 * volume);
  const double fdCenterOfMass = centerOfMassScalingExcessPressure(system);
  const double fdAtomic = atomicScalingExcessPressure(system);

  EXPECT_NEAR(molecularExcessPressure, fdCenterOfMass, 1e-6 * std::abs(fdCenterOfMass) + 1e-8)
      << "molecular virial disagrees with the center-of-mass scaling -dU/dV";
  EXPECT_NEAR(barostatExcessPressure, fdAtomic, 1e-6 * std::abs(fdAtomic) + 1e-8)
      << "barostat virial disagrees with the atomic scaling -dU/dV";

  // The distinction matters: the bonded virial of the strained chains is large.
  EXPECT_GT(std::abs(barostatExcessPressure - molecularExcessPressure), 0.1 * std::abs(molecularExcessPressure))
      << "the atomic and molecular virials hardly differ; the test configuration is not strained";

  // The off-diagonal elements follow the same conversion; the result is a symmetric tensor.
  EXPECT_NEAR(barostatVirial.ay, barostatVirial.bx, 1e-8 * std::abs(barostatVirial.trace()));
  EXPECT_NEAR(barostatVirial.az, barostatVirial.cx, 1e-8 * std::abs(barostatVirial.trace()));
  EXPECT_NEAR(barostatVirial.bz, barostatVirial.cy, 1e-8 * std::abs(barostatVirial.trace()));

  // Molecular coupling drives the centres of mass: its virial is the molecular one (symmetrized).
  const double3x3 molecularCoupling = computeBarostatVirial(system, molecularVirial, BarostatCoupling::Molecular);
  EXPECT_NEAR(molecularCoupling.trace(), molecularVirial.trace(), 1e-10 * std::abs(molecularVirial.trace()));
  EXPECT_NEAR(molecularCoupling.ay, 0.5 * (molecularVirial.ay + molecularVirial.bx),
              1e-10 * std::abs(molecularVirial.trace()));
}

// The kinetic partner of the barostat virial counts every flexible atom, so for a flexible liquid the
// pair (barostat virial, atomic kinetic energy) and the pair (molecular virial, molecular kinetic energy)
// are two estimators of the same pressure whose ideal-gas parts differ by (N_atoms - N_molecules) k T / V.
// The mixed pair is off by exactly that amount; make sure the mixed one is not what a flexible system
// would get: the kinetic virial of a flexible molecule must be the atomic one.
TEST(BAROSTAT_VIRIAL, kinetic_virial_of_flexible_molecules_is_atomic)
{
  RandomNumber random(5);
  System system = makeStrainedHeptaneBox(36.0, 2, random);
  std::span<AtomDynamics> dynamics = system.spanOfMoleculeDynamics();
  std::span<const Atom> atoms = std::as_const(system).spanOfMoleculeAtoms();
  double expected = 0.0;
  for (std::size_t i = 0; i < dynamics.size(); ++i)
  {
    dynamics[i].velocity = double3(random.uniform() - 0.5, random.uniform() - 0.5, random.uniform() - 0.5);
    const double mass = system.forceField.pseudoAtoms[atoms[i].type].mass;
    expected += mass * double3::dot(dynamics[i].velocity, dynamics[i].velocity);
  }
  EXPECT_NEAR(computeMolecularKineticVirial(system, BarostatCoupling::Atomic).trace(), expected, 1e-10 * expected);

  // Molecular coupling: the centre-of-mass momentum flux sum_I M_I V_I outer V_I, smaller than the atomic one.
  double expectedMolecular = 0.0;
  for (const Molecule &molecule : system.moleculeData)
  {
    double mass = 0.0;
    double3 momentum{};
    for (std::size_t k = 0; k < molecule.numberOfAtoms; ++k)
    {
      const double atomMass = system.forceField.pseudoAtoms[atoms[molecule.atomIndex + k].type].mass;
      mass += atomMass;
      momentum += atomMass * dynamics[molecule.atomIndex + k].velocity;
    }
    expectedMolecular += double3::dot(momentum, momentum) / mass;
  }
  EXPECT_NEAR(computeMolecularKineticVirial(system, BarostatCoupling::Molecular).trace(), expectedMolecular,
              1e-10 * expectedMolecular);
  EXPECT_LT(expectedMolecular, expected);
}

// A rigid molecule has a single coupled point, its center of mass: the barostat virial is the molecular one.
TEST(BAROSTAT_VIRIAL, rigid_molecules_leave_the_molecular_virial_unchanged)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  useSecondOrderTaylorShiftedLennardJones(forceField);
  forceField.useCharge = false;
  forceField.omitEwaldFourier = true;
  Component co2 = Component::makeCO2(forceField, 0, true);
  ASSERT_TRUE(co2.rigid);
  System system = System(forceField, SimulationBox(24.0, 24.0, 24.0), false, 300.0, 1e4, 1.0, {}, {co2}, {}, {10}, 5);

  RandomNumber random(3);
  std::span<Atom> atomData = system.spanOfMoleculeAtoms();
  for (Molecule &molecule : system.moleculeData)
  {
    molecule.centerOfMassPosition = 24.0 * double3(random.uniform(), random.uniform(), random.uniform());
    molecule.orientation = random.randomSimdQuatd();
    const double3x3 rotation = double3x3::buildRotationMatrixInverse(molecule.orientation);
    for (std::size_t k = 0; k < molecule.numberOfAtoms; ++k)
      atomData[molecule.atomIndex + k].position = molecule.centerOfMassPosition + rotation * co2.atoms[k].position;
  }
  computeGradients(system);

  const double3x3 molecularVirial = system.computeMolecularPressure().second;
  const double3x3 barostatVirial = computeBarostatVirial(system, molecularVirial, BarostatCoupling::Atomic);
  EXPECT_DOUBLE_EQ(barostatVirial.ax, molecularVirial.ax);
  EXPECT_DOUBLE_EQ(barostatVirial.by, molecularVirial.by);
  EXPECT_DOUBLE_EQ(barostatVirial.cz, molecularVirial.cz);
  EXPECT_DOUBLE_EQ(barostatVirial.ay, molecularVirial.ay);
  EXPECT_DOUBLE_EQ(barostatVirial.bz, molecularVirial.bz);
  EXPECT_DOUBLE_EQ(barostatVirial.az, molecularVirial.az);
}

TEST(BAROSTAT_VIRIAL, flexible_framework_bond_matches_atomic_scaling_finite_difference)
{
  ForceField forceField({{"A", false, 12.0, 0.0, 0.0, 6, false}}, {{0.0, 3.0}},
                        ForceField::MixingRule::Lorentz_Berthelot, 8.0, 8.0, 8.0, false, false, false);
  forceField.useCharge = false;
  forceField.omitEwaldFourier = true;
  const SimulationBox box(20.0, 20.0, 20.0);
  std::vector<Atom> atoms{Atom({5.0, 5.0, 5.0}, 0.0, 1.0, 0, 0, 0, 0, true),
                          Atom({6.4, 5.0, 5.0}, 0.0, 1.0, 0, 0, 0, 0, true)};
  Framework framework(forceField, "two-atom", box, 1, atoms, atoms, int3(1, 1, 1));
  framework.rigid = false;
  framework.intraMolecularPotentials.bonds = {BondPotential({0, 1}, BondType::Harmonic, {20000.0, 1.2})};
  framework.intraMolecularImageShifts.bonds = {{{int3{}, int3{}}}};

  System system(forceField, box, false, 300.0, 1.0e4, 1.0, framework, {}, {}, {}, 5);
  ASSERT_TRUE(system.framework->hasMobileAtoms());

  const double3x3 virial = system.computeMolecularPressure().second;
  const double volume = system.simulationBox.volume;
  const double excess = virial.trace() / (3.0 * volume);

  const double delta = 1e-5;
  const auto energyAtScale = [&](double scale)
  {
    SimulationBox scaledBox = system.simulationBox.scaled(scale);
    std::vector<Atom> scaled(system.spanOfFrameworkAtoms().begin(), system.spanOfFrameworkAtoms().end());
    for (Atom& atom : scaled) atom.position = scale * atom.position;
    return Interactions::computeFrameworkIntraMolecularEnergy(system.forceField, *system.framework, scaledBox, scaled)
        .potentialEnergy();
  };
  const double scalePlus = std::cbrt(1.0 + delta);
  const double scaleMinus = std::cbrt(1.0 - delta);
  const double fd =
      -(energyAtScale(scalePlus) - energyAtScale(scaleMinus)) /
      (system.simulationBox.scaled(scalePlus).volume - system.simulationBox.scaled(scaleMinus).volume);

  EXPECT_NEAR(excess, fd, 1e-5 * std::abs(fd) + 1e-8);
  EXPECT_GT(std::abs(fd), 1e-4);
}

TEST(BAROSTAT_VIRIAL, kinetic_virial_of_mixed_framework_uses_group_velocity)
{
  ForceField forceField({{"C", true, 12.0, 0.0, 0.0, 6, true}}, {{80.0, 3.0}},
                        ForceField::MixingRule::Lorentz_Berthelot, 5.0, 5.0, 5.0, false, false, false);
  const SimulationBox box(14.0, 14.0, 14.0);
  std::vector<Atom> atoms{Atom({4.0 / 14.0, 0.5, 0.5}, 0.0, 1.0, 0, 0, 0, 0, true),
                          Atom({5.4 / 14.0, 0.52, 0.5}, 0.0, 1.0, 0, 0, 0, 0, true),
                          Atom({6.8 / 14.0, 0.48, 0.5}, 0.0, 1.0, 0, 0, 0, 0, true),
                          Atom({8.2 / 14.0, 0.52, 0.5}, 0.0, 1.0, 0, 0, 0, 0, true),
                          Atom({9.6 / 14.0, 0.5, 0.5}, 0.0, 1.0, 0, 0, 0, 0, true),
                          Atom({11.0 / 14.0, 0.51, 0.5}, 0.0, 1.0, 0, 0, 0, 0, true)};
  Framework framework(forceField, "mixed-chain", box, 1, atoms, atoms, int3(1, 1, 1));
  const std::filesystem::path definition =
      std::filesystem::temp_directory_path() / "raspa3-mixed-framework-barostat.json";
  {
    std::ofstream stream(definition);
    stream << R"({"Type": "Flexible",
                  "Groups": [
                    {"Type": "Fixed", "Atoms": [0]},
                    {"Type": "Rigid", "Atoms": [1, 2, 3]},
                    {"Type": "Flexible", "Atoms": [4, 5]}
                  ]})";
  }
  framework.readFrameworkDefinition(forceField, definition.string());
  std::filesystem::remove(definition);

  System system(forceField, box, false, 300.0, 1.0e4, 1.0, {framework}, {}, {}, {}, 5);
  system.initializeFrameworkGroupData();
  ASSERT_EQ(system.frameworkGroupData.size(), 1u);

  for (AtomDynamics& dynamics : system.spanOfFrameworkDynamics()) dynamics.velocity = double3(0.0, 0.0, 0.0);
  system.frameworkGroupData[0].velocity = double3(0.3, -0.2, 0.4);
  system.spanOfFrameworkDynamics()[4].velocity = double3(0.1, 0.0, -0.05);
  system.spanOfFrameworkDynamics()[5].velocity = double3(-0.2, 0.15, 0.0);

  const double groupMass = system.framework->groups[1].mass;
  const double mass4 = system.forceField.pseudoAtoms[system.spanOfFrameworkAtoms()[4].type].mass;
  const double mass5 = system.forceField.pseudoAtoms[system.spanOfFrameworkAtoms()[5].type].mass;
  const double expected =
      groupMass * double3::dot(system.frameworkGroupData[0].velocity, system.frameworkGroupData[0].velocity) +
      mass4 * double3::dot(system.spanOfFrameworkDynamics()[4].velocity, system.spanOfFrameworkDynamics()[4].velocity) +
      mass5 * double3::dot(system.spanOfFrameworkDynamics()[5].velocity, system.spanOfFrameworkDynamics()[5].velocity);

  EXPECT_NEAR(computeMolecularKineticVirial(system, BarostatCoupling::Atomic).trace(), expected, 1e-10 * expected);
}

TEST(BAROSTAT_VIRIAL, mobile_framework_ewald_enters_the_pressure_tensor)
{
  ForceField forceField({{"C", true, 12.0, 0.4, 0.0, 6, true}}, {{0.0, 3.0}},
                        ForceField::MixingRule::Lorentz_Berthelot, 8.0, 8.0, 8.0, true, false, true);
  const SimulationBox box(16.0, 16.0, 16.0);
  std::vector<Atom> atoms{Atom({5.0, 8.0, 8.0}, 0.4, 1.0, 0, 0, 0, 0, true),
                          Atom({9.0, 8.0, 8.0}, -0.4, 1.0, 0, 0, 0, 0, true)};
  Framework framework(forceField, "charged-pair", box, 1, atoms, atoms, int3(1, 1, 1));
  framework.rigid = false;

  System system(forceField, box, false, 300.0, 1.0e4, 1.0, framework, {}, {}, {}, 5);
  system.forceField.initializeEwaldParameters(system.simulationBox);
  system.precomputeTotalRigidEnergy();
  ASSERT_TRUE(system.framework->hasMobileAtoms());

  const double3x3 virial = system.computeMolecularPressure().second;
  EXPECT_GT(std::abs(virial.trace()), 0.0);

  const auto ewaldEnergy = [&](double scale)
  {
    System copy = system;
    copy.simulationBox = system.simulationBox.scaled(scale);
    copy.forceField.initializeEwaldParameters(copy.simulationBox);
    copy.precomputeTotalRigidEnergy();
    for (Atom& atom : copy.spanOfFrameworkAtoms()) atom.position = scale * atom.position;
    return Integrators::updateGradients(copy.moleculeData, copy.spanOfMoleculeAtoms(), copy.spanOfMoleculeDynamics(),
                                        copy.spanOfFrameworkAtoms(), copy.forceField, copy.simulationBox,
                                        copy.components, copy.eik_x, copy.eik_y, copy.eik_z, copy.eik_xy, copy.trialEik,
                                        copy.fixedFrameworkStoredEik, copy.interpolationGrids,
                                        copy.numberOfMoleculesPerComponent, copy.framework,
                                        copy.spanOfFrameworkDynamics())
        .potentialEnergy();
  };
  const double delta = 1e-5;
  const double scalePlus = std::cbrt(1.0 + delta);
  const double scaleMinus = std::cbrt(1.0 - delta);
  const double fd = -(ewaldEnergy(scalePlus) - ewaldEnergy(scaleMinus)) /
                    (system.simulationBox.scaled(scalePlus).volume - system.simulationBox.scaled(scaleMinus).volume);
  EXPECT_NEAR(virial.trace() / (3.0 * system.simulationBox.volume), fd, 5e-5 * std::abs(fd) + 1e-8);
}

namespace
{
// A bonded, charged, flexible framework chain (bonds 0-1 and 1-2, so the pairs 0-1, 1-2 and 0-2 are excluded
// from the Coulomb pair sum): the exclusion correction of the excluded pairs is distance dependent and must
// enter the pressure tensor. Checks the analytic virial against the finite difference of the total energy of
// Integrators::updateGradients under the affine scaling of the framework atoms.
void checkBondedFrameworkChargeStrain(ForceField forceField)
{
  const SimulationBox box(16.0, 16.0, 16.0);
  // fractional unit-cell positions: a bent chain with 1.5 A bonds
  std::vector<Atom> atoms{Atom({5.0 / 16.0, 0.5, 0.5}, 0.4, 1.0, 0, 0, 0, 0, true),
                          Atom({6.5 / 16.0, 0.5, 0.5}, -0.8, 1.0, 0, 0, 0, 0, true),
                          Atom({7.3 / 16.0, 9.3 / 16.0, 0.5}, 0.4, 1.0, 0, 0, 0, 0, true)};
  Framework framework(forceField, "charged-chain", box, 1, atoms, atoms, int3(1, 1, 1));
  framework.rigid = false;
  framework.connectivityTable = ConnectivityTable(3);
  framework.connectivityTable[0, 1] = framework.connectivityTable[1, 0] = true;
  framework.connectivityTable[1, 2] = framework.connectivityTable[2, 1] = true;

  System system(forceField, box, false, 300.0, 1.0e4, 1.0, framework, {}, {}, {}, 5);
  system.forceField.initializeEwaldParameters(system.simulationBox);
  system.precomputeTotalRigidEnergy();
  ASSERT_TRUE(system.framework->hasMobileAtoms());

  const auto energyAtScale = [&](double scale)
  {
    System copy = system;
    copy.simulationBox = system.simulationBox.scaled(scale);
    copy.forceField.initializeEwaldParameters(copy.simulationBox);
    copy.precomputeTotalRigidEnergy();
    for (Atom& atom : copy.spanOfFrameworkAtoms()) atom.position = scale * atom.position;
    return Integrators::updateGradients(copy.moleculeData, copy.spanOfMoleculeAtoms(), copy.spanOfMoleculeDynamics(),
                                        copy.spanOfFrameworkAtoms(), copy.forceField, copy.simulationBox,
                                        copy.components, copy.eik_x, copy.eik_y, copy.eik_z, copy.eik_xy, copy.trialEik,
                                        copy.fixedFrameworkStoredEik, copy.interpolationGrids,
                                        copy.numberOfMoleculesPerComponent, copy.framework,
                                        copy.spanOfFrameworkDynamics())
        .potentialEnergy();
  };
  const double delta = 1e-5;
  const double scalePlus = std::cbrt(1.0 + delta);
  const double scaleMinus = std::cbrt(1.0 - delta);
  const double fd = -(energyAtScale(scalePlus) - energyAtScale(scaleMinus)) /
                    (system.simulationBox.scaled(scalePlus).volume - system.simulationBox.scaled(scaleMinus).volume);
  ASSERT_GT(std::abs(fd), 1e-6) << "the excluded pairs do not strain the energy; the test configuration is wrong";

  const double3x3 virial = system.computeMolecularPressure().second;
  EXPECT_NEAR(virial.trace() / (3.0 * system.simulationBox.volume), fd, 5e-5 * std::abs(fd) + 1e-8);
}
}  // namespace

TEST(BAROSTAT_VIRIAL, mobile_framework_ewald_exclusion_enters_the_pressure_tensor)
{
  ForceField forceField({{"C", true, 12.0, 0.4, 0.0, 6, true}}, {{0.0, 3.0}}, ForceField::MixingRule::Lorentz_Berthelot,
                        8.0, 8.0, 8.0, true, false, true);
  checkBondedFrameworkChargeStrain(forceField);
}

TEST(BAROSTAT_VIRIAL, mobile_framework_wolf_exclusion_enters_the_pressure_tensor)
{
  ForceField forceField({{"C", true, 12.0, 0.4, 0.0, 6, true}}, {{0.0, 3.0}}, ForceField::MixingRule::Lorentz_Berthelot,
                        8.0, 8.0, 8.0, true, false, true);
  forceField.chargeMethod = ForceField::ChargeMethod::Wolf;
  forceField.omitEwaldFourier = true;
  forceField.EwaldAlpha = 0.25;
  checkBondedFrameworkChargeStrain(forceField);
}

TEST(BAROSTAT_VIRIAL, propagate_cell_refreshes_fixed_framework_ewald)
{
  ForceField forceField({{"C", true, 12.0, 0.4, 0.0, 6, true}}, {{0.0, 3.0}},
                        ForceField::MixingRule::Lorentz_Berthelot, 8.0, 8.0, 8.0, true, false, true);
  const SimulationBox box(24.0, 24.0, 24.0);
  std::vector<Atom> atoms{Atom({5.0, 8.0, 8.0}, 0.4, 1.0, 0, 0, 0, 0, true),
                          Atom({9.0, 8.0, 8.0}, -0.4, 1.0, 0, 0, 0, 0, true)};
  Framework framework(forceField, "rigid-charged-pair", box, 1, atoms, atoms, int3(1, 1, 1));
  framework.rigid = true;
  System system(forceField, box, false, 300.0, 1.0e4, 1.0, framework, {}, {}, {}, 5);
  system.timeStep = 1.0e-3;
  system.thermobarostat =
      Thermobarostat(MolecularDynamicsEnsemble::NPT, CellMinimizationType::Isotropic, MonoclinicAngleType::Beta, 300.0,
                     system.pressure, system.timeStep, system.translationalDegreesOfFreedom, 3, 1, 1.0);
  RandomNumber random(4);
  system.thermobarostat->initialize(random);
  system.thermobarostat->logVolumeVelocity = 0.4;
  system.precomputeTotalGradients();
  system.precomputeTotalRigidEnergy();
  system.runningEnergies.translationalKineticEnergy = 0.0;
  system.runningEnergies.thermobarostatEnergy = system.thermobarostat->energy(system.simulationBox.volume);

  molecularDynamicsStep(system);

  const EwaldStructureFactor afterStep = system.fixedFrameworkStoredEik;
  system.precomputeTotalRigidEnergy();
  ASSERT_EQ(afterStep.size(), system.fixedFrameworkStoredEik.size());
  ASSERT_GT(afterStep.size(), 0uz);
  for (std::size_t k = 0; k < afterStep.size(); ++k)
  {
    EXPECT_NEAR(afterStep.sum[k].real(), system.fixedFrameworkStoredEik.sum[k].real(), 1e-10);
    EXPECT_NEAR(afterStep.sum[k].imag(), system.fixedFrameworkStoredEik.sum[k].imag(), 1e-10);
  }
}

// With "CutOff": "auto" the barostat moves the VDW cut-off with the cell at every step; the potential shifts and
// the tail corrections must be those of the new cut-off, not of the cut-off the force field was set up with.
TEST(BAROSTAT_VIRIAL, propagate_cell_follows_the_automatic_cutoff)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, true, false);
  forceField.useCharge = false;
  forceField.omitEwaldFourier = true;
  forceField.cutOffCoulomb = 8.0;
  forceField.cutOffFrameworkVDWAutomatic = true;
  forceField.cutOffMoleculeVDWAutomatic = true;

  Component methane =
      Component(forceField, "methane", 190.564, 4599200, 0.01142, {Atom({0, 0, 0}, 0.0, 1.0, 0, 6, 0, false, false)},
                {}, {}, 5, 21, {}, std::nullopt, false);
  const SimulationBox box(20.0, 20.0, 20.0);
  System system(forceField, box, false, 300.0, 1.0e4, 1.0, {}, {methane}, {}, {30}, 5);

  const double setupCutOff = system.forceField.cutOffMoleculeVDW;
  EXPECT_DOUBLE_EQ(setupCutOff, 10.0 - std::numeric_limits<double>::epsilon());
  const double setupShift = system.forceField(6, 6).shift;
  const double setupTail = system.forceField(6, 6).tailCorrectionEnergy;

  system.timeStep = 1.0e-3;
  system.thermobarostat =
      Thermobarostat(MolecularDynamicsEnsemble::NPT, CellMinimizationType::Isotropic, MonoclinicAngleType::Beta, 300.0,
                     system.pressure, system.timeStep, system.translationalDegreesOfFreedom, 3, 1, 1.0);
  RandomNumber random(4);
  system.thermobarostat->initialize(random);
  system.thermobarostat->logVolumeVelocity = 0.4;
  system.precomputeTotalGradients();
  system.precomputeTotalRigidEnergy();
  system.runningEnergies.translationalKineticEnergy = 0.0;
  system.runningEnergies.thermobarostatEnergy = system.thermobarostat->energy(system.simulationBox.volume);

  molecularDynamicsStep(system);

  // the cell changed and the cut-off followed it
  ASSERT_NE(system.simulationBox.lengthA, 20.0);
  const double3 widths = system.simulationBox.perpendicularWidths();
  const double expectedCutOff = 0.5 * std::min({widths.x, widths.y, widths.z}) - std::numeric_limits<double>::epsilon();
  EXPECT_DOUBLE_EQ(system.forceField.cutOffMoleculeVDW, expectedCutOff);
  EXPECT_DOUBLE_EQ(system.forceField.cutOffFrameworkVDW, expectedCutOff);
  EXPECT_NE(system.forceField.cutOffMoleculeVDW, setupCutOff);

  // the shifts and tail corrections of every pair belong to the new cut-off
  ForceField reference = system.forceField;
  reference.preComputeCutOffDependentParameters();
  for (std::size_t pair = 0; pair < reference.data.size(); ++pair)
  {
    EXPECT_DOUBLE_EQ(system.forceField.data[pair].shift, reference.data[pair].shift);
    EXPECT_DOUBLE_EQ(system.forceField.data[pair].tailCorrectionEnergy, reference.data[pair].tailCorrectionEnergy);
    EXPECT_DOUBLE_EQ(system.forceField.data[pair].tailCorrectionPressure, reference.data[pair].tailCorrectionPressure);
  }
  EXPECT_NE(system.forceField(6, 6).shift, setupShift);
  EXPECT_NE(system.forceField(6, 6).tailCorrectionEnergy, setupTail);
}

namespace
{
class TemporaryDefinition
{
 public:
  explicit TemporaryDefinition(std::string_view contents)
      : path(std::filesystem::temp_directory_path() /
             std::format("raspa3-barostat-mobile-framework-ewald-{}.json", static_cast<long long>(::getpid())))
  {
    // one file per process: the tests of this file run in parallel and used to overwrite and delete each other's
    std::ofstream(path) << contents;
  }
  ~TemporaryDefinition()
  {
    std::error_code ignored;
    std::filesystem::remove(path, ignored);
  }
  std::filesystem::path path;
};

// A five-atom bonded chain with charges 'charges' and a charged single-site adsorbate. 'mixed' makes the atoms 0 and
// 1 lab-fixed and the others flexible; otherwise the whole chain is flexible. No Lennard-Jones and no bonded terms.
System makeChargedChainSystem(bool mixed, const std::array<double, 5>& charges, double ionCharge)
{
  ForceField forceField({{"C", true, 12.0, 0.0, 0.0, 6, true}, {"M", false, 16.0, 0.0, 0.0, 8, false}},
                        {{0.0, 3.0}, {0.0, 3.0}}, ForceField::MixingRule::Lorentz_Berthelot, 7.5, 7.5, 7.5, true, false,
                        true);
  const SimulationBox box(16.0, 16.0, 16.0);
  std::vector<Atom> atoms;
  const std::array<double3, 5> positions{double3(5.0, 8.0, 8.0), double3(6.4, 8.3, 8.0), double3(7.8, 7.8, 8.0),
                                         double3(9.2, 8.3, 8.0), double3(10.6, 8.0, 8.0)};
  for (std::size_t i = 0; i != 5; ++i)
  {
    atoms.push_back(Atom(positions[i] / 16.0, charges[i], 1.0, 0, 0, 0, 0, true));
  }
  Framework framework(forceField, "charged-chain", box, 1, atoms, atoms, int3(1, 1, 1));

  const TemporaryDefinition definition(mixed ? R"({"Type": "Flexible",
                                                   "Groups": [{"Type": "Fixed", "Atoms": [0, 1]},
                                                              {"Type": "Flexible", "Atoms": [2, 3, 4]}],
                                                   "Intra14VanDerWaalsScalingValue": 0.0,
                                                   "Intra14ChargeChargeScalingValue": 0.0})"
                                             : R"({"Type": "Flexible",
                                                   "Intra14VanDerWaalsScalingValue": 0.0,
                                                   "Intra14ChargeChargeScalingValue": 0.0})");
  framework.readFrameworkDefinition(forceField, definition.path.string());
  EXPECT_EQ(framework.isMixed(), mixed);

  Component ion(forceField, "ion", 16.0, 1.0e6, 0.011, {Atom({0.0, 0.0, 0.0}, ionCharge, 1.0, 0, 1, 0, false, false)},
                ConnectivityTable(1), Potentials::IntraMolecularPotentials{}, 0, 0);
  System system(forceField, box, false, 300.0, 1.0e4, 1.0, {framework}, {ion}, {}, {1}, 5);
  system.spanOfMoleculeAtoms()[0].position = double3(9.0, 11.0, 9.5);
  system.forceField.initializeEwaldParameters(system.simulationBox);
  system.precomputeTotalRigidEnergy();
  return system;
}

RunningEnergy gradientEnergy(System& system)
{
  return Integrators::updateGradients(system.moleculeData, system.spanOfMoleculeAtoms(),
                                      system.spanOfMoleculeDynamics(), system.spanOfFrameworkAtoms(), system.forceField,
                                      system.simulationBox, system.components, system.eik_x, system.eik_y, system.eik_z,
                                      system.eik_xy, system.trialEik, system.fixedFrameworkStoredEik,
                                      system.interpolationGrids, system.numberOfMoleculesPerComponent, system.framework,
                                      system.spanOfFrameworkDynamics());
}

// The system with the cell and every atom scaled affinely.
System scaledSystem(const System& system, double scale)
{
  System copy = system;
  copy.simulationBox = system.simulationBox.scaled(scale);
  copy.forceField.initializeEwaldParameters(copy.simulationBox);
  for (Atom& atom : copy.spanOfFrameworkAtoms()) atom.position = scale * atom.position;
  for (Atom& atom : copy.spanOfMoleculeAtoms()) atom.position = scale * atom.position;
  copy.precomputeTotalRigidEnergy();
  return copy;
}

const std::array<double, 5> chainCharges{0.3, -0.3, 0.2, -0.4, 0.2};
}  // namespace

// The running energy of a Monte Carlo run and of the integrator are the same energy, for a mobile framework as well.
TEST(BAROSTAT_VIRIAL, total_energy_matches_the_gradient_energy_of_a_mobile_framework)
{
  for (const bool mixed : {false, true})
  {
    System system = makeChargedChainSystem(mixed, chainCharges, 0.4);
    const RunningEnergy total = system.computeTotalEnergies();
    const RunningEnergy gradient = gradientEnergy(system);
    EXPECT_NEAR(total.ewald_fourier, gradient.ewald_fourier, 1e-8) << "mixed " << mixed;
    EXPECT_NEAR(total.ewald_self, gradient.ewald_self, 1e-8) << "mixed " << mixed;
    EXPECT_NEAR(total.ewald_exclusion, gradient.ewald_exclusion, 1e-8) << "mixed " << mixed;
    EXPECT_NEAR(total.potentialEnergy(), gradient.potentialEnergy(), 1e-8) << "mixed " << mixed;
  }
}

// A Monte Carlo move of a charge next to a mobile framework: the incremental Fourier energy is the change of the total.
TEST(BAROSTAT_VIRIAL, ewald_move_difference_sees_a_mobile_framework)
{
  for (const bool mixed : {false, true})
  {
    System system = makeChargedChainSystem(mixed, chainCharges, 0.4);
    const RunningEnergy before = system.computeTotalEnergies();

    const std::vector<Atom> oldAtoms(system.spanOfMoleculeAtoms().begin(), system.spanOfMoleculeAtoms().end());
    std::vector<Atom> newAtoms = oldAtoms;
    newAtoms[0].position = double3(8.0, 9.5, 8.5);

    const RunningEnergy difference = Interactions::energyDifferenceEwaldFourier(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.trialEik, system.forceField,
        system.simulationBox, system.components, newAtoms, oldAtoms);

    system.spanOfMoleculeAtoms()[0].position = newAtoms[0].position;
    const RunningEnergy after = system.computeTotalEnergies();
    ASSERT_GT(std::abs(difference.ewald_fourier), 1e-3) << "mixed " << mixed;
    EXPECT_NEAR(difference.ewald_fourier, after.ewald_fourier - before.ewald_fourier, 1e-8) << "mixed " << mixed;
  }
}

// The framework-adsorbate entry of the energy status includes the mobile framework atoms.
TEST(BAROSTAT_VIRIAL, energy_status_holds_the_mobile_framework_fourier_energy)
{
  const auto fourier = [](const std::array<double, 5>& charges, double ionCharge)
  {
    System system = makeChargedChainSystem(false, charges, ionCharge);
    return gradientEnergy(system).ewald_fourier;
  };
  const double expected = fourier(chainCharges, 0.4) - fourier(chainCharges, 0.0) - fourier({}, 0.4);
  ASSERT_GT(std::abs(expected), 1e-3);

  System system = makeChargedChainSystem(false, chainCharges, 0.4);
  EnergyStatus status = system.computeMolecularPressure().first;
  EXPECT_NEAR(status.frameworkComponentEnergy(0, 0).CoulombicFourier.energy, expected, 1e-8);
}

// The bonded pair of two lab-fixed atoms is no part of any energy: its Fourier energy is subtracted with the stored
// structure factor, so an exclusion correction would be left over.
TEST(BAROSTAT_VIRIAL, fixed_framework_pairs_have_no_exclusion_energy)
{
  System system = makeChargedChainSystem(true, {0.3, -0.3, 0.0, 0.0, 0.0}, 0.0);
  EXPECT_DOUBLE_EQ(gradientEnergy(system).ewald_exclusion, 0.0);
  EXPECT_DOUBLE_EQ(system.computeTotalEnergies().ewald_exclusion, 0.0);
}

// The virial of a framework with lab-fixed and flexible atoms is the derivative of the energy it reports.
TEST(BAROSTAT_VIRIAL, mixed_framework_ewald_virial_matches_finite_difference)
{
  System system = makeChargedChainSystem(true, chainCharges, 0.4);
  const double delta = 1e-5;
  const double scalePlus = std::cbrt(1.0 + delta);
  const double scaleMinus = std::cbrt(1.0 - delta);
  System plus = scaledSystem(system, scalePlus);
  System minus = scaledSystem(system, scaleMinus);
  const double fd = -(gradientEnergy(plus).potentialEnergy() - gradientEnergy(minus).potentialEnergy()) /
                    (plus.simulationBox.volume - minus.simulationBox.volume);
  ASSERT_GT(std::abs(fd), 1e-6);

  const double3x3 virial = system.computeMolecularPressure().second;
  EXPECT_NEAR(virial.trace() / (3.0 * system.simulationBox.volume), fd, 5e-5 * std::abs(fd) + 1e-8);
}

// With the inter-molecular interactions omitted, the strain routine differentiates the framework-adsorbate energy
// of computeEwaldFourierEnergy, and neither reports nor differentiates the adsorbate-adsorbate terms.
TEST(BAROSTAT_VIRIAL, omit_inter_interactions_applies_to_the_ewald_strain_derivative)
{
  ForceField forceField = ForceField::makeZeoliteForceField(7.5, true, false, true);
  forceField.omitInterInteractions = true;
  forceField.automaticEwald = false;
  forceField.EwaldAlpha = 0.3;
  forceField.numberOfWaveVectors = int3(6, 6, 6);
  Framework framework = Framework::makeITQ29(forceField, int3(1, 1, 1));
  Component co2 = Component::makeCO2(forceField, 0, true);
  System system(forceField, std::nullopt, false, 300.0, 1.0e4, 1.0, {framework}, {co2}, {}, {2}, 5);
  std::span<Atom> molecules = system.spanOfMoleculeAtoms();
  molecules[0].position = double3(5.93355, 7.93355, 5.93355 + 1.149);
  molecules[1].position = double3(5.93355, 7.93355, 5.93355);
  molecules[2].position = double3(5.93355, 7.93355, 5.93355 - 1.149);
  molecules[3].position = double3(5.93355, 3.93355, 5.93355 + 1.149);
  molecules[4].position = double3(5.93355, 3.93355, 5.93355);
  molecules[5].position = double3(5.93355, 3.93355, 5.93355 - 1.149);
  system.precomputeTotalRigidEnergy();

  const auto ewaldEnergy = [](System& s)
  {
    const RunningEnergy energy = Interactions::computeEwaldFourierEnergy(
        s.eik_x, s.eik_y, s.eik_z, s.eik_xy, s.fixedFrameworkStoredEik, s.storedEik, s.forceField, s.simulationBox,
        s.components, s.numberOfMoleculesPerComponent, s.spanOfMoleculeAtoms(), s.netChargeFramework);
    return energy.ewald_fourier + energy.ewald_self + energy.ewald_exclusion;
  };

  std::vector<AtomDynamics> dynamics(system.spanOfMoleculeAtoms().size());
  auto [status, strain] = Interactions::computeEwaldFourierEnergyStrainDerivative(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
      system.forceField, system.simulationBox, system.framework, system.components, system.numberOfMoleculesPerComponent,
      system.spanOfMoleculeAtoms(), dynamics, system.netChargeFramework, system.netChargePerComponent);

  EXPECT_DOUBLE_EQ(status.componentEnergy(0, 0).CoulombicFourier.energy, 0.0);
  EXPECT_NEAR(status.frameworkComponentEnergy(0, 0).CoulombicFourier.energy, ewaldEnergy(system), 1e-8);

  const double delta = 1e-5;
  const double scalePlus = std::cbrt(1.0 + delta);
  const double scaleMinus = std::cbrt(1.0 - delta);
  System plus = scaledSystem(system, scalePlus);
  System minus = scaledSystem(system, scaleMinus);
  const double fd =
      (ewaldEnergy(plus) - ewaldEnergy(minus)) / (plus.simulationBox.volume - minus.simulationBox.volume);
  ASSERT_GT(std::abs(fd), 1e-6);
  EXPECT_NEAR(strain.trace() / (3.0 * system.simulationBox.volume), fd, 1e-4 * std::abs(fd) + 1e-8);

  // the force on an adsorbate atom comes from the framework only
  const double positionDelta = 1e-5;
  const double x = system.spanOfMoleculeAtoms()[0].position.x;
  system.spanOfMoleculeAtoms()[0].position.x = x + positionDelta;
  const double energyPlus = ewaldEnergy(system);
  system.spanOfMoleculeAtoms()[0].position.x = x - positionDelta;
  const double energyMinus = ewaldEnergy(system);
  system.spanOfMoleculeAtoms()[0].position.x = x;
  EXPECT_NEAR(dynamics[0].gradient.x, (energyPlus - energyMinus) / (2.0 * positionDelta), 1e-5);
}

// The Ewald strain derivative of a rigid framework differentiates the energy it reports: the framework has no energy
// with itself, so its own reciprocal-space term (which depends on alpha) must not enter the virial.
TEST(BAROSTAT_VIRIAL, rigid_framework_ewald_strain_matches_the_energy_derivative_for_every_alpha)
{
  for (const double alpha : {0.25, 0.3, 0.4})
  {
    ForceField forceField = ForceField::makeZeoliteForceField(7.5, true, false, true);
    forceField.automaticEwald = false;
    forceField.EwaldAlpha = alpha;
    forceField.numberOfWaveVectors = int3(10, 10, 10);
    Framework framework = Framework::makeITQ29(forceField, int3(1, 1, 1));
    Component co2 = Component::makeCO2(forceField, 0, true);
    System system(forceField, std::nullopt, false, 300.0, 1.0e4, 1.0, {framework}, {co2}, {}, {1}, 5);
    std::span<Atom> molecules = system.spanOfMoleculeAtoms();
    molecules[0].position = double3(5.93355, 7.93355, 5.93355 + 1.149);
    molecules[1].position = double3(5.93355, 7.93355, 5.93355);
    molecules[2].position = double3(5.93355, 7.93355, 5.93355 - 1.149);
    system.precomputeTotalRigidEnergy();

    const auto ewaldEnergy = [](System& s)
    {
      const RunningEnergy energy = Interactions::computeEwaldFourierEnergy(
          s.eik_x, s.eik_y, s.eik_z, s.eik_xy, s.fixedFrameworkStoredEik, s.storedEik, s.forceField, s.simulationBox,
          s.components, s.numberOfMoleculesPerComponent, s.spanOfMoleculeAtoms(), s.netChargeFramework);
      return energy.ewald_fourier + energy.ewald_self + energy.ewald_exclusion;
    };

    std::vector<AtomDynamics> dynamics(molecules.size());
    auto [status, strain] = Interactions::computeEwaldFourierEnergyStrainDerivative(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
        system.forceField, system.simulationBox, system.framework, system.components,
        system.numberOfMoleculesPerComponent, system.spanOfMoleculeAtoms(), dynamics, system.netChargeFramework,
        system.netChargePerComponent);

    const double delta = 1e-5;
    System plus = scaledSystem(system, std::cbrt(1.0 + delta));
    System minus = scaledSystem(system, std::cbrt(1.0 - delta));
    ASSERT_DOUBLE_EQ(plus.forceField.EwaldAlpha, alpha);
    const double fd =
        (ewaldEnergy(plus) - ewaldEnergy(minus)) / (plus.simulationBox.volume - minus.simulationBox.volume);
    ASSERT_GT(std::abs(fd), 1e-3);
    EXPECT_NEAR(strain.trace() / (3.0 * system.simulationBox.volume), fd, 2e-3 * std::abs(fd) + 1e-7) << "alpha " << alpha;
  }
}

// The tail correction of the virial holds the same pairs as the tail energy: framework-molecule and molecule-molecule,
// not framework-framework. An empty framework has no tail pressure, and the correction of a loaded one is the sum over
// those pairs.
TEST(BAROSTAT_VIRIAL, tail_pressure_excludes_framework_framework_pairs)
{
  const auto makeSystem = [](bool tailCorrections, std::size_t numberOfMolecules)
  {
    ForceField forceField = ForceField::makeZeoliteForceField(12.0, false, tailCorrections, false);
    Framework framework = Framework::makeITQ29(forceField, int3(2, 2, 2));
    Component methane = Component::makeMethane(forceField, 0);
    System system(forceField, std::nullopt, false, 300.0, 1.0e4, 1.0, {framework}, {methane}, {}, {numberOfMolecules}, 5);
    std::span<Atom> molecules = system.spanOfMoleculeAtoms();
    for (std::size_t i = 0; i < molecules.size(); ++i)
    {
      molecules[i].position = double3(5.93355 + 12.0 * static_cast<double>(i), 7.93355, 5.93355);
    }
    return system;
  };

  System empty = makeSystem(true, 0);
  ASSERT_NE(empty.forceField(0, 1).tailCorrectionPressure, 0.0);
  EXPECT_DOUBLE_EQ(empty.computeMolecularPressure().second.trace(), 0.0);

  System withTail = makeSystem(true, 2);
  System withoutTail = makeSystem(false, 2);

  const std::size_t numberOfTypes = withTail.forceField.numberOfPseudoAtoms;
  std::vector<double> frameworkCounts(numberOfTypes, 0.0);
  std::vector<double> moleculeCounts(numberOfTypes, 0.0);
  for (const Atom& atom : withTail.spanOfFrameworkAtoms()) frameworkCounts[atom.type] += atom.scalingVDW;
  for (const Atom& atom : withTail.spanOfMoleculeAtoms()) moleculeCounts[atom.type] += atom.scalingVDW;
  double expected = 0.0;
  for (std::size_t a = 0; a < numberOfTypes; ++a)
  {
    for (std::size_t b = 0; b < numberOfTypes; ++b)
    {
      expected += (moleculeCounts[a] * moleculeCounts[b] + 2.0 * frameworkCounts[a] * moleculeCounts[b]) *
                  withTail.forceField(a, b).tailCorrectionPressure;
    }
  }
  expected *= -2.0 * std::numbers::pi / (3.0 * withTail.simulationBox.volume);
  ASSERT_NE(expected, 0.0);

  const double difference =
      withTail.computeMolecularPressure().second.trace() - withoutTail.computeMolecularPressure().second.trace();
  EXPECT_NEAR(difference, 3.0 * expected, 1e-9 * std::abs(expected));
}

// The lab-fixed atoms of a mixed framework are a constant of the framework: the Ewald path leaves them out of its
// live atoms, so their self energy is no part of the system energy. The real-space charge methods summed the self
// energy of every framework atom, in the total energy, the gradient and the Hessian.
TEST(BAROSTAT_VIRIAL, real_space_self_energy_counts_the_mobile_framework_atoms_only)
{
  System system = makeChargedChainSystem(true, {0.3, -0.3, 0.2, -0.1, -0.1}, 0.4);
  ASSERT_TRUE(system.framework->isMixed());
  const std::size_t fixedAtoms = system.framework->numberOfFixedAtoms();
  ASSERT_EQ(fixedAtoms, 2uz);
  system.forceField.chargeMethod = ForceField::ChargeMethod::Wolf;

  const double prefactor = Units::CoulombicConversionFactor * Potentials::coulombSelfEnergyPrefactor(system.forceField);
  double expected = 0.0;
  for (const Atom& atom : system.spanOfFrameworkAtoms().subspan(fixedAtoms))
  {
    expected += prefactor * atom.scalingCoulomb * atom.charge * atom.scalingCoulomb * atom.charge;
  }
  for (const Atom& atom : system.spanOfMoleculeAtoms())
  {
    expected += prefactor * atom.scalingCoulomb * atom.charge * atom.scalingCoulomb * atom.charge;
  }
  ASSERT_NE(expected, 0.0);

  EXPECT_NEAR(system.computeTotalEnergies().ewald_self, expected, 1e-8 * std::abs(expected));
  EXPECT_NEAR(gradientEnergy(system).ewald_self, expected, 1e-8 * std::abs(expected));
}

// A caller without framework dynamics (dcTST hands an empty span) must still get the complete Ewald energies of a
// framework with mobile atoms. The gradient and the strain-derivative routines used to drop every mobile framework
// charge from the Fourier sum, the self term and the exclusion term in that case; only the framework forces may go.
TEST(BAROSTAT_VIRIAL, ewald_energies_do_not_depend_on_the_framework_dynamics_being_passed)
{
  System system = makeChargedChainSystem(false, {0.3, -0.3, 0.2, -0.1, -0.1}, 0.4);
  ASSERT_TRUE(system.framework->hasMobileAtoms());

  auto gradientEnergies = [&](std::span<AtomDynamics> frameworkDynamics)
  {
    for (AtomDynamics& dynamics : system.spanOfMoleculeDynamics()) dynamics.gradient = double3(0.0, 0.0, 0.0);
    return Interactions::computeEwaldFourierGradient(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.trialEik, system.fixedFrameworkStoredEik,
        system.forceField, system.simulationBox, system.components, system.numberOfMoleculesPerComponent,
        system.spanOfMoleculeAtoms(), system.spanOfMoleculeDynamics(), system.netChargeFramework, system.framework,
        system.spanOfFrameworkAtoms(), frameworkDynamics);
  };
  const RunningEnergy with = gradientEnergies(system.spanOfFrameworkDynamics());
  const RunningEnergy without = gradientEnergies({});
  ASSERT_NE(with.ewald_fourier, 0.0);
  EXPECT_NEAR(without.ewald_fourier, with.ewald_fourier, 1e-10 * std::abs(with.ewald_fourier));
  EXPECT_NEAR(without.ewald_self, with.ewald_self, 1e-10 * std::abs(with.ewald_self));
  EXPECT_NEAR(without.ewald_exclusion, with.ewald_exclusion, 1e-10 * std::max(1.0, std::abs(with.ewald_exclusion)));

  auto strain = [&](std::span<AtomDynamics> frameworkDynamics)
  {
    return Interactions::computeEwaldFourierEnergyStrainDerivative(
               system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik,
               system.storedEik, system.forceField, system.simulationBox, system.framework, system.components,
               system.numberOfMoleculesPerComponent, system.spanOfMoleculeAtoms(), system.spanOfMoleculeDynamics(),
               system.netChargeFramework, system.netChargePerComponent, system.spanOfFrameworkAtoms(),
               frameworkDynamics)
        .second;
  };
  const double3x3 strainWith = strain(system.spanOfFrameworkDynamics());
  const double3x3 strainWithout = strain({});
  auto near = [](double a, double b) { EXPECT_NEAR(a, b, 1e-10 * std::max(1.0, std::abs(b))); };
  near(strainWithout.ax, strainWith.ax);
  near(strainWithout.ay, strainWith.ay);
  near(strainWithout.az, strainWith.az);
  near(strainWithout.bx, strainWith.bx);
  near(strainWithout.by, strainWith.by);
  near(strainWithout.bz, strainWith.bz);
  near(strainWithout.cx, strainWith.cx);
  near(strainWithout.cy, strainWith.cy);
  near(strainWithout.cz, strainWith.cz);
  EXPECT_NE(strainWith.ax, 0.0);
}
