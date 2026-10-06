#include <gtest/gtest.h>

import std;

import int3;
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
import energy_status;
import running_energy;
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
