#include <gtest/gtest.h>

import std;

import double3;
import double3x3;
import int3;
import atom;
import molecule;
import simd_quatd;
import forcefield;
import component;
import framework;
import connectivity_table;
import intra_molecular_potentials;
import bond_potential;
import pseudo_atom;
import system;
import simulationbox;
import normal_modes;
import phonon_force_constants;
import phonon_dynamical_matrix;
import phonon_kpath;

namespace
{
Atom fractionalAtom(double3 position, std::size_t type)
{
  return Atom(position, 0.0, 1.0, 0, static_cast<std::uint16_t>(type), 0, 0, true);
}

void setRigidMoleculePositions(System& system, std::size_t moleculeIndex)
{
  Molecule& molecule = system.moleculeData[moleculeIndex];
  const Component& component = system.components[molecule.componentId];
  const double3x3 rotation = double3x3::buildRotationMatrixInverse(molecule.orientation);
  for (std::size_t localAtom = 0; localAtom < molecule.numberOfAtoms; ++localAtom)
  {
    system.spanOfMoleculeAtoms()[molecule.atomIndex + localAtom].position =
        molecule.centerOfMassPosition + rotation * component.atoms[localAtom].position;
  }
}

// Flexible one-dimensional framework chain (a single harmonic bond wrapped across the x boundary, so
// the force constants carry a nonzero-image block) plus a flexible single-site adsorbate. Neutral and
// non-polarizable, so the real-space force constants are complete.
System makeFlexibleChainSystem()
{
  ForceField forceField({{"C", false, 12.0, 0.0, 0.0, 6, false}, {"M", false, 16.0, 0.0, 0.0, 6, false}},
                        {{120.0, 3.4}, {100.0, 3.5}}, ForceField::MixingRule::Lorentz_Berthelot, 6.0, 6.0, 6.0, false,
                        false, false);
  const SimulationBox box(12.0, 12.0, 12.0);
  const std::size_t carbon = *forceField.findPseudoAtom("C");
  const std::size_t methaneLike = *forceField.findPseudoAtom("M");

  std::vector<Atom> frameworkAtoms{fractionalAtom({0.05, 0.5, 0.5}, carbon), fractionalAtom({0.95, 0.5, 0.5}, carbon)};
  Framework framework(forceField, "chain", box, 1, frameworkAtoms, frameworkAtoms, {1, 1, 1});
  framework.rigid = false;
  framework.intraMolecularPotentials.bonds.emplace_back(std::array<std::size_t, 2>{0, 1}, BondType::Harmonic,
                                                        std::vector<double>{1500.0, 1.2});
  framework.intraMolecularImageShifts.bonds.push_back({int3{0, 0, 0}, int3{-1, 0, 0}});

  Component molecule(forceField, "single-site", 16.0, 1.0e6, 0.011,
                     {Atom({0.0, 0.0, 0.0}, 0.0, 1.0, 0, static_cast<std::uint16_t>(methaneLike), 0, false, false)},
                     ConnectivityTable(1), Potentials::IntraMolecularPotentials{}, 0, 0);
  molecule.rigid = false;

  System system(forceField, box, false, 300.0, 1e4, 1.0, {framework}, {molecule}, {}, {1}, 5);
  system.spanOfMoleculeAtoms()[0].position = double3(0.6, 7.3, 6.0);  // ~1.3 A from framework atom 0
  return system;
}

// A charge-neutral gas of flexible single-atom ions (one +0.8 cation, two -0.4 anions) with Ewald
// summation enabled. Single-atom molecules have no intramolecular exclusions, so the reciprocal-space
// Ewald contribution to the dynamical matrix is complete: real-space force constants (van der Waals +
// erfc real Coulomb) plus the reciprocal Fourier term must reconstruct the full Gamma-point Hessian.
System makeChargedIonSystem()
{
  ForceField forceField({{"P", false, 15.0, 0.8, 0.0, 8, false}, {"N", false, 14.0, -0.4, 0.0, 8, false}},
                        {{60.0, 3.0}, {40.0, 3.2}}, ForceField::MixingRule::Lorentz_Berthelot, 11.0, 11.0, 11.0, true,
                        false, true);
  const SimulationBox box(25.0, 25.0, 25.0);

  Component cation(forceField, "cation", 100.0, 1e6, 0.2,
                   {Atom({0.0, 0.0, 0.0}, 0.8, 1.0, 0, 0, 0, false, false)}, ConnectivityTable(1),
                   Potentials::IntraMolecularPotentials{}, 0, 0);
  cation.rigid = false;
  Component anion(forceField, "anion", 100.0, 1e6, 0.2,
                  {Atom({0.0, 0.0, 0.0}, -0.4, 1.0, 0, 1, 0, false, false)}, ConnectivityTable(1),
                  Potentials::IntraMolecularPotentials{}, 0, 0);
  anion.rigid = false;

  System system(forceField, box, false, 300.0, 1e4, 1.0, {}, {cation, anion}, {}, {1, 2}, 5);
  std::span<Atom> atoms = system.spanOfMoleculeAtoms();
  atoms[0].position = double3(9.4, 11.2, 12.7);   // cation
  atoms[1].position = double3(13.1, 10.3, 9.8);   // anion
  atoms[2].position = double3(11.6, 14.0, 13.2);  // anion
  return system;
}

// Two flexible bent three-site molecules (charges -0.4, +0.8, -0.4, each net neutral) with Ewald.
// Intramolecular Coulomb is fully excluded, so the reciprocal Fourier term is corrected by the
// intramolecular exclusion (folded into the real-space force constants). Mirrors the charged-pair
// system of the Ewald Hessian finite-difference test.
System makeChargedMoleculePairSystem()
{
  ForceField forceField({{"P", false, 15.0, 0.8, 0.0, 8, false}, {"N", false, 14.0, -0.4, 0.0, 8, false}},
                        {{60.0, 3.0}, {40.0, 3.2}}, ForceField::MixingRule::Lorentz_Berthelot, 11.0, 11.0, 11.0, true,
                        false, true);
  const SimulationBox box(25.0, 25.0, 25.0);

  ConnectivityTable connectivityTable(3);
  connectivityTable[0, 1] = true;
  connectivityTable[1, 0] = true;
  connectivityTable[1, 2] = true;
  connectivityTable[2, 1] = true;

  Component component(forceField, "flexibleIon", 100.0, 1e6, 0.2,
                      {Atom({-0.9, 0.4, 0.1}, -0.4, 1.0, 0, 1, 0, false, false),
                       Atom({0.0, 0.0, 0.0}, 0.8, 1.0, 0, 0, 0, false, false),
                       Atom({1.0, 0.3, -0.2}, -0.4, 1.0, 0, 1, 0, false, false)},
                      connectivityTable, Potentials::IntraMolecularPotentials{}, 5, 21);
  component.rigid = false;

  System system(forceField, box, false, 300.0, 1e4, 1.0, {}, {component}, {}, {2}, 5);
  const std::array<double3, 6> positions = {double3(9.1, 10.4, 10.1),  double3(10.0, 10.0, 10.0),
                                            double3(11.0, 10.3, 9.8),  double3(12.6, 11.9, 10.7),
                                            double3(13.5, 11.5, 10.6), double3(14.4, 11.8, 11.3)};
  std::span<Atom> atoms = system.spanOfMoleculeAtoms();
  for (std::size_t i = 0; i < positions.size(); ++i) atoms[i].position = positions[i];
  return system;
}

Atom chargedFrameworkAtom(double3 position, double charge, std::size_t type)
{
  return Atom(position, charge, 1.0, 0, static_cast<std::uint16_t>(type), 0, 0, true);
}

// A charge-neutral flexible framework (two bonded atoms +0.4/-0.4 across the x boundary, with a harmonic
// bond) plus one charged single-atom adsorbate. The bonded framework pair is a Coulomb exclusion, so the
// framework branch of the exclusion correction must be folded into the force constants; the reciprocal
// Fourier term includes the flexible framework atoms as sites.
System makeChargedFrameworkSystem()
{
  ForceField forceField({{"P", false, 15.0, 0.8, 0.0, 8, false}, {"N", false, 14.0, -0.4, 0.0, 8, false}},
                        {{60.0, 3.0}, {40.0, 3.2}}, ForceField::MixingRule::Lorentz_Berthelot, 11.0, 11.0, 11.0, true,
                        false, true);
  const SimulationBox box(24.0, 24.0, 24.0);
  const std::size_t typeP = *forceField.findPseudoAtom("P");
  const std::size_t typeN = *forceField.findPseudoAtom("N");

  std::vector<Atom> frameworkAtoms{chargedFrameworkAtom({0.05, 0.5, 0.5}, 0.4, typeP),
                                   chargedFrameworkAtom({0.95, 0.5, 0.5}, -0.4, typeN)};
  Framework framework(forceField, "charged-chain", box, 1, frameworkAtoms, frameworkAtoms, {1, 1, 1});
  framework.rigid = false;
  framework.connectivityTable = ConnectivityTable(2);
  framework.connectivityTable[0, 1] = true;
  framework.connectivityTable[1, 0] = true;
  framework.intraMolecularPotentials.bonds.emplace_back(std::array<std::size_t, 2>{0, 1}, BondType::Harmonic,
                                                        std::vector<double>{1500.0, 2.4});
  framework.intraMolecularImageShifts.bonds.push_back({int3{0, 0, 0}, int3{-1, 0, 0}});

  Component molecule(forceField, "ion", 30.0, 1.0e6, 0.011,
                     {Atom({0.0, 0.0, 0.0}, -0.4, 1.0, 0, static_cast<std::uint16_t>(typeN), 0, false, false)},
                     ConnectivityTable(1), Potentials::IntraMolecularPotentials{}, 0, 0);
  molecule.rigid = false;

  System system(forceField, box, false, 300.0, 1e4, 1.0, {framework}, {molecule}, {}, {1}, 5);
  system.spanOfMoleculeAtoms()[0].position = double3(6.0, 14.0, 12.0);
  return system;
}

// A rigid charged framework (two fixed +0.4/-0.4 atoms) plus two flexible charged single-atom ions.
// The framework atoms carry no degrees of freedom, but their fixed charges shift each ion's self term
// through the static reciprocal structure factor.
System makeRigidChargedFrameworkSystem()
{
  ForceField forceField({{"P", false, 15.0, 0.8, 0.0, 8, false}, {"N", false, 14.0, -0.4, 0.0, 8, false}},
                        {{60.0, 3.0}, {40.0, 3.2}}, ForceField::MixingRule::Lorentz_Berthelot, 11.0, 11.0, 11.0, true,
                        false, true);
  const SimulationBox box(24.0, 24.0, 24.0);
  const std::size_t typeP = *forceField.findPseudoAtom("P");
  const std::size_t typeN = *forceField.findPseudoAtom("N");

  std::vector<Atom> frameworkAtoms{chargedFrameworkAtom({0.25, 0.5, 0.5}, 0.4, typeP),
                                   chargedFrameworkAtom({0.75, 0.5, 0.5}, -0.4, typeN)};
  Framework framework(forceField, "rigid-charged", box, 1, frameworkAtoms, frameworkAtoms, {1, 1, 1});
  framework.rigid = true;

  Component cation(forceField, "cation", 30.0, 1.0e6, 0.011,
                   {Atom({0.0, 0.0, 0.0}, 0.4, 1.0, 0, static_cast<std::uint16_t>(typeP), 0, false, false)},
                   ConnectivityTable(1), Potentials::IntraMolecularPotentials{}, 0, 0);
  cation.rigid = false;
  Component anion(forceField, "anion", 30.0, 1.0e6, 0.011,
                  {Atom({0.0, 0.0, 0.0}, -0.4, 1.0, 0, static_cast<std::uint16_t>(typeN), 0, false, false)},
                  ConnectivityTable(1), Potentials::IntraMolecularPotentials{}, 0, 0);
  anion.rigid = false;

  System system(forceField, box, false, 300.0, 1e4, 1.0, {framework}, {cation, anion}, {}, {1, 1}, 5);
  std::span<Atom> atoms = system.spanOfMoleculeAtoms();
  atoms[0].position = double3(6.0, 14.0, 12.0);
  atoms[1].position = double3(17.0, 9.5, 13.5);
  return system;
}

// Three rigid two-site "dumbbells" (linear, uncharged Lennard-Jones) in a periodic chain along x. Each
// molecule carries six generalized DOFs (center of mass + orientation); the linear geometry gives one
// zero-inertia rotation per molecule, so the phonon route must project the Cartesian force constants onto
// the center-of-mass/orientation coordinates and mass-weight with the inertia tensor. The nearest-neighbor
// coupling graph closes a loop with nonzero winding (0->1 and 1->2 in the home cell, 2->0 across the x
// boundary), so a single periodic image cannot be gauged away and the bands genuinely disperse.
System makeRigidMoleculeChainSystem()
{
  ForceField forceField({{"O", false, 16.0, 0.0, 0.0, 8, false}}, {{80.0, 3.0}},
                        ForceField::MixingRule::Lorentz_Berthelot, 4.0, 4.0, 4.0, false, false, false);
  const SimulationBox box(9.0, 9.0, 9.0);
  const std::size_t type = *forceField.findPseudoAtom("O");

  Component dumbbell(forceField, "dumbbell", 32.0, 1.0e6, 0.02,
                     {Atom({-0.6, 0.0, 0.0}, 0.0, 1.0, 0, static_cast<std::uint16_t>(type), 0, false, false),
                      Atom({0.6, 0.0, 0.0}, 0.0, 1.0, 0, static_cast<std::uint16_t>(type), 0, false, false)},
                     ConnectivityTable(2), Potentials::IntraMolecularPotentials{}, 0, 0);
  dumbbell.rigid = true;

  System system(forceField, box, false, 300.0, 1e4, 1.0, {}, {dumbbell}, {}, {3}, 5);
  // The molecules are inserted with random center-of-mass positions and orientations. Pin them to a
  // deterministic chain along x whose center-of-mass separations (~3 A) are all inside the 4 A cutoff,
  // so every nearest-neighbor pair interacts regardless of the internal body-frame orientation. The
  // loop closes across the x boundary in a 9 A box (0->1 and 1->2 in the home cell, 2->0 wrapping), so
  // the coupling graph encloses a nonzero lattice vector and the bands genuinely disperse.
  const std::array<double3, 3> centersOfMass = {double3(1.5, 4.5, 4.5), double3(4.5, 5.3, 4.0),
                                                double3(7.5, 3.8, 5.0)};
  for (std::size_t moleculeIndex = 0; moleculeIndex < 3; ++moleculeIndex)
  {
    system.moleculeData[moleculeIndex].centerOfMassPosition = centersOfMass[moleculeIndex];
    system.moleculeData[moleculeIndex].orientation = simd_quatd(0.0, 0.0, 0.0, 1.0);
    setRigidMoleculePositions(system, moleculeIndex);
  }
  return system;
}

// Three rigid, charge-neutral two-site "dipoles" (+0.5 / -0.5, linear Lennard-Jones) in a periodic chain
// along x with Ewald summation enabled. Each molecule keeps six generalized DOFs, the intramolecular
// charge pair is a rigid constant (fully excluded), and the reciprocal-space Ewald term must be projected
// onto the center-of-mass/orientation coordinates like the short-ranged part. The nearest-neighbor loop
// wraps across the x boundary of the 9 A box, so the bands genuinely disperse.
System makeChargedRigidMoleculeChainSystem()
{
  ForceField forceField({{"P", false, 16.0, 0.5, 0.0, 8, false}, {"N", false, 16.0, -0.5, 0.0, 8, false}},
                        {{80.0, 3.0}, {80.0, 3.0}}, ForceField::MixingRule::Lorentz_Berthelot, 4.0, 4.0, 4.0, false,
                        false, true);
  const SimulationBox box(9.0, 9.0, 9.0);
  const std::size_t typeP = *forceField.findPseudoAtom("P");
  const std::size_t typeN = *forceField.findPseudoAtom("N");

  Component dipole(forceField, "dipole", 32.0, 1.0e6, 0.02,
                   {Atom({-0.6, 0.0, 0.0}, 0.5, 1.0, 0, static_cast<std::uint16_t>(typeP), 0, false, false),
                    Atom({0.6, 0.0, 0.0}, -0.5, 1.0, 0, static_cast<std::uint16_t>(typeN), 0, false, false)},
                   ConnectivityTable(2), Potentials::IntraMolecularPotentials{}, 0, 0);
  dipole.rigid = true;

  System system(forceField, box, false, 300.0, 1e4, 1.0, {}, {dipole}, {}, {3}, 5);
  const std::array<double3, 3> centersOfMass = {double3(1.5, 4.5, 4.5), double3(4.5, 5.3, 4.0),
                                                double3(7.5, 3.8, 5.0)};
  for (std::size_t moleculeIndex = 0; moleculeIndex < 3; ++moleculeIndex)
  {
    system.moleculeData[moleculeIndex].centerOfMassPosition = centersOfMass[moleculeIndex];
    system.moleculeData[moleculeIndex].orientation = simd_quatd(0.0, 0.0, 0.0, 1.0);
    setRigidMoleculePositions(system, moleculeIndex);
  }
  return system;
}

// Pin a semi-flexible molecule to its component reference geometry translated by `offset` (identity
// orientation). The rigid-group geometry is then exactly rigid, consistent with the rigid-body
// regeneration used by the group-aware derivative machinery.
void setSemiFlexibleMoleculePositions(System& system, std::size_t moleculeIndex, double3 offset)
{
  Molecule& molecule = system.moleculeData[moleculeIndex];
  const Component& component = system.components[molecule.componentId];
  for (std::size_t localAtom = 0; localAtom < molecule.numberOfAtoms; ++localAtom)
  {
    system.spanOfMoleculeAtoms()[molecule.atomIndex + localAtom].position =
        offset + component.atoms[localAtom].position;
  }
  molecule.centerOfMassPosition = offset;
  molecule.orientation = simd_quatd(0.0, 0.0, 0.0, 1.0);
}

// Turn a four-site flexible component into a semi-flexible one: a rigid bent three-site core
// (nonlinear, so all three rotations carry inertia) plus one flexible tail atom attached through the
// harmonic junction bond 2-3. Each molecule carries 9 generalized DOFs (6 rigid-body + 3 Cartesian).
void makeComponentSemiFlexible(Component& component)
{
  component.groups = {MoleculeGroup(true, {0, 1, 2}), MoleculeGroup(false, {3})};
  component.atomGroupIds = {0, 0, 0, 1};
  component.computeGroupRigidProperties();
}

// Three semi-flexible molecules (rigid bent three-site core + flexible tail atom, uncharged
// Lennard-Jones) in a periodic chain along x. The nearest-neighbor coupling loop wraps the x boundary
// of the 9 A box, so a single periodic image cannot be gauged away and the bands genuinely disperse.
System makeSemiFlexibleMoleculeChainSystem()
{
  ForceField forceField({{"O", false, 16.0, 0.0, 0.0, 8, false}}, {{80.0, 3.0}},
                        ForceField::MixingRule::Lorentz_Berthelot, 4.0, 4.0, 4.0, false, false, false);
  const SimulationBox box(9.0, 9.0, 9.0);
  const std::uint16_t type = static_cast<std::uint16_t>(*forceField.findPseudoAtom("O"));

  ConnectivityTable connectivityTable(4);
  connectivityTable[0, 1] = true;
  connectivityTable[1, 0] = true;
  connectivityTable[1, 2] = true;
  connectivityTable[2, 1] = true;
  connectivityTable[2, 3] = true;
  connectivityTable[3, 2] = true;

  Potentials::IntraMolecularPotentials potentials;
  potentials.bonds.emplace_back(std::array<std::size_t, 2>{2, 3}, BondType::Harmonic,
                                std::vector<double>{1500.0, 0.85});

  Component semiFlexible(forceField, "semi-flexible", 64.0, 1.0e6, 0.02,
                         {Atom({-0.6, 0.25, 0.0}, 0.0, 1.0, 0, type, 0, false, false),
                          Atom({0.0, -0.35, 0.0}, 0.0, 1.0, 0, type, 0, false, false),
                          Atom({0.6, 0.25, 0.0}, 0.0, 1.0, 0, type, 0, false, false),
                          Atom({1.2, -0.25, 0.35}, 0.0, 1.0, 0, type, 0, false, false)},
                         connectivityTable, potentials, 0, 0);
  makeComponentSemiFlexible(semiFlexible);

  System system(forceField, box, false, 300.0, 1e4, 1.0, {}, {semiFlexible}, {}, {3}, 5);
  const std::array<double3, 3> centersOfMass = {double3(1.5, 4.5, 4.5), double3(4.5, 5.3, 4.0),
                                                double3(7.5, 3.8, 5.0)};
  for (std::size_t moleculeIndex = 0; moleculeIndex < 3; ++moleculeIndex)
  {
    setSemiFlexibleMoleculePositions(system, moleculeIndex, centersOfMass[moleculeIndex]);
  }
  return system;
}

// Charged variant of the semi-flexible chain: the rigid core carries -0.3/+0.6/-0.3 (net neutral, the
// tail is uncharged) with Ewald summation enabled. The intramolecular exclusion pairs inside the rigid
// core are constants of the motion; their Cartesian force-constant blocks must cancel against the
// gradient-curvature term after projection onto the group degrees of freedom.
System makeChargedSemiFlexibleMoleculeChainSystem()
{
  ForceField forceField({{"N", false, 16.0, -0.3, 0.0, 8, false},
                         {"P", false, 16.0, 0.6, 0.0, 8, false},
                         {"T", false, 16.0, 0.0, 0.0, 8, false}},
                        {{80.0, 3.0}, {80.0, 3.0}, {80.0, 3.0}}, ForceField::MixingRule::Lorentz_Berthelot, 4.0, 4.0,
                        4.0, false, false, true);
  const SimulationBox box(9.0, 9.0, 9.0);
  const std::uint16_t typeN = static_cast<std::uint16_t>(*forceField.findPseudoAtom("N"));
  const std::uint16_t typeP = static_cast<std::uint16_t>(*forceField.findPseudoAtom("P"));
  const std::uint16_t typeT = static_cast<std::uint16_t>(*forceField.findPseudoAtom("T"));

  ConnectivityTable connectivityTable(4);
  connectivityTable[0, 1] = true;
  connectivityTable[1, 0] = true;
  connectivityTable[1, 2] = true;
  connectivityTable[2, 1] = true;
  connectivityTable[2, 3] = true;
  connectivityTable[3, 2] = true;

  Potentials::IntraMolecularPotentials potentials;
  potentials.bonds.emplace_back(std::array<std::size_t, 2>{2, 3}, BondType::Harmonic,
                                std::vector<double>{1500.0, 0.85});

  Component semiFlexible(forceField, "semi-flexible-charged", 64.0, 1.0e6, 0.02,
                         {Atom({-0.6, 0.25, 0.0}, -0.3, 1.0, 0, typeN, 0, false, false),
                          Atom({0.0, -0.35, 0.0}, 0.6, 1.0, 0, typeP, 0, false, false),
                          Atom({0.6, 0.25, 0.0}, -0.3, 1.0, 0, typeN, 0, false, false),
                          Atom({1.2, -0.25, 0.35}, 0.0, 1.0, 0, typeT, 0, false, false)},
                         connectivityTable, potentials, 0, 0);
  makeComponentSemiFlexible(semiFlexible);

  System system(forceField, box, false, 300.0, 1e4, 1.0, {}, {semiFlexible}, {}, {3}, 5);
  const std::array<double3, 3> centersOfMass = {double3(1.5, 4.5, 4.5), double3(4.5, 5.3, 4.0),
                                                double3(7.5, 3.8, 5.0)};
  for (std::size_t moleculeIndex = 0; moleculeIndex < 3; ++moleculeIndex)
  {
    setSemiFlexibleMoleculePositions(system, moleculeIndex, centersOfMass[moleculeIndex]);
  }
  return system;
}
}  // namespace

TEST(phonon_dispersion, gamma_point_matches_normal_modes)
{
  System system = makeFlexibleChainSystem();

  const NormalModesResult modes = computeNormalModes(system);
  ASSERT_EQ(modes.numberOfModes, 9u);

  const std::vector<PhononModes> dispersion = computePhononDispersion(system, std::array<double3, 1>{double3(0, 0, 0)});
  ASSERT_EQ(dispersion.size(), 1u);
  ASSERT_EQ(dispersion[0].eigenvalues.size(), 9u);

  double scale = 0.0;
  for (const double value : modes.eigenvalues) scale = std::max(scale, std::abs(value));
  ASSERT_GT(scale, 0.0);

  // Both routines return ascending eigenvalues; at k = 0 the dynamical matrix is the mass-weighted
  // folded Hessian, so the squared frequencies must agree with the independent normal-mode analysis.
  for (std::size_t mode = 0; mode < modes.numberOfModes; ++mode)
  {
    EXPECT_NEAR(dispersion[0].eigenvalues[mode], modes.eigenvalues[mode], 1e-7 * scale + 1e-10)
        << "mode=" << mode;
  }
}

TEST(phonon_dispersion, gamma_point_has_three_acoustic_zero_modes)
{
  System system = makeFlexibleChainSystem();
  const ForceConstants forceConstants = computeRealSpaceForceConstants(system);
  const std::vector<double> inverseSqrtMass = phononInverseSqrtMasses(system);

  const PhononModes gamma = computePhononModes(forceConstants, inverseSqrtMass, double3(0.0, 0.0, 0.0));
  ASSERT_EQ(gamma.eigenvalues.size(), 9u);

  double scale = 0.0;
  for (const double value : gamma.eigenvalues) scale = std::max(scale, std::abs(value));
  ASSERT_GT(scale, 0.0);

  // Translational invariance of the force constants guarantees three zero-frequency acoustic modes
  // (uniform-translation eigenvectors), independent of the geometry. They need not be the lowest
  // eigenvalues away from a minimum, so count how many eigenvalues vanish.
  const auto zeroModes =
      std::ranges::count_if(gamma.eigenvalues, [&](double value) { return std::abs(value) < 1e-7 * scale; });
  EXPECT_EQ(zeroModes, 3);
}

TEST(phonon_dispersion, dynamical_matrix_is_hermitian_and_time_reversal_symmetric)
{
  System system = makeFlexibleChainSystem();
  const ForceConstants forceConstants = computeRealSpaceForceConstants(system);
  const std::vector<double> inverseSqrtMass = phononInverseSqrtMasses(system);
  const std::size_t dimension = forceConstants.dimension();

  const double3 k(0.3, 0.0, 0.0);
  const std::vector<std::complex<double>> matrix = computeDynamicalMatrix(forceConstants, inverseSqrtMass, k);
  ASSERT_EQ(matrix.size(), dimension * dimension);

  double maxHermitianError = 0.0;
  double maxMagnitude = 0.0;
  for (std::size_t row = 0; row < dimension; ++row)
  {
    for (std::size_t column = 0; column < dimension; ++column)
    {
      const std::complex<double> upper = matrix[row * dimension + column];
      const std::complex<double> lower = matrix[column * dimension + row];
      maxHermitianError = std::max(maxHermitianError, std::abs(upper - std::conj(lower)));
      maxMagnitude = std::max(maxMagnitude, std::abs(upper));
    }
  }
  ASSERT_GT(maxMagnitude, 0.0);
  EXPECT_LT(maxHermitianError, 1e-10 * maxMagnitude);

  // The nonzero image block makes the k-dependence real: at least one squared frequency must move.
  const PhononModes atK = computePhononModes(forceConstants, inverseSqrtMass, k);
  const PhononModes atGamma = computePhononModes(forceConstants, inverseSqrtMass, double3(0.0, 0.0, 0.0));
  double maxShift = 0.0;
  for (std::size_t mode = 0; mode < atK.eigenvalues.size(); ++mode)
  {
    maxShift = std::max(maxShift, std::abs(atK.eigenvalues[mode] - atGamma.eigenvalues[mode]));
  }
  EXPECT_GT(maxShift, 0.0);

  // Time-reversal symmetry: omega^2(k) = omega^2(-k).
  const PhononModes atMinusK = computePhononModes(forceConstants, inverseSqrtMass, double3(-0.3, 0.0, 0.0));
  double scale = 0.0;
  for (const double value : atK.eigenvalues) scale = std::max(scale, std::abs(value));
  for (std::size_t mode = 0; mode < atK.eigenvalues.size(); ++mode)
  {
    EXPECT_NEAR(atK.eigenvalues[mode], atMinusK.eigenvalues[mode], 1e-9 * scale + 1e-12) << "mode=" << mode;
  }
}

TEST(phonon_dispersion, ewald_reciprocal_matrix_is_hermitian_and_nonzero)
{
  System system = makeChargedIonSystem();
  const std::vector<double> inverseSqrtMass = phononInverseSqrtMasses(system);
  const std::size_t dimension = 3 * inverseSqrtMass.size();

  const double3 k(0.35, -0.1, 0.2);
  const std::vector<std::complex<double>> reciprocal =
      computeEwaldReciprocalDynamicalMatrix(system, inverseSqrtMass, k);
  ASSERT_EQ(reciprocal.size(), dimension * dimension);

  double maxHermitianError = 0.0;
  double maxMagnitude = 0.0;
  for (std::size_t row = 0; row < dimension; ++row)
  {
    for (std::size_t column = 0; column < dimension; ++column)
    {
      const std::complex<double> upper = reciprocal[row * dimension + column];
      const std::complex<double> lower = reciprocal[column * dimension + row];
      maxHermitianError = std::max(maxHermitianError, std::abs(upper - std::conj(lower)));
      maxMagnitude = std::max(maxMagnitude, std::abs(upper));
    }
  }
  ASSERT_GT(maxMagnitude, 0.0);  // the reciprocal term must actually contribute
  EXPECT_LT(maxHermitianError, 1e-10 * maxMagnitude);

  // At k = 0 the reciprocal contribution is real (its imaginary part cancels between +/- G).
  const std::vector<std::complex<double>> atGamma =
      computeEwaldReciprocalDynamicalMatrix(system, inverseSqrtMass, double3(0.0, 0.0, 0.0));
  double maxImaginary = 0.0;
  for (const std::complex<double>& value : atGamma) maxImaginary = std::max(maxImaginary, std::abs(value.imag()));
  EXPECT_LT(maxImaginary, 1e-10 * maxMagnitude);
}

TEST(phonon_dispersion, gamma_point_with_ewald_matches_normal_modes)
{
  System system = makeChargedIonSystem();

  // computeNormalModes mass-weights the full generalized Hessian from evaluateDerivatives, which
  // already contains the Ewald Fourier + self + exclusion terms. For single-atom ions the exclusion
  // and (position-independent) self terms vanish, so the phonon route (real-space force constants +
  // reciprocal Fourier) must reproduce the same Gamma-point squared frequencies.
  const NormalModesResult modes = computeNormalModes(system);
  ASSERT_EQ(modes.numberOfModes, 9u);

  const std::vector<PhononModes> dispersion = computePhononDispersion(system, std::array<double3, 1>{double3(0, 0, 0)});
  ASSERT_EQ(dispersion.size(), 1u);
  ASSERT_EQ(dispersion[0].eigenvalues.size(), 9u);

  double scale = 0.0;
  for (const double value : modes.eigenvalues) scale = std::max(scale, std::abs(value));
  ASSERT_GT(scale, 0.0);

  for (std::size_t mode = 0; mode < modes.numberOfModes; ++mode)
  {
    EXPECT_NEAR(dispersion[0].eigenvalues[mode], modes.eigenvalues[mode], 1e-6 * scale + 1e-10) << "mode=" << mode;
  }
}

TEST(phonon_dispersion, ewald_reciprocal_shifts_gamma_frequencies)
{
  System system = makeChargedIonSystem();
  const ForceConstants forceConstants = computeRealSpaceForceConstants(system);
  const std::vector<double> inverseSqrtMass = phononInverseSqrtMasses(system);

  // Real-space-only versus real-space + reciprocal Fourier at Gamma; the long-ranged Coulomb term
  // must change the spectrum.
  const PhononModes realSpaceOnly = computePhononModes(forceConstants, inverseSqrtMass, double3(0.0, 0.0, 0.0));
  const PhononModes withReciprocal =
      computePhononModes(system, forceConstants, inverseSqrtMass, double3(0.0, 0.0, 0.0));
  ASSERT_EQ(realSpaceOnly.eigenvalues.size(), withReciprocal.eigenvalues.size());

  double maxShift = 0.0;
  for (std::size_t mode = 0; mode < realSpaceOnly.eigenvalues.size(); ++mode)
  {
    maxShift = std::max(maxShift, std::abs(realSpaceOnly.eigenvalues[mode] - withReciprocal.eigenvalues[mode]));
  }
  EXPECT_GT(maxShift, 0.0);

  // The reciprocal term preserves the three acoustic (translational) zero modes at Gamma.
  double scale = 0.0;
  for (const double value : withReciprocal.eigenvalues) scale = std::max(scale, std::abs(value));
  ASSERT_GT(scale, 0.0);
  const auto zeroModes = std::ranges::count_if(withReciprocal.eigenvalues,
                                               [&](double value) { return std::abs(value) < 1e-7 * scale; });
  EXPECT_EQ(zeroModes, 3);
}

TEST(phonon_dispersion, gamma_point_with_ewald_matches_normal_modes_multi_atom_molecules)
{
  System system = makeChargedMoleculePairSystem();

  // Two three-site molecules -> 18 flexible DOFs. The intramolecular exclusion correction (folded into
  // the real-space force constants) together with the reciprocal Fourier term must reproduce the full
  // generalized Hessian used by computeNormalModes (which includes Ewald Fourier + self + exclusion).
  const NormalModesResult modes = computeNormalModes(system);
  ASSERT_EQ(modes.numberOfModes, 18u);

  const std::vector<PhononModes> dispersion = computePhononDispersion(system, std::array<double3, 1>{double3(0, 0, 0)});
  ASSERT_EQ(dispersion.size(), 1u);
  ASSERT_EQ(dispersion[0].eigenvalues.size(), 18u);

  double scale = 0.0;
  for (const double value : modes.eigenvalues) scale = std::max(scale, std::abs(value));
  ASSERT_GT(scale, 0.0);

  for (std::size_t mode = 0; mode < modes.numberOfModes; ++mode)
  {
    EXPECT_NEAR(dispersion[0].eigenvalues[mode], modes.eigenvalues[mode], 1e-6 * scale + 1e-10) << "mode=" << mode;
  }
}

TEST(phonon_dispersion, gamma_point_with_ewald_matches_normal_modes_rigid_charged_framework)
{
  System system = makeRigidChargedFrameworkSystem();

  // Two flexible ions in a fixed charged framework -> 6 flexible DOFs. The framework charges are static
  // (no rows), yet they enter each ion's reciprocal self term; the reconstruction must still match the
  // normal-mode spectrum. There is no acoustic sum rule here because the framework is fixed.
  const NormalModesResult modes = computeNormalModes(system);
  ASSERT_EQ(modes.numberOfModes, 6u);

  const std::vector<PhononModes> dispersion = computePhononDispersion(system, std::array<double3, 1>{double3(0, 0, 0)});
  ASSERT_EQ(dispersion.size(), 1u);
  ASSERT_EQ(dispersion[0].eigenvalues.size(), 6u);

  double scale = 0.0;
  for (const double value : modes.eigenvalues) scale = std::max(scale, std::abs(value));
  ASSERT_GT(scale, 0.0);

  for (std::size_t mode = 0; mode < modes.numberOfModes; ++mode)
  {
    EXPECT_NEAR(dispersion[0].eigenvalues[mode], modes.eigenvalues[mode], 1e-6 * scale + 1e-10) << "mode=" << mode;
  }
}

TEST(phonon_dispersion, gamma_point_matches_normal_modes_rigid_molecules)
{
  System system = makeRigidMoleculeChainSystem();

  // Three rigid dumbbells -> 18 generalized DOFs (3 x [center of mass + orientation]); the linear geometry
  // yields one zero-inertia rotation per molecule. computeNormalModes mass-weights the analytic
  // generalized Hessian directly; the phonon route must reproduce the same Gamma-point squared
  // frequencies by projecting the Cartesian force constants and adding the gradient-curvature term.
  const NormalModesResult modes = computeNormalModes(system);
  ASSERT_EQ(modes.numberOfModes, 18u);
  EXPECT_EQ(modes.discardedRotationalDofs, 3u);

  const std::vector<PhononModes> dispersion = computePhononDispersion(system, std::array<double3, 1>{double3(0, 0, 0)});
  ASSERT_EQ(dispersion.size(), 1u);
  ASSERT_EQ(dispersion[0].eigenvalues.size(), 18u);

  double scale = 0.0;
  for (const double value : modes.eigenvalues) scale = std::max(scale, std::abs(value));
  ASSERT_GT(scale, 0.0);

  for (std::size_t mode = 0; mode < modes.numberOfModes; ++mode)
  {
    EXPECT_NEAR(dispersion[0].eigenvalues[mode], modes.eigenvalues[mode], 1e-6 * scale + 1e-10) << "mode=" << mode;
  }
}

TEST(phonon_dispersion, rigid_molecule_dispersion_is_time_reversal_symmetric_and_k_dependent)
{
  System system = makeRigidMoleculeChainSystem();

  const std::array<double3, 3> kPath = {double3(0.0, 0.0, 0.0), double3(0.3, 0.0, 0.0), double3(-0.3, 0.0, 0.0)};
  const std::vector<PhononModes> dispersion = computePhononDispersion(system, kPath);
  ASSERT_EQ(dispersion.size(), 3u);
  for (const PhononModes& modes : dispersion) ASSERT_EQ(modes.eigenvalues.size(), 18u);

  double scale = 0.0;
  for (const double value : dispersion[1].eigenvalues) scale = std::max(scale, std::abs(value));
  ASSERT_GT(scale, 0.0);

  // The coupling loop encloses a nonzero lattice vector, so k must move at least one squared frequency.
  double maxShift = 0.0;
  for (std::size_t mode = 0; mode < 18u; ++mode)
  {
    maxShift = std::max(maxShift, std::abs(dispersion[1].eigenvalues[mode] - dispersion[0].eigenvalues[mode]));
  }
  EXPECT_GT(maxShift, 1e-6 * scale);

  // Time-reversal symmetry: omega^2(k) = omega^2(-k).
  for (std::size_t mode = 0; mode < 18u; ++mode)
  {
    EXPECT_NEAR(dispersion[1].eigenvalues[mode], dispersion[2].eigenvalues[mode], 1e-9 * scale + 1e-12) << "mode=" << mode;
  }
}

TEST(phonon_dispersion, gamma_point_matches_normal_modes_charged_rigid_molecules)
{
  System system = makeChargedRigidMoleculeChainSystem();

  // Three rigid charged dipoles -> 18 generalized DOFs. computeNormalModes mass-weights the analytic
  // generalized Hessian (real-space + Ewald Fourier + self + exclusion). The phonon route must reproduce
  // the same Gamma-point squared frequencies by adding the reciprocal-space Ewald matrix to the Cartesian
  // force constants, projecting onto the center-of-mass/orientation coordinates, and adding the
  // gradient-curvature term (whose driving gradient already includes the reciprocal force).
  const NormalModesResult modes = computeNormalModes(system);
  ASSERT_EQ(modes.numberOfModes, 18u);
  EXPECT_EQ(modes.discardedRotationalDofs, 3u);

  const std::vector<PhononModes> dispersion = computePhononDispersion(system, std::array<double3, 1>{double3(0, 0, 0)});
  ASSERT_EQ(dispersion.size(), 1u);
  ASSERT_EQ(dispersion[0].eigenvalues.size(), 18u);

  double scale = 0.0;
  for (const double value : modes.eigenvalues) scale = std::max(scale, std::abs(value));
  ASSERT_GT(scale, 0.0);

  for (std::size_t mode = 0; mode < modes.numberOfModes; ++mode)
  {
    EXPECT_NEAR(dispersion[0].eigenvalues[mode], modes.eigenvalues[mode], 1e-6 * scale + 1e-10) << "mode=" << mode;
  }
}

TEST(phonon_dispersion, charged_rigid_molecule_dispersion_is_time_reversal_symmetric_and_k_dependent)
{
  System system = makeChargedRigidMoleculeChainSystem();

  const std::array<double3, 3> kPath = {double3(0.0, 0.0, 0.0), double3(0.3, 0.0, 0.0), double3(-0.3, 0.0, 0.0)};
  const std::vector<PhononModes> dispersion = computePhononDispersion(system, kPath);
  ASSERT_EQ(dispersion.size(), 3u);
  for (const PhononModes& modes : dispersion) ASSERT_EQ(modes.eigenvalues.size(), 18u);

  double scale = 0.0;
  for (const double value : dispersion[1].eigenvalues) scale = std::max(scale, std::abs(value));
  ASSERT_GT(scale, 0.0);

  // Both the short-ranged force constants and the reciprocal-space term carry the wrapped coupling, so k
  // must move at least one squared frequency.
  double maxShift = 0.0;
  for (std::size_t mode = 0; mode < 18u; ++mode)
  {
    maxShift = std::max(maxShift, std::abs(dispersion[1].eigenvalues[mode] - dispersion[0].eigenvalues[mode]));
  }
  EXPECT_GT(maxShift, 1e-6 * scale);

  // Time-reversal symmetry: omega^2(k) = omega^2(-k).
  for (std::size_t mode = 0; mode < 18u; ++mode)
  {
    EXPECT_NEAR(dispersion[1].eigenvalues[mode], dispersion[2].eigenvalues[mode], 1e-9 * scale + 1e-12) << "mode=" << mode;
  }
}

TEST(phonon_dispersion, gamma_point_matches_normal_modes_semi_flexible_molecules)
{
  System system = makeSemiFlexibleMoleculeChainSystem();

  // Three semi-flexible molecules -> 27 generalized DOFs (per molecule: rigid-group center of mass +
  // orientation and one Cartesian tail atom). computeNormalModes mass-weights the group-aware analytic
  // generalized Hessian directly; the phonon route must reproduce the same Gamma-point squared
  // frequencies by projecting the Cartesian force constants of the de-grouped copy (which keeps the
  // junction bond stiffness) onto the group DOFs and adding the gradient-curvature term.
  const NormalModesResult modes = computeNormalModes(system);
  ASSERT_EQ(modes.numberOfModes, 27u);
  EXPECT_EQ(modes.discardedRotationalDofs, 0u);

  const std::vector<PhononModes> dispersion = computePhononDispersion(system, std::array<double3, 1>{double3(0, 0, 0)});
  ASSERT_EQ(dispersion.size(), 1u);
  ASSERT_EQ(dispersion[0].eigenvalues.size(), 27u);

  double scale = 0.0;
  for (const double value : modes.eigenvalues) scale = std::max(scale, std::abs(value));
  ASSERT_GT(scale, 0.0);

  for (std::size_t mode = 0; mode < modes.numberOfModes; ++mode)
  {
    EXPECT_NEAR(dispersion[0].eigenvalues[mode], modes.eigenvalues[mode], 1e-6 * scale + 1e-10) << "mode=" << mode;
  }
}

TEST(phonon_dispersion, semi_flexible_dispersion_is_time_reversal_symmetric_and_k_dependent)
{
  System system = makeSemiFlexibleMoleculeChainSystem();

  const std::array<double3, 3> kPath = {double3(0.0, 0.0, 0.0), double3(0.3, 0.0, 0.0), double3(-0.3, 0.0, 0.0)};
  const std::vector<PhononModes> dispersion = computePhononDispersion(system, kPath);
  ASSERT_EQ(dispersion.size(), 3u);
  for (const PhononModes& modes : dispersion) ASSERT_EQ(modes.eigenvalues.size(), 27u);

  double scale = 0.0;
  for (const double value : dispersion[1].eigenvalues) scale = std::max(scale, std::abs(value));
  ASSERT_GT(scale, 0.0);

  // The coupling loop encloses a nonzero lattice vector, so k must move at least one squared frequency.
  double maxShift = 0.0;
  for (std::size_t mode = 0; mode < 27u; ++mode)
  {
    maxShift = std::max(maxShift, std::abs(dispersion[1].eigenvalues[mode] - dispersion[0].eigenvalues[mode]));
  }
  EXPECT_GT(maxShift, 1e-6 * scale);

  // Time-reversal symmetry: omega^2(k) = omega^2(-k).
  for (std::size_t mode = 0; mode < 27u; ++mode)
  {
    EXPECT_NEAR(dispersion[1].eigenvalues[mode], dispersion[2].eigenvalues[mode], 1e-9 * scale + 1e-12) << "mode=" << mode;
  }
}

TEST(phonon_dispersion, gamma_point_matches_normal_modes_charged_semi_flexible_molecules)
{
  System system = makeChargedSemiFlexibleMoleculeChainSystem();

  // Charged semi-flexible molecules: the reciprocal-space Ewald matrix is added to the Cartesian force
  // constants before the projection onto the group DOFs. The intramolecular exclusion pairs inside the
  // rigid core are constants of the motion; their projected force-constant blocks must cancel against
  // the gradient-curvature term, reproducing the group-aware analytic Hessian of computeNormalModes.
  const NormalModesResult modes = computeNormalModes(system);
  ASSERT_EQ(modes.numberOfModes, 27u);
  EXPECT_EQ(modes.discardedRotationalDofs, 0u);

  const std::vector<PhononModes> dispersion = computePhononDispersion(system, std::array<double3, 1>{double3(0, 0, 0)});
  ASSERT_EQ(dispersion.size(), 1u);
  ASSERT_EQ(dispersion[0].eigenvalues.size(), 27u);

  double scale = 0.0;
  for (const double value : modes.eigenvalues) scale = std::max(scale, std::abs(value));
  ASSERT_GT(scale, 0.0);

  for (std::size_t mode = 0; mode < modes.numberOfModes; ++mode)
  {
    EXPECT_NEAR(dispersion[0].eigenvalues[mode], modes.eigenvalues[mode], 1e-6 * scale + 1e-10) << "mode=" << mode;
  }
}

TEST(phonon_dispersion, ewald_exclusion_correction_contributes_at_gamma)
{
  // Without the intramolecular exclusion correction the reciprocal Fourier term over-counts the
  // intramolecular charge pairs, so the Gamma spectrum would not match the normal-mode analysis. This
  // guards that computeRealSpaceForceConstants actually folds in the exclusion for multi-atom molecules
  // by checking the reciprocal term is nonzero and the reconstruction is complete (covered above).
  System system = makeChargedMoleculePairSystem();
  const std::vector<double> inverseSqrtMass = phononInverseSqrtMasses(system);
  const std::vector<std::complex<double>> reciprocal =
      computeEwaldReciprocalDynamicalMatrix(system, inverseSqrtMass, double3(0.0, 0.0, 0.0));
  double maxMagnitude = 0.0;
  for (const std::complex<double>& value : reciprocal) maxMagnitude = std::max(maxMagnitude, std::abs(value));
  EXPECT_GT(maxMagnitude, 0.0);
}

TEST(phonon_dispersion, kpath_samples_segments_with_labels_and_monotonic_coordinate)
{
  const SimulationBox box(12.0, 12.0, 12.0);
  const std::array<PhononPathNode, 3> nodes = {PhononPathNode{double3(0.0, 0.0, 0.0), "G"},
                                               PhononPathNode{double3(0.5, 0.0, 0.0), "X"},
                                               PhononPathNode{double3(0.5, 0.5, 0.0), "M"}};

  const std::size_t pointsPerSegment = 10;
  const std::vector<PhononKPoint> path = buildPhononKPath(nodes, pointsPerSegment, box);

  // One shared start node plus pointsPerSegment per segment.
  ASSERT_EQ(path.size(), 1u + 2u * pointsPerSegment);

  // Endpoints carry their labels; the path coordinate starts at zero and never decreases.
  EXPECT_EQ(path.front().label, "G");
  EXPECT_EQ(path[pointsPerSegment].label, "X");
  EXPECT_EQ(path.back().label, "M");
  EXPECT_DOUBLE_EQ(path.front().pathCoordinate, 0.0);
  for (std::size_t index = 1; index < path.size(); ++index)
  {
    EXPECT_GE(path[index].pathCoordinate, path[index - 1].pathCoordinate) << "index=" << index;
  }
  EXPECT_GT(path.back().pathCoordinate, 0.0);

  // Interior nodes are unlabeled and the fractional coordinate interpolates linearly along a segment.
  EXPECT_TRUE(path[3].label.empty());
  EXPECT_NEAR(path[5].kFractional.x, 0.25, 1e-12);  // halfway along G->X
  EXPECT_NEAR(path[5].kFractional.y, 0.0, 1e-12);
}

TEST(phonon_dispersion, dispersion_along_path_matches_pointwise_and_gamma_normal_modes)
{
  System system = makeFlexibleChainSystem();

  const std::array<PhononPathNode, 2> nodes = {PhononPathNode{double3(0.0, 0.0, 0.0), "G"},
                                               PhononPathNode{double3(0.5, 0.0, 0.0), "X"}};
  const std::size_t pointsPerSegment = 8;
  const PhononDispersionResult dispersion = computePhononDispersionAlongPath(system, nodes, pointsPerSegment);

  ASSERT_EQ(dispersion.path.size(), 1u + pointsPerSegment);
  ASSERT_EQ(dispersion.modes.size(), dispersion.path.size());
  for (const PhononModes& modes : dispersion.modes) ASSERT_EQ(modes.eigenvalues.size(), 9u);

  // The first path point is Gamma, which must reproduce the independent normal-mode spectrum.
  const NormalModesResult gamma = computeNormalModes(system);
  double scale = 0.0;
  for (const double value : gamma.eigenvalues) scale = std::max(scale, std::abs(value));
  ASSERT_GT(scale, 0.0);
  for (std::size_t mode = 0; mode < gamma.numberOfModes; ++mode)
  {
    EXPECT_NEAR(dispersion.modes.front().eigenvalues[mode], gamma.eigenvalues[mode], 1e-7 * scale + 1e-10)
        << "mode=" << mode;
  }

  // The along-path modes must agree with a direct pointwise evaluation at the same k-points.
  std::vector<double3> kPoints;
  for (const PhononKPoint& point : dispersion.path) kPoints.push_back(point.kFractional);
  const std::vector<PhononModes> pointwise = computePhononDispersion(system, kPoints);
  ASSERT_EQ(pointwise.size(), dispersion.modes.size());
  for (std::size_t index = 0; index < pointwise.size(); ++index)
  {
    for (std::size_t mode = 0; mode < 9u; ++mode)
    {
      EXPECT_DOUBLE_EQ(dispersion.modes[index].eigenvalues[mode], pointwise[index].eigenvalues[mode])
          << "index=" << index << " mode=" << mode;
    }
  }
}

TEST(phonon_dispersion, gamma_point_with_ewald_matches_normal_modes_charged_flexible_framework)
{
  System system = makeChargedFrameworkSystem();

  // Two flexible framework atoms + one adsorbate ion -> 9 flexible DOFs. The framework bonded pair is a
  // Coulomb exclusion, so the framework exclusion correction (folded into the force constants) plus the
  // reciprocal Fourier term must reproduce the full generalized Hessian used by computeNormalModes.
  const NormalModesResult modes = computeNormalModes(system);
  ASSERT_EQ(modes.numberOfModes, 9u);

  const std::vector<PhononModes> dispersion = computePhononDispersion(system, std::array<double3, 1>{double3(0, 0, 0)});
  ASSERT_EQ(dispersion.size(), 1u);
  ASSERT_EQ(dispersion[0].eigenvalues.size(), 9u);

  double scale = 0.0;
  for (const double value : modes.eigenvalues) scale = std::max(scale, std::abs(value));
  ASSERT_GT(scale, 0.0);

  for (std::size_t mode = 0; mode < modes.numberOfModes; ++mode)
  {
    EXPECT_NEAR(dispersion[0].eigenvalues[mode], modes.eigenvalues[mode], 1e-6 * scale + 1e-10) << "mode=" << mode;
  }
}
