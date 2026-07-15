#include <gtest/gtest.h>

import std;

import int3;
import double3;
import double3x3;
import simd_quatd;
import units;
import atom;
import molecule;
import connectivity_table;
import intra_molecular_potentials;
import pseudo_atom;
import vdwparameters;
import forcefield;
import framework;
import component;
import system;
import simulationbox;
import running_energy;
import interactions_polarization_derivatives;
import normal_modes;
import phonon_force_constants;
import phonon_dynamical_matrix;

namespace
{
// Two flexible, charge-neutral, polarizable three-site molecules in a box (no framework). The reciprocal
// field vanishes without a fixed framework structure factor, so the polarization Hessian is entirely
// real-space and exercises the many-body (three-center) Term A and the third-derivative Term B.
System makePolarizableBoxSystem()
{
  ForceField forceField({{"P", false, 15.0, 0.8, 1.5, 8, false}, {"N", false, 14.0, -0.4, 1.2, 8, false}},
                        {{60.0, 3.0}, {40.0, 3.2}}, ForceField::MixingRule::Lorentz_Berthelot, 11.0, 11.0, 11.0, true,
                        false, true);
  forceField.computePolarization = true;
  forceField.omitInterPolarization = false;
  forceField.omitInterInteractions = false;
  forceField.automaticEwald = false;
  forceField.EwaldAlpha = 0.3;
  forceField.numberOfWaveVectors = int3(8, 8, 8);

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

// A rigid charged framework (two fixed +0.4/-0.4 atoms) plus two flexible, polarizable single-atom ions
// with Ewald. The framework carries no degrees of freedom but produces both a real-space field and a
// static reciprocal (structure-factor) field on the ions, so this system exercises the reciprocal-space
// on-site polarization contribution in addition to the real-space terms.
Atom chargedFrameworkAtom(double3 position, double charge, std::size_t type)
{
  return Atom(position, charge, 1.0, 0, static_cast<std::uint16_t>(type), 0, 0, true);
}

System makePolarizableRigidFrameworkSystem()
{
  ForceField forceField({{"P", false, 15.0, 0.4, 1.5, 8, false}, {"N", false, 14.0, -0.4, 1.2, 8, false}},
                        {{60.0, 3.0}, {40.0, 3.2}}, ForceField::MixingRule::Lorentz_Berthelot, 11.0, 11.0, 11.0, true,
                        false, true);
  forceField.computePolarization = true;
  forceField.omitInterPolarization = false;
  forceField.omitInterInteractions = false;
  forceField.automaticEwald = false;
  forceField.EwaldAlpha = 0.3;
  forceField.numberOfWaveVectors = int3(8, 8, 8);

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
  atoms[1].position = double3(9.0, 11.5, 13.5);
  return system;
}

// Place a rigid molecule's atoms in the lab frame from its center of mass and orientation, matching the
// convention used by the rigid kinematics (space-frame offset s = R^{-1} * body-frame position).
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

// A bent, charge-neutral, polarizable three-site molecule (body frame centered on the mass-weighted
// center of mass: two -0.4 sites of mass 14 and one +0.8 site of mass 15). Non-collinear, so all three
// rotational inertias are non-zero and every orientation degree of freedom is exercised.
Component makeRigidPolarizableComponent(const ForceField& forceField, std::size_t typeP, std::size_t typeN)
{
  ConnectivityTable connectivityTable(3);
  connectivityTable[0, 1] = true;
  connectivityTable[1, 0] = true;
  connectivityTable[1, 2] = true;
  connectivityTable[2, 1] = true;

  Component component(forceField, "rigidIon", 43.0, 1.0e6, 0.02,
                      {Atom({-0.95, 0.30, 0.0}, -0.4, 1.0, 0, static_cast<std::uint16_t>(typeN), 0, false, false),
                       Atom({0.0, -0.56, 0.0}, 0.8, 1.0, 0, static_cast<std::uint16_t>(typeP), 0, false, false),
                       Atom({0.95, 0.30, 0.0}, -0.4, 1.0, 0, static_cast<std::uint16_t>(typeN), 0, false, false)},
                      connectivityTable, Potentials::IntraMolecularPotentials{}, 5, 21);
  component.rigid = true;
  return component;
}

// Two rigid polarizable three-site molecules in a box (no framework): the polarization Hessian is entirely
// real-space (many-body Term A and third-derivative Term B) and must be projected onto each molecule's
// center-of-mass and orientation degrees of freedom, including the gradient-curvature term.
System makePolarizableRigidMoleculeBoxSystem()
{
  ForceField forceField({{"P", false, 15.0, 0.8, 1.5, 8, false}, {"N", false, 14.0, -0.4, 1.2, 8, false}},
                        {{60.0, 3.0}, {40.0, 3.2}}, ForceField::MixingRule::Lorentz_Berthelot, 11.0, 11.0, 11.0, true,
                        false, true);
  forceField.computePolarization = true;
  forceField.omitInterPolarization = false;
  forceField.omitInterInteractions = false;
  forceField.automaticEwald = false;
  forceField.EwaldAlpha = 0.3;
  forceField.numberOfWaveVectors = int3(8, 8, 8);

  const SimulationBox box(25.0, 25.0, 25.0);
  const std::size_t typeP = *forceField.findPseudoAtom("P");
  const std::size_t typeN = *forceField.findPseudoAtom("N");

  Component component = makeRigidPolarizableComponent(forceField, typeP, typeN);

  System system(forceField, box, false, 300.0, 1e4, 1.0, {}, {component}, {}, {2}, 5);
  const std::array<double3, 2> centersOfMass = {double3(10.0, 10.0, 10.0), double3(13.5, 11.5, 10.6)};
  const std::array<simd_quatd, 2> orientations = {simd_quatd(0.0, 0.0, 0.0, 1.0),
                                                  simd_quatd(0.2, -0.1, 0.15, 0.9631)};
  for (std::size_t moleculeIndex = 0; moleculeIndex < 2; ++moleculeIndex)
  {
    system.moleculeData[moleculeIndex].centerOfMassPosition = centersOfMass[moleculeIndex];
    system.moleculeData[moleculeIndex].orientation = orientations[moleculeIndex];
    setRigidMoleculePositions(system, moleculeIndex);
  }
  return system;
}

// A rigid charged framework plus a rigid polarizable three-site molecule with Ewald: the static
// reciprocal (structure-factor) field of the framework adds an on-site reciprocal polarization
// contribution on top of the real-space terms, all projected onto the molecule's rigid-body coordinates.
System makePolarizableRigidMoleculeInFrameworkSystem()
{
  ForceField forceField({{"P", false, 15.0, 0.4, 1.5, 8, false}, {"N", false, 14.0, -0.4, 1.2, 8, false}},
                        {{60.0, 3.0}, {40.0, 3.2}}, ForceField::MixingRule::Lorentz_Berthelot, 11.0, 11.0, 11.0, true,
                        false, true);
  forceField.computePolarization = true;
  forceField.omitInterPolarization = false;
  forceField.omitInterInteractions = false;
  forceField.automaticEwald = false;
  forceField.EwaldAlpha = 0.3;
  forceField.numberOfWaveVectors = int3(8, 8, 8);

  const SimulationBox box(24.0, 24.0, 24.0);
  const std::size_t typeP = *forceField.findPseudoAtom("P");
  const std::size_t typeN = *forceField.findPseudoAtom("N");

  std::vector<Atom> frameworkAtoms{chargedFrameworkAtom({0.25, 0.5, 0.5}, 0.4, typeP),
                                   chargedFrameworkAtom({0.75, 0.5, 0.5}, -0.4, typeN)};
  Framework framework(forceField, "rigid-charged", box, 1, frameworkAtoms, frameworkAtoms, {1, 1, 1});
  framework.rigid = true;

  Component component = makeRigidPolarizableComponent(forceField, typeP, typeN);

  System system(forceField, box, false, 300.0, 1e4, 1.0, {framework}, {component}, {}, {1}, 5);
  system.moleculeData[0].centerOfMassPosition = double3(9.0, 11.5, 13.5);
  system.moleculeData[0].orientation = simd_quatd(0.15, 0.2, -0.1, 0.9581);
  setRigidMoleculePositions(system, 0);
  return system;
}

// A one-dimensional chain of flexible, polarizable ions along x used for supercell folding: each unit
// cell (length 8 A) holds a +0.5 cation near x = 1 and a -0.5 anion near x = 7. With a 3.9 A cutoff the
// only interactions are between neighbouring cells (cation<->anion across the boundary, ~2 A), so the
// bands genuinely disperse and the induced-dipole field couples an atom to its periodic images. Plain
// (non-Ewald) Coulomb keeps everything short-ranged, so a unit-cell dispersion sampled on the commensurate
// mesh must fold exactly onto the Gamma spectrum of the corresponding supercell.
System makePolarizableChain(std::size_t cells)
{
  ForceField forceField({{"P", false, 20.0, 0.5, 2.0, 8, false}, {"N", false, 20.0, -0.5, 2.0, 8, false}},
                        {{20.0, 3.0}, {20.0, 3.0}}, ForceField::MixingRule::Lorentz_Berthelot, 3.9, 3.9, 3.9, false,
                        false, true);
  forceField.computePolarization = true;
  forceField.omitInterPolarization = false;
  forceField.omitInterInteractions = false;
  forceField.chargeMethod = ForceField::ChargeMethod::Coulomb;  // plain truncated Coulomb: pure real-space
  forceField.automaticEwald = false;

  const double cellLength = 8.0;
  const SimulationBox box(cellLength * static_cast<double>(cells), 10.0, 10.0);
  const std::size_t typeP = *forceField.findPseudoAtom("P");
  const std::size_t typeN = *forceField.findPseudoAtom("N");

  Component cation(forceField, "cation", 20.0, 1.0e6, 0.011,
                   {Atom({0.0, 0.0, 0.0}, 0.5, 1.0, 0, static_cast<std::uint16_t>(typeP), 0, false, false)},
                   ConnectivityTable(1), Potentials::IntraMolecularPotentials{}, 0, 0);
  cation.rigid = false;
  Component anion(forceField, "anion", 20.0, 1.0e6, 0.011,
                  {Atom({0.0, 0.0, 0.0}, -0.5, 1.0, 0, static_cast<std::uint16_t>(typeN), 0, false, false)},
                  ConnectivityTable(1), Potentials::IntraMolecularPotentials{}, 0, 0);
  anion.rigid = false;

  System system(forceField, box, false, 300.0, 1e4, 1.0, {}, {cation, anion}, {},
                {cells, cells}, 5);
  std::span<Atom> atoms = system.spanOfMoleculeAtoms();
  for (std::size_t c = 0; c < cells; ++c)
  {
    atoms[c].position = double3(cellLength * static_cast<double>(c) + 1.0, 5.0, 5.0);          // cation
    atoms[cells + c].position = double3(cellLength * static_cast<double>(c) + 7.0, 5.0, 5.0);  // anion
  }
  return system;
}

std::vector<std::uint8_t> movableMask(const System& system)
{
  const std::size_t numberOfFrameworkAtoms = system.spanOfFrameworkAtoms().size();
  const std::size_t numberOfMoleculeAtoms = system.spanOfMoleculeAtoms().size();
  const bool flexibleFramework = system.framework.has_value() && !system.framework->rigid;
  std::vector<std::uint8_t> movable(numberOfFrameworkAtoms + numberOfMoleculeAtoms, 0);
  for (std::size_t atom = 0; atom < numberOfFrameworkAtoms; ++atom) movable[atom] = flexibleFramework ? 1 : 0;
  for (std::size_t atom = 0; atom < numberOfMoleculeAtoms; ++atom) movable[numberOfFrameworkAtoms + atom] = 1;
  return movable;
}

// Reference/writable access to a Cartesian coordinate of a global atom (framework atoms first).
double& coordinate(System& system, std::size_t global, std::size_t axis)
{
  const std::size_t numberOfFrameworkAtoms = system.spanOfFrameworkAtoms().size();
  double3& position = (global < numberOfFrameworkAtoms)
                          ? system.spanOfFrameworkAtoms()[global].position
                          : system.spanOfMoleculeAtoms()[global - numberOfFrameworkAtoms].position;
  return (&position.x)[axis];
}

double polarizationEnergy(System& system) { return system.computeTotalEnergies().polarization; }

std::vector<double3> analyticGradient(const System& system, std::span<const std::uint8_t> movable)
{
  return Interactions::computePolarizationDerivatives(system, movable).gradient;
}
}  // namespace

// The analytic polarization gradient must equal a central finite difference of the polarization energy.
TEST(polarization_phonon, analytic_gradient_matches_finite_difference_box)
{
  System system = makePolarizableBoxSystem();
  system.precomputeTotalRigidEnergy();

  const std::vector<std::uint8_t> movable = movableMask(system);
  const std::vector<double3> gradient = analyticGradient(system, movable);

  double scale = 0.0;
  for (const double3& g : gradient) scale = std::max({scale, std::abs(g.x), std::abs(g.y), std::abs(g.z)});
  ASSERT_GT(scale, 0.0);

  const double h = 1.0e-4;
  for (std::size_t global = 0; global < movable.size(); ++global)
  {
    if (!movable[global]) continue;
    for (std::size_t axis = 0; axis < 3; ++axis)
    {
      const double original = coordinate(system, global, axis);
      coordinate(system, global, axis) = original + h;
      const double plus = polarizationEnergy(system);
      coordinate(system, global, axis) = original - h;
      const double minus = polarizationEnergy(system);
      coordinate(system, global, axis) = original;
      polarizationEnergy(system);

      const double finiteDifference = (plus - minus) / (2.0 * h);
      EXPECT_NEAR(finiteDifference, (&gradient[global].x)[axis], 1.0e-4 * scale + 1.0e-7)
          << "global=" << global << " axis=" << axis;
    }
  }
}

TEST(polarization_phonon, analytic_gradient_matches_finite_difference_rigid_framework)
{
  System system = makePolarizableRigidFrameworkSystem();
  system.precomputeTotalRigidEnergy();

  const std::vector<std::uint8_t> movable = movableMask(system);
  const std::vector<double3> gradient = analyticGradient(system, movable);

  double scale = 0.0;
  for (const double3& g : gradient) scale = std::max({scale, std::abs(g.x), std::abs(g.y), std::abs(g.z)});
  ASSERT_GT(scale, 0.0);

  const double h = 1.0e-4;
  for (std::size_t global = 0; global < movable.size(); ++global)
  {
    if (!movable[global]) continue;
    for (std::size_t axis = 0; axis < 3; ++axis)
    {
      const double original = coordinate(system, global, axis);
      coordinate(system, global, axis) = original + h;
      const double plus = polarizationEnergy(system);
      coordinate(system, global, axis) = original - h;
      const double minus = polarizationEnergy(system);
      coordinate(system, global, axis) = original;
      polarizationEnergy(system);

      const double finiteDifference = (plus - minus) / (2.0 * h);
      EXPECT_NEAR(finiteDifference, (&gradient[global].x)[axis], 1.0e-4 * scale + 1.0e-7)
          << "global=" << global << " axis=" << axis;
    }
  }
}

// The analytic polarization Hessian must equal a central finite difference of the analytic gradient
// (which is itself validated against the energy). This checks the Term A / Term B block assembly,
// including the reciprocal on-site term for the rigid-framework system.
static void checkHessianAgainstGradientDerivative(System& system)
{
  system.precomputeTotalRigidEnergy();
  const std::vector<std::uint8_t> movable = movableMask(system);
  const Interactions::PolarizationDerivatives derivatives =
      Interactions::computePolarizationDerivatives(system, movable);

  double scale = 0.0;
  for (const auto& [pair, block] : derivatives.hessianBlocks)
    for (const double value : block) scale = std::max(scale, std::abs(value));
  ASSERT_GT(scale, 0.0);

  const double h = 1.0e-4;
  for (std::size_t globalJ = 0; globalJ < movable.size(); ++globalJ)
  {
    if (!movable[globalJ]) continue;
    for (std::size_t axisJ = 0; axisJ < 3; ++axisJ)
    {
      const double original = coordinate(system, globalJ, axisJ);
      coordinate(system, globalJ, axisJ) = original + h;
      const std::vector<double3> plus = analyticGradient(system, movable);
      coordinate(system, globalJ, axisJ) = original - h;
      const std::vector<double3> minus = analyticGradient(system, movable);
      coordinate(system, globalJ, axisJ) = original;

      for (std::size_t globalI = 0; globalI < movable.size(); ++globalI)
      {
        if (!movable[globalI]) continue;
        const auto it = derivatives.hessianBlocks.find({globalI, globalJ});
        for (std::size_t axisI = 0; axisI < 3; ++axisI)
        {
          const double finiteDifference =
              ((&plus[globalI].x)[axisI] - (&minus[globalI].x)[axisI]) / (2.0 * h);
          const double analytic =
              (it == derivatives.hessianBlocks.end()) ? 0.0 : it->second[axisI * 3 + axisJ];
          EXPECT_NEAR(finiteDifference, analytic, 1.0e-4 * scale + 1.0e-6)
              << "I=" << globalI << " axisI=" << axisI << " J=" << globalJ << " axisJ=" << axisJ;
        }
      }
    }
  }
}

TEST(polarization_phonon, analytic_hessian_matches_finite_difference_box)
{
  System system = makePolarizableBoxSystem();
  checkHessianAgainstGradientDerivative(system);
}

TEST(polarization_phonon, analytic_hessian_matches_finite_difference_rigid_framework)
{
  System system = makePolarizableRigidFrameworkSystem();
  checkHessianAgainstGradientDerivative(system);
}

// The Gamma-point phonon spectrum (real-space force constants + reciprocal dynamical matrix + folded
// polarization) must match the normal-mode spectrum (mass-weighted generalized Hessian with the
// polarization contribution added in evaluateDerivatives).
TEST(polarization_phonon, gamma_point_matches_normal_modes_box)
{
  System system = makePolarizableBoxSystem();

  const NormalModesResult modes = computeNormalModes(system);
  const std::vector<PhononModes> dispersion =
      computePhononDispersion(system, std::array<double3, 1>{double3(0.0, 0.0, 0.0)});
  ASSERT_EQ(dispersion.size(), 1u);
  ASSERT_EQ(dispersion[0].eigenvalues.size(), modes.eigenvalues.size());

  double scale = 0.0;
  for (const double value : modes.eigenvalues) scale = std::max(scale, std::abs(value));
  ASSERT_GT(scale, 0.0);

  for (std::size_t mode = 0; mode < modes.numberOfModes; ++mode)
  {
    EXPECT_NEAR(dispersion[0].eigenvalues[mode], modes.eigenvalues[mode], 1.0e-6 * scale + 1.0e-10)
        << "mode=" << mode;
  }
}

TEST(polarization_phonon, gamma_point_matches_normal_modes_rigid_framework)
{
  System system = makePolarizableRigidFrameworkSystem();

  const NormalModesResult modes = computeNormalModes(system);
  const std::vector<PhononModes> dispersion =
      computePhononDispersion(system, std::array<double3, 1>{double3(0.0, 0.0, 0.0)});
  ASSERT_EQ(dispersion.size(), 1u);
  ASSERT_EQ(dispersion[0].eigenvalues.size(), modes.eigenvalues.size());

  double scale = 0.0;
  for (const double value : modes.eigenvalues) scale = std::max(scale, std::abs(value));
  ASSERT_GT(scale, 0.0);

  for (std::size_t mode = 0; mode < modes.numberOfModes; ++mode)
  {
    EXPECT_NEAR(dispersion[0].eigenvalues[mode], modes.eigenvalues[mode], 1.0e-6 * scale + 1.0e-10)
        << "mode=" << mode;
  }
}

// Rigorous non-Gamma validation: the image-resolved polarization force constants carry the Bloch phase
// e^{2 pi i k.R}, so a unit-cell dispersion sampled on the commensurate mesh k = j/N must fold exactly
// onto the Gamma spectrum of the N-cell supercell (a standard, implementation-independent identity).
TEST(polarization_phonon, dispersion_folds_onto_supercell_gamma)
{
  constexpr std::size_t replicas = 3;

  const System unitCell = makePolarizableChain(1);
  std::vector<double3> mesh;
  for (std::size_t j = 0; j < replicas; ++j)
  {
    mesh.push_back(double3(static_cast<double>(j) / static_cast<double>(replicas), 0.0, 0.0));
  }
  const std::vector<PhononModes> unitDispersion = computePhononDispersion(unitCell, mesh);

  std::vector<double> foldedEigenvalues;
  for (const PhononModes& modes : unitDispersion)
  {
    foldedEigenvalues.insert(foldedEigenvalues.end(), modes.eigenvalues.begin(), modes.eigenvalues.end());
  }

  const System superCell = makePolarizableChain(replicas);
  const std::vector<PhononModes> superGamma =
      computePhononDispersion(superCell, std::array<double3, 1>{double3(0.0, 0.0, 0.0)});
  ASSERT_EQ(superGamma.size(), 1u);

  std::vector<double> superEigenvalues = superGamma[0].eigenvalues;
  ASSERT_EQ(foldedEigenvalues.size(), superEigenvalues.size());

  std::ranges::sort(foldedEigenvalues);
  std::ranges::sort(superEigenvalues);

  double scale = 0.0;
  for (const double value : superEigenvalues) scale = std::max(scale, std::abs(value));
  ASSERT_GT(scale, 0.0);

  for (std::size_t index = 0; index < superEigenvalues.size(); ++index)
  {
    EXPECT_NEAR(foldedEigenvalues[index], superEigenvalues[index], 1.0e-6 * scale + 1.0e-10) << "index=" << index;
  }
}

// The polarization contribution genuinely disperses: switching it off changes the off-Gamma spectrum,
// so the folding test above is actually exercising the induced-dipole force constants (not only the
// pairwise Coulomb/van der Waals chain).
TEST(polarization_phonon, polarization_changes_off_gamma_spectrum)
{
  const std::array<double3, 1> kPoint{double3(0.5, 0.0, 0.0)};

  System withPolarization = makePolarizableChain(1);
  System withoutPolarization = makePolarizableChain(1);
  withoutPolarization.forceField.computePolarization = false;

  const std::vector<PhononModes> polarized = computePhononDispersion(withPolarization, kPoint);
  const std::vector<PhononModes> unpolarized = computePhononDispersion(withoutPolarization, kPoint);
  ASSERT_EQ(polarized[0].eigenvalues.size(), unpolarized[0].eigenvalues.size());

  double maximumDifference = 0.0;
  for (std::size_t mode = 0; mode < polarized[0].eigenvalues.size(); ++mode)
  {
    maximumDifference =
        std::max(maximumDifference, std::abs(polarized[0].eigenvalues[mode] - unpolarized[0].eigenvalues[mode]));
  }
  EXPECT_GT(maximumDifference, 1.0e-6);
}

// A polarization dynamical matrix is Hermitian with real spectrum and obeys time-reversal symmetry
// D(k) ~ D(-k); the reciprocal (fixed-framework) on-site term is exercised here at a non-Gamma k.
TEST(polarization_phonon, dispersion_is_time_reversal_symmetric)
{
  System chain = makePolarizableChain(1);
  const double3 k(0.3, 0.0, 0.0);
  const std::vector<PhononModes> forward = computePhononDispersion(chain, std::array<double3, 1>{k});
  const std::vector<PhononModes> backward = computePhononDispersion(chain, std::array<double3, 1>{double3(-k.x, -k.y, -k.z)});
  ASSERT_EQ(forward[0].eigenvalues.size(), backward[0].eigenvalues.size());

  double scale = 0.0;
  for (const double value : forward[0].eigenvalues) scale = std::max(scale, std::abs(value));
  ASSERT_GT(scale, 0.0);
  for (std::size_t mode = 0; mode < forward[0].eigenvalues.size(); ++mode)
  {
    EXPECT_NEAR(forward[0].eigenvalues[mode], backward[0].eigenvalues[mode], 1.0e-8 * scale + 1.0e-10) << "mode=" << mode;
  }
}

// Rigid molecules also disperse with polarization now (the generalized-dispersion path projects the
// image-resolved Cartesian force constants at every k). Check it runs and is time-reversal symmetric.
TEST(polarization_phonon, non_gamma_runs_for_rigid_molecule)
{
  System system = makePolarizableRigidMoleculeBoxSystem();
  const double3 k(0.25, 0.1, 0.0);
  std::vector<PhononModes> forward;
  std::vector<PhononModes> backward;
  ASSERT_NO_THROW(forward = computePhononDispersion(system, std::array<double3, 1>{k}));
  ASSERT_NO_THROW(backward =
                      computePhononDispersion(system, std::array<double3, 1>{double3(-k.x, -k.y, -k.z)}));

  double scale = 0.0;
  for (const double value : forward[0].eigenvalues) scale = std::max(scale, std::abs(value));
  ASSERT_GT(scale, 0.0);
  for (std::size_t mode = 0; mode < forward[0].eigenvalues.size(); ++mode)
  {
    EXPECT_NEAR(forward[0].eigenvalues[mode], backward[0].eigenvalues[mode], 1.0e-8 * scale + 1.0e-10) << "mode=" << mode;
  }
}

// Rigid molecules: the Gamma-point phonon spectrum (generalized-dispersion path, which projects the
// Cartesian force constants and gradient onto the rigid-body coordinates) must match the normal-mode
// spectrum (evaluateDerivatives projects the polarization Cartesian Hessian and gradient-curvature term
// directly). The two rigid projections are implemented independently, so their agreement validates the
// rigid-molecule polarization Hessian.
TEST(polarization_phonon, gamma_point_matches_normal_modes_rigid_molecule_box)
{
  System system = makePolarizableRigidMoleculeBoxSystem();

  const NormalModesResult modes = computeNormalModes(system);
  const std::vector<PhononModes> dispersion =
      computePhononDispersion(system, std::array<double3, 1>{double3(0.0, 0.0, 0.0)});
  ASSERT_EQ(dispersion.size(), 1u);
  ASSERT_EQ(dispersion[0].eigenvalues.size(), modes.eigenvalues.size());

  double scale = 0.0;
  for (const double value : modes.eigenvalues) scale = std::max(scale, std::abs(value));
  ASSERT_GT(scale, 0.0);

  for (std::size_t mode = 0; mode < modes.numberOfModes; ++mode)
  {
    EXPECT_NEAR(dispersion[0].eigenvalues[mode], modes.eigenvalues[mode], 1.0e-6 * scale + 1.0e-10) << "mode=" << mode;
  }
}

TEST(polarization_phonon, gamma_point_matches_normal_modes_rigid_molecule_in_framework)
{
  System system = makePolarizableRigidMoleculeInFrameworkSystem();

  const NormalModesResult modes = computeNormalModes(system);
  const std::vector<PhononModes> dispersion =
      computePhononDispersion(system, std::array<double3, 1>{double3(0.0, 0.0, 0.0)});
  ASSERT_EQ(dispersion.size(), 1u);
  ASSERT_EQ(dispersion[0].eigenvalues.size(), modes.eigenvalues.size());

  double scale = 0.0;
  for (const double value : modes.eigenvalues) scale = std::max(scale, std::abs(value));
  ASSERT_GT(scale, 0.0);

  for (std::size_t mode = 0; mode < modes.numberOfModes; ++mode)
  {
    EXPECT_NEAR(dispersion[0].eigenvalues[mode], modes.eigenvalues[mode], 1.0e-6 * scale + 1.0e-10) << "mode=" << mode;
  }
}

