#include <gtest/gtest.h>

#include <algorithm>
#include <complex>
#include <cstddef>
#include <numbers>
#include <print>
#include <span>
#include <vector>

import int3;
import double3;
import double3x3;
import simd_quatd;
import factory;
import units;
import atom;
import pseudo_atom;
import vdwparameters;
import forcefield;
import framework;
import component;
import system;
import simulationbox;
import energy_factor;
import running_energy;
import interactions_intermolecular;
import interactions_framework_molecule;
import interactions_ewald;

TEST(rotation_matrix_reconstruction, Test_2_CO2_in_ITQ_29_2x2x2)
{
  ForceField forceField = TestFactories::makeDefaultFF(12.0, true, false, true);
  Framework f = TestFactories::makeITQ29(forceField, int3(2, 2, 2));
  Component c = TestFactories::makeCO2(forceField, 0, true);

  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {f}, {c}, {}, {40}, 5);

  std::span<Atom> atomData = system.spanOfMoleculeAtoms();

  for (std::size_t i = 0; i != 40; ++i)
  {
    double3 com = system.moleculeData[i].centerOfMassPosition;

    std::size_t atomIndex = system.moleculeData[i].atomIndex;
    std::size_t numberOfAtoms = system.moleculeData[i].numberOfAtoms;

    std::vector<double3> positions(numberOfAtoms);
    for (std::size_t j = 0; j != numberOfAtoms; ++j)
    {
      positions[j] = atomData[atomIndex + j].position;
    }

    std::vector<double3> reference_positions(numberOfAtoms);

    for (std::size_t j = 0; j != numberOfAtoms; ++j)
    {
      reference_positions[j] = c.atoms[j].position;
    }

    // compute rotation matrix that describes going from the space-fixed frame to the body-fixed frame
    double3x3 rotation_matrix = double3x3::computeRotationMatrix(system.moleculeData[i].centerOfMassPosition, positions,
                                                                 double3{0.0, 0.0, 0.0}, reference_positions);

    simd_quatd computed_quaternion = rotation_matrix.quaternion();

    double3x3 finalRotationMatrix = double3x3::buildRotationMatrixInverse(computed_quaternion);

    for (std::size_t j = 0; j != numberOfAtoms; ++j)
    {
      EXPECT_NEAR(atomData[atomIndex + j].position.x, (com + finalRotationMatrix * c.atoms[j].position).x, 1e-6);
      EXPECT_NEAR(atomData[atomIndex + j].position.y, (com + finalRotationMatrix * c.atoms[j].position).y, 1e-6);
      EXPECT_NEAR(atomData[atomIndex + j].position.z, (com + finalRotationMatrix * c.atoms[j].position).z, 1e-6);
    }
  }
}

TEST(rotation_matrix_reconstruction, Test_2_Water_in_ITQ_29_2x2x2)
{
  ForceField forceField = ForceField(
      {PseudoAtom("O", false, 15.9996, -0.84760, 0.0, 8, true), PseudoAtom("H", false, 1.0008, 0.42380, 0.0, 1, true)},
      {VDWParameters(78.19743111, 3.16555789), VDWParameters(0.0, 1.0)}, ForceField::MixingRule::Lorentz_Berthelot,
      10.0, 10.0, 10.0, false, false, true);

  SimulationBox box = SimulationBox(30.0, 30.0, 30.0, 100.0 * (std::numbers::pi / 180.0),
                                    95.0 * (std::numbers::pi / 180.0), 75.0 * (std::numbers::pi / 180.0));

  forceField.automaticEwald = false;
  forceField.EwaldAlpha = 0.2850;
  forceField.numberOfWaveVectors = int3(7, 7, 7);

  const double kvecx = forceField.numberOfWaveVectors.x * 2.0 * std::numbers::pi / box.cell.ax;
  const double kvecy = forceField.numberOfWaveVectors.y * 2.0 * std::numbers::pi / box.cell.by;
  const double kvecz = forceField.numberOfWaveVectors.z * 2.0 * std::numbers::pi / box.cell.cz;
  forceField.reciprocalCutOffSquared = std::max({kvecx * kvecx, kvecy * kvecy, kvecz * kvecz}) * 1.00001;

  Component c = Component(0, forceField, "H2O", 304.1282, 7377300.0, 0.22394,
                          {// double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type,
                           // uint8_t componentId, uint8_t groupId
                           Atom(double3(0.00000, -0.06461, 0.00000), -0.84760, 1.0, 0, 0, 0, false, false),
                           Atom(double3(0.81649, 0.51275, 0.00000), 0.42380, 1.0, 0, 1, 0, false, false),
                           Atom(double3(-0.81649, 0.51275, 0.00000), 0.42380, 1.0, 0, 1, 0, false, false)},
                          {}, {}, 5, 21);

  System system = System(0, forceField, box, 300.0, 1e4, 1.0, {}, {c}, {}, {400}, 5);

  std::span<Atom> atomData = system.spanOfMoleculeAtoms();

  for (std::size_t i = 0; i != 400; ++i)
  {
    double3 com = system.moleculeData[i].centerOfMassPosition;

    std::size_t atomIndex = system.moleculeData[i].atomIndex;
    std::size_t numberOfAtoms = system.moleculeData[i].numberOfAtoms;

    std::vector<double3> positions(numberOfAtoms);
    for (std::size_t j = 0; j != numberOfAtoms; ++j)
    {
      positions[j] = atomData[atomIndex + j].position;
    }

    std::vector<double3> reference_positions(numberOfAtoms);

    for (std::size_t j = 0; j != numberOfAtoms; ++j)
    {
      reference_positions[j] = c.atoms[j].position;
    }

    // compute rotation matrix that describes going from the space-fixed frame to the body-fixed frame
    double3x3 rotation_matrix = double3x3::computeRotationMatrix(system.moleculeData[i].centerOfMassPosition, positions,
                                                                 double3{0.0, 0.0, 0.0}, reference_positions);

    simd_quatd computed_quaternion = rotation_matrix.quaternion();

    double3x3 finalRotationMatrix = double3x3::buildRotationMatrixInverse(computed_quaternion);

    for (std::size_t j = 0; j != numberOfAtoms; ++j)
    {
      EXPECT_NEAR(atomData[atomIndex + j].position.x, (com + finalRotationMatrix * c.atoms[j].position).x, 1e-6);
      EXPECT_NEAR(atomData[atomIndex + j].position.y, (com + finalRotationMatrix * c.atoms[j].position).y, 1e-6);
      EXPECT_NEAR(atomData[atomIndex + j].position.z, (com + finalRotationMatrix * c.atoms[j].position).z, 1e-6);
    }
  }
}
