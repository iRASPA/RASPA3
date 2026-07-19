#include <gtest/gtest.h>

import std;

import int3;
import double3;
import simd_quatd;
import units;
import atom;
import pseudo_atom;
import vdwparameters;
import forcefield;
import framework;
import component;
import molecule;
import system;
import simulationbox;

// The combined translation-rotation smart MC move applies a translation followed by a rotation about
// the moved pivot in a single trial move. Its Metropolis-Hastings correction assumes that (1) the
// reverse proposal (-Delta r, -Delta phi) maps the trial state exactly back onto the old state, and
// (2) the order of the two operations does not matter (they commute because the rotation pivots about
// a material point of the molecule). These tests verify both properties for a rigid molecule.

static simd_quatd quaternionFromAngularDisplacement(double3 deltaPhi)
{
  const double angle = deltaPhi.length();
  if (angle < 1.0e-14)
  {
    return simd_quatd(0.0, 0.0, 0.0, 1.0);
  }
  return simd_quatd::fromAxisAngle(angle, deltaPhi * (1.0 / angle));
}

TEST(translation_rotation_smart_mc, reverse_proposal_recovers_old_state)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  Framework f = Framework::makeITQ29(forceField, int3(1, 1, 1));
  Component c = Component::makeCO2(forceField, 0, true);

  System system = System(forceField, std::nullopt, false, 300.0, 1e4, 1.0, {f}, {c}, {}, {1}, 5);

  std::span<Atom> atoms = system.spanOfMolecule(0, 0);
  const Molecule molecule = system.moleculeData[system.moleculeIndexOfComponent(0, 0)];
  const std::vector<Atom> originalAtoms(atoms.begin(), atoms.end());

  const double3 displacement(0.37, -0.21, 0.55);
  const double3 angularDisplacement(0.42, 0.13, -0.31);
  const simd_quatd rotation = quaternionFromAngularDisplacement(angularDisplacement);
  const simd_quatd inverseRotation = quaternionFromAngularDisplacement(-angularDisplacement);

  // Forward: translate, then rotate about the moved pivot.
  std::pair<Molecule, std::vector<Atom>> translated = c.translate(molecule, atoms, displacement);
  std::pair<Molecule, std::vector<Atom>> trial = c.rotate(translated.first, translated.second, rotation);

  // Reverse proposal generated at the trial state: translate by -Delta r, then rotate by the
  // quaternion formed from -Delta phi.
  std::pair<Molecule, std::vector<Atom>> reverseTranslated =
      c.translate(trial.first, trial.second, -displacement);
  std::pair<Molecule, std::vector<Atom>> recovered =
      c.rotate(reverseTranslated.first, reverseTranslated.second, inverseRotation);

  for (std::size_t i = 0; i != originalAtoms.size(); ++i)
  {
    EXPECT_NEAR(recovered.second[i].position.x, originalAtoms[i].position.x, 1e-10);
    EXPECT_NEAR(recovered.second[i].position.y, originalAtoms[i].position.y, 1e-10);
    EXPECT_NEAR(recovered.second[i].position.z, originalAtoms[i].position.z, 1e-10);
  }
  EXPECT_NEAR(recovered.first.centerOfMassPosition.x, molecule.centerOfMassPosition.x, 1e-10);
  EXPECT_NEAR(recovered.first.centerOfMassPosition.y, molecule.centerOfMassPosition.y, 1e-10);
  EXPECT_NEAR(recovered.first.centerOfMassPosition.z, molecule.centerOfMassPosition.z, 1e-10);
}

TEST(translation_rotation_smart_mc, translation_and_rotation_commute)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  Framework f = Framework::makeITQ29(forceField, int3(1, 1, 1));
  Component c = Component::makeCO2(forceField, 0, true);

  System system = System(forceField, std::nullopt, false, 300.0, 1e4, 1.0, {f}, {c}, {}, {1}, 5);

  std::span<Atom> atoms = system.spanOfMolecule(0, 0);
  const Molecule molecule = system.moleculeData[system.moleculeIndexOfComponent(0, 0)];

  const double3 displacement(-0.18, 0.44, 0.29);
  const double3 angularDisplacement(-0.25, 0.36, 0.17);
  const simd_quatd rotation = quaternionFromAngularDisplacement(angularDisplacement);

  // Translate then rotate.
  std::pair<Molecule, std::vector<Atom>> translated = c.translate(molecule, atoms, displacement);
  std::pair<Molecule, std::vector<Atom>> translateRotate = c.rotate(translated.first, translated.second, rotation);

  // Rotate then translate.
  std::pair<Molecule, std::vector<Atom>> rotated = c.rotate(molecule, atoms, rotation);
  std::pair<Molecule, std::vector<Atom>> rotateTranslate = c.translate(rotated.first, rotated.second, displacement);

  for (std::size_t i = 0; i != translateRotate.second.size(); ++i)
  {
    EXPECT_NEAR(translateRotate.second[i].position.x, rotateTranslate.second[i].position.x, 1e-10);
    EXPECT_NEAR(translateRotate.second[i].position.y, rotateTranslate.second[i].position.y, 1e-10);
    EXPECT_NEAR(translateRotate.second[i].position.z, rotateTranslate.second[i].position.z, 1e-10);
  }
}
