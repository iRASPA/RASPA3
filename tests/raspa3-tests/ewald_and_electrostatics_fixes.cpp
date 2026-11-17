#include <gtest/gtest.h>

import std;

import double3;
import units;
import atom;
import forcefield;
import component;
import system;
import simulationbox;
import potential_electrostatics;
import potential_pair_coulomb;
import potential_coulomb_real_space;

// Regression tests of the Ewald and electrostatics fixes: the electrostatic potential of a scaled pair is that of
// its energy, and a charged rigid molecule that reaches past half the box is refused.

namespace
{
Component rigidCO2(const ForceField& forceField, double halfLength)
{
  return Component(forceField, "CO2", 304.1282, 7377300.0, 0.22394,
                   {Atom({0, 0, halfLength}, -0.3256, 1.0, 0, 4, 0, false, false),
                    Atom({0, 0, 0.000}, 0.6512, 1.0, 0, 3, 0, false, false),
                    Atom({0, 0, -halfLength}, -0.3256, 1.0, 0, 4, 0, false, false)},
                   {}, {}, 5, 21, {}, std::nullopt, false);
}
}  // namespace

// The real-space Ewald energy of a pair with scaled charges is taken at r + delta (1 - s_A s_B), which keeps it
// finite for a switched-off site on top of another one. The electrostatic potential was taken at r: it was the
// potential of another energy, and 0/0 at r = 0.
TEST(EWALD_AND_ELECTROSTATICS, the_potential_of_a_scaled_pair_is_that_of_its_energy)
{
  const ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  ASSERT_EQ(forceField.chargeMethod, ForceField::ChargeMethod::Ewald);

  for (const double scalingA : {1.0, 0.6, 0.25})
  {
    for (const double r : {0.5, 2.0, 7.5})
    {
      const double scalingB = 1.0;
      const double chargeA = 0.7, chargeB = -0.4;
      // the energy of the pair is the scaled charge of A times the potential of B at A
      const double energy = Potentials::potentialCoulomb<0>(forceField, scalingA, scalingB, r, chargeA, chargeB).energy;
      const double potential =
          scalingB * chargeB * Potentials::potentialElectrostatics(forceField, 1.0, r, 1.0, scalingA * scalingB);
      EXPECT_NEAR(scalingA * chargeA * potential, energy, 1.0e-6 * std::abs(energy)) << scalingA << " " << r;
    }
  }

  // a switched-off site on top of another one: finite
  EXPECT_TRUE(std::isfinite(Potentials::potentialElectrostatics(forceField, 1.0, 0.0, 1.0, 0.0)));
  // the exclusion kernel of two sites on top of each other is the limit 2 alpha / sqrt(pi)
  EXPECT_NEAR(Potentials::ewaldExclusionFactors(forceField.EwaldAlpha, 1.0, 0.0).potential,
              2.0 * forceField.EwaldAlpha * std::numbers::inv_sqrtpi, 1.0e-12);
  EXPECT_TRUE(std::isfinite(Potentials::ewaldExclusionFactors(forceField.EwaldAlpha, 1.0, 0.0).firstDerivativeFactor));
}

// The exclusion of an intramolecular pair takes the direct interaction of the pair out of the Coulomb sum: the
// term at the separation of the two atoms in the molecule. At the minimum image it was, for a rigid molecule that
// reaches past half the box, the term of the nearest periodic image of the partner, which the orientation picks:
// the exclusion energy of a molecule depended on its orientation, with a cusp wherever the nearest image switched.
TEST(EWALD_AND_ELECTROSTATICS, the_exclusion_of_a_molecule_longer_than_half_the_box_is_at_its_own_separation)
{
  const ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  ASSERT_TRUE(forceField.useCharge);

  // the outer atoms are 16 apart, past half the box of 26
  const double halfLength = 8.0;
  System system(forceField, SimulationBox(26.0, 26.0, 26.0), false, 300.0, 1e4, 1.0, {}, {rigidCO2(forceField, halfLength)},
                {}, {1}, 5);
  ASSERT_EQ(system.spanOfMoleculeAtoms().size(), 3uz);

  const auto exclusionEnergy = [&](const double3& axis)
  {
    const double3 centre(13.0, 13.0, 13.0);
    const double3 direction = axis.normalized();
    std::span<Atom> atoms = system.spanOfMoleculeAtoms();
    atoms[0].position = centre + halfLength * direction;
    atoms[1].position = centre;
    atoms[2].position = centre - halfLength * direction;
    return system.computeTotalEnergies().ewald_exclusion;
  };

  // along a cell axis the outer pair is 16 apart and 10 from the image; along the body diagonal no component of
  // the separation reaches half the box
  const double alongAxis = exclusionEnergy(double3(0.0, 0.0, 1.0));
  const double alongDiagonal = exclusionEnergy(double3(1.0, 1.0, 1.0));
  EXPECT_NEAR(alongAxis, alongDiagonal, 1.0e-9 * std::abs(alongDiagonal));

  // -C sum q_i q_j erf(alpha r)/r over the three pairs of the molecule
  const double alpha = system.forceField.EwaldAlpha;
  const double qO = -0.3256, qC = 0.6512;
  const double expected =
      -Units::CoulombicConversionFactor * (2.0 * qO * qC * std::erf(alpha * halfLength) / halfLength +
                                           qO * qO * std::erf(alpha * 2.0 * halfLength) / (2.0 * halfLength));
  EXPECT_NEAR(alongAxis, expected, 1.0e-9 * std::abs(expected));
}
