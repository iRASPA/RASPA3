#include <gtest/gtest.h>

import std;

import double3;
import atom;
import molecule;
import forcefield;
import component;
import system;
import simulationbox;
import cross_links;

namespace
{
ForceField makeParticleForceField()
{
  return ForceField({{"X", false, 16.0, 0.0, 0.0, 6, false}}, {{120.0, 3.5}}, ForceField::MixingRule::Lorentz_Berthelot,
                    9.0, 9.0, 9.0, true, false, false);
}

// a box of single-site Lennard-Jones particles
System makeParticleSystem(const ForceField& forceField, std::size_t numberOfParticles)
{
  const Component particle(forceField, "particle", 190.0, 4.6e6, 0.0,
                           {Atom({0.0, 0.0, 0.0}, 0.0, 1.0, 0, 0, 0, false, false)}, {}, {}, 5, 21);
  return System(forceField, SimulationBox(22.0, 22.0, 22.0), false, 300.0, 1.0e5, 1.0, {}, {particle}, {},
                {numberOfParticles}, 5);
}
}  // namespace

// A probe pseudo-atom without a mass still gets a finite reference geometry; the mass is set afterwards.
TEST(system_bookkeeping, a_massless_spherical_probe_has_a_finite_geometry)
{
  const ForceField forceField({{"-", false, 0.0, 0.0, 0.0, 0, false}}, {{0.0, 1.0}},
                              ForceField::MixingRule::Lorentz_Berthelot, 9.0, 9.0, 9.0, false, false, false);
  const Component probe = Component::makeSphericalProbe(forceField, 0, "-");
  EXPECT_DOUBLE_EQ(probe.totalMass, 28.0134);
  // (the NaN of 0/0 prints as "nan"; isnan folds away under the fast-math flags of the project)
  EXPECT_EQ(std::to_string(probe.atoms[0].position.x).find("nan"), std::string::npos);
  ASSERT_FALSE(probe.fragmentGraph.fragments.empty());
  EXPECT_EQ(std::to_string(probe.fragmentGraph.fragments[0].centerOfMassReferencePosition.x).find("nan"),
            std::string::npos);
}

// Changing the angles of a rectangular box makes it triclinic, so the minimum image uses the full cell.
TEST(system_bookkeeping, setting_the_box_angles_updates_the_box_type)
{
  SimulationBox box(10.0, 10.0, 10.0);
  ASSERT_EQ(box.type, SimulationBox::Type::Rectangular);
  const double3 angles(0.5 * std::numbers::pi, 0.5 * std::numbers::pi, std::numbers::pi / 3.0);
  box.setBoxAngles(angles);
  EXPECT_EQ(box.type, SimulationBox::Type::Triclinic);

  const SimulationBox reference(10.0, 10.0, 10.0, angles.x, angles.y, angles.z);
  const double3 wrapped = box.applyPeriodicBoundaryConditions(double3(7.0, 2.0, 0.0));
  const double3 expected = reference.applyPeriodicBoundaryConditions(double3(7.0, 2.0, 0.0));
  EXPECT_NEAR(wrapped.x, expected.x, 1e-12);
  EXPECT_NEAR(wrapped.y, expected.y, 1e-12);
  EXPECT_NEAR(wrapped.z, expected.z, 1e-12);

  box.setBoxAngles(double3(0.5 * std::numbers::pi, 0.5 * std::numbers::pi, 0.5 * std::numbers::pi));
  EXPECT_EQ(box.type, SimulationBox::Type::Rectangular);
}

// A configuration exchange moves the fractional-molecule bookkeeping and every lambda coordinate with the atoms.
TEST(system_bookkeeping, exchanging_configurations_swaps_the_fractional_counts_and_every_lambda_bin)
{
  const ForceField forceField = makeParticleForceField();
  System a = makeParticleSystem(forceField, 4);
  System b = makeParticleSystem(forceField, 4);
  a.numberOfFractionalMoleculesPerComponent[0] = 1;
  b.numberOfFractionalMoleculesPerComponent[0] = 0;
  a.components[0].lambdaGibbs.currentBin = 7;
  b.components[0].lambdaGibbs.currentBin = 2;
  a.components[0].lambdaPairSwap.currentBin = 3;
  b.components[0].lambdaPairSwap.currentBin = 9;
  a.components[0].lambdaGroupSwapCB.currentBin = 5;
  b.components[0].lambdaGroupSwapCB.currentBin = 1;

  a.exchangeConfiguration(b);

  EXPECT_EQ(a.numberOfFractionalMoleculesPerComponent[0], 0uz);
  EXPECT_EQ(b.numberOfFractionalMoleculesPerComponent[0], 1uz);
  EXPECT_EQ(a.components[0].lambdaGibbs.currentBin, 2uz);
  EXPECT_EQ(b.components[0].lambdaGibbs.currentBin, 7uz);
  EXPECT_EQ(a.components[0].lambdaPairSwap.currentBin, 9uz);
  EXPECT_EQ(b.components[0].lambdaPairSwap.currentBin, 3uz);
  EXPECT_EQ(a.components[0].lambdaGroupSwapCB.currentBin, 1uz);
  EXPECT_EQ(b.components[0].lambdaGroupSwapCB.currentBin, 5uz);
}

// Moving a molecule that carries a cross-link renumbers the link sites instead of refusing the move.
TEST(system_bookkeeping, moving_a_cross_linked_molecule_renumbers_its_links)
{
  const ForceField forceField = makeParticleForceField();
  System system = makeParticleSystem(forceField, 6);
  system.crossLinks.addLink(CrossLink{.a = CrossLinkSite{.componentId = 0, .moleculeIndex = 1, .atomIndex = 0},
                                      .b = CrossLinkSite{.componentId = 0, .moleculeIndex = 3, .atomIndex = 0},
                                      .bondTypeId = 0});
  ASSERT_TRUE(system.crossLinks.moleculeIsLinked(0, 1));

  // molecule 1 moves behind molecule 4: molecules 2..4 shift down by one
  EXPECT_NO_THROW(system.moveMolecule(0, 1, 4));
  EXPECT_TRUE(system.crossLinks.moleculeIsLinked(0, 4));
  EXPECT_TRUE(system.crossLinks.moleculeIsLinked(0, 2));
  EXPECT_FALSE(system.crossLinks.moleculeIsLinked(0, 1));
  EXPECT_FALSE(system.crossLinks.moleculeIsLinked(0, 3));
}

// A single molecule whose centre-of-mass momentum is constrained has no translational degrees of freedom left; the
// status reports print a zero temperature rather than 0/0.
TEST(system_bookkeeping, md_status_reports_of_a_system_without_free_degrees_of_freedom_are_finite)
{
  const ForceField forceField = makeParticleForceField();
  System system = makeParticleSystem(forceField, 1);
  system.translationalCenterOfMassConstraint = system.translationalDegreesOfFreedom;
  for (const std::string& report :
       {system.writeEquilibrationStatusReportMD(1, 10), system.writeProductionStatusReportMD(1, 10)})
  {
    // the temperature lines only: the drift line of a fresh system has its own undefined reference
    std::istringstream stream(report);
    std::size_t temperatureLines = 0;
    for (std::string line; std::getline(stream, line);)
    {
      if (!line.contains("emperature")) continue;
      ++temperatureLines;
      EXPECT_EQ(line.find("nan"), std::string::npos) << line;
      EXPECT_EQ(line.find("inf"), std::string::npos) << line;
    }
    EXPECT_GE(temperatureLines, 3uz);
  }
}

// With a pinned lambda coordinate other than the grand-canonical one, the status line reports the pinned value.
TEST(system_bookkeeping, status_reports_print_the_pinned_lambda_coordinate)
{
  const ForceField forceField = makeParticleForceField();
  System system = makeParticleSystem(forceField, 4);
  Component& particle = system.components[0];
  particle.fixedLambdaBin = 10;
  particle.fixedLambdaCoordinate = Component::FixedLambdaCoordinate::PairSwap;
  particle.lambdaPairSwap.currentBin = 10;
  particle.lambdaPairSwap.computeDUdlambda = false;
  const std::string expected = std::format("lambda: {: g}", particle.lambdaPairSwap.lambdaValue());
  ASSERT_NE(particle.lambdaPairSwap.lambdaValue(), particle.lambdaGC.lambdaValue());

  EXPECT_NE(system.writeEquilibrationStatusReportMC(1, 10).find(expected), std::string::npos);
  EXPECT_NE(system.writeEquilibrationStatusReportMD(1, 10).find(expected), std::string::npos);
  EXPECT_NE(system.writeProductionStatusReportMC("", std::nullopt).find(expected), std::string::npos);
}
