#include <gtest/gtest.h>

#include <algorithm>
#include <complex>
#include <cstddef>
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
import energy_status;
import integrators;
import integrators_compute;
import integrators_update;
import interpolation_energy_grid;

TEST(rigid_gradient, Test_2_CO2_in_ITQ_29_2x2x2)
{
  ForceField forceField = TestFactories::makeDefaultFF(11.8, true, false, true);
  forceField.automaticEwald = false;
  forceField.EwaldAlpha = 0.25;
  forceField.numberOfWaveVectors = int3(8, 8, 8);
  Component c = TestFactories::makeCO2(forceField, 0, true);
  Framework f = TestFactories::makeITQ29(forceField, int3(2, 2, 2));

  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {f}, {c}, {2}, 5);

  std::span<Atom> atomPositions = system.spanOfMoleculeAtoms();
  atomPositions[0].position = double3(5.93355, 7.93355, 5.93355 + 1.149);
  atomPositions[1].position = double3(5.93355, 7.93355, 5.93355 + 0.0);
  atomPositions[2].position = double3(5.93355, 7.93355, 5.93355 - 1.149);
  atomPositions[3].position = double3(5.93355, 3.93355, 5.93355 + 1.149);
  atomPositions[4].position = double3(5.93355, 3.93355, 5.93355 + 0.0);
  atomPositions[5].position = double3(5.93355, 3.93355, 5.93355 - 1.149);

  system.moleculePositions[0].orientation = simd_quatd(0.0, 0.0, 0.0, 1.0);
  system.moleculePositions[1].orientation = simd_quatd(0.0, 0.0, 0.0, 1.0);

  for (Atom& atom : atomPositions)
  {
    atom.gradient = double3(0.0, 0.0, 0.0);
  }

  Integrators::updateGradients(system.spanOfMoleculeAtoms(), system.spanOfFrameworkAtoms(), system.forceField,
                               system.simulationBox, system.components, system.eik_x, system.eik_y, system.eik_z,
                               system.eik_xy, system.totalEik, system.fixedFrameworkStoredEik,
                               system.interpolationGrids, system.numberOfMoleculesPerComponent);
  Integrators::updateCenterOfMassAndQuaternionGradients(system.moleculePositions, system.spanOfMoleculeAtoms(),
                                                        system.components);

  EXPECT_NEAR(atomPositions[0].gradient.x, 0.000000000000, 1e-4);
  EXPECT_NEAR(atomPositions[0].gradient.y, 103.939706389550, 1e-4);
  EXPECT_NEAR(atomPositions[0].gradient.z, -18.706660418853, 1e-4);
  EXPECT_NEAR(atomPositions[1].gradient.x, 0.000000000000, 1e-4);
  EXPECT_NEAR(atomPositions[1].gradient.y, -574.210506496196, 1e-4);
  EXPECT_NEAR(atomPositions[1].gradient.z, 0.000000000000, 1e-4);
  EXPECT_NEAR(atomPositions[2].gradient.x, 0.000000000000, 1e-4);
  EXPECT_NEAR(atomPositions[2].gradient.y, 103.939706389549, 1e-4);
  EXPECT_NEAR(atomPositions[2].gradient.z, 18.706660418865, 1e-4);

  EXPECT_NEAR(atomPositions[3].gradient.x, 0.000000000000, 1e-4);
  EXPECT_NEAR(atomPositions[3].gradient.y, -103.939706389550, 1e-4);
  EXPECT_NEAR(atomPositions[3].gradient.z, -18.706660418853, 1e-4);
  EXPECT_NEAR(atomPositions[4].gradient.x, 0.000000000000, 1e-4);
  EXPECT_NEAR(atomPositions[4].gradient.y, 574.210506496196, 1e-4);
  EXPECT_NEAR(atomPositions[4].gradient.z, 0.000000000000, 1e-4);
  EXPECT_NEAR(atomPositions[5].gradient.x, 0.000000000000, 1e-4);
  EXPECT_NEAR(atomPositions[5].gradient.y, -103.939706389549, 1e-4);
  EXPECT_NEAR(atomPositions[5].gradient.z, 18.706660418865, 1e-4);

  EXPECT_NEAR(system.moleculePositions[0].gradient.x, 0.000000000000, 1e-4);
  EXPECT_NEAR(system.moleculePositions[0].gradient.y, -366.33109372, 1e-4);
  EXPECT_NEAR(system.moleculePositions[0].gradient.z, 0.000000000000, 1e-4);

  EXPECT_NEAR(system.moleculePositions[1].gradient.x, 0.000000000000, 1e-4);
  EXPECT_NEAR(system.moleculePositions[1].gradient.y, 366.33109372, 1e-4);
  EXPECT_NEAR(system.moleculePositions[1].gradient.z, 0.000000000000, 1e-4);

  EXPECT_NEAR(system.moleculePositions[0].orientationGradient.ix, 0.000000000000, 1e-6);
  EXPECT_NEAR(system.moleculePositions[0].orientationGradient.iy, 0.000000000000, 1e-6);
  EXPECT_NEAR(system.moleculePositions[0].orientationGradient.iz, 0.000000000000, 1e-6);
  EXPECT_NEAR(system.moleculePositions[0].orientationGradient.r, 0.000000000000, 1e-6);

  EXPECT_NEAR(system.moleculePositions[1].orientationGradient.ix, 0.000000000000, 1e-6);
  EXPECT_NEAR(system.moleculePositions[1].orientationGradient.iy, 0.000000000000, 1e-6);
  EXPECT_NEAR(system.moleculePositions[1].orientationGradient.iz, 0.000000000000, 1e-6);
  EXPECT_NEAR(system.moleculePositions[1].orientationGradient.r, 0.000000000000, 1e-6);
}

TEST(rigid_gradient, Test_2_CO2_in_ITQ_29_2x2x2_no_symmetry)
{
  ForceField forceField = TestFactories::makeDefaultFF(11.8, true, false, true);
  forceField.automaticEwald = false;
  forceField.EwaldAlpha = 0.25;
  forceField.numberOfWaveVectors = int3(8, 8, 8);
  Component c = TestFactories::makeCO2(forceField, 0, true);
  Framework f = TestFactories::makeITQ29(forceField, int3(2, 2, 2));

  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {f}, {c}, {2}, 5);

  std::span<Atom> atomPositions = system.spanOfMoleculeAtoms();
  atomPositions[0].position = double3(7.200346017629, 5.547320950450, 2.554032891860);
  atomPositions[1].position = double3(6.206010115191, 5.030510168430, 2.807811789148);
  atomPositions[2].position = double3(5.211674212752, 4.513699386409, 3.061590686436);
  atomPositions[3].position = double3(8.486922830523, 3.714749040941, 7.214116618614);
  atomPositions[4].position = double3(7.905757540326, 3.815476396960, 6.228062915605);
  atomPositions[5].position = double3(7.324592250129, 3.916203752979, 5.242009212595);

  system.moleculePositions[0].orientation = simd_quatd(0.44301474, 0.64356245, 0.55112425, 0.29296310);
  system.moleculePositions[1].orientation = simd_quatd(0.04909940, 0.26171955, 0.01333290, 0.96380203);

  for (Atom& atom : atomPositions)
  {
    atom.gradient = double3(0.0, 0.0, 0.0);
  }

  Integrators::updateGradients(system.spanOfMoleculeAtoms(), system.spanOfFrameworkAtoms(), system.forceField,
                               system.simulationBox, system.components, system.eik_x, system.eik_y, system.eik_z,
                               system.eik_xy, system.totalEik, system.fixedFrameworkStoredEik,
                               system.interpolationGrids, system.numberOfMoleculesPerComponent);
  Integrators::updateCenterOfMassAndQuaternionGradients(system.moleculePositions, system.spanOfMoleculeAtoms(),
                                                        system.components);

  EXPECT_NEAR(atomPositions[0].gradient.x, -401.918996903652, 1e-4);
  EXPECT_NEAR(atomPositions[0].gradient.y, -250.086986676442, 1e-4);
  EXPECT_NEAR(atomPositions[0].gradient.z, 76.393683260493, 1e-4);
  EXPECT_NEAR(atomPositions[1].gradient.x, -112.917384315368, 1e-4);
  EXPECT_NEAR(atomPositions[1].gradient.y, 135.881091839024, 1e-4);
  EXPECT_NEAR(atomPositions[1].gradient.z, 191.752321499793, 1e-4);
  EXPECT_NEAR(atomPositions[2].gradient.x, 649.466150035879, 1e-4);
  EXPECT_NEAR(atomPositions[2].gradient.y, 183.761020135682, 1e-4);
  EXPECT_NEAR(atomPositions[2].gradient.z, 542.179200465678, 1e-4);

  EXPECT_NEAR(atomPositions[3].gradient.x, -30.785378731879, 1e-4);
  EXPECT_NEAR(atomPositions[3].gradient.y, 80.096820989863, 1e-4);
  EXPECT_NEAR(atomPositions[3].gradient.z, 113.734455039940, 1e-4);
  EXPECT_NEAR(atomPositions[4].gradient.x, -40.135560076449, 1e-4);
  EXPECT_NEAR(atomPositions[4].gradient.y, 19.413856681000, 1e-4);
  EXPECT_NEAR(atomPositions[4].gradient.z, -302.257466181579, 1e-4);
  EXPECT_NEAR(atomPositions[5].gradient.x, -310.276701548065, 1e-4);
  EXPECT_NEAR(atomPositions[5].gradient.y, 181.549370345199, 1e-4);
  EXPECT_NEAR(atomPositions[5].gradient.z, -33.308103959230, 1e-4);

  EXPECT_NEAR(system.moleculePositions[0].gradient.x, 134.62976882, 1e-4);
  EXPECT_NEAR(system.moleculePositions[0].gradient.y, 69.55512530, 1e-4);
  EXPECT_NEAR(system.moleculePositions[0].gradient.z, 810.32520523, 1e-4);

  EXPECT_NEAR(system.moleculePositions[1].gradient.x, -381.19764036, 1e-4);
  EXPECT_NEAR(system.moleculePositions[1].gradient.y, 281.06004802, 1e-4);
  EXPECT_NEAR(system.moleculePositions[1].gradient.z, -221.83111510, 1e-4);

  EXPECT_NEAR(system.moleculePositions[0].orientationGradient.ix, 454.91960800, 1e-4);
  EXPECT_NEAR(system.moleculePositions[0].orientationGradient.iy, 913.61679756, 1e-4);
  EXPECT_NEAR(system.moleculePositions[0].orientationGradient.iz, -1032.71684255, 1e-4);
  EXPECT_NEAR(system.moleculePositions[0].orientationGradient.r, -752.14337829, 1e-4);

  EXPECT_NEAR(system.moleculePositions[1].orientationGradient.ix, 185.47928131, 1e-4);
  EXPECT_NEAR(system.moleculePositions[1].orientationGradient.iy, 361.21170189, 1e-4);
  EXPECT_NEAR(system.moleculePositions[1].orientationGradient.iz, -33.44656884, 1e-4);
  EXPECT_NEAR(system.moleculePositions[1].orientationGradient.r, -107.07297003, 1e-4);
}

TEST(rigid_gradient, Test_2_H2O_in_ITQ_29_2x2x2_no_symmetry)
{
  ForceField forceField = TestFactories::makeDefaultFF(11.8, true, false, true);
  forceField.automaticEwald = false;
  forceField.EwaldAlpha = 0.25;
  forceField.numberOfWaveVectors = int3(8, 8, 8);
  Component c = TestFactories::makeWater(forceField, 0, true);
  Framework f = TestFactories::makeITQ29(forceField, int3(2, 2, 2));

  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {f}, {c}, {2}, 5);

  std::span<Atom> atomPositions = system.spanOfMoleculeAtoms();
  atomPositions[0].position = double3(2.672476355982, 4.866114600102, 3.294027041381);
  atomPositions[1].position = double3(2.884778590882, 3.937807845649, 3.197041626529);
  atomPositions[2].position = double3(3.522761084928, 5.304109388792, 3.331554162278);
  atomPositions[3].position = double3(2.310905951524, 4.976954564602, 3.883078265715);
  atomPositions[4].position = double3(2.301056712725, 5.093499932362, 2.745991130810);
  atomPositions[5].position = double3(8.337288460238, 7.804054923450, 9.508718172285);
  atomPositions[6].position = double3(8.438085808732, 7.132544409030, 8.834073552088);
  atomPositions[7].position = double3(8.923408607748, 7.526719865362, 10.212834837504);
  atomPositions[8].position = double3(7.675123524814, 7.830783863478, 9.734158102333);
  atomPositions[9].position = double3(8.525606355084, 8.431855339730, 9.262947938533);

  double3 com0, com1;
  for (size_t i = 0; i < 5; ++i)
  {
    com0 += atomPositions[i].position * c.definedAtoms[i].second;
    com1 += atomPositions[i].position * c.definedAtoms[i].second;
  }
  system.moleculePositions[0].centerOfMassPosition = com0 / system.moleculePositions[0].invMass;
  system.moleculePositions[1].centerOfMassPosition = com1 / system.moleculePositions[0].invMass;

  atomPositions[0].velocity = double3(-2.56934959, 3.49026313, 2.46235469);
  atomPositions[1].velocity = double3(14.83854769, 8.31306952, -5.59358085);
  atomPositions[2].velocity = double3(-10.83853213, 19.76930089, -0.17504625);
  atomPositions[3].velocity = double3(-0.73607131, -7.74519202, 5.70179557);
  atomPositions[4].velocity = double3(-10.70666431, 0.16931747, 6.59935402);
  atomPositions[5].velocity = double3(0.95342220, -5.47957843, -2.32188856);
  atomPositions[6].velocity = double3(11.57668989, -11.62828553, 5.38545620);
  atomPositions[7].velocity = double3(10.25681953, 15.34872196, -1.86243447);
  atomPositions[8].velocity = double3(1.52395494, -16.52600113, 0.66358941);
  atomPositions[9].velocity = double3(-13.36285337, -4.55938179, -10.94095934);

  system.moleculePositions[0].orientation = simd_quatd(-0.34788397, 0.57840378, 0.41201608, 0.61210176);
  system.moleculePositions[1].orientation = simd_quatd(-0.07669892, 0.19685049, 0.29872787, 0.93066042);

  for (Atom& atom : atomPositions)
  {
    atom.gradient = double3(0.0, 0.0, 0.0);
  }

  RunningEnergy energy = system.computeTotalEnergies();
  EXPECT_NEAR(energy.frameworkMoleculeVDW * Units::EnergyToKelvin, 13286.65000939, 1e-4);
  EXPECT_NEAR(energy.frameworkMoleculeCharge * Units::EnergyToKelvin, 2340.77151436, 1e-4);
  EXPECT_NEAR(energy.moleculeMoleculeVDW * Units::EnergyToKelvin, -0.51506625, 1e-4);
  EXPECT_NEAR(energy.moleculeMoleculeCharge * Units::EnergyToKelvin, 1.94101152, 1e-4);
  EXPECT_NEAR((energy.ewald_fourier + energy.ewald_self + energy.ewald_exclusion) * Units::EnergyToKelvin,
              -578.26684244 + 18.14933840, 1e-4);
  EXPECT_NEAR(energy.potentialEnergy() * Units::EnergyToKelvin, 15068.729964973694, 1e-4);

  RunningEnergy force = Integrators::updateGradients(
      system.spanOfMoleculeAtoms(), system.spanOfFrameworkAtoms(), system.forceField, system.simulationBox,
      system.components, system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.totalEik,
      system.fixedFrameworkStoredEik, system.interpolationGrids, system.numberOfMoleculesPerComponent);
  Integrators::updateCenterOfMassAndQuaternionGradients(system.moleculePositions, system.spanOfMoleculeAtoms(),
                                                        system.components);
  Integrators::updateCenterOfMassAndQuaternionVelocities(system.moleculePositions, system.spanOfMoleculeAtoms(),
                                                         system.components);

  EXPECT_NEAR(force.frameworkMoleculeVDW * Units::EnergyToKelvin, 13286.65000939, 1e-4);
  EXPECT_NEAR(force.frameworkMoleculeCharge * Units::EnergyToKelvin, 2340.77151436, 1e-4);
  EXPECT_NEAR(force.moleculeMoleculeVDW * Units::EnergyToKelvin, -0.51506625, 1e-4);
  EXPECT_NEAR(force.moleculeMoleculeCharge * Units::EnergyToKelvin, 1.94101152, 1e-4);
  EXPECT_NEAR((energy.ewald_fourier + energy.ewald_self + energy.ewald_exclusion) * Units::EnergyToKelvin,
              -578.26684244 + 18.14933840, 1e-4);
  EXPECT_NEAR(energy.potentialEnergy() * Units::EnergyToKelvin, 15068.729964973694, 1e-4);

  EXPECT_NEAR(atomPositions[0].gradient.x, -10324.973732049806, 1e-4);
  EXPECT_NEAR(atomPositions[0].gradient.y, -4507.328164443007, 1e-4);
  EXPECT_NEAR(atomPositions[0].gradient.z, -11740.007010211551, 1e-4);
  EXPECT_NEAR(atomPositions[1].gradient.x, 63.952520769129, 1e-4);
  EXPECT_NEAR(atomPositions[1].gradient.y, 357.934385984066, 1e-4);
  EXPECT_NEAR(atomPositions[1].gradient.z, 229.902658234023, 1e-4);
  EXPECT_NEAR(atomPositions[2].gradient.x, -177.459324021042, 1e-4);
  EXPECT_NEAR(atomPositions[2].gradient.y, 13.702574609408, 1e-4);
  EXPECT_NEAR(atomPositions[2].gradient.z, -149.823840142605, 1e-4);
  EXPECT_NEAR(atomPositions[3].gradient.x, 40.721496369440, 1e-4);
  EXPECT_NEAR(atomPositions[3].gradient.y, 102.350781724646, 1e-4);
  EXPECT_NEAR(atomPositions[3].gradient.z, 544.715356017259, 1e-4);
  EXPECT_NEAR(atomPositions[4].gradient.x, 893.786993301498, 1e-4);
  EXPECT_NEAR(atomPositions[4].gradient.y, -462.571049611041, 1e-4);
  EXPECT_NEAR(atomPositions[4].gradient.z, 264.954428493690, 1e-4);
  EXPECT_NEAR(atomPositions[5].gradient.x, 21048.440564961642, 1e-4);
  EXPECT_NEAR(atomPositions[5].gradient.y, 7292.526006699147, 1e-4);
  EXPECT_NEAR(atomPositions[5].gradient.z, 38006.857957389766, 1e-4);
  EXPECT_NEAR(atomPositions[6].gradient.x, 144.018623861521, 1e-4);
  EXPECT_NEAR(atomPositions[6].gradient.y, -164.047506662763, 1e-4);
  EXPECT_NEAR(atomPositions[6].gradient.z, 122.652248106359, 1e-4);
  EXPECT_NEAR(atomPositions[7].gradient.x, 886.295335930291, 1e-4);
  EXPECT_NEAR(atomPositions[7].gradient.y, -1795.979021143564, 1e-4);
  EXPECT_NEAR(atomPositions[7].gradient.z, 5898.365905521006, 1e-4);
  EXPECT_NEAR(atomPositions[8].gradient.x, -65.898172944280, 1e-4);
  EXPECT_NEAR(atomPositions[8].gradient.y, -252.228439919443, 1e-4);
  EXPECT_NEAR(atomPositions[8].gradient.z, 73.787635477924, 1e-4);
  EXPECT_NEAR(atomPositions[9].gradient.x, 435.784773865589, 1e-4);
  EXPECT_NEAR(atomPositions[9].gradient.y, 424.010966892300, 1e-4);
  EXPECT_NEAR(atomPositions[9].gradient.z, 26.760057063541, 1e-4);

  EXPECT_NEAR(system.moleculePositions[0].gradient.x, -9503.97204563, 1e-4);
  EXPECT_NEAR(system.moleculePositions[0].gradient.y, -4495.91147174, 1e-4);
  EXPECT_NEAR(system.moleculePositions[0].gradient.z, -10850.25840761, 1e-4);

  EXPECT_NEAR(system.moleculePositions[1].gradient.x, 22448.64112567, 1e-4);
  EXPECT_NEAR(system.moleculePositions[1].gradient.y, 5504.28200587, 1e-4);
  EXPECT_NEAR(system.moleculePositions[1].gradient.z, 44128.42380356, 1e-4);

  EXPECT_NEAR(system.moleculePositions[0].velocity.x, -2.06127092, 1e-4);
  EXPECT_NEAR(system.moleculePositions[0].velocity.y, 4.66344746, 1e-4);
  EXPECT_NEAR(system.moleculePositions[0].velocity.z, 1.86784481, 1e-4);

  EXPECT_NEAR(system.moleculePositions[1].velocity.x, 2.06127092, 1e-4);
  EXPECT_NEAR(system.moleculePositions[1].velocity.y, -4.66344746, 1e-4);
  EXPECT_NEAR(system.moleculePositions[1].velocity.z, -1.86784481, 1e-4);

  EXPECT_NEAR(system.moleculePositions[0].orientationGradient.ix, -2303.05396459, 1e-4);
  EXPECT_NEAR(system.moleculePositions[0].orientationGradient.iy, -829.64716058, 1e-4);
  EXPECT_NEAR(system.moleculePositions[0].orientationGradient.iz, -544.25501509, 1e-4);
  EXPECT_NEAR(system.moleculePositions[0].orientationGradient.r, -158.60548902, 1e-4);

  EXPECT_NEAR(system.moleculePositions[1].orientationGradient.ix, 3659.18641846, 1e-4);
  EXPECT_NEAR(system.moleculePositions[1].orientationGradient.iy, -3391.54768162, 1e-4);
  EXPECT_NEAR(system.moleculePositions[1].orientationGradient.iz, -3436.34304600, 1e-4);
  EXPECT_NEAR(system.moleculePositions[1].orientationGradient.r, 2121.95001708, 1e-4);

  EXPECT_NEAR(system.moleculePositions[0].orientationMomentum.ix, -30.13618171, 1e-4);
  EXPECT_NEAR(system.moleculePositions[0].orientationMomentum.iy, -26.22576279, 1e-4);
  EXPECT_NEAR(system.moleculePositions[0].orientationMomentum.iz, 48.29543567, 1e-4);
  EXPECT_NEAR(system.moleculePositions[0].orientationMomentum.r, -24.85421796, 1e-4);

  EXPECT_NEAR(system.moleculePositions[1].orientationMomentum.ix, -52.66282302, 1e-4);
  EXPECT_NEAR(system.moleculePositions[1].orientationMomentum.iy, 8.54876925, 1e-4);
  EXPECT_NEAR(system.moleculePositions[1].orientationMomentum.iz, 27.14344361, 1e-4);
  EXPECT_NEAR(system.moleculePositions[1].orientationMomentum.r, -14.86096737, 1e-4);
}
