#include <gtest/gtest.h>

#include <algorithm>
#include <complex>
#include <cstddef>
#include <span>
#include <tuple>
#include <vector>

import int3;
import double3;
import double3x3;
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
import energy_status;
import running_energy;
import interactions_intermolecular;
import interactions_framework_molecule;
import interactions_ewald;

TEST(MC_strain_tensor, Test_fixed_10_CO2_in_Box_inter_only_VDW)
{
  double tolerance = 1e-4;

  ForceField forceField = TestFactories::makeDefaultFF(12.0, true, false, true);
  forceField.useCharge = false;
  forceField.omitEwaldFourier = true;
  Component c = TestFactories::makeCO2(forceField, 0, true);
  System system = System(0, forceField, SimulationBox(24.0, 24.0, 24.0), 300.0, 1e4, 1.0, {}, {c}, {10}, 5);

  std::span<Atom> atomPositions = system.spanOfMoleculeAtoms();

  atomPositions[0].position = double3(7.988074407572, 14.899584879375, 4.399987406643);
  atomPositions[1].position = double3(7.274874754254, 15.647465005356, 4.902205060325);
  atomPositions[2].position = double3(6.561675100936, 16.395345131337, 5.404422714007);
  atomPositions[3].position = double3(4.827191364635, 11.116729437613, 2.228559049200);
  atomPositions[4].position = double3(4.225319791126, 12.059757497190, 2.490569903656);
  atomPositions[5].position = double3(3.623448217617, 13.002785556767, 2.752580758111);
  atomPositions[6].position = double3(11.213207717903, 1.543908604534, 14.879720523734);
  atomPositions[7].position = double3(10.838296611500, 0.670996368314, 14.233448993994);
  atomPositions[8].position = double3(10.463385505097, -0.201915867907, 13.587177464255);
  atomPositions[9].position = double3(20.136034629402, 13.425705563836, 20.357126733679);
  atomPositions[10].position = double3(19.268912632143, 14.174691153792, 20.271563728193);
  atomPositions[11].position = double3(18.401790634884, 14.923676743748, 20.186000722707);
  atomPositions[12].position = double3(10.293429944185, 2.604978509902, 19.845212277981);
  atomPositions[13].position = double3(9.506941080551, 3.298537002688, 19.375517819354);
  atomPositions[14].position = double3(8.720452216918, 3.992095495474, 18.905823360727);
  atomPositions[15].position = double3(22.681356636431, 22.578271833871, 1.171765550468);
  atomPositions[16].position = double3(23.123101563696, 23.590426461847, 1.488949138160);
  atomPositions[17].position = double3(23.564846490961, 24.602581089824, 1.806132725853);
  atomPositions[18].position = double3(0.923319514763, 2.518001235362, 8.869299391713);
  atomPositions[19].position = double3(1.548565874314, 3.016318526023, 8.044103737521);
  atomPositions[20].position = double3(2.173812233865, 3.514635816685, 7.218908083330);
  atomPositions[21].position = double3(1.269009184104, 22.131774756102, 3.695331679629);
  atomPositions[22].position = double3(0.304888752117, 21.565726817790, 3.960402470471);
  atomPositions[23].position = double3(-0.659231679870, 20.999678879477, 4.225473261313);
  atomPositions[24].position = double3(13.843961437995, 22.253075838217, 21.744383226096);
  atomPositions[25].position = double3(12.854393374997, 21.899980534703, 22.209442099068);
  atomPositions[26].position = double3(11.864825311999, 21.546885231188, 22.674500972040);
  atomPositions[27].position = double3(1.635112578995, 18.974875387109, 5.600705491777);
  atomPositions[28].position = double3(1.273978603709, 19.816503820444, 6.294567749071);
  atomPositions[29].position = double3(0.912844628423, 20.658132253779, 6.988430006364);

  for (Atom& atom : atomPositions)
  {
    atom.gradient = double3(0.0, 0.0, 0.0);
  }

  RunningEnergy energy =
      Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomPositions);

  EXPECT_NEAR(energy.moleculeMoleculeVDW * Units::EnergyToKelvin, -1176.7918606686, tolerance);

  std::pair<EnergyStatus, double3x3> pressureInfo = system.computeMolecularPressure();
  pressureInfo.first.sumTotal();

  EXPECT_NEAR(pressureInfo.first.totalEnergy.energy * Units::EnergyToKelvin, -1176.7918606686, tolerance);
  EXPECT_NEAR(pressureInfo.second.ax, -267.419428, tolerance);
  EXPECT_NEAR(pressureInfo.second.ay, 118.444972, tolerance);
  EXPECT_NEAR(pressureInfo.second.az, -151.475714, tolerance);
  EXPECT_NEAR(pressureInfo.second.bx, 118.444972, tolerance);
  EXPECT_NEAR(pressureInfo.second.by, -846.854115, tolerance);
  EXPECT_NEAR(pressureInfo.second.bz, 228.464316, tolerance);
  EXPECT_NEAR(pressureInfo.second.cx, -151.475714, tolerance);
  EXPECT_NEAR(pressureInfo.second.cy, 228.464316, tolerance);
  EXPECT_NEAR(pressureInfo.second.cz, -684.541460, tolerance);
}

TEST(MC_strain_tensor, Test_fixed_10_CO2_in_Box_inter_no_fourier)
{
  double tolerance = 1e-4;

  ForceField forceField = TestFactories::makeDefaultFF(12.0, true, false, true);
  forceField.omitEwaldFourier = true;
  forceField.EwaldAlpha = 0.265058;
  forceField.numberOfWaveVectors = int3(7, 7, 7);
  Component c = TestFactories::makeCO2(forceField, 0, true);
  System system = System(0, forceField, SimulationBox(24.0, 24.0, 24.0), 300.0, 1e4, 1.0, {}, {c}, {10}, 5);

  std::span<Atom> atomPositions = system.spanOfMoleculeAtoms();

  atomPositions[0].position = double3(7.988074407572, 14.899584879375, 4.399987406643);
  atomPositions[1].position = double3(7.274874754254, 15.647465005356, 4.902205060325);
  atomPositions[2].position = double3(6.561675100936, 16.395345131337, 5.404422714007);
  atomPositions[3].position = double3(4.827191364635, 11.116729437613, 2.228559049200);
  atomPositions[4].position = double3(4.225319791126, 12.059757497190, 2.490569903656);
  atomPositions[5].position = double3(3.623448217617, 13.002785556767, 2.752580758111);
  atomPositions[6].position = double3(11.213207717903, 1.543908604534, 14.879720523734);
  atomPositions[7].position = double3(10.838296611500, 0.670996368314, 14.233448993994);
  atomPositions[8].position = double3(10.463385505097, -0.201915867907, 13.587177464255);
  atomPositions[9].position = double3(20.136034629402, 13.425705563836, 20.357126733679);
  atomPositions[10].position = double3(19.268912632143, 14.174691153792, 20.271563728193);
  atomPositions[11].position = double3(18.401790634884, 14.923676743748, 20.186000722707);
  atomPositions[12].position = double3(10.293429944185, 2.604978509902, 19.845212277981);
  atomPositions[13].position = double3(9.506941080551, 3.298537002688, 19.375517819354);
  atomPositions[14].position = double3(8.720452216918, 3.992095495474, 18.905823360727);
  atomPositions[15].position = double3(22.681356636431, 22.578271833871, 1.171765550468);
  atomPositions[16].position = double3(23.123101563696, 23.590426461847, 1.488949138160);
  atomPositions[17].position = double3(23.564846490961, 24.602581089824, 1.806132725853);
  atomPositions[18].position = double3(0.923319514763, 2.518001235362, 8.869299391713);
  atomPositions[19].position = double3(1.548565874314, 3.016318526023, 8.044103737521);
  atomPositions[20].position = double3(2.173812233865, 3.514635816685, 7.218908083330);
  atomPositions[21].position = double3(1.269009184104, 22.131774756102, 3.695331679629);
  atomPositions[22].position = double3(0.304888752117, 21.565726817790, 3.960402470471);
  atomPositions[23].position = double3(-0.659231679870, 20.999678879477, 4.225473261313);
  atomPositions[24].position = double3(13.843961437995, 22.253075838217, 21.744383226096);
  atomPositions[25].position = double3(12.854393374997, 21.899980534703, 22.209442099068);
  atomPositions[26].position = double3(11.864825311999, 21.546885231188, 22.674500972040);
  atomPositions[27].position = double3(1.635112578995, 18.974875387109, 5.600705491777);
  atomPositions[28].position = double3(1.273978603709, 19.816503820444, 6.294567749071);
  atomPositions[29].position = double3(0.912844628423, 20.658132253779, 6.988430006364);

  for (Atom& atom : atomPositions)
  {
    atom.gradient = double3(0.0, 0.0, 0.0);
  }

  RunningEnergy energy =
      Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomPositions);

  EXPECT_NEAR(energy.moleculeMoleculeVDW * Units::EnergyToKelvin, -1176.7918606686, tolerance);
  EXPECT_NEAR(energy.moleculeMoleculeCharge * Units::EnergyToKelvin, 355.02139790, tolerance);

  std::pair<EnergyStatus, double3x3> pressureInfo = system.computeMolecularPressure();
  pressureInfo.first.sumTotal();

  EXPECT_NEAR(pressureInfo.first.totalEnergy.energy * Units::EnergyToKelvin, -1176.7918606686 + 355.02139790,
              tolerance);
  EXPECT_NEAR(pressureInfo.second.ax, -338.451349, tolerance);
  EXPECT_NEAR(pressureInfo.second.ay, -78.549976, tolerance);
  EXPECT_NEAR(pressureInfo.second.az, -27.238158, tolerance);
  EXPECT_NEAR(pressureInfo.second.bx, -78.549976, tolerance);
  EXPECT_NEAR(pressureInfo.second.by, 38.059080, tolerance);
  EXPECT_NEAR(pressureInfo.second.bz, -532.756964, tolerance);
  EXPECT_NEAR(pressureInfo.second.cx, -27.238158, tolerance);
  EXPECT_NEAR(pressureInfo.second.cy, -532.756964, tolerance);
  EXPECT_NEAR(pressureInfo.second.cz, 20.703013, tolerance);
}

TEST(MC_strain_tensor, Test_fixed_10_CO2_in_Box_inter_ewald)
{
  double tolerance = 1e-4;

  ForceField forceField = TestFactories::makeDefaultFF(12.0, true, false, true);
  forceField.useCharge = true;
  forceField.omitEwaldFourier = false;
  forceField.automaticEwald = false;
  forceField.EwaldAlpha = 0.265058;
  forceField.numberOfWaveVectors = int3(7, 7, 7);

  Component c = TestFactories::makeCO2(forceField, 0, true);
  System system = System(0, forceField, SimulationBox(24.0, 24.0, 24.0), 300.0, 1e4, 1.0, {}, {c}, {10}, 5);

  std::span<Atom> atomPositions = system.spanOfMoleculeAtoms();

  atomPositions[0].position = double3(7.98807441, 14.89958488, 4.39998741);
  atomPositions[1].position = double3(7.27487475, 15.64746501, 4.90220506);
  atomPositions[2].position = double3(6.56167510, 16.39534513, 5.40442271);
  atomPositions[3].position = double3(4.82719136, 11.11672944, 2.22855905);
  atomPositions[4].position = double3(4.22531979, 12.05975750, 2.49056990);
  atomPositions[5].position = double3(3.62344822, 13.00278556, 2.75258076);
  atomPositions[6].position = double3(11.21320772, 1.54390860, 14.87972052);
  atomPositions[7].position = double3(10.83829661, 0.67099637, 14.23344899);
  atomPositions[8].position = double3(10.46338551, -0.20191587, 13.58717746);
  atomPositions[9].position = double3(20.13603463, 13.42570556, 20.35712673);
  atomPositions[10].position = double3(19.26891263, 14.17469115, 20.27156373);
  atomPositions[11].position = double3(18.40179063, 14.92367674, 20.18600072);
  atomPositions[12].position = double3(10.29342994, 2.60497851, 19.84521228);
  atomPositions[13].position = double3(9.50694108, 3.29853700, 19.37551782);
  atomPositions[14].position = double3(8.72045222, 3.99209550, 18.90582336);
  atomPositions[15].position = double3(22.68135664, 22.57827183, 1.17176555);
  atomPositions[16].position = double3(23.12310156, 23.59042646, 1.48894914);
  atomPositions[17].position = double3(23.56484649, 24.60258109, 1.80613273);
  atomPositions[18].position = double3(0.92331951, 2.51800124, 8.86929939);
  atomPositions[19].position = double3(1.54856587, 3.01631853, 8.04410374);
  atomPositions[20].position = double3(2.17381223, 3.51463582, 7.21890808);
  atomPositions[21].position = double3(1.26900918, 22.13177476, 3.69533168);
  atomPositions[22].position = double3(0.30488875, 21.56572682, 3.96040247);
  atomPositions[23].position = double3(-0.65923168, 20.99967888, 4.22547326);
  atomPositions[24].position = double3(13.84396144, 22.25307584, 21.74438323);
  atomPositions[25].position = double3(12.85439337, 21.89998053, 22.20944210);
  atomPositions[26].position = double3(11.86482531, 21.54688523, 22.67450097);
  atomPositions[27].position = double3(1.63511258, 18.97487539, 5.60070549);
  atomPositions[28].position = double3(1.27397860, 19.81650382, 6.29456775);
  atomPositions[29].position = double3(0.91284463, 20.65813225, 6.98843001);

  for (Atom& atom : atomPositions)
  {
    atom.gradient = double3(0.0, 0.0, 0.0);
  }

  system.precomputeTotalRigidEnergy();
  system.CoulombicFourierEnergySingleIon = Interactions::computeEwaldFourierEnergySingleIon(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.forceField, system.simulationBox,
      double3(0.0, 0.0, 0.0), 1.0);
  RunningEnergy energy = Interactions::computeEwaldFourierEnergy(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
      system.forceField, system.simulationBox, system.components, system.numberOfMoleculesPerComponent,
      system.spanOfMoleculeAtoms());

  EXPECT_NEAR(energy.ewaldFourier() * Units::EnergyToKelvin, 47.48747965, tolerance);

  std::pair<EnergyStatus, double3x3> pressureInfo = system.computeMolecularPressure();
  pressureInfo.first.sumTotal();

  EXPECT_NEAR(pressureInfo.first.totalEnergy.energy * Units::EnergyToKelvin,
              -1176.7918606686 + 355.02165559 + 47.48747965, tolerance);
  EXPECT_NEAR(pressureInfo.second.ax, -404.460729, tolerance);
  EXPECT_NEAR(pressureInfo.second.ay, -24.018427, tolerance);
  EXPECT_NEAR(pressureInfo.second.az, 7.832134, tolerance);
  EXPECT_NEAR(pressureInfo.second.bx, -24.018427, tolerance);
  EXPECT_NEAR(pressureInfo.second.by, 76.925477, tolerance);
  EXPECT_NEAR(pressureInfo.second.bz, -534.677553, tolerance);
  EXPECT_NEAR(pressureInfo.second.cx, 7.832134, tolerance);
  EXPECT_NEAR(pressureInfo.second.cy, -534.677553, tolerance);
  EXPECT_NEAR(pressureInfo.second.cz, 48.977581, tolerance);
}

TEST(MC_strain_tensor, Test_20_CH4_25x25x25_LJ)
{
  double delta = 1e-7;
  double tolerance = 1e-4;

  ForceField forceField = TestFactories::makeDefaultFF(12.0, true, false, false);
  Component c = TestFactories::makeMethane(forceField, 0);
  System system = System(0, forceField, SimulationBox(25.0, 25.0, 25.0), 300.0, 1e4, 1.0, {}, {c}, {20}, 5);

  std::span<Atom> moleculeAtomPositions = system.spanOfMoleculeAtoms();

  for (Atom& atom : moleculeAtomPositions)
  {
    atom.gradient = double3(0.0, 0.0, 0.0);
  }

  std::pair<EnergyStatus, double3x3> pressureInfo = Interactions::computeInterMolecularEnergyStrainDerivative(
      system.forceField, system.components, system.simulationBox, moleculeAtomPositions);

  std::vector<std::pair<double3x3, double>> strains{
      std::pair{double3x3{double3{delta, 0.0, 0.0}, double3{0.0, 0.0, 0.0}, double3{0.0, 0.0, 0.0}},
                pressureInfo.second.ax},
      std::pair{double3x3{double3{0.0, delta, 0.0}, double3{0.0, 0.0, 0.0}, double3{0.0, 0.0, 0.0}},
                pressureInfo.second.bx},
      std::pair{double3x3{double3{0.0, 0.0, delta}, double3{0.0, 0.0, 0.0}, double3{0.0, 0.0, 0.0}},
                pressureInfo.second.cx},
      std::pair{double3x3{double3{0.0, 0.0, 0.0}, double3{delta, 0.0, 0.0}, double3{0.0, 0.0, 0.0}},
                pressureInfo.second.ay},
      std::pair{double3x3{double3{0.0, 0.0, 0.0}, double3{0.0, delta, 0.0}, double3{0.0, 0.0, 0.0}},
                pressureInfo.second.by},
      std::pair{double3x3{double3{0.0, 0.0, 0.0}, double3{0.0, 0.0, delta}, double3{0.0, 0.0, 0.0}},
                pressureInfo.second.cy},
      std::pair{double3x3{double3{0.0, 0.0, 0.0}, double3{0.0, 0.0, 0.0}, double3{delta, 0.0, 0.0}},
                pressureInfo.second.az},
      std::pair{double3x3{double3{0.0, 0.0, 0.0}, double3{0.0, 0.0, 0.0}, double3{0.0, delta, 0.0}},
                pressureInfo.second.bz},
      std::pair{double3x3{double3{0.0, 0.0, 0.0}, double3{0.0, 0.0, 0.0}, double3{0.0, 0.0, delta}},
                pressureInfo.second.cz}};

  double3x3 inv = system.simulationBox.inverseCell;
  double3x3 identity{double3{1.0, 0.0, 0.0}, double3{0.0, 1.0, 0.0}, double3{0.0, 0.0, 1.0}};

  for (const std::pair<double3x3, double> strain : strains)
  {
    SimulationBox strainBox_forward2 =
        SimulationBox((identity + strain.first) * system.simulationBox.cell, SimulationBox::Type::Triclinic);
    std::vector<Atom> moleculeAtomPositions_forward2{};
    std::transform(moleculeAtomPositions.begin(), moleculeAtomPositions.end(),
                   std::back_inserter(moleculeAtomPositions_forward2),
                   [&strainBox_forward2, &inv](const Atom& m)
                   {
                     return Atom(strainBox_forward2.cell * (inv * m.position), m.charge, 1.0, m.moleculeId, m.type,
                                 m.componentId, m.groupId, m.isFractional);
                   });
    RunningEnergy EnergyForward2 = Interactions::computeInterMolecularEnergy(system.forceField, strainBox_forward2,
                                                                             moleculeAtomPositions_forward2);

    SimulationBox strainBox_forward1 =
        SimulationBox((identity + 0.5 * strain.first) * system.simulationBox.cell, SimulationBox::Type::Triclinic);
    std::vector<Atom> moleculeAtomPositions_forward1{};
    std::transform(moleculeAtomPositions.begin(), moleculeAtomPositions.end(),
                   std::back_inserter(moleculeAtomPositions_forward1),
                   [&strainBox_forward1, &inv](const Atom& m)
                   {
                     return Atom(strainBox_forward1.cell * (inv * m.position), m.charge, 1.0, m.moleculeId, m.type,
                                 m.componentId, m.groupId, m.isFractional);
                   });
    RunningEnergy EnergyForward1 = Interactions::computeInterMolecularEnergy(system.forceField, strainBox_forward1,
                                                                             moleculeAtomPositions_forward1);

    SimulationBox strainBox_backward1 =
        SimulationBox((identity - 0.5 * strain.first) * system.simulationBox.cell, SimulationBox::Type::Triclinic);
    std::vector<Atom> moleculeAtomPositions_backward1{};
    std::transform(moleculeAtomPositions.begin(), moleculeAtomPositions.end(),
                   std::back_inserter(moleculeAtomPositions_backward1),
                   [&strainBox_backward1, &inv](const Atom& m)
                   {
                     return Atom(strainBox_backward1.cell * (inv * m.position), m.charge, 1.0, m.moleculeId, m.type,
                                 m.componentId, m.groupId, m.isFractional);
                   });
    RunningEnergy EnergyBackward1 = Interactions::computeInterMolecularEnergy(system.forceField, strainBox_backward1,
                                                                              moleculeAtomPositions_backward1);

    SimulationBox strainBox_backward2 =
        SimulationBox((identity - strain.first) * system.simulationBox.cell, SimulationBox::Type::Triclinic);
    std::vector<Atom> moleculeAtomPositions_backward2{};
    std::transform(moleculeAtomPositions.begin(), moleculeAtomPositions.end(),
                   std::back_inserter(moleculeAtomPositions_backward2),
                   [&strainBox_backward2, &inv](const Atom& m)
                   {
                     return Atom(strainBox_backward2.cell * (inv * m.position), m.charge, 1.0, m.moleculeId, m.type,
                                 m.componentId, m.groupId, m.isFractional);
                   });
    RunningEnergy EnergyBackward2 = Interactions::computeInterMolecularEnergy(system.forceField, strainBox_backward2,
                                                                              moleculeAtomPositions_backward2);

    double strainDerivative = (-EnergyForward2.potentialEnergy() + 8.0 * EnergyForward1.potentialEnergy() -
                               8.0 * EnergyBackward1.potentialEnergy() + EnergyBackward2.potentialEnergy()) /
                              (6.0 * delta);

    EXPECT_NEAR(strainDerivative, strain.second, tolerance) << "Wrong strainDerivative";
  }
}

TEST(MC_strain_tensor, Test_20_Na_Cl_25x25x25_LJ_Real)
{
  double delta = 1e-7;
  double tolerance = 1e-4;

  ForceField forceField = TestFactories::makeDefaultFF(12.0, true, false, true);
  Component na = TestFactories::makeIon(forceField, 0, "Na", 6, 0.0);
  Component cl = TestFactories::makeIon(forceField, 1, "Cl", 7, 0.0);
  System system = System(0, forceField, SimulationBox(25.0, 25.0, 25.0), 300.0, 1e4, 1.0, {}, {na, cl}, {1, 1}, 5);

  std::span<Atom> moleculeAtomPositions = system.spanOfMoleculeAtoms();

  for (Atom& atom : moleculeAtomPositions)
  {
    atom.gradient = double3(0.0, 0.0, 0.0);
  }

  system.forceField.EwaldAlpha = 0.25;
  system.forceField.numberOfWaveVectors = int3(8, 8, 8);

  std::pair<EnergyStatus, double3x3> pressureInfo = Interactions::computeInterMolecularEnergyStrainDerivative(
      system.forceField, system.components, system.simulationBox, moleculeAtomPositions);
  pressureInfo.first.sumTotal();

  std::vector<std::pair<double3x3, double>> strains{
      std::pair{double3x3{double3{delta, 0.0, 0.0}, double3{0.0, 0.0, 0.0}, double3{0.0, 0.0, 0.0}},
                pressureInfo.second.ax},
      std::pair{double3x3{double3{0.0, delta, 0.0}, double3{0.0, 0.0, 0.0}, double3{0.0, 0.0, 0.0}},
                pressureInfo.second.bx},
      std::pair{double3x3{double3{0.0, 0.0, delta}, double3{0.0, 0.0, 0.0}, double3{0.0, 0.0, 0.0}},
                pressureInfo.second.cx},
      std::pair{double3x3{double3{0.0, 0.0, 0.0}, double3{delta, 0.0, 0.0}, double3{0.0, 0.0, 0.0}},
                pressureInfo.second.ay},
      std::pair{double3x3{double3{0.0, 0.0, 0.0}, double3{0.0, delta, 0.0}, double3{0.0, 0.0, 0.0}},
                pressureInfo.second.by},
      std::pair{double3x3{double3{0.0, 0.0, 0.0}, double3{0.0, 0.0, delta}, double3{0.0, 0.0, 0.0}},
                pressureInfo.second.cy},
      std::pair{double3x3{double3{0.0, 0.0, 0.0}, double3{0.0, 0.0, 0.0}, double3{delta, 0.0, 0.0}},
                pressureInfo.second.az},
      std::pair{double3x3{double3{0.0, 0.0, 0.0}, double3{0.0, 0.0, 0.0}, double3{0.0, delta, 0.0}},
                pressureInfo.second.bz},
      std::pair{double3x3{double3{0.0, 0.0, 0.0}, double3{0.0, 0.0, 0.0}, double3{0.0, 0.0, delta}},
                pressureInfo.second.cz}};

  double3x3 inv = system.simulationBox.inverseCell;
  double3x3 identity{double3{1.0, 0.0, 0.0}, double3{0.0, 1.0, 0.0}, double3{0.0, 0.0, 1.0}};

  for (const std::pair<double3x3, double> strain : strains)
  {
    SimulationBox strainBox_forward2 =
        SimulationBox((identity + strain.first) * system.simulationBox.cell, SimulationBox::Type::Triclinic);
    std::vector<Atom> moleculeAtomPositions_forward2{};
    std::transform(moleculeAtomPositions.begin(), moleculeAtomPositions.end(),
                   std::back_inserter(moleculeAtomPositions_forward2),
                   [&strainBox_forward2, &inv](const Atom& m)
                   {
                     return Atom(strainBox_forward2.cell * (inv * m.position), m.charge, 1.0, m.moleculeId, m.type,
                                 m.componentId, m.groupId, m.isFractional);
                   });
    RunningEnergy EnergyForward2 = Interactions::computeInterMolecularEnergy(system.forceField, strainBox_forward2,
                                                                             moleculeAtomPositions_forward2);

    SimulationBox strainBox_forward1 =
        SimulationBox((identity + 0.5 * strain.first) * system.simulationBox.cell, SimulationBox::Type::Triclinic);
    std::vector<Atom> moleculeAtomPositions_forward1{};
    std::transform(moleculeAtomPositions.begin(), moleculeAtomPositions.end(),
                   std::back_inserter(moleculeAtomPositions_forward1),
                   [&strainBox_forward1, &inv](const Atom& m)
                   {
                     return Atom(strainBox_forward1.cell * (inv * m.position), m.charge, 1.0, m.moleculeId, m.type,
                                 m.componentId, m.groupId, m.isFractional);
                   });
    RunningEnergy EnergyForward1 = Interactions::computeInterMolecularEnergy(system.forceField, strainBox_forward1,
                                                                             moleculeAtomPositions_forward1);

    SimulationBox strainBox_backward1 =
        SimulationBox((identity - 0.5 * strain.first) * system.simulationBox.cell, SimulationBox::Type::Triclinic);
    std::vector<Atom> moleculeAtomPositions_backward1{};
    std::transform(moleculeAtomPositions.begin(), moleculeAtomPositions.end(),
                   std::back_inserter(moleculeAtomPositions_backward1),
                   [&strainBox_backward1, &inv](const Atom& m)
                   {
                     return Atom(strainBox_backward1.cell * (inv * m.position), m.charge, 1.0, m.moleculeId, m.type,
                                 m.componentId, m.groupId, m.isFractional);
                   });
    RunningEnergy EnergyBackward1 = Interactions::computeInterMolecularEnergy(system.forceField, strainBox_backward1,
                                                                              moleculeAtomPositions_backward1);

    SimulationBox strainBox_backward2 =
        SimulationBox((identity - strain.first) * system.simulationBox.cell, SimulationBox::Type::Triclinic);
    std::vector<Atom> moleculeAtomPositions_backward2{};
    std::transform(moleculeAtomPositions.begin(), moleculeAtomPositions.end(),
                   std::back_inserter(moleculeAtomPositions_backward2),
                   [&strainBox_backward2, &inv](const Atom& m)
                   {
                     return Atom(strainBox_backward2.cell * (inv * m.position), m.charge, 1.0, m.moleculeId, m.type,
                                 m.componentId, m.groupId, m.isFractional);
                   });
    RunningEnergy EnergyBackward2 = Interactions::computeInterMolecularEnergy(system.forceField, strainBox_backward2,
                                                                              moleculeAtomPositions_backward2);

    double strainDerivative = (-EnergyForward2.potentialEnergy() + 8.0 * EnergyForward1.potentialEnergy() -
                               8.0 * EnergyBackward1.potentialEnergy() + EnergyBackward2.potentialEnergy()) /
                              (6.0 * delta);

    EXPECT_NEAR(strainDerivative, strain.second, tolerance) << "Wrong strainDerivative";
  }
}

TEST(MC_strain_tensor, Test_20_Na_Cl_in_Box_25x25x25_strain_derivative)
{
  double delta = 1e-4;
  double tolerance = 1e-3;

  ForceField forceField = TestFactories::makeDefaultFF(12.0, true, false, true);
  Component na = TestFactories::makeIon(forceField, 0, "Na", 6, 0.0);
  Component cl = TestFactories::makeIon(forceField, 1, "Cl", 7, 0.0);

  forceField.EwaldAlpha = 0.25;
  forceField.numberOfWaveVectors = int3(8, 8, 8);

  System system = System(0, forceField, SimulationBox(25.0, 25.0, 25.0), 300.0, 1e4, 1.0, {}, {na, cl}, {20, 20}, 5);

  std::span<Atom> moleculeAtomPositions = system.spanOfMoleculeAtoms();

  for (Atom& atom : moleculeAtomPositions)
  {
    atom.gradient = double3(0.0, 0.0, 0.0);
  }

  for (size_t i = 0; i < 20; ++i)
  {
    moleculeAtomPositions[i].charge = 1.0;
    system.atomPositions[i].charge = 1.0;
  }
  for (size_t i = 0; i < 20; ++i)
  {
    moleculeAtomPositions[i + 20].charge = -1.0;
    system.atomPositions[i + 20].charge = -1.0;
  }

  for (Atom& atom : moleculeAtomPositions)
  {
    atom.gradient = double3(0.0, 0.0, 0.0);
  }

  RunningEnergy energy, rigidenergy;

  system.precomputeTotalRigidEnergy();
  std::pair<EnergyStatus, double3x3> pressureInfo = Interactions::computeEwaldFourierEnergyStrainDerivative(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
      system.forceField, system.simulationBox, system.framework, system.components,
      system.numberOfMoleculesPerComponent, system.atomPositions, system.CoulombicFourierEnergySingleIon,
      system.netChargeFramework, system.netChargePerComponent);
  pressureInfo.first.sumTotal();

  std::vector<std::pair<double3x3, double>> strains{
      std::pair{double3x3{double3{delta, 0.0, 0.0}, double3{0.0, 0.0, 0.0}, double3{0.0, 0.0, 0.0}},
                pressureInfo.second.ax},
      std::pair{double3x3{double3{0.0, delta, 0.0}, double3{0.0, 0.0, 0.0}, double3{0.0, 0.0, 0.0}},
                pressureInfo.second.bx},
      std::pair{double3x3{double3{0.0, 0.0, delta}, double3{0.0, 0.0, 0.0}, double3{0.0, 0.0, 0.0}},
                pressureInfo.second.cx},
      std::pair{double3x3{double3{0.0, 0.0, 0.0}, double3{delta, 0.0, 0.0}, double3{0.0, 0.0, 0.0}},
                pressureInfo.second.ay},
      std::pair{double3x3{double3{0.0, 0.0, 0.0}, double3{0.0, delta, 0.0}, double3{0.0, 0.0, 0.0}},
                pressureInfo.second.by},
      std::pair{double3x3{double3{0.0, 0.0, 0.0}, double3{0.0, 0.0, delta}, double3{0.0, 0.0, 0.0}},
                pressureInfo.second.cy},
      std::pair{double3x3{double3{0.0, 0.0, 0.0}, double3{0.0, 0.0, 0.0}, double3{delta, 0.0, 0.0}},
                pressureInfo.second.az},
      std::pair{double3x3{double3{0.0, 0.0, 0.0}, double3{0.0, 0.0, 0.0}, double3{0.0, delta, 0.0}},
                pressureInfo.second.bz},
      std::pair{double3x3{double3{0.0, 0.0, 0.0}, double3{0.0, 0.0, 0.0}, double3{0.0, 0.0, delta}},
                pressureInfo.second.cz}};

  double3x3 inv = system.simulationBox.inverseCell;
  double3x3 identity{double3{1.0, 0.0, 0.0}, double3{0.0, 1.0, 0.0}, double3{0.0, 0.0, 1.0}};

  for (const std::pair<double3x3, double> strain : strains)
  {
    SimulationBox strainBox_forward2 =
        SimulationBox((identity + strain.first) * system.simulationBox.cell, SimulationBox::Type::Triclinic);
    std::vector<Atom> moleculeAtomPositions_forward2{};
    std::transform(moleculeAtomPositions.begin(), moleculeAtomPositions.end(),
                   std::back_inserter(moleculeAtomPositions_forward2),
                   [&strainBox_forward2, &inv](const Atom& m)
                   {
                     return Atom(strainBox_forward2.cell * (inv * m.position), m.charge, 1.0, m.moleculeId, m.type,
                                 m.componentId, m.groupId, m.isFractional);
                   });

    RunningEnergy EnergyForward2 = Interactions::computeEwaldFourierEnergy(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
        system.forceField, strainBox_forward2, system.components, system.numberOfMoleculesPerComponent,
        moleculeAtomPositions_forward2);

    SimulationBox strainBox_forward1 =
        SimulationBox((identity + 0.5 * strain.first) * system.simulationBox.cell, SimulationBox::Type::Triclinic);
    std::vector<Atom> moleculeAtomPositions_forward1{};
    std::transform(moleculeAtomPositions.begin(), moleculeAtomPositions.end(),
                   std::back_inserter(moleculeAtomPositions_forward1),
                   [&strainBox_forward1, &inv](const Atom& m)
                   {
                     return Atom(strainBox_forward1.cell * (inv * m.position), m.charge, 1.0, m.moleculeId, m.type,
                                 m.componentId, m.groupId, m.isFractional);
                   });

    RunningEnergy EnergyForward1 = Interactions::computeEwaldFourierEnergy(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
        system.forceField, strainBox_forward1, system.components, system.numberOfMoleculesPerComponent,
        moleculeAtomPositions_forward1);

    SimulationBox strainBox_backward1 =
        SimulationBox((identity - 0.5 * strain.first) * system.simulationBox.cell, SimulationBox::Type::Triclinic);
    std::vector<Atom> moleculeAtomPositions_backward1{};
    std::transform(moleculeAtomPositions.begin(), moleculeAtomPositions.end(),
                   std::back_inserter(moleculeAtomPositions_backward1),
                   [&strainBox_backward1, &inv](const Atom& m)
                   {
                     return Atom(strainBox_backward1.cell * (inv * m.position), m.charge, 1.0, m.moleculeId, m.type,
                                 m.componentId, m.groupId, m.isFractional);
                   });

    RunningEnergy EnergyBackward1 = Interactions::computeEwaldFourierEnergy(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
        system.forceField, strainBox_backward1, system.components, system.numberOfMoleculesPerComponent,
        moleculeAtomPositions_backward1);

    SimulationBox strainBox_backward2 =
        SimulationBox((identity - strain.first) * system.simulationBox.cell, SimulationBox::Type::Triclinic);
    std::vector<Atom> moleculeAtomPositions_backward2{};
    std::transform(moleculeAtomPositions.begin(), moleculeAtomPositions.end(),
                   std::back_inserter(moleculeAtomPositions_backward2),
                   [&strainBox_backward2, &inv](const Atom& m)
                   {
                     return Atom(strainBox_backward2.cell * (inv * m.position), m.charge, 1.0, m.moleculeId, m.type,
                                 m.componentId, m.groupId, m.isFractional);
                   });

    RunningEnergy EnergyBackward2 = Interactions::computeEwaldFourierEnergy(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
        system.forceField, strainBox_backward2, system.components, system.numberOfMoleculesPerComponent,
        moleculeAtomPositions_backward2);

    double strainDerivativeApproximation =
        (-EnergyForward2.potentialEnergy() + 8.0 * EnergyForward1.potentialEnergy() -
         8.0 * EnergyBackward1.potentialEnergy() + EnergyBackward2.potentialEnergy()) /
        (6.0 * delta);

    EXPECT_NEAR(strainDerivativeApproximation, strain.second, tolerance) << "Wrong strainDerivative";
  }
}

// TODO
/*
TEST(MC_strain_tensor, Test_10_CO2_in_Box_inter)
{
  double delta = 1e-7;
  double tolerance = 1e-4;

  ForceField forceField = ForceField(
      {
          PseudoAtom("Si", true, 28.0855, 2.05, 0.0, 14, false),
          PseudoAtom("O", true, 15.999, -1.025, 0.0, 8, false),
          PseudoAtom("CH4", false, 16.04246, 0.0, 0.0, 6, false),
          PseudoAtom("C_co2", false, 12.0, 0.6512, 0.0, 6, false),
          PseudoAtom("O_co2", false, 15.9994, -0.3256, 0.0, 8, false),
      },
      {VDWParameters(22.0, 2.30),
       VDWParameters(53.0, 3.3),
       VDWParameters(158.5, 3.72),
       VDWParameters(29.933, 2.745),
       VDWParameters(85.671, 3.017)},
      ForceField::MixingRule::Lorentz_Berthelot, 11.8, 11.8, 11.8, true, false, true);

  Component c = Component(
      0, forceField, "CO2", 304.1282, 7377300.0, 0.22394,
      {// double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type, uint8_t componentId,
       // uint8_t groupId
       Atom(double3(0.0, 0.0, 1.149), 0.0, 1.0, 0, 4, 1, 0),
       Atom(double3(0.0, 0.0, 0.0), 0.0, 1.0, 0, 3, 1, 0),
       Atom(double3(0.0, 0.0, -1.149), 0.0, 1.0, 0, 4, 1, 0)},
      5, 21);

  System system = System(0, forceField, SimulationBox(25.0, 25.0, 25.0), 300.0, 1e4, 1.0, {}, {c}, {10}, 5);

  std::span<Atom> moleculeAtomPositions = system.spanOfMoleculeAtoms();

  for (Atom& atom : moleculeAtomPositions)
  {
    atom.gradient = double3(0.0, 0.0, 0.0);
  }

  std::pair<EnergyStatus, double3x3> pressureInfo = Interactions::computeInterMolecularEnergyStrainDerivative(
      system.forceField, system.components, system.simulationBox, moleculeAtomPositions);

  std::vector<std::pair<double3x3, double>> strains{
      std::pair{double3x3{double3{delta, 0.0, 0.0}, double3{0.0, 0.0, 0.0}, double3{0.0, 0.0, 0.0}},
                pressureInfo.second.ax},
      std::pair{double3x3{double3{0.0, delta, 0.0}, double3{0.0, 0.0, 0.0}, double3{0.0, 0.0, 0.0}},
                pressureInfo.second.bx},
      std::pair{double3x3{double3{0.0, 0.0, delta}, double3{0.0, 0.0, 0.0}, double3{0.0, 0.0, 0.0}},
                pressureInfo.second.cx},
      std::pair{double3x3{double3{0.0, 0.0, 0.0}, double3{delta, 0.0, 0.0}, double3{0.0, 0.0, 0.0}},
                pressureInfo.second.ay},
      std::pair{double3x3{double3{0.0, 0.0, 0.0}, double3{0.0, delta, 0.0}, double3{0.0, 0.0, 0.0}},
                pressureInfo.second.by},
      std::pair{double3x3{double3{0.0, 0.0, 0.0}, double3{0.0, 0.0, delta}, double3{0.0, 0.0, 0.0}},
                pressureInfo.second.cy},
      std::pair{double3x3{double3{0.0, 0.0, 0.0}, double3{0.0, 0.0, 0.0}, double3{delta, 0.0, 0.0}},
                pressureInfo.second.az},
      std::pair{double3x3{double3{0.0, 0.0, 0.0}, double3{0.0, 0.0, 0.0}, double3{0.0, delta, 0.0}},
                pressureInfo.second.bz},
      std::pair{double3x3{double3{0.0, 0.0, 0.0}, double3{0.0, 0.0, 0.0}, double3{0.0, 0.0, delta}},
                pressureInfo.second.cz}};

  double3x3 inv = system.simulationBox.inverseCell;
  double3x3 identity{double3{1.0, 0.0, 0.0}, double3{0.0, 1.0, 0.0}, double3{0.0, 0.0, 1.0}};

  for (const std::pair<double3x3, double> strain : strains)
  {
    SimulationBox strainBox_forward2 =
        SimulationBox((identity + strain.first) * system.simulationBox.cell, SimulationBox::Type::Triclinic);
    std::vector<Atom> moleculeAtomPositions_forward2{};
    std::transform(moleculeAtomPositions.begin(), moleculeAtomPositions.end(),
                   std::back_inserter(moleculeAtomPositions_forward2),
                   [&strainBox_forward2, &inv](const Atom& m)
                   {
                     return Atom(strainBox_forward2.cell * (inv * m.position), m.charge, 1.0, m.moleculeId, m.type,
                                 m.componentId, m.groupId);
                   });
    RunningEnergy EnergyForward2 = Interactions::computeInterMolecularEnergy(system.forceField, strainBox_forward2,
                                                                             moleculeAtomPositions_forward2);

    SimulationBox strainBox_forward1 =
        SimulationBox((identity + 0.5 * strain.first) * system.simulationBox.cell, SimulationBox::Type::Triclinic);
    std::vector<Atom> moleculeAtomPositions_forward1{};
    std::transform(moleculeAtomPositions.begin(), moleculeAtomPositions.end(),
                   std::back_inserter(moleculeAtomPositions_forward1),
                   [&strainBox_forward1, &inv](const Atom& m)
                   {
                     return Atom(strainBox_forward1.cell * (inv * m.position), m.charge, 1.0, m.moleculeId, m.type,
                                 m.componentId, m.groupId);
                   });
    RunningEnergy EnergyForward1 = Interactions::computeInterMolecularEnergy(system.forceField, strainBox_forward1,
                                                                             moleculeAtomPositions_forward1);

    SimulationBox strainBox_backward1 =
        SimulationBox((identity - 0.5 * strain.first) * system.simulationBox.cell, SimulationBox::Type::Triclinic);
    std::vector<Atom> moleculeAtomPositions_backward1{};
    std::transform(moleculeAtomPositions.begin(), moleculeAtomPositions.end(),
                   std::back_inserter(moleculeAtomPositions_backward1),
                   [&strainBox_backward1, &inv](const Atom& m)
                   {
                     return Atom(strainBox_backward1.cell * (inv * m.position), m.charge, 1.0, m.moleculeId, m.type,
                                 m.componentId, m.groupId);
                   });
    RunningEnergy EnergyBackward1 = Interactions::computeInterMolecularEnergy(system.forceField, strainBox_backward1,
                                                                              moleculeAtomPositions_backward1);

    SimulationBox strainBox_backward2 =
        SimulationBox((identity - strain.first) * system.simulationBox.cell, SimulationBox::Type::Triclinic);
    std::vector<Atom> moleculeAtomPositions_backward2{};
    std::transform(moleculeAtomPositions.begin(), moleculeAtomPositions.end(),
                   std::back_inserter(moleculeAtomPositions_backward2),
                   [&strainBox_backward2, &inv](const Atom& m)
                   {
                     return Atom(strainBox_backward2.cell * (inv * m.position), m.charge, 1.0, m.moleculeId, m.type,
                                 m.componentId, m.groupId);
                   });
    RunningEnergy EnergyBackward2 = Interactions::computeInterMolecularEnergy(system.forceField, strainBox_backward2,
                                                                              moleculeAtomPositions_backward2);

    double strainDerivative = (-EnergyForward2.potentialEnergy() + 8.0 * EnergyForward1.potentialEnergy() -
                               8.0 * EnergyBackward1.potentialEnergy() + EnergyBackward2.potentialEnergy()) /
                              (6.0 * delta);

    //EXPECT_NEAR(strainDerivative, strain.second, tolerance) << "Wrong strainDerivative";
  }

}
*/
