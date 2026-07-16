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
import framework;
import component;
import system;
import simulationbox;
import energy_status;
import running_energy;
import randomnumbers;
import interactions_intermolecular;
import interactions_framework_molecule;
import interactions_ewald;
import interactions_polarization_derivatives;
import mc_moves_translation;
import connectivity_table;
import intra_molecular_potentials;

namespace
{
// Ten CO2 molecules (30 atoms) in a 24 A cubic box, positions taken from the existing pressure regression
// tests. Returned as a fully-built System so the intra-molecular exclusion of the finite-cutoff Coulomb
// methods is exercised (each CO2 has three charged sites within the Coulomb cutoff).
System makeTenCO2(ForceField forceField)
{
  Component c = Component::makeCO2(forceField, 0, true);
  System system = System(forceField, SimulationBox(24.0, 24.0, 24.0), false, 300.0, 1e4, 1.0, {}, {c}, {}, {10}, 5);

  std::span<Atom> atomData = system.spanOfMoleculeAtoms();
  const std::array<double3, 30> positions{
      double3(7.988074407572, 14.899584879375, 4.399987406643), double3(7.274874754254, 15.647465005356, 4.902205060325),
      double3(6.561675100936, 16.395345131337, 5.404422714007), double3(4.827191364635, 11.116729437613, 2.228559049200),
      double3(4.225319791126, 12.059757497190, 2.490569903656), double3(3.623448217617, 13.002785556767, 2.752580758111),
      double3(11.213207717903, 1.543908604534, 14.879720523734), double3(10.838296611500, 0.670996368314, 14.233448993994),
      double3(10.463385505097, -0.201915867907, 13.587177464255), double3(20.136034629402, 13.425705563836, 20.357126733679),
      double3(19.268912632143, 14.174691153792, 20.271563728193), double3(18.401790634884, 14.923676743748, 20.186000722707),
      double3(10.293429944185, 2.604978509902, 19.845212277981), double3(9.506941080551, 3.298537002688, 19.375517819354),
      double3(8.720452216918, 3.992095495474, 18.905823360727), double3(22.681356636431, 22.578271833871, 1.171765550468),
      double3(23.123101563696, 23.590426461847, 1.488949138160), double3(23.564846490961, 24.602581089824, 1.806132725853),
      double3(0.923319514763, 2.518001235362, 8.869299391713), double3(1.548565874314, 3.016318526023, 8.044103737521),
      double3(2.173812233865, 3.514635816685, 7.218908083330), double3(1.269009184104, 22.131774756102, 3.695331679629),
      double3(0.304888752117, 21.565726817790, 3.960402470471), double3(-0.659231679870, 20.999678879477, 4.225473261313),
      double3(13.843961437995, 22.253075838217, 21.744383226096), double3(12.854393374997, 21.899980534703, 22.209442099068),
      double3(11.864825311999, 21.546885231188, 22.674500972040), double3(1.635112578995, 18.974875387109, 5.600705491777),
      double3(1.273978603709, 19.816503820444, 6.294567749071), double3(0.912844628423, 20.658132253779, 6.988430006364)};
  for (std::size_t i = 0; i < positions.size(); ++i) atomData[i].position = positions[i];

  return system;
}

// Central 4th-order finite-difference of the (self + intra-molecular exclusion) energy of the finite-cutoff
// Coulomb method under an affine strain of ALL atom positions, matching the convention used by the
// Ewald strain-derivative regression tests.
double strainFiniteDifference(System& system, const double3x3& strainDirection, double delta)
{
  std::span<const Atom> atoms = system.spanOfMoleculeAtoms();
  double3x3 inv = system.simulationBox.inverseCell;
  double3x3 identity{double3{1.0, 0.0, 0.0}, double3{0.0, 1.0, 0.0}, double3{0.0, 0.0, 1.0}};

  auto energyAt = [&](double scale)
  {
    SimulationBox box((identity + scale * strainDirection) * system.simulationBox.cell, SimulationBox::Type::Triclinic);
    std::vector<Atom> strained;
    strained.reserve(atoms.size());
    std::transform(atoms.begin(), atoms.end(), std::back_inserter(strained), [&](const Atom& m)
                   { return Atom(box.cell * (inv * m.position), m.charge, 1.0, m.moleculeId, m.type, m.componentId,
                                 m.groupId, m.isFractional); });
    return Interactions::computeEwaldFourierEnergy(system.eik_x, system.eik_y, system.eik_z, system.eik_xy,
                                                   system.fixedFrameworkStoredEik, system.storedEik, system.forceField,
                                                   box, system.components, system.numberOfMoleculesPerComponent,
                                                   strained)
        .potentialEnergy();
  };

  return (-energyAt(2.0 * delta) + 8.0 * energyAt(delta) - 8.0 * energyAt(-delta) + energyAt(-2.0 * delta)) /
         (12.0 * delta);
}
}  // namespace

// The intra-molecular exclusion / completion term q_i q_j (V(r) - 1/r) of the Wolf method is distance
// dependent and must contribute to the strain derivative. This checks the full 3x3 strain tensor returned by
// computeEwaldFourierEnergyStrainDerivative against a finite-difference of computeEwaldFourierEnergy (the self
// term is strain-independent, so the finite difference isolates the exclusion contribution).
TEST(wolf_strain, exclusion_strain_matches_finite_difference)
{
  double delta = 1e-5;
  double tolerance = 1e-6;

  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  forceField.useCharge = true;
  forceField.chargeMethod = ForceField::ChargeMethod::Wolf;
  forceField.omitEwaldFourier = true;
  forceField.EwaldAlpha = 0.25;

  System system = makeTenCO2(forceField);

  std::span<AtomDynamics> moleculeDynamics = system.spanOfMoleculeDynamics();
  for (AtomDynamics& dyn : moleculeDynamics) dyn.gradient = double3(0.0, 0.0, 0.0);

  std::pair<EnergyStatus, double3x3> pressureInfo = Interactions::computeEwaldFourierEnergyStrainDerivative(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
      system.forceField, system.simulationBox, system.framework, system.components,
      system.numberOfMoleculesPerComponent, system.atomData, system.atomDynamics, system.netChargeFramework,
      system.netChargePerComponent);
  pressureInfo.first.sumTotal();

  // The energy decomposition must reproduce the Wolf self + exclusion energy of computeEwaldFourierEnergy.
  RunningEnergy reference = Interactions::computeEwaldFourierEnergy(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
      system.forceField, system.simulationBox, system.components, system.numberOfMoleculesPerComponent,
      system.spanOfMoleculeAtoms());
  EXPECT_NE(reference.ewaldFourier(), 0.0);
  EXPECT_NEAR(pressureInfo.first.totalEnergy.energy, reference.ewaldFourier(), 1e-8);

  const std::array<std::pair<double3x3, double>, 9> strains{
      std::pair{double3x3{double3{1, 0, 0}, double3{0, 0, 0}, double3{0, 0, 0}}, pressureInfo.second.ax},
      std::pair{double3x3{double3{0, 1, 0}, double3{0, 0, 0}, double3{0, 0, 0}}, pressureInfo.second.bx},
      std::pair{double3x3{double3{0, 0, 1}, double3{0, 0, 0}, double3{0, 0, 0}}, pressureInfo.second.cx},
      std::pair{double3x3{double3{0, 0, 0}, double3{1, 0, 0}, double3{0, 0, 0}}, pressureInfo.second.ay},
      std::pair{double3x3{double3{0, 0, 0}, double3{0, 1, 0}, double3{0, 0, 0}}, pressureInfo.second.by},
      std::pair{double3x3{double3{0, 0, 0}, double3{0, 0, 1}, double3{0, 0, 0}}, pressureInfo.second.cy},
      std::pair{double3x3{double3{0, 0, 0}, double3{0, 0, 0}, double3{1, 0, 0}}, pressureInfo.second.az},
      std::pair{double3x3{double3{0, 0, 0}, double3{0, 0, 0}, double3{0, 1, 0}}, pressureInfo.second.bz},
      std::pair{double3x3{double3{0, 0, 0}, double3{0, 0, 0}, double3{0, 0, 1}}, pressureInfo.second.cz}};

  for (const auto& [direction, analytic] : strains)
  {
    double numerical = strainFiniteDifference(system, direction, delta);
    EXPECT_NEAR(numerical, analytic, tolerance * std::max({1.0, std::abs(numerical), std::abs(analytic)}))
        << "Wolf exclusion strain derivative disagrees with finite difference";
  }
}

// Integration check: the Wolf molecular (center-of-mass based) excess pressure must equal the center-of-mass
// scaling -dU/dV used by the volume move. Because the self term is volume independent and the intra-molecular
// exclusion is invariant under rigid center-of-mass scaling, this validates that the new strain contributions
// are correctly removed by the atomic-to-molecular virial correction while the inter-molecular Wolf pair virial
// is retained.
TEST(wolf_strain, molecular_pressure_matches_com_scaling)
{
  double tolerance = 1e-5;

  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  forceField.useCharge = true;
  forceField.chargeMethod = ForceField::ChargeMethod::Wolf;
  forceField.omitEwaldFourier = true;
  forceField.EwaldAlpha = 0.25;
  for (VDWParameters& parameters : forceField.data)
  {
    if (parameters.type == VDWParameters::Type::LennardJones)
      parameters.type = VDWParameters::Type::LennardJonesSecondOrderTaylorShifted;
  }
  forceField.preComputeDerivedParameters();
  forceField.preComputePotentialShift();

  System system = makeTenCO2(forceField);

  std::pair<EnergyStatus, double3x3> pressureInfo = system.computeMolecularPressure();
  double volume = system.simulationBox.volume;
  double analyticExcessPressure = pressureInfo.second.trace() / (3.0 * volume);

  double delta = 1e-5;
  SimulationBox boxPlus = system.simulationBox.scaled(std::cbrt(1.0 + delta));
  SimulationBox boxMinus = system.simulationBox.scaled(std::cbrt(1.0 - delta));
  auto posPlus = system.scaledCenterOfMassPositions(system.simulationBox, boxPlus);
  auto posMinus = system.scaledCenterOfMassPositions(system.simulationBox, boxMinus);

  auto totalEnergy = [&](const SimulationBox& box, std::span<const Atom> atoms)
  {
    RunningEnergy inter = Interactions::computeInterMolecularEnergy(system.forceField, box, atoms) +
                          Interactions::computeInterMolecularTailEnergy(system.forceField, box, atoms);
    RunningEnergy ewald = Interactions::computeEwaldFourierEnergy(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
        system.forceField, box, system.components, system.numberOfMoleculesPerComponent, atoms);
    return (inter + ewald).potentialEnergy();
  };

  double dUdV = (totalEnergy(boxPlus, posPlus.second) - totalEnergy(boxMinus, posMinus.second)) /
                (boxPlus.volume - boxMinus.volume);

  EXPECT_NEAR(analyticExcessPressure, -dUdV, tolerance) << "Wolf molecular excess pressure disagrees with -dU/dV";
}

static double maxFieldDifference(std::span<const double3> a, std::span<const double3> b)
{
  double m = 0.0;
  for (std::size_t i = 0; i < a.size(); ++i) m = std::max(m, (a[i] - b[i]).length());
  return m;
}

// Polarization with the Wolf method is handled entirely in real space (there is no reciprocal field). The
// incrementally-maintained electric field must therefore stay consistent with a from-scratch rebuild across a
// sequence of accepted/rejected translation moves, and the running polarization energy must match a full
// recomputation.
TEST(wolf_polarization, stored_field_consistency_translation)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  forceField.useCharge = true;
  forceField.chargeMethod = ForceField::ChargeMethod::Wolf;
  forceField.omitEwaldFourier = true;
  forceField.EwaldAlpha = 0.25;
  forceField.computePolarization = true;
  forceField.omitInterPolarization = false;
  forceField.omitInterInteractions = false;

  Framework f = Framework::makeITQ29(forceField, int3(2, 2, 2));
  Component c = Component::makeCO2(forceField, 0, true);
  System system = System(forceField, std::nullopt, false, 300.0, 1e4, 1.0, {f}, {c}, {}, {10}, 5);

  RandomNumber random(42);

  RunningEnergy running = system.computeTotalEnergies();
  EXPECT_NE(running.polarization, 0.0);

  for (std::size_t step = 0; step < 400; ++step)
  {
    std::size_t selectedMolecule = system.randomMoleculeOfComponent(random, 0);
    std::optional<RunningEnergy> e = MC_Moves::translationMove(random, system, 0, selectedMolecule);
    if (e.has_value()) running += e.value();
  }

  std::span<double3> maintained = system.spanOfMoleculeElectricField();
  std::vector<double3> snapshot(maintained.begin(), maintained.end());

  system.computeTotalElectricField();
  std::span<double3> rebuilt = system.spanOfMoleculeElectricField();

  EXPECT_LT(maxFieldDifference(snapshot, rebuilt), 1e-8);

  RunningEnergy recomputed = system.computeTotalEnergies();
  EXPECT_NEAR(running.polarization - recomputed.polarization, 0.0, 1e-6);
  EXPECT_NEAR(running.potentialEnergy() - recomputed.potentialEnergy(), 0.0, 1e-6);
}

// Molecular excess pressure must include the polarization strain. With polarizable rigid molecules the
// analytic COM-based pressure (now with polarization cellGradient folded into the strain tensor) is
// checked against center-of-mass-scaling -dU/dV of the full energy including U_pol.
TEST(wolf_polarization, molecular_pressure_includes_polarization_strain)
{
  double tolerance = 1e-5;

  ForceField forceField({{"P", false, 15.0, 0.8, 1.5, 8, false}, {"N", false, 14.0, -0.4, 1.2, 8, false}},
                        {{60.0, 3.0}, {40.0, 3.2}}, ForceField::MixingRule::Lorentz_Berthelot, 9.0, 9.0, 9.0, true, false,
                        true);
  forceField.computePolarization = true;
  forceField.omitInterPolarization = false;
  forceField.omitInterInteractions = false;
  forceField.omitEwaldFourier = true;
  forceField.chargeMethod = ForceField::ChargeMethod::Wolf;
  forceField.EwaldAlpha = 0.25;
  for (VDWParameters& parameters : forceField.data)
  {
    if (parameters.type == VDWParameters::Type::LennardJones)
      parameters.type = VDWParameters::Type::LennardJonesSecondOrderTaylorShifted;
  }
  forceField.preComputeDerivedParameters();
  forceField.preComputePotentialShift();

  const std::size_t typeP = *forceField.findPseudoAtom("P");
  const std::size_t typeN = *forceField.findPseudoAtom("N");
  Component molecule =
      Component(forceField, "rigid-polar", 100.0, 1e6, 0.2,
                {Atom({-1.0, 0.0, 0.0}, -0.4, 1.0, 0, static_cast<std::uint16_t>(typeN), 0, false, false),
                 Atom({0.0, 0.0, 0.0}, 0.8, 1.0, 0, static_cast<std::uint16_t>(typeP), 0, false, false),
                 Atom({1.0, 0.0, 0.0}, -0.4, 1.0, 0, static_cast<std::uint16_t>(typeN), 0, false, false)},
                ConnectivityTable(3), Potentials::IntraMolecularPotentials{}, 5, 21);
  molecule.rigid = true;

  System system =
      System(forceField, SimulationBox(23.0, 21.0, 19.0), false, 300.0, 1.5e5, 1.0, {}, {molecule}, {}, {2}, 5);
  system.moleculeData[0].centerOfMassPosition = double3(7.3, 9.1, 8.4);
  system.moleculeData[1].centerOfMassPosition = double3(12.2, 11.4, 10.7);
  // Place atoms consistently with the stored COMs (rigid linear molecule along x).
  {
    std::span<Atom> atoms = system.spanOfMoleculeAtoms();
    atoms[0].position = double3(6.3, 9.1, 8.4);
    atoms[1].position = double3(7.3, 9.1, 8.4);
    atoms[2].position = double3(8.3, 9.1, 8.4);
    atoms[3].position = double3(11.2, 11.4, 10.7);
    atoms[4].position = double3(12.2, 11.4, 10.7);
    atoms[5].position = double3(13.2, 11.4, 10.7);
  }

  system.computeTotalElectricField();
  const double polarizationEnergy = system.computePolarizationEnergy().polarization;
  EXPECT_NE(polarizationEnergy, 0.0);

  std::pair<EnergyStatus, double3x3> pressureInfo = system.computeMolecularPressure();
  EXPECT_NEAR(pressureInfo.first.polarizationEnergy.energy, polarizationEnergy, 1e-10);

  const double analyticExcessPressure = pressureInfo.second.trace() / (3.0 * system.simulationBox.volume);

  // Without polarization the pressure must differ once U_pol is strain-dependent.
  System systemNoPol = system;
  systemNoPol.forceField.computePolarization = false;
  const double excessWithoutPolarization =
      systemNoPol.computeMolecularPressure().second.trace() / (3.0 * systemNoPol.simulationBox.volume);
  EXPECT_NE(analyticExcessPressure, excessWithoutPolarization);

  double delta = 1e-5;
  SimulationBox boxPlus = system.simulationBox.scaled(std::cbrt(1.0 + delta));
  SimulationBox boxMinus = system.simulationBox.scaled(std::cbrt(1.0 - delta));
  auto posPlus = system.scaledCenterOfMassPositions(system.simulationBox, boxPlus);
  auto posMinus = system.scaledCenterOfMassPositions(system.simulationBox, boxMinus);

  auto totalEnergy = [&](System& reference, const SimulationBox& box, std::span<const Atom> atoms)
  {
    RunningEnergy inter = Interactions::computeInterMolecularEnergy(reference.forceField, box, atoms) +
                          Interactions::computeInterMolecularTailEnergy(reference.forceField, box, atoms);
    RunningEnergy ewald = Interactions::computeEwaldFourierEnergy(
        reference.eik_x, reference.eik_y, reference.eik_z, reference.eik_xy, reference.fixedFrameworkStoredEik,
        reference.storedEik, reference.forceField, box, reference.components, reference.numberOfMoleculesPerComponent,
        atoms);

    System strained = reference;
    strained.simulationBox = box;
    std::span<Atom> strainedAtoms = strained.spanOfMoleculeAtoms();
    for (std::size_t i = 0; i < atoms.size(); ++i) strainedAtoms[i] = atoms[i];
    const std::size_t numberOfFrameworkAtoms = strained.spanOfFrameworkAtoms().size();
    const std::size_t numberOfMoleculeAtoms = strainedAtoms.size();
    std::vector<std::uint8_t> movable(numberOfFrameworkAtoms + numberOfMoleculeAtoms, 0);
    for (std::size_t atom = 0; atom < numberOfMoleculeAtoms; ++atom) movable[numberOfFrameworkAtoms + atom] = 1;
    const double polarization =
        Interactions::computePolarizationDerivatives(strained, movable, {}, false, true).energy;

    return (inter + ewald).potentialEnergy() + polarization;
  };

  double dUdV = (totalEnergy(system, boxPlus, posPlus.second) - totalEnergy(system, boxMinus, posMinus.second)) /
                (boxPlus.volume - boxMinus.volume);

  EXPECT_NEAR(analyticExcessPressure, -dUdV, tolerance)
      << "Polarizable molecular excess pressure disagrees with -dU/dV";
}
