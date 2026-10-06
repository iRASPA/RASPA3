#include <gtest/gtest.h>

import std;

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
import running_energy;
import randomnumbers;
import potential_pair_derivatives;
import potential_pair_vdw;
import potential_correction_vdw;
import potential_correction_pressure;
import integrators_update;
import spatial_decomposition_settings;
import spatial_decomposition_force_engine;
import spatial_decomposition_device_step;

// The switched Lennard-Jones truncations of the biomolecular force fields: the quintic potential switch of
// OpenMM / GROMACS ("TruncationMethod": "switched") and the CHARMM force switch ("force-switched"), both acting on
// [SwitchingDistance, cutoff]. The force field converts its Lennard-Jones pairs to the switched forms; the pair
// potential, the tail corrections, and every kernel of the spatial-decomposition engine (scalar, SIMD cluster,
// device) evaluate them.
namespace
{
constexpr double epsilonKelvin = 119.8;
constexpr double sigma = 3.405;
constexpr double cutOff = 10.0;
constexpr double switchingDistance = 8.0;

VDWParameters makeParameters(VDWParameters::Type type)
{
  VDWParameters parameters(epsilonKelvin, sigma, type);
  parameters.computeDerivedParameters(cutOff, 300.0, switchingDistance);
  return parameters;
}

double lennardJones(double r)
{
  const double s6 = std::pow(sigma / r, 6.0);
  return 4.0 * epsilonKelvin * Units::KelvinToEnergy * (s6 * s6 - s6);
}

double lennardJonesDerivative(double r)
{
  const double s6 = std::pow(sigma / r, 6.0);
  return 4.0 * epsilonKelvin * Units::KelvinToEnergy * (-12.0 * s6 * s6 + 6.0 * s6) / r;
}

// Simpson quadrature of f on [a, b]
double integrate(const std::function<double(double)>& f, double a, double b, std::size_t intervals = 20000)
{
  const double h = (b - a) / static_cast<double>(intervals);
  double sum = f(a) + f(b);
  for (std::size_t i = 1; i < intervals; ++i) sum += (i % 2 == 1 ? 4.0 : 2.0) * f(a + static_cast<double>(i) * h);
  return sum * h / 3.0;
}

ForceField makeArgonForceField(ForceField::TruncationMethod method)
{
  ForceField forceField({{"Ar", false, 39.948, 0.0, 0.0, 18, false}, {"Kr", false, 83.798, 0.0, 0.0, 36, false}},
                        {{epsilonKelvin, sigma}, {164.0, 3.63}}, ForceField::MixingRule::Lorentz_Berthelot, cutOff,
                        cutOff, cutOff, false, false, false);
  forceField.setTruncationMethod(method, switchingDistance);
  return forceField;
}

Component makeAtom(const ForceField& forceField, std::string name, std::uint16_t type)
{
  return Component(forceField, name, 150.0, 4.8e6, 0.0, {Atom({0.0, 0.0, 0.0}, 0.0, 1.0, 0, type, 0, false, false)},
                   {}, {}, 5, 21);
}
}  // namespace

TEST(vdw_switching, switched_forms_are_continuous_and_vanish_at_the_cutoff)
{
  for (VDWParameters::Type type :
       {VDWParameters::Type::LennardJonesSwitched, VDWParameters::Type::LennardJonesForceSwitched})
  {
    const VDWParameters parameters = makeParameters(type);
    EXPECT_EQ(parameters.parameters2.x, cutOff);
    EXPECT_EQ(parameters.parameters2.y, switchingDistance);

    // zero at and beyond the cutoff
    EXPECT_NEAR(parameters.potentialEnergyAtFullCoupling(cutOff * cutOff), 0.0, 1e-14);
    EXPECT_NEAR(parameters.radialDerivativeAtFullCoupling(cutOff), 0.0, 1e-14);
    EXPECT_EQ(parameters.potentialEnergyAtFullCoupling(11.0 * 11.0), 0.0);

    // continuous energy and force at the switching distance
    const double below = switchingDistance * (1.0 - 1e-9);
    const double above = switchingDistance * (1.0 + 1e-9);
    const double scale = std::abs(lennardJones(switchingDistance));
    EXPECT_NEAR(parameters.potentialEnergyAtFullCoupling(below * below),
                parameters.potentialEnergyAtFullCoupling(above * above), 1e-6 * scale);
    EXPECT_NEAR(parameters.radialDerivativeAtFullCoupling(below), parameters.radialDerivativeAtFullCoupling(above),
                1e-6 * scale);

    // the potential switch leaves the potential untouched below the switching distance and multiplies it by S(r)
    // above; the force switch is Lennard-Jones minus a constant below
    for (double r : {3.5, 5.0, 7.9})
    {
      const double energy = parameters.potentialEnergyAtFullCoupling(r * r);
      if (type == VDWParameters::Type::LennardJonesSwitched)
      {
        EXPECT_NEAR(energy, lennardJones(r), 1e-12 * std::abs(lennardJones(r)));
      }
      else
      {
        const double offset =
            4.0 * epsilonKelvin * Units::KelvinToEnergy *
            (std::pow(sigma, 12) / (std::pow(cutOff, 6) * std::pow(switchingDistance, 6)) -
             std::pow(sigma, 6) / (std::pow(cutOff, 3) * std::pow(switchingDistance, 3)));
        EXPECT_NEAR(energy, lennardJones(r) - offset, 1e-12 * std::abs(lennardJones(r)));
      }
      EXPECT_NEAR(parameters.radialDerivativeAtFullCoupling(r), lennardJonesDerivative(r),
                  1e-12 * std::abs(lennardJonesDerivative(r)));
    }
    if (type == VDWParameters::Type::LennardJonesSwitched)
    {
      const double r = 9.0;
      const double x = (r - switchingDistance) / (cutOff - switchingDistance);
      const double s = 1.0 - 10.0 * x * x * x + 15.0 * x * x * x * x - 6.0 * x * x * x * x * x;
      EXPECT_NEAR(parameters.potentialEnergyAtFullCoupling(r * r), lennardJones(r) * s,
                  1e-12 * std::abs(lennardJones(r)));
    }

    // the radial derivative against finite differences across the switching region
    for (double r : {8.3, 9.0, 9.7})
    {
      const double h = 1e-5;
      const double numerical = (parameters.potentialEnergyAtFullCoupling((r + h) * (r + h)) -
                                parameters.potentialEnergyAtFullCoupling((r - h) * (r - h))) /
                               (2.0 * h);
      EXPECT_NEAR(parameters.radialDerivativeAtFullCoupling(r), numerical, 1e-7 * std::abs(numerical) + 1e-12)
          << VDWParameters::nameOfType(type) << " r = " << r;
    }
  }
}

TEST(vdw_switching, tail_corrections_cover_the_switching_region)
{
  for (VDWParameters::Type type :
       {VDWParameters::Type::LennardJonesSwitched, VDWParameters::Type::LennardJonesForceSwitched})
  {
    const VDWParameters parameters = makeParameters(type);
    // the energy correction: Integrate[(U_LJ - U_switched) r^2, {r, rs, Infinity}] = the switched-away part on
    // [rs, rc] plus the plain tail beyond the cutoff (U_switched is zero there)
    // (numerically up to 400 Angstrom, the remaining -C6 / r^6 tail in closed form)
    const double farAway = 400.0;
    const double c6 = 4.0 * epsilonKelvin * Units::KelvinToEnergy * std::pow(sigma, 6);
    const double expectedEnergy =
        integrate([&](double r) { return (lennardJones(r) - parameters.potentialEnergyAtFullCoupling(r * r)) * r * r; },
                  switchingDistance, cutOff) +
        integrate([&](double r) { return lennardJones(r) * r * r; }, cutOff, farAway) -
        c6 / (3.0 * std::pow(farAway, 3));
    const double energy = Potentials::potentialCorrectionVDW(parameters, cutOff);
    EXPECT_NEAR(energy, expectedEnergy, 1e-8 * std::abs(expectedEnergy)) << VDWParameters::nameOfType(type);

    const double expectedPressure =
        integrate([&](double r)
                  { return (lennardJonesDerivative(r) - parameters.radialDerivativeAtFullCoupling(r)) * r * r * r; },
                  switchingDistance, cutOff) +
        integrate([&](double r) { return lennardJonesDerivative(r) * r * r * r; }, cutOff, farAway) +
        2.0 * c6 / std::pow(farAway, 3);
    const double pressure = Potentials::potentialCorrectionPressure(parameters, cutOff);
    EXPECT_NEAR(pressure, expectedPressure, 1e-8 * std::abs(expectedPressure)) << VDWParameters::nameOfType(type);

    // the removed part is attractive-dominated at these distances: the switched correction exceeds the plain one
    VDWParameters plain(epsilonKelvin, sigma);
    plain.computeDerivedParameters(cutOff, 300.0);
    EXPECT_LT(energy, Potentials::potentialCorrectionVDW(plain, cutOff));
  }
}

TEST(vdw_switching, force_field_converts_its_lennard_jones_pairs)
{
  for (ForceField::TruncationMethod method :
       {ForceField::TruncationMethod::Switched, ForceField::TruncationMethod::ForceSwitched})
  {
    const ForceField forceField = makeArgonForceField(method);
    const VDWParameters::Type expected = method == ForceField::TruncationMethod::Switched
                                             ? VDWParameters::Type::LennardJonesSwitched
                                             : VDWParameters::Type::LennardJonesForceSwitched;
    for (std::size_t i = 0; i < 2; ++i)
    {
      for (std::size_t j = 0; j < 2; ++j)
      {
        EXPECT_EQ(forceField(i, j).type, expected);
        EXPECT_EQ(forceField(i, j).parameters2.x, cutOff);
        EXPECT_EQ(forceField(i, j).parameters2.y, switchingDistance);
        EXPECT_EQ(forceField(i, j).shift, 0.0);
        EXPECT_FALSE(forceField.shiftPotentials[i * 2 + j]);
      }
    }
    // the mixed pair kept its Lorentz-Berthelot parameters
    EXPECT_NEAR(forceField(0, 1).parameters.x * Units::EnergyToKelvin, std::sqrt(epsilonKelvin * 164.0), 1e-10);
    EXPECT_NEAR(forceField(0, 1).parameters.y, 0.5 * (sigma + 3.63), 1e-12);

    // the pair potential of the force field is the switched form
    const double r = 9.0;
    const Potentials::PairDerivatives<1> value = Potentials::potentialVDW<1>(forceField, 1.0, 1.0, r * r, 0, 0);
    EXPECT_NEAR(value.energy, forceField(0, 0).potentialEnergyAtFullCoupling(r * r), 1e-12);
    EXPECT_NEAR(value.firstDerivativeFactor * r, forceField(0, 0).radialDerivativeAtFullCoupling(r), 1e-12);
    EXPECT_EQ(Potentials::potentialVDW<0>(forceField, 1.0, 1.0, cutOff * cutOff, 0, 0).energy, 0.0);

    // back to plain truncation
    ForceField truncated = forceField;
    truncated.setTruncationMethod(ForceField::TruncationMethod::Truncated);
    EXPECT_EQ(truncated(0, 1).type, VDWParameters::Type::LennardJones);
    ForceField shifted = forceField;
    shifted.setTruncationMethod(ForceField::TruncationMethod::Shifted);
    EXPECT_EQ(shifted(0, 0).type, VDWParameters::Type::LennardJones);
    EXPECT_NEAR(shifted(0, 0).shift, lennardJones(cutOff), 1e-14);
  }
}

TEST(vdw_switching, force_field_file_reads_truncation_method_and_switching_distance)
{
  const std::filesystem::path directory = std::filesystem::temp_directory_path() / "raspa_vdw_switching_test";
  std::filesystem::create_directories(directory);
  auto write = [&](std::string_view truncation, std::string_view switching)
  {
    std::ofstream file(directory / "force_field.json");
    file << std::format(R"({{
  "MixingRule": "Lorentz-Berthelot",
  "TruncationMethod": "{}",{}
  "TailCorrections": true,
  "CutOff": 10.0,
  "ChargeMethod": "None",
  "PseudoAtoms": [
    {{"name": "A", "framework": false, "element": "C", "mass": 14.0, "charge": 0.0}},
    {{"name": "B", "framework": false, "element": "O", "mass": 16.0, "charge": 0.0}}
  ],
  "SelfInteractions": [
    {{"name": "A", "type": "lennard-jones", "parameters": [80.0, 3.5], "parameters14": [20.0, 3.0]}},
    {{"name": "B", "type": "lennard-jones", "parameters": [50.0, 3.0]}}
  ]
}})",
                        truncation, switching);
  };

  write("force-switched", "\n  \"SwitchingDistance\": 8.5,");
  const ForceField forceSwitched((directory / "force_field.json").string());
  EXPECT_EQ(forceSwitched.truncationMethod, ForceField::TruncationMethod::ForceSwitched);
  EXPECT_EQ(forceSwitched.switchingDistance, 8.5);
  EXPECT_EQ(forceSwitched(0, 1).type, VDWParameters::Type::LennardJonesForceSwitched);
  EXPECT_EQ(forceSwitched(0, 1).parameters2.y, 8.5);
  // the 1-4 table is switched alike; the tail corrections are of the switched form
  ASSERT_TRUE(forceSwitched.hasPair14Parameters());
  EXPECT_EQ(forceSwitched.pair14(0, 0).type, VDWParameters::Type::LennardJonesForceSwitched);
  EXPECT_NEAR(forceSwitched.pair14(0, 0).parameters.x * Units::EnergyToKelvin, 20.0, 1e-10);
  EXPECT_NEAR(forceSwitched(0, 0).tailCorrectionEnergy,
              Potentials::potentialCorrectionVDW(forceSwitched(0, 0), 10.0), 1e-14);

  write("switched", "");
  const ForceField switched((directory / "force_field.json").string());
  EXPECT_EQ(switched.truncationMethod, ForceField::TruncationMethod::Switched);
  EXPECT_EQ(switched(1, 1).type, VDWParameters::Type::LennardJonesSwitched);
  EXPECT_EQ(switched(1, 1).parameters2.y, 8.0);  // the default: 2 Angstrom below the cutoff

  write("shifted", "");
  const ForceField shifted((directory / "force_field.json").string());
  EXPECT_EQ(shifted.truncationMethod, ForceField::TruncationMethod::Shifted);
  EXPECT_EQ(shifted(1, 1).type, VDWParameters::Type::LennardJones);
  EXPECT_TRUE(shifted.shiftPotentials[3]);

  write("rounded", "");
  EXPECT_THROW(ForceField((directory / "force_field.json").string()), std::runtime_error);
  write("switched", "\n  \"SwitchingDistance\": 12.0,");
  EXPECT_THROW(ForceField((directory / "force_field.json").string()), std::runtime_error);
  std::filesystem::remove_all(directory);
}

TEST(vdw_switching, spatial_decomposition_kernels_evaluate_the_switched_forms)
{
  RandomNumber random(31);
  for (ForceField::TruncationMethod method :
       {ForceField::TruncationMethod::Switched, ForceField::TruncationMethod::ForceSwitched})
  {
    const ForceField forceField = makeArgonForceField(method);
    const SimulationBox box(26.0, 26.0, 26.0);
    // a jittered lattice of 128 atoms of two kinds: many pairs fall in the switching region
    std::vector<double3> argon{}, krypton{};
    for (std::size_t ix = 0; ix < 4; ++ix)
      for (std::size_t iy = 0; iy < 4; ++iy)
        for (std::size_t iz = 0; iz < 4; ++iz)
        {
          const double3 position = 6.5 * double3(static_cast<double>(ix) + 0.5, static_cast<double>(iy) + 0.5,
                                                 static_cast<double>(iz) + 0.5) +
                                   1.5 * double3(random.uniform() - 0.5, random.uniform() - 0.5, random.uniform() - 0.5);
          ((ix + iy + iz) % 2 == 0 ? argon : krypton).push_back(position);
          const double3 second = position + double3(3.9, 0.0, 0.0);
          ((ix + iy + iz) % 2 == 0 ? krypton : argon).push_back(second);
        }
    System system = System(forceField, box, false, 300.0, 1e4, 1.0, {},
                           {makeAtom(forceField, "argon", 0), makeAtom(forceField, "krypton", 1)}, {argon, krypton},
                           {0, 0}, 5);
    ASSERT_EQ(system.spanOfMoleculeAtoms().size(), 128uz);

    const RunningEnergy reference = Integrators::updateGradients(
        system.moleculeData, system.spanOfMoleculeAtoms(), system.spanOfMoleculeDynamics(),
        system.spanOfFrameworkAtoms(), system.forceField, system.simulationBox, system.components, system.eik_x,
        system.eik_y, system.eik_z, system.eik_xy, system.trialEik, system.fixedFrameworkStoredEik,
        system.interpolationGrids, system.numberOfMoleculesPerComponent, system.framework,
        system.spanOfFrameworkDynamics(), &system.crossLinks);
    std::vector<double3> referenceGradient{};
    double gradientNorm = 0.0;
    for (const AtomDynamics& dynamics : system.spanOfMoleculeDynamics())
    {
      referenceGradient.push_back(dynamics.gradient);
      gradientNorm += double3::dot(dynamics.gradient, dynamics.gradient);
    }
    gradientNorm = std::sqrt(gradientNorm / 128.0);
    ASSERT_NE(reference.moleculeMoleculeVDW, 0.0);

    // the reference (the generic pair potential) against a direct sum of the switched form over the minimum
    // images
    double direct = 0.0;
    std::span<const Atom> atoms = system.spanOfMoleculeAtoms();
    for (std::size_t i = 0; i < atoms.size(); ++i)
    {
      for (std::size_t j = i + 1; j < atoms.size(); ++j)
      {
        const double3 dr = box.applyPeriodicBoundaryConditions(atoms[i].position - atoms[j].position);
        const double rr = double3::dot(dr, dr);
        if (rr < cutOff * cutOff) direct += forceField(atoms[i].type, atoms[j].type).potentialEnergyAtFullCoupling(rr);
      }
    }
    EXPECT_NEAR(reference.moleculeMoleculeVDW, direct, 1e-10 * std::abs(direct));

    auto check = [&](SpatialDecompositionSettings settings, double tolerance, std::string_view label)
    {
      SpatialDecompositionForceEngine engine(settings);
      engine.initialize(system);
      EXPECT_TRUE(engine.usesFastKernel()) << label;
      const RunningEnergy energy = engine.computeGradients(system, true);
      EXPECT_NEAR(energy.moleculeMoleculeVDW, reference.moleculeMoleculeVDW,
                  tolerance * std::abs(reference.moleculeMoleculeVDW))
          << label;
      std::size_t i = 0;
      double difference = 0.0;
      for (const AtomDynamics& dynamics : system.spanOfMoleculeDynamics())
      {
        difference += double3::dot(dynamics.gradient - referenceGradient[i], dynamics.gradient - referenceGradient[i]);
        ++i;
      }
      EXPECT_NEAR(std::sqrt(difference / 128.0), 0.0, tolerance * gradientNorm) << label;
    };

    const std::string name = ForceField::truncationMethodName(method);
    for (std::size_t threads : {1uz, 3uz})
    {
      SpatialDecompositionSettings scalar;
      scalar.numberOfThreads = threads;
      scalar.verletSkin = 1.0;
      check(scalar, 1e-10, std::format("{} scalar kernel, {} threads", name, threads));

      SpatialDecompositionSettings cluster;
      cluster.numberOfThreads = threads;
      cluster.verletSkin = 1.0;
      cluster.clusterKernelForDouble = true;
      check(cluster, 1e-10, std::format("{} cluster kernel (double), {} threads", name, threads));

      SpatialDecompositionSettings mixed;
      mixed.numberOfThreads = threads;
      mixed.verletSkin = 1.0;
      mixed.pairPrecision = PairPrecision::Mixed;
      check(mixed, 2e-5, std::format("{} cluster kernel (mixed), {} threads", name, threads));
    }

    for (PairDevice pairDevice : {PairDevice::OpenCL, PairDevice::Metal, PairDevice::CUDA})
    {
      if (!DeviceStep::available(pairDevice)) continue;
      SpatialDecompositionSettings settings;
      settings.numberOfThreads = 2;
      settings.verletSkin = 1.0;
      settings.pairDevice = pairDevice;
      check(settings, 2e-5, std::format("{} {}", name, pairDeviceName(pairDevice)));
    }
  }
}
