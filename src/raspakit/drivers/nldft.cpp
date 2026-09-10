module;

module nldft;

import std;

import stringutils;
import hardware_info;
import uint3;
import double3;
import double4;
import system;
import framework;
import forcefield;
import component;
import atom;
import units;
import input_reader;
import structure_input;
import crystal;
import pair_interactions;
import blocking_spheres;
import energy_shared_linear_probe;
import energy_shared_energy_backend;
import energy_backend;
import energy_opencl_backend;
import energy_shared_nldft;

namespace
{
constexpr std::size_t defaultMolecularOrientations = 128;
constexpr std::size_t defaultGridExtent = 128;

uint3 resolveGridSize(const System& system, const NLDFTParameters& parameters)
{
  if (parameters.gridSize.x > 0 && parameters.gridSize.y > 0 && parameters.gridSize.z > 0)
  {
    return parameters.gridSize;
  }
  if (system.forceField.numberOfVDWGridPoints.has_value())
  {
    return system.forceField.numberOfVDWGridPoints.value();
  }
  return uint3(defaultGridExtent, defaultGridExtent, defaultGridExtent);
}

LinearProbe makeProbe(const PairInteractions& interactions, const ForceField& forceField, const Component& component)
{
  if (std::optional<LinearProbe> named = LinearProbe::named(interactions, component.name); named.has_value())
  {
    return named.value();
  }

  if (component.definedAtoms.size() == 1)
  {
    const Atom& atom = component.definedAtoms.front().first;
    if (atom.type >= forceField.pseudoAtoms.size())
    {
      throw std::runtime_error(
          std::format("[NLDFT]: component '{}' has an out-of-range pseudo-atom type {}\n", component.name, atom.type));
    }
    const std::string& siteName = forceField.pseudoAtoms[atom.type].name;
    if (std::optional<LinearProbe> single = LinearProbe::singleSite(interactions, siteName); single.has_value())
    {
      return single.value();
    }
  }

  if (component.definedAtoms.size() >= 2)
  {
    double totalMass = 0.0;
    double3 centreOfMass{};
    for (const auto& [atom, mass] : component.definedAtoms)
    {
      centreOfMass += mass * atom.position;
      totalMass += mass;
    }
    if (!(totalMass > 0.0))
    {
      throw std::runtime_error(std::format("[NLDFT]: component '{}' has zero total mass\n", component.name));
    }
    centreOfMass = centreOfMass / totalMass;

    double3 axis{};
    double maxRadiusSquared = 0.0;
    for (const auto& [atom, mass] : component.definedAtoms)
    {
      const double3 delta = atom.position - centreOfMass;
      const double radiusSquared = delta.length_squared();
      if (radiusSquared > maxRadiusSquared)
      {
        maxRadiusSquared = radiusSquared;
        axis = delta;
      }
    }
    if (!(maxRadiusSquared > 0.0))
    {
      throw std::runtime_error(
          std::format("[NLDFT]: component '{}' is not linear (all sites coincide at the centre of mass)\n",
                      component.name));
    }
    axis = axis.normalized();

    LinearProbe probe;
    probe.name = component.name;
    probe.sites.reserve(component.definedAtoms.size());
    for (const auto& [atom, mass] : component.definedAtoms)
    {
      if (atom.type >= forceField.pseudoAtoms.size())
      {
        throw std::runtime_error(std::format("[NLDFT]: component '{}' has an out-of-range pseudo-atom type {}\n",
                                             component.name, atom.type));
      }
      const double offset = double3::dot(atom.position - centreOfMass, axis);
      probe.sites.push_back(
          LinearProbe::Site{atom.type, offset, atom.charge, forceField.pseudoAtoms[atom.type].name});
    }

    // Homogeneous diatomics (and CO2) look the same end-for-end; anything else needs a full sphere.
    probe.headTailSymmetric = true;
    for (const LinearProbe::Site& site : probe.sites)
    {
      const bool matched = std::ranges::any_of(probe.sites, [&](const LinearProbe::Site& other)
                                               {
                                                 return std::abs(other.offset + site.offset) < 1.0e-8 &&
                                                        other.type == site.type &&
                                                        std::abs(other.charge - site.charge) < 1.0e-12;
                                               });
      if (!matched)
      {
        probe.headTailSymmetric = false;
        break;
      }
    }
    return probe;
  }

  if (std::optional<LinearProbe> single = LinearProbe::singleSite(interactions, component.name); single.has_value())
  {
    return single.value();
  }

  throw std::runtime_error(std::format(
      "[NLDFT]: cannot build a linear probe from component '{}' (not a built-in shape, and no usable defined atoms)\n",
      component.name));
}

std::vector<BlockingSphere> blockingSpheresFromComponent(const Component& component)
{
  std::vector<BlockingSphere> spheres;
  spheres.reserve(component.blockingPockets.size());
  for (const double4& pocket : component.blockingPockets)
  {
    spheres.push_back(BlockingSphere{double3(pocket.x, pocket.y, pocket.z), pocket.w});
  }
  return spheres;
}
}  // namespace

NLDFT::NLDFT(InputReader& reader)
    : system(std::move(reader.systems.front())),
      parameters{.gridSize = reader.nldftGridSize,
                 .numberOfOrientations = reader.nldftNumberOfOrientations,
                 .useGPU = reader.nldftUseGPU,
                 .useElectrostatics = reader.nldftUseElectrostatics,
                 .relativePrecision = reader.nldftRelativePrecision}
{
  reader.systems.clear();

  if (!system.framework.has_value())
  {
    throw std::runtime_error("[NLDFT]: a framework is required (classical DFT needs an external field)\n");
  }
  if (system.components.empty())
  {
    throw std::runtime_error(
        "[NLDFT]: exactly one adsorbate component is required (it is the probe for V_ext / U(r, ω))\n");
  }
  if (system.components.size() != 1)
  {
    throw std::runtime_error(std::format(
        "[NLDFT]: exactly one adsorbate component is required, {} were declared\n", system.components.size()));
  }
}

NLDFT::NLDFT(System systemIn, NLDFTParameters parametersIn)
    : system(std::move(systemIn)), parameters(std::move(parametersIn))
{
  if (!system.framework.has_value())
  {
    throw std::runtime_error("[NLDFT]: a framework is required (classical DFT needs an external field)\n");
  }
  if (system.components.size() != 1)
  {
    throw std::runtime_error(std::format(
        "[NLDFT]: exactly one adsorbate component is required, {} were declared\n", system.components.size()));
  }
}

void NLDFT::run()
{
  setup();
  solve();
  output();
}

void NLDFT::setup()
{
  system.forceField.initializeAutomaticCutOff(system.simulationBox);
  system.forceField.initializeEwaldParameters(system.simulationBox);
  system.computeAutomaticBlockingPockets();

  if (!outputToFiles) return;

  std::filesystem::create_directories("output");
  const std::string fileName =
      std::format("output/output_{}_{}.s0.txt", system.temperature, system.input_pressure);
  stream = std::ofstream(fileName, std::ios::out);

  std::ostream out(stream.rdbuf());
  std::print(out, "{}", system.writeOutputHeader());
  std::print(out, "{}\n", HardwareInfo::writeInfo());
  std::print(out, "{}", Units::printStatus());
  std::print(out, "{}", system.writeSystemStatus());
  std::print(out, "{}", system.forceField.printPseudoAtomStatus());
  std::print(out, "{}", system.forceField.printForceFieldStatus());
  std::print(out, "{}", system.writeComponentStatus());
  std::print(out, "{}", system.writeNumberOfPseudoAtoms());
}

void NLDFT::solve()
{
  const Crystal crystal = StructureInput::makeCrystal(system.framework.value());
  const PairInteractions interactions = StructureInput::makeInteractions(system.forceField);
  const LinearProbe probe = makeProbe(interactions, system.forceField, system.components.front());
  const std::vector<BlockingSphere> spheres = blockingSpheresFromComponent(system.components.front());

  const uint3 gridSize = resolveGridSize(system, parameters);
  std::size_t orientations = parameters.numberOfOrientations;
  if (orientations == 0)
  {
    orientations = probe.sites.size() > 1 ? defaultMolecularOrientations : 1;
  }

  const bool useElectrostatics =
      parameters.useElectrostatics.value_or(system.forceField.useCharge) && probe.isCharged();
  const double relativePrecision = parameters.relativePrecision.value_or(system.forceField.EwaldPrecision);

  EnergyBackend backend = parameters.useGPU ? openCLEnergyBackend() : cpuEnergyBackend();

  if (outputToFiles)
  {
    std::ostream out(stream.rdbuf());
    std::print(out, "NLDFT probe: {} ({} site{}, {} orientation{})\n", probe.name, probe.sites.size(),
               probe.sites.size() == 1 ? "" : "s", orientations, orientations == 1 ? "" : "s");
    std::print(out, "NLDFT grid: {}x{}x{}, backend: {}, electrostatics: {}\n", gridSize.x, gridSize.y, gridSize.z,
               backend.name, useElectrostatics ? "on" : "off");
    if (!spheres.empty())
    {
      std::print(out, "NLDFT blocking spheres: {}\n", spheres.size());
    }
    std::print(out, "\n");
  }

  result.run(backend, interactions, crystal, probe, gridSize, orientations, system.temperature, spheres,
             useElectrostatics, relativePrecision);
}

void NLDFT::output()
{
  if (!outputToFiles) return;

  std::ostream out(stream.rdbuf());
  std::print(out, "NLDFT finished in {:.3f} s\n", result.seconds);
  std::print(out, "  mode: {}\n", result.molecular ? "ρ(r, ω)" : "spherical ρ(r)");
  std::print(out, "  model P0: {} Pa{}\n", result.bulk.saturationPressure,
             result.bulk.experimentalSaturation ? " (experimental N2 table; no bulk loop)" : "");
  std::print(out, "  Henry occupancy at experimental P0: {} / cell\n", result.henryMoleculesPerCell);
  std::print(out, "  BET gravimetric area: {} m^2/g (C = {}, n_m = {} / cell, window {}--{})\n",
             result.fit.gravimetricArea, result.fit.cConstant, result.fit.monolayerCapacity, result.fit.windowLow,
             result.fit.windowHigh);
  if (!result.isotherm.empty())
  {
    std::print(out, "  n(x -> 1): {} / cell\n", result.isotherm.back().moleculesPerCell);
  }
  if (result.unconverged > 0)
  {
    std::print(out, "  unconverged pressure points: {}\n", result.unconverged);
  }
  std::print(out, "\n");
}
