module;

module polymerization;

import std;

import archive;
import double3;
import atom;
import molecule;
import system;
import component;
import running_energy;
import input_reader;
import randomnumbers;
import json;
import polymer_links;
import polymerization_site;
import polymerization_settings;
import molecular_dynamics_spatial_decomposition;
import force_engine;

Polymerization::Polymerization(InputReader& reader) : md(reader), settings(reader.polymerizationSettings)
{
  // the overlay is evaluated on the host: the complete-step-on-device integrator and the device bonded kernels
  // (which bake the atom types of the components) are not used
  md.engineSettings.resident = false;
  md.engineSettings.deviceBonded = false;
  reactionsFinished.assign(md.systems.size(), 0);
  numberOfScans.assign(md.systems.size(), 0);
  activationCycle.assign(md.systems.size(), 0);

  for (std::size_t systemId = 0; const System& system : md.systems)
  {
    if (!system.crossLinks.bondTypes.empty())
    {
      throw std::runtime_error(std::format(
          "Polymerization: system {} defines 'CrossLinkBonds' (the equilibrium cross-link moves); they cannot be "
          "combined with the chain-growth polymerization\n",
          systemId));
    }
    bool anySite = false;
    for (const Component& component : system.components) anySite |= !component.polymerizationSites.empty();
    if (!anySite)
    {
      throw std::runtime_error(std::format(
          "Polymerization: no component of system {} declares 'PolymerizationSites' in its molecule file\n", systemId));
    }
    if (system.polymerLinks.templates.empty())
    {
      throw std::runtime_error(std::format(
          "Polymerization: system {} has no 'Polymerization' block with the bonded tables of the product\n",
          systemId));
    }
    ++systemId;
  }
}

void Polymerization::createOutputFiles()
{
  // restart: the MD driver's streams are reopened in append mode (a fresh run opens them in setup())
  md.createOutputFiles();
  openBondLogs();
}

void Polymerization::openBondLogs()
{
  if (!settings.writeBondLog || !md.outputToFiles) return;
  const std::ios::openmode mode =
      (md.simulationStage != MolecularDynamicsSpatialDecomposition::SimulationStage::Uninitialized) ? std::ios::app
                                                                                                      : std::ios::out;
  std::filesystem::create_directories("output");
  bondLogs.clear();
  for (std::size_t systemId = 0; systemId < md.systems.size(); ++systemId)
  {
    bondLogs.emplace_back(std::format("output/polymerization_bonds.s{}.txt", systemId), mode);
    if (mode == std::ios::out)
    {
      std::print(bondLogs.back(),
                 "# step cycle event site_a(component molecule site) site_b(component molecule site) distance[A] "
                 "conversion radicals bonds\n");
    }
  }
}

void Polymerization::writeSettings()
{
  if (!md.outputToFiles) return;
  for (std::size_t systemId = 0; const System& system : md.systems)
  {
    std::ostream stream(md.streams[systemId].rdbuf());
    const double scanTime = static_cast<double>(settings.scanEvery) * system.timeStep;
    std::print(stream, "Polymerization settings\n");
    std::print(stream, "========================================================================================================================\n");
    std::print(stream, "    reaction scan every:             {} MD steps ({:g} ps)\n", settings.scanEvery, scanTime);
    std::print(stream, "    reactions start at cycle:        {}\n", settings.startCycle);
    std::print(stream, "    target conversion:               {:g}\n", settings.targetConversion);
    std::print(stream, "    relaxation (force ramp):         {} MD steps\n", settings.relaxationSteps);
    std::print(stream, "    initiation:                      {} initial radicals, probability per site and scan {:g}\n",
               settings.initialRadicals,
               PolymerizationSettings::probabilityPerScan(settings.initiationProbability,
                                                          settings.initiationRateConstant, scanTime));
    std::print(stream, "    propagation:                     capture radius {:g} [A], probability per pair and scan {:g}, intramolecular {}\n",
               settings.propagationCaptureRadius,
               PolymerizationSettings::probabilityPerScan(settings.propagationProbability,
                                                          settings.propagationRateConstant, scanTime),
               settings.allowIntramolecular ? "allowed" : "forbidden");
    if (settings.terminationEnabled)
    {
      std::print(stream, "    termination (combination):       capture radius {:g} [A], probability per pair and scan {:g}\n",
                 settings.terminationCaptureRadius,
                 PolymerizationSettings::probabilityPerScan(settings.terminationProbability,
                                                            settings.terminationRateConstant, scanTime));
    }
    else
    {
      std::print(stream, "    termination:                     disabled\n");
    }
    std::print(stream, "{}", system.polymerLinks.templates.printStatus());
    for (std::size_t componentId = 0; const Component& component : system.components)
    {
      for (const PolymerizationSite& site : component.polymerizationSites)
      {
        std::print(stream, "    site of component {} ({}):        type '{}', head {} -> '{}', tail {} -> '{}'{}\n",
                   componentId, component.name, site.siteType, site.head, site.headProductType, site.tail,
                   site.tailProductType,
                   site.initiatedHeadType.has_value() ? std::format(", initiated head '{}'", *site.initiatedHeadType)
                                                      : std::string{});
      }
      ++componentId;
    }
    std::print(stream, "    engine:                          resident integrator and device bonded kernels disabled (host overlay)\n");
    std::print(stream, "\n\n");
    std::flush(stream);
    ++systemId;
  }
}

void Polymerization::run()
{
  md.afterSystemStep = [this](std::size_t systemId) { afterStep(systemId); };
  md.extraStatus = [this](std::size_t systemId) { return status(systemId); };
  md.checkpointWriter = [this]() { writeBinaryRestartFile(*this); };

  // resumed from a restart: the derived caches (layout, term lists) are not archived
  for (System& system : md.systems)
  {
    if (system.polymerLinks.activated)
    {
      system.polymerLinks.activate(system.components, system.numberOfMoleculesPerComponent);
      system.polymerLinks.relaxationSteps = settings.relaxationSteps;
      system.polymerLinks.rebuild(system.forceField, system.components, system.numberOfMoleculesPerComponent,
                                  system.spanOfMoleculeAtoms());
    }
  }

  using Stage = MolecularDynamicsSpatialDecomposition::SimulationStage;
  switch (md.simulationStage)
  {
    case Stage::Uninitialized:
      md.setup();
      openBondLogs();
      writeSettings();
      [[fallthrough]];
    case Stage::PreInitialization:
      md.preInitialize();
      [[fallthrough]];
    case Stage::Initialization:
      md.initialize();
      [[fallthrough]];
    case Stage::Equilibration:
      md.equilibrate();
      [[fallthrough]];
    case Stage::Production:
      md.production();
      break;
    default:
      throw std::runtime_error("Polymerization::run(): no resume dispatch for the checkpointed simulation stage");
  }

  if (md.outputToFiles)
  {
    for (std::size_t systemId = 0; System& system : md.systems)
    {
      std::ostream stream(md.streams[systemId].rdbuf());
      std::print(stream, "\n");
      std::print(stream, "Final polymerization state\n");
      std::print(stream, "========================================================================================================================\n");
      std::print(stream, "    reaction scans:                  {}\n", numberOfScans[systemId]);
      std::print(stream, "{}", system.polymerLinks.printStatus(system.simulationBox, system.components, system.spanOfMoleculeAtoms()));
      std::flush(stream);
      ++systemId;
    }
    writeTopology();
  }
  md.tearDown();
}

std::string Polymerization::status(std::size_t systemId)
{
  const System& system = md.systems[systemId];
  if (!system.polymerLinks.activated) return {};
  return system.polymerLinks.printStatus(system.simulationBox, system.components, system.spanOfMoleculeAtoms());
}

void Polymerization::afterStep(std::size_t systemId)
{
  if (md.simulationStage != MolecularDynamicsSpatialDecomposition::SimulationStage::Production) return;
  if (md.currentCycle < settings.startCycle) return;

  System& system = md.systems[systemId];
  PolymerLinkTable& table = system.polymerLinks;
  if (!table.activated) activate(systemId);

  ++table.currentStep;
  if (reactionsFinished[systemId]) return;
  if (table.currentStep % settings.scanEvery != 0) return;

  const std::size_t accepted = scan(systemId);
  ++numberOfScans[systemId];
  if (accepted > 0) topologyChanged(systemId);

  if (table.conversion() >= settings.targetConversion)
  {
    reactionsFinished[systemId] = 1;
    if (md.outputToFiles)
    {
      std::ostream stream(md.streams[systemId].rdbuf());
      std::print(stream, "Polymerization: target conversion {:g} reached at cycle {} (step {}); the reactions stop\n\n",
                 settings.targetConversion, md.currentCycle, table.currentStep);
    }
  }
}

void Polymerization::activate(std::size_t systemId)
{
  System& system = md.systems[systemId];
  PolymerLinkTable& table = system.polymerLinks;
  table.activate(system.components, system.numberOfMoleculesPerComponent);
  table.relaxationSteps = settings.relaxationSteps;
  activationCycle[systemId] = md.currentCycle;

  // the initial radicals: random monomer sites
  std::vector<std::size_t> monomers{};
  for (std::size_t i = 0; i < table.numberOfSites(); ++i)
  {
    if (table.state(i) == PolymerSiteState::Monomer) monomers.push_back(i);
  }
  const std::size_t count = std::min(settings.initialRadicals, monomers.size());
  for (std::size_t n = 0; n < count; ++n)
  {
    const std::size_t pick = n + static_cast<std::size_t>(md.random.uniform() * static_cast<double>(monomers.size() - n));
    std::swap(monomers[n], monomers[std::min(pick, monomers.size() - 1)]);
    table.initiate(monomers[n], system.forceField, system.components, system.spanOfMoleculeAtoms(),
                   std::span<Molecule>(system.moleculeData), table.currentStep);
    if (settings.writeBondLog && systemId < bondLogs.size())
    {
      const PolymerAtomRef r = table.siteRef(monomers[n]);
      std::print(bondLogs[systemId], "{} {} initiation {} {} {} - - - 0 {:.6f} {} {}\n", table.currentStep,
                 md.currentCycle, r.componentId, r.moleculeIndex, r.atomIndex, table.conversion(),
                 table.countSites(PolymerSiteState::Radical), table.numberOfBonds());
    }
  }
  if (md.outputToFiles)
  {
    std::ostream stream(md.streams[systemId].rdbuf());
    std::print(stream, "Polymerization: activated at cycle {} with {} sites and {} initial radicals\n\n",
               md.currentCycle, table.numberOfSites(), count);
    std::flush(stream);
  }
  if (count > 0) topologyChanged(systemId);
}

std::size_t Polymerization::scan(std::size_t systemId)
{
  System& system = md.systems[systemId];
  PolymerLinkTable& table = system.polymerLinks;
  const std::vector<Component>& components = system.components;
  std::span<Atom> atoms = system.spanOfMoleculeAtoms();
  const double scanTime = static_cast<double>(settings.scanEvery) * system.timeStep;
  const double probabilityPropagation = PolymerizationSettings::probabilityPerScan(
      settings.propagationProbability, settings.propagationRateConstant, scanTime);
  const double probabilityTermination = PolymerizationSettings::probabilityPerScan(
      settings.terminationProbability, settings.terminationRateConstant, scanTime);
  const double probabilityInitiation = PolymerizationSettings::probabilityPerScan(
      settings.initiationProbability, settings.initiationRateConstant, scanTime);

  struct Site
  {
    std::size_t index;
    PolymerAtomRef ref;  // the atom that reacts (the tail of a radical, the head of a monomer)
    double3 position;
  };
  std::vector<Site> radicals{};
  std::vector<Site> monomers{};
  for (std::size_t i = 0; i < table.numberOfSites(); ++i)
  {
    switch (table.state(i))
    {
      case PolymerSiteState::Radical:
      {
        const PolymerAtomRef tail = table.tailOf(i, components);
        radicals.push_back(Site{i, tail, atoms[table.layout.index(tail)].position});
        break;
      }
      case PolymerSiteState::Monomer:
      {
        const PolymerAtomRef head = table.headOf(i, components);
        monomers.push_back(Site{i, head, atoms[table.layout.index(head)].position});
        break;
      }
      case PolymerSiteState::Dead:
        break;
    }
  }

  struct Candidate
  {
    double distance;
    std::size_t a;
    std::size_t b;
    PolymerBondKind kind;
  };
  std::vector<Candidate> candidates{};
  const double propagationRadiusSquared = settings.propagationCaptureRadius * settings.propagationCaptureRadius;
  const double terminationRadiusSquared = settings.terminationCaptureRadius * settings.terminationCaptureRadius;
  if (!radicals.empty())
  {
    for (const Site& radical : radicals)
    {
      for (const Site& monomer : monomers)
      {
        if (!settings.allowIntramolecular && radical.ref.sameMolecule(monomer.ref)) continue;
        const double3 dr = system.simulationBox.applyPeriodicBoundaryConditions(radical.position - monomer.position);
        const double rr = double3::dot(dr, dr);
        if (rr < propagationRadiusSquared)
        {
          candidates.push_back(Candidate{std::sqrt(rr), radical.index, monomer.index, PolymerBondKind::Propagation});
        }
      }
    }
    if (settings.terminationEnabled)
    {
      for (std::size_t i = 0; i < radicals.size(); ++i)
      {
        for (std::size_t j = i + 1; j < radicals.size(); ++j)
        {
          if (!settings.allowIntramolecular && radicals[i].ref.sameMolecule(radicals[j].ref)) continue;
          const double3 dr =
              system.simulationBox.applyPeriodicBoundaryConditions(radicals[i].position - radicals[j].position);
          const double rr = double3::dot(dr, dr);
          if (rr < terminationRadiusSquared)
          {
            candidates.push_back(
                Candidate{std::sqrt(rr), radicals[i].index, radicals[j].index, PolymerBondKind::Termination});
          }
        }
      }
    }
  }
  // closest pairs first; ties broken by the site indices (a deterministic order for a given random stream)
  std::ranges::sort(candidates, [](const Candidate& x, const Candidate& y)
                    { return std::tie(x.distance, x.a, x.b) < std::tie(y.distance, y.a, y.b); });

  std::size_t accepted = 0;
  std::unordered_set<std::size_t> reacted{};
  auto log = [&](std::string_view event, std::size_t a, std::size_t b, double distance)
  {
    if (!settings.writeBondLog || systemId >= bondLogs.size()) return;
    const PolymerAtomRef ra = table.siteRef(a);
    const PolymerAtomRef rb = table.siteRef(b);
    std::print(bondLogs[systemId], "{} {} {} {} {} {} {} {} {} {:.4f} {:.6f} {} {}\n", table.currentStep,
               md.currentCycle, event, ra.componentId, ra.moleculeIndex, ra.atomIndex, rb.componentId,
               rb.moleculeIndex, rb.atomIndex, distance, table.conversion(),
               table.countSites(PolymerSiteState::Radical), table.numberOfBonds());
  };
  for (const Candidate& candidate : candidates)
  {
    if (accepted >= settings.maximumReactionsPerScan) break;
    if (reacted.contains(candidate.a) || reacted.contains(candidate.b)) continue;
    const double probability =
        candidate.kind == PolymerBondKind::Propagation ? probabilityPropagation : probabilityTermination;
    if (probability <= 0.0) continue;
    if (probability < 1.0 && md.random.uniform() >= probability) continue;

    if (candidate.kind == PolymerBondKind::Propagation)
    {
      table.propagate(candidate.a, candidate.b, system.forceField, components, atoms,
                      std::span<Molecule>(system.moleculeData), table.currentStep, candidate.distance);
      log("propagation", candidate.a, candidate.b, candidate.distance);
    }
    else
    {
      table.terminate(candidate.a, candidate.b, system.forceField, components, atoms,
                      std::span<Molecule>(system.moleculeData), table.currentStep, candidate.distance);
      log("termination", candidate.a, candidate.b, candidate.distance);
    }
    reacted.insert(candidate.a);
    reacted.insert(candidate.b);
    ++accepted;
  }

  // continuous initiation of the monomer sites that did not react in this scan
  if (probabilityInitiation > 0.0)
  {
    for (const Site& monomer : monomers)
    {
      if (accepted >= settings.maximumReactionsPerScan) break;
      if (reacted.contains(monomer.index)) continue;
      if (md.random.uniform() >= probabilityInitiation) continue;
      table.initiate(monomer.index, system.forceField, components, atoms, std::span<Molecule>(system.moleculeData),
                     table.currentStep);
      if (settings.writeBondLog && systemId < bondLogs.size())
      {
        const PolymerAtomRef r = table.siteRef(monomer.index);
        std::print(bondLogs[systemId], "{} {} initiation {} {} {} - - - 0 {:.6f} {} {}\n", table.currentStep,
                   md.currentCycle, r.componentId, r.moleculeIndex, r.atomIndex, table.conversion(),
                   table.countSites(PolymerSiteState::Radical), table.numberOfBonds());
      }
      reacted.insert(monomer.index);
      ++accepted;
    }
  }
  if (accepted > 0 && settings.writeBondLog && systemId < bondLogs.size()) std::flush(bondLogs[systemId]);
  return accepted;
}

void Polymerization::topologyChanged(std::size_t systemId)
{
  System& system = md.systems[systemId];
  PolymerLinkTable& table = system.polymerLinks;

  // the derived term lists follow the new bonds / types (an incomplete template table is reported here)
  table.rebuild(system.forceField, system.components, system.numberOfMoleculesPerComponent,
                system.spanOfMoleculeAtoms());

  // the engine: per-type sums, list rebuild (the device types follow); the driver: tail-correction counts
  md.engines[systemId].topologyChanged(system);
  if (systemId < md.tailCorrectionCounts.size())
  {
    md.tailCorrectionCounts[systemId].numberOfAtoms = std::numeric_limits<std::size_t>::max();
  }

  // forces and potential energies of the new topology; the kinetic, thermostat and barostat parts of the step
  // are kept, and the conserved-energy reference restarts (the potential is not continuous across a reaction)
  const RunningEnergy previous = system.runningEnergies;
  md.recomputeGradients(systemId);
  system.runningEnergies.translationalKineticEnergy = previous.translationalKineticEnergy;
  system.runningEnergies.rotationalKineticEnergy = previous.rotationalKineticEnergy;
  system.runningEnergies.NoseHooverEnergy = previous.NoseHooverEnergy;
  system.runningEnergies.thermobarostatEnergy = previous.thermobarostatEnergy;
  system.conservedEnergy = system.runningEnergies.conservedEnergy();
  system.referenceEnergy = system.conservedEnergy;
}

void Polymerization::writeTopology()
{
  if (!settings.writeTopology) return;
  std::filesystem::create_directories("output");
  for (std::size_t systemId = 0; const System& system : md.systems)
  {
    const PolymerLinkTable& table = system.polymerLinks;
    nlohmann::json topology{};
    topology["status"] = table.jsonStatus(system.simulationBox, system.components, system.spanOfMoleculeAtoms());
    nlohmann::json componentsJson = nlohmann::json::array();
    for (std::size_t componentId = 0; const Component& component : system.components)
    {
      nlohmann::json entry{};
      entry["id"] = componentId;
      entry["name"] = component.name;
      entry["numberOfMolecules"] = system.numberOfMoleculesPerComponent[componentId];
      entry["numberOfAtoms"] = component.atoms.size();
      nlohmann::json sites = nlohmann::json::array();
      for (const PolymerizationSite& site : component.polymerizationSites)
      {
        sites.push_back({{"type", site.siteType}, {"head", site.head}, {"tail", site.tail}});
      }
      entry["sites"] = sites;
      componentsJson.push_back(entry);
      ++componentId;
    }
    topology["components"] = componentsJson;
    nlohmann::json states = nlohmann::json::array();
    for (std::size_t i = 0; i < table.numberOfSites(); ++i)
    {
      const PolymerAtomRef r = table.siteRef(i);
      states.push_back({r.componentId, r.moleculeIndex, r.atomIndex, static_cast<int>(table.state(i))});
    }
    topology["siteStates"] = states;
    topology["bonds"] = table.jsonTopology();

    std::ofstream stream(std::format("output/polymer_topology.s{}.json", systemId));
    stream << topology.dump(2) << "\n";
    ++systemId;
  }
}

Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const Polymerization& p)
{
  archive << p.versionNumber;
  archive << p.md;
  archive << p.settings;
  archive << p.reactionsFinished;
  archive << p.numberOfScans;
  archive << p.activationCycle;
#if DEBUG_ARCHIVE
  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif
  return archive;
}

Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, Polymerization& p)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > p.versionNumber)
  {
    throw std::runtime_error(std::format("Invalid version reading 'Polymerization' at line {} in file {}\n",
                                         std::source_location::current().line(),
                                         std::source_location::current().file_name()));
  }
  archive >> p.md;
  archive >> p.settings;
  archive >> p.reactionsFinished;
  archive >> p.numberOfScans;
  archive >> p.activationCycle;
#if DEBUG_ARCHIVE
  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(
        std::format("Polymerization: Invalid magic number {} at the end of the restart data\n", magicNumber));
  }
#endif
  // the engine settings of the restart file are those of the input; the restrictions of this driver apply again
  p.md.engineSettings.resident = false;
  p.md.engineSettings.deviceBonded = false;
  return archive;
}
