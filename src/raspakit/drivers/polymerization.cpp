module;

module polymerization;

import std;

import archive;
import double3;
import double3x3;
import atom;
import molecule;
import system;
import component;
import running_energy;
import input_reader;
import randomnumbers;
import json;
import units;
import polymer_links;
import polymerization_site;
import polymerization_settings;
import property_radical_waiting_times;
import molecular_dynamics_spatial_decomposition;
import force_engine;
import lammps_io;
import lammps_topology;

namespace
{
using ReactionEvent = Polymerization::ReactionEvent;
using ReactionKind = Polymerization::ReactionKind;

/// Whether 'value' is NaN or infinite. The build uses -ffast-math, under which the compiler folds std::isfinite,
/// 'value != value' and even an integer test of the bit pattern to a constant; the volatile copy keeps the test.
bool nonFinite(double value)
{
  volatile std::uint64_t bits = std::bit_cast<std::uint64_t>(value);
  const std::uint64_t copy = bits;
  return ((copy >> 52) & 0x7FFu) == 0x7FFu;
}
bool nonFinite(const double3& value) { return nonFinite(value.x) || nonFinite(value.y) || nonFinite(value.z); }

/// Records the current types of the atoms of the molecule of 'site' in the event (before the reaction).
void snapshotTypes(const PolymerLinkTable& table, std::span<const Atom> atoms, std::size_t site, ReactionEvent& event)
{
  const PolymerAtomRef r = table.siteRef(site);
  const std::size_t start = table.layout.moleculeStart(r.componentId, r.moleculeIndex);
  const std::size_t end = start + table.layout.atomsPerMolecule[r.componentId];
  for (std::size_t g = start; g < end; ++g)
  {
    // (both sites of a reaction within one molecule list the molecule once)
    const bool present = std::ranges::any_of(event.retyped, [&](const auto& t) { return std::get<0>(t) == g; });
    if (!present) event.retyped.emplace_back(static_cast<std::uint32_t>(g), atoms[g].type, atoms[g].type);
  }
}

/// After the reaction: keeps the atoms whose type changed, with the type they have now.
void completeRetyping(std::span<const Atom> atoms, ReactionEvent& event)
{
  for (auto& [g, before, after] : event.retyped) after = atoms[g].type;
  std::erase_if(event.retyped, [](const auto& t) { return std::get<1>(t) == std::get<2>(t); });
}

/// The event record of the initiation of 'site' (to complete after the reaction).
ReactionEvent initiationEvent(const PolymerLinkTable& table, const std::vector<Component>& components,
                              std::span<const Atom> atoms, std::size_t site)
{
  ReactionEvent event{};
  event.kind = ReactionKind::Initiation;
  event.siteA = site;
  event.atomA = static_cast<std::uint32_t>(table.layout.index(table.tailOf(site, components)));
  snapshotTypes(table, atoms, site, event);
  return event;
}
}  // namespace

Polymerization::Polymerization(InputReader& reader) : md(reader), settings(reader.polymerizationSettings)
{
  reactionsFinished.assign(md.systems.size(), 0);
  numberOfScans.assign(md.systems.size(), 0);
  activationCycle.assign(md.systems.size(), 0);
  stereoRejections.assign(md.systems.size(), 0);
  orientationRejections.assign(md.systems.size(), 0);
  acceptedStrain.assign(md.systems.size(), 0.0);
  acceptedStrainCount.assign(md.systems.size(), 0);
  reactionsStoppedStep.assign(md.systems.size(), 0);
  eventFrames.assign(md.systems.size(), 0);

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
    // tacticity control needs the side group of every tail that has one (a tail with several substituents and no
    // 'SideGroup' is ambiguous; a tail without substituents has no stereocentre and is left alone)
    if (settings.tacticity != Tacticity::None)
    {
      for (const Component& component : system.components)
      {
        for (const PolymerizationSite& site : component.polymerizationSites)
        {
          if (site.sideGroup.has_value() || site.tail >= component.connectivityTable.numberOfBeads) continue;
          std::size_t substituents = 0;
          for (std::size_t neighbour : component.connectivityTable.findAllNeighbors(site.tail))
          {
            if (neighbour != site.head) ++substituents;
          }
          if (substituents > 1)
          {
            throw std::runtime_error(std::format(
                "Polymerization: 'Tacticity' is set but the tail {} of the '{}' site of component '{}' has {} "
                "substituents besides the head; give the one that defines the configuration with 'SideGroup'\n",
                site.tail, site.siteType, component.name, substituents));
          }
        }
      }
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
  if (!md.outputToFiles) return;
  const std::ios::openmode mode =
      (md.simulationStage != MolecularDynamicsSpatialDecomposition::SimulationStage::Uninitialized) ? std::ios::app
                                                                                                      : std::ios::out;
  std::filesystem::create_directories("output");
  bondLogs.clear();
  eventLogs.clear();
  for (std::size_t systemId = 0; systemId < md.systems.size(); ++systemId)
  {
    if (settings.writeEventFrames)
    {
      eventLogs.emplace_back(std::format("output/polymerization_events.s{}.jsonl", systemId), mode);
    }
    if (!settings.writeBondLog) continue;
    bondLogs.emplace_back(std::format("output/polymerization_bonds.s{}.txt", systemId), mode);
    if (mode == std::ios::out)
    {
      std::print(bondLogs.back(),
                 "# step cycle event site_a(component molecule site) site_b(component molecule site) distance[A] "
                 "conversion radicals bonds stereo strain[K] wait[ps]\n"
                 "# stereo: configuration (+/-) of the reacting tail(s) and the dyad with the previous unit (m meso, "
                 "r racemo, . none); '-' no stereocentre\n"
                 "# strain: angular strain of the new bond at the reaction (with an 'OrientationTemperature'); '-' not "
                 "evaluated\n"
                 "# wait: waiting time of the reacting radical since its last reaction (termination: both radicals a/b); "
                 "'-' initiation\n");
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
    if (settings.propagationOrientationTemperature.has_value() || settings.terminationOrientationTemperature.has_value())
    {
      auto describe = [](std::optional<double> temperature)
      { return temperature.has_value() ? std::format("T* = {:g} K", temperature.value()) : std::string("distance only"); };
      std::print(stream, "    orientation weight:              p = p0 exp(-dU_angular / k_B T*); propagation {}, termination {}\n",
                 describe(settings.propagationOrientationTemperature),
                 describe(settings.terminationOrientationTemperature));
    }
    else
    {
      std::print(stream, "    orientation weight:              none (distance criterion only)\n");
    }
    if (system.propertyRadicalWaitingTimes.has_value())
    {
      std::print(stream, "    waiting times:                   dt = t(next) - t(last reaction) per living radical, {} bins per decade, written every {} cycles\n",
                 system.propertyRadicalWaitingTimes->binsPerDecade, system.propertyRadicalWaitingTimes->writeEvery);
    }
    if (settings.tacticity != Tacticity::None)
    {
      std::print(stream, "    tacticity:                       {}{}, wrong-face candidates {}\n",
                 PolymerizationSettings::tacticityName(settings.tacticity),
                 settings.tacticity == Tacticity::Atactic
                     ? std::format(" (meso probability {:g})", settings.mesoProbability)
                     : std::string{},
                 settings.stereoSelective ? "skipped" : "accepted and driven by the restraint");
      std::print(stream, "    stereocentre restraint:          1/2 k (V - chi V_0)^2 with k/k_B {:g} [K], V_0 {:g}\n",
                 settings.chiralForceConstant, settings.chiralReferenceVolume);
    }
    else
    {
      std::print(stream, "    tacticity:                       not controlled (dyad statistics are reported)\n");
    }
    std::print(stream, "{}", system.polymerLinks.templates.printStatus());
    for (std::size_t componentId = 0; const Component& component : system.components)
    {
      for (const PolymerizationSite& site : component.polymerizationSites)
      {
        std::print(stream, "    site of component {} ({}):        type '{}', head {} -> '{}', tail {} -> '{}'{}{}\n",
                   componentId, component.name, site.siteType, site.head, site.headProductType, site.tail,
                   site.tailProductType,
                   site.initiatedHeadType.has_value() ? std::format(", initiated head '{}'", *site.initiatedHeadType)
                                                      : std::string{},
                   site.sideGroup.has_value() ? std::format(", side group {}", *site.sideGroup)
                                              : std::string(", no stereocentre"));
      }
      ++componentId;
    }
    std::print(stream,
               "    engine:                          polymer overlay on the device with the device bonded kernels "
               "(resident integrator allowed), on the host otherwise\n");
    if (settings.writeEventFrames)
    {
      std::print(stream,
                 "    event frames:                    LAMMPS data file of {} in polymerization_events/s{}/ "
                 "(polymer bonds in the topology), records in output/polymerization_events.s{}.jsonl\n",
                 settings.eventFrameEvery == 1 ? std::string("every reaction frame")
                                               : std::format("one in {} reaction frames", settings.eventFrameEvery),
                 systemId, systemId);
    }
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
      system.polymerLinks.chiralForceConstant = settings.chiralForceConstant * Units::KelvinToEnergy;
      system.polymerLinks.chiralReferenceVolume = settings.chiralReferenceVolume;
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
      if (settings.tacticity != Tacticity::None)
      {
        std::print(stream, "    wrong-face candidates skipped:   {}\n", stereoRejections[systemId]);
      }
      std::print(stream, "{}", orientationText(systemId));
      std::print(stream, "{}", system.polymerLinks.printStatus(system.simulationBox, system.components, system.spanOfMoleculeAtoms()));
      std::print(stream, "{}", waitingTimeText(systemId));
      writeWaitingTimes(systemId, true);
      std::flush(stream);
      writeTopology(systemId);
      ++systemId;
    }
  }
  md.tearDown();
}

std::string Polymerization::orientationText(std::size_t systemId) const
{
  if (!settings.propagationOrientationTemperature.has_value() && !settings.terminationOrientationTemperature.has_value())
  {
    return {};
  }
  const std::size_t count = acceptedStrainCount[systemId];
  return std::format(
      "    orientation weight:              {} candidates rejected, mean angular strain of the accepted reactions {:.1f} [K]\n",
      orientationRejections[systemId], count > 0 ? acceptedStrain[systemId] / static_cast<double>(count) : 0.0);
}

std::string Polymerization::status(std::size_t systemId)
{
  const System& system = md.systems[systemId];
  if (!system.polymerLinks.activated) return {};
  return system.polymerLinks.printStatus(system.simulationBox, system.components, system.spanOfMoleculeAtoms()) +
         orientationText(systemId) + waitingTimeText(systemId);
}

std::uint64_t Polymerization::waitingTimeReferenceStep(std::size_t systemId) const
{
  return reactionsStoppedStep[systemId] > 0 ? reactionsStoppedStep[systemId]
                                            : md.systems[systemId].polymerLinks.currentStep;
}

std::string Polymerization::waitingTimeText(std::size_t systemId) const
{
  const System& system = md.systems[systemId];
  if (!system.propertyRadicalWaitingTimes.has_value()) return {};
  return system.propertyRadicalWaitingTimes->printStatus(system.polymerLinks, system.timeStep,
                                                         waitingTimeReferenceStep(systemId));
}

void Polymerization::writeWaitingTimes(std::size_t systemId, bool force) const
{
  const System& system = md.systems[systemId];
  if (!system.propertyRadicalWaitingTimes.has_value() || !md.outputToFiles) return;
  system.propertyRadicalWaitingTimes->writeOutput(systemId, md.currentCycle, system.polymerLinks, system.timeStep,
                                                  waitingTimeReferenceStep(systemId), force);
}

void Polymerization::afterStep(std::size_t systemId)
{
  if (md.simulationStage != MolecularDynamicsSpatialDecomposition::SimulationStage::Production) return;
  if (md.currentCycle < settings.startCycle) return;

  System& system = md.systems[systemId];
  PolymerLinkTable& table = system.polymerLinks;
  if (!table.activated)
  {
    // the activation retypes atoms and re-evaluates the forces on the host: the host state must be current
    md.engines[systemId].downloadResidentState(system);
    activate(systemId);
  }

  ++table.currentStep;
  writeWaitingTimes(systemId, false);
  if (reactionsFinished[systemId]) return;
  if (table.currentStep % settings.scanEvery != 0) return;

  // the scan reads the positions (and a reaction modifies the atoms and re-evaluates the forces on the host):
  // with the resident integrator the device state is fetched first (a no-op when the host copy is current)
  md.engines[systemId].downloadResidentState(system);
  const std::vector<ReactionEvent> events = scan(systemId);
  ++numberOfScans[systemId];
  if (!events.empty())
  {
    topologyChanged(systemId);
    writeEventFrame(systemId, events);
  }

  if (table.conversion() >= settings.targetConversion)
  {
    reactionsFinished[systemId] = 1;
    reactionsStoppedStep[systemId] = table.currentStep;
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
  table.chiralForceConstant = settings.chiralForceConstant * Units::KelvinToEnergy;
  table.chiralReferenceVolume = settings.chiralReferenceVolume;
  activationCycle[systemId] = md.currentCycle;

  // the initial radicals: random monomer sites
  std::vector<std::size_t> monomers{};
  for (std::size_t i = 0; i < table.numberOfSites(); ++i)
  {
    if (table.state(i) == PolymerSiteState::Monomer) monomers.push_back(i);
  }
  const std::size_t count = std::min(settings.initialRadicals, monomers.size());
  std::vector<ReactionEvent> events{};
  for (std::size_t n = 0; n < count; ++n)
  {
    const std::size_t pick = n + static_cast<std::size_t>(md.random.uniform() * static_cast<double>(monomers.size() - n));
    std::swap(monomers[n], monomers[std::min(pick, monomers.size() - 1)]);
    ReactionEvent event = initiationEvent(table, system.components, system.spanOfMoleculeAtoms(), monomers[n]);
    table.initiate(monomers[n], system.forceField, system.components, system.spanOfMoleculeAtoms(),
                   std::span<Molecule>(system.moleculeData), table.currentStep);
    completeRetyping(system.spanOfMoleculeAtoms(), event);
    logEvent(systemId, event);
    events.push_back(std::move(event));
  }
  if (md.outputToFiles)
  {
    std::ostream stream(md.streams[systemId].rdbuf());
    std::print(stream, "Polymerization: activated at cycle {} with {} sites and {} initial radicals\n\n",
               md.currentCycle, table.numberOfSites(), count);
    std::flush(stream);
  }
  if (count > 0)
  {
    topologyChanged(systemId);
    writeEventFrame(systemId, events);
  }
}

std::int8_t Polymerization::nextChirality(std::int8_t previous)
{
  if (previous == 0) return 0;
  switch (settings.tacticity)
  {
    case Tacticity::None:
      return 0;
    case Tacticity::Isotactic:
      return previous;
    case Tacticity::Syndiotactic:
      return static_cast<std::int8_t>(-previous);
    case Tacticity::Atactic:
      return md.random.uniform() < settings.mesoProbability ? previous : static_cast<std::int8_t>(-previous);
  }
  return 0;
}

void Polymerization::logEvent(std::size_t systemId, const ReactionEvent& event)
{
  if (!settings.writeBondLog || systemId >= bondLogs.size()) return;
  const PolymerLinkTable& table = md.systems[systemId].polymerLinks;
  const PolymerAtomRef ra = table.siteRef(event.siteA);
  if (event.kind == ReactionKind::Initiation)
  {
    std::print(bondLogs[systemId], "{} {} initiation {} {} {} - - - 0 {:.6f} {} {} - - -\n", table.currentStep,
               md.currentCycle, ra.componentId, ra.moleculeIndex, ra.atomIndex, table.conversion(),
               table.countSites(PolymerSiteState::Radical), table.numberOfBonds());
    return;
  }
  const PolymerAtomRef rb = table.siteRef(event.siteB);
  auto waitText = [](const std::optional<double>& wait)
  { return wait.has_value() ? std::format("{:.4g}", wait.value()) : std::string("-"); };
  const std::string wait = event.kind == ReactionKind::Termination
                               ? std::format("{}/{}", waitText(event.waitingTime[0]), waitText(event.waitingTime[1]))
                               : waitText(event.waitingTime[0]);
  std::print(bondLogs[systemId], "{} {} {} {} {} {} {} {} {} {:.4f} {:.6f} {} {} {} {} {}\n", table.currentStep,
             md.currentCycle, ReactionEvent::name(event.kind), ra.componentId, ra.moleculeIndex, ra.atomIndex,
             rb.componentId, rb.moleculeIndex, rb.atomIndex, event.distance, table.conversion(),
             table.countSites(PolymerSiteState::Radical), table.numberOfBonds(), event.stereo,
             event.strain.has_value() ? std::format("{:.1f}", event.strain.value()) : std::string("-"), wait);
}

std::vector<Polymerization::ReactionEvent> Polymerization::scan(std::size_t systemId)
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
  // Diverged dynamics show up here as NaN coordinates: under -ffast-math a NaN passes the capture-radius test
  // and breaks the ordering of the candidate sort (which then reads past its range), so stop with a message.
  auto sitePosition = [&](const PolymerAtomRef& ref)
  {
    const double3 position = atoms[table.layout.index(ref)].position;
    if (nonFinite(position))
    {
      throw std::runtime_error(std::format(
          "[Polymerization]: the position of atom {} of molecule {} (component {}) is not finite at step {} "
          "(cycle {}); the dynamics diverged (check the time step, the relaxation of new bonds and the forces)\n",
          ref.atomIndex, ref.moleculeIndex, ref.componentId, table.currentStep, md.currentCycle));
    }
    return position;
  };
  for (std::size_t i = 0; i < table.numberOfSites(); ++i)
  {
    switch (table.state(i))
    {
      case PolymerSiteState::Radical:
      {
        const PolymerAtomRef tail = table.tailOf(i, components);
        radicals.push_back(Site{i, tail, sitePosition(tail)});
        break;
      }
      case PolymerSiteState::Monomer:
      {
        const PolymerAtomRef head = table.headOf(i, components);
        monomers.push_back(Site{i, head, sitePosition(head)});
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

  std::vector<ReactionEvent> events{};
  std::unordered_set<std::size_t> reacted{};
  double lastStrain = 0.0;  // the angular strain [K] of the reaction just accepted (for the log)

  // The orientation weight of a candidate: exp(-dU / k_B T*) with dU the angular strain of the new bond in the
  // current geometry (the monomer's head and tail typed as in the product), in Kelvin.
  const bool weighted = settings.propagationOrientationTemperature.has_value() ||
                        settings.terminationOrientationTemperature.has_value();
  if (weighted && !candidates.empty() && !table.derivedValid)
  {
    table.rebuild(system.forceField, components, system.numberOfMoleculesPerComponent, atoms);
  }
  auto angularStrainInKelvin = [&](const Candidate& candidate) -> double
  {
    std::map<std::uint32_t, std::uint16_t> retyped{};
    const std::uint32_t a = static_cast<std::uint32_t>(table.layout.index(table.tailOf(candidate.a, components)));
    std::uint32_t b{};
    if (candidate.kind == PolymerBondKind::Propagation)
    {
      const PolymerAtomRef head = table.headOf(candidate.b, components);
      const PolymerAtomRef tail = table.tailOf(candidate.b, components);
      const PolymerAtomRef r = table.siteRef(candidate.b);
      const PolymerizationSite& definition = components[r.componentId].polymerizationSites[r.atomIndex];
      b = static_cast<std::uint32_t>(table.layout.index(head));
      const std::optional<std::size_t> headType = system.forceField.findPseudoAtom(definition.headProductType);
      const std::optional<std::size_t> tailType = system.forceField.findPseudoAtom(definition.tailProductType);
      if (headType.has_value()) retyped[b] = static_cast<std::uint16_t>(headType.value());
      if (tailType.has_value())
      {
        retyped[static_cast<std::uint32_t>(table.layout.index(tail))] = static_cast<std::uint16_t>(tailType.value());
      }
    }
    else
    {
      b = static_cast<std::uint32_t>(table.layout.index(table.tailOf(candidate.b, components)));
    }
    return table.angularStrain(a, b, system.simulationBox, components, atoms, retyped) * Units::EnergyToKelvin;
  };
  // the acceptance draw: probability times the orientation weight; records the strain of an accepted reaction
  auto accept = [&](const Candidate& candidate, double probability) -> bool
  {
    const std::optional<double> orientationTemperature = candidate.kind == PolymerBondKind::Propagation
                                                             ? settings.propagationOrientationTemperature
                                                             : settings.terminationOrientationTemperature;
    double strain = 0.0;
    double weight = 1.0;
    if (orientationTemperature.has_value())
    {
      strain = angularStrainInKelvin(candidate);
      weight = PolymerizationSettings::orientationWeight(orientationTemperature, strain);
    }
    const double p = probability * weight;
    if (p < 1.0 && md.random.uniform() >= p)
    {
      if (orientationTemperature.has_value()) ++orientationRejections[systemId];
      return false;
    }
    if (orientationTemperature.has_value())
    {
      acceptedStrain[systemId] += strain;
      ++acceptedStrainCount[systemId];
    }
    lastStrain = strain;
    return true;
  };

  // The configuration the tail of radical 'site' takes when it bonds to 'partner': the one required by the
  // tacticity when set, otherwise the one its geometry gives (the face of approach); none when the required
  // and the geometric configuration differ and such candidates are skipped. 0: the tail is not a stereocentre.
  auto resolveChirality = [&](std::size_t site, const PolymerAtomRef& partner) -> std::optional<std::int8_t>
  {
    const std::uint32_t partnerAtom = static_cast<std::uint32_t>(table.layout.index(partner));
    const std::int8_t geometric =
        table.geometricChirality(site, partnerAtom, system.simulationBox, components, atoms);
    if (geometric == 0) return std::int8_t{0};
    const std::int8_t required = table.chirality(site);
    if (required == 0) return settings.tacticity == Tacticity::None ? std::int8_t{0} : geometric;
    if (settings.stereoSelective && required != geometric) return std::nullopt;
    return required;
  };
  // the dyad of the reacting tail with the centre before it in the chain, for the log
  auto stereoText = [&](std::size_t site, std::int8_t chirality) -> std::string
  {
    if (chirality == 0) return "-";
    const std::size_t previous = table.predecessorOf(site, components);
    const std::int8_t before = previous == PolymerLinkTable::noSite ? std::int8_t{0} : table.chirality(previous);
    const char dyad = before == 0 ? '.' : (before == chirality ? 'm' : 'r');
    return std::format("{}{}", chirality > 0 ? '+' : '-', dyad);
  };
  // the waiting time of a radical since its last reaction (its creation), in ps
  auto waitingTime = [&](std::size_t radical)
  { return static_cast<double>(table.currentStep - std::min(table.currentStep, table.siteStep[radical])) * system.timeStep; };
  // the event record of a bond-forming candidate, before the reaction (the retyped atoms are completed after it)
  auto bondEvent = [&](const Candidate& candidate, std::uint32_t atomB)
  {
    ReactionEvent event{};
    event.kind = candidate.kind == PolymerBondKind::Propagation ? ReactionKind::Propagation : ReactionKind::Termination;
    event.siteA = candidate.a;
    event.siteB = candidate.b;
    event.atomA = static_cast<std::uint32_t>(table.layout.index(table.tailOf(candidate.a, components)));
    event.atomB = atomB;
    event.distance = candidate.distance;
    event.bond = table.bonds.size();  // the index the new bond takes
    snapshotTypes(table, atoms, candidate.a, event);
    snapshotTypes(table, atoms, candidate.b, event);
    return event;
  };
  for (const Candidate& candidate : candidates)
  {
    if (events.size() >= settings.maximumReactionsPerScan) break;
    if (reacted.contains(candidate.a) || reacted.contains(candidate.b)) continue;
    const double probability =
        candidate.kind == PolymerBondKind::Propagation ? probabilityPropagation : probabilityTermination;
    if (probability <= 0.0) continue;

    if (candidate.kind == PolymerBondKind::Propagation)
    {
      const PolymerAtomRef head = table.headOf(candidate.b, components);
      const std::optional<std::int8_t> chiralityA = resolveChirality(candidate.a, head);
      if (!chiralityA.has_value())
      {
        ++stereoRejections[systemId];
        continue;
      }
      if (!accept(candidate, probability)) continue;
      ReactionEvent event = bondEvent(candidate, static_cast<std::uint32_t>(table.layout.index(head)));
      event.chirality = {chiralityA.value(), nextChirality(chiralityA.value())};
      event.stereo = stereoText(candidate.a, chiralityA.value());
      event.waitingTime[0] = waitingTime(candidate.a);
      if (weighted) event.strain = lastStrain;
      if (system.propertyRadicalWaitingTimes.has_value())
      {
        system.propertyRadicalWaitingTimes->addPropagation(event.waitingTime[0].value());
      }
      table.propagate(candidate.a, candidate.b, system.forceField, components, atoms,
                      std::span<Molecule>(system.moleculeData), table.currentStep, candidate.distance,
                      event.chirality[0], event.chirality[1]);
      completeRetyping(atoms, event);
      logEvent(systemId, event);
      events.push_back(std::move(event));
    }
    else
    {
      const PolymerAtomRef tailB = table.tailOf(candidate.b, components);
      const std::optional<std::int8_t> chiralityA = resolveChirality(candidate.a, tailB);
      const std::optional<std::int8_t> chiralityB =
          resolveChirality(candidate.b, table.tailOf(candidate.a, components));
      if (!chiralityA.has_value() || !chiralityB.has_value())
      {
        ++stereoRejections[systemId];
        continue;
      }
      if (!accept(candidate, probability)) continue;
      ReactionEvent event = bondEvent(candidate, static_cast<std::uint32_t>(table.layout.index(tailB)));
      event.chirality = {chiralityA.value(), chiralityB.value()};
      event.stereo = std::format("{}/{}", stereoText(candidate.a, chiralityA.value()),
                                 stereoText(candidate.b, chiralityB.value()));
      event.waitingTime = {waitingTime(candidate.a), waitingTime(candidate.b)};
      if (weighted) event.strain = lastStrain;
      if (system.propertyRadicalWaitingTimes.has_value())
      {
        system.propertyRadicalWaitingTimes->addTermination(event.waitingTime[0].value());
        system.propertyRadicalWaitingTimes->addTermination(event.waitingTime[1].value());
      }
      table.terminate(candidate.a, candidate.b, system.forceField, components, atoms,
                      std::span<Molecule>(system.moleculeData), table.currentStep, candidate.distance,
                      event.chirality[0], event.chirality[1]);
      completeRetyping(atoms, event);
      logEvent(systemId, event);
      events.push_back(std::move(event));
    }
    reacted.insert(candidate.a);
    reacted.insert(candidate.b);
  }

  // continuous initiation of the monomer sites that did not react in this scan
  if (probabilityInitiation > 0.0)
  {
    for (const Site& monomer : monomers)
    {
      if (events.size() >= settings.maximumReactionsPerScan) break;
      if (reacted.contains(monomer.index)) continue;
      if (md.random.uniform() >= probabilityInitiation) continue;
      ReactionEvent event = initiationEvent(table, components, atoms, monomer.index);
      table.initiate(monomer.index, system.forceField, components, atoms, std::span<Molecule>(system.moleculeData),
                     table.currentStep);
      completeRetyping(atoms, event);
      logEvent(systemId, event);
      events.push_back(std::move(event));
      reacted.insert(monomer.index);
    }
  }
  if (!events.empty() && settings.writeBondLog && systemId < bondLogs.size()) std::flush(bondLogs[systemId]);
  return events;
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

void Polymerization::writeTopology(std::size_t systemId)
{
  if (!settings.writeTopology || !md.outputToFiles) return;
  std::filesystem::create_directories("output");
  const System& system = md.systems[systemId];
  const PolymerLinkTable& table = system.polymerLinks;
  nlohmann::json topology{};
  topology["step"] = table.currentStep;
  topology["time"] = static_cast<double>(table.currentStep) * system.timeStep;
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
      nlohmann::json siteJson{{"type", site.siteType}, {"head", site.head}, {"tail", site.tail}};
      if (site.sideGroup.has_value()) siteJson["sideGroup"] = site.sideGroup.value();
      sites.push_back(siteJson);
    }
    entry["sites"] = sites;
    componentsJson.push_back(entry);
    ++componentId;
  }
  topology["components"] = componentsJson;
  topology["tacticity"] = PolymerizationSettings::tacticityName(settings.tacticity);
  nlohmann::json states = nlohmann::json::array();
  nlohmann::json chiralities = nlohmann::json::array();
  for (std::size_t i = 0; i < table.numberOfSites(); ++i)
  {
    const PolymerAtomRef r = table.siteRef(i);
    states.push_back({r.componentId, r.moleculeIndex, r.atomIndex, static_cast<int>(table.state(i))});
    chiralities.push_back(static_cast<int>(table.chirality(i)));
  }
  topology["siteStates"] = states;
  topology["siteChirality"] = chiralities;
  topology["bonds"] = table.jsonTopology();

  std::ofstream stream(std::format("output/polymer_topology.s{}.json", systemId));
  stream << topology.dump(2) << "\n";
}

void Polymerization::writeEventFrame(std::size_t systemId, std::span<const ReactionEvent> events)
{
  if (!settings.writeEventFrames || !md.outputToFiles || events.empty()) return;
  const System& system = md.systems[systemId];
  const PolymerLinkTable& table = system.polymerLinks;
  const std::size_t frame = eventFrames[systemId]++;

  std::string dataFile{};
  if (frame % settings.eventFrameEvery == 0)
  {
    const std::filesystem::path directory = std::format("polymerization_events/s{}", systemId);
    std::filesystem::create_directories(directory);
    const std::string stem = std::format("step_{:010d}", table.currentStep);
    std::size_t initiations = 0, propagations = 0, terminations = 0;
    for (const ReactionEvent& event : events)
    {
      initiations += event.kind == ReactionKind::Initiation;
      propagations += event.kind == ReactionKind::Propagation;
      terminations += event.kind == ReactionKind::Termination;
    }
    LAMMPS::ExportOptions options{};
    options.dataFile = stem + ".data";
    options.tableFile = stem + ".table";
    options.pairListFile = stem + ".pairs";
    options.comment = std::format(
        "polymerization event frame: step {}, time {:.6g} ps, cycle {}, {} initiation(s), {} propagation(s), {} "
        "termination(s), {} polymer bonds in total (see output/polymerization_events.s{}.jsonl)",
        table.currentStep, static_cast<double>(table.currentStep) * system.timeStep, md.currentCycle, initiations,
        propagations, terminations, table.numberOfBonds(), systemId);
    const LAMMPS::ExportFiles files = LAMMPS::exportSystem(
        system.components, system.atomData, system.atomDynamics, system.moleculeData, system.simulationBox,
        system.forceField, system.numberOfIntegerMoleculesPerComponent, system.framework, options, &table);
    std::ofstream(directory / options.dataFile) << files.data;
    std::ofstream(directory / (stem + ".in")) << files.input;
    if (!files.table.empty()) std::ofstream(directory / options.tableFile) << files.table;
    if (!files.pairList.empty()) std::ofstream(directory / options.pairListFile) << files.pairList;
    dataFile = (directory / options.dataFile).generic_string();
  }

  if (systemId < eventLogs.size())
  {
    eventLogs[systemId] << eventRecord(systemId, events, dataFile).dump() << "\n";
    std::flush(eventLogs[systemId]);
  }
  writeTopology(systemId);
}

nlohmann::json Polymerization::eventRecord(std::size_t systemId, std::span<const ReactionEvent> events,
                                           std::string_view dataFile) const
{
  const System& system = md.systems[systemId];
  const PolymerLinkTable& table = system.polymerLinks;
  const std::span<const Atom> atoms = system.spanOfMoleculeAtoms();
  // the LAMMPS atom id of the global atom g of the molecule span: the framework atoms come first, ids are 1-based
  const std::size_t idOffset = system.spanOfFrameworkAtoms().size() + 1;
  auto typeName = [&](std::uint16_t type) { return system.forceField.pseudoAtoms[type].name; };
  auto atomJson = [&](std::uint32_t g)
  {
    const PolymerAtomRef r = table.layout.ref(g);
    return nlohmann::json{{"id", g},
                          {"lammpsId", idOffset + g},
                          {"component", r.componentId},
                          {"molecule", r.moleculeIndex},
                          {"atom", r.atomIndex},
                          {"type", typeName(atoms[g].type)}};
  };
  auto optionalNumber = [](const std::optional<double>& value)
  { return value.has_value() ? nlohmann::json(value.value()) : nlohmann::json(nullptr); };

  nlohmann::json record{};
  record["step"] = table.currentStep;
  record["time"] = static_cast<double>(table.currentStep) * system.timeStep;
  record["cycle"] = md.currentCycle;
  record["frame"] = eventFrames[systemId] > 0 ? eventFrames[systemId] - 1 : 0;
  record["dataFile"] = dataFile.empty() ? nlohmann::json(nullptr) : nlohmann::json(std::string(dataFile));
  const double3x3& cell = system.simulationBox.cell;
  record["box"] = {{cell.ax, cell.ay, cell.az}, {cell.bx, cell.by, cell.bz}, {cell.cx, cell.cy, cell.cz}};

  nlohmann::json eventsJson = nlohmann::json::array();
  nlohmann::json retyped = nlohmann::json::array();
  nlohmann::json newBonds = nlohmann::json::array();
  for (const ReactionEvent& event : events)
  {
    nlohmann::json entry{};
    entry["kind"] = ReactionEvent::name(event.kind);
    entry["siteA"] = event.siteA;
    entry["a"] = atomJson(event.atomA);
    if (event.kind != ReactionKind::Initiation)
    {
      entry["siteB"] = event.siteB;
      entry["b"] = atomJson(event.atomB);
      entry["distance"] = event.distance;
      entry["chirality"] = {event.chirality[0], event.chirality[1]};
      entry["stereo"] = event.stereo;
      entry["strain"] = optionalNumber(event.strain);
      entry["waitingTime"] = {optionalNumber(event.waitingTime[0]), optionalNumber(event.waitingTime[1])};
      entry["bond"] = event.bond;
    }
    eventsJson.push_back(entry);
    for (const auto& [g, before, after] : event.retyped)
    {
      retyped.push_back({{"id", g}, {"lammpsId", idOffset + g}, {"from", typeName(before)}, {"to", typeName(after)}});
    }
    if (event.bond < table.bonds.size())
    {
      const PolymerBond& bond = table.bonds[event.bond];
      nlohmann::json bondJson{};
      bondJson["bond"] = event.bond;
      bondJson["kind"] = bond.kind == PolymerBondKind::Propagation ? "propagation" : "termination";
      bondJson["atoms"] = {table.layout.index(bond.a), table.layout.index(bond.b)};
      bondJson["formationDistance"] = bond.formationDistance;
      // the bonded terms the bond created (global atoms; valid after the rebuild of the table)
      nlohmann::json bends = nlohmann::json::array();
      nlohmann::json torsions = nlohmann::json::array();
      for (const PolymerLinkTable::LinkBendTerm& term : table.linkBends)
      {
        if (term.bondId == event.bond) bends.push_back(term.atoms);
      }
      for (const PolymerLinkTable::LinkTorsionTerm& term : table.linkTorsions)
      {
        if (term.bondId == event.bond) torsions.push_back(term.atoms);
      }
      bondJson["bends"] = bends;
      bondJson["torsions"] = torsions;
      newBonds.push_back(bondJson);
    }
  }
  record["events"] = eventsJson;
  record["retyped"] = retyped;
  record["newBonds"] = newBonds;
  record["totals"] = {{"bonds", table.numberOfBonds()},
                      {"monomers", table.countSites(PolymerSiteState::Monomer)},
                      {"radicals", table.countSites(PolymerSiteState::Radical)},
                      {"dead", table.countSites(PolymerSiteState::Dead)},
                      {"conversion", table.conversion()},
                      {"relaxationSteps", table.relaxationSteps},
                      {"rampsActive", table.rampsActive()}};
  return record;
}

Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const Polymerization& p)
{
  archive << p.versionNumber;
  archive << p.md;
  archive << p.settings;
  archive << p.reactionsFinished;
  archive << p.numberOfScans;
  archive << p.activationCycle;
  archive << p.stereoRejections;
  archive << p.orientationRejections;
  archive << p.acceptedStrain;
  archive << p.acceptedStrainCount;
  archive << p.reactionsStoppedStep;
  archive << p.eventFrames;
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
  if (versionNumber >= 2)
  {
    archive >> p.stereoRejections;
  }
  else
  {
    p.stereoRejections.assign(p.md.systems.size(), 0);
  }
  if (versionNumber >= 3)
  {
    archive >> p.orientationRejections;
    archive >> p.acceptedStrain;
    archive >> p.acceptedStrainCount;
  }
  else
  {
    p.orientationRejections.assign(p.md.systems.size(), 0);
    p.acceptedStrain.assign(p.md.systems.size(), 0.0);
    p.acceptedStrainCount.assign(p.md.systems.size(), 0);
  }
  if (versionNumber >= 4)
  {
    archive >> p.reactionsStoppedStep;
  }
  else
  {
    p.reactionsStoppedStep.assign(p.md.systems.size(), 0);
  }
  if (versionNumber >= 5)
  {
    archive >> p.eventFrames;
  }
  else
  {
    p.eventFrames.assign(p.md.systems.size(), 0);
  }
#if DEBUG_ARCHIVE
  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(
        std::format("Polymerization: Invalid magic number {} at the end of the restart data\n", magicNumber));
  }
#endif
  return archive;
}
