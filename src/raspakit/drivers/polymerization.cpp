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
import reactive_group;
import polymer_reactions;
import polymerization_settings;
import property_radical_waiting_times;
import molecular_dynamics_spatial_decomposition;
import force_engine;
import lammps_io;
import lammps_topology;

namespace
{
using ReactionEvent = Polymerization::ReactionEvent;

/// Whether 'value' is NaN or infinite. The build uses -ffast-math, under which the compiler folds std::isfinite,
/// 'value != value' and even an integer test of the bit pattern to a constant; the volatile copy keeps the test.
bool nonFinite(double value)
{
  volatile std::uint64_t bits = std::bit_cast<std::uint64_t>(value);
  const std::uint64_t copy = bits;
  return ((copy >> 52) & 0x7FFu) == 0x7FFu;
}
bool nonFinite(const double3& value) { return nonFinite(value.x) || nonFinite(value.y) || nonFinite(value.z); }

/// Records the current types of the atoms of the molecule of 'group' in the event (before the reaction).
void snapshotTypes(const PolymerLinkTable& table, std::span<const Atom> atoms, std::size_t group, ReactionEvent& event)
{
  const PolymerAtomRef r = table.groupRef(group);
  const std::size_t start = table.layout.moleculeStart(r.componentId, r.moleculeIndex);
  const std::size_t end = start + table.layout.atomsPerMolecule[r.componentId];
  for (std::size_t g = start; g < end; ++g)
  {
    // (both groups of a reaction within one molecule list the molecule once)
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

std::string optionalText(const std::optional<double>& value, std::string_view format)
{
  return value.has_value() ? std::vformat(format, std::make_format_args(value.value())) : std::string("-");
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
          "combined with the polymerization\n",
          systemId));
    }
    bool anyGroup = false;
    for (const Component& component : system.components) anyGroup |= !component.reactiveGroups.empty();
    if (!anyGroup)
    {
      throw std::runtime_error(std::format(
          "Polymerization: no component of system {} declares 'ReactiveGroups' in its molecule file\n", systemId));
    }
    if (system.polymerLinks.templates.empty())
    {
      throw std::runtime_error(std::format(
          "Polymerization: system {} has no 'Polymerization' block with the bonded tables of the product\n",
          systemId));
    }
    if (!system.polymerLinks.scheme.resolved || system.polymerLinks.scheme.empty())
    {
      throw std::runtime_error(std::format(
          "Polymerization: system {} has no 'Reactions' in its 'Polymerization' block\n", systemId));
    }
    if (settings.tacticity != Tacticity::None && !system.polymerLinks.scheme.hasStereo())
    {
      throw std::runtime_error(std::format(
          "Polymerization: 'Tacticity' is set but no reaction of system {} declares a stereocentre ('Stereo')\n",
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
                 "# step cycle reaction group_a(component molecule group) group_b(component molecule group) "
                 "distance[A] conversion active bonds stereo strain[K] wait[ps]\n"
                 "# group_b, distance: '-' for a unimolecular reaction\n"
                 "# stereo: per stereocentre of the reaction the configuration (+/-) and the dyad with the centre "
                 "before it in the chain (m meso, r racemo, . none); '-' no stereocentre\n"
                 "# strain: angular strain of the new bond(s) at the reaction (with an 'OrientationTemperature'); "
                 "'-' not evaluated\n"
                 "# wait: waiting time of the reacting group(s) in an active state since their last reaction "
                 "(a/b for a bimolecular reaction); '-' not active\n");
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
    std::print(stream, "    maximum reactions per scan:      {}\n",
               settings.maximumReactionsPerScan == std::numeric_limits<std::size_t>::max()
                   ? std::string("unlimited")
                   : std::to_string(settings.maximumReactionsPerScan));
    std::print(stream, "    relaxation (force ramp):         {} MD steps\n", settings.relaxationSteps);
    if (system.polymerLinks.scheme.orientationWeighted())
    {
      std::print(stream, "    orientation weight:              p = p0 exp(-dU_angular / k_B T*) for the reactions with a T*\n");
    }
    else
    {
      std::print(stream, "    orientation weight:              none (distance criterion only)\n");
    }
    if (system.propertyRadicalWaitingTimes.has_value())
    {
      std::print(stream, "    waiting times:                   dt = t(next) - t(last reaction) per active group, {} bins per decade, written every {} cycles\n",
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
    std::print(stream, "\n");
    for (std::size_t componentId = 0; const Component& component : system.components)
    {
      for (std::size_t g = 0; const ReactiveGroup& group : component.reactiveGroups)
      {
        std::string atoms{};
        for (const auto& [role, atom] : group.atoms)
        {
          atoms += std::format("{}{} {} ({})", atoms.empty() ? "" : ", ", role, atom,
                               system.forceField.pseudoAtoms[component.atoms[atom].type].name);
        }
        std::print(stream, "    group {} of component {} ({}):   type '{}', state '{}', atoms {}\n", g, componentId,
                   component.name, group.type, group.initialState, atoms);
        ++g;
      }
      ++componentId;
    }
    std::print(stream, "\n");
    std::print(stream, "{}", system.polymerLinks.scheme.printStatus(scanTime));
    std::print(stream, "{}", system.polymerLinks.templates.printStatus());
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
  if (!md.systems[systemId].polymerLinks.scheme.orientationWeighted()) return {};
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

  // the initial events of the unimolecular rules (the initial radicals): random eligible groups
  std::vector<ReactionEvent> events{};
  std::vector<std::size_t> applied(table.scheme.rules.size(), 0);
  for (std::size_t r = 0; r < table.scheme.rules.size(); ++r)
  {
    const ReactionRule& rule = table.scheme.rules[r];
    if (rule.bimolecular() || rule.initialEvents == 0) continue;
    std::vector<std::size_t> eligible{};
    for (std::size_t i = 0; i < table.numberOfGroups(); ++i)
    {
      if (table.typeOf(i) == rule.reactantType[0] && table.state(i) == rule.reactantState[0]) eligible.push_back(i);
    }
    const std::size_t count = std::min(rule.initialEvents, eligible.size());
    for (std::size_t n = 0; n < count; ++n)
    {
      const std::size_t pick =
          n + static_cast<std::size_t>(md.random.uniform() * static_cast<double>(eligible.size() - n));
      std::swap(eligible[n], eligible[std::min(pick, eligible.size() - 1)]);
      const std::size_t group = eligible[n];
      events.push_back(react(systemId, r, std::span<const std::size_t>(&group, 1), {}, {noAtom, noAtom}, 0.0, {}));
    }
    applied[r] = count;
  }
  if (md.outputToFiles)
  {
    std::ostream stream(md.streams[systemId].rdbuf());
    std::string initial{};
    for (std::size_t r = 0; r < applied.size(); ++r)
    {
      if (applied[r] == 0) continue;
      initial += std::format(", {} initial '{}' event(s)", applied[r], table.scheme.rules[r].name);
    }
    std::print(stream, "Polymerization: activated at cycle {} with {} groups{}\n\n", md.currentCycle,
               table.numberOfGroups(), initial);
    std::flush(stream);
  }
  if (!events.empty())
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

Polymerization::ReactionEvent Polymerization::react(std::size_t systemId, std::size_t ruleIndex,
                                                    std::span<const std::size_t> groups,
                                                    std::span<const std::int8_t> chiralities,
                                                    std::array<std::uint32_t, 2> captureAtoms, double distance,
                                                    std::optional<double> strain)
{
  System& system = md.systems[systemId];
  PolymerLinkTable& table = system.polymerLinks;
  const std::vector<Component>& components = system.components;
  std::span<Atom> atoms = system.spanOfMoleculeAtoms();
  const ReactionRule& rule = table.scheme.rules[ruleIndex];

  ReactionEvent event{};
  event.rule = static_cast<std::uint16_t>(ruleIndex);
  event.distance = distance;
  event.strain = strain;
  event.chirality.assign(chiralities.begin(), chiralities.end());
  for (std::size_t i = 0; i < groups.size(); ++i)
  {
    event.groups[i] = groups[i];
    event.atoms[i] = captureAtoms[i] != noAtom ? captureAtoms[i] : table.atomOf(groups[i], 0, components);
    // the waiting time of a reacting group in an active state since its last reaction, in ps
    if (table.active(groups[i]))
    {
      const std::uint64_t since = table.currentStep - std::min(table.currentStep, table.groupStep[groups[i]]);
      event.waitingTime[i] = static_cast<double>(since) * system.timeStep;
      if (system.propertyRadicalWaitingTimes.has_value())
      {
        system.propertyRadicalWaitingTimes->add(ruleIndex, event.waitingTime[i].value());
      }
    }
    snapshotTypes(table, atoms, groups[i], event);
  }

  // the stereo field of the log: per centre its configuration and the dyad with the centre before it in the chain
  // (from the geometry when the configuration is not controlled)
  std::vector<std::int8_t> nextChiralities(chiralities.size(), 0);
  if (!rule.stereoAtoms.empty())
  {
    std::string text{};
    for (std::size_t k = 0; k < rule.stereoAtoms.size(); ++k)
    {
      const std::array<std::uint32_t, 4> centreAtoms = table.stereoAtoms(rule, k, groups, components);
      const std::int8_t sign = chiralities[k] != 0
                                   ? chiralities[k]
                                   : table.geometricChirality(centreAtoms, system.simulationBox, atoms);
      const std::size_t centreGroup = groups[rule.stereoAtoms[k].centre.reactant];
      const std::size_t previous = table.previousGroup[centreGroup];
      std::int8_t before = 0;
      if (previous != PolymerLinkTable::noGroup)
      {
        const std::optional<std::size_t> centre = table.centreOf(previous);
        if (centre.has_value())
        {
          const StereoCentre& c = table.centres[centre.value()];
          before = c.chirality != 0 ? c.chirality : table.geometricChirality(c.atoms, system.simulationBox, atoms);
        }
      }
      const char dyad = before == 0 ? '.' : (before == sign ? 'm' : 'r');
      text += std::format("{}{}{}", k > 0 ? "/" : "", sign > 0 ? '+' : '-', dyad);
      if (rule.stereoAtoms[k].next.has_value()) nextChiralities[k] = nextChirality(chiralities[k]);
    }
    event.stereo = text;
  }

  const std::size_t brokenBefore = table.brokenBonds.size();
  event.bonds = table.apply(ruleIndex, groups, chiralities, nextChiralities, system.forceField, components, atoms,
                            std::span<Molecule>(system.moleculeData), system.simulationBox, table.currentStep);
  for (std::size_t i = brokenBefore; i < table.brokenBonds.size(); ++i) event.brokenBonds.push_back(i);
  completeRetyping(atoms, event);
  logEvent(systemId, event);
  return event;
}

void Polymerization::logEvent(std::size_t systemId, const ReactionEvent& event)
{
  if (!settings.writeBondLog || systemId >= bondLogs.size()) return;
  const PolymerLinkTable& table = md.systems[systemId].polymerLinks;
  const PolymerAtomRef ra = table.groupRef(event.groups[0]);
  const std::string b = event.bimolecular()
                            ? std::format("{} {} {}", table.groupRef(event.groups[1]).componentId,
                                          table.groupRef(event.groups[1]).moleculeIndex,
                                          table.groupRef(event.groups[1]).atomIndex)
                            : std::string("- - -");
  const std::string wait = event.bimolecular()
                               ? std::format("{}/{}", optionalText(event.waitingTime[0], "{:.4g}"),
                                             optionalText(event.waitingTime[1], "{:.4g}"))
                               : optionalText(event.waitingTime[0], "{:.4g}");
  std::print(bondLogs[systemId], "{} {} {} {} {} {} {} {} {:.6f} {} {} {} {} {}\n", table.currentStep, md.currentCycle,
             table.scheme.rules[event.rule].name, ra.componentId, ra.moleculeIndex, ra.atomIndex, b,
             event.bimolecular() ? std::format("{:.4f}", event.distance) : std::string("-"), table.conversion(),
             table.countActive(), table.numberOfBonds(), event.stereo, optionalText(event.strain, "{:.1f}"), wait);
}

std::vector<Polymerization::ReactionEvent> Polymerization::scan(std::size_t systemId)
{
  System& system = md.systems[systemId];
  PolymerLinkTable& table = system.polymerLinks;
  const ReactionScheme& scheme = table.scheme;
  const std::vector<Component>& components = system.components;
  std::span<Atom> atoms = system.spanOfMoleculeAtoms();
  const double scanTime = static_cast<double>(settings.scanEvery) * system.timeStep;

  // Diverged dynamics show up here as NaN coordinates: under -ffast-math a NaN passes the capture-radius test
  // and breaks the ordering of the candidate sort (which then reads past its range), so stop with a message.
  auto positionOf = [&](std::uint32_t atom)
  {
    const double3 position = atoms[atom].position;
    if (nonFinite(position))
    {
      const PolymerAtomRef ref = table.layout.ref(atom);
      throw std::runtime_error(std::format(
          "[Polymerization]: the position of atom {} of molecule {} (component {}) is not finite at step {} "
          "(cycle {}); the dynamics diverged (check the time step, the relaxation of new bonds and the forces)\n",
          ref.atomIndex, ref.moleculeIndex, ref.componentId, table.currentStep, md.currentCycle));
    }
    return position;
  };

  // ---- the bimolecular rules: every pair in the reactant states within the capture radius -----------------------
  struct Reactant
  {
    std::size_t group;
    std::uint32_t atom;  // the capture atom
    double3 position;
  };
  struct Candidate
  {
    double distance;
    std::uint16_t rule;
    std::size_t a;
    std::size_t b;
    std::uint32_t atomA;
    std::uint32_t atomB;
  };
  std::vector<Candidate> candidates{};
  for (std::size_t r = 0; r < scheme.rules.size(); ++r)
  {
    const ReactionRule& rule = scheme.rules[r];
    if (!rule.bimolecular() || rule.probabilityPerScan(scanTime) <= 0.0) continue;
    std::array<std::vector<Reactant>, 2> reactants{};
    for (std::size_t side = 0; side < (rule.symmetric() ? 1uz : 2uz); ++side)
    {
      for (std::size_t i = 0; i < table.numberOfGroups(); ++i)
      {
        if (table.typeOf(i) != rule.reactantType[side] || table.state(i) != rule.reactantState[side]) continue;
        const std::uint32_t atom = table.atomOf(i, rule.captureAtoms[side].role, components);
        reactants[side].push_back(Reactant{i, atom, positionOf(atom)});
      }
    }
    const double radiusSquared = rule.captureRadius * rule.captureRadius;
    auto consider = [&](const Reactant& x, const Reactant& y)
    {
      if (x.group == y.group) return;
      if (!rule.allowIntramolecular && table.groupRef(x.group).sameMolecule(table.groupRef(y.group))) return;
      const double3 dr = system.simulationBox.applyPeriodicBoundaryConditions(x.position - y.position);
      const double rr = double3::dot(dr, dr);
      if (rr < radiusSquared)
      {
        candidates.push_back(Candidate{std::sqrt(rr), static_cast<std::uint16_t>(r), x.group, y.group, x.atom, y.atom});
      }
    };
    if (rule.symmetric())
    {
      // Both reactants are drawn from the same list. With the same capture role on both sides the pair is
      // unordered (the roles are assigned at random on acceptance); with different roles both assignments are
      // distinct candidates with their own capture distance.
      const std::vector<Reactant>& list = reactants[0];
      const bool sameRole = rule.captureAtoms[0].role == rule.captureAtoms[1].role;
      for (std::size_t i = 0; i < list.size(); ++i)
      {
        for (std::size_t j = i + 1; j < list.size(); ++j)
        {
          if (sameRole)
          {
            consider(list[i], list[j]);
          }
          else
          {
            const std::uint32_t atomJ = table.atomOf(list[j].group, rule.captureAtoms[1].role, components);
            const std::uint32_t atomI = table.atomOf(list[i].group, rule.captureAtoms[1].role, components);
            consider(list[i], Reactant{list[j].group, atomJ, positionOf(atomJ)});
            consider(list[j], Reactant{list[i].group, atomI, positionOf(atomI)});
          }
        }
      }
    }
    else
    {
      for (const Reactant& x : reactants[0])
      {
        for (const Reactant& y : reactants[1]) consider(x, y);
      }
    }
  }
  // closest pairs first; ties broken by the rule and the group indices (a deterministic order for a given random
  // stream)
  std::ranges::sort(candidates, [](const Candidate& x, const Candidate& y)
                    { return std::tie(x.distance, x.rule, x.a, x.b) < std::tie(y.distance, y.rule, y.a, y.b); });

  std::vector<ReactionEvent> events{};
  std::unordered_set<std::size_t> reacted{};

  // The orientation weight of a candidate: exp(-dU / k_B T*) with dU the angular strain of the new bond(s) in the
  // current geometry (the atoms of the products typed as in the products), in Kelvin.
  if (scheme.orientationWeighted() && !candidates.empty() && !table.derivedValid)
  {
    table.rebuild(system.forceField, components, system.numberOfMoleculesPerComponent, atoms);
  }
  auto angularStrainInKelvin = [&](const ReactionRule& rule, std::span<const std::size_t> groups) -> double
  {
    std::map<std::uint32_t, std::uint16_t> retyped{};
    for (std::size_t i = 0; i < groups.size(); ++i)
    {
      for (const auto& [role, type] : rule.retype[i]) retyped[table.atomOf(groups[i], role, components)] = type;
    }
    double strain = 0.0;
    for (const std::array<RoleIndex, 2>& bond : rule.bondAtoms)
    {
      strain += table.angularStrain(table.atomOf(groups[bond[0].reactant], bond[0].role, components),
                                    table.atomOf(groups[bond[1].reactant], bond[1].role, components),
                                    system.simulationBox, components, atoms, retyped);
    }
    return strain * Units::EnergyToKelvin;
  };

  // The configurations the stereocentres of 'rule' take for 'groups': the one required by the tacticity when set
  // (inherited through 'Next'), otherwise the one the geometry gives (the face of approach), or none at all when
  // the tacticity is not controlled; nullopt when a required and a geometric configuration differ and such
  // candidates are skipped.
  auto resolveChiralities = [&](const ReactionRule& rule,
                                std::span<const std::size_t> groups) -> std::optional<std::vector<std::int8_t>>
  {
    std::vector<std::int8_t> result(rule.stereoAtoms.size(), 0);
    for (std::size_t k = 0; k < rule.stereoAtoms.size(); ++k)
    {
      const std::array<std::uint32_t, 4> centreAtoms = table.stereoAtoms(rule, k, groups, components);
      const std::int8_t geometric = table.geometricChirality(centreAtoms, system.simulationBox, atoms);
      const std::int8_t required = table.chirality(groups[rule.stereoAtoms[k].centre.reactant]);
      if (required == 0)
      {
        result[k] = settings.tacticity == Tacticity::None ? std::int8_t{0} : geometric;
        continue;
      }
      if (settings.stereoSelective && required != geometric) return std::nullopt;
      result[k] = required;
    }
    return result;
  };

  for (const Candidate& candidate : candidates)
  {
    if (events.size() >= settings.maximumReactionsPerScan) break;
    if (reacted.contains(candidate.a) || reacted.contains(candidate.b)) continue;
    const ReactionRule& rule = scheme.rules[candidate.rule];

    // a symmetric rule with one capture role: the roles of the pair are assigned at random
    std::array<std::size_t, 2> groups{candidate.a, candidate.b};
    std::array<std::uint32_t, 2> captureAtoms{candidate.atomA, candidate.atomB};
    if (rule.symmetric() && rule.captureAtoms[0].role == rule.captureAtoms[1].role && md.random.uniform() < 0.5)
    {
      std::swap(groups[0], groups[1]);
      std::swap(captureAtoms[0], captureAtoms[1]);
    }

    const std::optional<std::vector<std::int8_t>> chiralities = resolveChiralities(rule, groups);
    if (!chiralities.has_value())
    {
      ++stereoRejections[systemId];
      continue;
    }

    std::optional<double> strain{};
    double weight = 1.0;
    if (rule.orientationTemperature.has_value())
    {
      strain = angularStrainInKelvin(rule, groups);
      weight = rule.orientationWeight(strain.value());
    }
    const double p = rule.probabilityPerScan(scanTime) * weight;
    if (p < 1.0 && md.random.uniform() >= p)
    {
      if (rule.orientationTemperature.has_value()) ++orientationRejections[systemId];
      continue;
    }
    if (strain.has_value())
    {
      acceptedStrain[systemId] += strain.value();
      ++acceptedStrainCount[systemId];
    }

    events.push_back(
        react(systemId, candidate.rule, groups, chiralities.value(), captureAtoms, candidate.distance, strain));
    reacted.insert(candidate.a);
    reacted.insert(candidate.b);
  }

  // ---- the unimolecular rules: every group in the reactant state that did not react in this scan ---------------
  for (std::size_t r = 0; r < scheme.rules.size(); ++r)
  {
    const ReactionRule& rule = scheme.rules[r];
    if (rule.bimolecular()) continue;
    const double p = rule.probabilityPerScan(scanTime);
    if (p <= 0.0) continue;
    for (std::size_t i = 0; i < table.numberOfGroups(); ++i)
    {
      if (events.size() >= settings.maximumReactionsPerScan) break;
      if (table.typeOf(i) != rule.reactantType[0] || table.state(i) != rule.reactantState[0]) continue;
      if (reacted.contains(i)) continue;
      if (p < 1.0 && md.random.uniform() >= p) continue;
      events.push_back(react(systemId, r, std::span<const std::size_t>(&i, 1), {}, {noAtom, noAtom}, 0.0, {}));
      reacted.insert(i);
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
    nlohmann::json groups = nlohmann::json::array();
    for (const ReactiveGroup& group : component.reactiveGroups)
    {
      nlohmann::json atomsJson = nlohmann::json::object();
      for (const auto& [role, atom] : group.atoms) atomsJson[role] = atom;
      groups.push_back({{"type", group.type}, {"initialState", group.initialState}, {"atoms", atomsJson}});
    }
    entry["groups"] = groups;
    componentsJson.push_back(entry);
    ++componentId;
  }
  topology["components"] = componentsJson;
  topology["scheme"] = table.scheme.json();
  topology["tacticity"] = PolymerizationSettings::tacticityName(settings.tacticity);
  nlohmann::json states = nlohmann::json::array();
  nlohmann::json chiralities = nlohmann::json::array();
  for (std::size_t i = 0; i < table.numberOfGroups(); ++i)
  {
    const PolymerAtomRef r = table.groupRef(i);
    states.push_back({r.componentId, r.moleculeIndex, r.atomIndex, table.stateName(i)});
    chiralities.push_back(static_cast<int>(table.chirality(i)));
  }
  topology["groupStates"] = states;
  topology["groupChirality"] = chiralities;
  nlohmann::json centres = nlohmann::json::array();
  for (const StereoCentre& centre : table.centres)
  {
    centres.push_back({{"atoms", centre.atoms},
                       {"chirality", static_cast<int>(centre.chirality)},
                       {"bond", centre.bond},
                       {"group", centre.group}});
  }
  topology["stereoCentres"] = centres;
  topology["bonds"] = table.jsonTopology();
  topology["brokenBonds"] = table.jsonBrokenBonds();

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
    std::vector<std::size_t> counts(table.scheme.rules.size(), 0);
    for (const ReactionEvent& event : events) ++counts[event.rule];
    std::string summary{};
    for (std::size_t r = 0; r < counts.size(); ++r)
    {
      if (counts[r] > 0) summary += std::format("{} {}, ", counts[r], table.scheme.rules[r].name);
    }
    LAMMPS::ExportOptions options{};
    options.dataFile = stem + ".data";
    options.tableFile = stem + ".table";
    options.pairListFile = stem + ".pairs";
    options.comment = std::format(
        "polymerization event frame: step {}, time {:.6g} ps, cycle {}, {}{} polymer bonds in total (see "
        "output/polymerization_events.s{}.jsonl)",
        table.currentStep, static_cast<double>(table.currentStep) * system.timeStep, md.currentCycle, summary,
        table.numberOfBonds(), systemId);
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
  auto groupJson = [&](std::size_t group)
  {
    const PolymerAtomRef r = table.groupRef(group);
    return nlohmann::json{{"index", group},
                          {"component", r.componentId},
                          {"molecule", r.moleculeIndex},
                          {"group", r.atomIndex},
                          {"state", table.stateName(group)}};
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
  nlohmann::json brokenBonds = nlohmann::json::array();
  for (const ReactionEvent& event : events)
  {
    nlohmann::json entry{};
    entry["reaction"] = table.scheme.rules[event.rule].name;
    entry["a"] = groupJson(event.groups[0]);
    entry["atomA"] = atomJson(event.atoms[0]);
    if (event.bimolecular())
    {
      entry["b"] = groupJson(event.groups[1]);
      entry["atomB"] = atomJson(event.atoms[1]);
      entry["distance"] = event.distance;
    }
    entry["chirality"] = event.chirality;
    entry["stereo"] = event.stereo;
    entry["strain"] = optionalNumber(event.strain);
    entry["waitingTime"] = {optionalNumber(event.waitingTime[0]), optionalNumber(event.waitingTime[1])};
    entry["bonds"] = event.bonds;
    entry["brokenBonds"] = event.brokenBonds;
    eventsJson.push_back(entry);
    for (std::size_t bondIndex : event.brokenBonds)
    {
      if (bondIndex >= table.brokenBonds.size()) continue;
      const BrokenBond& bond = table.brokenBonds[bondIndex];
      const std::uint32_t a = static_cast<std::uint32_t>(table.layout.index(bond.a));
      const std::uint32_t b = static_cast<std::uint32_t>(table.layout.index(bond.b));
      nlohmann::json bondJson{};
      bondJson["bond"] = bondIndex;
      bondJson["reaction"] = table.scheme.rules[bond.reaction].name;
      bondJson["atoms"] = {a, b};
      bondJson["lammpsIds"] = {idOffset + a, idOffset + b};
      brokenBonds.push_back(bondJson);
    }
    for (const auto& [g, before, after] : event.retyped)
    {
      retyped.push_back({{"id", g}, {"lammpsId", idOffset + g}, {"from", typeName(before)}, {"to", typeName(after)}});
    }
    for (std::size_t bondIndex : event.bonds)
    {
      if (bondIndex >= table.bonds.size()) continue;
      const PolymerBond& bond = table.bonds[bondIndex];
      nlohmann::json bondJson{};
      bondJson["bond"] = bondIndex;
      bondJson["reaction"] = table.scheme.rules[bond.reaction].name;
      bondJson["atoms"] = {table.layout.index(bond.a), table.layout.index(bond.b)};
      bondJson["formationDistance"] = bond.formationDistance;
      // the bonded terms the bond created (global atoms; valid after the rebuild of the table)
      nlohmann::json bends = nlohmann::json::array();
      nlohmann::json torsions = nlohmann::json::array();
      for (const PolymerLinkTable::LinkBendTerm& term : table.linkBends)
      {
        if (term.bondId == bondIndex) bends.push_back(term.atoms);
      }
      for (const PolymerLinkTable::LinkTorsionTerm& term : table.linkTorsions)
      {
        if (term.bondId == bondIndex) torsions.push_back(term.atoms);
      }
      bondJson["bends"] = bends;
      bondJson["torsions"] = torsions;
      newBonds.push_back(bondJson);
    }
  }
  record["events"] = eventsJson;
  record["retyped"] = retyped;
  record["newBonds"] = newBonds;
  record["brokenBonds"] = brokenBonds;
  record["totals"] = {{"bonds", table.numberOfBonds()},
                      {"brokenBonds", table.numberOfBrokenBonds()},
                      {"groups", table.jsonGroupCounts()},
                      {"active", table.countActive()},
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
  archive >> p.stereoRejections;
  archive >> p.orientationRejections;
  archive >> p.acceptedStrain;
  archive >> p.acceptedStrainCount;
  archive >> p.reactionsStoppedStep;
  archive >> p.eventFrames;
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
