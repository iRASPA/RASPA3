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
import units;
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
  stereoRejections.assign(md.systems.size(), 0);
  orientationRejections.assign(md.systems.size(), 0);
  acceptedStrain.assign(md.systems.size(), 0.0);
  acceptedStrainCount.assign(md.systems.size(), 0);

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
                 "conversion radicals bonds stereo strain[K]\n"
                 "# stereo: configuration (+/-) of the reacting tail(s) and the dyad with the previous unit (m meso, "
                 "r racemo, . none); '-' no stereocentre\n"
                 "# strain: angular strain of the new bond at the reaction (with an 'OrientationTemperature'); '-' not "
                 "evaluated\n");
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
      std::flush(stream);
      ++systemId;
    }
    writeTopology();
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
         orientationText(systemId);
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
  for (std::size_t n = 0; n < count; ++n)
  {
    const std::size_t pick = n + static_cast<std::size_t>(md.random.uniform() * static_cast<double>(monomers.size() - n));
    std::swap(monomers[n], monomers[std::min(pick, monomers.size() - 1)]);
    table.initiate(monomers[n], system.forceField, system.components, system.spanOfMoleculeAtoms(),
                   std::span<Molecule>(system.moleculeData), table.currentStep);
    if (settings.writeBondLog && systemId < bondLogs.size())
    {
      const PolymerAtomRef r = table.siteRef(monomers[n]);
      std::print(bondLogs[systemId], "{} {} initiation {} {} {} - - - 0 {:.6f} {} {} - -\n", table.currentStep,
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
  auto log = [&](std::string_view event, std::size_t a, std::size_t b, double distance, std::string_view stereo)
  {
    if (!settings.writeBondLog || systemId >= bondLogs.size()) return;
    const PolymerAtomRef ra = table.siteRef(a);
    const PolymerAtomRef rb = table.siteRef(b);
    std::print(bondLogs[systemId], "{} {} {} {} {} {} {} {} {} {:.4f} {:.6f} {} {} {} {}\n", table.currentStep,
               md.currentCycle, event, ra.componentId, ra.moleculeIndex, ra.atomIndex, rb.componentId,
               rb.moleculeIndex, rb.atomIndex, distance, table.conversion(),
               table.countSites(PolymerSiteState::Radical), table.numberOfBonds(), stereo,
               weighted ? std::format("{:.1f}", lastStrain) : std::string("-"));
  };
  for (const Candidate& candidate : candidates)
  {
    if (accepted >= settings.maximumReactionsPerScan) break;
    if (reacted.contains(candidate.a) || reacted.contains(candidate.b)) continue;
    const double probability =
        candidate.kind == PolymerBondKind::Propagation ? probabilityPropagation : probabilityTermination;
    if (probability <= 0.0) continue;

    if (candidate.kind == PolymerBondKind::Propagation)
    {
      const std::optional<std::int8_t> chiralityA =
          resolveChirality(candidate.a, table.headOf(candidate.b, components));
      if (!chiralityA.has_value())
      {
        ++stereoRejections[systemId];
        continue;
      }
      if (!accept(candidate, probability)) continue;
      const std::string stereo = stereoText(candidate.a, chiralityA.value());
      table.propagate(candidate.a, candidate.b, system.forceField, components, atoms,
                      std::span<Molecule>(system.moleculeData), table.currentStep, candidate.distance,
                      chiralityA.value(), nextChirality(chiralityA.value()));
      log("propagation", candidate.a, candidate.b, candidate.distance, stereo);
    }
    else
    {
      const std::optional<std::int8_t> chiralityA =
          resolveChirality(candidate.a, table.tailOf(candidate.b, components));
      const std::optional<std::int8_t> chiralityB =
          resolveChirality(candidate.b, table.tailOf(candidate.a, components));
      if (!chiralityA.has_value() || !chiralityB.has_value())
      {
        ++stereoRejections[systemId];
        continue;
      }
      if (!accept(candidate, probability)) continue;
      const std::string stereo = std::format("{}/{}", stereoText(candidate.a, chiralityA.value()),
                                             stereoText(candidate.b, chiralityB.value()));
      table.terminate(candidate.a, candidate.b, system.forceField, components, atoms,
                      std::span<Molecule>(system.moleculeData), table.currentStep, candidate.distance,
                      chiralityA.value(), chiralityB.value());
      log("termination", candidate.a, candidate.b, candidate.distance, stereo);
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
        std::print(bondLogs[systemId], "{} {} initiation {} {} {} - - - 0 {:.6f} {} {} - -\n", table.currentStep,
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
  archive << p.stereoRejections;
  archive << p.orientationRejections;
  archive << p.acceptedStrain;
  archive << p.acceptedStrainCount;
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
