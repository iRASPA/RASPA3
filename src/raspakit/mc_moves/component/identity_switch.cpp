module;

module mc_moves_identity_switch;

import std;

import component;
import atom;
import molecule;
import double3;
import simd_quatd;
import simulationbox;
import cbmc;
import cbmc_chain_data;
import cbmc_growth_context;
import cbmc_interactions;
import cbmc_interactions_intermolecular;
import randomnumbers;
import system;
import running_energy;
import forcefield;
import interactions_framework_molecule;
import interactions_ewald;
import interactions_polarization;
import mc_moves_move_types;

std::optional<RunningEnergy> MC_Moves::identitySwitchMove(RandomNumber &random, System &system,
                                                          std::size_t selectedComponent)
{
  std::chrono::steady_clock::time_point time_begin, time_end;
  const Move::Types move = Move::Types::IdentitySwitchCBMC;

  Component &componentAData = system.components[selectedComponent];
  if (componentAData.identitySwitches.empty())
  {
    return std::nullopt;
  }

  const std::size_t componentA = selectedComponent;
  const std::size_t componentB =
      componentAData.identitySwitches[random.uniform_integer(0, componentAData.identitySwitches.size() - 1)];

  if (componentB >= system.components.size() || componentB == componentA)
  {
    return std::nullopt;
  }

  Component &componentBData = system.components[componentB];
  if (componentBData.type != componentAData.type)
  {
    return std::nullopt;
  }

  // Composition is conserved (one A and one B are exchanged), so picking one molecule of each
  // uniformly gives identical proposal probabilities for the forward and reverse moves; no
  // 50/50 role swap is needed for detailed balance.
  if (system.numberOfIntegerMoleculesPerComponent[componentA] == 0 ||
      system.numberOfIntegerMoleculesPerComponent[componentB] == 0)
  {
    return std::nullopt;
  }

  componentAData.mc_moves_statistics.addTrial(move);
  componentBData.mc_moves_statistics.addTrial(move);

  const std::size_t moleculeA = system.randomIntegerMoleculeOfComponent(random, componentA);
  const std::size_t moleculeB = system.randomIntegerMoleculeOfComponent(random, componentB);

  std::span<Atom> atomsA = system.spanOfMolecule(componentA, moleculeA);
  std::span<Atom> atomsB = system.spanOfMolecule(componentB, moleculeB);

  // copies: pristine old configurations (taken before any relabeling below)
  std::vector<Atom> oldA(atomsA.begin(), atomsA.end());
  std::vector<Atom> oldB(atomsB.begin(), atomsB.end());
  const Atom startingBeadA = oldA[componentAData.startingBead];
  const Atom startingBeadB = oldB[componentBData.startingBead];

  const std::size_t globalA = system.moleculeIndexOfComponent(componentA, moleculeA);
  const std::size_t globalB = system.moleculeIndexOfComponent(componentB, moleculeB);
  const std::make_signed_t<std::size_t> skipBothMolecules = static_cast<std::make_signed_t<std::size_t>>(globalA);

  // Determine cutoff distances based on whether dual cutoff is used.
  const double cutOffFrameworkVDW =
      system.forceField.useDualCutOff ? system.forceField.dualCutOff : system.forceField.cutOffFrameworkVDW;
  const double cutOffMoleculeVDW =
      system.forceField.useDualCutOff ? system.forceField.dualCutOff : system.forceField.cutOffMoleculeVDW;
  const double cutOffCoulomb =
      system.forceField.useDualCutOff ? system.forceField.dualCutOff : system.forceField.cutOffCoulomb;

  // Distinct trial ids that cannot collide with any existing molecule (global ids are 0..N-1).
  const std::size_t trialIdNewB = system.numberOfMolecules();
  const std::size_t trialIdNewA = system.numberOfMolecules() + 1;

  const CBMC::GrowContext growContext{system.hasExternalField,
                                      system.forceField,
                                      system.simulationBox,
                                      system.interpolationGrids,
                                      system.externalFieldInterpolationGrid,
                                      system.framework,
                                      system.spanOfFrameworkAtoms(),
                                      system.spanOfMoleculeAtoms(),
                                      system.beta,
                                      cutOffFrameworkVDW,
                                      cutOffMoleculeVDW,
                                      cutOffCoulomb};

  struct CBMCResult
  {
    ChainGrowData growNewB;     // new B grown at A's site
    ChainGrowData growNewA;     // new A grown at B's site
    ChainRetraceData retraceA;  // old A retraced at its site
    ChainRetraceData retraceB;  // old B retraced at its site
  };

  // Both molecules must be invisible to all four CBMC operations (grown/retraced against the
  // background "system minus {A, B}").  The CBMC helpers exclude molecules by global moleculeId:
  // the retraces self-exclude by id-matching and the grows take one skip id.  We therefore
  // temporarily relabel molecule B's atoms with molecule A's global id so that a single skip id
  // (and a single self-match) covers both molecules.  The relabeling is undone unconditionally
  // right after the CBMC phase; nothing inside the window reads moleculeId for any other purpose.
  for (Atom &atom : atomsB)
  {
    atom.moleculeId = static_cast<std::uint32_t>(globalA);
  }

  time_begin = std::chrono::steady_clock::now();
  std::optional<CBMCResult> cbmc = [&]() -> std::optional<CBMCResult>
  {
    std::optional<ChainGrowData> growNewB = CBMC::growMoleculeIdentityChangeInsertion(
        random, growContext, componentBData, componentB, trialIdNewB, startingBeadA, 1.0, 0, false, skipBothMolecules);
    if (!growNewB)
    {
      return std::nullopt;
    }
    if (system.insideBlockedPockets(componentBData,
                                    std::span<const Atom>(growNewB->atoms.begin(), growNewB->atoms.end())))
    {
      return std::nullopt;
    }

    std::optional<ChainGrowData> growNewA = CBMC::growMoleculeIdentityChangeInsertion(
        random, growContext, componentAData, componentA, trialIdNewA, startingBeadB, 1.0, 0, false, skipBothMolecules);
    if (!growNewA)
    {
      return std::nullopt;
    }
    if (system.insideBlockedPockets(componentAData,
                                    std::span<const Atom>(growNewA->atoms.begin(), growNewA->atoms.end())))
    {
      return std::nullopt;
    }

    // Retraces self-exclude by moleculeId; because of the relabeling both molecules are excluded.
    ChainRetraceData retraceA =
        CBMC::retraceMoleculeIdentityChangeDeletion(random, growContext, componentAData, atomsA);
    ChainRetraceData retraceB =
        CBMC::retraceMoleculeIdentityChangeDeletion(random, growContext, componentBData, atomsB);

    return CBMCResult{std::move(*growNewB), std::move(*growNewA), std::move(retraceA), std::move(retraceB)};
  }();

  // restore molecule B's global id unconditionally before anything else touches the system
  for (Atom &atom : atomsB)
  {
    atom.moleculeId = static_cast<std::uint32_t>(globalB);
  }
  time_end = std::chrono::steady_clock::now();
  componentAData.mc_moves_cputime[move][Move::Timing::NonEwald] += (time_end - time_begin);
  system.mc_moves_cputime[move][Move::Timing::NonEwald] += (time_end - time_begin);

  if (!cbmc)
  {
    return std::nullopt;
  }

  if (system.forceField.useDualCutOff)
  {
    // Dual cut-off scheme: correct both grown configurations from the inner cut-off to the full
    // cut-offs, using the same background (both exchanged molecules excluded) as the growth.
    std::optional<RunningEnergy> correctionNewB =
        CBMC::computeDualCutOffCorrection(growContext, componentBData, cbmc->growNewB.atoms, skipBothMolecules);
    std::optional<RunningEnergy> correctionNewA =
        CBMC::computeDualCutOffCorrection(growContext, componentAData, cbmc->growNewA.atoms, skipBothMolecules);
    if (!correctionNewB.has_value() || !correctionNewA.has_value())
    {
      return std::nullopt;
    }

    cbmc->growNewB.energies += correctionNewB.value();
    cbmc->growNewB.RosenbluthWeight *= std::exp(-system.beta * correctionNewB->potentialEnergy());
    cbmc->growNewA.energies += correctionNewA.value();
    cbmc->growNewA.RosenbluthWeight *= std::exp(-system.beta * correctionNewA->potentialEnergy());
  }

  componentAData.mc_moves_statistics.addConstructed(move);
  componentBData.mc_moves_statistics.addConstructed(move);

  // The direct A-B pair interaction is absent from all four CBMC weights (both molecules were
  // excluded from the background), so it enters the acceptance explicitly: new pair minus old pair.
  std::optional<RunningEnergy> newPairEnergy = CBMC::computeInterMolecularEnergy(
      system.forceField, system.simulationBox,
      std::span<const Atom>(cbmc->growNewA.atoms.begin(), cbmc->growNewA.atoms.end()), cutOffMoleculeVDW, cutOffCoulomb,
      std::span<Atom>(cbmc->growNewB.atoms.begin(), cbmc->growNewB.atoms.end()), -1, -1);
  if (!newPairEnergy)
  {
    return std::nullopt;
  }
  std::optional<RunningEnergy> oldPairEnergy = CBMC::computeInterMolecularEnergy(
      system.forceField, system.simulationBox, std::span<const Atom>(oldA.begin(), oldA.end()), cutOffMoleculeVDW,
      cutOffCoulomb, std::span<Atom>(oldB.begin(), oldB.end()), -1, -1);
  if (!oldPairEnergy)
  {
    // an existing configuration should never register as an overlap; reject defensively
    return std::nullopt;
  }

  // Combined Ewald Fourier difference for the simultaneous change of both molecules
  // (single call so the cross terms between the two exchanged molecules are handled exactly).
  std::vector<Atom> newAtoms;
  newAtoms.reserve(cbmc->growNewB.atoms.size() + cbmc->growNewA.atoms.size());
  newAtoms.insert(newAtoms.end(), cbmc->growNewB.atoms.begin(), cbmc->growNewB.atoms.end());
  newAtoms.insert(newAtoms.end(), cbmc->growNewA.atoms.begin(), cbmc->growNewA.atoms.end());

  std::vector<Atom> oldAtoms;
  oldAtoms.reserve(oldA.size() + oldB.size());
  oldAtoms.insert(oldAtoms.end(), oldA.begin(), oldA.end());
  oldAtoms.insert(oldAtoms.end(), oldB.begin(), oldB.end());

  time_begin = std::chrono::steady_clock::now();
  RunningEnergy energyFourierDifference = Interactions::energyDifferenceEwaldFourier(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.trialEik, system.forceField,
      system.simulationBox, newAtoms, oldAtoms, system.netCharge);
  time_end = std::chrono::steady_clock::now();
  componentAData.mc_moves_cputime[move][Move::Timing::Ewald] += (time_end - time_begin);
  system.mc_moves_cputime[move][Move::Timing::Ewald] += (time_end - time_begin);

  // Tail corrections cancel exactly: the multiset of atom types in the box is unchanged
  // (one A and one B on both sides of the move), for both molecule-molecule and
  // framework-molecule contributions.  Likewise the Ewald net-charge correction vanishes.

  RunningEnergy polarizationDifference;
  std::vector<double3> new_electric_field(newAtoms.size());
  std::vector<double3> old_electric_field(oldAtoms.size());
  if (system.forceField.computePolarization)
  {
    Interactions::computeFrameworkMoleculeElectricFieldDifference(system.forceField, system.simulationBox,
                                                                  system.spanOfFrameworkAtoms(), new_electric_field,
                                                                  old_electric_field, newAtoms, oldAtoms);

    Interactions::computeEwaldFourierElectricFieldDifference(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
        system.trialEik, system.forceField, system.simulationBox, new_electric_field, old_electric_field, newAtoms,
        oldAtoms);

    polarizationDifference = Interactions::computePolarizationEnergyDifference(system.forceField, new_electric_field,
                                                                               old_electric_field, newAtoms, oldAtoms);
  }

  const RunningEnergy pairEnergyDifference = newPairEnergy.value() - oldPairEnergy.value();

  // Canonical acceptance: composition is conserved, so fugacities, ideal-gas Rosenbluth weights and
  // the N_old/(N_new+1) factors of the semi-grand move all cancel identically.  What remains is the
  // Rosenbluth ratio of the two grows over the two retraces, corrected for the energy terms that the
  // CBMC weights do not contain (Fourier-space Ewald, polarization, and the explicit pair term).
  const double correctionFactor =
      std::exp(-system.beta * (energyFourierDifference.potentialEnergy() + polarizationDifference.potentialEnergy() +
                               pairEnergyDifference.potentialEnergy()));

  const double acceptanceProbability = correctionFactor *
                                       (cbmc->growNewB.RosenbluthWeight * cbmc->growNewA.RosenbluthWeight) /
                                       (cbmc->retraceA.RosenbluthWeight * cbmc->retraceB.RosenbluthWeight);

  if (random.uniform() < acceptanceProbability)
  {
    componentAData.mc_moves_statistics.addAccepted(move);
    componentBData.mc_moves_statistics.addAccepted(move);

    Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.trialEik);

    std::vector<Atom> acceptedAtomsNewB(cbmc->growNewB.atoms.begin(), cbmc->growNewB.atoms.end());
    for (Atom &atom : acceptedAtomsNewB)
    {
      atom.componentId = static_cast<std::uint8_t>(componentB);
    }
    Molecule acceptedMoleculeNewB = cbmc->growNewB.molecule;
    acceptedMoleculeNewB.componentId = componentB;

    std::vector<Atom> acceptedAtomsNewA(cbmc->growNewA.atoms.begin(), cbmc->growNewA.atoms.end());
    for (Atom &atom : acceptedAtomsNewA)
    {
      atom.componentId = static_cast<std::uint8_t>(componentA);
    }
    Molecule acceptedMoleculeNewA = cbmc->growNewA.molecule;
    acceptedMoleculeNewA.componentId = componentA;

    // Replace A first, then B.  Deleting a molecule of component A shifts atom storage but leaves
    // component B's per-component molecule index unchanged (insertions append at the end of a
    // component's block), so 'moleculeB' still refers to the original molecule; its span, however,
    // must be re-fetched because the underlying atom vector has been modified.
    if (system.forceField.computePolarization)
    {
      std::span<double3> electricFieldNewB = std::span(new_electric_field.begin(), cbmc->growNewB.atoms.size());
      std::span<double3> electricFieldNewA = std::span(
          new_electric_field.begin() + static_cast<std::vector<double3>::difference_type>(cbmc->growNewB.atoms.size()),
          cbmc->growNewA.atoms.size());

      system.deleteMolecule(componentA, moleculeA, atomsA);
      system.insertMoleculePolarization(componentB, acceptedMoleculeNewB, acceptedAtomsNewB, electricFieldNewB);

      std::span<Atom> refetchedAtomsB = system.spanOfMolecule(componentB, moleculeB);
      system.deleteMolecule(componentB, moleculeB, refetchedAtomsB);
      system.insertMoleculePolarization(componentA, acceptedMoleculeNewA, acceptedAtomsNewA, electricFieldNewA);
    }
    else
    {
      system.deleteMolecule(componentA, moleculeA, atomsA);
      system.insertMolecule(componentB, acceptedMoleculeNewB, acceptedAtomsNewB);

      std::span<Atom> refetchedAtomsB = system.spanOfMolecule(componentB, moleculeB);
      system.deleteMolecule(componentB, moleculeB, refetchedAtomsB);
      system.insertMolecule(componentA, acceptedMoleculeNewA, acceptedAtomsNewA);
    }

    const RunningEnergy energyDifference = (cbmc->growNewB.energies + cbmc->growNewA.energies) -
                                           (cbmc->retraceA.energies + cbmc->retraceB.energies) +
                                           energyFourierDifference + pairEnergyDifference + polarizationDifference;

    return energyDifference;
  }

  return std::nullopt;
}
