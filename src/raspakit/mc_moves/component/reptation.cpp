module;

module mc_moves_reptation;

import std;

import component;
import atom;
import molecule;
import double3;
import simulationbox;
import cbmc;
import randomnumbers;
import system;
import running_energy;
import forcefield;
import interactions_framework_molecule;
import interactions_intermolecular;
import interactions_ewald;
import interactions_external_field;
import interactions_polarization;
import mc_moves_move_types;

std::optional<RunningEnergy> MC_Moves::reptationMove(RandomNumber &random, System &system,
                                                     std::size_t selectedComponent, std::size_t selectedMolecule)
{
  std::span<Atom> molecule_atoms = system.spanOfMolecule(selectedComponent, selectedMolecule);
  Molecule &molecule = system.moleculeData[system.moleculeIndexOfComponent(selectedComponent, selectedMolecule)];

  std::chrono::steady_clock::time_point time_begin, time_end;
  Move::Types move = Move::Types::Reptation;
  Component &component = system.components[selectedComponent];

  // Restriction: the molecule must declare at least two validated repeat units ('RepeatUnits' in
  // the molecule JSON). This does not change during a simulation, so the attempts are not counted
  // as trials.
  const std::vector<std::vector<std::size_t>> &units = component.repeatUnits;
  if (units.size() < 2)
  {
    return std::nullopt;
  }

  component.mc_moves_statistics.addTrial(move);


  const std::size_t numberOfUnits = units.size();
  const std::size_t numberOfBeads = molecule_atoms.size();

  // Choose the reptation direction with equal probability: forward removes the head unit
  // (units[0]) and grows a new unit at the tail, backward is the mirror image. The reverse of a
  // forward step is a backward step, so the direction choice keeps the proposal symmetric.
  bool forward = random.uniform() < 0.5;
  const std::vector<std::size_t> &vacatedUnit = forward ? units.front() : units.back();
  const std::vector<std::size_t> &grownUnit = forward ? units.back() : units.front();

  // Placed sets for the CBMC growth plans (sorted, so the cached plans are shared between calls):
  // the retrace sees the current molecule with the vacated unit removable, the grow sees the
  // shifted molecule with the arriving unit's slots empty.
  std::vector<bool> isVacated(numberOfBeads, false);
  for (std::size_t atom : vacatedUnit) isVacated[atom] = true;
  std::vector<bool> isGrown(numberOfBeads, false);
  for (std::size_t atom : grownUnit) isGrown[atom] = true;
  std::vector<std::size_t> placedForRetrace{};
  std::vector<std::size_t> placedForGrow{};
  placedForRetrace.reserve(numberOfBeads - vacatedUnit.size());
  placedForGrow.reserve(numberOfBeads - grownUnit.size());
  for (std::size_t atom = 0; atom != numberOfBeads; ++atom)
  {
    if (!isVacated[atom]) placedForRetrace.push_back(atom);
    if (!isGrown[atom]) placedForGrow.push_back(atom);
  }

  // Construct the shifted fixed configuration: the surviving units move one block toward the
  // vacated end, so that after the move the molecule is again in canonical labeling. The
  // shift-periodicity validated at parse time guarantees that every atom keeps its type, charge
  // and topological role under this relabeling. The arriving unit's slots hold stale positions;
  // they are regrown below.
  std::vector<Atom> shiftedAtoms(molecule_atoms.begin(), molecule_atoms.end());
  if (forward)
  {
    for (std::size_t k = 1; k != numberOfUnits; ++k)
    {
      for (std::size_t j = 0; j != units[k].size(); ++j)
      {
        shiftedAtoms[units[k - 1][j]].position = molecule_atoms[units[k][j]].position;
      }
    }
  }
  else
  {
    for (std::size_t k = numberOfUnits - 1; k != 0; --k)
    {
      for (std::size_t j = 0; j != units[k].size(); ++j)
      {
        shiftedAtoms[units[k][j]].position = molecule_atoms[units[k - 1][j]].position;
      }
    }
  }

  // Grow the arriving unit attached to the shifted chain.
  const CBMC::GrowContext context = system.makeGrowContext();

  time_begin = std::chrono::steady_clock::now();
  std::optional<CBMC::GrowResult> growData = CBMC::regrowMolecule(
      random, context, component, molecule, shiftedAtoms,
      {.firstBead = CBMC::FirstBeadScheme::AlreadyPlaced, .beadsAlreadyPlaced = placedForGrow});
  time_end = std::chrono::steady_clock::now();
  component.mc_moves_cputime[move][Move::Timing::NonEwald] += (time_end - time_begin);
  system.mc_moves_cputime[move][Move::Timing::NonEwald] += (time_end - time_begin);

  if (!growData) return std::nullopt;

  std::span<const Atom> newMolecule = std::span(growData->atoms.begin(), growData->atoms.end());

  std::vector<Atom> old_molecule = std::vector(molecule_atoms.begin(), molecule_atoms.end());
  std::vector<double3> old_electric_field = std::vector<double3>(old_molecule.size());
  std::vector<double3> new_electric_field = std::vector<double3>(old_molecule.size());

  component.mc_moves_statistics.addConstructed(move);

  // Retrace the departing unit in the current configuration.
  time_begin = std::chrono::steady_clock::now();
  CBMC::RetraceResult retraceData = CBMC::retraceMolecule(
      random, context, component, molecule_atoms,
      {.firstBead = CBMC::FirstBeadScheme::AlreadyPlaced, .beadsAlreadyPlaced = placedForRetrace});
  time_end = std::chrono::steady_clock::now();
  component.mc_moves_cputime[move][Move::Timing::NonEwald] += (time_end - time_begin);
  system.mc_moves_cputime[move][Move::Timing::NonEwald] += (time_end - time_begin);

  // Compute the energy difference in the Fourier space due to Ewald summation. The relabeling of
  // the surviving units does not change the structure factors (identical types and charges per
  // slot), so the difference stems from the vacated and the grown unit only.
  time_begin = std::chrono::steady_clock::now();
  RunningEnergy energyFourierDifference = Interactions::energyDifferenceEwaldFourier(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.trialEik, system.forceField,
      system.simulationBox, newMolecule, molecule_atoms);
  time_end = std::chrono::steady_clock::now();
  component.mc_moves_cputime[move][Move::Timing::Ewald] += (time_end - time_begin);
  system.mc_moves_cputime[move][Move::Timing::Ewald] += (time_end - time_begin);

  // Dual cut-off scheme: correct the grown and retraced configurations from the inner cut-off to the
  // full cut-offs.
  if (!CBMC::applyDualCutOffCorrection(context, component, *growData) ||
      !CBMC::applyDualCutOffCorrection(context, component, old_molecule, retraceData))
  {
    return std::nullopt;
  }

  std::vector<double3> electricFieldNeighborDelta;
  RunningEnergy polarizationDifference;
  if (system.forceField.computePolarization)
  {
    Interactions::computeFrameworkMoleculeElectricFieldDifference(system.forceField, system.simulationBox,
                                                                  system.spanOfFrameworkAtoms(), new_electric_field,
                                                                  old_electric_field, growData->atoms, old_molecule);

    Interactions::computeEwaldFourierElectricFieldDifference(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
        system.trialEik, system.forceField, system.simulationBox, new_electric_field, old_electric_field,
        growData->atoms, old_molecule);

    // Molecule-molecule polarization (inter-molecular energy already handled via CBMC; discard returned energy).
    if (!system.forceField.omitInterPolarization)
    {
      electricFieldNeighborDelta.assign(system.spanOfMoleculeAtoms().size(), double3(0.0, 0.0, 0.0));
      [[maybe_unused]] std::optional<RunningEnergy> interPolarizationEnergy =
          Interactions::computeInterMolecularPolarizationElectricFieldDifference(
              system.forceField, system.simulationBox, electricFieldNeighborDelta, new_electric_field,
              old_electric_field, system.spanOfMoleculeAtoms(), growData->atoms, old_molecule);
    }

    polarizationDifference = Interactions::computePolarizationEnergyDifference(
        system.forceField, new_electric_field, old_electric_field, growData->atoms, old_molecule);

    if (!system.forceField.omitInterPolarization)
    {
      polarizationDifference += Interactions::computePolarizationEnergyNeighborDifference(
          system.forceField, system.spanOfMoleculeElectricField(), electricFieldNeighborDelta,
          system.spanOfMoleculeAtoms());
    }
  }

  // The grow and retrace plans belong to opposite chain ends and may partition the bonded terms
  // differently between the base samplers and the Rosenbluth weights (a branched unit grows its
  // side group as a sibling of the backbone at one end but sequentially at the other). The exact
  // base samplers make the base densities Boltzmann, but their normalizations are plan-dependent,
  // and the weight ratio below is only a valid acceptance up to the ratio of those normalizations:
  // correct by exp(logZ_growPlan - logZ_retracePlan). (For same-plan moves such as reinsertion this
  // factor is identically one.)
  const std::vector<CBMC::GrowStep> &growPlan = component.growthPlan(placedForGrow);
  const std::vector<CBMC::GrowStep> &retracePlan = component.growthPlan(placedForRetrace);
  const double logBaseNormalizationCorrection =
      CBMC::logBaseSamplerNormalization(system.beta, component, growPlan) -
      CBMC::logBaseSamplerNormalization(system.beta, component, retracePlan);

  // Metropolis acceptance with the configurational-bias weight ratio evaluated in log space.
  double logAcceptance =
      -system.beta * (energyFourierDifference.potentialEnergy() + polarizationDifference.potentialEnergy()) +
      growData->logRosenbluthWeight - retraceData.logRosenbluthWeight + logBaseNormalizationCorrection;

  if (random.uniform() < std::exp(logAcceptance))
  {
    component.mc_moves_statistics.addAccepted(move);

    Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.trialEik);
    std::copy(newMolecule.begin(), newMolecule.end(), molecule_atoms.begin());

    if (system.forceField.computePolarization && !system.forceField.omitInterPolarization)
    {
      std::span<double3> storedElectricField = system.spanOfMoleculeElectricField();
      for (std::size_t i = 0; i < storedElectricField.size(); ++i)
      {
        storedElectricField[i] += electricFieldNeighborDelta[i];
      }
    }

    if (system.forceField.computePolarization)
    {
      std::span<double3> electricFieldMolecule = system.spanElectricFieldOld(selectedComponent, selectedMolecule);
      std::copy(new_electric_field.begin(), new_electric_field.end(), electricFieldMolecule.begin());
    }

    molecule = growData->molecule;
    molecule.centerOfMassPosition = component.computeCenterOfMass(molecule_atoms);

    return (growData->energies - retraceData.energies) + energyFourierDifference + polarizationDifference;
  }

  return std::nullopt;
}
