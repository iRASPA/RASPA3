module;

export module cbmc;

import std;

import atom;
import molecule;
import double3x3;
import double3;
import randomnumbers;
import energy_status;
import energy_status_intra;
import energy_status_inter;
import running_energy;
import cbmc_chain_data;
import component;
import cbmc_util;
export import cbmc_growth_context;

// Error contract of the CBMC entry points below.
//
//  - A grow returns std::nullopt when the trial molecule can not be constructed (every trial of some
//    step overlaps or falls below 'minimumRosenbluthFactor', a recoil-growth dead end, ...). This is
//    an ordinary outcome of the move -- the caller counts it as a rejection.
//  - A retrace has no such outcome for a healthy simulation: the old configuration is an accepted
//    state and always has a weight. (The multiple-first-bead reinsertion retrace returns
//    std::nullopt only when the stored first bead overlaps against the modified background.)
//  - The functions throw std::runtime_error when the simulation state itself is inconsistent: no
//    growth step can be built for the requested placed set, the exact base sampler exhausts its
//    rejection budget, or a recoil-growth retrace finds the existing molecule overlapping (see
//    'retraceRecoilGrowthMoleculeChainDeletion'). No weight is defined in these cases and silently
//    continuing would corrupt the acceptance rule, so the error propagates to the driver, which
//    reports the message and stops. The entry points are therefore deliberately NOT noexcept.
export namespace CBMC
{
  // insertion
  [[nodiscard]] std::optional<ChainGrowData> growMoleculeSwapInsertion(RandomNumber &random, const GrowContext &context,
                                                                       Component &component,
                                                                       std::size_t selectedComponent,
                                                                       std::size_t selectedMolecule, double scaling,
                                                                       std::uint8_t groupId, bool isFractional);
  
  // deletion
  [[nodiscard]] ChainRetraceData retraceMoleculeSwapDeletion(RandomNumber &random, const GrowContext &context,
                                                             const Component &component,
                                                             std::span<Atom> molecule_atoms);
  
  // reinsertion grow
  [[nodiscard]] std::optional<ChainGrowData> growMoleculeReinsertion(RandomNumber &random, const GrowContext &context,
                                                                     Component &component,
                                                                     std::size_t selectedComponent, Molecule &molecule,
                                                                     std::span<Atom> molecule_atoms);
  
  // reinsertion retrace
  [[nodiscard]] std::optional<ChainRetraceData> retraceMoleculeReinsertion(
      RandomNumber &random, const GrowContext &context, const Component &component, Molecule &molecule,
      std::span<Atom> molecule_atoms, double storedR);
  
  // partial reinsertion grow
  [[nodiscard]] std::optional<ChainGrowData> growMoleculePartialReinsertion(
      RandomNumber &random, const GrowContext &context, Component &component, std::size_t selectedComponent,
      Molecule &molecule, std::span<Atom> molecule_atoms,
      const std::vector<std::size_t> &beadsAlreadyPlaced);
  
  // partial reinsertion retrace
  [[nodiscard]] ChainRetraceData retraceMoleculePartialReinsertion(
      RandomNumber &random, const GrowContext &context, const Component &component, Molecule &molecule,
      std::span<Atom> molecule_atoms, const std::vector<std::size_t> &beadsAlreadyPlaced);
  
  // identity change insertion
  [[nodiscard]] std::optional<ChainGrowData> growMoleculeIdentityChangeInsertion(
      RandomNumber &random, const GrowContext &context, Component &component, std::size_t selectedComponent,
      std::size_t selectedMolecule, const Atom &oldStartingBead, double scaling, std::uint8_t groupId,
      bool isFractional, std::make_signed_t<std::size_t> skipBackgroundMolecule = -1);
  
  // identity change deletion
  [[nodiscard]] ChainRetraceData retraceMoleculeIdentityChangeDeletion(RandomNumber &random, const GrowContext &context,
                                                                       const Component &component,
                                                                       std::span<Atom> molecule_atoms);
  
  // distance-biased ion-pair insertion: second molecule with fixed first-bead position
  [[nodiscard]] std::optional<ChainGrowData> growMoleculePairSecondSwapInsertion(
      RandomNumber &random, const GrowContext &context, Component &component, std::size_t selectedComponent,
      std::size_t selectedMolecule, double3 fixedFirstBeadPosition, double scaling, std::uint8_t groupId,
      bool isFractional);
  
  [[nodiscard]] ChainRetraceData retraceMoleculePairSecondSwapDeletion(const GrowContext &context,
                                                                       const Component &component,
                                                                       std::span<Atom> molecule_atoms);
}  // namespace CBMC
