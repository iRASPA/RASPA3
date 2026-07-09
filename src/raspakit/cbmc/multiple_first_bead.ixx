module;

export module cbmc_multiple_first_bead;

import std;

import atom;
import randomnumbers;
import cbmc_first_bead_data;
import framework;
import component;
import forcefield;
import simulationbox;
import interpolation_energy_grid;
import cbmc_util;
import cbmc_growth_context;

export namespace CBMC
{
[[nodiscard]] std::optional<FirstBeadData> growMoleculeMultipleFirstBeadSwapInsertion(
    RandomNumber &random, const GrowContext &context, const Component &component, const Atom &atom) noexcept;

[[nodiscard]] FirstBeadData retraceMultipleFirstBeadSwapDeletion(RandomNumber &random, const GrowContext &context,
                                                                 const Component &component, const Atom atom) noexcept;

[[nodiscard]] std::optional<FirstBeadData> growMultipleFirstBeadReinsertion(
    RandomNumber &random, const GrowContext &context, const Component &component, const Atom &atom,
    std::make_signed_t<std::size_t> skipBackgroundMolecule = -1) noexcept;

[[nodiscard]] std::optional<FirstBeadData> retraceMultipleFirstBeadReinsertion(
    RandomNumber &random, const GrowContext &context, const Component &component, const Atom &atom, double storedR,
    std::make_signed_t<std::size_t> skipBackgroundMolecule = -1) noexcept;

[[nodiscard]] std::optional<FirstBeadData> growMultipleFirstBeadPartialInsertion(
    const GrowContext &context, const Component &component, const Atom &atom,
    std::make_signed_t<std::size_t> skipBackgroundMolecule = -1) noexcept;

[[nodiscard]] FirstBeadData retraceMultipleFirstBeadPartialDeletion(const GrowContext &context,
                                                                    const Component &component,
                                                                    const Atom &atom) noexcept;

[[nodiscard]] std::optional<FirstBeadData> growFirstBeadAtFixedPosition(const GrowContext &context,
                                                                        const Component &component,
                                                                        const Atom &atom) noexcept;

[[nodiscard]] FirstBeadData retraceFirstBeadAtFixedPosition(const GrowContext &context, const Component &component,
                                                            const Atom atom) noexcept;
}  // namespace CBMC
