module;

#ifdef USE_LEGACY_HEADERS
#include <optional>
#include <span>
#include <vector>
#endif

export module cbmc_flexible_insertion;

#ifndef USE_LEGACY_HEADERS
import <vector>;
import <optional>;
import <span>;
#endif

import atom;
import double3x3;
import double3;
import randomnumbers;
import cbmc_chain_data;
import component;

// export [[nodiscard]] std::optional<ChainData>
// growFlexibleMoleculeSwapInsertion(const Component &component, RandomNumber &random, double cutOff,
//                                   double cutOffCoulomb, size_t selectedComponent, size_t selectedMolecule,
//                                   double scaling, std::vector<Atom> atoms) noexcept;
