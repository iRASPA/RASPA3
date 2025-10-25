module;

#ifdef USE_PRECOMPILED_HEADERS
#include "pch.h"
#endif

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <unordered_map>
#include <type_traits>
#endif

export module skrotationaloccurancetable;

#ifdef USE_STD_IMPORT
import std;
#endif

export struct SKRotationalOccuranceTable
{
  std::unordered_map<std::make_signed_t<std::size_t>, std::make_signed_t<std::size_t>> occurance{};

  SKRotationalOccuranceTable(std::make_signed_t<std::size_t> axis_6m, std::make_signed_t<std::size_t> axis_4m,
                             std::make_signed_t<std::size_t> axis_3m, std::make_signed_t<std::size_t> axis_2m,
                             std::make_signed_t<std::size_t> axis_1m, std::make_signed_t<std::size_t> axis_1,
                             std::make_signed_t<std::size_t> axis_2, std::make_signed_t<std::size_t> axis_3,
                             std::make_signed_t<std::size_t> axis_4, std::make_signed_t<std::size_t> axis_6);

  inline bool operator==(const SKRotationalOccuranceTable& b) const { return (this->occurance == b.occurance); }
};
