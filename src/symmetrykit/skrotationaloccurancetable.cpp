module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <map>
#include <type_traits>
#endif

module skrotationaloccurancetable;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import skrotationmatrix;

SKRotationalOccuranceTable::SKRotationalOccuranceTable(
    std::make_signed_t<std::size_t> axis_6m, std::make_signed_t<std::size_t> axis_4m,
    std::make_signed_t<std::size_t> axis_3m, std::make_signed_t<std::size_t> axis_2m,
    std::make_signed_t<std::size_t> axis_1m, std::make_signed_t<std::size_t> axis_1,
    std::make_signed_t<std::size_t> axis_2, std::make_signed_t<std::size_t> axis_3,
    std::make_signed_t<std::size_t> axis_4, std::make_signed_t<std::size_t> axis_6)
{
  occurance[static_cast<typename std::underlying_type<SKRotationMatrix::RotationType>::type>(
      SKRotationMatrix::RotationType::axis_6m)] = axis_6m;
  occurance[static_cast<typename std::underlying_type<SKRotationMatrix::RotationType>::type>(
      SKRotationMatrix::RotationType::axis_4m)] = axis_4m;
  occurance[static_cast<typename std::underlying_type<SKRotationMatrix::RotationType>::type>(
      SKRotationMatrix::RotationType::axis_3m)] = axis_3m;
  occurance[static_cast<typename std::underlying_type<SKRotationMatrix::RotationType>::type>(
      SKRotationMatrix::RotationType::axis_2m)] = axis_2m;
  occurance[static_cast<typename std::underlying_type<SKRotationMatrix::RotationType>::type>(
      SKRotationMatrix::RotationType::axis_1m)] = axis_1m;
  occurance[static_cast<typename std::underlying_type<SKRotationMatrix::RotationType>::type>(
      SKRotationMatrix::RotationType::axis_1)] = axis_1;
  occurance[static_cast<typename std::underlying_type<SKRotationMatrix::RotationType>::type>(
      SKRotationMatrix::RotationType::axis_2)] = axis_2;
  occurance[static_cast<typename std::underlying_type<SKRotationMatrix::RotationType>::type>(
      SKRotationMatrix::RotationType::axis_3)] = axis_3;
  occurance[static_cast<typename std::underlying_type<SKRotationMatrix::RotationType>::type>(
      SKRotationMatrix::RotationType::axis_4)] = axis_4;
  occurance[static_cast<typename std::underlying_type<SKRotationMatrix::RotationType>::type>(
      SKRotationMatrix::RotationType::axis_6)] = axis_6;
}
