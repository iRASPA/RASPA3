module;

module skrotationaloccurancetable;

import std;

import skrotationmatrix;

SKRotationalOccurrenceTable::SKRotationalOccurrenceTable(
    std::make_signed_t<std::size_t> axis_6m, std::make_signed_t<std::size_t> axis_4m,
    std::make_signed_t<std::size_t> axis_3m, std::make_signed_t<std::size_t> axis_2m,
    std::make_signed_t<std::size_t> axis_1m, std::make_signed_t<std::size_t> axis_1,
    std::make_signed_t<std::size_t> axis_2, std::make_signed_t<std::size_t> axis_3,
    std::make_signed_t<std::size_t> axis_4, std::make_signed_t<std::size_t> axis_6)
{
  // the slots of the types -6, -4, -3, -2, -1, 1, 2, 3, 4 and 6; the slots of -5, 0 and 5 stay zero
  occurrence = {axis_6m, 0, axis_4m, axis_3m, axis_2m, axis_1m, 0, axis_1, axis_2, axis_3, axis_4, 0, axis_6};
}
