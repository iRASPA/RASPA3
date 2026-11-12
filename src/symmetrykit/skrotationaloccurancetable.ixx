module;

export module skrotationaloccurancetable;

import std;

export struct SKRotationalOccurrenceTable
{
  // the count of each rotation type, in the slot of the type's value (-6 to 6) plus six
  std::array<std::make_signed_t<std::size_t>, 13> occurrence{};

  SKRotationalOccurrenceTable(std::make_signed_t<std::size_t> axis_6m, std::make_signed_t<std::size_t> axis_4m,
                             std::make_signed_t<std::size_t> axis_3m, std::make_signed_t<std::size_t> axis_2m,
                             std::make_signed_t<std::size_t> axis_1m, std::make_signed_t<std::size_t> axis_1,
                             std::make_signed_t<std::size_t> axis_2, std::make_signed_t<std::size_t> axis_3,
                             std::make_signed_t<std::size_t> axis_4, std::make_signed_t<std::size_t> axis_6);

  inline bool operator==(const SKRotationalOccurrenceTable& b) const { return (this->occurrence == b.occurrence); }
};
