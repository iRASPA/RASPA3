module;

export module skrotationaloccurancetable;

import std;

export struct SKRotationalOccurrenceTable
{
  std::unordered_map<std::make_signed_t<std::size_t>, std::make_signed_t<std::size_t>> occurrence{};

  SKRotationalOccurrenceTable(std::make_signed_t<std::size_t> axis_6m, std::make_signed_t<std::size_t> axis_4m,
                             std::make_signed_t<std::size_t> axis_3m, std::make_signed_t<std::size_t> axis_2m,
                             std::make_signed_t<std::size_t> axis_1m, std::make_signed_t<std::size_t> axis_1,
                             std::make_signed_t<std::size_t> axis_2, std::make_signed_t<std::size_t> axis_3,
                             std::make_signed_t<std::size_t> axis_4, std::make_signed_t<std::size_t> axis_6);

  inline bool operator==(const SKRotationalOccurrenceTable& b) const { return (this->occurrence == b.occurrence); }
};
