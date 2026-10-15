module;

export module polymerization_site;

import std;

import archive;

/**
 * \brief A polymerizable group of a molecule (a vinyl group of an acrylate, for example) for the
 * free-radical chain-growth model of the 'Polymerization' driver.
 *
 * Declared per component in the molecule file, separately from the equilibrium 'ReactiveSites' of
 * the cross-link moves:
 *
 *   "PolymerizationSites" : [
 *     {"Type" : "acrylate", "Head" : 0, "Tail" : 1, "HeadProductType" : "CH2_bb", "TailProductType" : "CH_bb"},
 *     {"Type" : "acrylate", "Head" : 13, "Tail" : 12, "HeadProductType" : "CH2_bb", "TailProductType" : "CH_bb"}
 *   ]
 *
 * 'Head' is the atom that a radical attacks (the CH2= of the double bond) and 'Tail' the atom that
 * carries the unpaired electron after the addition (the =CH-); the bond is formed between the tail
 * of the attacking radical and the head of the attacked monomer. When a site reacts its head and tail
 * change pseudo-atom type to the product types (the sp2 beads become the sp3 beads of the backbone),
 * which switches the pair interactions and selects the bonded parameters of the product from the
 * type-based tables of the 'Polymerization' block. The optional 'InitiatedHeadType' is the head type
 * of a site that is activated by an initiator (default: the head product type).
 */
export struct PolymerizationSite
{
  std::uint64_t versionNumber{1};

  std::string siteType{};
  std::size_t head{};
  std::size_t tail{};
  std::string headProductType{};
  std::string tailProductType{};
  std::optional<std::string> initiatedHeadType{};

  bool operator==(PolymerizationSite const &) const = default;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PolymerizationSite &s)
  {
    archive << s.versionNumber;
    archive << s.siteType;
    archive << s.head;
    archive << s.tail;
    archive << s.headProductType;
    archive << s.tailProductType;
    archive << s.initiatedHeadType;
    return archive;
  }
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PolymerizationSite &s)
  {
    std::uint64_t versionNumber;
    archive >> versionNumber;
    if (versionNumber > s.versionNumber)
    {
      throw std::runtime_error(std::format("Invalid version reading 'PolymerizationSite' at line {} in file {}\n",
                                           std::source_location::current().line(),
                                           std::source_location::current().file_name()));
    }
    archive >> s.siteType;
    archive >> s.head;
    archive >> s.tail;
    archive >> s.headProductType;
    archive >> s.tailProductType;
    archive >> s.initiatedHeadType;
    return archive;
  }
};
