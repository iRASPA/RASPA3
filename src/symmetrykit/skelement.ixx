module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <cstdlib>
#include <map>
#include <set>
#include <string>
#include <type_traits>
#include <vector>
#endif

export module skelement;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

export struct SKElement
{
  std::string _chemicalSymbol = std::string("Undefined");
  std::size_t _atomicNumber = 0;
  std::make_signed_t<std::size_t> _group = 0;
  std::size_t _period = 0;
  std::string _name = std::string("Undefined");
  double _mass = 1.0;
  double _atomRadius = 0.0;
  double _covalentRadius = 0.0;
  double _singleBondCovalentRadius = 0.0;
  double _doubleBondCovalentRadius = 0.0;
  double _tripleBondCovalentRadius = 0.0;
  double _VDWRadius = 1.0;
  std::vector<int> _possibleOxidationStates;
  std::size_t _oxidationState = 0;
  double _atomicPolarizability = 0.0;
  SKElement();
  SKElement(std::string string, std::size_t atomicNumber, std::make_signed_t<std::size_t> group, std::size_t period,
            std::string name, double mass, double atomRadius, double covalentRadius, double singleBondCovalentRadius,
            double doubleBondCovalentRadius, double tripleBondCovalentRadius, double vDWRadius,
            std::vector<int> possibleOxidationStates);
};

export struct PredefinedElements
{
  struct InsensitiveCompare
  {
    // case-independent (ci) compare_less binary function
    struct nocase_compare
    {
      bool operator()(const unsigned char& c1, const unsigned char& c2) const { return std::tolower(c1) < std::tolower(c2); }
    };

    bool operator()(const std::string& s1, const std::string& s2) const
    {
      return std::lexicographical_compare(s1.begin(), s1.end(),  // source range
                                          s2.begin(), s2.end(),  // dest range
                                          nocase_compare());     // comparison
    }
  };

  static std::vector<SKElement> predefinedElements;
  static std::map<std::string, std::size_t, PredefinedElements::InsensitiveCompare> atomicNumberData;

  static std::set<std::string> residueDefinitions;
  static std::map<std::string, std::string> residueDefinitionsElement;
  static std::map<std::string, std::string> residueDefinitionsType;
  static std::map<std::string, std::vector<std::string>> residueDefinitionsBondedAtoms;
};
