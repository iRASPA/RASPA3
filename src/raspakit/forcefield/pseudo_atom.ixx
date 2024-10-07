module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <fstream>
#include <iostream>
#include <optional>
#include <ostream>
#include <string>
#include <vector>
#endif

export module pseudo_atom;

#ifndef USE_LEGACY_HEADERS
import <vector>;
import <string>;
import <algorithm>;
import <iostream>;
import <ostream>;
import <fstream>;
import <optional>;
#endif

import archive;
import units;

export struct PseudoAtom
{
  PseudoAtom() {};
  PseudoAtom(std::string name, bool framework, double mass, double charge, double polarizability, size_t atomicNumber,
             bool printToPDB, std::string source = "")
      : name(name),
        framework(framework),
        mass(mass),
        charge(charge),
        polarizability(polarizability),
        atomicNumber(atomicNumber),
        printToPDB(printToPDB),
        source(source) {};
  uint64_t versionNumber{ 1 };

  std::string name{ "C" };
  bool framework{ false };
  double mass{ 1.0 };
  double charge{ 0.0 };
  double polarizability{ 0.0 };
  size_t atomicNumber{ 8 };
  size_t oxidationState{ 0 };
  bool printToPDB{true};
  std::string source{};

  bool operator==(const PseudoAtom &other) const;
  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PseudoAtom &a);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PseudoAtom &a);
};
