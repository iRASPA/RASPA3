module;

#ifdef USE_LEGACY_HEADERS
#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include <ostream>
#include <fstream>
#include <optional>
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

export struct PseudoAtom
{
  PseudoAtom() {};
  PseudoAtom(std::string name, double mass, double charge, size_t atomicNumber, bool printToPDB, std::string source = ""):
             name(name), mass(mass), charge(charge), atomicNumber(atomicNumber), printToPDB(printToPDB), source(source) {};
  uint64_t versionNumber{ 1 };
  std::string name{ "C" };
  double mass{ 1.0 };
  double charge{ 0.0 };
  size_t atomicNumber{ 8 };
  bool printToPDB{ true };
  std::string source{};

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PseudoAtom &a);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PseudoAtom &a);
};

