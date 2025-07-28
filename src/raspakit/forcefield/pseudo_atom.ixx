module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <optional>
#include <ostream>
#include <string>
#include <vector>
#endif

export module pseudo_atom;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import archive;
import units;

/**
 * \brief Represents a pseudo-atom in the simulation system.
 *
 * The PseudoAtom struct encapsulates the properties of an atom used in simulations,
 * including name, mass, charge, polarizability, atomic number, oxidation state, and other attributes.
 * It provides constructors for initializing pseudo-atoms and methods for serialization and comparison.
 */
export struct PseudoAtom
{
  /**
   * \brief Default constructor for the PseudoAtom struct.
   *
   * Initializes a PseudoAtom object with default values.
   */
  PseudoAtom() {};

  /**
   * \brief Constructs a PseudoAtom with specified parameters.
   *
   * Initializes a PseudoAtom with the provided name, framework flag, mass, charge,
   * polarizability, atomic number, printToPDB flag, and source.
   *
   * \param name The name of the pseudo-atom.
   * \param framework Indicates if the atom is part of the framework.
   * \param mass The mass of the pseudo-atom.
   * \param charge The electric charge of the pseudo-atom.
   * \param polarizability The polarizability of the pseudo-atom.
   * \param atomicNumber The atomic number of the pseudo-atom.
   * \param printToPDB Flag indicating whether to include the atom in PDB output.
   * \param source The source of the pseudo-atom data.
   */
  PseudoAtom(std::string name, bool framework, double mass, double charge, double polarizability,
             std::size_t atomicNumber, bool printToPDB, std::string source = "")
      : name(name),
        framework(framework),
        mass(mass),
        charge(charge),
        polarizability(polarizability),
        atomicNumber(atomicNumber),
        printToPDB(printToPDB),
        source(source) {};

  std::uint64_t versionNumber{1};  ///< Version number for serialization compatibility.

  std::string name{"C"};           ///< The name of the pseudo-atom.
  bool framework{false};           ///< Indicates if the atom is part of the framework.
  double mass{1.0};                ///< The mass of the pseudo-atom.
  double charge{0.0};              ///< The electric charge of the pseudo-atom.
  double polarizability{0.0};      ///< The polarizability of the pseudo-atom.
  std::size_t atomicNumber{8};     ///< The atomic number of the pseudo-atom.
  std::int64_t oxidationState{0};  ///< The oxidation state of the pseudo-atom.
  bool printToPDB{true};           ///< Flag indicating whether to include the atom in PDB output.
  std::string source{};            ///< The source of the pseudo-atom data.

  bool operator==(const PseudoAtom &other) const;
  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PseudoAtom &a);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PseudoAtom &a);
};
