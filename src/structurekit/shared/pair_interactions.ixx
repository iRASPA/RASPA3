module;

export module pair_interactions;

import std;

// What one atom of the structure does to one atom of the probe, for every pair of types that can meet, and
// nothing else that a force field happens to know.
//
// A force field is a large object with a mixing rule, a charge method, Ewald parameters, tail corrections and
// a temperature, and none of that is a question about a pore. The routines here ask three numbers of a pair
// and two cutoffs of the whole, so what they are given is those numbers, already mixed, already shifted, in a
// square table indexed by type. The mixing rule stays where it belongs, in whatever built this, and is
// applied once rather than being reapplied here in a way that might not agree.
//
// Everything is in the internal units of the caller. Sizes are in Ångström, strengths and shifts in whatever
// energy unit the report headers are later written in.
export struct PairParameters
{
  double sizeParameter{0.0};      // σ, Å
  double strengthParameter{0.0};  // ε
  double shift{0.0};              // what is subtracted so the potential reaches the cutoff at zero
};

export struct PairInteractions
{
  std::size_t numberOfTypes{0};

  // Per type, for the report headers and for the charge a probe site carries.
  std::vector<std::string> names;
  std::vector<double> charges;

  // Row-major over (type, type), so `numberOfTypes * numberOfTypes` long. Symmetric in practice, but stored
  // in full because that is what an indexed lookup in an inner loop wants.
  std::vector<PairParameters> parameters;

  double cutOffVDW{12.0};
  double cutOffCoulomb{12.0};

  // The pair of two types, and the pair of a type with itself. The second is what a probe's own size and
  // strength are read from, before any mixing.
  const PairParameters& operator()(std::size_t row, std::size_t column) const
  {
    return parameters[row * numberOfTypes + column];
  }
  const PairParameters& operator[](std::size_t row) const { return parameters[row * numberOfTypes + row]; }

  // The index of a named type, or nothing when the name is not one of them. This is how a probe named on a
  // command line becomes a row of the table.
  std::optional<std::size_t> findType(const std::string& name) const;
};
