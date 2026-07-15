module;

export module skbandpath;

import std;

import double3;
import double3x3;
import skdefinitions;

// High-symmetry Brillouin-zone path generation for band-structure and phonon-dispersion diagrams.
//
// The labelled k-points and recommended paths follow the crystallographic convention of
//   Y. Hinuma, G. Pizzi, Y. Kumagai, F. Oba, I. Tanaka,
//   "Band structure diagram paths based on crystallography",
//   Comp. Mat. Sci. 128, 140 (2017), DOI:10.1016/j.commatsci.2016.10.015 (the "HPKOT" paper).
//
// The k-point coordinates are expressed as fractional coordinates of the *primitive* reciprocal
// lattice basis of the crystallographic-standard primitive cell (see `primitiveReciprocalCell`).
// This implementation derives the tables directly from the article; it is a stand-alone
// re-implementation (not a port of any existing code base).

export struct SKKVector
{
  std::string label;      // e.g. "GAMMA", "X", "L_2"
  double3 coordinates{};  // fractional coordinates in the primitive reciprocal basis
};

export struct SKBandPath
{
  // Two-letter Bravais-lattice symbol, e.g. "cP", "oF", "mC", "hR", "aP".
  std::string bravaisLattice{};
  // Extended (case-specific) Bravais-lattice symbol from the HPKOT paper, e.g. "cP1", "oF3", "mC2".
  std::string extendedBravaisLattice{};

  bool hasInversionSymmetry{true};
  // True when the path was augmented with -k points (only when the structure has neither
  // inversion symmetry nor assumed time-reversal symmetry).
  bool timeReversalAugmented{false};

  // Labelled high-symmetry points (ordered as defined for the lattice).
  std::vector<SKKVector> kPoints{};
  // Recommended path as consecutive (start-label, end-label) segments. A break in the path is
  // indicated by consecutive segments that do not share an endpoint.
  std::vector<std::pair<std::string, std::string>> segments{};

  // Change-of-basis matrix P converting the conventional to the primitive cell:
  //   (a_p b_p c_p) = (a b c) P          (lattice vectors stored as columns)
  double3x3 primitiveTransformationMatrix{};
  // Linear map Q expressing the primitive cell in the conventional-cell basis: (a_p b_p c_p) = (a b c) Q.
  // Equals primitiveTransformationMatrix for all lattices except aP, where the cell is re-standardized.
  double3x3 conventionalToPrimitiveTransformation{};
  // Linear map expressing the analyzed cell (the unit cell passed to the structure-based overload) in the
  // conventional-cell basis: (analyzed cell) = (a b c) * conventionalToAnalyzedTransformation. Identity for
  // the parameter-based overload.
  double3x3 conventionalToAnalyzedTransformation{};
  // Crystallographic-standard conventional cell actually used (columns are lattice vectors).
  // For the triclinic (aP) case this is the re-standardized cell.
  double3x3 conventionalCell{};
  // Crystallographic-standard primitive cell (columns are lattice vectors), = conventionalCell * P.
  double3x3 primitiveCell{};
  // Reciprocal lattice of the primitive cell (columns are reciprocal vectors, includes the 2*pi factor),
  //   primitiveReciprocalCell = 2*pi * (primitiveCell^{-1})^T.
  // A labelled point with fractional coordinates f has Cartesian wavevector primitiveReciprocalCell * f.
  double3x3 primitiveReciprocalCell{};

  // Convenience lookup of a labelled point by name.
  std::optional<double3> coordinatesForLabel(std::string_view label) const;

  // Convert a point given as fractional coordinates in the primitive reciprocal basis into fractional
  // coordinates in the reciprocal basis of the analyzed cell (the cell passed to the structure-based
  // overload; for the parameter-based overload this is the conventional cell).
  double3 reciprocalFractionalInAnalyzedCell(const double3& primitiveReciprocalFractional) const;
};

// Return the HPKOT high-symmetry k-point path for a crystal, given its crystallographic-standard
// conventional cell (columns are lattice vectors), its space-group number (1..230) and whether its
// point group is centrosymmetric.
//
// The conventional cell must already be in the crystallographic-standard setting (e.g. as produced by
// SKSpaceGroup standardization); the axial ratios and the space-group number determine the extended
// Bravais-lattice case and hence the labelled points.
//
// `withTimeReversal == false` together with a non-centrosymmetric point group augments the path with
// the corresponding -k segments.
export std::optional<SKBandPath> SKBrillouinZonePath(double3x3 conventionalCell, std::size_t spaceGroupNumber,
                                                     bool hasInversionSymmetry, bool withTimeReversal = true,
                                                     double threshold = 1.0e-7);

// Convenience overload: detect the space group of the given structure and build the path.
// `unitCell` columns are lattice vectors; `atoms` are (fractional position, type, occupancy).
export std::optional<SKBandPath> SKBrillouinZonePath(
    double3x3 unitCell, std::vector<std::tuple<double3, std::size_t, double>> atoms, bool allowPartialOccupancies,
    double symmetryPrecision, bool withTimeReversal = true, double threshold = 1.0e-7);
