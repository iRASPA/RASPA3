module;

export module molecule_property_settings;

import std;

import archive;

// Per-component settings of the intra-molecular analyses. These are component options in the input
// ('Components' block): a component that carries one of these asks for the analysis, with its own
// sampling schedule and histogram sizing. The accumulators themselves live in the System (the
// 'Property...' classes), which read the settings of every component when they are constructed.

/// Histograms of the intra-molecular geometry (bonds, bends, torsions, end-to-end distance);
/// 'ComputeMoleculeProperties' in the input.
export struct MoleculePropertiesSettings
{
  std::uint64_t versionNumber{2};

  std::size_t sampleEvery{10};
  std::optional<std::size_t> writeEvery{5000};
  std::size_t numberOfBins{128};  ///< Bins of the bond, bend and torsion histograms (fixed ranges).
  double bondRange{4.0};          ///< Upper limit of the bond-length histogram [Angstrom].
  std::optional<double>
      endToEndRange{};  ///< Upper limit of the end-to-end histogram [Angstrom]; default: contour length.
  /// Bin width of the end-to-end histogram [Angstrom]; its range grows with the chain length, so the
  /// number of bins is ceil(range / width) rather than a fixed count.
  double endToEndBinWidth{0.25};

  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const MoleculePropertiesSettings& s);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, MoleculePropertiesSettings& s);
};

/// Gyration-tensor shape descriptors; 'ComputeMoleculeShape' in the input.
export struct MoleculeShapeSettings
{
  std::uint64_t versionNumber{2};

  std::size_t sampleEvery{10};
  std::optional<std::size_t> writeEvery{5000};
  std::size_t numberOfBins{128};  ///< Bins of the shape-anisotropy and prolateness histograms (fixed ranges).
  bool massWeighted{false};       ///< Mass weights instead of uniform bead weights.
  std::optional<double>
      radiusOfGyrationRange{};  ///< Upper limit of the Rg histogram [Angstrom]; default from the contour length.
  /// Bin width of the Rg histogram [Angstrom]; the number of bins is ceil(range / width).
  double radiusOfGyrationBinWidth{0.1};

  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const MoleculeShapeSettings& s);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, MoleculeShapeSettings& s);
};

/// Chain statistics along the backbone and the single-chain form factor; 'ComputeMoleculeBackbone'.
export struct MoleculeBackboneSettings
{
  std::uint64_t versionNumber{1};

  std::size_t sampleEvery{10};
  std::optional<std::size_t> writeEvery{5000};
  std::size_t numberOfWaveVectors{64};
  double waveVectorLowerLimit{0.01};  ///< [1/Angstrom]
  double waveVectorUpperLimit{5.0};   ///< [1/Angstrom]

  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const MoleculeBackboneSettings& s);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, MoleculeBackboneSettings& s);
};

/// End-to-end vector autocorrelation function (order-N); 'ComputeEndToEndACF'.
export struct EndToEndACFSettings
{
  std::uint64_t versionNumber{1};

  std::size_t sampleEvery{10};
  std::optional<std::size_t> writeEvery{5000};
  std::size_t numberOfBlockElements{25};  ///< Elements n per block of the order-N scheme.

  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const EndToEndACFSettings& s);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, EndToEndACFSettings& s);
};
