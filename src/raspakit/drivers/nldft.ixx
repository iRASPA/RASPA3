module;

export module nldft;

import std;

import uint3;
import system;
import input_reader;
import energy_shared_nldft;

/**
 * \brief Run-control settings for a classical nonlocal DFT (NLDFT) isotherm.
 *
 * The solver itself lives in structurekit (`NLDFTIsotherm`); this driver translates a raspakit
 * `System` into that API and writes the usual output header beside the isotherm file.
 */
export struct NLDFTParameters
{
  /// Energy / density grid size. Zero means use `ForceField::numberOfVDWGridPoints` when set, else 128³.
  uint3 gridSize{0, 0, 0};

  /// Orientations for ρ(r, ω). Zero means auto: 128 for a multi-site (linear) probe, 1 for a sphere.
  std::size_t numberOfOrientations{0};

  /// Prefer the OpenCL energy backend when available (molecular ρ(r, ω) still builds U(r, ω) on the CPU).
  bool useGPU{false};

  /// Include framework electrostatics when the probe carries charge (default: follow `ForceField::useCharge`).
  std::optional<bool> useElectrostatics{};

  /// Relative Ewald precision for the framework potential grid (default: `ForceField::EwaldPrecision`).
  std::optional<double> relativePrecision{};
};

/**
 * \brief Classical nonlocal DFT isotherm driver for a framework + adsorbate probe.
 *
 * Builds the external field V_ext(r) (spherical) or U(r, ω) (linear guest) on the framework energy
 * grid and solves the White Bear FMT + mean-field attraction grand potential, writing the isotherm
 * and Rouquerol/BET summary that `NLDFTIsotherm::run` already produces. Exactly one system with a
 * framework is required; the first component is the probe (built-in linear shapes such as N2, or a
 * single-site / defined-atom fallback).
 */
export struct NLDFT
{
  NLDFT() = delete;
  NLDFT(const NLDFT&) = delete;
  NLDFT& operator=(const NLDFT&) = delete;

  /**
   * \brief Takes ownership of the first system from the input reader.
   */
  explicit NLDFT(InputReader& reader);

  /**
   * \brief Programmatic construction without an input file.
   */
  NLDFT(System system, NLDFTParameters parameters = {});

  System system;
  NLDFTParameters parameters{};
  NLDFTIsotherm result{};

  std::ofstream stream;
  bool outputToFiles{true};

  void run();
  void setup();
  void solve();
  void output();
};
