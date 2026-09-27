module;

export module lammps_reader;

import std;

import json;

/**
 * \brief LAMMPS -> RASPA: reads a 'read_data' file (plus, optionally, the input script that names the styles)
 * and produces RASPA's force_field.json, one component JSON per distinct molecule template, and a
 * simulation.json skeleton.
 *
 * What is converted:
 *   - Masses, 'Pair Coeffs' / 'PairIJ Coeffs' (or pair_coeff lines of the script) -> pseudo-atoms and
 *     Lennard-Jones self / binary interactions. Atoms of one LAMMPS type that carry different charges
 *     become distinct RASPA pseudo-atoms (RASPA charges live on the pseudo-atom).
 *   - Bond / Angle / Dihedral Coeffs -> RASPA bond, bend and torsion definitions (index based), through
 *     LAMMPS::bondFromLammps etc. Styles come from the input script, from 'hybrid' prefixes in the
 *     coefficient lines, or default to harmonic / harmonic / nharmonic.
 *   - Molecules (grouped by molecule id) are compared by atom types, charges and topology; each distinct
 *     template becomes a component whose reference geometry is the first instance (unwrapped).
 *   - special_bonds 1-4 factors -> Intra14VanDerWaalsScalingValue / Intra14ChargeChargeScalingValue.
 *   - pair_style cut-offs, pair_modify mix/shift/tail, kspace_style -> force-field and simulation settings.
 *
 * Coordinates of all molecules are kept in 'positions' (Angstrom, unwrapped) so a caller can seed a RASPA
 * configuration; RASPA itself starts from the component templates.
 */
export namespace LAMMPS
{
struct ReadComponent
{
  std::string name{};
  nlohmann::json definition{};
  std::size_t count{};
  std::size_t atomsPerMolecule{};
};

struct ReadResult
{
  nlohmann::json forceField{};
  std::vector<ReadComponent> components{};
  nlohmann::json simulation{};
  std::array<double, 3> boxLengths{};
  std::array<double, 3> boxAngles{90.0, 90.0, 90.0};
  std::vector<std::array<double, 3>> positions{};  ///< all atoms in file order, unwrapped, Angstrom
  std::vector<std::string> warnings{};
};

ReadResult readDataFile(const std::filesystem::path &dataFile, const std::optional<std::filesystem::path> &inputScript);

/// Writes force_field.json, <component>.json files and simulation.json into 'directory'.
void writeRaspaInput(const ReadResult &result, const std::filesystem::path &directory);
}  // namespace LAMMPS
