module;

export module blocking_pockets;

import std;

import double4;
import json;
import atom;
import forcefield;
import framework;

// Where a component's blocking pockets come from when they are not written out by hand.
//
// A blocking pocket is a sphere a molecule is refused inside of, and the reason to have one is always the
// same: a cavity of the framework the molecule cannot get out of, and so cannot have got into, which the
// simulation would otherwise fill because insertion is not a journey. The sodalite cages of LTA are the
// standard example. Which cavities those are is a property of the framework and of the probe it is measured
// with, and nothing to do with the run, so it can be answered once from the CIF-file that was just read
// rather than looked up in a table and copied into the input by hand.
//
// That answer is the structural analysis's to give, and this is the engine's side of the question: a
// framework and a force field on the way in, the fractional centres and radii a Component holds on the way
// back. The route is the exact geometric one, which measures each pocket from its own surface and needs no
// grid, no sampling and no accelerator, so it costs the same on every machine and gives the same spheres.
export namespace BlockingPockets
{
/// The size of the nitrogen every structure is measured with, in Ångström. It is the probe the pore analyses
/// report against by default and the one the literature quotes accessible volumes and surface areas at, so a
/// pocket called closed here is a pocket closed by the same standard as everything else said about the
/// framework. A force field written for a simulation carries the guests it simulates and no probes, so this
/// is a length of its own rather than a pseudo-atom to be looked up; where the force field does define
/// 'probe-N2', that definition is preferred, since naming it is how a user overrides this.
inline constexpr double nitrogenProbeSizeParameter = 3.681;

/// The name looked for in the force field before falling back on the size above.
inline constexpr std::string_view nitrogenProbeName = "probe-N2";

/**
 * \struct Specification
 * \brief What an input file asked for: spheres it already knows, and whether more are to be worked out.
 */
struct Specification
{
  std::vector<double4> pockets{};  ///< Spheres listed in the input or read from a file.
  bool automatic{false};           ///< Whether the framework is to be measured once it has been read.
};

/**
 * \brief Reads a 'BlockingPockets' value in any of the forms the input accepts.
 *
 * The value is either a list of spheres, each of them the fractional positions s_x, s_y, s_z and a radius in
 * Ångström; or the string "auto", asking for the spheres to be worked out from the framework; or the name of
 * a `.block` file to read them from.
 *
 * The same value is accepted in the molecule definition file and in the 'Components' section of the
 * simulation input, and what the two say is added together.
 *
 * \param item The JSON value of the 'BlockingPockets' key.
 * \return The spheres it names and whether it also asked for the framework to be measured.
 *
 * \throws std::runtime_error If the value is neither a list of four-element arrays nor a string.
 */
Specification parse(const nlohmann::basic_json<nlohmann::raspa_map> &item);

/**
 * \brief Reads blocking pockets from a `.block` file.
 *
 * The format is the one the structural analysis writes: a count on the first line, then one line per sphere
 * holding the fractional centre s_x, s_y, s_z and a radius in Ångström. Comments are not allowed in it.
 *
 * \param fileName Path to the file; the `.block` extension is added when absent, and the file is looked up
 *        in the working directory and then in RASPA_DIR.
 * \return The spheres as (s_x, s_y, s_z, radius).
 *
 * \throws std::runtime_error If the file is not found, or holds fewer spheres than it says it does.
 */
std::vector<double4> readBlockingPocketFile(const std::string &fileName);

/**
 * \brief Computes the blocking pockets of a framework, as seen by a nitrogen probe.
 *
 * Inflates the framework atoms by the probe, splits the void into what the probe can reach and what it
 * cannot, and covers each unreachable pocket with a sphere at its centroid, of the lesser of the radius that
 * holds the pocket and the radius past which the sphere would reach a channel. Falls back to sampling the
 * void for the structures the surfaces cannot be measured on.
 *
 * The pockets are a property of the framework and of the probe, and so are the same for every component of
 * the system: what they say is that a cavity has no way out, which does not become false for a smaller guest
 * that has no way in either.
 *
 * \param framework The framework, one unit cell of which is measured.
 * \param forceField The force field the atom sizes are read from.
 * \return The spheres as (s_x, s_y, s_z, radius), fractional in the unit cell of the framework.
 *
 * \throws std::runtime_error If the force field gives the framework atoms no size to measure.
 */
std::vector<double4> compute(const Framework &framework, const ForceField &forceField);
}  // namespace BlockingPockets
