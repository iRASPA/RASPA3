module;

export module cbmc_util;

import std;

import atom;
import randomnumbers;

export namespace CBMC
{
/// Rosenbluth selection among trial directions given their log Boltzmann factors (-beta U): returns
/// the index of the selected trial, drawn with probability proportional to exp(logBoltzmannFactor).
std::size_t selectTrialPosition(RandomNumber &random, std::span<const double> logBoltzmannFactors);

/// Signed volume of the tetrahedron spanned by the four atoms of a chiral center; its sign is the
/// center's parity.
double chiralSignedVolume(const std::array<std::size_t, 4> &ids, std::span<const Atom> atoms);

/**
 * \brief Scoped scratch use of a few beads of a chain.
 *
 * The operator engine evaluates trial placements by writing the step's next-beads into the chain
 * itself (so every energy routine sees a complete, correctly indexed molecule) instead of copying the
 * whole chain per trial -- for a polymer of N beads that copy made every growth step O(N) and every
 * grow O(N^2). This guard saves the listed beads on construction and puts them back on destruction,
 * on every exit path including exceptions, so the chain is unchanged for the caller.
 */
class ScratchBeads
{
 public:
  ScratchBeads(std::vector<Atom> &chain, std::span<const std::size_t> beads) : chain_(chain), beads_(beads)
  {
    saved_.reserve(beads.size());
    for (std::size_t bead : beads) saved_.push_back(chain[bead]);
  }
  ~ScratchBeads() { restore(); }
  ScratchBeads(const ScratchBeads &) = delete;
  ScratchBeads &operator=(const ScratchBeads &) = delete;

  /// Puts the saved beads back (also done by the destructor).
  void restore() const
  {
    for (std::size_t k = 0; k != beads_.size(); ++k) chain_[beads_[k]] = saved_[k];
  }

 private:
  std::vector<Atom> &chain_;
  std::span<const std::size_t> beads_;
  std::vector<Atom> saved_{};
};
}  // namespace CBMC
