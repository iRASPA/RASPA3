module;

export module cbmc_util;

import std;

import atom;
import randomnumbers;
import cbmc_grow_step;

export namespace CBMC
{
/// Rosenbluth selection among trial directions given their log Boltzmann factors (-beta U): returns
/// the index of the selected trial, drawn with probability proportional to exp(logBoltzmannFactor).
std::size_t selectTrialPosition(RandomNumber &random, std::span<const double> logBoltzmannFactors);

/// Signed volume of the tetrahedron spanned by the four atoms of a chiral center; its sign is the
/// center's parity.
double chiralSignedVolume(const std::array<std::size_t, 4> &ids, std::span<const Atom> atoms);

/// Writes 'positions' (in 'step.nextBeads' order) into the step's next-beads of the chain.
void placeStepBeads(std::vector<Atom> &chainAtoms, const GrowStep &step, std::span<const Atom> positions);

/// The current atoms of the step's next-beads (in 'step.nextBeads' order): the existing configuration
/// a retrace pins as trial direction 0.
std::vector<Atom> stepBeadPositions(std::span<const Atom> chainAtoms, const GrowStep &step);

/**
 * \brief Scoped scratch use of a few beads of a chain.
 *
 * The operator engine evaluates trial placements by writing the step's next-beads into the chain
 * itself (so every energy routine sees a complete, correctly indexed molecule) instead of copying the
 * whole chain per trial -- for a polymer of N beads that copy made every growth step O(N) and every
 * grow O(N^2). This guard saves the listed beads on construction and puts them back on destruction,
 * on every exit path including exceptions, so the chain is unchanged for the caller. A caller that
 * wants to keep the placed beads (a completed growth level) calls 'keep()'.
 */
class ScratchBeads
{
 public:
  ScratchBeads(std::vector<Atom> &chain, std::span<const std::size_t> beads) : chain_(chain), beads_(beads)
  {
    saved_.reserve(beads.size());
    for (std::size_t bead : beads) saved_.push_back(chain[bead]);
  }
  ScratchBeads(std::vector<Atom> &chain, const GrowStep &step) : ScratchBeads(chain, step.nextBeads) {}
  ~ScratchBeads()
  {
    if (!kept_) restore();
  }
  ScratchBeads(const ScratchBeads &) = delete;
  ScratchBeads &operator=(const ScratchBeads &) = delete;

  /// Puts the saved beads back (also done by the destructor unless kept).
  void restore() const
  {
    for (std::size_t k = 0; k != beads_.size(); ++k) chain_[beads_[k]] = saved_[k];
  }

  /// Leaves the beads as they are now on destruction.
  void keep() noexcept { kept_ = true; }

 private:
  std::vector<Atom> &chain_;
  std::span<const std::size_t> beads_;
  std::vector<Atom> saved_{};
  bool kept_{false};
};
}  // namespace CBMC
