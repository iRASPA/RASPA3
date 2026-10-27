module;

export module mc_moves_multi_trial_common;

import std;

import atom;
import running_energy;
import system;
import cbmc_multi_trial;

// Shared by the multiple-try moves (swap, reinsertion, Gibbs swap): the grower they run on and the
// long-range (Ewald Fourier + tail) evaluation they fold into every trial's weight.
export namespace MC_Moves
{
/// The grower of the system: the device grower when the system has a device environment, the CPU
/// reference grower otherwise.
[[nodiscard]] const CBMC::MultiTrialGrower& multiTrialGrower(const System& system);

/// The long-range contributions of exchanging 'removed' for 'inserted' in 'system': the Ewald Fourier
/// difference (which also leaves the trial structure factor of the exchange in 'system.trialEik') and,
/// when 'includeTail' is set, the tail-correction difference. Returned separately so the move can add them
/// to the running energies on acceptance.
struct LongRangeDifference
{
  RunningEnergy fourier;
  RunningEnergy tail;
  [[nodiscard]] double potentialEnergy() const { return fourier.potentialEnergy() + tail.potentialEnergy(); }
  [[nodiscard]] RunningEnergy total() const { return fourier + tail; }
};
[[nodiscard]] LongRangeDifference longRangeDifference(System& system, std::span<const Atom> inserted,
                                                      std::span<const Atom> removed, bool includeTail);

/// The correction callback of a trial set that exchanges 'removed' for a trial in 'system': the long-range
/// energy dU(trial, removed) - 'baseline'. Every weight of a move must be measured from one common state:
/// the insertion and the reinsertion measure from the current state ('baseline' = 0; the old configuration
/// of a reinsertion then evaluates to 0), the deletion from the state without the molecule
/// ('baseline' = dU({}, removed); the molecule itself then evaluates to -baseline, its own contribution).
/// The existing molecule is recognised by identity of the span and evaluates to -baseline exactly.
[[nodiscard]] CBMC::LongRangeCorrection longRangeCorrection(System& system, std::span<const Atom> removed,
                                                            bool includeTail, double baseline);

/// Throws when the force field asks for something the multiple-try moves do not support (polarization).
void checkMultiTrialSupported(const System& system, std::string_view moveName);
}  // namespace MC_Moves
