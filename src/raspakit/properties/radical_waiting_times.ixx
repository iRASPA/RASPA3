module;

export module property_radical_waiting_times;

import std;

import archive;
import polymer_links;

/**
 * \brief Waiting-time distribution of the active groups (the living chain ends) of a polymerization
 * ('ComputeRadicalWaitingTimes').
 *
 * For every group in an active state of the ReactionScheme (a radical of a chain growth, a partly reacted
 * amine of a step growth) the time between its successive reactions, dt = t(next reaction) - t(last reaction)
 * (for the first reaction: since it entered the active state), is collected when it reacts, in a log-spaced
 * histogram from 1e-3 to 1e6 ps with 'binsPerDecade' bins per decade; separately per reaction rule that closes
 * the interval. The open intervals, the ages of the groups that are still active (from their last reaction up to
 * the current step, or up to the step at which the reactions stopped), are built from the PolymerLinkTable when
 * the property reports.
 *
 * A single-exponential waiting time is a well-defined propagation rate; a stretched tail, a second slow
 * component, or a pile-up of old active groups that no longer react is the signature of radical trapping. The
 * property measures, it does not impose: the reaction scan treats every active group as reactive.
 *
 * The samples are added by the Polymerization driver at the reaction (the only moment the closed interval is
 * known), not by a sampling interval. Output: 'radical_waiting_times/radical_waiting_times.s<i>.txt' every
 * 'writeEvery' cycles and at the end of the run.
 */
export struct PropertyRadicalWaitingTimes
{
  static constexpr double minimumExponent{-3.0};
  static constexpr double maximumExponent{6.0};

  /// A log-spaced histogram of times [ps] with its running sums.
  struct Histogram
  {
    std::size_t binsPerDecade{10};
    std::vector<std::uint64_t> counts{};
    std::uint64_t count{0};
    double sum{0.0};
    double sumOfLogarithms{0.0};
    double maximum{0.0};

    explicit Histogram(std::size_t bins = 10)
        : binsPerDecade(bins), counts(static_cast<std::size_t>(maximumExponent - minimumExponent) * bins, 0)
    {
    }

    std::size_t numberOfBins() const { return counts.size(); }
    double lowerEdge(std::size_t bin) const
    {
      return std::pow(10.0, minimumExponent + static_cast<double>(bin) / static_cast<double>(binsPerDecade));
    }
    double upperEdge(std::size_t bin) const { return lowerEdge(bin + 1); }
    /// The width of a bin in decades of time.
    double binWidth() const { return 1.0 / static_cast<double>(binsPerDecade); }
    std::size_t binOf(double dt) const;
    void add(double dt);
    double mean() const { return count > 0 ? sum / static_cast<double>(count) : 0.0; }
    double geometricMean() const
    {
      return count > 0 ? std::exp(sumOfLogarithms / static_cast<double>(count)) : 0.0;
    }
    /// The median from the histogram (the geometric centre of the bin containing it).
    double median() const;
    /// The number of values in bins at or above 'threshold'.
    std::uint64_t countAbove(double threshold) const;
    /// The probability density per decade of time of a bin.
    double density(std::size_t bin) const
    {
      return count > 0 ? static_cast<double>(counts[bin]) / (static_cast<double>(count) * binWidth()) : 0.0;
    }
    Histogram &operator+=(const Histogram &other);

    friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const Histogram &h);
    friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, Histogram &h);
  };

  PropertyRadicalWaitingTimes() {}
  PropertyRadicalWaitingTimes(std::size_t binsPerDecade, std::size_t writeEvery)
      : binsPerDecade(binsPerDecade), writeEvery(writeEvery)
  {
  }

  std::uint64_t versionNumber{2};

  std::size_t binsPerDecade{10};
  std::size_t writeEvery{5000};
  std::vector<Histogram> byReaction{};  ///< Closed intervals, per reaction rule that ended them.

  /// Adds the waiting time [ps] of an active group that reacts now through rule 'reaction'.
  void add(std::size_t reaction, double dt);

  /// The closed intervals of all reactions together.
  Histogram closed() const;
  /// The open intervals: the ages [ps] of the active groups of 'table' at 'referenceStep'.
  Histogram open(const PolymerLinkTable &table, double timeStep, std::uint64_t referenceStep) const;

  /// Writes the output file when 'currentCycle' is a multiple of 'writeEvery' (or always with 'force').
  void writeOutput(std::size_t systemId, std::size_t currentCycle, const PolymerLinkTable &table, double timeStep,
                   std::uint64_t referenceStep, bool force = false) const;

  /// The two lines of the status report (the closed intervals and the ages of the active groups).
  std::string printStatus(const PolymerLinkTable &table, double timeStep, std::uint64_t referenceStep) const;
  std::string printSettings() const;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PropertyRadicalWaitingTimes &p);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyRadicalWaitingTimes &p);
};
