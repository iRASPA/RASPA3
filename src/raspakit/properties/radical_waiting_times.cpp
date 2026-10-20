module;

module property_radical_waiting_times;

import std;

import archive;
import polymer_links;

// ---- Histogram ------------------------------------------------------------------------------------------------------

std::size_t PropertyRadicalWaitingTimes::Histogram::binOf(double dt) const
{
  if (!(dt > 0.0)) return 0;
  const double x = (std::log10(dt) - minimumExponent) * static_cast<double>(binsPerDecade);
  if (x < 0.0) return 0;
  return std::min(static_cast<std::size_t>(x), counts.size() - 1);
}

void PropertyRadicalWaitingTimes::Histogram::add(double dt)
{
  ++counts[binOf(dt)];
  ++count;
  sum += dt;
  sumOfLogarithms += std::log(std::max(dt, 1e-300));
  maximum = std::max(maximum, dt);
}

double PropertyRadicalWaitingTimes::Histogram::median() const
{
  if (count == 0) return 0.0;
  std::uint64_t cumulative = 0;
  for (std::size_t bin = 0; bin < counts.size(); ++bin)
  {
    cumulative += counts[bin];
    if (2 * cumulative >= count) return std::sqrt(lowerEdge(bin) * upperEdge(bin));
  }
  return upperEdge(counts.size() - 1);
}

std::uint64_t PropertyRadicalWaitingTimes::Histogram::countAbove(double threshold) const
{
  std::uint64_t total = 0;
  for (std::size_t bin = 0; bin < counts.size(); ++bin)
  {
    if (lowerEdge(bin) >= threshold) total += counts[bin];
  }
  return total;
}

PropertyRadicalWaitingTimes::Histogram &PropertyRadicalWaitingTimes::Histogram::operator+=(const Histogram &other)
{
  for (std::size_t bin = 0; bin < std::min(counts.size(), other.counts.size()); ++bin) counts[bin] += other.counts[bin];
  count += other.count;
  sum += other.sum;
  sumOfLogarithms += other.sumOfLogarithms;
  maximum = std::max(maximum, other.maximum);
  return *this;
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PropertyRadicalWaitingTimes::Histogram &h)
{
  archive << h.binsPerDecade;
  archive << h.counts;
  archive << h.count;
  archive << h.sum;
  archive << h.sumOfLogarithms;
  archive << h.maximum;
  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyRadicalWaitingTimes::Histogram &h)
{
  archive >> h.binsPerDecade;
  archive >> h.counts;
  archive >> h.count;
  archive >> h.sum;
  archive >> h.sumOfLogarithms;
  archive >> h.maximum;
  return archive;
}

// ---- property -------------------------------------------------------------------------------------------------------

PropertyRadicalWaitingTimes::Histogram PropertyRadicalWaitingTimes::closed() const
{
  Histogram all = propagation;
  all += termination;
  return all;
}

PropertyRadicalWaitingTimes::Histogram PropertyRadicalWaitingTimes::open(const PolymerLinkTable &table,
                                                                         double timeStep,
                                                                         std::uint64_t referenceStep) const
{
  Histogram ages(binsPerDecade);
  if (!table.activated) return ages;
  for (std::size_t i = 0; i < table.numberOfSites(); ++i)
  {
    if (table.state(i) != PolymerSiteState::Radical) continue;
    const std::uint64_t since = table.siteStep[i] < referenceStep ? referenceStep - table.siteStep[i] : 0;
    ages.add(static_cast<double>(since) * timeStep);
  }
  return ages;
}

void PropertyRadicalWaitingTimes::writeOutput(std::size_t systemId, std::size_t currentCycle,
                                              const PolymerLinkTable &table, double timeStep,
                                              std::uint64_t referenceStep, bool force) const
{
  if (!force && currentCycle % writeEvery != 0uz) return;
  if (!table.activated) return;

  const Histogram all = closed();
  const Histogram ages = open(table, timeStep, referenceStep);

  std::filesystem::create_directory("radical_waiting_times");
  std::ofstream stream(std::format("radical_waiting_times/radical_waiting_times.s{}.txt", systemId));

  std::print(stream, "# waiting times of the living radicals: dt = t(next reaction) - t(last reaction) of the radical that reacts\n");
  std::print(stream, "# closed intervals end in a propagation (the adding radical) or a termination (both radicals); the open\n");
  std::print(stream, "# intervals are the ages of the radicals still living at the time of writing (from their last reaction,\n");
  std::print(stream, "# or their initiation). A single exponential is a well-defined rate; a stretched tail, a second slow\n");
  std::print(stream, "# component, or a pile-up of old living radicals is the trapping signature.\n");
  std::print(stream, "# step {}, time {:.6g} [ps], conversion {:.6f}, living radicals {}, reactions {}\n", table.currentStep,
             static_cast<double>(table.currentStep) * timeStep, table.conversion(), ages.count,
             referenceStep < table.currentStep ? std::format("stopped at step {}", referenceStep)
                                               : std::string("running"));
  std::print(stream, "# closed: {} (propagation {}, termination {}), mean {:.6g} [ps], geometric mean {:.6g} [ps], max {:.6g} [ps]\n",
             all.count, propagation.count, termination.count, all.mean(), all.geometricMean(), all.maximum);
  std::print(stream, "# open:   {}, mean {:.6g} [ps], geometric mean {:.6g} [ps], max {:.6g} [ps]\n", ages.count,
             ages.mean(), ages.geometricMean(), ages.maximum);
  std::print(stream, "# log-spaced bins ({} per decade); the densities are per decade of time: count / (N * {:.6g})\n",
             binsPerDecade, all.binWidth());
  std::print(stream, "# column 1: lower edge of the bin [ps]\n");
  std::print(stream, "# column 2: upper edge of the bin [ps]\n");
  std::print(stream, "# column 3: closed intervals ended by a propagation [-]\n");
  std::print(stream, "# column 4: closed intervals ended by a termination [-]\n");
  std::print(stream, "# column 5: closed intervals, total [-]\n");
  std::print(stream, "# column 6: open intervals (ages of the living radicals) [-]\n");
  std::print(stream, "# column 7: density of the closed intervals [1/decade]\n");
  std::print(stream, "# column 8: density of the open intervals [1/decade]\n");
  std::print(stream, "# column 9: cumulative distribution of the closed intervals [-]\n");
  std::print(stream, "# column 10: cumulative distribution of the open intervals [-]\n");

  std::uint64_t closedCumulative = 0, openCumulative = 0;
  for (std::size_t bin = 0; bin < all.numberOfBins(); ++bin)
  {
    closedCumulative += all.counts[bin];
    openCumulative += ages.counts[bin];
    if (all.counts[bin] == 0 && ages.counts[bin] == 0) continue;
    std::print(stream, "{:.6g} {:.6g} {} {} {} {} {:.6g} {:.6g} {:.6f} {:.6f}\n", all.lowerEdge(bin), all.upperEdge(bin),
               propagation.counts[bin], termination.counts[bin], all.counts[bin], ages.counts[bin], all.density(bin),
               ages.density(bin),
               all.count > 0 ? static_cast<double>(closedCumulative) / static_cast<double>(all.count) : 0.0,
               ages.count > 0 ? static_cast<double>(openCumulative) / static_cast<double>(ages.count) : 0.0);
  }
}

std::string PropertyRadicalWaitingTimes::printStatus(const PolymerLinkTable &table, double timeStep,
                                                     std::uint64_t referenceStep) const
{
  if (!table.activated) return {};
  const Histogram all = closed();
  const Histogram ages = open(table, timeStep, referenceStep);

  std::string text{};
  text += std::format(
      "    waiting times:                   {} closed intervals (propagation {}, termination {}), mean {:.4g} [ps], "
      "median {:.3g} [ps], max {:.4g} [ps]\n",
      all.count, propagation.count, termination.count, all.mean(), all.median(), all.maximum);
  // the living radicals older than ten times the mean closed interval are the candidates for trapping
  const double threshold = all.count > 0 ? 10.0 * all.mean() : 0.0;
  text += std::format(
      "    living radicals:                 {} with mean age {:.4g} [ps], median {:.3g} [ps], oldest {:.4g} [ps]{}\n",
      ages.count, ages.mean(), ages.median(), ages.maximum,
      all.count > 0 ? std::format(", {} older than 10x the mean waiting time ({:.3g} [ps])", ages.countAbove(threshold),
                                  threshold)
                    : std::string{});
  return text;
}

std::string PropertyRadicalWaitingTimes::printSettings() const
{
  std::ostringstream stream;
  std::print(stream, "Radical waiting times:\n");
  std::print(stream, "    dt = t(next reaction) - t(last reaction) of every living radical, sampled at the reaction\n");
  std::print(stream, "    bins per decade: {} (1e{:g} to 1e{:g} ps)\n", binsPerDecade, minimumExponent, maximumExponent);
  std::print(stream, "    write every: {}\n", writeEvery);
  std::print(stream, "\n");
  return stream.str();
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PropertyRadicalWaitingTimes &p)
{
  archive << p.versionNumber;
  archive << p.binsPerDecade;
  archive << p.writeEvery;
  archive << p.propagation;
  archive << p.termination;
#if DEBUG_ARCHIVE
  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif
  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyRadicalWaitingTimes &p)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > p.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'PropertyRadicalWaitingTimes' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }
  archive >> p.binsPerDecade;
  archive >> p.writeEvery;
  archive >> p.propagation;
  archive >> p.termination;
#if DEBUG_ARCHIVE
  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("PropertyRadicalWaitingTimes: Error in binary restart\n"));
  }
#endif
  return archive;
}
