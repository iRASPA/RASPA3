module;

export module rejection_envelope;

import std;

export namespace Potentials
{
/// Identifies the one-dimensional density a rejection sampler draws from: the kind of term (bond, Urey-Bradley,
/// bend), its type and its parameters.
struct PotentialTag
{
  int kind{};
  int type{};
  std::array<double, 8> parameters{};

  PotentialTag(int kindValue, int typeValue, std::span<const double> parameterValues) : kind(kindValue), type(typeValue)
  {
    for (std::size_t i = 0; i < parameterValues.size() && i < parameters.size(); ++i)
      parameters[i] = parameterValues[i];
  }
};

/// The envelope of a rejection sampler (the headroom-scaled maximum of its density over a grid), remembered per
/// thread for the density, the temperature and the interval it was estimated on. The grid estimate costs a
/// thousand evaluations of the potential; the samplers used to repeat it on every single draw, i.e. per bead of
/// every CBMC trial, although it depends on nothing that changes between draws. 'compute' is called on a miss.
template <typename Compute>
double cachedEnvelope(const PotentialTag& tag, double beta, double lower, double upper, Compute&& compute)
{
  struct Entry
  {
    bool used{false};
    PotentialTag tag{0, 0, {}};
    double beta{}, lower{}, upper{};
    double envelope{};
  };
  thread_local std::array<Entry, 64> cache{};
  thread_local std::size_t next = 0;
  for (const Entry& entry : cache)
  {
    if (entry.used && entry.tag.kind == tag.kind && entry.tag.type == tag.type &&
        entry.tag.parameters == tag.parameters && entry.beta == beta && entry.lower == lower && entry.upper == upper)
    {
      return entry.envelope;
    }
  }
  Entry& entry = cache[next];
  next = (next + 1) % cache.size();
  entry = Entry{true, tag, beta, lower, upper, compute()};
  return entry.envelope;
}
}  // namespace Potentials
