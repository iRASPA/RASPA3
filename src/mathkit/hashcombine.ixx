module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#endif

export module hashcombine;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

// https://stackoverflow.com/questions/35985960/c-why-is-boosthash-combine-the-best-way-to-combine-hash-values/50978188

export inline void hash_combine([[maybe_unused]] std::size_t& seed) {}

export template <typename T, typename... Rest>
inline void hash_combine(std::size_t& seed, const T& v, Rest... rest)
{
  std::hash<T> hasher;
  seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
  hash_combine(seed, rest...);
}
