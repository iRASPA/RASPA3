module;

export module finite_value;

import std;

/// Whether 'value' is a finite number (neither NaN nor infinite), read off its exponent bits through a volatile
/// copy. Every kit is compiled with -ffast-math, under which std::isfinite is true of everything, a comparison with
/// a NaN is whatever the optimiser makes of it, and even an integer test of the bit pattern of a value the compiler
/// holds to be finite is folded away; the volatile read is not.
export inline bool isFiniteValue(double value)
{
  volatile std::uint64_t bits = std::bit_cast<std::uint64_t>(value);
  const std::uint64_t copy = bits;
  return ((copy >> 52) & 0x7FFu) != 0x7FFu;
}
