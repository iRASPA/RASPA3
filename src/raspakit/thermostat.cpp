module;

#ifdef USE_LEGACY_HEADERS
#if defined(__has_include) && __has_include(<format>)
#include <format>
#endif
#include <vector>
#include <tuple>
#if defined(__has_include) && __has_include(<print>)
  #include <print>
#endif
#endif

module thermostat;

#ifndef USE_LEGACY_HEADERS
import <vector>;
import <tuple>;
#if defined(__has_include) && __has_include(<print>)
  import <print>;
#endif
#endif

#if !(defined(__has_include) && __has_include(<print>))
  import print;
#endif

import archive;


