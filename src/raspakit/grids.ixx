module;

#ifdef USE_LEGACY_HEADERS
#include <cmath>
#include <cstddef>
#include <istream>
#include <ostream>
#include <fstream>
#include <sstream>
#include <type_traits>
#if defined(__has_include) && __has_include(<print>)
  #include <print>
#endif
#endif

#if defined(_WIN32)
  import <cassert>;
#else
  #include <assert.h>
#endif

export module grids;

#ifndef USE_LEGACY_HEADERS
import <cmath>;
import <cstddef>;
import <istream>;
import <ostream>;
import <fstream>;
import <sstream>;
import <type_traits>;
#if defined(__has_include) && __has_include(<print>)
  import <print>;
#endif
#endif

#if !(defined(__has_include) && __has_include(<print>))
  import print;
#endif

import archive;
import double3;
import stringutils;

import scaling;


export struct Grid
{
};
