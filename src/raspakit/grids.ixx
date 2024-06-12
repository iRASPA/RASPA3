module;

#ifdef USE_LEGACY_HEADERS
#include <cmath>
#include <cstddef>
#include <istream>
#include <ostream>
#include <fstream>
#include <sstream>
#include <type_traits>
#include <print>
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
import <print>;
#endif


import archive;
import double3;
import stringutils;

import scaling;


export struct Grid
{
};
