module;

#ifdef USE_LEGACY_HEADERS
#include <iostream>
#include <fstream>
#include <string_view>
#include <string>
#include <cstdio>
#include <utility>
#include <version>
#if defined(__has_include) && __has_include(<format>)
#include <format>
#endif
#endif

#if defined(__has_include) && __has_include(<print>) && __has_include(<format>)
    #include <format>
    #define BWP_FMT_LIB "std"
    #define BWP_FMTNS std
#else
    #include <utility>
    #define FMT_HEADER_ONLY
    #include "fmt/core.h"
    #include "fmt/format.h"
    //#include "fmt/xchar.h"
    #include "fmt/ostream.h"
    #include "fmt/printf.h"
    #include "fmt/chrono.h"
    #define BWP_FMT_LIB "libfmt"
    #define BWP_FMTNS fmt
#endif // __cpp_lib_format

export module print;

#ifndef USE_LEGACY_HEADERS
import <iostream>;
import <fstream>;
import <string_view>;
import <string>;
import <cstdio>;
import <utility>;
import <version>;
import <format>;
#endif


export namespace std 
{
  template<typename... Args> std::string format(const char* str_fmt, Args&&... args) {
      return BWP_FMTNS::vformat(str_fmt, BWP_FMTNS::make_format_args(args...));
  }

  // send to FILE*
  template<typename... Args> constexpr void print(FILE* fdest, const std::string_view str_fmt, Args&&... args) {
      try
      {
         fputs(BWP_FMTNS::vformat(str_fmt, BWP_FMTNS::make_format_args(args...)).c_str(), fdest);
      }
      catch (std::exception const& e)
      {
          std::cerr << e.what();
      }
  }

  // send to ostream
  template<typename... Args> constexpr void print(std::ostream & ostream_dest, const std::string_view str_fmt, Args&&... args) {
      ostream_dest << BWP_FMTNS::vformat(str_fmt, BWP_FMTNS::make_format_args(args...));
  }

  // send to ofstream
  template<typename... Args> constexpr void print(std::ofstream & ofstream_dest, const std::string_view str_fmt, Args&&... args) {
      ofstream_dest << BWP_FMTNS::vformat(str_fmt, BWP_FMTNS::make_format_args(args...));
  }
}
