export module print;

import <iostream>;
import <fstream>;
import <string_view>;
import <string>;
import <cstdio>;
import <utility>;
import <version>;
import <format>;


#if defined(WIN32) || defined(_WIN32) || defined(_WIN64) || defined(__WIN32__) || defined(__NT__)
    import <format>;
    #define BWP_FMT_LIB "std"
    #define BWP_FMTNS std
#elif defined (__cpp_lib_format)
    #include <format>
    #define BWP_FMT_LIB "std"
    #define BWP_FMTNS std
#else
    #define FMT_HEADER_ONLY
    #include "fmt/core.h"
    #include "fmt/format.h"
    #include "fmt/xchar.h"
    #include "fmt/ostream.h"
    #include "fmt/printf.h"
    #include "fmt/chrono.h"
    #define BWP_FMT_LIB "libfmt"
    #define BWP_FMTNS fmt
#endif // __cpp_lib_format

export namespace std 
{
    template<typename... Args> std::string print(const char* str_fmt, Args&&... args) {
        return BWP_FMTNS::vformat(str_fmt, BWP_FMTNS::make_format_args(args...));
    }

    // default to stdout
    //template<typename... Args> constexpr void print(const std::string_view str_fmt, Args&&... args) {
    //    fputs(std::vformat(str_fmt, std::make_format_args(args...)).c_str(), stdout);
    //}

    
    // send to FILE*
    template<typename... Args> constexpr void print(FILE* fdest, const std::string_view str_fmt, Args&&... args) {
        fputs(BWP_FMTNS::vformat(str_fmt, BWP_FMTNS::make_format_args(args...)).c_str(), fdest);
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
