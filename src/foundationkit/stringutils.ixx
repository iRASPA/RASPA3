module;

#ifdef USE_PRECOMPILED_HEADERS
#include "pch.h"
#endif

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <cctype>
#include <cstddef>
#include <format>
#include <locale>
#include <print>
#include <string>
#include <type_traits>
#endif

export module stringutils;

#ifdef USE_STD_IMPORT
import std;
import std.compat;
#endif

export inline std::string simplified(std::string a) { return a; };

export inline std::string toLower(const std::string &a)
{
  std::string data = a;
  std::transform(data.begin(), data.end(), data.begin(),
    [](unsigned char c){ return std::tolower(c); });
  return data;
}



export inline bool caseInSensStringCompare(const std::string& lhs, const std::string& rhs)
{
  return lhs.size() == rhs.size() && std::equal(lhs.begin(), lhs.end(), rhs.begin(),
                                                [](auto a, auto b) { return std::tolower(a) == std::tolower(b); });
}

export struct caseInsensitiveComparator
{
  struct nocase_compare
  {
    bool operator()(const unsigned char& c1, const unsigned char& c2) const { return tolower(c1) < tolower(c2); }
  };
  bool operator()(const std::string& lhs, const std::string& rhs) const
  {
    return std::lexicographical_compare(lhs.begin(), lhs.end(), rhs.begin(), rhs.end(), nocase_compare());
  }
};

export inline bool startsWith(const std::string& str, const std::string& prefix)
{
  return str.size() >= prefix.size() && str.substr(0, prefix.size()) == prefix;
}

export inline std::string trim(const std::string& s)
{
  auto start = s.begin();
  while (start != s.end() && std::isspace(*start))
  {
    start++;
  }

  auto end = s.end();
  do
  {
    end--;
  } while (std::distance(start, end) > 0 && std::isspace(*end));

  return std::string(start, end + 1);
}

export inline std::string addExtension(const std::string& fileName, const std::string& extension)
{
  if (fileName.length() >= extension.length() &&
      fileName.compare(fileName.length() - extension.length(), extension.length(), extension) == 0)
  {
    return fileName;
  }
  else
  {
    return fileName + extension;
  }
}
