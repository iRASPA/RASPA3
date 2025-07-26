module;

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

#ifndef USE_LEGACY_HEADERS
import std;
#endif

export inline bool caseInSensStringCompare(const std::string& str1, const std::string& str2)
{
  return str1.size() == str2.size() && std::equal(str1.begin(), str1.end(), str2.begin(),
                                                  [](auto a, auto b) { return std::tolower(a) == std::tolower(b); });
}

export struct caseInsensitiveComparator
{
  bool operator()(const std::string &lhs, const std::string &rhs) const
  {
#if defined(_WIN32)
    return _stricmp(lhs.c_str(), rhs.c_str()) < 0;
#else
    return strcasecmp(lhs.c_str(), rhs.c_str()) < 0;
#endif
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
