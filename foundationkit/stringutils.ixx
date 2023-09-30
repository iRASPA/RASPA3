export module stringutils;

import <string>;
import <locale>;
import <algorithm>;
import <cctype>;
import <print>;
import <format>;

export namespace std
{
  template<typename... Args> constexpr void print(std::ostream & os, const std::string_view str_fmt, Args&&... args)
  {
    os << std::vformat(str_fmt, std::make_format_args(args...));
  }

  template<typename... Args> constexpr void print(std::ofstream & ofs, const std::string_view str_fmt, Args&&... args)
  {
    ofs << std::vformat(str_fmt, std::make_format_args(args...));
  }
}


export inline bool caseInSensStringCompare(const std::string& str1, const std::string& str2)
{
  return str1.size() == str2.size() && std::equal(str1.begin(), str1.end(), str2.begin(), [](auto a, auto b) {return std::tolower(a) == std::tolower(b); });
}

export inline bool startsWith(const std::string &str, const std::string &prefix) {
    return str.size() >= prefix.size() && str.substr(0, prefix.size()) == prefix;
}

export inline std::string trim(const std::string& s)
{
  auto start = s.begin();
  while (start != s.end() && std::isspace(*start)) {
    start++;
  }

  auto end = s.end();
  do {
    end--;
  } while (std::distance(start, end) > 0 && std::isspace(*end));

  return std::string(start, end + 1);
}
