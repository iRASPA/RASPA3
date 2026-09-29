module;

export module stringutils;

import std;
import std.compat;

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

/// Maximum number of columns of a line in the output file.
export inline constexpr std::size_t outputLineWidth = 120;

/// Number of display columns of a UTF-8 string (counts code points, not bytes).
export inline std::size_t displayWidth(std::string_view s)
{
  std::size_t width = 0;
  for (unsigned char c : s)
  {
    if ((c & 0xC0) != 0x80) ++width;
  }
  return width;
}

/// Word-wraps `text` so that every emitted line, including its indent, fits within `width` columns.
/// The first line is prefixed with `firstIndent`, subsequent lines with `continuationIndent`.
/// Words are separated on white space, except that a bracketed unit such as "[K]" stays attached to the
/// value preceding it. A single word longer than the available width is emitted on its own line.
/// Every emitted line is terminated with a newline.
export inline std::string wrapText(std::string_view text, std::string_view firstIndent,
                                   std::string_view continuationIndent, std::size_t width = outputLineWidth)
{
  std::vector<std::string> words;
  std::size_t position = 0;
  while (position < text.size())
  {
    std::size_t begin = text.find_first_not_of(" \t\n\r", position);
    if (begin == std::string_view::npos) break;
    std::size_t end = text.find_first_of(" \t\n\r", begin);
    if (end == std::string_view::npos) end = text.size();
    std::string_view word = text.substr(begin, end - begin);
    if (!words.empty() && word.starts_with('['))
    {
      words.back() += ' ';
      words.back() += word;
    }
    else
    {
      words.emplace_back(word);
    }
    position = end;
  }

  std::string result;
  std::string line(firstIndent);
  std::size_t lineWidth = displayWidth(firstIndent);
  bool lineHasWord = false;
  for (const std::string& word : words)
  {
    std::size_t wordWidth = displayWidth(word);
    if (lineHasWord && lineWidth + 1 + wordWidth > width)
    {
      result += line;
      result += '\n';
      line = std::string(continuationIndent);
      lineWidth = displayWidth(continuationIndent);
      lineHasWord = false;
    }
    if (lineHasWord)
    {
      line += ' ';
      ++lineWidth;
    }
    line += word;
    lineWidth += wordWidth;
    lineHasWord = true;
  }
  result += line;
  result += '\n';
  return result;
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

export std::string readFileContent(const std::string &fileName, const std::string &extension)
{
  std::string file_name_string = addExtension(fileName, extension);

  std::filesystem::path path = std::filesystem::path(file_name_string);
  if (!std::filesystem::exists(path))
  {
    if (const char* env_p = std::getenv("RASPA_DIR"))
    {
      path = std::filesystem::path(env_p) / file_name_string;
    }
  }

  if (!std::filesystem::exists(path))
  {
    throw std::runtime_error(
        std::format("File '{}' not found (also not in 'RASPA_DIR')\n", file_name_string));
  }

  std::ifstream t(path);
  std::string file_content((std::istreambuf_iterator<char>(t)), std::istreambuf_iterator<char>());

  return file_content;
}

