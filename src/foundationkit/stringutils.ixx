module;

export module stringutils;

import std;
import std.compat;

/// The string without leading and trailing whitespace and with every run of internal whitespace collapsed to one
/// space (the semantics of Qt's QString::simplified); it used to return its argument unchanged.
export inline std::string simplified(const std::string &a)
{
  std::string result;
  result.reserve(a.size());
  bool pendingSpace = false;
  for (unsigned char c : a)
  {
    if (std::isspace(c))
    {
      pendingSpace = !result.empty();
    }
    else
    {
      if (pendingSpace) result.push_back(' ');
      pendingSpace = false;
      result.push_back(static_cast<char>(c));
    }
  }
  return result;
}

export inline std::string toLower(const std::string &a)
{
  std::string data = a;
  std::transform(data.begin(), data.end(), data.begin(),
    [](unsigned char c){ return std::tolower(c); });
  return data;
}



export inline bool caseInSensStringCompare(const std::string& lhs, const std::string& rhs)
{
  // std::tolower takes an int that must be representable as unsigned char: a plain (signed) char above 0x7f,
  // as in UTF-8 text, would be negative
  return lhs.size() == rhs.size() &&
         std::equal(lhs.begin(), lhs.end(), rhs.begin(),
                    [](unsigned char a, unsigned char b) { return std::tolower(a) == std::tolower(b); });
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
  return str.size() >= prefix.size() && str.compare(0, prefix.size(), prefix) == 0;
}

export inline std::string trim(const std::string& s)
{
  const auto isSpace = [](char c) { return std::isspace(static_cast<unsigned char>(c)) != 0; };
  const auto start = std::find_if_not(s.begin(), s.end(), isSpace);
  const auto end = std::find_if_not(s.rbegin(), s.rend(), isSpace).base();
  // an empty or all-whitespace string has nothing left
  return start < end ? std::string(start, end) : std::string();
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
  // a directory exists too (an empty name with RASPA_DIR set resolves to that directory), and a stream over it
  // reads as empty
  if (!std::filesystem::is_regular_file(path))
  {
    throw std::runtime_error(std::format("'{}' is not a file\n", path.string()));
  }

  std::ifstream t(path, std::ios::binary);
  if (!t.is_open())
  {
    throw std::runtime_error(std::format("File '{}' cannot be opened for reading\n", path.string()));
  }
  // sized from the file and read in one call, instead of appended byte by byte through a stream iterator
  std::string file_content(static_cast<std::size_t>(std::filesystem::file_size(path)), '\0');
  t.read(file_content.data(), static_cast<std::streamsize>(file_content.size()));
  file_content.resize(static_cast<std::size_t>(t.gcount()));

  return file_content;
}

