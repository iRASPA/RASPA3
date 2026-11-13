module;

module scanner;

import std;

import characterset;

Scanner::Scanner(const std::string& content, const CharacterSet& charactersToBeSkipped)
    : _charactersToBeSkipped(charactersToBeSkipped)
{
  _string = content;
  _scanLocation = _string.cbegin();
}

std::string::const_iterator Scanner::find_first_not_of(const std::string& chars, const std::string& text,
                                                       std::string::const_iterator location)
{
  std::string::const_iterator it = location;
  std::string::const_iterator end = text.cend();
  while (it != end)
  {
    if (chars.contains(*it) == false) return it;
    ++it;
  }

  return text.cend();
}

std::string::const_iterator Scanner::find_first_of(const std::string& chars, const std::string& text,
                                                   std::string::const_iterator location)
{
  std::string::const_iterator it = location;
  std::string::const_iterator end = text.cend();
  while (it != end)
  {
    if (chars.contains(*it)) return it;
    ++it;
  }

  return text.cend();
}

bool Scanner::scanCharacters(const CharacterSet& set, std::string& into)
{
  std::string::const_iterator found = find_first_not_of(set.string(), _string, _scanLocation);

  // the characters of the set up to the next other character, or up to the end of the string; none at all (the
  // next character is not in the set, or the string is done) is no match, as for the other scans
  if (found == _scanLocation || _scanLocation >= _string.cend())
  {
    into = std::string("");
    return false;
  }
  into = std::string(_scanLocation, found);
  _scanLocation = found;
  return true;
}

bool Scanner::scanLine(std::string& into)
{
  // the set is the same for every line: made once
  static const CharacterSet newline = CharacterSet::newlineCharacter();
  std::string::const_iterator newlineLocation =
      find_first_of(newline.string(), _string, _scanLocation);

  // nothing left
  if (_scanLocation >= _string.cend())
  {
    into = std::string("");
    return false;
  }

  // the last line need not end with a newline
  if (newlineLocation >= _string.cend())
  {
    into = std::string(_scanLocation, _string.cend());
    _scanLocation = _string.cend();
    return true;
  }

  into = std::string(_scanLocation, newlineLocation);
  _scanLocation = newlineLocation + 1;
  return true;
}

bool Scanner::scanUpToCharacters(const CharacterSet& set, std::string& into)
{
  std::string::const_iterator found = find_first_not_of(_charactersToBeSkipped.string(), _string, _scanLocation);

  if (found >= _string.cend())
  {
    _scanLocation = _string.cend();
    into = std::string("");
    return false;
  }

  _scanLocation = found;

  // up to the next character of the set, or up to the end of the string (the last token need not be followed by
  // a delimiter)
  found = find_first_of(set.string(), _string, _scanLocation);
  into = std::string(_scanLocation, found);
  _scanLocation = found;
  return true;
}

bool Scanner::isAtEnd() { return _scanLocation >= _string.cend(); }

bool Scanner::scanDouble(double& value)
{
  std::string::const_iterator found = find_first_not_of(_charactersToBeSkipped.string(), _string, _scanLocation);

  if (found >= _string.cend())
  {
    _scanLocation = _string.cend();
    return false;
  }

  _scanLocation = found;

  // the token runs up to the next skipped character or to the end of the string
  found = find_first_of(_charactersToBeSkipped.string(), _string, _scanLocation);
  // parsed in place; the library parser is the fallback for the forms 'from_chars' does not take (a leading plus,
  // a leading space), so what was accepted before is accepted still
  const char* first = std::to_address(_scanLocation);
  const char* last = std::to_address(found);
  _scanLocation = found;
  if (const std::from_chars_result result = std::from_chars(first, last, value);
      result.ec == std::errc{} && result.ptr == last)
  {
    return true;
  }
  std::string into = std::string(first, last);
  try
  {
    value = std::stod(into);
  }
  catch ([[maybe_unused]] std::invalid_argument const& e)
  {
    return false;
  }
  catch ([[maybe_unused]] std::out_of_range const& e)
  {
    return false;
  }
  return true;
}

bool Scanner::scanInt(int& value)
{
  std::string::const_iterator found = find_first_not_of(_charactersToBeSkipped.string(), _string, _scanLocation);

  if (found >= _string.cend())
  {
    _scanLocation = _string.cend();
    return false;
  }

  _scanLocation = found;

  // the token runs up to the next skipped character or to the end of the string
  found = find_first_of(_charactersToBeSkipped.string(), _string, _scanLocation);
  // parsed in place; the library parser is the fallback for the forms 'from_chars' does not take (a leading plus,
  // a leading space), so what was accepted before is accepted still
  const char* first = std::to_address(_scanLocation);
  const char* last = std::to_address(found);
  _scanLocation = found;
  if (const std::from_chars_result result = std::from_chars(first, last, value);
      result.ec == std::errc{} && result.ptr == last)
  {
    return true;
  }
  std::string into = std::string(first, last);
  try
  {
    value = std::stoi(into);
  }
  catch ([[maybe_unused]] std::invalid_argument const& e)
  {
    return false;
  }
  catch ([[maybe_unused]] std::out_of_range const& e)
  {
    return false;
  }
  return true;
}
