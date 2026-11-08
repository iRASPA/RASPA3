module;

module scanner;

import std;

import characterset;

Scanner::Scanner(const std::string& content, CharacterSet charactersToBeSkipped)
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

bool Scanner::scanCharacters(CharacterSet set, std::string& into)
{
  std::string::const_iterator found = find_first_not_of(set.string(), _string, _scanLocation);

  // the characters of the set up to the next other character, or up to the end of the string
  if (found >= _string.cend() && _scanLocation >= _string.cend())
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
  std::string::const_iterator newlineLocation =
      find_first_of(CharacterSet::newlineCharacter().string(), _string, _scanLocation);

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

bool Scanner::scanUpToCharacters(CharacterSet set, std::string& into)
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
  std::string into = std::string(_scanLocation, found);
  _scanLocation = found;
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
  std::string into = std::string(_scanLocation, found);
  _scanLocation = found;
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
