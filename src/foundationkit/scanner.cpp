module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>
#endif

module scanner;

#ifndef USE_LEGACY_HEADERS
import <string>;
import <vector>;
import <iostream>;
import <sstream>;
import <algorithm>;
import <iterator>;
import <iomanip>;
import <string_view>;
import <type_traits>;
#endif

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

  if (found >= _string.cend())
  {
    _scanLocation = _string.cend();
    into = std::string("");
    return false;
  }

  if (found < _string.cend())
  {
    into = std::string(_scanLocation, found);
    // into = std::string(_scanLocation, int(found - _scanLocation));
    _scanLocation = found;
    return true;
  }
  _scanLocation = _string.cend();
  into = std::string("");
  return false;
}

bool Scanner::scanLine(std::string& into)
{
  std::string::const_iterator newlineLocation =
      find_first_of(CharacterSet::newlineCharacter().string(), _string, _scanLocation);

  // not found
  if (newlineLocation >= _string.cend())
  {
    _scanLocation = _string.cend();
    into = std::string("");
    return false;
  }

  // found
  if (newlineLocation < _string.cend())
  {
    // into = std::string(_scanLocation, newlineLocation - _scanLocation);
    into = std::string(_scanLocation, newlineLocation);
    _scanLocation = newlineLocation + 1;
    return true;
  }

  // found newline as the last character
  _scanLocation = _string.cend();
  into = std::string("");
  return false;
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

  found = find_first_of(set.string(), _string, _scanLocation);

  if (found < _string.cend())
  {
    // into = std::string(_scanLocation, found - _scanLocation);
    into = std::string(_scanLocation, found);
    _scanLocation = found;
    return true;
  }
  _scanLocation = _string.cend();
  into = std::string("");
  return false;
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

  found = find_first_of(_charactersToBeSkipped.string(), _string, _scanLocation);

  if (found < _string.cend())
  {
    // std::string into = std::string(_scanLocation, found - _scanLocation);
    std::string into = std::string(_scanLocation, found);
    _scanLocation = found;
    bool success = false;
    // value = into.toDouble(&success);
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
    return success;
  }
  _scanLocation = _string.cend();
  return false;
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

  found = find_first_of(_charactersToBeSkipped.string(), _string, _scanLocation);

  if (found < _string.cend())
  {
    // std::string into = std::string(_scanLocation, found - _scanLocation);
    std::string into = std::string(_scanLocation, found);
    _scanLocation = found;
    bool success = false;
    // value = into.toInt(&success);
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
    return success;
  }
  _scanLocation = _string.cend();
  return false;
}
