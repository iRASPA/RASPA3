module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>
#endif

export module scanner;

#ifndef USE_LEGACY_HEADERS
import <string>;
import <vector>;
import <iostream>;
import <sstream>;
import <algorithm>;
import <iterator>;
import <iomanip>;
#endif

import characterset;

export class Scanner
{
 public:
  Scanner(const std::string& content, CharacterSet charactersToBeSkipped);
  std::string string() const { return _string; }
  std::string::const_iterator scanLocation() const { return _scanLocation; }
  void setScanLocation(std::string::const_iterator location) { _scanLocation = location; }
  bool scanCharacters(CharacterSet set, std::string& into);
  bool scanLine(std::string& into);
  bool scanUpToCharacters(CharacterSet set, std::string& into);
  bool isAtEnd();
  bool scanDouble(double& value);
  bool scanInt(int& value);
  const std::string displayName() { return _displayName; }

 private:
  CharacterSet _charactersToBeSkipped;
  std::string _displayName;
  std::string _string;
  std::string::const_iterator _scanLocation;

  std::string::const_iterator find_first_not_of(const std::string& chars, const std::string& text,
                                                std::string::const_iterator location);
  std::string::const_iterator find_first_of(const std::string& chars, const std::string& text,
                                            std::string::const_iterator location);
};
