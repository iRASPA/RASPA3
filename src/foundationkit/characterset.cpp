module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <string>
#endif

module characterset;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

CharacterSet::CharacterSet() {}

CharacterSet::CharacterSet(std::string chars) { _string = chars; }

// \x0A   \n    Linefeed (LF)
// \x0B   \e    Escape (ESC)
// \x0C   \f    Formfeed (FF)
// \x0D   \r    Carriage return (CR)
// \x85

CharacterSet CharacterSet::newlineCharacter() { return CharacterSet(std::string("\x0A", 1)); }

CharacterSet CharacterSet::newlineCharacterSet() { return CharacterSet(std::string("\x0A\x0B\x0C\x0D\x85", 5)); }

CharacterSet CharacterSet::whitespaceAndNewlineCharacterSet()
{
  return CharacterSet(std::string("\x0A\x0B\x0C\x0D\x85\x0D\x85\x20\x09\xA0", 8));
}

CharacterSet CharacterSet::whitespaceCharacterSet() { return CharacterSet(std::string("\x20\x09\xA0", 3)); }
