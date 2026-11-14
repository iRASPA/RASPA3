module;

module characterset;

import std;

CharacterSet::CharacterSet() {}

CharacterSet::CharacterSet(std::string chars) { _string = chars; }

// \x0A   \n    Linefeed (LF)
// \x0B   \v    Vertical tab (VT)
// \x0C   \f    Formfeed (FF)
// \x0D   \r    Carriage return (CR)
//
// The sets are matched byte by byte on UTF-8 text: the Latin-1 bytes NEL (\x85) and NBSP (\xA0) are continuation
// bytes there (the 0xC3 0x85 of 'Å'), and are not part of the sets.

CharacterSet CharacterSet::newlineCharacter() { return CharacterSet(std::string("\x0A", 1)); }

CharacterSet CharacterSet::newlineCharacterSet() { return CharacterSet(std::string("\x0A\x0B\x0C\x0D", 4)); }

CharacterSet CharacterSet::whitespaceAndNewlineCharacterSet()
{
  // the newline set and the whitespace set together
  return CharacterSet(std::string("\x0A\x0B\x0C\x0D\x20\x09", 6));
}

CharacterSet CharacterSet::whitespaceCharacterSet() { return CharacterSet(std::string("\x20\x09", 2)); }
