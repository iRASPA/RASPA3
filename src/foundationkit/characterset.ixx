module;

#ifdef USE_PRECOMPILED_HEADERS
#include "pch.h"
#endif

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <string>
#include <vector>
#endif

export module characterset;

#ifdef USE_STD_IMPORT
import std;
#endif

export class CharacterSet
{
 public:
  CharacterSet();
  CharacterSet(std::string chars);
  std::string string() { return _string; }
  static CharacterSet newlineCharacter();
  static CharacterSet newlineCharacterSet();
  static CharacterSet whitespaceAndNewlineCharacterSet();
  static CharacterSet whitespaceCharacterSet();

 private:
  std::string _string;
};
