module;

#ifdef USE_LEGACY_HEADERS
#include <string>
#include <vector>
#endif

export module characterset;

#ifndef USE_LEGACY_HEADERS
import <string>;
import <vector>;
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
