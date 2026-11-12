module;

export module characterset;

import std;

export class CharacterSet
{
 public:
  CharacterSet();
  CharacterSet(std::string chars);
  const std::string& string() const { return _string; }
  static CharacterSet newlineCharacter();
  static CharacterSet newlineCharacterSet();
  static CharacterSet whitespaceAndNewlineCharacterSet();
  static CharacterSet whitespaceCharacterSet();

 private:
  std::string _string;
};
