export module characterset;

import <string>;
import <vector>;

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