export module skonethirdseitzmatrix;

import <string>;

export class SKOneThirdSeitzMatrix
{
public:
	SKOneThirdSeitzMatrix(std::string text, uint8_t encoding, int8_t r1, int8_t r2, int8_t r3, int8_t t);
	std::string text() { return _text; }
	uint8_t encoding() { return _encoding; }
	int8_t r1() { return _r1; }
	int8_t r2() { return _r2; }
	int8_t r3() { return _r3; }
	int8_t t() { return _t; }
private:
	std::string _text;
	uint8_t _encoding;
	int8_t _r1;
	int8_t _r2;
	int8_t _r3;
	int8_t _t;
};