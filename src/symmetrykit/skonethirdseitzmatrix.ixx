module;

export module skonethirdseitzmatrix;

import std;

export class SKOneThirdSeitzMatrix
{
 public:
  SKOneThirdSeitzMatrix(std::string text, std::uint8_t encoding, std::int8_t r1, std::int8_t r2, std::int8_t r3,
                        std::int8_t t);
  std::string text() { return _text; }
  std::uint8_t encoding() { return _encoding; }
  std::int8_t r1() { return _r1; }
  std::int8_t r2() { return _r2; }
  std::int8_t r3() { return _r3; }
  std::int8_t t() { return _t; }

 private:
  std::string _text;
  std::uint8_t _encoding;
  std::int8_t _r1;
  std::int8_t _r2;
  std::int8_t _r3;
  std::int8_t _t;
};
