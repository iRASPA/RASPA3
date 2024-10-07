#include <gtest/gtest.h>

#include <bit>
#include <complex>
#include <fstream>
#include <iostream>
#include <print>

/*
import archive;
struct teststruct
{
  size_t a{0};
  int b{0};
  int c{0};
  std::vector<int> d{};
  size_t e{0};

  bool operator==(const teststruct& rhs) const {
    return
       a == rhs.a
       && b == rhs.b
       && c == rhs.c
       && d == rhs.d
       && e == rhs.e
    ;
}

  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& stream, const teststruct& mc);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& stream, teststruct& mc);
};

Archive<std::ofstream>& operator<<(Archive<std::ofstream>& stream, const teststruct& mc)
{
  stream << mc.a;
  stream << mc.b;
  stream << mc.c;
  stream << mc.d;
  stream << mc.e;
  return stream;
}

Archive<std::ifstream>& operator>>(Archive<std::ifstream>& stream, teststruct& mc)
{
  stream >> mc.a;
  stream >> mc.b;
  stream >> mc.c;
  stream >> mc.d;
  stream >> mc.e;
  return stream;
}


TEST(archive, Test_archiving)
{
  std::ofstream ofile("test.bin", std::ios::binary);
  Archive<std::ofstream> archive(ofile);

  teststruct a;
  a.a = 100;
  a.b = 3;
  a.c = 4;
  a.d = std::vector<int>{5,6,7};
  a.e = 333;
  archive << a;

  ofile.close();

  std::ifstream ifile("test.bin", std::ios::binary);
  if(!ifile.is_open())
  {
      std::cout << "ERROR file doesn't exist.." << std::endl;
      return;
  }

  Archive<std::ifstream> archive2(ifile);
  teststruct b{};
  archive2 >> b;
  ifile.close();

  EXPECT_EQ(a, b);
}
*/
