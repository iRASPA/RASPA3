#include <gtest/gtest.h>
#include <iostream>

int main(int argc, char** argv)
{
  std::cout << "start" << std::endl;
  testing::InitGoogleTest(&argc, argv);
  std::cout << "here" << std::endl;
  return RUN_ALL_TESTS();
}
