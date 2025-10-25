#ifdef USE_LEGACY_HEADERS
#include <gtest/gtest.h>

#include <iostream>
#endif

int main(int argc, char** argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
