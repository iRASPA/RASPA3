#include <gtest/gtest.h>

import std;

import getopt;

// An option that requires a value takes the next argument unless it is an option itself: '--grid-size --gpu' took
// '--gpu' for the grid size and dropped the option. The missing value is reported instead; a value that looks like
// an option but is none (a negative number) is still passed on.
TEST(getopt, a_required_argument_does_not_take_the_next_option)
{
  const auto parse = [](std::vector<std::string> storage, std::vector<std::string>& values, bool& gpu)
  {
    std::vector<char*> argv{};
    for (std::string& argument : storage) argv.push_back(argument.data());
    argparser::argparser parser(static_cast<int>(argv.size()), argv.data());
    parser.reg({"--grid-size"}, "size", argparser::argparser::required_argument, "",
               [&](const std::string& value) { values.push_back(value); });
    parser.reg({"--gpu"}, argparser::argparser::no_argument, "", [&](const std::string&) { gpu = true; });
    parser();
  };

  std::vector<std::string> values;
  bool gpu = false;
  EXPECT_THROW(parse({"raspa3", "--grid-size", "--gpu"}, values, gpu), argparser::argument_required_exception);
  EXPECT_TRUE(values.empty());

  values.clear();
  gpu = false;
  parse({"raspa3", "--grid-size", "-1", "--gpu"}, values, gpu);
  ASSERT_EQ(values.size(), 1uz);
  EXPECT_EQ(values[0], "-1");
  EXPECT_TRUE(gpu);
}
