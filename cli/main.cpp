#include <locale.h>

import std;

import opencl;
import commandline;

int main(int argc, char* argv[])
{
  using namespace std::literals;

  setlocale(LC_ALL, "en-US");

  OpenCL::initialize();

  try
  {
    CommandLine::run(argc, argv);
  }
  catch (const std::exception& e)
  {
    std::cerr << e.what();
    std::exit(-1);
  }
}
