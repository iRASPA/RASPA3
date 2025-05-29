module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <exception>
#include <iostream>
#include <fstream>
#include <vector>
#include <span>
#include <deque>
#include <optional>
#include <semaphore>
#include <mutex>
#include <complex>
#include <locale>
#include <ranges>
#include <string_view>
#include <filesystem>
#include <bitset>
#endif

module commandline;

#ifndef USE_LEGACY_HEADERS
import <cstddef>;
import <exception>;
import <iostream>;
import <fstream>;
import <vector>;
import <span>;
import <deque>;
import <optional>;
import <semaphore>;
import <mutex>;
import <complex>;
import <locale>;
import <string_view>;
#endif

import archive;
import int3;
import threadpool;
import input_reader;
import framework;
import forcefield;
import opencl;
import mc_void_fraction;
import mc_surface_area;
import mc_pore_size_distribution;
import mc_opencl_void_fraction;
import mc_opencl_surface_area;
import mc_opencl_pore_size_distribution;
import energy_opencl_void_fraction;
import energy_opencl_surface_area;
import energy_void_fraction;
import energy_surface_area;
import getopt;
#ifdef BUILD_LIBTORCH
import libtorch_test;
#endif


ForceField CommandLine::defaultForceFieldZeolite(double rc, bool shifted, bool tailCorrections, bool useEwald)
{
  return ForceField({{"Si", true, 28.0855, 2.05, 0.0, 14, false},
                     {"O", true, 15.999, -1.025, 0.0, 8, false},
                     {"He", false, 4.002602, 0.0, 0.0, 2, false},
                     {"Ar", false, 39.948, 0.0, 0.0, 18, false},
                     {"CH4", false, 16.04246, 0.0, 0.0, 6, false},
                     {"C_co2", false, 12.0, 0.6512, 0.2, 6, false},
                     {"O_co2", false, 15.9994, -0.3256, 0.1, 8, false}},
                    {{22.0, 2.30},
                     {53.0, 3.30},
                     {10.9, 2.64},
                     {124.070, 3.38},
                     {158.5, 3.72},
                     {29.933, 2.745},
                     {85.671, 3.017}},
                    ForceField::MixingRule::Lorentz_Berthelot, rc, rc, rc, shifted, tailCorrections, useEwald);
}

ForceField CommandLine::defaultForceFieldMOF(double rc, bool shifted, bool tailCorrections, bool useEwald)
{
  return ForceField({{"O",  false, 15.999,    0.0, 0.0,  8, false},
                     {"N",  false, 14.0067,   0.0, 0.0,  7, false},
                     {"C",  false, 12.011,    0.0, 0.0,  6, false},
                     {"F",  false, 18.998403, 0.0, 0.0,  9, false},
                     {"B",  false, 10.811,    0.0, 0.0,  5, false},
                     {"P",  false, 30.973762, 0.0, 0.0, 15, false},
                     {"S",  false, 32.065,    0.0, 0.0, 16, false},
                     {"Cl", false, 35.453,    0.0, 0.0, 17, false},
                     {"Br", false, 79.904,    0.0, 0.0, 35, false},
                     {"H",  false, 1.00784,   0.0, 0.0, 1,  false},
                     {"Al", false, 26.981539, 0.0, 0.0, 13, false},
                     {"Si", false, 28.0855,   0.0, 0.0, 14, false},
                     {"Zn", false, 65.38,     0.0, 0.0, 30, false},
                     {"Be", false, 9.012182,  0.0, 0.0, 4,  false},
                     {"Cr", false, 51.9961,   0.0, 0.0, 24, false},
                     {"Fe", false, 55.845,    0.0, 0.0, 26, false},
                     {"Mn", false, 54.938044, 0.0, 0.0, 25, false},
                     {"Cu", false, 63.546,    0.0, 0.0, 29, false},
                     {"Co", false, 58.933195, 0.0, 0.0, 27, false},
                     {"Ga", false, 72.64,     0.0, 0.0, 32, false},
                     {"Ti", false, 47.867,    0.0, 0.0, 22, false},
                     {"Sc", false, 44.955912, 0.0, 0.0, 21, false},
                     {"V",  false, 50.9415,   0.0, 0.0, 23, false},
                     {"Ni", false, 58.6934,   0.0, 0.0, 28, false},
                     {"Zr", false, 91.224,    0.0, 0.0, 40, false},
                     {"Mg", false, 24.305,    0.0, 0.0, 12, false},
                     {"Ne", false, 20.1797,   0.0, 0.0, 10, false},
                     {"Ag", false, 107.8682,  0.0, 0.0, 47, false},
                     {"In", false, 114.818 ,  0.0, 0.0, 49, false},
                     {"Cd", false, 112.41,    0.0, 0.0, 48, false},
                     {"Sb", false, 121.76,    0.0, 0.0, 51, false},
                     {"Te", false, 127.6,     0.0, 0.0, 52, false},
                     {"He", false, 4.002602,  0.0, 0.0, 2, false},
                     {"Ar", false, 39.948,    0.0, 0.0, 18, false},
                     {"CH4", false, 16.04246, 0.0, 0.0, 6, false},
                     {"C_co2", false, 12.0, 0.6512, 0.2, 6, false},
                     {"O_co2", false, 15.9994, -0.3256, 0.1, 8, false}},
                    //{{ 48.1581,  3.03315},
                    {{ 53.0,       3.30},             // O
                     { 38.9492,    3.26256},          // N
                     { 47.8562,    3.47299},          // C
                     { 36.4834,    3.0932 },          // F
                     { 47.8058,    3.58141},          // B
                     {161.03,      3.69723},          // P
                     {173.107,     3.59032},          // S
                     {142.562,     3.51932},          // Cl
                     {186.191,     3.51905},          // Br
                     {  7.64893,   2.84642},          // H
                     {155.998,     3.91105},          // Al
                     {155.998,     3.80414},          // Si
                     { 62.3988584, 2.461553158},      // Zn
                     { 42.7736,    2.44552},          // Be
                     {  7.54829,   2.69319},          // Cr
                     {  6.54185,   2.5943 },          // Fe
                     {  6.54185,   2.63795},          // Mn
                     {  2.5161,    3.11369},          // Cu
                     {  7.04507,   2.55866},          // Co
                     {208.836,     3.90481},          // Ga
                     {  8.55473,   2.8286 },          // Ti
                     {  9.56117,   2.93551},          // Sc
                     {  8.05151,   2.80099},          // V
                     {  7.54829,   2.52481},          // Ni
                     { 34.7221,    2.78317},          // Zr
                     { 55.8574,    2.69141},          // Mg
                     { 21.1352,    2.88918},          // Ne
                     { 18.1159,    2.80455},          // Ag
                     {301.428,     3.97608},          // In
                     {114.734,     2.53728},          // Cd
                     {225.946,     3.93777},          // Sb
                     {200.281,     3.98232},          // Te
                     { 10.9,       2.64},             // He
                     {124.070,     3.38},             // Ar
                     {158.5,       3.72},             // CH4
                     { 29.933,     2.745},            // C_co2
                     { 85.671,     3.017}},           // O_co2
                    ForceField::MixingRule::Lorentz_Berthelot, rc, rc, rc, shifted, tailCorrections, useEwald);
}

void CommandLine::run(int argc, char* argv[])
{
  // variables modified by command-line switches
  unsigned long num_threads = 1;
  bool use_gridbased_methods{false};
  bool use_monte_carlo_methods{false};
  bool use_energy_methods{false};
  bool use_cpu{false};
  bool use_gpu{false};
  std::bitset<CommandLine::Flag::Last> state;
  std::string input_files;
  size_t number_of_iterations{10000};
  std::optional<ForceField> forceField;
  bool is_zeolite{false};
  bool is_mof{true};

  // definition of command-line switches
  using argparser = argparser::argparser;
  argparser opt(argc, argv);
  opt
    .info("raspa3", argv[0])
    // register command-line switches that execute
    // the built-in help display
    .help({"-h", "--help"}, "Display this help")
    // register command-line switches that modify
    // a variable according to the argument given
    .reg({"-N", "--number-of-iterations"},
         "NUM_ITERATIONS", // will be used in help to illustrate the argument
         argparser::required_argument,
         "Set number of iterations", // will be displayed in help
         [&number_of_iterations](std::string const &arg)
         {
            number_of_iterations = std::stoul(arg);
         })
    .reg({"-p", "--threads"},
         "NUM_THREADS", // will be used in help to illustrate the argument
         argparser::required_argument,
         "Set number of threads", // will be displayed in help
         [&num_threads](std::string const &arg)
         {
            num_threads = std::stoul(arg);
         })
    .reg({"-f", "--force-field"},
         "FILE_NAME", // will be used in help to illustrate the argument
         argparser::required_argument,
         "Sets the force field", // will be displayed in help
         [](std::string const &arg)
         {
           //forceField = ForceField();
           std::cout << "filename: " << arg << std::endl;
         })
    // register command-line switch that doesn't
    // expect an argument
    .reg({"-t", "--opencl-test"},
         argparser::no_argument,
         "Show detected opencl devices",
         [](std::string const &)
         {
            std::cout << OpenCL::printBestOpenCLDevice();
            exit(0);
         })
    .reg({"--zeolite"},
         argparser::no_argument,
         "Use generic zeolite model (TraPPE zeo)",
         [&is_zeolite](std::string const &)
         {
           is_zeolite = true;
         })
    .reg({"--mof"},
         argparser::no_argument,
         "Use generic MOF model (TraPPE zeo)",
         [&is_mof](std::string const &)
         {
           is_mof = true;
         })
    .reg({"--grids"},
         argparser::no_argument,
         "Use grid-based methods",
         [&use_gridbased_methods](std::string const &)
         {
           use_gridbased_methods = true;
         })
    .reg({"--cpu"},
         argparser::no_argument,
         "Compute on the gpu",
         [&use_cpu](std::string const &)
         {
           use_cpu = true;
         })
    .reg({"--gpu"},
         argparser::no_argument,
         "Compute on the gpu",
         [&use_gpu](std::string const &)
         {
           use_gpu = true;
         })
    .reg({"--monte-carlo"},
         argparser::no_argument,
         "Use Monte Carlo-based methods",
         [&use_monte_carlo_methods](std::string const &)
         {
           use_monte_carlo_methods = true;
         })
    .reg({"--energy"},
         argparser::no_argument,
         "Use energy-based methods",
         [&use_energy_methods](std::string const &)
         {
           use_energy_methods = true;
         })
    .reg({"-s", "--surface-area"},
         argparser::no_argument,
         "Compute surface area",
         [&state](std::string const &)
         {
           state.set(Flag::SurfaceArea);
         })
    .reg({"-v", "--void-fraction"},
         argparser::no_argument,
         "Compute void fraction",
         [&state](std::string const &)
         {
           state.set(Flag::VoidFraction);
         })
    .reg({"-p", "--pore-size-distribution"},
         argparser::no_argument,
         "Compute pore size distribution",
         [&state](std::string const &)
         {
           state.set(Flag::PSD);
         })
    // register positional arguments
    .pos("INPUT_CIF_FILE", // will be used in help to illustrate the argument
         "Set CIF file to read", // will be displayed in help
         [&input_files](std::string const &arg)
         { 
           input_files = arg;
         });

  try
  {
    // parse command-line and execute handlers accordingly
    opt();
  }
  // user requested help, help was shown,
  // everything's fine, so exit smoothly
  catch (::argparser::help_requested_exception const &)
  {
    exit(0);
  }
  // handle cases where a required argument was not given
  catch (::argparser::argument_required_exception const &e)
  {
    std::cerr << "\u001b[31;1mERROR: " << e.what() << "\u001b[0m\n";
    exit(-1);
  }
  // handle cases where an unknown switch was given
  catch (::argparser::unknown_option_exception const &e)
  {
    std::cerr << "\u001b[31;1mERROR: " << e.what() << "\u001b[0m\n";
    opt.display_help();
    exit(-2);
  }
  // handle all other exceptions
  catch (std::exception const &e)
  {
    std::cerr << "\u001b[31;1mERROR: " << e.what() << "\u001b[0m\n";
    exit(-3);
  }

  if(!use_cpu && !use_gpu) use_cpu = true;

  // start running
  std::vector<std::string> filenames = std::ranges::to<std::vector<std::string>>(std::string_view(input_files) | std::ranges::views::split(' '));
  for(const std::string &filename : filenames)
  {
    std::cout << "computing structure " << filename << std::endl;

    std::filesystem::path frameworkPathfile = std::filesystem::path(filename);
    if (!std::filesystem::exists(frameworkPathfile))
    {
      throw std::runtime_error(std::format("File '{}' not found\n", filename));
    }

    std::string stem = std::filesystem::path(filename).stem();

    // if no force-field specified, use the default one
    if(is_zeolite)
    {
      forceField = forceField.value_or(defaultForceFieldZeolite(12.0, false, false, false));
    }
    else
    {
      forceField = forceField.value_or(defaultForceFieldMOF(12.0, false, false, false));
    }

    Framework framework = Framework(0, forceField.value(), stem, filename, std::nullopt, Framework::UseChargesFrom::CIF_File);

    if(state.test(CommandLine::Flag::State::SurfaceArea))
    {
      std::cout << "Compute surface area" << std::endl;

      if(use_monte_carlo_methods)
      {
        if(use_cpu)
        {
          MC_SurfaceArea sa;
          sa.run(forceField.value(), framework, 1.0, "Ar", number_of_iterations);
        }

        if(use_gpu)
        {
          EnergyOpenCLSurfaceArea sa;
          sa.run(forceField.value(), framework);
        }
      }

      if(use_energy_methods)
      {
        if(use_cpu)
        {
          EnergySurfaceArea sa;
          sa.run(forceField.value(), framework);
        }

        if(use_gpu)
        {
          EnergyOpenCLSurfaceArea sa;
          sa.run(forceField.value(), framework);
        }
      }
    }

    if(state.test(CommandLine::Flag::State::VoidFraction))
    {
      std::cout << "Compute void fraction" << std::endl;

      if(use_monte_carlo_methods)
      {
        if(use_cpu)
        {
          MC_VoidFraction vf;
          vf.run(forceField.value(), framework, number_of_iterations);
        }

        if(use_gpu)
        {
          EnergyOpenCLVoidFraction vf;
          vf.run(forceField.value(), framework);
        }
      }

      if(use_energy_methods)
      {
        if(use_cpu)
        {
          EnergyVoidFraction vf;
          vf.run(forceField.value(), framework);
        }

        if(use_gpu)
        {
          EnergyOpenCLVoidFraction vf;
          vf.run(forceField.value(), framework);
        }
      }
    }

    if(state.test(CommandLine::Flag::State::PSD))
    {
      std::cout << "Compute PSD" << std::endl;
      if(use_monte_carlo_methods)
      {
        if(use_cpu)
        {
          MC_PoreSizeDistribution psd(1000);
          psd.run(forceField.value(), framework, 1.0, number_of_iterations);
        }

        if(use_gpu)
        {
          MC_OpenCL_PoreSizeDistribution psd(1000);
          psd.run(forceField.value(), framework, 1.0, number_of_iterations);
        }
      }

      if(use_energy_methods)
      {
        std::print("TODO: not implemented yet\n");
      }

    }
  }
}
