module;

#ifdef USE_PRECOMPILED_HEADERS
#include "pch.h"
#include "mdspanwrapper.h"
#endif

#ifdef USE_LEGACY_HEADERS
#include <bitset>
#include <complex>
#include <cstddef>
#include <deque>
#include <exception>
#include <filesystem>
#include <format>
#include <fstream>
#include <iostream>
#include <locale>
#include <mutex>
#include <optional>
#include <print>
#include <ranges>
#include <semaphore>
#include <span>
#include <string_view>
#include <vector>
#include <tuple>
#include "mdspanwrapper.h"
#endif

module commandline;

#ifdef USE_STD_IMPORT
import std;
#endif

import archive;
import int3;
import uint3;
import double3;
import threadpool;
import hdf5;
import input_reader;
import framework;
import vdwparameters;
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
import integration_surface_area;
import integration_opencl_surface_area;
import getopt;
import tessellation;
import interpolation_energy_grid;
import pore_size_distribution_ban_vlugt;
#ifdef BUILD_LIBTORCH
import libtorch_test;
#endif
#if !(defined(__has_include) && __has_include(<mdspan>))
//import mdspan;
#endif

ForceField CommandLine::defaultForceFieldZeolite(double rc, bool shifted, bool tailCorrections, bool useEwald)
{
  return ForceField(
      {{"-", false, 0.0, 0.0, 0.0, 0, false},
       {"Si", true, 28.0855, 2.05, 0.0, 14, false},
       {"O", true, 15.999, -1.025, 0.0, 8, false},
       {"He", false, 4.002602, 0.0, 0.0, 2, false},
       {"Ar", false, 39.948, 0.0, 0.0, 18, false},
       {"CH4", false, 16.04246, 0.0, 0.0, 6, false},
       {"C_co2", false, 12.0, 0.6512, 0.2, 6, false},
       {"O_co2", false, 15.9994, -0.3256, 0.1, 8, false},
       {"N", false, 14.00674, 0.0, 0.0, 6, false}},
      {{1.0, 1.0}, 
      {22.0, 2.30}, 
      {53.0, 3.30}, 
      {10.9, 2.64},
      {124.070, 3.38},
      {158.5, 3.72},
      {29.933, 2.745},
      {85.671, 3.017},
      {91.5, 3.681}},
      ForceField::MixingRule::Lorentz_Berthelot, rc, rc, rc, shifted, tailCorrections, useEwald);
}

ForceField CommandLine::defaultForceFieldMOF(double rc, bool shifted, bool tailCorrections, bool useEwald)
{
  return ForceField({{"-", false, 0.0, 0.0, 0.0, 0, false},
                     {"O",  true, 15.999, 0.0, 0.0, 8, false},         {"N",  true, 14.0067, 0.0, 0.0, 7, false},
                     {"C",  true, 12.011, 0.0, 0.0, 6, false},         {"F",  true, 18.998403, 0.0, 0.0, 9, false},
                     {"B",  true, 10.811, 0.0, 0.0, 5, false},         {"P",  true, 30.973762, 0.0, 0.0, 15, false},
                     {"S",  true, 32.065, 0.0, 0.0, 16, false},        {"Cl", true, 35.453, 0.0, 0.0, 17, false},
                     {"Br", true, 79.904, 0.0, 0.0, 35, false},        {"H",  true, 1.00784, 0.0, 0.0, 1, false},
                     {"Al", true, 26.981539, 0.0, 0.0, 13, false},     {"Si", true, 28.0855, 0.0, 0.0, 14, false},
                     {"Zn", true, 65.38, 0.0, 0.0, 30, false},         {"Be", true, 9.012182, 0.0, 0.0, 4, false},
                     {"Cr", true, 51.9961, 0.0, 0.0, 24, false},       {"Fe", true, 55.845, 0.0, 0.0, 26, false},
                     {"Mn", true, 54.938044, 0.0, 0.0, 25, false},     {"Cu", true, 63.546, 0.0, 0.0, 29, false},
                     {"Co", true, 58.933195, 0.0, 0.0, 27, false},     {"Ga", true, 72.64, 0.0, 0.0, 32, false},
                     {"Ti", true, 47.867, 0.0, 0.0, 22, false},        {"Sc", true, 44.955912, 0.0, 0.0, 21, false},
                     {"V",  true, 50.9415, 0.0, 0.0, 23, false},       {"Ni", true, 58.6934, 0.0, 0.0, 28, false},
                     {"Zr", true, 91.224, 0.0, 0.0, 40, false},        {"Mg", true, 24.305, 0.0, 0.0, 12, false},
                     {"Ne", true, 20.1797, 0.0, 0.0, 10, false},       {"Ag", true, 107.8682, 0.0, 0.0, 47, false},
                     {"In", true, 114.818, 0.0, 0.0, 49, false},       {"Cd", true, 112.41, 0.0, 0.0, 48, false},
                     {"Sb", true, 121.76, 0.0, 0.0, 51, false},        {"Te", true, 127.6, 0.0, 0.0, 52, false},
                     {"He", false, 4.002602, 0.0, 0.0, 2, false},      {"Ar", true, 39.948, 0.0, 0.0, 18, false},
                     {"CH4", false, 16.04246, 0.0, 0.0, 6, false},     {"C_co2", false, 12.0, 0.6512, 0.2, 6, false},
                     {"O_co2", false, 15.9994, -0.3256, 0.1, 8, false}},
                    //{{ 48.1581,  3.03315},
                    {{1.0, 1.0},                 // custom
                     {48.1581, 3.03315},         // O
                     {38.9492, 3.26256},         // N
                     {47.8562, 3.47299},         // C
                     {36.4834, 3.0932},          // F
                     {47.8058, 3.58141},         // B
                     {161.03, 3.69723},          // P
                     {173.107, 3.59032},         // S
                     {142.562, 3.51932},         // Cl
                     {186.191, 3.51905},         // Br
                     {7.64893, 2.84642},         // H
                     {155.998, 3.91105},         // Al
                     {155.998, 3.80414},         // Si
                     {27.677074, 4.04},          // Zn
                     {42.7736, 2.44552},         // Be
                     {7.54829, 2.69319},         // Cr
                     {6.54185, 2.5943},          // Fe
                     {6.54185, 2.63795},         // Mn
                     {2.5161, 3.11369},          // Cu
                     {7.04507, 2.55866},         // Co
                     {208.836, 3.90481},         // Ga
                     {8.55473, 2.8286},          // Ti
                     {9.56117, 2.93551},         // Sc
                     {8.05151, 2.80099},         // V
                     {7.54829, 2.52481},         // Ni
                     {34.7221, 2.78317},         // Zr
                     {55.8574, 2.69141},         // Mg
                     {21.1352, 2.88918},         // Ne
                     {18.1159, 2.80455},         // Ag
                     {301.428, 3.97608},         // In
                     {114.734, 2.53728},         // Cd
                     {225.946, 3.93777},         // Sb
                     {200.281, 3.98232},         // Te
                     {10.9, 2.64},               // He
                     {124.070, 3.38},            // Ar
                     {158.5, 3.72},              // CH4
                     {29.933, 2.745},            // C_co2
                     {85.671, 3.017}},           // O_co2
                    ForceField::MixingRule::Lorentz_Berthelot, rc, rc, rc, shifted, tailCorrections, useEwald);
}

void CommandLine::run(int argc, char *argv[])
{
  // variables modified by command-line switches
  unsigned long num_threads = 1;
  bool use_gridbased_methods{false};
  bool use_geometric_methods{false};
  bool use_monte_carlo_methods{false};
  bool use_energy_methods{false};
  bool use_integration_methods{false};
  bool use_cpu{false};
  bool use_gpu{false};
  std::bitset<CommandLine::State::Last> state;
  std::string input_files;
  std::optional<std::size_t> number_of_iterations{};
  std::optional<std::size_t> number_of_inner_steps{};
  std::optional<double> minimum_range{};
  std::optional<double> maximum_range{};
  std::optional<std::string> probe_atom_name{};
  std::optional<double> probe_size{};
  std::optional<double> probe_strength{};
  double well_depth_factor{ 1.0 };
  double iso_value{ 0.0 };
  std::optional<ForceField> forceField;
  bool is_zeolite{false};
  bool is_mof{true};
  uint3 gridSize{128, 128, 128};
  std::optional<std::size_t> number_of_slices{ };
  std::vector<std::size_t> pseudoAtomsGrid;
  ForceField::InterpolationScheme order{ForceField::InterpolationScheme::Tricubic};
  ForceField::InterpolationGridType gridType{ForceField::InterpolationGridType::LennardJones};

  // definition of command-line switches
  using argparser = argparser::argparser;
  argparser opt(argc, argv);
  opt.info("raspa3", argv[0])
      // register command-line switches that execute
      // the built-in help display
      .help({"-h", "--help"}, "Display this help")
      // register command-line switches that modify
      // a variable according to the argument given
      .reg({"-N", "--number-of-iterations"},
           "NUM_ITERATIONS",  // will be used in help to illustrate the argument
           argparser::required_argument,
           "Set number of iterations",  // will be displayed in help
           [&number_of_iterations](std::string const &arg) { std::print("arg: {}\n",arg); number_of_iterations = std::stoul(arg); std::print("arg: {}\n",arg); })
      .reg({"-M", "--number-of-inner-steps"},
           "NUM_INNER_ITERATIONS",  // will be used in help to illustrate the argument
           argparser::required_argument,
           "Set number of inner steps",  // will be displayed in help
           [&number_of_inner_steps](std::string const &arg) { std::print("arg: {}\n",arg); number_of_inner_steps = std::stoul(arg); std::print("arg: {}\n",arg); })
      .reg({"--threads"},
           "NUM_THREADS",  // will be used in help to illustrate the argument
           argparser::required_argument,
           "Set number of threads",  // will be displayed in help
           [&num_threads](std::string const &arg) { std::print("arg: {}\n",arg); num_threads = std::stoul(arg); })
      .reg({"-f", "--force-field"},
           "FILE_NAME",  // will be used in help to illustrate the argument
           argparser::required_argument,
           "Sets the force field",  // will be displayed in help
           [&forceField](std::string const &arg)
           {
             forceField = ForceField(arg);
             std::cout << "filename: " << arg << std::endl;
           })
      // register command-line switch that doesn't
      // expect an argument
      .reg({"-t", "--opencl-test"}, argparser::no_argument, "Show detected opencl devices",
           [](std::string const &)
           {
             std::cout << OpenCL::printBestOpenCLDevice();
             std::exit(0);
           })
      .reg({"-s", "--surface-area"}, argparser::no_argument, "Compute surface area",
           [&state](std::string const &) { state.set(State::SurfaceArea); })
      .reg({"-v", "--void-fraction"}, argparser::no_argument, "Compute void fraction",
           [&state](std::string const &) { state.set(State::VoidFraction); })
      .reg({"-p", "--pore-size-distribution"}, argparser::no_argument, "Compute pore size distribution",
           [&state](std::string const &) { state.set(State::PSD); })
      .reg({"--pore-size-distribution-ban-vlugt"}, argparser::no_argument,
           "Use pore size distribution method from Ban, Vlugt paper",
           [&state](std::string const &){state.set(State::PSD_BV);})
      .reg({"--tessellation"}, argparser::no_argument, "Use tessellation method",
           [&state](std::string const &) { state.set(State::TessellationComputation); })
      .reg({"--zeolite"}, argparser::no_argument, "Use generic zeolite model (TraPPE zeo)",
           [&is_zeolite](std::string const &) { is_zeolite = true; })
      .reg({"--mof"}, argparser::no_argument, "Use generic MOF model (TraPPE zeo)",
           [&is_mof](std::string const &) { is_mof = true; })
      .reg({"--grids"}, argparser::no_argument, "Use grid-based methods",
           [&use_gridbased_methods](std::string const &) { use_gridbased_methods = true; })
      .reg({"--geometric"}, argparser::no_argument, "Use geometric methods",
           [&use_geometric_methods](std::string const &) { use_geometric_methods = true; })
      .reg({"--energy"}, argparser::no_argument, "Use energy-based methods",
           [&use_energy_methods](std::string const &) { use_energy_methods = true; })
      .reg({"--monte-carlo"}, argparser::no_argument, "Use Monte Carlo algorithm",
           [&use_monte_carlo_methods](std::string const &) { use_monte_carlo_methods = true; })
      .reg({"--integration"}, argparser::no_argument, "Use integration algorithm",
           [&use_integration_methods](std::string const &) { use_integration_methods = true; })
      .reg({"--128"}, argparser::no_argument, "Use low-accuracy 128x128x128 grid",
           [&gridSize](std::string const &) { gridSize = uint3(128, 128, 128); })
      .reg({"--256"}, argparser::no_argument, "Use medium-accuracy 256x256x256 grid",
           [&gridSize](std::string const &) { gridSize = uint3(256, 256, 256); })
      .reg({"--512"}, argparser::no_argument, "Use high-accuracy 512x512x512 grid",
           [&gridSize](std::string const &) { gridSize = uint3(512, 512, 512); })
      .reg({"--grid-size"}, argparser::required_argument, "Set size of grid: 'X Y Z'",
           [&gridSize](std::string const &arg)
           {
             std::istringstream iss(arg);
             std::size_t x, y, z;
             if (!(iss >> x >> y >> z))
             {
               throw std::runtime_error(
                   "Invalid --grid-size: expected three integers separated by spaces wrapped in quotation marks");
             }
             gridSize = uint3(x, y, z);
           })
      .reg({"--minimum-range"},
           argparser::required_argument,
           "Minimum range", 
           [&minimum_range](std::string const &arg) { minimum_range = std::stod(arg); })
      .reg({"--maximum-range"},
           argparser::required_argument,
           "Maximum range", 
           [&maximum_range](std::string const &arg) { maximum_range = std::stod(arg); })
      .reg({"--probe-atom-name"},
           argparser::required_argument,
           "The name of the probe atom", 
           [&probe_atom_name](std::string const &arg) { probe_atom_name = arg; })
      .reg({"--probe-size-parameter"},
           argparser::required_argument,
           "The size of the probe atom", 
           [&probe_size](std::string const &arg) { probe_size = std::stod(arg); })
      .reg({"--probe-strength-parameter"},
           argparser::required_argument,
           "The strength of the probe atom", 
           [&probe_strength](std::string const &arg) { probe_strength = std::stod(arg); })
      .reg({"--use-well-depth-as-size"},
           argparser::no_argument,
           "Uses the well-depth location (1.1225 times sigma) as the size of the atom", 
           [&well_depth_factor](std::string const &) { well_depth_factor = std::pow(2.0, 1.0 / 6.0); })
      .reg({"--iso-value"},
           argparser::required_argument,
           "Sets the isovalue for the iso-surface energy surfaces", 
           [&iso_value](std::string const &arg) { iso_value = std::stod(arg); })
      .reg({"--number-of-slices"}, argparser::required_argument, "Set number of slices",
           [&number_of_slices](std::string const &arg)
           {
             std::istringstream iss(arg);
             std::size_t x;
             if (!(iss >> x))
             {
               throw std::runtime_error(
                   "Invalid --grid-size: expected three integers separated by spaces wrapped in quotation marks");
             }
             number_of_slices = x;
           })
      .reg({"-e", "--energy-grid-computation"}, argparser::required_argument,
           "Compute Energy grid for given list of pseudo atoms",
           [&state, &pseudoAtomsGrid](std::string const &arg)
           {
             state.set(State::EnergyGrid);
             std::istringstream iss(arg);
             std::size_t idx;
             while (iss >> idx)
             {
               pseudoAtomsGrid.push_back(idx);
             }
             if (pseudoAtomsGrid.empty())
             {
               throw std::runtime_error("Invalid Energy Grid Computation, list pseudo atom indices with '0 1 2 ...'");
             }
           })
      .reg({"--tricubic"}, argparser::no_argument, "Set interpolation scheme to tricubic",
           [&order](std::string const &) { order = ForceField::InterpolationScheme::Tricubic; })
      .reg({"--triquintic"}, argparser::no_argument, "Set interpolation scheme to triquintic",
           [&order](std::string const &) { order = ForceField::InterpolationScheme::Triquintic; })
      .reg({"--Lennard-Jones"}, argparser::no_argument, "Set interpolation energy grid to Lennard-Jones",
           [&gridType](std::string const &) { gridType = ForceField::InterpolationGridType::LennardJones; })
      .reg({"--Ewald"}, argparser::no_argument, "Set interpolation energy grid to Ewald",
           [&gridType](std::string const &) { gridType = ForceField::InterpolationGridType::EwaldReal; })
      .reg({"--cpu"}, argparser::no_argument, "Compute on the gpu", [&use_cpu](std::string const &) { use_cpu = true; })
      .reg({"--gpu"}, argparser::no_argument, "Compute on the gpu", [&use_gpu](std::string const &) { use_gpu = true; })
      // register positional arguments
      .pos("INPUT_CIF_FILE",        // will be used in help to illustrate the argument
           "Set CIF file to read",  // will be displayed in help
           [&input_files](std::string const &arg) { input_files = arg; });

  try
  {
    // parse command-line and execute handlers accordingly
    opt();
  }
  // user requested help, help was shown,
  // everything's fine, so exit smoothly
  catch (::argparser::help_requested_exception const &)
  {
    std::exit(0);
  }
  // handle cases where a required argument was not given
  catch (::argparser::argument_required_exception const &e)
  {
    std::cerr << "\u001b[31;1mERROR: " << e.what() << "\u001b[0m\n";
    std::exit(-1);
  }
  // handle cases where an unknown switch was given
  catch (::argparser::unknown_option_exception const &e)
  {
    std::cerr << "\u001b[31;1mERROR: " << e.what() << "\u001b[0m\n";
    opt.display_help();
    std::exit(-2);
  }
  // handle all other exceptions
  catch (std::exception const &e)
  {
    std::cerr << "\u001b[31;1mERROR: " << e.what() << "\u001b[0m\n";
    std::exit(-3);
  }

  if (!use_cpu && !use_gpu) use_cpu = true;

  // start running
  std::vector<std::string> filenames =
      std::ranges::to<std::vector<std::string>>(std::string_view(input_files) | std::ranges::views::split(' '));
  for (const std::string &filename : filenames)
  {
    std::cout << "computing structure " << filename << std::endl;

    std::filesystem::path frameworkPathfile = std::filesystem::path(filename);
    if (!std::filesystem::exists(frameworkPathfile))
    {
      throw std::runtime_error(std::format("File '{}' not found\n", filename));
    }

    std::string stem = std::filesystem::path(filename).stem().string();

    // if no force-field specified, use the default one
    if (!forceField.has_value())
    {
      if (is_zeolite)
      {
        forceField = forceField.value_or(defaultForceFieldZeolite(12.0, false, false, false));
      }
      else
      {
        forceField = forceField.value_or(defaultForceFieldMOF(12.0, false, false, false));
      }
    }

    // Handle custom probe size
    if(probe_size.has_value())
    {
      forceField->data.front() = VDWParameters(probe_strength.value_or(1.0), probe_size.value());
      forceField->applyMixingRule();
      forceField->preComputePotentialShift();
      forceField->preComputeTailCorrection();
      probe_atom_name = "-";
    }

    Framework framework =
        Framework(0, forceField.value(), stem, filename, std::nullopt, Framework::UseChargesFrom::CIF_File);

    if (state.test(CommandLine::State::TessellationComputation))
    {
      std::cout << "Compute tesselation" << std::endl;

      if (use_gridbased_methods)
      {
        if (use_cpu)
        {
        }

        if (use_gpu)
        {
          Tessellation s(gridSize);
          s.run(forceField.value(), framework);
        }
      }
    }

    if(state.test(CommandLine::State::PSD_BV))
    {
      std::cout << "Compute pore size distribution using the method from the Ban, Vlugt paper" << std::endl;

      if(use_gridbased_methods)
      {
        if(use_cpu)
        {
        }

        if(use_gpu)
        {
          BanVlugtPoreSizeDistribution s(gridSize);
          s.run(forceField.value(), framework);
        }
      }
    }

    if (state.test(CommandLine::State::SurfaceArea))
    {
      std::cout << "Compute surface area" << std::endl;

      if (use_geometric_methods)
      {
        if (use_monte_carlo_methods)
        {
          if (use_cpu)
          {
            MC_SurfaceArea sa;

            sa.run(forceField.value(), framework, well_depth_factor, probe_atom_name.value_or("Ar"), number_of_iterations, number_of_inner_steps);
          }

          if (use_gpu)
          {
            MC_OpenCL_SurfaceArea sa;
            sa.run(forceField.value(), framework, well_depth_factor, probe_atom_name.value_or("Ar"), number_of_iterations, number_of_inner_steps);
          }
        }

        if (use_integration_methods)
        {
          if (use_cpu)
          {
            Integration_SurfaceArea sa;
            sa.run(forceField.value(), framework, well_depth_factor, probe_atom_name.value_or("Ar"), number_of_slices);
          }

          if (use_gpu)
          {
            Integration_OpenCL_SurfaceArea sa;
            sa.run(forceField.value(), framework, well_depth_factor, probe_atom_name.value_or("Ar"), number_of_slices);
          }
        }
      }

      if (use_energy_methods)
      {
        if (use_cpu)
        {
          EnergySurfaceArea sa;
          sa.run(forceField.value(), framework, iso_value, probe_atom_name.value_or("Ar"));
        }

        if (use_gpu)
        {
          EnergyOpenCLSurfaceArea sa;
          sa.run(forceField.value(), framework, iso_value, probe_atom_name.value_or("Ar"), gridSize);
        }
      }
    }

    if (state.test(CommandLine::State::VoidFraction))
    {
      std::cout << "Compute void fraction" << std::endl;

      if (use_monte_carlo_methods)
      {
        if (use_cpu)
        {
          MC_VoidFraction vf;
          vf.run(forceField.value(), framework, probe_atom_name.value_or("He"), number_of_iterations);
        }

        if (use_gpu)
        {
          EnergyOpenCLVoidFraction vf;
          vf.run(forceField.value(), framework);
        }
      }

      if (use_energy_methods)
      {
        if (use_cpu)
        {
          EnergyVoidFraction vf;
          vf.run(forceField.value(), framework);
        }

        if (use_gpu)
        {
          EnergyOpenCLVoidFraction vf;
          vf.run(forceField.value(), framework);
        }
      }
    }

    if (state.test(CommandLine::State::PSD))
    {
      if (use_monte_carlo_methods)
      {
        if (use_cpu)
        {
          MC_PoreSizeDistribution psd(1000);
          psd.run(forceField.value(), framework, 1.0, number_of_iterations, number_of_inner_steps, maximum_range);
        }

        if (use_gpu)
        {
          MC_OpenCL_PoreSizeDistribution psd(1000);
          psd.run(forceField.value(), framework, 1.0, number_of_iterations, number_of_inner_steps, maximum_range);
        }
      }

      if (use_energy_methods)
      {
        std::print("TODO: not implemented yet\n");
      }
    }

    if (state.test(CommandLine::State::EnergyGrid))
    {
      HDF5Writer h5File(
          std::format("{}_{}_{}.h5", framework.name,
                      (order == ForceField::InterpolationScheme::Tricubic) ? "Tricubic" : "Triquintic",
                      (gridType == ForceField::InterpolationGridType::LennardJones) ? "LennardJones" : "Ewald"));
      for (std::size_t pseudoAtomIdx : pseudoAtomsGrid)
      {
        InterpolationEnergyGrid grid(framework.simulationBox, double3{},
                                     gridSize, order);
        std::ostream nullStream{nullptr};
        grid.makeFrameworkInterpolationGrid(nullStream, gridType, forceField.value(), framework,
                                   (gridType == ForceField::InterpolationGridType::LennardJones)
                                       ? forceField->cutOffFrameworkVDW
                                       : forceField->cutOffCoulomb,
                                   pseudoAtomIdx);

        // unpreferable copy implementation to reorder fortran layout of grid data to c layout
        std::vector<double> reordered(grid.data.size());
        std::mdspan<double, std::dextents<std::size_t, 4>, std::layout_left> spanFortran(
            grid.data.data(), std::to_underlying(order), gridSize.x + 1, gridSize.y + 1, gridSize.z + 1);
        std::mdspan<double, std::dextents<std::size_t, 4>> spanC(reordered.data(), std::to_underlying(order),
                                                                 gridSize.x + 1, gridSize.y + 1, gridSize.z + 1);

        for (std::size_t v = 0; v < std::to_underlying(order); v++)
        {
          for (std::size_t ix = 0; ix < static_cast<std::size_t>(gridSize.x + 1); ix++)
          {
            for (std::size_t iy = 0; iy < static_cast<std::size_t>(gridSize.y + 1); iy++)
            {
              for (std::size_t iz = 0; iz < static_cast<std::size_t>(gridSize.z + 1); iz++)
              {
                spanC[v, ix, iy, iz] = spanFortran[v, ix, iy, iz];
              }
            }
          }
        }

        h5File.createDataset<double>(
            "/", forceField->pseudoAtoms[pseudoAtomIdx].name,
            {std::to_underlying(order), static_cast<std::size_t>(gridSize.x + 1),
             static_cast<std::size_t>(gridSize.y + 1), static_cast<std::size_t>(gridSize.z + 1)},
            {});
        h5File.writeVector<double>("/", forceField->pseudoAtoms[pseudoAtomIdx].name, reordered);
      }
    }
  }
}
