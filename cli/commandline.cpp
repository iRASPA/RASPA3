module;

module commandline;

import std;

import archive;
import int3;
import uint3;
import double3;
import threadpool;
import stringutils;
import hdf5;
import input_reader;
import framework;
import vdwparameters;
import forcefield;
import cif_reader;
import atom;
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
import mdspan;
#endif


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
  uint3 gridSize{128, 128, 128};
  std::optional<std::size_t> number_of_slices{ };
  std::vector<std::size_t> pseudoAtomsGrid;
  ForceField::InterpolationScheme order{ForceField::InterpolationScheme::Tricubic};
  ForceField::InterpolationGridType gridType{ForceField::InterpolationGridType::LennardJones};
  std::optional<ForceField> forceField{};
  Framework framework{};


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
           [&forceField](std::string const &) { forceField = ForceField::makeZeoliteForceField(12.0, true, false, false);; })
      .reg({"--mof"}, argparser::no_argument, "Use generic MOF model (TraPPE zeo)",
           [&forceField](std::string const &) { forceField = ForceField::makeMetalOrganicFrameworkForceField(12.0, true, false, false); })
      .reg({"--zeo++"}, argparser::no_argument, "Use zeo++ radii",
           [&forceField](std::string const &) { forceField = ForceField::makeZeoPlusPlusForceField(12.0, true, false, false); })
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
           [&probe_atom_name](std::string const &arg) { probe_atom_name = "probe-" + arg; })
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

    const std::string file_content = readFileContent(stem, ".cif");

    if(forceField.has_value())
    {
      // Case: force field has been manually selected
      if(const auto cif = CIFReader::readCIFString(file_content, forceField.value(), CIFReader::UseChargesFrom::CIF_File); cif.has_value())
      {
        auto [simulation_box, space_group_hall_symbol, defined_atoms, fractional_atoms_unit_cell] = cif.value();
        framework = Framework(forceField.value(), stem, simulation_box, space_group_hall_symbol,
                              defined_atoms, fractional_atoms_unit_cell, {1, 1, 1});
      }
      else if (cif.error() == CIFReader::ParseError::invalidForceField)
      {
        std::print("invalid forcefield\n");
        std::exit(-1);
      }
    }
    else
    {
      // Case: auto-detect force field
      // First: try with zeolite force field
      ForceField trial_zeolite_force_field = ForceField::makeZeoliteForceField(12.0, true, false, false);
      if(const auto zeolite_cif = CIFReader::readCIFString(file_content, trial_zeolite_force_field, CIFReader::UseChargesFrom::CIF_File); zeolite_cif.has_value())
      {
        forceField = trial_zeolite_force_field;
        auto [simulation_box, space_group_hall_symbol, defined_atoms, fractional_atoms_unit_cell] = zeolite_cif.value();
        framework = Framework(forceField.value(), stem, simulation_box, space_group_hall_symbol,
                              defined_atoms, fractional_atoms_unit_cell, {1, 1, 1});
      }
      else if (zeolite_cif.error() == CIFReader::ParseError::invalidForceField)
      {
        // Second: try with general MOF force field
        ForceField trial_mof_force_field = ForceField::makeMetalOrganicFrameworkForceField(12.0, true, false, false);
        if(const auto mof_cif = CIFReader::readCIFString(file_content, trial_mof_force_field, CIFReader::UseChargesFrom::CIF_File); mof_cif.has_value())
        {
          forceField = trial_mof_force_field;
          auto [simulation_box, space_group_hall_symbol, defined_atoms, fractional_atoms_unit_cell] = mof_cif.value();
          framework = Framework(forceField.value(), stem, simulation_box, space_group_hall_symbol,
                                defined_atoms, fractional_atoms_unit_cell, {1, 1, 1});
        }
        else if (mof_cif.error() == CIFReader::ParseError::invalidForceField)
        {
        }
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

            sa.run(forceField.value(), framework, well_depth_factor, probe_atom_name.value_or("probe-N2"), number_of_iterations, number_of_inner_steps);
          }

          if (use_gpu)
          {
            MC_OpenCL_SurfaceArea sa;
            sa.run(forceField.value(), framework, well_depth_factor, probe_atom_name.value_or("probe-N2"), number_of_iterations, number_of_inner_steps);
          }
        }

        if (use_integration_methods)
        {
          if (use_cpu)
          {
            Integration_SurfaceArea sa;
            sa.run(forceField.value(), framework, well_depth_factor, probe_atom_name.value_or("probe-N2"), number_of_slices);
          }

          if (use_gpu)
          {
            Integration_OpenCL_SurfaceArea sa;
            sa.run(forceField.value(), framework, well_depth_factor, probe_atom_name.value_or("probe-N2"), number_of_slices);
          }
        }
      }

      if (use_energy_methods)
      {
        if (use_cpu)
        {
          EnergySurfaceArea sa;
          sa.run(forceField.value(), framework, iso_value, probe_atom_name.value_or("probe-N2"));
        }

        if (use_gpu)
        {
          EnergyOpenCLSurfaceArea sa;
          sa.run(forceField.value(), framework, iso_value, probe_atom_name.value_or("probe-N2"), gridSize);
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
          //vf.run(forceField.value(), framework, well_depth_factor, probe_atom_name.value_or("probe-He"), number_of_iterations, number_of_inner_steps);
        }

        if (use_gpu)
        {
          EnergyOpenCLVoidFraction vf;
          //vf.run(forceField.value(), framework);
        }
      }

      if (use_energy_methods)
      {
        if (use_cpu)
        {
          EnergyVoidFraction vf;
          vf.run(forceField.value(), framework, probe_atom_name.value_or("probe-He"), number_of_iterations, number_of_inner_steps);
        }

        if (use_gpu)
        {
          EnergyOpenCLVoidFraction vf;
          vf.run(forceField.value(), framework, probe_atom_name.value_or("probe-He"), gridSize);
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
          psd.run(forceField.value(), framework, well_depth_factor, number_of_iterations, number_of_inner_steps, maximum_range);
        }

        if (use_gpu)
        {
          MC_OpenCL_PoreSizeDistribution psd(1000);
          psd.run(forceField.value(), framework, well_depth_factor, number_of_iterations, number_of_inner_steps, maximum_range);
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
