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
import structure_input;
import crystal;
import pair_interactions;
import sampled_structure;
import mc_void_fraction;
import mc_surface_area;
import mc_pore_size_distribution;
import mc_pore_diameters;
import mc_pore_volume;
import mc_channels;
import mc_window_shape;
import mc_blocking_pockets;
import sampled_roadmap;
import sampling_backend;
import mc_backend;
import mc_opencl_backend;
import mc_opencl_void_fraction;
import mc_opencl_surface_area;
import mc_opencl_pore_size_distribution;
import energy_opencl_void_fraction;
import energy_opencl_surface_area;
import energy_void_fraction;
import integration_surface_area;
import integration_opencl_surface_area;
import getopt;
import interpolation_energy_grid;
import pore_size_distribution_ban_vlugt;
import opencl_clearance_grid;
import grid_connected_components;
import opencl_pore_analysis;
import opencl_void_fraction;
import opencl_surface_area;
import opencl_pore_size_distribution;
import opencl_blocking_spheres;
import energy_backend;
import energy_opencl_backend;
import energy_opencl_probe_energy_grid;
import energy_shared_energy_barrier;
import energy_shared_linear_probe;
import energy_shared_molecular_energy_barrier;
import energy_shared_molecular_void_fraction;
import energy_shared_molecular_surface_area;
import energy_shared_pore_size_distribution;
import energy_shared_pore_analysis;
import energy_shared_pore_volume;
import energy_shared_blocking_spheres;
import energy_shared_tessellation;
import units;
import voronoi_pore_diameters;
import voronoi_channels;
import voronoi_surface_area;
import voronoi_accessible_volume;
import voronoi_blocking_spheres;
import voronoi_pore_size_distribution;
import apollonius_pore_analysis;
import apollonius_pore_size_distribution;
import apollonius_surface_area;
import apollonius_accessible_volume;
import apollonius_blocking_spheres;
import brute_force_validation;
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
  bool use_voronoi{false};
  bool use_apollonius{false};
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
  double temperature{ 298.0 };
  double blocking_threshold{ 30.0 };
  double brute_force_spacing{ 0.15 };
  std::size_t brute_force_samples{ 20000 };
  std::size_t brute_force_points{ 4000000 };
  bool brute_force_skip_excluded{ false };
  std::optional<std::string> molecule_name{};
  std::size_t number_of_orientations{ 128 };
  CIFReader::UseChargesFrom charges_from{CIFReader::UseChargesFrom::CIF_File};
  uint3 gridSize{128, 128, 128};
  std::optional<std::size_t> number_of_slices{ };
  std::optional<std::size_t> number_of_bins{ };
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
           [&number_of_iterations](std::string const &arg) { number_of_iterations = std::stoul(arg); })
      .reg({"-M", "--number-of-inner-steps"},
           "NUM_INNER_ITERATIONS",  // will be used in help to illustrate the argument
           argparser::required_argument,
           "Set number of inner steps",  // will be displayed in help
           [&number_of_inner_steps](std::string const &arg) { number_of_inner_steps = std::stoul(arg); })
      .reg({"--threads"},
           "NUM_THREADS",  // will be used in help to illustrate the argument
           argparser::required_argument,
           "Set number of threads the analyses may run on (default 1, one thread)",  // will be displayed in help
           [&num_threads](std::string const &arg) { num_threads = std::stoul(arg); })
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
      .reg({"--pore-analysis"}, argparser::no_argument,
           "Compute the pore diameters Di, Df and Dif and the channel/pocket analysis (as zeo++ -res and -chan)",
           [&state](std::string const &) { state.set(State::PoreAnalysis); })
      .reg({"--blocking-spheres"}, argparser::no_argument,
           "Compute spheres covering the pockets the probe cannot reach, in RASPA .block format",
           [&state](std::string const &) { state.set(State::BlockingSpheres); })
      .reg({"--energy-barrier"}, argparser::no_argument,
           "Compute the lowest energy at which the probe percolates, the energetic counterpart of Df",
           [&state](std::string const &) { state.set(State::PercolationBarrier); })
      .reg({"--brute-force"}, argparser::no_argument,
           "Work the exact geometric properties out again from the atom positions and radii alone, on a grid "
           "and with random points, and report the two side by side",
           [&state](std::string const &) { state.set(State::BruteForce); })
      .reg({"--brute-force-spacing"}, argparser::required_argument,
           "Spacing of the grid the brute-force check floods the void on [Å], default 0.15. What Df and Dif "
           "can resolve is set by this",
           [&brute_force_spacing](std::string const &arg) { brute_force_spacing = std::stod(arg); })
      .reg({"--brute-force-samples"}, argparser::required_argument,
           "Directions thrown at each atom by the brute-force check, default 20000. The surface areas converge "
           "as one over the root of this",
           [&brute_force_samples](std::string const &arg) { brute_force_samples = std::stoul(arg); })
      .reg({"--brute-force-points"}, argparser::required_argument,
           "Points thrown at the cell by the brute-force check, default 4000000",
           [&brute_force_points](std::string const &arg) { brute_force_points = std::stoul(arg); })
      .reg({"--brute-force-skip-excluded"}, argparser::no_argument,
           "Leave the convex/saddle/concave decomposition out of the brute-force check, which is the slowest "
           "part of it",
           [&brute_force_skip_excluded](std::string const &) { brute_force_skip_excluded = true; })
      .reg({"--temperature"}, argparser::required_argument,
           "Set the temperature the energy barrier is reported against [K], default 298",
           [&temperature](std::string const &arg) { temperature = std::stod(arg); })
      .reg({"--molecule"}, argparser::required_argument,
           "Use a rigid molecule with a shape of its own, such as CO2, for the energy barrier, the void "
           "fraction or the surface area",
           [&molecule_name](std::string const &arg) { molecule_name = arg; })
      .reg({"--orientations"}, argparser::required_argument,
           "Number of orientations to average a molecular property over, default 128. An area needs more of "
           "them than a barrier does",
           [&number_of_orientations](std::string const &arg) { number_of_orientations = std::stoul(arg); })
      .reg({"--charges-from"}, argparser::required_argument,
           "Where framework charges come from: 'cif' (default) or 'pseudo-atoms'. A CIF without a charge "
           "column leaves every atom neutral, so a charged guest needs 'pseudo-atoms' to feel anything",
           [&charges_from](std::string const &arg)
           {
             if (arg == "cif") charges_from = CIFReader::UseChargesFrom::CIF_File;
             else if (arg == "pseudo-atoms") charges_from = CIFReader::UseChargesFrom::PseudoAtoms;
             else throw std::runtime_error(std::format("Unknown charge source '{}', expected 'cif' or 'pseudo-atoms'\n", arg));
           })
      .reg({"--voronoi"}, argparser::no_argument,
           "Use the Voronoi (radical) pore network for the geometric analyses",
           [&use_voronoi](std::string const &) { use_voronoi = true; })
      .reg({"-apollonius", "--apollonius"}, argparser::no_argument,
           "Use the Apollonius diagram instead of the Voronoi network for the geometric analyses (on its own, "
           "computes the pore analysis)",
           [&use_apollonius](std::string const &) { use_apollonius = true; })
      .reg({"--tessellation"}, argparser::no_argument, "Use tessellation method",
           [&state](std::string const &) { state.set(State::TessellationComputation); })
      .reg({"--zeolite"}, argparser::no_argument, "Use generic zeolite model (TraPPE zeo)",
           [&forceField](std::string const &)
           { forceField = ForceField::makeZeoliteForceField(12.0, true, false, false); })
      .reg({"--mof"}, argparser::no_argument, "Use generic MOF model (TraPPE zeo)",
           [&forceField](std::string const &)
           { forceField = ForceField::makeMetalOrganicFrameworkForceField(12.0, true, false, false); })
      .reg({"--zeo++"}, argparser::no_argument, "Use zeo++ radii",
           [&forceField](std::string const &)
           { forceField = ForceField::makeZeoPlusPlusForceField(12.0, true, false, false); })
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
      .reg({"--blocking-threshold"},
           argparser::required_argument,
           "How many kT a region must cost a molecule to leave before the energy route blocks it (default 30, "
           "which holds it about ten seconds at a lattice frequency of 1e12 per second)",
           [&blocking_threshold](std::string const &arg) { blocking_threshold = std::stod(arg); })
      .reg({"--number-of-slices"}, argparser::required_argument,
           "Set number of slices: of the energy integration, or of each smooth piece of the exact surface area",
           [&number_of_slices](std::string const &arg)
           {
             std::istringstream iss(arg);
             std::size_t x;
             if (!(iss >> x))
             {
               throw std::runtime_error("Invalid --number-of-slices: expected a single integer");
             }
             number_of_slices = x;
           })
      .reg({"--number-of-bins"}, argparser::required_argument,
           "Set the number of pore sizes the pore-size distribution is evaluated at",
           [&number_of_bins](std::string const &arg)
           {
             std::istringstream iss(arg);
             std::size_t x;
             if (!(iss >> x))
             {
               throw std::runtime_error("Invalid --number-of-bins: expected a single integer");
             }
             number_of_bins = x;
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

  // The energy-based properties are the same arithmetic on a field whichever machine filled the field in, so
  // the backend is chosen once here and handed down. Asking for the other one changes nothing but the
  // precision of the field and the time it takes to build.
  auto energyBackend = [&] { return use_gpu ? openCLEnergyBackend() : cpuEnergyBackend(); };

  // Serial unless threads were asked for. An analysis that can use them looks the pool up and finds it empty
  // otherwise, and takes the route it has always taken; see `exact_parallel` for what the threaded routes are
  // and are not promised to agree with.
  ThreadPool::ThreadPool<ThreadPool::details::default_function_type, std::jthread>::instance().init(
      num_threads, num_threads > 1 ? ThreadPool::ThreadingType::ThreadPool : ThreadPool::ThreadingType::Serial);

  // A pore network is only of use to an analysis that reads one, so asking for the Apollonius diagram
  // without saying what to do with it is asking for the pore analysis.
  if ((use_apollonius || use_voronoi) && !state.test(State::SurfaceArea) && !state.test(State::VoidFraction) &&
      !state.test(State::BlockingSpheres))
  {
    state.set(State::PoreAnalysis);
  }

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
      if(const auto cif = CIFReader::readCIFString(file_content, forceField.value(),
                                                   charges_from); cif.has_value())
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
      if(const auto zeolite_cif = CIFReader::readCIFString(file_content, trial_zeolite_force_field,
                                                           charges_from);
         zeolite_cif.has_value())
      {
        forceField = trial_zeolite_force_field;
        auto [simulation_box, space_group_hall_symbol, defined_atoms, fractional_atoms_unit_cell] = zeolite_cif.value();
        framework = Framework(forceField.value(), stem, simulation_box, space_group_hall_symbol,
                              defined_atoms, fractional_atoms_unit_cell, {1, 1, 1});
      }
      else if (zeolite_cif.error() == CIFReader::ParseError::invalidForceField)
      {
        // Second: try with general MOF force field
        ForceField trial_mof_force_field =
            ForceField::makeMetalOrganicFrameworkForceField(12.0, true, false, false);
        if(const auto mof_cif = CIFReader::readCIFString(file_content, trial_mof_force_field,
                                                        CIFReader::UseChargesFrom::CIF_File); mof_cif.has_value())
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

    // The structural analysis is a library of its own, built on the same foundations as the simulation
    // engine rather than on the engine, so it is handed neither a framework nor a force field. This is where
    // the two become what it does take: a cell with a set of typed and charged atom centres in it, and a
    // table of what a pair of types does to one another. It happens once, after any custom probe has been
    // put into the force field above, so that everything below measures the same structure.
    Crystal crystal = StructureInput::makeCrystal(framework);
    PairInteractions interactions = StructureInput::makeInteractions(forceField.value());

    // The probe the geometric analyses sample with. Not every force field defines every probe, a
    // custom one read from a file need define none of them, so an unchosen default falls back to
    // nitrogen rather than aborting the run over a name the user never asked for. A probe named on
    // the command line is passed through, and the analysis says so if it does not exist.
    auto geometricProbe = [&](const std::string &preferred) -> std::string
    {
      if (probe_atom_name.has_value()) return probe_atom_name.value();
      if (forceField->findPseudoAtom(preferred).has_value()) return preferred;
      return "probe-N2";
    };

    // The molecule the energy routes answer for. A name given with --molecule is looked for among the ones
    // with a shape of their own first and among the pseudo-atoms after, so that a single site can be sent
    // the long way round and held against the direct route; with no name it is the probe atom, as one site.
    auto energyMolecule = [&](const std::string &what, const std::string &preferred) -> LinearProbe
    {
      std::string name = molecule_name.value_or(geometricProbe(preferred));
      std::optional<LinearProbe> molecule = LinearProbe::named(interactions, name);
      if (!molecule.has_value()) molecule = LinearProbe::singleSite(interactions, name);
      if (!molecule.has_value())
      {
        throw std::runtime_error(std::format("Unknown molecule '{}' for the {}\n", name, what));
      }
      return molecule.value();
    };

    // A molecule with a shape has to be averaged over its orientations; a single site has only the one.
    auto energyOrientations = [&] { return molecule_name.has_value() ? number_of_orientations : std::size_t{1}; };

    // The sampling routines go further still and ask for no types at all, only a cell, a set of atom centres
    // and the radius of each atom's contact sphere. This is where the mixing rule becomes those radii, once
    // and in one place rather than repeated in each of them.
    auto sampledStructure = [&](std::vector<double> radii) -> SampledStructure
    {
      return SampledStructure{.name = crystal.name,
                              .spaceGroupHallNumber = crystal.spaceGroupHallNumber,
                              .unitCell = crystal.unitCell,
                              .positions = crystal.cartesianPositions(),
                              .radii = std::move(radii),
                              .mass = crystal.mass};
    };

    // For the routines that measure the wall: a point is in contact when the probe's centre is a mixed
    // diameter from the atom's, so that is the radius of each atom's sphere. An unknown probe is refused
    // here rather than after the structure has been assembled.
    auto sampledSurface = [&](const std::string &probeName) -> std::pair<SampledStructure, SampledProbe>
    {
      std::optional<std::size_t> probeType = interactions.findType(probeName);
      if (!probeType.has_value())
      {
        throw std::runtime_error(std::format("Unknown probe atom '{}' for the surface area\n", probeName));
      }

      std::vector<double> radii;
      radii.reserve(crystal.size());
      for (const CrystalAtom &atom : crystal.atoms)
      {
        radii.push_back(well_depth_factor * interactions(probeType.value(), atom.type).sizeParameter);
      }

      SampledProbe probe{.name = probeName,
                         .sizeParameter = interactions[probeType.value()].sizeParameter,
                         .wellDepthFactor = well_depth_factor};

      return {sampledStructure(std::move(radii)), probe};
    };

    // For the routines that measure the void: the probe is a point there and what it is kept out of is the
    // atom's own sphere, so the radius is half the atom's own diameter.
    auto sampledVoid = [&]() -> SampledStructure
    {
      std::vector<double> radii;
      radii.reserve(crystal.size());
      for (const CrystalAtom &atom : crystal.atoms)
      {
        radii.push_back(0.5 * well_depth_factor * interactions[atom.type].sizeParameter);
      }

      return sampledStructure(std::move(radii));
    };

    // And for the routines that ask what a probe of a size of its own can reach. The probe's centre is kept
    // the two radii apart from the atom's centre, which under the Lorentz-Berthelot rule is the mixed
    // diameter of the pair -- the same inflation the diagram routes apply to their atoms before building
    // anything, so the two are measuring the same void.
    auto sampledVoidFor = [&](const std::string &probeName) -> SampledStructure
    {
      std::optional<std::size_t> probeType = interactions.findType(probeName);
      if (!probeType.has_value())
      {
        throw std::runtime_error(std::format("Unknown probe atom '{}' for the sampled void\n", probeName));
      }

      std::vector<double> radii;
      radii.reserve(crystal.size());
      for (const CrystalAtom &atom : crystal.atoms)
      {
        radii.push_back(well_depth_factor * interactions(probeType.value(), atom.type).sizeParameter);
      }

      return sampledStructure(std::move(radii));
    };

    if (state.test(CommandLine::State::TessellationComputation))
    {
      std::cout << "Compute tesselation" << std::endl;

      if (use_gridbased_methods)
      {
        if (use_gpu)
        {
          // The clearance field already names the nearest atom at every point, that being what its distance is
          // measured to, so the tessellation is a report on the field rather than a computation of its own.
          ClearanceGrid grid = ClearanceGrid::compute(interactions, crystal, gridSize);
          grid.writeTessellation(crystal);
        }
      }

      if (use_energy_methods)
      {
        // The same division, by strongest attraction rather than by nearest surface. It depends on the probe
        // in a way the geometric one does not, so which probe was used is part of the answer and goes in the
        // name of the file.
        EnergyTessellation tessellation;
        if (molecule_name.has_value())
        {
          tessellation.run(energyBackend(), interactions, crystal,
                           energyMolecule("tessellation", "probe-He"), iso_value, gridSize, energyOrientations(),
                           temperature);
        }
        else
        {
          tessellation.run(energyBackend(), interactions, crystal, geometricProbe("probe-He"), iso_value,
                           gridSize);
        }
        std::cout << "iso-surface " << tessellation.totalGravimetricArea << " m^2/g divided among "
                  << tessellation.atoms.size() << " atoms, " << tessellation.undecidedArea
                  << " A^2 of it in no atom's cell" << std::endl;
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

      if (use_gridbased_methods && use_gpu)
      {
        std::cout << "Compute the accessible surface area from the clearance grid" << std::endl;

        GridSurfaceArea sa;
        sa.run(interactions, crystal, geometricProbe("probe-N2"), gridSize);
      }

      if (use_apollonius)
      {
        // The area itself is measured rather than sampled unless the sampled estimate is asked for by
        // name, there being no reason to prefer a statistical answer to an exact one at the same cost.
        ApolloniusSurfaceArea sa;
        sa.run(interactions, crystal, geometricProbe("probe-N2"),
               use_monte_carlo_methods ? ApolloniusSurfaceArea::Method::Sampled
                                       : ApolloniusSurfaceArea::Method::Exact,
               number_of_iterations, number_of_slices);
      }
      else if (use_voronoi)
      {
        // The area is measured against the radical network too: it is the same union of atoms and comes
        // to the same total, and what the network then decides is only how that total divides.
        VoronoiSurfaceArea sa;
        sa.run(interactions, crystal, geometricProbe("probe-N2"),
               use_monte_carlo_methods ? VoronoiSurfaceArea::Method::Sampled
                                       : VoronoiSurfaceArea::Method::Exact,
               number_of_iterations, number_of_slices);
      }

      if (use_geometric_methods)
      {
        auto [structure, probe] = sampledSurface(probe_atom_name.value_or("probe-N2"));

        if (use_monte_carlo_methods)
        {
          if (use_cpu)
          {
            MC_SurfaceArea sa;
            sa.run(structure, probe, number_of_iterations, number_of_inner_steps);
          }

          if (use_gpu)
          {
            MC_OpenCL_SurfaceArea sa;
            sa.run(structure, probe, number_of_iterations, number_of_inner_steps);
          }
        }

        if (use_integration_methods)
        {
          if (use_cpu)
          {
            Integration_SurfaceArea sa;
            sa.run(structure, probe, number_of_slices);
          }

          if (use_gpu)
          {
            Integration_OpenCL_SurfaceArea sa;
            sa.run(structure, probe, number_of_slices);
          }
        }
      }

      if (use_energy_methods)
      {
        // As with the void fraction: a molecule with a shape of its own goes through the orientational
        // landscape, anything else is a single site and takes the cheaper route.
        std::optional<LinearProbe> molecule{};
        if (molecule_name.has_value())
        {
          molecule = LinearProbe::named(interactions, molecule_name.value());
          if (!molecule.has_value()) molecule = LinearProbe::singleSite(interactions, molecule_name.value());
          if (!molecule.has_value())
          {
            throw std::runtime_error(
                std::format("Unknown molecule '{}' for the surface area\n", molecule_name.value()));
          }
        }

        if (molecule.has_value())
        {
          MolecularSurfaceArea sa;
          sa.run(energyBackend(), interactions, crystal, molecule.value(), iso_value, gridSize,
                 number_of_orientations, temperature);
          std::cout << "free-energy surface: " << sa.gravimetricArea << " m^2/g, minimum-energy surface: "
                    << sa.minimumEnergyGravimetricArea << " m^2/g (orientational excess "
                    << sa.orientationalExcess() << ")" << std::endl;
        }
        else if (use_gpu)
        {
          EnergyOpenCLSurfaceArea sa;
          sa.run(interactions, crystal, iso_value, probe_atom_name.value_or("probe-N2"), gridSize);
        }
        else
        {
          // The single-site route has no processor version of its own, so the probe is sent round as a
          // molecule of one site. The two agree exactly, which is what makes the substitution honest.
          std::string probe = probe_atom_name.value_or("probe-N2");
          std::optional<LinearProbe> single = LinearProbe::singleSite(interactions, probe);
          if (!single.has_value())
          {
            throw std::runtime_error(std::format("Unknown probe atom '{}' for the surface area\n", probe));
          }

          MolecularSurfaceArea sa;
          sa.run(energyBackend(), interactions, crystal, single.value(), iso_value, gridSize, 1, temperature);
          std::cout << "surface: " << sa.gravimetricArea << " m^2/g" << std::endl;
        }
      }
    }

    if (state.test(CommandLine::State::VoidFraction))
    {
      std::cout << "Compute void fraction" << std::endl;

      // The geometric void fraction is the probe-accessible void, split from the void closed off in
      // pockets, which is what zeo++ reports as the accessible volume.
      if (use_gridbased_methods && use_gpu)
      {
        std::cout << "Compute the void volume from the clearance grid" << std::endl;

        GridVoidFraction vf;
        vf.run(interactions, crystal, geometricProbe("probe-He"), gridSize);
      }

      if (use_apollonius)
      {
        // Both how much void there is and how it divides between channels and pockets are measured
        // rather than sampled, unless the sampled estimate is asked for by name. The division falls
        // back on sampling only where the surface cannot supply it, and says so when it does.
        ApolloniusAccessibleVolume av;
        av.run(interactions, crystal, geometricProbe("probe-He"),
               use_monte_carlo_methods ? ApolloniusAccessibleVolume::Method::Sampled
                                       : ApolloniusAccessibleVolume::Method::Exact,
               number_of_iterations, number_of_slices);
      }
      else if (use_voronoi)
      {
        VoronoiAccessibleVolume av;
        av.run(interactions, crystal, geometricProbe("probe-He"),
               use_monte_carlo_methods ? VoronoiAccessibleVolume::Method::Sampled
                                       : VoronoiAccessibleVolume::Method::Exact,
               number_of_iterations, number_of_slices);
      }

      if (use_monte_carlo_methods && !use_energy_methods)
      {
        // Throwing points at the cell and counting the ones that land clear of every atom. The fraction
        // needs nothing else and is unbiased; splitting it into the part a probe can reach and the part
        // shut away needs to know what is connected to what, which is what the roadmap is for.
        std::cout << "Compute the void fraction and pore volume by sampling the void" << std::endl;

        SamplingBackend sampling = use_gpu ? samplingBackendOpenCL() : samplingBackendCPU();
        SampledStructure probedVoid = sampledVoidFor(geometricProbe("probe-He"));

        MC_PoreVolume volume;
        volume.run(probedVoid,
                   SampledRoadmap::build(probedVoid, sampling, number_of_iterations, number_of_inner_steps));

        std::cout << "void fraction " << volume.voidFraction << " +/- " << volume.voidFractionError << ", of which "
                  << volume.accessibleVolumeFraction << " accessible and " << volume.inaccessibleVolumeFraction
                  << " shut in " << volume.numberOfPockets << " pocket(s)" << std::endl;
        std::cout << "pore volume " << volume.gravimetricVoidVolume << " cm^3/g, of which "
                  << volume.gravimetricAccessibleVolume << " cm^3/g a probe can reach" << std::endl;
      }

      if (use_energy_methods && use_monte_carlo_methods)
      {
        // Widom insertion samples the same average the grid route sums, at random points rather than on a
        // lattice, and --monte-carlo asks for the sampled estimate here as it does for the geometric
        // properties. It shares no arithmetic with the field routes, so its agreeing with them is worth
        // something.
        EnergyVoidFraction vf;
        vf.run(interactions, crystal, probe_atom_name.value_or("probe-He"), number_of_iterations,
               number_of_inner_steps);
      }
      else if (use_energy_methods)
      {
        // A molecule with a shape of its own goes through the orientational landscape; anything else is a
        // single site and takes the cheaper route, the two agreeing when the molecule named is one site.
        std::optional<LinearProbe> molecule{};
        if (molecule_name.has_value())
        {
          molecule = LinearProbe::named(interactions, molecule_name.value());
          if (!molecule.has_value()) molecule = LinearProbe::singleSite(interactions, molecule_name.value());
          if (!molecule.has_value())
          {
            throw std::runtime_error(std::format("Unknown molecule '{}' for the void fraction\n",
                                                 molecule_name.value()));
          }
        }

        if (molecule.has_value())
        {
          MolecularVoidFraction vf;
          vf.run(energyBackend(), interactions, crystal, molecule.value(), gridSize, number_of_orientations,
                 temperature);
          std::cout << "average of exp(-A/kT): " << vf.boltzmannAverage
                    << (vf.readsAsFraction ? " (reads as a void fraction)"
                                           : " (above one, so a Henry coefficient rather than a fraction)")
                    << ", K_H = " << vf.henryCoefficient << " mol/kg/Pa" << std::endl;
        }
        else if (use_gpu)
        {
          EnergyOpenCLVoidFraction vf;
          vf.run(interactions, crystal, probe_atom_name.value_or("probe-He"), gridSize, temperature);
        }
        else
        {
          // The single-site route has no processor version of its own, so the probe is sent round as a
          // molecule of one site. The two agree exactly, which is what makes the substitution honest.
          std::string probe = probe_atom_name.value_or("probe-He");
          std::optional<LinearProbe> single = LinearProbe::singleSite(interactions, probe);
          if (!single.has_value())
          {
            throw std::runtime_error(std::format("Unknown probe atom '{}' for the void fraction\n", probe));
          }

          MolecularVoidFraction vf;
          vf.run(energyBackend(), interactions, crystal, single.value(), gridSize, 1, temperature);
          std::cout << "average of exp(-U/kT): " << vf.boltzmannAverage
                    << (vf.readsAsFraction ? " (reads as a void fraction)"
                                           : " (above one, so a Henry coefficient rather than a fraction)")
                    << ", K_H = " << vf.henryCoefficient << " mol/kg/Pa" << std::endl;
        }

        // The average above is one reading of how much room there is and it is not a volume. Counting the
        // region the molecule can occupy is the other, and it is the one that comes out in Å³ and cm³/g and
        // can be set beside a geometric table. Both are read off the same landscape, which is built once.
        EnergyPoreVolume pv;
        pv.run(energyBackend(), interactions, crystal, energyMolecule("pore volume", "probe-He"), iso_value,
               gridSize, energyOrientations(), temperature, blocking_threshold);
        std::cout << "pore volume " << pv.gravimetricVoidVolume << " cm^3/g, of which "
                  << pv.gravimetricReachableVolume << " cm^3/g the molecule can reach; by Boltzmann weight "
                  << pv.gravimetricBoltzmannVolume << " cm^3/g" << std::endl;
      }
    }

    if (state.test(CommandLine::State::PSD))
    {
      // The distribution itself is a closed form over the surface of the framework, so it is evaluated rather
      // than sampled unless the sampled estimate is asked for by name. Which diagram is named decides only how
      // the curve is divided between the void a probe can reach and the void it cannot.
      //
      // The probe named is the one the accessible distribution is reported for, beside the distribution of the
      // whole of the void. Helium by default, as for the accessible volume: it is the molecule a void volume is
      // measured with, and being the smallest of the probes it is the one that separates the pores a molecule
      // cannot enter at all from the pores it merely finds narrow. A larger probe answers a narrower question,
      // and answers nothing whatever in a framework whose windows it cannot pass.
      if (use_gridbased_methods && use_gpu)
      {
        std::cout << "Compute the pore-size distribution from the clearance grid" << std::endl;

        GridPoreSizeDistribution psd;
        psd.run(interactions, crystal, geometricProbe("probe-He"), gridSize, maximum_range, number_of_bins);
      }

      if (use_apollonius)
      {
        std::cout << "Compute the pore-size distribution from the Apollonius diagram" << std::endl;
        ApolloniusPoreSizeDistribution psd;
        psd.run(interactions, crystal, geometricProbe("probe-He"), maximum_range, number_of_bins,
                number_of_slices.value_or(1));
      }
      else if (use_voronoi)
      {
        std::cout << "Compute the pore-size distribution from the radical (Voronoi) network" << std::endl;
        VoronoiPoreSizeDistribution psd;
        psd.run(interactions, crystal, geometricProbe("probe-He"), maximum_range, number_of_bins,
                number_of_slices.value_or(1));
      }

      if (use_monte_carlo_methods)
      {
        SampledStructure structure = sampledVoid();

        if (use_cpu)
        {
          MC_PoreSizeDistribution psd(1000);
          psd.run(structure, number_of_iterations, number_of_inner_steps, maximum_range);
        }

        if (use_gpu)
        {
          MC_OpenCL_PoreSizeDistribution psd(1000);
          psd.run(structure, number_of_iterations, number_of_inner_steps, maximum_range);
        }
      }

      if (use_energy_methods)
      {
        // As with the surface area: a molecule with a shape of its own goes through the orientational
        // landscape, and anything else is sent round as a molecule of one site, which the landscape
        // reproduces exactly.
        std::string probe = molecule_name.value_or(probe_atom_name.value_or("probe-N2"));
        std::optional<LinearProbe> molecule = LinearProbe::named(interactions, probe);
        if (!molecule.has_value()) molecule = LinearProbe::singleSite(interactions, probe);
        if (!molecule.has_value())
        {
          throw std::runtime_error(
              std::format("Unknown molecule '{}' for the pore-size distribution\n", probe));
        }

        std::size_t orientations = molecule_name.has_value() ? number_of_orientations : 1;

        EnergyPoreSizeDistribution psd;
        psd.run(energyBackend(), interactions, crystal, molecule.value(), iso_value, gridSize, orientations,
                temperature, blocking_threshold, true, 1e-6, maximum_range, number_of_bins);
        std::cout << "largest sphere the void holds: " << psd.largestDiameter << " A, void fraction at this level "
                  << psd.voidFraction << ", running in " << psd.dimensionality << " directions" << std::endl;
        std::cout << "mean pore per unit volume " << psd.volumetricMeanDiameter << " A, per molecule "
                  << psd.occupancyMeanDiameter << " A, with " << psd.reachableOccupancyFraction
                  << " of the molecules in reachable void" << std::endl;
      }
    }

    if (state.test(CommandLine::State::PoreAnalysis))
    {
      std::string probe = geometricProbe("probe-N2");

      bool onTheGrid = use_gridbased_methods && use_gpu;

      if (onTheGrid)
      {
        std::cout << "Compute pore diameters and channels from the clearance grid" << std::endl;

        GridPoreAnalysis analysis;
        analysis.run(interactions, crystal, probe, gridSize);
      }

      if (use_energy_methods)
      {
        LinearProbe molecule = energyMolecule("pore analysis", "probe-N2");
        std::cout << "Compute pore diameters and channels from the energy landscape for " << molecule.name
                  << std::endl;

        EnergyPoreAnalysis analysis;
        analysis.run(energyBackend(), interactions, crystal, molecule, iso_value, gridSize,
                     energyOrientations(), temperature);

        // The levels are the answers with no parameter in them, so they are the ones worth putting on the
        // screen; the lengths need the iso-value and are left to the file.
        std::cout << "Di " << analysis.levels.deepestWell * Units::EnergyToKelvin << " K";
        if (analysis.levels.percolates)
        {
          std::cout << ", Df " << analysis.levels.percolationBarrier * Units::EnergyToKelvin << " K ("
                    << analysis.levels.percolationBarrier / (Units::KB * temperature) << " kT), Dif "
                    << analysis.levels.deepestWellOnPath * Units::EnergyToKelvin << " K, running in "
                    << analysis.levels.dimensionalityWithin(Units::KB * temperature).first << " direction(s)";
        }
        else
        {
          std::cout << ", and no path through the crystal at any energy";
        }
        std::cout << std::endl;
        std::cout << "on the contour: Di " << analysis.diameters.includedSphereDiameter << " A, Df "
                  << analysis.diameters.freeSphereDiameter << " A, Dif "
                  << analysis.diameters.includedAlongFreePathDiameter << " A" << std::endl;
      }

      if (use_monte_carlo_methods)
      {
        std::cout << "Compute pore diameters, channels and window shapes by sampling the void" << std::endl;

        // The diameters and the windows are measured with a point probe, as the diagram routes measure
        // theirs; the channels are what the named probe can travel, so they are sampled against its own
        // contact radii. Each of the two roadmaps is built once and read three times over.
        SamplingBackend sampling = use_gpu ? samplingBackendOpenCL() : samplingBackendCPU();

        SampledStructure pointProbe = sampledVoid();
        SampledRoadmap voidMap =
            SampledRoadmap::build(pointProbe, sampling, number_of_iterations, number_of_inner_steps);

        MC_PoreDiameters diameters;
        diameters.run(pointProbe, voidMap);

        MC_WindowShape windows;
        windows.run(pointProbe, voidMap);

        SampledStructure probedVoid = sampledVoidFor(probe);
        MC_Channels channels;
        channels.run(probedVoid,
                     SampledRoadmap::build(probedVoid, sampling, number_of_iterations, number_of_inner_steps));

        // Every one of these is a lower bound, so it is worth saying what it was read from: a run whose
        // sample is too small says so by finding too few points in the void, not by scattering.
        std::cout << "Di " << diameters.result.includedSphereDiameter << " A";
        if (diameters.percolates)
        {
          std::cout << ", Df " << diameters.result.freeSphereDiameter << " A, Dif "
                    << diameters.result.includedAlongFreePathDiameter << " A";
        }
        else
        {
          std::cout << ", and no path through the crystal was found";
        }
        std::cout << std::endl;
        std::cout << "from " << voidMap.numberOfVoidSamples << " points in the void of " << voidMap.numberOfSamples
                  << " thrown (void fraction " << voidMap.voidFraction << "), " << voidMap.numberOfLinks
                  << " links, " << voidMap.numberOfPocketCentres << " pocket centres, on " << voidMap.backend
                  << std::endl;

        if (windows.freeSphere.measured)
        {
          std::cout << "Df window: free width " << windows.freeSphere.smallestFreeWidth << " - "
                    << windows.freeSphere.largestFreeWidth << " A, ellipse " << windows.freeSphere.minorAxis
                    << " x " << windows.freeSphere.majorAxis << " A, ring of "
                    << windows.freeSphere.boundingAtoms << " atoms" << std::endl;
        }

        std::cout << channels.numberOfChannels << " channel(s) and " << channels.numberOfPockets
                  << " pocket(s) for " << probe << ", running in " << channels.dimensionality
                  << " direction(s), holding " << channels.channelShareOfVoid << " of the void" << std::endl;
      }

      if (use_apollonius)
      {
        std::cout << "Compute pore diameters and channels from the Apollonius diagram" << std::endl;

        ApolloniusPoreAnalysis analysis;
        analysis.run(interactions, crystal, probe);
      }
      else if (!onTheGrid && (use_voronoi || !(use_energy_methods || use_monte_carlo_methods)))
      {
        // The Voronoi route is what a run that named no route at all gets. Asking for the energy one counts
        // as naming a route, so it is not also given this one behind its back.
        std::cout << "Compute pore diameters and channels from the Voronoi network" << std::endl;

        VoronoiPoreDiameters diameters;
        diameters.run(interactions, crystal);

        // Both analyses read the same network, and building it costs more than either of them.
        VoronoiChannels channels;
        channels.run(interactions, crystal, probe, diameters.network);
      }
    }

    if (state.test(CommandLine::State::BlockingSpheres))
    {
      std::string probe = geometricProbe("probe-N2");

      // The spheres come from the surfaces of the pockets themselves; the network is there for the cluster the
      // surfaces cannot place and for the sampled fallback, so it is worth saying which of the two was used.
      auto report = [](std::size_t spheres, bool measured, const std::string &reason)
      {
        std::cout << spheres << (measured ? " spheres, one per pocket, measured from its own surface"
                                          : " spheres, sampled: " + reason)
                  << std::endl;
      };

      bool onTheGrid = use_gridbased_methods && use_gpu;

      if (onTheGrid)
      {
        std::cout << "Compute blocking spheres from the clearance grid" << std::endl;

        GridBlockingSpheres blocks;
        blocks.run(interactions, crystal, probe, gridSize);
        std::cout << blocks.spheres.size() << " spheres, covering " << blocks.numberOfPockets << " pockets"
                  << std::endl;
      }

      if (use_energy_methods)
      {
        LinearProbe molecule = energyMolecule("blocking spheres", "probe-N2");
        std::cout << "Compute blocking spheres from the energy landscape for " << molecule.name << std::endl;

        EnergyBlockingSpheres blocks;
        blocks.run(energyBackend(), interactions, crystal, molecule, iso_value, gridSize,
                   energyOrientations(), temperature, blocking_threshold);

        std::cout << blocks.spheres.size() << " spheres over " << blocks.numberOfBlockedCavities
                  << " pieces too deep to leave, at " << blocking_threshold << " kT" << std::endl;
        if (blocks.numberOfLeakyCavities > 0)
        {
          std::cout << blocks.numberOfLeakyCavities << " further piece(s), holding " << blocks.leakyFraction
                    << " of the cell, do not run anywhere but are left open: the molecule gets out of them"
                    << std::endl;
        }
      }

      if (use_monte_carlo_methods)
      {
        // The pockets are the pieces of the roadmap that run nowhere, and the spheres are grown over the
        // points in them. Unlike the surface routes this one can miss a cavity no point landed in, so what
        // it found is printed beside what it wrote.
        std::cout << "Compute blocking spheres by sampling the void" << std::endl;

        SamplingBackend sampling = use_gpu ? samplingBackendOpenCL() : samplingBackendCPU();
        SampledStructure probedVoid = sampledVoidFor(probe);

        MC_BlockingPockets blocks;
        blocks.run(probedVoid,
                   SampledRoadmap::build(probedVoid, sampling, number_of_iterations, number_of_inner_steps));

        std::cout << blocks.spheres.size() << " spheres over " << blocks.numberOfCoveredPockets << " of "
                  << blocks.numberOfPockets << " pocket(s), holding " << blocks.blockedFractionOfVoid
                  << " of the sampled void" << std::endl;
      }

      if (use_apollonius)
      {
        std::cout << "Compute blocking spheres from the Apollonius diagram" << std::endl;

        ApolloniusBlockingSpheres blocks;
        blocks.run(interactions, crystal, probe);
        report(blocks.spheres.size(), blocks.measured, blocks.fallbackReason);
      }
      else if (!onTheGrid && !use_monte_carlo_methods && (use_voronoi || !use_energy_methods))
      {
        // As with the pore analysis: this is the route a run that named none gets, and asking for the energy
        // one is naming one.
        std::cout << "Compute blocking spheres from the Voronoi network" << std::endl;

        VoronoiBlockingSpheres blocks;
        blocks.run(interactions, crystal, probe);
        report(blocks.spheres.size(), blocks.measured, blocks.fallbackReason);
      }
    }

    if (state.test(CommandLine::State::BruteForce))
    {
      std::cout << "Check the exact geometry against brute force" << std::endl;

      BruteForceSettings settings{.spacing = brute_force_spacing,
                                  .samplesPerAtom = brute_force_samples,
                                  .volumePoints = brute_force_points,
                                  .subdivisions = number_of_slices.value_or(1),
                                  .skipSolventExcluded = brute_force_skip_excluded};

      // The same probes the routes being checked use by default, so that what comes out is comparable with
      // what those routes wrote rather than with a third convention invented here.
      BruteForceValidation validation;
      validation.run(interactions, crystal, geometricProbe("probe-N2"), geometricProbe("probe-He"),
                     settings);

      std::size_t disagreements = validation.numberOfDisagreements();
      std::cout << validation.checks.size() - disagreements << " of " << validation.checks.size()
                << " checks agree";
      if (disagreements > 0)
      {
        std::cout << "; see " << framework.name << ".brute-force.txt for which";
      }
      std::cout << std::endl;
    }

    if (state.test(CommandLine::State::PercolationBarrier))
    {
      // The barrier is a property of a molecule in a framework rather than of the framework alone, so the
      // probe is the one the energy routes use rather than the geometric stand-in, and it is not swapped for
      // a bare radius when the run is asked for zeo++ conventions.
      // A molecule named here is one with a shape of its own, and it is carried through the orientational
      // route. Anything else is a single site and goes the cheaper way, the two agreeing exactly when the
      // molecule named happens to be a single site.
      // A name is looked for among the molecules with a shape first, and failing that among the pseudo-atoms,
      // so that a single site can be sent the long way round. That is not only a convenience: a single site
      // has to come back from the orientational route with exactly what the direct route gives, and being
      // able to ask for it is what makes that comparable.
      std::optional<LinearProbe> molecule{};
      if (molecule_name.has_value())
      {
        molecule = LinearProbe::named(interactions, molecule_name.value());
        if (!molecule.has_value()) molecule = LinearProbe::singleSite(interactions, molecule_name.value());
      }

      if (molecule_name.has_value() && !molecule.has_value())
      {
        throw std::runtime_error(
            std::format("Unknown molecule '{}' for the energy barrier; the ones with a shape of their own are {}, "
                        "and any pseudo-atom may be named with --probe-atom-name instead\n",
                        molecule_name.value(), [] {
                          std::string names;
                          for (const std::string &name : LinearProbe::builtInNames())
                            names += (names.empty() ? "" : ", ") + name;
                          return names;
                        }()));
      }

      if (molecule.has_value())
      {
        std::cout << "Compute the percolation barrier for " << molecule->name << " over " << number_of_orientations
                  << " orientations" << std::endl;

        MolecularEnergyBarrier barrier;
        barrier.run(energyBackend(), interactions, crystal, molecule.value(), gridSize,
                    number_of_orientations, temperature);

        if (barrier.fromFreeEnergy.percolates)
        {
          std::cout << "barrier " << barrier.fromFreeEnergy.percolationBarrier * Units::EnergyToKelvin << " K ("
                    << barrier.fromFreeEnergy.percolationBarrier / (Units::KB * temperature) << " kT at "
                    << temperature << " K), running in "
                    << barrier.fromFreeEnergy.dimensionalityWithin(Units::KB * temperature).first
                    << " direction(s); orientation costs "
                    << barrier.orientationalPenalty * Units::EnergyToKelvin << " K of it" << std::endl;
        }
        if (barrier.grid.chargesIgnored)
        {
          std::cout << "warning: " << molecule->name
                    << " carries partial charges that this landscape does not act on" << std::endl;
        }
        else if (barrier.grid.chargesIncluded && barrier.potential.largestFrameworkCharge == 0.0)
        {
          std::cout << "warning: every framework atom is neutral, so the electrostatics did nothing; pass "
                       "--charges-from pseudo-atoms to take charges from the force field"
                    << std::endl;
        }
      }
      else
      {
        std::string probe = probe_atom_name.value_or("probe-N2");

        std::cout << "Compute the percolation barrier from the probe energy grid" << std::endl;

        EnergyBarrier barrier;
        barrier.run(energyBackend(), interactions, crystal, probe, gridSize, temperature);

        if (barrier.percolates)
        {
          std::cout << "barrier " << barrier.percolationBarrier * Units::EnergyToKelvin << " K ("
                    << barrier.percolationBarrier / (Units::KB * temperature) << " kT at " << temperature
                    << " K), running in " << barrier.dimensionalityWithin(Units::KB * temperature).first
                    << " direction(s)" << std::endl;
        }
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
