module;

export module run_simulation;

import std;

import input_reader;

/**
 * \brief Runs the simulation described by a parsed input file.
 *
 * Constructs the driver belonging to the simulation type in the input, restores it from the binary
 * restart file when the input asks for it, and runs it to completion. Which drivers have to reopen
 * their output files and rebuild their interpolation grids after such a restart is a property of the
 * drivers themselves, so that knowledge lives here rather than in each program that starts a run.
 *
 * The thread pool is initialized from the number of threads and the threading type in the input, so
 * callers do not have to do so themselves.
 *
 * \param inputReader The parsed input file. The driver takes over its systems.
 */
export void runSimulation(InputReader& inputReader);
