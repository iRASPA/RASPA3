module;

#include <csignal>

export module graceful_shutdown;

import std;

// Cooperative shutdown on SIGTERM/SIGINT (and SIGUSR1 on POSIX, commonly sent by HPC schedulers
// ahead of the walltime kill): the signal handler only sets a flag; the drivers poll it at cycle
// boundaries (where the simulation state is consistent), write a final binary restart file, and
// exit cleanly. This turns a scheduler kill into a checkpoint at the kill point instead of losing
// up to 'writeBinaryRestartEvery' cycles of work.

namespace
{
std::atomic<bool> stopRequestedFlag{false};

extern "C" void gracefulShutdownSignalHandler(int)
{
  // only async-signal-safe operations are allowed here
  stopRequestedFlag.store(true, std::memory_order_relaxed);
}
}  // namespace

export namespace GracefulShutdown
{
/**
 * \brief Installs the shutdown signal handlers; called once from main.
 */
void install()
{
  std::signal(SIGTERM, gracefulShutdownSignalHandler);
  std::signal(SIGINT, gracefulShutdownSignalHandler);
#if defined(SIGUSR1)
  std::signal(SIGUSR1, gracefulShutdownSignalHandler);
#endif
}

/**
 * \brief Returns whether a shutdown signal was received; polled by the drivers at cycle
 *        boundaries.
 */
bool requested() { return stopRequestedFlag.load(std::memory_order_relaxed); }

/**
 * \brief Reports the shutdown and exits; called after the driver has written its restart file.
 */
[[noreturn]] void exitAfterCheckpoint()
{
  std::cerr << "Termination signal received: binary restart file written, exiting\n";
  std::exit(0);
}
}  // namespace GracefulShutdown
