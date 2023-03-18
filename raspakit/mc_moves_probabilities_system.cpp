module mc_moves_probabilities_system;

import <string>;
import <sstream>;

import double3;

import print;

import move_statistics;


void MCMoveProbabilitiesSystem::clear()
{
  statistics_VolumeMove.clear();
  statistics_GibbsVolumeMove.clear();
  statistics_GibbsSwapMove_CBMC.clear();
  statistics_GibbsSwapMove_CFCMC.clear();
  statistics_GibbsSwapMove_CFCMC_CBMC.clear();
}

void MCMoveProbabilitiesSystem::optimizeAcceptance()
{
  statistics_VolumeMove.optimizeAcceptance();
  statistics_GibbsVolumeMove.optimizeAcceptance();
}

void MCMoveProbabilitiesSystem::clearTimingStatistics()
{
  cpuTime_VolumeMove = std::chrono::duration<double>(0.0);
  cpuTime_GibbsVolumeMove = std::chrono::duration<double>(0.0);
  cpuTime_GibbsSwapMove_CBMC = std::chrono::duration<double>(0.0);
  cpuTime_GibbsSwapLambdaMove_CFCMC = std::chrono::duration<double>(0.0);
  cpuTime_GibbsSwapLambdaMove_CFCMC_CBMC = std::chrono::duration<double>(0.0);

  cpuTime_VolumeMove_NonEwald = std::chrono::duration<double>(0.0);
  cpuTime_GibbsVolumeMove_NonEwald = std::chrono::duration<double>(0.0);
  cpuTime_GibbsSwapMove_CBMC_NonEwald = std::chrono::duration<double>(0.0);
  cpuTime_GibbsSwapLambdaMove_CFCMC_NonEwald = std::chrono::duration<double>(0.0);
  cpuTime_GibbsSwapLambdaMove_CFCMC_CBMC_NonEwald = std::chrono::duration<double>(0.0);

  cpuTime_VolumeMove_Ewald = std::chrono::duration<double>(0.0);
  cpuTime_GibbsVolumeMove_Ewald = std::chrono::duration<double>(0.0);
  cpuTime_GibbsSwapMove_CBMC_Ewald = std::chrono::duration<double>(0.0);
  cpuTime_GibbsSwapLambdaMove_CFCMC_Ewald = std::chrono::duration<double>(0.0);
  cpuTime_GibbsSwapLambdaMove_CFCMC_CBMC_Ewald = std::chrono::duration<double>(0.0);
}

const std::string MCMoveProbabilitiesSystem::writeMCMoveCPUTimeStatistics() const
{
  std::ostringstream stream;
  if (cpuTime_VolumeMove.count() > 0.0)
  {
    std::print(stream, "    Volume move:            {:14f} [s]\n", cpuTime_VolumeMove.count());
    std::print(stream, "        Non-Ewald:              {:14f} [s]\n", cpuTime_VolumeMove_NonEwald.count());
    std::print(stream, "        Ewald:                  {:14f} [s]\n", cpuTime_VolumeMove_Ewald.count());
    std::print(stream, "        Overhead:               {:14f} [s]\n", cpuTime_VolumeMove.count() - cpuTime_VolumeMove_NonEwald.count() - cpuTime_VolumeMove_Ewald.count());
  }

  if (cpuTime_GibbsVolumeMove.count() > 0.0)
  {
    std::print(stream, "    Gibbs Volume:           {:14f} [s]\n", cpuTime_GibbsVolumeMove.count());
    std::print(stream, "        Non-Ewald:              {:14f} [s]\n", cpuTime_GibbsVolumeMove_NonEwald.count());
    std::print(stream, "        Ewald:                  {:14f} [s]\n", cpuTime_GibbsVolumeMove_Ewald.count());
    std::print(stream, "        Overhead:               {:14f} [s]\n", cpuTime_GibbsVolumeMove.count() - cpuTime_GibbsVolumeMove_NonEwald.count() - cpuTime_GibbsVolumeMove_Ewald.count());
  }
  if (cpuTime_GibbsSwapMove_CBMC.count() > 0.0)
  {
    std::print(stream, "    Gibbs swap (CBMC):      {:14f} [s]\n", cpuTime_GibbsSwapMove_CBMC.count());
    std::print(stream, "        Non-Ewald:              {:14f} [s]\n", cpuTime_GibbsSwapMove_CBMC_NonEwald.count());
    std::print(stream, "        Ewald:                  {:14f} [s]\n", cpuTime_GibbsSwapMove_CBMC_Ewald.count());
    std::print(stream, "        Overhead:               {:14f} [s]\n", cpuTime_GibbsSwapMove_CBMC.count() - cpuTime_GibbsSwapMove_CBMC_NonEwald.count() - cpuTime_GibbsSwapMove_CBMC_Ewald.count());
  }
  if (cpuTime_GibbsSwapLambdaMove_CFCMC.count() > 0.0)
  {
    std::print(stream, "    Gibbs swap (CFCMC):     {:14f} [s]\n", cpuTime_GibbsSwapLambdaMove_CFCMC.count());
    std::print(stream, "        Non-Ewald:              {:14f} [s]\n", cpuTime_GibbsSwapLambdaMove_CFCMC_NonEwald.count());
    std::print(stream, "        Ewald:                  {:14f} [s]\n", cpuTime_GibbsSwapLambdaMove_CFCMC_Ewald.count());
    std::print(stream, "        Overhead:               {:14f} [s]\n", cpuTime_GibbsSwapLambdaMove_CFCMC.count() - cpuTime_GibbsSwapLambdaMove_CFCMC_NonEwald.count() - cpuTime_GibbsSwapLambdaMove_CFCMC_Ewald.count());
  }
  if (cpuTime_GibbsSwapLambdaMove_CFCMC_CBMC.count() > 0.0)
  {
    std::print(stream, "    Gibbs swap (CB/CFCMC):  {:14f} [s]\n", cpuTime_GibbsSwapLambdaMove_CFCMC_CBMC.count());
    std::print(stream, "        Non-Ewald:              {:14f} [s]\n", cpuTime_GibbsSwapLambdaMove_CFCMC_CBMC_NonEwald.count());
    std::print(stream, "        Ewald:                  {:14f} [s]\n", cpuTime_GibbsSwapLambdaMove_CFCMC_CBMC_Ewald.count());
    std::print(stream, "        Overhead:               {:14f} [s]\n", cpuTime_GibbsSwapLambdaMove_CFCMC_CBMC.count() - cpuTime_GibbsSwapLambdaMove_CFCMC_CBMC_NonEwald.count() - cpuTime_GibbsSwapLambdaMove_CFCMC_CBMC_Ewald.count());
  }

  return stream.str();
}