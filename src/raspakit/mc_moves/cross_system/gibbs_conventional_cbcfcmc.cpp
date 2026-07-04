module;

module mc_moves_gibbs_conventional_cbcfcmc;

import std;

import randomnumbers;
import running_energy;
import system;
import mc_moves_gibbs_conventional_common;

std::optional<std::pair<RunningEnergy, RunningEnergy>> MC_Moves::GibbsConventionalCBCFCMCMove(
    RandomNumber& random, System& systemA, System& systemB, std::size_t selectedComponent)
{
  return GibbsConventionalCommon::gibbsConventionalMove(random, systemA, systemB, selectedComponent, true);
}
