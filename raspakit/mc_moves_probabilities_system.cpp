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
}

void MCMoveProbabilitiesSystem::optimizeAcceptance()
{
  statistics_VolumeMove.optimizeAcceptance(0.01, 1.5);
  statistics_GibbsVolumeMove.optimizeAcceptance(0.01, 1.5);
}
