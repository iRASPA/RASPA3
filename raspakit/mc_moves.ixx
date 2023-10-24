export module mc_moves;

import <vector>;
import <tuple>;
import <optional>;
import <span>;
import <fstream>;

import archive;
import component;
import atom;
import double3;
import double3x3;
import simulationbox;
import randomnumbers;
import system;
import energy_status;
import running_energy;


export struct MC_Moves
{
  MC_Moves() {}

  bool operator==(MC_Moves const&) const = default;

  uint64_t versionNumber{ 1 }; 
  double energyOverlapCriteria = 1e6;
  bool useDualCutOff{ false };

  void performRandomMove(RandomNumber &random, System& selectedSystem, System& selectedSecondystem, size_t selectedComponent, size_t &fractionalMoleculeSystem);
  void performRandomMoveProduction(RandomNumber &random, System& selectedSystem, System& selectedSecondSystem, size_t selectedComponent, size_t &fractionalMoleculeSystem, size_t currentBlock);

  std::optional<RunningEnergy> translationMove(RandomNumber &random, System& system, size_t selectedComponent, std::span<Atom> molecule) const;
  std::optional<RunningEnergy> randomTranslationMove(RandomNumber &random, System& system, size_t selectedComponent, std::span<Atom> molecule);

  std::optional<RunningEnergy> rotationMove(RandomNumber &random, System& system, size_t selectedComponent, std::span<Atom> molecule);
  std::optional<RunningEnergy> randomRotationMove(RandomNumber &random, System& system, size_t selectedComponent, std::span<Atom> molecule);

  std::optional<RunningEnergy> identityChangeMove(RandomNumber &random, System& system, size_t selectedComponent, size_t selectedMolecule, std::span<Atom> atoms);
  std::optional<RunningEnergy> reinsertionMove(RandomNumber &random, System& system, size_t selectedComponent, size_t selectedMolecule, std::span<Atom> atoms);

  std::pair<std::optional<RunningEnergy>, double3> insertionMove(RandomNumber &random, System& system, size_t selectedComponent);
  std::pair<std::optional<RunningEnergy>, double3> deletionMove(RandomNumber &random, System& system, size_t selectedComponent, size_t selectedMolecule);  
  std::pair<std::optional<RunningEnergy>, double3> swapMove_CFCMC(RandomNumber &random, System& system, size_t selectedComponent, size_t selectedMolecule, bool insertionDisabled = false, bool deletionDisabled = false);
  std::pair<std::optional<RunningEnergy>, double3> swapMove_CFCMC_CBMC(RandomNumber &random, System& system, size_t selectedComponent, size_t selectedMolecule, bool insertionDisabled = false, bool deletionDisabled = false);

  std::optional<double> WidomMove(RandomNumber &random, System& system, size_t selectedComponent);

  std::optional<std::pair<RunningEnergy, RunningEnergy>> GibbsSwapMove_CBMC(RandomNumber &random, System& systemA, System& systemB, size_t selectedComponent);
  std::optional<std::pair<RunningEnergy, RunningEnergy>> GibbsSwapMove_CFCMC(RandomNumber &random, System& systemA, System& systemB, size_t selectedComponent, size_t &fractionalMoleculeSystem);
  std::optional<std::pair<RunningEnergy, RunningEnergy>> GibbsSwapMove_CFCMC_CBMC(RandomNumber &random, System& systemA, System& systemB, size_t selectedComponent, size_t &fractionalMoleculeSystem);

  std::optional<RunningEnergy> volumeMove(RandomNumber &random, System &system) const;
  std::optional<std::pair<RunningEnergy, RunningEnergy>> GibbsVolumeMove(RandomNumber &random, System &systemA, System &systemB) const;

  std::optional<RunningEnergy> reactionMove(RandomNumber &random, System& system, const std::vector<size_t> reactantStoichiometry, const std::vector<size_t> productStoichiometry) const;
  std::optional<RunningEnergy> reactionMove_CFCMC(RandomNumber &random, System& system, const std::vector<size_t> reactantStoichiometry, const std::vector<size_t> productStoichiometry) const;
  std::optional<RunningEnergy> reactionMove_CFCMC_CBMC(RandomNumber &random, System& system, const std::vector<size_t> reactantStoichiometry, const std::vector<size_t> productStoichiometry) const;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const MC_Moves &mc);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, MC_Moves &mc);
};
