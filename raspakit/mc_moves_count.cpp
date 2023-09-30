module mc_moves_count;

import <chrono>;
import <string>;
import <sstream>;
import <print>;

import double3;
import stringutils;


void MCMoveCount::clearCountStatistics()
{
  translationMove = size_t{0};
  randomTranslationMove = size_t{0};
  rotationMove = size_t{0};
  randomRotationMove = size_t{0};
  reinsertionMoveCBMC = size_t{0};
  swapInsertionMoveCBMC = size_t{0};
  swapDeletionMoveCBMC = size_t{0};
  swapLambdaMoveCFCMC = size_t{0};
  swapLambdaMoveCBCFCMC = size_t{0};
  GibbsSwapMoveCBMC = size_t{0};
  GibbsSwapLambdaMoveCFCMC = size_t{0};
  WidomMoveCBMC = size_t{0};
  WidomMoveCFCMC = size_t{0};
  WidomMoveCBCFCMC = size_t{0};
  volumeMove = size_t{0};
  GibbsVolumeMove = size_t{0};
}

const std::string MCMoveCount::writeMCMoveCountSystemStatistics(size_t countTotal) const
{
  std::ostringstream stream;

  std::print(stream, "\n");
  std::print(stream, "Volume:                          {:14f} [-]\n", 100.0 * static_cast<double>(volumeMove) / static_cast<double>(countTotal));
  std::print(stream, "Gibbs Volume:                    {:14f} [-]\n", 100.0 * static_cast<double>(GibbsVolumeMove) / static_cast<double>(countTotal));

  std::print(stream, "\n\n");

  return stream.str();
}

const std::string MCMoveCount::writeMCMoveCountComponentStatistics(size_t countTotal, size_t componentId, const std::string &componentName) const
{
  std::ostringstream stream;

  std::print(stream, "Component {} {}\n", componentId, componentName);

  std::print(stream, "\n");
  std::print(stream, "    Translation:                 {:14f} [%]\n", 100.0 * static_cast<double>(translationMove) / static_cast<double>(countTotal));
  std::print(stream, "    Random Translation:          {:14f} [%]\n", 100.0 * static_cast<double>(randomTranslationMove) / static_cast<double>(countTotal));
  std::print(stream, "    Rotation:                    {:14f} [%]\n", 100.0 * static_cast<double>(rotationMove) / static_cast<double>(countTotal));
  std::print(stream, "    Random Rotation:             {:14f} [%]\n", 100.0 * static_cast<double>(randomRotationMove) / static_cast<double>(countTotal));
  std::print(stream, "    Reinsertion CBMC:            {:14f} [%]\n", 100.0 * static_cast<double>(reinsertionMoveCBMC) / static_cast<double>(countTotal));
  std::print(stream, "    Insertion CBMC:              {:14f} [%]\n", 100.0 * static_cast<double>(swapInsertionMoveCBMC) / static_cast<double>(countTotal));
  std::print(stream, "    Deletion CBMC:               {:14f} [%]\n", 100.0 * static_cast<double>(swapDeletionMoveCBMC) / static_cast<double>(countTotal));
  std::print(stream, "    Insertion/deletion CFCMC:    {:14f} [%]\n", 100.0 * static_cast<double>(swapLambdaMoveCFCMC) / static_cast<double>(countTotal));
  std::print(stream, "    Insertion/deletion CB/CFCMC: {:14f} [%]\n", 100.0 * static_cast<double>(swapLambdaMoveCBCFCMC) / static_cast<double>(countTotal));
  std::print(stream, "    Gibbs-swap CBMC:             {:14f} [%]\n", 100.0 * static_cast<double>(GibbsSwapMoveCBMC) / static_cast<double>(countTotal));
  std::print(stream, "    Gibbs-swap CFCMC:            {:14f} [%]\n", 100.0 * static_cast<double>(GibbsSwapLambdaMoveCFCMC) / static_cast<double>(countTotal));
  std::print(stream, "    Widom CBMC:                  {:14f} [%]\n", 100.0 * static_cast<double>(WidomMoveCBMC) / static_cast<double>(countTotal));
  std::print(stream, "    Widom CFCMC:                 {:14f} [%]\n", 100.0 * static_cast<double>(WidomMoveCFCMC) / static_cast<double>(countTotal));
  std::print(stream, "    Widom CB/CFCMC:              {:14f} [%]\n", 100.0 * static_cast<double>(WidomMoveCBCFCMC) / static_cast<double>(countTotal));

  std::print(stream, "\n\n");

  return stream.str();
}

const std::string MCMoveCount::writeMCMoveCountAllSystemStatistics(size_t countTotal) const
{
  std::ostringstream stream;

  std::print(stream, "Translation:                 {:14f} [%]\n", 100.0 * static_cast<double>(translationMove) / static_cast<double>(countTotal));
  std::print(stream, "Random Translation:          {:14f} [%]\n", 100.0 * static_cast<double>(randomTranslationMove) / static_cast<double>(countTotal));
  std::print(stream, "Rotation:                    {:14f} [%]\n", 100.0 * static_cast<double>(rotationMove) / static_cast<double>(countTotal));
  std::print(stream, "Random rotation:             {:14f} [%]\n", 100.0 * static_cast<double>(randomRotationMove) / static_cast<double>(countTotal));
  std::print(stream, "Reinsertion CBMC:            {:14f} [%]\n", 100.0 * static_cast<double>(reinsertionMoveCBMC) / static_cast<double>(countTotal));
  std::print(stream, "Insertion CBMC:              {:14f} [%]\n", 100.0 * static_cast<double>(swapInsertionMoveCBMC) / static_cast<double>(countTotal));
  std::print(stream, "Deletion CBMC:               {:14f} [%]\n", 100.0 * static_cast<double>(swapDeletionMoveCBMC) / static_cast<double>(countTotal));
  std::print(stream, "Insertion/deletion CFCMC:    {:14f} [%]\n", 100.0 * static_cast<double>(swapLambdaMoveCFCMC) / static_cast<double>(countTotal));
  std::print(stream, "Insertion/deletion CB/CFCMC: {:14f} [%]\n", 100.0 * static_cast<double>(swapLambdaMoveCBCFCMC) / static_cast<double>(countTotal));
  std::print(stream, "Gibbs-swap CBMC:             {:14f} [%]\n", 100.0 * static_cast<double>(GibbsSwapMoveCBMC) / static_cast<double>(countTotal));
  std::print(stream, "Gibbs-swap CFCMC:            {:14f} [%]\n", 100.0 * static_cast<double>(GibbsSwapLambdaMoveCFCMC) / static_cast<double>(countTotal));
  std::print(stream, "Volume:                      {:14f} [%]\n", 100.0 * static_cast<double>(volumeMove) / static_cast<double>(countTotal));
  std::print(stream, "Widom CBMC:                  {:14f} [%]\n", 100.0 * static_cast<double>(WidomMoveCBMC) / static_cast<double>(countTotal));
  std::print(stream, "Widom CFCMC:                 {:14f} [%]\n", 100.0 * static_cast<double>(WidomMoveCFCMC) / static_cast<double>(countTotal));
  std::print(stream, "Widom CB/CFCMC:              {:14f} [%]\n", 100.0 * static_cast<double>(WidomMoveCBCFCMC) / static_cast<double>(countTotal));
  std::print(stream, "Gibbs Volume:                {:14f} [%]\n", 100.0 * static_cast<double>(GibbsVolumeMove) / static_cast<double>(countTotal));

  std::print(stream, "\n");
  std::print(stream, "Production mc-move count:    {:14d} [-]\n", countTotal);
  std::print(stream, "                All summed:  {:14d} [-]\n", total());
  std::print(stream, "                difference:  {:14d} [-]\n", countTotal - total());

  std::print(stream, "\n\n");

  return stream.str();
}
