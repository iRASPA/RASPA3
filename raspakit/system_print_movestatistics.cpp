module;

module system;

import <iostream>;
import <random>;
import <sstream>;
import <cmath>;
import <optional>;

import double3;
import component;
import print;
import property_lambda_probability_histogram;
import simulationbox;
import units;
import loadings;
import averages;
import property_lambda_probability_histogram;
import property_widom;
import property_loading;
import move_statistics;
import mc_moves_probabilities_system;
import mc_moves_probabilities_particles;


inline std::string formatMoveStatistics(const std::string name, const MoveStatistics<double>& move)
{
  std::ostringstream stream;

  std::print(stream, "{} total:        {:10}\n", name, move.counts);
  std::print(stream, "{} constructed:  {:10}\n", name, move.constructed);
  std::print(stream, "{} accepted:     {:10}\n", name, move.accepted);
  std::print(stream, "{} fraction:     {:10f}\n", name, move.accepted / std::max(1.0, double(move.counts)));
  std::print(stream, "{} max-change:   {:10f}\n\n", name, move.maxChange);

  return stream.str();
}

inline std::string formatMoveStatistics(const std::string name, const MoveStatistics<double3> &move)
{
  std::ostringstream stream;

  std::print(stream, "{} total:        {:10} {:10} {:10}\n", name, move.counts.x, move.counts.y, move.counts.z);
  std::print(stream, "{} constructed:  {:10} {:10} {:10}\n", name, move.constructed.x, move.constructed.y, move.constructed.z);
  std::print(stream, "{} accepted:     {:10} {:10} {:10}\n", name, move.accepted.x, move.accepted.y, move.accepted.z);
  std::print(stream, "{} fraction:     {:10f} {:10f} {:10f}\n", name, move.accepted.x / std::max(1.0, double(move.counts.x)),
                 move.accepted.y / std::max(1.0, double(move.counts.y)), move.accepted.z / std::max(1.0, double(move.counts.z)));
  std::print(stream, "{} max-change:   {:10f} {:10f} {:10f}\n\n", name, move.maxChange.x, move.maxChange.y, move.maxChange.z);

  return stream.str();
}

std::string System::writeMCMoveStatistics() const
{
  std::ostringstream stream;

  if (mc_moves_probabilities.probabilityVolumeMove > 0.0) std::print(stream, formatMoveStatistics( "Volume", mc_moves_probabilities.statistics_VolumeMove));
  if (mc_moves_probabilities.probabilityGibbsVolumeMove > 0.0) std::print(stream, formatMoveStatistics("Gibbs Volume", mc_moves_probabilities.statistics_GibbsVolumeMove));

  for (size_t componentId = 0; const Component& component: components)
  {
    std::print(stream,"Component {} [{}]\n", componentId, component.name);

    std::print(stream, component.mc_moves_probabilities.writeMCMoveStatistics());

    if(component.hasFractionalMolecule)
    {
      double imposedChemicalPotential = std::log(beta * component.molFraction * pressure) / beta;
      double imposedFugacity = component.molFraction * pressure;

      std::print(stream, component.lambdaGC.writeAveragesStatistics(beta, imposedChemicalPotential, imposedFugacity));
      std::print(stream, component.lambdaGC.writeDUdLambdaStatistics(beta, imposedChemicalPotential, imposedFugacity));
    }

    if(component.mc_moves_probabilities.probabilityWidomMove > 0.0)
    {
      double imposedChemicalPotential = std::log(beta * component.molFraction * pressure) / beta;
      double imposedFugacity = component.molFraction * pressure;
      std::print(stream, component.averageRosenbluthWeights.writeAveragesStatistics(beta, imposedChemicalPotential, imposedFugacity));
    }

    ++componentId;
  }

 
  std::print(stream, "\n\n");

  return stream.str();
}
