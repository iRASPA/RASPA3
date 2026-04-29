module;

export module property_vacf;

import std;

import archive;
import double3;
import double4;
import atom;
import molecule;
import molecule;
import simulationbox;
import forcefield;
import component;
import velocity_autocorrelation_function_data;

// Computes Velocity Auto-Correlation Function (VASF)

export struct PropertyVelocityAutoCorrelationFunction
{
  PropertyVelocityAutoCorrelationFunction() {};

  PropertyVelocityAutoCorrelationFunction(std::size_t numberOfBuffersVACF,
                                          std::size_t bufferLengthVACF,
                                          std::size_t sampleEvery, 
                                          std::optional<std::size_t> writeEvery):
        numberOfBuffersVACF(numberOfBuffersVACF),
        bufferLengthVACF(bufferLengthVACF),
        sampleEvery(sampleEvery),
        writeEvery(writeEvery)
  {
  }

  PropertyVelocityAutoCorrelationFunction(const std::vector<std::size_t> &numberOfMoleculesPerComponent,
                                          std::size_t numberOfParticles,
                                          double timeStep,
                                          std::size_t numberOfBuffersVACF, 
                                          std::size_t bufferLengthVACF,
                                          std::size_t sampleEvery, 
                                          std::optional<std::size_t> writeEvery):
        numberOfMoleculesPerComponent(numberOfMoleculesPerComponent),
        numberOfComponents(numberOfMoleculesPerComponent.size()),
        numberOfParticles(numberOfParticles),
        timeStep(timeStep),
        numberOfBuffersVACF(numberOfBuffersVACF),
        bufferLengthVACF(bufferLengthVACF),
        sampleEvery(sampleEvery),
        writeEvery(writeEvery),
        originVACF(numberOfBuffersVACF, std::vector<double3>(numberOfParticles, double3(0.0, 0.0, 0.0))),
        originOnsagerVACF(numberOfBuffersVACF, std::vector<double3>(numberOfComponents, double3(0.0, 0.0, 0.0))),
        acfVACF(numberOfBuffersVACF,
                std::vector<std::vector<double4>>(numberOfComponents,
                                                  std::vector<double4>(bufferLengthVACF, double4(0.0, 0.0, 0.0, 0.0)))),
        acfOnsagerVACF(
            numberOfBuffersVACF,
            std::vector<std::vector<std::vector<double4>>>(
                numberOfComponents,
                std::vector<std::vector<double4>>(
                    numberOfComponents, std::vector<double4>(bufferLengthVACF, double4(0.0, 0.0, 0.0, 0.0))))),
        accumulatedAcfVACF(numberOfComponents, std::vector<double4>(bufferLengthVACF, double4(0.0, 0.0, 0.0, 0.0))),
        accumulatedAcfOnsagerVACF(
            numberOfBuffersVACF,
            std::vector<std::vector<double4>>(numberOfComponents,
                                              std::vector<double4>(bufferLengthVACF, double4(0.0, 0.0, 0.0, 0.0)))),
        countVACF(numberOfBuffersVACF),
        countAccumulatedVACF(0),
        sumVel(numberOfComponents)
  {
    // trick to space the origins evenly (see for example Rapaport 2004)
    for (std::size_t currentBuffer = 0; currentBuffer < numberOfBuffersVACF; ++currentBuffer)
    {
      countVACF[currentBuffer] = -static_cast<std::int64_t>(currentBuffer * bufferLengthVACF / numberOfBuffersVACF);
    }
  }

  std::uint64_t versionNumber{1};

  std::vector<std::size_t> numberOfMoleculesPerComponent;
  std::size_t numberOfComponents;
  std::size_t numberOfParticles;
  double timeStep;

  std::size_t numberOfBuffersVACF;
  std::size_t bufferLengthVACF;

  std::size_t sampleEvery;
  std::optional<std::size_t> writeEvery;

  std::vector<std::vector<double3>> originVACF;
  std::vector<std::vector<double3>> originOnsagerVACF;
  std::vector<std::vector<std::vector<double4>>> acfVACF;
  std::vector<std::vector<std::vector<std::vector<double4>>>> acfOnsagerVACF;

  std::vector<std::vector<double4>> accumulatedAcfVACF;
  std::vector<std::vector<std::vector<double4>>> accumulatedAcfOnsagerVACF;

  std::vector<std::int64_t> countVACF;
  std::size_t countAccumulatedVACF;
  std::vector<double3> sumVel;

  void addSample(std::size_t currentCycle, std::vector<Molecule> &molecules);

  std::vector<std::vector<VelocityAutoCorrelationFunctionData>> result();

  void writeOutput(std::size_t systemId, const std::vector<Component> &components, std::size_t currentCycle);

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive,
                                            const PropertyVelocityAutoCorrelationFunction &msd);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive,
                                            PropertyVelocityAutoCorrelationFunction &msd);
};
