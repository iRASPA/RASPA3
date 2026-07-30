module;

module apollonius_accessibility;

import std;

import double3;
import simulationbox;
import apollonius_network;
import pore_accessibility;

ApolloniusAccessibility ApolloniusAccessibility::create(const SimulationBox& simulationBox,
                                                        const std::vector<double3>& fractionalPositions,
                                                        const std::vector<double>& radii, double probeRadius)
{
  std::vector<double> inflatedRadii(radii.size());
  for (std::size_t i = 0; i < radii.size(); ++i) inflatedRadii[i] = radii[i] + probeRadius;

  ApolloniusAccessibility result;
  result.diagram = ApolloniusPoreNetwork::create(simulationBox, fractionalPositions, inflatedRadii);
  result.accessibility = PoreAccessibility::createFromNetwork(std::move(result.diagram.network),
                                                                 PoreAccessibility::Metric::Clearance, probeRadius);
  return result;
}
