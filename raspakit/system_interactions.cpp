module;

module system;

import atom;
import energy_factor;
import energy_status;
import energy_status_intra;
import energy_status_inter;
import running_energy;
import component;
import double3;
import double3x3;
import forcefield;
import simulationbox;
import units;
import potential_energy_vdw;
import potential_energy_coulomb;
import interactions_framework_molecule;
import interactions_intermolecular;

import <iomanip>;
import <iostream>;
import <algorithm>;
import <vector>;
import <span>;
import <optional>;
import <cmath>;

// system_interactions.cpp

inline std::pair<EnergyStatus, double3x3> pair_acc(const std::pair<EnergyStatus, double3x3> &lhs, const std::pair<EnergyStatus, double3x3> &rhs)
{
  return std::make_pair(lhs.first + rhs.first, lhs.second + rhs.second);
}

void System::precomputeTotalRigidEnergy() noexcept
{
  rigidEnergies.zero();
  computeEwaldFourierRigidEnergy(simulationBox, rigidEnergies);
}

void System::recomputeTotalEnergies() noexcept
{
  runningEnergies.zero();

  std::span<const Atom> frameworkAtomPositions = spanOfFrameworkAtoms();
  std::span<const Atom> moleculeAtomPositions = spanOfMoleculeAtoms();

  Interactions::computeFrameworkMoleculeEnergy(forceField, simulationBox, frameworkAtomPositions, moleculeAtomPositions, runningEnergies);
  Interactions::computeInterMolecularEnergy(forceField, simulationBox, moleculeAtomPositions, runningEnergies);

  Interactions::computeFrameworkMoleculeTailEnergy(forceField, simulationBox, frameworkAtomPositions, moleculeAtomPositions, runningEnergies);
  Interactions::computeInterMolecularTailEnergy(forceField, simulationBox,moleculeAtomPositions, runningEnergies);

  computeEwaldFourierEnergy(simulationBox, runningEnergies);

  // correct for the energy of rigid parts
  runningEnergies -= rigidEnergies;
}


RunningEnergy System::computeTotalEnergies() noexcept
{
  RunningEnergy runningEnergy{};

  std::span<const Atom> frameworkAtomPositions = spanOfFrameworkAtoms();
  std::span<const Atom> moleculeAtomPositions = spanOfMoleculeAtoms();

  Interactions::computeFrameworkMoleculeEnergy(forceField, simulationBox, frameworkAtomPositions, moleculeAtomPositions, runningEnergy);
  Interactions::computeInterMolecularEnergy(forceField, simulationBox, moleculeAtomPositions, runningEnergy);

  Interactions::computeFrameworkMoleculeTailEnergy(forceField, simulationBox, frameworkAtomPositions, moleculeAtomPositions, runningEnergy);
  Interactions::computeInterMolecularTailEnergy(forceField, simulationBox, moleculeAtomPositions, runningEnergy);

  computeEwaldFourierEnergy(simulationBox, runningEnergy);

  // correct for the energy of rigid parts
  runningEnergy -= rigidEnergies;

  return runningEnergy;
}

void System::computeTotalGradients() noexcept
{
}

std::pair<EnergyStatus, double3x3> System::computeMolecularPressure() noexcept
{
  for(Atom& atom: atomPositions)
  {
    atom.gradient = double3(0.0, 0.0, 0.0);
  }
  std::pair<EnergyStatus, double3x3> pressureInfo = computeFrameworkMoleculeEnergyStrainDerivative();

  pressureInfo = pair_acc(pressureInfo, computeInterMolecularEnergyStrainDerivative());
  pressureInfo = pair_acc(pressureInfo, computeEwaldFourierEnergyStrainDerivative());

  pressureInfo.first.sumTotal();

  // Correct rigid molecule contribution using the constraints forces
  double3x3 correctionTerm;
  for(size_t componentId = 0; componentId < components.size(); ++componentId)
  {
    if(components[componentId].rigid)
    {
      for(size_t i = 0; i < numberOfMoleculesPerComponent[componentId]; ++i)
      {
        std::span<Atom> span = spanOfMolecule(componentId, i);

        double totalMass = 0.0;
        double3 com(0.0, 0.0, 0.0);
        for(const Atom& atom: span)
        {
          double mass = forceField.pseudoAtoms[static_cast<size_t>(atom.type)].mass;
          com += mass * atom.position;
          totalMass += mass;
        }
        com = com / totalMass;

        for(const Atom& atom: span)
        {
          correctionTerm.ax += (atom.position.x - com.x) * atom.gradient.x;
          correctionTerm.ay += (atom.position.x - com.x) * atom.gradient.y;
          correctionTerm.az += (atom.position.x - com.x) * atom.gradient.z;

          correctionTerm.bx += (atom.position.y - com.y) * atom.gradient.x;
          correctionTerm.by += (atom.position.y - com.y) * atom.gradient.y;
          correctionTerm.bz += (atom.position.y - com.y) * atom.gradient.z;

          correctionTerm.cx += (atom.position.z - com.z) * atom.gradient.x;
          correctionTerm.cy += (atom.position.z - com.z) * atom.gradient.y;
          correctionTerm.cz += (atom.position.z - com.z) * atom.gradient.z;
        }
      }
    }
  }
  pressureInfo.second = -(pressureInfo.second - correctionTerm);

  return pressureInfo;
}

