module;

module system;

import atom;
import energy_factor;
import energy_status;
import energy_status_intra;
import energy_status_inter;
import component;
import double3;
import double3x3;
import forcefield;
import simulationbox;

import <iomanip>;
import <iostream>;
import <algorithm>;
import <vector>;
import <span>;
import <optional>;

void System::computeTotalEnergies() noexcept
{
  runningEnergies.zero();

  computeFrameworkMoleculeEnergy();
  computeInterMolecularEnergy();

  computeTailCorrectionVDWEnergy();

  computeEwaldFourierEnergy();

  for(size_t l = 0; l != components.size(); ++l)
  {
    if(components[l].type == Component::Type::Framework)
    {
       runningEnergies(l,l).VanDerWaals = EnergyFactor(0.0, 0.0);
       runningEnergies(l,l).VanDerWaalsTailCorrection = EnergyFactor(0.0, 0.0);
       runningEnergies(l,l).CoulombicReal = EnergyFactor(0.0, 0.0);
       runningEnergies(l,l).CoulombicFourier = EnergyFactor(0.0, 0.0);
    }
  }

  runningEnergies.sumTotal();
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
    std::pair<EnergyStatus, double3x3> pressureInfo = computeInterMolecularMolecularPressure();

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
    pressureInfo.second -= correctionTerm;

    return pressureInfo;
}

[[nodiscard]] const std::vector<std::pair<Atom, EnergyStatus>> System::computeExternalNonOverlappingEnergies(std::vector<Atom>& trialPositions) const noexcept
{
	std::vector<std::pair<Atom, EnergyStatus>> energies{};
	for (auto it = trialPositions.begin(); it != trialPositions.end(); ++it)
	{
        // skip trial-positions that have an overlap in inter-molecular energy
		std::optional<EnergyStatus> interEnergy = computeInterMolecularEnergy({ it,1 });
        if(!interEnergy.has_value()) continue;

        // skip trial-positions that have an overlap in framework-molecule energy
        std::optional<EnergyStatus> frameworkEnergy = computeFrameworkMoleculeEnergy({ it,1 });
        if(!frameworkEnergy.has_value()) continue;

	    energies.push_back(std::make_pair(*it, interEnergy.value() + frameworkEnergy.value()));
	}
	return energies;
}

const std::vector<std::pair<std::vector<Atom>,EnergyStatus>> System::computeExternalNonOverlappingEnergies(std::vector<std::vector<Atom>>& trialPositionSets, std::make_signed_t<std::size_t> skip) const noexcept
{
	std::vector<std::pair<std::vector<Atom>,EnergyStatus>> energies{};

	for (std::vector<Atom> trialPositionSet : trialPositionSets)
	{
        std::span<Atom> span = std::span<Atom>(trialPositionSet.begin(), trialPositionSet.end());

		std::optional<EnergyStatus> interEnergy = computeInterMolecularEnergy(span, skip);
        if(!interEnergy.has_value()) continue;

		std::optional<EnergyStatus> frameworkEnergy = computeFrameworkMoleculeEnergy(span, skip);
        if(!frameworkEnergy.has_value()) continue;

		energies.push_back(std::make_pair(trialPositionSet, interEnergy.value() + frameworkEnergy.value()));
	}
	return energies;
}
