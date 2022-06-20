module;

module system;

import atom;
import energy_status;
import energy_status_intra;
import energy_status_inter;
import component;

import <iomanip>;
import <iostream>;
import <algorithm>;
import <vector>;
import <span>;
import <optional>;

void System::computeTotalEnergies() noexcept
{
  runningEnergies.zero();

  computeFrameworkMoleculeVDWEnergy();
  computeInterMolecularVDWEnergy();

  computeTailCorrectionVDWEnergy();

  computeEwaldFourierEnergy();

  for(size_t l = 0; l != components.size(); ++l)
  {
    if(components[l].type == Component::Type::Framework)
    {
       runningEnergies(l,l).VanDerWaals = 0.0;
       runningEnergies(l,l).VanDerWaalsTailCorrection = 0.0;
       runningEnergies(l,l).CoulombicReal = 0.0;
       runningEnergies(l,l).CoulombicFourier = 0.0;
    }
  }


  runningEnergies.sumTotal();
}

[[nodiscard]] const std::vector<std::pair<Atom, EnergyStatus>> System::computeExternalNonOverlappingEnergies(std::vector<Atom>& trialPositions) const noexcept
{
	std::vector<std::pair<Atom, EnergyStatus>> energies{};
	for (auto it = trialPositions.begin(); it != trialPositions.end(); ++it)
	{
		std::optional<EnergyStatus> interEnergy = computeInterMolecularVDWEnergy({ it,1 });
        if(!interEnergy.has_value()) continue;

        std::optional<EnergyStatus> frameworkEnergy = computeFrameworkMoleculeVDWEnergy({ it,1 });
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

		std::optional<EnergyStatus> interEnergy = computeInterMolecularVDWEnergy(span, skip);
        if(!interEnergy.has_value()) continue;

		std::optional<EnergyStatus> frameworkEnergy = computeFrameworkMoleculeVDWEnergy(span, skip);
        if(!frameworkEnergy.has_value()) continue;

		energies.push_back(std::make_pair(trialPositionSet, interEnergy.value() + frameworkEnergy.value()));
	}
	return energies;
}
