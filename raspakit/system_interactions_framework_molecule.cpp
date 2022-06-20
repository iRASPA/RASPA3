module;

module system;

import energy_status;
import potential_energy_vdw;
import potential_energy_coulomb;
import simulationbox;
import double3;
import double3x3;
import forcefield;
import atom;
import energy_status_inter;
import units;

import <optional>;
import <iostream>;
import <algorithm>;
import <vector>;
import <span>;
import <cmath>;

void System::computeFrameworkMoleculeVDWEnergy() noexcept
{
	double3 dr, posA, posB, f;
	double rr, energy;

    const double cutOffVDWSquared = forceField.cutOff * forceField.cutOff;
    const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;
    const double prefactor = Units::CoulombicConversionFactor;

    std::span<const Atom> frameworkAtoms = spanOfFrameworkAtoms();
    std::span<const Atom> moleculeAtoms = spanOfMoleculeAtoms();
	if (moleculeAtoms.empty()) return;

	for (std::span<const Atom>::iterator it1 = frameworkAtoms.begin(); it1 != frameworkAtoms.end(); ++it1)
	{
		posA = it1->position;
		size_t compA = static_cast<size_t>(it1->componentId);
		size_t typeA = static_cast<size_t>(it1->type);
		double scaleA = it1->scalingVDW;
		double chargeA = it1->charge;
		for (std::span<const Atom>::iterator it2 = moleculeAtoms.begin(); it2 != moleculeAtoms.end(); ++it2)
		{
			size_t compB = static_cast<size_t>(it2->componentId);

			posB = it2->position;
			size_t typeB = static_cast<size_t>(it2->type);
			double scaleB = it2->scalingVDW;
			double chargeB = it2->charge;

			dr = posA - posB;
			dr = simulationBox.applyPeriodicBoundaryConditions(dr);
			rr = double3::dot(dr, dr);

			if (rr < cutOffVDWSquared)
			{
				double scaling = scaleA * scaleB;
				EnergyFactor energyFactor = potentialVDWEnergy(forceField, scaling, rr, typeA, typeB);

				runningEnergies(compA, compB).VanDerWaals += 0.5 * energyFactor.energy;
				runningEnergies(compB, compA).VanDerWaals += 0.5 * energyFactor.energy;
                runningEnergies.dUdlambda += energyFactor.dUdlambda;
			}
            if (rr < cutOffChargeSquared)
            {
                double r = std::sqrt(rr);
                double scaling = it1->scalingCoulomb * it2->scalingCoulomb;
                energy = prefactor * potentialCoulombEnergy(forceField, scaling, r, chargeA, chargeB);

                runningEnergies(compA, compB).CoulombicReal += 0.5 * energy;
                runningEnergies(compB, compA).CoulombicReal += 0.5 * energy;
            }
		}
	}
}

[[nodiscard]] std::optional<EnergyStatus> System::computeFrameworkMoleculeVDWEnergy(std::span<Atom> atoms, std::make_signed_t<std::size_t> skip) const noexcept
{
	double3 dr, s, t;
	double rr, energy;

	EnergyStatus energySum(components.size());

    [[maybe_unused]] const double overlapCriteria = forceField.overlapCriteria;
    const double cutOffVDWSquared = forceField.cutOff * forceField.cutOff;
    const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;
    const double prefactor = Units::CoulombicConversionFactor;

    std::span<const Atom> frameworkAtoms = spanOfFrameworkAtoms();

    for (std::span<const Atom>::iterator it1 = frameworkAtoms.begin(); it1 != frameworkAtoms.end(); ++it1)
    {
	   size_t compA = static_cast<size_t>(it1->componentId);
	   double3 posA = it1->position;
	   size_t typeA = static_cast<size_t>(it1->type);
	   double scaleA = it1->scalingVDW;
	   double chargeA = it1->charge;
       double scalingCoulombA = it1->scalingCoulomb;

	   for (int index = 0; const Atom& atom : atoms)
	   {
		    if (index != skip)
		    {
		    	double3 posB = atom.position;
		    	size_t compB = static_cast<size_t>(atom.componentId);
		    	size_t typeB = static_cast<size_t>(atom.type);
		    	double scaleB = atom.scalingVDW;
		    	double chargeB = atom.charge;
		    	double scalingCoulombB = atom.scalingCoulomb;

				dr = posA - posB;
				dr = simulationBox.applyPeriodicBoundaryConditions(dr);
				rr = double3::dot(dr, dr);

				if (rr < cutOffVDWSquared)
				{
					double scaling = scaleA * scaleB;
					EnergyFactor energyFactor = potentialVDWEnergy(forceField, scaling, rr, typeA, typeB);

                    if (energyFactor.energy > overlapCriteria) return std::nullopt;
					energySum(compA, compB).VanDerWaals += 0.5 * energyFactor.energy;
					energySum(compB, compA).VanDerWaals += 0.5 * energyFactor.energy;
                    energySum.dUdlambda += energyFactor.dUdlambda;
				}
                if (rr < cutOffChargeSquared)
                {
                    double r = std::sqrt(rr);
                    double scaling = scalingCoulombA * scalingCoulombB;
                    energy = prefactor * potentialCoulombEnergy(forceField, scaling, r, chargeA, chargeB);

                    energySum(compA, compB).CoulombicReal += 0.5 * energy;
                    energySum(compB, compA).CoulombicReal += 0.5 * energy;
                }
			}
		    ++index;
		}
	}

	energySum.sumTotal();
	return energySum;
}


[[nodiscard]] std::optional<EnergyStatus> System::computeFrameworkMoleculeEnergyDifference(std::span<const Atom> newatoms, std::span<const Atom> oldatoms) const noexcept
{
	double3 dr, s, t;
	double rr, energy;

	EnergyStatus energySum(components.size());

    const double overlapCriteria = forceField.overlapCriteria;
    const double cutOffVDWSquared = forceField.cutOff * forceField.cutOff;
    const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;
    const double prefactor = Units::CoulombicConversionFactor;

    std::span<const Atom> frameworkAtoms = spanOfFrameworkAtoms();

    for (std::span<const Atom>::iterator it1 = frameworkAtoms.begin(); it1 != frameworkAtoms.end(); ++it1)
    {
    	size_t compA = static_cast<size_t>(it1->componentId);
    	double3 posA = it1->position;
    	size_t typeA = static_cast<size_t>(it1->type);
    	double scaleA = it1->scalingVDW;
    	double chargeA = it1->charge;
        double scalingCoulombA = it1->scalingCoulomb;

    	for (const Atom& atom : newatoms)
     	{
      	  double3 posB = atom.position;
      	  size_t compB = static_cast<size_t>(atom.componentId);
      	  size_t typeB = static_cast<size_t>(atom.type);
      	  double scaleB = atom.scalingVDW;
      	  double chargeB = atom.charge;
          double scalingCoulombB = atom.scalingCoulomb;

		  dr = posA - posB;
		  dr = simulationBox.applyPeriodicBoundaryConditions(dr);
		  rr = double3::dot(dr, dr);

		  if (rr < cutOffVDWSquared)
		  {
		  	double scaling = scaleA * scaleB;
		  	EnergyFactor energyFactor = potentialVDWEnergy(forceField, scaling, rr, typeA, typeB);
            if (energyFactor.energy > overlapCriteria) return std::nullopt;

		  	energySum(compA, compB).VanDerWaals += 0.5 * energyFactor.energy;
		  	energySum(compB, compA).VanDerWaals += 0.5 * energyFactor.energy;
            energySum.dUdlambda += energyFactor.dUdlambda;
		  }
          if (rr < cutOffChargeSquared)
          {
              double r = std::sqrt(rr);
              double scaling = scalingCoulombA * scalingCoulombB;
              energy = prefactor * potentialCoulombEnergy(forceField, scaling, r, chargeA, chargeB);

              energySum(compA, compB).CoulombicReal += 0.5 * energy;
              energySum(compB, compA).CoulombicReal += 0.5 * energy;
          }
		}

    	for (const Atom& atom : oldatoms)
     	{
      	  double3 posB = atom.position;
      	  size_t compB = static_cast<size_t>(atom.componentId);
      	  size_t typeB = static_cast<size_t>(atom.type);
      	  double scaleB = atom.scalingVDW;
      	  double chargeB = atom.charge;
          double scalingCoulombB = atom.scalingCoulomb;

		  dr = posA - posB;
		  dr = simulationBox.applyPeriodicBoundaryConditions(dr);
		  rr = double3::dot(dr, dr);

		  if (rr < cutOffVDWSquared)
		  {
		  	double scaling = scaleA * scaleB;
		  	EnergyFactor energyFactor = potentialVDWEnergy(forceField, scaling, rr, typeA, typeB);

                 //if (energyFactor.energy > overlapCriteria) return std::nullopt;
		  	energySum(compA, compB).VanDerWaals -= 0.5 * energyFactor.energy;
		  	energySum(compB, compA).VanDerWaals -= 0.5 * energyFactor.energy;
            energySum.dUdlambda -= energyFactor.dUdlambda;
		  }
          if (rr < cutOffChargeSquared)
          {
              double r = std::sqrt(rr);
              double scaling = scalingCoulombA * scalingCoulombB;
              energy = prefactor * potentialCoulombEnergy(forceField, scaling, r, chargeA, chargeB);

              energySum(compA, compB).CoulombicReal -= 0.5 * energy;
              energySum(compB, compA).CoulombicReal -= 0.5 * energy;
          }
		}
	}

	energySum.sumTotal();
	return energySum;
}


