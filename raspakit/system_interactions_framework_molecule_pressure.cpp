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
import energy_factor;
import energy_status_inter;
import running_energy;
import units;
import threadpool;
import threading;

import <optional>;
import <iostream>;
import <algorithm>;
import <vector>;
import <span>;
import <cmath>;
import <thread>;
import <future>;
import <deque>;
import <semaphore>;


[[nodiscard]] std::pair<EnergyStatus, double3x3> System::computeFrameworkMoleculeEnergyStrainDerivative() noexcept
{
	double3 dr, posA, posB, f;
	double rr;

    double3x3 strainDerivative;
    EnergyStatus energy(components.size());

    const double cutOffVDWSquared = forceField.cutOffVDW * forceField.cutOffVDW;
    const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;
    const double prefactor = Units::CoulombicConversionFactor;

    std::span<const Atom> frameworkAtoms = spanOfFrameworkAtoms();
    std::span<const Atom> moleculeAtoms = spanOfMoleculeAtoms();
	if (moleculeAtoms.empty()) return std::make_pair(energy, strainDerivative);

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

				energy(compA, compB).VanDerWaals += 0.5 * energyFactor;
				energy(compB, compA).VanDerWaals += 0.5 * energyFactor;
			}
            if (!noCharges && rr < cutOffChargeSquared)
            {
                double r = std::sqrt(rr);
                double scaling = it1->scalingCoulomb * it2->scalingCoulomb;
                EnergyFactor energyFactor = prefactor * potentialCoulombEnergy(forceField, scaling, r, chargeA, chargeB);

                energy(compA, compB).CoulombicReal += 0.5 * energyFactor;
                energy(compB, compA).CoulombicReal += 0.5 * energyFactor;
            }
		}
	}

    return std::make_pair(energy, strainDerivative);
}