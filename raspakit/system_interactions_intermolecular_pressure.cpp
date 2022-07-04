module;

module system;

import energy_status;
import potential_energy_vdw;
import potential_gradient_vdw;
import potential_energy_coulomb;
import potential_gradient_coulomb;
import simulationbox;
import double3;
import double3x3;
import forcefield;
import atom;
import energy_factor;
import force_factor;
import energy_status_inter;
import units;

import <iostream>;
import <algorithm>;
import <vector>;
import <span>;
import <cmath>;
import <optional>;


std::pair<EnergyStatus, double3x3> System::computeInterMolecularMolecularPressure() noexcept
{
	double3 dr, posA, posB;
	double rr;

    const double cutOffVDWSquared = forceField.cutOff * forceField.cutOff;
    const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;
    const double prefactor = Units::CoulombicConversionFactor;

    EnergyStatus energy(components.size());
    double3x3 strainDerivativeTensor;

	std::span<Atom> moleculeAtoms = spanOfMoleculeAtoms();
	if (moleculeAtoms.empty()) return {energy, strainDerivativeTensor};

	for (std::span<Atom>::iterator it1 = moleculeAtoms.begin(); it1 != moleculeAtoms.end() - 1; ++it1)
	{
		posA = it1->position;
		size_t molA = static_cast<size_t>(it1->moleculeId);
		size_t compA = static_cast<size_t>(it1->componentId);
		size_t typeA = static_cast<size_t>(it1->type);
		double scaleA = it1->scalingVDW;
		double chargeA = it1->charge;
		for (std::span<Atom>::iterator it2 = it1 + 1; it2 != moleculeAtoms.end(); ++it2)
		{
			size_t molB = static_cast<size_t>(it2->moleculeId);
			size_t compB = static_cast<size_t>(it2->componentId);

			// skip interactions within the same molecule
			if (!((compA == compB) && (molA == molB)))
			{
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
					ForceFactor forceFactor = potentialVDWGradient(forceField, scaling, rr, typeA, typeB);
					
					energy(compA, compB).VanDerWaals += 0.5 * EnergyFactor(forceFactor.energy, 0.0);
					energy(compB, compA).VanDerWaals += 0.5 * EnergyFactor(forceFactor.energy, 0.0);

                    const double3 f = forceFactor.forceFactor * dr;

                    it1->gradient += f;
                    it2->gradient -= f;

                    strainDerivativeTensor.ax += f.x*dr.x;
                    strainDerivativeTensor.bx += f.y*dr.x;
                    strainDerivativeTensor.cx += f.z*dr.x;

                    strainDerivativeTensor.ay += f.x*dr.y;
                    strainDerivativeTensor.by += f.y*dr.y;
                    strainDerivativeTensor.cy += f.z*dr.y;

                    strainDerivativeTensor.az += f.x*dr.z;
                    strainDerivativeTensor.bz += f.y*dr.z;
                    strainDerivativeTensor.cz += f.z*dr.z;
				}
                if (!noCharges && rr < cutOffChargeSquared)
                {
                    double r = std::sqrt(rr);
                    double scaling = it1->scalingCoulomb * it2->scalingCoulomb;
                    ForceFactor energyFactor = prefactor * potentialCoulombGradient(forceField, scaling, r, chargeA, chargeB);

                    energy(compA, compB).CoulombicReal += 0.5 * EnergyFactor(energyFactor.energy, 0);
                    energy(compB, compA).CoulombicReal += 0.5 * EnergyFactor(energyFactor.energy, 0);

                    const double3 f = energyFactor.forceFactor * dr;

                    it1->gradient += f;
                    it2->gradient -= f;

                    strainDerivativeTensor.ax += f.x*dr.x;
                    strainDerivativeTensor.bx += f.y*dr.x;
                    strainDerivativeTensor.cx += f.z*dr.x;

                    strainDerivativeTensor.ay += f.x*dr.y;
                    strainDerivativeTensor.by += f.y*dr.y;
                    strainDerivativeTensor.cy += f.z*dr.y;

                    strainDerivativeTensor.az += f.x*dr.z;
                    strainDerivativeTensor.bz += f.y*dr.z;
                    strainDerivativeTensor.cz += f.z*dr.z;
                }
			}
		}
	}



    return {energy, strainDerivativeTensor};
}

