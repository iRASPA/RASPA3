module;

export module energy_interactions_intermolecular;

import energy_status;
import potential_energy_vdw;
//import potential_gradient_vdw;
import potential_energy_coulomb;
//import potential_gradient_coulomb;
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
import <execution>;

export EnergyStatus computeInterMolecularEnergy(size_t numberOfComponents, const ForceField &forceField, const SimulationBox &simulationBox, const std::span<const Atom> &moleculeAtoms) noexcept
{
	double3 dr, posA, posB, f;
	double rr;

    const double cutOffVDWSquared = forceField.cutOff * forceField.cutOff;
    const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;
    const double prefactor = Units::CoulombicConversionFactor;

    EnergyStatus energies(numberOfComponents);

	if (moleculeAtoms.empty()) return energies;

    for (std::span<const Atom>::iterator it1 = moleculeAtoms.begin(); it1 != moleculeAtoms.end() - 1; ++it1)
	{
		posA = it1->position;
		size_t molA = static_cast<size_t>(it1->moleculeId);
		size_t compA = static_cast<size_t>(it1->componentId);
		size_t typeA = static_cast<size_t>(it1->type);
		double scaleA = it1->scalingVDW;
		double chargeA = it1->charge;
		for (std::span<const Atom>::iterator it2 = it1 + 1; it2 != moleculeAtoms.end(); ++it2)
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
					EnergyFactor energyFactor = potentialVDWEnergy(forceField, scaling, rr, typeA, typeB);
					
					energies(compA, compB).VanDerWaals += 0.5 * energyFactor;
					energies(compB, compA).VanDerWaals += 0.5 * energyFactor;
                    energies.dUdlambda += energyFactor.dUdlambda;
				}
                // FIX: noCharges
                if (rr < cutOffChargeSquared)
                {
                    double r = std::sqrt(rr);
                    double scaling = it1->scalingCoulomb * it2->scalingCoulomb;
                    EnergyFactor energyFactor = prefactor * potentialCoulombEnergy(forceField, scaling, r, chargeA, chargeB);

                    energies(compA, compB).CoulombicReal += 0.5 * energyFactor;
                    energies(compB, compA).CoulombicReal += 0.5 * energyFactor;

                }
			}
		}
	}

    return energies;
}

/*
// https://livebook.manning.com/book/c-plus-plus-concurrency-in-action-second-edition/chapter-2/99
template<typename Iterator,typename T>
struct accumulate_block
{
    void operator()(Iterator first,Iterator last,T& result)
    {
        result=std::accumulate(first,last,result);
    }
};
template<typename Iterator,typename T>
T parallel_accumulate(Iterator first,Iterator last,T init)
{
    unsigned long const length=std::distance(first,last);
    if(!length)
        return init;
    unsigned long const min_per_thread=25;
    unsigned long const max_threads=
        (length+min_per_thread-1)/min_per_thread;
    unsigned long const hardware_threads=
        std::thread::hardware_concurrency();
    unsigned long const num_threads=
        std::min(hardware_threads!=0?hardware_threads:2,max_threads);
    unsigned long const block_size=length/num_threads;
    std::vector<T> results(num_threads);
    std::vector<std::thread>  threads(num_threads-1);
    Iterator block_start=first;
    for(unsigned long i=0;i<(num_threads-1);++i)
    {
        Iterator block_end=block_start;
        std::advance(block_end,block_size);
        threads[i]=std::thread(
            accumulate_block<Iterator,T>(),
            block_start,block_end,std::ref(results[i]));
        block_start=block_end;
    }
    accumulate_block<Iterator,T>()(
        block_start,last,results[num_threads-1]);

    for(auto& entry: threads)
           entry.join();
    return std::accumulate(results.begin(),results.end(),init);
}
*/
