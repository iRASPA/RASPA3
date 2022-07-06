module;

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
import energy_status_intra;
import units;
import threadpool;

import <optional>;
import <iostream>;
import <algorithm>;
import <vector>;
import <span>;
import <cmath>;
import <tuple>;
import <thread>;
import <future>;

module energy_interactions_framework_molecule;

inline EnergyStatus computeFrameworkMoleculeEnergyParallel(const size_t numberOfComponents, const ForceField &forceField, const SimulationBox &simulationBox, std::span<Atom> atoms, std::make_signed_t<std::size_t> skip, std::span<const Atom>::iterator startIterator, std::span<const Atom>::iterator endIterator)
{
    double3 dr, s, t;
    double rr;

    [[maybe_unused]] const double overlapCriteria = forceField.overlapCriteria;
    const double cutOffVDWSquared = forceField.cutOff * forceField.cutOff;
    const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;
    const double prefactor = Units::CoulombicConversionFactor;

    EnergyStatus energySum(numberOfComponents);

    for (std::span<const Atom>::iterator it1 = startIterator; it1 != endIterator; ++it1)
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

                    //if (energyFactor.energy > overlapCriteria) 
                    //{
                    //    //energySum.second = true;
                    //    return;
                    //}
                    energySum(compA, compB).VanDerWaals += 0.5 * energyFactor;
                    energySum(compB, compA).VanDerWaals += 0.5 * energyFactor;
                    energySum.dUdlambda += energyFactor.dUdlambda;
                }
                //if (!noCharges && rr < cutOffChargeSquared)
                if ( rr < cutOffChargeSquared)
                {
                    double r = std::sqrt(rr);
                    double scaling = scalingCoulombA * scalingCoulombB;
                    EnergyFactor energy = prefactor * potentialCoulombEnergy(forceField, scaling, r, chargeA, chargeB);

                    energySum(compA, compB).CoulombicReal += 0.5 * energy;
                    energySum(compB, compA).CoulombicReal += 0.5 * energy;
                }
            }
            ++index;
        }
    }

    return energySum;
}

std::optional<EnergyStatus> computeFrameworkMoleculeEnergyPar([[maybe_unused]] size_t numberOfComponents, [[maybe_unused]] const ForceField &forceField, [[maybe_unused]] const SimulationBox &simulationBox, [[maybe_unused]] std::span<Atom> &atoms, [[maybe_unused]] std::make_signed_t<std::size_t> skip, [[maybe_unused]] std::span<const Atom> &frameworkAtoms) noexcept
{
    ThreadPool &pool = ThreadPool::instance();
    const size_t numberOfHelperThreads = pool.get_thread_count();

    std::vector<std::future<EnergyStatus>> threads(numberOfHelperThreads);

    size_t const block_size = frameworkAtoms.size() / (numberOfHelperThreads + 1);

    std::span<const Atom>::iterator block_start=frameworkAtoms.begin();
    for(size_t i = 0 ; i < numberOfHelperThreads; ++i)
    {
        std::span<const Atom>::iterator block_end=block_start;
        std::advance(block_end,block_size);

        threads[i] = pool.submit(computeFrameworkMoleculeEnergyParallel,
            numberOfComponents, forceField, simulationBox, atoms, skip, block_start, block_end); 
        block_start=block_end;
    }
    EnergyStatus energy = computeFrameworkMoleculeEnergyParallel(numberOfComponents, forceField, simulationBox, atoms, skip, block_start, frameworkAtoms.end());

    for(size_t i=0; i < threads.size(); ++i)
    {
      energy += threads[i].get();
    }

    energy.sumTotal();

    return energy;
}

