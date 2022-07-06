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
import threadpool;

import <iostream>;
import <algorithm>;
import <vector>;
import <span>;
import <cmath>;
import <optional>;
import <thread>;
import <future>;

void System::computeInterMolecularEnergy() noexcept
{
	double3 dr, posA, posB, f;
	double rr;

    const double cutOffVDWSquared = forceField.cutOff * forceField.cutOff;
    const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;
    const double prefactor = Units::CoulombicConversionFactor;


	std::span<const Atom> moleculeAtoms = spanOfMoleculeAtoms();
	if (moleculeAtoms.empty()) return;

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
					
					runningEnergies(compA, compB).VanDerWaals += 0.5 * energyFactor;
					runningEnergies(compB, compA).VanDerWaals += 0.5 * energyFactor;
                    runningEnergies.dUdlambda += energyFactor.dUdlambda;
				}
                if (!noCharges && rr < cutOffChargeSquared)
                {
                    double r = std::sqrt(rr);
                    double scaling = it1->scalingCoulomb * it2->scalingCoulomb;
                    EnergyFactor energyFactor = prefactor * potentialCoulombEnergy(forceField, scaling, r, chargeA, chargeB);

                    runningEnergies(compA, compB).CoulombicReal += 0.5 * energyFactor;
                    runningEnergies(compB, compA).CoulombicReal += 0.5 * energyFactor;

                }
			}
		}
	}
}

[[nodiscard]] std::optional<EnergyStatus> System::computeInterMolecularSpanEnergy(std::span<const Atom>::iterator startIterator, std::span<const Atom>::iterator endIterator, std::span<Atom> atoms, std::make_signed_t<std::size_t> skip) const noexcept
{
	double3 dr, s, t;
	double rr;

	EnergyStatus energySum(components.size());

    const double overlapCriteria = forceField.overlapCriteria;
    const double cutOffVDWSquared = forceField.cutOff * forceField.cutOff;
    const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;
    const double prefactor = Units::CoulombicConversionFactor;

	for (std::span<const Atom>::iterator it1 = startIterator; it1 != endIterator; ++it1)
	{
    	size_t molA = static_cast<size_t>(it1->moleculeId);
		size_t compA = static_cast<size_t>(it1->componentId);
	    double3 posA = it1->position;
	    size_t typeA = static_cast<size_t>(it1->type);
	    double chargeA = it1->charge;
	    double scaleA = it1->scalingVDW;
	    double scalingCoulombA = it1->scalingCoulomb;

     	for (int index = 0; const Atom& atom : atoms)
	    {
		    if (index != skip)
		    {
			    double3 posB = atom.position;
			    size_t compB = static_cast<size_t>(atom.componentId);
			    size_t molB = static_cast<size_t>(atom.moleculeId);
			    size_t typeB = static_cast<size_t>(atom.type);
			    double scaleB = atom.scalingVDW;
			    double chargeB = atom.charge;
	            double scalingCoulombB = atom.scalingCoulomb;

				if (!(compA == compB && molA == molB))
				{
					dr = posA - posB;
					dr = simulationBox.applyPeriodicBoundaryConditions(dr);
					rr = double3::dot(dr, dr);

					if (rr < cutOffVDWSquared)
					{
						double scaling = scaleA * scaleB;
						EnergyFactor energyFactor = potentialVDWEnergy(forceField, scaling, rr, typeA, typeB);
                        if (energyFactor.energy > overlapCriteria) return std::nullopt;

						energySum(compA, compB).VanDerWaals += 0.5 * energyFactor;
						energySum(compB, compA).VanDerWaals += 0.5 * energyFactor;
                        energySum.dUdlambda += energyFactor.dUdlambda;
					}
                    if (!noCharges && rr < cutOffChargeSquared)
                    {
                        double r = std::sqrt(rr);
                        double scaling = scalingCoulombA * scalingCoulombB;
                        EnergyFactor energy = prefactor * potentialCoulombEnergy(forceField, scaling, r, chargeA, chargeB);

                        energySum(compA, compB).CoulombicReal += 0.5 * energy;
                        energySum(compB, compA).CoulombicReal += 0.5 * energy;
                    }
				}
			}
		    ++index;
		}
	}
	
	//energySum.sumTotal();
	return energySum;
}

[[nodiscard]] std::optional<EnergyStatus> System::computeInterMolecularEnergy(std::span<Atom> atoms, std::make_signed_t<std::size_t> skip) const noexcept
{
	double3 dr, s, t;
	double rr;

	EnergyStatus energySum(components.size());

    std::span<const Atom> moleculeAtoms = spanOfMoleculeAtoms();

    const double overlapCriteria = forceField.overlapCriteria;
    const double cutOffVDWSquared = forceField.cutOff * forceField.cutOff;
    const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;
    const double prefactor = Units::CoulombicConversionFactor;

	for (std::span<const Atom>::iterator it1 = moleculeAtoms.begin(); it1 != moleculeAtoms.end(); ++it1)
	{
    	size_t molA = static_cast<size_t>(it1->moleculeId);
		size_t compA = static_cast<size_t>(it1->componentId);
	    double3 posA = it1->position;
	    size_t typeA = static_cast<size_t>(it1->type);
	    double chargeA = it1->charge;
	    double scaleA = it1->scalingVDW;
	    double scalingCoulombA = it1->scalingCoulomb;

     	for (int index = 0; const Atom& atom : atoms)
	    {
		    if (index != skip)
		    {
			    double3 posB = atom.position;
			    size_t compB = static_cast<size_t>(atom.componentId);
			    size_t molB = static_cast<size_t>(atom.moleculeId);
			    size_t typeB = static_cast<size_t>(atom.type);
			    double scaleB = atom.scalingVDW;
			    double chargeB = atom.charge;
	            double scalingCoulombB = atom.scalingCoulomb;

				if (!(compA == compB && molA == molB))
				{
					dr = posA - posB;
					dr = simulationBox.applyPeriodicBoundaryConditions(dr);
					rr = double3::dot(dr, dr);

					if (rr < cutOffVDWSquared)
					{
						double scaling = scaleA * scaleB;
						EnergyFactor energyFactor = potentialVDWEnergy(forceField, scaling, rr, typeA, typeB);
                        if (energyFactor.energy > overlapCriteria) return std::nullopt;

						energySum(compA, compB).VanDerWaals += 0.5 * energyFactor;
						energySum(compB, compA).VanDerWaals += 0.5 * energyFactor;
                        energySum.dUdlambda += energyFactor.dUdlambda;
					}
                    if (!noCharges && rr < cutOffChargeSquared)
                    {
                        double r = std::sqrt(rr);
                        double scaling = scalingCoulombA * scalingCoulombB;
                        EnergyFactor energy = prefactor * potentialCoulombEnergy(forceField, scaling, r, chargeA, chargeB);

                        energySum(compA, compB).CoulombicReal += 0.5 * energy;
                        energySum(compB, compA).CoulombicReal += 0.5 * energy;
                    }
				}
			}
		    ++index;
		}
	}
	
	energySum.sumTotal();
	return energySum;
    
    //std::span<const Atom> moleculeAtoms = spanOfMoleculeAtoms();
    //return computeInterMolecularSpanEnergy(moleculeAtoms.begin(), moleculeAtoms.end(), atoms, skip);

    /*
  ThreadPool &pool = ThreadPool::instance();
  const size_t numberOfHelperThreads = pool.get_thread_count();

  std::vector<std::future<std::optional<EnergyStatus>>> threads(numberOfHelperThreads);

  std::span<const Atom> moleculeAtoms = spanOfMoleculeAtoms();
  size_t const block_size = moleculeAtoms.size() / (numberOfHelperThreads + 1);

  std::span<const Atom>::iterator block_start = moleculeAtoms.begin();
  for(size_t i = 0 ; i < numberOfHelperThreads; ++i)
  {
      std::span<const Atom>::iterator block_end = block_start;
      std::advance(block_end,block_size);

      threads[i] = pool.submit(std::bind(&System::computeInterMolecularSpanEnergy, this, block_start, block_end, atoms, skip));
      block_start=block_end;
  }
  std::optional<EnergyStatus> energy = computeInterMolecularSpanEnergy(block_start, moleculeAtoms.end(), atoms, skip);


  for(std::future<std::optional<EnergyStatus>> &future : threads)
  {
      future.wait();
  }

  if(!energy.has_value()) return std::nullopt;

  for(size_t i=0; i < threads.size(); ++i)
  {
    std::optional<EnergyStatus> energyContribution = threads[i].get();
    if(!energyContribution.has_value()) return std::nullopt;
    energy = energy.value() + energyContribution.value();
  }

  energy->sumTotal();

  return energy;
*/
}

[[nodiscard]] std::optional<EnergyStatus> System::computeInterMolecularEnergyDifference(std::span<const Atom> newatoms, std::span<const Atom> oldatoms) const noexcept
{
	double3 dr, s, t;
	double rr;

	EnergyStatus energySum(components.size());

    const double overlapCriteria = forceField.overlapCriteria;
    const double cutOffVDWSquared = forceField.cutOff * forceField.cutOff;
    const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;
    const double prefactor = Units::CoulombicConversionFactor;

    std::span<const Atom> moleculeAtoms = spanOfMoleculeAtoms();

	for (std::span<const Atom>::iterator it1 = moleculeAtoms.begin(); it1 != moleculeAtoms.end(); ++it1)
	{
		size_t molA = static_cast<size_t>(it1->moleculeId);
		size_t compA = static_cast<size_t>(it1->componentId);
		double3 posA = it1->position;
		size_t typeA = static_cast<size_t>(it1->type);
		double scaleA = it1->scalingVDW;
		double chargeA = it1->charge;
        double scalingCoulombA = it1->scalingCoulomb;

    	for (const Atom& atom : newatoms)
    	{
			size_t compB = static_cast<size_t>(atom.componentId);
		    size_t molB = static_cast<size_t>(atom.moleculeId);

			if (!(compA == compB && molA == molB))
			{
	     	     double3 posB = atom.position;
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

					energySum(compA, compB).VanDerWaals += 0.5 * energyFactor;
					energySum(compB, compA).VanDerWaals += 0.5 * energyFactor;
                    energySum.dUdlambda += energyFactor.dUdlambda;
				}
                if (!noCharges && rr < cutOffChargeSquared)
                {
                    double r = std::sqrt(rr);
                    double scaling = scalingCoulombA * scalingCoulombB;
                    EnergyFactor energy = prefactor * potentialCoulombEnergy(forceField, scaling, r, chargeA, chargeB);

                    energySum(compA, compB).CoulombicReal += 0.5 * energy;
                    energySum(compB, compA).CoulombicReal += 0.5 * energy;
                }
			}
		}

    	for (const Atom& atom : oldatoms)
    	{
			size_t compB = static_cast<size_t>(atom.componentId);
		    size_t molB = static_cast<size_t>(atom.moleculeId);

			if (!(compA == compB && molA == molB))
			{
	     	     double3 posB = atom.position;
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

					energySum(compA, compB).VanDerWaals -= 0.5 * energyFactor;
					energySum(compB, compA).VanDerWaals -= 0.5 * energyFactor;
                    energySum.dUdlambda -= energyFactor.dUdlambda;
				}
                if (!noCharges && rr < cutOffChargeSquared)
                {
                    double r = std::sqrt(rr);
                    double scaling = scalingCoulombA * scalingCoulombB;
                    EnergyFactor energy = prefactor * potentialCoulombEnergy(forceField, scaling, r, chargeA, chargeB);

                    energySum(compA, compB).CoulombicReal -= 0.5 * energy;
                    energySum(compB, compA).CoulombicReal -= 0.5 * energy;
                }
			}
		}
	}
	
	energySum.sumTotal();
	return energySum;
}
