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

void System::computeFrameworkMoleculeEnergy(const SimulationBox &box, std::span<const Atom> frameworkAtomPositions, std::span<const Atom> moleculeAtomPositions, RunningEnergy &energyStatus) noexcept
{
  double3 dr, posA, posB, f;
  double rr;

  const double cutOffVDWSquared = forceField.cutOff * forceField.cutOff;
  const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;
  const double prefactor = Units::CoulombicConversionFactor;

  if (moleculeAtomPositions.empty()) return;

  for (std::span<const Atom>::iterator it1 = frameworkAtomPositions.begin(); it1 != frameworkAtomPositions.end(); ++it1)
  {
    posA = it1->position;
    size_t typeA = static_cast<size_t>(it1->type);
    double scaleA = it1->scalingVDW;
    double chargeA = it1->charge;
    for (std::span<const Atom>::iterator it2 = moleculeAtomPositions.begin(); it2 != moleculeAtomPositions.end(); ++it2)
    {
      posB = it2->position;
      size_t typeB = static_cast<size_t>(it2->type);
      double scaleB = it2->scalingVDW;
      double chargeB = it2->charge;

      dr = posA - posB;
      dr = box.applyPeriodicBoundaryConditions(dr);
      rr = double3::dot(dr, dr);

      if (rr < cutOffVDWSquared)
      {
        double scaling = scaleA * scaleB;
        EnergyFactor energyFactor = potentialVDWEnergy(forceField, scaling, rr, typeA, typeB);

        energyStatus.frameworkMoleculeVDW += energyFactor.energy;
        energyStatus.dUdlambda += energyFactor.dUdlambda;
      }
      if (!noCharges && rr < cutOffChargeSquared)
      {
        double r = std::sqrt(rr);
        double scaling = it1->scalingCoulomb * it2->scalingCoulomb;
        EnergyFactor energyFactor = prefactor * potentialCoulombEnergy(forceField, scaling, r, chargeA, chargeB);

        energyStatus.frameworkMoleculeCharge += energyFactor.energy;
        energyStatus.dUdlambda += energyFactor.dUdlambda;
      }
    }
  }
}

[[nodiscard]] RunningEnergy System::computeFrameworkSpanMoleculeEnergy(std::span<const Atom>::iterator startIterator, std::span<const Atom>::iterator endIterator, std::span<Atom> atoms, std::make_signed_t<std::size_t> skip, std::atomic_flag & cancel) const
{
  RunningEnergy energySum;

  [[maybe_unused]] const double overlapCriteria = forceField.overlapCriteria;
  const double cutOffVDWSquared = forceField.cutOff * forceField.cutOff;
  const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;
  const double prefactor = Units::CoulombicConversionFactor;

  for (std::span<const Atom>::iterator it1 = startIterator; it1 != endIterator; ++it1)
  {
    if(cancel.test()) return energySum;
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
        size_t typeB = static_cast<size_t>(atom.type);
        double scaleB = atom.scalingVDW;
        double chargeB = atom.charge;
        double scalingCoulombB = atom.scalingCoulomb;

        double3 dr = posA - posB;
        dr = simulationBox.applyPeriodicBoundaryConditions(dr);
        double rr = double3::dot(dr, dr);

        if (rr < cutOffVDWSquared)
        {
          double scaling = scaleA * scaleB;
          EnergyFactor energyFactor = potentialVDWEnergy(forceField, scaling, rr, typeA, typeB);

          if (energyFactor.energy > overlapCriteria) 
          {
            cancel.test_and_set();
            return energySum;
          }
          energySum.frameworkMoleculeVDW += energyFactor.energy;
          energySum.dUdlambda += energyFactor.dUdlambda;
        }
        if (!noCharges && rr < cutOffChargeSquared)
        {
          double r = std::sqrt(rr);
          double scaling = scalingCoulombA * scalingCoulombB;
          EnergyFactor energyFactor = prefactor * potentialCoulombEnergy(forceField, scaling, r, chargeA, chargeB);

          energySum.frameworkMoleculeCharge += energyFactor.energy;
          energySum.dUdlambda += energyFactor.dUdlambda;
        }
      }
      ++index;
    }
  }

  return energySum;
}

[[nodiscard]] std::optional<RunningEnergy> System::computeFrameworkMoleculeEnergy(std::span<Atom> atoms, std::make_signed_t<std::size_t> skip) const noexcept
{
  ThreadPool &pool = ThreadPool::instance();

  std::span<const Atom> frameworkAtoms = spanOfFrameworkAtoms();

  [[maybe_unused]] const double overlapCriteria = forceField.overlapCriteria;
  const double cutOffVDWSquared = forceField.cutOff * forceField.cutOff;
  const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;
  const double prefactor = Units::CoulombicConversionFactor;

  std::atomic_flag cancel{false};

  switch(pool.threadingType)
  {
    case ThreadPool::ThreadingType::Serial:
    {
      RunningEnergy energySum;
      for (std::span<const Atom>::iterator it1 = frameworkAtoms.begin(); it1 != frameworkAtoms.end(); ++it1)
      {
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
            size_t typeB = static_cast<size_t>(atom.type);
            double scaleB = atom.scalingVDW;
            double chargeB = atom.charge;
            double scalingCoulombB = atom.scalingCoulomb;

            double3 dr = posA - posB;
            dr = simulationBox.applyPeriodicBoundaryConditions(dr);
            double rr = double3::dot(dr, dr);

            if (rr < cutOffVDWSquared)
            {
              double scaling = scaleA * scaleB;
              EnergyFactor energyFactor = potentialVDWEnergy(forceField, scaling, rr, typeA, typeB);
              if (energyFactor.energy > overlapCriteria)
              {
                return std::nullopt;
              }
              energySum.frameworkMoleculeVDW += energyFactor.energy;
              energySum.dUdlambda += energyFactor.dUdlambda;
            }
            if (!noCharges && rr < cutOffChargeSquared)
            {
              double r = std::sqrt(rr);
              double scaling = scalingCoulombA * scalingCoulombB;
              EnergyFactor energyFactor = prefactor * potentialCoulombEnergy(forceField, scaling, r, chargeA, chargeB);

              energySum.frameworkMoleculeCharge += energyFactor.energy;
              energySum.dUdlambda += energyFactor.dUdlambda;
            }
          }
          ++index;
        }
      }
      return energySum;
    }
    case ThreadPool::ThreadingType::OpenMP:
    {
      RunningEnergy energySum;
      #pragma omp parallel for reduction(+:energySum)
      for (std::span<const Atom>::iterator it1 = frameworkAtoms.begin(); it1 != frameworkAtoms.end(); ++it1)
      {
        if(!cancel.test())
        {
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
              size_t typeB = static_cast<size_t>(atom.type);
              double scaleB = atom.scalingVDW;
              double chargeB = atom.charge;
              double scalingCoulombB = atom.scalingCoulomb;

              double3 dr = posA - posB;
              dr = simulationBox.applyPeriodicBoundaryConditions(dr);
              double rr = double3::dot(dr, dr);

              if (rr < cutOffVDWSquared)
              {
                double scaling = scaleA * scaleB;
                EnergyFactor energyFactor = potentialVDWEnergy(forceField, scaling, rr, typeA, typeB);
                if (energyFactor.energy > overlapCriteria)
                {
                  cancel.test_and_set();
                }
                energySum.frameworkMoleculeVDW += energyFactor.energy;
                energySum.dUdlambda += energyFactor.dUdlambda;
              }
              if (!noCharges && rr < cutOffChargeSquared)
              {
                double r = std::sqrt(rr);
                double scaling = scalingCoulombA * scalingCoulombB;
                EnergyFactor energyFactor = prefactor * potentialCoulombEnergy(forceField, scaling, r, chargeA, chargeB);

                energySum.frameworkMoleculeCharge += energyFactor.energy;
                energySum.dUdlambda += energyFactor.dUdlambda;
              }
            }
            ++index;
          }
        }
      }
      if(cancel.test()) return std::nullopt;
      return energySum;
    }
    case ThreadPool::ThreadingType::ThreadPool:
    {
      const size_t numberOfHelperThreads = pool.getThreadCount();

      std::vector<std::future<RunningEnergy>> threads(numberOfHelperThreads);

      size_t const block_size = frameworkAtoms.size() / (numberOfHelperThreads + 1);


      std::span<const Atom>::iterator block_start = frameworkAtoms.begin();
      for(size_t i = 0 ; i != numberOfHelperThreads; ++i)
      {
        std::span<const Atom>::iterator block_end = block_start;
        std::advance(block_end,block_size);

        threads[i] = pool.enqueue(std::bind(&System::computeFrameworkSpanMoleculeEnergy, this, block_start, block_end, 
                                 atoms, skip, std::ref(cancel)));

        block_start=block_end;
      }
      RunningEnergy energy = computeFrameworkSpanMoleculeEnergy(block_start, frameworkAtoms.end(), atoms, skip, cancel);

      for(size_t i = 0; i != numberOfHelperThreads; ++i)
      {
        energy += threads[i].get();
      }
      if(cancel.test()) return std::nullopt;

      return energy;
    }
    case ThreadPool::ThreadingType::GPU_Offload:
    {
      RunningEnergy energySum;
      for (std::span<const Atom>::iterator it1 = frameworkAtoms.begin(); it1 != frameworkAtoms.end(); ++it1)
      {
        if(!cancel.test())
        {
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
              size_t typeB = static_cast<size_t>(atom.type);
              double scaleB = atom.scalingVDW;
              double chargeB = atom.charge;
              double scalingCoulombB = atom.scalingCoulomb;

              double3 dr = posA - posB;
              dr = simulationBox.applyPeriodicBoundaryConditions(dr);
              double rr = double3::dot(dr, dr);

              if (rr < cutOffVDWSquared)
              {
                double scaling = scaleA * scaleB;
                EnergyFactor energyFactor = potentialVDWEnergy(forceField, scaling, rr, typeA, typeB);
                if (energyFactor.energy > overlapCriteria)
                {
                  cancel.test_and_set();
                }
                energySum.frameworkMoleculeVDW += energyFactor.energy;
                energySum.dUdlambda += energyFactor.dUdlambda;
              }
              if (!noCharges && rr < cutOffChargeSquared)
              {
                double r = std::sqrt(rr);
                double scaling = scalingCoulombA * scalingCoulombB;
                EnergyFactor energyFactor = prefactor * potentialCoulombEnergy(forceField, scaling, r, chargeA, chargeB);

                energySum.frameworkMoleculeCharge += energyFactor.energy;
                energySum.dUdlambda += energyFactor.dUdlambda;
              }
            }
            ++index;
          }
        }
      }
      if(cancel.test()) return std::nullopt;
      return energySum;
    }
  }
}


[[nodiscard]] std::optional<RunningEnergy> System::computeFrameworkMoleculeEnergyDifference(std::span<const Atom> newatoms, std::span<const Atom> oldatoms) const noexcept
{
  double3 dr, s, t;
  double rr;

  RunningEnergy energySum{};

  const double overlapCriteria = forceField.overlapCriteria;
  const double cutOffVDWSquared = forceField.cutOff * forceField.cutOff;
  const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;
  const double prefactor = Units::CoulombicConversionFactor;

  std::span<const Atom> frameworkAtoms = spanOfFrameworkAtoms();

  for (std::span<const Atom>::iterator it1 = frameworkAtoms.begin(); it1 != frameworkAtoms.end(); ++it1)
  {
    double3 posA = it1->position;
    size_t typeA = static_cast<size_t>(it1->type);
    double scaleA = it1->scalingVDW;
    double chargeA = it1->charge;
    double scalingCoulombA = it1->scalingCoulomb;

    for (const Atom& atom : newatoms)
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

        energySum.frameworkMoleculeVDW +=  energyFactor.energy;
        energySum.dUdlambda += energyFactor.dUdlambda;
      }
      if (!noCharges && rr < cutOffChargeSquared)
      {
        double r = std::sqrt(rr);
        double scaling = scalingCoulombA * scalingCoulombB;
        EnergyFactor energyFactor = prefactor * potentialCoulombEnergy(forceField, scaling, r, chargeA, chargeB);

        energySum.frameworkMoleculeCharge += energyFactor.energy;
        energySum.dUdlambda += energyFactor.dUdlambda;
      }
    }

    for (const Atom& atom : oldatoms)
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

        energySum.frameworkMoleculeVDW -= energyFactor.energy;
        energySum.dUdlambda -= energyFactor.dUdlambda;
      }
      if (!noCharges && rr < cutOffChargeSquared)
      {
        double r = std::sqrt(rr);
        double scaling = scalingCoulombA * scalingCoulombB;
        EnergyFactor energyFactor = prefactor * potentialCoulombEnergy(forceField, scaling, r, chargeA, chargeB);

        energySum.frameworkMoleculeCharge -= energyFactor.energy;
        energySum.dUdlambda -= energyFactor.dUdlambda;
      }
    }
  }

  return energySum;
}

