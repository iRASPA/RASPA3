module;

module system;

import <vector>;
import <numbers>;
import <iostream>;

import forcefield;
import energy_factor;
import potential_correction_vdw;
import energy_status;
import energy_status_inter;
import component;
import atom;
import simulationbox;
import units;

void System::computeTailCorrectionVDWEnergy() noexcept
{
  for(size_t k = 0; k != components.size(); ++k)
  {
    for(size_t l = 0; l != components.size(); ++l)
    {
      double energy = 0.0;
      //double pressure = 0.0;
      for(size_t i = 0; i != forceField.pseudoAtoms.size(); ++i)
      {
        for(size_t j = 0; j != forceField.pseudoAtoms.size(); ++j)
        {
          if(forceField(i,j).tailCorrection)
          {
            energy += 2.0 * std::numbers::pi * static_cast<double>(numberOfPseudoAtoms[k][i]) * 
                                               static_cast<double>(numberOfPseudoAtoms[l][j]) * 
                                               potentialCorrectionVDW(forceField,i,j);
            //pressure-=(2.0/3.0)*M_PI*(NumberOfPseudoAtomsType[CurrentSystem][i]-NumberOfFractionalPseudoAtomsType[CurrentSystem][i])*
            //                         (NumberOfPseudoAtomsType[CurrentSystem][j]-NumberOfFractionalPseudoAtomsType[CurrentSystem][j])*
            //                         PotentialCorrectionPressure(i,j,CutOffVDW);
          }
          else if(!forceField(i,j).shiftPotential)
          {
            // impulsive correction
            //pressure+=(2.0/3.0)*M_PI*(NumberOfPseudoAtomsType[CurrentSystem][i]-NumberOfFractionalPseudoAtomsType[CurrentSystem][i])*
            //                         (NumberOfPseudoAtomsType[CurrentSystem][j]-NumberOfFractionalPseudoAtomsType[CurrentSystem][j])*
            //                         CUBE(CutOffVDW)*PotentialValue(i,j,CutOffVDWSquared,1.0);
          }
        }
      }
      runningEnergies(k, l).VanDerWaalsTailCorrection = EnergyFactor(energy / simulationBox.volume, 0.0);
    }
  }
}

[[nodiscard]] EnergyStatus System::computeTailCorrectionVDWOldEnergy() const noexcept
{
  EnergyStatus energySum(components.size());

  for(size_t k = 0; k != components.size(); ++k)
  {
    for(size_t l = 0; l != components.size(); ++l)
    {
      double energy = 0.0;
      for(size_t i = 0; i != forceField.pseudoAtoms.size(); ++i)
      {
        for(size_t j = 0; j != forceField.pseudoAtoms.size(); ++j)
        {
          if(forceField(i,j).tailCorrection)
          {
            energy += 2.0 * std::numbers::pi * static_cast<double>(numberOfPseudoAtoms[k][i]) * 
                                               static_cast<double>(numberOfPseudoAtoms[l][j]) * 
                                               potentialCorrectionVDW(forceField,i,j);
          }
        }
      }
      energySum(k, l).VanDerWaalsTailCorrection = EnergyFactor(energy / simulationBox.volume, 0.0);
    }
  }
  energySum.sumTotal();
  return energySum;
}

[[nodiscard]] EnergyStatus System::computeTailCorrectionVDWAddEnergy(size_t selectedComponent) const noexcept
{
  EnergyStatus energySum(components.size());

  std::vector<std::vector<size_t>> newNumberOfPseudoAtoms(numberOfPseudoAtoms);

  for(const Atom& atom: components[selectedComponent].atoms)
  {
    size_t type = static_cast<size_t>(atom.type);
    newNumberOfPseudoAtoms[selectedComponent][type] += 1;
  }

  for(size_t k = 0; k != components.size(); ++k)
  {
    for(size_t l = 0; l != components.size(); ++l)
    {
      double energy = 0.0;
      for(size_t i = 0; i != forceField.pseudoAtoms.size(); ++i)
      {
        for(size_t j = 0; j != forceField.pseudoAtoms.size(); ++j)
        {
          if(forceField(i,j).tailCorrection)
          {
            energy += 2.0 * std::numbers::pi * static_cast<double>(newNumberOfPseudoAtoms[k][i]) * 
                                               static_cast<double>(newNumberOfPseudoAtoms[l][j]) * 
                                               potentialCorrectionVDW(forceField,i,j);
          }
        }
      }
      energySum(k, l).VanDerWaalsTailCorrection = EnergyFactor(energy / simulationBox.volume, 0.0);
    }
  }
  energySum.sumTotal();

  return energySum;
}

[[nodiscard]] EnergyStatus System::computeTailCorrectionVDWRemoveEnergy(size_t selectedComponent) const noexcept
{
  EnergyStatus energySum(components.size());

  std::vector<std::vector<size_t>> newNumberOfPseudoAtoms(numberOfPseudoAtoms);

  for(const Atom& atom: components[selectedComponent].atoms)
  {
    size_t type = static_cast<size_t>(atom.type);
    newNumberOfPseudoAtoms[selectedComponent][type] -= 1;
  }

  for(size_t k = 0; k != components.size(); ++k)
  {
    for(size_t l = 0; l != components.size(); ++l)
    {
      double energy = 0.0;
      for(size_t i = 0; i != forceField.pseudoAtoms.size(); ++i)
      {
        for(size_t j = 0; j != forceField.pseudoAtoms.size(); ++j)
        {
          if(forceField(i,j).tailCorrection)
          {
            energy += 2.0 * std::numbers::pi * static_cast<double>(newNumberOfPseudoAtoms[k][i]) * 
                                               static_cast<double>(newNumberOfPseudoAtoms[l][j]) * 
                                               potentialCorrectionVDW(forceField,i,j);
          }
        }
      }
      energySum(k, l).VanDerWaalsTailCorrection = EnergyFactor(energy / simulationBox.volume, 0.0);
    }
  }
  energySum.sumTotal();

  return energySum;
}
