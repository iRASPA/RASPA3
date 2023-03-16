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


std::pair<EnergyStatus, double3x3> System::computeInterMolecularEnergyStrainDerivative() noexcept
{
  double3 dr, posA, posB;
  double rr;

  const double cutOffVDWSquared = forceField.cutOffVDW * forceField.cutOffVDW;
  const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;

  EnergyStatus energy(components.size());
  double3x3 strainDerivativeTensor{};

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

          const double3 g = forceFactor.forceFactor * dr;

          it1->gradient += g;
          it2->gradient -= g;

          strainDerivativeTensor.ax += g.x*dr.x;
          strainDerivativeTensor.bx += g.y*dr.x;
          strainDerivativeTensor.cx += g.z*dr.x;
  
          strainDerivativeTensor.ay += g.x*dr.y;
          strainDerivativeTensor.by += g.y*dr.y;
          strainDerivativeTensor.cy += g.z*dr.y;
 
          strainDerivativeTensor.az += g.x*dr.z;
          strainDerivativeTensor.bz += g.y*dr.z;
          strainDerivativeTensor.cz += g.z*dr.z;
        }
        if (!noCharges && rr < cutOffChargeSquared)
        {
          double r = std::sqrt(rr);
          double scaling = it1->scalingCoulomb * it2->scalingCoulomb;
          ForceFactor energyFactor = potentialCoulombGradient(forceField, scaling, r, chargeA, chargeB);

          energy(compA, compB).CoulombicReal += 0.5 * EnergyFactor(energyFactor.energy, 0);
          energy(compB, compA).CoulombicReal += 0.5 * EnergyFactor(energyFactor.energy, 0);

          const double3 g = energyFactor.forceFactor * dr;

          it1->gradient += g;
          it2->gradient -= g;

          strainDerivativeTensor.ax += g.x*dr.x;
          strainDerivativeTensor.bx += g.y*dr.x;
          strainDerivativeTensor.cx += g.z*dr.x;

          strainDerivativeTensor.ay += g.x*dr.y;
          strainDerivativeTensor.by += g.y*dr.y;
          strainDerivativeTensor.cy += g.z*dr.y;

          strainDerivativeTensor.az += g.x*dr.z;
          strainDerivativeTensor.bz += g.y*dr.z;
          strainDerivativeTensor.cz += g.z*dr.z;
        }
      }
    }
  }

  return {energy, strainDerivativeTensor};
}

