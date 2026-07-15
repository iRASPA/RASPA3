module;

module cbmc_interactions_intermolecular;

import std;

import energy_status;
import potential_pair_derivatives;
import potential_pair_vdw;
import potential_pair_coulomb;
import potential_correction_vdw;
import simulationbox;
import double3;
import double3x3;
import forcefield;
import atom;
import energy_status_inter;
import running_energy;
import units;
import threadpool;

[[nodiscard]] std::optional<RunningEnergy> CBMC::computeInterMolecularEnergy(
    const ForceField &forceField, const SimulationBox &simulationBox, std::span<const Atom> moleculeAtoms,
    double cutOffVDW, double cutOffCoulomb, std::span<Atom> atoms, std::make_signed_t<std::size_t> skip,
    std::make_signed_t<std::size_t> skipBackgroundMolecule) noexcept
{
  double3 dr, s, t;
  double rr;

  RunningEnergy energySum;

  bool useCharge = forceField.useCharge;

  const double overlapCriteria = forceField.energyOverlapCriteria;
  const double cutOffVDWSquared = cutOffVDW * cutOffVDW;
  const double cutOffChargeSquared = cutOffCoulomb * cutOffCoulomb;

  for (std::span<const Atom>::iterator it1 = moleculeAtoms.begin(); it1 != moleculeAtoms.end(); ++it1)
  {
    std::size_t molA = static_cast<std::size_t>(it1->moleculeId);
    if (skipBackgroundMolecule >= 0 && molA == static_cast<std::size_t>(skipBackgroundMolecule))
    {
      continue;
    }
    double3 posA = it1->position;
    std::size_t typeA = static_cast<std::size_t>(it1->type);
    std::uint8_t groupIdA = it1->groupId;
    double scalingVDWA = it1->scalingVDW;
    double scalingCoulombA = it1->scalingCoulomb;
    double chargeA = it1->charge;

    for (int index = 0; const Atom &atom : atoms)
    {
      if (index != skip)
      {
        double3 posB = atom.position;
        std::size_t molB = static_cast<std::size_t>(atom.moleculeId);
        std::size_t typeB = static_cast<std::size_t>(atom.type);
        std::uint8_t groupIdB = atom.groupId;
        double scalingVDWB = atom.scalingVDW;
        double scalingCoulombB = atom.scalingCoulomb;
        double chargeB = atom.charge;

        if (molA != molB)
        {
          dr = posA - posB;
          dr = simulationBox.applyPeriodicBoundaryConditions(dr);
          rr = double3::dot(dr, dr);

          if (rr < cutOffVDWSquared)
          {
            Potentials::PairDerivatives<0> energyFactor = Potentials::potentialVDW<0>(
                forceField, scalingVDWA, scalingVDWB, rr, typeA, typeB);
            if (energyFactor.energy > overlapCriteria) return std::nullopt;

            energySum.moleculeMoleculeVDW += energyFactor.energy;
            energySum.addDudlambdaVDW(groupIdA, groupIdB, scalingVDWA, scalingVDWB, energyFactor.dUdlambda);
          }
          if (useCharge && rr < cutOffChargeSquared)
          {
            double r = std::sqrt(rr);
            Potentials::PairDerivatives<0> energyFactor = Potentials::potentialCoulomb<0>(
                forceField, scalingCoulombA, scalingCoulombB, r, chargeA, chargeB);

            energySum.moleculeMoleculeCharge += energyFactor.energy;
            energySum.addDudlambdaCharge(groupIdA, groupIdB, scalingCoulombA, scalingCoulombB, energyFactor.dUdlambda);
          }
        }
      }
      ++index;
    }
  }

  return energySum;
}
