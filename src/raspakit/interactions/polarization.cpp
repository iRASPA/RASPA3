module;

module interactions_polarization;

import std;

import int3;
import double3;
import double3x3;
import atom;
import potential_pair_derivatives;
import potential_pair_coulomb;
import simulationbox;
import energy_status;
import energy_status_inter;
import units;
import running_energy;
import framework;
import component;
import forcefield;

// The polarization coupling of an atom is scaled by its Coulomb scaling factor: a fractional (CFCMC)
// molecule decouples smoothly from the field as lambda decreases and feels no polarization at all below
// lambda = 1/2 -- where the VDW interaction is still (partially) off and the atom could approach a charge
// singularity where |E|^2 diverges. For integer molecules scalingCoulomb == 1 and nothing changes.
RunningEnergy Interactions::computePolarizationEnergyDifference(const ForceField &forceField,
                                                                std::span<double3> electricFieldNew,
                                                                std::span<double3> electricField,
                                                                std::span<Atom> moleculeAtomPositions)
{
  RunningEnergy energy;

  for (std::size_t i = 0; i < moleculeAtomPositions.size(); ++i)
  {
    std::size_t type = moleculeAtomPositions[i].type;
    double polarizability = moleculeAtomPositions[i].scalingCoulomb * forceField.pseudoAtoms[type].polarizability /
                            Units::CoulombicConversionFactor;
    energy.polarization -= 0.5 * polarizability * double3::dot(electricFieldNew[i], electricFieldNew[i]);
  }

  for (std::size_t i = 0; i < moleculeAtomPositions.size(); ++i)
  {
    std::size_t type = moleculeAtomPositions[i].type;
    double polarizability = moleculeAtomPositions[i].scalingCoulomb * forceField.pseudoAtoms[type].polarizability /
                            Units::CoulombicConversionFactor;
    energy.polarization += 0.5 * polarizability * double3::dot(electricField[i], electricField[i]);
  }

  return energy;
}

RunningEnergy Interactions::computePolarizationEnergyDifference(const ForceField &forceField,
                                                                std::span<double3> electricFieldNew,
                                                                std::span<double3> electricField,
                                                                std::span<Atom> moleculeAtomPositionsNew,
                                                                std::span<Atom> moleculeAtomPositionsOld)
{
  RunningEnergy energy;

  for (std::size_t i = 0; i < moleculeAtomPositionsNew.size(); ++i)
  {
    std::size_t type = moleculeAtomPositionsNew[i].type;
    double polarizability = moleculeAtomPositionsNew[i].scalingCoulomb * forceField.pseudoAtoms[type].polarizability /
                            Units::CoulombicConversionFactor;
    energy.polarization -= 0.5 * polarizability * double3::dot(electricFieldNew[i], electricFieldNew[i]);
  }

  for (std::size_t i = 0; i < moleculeAtomPositionsOld.size(); ++i)
  {
    std::size_t type = moleculeAtomPositionsOld[i].type;
    double polarizability = moleculeAtomPositionsOld[i].scalingCoulomb * forceField.pseudoAtoms[type].polarizability /
                            Units::CoulombicConversionFactor;
    energy.polarization += 0.5 * polarizability * double3::dot(electricField[i], electricField[i]);
  }

  return energy;
}

RunningEnergy Interactions::computePolarizationEnergyNeighborDifference(const ForceField &forceField,
                                                                       std::span<const double3> electricField,
                                                                       std::span<const double3> electricFieldDelta,
                                                                       std::span<const Atom> moleculeAtoms)
{
  RunningEnergy energy;

  for (std::size_t i = 0; i < moleculeAtoms.size(); ++i)
  {
    std::size_t type = moleculeAtoms[i].type;
    double polarizability = moleculeAtoms[i].scalingCoulomb * forceField.pseudoAtoms[type].polarizability /
                            Units::CoulombicConversionFactor;
    double3 electricFieldOld = electricField[i];
    double3 electricFieldNew = electricFieldOld + electricFieldDelta[i];
    energy.polarization -= 0.5 * polarizability *
                           (double3::dot(electricFieldNew, electricFieldNew) -
                            double3::dot(electricFieldOld, electricFieldOld));
  }

  return energy;
}

double Interactions::computePolarizationDUdlambda(const ForceField &forceField, const SimulationBox &simulationBox,
                                                  std::span<const Atom> moleculeAtoms,
                                                  std::span<const double3> electricField, std::size_t groupId)
{
  if (!forceField.computePolarization || groupId == 0) return 0.0;

  double value = 0.0;

  // Coupling term of the group's own atoms: d/ds of -1/2 s alpha |E|^2. The field an atom feels
  // excludes its own molecule, so it does not depend on the atom's own coupling.
  for (std::size_t i = 0; i < moleculeAtoms.size(); ++i)
  {
    if (static_cast<std::size_t>(moleculeAtoms[i].groupId) != groupId) continue;
    std::size_t type = static_cast<std::size_t>(moleculeAtoms[i].type);
    double polarizability = forceField.pseudoAtoms[type].polarizability / Units::CoulombicConversionFactor;
    value -= 0.5 * polarizability * double3::dot(electricField[i], electricField[i]);
  }

  // Field term: the group's atoms produce a real-space field at other molecules that is proportional
  // to the group's Coulomb scaling. Only present when molecule-molecule polarization is active; the
  // reciprocal part of the stored field is built from the fixed framework only and is lambda-independent.
  if (!forceField.useCharge || forceField.omitInterInteractions || forceField.omitInterPolarization) return value;

  std::vector<std::size_t> groupAtoms;
  for (std::size_t i = 0; i < moleculeAtoms.size(); ++i)
  {
    if (static_cast<std::size_t>(moleculeAtoms[i].groupId) == groupId) groupAtoms.push_back(i);
  }
  if (groupAtoms.empty()) return value;

  const double cutOffChargeSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;

  for (std::size_t j = 0; j < moleculeAtoms.size(); ++j)
  {
    std::size_t type = static_cast<std::size_t>(moleculeAtoms[j].type);
    double polarizability = moleculeAtoms[j].scalingCoulomb * forceField.pseudoAtoms[type].polarizability /
                            Units::CoulombicConversionFactor;
    if (polarizability == 0.0) continue;

    double3 unscaledGroupField{};
    for (std::size_t s : groupAtoms)
    {
      if (moleculeAtoms[s].moleculeId == moleculeAtoms[j].moleculeId) continue;

      double3 dr = moleculeAtoms[j].position - moleculeAtoms[s].position;
      dr = simulationBox.applyPeriodicBoundaryConditions(dr);
      double rr = double3::dot(dr, dr);
      if (rr < cutOffChargeSquared)
      {
        double r = std::sqrt(rr);
        Potentials::PairDerivatives<1> gradient = Potentials::potentialCoulomb<1>(forceField, 1.0, 1.0, r, 1.0, 1.0);
        unscaledGroupField -= moleculeAtoms[s].charge * gradient.firstDerivativeFactor * dr;
      }
    }
    value -= polarizability * double3::dot(electricField[j], unscaledGroupField);
  }

  return value;
}
