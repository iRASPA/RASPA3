#include <gtest/gtest.h>

import std;

import double3;
import units;
import atom;
import forcefield;
import component;
import system;
import simulationbox;
import running_energy;
import interactions_ewald;
import potential_coulomb_real_space;

// 'energyDifferenceEwaldFourier' is called with a span holding a single molecule by every ordinary
// Monte Carlo move, and with a span holding several molecules at once by the reaction moves and by
// the Gibbs CB/CFCMC flag swap. The multi-molecule call has to agree with the sequence of
// single-molecule calls that the chaining callers perform, which requires the intramolecular
// exclusion to be matched on moleculeId instead of being taken over every pair in the span.

namespace
{
using EikVector = std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>>;

// Exclusion terms of the pairs whose two atoms belong to different molecules: precisely the spurious
// contribution that appears when a multi-molecule span is treated as one molecule. Used to assert
// that the comparison below is actually sensitive to the difference.
double crossMoleculeExclusion(const ForceField& forceField, const SimulationBox& box, std::span<const Atom> atoms)
{
  const double cutOffSquared = forceField.cutOffCoulomb * forceField.cutOffCoulomb;
  double energy = 0.0;
  for (std::size_t i = 0; i + 1 < atoms.size(); ++i)
  {
    for (std::size_t j = i + 1; j != atoms.size(); ++j)
    {
      if (atoms[i].moleculeId == atoms[j].moleculeId) continue;

      double3 dr = box.applyPeriodicBoundaryConditions(atoms[i].position - atoms[j].position);
      double rr = double3::dot(dr, dr);
      double r = std::sqrt(rr);
      double scaling = atoms[i].scalingCoulomb * atoms[j].scalingCoulomb;
      double prefactor = Units::CoulombicConversionFactor * atoms[i].charge * atoms[j].charge;

      if (forceField.usesEwaldFourier())
      {
        energy += scaling * prefactor * Potentials::ewaldExclusionFactors(forceField.EwaldAlpha, scaling, r).potential;
      }
      else if (rr < cutOffSquared)
      {
        energy += scaling * prefactor * (Potentials::coulombRealSpaceFactors(forceField, r).potential - 1.0 / r);
      }
    }
  }
  return energy;
}

// The reference treatment: replace one molecule at a time, threading the structure factor of each
// call into the next, as the chaining callers (pair swap, group swap, Gibbs swap) do.
RunningEnergy chainedDifference(System& system, const std::vector<std::vector<Atom>>& newMolecules,
                                const std::vector<std::vector<Atom>>& oldMolecules)
{
  RunningEnergy total{};
  EikVector working = system.storedEik;
  EikVector trial;
  for (std::size_t k = 0; k != oldMolecules.size(); ++k)
  {
    total += Interactions::energyDifferenceEwaldFourier(system.eik_x, system.eik_y, system.eik_z, system.eik_xy,
                                                        working, trial, system.forceField, system.simulationBox,
                                                        newMolecules[k], oldMolecules[k]);
    working = trial;
  }
  return total;
}

// Three CO2 molecules: the first two are the ones the move changes, the third stays put so that the
// structure factor the difference is taken against is not trivial. The second molecule is fractional
// so that the per-group dU/dlambda is exercised as well.
System makeThreeMoleculeSystem(ForceField::ChargeMethod chargeMethod)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  forceField.automaticEwald = false;
  forceField.EwaldAlpha = 0.25;
  forceField.chargeMethod = chargeMethod;

  Component c = Component::makeCO2(forceField, 0, true);
  System system = System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e4, 1.0, {}, {c}, {}, {3}, 5);

  std::span<Atom> atomData = system.spanOfMoleculeAtoms();
  const std::array<double3, 3> centers{double3(-2.0, 0.0, 0.0), double3(2.0, 0.0, 0.0), double3(0.0, 8.0, 0.0)};
  for (std::size_t m = 0; m != centers.size(); ++m)
  {
    atomData[3 * m + 0].position = centers[m] + double3(0.0, 0.0, 1.149);
    atomData[3 * m + 1].position = centers[m];
    atomData[3 * m + 2].position = centers[m] + double3(0.0, 0.0, -1.149);
  }

  for (std::size_t i = 3; i != 6; ++i) atomData[i].setScalingToFractional(0.4, 1);

  (void)Interactions::computeEwaldFourierEnergy(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
      system.forceField, system.simulationBox, system.components, system.numberOfMoleculesPerComponent,
      system.spanOfMoleculeAtoms());

  return system;
}

// The trial state: both molecules move, and the fractional one also changes its coupling.
std::vector<Atom> trialStateOf(std::span<const Atom> molecule, double3 displacement, std::optional<double> lambda)
{
  std::vector<Atom> atoms(molecule.begin(), molecule.end());
  for (Atom& atom : atoms)
  {
    atom.position += displacement;
    if (lambda.has_value()) atom.setScalingToFractional(lambda.value(), 1);
  }
  return atoms;
}

void expectCombinedMatchesChained(ForceField::ChargeMethod chargeMethod)
{
  System system = makeThreeMoleculeSystem(chargeMethod);
  std::span<const Atom> atomData = system.spanOfMoleculeAtoms();

  const std::vector<std::vector<Atom>> oldMolecules{
      std::vector<Atom>(atomData.begin(), atomData.begin() + 3),
      std::vector<Atom>(atomData.begin() + 3, atomData.begin() + 6)};
  const std::vector<std::vector<Atom>> newMolecules{
      trialStateOf(oldMolecules[0], double3(0.3, -0.2, 0.1), std::nullopt),
      trialStateOf(oldMolecules[1], double3(-0.15, 0.25, -0.2), 0.7)};

  std::vector<Atom> oldAtoms;
  std::vector<Atom> newAtoms;
  for (std::size_t k = 0; k != oldMolecules.size(); ++k)
  {
    oldAtoms.insert(oldAtoms.end(), oldMolecules[k].begin(), oldMolecules[k].end());
    newAtoms.insert(newAtoms.end(), newMolecules[k].begin(), newMolecules[k].end());
  }

  // The two molecules sit 4 Angstrom apart, so treating the span as a single molecule would add a
  // clearly non-zero spurious exclusion; without this the comparison below could pass vacuously.
  const double spurious = crossMoleculeExclusion(system.forceField, system.simulationBox, oldAtoms) -
                          crossMoleculeExclusion(system.forceField, system.simulationBox, newAtoms);
  EXPECT_GT(std::abs(spurious), 1.0);

  const RunningEnergy combined = Interactions::energyDifferenceEwaldFourier(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.trialEik, system.forceField,
      system.simulationBox, newAtoms, oldAtoms);

  const RunningEnergy chained = chainedDifference(system, newMolecules, oldMolecules);

  EXPECT_NEAR(combined.ewald_exclusion, chained.ewald_exclusion, 1e-8);
  EXPECT_NEAR(combined.ewald_fourier, chained.ewald_fourier, 1e-8);
  EXPECT_NEAR(combined.ewald_self, chained.ewald_self, 1e-8);
  EXPECT_NEAR(combined.dudlambdaEwald[0], chained.dudlambdaEwald[0], 1e-8);
  EXPECT_NEAR(combined.potentialEnergy(), chained.potentialEnergy(), 1e-8);
}
}  // namespace

TEST(ewald_exclusion_multi_molecule, ewald_combined_span_matches_chained_single_molecule_calls)
{
  expectCombinedMatchesChained(ForceField::ChargeMethod::Ewald);
}

TEST(ewald_exclusion_multi_molecule, wolf_combined_span_matches_chained_single_molecule_calls)
{
  expectCombinedMatchesChained(ForceField::ChargeMethod::Wolf);
}
