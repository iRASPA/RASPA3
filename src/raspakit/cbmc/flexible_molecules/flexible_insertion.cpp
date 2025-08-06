module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <optional>
#include <span>
#include <tuple>
#include <type_traits>
#include <vector>
#endif

module cbmc_flexible_insertion;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import randomnumbers;
import units;
import component;
import molecule;
import atom;
import double3;
import simd_quatd;
import double3x3;
import simulationbox;
import energy_status;
import cbmc_first_bead_data;
import cbmc_chain_data;
import cbmc_util;
import cbmc_multiple_first_bead;
import cbmc_interactions;
import cbmc_generate_trialorientations_sphere;
import cbmc_generate_trialorientations_mc;
import forcefield;
import energy_factor;
import running_energy;
import framework;
import component;
import interpolation_energy_grid;
import connectivity_table;
import internal_potentials;
import bond_potential;

// atoms is a recentered copy of the molecule (recentered around the starting bead)
[[nodiscard]] std::optional<ChainData> CBMC::growFlexibleMoleculeSwapInsertion(
    RandomNumber &random, Component &component, bool hasExternalField, const std::vector<Component> &components,
    const ForceField &forceField, const SimulationBox &simulationBox,
    const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> moleculeAtoms, double beta, double cutOffFrameworkVDW, double cutOffMoleculeVDW,
    double cutOffCoulomb, std::size_t selectedComponent, std::size_t selectedMolecule, double scaling, bool groupId,
    bool isFractional, std::size_t numberOfTrialDirections)
{
  std::size_t startingBead = components[selectedComponent].startingBead;
  Atom firstBead = components[selectedComponent].atoms[startingBead];
  firstBead.moleculeId = static_cast<std::uint32_t>(selectedMolecule);
  firstBead.groupId = groupId;
  firstBead.isFractional = isFractional;
  firstBead.setScaling(scaling);

  std::optional<FirstBeadData> const firstBeadData = CBMC::growMoleculeMultipleFirstBeadSwapInsertion(
      random, component, hasExternalField, forceField, simulationBox, interpolationGrids, framework, frameworkAtoms,
      moleculeAtoms, beta, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb, firstBead, numberOfTrialDirections);

  if (!firstBeadData) return std::nullopt;

  if (components[selectedComponent].atoms.size() == 1)
  {
    return ChainData(Molecule(double3(firstBeadData->atom.position), simd_quatd(0.0, 0.0, 0.0, 1.0),
                              component.totalMass, component.componentId, component.definedAtoms.size()),
                     {firstBeadData->atom}, firstBeadData->energies, firstBeadData->RosenbluthWeight, 0.0);
  }

  // place the molecule centered around the first bead at 'firstBeadData->atom.position'
  std::vector<Atom> atoms = components[selectedComponent].atoms;
  atoms[0] = firstBead;
  std::for_each(atoms.begin(), atoms.end(),
                [&](Atom &atom)
                {
                  atom.position +=
                      firstBeadData->atom.position - components[selectedComponent].atoms[startingBead].position;
                  atom.moleculeId = static_cast<std::uint32_t>(selectedMolecule);
                  atom.groupId = groupId;
                  atom.isFractional = isFractional;
                  atom.setScaling(scaling);
                });

  std::optional<ChainData> const flexibleChainData = CBMC::growFlexibleMoleculeChainInsertion(
      random, component, hasExternalField, forceField, simulationBox, interpolationGrids, framework, frameworkAtoms,
      moleculeAtoms, beta, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb, startingBead, atoms,
      numberOfTrialDirections, selectedMolecule, scaling, groupId, isFractional, components, selectedComponent);

  if (!flexibleChainData) return std::nullopt;

  return ChainData(flexibleChainData->molecule, flexibleChainData->atom,
                   firstBeadData->energies + flexibleChainData->energies,
                   firstBeadData->RosenbluthWeight * flexibleChainData->RosenbluthWeight, 0.0);
}

[[nodiscard]] std::optional<ChainData> CBMC::growFlexibleMoleculeChainInsertion(
    RandomNumber &random, Component &component, bool hasExternalField, const ForceField &forceField,
    const SimulationBox &simulationBox, const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> moleculeAtoms, double beta, double cutOffFrameworkVDW, double cutOffMoleculeVDW,
    double cutOffCoulomb, std::size_t startingBead, std::vector<Atom> molecule, std::size_t numberOfTrialDirections,
    std::size_t selectedMolecule, double scaling, bool groupId, bool isFractional,
    const std::vector<Component> &components, std::size_t selectedComponent)
{
  std::size_t numberOfBeads = component.connectivityTable.numberOfBeads;
  std::vector<std::vector<Atom>> trialPositions{};
  trialPositions.reserve(numberOfTrialDirections);
  std::vector<Atom> chain_atoms(molecule.begin(), molecule.end());

  std::vector<std::size_t> beads_already_placed{ component.startingBead };
  const Atom first_bead = molecule[component.startingBead];
  chain_atoms[component.startingBead] = first_bead;

  double chain_rosen_bluth_weight = 1.0;
  RunningEnergy chain_external_energies{};

  //for(std::size_t i = 0; i != molecule.size(); ++i)
  //{
  //  std::print("atoms: {}: {}\n", i, molecule[i]);
  //}
  //for(std::size_t i = 0; i != chain_atoms.size(); ++i)
  //{
  //  std::print("chain_atoms: {}: {}\n", i, chain_atoms[i]);
  //}

  do
  {
    auto [previous_bead, current_bead, nextBeads] = component.connectivityTable.nextBeads(beads_already_placed);

    Potentials::InternalPotentials internalPotentials = 
      component.internalPotentials.filteredInteractions(numberOfBeads, beads_already_placed, nextBeads);

    if(!previous_bead.has_value())
    {
      // case: growing a single bond with no previous beads
      //       for example: dimer, or starting in the middle of a linear chain
      
      if(internalPotentials.bonds.size() != 1)
      {
        throw std::runtime_error(std::format("[CBMC]: multiple bonds detected in 'growFlexibleMoleculeChainInsertion'\n"));
      }

      const BondPotential bond = internalPotentials.bonds.front();

      Atom argument = chain_atoms[nextBeads[0]];
      argument.position = first_bead.position;

      trialPositions = generateTrialOrientationsSimpleSphere(random, beta,  argument, numberOfTrialDirections, bond);
      //for(const std::vector<Atom> &t : trialPositions)
      //{
      //  for(const Atom &s : t)
      //  {
      //    std::print("trial_orientations sphere: {}\n", s);
      //  }
      //}
    }
    else
    {

      std::vector<Atom> trial_orientations = generateTrialOrientationsMonteCarloScheme(random, beta, 
                                                    component, chain_atoms, 
                                                    previous_bead.value(), current_bead, nextBeads,
                                                    numberOfTrialDirections, internalPotentials);

      //for(const Atom &s : trial_orientations)
      //{
      //  std::print("trial_orientations: {}\n",s);
      //}

      // Rotate around to obtain the other trial-orientations
      trialPositions = {};
      double3 last_bond_vector = (chain_atoms[previous_bead.value()].position - chain_atoms[current_bead].position).normalized();

      for(std::size_t i = 0; i != numberOfTrialDirections; ++i)
      {
        double random_angle = (2.0 * random.uniform() - 1.0) * std::numbers::pi;

        // bond_vector.rotateAroundAxis(last_bond_vector, random_angle)
        std::vector<Atom> rotated_atoms = trial_orientations;
        std::for_each(rotated_atoms.begin(), rotated_atoms.end(), [&](Atom &atom)
          {
            atom.position = chain_atoms[current_bead].position + last_bond_vector.rotateAroundAxis(atom.position - chain_atoms[current_bead].position, random_angle);
          });
        trialPositions.push_back(rotated_atoms);
      }
    }    

    // compute the external-energies for the next-beads
    const std::vector<std::pair<std::vector<Atom>, RunningEnergy>> externalEnergies = 
      CBMC::computeExternalNonOverlappingEnergies(component, hasExternalField, forceField, simulationBox,
                          interpolationGrids, framework, frameworkAtoms, moleculeAtoms,
                          cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb, trialPositions, -1);


    if (externalEnergies.empty()) return std::nullopt;


    std::vector<double> logBoltmannFactors{};
    std::transform(externalEnergies.begin(), externalEnergies.end(), std::back_inserter(logBoltmannFactors),
                   [&](const std::pair<std::vector<Atom>, RunningEnergy> &v)
                   { return -beta * std::get<1>(v).potentialEnergy(); });

    std::size_t selected = CBMC::selectTrialPosition(random, logBoltmannFactors);

    double rosen_bluth_weight = std::accumulate(logBoltmannFactors.begin(), logBoltmannFactors.end(), 0.0,
                                              [](const double &acc, const double &logBoltmannFactor)
                                              { return acc + std::exp(logBoltmannFactor); });

    if (rosen_bluth_weight < forceField.minimumRosenbluthFactor) return std::nullopt;

    chain_rosen_bluth_weight *= rosen_bluth_weight / static_cast<double>(numberOfTrialDirections);

    const std::vector<Atom> &selected_trial_atoms = externalEnergies[selected].first;
    chain_external_energies += externalEnergies[selected].second;

    // add 'nextBeads' to 'beads_already_placed'
    beads_already_placed.insert(beads_already_placed.end(), nextBeads.begin(), nextBeads.end());

    // add the selected atoms 
    for(std::size_t i = 0; i != nextBeads.size(); ++i)
    {
      std::size_t index = nextBeads[i];
      chain_atoms[index] = selected_trial_atoms[i];
    }

  } while(beads_already_placed.size() < component.connectivityTable.numberOfBeads);

  // recompute all the internal interactions
  RunningEnergy internal_energies = component.internalPotentials.computeInternalEnergies(chain_atoms);

  //for(std::size_t i = 0; i != chain_atoms.size(); ++i)
  //{
  //  std::print("chain_atoms: {}: {}\n", i, chain_atoms[i]);
  //}
  // std::print("\n");

  //std::print("lengths: {} {}\n", (chain_atoms[1].position - chain_atoms[0].position).length(),
  //                               (chain_atoms[2].position - chain_atoms[1].position).length());

  //double check_angle = double3::angle(chain_atoms[0].position, chain_atoms[1].position, chain_atoms[2].position);
  //std::print("Angle CHECK: {}\n", check_angle * Units::RadiansToDegrees);


  //std::exit(0);

  return ChainData({}, chain_atoms,  chain_external_energies + internal_energies, 
                   chain_rosen_bluth_weight, 0.0);
}
