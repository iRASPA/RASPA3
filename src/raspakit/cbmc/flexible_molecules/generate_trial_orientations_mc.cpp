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
    
module cbmc_generate_trialorientations_mc;
    
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
import forcefield;  
import energy_factor;
import running_energy;
import framework;
import component;
import interpolation_energy_grid;
import bond_potential;
import intra_molecular_potentials;
import cbmc_move_statistics;

std::vector<Atom> CBMC::generateTrialOrientationsMonteCarloScheme(RandomNumber &random, double beta,
            Component &component,
            const std::vector<Atom> molecule_atoms, std::size_t previousBead,
            std::size_t currentBead, std::vector<std::size_t> nextBeads, std::size_t numberOfTrialDirections,
            const Potentials::IntraMolecularPotentials &intraMolecularInteractions)
{
  std::vector<Atom> chain_atoms(molecule_atoms);

  // Calculate initial positions for the next beads
  if(!component.grownAtoms.empty())
  {
    // Restore from previous configuration stored in component
    double3 last_bond_vector = chain_atoms[previousBead].position - chain_atoms[currentBead].position;
    double3 v = component.grownAtoms[previousBead].position - component.grownAtoms[currentBead].position;

    double3x3 rotation_matrix = double3x3::computeRotationMatrix(v, last_bond_vector);

    for(std::size_t i = 0; i != nextBeads.size(); ++i)
    {
      std::optional<BondPotential> bond = intraMolecularInteractions.findBondPotential(currentBead, nextBeads[i]);

      double3 va = component.grownAtoms[nextBeads[i]].position - component.grownAtoms[currentBead].position;

      double3 vec = rotation_matrix * va;

      chain_atoms[nextBeads[i]].position = chain_atoms[currentBead].position + vec;
    }
  }
  else
  {
    // Generate positions from scratch
    for(std::size_t i = 0; i != nextBeads.size(); ++i)
    {
      std::optional<BondPotential> bond = intraMolecularInteractions.findBondPotential(currentBead, nextBeads[i]);

      double bond_length = bond.transform([&random, &beta](const BondPotential &bond) { return bond.generateBondLength(random, beta);}).value_or(1.54);
      double angle = intraMolecularInteractions.bends.empty() ? 120.0 * Units::DegreesToRadians :
                         intraMolecularInteractions.bends.front().generateBendAngle(random, beta);

      double3 last_bond_vector = (chain_atoms[previousBead].position - chain_atoms[currentBead].position).normalized();
      double3 vec = random.randomVectorOnCone(last_bond_vector, angle);

      chain_atoms[nextBeads[i]].position = chain_atoms[currentBead].position + bond_length * vec;
    }
  }

  std::size_t number_of_trials = 600;

  std::vector<double> move_probabilities{};
  if( intraMolecularInteractions.bends.size() == 1 )
  {
    move_probabilities = { 0.6, 0.4, 0.0 };
  }
  else
  {
    move_probabilities = { 0.4, 0.4, 0.2 };
  }

  // pick the move accoording to the prescribed weights
  MoveType move_type = MoveType(random.categoricalDistribution(move_probabilities));

  RunningEnergy current_bond_energy = intraMolecularInteractions.calculateBondSmallMCEnergies(chain_atoms);
  RunningEnergy current_bend_energy = intraMolecularInteractions.calculateBendSmallMCEnergies(chain_atoms);

  for(std::size_t trial = 0; trial != number_of_trials; ++trial)
  {
    switch(move_type)
    {
      case MoveType::BondLengthChange: 
      {
        // Move: change bond-length
        // bond-energies are affected by only bond-terms
        // bond-bond needs to be handled after growing the molecule and correcting for it

        component.cbmc_moves_statistics[currentBead].bondLengthChange.counts += 1;
        component.cbmc_moves_statistics[currentBead].bondLengthChange.totalCounts += 1;

        std::size_t selected_next_bead = nextBeads[random.uniform_integer(0, nextBeads.size() - 1)];
        double3 current_bond_vector = chain_atoms[selected_next_bead].position - chain_atoms[currentBead].position;
        double current_bond_length = current_bond_vector.length();

        double3 saved_current_position = chain_atoms[selected_next_bead].position;
        double max_change = component.cbmc_moves_statistics[currentBead].bondLengthChange.maxChange;
        double new_bond_length = current_bond_length + (2.0 * random.uniform() - 1.0) * max_change;
        chain_atoms[selected_next_bead].position = chain_atoms[currentBead].position + (new_bond_length / current_bond_length) * current_bond_vector;

        RunningEnergy new_bond_energy = intraMolecularInteractions.calculateBondSmallMCEnergies(chain_atoms);

        double energy_difference = (new_bond_energy - current_bond_energy).potentialEnergy();

        component.cbmc_moves_statistics[currentBead].bondLengthChange.constructed += 1;
        component.cbmc_moves_statistics[currentBead].bondLengthChange.totalConstructed += 1;

        // compute energy
        if(random.uniform() < (new_bond_length / current_bond_length) * (new_bond_length / current_bond_length) * std::exp(-beta * energy_difference))
        {
          // update component statistics
          current_bond_energy = new_bond_energy;

          component.cbmc_moves_statistics[currentBead].bondLengthChange.accepted += 1;
          component.cbmc_moves_statistics[currentBead].bondLengthChange.totalAccepted += 1;
        }
        else
        {
          chain_atoms[selected_next_bead].position = saved_current_position;
        }
        break;
      }
      case MoveType::BendAngleChange: 
      {
        // Move: change bend angle of one chosen bead
        // bend-energies are affected by: 1) bend-terms, 2) urey-bradley 3) bond-bend, 4) bend-bend, 5) inversion-bend, 6) out-of-plane bend
        
        component.cbmc_moves_statistics[currentBead].bendAngleChange.counts += 1;
        component.cbmc_moves_statistics[currentBead].bendAngleChange.totalCounts += 1;

        std::size_t selected_next_bead = nextBeads[random.uniform_integer(0, nextBeads.size() - 1)];

        double current_angle = double3::angle(chain_atoms[previousBead].position, 
                                              chain_atoms[currentBead].position,
                                              chain_atoms[selected_next_bead].position);
        double current_sinus = std::sin(current_angle);

        double3 last_bond_vector = (chain_atoms[previousBead].position - chain_atoms[currentBead].position).normalized();
        double3 negative_bond_vector = chain_atoms[currentBead].position - chain_atoms[selected_next_bead].position;
        double3 perpendicular_vector = double3::perpendicular(negative_bond_vector, last_bond_vector);

        double max_change = component.cbmc_moves_statistics[currentBead].bendAngleChange.maxChange;
        double random_angle = (2.0 * random.uniform() - 1.0) * max_change;
        double3 bond_vector = chain_atoms[selected_next_bead].position - chain_atoms[currentBead].position;
        double3 new_vec = perpendicular_vector.rotateAroundAxis(bond_vector, random_angle);

        double3 saved_current_position = chain_atoms[selected_next_bead].position;
        chain_atoms[selected_next_bead].position = chain_atoms[currentBead].position + new_vec;

        double new_angle = double3::angle(chain_atoms[previousBead].position,
                                          chain_atoms[currentBead].position,
                                          chain_atoms[selected_next_bead].position);
        double new_sinus = std::sin(new_angle);

        RunningEnergy new_bend_energy = intraMolecularInteractions.calculateBendSmallMCEnergies(chain_atoms);

        double energy_difference = (new_bend_energy - current_bend_energy).potentialEnergy();

        //std::print("random angle: {} new: {} {} old {} {} delta: {}\n", random_angle * Units::RadiansToDegrees,
        //                new_angle * Units::RadiansToDegrees, new_bend_energy.potentialEnergy(), 
        //                current_angle * Units::RadiansToDegrees, current_bend_energy.potentialEnergy(),
        //               energy_difference);
        
        component.cbmc_moves_statistics[currentBead].bendAngleChange.constructed += 1;
        component.cbmc_moves_statistics[currentBead].bendAngleChange.totalConstructed += 1;

        // compute energy
        if(random.uniform() < (new_sinus / current_sinus) * std::exp(-beta * energy_difference))
        {
          // update component statistics
          current_bend_energy = new_bend_energy;

          component.cbmc_moves_statistics[currentBead].bendAngleChange.accepted += 1;
          component.cbmc_moves_statistics[currentBead].bendAngleChange.totalAccepted += 1;
        }
        else
        {
          chain_atoms[selected_next_bead].position = saved_current_position;
        }
        break;
      }
      case MoveType::ConeChange: 
      {
        // Move: rotate one position on the cone while keeping the cone angle fixed
        // bend-energies are affected by: 1) bend-terms, 2) bond-bend, 3) bend-bend, 4) inversion-bend, 5) out-of-plane bend
        
        component.cbmc_moves_statistics[currentBead].conePositionChange.counts += 1;
        component.cbmc_moves_statistics[currentBead].conePositionChange.totalCounts += 1;

        std::size_t selected_next_bead = nextBeads[random.uniform_integer(0, nextBeads.size() - 1)];

        double3 last_bond_vector = (chain_atoms[previousBead].position - chain_atoms[currentBead].position).normalized();
        double3 bond_vector = chain_atoms[selected_next_bead].position - chain_atoms[currentBead].position;

        double max_change = component.cbmc_moves_statistics[currentBead].conePositionChange.maxChange;
        double random_angle = (2.0 * random.uniform() - 1.0) * max_change;
        double3 new_vec = last_bond_vector.rotateAroundAxis(bond_vector, random_angle);

        double3 saved_current_position = chain_atoms[selected_next_bead].position;
        chain_atoms[selected_next_bead].position = chain_atoms[currentBead].position + new_vec;

        double new_angle = double3::angle(chain_atoms[previousBead].position,
                                                   chain_atoms[currentBead].position,
                                                   chain_atoms[selected_next_bead].position);

        RunningEnergy new_bend_energy = intraMolecularInteractions.calculateBendSmallMCEnergies(chain_atoms);

        double energy_difference = (new_bend_energy - current_bend_energy).potentialEnergy();

        component.cbmc_moves_statistics[currentBead].conePositionChange.constructed += 1;
        component.cbmc_moves_statistics[currentBead].conePositionChange.totalConstructed += 1;

        // compute energy
        if(random.uniform() < std::exp(-beta * energy_difference))
        {
          // update component statistics
          current_bend_energy = new_bend_energy;

          component.cbmc_moves_statistics[currentBead].conePositionChange.accepted += 1;
          component.cbmc_moves_statistics[currentBead].conePositionChange.totalAccepted += 1;
        }
        else
        {
          chain_atoms[selected_next_bead].position = saved_current_position;
        }
        break;
      }
      default:
        std::unreachable();
    }
  }

  std::vector<Atom> next_bead_atoms(nextBeads.size());
  for(std::size_t i = 0; i != nextBeads.size(); ++i)
  {
    next_bead_atoms[i] = chain_atoms[nextBeads[i]];
  }


  return next_bead_atoms;
}
