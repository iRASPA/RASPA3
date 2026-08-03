module;

module energy_shared_molecular_void_fraction;

import std;

import uint3;
import double3;
import skspacegroupdatabase;
import crystal;
import pair_interactions;
import units;
import energy_shared_linear_probe;
import energy_shared_energy_backend;
import energy_shared_molecular_energy_grid;
import energy_shared_electrostatic_potential_grid;


MolecularVoidFraction::MolecularVoidFraction() {}


MolecularVoidFraction::~MolecularVoidFraction() {}


MolecularVoidFraction MolecularVoidFraction::fromGrid(const MolecularEnergyGrid &grid, const Crystal &framework,
                                                      double temperature)
{
  std::chrono::steady_clock::time_point time_begin = std::chrono::steady_clock::now();

  MolecularVoidFraction result;
  result.grid = grid;
  result.temperature = temperature;
  result.numberOfVoxels = grid.numberOfVoxels();

  if (grid.numberOfVoxels() == 0) return result;

  const double beta = 1.0 / (Units::KB * temperature);

  // Summed in double although the field is single, because the terms span many orders of magnitude: a point
  // in the wall contributes nothing at all and a deep site contributes a great deal, and adding those in
  // single precision would lose the small ones entirely. A point inside an atom sits at the ceiling and its
  // term underflows to zero, which is the right answer for a place the molecule cannot be.
  double total = 0.0;
  for (float freeEnergy : grid.freeEnergy) total += std::exp(-beta * static_cast<double>(freeEnergy));

  result.boltzmannAverage = total / static_cast<double>(grid.numberOfVoxels());
  result.readsAsFraction = result.boltzmannAverage <= 1.0;

  if (result.boltzmannAverage > 0.0)
  {
    result.excessChemicalPotential = -Units::KB * temperature * std::log(result.boltzmannAverage);
  }

  // The framework's density is what turns an average per unit volume into one per unit mass, which is the
  // form an isotherm is measured in.
  double density = 1e-3 * framework.mass /
                   (framework.unitCell.volume * Units::Angstrom * Units::Angstrom * Units::Angstrom *
                    Units::AvogadroConstant);
  if (density > 0.0)
  {
    result.henryCoefficient = result.boltzmannAverage / (Units::MolarGasConstant * temperature * density);
  }

  std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - time_begin;
  result.seconds = elapsed.count();

  return result;
}


void MolecularVoidFraction::run(const EnergyBackend &backend, const PairInteractions &interactions, const Crystal &framework, const LinearProbe &probe,
                                uint3 gridSize, std::size_t numberOfOrientations, double temperature,
                                bool useElectrostatics, double relativePrecision)
{
  MolecularField field = buildMolecularField(backend, interactions, framework, probe, gridSize, numberOfOrientations,
                                             temperature, useElectrostatics, relativePrecision);

  *this = MolecularVoidFraction::fromGrid(field.grid, framework, temperature);
  this->potential = field.potential;

  if (this->numberOfVoxels == 0) return;

  double3 spacing = this->grid.spacing();
  double density = 1e-3 * framework.mass /
                   (framework.unitCell.volume * Units::Angstrom * Units::Angstrom * Units::Angstrom *
                    Units::AvogadroConstant);

  std::ofstream myfile;
  myfile.open(framework.name + "." + probe.name + ".energy.vf." + this->grid.backend + ".txt");
  std::print(myfile, "# Energy-based void fraction for a rigid linear molecule\n");
  std::print(myfile, "# Crystal: {}\n", framework.name);
  std::print(myfile, "# Space-group HM-symbol: {}\n",
             SKSpaceGroupDataBase::spaceGroupData[framework.spaceGroupHallNumber].HMString());
  std::print(myfile, "# Number of framework atoms: {}\n", framework.atoms.size());
  std::print(myfile, "# Crystal volume: {} [Å³]\n", framework.unitCell.volume);
  std::print(myfile, "# Crystal mass: {} [g/mol]\n", framework.mass);
  std::print(myfile, "# Crystal density: {} [kg/m³]\n", density);
  std::print(myfile, "# Molecule: {}, {} sites, {:.4f} [Å] end to end\n", probe.name, probe.sites.size(),
             probe.length());
  std::print(myfile, "# Orientations sampled: {} over the {}\n", this->grid.numberOfOrientations,
             this->grid.overHemisphere ? "hemisphere, the molecule being the same end for end" : "whole sphere");
  std::print(myfile, "# Grid: {} x {} x {} points, spacing {:.5f} x {:.5f} x {:.5f} [Å], far face left out\n",
             this->grid.gridSize.x, this->grid.gridSize.y, this->grid.gridSize.z, spacing.x, spacing.y, spacing.z);
  std::print(myfile, "# Temperature: {} [K]\n", temperature);
  std::print(myfile, "# Cutoff: {} [Å]\n", this->grid.cutOff);

  if (this->grid.chargesIncluded)
  {
    std::print(myfile, "# Electrostatics: Ewald, split at alpha = {:.5f} [1/Å], {} wave vectors\n",
               this->potential.alpha, this->potential.numberOfWaveVectors);
  }
  if (this->grid.chargesIgnored)
  {
    std::print(myfile, "#\n");
    std::print(myfile, "# WARNING: this molecule carries partial charges and they have not been acted on.\n");
  }

  std::print(myfile, "# Timing ({}): {} [s] for the landscape\n", this->grid.backend, this->grid.seconds);
  std::print(myfile, "# CPU Timing: {} [s] for the average\n", this->seconds);
  std::print(myfile, "#\n");
  std::print(myfile, "# The average of exp(-U/kT) over positions and over orientations both. Averaging over\n");
  std::print(myfile, "# orientations is what the free energy already is, so this is the spatial average of\n");
  std::print(myfile, "# exp(-A/kT) over the very field the barrier is read from, and the two agree by\n");
  std::print(myfile, "# construction rather than by coincidence.\n");
  std::print(myfile, "\n");

  std::print(myfile, "{:36} {:16.8f}\n", "Average of exp(-A/kT):", this->boltzmannAverage);
  std::print(myfile, "{:36} {:16.4f} [K] {:12.4f} [kJ/mol]\n", "Excess chemical potential:",
             this->excessChemicalPotential * Units::EnergyToKelvin,
             this->excessChemicalPotential * Units::EnergyToKJPerMol);
  std::print(myfile, "{:36} {:16.6e} [mol/kg/Pa]\n", "Henry coefficient:", this->henryCoefficient);
  std::print(myfile, "\n");

  if (this->readsAsFraction)
  {
    std::print(myfile, "# The average is at or below one, so it may be read as a void fraction: the share of\n");
    std::print(myfile, "# the cell the molecule finds open to it, weighted by how welcome it is there.\n");
  }
  else
  {
    std::print(myfile, "# The average is above one, so it is NOT a void fraction and should not be quoted as\n");
    std::print(myfile, "# one. A molecule that binds is found in the framework more often than in the same\n");
    std::print(myfile, "# volume of gas, and the average measures that excess rather than any share of space.\n");
    std::print(myfile, "# Only a probe barely held at all, helium being the usual choice, gives an average that\n");
    std::print(myfile, "# is a fraction. What the number means here is the Henry coefficient above it.\n");
  }

  myfile.close();
}
