module;

module energy_shared_molecular_energy_grid;

import std;

import double3;
import unit_cell;

MolecularEnergyGrid::MolecularEnergyGrid() {}

MolecularEnergyGrid::~MolecularEnergyGrid() {}

double3 MolecularEnergyGrid::fractionalPosition(std::size_t i, std::size_t j, std::size_t k) const
{
  return double3(static_cast<double>(i) / static_cast<double>(this->gridSize.x),
                 static_cast<double>(j) / static_cast<double>(this->gridSize.y),
                 static_cast<double>(k) / static_cast<double>(this->gridSize.z));
}

double3 MolecularEnergyGrid::spacing() const
{
  double3 a = this->unitCell.cell[0];
  double3 b = this->unitCell.cell[1];
  double3 c = this->unitCell.cell[2];

  return double3(a.length() / static_cast<double>(this->gridSize.x),
                 b.length() / static_cast<double>(this->gridSize.y),
                 c.length() / static_cast<double>(this->gridSize.z));
}

double MolecularEnergyGrid::deepestFreeEnergy() const
{
  if (this->freeEnergy.empty()) return 0.0;
  return static_cast<double>(*std::ranges::min_element(this->freeEnergy));
}

double MolecularEnergyGrid::deepestMinimumEnergy() const
{
  if (this->minimumEnergy.empty()) return 0.0;
  return static_cast<double>(*std::ranges::min_element(this->minimumEnergy));
}
