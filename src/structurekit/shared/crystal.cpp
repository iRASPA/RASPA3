module;

module crystal;

import std;

import double2;
import double3;

import pair_interactions;

std::vector<double3> Crystal::cartesianPositions() const
{
  std::vector<double3> positions;
  positions.reserve(this->atoms.size());

  for (const CrystalAtom& atom : this->atoms)
  {
    positions.push_back(atom.position);
  }

  return positions;
}

std::vector<double2> Crystal::lennardJonesParameters(const PairInteractions& interactions) const
{
  std::vector<double2> parameters;
  parameters.reserve(this->atoms.size());

  for (const CrystalAtom& atom : this->atoms)
  {
    const PairParameters& self = interactions[atom.type];
    parameters.push_back(double2(self.strengthParameter, self.sizeParameter));
  }

  return parameters;
}
