export module cbmc_first_bead_data;

import <vector>;

import atom;
import double3x3;
import double3;
import randomnumbers;
import running_energy;


export struct FirstBeadData
{
  Atom atom;
  RunningEnergy energies;
  double RosenbluthWeight;
  double storedR;

  FirstBeadData(Atom atom, RunningEnergy energies, double RosenbluthWeight, double storedR) noexcept :
      atom(atom), energies(energies), RosenbluthWeight(RosenbluthWeight), storedR(storedR)
  {
  }

  FirstBeadData() noexcept = delete;
  FirstBeadData(const FirstBeadData& a) noexcept = default;
  FirstBeadData& operator=(const FirstBeadData& a) noexcept = default;
  FirstBeadData(FirstBeadData&& a) noexcept = default;
  FirstBeadData& operator=(FirstBeadData&& a) noexcept = default;
  ~FirstBeadData() noexcept = default;
};

