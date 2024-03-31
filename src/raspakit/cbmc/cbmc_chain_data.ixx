export module cbmc_chain_data;

import <vector>;

import atom;
import double3x3;
import double3;
import randomnumbers;
import running_energy;


export struct ChainData
{
  std::vector<Atom> atom;
  RunningEnergy energies;
  double RosenbluthWeight;
  double storedR;

  ChainData(std::vector<Atom> atom, RunningEnergy energies, double RosenbluthWeight, double storedR) noexcept :
      atom(atom), energies(energies), RosenbluthWeight(RosenbluthWeight), storedR(storedR)
  {
  }
  ChainData() = delete;
  ChainData(const ChainData& a) noexcept = default;
  ChainData& operator=(const ChainData& a) noexcept = default;
  ChainData(ChainData&& a) noexcept = default;
  ChainData& operator=(ChainData&& a) noexcept = default;
  ~ChainData() noexcept = default;
};
