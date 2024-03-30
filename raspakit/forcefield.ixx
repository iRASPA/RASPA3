export module forcefield;

import <vector>;
import <string>;
import <algorithm>;
import <iostream>;
import <ostream>;
import <fstream>;
import <optional>;

import archive;
import double4;
import double3;
import int3;


export struct VDWParameters
{
  enum class Type : int
  {
    LennardJones = 0,
    BuckingHam = 1,
    Morse = 2,
    FeynmannHibbs = 3,
    MM3 = 4,
    BornHugginsMeyer = 5
  };

  double4 parameters;      // for LJ: epsilon, sigma, for Buckingham: 3 parameters
  double shift;
  double tailCorrectionEnergy;
  Type type{ 0 };

  VDWParameters(): parameters(double4(158.5/ 1.2027242847,3.72,0.0,0.0)), shift(-0.56217796) {}

  VDWParameters(double epsilon, double sigma) : parameters(double4(epsilon, sigma, 0.0, 0.0)), shift(0.0), 
                                                tailCorrectionEnergy(0.0), type(Type::LennardJones)
  {
  }

  bool operator==(VDWParameters const&) const = default;

  void computeShiftAtCutOff(double cutOff)
  {
    shift = 0.0;
    double scaling = 1.0;
    double arg1 = parameters.x;
    double arg2 = parameters.y * parameters.y;
    double rr = cutOff * cutOff;
    double temp = (rr / arg2);
    double rri3 = 1.0 / ((temp * temp * temp) + 0.5 * (1.0 - scaling) * (1.0 - scaling));
    shift = scaling * (4.0 * arg1 * (rri3 * (rri3 - 1.0)));
  }

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const VDWParameters &p);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, VDWParameters &p);
};

export struct PseudoAtom
{
  PseudoAtom() {};
  PseudoAtom(std::string name, double mass, double charge, size_t atomicNumber, bool printToPDB):
             name(name), mass(mass), charge(charge), atomicNumber(atomicNumber), printToPDB(printToPDB) {};
  uint64_t versionNumber{ 1 };
  std::string name{ "C" };
  double mass{ 1.0 };
  double charge{ 0.0 };
  size_t atomicNumber{ 8 };
  bool printToPDB{ true };

  bool operator==(PseudoAtom const&) const = default;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PseudoAtom &a);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PseudoAtom &a);
};

export struct ForceField
{
  enum class ChargeMethod : int
  {
    Ewald = 0,
    Coulomb = 1,
    Wolf = 2,
    ModifiedWolf = 3
  };

  enum class MixingRule : int
  {
      Lorentz_Berthelot = 0
  };
  
  uint64_t versionNumber{ 1 };

  // 2D-vector, size numberOfPseudoAtoms squared
  std::vector<VDWParameters> data{};
  std::vector<bool> shiftPotentials{};
  std::vector<bool> tailCorrections{};
  double cutOffVDW{ 12.0 };
  double cutOffCoulomb{ 12.0 };
  double dualCutOff{ 6.0 };
  
  size_t numberOfPseudoAtoms{ 0 };
  std::vector<PseudoAtom> pseudoAtoms{};

  ChargeMethod chargeMethod { ChargeMethod::Ewald};

  double overlapCriteria{ 1e5 };

  double EwaldPrecision{ 1e-6 };
  double EwaldAlpha{ 0.265058 };
  int3 numberOfWaveVectors{ 8, 8, 8 };
  bool automaticEwald{ true };

  bool noCharges{ false };
  bool omitEwaldFourier{ false };

  double minimumRosenbluthFactor{ 1e-150 };
  double energyOverlapCriteria = 1e6;
  bool useDualCutOff{ false };

  ForceField(std::vector<PseudoAtom> pseudoAtoms, std::vector<VDWParameters> parameters, MixingRule mixingRule, 
             double cutOff, bool shifted, bool tailCorrecions) noexcept(false);
  ForceField(size_t systemId) noexcept(false);

  ForceField() noexcept = default;
  ForceField(const ForceField &a) noexcept = default;
  ForceField& operator=(const ForceField& a) noexcept = default;
  ForceField(ForceField&& a) noexcept = default;
  ForceField& operator=(ForceField&& a) noexcept = default;
  ~ForceField() noexcept = default;

  VDWParameters& operator() (size_t row, size_t col) { return data[row * numberOfPseudoAtoms + col]; }
  const VDWParameters&  operator() (size_t row, size_t col) const { return data[row * numberOfPseudoAtoms + col]; }

  void preComputeTailCorrection();

  void ReadPseudoAtoms(std::string pseudoAtomsFileName) noexcept(false);
  void ReadForceFieldMixing(std::string pseudoAtomsFileName) noexcept(false);

  std::string printPseudoAtomStatus() const;
  std::string printForceFieldStatus() const;

  std::optional<size_t> findPseudoAtom(const std::string &name) const;

  void initializeEwaldParameters(double3 perpendicularWidths);

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const ForceField &f);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, ForceField &f);

  std::string repr() const {return printPseudoAtomStatus() + "\n" + printForceFieldStatus();}
};
