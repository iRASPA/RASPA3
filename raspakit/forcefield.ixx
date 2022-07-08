export module forcefield;

import double4;
import int3;
import <vector>;
import <string>;
import <algorithm>;
import <iostream>;

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
    double reserved;
    double reserved2;
    Type type{ 0 };
    bool tailCorrection;
    bool shiftPotential;

    VDWParameters(): parameters(double4(158.5/ 1.2027242847,3.72,0.0,0.0)), shift(-0.56217796) {}
    VDWParameters(double epsilon, double sigma) : parameters(double4(epsilon, sigma, 0.0, 0.0)), type(Type::LennardJones), tailCorrection(false), shiftPotential(true)
    {
        double scaling = 1.0;
        double arg1 = epsilon;
        double arg2 = sigma * sigma;
        double rr = 12.0 * 12.0;
        double temp = (rr / arg2);
        double rri3 = 1.0 / ((temp * temp * temp) + 0.5 * (1.0 - scaling) * (1.0 - scaling));
        shift = scaling * (4.0 * arg1 * (rri3 * (rri3 - 1.0)));
    }
};

export struct PseudoAtom
{
    std::string name{ "C" };
    double mass{ 1.0 };
    double charge{ 0.0 };
    size_t atomicNumber{ 8 };
    bool printToPDB{ true };
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
    
    // 2D-vector, size numberOfPseudoAtoms squared
    std::vector<VDWParameters> data{};
    double cutOff{ 12.0 };
    double dualCutOff{ 5.0 };
    double cutOffCoulomb{ 12.0 };
    double alpha { 0.265058 };
    int3 numberOfWaveVectors{8, 8, 8};
    size_t numberOfPseudoAtoms{ 0 };
    std::vector< PseudoAtom> pseudoAtoms;

    ChargeMethod chargeMethod { ChargeMethod::Ewald};

    double overlapCriteria{ 1e5 };

    ForceField(std::string pseudoAtomsFileName, std::string forceFieldmixingFileName, std::string forceFieldOverwriteFileName) noexcept(false);
    
    ForceField(size_t size) : data(size* size)
    {
        std::fill(data.begin(), data.end(), VDWParameters(158.5/ 1.2027242847,3.72));
    }

    VDWParameters& operator() (size_t row, size_t col) { return data[row * numberOfPseudoAtoms + col]; }
    const VDWParameters&  operator() (size_t row, size_t col) const { return data[row * numberOfPseudoAtoms + col]; }

    void ReadPseudoAtoms(std::string pseudoAtomsFileName) noexcept(false);
    void ReadForceFieldMixing(std::string pseudoAtomsFileName) noexcept(false);

    std::string printPseudoAtomStatus();
    std::string printForceFieldStatus();
};
