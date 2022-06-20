export module skspacegroup;

import <vector>;
import <tuple>;
import <array>;
import <string>;
import <optional>;
import <unordered_set>;
import mathkit;
import skseitzmatrix;
import skpointgroup;
import skpointsymmetryset;
import sksymmetryoperationset;
import skspacegroupsetting;
import skrotationalchangeofbasis;
import skspacegroupdatabase;
import sksymmetrycell;

export class SKSpaceGroup
{
public:
    struct FoundSpaceGroupInfo
    {
        int HallNumber;
        double3 origin;
        SKSymmetryCell cell;
        SKRotationalChangeOfBasis changeOfBasis;
        double3x3 transformationMatrix;
        double3x3 rotationMatrix;
        std::vector<std::tuple<double3, int, double>> atoms;
        std::vector<std::tuple<double3, int, double>> asymmetricAtoms;
    };

    struct FoundNiggliCellInfo
    {
        int HallNumber;
        SKSymmetryCell cell;
        std::vector<std::tuple<double3, int, double>> atoms;
    };

    struct FoundPrimitiveCellInfo
    {
        SKSymmetryCell cell;
        std::vector<std::tuple<double3, int, double>> atoms;
    };

    SKSpaceGroup(int HallNumber);
    std::vector<double3> listOfSymmetricPositions(double3 pos);
    const SKSpaceGroupSetting& spaceGroupSetting() const { return _spaceGroupSetting; }

    static std::vector<std::string> latticeTranslationStrings(int HallNumber);
    static std::string inversionCenterString(int HallNumber);
    static std::optional<int> HallNumberFromHMString(std::string inputString);
    static std::optional<int> HallNumberFromSpaceGroupNumber(int);
    static std::optional<int> HallNumber(std::string inputString);
    static std::optional<FoundPrimitiveCellInfo> SKFindPrimitive(double3x3 unitCell, std::vector<std::tuple<double3, int, double>> atoms, bool allowPartialOccupancies, double symmetryPrecision);
    static std::optional<FoundNiggliCellInfo> findNiggliCell(double3x3 unitCell, std::vector<std::tuple<double3, int, double> > atoms, bool allowPartialOccupancies, double symmetryPrecision);
    static std::optional<FoundSpaceGroupInfo> findSpaceGroup(double3x3 unitCell, std::vector<std::tuple<double3, int, double> > atoms, bool allowPartialOccupancies, double symmetryPrecision);

    static SKSymmetryOperationSet findSpaceGroupSymmetry(union double3x3 unitCell, std::vector<std::tuple<double3, int, double>> reducedAtoms, std::vector<std::tuple<double3, int, double>> atoms, SKPointSymmetrySet latticeSymmetries, bool allowPartialOccupancies, double symmetryPrecision);
    static std::optional<std::pair<double3, SKRotationalChangeOfBasis>> matchSpaceGroup(int HallNumber, double3x3 lattice, Centring entering, std::vector<SKSeitzMatrix> seitzMatrices, double symmetryPrecision);
    static std::optional<double3> getOriginShift(int HallNumber, Centring centering, SKRotationalChangeOfBasis changeOfBasis, std::vector<SKSeitzMatrix> seitzMatrices, double symmetryPrecision);
private:
    SKSpaceGroupSetting _spaceGroupSetting = SKSpaceGroupDataBase::spaceGroupData[1];

    static bool matchSpacegroup(std::string spaceSearchGroupString, std::string storedSpaceGroupString);
};