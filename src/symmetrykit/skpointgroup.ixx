module;

#ifdef USE_LEGACY_HEADERS
#include <string>
#include <vector>
#include <map>
#include <cstdint>
#include <optional>
#endif

#if defined(_WIN32)
  import <cassert>;
#else
  #include <assert.h>
#endif

export module skpointgroup;

#ifndef USE_LEGACY_HEADERS
import <string>;
import <vector>;
import <map>;
import <cstdint>;
import <optional>;
#endif

import skdefinitions;
import skrotationaloccurancetable;
import skpointsymmetryset;
import sktransformationmatrix;
import skrotationmatrix;

import int3;
import int3x3;
import double3;
import double3x3;

export enum class Laue : size_t
{
    none = 0, laue_1 = 1, laue_2m = 2, laue_mmm = 3, laue_4m = 4, laue_4mmm = 5, laue_3 = 6, laue_3m = 7, laue_6m = 8, laue_6mmm = 9, laue_m3 = 10, laue_m3m = 11
};


export class SKPointGroup
{
public:

    SKPointGroup(SKRotationalOccuranceTable table, size_t number, std::string symbol, std::string schoenflies, Holohedry holohedry, Laue laue, bool centrosymmetric, bool enantiomorphic);
    SKPointGroup(SKPointSymmetrySet set);
    static std::vector<SKPointGroup> pointGroupData;

    size_t number() { return _number; }
    Holohedry holohedry() const { return _holohedry; }
    std::string holohedryString() const;
    std::string LaueString() const;
    std::string symbol() { return _symbol; }
    std::string schoenflies() { return _schoenflies; }
    bool centrosymmetric() { return _centrosymmetric; }
    bool enantiomorphic() { return _enantiomorphic; }

    Laue laue() const { return _laue; }
    size_t number() const { return _number; }
    Centring computeCentering(SKTransformationMatrix basis);

    SKTransformationMatrix computeBasisCorrection(SKTransformationMatrix basis, Centring& centering);
    const std::optional<SKTransformationMatrix> constructAxes(std::vector<SKRotationMatrix> rotations) const;
    const SKRotationalOccuranceTable& table() const { return _table; }
private:
    SKRotationalOccuranceTable _table;
    size_t _number = 0;
    std::string _symbol = "";
    std::string _schoenflies = "";
    Holohedry _holohedry = Holohedry::none;
    Laue _laue = Laue::none;
    bool _centrosymmetric = false;
    bool _enantiomorphic = false;

    static std::map<Laue, size_t> rotationTypeForBasis;
};
