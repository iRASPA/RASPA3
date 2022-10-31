module;

module skparser;

import <cstdlib>;
import <vector>;
import <tuple>;
import <memory>;
import <numbers>;
import double3;

import skstructure;
import skasymmetricatom;

SKParser::SKParser() : _a(20.0), _b(20.0), _c(20.0), _alpha(90.0 * std::numbers::pi / 180.0), _beta(90.0 * std::numbers::pi / 180.0), _gamma(90.0 * std::numbers::pi / 180.0)
{

}

SKParser::~SKParser()
{
    // Compulsory virtual destructor definition
}


std::vector<std::vector<std::shared_ptr<SKStructure>>> SKParser::movies()
{
    return _movies;
}

std::vector<std::tuple<double3, int, double> > SKParser::firstTestFrame()
{
    std::vector<std::tuple<double3, int, double> > atoms{};

    for (const std::vector<std::shared_ptr<SKStructure>>& movie : _movies)
    {
        for (const std::shared_ptr<SKStructure>& structure : movie)
        {
            for (const std::shared_ptr<SKAsymmetricAtom>& atom : structure->atoms)
            {
                std::tuple<double3, int, double> atomTuple = std::make_tuple<double3, int, double>(atom->position(), atom->elementIdentifier(), 1.0);
                atoms.push_back(atomTuple);
            }
        }
    }
    return atoms;
}

