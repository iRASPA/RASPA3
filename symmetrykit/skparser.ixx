export module skparser;

import <cstdlib>;
import <vector>;
import <tuple>;
import <memory>;
import <numbers>;
import double3;
import skstructure;

export class SKParser
{
public:
    enum class ImportType : int64_t
    {
        asSeperateProjects = 0,
        asSingleProject = 1,
        asMovieFrames = 2
    };
    SKParser();
    virtual ~SKParser();
    virtual void startParsing() noexcept(false) = 0;
    std::vector<std::vector<std::shared_ptr<SKStructure>>> movies();

    std::vector<std::tuple<double3, int, double> > firstTestFrame();
protected:
    double _a, _b, _c;
    double _alpha, _beta, _gamma;
    std::vector<std::vector<std::shared_ptr<SKStructure>>> _movies;
};
