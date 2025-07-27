module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <cstdlib>
#include <memory>
#include <numbers>
#include <tuple>
#include <vector>
#endif

export module skparser;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import double3;
import skstructure;

export class SKParser
{
 public:
  enum class ImportType : std::int64_t
  {
    asSeperateProjects = 0,
    asSingleProject = 1,
    asMovieFrames = 2
  };
  SKParser();
  virtual ~SKParser();
  virtual void startParsing() noexcept(false) = 0;
  std::vector<std::vector<std::shared_ptr<SKStructure>>> movies();

  std::vector<std::tuple<double3, std::size_t, double>> firstTestFrame();

 protected:
  double _a, _b, _c;
  double _alpha, _beta, _gamma;
  std::vector<std::vector<std::shared_ptr<SKStructure>>> _movies;
};
