module;

export module energy_shared_linear_probe;

import std;

import double3;
import pair_interactions;

// A rigid linear probe molecule: a few interaction sites strung along one axis, at fixed distances from the
// centre of mass. Carbon dioxide and nitrogen are of this kind, and so is any diatomic.
//
// The single-site probes the other routes use are the degenerate case, one site at the centre, and they go
// through the same code. That is worth having: a one-site probe carried through the orientational machinery
// has to give back exactly what the single-site route gives, which is the sharpest test available of the
// machinery being right.
export struct LinearProbe
{
  struct Site
  {
    // The pseudo-atom this site is, so that the force field's own mixing rule applies to it.
    std::size_t type{0};

    // Where it sits along the axis, in Å, measured from the centre of mass. Negative is the other end.
    double offset{0.0};

    // Its partial charge, in units of the elementary charge.
    double charge{0.0};

    std::string name;
  };

  std::string name;
  std::vector<Site> sites;

  // Whether turning the molecule end for end leaves it unchanged, as it does for carbon dioxide and
  // nitrogen. When it does, orientations come in pairs that give the same energy, and only half the sphere
  // need be sampled.
  bool headTailSymmetric{true};

  bool isCharged() const
  {
    return std::ranges::any_of(this->sites, [](const Site &site) { return site.charge != 0.0; });
  }

  double length() const
  {
    if (this->sites.empty()) return 0.0;
    auto [low, high] = std::ranges::minmax(this->sites | std::views::transform(&Site::offset));
    return high - low;
  }

  // The probes that can be named on the command line. A single-site probe of any pseudo-atom is made on
  // demand, so only the ones with a shape of their own are listed here.
  static std::optional<LinearProbe> named(const PairInteractions &interactions, const std::string &name);

  // One site at the centre, which is what the other routes mean by a probe.
  static std::optional<LinearProbe> singleSite(const PairInteractions &interactions, const std::string &pseudoAtomName);

  static std::vector<std::string> builtInNames();
};

// A set of directions the molecule's axis may point along, spread as evenly over the sphere as a set of that
// size can be. `overHemisphere` halves the sphere for a molecule that is the same end for end, which is not
// an approximation but a removal of exact duplicates.
//
// The points are the spherical Fibonacci set: taking the height in equal steps makes them equal-area, so a
// plain average over them is the average over the sphere, with no weights to carry.
export std::vector<double3> orientationSet(std::size_t count, bool overHemisphere);
