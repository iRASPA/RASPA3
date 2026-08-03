module;

module energy_shared_linear_probe;

import std;

import double3;
import pair_interactions;

namespace
{
// The bond lengths of the built-in molecules. The rest of each molecule, the sizes, strengths and charges,
// is taken from the force field's own pseudo-atoms, so that a probe here is the same molecule a simulation
// with this force field would use and no parameter is written down twice.
struct BuiltIn
{
  std::string name;
  std::vector<std::pair<std::string, double>> sites;  // pseudo-atom name and its offset along the axis [Å]
  bool headTailSymmetric;
};

const std::vector<BuiltIn> &builtIns()
{
  // The geometries are the ones this project's own molecule definitions use, so that a probe here is the
  // molecule a simulation with the same force field would insert. Only the shape is written down; every
  // size, strength and charge is read from the pseudo-atoms at run time.
  //
  // Nitrogen and oxygen carry a massless site at the centre holding a compensating positive charge, which
  // is how a quadrupole is put on a homonuclear diatomic without giving it a dipole it does not have. That
  // site has no dispersion, and the force fields say so by giving it an interaction of type `none`; it
  // costs nothing here beyond a term in the electrostatic sum, and leaving it out would throw away the
  // whole of the quadrupole, which is most of what nitrogen does in a zeolite.
  static const std::vector<BuiltIn> table{
      // The two oxygens 1.149 Å either side of the carbon, which is the geometry that goes with the C_co2
      // and O_co2 parameters the force fields carry.
      {"CO2", {{"O_co2", -1.149}, {"C_co2", 0.0}, {"O_co2", 1.149}}, true},

      // A bond of 1.10 Å, as in this project's own n2.json.
      {"N2", {{"N_n2", -0.550}, {"N_com", 0.0}, {"N_n2", 0.550}}, true},

      // A bond of 1.21 Å, which is the measured one and is what both of the Martin-Calvo parameterisations
      // of oxygen are built on, so the shape does not turn on which of them a force field carries.
      {"O2", {{"O_o2", -0.605}, {"O_com", 0.0}, {"O_o2", 0.605}}, true},

      // Hydrogen is put together the other way about from the two above: Darkrim and Levesque put the
      // dispersion at the centre and nothing but charge on the two protons, 0.371 Å out. Nothing here has
      // to know that, the sites being read from the force field either way, but a force field carrying the
      // united-atom hydrogen instead will not have H_com and the molecule will not be built, which is the
      // right outcome: the two are different models and neither is a substitute for the other.
      {"H2", {{"H_h2", -0.371}, {"H_com", 0.0}, {"H_h2", 0.371}}, true},

      // Nitric oxide is the one here whose two ends are not alike, so the centre of mass is not the middle
      // of the bond: over 1.15 Å the heavier oxygen pulls it to 0.613 Å from the nitrogen. Turning the
      // molecule end for end is a real change to it, and `false` says so, which costs a full sphere of
      // orientations instead of a half.
      {"NO", {{"N_no", -0.6132}, {"O_no", 0.5368}}, false}};
  return table;
}
}  // namespace


std::optional<LinearProbe> LinearProbe::named(const PairInteractions &interactions, const std::string &name)
{
  for (const BuiltIn &candidate : builtIns())
  {
    if (candidate.name != name) continue;

    LinearProbe probe;
    probe.name = candidate.name;
    probe.headTailSymmetric = candidate.headTailSymmetric;
    for (const auto &[siteName, offset] : candidate.sites)
    {
      std::optional<std::size_t> type = interactions.findType(siteName);
      if (!type.has_value()) return std::nullopt;

      probe.sites.push_back(
          LinearProbe::Site{type.value(), offset, interactions.charges[type.value()], siteName});
    }
    return probe;
  }
  return std::nullopt;
}


std::optional<LinearProbe> LinearProbe::singleSite(const PairInteractions &interactions, const std::string &pseudoAtomName)
{
  std::optional<std::size_t> type = interactions.findType(pseudoAtomName);
  if (!type.has_value()) return std::nullopt;

  LinearProbe probe;
  probe.name = pseudoAtomName;
  probe.headTailSymmetric = true;
  probe.sites.push_back(
      LinearProbe::Site{type.value(), 0.0, interactions.charges[type.value()], pseudoAtomName});
  return probe;
}


std::vector<std::string> LinearProbe::builtInNames()
{
  std::vector<std::string> names;
  for (const BuiltIn &candidate : builtIns()) names.push_back(candidate.name);
  return names;
}


std::vector<double3> orientationSet(std::size_t count, bool overHemisphere)
{
  std::vector<double3> directions;
  if (count == 0) return directions;
  directions.reserve(count);

  // Stepping the height in equal intervals makes the points equal-area, since the area of a zone of the
  // sphere depends only on the height it spans. Turning by the golden angle each step keeps them from
  // lining up into the spiral arms an ordinary angle would produce.
  const double goldenAngle = 2.0 * std::numbers::pi / std::numbers::phi;
  const double span = overHemisphere ? 1.0 : 2.0;

  for (std::size_t i = 0; i < count; ++i)
  {
    double height = 1.0 - span * (static_cast<double>(i) + 0.5) / static_cast<double>(count);
    double radius = std::sqrt(std::max(0.0, 1.0 - height * height));
    double angle = goldenAngle * static_cast<double>(i);

    directions.push_back(double3(radius * std::cos(angle), radius * std::sin(angle), height));
  }

  return directions;
}
