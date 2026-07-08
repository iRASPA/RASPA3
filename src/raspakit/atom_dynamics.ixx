module;

export module atom_dynamics;

import std;

import archive;
import double3;
import json;

/**
 * \brief Cold per-atom molecular-dynamics fields, stored parallel to Atom.
 *
 * The velocity and gradient (force) of an atom are not read in the energy inner
 * loops, so they are split out of Atom into this separate array to keep Atom a
 * single 64-byte cache line. An AtomDynamics record is stored at the same index
 * as its Atom in the parallel storage; see System::atomData / System::atomDynamics.
 * AtomDynamics is exactly 64 bytes (2 times double4).
 */
export struct AtomDynamics
{
  double3 velocity{};  ///< The velocity of the atom.
  double3 gradient{};  ///< The gradient (force) acting on the atom.

  AtomDynamics() noexcept = default;
  AtomDynamics(double3 velocity, double3 gradient) noexcept : velocity(velocity), gradient(gradient) {}

  bool operator==(AtomDynamics const&) const = default;

  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const AtomDynamics& atom);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, AtomDynamics& atom);
  friend void to_json(nlohmann::json&, const AtomDynamics&);
  friend void from_json(const nlohmann::json&, AtomDynamics&);

  inline std::string repr() const
  {
    std::ostringstream stream;
    std::print(stream, "velocity: ({}, {}, {}), gradient: ({}, {}, {})\n", velocity.x, velocity.y, velocity.z,
               gradient.x, gradient.y, gradient.z);
    return stream.str();
  }
};

// velocity (double4) + gradient (double4) = 2x32 = 64 bytes
static_assert(sizeof(AtomDynamics) == 64, "struct AtomDynamics size is not 64");

export void to_json(nlohmann::json& j, const AtomDynamics& a)
{
  j = nlohmann::json{{"velocity", a.velocity}, {"gradient", a.gradient}};
}

export void from_json(const nlohmann::json& j, AtomDynamics& a)
{
  j.at("velocity").get_to(a.velocity);
  j.at("gradient").get_to(a.gradient);
}
