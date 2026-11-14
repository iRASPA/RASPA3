module;

export module minimization_torsion_hessian_geometry;

import std;

import double3;

export namespace Minimization
{
/**
 * Geometry and derivative factors of a dihedral A-B-C-D needed to assemble its Hessian.
 *
 * Follows RASPA2 conventions: DF = dU/dCosPhi, DDF = d²U/dCosPhi²,
 * Dcb, dr and ds are unit vectors, d = (Dab.Dcb)/rbc and e = (Ddc.Dcb)/rbc.
 */
struct TorsionHessianGeometry
{
  double DF{};
  double DDF{};
  double3 dtA{};
  double3 dtB{};
  double3 dtC{};
  double3 dtD{};
  double3 Dab{};
  double3 Dcb{};  // normalized
  double3 Ddc{};
  double rbc{};
  double3 dr{};  // normalized
  double3 ds{};  // normalized
  double r{};
  double s{};
  double d{};
  double e{};
  double cosPhi{};
  /// The Cartesian Hessian itself (coordinates A.x .. D.z), for a potential written in phi rather than cos(phi):
  /// the cos(phi) chain rule (DF, DDF) is singular at cis and trans. When set, it replaces the DF/DDF assembly.
  std::optional<std::array<std::array<double, 12>, 12>> cartesianHessian{};
};
}  // namespace Minimization
