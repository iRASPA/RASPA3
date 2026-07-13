module;

export module minimization_bend_hessian_geometry;

import double3;

export namespace Minimization
{
struct BendHessianGeometry
{
  double DF{};
  double DDF{};
  double3 dtA{};
  double3 dtB{};
  double3 dtC{};
  double3 Rab{};
  double3 Rbc{};
  double cosTheta{};
  double rab{};
  double rbc{};
};
}  // namespace Minimization
