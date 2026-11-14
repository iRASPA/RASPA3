module;

module interactions_external_field;

import std;

import int3;
import double3;
import double4;
import double3x3;
import atom;
import simulationbox;
import energy_status;
import energy_status_inter;
import units;
import potential_pair_derivatives;
import running_energy;
import component;
import forcefield;
import tricubic_derivatives_external_field;
import triquintic_derivatives_external_field;
import interpolation_energy_grid;

namespace
{
/// The external field at one atom. Like the framework interaction, the field couples to the fractional molecules
/// through the VDW scaling: U = scalingVDW * U_ext(r) for the smooth fields (analytic or tabulated), so that the
/// analytic evaluation of a field and its interpolation grid agree, and a switched-off molecule (scaling zero) has
/// no field energy at all, as the serial Gibbs CFCMC bookkeeping of inactive fractional molecules assumes.
/// 'energy' is that scaled energy, 'dUdlambda' the derivative factor dU/d(scalingVDW) = U_ext (see
/// PairDerivatives), and 'inside' is false when a coupled atom lies outside a confining geometry (cylinder /
/// rectangular channel through the cell centre). The callers turn 'outside' into either an overlap rejection
/// (trial positions, like the CBMC growth) or the overlap sentinel energy 10 * energyOverlapCriteria (the total
/// energy and positions that are already part of the configuration), so that a molecule that starts outside the
/// channel is not trapped there: a move that brings it inside lowers the energy by the sentinel and is accepted.
struct ExternalFieldValue
{
  double energy{0.0};
  double dUdlambda{0.0};
  bool inside{true};
};

// Whether 'posA' lies within 'radius' of the axis through the cell centre along 'axisBegin -> axisEnd'.
bool insideCylinder(const double3& posA, const double3& axisBegin, const double3& axisEnd, double radius)
{
  const double3 v = (axisEnd - axisBegin).normalized();
  const double3 w = double3::cross(posA - axisBegin, v);
  return w.length_squared() < radius * radius;
}

// Whether 'posA' lies within the rectangular channel along 'axisBegin -> axisEnd'. 'first' and 'second' index
// the two transverse components of w = (posA - axisBegin) x v that are compared against the channel width
// (externalFieldGeometryParameters.x) and height (externalFieldGeometryParameters.y). The same cyclic choice of
// components is made for the three axes: (y, z) for a channel along x, (z, x) along y, and (x, y) along z, so
// that the width is the extent along z, x, and y, and the height the extent along y, z, and x, respectively.
bool insideRectangle(const double3& posA, const double3& axisBegin, const double3& axisEnd, std::size_t first,
                     std::size_t second, const double4& externalFieldGeometryParameters)
{
  const double3 v = (axisEnd - axisBegin).normalized();
  const double3 w = double3::cross(posA - axisBegin, v);
  return std::abs(w[first]) < externalFieldGeometryParameters.x &&
         std::abs(w[second]) < externalFieldGeometryParameters.y;
}

ExternalFieldValue externalFieldValue(const ForceField& forceField, const SimulationBox& simulationBox,
                                      const std::optional<InterpolationEnergyGrid>& externalFieldInterpolationGrid,
                                      const Atom& atom)
{
  const double3 posA = atom.position;
  const double scalingVDWA = atom.scalingVDW;
  const double4 externalFieldGeometryParameters = forceField.externalFieldGeometryParameters;

  ExternalFieldValue value{};

  // a decoupled atom (an inactive fractional molecule) feels neither the field nor the confining geometry
  if (scalingVDWA == 0.0) return value;

  if (externalFieldInterpolationGrid.has_value())
  {
    value.dUdlambda = externalFieldInterpolationGrid->interpolate(posA);
    value.energy = scalingVDWA * value.dUdlambda;
    return value;
  }

  const double fieldEnergy = Interactions::analyticExternalFieldEnergy(forceField, posA);
  // the confining geometry is periodic: test the position wrapped into the unit cell, so that a molecule that
  // drifted out of the cell, or a bead of a chain overhanging its edge, is tested against its own image of the
  // channel (as the interpolation grid does) and not against the channel of the home cell
  const double3 wrappedA = simulationBox.cell * double3::fract(simulationBox.inverseCell * posA);
  switch (forceField.potentialEnergySurfaceType)
  {
    case ForceField::PotentialEnergySurfaceType::None:
      break;
    case ForceField::PotentialEnergySurfaceType::GridFile:
      // A field read from file is only available through the interpolation grid (checked by the force-field
      // reader: 'GridFile' requires 'UseExternalFieldGrid').
      break;
    case ForceField::PotentialEnergySurfaceType::CylinderX:
      value.inside = insideCylinder(wrappedA, simulationBox.cell * double3{0.0, 0.5, 0.5},
                                    simulationBox.cell * double3{1.0, 0.5, 0.5}, externalFieldGeometryParameters.x);
      break;
    case ForceField::PotentialEnergySurfaceType::CylinderY:
      value.inside = insideCylinder(wrappedA, simulationBox.cell * double3{0.5, 0.0, 0.5},
                                    simulationBox.cell * double3{0.5, 1.0, 0.5}, externalFieldGeometryParameters.x);
      break;
    case ForceField::PotentialEnergySurfaceType::CylinderZ:
      value.inside = insideCylinder(wrappedA, simulationBox.cell * double3{0.5, 0.5, 0.0},
                                    simulationBox.cell * double3{0.5, 0.5, 1.0}, externalFieldGeometryParameters.x);
      break;
    case ForceField::PotentialEnergySurfaceType::RectangleX:
      value.inside =
          insideRectangle(wrappedA, simulationBox.cell * double3{0.0, 0.5, 0.5},
                          simulationBox.cell * double3{1.0, 0.5, 0.5}, 1, 2, externalFieldGeometryParameters);
      break;
    case ForceField::PotentialEnergySurfaceType::RectangleY:
      value.inside =
          insideRectangle(wrappedA, simulationBox.cell * double3{0.5, 0.0, 0.5},
                          simulationBox.cell * double3{0.5, 1.0, 0.5}, 2, 0, externalFieldGeometryParameters);
      break;
    case ForceField::PotentialEnergySurfaceType::RectangleZ:
      value.inside =
          insideRectangle(wrappedA, simulationBox.cell * double3{0.5, 0.5, 0.0},
                          simulationBox.cell * double3{0.5, 0.5, 1.0}, 0, 1, externalFieldGeometryParameters);
      break;
    default:
      break;
  }

  if (!value.inside)
  {
    // the hard wall of a confining geometry is not softened with lambda
    value.energy = 10.0 * forceField.energyOverlapCriteria;
    return value;
  }

  value.dUdlambda = fieldEnergy;
  value.energy = scalingVDWA * fieldEnergy;
  return value;
}
}  // namespace

double Interactions::analyticExternalFieldEnergy(const ForceField& forceField, const double3& position)
{
  const double3 pos = position + forceField.potentialEnergySurfaceOrigin;

  switch (forceField.potentialEnergySurfaceType)
  {
    case ForceField::PotentialEnergySurfaceType::SecondOrderPolynomialTestFunction:
      return std::get<0>(Potentials::potentialTricubicSecondOrderPolynomialTestFunction(pos));
    case ForceField::PotentialEnergySurfaceType::ThirdOrderPolynomialTestFunction:
      return std::get<0>(Potentials::potentialTricubicThirdOrderPolynomialTestFunction(pos));
    case ForceField::PotentialEnergySurfaceType::FourthOrderPolynomialTestFunction:
      return std::get<0>(Potentials::potentialTricubicFourthOrderPolynomialTestFunction(pos));
    case ForceField::PotentialEnergySurfaceType::FifthOrderPolynomialTestFunction:
      return std::get<0>(Potentials::potentialTricubicFifthOrderPolynomialTestFunction(pos));
    case ForceField::PotentialEnergySurfaceType::SixthOrderPolynomialTestFunction:
      return std::get<0>(Potentials::potentialTricubicSixthOrderPolynomialTestFunction(pos));
    case ForceField::PotentialEnergySurfaceType::ExponentialNonPolynomialTestFunction:
      return std::get<0>(Potentials::potentialTricubicExponentialNonPolynomialTestFunction(pos));
    case ForceField::PotentialEnergySurfaceType::MullerBrown:
    {
      const double4 A{-200.0, -100.0, -170.0, 15.0};
      const double4 a{-1.0, -1.0, -6.5, 0.7};
      const double4 x0{1.0, 0.0, -0.5, -1.0};
      const double4 b{0.0, 0.0, 11.0, 0.6};
      const double4 y0{0.0, 0.5, 1.5, 1.0};
      const double4 c{-10.0, -10.0, -6.5, 0.7};
      double energy{0.0};
      for (std::size_t i = 0; i < 4; ++i)
      {
        energy += A[i] * std::exp(a[i] * (pos.x - x0[i]) * (pos.x - x0[i]) + b[i] * (pos.x - x0[i]) * (pos.y - y0[i]) +
                                  c[i] * (pos.y - y0[i]) * (pos.y - y0[i]));
      }
      return energy;
    }
    default:
      return 0.0;
  }
}

void Interactions::computeExternalFieldEnergy(
    bool hasExternalField, [[maybe_unused]] const ForceField& forceField,
    [[maybe_unused]] const SimulationBox& simulationBox, [[maybe_unused]] std::span<const Atom> moleculeAtoms,
    [[maybe_unused]] RunningEnergy& energyStatus,
    [[maybe_unused]] const std::optional<InterpolationEnergyGrid>& externalFieldInterpolationGrid) noexcept
{
  if (hasExternalField)
  {
    for (const Atom& atom : moleculeAtoms)
    {
      const ExternalFieldValue value =
          externalFieldValue(forceField, simulationBox, externalFieldInterpolationGrid, atom);

      energyStatus.externalFieldVDW += value.energy;
      energyStatus.addDudlambdaVDW(atom.groupId, 0, 1.0, 1.0, value.dUdlambda);
    }
  }
}

void Interactions::computeExternalFieldTailEnergy(bool hasExternalField, [[maybe_unused]] const ForceField& forceField,
                                                  [[maybe_unused]] const SimulationBox& simulationBox,
                                                  [[maybe_unused]] std::span<const Atom> moleculeAtoms,
                                                  [[maybe_unused]] RunningEnergy& energyStatus) noexcept
{
  if (hasExternalField)
  {
  }
}

[[nodiscard]] std::optional<RunningEnergy> Interactions::computeExternalFieldEnergyDifference(
    bool hasExternalField, [[maybe_unused]] const ForceField& forceField,
    [[maybe_unused]] const SimulationBox& simulationBox,
    [[maybe_unused]] const std::optional<InterpolationEnergyGrid>& externalFieldInterpolationGrid,
    [[maybe_unused]] std::span<const Atom> newatoms, [[maybe_unused]] std::span<const Atom> oldatoms) noexcept
{
  RunningEnergy energySum;

  if (hasExternalField)
  {
    // Trial positions outside a confining geometry are an overlap (as in the CBMC growth).
    for (const Atom& atom : newatoms)
    {
      const ExternalFieldValue value =
          externalFieldValue(forceField, simulationBox, externalFieldInterpolationGrid, atom);
      if (!value.inside) return std::nullopt;

      energySum.externalFieldVDW += value.energy;
      energySum.addDudlambdaVDW(atom.groupId, 0, 1.0, 1.0, value.dUdlambda);
    }

    // The old positions are part of the configuration: outside the geometry they carry the sentinel energy that
    // the total energy assigned to them, so the difference stays consistent with the running total and a molecule
    // that starts outside can be moved back in.
    for (const Atom& atom : oldatoms)
    {
      const ExternalFieldValue value =
          externalFieldValue(forceField, simulationBox, externalFieldInterpolationGrid, atom);

      energySum.externalFieldVDW -= value.energy;
      energySum.addDudlambdaVDW(atom.groupId, 0, 1.0, 1.0, -value.dUdlambda);
    }
  }

  return energySum;
}
