module;

module interactions_hessian_intramolecular;

import std;

import double3;
import double3x3;
import int3;
import units;
import forcefield;
import framework;
import simulationbox;
import intra_molecular_potentials;
import bond_potential;
import bend_potential;
import urey_bradley_potential;
import van_der_waals_potential;
import coulomb_potential;
import torsion_potential;
import component;
import bend_potential_gradient_hessian_strain;
import torsion_potential_gradient_hessian_strain;
import distance_potential_gradient_hessian_strain;
import minimization_hessian_scatter;
import minimization_rigid_kinematics;
import minimization_cell_layout;
import potential_pair_derivatives;
import potential_pair_vdw;
import potential_pair_coulomb;

namespace
{
// UreyBradleyType enumerators are the BondType enumerators shifted by one (no 'None' entry).
BondType bondTypeFromUreyBradley(UreyBradleyType type)
{
  return static_cast<BondType>(static_cast<std::size_t>(type) + 1);
}

// Under the molecular (rigid-body) strain convention, atoms inside a rigid group move with the
// group's center of mass: their internal offset (position - group CoM) does not scale with the
// cell. The intramolecular kernels accumulate the atomic virial (arm = atom position), so for
// semi-flexible molecules the non-scaling offset contribution g_i (x) (pos_i - com_i) of every
// rigid-group atom of the term has to be removed from the strain gradient.
template <std::size_t N>
void removeRigidOffsetStrainGradient(GeneralizedHessian& hessian, const MinimizationDofLayout& layout,
                                     const Minimization::RigidDerivativeCache& rigidCache, std::size_t moleculeIndex,
                                     std::span<const Atom> moleculeAtoms, const std::array<std::size_t, N>& atomIds,
                                     const std::array<double3, N>& gradients)
{
  for (std::size_t i = 0; i < N; ++i)
  {
    const std::size_t localAtom = atomIds[i];
    if (!layout.atomRigidComDof(moleculeIndex, localAtom).has_value())
    {
      continue;
    }
    const double3 offset = moleculeAtoms[localAtom].position - rigidCache.bodyCenterOfMass(moleculeIndex, localAtom);
    const double3& gradient = gradients[i];
    hessian.strainGradient().ax -= offset.x * gradient.x;
    hessian.strainGradient().bx -= offset.y * gradient.x;
    hessian.strainGradient().cx -= offset.z * gradient.x;
    hessian.strainGradient().ay -= offset.x * gradient.y;
    hessian.strainGradient().by -= offset.y * gradient.y;
    hessian.strainGradient().cy -= offset.z * gradient.y;
    hessian.strainGradient().az -= offset.x * gradient.z;
    hessian.strainGradient().bz -= offset.y * gradient.z;
    hessian.strainGradient().cz -= offset.z * gradient.z;
  }
}

template <std::size_t N>
bool hasPeriodicShift(const std::array<int3, N>& shifts)
{
  return std::ranges::any_of(shifts, [](const int3& shift) { return shift.x != 0 || shift.y != 0 || shift.z != 0; });
}

template <std::size_t N>
void addFrameworkTermCellCorrection(GeneralizedHessian& hessian, const MinimizationDofLayout& layout,
                                    const CellMinimizationLayout& cellLayout, const std::array<std::size_t, N>& atomIds,
                                    const std::array<double3, N>& positions,
                                    const std::array<double3, N>& shiftedPositions,
                                    const std::array<double3, N>& gradients, const GeneralizedHessian& localHessian)
{
  if (cellLayout.empty()) return;
  std::array<double3, N> shift{};
  bool hasPeriodicShift = false;
  for (std::size_t i = 0; i < N; ++i)
  {
    shift[i] = shiftedPositions[i] - positions[i];
    hasPeriodicShift = hasPeriodicShift || double3::dot(shift[i], shift[i]) > 0.0;
  }
  if (!hasPeriodicShift) return;

  for (std::size_t a = 0; a < cellLayout.size(); ++a)
  {
    const std::size_t cellA = *layout.cellDof(a);
    std::array<double3, N> rawDirection{};
    std::array<double3, N> exactDirection{};
    std::array<double3, N> correctionDirection{};
    for (std::size_t i = 0; i < N; ++i)
    {
      rawDirection[i] = cellLayout.bases[a] * positions[i];
      exactDirection[i] = cellLayout.bases[a] * shiftedPositions[i];
      correctionDirection[i] = exactDirection[i] - rawDirection[i];
    }

    for (std::size_t i = 0; i < N; ++i)
    {
      const std::size_t globalBase = *layout.frameworkAtomDof(atomIds[i], MinimizationDofAxis::X);
      for (std::size_t axisI = 0; axisI < 3; ++axisI)
      {
        double value = 0.0;
        const std::size_t localI = 3 * i + axisI;
        for (std::size_t j = 0; j < N; ++j)
        {
          for (std::size_t axisJ = 0; axisJ < 3; ++axisJ)
          {
            value += localHessian(localI, 3 * j + axisJ) * (&correctionDirection[j].x)[axisJ];
          }
        }
        hessian.add(globalBase + axisI, cellA, value);
        hessian.add(cellA, globalBase + axisI, value);
      }
    }

    for (std::size_t b = 0; b < cellLayout.size(); ++b)
    {
      double value = 0.0;
      for (std::size_t i = 0; i < N; ++i)
      {
        const double3 rawI = rawDirection[i];
        const double3 exactI = exactDirection[i];
        value += double3::dot(gradients[i], cellStrainSecondDerivative(cellLayout, a, b) * shift[i]);
        for (std::size_t j = 0; j < N; ++j)
        {
          const double3 rawJ = cellLayout.bases[b] * positions[j];
          const double3 exactJ = cellLayout.bases[b] * shiftedPositions[j];
          for (std::size_t axisI = 0; axisI < 3; ++axisI)
          {
            for (std::size_t axisJ = 0; axisJ < 3; ++axisJ)
            {
              value += ((&exactI.x)[axisI] * (&exactJ.x)[axisJ] - (&rawI.x)[axisI] * (&rawJ.x)[axisJ]) *
                       localHessian(3 * i + axisI, 3 * j + axisJ);
            }
          }
        }
      }
      hessian.add(cellA, *layout.cellDof(b), value);
    }
  }
}
}  // namespace

RunningEnergy Interactions::computeFrameworkIntraMolecularHessian(
    const ForceField& forceField, const Framework& framework, const SimulationBox& simulationBox,
    std::span<const Atom> atoms, const MinimizationDofLayout& layout, GeneralizedHessian& hessian,
    std::span<AtomDynamics> dynamics, const CellMinimizationLayout& cellLayout)
{
  RunningEnergy energies{};
  if (framework.rigid) return energies;

  const auto shiftedPosition = [&](std::size_t atom, const int3& shift)
  {
    return atoms[atom].position + simulationBox.cell * double3(static_cast<double>(shift.x),
                                                               static_cast<double>(shift.y),
                                                               static_cast<double>(shift.z));
  };
  const auto dofBase = [&](std::size_t atom) { return *layout.frameworkAtomDof(atom, MinimizationDofAxis::X); };

  const auto accumulateRadial = [&](std::size_t A, std::size_t B, const std::array<int3, 2>& shifts, double energy,
                                    double f1, double f2, double& energyAccumulator)
  {
    const double3 positionA = shiftedPosition(A, shifts[0]);
    const double3 positionB = shiftedPosition(B, shifts[1]);
    const double3 dr = positionA - positionB;
    const double3 gradientA = f1 * dr;
    energyAccumulator += energy;
    dynamics[A].gradient += gradientA;
    dynamics[B].gradient -= gradientA;
    hessian.strainGradient().ax += dr.x * gradientA.x;
    hessian.strainGradient().bx += dr.y * gradientA.x;
    hessian.strainGradient().cx += dr.z * gradientA.x;
    hessian.strainGradient().ay += dr.x * gradientA.y;
    hessian.strainGradient().by += dr.y * gradientA.y;
    hessian.strainGradient().cy += dr.z * gradientA.y;
    hessian.strainGradient().az += dr.x * gradientA.z;
    hessian.strainGradient().bz += dr.y * gradientA.z;
    hessian.strainGradient().cz += dr.z * gradientA.z;
    Minimization::scatterAtomicPositionPositionByDof(hessian, dofBase(A), dofBase(B), f1, f2, dr);
    if (!cellLayout.empty() && hasPeriodicShift(shifts))
    {
      GeneralizedHessian localHessian(6, 0);
      Minimization::scatterAtomicPositionPositionByDof(localHessian, 0, 3, f1, f2, dr);
      addFrameworkTermCellCorrection(hessian, layout, cellLayout, std::array<std::size_t, 2>{A, B},
                                     std::array<double3, 2>{atoms[A].position, atoms[B].position},
                                     std::array<double3, 2>{positionA, positionB},
                                     std::array<double3, 2>{gradientA, -gradientA}, localHessian);
    }
  };

  const auto& potentials = framework.intraMolecularPotentials;
  const auto& images = framework.intraMolecularImageShifts;
  std::set<std::array<std::size_t, 2>> pairs12And13;
  std::set<std::array<std::size_t, 2>> pairs14;
  if (!framework.connectivityTable.table.empty())
  {
    for (const std::array<std::size_t, 2>& bond : framework.connectivityTable.findAllBonds())
    {
      pairs12And13.insert({std::min(bond[0], bond[1]), std::max(bond[0], bond[1])});
    }
    for (const std::array<std::size_t, 3>& bend : framework.connectivityTable.findAllBends())
    {
      pairs12And13.insert({std::min(bend[0], bend[2]), std::max(bend[0], bend[2])});
    }
    for (const std::array<std::size_t, 4>& torsion : framework.connectivityTable.findAllTorsions())
    {
      const std::array<std::size_t, 2> pair{std::min(torsion[0], torsion[3]), std::max(torsion[0], torsion[3])};
      if (!pairs12And13.contains(pair)) pairs14.insert(pair);
    }
  }
  for (std::size_t index = 0; index < potentials.bonds.size(); ++index)
  {
    const BondPotential& bond = potentials.bonds[index];
    const std::size_t A = bond.identifiers[0];
    const std::size_t B = bond.identifiers[1];
    const double3 positionA = shiftedPosition(A, images.bonds[index][0]);
    const double3 positionB = shiftedPosition(B, images.bonds[index][1]);
    auto [energy, gradient, strain, f1, f2] = bond.potentialEnergyGradientHessianStrain(positionA, positionB);
    energies.bond += energy;
    dynamics[A].gradient += gradient[0];
    dynamics[B].gradient += gradient[1];
    hessian.strainGradient() += strain;
    Minimization::scatterAtomicPositionPositionByDof(hessian, dofBase(A), dofBase(B), f1, f2, positionA - positionB);
    if (!cellLayout.empty() && hasPeriodicShift(images.bonds[index]))
    {
      GeneralizedHessian localHessian(6, 0);
      Minimization::scatterAtomicPositionPositionByDof(localHessian, 0, 3, f1, f2, positionA - positionB);
      addFrameworkTermCellCorrection(hessian, layout, cellLayout, std::array<std::size_t, 2>{A, B},
                                     std::array<double3, 2>{atoms[A].position, atoms[B].position},
                                     std::array<double3, 2>{positionA, positionB}, gradient, localHessian);
    }
  }

  for (std::size_t index = 0; index < potentials.bends.size(); ++index)
  {
    const BendPotential& bend = potentials.bends[index];
    const auto& ids = bend.identifiers;
    const auto& shifts = images.bends[index];
    const std::array<double3, 3> shifted = {shiftedPosition(ids[0], shifts[0]), shiftedPosition(ids[1], shifts[1]),
                                            shiftedPosition(ids[2], shifts[2])};
    auto [energy, gradient, strain, geometry] = Potentials::Internal::bendPotentialEnergyGradientHessianStrain(
        bend.type, bend.parameters, shifted[0], shifted[1], shifted[2]);
    energies.bend += energy;
    for (std::size_t i = 0; i < 3; ++i) dynamics[ids[i]].gradient += gradient[i];
    hessian.strainGradient() += strain;
    Minimization::scatterBendHessianByDof(hessian, {dofBase(ids[0]), dofBase(ids[1]), dofBase(ids[2])}, geometry);
    if (!cellLayout.empty() && hasPeriodicShift(shifts))
    {
      GeneralizedHessian localHessian(9, 0);
      Minimization::scatterBendHessianByDof(localHessian, {0, 3, 6}, geometry);
      addFrameworkTermCellCorrection(
          hessian, layout, cellLayout, ids,
          std::array<double3, 3>{atoms[ids[0]].position, atoms[ids[1]].position, atoms[ids[2]].position}, shifted,
          gradient, localHessian);
    }
  }

  const auto accumulateTorsions = [&](const std::vector<TorsionPotential>& torsions,
                                      const std::vector<std::array<int3, 4>>& termImages,
                                      double RunningEnergy::* energyMember)
  {
    for (std::size_t index = 0; index < torsions.size(); ++index)
    {
      const TorsionPotential& torsion = torsions[index];
      const auto& ids = torsion.identifiers;
      const auto& shifts = termImages[index];
      const std::array<double3, 4> shifted = {shiftedPosition(ids[0], shifts[0]), shiftedPosition(ids[1], shifts[1]),
                                              shiftedPosition(ids[2], shifts[2]), shiftedPosition(ids[3], shifts[3])};
      auto [energy, gradient, strain, geometry] = Potentials::Internal::torsionPotentialEnergyGradientHessianStrain(
          torsion.type, torsion.parameters, shifted[0], shifted[1], shifted[2], shifted[3]);
      energies.*energyMember += energy;
      for (std::size_t i = 0; i < 4; ++i) dynamics[ids[i]].gradient += gradient[i];
      hessian.strainGradient() += strain;
      Minimization::scatterTorsionHessianByDof(
          hessian, {dofBase(ids[0]), dofBase(ids[1]), dofBase(ids[2]), dofBase(ids[3])}, geometry);
      if (!cellLayout.empty() && hasPeriodicShift(shifts))
      {
        GeneralizedHessian localHessian(12, 0);
        Minimization::scatterTorsionHessianByDof(localHessian, {0, 3, 6, 9}, geometry);
        addFrameworkTermCellCorrection(hessian, layout, cellLayout, ids,
                                       std::array<double3, 4>{atoms[ids[0]].position, atoms[ids[1]].position,
                                                              atoms[ids[2]].position, atoms[ids[3]].position},
                                       shifted, gradient, localHessian);
      }
    }
  };
  accumulateTorsions(potentials.torsions, images.torsions, &RunningEnergy::torsion);
  accumulateTorsions(potentials.improperTorsions, images.improperTorsions, &RunningEnergy::improperTorsion);

  for (std::size_t index = 0; index < potentials.vanDerWaals.size(); ++index)
  {
    const VanDerWaalsPotential& potential = potentials.vanDerWaals[index];
    const std::size_t A = potential.identifiers[0];
    const std::size_t B = potential.identifiers[1];
    const double3 dr =
        shiftedPosition(A, images.vanDerWaals[index][0]) - shiftedPosition(B, images.vanDerWaals[index][1]);
    const double rr = double3::dot(dr, dr);
    const bool isScaled14 =
        pairs14.contains({std::min(A, B), std::max(A, B)}) && potential.scaling != 1.0;
    if (isScaled14 || rr < forceField.cutOffFrameworkVDW * forceField.cutOffFrameworkVDW)
    {
      const Potentials::PairDerivatives<2> factors = Potentials::potentialVDW<2>(
          forceField, 1.0, 1.0, rr, static_cast<std::size_t>(atoms[A].type), static_cast<std::size_t>(atoms[B].type));
      accumulateRadial(A, B, images.vanDerWaals[index], potential.scaling * factors.energy,
                       potential.scaling * factors.firstDerivativeFactor,
                       potential.scaling * factors.secondDerivativeFactor, energies.intraVDW);
    }
  }

  for (std::size_t index = 0; index < potentials.coulombs.size(); ++index)
  {
    const CoulombPotential& potential = potentials.coulombs[index];
    const std::size_t A = potential.identifiers[0];
    const std::size_t B = potential.identifiers[1];
    const double3 dr = shiftedPosition(A, images.coulombs[index][0]) - shiftedPosition(B, images.coulombs[index][1]);
    const double rr = double3::dot(dr, dr);
    if (pairs14.contains({std::min(A, B), std::max(A, B)}) && potential.scaling != 1.0)
    {
      const double r = std::sqrt(rr);
      const double k = potential.scaling * Units::CoulombicConversionFactor * potential.chargeA * potential.chargeB;
      accumulateRadial(A, B, images.coulombs[index], k / r, -k / (r * rr), 3.0 * k / (r * rr * rr),
                       energies.intraCoul);
    }
    else if (rr < forceField.cutOffCoulomb * forceField.cutOffCoulomb)
    {
      const double r = std::sqrt(rr);
      const Potentials::PairDerivatives<2> factors =
          Potentials::potentialCoulomb<2>(forceField, potential.scaling, 1.0, r, potential.chargeA, potential.chargeB);
      accumulateRadial(A, B, images.coulombs[index], factors.energy, factors.firstDerivativeFactor,
                       factors.secondDerivativeFactor, energies.intraCoul);
    }
  }

  return energies;
}

RunningEnergy Interactions::computeIntraMolecularBondHessian(
    std::span<const Molecule> moleculeData, std::span<const Atom> atoms, std::span<const Component> components,
    const MinimizationDofLayout& layout, GeneralizedHessian& hessian, std::span<AtomDynamics> dynamics)
{
  RunningEnergy energies{};

  const Minimization::RigidDerivativeCache rigidCache =
      Minimization::RigidDerivativeCache::build(moleculeData, components, atoms);

  for (std::size_t moleculeIndex = 0; moleculeIndex < moleculeData.size(); ++moleculeIndex)
  {
    const Molecule& molecule = moleculeData[moleculeIndex];
    if (components[molecule.componentId].rigid)
    {
      continue;
    }
    const bool semiFlexible = components[molecule.componentId].isSemiFlexible();

    const Potentials::IntraMolecularPotentials& potentials = components[molecule.componentId].intraMolecularPotentials;
    std::span<const Atom> atom_molecule_span = {&atoms[molecule.atomIndex], molecule.numberOfAtoms};
    std::span<AtomDynamics> dynamics_molecule_span = {&dynamics[molecule.atomIndex], molecule.numberOfAtoms};

    for (const BondPotential& bond : potentials.bonds)
    {
      const std::size_t A = bond.identifiers[0];
      const std::size_t B = bond.identifiers[1];
      auto [energy, gradient, strain, f1, f2] =
          bond.potentialEnergyGradientHessianStrain(atom_molecule_span[A].position, atom_molecule_span[B].position);
      energies.bond += energy;
      dynamics_molecule_span[A].gradient += gradient[0];
      dynamics_molecule_span[B].gradient += gradient[1];
      hessian.strainGradient() += strain;

      const double3 dr = atom_molecule_span[A].position - atom_molecule_span[B].position;
      if (semiFlexible)
      {
        removeRigidOffsetStrainGradient<2>(hessian, layout, rigidCache, moleculeIndex, atom_molecule_span, {A, B},
                                           {gradient[0], gradient[1]});
        Minimization::scatterRadialHessianSites(
            hessian, {Minimization::makeHessianSite(layout, rigidCache, moleculeIndex, A),
                      Minimization::makeHessianSite(layout, rigidCache, moleculeIndex, B)},
            {gradient[0], gradient[1]}, f1, f2, dr);
        continue;
      }
      Minimization::scatterAtomicPositionPosition(hessian, layout, moleculeIndex, A, moleculeIndex, B, f1, f2, dr);
      if (hessian.numStrain() == 1)
      {
        Minimization::scatterAtomicPositionStrainIsotropic(hessian, layout, moleculeIndex, A, moleculeIndex, B, f1, f2,
                                                           dr);
        Minimization::scatterAtomicStrainStrainIsotropic(hessian, f1, f2, dr, atom_molecule_span[A].position,
                                                         atom_molecule_span[A].position, atom_molecule_span[B].position,
                                                         atom_molecule_span[B].position, false, false);
      }
    }
  }

  return energies;
}

RunningEnergy Interactions::computeIntraMolecularBendHessian(
    std::span<const Molecule> moleculeData, std::span<const Atom> atoms, std::span<const Component> components,
    const MinimizationDofLayout& layout, GeneralizedHessian& hessian, std::span<AtomDynamics> dynamics)
{
  RunningEnergy energies{};

  const Minimization::RigidDerivativeCache rigidCache =
      Minimization::RigidDerivativeCache::build(moleculeData, components, atoms);

  for (std::size_t moleculeIndex = 0; moleculeIndex < moleculeData.size(); ++moleculeIndex)
  {
    const Molecule& molecule = moleculeData[moleculeIndex];
    if (components[molecule.componentId].rigid)
    {
      continue;
    }
    const bool semiFlexible = components[molecule.componentId].isSemiFlexible();

    const Potentials::IntraMolecularPotentials& potentials = components[molecule.componentId].intraMolecularPotentials;
    std::span<const Atom> atom_molecule_span = {&atoms[molecule.atomIndex], molecule.numberOfAtoms};
    std::span<AtomDynamics> dynamics_molecule_span = {&dynamics[molecule.atomIndex], molecule.numberOfAtoms};

    for (const BendPotential& bend : potentials.bends)
    {
      const std::size_t A = bend.identifiers[0];
      const std::size_t B = bend.identifiers[1];
      const std::size_t C = bend.identifiers[2];
      auto [energy, gradient, strain, geometry] = Potentials::Internal::bendPotentialEnergyGradientHessianStrain(
          bend.type, bend.parameters, atom_molecule_span[A].position, atom_molecule_span[B].position,
          atom_molecule_span[C].position);
      energies.bend += energy;
      dynamics_molecule_span[A].gradient += gradient[0];
      dynamics_molecule_span[B].gradient += gradient[1];
      dynamics_molecule_span[C].gradient += gradient[2];
      hessian.strainGradient() += strain;
      if (semiFlexible)
      {
        removeRigidOffsetStrainGradient<3>(hessian, layout, rigidCache, moleculeIndex, atom_molecule_span, {A, B, C},
                                           {gradient[0], gradient[1], gradient[2]});
        Minimization::scatterBendHessianSites(
            hessian,
            {Minimization::makeHessianSite(layout, rigidCache, moleculeIndex, A),
             Minimization::makeHessianSite(layout, rigidCache, moleculeIndex, B),
             Minimization::makeHessianSite(layout, rigidCache, moleculeIndex, C)},
            {gradient[0], gradient[1], gradient[2]}, geometry);
        continue;
      }
      Minimization::scatterBendHessian(hessian, layout, moleculeIndex, A, B, C, geometry);
    }
  }

  return energies;
}

RunningEnergy Interactions::computeIntraMolecularUreyBradleyHessian(
    std::span<const Molecule> moleculeData, std::span<const Atom> atoms, std::span<const Component> components,
    const MinimizationDofLayout& layout, GeneralizedHessian& hessian, std::span<AtomDynamics> dynamics)
{
  RunningEnergy energies{};

  const Minimization::RigidDerivativeCache rigidCache =
      Minimization::RigidDerivativeCache::build(moleculeData, components, atoms);

  for (std::size_t moleculeIndex = 0; moleculeIndex < moleculeData.size(); ++moleculeIndex)
  {
    const Molecule& molecule = moleculeData[moleculeIndex];
    if (components[molecule.componentId].rigid)
    {
      continue;
    }
    const bool semiFlexible = components[molecule.componentId].isSemiFlexible();

    const Potentials::IntraMolecularPotentials& potentials = components[molecule.componentId].intraMolecularPotentials;
    std::span<const Atom> atom_molecule_span = {&atoms[molecule.atomIndex], molecule.numberOfAtoms};
    std::span<AtomDynamics> dynamics_molecule_span = {&dynamics[molecule.atomIndex], molecule.numberOfAtoms};

    for (const UreyBradleyPotential& ureyBradley : potentials.ureyBradleys)
    {
      const std::size_t A = ureyBradley.identifiers[0];
      const std::size_t B = ureyBradley.identifiers[1];
      auto [energy, gradient, strain, f1, f2] = Potentials::Internal::distancePotentialEnergyGradientHessianStrain(
          bondTypeFromUreyBradley(ureyBradley.type), ureyBradley.parameters, atom_molecule_span[A].position,
          atom_molecule_span[B].position);
      energies.ureyBradley += energy;
      dynamics_molecule_span[A].gradient += gradient[0];
      dynamics_molecule_span[B].gradient += gradient[1];
      hessian.strainGradient() += strain;

      const double3 dr = atom_molecule_span[A].position - atom_molecule_span[B].position;
      if (semiFlexible)
      {
        removeRigidOffsetStrainGradient<2>(hessian, layout, rigidCache, moleculeIndex, atom_molecule_span, {A, B},
                                           {gradient[0], gradient[1]});
        Minimization::scatterRadialHessianSites(
            hessian, {Minimization::makeHessianSite(layout, rigidCache, moleculeIndex, A),
                      Minimization::makeHessianSite(layout, rigidCache, moleculeIndex, B)},
            {gradient[0], gradient[1]}, f1, f2, dr);
        continue;
      }
      Minimization::scatterAtomicPositionPosition(hessian, layout, moleculeIndex, A, moleculeIndex, B, f1, f2, dr);
      if (hessian.numStrain() == 1)
      {
        Minimization::scatterAtomicPositionStrainIsotropic(hessian, layout, moleculeIndex, A, moleculeIndex, B, f1, f2,
                                                           dr);
        Minimization::scatterAtomicStrainStrainIsotropic(hessian, f1, f2, dr, atom_molecule_span[A].position,
                                                         atom_molecule_span[A].position, atom_molecule_span[B].position,
                                                         atom_molecule_span[B].position, false, false);
      }
    }
  }

  return energies;
}

RunningEnergy Interactions::computeIntraMolecularVanDerWaalsHessian(
    std::span<const Molecule> moleculeData, std::span<const Atom> atoms, std::span<const Component> components,
    const MinimizationDofLayout& layout, GeneralizedHessian& hessian, std::span<AtomDynamics> dynamics)
{
  RunningEnergy energies{};

  const Minimization::RigidDerivativeCache rigidCache =
      Minimization::RigidDerivativeCache::build(moleculeData, components, atoms);

  for (std::size_t moleculeIndex = 0; moleculeIndex < moleculeData.size(); ++moleculeIndex)
  {
    const Molecule& molecule = moleculeData[moleculeIndex];
    if (components[molecule.componentId].rigid)
    {
      continue;
    }
    const bool semiFlexible = components[molecule.componentId].isSemiFlexible();

    const Potentials::IntraMolecularPotentials& potentials = components[molecule.componentId].intraMolecularPotentials;
    std::span<const Atom> atom_molecule_span = {&atoms[molecule.atomIndex], molecule.numberOfAtoms};
    std::span<AtomDynamics> dynamics_molecule_span = {&dynamics[molecule.atomIndex], molecule.numberOfAtoms};

    for (const VanDerWaalsPotential& vanDerWaals : potentials.vanDerWaals)
    {
      const std::size_t A = vanDerWaals.identifiers[0];
      const std::size_t B = vanDerWaals.identifiers[1];
      const double3 dr = atom_molecule_span[A].position - atom_molecule_span[B].position;
      const double rr = double3::dot(dr, dr);

      // Lennard-Jones: 4*eps*((sigma/r)^12 - (sigma/r)^6), parameters = {eps, sigma}.
      const double sigmaOverR2 = (vanDerWaals.parameters[1] * vanDerWaals.parameters[1]) / rr;
      const double t = sigmaOverR2 * sigmaOverR2 * sigmaOverR2;
      const double prefactor = vanDerWaals.scaling * vanDerWaals.parameters[0];
      const double energy = 4.0 * prefactor * t * (t - 1.0);
      const double f1 = 24.0 * prefactor * t * (1.0 - 2.0 * t) / rr;
      const double f2 = 96.0 * prefactor * t * (7.0 * t - 2.0) / (rr * rr);

      energies.intraVDW += energy;
      const double3 gradientA = f1 * dr;
      dynamics_molecule_span[A].gradient += gradientA;
      dynamics_molecule_span[B].gradient -= gradientA;

      double3x3 strain{};
      strain.ax = dr.x * gradientA.x;
      strain.bx = dr.y * gradientA.x;
      strain.cx = dr.z * gradientA.x;
      strain.ay = dr.x * gradientA.y;
      strain.by = dr.y * gradientA.y;
      strain.cy = dr.z * gradientA.y;
      strain.az = dr.x * gradientA.z;
      strain.bz = dr.y * gradientA.z;
      strain.cz = dr.z * gradientA.z;
      hessian.strainGradient() += strain;

      if (semiFlexible)
      {
        removeRigidOffsetStrainGradient<2>(hessian, layout, rigidCache, moleculeIndex, atom_molecule_span, {A, B},
                                           {gradientA, -1.0 * gradientA});
        Minimization::scatterRadialHessianSites(
            hessian, {Minimization::makeHessianSite(layout, rigidCache, moleculeIndex, A),
                      Minimization::makeHessianSite(layout, rigidCache, moleculeIndex, B)},
            {gradientA, -1.0 * gradientA}, f1, f2, dr);
        continue;
      }
      Minimization::scatterAtomicPositionPosition(hessian, layout, moleculeIndex, A, moleculeIndex, B, f1, f2, dr);
      if (hessian.numStrain() == 1)
      {
        Minimization::scatterAtomicPositionStrainIsotropic(hessian, layout, moleculeIndex, A, moleculeIndex, B, f1, f2,
                                                           dr);
        Minimization::scatterAtomicStrainStrainIsotropic(hessian, f1, f2, dr, atom_molecule_span[A].position,
                                                         atom_molecule_span[A].position, atom_molecule_span[B].position,
                                                         atom_molecule_span[B].position, false, false);
      }
    }
  }

  return energies;
}

RunningEnergy Interactions::computeIntraMolecularCoulombHessian(
    std::span<const Molecule> moleculeData, std::span<const Atom> atoms, std::span<const Component> components,
    const MinimizationDofLayout& layout, GeneralizedHessian& hessian, std::span<AtomDynamics> dynamics)
{
  RunningEnergy energies{};

  const Minimization::RigidDerivativeCache rigidCache =
      Minimization::RigidDerivativeCache::build(moleculeData, components, atoms);

  for (std::size_t moleculeIndex = 0; moleculeIndex < moleculeData.size(); ++moleculeIndex)
  {
    const Molecule& molecule = moleculeData[moleculeIndex];
    if (components[molecule.componentId].rigid)
    {
      continue;
    }
    const bool semiFlexible = components[molecule.componentId].isSemiFlexible();

    const Potentials::IntraMolecularPotentials& potentials = components[molecule.componentId].intraMolecularPotentials;
    std::span<const Atom> atom_molecule_span = {&atoms[molecule.atomIndex], molecule.numberOfAtoms};
    std::span<AtomDynamics> dynamics_molecule_span = {&dynamics[molecule.atomIndex], molecule.numberOfAtoms};

    for (const CoulombPotential& coulomb : potentials.coulombs)
    {
      const std::size_t A = coulomb.identifiers[0];
      const std::size_t B = coulomb.identifiers[1];
      const double3 dr = atom_molecule_span[A].position - atom_molecule_span[B].position;
      const double rr = double3::dot(dr, dr);
      const double r = std::sqrt(rr);

      // U = k/r with k = scaling * conversion * qA * qB.
      const double k = coulomb.scaling * Units::CoulombicConversionFactor * coulomb.chargeA * coulomb.chargeB;
      const double energy = k / r;
      const double f1 = -k / (r * rr);
      const double f2 = 3.0 * k / (r * rr * rr);

      energies.intraCoul += energy;
      const double3 gradientA = f1 * dr;
      dynamics_molecule_span[A].gradient += gradientA;
      dynamics_molecule_span[B].gradient -= gradientA;

      double3x3 strain{};
      strain.ax = dr.x * gradientA.x;
      strain.bx = dr.y * gradientA.x;
      strain.cx = dr.z * gradientA.x;
      strain.ay = dr.x * gradientA.y;
      strain.by = dr.y * gradientA.y;
      strain.cy = dr.z * gradientA.y;
      strain.az = dr.x * gradientA.z;
      strain.bz = dr.y * gradientA.z;
      strain.cz = dr.z * gradientA.z;
      hessian.strainGradient() += strain;

      if (semiFlexible)
      {
        removeRigidOffsetStrainGradient<2>(hessian, layout, rigidCache, moleculeIndex, atom_molecule_span, {A, B},
                                           {gradientA, -1.0 * gradientA});
        Minimization::scatterRadialHessianSites(
            hessian, {Minimization::makeHessianSite(layout, rigidCache, moleculeIndex, A),
                      Minimization::makeHessianSite(layout, rigidCache, moleculeIndex, B)},
            {gradientA, -1.0 * gradientA}, f1, f2, dr);
        continue;
      }
      Minimization::scatterAtomicPositionPosition(hessian, layout, moleculeIndex, A, moleculeIndex, B, f1, f2, dr);
      if (hessian.numStrain() == 1)
      {
        Minimization::scatterAtomicPositionStrainIsotropic(hessian, layout, moleculeIndex, A, moleculeIndex, B, f1, f2,
                                                           dr);
        Minimization::scatterAtomicStrainStrainIsotropic(hessian, f1, f2, dr, atom_molecule_span[A].position,
                                                         atom_molecule_span[A].position, atom_molecule_span[B].position,
                                                         atom_molecule_span[B].position, false, false);
      }
    }
  }

  return energies;
}

RunningEnergy Interactions::computeIntraMolecularTorsionHessian(
    std::span<const Molecule> moleculeData, std::span<const Atom> atoms, std::span<const Component> components,
    const MinimizationDofLayout& layout, GeneralizedHessian& hessian, std::span<AtomDynamics> dynamics)
{
  RunningEnergy energies{};

  const Minimization::RigidDerivativeCache rigidCache =
      Minimization::RigidDerivativeCache::build(moleculeData, components, atoms);

  for (std::size_t moleculeIndex = 0; moleculeIndex < moleculeData.size(); ++moleculeIndex)
  {
    const Molecule& molecule = moleculeData[moleculeIndex];
    if (components[molecule.componentId].rigid)
    {
      continue;
    }
    const bool semiFlexible = components[molecule.componentId].isSemiFlexible();

    const Potentials::IntraMolecularPotentials& potentials = components[molecule.componentId].intraMolecularPotentials;
    std::span<const Atom> atom_molecule_span = {&atoms[molecule.atomIndex], molecule.numberOfAtoms};
    std::span<AtomDynamics> dynamics_molecule_span = {&dynamics[molecule.atomIndex], molecule.numberOfAtoms};

    auto accumulateTorsion = [&](const TorsionPotential& torsion, double& energyAccumulator)
    {
      const std::size_t A = torsion.identifiers[0];
      const std::size_t B = torsion.identifiers[1];
      const std::size_t C = torsion.identifiers[2];
      const std::size_t D = torsion.identifiers[3];
      auto [energy, gradient, strain, geometry] = Potentials::Internal::torsionPotentialEnergyGradientHessianStrain(
          torsion.type, torsion.parameters, atom_molecule_span[A].position, atom_molecule_span[B].position,
          atom_molecule_span[C].position, atom_molecule_span[D].position);
      energyAccumulator += energy;
      dynamics_molecule_span[A].gradient += gradient[0];
      dynamics_molecule_span[B].gradient += gradient[1];
      dynamics_molecule_span[C].gradient += gradient[2];
      dynamics_molecule_span[D].gradient += gradient[3];
      hessian.strainGradient() += strain;
      if (semiFlexible)
      {
        removeRigidOffsetStrainGradient<4>(hessian, layout, rigidCache, moleculeIndex, atom_molecule_span,
                                           {A, B, C, D}, {gradient[0], gradient[1], gradient[2], gradient[3]});
        Minimization::scatterTorsionHessianSites(
            hessian,
            {Minimization::makeHessianSite(layout, rigidCache, moleculeIndex, A),
             Minimization::makeHessianSite(layout, rigidCache, moleculeIndex, B),
             Minimization::makeHessianSite(layout, rigidCache, moleculeIndex, C),
             Minimization::makeHessianSite(layout, rigidCache, moleculeIndex, D)},
            {gradient[0], gradient[1], gradient[2], gradient[3]}, geometry);
        return;
      }
      Minimization::scatterTorsionHessian(hessian, layout, moleculeIndex, A, B, C, D, geometry);
    };

    for (const TorsionPotential& torsion : potentials.torsions)
    {
      accumulateTorsion(torsion, energies.torsion);
    }
    for (const TorsionPotential& improperTorsion : potentials.improperTorsions)
    {
      accumulateTorsion(improperTorsion, energies.improperTorsion);
    }
  }

  return energies;
}
