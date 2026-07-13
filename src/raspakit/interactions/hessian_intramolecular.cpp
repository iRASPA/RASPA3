module;

module interactions_hessian_intramolecular;

import std;

import double3;
import double3x3;
import units;
import intra_molecular_potentials;
import bond_potential;
import bend_potential;
import urey_bradley_potential;
import van_der_waals_potential;
import coulomb_potential;
import torsion_potential;
import bend_potential_gradient_hessian_strain;
import torsion_potential_gradient_hessian_strain;
import distance_potential_gradient_hessian_strain;
import minimization_hessian_scatter;

namespace
{
// UreyBradleyType enumerators are the BondType enumerators shifted by one (no 'None' entry).
BondType bondTypeFromUreyBradley(UreyBradleyType type)
{
  return static_cast<BondType>(static_cast<std::size_t>(type) + 1);
}
}  // namespace

RunningEnergy Interactions::computeIntraMolecularBondHessian(std::span<const Molecule> moleculeData,
                                                           std::span<const Atom> atoms,
                                                           std::span<const Component> components,
                                                           const MinimizationDofLayout &layout,
                                                           GeneralizedHessian &hessian,
                                                           std::span<AtomDynamics> dynamics)
{
  RunningEnergy energies{};

  for (std::size_t moleculeIndex = 0; moleculeIndex < moleculeData.size(); ++moleculeIndex)
  {
    const Molecule &molecule = moleculeData[moleculeIndex];
    if (components[molecule.componentId].rigid)
    {
      continue;
    }

    const Potentials::IntraMolecularPotentials &potentials = components[molecule.componentId].intraMolecularPotentials;
    std::span<const Atom> atom_molecule_span = {&atoms[molecule.atomIndex], molecule.numberOfAtoms};
    std::span<AtomDynamics> dynamics_molecule_span = {&dynamics[molecule.atomIndex], molecule.numberOfAtoms};

    for (const BondPotential &bond : potentials.bonds)
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
      Minimization::scatterAtomicPositionPosition(hessian, layout, moleculeIndex, A, moleculeIndex, B, f1, f2, dr);
      if (hessian.numStrain() == 1)
      {
        Minimization::scatterAtomicPositionStrainIsotropic(hessian, layout, moleculeIndex, A, moleculeIndex, B, f1, f2,
                                                           dr);
        Minimization::scatterAtomicStrainStrainIsotropic(
            hessian, f1, f2, dr, atom_molecule_span[A].position, atom_molecule_span[A].position,
            atom_molecule_span[B].position, atom_molecule_span[B].position, false, false);
      }
    }
  }

  return energies;
}

RunningEnergy Interactions::computeIntraMolecularBendHessian(std::span<const Molecule> moleculeData,
                                                             std::span<const Atom> atoms,
                                                             std::span<const Component> components,
                                                             const MinimizationDofLayout &layout,
                                                             GeneralizedHessian &hessian,
                                                             std::span<AtomDynamics> dynamics)
{
  RunningEnergy energies{};

  for (std::size_t moleculeIndex = 0; moleculeIndex < moleculeData.size(); ++moleculeIndex)
  {
    const Molecule &molecule = moleculeData[moleculeIndex];
    if (components[molecule.componentId].rigid)
    {
      continue;
    }

    const Potentials::IntraMolecularPotentials &potentials = components[molecule.componentId].intraMolecularPotentials;
    std::span<const Atom> atom_molecule_span = {&atoms[molecule.atomIndex], molecule.numberOfAtoms};
    std::span<AtomDynamics> dynamics_molecule_span = {&dynamics[molecule.atomIndex], molecule.numberOfAtoms};

    for (const BendPotential &bend : potentials.bends)
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
      Minimization::scatterBendHessian(hessian, layout, moleculeIndex, A, B, C, geometry);
    }
  }

  return energies;
}

RunningEnergy Interactions::computeIntraMolecularUreyBradleyHessian(std::span<const Molecule> moleculeData,
                                                                    std::span<const Atom> atoms,
                                                                    std::span<const Component> components,
                                                                    const MinimizationDofLayout &layout,
                                                                    GeneralizedHessian &hessian,
                                                                    std::span<AtomDynamics> dynamics)
{
  RunningEnergy energies{};

  for (std::size_t moleculeIndex = 0; moleculeIndex < moleculeData.size(); ++moleculeIndex)
  {
    const Molecule &molecule = moleculeData[moleculeIndex];
    if (components[molecule.componentId].rigid)
    {
      continue;
    }

    const Potentials::IntraMolecularPotentials &potentials = components[molecule.componentId].intraMolecularPotentials;
    std::span<const Atom> atom_molecule_span = {&atoms[molecule.atomIndex], molecule.numberOfAtoms};
    std::span<AtomDynamics> dynamics_molecule_span = {&dynamics[molecule.atomIndex], molecule.numberOfAtoms};

    for (const UreyBradleyPotential &ureyBradley : potentials.ureyBradleys)
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
      Minimization::scatterAtomicPositionPosition(hessian, layout, moleculeIndex, A, moleculeIndex, B, f1, f2, dr);
      if (hessian.numStrain() == 1)
      {
        Minimization::scatterAtomicPositionStrainIsotropic(hessian, layout, moleculeIndex, A, moleculeIndex, B, f1, f2,
                                                           dr);
        Minimization::scatterAtomicStrainStrainIsotropic(
            hessian, f1, f2, dr, atom_molecule_span[A].position, atom_molecule_span[A].position,
            atom_molecule_span[B].position, atom_molecule_span[B].position, false, false);
      }
    }
  }

  return energies;
}

RunningEnergy Interactions::computeIntraMolecularVanDerWaalsHessian(std::span<const Molecule> moleculeData,
                                                                    std::span<const Atom> atoms,
                                                                    std::span<const Component> components,
                                                                    const MinimizationDofLayout &layout,
                                                                    GeneralizedHessian &hessian,
                                                                    std::span<AtomDynamics> dynamics)
{
  RunningEnergy energies{};

  for (std::size_t moleculeIndex = 0; moleculeIndex < moleculeData.size(); ++moleculeIndex)
  {
    const Molecule &molecule = moleculeData[moleculeIndex];
    if (components[molecule.componentId].rigid)
    {
      continue;
    }

    const Potentials::IntraMolecularPotentials &potentials = components[molecule.componentId].intraMolecularPotentials;
    std::span<const Atom> atom_molecule_span = {&atoms[molecule.atomIndex], molecule.numberOfAtoms};
    std::span<AtomDynamics> dynamics_molecule_span = {&dynamics[molecule.atomIndex], molecule.numberOfAtoms};

    for (const VanDerWaalsPotential &vanDerWaals : potentials.vanDerWaals)
    {
      const std::size_t A = vanDerWaals.identifiers[0];
      const std::size_t B = vanDerWaals.identifiers[1];
      const double3 dr = atom_molecule_span[A].position - atom_molecule_span[B].position;
      const double rr = double3::dot(dr, dr);
      const double r = std::sqrt(rr);

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

      Minimization::scatterAtomicPositionPosition(hessian, layout, moleculeIndex, A, moleculeIndex, B, f1, f2, dr);
      if (hessian.numStrain() == 1)
      {
        Minimization::scatterAtomicPositionStrainIsotropic(hessian, layout, moleculeIndex, A, moleculeIndex, B, f1, f2,
                                                           dr);
        Minimization::scatterAtomicStrainStrainIsotropic(
            hessian, f1, f2, dr, atom_molecule_span[A].position, atom_molecule_span[A].position,
            atom_molecule_span[B].position, atom_molecule_span[B].position, false, false);
      }
    }
  }

  return energies;
}

RunningEnergy Interactions::computeIntraMolecularCoulombHessian(std::span<const Molecule> moleculeData,
                                                                std::span<const Atom> atoms,
                                                                std::span<const Component> components,
                                                                const MinimizationDofLayout &layout,
                                                                GeneralizedHessian &hessian,
                                                                std::span<AtomDynamics> dynamics)
{
  RunningEnergy energies{};

  for (std::size_t moleculeIndex = 0; moleculeIndex < moleculeData.size(); ++moleculeIndex)
  {
    const Molecule &molecule = moleculeData[moleculeIndex];
    if (components[molecule.componentId].rigid)
    {
      continue;
    }

    const Potentials::IntraMolecularPotentials &potentials = components[molecule.componentId].intraMolecularPotentials;
    std::span<const Atom> atom_molecule_span = {&atoms[molecule.atomIndex], molecule.numberOfAtoms};
    std::span<AtomDynamics> dynamics_molecule_span = {&dynamics[molecule.atomIndex], molecule.numberOfAtoms};

    for (const CoulombPotential &coulomb : potentials.coulombs)
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

      Minimization::scatterAtomicPositionPosition(hessian, layout, moleculeIndex, A, moleculeIndex, B, f1, f2, dr);
      if (hessian.numStrain() == 1)
      {
        Minimization::scatterAtomicPositionStrainIsotropic(hessian, layout, moleculeIndex, A, moleculeIndex, B, f1, f2,
                                                           dr);
        Minimization::scatterAtomicStrainStrainIsotropic(
            hessian, f1, f2, dr, atom_molecule_span[A].position, atom_molecule_span[A].position,
            atom_molecule_span[B].position, atom_molecule_span[B].position, false, false);
      }
    }
  }

  return energies;
}

RunningEnergy Interactions::computeIntraMolecularTorsionHessian(std::span<const Molecule> moleculeData,
                                                                std::span<const Atom> atoms,
                                                                std::span<const Component> components,
                                                                const MinimizationDofLayout &layout,
                                                                GeneralizedHessian &hessian,
                                                                std::span<AtomDynamics> dynamics)
{
  RunningEnergy energies{};

  for (std::size_t moleculeIndex = 0; moleculeIndex < moleculeData.size(); ++moleculeIndex)
  {
    const Molecule &molecule = moleculeData[moleculeIndex];
    if (components[molecule.componentId].rigid)
    {
      continue;
    }

    const Potentials::IntraMolecularPotentials &potentials = components[molecule.componentId].intraMolecularPotentials;
    std::span<const Atom> atom_molecule_span = {&atoms[molecule.atomIndex], molecule.numberOfAtoms};
    std::span<AtomDynamics> dynamics_molecule_span = {&dynamics[molecule.atomIndex], molecule.numberOfAtoms};

    auto accumulateTorsion = [&](const TorsionPotential &torsion, double &energyAccumulator)
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
      Minimization::scatterTorsionHessian(hessian, layout, moleculeIndex, A, B, C, D, geometry);
    };

    for (const TorsionPotential &torsion : potentials.torsions)
    {
      accumulateTorsion(torsion, energies.torsion);
    }
    for (const TorsionPotential &improperTorsion : potentials.improperTorsions)
    {
      accumulateTorsion(improperTorsion, energies.improperTorsion);
    }
  }

  return energies;
}
