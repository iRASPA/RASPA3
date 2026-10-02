module;

module spatial_decomposition_device_bonded_topology;

import std;

import atom;
import molecule;
import component;
import forcefield;
import system;
import bond_potential;
import bend_potential;
import torsion_potential;
import van_der_waals_potential;
import coulomb_potential;
import units;
import intra_molecular_potentials;

bool BondedTopology::supports(const System& system, std::string& reason)
{
  for (const Component& component : system.components)
  {
    const Potentials::IntraMolecularPotentials& potentials = component.intraMolecularPotentials;
    const char* unsupported = nullptr;
    if (!potentials.ureyBradleys.empty())
      unsupported = "Urey-Bradley";
    else if (!potentials.inversionBends.empty())
      unsupported = "inversion-bend";
    else if (!potentials.outOfPlaneBends.empty())
      unsupported = "out-of-plane-bend";
    else if (!potentials.bondBonds.empty())
      unsupported = "bond-bond";
    else if (!potentials.bondBends.empty())
      unsupported = "bond-bend";
    else if (!potentials.bondTorsions.empty())
      unsupported = "bond-torsion";
    else if (!potentials.bendBends.empty())
      unsupported = "bend-bend";
    else if (!potentials.bendTorsions.empty())
      unsupported = "bend-torsion";
    if (unsupported)
    {
      reason = std::format("{} terms (component '{}')", unsupported, component.name);
      return false;
    }
    if (component.atoms.size() > maximumAtomsPerMolecule)
    {
      reason =
          std::format("molecules with more than {} atoms (component '{}')", maximumAtomsPerMolecule, component.name);
      return false;
    }
  }
  for (const Atom& atom : system.spanOfMoleculeAtoms())
  {
    if (atom.groupId != 0)
    {
      reason = "dU/dlambda group atoms";
      return false;
    }
  }
  return true;
}

void BondedTopology::build(const System& system)
{
  terms.clear();
  gradientOffset.clear();
  atomGradientStart.clear();
  atomGradients.clear();
  instanceMolecule.clear();
  molecules.clear();
  massOfAtom.clear();

  // the terms of every component, flattened and grouped by kind, with the gradient slots per component atom (CSR)
  std::vector<std::uint32_t> componentAtomOffset(system.components.size(), 0);
  std::vector<std::uint32_t> componentTermOffset(system.components.size(), 0);
  std::vector<std::uint32_t> componentTerms(system.components.size(), 0);
  std::vector<std::uint32_t> componentGradients(system.components.size(), 0);
  std::vector<std::vector<std::uint32_t>> referencesPerAtom;
  for (std::size_t c = 0; c < system.components.size(); ++c)
  {
    const Component& component = system.components[c];
    const Potentials::IntraMolecularPotentials& potentials = component.intraMolecularPotentials;
    componentAtomOffset[c] = static_cast<std::uint32_t>(referencesPerAtom.size());
    componentTermOffset[c] = static_cast<std::uint32_t>(terms.size());
    const std::size_t atomsInComponent = component.atoms.size();
    referencesPerAtom.resize(referencesPerAtom.size() + atomsInComponent);
    std::uint32_t gradients = 0;
    auto addTerm = [&](std::uint32_t kind, std::size_t type, std::span<const std::size_t> identifiers,
                       std::span<const double> values)
    {
      Term term{};
      term.kind = kind;
      term.type = static_cast<std::uint32_t>(type);
      for (std::size_t k = 0; k < identifiers.size(); ++k) term.atoms[k] = static_cast<std::uint32_t>(identifiers[k]);
      for (std::size_t k = 0; k < std::min<std::size_t>(values.size(), 6); ++k)
      {
        term.parameters[k] = static_cast<float>(values[k]);
      }
      terms.push_back(term);
      gradientOffset.push_back(gradients);
      for (std::uint32_t role = 0; role < identifiers.size(); ++role)
      {
        referencesPerAtom[componentAtomOffset[c] + identifiers[role]].push_back(gradients + role);
      }
      gradients += static_cast<std::uint32_t>(identifiers.size());
    };
    for (const BondPotential& bond : potentials.bonds)
    {
      addTerm(0, std::to_underlying(bond.type), bond.identifiers, bond.parameters);
    }
    for (const BendPotential& bend : potentials.bends)
    {
      addTerm(1, std::to_underlying(bend.type), bend.identifiers, bend.parameters);
    }
    for (const TorsionPotential& torsion : potentials.torsions)
    {
      addTerm(2, std::to_underlying(torsion.type), torsion.identifiers, torsion.parameters);
    }
    for (const TorsionPotential& torsion : potentials.improperTorsions)
    {
      addTerm(3, std::to_underlying(torsion.type), torsion.identifiers, torsion.parameters);
    }
    for (const VanDerWaalsPotential& pair : potentials.vanDerWaals)
    {
      const double values[2] = {pair.scaling * 4.0 * pair.parameters[0], pair.parameters[1] * pair.parameters[1]};
      addTerm(4, std::to_underlying(pair.type), pair.identifiers, values);
    }
    for (const CoulombPotential& pair : potentials.coulombs)
    {
      const double values[1] = {pair.scaling * Units::CoulombicConversionFactor * pair.chargeA * pair.chargeB};
      addTerm(5, std::to_underlying(pair.type), pair.identifiers, values);
    }
    componentTerms[c] = static_cast<std::uint32_t>(terms.size()) - componentTermOffset[c];
    componentGradients[c] = gradients;
  }
  atomGradientStart.reserve(referencesPerAtom.size() + 1);
  atomGradientStart.push_back(0);
  for (const std::vector<std::uint32_t>& references : referencesPerAtom)
  {
    atomGradients.insert(atomGradients.end(), references.begin(), references.end());
    atomGradientStart.push_back(static_cast<std::uint32_t>(atomGradients.size()));
  }

  // the molecules: their term instances are consecutive (ordered by kind within a molecule), their gradient slots
  // as well
  std::size_t instances = 0;
  numberOfGradients = 0;
  molecules.reserve(system.moleculeData.size());
  for (std::size_t m = 0; m < system.moleculeData.size(); ++m)
  {
    const Molecule& molecule = system.moleculeData[m];
    const std::size_t c = molecule.componentId;
    MoleculeInfo info{};
    info.firstAtom = static_cast<std::uint32_t>(molecule.atomIndex);
    info.numberOfAtoms = static_cast<std::uint32_t>(molecule.numberOfAtoms);
    info.atomOffset = componentAtomOffset[c];
    info.termOffset = componentTermOffset[c];
    info.instanceBase = static_cast<std::uint32_t>(instances);
    info.gradientBase = static_cast<std::uint32_t>(numberOfGradients);
    molecules.push_back(info);
    instanceMolecule.insert(instanceMolecule.end(), componentTerms[c], static_cast<std::uint32_t>(m));
    instances += componentTerms[c];
    numberOfGradients += componentGradients[c];
  }
  if (instances > std::numeric_limits<std::uint32_t>::max() || numberOfGradients > std::numeric_limits<std::uint32_t>::max())
  {
    throw std::runtime_error("[Device bonded kernel]: too many bonded term instances for 32-bit indices\n");
  }
  numberOfInstances = instances;

  const std::span<const Atom> atoms = system.spanOfMoleculeAtoms();
  massOfAtom.reserve(atoms.size());
  chargeSquaredSum = 0.0;
  for (const Atom& atom : atoms)
  {
    massOfAtom.push_back(static_cast<float>(system.forceField.pseudoAtoms[atom.type].mass));
    chargeSquaredSum += atom.charge * atom.charge;
  }
}

void BondedTopology::layout(std::span<const std::uint32_t> slotOfSorted, std::span<const std::uint32_t> originalToSorted,
                            std::size_t slots, std::vector<std::uint32_t>& slotMolecule,
                            std::vector<std::uint32_t>& referenceOfSorted) const
{
  slotMolecule.assign(slots, noAtom);
  referenceOfSorted.assign(originalToSorted.size(), 0);
  for (std::size_t m = 0; m < molecules.size(); ++m)
  {
    const MoleculeInfo& molecule = molecules[m];
    const std::uint32_t reference = originalToSorted[molecule.firstAtom];
    for (std::uint32_t b = 0; b < molecule.numberOfAtoms; ++b)
    {
      const std::uint32_t sorted = originalToSorted[molecule.firstAtom + b];
      slotMolecule[slotOfSorted[sorted]] = (static_cast<std::uint32_t>(m) << 8) | b;
      referenceOfSorted[sorted] = reference;
    }
  }
}
