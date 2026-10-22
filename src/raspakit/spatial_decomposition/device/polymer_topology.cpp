module;

module spatial_decomposition_device_polymer_topology;

import std;

import atom;
import molecule;
import component;
import forcefield;
import vdwparameters;
import system;
import polymer_links;
import bond_potential;
import bend_potential;
import torsion_potential;
import intra_molecular_exclusions;

namespace
{
using Term = PolymerOverlayTopology::Term;
using PairState = PolymerLinkTable::PairState;
using PairScalingInfo = PolymerLinkTable::PairScalingInfo;

/// The (4 epsilon, sigma^6, shift) of a Lennard-Jones pair as the device kernels take it (the force-switched form
/// with the constant it adds below the switching distance, see DeviceStep / BondedTopology); zeros when the pair
/// potential is not a Lennard-Jones form (the overlay then has no van der Waals term for the pair).
std::array<double, 3> lennardJonesTriple(const VDWParameters& parameters)
{
  if (parameters.type != VDWParameters::Type::LennardJones &&
      parameters.type != VDWParameters::Type::LennardJonesSwitched &&
      parameters.type != VDWParameters::Type::LennardJonesForceSwitched)
  {
    return {0.0, 0.0, 0.0};
  }
  const double epsilon4 = 4.0 * parameters.parameters.x;
  const double sigma2 = parameters.parameters.y * parameters.parameters.y;
  const double sigma6 = sigma2 * sigma2 * sigma2;
  double shift = parameters.shift;
  if (parameters.type == VDWParameters::Type::LennardJonesForceSwitched)
  {
    const double rc = parameters.parameters2.x;
    const double invRc3 = 1.0 / (rc * rc * rc);
    const double invRs3 = invRc3 / parameters.parameters2.z;
    shift = epsilon4 * sigma6 * (sigma6 * invRc3 * invRc3 * invRs3 * invRs3 - invRc3 * invRs3);
  }
  return {epsilon4, sigma6, shift};
}

/// The parameters of a pair in one treatment: full (the regular table), scaled (the pair table of the scaling,
/// times its van der Waals scaling, with the Coulomb scaling in slot 3) or excluded (none).
void pairParameters(const ForceField& forceField, PairState state, const PairScalingInfo& scaling, std::size_t typeA,
                    std::size_t typeB, float* out)
{
  for (std::size_t k = 0; k < 6; ++k) out[k] = 0.0f;
  switch (state)
  {
    case PairState::Full:
    {
      const std::array<double, 3> lj = lennardJonesTriple(forceField(typeA, typeB));
      out[0] = static_cast<float>(lj[0]);
      out[1] = static_cast<float>(lj[1]);
      out[2] = static_cast<float>(lj[2]);
      out[3] = 1.0f;
      break;
    }
    case PairState::Scaled:
    {
      const std::array<double, 3> lj = lennardJonesTriple(forceField.pair(typeA, typeB, scaling.pair14));
      out[0] = static_cast<float>(scaling.scalingVDW * lj[0]);
      out[1] = static_cast<float>(lj[1]);
      out[2] = static_cast<float>(scaling.scalingVDW * lj[2]);
      out[3] = static_cast<float>(scaling.scalingCoulomb);
      break;
    }
    case PairState::Excluded:
      break;
  }
}

template <std::size_t N>
void copyParameters(const std::array<double, N>& values, float* out)
{
  for (std::size_t k = 0; k < 6; ++k) out[k] = k < N ? static_cast<float>(values[k]) : 0.0f;
}
}  // namespace

void PolymerOverlayTopology::build(const System& system, const PolymerLinkTable& table)
{
  terms.clear();
  atoms.clear();
  refStart.clear();
  refLinkEnd.clear();
  refs.clear();
  rampSteps.clear();
  numberOfRetypeCorrections = 0;
  if (!table.derivedValid)
  {
    throw std::runtime_error("[Polymerization]: the derived topology of the polymer bonds is out of date\n");
  }
  const ForceField& forceField = system.forceField;
  const std::span<const Atom> systemAtoms = system.spanOfMoleculeAtoms();

  std::map<std::uint64_t, std::uint32_t> rampOfStep{};
  auto rampIndex = [&](std::uint64_t step)
  {
    auto [it, inserted] = rampOfStep.try_emplace(step, static_cast<std::uint32_t>(rampSteps.size()));
    if (inserted) rampSteps.push_back(step);
    return it->second;
  };
  auto stepOfBond = [&](std::uint32_t bondId) { return table.bonds[bondId].formedAtStep; };
  auto add = [&](Term term) { terms.push_back(term); };

  for (const PolymerLinkTable::LinkBondTerm& link : table.linkBonds)
  {
    Term term{};
    term.kind = 0;
    term.type = static_cast<std::uint32_t>(std::to_underlying(link.potential.type));
    term.ramp = rampIndex(stepOfBond(link.bondId));
    term.atoms[0] = link.atoms[0];
    term.atoms[1] = link.atoms[1];
    term.flags = interMolecular;
    copyParameters(link.potential.parameters, term.parameters);
    add(term);
  }
  for (const PolymerLinkTable::LinkBendTerm& link : table.linkBends)
  {
    Term term{};
    term.kind = 1;
    term.type = static_cast<std::uint32_t>(std::to_underlying(link.potential.type));
    term.ramp = rampIndex(stepOfBond(link.bondId));
    for (std::size_t i = 0; i < 3; ++i) term.atoms[i] = link.atoms[i];
    term.flags = interMolecular;
    copyParameters(link.potential.parameters, term.parameters);
    add(term);
  }
  for (const PolymerLinkTable::LinkTorsionTerm& link : table.linkTorsions)
  {
    Term term{};
    term.kind = 2;
    term.type = static_cast<std::uint32_t>(std::to_underlying(link.potential.type));
    term.ramp = rampIndex(stepOfBond(link.bondId));
    for (std::size_t i = 0; i < 4; ++i) term.atoms[i] = link.atoms[i];
    term.flags = interMolecular;
    copyParameters(link.potential.parameters, term.parameters);
    add(term);
  }
  if (table.chiralForceConstant > 0.0)
  {
    for (const PolymerLinkTable::StereoCentreTerm& centre : table.stereoCentres)
    {
      Term term{};
      term.kind = 7;
      term.ramp = rampIndex(stepOfBond(centre.bondId));
      for (std::size_t i = 0; i < 4; ++i) term.atoms[i] = centre.atoms[i];
      term.flags = interMolecular;
      term.parameters[0] = static_cast<float>(table.chiralForceConstant);
      term.parameters[1] = static_cast<float>(static_cast<double>(centre.chirality) * table.chiralReferenceVolume);
      add(term);
    }
  }
  for (const PolymerLinkTable::PairCorrection& pair : table.pairCorrections)
  {
    Term term{};
    term.kind = 6;
    term.type = static_cast<std::uint32_t>(std::to_underlying(pair.wanted));
    term.type2 = static_cast<std::uint32_t>(std::to_underlying(pair.counted));
    term.ramp = rampIndex(stepOfBond(pair.bondId));
    term.atoms[0] = pair.atoms[0];
    term.atoms[1] = pair.atoms[1];
    term.flags = interMolecular;
    const std::size_t typeA = static_cast<std::size_t>(systemAtoms[pair.atoms[0]].type);
    const std::size_t typeB = static_cast<std::size_t>(systemAtoms[pair.atoms[1]].type);
    pairParameters(forceField, pair.wanted, pair.wantedScaling, typeA, typeB, term.parameters);
    pairParameters(forceField, pair.counted, pair.countedScaling, typeA, typeB, term.parameters2);
    add(term);
  }
  for (const PolymerLinkTable::OverrideBondTerm& over : table.overrideBonds)
  {
    Term term{};
    term.kind = 3;
    term.type = static_cast<std::uint32_t>(std::to_underlying(over.replacement.type));
    term.type2 = static_cast<std::uint32_t>(std::to_underlying(over.original.type));
    term.ramp = rampIndex(over.step);
    term.atoms[0] = over.atoms[0];
    term.atoms[1] = over.atoms[1];
    copyParameters(over.replacement.parameters, term.parameters);
    copyParameters(over.original.parameters, term.parameters2);
    add(term);
  }
  for (const PolymerLinkTable::OverrideBendTerm& over : table.overrideBends)
  {
    Term term{};
    term.kind = 4;
    term.type = static_cast<std::uint32_t>(std::to_underlying(over.replacement.type));
    term.type2 = static_cast<std::uint32_t>(std::to_underlying(over.original.type));
    term.ramp = rampIndex(over.step);
    for (std::size_t i = 0; i < 3; ++i) term.atoms[i] = over.atoms[i];
    copyParameters(over.replacement.parameters, term.parameters);
    copyParameters(over.original.parameters, term.parameters2);
    add(term);
  }
  for (const PolymerLinkTable::OverrideTorsionTerm& over : table.overrideTorsions)
  {
    Term term{};
    term.kind = 5;
    term.type = static_cast<std::uint32_t>(std::to_underlying(over.replacement.type));
    term.type2 = static_cast<std::uint32_t>(std::to_underlying(over.original.type));
    term.ramp = rampIndex(over.step);
    for (std::size_t i = 0; i < 4; ++i) term.atoms[i] = over.atoms[i];
    copyParameters(over.replacement.parameters, term.parameters);
    copyParameters(over.original.parameters, term.parameters2);
    add(term);
  }

  // The scaled pairs of the retyped molecules: the device bonded kernels hold the component's types, the host
  // follows the atoms (the retyping is not ramped on the host, so neither is this: the always-on ramp)
  {
    std::uint32_t fullRamp = std::numeric_limits<std::uint32_t>::max();
    for (const Molecule& molecule : system.moleculeData)
    {
      const Component& component = system.components[molecule.componentId];
      const std::span<const Atom> moleculeAtoms = systemAtoms.subspan(molecule.atomIndex, molecule.numberOfAtoms);
      bool retyped = false;
      for (std::size_t b = 0; b < molecule.numberOfAtoms && !retyped; ++b)
      {
        retyped = moleculeAtoms[b].type != component.atoms[b].type;
      }
      if (!retyped) continue;
      for (const IntraMolecularExclusions::ScaledPair& pair : component.intraMolecularPotentials.exclusions.scaledPairs)
      {
        const std::size_t currentA = static_cast<std::size_t>(moleculeAtoms[pair.atomA].type);
        const std::size_t currentB = static_cast<std::size_t>(moleculeAtoms[pair.atomB].type);
        const std::size_t originalA = static_cast<std::size_t>(component.atoms[pair.atomA].type);
        const std::size_t originalB = static_cast<std::size_t>(component.atoms[pair.atomB].type);
        if (currentA == originalA && currentB == originalB) continue;
        if (pair.scalingVDW == 0.0) continue;
        const PairScalingInfo scaling{pair.scalingVDW, pair.scalingCoulomb, pair.pair14};
        Term term{};
        term.kind = 6;
        term.type = static_cast<std::uint32_t>(std::to_underlying(PairState::Scaled));
        term.type2 = term.type;
        term.atoms[0] = static_cast<std::uint32_t>(molecule.atomIndex + pair.atomA);
        term.atoms[1] = static_cast<std::uint32_t>(molecule.atomIndex + pair.atomB);
        term.flags = intramolecularPair;
        pairParameters(forceField, PairState::Scaled, scaling, currentA, currentB, term.parameters);
        pairParameters(forceField, PairState::Scaled, scaling, originalA, originalB, term.parameters2);
        // the Coulomb parts cancel exactly (same charges, same scaling): leave them out
        term.parameters[3] = 0.0f;
        term.parameters2[3] = 0.0f;
        if (term.parameters[0] == term.parameters2[0] && term.parameters[1] == term.parameters2[1] &&
            term.parameters[2] == term.parameters2[2])
        {
          continue;
        }
        if (fullRamp == std::numeric_limits<std::uint32_t>::max()) fullRamp = rampIndex(alwaysOn);
        term.ramp = fullRamp;
        add(term);
        ++numberOfRetypeCorrections;
      }
    }
  }

  if (terms.size() * 4 > std::numeric_limits<std::uint32_t>::max())
  {
    throw std::runtime_error("[Polymerization]: too many overlay terms for 32-bit gradient indices\n");
  }

  // the gradient slots per participating atom, the inter-molecular ones first
  std::map<std::uint32_t, std::pair<std::vector<std::uint32_t>, std::vector<std::uint32_t>>> perAtom{};
  for (std::size_t t = 0; t < terms.size(); ++t)
  {
    const Term& term = terms[t];
    const std::size_t count = (term.kind == 0 || term.kind == 3 || term.kind == 6) ? 2
                              : (term.kind == 1 || term.kind == 4)                 ? 3
                                                                                   : 4;
    for (std::size_t role = 0; role < count; ++role)
    {
      auto& lists = perAtom[term.atoms[role]];
      const std::uint32_t slot = static_cast<std::uint32_t>(4 * t + role);
      if ((term.flags & interMolecular) != 0)
        lists.first.push_back(slot);
      else
        lists.second.push_back(slot);
    }
  }
  atoms.reserve(perAtom.size());
  refStart.reserve(perAtom.size() + 1);
  refLinkEnd.reserve(perAtom.size());
  refStart.push_back(0);
  for (const auto& [atom, lists] : perAtom)
  {
    atoms.push_back(atom);
    refs.insert(refs.end(), lists.first.begin(), lists.first.end());
    refLinkEnd.push_back(static_cast<std::uint32_t>(refs.size()));
    refs.insert(refs.end(), lists.second.begin(), lists.second.end());
    refStart.push_back(static_cast<std::uint32_t>(refs.size()));
  }
}

void PolymerOverlayTopology::lambdas(const PolymerLinkTable& table, std::vector<float>& out) const
{
  out.resize(rampSteps.size());
  for (std::size_t k = 0; k < rampSteps.size(); ++k)
  {
    out[k] = rampSteps[k] == alwaysOn ? 1.0f : static_cast<float>(table.rampFactor(rampSteps[k]));
  }
}
