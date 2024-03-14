export module cbmc_multiple_first_bead;

import <optional>;
import <span>;

import atom;
import randomnumbers;
import cbmc_first_bead_data;
import forcefield;
import simulationbox;

export namespace CBMC
{
  [[nodiscard]] std::optional<FirstBeadData> 
  growMoleculeMultipleFirstBeadSwapInsertion(RandomNumber &random, bool hasExternalField, const ForceField &forceField, 
                                             const SimulationBox &simulationBox, std::span<const Atom> frameworkAtoms, 
                                             std::span<const Atom> moleculeAtoms, double beta, double cutOff, 
                                             double cutOffCoulomb, const Atom& atom, size_t numberOfTrialDirections) 
                                             noexcept;

  [[nodiscard]] FirstBeadData 
  retraceRigidMultipleFirstBeadSwapDeletion(RandomNumber &random, bool hasExternalField, const ForceField &forcefield, 
                                            const SimulationBox &simulationBox, std::span<const Atom> frameworkAtoms, 
                                            std::span<const Atom> moleculeAtoms, double beta, double cutOff, 
                                            double cutOffCoulomb, const Atom& atom, double scaling, double storedR, 
                                            size_t numberOfTrialDirections) noexcept;

  [[nodiscard]] std::optional<FirstBeadData> 
  growRigidMultipleFirstBeadReinsertion(RandomNumber &random, bool hasExternalField, const ForceField &forceField, 
                                        const SimulationBox &simulationBox, std::span<const Atom> frameworkAtoms, 
                                        std::span<const Atom> moleculeAtoms, double beta, double cutOff, 
                                        double cutOffCoulomb, const Atom& atom, size_t numberOfTrialDirections) 
                                        noexcept;

  [[nodiscard]] FirstBeadData 
  retraceRigidMultipleFirstBeadReinsertion(RandomNumber &random, bool hasExternalField, const ForceField &forceField, 
                                           const SimulationBox &simulationBox, std::span<const Atom> frameworkAtoms, 
                                           std::span<const Atom> moleculeAtoms, double beta, double cutOff, 
                                           double cutOffCoulomb, const Atom& atom, double storedR, 
                                           size_t numberOfTrialDirections) noexcept;
}
