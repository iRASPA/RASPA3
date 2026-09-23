module;

export module cbmc_growth_context;

import atom;
import forcefield;
import framework;
import simulationbox;
import interpolation_energy_grid;

import std;

export namespace CBMC
{
  struct GrowContext
  {
    bool hasExternalField;
    const ForceField &forceField;
    const SimulationBox &simulationBox;
    const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids;
    const std::optional<InterpolationEnergyGrid> &externalFieldInterpolationGrid;
    const std::optional<Framework> &framework;
    std::span<const Atom> frameworkAtoms;
    std::span<const Atom> moleculeAtoms;
    double beta;
    double cutOffFrameworkVDW;
    double cutOffMoleculeVDW;
    double cutOffCoulomb;
    // Scale factor of the intramolecular non-bonded (intra van der Waals + intra Coulomb) energy in
    // the trial-selection weights of chain growth and retracing. Used by NCMC reinsertion to grow
    // and retrace chains with softened intramolecular strain; 1.0 (the default, left untouched by
    // the aggregate initialization of every other caller) recovers standard CBMC.
    double intraScaling{1.0};
  };
}
