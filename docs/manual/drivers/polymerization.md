# Polymerization
\page driver_polymerization Polymerization

`"SimulationType" : "Polymerization"`

The polymerization driver grows a polymer, or cures a network, during a
molecular-dynamics run: reactive groups that come close to each other are
bonded according to reaction rules given in the input. It is built on
[spatial-decomposition molecular dynamics](\ref driver_molecular_dynamics_spatial_decomposition),
all of whose settings apply, because the systems of interest are large boxes
of monomers.

## What kind of model this is

Classical force fields have a fixed set of bonds, so a chemical reaction
cannot happen by itself. A quantum-mechanical description of bond formation
is far too expensive for the tens of thousands of atoms and the long times of
a curing process. *Reactive molecular dynamics with templates* takes a middle
road: the dynamics is classical, and at intervals the bonding topology is
changed by rule. The rules do not describe *how* a bond forms. They describe
*when* one may form (a distance criterion, and a probability) and *what the
force field is afterwards*.

The model therefore gives the structure that results from the reactions: the
chain-length distribution, the network topology, the gel point, the density
and shrinkage of the cured material, and from the cured structure its
mechanical and thermal properties. The rate constants are input, not output.

## Reactive groups and states

A monomer carries one or more *reactive groups*. A group has a type (for
instance an acrylate, an epoxide or an amine), a set of named atoms that play
a role in the reactions, and a *state*, which starts as `monomer`. A reaction
is a transition between states:

- **unimolecular**, of one group: an initiation, monomer \f$\to\f$ radical;
- **bimolecular**, of a pair of groups whose capture atoms are within a
  radius: a propagation (a radical bonds to a monomer, which becomes the new
  radical), or a termination by combination (two radicals bond).

Chain growth and step growth are written in the same language. A free-radical
polymerization is the set initiation, propagation, termination. An
epoxy-amine cure is an epoxide reacting with a primary amine, which becomes a
secondary amine and can react once more.

## The reaction step

Every `"ScanEvery"` steps of the dynamics the groups are scanned.

1. For every bimolecular rule, the pairs of groups in its reactant states
   whose capture atoms are within the capture radius are candidates.
2. The candidates of all rules are pooled and tried in order of increasing
   distance, at most one reaction per group per scan. Each is accepted with
   the probability of its rule.
3. Every group in the reactant state of a unimolecular rule reacts with the
   probability of that rule.

The probability can be given directly or as a rate constant \f$k\f$. A
first-order process with rate \f$k\f$ has a probability of having happened
within the scan interval \f$\Delta t_\text{scan}\f$ of

\f[
p = 1 - e^{-k\,\Delta t_\text{scan}} ,
\f]

so that a rate given this way does not change when the scan interval does.

**The time scale.** As in every method of this kind, real polymerizations
take seconds to hours while a simulation covers nanoseconds. The reaction
probabilities that are used in practice are therefore far higher than the
physical ones, and the simulated time is not the physical reaction time. The
sensitivity of the final structure to the probabilities and the scan interval
is something to test.

## The orientation weight

A distance criterion alone accepts a pair whatever the relative orientation
of the two groups, and then pulls the new bond into its equilibrium geometry.
A rule can be given an orientation temperature \f$T^*\f$; the acceptance
probability of a candidate is then multiplied by

\f[
\exp\!\left(-\frac{\Delta U_\text{ang}}{k_B T^*}\right) ,
\f]

where \f$\Delta U_\text{ang}\f$ is the angular strain the new bond would have
in the current geometry: the sum of the bend and torsion energies through
the new bond, each relative to its minimum. A pair that is already close to
the product geometry is favoured. \f$T^*\f$ is a selectivity parameter, not
the physical temperature; the strain of a random approach is thousands of
kelvin, so useful values are of that order.

## Changing the force field

An accepted reaction does three things.

- **Retyping.** The atoms named in the rule get new types, and with them new
  charges and Lennard-Jones parameters: an \f$sp^2\f$ carbon of the monomer
  becomes a backbone carbon of the polymer.
- **Bonds.** The bonds of the rule are formed, and bonds may be broken (the
  ring bond of an epoxide that opens).
- **Bonded terms from type-based tables.** The bond, the bends and the
  torsions through a new bond are looked up, by the *current* types of the
  atoms, in tables of the product. Terms of the monomer that contain a
  retyped atom are replaced by the entry for the new types. The non-bonded
  interaction of pairs that have become 1-2 or 1-3 neighbours is removed, and
  1-4 pairs are scaled.

Because the terms are found by type, the tables are those of an ordinary
molecule file of the polymer. No template of the full reacting neighbourhood
has to be drawn by hand, and a missing entry is reported with the tuple of
types to add.

The molecules of the system never change. The polymer is the graph of bonds
between them, and that graph is the result: it is written at the end, and
every reaction is logged.

## Relaxation and energy

At the moment of reaction the two atoms are up to a capture radius apart,
much further than a bond length. Switching the bond on at once would give an
enormous force. All new or changed terms are therefore switched on linearly,
with a factor that rises from 0 to 1 over `"RelaxationSteps"` steps.

During the ramp the Hamiltonian depends on time, so the conserved quantity of
the dynamics is not conserved: it changes by the work of the ramp, which is
the heat of reaction of the model and is large (of the order of
\f$10^4\f$-\f$10^5\f$ K per bond). The thermostat removes it. The output
keeps the two apart: the change of the conserved energy during ramps, plus
the jump of the potential at the reactions, is reported as ramp work, so that
the remaining drift still measures the quality of the integration.

## Kinetics that can be read off

- **Conversion**: the fraction of reactive groups that have reacted. The
  reactions stop at the target conversion; the dynamics continues, so that
  the cured structure equilibrates.
- **Active groups**: groups in a state that only a reaction can produce, the
  living chain ends. Their number is reported, and the distribution of their
  waiting times between reactions is collected.
- **The topology**: chain lengths, cycles, cross-link density, the gel point.

## Tacticity

A reaction that creates a stereocentre can be given a rule for its
configuration, with a probability of a meso placement relative to the
previous stereocentre of the chain, so that isotactic, syndiotactic and
atactic chains can be grown.

## Stages

The initialization (Monte Carlo) and the equilibration (molecular dynamics)
run on the unreacted monomers. The reactions are active during production
only.

## Limits

- The reactions are irreversible and their rates are input.
- The method is separate from the equilibrium cross-link Monte Carlo moves
  and cannot be combined with them.
- The limits of
  [spatial-decomposition molecular dynamics](\ref driver_molecular_dynamics_spatial_decomposition)
  apply: molecules in a box, without framework.

## References

- J.R. Gissinger, B.D. Jensen and K.E. Wise, *Modeling chemical reactions in
  classical molecular dynamics simulations*, Polymer **128**, 211 (2017).
- J.R. Gissinger, B.D. Jensen and K.E. Wise, *REACTER: A Heuristic Method for
  Reactive Molecular Dynamics*, Macromolecules **53**, 9953 (2020).
- V. Varshney, S.S. Patnaik, A.K. Roy and B.L. Farmer, *A Molecular Dynamics
  Study of Epoxy-Based Networks: Cross-Linking Procedure and Prediction of
  Molecular and Material Properties*, Macromolecules **41**, 6837 (2008).
- P.J. Flory, *Principles of Polymer Chemistry*, Cornell University Press
  (1953).
