# Polymers
\page polymers Polymers

A polymer is a molecule made of many repeat units, and a network is one
molecule that fills the box. `RASPA3` handles them in three ways, which this
page describes together because they share their concepts:

| What | How it is represented | Method | Section |
|---|---|---|---|
| A chain of given length | One molecule, defined in a molecule file | Monte Carlo with moves made for chains, or molecular dynamics | [Chains](#polymers-chains) |
| A reversible network | Molecules plus a table of bonds between them, at equilibrium | Monte Carlo moves that form, break and rewire the bonds | [Cross-links](#polymers-cross-links) |
| A polymerization or a cure | Monomers plus a table of the bonds that reactions formed, irreversible | Molecular dynamics with reaction rules | [Polymerization](#polymers-polymerization) |

In the last two the polymer is not a molecule. The molecules stay what their
molecule files say they are, and the bonds between them are kept in an
[overlay](#polymers-overlay). The largest part of this page is the
[reaction input](#polymers-reaction-input) of the polymerization, explained
line by line for the curing of hexanediol diacrylate (HDDA).

## Chains as molecules {#polymers-chains}

A chain of known length is an ordinary flexible molecule: the molecule file
lists its atoms, its connectivity and its [bonded terms](\ref potentials).
The file `pHDDA-20.json` of the examples, a chain of twenty HDDA units, is of
this kind. Three keys of the molecule file exist for chains:

| Key | Meaning |
|---|---|
| `"EndToEndAtoms" : [i, j]` | The two ends of the chain. They define the backbone (the path between them) for the moves and the end-to-end vector for the properties. Default: the two atoms farthest apart along the bonds. |
| `"RepeatUnits" : [[...], [...]]` | The atoms of every repeat unit, in order along the chain. Needed by the reptation move. |
| `"ChiralCenters"` | The stereocentres of the chain. The moves that change the internal geometry reject a trial that inverts one. |

### Moves for chains {#polymers-chain-moves}

A long chain in a dense system cannot be sampled by regrowing it: the
probability that a new conformation fits decreases exponentially with its
length. The moves below change a chain locally, or change which atoms are
connected, and are set per component with the keys of the
[commands](\ref commands) page.

| Move | Key | What moves | What is conserved |
|---|---|---|---|
| Bead displacement | `"BeadDisplacementProbability"` | One bead, along one axis | Nothing: it is the move that relaxes flexible bond lengths and bend angles |
| Bead flip | `"BeadFlipProbability"` | One bead, rotated about the axis through its two neighbours; a chain end about its neighbour | Bond lengths |
| Crankshaft | `"CrankshaftProbability"` | A segment of a few atoms between two anchors, rotated about the axis through the anchors | Bond lengths; the chain ends stay in place |
| Pivot | `"PivotProbability"`, `"PivotCBMCProbability"` | The smaller of the two parts on either side of a bond, rotated about that bond | Bond lengths and bend angles |
| Concerted rotation | `"ConcertedRotationProbability"` | Four atoms inside a window of eight backbone atoms, by solving a loop closure | Bond lengths and bend angles; everything outside the window |
| Reptation | `"ReptationProbability"` | A repeat unit is removed at one end and grown at the other | The chain length |
| Double bridging | `"DoubleBridgingProbability"` | Two chains exchange their tails: a trimer is cut out of each and the head of one is bridged to the tail of the other | Every chain length, bond lengths and bend angles |
| Intramolecular double rebridging | `"IntramolecularDoubleRebridgingProbability"` | The segment between two sites of one chain is reversed | Chain length, bond lengths and bend angles |
| Reinsertion, partial reinsertion | `"ReinsertionProbability"`, `"PartialReinsertionProbability"` | The whole chain, or a part of it, regrown with configurational bias or recoil growth | |

How they complement each other:

- The **local moves** (bead displacement, bead flip, crankshaft, concerted
  rotation) are accepted in a dense melt because they move few atoms.
  Concerted rotation is the one of choice for long chains: it is strictly
  local and keeps all bond lengths and angles. It does not move the first and
  last two backbone atoms, so it needs a move for the chain ends next to it.
- The **pivot** changes the shape of a whole chain in one step. It is very
  efficient for an isolated chain or a dilute solution, and fails in a melt,
  where the rotated part overlaps with its surroundings.
- The **connectivity-altering moves** (double bridging, double rebridging)
  move few atoms but change which atoms belong to which chain, so the
  end-to-end vectors of two chains change at once. They decorrelate the
  large-scale conformation of a melt orders of magnitude faster than local
  moves.
- **Reptation** moves a chain along its own contour, and is effective as
  long as a chain end can find room.

The concerted rotation and the bridging moves solve a geometric closure
problem, which can have several solutions or none. Their acceptance rules
contain the ratio of the numbers of solutions and the ratio of the Jacobians
of the closure, besides the Boltzmann factor.

The conformation of the chains is measured with the properties described on
the [properties](\ref properties) page: the end-to-end distance and its
autocorrelation function, whose relaxation time tells how far apart two
independent samples are, and the gyration tensor.

## The overlay {#polymers-overlay}

A network that spans the box is, chemically, one molecule with as many atoms
as the system. Representing it as one would make every data structure that
assumes small, identical molecules useless, and every reaction would have to
merge two molecule definitions into a new one.

`RASPA3` does not do that. **The molecules never change.** Each keeps the
atoms, the connectivity and the bonded terms of its component. The bonds
between molecules are kept in a table owned by the system, and that table is
turned into a set of extra energy terms: the *overlay*.

The energy is evaluated in two layers. The regular routines compute the
energy as if no bond between molecules existed: the bonded terms within each
molecule, and the non-bonded interaction of all pairs of atoms of different
molecules. The overlay then adds what the bonds change,

\f[
U = U_\text{molecules}
+ \sum_{\text{links}} \left[ u_\text{bond} + \sum u_\text{bend} + \sum u_\text{torsion} \right]
- \sum_{\text{pairs across a link}} \left[ u_\text{counted}(r_{ij}) - u_\text{wanted}(r_{ij}) \right] .
\f]

The last sum is the bookkeeping that makes two linked molecules have exactly
the energy of one molecule with that bond. Two atoms of different molecules
interact at full strength in the regular sums. Once a bond connects them,
they are a 1-2 pair, and the atoms next to them are 1-3 and 1-4 pairs, which
a [force field excludes or scales](\ref potential_intramolecular_nonbonded).
The pair routines are left untouched; the overlay subtracts exactly what they
counted (\f$u_\text{counted}\f$) and adds what the topology requires
(\f$u_\text{wanted}\f$: nothing for an excluded pair, the scaled interaction
for a 1-4 pair). With the Ewald summation the correction includes the
exclusion term of the Fourier part, as for a pair within a molecule.

What this buys:

- The caches and tables per component stay valid, and all molecules of a
  component remain identical in layout.
- A percolating network needs no giant molecule.
- Forming or breaking a bond is an edit of a table, followed by a rebuild of
  the overlay's term lists.
- The graph of bonds *is* the result of a polymerization, and it is available
  as such: chain lengths, clusters, cycles and the gel point are properties
  of the table.

There are two overlays. They share the idea and differ in what a bond brings
with it:

| | Cross-links | Polymer bonds |
|---|---|---|
| Table | `CrossLinkTable` | `PolymerLinkTable` |
| Used by | Monte Carlo (and molecular dynamics) at equilibrium | The polymerization driver |
| Bonds are | Reversible: formed, broken and rewired by moves | Irreversible: formed by reaction rules |
| A bond brings | The bond, optional bends at the junction, 1-2 and 1-3 exclusions, a formation energy | The complete bonded topology of a real bond: bond, bends, torsions, 1-2 and 1-3 exclusions, scaled 1-4 pairs |
| Atom types | Unchanged | Reacted atoms change type |
| Bonds of the molecule | Unchanged | Can be broken (ring opening) |
| Switched on | At once | With a ramp over a number of steps |

The two cannot be combined in one system. Both report their bonded energy in
the `cross-link` slot of the energy output.

## Reversible cross-links {#polymers-cross-links}

Cross-links are bonds between molecules that exist at equilibrium: they form
and break, and the simulation samples the distribution over topologies. They
describe associating polymers, reversible gels, and vitrimers, in which bonds
exchange partners while their number stays constant.

### Input

A molecule declares which of its atoms can carry a link, in its molecule
file:

```json
"ReactiveSites" : [
  [0, "X"],
  {"Atom" : 3, "Type" : "X", "Valence" : 2, "LinkBend" : ["HARMONIC", [500.0, 120.0]]}
]
```

An entry is `[atom, "type"]`, `[atom, "type", valence]`, or an object. The
*type* names the chemistry of the site. The *valence* is the largest number
of links the site may carry (default 1). A site with a valence above one can
be given a `"LinkBend"`, a bend on the angle between two of its links.

The system declares which site types can bond, and with what potential:

```json
"CrossLinkBonds" : [
  {
    "Sites" : ["X", "X"],
    "Bond" : ["HARMONIC", [1000.0, 3.8]],
    "JunctionBend" : ["HARMONIC", [500.0, 120.0]],
    "CaptureRadius" : 6.5,
    "FormationEnergy" : -1000.0
  }
]
```

| Key | Meaning |
|---|---|
| `"Sites"` | The two site types the bond joins. |
| `"Bond"` | The [bond potential](\ref potential_bonds) of the link. |
| `"JunctionBend"` | Optional [bend](\ref potential_bends) on every angle neighbour-site-partner across the link. |
| `"CaptureRadius"` | The largest distance [Å] at which the moves propose a link. It affects the efficiency, not the distribution that is sampled. |
| `"FormationEnergy"` | A constant energy [K] per link. A negative value favours bonding; it sets the association constant and thereby the degree of cross-linking. |

`"InitialCrossLinks"` starts the simulation from given links. The links are
part of the restart file.

### The energy of a link

\f[
u_\text{link} = u_\text{bond}(r_{AB}) + \epsilon_\text{form}
+ \sum_{n_A} u_\text{junction}(\theta_{n_A A B})
+ \sum_{n_B} u_\text{junction}(\theta_{A B n_B})
- \sum_{\text{1-2, 1-3 pairs}} u_\text{pair}(r_{ij}) ,
\f]

with \f$n_A\f$ the neighbours of site \f$A\f$ within its own molecule. The
excluded pairs are the two sites themselves and the 1-3 pairs across the
link. A pair that several links reach is excluded once.

### The topology moves {#polymers-topology-moves}

Three moves change the table. They are set per system.

**Formation and scission** (`"CrossLinkFormationProbability"`). Half of the
attempts form a link, half break one. Formation picks a free site \f$i\f$
among the \f$N_F\f$ free sites and a partner \f$j\f$ among its \f$n_i\f$
candidates: free sites of a compatible type on another molecule, within the
capture radius. The same pair could have been proposed starting from
\f$j\f$, so with \f$N_b\f$ the number of links,

\f[
P_\text{acc}(\text{form}) = \min\left[1,\;
\frac{N_F}{(N_b + 1)\left(\frac{1}{n_i} + \frac{1}{n_j}\right)}\,
e^{-\beta\,\Delta U}\right] .
\f]

Scission picks one of the \f$N_b\f$ links. With the counts \f$N_F'\f$,
\f$n_i'\f$, \f$n_j'\f$ of the state without the link,

\f[
P_\text{acc}(\text{break}) = \min\left[1,\;
\frac{N_b\left(\frac{1}{n_i'} + \frac{1}{n_j'}\right)}{N_F'}\,
e^{-\beta\,\Delta U}\right] .
\f]

A link whose sites are farther apart than the capture radius is never broken,
because the formation move could not propose it. \f$\Delta U\f$ contains the
formation energy, so this pair of moves samples the equilibrium number of
links.

**Swap** (`"CrossLinkSwapProbability"`). A link \f$(i, j)\f$ and one of its
ends, the pivot \f$i\f$, are picked; the other end moves from \f$j\f$ to a
free candidate \f$k\f$ of \f$i\f$. The number of links is conserved. With
the link counted as absent the set of candidates is the same before and
after, so the proposal is symmetric and the move is accepted with

\f[
P_\text{acc} = \min\left[1,\; e^{-\beta\,\Delta U}\right] .
\f]

**Exchange** (`"CrossLinkExchangeProbability"`). Two links trade partners,

\f[
(a, b) + (c, d) \;\rightarrow\; (a, c) + (b, d) .
\f]

A link and its pivot \f$a\f$ are picked; \f$c\f$ is drawn among the
*linked* sites near \f$a\f$, and \f$(c, d)\f$ among the links of \f$c\f$.
Neither the number of links nor the number of links of any site changes.
With \f$n_\text{links}\f$ the number of links a site carries,

\f[
P_\text{acc} = \min\left[1,\;
\frac{n_\text{links}(c)}{n_\text{links}(b)}\,
e^{-\beta\left[u(a,c) + u(b,d) - u(a,b) - u(c,d)\right]}\right] ,
\f]

which is the plain Metropolis rule for sites of valence one. This is the only
move that rewires a network whose sites are all occupied, where the swap
finds no free site: the exchange reaction of a vitrimer.

### Moving linked molecules

A linked molecule is still translated, rotated and changed in shape by the
ordinary moves; the change of the link energy is part of the acceptance.
Reinsertion and partial reinsertion regrow a linked molecule with its linked
atoms kept in place. The partner side of every link is frozen into a
*tether*: a copy of the partner's site and of its neighbours, against which
the bond, the junction bends and the exclusions of the link are evaluated for
each trial position. These enter the Rosenbluth weights, so the geometry of
the junction is sampled correctly. A part of the molecule between two linked
sites is regrown with fixed end points.

A molecule that carries a link cannot be deleted; its links have to be broken
first. Components with reactive sites cannot use the fractional-molecule
moves, identity changes, pair and group swaps, reptation, double bridging,
the Gibbs swap moves or the reaction-ensemble moves, and cannot be combined
with polarization.

### The table {#polymers-cross-link-table}

`CrossLinkTable` holds

| Member | Content |
|---|---|
| `bondTypes` | The entries of `"CrossLinkBonds"` |
| `links` | The links, each two sites (component, molecule, atom) and the index of its bond type |
| `linksPerMolecule` | Derived: per molecule the links it takes part in |
| `linkCountPerSite` | Derived: per site the number of links it carries, compared with the valence |

Only the first two are stored in a restart file; the indices are rebuilt.
Because a site refers to its molecule by index, the table is renumbered when
a molecule is inserted, deleted or moved in the list. In parallel tempering
the links are exchanged along with the configurations.

## Polymerization {#polymers-polymerization}

The polymerization driver forms bonds during a molecular-dynamics run,
according to reaction rules. The method, the reaction step, the orientation
weight and the ramp are explained on the page of the
[polymerization driver](\ref driver_polymerization). Here follow the data
structure and the input.

### The polymer overlay {#polymers-link-table}

`PolymerLinkTable` has two parts. The **authoritative state** is the record
of what happened, and is what a restart file stores:

| Member | Content |
|---|---|
| `scheme` | The reaction scheme: the group types with their roles and states, and the rules |
| `templates` | The bond, bend and torsion parameters of the product, by pseudo-atom types |
| `bonds` | The bonds formed: the two atoms (component, molecule, atom), the rule, the step, the distance at formation |
| `brokenBonds` | The bonds of a molecule that a reaction broke |
| `retypings` | The log of type changes: atom, type before, type after, step |
| `centres` | The stereocentres created: the centre, its three substituents, the configuration held |
| `groupState` and related arrays | Per reactive group: its state, the step of its last change, the configuration its next stereocentre must take, the group before it in the chain |
| `radicals`, `groupRadical`, `groupLastRadical` | The life of every radical: where and when it was born, the group it sits on, when it died, the number of propagations |

The current types of the atoms are stored with the atoms themselves.

The **derived term lists** are what the energy and force routines evaluate.
They are rebuilt from the state after every reaction and after a restart:

| List | Terms | Energy added |
|---|---|---|
| `linkBonds`, `linkBends`, `linkTorsions` | The new bond \f$A\!-\!B\f$, the bends \f$X\!-\!A\!-\!B\f$ and \f$A\!-\!B\!-\!Y\f$, and the torsions \f$W\!-\!X\!-\!A\!-\!B\f$, \f$X\!-\!A\!-\!B\!-\!Y\f$ and \f$A\!-\!B\!-\!Y\!-\!Z\f$ through it | \f$\lambda\,u\f$ |
| `overrideBonds`, `overrideBends`, `overrideTorsions` | Terms of a molecule that contain a retyped atom, or that run through a broken bond | \f$\lambda\,(u_\text{replacement} - u_\text{original})\f$ |
| `pairCorrections` | Pairs of atoms whose non-bonded treatment changed | \f$\lambda\,(u_\text{wanted} - u_\text{counted})\f$ |
| `stereoCentres` | The restraints that hold the stereocentres | \f$\lambda\,\frac{1}{2}k\,(\hat V - \chi \hat V_0)^2\f$ |

\f$\lambda\f$ is the ramp factor of the term. It rises linearly from 0 to 1
over `"RelaxationSteps"` steps after the reaction that created the term, so a
bond that was captured at 4 Å pulls its atoms together gradually.

Some details of the bookkeeping:

- **Overrides.** The regular routines still evaluate the molecule with the
  parameters of its component. A term with a retyped atom is corrected by
  adding the difference to the product parameters. A term through a broken
  bond gets the zero potential as replacement, so the difference removes it.
  A term that is retyped twice (the second reaction of an amine) goes from
  the component's parameters to those of the first product with the ramp of
  the first reaction, and on to the second with its own ramp.
- **Pair corrections.** Each pair has a state *counted* (what the regular
  sums do: full, scaled or excluded) and a state *wanted* (what the topology
  requires). The correction is the difference, with the 1-4 scaling of the
  component.
- **Unlisted pairs.** The Lennard-Jones energy of two bonded atoms at a bond
  length is enormous. Computing it in the pair sum and subtracting it again
  in the overlay would leave a rounding error that is a net force, in
  particular in the single precision of a GPU. The pairs the topology
  excludes are therefore removed from the neighbour lists (`unlisted`), and
  the overlay supplies the fading remainder \f$(1 - \lambda)\,u_\text{counted}\f$
  during the ramp.
- **Neighbours.** `linkNeighbours` and `brokenPairs` give the current
  connectivity: that of the component, without the broken bonds, plus the
  polymer bonds. The bends, torsions and exclusions of a new bond are found
  by walking this graph, which is why a later reaction next to an earlier one
  gets the right terms.

### The reaction input, with HDDA as example {#polymers-reaction-input}

The example is `examples/polymers/8_polymerization_hdda_in_box`: 2690
molecules of hexanediol diacrylate are cured by free-radical polymerization
into a network, following A. Torres-Knoop *et al.* (2018).

The input has three parts, in three files:

| Part | File | Says |
|---|---|---|
| The reactive groups | the molecule file `HDDA.json` | *where* a molecule can react: which atoms play which role |
| The pseudo-atoms | `force_field.json` | the types of the atoms before and after reaction |
| The reaction rules and product tables | the `"Polymerization"` block of the system in `simulation.json` | *what* happens, *when*, and the force field afterwards |

A molecule file contains no chemistry. It only names atoms.

#### 1. The monomer and its reactive groups

HDDA is modelled with 16 united-atom sites. Numbered as in the molecule file:

```
   14 O                                              O 15
      ‖                                              ‖
 0    1    2    3    4    5    6    7    8    9   10   11   12   13
CH2 = CH - C  - O  - CH2- CH2- CH2- CH2- CH2- CH2- O  - C  - CH = CH2
head tail side                                         side tail head
```

Atoms 0 and 13 have the type `CH2_sp2`, atoms 1 and 12 `CH_sp2`, atoms 2 and
11 `C_CO`. The two double bonds are the two reactive groups:

```json
"ReactiveGroups" : [
  {"Type" : "acrylate", "State" : "monomer", "Atoms" : {"head" : 0,  "tail" : 1,  "side" : 2}},
  {"Type" : "acrylate", "State" : "monomer", "Atoms" : {"head" : 13, "tail" : 12, "side" : 11}}
]
```

| Key | Meaning |
|---|---|
| `"Type"` | The kind of group. The rules refer to it. All groups of a type must have the same roles. |
| `"State"` | The state the group starts in (default `"monomer"`). |
| `"Atoms"` | The *roles*: names for the atoms of the group that the rules need, each with the index of the atom in the molecule. |

The roles of an acrylate group are chosen after what the atoms do in the
reaction. The `head` is the \f$\text{CH}_2\!=\f$ carbon that a radical
attacks. The `tail` is the \f$=\!\text{CH}-\f$ carbon that carries the
radical afterwards. The `side` is the ester carbon on the tail; no reaction
changes it, but it is needed to define the configuration of the tail when
that becomes a stereocentre.

Because a molecule has two groups, it can become part of two chains. That is
what makes HDDA a cross-linker and the product a network instead of a set of
linear chains.

#### 2. The pseudo-atoms of the product

When the double bond reacts, its two \f$sp^2\f$ carbons become \f$sp^3\f$
carbons of the backbone. In a united-atom force field those are different
pseudo-atoms, so `force_field.json` defines both:

| Atom | Before | \f$\epsilon/k_B\f$ [K], \f$\sigma\f$ [Å] | After | \f$\epsilon/k_B\f$ [K], \f$\sigma\f$ [Å] |
|---|---|---|---|---|
| head | `CH2_sp2` | 85.0, 3.675 | `CH2_bb` | 46.0, 3.95 |
| tail | `CH_sp2` | 52.0, 3.71 | `CH_bb` | 10.0, 4.68 |

No molecule contains `CH2_bb` or `CH_bb` at the start. They appear through
the `"Retype"` of the rules.

#### 3. States

The rules are transitions between the states of a group. The example uses
four states of the type `acrylate`:

| State | The group is |
|---|---|
| `monomer` | an unreacted double bond |
| `radical_tail` | a chain end with the radical on the `tail` carbon, \f$-\text{CH}^\bullet-\text{COOR}\f$: the normal, secondary radical |
| `radical_head` | a chain end with the radical on the `head` carbon, \f$-\text{CH}_2^\bullet\f$: the primary radical left by a mis-addition |
| `dead` | reacted, part of a chain, no radical |

States are not declared anywhere: a state exists because a group starts in
it or a rule names it. A state is *active* when a bimolecular rule has it as
reactant and no group starts in it. Here those are the two radical states.
The groups in active states are the living chain ends; they are counted, and
each carries a radical with an identity that is followed through the
reactions.

#### 4. The anatomy of a rule

The propagation rule of the example, key by key:

```json
{
  "Name" : "propagation",
  "Reactants" : [{"Label" : "R", "Type" : "acrylate", "State" : "radical_tail"},
                 {"Label" : "M", "Type" : "acrylate", "State" : "monomer"}],
  "Capture" : {"Atoms" : ["R.tail", "M.head"], "Radius" : 4.0},
  "FormBonds" : [["R.tail", "M.head"]],
  "Products" : [{"Label" : "R", "State" : "dead"},
                {"Label" : "M", "State" : "radical_tail",
                 "Retype" : {"head" : "CH2_bb", "tail" : "CH_bb"}}],
  "Stereo" : {"Centre" : "R.tail", "Substituents" : ["R.head", "M.head", "R.side"], "Next" : "M"},
  "Probability" : 0.5,
  "AllowIntramolecular" : true
}
```

| Key | In this rule | In general |
|---|---|---|
| `"Name"` | `propagation` | Appears in the logs and the status report. |
| `"Reactants"` | A group `R` that is a tail radical, and a group `M` that is a monomer | One reactant makes a unimolecular rule, two a bimolecular one. Each gives the group type and the state a group must be in. The *label* is a name for that reactant within the rule (default `A` and `B`). |
| `"Capture"` | The distance between the tail of `R` and the head of `M` must be below 4 Å | The two atoms whose distance decides whether a pair is a candidate, written `Label.role`. Required for a bimolecular rule. Radius default 4 Å. |
| `"FormBonds"` | A bond between the tail of `R` and the head of `M` | The bonds formed, each between the two reactants. |
| `"Products"` | `R` becomes `dead`. `M` becomes `radical_tail`, and its head and tail become backbone atoms | Per reactant, matched by label: the new state, and `"Retype"`, the roles whose atom changes pseudo-atom type. A reactant without an entry is unchanged. |
| `"Stereo"` | The tail of `R` becomes a stereocentre | See [below](#polymers-stereo). |
| `"Probability"` | Half of the candidates react at each scan | The probability per candidate per scan (default 1). Alternatively `"RateConstant"` in ps⁻¹. |
| `"AllowIntramolecular"` | The two groups of one molecule may react with each other | With `false`, they may not (no primary cycles). Default `true`. |

Further keys that this rule does not use:

| Key | Meaning |
|---|---|
| `"BreakBonds" : [["E.carbon", "E.oxygen"]]` | Bonds of the molecule, between two roles of one reactant, that the reaction breaks: the ring opening of an epoxide. |
| `"OrientationTemperature" : T` | Multiplies the probability by \f$\exp(-\Delta U_\text{ang}/k_B T)\f$, with \f$\Delta U_\text{ang}\f$ the angular strain the new bond would have. See the [driver page](\ref driver_polymerization). |
| `"InitialEvents" : n` | For a unimolecular rule: applied to \f$n\f$ random groups when the reactions start. |

What happens when the rule fires on a pair is the chemistry

\f[
\sim\!\text{CH}_2\!-\!\text{CH}^\bullet(\text{COOR}) \;+\; \text{CH}_2\!=\!\text{CH}(\text{COOR})
\;\longrightarrow\;
\sim\!\text{CH}_2\!-\!\text{CH}(\text{COOR})\!-\!\text{CH}_2\!-\!\text{CH}^\bullet(\text{COOR}) .
\f]

The radical has moved from `R` to `M`. In the bookkeeping, the radical of `R`
is handed to `M`, because `R` leaves an active state and `M` enters one in
the same reaction. `M` is then a `radical_tail` and is the `R` of the next
propagation: the chain grows.

Note that the retyping is on `M` only. `R` was already retyped by the
reaction that made it a radical.

#### 5. The eight rules of the example

**Initiation.**

```json
{
  "Name" : "initiation",
  "Reactants" : [{"Label" : "M", "Type" : "acrylate", "State" : "monomer"}],
  "Products" : [{"Label" : "M", "State" : "radical_tail", "Retype" : {"head" : "CH2_bb", "tail" : "CH_bb"}}],
  "Probability" : 0.0,
  "InitialEvents" : 135
}
```

A unimolecular rule: a monomer group becomes a tail radical, and its two
carbons are retyped. The initiator is not modelled as a molecule; the group
is treated as if an initiator fragment (an implicit hydrogen in the
united-atom picture) had added to its head. With `"InitialEvents" : 135` the
rule is applied to 135 random groups at the start of the reactions, 5 % of
the number of monomers. With `"Probability" : 0.0` it never fires afterwards:
all radicals exist from the start, as in the reference study.

**Propagation, four rules.** A radical can be on the tail or on the head,
and it can add to the head or to the tail of a monomer:

| Rule | Radical on | Adds to | Bond | New radical on | Probability |
|---|---|---|---|---|---|
| `propagation` | `R.tail` | `M.head` | head-to-tail, the regular one | `M.tail` | 0.5 |
| `propagation_tail_tail` | `R.tail` | `M.tail` | a mis-addition | `M.head` | 0.01 |
| `propagation_head_head` | `R.head` | `M.head` | returns to the regular radical | `M.tail` | 0.5 |
| `propagation_head_tail` | `R.head` | `M.tail` | | `M.head` | 0.01 |

The four differ in the reactant state of `R`, in the capture pair and the
bond, and in the product state of `M`. Addition to the tail carbon is
disfavoured because it gives the less stable primary radical; the rules
express that with a probability fifty times lower. Each rule has its own
capture pair. When both carbons of a monomer are within reach of a radical,
the closer pair is tried first with the probability of its rule, and the
other gets its own trial if that fails.

**Termination by combination, three rules.** Two radicals bond and both
groups become `dead`:

| Rule | Reactants | Bond |
|---|---|---|
| `termination` | `radical_tail` + `radical_tail` | `A.tail`–`B.tail` |
| `termination_head_tail` | `radical_head` + `radical_tail` | `A.head`–`B.tail` |
| `termination_head_head` | `radical_head` + `radical_head` | `A.head`–`B.head` |

No atoms are retyped: both groups were retyped when they became radicals.
When the two reactants have the same type and state the pair has no order,
and the labels `A` and `B` are assigned at random.

#### 6. Stereocentres {#polymers-stereo}

When a tail carbon gets its second backbone bond it has three different
substituents besides its implicit hydrogen, and is a stereocentre. A
`"Stereo"` entry declares it:

```json
"Stereo" : {"Centre" : "R.tail", "Substituents" : ["R.head", "M.head", "R.side"], "Next" : "M"}
```

The configuration of the centre is the sign of the normalised volume spanned
by the three bond vectors \f$\mathbf{b}_1, \mathbf{b}_2, \mathbf{b}_3\f$ from
the centre to the substituents, in the order listed,

\f[
\hat V = \frac{(\mathbf{b}_1 \times \mathbf{b}_2)\cdot\mathbf{b}_3}
{b_1\,b_2\,b_3} .
\f]

The order (previous unit, next unit, side group) follows the chain, so that
two consecutive centres with the same sign form a *meso* dyad and with
opposite signs a *racemo* dyad. `"Next"` names the reactant whose own centre
comes next along the chain, so that the dyads can be paired. A rule that
creates two centres (the tail-to-tail additions, and termination) lists two
entries.

The HDDA example sets no `"Tacticity"`, so the configurations are recorded
and the dyad and triad fractions reported, but nothing is selected or
restrained. With `"Tacticity"` the configuration of each new centre is
prescribed and held with a restraint on \f$\hat V\f$; see example 10 of the
[polymer examples](\ref examples_polymers).

#### 7. The conversion

```json
"Conversion" : [{"Type" : "acrylate", "States" : ["radical_tail", "radical_head", "dead"]}]
```

The conversion is the fraction of the groups of the listed types that are in
the listed states: here, the fraction of double bonds that have reacted.
Without the key it is the fraction of all groups that have left their initial
state, which is the same for this example. The reactions stop when the
conversion reaches `"TargetConversion"`.

#### 8. The product tables

A reaction creates bonded terms that no monomer has. Their parameters are
given in the `"Polymerization"` block, in the format of a molecule file and
by pseudo-atom *type*:

```json
"Bonds" : [
  [["CH2_bb", "CH_bb"], "HARMONIC", [301931.64, 1.54]],
  [["CH_bb", "CH_bb"], "HARMONIC", [301931.64, 1.54]],
  [["CH2_bb", "CH2_bb"], "HARMONIC", [301931.64, 1.54]],
  [["CH_bb", "C_CO"], "HARMONIC", [301931.64, 1.52]]
],
"Bends" : [
  [["CH2_bb", "CH_bb", "CH2_bb"], "HARMONIC", [62500.0, 114.0]],
  [["CH2_bb", "CH_bb", "C_CO"], "HARMONIC", [62500.0, 114.0]],
  [["CH_bb", "CH2_bb", "CH_bb"], "HARMONIC", [62500.0, 114.0]]
],
"Torsions" : [
  [["CH2_bb", "CH_bb", "CH2_bb", "CH_bb"], "TRAPPE_EXTENDED", [1078.16, 355.03, 68.19, 791.32, 0.0]],
  [["CH2_bb", "CH_bb", "C_CO", "O_CO"], "TRAPPE_EXTENDED", [823.03, -47.91, -773.13, -1.99, 0.0]]
]
```

(an excerpt; the example lists 6 bonds, 11 bends and 18 torsions). The lookup
is by the types the atoms have *after* the reaction, in either direction.
`"TemplateFile" : "pHDDA-20.json"` reads the tables from a molecule file of
the polymer instead, and inline entries are added to it.

Follow one regular propagation to see which entries are used. The tail of
`R` (type `CH_bb`) bonds to the head of `M` (retyped to `CH2_bb`):

| Created | Atoms | Table entry |
|---|---|---|
| Bond | `R.tail`–`M.head` | `CH_bb`–`CH2_bb` |
| Bends at `R.tail` | `R.head`–`R.tail`–`M.head` | `CH2_bb`–`CH_bb`–`CH2_bb` |
| | `R.side`–`R.tail`–`M.head` | `C_CO`–`CH_bb`–`CH2_bb` |
| Bend at `M.head` | `R.tail`–`M.head`–`M.tail` | `CH_bb`–`CH2_bb`–`CH_bb` |
| Torsions about the new bond | `R.head`–`R.tail`–`M.head`–`M.tail` | `CH2_bb`–`CH_bb`–`CH2_bb`–`CH_bb` |
| | `R.side`–`R.tail`–`M.head`–`M.tail` | `C_CO`–`CH_bb`–`CH2_bb`–`CH_bb` |
| Torsions ending on the new bond | O–`R.side`–`R.tail`–`M.head` (twice) | `O_CO`–`C_CO`–`CH_bb`–`CH2_bb`, `O_es`–`C_CO`–`CH_bb`–`CH2_bb` |
| | previous unit–`R.head`–`R.tail`–`M.head` | `CH_bb`–`CH2_bb`–`CH_bb`–`CH2_bb` |
| | `R.tail`–`M.head`–`M.tail`–`M.side` | `CH_bb`–`CH2_bb`–`CH_bb`–`C_CO` |

The same reaction retypes the head and tail of `M`, so the terms of that
molecule that contain them are overridden:

| Term of the monomer | Becomes |
|---|---|
| Bond `CH2_sp2`–`CH_sp2`, 1.33 Å | `CH2_bb`–`CH_bb`, 1.54 Å: the double bond becomes a single bond |
| Bond `CH_sp2`–`C_CO` | `CH_bb`–`C_CO` |
| Bend `CH2_sp2`–`CH_sp2`–`C_CO`, 119.7° | `CH2_bb`–`CH_bb`–`C_CO`, 114° |
| Bends `CH_sp2`–`C_CO`–`O_CO` and `CH_sp2`–`C_CO`–`O_es` | the same with `CH_bb` |
| Torsions `CH2_sp2`–`CH_sp2`–`C_CO`–O and `CH2_O`–`O_es`–`C_CO`–`CH_sp2` | the same with `CH2_bb` and `CH_bb` |

And the non-bonded pairs across the new bond change: `R.tail`–`M.head` is
now a 1-2 pair, and `R.head`–`M.head`, `R.side`–`M.head` and
`R.tail`–`M.tail` are 1-3 pairs, all excluded; the pairs three bonds apart
are scaled as the 1-4 pairs of the component.

The entries with two equal backbone types (`CH_bb`–`CH_bb`,
`CH2_bb`–`CH2_bb`, and the bends and torsions around them) are not in a
regular head-to-tail chain. They are needed by the terminations and by the
mis-additions. **Every bond a rule can form needs its terms.** A missing
entry is an error at the moment of the reaction, and the message names the
tuple of types to add.

#### 9. The settings of the block

```json
"Polymerization" : {
  "ScanEvery" : 10000,
  "StartCycle" : 0,
  "TargetConversion" : 1.0,
  "RelaxationSteps" : 2000,
  "WriteEventFrames" : true,
  "EventFrameEvery" : 1,
  "Reactions" : [ ... ],
  "Conversion" : [ ... ],
  "Bonds" : [ ... ], "Bends" : [ ... ], "Torsions" : [ ... ]
}
```

| Key | Default | Meaning |
|---|---|---|
| `"ScanEvery"` | 10 | Steps between two scans for reactions. With a time step of 1 fs, 10000 steps is the 10 ps of the reference study. |
| `"StartCycle"` | 0 | The production step at which the reactions start. |
| `"TargetConversion"` | 1.0 | The conversion at which the reactions stop. The dynamics continues. |
| `"RelaxationSteps"` | 1000 | The length of the ramp. |
| `"MaximumReactionsPerScan"` | unlimited | A cap on the reactions of one scan. |
| `"Tacticity"` | none | The configuration of the stereocentres. |
| `"TemplateFile"` | | A molecule file of the polymer to take the product tables from. |
| `"WriteBondLog"` | true | One line per reaction in `output/polymerization_bonds.s<i>.txt`. |
| `"WriteTopology"` | true | The groups, states, bonds and radicals in `output/polymer_topology.s<i>.json`. |
| `"WriteEventFrames"` | false | The configuration at every scan with a reaction, as a LAMMPS data file with the polymer bonds, and a JSON record per frame. |
| `"EventFrameEvery"` | 1 | Write the data file of every \f$n\f$-th frame only. |

All keys, and the format of the output files, are listed on the
[commands](\ref commands) page.

#### 10. Probabilities and rates

A probability is per candidate and per scan, so it is tied to
`"ScanEvery"`. A rate constant \f$k\f$ is not:

\f[
p = 1 - e^{-k\,\Delta t_\text{scan}}, \qquad
\Delta t_\text{scan} = \texttt{ScanEvery} \times \texttt{TimeStep} .
\f]

Example 9 (`9_polymerization_hdda_in_box_realistic`) is the same system with
slower kinetics, and shows the other options:

| | Example 8 | Example 9 |
|---|---|---|
| Temperature | 600 K | 300 K |
| Initiation | 135 radicals at the start, none later | 5 at the start, then `"RateConstant" : 4e-07` ps⁻¹, about \f$4\times10^{-6}\f$ per group per scan |
| Propagation probability | 0.5 (mis-additions 0.01) | 0.05 (mis-additions 0.001) |
| Orientation | distance criterion only | `"OrientationTemperature" : 2000.0` |
| Target conversion | 1.0 | 0.9 |

In example 9 radicals are created throughout the run, chains start at
different times, and a pair must approach in a geometry close to that of the
product to react.

#### 11. What to read in the output

The status report gives, at every print interval, the number of groups per
state, the number of active groups, the conversion, the number of events per
rule, the dyad statistics, and the network statistics of the graph whose
nodes are the molecules and whose edges are the polymer bonds: the number of
clusters, the size of the largest one (the gel fraction), the cycles, and
whether a cluster is connected to its own periodic image. The first
appearance of such a spanning cluster is the gel point of the box.

`"ComputeRadicalWaitingTimes"` collects the time every radical waits between
two reactions; a radical that stops reacting is trapped.
`"ComputePolymerChains"` measures the size and shape of the chains while
they grow, separately for every chain length and only after the force ramp of
their last reaction. Repeated runs at the same temperature are then combined
per length into the structure of the polymer as it forms. See the
[properties](\ref properties) page.

### Writing the rules for another chemistry

1. Identify the reactive groups and give the atoms that take part a role.
2. List the states a group goes through. Use a separate state wherever the
   next reaction differs: the HDDA example needs two radical states because a
   head radical and a tail radical bond through different atoms.
3. Write one rule per transition. A group that reacts more than once (an
   amine with two hydrogens) gets a chain of states and one rule per step.
4. Add pseudo-atoms for the atoms whose bonding changes, and retype them in
   the rules.
5. Give the product tables. Running the simulation reports missing entries
   by name.

A step-growth chemistry is written in the same way. Example 11, an
epoxy-amine cure, uses two group types, the states primary, secondary and
tertiary of the amine, and `"BreakBonds"` for the opening of the epoxide
ring. No atom is ever removed: a proton that moves must be an implicit
hydrogen of a united-atom site, carried by the retyping.

## Limits

- The polymerization is irreversible and its rates are input. It runs with
  [spatial-decomposition molecular dynamics](\ref driver_molecular_dynamics_spatial_decomposition):
  molecules in a box, without a framework.
- Only bonds, bends and torsions are created and overridden. A bond that
  takes part in a Urey-Bradley, inversion, improper, cross or CMAP term
  cannot be broken.
- Cross-links and polymerization cannot be used in the same system.
- The number of molecules is fixed once the reactions have started.

## References

- A. Torres-Knoop, I. Kryven, V. Schamboeck and P.D. Iedema, *Modeling the
  free-radical polymerization of hexanediol diacrylate (HDDA): a molecular
  dynamics and graph theory approach*, Soft Matter **14**, 3404 (2018).
- J.R. Gissinger, B.D. Jensen and K.E. Wise, *Modeling chemical reactions in
  classical molecular dynamics simulations*, Polymer **128**, 211 (2017).
- L.R. Dodd, T.D. Boone and D.N. Theodorou, *A concerted rotation algorithm
  for atomistic Monte Carlo simulation of polymer melts and glasses*,
  Mol. Phys. **78**, 961 (1993).
- N.C. Karayiannis, V.G. Mavrantzas and D.N. Theodorou, *A Novel Monte Carlo
  Scheme for the Rapid Equilibration of Atomistic Model Polymer Systems of
  Precisely Defined Molecular Architecture*, Phys. Rev. Lett. **88**, 105503
  (2002).
- F. Smallenburg and F. Sciortino, *Liquids more stable than crystals in
  particles with limited valence and flexible bonds*, Nature Phys. **9**, 554
  (2013).
- K.A. Maerzke, N.E. Schultz, R.B. Ross and J.I. Siepmann, *TraPPE-UA Force
  Field for Acrylates and Monte Carlo Simulations for Their Mixtures with
  Alkanes and Alcohols*, J. Phys. Chem. B **113**, 6415 (2009).
