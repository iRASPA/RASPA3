# Examples Polymers
\page examples_polymers Examples Polymers


## Table of Contents
1. [Monte Carlo: cross-linked telechelic chains in box (reversible network)](#Example_polymers_1)
2. [Polymerization: free-radical curing of HDDA in box](#Example_polymers_4)
3. [Polymerization: isotactic poly(methyl acrylate) chains in box](#Example_polymers_5)


#### Monte Carlo: cross-linked telechelic chains in box (reversible network) <a name="Example_polymers_1"></a>

40 united-atom butane-like chains in a \f$25 \times 25 \times 25\f$ &Aring; box at
400 K whose two end beads are *reactive sites*: any two ends of different
molecules may be joined by a reversible inter-molecular bond, so that the
chains associate into dimers, longer chains, rings and a network whose
connectivity fluctuates at equilibrium. It demonstrates the cross-link
machinery: the topology moves that form, break and rewire the links, and the
regrowth of linked molecules with their linked sites kept in place.

The molecule declares its reactive atoms with `"ReactiveSites"` (site type
`"X"`, valence 1: each end carries at most one link); everything else in
`telechelic.json` is the TraPPE butane of basic example 4:

```json
  "ReactiveSites" : [
    [0, "X"],
    [3, "X"]
  ]
```

The bond that can form between two `"X"` sites is declared per system with
`"CrossLinkBonds"`: a soft harmonic bond centred at the Lennard-Jones contact
distance, a junction bend on the angles CH2–CH3···CH3' across the link, the
capture radius within which the formation move proposes a link (a proposal
parameter only) and a formation energy that sets the association constant
and thereby the degree of cross-linking. The three topology moves are enabled
per system with `"CrossLinkFormationProbability"` (formation/scission),
`"CrossLinkSwapProbability"` (one end of a link moves to a free site, keeping
the number of links constant) and `"CrossLinkExchangeProbability"` (two links
trade partners, (a,b)+(c,d) → (a,c)+(b,d), keeping every site's number of
links constant; the move that keeps a saturated network mobile). Linked
molecules are still translated, rotated and regrown: the
reinsertion and partial-reinsertion moves keep the linked end beads fixed and
regrow the rest of the chain with the link's bond, junction bends and
exclusion corrections entering the Rosenbluth weights (a fixed-endpoint
regrowth when both ends are linked).

```json
{
  "SimulationType" : "MonteCarlo",
  "NumberOfProductionCycles" : 10000,
  "NumberOfInitializationCycles" : 2000,
  "PrintEvery" : 1000,

  "Systems" :
  [
    {
      "Type" : "Box",
      "BoxLengths" : [25.0, 25.0, 25.0],
      "ExternalTemperature" : 400.0,
      "ChargeMethod" : "None",
      "CrossLinkSwapProbability" : 0.5,
      "CrossLinkFormationProbability" : 1.0,
      "CrossLinkExchangeProbability" : 0.5,
      "CrossLinkBonds" : [
        {
          "Sites" : ["X", "X"],
          "Bond" : ["HARMONIC", [1000.0, 3.8]],
          "JunctionBend" : ["HARMONIC", [500.0, 120.0]],
          "CaptureRadius" : 6.5,
          "FormationEnergy" : -1000.0
        }
      ]
    }
  ],

  "Components" :
  [
    {
      "Name" : "telechelic",
      "TranslationProbability" : 1.0,
      "RotationProbability" : 1.0,
      "ReinsertionProbability" : 1.0,
      "PartialReinsertionProbability" : 1.0,
      "CreateNumberOfMolecules" : 40
    }
  ]
}
```

The status reports list the current `Number of cross-links` next to the
number of molecules (here about 28 of the 40 possible links, fluctuating
between roughly 24 and 31), the energy breakdown carries the bonded link
terms in the `cross-link` slot, and the move statistics report the
acceptance of the formation/scission, swap and exchange moves alongside the
ordinary moves. Making `"FormationEnergy"` more negative drives the system
towards a fully linked network (at `-2500` K essentially all 40 links are
present); there the swap move finds hardly any free site to relocate to,
while the exchange move needs none and keeps rewiring the network (its
acceptance is modest, a few percent, because the second new bond joins two
sites whose separation is not controlled by the capture radius, but the move
costs only four bonded-energy evaluations). Setting it to zero leaves
the association to the bond well alone. The links are stored in the restart file, so a run can be
continued from a linked state.


#### Polymerization: free-radical curing of HDDA in box <a name="Example_polymers_4"></a>

150 hexanediol diacrylate (HDDA) monomers in a \f$38 \times 38 \times 38\f$
&Aring; box at 300 K (about the liquid density) are cured into a cross-linked
poly-HDDA network with `"SimulationType" : "Polymerization"`: a
spatial-decomposition MD run in which, during the production stage, the vinyl
groups react by stochastic free-radical chain growth. The monomer file
`HDDA.json` is the 16-bead TraPPE-UA acrylate of the pHDDA examples (the
polymer file `pHDDA-10.json` of
`examples/polymers/2_mc_end_to_end_distance_phdda_20_in_box` is a *pre-built*
chain of ten such units); its two double bonds are declared as reactive sites
with the pseudo-atoms they become in the polymer backbone:

```json
  "PolymerizationSites" : [
    {"Type" : "acrylate", "Head" : 0,  "Tail" : 1,  "HeadProductType" : "CH2_bb", "TailProductType" : "CH_bb"},
    {"Type" : "acrylate", "Head" : 13, "Tail" : 12, "HeadProductType" : "CH2_bb", "TailProductType" : "CH_bb"}
  ]
```

The `"Polymerization"` block of the system sets the kinetics and gives the
bonded tables of the product by pseudo-atom type (the `CH2_bb`/`CH_bb` entries
of the pHDDA-10 molecule file, extended with the `CH_bb–CH_bb` terms that a
termination by combination needs). Eight random sites become radicals at the
start of production; every 10 MD steps (10 fs) each radical tail within 4
&Aring; of a monomer head adds to it with probability 0.05, and two radical
tails within 4 &Aring; combine with the same probability. An accepted reaction
retypes the head and tail to `CH2_bb`/`CH_bb` (which changes their
Lennard-Jones parameters and charges), adds the new bond with its bends,
torsions and 1-2/1-3 exclusions and scaled 1-4 pairs, and ramps the new
terms in over 2 ps so that the capture does not heat the system (the
thermostat holds 300 K within a few kelvin). The reactions stop at 90%
conversion and the dynamics continues on the cured network.

```json
{
  "SimulationType" : "Polymerization",
  "NumberOfThreads" : 4,
  "NumberOfInitializationCycles" : 200,
  "NumberOfEquilibrationCycles" : 20000,
  "NumberOfProductionCycles" : 2000000,
  "PrintEvery" : 10000,

  "Systems" : [
    {
      "Type" : "Box",
      "BoxLengths" : [38.0, 38.0, 38.0],
      "ExternalTemperature" : 300.0,
      "ExternalPressure" : 100000.0,
      "Ensemble" : "NPT",
      "TimeStep" : 0.001,
      "ChargeMethod" : "Ewald",
      "Polymerization" : {
        "ScanEvery" : 10,
        "TargetConversion" : 0.9,
        "RelaxationSteps" : 2000,
        "Initiation" :  {"InitialRadicals" : 8, "Probability" : 0.0},
        "Propagation" : {"CaptureRadius" : 4.0, "Probability" : 0.05, "AllowIntramolecular" : true},
        "Termination" : {"Enabled" : true, "CaptureRadius" : 4.0, "Probability" : 0.05},
        "Bonds" : [
          [["CH2_bb", "CH_bb"], "HARMONIC", [301931.64, 1.54]],
          [["CH_bb", "CH_bb"], "HARMONIC", [301931.64, 1.54]],
          [["CH_bb", "C_CO"], "HARMONIC", [301931.64, 1.52]]
        ],
        "Bends" : [ "..." ],
        "Torsions" : [ "..." ]
      }
    }
  ],

  "Components" : [
    {
      "Name" : "HDDA",
      "TranslationProbability" : 1.0,
      "RotationProbability" : 1.0,
      "ReinsertionProbability" : 0.1,
      "PartialReinsertionProbability" : 5.0,
      "CrankshaftProbability" : 1.0,
      "BeadDisplacementProbability" : 1.0,
      "CreateNumberOfMolecules" : 150
    }
  ]
}
```

The status report of every `"PrintEvery"` cycles adds a `Polymerization
status` block with the conversion, the numbers of monomer, radical and dead
sites, the numbers of initiations, propagations and terminations, and the
network statistics: the molecules connected into polymer, the number of
clusters, the size of the largest one (its fraction of all molecules is the
gel fraction), the number of bonds that closed a cycle within one molecule
(primary cycles, the two acrylates of one HDDA reacting with the same chain),
the cycle rank of the molecule graph and the number of independent loops per
box of the unfolded network (the cycle rank minus one per spanning cluster;
the effective cross-links), and the percolation analysis: whether a cluster
is bonded to its own periodic image through the boundary (a *winding loop*
of the molecule graph), in how many directions, and the fraction of the
molecules in such spanning clusters. The appearance of the first spanning
cluster marks the gel point of the finite box. The bonded energy of the polymer bonds is
reported in the `cross-link` slot of the energy breakdown. Every reaction is
appended to `output/polymerization_bonds.s0.txt` (step, cycle, event, the two
sites as component/molecule/atom, the distance at which the bond formed and
the conversion), and `output/polymer_topology.s0.json` holds the final site
states and the bond list, from which the polymer graph can be rebuilt. The
polymer bonds, the retyped atoms and the reaction state are stored in the
binary restart file, so a curing run can be continued. Within 12 ps of a
test run the conversion reaches about 30% with one gel cluster of 68 of the
150 molecules; a complete cure takes of the order of a nanosecond, which is
why the production stage is long and the reaction probabilities are small
(in the real system the propagation rate constant of acrylates is of order
\f$10^4\f$ L mol⁻¹ s⁻¹, far slower than any MD run, so the rate here is a
compromise that keeps the growth diffusion-limited and lets the network relax
between reactions).


#### Polymerization: isotactic poly(methyl acrylate) chains in box <a name="Example_polymers_5"></a>

200 methyl acrylate monomers (the six-bead TraPPE-UA ester
`CH2_sp2=CH_sp2–C_CO(=O_CO)–O_es–CH3_O`) in a \f$32 \times 32 \times 32\f$
&Aring; box at 300 K and 1 bar (the box relaxes to the liquid density, about
31 &Aring;) are polymerized into linear chains with a
controlled tacticity. The monomer has one vinyl group; its site declares the
ester carbon as the side group that defines the configuration of the
backbone CH (it is the only substituent of the tail besides the head, so the
key could be left out):

```json
  "PolymerizationSites" : [
    {"Type" : "acrylate", "Head" : 0, "Tail" : 1, "SideGroup" : 2, "HeadProductType" : "CH2_bb", "TailProductType" : "CH_bb"}
  ]
```

Four random sites are initiated at the start of production and termination
is switched off, so four chains grow until 90% of the monomers have reacted
(about 45 units each). The `"Tacticity"` block asks for isotactic chains:

```json
      "Polymerization" : {
        "ScanEvery" : 10,
        "TargetConversion" : 0.9,
        "RelaxationSteps" : 2000,
        "Initiation" :  {"InitialRadicals" : 4, "Probability" : 0.0},
        "Propagation" : {"CaptureRadius" : 4.0, "Probability" : 0.1},
        "Termination" : {"Enabled" : false},
        "Tacticity" : {
          "Type" : "Isotactic",
          "StereoSelective" : true,
          "ForceConstant" : 40000.0,
          "ReferenceVolume" : 0.6
        },
        "Bonds" : [ "..." ], "Bends" : [ "..." ], "Torsions" : [ "..." ]
      }
```

The backbone CH of a reacted unit is a stereocentre (head, side group,
polymer partner, and the implicit hydrogen of the united-atom bead) whose
configuration is the sign of the chiral volume of its three explicit
substituents, ordered along the chain. The first unit of a chain is free;
when its radical adds to a monomer the face from which the monomer
approaches decides the configuration, and from then on every new radical
carries the configuration it has to take: the same sign as the centre
before it (a meso dyad) for isotactic chains, the opposite sign (racemo) for
syndiotactic ones, and a Bernoullian choice with `"MesoProbability"` for
atactic ones. With `"StereoSelective"` a radical–monomer pair whose geometry
would give the wrong configuration is skipped in that scan (the monomer has
to come around to the other face; the number of such skipped candidates is
reported at the end), so the selection costs reaction rate but no energy.
Once reacted, every centre is held in its configuration by a harmonic
restraint on the normalised chiral volume, ramped in with the bond, because
the united-atom CH has no hydrogen that would otherwise keep it from
inverting through the planar arrangement (the three bends alone give a
barrier of only a few \f$k_B T\f$). The `Polymerization status` block
reports the stereocentres, the meso/racemo dyad counts and the meso
fraction, the triad fractions mm/mr/rr, and how many restrained centres are
momentarily inverted; for this input the meso fraction is 1 (mm = 1) by
construction, and the same input with `"Type" : "Syndiotactic"` gives meso
fraction 0 (rr = 1), with `"Type" : "Atactic"` about 0.5 (mm : mr : rr
≈ 1 : 2 : 1). Every reaction in `output/polymerization_bonds.s0.txt`
carries the configuration of the reacting tail and its dyad (`+m`, `-m`,
`+.` for the first centre of a chain), and `output/polymer_topology.s0.json`
holds the configuration of every site next to its state, so the
stereosequence of every chain can be read off. Leaving out the
`"Tacticity"` block polymerizes without selection or restraints; the dyad
statistics are then still reported and show the natural (near 1 : 1)
outcome of the capture geometry.
