# Examples Polymers
\page examples_polymers Examples Polymers


## Table of Contents
- [1. Monte Carlo: cross-linked telechelic chains in box (reversible network)](#Example_polymers_1)
- [8. Polymerization: free-radical curing of HDDA in box](#Example_polymers_8)
- [9. Polymerization: free-radical curing of HDDA in box at realistic rates](#Example_polymers_9)
- [10. Polymerization: isotactic poly(methyl acrylate) chains in box](#Example_polymers_10)
- [11. Polymerization: epoxy–amine curing with ring opening in box](#Example_polymers_11)


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


#### Polymerization: free-radical curing of HDDA in box <a name="Example_polymers_8"></a>

Run from `examples/polymers/8_polymerization_hdda_in_box`.
2690 hexanediol diacrylate (HDDA) monomers in a \f$111.4 \times 111.4 \times 111.4\f$
&Aring; box at 600 K and 1 bar (NPT) are cured into a cross-linked
poly-HDDA network with `"SimulationType" : "Polymerization"`, following the
protocol of A. Torres-Knoop, I. Kryven, V. Schamboeck and P. D. Iedema,
*Soft Matter* **14**, 3404 (2018) (TraPPE-UA acrylates with harmonic bonds,
curing at 600 K in NPT, Nos&eacute;&ndash;Hoover coupling 1 ps, 5 % of the
monomers initially radical, capture radius 4 &Aring;, reaction probability 0.5,
proximity check every 10 ps): a
spatial-decomposition MD run (`"PairDevice" : "GPU"`) in which, during the production stage, the vinyl
groups react by stochastic free-radical chain growth. The monomer file
`HDDA.json` is the 16-bead TraPPE-UA acrylate of the pHDDA examples (the
polymer file `pHDDA-20.json` of
`examples/polymers/2_mc_end_to_end_distance_phdda_20_in_box` is a *pre-built*
chain of twenty such units); its two double bonds are declared as reactive
groups of type `acrylate`, each with its `head` (the `CH2=` a radical adds
to), `tail` (the `=CH–` that carries the radical afterwards) and `side`
(the ester carbon on the tail) atoms named by role:

```json
  "ReactiveGroups" : [
    {"Type" : "acrylate", "State" : "monomer", "Atoms" : {"head" : 0,  "tail" : 1,  "side" : 2}},
    {"Type" : "acrylate", "State" : "monomer", "Atoms" : {"head" : 13, "tail" : 12, "side" : 11}}
  ]
```

The `"Polymerization"` block of the system spells out the chemistry as
reaction rules between the states of the groups and gives the bonded tables
of the product by pseudo-atom type (the `CH2_bb`/`CH_bb` entries of the
pHDDA-20 molecule file, extended with the `CH_bb–CH_bb` and `CH2_bb–CH2_bb`
terms that a termination by combination and the head-to-head linkages need).
The radical of a group is located by its state: `radical_tail` is the usual
secondary radical on the `tail` carbon (`–CH•–COOR`), `radical_head` the
primary radical on the `head` carbon (`–CH2•`) that a mis-addition leaves
behind. The *initiation* rule turns a `monomer`
group into a `radical_tail` and retypes its head and tail to the backbone
pseudo-atoms `CH2_bb`/`CH_bb`; it is applied to 135 random groups at the
start of production (`"InitialEvents"`) and never afterwards
(`"Probability" : 0.0`). The *propagation* rule pairs a `radical_tail` group
`R` with a `monomer` group `M` whose `R.tail`–`M.head` distance is within 4
&Aring;, forms that bond, leaves `R` as `dead` and makes `M` the new
`radical_tail` (retyped), with probability 0.5 per scan every 10000 MD steps
(10 ps): the regular head-to-tail addition. Three further propagation rules
give the regiochemistry: `propagation_tail_tail` adds `R.tail` to `M.tail`
instead and leaves the radical on `M.head` (`radical_head`), with probability
0.01, i.e. 2 % of the regular rate when both monomer carbons are in range;
`propagation_head_head` lets a `radical_head` add to `M.head` and returns
the radical to a tail (probability 0.5, a tail-to-tail linkage in the
chain); `propagation_head_tail` lets it add to `M.tail` and keeps a head
radical (0.01). Each rule has its own capture pair, so for a radical with
both carbons of a monomer in range the closer pair is tried first with its
rule's probability and the other gets its own trial when that fails. The
*termination* rules bond two radicals within 4 &Aring; with probability 0.5
and leave both `dead`, one rule per combination of radical states
(`tail–tail`, `head–tail`, `head–head`). The conversion is
the fraction of acrylate groups in the states `radical_tail`, `radical_head`
or `dead`. An
accepted reaction retypes the atoms (which changes their Lennard-Jones
parameters and charges), adds the new bond with its bends, torsions and
1-2/1-3 exclusions and scaled 1-4 pairs, and ramps the new terms in over 2 ps
so that the capture does not heat the system. The reactions run to full
conversion and the dynamics continues on the cured network. Example 9 is the
same monomers with the slow, orientation-weighted kinetics.

```json
{
  "SimulationType" : "Polymerization",
  "NumberOfThreads" : 8,
  "PairDevice" : "GPU",
  "NumberOfInitializationCycles" : 50,
  "NumberOfEquilibrationCycles" : 500000,
  "NumberOfProductionCycles" : 10000000,
  "PrintEvery" : 10000,

  "Systems" : [
    {
      "Type" : "Box",
      "BoxLengths" : [111.4, 111.4, 111.4],
      "ExternalTemperature" : 600.0,
      "ExternalPressure" : 100000.0,
      "Ensemble" : "NPT",
      "TimeScaleParameterThermostat" : 1.0,
      "TimeStep" : 0.001,
      "ChargeMethod" : "Ewald",
      "Polymerization" : {
        "ScanEvery" : 10000,
        "TargetConversion" : 1.0,
        "RelaxationSteps" : 2000,
        "WriteEventFrames" : true,
        "EventFrameEvery" : 1,
        "Reactions" : [
          {
            "Name" : "initiation",
            "Reactants" : [{"Label" : "M", "Type" : "acrylate", "State" : "monomer"}],
            "Products" : [{"Label" : "M", "State" : "radical_tail", "Retype" : {"head" : "CH2_bb", "tail" : "CH_bb"}}],
            "Probability" : 0.0,
            "InitialEvents" : 135
          },
          {
            "Name" : "propagation",
            "Reactants" : [{"Label" : "R", "Type" : "acrylate", "State" : "radical_tail"},
                           {"Label" : "M", "Type" : "acrylate", "State" : "monomer"}],
            "Capture" : {"Atoms" : ["R.tail", "M.head"], "Radius" : 4.0},
            "FormBonds" : [["R.tail", "M.head"]],
            "Products" : [{"Label" : "R", "State" : "dead"},
                          {"Label" : "M", "State" : "radical_tail", "Retype" : {"head" : "CH2_bb", "tail" : "CH_bb"}}],
            "Stereo" : {"Centre" : "R.tail", "Substituents" : ["R.head", "M.head", "R.side"], "Next" : "M"},
            "Probability" : 0.5,
            "AllowIntramolecular" : true
          },
          {
            "Name" : "propagation_tail_tail",
            "Reactants" : [{"Label" : "R", "Type" : "acrylate", "State" : "radical_tail"},
                           {"Label" : "M", "Type" : "acrylate", "State" : "monomer"}],
            "Capture" : {"Atoms" : ["R.tail", "M.tail"], "Radius" : 4.0},
            "FormBonds" : [["R.tail", "M.tail"]],
            "Products" : [{"Label" : "R", "State" : "dead"},
                          {"Label" : "M", "State" : "radical_head", "Retype" : {"head" : "CH2_bb", "tail" : "CH_bb"}}],
            "Stereo" : [{"Centre" : "R.tail", "Substituents" : ["R.head", "M.tail", "R.side"]},
                        {"Centre" : "M.tail", "Substituents" : ["M.head", "R.tail", "M.side"]}],
            "Probability" : 0.01,
            "AllowIntramolecular" : true
          },
          {
            "Name" : "propagation_head_head",
            "Reactants" : [{"Label" : "R", "Type" : "acrylate", "State" : "radical_head"},
                           {"Label" : "M", "Type" : "acrylate", "State" : "monomer"}],
            "Capture" : {"Atoms" : ["R.head", "M.head"], "Radius" : 4.0},
            "FormBonds" : [["R.head", "M.head"]],
            "Products" : [{"Label" : "R", "State" : "dead"},
                          {"Label" : "M", "State" : "radical_tail", "Retype" : {"head" : "CH2_bb", "tail" : "CH_bb"}}],
            "Probability" : 0.5,
            "AllowIntramolecular" : true
          },
          {
            "Name" : "propagation_head_tail",
            "Reactants" : [{"Label" : "R", "Type" : "acrylate", "State" : "radical_head"},
                           {"Label" : "M", "Type" : "acrylate", "State" : "monomer"}],
            "Capture" : {"Atoms" : ["R.head", "M.tail"], "Radius" : 4.0},
            "FormBonds" : [["R.head", "M.tail"]],
            "Products" : [{"Label" : "R", "State" : "dead"},
                          {"Label" : "M", "State" : "radical_head", "Retype" : {"head" : "CH2_bb", "tail" : "CH_bb"}}],
            "Stereo" : {"Centre" : "M.tail", "Substituents" : ["M.head", "R.head", "M.side"]},
            "Probability" : 0.01,
            "AllowIntramolecular" : true
          },
          {
            "Name" : "termination",
            "Reactants" : [{"Label" : "A", "Type" : "acrylate", "State" : "radical_tail"},
                           {"Label" : "B", "Type" : "acrylate", "State" : "radical_tail"}],
            "Capture" : {"Atoms" : ["A.tail", "B.tail"], "Radius" : 4.0},
            "FormBonds" : [["A.tail", "B.tail"]],
            "Products" : [{"Label" : "A", "State" : "dead"}, {"Label" : "B", "State" : "dead"}],
            "Stereo" : [{"Centre" : "A.tail", "Substituents" : ["A.head", "B.tail", "A.side"]},
                        {"Centre" : "B.tail", "Substituents" : ["B.head", "A.tail", "B.side"]}],
            "Probability" : 0.5
          },
          {
            "Name" : "termination_head_tail",
            "Reactants" : [{"Label" : "A", "Type" : "acrylate", "State" : "radical_head"},
                           {"Label" : "B", "Type" : "acrylate", "State" : "radical_tail"}],
            "Capture" : {"Atoms" : ["A.head", "B.tail"], "Radius" : 4.0},
            "FormBonds" : [["A.head", "B.tail"]],
            "Products" : [{"Label" : "A", "State" : "dead"}, {"Label" : "B", "State" : "dead"}],
            "Stereo" : {"Centre" : "B.tail", "Substituents" : ["B.head", "A.head", "B.side"]},
            "Probability" : 0.5
          },
          {
            "Name" : "termination_head_head",
            "Reactants" : [{"Label" : "A", "Type" : "acrylate", "State" : "radical_head"},
                           {"Label" : "B", "Type" : "acrylate", "State" : "radical_head"}],
            "Capture" : {"Atoms" : ["A.head", "B.head"], "Radius" : 4.0},
            "FormBonds" : [["A.head", "B.head"]],
            "Products" : [{"Label" : "A", "State" : "dead"}, {"Label" : "B", "State" : "dead"}],
            "Probability" : 0.5
          }
        ],
        "Conversion" : [{"Type" : "acrylate", "States" : ["radical_tail", "radical_head", "dead"]}],
        "Bonds" : [
          [["CH2_bb", "CH_bb"], "HARMONIC", [301931.64, 1.54]],
          [["CH_bb", "CH_bb"], "HARMONIC", [301931.64, 1.54]],
          [["CH2_bb", "CH2_bb"], "HARMONIC", [301931.64, 1.54]],
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
      "CreateNumberOfMolecules" : 2690
    }
  ]
}
```

The status report of every `"PrintEvery"` cycles adds a `Polymerization
status` block with the numbers of groups per state (monomer, radical, dead),
the active groups (the living radicals), the conversion, the number of
events per reaction rule, and the network statistics: the molecules connected into polymer, the number of
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
appended to `output/polymerization_bonds.s0.txt` (step, cycle, the rule, the
two groups as component/molecule/group, the distance at which the bond formed,
the conversion and the radical id of each reacting group), and
`output/polymer_topology.s0.json` holds the reaction scheme, the final group
states and the bond list, from which the polymer graph can be rebuilt. Every
living radical carries an id from its initiation (`41*` in the log) through
every propagation (`17/17`: handed from the radical group to the added
monomer) to its termination (`17x/23x`), so the reactions of one chain end
can be followed through the log; the `radicals` table of the topology file
summarises each (step and rule of birth and death, the group it sits or died
on with that group's state, which says whether the radical is on the tail or
the head carbon, the number of propagations, which is the kinetic chain
length), and `groupLastRadical` names for every group the radical that
passed through it, i.e. the chain a dead group belongs to. The `events`
counts of the status block (and the rule names in the log) give the
regiochemistry of the cure: the number of `propagation_tail_tail` and
`propagation_head_tail` events over all propagations is the fraction of
head-to-head linkages. With
`"WriteEventFrames"` (here a data file for every reaction, `"EventFrameEvery": 1`) the
configuration at the moment of a reaction is written as a LAMMPS data file to
`polymerization_events/s0/step_<step>.data`, with the polymer bonds and their
bends and torsions in its topology and the reacted atoms retyped, so the
growing network can be followed in OVITO or VMD (both read the bonds of a
data file); the reactions of every frame, with the atoms before and after
and the new bonded terms, are the records of
`output/polymerization_events.s0.jsonl`. With
`"ComputeRadicalWaitingTimes"` the waiting time of every living radical
between its successive reactions, and the ages of the radicals still living,
are collected in log-spaced histograms in
`radical_waiting_times/radical_waiting_times.s0.txt`; a stretched tail or a
pile-up of old living radicals in that distribution is the signature of
radical trapping. The
polymer bonds, the retyped atoms and the reaction state are stored in the
binary restart file, so a curing run can be continued. A complete cure takes
of the order of a nanosecond even at these elevated probabilities, which is
why the production stage is long (in the real system the propagation rate
constant of acrylates is of order \f$10^4\f$ L mol⁻¹ s⁻¹, far slower than any
MD run). Example 9 keeps the scan frequent and the probability small, which
is the setup that leaves time for the network to relax between reactions.


#### Polymerization: free-radical curing of HDDA in box at realistic rates <a name="Example_polymers_9"></a>

The same 2690 HDDA monomers, monomer file and product tables as example 8,
run from `examples/polymers/9_polymerization_hdda_in_box_realistic`. The box
is \f$125 \times 125 \times 125\f$ &Aring; at 300 K and 1 bar (`"Ensemble" :
"NPT"`). Five groups start as radicals and further radicals are created by
the initiation rule with `"RateConstant" : 4e-7`. Propagation and
termination are attempted every 10 MD steps (10 fs) inside a 4 &Aring;
capture radius with probability 0.05 for the regular head-to-tail addition
and the terminations and 0.001 for the two mis-additions that leave a head
radical (`propagation_tail_tail`, `propagation_head_tail`), and an
orientation weight (`"OrientationTemperature" : 2000` on every bimolecular
rule) rejects poorly aligned approaches. The reactions stop at 90%
conversion. Radical waiting times are collected, as in example 8. The rules
are those of example 8 with the kinetic keys changed:

```json
      "Polymerization" : {
        "ScanEvery" : 10,
        "TargetConversion" : 0.9,
        "RelaxationSteps" : 2000,
        "WriteEventFrames" : true,
        "EventFrameEvery" : 1,
        "Reactions" : [
          {"Name" : "initiation", "...", "RateConstant" : 4e-7, "InitialEvents" : 5},
          {"Name" : "propagation", "...", "Probability" : 0.05, "AllowIntramolecular" : true, "OrientationTemperature" : 2000.0},
          {"Name" : "propagation_tail_tail", "...", "Probability" : 0.001, "AllowIntramolecular" : true, "OrientationTemperature" : 2000.0},
          {"Name" : "propagation_head_head", "...", "Probability" : 0.05, "AllowIntramolecular" : true, "OrientationTemperature" : 2000.0},
          {"Name" : "propagation_head_tail", "...", "Probability" : 0.001, "AllowIntramolecular" : true, "OrientationTemperature" : 2000.0},
          {"Name" : "termination", "...", "Probability" : 0.05, "OrientationTemperature" : 2000.0},
          {"Name" : "termination_head_tail", "...", "Probability" : 0.05, "OrientationTemperature" : 2000.0},
          {"Name" : "termination_head_head", "...", "Probability" : 0.05, "OrientationTemperature" : 2000.0}
        ],
        "Conversion" : [{"Type" : "acrylate", "States" : ["radical_tail", "radical_head", "dead"]}]
      }
```


#### Polymerization: isotactic poly(methyl acrylate) chains in box <a name="Example_polymers_10"></a>

Run from `examples/polymers/10_polymerization_isotactic_methyl_acrylate_in_box`.
200 methyl acrylate monomers (the six-bead TraPPE-UA ester
`CH2_sp2=CH_sp2–C_CO(=O_CO)–O_es–CH3_O`) in a \f$32 \times 32 \times 32\f$
&Aring; box at 300 K and 1 bar (the box relaxes to the liquid density, about
31 &Aring;) are polymerized into linear chains with a
controlled tacticity. The monomer has one vinyl group; its `side` role names
the ester carbon as the substituent that defines the configuration of the
backbone CH:

```json
  "ReactiveGroups" : [
    {"Type" : "acrylate", "State" : "monomer", "Atoms" : {"head" : 0, "tail" : 1, "side" : 2}}
  ]
```

Four random groups are initiated at the start of production and there is no
termination rule, so four chains grow until 90% of the monomers have reacted
(about 45 units each). The `"Stereo"` entry of the propagation rule declares
the stereocentre the reaction creates — the radical's tail with its
substituents ordered along the chain (`R.head`, the new partner `M.head`, and
`R.side`) — and names the monomer as the carrier of the next centre
(`"Next" : "M"`); the `"Tacticity"` block asks for isotactic chains:

```json
      "Polymerization" : {
        "ScanEvery" : 10,
        "TargetConversion" : 0.9,
        "RelaxationSteps" : 2000,
        "WriteEventFrames" : true,
        "EventFrameEvery" : 5,
        "Reactions" : [
          {
            "Name" : "initiation",
            "Reactants" : [{"Label" : "M", "Type" : "acrylate", "State" : "monomer"}],
            "Products" : [{"Label" : "M", "State" : "radical", "Retype" : {"head" : "CH2_bb", "tail" : "CH_bb"}}],
            "Probability" : 0.0,
            "InitialEvents" : 4
          },
          {
            "Name" : "propagation",
            "Reactants" : [{"Label" : "R", "Type" : "acrylate", "State" : "radical"},
                           {"Label" : "M", "Type" : "acrylate", "State" : "monomer"}],
            "Capture" : {"Atoms" : ["R.tail", "M.head"], "Radius" : 4.0},
            "FormBonds" : [["R.tail", "M.head"]],
            "Products" : [{"Label" : "R", "State" : "dead"},
                          {"Label" : "M", "State" : "radical", "Retype" : {"head" : "CH2_bb", "tail" : "CH_bb"}}],
            "Stereo" : {"Centre" : "R.tail", "Substituents" : ["R.head", "M.head", "R.side"], "Next" : "M"},
            "Probability" : 0.1
          }
        ],
        "Conversion" : [{"Type" : "acrylate", "States" : ["radical", "dead"]}],
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
holds every stereocentre with its configuration next to the group states, so
the stereosequence of every chain can be read off. With `"WriteEventFrames"`
every propagation is also recorded in `output/polymerization_events.s0.jsonl`
(the atoms before and after, the configuration taken, the dyad, the new
bonded terms) and one frame in five is written as a LAMMPS data file to
`polymerization_events/s0/`, with the chains as bonded topology; the
stereocentre restraints have no LAMMPS equivalent and are not part of those
files (the input script notes it). Leaving out the
`"Tacticity"` block polymerizes without selection or restraints; the dyad
statistics are then still reported and show the natural (near 1 : 1)
outcome of the capture geometry.

#### Polymerization: epoxy–amine curing with ring opening in box <a name="Example_polymers_11"></a>

Run from `examples/polymers/11_polymerization_epoxy_amine_in_box`.
120 diglycidyl ether molecules (a nine-bead united-atom model with an epoxide
ring `CH2_ep–CH_ep–O_ep` at each end of a `CH2_O–O_eth–CH2_O` ether linker)
and 60 ethylenediamine molecules (`N_p–CH2_N–CH2_N–N_p`, the amine groups as
single beads with their hydrogens folded in) in a
\f$32 \times 32 \times 32\f$ &Aring; box at 350 K and 1 bar are cured into a
network by the step-growth addition of the amines to the epoxides. Unlike
the chain-growth examples, the reaction here *breaks* a bond of the monomer:
the amine nitrogen attacks the terminal CH2 of the ring, the CH2–O ring bond
opens and the oxygen becomes a hydroxyl. Each diglycidyl ether carries two
`epoxide` groups whose roles name the three ring atoms, each diamine two
`amine` groups:

```json
  "ReactiveGroups" : [
    {"Type" : "epoxide", "State" : "closed", "Atoms" : {"CH2" : 0, "CH" : 1, "O" : 2}},
    {"Type" : "epoxide", "State" : "closed", "Atoms" : {"CH2" : 7, "CH" : 6, "O" : 8}}
  ]
```

```json
  "ReactiveGroups" : [
    {"Type" : "amine", "State" : "primary", "Atoms" : {"N" : 0}},
    {"Type" : "amine", "State" : "primary", "Atoms" : {"N" : 3}}
  ]
```

A primary amine adds to one ring and becomes secondary, a secondary amine
adds to a second ring and becomes tertiary (the two reactive hydrogens of an
NH2); both rules form the CH2–N bond, break the CH2–O ring bond
(`"BreakBonds"`) and retype the opened ring (`CH2_N`, the hydroxyl-bearing
`CH_OH` and the hydroxyl `O_OH`) and the nitrogen (`NH`, `N_t`):

```json
      "Polymerization" : {
        "ScanEvery" : 2000,
        "TargetConversion" : 0.9,
        "RelaxationSteps" : 2000,
        "WriteEventFrames" : true,
        "EventFrameEvery" : 5,
        "Reactions" : [
          {
            "Name" : "primary_addition",
            "Reactants" : [{"Label" : "E", "Type" : "epoxide", "State" : "closed"},
                           {"Label" : "A", "Type" : "amine", "State" : "primary"}],
            "Capture" : {"Atoms" : ["E.CH2", "A.N"], "Radius" : 3.5},
            "FormBonds" : [["E.CH2", "A.N"]],
            "BreakBonds" : [["E.CH2", "E.O"]],
            "Products" : [{"Label" : "E", "State" : "opened", "Retype" : {"CH2" : "CH2_N", "CH" : "CH_OH", "O" : "O_OH"}},
                          {"Label" : "A", "State" : "secondary", "Retype" : {"N" : "NH"}}],
            "Probability" : 0.5
          },
          {
            "Name" : "secondary_addition",
            "Reactants" : [{"Label" : "E", "Type" : "epoxide", "State" : "closed"},
                           {"Label" : "A", "Type" : "amine", "State" : "secondary"}],
            "Capture" : {"Atoms" : ["E.CH2", "A.N"], "Radius" : 3.5},
            "FormBonds" : [["E.CH2", "A.N"]],
            "BreakBonds" : [["E.CH2", "E.O"]],
            "Products" : [{"Label" : "E", "State" : "opened", "Retype" : {"CH2" : "CH2_N", "CH" : "CH_OH", "O" : "O_OH"}},
                          {"Label" : "A", "State" : "tertiary", "Retype" : {"N" : "N_t"}}],
            "Probability" : 0.25
          }
        ],
        "Conversion" : [{"Type" : "epoxide", "States" : ["opened"]}],
        "Bonds" : [ "..." ], "Bends" : [ "..." ], "Torsions" : [ "..." ]
      }
```

A broken bond is ramped out as a formed bond is ramped in: over
`"RelaxationSteps"` the ring bond and the ring bends and torsions of the
monomer that run through it go to zero (no table entry is needed for them),
while the new CH2–N bond with its bends and torsions and the retyped terms
of the opened ring (the `CH2_N–CH_OH–O_OH` bend takes the place of the 60°
ring bend) come in from the `"Bonds"`, `"Bends"` and `"Torsions"` tables.
The polymer topology follows the opened ring: the oxygen is three bonds from
the nitrogen (a scaled 1-4 pair), no longer two. Every group reacts at most
as its state allows, so the scan runs once per ramp (`"ScanEvery"` equal to
`"RelaxationSteps"`) and every bond is relaxed before the next scan; with
two reactive hydrogens per amine group and an epoxide on both ends of the
ether, the diamine is tetrafunctional and the product a network. The second
addition to an amine does not disturb the first: the bonded terms and pair
corrections of the earlier bond keep their state, and only what the new
reaction changes (the `N_t` parameters of the nitrogen's terms, the new
bond) is ramped. The `Polymerization status` block counts the `epoxide`
groups per state (`closed`, `opened`), the `amine` groups (`primary`,
`secondary`, `tertiary`), the polymer bonds and the broken bonds, the
clusters and the percolation; `output/polymer_topology.s0.json` lists the
broken bonds next to the formed ones, and the event frames in
`polymerization_events/s0/` (LAMMPS data files) carry the opened rings with
the ring bonds, angles and dihedrals through the broken bond left out (the
input script notes how many).

The ring torsions deserve a note: the connectivity of a three-membered ring
produces, besides the proper torsions around the ring bonds, the degenerate
torsions whose first and last atom coincide (`O_ep–CH2_ep–CH_ep–O_ep` and
the like); the molecule file lists these with a zero-parameter TraPPE
torsion, since every torsion of the connectivity has to be defined.
