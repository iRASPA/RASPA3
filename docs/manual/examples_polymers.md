# Examples Polymers
\page examples_polymers Examples Polymers


## Table of Contents
1. [Monte Carlo: cross-linked telechelic chains in box (reversible network)](#Example_polymers_1)


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
and thereby the degree of cross-linking. The two topology moves are enabled
per system with `"CrossLinkFormationProbability"` (formation/scission) and
`"CrossLinkSwapProbability"` (rewiring two links, keeping the number of links
constant). Linked molecules are still translated, rotated and regrown: the
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
acceptance of the formation/scission and swap moves alongside the ordinary
moves. Making `"FormationEnergy"` more negative drives the system towards a
fully linked network; setting it to zero leaves the association to the bond
well alone. The links are stored in the restart file, so a run can be
continued from a linked state.
