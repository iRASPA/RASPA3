# Monte Carlo moves
\page mc_moves Monte Carlo moves

This page describes all Monte Carlo moves implemented in `RASPA3` (see
`src/raspakit/mc_moves/move_types.ixx`). For every move it lists what the move
does, the detailed algorithmic steps, the acceptance rule, the type of system
and ensemble for which the move is useful, and a reference to the original
paper. The input keyword that activates the move is given in each section; see
the [Commands](\ref commands) page for the full input documentation.

## Table of Contents

<!-- TOC -->
* [Move selection](#move-selection)
* [Notation](#notation)
* [Single-molecule moves](#single-molecule-moves)
  * [Translation](#translation)
  * [Random translation](#random-translation)
  * [Rotation](#rotation)
  * [Random rotation](#random-rotation)
  * [Force-bias translation (smart MC)](#force-bias-translation)
* [Configurational-bias moves](#configurational-bias-moves)
  * [Reinsertion (CBMC)](#reinsertion-cbmc)
  * [Partial reinsertion (CBMC)](#partial-reinsertion-cbmc)
  * [Identity change (CBMC)](#identity-change-cbmc)
* [Insertion and deletion moves](#insertion-and-deletion-moves)
  * [Swap (conventional)](#swap-conventional)
  * [Swap (CBMC)](#swap-cbmc)
  * [Swap (CFCMC)](#swap-cfcmc)
  * [Swap (CB/CFCMC)](#swap-cbcfcmc)
  * [Pair swap (conventional)](#pair-swap-conventional)
  * [Pair swap (CBMC)](#pair-swap-cbmc)
  * [Pair swap (CFCMC)](#pair-swap-cfcmc)
  * [Pair swap (CB/CFCMC)](#pair-swap-cbcfcmc)
* [System moves](#system-moves)
  * [Volume change](#volume-change)
  * [Anisotropic volume change](#anisotropic-volume-change)
  * [Hybrid MC](#hybrid-mc)
  * [Force-bias translation of all molecules](#force-bias-translation-all)
* [Gibbs ensemble moves](#gibbs-ensemble-moves)
  * [Gibbs volume change](#gibbs-volume-change)
  * [Gibbs swap (CBMC)](#gibbs-swap-cbmc)
  * [Gibbs swap (CFCMC, serial)](#gibbs-swap-cfcmc)
  * [Gibbs swap (CB/CFCMC, serial)](#gibbs-swap-cbcfcmc)
  * [Gibbs conventional CFCMC (parallel)](#gibbs-conventional-cfcmc)
  * [Gibbs conventional CB/CFCMC (parallel)](#gibbs-conventional-cbcfcmc)
  * [Gibbs identity change (CBMC)](#gibbs-identity-change-cbmc)
* [Chemical-potential measurement moves](#chemical-potential-measurement-moves)
  * [Widom insertion](#widom-insertion)
  * [Widom CFCMC](#widom-cfcmc)
  * [Widom CB/CFCMC](#widom-cbcfcmc)
* [Parallel tempering](#parallel-tempering)
* [Reaction ensemble moves](#reaction-ensemble-moves)
  * [Reaction (CBMC)](#reaction-cbmc)
  * [Reaction conventional CFCMC (parallel)](#reaction-conventional-cfcmc)
  * [Reaction conventional CB/CFCMC (parallel)](#reaction-conventional-cbcfcmc)
  * [Reaction CFCMC (serial)](#reaction-cfcmc)
  * [Reaction CB/CFCMC (serial)](#reaction-cbcfcmc)
<!-- TOC -->

----------------------------------------------------------------------------------

## Move selection <a name="move-selection"></a>

A Monte Carlo cycle in `RASPA3` consists of a number of Monte Carlo steps (at
minimum the number of molecules present). Each step proceeds as:

- select a random system (box),
- select a random component of that system,
- sample a move type from the (normalized) move-probability table of that
  component; the probabilities are set in the input file with the keywords
  listed below,
- if the move acts on a single molecule: select a random molecule of that
  component in the chosen box,
- perform the move and accept or reject the trial configuration according to
  the acceptance rule of the move.

System moves (volume changes, hybrid MC, parallel tempering, reactions) are
also sampled through this mechanism but ignore the selected molecule and act on
the box(es) as a whole. The maximum displacement/rotation angle/volume change
of the adaptive moves is optimized during initialization and equilibration to
reach an acceptance ratio of approximately 50%.

## Notation <a name="notation"></a>

- \f$\beta = 1/(k_B T)\f$: inverse temperature,
- \f$\Delta U = U(\mathbf{n}) - U(\mathbf{o})\f$: potential-energy difference between the new
  (trial) configuration \f$\mathbf{n}\f$ and the old configuration \f$\mathbf{o}\f$,
- \f$N\f$: number of (integer) molecules of the selected component in the box,
- \f$V\f$: volume of the box,
- \f$f = \gamma\, x\, p\f$: fugacity of the component (fugacity coefficient
  \f$\gamma\f$, mole fraction \f$x\f$, pressure \f$p\f$),
- \f$W\f$: Rosenbluth weight of a CBMC-grown chain, \f$W^\text{IG}\f$ the ideal-gas
  Rosenbluth weight of an isolated molecule,
- \f$\lambda \in [0,1]\f$: coupling parameter of the fractional molecule in CFCMC
  moves; \f$\eta(\lambda)\f$ the corresponding bias factor (obtained by Wang-Landau
  weighting during equilibration),
- all acceptance rules are of the Metropolis form
  \f$\text{acc}(\mathbf{o} \rightarrow \mathbf{n}) = \min(1, \text{expression})\f$.

Trial configurations that place a molecule inside a blocked pocket are rejected
immediately.

----------------------------------------------------------------------------------

## Single-molecule moves <a name="single-molecule-moves"></a>

### Translation <a name="translation"></a>

Keyword: `"TranslationProbability"` (component move)

Displaces a single rigid-body molecule over a small random distance along one
of the Cartesian axes, leaving its orientation and internal conformation
unchanged. The maximum displacement is adjusted per component and per axis
during the run to obtain approximately 50% acceptance.

Steps:

- select a random box,
- select a random component,
- select a random molecule of that component in the chosen box,
- select a random Cartesian direction (x, y, or z),
- draw a random displacement \f$\Delta r\f$ uniformly from
  \f$[-\Delta_\text{max}, +\Delta_\text{max}]\f$ along the selected direction and translate the
  whole molecule,
- compute the energy difference \f$\Delta U\f$ (framework-molecule,
  molecule-molecule, Ewald Fourier, external field),
- accept or reject with the Metropolis rule.

Acceptance rule:

\f[
\text{acc}(\mathbf{o} \rightarrow \mathbf{n}) = \min\left(1, e^{-\beta \Delta U}\right)
\f]

Useful for: virtually every simulation; it is the basic move for thermal
equilibration of adsorbates in frameworks and of molecules in liquid or gas
boxes. Inexpensive because only an energy difference of a single molecule is
required.

Reference: N. Metropolis, A.W. Rosenbluth, M.N. Rosenbluth, A.H. Teller, and
E. Teller, *Equation of State Calculations by Fast Computing Machines*,
J. Chem. Phys. **21**, 1087-1092 (1953).
<https://doi.org/10.1063/1.1699114>

### Random translation <a name="random-translation"></a>

Keyword: `"RandomTranslationProbability"` (component move)

Translates a single molecule to an arbitrary position along one axis of the
box: the displacement is drawn over the full box length instead of a small
tuned maximum displacement. In contrast to reinsertion, the orientation and
internal conformation are kept, and no biasing is used.

Steps:

- select a random box,
- select a random component,
- select a random molecule of that component in the chosen box,
- select a random Cartesian axis of the cell,
- draw a fractional coordinate \f$s \in [0,1)\f$ along that axis and translate the
  molecule by the corresponding cell vector fraction,
- compute \f$\Delta U\f$ and apply the Metropolis rule.

Acceptance rule:

\f[
\text{acc}(\mathbf{o} \rightarrow \mathbf{n}) = \min\left(1, e^{-\beta \Delta U}\right)
\f]

Useful for: systems where molecules are trapped in local free-energy minima
(e.g. strongly adsorbing sites or channels in nanoporous materials) and normal
translations diffuse too slowly through phase space. Acceptance is usually low,
so it complements (not replaces) the regular translation move.

Reference: D. Dubbeldam, A. Torres-Knoop, and K.S. Walton, *On the Inner
Workings of Monte Carlo Codes*, Mol. Simul. **39**, 1253-1292 (2013).
<https://doi.org/10.1080/08927022.2013.819102>

### Rotation <a name="rotation"></a>

Keyword: `"RotationProbability"` (component move)

Rotates a single molecule over a small random angle around one of the Cartesian
axes through its rotation center (the starting bead). The maximum angle is
adjusted during the run towards 50% acceptance.

Steps:

- select a random box,
- select a random component,
- select a random molecule of that component in the chosen box,
- select a random rotation axis (x, y, or z),
- draw a random angle \f$\Delta\theta\f$ uniformly from
  \f$[-\theta_\text{max}, +\theta_\text{max}]\f$ and rotate the molecule rigidly around
  the axis,
- compute \f$\Delta U\f$ and apply the Metropolis rule.

Acceptance rule:

\f[
\text{acc}(\mathbf{o} \rightarrow \mathbf{n}) = \min\left(1, e^{-\beta \Delta U}\right)
\f]

Useful for: all systems containing non-spherical (multi-atom) molecules;
essential for equilibrating orientational degrees of freedom of rigid
molecules such as CO<sub>2</sub>, N<sub>2</sub>, water, or benzene.

Reference: J.A. Barker and R.O. Watts, *Structure of water; A Monte Carlo
calculation*, Chem. Phys. Lett. **3**, 144-145 (1969).
<https://doi.org/10.1016/0009-2614(69)80119-3>

### Random rotation <a name="random-rotation"></a>

Keyword: `"RandomRotationProbability"` (component move)

Gives a single molecule a completely new random orientation, drawn uniformly on
the rotation group (random vector on the unit sphere plus random angle),
independent of the current orientation.

Steps:

- select a random box,
- select a random component,
- select a random molecule of that component in the chosen box,
- generate a random new orientation (uniform random rotation),
- rotate the molecule around its rotation center to that orientation,
- compute \f$\Delta U\f$ and apply the Metropolis rule.

Acceptance rule:

\f[
\text{acc}(\mathbf{o} \rightarrow \mathbf{n}) = \min\left(1, e^{-\beta \Delta U}\right)
\f]

Useful for: molecules whose reorientation is severely hindered (e.g. dipolar
molecules in strong electrostatic fields, molecules in tight adsorption sites)
where sequential small rotations mix slowly.

Reference: D. Dubbeldam, A. Torres-Knoop, and K.S. Walton, *On the Inner
Workings of Monte Carlo Codes*, Mol. Simul. **39**, 1253-1292 (2013).
<https://doi.org/10.1080/08927022.2013.819102>

### Force-bias translation (smart MC) <a name="force-bias-translation"></a>

Keyword: `"ForceBiasTranslationProbability"` (component move)

A "smart Monte Carlo" displacement: the trial translation of a single molecule
is biased in the direction of the force acting on it, so that trial moves
preferentially go downhill in energy. The asymmetry of the proposal is removed
exactly by a Metropolis-Hastings correction, so the sampled distribution
remains the correct Boltzmann distribution.

Steps:

- select a random box,
- select a random component,
- select a random molecule of that component in the chosen box,
- compute the total (center-of-mass) force \f$\mathbf{F}_\mathbf{o}\f$ on the molecule in the
  current configuration,
- draw a Gaussian random vector \f$\boldsymbol{\xi}\f$ and displace the molecule rigidly over
  \f[
  \Delta\mathbf{r} = b\,\mathbf{F}_\mathbf{o} + \sigma\,\boldsymbol{\xi},
  \qquad b = \frac{1}{2}\beta\sigma^2 ,
  \f]
  where \f$\sigma\f$ is the tuned maximum-displacement parameter of the move,
- compute the energy difference \f$\Delta U\f$ and the force \f$\mathbf{F}_\mathbf{n}\f$ in the trial
  configuration,
- accept or reject including the Hastings correction for the asymmetric
  proposal.

Acceptance rule:

\f[
\text{acc}(\mathbf{o} \rightarrow \mathbf{n}) = \min\left(1,
   e^{-\beta \Delta U}\,
   \exp\left[\frac{\left|\Delta\mathbf{r} - b\mathbf{F}_\mathbf{o}\right|^2
                 - \left|\Delta\mathbf{r} + b\mathbf{F}_\mathbf{n}\right|^2}{2\sigma^2}\right]
   \right)
\f]

Useful for: dense systems (liquids, high loadings) where unbiased random
displacements have low acceptance; the force bias allows larger effective step
sizes at the cost of force evaluations.

References: C. Pangali, M. Rao, and B.J. Berne, *On a novel Monte Carlo scheme
for simulating water and aqueous solutions*, Chem. Phys. Lett. **55**, 413-417
(1978). <https://doi.org/10.1016/0009-2614(78)84003-2>;
P.J. Rossky, J.D. Doll, and H.L. Friedman, *Brownian dynamics as smart Monte
Carlo simulation*, J. Chem. Phys. **69**, 4628-4633 (1978).
<https://doi.org/10.1063/1.436415>

----------------------------------------------------------------------------------

## Configurational-bias moves <a name="configurational-bias-moves"></a>

### Reinsertion (CBMC) <a name="reinsertion-cbmc"></a>

Keyword: `"ReinsertionProbability"` (component move)

Removes a molecule from its current position and regrows it at a random
position in the box using Configurational-Bias Monte Carlo (CBMC). Multiple
trial positions are generated for the first bead and for every subsequent
segment; one is selected according to its Boltzmann weight. The old
configuration is "retraced" to compute its Rosenbluth weight, which enters the
acceptance rule and removes the bias exactly.

Steps:

- select a random box,
- select a random component,
- select a random molecule of that component in the chosen box,
- grow a new configuration of the molecule at a random position: place \f$k\f$
  trial positions for the first bead, select one with probability proportional
  to its Boltzmann factor, then grow the remaining segments segment-by-segment,
  each time generating \f$k\f$ trial orientations and selecting one according to
  its Boltzmann factor; accumulate the new Rosenbluth weight \f$W(\mathbf{n})\f$,
- retrace the old configuration with the same procedure (using the actual old
  positions as one of the trial sets) to obtain the old Rosenbluth weight
  \f$W(\mathbf{o})\f$,
- compute the remaining (Ewald Fourier, tail-correction) energy differences not
  included in the CBMC growth,
- accept or reject based on the ratio of Rosenbluth weights.

Acceptance rule:

\f[
\text{acc}(\mathbf{o} \rightarrow \mathbf{n}) = \min\left(1, \frac{W(\mathbf{n})}{W(\mathbf{o})}\,
     e^{-\beta \Delta U_\text{non-CBMC}}\right)
\f]

where \f$\Delta U_\text{non-CBMC}\f$ contains the energy contributions (Ewald Fourier and
tail corrections) not accounted for during the biased growth.

Useful for: flexible molecules (alkanes, alcohols, etc.), where it is often the
only efficient way to change the internal conformation; also useful for rigid
molecules to jump between adsorption sites without diffusion. Highly
recommended for chain molecules in confinement.

References: J.I. Siepmann and D. Frenkel, *Configurational bias Monte Carlo: a
new sampling scheme for flexible chains*, Mol. Phys. **75**, 59-70 (1992).
<https://doi.org/10.1080/00268979200100061>;
D. Frenkel, G.C.A.M. Mooij, and B. Smit, *Novel scheme to study structural and
thermal properties of continuously deformable molecules*, J. Phys.: Condens.
Matter **4**, 3053-3076 (1992). <https://doi.org/10.1088/0953-8984/4/12/006>

### Partial reinsertion (CBMC) <a name="partial-reinsertion-cbmc"></a>

Keyword: `"PartialReinsertionProbability"` (component move)

Regrows only a part of the molecule with CBMC while a chosen subset of atoms
remains fixed in place. The lists of fixed atoms per configuration are defined
in the component file (`partialReinsertionFixedAtoms`).

Steps:

- select a random box,
- select a random component,
- select a random molecule of that component in the chosen box,
- select at random one of the defined "fixed atom" configurations of the
  component,
- regrow all non-fixed atoms with CBMC starting from the fixed atoms,
  accumulating the new Rosenbluth weight \f$W(\mathbf{n})\f$,
- retrace the old positions of the regrown atoms to obtain \f$W(\mathbf{o})\f$,
- accept or reject based on the ratio of Rosenbluth weights.

Acceptance rule:

\f[
\text{acc}(\mathbf{o} \rightarrow \mathbf{n}) = \min\left(1, \frac{W(\mathbf{n})}{W(\mathbf{o})}\,
     e^{-\beta \Delta U_\text{non-CBMC}}\right)
\f]

Useful for: large or strongly adsorbed flexible molecules where a full
reinsertion has a very low acceptance; e.g. regrowing side chains or tail
segments while keeping an anchoring group in place.

Reference: T.J.H. Vlugt, R. Krishna, and B. Smit, *Molecular Simulations of
Adsorption Isotherms for Linear and Branched Alkanes and Their Mixtures in
Silicalite*, J. Phys. Chem. B **103**, 1102-1118 (1999).
<https://doi.org/10.1021/jp982736c>

### Identity change (CBMC) <a name="identity-change-cbmc"></a>

Keyword: `"IdentityChangeProbability"` (component move), together with
`"IdentityChangesList"`

Transforms a molecule of one component into a molecule of another component at
(approximately) the same position, using CBMC to grow the new identity. This
is the natural move of the semi-grand ensemble and greatly accelerates
equilibration of mixture compositions.

Steps:

- select a random box,
- select a random component,
- select at random a partner component from the identity-change list,
- select a random molecule of the old component in the chosen box,
- grow a molecule of the new component with CBMC, placing its first bead at the
  position of the first bead of the old molecule; accumulate \f$W(\mathbf{n})\f$,
- retrace the old molecule to obtain \f$W(\mathbf{o})\f$,
- delete the old molecule and insert the new one upon acceptance.

Acceptance rule (change of a molecule of component \f$A\f$ into component \f$B\f$):

\f[
\text{acc}(\mathbf{o} \rightarrow \mathbf{n}) = \min\left(1,
   \frac{W(\mathbf{n})}{W(\mathbf{o})}\,
   \frac{f_B\, N_A}{f_A\, (N_B + 1)}\,
   e^{-\beta \Delta U_\text{non-CBMC}} \right)
\f]

Useful for: adsorption of mixtures (e.g. binary/multicomponent isotherms in
zeolites and MOFs), where direct swaps of the minority component are rare; the
identity change equilibrates the composition far more efficiently.

References: D.A. Kofke and E.D. Glandt, *Monte Carlo simulation of
multicomponent equilibria in a semigrand canonical ensemble*, Mol. Phys.
**64**, 1105-1131 (1988). <https://doi.org/10.1080/00268978800100743>;
M.G. Martin and J.I. Siepmann, *Predicting multicomponent phase equilibria and
free energies of transfer for alkanes by molecular simulation*, J. Am. Chem.
Soc. **119**, 8921-8924 (1997). <https://doi.org/10.1021/ja964218q>

----------------------------------------------------------------------------------

## Insertion and deletion moves <a name="insertion-and-deletion-moves"></a>

These moves impose chemical equilibrium between the box and an imaginary
reservoir at the imposed fugacity \f$f\f$; they are the defining moves of grand
canonical (\f$\mu VT\f$) Monte Carlo.

### Swap (conventional) <a name="swap-conventional"></a>

Keyword: `"SwapConventionalProbability"` (component move)

Inserts or deletes a molecule at a random position without any biasing.

Steps:

- select a random box,
- select a random component,
- with 50% probability attempt an insertion, otherwise a deletion,
- insertion: place the molecule at a uniformly random position with a random
  orientation, compute \f$\Delta U\f$,
- deletion: select a random (integer) molecule of the component in the chosen
  box and compute the energy of removing it,
- accept or reject.

Acceptance rules:

\f[
\text{acc}_\text{ins} = \min\left(1, \frac{\beta f V}{N+1}\, e^{-\beta \Delta U} \right),
\qquad
\text{acc}_\text{del} = \min\left(1, \frac{N}{\beta f V}\, e^{-\beta \Delta U} \right)
\f]

Useful for: grand-canonical adsorption of small rigid molecules (CH<sub>4</sub>,
CO<sub>2</sub>, N<sub>2</sub>, noble gases) at low to moderate densities. For chain
molecules or dense systems use the CBMC or CFCMC variants instead.

Reference: D.J. Adams, *Grand canonical ensemble Monte Carlo for a
Lennard-Jones fluid*, Mol. Phys. **29**, 307-311 (1975).
<https://doi.org/10.1080/00268977500100221>

### Swap (CBMC) <a name="swap-cbmc"></a>

Keyword: `"SwapProbability"` (component move)

Grand-canonical insertion/deletion using configurational-bias growth: molecules
are grown atom-by-atom with multiple trial directions, which raises the
insertion acceptance by orders of magnitude for chain molecules and in
confined systems.

Steps:

- select a random box,
- select a random component,
- with 50% probability attempt an insertion, otherwise a deletion,
- insertion: grow a new molecule with CBMC at a random position (multiple trial
  first beads and trial segment orientations) and accumulate the Rosenbluth
  weight \f$W(\mathbf{n})\f$,
- deletion: select a random (integer) molecule of the component in the chosen
  box and retrace it with CBMC to obtain \f$W(\mathbf{o})\f$,
- add the non-CBMC energy contributions (Ewald Fourier, tail corrections),
- accept or reject.

Acceptance rules (with the ideal-gas Rosenbluth weight \f$W^\text{IG}\f$ of the
isolated molecule as normalization):

\f[
\text{acc}_\text{ins} = \min\left(1, \frac{\beta f V}{N+1}\,
      \frac{W(\mathbf{n})}{W^\text{IG}}\, e^{-\beta \Delta U_\text{non-CBMC}} \right),
\qquad
\text{acc}_\text{del} = \min\left(1, \frac{N}{\beta f V}\,
      \frac{W^\text{IG}}{W(\mathbf{o})}\, e^{-\beta \Delta U_\text{non-CBMC}} \right)
\f]

Useful for: adsorption isotherms of flexible/chain molecules (alkanes,
alcohols) in nanoporous materials; the standard workhorse for GCMC of anything
larger than a few atoms.

References: B. Smit, *Grand canonical Monte Carlo simulations of chain
molecules: adsorption isotherms of alkanes in zeolites*, Mol. Phys. **85**,
153-172 (1995). <https://doi.org/10.1080/00268979500101011>;
T.J.H. Vlugt, R. Krishna, and B. Smit, J. Phys. Chem. B **103**, 1102-1118
(1999). <https://doi.org/10.1021/jp982736c>

### Swap (CFCMC) <a name="swap-cfcmc"></a>

Keyword: `"CFCMC_SwapProbability"` (component move)

Continuous Fractional Component Monte Carlo insertion/deletion. Each component
carries one *fractional* molecule whose interactions with the surroundings are
scaled by a coupling parameter \f$\lambda \in [0,1]\f$ (\f$\lambda=0\f$: ideal gas /
non-interacting, \f$\lambda=1\f$: fully interacting). Instead of inserting a whole
molecule at once, the move gradually changes \f$\lambda\f$; a biasing potential
\f$\eta(\lambda)\f$ (measured with Wang-Landau updating during equilibration)
flattens the \f$\lambda\f$ histogram so that insertions and deletions become
possible even in dense systems.

Steps:

- select a random box,
- select a random component,
- draw a random change \f$\Delta\lambda\f$ of the coupling parameter of the fractional
  molecule (the maximum change is tuned to 50% acceptance); three cases arise:
- if \f$\lambda + \Delta\lambda\f$ stays within \f$[0,1]\f$ — **\f$\lambda\f$-change**: only the scaling of
  the fractional molecule changes; compute \f$\Delta U\f$ and accept with the biased
  Metropolis rule,
- if \f$\lambda + \Delta\lambda > 1\f$ — **insertion**: the fractional molecule becomes a fully
  coupled integer molecule (\f$\lambda \rightarrow 1\f$), and a new fractional molecule is
  inserted at a random position with the remainder \f$\lambda_\mathbf{n} = \lambda + \Delta\lambda - 1\f$,
- if \f$\lambda + \Delta\lambda < 0\f$ — **deletion**: the current fractional molecule is removed
  (\f$\lambda \rightarrow 0\f$), and a randomly chosen integer molecule of the same component
  becomes the new fractional molecule with \f$\lambda_\mathbf{n} = \lambda + \Delta\lambda + 1\f$,
- accept or reject with the appropriate rule below.

Acceptance rules (\f$\Delta\eta = \eta(\lambda_\mathbf{n}) - \eta(\lambda_\mathbf{o})\f$):

\f[
\text{acc}_{\lambda} = \min\left(1, e^{-\beta \Delta U + \Delta\eta}\right),
\qquad
\text{acc}_\text{ins} = \min\left(1, \frac{\beta f V}{N+1}\,
   e^{-\beta \Delta U + \Delta\eta}\right),
\qquad
\text{acc}_\text{del} = \min\left(1, \frac{N}{\beta f V}\,
   e^{-\beta \Delta U + \Delta\eta}\right)
\f]

Useful for: open-ensemble simulations of dense phases where direct insertions
essentially never succeed: high-loading adsorption, liquids, ionic liquids. As
a by-product the chemical potential is obtained directly from the \f$\lambda\f$
histogram.

Reference: W. Shi and E.J. Maginn, *Continuous Fractional Component Monte
Carlo: An Adaptive Biasing Method for Open System Atomistic Simulations*,
J. Chem. Theory Comput. **3**, 1451-1463 (2007).
<https://doi.org/10.1021/ct7000039>

### Swap (CB/CFCMC) <a name="swap-cbcfcmc"></a>

Keyword: `"CFCMC_CBMC_SwapProbability"` (component move)

Hybrid of CBMC and CFCMC: the insertion/deletion of the fractional molecule is
performed with configurational-bias growth (multiple trial positions and
orientations), while the coupling strength is still changed gradually through
\f$\lambda\f$. This combines the orientational/conformational pre-sampling of CBMC
with the smooth coupling of CFCMC.

Steps:

- select a random box,
- select a random component,
- draw a random change \f$\Delta\lambda\f$ of the coupling parameter of the fractional
  molecule; three cases arise:
- if \f$\lambda + \Delta\lambda\f$ stays within \f$[0,1]\f$ — **\f$\lambda\f$-change**: rescale the
  interactions of the fractional molecule and compute \f$\Delta U\f$,
- if \f$\lambda + \Delta\lambda > 1\f$ — **insertion**: the fractional molecule becomes a fully
  coupled integer molecule (\f$\lambda \rightarrow 1\f$), and a new fractional molecule with
  \f$\lambda_\mathbf{n} = \lambda + \Delta\lambda - 1\f$ is grown with CBMC (multiple trial positions and
  orientations), accumulating the Rosenbluth weight \f$W(\mathbf{n})\f$,
- if \f$\lambda + \Delta\lambda < 0\f$ — **deletion**: the fractional molecule is retraced with
  CBMC (Rosenbluth weight \f$W(\mathbf{o})\f$) and removed, and a randomly selected
  integer molecule becomes the new fractional molecule with
  \f$\lambda_\mathbf{n} = \lambda + \Delta\lambda + 1\f$,
- accept or reject.

Acceptance rules (\f$\Delta\eta = \eta(\lambda_\mathbf{n}) - \eta(\lambda_\mathbf{o})\f$):

\f[
\text{acc}_{\lambda} = \min\left(1, e^{-\beta \Delta U + \Delta\eta}\right)
\f]

\f[
\text{acc}_\text{ins} = \min\left(1, \frac{\beta f V}{N+1}\,
   \frac{W(\mathbf{n})}{W^\text{IG}}\,
   e^{-\beta \Delta U + \Delta\eta}\right),
\qquad
\text{acc}_\text{del} = \min\left(1, \frac{N}{\beta f V}\,
   \frac{W^\text{IG}}{W(\mathbf{o})}\,
   e^{-\beta \Delta U + \Delta\eta}\right)
\f]

Useful for: open-ensemble simulations of *flexible* molecules in dense phases
(long alkanes at high loading, chain molecules in liquids), where both the
conformational bias and the gradual coupling are needed.

Reference: A. Torres-Knoop, S.P. Balaji, T.J.H. Vlugt, and D. Dubbeldam, *A
Comparison of Advanced Monte Carlo Methods for Open Systems: CFCMC vs CBMC*,
J. Chem. Theory Comput. **10**, 942-952 (2014).
<https://doi.org/10.1021/ct4009766>

### Pair swap (conventional) <a name="pair-swap-conventional"></a>

Keyword: `"PairSwapConventionalProbability"` (component move), together with
the pair-partner component and `"MaximumPairDistance"` \f$R_\text{max}\f$

Inserts or deletes a *pair* of molecules (the selected component \f$A\f$ together
with its defined partner component \f$B\f$) in a single move. The first molecule
is placed at a random position; the second molecule is placed with its first
bead uniformly inside a sphere of radius \f$R_\text{max}\f$ around the first bead of
molecule \f$A\f$ (the distance is drawn as \f$r = R_\text{max}\sqrt[3]{u}\f$, uniform in
the sphere volume). Keeping the pair together at insertion is essential for
oppositely charged species.

Steps:

- select a random box,
- select a random component \f$A\f$ (the pair partner \f$B\f$ is defined in the
  input),
- with 50% probability attempt a pair insertion, otherwise a pair deletion,
- insertion: grow molecule \f$A\f$ at a random position; draw a distance
  \f$r = R_\text{max}\sqrt[3]{u}\f$ and a random direction on the unit sphere, and grow
  molecule \f$B\f$ with its first bead at that position,
- deletion: select a random integer molecule of component \f$A\f$ and its paired
  molecule of component \f$B\f$, and retrace both,
- compute the energy differences and accept or reject.

Acceptance rules:

\f[
\text{acc}_\text{ins} = \min\left(1,
   \frac{\beta f_A f_B V}{(N_A+1)(N_B+1)}\,
   \frac{W_A(\mathbf{n})}{W_A^\text{IG}}\,\frac{W_B(\mathbf{n})}{W_B^\text{IG}}\,
   e^{-\beta \Delta U_\text{non-CBMC}}\right)
\f]

\f[
\text{acc}_\text{del} = \min\left(1,
   \frac{N_A N_B}{\beta f_A f_B V}\,
   \frac{W_A^\text{IG}}{W_A(\mathbf{o})}\,\frac{W_B^\text{IG}}{W_B(\mathbf{o})}\,
   e^{-\beta \Delta U_\text{non-CBMC}}\right)
\f]

A single volume factor \f$V\f$ appears because molecule \f$B\f$ is confined to the
pair sphere around molecule \f$A\f$ rather than inserted anywhere in the box.

Useful for: systems where single-molecule swaps would break a constraint that
must hold at all times, most notably charge neutrality when inserting ionic
species (salt ion pairs, charged adsorbates with counter-ions).

Reference: G. Orkoulas and A.Z. Panagiotopoulos, *Free energy and phase
equilibria for the restricted primitive model of ionic fluids from Monte Carlo
simulations*, J. Chem. Phys. **101**, 1452-1459 (1994).
<https://doi.org/10.1063/1.467770>

### Pair swap (CBMC) <a name="pair-swap-cbmc"></a>

Keyword: `"PairSwapProbability"` (component move)

The configurational-bias variant of the pair swap: both molecules of the pair
are grown with CBMC (multiple trial positions and orientations). The pair
distance is drawn uniformly in \f$r \in [0, R_\text{max}]\f$, which introduces an
explicit distance-bias factor \f$b(r) = 3r^2/R_\text{max}^2\f$ in the acceptance rule
(the reverse of the \f$r^2\f$ volume element of the spherical shell).

Steps:

- select a random box,
- select a random component \f$A\f$ (the pair partner \f$B\f$ is defined in the
  input),
- with 50% probability attempt a pair insertion, otherwise a pair deletion,
- insertion: grow molecule \f$A\f$ with CBMC at a random position; draw
  \f$r = R_\text{max}\, u\f$ and a random direction, and grow molecule \f$B\f$ with CBMC
  with its first bead fixed at that position,
- deletion: select a random integer molecule of \f$A\f$ and its paired molecule of
  \f$B\f$ and retrace both with CBMC,
- accept or reject.

Acceptance rules:

\f[
\text{acc}_\text{ins} = \min\left(1,
   \frac{\beta f_A f_B V}{(N_A+1)(N_B+1)}\;
   \frac{3 r^2}{R_\text{max}^2}\;
   \frac{W_A(\mathbf{n})}{W_A^\text{IG}}\,\frac{W_B(\mathbf{n})}{W_B^\text{IG}}\,
   e^{-\beta \Delta U_\text{non-CBMC}}\right)
\f]

\f[
\text{acc}_\text{del} = \min\left(1,
   \frac{N_A N_B}{\beta f_A f_B V}\;
   \frac{R_\text{max}^2}{3 r^2}\;
   \frac{W_A^\text{IG}}{W_A(\mathbf{o})}\,\frac{W_B^\text{IG}}{W_B(\mathbf{o})}\,
   e^{-\beta \Delta U_\text{non-CBMC}}\right)
\f]

Useful for: pair insertions of *flexible* charged species (e.g. molecular ions
with several beads) where the conventional pair placement has too low an
acceptance.

Reference: G. Orkoulas and A.Z. Panagiotopoulos, J. Chem. Phys. **101**,
1452-1459 (1994). <https://doi.org/10.1063/1.467770>

### Pair swap (CFCMC) <a name="pair-swap-cfcmc"></a>

Keyword: `"CFCMC_PairSwapProbability"` (component move)

The CFCMC variant of the pair swap: the pair of fractional molecules (one of
component \f$A\f$, one of component \f$B\f$) is coupled through a *single, common*
coupling parameter \f$\lambda\f$, so both molecules appear and disappear together and
the system stays charge neutral at every value of \f$\lambda\f$.

Steps:

- select a random box,
- select a random component \f$A\f$ (the pair partner \f$B\f$ is defined in the
  input),
- draw a random change \f$\Delta\lambda\f$ of the common coupling parameter; three cases
  arise:
- if \f$\lambda + \Delta\lambda\f$ stays within \f$[0,1]\f$ — **\f$\lambda\f$-change**: rescale both
  fractional molecules simultaneously,
- if \f$\lambda + \Delta\lambda > 1\f$ — **pair insertion**: both fractional molecules become
  integer molecules and a new fractional pair is inserted at random positions
  with \f$\lambda_\mathbf{n} = \lambda + \Delta\lambda - 1\f$,
- if \f$\lambda + \Delta\lambda < 0\f$ — **pair deletion**: the fractional pair is removed and
  a randomly selected integer molecule of each component becomes the new
  fractional pair with \f$\lambda_\mathbf{n} = \lambda + \Delta\lambda + 1\f$,
- accept or reject.

Acceptance rules (\f$\Delta\eta = \eta(\lambda_\mathbf{n}) - \eta(\lambda_\mathbf{o})\f$):

\f[
\text{acc}_{\lambda} = \min\left(1, e^{-\beta \Delta U + \Delta\eta}\right)
\f]

\f[
\text{acc}_\text{ins} = \min\left(1,
   \frac{\beta f_A V}{N_A+1}\,\frac{\beta f_B V}{N_B+1}\,
   e^{-\beta \Delta U + \Delta\eta}\right),
\qquad
\text{acc}_\text{del} = \min\left(1,
   \frac{N_A}{\beta f_A V}\,\frac{N_B}{\beta f_B V}\,
   e^{-\beta \Delta U + \Delta\eta}\right)
\f]

Useful for: open-ensemble simulations of ionic species in dense phases (salts
in liquids, ionic adsorbates at high loading) where both charge neutrality and
gradual coupling are required.

References: W. Shi and E.J. Maginn, J. Chem. Theory Comput. **3**, 1451-1463
(2007). <https://doi.org/10.1021/ct7000039>; G. Orkoulas and A.Z.
Panagiotopoulos, J. Chem. Phys. **101**, 1452-1459 (1994).
<https://doi.org/10.1063/1.467770>

### Pair swap (CB/CFCMC) <a name="pair-swap-cbcfcmc"></a>

Keyword: `"CFCMC_CBMC_PairSwapProbability"` (component move)

The CFCMC pair swap with configurational-bias growth: in the pair-insertion
and pair-deletion branches the fractional molecules are grown/retraced with
CBMC, adding the Rosenbluth-weight ratios of both molecules to the acceptance
rules.

Acceptance rules:

\f[
\text{acc}_{\lambda} = \min\left(1, e^{-\beta \Delta U + \Delta\eta}\right)
\f]

\f[
\text{acc}_\text{ins} = \min\left(1,
   \frac{\beta f_A V}{N_A+1}\,\frac{\beta f_B V}{N_B+1}\,
   \frac{W_A(\mathbf{n})}{W_A^\text{IG}}\,\frac{W_B(\mathbf{n})}{W_B^\text{IG}}\,
   e^{-\beta \Delta U + \Delta\eta}\right)
\f]

\f[
\text{acc}_\text{del} = \min\left(1,
   \frac{N_A}{\beta f_A V}\,\frac{N_B}{\beta f_B V}\,
   \frac{W_A^\text{IG}}{W_A(\mathbf{o})}\,\frac{W_B^\text{IG}}{W_B(\mathbf{o})}\,
   e^{-\beta \Delta U + \Delta\eta}\right)
\f]

Useful for: flexible ionic species in dense phases, combining charge-neutral
pair transfer, gradual coupling, and conformational biasing.

References: A. Torres-Knoop, S.P. Balaji, T.J.H. Vlugt, and D. Dubbeldam,
J. Chem. Theory Comput. **10**, 942-952 (2014).
<https://doi.org/10.1021/ct4009766>; W. Shi and E.J. Maginn, J. Chem. Theory
Comput. **3**, 1451-1463 (2007). <https://doi.org/10.1021/ct7000039>

----------------------------------------------------------------------------------

## System moves <a name="system-moves"></a>

### Volume change <a name="volume-change"></a>

Keyword: `"VolumeMoveProbability"` (system move)

Isotropic change of the box volume at constant external pressure — the basic
move of the \f$NpT\f$ ensemble. The volume change is sampled in \f$\ln V\f$. Rigid
molecules are scaled by their center of mass (internal geometry preserved);
framework and flexible molecules are scaled atom-by-atom.

Steps:

- select a random box,
- draw a new volume from
  \f$\ln V_\mathbf{n} = \ln V_\mathbf{o} + \Delta_\text{max}\,(2u - 1)\f$ with \f$u\f$ uniform in \f$[0,1)\f$,
- scale the box vectors isotropically by \f$(V_\mathbf{n}/V_\mathbf{o})^{1/3}\f$ and scale all
  molecule positions (center-of-mass scaling for rigid molecules),
- recompute the total energy of the scaled configuration,
- accept or reject.

Acceptance rule (for sampling in \f$\ln V\f$):

\f[
\text{acc}(\mathbf{o} \rightarrow \mathbf{n}) = \min\left(1,
  \left(\frac{V_\mathbf{n}}{V_\mathbf{o}}\right)^{N+1}
  e^{-\beta\left[\Delta U + p\,(V_\mathbf{n} - V_\mathbf{o})\right]}\right)
\f]

Useful for: \f$NpT\f$ simulations of bulk fluids (densities, equations of state)
and for equilibrating framework densities of flexible frameworks. Expensive:
requires a full energy recalculation.

References: I.R. McDonald, *NpT-ensemble Monte Carlo calculations for binary
liquid mixtures*, Mol. Phys. **23**, 41-58 (1972).
<https://doi.org/10.1080/00268977200100031>;
R. Eppenga and D. Frenkel, *Monte Carlo study of the isotropic and nematic
phases of infinitely thin hard platelets*, Mol. Phys. **52**, 1303-1334 (1984).
<https://doi.org/10.1080/00268978400101951>

### Anisotropic volume change <a name="anisotropic-volume-change"></a>

Keyword: `"AnisotropicVolumeMoveProbability"` (system move)

Changes the three box lengths independently (the box shape can change while
the pressure tensor is imposed per diagonal component). Molecules are
translated rigidly with their center-of-mass fractional position kept fixed,
so intramolecular energies are unchanged.

Steps:

- select a random box,
- draw three independent scale factors \f$s_\alpha = e^{\Delta_{\text{max},\alpha}(2u_\alpha-1)}\f$
  for \f$\alpha = x,y,z\f$,
- scale the box vectors and the molecular center-of-mass positions accordingly,
- recompute the total energy,
- accept or reject using the volume ratio
  \f$V_\mathbf{n}/V_\mathbf{o} = s_x s_y s_z\f$ and the pressure-tensor work
  \f$\Delta W = V_\mathbf{o} \sum_\alpha p_{\alpha\alpha}\,(s_\alpha - 1)\f$.

Acceptance rule:

\f[
\text{acc}(\mathbf{o} \rightarrow \mathbf{n}) = \min\left(1,
  \left(\frac{V_\mathbf{n}}{V_\mathbf{o}}\right)^{N+1}
  e^{-\beta\left[\Delta U + \Delta W\right]}\right)
\f]

Useful for: solids and frameworks where the cell shape must relax (structural
transitions, breathing MOFs), and for fluids in anisotropic confinement where
the box aspect ratio matters.

Reference: S. Yashonath and C.N.R. Rao, *A Monte Carlo study of crystal
structure transformations*, Mol. Phys. **54**, 245-251 (1985).
<https://doi.org/10.1080/00268978500100201>

### Hybrid MC <a name="hybrid-mc"></a>

Keyword: `"HybridMCProbability"` (system move), with
`"HybridMCMoveNumberOfSteps"` for the trajectory length

Performs a short molecular-dynamics trajectory as a single collective Monte
Carlo move: all molecules move simultaneously, giving efficient collective
relaxation while the Metropolis step removes integration errors.

Steps:

- select a random box,
- draw fresh velocities for all molecules from the Maxwell-Boltzmann
  distribution at the system temperature,
- integrate Newton's equations of motion (velocity Verlet, NVE) for the chosen
  number of MD steps,
- measure the drift of the conserved (total) energy
  \f$\Delta E = E_\text{end} - E_\text{start}\f$ over the trajectory,
- accept the final configuration or restore the initial one.

Acceptance rule (as implemented, using the absolute conserved-energy drift):

\f[
\text{acc}(\mathbf{o} \rightarrow \mathbf{n}) = \min\left(1, e^{-\beta\, \left|\Delta E\right|}\right)
\f]

Useful for: dense fluids and flexible frameworks where single-molecule moves
relax the structure slowly; also the standard way to include framework
flexibility and collective motions in an MC simulation.

References: S. Duane, A.D. Kennedy, B.J. Pendleton, and D. Roweth, *Hybrid
Monte Carlo*, Phys. Lett. B **195**, 216-222 (1987).
<https://doi.org/10.1016/0370-2693(87)91197-X>;
B. Mehlig, D.W. Heermann, and B.M. Forrest, *Hybrid Monte Carlo method for
condensed-matter systems*, Phys. Rev. B **45**, 679-685 (1992).
<https://doi.org/10.1103/PhysRevB.45.679>

### Force-bias translation of all molecules <a name="force-bias-translation-all"></a>

Keyword: `"ForceBiasTranslationAllProbability"` (system move)

The collective variant of the [force-bias translation](#force-bias-translation):
*all* molecules in the box are displaced simultaneously, each biased along the
force acting on it, with a global Metropolis-Hastings correction. It behaves
like a single Brownian-dynamics step wrapped in an exact acceptance rule.

Steps:

- select a random box,
- compute the center-of-mass forces on all molecules,
- displace every molecule by \f$\Delta\mathbf{r}_i = b\,\mathbf{F}_{i,\mathbf{o}} + \sigma\,\boldsymbol{\xi}_i\f$
  with \f$b = \frac{1}{2}\beta\sigma^2\f$ and independent Gaussian vectors \f$\boldsymbol{\xi}_i\f$,
- recompute the energy and the forces in the trial configuration,
- accept or reject all displacements together with the Hastings-corrected
  Metropolis rule (the product over molecules of the single-molecule bias
  factors).

Acceptance rule:

\f[
\text{acc}(\mathbf{o} \rightarrow \mathbf{n}) = \min\left(1,
  e^{-\beta \Delta U}\,
  \exp\left[\sum_i \frac{\left|\Delta\mathbf{r}_i - b\mathbf{F}_{i,\mathbf{o}}\right|^2
     - \left|\Delta\mathbf{r}_i + b\mathbf{F}_{i,\mathbf{n}}\right|^2}{2\sigma^2}\right]\right)
\f]

Useful for: dense liquids and high loadings, as a cheaper alternative to
hybrid MC for collective translational relaxation.

Reference: P.J. Rossky, J.D. Doll, and H.L. Friedman, *Brownian dynamics as
smart Monte Carlo simulation*, J. Chem. Phys. **69**, 4628-4633 (1978).
<https://doi.org/10.1063/1.436415>

----------------------------------------------------------------------------------

## Gibbs ensemble moves <a name="gibbs-ensemble-moves"></a>

The Gibbs ensemble simulates two coupled boxes (typically a liquid and a vapor
phase) that exchange volume and molecules so that both phases are in mutual
thermodynamic equilibrium without an explicit interface.

### Gibbs volume change <a name="gibbs-volume-change"></a>

Keyword: `"GibbsVolumeMoveProbability"` (system move)

Exchanges volume between the two boxes at constant total volume
\f$V = V_A + V_B\f$; the change is sampled in \f$\ln(V_A/V_B)\f$.

Steps:

- draw \f$\ln(V_{A,\mathbf{n}}/V_{B,\mathbf{n}}) = \ln(V_{A,\mathbf{o}}/V_{B,\mathbf{o}}) + \Delta_\text{max}(2u-1)\f$,
- compute the two new volumes with \f$V_{A,\mathbf{n}} + V_{B,\mathbf{n}} = V\f$,
- scale both boxes isotropically and scale the molecule positions by their
  centers of mass,
- recompute the total energies of both boxes,
- accept or reject both scalings together.

Acceptance rule:

\f[
\text{acc}(\mathbf{o} \rightarrow \mathbf{n}) = \min\left(1,
  \left(\frac{V_{A,\mathbf{n}}}{V_{A,\mathbf{o}}}\right)^{N_A+1}
  \left(\frac{V_{B,\mathbf{n}}}{V_{B,\mathbf{o}}}\right)^{N_B+1}
  e^{-\beta\left[\Delta U_A + \Delta U_B\right]}\right)
\f]

Useful for: vapor-liquid equilibrium (VLE) calculations of pure components and
mixtures; mandatory in any Gibbs ensemble simulation (together with a Gibbs
swap move).

Reference: A.Z. Panagiotopoulos, *Direct determination of phase coexistence
properties of fluids by Monte Carlo simulation in a new ensemble*, Mol. Phys.
**61**, 813-826 (1987). <https://doi.org/10.1080/00268978700101491>

### Gibbs swap (CBMC) <a name="gibbs-swap-cbmc"></a>

Keyword: `"GibbsSwapCBMCProbability"` (component move)

Transfers one molecule from one box to the other using CBMC growth in the
receiving box and CBMC retracing in the donating box. This move equalizes the
chemical potential of the component between the two phases.

Steps:

- choose the transfer direction (box \f$A \rightarrow B\f$ or \f$B \rightarrow A\f$, 50% each),
- select a random molecule of the selected component in the donating box,
- grow a trial molecule with CBMC at a random position in the receiving box and
  accumulate \f$W(\mathbf{n})\f$,
- retrace the molecule in the donating box to obtain \f$W(\mathbf{o})\f$,
- accept or reject the simultaneous insertion and deletion.

Acceptance rule (transfer from box \f$B\f$ to box \f$A\f$):

\f[
\text{acc}(\mathbf{o} \rightarrow \mathbf{n}) = \min\left(1,
  \frac{W(\mathbf{n})\, N_B\, V_A}{W(\mathbf{o})\, (N_A+1)\, V_B}\,
  e^{-\beta \Delta U_\text{non-CBMC}}\right)
\f]

Useful for: VLE and liquid-liquid equilibria of molecular fluids and mixtures;
the CBMC growth keeps acceptance workable for chain molecules.

References: A.Z. Panagiotopoulos, N. Quirke, M. Stapleton, and D.J. Tildesley,
*Phase equilibria by simulation in the Gibbs ensemble: alternative derivation,
generalization and application to mixture and membrane equilibria*, Mol. Phys.
**63**, 527-545 (1988). <https://doi.org/10.1080/00268978800100361>;
G.C.A.M. Mooij, D. Frenkel, and B. Smit, *Direct simulation of phase
equilibria of chain molecules*, J. Phys.: Condens. Matter **4**, L255-L259
(1992). <https://doi.org/10.1088/0953-8984/4/16/001>

### Gibbs swap (CFCMC, serial) <a name="gibbs-swap-cfcmc"></a>

Keyword: `"GibbsSwapCFCMCProbability"` (component move)

Serial CFCMC version of the Gibbs swap: there is a *single* fractional molecule
in the whole two-box system, and the coupling parameter \f$\lambda\f$ interpolates the
molecule between the two boxes. Increasing \f$\lambda\f$ couples the fractional
molecule in one box while decoupling it in the other, so molecule transfer
proceeds gradually instead of in one step.

Steps (box \f$A\f$ is the box currently holding the fractional molecule; each
box keeps its own bias function \f$\eta_A(\lambda)\f$, \f$\eta_B(\lambda)\f$):

- identify the box \f$A\f$ holding the fractional molecule,
- select one of three sub-moves at random: with 25% probability a **swap**
  sub-move, with 25% probability a **transfer** sub-move, and with 50%
  probability a **\f$\lambda\f$-change** sub-move,
- **swap** sub-move: the fractional molecule in box \f$A\f$ becomes a whole
  (integer) molecule at its current position, and a randomly selected integer
  molecule in box \f$B\f$ becomes the new fractional molecule at the *same*
  \f$\lambda\f$; compute the energy differences in both boxes,
- **transfer** sub-move: the fractional molecule is moved (at the same \f$\lambda\f$)
  from box \f$A\f$ to a random position in box \f$B\f$; compute the energy
  differences in both boxes,
- **\f$\lambda\f$-change** sub-move: draw a new \f$\lambda\f$ bin and rescale the fractional
  molecule inside box \f$A\f$,
- accept or reject with the corresponding rule below.

Acceptance rules (integer-molecule counts \f$N_A\f$, \f$N_B\f$; box volumes \f$V_A\f$,
\f$V_B\f$):

\f[
\text{acc}_\text{swap} = \min\left(1,
  \frac{N_B}{N_A+1}\,
  e^{-\beta\left[\Delta U_A + \Delta U_B\right]
     + \eta_B(\lambda) - \eta_A(\lambda)}\right)
\f]

\f[
\text{acc}_\text{transfer} = \min\left(1,
  \frac{V_B}{V_A}\,
  e^{-\beta\left[\Delta U_A + \Delta U_B\right]
     + \eta_B(\lambda) - \eta_A(\lambda)}\right)
\f]

\f[
\text{acc}_{\lambda} = \min\left(1,
  e^{-\beta \Delta U + \eta_A(\lambda_\mathbf{n}) - \eta_A(\lambda_\mathbf{o})}\right)
\f]

Useful for: VLE of dense and strongly interacting fluids where direct Gibbs
swaps fail (e.g. water, ionic liquids, high-density liquids near the triple
point). Also yields the chemical potential of both phases directly.

Reference: A. Poursaeidesfahani, A. Torres-Knoop, D. Dubbeldam, and T.J.H.
Vlugt, *Direct Free Energy Calculation in the Continuous Fractional Component
Gibbs Ensemble*, J. Chem. Theory Comput. **12**, 1481-1490 (2016).
<https://doi.org/10.1021/acs.jctc.5b01230>

### Gibbs swap (CB/CFCMC, serial) <a name="gibbs-swap-cbcfcmc"></a>

Keyword: `"GibbsSwapCBCFCMCProbability"` (component move)

The serial CFCMC Gibbs swap combined with CBMC: when the fractional molecule
hops between boxes it is grown/retraced with configurational-bias sampling, so
the Rosenbluth-weight ratio enters the acceptance rule in addition to the
CFCMC bias terms.

Steps: as for [Gibbs swap (CFCMC)](#gibbs-swap-cfcmc) — the same three
sub-moves (25% swap, 25% transfer, 50% \f$\lambda\f$-change) — but in the swap and
transfer sub-moves the new fractional molecule is grown with CBMC in the
receiving box (Rosenbluth weight \f$W(\mathbf{n})\f$) and the old fractional molecule
is retraced with CBMC in the donating box (Rosenbluth weight \f$W(\mathbf{o})\f$).

Acceptance rules:

\f[
\text{acc}_\text{swap} = \min\left(1,
  \frac{N_B}{N_A+1}\,
  \frac{W(\mathbf{n})}{W(\mathbf{o})}\,
  e^{-\beta\left[\Delta U_A + \Delta U_B\right]
     + \eta_B(\lambda) - \eta_A(\lambda)}\right)
\f]

\f[
\text{acc}_\text{transfer} = \min\left(1,
  \frac{V_B}{V_A}\,
  \frac{W(\mathbf{n})}{W(\mathbf{o})}\,
  e^{-\beta\left[\Delta U_A + \Delta U_B\right]
     + \eta_B(\lambda) - \eta_A(\lambda)}\right)
\f]

\f[
\text{acc}_{\lambda} = \min\left(1,
  e^{-\beta \Delta U + \eta_A(\lambda_\mathbf{n}) - \eta_A(\lambda_\mathbf{o})}\right)
\f]

Useful for: Gibbs-ensemble phase equilibria of flexible chain molecules in
dense phases.

Reference: A. Torres-Knoop, S.P. Balaji, T.J.H. Vlugt, and D. Dubbeldam,
J. Chem. Theory Comput. **10**, 942-952 (2014).
<https://doi.org/10.1021/ct4009766>

### Gibbs conventional CFCMC (parallel) <a name="gibbs-conventional-cfcmc"></a>

Keyword: `"GibbsConventionalCFCMCProbability"` (component move)

Parallel CFCMC version of the Gibbs ensemble in the spirit of Shi and Maginn:
*each box* carries its own fractional molecule of the component, with coupling
parameters \f$\lambda_A\f$ and \f$\lambda_B\f$ and bias functions \f$\eta_A\f$, \f$\eta_B\f$. A single
random change \f$\nu\f$ moves the two coupling parameters in opposite directions,
\f$\lambda_A \rightarrow \lambda_A + \nu\f$ and \f$\lambda_B \rightarrow \lambda_B - \nu\f$, so that increasing the coupling
in one box decreases it in the other. When a \f$\lambda\f$ crosses a boundary the
move becomes a molecule transfer between the boxes.

Steps:

- draw a random coupled change \f$\nu\f$ and compute
  \f$\lambda_{A,\mathbf{n}} = \lambda_A + \nu\f$, \f$\lambda_{B,\mathbf{n}} = \lambda_B - \nu\f$,
- if both stay within \f$[0,1]\f$ — **coupled \f$\lambda\f$-change**: rescale both
  fractional molecules simultaneously,
- if \f$\lambda_{A,\mathbf{n}} > 1\f$ — **transfer \f$B \rightarrow A\f$**: the fractional molecule of
  box \f$A\f$ becomes a whole molecule at its position; a new fractional molecule
  is inserted at a random position in box \f$A\f$ with the wrapped
  \f$\lambda_{A,\mathbf{n}} - 1\f$; the fractional molecule of box \f$B\f$ is removed, and a
  randomly selected whole molecule of box \f$B\f$ becomes the new fractional
  molecule with the wrapped \f$\lambda_{B,\mathbf{n}} + 1\f$,
- if \f$\lambda_{A,\mathbf{n}} < 0\f$ — **transfer \f$A \rightarrow B\f$**: the mirror image,
- accept or reject.

Acceptance rules (\f$\Delta\eta_A = \eta_A(\lambda_{A,\mathbf{n}}) - \eta_A(\lambda_{A,\mathbf{o}})\f$ and
likewise for \f$B\f$; \f$N_A\f$, \f$N_B\f$ are the total molecule counts of the component
including the fractional molecules):

\f[
\text{acc}_{\lambda} = \min\left(1,
  e^{-\beta\left[\Delta U_A + \Delta U_B\right]
     + \Delta\eta_A + \Delta\eta_B}\right)
\f]

\f[
\text{acc}_{B \rightarrow A} = \min\left(1,
  \frac{(N_B - 1)\, V_A}{N_A\, V_B}\,
  e^{-\beta\left[\Delta U_A + \Delta U_B\right]
     + \Delta\eta_A + \Delta\eta_B}\right),
\qquad
\text{acc}_{A \rightarrow B} = \min\left(1,
  \frac{(N_A - 1)\, V_B}{N_B\, V_A}\,
  e^{-\beta\left[\Delta U_A + \Delta U_B\right]
     + \Delta\eta_A + \Delta\eta_B}\right)
\f]

Useful for: Gibbs-ensemble VLE of dense systems while keeping the two boxes
decoupled in their \f$\lambda\f$ administration (a fractional molecule in each phase);
an alternative to the serial scheme, at the cost of slightly more bookkeeping.

Reference: W. Shi and E.J. Maginn, *Improvement in molecule exchange
efficiency in Gibbs ensemble Monte Carlo: Development and implementation of
the continuous fractional component move*, J. Comput. Chem. **29**, 2520-2530
(2008). <https://doi.org/10.1002/jcc.20977>

### Gibbs conventional CB/CFCMC (parallel) <a name="gibbs-conventional-cbcfcmc"></a>

Keyword: `"GibbsConventionalCBCFCMCProbability"` (component move)

The parallel CFCMC Gibbs move with CBMC growth/retracing of the fractional
molecules during the transfer sub-moves: the new fractional molecule of the
receiving box is grown with configurational bias (Rosenbluth weight
\f$W(\mathbf{n})\f$) and the removed fractional molecule of the donating box is
retraced (Rosenbluth weight \f$W(\mathbf{o})\f$). The steps are identical to
[Gibbs conventional CFCMC](#gibbs-conventional-cfcmc).

Acceptance rules:

\f[
\text{acc}_{\lambda} = \min\left(1,
  e^{-\beta\left[\Delta U_A + \Delta U_B\right]
     + \Delta\eta_A + \Delta\eta_B}\right)
\f]

\f[
\text{acc}_{B \rightarrow A} = \min\left(1,
  \frac{(N_B - 1)\, V_A}{N_A\, V_B}\,
  \frac{W(\mathbf{n})}{W(\mathbf{o})}\,
  e^{-\beta\left[\Delta U_A + \Delta U_B\right]
     + \Delta\eta_A + \Delta\eta_B}\right),
\qquad
\text{acc}_{A \rightarrow B} = \min\left(1,
  \frac{(N_A - 1)\, V_B}{N_B\, V_A}\,
  \frac{W(\mathbf{n})}{W(\mathbf{o})}\,
  e^{-\beta\left[\Delta U_A + \Delta U_B\right]
     + \Delta\eta_A + \Delta\eta_B}\right)
\f]

Useful for: parallel-scheme Gibbs VLE of flexible molecules.

References: W. Shi and E.J. Maginn, J. Comput. Chem. **29**, 2520-2530 (2008).
<https://doi.org/10.1002/jcc.20977>; A. Torres-Knoop et al., J. Chem. Theory
Comput. **10**, 942-952 (2014). <https://doi.org/10.1021/ct4009766>

### Gibbs identity change (CBMC) <a name="gibbs-identity-change-cbmc"></a>

Keyword: `"GibbsIdentityChangeProbability"` (component move), together with
`"GibbsIdentityChangesList"`

Simultaneously changes the identity of a molecule in one box and performs the
inverse change in the other box: a molecule of component \f$A\f$ in box I becomes
component \f$B\f$ while a molecule of component \f$B\f$ in box II becomes component
\f$A\f$. Both regrowths use CBMC. The total number of molecules of each component
is conserved, which makes the move ideal for closed multicomponent Gibbs
simulations.

Steps:

- select the pair of components from the identity-change list and the direction
  of the exchange,
- select a random molecule of component \f$A\f$ in box I and of component \f$B\f$
  in box II,
- regrow (CBMC) a molecule of component \f$B\f$ at the position of the first bead
  of the selected molecule in box I, and a molecule of component \f$A\f$ in
  box II,
- retrace both old molecules for their Rosenbluth weights,
- accept or reject the double transformation as a single move.

Acceptance rule:

\f[
\text{acc}(\mathbf{o} \rightarrow \mathbf{n}) = \min\left(1,
   \frac{W_\text{I}(\mathbf{n})\,W_\text{II}(\mathbf{n})}{W_\text{I}(\mathbf{o})\,W_\text{II}(\mathbf{o})}\,
   \frac{N_{A,\text{I}}\; N_{B,\text{II}}}{(N_{B,\text{I}}+1)(N_{A,\text{II}}+1)}\,
   e^{-\beta \Delta U_\text{non-CBMC}}\right)
\f]

Useful for: multicomponent VLE where compositions of the phases equilibrate
slowly; particularly effective for molecules of similar size (e.g. isomer
mixtures).

Reference: M.G. Martin and J.I. Siepmann, *Predicting multicomponent phase
equilibria and free energies of transfer for alkanes by molecular simulation*,
J. Am. Chem. Soc. **119**, 8921-8924 (1997).
<https://doi.org/10.1021/ja964218q>

----------------------------------------------------------------------------------

## Chemical-potential measurement moves <a name="chemical-potential-measurement-moves"></a>

### Widom insertion <a name="widom-insertion"></a>

Keyword: `"WidomProbability"` (component move)

A *measurement* move: a ghost molecule is grown at a random position with CBMC,
its Rosenbluth weight is recorded, and the molecule is always removed again.
The configuration of the system never changes. The excess chemical potential
follows from the average Rosenbluth weight; in adsorption studies this gives
direct access to the Henry coefficient and the heat of adsorption at infinite
dilution.

Steps:

- select a random box,
- select a random component,
- grow a ghost molecule at a random position with CBMC and record
  \f$W/W^\text{IG}\f$,
- remove the ghost molecule (the system is unchanged),
- accumulate the average.

There is no acceptance rule; the measured quantity is

\f[
\mu_\text{ex} = -k_B T \ln \left\langle \frac{W}{W^\text{IG}} \right\rangle
\f]

Useful for: chemical potentials of dilute species, Henry coefficients in
frameworks, consistency checks of open-ensemble simulations. Becomes
inefficient at high densities (use the CFCMC variants there).

Reference: B. Widom, *Some Topics in the Theory of Fluids*, J. Chem. Phys.
**39**, 2808-2812 (1963). <https://doi.org/10.1063/1.1734110>

### Widom CFCMC <a name="widom-cfcmc"></a>

Keyword: `"CFCMC_WidomProbability"` (component move)

Measures the chemical potential with the CFCMC machinery while keeping the
number of integer molecules fixed: the fractional molecule performs
\f$\lambda\f$-changes and reinsertions restricted so that no net molecule creation
occurs. The chemical potential follows from the ratio of the probabilities of
the \f$\lambda \rightarrow 1\f$ and \f$\lambda \rightarrow 0\f$ ends of the sampled \f$\lambda\f$ histogram,

\f[
\mu_\text{ex} = -k_B T \ln \frac{p(\lambda \rightarrow 1)}{p(\lambda \rightarrow 0)}
\f]

Steps:

- select a random box,
- select a random component,
- draw a random change \f$\Delta\lambda\f$ of the coupling parameter of the fractional
  molecule,
- if \f$\lambda + \Delta\lambda\f$ stays within \f$[0,1]\f$ — **\f$\lambda\f$-change**: rescale the
  fractional molecule,
- if \f$\lambda + \Delta\lambda\f$ crosses a boundary — **reinsertion**: the fractional
  molecule is moved to a random position with the wrapped \f$\lambda\f$ value; no
  integer molecule is created or destroyed,
- accept or reject and record the \f$\lambda\f$ histogram.

Acceptance rules (\f$\Delta\eta = \eta(\lambda_\mathbf{n}) - \eta(\lambda_\mathbf{o})\f$):

\f[
\text{acc}_{\lambda} = \min\left(1, e^{-\beta \Delta U + \Delta\eta}\right),
\qquad
\text{acc}_\text{reins} = \min\left(1, e^{-\beta \Delta U + \Delta\eta}\right)
\f]

Useful for: chemical potentials in dense phases (liquids, high loadings) where
plain Widom insertion fails.

Reference: W. Shi and E.J. Maginn, J. Chem. Theory Comput. **3**, 1451-1463
(2007). <https://doi.org/10.1021/ct7000039>

### Widom CB/CFCMC <a name="widom-cbcfcmc"></a>

Keyword: `"CFCMC_CBMC_WidomProbability"` (component move)

The CFCMC Widom measurement with configurational-bias growth of the fractional
molecule, for flexible molecules in dense phases. The steps are identical to
[Widom CFCMC](#widom-cfcmc), but in the reinsertion sub-move the fractional
molecule is grown with CBMC at the new position (Rosenbluth weight \f$W(\mathbf{n})\f$)
and retraced at the old position (Rosenbluth weight \f$W(\mathbf{o})\f$).

Acceptance rules:

\f[
\text{acc}_{\lambda} = \min\left(1, e^{-\beta \Delta U + \Delta\eta}\right),
\qquad
\text{acc}_\text{reins} = \min\left(1,
   \frac{W(\mathbf{n})}{W(\mathbf{o})}\,
   e^{-\beta \Delta U_\text{non-CBMC} + \Delta\eta}\right)
\f]

Reference: A. Torres-Knoop, S.P. Balaji, T.J.H. Vlugt, and D. Dubbeldam,
J. Chem. Theory Comput. **10**, 942-952 (2014).
<https://doi.org/10.1021/ct4009766>

----------------------------------------------------------------------------------

## Parallel tempering <a name="parallel-tempering"></a>

Keyword: `"ParallelTemperingSwapProbability"` (system move)

Attempts to exchange the complete configurations (all positions, box shapes,
and molecule counts) of two systems that are simulated in parallel at
different state points (temperature, pressure, or even force field). Replicas
trapped in low-temperature/high-density local minima are released by passing
through the more mobile replicas, dramatically improving ergodicity.

Steps:

- select a pair of systems (replicas),
- if the force fields of the two systems differ, recompute the energy of each
  configuration in the other system's Hamiltonian,
- evaluate the exchange acceptance factor (below),
- upon acceptance, swap the configurations, simulation boxes, molecule counts,
  and running energies of the two systems.

Acceptance rule (equal Hamiltonians, different temperatures):

\f[
\text{acc}(\mathbf{o} \rightarrow \mathbf{n}) = \min\left(1,
  e^{(\beta_B - \beta_A)\,(U_B - U_A)}\right)
\f]

For systems at different pressures the factor
\f$(p_B/p_A)^{\,N_B - N_A}\f$ multiplies the ratio (hyper-parallel tempering for
open systems), and for different Hamiltonians the full cross-energy form
\f$\exp\left[-\beta_A (U_A(\mathbf{x}_B) - U_A(\mathbf{x}_A)) - \beta_B (U_B(\mathbf{x}_A) - U_B(\mathbf{x}_B))\right]\f$
is used.

Useful for: systems with rugged free-energy landscapes: low-temperature
adsorption, strongly binding sites, phase transitions, and whenever a series of
state points (an isotherm or isobar) is computed anyway — the replicas then
improve each other's sampling at no extra cost.

References: G. Yan and J.J. de Pablo, *Hyper-parallel tempering Monte Carlo:
Application to the Lennard-Jones fluid and the restricted primitive model*,
J. Chem. Phys. **111**, 9509-9516 (1999). <https://doi.org/10.1063/1.480282>;
K. Hukushima and K. Nemoto, *Exchange Monte Carlo Method and Application to
Spin Glass Simulations*, J. Phys. Soc. Jpn. **65**, 1604-1608 (1996).
<https://doi.org/10.1143/JPSJ.65.1604>

----------------------------------------------------------------------------------

## Reaction ensemble moves <a name="reaction-ensemble-moves"></a>

Reaction-ensemble Monte Carlo (RxMC) samples chemical reaction equilibria
without simulating the reaction dynamics: a reaction step deletes the reactant
molecules and inserts the product molecules according to the stoichiometry,
with the ideal-gas partition functions of the species entering the acceptance
rule. Reactions are defined in the input via `"Reactions"` with stoichiometric
coefficients \f$\nu_i\f$ (negative for reactants, positive for products).

### Reaction (CBMC) <a name="reaction-cbmc"></a>

Keyword: `"ReactionCBMCProbability"` (system move)

Performs one reaction step in the forward or backward direction using CBMC to
insert the products and retrace the reactants.

Steps:

- select a random box,
- select a random reaction from the reaction list and a random direction
  (forward/backward),
- select random reactant molecules according to the stoichiometry and retrace
  them with CBMC,
- grow the product molecules at random positions with CBMC,
- accept or reject the combined deletion/insertion.

Acceptance rule (forward direction, ideal-gas partition functions per unit
volume \f$q_i\f$):

\f[
\text{acc}(\mathbf{o} \rightarrow \mathbf{n}) = \min\left(1,
   \prod_i \left(q_i V\right)^{\nu_i}\,
   \prod_i \frac{N_i!}{(N_i + \nu_i)!}\;
   \frac{\prod_\text{products} W(\mathbf{n})/W^\text{IG}}{\prod_\text{reactants} W(\mathbf{o})/W^\text{IG}}\;
   e^{-\beta \Delta U_\text{non-CBMC}}\right)
\f]

Useful for: chemical equilibria (e.g. ammonia synthesis, esterification,
propene metathesis) in bulk phases and inside nanoporous catalysts, including
the shift of reaction equilibria by confinement.

References: W.R. Smith and B. Triska, *The reaction ensemble method for the
computer simulation of chemical and phase equilibria*, J. Chem. Phys. **100**,
3019-3027 (1994). <https://doi.org/10.1063/1.466443>;
J.K. Johnson, A.Z. Panagiotopoulos, and K.E. Gubbins, *Reactive canonical Monte
Carlo: a new simulation technique for reacting or associating fluids*, Mol.
Phys. **81**, 717-733 (1994). <https://doi.org/10.1080/00268979400100481>

### Reaction conventional CFCMC (parallel) <a name="reaction-conventional-cfcmc"></a>

Keyword: `"ReactionConventionalCFCMCProbability"` (system move)

Parallel Rx/CFC version of the reaction step: the system permanently carries
fractional molecules of *both* the reactants and the products of the reaction,
coupled through a single parameter \f$\lambda\f$: the products are scaled with
\f$\lambda\f$ and the reactants with \f$1 - \lambda\f$, so the reaction proceeds gradually as
\f$\lambda\f$ moves from 0 to 1. A Wang-Landau bias \f$\eta(\lambda)\f$ flattens the reaction
path. When \f$\lambda\f$ crosses a boundary, an actual chemical reaction is performed
with the wrapped \f$\lambda\f$.

Steps:

- select a random box,
- select a random reaction from the reaction list,
- draw \f$\lambda_\mathbf{n} = \lambda_\mathbf{o} + \Delta\lambda\f$,
- if \f$\lambda_\mathbf{n}\f$ stays within \f$[0,1]\f$ — **\f$\lambda\f$-change**: rescale all fractional
  molecules of the reaction (reactants to \f$1-\lambda_\mathbf{n}\f$, products to \f$\lambda_\mathbf{n}\f$)
  and compute \f$\Delta U\f$,
- if \f$\lambda_\mathbf{n} > 1\f$ — **forward reaction** at the wrapped
  \f$\lambda' = \lambda_\mathbf{n} - 1\f$: the (nearly fully coupled) product fractional molecules
  become whole molecules in place; randomly selected whole reactant molecules
  become the new reactant fractional molecules at \f$\lambda'\f$; the nearly decoupled
  reactant fractional molecules are deleted (retraced, Rosenbluth weight
  \f$W(\mathbf{o})\f$); new nearly decoupled product fractional molecules are grown
  (Rosenbluth weight \f$W(\mathbf{n})\f$),
- if \f$\lambda_\mathbf{n} < 0\f$ — **backward reaction**: the mirror image,
- accept or reject.

Acceptance rules (\f$q_i\f$ the ideal-gas partition function of component \f$i\f$,
\f$\nu_i\f$ signed stoichiometric coefficients of the attempted direction, \f$N_i\f$
whole-molecule counts, \f$\Delta\eta = \eta(\lambda_\mathbf{n}) - \eta(\lambda_\mathbf{o})\f$):

\f[
\text{acc}_{\lambda} = \min\left(1, e^{-\beta \Delta U + \Delta\eta}\right)
\f]

\f[
\text{acc}_\text{rxn} = \min\left(1,
   \frac{W(\mathbf{n})/W^\text{IG}(\mathbf{n})}{W(\mathbf{o})/W^\text{IG}(\mathbf{o})}\,
   \prod_i \left(q_i V\right)^{\nu_i}
   \prod_i \frac{N_i!}{(N_i + \nu_i)!}\;
   e^{-\beta \Delta U + \Delta\eta}\right)
\f]

Useful for: reaction equilibria in dense phases (liquid-phase reactions,
reactions at high loading in nanoporous materials) where direct
deletion/insertion of whole molecules fails.

References: T.W. Rosch and E.J. Maginn, *Reaction Ensemble Monte Carlo
Simulation of Complex Molecular Systems*, J. Chem. Theory Comput. **7**,
269-279 (2011). <https://doi.org/10.1021/ct100615j>;
A. Poursaeidesfahani, R. Hens, A. Rahbari, M. Ramdin, D. Dubbeldam, and
T.J.H. Vlugt, *Efficient Application of Continuous Fractional Component Monte
Carlo in the Reaction Ensemble*, J. Chem. Theory Comput. **13**, 4452-4466
(2017). <https://doi.org/10.1021/acs.jctc.7b00092>

### Reaction conventional CB/CFCMC (parallel) <a name="reaction-conventional-cbcfcmc"></a>

Keyword: `"ReactionConventionalCBCFCMCProbability"` (system move)

The parallel Rx/CFC reaction move with configurational-bias growth of the
fractional molecules: in the boundary-crossing (reaction) branch the new
nearly decoupled fractional molecules are grown with CBMC using multiple trial
positions and orientations, and the deleted fractional molecules are retraced
with CBMC. The steps and acceptance rules are those of
[Reaction conventional CFCMC](#reaction-conventional-cfcmc):

\f[
\text{acc}_{\lambda} = \min\left(1, e^{-\beta \Delta U + \Delta\eta}\right)
\f]

\f[
\text{acc}_\text{rxn} = \min\left(1,
   \frac{W(\mathbf{n})/W^\text{IG}(\mathbf{n})}{W(\mathbf{o})/W^\text{IG}(\mathbf{o})}\,
   \prod_i \left(q_i V\right)^{\nu_i}
   \prod_i \frac{N_i!}{(N_i + \nu_i)!}\;
   e^{-\beta \Delta U + \Delta\eta}\right)
\f]

with \f$W(\mathbf{n})\f$ and \f$W(\mathbf{o})\f$ the CBMC Rosenbluth weights of the grown and
retraced fractional molecules.

Useful for: reaction equilibria of flexible molecules in dense phases.

References: T.W. Rosch and E.J. Maginn, J. Chem. Theory Comput. **7**, 269-279
(2011). <https://doi.org/10.1021/ct100615j>; A. Torres-Knoop et al., J. Chem.
Theory Comput. **10**, 942-952 (2014). <https://doi.org/10.1021/ct4009766>

### Reaction CFCMC (serial) <a name="reaction-cfcmc"></a>

Keyword: `"ReactionCFCMCProbability"` (system move)

Serial Rx/CFC version of the reaction move: only *one side* of the reaction is
fractional at a time (indicator \f$\delta\f$: reactant side or product side), coupled
through a single \f$\lambda \in [0,1]\f$ with separate bias functions
\f$\eta_\text{reac}(\lambda)\f$ and \f$\eta_\text{prod}(\lambda)\f$ for the two sides. Three
sub-moves change \f$\lambda\f$, flip the fractional side at fixed \f$\lambda\f$, or perform
the reaction on whole molecules.

Steps:

- select a random box,
- select a random reaction from the reaction list,
- select a sub-move: with 50% probability a **\f$\lambda\f$-change**; otherwise, below
  the \f$\lambda\f$ switch point a **fractional reaction**, above it a
  **whole-molecule reaction**,
- **\f$\lambda\f$-change**: draw \f$\lambda_\mathbf{n} = \lambda_\mathbf{o} + \Delta\lambda\f$ (rejected if outside
  \f$[0,1]\f$) and rescale the fractional molecules of the active side,
- **fractional reaction**: delete the fractional molecules of the active side
  (retraced, Rosenbluth weight \f$W(\mathbf{o})\f$) and grow fractional molecules of the
  other side at the *same* \f$\lambda\f$ (Rosenbluth weight \f$W(\mathbf{n})\f$); \f$\delta\f$ flips,
  the whole-molecule counts do not change,
- **whole-molecule reaction**: the fractional molecules of the active side
  become whole molecules in place, and randomly selected whole molecules of
  the other side become the new fractional molecules at the same \f$\lambda\f$;
  \f$\delta\f$ flips, all positions and \f$\lambda\f$ stay the same,
- accept or reject.

Acceptance rules (\f$\nu_i^\text{lost}\f$/\f$\nu_i^\text{gained}\f$ the stoichiometric
coefficients of the disappearing/appearing fractional side, \f$\bar\nu_i\f$ the
signed change in the whole-molecule count of component \f$i\f$, and
\f$\Delta\eta\f$ the difference of the side bias functions evaluated at the current
\f$\lambda\f$):

\f[
\text{acc}_{\lambda} = \min\left(1,
   e^{-\beta \Delta U + \eta_\delta(\lambda_\mathbf{n}) - \eta_\delta(\lambda_\mathbf{o})}\right)
\f]

\f[
\text{acc}_\text{frac} = \min\left(1,
   \frac{W(\mathbf{n})/W^\text{IG}(\mathbf{n})}{W(\mathbf{o})/W^\text{IG}(\mathbf{o})}\,
   \prod_i \left(q_i V\right)^{\nu_i^\text{gained} - \nu_i^\text{lost}}\,
   e^{-\beta \Delta U + \Delta\eta}\right)
\f]

\f[
\text{acc}_\text{whole} = \min\left(1,
   \prod_i \frac{N_i!}{(N_i + \bar\nu_i)!}\;
   e^{-\beta \Delta U + \Delta\eta}\right)
\f]

Note that the fractional reaction carries the \f$(q_i V)^{\nu_i}\f$ factors but no
factorials (whole-molecule counts unchanged), while the whole-molecule
reaction carries the factorials but no \f$(q_i V)^{\nu_i}\f$ factors (the total
number of molecules per component, whole plus fractional, is unchanged).

Useful for: the same problems as the parallel scheme; the serial scheme keeps
only one fractional set at a time, which simplifies the simulation and
improves efficiency for many systems.

References: T.W. Rosch and E.J. Maginn, J. Chem. Theory Comput. **7**, 269-279
(2011). <https://doi.org/10.1021/ct100615j>;
A. Poursaeidesfahani, R. Hens, A. Rahbari, M. Ramdin, D. Dubbeldam, and
T.J.H. Vlugt, J. Chem. Theory Comput. **13**, 4452-4466 (2017).
<https://doi.org/10.1021/acs.jctc.7b00092>

### Reaction CB/CFCMC (serial) <a name="reaction-cbcfcmc"></a>

Keyword: `"ReactionCBCFCMCProbability"` (system move)

The serial Rx/CFC reaction move with configurational-bias growth: in the
fractional-reaction sub-move the new fractional molecules are grown with CBMC
(multiple trial positions and orientations) and the removed ones are retraced
with CBMC. The steps and acceptance rules are those of
[Reaction CFCMC (serial)](#reaction-cfcmc):

\f[
\text{acc}_{\lambda} = \min\left(1,
   e^{-\beta \Delta U + \eta_\delta(\lambda_\mathbf{n}) - \eta_\delta(\lambda_\mathbf{o})}\right)
\f]

\f[
\text{acc}_\text{frac} = \min\left(1,
   \frac{W(\mathbf{n})/W^\text{IG}(\mathbf{n})}{W(\mathbf{o})/W^\text{IG}(\mathbf{o})}\,
   \prod_i \left(q_i V\right)^{\nu_i^\text{gained} - \nu_i^\text{lost}}\,
   e^{-\beta \Delta U + \Delta\eta}\right)
\f]

\f[
\text{acc}_\text{whole} = \min\left(1,
   \prod_i \frac{N_i!}{(N_i + \bar\nu_i)!}\;
   e^{-\beta \Delta U + \Delta\eta}\right)
\f]

with \f$W(\mathbf{n})\f$ and \f$W(\mathbf{o})\f$ the CBMC Rosenbluth weights.

Useful for: reaction equilibria of flexible molecules in dense phases with the
serial (single fractional set) bookkeeping.

References: T.W. Rosch and E.J. Maginn, J. Chem. Theory Comput. **7**, 269-279
(2011). <https://doi.org/10.1021/ct100615j>; A. Torres-Knoop et al., J. Chem.
Theory Comput. **10**, 942-952 (2014). <https://doi.org/10.1021/ct4009766>

----------------------------------------------------------------------------------

For the derivations of these acceptance rules and implementation details, see:

> D. Dubbeldam, A. Torres-Knoop, and K.S. Walton, *On the Inner Workings of
> Monte Carlo Codes*, Mol. Simul. **39**, 1253-1292 (2013).
> <http://dx.doi.org/10.1080/08927022.2013.819102>
