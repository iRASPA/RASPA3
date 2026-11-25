# Properties
\page properties Properties

This page lists everything `RASPA3` measures during a simulation: what the
quantity is, how it is defined, how it is sampled, and where the result is
written. Some properties are always computed; most are switched on with a key
in the input. The keys themselves, with their defaults, are documented on the
[Commands](\ref commands) page; the methods that produce whole isotherms and
free-energy profiles are on the [Drivers](\ref drivers) pages.

## Table of Contents

<!-- TOC -->
* [How properties are sampled](#how-properties-are-sampled)
* [Overview](#property-overview)
* [Thermodynamic averages](#thermodynamic-averages)
  * [Energies](#property-energies)
  * [Temperature](#property-temperature)
  * [Pressure](#property-pressure)
  * [Simulation box, volume and density](#property-box)
  * [Loading](#property-loading)
  * [Enthalpy of adsorption](#property-enthalpy)
  * [Partial molar properties](#property-partial-molar)
* [Chemical potential and free energy](#chemical-potential-and-free-energy)
  * [Rosenbluth weight, Henry coefficient and chemical potential](#property-widom)
  * [The λ histogram](#property-lambda-histogram)
  * [The derivative ∂U/∂λ](#property-dudlambda)
  * [The macrostate distribution](#property-lnpi)
  * [Reweighted isotherms, coexistence and BET area](#property-isotherms)
  * [Free-energy profile and transmission coefficient](#property-dctst)
* [Distributions](#distributions)
  * [Histogram of the energy](#property-energy-histogram)
  * [Histogram of the number of molecules](#property-number-histogram)
  * [Density grids](#property-density-grid)
  * [Radial distribution function, conventional](#property-conventional-rdf)
  * [Radial distribution function, force-based](#property-rdf)
* [Evolution during the run](#evolution-during-the-run)
* [Dynamics](#dynamics)
  * [Mean-squared displacement](#property-msd)
  * [Velocity autocorrelation function](#property-vacf)
* [Conformation of molecules](#conformation-of-molecules)
  * [Bond, bend and torsion histograms](#property-molecule-properties)
  * [Shape: the gyration tensor](#property-molecule-shape)
  * [Backbone: chain statistics](#property-molecule-backbone)
  * [End-to-end autocorrelation and relaxation time](#property-end-to-end-acf)
* [Mechanical properties](#mechanical-properties)
  * [Elastic constants from fluctuations](#property-elastic-fluctuation)
  * [Elastic constants, normal modes and phonons at a minimum](#property-minimization)
* [Reaction kinetics](#reaction-kinetics)
  * [Radical waiting times](#property-radical-waiting-times)
* [Configurations](#configurations)
<!-- TOC -->

----------------------------------------------------------------------------------

## How properties are sampled {#how-properties-are-sampled}

**When.** Properties are measured during the production stage only. A
property that is switched on has a *sampling* interval (`"Sample...Every"`)
and a *writing* interval (`"Write...Every"`), both in cycles; a cycle is a
sweep of Monte Carlo moves or one time step of molecular dynamics. The output
files are rewritten at the writing interval, so the convergence of a property
can be followed while the simulation runs.

**Errors.** The production run is divided into `"NumberOfBlocks"` contiguous
blocks. An average is the mean of its block averages \f$A_b\f$, and its error
is the 95% confidence interval of that mean,

\f[
\delta A = t_{0.975,\,n-1}\,\frac{s}{\sqrt{n}}, \qquad
s^2 = \frac{1}{n-1}\sum_{b=1}^{n}(A_b - \bar A)^2 ,
\f]

with \f$n\f$ the number of blocks and \f$t\f$ the Student factor. The blocks
must be long compared to the correlation time of the property for this to be
a valid error. Quantities that are *functions* of averages (a ratio of
fluctuations, for instance) are formed per block from the block sums and then
averaged, so that they have errors of the same kind.

**Biased sampling.** With a fractional molecule the simulation performs a
biased walk over the coupling parameter \f$\lambda\f$. Every sample then
carries the weight that removes the bias,
\f$w = e^{-\eta(\lambda)}\f$, and an average is
\f$\langle A\rangle = \sum_s w_s A_s/\sum_s w_s\f$, so that the results are
those of the unbiased ensemble. The density grids and the radial distribution
functions, among others, are accumulated this way.

**Which molecules count.** A fractional molecule is not a molecule of the
system. The loadings count the whole (integer) molecules.

**Units.** Energies are reported in kelvin (\f$U/k_B\f$) unless stated
otherwise, lengths in ångström, times in picoseconds, pressures in pascal.

## Overview {#property-overview}

| Property | Switched on with | Engine | Output |
|---|---|---|---|
| [Energies](#property-energies) | always | MC, MD | output file |
| [Temperature](#property-temperature) | always | MD | output file |
| [Pressure](#property-pressure) | always | MC, MD | output file |
| [Box, volume, density](#property-box) | always | MC, MD | output file |
| [Loading](#property-loading) | always | MC, MD | output file |
| [Enthalpy of adsorption](#property-enthalpy) | always (open ensembles) | MC | output file |
| [Partial molar properties](#property-partial-molar) | always (open ensembles) | MC | output file |
| [Rosenbluth weight, Henry coefficient](#property-widom) | `"WidomProbability"` | MC | output file |
| [λ histogram](#property-lambda-histogram) | a CFCMC move | MC | output file |
| [∂U/∂λ](#property-dudlambda) | `"ThermodynamicIntegration"` | MC | output file |
| [Macrostate distribution](#property-lnpi) | TMMC drivers | MC | `tmmc`, `output` |
| [Energy histogram](#property-energy-histogram) | `"ComputeEnergyHistogram"` | MC, MD | `energy_histogram` |
| [Number-of-molecules histogram](#property-number-histogram) | `"ComputeNumberOfMoleculesHistogram"` | MC | `number_of_molecules_histogram` |
| [Density grid](#property-density-grid) | `"ComputeDensityGrid"` | MC, MD | `density_grids` |
| [RDF, conventional](#property-conventional-rdf) | `"ComputeConventionalRDF"` | MC, MD | `conventional_rdf` |
| [RDF, force-based](#property-rdf) | `"ComputeRDF"` | MC, MD | `rdf` |
| [Evolution of N and V](#evolution-during-the-run) | `"ComputeNumberOfMoleculesEvolution"`, `"ComputeVolumeEvolution"` | MC, MD | `number_of_molecules_evolution`, `volume_evolution` |
| [Mean-squared displacement](#property-msd) | `"ComputeMSD"` | MD | `msd` |
| [Velocity autocorrelation](#property-vacf) | `"ComputeVACF"` | MD | `vacf` |
| [Bond, bend, torsion histograms](#property-molecule-properties) | `"ComputeMoleculeProperties"` | MC, MD | `molecule_properties` |
| [Molecule shape](#property-molecule-shape) | `"ComputeMoleculeShape"` | MC, MD | `molecule_shape` |
| [Backbone statistics](#property-molecule-backbone) | `"ComputeMoleculeBackbone"` | MC, MD | `molecule_backbone` |
| [End-to-end autocorrelation](#property-end-to-end-acf) | `"ComputeEndToEndACF"` | MC, MD | `end_to_end_acf` |
| [Elastic constants (fluctuations)](#property-elastic-fluctuation) | `"ComputeElasticConstantsFromFluctuations"` | MD | output file |
| [Elastic constants, modes, phonons](#property-minimization) | `"ComputeElasticConstants"`, `"ComputeNormalModes"`, `"ComputePhononDispersion"` | minimization | output file |
| [Radical waiting times](#property-radical-waiting-times) | `"ComputeRadicalWaitingTimes"` | polymerization | `radical_waiting_times` |
| [Movies and data files](#configurations) | `"OutputPDBMovie"`, `"WriteLammpsData"` | MC, MD | `movies`, `lammps` |

The first group of keys belongs in the system block, the conformation
analyses in the component block.

----------------------------------------------------------------------------------

## Thermodynamic averages {#thermodynamic-averages}

These are computed in every run and reported in the output file, with their
errors, at the end of the production stage.

### Energies {#property-energies}

The potential energy is averaged in total and per contribution: the
interaction of the molecules with the framework and with each other, each
split into van der Waals, real-space Coulomb and reciprocal-space (Ewald)
parts; the tail correction; the polarization energy when it is switched on;
the external field; and the intramolecular energy (bonds, bends, torsions,
and the intramolecular non-bonded terms). The interactions between molecules
are also given per pair of components.

In molecular dynamics the kinetic energy (translational and rotational), the
energy of the thermostat and barostat, and the conserved quantity of the
integrator are reported as well.

*Check.* At the end of a Monte Carlo run the running energies, which were
updated with the energy difference of every accepted move, are compared with
the energies computed from scratch. In molecular dynamics the drift of the
conserved quantity plays this role.

### Temperature {#property-temperature}

In molecular dynamics the kinetic temperature follows from equipartition,

\f[
T = \frac{2\,\langle K\rangle}{N_f\,k_B} ,
\f]

with \f$K\f$ the kinetic energy and \f$N_f\f$ the number of degrees of
freedom (constraints are subtracted). It is
reported in total and separately for the translational degrees of freedom
(centres of mass and flexible atoms) and the rotational degrees of freedom of
rigid molecules. The separate temperatures must agree; when they do not, the
time step is too large or the equilibration too short.

### Pressure {#property-pressure}

The pressure is the sum of an ideal-gas and an excess part,

\f[
p = \rho\,k_B T + p^\text{ex}, \qquad
p^\text{ex} = -\left\langle \frac{\partial U}{\partial V}\right\rangle
= \frac{1}{3V}\left\langle \sum_{i<j} \mathbf{r}_{ij}\cdot\mathbf{F}_{ij}\right\rangle ,
\f]

where the sum is the virial. The full pressure *tensor*,
\f$p_{\alpha\beta} = \rho k_BT\,\delta_{\alpha\beta}
+ \frac{1}{V}\langle\sum_{i<j} r_{ij,\alpha}F_{ij,\beta}\rangle\f$, is
reported as well; its diagonal elements are equal in an isotropic fluid.

The reported pressure is the *molecular* pressure: \f$\rho\f$ is the number
density of molecules and the virial is that of the forces between the centres
of mass. The atomic form (density of atoms, virial of all forces including
the bonded ones) has the same average but different fluctuations. The
molecular form is the one a Monte Carlo volume move samples. The long-range
(Ewald) part and the tail correction are included.

*Note.* For an adsorbate inside a rigid framework this is the pressure of the
adsorbed phase in a fixed external field; it is not the pressure of the gas
reservoir, which is the input of the simulation.

### Simulation box, volume and density {#property-box}

When the box can change (a Monte Carlo volume move, a barostat, the Gibbs
ensemble) the averages of the box lengths, the angles and the volume are
reported, with the number density and the mass density that follow from them.
In the isothermal-isobaric ensemble the average density at the imposed
pressure is the equation of state of the fluid.

### Loading {#property-loading}

The *absolute* loading is the average number of molecules of each component
in the system, \f$\langle N_i\rangle\f$. It is reported in molecules per
simulation cell, molecules per unit cell, mol per kilogram of framework and
milligram per gram of framework.

Experiments measure the *excess* loading: the amount adsorbed beyond what
would be present in the pore volume at the density of the bulk gas,

\f[
n^\text{ex}_i = n^\text{abs}_i - \rho^\text{gas}_i\,V_\text{pore},
\qquad V_\text{pore} = \xi\,V ,
\f]

with \f$\xi\f$ the helium void fraction of the framework
(`"HeliumVoidFraction"`, from a separate calculation) and
\f$\rho^\text{gas}\f$ the density of the gas phase at the pressure and
temperature of the simulation, from the Peng-Robinson equation of state. At
low pressure the two loadings are nearly equal; at high pressure the excess
loading goes through a maximum and decreases, while the absolute loading
keeps rising.

A series of loadings at increasing pressure is the adsorption isotherm.

### Enthalpy of adsorption {#property-enthalpy}

The enthalpy of adsorption (minus the isosteric heat) is the enthalpy change
when a molecule goes from the gas phase to the adsorbed phase. In a
grand-canonical simulation it follows from the fluctuations of the energy and
the number of molecules, without a second simulation:

\f[
\Delta H = \frac{\langle U N\rangle - \langle U\rangle\langle N\rangle}
{\langle N^2\rangle - \langle N\rangle^2} - \langle U^\text{g}\rangle - k_B T .
\f]

The first term is \f$(\partial U/\partial N)_{V,T}\f$, the energy gained per
added molecule. \f$\langle U^\text{g}\rangle\f$ is the average intramolecular
energy of an isolated molecule in the ideal gas, which is zero for a rigid
molecule; without it the enthalpy of a flexible molecule would contain its
internal energy. The term \f$k_BT\f$ is the \f$pV\f$ work of the ideal gas.

For a mixture the derivative with respect to the number of molecules of
component \f$i\f$ at constant numbers of the others requires the full
covariance matrix of the particle numbers,

\f[
\Delta H_i = \sum_j \left[\langle U N_j\rangle - \langle U\rangle\langle N_j\rangle\right]
\left[\mathrm{Cov}(N)^{-1}\right]_{ji} - \langle U^\text{g}_i\rangle - k_B T ,
\qquad
\mathrm{Cov}(N)_{ij} = \langle N_i N_j\rangle - \langle N_i\rangle\langle N_j\rangle .
\f]

The total enthalpy of adsorption is the average of the \f$\Delta H_i\f$
weighted with the mole fractions of the adsorbed phase.

*Notes.* The formula needs fluctuations in the number of molecules, so it
exists in open ensembles only; in a closed ensemble the covariance matrix is
singular and no value is reported. It converges slowly at high loading, where
the fluctuations of \f$N\f$ become small. At zero loading the Widom route
below is the better one.

### Partial molar properties {#property-partial-molar}

The same fluctuation approach gives the partial molar internal energy and the
partial molar volume of every component that is exchanged with a reservoir,

\f[
\bar U_i = \left(\frac{\partial U}{\partial N_i}\right)_{T,V,N_{j\ne i}}
= \sum_j \left[\langle N_j U\rangle - \langle N_j\rangle\langle U\rangle\right]\left[\mathrm{Cov}(N)^{-1}\right]_{ji},
\f]

\f[
\bar V_i = \left(\frac{\partial V}{\partial N_i}\right)_{T,p,N_{j\ne i}}
= \sum_j \left[\langle N_j V\rangle - \langle N_j\rangle\langle V\rangle\right]\left[\mathrm{Cov}(N)^{-1}\right]_{ji} .
\f]

The partial molar volume requires an ensemble in which both the volume and
the numbers of molecules fluctuate.

----------------------------------------------------------------------------------

## Chemical potential and free energy {#chemical-potential-and-free-energy}

### Rosenbluth weight, Henry coefficient and chemical potential {#property-widom}

Switched on with `"WidomProbability"` (a component move).

The Widom move inserts a ghost molecule with configurational bias, records
the Rosenbluth weight \f$W\f$ of the growth, and removes the molecule again.
The average weight, relative to that of an isolated molecule in the ideal gas
\f$\langle W^\text{IG}\rangle\f$ (`"IdealGasRosenbluthWeight"`, which is one
for a rigid molecule), gives the excess chemical potential,

\f[
\mu^\text{ex} = -k_B T\,\ln\frac{\langle W\rangle}{\langle W^\text{IG}\rangle} .
\f]

In an *empty* framework the same ratio is the Henry coefficient, the slope of
the isotherm at zero pressure,

\f[
K_H = \frac{1}{R\,T\,\rho_f}\,\frac{\langle W\rangle}{\langle W^\text{IG}\rangle} ,
\f]

with \f$\rho_f\f$ the density of the framework; \f$n = K_H\,p\f$ at low
pressure. Ratios of Henry coefficients are the adsorption selectivities at
infinite dilution.

*Notes.* The method fails at high density, where almost every insertion
overlaps and the average is dominated by rare events. For flexible molecules
the ideal-gas Rosenbluth weight must be computed first, at the same
temperature. In the Gibbs ensemble the corresponding quantity is measured
during the transfer moves.

### The \f$\lambda\f$ histogram {#property-lambda-histogram}

Available with every continuous-fractional-component move.

The fractional molecule performs a random walk over its coupling parameter
\f$\lambda\f$, flattened by the bias \f$\eta(\lambda)\f$ that is built during
equilibration. The observed histogram, corrected for the bias, is the
probability distribution \f$p(\lambda)\f$ of the unbiased system, and its two
ends give the excess chemical potential,

\f[
\mu^\text{ex} = -k_B T\,\ln\frac{p(\lambda = 1)}{p(\lambda = 0)} ,
\f]

since \f$\lambda = 0\f$ is a molecule that does not interact and
\f$\lambda = 1\f$ a molecule that interacts fully. The output lists the
histogram, the bias and the resulting free energy per \f$\lambda\f$ bin. A
flat observed histogram means the bias is good; a histogram with unvisited
bins means the equilibration was too short.

### The derivative \f$\partial U/\partial\lambda\f$ {#property-dudlambda}

With `"ThermodynamicIntegration"` on a component, the derivative of the
energy with respect to the coupling parameter is accumulated per
\f$\lambda\f$ bin. Its average is the integrand of thermodynamic integration,

\f[
\mu^\text{ex} = \int_0^1 \left\langle\frac{\partial U}{\partial\lambda}\right\rangle_\lambda d\lambda .
\f]

See the pages on [thermodynamic integration](\ref driver_thermodynamic_integration)
and its [parallel version](\ref driver_parallel_thermodynamic_integration),
which also perform the integration.

### The macrostate distribution {#property-lnpi}

The transition-matrix drivers record the collection matrix of acceptance
probabilities and derive from it \f$\ln\Pi(N)\f$, the logarithm of the
probability of every number of molecules. It is written with its error and
the histogram of visits. See
[transition-matrix Monte Carlo](\ref driver_transition_matrix) and
[parallel TMMC](\ref driver_parallel_tmmc).

### Reweighted isotherms, coexistence and BET area {#property-isotherms}

The [reweighted-histogram](\ref driver_reweighted_histogram) and
[parallel TMMC](\ref driver_parallel_tmmc) drivers turn their distributions
into derived properties:

- **isotherms** on a fine grid of pressures, at one or several temperatures,
  with the adsorption and desorption branches where the distribution is
  bimodal;
- **vapour-liquid coexistence** of a fluid in a box: the saturation pressure
  and the saturated vapour and liquid densities, from the equal-area
  criterion;
- the **BET surface area** (`"ComputeBET"`): the BET equation fitted to the
  isotherm over the pressure range selected by the Rouquerol consistency
  criteria, which gives the monolayer capacity \f$n_m\f$ and the constant
  \f$C\f$; the area is \f$n_m\f$ times the cross-section of the adsorbate
  molecule.

The [NLDFT driver](\ref driver_nldft) reports an isotherm and a BET area from
density functional theory.

### Free-energy profile and transmission coefficient {#property-dctst}

The dcTST drivers measure the free-energy profile \f$\beta F(q)\f$ of a
molecule along a hopping coordinate, the probability at the barrier, the
transmission coefficients and the resulting hopping rates and diffusion
coefficients. See the pages on the
[free-energy profile](\ref driver_dctst_free_energy) and the
[transmission coefficient](\ref driver_dctst).

----------------------------------------------------------------------------------

## Distributions {#distributions}

### Histogram of the energy {#property-energy-histogram}

`"ComputeEnergyHistogram"` (system). Output directory: `energy_histogram`.

The probability distribution of the potential energy of the system, in total
and for its van der Waals, Coulomb and polarization parts, between the given
lower and upper limits. In the canonical ensemble its width is the heat
capacity,

\f[
C_V^\text{ex} = \frac{\langle U^2\rangle - \langle U\rangle^2}{k_B T^2} ,
\f]

and two peaks signal coexisting states. In adsorption the distribution shows
the adsorption sites: molecules at strong sites and at weak sites give
separate contributions.

### Histogram of the number of molecules {#property-number-histogram}

`"ComputeNumberOfMoleculesHistogram"` (system). Output directory:
`number_of_molecules_histogram`.

The probability \f$P(N)\f$ of finding \f$N\f$ molecules of each component, in
ensembles where the number fluctuates. Its mean is the loading and its
variance is related to the slope of the isotherm,

\f[
\langle N^2\rangle - \langle N\rangle^2 = k_B T\left(\frac{\partial\langle N\rangle}{\partial\mu}\right)_{V,T} ,
\f]

so a steep isotherm has a wide distribution. A distribution with two peaks
means the system switches between two states, as at capillary condensation;
the averages of such a run depend on how often it switched.

### Density grids {#property-density-grid}

`"ComputeDensityGrid"` (system). Output directory: `density_grids`.

The simulation cell is divided into voxels and the positions of the atoms are
accumulated, giving the three-dimensional average density
\f$\rho(\mathbf{r})\f$ of the adsorbate. It shows *where* molecules adsorb.

- **Normalization.** `"Max"` scales the grid to its maximum, for
  visualization. `"NumberDensity"` gives the number density, whose integral
  over the cell is the loading.
- **Binning.** With `"Standard"` binning an atom contributes to the voxel it
  is in. With `"Equitable"` binning it is spread over the neighbouring voxels
  according to its position, which gives a smoother grid at the same number
  of samples, particularly on fine grids.
- **Per atom type.** With a list of pseudo-atoms a separate grid is made for
  each, which resolves the orientation of the molecules: the carbon and the
  oxygen atoms of CO<sub>2</sub>, for instance.

The grids are written as Gaussian cube files, one per component (and per
listed pseudo-atom), which visualization programs show as isosurfaces
together with the framework.

### Radial distribution function, conventional {#property-conventional-rdf}

`"ComputeConventionalRDF"` (system). Output directory: `conventional_rdf`.

The radial distribution function \f$g_{AB}(r)\f$ is the density of atoms of
type \f$B\f$ at a distance \f$r\f$ from an atom of type \f$A\f$, relative to
the average density. It is measured by counting pairs in shells,

\f[
g_{AB}(r) = \frac{V}{N_A N_B}\;
\frac{\left\langle \sum_{i\in A}\sum_{j\in B} \delta(r - r_{ij})\right\rangle}{4\pi r^2} ,
\f]

with the delta function replaced by a bin of width \f$\Delta r\f$. It tends
to one at large distance in a homogeneous fluid. Its first peak gives the
nearest-neighbour distance and its integral the coordination number,

\f[
n_{AB}(r) = 4\pi\rho_B \int_0^{r} g_{AB}(r')\,r'^2\,dr' .
\f]

The range cannot exceed half the smallest perpendicular width of the box.
The statistical noise of a bin grows as the bin is made narrower, because
fewer pairs fall into it.

### Radial distribution function, force-based {#property-rdf}

`"ComputeRDF"` (system). Output directory: `rdf`.

The force-based estimator of Borgis *et al.* removes the trade-off between
resolution and noise of the histogram method. It uses the identity that the
derivative of \f$g(r)\f$ is determined by the mean force between a pair at
that distance. Integrating that identity gives \f$g(r)\f$ as an average over
*all* pairs beyond (or within) the distance \f$r\f$, each weighted with the
relative force along the line connecting the two atoms,

\f[
g_{AB}(r) \;\propto\; \beta\left\langle \sum_{i\in A}\sum_{j\in B}
\frac{(\mathbf{F}_j - \mathbf{F}_i)\cdot\hat{\mathbf{r}}_{ij}}{2\,r_{ij}^2}\;
\theta(r - r_{ij})\right\rangle ,
\f]

with \f$\theta\f$ the step function. Every pair contributes to every bin
beyond its distance instead of to one bin, so the noise does not grow when
the bins are made narrower, and the result is smooth at any resolution.

In molecular dynamics the forces of the integrator are reused. In Monte Carlo
a full evaluation of the forces is made at each sample, intramolecular forces
included, so flexible molecules are treated correctly; the sampling interval
should reflect that cost. The estimator is normalized by the average number
of pairs, so that it remains valid when the number of molecules or the volume
fluctuates. It is refused with an external field or with polarization, whose
forces are not evaluated; the conventional function is available there.

----------------------------------------------------------------------------------

## Evolution during the run {#evolution-during-the-run}

Three properties record a quantity as a function of the cycle number, over
all stages of the run, rather than its average:

- `"ComputeNumberOfMoleculesEvolution"`: the number of molecules of each
  component;
- `"ComputeVolumeEvolution"`: the volume of the box;
- the conserved energy of molecular dynamics, with the terms of the running
  energy. This one has no input key; it is switched on from the Python
  interface (`set_property_conserved_energy_evolution`).

They are written to the directories `number_of_molecules_evolution`,
`volume_evolution` and `conserved_energy`. They answer the question that averages cannot: *has the system
equilibrated?* The number of molecules during the initialization of an
adsorption run rises and levels off; production should start after the
plateau is reached. The same holds for the volume in a constant-pressure run.
The conserved energy of molecular dynamics should fluctuate about a constant
value; its drift measures the quality of the integration.

----------------------------------------------------------------------------------

## Dynamics {#dynamics}

These properties require molecular dynamics and a fixed number of molecules.

### Mean-squared displacement {#property-msd}

`"ComputeMSD"` (system). Output directory: `msd`.

The **self** mean-squared displacement of a component is the average over its
molecules and over time origins of the squared displacement of the centre of
mass,

\f[
\text{MSD}(t) = \left\langle \left|\mathbf{r}_i(t_0 + t) - \mathbf{r}_i(t_0)\right|^2 \right\rangle ,
\f]

written in total and per direction. At short times the motion is ballistic
(\f$\propto t^2\f$), at intermediate times a molecule may be confined by its
neighbours or its cage (a plateau), and at long times it is diffusive
(\f$\propto t\f$). The self-diffusion coefficient is the slope in the
diffusive regime (the Einstein relation),

\f[
D_S = \lim_{t\to\infty}\frac{\text{MSD}(t)}{2\,d\,t} ,
\f]

with \f$d\f$ the number of dimensions over which the displacement is taken:
three for the total, one for a single direction. On a log-log plot the
diffusive regime has slope one; only that part may be fitted.

The **collective** (Onsager) displacements are written per pair of
components,

\f[
\text{MSD}_{ij}(t) = \frac{1}{N}\left\langle \Delta\mathbf{R}_i(t)\cdot\Delta\mathbf{R}_j(t)\right\rangle,
\qquad \Delta\mathbf{R}_i(t) = \sum_{k\in i}\left[\mathbf{r}_k(t_0+t) - \mathbf{r}_k(t_0)\right] ,
\f]

with \f$N\f$ the total number of molecules. Their slopes are the Onsager
coefficients, from which the Maxwell-Stefan and the Fick diffusion
coefficients follow. The collective functions are sums over molecules, not
averages, so they have far poorer statistics than the self function: one
sample per time origin instead of one per molecule.

*The order-\f$n\f$ algorithm.* Storing every position at every time origin
would need memory proportional to the length of the run. The displacements
are instead kept in blocks of increasing coarseness: the first block holds
the last `"NumberOfBlockElementsMSD"` samples, the second holds sums over
that many samples, and so on. This covers times from one sampling interval to
the length of the run, with points that are equidistant on a logarithmic
scale, at a memory cost that grows with the logarithm of the run length.

### Velocity autocorrelation function {#property-vacf}

`"ComputeVACF"` (system). Output directory: `vacf`.

The velocity autocorrelation function,

\f[
C_v(t) = \left\langle \mathbf{v}_i(t_0 + t)\cdot\mathbf{v}_i(t_0)\right\rangle ,
\f]

describes how long a molecule remembers its velocity. It starts at
\f$d\,k_BT/m\f$, decays through collisions, and in a dense fluid or a cage
becomes negative when molecules rebound. Its time integral is the
self-diffusion coefficient (the Green-Kubo relation),

\f[
D_S = \frac{1}{d}\int_0^\infty C_v(t)\,dt ,
\f]

equivalent to the Einstein relation. The collective functions
\f$\langle\mathbf{V}_i(t)\cdot\mathbf{V}_j(0)\rangle/N\f$, with
\f$\mathbf{V}_i\f$ the sum of the velocities of the molecules of component
\f$i\f$, give the Onsager coefficients in the same way. The Fourier transform
of the function is the vibrational density of states.

The function is accumulated with several staggered buffers: time origins are
started at regular intervals, so that `"NumberOfBuffersVACF"` origins are in
use at any moment, each followed for `"BufferLengthVACF"` samples. The
sampling interval must resolve the fastest motion, so it is short compared
to that of the mean-squared displacement.

*Choosing between the two.* The mean-squared displacement is the robust route
to a diffusion coefficient. The integral of the autocorrelation function
converges slowly because of its long-time tail, but the function itself
shows the mechanism of the motion.

----------------------------------------------------------------------------------

## Conformation of molecules {#conformation-of-molecules}

These analyses are switched on per component and are available in Monte
Carlo and molecular dynamics. Their accumulated data are stored in the binary
restart files.

### Bond, bend and torsion histograms {#property-molecule-properties}

`"ComputeMoleculeProperties"` (component). Output directory:
`molecule_properties`.

Probability histograms of the internal coordinates of a molecule: every bond
length, every bend angle between 0° and 180°, and every torsion angle between
−180° and 180°. When the component has end-to-end atoms the distribution of
the end-to-end distance is added.

For an isolated coordinate with potential \f$u(x)\f$ the distribution is the
Boltzmann factor times the volume element of the coordinate,

\f[
P(l) \propto l^2\,e^{-\beta u(l)}, \qquad
P(\theta) \propto \sin\theta\; e^{-\beta u(\theta)}, \qquad
P(\phi) \propto e^{-\beta u(\phi)} ,
\f]

which is the test of a correct sampling of flexible molecules. Deviations in
a dense phase or in a pore show how the environment changes the conformation:
a shift of the trans-gauche populations of a chain in a narrow channel, for
instance.

### Shape: the gyration tensor {#property-molecule-shape}

`"ComputeMoleculeShape"` (component). Output directory: `molecule_shape`.

The shape of a molecule is characterized by its gyration tensor,

\f[
S_{\alpha\beta} = \sum_i w_i\,(r_{i\alpha} - r_{\text{cm},\alpha})(r_{i\beta} - r_{\text{cm},\beta}) ,
\f]

with uniform weights \f$w_i = 1/N\f$ (the convention of polymer physics) or
mass weights. Its eigenvalues \f$\lambda_1 \ge \lambda_2 \ge \lambda_3\f$ are
the squared extents of the molecule along its principal axes, and from them
follow

| Descriptor | Definition | Meaning |
|---|---|---|
| radius of gyration | \f$R_g^2 = \lambda_1 + \lambda_2 + \lambda_3\f$ | overall size |
| asphericity | \f$b = \lambda_1 - \frac{1}{2}(\lambda_2 + \lambda_3)\f$ | 0 for a sphere |
| acylindricity | \f$c = \lambda_2 - \lambda_3\f$ | 0 for a cylinder |
| relative shape anisotropy | \f$\kappa^2 = (b^2 + \frac{3}{4}c^2)/R_g^4\f$ | 0 sphere, 1 rod |
| prolateness | \f$S = 27\prod_i(\lambda_i - \bar\lambda)/R_g^6\f$ | \f$-\frac{1}{4}\f$ oblate, 0 sphere, 2 prolate |
| hydrodynamic radius | \f$R_h^{-1} = \langle N^{-2}\sum_{i\ne j} 1/r_{ij}\rangle\f$ | Kirkwood approximation |

with \f$\bar\lambda = R_g^2/3\f$. The ratio \f$\sqrt{\langle R_g^2\rangle}/R_h\f$
is 1.5 for a Gaussian chain and 0.78 for a hard sphere, and
\f$\langle R^2\rangle/\langle R_g^2\rangle\f$ (with \f$R\f$ the end-to-end
distance) is 6 for an ideal chain.

Averages of ratios and ratios of averages differ. Both are reported: the
per-molecule averages \f$\langle\kappa^2\rangle\f$ and \f$\langle S\rangle\f$,
and the ensemble-ratio forms of Theodorou and Suter.

The *ensemble-averaged* tensor \f$\langle S_{\alpha\beta}\rangle\f$ in the
laboratory frame is written too. For an isotropic fluid its three eigenvalues
are equal; in a framework its principal axes show along which direction of
the box the chains are aligned. For components with repeat units the
descriptors are also given per monomer along the chain, which shows, for
instance, whether end monomers are more extended than interior ones.

### Backbone: chain statistics {#property-molecule-backbone}

`"ComputeMoleculeBackbone"` (component). Output directory:
`molecule_backbone`.

The backbone is the shortest path of bonds between the two ends of a chain.
With bond vectors \f$\mathbf{b}_i = \mathbf{r}_{i+1} - \mathbf{r}_i\f$ along
it, the following are sampled.

- **Internal distances** \f$\langle r^2(k)\rangle\f$ between beads \f$k\f$
  bonds apart. For an ideal chain \f$\langle r^2(k)\rangle/(k\langle l\rangle^2)\f$
  is constant; its rise at small \f$k\f$ shows the local stiffness.
- **Bond correlation** \f$C(k) = \langle\hat{\mathbf{b}}_i\cdot\hat{\mathbf{b}}_{i+k}\rangle\f$,
  which decays as \f$e^{-k\langle l\rangle/l_p}\f$ with the persistence length
  \f$l_p\f$ for a worm-like chain.
- **The form factor**, the scattering function of a single chain,

\f[
P(q) = \left\langle \frac{1}{N^2}\sum_{i,j}\frac{\sin(q\,r_{ij})}{q\,r_{ij}}\right\rangle ,
\f]

  on a logarithmic grid of wave vectors, next to the Debye function of a
  Gaussian chain with the same radius of gyration. It is what small-angle
  scattering measures.

The summary derives from them

| Quantity | Definition |
|---|---|
| characteristic ratio | \f$C_N = \langle R^2\rangle/((N_b - 1)\langle l\rangle^2)\f$ |
| Kuhn length | \f$b_K = \langle R^2\rangle/R_\text{max}\f$, with \f$R_\text{max}\f$ the contour length |
| persistence length | from the projection of the chain on its first bond, and from a fit to \f$\ln C(k)\f$ |
| Flory exponent \f$\nu\f$ | the slope of \f$\ln\langle r^2(k)\rangle\f$ against \f$\ln k\f$, divided by two |

The Flory exponent distinguishes a collapsed chain (\f$\nu = 1/3\f$), an
ideal chain (\f$1/2\f$), a swollen chain in a good solvent (\f$\approx 0.588\f$)
and a rod (\f$1\f$). It needs a chain long enough to have a scaling regime.

### End-to-end autocorrelation and relaxation time {#property-end-to-end-acf}

`"ComputeEndToEndACF"` (component). Output directory: `end_to_end_acf`.

The autocorrelation function of the end-to-end vector \f$\mathbf{R}\f$ of a
chain,

\f[
C(t) = \left\langle \mathbf{R}(0)\cdot\mathbf{R}(t)\right\rangle ,
\f]

decays on the longest relaxation time of the chain, \f$\tau_R\f$. That time
is the spacing of statistically independent samples of the conformation: a
run of length \f$T\f$ contains about \f$T/(2\tau_R)\f$ independent samples
per chain. It is therefore the quantity that tells whether a simulation of a
polymer is long enough. Three estimates of \f$\tau_R\f$ are reported: the
integral of \f$C(t)/C(0)\f$ up to its first zero, the time at which it falls
to \f$1/e\f$, and the decay time of an exponential fitted to its tail.

The vector function also decays when a chain of unchanged shape merely
tumbles. Two rotation-invariant functions are therefore accumulated on the
same time lags: the autocorrelation of the fluctuations of the squared
end-to-end distance and of the squared radius of gyration,

\f[
C_{R^2}(t) = \frac{\left\langle \delta R^2(0)\,\delta R^2(t)\right\rangle}{\left\langle(\delta R^2)^2\right\rangle},
\qquad \delta R^2 = R^2 - \langle R^2\rangle .
\f]

These decay only when the chain changes shape. A vector relaxation time much
shorter than the conformational ones means the chain rotates while its fold
persists; comparable times mean the conformation itself relaxes. For a
Gaussian chain in the Rouse model \f$C_{R^2}\f$ decays about twice as fast as
the vector function.

The function is accumulated with the order-\f$n\f$ scheme of the mean-squared
displacement. In molecular dynamics the time lag is in picoseconds. In Monte
Carlo it is in cycles, and the function then measures how efficiently the
moves decorrelate the conformations, which is a way of comparing moves.

----------------------------------------------------------------------------------

## Mechanical properties {#mechanical-properties}

### Elastic constants from fluctuations {#property-elastic-fluctuation}

`"ComputeElasticConstantsFromFluctuations"` (system), in molecular dynamics
at constant volume and temperature.

At finite temperature the isothermal elastic tensor is the second derivative
of the free energy with respect to strain. In a simulation at fixed cell it
is the sum of three terms,

\f[
C_{\alpha\beta\gamma\delta} = \left\langle C^\text{Born}_{\alpha\beta\gamma\delta}\right\rangle
- \frac{V}{k_B T}\left[\langle\sigma_{\alpha\beta}\,\sigma_{\gamma\delta}\rangle
- \langle\sigma_{\alpha\beta}\rangle\langle\sigma_{\gamma\delta}\rangle\right]
+ \frac{N k_B T}{V}\left(\delta_{\alpha\gamma}\delta_{\beta\delta} + \delta_{\alpha\delta}\delta_{\beta\gamma}\right) .
\f]

The **Born term** is the average second derivative of the potential energy
with respect to a homogeneous strain: the stiffness if all atoms followed the
deformation affinely. The **fluctuation term** is the covariance of the
stress; it accounts for the relaxation of the atoms and always softens the
material. The **kinetic term** is the ideal-gas contribution. The output
gives the three terms separately, the stiffness matrix in Voigt notation with
its error, the Born stability eigenvalues and the bulk, shear and Young's
moduli.

*Notes.* The fluctuation term converges slowly: the blocks must be much
longer than the correlation time of the stress, and the result should be
checked against the length of the run, the sampling interval and the size of
the system. The Born term needs the second derivatives of the energy, which
cost considerably more than the forces; the sampling interval is set
separately for that reason. The cell should first be equilibrated at the
temperature and mean pressure of interest.

### Elastic constants, normal modes and phonons at a minimum {#property-minimization}

The [minimization driver](\ref driver_minimization) computes the
zero-temperature counterparts at an energy minimum: the relaxed elastic
tensor (`"ComputeElasticConstants"`), the normal-mode frequencies
(`"ComputeNormalModes"`) and the phonon band structure
(`"ComputePhononDispersion"`). The static elastic tensor is the
zero-temperature limit of the fluctuation formula, in which the stress
covariance becomes the relaxation term.

----------------------------------------------------------------------------------

## Reaction kinetics {#reaction-kinetics}

### Radical waiting times {#property-radical-waiting-times}

`"ComputeRadicalWaitingTimes"` (system), with the
[polymerization driver](\ref driver_polymerization). Output directory:
`radical_waiting_times`.

For every active group (a living chain end) the time between two successive
reactions is recorded at the moment it reacts, in a histogram with
logarithmically spaced bins from \f$10^{-3}\f$ to \f$10^{6}\f$ ps, separately
per reaction rule. The ages of the groups that are still active when the
property is written, the *open* intervals, are given as a second histogram.

For a first-order process with rate \f$k\f$ the waiting times are
exponentially distributed,

\f[
P(\Delta t) = k\,e^{-k\,\Delta t} ,
\f]

with mean \f$1/k\f$. A single exponential therefore means a well-defined
propagation rate. A stretched tail, a second slow component, or an
accumulation of old active groups that no longer react is the signature of
*trapping*: radicals buried in the network that monomers cannot reach. The
conversion at which this sets in is read from successive output files.

----------------------------------------------------------------------------------

## Configurations {#configurations}

Two options write the configuration itself, for visualization or for use in
another program:

- `"OutputPDBMovie"` writes snapshots every `"SampleMovieEvery"` cycles to
  the directory `movies`, as a PDB file with one model per snapshot
  (`movies/movie.s<system>.pdb`). The positions can be wrapped into the box or
  left unwrapped.
- `"WriteLammpsData"` writes the configuration as a LAMMPS data file, with
  its topology and tables, to the directory `lammps`.

Movies of the normal modes are written by the
[minimization driver](\ref driver_minimization).

----------------------------------------------------------------------------------

## References

- D. Frenkel and B. Smit, *Understanding Molecular Simulation*, 3rd edition,
  Academic Press (2023).
- M.P. Allen and D.J. Tildesley, *Computer Simulation of Liquids*, 2nd
  edition, Oxford University Press (2017).
- H. Flyvbjerg and H.G. Petersen, *Error estimates on averages of correlated
  data*, J. Chem. Phys. **91**, 461 (1989).
- T.J.H. Vlugt, E. García-Pérez, D. Dubbeldam, S. Ban and S. Calero,
  *Computing the Heat of Adsorption using Molecular Simulations: The Effect of
  Strong Coulombic Interactions*, J. Chem. Theory Comput. **4**, 1107 (2008).
- B. Widom, *Some Topics in the Theory of Fluids*,
  J. Chem. Phys. **39**, 2808 (1963).
- D. Borgis, R. Assaraf, B. Rotenberg and R. Vuilleumier, *Computation of pair
  distribution functions and three-dimensional densities with a reduced
  variance principle*, Mol. Phys. **111**, 3486 (2013).
- D. Dubbeldam, D.C. Ford, D.E. Ellis and R.Q. Snurr, *A new perspective on the
  order-n algorithm for computing correlation functions*,
  Mol. Simul. **35**, 1084 (2009).
- R. Krishna and J.M. van Baten, *Diffusion of Alkane Mixtures in Zeolites:
  Validating the Maxwell-Stefan Formulation Using MD Simulations*,
  J. Phys. Chem. B **109**, 6386 (2005).
- D.N. Theodorou and U.W. Suter, *Shape of unperturbed linear polymers:
  polypropylene*, Macromolecules **18**, 1206 (1985).
- M. Rubinstein and R.H. Colby, *Polymer Physics*, Oxford University Press
  (2003).
- D.R. Squire, A.C. Holt and W.G. Hoover, *Isothermal elastic constants for
  argon. Theory and Monte Carlo calculations*, Physica **42**, 388 (1969).
- J.F. Lutsko, *Generalized expressions for the calculation of elastic
  constants by computer simulation*, J. Appl. Phys. **65**, 2991 (1989).
- J. Rouquerol, P. Llewellyn and F. Rouquerol, *Is the BET equation
  applicable to microporous adsorbents?*,
  Stud. Surf. Sci. Catal. **160**, 49 (2007).
