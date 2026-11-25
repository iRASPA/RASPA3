# Nonlocal density functional theory
\page driver_nldft Nonlocal density functional theory

`"SimulationType" : "NLDFT"`

Classical density functional theory computes an adsorption isotherm without
sampling: it finds the equilibrium *density profile* of the adsorbate
directly, by minimizing a free-energy functional. The driver does this on the
energy grid of the actual framework, for nitrogen-like probes, and reports
the isotherm together with a BET analysis of it.

## Forward instead of inverse

In experimental porosimetry density functional theory is used for the
*inverse* problem: an isotherm is measured, the structure is unknown, and a
library of isotherms of ideal slit or cylinder pores (a kernel) is fitted to
it to obtain a pore-size distribution. Here the structure is known exactly,
so the *forward* problem is solved: the external potential is the field of
the framework itself, and the isotherm is predicted.

## The grand potential functional

The central theorem of density functional theory states that the grand
potential is a functional of the one-body density \f$\rho(\mathbf{r})\f$ that
is minimal at the equilibrium density, where its value is the grand potential
of the system. The functional is

\f[
\Omega[\rho] = F_\text{id}[\rho] + F_\text{hs}[\rho] + F_\text{att}[\rho]
+ \int \rho(\mathbf{r})\left[V_\text{ext}(\mathbf{r}) - \mu\right] d\mathbf{r} .
\f]

The four terms are:

- **The ideal gas**, which is known exactly,

\f[
\beta F_\text{id}[\rho] = \int \rho(\mathbf{r})\left[\ln\!\left(\Lambda^3\rho(\mathbf{r})\right) - 1\right] d\mathbf{r} .
\f]

- **The hard-sphere repulsion**, by fundamental measure theory (below).
- **The attraction**, in mean field,

\f[
F_\text{att}[\rho] = \frac{1}{2}\iint \rho(\mathbf{r})\,\rho(\mathbf{r}')\,
u_\text{att}(|\mathbf{r} - \mathbf{r}'|)\,d\mathbf{r}\,d\mathbf{r}' ,
\f]

  with \f$u_\text{att}\f$ the attractive part of the Lennard-Jones potential.
  For the spherical probe this is the Weeks-Chandler-Andersen division
  (constant at the well depth inside the minimum, the Lennard-Jones tail
  outside); for the nitrogen molecule it is the 12-6 potential between the
  two sites from \f$r = \sigma\f$ outward.
- **The external field** of the framework, its dispersion and electrostatic
  potential on the grid, and the chemical potential of the reservoir.

"Nonlocal" refers to the repulsive term. A *local* functional, which uses the
free energy of the uniform fluid at the density of each point, cannot
describe the layering of a fluid against a wall, where the local density
exceeds that of close packing. A nonlocal functional uses densities averaged
over the size of a molecule.

## Fundamental measure theory

In Rosenfeld's fundamental measure theory the hard-sphere free energy is

\f[
\beta F_\text{hs}[\rho] = \int \Phi\!\left(\{n_\alpha(\mathbf{r})\}\right) d\mathbf{r} ,
\f]

a local function of *weighted densities*,

\f[
n_\alpha(\mathbf{r}) = \int \rho(\mathbf{r}')\,\omega_\alpha(\mathbf{r} - \mathbf{r}')\,d\mathbf{r}' ,
\f]

which are convolutions of the density with weight functions that describe the
geometry of a sphere: its volume (\f$n_3\f$, the local packing fraction), its
surface (\f$n_2\f$), and a vector weight for the surface normal
(\f$\mathbf{n}_2\f$). The driver uses the *White Bear* form of \f$\Phi\f$,
whose uniform limit is the accurate Carnahan-Starling equation of state, with
the third term in the Rosenfeld-Tarazona form, which keeps the free-energy
density non-negative and the functional bounded below.

## The equation for the density

Setting the functional derivative of \f$\Omega\f$ to zero gives the
Euler-Lagrange equation,

\f[
\rho(\mathbf{r}) = \rho_b\,\exp\!\left\{-\beta\left[V_\text{ext}(\mathbf{r})
+ \mu_\text{hs}[\rho](\mathbf{r}) - \mu_\text{hs,b}
+ \int \rho(\mathbf{r}')\,u_\text{att}(|\mathbf{r}-\mathbf{r}'|)\,d\mathbf{r}' - \rho_b\,a\right]\right\} ,
\f]

with \f$\rho_b\f$ the density of the bulk fluid at the same chemical
potential and \f$a = \int u_\text{att}\,d\mathbf{r}\f$. It reads as a
Boltzmann factor of an effective potential: the external field, plus the
excess chemical potential of the hard spheres and the mean-field attraction
at that point, each relative to the bulk.

Two limits are instructive. At low pressure the terms in \f$\rho\f$ vanish
and \f$\rho(\mathbf{r}) = \rho_b\,e^{-\beta V_\text{ext}(\mathbf{r})}\f$: the
Henry regime, identical to what Monte Carlo gives. At high pressure the
hard-sphere term is what stops the density from growing in the wells, and
produces the layering and the filling of the pore.

## Numerical solution

The density lives on the grid of the framework's energy grid, periodic in the
cell. All convolutions, those of the weighted densities and that of the
attractive tail, are products in Fourier space, so one evaluation of the
functional costs a fixed number of fast Fourier transforms, for any cell
shape.

The equation is not iterated directly, which is unstable at high density.
Instead \f$\Omega\f$ is minimized with a quasi-Newton method (L-BFGS) in the
variable \f$\ln\rho\f$, which keeps the density positive, and every step is
required to lower \f$\Omega\f$. The weight functions are laid onto the grid
as non-negative kernels in real space, so that the weighted densities cannot
become negative or exceed their bounds on any grid.

## The molecular density

A spherical probe has a density \f$\rho(\mathbf{r})\f$. For a linear molecule
such as nitrogen the density also depends on the orientation,
\f$\rho(\mathbf{r}, \omega)\f$, and the external field becomes the energy of
the molecule at position \f$\mathbf{r}\f$ with orientation \f$\omega\f$. The
ideal term is then that of an ideal rotator, and the low-pressure limit is

\f[
\rho(\mathbf{r}) = \rho_b\,\left\langle e^{-\beta U(\mathbf{r},\omega)}\right\rangle_\omega ,
\f]

the same orientational average that defines the Henry coefficient in Monte
Carlo. This matters in narrow pores, where a molecule fits only in certain
orientations: a spherical probe cannot represent that.

## The bulk fluid and the relative pressure

The same functional at uniform density is an equation of state of the bulk
fluid. Its vapour-liquid coexistence at the temperature of the run is found
once, by a Maxwell construction, and gives the saturation pressure \f$P_0\f$
*of the model*. The relative pressure \f$x = P/P_0\f$ is measured against
that value, not against the experimental one. This is essential for
consistency: only with the model's own \f$P_0\f$ does \f$x \to 1\f$ mean
condensation of the bulk fluid in the model.

## Sweeps and hysteresis

The pressure is swept upward, each point starting from the converged profile
of the previous one. Where that sweep fills a pore, the pressure is also
swept back down from the filled profile. Inside a hysteresis loop both
profiles are local minima of \f$\Omega\f$, metastable in opposite directions.
The reported isotherm is the *equilibrium* one: at each pressure, the profile
with the lower grand potential, which is the same criterion that locates the
bulk coexistence.

Layering, micropore filling and capillary condensation come out of the real
geometry; nothing is assembled from ideal pores.

## BET analysis

The BET surface area is obtained by fitting the BET equation to the computed
isotherm over the pressure range selected by the Rouquerol consistency
criteria, the procedure applied to a measured or simulated isotherm, so that
the number can be compared with experiment on equal terms.

## What the method does and does not include

- **No sampling.** The result is the minimum of a functional, so the isotherm
  has no statistical noise; its error is that of the functional.
- **Mean-field attraction.** Correlations due to the attraction are
  neglected, which is the main approximation of the method. It is not a
  replacement for Monte Carlo for strongly associating fluids.
- **A model fluid.** The fluid is a hard sphere with a Lennard-Jones tail (or
  a fused dumbbell for nitrogen), with parameters fitted to the bulk fluid.
  The interactions between adsorbate molecules do not include electrostatics.
- **Deep wells.** For a spherical probe the external potential is clipped
  from below, since a spherical density cannot represent the orientational
  confinement that makes such wells narrow.
- Exactly one system with a framework; the first component is the probe.

## References

- R. Evans, *The nature of the liquid-vapour interface and other topics in
  the statistical mechanics of non-uniform, classical fluids*,
  Adv. Phys. **28**, 143 (1979).
- Y. Rosenfeld, *Free-energy model for the inhomogeneous hard-sphere fluid
  mixture and density-functional theory of freezing*,
  Phys. Rev. Lett. **63**, 980 (1989).
- R. Roth, R. Evans, A. Lang and G. Kahl, *Fundamental measure theory for
  hard-sphere mixtures revisited: the White Bear version*,
  J. Phys.: Condens. Matter **14**, 12063 (2002).
- P. Tarazona, *Density Functional for Hard Sphere Crystals: A Fundamental
  Measure Approach*, Phys. Rev. Lett. **84**, 694 (2000).
- P.I. Ravikovitch, A. Vishnyakov, R. Russo and A.V. Neimark, *Unified
  Approach to Pore Size Characterization of Microporous Carbonaceous
  Materials from N2, Ar, and CO2 Adsorption Isotherms*,
  Langmuir **16**, 2311 (2000).
- J.D. Weeks, D. Chandler and H.C. Andersen, *Role of Repulsive Forces in
  Determining the Equilibrium Structure of Simple Liquids*,
  J. Chem. Phys. **54**, 5237 (1971).
- J. Rouquerol, P. Llewellyn and F. Rouquerol, *Is the BET equation
  applicable to microporous adsorbents?*,
  Stud. Surf. Sci. Catal. **160**, 49 (2007).
