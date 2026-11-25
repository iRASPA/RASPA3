# Molecular dynamics
\page driver_molecular_dynamics Molecular dynamics

`"SimulationType" : "MolecularDynamics"`

The molecular-dynamics driver integrates the classical equations of motion.
Next to the equilibrium averages that Monte Carlo gives, it gives the
quantities that depend on time: mean-squared displacements and diffusion
coefficients, velocity autocorrelation functions, and relaxation times.

## Equations of motion

For atoms with positions \f$\mathbf{r}_i\f$, momenta \f$\mathbf{p}_i\f$ and
masses \f$m_i\f$ Hamilton's equations are

\f[
\dot{\mathbf{r}}_i = \frac{\mathbf{p}_i}{m_i}, \qquad
\dot{\mathbf{p}}_i = \mathbf{F}_i = -\frac{\partial U}{\partial \mathbf{r}_i} .
\f]

They conserve the energy \f$H = \sum_i p_i^2/2m_i + U\f$, so a trajectory
samples the micro-canonical (\f$NVE\f$) ensemble. Time averages over a long
trajectory are taken as ensemble averages (the ergodic hypothesis).

## The integrator

The equations are integrated with the velocity-Verlet scheme. Written with
the Liouville operator \f$iL = iL_r + iL_p\f$, one time step \f$\Delta t\f$ is
the symmetric (Trotter) factorization

\f[
e^{iL\Delta t} \approx e^{iL_p \Delta t/2}\; e^{iL_r \Delta t}\; e^{iL_p \Delta t/2},
\f]

that is: half a step of the momenta with the present forces, a full step of
the positions, new forces, and the second half step of the momenta. Each
factor is an exact flow of a part of the Hamiltonian, so the scheme is
*symplectic* and time reversible. A symplectic integrator does not conserve
\f$H\f$ exactly, but it conserves a nearby "shadow" Hamiltonian, so the
energy fluctuates around a constant with an amplitude of order
\f$\Delta t^2\f$ instead of drifting. The absence of drift is the practical
test of a time step.

**Rigid molecules.** A rigid molecule has the position of its centre of mass
and an orientation, stored as a unit quaternion, as degrees of freedom, with
the linear and angular momentum conjugate to them. The rotational part is
integrated with the symplectic quaternion scheme of Miller *et al.*
(NO_SQUISH), which splits the free rotation of an asymmetric body into exactly
solvable rotations about the body axes.

**Constraints.** Bonds, and bends between two fixed bonds, that are marked as
fixed in a flexible molecule are held with SHAKE (positions) and RATTLE
(velocities): after the unconstrained step the positions are corrected along
the bonds until every constrained distance is restored, and the components of
the velocities along the bonds are removed. Every independent constraint
removes one degree of freedom from the temperature.

## Thermostats

To sample at a temperature the system is coupled to a heat bath.

**Nosé-Hoover chain.** The momenta get a friction term whose coefficient is
itself a dynamical variable, driven by the difference between the
instantaneous kinetic energy and its target,

\f[
\dot{\mathbf{p}}_i = \mathbf{F}_i - \frac{p_{\xi_1}}{Q_1}\,\mathbf{p}_i, \qquad
\dot{p}_{\xi_1} = \sum_i \frac{p_i^2}{m_i} - N_f k_B T - \frac{p_{\xi_2}}{Q_2}\,p_{\xi_1} ,
\f]

and every thermostat variable is in turn thermostatted by the next one of the
chain. The dynamics is deterministic and time reversible and generates the
canonical distribution; the chain repairs the lack of ergodicity of a single
Nosé-Hoover thermostat for small or stiff systems. The thermostat masses
follow from one time scale, \f$Q = N_f k_B T \tau^2\f$. The chain is
integrated with a Yoshida-Suzuki factorization inside multiple small steps.
The conserved quantity is the energy of the system plus that of the bath.

**Langevin.** Alternatively every degree of freedom gets a friction and a
random force that satisfy the fluctuation-dissipation theorem,

\f[
\dot{\mathbf{p}}_i = \mathbf{F}_i - \gamma\,\mathbf{p}_i
+ \sqrt{2\gamma m_i k_B T}\;\boldsymbol{\eta}_i(t) ,
\f]

with \f$\boldsymbol{\eta}\f$ white noise. It thermalizes every degree of
freedom locally, at the price of perturbing the dynamics on the time scale
\f$1/\gamma\f$. The heat exchanged with the bath is tracked, so that a
conserved quantity is still available as a check.

## Barostats

At constant pressure the volume becomes a dynamical variable with a mass. In
the Martyna-Tobias-Klein formulation the positions and momenta are scaled by
the rate of change of the logarithm of the volume, which is driven by the
difference between the instantaneous and the external pressure,

\f[
\dot{\mathbf{r}}_i = \frac{\mathbf{p}_i}{m_i} + \frac{p_\epsilon}{W}\,\mathbf{r}_i, \qquad
\dot V = \frac{d\,V\,p_\epsilon}{W}, \qquad
\dot p_\epsilon = d\,V\,(P_\text{inst} - P_\text{ext}) + \frac{d}{N_f}\sum_i \frac{p_i^2}{m_i} .
\f]

`"NPT"` scales the three box lengths together. `"NPTPR"` (Parrinello-Rahman
in the Martyna-Tobias-Klein form) makes the cell matrix dynamical, driven by
the pressure *tensor*, so that the shape of the box can change too; the
`"CellType"` restricts which components may move. The conserved quantity is
the extended enthalpy. The instantaneous pressure is that of the points the
barostat scales: the centres of mass of rigid molecules, and the atoms of
flexible ones.

## Open ensembles

In the `"MuVT"`, `"MuPT"` and `"MuPTPR"` ensembles molecules are exchanged
with a reservoir during the dynamics: a configurational-bias insertion or
deletion is attempted at regular intervals, with the grand-canonical
acceptance rule, between ordinary time steps of the corresponding closed
ensemble. An inserted molecule gets a velocity from the Maxwell-Boltzmann
distribution.

## Stages

Initialization uses Monte Carlo cycles, which relax a starting configuration
far more robustly than dynamics does. At the start of the equilibration the
velocities are drawn from the Maxwell-Boltzmann distribution at the external
temperature, and from then on a cycle is one time step.

## Transport properties

**Self-diffusion** follows from the mean-squared displacement of the
individual molecules (the Einstein relation), or equivalently from the
integral of the velocity autocorrelation function (Green-Kubo),

\f[
D_S = \lim_{t\to\infty} \frac{1}{2 d\,t}\,
\left\langle \left|\mathbf{r}_i(t) - \mathbf{r}_i(0)\right|^2 \right\rangle
= \frac{1}{d}\int_0^\infty \langle \mathbf{v}_i(t)\cdot\mathbf{v}_i(0)\rangle\,dt ,
\f]

with \f$d\f$ the dimension. **Collective** (corrected, Maxwell-Stefan)
diffusion uses the displacement of the centre of mass of all molecules of a
component,

\f[
D_C = \lim_{t\to\infty} \frac{1}{2 d\,N\,t}\,
\left\langle \left|\sum_{i=1}^{N} \left[\mathbf{r}_i(t) - \mathbf{r}_i(0)\right]\right|^2 \right\rangle .
\f]

The mean-squared displacement is sampled with an order-\f$n\f$ scheme: blocks
of increasing length store coarse-grained displacements, so that times from
one sampling interval up to the length of the run are covered with a memory
that grows with the logarithm of the run length. Only the linear regime of
the mean-squared displacement, after the ballistic regime and the cage
regime, gives a diffusion coefficient.

## Things to check

- **Drift of the conserved quantity.** It should fluctuate without trend. A
  drift means the time step is too large, the forces are discontinuous (a
  truncated instead of a shifted potential), or the constraints are not
  converged.
- **Temperature.** The average kinetic temperature must equal the external
  temperature; the translational and rotational temperatures are reported
  separately and must agree.
- **Rare events.** When molecules cross a barrier less than a few times
  during the run, the diffusion coefficient from the mean-squared displacement
  is not converged. The [dcTST drivers](\ref driver_dctst) are made for that
  regime.

## References

- L. Verlet, *Computer "Experiments" on Classical Fluids. I*,
  Phys. Rev. **159**, 98 (1967).
- M. Tuckerman, B.J. Berne and G.J. Martyna, *Reversible multiple time scale
  molecular dynamics*, J. Chem. Phys. **97**, 1990 (1992).
- G.J. Martyna, M.L. Klein and M. Tuckerman, *Nosé-Hoover chains: The canonical
  ensemble via continuous dynamics*, J. Chem. Phys. **97**, 2635 (1992).
- G.J. Martyna, D.J. Tobias and M.L. Klein, *Constant pressure molecular
  dynamics algorithms*, J. Chem. Phys. **101**, 4177 (1994).
- G.J. Martyna, M.E. Tuckerman, D.J. Tobias and M.L. Klein, *Explicit
  reversible integrators for extended systems dynamics*,
  Mol. Phys. **87**, 1117 (1996).
- T.F. Miller III, M. Eleftheriou, P. Pattnaik, A. Ndirango, D. Newns and
  G.J. Martyna, *Symplectic quaternion scheme for biophysical molecular
  dynamics*, J. Chem. Phys. **116**, 8649 (2002).
- J.-P. Ryckaert, G. Ciccotti and H.J.C. Berendsen, *Numerical integration of
  the Cartesian equations of motion of a system with constraints*,
  J. Comput. Phys. **23**, 327 (1977).
- H.C. Andersen, *Rattle: A "velocity" version of the Shake algorithm*,
  J. Comput. Phys. **52**, 24 (1983).
- D. Dubbeldam, D.C. Ford, D.E. Ellis and R.Q. Snurr, *A new perspective on the
  order-n algorithm for computing correlation functions*,
  Mol. Simul. **35**, 1084 (2009).
- D. Dubbeldam and R.Q. Snurr, *Recent Developments in the Molecular Modeling
  of Diffusion in Nanoporous Materials*, Mol. Simul. **33**, 305 (2007).
