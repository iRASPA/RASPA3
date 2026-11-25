# Replica-exchange molecular dynamics
\page driver_parallel_tempering_molecular_dynamics Replica-exchange molecular dynamics

`"SimulationType" : "ParallelTemperingMolecularDynamics"`
(aliases `"ReplicaExchangeMolecularDynamics"`, `"ParallelTemperingMD"`)

Replica-exchange molecular dynamics (REMD) is
[parallel tempering](\ref driver_parallel_tempering) in which the replicas are
propagated by [molecular dynamics](\ref driver_molecular_dynamics) instead of
Monte Carlo. It is the method of choice when the slow degrees of freedom are
collective and continuous, such as the conformation of a polymer chain, for
which dynamics moves all atoms together where single-molecule Monte Carlo
moves do not.

## The extended ensemble with momenta

A replica in the canonical ensemble has the weight
\f$e^{-\beta[K(\mathbf{p}) + U(\mathbf{r})]}\f$ over positions *and* momenta.
Exchanging the full phase-space points of replicas \f$A\f$ and \f$B\f$ would
give an acceptance that contains the kinetic energies as well. Sugita and
Okamoto showed that this can be avoided: when the momenta that travel with a
configuration are rescaled with the square root of the temperature ratio,

\f[
\mathbf{p}^\text{new} = \sqrt{\frac{T_\text{new}}{T_\text{old}}}\;\mathbf{p}^\text{old} ,
\f]

the kinetic energy changes by the factor \f$T_\text{new}/T_\text{old}\f$, so
\f$\beta_\text{new} K^\text{new} = \beta_\text{old} K^\text{old}\f$ and the
kinetic parts of the two Boltzmann factors cancel. The acceptance rule is
then that of Monte Carlo parallel tempering, on the potential energies only,

\f[
\text{acc} = \min\left(1,\; e^{(\beta_B - \beta_A)(U_B - U_A)}\right) .
\f]

The rescaling also means that the kinetic temperature of each replica is the
same before and after a swap: the configuration arrives with velocities that
already belong to its new temperature, and no re-equilibration is needed.

## What stays and what travels

- **Travels with the configuration:** the positions, the (rescaled) momenta,
  and for rigid molecules the orientations and angular momenta.
- **Stays with the replica:** the temperature, and the Nosé-Hoover chain. The
  thermostat variables describe the heat bath at that temperature, not the
  system, so they are not exchanged.

The conserved quantity of the extended system jumps at an accepted swap by
construction. Its reference value is therefore reset at a swap, so that the
reported drift measures the quality of the integration between swaps only.

## Why a thermostat is required

The acceptance rule is derived for replicas that sample the canonical
distribution at their temperatures. The driver therefore requires the
`"NVT"` ensemble. A micro-canonical replica has no temperature in that sense,
and the barostat and particle-exchange ensembles would let the volume or the
number of molecules differ between replicas, which the swap does not account
for.

## Stages and the swap interval

As in molecular dynamics, the initialization stages are Monte Carlo cycles.
At the start of the equilibration every replica gets Maxwell-Boltzmann
velocities at its own temperature, and from then on a cycle is one time step.
Swaps are attempted every `"ParallelTemperingSwapEvery"` cycles between
alternating sets of neighbour pairs.

For Monte Carlo a swap every ten cycles is reasonable, because a cycle
decorrelates the configuration appreciably. One MD step does not. After a
swap a replica needs time to lose the memory of the exchange and for its
potential energy to relax to the new temperature; attempting swaps much more
often than that only moves the same configurations back and forth. An
interval of a hundred to a thousand steps is usual.

## Diagnostics

The combined output has the same per-pair acceptance table and round-trip
statistics as [parallel tempering](\ref driver_parallel_tempering) (round
trips, mean round-trip time, the up-moving fraction \f$f(T)\f$), and in
addition the kinetic temperature and the drift of the conserved energy of
every replica. A replica whose kinetic temperature differs from its target
indicates a problem with the time step or the thermostat, not with the
exchange.

The time-dependent properties of a replica (mean-squared displacement,
velocity autocorrelation) are interrupted at every accepted swap, because the
configuration is replaced. They are meaningful only on time scales shorter
than the swap interval. The equilibrium averages are unaffected.

## Limits

- All replicas must have the same Hamiltonian and topology; systems with
  reactions or with pair, group or Gibbs fractional molecules are refused.
- The number of replicas needed grows with the square root of the number of
  degrees of freedom, as for any tempering method. For a solvated chain most
  of those degrees of freedom belong to the solvent.
- Binary restart files are written at the synchronization points and can be
  used to continue a run.

## References

- Y. Sugita and Y. Okamoto, *Replica-exchange molecular dynamics method for
  protein folding*, Chem. Phys. Lett. **314**, 141 (1999).
- U.H.E. Hansmann, *Parallel tempering algorithm for conformational studies of
  biological molecules*, Chem. Phys. Lett. **281**, 140 (1997).
- D.J. Earl and M.W. Deem, *Parallel tempering: Theory, applications, and new
  perspectives*, Phys. Chem. Chem. Phys. **7**, 3910 (2005).
- H.G. Katzgraber, S. Trebst, D.A. Huse and M. Troyer, *Feedback-optimized
  parallel tempering Monte Carlo*, J. Stat. Mech. P03018 (2006).
