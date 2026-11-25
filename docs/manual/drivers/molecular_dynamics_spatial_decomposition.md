# Spatial-decomposition molecular dynamics
\page driver_molecular_dynamics_spatial_decomposition Spatial-decomposition molecular dynamics

`"SimulationType" : "MolecularDynamicsSpatialDecomposition"`
(aliases `"SpatialDecompositionMolecularDynamics"`, `"MolecularDynamicsSD"`)

This driver integrates the same equations of motion as
[molecular dynamics](\ref driver_molecular_dynamics), with the same
integrators, thermostats, barostats and constraints. What differs is how the
forces are computed: for large boxes of molecules (tens of thousands of atoms
and more) the all-pairs loop and the Ewald sum over wave vectors are replaced
by methods whose cost grows linearly, or nearly so, with the number of atoms,
and the work is divided over threads by dividing space.

## Why the standard force routine does not scale

The force on an atom is a sum over all atoms within the cutoff \f$r_c\f$.
Finding them by testing every pair costs \f$O(N^2)\f$. The Ewald sum adds a
reciprocal-space part,

\f[
U_\text{rec} = \frac{1}{2V\epsilon_0}\sum_{\mathbf{k}\neq 0}
\frac{e^{-k^2/4\alpha^2}}{k^2}\,\left|S(\mathbf{k})\right|^2, \qquad
S(\mathbf{k}) = \sum_{j} q_j\, e^{i\mathbf{k}\cdot\mathbf{r}_j} ,
\f]

whose number of wave vectors grows with the volume at a fixed accuracy, so
that it costs \f$O(N^{3/2})\f$ at best. Both become the bottleneck long before
the integration itself does.

## Cell lists and Verlet lists

The box is divided into cells. An atom can only interact with atoms in its own
cell and in the cells near it, so the candidate pairs are found in
\f$O(N)\f$. The driver uses cells of a third of the search radius, which
makes the searched region follow the cutoff sphere more closely than cells of
the full radius would.

The pairs that are found are stored in a *Verlet list*, built with a radius
\f$r_c + r_\text{skin}\f$. The list stays valid as long as no atom has moved
more than half the skin, because two atoms outside the list radius cannot have
come within the cutoff before that. The lists are therefore rebuilt only when
needed, not every step; the skin trades the cost of rebuilding against the
length of the lists.

## Domain decomposition

The box is divided into as many sub-domains as there are threads, on a
\f$p_x \times p_y \times p_z\f$ grid chosen to minimize the surface between
them. The cut planes are placed on the atom positions at every rebuild, so
that every thread owns the same number of atoms even when the density is not
uniform.

Each thread works on a private, compact copy of its own atoms and of the
*ghost images* it needs: the atoms of neighbouring sub-domains, and periodic
images, that lie within the list radius of its domain. The periodic shift of
a ghost is resolved when the list is built, so the inner loop over pairs
needs no minimum-image operation. Every pair is evaluated once, by one
thread, with Newton's third law; the forces that end up on ghosts are
collected by the threads that own those atoms after the pair phase. The
threads meet only at these collection points, which is why the method scales:
the amount of shared data is proportional to the surface of a sub-domain, the
amount of work to its volume.

## Particle-mesh Ewald

The Ewald method splits the Coulomb interaction into a short-ranged part,
summed in real space with the other pair interactions, and a smooth
long-ranged part, summed in reciprocal space. Only the second sum is
replaced. In the smooth particle-mesh method the structure factor
\f$S(\mathbf{k})\f$ is not summed atom by atom but obtained from a mesh:

1. the charges are spread onto a regular mesh with cardinal B-splines of a
   chosen order (every charge contributes to the mesh points around it),
2. the mesh is Fourier transformed with a fast Fourier transform,
3. in reciprocal space it is multiplied by the Ewald influence function,
   corrected for the smoothing of the B-splines,
4. the transform back gives the potential on the mesh, and the forces on the
   atoms follow by differentiating the same B-splines.

The cost is \f$O(N)\f$ for the spreading and interpolation plus
\f$O(M \log M)\f$ for the transforms on \f$M\f$ mesh points. The splitting
parameter \f$\alpha\f$ is that of the Ewald sum of the force field, so the
real-space, self and exclusion terms are unchanged, and the result converges
to the exact Ewald energy and forces as the mesh is refined and the spline
order raised. Because the forces are the exact derivatives of the mesh
energy, the energy is conserved to the accuracy of the integrator.

## The real-space kernel

For plain 12-6 Lennard-Jones pairs with Ewald (or no) electrostatics a
specialised kernel is used. The Ewald real-space term
\f$\mathrm{erfc}(\alpha r)/r\f$ is tabulated as a cubic Hermite spline in
\f$r^2\f$, so the kernel needs neither a square root nor a division nor the
error function per pair, and the force is the exact derivative of the
interpolated energy. Other potentials use the general kernels of the rest of
the code.

## Pressure

For the barostats the pressure tensor is needed every step. It is assembled
by the same threads from the pair, reciprocal, exclusion and tail
contributions and converted to the molecular (centre-of-mass) virial, so the
constant-pressure ensembles do not need a separate \f$O(N^2)\f$ evaluation.

## What the run checks

At the start of the equilibration and of the production stage the energy and
the forces of this engine are compared once with those of the exact code, and
the differences are printed. This is the direct test that the mesh spacing
and the spline order are fine enough. The output also lists the cell grid,
the layout of the sub-domains with their numbers of atoms and pairs, and the
time spent per phase.

## Limits

- Every sub-domain must be at least one cell wide in every direction, so the
  number of threads must fit the box.
- The lists use the minimum-image convention, so
  \f$r_c + r_\text{skin} \le L^\perp_\text{min}/2\f$. The automatic Coulomb
  cutoff is exactly half the box; set the cutoffs explicitly.
- The driver handles molecules in a box. Frameworks, external fields,
  polarization, fractional molecules and the particle-exchange ensembles are
  refused, with a message that points to
  [molecular dynamics](\ref driver_molecular_dynamics).

## References

- R.W. Hockney and J.W. Eastwood, *Computer Simulation Using Particles*,
  Adam Hilger (1988).
- U. Essmann, L. Perera, M.L. Berkowitz, T. Darden, H. Lee and L.G. Pedersen,
  *A smooth particle mesh Ewald method*, J. Chem. Phys. **103**, 8577 (1995).
- S. Plimpton, *Fast Parallel Algorithms for Short-Range Molecular Dynamics*,
  J. Comput. Phys. **117**, 1 (1995).
- L. Verlet, *Computer "Experiments" on Classical Fluids. I*,
  Phys. Rev. **159**, 98 (1967).
- M. Deserno and C. Holm, *How to mesh up Ewald sums. I*,
  J. Chem. Phys. **109**, 7678 (1998).
