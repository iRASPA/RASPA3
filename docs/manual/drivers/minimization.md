# Minimization
\page driver_minimization Minimization

`"SimulationType" : "Minimization"`

The minimization driver moves the initial configuration to the nearest local
minimum of the potential energy. At the minimum it can compute the quantities
that follow from the second derivatives of the energy: the elastic constants,
the normal modes and the phonon dispersion. A minimum is a zero-temperature
structure; nothing is sampled.

## The problem

Let \f$\mathbf{x}\f$ collect all degrees of freedom: the positions of the
atoms of flexible molecules, the centres of mass and orientations of rigid
molecules, and, when the cell is allowed to change, the strain of the cell.
Near a point \f$\mathbf{x}\f$ the energy is

\f[
U(\mathbf{x} + \Delta\mathbf{x}) \approx U(\mathbf{x})
+ \mathbf{g}^T \Delta\mathbf{x}
+ \frac{1}{2}\,\Delta\mathbf{x}^T \mathbf{H}\,\Delta\mathbf{x} ,
\f]

with the gradient \f$\mathbf{g}\f$ and the Hessian \f$\mathbf{H}\f$ in these
generalized coordinates. A minimum has \f$\mathbf{g} = 0\f$ and a positive
definite Hessian (apart from the zero modes of overall translation).

The Newton-Raphson step \f$\Delta\mathbf{x} = -\mathbf{H}^{-1}\mathbf{g}\f$
goes to the stationary point of the quadratic model. It converges
quadratically close to a minimum, but far from one the Hessian has negative
eigenvalues and the step walks uphill along them, towards a saddle point.

## Eigenvector following

The driver uses Baker's eigenvector-following method, which keeps the
second-order convergence and repairs the step. The Hessian is diagonalized,
\f$\mathbf{H}\mathbf{v}_i = b_i \mathbf{v}_i\f$, and the gradient is expressed
in the eigenvectors, \f$F_i = \mathbf{v}_i^T\mathbf{g}\f$. The step is a
shifted Newton step,

\f[
\Delta\mathbf{x} = -\sum_i \frac{F_i}{b_i - \lambda}\,\mathbf{v}_i ,
\f]

where the shift \f$\lambda\f$ comes from the *rational function
optimization* model: it is the lowest eigenvalue of the Hessian augmented
with the gradient, i.e. the lowest solution of

\f[
\lambda = \sum_i \frac{F_i^2}{\lambda - b_i} .
\f]

That solution is negative and lies below the lowest eigenvalue, so every
denominator \f$b_i - \lambda\f$ is positive and the step goes downhill along
*every* mode, including those with negative curvature. As the minimum is
approached \f$\lambda \to 0\f$ and the step becomes the Newton step. Modes
with an eigenvalue below a small threshold (the translations, and soft modes)
are left out, and the step is scaled down to a maximum length when the
quadratic model cannot be trusted that far.

Because the full Hessian is computed analytically and diagonalized every
step, the cost per step grows with the cube of the number of degrees of
freedom. In return the method converges in few steps, also for the soft,
strongly coupled modes of framework materials where first-derivative methods
stall, and it *knows* that the point it ends at is a minimum: all eigenvalues
are positive.

## Variable cell

With a variable cell the strain is part of \f$\mathbf{x}\f$, and its gradient
is the difference between the internal stress and the applied pressure. The
cell degrees of freedom are taken logarithmically, so that a step multiplies
the cell and cannot make a box length negative. The minimized quantity is
then the enthalpy \f$U + pV\f$.

## Elastic constants

At the minimum the energy is expanded in a homogeneous strain
\f$\boldsymbol{\epsilon}\f$ and the internal displacements
\f$\mathbf{u}\f$ of the atoms,

\f[
\Delta U = \frac{1}{2}\,\mathbf{u}^T \mathbf{H}_{uu}\,\mathbf{u}
+ \mathbf{u}^T \mathbf{H}_{u\epsilon}\,\boldsymbol{\epsilon}
+ \frac{1}{2}\,\boldsymbol{\epsilon}^T \mathbf{H}_{\epsilon\epsilon}\,\boldsymbol{\epsilon} .
\f]

When the cell is strained the atoms relax to
\f$\mathbf{u} = -\mathbf{H}_{uu}^{-1}\mathbf{H}_{u\epsilon}\boldsymbol{\epsilon}\f$.
Eliminating them gives the relaxed elastic tensor,

\f[
\mathbf{C} = \frac{1}{V}\left[\mathbf{H}_{\epsilon\epsilon}
- \mathbf{H}_{\epsilon u}\,\mathbf{H}_{uu}^{-1}\,\mathbf{H}_{u\epsilon}\right]
+ \text{pressure terms} .
\f]

The first term is the *Born* term, the response to an affine deformation in
which every atom follows the strain. The second is the *relaxation* term; it
is always a softening, and for flexible frameworks it is large. The inverse
of \f$\mathbf{H}_{uu}\f$ is taken on the space without the translational zero
modes. From the \f$6\times6\f$ stiffness matrix in Voigt notation follow the
compliance matrix \f$\mathbf{S} = \mathbf{C}^{-1}\f$, the Born stability
criterion (all eigenvalues of \f$\mathbf{C}\f$ positive) and the bulk, shear
and Young's moduli in the Voigt (uniform strain), Reuss (uniform stress) and
Hill (their mean) averages for a polycrystal.

## Normal modes

Small oscillations about the minimum follow
\f$\mathbf{M}\ddot{\mathbf{x}} = -\mathbf{H}\mathbf{x}\f$. With the
mass-weighted Hessian
\f$\tilde{\mathbf{H}} = \mathbf{M}^{-1/2}\mathbf{H}\mathbf{M}^{-1/2}\f$ this is
an eigenvalue problem,

\f[
\tilde{\mathbf{H}}\,\mathbf{e}_k = \omega_k^2\,\mathbf{e}_k ,
\f]

whose eigenvalues are the squared frequencies and whose eigenvectors are the
displacement patterns of the modes. For rigid molecules the mass matrix holds
the molecular mass for the centre of mass and the inertia tensor for the
orientation. A negative eigenvalue is an imaginary frequency: the structure
is not a minimum along that mode.

## Phonon dispersion

A periodic crystal has modes with a wave vector \f$\mathbf{k}\f$, in which the
displacement in the cell at lattice vector \f$\mathbf{R}\f$ carries the phase
\f$e^{i\mathbf{k}\cdot\mathbf{R}}\f$. They are the eigenvectors of the
*dynamical matrix*, the Fourier transform of the force constants between an
atom and the periodic images of the others,

\f[
D_{ij}(\mathbf{k}) = \frac{1}{\sqrt{m_i m_j}}
\sum_{\mathbf{R}} \Phi_{ij}(\mathbf{R})\, e^{i\mathbf{k}\cdot\mathbf{R}} .
\f]

Diagonalizing \f$\mathbf{D}(\mathbf{k})\f$ along a path between
high-symmetry points gives the band structure \f$\omega_k(\mathbf{k})\f$. At
\f$\mathbf{k} = 0\f$ it reduces to the normal-mode problem. For charged
systems the reciprocal-space part of the Ewald sum contributes its own
\f$\mathbf{k}\f$-dependent term.

## Things to check

- The run must end with a gradient below the tolerance *and* without negative
  eigenvalues; otherwise the elastic constants and frequencies are not those
  of a stable structure.
- A minimum is local. Different starting configurations of a flexible
  framework, or of adsorbed molecules, can end in different minima.
- Potentials must be twice differentiable at the cutoff for the Hessian to be
  meaningful; a truncated potential gives discontinuous forces.

## References

- J. Baker, *An algorithm for the location of transition states*,
  J. Comput. Chem. **7**, 385 (1986).
- A. Banerjee, N. Adams, J. Simons and R. Shepard, *Search for stationary
  points on surfaces*, J. Phys. Chem. **89**, 52 (1985).
- D. Dubbeldam, R. Krishna and R.Q. Snurr, *Method for Analyzing Structural
  Changes of Flexible Metal-Organic Frameworks Induced by Adsorbates*,
  J. Phys. Chem. C **113**, 19317 (2009).
- M. Born and K. Huang, *Dynamical Theory of Crystal Lattices*, Oxford
  University Press (1954).
- J.F. Lutsko, *Generalized expressions for the calculation of elastic
  constants by computer simulation*, J. Appl. Phys. **65**, 2991 (1989).
- R. Hill, *The Elastic Behaviour of a Crystalline Aggregate*,
  Proc. Phys. Soc. A **65**, 349 (1952).
