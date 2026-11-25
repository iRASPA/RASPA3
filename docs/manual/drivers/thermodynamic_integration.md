# Thermodynamic integration
\page driver_thermodynamic_integration Thermodynamic integration

`"SimulationType" : "ThermodynamicIntegration"`

Thermodynamic integration computes a free-energy difference, here the excess
chemical potential of a molecule, as the integral of an ordinary ensemble
average along a path that switches the interactions of that molecule on. This
driver computes one point of the integrand: the average of
\f$\partial U/\partial\lambda\f$ at one fixed value of the coupling parameter.
One run per value of \f$\lambda\f$ gives the whole curve; the
[parallel driver](\ref driver_parallel_thermodynamic_integration) does all
values in one run.

## Why a path is needed

A free energy is not an average of a function of the configuration. It is the
logarithm of a partition function,

\f[
F = -k_B T \ln Q, \qquad Q = \int e^{-\beta U(\mathbf{x})}\,d\mathbf{x} ,
\f]

which a simulation cannot measure directly. Its *derivative* with respect to
a parameter of the energy, however, is an average. Let the energy depend on a
coupling parameter \f$\lambda\f$. Then

\f[
\frac{\partial F}{\partial\lambda}
= -k_B T\,\frac{1}{Q}\frac{\partial Q}{\partial\lambda}
= \frac{\int \dfrac{\partial U}{\partial\lambda}\, e^{-\beta U}\,d\mathbf{x}}{\int e^{-\beta U}\,d\mathbf{x}}
= \left\langle \frac{\partial U}{\partial\lambda}\right\rangle_\lambda ,
\f]

and the free-energy difference between two states is the integral of that
average along any path that connects them (Kirkwood, 1935),

\f[
\Delta F = F(1) - F(0) = \int_0^1 \left\langle \frac{\partial U}{\partial\lambda}\right\rangle_\lambda d\lambda .
\f]

The free energy is a state function, so the result does not depend on the
path. The path may be unphysical; it only has to be reversible, meaning that
the system stays in equilibrium at every \f$\lambda\f$.

## The path: a fractional molecule

The path used here is that of the continuous-fractional-component moves. One
molecule is *fractional*: its interactions with all other molecules and with
the framework are scaled with \f$\lambda\f$. At \f$\lambda = 0\f$ it does not
interact (an ideal-gas molecule in the volume), at \f$\lambda = 1\f$ it is an
ordinary molecule. The free-energy difference between the two ends is exactly
the excess chemical potential,

\f[
\mu^\text{ex} = \int_0^1 \left\langle \frac{\partial U}{\partial\lambda}\right\rangle_\lambda d\lambda ,
\f]

the work of bringing one molecule from the ideal gas into the system.

The scaling is not linear. With a linearly scaled Lennard-Jones potential the
repulsive core at \f$r\to 0\f$ stays infinite for every \f$\lambda > 0\f$,
and \f$\langle\partial U/\partial\lambda\rangle\f$ diverges as
\f$\lambda\to 0\f$. A *soft-core* scaling shifts the distance inside the
potential so that the core has a finite height that grows with
\f$\lambda\f$; the integrand is then finite and smooth over the whole
interval. The charges of the fractional molecule are scaled as well.

Which coupling is pinned follows from the definition of the component: a
single fractional molecule, both molecules of an ion pair with a common
\f$\lambda\f$, or all members of a group. The pair and group couplings keep
the system electrically neutral at every \f$\lambda\f$, which is what makes
the excess chemical potential of a salt computable.

## What the driver does

In the CFCMC moves \f$\lambda\f$ is a variable that performs a random walk
under a Wang-Landau bias. Here it is not: the fractional molecule is created
at the requested \f$\lambda\f$-bin,
\f$\lambda = k/(n_\text{bins}-1)\f$, and stays there. There are no
\f$\lambda\f$-changing moves and no bias. All other moves act as usual, on
the fractional molecule too, so the system equilibrates around it.

During production the derivative of the energy with respect to
\f$\lambda\f$ is evaluated analytically, together with the energy, and
averaged. Its error follows from the block averages.

## Completing the integral

With the averages at all \f$\lambda\f$-bins the integral is a quadrature of a
tabulated function on an equidistant grid. The shape of the integrand
determines how many points are needed: it is usually steepest near
\f$\lambda = 0\f$, where the molecule first excludes volume, and smooth
elsewhere. A calculation should show that the result does not change when
points are added.

The error of the integral follows from the errors of the points, which are
independent because they come from separate runs:
\f$\delta^2(\Delta F) = \sum_k w_k^2\,\delta_k^2\f$ with \f$w_k\f$ the weights
of the quadrature rule.

## Compared with the other routes to a chemical potential

- **Widom insertion** measures
  \f$\mu^\text{ex} = -k_BT\ln\langle e^{-\beta\Delta U}\rangle\f$ from
  attempted insertions of a whole molecule. It fails in dense systems, where
  almost every insertion overlaps.
- **The CFCMC \f$\lambda\f$ histogram** gives the same free energy from the
  probabilities of \f$\lambda = 0\f$ and \f$\lambda = 1\f$ in a walk over
  \f$\lambda\f$. It needs a good bias to visit both ends.
- **Thermodynamic integration** needs neither insertions nor a bias. Every
  point is an ordinary equilibrium average with an ordinary error bar, and
  the points that matter most can be run longest. Its systematic error is the
  quadrature error.

## Things to check

- The system must be equilibrated at each \f$\lambda\f$ separately; the
  environment of a nearly decoupled molecule differs from that of a coupled
  one.
- The curve \f$\langle\partial U/\partial\lambda\rangle(\lambda)\f$ should be
  inspected, not only its integral. A kink or a jump signals a transition
  along the path (the molecule moving to a different site, for instance),
  where the sampling at neighbouring \f$\lambda\f$ values is not connected.

## References

- J.G. Kirkwood, *Statistical Mechanics of Fluid Mixtures*,
  J. Chem. Phys. **3**, 300 (1935).
- D. Frenkel and B. Smit, *Understanding Molecular Simulation*, 3rd edition,
  Academic Press (2023), chapter on free-energy calculations.
- T.C. Beutler, A.E. Mark, R.C. van Schaik, P.R. Gerber and
  W.F. van Gunsteren, *Avoiding singularities and numerical instabilities in
  free energy calculations based on molecular simulations*,
  Chem. Phys. Lett. **222**, 529 (1994).
- W. Shi and E.J. Maginn, *Continuous Fractional Component Monte Carlo: An
  Adaptive Biasing Method for Open System Atomistic Simulations*,
  J. Chem. Theory Comput. **3**, 1451 (2007).
- A. Rahbari, R. Hens, M. Ramdin, O.A. Moultos, D. Dubbeldam and
  T.J.H. Vlugt, *Recent advances in the continuous fractional component Monte
  Carlo methodology*, Mol. Simul. **47**, 804 (2021).
