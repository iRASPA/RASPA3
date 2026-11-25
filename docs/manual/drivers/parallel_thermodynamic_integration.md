# Parallel thermodynamic integration
\page driver_parallel_thermodynamic_integration Parallel thermodynamic integration

`"SimulationType" : "ParallelThermodynamicIntegration"`

This driver computes the whole curve
\f$\langle\partial U/\partial\lambda\rangle(\lambda)\f$ of
[thermodynamic integration](\ref driver_thermodynamic_integration), and its
integral, the excess chemical potential, in one run. One replica is made per
\f$\lambda\f$-bin, every replica runs in its own thread, and the replicas
exchange their values of \f$\lambda\f$. The theory of the integral itself is
on the page of the serial driver; this page explains the exchange and the
quadrature.

## Replicas along the coupling parameter

Replica \f$k\f$ holds a fractional molecule (or pair, or group) pinned at
\f$\lambda_k = k/(n_\text{bins}-1)\f$ and samples the canonical distribution
of the Hamiltonian \f$U(\mathbf{x};\lambda_k)\f$. All replicas have the same
temperature; they differ in the Hamiltonian. The set of replicas has the
product weight

\f[
\pi(\mathbf{x}_0,\dots,\mathbf{x}_{n-1}) \propto \prod_{k} e^{-\beta\,U(\mathbf{x}_k;\,\lambda_k)} .
\f]

## Hamiltonian replica exchange

A separate simulation at each \f$\lambda\f$ has a weakness: the environment
of a nearly decoupled molecule relaxes slowly at large \f$\lambda\f$, and a
molecule trapped at one site cannot reach another. The exchange move removes
it. Two replicas \f$A\f$ and \f$B\f$ at neighbouring bins swap their
\f$\lambda\f$ values: configuration \f$\mathbf{x}_A\f$ continues at
\f$\lambda_B\f$ and \f$\mathbf{x}_B\f$ at \f$\lambda_A\f$. The ratio of the
weights is

\f[
\frac{\pi(\text{new})}{\pi(\text{old})}
= \frac{e^{-\beta U(\mathbf{x}_A;\lambda_B)}\; e^{-\beta U(\mathbf{x}_B;\lambda_A)}}
{e^{-\beta U(\mathbf{x}_A;\lambda_A)}\; e^{-\beta U(\mathbf{x}_B;\lambda_B)}} ,
\f]

so with \f$\Delta U_A = U(\mathbf{x}_A;\lambda_B) - U(\mathbf{x}_A;\lambda_A)\f$
and \f$\Delta U_B = U(\mathbf{x}_B;\lambda_A) - U(\mathbf{x}_B;\lambda_B)\f$
the acceptance rule is

\f[
\text{acc} = \min\left(1,\; e^{-\beta(\Delta U_A + \Delta U_B)}\right) .
\f]

In contrast to temperature exchange, the energies after the swap are not
known in advance: each configuration must be evaluated at the other
\f$\lambda\f$. Only the interactions of the fractional molecule depend on
\f$\lambda\f$, so this costs the energy of one molecule, not of the system.

Through the exchanges a configuration performs a random walk in
\f$\lambda\f$. At small \f$\lambda\f$ the fractional molecule passes through
other molecules and through barriers; when the configuration returns to large
\f$\lambda\f$ it arrives with the molecule somewhere else. This is the same
mechanism by which the CFCMC moves equilibrate dense systems, obtained here
without a bias.

Because the exchanges only permute the \f$\lambda\f$ values over the
replicas, every \f$\lambda\f$-bin is occupied by exactly one replica at all
times, and every point of the curve is sampled in every cycle.

## Bookkeeping: averages belong to a \f$\lambda\f$

An ensemble average is a property of a Hamiltonian, not of a thread. When two
replicas exchange \f$\lambda\f$, the accumulated averages (energies,
pressure, loadings, histograms, and \f$\partial U/\partial\lambda\f$) go with
the \f$\lambda\f$. What a replica reports at the end are therefore the
averages of the \f$\lambda\f$-bin it holds, sampled by all the configurations
that passed through that bin.

The number of Monte Carlo steps per production cycle is fixed at the start of
production, from the mean number of molecules of the replicas, so that all
bins are sampled at the same rate.

## The integral

At the end the averages of all bins are joined into one curve on an
equidistant grid and integrated from 0 to 1. Three rules are reported:

- a **natural cubic spline** through the points, integrated exactly: the
  recommended value;
- **composite Simpson** (the 1/3 rule, with a 3/8 segment at the end when the
  number of intervals is odd);
- the **trapezoidal rule**, as a baseline.

Gaussian quadrature, which would need fewer points for the same accuracy,
does not apply: it requires the integrand at its own non-equidistant nodes,
and the data are on fixed bins.

The comparison of the three is the test of the *discretization*. When they
agree to much better than the statistical error, the \f$\lambda\f$ grid
resolves the curve and the error of the result is the sampling error. When
they differ by more than the sampling error, the curve has more curvature
than the grid can represent, and more bins are needed; running longer would
not help.

The *statistical* error is obtained by integrating each block of the
production run separately with the same rule: the spread of the block-wise
integrals gives the confidence interval. This is correct even though the bins
are correlated by the exchanges, because the correlation is contained within
each block.

During the run the current curve and the running spline integral are written
to a file that is overwritten at regular intervals, so the convergence can be
followed.

## Things to check

- **Exchange acceptance per pair.** A pair of bins with a low acceptance
  splits the ladder into two parts that do not communicate. The overlap is
  usually poorest near \f$\lambda = 0\f$; more bins help there.
- **The shape of the curve.** It should be smooth. A feature that is one bin
  wide is not resolved.
- **The three quadratures**, as described above.

## Limits

- Exactly one system, with exactly one component marked for thermodynamic
  integration.
- The number of threads is the number of \f$\lambda\f$-bins.
- Binary restart files are not supported by this driver.

## References

- J.G. Kirkwood, *Statistical Mechanics of Fluid Mixtures*,
  J. Chem. Phys. **3**, 300 (1935).
- H. Fukunishi, O. Watanabe and S. Takada, *On the Hamiltonian replica
  exchange method for efficient sampling of biomolecular systems: Application
  to protein structure prediction*, J. Chem. Phys. **116**, 9058 (2002).
- Y. Sugita, A. Kitao and Y. Okamoto, *Multidimensional replica-exchange
  method for free-energy calculations*, J. Chem. Phys. **113**, 6042 (2000).
- A. Rahbari, R. Hens, M. Ramdin, O.A. Moultos, D. Dubbeldam and
  T.J.H. Vlugt, *Recent advances in the continuous fractional component Monte
  Carlo methodology*, Mol. Simul. **47**, 804 (2021).
