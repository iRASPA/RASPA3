# Monte Carlo
\page driver_monte_carlo Monte Carlo

`"SimulationType" : "MonteCarlo"`

The Monte Carlo driver samples the equilibrium distribution of a system by a
Markov chain of trial moves. It gives thermodynamic averages (energies,
loadings, densities, pressures, heats of adsorption) but no dynamics: the
order of the configurations has no physical meaning.

## What is sampled

In statistical mechanics an average is an integral over all configurations
\f$\mathbf{x}\f$, weighted with the probability density of the ensemble,

\f[
\langle A\rangle = \int A(\mathbf{x})\,\pi(\mathbf{x})\,d\mathbf{x} .
\f]

The weights of the ensembles used in `RASPA3` are, up to normalization,

| Ensemble | Fixed | \f$\pi(\mathbf{x})\f$ proportional to |
|---|---|---|
| canonical | \f$N, V, T\f$ | \f$e^{-\beta U}\f$ |
| isothermal-isobaric | \f$N, p, T\f$ | \f$V^{N} e^{-\beta (U + pV)}\f$ |
| grand-canonical | \f$\mu, V, T\f$ | \f$\dfrac{(\beta f V)^{N}}{N!}\, e^{-\beta U}\f$ |
| Gibbs | \f$N, V, T\f$ of two boxes together | \f$\dfrac{V_A^{N_A} V_B^{N_B}}{N_A!\,N_B!}\, e^{-\beta (U_A + U_B)}\f$ |

Here \f$f\f$ is the fugacity of the reservoir, written this way because the
ideal-gas factors of the chemical potential then cancel against the momentum
integrals. The integral cannot be evaluated on a grid, because almost all
configurations have a vanishing weight. Monte Carlo replaces it by an average
over configurations that are *generated with* probability \f$\pi\f$,

\f[
\langle A\rangle \approx \frac{1}{M}\sum_{i=1}^{M} A(\mathbf{x}_i) .
\f]

## The Markov chain

The configurations are generated one from another. From the old configuration
\f$\mathbf{o}\f$ a trial configuration \f$\mathbf{n}\f$ is proposed with
probability \f$\alpha(\mathbf{o}\to\mathbf{n})\f$ and accepted with
probability \f$\text{acc}(\mathbf{o}\to\mathbf{n})\f$; when it is rejected the
old configuration is counted again. The chain converges to \f$\pi\f$ when it
can reach every configuration and leaves \f$\pi\f$ unchanged. The second
condition is imposed in its strict form, *detailed balance*: the flow from
\f$\mathbf{o}\f$ to \f$\mathbf{n}\f$ equals the reverse flow,

\f[
\pi(\mathbf{o})\,\alpha(\mathbf{o}\to\mathbf{n})\,\text{acc}(\mathbf{o}\to\mathbf{n})
= \pi(\mathbf{n})\,\alpha(\mathbf{n}\to\mathbf{o})\,\text{acc}(\mathbf{n}\to\mathbf{o}) .
\f]

The Metropolis-Hastings choice satisfies it,

\f[
\text{acc}(\mathbf{o}\to\mathbf{n}) = \min\left(1,
\frac{\pi(\mathbf{n})\,\alpha(\mathbf{n}\to\mathbf{o})}{\pi(\mathbf{o})\,\alpha(\mathbf{o}\to\mathbf{n})}\right) .
\f]

For a symmetric proposal, such as a random displacement, this is the familiar
\f$\min(1, e^{-\beta\Delta U})\f$. For a biased proposal the ratio of the
proposal probabilities removes the bias exactly; this is how configurational
bias, smart Monte Carlo and the continuous-fractional-component moves work.
Every move of the [Monte Carlo moves](\ref mc_moves) page is one choice of
\f$\alpha\f$ with the acceptance rule that follows from it.

The ensemble is not given in the input. It follows from the moves that have a
non-zero probability: translation, rotation and reinsertion alone sample the
canonical ensemble, a volume move adds the pressure, a swap move adds the
chemical potential, and the Gibbs moves couple two boxes.

## Cycles and stages

A step is one trial move. The move and the molecule it acts on are drawn at
random, with the probabilities of the input, so that the chain satisfies
detailed balance as a whole and not only move by move. A *cycle* is
\f$\max(20, N)\f$ steps, where \f$N\f$ is the number of molecules, so that on
average every molecule is tried once per cycle. During production the number
of steps per cycle is fixed (the average number of molecules over the second
half of the equilibration): in an open ensemble a cycle whose length followed
the current number of molecules would sample the states of high loading less
often per step, and bias the averages.

The run has up to three stages:

- **Initialization.** Moves are made and nothing is measured. The maximum
  displacements, rotation angles and volume changes are tuned towards an
  acceptance of about one half.
- **Equilibration.** As initialization, and in addition the biasing functions
  of the fractional molecules are built with the Wang-Landau scheme: the bias
  of the visited \f$\lambda\f$ bin is lowered by a factor that is reduced
  whenever the histogram over \f$\lambda\f$ is flat enough. At the end the
  bias is frozen, which is what makes the production run a proper Markov
  chain.
- **Production.** The tuned parameters and the biases are fixed, and the
  properties are sampled once per cycle.

## Averages and their errors

Successive configurations are correlated, so the standard deviation of the
samples underestimates the error of the mean. The production run is therefore
divided into `"NumberOfBlocks"` contiguous blocks; the blocks are long
compared to the correlation time, so their averages \f$A_b\f$ are nearly
independent. The reported value is the mean of the block averages and the
error is the 95% confidence interval of that mean,

\f[
\delta A = t_{0.975,\,n-1}\;\frac{s}{\sqrt{n}}, \qquad
s^2 = \frac{1}{n-1}\sum_{b=1}^{n}\left(A_b - \bar A\right)^2 ,
\f]

with \f$n\f$ the number of blocks and \f$t\f$ the Student factor. Quantities
that are ratios of averages, such as the heat of adsorption from energy and
particle-number fluctuations,

\f[
-\Delta H = k_B T - \frac{\langle U N\rangle - \langle U\rangle\langle N\rangle}
{\langle N^2\rangle - \langle N\rangle^2} ,
\f]

are formed per block from the block sums and then averaged, so their errors
come out of the same procedure.

## Things to check

- **Energy drift.** The running energies are updated with the energy
  difference of every accepted move. At the end they are compared with the
  energy computed from scratch; a difference beyond round-off means that a
  move does not compute the energy difference it should.
- **Equilibration.** The averages are only meaningful when the production run
  starts from an equilibrated system. The evolution of the energy and of the
  number of molecules shows whether it did.
- **Acceptance.** A move that is almost never accepted contributes nothing.
  Insertions in dense systems are the usual case; the configurational-bias
  and fractional-molecule variants exist for them.

## References

- N. Metropolis, A.W. Rosenbluth, M.N. Rosenbluth, A.H. Teller and E. Teller,
  *Equation of State Calculations by Fast Computing Machines*,
  J. Chem. Phys. **21**, 1087 (1953).
- W.K. Hastings, *Monte Carlo sampling methods using Markov chains and their
  applications*, Biometrika **57**, 97 (1970).
- D. Frenkel and B. Smit, *Understanding Molecular Simulation*, 3rd edition,
  Academic Press (2023).
- D. Dubbeldam, A. Torres-Knoop and K.S. Walton, *On the Inner Workings of
  Monte Carlo Codes*, Mol. Simul. **39**, 1253 (2013).
- F. Wang and D.P. Landau, *Efficient, Multiple-Range Random Walk Algorithm to
  Calculate the Density of States*, Phys. Rev. Lett. **86**, 2050 (2001).
- H. Flyvbjerg and H.G. Petersen, *Error estimates on averages of correlated
  data*, J. Chem. Phys. **91**, 461 (1989).
