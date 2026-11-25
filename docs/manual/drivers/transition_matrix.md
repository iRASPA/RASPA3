# Transition-matrix Monte Carlo
\page driver_transition_matrix Transition-matrix Monte Carlo

`"SimulationType" : "MonteCarloTransitionMatrix"`

Transition-matrix Monte Carlo (TMMC) computes the probability \f$\Pi(N)\f$ of
finding \f$N\f$ molecules in the system, over a whole range of \f$N\f$, in a
single grand-canonical simulation. From \f$\Pi(N)\f$ follow the complete
adsorption isotherm at that temperature, the free energy as a function of
loading, and, where the distribution has two peaks, the conditions of phase
coexistence and the two branches of a hysteresis loop.

## Macrostates

The number of molecules \f$N\f$ of the adsorbate is the *macrostate*. In the
grand-canonical ensemble its probability is

\f[
\Pi(N;\,\beta, f) = \frac{(\beta f)^N\, Q(N, V, T)}{\Xi(\beta, f)} ,
\f]

with \f$Q\f$ the canonical partition function (including \f$V^N/N!\f$). An
ordinary simulation measures \f$\Pi(N)\f$ as a histogram of visits, and so
only knows it where \f$N\f$ is visited often: a few molecules around the
average. TMMC obtains it from a different observable, which is available for
*every* attempted move.

## The collection matrix

Every insertion or deletion attempt from a configuration with \f$N\f$
molecules has an acceptance probability
\f$a = \min(1, \dots)\f$ that the simulation computes anyway. Instead of only
using it to decide the move, TMMC records it. For an attempted insertion,

\f[
C(N \to N+1) \mathrel{+}= a, \qquad C(N \to N) \mathrel{+}= 1 - a ,
\f]

and likewise for a deletion with \f$N-1\f$. The recorded value is the *unbiased* acceptance probability,
whether the move was then accepted or not, and whatever bias was used to make
that decision. The matrix \f$C\f$ is tridiagonal, because a move changes
\f$N\f$ by at most one.

Normalizing each row gives an estimate of the probability that a step from
macrostate \f$N\f$ leads to macrostate \f$M\f$,

\f[
P(N \to M) = \frac{C(N \to M)}{\sum_{K} C(N \to K)} .
\f]

Recording the probability instead of the outcome is what makes the method
efficient: an insertion that is accepted once in a million attempts
contributes \f$10^{-6}\f$ at every attempt, instead of one count after a
million.

## From transition probabilities to the distribution

The Markov chain satisfies detailed balance between macrostates,

\f[
\Pi(N)\,P(N \to N+1) = \Pi(N+1)\,P(N+1 \to N) ,
\f]

so the distribution follows by walking up the chain of macrostates,

\f[
\ln\Pi(N+1) = \ln\Pi(N) + \ln\frac{P(N \to N+1)}{P(N+1 \to N)} .
\f]

Each link needs the transitions between two neighbouring macrostates only.
Errors in distant parts of the range do not affect each other, and the
relative accuracy of \f$\Pi\f$ is uniform over many orders of magnitude,
which a visit histogram cannot achieve.

## The bias

For the recursion every macrostate of the range must be visited, including
those that are improbable at the simulated fugacity. The moves that change
\f$N\f$ are therefore accepted with a *biased* rule,

\f[
\text{acc}_\text{biased}(N \to M) = \min\left(1,\;
\frac{\pi(\mathbf{n})\,\alpha(\mathbf{n}\to\mathbf{o})}{\pi(\mathbf{o})\,\alpha(\mathbf{o}\to\mathbf{n})}\;
e^{\,\eta(M) - \eta(N)}\right), \qquad \eta(N) = -\ln\Pi(N) ,
\f]

which makes all macrostates equally probable: a flat histogram. The walk then
diffuses freely over the whole range, across free-energy barriers of any
height.

The bias needs \f$\Pi\f$, which is what the simulation is computing. The two
are obtained together: every `"TMMCUpdateEvery"` steps the bias is rebuilt
from the current collection matrix. This is safe because the collection
matrix stores unbiased probabilities, so its entries stay valid when the bias
changes; a poor early bias only costs efficiency, never correctness. That is
the essential difference with the Wang-Landau method, in which the bias *is*
the estimate and the changing bias breaks detailed balance until its updates
have become negligible.

At the very start the collection matrix is empty and the bias is flat over
the states that were never visited. A Wang-Landau bias on the macrostate can
be switched on for the initialization and equilibration stages to drive the
walk through the range for the first time; in production the
transition-matrix bias takes over.

## Reweighting

The distribution is computed at the fugacity \f$f_\text{ref}\f$ of the
simulation, but it is known at every other fugacity exactly, because
\f$\Pi(N)\propto(\beta f)^N Q(N)\f$,

\f[
\ln\Pi(N;\,f) = \ln\Pi(N;\,f_\text{ref}) + N\ln\frac{f}{f_\text{ref}} + \text{constant} .
\f]

The average loading at any pressure follows as
\f$\langle N\rangle = \sum_N N\,\Pi(N;f)\f$, which gives the whole isotherm
at this temperature from one simulation. Reweighting in temperature is not
possible: the energy is not recorded.

## What the distribution tells

- **Free energy.** \f$-k_BT\ln\Pi(N)\f$ is the free energy of the system as a
  function of the number of molecules, at that chemical potential.
- **Phase coexistence.** A first-order transition (vapour-liquid in a box,
  capillary condensation in a pore) shows as two peaks. At coexistence the
  two peaks have equal area. The barrier between them is the free energy of
  the interface.
- **Metastable branches.** Averaging \f$N\f$ over the low-density peak alone
  gives the adsorption branch of a hysteresis loop, and over the high-density
  peak alone the desorption branch.

## Fractional molecules

When the insertion and deletion are done with a fractional molecule over
several \f$\lambda\f$ bins, the macrostate is the pair \f$(N, \lambda)\f$ and
the chain of macrostates runs through all \f$\lambda\f$ bins between
consecutive integer states. The recursion is unchanged; the reported
\f$\Pi(N)\f$ is that of the integer states.

## Things to check

- **Coverage.** Every macrostate of the range must be visited, in both
  directions. The statistics file numbers the *segments*: runs of states
  joined by links with accepted transitions both ways. Across a link without
  them \f$\ln\Pi\f$ is carried flat to keep the bias continuous; that step is
  not an estimate. Differences of \f$\ln\Pi\f$ mean something only within one
  segment.
- **The upper bound.** When the average loading approaches the maximum of the
  macrostate range, the isotherm saturates artificially. The range must
  extend beyond the highest loading of interest.
- **Round trips.** As with any flat-histogram method, the walk must travel
  from one end of the range to the other many times.

## References

- M. Fitzgerald, R.R. Picard and R.N. Silver, *Canonical transition
  probabilities for adaptive Metropolis simulation*,
  Europhys. Lett. **46**, 282 (1999).
- J.R. Errington, *Direct calculation of liquid-vapor phase equilibria from
  transition matrix Monte Carlo simulation*, J. Chem. Phys. **118**, 9915
  (2003).
- V.K. Shen and J.R. Errington, *Metastability and Instability in the
  Lennard-Jones Fluid Investigated by Transition-Matrix Monte Carlo*,
  J. Phys. Chem. B **108**, 19595 (2004).
- D.W. Siderius and V.K. Shen, *Use of the Grand Canonical Transition-Matrix
  Monte Carlo Method to Model Gas Adsorption in Porous Materials*,
  J. Phys. Chem. C **117**, 5861 (2013).
- F. Wang and D.P. Landau, *Efficient, Multiple-Range Random Walk Algorithm to
  Calculate the Density of States*, Phys. Rev. Lett. **86**, 2050 (2001).
