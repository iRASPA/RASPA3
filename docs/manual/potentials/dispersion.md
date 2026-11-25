# Dispersion and repulsion
\page potential_dispersion Dispersion and repulsion

The *van der Waals* interaction between two atoms that are not bonded has two
parts. At short distance the electron clouds overlap and repel each other
strongly (Pauli repulsion). At larger distance the atoms attract: the
fluctuating charge distribution of one atom induces a dipole in the other,
and the two correlated dipoles lower the energy. This is the London
dispersion force, and to leading order it falls off as \f$-C_6/r^6\f$. Higher
multipoles add terms in \f$r^{-8}\f$ and \f$r^{-10}\f$.

In a classical force field both parts are represented by an effective pair
potential \f$u(r)\f$ between interaction sites,

\f[
U_\text{vdW} = \sum_{i<j} u_{ij}(r_{ij}) ,
\f]

whose parameters are fitted, usually to the properties of fluids or to
adsorption data. In adsorption this term usually dominates the interaction of
a molecule with a framework.

## Where it is defined

The interactions are given in `force_field.json`. `"SelfInteractions"` gives
the potential between two sites of the same type, and `"BinaryInteractions"`
overrides the potential between two particular types:

```json
"SelfInteractions" : [
  {"name" : "CH4", "type" : "lennard-jones", "parameters" : [158.5, 3.72], "source" : "..."}
],
"BinaryInteractions" : [
  {"names" : ["CH4", "O"], "type" : "lennard-jones", "parameters" : [115.0, 3.47], "source" : "..."}
]
```

Pairs that are not listed are generated from the self interactions with the
[mixing rule](#potential-dispersion-mixing).

## The Lennard-Jones potential

The standard form is the 12-6 Lennard-Jones potential,

\f[
u(r) = 4\epsilon\left[\left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^{6}\right] .
\f]

The \f$r^{-6}\f$ term is the dispersion; the \f$r^{-12}\f$ repulsion has no
physical basis beyond being steep and cheap to compute (the square of the
attractive term). \f$\sigma\f$ is the distance at which the potential is
zero, the minimum lies at \f$r_\text{min} = 2^{1/6}\sigma\f$, and
\f$\epsilon\f$ is its depth.

## Available forms {#potential-dispersion-forms}

`"type"` in the force-field file, the functional form, and the parameters in
the order they are given. Energies are in K and lengths in Å.

| `"type"` | \f$u(r)\f$ | Parameters |
|---|---|---|
| `"lennard-jones"` | \f$4p_0\left[(p_1/r)^{12} - (p_1/r)^{6}\right]\f$ | \f$\epsilon\f$ [K], \f$\sigma\f$ [Å] |
| `"lennard-jones-shifted-force"` | Lennard-Jones with the energy *and* the force shifted to zero at the cutoff | \f$\epsilon\f$, \f$\sigma\f$ |
| `"lennard-jones-second-order-taylor-shifted"` | Lennard-Jones minus its second-order Taylor expansion at the cutoff, so that \f$u\f$, \f$u'\f$ and \f$u''\f$ vanish there | \f$\epsilon\f$, \f$\sigma\f$ |
| `"lennard-jones-switched"` | Lennard-Jones times a quintic switching function between the switching distance and the cutoff | \f$\epsilon\f$, \f$\sigma\f$ |
| `"lennard-jones-force-switched"` | The CHARMM force switch: the force goes smoothly to zero between the switching distance and the cutoff | \f$\epsilon\f$, \f$\sigma\f$ |
| `"wca"` | \f$4p_0\left[(p_1/r)^{12} - (p_1/r)^{6}\right] + p_0\f$ for \f$r < 2^{1/6}p_1\f$, zero beyond | \f$\epsilon\f$, \f$\sigma\f$ |
| `"feynman-hibbs-lennard-jones"` | Lennard-Jones plus the quantum correction \f$\dfrac{\hbar^2}{24\mu k_BT}\,4p_0\,\dfrac{132(p_1/r)^{12} - 30(p_1/r)^{6}}{r^2}\f$ | \f$\epsilon\f$, \f$\sigma\f$, reduced mass \f$\mu\f$ [amu] |
| `"12-6"` | \f$p_0/r^{12} - p_1/r^{6}\f$ | [K Å¹²], [K Å⁶] |
| `"12-6-2-0"` | \f$p_0/r^{12} + p_1/r^{6} + p_2/r^{2} + p_3\f$ | [K Å¹²], [K Å⁶], [K Å²], [K] |
| `"mie"` | \f$p_0/r^{p_1} - p_2/r^{p_3}\f$ | [K Å^p1], repulsive exponent, [K Å^p3], attractive exponent |
| `"cff-9-6"` | \f$p_0/r^{9} - p_1/r^{6}\f$ | [K Å⁹], [K Å⁶] |
| `"cff-eps-sigma"` | \f$p_0\left[2(p_1/r)^{9} - 3(p_1/r)^{6}\right]\f$ | \f$\epsilon\f$, \f$\sigma\f$ |
| `"buckingham"` | \f$p_0\,e^{-p_1 r} - p_2/r^{6}\f$ | \f$A\f$ [K], \f$b\f$ [Å⁻¹], \f$C_6\f$ [K Å⁶] |
| `"born-huggins-meyer"` | \f$p_0\,e^{p_1(p_2 - r)} - p_3/r^{6} - p_4/r^{8}\f$ | \f$A\f$ [K], \f$b\f$ [Å⁻¹], \f$\sigma\f$ [Å], \f$C_6\f$, \f$C_8\f$ |
| `"generic"` | \f$p_0\,e^{-p_1 r} - p_2/r^{4} - p_3/r^{6} - p_4/r^{8} - p_5/r^{10}\f$ | \f$A\f$, \f$b\f$, \f$C_4\f$, \f$C_6\f$, \f$C_8\f$, \f$C_{10}\f$ |
| `"pellenq-nicholson"` | \f$p_0\,e^{-p_1 r} - f_6\,p_2/r^{6} - f_8\,p_3/r^{8} - f_{10}\,p_4/r^{10}\f$, with damping functions \f$f_n\f$ | \f$A\f$, \f$b\f$, \f$C_6\f$, \f$C_8\f$, \f$C_{10}\f$ |
| `"hydrated-ion-water"` | \f$p_0\,e^{-p_1 r} - p_2/r^{4} - p_3/r^{6} - p_4/r^{12}\f$ | \f$A\f$, \f$b\f$, \f$C_4\f$, \f$C_6\f$, \f$C_{12}\f$ |
| `"matsuoka-clementi-yoshimine"` | \f$p_0\,e^{-p_1 r} + p_2\,e^{-p_3 r}\f$ | \f$A\f$ [K], \f$a\f$ [Å⁻¹], \f$B\f$ [K], \f$b\f$ [Å⁻¹] |
| `"morse"` | \f$p_0\left\{\left[1 - e^{-p_1(r - p_2)}\right]^2 - 1\right\}\f$ | \f$D\f$ [K], \f$a\f$ [Å⁻¹], \f$r_0\f$ [Å] |
| `"mm3"` | \f$p_0\left[1.84\times10^{5}\,e^{-12/P} - 2.25\,P^{6}\right]\f$ with \f$P = p_1/r\f$; \f$192.27\,p_0\,P^2\f$ for \f$P > 3.02\f$ | \f$\epsilon\f$, \f$\sigma\f$ |
| `"repulsive-harmonic"` | a purely repulsive harmonic wall | \f$k\f$ [K], \f$r_0\f$ [Å] |
| `"none"` | no interaction | |

Several types have alternative spellings (for example `"lennardjones"`,
`"mcy"`, `"pcff"` and `"class2"` for `"cff-eps-sigma"`).

**Which to use.** The form is part of the force field: parameters fitted with
one form are not valid with another. Some remarks on the families:

- *Exponential repulsion* (Buckingham and the forms derived from it) follows
  the decay of the electron density more closely than \f$r^{-12}\f$. The plain
  Buckingham potential turns over and goes to \f$-\infty\f$ at very short
  distance, which the simulation has to guard against.
- *Softer repulsion* (9-6, Mie) is used in class-II and coarse-grained force
  fields. The Mie potential makes the two exponents parameters.
- *Feynman-Hibbs* adds the leading quantum correction for light molecules
  (hydrogen, helium, neon) at low temperature. It depends on the temperature,
  through the thermal de Broglie wavelength, so the potential changes with the
  state point; methods that reuse energies at other temperatures refuse it.
- *Weeks-Chandler-Andersen* is the repulsive part of Lennard-Jones only, the
  reference fluid of perturbation theory.

## Mixing rules {#potential-dispersion-mixing}

A force field with \f$n\f$ types of sites has \f$n(n+1)/2\f$ pair
interactions. Usually only the \f$n\f$ self interactions are fitted, and the
cross interactions follow from a *mixing rule*, `"MixingRule"`:

| Rule | \f$\epsilon_{ij}\f$ | \f$\sigma_{ij}\f$ |
|---|---|---|
| `"Lorentz-Berthelot"` | \f$\sqrt{\epsilon_i\epsilon_j}\f$ | \f$(\sigma_i + \sigma_j)/2\f$ |
| `"Jorgensen"` | \f$\sqrt{\epsilon_i\epsilon_j}\f$ | \f$\sqrt{\sigma_i\sigma_j}\f$ |
| `"SixthPower"` (Waldman-Hagler) | \f$2\sqrt{\epsilon_i\epsilon_j}\,\dfrac{\sigma_i^3\sigma_j^3}{\sigma_i^6 + \sigma_j^6}\f$ | \f$\left[(\sigma_i^6 + \sigma_j^6)/2\right]^{1/6}\f$ |

The geometric mean of the well depths follows from London's formula for the
dispersion between unlike atoms; the arithmetic mean of the sizes is exact for
hard spheres. Neither is exact for real atoms, and the rule is part of the
force field: TraPPE and most zeolite force fields use Lorentz-Berthelot, OPLS
uses the geometric mean for both.

A rule applies when both sites use the same form with a strength and a size
parameter (Lennard-Jones and its variants, WCA, MM3, the CFF
\f$\epsilon\f$-\f$\sigma\f$ form and Feynman-Hibbs). Any other pair, such as
two different forms or Buckingham, Morse or Mie, has *no* interaction unless
`"BinaryInteractions"` gives one.

A binary interaction always takes precedence over the mixing rule. This is
how force fields that are fitted per pair (the adsorbate-oxygen interactions
of zeolite force fields, for instance) are entered.

## The cutoff {#potential-dispersion-cutoff}

The interaction is evaluated up to a cutoff distance \f$r_c\f$ only, with the
nearest periodic image of each atom. Separate cutoffs can be given for the
interaction with the framework and between molecules. The minimum-image
convention requires the cutoff to be at most half the smallest perpendicular
width of the box; for a framework the number of unit cells is chosen to
satisfy this.

What happens *at* the cutoff is set with `"TruncationMethod"`:

- `"truncated"`: \f$u(r) = 0\f$ for \f$r > r_c\f$. The energy has a small
  jump at the cutoff, so the force is infinite there. This is harmless in
  Monte Carlo, which uses energies only, but unsuitable for molecular
  dynamics.
- `"shifted"`: \f$u(r) - u(r_c)\f$. The energy is continuous; the force still
  has a jump.
- `"switched"`: the potential is multiplied by
  \f$S(x) = 1 - 10x^3 + 15x^4 - 6x^5\f$ with \f$x = (r - r_s)/(r_c - r_s)\f$
  between the switching distance \f$r_s\f$ and the cutoff, so the energy and
  the force go smoothly to zero.
- `"force-switched"`: the *force* is switched smoothly to zero, and the
  potential is its integral.

The cutoff and the truncation are part of the force field. A force field
fitted with a truncated-and-shifted potential at 12 Å gives different
adsorption with a 14 Å cutoff and tail corrections.

## Tail corrections {#potential-dispersion-tail}

Beyond the cutoff the dispersion is not zero. Assuming the fluid is uniform
there, \f$g(r) = 1\f$, the missing energy is

\f[
U_\text{tail} = \frac{1}{2}\sum_{i}\sum_{j} N_i\,\rho_j\; 4\pi\int_{r_c}^{\infty} u_{ij}(r)\,r^2\,dr ,
\f]

which for Lennard-Jones is

\f[
4\pi\int_{r_c}^{\infty} u(r)\,r^2\,dr
= \frac{16\pi}{3}\,\epsilon\,\sigma^3\left[\frac{1}{3}\left(\frac{\sigma}{r_c}\right)^{9} - \left(\frac{\sigma}{r_c}\right)^{3}\right] .
\f]

A corresponding term is added to the pressure. The correction is switched on
with `"TailCorrections"`. It is appropriate for a bulk fluid, where the
density is uniform beyond the cutoff. In a framework the assumption of a
uniform density is poor, and whether to use it follows the convention of the
force field.

## The 1-4 parameters

Some force fields (CHARMM) use different Lennard-Jones parameters for the
1-4 pairs within a molecule. They are given as `"parameters14"` next to the
regular parameters. See
[non-bonded interactions within a molecule](\ref potential_intramolecular_nonbonded).

## Scaling for a fractional molecule

The interaction of a fractional molecule is scaled with its coupling
parameter \f$\lambda\f$. Multiplying the potential by \f$\lambda\f$ is not
enough: the repulsive wall would remain infinite for any \f$\lambda > 0\f$,
and the molecule could never overlap with another. A *soft-core* form is used
instead, in which the distance in the potential is shifted,

\f[
u_\lambda(r) = \lambda\,4\epsilon\left[\frac{1}{\left[\alpha(1-\lambda)^2 + (r/\sigma)^6\right]^2}
- \frac{1}{\alpha(1-\lambda)^2 + (r/\sigma)^6}\right] ,
\f]

with \f$\alpha = 1/2\f$, so that the potential is finite at \f$r = 0\f$ for
\f$\lambda < 1\f$ and becomes the ordinary potential at \f$\lambda = 1\f$.

## Interpolation grids

The interaction of an atom with a rigid framework depends on its position
only, so it can be tabulated on a grid once and interpolated during the
simulation, which is much faster than summing over the framework atoms. See
the grid options on the [Commands](\ref commands) page.

## References

- J.E. Lennard-Jones, *Cohesion*, Proc. Phys. Soc. **43**, 461 (1931).
- F. London, *The general theory of molecular forces*,
  Trans. Faraday Soc. **33**, 8 (1937).
- R.A. Buckingham, *The classical equation of state of gaseous helium, neon
  and argon*, Proc. R. Soc. Lond. A **168**, 264 (1938).
- J.D. Weeks, D. Chandler and H.C. Andersen, *Role of Repulsive Forces in
  Determining the Equilibrium Structure of Simple Liquids*,
  J. Chem. Phys. **54**, 5237 (1971).
- R.P. Feynman and A.R. Hibbs, *Quantum Mechanics and Path Integrals*,
  McGraw-Hill (1965).
- M. Waldman and A.T. Hagler, *New combining rules for rare gas van der Waals
  parameters*, J. Comput. Chem. **14**, 1077 (1993).
- P.J. Steinbach and B.R. Brooks, *New spherical-cutoff methods for
  long-range forces in macromolecular simulation*,
  J. Comput. Chem. **15**, 667 (1994).
- T.C. Beutler, A.E. Mark, R.C. van Schaik, P.R. Gerber and
  W.F. van Gunsteren, *Avoiding singularities and numerical instabilities in
  free energy calculations based on molecular simulations*,
  Chem. Phys. Lett. **222**, 529 (1994).
- D. Dubbeldam, K.S. Walton, T.J.H. Vlugt and S. Calero, *Design,
  Parameterization, and Implementation of Atomic Force Fields for Adsorption
  in Nanoporous Materials*, Adv. Theory Simul. **2**, 1900135 (2019).
