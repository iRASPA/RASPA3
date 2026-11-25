# Bends
\page potential_bends Bends

A bend term describes the energy of changing the angle \f$\theta\f$ between
two bonds that share an atom. For atoms A-B-C, with B the central atom,

\f[
\cos\theta = \frac{\mathbf{r}_{BA}\cdot\mathbf{r}_{BC}}{r_{BA}\,r_{BC}} ,
\qquad
U_\text{bend} = \sum_{\text{angles}} u(\theta) .
\f]

## The physics

The preferred angle reflects the hybridization of the central atom: close to
109.5° for a tetrahedral \f$sp^3\f$ carbon, 120° for a planar \f$sp^2\f$
centre, 180° for a linear \f$sp\f$ centre. Bending is softer than stretching:
typical force constants are of the order of \f$5\times10^4\f$ K/rad², so an
angle fluctuates by several degrees at room temperature, and bending
vibrations have periods of 20-50 fs. Unlike a bond length, the angle is a
coordinate that a simulation of a flexible molecule has to sample.

## In the molecule file

Each entry of `"Bends"` names the three pseudo-atom types (the central atom
in the middle), the type of potential and its parameters:

```json
"Bends" : [
  [["CH3", "CH2", "CH2"], "HARMONIC", [62500.0, 114.0]]
]
```

Reference angles are given in degrees. Force constants are per radian
squared.

## Available forms

| Type | \f$u(\theta)\f$ | Parameters |
|---|---|---|
| `"HARMONIC"` | \f$\frac{1}{2}p_0\,(\theta - p_1)^2\f$ | \f$k\f$ [K/rad²], \f$\theta_0\f$ [°] |
| `"HARMONIC_COSINE"` | \f$\frac{1}{2}p_0\,(\cos\theta - \cos p_1)^2\f$ | \f$k\f$ [K], \f$\theta_0\f$ [°] |
| `"COSINE"` | \f$p_0\left[1 + \cos(p_1\theta - p_2)\right]\f$ | \f$k\f$ [K], multiplicity, phase [°] |
| `"QUARTIC"` | \f$\frac{1}{2}p_0(\theta-p_1)^2 + \frac{1}{3}p_2(\theta-p_1)^3 + \frac{1}{4}p_3(\theta-p_1)^4\f$ | [K/rad²], \f$\theta_0\f$ [°], [K/rad³], [K/rad⁴] |
| `"CFF_QUARTIC"` | \f$p_0(\theta-p_1)^2 + p_2(\theta-p_1)^3 + p_3(\theta-p_1)^4\f$ | [K/rad²], \f$\theta_0\f$ [°], [K/rad³], [K/rad⁴] |
| `"MM3"` | \f$p_0\,\Delta^2\left[1 - 0.014\Delta + 5.6\times10^{-5}\Delta^2 - 7\times10^{-7}\Delta^3 + 2.2\times10^{-8}\Delta^4\right]\f$, \f$\Delta = \theta - p_1\f$ | [mdyne Å/rad²], \f$\theta_0\f$ [°] |
| `"TAFIPOLSKY"` | \f$\frac{1}{2}p_0\,(1 + \cos\theta)(1 + \cos 2\theta)\f$ | \f$k\f$ [K] |
| `"CORE_SHELL"` | \f$\frac{1}{2}p_0\,(\theta - p_1)^2\f$ | \f$k\f$ [K/rad²], \f$\theta_0\f$ [°] |
| `"FIXED"` | a constraint: the angle is held at its value | |
| `"RIGID"` | the angle is held at the given value | \f$\theta_0\f$ [°] |

Notes on the forms:

- **Harmonic.** The standard form, with the factor \f$\frac{1}{2}\f$ included;
  force fields that write \f$k(\theta-\theta_0)^2\f$ quote half the value.
- **Harmonic cosine.** Harmonic in \f$\cos\theta\f$ instead of
  \f$\theta\f$. Near the minimum the two are related by
  \f$k_\theta = k_{\cos}\sin^2\theta_0\f$. The cosine form has no
  singularity in its derivative at \f$\theta = 180^\circ\f$, where
  \f$d\theta/d\cos\theta\f$ diverges, which makes it the better choice for
  angles that can become linear. It is the form used by DREIDING.
- **Cosine.** A periodic form for centres with several equivalent minima, as
  in octahedral or square-planar coordination around a metal, where the angle
  can be 90° or 180°. This is the bend of UFF.
- **Quartic forms and MM3.** Anharmonic corrections for accurate vibrational
  frequencies and for large distortions, as in class-II force fields.
- **Tafipolsky.** A form for the coordination geometry of metal centres in
  metal-organic frameworks.

## Sampling an angle

In configurational-bias growth a new atom is placed on a cone around the
previous bond. The angle is drawn from

\f[
P(\theta) \propto \sin\theta\; e^{-\beta u(\theta)} ,
\f]

where \f$\sin\theta\f$ is the solid-angle element: there are more directions
at 90° from an axis than close to it. A stiff harmonic bend is sampled from
its Gaussian Boltzmann factor directly. The factor \f$\sin\theta\f$ is also
why a linear angle is never exactly linear on average: the weight of
\f$\theta = 180^\circ\f$ itself is zero.

When two new atoms are attached to the same atom (a branch), their bends and
the bend between them are not independent, and the two positions have to be
generated together.

## Fixed angles

A `"FIXED"` bend between two fixed bonds makes the three atoms a rigid unit.
In molecular dynamics it is imposed as a constraint on the 1-3 distance.

## References

- A.K. Rappé, C.J. Casewit, K.S. Colwell, W.A. Goddard III and W.M. Skiff,
  *UFF, a full periodic table force field for molecular mechanics and
  molecular dynamics simulations*, J. Am. Chem. Soc. **114**, 10024 (1992).
- S.L. Mayo, B.D. Olafson and W.A. Goddard III, *DREIDING: a generic force
  field for molecular simulations*, J. Phys. Chem. **94**, 8897 (1990).
- N.L. Allinger, Y.H. Yuh and J.-H. Lii, *Molecular mechanics. The MM3 force
  field for hydrocarbons. 1*, J. Am. Chem. Soc. **111**, 8551 (1989).
- M. Tafipolsky, S. Amirjalayer and R. Schmid, *Ab initio parametrized MM3
  force field for the metal-organic framework MOF-5*,
  J. Comput. Chem. **28**, 1169 (2007).
- M.G. Martin and J.I. Siepmann, *Transferable Potentials for Phase
  Equilibria. 1. United-Atom Description of n-Alkanes*,
  J. Phys. Chem. B **102**, 2569 (1998).
