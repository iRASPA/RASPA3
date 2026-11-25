# Bonds
\page potential_bonds Bonds

A chemical bond holds two atoms at a preferred distance. Stretching or
compressing it costs energy, and the bond term describes that cost as a
function of the distance \f$r\f$ between the two atoms,

\f[
U_\text{bond} = \sum_{\text{bonds}} u(r) .
\f]

## The physics

The true energy of a bond as a function of its length is an asymmetric well:
steep on compression, where the nuclei repel, and flattening out on
extension, towards the dissociation energy. Near the minimum any well is a
parabola. At ordinary temperatures a bond explores only a few hundredths of an
ångström around its equilibrium length, so the harmonic approximation,

\f[
u(r) = \frac{1}{2}\,k\,(r - r_0)^2 ,
\f]

is adequate for structure and thermodynamics. The force constant is large:
typical values of \f$k/k_B\f$ are of the order of \f$10^5\f$ K/Å², which
gives vibrations with periods of about 10 fs. The anharmonic forms matter
when the vibrational frequencies themselves are of interest, when bonds are
strongly stretched, or when a force field couples the bond to other
coordinates.

## In the molecule file

Each entry of `"Bonds"` names the two pseudo-atom types, the type of
potential and its parameters:

```json
"Bonds" : [
  [["CH3", "CH2"], "HARMONIC", [96500.0, 1.54]]
]
```

The entry applies to every bond of the `"Connectivity"` between atoms of
those two types.

## Available forms

| Type | \f$u(r)\f$ | Parameters |
|---|---|---|
| `"HARMONIC"` | \f$\frac{1}{2}p_0\,(r - p_1)^2\f$ | \f$k\f$ [K/Å²], \f$r_0\f$ [Å] |
| `"MORSE"` | \f$p_0\left\{\left[1 - e^{-p_1(r - p_2)}\right]^2 - 1\right\}\f$ | \f$D\f$ [K], \f$a\f$ [Å⁻¹], \f$r_0\f$ [Å] |
| `"QUARTIC"` | \f$\frac{1}{2}p_0(r-p_1)^2 + \frac{1}{3}p_2(r-p_1)^3 + \frac{1}{4}p_3(r-p_1)^4\f$ | [K/Å²], \f$r_0\f$ [Å], [K/Å³], [K/Å⁴] |
| `"CFF_QUARTIC"` | \f$p_0(r-p_1)^2 + p_2(r-p_1)^3 + p_3(r-p_1)^4\f$ | [K/Å²], \f$r_0\f$ [Å], [K/Å³], [K/Å⁴] |
| `"MM3"` | \f$p_0(r-p_1)^2\left[1 - 2.55(r-p_1) + \frac{7}{12}\,2.55^2\,(r-p_1)^2\right]\f$ | [mdyne/Å per molecule], \f$r_0\f$ [Å] |
| `"RESTRAINED_HARMONIC"` | \f$\frac{1}{2}p_0(r-p_1)^2\f$ for \f$\lvert r-p_1\rvert \le p_2\f$; linear beyond: \f$\frac{1}{2}p_0p_2^2 + p_0p_2(\lvert r-p_1\rvert - p_2)\f$ | [K/Å²], \f$r_0\f$ [Å], [Å] |
| `"CORE_SHELL_SPRING"` | \f$\frac{1}{2}p_0\,r^2\f$ | [K/Å²] |
| `"LENNARD_JONES"` | \f$4p_0\left[(p_1/r)^{12} - (p_1/r)^{6}\right]\f$ | \f$\epsilon\f$ [K], \f$\sigma\f$ [Å] |
| `"LJ_12_6"` | \f$p_0/r^{12} - p_1/r^{6}\f$ | [K Å¹²], [K Å⁶] |
| `"BUCKINGHAM"` | \f$p_0\,e^{-p_1 r} - p_2/r^{6}\f$ | [K], [Å⁻¹], [K Å⁶] |
| `"FIXED"` | a constraint: the distance is held at its value | |
| `"NONE"` | no bond term | |

Notes on the forms:

- **Harmonic.** Note the factor \f$\frac{1}{2}\f$: force fields that write
  \f$k(r-r_0)^2\f$ (AMBER, CHARMM, CFF) quote a force constant half as large.
- **Morse.** \f$D\f$ is the depth of the well and \f$a\f$ its inverse width;
  the harmonic force constant at the minimum is \f$2Da^2\f$. The energy tends
  to zero at large distance, so the bond can dissociate.
- **Quartic forms.** The cubic term makes the well asymmetric, as a real bond
  is; the quartic term keeps the energy from going to \f$-\infty\f$ on
  extension. The two forms differ in their prefactors only, following the
  conventions of the force fields they come from.
- **Restrained harmonic.** Harmonic near the minimum and linear beyond a
  distance, so that the force is bounded. It is meant for restraints.
- **Core-shell spring.** The spring between the core and the shell of a
  polarizable ion in the shell model; its rest length is zero.
- **Non-bonded forms as a bond.** A Lennard-Jones or Buckingham term between
  two bonded atoms describes a bond that is weak or partly ionic, as between a
  metal and a ligand.

## Fixed bonds

Because bond vibrations are fast and of small amplitude, they are often
removed: the bond is given the type `"FIXED"` and its length is a constraint.
This is the choice of the united-atom alkane models (TraPPE), and it has two
advantages. In molecular dynamics the time step can be several times larger,
since it no longer has to resolve the bond vibration; the constraint is
imposed with SHAKE and RATTLE. In Monte Carlo, a molecule is grown and moved
with the bond length fixed, and no moves are spent on sampling a coordinate
that hardly varies.

A stiff harmonic bond and a fixed bond are not exactly equivalent: a
constraint removes a degree of freedom from the kinetic energy and changes
the distribution of the other coordinates slightly. The difference is small
for most properties.

## Sampling a bond length

In configurational-bias growth the length of a new bond is drawn from its own
Boltzmann distribution,

\f[
P(r) \propto r^2\,e^{-\beta u(r)} ,
\f]

where \f$r^2\f$ is the volume element of the sphere on which the new atom is
placed. Because the bond is stiff and independent of the other terms to a
good approximation, this is efficient, and the non-bonded energy then only
has to select among orientations.

## References

- P.M. Morse, *Diatomic Molecules According to the Wave Mechanics. II.
  Vibrational Levels*, Phys. Rev. **34**, 57 (1929).
- N.L. Allinger, Y.H. Yuh and J.-H. Lii, *Molecular mechanics. The MM3 force
  field for hydrocarbons. 1*, J. Am. Chem. Soc. **111**, 8551 (1989).
- J.R. Maple, M.-J. Hwang, T.P. Stockfisch, U. Dinur, M. Waldman, C.S. Ewig
  and A.T. Hagler, *Derivation of class II force fields. I. Methodology and
  quantum force field for the alkyl functional group and alkane molecules*,
  J. Comput. Chem. **15**, 162 (1994).
- J.-P. Ryckaert, G. Ciccotti and H.J.C. Berendsen, *Numerical integration of
  the Cartesian equations of motion of a system with constraints*,
  J. Comput. Phys. **23**, 327 (1977).
- B.G. Dick and A.W. Overhauser, *Theory of the Dielectric Constants of
  Alkali Halide Crystals*, Phys. Rev. **112**, 90 (1958).
