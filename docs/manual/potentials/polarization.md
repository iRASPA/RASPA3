# Polarization
\page potential_polarization Polarization

Fixed partial charges describe the *permanent* charge distribution of a
molecule. In an electric field the electron cloud of an atom is also
distorted, which gives it an *induced* dipole and lowers the energy. This is
polarization (induction). It matters where the fields are strong: near
cations and open metal sites of a framework, and for molecules with a large
polarizability and no permanent multipole, such as noble gases and alkanes,
whose interaction with a charged framework is otherwise missed entirely.

## The model

Each site \f$i\f$ has an isotropic polarizability \f$\alpha_i\f$. In the
electric field \f$\mathbf{E}_i\f$ at its position it acquires the dipole

\f[
\boldsymbol{\mu}_i = \alpha_i\,\mathbf{E}_i ,
\f]

and the energy of the induced dipole in the field, including the work needed
to create it, is

\f[
U_\text{pol} = -\frac{1}{2}\sum_i \alpha_i\,\left|\mathbf{E}_i\right|^2 .
\f]

The field is that of the *permanent charges* of the other molecules and of
the framework,

\f[
\mathbf{E}_i = \frac{1}{4\pi\epsilon_0}\sum_{j} q_j\,\frac{\mathbf{r}_{ij}}{r_{ij}^3} .
\f]

## What is and is not included

The induced dipoles respond to the field of the charges. The field of the
induced dipoles themselves on each other (*back-polarization*), which would
require solving for the dipoles self-consistently at every step, is not
included. This first-order treatment captures most of the induction energy in
adsorption, where the field of the framework dominates, at a small fraction
of the cost.

## Many-body character

The energy is quadratic in the field, and the field is a sum over charges, so
\f$|\mathbf{E}_i|^2\f$ contains products of contributions of *two* other
charges. Polarization is therefore not a sum of pair interactions: the energy
of a molecule depends on where all the others are, in a way that cannot be
tabulated per pair. Two practical consequences follow.

- The interaction with the framework cannot be taken from a pre-computed
  energy grid, because the electric field at each atom is needed. With
  polarization the framework interaction is evaluated pair by pair.
- A Monte Carlo move of one molecule changes the field at every other
  molecule. The moves that rely on the energy being a sum of independent
  contributions are restricted or not available with polarization; the
  [Monte Carlo moves](\ref mc_moves) page says so per move.

## Input

Polarization is switched on with `"ComputePolarization"` in the force field.
The polarizability of a site is a property of its pseudo-atom. The
polarization energy is reported separately in the energy averages and in the
energy histogram.

## Limits

- The gradient of the polarization energy is not evaluated by all force
  routines; properties and methods that need those forces (the force-based
  radial distribution function, the stress-fluctuation elastic constants, the
  spatial-decomposition engine) refuse to run with polarization.
- Force fields are fitted with or without explicit polarization. Adding it to
  a force field whose dispersion parameters were fitted without it counts the
  induction twice.

## References

- J. Applequist, J.R. Carl and K.-K. Fung, *Atom dipole interaction model for
  molecular polarizability*, J. Am. Chem. Soc. **94**, 2952 (1972).
- B.T. Thole, *Molecular polarizabilities calculated with a modified dipole
  interaction*, Chem. Phys. **59**, 341 (1981).
- T.M. Becker, D. Dubbeldam, L.-C. Lin and T.J.H. Vlugt, *Investigating
  polarization effects of CO2 adsorption in MgMOF-74*,
  J. Comput. Sci. **15**, 86 (2016).
- T.M. Becker, J. Heinen, D. Dubbeldam, L.-C. Lin and T.J.H. Vlugt,
  *Polarizable Force Fields for CO2 and CH4 Adsorption in M-MOF-74*,
  J. Phys. Chem. C **121**, 4659 (2017).
