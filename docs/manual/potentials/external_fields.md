# External fields
\page potential_external_fields External fields

An external field is a potential that acts on each atom separately and
depends only on where the atom is,

\f[
U_\text{ext} = \sum_i u_\text{ext}(\mathbf{r}_i) .
\f]

It represents an environment that is not modelled atom by atom: a pore of
ideal shape, a wall, or a field computed by another method.

## The physics

A framework of atoms is itself an external field for the adsorbate when it is
rigid: the dispersion and electrostatic interaction of a molecule with all
framework atoms is a function of its position only. An explicit external
field replaces that sum by a given function. There are two reasons to do so.

- **Idealized confinement.** Much of the theory of fluids in pores is
  formulated for simple geometries: slits, cylinders, channels. Simulating a
  fluid in such a geometry separates the effect of confinement itself
  (layering, shifted phase transitions, single-file diffusion) from the
  chemical detail of a particular material, and allows direct comparison
  with theory.
- **A field from elsewhere.** A potential-energy surface from a
  quantum-chemical calculation, or from any other source, can be used as it
  is, tabulated on a grid.

## Confining geometries

A confining geometry is a region of the box inside which the field is zero
and outside which a molecule cannot be: a *hard wall*.

\f[
u_\text{ext}(\mathbf{r}) =
\begin{cases} 0 & \text{inside} \\ \infty & \text{outside} \end{cases}
\f]

| `"ExternalFieldPotentialEnergySurface"` | Region | Parameters |
|---|---|---|
| `"CylinderX"`, `"CylinderY"`, `"CylinderZ"` | a cylinder along the given axis, through the centre of the cell | `"ExternalFieldCylinderRadius"` [Å] |
| `"RectangleX"`, `"RectangleY"`, `"RectangleZ"` | a rectangular channel along the given axis | `"ExternalFieldRectangularChannelWidth"`, `"ExternalFieldRectangularChannelHeight"` [Å], half-widths from the axis |

A trial position outside the region is an overlap and is rejected. The wall
is tested exactly, not from a table. Because the box is periodic along the
axis of the channel, the channel is infinitely long.

A hard wall acts on the positions of the *atoms*. A molecule can therefore
approach the wall until one of its atoms touches it, and the accessible
volume for the centre of a molecule is smaller than the channel by the size
of the molecule.

## Tabulated fields

A field given on a grid is read from a cube file
(`"ExternalFieldPotentialEnergySurface" : "GridFile"` with
`"ExternalFieldGridFileName"`) and interpolated during the simulation. A
tabulated field is an energy surface, so it has finite values and gradients
everywhere, and can be used in molecular dynamics.

A tabulated field couples linearly to a fractional molecule,

\f[
U = \lambda\,U_\text{ext} ,
\f]

with the corresponding term in \f$\partial U/\partial\lambda\f$.

## In the input

The external field is switched on per system with `"ExternalField" : true`
in `simulation.json`. Its definition is in `force_field.json`:

```json
{
  "UseExternalFieldGrid" : false,
  "ExternalFieldPotentialEnergySurface" : "CylinderZ",
  "ExternalFieldCylinderRadius" : 5.0
}
```

The energy of the field is reported as a separate contribution in the energy
averages. Worked examples (water in a cylindrical and in a rectangular
channel, and a field on a grid) are on the
[external-field examples](\ref examples_external_fields) page.

## Limits

- A hard wall has no gradient, so the confining geometries are meant for
  Monte Carlo. Several properties and drivers that need forces refuse to run with
  an external field: the force-based radial distribution function, the
  elastic constants from stress fluctuations, the spatial-decomposition
  engine, and the exchange of configurations in parallel tempering.
- The pressure of a fluid in an external field is the pressure of a
  confined, inhomogeneous fluid. It is a tensor whose components along and
  across the channel differ, and it is not the pressure of the bulk fluid in
  equilibrium with it.

## References

- J.S. Rowlinson and B. Widom, *Molecular Theory of Capillarity*, Oxford
  University Press (1982).
- L.D. Gelb, K.E. Gubbins, R. Radhakrishnan and M. Sliwinska-Bartkowiak,
  *Phase separation in confined systems*, Rep. Prog. Phys. **62**, 1573
  (1999).
- R. Evans, *Fluids adsorbed in narrow pores: phase equilibria and
  structure*, J. Phys.: Condens. Matter **2**, 8989 (1990).
