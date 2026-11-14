# Examples External Fields
\page examples_external_fields Examples External Fields

These examples confine water to a channel by an external potential rather than
by a framework. The channel geometry is set in `force_field.json` with
`"ExternalFieldPotentialEnergySurface"`, and the system enables the field with
`"ExternalField" : true`.

## Table of Contents
1. [Monte Carlo: water in a cylindrical channel](#Example_external_1)
2. [Monte Carlo: water in a square channel](#Example_external_2)
3. [Monte Carlo: water in a cylindrical channel with grid interpolation](#Example_external_3)

----------------------------------------------------------------------------------

#### Monte Carlo: water in a cylindrical channel <a name="Example_external_1"></a>

100 water molecules in a \f$24.83 \times 24.83 \times 24.83\f$ &Aring; box at
298 K, confined to a cylinder of radius 5 &Aring; along z. The external field is
evaluated analytically (`"UseExternalFieldGrid" : false`).

Run from `examples/external_fields/1_mc_water_in_cylinder_channel`.

`simulation.json`:

```json
{
  "SimulationType" : "MonteCarlo",
  "NumberOfProductionCycles" : 1000,
  "NumberOfInitializationCycles" : 1000,
  "PrintEvery" : 100,

  "Systems" :
  [
    {
      "Type" : "Box",
      "BoxLengths" : [24.83, 24.83, 24.83],
      "ExternalTemperature" : 298.0,
      "ExternalField" : true,
      "ChargeMethod" : "Ewald",
      "OutputPDBMovie" : false,
      "SampleMovieEvery" : 10
    }
  ],

  "Components" :
  [
    {
      "Name" : "water",
      "TranslationProbability" : 0.5,
      "RotationProbability" : 0.5,
      "ReinsertionProbability" : 1.0,
      "CreateNumberOfMolecules" : 100
    }
  ]
}
```

`force_field.json` (channel settings):

```json
{
  "UseExternalFieldGrid" : false,
  "ExternalFieldPotentialEnergySurface" : "CylinderZ",
  "ExternalFieldCylinderRadius" : 5.0
}
```

#### Monte Carlo: water in a square channel <a name="Example_external_2"></a>

The same water Monte Carlo setup, confined to a rectangular channel of width
and height 5 &Aring; along z (`"RectangleZ"`). The width and the height are
half-widths measured from the channel axis through the centre of the cell: for
`"RectangleX"` the width is the extent along z and the height along y, for
`"RectangleY"` along x and z, and for `"RectangleZ"` along y and x. Inside the
channel (and inside the cylinder of `"CylinderX/Y/Z"`) the field is zero;
outside it is an overlap: trial positions outside are rejected, and a molecule
that starts outside (e.g. from a restart file) carries the overlap energy until
a move brings it back in.

Run from `examples/external_fields/2_mc_water_in_square_channel`.

```json
{
  "UseExternalFieldGrid" : false,
  "ExternalFieldPotentialEnergySurface" : "RectangleZ",
  "ExternalFieldRectangularChannelWidth" : 5.0,
  "ExternalFieldRectangularChannelHeight" : 5.0
}
```

#### Monte Carlo: water in a cylindrical channel with grid interpolation <a name="Example_external_3"></a>

The cylindrical channel of example 1, with the external field tabulated on a
\f$128^3\f$ grid and interpolated during the run (`"UseExternalFieldGrid" : true`).
A tabulated field couples linearly to the CFCMC fractional molecules
(\f$U = \lambda_\mathrm{VDW} U_\mathrm{ext}\f$, with the matching
\f$\partial U/\partial\lambda\f$ in the dU/dλ estimators). A field read from a
cube file (`"ExternalFieldPotentialEnergySurface" : "GridFile"`,
`"ExternalFieldGridFileName"`) exists only as this grid and therefore requires
`"UseExternalFieldGrid" : true`. A confining geometry (`CylinderX/Y/Z`, `RectangleX/Y/Z`) is a
hard wall, not an energy surface, and is never tabulated: with such a geometry
`"UseExternalFieldGrid"` has no effect and the wall is tested exactly. A periodic
cube file (N voxels a step L/N apart, as CP2K writes it) gets the plane at the
cell boundary as the copy of the first one; a file whose N points span the cell
end to end (as RASPA writes it) is read as it is.

Run from `examples/external_fields/3_mc_water_in_cylinder_channel_grids`.

```json
{
  "UseExternalFieldGrid" : true,
  "ExternalFieldPotentialEnergySurface" : "CylinderZ",
  "ExternalFieldCylinderRadius" : 5.0,
  "NumberOfExternalFieldGridPoints" : [128, 128, 128],
  "SpacingVDWGrid" : 0.15,
  "SpacingCoulombGrid" : 0.15,
  "NumberOfGridTestPoints" : 100000,
  "InterpolationScheme" : 3
}
```
