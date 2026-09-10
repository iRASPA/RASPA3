module;

module energy_shared_well_surface;

import std;

import uint3;
import double3;
import skspacegroupdatabase;
import crystal;
import pair_interactions;
import units;
import blocking_spheres;
import energy_shared_linear_probe;
import energy_shared_electrostatic_potential_grid;
import energy_shared_energy_backend;
import energy_shared_probe_energy_grid;
import energy_shared_blocking_mask;

// Everything that has to know about Kelvin, kept apart from the arithmetic that does not.
//
// The measurement next door works in the units the field is held in and so needs nothing from the unit system;
// a temperature in Kelvin is a thing a *report* and a *driver* deal in. Splitting them is not tidiness: the
// conversion factors are mutable statics belonging to the engine, and a translation unit that touches them
// drags that dependency along with it. Keeping them out of the measurement is what lets the measurement be
// linked, and tested, against the structural library alone.

void writeWellSurface(std::ostream &stream, const WellSurface &surface)
{
  const double toKelvin = Units::EnergyToKelvin;

  std::print(stream, "#\n");
  std::print(stream, "# THE WELL SURFACE: WHERE THE MOLECULE ACTUALLY SITS\n");
  std::print(stream, "#\n");
  std::print(stream, "# The iso-surface of the energy is the inner turning point. Zero energy is where the\n");
  std::print(stream, "# repulsion balances the attraction on the way in, so it is the closest a molecule with\n");
  std::print(stream, "# nothing to spare gets before it is thrown back; it is not where it sits. It sits at the\n");
  std::print(stream, "# bottom of the well, a little further out, and that is the surface measured here.\n");
  std::print(stream, "#\n");
  std::print(stream, "# Construction. Topology comes from a distance field, geometry from the energy. d(x) is the\n");
  std::print(stream, "# additively weighted (Apollonius) distance to the framework, zero on the probe-contact\n");
  std::print(stream, "# offset surface. Marching cubes is run on max(-d, s (U - iso)), whose zero set is the\n");
  std::print(stream, "# boundary of the region the probe can occupy with a well at least iso deep. Being the\n");
  std::print(stream, "# boundary of a region, the mesh is watertight and single-sheeted: no interior membranes,\n");
  std::print(stream, "# no domes across intersections, no flaps at sheet junctions. Each vertex is then slid\n");
  std::print(stream, "# along the ray into the wall onto the 1D minimum of the analytic energy, which is the\n");
  std::print(stream, "# true multi-atom well floor.\n");
  std::print(stream, "#\n");
  std::print(stream, "# This replaces walking out from the zero-energy surface along the wall normal. That walk\n");
  std::print(stream, "# is the right physics on a smooth isolated wall and the wrong mesh on a real framework:\n");
  std::print(stream, "# the crease set of the energy contains one-sided sheets, and no local quantity separates\n");
  std::print(stream, "# them from the wall sheet. The contact surface cannot produce those, because d has exactly\n");
  std::print(stream, "# one zero crossing along any ray into a wall.\n");
  std::print(stream, "#\n");
  std::print(stream, "# NARROW CHANNELS. Where the pore is narrower than the probe's contact diameter there is no\n");
  std::print(stream, "# sheet to measure: the transverse minima have merged onto the channel axis and the well is\n");
  std::print(stream, "# a line, not a surface. That region is reported separately as the filament, a thin tube\n");
  std::print(stream, "# along the medial set. The ridge of that tube is the Gelb-Gubbins inscribed radius of\n");
  std::print(stream, "# the energy iso (distance to U = iso), the same field the energy PSD starts from, not\n");
  std::print(stream, "# the distance transform of the reliability overlay. A 1-D file packs along the graph\n");
  std::print(stream, "# length of that ridge (over the liquid diameter); a 2-D slit packs on the planar ridge\n");
  std::print(stream, "# (area / 16.2 Å²). A² / (4π V) of the whole tube is a cylinder check, not the packing\n");
  std::print(stream, "# of a mixed 1-D/2-D well. A tight-fitting probe therefore reports less sheet area than\n");
  std::print(stream, "# rolling it over the wall would, which is the honest answer for a surface that is not\n");
  std::print(stream, "# there.\n");
  std::print(stream, "#\n");
  std::print(stream, "# Blocking pockets close the surface off around themselves, which takes an inaccessible\n");
  std::print(stream, "# cage's own internal sheet out of the total.\n");
  std::print(stream, "#\n");
  std::print(stream, "# Trim isovalue: {:.4f} [K] (internal {:.6g}), energy scale {:.6g} [internal / Å]\n",
             surface.isoValue * toKelvin, surface.isoValue, surface.energyScale);
  std::print(stream, "# Triangles: {} on the sheet ({} discarded as implausibly large), {} vertices left on the\n",
             surface.numberOfTriangles, surface.numberOfRejectedTriangles, surface.numberOfTrimmedVertices);
  std::print(stream, "# trim cap; {} triangles on the filament after specks below {:.1f} [Å²] were removed\n",
             surface.numberOfFilamentTriangles, wellFilamentMinimumArea);
  std::print(stream, "# Timing: {} [s] for the sheet, the refinement and the filament\n", surface.seconds);
  std::print(stream, "\n");

  std::print(stream, "{:34} {:14.4f} [Å²] {:12.4f} [m²/g] {:12.4f} [m²/cm³]\n",
             "Well surface:", surface.area, surface.gravimetricArea, surface.volumetricArea);
  std::print(stream, "{:34} {:14.4f} [Å²]\n", "Filament area:", surface.filamentArea);
  std::print(stream, "{:34} {:14.4f} [Å³]\n", "Filament volume:", surface.filamentVolume);
  std::print(stream, "{:34} {:14.4f} [Å]\n", "Cylinder A²/(4πV):", surface.filamentLength);
  std::print(stream, "{:34} {:14.4f} [Å]\n", "1-D file (ridge):", surface.filamentRidgeLength);
  std::print(stream, "{:34} {:14.4f} [Å²]\n", "2-D midplane (ridge):", surface.filamentMedialArea);
  std::print(stream, "{:34} {:14.4f} [K]\n", "Mean well depth:", surface.meanDepth * toKelvin);
  std::print(stream, "{:34} {:14.4f} [K]\n", "Deepest well:", surface.deepestWell * toKelvin);
  std::print(stream, "\n");

  std::print(stream, "# The Boltzmann-weighted area is the integral of exp(-U_min/kT) over the well surface at\n");
  std::print(stream, "# {:.2f} [K]. On a level set there is no such quantity to be had --- the energy is the same\n",
             surface.thermalEnergy * toKelvin);
  std::print(stream, "# on all of it by construction, so the weight is one constant that divides straight back\n");
  std::print(stream, "# out. Here the depth varies over the surface, and the variation is physical rather than\n");
  std::print(stream, "# discretization: a shallow well on a convex cap facing a wide pore counts for little, a\n");
  std::print(stream, "# deep one in a corner touched on several sides at once counts for a great deal.\n");
  std::print(stream, "\n");
  std::print(stream, "{:34} {:14.4f} [Å²] {:12.4f} [m²/g]\n",
             "Boltzmann-weighted area:", surface.weightedArea, surface.gravimetricWeightedArea);
  std::print(stream, "{:34} {:14.5f} [-]\n", "Mean weight:", surface.enhancement());
}


MolecularWellSurface::MolecularWellSurface() {}

MolecularWellSurface::~MolecularWellSurface() {}

void MolecularWellSurface::run(const EnergyBackend &backend, const PairInteractions &interactions,
                               const Crystal &framework, const LinearProbe &probe, double level, uint3 gridSize,
                               std::size_t numberOfOrientations, double temperature,
                               std::span<const BlockingSphere> blockingSpheres, double energyScale,
                               double blockedEnergyPerAngstrom, double ceiling, bool useElectrostatics,
                               double relativePrecision)
{
  // 0.001 Å per kelvin of the field, expressed per internal energy unit.
  if (energyScale == 0.0) energyScale = wellEnergyScalePerKelvin * Units::EnergyToKelvin;
  if (blockedEnergyPerAngstrom == 0.0)
  {
    blockedEnergyPerAngstrom = blockedEnergyPerAngstromInKelvin * Units::KelvinToEnergy;
  }
  if (ceiling == 0.0) ceiling = probeEnergyCeilingInKelvin * Units::KelvinToEnergy;

  // A charged probe acts on the framework's potential: the smooth far half of the Ewald sum comes from the
  // backend as a grid, and the near half is summed inside the field alongside the dispersion, with the same
  // conversion factor the grid was built with.
  ElectrostaticPotentialGrid potential;
  bool wantsElectrostatics = useElectrostatics && probe.isCharged();
  if (wantsElectrostatics)
  {
    potential = backend.electrostaticPotentialGrid(interactions, framework, gridSize, relativePrecision);
  }
  const ElectrostaticPotentialGrid *potentialOrNothing = wantsElectrostatics ? &potential : nullptr;

  // The one temperature serves twice: as the kT of the orientational free energy the field is made of, and
  // as the Boltzmann weight over the finished sheet.
  this->isoValue = level;
  const double thermalEnergy = temperature * Units::KelvinToEnergy;
  this->field = backend.wellField
                    ? backend.wellField(interactions, framework, probe, gridSize, numberOfOrientations, thermalEnergy,
                                        blockingSpheres, blockedEnergyPerAngstrom, ceiling, potentialOrNothing,
                                        Units::CoulombicConversionFactor)
                    : computeWellField(interactions, framework, probe, gridSize, numberOfOrientations, thermalEnergy,
                                       blockingSpheres, blockedEnergyPerAngstrom, ceiling, potentialOrNothing,
                                       Units::CoulombicConversionFactor);

  if (this->field.numberOfVoxels() == 0) return;

  this->surface = wellSurfaceOfField(framework, interactions, this->field, level, energyScale,
                                     temperature * Units::KelvinToEnergy, blockingSpheres, &backend,
                                     potentialOrNothing, Units::CoulombicConversionFactor);

  double3 spacing = double3(this->field.unitCell.cell[0].length() / static_cast<double>(this->field.gridSize.x),
                            this->field.unitCell.cell[1].length() / static_cast<double>(this->field.gridSize.y),
                            this->field.unitCell.cell[2].length() / static_cast<double>(this->field.gridSize.z));
  double density = 1e-3 * framework.mass / (framework.unitCell.volume * Units::Angstrom * Units::Angstrom *
                                            Units::Angstrom * Units::AvogadroConstant);

  std::ofstream myfile;
  myfile.open(framework.name + "." + probe.name + ".energy.wells.txt");
  std::print(myfile, "# The adsorption surface: the well-floor contact sheet of the probe\n");
  std::print(myfile, "# Crystal: {}\n", framework.name);
  std::print(myfile, "# Space-group HM-symbol: {}\n",
             SKSpaceGroupDataBase::spaceGroupData[framework.spaceGroupHallNumber].HMString());
  std::print(myfile, "# Number of framework atoms: {}\n", framework.atoms.size());
  std::print(myfile, "# Crystal volume: {} [Å³]\n", framework.unitCell.volume);
  std::print(myfile, "# Crystal mass: {} [g/mol]\n", framework.mass);
  std::print(myfile, "# Crystal density: {} [kg/m³]\n", density);
  if (probe.sites.size() == 1)
  {
    std::print(myfile, "# Probe: {} (single site)\n", probe.name);
  }
  else
  {
    std::print(myfile, "# Probe: {} ({} sites over {:.3f} Å, {} orientations). The energy is the orientational\n",
               probe.name, probe.sites.size(), probe.length(), this->field.numberOfOrientations);
    std::print(myfile, "# free energy -kT ln <exp(-U/kT)>, so the turning a narrow pore denies the molecule is\n");
    std::print(myfile, "# paid for, and the contact distance is that of the best-fitting way round\n");
  }
  if (this->field.chargesIncluded)
  {
    std::print(myfile, "# The probe's partial charges act on the framework's potential: Ewald alpha {:.4f} [1/Å],\n",
               this->field.ewaldAlpha);
    std::print(myfile, "# {} wave vectors, the near half summed pair by pair alongside the dispersion; the charges\n",
               this->field.numberOfWaveVectors);
    std::print(myfile, "# deepen and turn the wells but take no part in the contact distance\n");
  }
  else if (this->field.chargesIgnored)
  {
    std::print(myfile, "# The probe's partial charges take no part here: the well field is dispersion only\n");
  }
  std::print(myfile, "# Grid: {} x {} x {} points, spacing {:.5f} x {:.5f} x {:.5f} [Å], far face left out\n",
             this->field.gridSize.x, this->field.gridSize.y, this->field.gridSize.z, spacing.x, spacing.y, spacing.z);
  std::print(myfile, "# Temperature: {} [K] (the orientational average and the Boltzmann weight)\n", temperature);
  std::print(myfile, "# Cutoff: {} [Å]\n", this->field.cutOff);
  std::print(myfile, "# Blocking spheres: {}\n", blockingSpheres.size());
  std::print(myfile, "# Timing: {} [s] for the well field\n", this->field.seconds);

  writeWellSurface(myfile, this->surface);
  myfile.close();
}
