module;

module energy_shared_bet_surface_area;

import std;

import uint3;
import double3;
import double3x3;
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
import energy_shared_well_surface;

void writeIsothermDerivedProperties(std::ostream &stream, const BETSurfaceArea &bet)
{
  if (bet.plateauReading)
  {
    std::print(stream, "# No window in this isotherm is a BET line: the pore fills within about a decade of x,\n");
    std::print(stream, "# so Rouquerol's ceiling sits at the top of the filling step and x/(n(1-x)) falls rather\n");
    std::print(stream, "# than rises everywhere below it. Read as the Type I isotherm it is instead: the monolayer\n");
    std::print(stream, "# is the plateau, and C is the Henry slope n/x divided by it. Neither depends on a window.\n");
  }
  if (!(bet.microporeVolume > 0.0)) return;

  std::print(stream, "{:34} {:14.4f} [mL / g] (Gurvich, {:.4f} molecules / cell as liquid)\n",
             "Micropore volume:", bet.microporeVolume, bet.saturationLoading);
  if (bet.tPlotNumberOfPoints >= 3)
  {
    std::print(stream, "{:34} {:14.4f} [mL / g]\n", "  t-plot micropore volume:", bet.tPlotMicroporeVolume);
    std::print(stream, "{:34} {:14.4f} [m²/g] ({} points, r² = {:.5f})\n",
               "  t-plot mesopore + external:", bet.tPlotExternalArea, bet.tPlotNumberOfPoints, bet.tPlotRSquared);
    std::print(stream, "# The t-plot compares the isotherm against the film a non-porous silica would have grown at\n");
    std::print(stream, "# the same relative pressure, over 3.2 to 4.2 Å of thickness: the intercept is read as the\n");
    std::print(stream, "# micropore volume and the slope as the surface a film is still free to grow on. A periodic\n");
    std::print(stream, "# crystal has no external surface, so that slope is not a surface but the bias of the\n");
    std::print(stream, "# reading itself, which experimental work has to calibrate away against reference solids.\n");
  }
}

void writeBETSurfaceArea(std::ostream &stream, const BETSurfaceArea &bet)
{
  std::print(stream, "#\n");
  std::print(stream, "# PREDICTED BET SURFACE AREA\n");
  std::print(stream, "#\n");
  std::print(stream, "# Not a geometric area. BET is a fit to an isotherm; the isotherm here is a Fowler–\n");
  std::print(stream, "# Guggenheim / Bragg–Williams mean-field lattice gas over the molecules the field will hold.\n");
  std::print(stream, "# Those are counted rather than measured: the attractive region is filled deepest-first,\n");
  std::print(stream, "# refusing to place a molecule within v_L^{{1/3}} of one already there, and each molecule so\n");
  std::print(stream, "# placed is one site carrying the energy of the point it took. A wide pore fills at liquid\n");
  std::print(stream, "# density because that spacing is the liquid's, and a channel too narrow for two abreast\n");
  std::print(stream, "# takes a single file, neither of which has to be decided in advance. Nitrogen attracts\n");
  std::print(stream, "# itself with W = -2 q_L, which makes bulk liquid condense at P0 by construction; the site\n");
  std::print(stream, "# occupancy is the Maxwell-construction minimiser of its grand potential, so deep wells fill\n");
  std::print(stream, "# first and shallow corners last. Packing dV / v_L over the attractive voxels is not used: a\n");
  std::print(stream, "# molecule's centre in an ultramicroporous channel is confined to a filament of a few Å³ and\n");
  std::print(stream, "# that reading undercounts it by an order of magnitude. The x -> 0 slope is anchored to the\n");
  std::print(stream, "# exact grid Henry coefficient (P0 V / kT) <exp(-U/kT)>.\n");
  std::print(stream, "#\n");
  std::print(stream, "# The Rouquerol consistency criteria pick the linear window, as they would on a measured or\n");
  std::print(stream, "# GCMC isotherm: the ceiling is the largest x to which n(1-x) still rises (at most the\n");
  std::print(stream, "# conventional 0.30) and the intercept must be positive. In a micropore the rise stops where\n");
  std::print(stream, "# the pores fill and the window follows it below the free-surface convention (Walton and\n");
  std::print(stream, "# Snurr, JACS 129, 8552 (2007)); the fitted n_m then can read the filling capacity rather than\n");
  std::print(stream, "# a geometric monolayer, which is the experimental artifact of BET on ultramicroporous\n");
  std::print(stream, "# materials and is left in, not corrected away.\n");
  std::print(stream, "#\n");
  std::print(stream, "# Nitrogen cross-section: {:.1f} [Å²]; liquid volume: {:.1f} [Å³]; linear spacing: {:.3f} [Å]; q_L: {:.2f} [kJ/mol]\n",
             nitrogenCrossSection, nitrogenLiquidVolume, nitrogenLinearSpacing, nitrogenHeatOfLiquefactionKJMol);
  std::print(stream, "# Fit window: {:.4g} -- {:.4g} in P/P0, r² = {:.5f}\n", bet.windowLow, bet.windowHigh,
             bet.rSquared);
  std::print(stream, "# Henry match prefactor f: {:.5g}\n", bet.henryAnchor);
  std::print(stream, "# Timing: {} [s] for the well field, the sheet and the fit\n", bet.seconds);
  std::print(stream, "\n");

  std::print(stream, "{:34} {:14.4f} [m²/g] {:12.4f} [m²/cm³]\n",
             "BET area:", bet.gravimetricArea, bet.volumetricArea);
  std::print(stream, "{:34} {:14.4f} [molecules / cell]\n", "Monolayer capacity n_m:", bet.monolayerCapacity);
  std::print(stream, "{:34} {:14.4f} [-]\n", "BET C constant:", bet.cConstant);
  std::print(stream, "{:34} {:14.4f} [molecules / cell]\n", "Saturation capacity (x -> 1):", bet.saturationCapacity);
  writeIsothermDerivedProperties(stream, bet);
  std::print(stream, "\n");
  std::print(stream, "{:34} {:14.4f} [molecules / cell]\n", "Packed into the field:", bet.fieldCapacity);
  std::print(stream, "{:34} {:14.4f} [neighbours] against {:.2f} in the bulk liquid\n",
             "Mean coordination:", bet.fieldCoordination, nitrogenLiquidCoordination);
  std::print(stream, "# That is the capacity the isotherm was built on. The rows below are the same pore read off\n");
  std::print(stream, "# the meshes drawn through the field instead --- a wall to spread a monolayer over, a file\n");
  std::print(stream, "# to thread, a midplane to tile --- which is a reading of the shape of the pore and not of\n");
  std::print(stream, "# what fits in it. They are the surface the well field describes; they no longer set the\n");
  std::print(stream, "# capacity. The fitted n_m need not equal either: heterogeneous condensation and the\n");
  std::print(stream, "# Rouquerol window can read more or less than the geometric monolayer.\n");
  std::print(stream, "{:34} {:14.4f} [molecules / cell] {:12.4f} [m²/g]\n",
             "Contact sheet:", bet.sheetCapacity, bet.sheetGravimetricArea);
  std::print(stream, "{:34} {:14.4f} [molecules / cell] {:12.4f} [m²/g]\n",
             "1-D file:", bet.filamentFileCapacity, bet.filamentFileGravimetricArea);
  std::print(stream, "{:34} {:14.4f} [molecules / cell] {:12.4f} [m²/g]\n",
             "2-D midplane:", bet.filamentMidplaneCapacity, bet.filamentMidplaneGravimetricArea);
  if (bet.filamentBlobCapacity > 0.0)
  {
    std::print(stream, "{:34} {:14.4f} [molecules / cell]\n", "Filament blob (V/v_L):", bet.filamentBlobCapacity);
  }
  std::print(stream, "{:34} {:14.4f} [molecules / cell]\n", "Filament capacity:", bet.filamentCapacity);
  if (bet.unsheetedFilamentShare > 0.0)
  {
    std::print(stream, "{:34} {:14.4f} [-]\n", "Filament not drawn as sheet:", bet.unsheetedFilamentShare);
    std::print(stream, "# The contact sheet and the filament overlay bound the same merged well, so the share of\n");
    std::print(stream, "# that boundary the sheet has drawn is a share of the pore already counted and only the\n");
    std::print(stream, "# rest is the filament's. A channel wide enough for a sheet keeps almost none of its\n");
    std::print(stream, "# medial curve; one that has closed over into a tube keeps nearly all of it.\n");
  }

  // The model isotherm itself, so the window can be judged the way the papers ask it to be reported:
  // n(1-x) is Rouquerol's first-criterion ordinate, x/(n(1-x)) the BET-plot ordinate.
  std::print(stream, "\n");
  std::print(stream, "# {:>12} {:>16} {:>16} {:>16}\n", "x = P/P0", "n [molec/cell]", "n(1-x)", "x/(n(1-x))");
  for (const IsothermPoint &point : bet.isotherm)
  {
    double g = point.moleculesPerCell * (1.0 - point.relativePressure);
    double y = (g > 0.0) ? point.relativePressure / g : 0.0;
    std::print(stream, "  {:>12.6g} {:>16.8g} {:>16.8g} {:>16.8g}\n", point.relativePressure, point.moleculesPerCell,
               g, y);
  }
}


void BETSurfaceArea::run(const EnergyBackend &backend, const PairInteractions &interactions, const Crystal &framework,
                         const LinearProbe &probe, double level, uint3 gridSize, std::size_t numberOfOrientations,
                         double temperature, std::span<const BlockingSphere> blockingSpheres, double energyScale,
                         double blockedEnergyPerAngstrom, double ceiling, bool useElectrostatics,
                         double relativePrecision)
{
  std::chrono::steady_clock::time_point time_begin = std::chrono::steady_clock::now();

  // 0.001 Å per kelvin of the field, expressed per internal energy unit.
  if (energyScale == 0.0) energyScale = wellEnergyScalePerKelvin * Units::EnergyToKelvin;
  if (blockedEnergyPerAngstrom == 0.0)
  {
    blockedEnergyPerAngstrom = blockedEnergyPerAngstromInKelvin * Units::KelvinToEnergy;
  }
  if (ceiling == 0.0) ceiling = probeEnergyCeilingInKelvin * Units::KelvinToEnergy;
  if (temperature <= 0.0) temperature = nitrogenBETTemperature;

  // A charged probe acts on the framework's potential, exactly as in the well surface it stands on.
  ElectrostaticPotentialGrid potential;
  bool wantsElectrostatics = useElectrostatics && probe.isCharged();
  if (wantsElectrostatics)
  {
    potential = backend.electrostaticPotentialGrid(interactions, framework, gridSize, relativePrecision);
  }
  const ElectrostaticPotentialGrid *potentialOrNothing = wantsElectrostatics ? &potential : nullptr;

  // For a molecule the field is the orientational free energy, which is exactly what the Henry limit of a
  // rigid rotor wants: (V/kT) <exp(-U/kT)> over positions and orientations is <exp(-F/kT)> over positions.
  const double thermalEnergy = temperature * Units::KelvinToEnergy;
  this->field = backend.wellField
                    ? backend.wellField(interactions, framework, probe, gridSize, numberOfOrientations, thermalEnergy,
                                        blockingSpheres, blockedEnergyPerAngstrom, ceiling, potentialOrNothing,
                                        Units::CoulombicConversionFactor)
                    : computeWellField(interactions, framework, probe, gridSize, numberOfOrientations, thermalEnergy,
                                       blockingSpheres, blockedEnergyPerAngstrom, ceiling, potentialOrNothing,
                                       Units::CoulombicConversionFactor);
  this->surface = wellSurfaceOfField(framework, interactions, this->field, level, energyScale,
                                     temperature * Units::KelvinToEnergy, blockingSpheres, &backend,
                                     potentialOrNothing, Units::CoulombicConversionFactor);

  const double beta = 1.0 / (Units::KB * temperature);
  double boltzmann = 0.0;
  if (!this->field.energy.empty())
  {
    for (float energy : this->field.energy) boltzmann += std::exp(-beta * static_cast<double>(energy));
    boltzmann /= static_cast<double>(this->field.energy.size());
  }

  // n = (P0 V / kT) <exp(-U/kT)>, the Henry occupancy of the cell at saturation pressure, in molecules.
  const double volumeSI = framework.unitCell.volume * Units::Angstrom * Units::Angstrom * Units::Angstrom;
  const double kT = Units::BoltzmannConstant * temperature;
  const double henryMoleculesPerCell = (kT > 0.0) ? nitrogenSaturationPressure * volumeSI / kT * boltzmann : 0.0;

  const double heatOfLiquefaction = nitrogenHeatOfLiquefactionKJMol / Units::EnergyToKJPerMol;

  // How many monolayers the pore actually holds: diagnostic only (the lattice-gas isotherm is monolayer
  // with lateral interactions; micropore type-I shape comes from the filling step itself).
  this->multilayerRoom = 0.0;
  const double sheetSites = this->surface.area / nitrogenCrossSection;

  if (sheetSites > 0.0 && !this->field.energy.empty())
  {
    const float trim = static_cast<float>(this->surface.isoValue);
    std::size_t poreVoxels = 0;
    for (float energy : this->field.energy)
    {
      if (energy < trim) ++poreVoxels;
    }
    const double poreVolume =
        framework.unitCell.volume * static_cast<double>(poreVoxels) / static_cast<double>(this->field.energy.size());
    this->multilayerRoom =
        std::max(1.0, (poreVolume - this->surface.filamentVolume) / (sheetSites * nitrogenLiquidVolume));
  }

  // The contact sheet and the filament overlay bound the same merged well, so the sheet's share of that
  // boundary is a share of the pore it has already counted and the filament keeps only the rest.
  this->unsheetedFilamentShare =
      (this->surface.filamentArea > 0.0)
          ? std::max(0.0, 1.0 - std::min(1.0, this->surface.area / this->surface.filamentArea))
          : 0.0;
  GeometricSites mesh =
      geometricAdsorptionSites(this->surface.patches, this->surface.filamentVoxels, this->surface.filamentArea);
  this->sheetCapacity = mesh.sheetCapacity;
  this->filamentFileCapacity = mesh.filamentFileCapacity;
  this->filamentMidplaneCapacity = mesh.filamentMidplaneCapacity;
  this->filamentBlobCapacity = mesh.filamentBlobCapacity;
  this->filamentCapacity = mesh.filamentCapacity;

  // The isotherm is run on the field's own reading of the capacity. The mesh capacities above stay in the
  // report as the surface they describe, but they are a reading of the shape of the pore and not of what
  // fits in it, and that is where the two of them parted company with a simulated isotherm: mordenite 43%
  // under, beta 28% over.
  GeometricSites geometry = fieldAdsorptionSites(this->field);
  this->fieldCapacity = static_cast<double>(geometry.sites.size());
  std::size_t pairs = 0uz;
  for (const std::vector<std::uint32_t> &list : geometry.neighbours) pairs += list.size();
  this->fieldCoordination =
      geometry.neighbours.empty() ? 0.0 : static_cast<double>(pairs) / static_cast<double>(geometry.neighbours.size());

  constexpr double angstromSquaredToSquareMetrePerMol = 6.0221419947e3;
  auto gravimetricOf = [&](double molecules)
  {
    return (framework.mass > 0.0)
               ? molecules * nitrogenCrossSection * angstromSquaredToSquareMetrePerMol / framework.mass
               : 0.0;
  };
  this->sheetGravimetricArea = gravimetricOf(mesh.sheetCapacity);
  this->filamentFileGravimetricArea = gravimetricOf(mesh.filamentFileCapacity);
  this->filamentMidplaneGravimetricArea = gravimetricOf(mesh.filamentMidplaneCapacity);

  LatticeIsotherm lattice = latticeGasIsotherm(geometry.sites, henryMoleculesPerCell,
                                               temperature * Units::KelvinToEnergy, heatOfLiquefaction,
                                               geometry.neighbours);
  BETSurfaceArea latticeFit =
      BETSurfaceArea::fromIsotherm(std::move(lattice.isotherm), framework.mass, framework.unitCell.volume,
                                   nitrogenCrossSection, nitrogenLiquidVolume);

  this->gravimetricArea = latticeFit.gravimetricArea;
  this->volumetricArea = latticeFit.volumetricArea;
  this->monolayerCapacity = latticeFit.monolayerCapacity;
  this->cConstant = latticeFit.cConstant;
  this->windowLow = latticeFit.windowLow;
  this->windowHigh = latticeFit.windowHigh;
  this->rSquared = latticeFit.rSquared;
  this->henryAnchor = lattice.henryPrefactor;
  this->saturationCapacity = lattice.saturationCapacity;
  this->saturationLoading = latticeFit.saturationLoading;
  this->microporeVolume = latticeFit.microporeVolume;
  this->tPlotMicroporeVolume = latticeFit.tPlotMicroporeVolume;
  this->tPlotExternalArea = latticeFit.tPlotExternalArea;
  this->tPlotRSquared = latticeFit.tPlotRSquared;
  this->tPlotNumberOfPoints = latticeFit.tPlotNumberOfPoints;
  this->plateauReading = latticeFit.plateauReading;
  this->isotherm = std::move(latticeFit.isotherm);

  std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - time_begin;
  this->seconds = elapsed.count();

  double3 spacing = double3(this->field.unitCell.cell[0].length() / static_cast<double>(this->field.gridSize.x),
                            this->field.unitCell.cell[1].length() / static_cast<double>(this->field.gridSize.y),
                            this->field.unitCell.cell[2].length() / static_cast<double>(this->field.gridSize.z));
  double density = 1e-3 * framework.mass / (framework.unitCell.volume * Units::Angstrom * Units::Angstrom *
                                            Units::Angstrom * Units::AvogadroConstant);

  std::ofstream myfile;
  myfile.open(framework.name + "." + probe.name + ".energy.bet.txt");
  std::print(myfile, "# Predicted BET surface area from the well-surface contact sheet\n");
  std::print(myfile, "# Crystal: {}\n", framework.name);
  std::print(myfile, "# Space-group HM-symbol: {}\n",
             SKSpaceGroupDataBase::spaceGroupData[framework.spaceGroupHallNumber].HMString());
  std::print(myfile, "# Number of framework atoms: {}\n", framework.atoms.size());
  std::print(myfile, "# Crystal volume: {} [Å³]\n", framework.unitCell.volume);
  std::print(myfile, "# Crystal mass: {} [g/mol]\n", framework.mass);
  std::print(myfile, "# Crystal density: {} [kg/m³]\n", density);
  if (probe.sites.size() == 1)
  {
    std::print(myfile, "# Probe: {}\n", probe.name);
  }
  else
  {
    std::print(myfile, "# Probe: {} ({} sites over {:.3f} Å, {} orientations)\n", probe.name, probe.sites.size(),
               probe.length(), this->field.numberOfOrientations);
  }
  if (this->field.chargesIncluded)
  {
    std::print(myfile, "# The probe's partial charges act on the framework's potential (Ewald alpha {:.4f} [1/Å],\n",
               this->field.ewaldAlpha);
    std::print(myfile, "# {} wave vectors); they take no part in the contact distance\n",
               this->field.numberOfWaveVectors);
  }
  else if (this->field.chargesIgnored)
  {
    std::print(myfile, "# The probe's partial charges take no part here: the well field is dispersion only\n");
  }
  std::print(myfile, "# Grid: {} x {} x {} points, spacing {:.5f} x {:.5f} x {:.5f} [Å]\n", this->field.gridSize.x,
             this->field.gridSize.y, this->field.gridSize.z, spacing.x, spacing.y, spacing.z);
  std::print(myfile, "# Temperature: {} [K]\n", temperature);
  std::print(myfile, "# Blocking spheres: {}\n", blockingSpheres.size());
  std::print(myfile, "# Henry occupancy at P0: {:.6g} [molecules / cell]; <exp(-U/kT)> = {:.6g}\n",
             henryMoleculesPerCell, boltzmann);

  // Where the Henry average comes from. <exp(-U/kT)> is a mean over voxels, so a small number of very deep
  // ones can carry all of it while the contact sheet the BET line is fitted to knows nothing about them.
  // The table says how much of the average survives when wells past a given depth are discarded, and the
  // listing names the framework atom each of the deepest voxels sits on, which is what tells a genuine
  // adsorption site apart from a singularity in the potential.
  if (!this->field.energy.empty())
  {
    const std::size_t numberOfVoxels = this->field.energy.size();
    const double totalWeight = boltzmann * static_cast<double>(numberOfVoxels);

    std::print(myfile, "#\n");
    std::print(myfile, "# Depth profile of the Henry average:\n");
    std::print(myfile, "# {:>10}  {:>12}  {:>14}  {:>12}\n", "deeper than", "voxels", "fraction of", "<exp> without");
    std::print(myfile, "# {:>10}  {:>12}  {:>14}  {:>12}\n", "[kT]", "[-]", "<exp(-U/kT)>", "them");
    for (const double threshold : {10.0, 20.0, 30.0, 40.0, 60.0, 80.0, 100.0})
    {
      std::size_t count = 0uz;
      double weight = 0.0;
      for (float energy : this->field.energy)
      {
        const double betaU = beta * static_cast<double>(energy);
        if (betaU < -threshold)
        {
          ++count;
          weight += std::exp(-betaU);
        }
      }
      const double remaining = (totalWeight - weight) / static_cast<double>(numberOfVoxels);
      std::print(myfile, "# {:>10.0f}  {:>12d}  {:>14.6g}  {:>12.6g}\n", threshold, count,
                 (totalWeight > 0.0) ? weight / totalWeight : 0.0, remaining);
    }

    // The deepest handful, each placed against the framework atom nearest to it.
    constexpr std::size_t numberOfDeepestToList = 8uz;
    std::vector<std::size_t> order(numberOfVoxels);
    std::iota(order.begin(), order.end(), 0uz);
    const std::size_t listed = std::min(numberOfDeepestToList, numberOfVoxels);
    std::partial_sort(order.begin(), order.begin() + static_cast<std::ptrdiff_t>(listed), order.end(),
                      [&](std::size_t a, std::size_t b) { return this->field.energy[a] < this->field.energy[b]; });

    const uint3 size = this->field.gridSize;
    const double3x3 cell = framework.unitCell.cell;
    std::print(myfile, "#\n");
    std::print(myfile, "# The {} deepest voxels and the framework atom each is nearest to:\n", listed);
    std::print(myfile, "# {:>8} {:>8} {:>8}  {:>14}  {:>10}  {:>8}  {:>10}\n", "frac a", "frac b", "frac c", "U [K]",
               "U/kT", "atom", "r [Å]");
    for (std::size_t rank = 0; rank < listed; ++rank)
    {
      const std::size_t voxel = order[rank];
      const std::size_t ix = voxel % size.x;
      const std::size_t iy = (voxel / size.x) % size.y;
      const std::size_t iz = voxel / (static_cast<std::size_t>(size.x) * size.y);
      const double3 point(static_cast<double>(ix) / static_cast<double>(size.x),
                          static_cast<double>(iy) / static_cast<double>(size.y),
                          static_cast<double>(iz) / static_cast<double>(size.z));

      double nearestSquared = std::numeric_limits<double>::max();
      std::size_t nearestAtom = 0uz;
      for (std::size_t iatom = 0; iatom < framework.fractionalPositions.size(); ++iatom)
      {
        double3 ds = point - framework.fractionalPositions[iatom];
        ds.x -= std::rint(ds.x);
        ds.y -= std::rint(ds.y);
        ds.z -= std::rint(ds.z);
        const double3 dr = cell * ds;
        const double rr = double3::dot(dr, dr);
        if (rr < nearestSquared)
        {
          nearestSquared = rr;
          nearestAtom = iatom;
        }
      }

      const double energyKelvin = static_cast<double>(this->field.energy[voxel]) * Units::EnergyToKelvin;
      std::print(myfile, "# {:>8.5f} {:>8.5f} {:>8.5f}  {:>14.6g}  {:>10.4g}  {:>8d}  {:>10.4f}\n", point.x, point.y,
                 point.z, energyKelvin, beta * static_cast<double>(this->field.energy[voxel]),
                 framework.atoms[nearestAtom].type, std::sqrt(nearestSquared));
    }
  }

  // What the one constant of the packing is worth. The capacity is the count of molecules the field will
  // take at a given closest approach, and how fast it falls with that distance says which limit the pore
  // is in: a single file loses molecules as 1/d, a filled cage as 1/d³. The field is already built, so
  // the whole curve costs no more than the packing itself.
  {
    std::print(myfile, "#\n");
    std::print(myfile, "# Capacity against the closest approach allowed between two molecules:\n");
    std::print(myfile, "# {:>10}  {:>14}  {:>16}\n", "d [Å]", "molecules/cell", "mean neighbours");
    for (const double distance : {3.2, 3.4, 3.5, 3.6, nitrogenLinearSpacing, 4.1, 4.4})
    {
      const GeometricSites trial = fieldAdsorptionSites(this->field, distance);
      std::size_t trialPairs = 0uz;
      for (const std::vector<std::uint32_t> &list : trial.neighbours) trialPairs += list.size();
      const double coordination =
          trial.neighbours.empty() ? 0.0
                                   : static_cast<double>(trialPairs) / static_cast<double>(trial.neighbours.size());
      std::print(myfile, "# {:>10.3f}  {:>14d}  {:>16.3f}\n", distance, trial.sites.size(), coordination);
    }
  }

  writeBETSurfaceArea(myfile, *this);
  writeWellSurface(myfile, this->surface);
  myfile.close();
}
