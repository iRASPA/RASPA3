module;

module energy_shared_nldft;

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
import energy_shared_well_field;
import energy_shared_well_surface;
import energy_shared_bet_surface_area;

NLDFTIsotherm::NLDFTIsotherm() {}

NLDFTIsotherm::~NLDFTIsotherm() {}

void NLDFTIsotherm::run(const EnergyBackend &backend, const PairInteractions &interactions, const Crystal &framework,
                        const LinearProbe &probe, uint3 gridSize, std::size_t numberOfOrientations,
                        double temperature, std::span<const BlockingSphere> blockingSpheres, bool useElectrostatics,
                        double relativePrecision)
{
  std::chrono::steady_clock::time_point time_begin = std::chrono::steady_clock::now();

  if (temperature <= 0.0) temperature = nitrogenBETTemperature;
  const double blockedEnergyPerAngstrom = blockedEnergyPerAngstromInKelvin * Units::KelvinToEnergy;
  const double ceiling = probeEnergyCeilingInKelvin * Units::KelvinToEnergy;

  // A charged probe acts on the framework's potential, exactly as in the well field.
  ElectrostaticPotentialGrid potential;
  bool wantsElectrostatics = useElectrostatics && probe.isCharged();
  if (wantsElectrostatics)
  {
    potential = backend.electrostaticPotentialGrid(interactions, framework, gridSize, relativePrecision);
  }
  const ElectrostaticPotentialGrid *potentialOrNothing = wantsElectrostatics ? &potential : nullptr;

  const double thermalEnergy = temperature * Units::KelvinToEnergy;
  // A linear guest (`--nldft --molecule N2`) needs U(r, ω) on the CPU field: GPU well-field kernels
  // keep only the Helmholtz average. A spherical probe keeps the existing backend (GPU if present).
  const bool wantsMolecular = probe.sites.size() > 1 && numberOfOrientations > 1;
  if (wantsMolecular || !backend.wellField)
  {
    this->field = computeWellField(interactions, framework, probe, gridSize, numberOfOrientations, thermalEnergy,
                                   blockingSpheres, blockedEnergyPerAngstrom, ceiling, potentialOrNothing,
                                   Units::CoulombicConversionFactor);
  }
  else
  {
    this->field = backend.wellField(interactions, framework, probe, gridSize, numberOfOrientations, thermalEnergy,
                                    blockingSpheres, blockedEnergyPerAngstrom, ceiling, potentialOrNothing,
                                    Units::CoulombicConversionFactor);
  }

  std::vector<float> energyKelvin(this->field.energy.size());
  for (std::size_t i = 0; i < energyKelvin.size(); ++i)
  {
    energyKelvin[i] = static_cast<float>(static_cast<double>(this->field.energy[i]) * Units::EnergyToKelvin);
  }

  const bool molecular =
      wantsMolecular && this->field.orientationEnergy.size() == energyKelvin.size() * numberOfOrientations;
  this->molecular = molecular;
  this->numberOfOrientations = molecular ? numberOfOrientations : 1;

  NLDFTOptions options;
  options.temperature = temperature;
  options.headTailSymmetric = probe.headTailSymmetric;
  if (molecular)
  {
    // Guest–guest TraPPE N–N (or the force field's dispersing sites): full 12-6 from r = σ, not WCA.
    // Hard-sphere FMT uses the same dispersing sites as a fused dumbbell at TraPPE σ — matching MC
    // packing of the N atoms — not the older Ravikovitch COM sphere (d = 3.575 Å). N_com has no LJ
    // and is skipped for both FMT and the attractive mean field. Guest–guest Coulomb is left out of
    // F_att: a cutoff 1/r mean field of the partial charges is monopole-dominated and wrong for a
    // neutral quadrupole; orientation-averaged QQ (~1/r^5) is not implemented.
    for (const auto &site : probe.sites)
    {
      const PairParameters &self = interactions[site.type];
      if (!(self.sizeParameter > 0.0) || !(self.strengthParameter > 0.0)) continue;
      if (options.attractiveSiteOffsets.empty())
      {
        options.attractiveSigma = self.sizeParameter;
        options.attractiveEpsilon = self.strengthParameter * Units::EnergyToKelvin;
      }
      else if (std::abs(self.sizeParameter - options.attractiveSigma) > 1.0e-8 ||
               std::abs(self.strengthParameter * Units::EnergyToKelvin - options.attractiveEpsilon) > 1.0e-8)
      {
        continue;
      }
      options.attractiveSiteOffsets.push_back(site.offset);
    }
    if (!options.attractiveSiteOffsets.empty())
    {
      options.dumbbellSiteOffsets = options.attractiveSiteOffsets;
      options.dumbbellSiteDiameter = options.attractiveSigma;
      options.hardSphereDiameter = options.attractiveSigma;
      options.sigma = options.attractiveSigma;
      options.cutoff = nldftCutoffFactor * options.attractiveSigma;
    }
  }
  const float ceilingKelvin = static_cast<float>(probeEnergyCeilingInKelvin);
  NLDFTExternalFieldStats fieldStats{};
  if (!molecular)
  {
    fieldStats = nldftMaskAndClipExternalPotential(energyKelvin, temperature, ceilingKelvin, options.wellFloorInKT);
  }

  std::vector<float> orientationKelvin;
  if (molecular)
  {
    orientationKelvin.resize(this->field.orientationEnergy.size());
    for (std::size_t i = 0; i < orientationKelvin.size(); ++i)
    {
      orientationKelvin[i] =
          static_cast<float>(static_cast<double>(this->field.orientationEnergy[i]) * Units::EnergyToKelvin);
    }
  }

  // The exact grid Henry limit on the unclipped U(r, ω). This is the GCMC ⟨e^{-βU}⟩ on the same grid.
  double boltzmann = 0.0;
  if (molecular && !orientationKelvin.empty())
  {
    const std::size_t nVoxels = energyKelvin.size();
    const std::size_t M = this->numberOfOrientations;
    for (std::size_t i = 0; i < nVoxels; ++i)
    {
      double sum = 0.0;
      for (std::size_t o = 0; o < M; ++o)
        sum += std::exp(-static_cast<double>(orientationKelvin[i * M + o]) / temperature);
      boltzmann += sum / static_cast<double>(M);
    }
    boltzmann /= static_cast<double>(nVoxels);
  }
  else if (!energyKelvin.empty())
  {
    for (float energy : energyKelvin) boltzmann += std::exp(-static_cast<double>(energy) / temperature);
    boltzmann /= static_cast<double>(energyKelvin.size());
  }
  const double volumeSI = framework.unitCell.volume * Units::Angstrom * Units::Angstrom * Units::Angstrom;
  const double kT = Units::BoltzmannConstant * temperature;
  this->henryMoleculesPerCell = (kT > 0.0) ? nitrogenSaturationPressure * volumeSI / kT * boltzmann : 0.0;

  struct DeepPose
  {
    double kelvin{0.0};
    std::size_t voxel{0};
    std::size_t orientation{0};
  };
  DeepPose deepest[3]{};
  std::size_t nCeiling = 0, nBelow80 = 0, nBelow40 = 0, nBelow20 = 0, nBelow12 = 0;
  std::size_t voxelsAttractive = 0, voxelsOpen = 0;
  double boltzmannBelow40 = 0.0, boltzmannCapped12 = 0.0;
  if (molecular && !orientationKelvin.empty())
  {
    const std::size_t nVoxels = energyKelvin.size();
    const std::size_t M = this->numberOfOrientations;
    const double kTkelvin = temperature;
    const double cap12 = 12.0 * kTkelvin;
    const float halfCeiling = 0.5f * ceilingKelvin;
    for (std::size_t i = 0; i < nVoxels; ++i)
    {
      bool open = false;
      bool attractive = false;
      double voxelBelow40 = 0.0;
      double voxelCap12 = 0.0;
      for (std::size_t o = 0; o < M; ++o)
      {
        const double u = static_cast<double>(orientationKelvin[i * M + o]);
        const double w = std::exp(-u / kTkelvin);
        voxelCap12 += std::exp(-std::max(u, -cap12) / kTkelvin);
        if (u >= halfCeiling)
        {
          ++nCeiling;
        }
        else
        {
          open = true;
          if (u < 0.0) attractive = true;
          if (u < -80.0 * kTkelvin) ++nBelow80;
          if (u < -40.0 * kTkelvin)
          {
            ++nBelow40;
            voxelBelow40 += w;
          }
          if (u < -20.0 * kTkelvin) ++nBelow20;
          if (u < -12.0 * kTkelvin) ++nBelow12;
          if (u < deepest[2].kelvin)
          {
            DeepPose pose{u, i, o};
            if (pose.kelvin < deepest[0].kelvin)
            {
              deepest[2] = deepest[1];
              deepest[1] = deepest[0];
              deepest[0] = pose;
            }
            else if (pose.kelvin < deepest[1].kelvin)
            {
              deepest[2] = deepest[1];
              deepest[1] = pose;
            }
            else
            {
              deepest[2] = pose;
            }
          }
        }
      }
      if (open) ++voxelsOpen;
      if (attractive) ++voxelsAttractive;
      boltzmannBelow40 += voxelBelow40 / static_cast<double>(M);
      boltzmannCapped12 += voxelCap12 / static_cast<double>(M);
    }
    boltzmannBelow40 /= static_cast<double>(nVoxels);
    boltzmannCapped12 /= static_cast<double>(nVoxels);
  }

  // Pressure window from the grid Henry, same rule as WHAM/TMMC: below the fill step for strong
  // binders, up to the isotherm's own P0. Preview coexistence so x uses the same P0 the solve will.
  {
    const NLDFTBulk bulkPreview = nldftBulkCoexistence(options);
    nldftApplyHenryPressureWindow(options, this->henryMoleculesPerCell, framework.unitCell.volume,
                                  bulkPreview.saturationPressure);
  }

  NLDFTGridIsotherm grid =
      molecular ? nldftIsothermOnGrid(energyKelvin, this->field.gridSize, this->field.unitCell, options,
                                      orientationKelvin, this->numberOfOrientations)
                : nldftIsothermOnGrid(energyKelvin, this->field.gridSize, this->field.unitCell, options);
  this->bulk = grid.bulk;
  this->unconverged = grid.unconverged;
  this->isotherm = std::move(grid.isotherm);

  this->fit = BETSurfaceArea::fromIsotherm(this->isotherm, framework.mass, framework.unitCell.volume,
                                           nitrogenCrossSection, nitrogenLiquidVolume);

  std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - time_begin;
  this->seconds = elapsed.count();

  double3 spacing = double3(this->field.unitCell.cell[0].length() / static_cast<double>(this->field.gridSize.x),
                            this->field.unitCell.cell[1].length() / static_cast<double>(this->field.gridSize.y),
                            this->field.unitCell.cell[2].length() / static_cast<double>(this->field.gridSize.z));
  double density = 1e-3 * framework.mass / (framework.unitCell.volume * Units::Angstrom * Units::Angstrom *
                                            Units::Angstrom * Units::AvogadroConstant);

  std::ofstream myfile;
  myfile.open(framework.name + "." + probe.name + ".energy.nldft.txt");
  std::print(myfile, "# Predicted nitrogen isotherm by classical density functional theory on the energy grid\n");
  std::print(myfile, "# Crystal: {}\n", framework.name);
  std::print(myfile, "# Space-group HM-symbol: {}\n",
             SKSpaceGroupDataBase::spaceGroupData[framework.spaceGroupHallNumber].HMString());
  std::print(myfile, "# Number of framework atoms: {}\n", framework.atoms.size());
  std::print(myfile, "# Crystal volume: {} [Å³]\n", framework.unitCell.volume);
  std::print(myfile, "# Crystal mass: {} [g/mol]\n", framework.mass);
  std::print(myfile, "# Crystal density: {} [kg/m³]\n", density);
  if (probe.sites.size() == 1)
  {
    std::print(myfile, "# Probe (external field): {}\n", probe.name);
  }
  else
  {
    std::print(myfile, "# Probe (external field): {} ({} sites over {:.3f} Å, {} orientations)\n", probe.name,
               probe.sites.size(), probe.length(), this->field.numberOfOrientations);
  }
  if (this->field.chargesIncluded)
  {
    std::print(myfile, "# The probe's partial charges act on the framework's potential (Ewald alpha {:.4f} [1/Å],\n",
               this->field.ewaldAlpha);
    std::print(myfile, "# {} wave vectors)\n", this->field.numberOfWaveVectors);
  }
  else if (this->field.chargesIgnored)
  {
    std::print(myfile, "# The probe's partial charges take no part here: the external field is dispersion only\n");
  }
  std::print(myfile, "# Grid: {} x {} x {} points, spacing {:.5f} x {:.5f} x {:.5f} [Å]\n", this->field.gridSize.x,
             this->field.gridSize.y, this->field.gridSize.z, spacing.x, spacing.y, spacing.z);
  std::print(myfile, "# Temperature: {} [K]\n", temperature);
  std::print(myfile, "# Blocking spheres: {}\n", blockingSpheres.size());
  if (this->molecular)
  {
    if (!options.dumbbellSiteOffsets.empty())
    {
      std::print(myfile,
                 "# Density: ρ(r, ω) for the linear guest. White Bear FMT is a fused dumbbell at TraPPE σ "
                 "({:.3f} Å, {} sites).\n",
                 options.dumbbellSiteDiameter, options.dumbbellSiteOffsets.size());
    }
    else
    {
      std::print(myfile, "# Density: ρ(r, ω) for the linear guest. White Bear FMT stays on the COM (d = {:.3f} Å).\n",
                 nldftHardSphereDiameter);
    }
    if (!options.attractiveSiteOffsets.empty())
    {
      std::print(myfile,
                 "# Guest–guest attraction is TraPPE site-site 12-6 from r = σ ({:.3f} Å, ε/k {:.2f} K) at {} "
                 "dispersing sites, not the Ravikovitch WCA tail.\n",
                 options.attractiveSigma, options.attractiveEpsilon, options.attractiveSiteOffsets.size());
      std::print(myfile,
                 "# Guest–guest Coulomb (quadrupole) is not in the mean-field attraction: a cutoff 1/r "
                 "kernel is monopole-dominated for partial charges.\n");
    }
    else
    {
      std::print(myfile, "# WCA attraction stays on the COM density.\n");
    }
    std::print(myfile, "# ρ_n = (1/M) Σ_ω ρ(r, ω), M = {}.\n", this->numberOfOrientations);
    std::print(myfile, "# U(r, ω) is unclipped in both Henry and the Euler-Lagrange map; overlap poses are dropped.\n");
  }
  else
  {
    std::print(myfile, "# Density: spherical ρ(r). V_ext is the Helmholtz field, clipped at {:.3g} kT.\n",
               options.wellFloorInKT);
    if (wantsMolecular && !this->molecular)
    {
      std::print(myfile, "# WARNING: U(r, ω) did not fit in memory; fell back to spherical ρ on the Helmholtz field\n");
    }
  }
  std::print(myfile, "#\n");
  std::print(myfile, "# The grand potential is the standard NLDFT functional: fundamental-measure hard spheres\n");
  std::print(myfile, "# (White Bear; the uniform limit is Carnahan-Starling) plus a mean-field attraction,\n");
  std::print(myfile, "# in the framework's own dispersion+Ewald field --- no slit or cylinder wall model. The\n");
  std::print(myfile, "# Euler-Lagrange equation is solved by Picard iteration, all convolutions by FFT over the\n");
  std::print(myfile, "# periodic cell; each pressure is warm-started from the one below it, which keeps the\n");
  std::print(myfile, "# profile on the adsorption branch through condensation.\n");
  std::print(myfile, "#\n");
  if (this->molecular && !options.attractiveSiteOffsets.empty())
  {
    if (!options.dumbbellSiteOffsets.empty())
    {
      std::print(myfile,
                 "# Fluid: fused-dumbbell FMT d {:.3f} [Å] at {} sites; TraPPE N–N 12-6 from σ {:.3f} [Å], "
                 "epsilon/k {:.2f} [K], cutoff {:.3f} [Å]\n",
                 options.dumbbellSiteDiameter, options.dumbbellSiteOffsets.size(), options.attractiveSigma,
                 options.attractiveEpsilon, options.cutoff);
    }
    else
    {
      std::print(myfile,
                 "# Fluid: COM FMT d {:.3f} [Å]; TraPPE N–N 12-6 from σ {:.3f} [Å], epsilon/k {:.2f} [K] at {} sites, "
                 "cutoff {:.3f} [Å]\n",
                 nldftHardSphereDiameter, options.attractiveSigma, options.attractiveEpsilon,
                 options.attractiveSiteOffsets.size(), options.cutoff);
    }
  }
  else if (this->molecular && !options.dumbbellSiteOffsets.empty())
  {
    std::print(myfile, "# Fluid: WCA sigma {:.3f} [Å], epsilon/k {:.2f} [K]; dumbbell FMT d {:.3f} [Å] at {} sites, cutoff {:.3f} [Å]\n",
               nldftSigma, nldftEpsilonKelvin, options.dumbbellSiteDiameter, options.dumbbellSiteOffsets.size(),
               nldftCutoffFactor * nldftSigma);
  }
  else
  {
    std::print(myfile, "# Fluid: sigma {:.3f} [Å], epsilon/k {:.2f} [K], hard-sphere d {:.3f} [Å], cutoff {:.3f} [Å]\n",
               nldftSigma, nldftEpsilonKelvin, nldftHardSphereDiameter, nldftCutoffFactor * nldftSigma);
  }
  std::print(myfile, "# Model bulk at {} [K]: P0 {:.6g} [Pa] ({:.4g} atm), gas {:.6g} [1/Å³], liquid {:.6g} [1/Å³]\n",
             temperature, this->bulk.saturationPressure, this->bulk.saturationPressure / 101325.0,
             this->bulk.gasDensity, this->bulk.liquidDensity);
  if (this->bulk.experimentalSaturation)
  {
    std::print(myfile, "# The uniform fluid has no atmosphere-like spinodal at this T; x = P / P0 uses experimental N2 P0.\n");
  }
  else
  {
    std::print(myfile, "# x = P / P0 is measured against the model fluid's own P0, as kernel builders do\n");
  }
  std::print(myfile, "# Exact grid Henry occupancy at experimental P0: {:.6g} [molecules / cell]\n",
             this->henryMoleculesPerCell);
  {
    const double henryCoefficient =
        (this->henryMoleculesPerCell > 0.0) ? this->henryMoleculesPerCell / nitrogenSaturationPressure : 0.0;
    const double packingCapacity = framework.unitCell.volume / nitrogenLiquidVolume;
    std::print(myfile,
               "# Pressure window: x = {:.6g} -- {:.6g} (P = {:.6g} -- {:.6g} Pa), Henry-based like WHAM/TMMC\n",
               options.xLow, options.xHigh, options.xLow * this->bulk.saturationPressure,
               options.xHigh * this->bulk.saturationPressure);
    std::print(myfile,
               "#   K_H {:.6g} /cell/Pa, packing {:.4g} /cell, bottom where Henry holds {:g} of packing "
               "(floor {:.4g} Pa)\n",
               henryCoefficient, packingCapacity, nldftPressureWindowFillFraction, nldftPressureWindowFloorPa);
  }
  if (molecular && !orientationKelvin.empty())
  {
    const double henryCapped12 =
        (kT > 0.0) ? nitrogenSaturationPressure * volumeSI / kT * boltzmannCapped12 : 0.0;
    const double henryFromHoles =
        (kT > 0.0) ? nitrogenSaturationPressure * volumeSI / kT * boltzmannBelow40 : 0.0;
    const std::size_t nStates = energyKelvin.size() * this->numberOfOrientations;
    std::print(myfile, "# U(r, ω): {} / {} states at the overlap ceiling, {} voxels with a free orientation,\n",
               nCeiling, nStates, voxelsOpen);
    std::print(myfile, "# {} voxels with any attractive ω, {} / {} / {} / {} states below −12 / −20 / −40 / −80 kT\n",
               voxelsAttractive, nBelow12, nBelow20, nBelow40, nBelow80);
    std::print(myfile, "# Henry if wells deeper than 12 kT are clipped: {:.6g} / cell\n", henryCapped12);
    std::print(myfile, "# Henry from states below −40 kT alone: {:.6g} / cell ({:.3g} of the total)\n",
               henryFromHoles, (boltzmann > 0.0) ? boltzmannBelow40 / boltzmann : 0.0);
    auto describePose = [&](const DeepPose &pose)
    {
      if (!(pose.kelvin < 0.0)) return;
      const uint3 g = this->field.gridSize;
      const std::size_t nx = g.x;
      const std::size_t ny = g.y;
      const std::size_t ix = pose.voxel % nx;
      const std::size_t iy = (pose.voxel / nx) % ny;
      const std::size_t iz = pose.voxel / (nx * ny);
      const double3 frac(static_cast<double>(ix) / static_cast<double>(nx),
                         static_cast<double>(iy) / static_cast<double>(ny),
                         static_cast<double>(iz) / static_cast<double>(g.z));
      const double3 com = this->field.unitCell.cell * frac;
      const std::vector<double3> axes = orientationSet(this->numberOfOrientations, probe.headTailSymmetric);
      const double3 axis = axes[pose.orientation];
      std::print(myfile, "# Pose U={:.5g} K ({:.3g} kT) voxel ({}, {}, {}) ω={}  well-distance={:.3f} Å\n",
                 pose.kelvin, pose.kelvin / temperature, ix, iy, iz, pose.orientation,
                 static_cast<double>(this->field.distance[pose.voxel]));
      std::print(myfile, "#   COM ({:.3f}, {:.3f}, {:.3f}) Å  axis ({:.3f}, {:.3f}, {:.3f})\n", com.x, com.y, com.z,
                 axis.x, axis.y, axis.z);
      for (const LinearProbe::Site &site : probe.sites)
      {
        const double3 pos = com + site.offset * axis;
        double bestR = 1.0e10;
        std::size_t bestAtom = 0;
        for (std::size_t a = 0; a < framework.atoms.size(); ++a)
        {
          const double3 dr =
              this->field.unitCell.applyPeriodicBoundaryConditions(pos - framework.atoms[a].position);
          const double r = dr.length();
          if (r < bestR)
          {
            bestR = r;
            bestAtom = a;
          }
        }
        const CrystalAtom &atom = framework.atoms[bestAtom];
        const PairParameters &pair = interactions(site.type, atom.type);
        const double rmin = wellContactPrefactor * pair.sizeParameter;
        const std::string typeName =
            atom.type < interactions.names.size() ? interactions.names[atom.type] : "?";
        std::print(myfile,
                   "#   {:<6} q={:+.3f}  r={:.3f} Å to {} q={:+.3f}  σ={:.3f} rmin={:.3f} clearance={:+.3f} Å\n",
                   site.name, site.charge, bestR, typeName, atom.charge, pair.sizeParameter, rmin, bestR - rmin);
      }
    };
    for (const DeepPose &pose : deepest) describePose(pose);
  }
  if (fieldStats.clippedVoxels > 0)
  {
    std::print(myfile, "# V_ext: {} voxels clipped to the well floor of {:.4g} K ({:.3g} kT)\n",
               fieldStats.clippedVoxels, fieldStats.wellFloorKelvin, options.wellFloorInKT);
  }
  if (this->unconverged > 0)
  {
    std::print(myfile, "# WARNING: {} pressure points ran out of iterations before converging\n", this->unconverged);
  }
  std::print(myfile, "# Timing: {} [s] for the field and the {} pressure points\n", this->seconds,
             this->isotherm.size());
  std::print(myfile, "\n");

  std::print(myfile, "# The Rouquerol/BET read of this isotherm, fitted exactly as it would be to a measured or\n");
  std::print(myfile, "# GCMC isotherm.\n");
  std::print(myfile, "{:34} {:14.4f} [m²/g] {:12.4f} [m²/cm³]\n",
             "BET area:", this->fit.gravimetricArea, this->fit.volumetricArea);
  std::print(myfile, "{:34} {:14.4f} [molecules / cell]\n", "Monolayer capacity n_m:", this->fit.monolayerCapacity);
  std::print(myfile, "{:34} {:14.4f} [-]\n", "BET C constant:", this->fit.cConstant);
  std::print(myfile, "# Fit window: {:.4g} -- {:.4g} in P/P0, r² = {:.5f}\n", this->fit.windowLow,
             this->fit.windowHigh, this->fit.rSquared);
  writeIsothermDerivedProperties(myfile, this->fit);
  std::print(myfile, "\n");

  std::print(myfile, "# {:>12} {:>16} {:>16} {:>16}\n", "x = P/P0", "n [molec/cell]", "n(1-x)", "x/(n(1-x))");
  for (const IsothermPoint &point : this->isotherm)
  {
    double g = point.moleculesPerCell * (1.0 - point.relativePressure);
    double y = (g > 0.0) ? point.relativePressure / g : 0.0;
    std::print(myfile, "  {:>12.6g} {:>16.8g} {:>16.8g} {:>16.8g}\n", point.relativePressure, point.moleculesPerCell,
               g, y);
  }
  myfile.close();
}
