module;

module system;

import std;

import randomnumbers;
import int3;
import uint3;
import double3;
import double3x3;
import atom;
import framework;
import forcefield;
import simulationbox;
import units;
import interpolation_energy_grid;
import interactions_framework_molecule;
import interactions_framework_molecule_grid;
import interactions_intermolecular;
import interactions_ewald;
import interactions_external_field;
import interactions_external_field_grid;

// System interpolation grids: external-field and framework energy grid setup/validation.

void System::createExternalFieldInterpolationGrid(std::ostream& stream, std::size_t systemId)
{
  // use local random-number generator (so that it does not interfere with a binary-restart)
  RandomNumber random{std::nullopt};

  std::size_t numberOfGridTestPoints = forceField.numberOfGridTestPoints;

  if(hasExternalField)
  {
    if(forceField.useExternalFieldGrid)
    {
      // int3 numberOfExternalFieldGridPoints  = forceField.numberOfExternalFieldGridPoints;
      uint3 numberOfExternalFieldGridPoints = InterpolationEnergyGrid::parseExternalFieldGridDimensions(forceField.externalFieldGridFileName);

      externalFieldInterpolationGrid = InterpolationEnergyGrid(simulationBox, forceField.potentialEnergySurfaceOrigin,
                                                      numberOfExternalFieldGridPoints, forceField.interpolationScheme);

      std::print(stream, "Generating an external field interpolation grid ({}x{}x{})\n",
                 externalFieldInterpolationGrid->numberOfGridPoints.x, 
                 externalFieldInterpolationGrid->numberOfGridPoints.y, 
                 externalFieldInterpolationGrid->numberOfGridPoints.z);
      std::print(stream, "===============================================================================\n");
      externalFieldInterpolationGrid->makeExternalFieldInterpolationGrid(stream, forceField, simulationBox);

      double count{};
      double summed_errors{};
      double boltzmann_weighted_difference_squared_summed{};
      double boltzmann_weighted_full_squared_summed{};

      for (std::size_t i = 0; i < numberOfGridTestPoints; ++i)
      {
        // generate random position in super cell
        double3 s = double3(random.uniform(), random.uniform(), random.uniform());
        double3 pos = simulationBox.cell * s;

        auto [interpolated_energy, interpolated_gradient, interpolated_hessian] =
            externalFieldInterpolationGrid->interpolateHessian(pos);

        // convert to Kelvin
        interpolated_energy *= Units::EnergyToKelvin;

        interpolated_gradient.x *= Units::EnergyToKelvin;
        interpolated_gradient.y *= Units::EnergyToKelvin;
        interpolated_gradient.z *= Units::EnergyToKelvin;

        interpolated_hessian.ax *= Units::EnergyToKelvin;
        interpolated_hessian.ay *= Units::EnergyToKelvin;
        interpolated_hessian.az *= Units::EnergyToKelvin;
        interpolated_hessian.bx *= Units::EnergyToKelvin;
        interpolated_hessian.by *= Units::EnergyToKelvin;
        interpolated_hessian.bz *= Units::EnergyToKelvin;
        interpolated_hessian.cx *= Units::EnergyToKelvin;
        interpolated_hessian.cy *= Units::EnergyToKelvin;
        interpolated_hessian.cz *= Units::EnergyToKelvin;

        double analytical_energy{};
        double3 analytical_gradient{};
        double3x3 analytical_hessian{};
        switch (forceField.interpolationScheme)
        {
          case ForceField::InterpolationScheme::Polynomial:
          {
            std::array<double, 8> analytical_polynomial = Interactions::calculateTricubicFractionalAtPositionExternalField(
                 forceField, simulationBox, pos + forceField.potentialEnergySurfaceOrigin);

            analytical_energy = analytical_polynomial[0] * Units::EnergyToKelvin;
          }
          break;
          case ForceField::InterpolationScheme::Tricubic:
          {
            std::array<double, 8> analytical_tricubic = Interactions::calculateTricubicFractionalAtPositionExternalField(
                 forceField, simulationBox, pos + forceField.potentialEnergySurfaceOrigin);

            analytical_energy = analytical_tricubic[0] * Units::EnergyToKelvin;

            // convert gradient from fractional to Cartesian
            analytical_gradient =
                framework->simulationBox.inverseCell.transpose() *
                double3(analytical_tricubic[1], analytical_tricubic[2], analytical_tricubic[3]);

            analytical_gradient.x *= Units::EnergyToKelvin;
            analytical_gradient.y *= Units::EnergyToKelvin;
            analytical_gradient.z *= Units::EnergyToKelvin;
          }
          break;
          case ForceField::InterpolationScheme::Triquintic:
          {
            std::array<double, 27> analytical_triquintic = Interactions::calculateTriquinticFractionalAtPositionExternalField(
                 forceField, simulationBox, pos + forceField.potentialEnergySurfaceOrigin);

            analytical_energy = analytical_triquintic[0] * Units::EnergyToKelvin;

            // convert gradient from fractional to Cartesian
            analytical_gradient =
                framework->simulationBox.inverseCell.transpose() *
                double3(analytical_triquintic[1], analytical_triquintic[2], analytical_triquintic[3]);

            analytical_gradient.x *= Units::EnergyToKelvin;
            analytical_gradient.y *= Units::EnergyToKelvin;
            analytical_gradient.z *= Units::EnergyToKelvin;

            double3x3 hessian =
                double3x3(analytical_triquintic[4], analytical_triquintic[5], analytical_triquintic[6],
                          analytical_triquintic[5], analytical_triquintic[7], analytical_triquintic[8],
                          analytical_triquintic[6], analytical_triquintic[8], analytical_triquintic[9]);
            analytical_hessian =
                framework->simulationBox.inverseCell.transpose() * hessian * framework->simulationBox.inverseCell;

            analytical_hessian.ax *= Units::EnergyToKelvin;
            analytical_hessian.ay *= Units::EnergyToKelvin;
            analytical_hessian.az *= Units::EnergyToKelvin;
            analytical_hessian.bx *= Units::EnergyToKelvin;
            analytical_hessian.by *= Units::EnergyToKelvin;
            analytical_hessian.bz *= Units::EnergyToKelvin;
            analytical_hessian.cx *= Units::EnergyToKelvin;
            analytical_hessian.cy *= Units::EnergyToKelvin;
            analytical_hessian.cz *= Units::EnergyToKelvin;
          }
          break;
        }

        double boltzmann_weight_value = std::exp(-beta * analytical_energy);
        double difference = analytical_energy - interpolated_energy;

        summed_errors += std::fabs(difference);
        count += 1.0;

        boltzmann_weighted_difference_squared_summed += difference * difference * boltzmann_weight_value;
        boltzmann_weighted_full_squared_summed += boltzmann_weight_value;
      }

      std::print(stream, "Testing external-field interpolation grid ({}x{}x{})\n", numberOfExternalFieldGridPoints.x,
                 numberOfExternalFieldGridPoints.y, numberOfExternalFieldGridPoints.z);
      std::print(stream, "-------------------------------------------------------------------------------\n");
      std::print(stream, "(Using {} points for testing)\n\n", numberOfGridTestPoints);

      std::print(stream, "Absolute error energy:                  {}\n\n", summed_errors / count);
      std::print(stream, "Boltzmann weighted error energy:        {}\n\n",
                        std::sqrt(boltzmann_weighted_difference_squared_summed / 
                                  boltzmann_weighted_full_squared_summed));
    }
  }

  if(forceField.writeExternalFieldInterpolationGrid)
  {
    externalFieldInterpolationGrid->writeOutput(systemId, simulationBox, forceField);
  }
}

void System::createFrameworkInterpolationGrids(std::ostream& stream)
{
  // use local random-number generator (so that it does not interfere with a binary-restart)
  RandomNumber random{std::nullopt};

  std::size_t numberOfGridTestPoints = forceField.numberOfGridTestPoints;

  if (framework.has_value())
  {
    uint3 numberOfCoulombGridPoints{};
    if (forceField.numberOfVDWGridPoints.has_value())
    {
      numberOfCoulombGridPoints = forceField.numberOfCoulombGridPoints.value();
    }
    else
    {
      const double3 perpendicular_widths = framework->simulationBox.perpendicularWidths();
      numberOfCoulombGridPoints.x =
          static_cast<std::size_t>(perpendicular_widths.x / forceField.spacingCoulombGrid + 0.5);
      numberOfCoulombGridPoints.y =
          static_cast<std::size_t>(perpendicular_widths.y / forceField.spacingCoulombGrid + 0.5);
      numberOfCoulombGridPoints.z =
          static_cast<std::size_t>(perpendicular_widths.z / forceField.spacingCoulombGrid + 0.5);
    }

    // also create a Charge grid when needed
    if (!forceField.gridPseudoAtomIndices.empty())
    {
      std::print(stream, "Generating an Ewald Real interpolation grid ({}x{}x{}) for a unit charge\n",
                 numberOfCoulombGridPoints.x, numberOfCoulombGridPoints.y, numberOfCoulombGridPoints.z);
      std::print(stream, "===============================================================================\n");

      interpolationGrids.back() =
          InterpolationEnergyGrid(framework->simulationBox, forceField.potentialEnergySurfaceOrigin,
                                  numberOfCoulombGridPoints, forceField.interpolationScheme);
      interpolationGrids.back()->makeFrameworkInterpolationGrid(stream, ForceField::InterpolationGridType::EwaldReal, forceField,
                                                       framework.value(), forceField.cutOffCoulomb, 0);
    }

    uint3 numberOfVDWGridPoints{};
    if (forceField.numberOfVDWGridPoints.has_value())
    {
      numberOfVDWGridPoints = forceField.numberOfVDWGridPoints.value();
    }
    else
    {
      const double3 perpendicular_widths = framework->simulationBox.perpendicularWidths();
      numberOfVDWGridPoints.x = static_cast<std::size_t>(perpendicular_widths.x / forceField.spacingVDWGrid + 0.5);
      numberOfVDWGridPoints.y = static_cast<std::size_t>(perpendicular_widths.y / forceField.spacingVDWGrid + 0.5);
      numberOfVDWGridPoints.z = static_cast<std::size_t>(perpendicular_widths.z / forceField.spacingVDWGrid + 0.5);
    }

    for (const std::size_t& index : forceField.gridPseudoAtomIndices)
    {
      std::print(stream, "Generating an VDW interpolation grid ({}x{}x{}) for {}\n", numberOfVDWGridPoints.x,
                 numberOfVDWGridPoints.y, numberOfVDWGridPoints.z, forceField.pseudoAtoms[index].name);
      std::print(stream, "===============================================================================\n");

      interpolationGrids[index] =
          InterpolationEnergyGrid(framework->simulationBox, forceField.potentialEnergySurfaceOrigin, numberOfVDWGridPoints, forceField.interpolationScheme);
      interpolationGrids[index]->makeFrameworkInterpolationGrid(stream, ForceField::InterpolationGridType::LennardJones,
                                                       forceField, framework.value(), forceField.cutOffFrameworkVDW,
                                                       index);

      double boltzmann_weight_summed_vdw{};
      double boltzmann_weighted_energy_full_summed_vdw{};
      double boltzmann_weighted_energy_interpolated_summed_vdw{};
      double boltzmann_weighted_difference_squared_summed_vdw{};
      double boltzmann_weighted_full_squared_summed_vdw{};

      double3 boltzmann_weighted_gradient_full_summed_vdw{};
      double3 boltzmann_weighted_gradient_interpolated_summed_vdw{};
      double3 boltzmann_weighted_difference_squared_summed_vdw_gradient{};
      double3 boltzmann_weighted_full_squared_summed_vdw_gradient{};

      double3x3 boltzmann_weighted_hessian_full_summed_vdw{};
      double3x3 boltzmann_weighted_hessian_interpolated_summed_vdw{};
      double3x3 boltzmann_weighted_difference_squared_summed_vdw_hessian{};
      double3x3 boltzmann_weighted_full_squared_summed_vdw_hessian{};

      double boltzmann_weight_summed_real_ewald{};
      double boltzmann_weighted_energy_full_summed_real_ewald{};
      double boltzmann_weighted_energy_interpolated_summed_real_ewald{};
      double boltzmann_weighted_difference_squared_summed_real_ewald{};
      double boltzmann_weighted_full_squared_summed_real_ewald{};

      double3 boltzmann_weighted_gradient_full_summed_real_ewald{};
      double3 boltzmann_weighted_gradient_interpolated_summed_real_ewald{};
      double3 boltzmann_weighted_difference_squared_summed_real_ewald_gradient{};
      double3 boltzmann_weighted_full_squared_summed_real_ewald_gradient{};

      double3x3 boltzmann_weighted_hessian_full_summed_real_ewald{};
      double3x3 boltzmann_weighted_hessian_interpolated_summed_real_ewald{};
      double3x3 boltzmann_weighted_difference_squared_summed_real_ewald_hessian{};
      double3x3 boltzmann_weighted_full_squared_summed_real_ewald_hessian{};

      for (std::size_t i = 0; i < numberOfGridTestPoints; ++i)
      {
        // generate random position in super cell
        double3 s = double3(random.uniform(), random.uniform(), random.uniform());
        double3 pos = simulationBox.cell * s;

        auto [interpolated_energy_vdw, interpolated_gradient_vdw, interpolated_hessian_vdw] =
            interpolationGrids[index]->interpolateHessian(pos);

        // convert to Kelvin
        interpolated_energy_vdw *= Units::EnergyToKelvin;
        interpolated_gradient_vdw.x *= Units::EnergyToKelvin;
        interpolated_gradient_vdw.y *= Units::EnergyToKelvin;
        interpolated_gradient_vdw.z *= Units::EnergyToKelvin;
        interpolated_hessian_vdw.ax *= Units::EnergyToKelvin;
        interpolated_hessian_vdw.ay *= Units::EnergyToKelvin;
        interpolated_hessian_vdw.az *= Units::EnergyToKelvin;
        interpolated_hessian_vdw.bx *= Units::EnergyToKelvin;
        interpolated_hessian_vdw.by *= Units::EnergyToKelvin;
        interpolated_hessian_vdw.bz *= Units::EnergyToKelvin;
        interpolated_hessian_vdw.cx *= Units::EnergyToKelvin;
        interpolated_hessian_vdw.cy *= Units::EnergyToKelvin;
        interpolated_hessian_vdw.cz *= Units::EnergyToKelvin;

        double analytical_vdw_energy{};
        double3 analytical_vdw_gradient{};
        double3x3 analytical_vdw_hessian{};
        switch (forceField.interpolationScheme)
        {
          case ForceField::InterpolationScheme::Polynomial:
          {
            double analytical_vdw_quintic =
                Interactions::calculateEnergyAtPosition(ForceField::InterpolationGridType::LennardJones, forceField,
                                                        simulationBox, pos, index, spanOfFrameworkAtoms());

            analytical_vdw_energy = analytical_vdw_quintic * Units::EnergyToKelvin;
          }
          break;
          case ForceField::InterpolationScheme::Tricubic:
          {
            std::array<double, 8> analytical_vdw_tricubic = Interactions::calculateTricubicFractionalAtPosition(
                ForceField::InterpolationGridType::LennardJones, forceField, simulationBox, pos, index,
                framework->simulationBox, spanOfFrameworkAtoms());

            analytical_vdw_energy = analytical_vdw_tricubic[0] * Units::EnergyToKelvin;

            // convert gradient from fractional to Cartesian
            analytical_vdw_gradient =
                framework->simulationBox.inverseCell.transpose() *
                double3(analytical_vdw_tricubic[1], analytical_vdw_tricubic[2], analytical_vdw_tricubic[3]);

            analytical_vdw_gradient.x *= Units::EnergyToKelvin;
            analytical_vdw_gradient.y *= Units::EnergyToKelvin;
            analytical_vdw_gradient.z *= Units::EnergyToKelvin;
          }
          break;
          case ForceField::InterpolationScheme::Triquintic:
          {
            std::array<double, 27> analytical_vdw_triquintic = Interactions::calculateTriquinticFractionalAtPosition(
                ForceField::InterpolationGridType::LennardJones, forceField, simulationBox, pos, index,
                framework->simulationBox, spanOfFrameworkAtoms());

            analytical_vdw_energy = analytical_vdw_triquintic[0] * Units::EnergyToKelvin;

            // convert gradient from fractional to Cartesian
            analytical_vdw_gradient =
                framework->simulationBox.inverseCell.transpose() *
                double3(analytical_vdw_triquintic[1], analytical_vdw_triquintic[2], analytical_vdw_triquintic[3]);

            analytical_vdw_gradient.x *= Units::EnergyToKelvin;
            analytical_vdw_gradient.y *= Units::EnergyToKelvin;
            analytical_vdw_gradient.z *= Units::EnergyToKelvin;

            double3x3 hessian =
                double3x3(analytical_vdw_triquintic[4], analytical_vdw_triquintic[5], analytical_vdw_triquintic[6],
                          analytical_vdw_triquintic[5], analytical_vdw_triquintic[7], analytical_vdw_triquintic[8],
                          analytical_vdw_triquintic[6], analytical_vdw_triquintic[8], analytical_vdw_triquintic[9]);
            analytical_vdw_hessian =
                framework->simulationBox.inverseCell.transpose() * hessian * framework->simulationBox.inverseCell;

            analytical_vdw_hessian.ax *= Units::EnergyToKelvin;
            analytical_vdw_hessian.ay *= Units::EnergyToKelvin;
            analytical_vdw_hessian.az *= Units::EnergyToKelvin;
            analytical_vdw_hessian.bx *= Units::EnergyToKelvin;
            analytical_vdw_hessian.by *= Units::EnergyToKelvin;
            analytical_vdw_hessian.bz *= Units::EnergyToKelvin;
            analytical_vdw_hessian.cx *= Units::EnergyToKelvin;
            analytical_vdw_hessian.cy *= Units::EnergyToKelvin;
            analytical_vdw_hessian.cz *= Units::EnergyToKelvin;
          }
          break;
        }

        double boltzmann_weight_vdw = std::exp(-beta * analytical_vdw_energy);

        boltzmann_weight_summed_vdw += boltzmann_weight_vdw;

        boltzmann_weighted_energy_interpolated_summed_vdw += interpolated_energy_vdw * boltzmann_weight_vdw;
        boltzmann_weighted_energy_full_summed_vdw += analytical_vdw_energy * boltzmann_weight_vdw;
        boltzmann_weighted_difference_squared_summed_vdw += (analytical_vdw_energy - interpolated_energy_vdw) *
                                                            (analytical_vdw_energy - interpolated_energy_vdw) *
                                                            boltzmann_weight_vdw;
        boltzmann_weighted_full_squared_summed_vdw +=
            analytical_vdw_energy * analytical_vdw_energy * boltzmann_weight_vdw;

        // gradient
        boltzmann_weighted_gradient_interpolated_summed_vdw.x += interpolated_gradient_vdw.x * boltzmann_weight_vdw;
        boltzmann_weighted_gradient_interpolated_summed_vdw.y += interpolated_gradient_vdw.y * boltzmann_weight_vdw;
        boltzmann_weighted_gradient_interpolated_summed_vdw.z += interpolated_gradient_vdw.z * boltzmann_weight_vdw;

        boltzmann_weighted_gradient_full_summed_vdw.x += analytical_vdw_gradient.x * boltzmann_weight_vdw;
        boltzmann_weighted_gradient_full_summed_vdw.y += analytical_vdw_gradient.y * boltzmann_weight_vdw;
        boltzmann_weighted_gradient_full_summed_vdw.z += analytical_vdw_gradient.z * boltzmann_weight_vdw;

        boltzmann_weighted_difference_squared_summed_vdw_gradient.x +=
            (analytical_vdw_gradient.x - interpolated_gradient_vdw.x) *
            (analytical_vdw_gradient.x - interpolated_gradient_vdw.x) * boltzmann_weight_vdw;
        boltzmann_weighted_difference_squared_summed_vdw_gradient.y +=
            (analytical_vdw_gradient.y - interpolated_gradient_vdw.y) *
            (analytical_vdw_gradient.y - interpolated_gradient_vdw.y) * boltzmann_weight_vdw;
        boltzmann_weighted_difference_squared_summed_vdw_gradient.z +=
            (analytical_vdw_gradient.z - interpolated_gradient_vdw.z) *
            (analytical_vdw_gradient.z - interpolated_gradient_vdw.z) * boltzmann_weight_vdw;

        boltzmann_weighted_full_squared_summed_vdw_gradient.x +=
            analytical_vdw_gradient.x * analytical_vdw_gradient.x * boltzmann_weight_vdw;
        boltzmann_weighted_full_squared_summed_vdw_gradient.y +=
            analytical_vdw_gradient.y * analytical_vdw_gradient.y * boltzmann_weight_vdw;
        boltzmann_weighted_full_squared_summed_vdw_gradient.z +=
            analytical_vdw_gradient.z * analytical_vdw_gradient.z * boltzmann_weight_vdw;

        // Hessian
        boltzmann_weighted_hessian_interpolated_summed_vdw.ax += interpolated_hessian_vdw.ax * boltzmann_weight_vdw;
        boltzmann_weighted_hessian_interpolated_summed_vdw.ay += interpolated_hessian_vdw.ay * boltzmann_weight_vdw;
        boltzmann_weighted_hessian_interpolated_summed_vdw.az += interpolated_hessian_vdw.az * boltzmann_weight_vdw;
        boltzmann_weighted_hessian_interpolated_summed_vdw.by += interpolated_hessian_vdw.by * boltzmann_weight_vdw;
        boltzmann_weighted_hessian_interpolated_summed_vdw.bz += interpolated_hessian_vdw.bz * boltzmann_weight_vdw;
        boltzmann_weighted_hessian_interpolated_summed_vdw.cz += interpolated_hessian_vdw.cz * boltzmann_weight_vdw;

        boltzmann_weighted_hessian_full_summed_vdw.ax += analytical_vdw_hessian.ax * boltzmann_weight_vdw;
        boltzmann_weighted_hessian_full_summed_vdw.ay += analytical_vdw_hessian.ay * boltzmann_weight_vdw;
        boltzmann_weighted_hessian_full_summed_vdw.az += analytical_vdw_hessian.az * boltzmann_weight_vdw;
        boltzmann_weighted_hessian_full_summed_vdw.by += analytical_vdw_hessian.by * boltzmann_weight_vdw;
        boltzmann_weighted_hessian_full_summed_vdw.bz += analytical_vdw_hessian.bz * boltzmann_weight_vdw;
        boltzmann_weighted_hessian_full_summed_vdw.cz += analytical_vdw_hessian.cz * boltzmann_weight_vdw;

        boltzmann_weighted_difference_squared_summed_vdw_hessian.ax +=
            (analytical_vdw_hessian.ax - interpolated_hessian_vdw.ax) *
            (analytical_vdw_hessian.ax - interpolated_hessian_vdw.ax) * boltzmann_weight_vdw;
        boltzmann_weighted_difference_squared_summed_vdw_hessian.ay +=
            (analytical_vdw_hessian.ay - interpolated_hessian_vdw.ay) *
            (analytical_vdw_hessian.ay - interpolated_hessian_vdw.ay) * boltzmann_weight_vdw;
        boltzmann_weighted_difference_squared_summed_vdw_hessian.az +=
            (analytical_vdw_hessian.az - interpolated_hessian_vdw.az) *
            (analytical_vdw_hessian.az - interpolated_hessian_vdw.az) * boltzmann_weight_vdw;
        boltzmann_weighted_difference_squared_summed_vdw_hessian.by +=
            (analytical_vdw_hessian.by - interpolated_hessian_vdw.by) *
            (analytical_vdw_hessian.by - interpolated_hessian_vdw.by) * boltzmann_weight_vdw;
        boltzmann_weighted_difference_squared_summed_vdw_hessian.bz +=
            (analytical_vdw_hessian.bz - interpolated_hessian_vdw.bz) *
            (analytical_vdw_hessian.bz - interpolated_hessian_vdw.bz) * boltzmann_weight_vdw;
        boltzmann_weighted_difference_squared_summed_vdw_hessian.cz +=
            (analytical_vdw_hessian.cz - interpolated_hessian_vdw.cz) *
            (analytical_vdw_hessian.cz - interpolated_hessian_vdw.cz) * boltzmann_weight_vdw;

        boltzmann_weighted_full_squared_summed_vdw_hessian.ax +=
            analytical_vdw_hessian.ax * analytical_vdw_hessian.ax * boltzmann_weight_vdw;
        boltzmann_weighted_full_squared_summed_vdw_hessian.ay +=
            analytical_vdw_hessian.ay * analytical_vdw_hessian.ay * boltzmann_weight_vdw;
        boltzmann_weighted_full_squared_summed_vdw_hessian.az +=
            analytical_vdw_hessian.az * analytical_vdw_hessian.az * boltzmann_weight_vdw;
        boltzmann_weighted_full_squared_summed_vdw_hessian.by +=
            analytical_vdw_hessian.by * analytical_vdw_hessian.by * boltzmann_weight_vdw;
        boltzmann_weighted_full_squared_summed_vdw_hessian.bz +=
            analytical_vdw_hessian.bz * analytical_vdw_hessian.bz * boltzmann_weight_vdw;
        boltzmann_weighted_full_squared_summed_vdw_hessian.cz +=
            analytical_vdw_hessian.cz * analytical_vdw_hessian.cz * boltzmann_weight_vdw;

        // test charges when no VDW overlap detected (putting a unit charge on top of negative framework atom can lead
        // to very negative energies)
        if (analytical_vdw_energy < forceField.energyOverlapCriteria)
        {
          double charge = forceField.pseudoAtoms[index].charge;
          auto [interpolated_value_real_ewald, interpolated_gradient_real_ewald, interpolated_hessian_real_ewald] =
              interpolationGrids.back()->interpolateHessian(pos);

          interpolated_value_real_ewald *= (charge * Units::EnergyToKelvin);
          interpolated_gradient_real_ewald.x *= (charge * Units::EnergyToKelvin);
          interpolated_gradient_real_ewald.y *= (charge * Units::EnergyToKelvin);
          interpolated_gradient_real_ewald.z *= (charge * Units::EnergyToKelvin);
          interpolated_hessian_real_ewald.ax *= (charge * Units::EnergyToKelvin);
          interpolated_hessian_real_ewald.ay *= (charge * Units::EnergyToKelvin);
          interpolated_hessian_real_ewald.az *= (charge * Units::EnergyToKelvin);
          interpolated_hessian_real_ewald.bx *= (charge * Units::EnergyToKelvin);
          interpolated_hessian_real_ewald.by *= (charge * Units::EnergyToKelvin);
          interpolated_hessian_real_ewald.bz *= (charge * Units::EnergyToKelvin);
          interpolated_hessian_real_ewald.cx *= (charge * Units::EnergyToKelvin);
          interpolated_hessian_real_ewald.cy *= (charge * Units::EnergyToKelvin);
          interpolated_hessian_real_ewald.cz *= (charge * Units::EnergyToKelvin);

          double analytical_real_ewald_energy{};
          double3 analytical_real_ewald_gradient{};
          double3x3 analytical_real_ewald_hessian{};
          switch (forceField.interpolationScheme)
          {
            case ForceField::InterpolationScheme::Polynomial:
            {
              double analytical_real_ewald_tricubic =
                  Interactions::calculateEnergyAtPosition(ForceField::InterpolationGridType::EwaldReal, forceField,
                                                          simulationBox, pos, index, spanOfFrameworkAtoms());

              analytical_real_ewald_energy = charge * analytical_real_ewald_tricubic * Units::EnergyToKelvin;
            }
            break;
            case ForceField::InterpolationScheme::Tricubic:
            {
              std::array<double, 8> analytical_real_ewald_tricubic =
                  Interactions::calculateTricubicFractionalAtPosition(ForceField::InterpolationGridType::EwaldReal,
                                                                      forceField, simulationBox, pos, index,
                                                                      framework->simulationBox, spanOfFrameworkAtoms());

              analytical_real_ewald_energy = charge * analytical_real_ewald_tricubic[0] * Units::EnergyToKelvin;

              // convert gradient from fractional to Cartesian
              analytical_real_ewald_gradient =
                  framework->simulationBox.inverseCell.transpose() * double3(analytical_real_ewald_tricubic[1],
                                                                             analytical_real_ewald_tricubic[2],
                                                                             analytical_real_ewald_tricubic[3]);
              analytical_real_ewald_gradient.x *= (charge * Units::EnergyToKelvin);
              analytical_real_ewald_gradient.y *= (charge * Units::EnergyToKelvin);
              analytical_real_ewald_gradient.z *= (charge * Units::EnergyToKelvin);
            }
            break;
            case ForceField::InterpolationScheme::Triquintic:
            {
              std::array<double, 27> analytical_real_ewald_triquintic =
                  Interactions::calculateTriquinticFractionalAtPosition(
                      ForceField::InterpolationGridType::EwaldReal, forceField, simulationBox, pos, index,
                      framework->simulationBox, spanOfFrameworkAtoms());

              analytical_real_ewald_energy = charge * analytical_real_ewald_triquintic[0] * Units::EnergyToKelvin;

              // convert gradient from fractional to Cartesian
              analytical_real_ewald_gradient =
                  framework->simulationBox.inverseCell.transpose() * double3(analytical_real_ewald_triquintic[1],
                                                                             analytical_real_ewald_triquintic[2],
                                                                             analytical_real_ewald_triquintic[3]);

              analytical_real_ewald_gradient.x *= (charge * Units::EnergyToKelvin);
              analytical_real_ewald_gradient.y *= (charge * Units::EnergyToKelvin);
              analytical_real_ewald_gradient.z *= (charge * Units::EnergyToKelvin);

              double3x3 hessian = double3x3(analytical_real_ewald_triquintic[4], analytical_real_ewald_triquintic[5],
                                            analytical_real_ewald_triquintic[6], analytical_real_ewald_triquintic[5],
                                            analytical_real_ewald_triquintic[7], analytical_real_ewald_triquintic[8],
                                            analytical_real_ewald_triquintic[6], analytical_real_ewald_triquintic[8],
                                            analytical_real_ewald_triquintic[9]);
              analytical_real_ewald_hessian =
                  framework->simulationBox.inverseCell.transpose() * hessian * framework->simulationBox.inverseCell;

              analytical_real_ewald_hessian.ax *= (charge * Units::EnergyToKelvin);
              analytical_real_ewald_hessian.ay *= (charge * Units::EnergyToKelvin);
              analytical_real_ewald_hessian.az *= (charge * Units::EnergyToKelvin);
              analytical_real_ewald_hessian.bx *= (charge * Units::EnergyToKelvin);
              analytical_real_ewald_hessian.by *= (charge * Units::EnergyToKelvin);
              analytical_real_ewald_hessian.bz *= (charge * Units::EnergyToKelvin);
              analytical_real_ewald_hessian.cx *= (charge * Units::EnergyToKelvin);
              analytical_real_ewald_hessian.cy *= (charge * Units::EnergyToKelvin);
              analytical_real_ewald_hessian.cz *= (charge * Units::EnergyToKelvin);
            }
            break;
          }

          double boltzmann_weight_real_ewald = std::exp(-beta * analytical_vdw_energy);

          boltzmann_weight_summed_real_ewald += boltzmann_weight_real_ewald;

          // energy
          boltzmann_weighted_energy_interpolated_summed_real_ewald +=
              interpolated_value_real_ewald * boltzmann_weight_real_ewald;
          boltzmann_weighted_energy_full_summed_real_ewald +=
              analytical_real_ewald_energy * boltzmann_weight_real_ewald;
          boltzmann_weighted_difference_squared_summed_real_ewald +=
              (analytical_real_ewald_energy - interpolated_value_real_ewald) *
              (analytical_real_ewald_energy - interpolated_value_real_ewald) * boltzmann_weight_real_ewald;
          boltzmann_weighted_full_squared_summed_real_ewald +=
              analytical_real_ewald_energy * analytical_real_ewald_energy * boltzmann_weight_real_ewald;

          // gradients
          boltzmann_weighted_gradient_interpolated_summed_real_ewald.x +=
              interpolated_gradient_real_ewald.x * boltzmann_weight_real_ewald;
          boltzmann_weighted_gradient_interpolated_summed_real_ewald.y +=
              interpolated_gradient_real_ewald.y * boltzmann_weight_real_ewald;
          boltzmann_weighted_gradient_interpolated_summed_real_ewald.z +=
              interpolated_gradient_real_ewald.z * boltzmann_weight_real_ewald;

          boltzmann_weighted_gradient_full_summed_real_ewald.x +=
              analytical_real_ewald_gradient.x * boltzmann_weight_real_ewald;
          boltzmann_weighted_gradient_full_summed_real_ewald.y +=
              analytical_real_ewald_gradient.y * boltzmann_weight_real_ewald;
          boltzmann_weighted_gradient_full_summed_real_ewald.z +=
              analytical_real_ewald_gradient.z * boltzmann_weight_real_ewald;

          boltzmann_weighted_difference_squared_summed_real_ewald_gradient.x +=
              (analytical_real_ewald_gradient.x - interpolated_gradient_real_ewald.x) *
              (analytical_real_ewald_gradient.x - interpolated_gradient_real_ewald.x) * boltzmann_weight_real_ewald;
          boltzmann_weighted_difference_squared_summed_real_ewald_gradient.y +=
              (analytical_real_ewald_gradient.y - interpolated_gradient_real_ewald.y) *
              (analytical_real_ewald_gradient.y - interpolated_gradient_real_ewald.y) * boltzmann_weight_real_ewald;
          boltzmann_weighted_difference_squared_summed_real_ewald_gradient.z +=
              (analytical_real_ewald_gradient.z - interpolated_gradient_real_ewald.z) *
              (analytical_real_ewald_gradient.z - interpolated_gradient_real_ewald.z) * boltzmann_weight_real_ewald;

          boltzmann_weighted_full_squared_summed_real_ewald_gradient.x +=
              analytical_real_ewald_gradient.x * analytical_real_ewald_gradient.x * boltzmann_weight_real_ewald;
          boltzmann_weighted_full_squared_summed_real_ewald_gradient.y +=
              analytical_real_ewald_gradient.y * analytical_real_ewald_gradient.y * boltzmann_weight_real_ewald;
          boltzmann_weighted_full_squared_summed_real_ewald_gradient.z +=
              analytical_real_ewald_gradient.z * analytical_real_ewald_gradient.z * boltzmann_weight_real_ewald;

          // Hessian
          boltzmann_weighted_hessian_interpolated_summed_real_ewald.ax +=
              interpolated_hessian_real_ewald.ax * boltzmann_weight_real_ewald;
          boltzmann_weighted_hessian_interpolated_summed_real_ewald.ay +=
              interpolated_hessian_real_ewald.ay * boltzmann_weight_real_ewald;
          boltzmann_weighted_hessian_interpolated_summed_real_ewald.az +=
              interpolated_hessian_real_ewald.az * boltzmann_weight_real_ewald;
          boltzmann_weighted_hessian_interpolated_summed_real_ewald.by +=
              interpolated_hessian_real_ewald.by * boltzmann_weight_real_ewald;
          boltzmann_weighted_hessian_interpolated_summed_real_ewald.bz +=
              interpolated_hessian_real_ewald.bz * boltzmann_weight_real_ewald;
          boltzmann_weighted_hessian_interpolated_summed_real_ewald.cz +=
              interpolated_hessian_real_ewald.cz * boltzmann_weight_real_ewald;

          boltzmann_weighted_hessian_full_summed_real_ewald.ax +=
              analytical_real_ewald_hessian.ax * boltzmann_weight_real_ewald;
          boltzmann_weighted_hessian_full_summed_real_ewald.ay +=
              analytical_real_ewald_hessian.ay * boltzmann_weight_real_ewald;
          boltzmann_weighted_hessian_full_summed_real_ewald.az +=
              analytical_real_ewald_hessian.az * boltzmann_weight_real_ewald;
          boltzmann_weighted_hessian_full_summed_real_ewald.by +=
              analytical_real_ewald_hessian.by * boltzmann_weight_real_ewald;
          boltzmann_weighted_hessian_full_summed_real_ewald.bz +=
              analytical_real_ewald_hessian.bz * boltzmann_weight_real_ewald;
          boltzmann_weighted_hessian_full_summed_real_ewald.cz +=
              analytical_real_ewald_hessian.cz * boltzmann_weight_real_ewald;

          boltzmann_weighted_difference_squared_summed_real_ewald_hessian.ax +=
              (analytical_real_ewald_hessian.ax - interpolated_hessian_real_ewald.ax) *
              (analytical_real_ewald_hessian.ax - interpolated_hessian_real_ewald.ax) * boltzmann_weight_real_ewald;
          boltzmann_weighted_difference_squared_summed_real_ewald_hessian.ay +=
              (analytical_real_ewald_hessian.ay - interpolated_hessian_real_ewald.ay) *
              (analytical_real_ewald_hessian.ay - interpolated_hessian_real_ewald.ay) * boltzmann_weight_real_ewald;
          boltzmann_weighted_difference_squared_summed_real_ewald_hessian.az +=
              (analytical_real_ewald_hessian.az - interpolated_hessian_real_ewald.az) *
              (analytical_real_ewald_hessian.az - interpolated_hessian_real_ewald.az) * boltzmann_weight_real_ewald;
          boltzmann_weighted_difference_squared_summed_real_ewald_hessian.by +=
              (analytical_real_ewald_hessian.by - interpolated_hessian_real_ewald.by) *
              (analytical_real_ewald_hessian.by - interpolated_hessian_real_ewald.by) * boltzmann_weight_real_ewald;
          boltzmann_weighted_difference_squared_summed_real_ewald_hessian.bz +=
              (analytical_real_ewald_hessian.bz - interpolated_hessian_real_ewald.bz) *
              (analytical_real_ewald_hessian.bz - interpolated_hessian_real_ewald.bz) * boltzmann_weight_real_ewald;
          boltzmann_weighted_difference_squared_summed_real_ewald_hessian.cz +=
              (analytical_real_ewald_hessian.cz - interpolated_hessian_real_ewald.cz) *
              (analytical_real_ewald_hessian.cz - interpolated_hessian_real_ewald.cz) * boltzmann_weight_real_ewald;

          boltzmann_weighted_full_squared_summed_real_ewald_hessian.ax +=
              analytical_real_ewald_hessian.ax * analytical_real_ewald_hessian.ax * boltzmann_weight_real_ewald;
          boltzmann_weighted_full_squared_summed_real_ewald_hessian.ay +=
              analytical_real_ewald_hessian.ay * analytical_real_ewald_hessian.ay * boltzmann_weight_real_ewald;
          boltzmann_weighted_full_squared_summed_real_ewald_hessian.az +=
              analytical_real_ewald_hessian.az * analytical_real_ewald_hessian.az * boltzmann_weight_real_ewald;
          boltzmann_weighted_full_squared_summed_real_ewald_hessian.by +=
              analytical_real_ewald_hessian.by * analytical_real_ewald_hessian.by * boltzmann_weight_real_ewald;
          boltzmann_weighted_full_squared_summed_real_ewald_hessian.bz +=
              analytical_real_ewald_hessian.bz * analytical_real_ewald_hessian.bz * boltzmann_weight_real_ewald;
          boltzmann_weighted_full_squared_summed_real_ewald_hessian.cz +=
              analytical_real_ewald_hessian.cz * analytical_real_ewald_hessian.cz * boltzmann_weight_real_ewald;
        }
      }

      std::print(stream, "Testing VDW interpolation grid ({}x{}x{}) for {}\n", numberOfVDWGridPoints.x,
                 numberOfVDWGridPoints.y, numberOfVDWGridPoints.z, forceField.pseudoAtoms[index].name);
      std::print(stream, "-------------------------------------------------------------------------------\n");
      std::print(stream, "(Using {} points for testing)\n\n", numberOfGridTestPoints);

      std::print(stream, "Boltzmann average energy VDW (table):      {}\n",
                 boltzmann_weighted_energy_interpolated_summed_vdw / boltzmann_weight_summed_vdw);
      std::print(stream, "Boltzmann average energy VDW (full):       {}\n",
                 boltzmann_weighted_energy_full_summed_vdw / boltzmann_weight_summed_vdw);
      std::print(
          stream, "Boltzmann relative error:                  {}\n\n",
          std::sqrt(boltzmann_weighted_difference_squared_summed_vdw / boltzmann_weighted_full_squared_summed_vdw));

      if (forceField.interpolationScheme == ForceField::InterpolationScheme::Tricubic ||
          forceField.interpolationScheme == ForceField::InterpolationScheme::Triquintic)
      {
        std::print(stream, "Boltzmann average gradient(x) VDW (table): {}\n",
                   boltzmann_weighted_gradient_interpolated_summed_vdw.x / boltzmann_weight_summed_vdw);
        std::print(stream, "Boltzmann average gradient(x) VDW (full):  {}\n",
                   boltzmann_weighted_gradient_full_summed_vdw.x / boltzmann_weight_summed_vdw);
        std::print(stream, "Boltzmann relative error:                  {}\n\n",
                   std::sqrt(boltzmann_weighted_difference_squared_summed_vdw_gradient.x /
                             boltzmann_weighted_full_squared_summed_vdw_gradient.x));

        std::print(stream, "Boltzmann average gradient(y) VDW (table): {}\n",
                   boltzmann_weighted_gradient_interpolated_summed_vdw.y / boltzmann_weight_summed_vdw);
        std::print(stream, "Boltzmann average gradient(y) VDW (full):  {}\n",
                   boltzmann_weighted_gradient_full_summed_vdw.y / boltzmann_weight_summed_vdw);
        std::print(stream, "Boltzmann relative error:                  {}\n\n",
                   std::sqrt(boltzmann_weighted_difference_squared_summed_vdw_gradient.y /
                             boltzmann_weighted_full_squared_summed_vdw_gradient.y));

        std::print(stream, "Boltzmann average gradient(z) VDW (table): {}\n",
                   boltzmann_weighted_gradient_interpolated_summed_vdw.z / boltzmann_weight_summed_vdw);
        std::print(stream, "Boltzmann average gradient(z) VDW (full):  {}\n",
                   boltzmann_weighted_gradient_full_summed_vdw.z / boltzmann_weight_summed_vdw);
        std::print(stream, "Boltzmann relative error:                  {}\n\n",
                   std::sqrt(boltzmann_weighted_difference_squared_summed_vdw_gradient.z /
                             boltzmann_weighted_full_squared_summed_vdw_gradient.z));
      }

      if (forceField.interpolationScheme == ForceField::InterpolationScheme::Triquintic)
      {
        std::print(stream, "Boltzmann average hessian(ax) VDW (table): {}\n",
                   boltzmann_weighted_hessian_interpolated_summed_vdw.ax / boltzmann_weight_summed_vdw);
        std::print(stream, "Boltzmann average hessian(ax) VDW (full):  {}\n",
                   boltzmann_weighted_hessian_full_summed_vdw.ax / boltzmann_weight_summed_vdw);
        std::print(stream, "Boltzmann relative error:                  {}\n\n",
                   std::sqrt(boltzmann_weighted_difference_squared_summed_vdw_hessian.ax /
                             boltzmann_weighted_full_squared_summed_vdw_hessian.ax));

        std::print(stream, "Boltzmann average hessian(ay) VDW (table): {}\n",
                   boltzmann_weighted_hessian_interpolated_summed_vdw.ay / boltzmann_weight_summed_vdw);
        std::print(stream, "Boltzmann average hessian(ay) VDW (full):  {}\n",
                   boltzmann_weighted_hessian_full_summed_vdw.ay / boltzmann_weight_summed_vdw);
        std::print(stream, "Boltzmann relative error:                  {}\n\n",
                   std::sqrt(boltzmann_weighted_difference_squared_summed_vdw_hessian.ay /
                             boltzmann_weighted_full_squared_summed_vdw_hessian.ay));

        std::print(stream, "Boltzmann average hessian(az) VDW (table): {}\n",
                   boltzmann_weighted_hessian_interpolated_summed_vdw.az / boltzmann_weight_summed_vdw);
        std::print(stream, "Boltzmann average hessian(az) VDW (full):  {}\n",
                   boltzmann_weighted_hessian_full_summed_vdw.az / boltzmann_weight_summed_vdw);
        std::print(stream, "Boltzmann relative error:                  {}\n\n",
                   std::sqrt(boltzmann_weighted_difference_squared_summed_vdw_hessian.az /
                             boltzmann_weighted_full_squared_summed_vdw_hessian.az));

        std::print(stream, "Boltzmann average hessian(by) VDW (table): {}\n",
                   boltzmann_weighted_hessian_interpolated_summed_vdw.by / boltzmann_weight_summed_vdw);
        std::print(stream, "Boltzmann average hessian(by) VDW (full):  {}\n",
                   boltzmann_weighted_hessian_full_summed_vdw.by / boltzmann_weight_summed_vdw);
        std::print(stream, "Boltzmann relative error:                  {}\n\n",
                   std::sqrt(boltzmann_weighted_difference_squared_summed_vdw_hessian.by /
                             boltzmann_weighted_full_squared_summed_vdw_hessian.by));

        std::print(stream, "Boltzmann average hessian(bz) VDW (table): {}\n",
                   boltzmann_weighted_hessian_interpolated_summed_vdw.bz / boltzmann_weight_summed_vdw);
        std::print(stream, "Boltzmann average hessian(bz) VDW (full):  {}\n",
                   boltzmann_weighted_hessian_full_summed_vdw.bz / boltzmann_weight_summed_vdw);
        std::print(stream, "Boltzmann relative error:                  {}\n\n",
                   std::sqrt(boltzmann_weighted_difference_squared_summed_vdw_hessian.bz /
                             boltzmann_weighted_full_squared_summed_vdw_hessian.bz));

        std::print(stream, "Boltzmann average hessian(cz) VDW (table): {}\n",
                   boltzmann_weighted_hessian_interpolated_summed_vdw.cz / boltzmann_weight_summed_vdw);
        std::print(stream, "Boltzmann average hessian(cz) VDW (full):  {}\n",
                   boltzmann_weighted_hessian_full_summed_vdw.cz / boltzmann_weight_summed_vdw);
        std::print(stream, "Boltzmann relative error:                  {}\n\n",
                   std::sqrt(boltzmann_weighted_difference_squared_summed_vdw_hessian.cz /
                             boltzmann_weighted_full_squared_summed_vdw_hessian.cz));
      }

      std::print(stream, "Testing Coulomb interpolation grid ({}x{}x{}) for {}\n", numberOfCoulombGridPoints.x,
                 numberOfCoulombGridPoints.y, numberOfCoulombGridPoints.z, forceField.pseudoAtoms[index].name);
      std::print(stream, "-------------------------------------------------------------------------------\n");
      std::print(stream, "(Using {} points for testing)\n\n", numberOfGridTestPoints);

      std::print(stream, "Boltzmann average energy Real Ewald (table):      {}\n",
                 boltzmann_weighted_energy_interpolated_summed_real_ewald / boltzmann_weight_summed_real_ewald);
      std::print(stream, "Boltzmann average energy Real Ewald (full):       {}\n",
                 boltzmann_weighted_energy_full_summed_real_ewald / boltzmann_weight_summed_real_ewald);
      std::print(stream, "Boltzmann relative error:                         {}\n\n",
                 std::sqrt(boltzmann_weighted_difference_squared_summed_real_ewald /
                           boltzmann_weighted_full_squared_summed_real_ewald));

      if (forceField.interpolationScheme == ForceField::InterpolationScheme::Tricubic ||
          forceField.interpolationScheme == ForceField::InterpolationScheme::Triquintic)
      {
        std::print(stream, "Boltzmann average gradient(x) Real Ewald (table): {}\n",
                   boltzmann_weighted_gradient_interpolated_summed_real_ewald.x / boltzmann_weight_summed_real_ewald);
        std::print(stream, "Boltzmann average gradient(x) Real Ewald (full):  {}\n",
                   boltzmann_weighted_gradient_full_summed_real_ewald.x / boltzmann_weight_summed_real_ewald);
        std::print(stream, "Boltzmann relative error:                         {}\n\n",
                   std::sqrt(boltzmann_weighted_difference_squared_summed_real_ewald_gradient.x /
                             boltzmann_weighted_full_squared_summed_real_ewald_gradient.x));

        std::print(stream, "Boltzmann average gradient(y) Real Ewald (table): {}\n",
                   boltzmann_weighted_gradient_interpolated_summed_real_ewald.y / boltzmann_weight_summed_real_ewald);
        std::print(stream, "Boltzmann average gradient(y) Real Ewald (full):  {}\n",
                   boltzmann_weighted_gradient_full_summed_real_ewald.y / boltzmann_weight_summed_real_ewald);
        std::print(stream, "Boltzmann relative error:                         {}\n\n",
                   std::sqrt(boltzmann_weighted_difference_squared_summed_real_ewald_gradient.y /
                             boltzmann_weighted_full_squared_summed_real_ewald_gradient.y));

        std::print(stream, "Boltzmann average gradient(z) Real Ewald (table): {}\n",
                   boltzmann_weighted_gradient_interpolated_summed_real_ewald.z / boltzmann_weight_summed_real_ewald);
        std::print(stream, "Boltzmann average gradient(z) Real Ewald (full):  {}\n",
                   boltzmann_weighted_gradient_full_summed_real_ewald.z / boltzmann_weight_summed_real_ewald);
        std::print(stream, "Boltzmann relative error:                         {}\n\n",
                   std::sqrt(boltzmann_weighted_difference_squared_summed_real_ewald_gradient.z /
                             boltzmann_weighted_full_squared_summed_real_ewald_gradient.z));
      }

      if (forceField.interpolationScheme == ForceField::InterpolationScheme::Triquintic)
      {
        std::print(stream, "Boltzmann average Hessian(ax) Real Ewald (table): {}\n",
                   boltzmann_weighted_hessian_interpolated_summed_real_ewald.ax / boltzmann_weight_summed_real_ewald);
        std::print(stream, "Boltzmann average Hessian(ax) Real Ewald (full):  {}\n",
                   boltzmann_weighted_hessian_full_summed_real_ewald.ax / boltzmann_weight_summed_real_ewald);
        std::print(stream, "Boltzmann relative error:                         {}\n\n",
                   std::sqrt(boltzmann_weighted_difference_squared_summed_real_ewald_hessian.ax /
                             boltzmann_weighted_full_squared_summed_real_ewald_hessian.ax));

        std::print(stream, "Boltzmann average Hessian(ay) Real Ewald (table): {}\n",
                   boltzmann_weighted_hessian_interpolated_summed_real_ewald.ay / boltzmann_weight_summed_real_ewald);
        std::print(stream, "Boltzmann average Hessian(ay) Real Ewald (full):  {}\n",
                   boltzmann_weighted_hessian_full_summed_real_ewald.ay / boltzmann_weight_summed_real_ewald);
        std::print(stream, "Boltzmann relative error:                         {}\n\n",
                   std::sqrt(boltzmann_weighted_difference_squared_summed_real_ewald_hessian.ay /
                             boltzmann_weighted_full_squared_summed_real_ewald_hessian.ay));

        std::print(stream, "Boltzmann average Hessian(az) Real Ewald (table): {}\n",
                   boltzmann_weighted_hessian_interpolated_summed_real_ewald.az / boltzmann_weight_summed_real_ewald);
        std::print(stream, "Boltzmann average Hessian(az) Real Ewald (full):  {}\n",
                   boltzmann_weighted_hessian_full_summed_real_ewald.az / boltzmann_weight_summed_real_ewald);
        std::print(stream, "Boltzmann relative error:                         {}\n\n",
                   std::sqrt(boltzmann_weighted_difference_squared_summed_real_ewald_hessian.az /
                             boltzmann_weighted_full_squared_summed_real_ewald_hessian.az));

        std::print(stream, "Boltzmann average Hessian(by) Real Ewald (table): {}\n",
                   boltzmann_weighted_hessian_interpolated_summed_real_ewald.by / boltzmann_weight_summed_real_ewald);
        std::print(stream, "Boltzmann average Hessian(by) Real Ewald (full):  {}\n",
                   boltzmann_weighted_hessian_full_summed_real_ewald.by / boltzmann_weight_summed_real_ewald);
        std::print(stream, "Boltzmann relative error:                         {}\n\n",
                   std::sqrt(boltzmann_weighted_difference_squared_summed_real_ewald_hessian.by /
                             boltzmann_weighted_full_squared_summed_real_ewald_hessian.by));

        std::print(stream, "Boltzmann average Hessian(bz) Real Ewald (table): {}\n",
                   boltzmann_weighted_hessian_interpolated_summed_real_ewald.bz / boltzmann_weight_summed_real_ewald);
        std::print(stream, "Boltzmann average Hessian(bz) Real Ewald (full):  {}\n",
                   boltzmann_weighted_hessian_full_summed_real_ewald.bz / boltzmann_weight_summed_real_ewald);
        std::print(stream, "Boltzmann relative error:                         {}\n\n",
                   std::sqrt(boltzmann_weighted_difference_squared_summed_real_ewald_hessian.bz /
                             boltzmann_weighted_full_squared_summed_real_ewald_hessian.bz));

        std::print(stream, "Boltzmann average Hessian(cz) Real Ewald (table): {}\n",
                   boltzmann_weighted_hessian_interpolated_summed_real_ewald.cz / boltzmann_weight_summed_real_ewald);
        std::print(stream, "Boltzmann average Hessian(cz) Real Ewald (full):  {}\n",
                   boltzmann_weighted_hessian_full_summed_real_ewald.cz / boltzmann_weight_summed_real_ewald);
        std::print(stream, "Boltzmann relative error:                         {}\n\n",
                   std::sqrt(boltzmann_weighted_difference_squared_summed_real_ewald_hessian.cz /
                             boltzmann_weighted_full_squared_summed_real_ewald_hessian.cz));
      }

      std::print(stream, "\n");
      std::flush(stream);
    }
  }
}
