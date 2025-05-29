module;
     
#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>  
#include <chrono> 
#include <complex>
#include <cstddef>
#include <cmath>
#include <algorithm>
#include <exception>
#include <format>
#include <fstream>
#include <istream>
#include <map>                                      
#include <ostream>                                  
#include <print>                                    
#include <source_location>
#include <sstream>
#include <tuple>          
#include <utility>
#include <vector>
#include <string>
#if defined(__has_include) && __has_include(<mdspan>)
#include <mdspan>
#endif
#endif
  
module energy_void_fraction;

import opencl;
import float4;
import double4;
import forcefield;
import framework;
import units;

EnergyVoidFraction::EnergyVoidFraction()
{
};


EnergyVoidFraction::~EnergyVoidFraction()
{
}


void EnergyVoidFraction::run(const ForceField &forceField, const Framework &framework)
{
  int3 grid_size = int3(128,128,128);
  double2 probeParameter = double2(10.9 * Units::KelvinToEnergy, 2.64);
  double cutoff = forceField.cutOffFrameworkVDW;
  double3x3 unitCell = framework.simulationBox.cell;
  int3 numberOfReplicas = framework.simulationBox.smallestNumberOfUnitCellsForMinimumImagesConvention(cutoff);
  std::vector<double3> positions = framework.fractionalAtomPositionsUnitCell();
  std::vector<double2> potentialParameters = framework.atomUnitCellLennardJonesPotentialParameters(forceField);
  std::chrono::system_clock::time_point time_begin, time_end;

  time_begin = std::chrono::system_clock::now();

  size_t numberOfAtoms = positions.size();
  int temp = grid_size.x*grid_size.y*grid_size.z;

  std::vector<double> outputData = std::vector<double>(static_cast<size_t>(temp));

  double3 correction = double3(1.0/double(numberOfReplicas.x), 1.0/double(numberOfReplicas.y), 1.0/double(numberOfReplicas.z));

  double3x3 replicaCell = double3x3(double(numberOfReplicas.x) * unitCell[0],
                                    double(numberOfReplicas.y) * unitCell[1],
                                    double(numberOfReplicas.z) * unitCell[2]);

  int totalNumberOfReplicas = numberOfReplicas.x * numberOfReplicas.y * numberOfReplicas.z;
  std::vector<double3> replicaVector(static_cast<size_t>(totalNumberOfReplicas));
  size_t index = 0;
  for(size_t i=0; i<static_cast<size_t>(numberOfReplicas.x); i++)
  {
    for(size_t j=0; j<static_cast<size_t>(numberOfReplicas.y); j++)
    {
      for(size_t k=0; k<static_cast<size_t>(numberOfReplicas.z); k++)
      {
        replicaVector[index] = double3((double(i)/double(numberOfReplicas.x)),
                                       (double(j)/double(numberOfReplicas.y)),
                                       (double(k)/double(numberOfReplicas.z)));
        index += 1;
      }
    }
  }

  for(size_t z=0; z<static_cast<size_t>(grid_size.z);z++)
  {
    for(size_t y=0; y<static_cast<size_t>(grid_size.y); y++)
    {
      for(size_t x=0 ; x<static_cast<size_t>(grid_size.x); x++)
      {
        double3 gridPosition = correction * double3(double(x)/double(grid_size.x-1),double(y)/double(grid_size.y-1),double(z)/double(grid_size.z-1));

        double value = 0.0;
        for(size_t i=0 ; i<numberOfAtoms; i++)
        {
          double3 position = correction * positions[i];
          double2 currentPotentialParameters = potentialParameters[i];

          // use 4 x epsilon for a probe epsilon of unity
          double epsilon = 4.0*sqrt(currentPotentialParameters.x * probeParameter.x);

          // mixing rule for the atom and the probe
          double sigma = 0.5 * (currentPotentialParameters.y + probeParameter.y);

          for(size_t j=0;j<static_cast<size_t>(totalNumberOfReplicas);j++)
          {
            double3 ds = gridPosition - position - replicaVector[j];
            ds.x -= std::rint(ds.x);
            ds.y -= std::rint(ds.y);
            ds.z -= std::rint(ds.z);
            double3 dr = replicaCell * ds;

            double rr = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;
            if (rr<12.0*12.0)
            {
              double sigma2rr = sigma*sigma/rr;
              double rri3 = sigma2rr * sigma2rr * sigma2rr;
              value += epsilon*(rri3*(rri3-1.0));
            }
          }
        }

        outputData[x + y* static_cast<size_t>(grid_size.x) + z * static_cast<size_t>(grid_size.x * grid_size.y)] += std::min(value,10000000.0);
      }
    }
  }

  double sumBoltzmannWeight = 0.0;
  for(const double &value: outputData)
  {
    sumBoltzmannWeight += exp(-(1.0/(Units::KB*298.0) * value)); 
  }

  time_end = std::chrono::system_clock::now();

  std::chrono::duration<double> timing = time_end - time_begin;

  std::ofstream myfile;
  myfile.open(framework.name + ".energy.vf.cpu.txt");
  std::print(myfile, "# Void-fraction using energy-based method\n");
  std::print(myfile, "# CPU Timing: {} [s]\n", timing.count());
  myfile << sumBoltzmannWeight/static_cast<double>(grid_size.x * grid_size.y * grid_size.z)  << std::endl;
  myfile.close();
}
