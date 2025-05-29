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
  
module energy_surface_area;

import opencl;
import float4;
import double4;
import marching_cubes;
import forcefield;
import framework;
import units;

EnergySurfaceArea::EnergySurfaceArea()
{
};


void EnergySurfaceArea::run(const ForceField &forceField, const Framework &framework)
{
  double isoValue = 0.0;
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
  for(int i=0; i<numberOfReplicas.x; i++)
  {
    for(int j=0; j<numberOfReplicas.y; j++)
    {
      for(int k=0; k<numberOfReplicas.z; k++)
      {
        replicaVector[index] = double3((double(i)/double(numberOfReplicas.x)),
                                       (double(j)/double(numberOfReplicas.y)),
                                       (double(k)/double(numberOfReplicas.z)));
        index += 1;
      }
    }
  }

  for(size_t z = 0; z < static_cast<size_t>(grid_size.z); z++)
  {
    for(size_t y = 0; y < static_cast<size_t>(grid_size.y); y++)
    {
      for(size_t x = 0 ; x < static_cast<size_t>(grid_size.x); x++)
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
              value += epsilon * (rri3 * (rri3 - 1.0));
            }
          }
        }

        outputData[x + y* static_cast<size_t>(grid_size.x) + z * static_cast<size_t>(grid_size.x * grid_size.y)] += std::min(value,10000000.0);
      }
    }
  }


  int largestSize = std::max({grid_size.x,grid_size.y,grid_size.z});
  int powerOfTwo = 1;
  while(largestSize > std::pow(2,powerOfTwo))
  {
    powerOfTwo += 1;
  }
  int size = static_cast<int>(std::pow(2,powerOfTwo)); // the encompassing size to use as textures (size 16,32,64,128,256, and 512 are supported).

  // brute force implementation
  // Create marching cubes object
  MarchingCubes cube(grid_size.x,grid_size.y,grid_size.z);

  // Initiate the cube
  cube.init_all();

  // Set the data
  for(size_t i=0; i < static_cast<size_t>(grid_size.x); i++)
  {
    for(size_t j=0; j < static_cast<size_t>(grid_size.y); j++)
    {
      for(size_t k=0; k < static_cast<size_t>(grid_size.z); k++)
      {
        double value = outputData[i + static_cast<size_t>(size) * j + k * static_cast<size_t>(size * size)];
        cube.set_data(value, i, j, k);
      }
    }
  }

  cube.run( isoValue );

  size_t numberOfTriangles = cube.ntrigs();
  std::vector<float4> triangleData{};
  triangleData.reserve(3*3*numberOfTriangles);

    // Fetch the info
  for(size_t i=0; i < cube.ntrigs(); i++)
  {
    const Triangle* tri = cube.trig(static_cast<int>(i));

    const Vertex* vertex1 = cube.vert(tri->v1);
    triangleData.push_back(double4(vertex1->x/grid_size.x,vertex1->y/grid_size.y,vertex1->z/grid_size.z,1.0));
    triangleData.push_back(double4(vertex1->nx,vertex1->ny,vertex1->nz,0.0));
    triangleData.push_back(double4(0.0,0.0,0.0,0.0));

    const Vertex* vertex2 = cube.vert(tri->v2);
    triangleData.push_back(double4(vertex2->x/grid_size.x,vertex2->y/grid_size.y,vertex2->z/grid_size.z,1.0));
    triangleData.push_back(double4(vertex2->nx,vertex2->ny,vertex2->nz,0.0));
    triangleData.push_back(double4(0.0,0.0,0.0,0.0));

    const Vertex* vertex3 = cube.vert(tri->v3);
    triangleData.push_back(double4(vertex3->x/grid_size.x,vertex3->y/grid_size.y,vertex3->z/grid_size.z,1.0));
    triangleData.push_back(double4(vertex3->nx,vertex3->ny,vertex3->nz,0.0));
    triangleData.push_back(double4(0.0,0.0,0.0,0.0));
  }

  double accumulated_surface_area=0.0;
  for(size_t i=0; i<triangleData.size(); i+=9)
  {
    double3 p1 = unitCell * double3(static_cast<double>(triangleData[i].x), static_cast<double>(triangleData[i].y), static_cast<double>(triangleData[i].z));
    double3 p2 = unitCell * double3(static_cast<double>(triangleData[i+3].x), static_cast<double>(triangleData[i+3].y), static_cast<double>(triangleData[i+3].z));
    double3 p3 = unitCell * double3(static_cast<double>(triangleData[i+6].x) , static_cast<double>(triangleData[i+6].y), static_cast<double>(triangleData[i+6].z));
    double3 v = double3::cross(p2-p1, p3-p1);
    double area = 0.5 * v.length();
    if(std::isfinite(area) && std::fabs(area) < 1.0 )
    {
      accumulated_surface_area += area;
    }
  }

  time_end = std::chrono::system_clock::now();

  std::chrono::duration<double> timing = time_end - time_begin;

  std::ofstream myfile;
  myfile.open(framework.name + ".energy.sa.cpu.txt");
  std::print(myfile, "# Surface area using energy-based method\n");
  std::print(myfile, "# CPU Timing: {} [s]\n", timing.count());
  myfile << accumulated_surface_area << " [A^2]" << std::endl;
  myfile << accumulated_surface_area * Units::Angstrom * Units::Angstrom * Units::AvogadroConstant / framework.unitCellMass << " [m^2/g]" << std::endl;
  myfile << 1.0e4 * accumulated_surface_area / framework.simulationBox.volume << " [m^2/cm^3]" << std::endl;
  myfile.close();
}
