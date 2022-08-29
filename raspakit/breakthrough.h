#include <cstddef>
#include <vector>

#include "component.h"
#include "inputreader.h"
#include "iast.h"

const double R=8.31446261815324;

struct Breakthrough
{
  Breakthrough(const InputReader &inputreader);

  void print();
  void initialize();
  void run();

  std::vector<Component> components;
  size_t Ncomp;      // number of components
  size_t Ngrid;      // number of grid points

  size_t printEvery; // print time step to the screen every printEvery steps

  double T;          // absolute temperature [K]
  double p_total;    // total pressure column [Pa]
  double dptdx;      // pressure gradient [N/m3]
  double epsilon;    // void-fraction of the column [-]
  double rho_p;      // particle density [kg/m3]
  double u_in;       // initial velocity at the begin of the column [m/s]
  
  double dx;         // spacing in spatial direction
  double L;          // length of the column
  double dt;         // timestep integration
  size_t Nsteps;     // total number of steps
  size_t Niter;      // number of inner iterations
  double Pi;
  std::vector<double> prefactor;

  std::vector<double> Yi;     // ideal gas mol-fraction for each component
  std::vector<double> Xi;     // adsorbed mol-fraction for each component
  std::vector<double> Pi0;
  std::vector<double> Ni;     // number of molecules for each component
  std::vector<double> U;      // superficial gas velocity along the column
  std::vector<double> Unew; 
  std::vector<double> Pt;     // total pressure along the column
  std::vector<double> Ptinit; // initial pressure along the column

  // vector of (Ngrid + 1) * Ncomp, component data contiguous
  std::vector<double> P;        // partial pressure at every grid point for each component
  std::vector<double> Pnew;
  std::vector<double> Q;        // volume-averaged adsorption amount at every grid point for each component
  std::vector<double> Qnew;
  std::vector<double> Qeq;      // equilibrium adsorption amount at every grid point for each component
  std::vector<double> Qeqnew;
  std::vector<double> Dpdt;     // derivative of P with respect to time
  std::vector<double> Dpdtnew;
  std::vector<double> Dqdt;     // derivative of Q with respect to time
  std::vector<double> Dqdtnew;
  std::vector<double> Polditer;
  std::vector<double> Qolditer;

  enum class TimeIntegrationScheme
  {
    CranckNicolson = 0,
    BackwardEuler = 1
  };
  enum class AdvectiveIntegrationScheme
  {
    Upwind = 0,
    CentralDifferencing = 1
  };
  enum class VelocityIntegrationScheme
  {
    Upwind = 0,
    CranckNicolson = 1,
    CentralDifferencing = 2
  };

  TimeIntegrationScheme timeIntegrationScheme;
  AdvectiveIntegrationScheme advectiveIntegrationScheme;
  VelocityIntegrationScheme velocityIntegrationScheme;
  IastScheme iastScheme;

  size_t numberOfIASTSteps{ 0 };

  void derivatives(std::vector<double> &Dqdt,
                   std::vector<double> &Dpdt,
                   const std::vector<double> &Qeq,
                   const std::vector<double> &Q,
                   const std::vector<double> &U,
                   const std::vector<double> &P);

  void velocity(std::vector<double> &U,
                const std::vector<double> &Q,
                const std::vector<double> &Qeq,
                const std::vector<double> &P,
                const std::vector<double> &Pt);
};
