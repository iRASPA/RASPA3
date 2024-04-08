module;

#ifdef USE_LEGACY_HEADERS
#include <vector>
#include <span>
#include <tuple>
#include <string>
#include <fstream>
#if defined(__has_include) && __has_include(<mdspan>)
  #include <mdspan>
#endif
#endif

export module breakthrough;

#ifndef USE_LEGACY_HEADERS
import <vector>;
import <span>;
import <tuple>;
import <string>;
import <fstream>;
#if defined(__has_include) && __has_include(<mdspan>)
  import <mdspan>;
#endif
#endif

#if !(defined(__has_include) && __has_include(<mdspan>))
  import mdspan;
#endif

import input_reader;
import component;
import system;
import mixture_prediction;

export struct Breakthrough
{
  public:
    Breakthrough(System &system);

    void print() const;

    std::string writeHeader();
    void run(std::ostream &stream);

    void createPlotScript();
    void createMovieScripts();

  private:
    const System &system;
    const std::string displayName;
    const std::vector<Component> components;
    //size_t carrierGasComponent{ 0 }; 
    size_t Ncomp;      // number of components
    size_t Ngrid;      // number of grid points

    size_t printEvery; // print time step to the screen every printEvery steps
    size_t writeEvery; // write data to files every writeEvery steps

    double T;          // absolute temperature [K]
    double p_total;    // total pressure column [Pa]
    double dptdx;      // pressure gradient [N/m3]
    double epsilon;    // void-fraction of the column [-]
    double rho_p;      // particle density [kg/m3]
    double v_in;       // interstitial velocity at the begin of the column [m/s]
    
    double L;          // length of the column
    double dx;         // spacing in spatial direction
    double dt;         // timestep integration
    size_t Nsteps;     // total number of steps
    bool autoSteps;    // use automatic number of steps
    bool pulse;        // pulsed inlet condition for breakthrough
    double tpulse;     // pulse time
    MixturePrediction mixture;
    std::pair<size_t, size_t> iastPerformance{ 0, 0 };
    
    // vector of size 'Ncomp'
    std::vector<double> prefactor;
    std::vector<double> Yi;        // ideal gas mol-fraction for each component
    std::vector<double> Xi;        // adsorbed mol-fraction for each component
    std::vector<double> Ni;        // number of molecules for each component

    // vector of size '(Ngrid + 1)'
    std::vector<double> V;         // interstitial gas velocity along the column
    std::vector<double> Vnew; 
    std::vector<double> Pt;        // total pressure along the column

    // vector of size '(Ngrid + 1) * Ncomp', for each grid point, data per component (contiguous)
    std::vector<double> P_vector;         // partial pressure at every grid point for each component
    #if defined(__has_include) && __has_include(<mdspan>)
      std::mdspan<double, std::dextents<size_t, 2>> P;
    #else
      std::experimental::mdspan<double, std::experimental::dextents<size_t, 2>> P;
    #endif
    std::vector<double> Pnew_vector;
    #if defined(__has_include) && __has_include(<mdspan>)
      std::mdspan<double, std::dextents<size_t, 2>> Pnew;
    #else
      std::experimental::mdspan<double, std::experimental::dextents<size_t, 2>> Pnew;
    #endif
    std::vector<double> Q_vector;         // volume-averaged adsorption amount at every grid point for each component
    #if defined(__has_include) && __has_include(<mdspan>)
      std::mdspan<double, std::dextents<size_t, 2>> Q;
    #else
      std::experimental::mdspan<double, std::experimental::dextents<size_t, 2>> Q;
    #endif
    std::vector<double> Qnew_vector;
    #if defined(__has_include) && __has_include(<mdspan>)
      std::mdspan<double, std::dextents<size_t, 2>> Qnew;
    #else
      std::experimental::mdspan<double, std::experimental::dextents<size_t, 2>> Qnew;
    #endif
    std::vector<double> Qeq_vector;       // equilibrium adsorption amount at every grid point for each component
    #if defined(__has_include) && __has_include(<mdspan>)
      std::mdspan<double, std::dextents<size_t, 2>> Qeq;
    #else
      std::experimental::mdspan<double, std::experimental::dextents<size_t, 2>> Qeq;
    #endif
    std::vector<double> Qeqnew_vector;
    #if defined(__has_include) && __has_include(<mdspan>)
      std::mdspan<double, std::dextents<size_t, 2>> Qeqnew;
    #else
      std::experimental::mdspan<double, std::experimental::dextents<size_t, 2>> Qeqnew;
    #endif
    std::vector<double> Dpdt_vector;      // derivative of P with respect to time
    #if defined(__has_include) && __has_include(<mdspan>)
      std::mdspan<double, std::dextents<size_t, 2>> Dpdt;
    #else
      std::experimental::mdspan<double, std::experimental::dextents<size_t, 2>> Dpdt;
    #endif
    std::vector<double> Dpdtnew_vector;
    #if defined(__has_include) && __has_include(<mdspan>)
      std::mdspan<double, std::dextents<size_t, 2>> Dpdtnew;
    #else
      std::experimental::mdspan<double, std::experimental::dextents<size_t, 2>> Dpdtnew;
    #endif
    std::vector<double> Dqdt_vector;      // derivative of Q with respect to time
    #if defined(__has_include) && __has_include(<mdspan>)
      std::mdspan<double, std::dextents<size_t, 2>> Dqdt;
    #else
      std::experimental::mdspan<double, std::experimental::dextents<size_t, 2>> Dqdt;
    #endif
    std::vector<double> Dqdtnew_vector;
    #if defined(__has_include) && __has_include(<mdspan>)
      std::mdspan<double, std::dextents<size_t, 2>> Dqdtnew;
    #else
      std::experimental::mdspan<double, std::experimental::dextents<size_t, 2>> Dqdtnew;
    #endif
    std::vector<double> cachedP0;  // cached hypothetical pressure
    std::vector<double> cachedPsi; // cached reduced grand potential over the column


    enum class IntegrationScheme
    {
      SSP_RK = 0,
      Iterative = 1
    };

    #if defined(__has_include) && __has_include(<mdspan>)
      void computeFirstDerivatives(std::mdspan<double, std::dextents<size_t, 2>> &dqdt,
                                   std::mdspan<double, std::dextents<size_t, 2>> &dpdt,
                                   const std::mdspan<double, std::dextents<size_t, 2>> &q_eq,
                                   const std::mdspan<double, std::dextents<size_t, 2>> &q,
                                   const std::vector<double> &v,
                                   const std::mdspan<double, std::dextents<size_t, 2>> &p);
    #else
      void computeFirstDerivatives(std::experimental::mdspan<double, std::experimental::dextents<size_t, 2>> &dqdt,
                                   std::experimental::mdspan<double, std::experimental::dextents<size_t, 2>> &dpdt,
                                   const std::experimental::mdspan<double, std::experimental::dextents<size_t, 2>> &q_eq,
                                   const std::experimental::mdspan<double, std::experimental::dextents<size_t, 2>> &q,
                                   const std::vector<double> &v,
                                   const std::experimental::mdspan<double, std::experimental::dextents<size_t, 2>> &p);
    #endif

    void computeEquilibriumLoadings();

    void computeVelocity();

    void createMovieScriptColumnV();
    void createMovieScriptColumnPt();
    void createMovieScriptColumnQ();
    void createMovieScriptColumnQeq();
    void createMovieScriptColumnP();
    void createMovieScriptColumnDpdt();
    void createMovieScriptColumnDqdt();
    void createMovieScriptColumnPnormalized();
};
