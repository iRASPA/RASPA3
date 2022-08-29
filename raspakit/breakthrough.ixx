export module breakthrough;

import <vector>;
import <span>;
import <tuple>;
import <string>;
import <fstream>;

import input_reader;
import component;
import system;
import mixture_prediction;

export struct Breakthrough
{
  public:
    Breakthrough(System &system);

    void writeHeader(std::ostream &stream);
    void run(std::ostream &stream);

    void createPlotScript(std::string directoryName);
    void createMovieScripts(std::string directoryName);

  private:
    const System &system;
    std::string displayName;
    std::span<const Component> components;
    size_t carrierGasComponent{ 0 }; 
    size_t Ncomp;      // number of components
    size_t Ngrid;      // number of grid points

    size_t printEvery; // print time step to the screen every printEvery steps
    size_t writeEvery; // write data to files every writeEvery steps

    double T;
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
    std::vector<double> P;         // partial pressure at every grid point for each component
    std::vector<double> Pnew;
    std::vector<double> Q;         // volume-averaged adsorption amount at every grid point for each component
    std::vector<double> Qnew;
    std::vector<double> Qeq;       // equilibrium adsorption amount at every grid point for each component
    std::vector<double> Qeqnew;
    std::vector<double> Dpdt;      // derivative of P with respect to time
    std::vector<double> Dpdtnew;
    std::vector<double> Dqdt;      // derivative of Q with respect to time
    std::vector<double> Dqdtnew;
    std::vector<double> cachedP0;  // cached hypothetical pressure
    std::vector<double> cachedPsi; // cached reduced grand potential over the column

    std::pair<size_t, size_t> iastPerformance{ 0, 0 };

    void computeFirstDerivatives(std::vector<double> &dqdt,
                                 std::vector<double> &dpdt,
                                 const std::vector<double> &q_eq,
                                 const std::vector<double> &q,
                                 const std::vector<double> &v,
                                 const std::vector<double> &p);

    void computeEquilibriumLoadings(MixturePrediction &mixture);

    void computeVelocity();

    void createMovieScriptColumnV(std::string directoryName);
    void createMovieScriptColumnPt(std::string directoryName);
    void createMovieScriptColumnQ(std::string directoryName);
    void createMovieScriptColumnQeq(std::string directoryName);
    void createMovieScriptColumnP(std::string directoryName);
    void createMovieScriptColumnDpdt(std::string directoryName);
    void createMovieScriptColumnDqdt(std::string directoryName);
    void createMovieScriptColumnPnormalized(std::string directoryName);
};
