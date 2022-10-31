module;

module breakthrough;

import <vector>;
import <span>;
import <cmath>;
import <string>;
import <iostream>;
import <fstream>;
import <limits>;
import <filesystem>;
import <algorithm>;
import <numeric>;
import <sstream>;
import <chrono>;

import print;
import input_reader;
import component;
import system;
import simulationbox;
import mixture_prediction;

const double R=8.31446261815324;

inline double maxVectorDifference(const std::vector<double> &v, const std::vector<double> &w)
{
  if(v.empty() || w.empty()) return 0.0;
  if(v.size() != w.size()) throw std::runtime_error("Error: unequal vector size\n");

  double max = std::abs(v[0] - w[0]);
  for(size_t i = 1; i < v.size(); ++i)
  {
    double temp = std::abs(v[i] - w[i]);
    if(temp > max) max = temp;
  }
  return max;
}


// allow std::pairs to be added
template <typename T,typename U>                                                   
std::pair<T,U> operator+(const std::pair<T,U> & l,const std::pair<T,U> & r) {   
    return {l.first+r.first,l.second+r.second};
}
template <typename T, typename U>
std::pair<T,U> &operator+=(std::pair<T,U> & l, const std::pair<T,U> & r) {   
    l.first += r.first;
    l.second += r.second;
    return l;
}

Breakthrough::Breakthrough(System &system):
    system(system),
    displayName(system.components.front().name),
    components(system.spanOfAdsorbateComponents()),
    Ncomp(components.size()),
    Ngrid(system.columnNumberOfGridPoints),
    printEvery(10000),
    writeEvery(5000),
    T(system.simulationBox.temperature),
    p_total(system.columnTotalPressure),
    dptdx(system.columnPressureGradient),
    epsilon(system.columnVoidFraction),
    rho_p(system.columnParticleDensity),
    v_in(system.columnEntranceVelocity),
    L(system.columnLength),
    dx(L / static_cast<double>(Ngrid)),
    dt(system.columnTimeStep),
    Nsteps(system.columnNumberOfTimeSteps),
    autoSteps(system.columnAutoNumberOfTimeSteps),
    prefactor(Ncomp),
    Yi(Ncomp),
    Xi(Ncomp),
    Ni(Ncomp),
    V(Ngrid+1),
    Vnew(Ngrid+1),
    Pt(Ngrid+1),
    P((Ngrid + 1) * Ncomp),
    Pnew((Ngrid + 1) * Ncomp),
    Q((Ngrid + 1) * Ncomp),
    Qnew((Ngrid + 1) * Ncomp),
    Qeq((Ngrid + 1) * Ncomp),
    Qeqnew((Ngrid + 1) * Ncomp),
    Dpdt((Ngrid + 1) * Ncomp),
    Dpdtnew((Ngrid + 1) * Ncomp),
    Dqdt((Ngrid + 1) * Ncomp),
    Dqdtnew((Ngrid + 1) * Ncomp),
    cachedP0((Ngrid + 1) * Ncomp),
    cachedPsi(Ngrid + 1)
{
  auto itr = std::find_if(components.begin(), components.end(),
       [&](const Component &c) {
          return c.isCarrierGas;
       });

  if (itr == components.end())
  {
    throw std::runtime_error("Error [Breakthrough]: no carrier gas component present");
  }
  carrierGasComponent = static_cast<size_t>(std::distance(components.begin(), itr));
}

void Breakthrough::writeHeader(std::ostream &stream)
{
  std::print(stream, "\nColumn properties\n");
  std::print(stream, "=======================================================\n");
  std::print(stream, "Display-name:                          {}\n", displayName);
  std::print(stream, "Temperature:                           {} [K]\n", T);
  std::print(stream, "Column length:                         {} [m]\n", L);
  std::print(stream, "Column void-fraction:                  {} [-]\n", epsilon);
  std::print(stream, "Particle density:                      {} [kg/m^3]\n", rho_p);
  std::print(stream, "Total pressure:                        {} [Pa]\n", p_total);
  std::print(stream, "Pressure gradient:                     {} [Pa/m]\n", dptdx);
  std::print(stream, "Column entrance interstitial velocity: {} [m/s]\n", v_in);
  std::print(stream, "\n\n");

  std::print("Breakthrough settings\n");
  std::print("=======================================================\n");
  std::print("Number of time steps:          {}\n", Nsteps);
  std::print("Print every step:              {}\n", printEvery);
  std::print("Write data every step:         {}\n", writeEvery);
  std::print("\n\n");

  std::print(stream, "Integration details\n");
  std::print(stream, "=======================================================\n");
  std::print(stream, "Time step:                     {} [s]\n", dt);
  std::print(stream, "Number of column grid points:  {} [-]\n", Ngrid);
  std::print(stream, "Column spacing:                {} [m]\n", dx);
  std::print(stream, "\n\n");

  std::print("Component data\n");
  std::print("=======================================================\n");
  for(size_t i = 0; i < Ncomp; ++i)
  {
    components[i].printBreakthroughStatus(stream);
    std::print(stream, "\n");
  }
}

void Breakthrough::run(std::ostream &stream)
{
  MixturePrediction mixture(system, 0.0, 0.0, 100, MixturePrediction::PressureScale::Log);

  std::string directoryNameString = std::print("Breakthrough/System_{}/", system.systemId);;
  std::filesystem::path directoryName{ directoryNameString };
  std::filesystem::create_directories(directoryName);

  // create the output files
  std::vector<std::ofstream> streams;
  for (const Component &component: components)
  {
    std::string fileName = directoryNameString + "component_" + std::to_string(component.componentId) + "_" + component.name + ".data";
    streams.emplace_back(std::ofstream{ fileName });
  }

  std::ofstream movieStream(directoryNameString + "column.data");

  createPlotScript(directoryNameString);
  createMovieScripts(directoryNameString);

  // precomputed factor for mass transfer
  for(size_t j = 0; j < components.size(); ++j)
  {
    prefactor[j] = R * T * ((1.0 - epsilon) / epsilon) * rho_p * components[j].massTransferCoefficient;
  }

  // set P and Q to zero
  std::fill(P.begin(), P.end(), 0.0);
  std::fill(Q.begin(), Q.end(), 0.0);

  // initial pressure along the column
  std::vector<double> pt_init(Ngrid + 1);

  // set the initial total pressure along the column assuming the pressure gradient is constant
  for(size_t i = 0; i < Ngrid + 1; ++i)
  {
    pt_init[i] = p_total + dptdx * static_cast<double>(i) * dx;
  }

  // initialize the interstitial gas velocity in the column
  for(size_t i = 0; i < Ngrid + 1; ++i)
  {
    V[i] = v_in * p_total / pt_init[i];
  }

  // set the partial pressure of the carrier gas to the total initial pressure
  // for the column except for the entrance (i=0)
  for(size_t i = 1; i < Ngrid + 1; ++i)
  {
    P[i * Ncomp + carrierGasComponent] = pt_init[i];
  }

  // at the column entrance, the mol-fractions of the components in the gas phase are fixed
  // the partial pressures of the components at the entrance are the mol-fractions times the 
  // total pressure
  for(size_t j = 0; j < Ncomp; ++j)
  {
    P[0 * Ncomp + j] = p_total * components[j].molFraction;
  }

  // at the entrance: mol-fractions Yi are the gas-phase mol-fractions
  // for the column: the initial mol-fraction of the carrier-gas is 1, and 0 for the other components
  //
  // the K of the carrier gas is chosen as zero 
  // so Qeq is zero for all components in the column after the entrance
  // only the values for Yi at the entrance are effected by adsorption
  // Note: for non-zero K of the carrier gas the IAST does not converge
  for(size_t i = 0; i < Ngrid + 1; ++i)
  {
    double sum = 0.0;
    for(size_t j = 0; j < Ncomp; ++j)
    {
      Yi[j] = std::max(P[i * Ncomp + j] / pt_init[i], 0.0);
      sum += Yi[j];
    }
    for(size_t j = 0; j < Ncomp; ++j)
    {
      Yi[j] /= sum;
    }

    iastPerformance += mixture.predictMixture(Yi, pt_init[i], Xi, Ni, &cachedP0[i * Ncomp], cachedPsi[i]);

    for(size_t j = 0; j < Ncomp; ++j)
    {
      Qeq[i * Ncomp + j] = Ni[j];
    }
  }

  for(size_t i = 0; i < Ngrid + 1; ++i)
  {
    Pt[i] = 0.0;
    for(size_t j = 0; j < Ncomp; ++j)
    {
      Pt[i] += P[i * Ncomp + j];
    }
  }


  for(size_t step = 0; (step < Nsteps || autoSteps); ++step)
  {
    double t = static_cast<double>(step) * dt;

    if(step % writeEvery == 0)
    {
      // write breakthrough output to files
      // column 1: dimensionless time
      // column 2: time [minutes]
      // column 3: normalized partial pressure
      for(size_t j = 0; j < Ncomp; ++j)
      {
        std::print(streams[j], "{} {} {}\n", t * v_in / L,  t/60.0, 
            P[Ngrid * Ncomp + j] / (p_total * components[j].molFraction));
      }

      for(size_t i = 0; i < Ngrid + 1; ++i)
      {
        std::print(movieStream, "{} {} {}", static_cast<double>(i) * dx, V[i], Pt[i]);
        for(size_t j = 0; j < Ncomp; ++j)
        {
          std::print(movieStream, " {} {} {} {} {} {}", Q[i * Ncomp + j], Qeq[i * Ncomp + j], P[i * Ncomp + j],
              P[i * Ncomp + j] / (p_total * components[j].molFraction), Dpdt[i * Ncomp + j], Dqdt[i * Ncomp + j]);
        }
        std::print(movieStream, "\n");
      }
      std::print(movieStream, "\n\n");
    }


    if(step % printEvery == 0)
    {
      std::print(stream, "Timestep {}, time: {} [s]\n", std::to_string(step), std::to_string(t));
      std::print(stream, "    Average number of IAST steps: {}\n", std::to_string(static_cast<double>(iastPerformance.first)/
                   static_cast<double>(iastPerformance.second)));
    }

    // check if we can set the expected end-time based on 10% longer time than when all 
    // adorbed mol-fractions are smaller than 1% of unity
    if(autoSteps)
    {
      double tolerance = 0.0;
      for(size_t i = 0; i < Ngrid + 1; ++i)
      {
        for(size_t j = 0; j < Ncomp; ++j)
        {
          tolerance = std::max(tolerance, std::abs((P[i * Ncomp + j] / (p_total * components[j].molFraction)) - 1.0));
        }
      }
      // consider 1% as being visibily indistinguishable from 'converged'
      // use a 10% longer time for display purposes
      if(tolerance < 0.01)
      {
        std::print(stream, "\nConvergence criteria reached, running 10% longer\n\n");
        Nsteps = static_cast<size_t>(1.1 * static_cast<double>(step));
        autoSteps = false;
      }
    }

    // SSP-RK Step 1
    // ======================================================================

    // calculate the derivatives Dq/dt and Dp/dt based on Qeq, Q, V, and P
    computeFirstDerivatives(Dqdt,Dpdt,Qeq,Q,V,P);

    // Dqdt and Dpdt are calculated at old time step
    // make estimate for the new loadings and new gas phase partial pressures
    // first iteration is made using the Explicit Euler scheme
    for(size_t i = 0; i < Ngrid + 1; ++i)
    {
      for(size_t j = 0; j < Ncomp; ++j)
      {
         Qnew[i * Ncomp + j] = Q[i * Ncomp + j] + dt * Dqdt[i * Ncomp + j];
         Pnew[i * Ncomp + j] = P[i * Ncomp + j] + dt * Dpdt[i * Ncomp + j];
      }
    }

    computeEquilibriumLoadings(mixture);

    computeVelocity();

    // SSP-RK Step 2
    // ======================================================================
    
    // calculate new derivatives at new (current) timestep
    // calculate the derivatives Dq/dt and Dp/dt based on Qeq, Q, V, and P at new (current) timestep
    computeFirstDerivatives(Dqdtnew,Dpdtnew,Qeqnew,Qnew,Vnew,Pnew);

    for(size_t i = 0; i < Ngrid + 1; ++i)
    {
      for(size_t j = 0; j < Ncomp; ++j)
      {
         Qnew[i * Ncomp + j] = 0.75 * Q[i * Ncomp + j] + 0.25 * Qnew[i * Ncomp + j] +
                               0.25 * dt * Dqdtnew[i * Ncomp + j];
         Pnew[i * Ncomp + j] = 0.75 * P[i * Ncomp + j] + 0.25 * Pnew[i * Ncomp + j] +
                               0.25 * dt * Dpdtnew[i * Ncomp + j];
      }
    }

    computeEquilibriumLoadings(mixture);

    computeVelocity();

    // SSP-RK Step 3
    // ======================================================================

    for(size_t i = 0; i < Ngrid + 1; ++i)
    {
      for(size_t j = 0; j < Ncomp; ++j)
      {
         Qnew[i * Ncomp + j] = (1.0/3.0) * Q[i * Ncomp + j] + (2.0/3.0) * Qnew[i * Ncomp + j] +
                               (2.0/3.0) * dt * Dqdtnew[i * Ncomp + j];
         Pnew[i * Ncomp + j] = (1.0/3.0) * P[i * Ncomp + j] + (2.0/3.0) * Pnew[i * Ncomp + j] +
                               (2.0/3.0) * dt * Dpdtnew[i * Ncomp + j];
      }
    }

    computeEquilibriumLoadings(mixture);

    computeVelocity();

    // update to the new time step
    std::copy(Qnew.begin(), Qnew.end(), Q.begin());
    std::copy(Pnew.begin(), Pnew.end(), P.begin());
    std::copy(Qeqnew.begin(), Qeqnew.end(), Qeq.begin());
    std::copy(Vnew.begin(), Vnew.end(), V.begin());
  }

  std::print(stream, "Final timestep {}, time: {}\n", std::to_string(Nsteps), std::to_string(dt * static_cast<double>(Nsteps)));
}

void Breakthrough::computeEquilibriumLoadings(MixturePrediction &mixture)
{
    // calculate new equilibrium loadings Qeqnew corresponding to the new timestep
  for(size_t i = 0; i < Ngrid + 1; ++i)
  {
    // estimation of total pressure Pt at each grid point from partial pressures
    Pt[i] = 0.0;
    for(size_t j = 0; j < Ncomp; ++j)
    {
      Pt[i] += Pnew[i * Ncomp + j];
    }

    // compute gas-phase mol-fractions
    // force the gas-phase mol-fractions to be positive and normalized
    double sum = 0.0;
    for(size_t j = 0; j < Ncomp; ++j)
    {
      Yi[j] = std::max(Pnew[i * Ncomp + j], 0.0);
      sum += Yi[j];
    }
    for(size_t j = 0; j < Ncomp; ++j)
    {
      Yi[j] /= sum;
    }


    // use Yi and Pt[i] to compute the loadings in the adsorption mixture
    iastPerformance += mixture.predictMixture(Yi, Pt[i], Xi, Ni, &cachedP0[i * Ncomp], cachedPsi[i]);

    for(size_t j = 0; j < Ncomp; ++j)
    {
      Qeqnew[i * Ncomp + j] = Ni[j];
    }
  }

  // check the total pressure at the outlet, it should not be negative
  if (Pt[0] + dptdx * L < 0.0)
  {
    throw std::runtime_error("Error: pressure gradient is too large (negative outlet pressure)\n");
  }
}


// calculate the derivatives Dq/dt and Dp/dt along the column
void Breakthrough::computeFirstDerivatives(std::vector<double> &dqdt,
                                           std::vector<double> &dpdt,
                                           const std::vector<double> &q_eq,
                                           const std::vector<double> &q,
                                           const std::vector<double> &v,
                                           const std::vector<double> &p)
{
  double idx = 1.0 / dx;
  double idx2 = 1.0 / (dx * dx);

  // first gridpoint
  for(size_t j = 0; j < Ncomp; ++j)
  {
    dqdt[0 * Ncomp + j] = components[j].massTransferCoefficient * (q_eq[0 * Ncomp + j] - q[0 * Ncomp + j]);
    dpdt[0 * Ncomp + j] = 0.0;
  }

  // middle gridpoints
  for(size_t i = 1; i < Ngrid; i++)
  {
    for(size_t j = 0; j < Ncomp; ++j)
    {
      dqdt[i * Ncomp + j] = components[j].massTransferCoefficient * (q_eq[i * Ncomp + j] - q[i * Ncomp + j]);
      dpdt[i * Ncomp + j] = (v[i - 1] * p[(i - 1) * Ncomp + j] - v[i] * p[i * Ncomp + j]) * idx
                            + components[j].axialDispersionCoefficient * (p[(i + 1) * Ncomp + j] - 2.0 * p[i * Ncomp + j] + p[(i - 1) * Ncomp + j]) * idx2
                            - prefactor[j] * (q_eq[i * Ncomp + j] - q[i * Ncomp + j]);
    }
  }

  // last gridpoint
  for(size_t j = 0; j < Ncomp; ++j)
  {
    dqdt[Ngrid * Ncomp + j] = components[j].massTransferCoefficient * (q_eq[Ngrid * Ncomp + j] - q[Ngrid * Ncomp + j]);
    dpdt[Ngrid * Ncomp + j] = (v[Ngrid - 1] * p[(Ngrid - 1) * Ncomp + j] - v[Ngrid] * P[Ngrid * Ncomp + j]) * idx
                              + components[j].axialDispersionCoefficient * (p[(Ngrid - 1) * Ncomp + j] - p[Ngrid * Ncomp + j]) * idx2
                              - prefactor[j] * (q_eq[Ngrid * Ncomp + j] - q[Ngrid * Ncomp + j]);
  }

}

// calculate new velocity Vnew from Qnew, Qeqnew, Pnew, Pt
void Breakthrough::computeVelocity()
{
  double idx2 = 1.0 / (dx * dx);

  // first grid point
  Vnew[0] = v_in;
 
  // middle gridpoints
  for(size_t i = 1; i < Ngrid; ++i)  
  {
    // sum = derivative at the actual gridpoint i
    double sum = 0.0;
    for(size_t j = 0; j < Ncomp; ++j)
    {
      sum = sum - prefactor[j] * (Qeqnew[i * Ncomp + j] - Qnew[i * Ncomp + j]) +
            components[j].axialDispersionCoefficient * (Pnew[(i - 1) * Ncomp + j] - 2.0 * Pnew[i * Ncomp + j] + Pnew[(i + 1) * Ncomp + j]) * idx2;
    }
  
    // explicit version
    Vnew[i] = Vnew[i - 1] + dx * (sum - Vnew[i - 1] * dptdx) / Pt[i];
  }
  
  // last grid point
  double sum = 0.0;
  for(size_t j = 0; j < Ncomp; ++j)
  {
    sum = sum - prefactor[j] * (Qeqnew[Ngrid * Ncomp + j] - Qnew[Ngrid * Ncomp + j]) +
          components[j].axialDispersionCoefficient * (Pnew[(Ngrid - 1) * Ncomp + j] - Pnew[Ngrid * Ncomp + j]) * idx2;
  }
  
  // explicit version
  Vnew[Ngrid] = Vnew[Ngrid-1] + dx * (sum - Vnew[Ngrid - 1] * dptdx) / Pt[Ngrid];
}

void Breakthrough::createPlotScript(std::string directoryName)
{
  std::ofstream stream(directoryName + "plot_breakthrough");

  std::print(stream, "set encoding utf8\n");
  std::print(stream, "set xlabel 'Dimensionless time, {{/Helvetica-Italic τ}}={{/Helvetica-Italic tv/L}} / [-]' font 'Helvetica,18'\n");
  std::print(stream, "set ylabel 'Concentration exit gas, {{/Helvetica-Italic c}}_i/{{/Helvetica-Italic c}}_{{i,0}} / [-]' offset 0.0,0 font 'Helvetica,18'\n");
  std::print(stream, "set bmargin 4\n");
  std::print(stream, "set yrange[0:]\n");

  std::print(stream, "set key outside top center horizontal samplen 2.5 height 0.5 spacing 1.5 font 'Helvetica, 10'\n");
  std::print(stream, "set key title '{} {{/:Italic T}}={} K, {{/:Italic p_t}}={} kPa'\n", displayName, T, p_total * 1e-3);

  std::print(stream, "set output 'breakthrough_dimensionless.pdf'\n");
  std::print(stream, "set term pdf color solid\n");

  std::print(stream, "ev=1\n");
  std::print(stream, "plot \\\n");
  for (size_t i = 0; i < Ncomp; i++)
  {
    std::string fileName = "component_" + std::to_string(components[i].componentId) + "_" + components[i].name + ".data";
    std::print(stream, "    '{}' us ($1):($3) every ev title '{} (y_i={})' with li lw 2{}\n",
               fileName, components[i].name, components[i].molFraction, (i < Ncomp - 1 ? ",\\" : ""));
  }
  std::print(stream, "set output 'breakthrough.pdf'\n");
  std::print(stream, "set xlabel 'Time, {{/Helvetica-Italic t}} / [min.]' font 'Helvetica,18'\n");
  std::print(stream, "plot \\\n");
  for (size_t i = 0; i < Ncomp; i++)
  {
    std::string fileName = "component_" + std::to_string(components[i].componentId) + "_" + components[i].name + ".data";
    std::print(stream, "    '{}' us ($2):($3) every ev title '{} (y_i={})' with li lw 2{}\n",
               fileName, components[i].name, components[i].molFraction, (i < Ncomp - 1 ? ",\\" : ""));
  }
}

void Breakthrough::createMovieScripts(std::string directoryName)
{
  std::ofstream makeMovieStream(directoryName + "make_movies");
  std::print(makeMovieStream, "./make_movie_V \"$@\"\n");
  std::print(makeMovieStream, "./make_movie_Pt \"$@\"\n");
  std::print(makeMovieStream, "./make_movie_Q \"$@\"\n");
  std::print(makeMovieStream, "./make_movie_Qeq \"$@\"\n");
  std::print(makeMovieStream, "./make_movie_P \"$@\"\n");
  std::print(makeMovieStream, "./make_movie_Pnorm \"$@\"\n");
  std::print(makeMovieStream, "./make_movie_Dpdt \"$@\"\n");
  std::print(makeMovieStream, "./make_movie_Dqdt \"$@\"\n");

  std::filesystem::path path{directoryName + "make_movies"};
  std::filesystem::permissions(path, std::filesystem::perms::owner_exec, std::filesystem::perm_options::add);

  createMovieScriptColumnV(directoryName);
  createMovieScriptColumnPt(directoryName);
  createMovieScriptColumnQ(directoryName);
  createMovieScriptColumnQeq(directoryName);
  createMovieScriptColumnP(directoryName);
  createMovieScriptColumnDpdt(directoryName);
  createMovieScriptColumnDqdt(directoryName);
  createMovieScriptColumnPnormalized(directoryName);
}

// -crf 18: the range of the CRF scale is 0–51, where 0 is lossless, 23 is the default, 
//          and 51 is worst quality possible; 18 is visually lossless or nearly so.
// -pix_fmt yuv420p: needed on apple devices
std::string movieScriptTemplate(std::string s)
{
  std::ostringstream stream;

  std::print(stream, "rm -f column_movie_{}.mp4\n", s);
  std::print(stream, "every=1\n");
  std::print(stream, "format='-c:v libx265 -tag:v hvc1'\n");
  std::print(stream, "width=800\n");
  std::print(stream, "height=600\n");
  std::print(stream, "quality=18\n");
  std::print(stream, "while getopts e:w:h:q:l flag\n");
  std::print(stream, "do\n");
  std::print(stream, "    case \"${{flag}}\" in\n");
  std::print(stream, "        e) every=${{OPTARG}};;\n");
  std::print(stream, "        w) width=${{OPTARG}};;\n");
  std::print(stream, "        h) height=${{OPTARG}};;\n");
  std::print(stream, "        q) quality=${{OPTARG}};;\n");
  std::print(stream, "        l) format=\"-c:v libx264\";;\n");
  std::print(stream, "    esac\n");
  std::print(stream, "done\n");
  std::print(stream, "gnuplot -c plot_column_{} $every $width $height | ffmpeg -f png_pipe -s:v \"${{width}},${{height}}\" -i pipe: $format -pix_fmt yuv420p -crf $quality -c:a aac column_movie_{}.mp4\n", s, s);
  return stream.str();
}

void Breakthrough::createMovieScriptColumnV(std::string directoryName)
{
  std::ofstream makeMovieStream(directoryName + "make_movie_V");
  makeMovieStream << movieScriptTemplate("V");

  std::filesystem::path path{directoryName + "make_movie_V"};
  std::filesystem::permissions(path, std::filesystem::perms::owner_exec, std::filesystem::perm_options::add);

  std::ofstream stream(directoryName + "plot_column_V");

  std::print(stream, "set encoding utf8\n");
  std::print(stream, "set terminal pngcairo size ARG2,ARG3 enhanced font 'Helvetica,10'\n");
  std::print(stream, "set xlabel 'Adsorber position / [m]' font 'Helvetica,18'\n");
  std::print(stream, "set ylabel 'Interstitial velocity, {{/Helvetica-Italic v}} / [m/s]' offset 0.0,0 font 'Helvetica,18'\n");
  std::print(stream, "set key outside top center horizontal samplen 2.5 height 0.5 spacing 1.5 font 'Helvetica, 10'\n");
  std::print(stream, "set bmargin 4\n");
  std::print(stream, "set title '{} {{/:Italic T}}={} K, {{/:Italic p_t}}={} kPa'\n", displayName, T, p_total * 1e-3);
  std::print(stream, "stats 'column.data' us 2 nooutput\n");
  std::print(stream, "max=STATS_max\n");
  std::print(stream, "stats 'column.data' us 1 nooutput\n");
  std::print(stream, "set xrange[0:STATS_max]\n");
  std::print(stream, "set yrange[0:1.1*max]\n");
  std::print(stream, "ev=int(ARG1)\n");
  std::print(stream, "do for [i=0:int((STATS_blocks-2)/ev)] {{\n");
  std::print(stream, "  plot \\\n");
  std::print(stream, "    'column.data' us 1:2 index ev*i notitle with li lw 2,\\\n");
  std::print(stream, "    'column.data' us 1:2 index ev*i notitle with po\n");
  std::print(stream, "}}\n");
}


void Breakthrough::createMovieScriptColumnPt(std::string directoryName)
{
  std::ofstream makeMovieStream(directoryName + "make_movie_Pt");
  makeMovieStream << movieScriptTemplate("Pt");

  std::filesystem::path path{directoryName + "make_movie_Pt"};
  std::filesystem::permissions(path, std::filesystem::perms::owner_exec, std::filesystem::perm_options::add);

  std::ofstream stream(directoryName + "plot_column_Pt");

  std::print(stream, "set encoding utf8\n");
  std::print(stream, "set terminal pngcairo size ARG2,ARG3 enhanced font 'Helvetica,10'\n");
  std::print(stream, "set xlabel 'Adsorber position / [m]' font 'Helvetica,18'\n");
  std::print(stream, "set ylabel 'Total Pressure, {{/Helvetica-Italic p_t}} / [Pa]' offset 0.0,0 font 'Helvetica,18'\n");
  std::print(stream, "set key outside top center horizontal samplen 2.5 height 0.5 spacing 1.5 font 'Helvetica, 10'\n");
  std::print(stream, "set bmargin 4\n");
  std::print(stream, "set title '{} {{/:Italic T}}={} K, {{/:Italic p_t}}={} kPa'\n", displayName, T, p_total * 1e-3);
  std::print(stream, "stats 'column.data' us 3 nooutput\n");
  std::print(stream, "max=STATS_max\n");
  std::print(stream, "stats 'column.data' us 1 nooutput\n");
  std::print(stream, "set xrange[0:STATS_max]\n");
  std::print(stream, "set yrange[0:1.1*max]\n");
  std::print(stream, "ev=int(ARG1)\n");
  std::print(stream, "do for [i=0:int((STATS_blocks-2)/ev)] {{\n");
  std::print(stream, "  plot \\\n");
  std::print(stream, "    'column.data' us 1:3 index ev*i notitle with li lw 2,\\\n");
  std::print(stream, "    'column.data' us 1:3 index ev*i notitle with po\n");
  std::print(stream, "}}\n");
}

void Breakthrough::createMovieScriptColumnQ(std::string directoryName)
{
  std::ofstream makeMovieStream(directoryName + "make_movie_Q");
  makeMovieStream << movieScriptTemplate("Q");

  std::filesystem::path path{directoryName + "make_movie_Q"};
  std::filesystem::permissions(path, std::filesystem::perms::owner_exec, std::filesystem::perm_options::add);

  std::ofstream stream(directoryName + "plot_column_Q");

  std::print(stream, "set encoding utf8\n");
  std::print(stream, "set terminal pngcairo size ARG2,ARG3 enhanced font 'Helvetica,10'\n");
  std::print(stream, "set xlabel 'Adsorber position / [m]' font 'Helvetica,18'\n");
  std::print(stream, "set ylabel 'Concentration, {{/Helvetica-Italic c}}_i / [mol/kg]' offset 0.0,0 font 'Helvetica,18'\n");
  std::print(stream, "set key outside top center horizontal samplen 2.5 height 0.5 spacing 1.5 font 'Helvetica, 10'\n");
  std::print(stream, "set bmargin 4\n");
  std::print(stream, "set key title '{} {{/:Italic T}}={} K, {{/:Italic p_t}}={} kPa'\n", displayName, T, p_total * 1e-3);
  std::print(stream, "stats 'column.data' nooutput\n");
  std::print(stream, "max = 0.0;\n");
  std::print(stream, "do for [i=4:STATS_columns:6] {{\n");
  std::print(stream, "  stats 'column.data' us i nooutput\n");
  std::print(stream, "  if (max<STATS_max) {{\n");
  std::print(stream, "    max=STATS_max\n");
  std::print(stream, "  }}\n");
  std::print(stream, "}}\n");
  std::print(stream, "stats 'column.data' us 1 nooutput\n");
  std::print(stream, "set xrange[0:STATS_max]\n");
  std::print(stream, "set yrange[0:1.1*max]\n");
  std::print(stream, "ev=int(ARG1)\n");
  std::print(stream, "do for [i=0:int((STATS_blocks-2)/ev)] {{\n");
  std::print(stream, "  plot \\\n");
  for (size_t i = 0; i < Ncomp; i++)
  {
    std::print(stream, "    'column.data' us 1:{} index ev*i notitle with li lw 2,\\\n", std::to_string(4 + i * 6));
  }
  for (size_t i = 0; i < Ncomp; i++)
  {
    std::print(stream, "    'column.data' us 1:{} index ev*i title '{} (y_i={})' with po{}\n",
               std::to_string(4 + i * 6), components[i].name, components[i].molFraction, (i < Ncomp - 1 ? ",\\" : ""));
  }
  std::print(stream, "}}\n");
}

void Breakthrough::createMovieScriptColumnQeq(std::string directoryName)
{
  std::ofstream makeMovieStream(directoryName + "make_movie_Qeq");
  makeMovieStream << movieScriptTemplate("Qeq");

  std::filesystem::path path{directoryName + "make_movie_Qeq"};
  std::filesystem::permissions(path, std::filesystem::perms::owner_exec, std::filesystem::perm_options::add);

  std::ofstream stream(directoryName + "plot_column_Qeq");

  std::print(stream, "set encoding utf8\n");
  std::print(stream, "set terminal pngcairo size ARG2,ARG3 enhanced font 'Helvetica,10'\n");
  std::print(stream, "set xlabel 'Adsorber position / [m]' font 'Helvetica,18'\n");
  std::print(stream, "set ylabel 'Concentration, {{/Helvetica-Italic c}}_i / [mol/kg]' offset 0.0,0 font 'Helvetica,18'\n");
  std::print(stream, "set key outside top center horizontal samplen 2.5 height 0.5 spacing 1.5 font 'Helvetica, 10'\n");
  std::print(stream, "set bmargin 4\n");
  std::print(stream, "set key title '{} {{/:Italic T}}={} K, {{/:Italic p_t}}={} kPa'\n", displayName, T, p_total * 1e-3);
  std::print(stream, "stats 'column.data' nooutput\n");
  std::print(stream, "max = 0.0;\n");
  std::print(stream, "do for [i=5:STATS_columns:6] {{\n");
  std::print(stream, "  stats 'column.data' us i nooutput\n");
  std::print(stream, "  if (max<STATS_max) {{\n");
  std::print(stream, "    max=STATS_max\n");
  std::print(stream, "  }}\n");
  std::print(stream, "}}\n");
  std::print(stream, "stats 'column.data' us 1 nooutput\n");
  std::print(stream, "set xrange[0:STATS_max]\n");
  std::print(stream, "set yrange[0:1.1*max]\n");
  std::print(stream, "ev=int(ARG1)\n");
  std::print(stream, "do for [i=0:int((STATS_blocks-2)/ev)] {{\n");
  std::print(stream, "  plot \\\n");
  for (size_t i = 0; i < Ncomp; i++)
  {
    std::print(stream, "    'column.data' us 1:{} index ev*i notitle with li lw 2,\\\n", std::to_string(5 + i * 6));
  }
  for (size_t i = 0; i < Ncomp; i++)
  {
    std::print(stream, "    'column.data' us 1:{} index ev*i title '{} (y_i={})' with po{}\n",
               std::to_string(5 + i * 6), components[i].name, components[i].molFraction, (i < Ncomp - 1 ? ",\\" : ""));
  }
  std::print(stream, "}}\n");
}

void Breakthrough::createMovieScriptColumnP(std::string directoryName)
{
  std::ofstream makeMovieStream(directoryName + "make_movie_P");
  makeMovieStream << movieScriptTemplate("P");

  std::filesystem::path path{directoryName + "make_movie_P"};
  std::filesystem::permissions(path, std::filesystem::perms::owner_exec, std::filesystem::perm_options::add);

  std::ofstream stream(directoryName + "plot_column_P");

  std::print(stream, "set encoding utf8\n");
  std::print(stream, "set terminal pngcairo size ARG2,ARG3 enhanced font 'Helvetica,10'\n");
  std::print(stream, "set xlabel 'Adsorber position / [m]' font 'Helvetica,18'\n");
  std::print(stream, "set ylabel 'Partial pressure, {{/Helvetica-Italic p}}_i / [Pa]' offset 0.0,0 font 'Helvetica,18'\n");
  std::print(stream, "set key outside top center horizontal samplen 2.5 height 0.5 spacing 1.5 font 'Helvetica, 10'\n");
  std::print(stream, "set bmargin 4\n");
  std::print(stream, "set key title '{} {{/:Italic T}}={} K, {{/:Italic p_t}}={} kPa'\n", displayName, T, p_total * 1e-3);
  std::print(stream, "stats 'column.data' nooutput\n");
  std::print(stream, "max = 0.0;\n");
  std::print(stream, "do for [i=6:STATS_columns:6] {{\n");
  std::print(stream, "  stats 'column.data' us i nooutput\n");
  std::print(stream, "  if (max<STATS_max) {{\n");
  std::print(stream, "    max=STATS_max\n");
  std::print(stream, "  }}\n");
  std::print(stream, "}}\n");
  std::print(stream, "stats 'column.data' us 1 nooutput\n");
  std::print(stream, "set xrange[0:STATS_max]\n");
  std::print(stream, "set yrange[0:1.1*max]\n");
  std::print(stream, "ev=int(ARG1)\n");
  std::print(stream, "do for [i=0:int((STATS_blocks-2)/ev)] {{\n");
  std::print(stream, "  plot \\\n");
  for (size_t i = 0; i < Ncomp; i++)
  {
    std::print(stream, "    'column.data' us 1:{} index ev*i notitle with li lw 2,\\\n", std::to_string(6 + i * 6));
  }
  for (size_t i = 0; i < Ncomp; i++)
  {
    std::print(stream, "    'column.data' us 1:{} index ev*i title '{} (y_i={})' with po{}\n",
               std::to_string(6 + i * 6), components[i].name, components[i].molFraction, (i < Ncomp - 1 ? ",\\" : ""));;
  }
  std::print(stream, "}}\n");
}

void Breakthrough::createMovieScriptColumnPnormalized(std::string directoryName)
{
  std::ofstream makeMovieStream(directoryName + "make_movie_Pnorm");
  makeMovieStream << movieScriptTemplate("Pnorm");

  std::filesystem::path path{directoryName + "make_movie_Pnorm"};
  std::filesystem::permissions(path, std::filesystem::perms::owner_exec, std::filesystem::perm_options::add);

  std::ofstream stream(directoryName + "plot_column_Pnorm");

  std::print(stream, "set encoding utf8\n");
  std::print(stream, "set terminal pngcairo size ARG2,ARG3 enhanced font 'Helvetica,10'\n");
  std::print(stream, "set xlabel 'Adsorber position / [m]' font 'Helvetica,18'\n");
  std::print(stream, "set ylabel 'Partial pressure, {{/Helvetica-Italic p}}_i / [-]' offset 0.0,0 font 'Helvetica,18'\n");
  std::print(stream, "set key outside top center horizontal samplen 2.5 height 0.5 spacing 1.5 font 'Helvetica, 10'\n");
  std::print(stream, "set bmargin 4\n");
  std::print(stream, "set key title '{} {{/:Italic T}}={} K, {{/:Italic p_t}}={} kPa'\n", displayName, T, p_total * 1e-3);
  std::print(stream, "stats 'column.data' nooutput\n");
  std::print(stream, "max = 0.0;\n");
  std::print(stream, "do for [i=7:STATS_columns:6] {{\n");
  std::print(stream, "  stats 'column.data' us i nooutput\n");
  std::print(stream, "  if (max<STATS_max) {{\n");
  std::print(stream, "    max=STATS_max\n");
  std::print(stream, "  }}\n");
  std::print(stream, "}}\n");
  std::print(stream, "stats 'column.data' us 1 nooutput\n");
  std::print(stream, "set xrange[0:STATS_max]\n");
  std::print(stream, "set yrange[0:1.1*max]\n");
  std::print(stream, "ev=int(ARG1)\n");
  std::print(stream, "do for [i=0:int((STATS_blocks-2)/ev)] {{\n");
  std::print(stream, "  plot \\\n");
  for (size_t i = 0; i < Ncomp; i++)
  {
    std::print(stream, "    'column.data' us 1:{} index ev*i notitle with li lw 2,\\\n", std::to_string(7 + i * 6));
  }
  for (size_t i = 0; i < Ncomp; i++)
  {
    std::print(stream, "    'column.data' us 1:{} index ev*i title '{} (y_i={})' with po{}\n",
               std::to_string(7 + i * 6), components[i].name, components[i].molFraction, (i < Ncomp - 1 ? ",\\" : ""));
  }
  std::print(stream, "}}\n");
}

void Breakthrough::createMovieScriptColumnDpdt(std::string directoryName)
{
  std::ofstream makeMovieStream(directoryName + "make_movie_Dpdt");
  makeMovieStream << movieScriptTemplate("Dpdt");

  std::filesystem::path path{directoryName + "make_movie_Dpdt"};
  std::filesystem::permissions(path, std::filesystem::perms::owner_exec, std::filesystem::perm_options::add);

  std::ofstream stream(directoryName + "plot_column_Dpdt");

  std::print(stream, "set encoding utf8\n");
  std::print(stream, "set terminal pngcairo size ARG2,ARG3 enhanced font 'Helvetica,10'\n");
  std::print(stream, "set xlabel 'Adsorber position / [m]' font 'Helvetica,18'\n");
  std::print(stream, "set ylabel 'Pressure derivative, {{/Helvetica-Italic dp_/dt}} / [Pa/s]' offset 0.0,0 font 'Helvetica,18'\n");
  std::print(stream, "set key outside top center horizontal samplen 2.5 height 0.5 spacing 1.5 font 'Helvetica, 10'\n");
  std::print(stream, "set bmargin 4\n");
  std::print(stream, "set key title '{} {{/:Italic T}}={} K, {{/:Italic p_t}}={} kPa'\n", displayName, T, p_total * 1e-3);
  std::print(stream, "stats 'column.data' nooutput\n");
  std::print(stream, "max = -1e10;\n");
  std::print(stream, "min = 1e10;\n");
  std::print(stream, "do for [i=8:STATS_columns:6] {{\n");
  std::print(stream, "  stats 'column.data' us i nooutput\n");
  std::print(stream, "  if (STATS_max>max) {{\n");
  std::print(stream, "    max=STATS_max\n");
  std::print(stream, "  }}\n");
  std::print(stream, "  if (STATS_min<min) {{\n");
  std::print(stream, "    min=STATS_min\n");
  std::print(stream, "  }}\n");
  std::print(stream, "}}\n");
  std::print(stream, "stats 'column.data' us 1 nooutput\n");
  std::print(stream, "set xrange[0:STATS_max]\n");
  std::print(stream, "set yrange[1.1*min:1.1*max]\n");
  std::print(stream, "ev=int(ARG1)\n");
  std::print(stream, "do for [i=0:int((STATS_blocks-2)/ev)] {{\n");
  std::print(stream, "  plot \\\n");
  for (size_t i = 0; i < Ncomp; i++)
  {
    std::print(stream, "    'column.data' us 1:{} index ev*i notitle with li lw 2,\\\n", std::to_string(8 + i * 6));
  }
  for (size_t i = 0; i < Ncomp; i++)
  {
    std::print(stream, "    'column.data' us 1:{} index ev*i title '{} (y_i={})' with po{}\n",
               std::to_string(8 + i * 6), components[i].name, components[i].molFraction, (i < Ncomp - 1 ? ",\\" : ""));
  }
  std::print(stream, "}}\n");
}

void Breakthrough::createMovieScriptColumnDqdt(std::string directoryName)
{
  std::ofstream makeMovieStream(directoryName + "make_movie_Dqdt");
  makeMovieStream << movieScriptTemplate("Dqdt");

  std::filesystem::path path{directoryName + "make_movie_Dqdt"};
  std::filesystem::permissions(path, std::filesystem::perms::owner_exec, std::filesystem::perm_options::add);

  std::ofstream stream(directoryName + "plot_column_Dqdt");

  std::print(stream, "set encoding utf8\n");
  std::print(stream, "set terminal pngcairo size ARG2,ARG3 enhanced font 'Helvetica,10'\n");
  std::print(stream, "set xlabel 'Adsorber position / [m]' font 'Helvetica,18'\n");
  std::print(stream, "set ylabel 'Loading derivative, {{/Helvetica-Italic dq_i/dt}} / [mol/kg/s]' offset 0.0,0 font 'Helvetica,18'\n");
  std::print(stream, "set key outside top center horizontal samplen 2.5 height 0.5 spacing 1.5 font 'Helvetica, 10'\n");
  std::print(stream, "set bmargin 4\n");
  std::print(stream, "set key title '{} {{/:Italic T}}={} K, {{/:Italic p_t}}={} kPa'\n", displayName, T, p_total * 1e-3);
  std::print(stream, "stats 'column.data' nooutput\n");
  std::print(stream, "max = -1e10;\n");
  std::print(stream, "min = 1e10;\n");
  std::print(stream, "min = 10000000000000.0;\n");
  std::print(stream, "do for [i=9:STATS_columns:6] {{\n");
  std::print(stream, "  stats 'column.data' us i nooutput\n");
  std::print(stream, "  if (STATS_max>max) {{\n");
  std::print(stream, "    max=STATS_max\n");
  std::print(stream, "  }}\n");
  std::print(stream, "  if (STATS_min<min) {{\n");
  std::print(stream, "    min=STATS_min\n");
  std::print(stream, "  }}\n");
  std::print(stream, "}}\n");
  std::print(stream, "stats 'column.data' us 1 nooutput\n");
  std::print(stream, "set xrange[0:STATS_max]\n");
  std::print(stream, "set yrange[1.1*min:1.1*max]\n");
  std::print(stream, "ev=int(ARG1)\n");
  std::print(stream, "do for [i=0:int((STATS_blocks-2)/ev)] {{\n");
  std::print(stream, "  plot \\\n");
  for (size_t i = 0; i < Ncomp; i++)
  {
    std::print(stream, "    'column.data' us 1:{} index ev*i notitle with li lw 2,\\\n", std::to_string(9 + i * 6));
  }
  for (size_t i = 0; i < Ncomp; i++)
  {
  std::print(stream, "    'column.data' us 1:{} index ev*i title '{} (y_i={})' with po{}\n",
             std::to_string(9 + i * 6), components[i].name, components[i].molFraction, (i < Ncomp - 1 ? ",\\" : ""));
  }
  std::print(stream, "}}\n");
}
