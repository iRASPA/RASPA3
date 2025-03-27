module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <print>
#include <span>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>
#if defined(__has_include) && __has_include(<mdspan>)
#include <mdspan>
#endif
#endif

module breakthrough;

#ifndef USE_LEGACY_HEADERS
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
import <type_traits>;
import <print>;
import <mdspan>;
#endif

import stringutils;
import input_reader;
import component;
import system;
import simulationbox;
import mixture_prediction;
#if !(defined(__has_include) && __has_include(<mdspan>))
import mdspan;
#endif

// TODO: move std::span to std::mdarray in C++26
// TODO: move cachedP0 and cachedPsi to submdspan in C++26

const double R = 8.31446261815324;

inline double maxVectorDifference(const std::vector<double> &v, const std::vector<double> &w)
{
  if (v.empty() || w.empty()) return 0.0;
  if (v.size() != w.size()) throw std::runtime_error("Error: unequal vector size\n");

  double max = std::abs(v[0] - w[0]);
  for (size_t i = 1; i < v.size(); ++i)
  {
    double temp = std::abs(v[i] - w[i]);
    if (temp > max) max = temp;
  }
  return max;
}

// allow std::pairs to be added
template <typename T, typename U>
std::pair<T, U> operator+(const std::pair<T, U> &l, const std::pair<T, U> &r)
{
  return {l.first + r.first, l.second + r.second};
}
template <typename T, typename U>
std::pair<T, U> &operator+=(std::pair<T, U> &l, const std::pair<T, U> &r)
{
  l.first += r.first;
  l.second += r.second;
  return l;
}

Breakthrough::Breakthrough(System &system)
    : system(system),
      displayName(system.components.front().name),
      components(system.components),
      Ncomp(components.size()),
      Ngrid(system.columnNumberOfGridPoints),
      printEvery(10000),
      writeEvery(5000),
      T(system.temperature),
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
      mixture(system),
      prefactor(Ncomp),
      Yi(Ncomp),
      Xi(Ncomp),
      Ni(Ncomp),
      V(Ngrid + 1),
      Vnew(Ngrid + 1),
      Pt(Ngrid + 1),
      P_vector((Ngrid + 1) * Ncomp),
      P(P_vector.data(), Ngrid + 1, Ncomp),
      Pnew_vector((Ngrid + 1) * Ncomp),
      Pnew(Pnew_vector.data(), Ngrid + 1, Ncomp),
      Q_vector((Ngrid + 1) * Ncomp),
      Q(Q_vector.data(), Ngrid + 1, Ncomp),
      Qnew_vector((Ngrid + 1) * Ncomp),
      Qnew(Qnew_vector.data(), Ngrid + 1, Ncomp),
      Qeq_vector((Ngrid + 1) * Ncomp),
      Qeq(Qeq_vector.data(), Ngrid + 1, Ncomp),
      Qeqnew_vector((Ngrid + 1) * Ncomp),
      Qeqnew(Qeqnew_vector.data(), Ngrid + 1, Ncomp),
      Dpdt_vector((Ngrid + 1) * Ncomp),
      Dpdt(Dpdt_vector.data(), Ngrid + 1, Ncomp),
      Dpdtnew_vector((Ngrid + 1) * Ncomp),
      Dpdtnew(Dpdtnew_vector.data(), Ngrid + 1, Ncomp),
      Dqdt_vector((Ngrid + 1) * Ncomp),
      Dqdt(Dqdt_vector.data(), Ngrid + 1, Ncomp),
      Dqdtnew_vector((Ngrid + 1) * Ncomp),
      Dqdtnew(Dqdtnew_vector.data(), Ngrid + 1, Ncomp),
      cachedP0((Ngrid + 1) * Ncomp * system.maxIsothermTerms),
      cachedPsi((Ngrid + 1) * system.maxIsothermTerms)
{
  // precomputed factor for mass transfer
  for (size_t j = 0; j < Ncomp; ++j)
  {
    prefactor[j] = R * T * ((1.0 - epsilon) / epsilon) * rho_p * components[j].massTransferCoefficient;
  }

  // set P and Q to zero
  std::fill(P_vector.begin(), P_vector.end(), 0.0);
  std::fill(Q_vector.begin(), Q_vector.end(), 0.0);

  // initial pressure along the column
  std::vector<double> pt_init(Ngrid + 1);

  // set the initial total pressure along the column assuming the pressure gradient is constant
  for (size_t i = 0; i < Ngrid + 1; ++i)
  {
    pt_init[i] = p_total + dptdx * static_cast<double>(i) * dx;
  }

  // initialize the interstitial gas velocity in the column
  for (size_t i = 0; i <= Ngrid; ++i)
  {
    V[i] = v_in * p_total / pt_init[i];
  }

  // set the partial pressure of the carrier gas to the total initial pressure
  // for the column except for the entrance (i=0)
  for (size_t i = 1; i <= Ngrid; ++i)
  {
    P[i, system.carrierGasComponent] = pt_init[i];
  }

  // auto st = std::experimental::mdspan(P.data(), 2, 6);
  // std::experimental::mdspan<double, std::experimental::extents<size_t, 20, 10>> m(P.data(), 16, 32);

  // at the column entrance, the mol-fractions of the components in the gas phase are fixed
  // the partial pressures of the components at the entrance are the mol-fractions times the
  // total pressure
  for (size_t j = 0; j < Ncomp; ++j)
  {
    P[0, j] = p_total * components[j].molFraction;
  }

  // at the entrance: mol-fractions Yi are the gas-phase mol-fractions
  // for the column: the initial mol-fraction of the carrier-gas is 1, and 0 for the other components
  //
  // the K of the carrier gas is chosen as zero
  // so Qeq is zero for all components in the column after the entrance
  // only the values for Yi at the entrance are effected by adsorption
  for (size_t i = 0; i < Ngrid + 1; ++i)
  {
    double sum = 0.0;
    for (size_t j = 0; j < Ncomp; ++j)
    {
      Yi[j] = std::max(P[i, j] / pt_init[i], 0.0);
      sum += Yi[j];
    }
    for (size_t j = 0; j < Ncomp; ++j)
    {
      Yi[j] /= sum;
    }

    iastPerformance += mixture.predictMixture(Yi, pt_init[i], Xi, Ni, &cachedP0[i * Ncomp * system.maxIsothermTerms],
                                              &cachedPsi[i * system.maxIsothermTerms]);

    for (size_t j = 0; j < Ncomp; ++j)
    {
      Qeq[i, j] = Ni[j];
    }
  }

  for (size_t i = 0; i < Ngrid + 1; ++i)
  {
    Pt[i] = 0.0;
    for (size_t j = 0; j < Ncomp; ++j)
    {
      Pt[i] += std::max(0.0, P[i, j]);
    }
  }
}

std::string Breakthrough::writeHeader()
{
  std::ostringstream stream;

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

  std::print(stream, "Breakthrough settings\n");
  std::print(stream, "=======================================================\n");
  std::print(stream, "Number of time steps:          {}\n", Nsteps);
  std::print(stream, "Print every step:              {}\n", printEvery);
  std::print(stream, "Write data every step:         {}\n", writeEvery);
  std::print(stream, "\n\n");

  std::print(stream, "Integration details\n");
  std::print(stream, "=======================================================\n");
  std::print(stream, "Time step:                     {} [s]\n", dt);
  std::print(stream, "Number of column grid points:  {} [-]\n", Ngrid);
  std::print(stream, "Column spacing:                {} [m]\n", dx);
  std::print(stream, "\n\n");

  std::print(stream, "Component data\n");
  std::print(stream, "=======================================================\n");
  for (size_t i = 0; i < Ncomp; ++i)
  {
    // std::print(stream, components[i].printBreakthroughStatus());
    std::print(stream, "\n");
  }

  return stream.str();
}

void Breakthrough::run(std::ostream &stream)
{
  std::filesystem::create_directory("Breakthrough");
  std::filesystem::create_directory(std::format("Breakthrough/System_{}", system.systemId));

  // create the output files
  std::vector<std::ofstream> streams;
  for (size_t i = 0; i < Ncomp; i++)
  {
    std::string fileName = std::format("Breakthrough/System_{}/component_{}_{}.txt", system.systemId, std::to_string(i),
                                       components[i].name);
    streams.emplace_back(std::ofstream{fileName});
  }

  std::ofstream movieStream(std::format("Breakthrough/System_{}/column.txt", system.systemId));

  for (size_t step = 0; (step < Nsteps || autoSteps); ++step)
  {
    double t = static_cast<double>(step) * dt;

    // pulse boundary condition
    if (pulse == true)
    {
      if (t > tpulse)
      {
        for (size_t j = 0; j < Ncomp; ++j)
        {
          if (j == system.carrierGasComponent)
          {
            P[0, j] = p_total;
          }
          else
          {
            P[0, j] = 0.0;
          }
        }
      }
    }

    if (step % writeEvery == 0)
    {
      // write breakthrough output to files
      // column 1: dimensionless time
      // column 2: time [minutes]
      // column 3: normalized partial pressure
      for (size_t j = 0; j < Ncomp; ++j)
      {
        streams[j] << t * v_in / L << " " << t / 60.0 << " "
                   << P[Ngrid, j] / ((p_total + dptdx * L) * components[j].molFraction) << std::endl;
      }

      size_t column_nr = 1;
      movieStream << "# column " << column_nr++ << ": z  (column position)" << std::endl;
      movieStream << "# column " << column_nr++ << ": V  (velocity)" << std::endl;
      movieStream << "# column " << column_nr++ << ": Pt (total pressure)" << std::endl;
      for (size_t j = 0; j < Ncomp; ++j)
      {
        movieStream << "# column " << column_nr++ << ": component " << j << " Q     (loading)\n";
        movieStream << "# column " << column_nr++ << ": component " << j << " Qeq   (equilibrium loading)\n";
        movieStream << "# column " << column_nr++ << ": component " << j << " P     (partial pressure)\n";
        movieStream << "# column " << column_nr++ << ": component " << j << " Pnorm (normalized partial pressure)\n";
        movieStream << "# column " << column_nr++ << ": component " << j << " Dpdt  (derivative P with t)\n";
        movieStream << "# column " << column_nr++ << ": component " << j << " Dqdt  (derivative Q with t)\n";
      }

      for (size_t i = 0; i < Ngrid + 1; ++i)
      {
        movieStream << static_cast<double>(i) * dx << " ";
        movieStream << V[i] << " ";
        movieStream << Pt[i] << " ";
        for (size_t j = 0; j < Ncomp; ++j)
        {
          movieStream << Q[i, j] << " " << Qeq[i, j] << " " << P[i, j] << " "
                      << P[i, j] / (Pt[i] * components[j].molFraction) << " " << Dpdt[i, j] << " " << Dqdt[i, j] << " ";
        }
        movieStream << "\n";
      }
      movieStream << "\n\n";
    }

    if (step % printEvery == 0)
    {
      std::print(stream, "Timestep {}, time: {} [s]\n", std::to_string(step), std::to_string(t));
      std::print(
          stream, "    Average number of mixture-prediction steps: {}\n",
          std::to_string(static_cast<double>(iastPerformance.first) / static_cast<double>(iastPerformance.second)));
    }

    // check if we can set the expected end-time based on 10% longer time than when all
    // adorbed mol-fractions are smaller than 1% of unity
    if (autoSteps)
    {
      double tolerance = 0.0;
      for (size_t j = 0; j < Ncomp; ++j)
      {
        tolerance =
            std::max(tolerance, std::abs((P[Ngrid, j] / ((p_total + dptdx * L) * components[j].molFraction)) - 1.0));
      }

      // consider 1% as being visibily indistinguishable from 'converged'
      // use a 10% longer time for display purposes
      if (tolerance < 0.01)
      {
        std::print(stream, "\nConvergence criteria reached, running 10% longer\n\n");
        Nsteps = static_cast<size_t>(1.1 * static_cast<double>(step));
        autoSteps = false;
      }
    }

    // SSP-RK Step 1
    // ======================================================================

    // calculate the derivatives Dq/dt and Dp/dt based on Qeq, Q, V, and P
    computeFirstDerivatives(Dqdt, Dpdt, Qeq, Q, V, P);

    // Dqdt and Dpdt are calculated at old time step
    // make estimate for the new loadings and new gas phase partial pressures
    // first iteration is made using the Explicit Euler scheme
    for (size_t i = 0; i < Ngrid + 1; ++i)
    {
      for (size_t j = 0; j < Ncomp; ++j)
      {
        Qnew[i, j] = Q[i, j] + dt * Dqdt[i, j];
        Pnew[i, j] = P[i, j] + dt * Dpdt[i, j];
      }
    }

    computeEquilibriumLoadings();

    computeVelocity();

    // SSP-RK Step 2
    // ======================================================================

    // calculate new derivatives at new (current) timestep
    // calculate the derivatives Dq/dt and Dp/dt based on Qeq, Q, V, and P at new (current) timestep
    computeFirstDerivatives(Dqdtnew, Dpdtnew, Qeqnew, Qnew, Vnew, Pnew);

    for (size_t i = 0; i < Ngrid + 1; ++i)
    {
      for (size_t j = 0; j < Ncomp; ++j)
      {
        Qnew[i, j] = 0.75 * Q[i, j] + 0.25 * Qnew[i, j] + 0.25 * dt * Dqdtnew[i, j];
        Pnew[i, j] = 0.75 * P[i, j] + 0.25 * Pnew[i, j] + 0.25 * dt * Dpdtnew[i, j];
      }
    }

    computeEquilibriumLoadings();

    computeVelocity();

    // SSP-RK Step 3
    // ======================================================================

    // calculate new derivatives at new (current) timestep
    // calculate the derivatives Dq/dt and Dp/dt based on Qeq, Q, V, and P at new (current) timestep
    computeFirstDerivatives(Dqdtnew, Dpdtnew, Qeqnew, Qnew, Vnew, Pnew);

    for (size_t i = 0; i < Ngrid + 1; ++i)
    {
      for (size_t j = 0; j < Ncomp; ++j)
      {
        Qnew[i, j] = (1.0 / 3.0) * Q[i, j] + (2.0 / 3.0) * Qnew[i, j] + (2.0 / 3.0) * dt * Dqdtnew[i, j];
        Pnew[i, j] = (1.0 / 3.0) * P[i, j] + (2.0 / 3.0) * Pnew[i, j] + (2.0 / 3.0) * dt * Dpdtnew[i, j];
      }
    }

    computeEquilibriumLoadings();

    computeVelocity();

    // update to the new time step
    std::copy(Qnew_vector.begin(), Qnew_vector.end(), Q_vector.begin());
    std::copy(Pnew_vector.begin(), Pnew_vector.end(), P_vector.begin());
    std::copy(Qeqnew_vector.begin(), Qeqnew_vector.end(), Qeq_vector.begin());
    std::copy(Vnew.begin(), Vnew.end(), V.begin());
  }

  std::cout << "Final timestep " + std::to_string(Nsteps) +
                   ", time: " + std::to_string(dt * static_cast<double>(Nsteps)) + " [s]"
            << std::endl;
}

void Breakthrough::computeEquilibriumLoadings()
{
  // calculate new equilibrium loadings Qeqnew corresponding to the new timestep
  for (size_t i = 0; i < Ngrid + 1; ++i)
  {
    // estimation of total pressure Pt at each grid point from partial pressures
    Pt[i] = 0.0;
    for (size_t j = 0; j < Ncomp; ++j)
    {
      Pt[i] += std::max(0.0, Pnew[i, j]);
    }

    // compute gas-phase mol-fractions
    // force the gas-phase mol-fractions to be positive and normalized
    double sum = 0.0;
    for (size_t j = 0; j < Ncomp; ++j)
    {
      Yi[j] = std::max(Pnew[i, j], 0.0);
      sum += Yi[j];
    }
    for (size_t j = 0; j < Ncomp; ++j)
    {
      Yi[j] /= sum;
    }

    // use Yi and Pt[i] to compute the loadings in the adsorption mixture via mixture prediction
    iastPerformance += mixture.predictMixture(Yi, Pt[i], Xi, Ni, &cachedP0[i * Ncomp * system.maxIsothermTerms],
                                              &cachedPsi[i * system.maxIsothermTerms]);

    for (size_t j = 0; j < Ncomp; ++j)
    {
      Qeqnew[i, j] = Ni[j];
    }
  }

  // check the total pressure at the outlet, it should not be negative
  if (Pt[0] + dptdx * L < 0.0)
  {
    throw std::runtime_error("Error: pressure gradient is too large (negative outlet pressure)\n");
  }
}

// calculate the derivatives Dq/dt and Dp/dt along the column
void Breakthrough::computeFirstDerivatives(std::mdspan<double, std::dextents<size_t, 2>> &dqdt,
                                           std::mdspan<double, std::dextents<size_t, 2>> &dpdt,
                                           const std::mdspan<double, std::dextents<size_t, 2>> &q_eq,
                                           const std::mdspan<double, std::dextents<size_t, 2>> &q,
                                           const std::vector<double> &v,
                                           const std::mdspan<double, std::dextents<size_t, 2>> &p)
{
  double idx = 1.0 / dx;
  double idx2 = 1.0 / (dx * dx);

  // first gridpoint
  for (size_t j = 0; j < Ncomp; ++j)
  {
    dqdt[0, j] = components[j].massTransferCoefficient * (q_eq[0, j] - q[0, j]);
    dpdt[0, j] = 0.0;
  }

  // middle gridpoints
  for (size_t i = 1; i < Ngrid; i++)
  {
    for (size_t j = 0; j < Ncomp; ++j)
    {
      dqdt[i, j] = components[j].massTransferCoefficient * (q_eq[i, j] - q[i, j]);
      dpdt[i, j] = (v[i - 1] * p[i - 1, j] - v[i] * p[i, j]) * idx +
                   components[j].axialDispersionCoefficient * (p[i + 1, j] - 2.0 * p[i, j] + p[i - 1, j]) * idx2 -
                   prefactor[j] * (q_eq[i, j] - q[i, j]);
    }
  }

  // last gridpoint
  for (size_t j = 0; j < Ncomp; ++j)
  {
    dqdt[Ngrid, j] = components[j].massTransferCoefficient * (q_eq[Ngrid, j] - q[Ngrid, j]);
    dpdt[Ngrid, j] = (v[Ngrid - 1] * p[Ngrid - 1, j] - v[Ngrid] * p[Ngrid, j]) * idx +
                     components[j].axialDispersionCoefficient * (p[Ngrid - 1, j] - p[Ngrid, j]) * idx2 -
                     prefactor[j] * (q_eq[Ngrid, j] - q[Ngrid, j]);
  }
}

// calculate new velocity Vnew from Qnew, Qeqnew, Pnew, Pt
void Breakthrough::computeVelocity()
{
  double idx2 = 1.0 / (dx * dx);

  // first grid point
  Vnew[0] = v_in;

  // middle gridpoints
  for (size_t i = 1; i < Ngrid; ++i)
  {
    // sum = derivative at the actual gridpoint i
    double sum = 0.0;
    for (size_t j = 0; j < Ncomp; ++j)
    {
      sum = sum - prefactor[j] * (Qeqnew[i, j] - Qnew[i, j]) +
            components[j].axialDispersionCoefficient * (Pnew[i - 1, j] - 2.0 * Pnew[i, j] + Pnew[i + 1, j]) * idx2;
    }

    // explicit version
    Vnew[i] = Vnew[i - 1] + dx * (sum - Vnew[i - 1] * dptdx) / Pt[i];
  }

  // last grid point
  double sum = 0.0;
  for (size_t j = 0; j < Ncomp; ++j)
  {
    sum = sum - prefactor[j] * (Qeqnew[Ngrid, j] - Qnew[Ngrid, j]) +
          components[j].axialDispersionCoefficient * (Pnew[Ngrid - 1, j] - Pnew[Ngrid, j]) * idx2;
  }

  // explicit version
  Vnew[Ngrid] = Vnew[Ngrid - 1] + dx * (sum - Vnew[Ngrid - 1] * dptdx) / Pt[Ngrid];
}

void Breakthrough::print() const
{
  std::cout << "Column properties\n";
  std::cout << "=======================================================\n";
  std::cout << "Display-name:                          " << displayName << "\n";
  std::cout << "Temperature:                           " << T << " [K]\n";
  std::cout << "Column length:                         " << L << " [m]\n";
  std::cout << "Column void-fraction:                  " << epsilon << " [-]\n";
  std::cout << "Particle density:                      " << rho_p << " [kg/m^3]\n";
  std::cout << "Total pressure:                        " << p_total << " [Pa]\n";
  std::cout << "Pressure gradient:                     " << dptdx << " [Pa/m]\n";
  std::cout << "Column entrance interstitial velocity: " << v_in << " [m/s]\n";
  std::cout << "\n\n";

  std::cout << "Breakthrough settings\n";
  std::cout << "=======================================================\n";
  std::cout << "Number of time steps:          " << Nsteps << "\n";
  std::cout << "Print every step:              " << printEvery << "\n";
  std::cout << "Write data every step:         " << writeEvery << "\n";
  std::cout << "\n\n";

  std::cout << "Integration details\n";
  std::cout << "=======================================================\n";
  std::cout << "Time step:                     " << dt << " [s]\n";
  std::cout << "Number of column grid points:  " << Ngrid << "\n";
  std::cout << "Column spacing:                " << dx << " [m]\n";
  std::cout << "\n\n";

  std::cout << "Component data\n";
  std::cout << "=======================================================\n";
  std::cout << "maximum isotherm terms:        " << system.maxIsothermTerms << "\n";
  for (size_t i = 0; i < Ncomp; ++i)
  {
    // std::cout << components[i].print(i);
    std::cout << "\n";
  }
}

void Breakthrough::createPlotScript()
{
  std::filesystem::create_directory("Breakthrough");
  std::filesystem::create_directory(std::format("Breakthrough/System_{}", system.systemId));

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
  std::ofstream stream_graphs(std::format("Breakthrough/System_{}/make_graphs.bat", system.systemId));
  stream_graphs << "set PATH=%PATH%;C:\\Program Files\\gnuplot\\bin;"
                << "C:\\Program Files\\ffmpeg-master-latest-win64-gpl\\bin;C:\\Program Files\\ffmpeg\\bin\n";
  stream_graphs << "gnuplot.exe plot_breakthrough\n";

  std::filesystem::path path{std::format("Breakthrough/System_{}/make_graphs.bat", system.systemId)};
  std::filesystem::permissions(path, std::filesystem::perms::owner_exec, std::filesystem::perm_options::add);
#else
  std::ofstream stream_graphs(std::format("Breakthrough/System_{}/make_graphs", system.systemId));
  stream_graphs << "#!/bin/sh\n";
  stream_graphs << "cd -- \"$(dirname \"$0\")\"\n";
  stream_graphs << "gnuplot plot_breakthrough\n";

  std::filesystem::path path{std::format("Breakthrough/System_{}/make_graphs", system.systemId)};
  std::filesystem::permissions(path, std::filesystem::perms::owner_exec, std::filesystem::perm_options::add);
#endif

  std::ofstream stream(std::format("Breakthrough/System_{}/plot_breakthrough", system.systemId));
  stream << "set encoding utf8\n";
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
  stream << "set xlabel 'Dimensionless time, {/Arial-Italic τ}={/Arial-Italic tv/L} / [-]' font \"Arial,14\"\n";
  stream << "set ylabel 'Concentration exit gas, {/Arial-Italic c}_i/"
         << "{/Arial-Italic c}_{i,0} / [-]' offset 0.0,0 font \"Arial,14\"\n";
  stream << "set key outside top center horizontal samplen 2.5 height 0.5 spacing 1.5 font 'Arial, 10'\n";
#else
  stream << "set xlabel 'Dimensionless time, {/Helvetica-Italic τ}={/Helvetica-Italic tv/L} / [-]' "
            "font \"Helvetica,18\"\n";
  stream << "set ylabel 'Concentration exit gas, {/Helvetica-Italic c}_i/"
            "{/Helvetica-Italic c}_{i,0} / [-]' offset 0.0,0 font \"Helvetica,18\"\n";
  stream << "set key outside top center horizontal samplen 2.5 height 0.5 spacing 1.5 font 'Helvetica, 10'\n";
#endif
  stream << "set bmargin 4\n";
  stream << "set yrange[0:]\n";

  stream << "set key title '" << displayName << " {/:Italic T}=" << T << " K, {/:Italic p_t}=" << p_total * 1e-3
         << " kPa'\n";

  stream << "set output 'breakthrough_dimensionless.pdf'\n";
  stream << "set term pdf color solid\n";

  // colorscheme from book 'gnuplot in action', listing 12.7
  stream << "set linetype 1 pt 5 ps 1 lw 4 lc rgb '0xee0000'\n";
  stream << "set linetype 2 pt 7 ps 1 lw 4 lc rgb '0x008b00'\n";
  stream << "set linetype 3 pt 9 ps 1 lw 4 lc rgb '0x0000cd'\n";
  stream << "set linetype 4 pt 11 ps 1 lw 4 lc rgb '0xff3fb3'\n";
  stream << "set linetype 5 pt 13 ps 1 lw 4 lc rgb '0x00cdcd'\n";
  stream << "set linetype 6 pt 15 ps 1 lw 4 lc rgb '0xcd9b1d'\n";
  stream << "set linetype 7 pt  4 ps 1 lw 4 lc rgb '0x8968ed'\n";
  stream << "set linetype 8 pt  6 ps 1 lw 4 lc rgb '0x8b8b83'\n";
  stream << "set linetype 9 pt  8 ps 1 lw 4 lc rgb '0x00bb00'\n";
  stream << "set linetype 10 pt 10 ps 1 lw 4 lc rgb '0x1e90ff'\n";
  stream << "set linetype 11 pt 12 ps 1 lw 4 lc rgb '0x8b2500'\n";
  stream << "set linetype 12 pt 14 ps 1 lw 4 lc rgb '0x000000'\n";

  stream << "ev=1\n";
  stream << "plot \\\n";
  for (size_t i = 0; i < Ncomp; i++)
  {
    std::string fileName = "component_" + std::to_string(i) + "_" + components[i].name + ".txt";
    stream << "    " << "\"" << fileName << "\"" << " us ($1):($3) every ev" << " title \"" << components[i].name
           << " (y_i=" << components[i].molFraction << ")\"" << " with li lt " << i + 1 << (i < Ncomp - 1 ? ",\\" : "")
           << "\n";
  }
  stream << "set output 'breakthrough.pdf'\n";
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
  stream << "set xlabel 'Time, {/Arial-Italic t} / [min.]' font \"Arial,14\"\n";
#else
  stream << "set xlabel 'Time, {/Helvetica-Italic t} / [min.]' font \"Helvetica,18\"\n";
#endif
  stream << "plot \\\n";
  for (size_t i = 0; i < Ncomp; i++)
  {
    std::string fileName = "component_" + std::to_string(i) + "_" + components[i].name + ".txt";
    stream << "    " << "\"" << fileName << "\"" << " us ($2):($3) every ev" << " title \"" << components[i].name
           << " (y_i=" << components[i].molFraction << ")\"" << " with li lt " << i + 1 << (i < Ncomp - 1 ? ",\\" : "")
           << "\n";
  }
}

void Breakthrough::createMovieScripts()
{
  std::filesystem::create_directory("Breakthrough");
  std::filesystem::create_directory(std::format("Breakthrough/System_{}", system.systemId));

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
  std::ofstream makeMovieStream(std::format("Breakthrough/System_{}/make_movies.bat", system.systemId));
  makeMovieStream << "CALL make_movie_V.bat %1 %2 %3 %4\n";
  makeMovieStream << "CALL make_movie_Pt.bat %1 %2 %3 %4\n";
  makeMovieStream << "CALL make_movie_Q.bat %1 %2 %3 %4\n";
  makeMovieStream << "CALL make_movie_Qeq.bat %1 %2 %3 %4\n";
  makeMovieStream << "CALL make_movie_P.bat %1 %2 %3 %4\n";
  makeMovieStream << "CALL make_movie_Pnorm.bat %1 %2 %3 %4\n";
  makeMovieStream << "CALL make_movie_Dpdt.bat %1 %2 %3 %4\n";
  makeMovieStream << "CALL make_movie_Dqdt.bat %1 %2 %3 %4\n";

  std::filesystem::path path{std::format("Breakthrough/System_{}/make_movies.bat", system.systemId)};
  std::filesystem::permissions(path, std::filesystem::perms::owner_exec, std::filesystem::perm_options::add);
#else
  std::ofstream makeMovieStream(std::format("Breakthrough/System_{}/make_movies", system.systemId));
  makeMovieStream << "#!/bin/sh\n";
  makeMovieStream << "cd -- \"$(dirname \"$0\")\"\n";
  makeMovieStream << "./make_movie_V \"$@\"\n";
  makeMovieStream << "./make_movie_Pt \"$@\"\n";
  makeMovieStream << "./make_movie_Q \"$@\"\n";
  makeMovieStream << "./make_movie_Qeq \"$@\"\n";
  makeMovieStream << "./make_movie_P \"$@\"\n";
  makeMovieStream << "./make_movie_Pnorm \"$@\"\n";
  makeMovieStream << "./make_movie_Dpdt \"$@\"\n";
  makeMovieStream << "./make_movie_Dqdt \"$@\"\n";

  std::filesystem::path path{std::format("Breakthrough/System_{}/make_movies", system.systemId)};
  std::filesystem::permissions(path, std::filesystem::perms::owner_exec, std::filesystem::perm_options::add);
#endif

  createMovieScriptColumnV();
  createMovieScriptColumnPt();
  createMovieScriptColumnQ();
  createMovieScriptColumnQeq();
  createMovieScriptColumnP();
  createMovieScriptColumnDpdt();
  createMovieScriptColumnDqdt();
  createMovieScriptColumnPnormalized();
}

// -crf 18: the range of the CRF scale is 0–51, where 0 is lossless, 23 is the default,
//          and 51 is worst quality possible; 18 is visually lossless or nearly so.
// -pix_fmt yuv420p: needed on apple devices
std::string movieScriptTemplate(std::string s)
{
  std::ostringstream stream;

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
  stream << "del column_movie_" << s << ".mp4\n";
  stream << "set /A argVec[1]=1\n";
  stream << "set /A argVec[2]=1200\n";
  stream << "set /A argVec[3]=800\n";
  stream << "set /A argVec[4]=18\n";
  stream << "setlocal enabledelayedexpansion\n";
  stream << "set argCount=0\n";
  stream << "for %%x in (%*) do (\n";
  stream << "   set /A argCount+=1\n";
  stream << "   set \"argVec[!argCount!]=%%~x\"'n";
  stream << ")\n";
  stream << "set PATH=%PATH%;C:\\Program Files\\gnuplot\\bin;"
         << "C:\\Program Files\\ffmpeg-master-latest-win64-gpl\\bin;C:\\Program Files\\ffmpeg\\bin\n";
  stream << "gnuplot.exe -c plot_column_" << s
         << " %argVec[1]% %argVec[2]% %argVec[3]% | ffmpeg.exe "
            "-f png_pipe -s:v \"%argVec[2]%,%argVec[3]%\" -i pipe: -c:v libx264 -pix_fmt yuv420p "
            "-crf %argVec[4]% -c:a aac column_movie_"
         << s + ".mp4\n";
#else
  stream << "rm -f " << "column_movie_" << s << ".mp4\n";
  stream << "every=1\n";
  stream << "format=\"-c:v libx265 -tag:v hvc1\"\n";
  stream << "width=1200\n";
  stream << "height=800\n";
  stream << "quality=18\n";
  stream << "while getopts e:w:h:q:l flag\n";
  stream << "do\n";
  stream << "    case \"${flag}\" in\n";
  stream << "        e) every=${OPTARG};;\n";
  stream << "        w) width=${OPTARG};;\n";
  stream << "        h) height=${OPTARG};;\n";
  stream << "        q) quality=${OPTARG};;\n";
  stream << "        l) format=\"-c:v libx264\";;\n";
  stream << "    esac\n";
  stream << "done\n";
  stream << "gnuplot -c plot_column_" << s
         << " $every $width $height | ffmpeg -f png_pipe "
            "-s:v \"${width},${height}\" -i pipe: $format -pix_fmt yuv420p -crf $quality "
            "-c:a aac column_movie_"
         << s + ".mp4\n";
#endif
  return stream.str();
}

void Breakthrough::createMovieScriptColumnV()
{
  std::filesystem::create_directory("Breakthrough");
  std::filesystem::create_directory(std::format("Breakthrough/System_{}", system.systemId));

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
  std::ofstream makeMovieStream(std::format("Breakthrough/System_{}/make_movie_V.bat", system.systemId));

  std::filesystem::path path{std::format("Breakthrough/System_{}/make_movie_V", system.systemId)};
  std::filesystem::permissions(path, std::filesystem::perms::owner_exec, std::filesystem::perm_options::add);
#else
  std::ofstream makeMovieStream(std::format("Breakthrough/System_{}/make_movie_V", system.systemId));
  makeMovieStream << "#!/bin/sh\n";
  makeMovieStream << "cd -- \"$(dirname \"$0\")\"\n";

  std::filesystem::path path{std::format("Breakthrough/System_{}/make_movie_V", system.systemId)};
  std::filesystem::permissions(path, std::filesystem::perms::owner_exec, std::filesystem::perm_options::add);
#endif
  makeMovieStream << movieScriptTemplate("V");

  std::ofstream stream(std::format("Breakthrough/System_{}/plot_column_V", system.systemId));

  stream << "set encoding utf8\n";
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
  stream << "set terminal pngcairo size ARG2,ARG3 enhanced font 'Arial,10'\n";
  stream << "set xlabel 'Adsorber position / [m]' font 'Arial,14'\n";
  stream << "set ylabel 'Interstitial velocity, {/Arial-Italic v} / [m/s]' offset 0.0,0 font 'Arial,14'\n";
  stream << "set key outside top center horizontal samplen 2.5 height 0.5 spacing 1.5 font 'Arial, 10'\n";
#else
  stream << "set terminal pngcairo size ARG2,ARG3 enhanced font 'Helvetica,10'\n";
  stream << "set xlabel 'Adsorber position / [m]' font 'Helvetica,18'\n";
  stream << "set ylabel 'Interstitial velocity, {/Helvetica-Italic v} / [m/s]' offset 0.0,0 font 'Helvetica,18'\n";
  stream << "set key outside top center horizontal samplen 2.5 height 0.5 spacing 1.5 font 'Helvetica, 10'\n";
#endif

  // colorscheme from book 'gnuplot in action', listing 12.7
  stream << "set linetype 1 pt 5 ps 1 lw 4 lc rgb '0xee0000'\n";
  stream << "set linetype 2 pt 7 ps 1 lw 4 lc rgb '0x008b00'\n";
  stream << "set linetype 3 pt 9 ps 1 lw 4 lc rgb '0x0000cd'\n";
  stream << "set linetype 4 pt 11 ps 1 lw 4 lc rgb '0xff3fb3'\n";
  stream << "set linetype 5 pt 13 ps 1 lw 4 lc rgb '0x00cdcd'\n";
  stream << "set linetype 6 pt 15 ps 1 lw 4 lc rgb '0xcd9b1d'\n";
  stream << "set linetype 7 pt  4 ps 1 lw 4 lc rgb '0x8968ed'\n";
  stream << "set linetype 8 pt  6 ps 1 lw 4 lc rgb '0x8b8b83'\n";
  stream << "set linetype 9 pt  8 ps 1 lw 4 lc rgb '0x00bb00'\n";
  stream << "set linetype 10 pt 10 ps 1 lw 4 lc rgb '0x1e90ff'\n";
  stream << "set linetype 11 pt 12 ps 1 lw 4 lc rgb '0x8b2500'\n";
  stream << "set linetype 12 pt 14 ps 1 lw 4 lc rgb '0x000000'\n";

  stream << "set bmargin 4\n";
  stream << "set title '" << displayName << " {/:Italic T}=" << T << " K, {/:Italic p_t}=" << p_total * 1e-3
         << " kPa'\n";
  stream << "stats 'column.txt' us 2 nooutput\n";
  stream << "max=STATS_max\n";
  stream << "stats 'column.txt' us 1 nooutput\n";
  stream << "set xrange[0:STATS_max]\n";
  stream << "set yrange[0:1.1*max]\n";
  stream << "ev=int(ARG1)\n";
  stream << "do for [i=0:int((STATS_blocks-2)/ev)] {\n";
  stream << "  plot \\\n";
  stream << "    " << "'column.txt'" << " us 1:2 index ev*i notitle with li lt 1,\\\n";
  stream << "    " << "'column.txt'" << " us 1:2 index ev*i notitle with po lt 1\n";
  stream << "}\n";
}

void Breakthrough::createMovieScriptColumnPt()
{
  std::filesystem::create_directory("Breakthrough");
  std::filesystem::create_directory(std::format("Breakthrough/System_{}", system.systemId));

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
  std::ofstream makeMovieStream(std::format("Breakthrough/System_{}/make_movie_Pt.bat", system.systemId));

  std::filesystem::path path{std::format("Breakthrough/System_{}/make_movie_Pt.bat", system.systemId)};
  std::filesystem::permissions(path, std::filesystem::perms::owner_exec, std::filesystem::perm_options::add);
#else
  std::ofstream makeMovieStream(std::format("Breakthrough/System_{}/make_movie_Pt.bat", system.systemId));
  makeMovieStream << "#!/bin/sh\n";
  makeMovieStream << "cd -- \"$(dirname \"$0\")\"\n";

  std::filesystem::path path{std::format("Breakthrough/System_{}/make_movie_Pt.bat", system.systemId)};
  std::filesystem::permissions(path, std::filesystem::perms::owner_exec, std::filesystem::perm_options::add);
#endif
  makeMovieStream << movieScriptTemplate("Pt");

  std::ofstream stream(std::format("Breakthrough/System_{}/plot_column_Pt", system.systemId));

  stream << "set encoding utf8\n";
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
  stream << "set terminal pngcairo size ARG2,ARG3 enhanced font 'Arial,10'\n";
  stream << "set xlabel 'Adsorber position / [m]' font 'Arial,14'\n";
  stream << "set ylabel 'Total Pressure, {/Arial-Italic p_t} / [Pa]' offset 0.0,0 font 'Arial,14'\n";
  stream << "set key outside top center horizontal samplen 2.5 height 0.5 spacing 1.5 font 'Arial, 10'\n";
#else
  stream << "set terminal pngcairo size ARG2,ARG3 enhanced font 'Helvetica,10'\n";
  stream << "set xlabel 'Adsorber position / [m]' font 'Helvetica,18'\n";
  stream << "set ylabel 'Total Pressure, {/Helvetica-Italic p_t} / [Pa]' offset 0.0,0 font 'Helvetica,18'\n";
  stream << "set key outside top center horizontal samplen 2.5 height 0.5 spacing 1.5 font 'Helvetica, 10'\n";
#endif

  // colorscheme from book 'gnuplot in action', listing 12.7
  stream << "set linetype 1 pt 5 ps 1 lw 4 lc rgb '0xee0000'\n";
  stream << "set linetype 2 pt 7 ps 1 lw 4 lc rgb '0x008b00'\n";
  stream << "set linetype 3 pt 9 ps 1 lw 4 lc rgb '0x0000cd'\n";
  stream << "set linetype 4 pt 11 ps 1 lw 4 lc rgb '0xff3fb3'\n";
  stream << "set linetype 5 pt 13 ps 1 lw 4 lc rgb '0x00cdcd'\n";
  stream << "set linetype 6 pt 15 ps 1 lw 4 lc rgb '0xcd9b1d'\n";
  stream << "set linetype 7 pt  4 ps 1 lw 4 lc rgb '0x8968ed'\n";
  stream << "set linetype 8 pt  6 ps 1 lw 4 lc rgb '0x8b8b83'\n";
  stream << "set linetype 9 pt  8 ps 1 lw 4 lc rgb '0x00bb00'\n";
  stream << "set linetype 10 pt 10 ps 1 lw 4 lc rgb '0x1e90ff'\n";
  stream << "set linetype 11 pt 12 ps 1 lw 4 lc rgb '0x8b2500'\n";
  stream << "set linetype 12 pt 14 ps 1 lw 4 lc rgb '0x000000'\n";

  stream << "set bmargin 4\n";
  stream << "set title '" << displayName << " {/:Italic T}=" << T << " K, {/:Italic p_t}=" << p_total * 1e-3
         << " kPa'\n";
  stream << "stats 'column.txt' us 3 nooutput\n";
  stream << "max=STATS_max\n";
  stream << "stats 'column.txt' us 1 nooutput\n";
  stream << "set xrange[0:STATS_max]\n";
  stream << "set yrange[0:1.1*max]\n";
  stream << "ev=int(ARG1)\n";
  stream << "do for [i=0:int((STATS_blocks-2)/ev)] {\n";
  stream << "  plot \\\n";
  stream << "    " << "'column.txt'" << " us 1:3 index ev*i notitle with li lt 1,\\\n";
  stream << "    " << "'column.txt'" << " us 1:3 index ev*i notitle with po lt 1\n";
  stream << "}\n";
}

void Breakthrough::createMovieScriptColumnQ()
{
  std::filesystem::create_directory("Breakthrough");
  std::filesystem::create_directory(std::format("Breakthrough/System_{}", system.systemId));

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
  std::ofstream makeMovieStream(std::format("Breakthrough/System_{}/make_movie_Q.bat", system.systemId));

  std::filesystem::path path{std::format("Breakthrough/System_{}/make_movie_Q.bat", system.systemId)};
  std::filesystem::permissions(path, std::filesystem::perms::owner_exec, std::filesystem::perm_options::add);
#else
  std::ofstream makeMovieStream(std::format("Breakthrough/System_{}/make_movie_Q", system.systemId));
  makeMovieStream << "#!/bin/sh\n";
  makeMovieStream << "cd -- \"$(dirname \"$0\")\"\n";

  std::filesystem::path path{std::format("Breakthrough/System_{}/make_movie_Q", system.systemId)};
  std::filesystem::permissions(path, std::filesystem::perms::owner_exec, std::filesystem::perm_options::add);
#endif
  makeMovieStream << movieScriptTemplate("Q");

  std::ofstream stream(std::format("Breakthrough/System_{}/plot_column_Q", system.systemId));

  stream << "set encoding utf8\n";
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
  stream << "set terminal pngcairo size ARG2,ARG3 enhanced font 'Arial,10'\n";
  stream << "set xlabel 'Adsorber position / [m]' font 'Arial,14'\n";
  stream << "set ylabel 'Concentration, {/Arial-Italic c}_i / [mol/kg]' offset 0.0,0 font 'Arial,14'\n";
  stream << "set key outside top center horizontal samplen 2.5 height 0.5 spacing 1.5 font 'Arial, 10'\n";
#else
  stream << "set terminal pngcairo size ARG2,ARG3 enhanced font 'Helvetica,10'\n";
  stream << "set xlabel 'Adsorber position / [m]' font 'Helvetica,18'\n";
  stream << "set ylabel 'Concentration, {/Helvetica-Italic c}_i / [mol/kg]' offset 0.0,0 font 'Helvetica,18'\n";
  stream << "set key outside top center horizontal samplen 2.5 height 0.5 spacing 1.5 font 'Helvetica, 10'\n";
#endif

  // colorscheme from book 'gnuplot in action', listing 12.7
  stream << "set linetype 1 pt 5 ps 1 lw 4 lc rgb '0xee0000'\n";
  stream << "set linetype 2 pt 7 ps 1 lw 4 lc rgb '0x008b00'\n";
  stream << "set linetype 3 pt 9 ps 1 lw 4 lc rgb '0x0000cd'\n";
  stream << "set linetype 4 pt 11 ps 1 lw 4 lc rgb '0xff3fb3'\n";
  stream << "set linetype 5 pt 13 ps 1 lw 4 lc rgb '0x00cdcd'\n";
  stream << "set linetype 6 pt 15 ps 1 lw 4 lc rgb '0xcd9b1d'\n";
  stream << "set linetype 7 pt  4 ps 1 lw 4 lc rgb '0x8968ed'\n";
  stream << "set linetype 8 pt  6 ps 1 lw 4 lc rgb '0x8b8b83'\n";
  stream << "set linetype 9 pt  8 ps 1 lw 4 lc rgb '0x00bb00'\n";
  stream << "set linetype 10 pt 10 ps 1 lw 4 lc rgb '0x1e90ff'\n";
  stream << "set linetype 11 pt 12 ps 1 lw 4 lc rgb '0x8b2500'\n";
  stream << "set linetype 12 pt 14 ps 1 lw 4 lc rgb '0x000000'\n";

  stream << "set bmargin 4\n";
  stream << "set key title '" << displayName << " {/:Italic T}=" << T << " K, {/:Italic p_t}=" << p_total * 1e-3
         << " kPa'\n";
  stream << "stats 'column.txt' nooutput\n";
  stream << "max = 0.0;\n";
  stream << "do for [i=4:STATS_columns:6] {\n";
  stream << "  stats 'column.txt' us i nooutput\n";
  stream << "  if (max<STATS_max) {\n";
  stream << "    max=STATS_max\n";
  stream << "  }\n";
  stream << "}\n";
  stream << "stats 'column.txt' us 1 nooutput\n";
  stream << "set xrange[0:STATS_max]\n";
  stream << "set yrange[0:1.1*max]\n";
  stream << "ev=int(ARG1)\n";
  stream << "do for [i=0:int((STATS_blocks-2)/ev)] {\n";
  stream << "  plot \\\n";
  for (size_t i = 0; i < Ncomp; i++)
  {
    stream << "    " << "'column.txt'" << " us 1:" << std::to_string(4 + i * 6) << " index ev*i notitle "
           << " with li lt " << i + 1 << ",\\\n";
  }
  for (size_t i = 0; i < Ncomp; i++)
  {
    stream << "    " << "'column.txt'" << " us 1:" << std::to_string(4 + i * 6) << " index ev*i title '"
           << components[i].name << " (y_i=" << components[i].molFraction << ")'" << " with po lt " << i + 1
           << (i < Ncomp - 1 ? ",\\" : "") << "\n";
  }
  stream << "}\n";
}

void Breakthrough::createMovieScriptColumnQeq()
{
  std::filesystem::create_directory("Breakthrough");
  std::filesystem::create_directory(std::format("Breakthrough/System_{}", system.systemId));

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
  std::ofstream makeMovieStream(std::format("Breakthrough/System_{}/make_movie_Qeq.bat", system.systemId));

  std::filesystem::path path{std::format("Breakthrough/System_{}/make_movie_Qeq.bat", system.systemId)};
  std::filesystem::permissions(path, std::filesystem::perms::owner_exec, std::filesystem::perm_options::add);
#else
  std::ofstream makeMovieStream(std::format("Breakthrough/System_{}/make_movie_Qeq", system.systemId));
  makeMovieStream << "#!/bin/sh\n";
  makeMovieStream << "cd -- \"$(dirname \"$0\")\"\n";

  std::filesystem::path path{std::format("Breakthrough/System_{}/make_movie_Qeq", system.systemId)};
  std::filesystem::permissions(path, std::filesystem::perms::owner_exec, std::filesystem::perm_options::add);
#endif
  makeMovieStream << movieScriptTemplate("Qeq");

  std::ofstream stream(std::format("Breakthrough/System_{}/plot_column_Qeq", system.systemId));

  stream << "set encoding utf8\n";
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
  stream << "set terminal pngcairo size ARG2,ARG3 enhanced font 'Arial,10'\n";
  stream << "set xlabel 'Adsorber position / [m]' font 'Arial,14'\n";
  stream << "set ylabel 'Concentration, {/Arial-Italic c}_i / [mol/kg]' offset 0.0,0 font 'Arial,14'\n";
  stream << "set key outside top center horizontal samplen 2.5 height 0.5 spacing 1.5 font 'Arial, 10'\n";
#else
  stream << "set terminal pngcairo size ARG2,ARG3 enhanced font 'Helvetica,10'\n";
  stream << "set xlabel 'Adsorber position / [m]' font 'Helvetica,18'\n";
  stream << "set ylabel 'Concentration, {/Helvetica-Italic c}_i / [mol/kg]' offset 0.0,0 font 'Helvetica,18'\n";
  stream << "set key outside top center horizontal samplen 2.5 height 0.5 spacing 1.5 font 'Helvetica, 10'\n";
#endif

  // colorscheme from book 'gnuplot in action', listing 12.7
  stream << "set linetype 1 pt 5 ps 1 lw 4 lc rgb '0xee0000'\n";
  stream << "set linetype 2 pt 7 ps 1 lw 4 lc rgb '0x008b00'\n";
  stream << "set linetype 3 pt 9 ps 1 lw 4 lc rgb '0x0000cd'\n";
  stream << "set linetype 4 pt 11 ps 1 lw 4 lc rgb '0xff3fb3'\n";
  stream << "set linetype 5 pt 13 ps 1 lw 4 lc rgb '0x00cdcd'\n";
  stream << "set linetype 6 pt 15 ps 1 lw 4 lc rgb '0xcd9b1d'\n";
  stream << "set linetype 7 pt  4 ps 1 lw 4 lc rgb '0x8968ed'\n";
  stream << "set linetype 8 pt  6 ps 1 lw 4 lc rgb '0x8b8b83'\n";
  stream << "set linetype 9 pt  8 ps 1 lw 4 lc rgb '0x00bb00'\n";
  stream << "set linetype 10 pt 10 ps 1 lw 4 lc rgb '0x1e90ff'\n";
  stream << "set linetype 11 pt 12 ps 1 lw 4 lc rgb '0x8b2500'\n";
  stream << "set linetype 12 pt 14 ps 1 lw 4 lc rgb '0x000000'\n";

  stream << "set bmargin 4\n";
  stream << "set key title '" << displayName << " {/:Italic T}=" << T << " K, {/:Italic p_t}=" << p_total * 1e-3
         << " kPa'\n";
  stream << "stats 'column.txt' nooutput\n";
  stream << "max = 0.0;\n";
  stream << "do for [i=5:STATS_columns:6] {\n";
  stream << "  stats 'column.txt' us i nooutput\n";
  stream << "  if (max<STATS_max) {\n";
  stream << "    max=STATS_max\n";
  stream << "  }\n";
  stream << "}\n";
  stream << "stats 'column.txt' us 1 nooutput\n";
  stream << "set xrange[0:STATS_max]\n";
  stream << "set yrange[0:1.1*max]\n";
  stream << "ev=int(ARG1)\n";
  stream << "do for [i=0:int((STATS_blocks-2)/ev)] {\n";
  stream << "  plot \\\n";
  for (size_t i = 0; i < Ncomp; i++)
  {
    stream << "    " << "'column.txt'" << " us 1:" << std::to_string(5 + i * 6) << " index ev*i notitle "
           << " with li lt " << i + 1 << ",\\\n";
  }
  for (size_t i = 0; i < Ncomp; i++)
  {
    stream << "    " << "'column.txt'" << " us 1:" << std::to_string(5 + i * 6) << " index ev*i title '"
           << components[i].name << " (y_i=" << components[i].molFraction << ")'" << " with po lt " << i + 1
           << (i < Ncomp - 1 ? ",\\" : "") << "\n";
  }
  stream << "}\n";
}

void Breakthrough::createMovieScriptColumnP()
{
  std::filesystem::create_directory("Breakthrough");
  std::filesystem::create_directory(std::format("Breakthrough/System_{}", system.systemId));

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
  std::ofstream makeMovieStream(std::format("Breakthrough/System_{}/make_movie_P.bat", system.systemId));

  std::filesystem::path path{std::format("Breakthrough/System_{}/make_movie_P.bat", system.systemId)};
  std::filesystem::permissions(path, std::filesystem::perms::owner_exec, std::filesystem::perm_options::add);
#else
  std::ofstream makeMovieStream(std::format("Breakthrough/System_{}/make_movie_P", system.systemId));
  makeMovieStream << "#!/bin/sh\n";
  makeMovieStream << "cd -- \"$(dirname \"$0\")\"\n";

  std::filesystem::path path{std::format("Breakthrough/System_{}/make_movie_P", system.systemId)};
  std::filesystem::permissions(path, std::filesystem::perms::owner_exec, std::filesystem::perm_options::add);
#endif
  makeMovieStream << movieScriptTemplate("P");

  std::ofstream stream(std::format("Breakthrough/System_{}/plot_column_P", system.systemId));

  stream << "set encoding utf8\n";
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
  stream << "set terminal pngcairo size ARG2,ARG3 enhanced font 'Arial,10'\n";
  stream << "set xlabel 'Adsorber position / [m]' font 'Arial,14'\n";
  stream << "set ylabel 'Partial pressure, {/Arial-Italic p}_i / [Pa]' offset 0.0,0 font 'Arial,14'\n";
  stream << "set key outside top center horizontal samplen 2.5 height 0.5 spacing 1.5 font 'Arial, 10'\n";
#else
  stream << "set terminal pngcairo size ARG2,ARG3 enhanced font 'Helvetica,10'\n";
  stream << "set xlabel 'Adsorber position / [m]' font 'Helvetica,18'\n";
  stream << "set ylabel 'Partial pressure, {/Helvetica-Italic p}_i / [Pa]' offset 0.0,0 font 'Helvetica,18'\n";
  stream << "set key outside top center horizontal samplen 2.5 height 0.5 spacing 1.5 font 'Helvetica, 10'\n";
#endif

  // colorscheme from book 'gnuplot in action', listing 12.7
  stream << "set linetype 1 pt 5 ps 1 lw 4 lc rgb '0xee0000'\n";
  stream << "set linetype 2 pt 7 ps 1 lw 4 lc rgb '0x008b00'\n";
  stream << "set linetype 3 pt 9 ps 1 lw 4 lc rgb '0x0000cd'\n";
  stream << "set linetype 4 pt 11 ps 1 lw 4 lc rgb '0xff3fb3'\n";
  stream << "set linetype 5 pt 13 ps 1 lw 4 lc rgb '0x00cdcd'\n";
  stream << "set linetype 6 pt 15 ps 1 lw 4 lc rgb '0xcd9b1d'\n";
  stream << "set linetype 7 pt  4 ps 1 lw 4 lc rgb '0x8968ed'\n";
  stream << "set linetype 8 pt  6 ps 1 lw 4 lc rgb '0x8b8b83'\n";
  stream << "set linetype 9 pt  8 ps 1 lw 4 lc rgb '0x00bb00'\n";
  stream << "set linetype 10 pt 10 ps 1 lw 4 lc rgb '0x1e90ff'\n";
  stream << "set linetype 11 pt 12 ps 1 lw 4 lc rgb '0x8b2500'\n";
  stream << "set linetype 12 pt 14 ps 1 lw 4 lc rgb '0x000000'\n";

  stream << "set bmargin 4\n";
  stream << "set key title '" << displayName << " {/:Italic T}=" << T << " K, {/:Italic p_t}=" << p_total * 1e-3
         << " kPa'\n";
  stream << "stats 'column.txt' nooutput\n";
  stream << "max = 0.0;\n";
  stream << "do for [i=6:STATS_columns:6] {\n";
  stream << "  stats 'column.txt' us i nooutput\n";
  stream << "  if (max<STATS_max) {\n";
  stream << "    max=STATS_max\n";
  stream << "  }\n";
  stream << "}\n";
  stream << "stats 'column.txt' us 1 nooutput\n";
  stream << "set xrange[0:STATS_max]\n";
  stream << "set yrange[0:1.1*max]\n";
  stream << "ev=int(ARG1)\n";
  stream << "do for [i=0:int((STATS_blocks-2)/ev)] {\n";
  stream << "  plot \\\n";
  for (size_t i = 0; i < Ncomp; i++)
  {
    stream << "    " << "'column.txt'" << " us 1:" << std::to_string(6 + i * 6) << " index ev*i notitle "
           << " with li lt " << i + 1 << ",\\\n";
  }
  for (size_t i = 0; i < Ncomp; i++)
  {
    stream << "    " << "'column.txt'" << " us 1:" << std::to_string(6 + i * 6) << " index ev*i title '"
           << components[i].name << " (y_i=" << components[i].molFraction << ")'" << " with po lt " << i + 1
           << (i < Ncomp - 1 ? ",\\" : "") << "\n";
  }
  stream << "}\n";
}

void Breakthrough::createMovieScriptColumnPnormalized()
{
  std::filesystem::create_directory("Breakthrough");
  std::filesystem::create_directory(std::format("Breakthrough/System_{}", system.systemId));

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
  std::ofstream makeMovieStream(std::format("Breakthrough/System_{}/make_movie_Pnorm.bat", system.systemId));

  std::filesystem::path path{std::format("Breakthrough/System_{}/make_movie_Pnorm.bat", system.systemId)};
  std::filesystem::permissions(path, std::filesystem::perms::owner_exec, std::filesystem::perm_options::add);
#else
  std::ofstream makeMovieStream(std::format("Breakthrough/System_{}/make_movie_Pnorm", system.systemId));
  makeMovieStream << "#!/bin/sh\n";
  makeMovieStream << "cd -- \"$(dirname \"$0\")\"\n";

  std::filesystem::path path{std::format("Breakthrough/System_{}/make_movie_Pnorm", system.systemId)};
  std::filesystem::permissions(path, std::filesystem::perms::owner_exec, std::filesystem::perm_options::add);
#endif
  makeMovieStream << movieScriptTemplate("Pnorm");

  std::ofstream stream(std::format("Breakthrough/System_{}/plot_column_Pnorm", system.systemId));

  stream << "set encoding utf8\n";
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
  stream << "set terminal pngcairo size ARG2,ARG3 enhanced font 'Arial,10'\n";
  stream << "set xlabel 'Adsorber position / [m]' font 'Arial,14'\n";
  stream << "set ylabel 'Partial pressure, {/Arial-Italic p}_i / [-]' offset 0.0,0 font 'Arial,14'\n";
  stream << "set key outside top center horizontal samplen 2.5 height 0.5 spacing 1.5 font 'Arial, 10'\n";
#else
  stream << "set terminal pngcairo size ARG2,ARG3 enhanced font 'Helvetica,10'\n";
  stream << "set xlabel 'Adsorber position / [m]' font 'Helvetica,18'\n";
  stream << "set ylabel 'Partial pressure, {/Helvetica-Italic p}_i / [-]' offset 0.0,0 font 'Helvetica,18'\n";
  stream << "set key outside top center horizontal samplen 2.5 height 0.5 spacing 1.5 font 'Helvetica, 10'\n";
#endif

  // colorscheme from book 'gnuplot in action', listing 12.7
  stream << "set linetype 1 pt 5 ps 1 lw 4 lc rgb '0xee0000'\n";
  stream << "set linetype 2 pt 7 ps 1 lw 4 lc rgb '0x008b00'\n";
  stream << "set linetype 3 pt 9 ps 1 lw 4 lc rgb '0x0000cd'\n";
  stream << "set linetype 4 pt 11 ps 1 lw 4 lc rgb '0xff3fb3'\n";
  stream << "set linetype 5 pt 13 ps 1 lw 4 lc rgb '0x00cdcd'\n";
  stream << "set linetype 6 pt 15 ps 1 lw 4 lc rgb '0xcd9b1d'\n";
  stream << "set linetype 7 pt  4 ps 1 lw 4 lc rgb '0x8968ed'\n";
  stream << "set linetype 8 pt  6 ps 1 lw 4 lc rgb '0x8b8b83'\n";
  stream << "set linetype 9 pt  8 ps 1 lw 4 lc rgb '0x00bb00'\n";
  stream << "set linetype 10 pt 10 ps 1 lw 4 lc rgb '0x1e90ff'\n";
  stream << "set linetype 11 pt 12 ps 1 lw 4 lc rgb '0x8b2500'\n";
  stream << "set linetype 12 pt 14 ps 1 lw 4 lc rgb '0x000000'\n";

  stream << "set bmargin 4\n";
  stream << "set key title '" << displayName << " {/:Italic T}=" << T << " K, {/:Italic p_t}=" << p_total * 1e-3
         << " kPa'\n";
  stream << "stats 'column.txt' nooutput\n";
  stream << "max = 0.0;\n";
  stream << "do for [i=7:STATS_columns:6] {\n";
  stream << "  stats 'column.txt' us i nooutput\n";
  stream << "  if (max<STATS_max) {\n";
  stream << "    max=STATS_max\n";
  stream << "  }\n";
  stream << "}\n";
  stream << "stats 'column.txt' us 1 nooutput\n";
  stream << "set xrange[0:STATS_max]\n";
  stream << "set yrange[0:1.1*max]\n";
  stream << "ev=int(ARG1)\n";
  stream << "do for [i=0:int((STATS_blocks-2)/ev)] {\n";
  stream << "  plot \\\n";
  for (size_t i = 0; i < Ncomp; i++)
  {
    stream << "    " << "'column.txt'" << " us 1:" << std::to_string(7 + i * 6) << " index ev*i notitle "
           << " with li lt " << i + 1 << ",\\\n";
  }
  for (size_t i = 0; i < Ncomp; i++)
  {
    stream << "    " << "'column.txt'" << " us 1:" << std::to_string(7 + i * 6) << " index ev*i title '"
           << components[i].name << " (y_i=" << components[i].molFraction << ")'" << " with po lt " << i + 1
           << (i < Ncomp - 1 ? ",\\" : "") << "\n";
  }
  stream << "}\n";
}

void Breakthrough::createMovieScriptColumnDpdt()
{
  std::filesystem::create_directory("Breakthrough");
  std::filesystem::create_directory(std::format("Breakthrough/System_{}", system.systemId));

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
  std::ofstream makeMovieStream(std::format("Breakthrough/System_{}/make_movie_Dpdt.bat", system.systemId));

  std::filesystem::path path{std::format("Breakthrough/System_{}/make_movie_Dpdt.bat", system.systemId)};
  std::filesystem::permissions(path, std::filesystem::perms::owner_exec, std::filesystem::perm_options::add);
#else
  std::ofstream makeMovieStream(std::format("Breakthrough/System_{}/make_movie_Dpdt", system.systemId));
  makeMovieStream << "#!/bin/sh\n";
  makeMovieStream << "cd -- \"$(dirname \"$0\")\"\n";

  std::filesystem::path path{std::format("Breakthrough/System_{}/make_movie_Dpdt", system.systemId)};
  std::filesystem::permissions(path, std::filesystem::perms::owner_exec, std::filesystem::perm_options::add);
#endif
  makeMovieStream << movieScriptTemplate("Dpdt");

  std::ofstream stream(std::format("Breakthrough/System_{}/plot_column_Dpdt", system.systemId));

  stream << "set encoding utf8\n";
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
  stream << "set terminal pngcairo size ARG2,ARG3 enhanced font 'Arial,10'\n";
  stream << "set xlabel 'Adsorber position / [m]' font 'Arial,14'\n";
  stream << "set ylabel 'Pressure derivative, {/Arial-Italic dp_/dt} / [Pa/s]' offset 0.0,0 font 'Arial,14'\n";
  stream << "set key outside top center horizontal samplen 2.5 height 0.5 spacing 1.5 font 'Arial, 10'\n";
#else
  stream << "set terminal pngcairo size ARG2,ARG3 enhanced font 'Helvetica,10'\n";
  stream << "set xlabel 'Adsorber position / [m]' font 'Helvetica,18'\n";
  stream << "set ylabel 'Pressure derivative, {/Helvetica-Italic dp_/dt} / [Pa/s]' "
            "offset 0.0,0 font 'Helvetica,18'\n";
  stream << "set key outside top center horizontal samplen 2.5 height 0.5 spacing 1.5 font 'Helvetica, 10'\n";
#endif

  // colorscheme from book 'gnuplot in action', listing 12.7
  stream << "set linetype 1 pt 5 ps 1 lw 4 lc rgb '0xee0000'\n";
  stream << "set linetype 2 pt 7 ps 1 lw 4 lc rgb '0x008b00'\n";
  stream << "set linetype 3 pt 9 ps 1 lw 4 lc rgb '0x0000cd'\n";
  stream << "set linetype 4 pt 11 ps 1 lw 4 lc rgb '0xff3fb3'\n";
  stream << "set linetype 5 pt 13 ps 1 lw 4 lc rgb '0x00cdcd'\n";
  stream << "set linetype 6 pt 15 ps 1 lw 4 lc rgb '0xcd9b1d'\n";
  stream << "set linetype 7 pt  4 ps 1 lw 4 lc rgb '0x8968ed'\n";
  stream << "set linetype 8 pt  6 ps 1 lw 4 lc rgb '0x8b8b83'\n";
  stream << "set linetype 9 pt  8 ps 1 lw 4 lc rgb '0x00bb00'\n";
  stream << "set linetype 10 pt 10 ps 1 lw 4 lc rgb '0x1e90ff'\n";
  stream << "set linetype 11 pt 12 ps 1 lw 4 lc rgb '0x8b2500'\n";
  stream << "set linetype 12 pt 14 ps 1 lw 4 lc rgb '0x000000'\n";

  stream << "set bmargin 4\n";
  stream << "set key title '" << displayName << " {/:Italic T}=" << T << " K, {/:Italic p_t}=" << p_total * 1e-3
         << " kPa'\n";
  stream << "stats 'column.txt' nooutput\n";
  stream << "max = -1e10;\n";
  stream << "min = 1e10;\n";
  stream << "do for [i=8:STATS_columns:6] {\n";
  stream << "  stats 'column.txt' us i nooutput\n";
  stream << "  if (STATS_max>max) {\n";
  stream << "    max=STATS_max\n";
  stream << "  }\n";
  stream << "  if (STATS_min<min) {\n";
  stream << "    min=STATS_min\n";
  stream << "  }\n";
  stream << "}\n";
  stream << "stats 'column.txt' us 1 nooutput\n";
  stream << "set xrange[0:STATS_max]\n";
  stream << "set yrange[1.1*min:1.1*max]\n";
  stream << "ev=int(ARG1)\n";
  stream << "do for [i=0:int((STATS_blocks-2)/ev)] {\n";
  stream << "  plot \\\n";
  for (size_t i = 0; i < Ncomp; i++)
  {
    stream << "    " << "'column.txt'" << " us 1:" << std::to_string(8 + i * 6) << " index ev*i notitle "
           << " with li lt " << i + 1 << ",\\\n";
  }
  for (size_t i = 0; i < Ncomp; i++)
  {
    stream << "    " << "'column.txt'" << " us 1:" << std::to_string(8 + i * 6) << " index ev*i title '"
           << components[i].name << " (y_i=" << components[i].molFraction << ")'" << " with po lt " << i + 1
           << (i < Ncomp - 1 ? ",\\" : "") << "\n";
  }
  stream << "}\n";
}

void Breakthrough::createMovieScriptColumnDqdt()
{
  std::filesystem::create_directory("Breakthrough");
  std::filesystem::create_directory(std::format("Breakthrough/System_{}", system.systemId));

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
  std::ofstream makeMovieStream(std::format("Breakthrough/System_{}/make_movie_Dqdt.bat", system.systemId));

  std::filesystem::path path{std::format("Breakthrough/System_{}/make_movie_Dqdt.bat", system.systemId)};
  std::filesystem::permissions(path, std::filesystem::perms::owner_exec, std::filesystem::perm_options::add);
#else
  std::ofstream makeMovieStream(std::format("Breakthrough/System_{}/make_movie_Dqdt", system.systemId));
  makeMovieStream << "#!/bin/sh\n";
  makeMovieStream << "cd -- \"$(dirname \"$0\")\"\n";

  std::filesystem::path path{std::format("Breakthrough/System_{}/make_movie_Dqdt", system.systemId)};
  std::filesystem::permissions(path, std::filesystem::perms::owner_exec, std::filesystem::perm_options::add);
#endif
  makeMovieStream << movieScriptTemplate("Dqdt");

  std::ofstream stream(std::format("Breakthrough/System_{}/plot_column_Dqdt", system.systemId));

  stream << "set encoding utf8\n";
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
  stream << "set terminal pngcairo size ARG2,ARG3 enhanced font 'Arial,10'\n";
  stream << "set xlabel 'Adsorber position / [m]' font 'Arial,14'\n";
  stream << "set ylabel 'Loading derivative, {/Arial-Italic dq_i/dt} / [mol/kg/s]' offset 0.0,0 font 'Arial,14'\n";
  stream << "set key outside top center horizontal samplen 2.5 height 0.5 spacing 1.5 font 'Arial, 10'\n";
#else
  stream << "set terminal pngcairo size ARG2,ARG3 enhanced font 'Helvetica,10'\n";
  stream << "set xlabel 'Adsorber position / [m]' font 'Helvetica,18'\n";
  stream << "set ylabel 'Loading derivative, {/Helvetica-Italic dq_i/dt} / [mol/kg/s]' "
            "offset 0.0,0 font 'Helvetica,18'\n";
  stream << "set key outside top center horizontal samplen 2.5 height 0.5 spacing 1.5 font 'Helvetica, 10'\n";
#endif

  // colorscheme from book 'gnuplot in action', listing 12.7
  stream << "set linetype 1 pt 5 ps 1 lw 4 lc rgb '0xee0000'\n";
  stream << "set linetype 2 pt 7 ps 1 lw 4 lc rgb '0x008b00'\n";
  stream << "set linetype 3 pt 9 ps 1 lw 4 lc rgb '0x0000cd'\n";
  stream << "set linetype 4 pt 11 ps 1 lw 4 lc rgb '0xff3fb3'\n";
  stream << "set linetype 5 pt 13 ps 1 lw 4 lc rgb '0x00cdcd'\n";
  stream << "set linetype 6 pt 15 ps 1 lw 4 lc rgb '0xcd9b1d'\n";
  stream << "set linetype 7 pt  4 ps 1 lw 4 lc rgb '0x8968ed'\n";
  stream << "set linetype 8 pt  6 ps 1 lw 4 lc rgb '0x8b8b83'\n";
  stream << "set linetype 9 pt  8 ps 1 lw 4 lc rgb '0x00bb00'\n";
  stream << "set linetype 10 pt 10 ps 1 lw 4 lc rgb '0x1e90ff'\n";
  stream << "set linetype 11 pt 12 ps 1 lw 4 lc rgb '0x8b2500'\n";
  stream << "set linetype 12 pt 14 ps 1 lw 4 lc rgb '0x000000'\n";

  stream << "set bmargin 4\n";
  stream << "set key title '" << displayName << " {/:Italic T}=" << T << " K, {/:Italic p_t}=" << p_total * 1e-3
         << " kPa'\n";
  stream << "stats 'column.txt' nooutput\n";
  stream << "max = -1e10;\n";
  stream << "min = 1e10;\n";
  stream << "min = 10000000000000.0;\n";
  stream << "do for [i=9:STATS_columns:6] {\n";
  stream << "  stats 'column.txt' us i nooutput\n";
  stream << "  if (STATS_max>max) {\n";
  stream << "    max=STATS_max\n";
  stream << "  }\n";
  stream << "  if (STATS_min<min) {\n";
  stream << "    min=STATS_min\n";
  stream << "  }\n";
  stream << "}\n";
  stream << "stats 'column.txt' us 1 nooutput\n";
  stream << "set xrange[0:STATS_max]\n";
  stream << "set yrange[1.1*min:1.1*max]\n";
  stream << "ev=int(ARG1)\n";
  stream << "do for [i=0:int((STATS_blocks-2)/ev)] {\n";
  stream << "  plot \\\n";
  for (size_t i = 0; i < Ncomp; i++)
  {
    stream << "    " << "'column.txt'" << " us 1:" << std::to_string(9 + i * 6) << " index ev*i notitle "
           << " with li lt " << i + 1 << ",\\\n";
  }
  for (size_t i = 0; i < Ncomp; i++)
  {
    stream << "    " << "'column.txt'" << " us 1:" << std::to_string(9 + i * 6) << " index ev*i title '"
           << components[i].name << " (y_i=" << components[i].molFraction << ")'" << " with po lt " << i + 1
           << (i < Ncomp - 1 ? ",\\" : "") << "\n";
  }
  stream << "}\n";
}
