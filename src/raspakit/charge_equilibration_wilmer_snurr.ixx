module;

#ifdef USE_PRECOMPILED_HEADERS
#include "pch.h"
#endif

#ifdef USE_LEGACY_HEADERS
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <string>
#include <tuple>
#include <vector>
#include <array>
#include <string_view>
#endif

export module charge_equilibration_wilmer_snurr;

#ifdef USE_STD_IMPORT
import std;
#endif

import atom;
import forcefield;
import simulationbox;

using namespace std::literals;

// "Towards rapid computational screening of metal-organic frameworks for carbon
//   dioxide capture: Calculation of framework charges via charge equilibration"
// C.E. Wilmer and R.Q. Snurr
// Chemical Engineering Journal, 171(3), 775-781
//
// "An Extended Charge Equilibration Method",
// Christopher E. Wilmer, Ki Chul Kim, and Randall Q. Snurr,
// The journal of Physical Chemistry Letters, 3, 2506-2511, 2012

export namespace ChargeEquilibration
{
enum class Type : std::size_t
{
  NonPeriodic = 0,
  Periodic = 1,
  PeriodicDirectSum = 2,
  PeriodicEwaldSum = 3
};

void computeChargeEquilibration(const ForceField &forceField, const SimulationBox &box, std::span<Atom> frameworkAtoms,
                                ChargeEquilibration::Type type);

struct ChargeEquilbrationElementData
{
  std::string_view label;
  double electronAffinity;               // in units of [eV]
  std::array<double, 10> ionizationEnergy;  // in units of [eV]
};

//{ "H", 0.7520812,  {13.598443}},
//  {"Zn", 0.0,        {9.394199, 17.96439, 39.723, 59.4, 82.6, 108, 134, 174}},
//  { "O", 1.4611135,  {13.61805, 35.1211, 54.9355, 77.41353, 113.8990, 138.1197, 739.29, 871.4101}},
//  { "C", 1.262119,   {11.26030, 24.3833, 47.8878, 64.4939, 392.087, 489.99334}},

constexpr std::array<ChargeEquilbrationElementData, 1> chargeEquilbrationElementData = {
    {""sv, 0.0, std::array<double, 10>{}}};

  /*
    {""sv, 0.0, {}},
    {"H"sv, -2.4172, {11.4732}},
    {"He"sv, 0.00000, {24.587387, 54.417760}},
    {"Li"sv, 0.618049, {5.391719, 75.6400, 122.45429}},
    {"Be"sv, 0.00000, {9.32270, 18.21114, 153.89661, 217.71865}},
    {"B"sv, 0.279730, {8.29802, 25.1548, 37.93064, 259.37521, 340.22580}},
    {"C"sv, 1.26212, {11.26000, 24.38300, 47.88700, 64.49200, 392.07700, 489.98100}},
    {"N"sv, 0.000000, {14.5341, 29.6013, 47.44924, 77.4735, 97.8902, 552.0718, 667.046}},
    {"O"sv, 1.46111, {13.61800, 35.11600, 54.93400, 77.41200, 113.89600, 138.11600, 739.31500, 871.38700}},
    {"F"sv, 3.4011895, {17.4228, 34.9708, 62.7084, 87.1398, 114.2428, 157.1651, 185.186, 953.9112}},
    {"Ne"sv, 0.000000, {21.56454, 40.96296, 63.45, 97.12, 126.21, 157.93, 207.2759, 239.0989}},
    {"Na"sv, 0.547926, {5.139076, 47.2864, 71.6200, 98.91, 138.40, 172.18, 208.50, 264.25}},
    {"Mg"sv, 0.000000, {7.646235, 15.03527, 80.1437, 109.2655, 141.27, 186.76, 225.02, 265.96}},
    {"Al"sv, 0.43283, {5.985768, 18.82855, 28.44765, 119.992, 153.825, 190.49, 241.76, 284.66}},
    {"Si"sv, 1.3895213, {8.15168, 16.34584, 33.49302, 45.14181, 166.767, 205.27, 246.5, 303.54}},
    {"P"sv, 0.7465, {10.48669, 19.7695, 30.2027, 51.4439, 65.0251, 220.421, 263.57, 309.60}},
    {"S"sv, 2.07710418, {10.36001, 23.33788, 34.79, 47.222, 72.5945, 88.0530, 280.948, 328.75}},
    {"Cl"sv, 3.612724, {12.96763, 23.8136, 39.61, 53.4652, 67.8, 97.03, 114.1958, 348.28}},
    {"Ar"sv, 0.000000, {15.759610, 27.62966, 40.64, 59.81, 75.02, 91.009, 124.323, 143.460}},
    {"K"sv, 0.50147, {4.3406633, 31.63, 45.806, 60.91, 82.66, 99.4, 117.56, 154.88}},
    {"Ca"sv, 0.02455, {6.11316, 11.87172, 50.9131, 67.27, 84.50, 108.78, 127.2, 147.24}},
    {"Sc"sv, 0.188, {6.56149, 12.79977, 24.75666, 73.4894, 91.65, 110.68, 138.0, 158.1}},
    {"Ti"sv, 0.079, {6.82812, 13.5755, 27.4917, 43.2672, 99.30, 119.53, 140.8, 170.4}},
    {"V"sv, 0.525, {6.74619, 14.618, 29.311, 46.709, 65.2817, 128.13, 150.6, 173.4}},
    {"Cr"sv, 0.666, {6.76651, 16.4857, 30.96, 49.16, 69.46, 90.6349, 160.18, 184.7}},
    {"Mn"sv, 0.000000, {7.43402, 15.6400, 33.668, 51.2, 72.4, 95.6, 119.203, 194.5}},
    {"Fe"sv, 0.151, {7.9024, 16.1877, 30.652, 54.8, 75.0, 99.1, 124.98, 151.06}},
    {"Co"sv, 0.662, {7.88101, 17.084, 33.50, 51.3, 79.5, 102.0, 128.9, 157.8}},
    {"Ni"sv, 1.156, {7.6398, 18.16884, 35.19, 54.9, 76.06, 108, 133, 162}},
    {"Cu"sv, 1.235, {7.72638, 20.2924, 36.841, 57.38, 79.8, 103, 139, 166}},
    {"Zn"sv, 0.0, {9.39400, 17.96400, 39.72200, 59.4000, 82.60000, 108.00000, 134.00000, 174.00000}},
    {"Ga"sv, 0.43, {5.999301, 20.51515, 30.7258, 63.241, 86.01, 112.7, 140.9, 169.9}},
    {"Ge"sv, 1.232712, {7.89943, 15.93461, 34.2241, 45.7131, 93.5}},
    {"As"sv, 0.804, {9.7886, 18.5892, 28.351, 50.13, 62.63, 127.6}},
    {"Se"sv, 2.02067, {9.75239, 21.19, 30.8204, 42.9450, 68.3, 81.7, 155.4}},
    {"Br"sv, 3.363588, {11.8138, 21.591, 36, 47.3, 59.7, 88.6, 103.0, 192.8}},
    {"Kr"sv, 0.000000, {13.99961, 24.35984, 36.950, 52.5, 64.7, 78.5, 111.0, 125.802}},
    {"Rb"sv, 0.48592, {4.177128, 27.2895, 40, 52.6, 71.0, 84.4, 99.2, 136}},
    {"Sr"sv, 0.048, {5.69485, 11.0301, 42.89, 57, 71.6, 90.8, 106, 122.3}},
    {"Y"sv, 0.307, {6.2173, 12.224, 20.52, 60.597, 77.0, 93.0, 116, 129}},
    {"Zr"sv, 0.426, {6.63390, 13.1, 22.99, 34.34, 80.348}},
    {"Nb"sv, 0.916, {6.75885, 14.0, 25.04, 38.3, 50.55, 102.057, 125}},
    {"Mo"sv, 0.748, {7.09243, 16.16, 27.13, 46.4, 54.49, 68.8276, 125.664, 143.6}},
    {"Tc"sv, 0.55, {7.28, 15.26, 29.54}},
    {"Ru"sv, 1.05, {7.36050, 16.76, 28.47}},
    {"Rh"sv, 1.137, {7.45890, 18.08, 31.06}},
    {"Pd"sv, 0.562, {8.3369, 19.43, 32.93}},
    {"Ag"sv, 1.302, {7.57623, 21.47746, 34.83}},
    {"Cd"sv, 0.000000, {8.99382, 16.90831, 37.48}},
    {"In"sv, 0.3, {5.78636, 18.8703, 28.03, 54}},
    {"Sn"sv, 1.112067, {7.34392, 14.6322, 30.50260, 40.73502, 72.28}},
    {"Sb"sv, 1.046, {8.60839, 16.63, 25.3, 44.2, 56, 108}},
    {"Te"sv, 1.970876, {9.0096, 18.6, 27.96, 37.41, 58.75, 70.7, 137}},
    {"I"sv, 3.0590463, {10.45126, 19.1313, 33}},
    {"Xe"sv, 0.000000, {12.12984, 20.9750, 32.1230}},
    {"Cs"sv, 0.471626, {3.893905, 23.15744}},
    {"Ba"sv, 0.14462, {5.211664, 10.00383}},
    {"La"sv, 0.47, {5.5769, 11.059, 19.1773, 49.950, 61.6}},
    {"Ce"sv, 0.65, {5.5387, 10.85, 20.198, 36.758, 65.55, 77.6}},
    {"Pr"sv, 0.962, {5.5473, 10.55, 21.624, 38.98, 57.53}},
    {"Nd"sv, 1.916, {5.5250, 10.72, 22.1, 40.4}},
    {"Pm"sv, 0.000000, {5.582, 10.9, 22.3, 41.1}},
    {"Sm"sv, 0.000000, {5.6437, 11.07, 23.4, 41.4}},
    {"Eu"sv, 0.864, {5.67038, 11.25, 24.92, 42.7}},
    {"Gd"sv, 0.000000, {6.14980, 12.09, 20.63, 44.0}},
    {"Tb"sv, 1.165, {5.8638, 11.52, 21.91, 39.79}},
    {"Dy"sv, 0.000000, {5.9389, 11.67, 22.8, 41.47}},
    {"Ho"sv, 0.000000, {6.0215, 11.8, 22.84, 42.5}},
    {"Er"sv, 0.000000, {6.1077, 11.93, 22.74, 42.7}},
    {"Tm"sv, 1.029, {6.18431, 12.05, 23.68, 42.7}},
    {"Yb"sv, -0.020, {6.25416, 12.176, 25.05, 43.56}},
    {"Lu"sv, 0.34, {5.42586, 13.9, 20.9594, 45.25, 66.8}},
    {"Hf"sv, 0.017, {6.82507, 15, 23.3, 33.33}},
    {"Ta"sv, 0.322000, {7.54957}},
    {"W"sv, 0.81626, {7.86403, 16.1}},
    {"Re"sv, 0.15, {7.83352}},
    {"Os"sv, 1.1, {8.43823}},
    {"Ir"sv, 1.5638, {8.96702}},
    {"Pt"sv, 2.128, {8.9588, 18.563}},
    {"Au"sv, 2.30863, {9.22553, 20.2}},
    {"Hg"sv, 0.000000, {10.4375, 18.7568, 34.2}},
    {"Tl"sv, 0.377, {6.108194, 20.4283, 29.83}},
    {"Pb"sv, 0.364, {7.41663, 15.03248, 31.9373, 42.32, 68.8}},
    {"Bi"sv, 0.942362, {7.2855, 16.703, 25.56, 45.3, 56.0, 88.3}},
    {"Po"sv, 1.9, {8.414}}};
*/

// X is the 'electronegativity' or 'Mulliken electronegativity'
// "Towards rapid computational screening of metal-organic frameworks for carbon
//   dioxide capture: Calculation of framework charges via charge equilibration"
// Chemical Engineering Journal, 171(3), 775-781
// C.E. Wilmer and R.Q. Snurr
// X_0 = (\partial dE/\partial Q)_{Q=0}  = 0.5 * (IP + EA)
// the electron affinity can be thought of as the "zero-th" ionization energy
inline double referenceTableX(const ForceField &forceField, std::size_t pseudo_atom_type, std::size_t index)
{
  double electron_affinity, ionization_potential;

  std::size_t oxidationState =
      static_cast<std::size_t>(std::abs(forceField.pseudoAtoms[pseudo_atom_type].oxidationState));
  if (oxidationState == 0)
  {
    ionization_potential = chargeEquilbrationElementData[index].ionizationEnergy[0];
    electron_affinity = chargeEquilbrationElementData[index].electronAffinity;
  }
  else
  {
    ionization_potential = chargeEquilbrationElementData[index].ionizationEnergy[oxidationState];
    electron_affinity = chargeEquilbrationElementData[index].ionizationEnergy[oxidationState - 1];
  }

  return 0.5 * (ionization_potential + electron_affinity);
}

// J is the 'idempotential' of the 'chemical hardness' parameter of the atom
// J_0 = \partial d^2E/\partial Q^2  = IP - EA
// the electron affinity can be thought of as the "zero-th" ionization energy
inline double referenceTableJ(const ForceField &forceField, std::size_t pseudo_atom_type, std::size_t index)
{
  double electron_affinity, ionization_potential;

  std::size_t oxidationState =
      static_cast<std::size_t>(std::abs(forceField.pseudoAtoms[pseudo_atom_type].oxidationState));
  if (oxidationState == 0)
  {
    ionization_potential = chargeEquilbrationElementData[index].ionizationEnergy[0];
    electron_affinity = chargeEquilbrationElementData[index].electronAffinity;
  }
  else
  {
    // TODO: check index and vector-size

    ionization_potential = chargeEquilbrationElementData[index].ionizationEnergy[oxidationState];
    electron_affinity = chargeEquilbrationElementData[index].ionizationEnergy[oxidationState - 1];
  }

  return (ionization_potential - electron_affinity);
}

inline double referenceTableXc(const ForceField &forceField, std::size_t pseudo_atom_type, std::size_t index)
{
  double electron_affinity, ionization_potential;

  std::size_t oxidationState =
      static_cast<std::size_t>(std::abs(forceField.pseudoAtoms[pseudo_atom_type].oxidationState));
  if (oxidationState == 0)
  {
    electron_affinity = chargeEquilbrationElementData[index].electronAffinity;
    ionization_potential = chargeEquilbrationElementData[index].ionizationEnergy[0];
  }
  else
  {
    electron_affinity = chargeEquilbrationElementData[index].ionizationEnergy[oxidationState - 1];
    ionization_potential = chargeEquilbrationElementData[index].ionizationEnergy[oxidationState];
  }

  return (ionization_potential - electron_affinity) * static_cast<double>(oxidationState);
}
}  // namespace ChargeEquilibration
