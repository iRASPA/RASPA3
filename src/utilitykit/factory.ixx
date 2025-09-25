module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <complex>
#include <cstddef>
#include <numbers>
#include <span>
#include <vector>
#endif

export module factory;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import int3;
import double3;
import double3x3;

import units;
import atom;
import pseudo_atom;
import vdwparameters;
import forcefield;
import framework;
import component;
import system;
import simulationbox;
import energy_factor;
import running_energy;
import gradient_factor;
import energy_status;
import units;
import interactions_intermolecular;
import interactions_framework_molecule;
import interactions_ewald;

export namespace TestFactories
{
inline ForceField makeDefaultFF(double rc = 12.0, bool shifted = true, bool tailCorrections = false,
                                bool useEwald = false)
{
  return ForceField({{"Si", true, 28.0855, 2.05, 0.0, 14, false},
                     {"O", true, 15.999, -1.025, 0.0, 8, false},
                     {"CH4", false, 16.04246, 0.0, 0.0, 6, false},
                     {"C_co2", false, 12.0, 0.6512, 0.2, 6, false},
                     {"O_co2", false, 15.9994, -0.3256, 0.1, 8, false},
                     {"Na+", false, 12.0, 1.0, 0.0, 6, false},
                     {"Cl-", false, 15.9994, -1.0, 0.0, 8, false},
                     {"Ow", false, 15.9996, 0.0, 0.0, 8, false},
                     {"Hw", false, 1.0008, 0.241, 0.0, 1, false},
                     {"Lw", false, 0.0, -0.241, 0.0, 0, false}},

                    {{22.0, 2.30},
                     {53.0, 3.30},
                     {158.5, 3.72},
                     {29.933, 2.745},
                     {85.671, 3.017},
                     {15.0966, 2.65755},
                     {142.562, 3.51932},
                     {89.633, 3.097},
                     {0.0, 1.0},
                     {0.0, 1.0}},
                    ForceField::MixingRule::Lorentz_Berthelot, rc, rc, rc, shifted, tailCorrections, useEwald);
}

inline Component makeMethane(const ForceField& ff, std::uint8_t id = 0)
{
  return Component(id, ff, "methane", 190.564, 45599200, 0.01142, {Atom({0, 0, 0}, 0.0, 1.0, 0, 2, id, false, false)},
                   {}, {}, 5, 21);
}

inline Component makeCO2(const ForceField& ff, std::uint8_t id = 0, bool useCharges = false)
{
  const double qC = useCharges ? 0.6512 : 0.0;
  const double qO = useCharges ? -0.3256 : 0.0;

  return Component(
      id, ff, "CO2", 304.1282, 7377300.0, 0.22394,
      {Atom({0, 0, 1.149}, qO, 1.0, 0, 4, id, false, false), Atom({0, 0, 0.000}, qC, 1.0, 0, 3, id, false, false),
       Atom({0, 0, -1.149}, qO, 1.0, 0, 4, id, false, false)},
      {}, {}, 5, 21);
}

inline Component makeWater(const ForceField& ff, std::uint8_t id = 0, bool useCharges = false)
{
  const double qh = useCharges ? 0.241 : 0.0;
  const double ql = useCharges ? -0.241 : 0.0;

  return Component(
      id, ff, "water", 0.0, 0.0, 0.0,
      {Atom(double3(0.0, 0.0, 0.0), 0.0, 1.0, 0, 7, id, false, false),
       Atom(double3(-0.75695032726366118157, 0.0, -0.58588227661829499395), qh, 1.0, 0, 8, id, false, false),
       Atom(double3(0.75695032726366118157, 0.0, -0.58588227661829499395), qh, 1.0, 0, 8, id, false, false),
       Atom(double3(0.0, -0.57154330164408200866, 0.40415127656087122858), ql, 1.0, 0, 9, id, false, false),
       Atom(double3(0.0, 0.57154330164408200866, 0.40415127656087122858), ql, 1.0, 0, 9, id, false, false)},
      {}, {}, 5, 21);
}

inline Component makeIon(const ForceField& ff, std::uint8_t id, std::string_view name, std::uint16_t type, double q)
{
  return Component(id, ff, std::string{name}, 0.0, 0.0, 0.0, {Atom({0, 0, 0}, q, 1.0, 0, type, id, false, false)}, {},
                   {}, 5, 21);
}

inline Framework makeFAU(const ForceField& ff, int3 replicate = {1, 1, 1})
{
  return Framework(0, ff, "FAU", SimulationBox(24.2576, 24.2576, 24.2576), 526,
                   {Atom({-0.05392, 0.1253, 0.03589}, 2.05, 1, 0, 0, 0, false, false),
                    Atom({0, -0.10623, 0.10623}, -1.025, 1, 0, 1, 0, false, false),
                    Atom({-0.00323, -0.00323, 0.14066}, -1.025, 1, 0, 1, 0, false, false),
                    Atom({0.0757, 0.0757, -0.03577}, -1.025, 1, 0, 1, 0, false, false),
                    Atom({0.07063, 0.07063, 0.32115}, -1.025, 1, 0, 1, 0, false, false)},
                   replicate);
}

inline Framework makeITQ29(const ForceField& ff, int3 replicate = {1, 1, 1})
{
  return Framework(0, ff, "ITQ-29", SimulationBox(11.8671, 11.8671, 11.8671), 517,

                   {Atom({0.3683, 0.1847, 0}, 2.05, 1, 0, 0, 0, false, false),
                    Atom({0.5, 0.2179, 0}, -1.025, 1, 0, 1, 0, false, false),
                    Atom({0.2939, 0.2939, 0}, -1.025, 1, 0, 1, 0, false, false),
                    Atom({0.3429, 0.1098, 0.1098}, -1.025, 1, 0, 1, 0, false, false)},
                   replicate);
}

inline Framework makeMFI_Si(const ForceField& ff, int3 replicate = {1, 1, 1})
{
  return Framework(0, ff, "MFI_SI", SimulationBox(20.022, 19.899, 13.383), 292,
                   {
                       Atom(double3(0.42238, 0.0565, -0.33598), 2.05, 1.0, 0, 0, 0, false, false),
                       Atom(double3(0.30716, 0.02772, -0.1893), 2.05, 1.0, 0, 0, 0, false, false),
                       Atom(double3(0.27911, 0.06127, 0.0312), 2.05, 1.0, 0, 0, 0, false, false),
                       Atom(double3(0.12215, 0.06298, 0.0267), 2.05, 1.0, 0, 0, 0, false, false),
                       Atom(double3(0.07128, 0.02722, -0.18551), 2.05, 1.0, 0, 0, 0, false, false),
                       Atom(double3(0.18641, 0.05896, -0.32818), 2.05, 1.0, 0, 0, 0, false, false),
                       Atom(double3(0.42265, -0.1725, -0.32718), 2.05, 1.0, 0, 0, 0, false, false),
                       Atom(double3(0.30778, -0.13016, -0.18548), 2.05, 1.0, 0, 0, 0, false, false),
                       Atom(double3(0.27554, -0.17279, 0.03109), 2.05, 1.0, 0, 0, 0, false, false),
                       Atom(double3(0.12058, -0.1731, 0.02979), 2.05, 1.0, 0, 0, 0, false, false),
                       Atom(double3(0.07044, -0.13037, -0.182), 2.05, 1.0, 0, 0, 0, false, false),
                       Atom(double3(0.18706, -0.17327, -0.31933), 2.05, 1.0, 0, 0, 0, false, false),
                       Atom(double3(0.3726, 0.0534, -0.2442), -1.025, 1.0, 0, 1, 0, false, false),
                       Atom(double3(0.3084, 0.0587, -0.0789), -1.025, 1.0, 0, 1, 0, false, false),
                       Atom(double3(0.2007, 0.0592, 0.0289), -1.025, 1.0, 0, 1, 0, false, false),
                       Atom(double3(0.0969, 0.0611, -0.0856), -1.025, 1.0, 0, 1, 0, false, false),
                       Atom(double3(0.1149, 0.0541, -0.2763), -1.025, 1.0, 0, 1, 0, false, false),
                       Atom(double3(0.2435, 0.0553, -0.246), -1.025, 1.0, 0, 1, 0, false, false),
                       Atom(double3(0.3742, -0.1561, -0.2372), -1.025, 1.0, 0, 1, 0, false, false),
                       Atom(double3(0.3085, -0.1552, -0.0728), -1.025, 1.0, 0, 1, 0, false, false),
                       Atom(double3(0.198, -0.1554, 0.0288), -1.025, 1.0, 0, 1, 0, false, false),
                       Atom(double3(0.091, -0.1614, -0.0777), -1.025, 1.0, 0, 1, 0, false, false),
                       Atom(double3(0.1169, -0.1578, -0.2694), -1.025, 1.0, 0, 1, 0, false, false),
                       Atom(double3(0.2448, -0.1594, -0.2422), -1.025, 1.0, 0, 1, 0, false, false),
                       Atom(double3(0.3047, -0.051, -0.1866), -1.025, 1.0, 0, 1, 0, false, false),
                       Atom(double3(0.0768, -0.0519, -0.1769), -1.025, 1.0, 0, 1, 0, false, false),
                       Atom(double3(0.4161, 0.1276, -0.3896), -1.025, 1.0, 0, 1, 0, false, false),
                       Atom(double3(0.4086, -0.0017, -0.4136), -1.025, 1.0, 0, 1, 0, false, false),
                       Atom(double3(0.402, -0.1314, -0.4239), -1.025, 1.0, 0, 1, 0, false, false),
                       Atom(double3(0.1886, 0.1298, -0.3836), -1.025, 1.0, 0, 1, 0, false, false),
                       Atom(double3(0.194, 0.0007, -0.4082), -1.025, 1.0, 0, 1, 0, false, false),
                       Atom(double3(0.1951, -0.1291, -0.419), -1.025, 1.0, 0, 1, 0, false, false),
                       Atom(double3(-0.0037, 0.0502, -0.208), -1.025, 1.0, 0, 1, 0, false, false),
                       Atom(double3(-0.004, -0.1528, -0.2078), -1.025, 1.0, 0, 1, 0, false, false),
                       Atom(double3(0.4192, -0.25, -0.354), -1.025, 1.0, 0, 1, 0, false, false),
                       Atom(double3(0.1884, -0.25, -0.3538), -1.025, 1.0, 0, 1, 0, false, false),
                       Atom(double3(0.2883, -0.25, 0.0579), -1.025, 1.0, 0, 1, 0, false, false),
                       Atom(double3(0.1085, -0.25, 0.0611), -1.025, 1.0, 0, 1, 0, false, false),
                   },
                   replicate);
}

inline Framework makeCHA(const ForceField& ff, int3 replicate = {1, 1, 1})
{
  return Framework(0, ff, "CHA",
                   SimulationBox(9.459, 9.459, 9.459, 94.07 * std::numbers::pi / 180.0,
                                 94.07 * std::numbers::pi / 180.0, 94.07 * std::numbers::pi / 180.0),
                   1,
                   {// double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type,
                    // uint8_t componentId, uint8_t groupId
                    Atom(double3(0.10330000, 0.33310000, 0.87430000), 2.05, 1.0, 0, 0, 0, false, false),
                    Atom(double3(0.12570000, 0.89670000, 0.66690000), 2.05, 1.0, 0, 0, 0, false, false),
                    Atom(double3(0.66690000, 0.89670000, 0.12570000), 2.05, 1.0, 0, 0, 0, false, false),
                    Atom(double3(0.10330000, 0.87430000, 0.33310000), 2.05, 1.0, 0, 0, 0, false, false),
                    Atom(double3(0.87430000, 0.33310000, 0.10330000), 2.05, 1.0, 0, 0, 0, false, false),
                    Atom(double3(0.33310000, 0.87430000, 0.10330000), 2.05, 1.0, 0, 0, 0, false, false),
                    Atom(double3(0.12570000, 0.66690000, 0.89670000), 2.05, 1.0, 0, 0, 0, false, false),
                    Atom(double3(0.33310000, 0.10330000, 0.87430000), 2.05, 1.0, 0, 0, 0, false, false),
                    Atom(double3(0.66690000, 0.12570000, 0.89670000), 2.05, 1.0, 0, 0, 0, false, false),
                    Atom(double3(0.87430000, 0.10330000, 0.33310000), 2.05, 1.0, 0, 0, 0, false, false),
                    Atom(double3(0.89670000, 0.12570000, 0.66690000), 2.05, 1.0, 0, 0, 0, false, false),
                    Atom(double3(0.89670000, 0.66690000, 0.12570000), 2.05, 1.0, 0, 0, 0, false, false),
                    Atom(double3(0.00000000, 0.73350000, 0.26650000), -1.025, 1.0, 0, 1, 0, false, false),
                    Atom(double3(0.26650000, 0.73350000, 1.00000000), -1.025, 1.0, 0, 1, 0, false, false),
                    Atom(double3(0.26650000, 1.00000000, 0.73350000), -1.025, 1.0, 0, 1, 0, false, false),
                    Atom(double3(0.73350000, 1.00000000, 0.26650000), -1.025, 1.0, 0, 1, 0, false, false),
                    Atom(double3(1.00000000, 0.26650000, 0.73350000), -1.025, 1.0, 0, 1, 0, false, false),
                    Atom(double3(0.73350000, 0.26650000, 1.00000000), -1.025, 1.0, 0, 1, 0, false, false),
                    Atom(double3(0.15060000, 0.84940000, 0.50000000), -1.025, 1.0, 0, 1, 0, false, false),
                    Atom(double3(0.50000000, 0.84940000, 0.15060000), -1.025, 1.0, 0, 1, 0, false, false),
                    Atom(double3(0.15060000, 0.50000000, 0.84940000), -1.025, 1.0, 0, 1, 0, false, false),
                    Atom(double3(0.84940000, 0.15060000, 0.50000000), -1.025, 1.0, 0, 1, 0, false, false),
                    Atom(double3(0.84940000, 0.50000000, 0.15060000), -1.025, 1.0, 0, 1, 0, false, false),
                    Atom(double3(0.50000000, 0.15060000, 0.84940000), -1.025, 1.0, 0, 1, 0, false, false),
                    Atom(double3(0.25030000, 0.89300000, 0.25030000), -1.025, 1.0, 0, 1, 0, false, false),
                    Atom(double3(0.10700000, 0.74970000, 0.74970000), -1.025, 1.0, 0, 1, 0, false, false),
                    Atom(double3(0.25030000, 0.25030000, 0.89300000), -1.025, 1.0, 0, 1, 0, false, false),
                    Atom(double3(0.89300000, 0.25030000, 0.25030000), -1.025, 1.0, 0, 1, 0, false, false),
                    Atom(double3(0.74970000, 0.10700000, 0.74970000), -1.025, 1.0, 0, 1, 0, false, false),
                    Atom(double3(0.74970000, 0.74970000, 0.10700000), -1.025, 1.0, 0, 1, 0, false, false),
                    Atom(double3(0.02040000, 0.31930000, 0.02040000), -1.025, 1.0, 0, 1, 0, false, false),
                    Atom(double3(0.68070000, 0.97960000, 0.97960000), -1.025, 1.0, 0, 1, 0, false, false),
                    Atom(double3(0.02040000, 0.02040000, 0.31930000), -1.025, 1.0, 0, 1, 0, false, false),
                    Atom(double3(0.31930000, 0.02040000, 0.02040000), -1.025, 1.0, 0, 1, 0, false, false),
                    Atom(double3(0.97960000, 0.68070000, 0.97960000), -1.025, 1.0, 0, 1, 0, false, false),
                    Atom(double3(0.97960000, 0.97960000, 0.68070000), -1.025, 1.0, 0, 1, 0, false, false)},
                   replicate);
}

}  // namespace TestFactories
