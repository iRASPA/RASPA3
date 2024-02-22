from RaspaKit import *

atomTypes = [PseudoAtom("Si",    28.0855,   2.05,  14, False), \
             PseudoAtom("O",     15.999,   -1.025,  8, False), \
             PseudoAtom("CH4",   16.04246,  0.0,    6, False), \
             PseudoAtom("C_co2", 12.0,      0.6512, 6, False), \
             PseudoAtom("O_co2", 15.9994,  -0.3256, 8, False)]
parameters = [VDWParameters(22.0 / 1.2027242847, 2.30), \
              VDWParameters(53.0 / 1.2027242847, 3.3), \
              VDWParameters(158.5 / 1.2027242847, 3.72), \
              VDWParameters(29.933 / 1.2027242847, 2.745), \
              VDWParameters(85.671 / 1.2027242847, 3.017)]

force_field = ForceField(atomTypes, parameters, ForceField.MixingRule.Lorentz_Berthelot, 12.0, True, False)

print(force_field)

framework = Component(0, "ITQ-29", 46144.748974602669, SimulationBox(11.8671, 11.8671, 11.8671), \
                      517, \
                      [ \
                        Atom(double3(0.3683, 0.1847, 0),       2.05,  1.0, 0, 0, 0), \
                        Atom(double3(0.5,    0.2179, 0),      -1.025, 1.0, 1, 0, 0), \
                        Atom(double3(0.2939, 0.2939, 0),      -1.025, 1.0, 1, 0, 0), \
                        Atom(double3(0.3429, 0.1098, 0.1098), -1.025, 1.0, 1, 0, 0) \
                      ], \
                      int3(1, 1, 1), 5, 21)

component = Component(1, \
                      "CO2", \
                      43.9988, \
                      SimulationBox(0.0, 0.0, 0.0), \
                      304.1282, 7377300.0, 0.22394, \
                      [ \
                         Atom(double3(0.0, 0.0,  1.149), -0.3256, 1.0, 4, 1, 0), \
                         Atom(double3(0.0, 0.0,  0.0  ),  0.6512, 1.0, 3, 1, 0), \
                         Atom(double3(0.0, 0.0, -1.149), -0.3256, 1.0, 4, 1, 0) \
                      ], 5, 21)


system = System(0, 300.0, 1e4, force_field, [ framework, component ], [ 0, 2 ], 5)

system.atomPositions[72].position = double3(5.93355, 7.93355, 5.93355 + 1.149);
system.atomPositions[73].position = double3(5.93355, 7.93355, 5.93355 + 0.0);
system.atomPositions[74].position = double3(5.93355, 7.93355, 5.93355 - 1.149);
system.atomPositions[75].position = double3(5.93355, 3.93355, 5.93355 + 1.149);
system.atomPositions[76].position = double3(5.93355, 3.93355, 5.93355 + 0.0);
system.atomPositions[77].position = double3(5.93355, 3.93355, 5.93355 - 1.149);

#print(system.atomPositions)

energy = system.computeTotalEnergies()
print(1.2027242847 * energy.moleculeMoleculeVDW)
