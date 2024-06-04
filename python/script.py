from _raspa import *

atomTypes = [
    PseudoAtom("Si", 28.0855, 2.05, 14, False),
    PseudoAtom("O", 15.999, -1.025, 8, False),
    PseudoAtom("CH4", 16.04246, 0.0, 6, False),
    PseudoAtom("C_co2", 12.0, 0.6512, 6, False),
    PseudoAtom("O_co2", 15.9994, -0.3256, 8, False)
]
parameters = [
    VDWParameters(22.0, 2.30),
    VDWParameters(53.0, 3.3),
    VDWParameters(158.5, 3.72),
    VDWParameters(29.933, 2.745),
    VDWParameters(85.671, 3.017)
]

force_field = ForceField(atomTypes, parameters, ForceField.MixingRule.Lorentz_Berthelot, 12.0, True, False)

print(force_field)

framework = Framework(
    0,
    force_field,
    "ITQ-29",
    SimulationBox(11.8671, 11.8671, 11.8671),
    517,
    [
        Atom(double3(0.3683, 0.1847, 0), 2.05, 1.0, 0, 0, 0, 0),
        Atom(double3(0.5, 0.2179, 0), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.2939, 0.2939, 0), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.3429, 0.1098, 0.1098), -1.025, 1.0, 0, 1, 0, 0)
    ],
    int3(1, 1, 1),
)

component = Component(
    0,
    force_field,
    "CO2",
    304.1282,
    7377300.0,
    0.22394,
    [
        Atom(double3(0.0, 0.0, 1.149), -0.3256, 1.0, 0, 4, 1, 0),
        Atom(double3(0.0, 0.0, 0.0), 0.6512, 1.0, 0, 3, 1, 0),
        Atom(double3(0.0, 0.0, -1.149), -0.3256, 1.0, 0, 4, 1, 0)
    ],
    5,
    21,
    MCMoveProbabilitiesParticles(probabilityTranslationMove=1.0, probabilityRotationMove=1.0),
)

system = System(0, None, 300.0, 1e4, force_field, [framework], [component], [2], 5)

energy = system.computeTotalEnergies()
print(energy.moleculeMoleculeVDW)
