import enum
import collections.abc
import sys
import platform
import ctypes
import subprocess

def _get_cpu_feature_level():
    machine = platform.machine().lower()
    if "x86_64" not in machine and "amd64" not in machine:
        return "base"

    try:
        if platform.system() == "Windows":
            # PF_AVX512F_INSTRUCTIONS_AVAILABLE = 41
            # PF_AVX2_INSTRUCTIONS_AVAILABLE = 40
            if ctypes.windll.kernel32.IsProcessorFeaturePresent(41):
                return "avx512"
            if ctypes.windll.kernel32.IsProcessorFeaturePresent(40):
                return "avx2"
        
        elif platform.system() == "Linux":
            with open("/proc/cpuinfo", "r") as f:
                content = f.read()
                if "avx512f" in content: return "avx512"
                if "avx2" in content: return "avx2"
        
        elif platform.system() == "Darwin":
            cmd = ["sysctl", "-n", "machdep.cpu.leaf7_features"]
            output = subprocess.check_output(cmd).decode().lower()
            if "avx512f" in output: return "avx512"
            if "avx2" in output: return "avx2"
            
    except Exception:
        pass
    return "base"

# Specialized Loader
level = _get_cpu_feature_level()
try:
    if level == "avx512":
        from raspalib.raspalib_avx512 import *
    elif level == "avx2":
        from raspalib.raspalib_avx2 import *
    else:
        from raspalib.raspalib import *
except ImportError:
    # Fallback if a specific optimized binary wasn't shipped
    from raspalib.raspalib import *

from raspalib.version import __version__
from raspalib.atom import *
from raspalib.pseudo_atom import *
from raspalib.vdw_parameters import *
from raspalib.force_field import *
from raspalib.mc_move_probabilities import *
from raspalib.connectivity_table import *
from raspalib.intra_molecular_potentials import *
from raspalib.move import *
from raspalib.move_statistics_double import *
from raspalib.move_statistics_double3 import *
from raspalib.mc_move_statistics import *
from raspalib.widom_data import *
from raspalib.loading_data import *
from raspalib.pressure_data import *
from raspalib.sample_movie import *
from raspalib.energy_factor import *
from raspalib.energy_status import *
from raspalib.average_energy_type import *
from raspalib.property_energy import *
from raspalib.property_energy_histogram import *
from raspalib.property_loading import *
from raspalib.property_number_of_molecules_histogram import *
from raspalib.enthalpy_of_adsorption_data import *
from raspalib.property_enthalpy import *
from raspalib.property_pressure import *
from raspalib.property_number_of_molecules_evolution import *
from raspalib.property_volume_evolution import *
from raspalib.running_energy import *
from raspalib.property_conserved_energy_evolution import *
from raspalib.property_density_grid import *
from raspalib.property_lambda_probability_histogram import *
from raspalib.property_widom import *
from raspalib.property_conventional_radial_distribution_function import *
from raspalib.property_mean_squared_displacement import *
from raspalib.property_velocity_autocorrelation_function import *
from raspalib.thermostat import *
from raspalib.simulation_box import *
from raspalib.component import *
from raspalib.framework import *
from raspalib.system import *
from raspalib.monte_carlo import *
from raspalib.molecular_dynamics import *
from raspalib.equation_of_state import *

__all__ = [
    "PseudoAtom",
    "VDWParameters",
    "Atom",
    "ForceField",
    "MCMoveProbabilities",
    "MCMoveStatistics",
    "IntraMolecularPotentials",
    "Component",
    "SimulationBox",
    "Framework",
    "System",
    "MonteCarlo",
    "LoadingData",
    "SampleMovie",
    "EnergyFactor",
    "EnergyStatus",
    "AverageEnergyType",
    "PropertyEnergy",
    "PropertyEnergyHistogram",
    "PropertyNumberOfMoleculesHistogram",
    "PropertyLoading",
    "PropertyNumberOfMoleculesEvolution",
    "PropertyVolumeEvolution",
    "PropertyDensityGrid",
    "PropertyLambdaProbabilityHistogram",
    "Move",
    "MoveStatisticsDouble3",
    "WidomData",
    "EquationOfState",
    "MolecularDynamics",
    "RunningEnergy",
    "PropertyConservedEnergyEvolution",
    "Thermostat",
    "PropertyConventionalRadialDistributionFunction",
    "PropertyMeanSquaredDisplacement",
    "PropertyVelocityAutoCorrelationFunction"
]

