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
     "RunningEnergy",
     "System",
     "MonteCarlo",
     "MolecularDynamics",
     "PropertyConservedEnergyEvolution",
     "Thermostat",
     "PropertyConventionalRadialDistributionFunction",
     "PropertyMeanSquaredDisplacement",
     "PropertyVelocityAutoCorrelationFunction",
     "PropertyEnthalpy"
]

