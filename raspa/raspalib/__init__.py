import enum
import collections.abc
import os
import platform
import ctypes
import ctypes.util
import warnings

def _has_opencl():
    # 1. Standard search (Windows/macOS/Linux dev)
    if ctypes.util.find_library("OpenCL"):
        return True

    # 2. Linux-specific deep search for runtime versions
    if platform.system() == "Linux":
        patterns = [
            "libOpenCL.so.1", # Try loading by name via OS loader
            "/usr/lib/x86_64-linux-gnu/libOpenCL.so.1",
            "/usr/lib64/libOpenCL.so.1",
            "/usr/lib/libOpenCL.so.1"
        ]
        for p in patterns:
            if p.startswith("/"): # Absolute path check
                if os.path.exists(p): return True
            else: # Name-only check using OS dynamic loader
                try:
                    ctypes.CDLL(p)
                    return True
                except Exception:
                    continue
    return False


__has_opencl__ = _has_opencl()

if not __has_opencl__ and os.environ.get("RASPALIB_QUIET") != "1":
    distro_cmd = "Please install OpenCL drivers."
    sys_name = platform.system()
    
    if sys_name == "Windows":
        distro_cmd = "Please update your GPU drivers (NVIDIA, AMD, or Intel) to install OpenCL.dll."
    elif sys_name == "Linux":
        try:
            info = platform.freedesktop_os_release()
            d_id = info.get("ID", "").lower()
            if "ubuntu" in d_id or "debian" in d_id:
                distro_cmd = "Run: 'sudo apt install ocl-icd-libopencl1'"
            elif "arch" in d_id:
                distro_cmd = "Run: 'sudo pacman -S ocl-icd'"
            elif "fedora" in d_id or "rhel" in d_id or "centos" in d_id:
                distro_cmd = "Run: 'sudo dnf install ocl-icd'"
        except:
            pass
    elif sys_name == "Darwin":
        distro_cmd = "OpenCL should be included with macOS. Check your system updates."

    warnings.warn(
        f"OpenCL library not found. GPU acceleration is non-functional. {distro_cmd}",
        RuntimeWarning
    )

from . import _raspalib as _module

globals().update({
    name: export 
    for name, export in _module.__dict__.items() 
    if not name.startswith('_')
})
__file_imported__ = _module.__file__
del _module

def print_debug_info():
    print(f"Raspalib Info:")
    print(f"  Binary Path:     {__file_imported__}")
    print(f"  OpenCL detected: {__has_opencl__}")


from raspalib.version import __version__

__all__ = [
    "PseudoAtom",
    "VDWParameters",
    "Atom",
    "ForceField",
    "MCMoveProbabilities",
    "MCMoveStatistics",
    "ConnectivityTable",
    "IntraMolecularPotentials",
    "Component",
    "SimulationBox",
    "Framework",
    "System",
    "MonteCarlo",
    "LoadingData",
    "SampleMovie",
    "EnergyDuDlambda",
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
    "PropertyWidom",
    "PropertyPressure",
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
    "PropertyVelocityAutoCorrelationFunction",
    "PropertyEnthalpy",
    "EnthalpyOfAdsorptionData"
]

