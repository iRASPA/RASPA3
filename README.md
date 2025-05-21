<h1 align="center">
  <br>
  <a href="https://github.com/iRASPA/raspa3"><img src="https://avatars.githubusercontent.com/u/46400041?s=200&v=4" alt="Markdownify" width="200"></a>
  <br>
  RASPA3
  <br>
</h1>

<h4 align="center">This software is a general purpose classical simulation package.</h4>

<p align="center">It has been developed at the University of Amsterdam (Amsterdam, The Netherlands) during 2022/2024 in active collaboration with Eindhoven University of Technology (Eindhoven, Netherlands), Delft University of Technology (Delft, The Netherlands), and Northwestern University (Evanston, USA).</p>

<p align="center">
  <a href="https://github.com/iRASPA/raspa3/releases">
<img alt="GitHub Actions Workflow Status" src="https://img.shields.io/github/actions/workflow/status/iRASPA/RASPA3/create-release.yml">
  </a>
  <a href="https://github.com/iRASPA/raspa3/issues"><img alt="GitHub Issues or Pull Requests" src="https://img.shields.io/github/issues/iRASPA/raspa3">
</a>
 <a href="https://iRASPA.github.io/RASPA3"><img alt="Documentation" src="https://img.shields.io/badge/documentation-blue">
 </a>
 <a href="https://github.com/iRASPA/RASPA3/actions/workflows/unit-test.yml"><img alt="Unittests" src="https://github.com/iraspa/raspa3/actions/workflows/unit-test.yml/badge.svg">
 </a>
</p>

<p align="center">
  •
  <a href="#authors">Authors</a> •
  <a href="#contributors">Contributors</a> •
  <a href="#running">Running</a> •
  <a href="#python">Python</a> •
  <a href="#dependencies">Dependencies</a> •
  <a href="#installation-guide">Installation Guide</a> •
</p>

# Authors

Drs. Youri Ran, University of Amsterdam<br>
Drs. Shrinjay Sharma, Delft University of Technology<br>
Dr. Salvador R.G. Balestra, Universidad Pablo de Olavide<br>
Drs. Zhao Li, Northwestern University<br>
Prof. Sofia Calero,  Eindhoven University of Technology<br>
Prof. Thijs Vlugt, Delft University of Technology<br>
Prof. Randall Q. Snurr, Northwestern University<br>
Dr. David Dubbeldam, University of Amsterdam

# Contributors

Alvaro Vazquez Mayagoitia, Argonne National Lab, contribution to openmp-implementation discussion<br>
Anserme, better README.md and packaging

# Citing RASPA3

Y.A. Ran, S. Sharma, S.R.G. Balestra, Z. Li, S. Calero, T.J.H. Vlugt, R.Q. Snurr, D. Dubbeldam, _"RASPA3: A Monte Carlo code for computing adsorption and diffusion in nanoporous materials and thermodynamics properties of fluids"_, **2024**, J. Chem. Phys., 161, 114106, [DOI](https://doi.org/10.1063/5.0226249)

# Running

cd examples/basic/1_mc_methane_in_box<br>
./run


# Implemented so far
- rigid molecules and rigid framework
- multiple systems, box or framework
- grand canonical ensemble (CBMC, CFCMC, and CB/CFCMC)
- Gibbs ensemble (CBMC, CFCMC, and CB/CFCMC)
- Monte Carlo NPT ensemble
- transition matrix Monte Carlo
- Molecular Dynamics NVT ensemble (Nose-Hoover thermostat)
- binary restart
- blocking pockets
- PDB-movies, energy histograms, number of molecule histograms
- RDF, MSD-order-N, VACF, density grids (cube files)
- charge equilibration
- MC/MD hybrid move
- tail-corrections for CFCMC
- grids for rigid frameworks

# Todo-list
- restart-file
- flexible molecules
- flexible frameworks
- reaction ensemble
- identity change
- polarization
- cell-lists for rigid frameworks
- partial molar volumes
- zeo++-type calculations
- partial molar enthalpies/volumes
- optimization
- elastic constants
- HDF5 property writing


# Installation Guide

RASPA3 makes use of modern C++ and requires C++ version 23. This is only included with later versions of compilers and the code is only tested using LLVM & Clang version 18. All prebuilt versions are built with static linking, making it unnecessary for users to install dependencies. 

For contributors or others who need to build from scratch, we recommend loading the supplied Dockerfiles, which include all the dependencies. 

For python users we currently only offer installation via pip, that includes building from source, requiring users to have installed all dependencies.

<!-- TOC -->
* [Prebuilt installation](#download-prebuilt-installation)
* [Building from source](#build-from-source)
    * [Presets](#presets)
    * [Linux](#linux)
      * [Ubuntu 24](#ubuntu-24)
      * [Fedora 40](#fedora-40)
    * [macOS](#macos)
    * [Python](#python)
<!-- TOC -->

## Download prebuilt installation

In Github, on the right side of the page, you will find the releases section. Select your OS and use the installer to install RASPA3.
We provide packages for:<br>
- arm64<br>
- Intel/AMD core-avx2<br>
- Intel skylake-avx512<br>
<br>
Use the 'core-avx2' version for Intel and AMD cpu's that support the avx2 instruction set:<br>
- Intel: Haswell processors 2013 and newer<br>
- AMD: Excavator processors 2015 and newer<br>
<br>
The 'skylake-avx512' binary can provide a roughly 25% increase in speed compared to the 'core-avx2' on Intel cpu's that support avx512.

## Build from source

### Dependencies

- cmake 3.28 or higher<br>
- ninja 1.11 or higher<br>
- llvm-18 or higher<br>
- python3 and pybind11<br>
- blas and lapack (64-bit integers)
- openmp
- hdf5

### Installing dependencies with Conda

The most practical method of compiling raspa from source is through installing the dependencies with conda. This is readily achieved through:

```bash
conda env create -f env.yml
conda acitvate raspa
cmake --preset=linux_conda   # or -preset=mac_conda
ninja -C build
```

### presets

cmake --list-presets<br>
Available configure presets:<br>
> "windows_conda_raspa3"<br>
> "mac_conda"<br>
> "mac_conda_raspa3"<br>
> "linux_conda"<br>
> "linux_conda_raspa3"<br>
> "windows_conda_raspalib"<br>
> "mac_conda_raspalib"<br>
> "linux_conda_raspalib"<br>
> "windows-arm64-package"<br>
> "windows-x64-core-avx2-package"<br>
> "windows-x64-skylake-avx512-package"<br>
> "macos-x64-core-avx2"<br>
> "macos-x64-core-avx2-package"<br>
> "macos-x64-skylake-avx512"<br>
> "macos-x64-skylake-avx512-package"<br>
> "macos-x64-debug"<br>
> "macos-x64-profile"<br>
> "macos-apple-silicon"<br>
> "macos-apple-silicon-package"<br>
> "macos-apple-silicon-profile"<br>
> "macos-apple-silicon-debug"<br>
> "linux-x86_64"<br>
> "linux-x86_64-carbon"<br>
> "linux-x86_64-core-avx2-opensuse-leap-15.2"<br>
> "linux-x86_64-skylake-avx512-opensuse-leap-15.2"<br>
> "linux-x86_64-core-avx2-opensuse-leap-15.3"<br>
> "linux-x86_64-skylake-avx512-opensuse-leap-15.3"<br>
> "linux-x86_64-core-avx2-opensuse-leap-15.4"<br>
> "linux-x86_64-skylake-avx512-opensuse-leap-15.4"<br>
> "linux-x86_64-core-avx2-opensuse-leap-15.5"<br>
> "linux-x86_64-skylake-avx512-opensuse-leap-15.5"<br>
> "linux-x86_64-core-avx2-opensuse-leap-15.6"<br>
> "linux-x86_64-skylake-avx512-opensuse-leap-15.6"<br>
> "linux-x86_64-core-avx2-opensuse-tumbleweed"<br>
> "linux-x86_64-skylake-avx512-opensuse-tumbleweed"<br>
> "linux-x86_64-core-avx2-archlinux"<br>
> "linux-x86_64-skylake-avx512-archlinux"<br>
> "linux-x86_64-core-avx2-redhat-6"<br>
> "linux-x86_64-skylake-avx512-redhat-6"<br>
> "linux-x86_64-core-avx2-redhat-7"<br>
> "linux-x86_64-skylake-avx512-redhat-7"<br>
> "linux-x86_64-core-avx2-redhat-8"<br>
> "linux-x86_64-skylake-avx512-redhat-8"<br>
> "linux-x86_64-core-avx2-redhat-9"<br>
> "linux-x86_64-skylake-avx512-redhat-9"<br>
> "linux-x86_64-core-avx2-debian-12"<br>
> "linux-x86_64-skylake-avx512-debian-12"<br>
> "linux-x86_64-core-avx2-debian-11"<br>
> "linux-x86_64-skylake-avx512-debian-11"<br>
> "linux-x86_64-core-avx2-debian-10"<br>
> "linux-x86_64-skylake-avx512-debian-10"<br>
> "linux-x86_64-core-avx2-ubuntu-24"<br>
> "linux-x86_64-skylake-avx512-ubuntu-24"<br>
> "linux-x86_64-core-avx2-ubuntu-22"<br>
> "linux-x86_64-skylake-avx512-ubuntu-22"<br>
> "linux-x86_64-core-avx2-ubuntu-20"<br>
> "linux-x86_64-skylake-avx512-ubuntu-20"<br>
> "linux-x86_64-core-avx2-fedora-35"<br>
> "linux-x86_64-skylake-avx512-fedora-35"<br>
> "linux-x86_64-core-avx2-fedora-36"<br>
> "linux-x86_64-skylake-avx512-fedora-36"<br>
> "linux-x86_64-core-avx2-fedora-37"<br>
> "linux-x86_64-skylake-avx512-fedora-37"<br>
> "linux-x86_64-core-avx2-fedora-38"<br>
> "linux-x86_64-skylake-avx512-fedora-38"<br>
> "linux-x86_64-core-avx2-fedora-39"<br>
> "linux-x86_64-skylake-avx512-fedora-39"<br>
> "linux-x86_64-core-avx2-fedora-40"<br>
> "linux-x86_64-skylake-avx512-fedora-40"<br>
> "linux-x86_64-core-avx2-fedora-41"<br>
> "linux-x86_64-skylake-avx512-fedora-41"<br>
> "linux-aarch64"<br>
> "linux-aarch64-ubuntu-24"<br>
> "linux-aarch64-redhat-9"

### linux

#### Ubuntu 24

apt-get install -y --no-install-recommends git ca-certificates cmake ninja-build<br>
apt-get install -y --no-install-recommends llvm lld clang clang-tools clang-tidy libc++-dev libc++abi-dev libomp-dev libclang-rt-dev<br>
apt-get install -y --no-install-recommends python3 pybind11-dev python3-pybind11 python3-dev<br>
apt-get install -y --no-install-recommends liblapack64-dev libblas64-dev<br>
cmake -B build --preset=linux-ubuntu-24<br>
ninja -C build<br>
ninja -C build install<br>
ctest --test-dir build/tests/raspakit-tests --verbose

#### Fedora 40

dnf install -y wget git rpm-build<br>
dnf install -y llvm lld cmake clang clang-tools-extra ninja-build<br>
dnf install -y libomp-devel libcxx libcxxabi libcxx-devel libcxxabi-devel libcxx-static libcxxabi-static<br>
dnf install -y lapack-devel lapack64 blas64<br>
dnf install -y python3 python3-devel python3-pybind11<br>
dnf install -y pybind11-devel<br>
cmake -B --preset linux-fedora-40<br>
ninja -C build<br>
ninja -C build install<br>
ctest --test-dir build/tests/raspakit-tests --verbose

### macOS
brew install llvm lld libomp hdf5 libaec ninja cmake doxygen graphviz lapack pybind11<br>
(make sure '\`brew --prefix hdf5\`/bin'  for hdf5 is in your path)<br>
cmake -B --preset macos-apple-silicon<br>
(or cmake -B --preset macos-x64)<br>
ninja -C build<br>
ninja -C build install<br>
ctest --test-dir build/tests/raspakit-tests --verbose


### python

This package can also be built as a library for python. To build the python package the pip packaging system can be used. Note that due to compilation of the full package this might take a few minutes. To install, run the following command:

```bash
export CMAKE_PRESET="--preset=python"
pip install .
```

This will install the package to the current python environment.

We strongly advise users to use the CMakePresets.json preset for their given system. For building the python package with a given preset change the following line to reflect the given preset:

```
export CMAKE_PRESET="--preset=macos-apple-silicon"
pip install .
```
