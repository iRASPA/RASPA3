# Compiling `RASPA`
\page compiling Compiling

## Requirements

`CMake 3.28` and later support `C++` modules. `C++20` named modules are
now supported by the `Ninja Generator 1.11` or newer in combination with
the `LLVM/Clang 16.0` and newer, `MSVC toolset 14.34` and newer, or
`GCC 14` and newer. We recommend LLVM/Clang 18 or higher for compiling
`RASPA3`. This version has support for `std::format`, `std::print`, and
`std::jthread`. The output-files of `RASPA3` are in `UTF-8` encoding and
contain unicode characters (e.g. for Å). `RASPA3` depends on

-   `C++23` compliant compiler
-   `Cmake 3.28`
-   `Ninja Generator 1.11`
-   `Openmp`
-   `hdf5`
-   `lapack` and `blas` (64-bit integers)
-   `pybind11`

## `RASPA` from `git`

Working with '`git`' and a remote repository means that you will have to
distinguish between two locations of the code:

1.  The repository (visible to everyone)

2.  your local copy (only visible to you)

To check-out the code for the first time do:

         git clone https://github.com/iraspa/RASPA3

After that, you can update the code by using

         git pull

## Compiling `RASPA`

Using `cmake --list-presets` you can see the list of `CMake` presets.
These include

-  `python`
-  `windows_conda_raspa3`
-  `mac_conda_raspa3`
-  `linux_conda_raspa3`
-  `windows_conda_raspalib`
-  `mac_conda_raspalib`
-  `linux_conda_raspalib`
-  `windows-arm64-package`
-  `windows-x64-core-avx2-package`
-  `windows-x64-skylake-avx512-package`
-  `macos-universal`
-  `macos-universal-package`
-  `macos-x64-core-avx2`
-  `macos-x64-core-avx2-package`
-  `macos-x64-skylake-avx512`
-  `macos-x64-skylake-avx512-package`
-  `macos-x64-debug`
-  `macos-x64-profile`
-  `macos-apple-silicon`
-  `macos-apple-silicon-package`
-  `macos-apple-silicon-profile`
-  `macos-apple-silicon-debug`
-  `linux-x86_64`
-  `linux-x86_64-carbon`
-  `linux-x86_64-core-avx2-opensuse-leap-15.2`
-  `linux-x86_64-skylake-avx512-opensuse-leap-15.2`
-  `linux-x86_64-core-avx2-opensuse-leap-15.3`
-  `linux-x86_64-skylake-avx512-opensuse-leap-15.3`
-  `linux-x86_64-core-avx2-opensuse-leap-15.4`
-  `linux-x86_64-skylake-avx512-opensuse-leap-15.4`
-  `linux-x86_64-core-avx2-opensuse-leap-15.5`
-  `linux-x86_64-skylake-avx512-opensuse-leap-15.5`
-  `linux-x86_64-core-avx2-opensuse-leap-15.6`
-  `linux-x86_64-skylake-avx512-opensuse-leap-15.6`
-  `linux-x86_64-core-avx2-opensuse-tumbleweed`
-  `linux-x86_64-skylake-avx512-opensuse-tumbleweed`
-  `linux-x86_64-core-avx2-archlinux`
-  `linux-x86_64-skylake-avx512-archlinux`
-  `linux-x86_64-core-avx2-redhat-6`
-  `linux-x86_64-skylake-avx512-redhat-6`
-  `linux-x86_64-core-avx2-redhat-7`
-  `linux-x86_64-skylake-avx512-redhat-7`
-  `linux-x86_64-core-avx2-redhat-8`
-  `linux-x86_64-skylake-avx512-redhat-8`
-  `linux-x86_64-core-avx2-redhat-9`
-  `linux-x86_64-skylake-avx512-redhat-9`
-  `linux-x86_64-core-avx2-debian-12`
-  `linux-x86_64-skylake-avx512-debian-12`
-  `linux-x86_64-core-avx2-debian-11`
-  `linux-x86_64-skylake-avx512-debian-11`
-  `linux-x86_64-core-avx2-debian-10`
-  `linux-x86_64-skylake-avx512-debian-10`
-  `linux-x86_64-core-avx2-ubuntu-24`
-  `linux-x86_64-skylake-avx512-ubuntu-24`
-  `linux-x86_64-core-avx2-ubuntu-22`
-  `linux-x86_64-skylake-avx512-ubuntu-22`
-  `linux-x86_64-core-avx2-ubuntu-20`
-  `linux-x86_64-skylake-avx512-ubuntu-20`
-  `linux-x86_64-core-avx2-fedora-35`
-  `linux-x86_64-skylake-avx512-fedora-35`
-  `linux-x86_64-core-avx2-fedora-36`
-  `linux-x86_64-skylake-avx512-fedora-36`
-  `linux-x86_64-core-avx2-fedora-37`
-  `linux-x86_64-skylake-avx512-fedora-37`
-  `linux-x86_64-core-avx2-fedora-38`
-  `linux-x86_64-skylake-avx512-fedora-38`
-  `linux-x86_64-core-avx2-fedora-39`
-  `linux-x86_64-skylake-avx512-fedora-39`
-  `linux-x86_64-core-avx2-fedora-40`
-  `linux-x86_64-skylake-avx512-fedora-40`
-  `linux-aarch64`
-  `linux-aarch64-ubuntu-22`
-  `linux-aarch64-ubuntu-24`
-  `linux-aarch64-redhat-9`

Presets predetermine for your system what the cmake options are. For example, the paths to compilers and libraries are given. Also, the CXX flags determining the optimization level and system architecture are preset. Pick your preset, depending on your system:
- OS base (linux / mac / windows)
- Architecture (x86_64 for intel, aarch64 for arm or apple-silicon)
- Instruction set (skylake avx512 offers improved performance only on intel chips)
- Linux distribution


Installation on recent `macOS` computers is then accomplished for
example by

    cmake -B build --preset macos-apple-silicon
    ninja -C build install

or on `Linux`

    cmake -B build --preset linux-x86_64-core-avx2-ubuntu-22
    ninja -C build install

afterwhich the unit tests can be run

    ctest --test-dir build/tests --verbose

Packages can be created with

    ninja -C build package

`Doxygen` code documentation can be created with

    ninja -C build documentation

### Installing required packages on Ubuntu-24 or higher

    apt-get install -y git ca-certificates cmake ninja-build
    apt-get install -y llvm lld clang clang-tools clang-tidy
    apt-get install -y libc++-dev libc++abi-dev libomp-dev libclang-rt-dev
    apt-get install -y python3 pybind11-dev python3-pybind11 python3-dev
    apt-get install -y liblapack64-dev libblas64-dev
    apt-get install -y libhdf5-cpp-103-1t64 libhdf5-dev python3-h5py python3-tables

### Installing required packages on Fedora-40 or higher

    dnf install -y wget git rpm-build
    dnf install -y llvm lld cmake clang clang-tools-extra ninja-build
    dnf install -y libcxx libcxxabi libcxx-devel libcxxabi-devel 
    dnf install -y libomp-devel libcxx-static libcxxabi-static
    dnf install -y lapack-devel lapack64 blas64
    dnf install -y python3 python3-devel python3-pybind11
    dnf install -y pybind11-devel
    dnf install -y hdf5 hdf5-devel hdf5-static python3-h5py python3-tables

### Installing required packages on `macOS` using `Homebrew`

    dnf install -y wget git rpm-build
    dnf install -y llvm lld cmake clang clang-tools-extra ninja-build
    dnf install -y libcxx libcxxabi libcxx-devel libcxxabi-devel
    dnf install -y libomp-devel libcxx-static libcxxabi-static
    dnf install -y lapack-devel lapack64 blas64
    dnf install -y python3 python3-devel python3-pybind11
    dnf install -y pybind11-devel hdf5