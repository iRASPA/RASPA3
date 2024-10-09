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

-   `macos-x64`
-   `macos-x64-debug`
-   `macos-apple-silicon`
-   `macos-apple-silicon-debug`
-   `windows-x64`
-   `windows-arm64`
-   `linux`
-   `linux-opensuse-leap-15.2`
-   `linux-opensuse-leap-15.3`
-   `linux-opensuse-leap-15.4`
-   `linux-opensuse-leap-15.5`
-   `linux-opensuse-tumbleweed`
-   `linux-archlinux`
-   `linux-redhat-6`
-   `linux-redhat-7`
-   `linux-redhat-8`
-   `linux-redhat-9`
-   `linux-debian-12`
-   `linux-debian-11`
-   `linux-debian-10`
-   `linux-ubuntu-24`
-   `linux-ubuntu-22`
-   `linux-ubuntu-20`
-   `linux-fedora-35`
-   `linux-fedora-36`
-   `linux-fedora-37`
-   `linux-fedora-38`
-   `linux-fedora-39`
-   `linux-fedora-40`

Installation on recent `macOS` computers is then accomplished for
example by

    cmake -B build --preset macos-apple-silicon
    ninja -C build install

or on `Linux`

    cmake -B build --preset linux-ubuntu-22
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