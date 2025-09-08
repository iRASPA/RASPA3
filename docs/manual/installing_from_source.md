# Installing `RASPA` from Source
\page installing_from_source Installing from Source


## 1. Building from Source

### 1.1 Prerequisites

| Package                | Min. version |
| ---------------------- | ------------ |
| clang/clang++          | 18           |
| cmake                  | 3.28         |
| ninja                  | 1.11         |
| python                 | 3.11         |
| pybind11               | 2.12         |
| BLAS + LAPACK          | 64‑bit ints  |
| HDF5                   | 1.12         |
| OpenMP · OpenCL        | latest       |

### 1.2 Build with Conda (recommended)

```bash
git clone https://github.com/raspa3/raspa3.git
cd raspa3
conda env create -f env.yml
conda activate raspa
cmake --preset=linux_conda  # or mac_conda / windows_conda_raspa3
ninja -C build
ninja -C build install      # optional
```

**Note**: remove the package ocl-icd from env.yml for Mac and Windows installations, as OpenCL is only available for Linux. For MacOS using the system OpenCL is recommended.

**Note**: it is often necessary to add the conda installed libraries raspa3 links against to your dynamic library path using `export LD_LIBRARY_PATH=${CONDA_PREFIX}/lib:$LD_LIBRARY_PATH`.

### 1.3 Linux examples

#### Ubuntu 24.04

```bash
sudo apt install -y git ca-certificates cmake ninja-build \
  llvm lld clang clang-tidy libc++-dev libc++abi-dev \
  libomp-dev libclang-rt-dev python3 python3-dev pybind11-dev \
  liblapack64-dev libblas64-dev

git clone https://github.com/raspa3/raspa3.git
cd raspa3
cmake -B build --preset=linux-x86_64-core-avx2-ubuntu-24
ninja -C build
ninja -C build install
```

#### Fedora 40

```bash
sudo dnf install -y git llvm lld cmake clang clang-tools-extra ninja-build \
  libomp-devel libcxx libcxx-devel libcxxabi libcxxabi-devel \
  lapack64 blas64 python3 python3-devel pybind11-devel

git clone https://github.com/raspa3/raspa3.git
cd raspa3
cmake -B build --preset=linux-x86_64-core-avx2-fedora-40
ninja -C build
ninja -C build install
```

### 1.4 macOS (Homebrew)

```zsh
brew install llvm lld libomp hdf5 libaec ninja cmake doxygen graphviz lapack pybind11
export PATH="$(brew --prefix llvm)/bin:$PATH"  # clang‑18
export HDF5_ROOT="$(brew --prefix hdf5)"
git clone https://github.com/raspa3/raspa3.git
cd raspa3
cmake -B build --preset=macos-apple-silicon  # or macos-x64-core-avx2
ninja -C build
ninja -C build install
```

---

## 2. CMake Presets

Run `cmake --list-presets` to see every preset.
Common shortcuts:

| Preset                                 | Target                |
| -------------------------------------- | --------------------- |
| linux\_conda                           | Generic Linux + Conda |
| linux-x86\_64-core-avx2-ubuntu-24      | Ubuntu 24, AVX2       |
| linux-x86\_64-skylake-avx512-fedora-40 | Fedora 40, AVX‑512    |
| mac\_conda                             | macOS + Conda         |
| macos-apple-silicon                    | macOS Apple Silicon   |
| windows\_conda\_raspa3                 | Windows + Conda       |

---

