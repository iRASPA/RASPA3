# RASPA 3 – Installation Guide

RASPA 3 is a modern Monte‑Carlo simulation package that uses C++23.
Pre‑built binaries are **statically linked** and therefore do not require any external runtime libraries.
If you only want to *use* RASPA 3, download a pre‑built package.
If you want to *develop* or *debug* RASPA 3, follow the build‑from‑source instructions.

> **Tested toolchain:** LLVM 18 / Clang 18

---

## 1. Quick Start (Pre‑built Packages)

All official releases are published on the **GitHub “Releases” tab**.
Choose the file that matches your operating system, CPU architecture and AVX capability.

| Platform        | File suffix | Example                                           |
| --------------- | ----------- | ------------------------------------------------- |
| Conda           | –           | install with `conda install raspa3 raspalib`      |
| Debian / Ubuntu | `.deb`      | `raspa‑{VERSION}‑{OS}‑{ARCH}‑{AVX}.deb`      |
| RPM‑based Linux | `.rpm`      | `raspa‑{VERSION}‑{OS}‑{ARCH}‑{AVX}.rpm` |
| macOS           | `.pkg`      | `raspa‑{VERSION}‑mac‑{ARCH}.pkg`                      |
| Windows         | `.exe`      | `raspa‑{VERSION}‑windows‑{ARCH}‑{AVX}.exe`       |

### 1.1 Conda (all platforms)

```bash
conda install -c conda-forge raspa3 raspalib
```

### 1.2 Debian / Ubuntu (.deb)

```bash
apt install ./raspa3-<VERSION>-ubuntu<20|22|24>-<ARCH>-<AVX>.deb
```

### 1.3 RPM‑based Linux (.rpm)

```bash
rpm -i raspa3-<VERSION>-<OS>-<ARCH>-<AVX>.rpm
```

*To unpack without installing:* `rpm2cpio raspa3-*.rpm | cpio -idmv`

### 1.4 macOS (.pkg)

Pick either the arm64 version (Apple silicon) or the x86_64 version (intel chips). Double‑click the `.pkg` and follow the wizard.

### 1.5 Windows (.exe)

Run the `.exe`; the wizard adds `%RASPA3_HOME%\bin` to your `%PATH%`.

### 1.6 Choosing AVX2 vs AVX‑512

| Binary         | Use when                                             |
| -------------- | ---------------------------------------------------- |
| core‑avx2      | Intel ≥ Haswell (2013) or any AMD ≥ Excavator (2015) |
| skylake‑avx512 | Intel CPUs with AVX‑512 (≈ 25 % faster)              |

---

## 2. Building from Source

### 2.1 Prerequisites

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

### 2.2 Build with Conda (recommended)

```bash
git clone https://github.com/raspa3/raspa3.git
cd raspa3
conda env create -f env.yml
conda activate raspa
cmake --preset=linux_conda  # or mac_conda / windows_conda_raspa3
ninja -C build
ninja -C build install      # optional
```

**Note**: remove the package ocl-icd from env.yml for Mac and Windows installations, as OpenCL is only available for Linux. For Max using the system OpenCL is recommended.

### 2.3 Linux examples

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

### 2.4 macOS (Homebrew)

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

## 3. CMake Presets

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

Need help?  Open an issue on GitHub or ask on the discussion forum.
