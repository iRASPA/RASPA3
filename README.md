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
<img alt="GitHub Actions Workflow Status" src="https://img.shields.io/github/actions/workflow/status/iRASPA/raspa3/github-actions-create-release.yml">
  </a>
  <a href="https://github.com/iRASPA/raspa3/issues"><img alt="GitHub Issues or Pull Requests" src="https://img.shields.io/github/issues/iRASPA/raspa3">
</a>
</p>

<p align="center">
  <a href="#authors">Authors</a> •
  <a href="#contributors">Contributors</a> •
  <a href="#installation-guide">Installation Guide</a> •
  <a href="#how-to-use">How To Use</a> 
</p>

# Authors

* Drs. Youri Ran, University of Amsterdam
* Drs. Shrinjay Shoarma, Delft University of Technology
* Drs. Zhao Li, Northwestern University
* Dr. David Dubbeldam, University of Amsterdam
* Prof. Sofia Calero, Eindhoven University of Technology
* Prof. Thijs Vlugt, Delft University of Technology
* Prof. Randall Q. Snurr, Northwestern University

# Contributors

* Alvaro Vazquez Mayagoitia, Argonne National Lab, contribution to openmp-implementation

# Installation Guide

<!-- TOC -->
* [Installation Guide](#installation-guide)
    * [Linux](#linux)
      * [Ubuntu 24 (fully supported)](#ubuntu-24-fully-supported)
      * [Ubuntu 22.04 and LLVM-17](#ubuntu-2204-and-llvm-17)
      * [Redhat/Centos 8](#redhatcentos-8)
      * [Redhat 9](#redhat-9)
      * [carbonOS](#carbonos)
      * [Make rpm](#make-rpm)
    * [Windows](#windows)
    * [macOS](#macos)
    * [Python](#python)
<!-- TOC -->

### Linux

#### Ubuntu 24 (fully supported)

1. First, you need to obtain the latest LLVM toochain. Simply use the following script to easily install the LLVM into
   the system.

* Download the script

   ```bash
   wget https://apt.llvm.org/llvm.sh
   ```

* Make it executable

   ```bash
   chmod +x llvm.sh
   ```

* Install LLVM

   ```bash
   sudo ./llvm.sh all
   ```

* Cleaning

   ```bash
   rm -f ./llvm.sh
   ```

2. Install dependencies.

* Requried

   ```bash
   sudo apt install pybind11-dev python3-pybind11 doxygen graphviz git cmake ninja-build
   ```

* For BLAS and LAPACK (optional but **RECOMMEND**）

   ```bash
   sudo apt install libblas64-dev liblapack64-dev
   ```

3. Clone the project

   ```bash
   git clone https://github.com/iRASPA/raspa3.git
   ```
4. Compile

* Chaning working directory

   ```bash
   cd raspa3
   ```

* Configure Makefiles

   ```bash
   cmake -B build -DCMAKE_CXX_COMPILER=clang++-18 -DCMAKE_C_COMPILER=clang-18 -DCMAKE_CXX_FLAGS="-Wno-unused-but-set-variable" -DCMAKE_CXX_STANDARD_INCLUDE_DIRECTORIES="/usr/lib/llvm-18/include/c++/v1;/usr/lib/llvm-18/include/x86_64-pc-linux-gnu/c++/v1" -DCMAKE_LIBRARY_PATH="/usr/lib/llvm-18/lib" --preset linux-ubuntu-24
   ```

* Build

   ```bash
   ninja -C build package
   ```

5. Install

   ```bash
   sudo dpkg -i build/raspa_3.0.0_amd64.deb
   ```

#### Ubuntu 22.04 and LLVM-17

```bash
sudo apt install build-essential git cmake
sudo add-apt-repository 'deb [http://apt.llvm.org/jammy/](http://apt.llvm.org/jammy/) llvm-toolchain-jammy-17 main'
sudo add-apt-repository 'deb-src [http://apt.llvm.org/jammy/](http://apt.llvm.org/jammy/) llvm-toolchain-jammy-17 main'
sudo apt install libllvm-17-ocaml-dev libllvm17 llvm-17 llvm-17-dev llvm-17-runtime
sudo apt install libc++-17-dev libc++abi-17-dev libomp-17-dev pybind11-dev python3-pybind11 liblapack-dev clang-tools-17
cd src
make -f makefile-manual
make -f makefile-manual install
```

#### Redhat/Centos 8

```bash
sudo dnf group install "Development Tools"
sudo dnf install llvm-toolset
sudo dnf install llvm-devel clang-devel
sudo dnf install lldb python3-lit
sudo dnf install cmake
#sudo dnf install fmt-devel
sudo dnf --enablerepo=devel install python3-devel python3.11-pybind11-devel
sudo dnf --enablerepo=devel install pybind11-devel
sudo dnf --enablerepo=devel install lapack-devel
```

> Note: no libc++-devel for Redhat-8

#### Redhat 9

```bash
sudo dnf group install "Development Tools"
sudo dnf install llvm-toolset
sudo dnf install llvm-devel clang-devel
sudo dnf install lldb python3-lit
sudo dnf install cmake
sudo dnf install fmt-devel
sudo dnf --enablerepo=devel install python3-devel python3-pybind11
sudo dnf --enablerepo=devel install pybind11-devel
sudo dnf --enablerepo=devel install lapack-devel
sudo dnf install [https://kojipkgs.fedoraproject.org//packages/libcxx/16.0.6/1.fc38/x86_64/libcxx-16.0.6-1.fc38.x86_64.rpm](https://kojipkgs.fedoraproject.org//packages/libcxx/16.0.6/1.fc38/x86_64/libcxx-16.0.6-1.fc38.x86_64.rpm)
sudo dnf install [https://kojipkgs.fedoraproject.org//packages/libcxx/16.0.6/1.fc38/x86_64/libcxx-devel-16.0.6-1.fc38.x86_64.rpm](https://kojipkgs.fedoraproject.org//packages/libcxx/16.0.6/1.fc38/x86_64/libcxx-devel-16.0.6-1.fc38.x86_64.rpm)
sudo dnf install [https://kojipkgs.fedoraproject.org//packages/libcxx/16.0.6/1.fc38/x86_64/libcxx-static-16.0.6-1.fc38.x86_64.rpm](https://kojipkgs.fedoraproject.org//packages/libcxx/16.0.6/1.fc38/x86_64/libcxx-static-16.0.6-1.fc38.x86_64.rpm)
sudo dnf install [https://kojipkgs.fedoraproject.org//packages/libcxx/16.0.6/1.fc38/x86_64/libcxxabi-16.0.6-1.fc38.x86_64.rpm](https://kojipkgs.fedoraproject.org//packages/libcxx/16.0.6/1.fc38/x86_64/libcxxabi-16.0.6-1.fc38.x86_64.rpm)
sudo dnf install [https://kojipkgs.fedoraproject.org//packages/libcxx/16.0.6/1.fc38/x86_64/libcxxabi-devel-16.0.6-1.fc38.x86_64.rpm](https://kojipkgs.fedoraproject.org//packages/libcxx/16.0.6/1.fc38/x86_64/libcxxabi-devel-16.0.6-1.fc38.x86_64.rpm)
sudo dnf install [https://kojipkgs.fedoraproject.org//packages/libcxx/16.0.6/1.fc38/x86_64/libcxxabi-static-16.0.6-1.fc38.x86_64.rpm](https://kojipkgs.fedoraproject.org//packages/libcxx/16.0.6/1.fc38/x86_64/libcxxabi-static-16.0.6-1.fc38.x86_64.rpm)
```

> Note: no libc++-devel for Redhat-9, the one of Fedora 38 works

#### carbonOS

```bash
export PATH=${HOME}/software/bin/:${PATH}
mkdir build
${HOME}/software/bin/cmake -B build -GNinja -DCMAKE_INSTALL_PREFIX=${HOME}/raspa3 -DCMAKE_CC_COMPILER=${HOME}/software/bin/clang -DCMAKE_CXX_COMPILER=${HOME}/software/bin/clang++ -DCMAKE_BUILD_WITH_INSTALL_RPATH=ON -DCMAKE_MAKE_PROGRAM=${HOME}/software/bin/ninja -DCMAKE_BUILD_TYPE=Release .
${HOME}/software/bin/cmake --build build (or: ninja -C build -v)
${HOME}/software/bin/cmake --install build --config Release
${HOME}/software/bin/ctest --test-dir build/tests --verbose
${HOME}/software/bin/ctest --test-dir build/tests/raspakit-tests --verbose
```

#### Make rpm

```bash
cmake -B build -GNinja -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_BUILD_WITH_INSTALL_RPATH=ON -DCMAKE_BUILD_TYPE=Release -DBUILD_RPM_PACKAGE=ON . ninja -C build ninja -C build package rpm -qlp /home/dubbelda/source/raspa3/build/raspa3-3.0.0-1.el7.centos.x86_64.rpm
```

### Windows

Visual studio 2022 project file included. Showstopper: msvc does not support the multidimensional subscript operator

### macOS

```bash
export SDKROOT=$(xcrun --sdk macosx --show-sdk-path)
```

### Python

```bash
export PYTHONPATH=$PYTHONPATH:/usr/share/raspa3/lib:/usr/local/share/raspa3/lib 
cd python 
python3 script.py
```

# How To Use

1. Set environment variables

* For bash

   ```bash
   echo 'export RASPA_DIR=/usr/share/raspa3' >> ~/.bashrc
   ```

* For zsh

   ```bash
   echo 'export RASPA_DIR=/usr/share/raspa3' >> ~/.zshrc
   ```

2. Run

   ```bash
   raspa
   ```


