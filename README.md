RASPA3
======

This software is a general purpose classical simulation package. It has been developed at
the University of Amsterdam (Amsterdam, The Netherlands) during 2022/2024 in active collaboration
with Eindhoven University of Technology (Eindhoven, Netherlands), Delft University of
Technology (Delft, The Netherlands), and Northwestern University (Evanston, USA).

Authors
=======
Drs. Youri Ran, University of Amsterdam<br>
Drs. Shrinjay Shoarma, Delft University of Technology<br>
Drs. Zhao Li, Northwestern University<br>
Dr. David Dubbeldam, University of Amsterdam<br>
Prof. Sofia Calero,  Eindhoven University of Technology<br>
Prof. Thijs Vlugt, Delft University of Technology<br>
Prof. Randall Q. Snurr, Northwestern University

Contributors
============
Alvaro Vazquez Mayagoitia, Argonne National Lab, contribution to openmp-implementation

Mac
=======
export SDKROOT=$(xcrun --sdk macosx --show-sdk-path)

Windows
=======
Visual studio 2022 project file included.
Showstopper: msvc does not support the multidimensional subscript operator

Linux
=====
Edit 'makefile'<br>
  CXX=clang++<br>
<br>
Set your compiler, e.g.<br>
  CXX=clang++-16<br>
or<br>
  CXX=clang++-17<br>

Ubuntu
======
sudo apt install build-essential git cmake<br>
wget -O - https://apt.llvm.org/llvm-snapshot.gpg.key | sudo apt-key add -<br>
wget -qO- https://apt.llvm.org/llvm-snapshot.gpg.key | sudo tee /etc/apt/trusted.gpg.d/apt.llvm.org.asc<br>
sudo apt install libc++-17-dev libc++abi-17-dev libomp-17-dev pybind11-dev python3-pybind11 liblapack-dev<br>

Redhat/Centos 8
===============
sudo dnf group install "Development Tools"<br>
sudo dnf install llvm-toolset<br>
sudo dnf install llvm-devel clang-devel<br>
sudo dnf install lldb python3-lit<br>
sudo dnf install cmake<br>
#sudo dnf install fmt-devel<br>
sudo dnf --enablerepo=devel install python3-devel python3.11-pybind11-devel<br>
sudo dnf --enablerepo=devel install pybind11-devel<br>
sudo dnf --enablerepo=devel install lapack-devel<br>
Note: no libc++-devel for Redhat-8<br>

Redhat 9
===============
sudo dnf group install "Development Tools"<br>
sudo dnf install llvm-toolset<br>
sudo dnf install llvm-devel clang-devel<br>
sudo dnf install lldb python3-lit<br>
sudo dnf install cmake<br>
sudo dnf install fmt-devel<br>
sudo dnf --enablerepo=devel install python3-devel python3-pybind11<br>
sudo dnf --enablerepo=devel install pybind11-devel<br>
sudo dnf --enablerepo=devel install lapack-devel<br>
Note: no libc++-devel for Redhat-9, the one of Fedora 38 works<br>
wget https://kojipkgs.fedoraproject.org//packages/libcxx/16.0.6/1.fc38/x86_64/libcxx-16.0.6-1.fc38.x86_64.rpm<br>
wget https://kojipkgs.fedoraproject.org//packages/libcxx/16.0.6/1.fc38/x86_64/libcxx-devel-16.0.6-1.fc38.x86_64.rpm<br>
wget https://kojipkgs.fedoraproject.org//packages/libcxx/16.0.6/1.fc38/x86_64/libcxx-static-16.0.6-1.fc38.x86_64.rpm<br>
wget https://kojipkgs.fedoraproject.org//packages/libcxx/16.0.6/1.fc38/x86_64/libcxxabi-16.0.6-1.fc38.x86_64.rpm<br>
wget https://kojipkgs.fedoraproject.org//packages/libcxx/16.0.6/1.fc38/x86_64/libcxxabi-devel-16.0.6-1.fc38.x86_64.rpm<br>
wget https://kojipkgs.fedoraproject.org//packages/libcxx/16.0.6/1.fc38/x86_64/libcxxabi-static-16.0.6-1.fc38.x86_64.rpm<br>
sudo dnf install libcxx*.rpm

Compilation
===========
for a parallel build using 16 cores use:<br>
make -j 16 release<br>
make -j 16 debug<br>
make -j 16 shared<br>
<br>
make<br>
defaults to<br>
make release<br>
<br>
to clean: make clean<br>

Running
=======
cd examples/basic/1_mc_methane_in_box<br>
./run

CMake (TODO)
============
export CXX=/usr/local/opt/llvm@17/bin/clang++
export CC=/usr/local/opt/llvm@17/bin/clang
export CXX=g++-13 
export CC=gcc-13 
mkdir build
cd build
cmake .. -GNinja
cmake --build . -v
cmake --build build --target test
cmake --build build --target docs

From Source
===========
mkdir ~/software
mkdir ~/source
cd ~/source
wget https://github.com/Kitware/CMake/releases/download/v3.29.1/cmake-3.29.1.tar.gz
tar -zxvf cmake-3.29.1.tar.gz
cd cmake-3.29.1
./bootstrap --prefix=${HOME}/software
make
make install
cd ~/source

wget https://github.com/llvm/llvm-project/archive/refs/tags/llvmorg-17.0.6.tar.gz
tar -zxvf llvmorg-17.0.6.tar.gz
cd llvm-project-llvmorg-17.0.6 
mkdir build
~/software/bin/cmake -G Ninja -S llvm -B build -DCMAKE_INSTALL_PREFIX=${HOME}/software -DCMAKE_BUILD_TYPE=Release  -DLLVM_ENABLE_PROJECTS="clang;openmp;lld" -DLLVM_ENABLE_RUNTIMES="compiler-rt;libcxx;libcxxabi;libunwind"  -DLLVM_TARGETS_TO_BUILD=X86 -DLLVM_RUNTIME_TARGETS="x86_64-unknown-linux-gnu" -DLLVM_BINUTILS_INCDIR=/usr/include
ninja -C build
ninja -C build install
~/software/bin/cmake -G Ninja -S runtimes -B build -DCMAKE_INSTALL_PREFIX=${HOME}/software -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=${HOME}/software/bin/clang -DCMAKE_CXX_COMPILER=${HOME}/software/bin/clang++ -DLLVM_ENABLE_RUNTIMES="libcxx;libcxxabi;libunwind"
rm -rf build
mkdir build
ninja -C build 
ninja -C build install

wget https://github.com/pybind/pybind11/archive/refs/tags/v2.12.0.tar.gz
tar -zxvf v2.12.0.tar.gz
cd pybind11-2.12.0
mkdir build
cmake3 -B build -DCMAKE_INSTALL_PREFIX=${HOME}/software .
make -j 16
make install

wget https://github.com/ninja-build/ninja/archive/refs/tags/v1.11.1.tar.gz
tar -zxvf v1.11.1.tar.gz 
cd ninja-1.11.1/

use run file:

#! /bin/sh -f
export RASPA_DIR=${HOME}/raspa3
export LD_LIBRARY_PATH=${HOME}/software/lib
../../../src/raspa3.exe


Compiling RASPA3
================
mkdir build
cmake -B build -GNinja -DCMAKE_INSTALL_PREFIX=${HOME}/raspa3 -DCMAKE_CXX_COMPILER=/usr/local/opt/llvm@17/bin/clang++ -DCMAKE_BUILD_WITH_INSTALL_RPATH=ON -DCMAKE_BUILD_TYPE=Release .
cmake --build build  (or: ninja -C build -v)
cmake --install build --config Release
ctest --test-dir build/tests --verbose
ctest --test-dir build/tests/raspakit-tests --verbose
