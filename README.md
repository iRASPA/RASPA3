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

LLVM general
============
wget https://apt.llvm.org/llvm.sh<br>
chmod u+x llvm.sh<br>
sudo ./llvm.sh 17<br>

Ubuntu 22.04 and LLVM-17
========================
sudo apt install build-essential git cmake<br>
sudo add-apt-repository 'deb http://apt.llvm.org/jammy/ llvm-toolchain-jammy-17 main'<br>
sudo add-apt-repository 'deb-src http://apt.llvm.org/jammy/ llvm-toolchain-jammy-17 main'<br>
sudo apt install libllvm-17-ocaml-dev libllvm17 llvm-17 llvm-17-dev llvm-17-runtime<br>
sudo apt install libc++-17-dev libc++abi-17-dev libomp-17-dev pybind11-dev python3-pybind11 liblapack-dev clang-tools-17<br>
cd src<br>
make -f makefile-manual<br>
make -f makefile-manual install<br>

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
sudo dnf install https://kojipkgs.fedoraproject.org//packages/libcxx/16.0.6/1.fc38/x86_64/libcxx-16.0.6-1.fc38.x86_64.rpm<br>
sudo dnf install https://kojipkgs.fedoraproject.org//packages/libcxx/16.0.6/1.fc38/x86_64/libcxx-devel-16.0.6-1.fc38.x86_64.rpm<br>
sudo dnf install https://kojipkgs.fedoraproject.org//packages/libcxx/16.0.6/1.fc38/x86_64/libcxx-static-16.0.6-1.fc38.x86_64.rpm<br>
sudo dnf install https://kojipkgs.fedoraproject.org//packages/libcxx/16.0.6/1.fc38/x86_64/libcxxabi-16.0.6-1.fc38.x86_64.rpm<br>
sudo dnf install https://kojipkgs.fedoraproject.org//packages/libcxx/16.0.6/1.fc38/x86_64/libcxxabi-devel-16.0.6-1.fc38.x86_64.rpm<br>
sudo dnf install https://kojipkgs.fedoraproject.org//packages/libcxx/16.0.6/1.fc38/x86_64/libcxxabi-static-16.0.6-1.fc38.x86_64.rpm<br>

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

Ubuntu 24 (fully supported)
===========================
sudo apt install build-essential git cmake ninja-build<br>
sudo apt install llvm-17 clang-17 clang-tools-17 clang-tidy-17 libc++-17-dev libc++abi-17-dev libomp-17-dev<br>
sudo apt pybind11-dev python3-pybind11<br>
sudo apt liblapack-dev<br>
sudo apt doxygen graphviz<br>
cmake -B build -GNinja -DCMAKE_INSTALL_PREFIX=${HOME}/raspa3 -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_BUILD_WITH_INSTALL_RPATH=ON -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON .<br>
cmake --build build  (or: ninja -C build -v)<br>
cmake --install build --config Release<br>
ctest --test-dir build/tests --verbose<br>
ctest --test-dir build/tests/raspakit-tests --verbose<br>

CMake (TODO)
============
export CXX=/usr/local/opt/llvm@17/bin/clang++<br>
export CC=/usr/local/opt/llvm@17/bin/clang<br>
export CXX=g++-13<br>
export CC=gcc-13 <br>
mkdir build<br>
cd build<br>
cmake .. -GNinja<br>
cmake --build . -v<br>
cmake --build build --target test<br>
cmake --build build --target docs<br>

From Source
===========
mkdir ~/software<br>
mkdir ~/source<br>
cd ~/source<br>
wget https://github.com/Kitware/CMake/releases/download/v3.29.1/cmake-3.29.1.tar.gz<br>
tar -zxvf cmake-3.29.1.tar.gz<br>
cd cmake-3.29.1<br>
./bootstrap --prefix=${HOME}/software<br>
make<br>
make install<br>
cd ~/source<br>
<br>
wget https://github.com/llvm/llvm-project/archive/refs/tags/llvmorg-17.0.6.tar.gz<br>
tar -zxvf llvmorg-17.0.6.tar.gz<br>
cd llvm-project-llvmorg-17.0.6 <br>
mkdir build<br>
~/software/bin/cmake -G Ninja -S llvm -B build -DCMAKE_INSTALL_PREFIX=${HOME}/software -DCMAKE_BUILD_TYPE=Release  -DLLVM_ENABLE_PROJECTS="clang;openmp;lld" -DLLVM_ENABLE_RUNTIMES="compiler-rt;libcxx;libcxxabi;libunwind"  -DLLVM_TARGETS_TO_BUILD=X86 -DLLVM_RUNTIME_TARGETS="x86_64-unknown-linux-gnu" -DLLVM_BINUTILS_INCDIR=/usr/include<br>
ninja -C build<br>
ninja -C build install<br>
~/software/bin/cmake -G Ninja -S runtimes -B build -DCMAKE_INSTALL_PREFIX=${HOME}/software -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=${HOME}/software/bin/clang -DCMAKE_CXX_COMPILER=${HOME}/software/bin/clang++ -DLLVM_ENABLE_RUNTIMES="libcxx;libcxxabi;libunwind"<br>
rm -rf build<br>
mkdir build<br>
ninja -C build<br>
ninja -C build install<br>
<br>
wget https://github.com/pybind/pybind11/archive/refs/tags/v2.12.0.tar.gz<br>
tar -zxvf v2.12.0.tar.gz<br>
cd pybind11-2.12.0<br>
mkdir build<br>
cmake3 -B build -DCMAKE_INSTALL_PREFIX=${HOME}/software .<br>
make -j 16<br>
make install<br>

wget https://github.com/ninja-build/ninja/archive/refs/tags/v1.11.1.tar.gz<br>
tar -zxvf v1.11.1.tar.gz<br>
cd ninja-1.11.1/<br>
<br>
use run file:<br>
<br>
#! /bin/sh -f<br>
export RASPA_DIR=${HOME}/raspa3<br>
export LD_LIBRARY_PATH=${HOME}/software/lib<br>
../../../src/raspa3.exe<br>


Compiling RASPA3
================
mkdir build<br>
cmake -B build -GNinja -DCMAKE_INSTALL_PREFIX=${HOME}/raspa3 -DCMAKE_CXX_COMPILER=/usr/local/opt/llvm@17/bin/clang++ -DCMAKE_BUILD_WITH_INSTALL_RPATH=ON -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON .<br>
cmake --build build  (or: ninja -C build -v)<br>
cmake --install build --config Release<br>
ctest --test-dir build/tests --verbose<br>
ctest --test-dir build/tests/raspakit-tests --verbose<br>

Carbon
======
export PATH=${HOME}/software/bin/:${PATH}<br>
mkdir build<br>
${HOME}/software/bin/cmake -B build -GNinja -DCMAKE_INSTALL_PREFIX=${HOME}/raspa3 -DCMAKE_CC_COMPILER=${HOME}/software/bin/clang  -DCMAKE_CXX_COMPILER=${HOME}/software/bin/clang++ -DCMAKE_BUILD_WITH_INSTALL_RPATH=ON -DCMAKE_MAKE_PROGRAM=${HOME}/software/bin/ninja -DCMAKE_BUILD_TYPE=Release .<br>
${HOME}/software/bin/cmake --build build  (or: ninja -C build -v)<br>
${HOME}/software/bin/cmake --install build --config Release<br>
${HOME}/software/bin/ctest --test-dir build/tests --verbose<br>
${HOME}/software/bin/ctest --test-dir build/tests/raspakit-tests --verbose<br>

Make rpm
========
cmake -B build -GNinja  -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_BUILD_WITH_INSTALL_RPATH=ON -DCMAKE_BUILD_TYPE=Release -DBUILD_PACKAGE=ON .
ninja -C build
ninja -C build package
rpm -qlp /home/dubbelda/source/raspa3/build/raspa3-3.0.0-1.el7.centos.x86_64.rpm<br>

