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

Windows
=======
Visual studio 2022 project file included.

Ubuntu
======
sudo apt install build-essential git cmake
wget -O - https://apt.llvm.org/llvm-snapshot.gpg.key | sudo apt-key add -
wget -qO- https://apt.llvm.org/llvm-snapshot.gpg.key | sudo tee /etc/apt/trusted.gpg.d/apt.llvm.org.asc
sudo apt install libc++-17-dev libc++abi-17-dev libomp-17-dev pybind11-dev python3-pybind11 liblapack-dev

Redhat/Centos 8
===============
yum groupinstall 'Development Tools'
yum module install llvm-toolset
yum install llvm-devel clang-devel
yum install lldb python3-lit
yum install cmake
yum install lapack
Showstopper: no libc++ for Redhat-8

Redhat 9
===============
sudo dnf group install "Development Tools"
sudo dnf install llvm-toolset
sudo dnf install llvm-devel clang-devel
sudo dnf install lldb python3-lit
sudo dnf install cmake
sudo dnf --enablerepo=devel install python3-pybind11
sudo dnf --enablerepo=devel install pybind11-devel
sudo dnf --enablerepo=devel install lapack-devel
Note: no libc++-devel for Redhat-9
wget https://kojipkgs.fedoraproject.org//packages/libcxx/16.0.6/1.fc38/x86_64/libcxx-16.0.6-1.fc38.x86_64.rpm
wget https://kojipkgs.fedoraproject.org//packages/libcxx/16.0.6/1.fc38/x86_64/libcxx-devel-16.0.6-1.fc38.x86_64.rpm
wget https://kojipkgs.fedoraproject.org//packages/libcxx/16.0.6/1.fc38/x86_64/libcxx-static-16.0.6-1.fc38.x86_64.rpm
wget https://kojipkgs.fedoraproject.org//packages/libcxx/16.0.6/1.fc38/x86_64/libcxxabi-16.0.6-1.fc38.x86_64.rpm
wget https://kojipkgs.fedoraproject.org//packages/libcxx/16.0.6/1.fc38/x86_64/libcxxabi-devel-16.0.6-1.fc38.x86_64.rpm
wget https://kojipkgs.fedoraproject.org//packages/libcxx/16.0.6/1.fc38/x86_64/libcxxabi-static-16.0.6-1.fc38.x86_64.rpm
sudo dnf install libcxx*.rpm

Compilation
===========
make -j 32
or
make -j 32 debug

to clean: make clean

Running
=======
cd "raspa3-tests/1"
./run


