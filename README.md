RASPA3
======

This software is a general purpose classical simulation package. It has been developed at
the University of Amsterdam (Amsterdam, The Netherlands) during 2022 in active collaboration
with Eindhoven University of Technology (Eindhoven, Netherlands), Delft University of
Technology (Delft, The Netherlands), and Northwestern University (Evanston, USA).

Authors
=======
Dr. David Dubbeldam, University of Amsterdam<br>
Prof. Sofia Calero,  Eindhoven University of Technology<br>
Prof. Thijs Vlugt, Delft University of Technology<br>
Prof. Randall Q. Snurr, Northwestern University

Contributors
============
Zhao Li, Northwestern University, contribution to openmp-implementation
Alvaro Vazquez Mayagoitia, Argonne National Lab, contribution to openmp-implementation
Shrinjay, Delft University of Technology, contribution to the breakthrough implementation

Windows
=======
Visual studio 2022 project file included.

Ubuntu 22.04
============
sudo apt install build-essential git clang-14 libc++-dev libc++abi-dev gfortran liblapack-dev libfmt-dev libomp-dev

Compilation
===========
make -j 32

to clean: make clean

Running
=======
cd "raspa3-tests/1"
./run


