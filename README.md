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

Windows
=======
Visual studio 2022 project file included.

Ubuntu 22.04
============
sudo apt install build-essential git clang-14 libc++-dev libc++abi-dev gfortran liblapack-dev libfmt-dev

Intel Threading Building Blocks
===============================
sudo apt install libtbb-dev
link with: -ltbb

use as:
#include <execution>
#include <algorithm>
std::sort(std::execution::par_unseq, input.begin(), input.end());


Modify libfmt
=============
edit:
/usr/include/fmt/core.h
change line 1358 from: return to_string_view(val)
to: return basic_string_view<char>(val);

(This should be solved in 'clang-15' when 'std::print' is fully implemented.

Compilation
===========
make -f makefile-llvm all

to clean: make -f makefile-llvm allclean

Running
=======
cd "raspa3-tests/1"
./run


