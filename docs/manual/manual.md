# Introduction
\page manual_main Manual

## Design philosophy

`RASPA3` is a molecular simulation code for computing adsorption and
diffusion in nanoporous materials, and thermodynamic and transport
properties of fluids. It implements force field based classical Monte
Carlo/Molecular Dynamics in various ensembles. `RASPA3` can be used for
the simulation of molecules in gases, fluids, zeolites,
aluminosilicates, metal-organic frameworks, and carbon nanotubes. One of
the main applications of `RASPA3` is to compute adsorption isotherms and
to understand the atomic-level mechanism of adsorption and separations.

`RASPA3` is redesigned and rewritten from the ground up in `C++23`,
based on the following ideas:

-   Composition and value semantics\
    Major improvements in efficiency result from the implementation of
    `C++` code based on composition and value semantics (in contrast to
    inheritance and reference semantics). Composition is a structural
    approach in which complex objects are built from simple building
    blocks. Value semantics avoid sharing mutable state and upholds the
    independence of values to support local reasoning. References are
    only used implicitly, at function boundaries, and are never stored.
    This style of coding has many advantages. Avoiding complex
    inheritance hierarchies and reducing the need for virtual functions
    lead to major code simplifications. It encourages code clarity and
    local reasoning. With composition, components are directly included
    in an object, removing the necessity for pointers, manual memory
    allocations, and indirections. This leads to straightforward memory
    patterns and access safety. Lifetimes of objects are coupled with
    composition, eliminating the need for explicit memory management.
    These lead to better compiler optimization and therefore performance
    of the code.

-   Correctness and accuracy\
    For all the techniques and algorithms available in `RASPA3` we have
    implemented the 'best' ones available in literature. For example,
    `RASPA3` uses Configurational-Bias Monte-Carlo, it uses the Ewald
    summation for electrostatics, molecular dynamics is based on
    'symplectic' integrators, all Monte-Carlo moves obey detailed
    balance etc. Unit tests test the smallest functional units of code
    and prevent developers breaking existing code when refactoring or
    adding new code. The unit tests in `RASPA3` are arranged
    hierarchically. Atoms are placed at several distances and the
    Lennard-Jones potentials are tested and compared with the
    analytically computed energy. Then molecules are placed within a
    zeolite framework at predetermined positions and compared to the
    energies computed with `RASPA2`. The various energy routines used in
    biased-sampling are tested by comparing to the general routines. The
    various routines for the gradients are tested by comparing the
    computed gradients to the values computed by finite difference
    schemes based on the energy. Likewise, the strain-derivative tensor
    (related to the pressure) is tested by comparing to the values
    computed by finite difference schemes based on the energy of a
    strained cell.

-   Input made easy\
    The input of `RASPA3` has been changed to `JSON` format. `JSON`
    stands for `JavaScript Object Notation`. It is a lightweight
    data-interchange format in plain text written in JavaScript object
    notation. It is language independent. The requirements for the input
    files is kept as minimal as possible. Only for more advanced options
    extra commands in the input file are needed. Also the format of the
    input is straightforward. Default settings are usually the best
    ones. Fugacity coefficients and excess adsorption are automatically
    computed.

-   Output made easy\
    Similarly, the `JSON` format is now also used for the output. Where
    previously all initial data (e.g. hardware info, unit conversion
    factors) and statistics (e.g. `CPU` timings, `MC` move statistics)
    were written to a text file, these are now written as a nested
    dictionary to a `JSON` file.

-   Integrated simulation environment\
    The `RASPA3` `C++` simulation engine is made available to `Python`
    via a `Pybind11` interface. For use in `Python`, `RASPA3` is built
    as a shared library, allowing its functions to be used by `Python`
    users and enabling seamless interactions with the simulation
    routines. Through the API, the library allows for invocation of
    `RASPA3`'s simulation routines from `Python` scripts, calling the
    same simulation routines as via the `JSON` input. This execution
    directly from `Python` enables the ability to prototype simulation
    settings and incorporate `RASPA3` into existing workflows. Extension
    and modification of the code is relatively straightforward.

## Output from `RASPA`

`RASPA` generates output from the simulation. Some data is just
information on the status, while other data are written because you
specifically asked the program to compute it for you. The output is
written to be used with other programs like:

-   `gnuplot`

-   `iRASPA`

-   `VMD`

## Citing `RASPA`

If you are using `RASPA` and would like to cite it in your journal
articles or book-chapters, then for `RASPA`:

> Y.A. Ran, S. Sharma, S.R.G. Balestra, Z. Li, S. Calero, T.J.H Vlugt,
> R.Q. Snurr, and D. Dubbeldam RASPA3: A Monte Carlo Code for Computing
> Adsorption and Diffusion in Nanoporous Materials and Thermodynamics
> Properties of Fluids, 
> J. Chem. Phys. 2024, 161, 114106
> <https://doi.org/10.1063/5.0226249>

> D. Dubbeldam, S. Calero, D.E. Ellis, and R.Q. Snurr, RASPA: Molecular
> Simulation Software for Adsorption and Diffusion in Flexible
> Nanoporous Materials,
> 2015
> <http://dx.doi.org/10.1080/08927022.2015.1010082>

For the inner workings of Monte Carlo codes:

> D. Dubbeldam, A. Torres-Knoop, and K.S. Walton, On the Inner Workings
> of Monte Carlo Codes, <http://dx.doi.org/10.1080/08927022.2013.819102>
> , 39(14-15), 1253-1292, 2013.

For the description of Molecular Dynamics and diffusion:

> D. Dubbeldam and R.Q. Snurr, Recent Developments in the Molecular
> Modeling of Diffusion in Nanoporous Materials,
> <http://dx.doi.org/10.1080/08927020601156418>, , 33(4-5), 305-325,
> 2007.

For the description of the implementation of force fields:

> D. Dubbeldam, K.S. Walton, T.J.H. Vlugt, and S. Calero, Design,
> Parameterization, and Implementation of Atomic Force Fields for
> Adsorption in Nanoporous Materials,
> <https://doi.org/10.1002/adts.201900135>, , 2(11), 1900135, 2019.

For the implementation of the grid interpolation methods:

> Y.A. Ran, J. Tapia, S. Sharma, P. Bai, S. Calero, T.J.H. Vlugt, D. Dubbeldam
> Implementation of energy and force grids in molecular simulations of porous 
> materials.
> Molecular Physics, 2025, e2542486
> <https://doi.org/10.1080/00268976.2025.2542486>