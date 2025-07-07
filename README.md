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
<img alt="GitHub Actions Workflow Status" src="https://img.shields.io/github/actions/workflow/status/iRASPA/RASPA3/create-release.yml">
  </a>
  <a href="https://github.com/iRASPA/raspa3/issues"><img alt="GitHub Issues or Pull Requests" src="https://img.shields.io/github/issues/iRASPA/raspa3">
</a>
 <a href="https://iRASPA.github.io/RASPA3"><img alt="Documentation" src="https://img.shields.io/badge/documentation-blue">
 </a>
 <a href="https://github.com/iRASPA/RASPA3/actions/workflows/unit-test.yml"><img alt="Unittests" src="https://github.com/iraspa/raspa3/actions/workflows/unit-test.yml/badge.svg">
 </a>
</p>

<p align="center">
  •
  <a href="#authors">Authors</a> •
  <a href="#contributors">Contributors</a> •
  <a href="#running">Running</a> •
  <a href="#python">Python</a> •
  <a href="#dependencies">Dependencies</a> •
  <a href="#installation-guide">Installation Guide</a> •
</p>

# Authors

Drs. Youri Ran, University of Amsterdam<br>
Drs. Shrinjay Sharma, Delft University of Technology<br>
Dr. Salvador R.G. Balestra, Universidad Pablo de Olavide<br>
Drs. Zhao Li, Northwestern University<br>
Prof. Sofia Calero,  Eindhoven University of Technology<br>
Prof. Thijs Vlugt, Delft University of Technology<br>
Prof. Randall Q. Snurr, Northwestern University<br>
Dr. David Dubbeldam, University of Amsterdam

# Contributors

Alvaro Vazquez Mayagoitia, Argonne National Lab, contribution to openmp-implementation discussion<br>
Anserme, better README.md and packaging

# Citing RASPA3

Y.A. Ran, S. Sharma, S.R.G. Balestra, Z. Li, S. Calero, T.J.H. Vlugt, R.Q. Snurr, D. Dubbeldam, _"RASPA3: A Monte Carlo code for computing adsorption and diffusion in nanoporous materials and thermodynamics properties of fluids"_, **2024**, J. Chem. Phys., 161, 114106, [DOI](https://doi.org/10.1063/5.0226249)

# Running

```
cd examples/basic/1_mc_methane_in_box
./run
```

# Implemented so far
- rigid molecules and rigid framework
- multiple systems, box or framework
- grand canonical ensemble (CBMC, CFCMC, and CB/CFCMC)
- Gibbs ensemble (CBMC, CFCMC, and CB/CFCMC)
- Monte Carlo NPT ensemble
- transition matrix Monte Carlo
- Molecular Dynamics NVT ensemble (Nose-Hoover thermostat)
- binary restart
- blocking pockets
- PDB-movies, energy histograms, number of molecule histograms
- RDF, MSD-order-N, VACF, density grids (cube files)
- charge equilibration
- MC/MD hybrid move
- tail-corrections for CFCMC
- grids for rigid frameworks

# Todo-list
- restart-file
- flexible molecules
- flexible frameworks
- reaction ensemble
- identity change
- polarization
- cell-lists for rigid frameworks
- partial molar volumes
- zeo++-type calculations
- partial molar enthalpies/volumes
- optimization
- elastic constants
- HDF5 property writing
