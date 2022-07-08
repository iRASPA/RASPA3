 1) Change to return type to 'double' and compute energy with pressure-routine.
 2) compute pressure (Ewald)
 3) volume move (NPT)
 4) rescale MC moves max-change 
 5) identity change
 6) Gibbs ensemble
 7) Molecular Dynamics
 8) read framework
 9) read molecule
10) make routines cutoff-dependence for dual-cutoff method
11) reaction ensemble
12) polarization
13) cell-lists for rigid frameworks
14) grids for rigid frameworks
15) partial molar volumes

fix
===
1) density in kg/m^3 (multiply by mass)


long-term:
==========
1) gcc makefiles, but gcc-11 and gcc-12 not working so far (bugs in module implementation).



