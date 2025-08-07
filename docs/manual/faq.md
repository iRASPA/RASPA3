# FAQ
\page FAQ FAQ

If your issue is not in here, feel free to open a new issue on [Github](https://www.github.com/iRASPA/RASPA3/issues).

### Error while loading libraries
This error can occur when raspa can not find the libraries it is linked against, such as BLAS and LAPACK, which are used for linear algebra. Finding all linked libraries for the executable  can be done with `ldd /path/to/raspa3` (Linux) or `otool -L /path/to/raspa3` (MacOS). If there are unlinked libraries, locate them on your system and add them to your path using `export LD_LIBRARY_PATH=/path/to/lib:$LD_LIBRARY_PATH`. 

```
raspa3: error while loading shared libraries: libopenblas.so.0: cannot open shared object file: No such file or directory
```

### Excess loading is negative
This usually happens when computing an isotherm and the next pressure is above the vapor pressure. The boundary from gas to liquid adsorption has been crossed and the amount of excess molecules increases by orders of magnitude. There is a reason why experimental gas-phase isotherms are of finite range, they usually stop at the vapor pressure. Also, if the pressure is very high the fluid outside the crystal is compressed more and more while the loading inside the crystal remains the same (at maximum loading). Hence, excess adsorption evantually becomes negative.

### Large drift in Monte Carlo energies 
This should _not_ happen and signals an error in (one of) the Monte Carlo routines. During the Monte Carlo simulations, the running-energies are stored. These are starting energy, and all the added energy _differences_. At the final stage, the energy is recomputed again, and these should match within an error of about $10^{−5}$ or lower. If you have added your own MC move, check whether you have properly added the energy differences to the running energies.

### Energy is not conserved in molecular dynamics
Usually, this happens because the time step is too large. Also, at initialization, the system can be far from equilibrated and a smaller time step is needed.

### RASPA ”hangs” at initialization
Put ’CreateNumberOfMolecules 0’ and check if that solves the problem. If so, then you have tried to add too many molecules in the system (i.e. more than actually fit in the system). For systems without a framework, one can also increase the size of the box.

### Mean-square displacement is not linear
There are several known causes:
- Your system is one-dimensional and particles are unable to pass each other. This is known as ’single- file-diffusion’ and the mean square displacement is propertial to the square root of time,
- You did not simulate long enough. In some systems it can take up to several nanoseoconds before the msd becomes linear in time,
- You forgot to specify interactions between the molecules and they are not interacting.

### Strange behaviour when using cations
The problem could be related to CBMC of net-charged molecules. The Rosenbluth weights can become very large or very small, because the energy difference when displacing an ion is large. This can lead to numerical problems for ratio’s of combination of small/large, for example in the reinsertion move. To see whether this is the cause use ’RandomTranslationProbability’ and set ’ReinsertionProbability’ to zero.

### Parallel tempering does not work for systems with cations 
The problem is physical, it is just that all the energy distributions are more ’spiked’ and overlap between the energies of the systems is rare. The only solution is to use more systems (and smaller temperature differences) to increase the acceptance rates of swapped between neighboring system. The same problem happens when one increases the system size.

### Results do not match data from literature 
Common reasons include difference in simulation length, system size, cutoff value, tail-corrections vs shifted potentials, and handling of electrostatics. For adsorption, is the crystal structure that you used the same? Another very often made error, is comparing against different units, or an error in the conversion of units.