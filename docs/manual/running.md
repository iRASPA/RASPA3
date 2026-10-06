# Running `RASPA`
\page running Running

## Input file

An input-file describing the type of simulation and the parameters. In the same directory as the 'run'-file, there needs to be a file called `simulation.json`. An example file is:
```json
{
    "SimulationType" : "MonteCarlo",
    "NumberOfProductionCycles" : 100000,
    "NumberOfInitializationCycles" : 1000,
    "NumberOfEquilibrationCycles" : 10000,
    "PrintEvery" : 1000,

    "Systems" :
    [
    {
        "Type" : "Box",
        "BoxLengths" : [30.0, 30.0, 30.0],
        "ExternalTemperature" : 300.0,
        "ChargeMethod" : "None",
        "OutputPDBMovie" : true,
        "SampleMovieEvery" : 10
    }
    ],

    "Components" :
    [
    {
        "Name" : "methane",
        "MoleculeDefinition" : "ExampleDefinitions",
        "TranslationProbability" : 1.0,
        "CreateNumberOfMolecules" : 100
    }
    ]
}
```
This tells `RASPA` to run a Monte-Carlo simulation of 100 methane
molecules in a $30\times30\times30$ Å cubic box (with 90$^\circ$
angles) at 300 Kelvin. It will start with 1000 cycles to initialize
the system, 10000 cycles to equilibrate the system, and will use
100000 cycle to obtain thermodynamic properties of interest. Every
1000 cycles a status-report is printed to the output. The
Monte-Carlo program will use only the 'translation move' where a
particle is given a random translation and the move is accepted or
rejected based on the energy difference.

Further settings for writing input files can be found in the ![commands](docs/manual/commands.md) section.

----------------------------------------------------------------------------------

## Importing topologies from other codes

`raspa3` can write a complete set of input files from the topology files of other simulation
packages and exit. The converters produce `force_field.json`, one `<component>.json` per distinct
molecule and a `simulation.json` skeleton in the directory given with `--output-dir`; edit the
skeleton (number of cycles, ensemble, moves) before running. Warnings about terms that could not be
converted are listed in `conversion_warnings.txt`.

### AMBER (`prmtop` / `inpcrd`)

    raspa3 --from-prmtop system.prmtop --inpcrd system.inpcrd --output-dir raspa-input

| Option | Meaning |
|---|---|
| `--from-prmtop FILE` | AMBER7-format topology (`%VERSION`, `%FLAG` sections) |
| `--inpcrd FILE` | ASCII coordinate file (`inpcrd` or `rst7`; NetCDF restarts must be converted with `cpptraj` first). Positions and the box are written to `restart.json`, referenced from `simulation.json` through `RestartFileName`, so the simulation starts from the AMBER configuration |
| `--cutoff R` | van der Waals and real-space Coulomb cut-off in Å (default 10) |
| `--truncation M` | van der Waals truncation written to `force_field.json`: `truncated` (default), `shifted`, `switched` (the OpenMM/GROMACS potential switch) or `force-switched` (the CHARMM `vfswitch` force switch) |
| `--switching-distance R` | where the switching of `switched`/`force-switched` starts, in Å (default: 2 Å below the cut-off) |
| `--output-dir DIR` | output directory (default `raspa-from-prmtop`) |

The conversion:

- One pseudo-atom per distinct AMBER atom type (`AMBER_ATOM_TYPE`), with Lennard-Jones parameters
  from `LENNARD_JONES_ACOEF`/`BCOEF` ($\epsilon = B^2/4A$, $\sigma = (A/B)^{1/6}$). Pairs that deviate
  from Lorentz-Berthelot mixing (NBFIX-style ion parameters) are written as `BinaryInteractions`.
  Charges are per atom (`CHARGE` divided by 18.2223; by $\sqrt{332.0716}$ for a CHAMBER topology, which
  uses the CHARMM Coulomb constant). A CHAMBER topology (CHARMM force field converted with ParmEd) carries a
  separate Lennard-Jones table for the 1-4 pairs (`LENNARD_JONES_14_ACOEF`/`BCOEF`); it is written as
  `parameters14` of the self and binary interactions.
- Bonds and angles become `HARMONIC` terms ($p_0 = 2k$, since AMBER writes $k(r-r_0)^2$ and RASPA
  $\tfrac{1}{2}p_0(r-p_1)^2$). All dihedral terms of a bonded quadruple with a phase of 0 or 180
  degrees and periodicity $\leq 5$ are summed into one `POLYNOMIAL` in $\cos\phi$ (exact, including the
  constant); terms with other phases or higher periodicities, and the AMBER impropers, are written as
  `CVFF` entries under `ImproperTorsions`. `SCEE`/`SCNB` become `Intra14ChargeChargeScalingValue` and
  `Intra14VanDerWaalsScalingValue`.
- CMAP corrections (ff19SB `CMAP_*`, CHAMBER `CHARMM_CMAP_*`): every map becomes an entry of `CMAPs` in
  `force_field.json` (named `CMAP_1`, `CMAP_2`, ...; grid in K) and every five-atom term an entry of
  `CMAPTorsions` in the component file.
- Molecules are the connected components of the bond graph, grouped into components by their types,
  charges and parametrised topology. A single-residue molecule is named after its residue (`WAT`,
  `Na_plus`, `Cl-`), a chain of amino acids `protein`, a chain of nucleotides `nucleic_acid`.
- Water (residues `WAT`, `HOH`, `TIP3`, `TP3`, `SOL`) and single atoms are rigid components. Water takes the
  ideal geometry from its equilibrium bond lengths (and angle, or the H-H bond of a SETTLE topology), and is
  integrated as a rigid body (centre of mass and quaternion); the AMBER coordinates are fitted onto that
  geometry when the configuration is read.
- With a box (`inpcrd` box line or `BOX_DIMENSIONS`) the system is periodic with Ewald summation. Without
  one, a vacuum box is built around the molecules with direct Coulomb and a cut-off spanning every atom pair.

Single-point energies of a converted ff14SB peptide (ACE-ALA-PRO-PHE-NME) agree with OpenMM to $10^{-5}$ kcal/mol
per term in vacuum, and to $4\times10^{-6}$ relative (the Ewald/PME precision) for the same peptide in 437 TIP3P
waters with NaCl; with `--truncation switched` the Lennard-Jones energy of the solvated system agrees with
OpenMM's switching function to $10^{-7}$ relative. The CMAP energy of an ff19SB peptide (ACE-ALA-GLY-ALA-NME,
three CMAP terms) agrees with OpenMM's `CMAPTorsionForce` to $2\times10^{-7}$ kcal/mol (the precision of the
grid in the prmtop), and the 1-4 Lennard-Jones table of a CHAMBER topology to $10^{-6}$ kcal/mol.

Not converted: 10-12 hydrogen-bond terms, polarizabilities and GB radii. The Coulomb constant
differs slightly between AMBER (332.0522 kcal Å/mol) and RASPA's CODATA value (332.0637 kcal Å/mol), a
relative difference of $3.5\times10^{-5}$ in the electrostatic energies.

### LAMMPS (`read_data`)

    raspa3 --from-lammps data.lammps --lammps-input in.lammps --output-dir raspa-input

The data file provides masses, pair coefficients, bonds, angles, dihedrals and molecules; the input
script (optional) names the styles, `special_bonds` and `pair_modify` settings.

----------------------------------------------------------------------------------

## RASPA_DIR

The RASPA_DIR should be linked, as many of the default .cif and force field files are given there. Make sure the the environment variable `RASPA_DIR` is set up correctly or run using the run file, which automatically includes the right path after building.

`run` file:
```
#! /bin/sh -f
export RASPA_DIR=/usr/share/raspa3
/usr/bin/raspa3
```
This type of file is know as a '`shell script`'. `RASPA` needs the
variable '`RASPA_DIR`' to be set in order to look up the molecules,
frameworks, etc. The scripts sets the variable and runs `RASPA`.
`RASPA` can then be run from any directory you would like.

----------------------------------------------------------------------------------

## Run job on cluster

In order to run it on a cluster using a queuing system one needs an
additional file '`bsub.job`' (arbitrary name)

-   '`gridengine`'

             #!/bin/bash
             # Serial sample script for Grid Engine
             # Replace items enclosed by {}
             #$ -S /bin/bash
             #$ -N Test
             #$ -V
             #$ -cwd
             echo $PBS_JOBID > jobid
             export RASPA_DIR=/usr/share/raspa3
             /usr/bin/raspa3

    The job can be submitted using '`qsub bsub.job`'.

-   '`torque`'

             #!/bin/bash
             #PBS -N Test
             #PBS -o pbs.out
             #PBS -e pbs.err
             #PBS -r n
             #PBS -V
             #PBS -mba
             cd $PBS_O_WORKDIR
             echo $PBS_JOBID > jobid
             export RASPA_DIR=/usr/share/raspa3
             /usr/bin/raspa3

    The job can be submitted using '`qsub bsub.job`'.

-   '`slurm`'

          #!/bin/bash 
          #SBATCH -N 1
          #SBATCH --job-name=Test
          #SBATCH --export=ALL
          echo $SLURM_JOBID > jobid
          valhost=$SLURM_JOB_NODELIST
          echo $valhost > hostname
          module load slurm
          export RASPA_DIR=/usr/share/raspa3
          /usr/bin/raspa3

    The job can be submitted using '`sbatch bsub.job`'.
