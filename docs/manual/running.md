# Running `RASPA`
\page running Running

## Input file

An input-file describing the type of simulation and the parameters. In the same directory as the 'run'-file, there needs to be a file called `simulation.json`. An example file is:
```json
{
    "SimulationType" : "MonteCarlo",
    "NumberOfCycles" : 100000,
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
