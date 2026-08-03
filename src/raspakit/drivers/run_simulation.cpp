module;

module run_simulation;

import std;

import archive;
import threadpool;
import input_reader;
import monte_carlo;
import monte_carlo_transition_matrix;
import molecular_dynamics;
import minimization;
import thermodynamic_integration;
import parallel_thermodynamic_integration;
import parallel_tempering;
import hyper_parallel_tempering;
import reweighted_histogram;
import parallel_tmmc;

void runSimulation(InputReader& inputReader)
{
  auto& pool = ThreadPool::ThreadPool<ThreadPool::details::default_function_type>::instance();
  pool.init(inputReader.numberOfThreads, inputReader.threadingType);

  switch (inputReader.simulationType)
  {
    case InputReader::SimulationType::MonteCarlo:
    {
      MonteCarlo mc(inputReader);
      if (inputReader.restartFromBinary)
      {
        readBinaryRestartFile(mc, inputReader.restartFromBinaryFileName);
        // the output files are opened in append mode (the header is already in the file), and
        // the interpolation grids are not stored in the restart file and must be rebuilt
        mc.createOutputFiles();
        mc.createInterpolationGrids();
      }

      mc.run();
      break;
    }
    case InputReader::SimulationType::MonteCarloTransitionMatrix:
    {
      MonteCarloTransitionMatrix mc(inputReader);
      if (inputReader.restartFromBinary)
      {
        readBinaryRestartFile(mc, inputReader.restartFromBinaryFileName);
        // the resumed stage jumps past the header block, so the output streams must exist
        mc.createOutputFiles();
      }
      mc.run();
      break;
    }
    case InputReader::SimulationType::MolecularDynamics:
    {
      MolecularDynamics md(inputReader);
      if (inputReader.restartFromBinary)
      {
        readBinaryRestartFile(md, inputReader.restartFromBinaryFileName);
        md.createOutputFiles();
        // the grids are not stored in the restart file and must be rebuilt
        md.createInterpolationGrids();
      }

      md.run();
      break;
    }
    case InputReader::SimulationType::Minimization:
    {
      Minimization minimization(inputReader);
      if (inputReader.restartFromBinary)
      {
        readBinaryRestartFile(minimization, inputReader.restartFromBinaryFileName);
      }
      minimization.run();
      break;
    }
    case InputReader::SimulationType::ThermodynamicIntegration:
    {
      ThermodynamicIntegration ti(inputReader);
      ti.run();
      break;
    }
    case InputReader::SimulationType::ParallelThermodynamicIntegration:
    {
      ParallelThermodynamicIntegration parallel_ti(inputReader);
      if (inputReader.restartFromBinary)
      {
        readBinaryRestartFile(parallel_ti, inputReader.restartFromBinaryFileName);
      }
      parallel_ti.run();
      break;
    }
    case InputReader::SimulationType::ParallelTempering:
    {
      ParallelTempering parallel_tempering(inputReader);
      if (inputReader.restartFromBinary)
      {
        readBinaryRestartFile(parallel_tempering, inputReader.restartFromBinaryFileName);
      }
      parallel_tempering.run();
      break;
    }
    case InputReader::SimulationType::HyperParallelTempering:
    {
      HyperParallelTempering hyper_parallel_tempering(inputReader);
      if (inputReader.restartFromBinary)
      {
        readBinaryRestartFile(hyper_parallel_tempering, inputReader.restartFromBinaryFileName);
      }
      hyper_parallel_tempering.run();
      break;
    }
    case InputReader::SimulationType::ReweightedHistogram:
    {
      ReweightedHistogram reweighted_histogram(inputReader);
      if (inputReader.restartFromBinary)
      {
        readBinaryRestartFile(reweighted_histogram, inputReader.restartFromBinaryFileName);
      }
      reweighted_histogram.run();
      break;
    }
    case InputReader::SimulationType::ParallelTMMC:
    {
      ParallelTMMC parallel_tmmc(inputReader);
      if (inputReader.restartFromBinary)
      {
        readBinaryRestartFile(parallel_tmmc, inputReader.restartFromBinaryFileName);
      }
      parallel_tmmc.run();
      break;
    }
    default:
      break;
  }
}
