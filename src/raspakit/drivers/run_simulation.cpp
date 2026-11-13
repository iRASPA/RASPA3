module;

module run_simulation;

import std;
import system;

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
import parallel_tempering_molecular_dynamics;
import molecular_dynamics_spatial_decomposition;
import hyper_parallel_tempering;
import reweighted_histogram;
import parallel_tmmc;
import nldft;
import dctst;
import dctst_free_energy;
import polymerization;

void runSimulation(InputReader& inputReader)
{
  auto& pool = ThreadPool::ThreadPool<ThreadPool::details::default_function_type>::instance();
  pool.init(inputReader.numberOfThreads, inputReader.threadingType);

  // The device of the multiple-try CBMC grower is a property of the machine running the job, not of the state: it
  // is not stored in a binary restart (a resumed run may lack the device), and the systems read from one come back
  // with the default. The choice of the input the run was started with is put back after reading.
  auto applyInputDevices = [&](std::vector<System>& systems)
  {
    for (System& system : systems) system.cbmcDevice = inputReader.cbmcDevice;
  };

  switch (inputReader.simulationType)
  {
    case InputReader::SimulationType::MonteCarlo:
    {
      MonteCarlo mc(inputReader);
      if (inputReader.restartFromBinary)
      {
        readBinaryRestartFile(mc, inputReader.restartFromBinaryFileName);
        applyInputDevices(mc.systems);
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
        applyInputDevices(mc.systems);
        // the resumed stage jumps past setup(), so the output streams must exist, and the
        // interpolation grids are not stored in the restart file and must be rebuilt
        mc.createOutputFiles();
        mc.createInterpolationGrids();
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
        applyInputDevices(md.systems);
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
        applyInputDevices(minimization.systems);
      }
      minimization.run();
      break;
    }
    case InputReader::SimulationType::ThermodynamicIntegration:
    {
      // this driver has no restart file format: continuing from scratch without a word would look like a resume
      if (inputReader.restartFromBinary)
      {
        throw std::runtime_error("[Input reader]: 'ThermodynamicIntegration' does not support 'RestartFromBinary'\n");
      }
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
        applyInputDevices(parallel_ti.systems);
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
        applyInputDevices(parallel_tempering.systems);
      }
      parallel_tempering.run();
      break;
    }
    case InputReader::SimulationType::ParallelTemperingMolecularDynamics:
    {
      ParallelTemperingMolecularDynamics parallel_tempering_md(inputReader);
      if (inputReader.restartFromBinary)
      {
        readBinaryRestartFile(parallel_tempering_md, inputReader.restartFromBinaryFileName);
        applyInputDevices(parallel_tempering_md.systems);
      }
      parallel_tempering_md.run();
      break;
    }
    case InputReader::SimulationType::MolecularDynamicsSpatialDecomposition:
    {
      MolecularDynamicsSpatialDecomposition md(inputReader);
      if (inputReader.restartFromBinary)
      {
        readBinaryRestartFile(md, inputReader.restartFromBinaryFileName);
        applyInputDevices(md.systems);
        md.createOutputFiles();
        // the force engines are not stored in the restart file; the resumed stage rebuilds them
      }
      md.run();
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
    case InputReader::SimulationType::NLDFT:
    {
      if (inputReader.restartFromBinary)
      {
        throw std::runtime_error("[Input reader]: 'NLDFT' does not support 'RestartFromBinary'\n");
      }
      NLDFT nldft(inputReader);
      nldft.run();
      break;
    }
    case InputReader::SimulationType::DynamicallyCorrectedTST:
    {
      if (inputReader.restartFromBinary)
      {
        throw std::runtime_error("[Input reader]: 'DynamicallyCorrectedTST' does not support 'RestartFromBinary'\n");
      }
      DynamicallyCorrectedTST dctst(inputReader);
      dctst.run();
      break;
    }
    case InputReader::SimulationType::DynamicallyCorrectedTSTFreeEnergy:
    {
      if (inputReader.restartFromBinary)
      {
        throw std::runtime_error(
            "[Input reader]: 'DynamicallyCorrectedTSTFreeEnergy' does not support 'RestartFromBinary'\n");
      }
      DynamicallyCorrectedTSTFreeEnergy dctst(inputReader);
      dctst.run();
      break;
    }
    case InputReader::SimulationType::Polymerization:
    {
      Polymerization polymerization(inputReader);
      if (inputReader.restartFromBinary)
      {
        readBinaryRestartFile(polymerization, inputReader.restartFromBinaryFileName);
        polymerization.createOutputFiles();
      }
      polymerization.run();
      break;
    }
    default:
      break;
  }
}
