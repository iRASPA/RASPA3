export module monte_carlo;

import randomnumbers;
import threadpool;

import averages;
import system;
import mc_moves;
import input_reader;
import energy_status;

import <vector>;
import <iostream>;
import <fstream>;
import <chrono>;

export struct MonteCarlo
{
    enum class WeightingMethod : size_t
    {
        LambdaZero = 0,
        AllLambdas = 1
    };

    MonteCarlo(InputReader& reader) noexcept;

    void run();
    void initialize();
    void equilibrate();
    void production();
    void output();
    void cleanup();

    System& randomSystem();

    size_t numberOfCycles;
    size_t numberOfInitializationCycles;
    size_t numberOfEquilibrationCycles;
    size_t optimizeMCMovesEvery{ 5000 };
    size_t printEvery;

    std::vector<System> systems;
    size_t fractionalMoleculeSystem{ 0 };   // the system where the fractional molecule is located

    std::vector<std::ofstream> streams;

    BlockErrorEstimation estimation;

    MC_Moves particleMoves;

    std::chrono::system_clock::time_point t1;
    std::chrono::system_clock::time_point t2;
};
