export module sample_movies;

import atom;
import simulationbox;
import forcefield;

import <vector>;
import <numeric>;
import <fstream>;
import <utility>;
import <string>;


export struct SampleMovie
{
    SampleMovie(size_t systemId, const ForceField &forceField, const SimulationBox& simulationBox, 
                const std::vector<Atom>& atomPositions);
    //SampleMovie() noexcept = default;
    //SampleMovie(const SampleMovie& a) noexcept = default;
    //SampleMovie& operator=(const SampleMovie& a) noexcept = default;
    SampleMovie(SampleMovie&& a) noexcept;
    //SampleMovie& operator=(SampleMovie&& a) noexcept = default;
    //~SampleMovie() = default;

    const size_t systemId;
    const ForceField& forceField;
    const SimulationBox& simulationBox; 
    const std::vector<Atom>& atomPositions;

    void initialize();
    void update(size_t cycle);
    void closeOutputFile();

    size_t writeEvery{ 5000 };
    bool sample{ false };

    int modelNumber{ 1 };

    //std::ofstream outputFile{};
};
