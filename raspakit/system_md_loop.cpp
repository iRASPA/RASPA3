module;

module system;

import <iostream>;

// system_md_loop.cpp

void System::MD_Loop()
{
    /*
    energyBookkeeping.initialize();

    EnergyBookkeeping energies = integrate();

    double totalEnergy = energies.sum();

    std::cout << "totalEnergy: " << totalEnergy;

    totalEnergy = 0.0;
        for (int i = 0; i < 2000; i++)
        {
            //totalEnergy += interMolecularInteractions.computeEnergy(0, i);

        }
    std::cout << "totalEnergy: " << totalEnergy/2.0;

    //std::optional<std::pair<double, std::vector<double3>>> rosenbluthWeight = cbmc.growMoleculeInsertion(0, 1.0);
    

    for (int i = 0; i < 100000; i++)
    {
        EnergyBookkeeping energies = integrate();

        double totalEnergy = energies.sum();

        double kineticEnergy = 0.0;
        for (const double3& vel : adsorbateAtomVelocities)
        {
            kineticEnergy += 0.5 * (vel.x * vel.x + vel.y * vel.y + vel.z * vel.z);
        }

        if (i % 1000 == 0)
        {
            std::cout << "totalEnergy: " << totalEnergy + kineticEnergy << " " << totalEnergy << " " << kineticEnergy << std::endl;
        }
    }*/
};
