module;

module system;

// system_integration.cpp

/*
[[nodiscard]] EnergyBookkeeping System::integrate() noexcept
{
    for (size_t i = 0; i < adsorbateAtomVelocities.size(); ++i)
    {
        adsorbateAtomVelocities[i] += 0.5 * adsorbateAtomForces[i] * timeStep;
        adsorbateAtomPositions[i].position += adsorbateAtomVelocities[i] * timeStep;
    }

    EnergyBookkeeping energies = computeTotalForces();

    for (size_t i = 0; i < adsorbateAtomVelocities.size(); ++i)
    {
        adsorbateAtomVelocities[i] += 0.5 * adsorbateAtomForces[i] * timeStep;
    }

    return energies;
}
*/
