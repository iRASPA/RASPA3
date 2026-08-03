module;

export module energy_backend;

import energy_shared_energy_backend;

// The processor's answer to the four field builders. Handing this to a driver is the whole of choosing the
// CPU route; nothing else in the energy-based properties knows which backend it is running on.
export EnergyBackend cpuEnergyBackend();
