module;

export module energy_opencl_backend;

import energy_shared_energy_backend;

// The GPU's answer to the four field builders. Handing this to a driver is the whole of choosing the GPU
// route; nothing else in the energy-based properties knows which backend it is running on.
export EnergyBackend openCLEnergyBackend();
