module;

export module mc_backend;

import std;

import sampling_backend;

// The four sampling primitives on the processor, spread over its cores with OpenMP. Each is a plain loop
// over the cases, and each case a plain loop over the atoms of the cell; there is nothing in here that the
// GPU version does differently except where it runs.
export SamplingBackend samplingBackendCPU();
