module;

export module mc_opencl_backend;

import std;

import sampling_backend;

// The four sampling primitives on the GPU.
//
// What they compute is defined in `sampling_backend` and is the same on either side; what is worth saying
// here is why these four and not the whole construction. Building a roadmap is a graph algorithm wrapped
// around an enormous amount of arithmetic about spheres. The graph part -- binning the samples, walking the
// components, ranking the hops, deduplicating the pockets -- is branchy, serial and cheap, and putting it
// on a device would cost more in staging than it could ever save. The arithmetic part is a loop over the
// atoms of the cell run some hundreds of thousands of times over independent cases, which is the shape of
// problem a device exists for. So the two are split along that line and only the second crosses over.
//
// Each call stages its cases into device memory, runs one work item per case and reads the answers back.
// The atoms are staged once per call rather than once per case, since every case reads all of them; on a
// cell of a few hundred atoms they sit in cache and the kernels run at the speed of the arithmetic.
//
// The device works in single precision, which is what it is fast at, and the answers come back as doubles.
// That sets a floor on how finely a clearance can be told, of order a part in ten million of the cell,
// which is some four orders below the convergence of the sample itself: two runs of the same size on the
// two backends differ by far more than this because they are different random streams, and neither differs
// from the exact answer by anything like as little.
export SamplingBackend samplingBackendOpenCL();

// Whether there is a device to run on. A run that asks for the GPU without one gets told so rather than
// getting the processor's answer under the GPU's name, the two being distinguishable in the reports.
export bool samplingOpenCLAvailable();
