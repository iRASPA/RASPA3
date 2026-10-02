module;

#define CL_TARGET_OPENCL_VERSION 120
#ifdef __APPLE__
#include <OpenCL/cl.h>
#elif _WIN32
#include <CL/cl.h>
#else
#include <CL/opencl.h>
#endif

export module spatial_decomposition_opencl_bonded;

import std;

import double3;
import double3x3;
import spatial_decomposition_opencl_handles;
import spatial_decomposition_device_kernels;
import spatial_decomposition_device_backend;
import spatial_decomposition_device_bonded_topology;

/**
 * \brief The per-molecule terms on the OpenCL device: Ewald self and exclusion corrections, bonds, bends,
 * torsions and improper torsions, intramolecular Lennard-Jones and Coulomb pairs, and the atomic-to-molecular
 * virial correction.
 *
 * Two kernels (the shared bonded_kernel_source.cpp compiled with the OpenCL dialect header): one work-item per
 * term instance evaluates the term once and writes the gradient on each of its atoms to a per-instance slot; one
 * work-item per slot (atom) evaluates the exclusion pairs and the virial correction, gathers the gradients of its
 * terms and adds the sum to the device force of the slot. Both leave per-group partial sums of the energies, the
 * exclusion strain derivative and the virial correction, which the host sums in double. The topology
 * (BondedTopology, built by the neutral DeviceStep) is uploaded once; the slot layout (molecule and index of
 * every slot) after every list build. The positions of the terms are the float positions relative to the first
 * atom of the molecule, indexed by the atom's index in the system, packed by the host.
 */
export class OpenCLBonded
{
 public:
  OpenCLBonded() = default;
  ~OpenCLBonded() = default;
  OpenCLBonded(const OpenCLBonded&) = delete;
  OpenCLBonded& operator=(const OpenCLBonded&) = delete;
  OpenCLBonded(OpenCLBonded&&) noexcept = default;
  OpenCLBonded& operator=(OpenCLBonded&&) noexcept = default;

  void initialize(cl_context context, cl_device_id device);
  bool initialized() const { return atomKernel.get() != nullptr; }

  /// Uploads the terms of the components (per atom), the molecule table and the masses (blocking, shared queue).
  void setTopology(const BondedTopology& topology);
  /// The Ewald parameters (`useCharge`: self and exclusion corrections on).
  void setParameters(double alpha, double conversionFactor, bool useCharge);
  /// Follows a change of the Ewald alpha.
  void setAlpha(double alpha);

  /// The slot layout of a list build (BondedTopology::layout); enqueues the upload on `queue` (non-blocking; the
  /// span must stay valid until the queue completed it).
  void setLayout(cl_command_queue queue, std::span<const std::uint32_t> slotMolecule);

  /// Enqueues the kernels (after the pair kernel and the mesh interpolation wrote `force`); `relative` holds the
  /// relative position and charge of every atom in the system order.
  void enqueue(cl_command_queue queue, cl_mem relative, cl_mem force);
  /// Enqueues the read-back of the partial sums; `event` receives the read's event.
  void enqueueRead(cl_command_queue queue, cl_event* event);

  /// After the read completed: the sums over all molecules.
  DeviceBondedResults collect() const;

  std::string status() const;

 private:
  /// Mirrors the BondedParameters struct of the kernel source.
  struct Parameters
  {
    float alpha{0.0f};
    float twoAlphaOverSqrtPi{0.0f};
    float selfPrefactor{0.0f};
    float coulombFactor{0.0f};
    std::uint32_t useCharge{0};
    std::uint32_t numberOfSlots{0};
    std::uint32_t numberOfInstances{0};
    std::uint32_t atomPartialOffset{0};
  };
  static_assert(sizeof(Parameters) == 32);

  cl_context context{nullptr};
  cl_device_id device{nullptr};
  OpenCLDevice::ProgramHandle program{};
  OpenCLDevice::KernelHandle termKernel{}, atomKernel{};
  OpenCLDevice::MemHandle slotMoleculeBuffer{}, moleculeInfoBuffer{}, instanceMoleculeBuffer{}, massBuffer{};
  OpenCLDevice::MemHandle termsBuffer{}, gradientOffsetBuffer{}, atomGradientStartBuffer{}, atomGradientsBuffer{};
  OpenCLDevice::MemHandle termGradientBuffer{}, parameterBuffer{}, partialBuffer{};
  std::size_t slotCapacity{0}, partialCapacity{0};

  Parameters parameters{};
  bool parametersChanged{true};
  double alphaValue{0.0}, conversionFactor{0.0}, chargeSquaredSum{0.0};
  std::size_t numberOfMolecules{0}, numberOfTerms{0}, termGroups{0};

  std::vector<float> hostPartials{};
  std::size_t atomGroups{0};
};
