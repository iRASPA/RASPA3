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
import system;
import spatial_decomposition_opencl_handles;

/// The OpenCL C source of the bonded kernels (bonded_kernel_source.cpp).
extern const char* const openclBondedKernelSource;

/**
 * \brief The per-molecule terms on the OpenCL device: Ewald self and exclusion corrections, bonds, bends,
 * torsions and improper torsions, intramolecular Lennard-Jones and Coulomb pairs, and the atomic-to-molecular
 * virial correction.
 *
 * One work-item per slot gathers the terms acting on its atom (see bonded_kernel_source.cpp), adds the gradient
 * to the device force of the slot and leaves per-group partial sums of the energies, the exclusion strain
 * derivative and the virial correction, which the host sums in double. The topology (the terms of every
 * component, listed per atom) is uploaded once; the slot layout (molecule and index of every slot, slot of every
 * atom) after every list build. The positions of the terms are the float positions relative to the first atom of
 * the molecule, packed by the host (OpenCLPairKernel::packPositions) with the reference atoms this class provides.
 *
 * Covers all BondType, BendType and TorsionType potentials and the intramolecular Lennard-Jones and Coulomb pairs;
 * `supports` reports the components that use other intramolecular terms (Urey-Bradley, inversion bends, cross
 * terms), molecules with more than 256 atoms and lambda-group atoms, for which the engine keeps the bonded work on
 * the host.
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

  /// Whether the device kernels cover the intramolecular terms of the system; on false, `reason` names the first
  /// unsupported feature.
  static bool supports(const System& system, std::string& reason);

  void initialize(cl_context context, cl_device_id device);
  bool initialized() const { return kernel.get() != nullptr; }

  /// Uploads the terms of the components (per atom), the molecule table and the masses of the pseudo-atom types.
  void setTopology(const System& system);
  /// The Ewald parameters (`useCharge`: self and exclusion corrections on).
  void setParameters(double alpha, double conversionFactor, bool useCharge);
  /// Follows a change of the Ewald alpha.
  void setAlpha(double alpha);

  /// Maps the slot layout of a list build: `slotOfSorted` per cell-list (sorted) atom and `originalToSorted` per
  /// atom in the system order. Enqueues the uploads on `queue` (non-blocking; the host arrays are members).
  void setLayout(cl_command_queue queue, std::span<const std::uint32_t> slotOfSorted,
                 std::span<const std::uint32_t> originalToSorted, std::size_t slots);

  /// Per sorted atom: the sorted index of the first atom of its molecule (the reference of the relative positions).
  std::span<const std::uint32_t> referenceAtoms() const { return referenceOfSorted; }

  /// Enqueues the kernel over the slots (after the pair kernel and the mesh interpolation wrote `force`).
  void enqueue(cl_command_queue queue, cl_mem position, cl_mem relative, cl_mem typeOf, cl_mem force);
  /// Enqueues the read-back of the partial sums; `event` receives the read's event.
  void enqueueRead(cl_command_queue queue, cl_event* event);

  struct Results
  {
    /// The self energy is evaluated in double on the host (it depends on the charges and alpha only); the
    /// device evaluates self + exclusion without their mutual cancellation and `exclusion` is the difference.
    double self{0.0}, exclusion{0.0}, bond{0.0}, bend{0.0}, torsion{0.0}, improperTorsion{0.0};
    double intraVDW{0.0}, intraCoulomb{0.0};
    double3x3 exclusionStrain{};
    double3x3 correction{};
  };
  /// After the read completed: the sums over all molecules.
  Results collect() const;

  std::size_t numberOfTerms() const { return terms.size(); }
  std::string status() const;

 private:
  /// Mirrors the Term struct of the kernel source.
  struct Term
  {
    std::uint32_t kind{0};
    std::uint32_t type{0};
    std::uint32_t atoms[4]{};
    float parameters[6]{};
  };
  static_assert(sizeof(Term) == 48);

  /// Mirrors the BondedParameters struct of the kernel source.
  struct Parameters
  {
    float alpha{0.0f};
    float twoAlphaOverSqrtPi{0.0f};
    float selfPrefactor{0.0f};
    float coulombFactor{0.0f};
    std::uint32_t useCharge{0};
    std::uint32_t numberOfSlots{0};
    std::uint32_t padding[2]{};
  };
  static_assert(sizeof(Parameters) == 32);

  struct MoleculeInfo
  {
    std::uint32_t firstAtom{0}, numberOfAtoms{0}, atomOffset{0}, padding{0};
  };

  cl_context context{nullptr};
  cl_device_id device{nullptr};
  OpenCLDevice::ProgramHandle program{};
  OpenCLDevice::KernelHandle kernel{};
  OpenCLDevice::MemHandle slotMoleculeBuffer{}, slotOfOriginalBuffer{}, moleculeInfoBuffer{}, atomTermStartBuffer{};
  OpenCLDevice::MemHandle atomTermsBuffer{}, termsBuffer{}, massBuffer{}, parameterBuffer{}, partialBuffer{};
  std::size_t slotCapacity{0}, atomCapacity{0}, partialCapacity{0};

  Parameters parameters{};
  bool parametersChanged{true};
  double alphaValue{0.0}, conversionFactor{0.0}, chargeSquaredSum{0.0};
  std::vector<Term> terms{};
  std::vector<std::uint32_t> atomTermStart{}, atomTerms{};
  std::vector<MoleculeInfo> molecules{};
  std::vector<float> masses{};

  std::vector<std::uint32_t> slotMolecule{}, slotOfOriginal{}, referenceOfSorted{};
  std::vector<float> hostPartials{};
  std::size_t groups{0};
};
