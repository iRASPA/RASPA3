module;

export module spatial_decomposition_device_bonded;

import std;

import double3;
import double3x3;
import spatial_decomposition_device_context;
import spatial_decomposition_device_bonded_topology;

/// The per-molecule results of a step (bonded chain), summed over all molecules.
export struct DeviceBondedResults
{
  /// The self energy is evaluated in double on the host (it depends on the charges and alpha only); the device
  /// evaluates the exclusion of the excluded pairs in a reduced form (the pair term minus its small-r limit) and
  /// the host adds the constant that completes it.
  double self{0.0}, exclusion{0.0}, bond{0.0}, bend{0.0}, torsion{0.0}, improperTorsion{0.0};
  double intraVDW{0.0}, intraCoulomb{0.0};
  double3x3 exclusionStrain{};
  double3x3 correction{};
};

/**
 * \brief The per-molecule terms on the device: Ewald self and exclusion corrections, bonds, bends, torsions and
 * improper torsions, intramolecular Lennard-Jones and Coulomb pairs, and the atomic-to-molecular virial
 * correction.
 *
 * Two kernels (the shared bonded_kernel_source.cpp compiled by the DeviceContext): one work-item per term
 * instance evaluates the term once and writes the gradient on each of its atoms to a per-instance slot; one
 * work-item per slot (atom) evaluates the exclusion pairs and the virial correction, gathers the gradients of its
 * terms and adds the sum to the device force of the slot. Both leave per-group partial sums of the energies, the
 * exclusion strain derivative and the virial correction, which the host sums in double. The topology
 * (BondedTopology, built by DeviceStep) is uploaded once; the slot layout (molecule and index of every slot)
 * after every list build. The positions of the terms are the float positions relative to the first atom of the
 * molecule, indexed by the atom's index in the system, packed by the host.
 */
export class DeviceBonded
{
 public:
  DeviceBonded() = default;
  ~DeviceBonded() = default;
  DeviceBonded(const DeviceBonded&) = delete;
  DeviceBonded& operator=(const DeviceBonded&) = delete;
  DeviceBonded(DeviceBonded&&) noexcept = default;
  DeviceBonded& operator=(DeviceBonded&&) noexcept = default;

  void initialize(DeviceContext& context);
  bool initialized() const { return static_cast<bool>(atomKernel); }

  /// Uploads the terms of the components (per atom), the molecule table and the masses (blocking).
  void setTopology(const BondedTopology& topology);
  /// The Ewald parameters (`useCharge`: self and exclusion corrections and the intramolecular Coulomb pairs on)
  /// and the cutoffs of the scaled intramolecular pairs.
  void setParameters(double alpha, double conversionFactor, bool useCharge, double cutOffVDW, double cutOffCharge);
  /// Follows a change of the Ewald alpha.
  void setAlpha(double alpha);

  /// The slot layout of a list build (BondedTopology::layout); enqueues the upload (non-blocking; the span must
  /// stay valid until the stream completed it).
  void setLayout(std::span<const std::uint32_t> slotMolecule);

  /// Enqueues the kernels (after the pair kernel and the mesh interpolation wrote `force`); `relative` holds the
  /// relative position and charge of every atom in the system order.
  void enqueue(DeviceBuffer relative, DeviceBuffer force);
  /// Enqueues the read-back of the partial sums (complete with a mark / finish of the context).
  void enqueueRead();

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
    float cutOffVDWSquared{0.0f};
    float cutOffChargeSquared{0.0f};
    std::uint32_t padding[2]{};
  };
  static_assert(sizeof(Parameters) == 48);

  DeviceContext* context{nullptr};
  DeviceKernel termKernel{}, atomKernel{};
  DeviceBufferOwner slotMoleculeBuffer{}, moleculeInfoBuffer{}, instanceMoleculeBuffer{}, massBuffer{};
  DeviceBufferOwner termsBuffer{}, gradientOffsetBuffer{}, atomGradientStartBuffer{}, atomGradientsBuffer{};
  DeviceBufferOwner exclusionStartBuffer{}, exclusionPartnersBuffer{};
  DeviceBufferOwner termGradientBuffer{}, parameterBuffer{}, partialBuffer{};
  std::size_t slotCapacity{0}, partialCapacity{0};

  Parameters parameters{};
  bool parametersChanged{true};
  double alphaValue{0.0}, conversionFactor{0.0}, chargeSquaredSum{0.0}, exclusionChargeProductSum{0.0};
  std::size_t numberOfMolecules{0}, numberOfTerms{0}, termGroups{0};

  std::vector<float> hostPartials{};
  std::size_t atomGroups{0};
};
