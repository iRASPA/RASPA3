module;

export module spatial_decomposition_device_polymer;

import std;

import double3x3;
import spatial_decomposition_device_context;
import spatial_decomposition_device_polymer_topology;

/// The sums of the polymer-bond overlay over a step (the slots of RunningEnergy it adds to, the strain derivative
/// and the atomic-to-molecular virial correction of its inter-molecular terms).
export struct DevicePolymerResults
{
  double crossLink{0.0}, bond{0.0}, bend{0.0}, torsion{0.0};
  double moleculeMoleculeVDW{0.0}, moleculeMoleculeCharge{0.0}, ewaldExclusion{0.0}, intraVDW{0.0};
  double3x3 strain{};
  double3x3 correction{};
};

/// The box, cutoff, Ewald and Lennard-Jones switching parameters the overlay kernels share with the pair kernel.
export struct DevicePolymerPairSetup
{
  float cell[9]{};
  float inverseCell[9]{};
  float cutOffVDWSquared{0.0f};
  float cutOffChargeSquared{0.0f};
  float alpha{0.0f};
  float coulombFactor{0.0f};
  std::uint32_t useCharge{0};
  std::uint32_t orthorhombic{0};
  std::uint32_t switchMode{0};
  float switchDistanceSquared{0.0f};
  float switchDistance{0.0f};
  float switchInverseWidth{0.0f};
  float switchInverseCutOff3{0.0f};
  float switchA12{1.0f};
  float switchA6{1.0f};
  bool operator==(const DevicePolymerPairSetup&) const = default;
};

/**
 * \brief The polymer-bond overlay of the Polymerization driver on the device: the terms of PolymerOverlayTopology
 * evaluated by the two kernels of polymer_kernel_source.cpp after the bonded kernels, their gradients added to the
 * device forces, their energies, strain derivative and virial correction reduced per work-group and summed on the
 * host in double.
 *
 * The topology is uploaded after every reaction (setTopology), the ramp factors whenever they change (setLambdas),
 * the shared pair parameters when the box or the force field changes (setPairSetup). The kernels read the wrapped
 * slot positions of the pair lists through the slot of every atom (uploaded by DeviceStep at every list build),
 * and the relative positions and molecule centers of the bonded kernels for the virial correction.
 */
export class DevicePolymerOverlay
{
 public:
  DevicePolymerOverlay() = default;
  ~DevicePolymerOverlay() = default;
  DevicePolymerOverlay(const DevicePolymerOverlay&) = delete;
  DevicePolymerOverlay& operator=(const DevicePolymerOverlay&) = delete;
  DevicePolymerOverlay(DevicePolymerOverlay&&) noexcept = default;
  DevicePolymerOverlay& operator=(DevicePolymerOverlay&&) noexcept = default;

  void initialize(DeviceContext& context);
  bool initialized() const { return static_cast<bool>(termKernel); }

  /// Uploads the terms, the participating atoms and their gradient-slot lists (blocking).
  void setTopology(const PolymerOverlayTopology& topology);
  /// Writes the ramp factors (one per ramp of the topology; the write completes before the next launch).
  void setLambdas(std::span<const float> lambdas);
  /// The parameters shared with the pair kernel (written at the next enqueue when they changed).
  void setPairSetup(const DevicePolymerPairSetup& setup);

  std::size_t numberOfTerms() const { return terms; }
  std::size_t numberOfAtoms() const { return atoms; }
  std::size_t numberOfRamps() const { return ramps; }

  /// Enqueues the kernels (after the bonded kernels wrote `force` and `moleculeCenter`).
  void enqueue(DeviceBuffer position, DeviceBuffer slotOfAtom, DeviceBuffer relative, DeviceBuffer moleculeOfAtom,
               DeviceBuffer moleculeCenter, DeviceBuffer force);
  /// Enqueues the read-back of the partial sums (complete with a mark / finish of the context).
  void enqueueRead();
  /// After the read completed: the sums over all terms.
  DevicePolymerResults collect() const;

  std::string status() const;

 private:
  /// Mirrors the PolymerParameters struct of the kernel source.
  struct Parameters
  {
    float cell[9]{};
    float inverseCell[9]{};
    float cutOffVDWSquared{0.0f};
    float cutOffChargeSquared{0.0f};
    float alpha{0.0f};
    float alphaSquared{0.0f};
    float alphaOverSqrtPi{0.0f};
    float twoAlphaOverSqrtPi{0.0f};
    float coulombFactor{0.0f};
    std::uint32_t useCharge{0};
    std::uint32_t orthorhombic{0};
    std::uint32_t numberOfTerms{0};
    std::uint32_t numberOfAtoms{0};
    std::uint32_t switchMode{0};
    float switchDistanceSquared{0.0f};
    float switchDistance{0.0f};
    float switchInverseWidth{0.0f};
    float switchInverseCutOff3{0.0f};
    float switchA12{1.0f};
    float switchA6{1.0f};
  };
  static_assert(sizeof(Parameters) == 144);

  DeviceContext* context{nullptr};
  DeviceKernel termKernel{}, atomKernel{};
  DeviceBufferOwner termsBuffer{}, atomsBuffer{}, refStartBuffer{}, refLinkEndBuffer{}, refsBuffer{};
  DeviceBufferOwner lambdaBuffer{}, termGradientBuffer{}, parameterBuffer{};
  DeviceBufferOwner termPartialBuffer{}, atomPartialBuffer{};  ///< Group partials of the two kernels.
  std::size_t lambdaCapacity{0}, termPartialCapacity{0}, atomPartialCapacity{0};

  Parameters parameters{};
  DevicePolymerPairSetup setup{};
  bool setupApplied{false};
  bool parametersChanged{true};
  std::size_t terms{0}, atoms{0}, ramps{0}, termGroups{0}, atomGroups{0};
  std::array<std::vector<float>, 2> lambdaHost{};
  std::size_t lambdaTurn{0};
  std::vector<float> hostPartials{};  ///< The term partials followed by the atom partials.
};
