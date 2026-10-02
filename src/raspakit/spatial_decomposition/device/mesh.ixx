module;

export module spatial_decomposition_device_mesh;

import std;

import double3;
import double3x3;
import int3;
import simulationbox;
import spatial_decomposition_device_context;

/**
 * \brief The particle-mesh Ewald sum (PPPM) on the device, over the slots of the pair kernel.
 *
 * The device transcription of PPPM in single precision: charge spreading with cardinal B-splines into a
 * fixed-point mesh (integer atomics: not every device has floating-point atomics, and the integer sum is order
 * independent), a real-to-complex 3D FFT of the mesh (own Stockham radix 2/3/4/5 kernels over tiles of lines in
 * local memory; the real z lines are packed into half-length complex transforms, so only the half spectrum
 * kz = 0..Kz/2 is computed; the mesh sizes are 2^a 3^b 5^c by construction and Kz is made even), the influence
 * function computed on the fly per wave vector (so a cell change costs nothing), the inverse FFT and the
 * interpolation of the gradients, which are added to the device force of the slots. The reciprocal energy, its
 * strain derivative and the single-ion sums of the net-charge correction are reduced per work-group on the device
 * and summed in double on the host.
 *
 * The kernels are the shared mesh_kernel_source.cpp, compiled by the DeviceContext of DeviceStep, which owns this
 * object and provides the position / force buffers; the chain is enqueued after the pair kernel of the step (the
 * interpolation adds to the pair forces).
 */
export class DeviceMesh
{
 public:
  DeviceMesh() = default;
  ~DeviceMesh() = default;
  DeviceMesh(const DeviceMesh&) = delete;
  DeviceMesh& operator=(const DeviceMesh&) = delete;
  DeviceMesh(DeviceMesh&&) noexcept = default;
  DeviceMesh& operator=(DeviceMesh&&) noexcept = default;

  /// Builds the kernels for the given B-spline order (3 to 7).
  void initialize(DeviceContext& context, std::size_t order);
  bool initialized() const { return static_cast<bool>(spreadKernel); }

  /// Fixes the mesh, alpha and the Coulomb conversion factor: allocates the mesh buffers, plans the FFTs
  /// (twiddle tables) and computes the B-spline moduli. An odd z size is raised to the next even FFT-friendly size
  /// (see meshSize()). Throws when an axis exceeds the local memory of the device.
  void setup(int3 mesh, double alpha, double conversionFactor);

  /// Number of slots (positions) to spread and interpolate.
  void setSlots(std::size_t slots);
  /// Follows a change of the Ewald alpha.
  void setAlpha(double alphaValue);
  /// Cell-dependent parameters (inverse cell, volume); written to the device at the next enqueue when changed.
  void updateBox(const SimulationBox& box);

  /// Enqueues the step: spreading, forward FFT, influence function, inverse FFT and (with `interpolate`) the
  /// interpolation of the gradients into `force`. Non-blocking.
  void enqueue(DeviceBuffer position, DeviceBuffer force, bool interpolate);
  /// Enqueues the read-back of the per-group partial sums (complete with a mark / finish of the context).
  void enqueueRead();
  /// Times the stages of one step synchronously (the fastest of a few runs each; the mesh and `force` are left
  /// modified) and keeps the result for the status report.
  void profile(DeviceBuffer position, DeviceBuffer force);
  /// After the read completed: the reciprocal energy, its strain derivative, the single-ion sum and tensor.
  void collect(double& energy, double3x3& strain, double& singleIonSum, double3x3& singleIonStrain) const;

  int3 meshSize() const { return mesh; }
  std::size_t interpolationOrder() const { return order; }
  std::string status() const;

 private:
  /// Mirrors the MeshParameters struct of the kernel source.
  struct Parameters
  {
    float inverseCell[9]{};
    float scale{0.0f};
    float inverseScale{0.0f};
    float prefactor{0.0f};
    float alphaFactor{0.0f};
    float inverseFourAlphaSquared{0.0f};
    std::uint32_t meshX{0}, meshY{0}, meshZ{0};
    std::uint32_t numberOfSlots{0};
    std::uint32_t padding[6]{};
  };
  static_assert(sizeof(Parameters) == 96);

  /// The 1D transforms along one axis: the line geometry and the Stockham stage plan. For the real z axis N is the
  /// half length Kz / 2 and the local buffers hold `localLength` = N + 1 values per line (the half spectrum).
  struct AxisPlan
  {
    std::uint32_t N{0}, localLength{0}, radixCode{0}, stages{0};
    std::uint32_t axisStride{0}, lineStride{0}, innerCount{0}, outerStride{0}, tile{1}, tileShift{0};
    std::uint32_t tilesPerOuter{1};  // the tile (lines per work-group) is 2^tileShift
    std::size_t groups{0}, groupSize{64};
    DeviceBufferOwner twiddle{};
  };

  DeviceContext* context{nullptr};
  DeviceKernel spreadKernel{}, realForwardKernel{}, realBackwardKernel{}, fftKernel{}, influenceKernel{},
      interpolateKernel{};
  /// Fixed-point charge mesh, the half spectrum (complex), the real potential mesh.
  DeviceBufferOwner meshBuffer{}, dataBuffer{}, potentialBuffer{}, parameterBuffer{}, partialBuffer{};
  DeviceBufferOwner moduliX{}, moduliY{}, moduliZ{};
  /// exp(-2 pi i k / Kz), k = 0..Kz/2: the twiddles of the real-to-complex unpacking.
  DeviceBufferOwner halfTwiddle{};
  AxisPlan planX{}, planY{}, planZ{};
  std::size_t localMemory{32768};
  std::size_t maxGroupSize{256};

  int3 mesh{0, 0, 0};
  std::size_t order{5};
  double alpha{0.0};
  double conversionFactor{1.0};
  Parameters parameters{};
  bool parametersChanged{true};
  std::size_t influenceGroups{0};
  std::uint32_t lines{0};  ///< z lines of the real transforms (Kx Ky)
  std::vector<float> hostPartials{};
  /// Sampled stage times (seconds): spreading, forward FFTs (with the conversion), influence, backward FFTs,
  /// interpolation.
  std::array<double, 5> stageTimes{};
  bool profiled{false};

  void writeParameters(bool blocking);
  void enqueueSpread(DeviceBuffer position);
  void enqueueForwardTransforms();
  void enqueueInfluence();
  void enqueueBackwardTransforms();
  void enqueueInterpolate(DeviceBuffer position, DeviceBuffer force);

  void planAxis(AxisPlan& plan, std::uint32_t N, std::uint32_t localLength, std::uint32_t axisStride,
                std::uint32_t lineStride, std::uint32_t innerCount, std::uint32_t outerStride,
                std::uint32_t outerCount);
  void enqueueTransform(const AxisPlan& plan, float sign);
  void enqueueRealTransform(bool forward);
  std::size_t meshPoints() const
  {
    return static_cast<std::size_t>(mesh.x) * static_cast<std::size_t>(mesh.y) * static_cast<std::size_t>(mesh.z);
  }
  /// Points of the half spectrum, Kx Ky (Kz / 2 + 1).
  std::size_t spectrumPoints() const
  {
    return static_cast<std::size_t>(mesh.x) * static_cast<std::size_t>(mesh.y) *
           (static_cast<std::size_t>(mesh.z) / 2 + 1);
  }
};
