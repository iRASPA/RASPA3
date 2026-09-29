module;

export module spatial_decomposition_pppm;

import std;

import int3;
import double3;
import double3x3;
import simulationbox;

/**
 * \brief Smooth particle-mesh Ewald (SPME) evaluation of the reciprocal-space Coulomb sum.
 *
 * Replaces the direct k-space loop of the Ewald sum (O(N k_max^3)) by charge assignment with cardinal B-splines
 * of order p onto a K_x x K_y x K_z mesh, one real-to-complex FFT, a multiplication with the influence function
 *   G(m) = C (2 pi / V) exp(-k^2 / 4 alpha^2) / k^2 * |b_x(m_x) b_y(m_y) b_z(m_z)|^2,   k = 2 pi (m . inverseCell),
 * one complex-to-real FFT back, and interpolation of the resulting potential mesh with the analytic derivatives
 * of the same B-splines (Essmann et al., J. Chem. Phys. 103, 8577 (1995)). The reciprocal energy is
 *   E = sum_m G(m) |F(Q)(m)|^2,
 * evaluated in the half spectrum of the real transform, and the gradient on atom i is
 *   dE/dr_i = 2 q_i sum_nodes phi(node) d(M_x M_y M_z)/dr_i.
 * The same Ewald alpha as the force field's real-space part is used, so the self and intramolecular exclusion
 * terms of the direct Ewald code apply unchanged (Interactions::addChargeSelfEnergy and
 * addIntraMolecularChargeExclusionGradient), as does the net-charge correction, for which
 * singleIonFourierSum() provides the mesh analogue of the direct sum's single-ion term.
 *
 * The mesh is shared work: every thread spreads the atoms it owns into a private mesh copy (clearMesh, spread),
 * the copies are summed in parallel slabs (reduceMeshes), thread 0 performs the transforms (solve, using FFTW's
 * own threads), and every thread interpolates the gradients of its own atoms (interpolate). Triclinic cells are
 * handled in fractional coordinates. The mesh dimensions are fixed at initialization; the influence function is
 * recomputed when the cell changes (NPT).
 */
export class PPPM
{
 public:
  PPPM() = default;
  ~PPPM();
  PPPM(const PPPM&) = delete;
  PPPM& operator=(const PPPM&) = delete;

  /// Chooses the mesh (smallest 2^a 3^b 5^c sizes at or below `meshSpacing` along each cell vector), plans the
  /// FFTs and computes the B-spline moduli and the influence function for `box`.
  void initialize(const SimulationBox& box, double alpha, double meshSpacing, std::size_t order,
                  std::size_t numberOfThreads, double coulombConversionFactor);

  bool initialized() const { return forwardPlan != nullptr; }
  int3 meshSize() const { return mesh; }
  std::size_t interpolationOrder() const { return order; }
  std::size_t numberOfMeshCopies() const { return meshCopies.size(); }

  /// Recomputes the influence function when the cell or alpha differ from what it was built for. Thread 0 only.
  void updateBox(const SimulationBox& box, double alphaValue);

  void clearMesh(std::size_t thread);

  /// Spreads the charges of the given atoms (sorted-order indices into the SoA arrays) into mesh copy `thread`.
  void spread(std::size_t thread, std::span<const std::uint32_t> atoms, const double* x, const double* y,
              const double* z, const double* charge, const double* scalingCoulomb);

  /// Sums all mesh copies into copy 0 for the x-slab of this thread out of `numberOfThreads`.
  void reduceMeshes(std::size_t thread, std::size_t numberOfThreads);

  /// Forward transform of copy 0, influence-function multiplication and inverse transform into the potential
  /// mesh. Returns the reciprocal energy (energy units of the force field). With `withVirial` the strain
  /// derivative of the reciprocal energy is accumulated as well (see reciprocalStrainDerivative). Thread 0 only.
  double solve(bool withVirial);

  /// Strain derivative sum_m G(m) |F(Q)(m)|^2 [I - 2 (1/k^2 + 1/(4 alpha^2)) k k^T] of the last solve(true).
  const double3x3& reciprocalStrainTensor() const { return reciprocalStrain; }

  /// The same tensor for a unit point charge without the B-spline moduli, sum_m bare(m) [I - 2 (...) k k^T];
  /// multiplied by the squared net charge it is the strain response of the net-charge correction.
  const double3x3& singleIonStrainTensor() const { return ionStrain; }

  /// Adds the reciprocal gradient dE/dr of the given atoms to fx/fy/fz (sorted-order indexing).
  void interpolate(std::span<const std::uint32_t> atoms, const double* x, const double* y, const double* z,
                   const double* charge, const double* scalingCoulomb, double* fx, double* fy, double* fz) const;

  /// C sum_{m != 0} (2 pi / V) exp(-k^2 / 4 alpha^2) / k^2 over the mesh wave vectors (all m, without the B-spline
  /// moduli): the single-ion term of the net-charge correction.
  double singleIonFourierSum() const { return singleIonSum; }

  std::string status() const;

  /// Smallest integer >= n whose only prime factors are 2, 3 and 5 (FFT-friendly).
  static std::size_t nextFFTFriendly(std::size_t n);

  /// Cardinal B-spline weights M_p(w + j), j = 0..p-1, and their derivatives for the fractional offset w in [0, 1).
  static void bsplineWeights(std::size_t order, double w, std::span<double> weights, std::span<double> derivatives);

 private:
  int3 mesh{0, 0, 0};
  std::size_t order{5};
  double alpha{0.0};
  double conversionFactor{1.0};
  double3x3 inverseCellAtBuild{};
  double3x3 inverseCell{};
  double volume{0.0};

  std::vector<double*> meshCopies{};  ///< fftw_malloc'ed real meshes, one per thread; copy 0 is the FFT input.
  double* potential{nullptr};         ///< Real output of the inverse transform.
  void* spectrum{nullptr};            ///< fftw_complex half spectrum.
  void* forwardPlan{nullptr};
  void* backwardPlan{nullptr};

  std::vector<double> influence{};  ///< G(m) on the half spectrum.
  std::vector<double> bsplineModulusX{}, bsplineModulusY{}, bsplineModulusZ{};
  double singleIonSum{0.0};
  double3x3 ionStrain{};
  double3x3 reciprocalStrain{};

  std::size_t realSize() const;
  std::size_t complexSize() const;
  void computeBsplineModuli();
  void computeInfluenceFunction(const SimulationBox& box);
  void release();
};
