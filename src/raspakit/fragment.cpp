module;

module fragment;

import std;

import archive;
import double3;
import double3x3;
import simd_quatd;
import atom;
import molecule;

void Fragment::computeRigidProperties(std::span<const double3> referencePositions, std::span<const double> masses)
{
  mass = 0.0;
  centerOfMassReferencePosition = double3{};
  inertiaVector = double3{};
  inverseInertiaVector = double3{};
  bodyFixedPositions.clear();
  rotationalDegreesOfFreedom = 0;
  shapeType = 2;  // point

  // Cache the per-atom masses so the fragment is self-contained (used by 'deriveState').
  atomMasses.clear();
  atomMasses.reserve(atoms.size());
  for (std::size_t atom : atoms) atomMasses.push_back(masses[atom]);

  // Fragment center of mass in the molecule reference frame ('referencePositions' is the reference
  // geometry the CBMC growth places rigid fragments from, already in the whole-molecule frame).
  double3 com{};
  double total_mass{};
  for (std::size_t atom : atoms)
  {
    com += masses[atom] * referencePositions[atom];
    total_mass += masses[atom];
  }
  com = com / total_mass;
  mass = total_mass;
  centerOfMassReferencePosition = com;

  // A single-atom fragment is a point mass: no orientation, no rotational degrees of freedom.
  if (atoms.size() == 1)
  {
    bodyFixedPositions.push_back(double3{});
    return;
  }

  // Inertia tensor about the fragment center of mass.
  double3x3 inertiaTensor{};
  for (std::size_t atom : atoms)
  {
    double m = masses[atom];
    double3 dr = referencePositions[atom] - com;
    inertiaTensor.ax += m * dr.x * dr.x;
    inertiaTensor.bx += m * dr.y * dr.x;
    inertiaTensor.cx += m * dr.z * dr.x;
    inertiaTensor.ay += m * dr.x * dr.y;
    inertiaTensor.by += m * dr.y * dr.y;
    inertiaTensor.cy += m * dr.z * dr.y;
    inertiaTensor.az += m * dr.x * dr.z;
    inertiaTensor.bz += m * dr.y * dr.z;
    inertiaTensor.cz += m * dr.z * dr.z;
  }

  // The body frame is the frame in which the inertia tensor is diagonal.
  double3 eigenvalues{};
  double3x3 eigenvectors{};
  inertiaTensor.EigenSystemSymmetric(eigenvalues, eigenvectors);

  bodyFixedPositions.reserve(atoms.size());
  for (std::size_t atom : atoms)
  {
    double3 dr = referencePositions[atom] - com;
    double3 pos = eigenvectors.transpose() * dr;
    if (std::abs(pos.x) < 1e-8) pos.x = 0.0;
    if (std::abs(pos.y) < 1e-8) pos.y = 0.0;
    if (std::abs(pos.z) < 1e-8) pos.z = 0.0;
    bodyFixedPositions.push_back(pos);

    inertiaVector.x += masses[atom] * (pos.y * pos.y + pos.z * pos.z);
    inertiaVector.y += masses[atom] * (pos.x * pos.x + pos.z * pos.z);
    inertiaVector.z += masses[atom] * (pos.x * pos.x + pos.y * pos.y);
  }

  double rotlim = std::max(1.0e-2, inertiaVector.x + inertiaVector.y + inertiaVector.z) * 1.0e-5;
  inverseInertiaVector.x = (inertiaVector.x < rotlim) ? 0.0 : 1.0 / inertiaVector.x;
  inverseInertiaVector.y = (inertiaVector.y < rotlim) ? 0.0 : 1.0 / inertiaVector.y;
  inverseInertiaVector.z = (inertiaVector.z < rotlim) ? 0.0 : 1.0 / inertiaVector.z;

  double rotall = inertiaVector.x + inertiaVector.y + inertiaVector.z > 1.0e-5
                      ? inertiaVector.x + inertiaVector.y + inertiaVector.z
                      : 1.0;
  std::size_t zeroModes = 0;
  if (inertiaVector.x / rotall < 1.0e-5) ++zeroModes;
  if (inertiaVector.y / rotall < 1.0e-5) ++zeroModes;
  if (inertiaVector.z / rotall < 1.0e-5) ++zeroModes;

  // zeroModes: 0 -> nonlinear (3 rot DOF), 1 -> linear (2 rot DOF), 3 -> single point (0 rot DOF).
  rotationalDegreesOfFreedom = 3 - std::min<std::size_t>(zeroModes, 3);
  shapeType = (zeroModes == 0) ? 0 : (zeroModes >= 3 ? 2 : 1);
}

void Fragment::regenerateAtoms(const GroupState &state, std::span<Atom> moleculeAtoms) const
{
  double3x3 rotation = double3x3::buildRotationMatrixInverse(state.orientation);
  for (std::size_t k = 0; k != atoms.size(); ++k)
  {
    moleculeAtoms[atoms[k]].position = state.centerOfMassPosition + rotation * bodyFixedPositions[k];
  }
}

GroupState Fragment::deriveState(std::span<const Atom> moleculeAtoms) const
{
  GroupState state{};

  double3 com{};
  for (std::size_t k = 0; k != atoms.size(); ++k)
  {
    com += atomMasses[k] * moleculeAtoms[atoms[k]].position;
  }
  state.centerOfMassPosition = com / mass;

  if (atoms.size() == 1)
  {
    state.orientation = simd_quatd(0.0, 0.0, 0.0, 1.0);
    return state;
  }

  // Orientation by an orthogonal Procrustes fit of the body-fixed positions onto the laboratory
  // positions: 'computeRotationMatrix' returns R with p_k - com = R * bodyFixed_k. The stored
  // quaternion follows the fully rigid molecule convention p = com + buildRotationMatrixInverse(q) *
  // bodyFixed, so buildRotationMatrixInverse(q) = R and buildRotationMatrix(q) = R^T.
  std::vector<double3> bodyPositions(bodyFixedPositions.begin(), bodyFixedPositions.end());
  std::vector<double3> labPositions(atoms.size());
  for (std::size_t k = 0; k != atoms.size(); ++k)
  {
    labPositions[k] = moleculeAtoms[atoms[k]].position;
  }

  double3x3 rotation = double3x3::computeRotationMatrix(double3(0.0, 0.0, 0.0), bodyPositions,
                                                        state.centerOfMassPosition, labPositions);
  double3x3 rotationTranspose = rotation.transpose();
  state.orientation = rotationTranspose.quaternion().normalized();
  return state;
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const Fragment &f)
{
  archive << f.versionNumber;
  archive << f.atoms;
  archive << f.atomMasses;
  archive << f.mass;
  archive << f.centerOfMassReferencePosition;
  archive << f.inertiaVector;
  archive << f.inverseInertiaVector;
  archive << f.bodyFixedPositions;
  archive << f.rotationalDegreesOfFreedom;
  archive << f.shapeType;
  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, Fragment &f)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > f.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'Fragment' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }
  archive >> f.atoms;
  archive >> f.atomMasses;
  archive >> f.mass;
  archive >> f.centerOfMassReferencePosition;
  archive >> f.inertiaVector;
  archive >> f.inverseInertiaVector;
  archive >> f.bodyFixedPositions;
  archive >> f.rotationalDegreesOfFreedom;
  archive >> f.shapeType;
  return archive;
}
