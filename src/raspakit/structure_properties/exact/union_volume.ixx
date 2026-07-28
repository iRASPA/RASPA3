module;

export module exact_union_volume;

import std;

import double2;
import voronoi_accessibility;
import exact_surface_patches;

// The volume of the union of the probe-inflated atoms, measured rather than sampled.
//
// The void fraction of a framework is what is left of the cell once that union is taken away, so this
// is the geometric void fraction with the sampling removed from it. Two facts make it a finite sum of
// closed-form pieces.
//
// The first is that a cell of the diagram sees only its own atom. If a point lies in the cell of atom i
// and inside any inflated atom at all, then the clearance there is at most zero and atom i is the atom
// that attains it, so the point lies inside atom i as well. The union is therefore cut into one piece
// per atom, the part of that atom inside its own cell, with nothing missed and nothing counted twice.
//
// The second is the divergence theorem applied to one such piece with the field x - p_i, whose
// divergence is three. The piece is bounded by part of the sphere of atom i and by flat faces, and on
// each the field's normal component is a constant: the radius on the sphere, and the distance from the
// centre to the plane on a face. So
//
//   3 |B_i cap C_i|  =  R_i |S_i|  +  sum_j h_ij |F_ij| ,
//
// with S_i the exposed patch of the sphere and F_ij the part of the face on the plane between i and j
// that lies inside both the ball and the cell. The faces are flat for the radical, or power, cells
// rather than the Apollonius ones, so it is those that are used here; the two agree on which atom a
// point of the surface belongs to, so the patches S_i are the same either way and are the ones the
// surface-area code already measures. A volume needs no more than a partition in which each cell sees
// its own atom, which the radical cells provide and which is why they serve here although they are the
// wrong construction for the pore diameters.
//
// Each face is then a disc cut by half planes, whose area is elementary, so the only quadrature
// anywhere in the result is the latitude integral behind |S_i|.
export double unionOfBallsVolume(const VoronoiAccessibility& accessibility, std::size_t subdivisions = 1);

// The same, from a sweep of the patches already done. The spheres enter only through their radius
// weighted area, so a caller that has measured the surface, and had to classify it for other reasons,
// need not sweep a second time.
export double unionOfBallsVolume(const VoronoiAccessibility& accessibility, const ExactSurfaceAreaSample& patches);

// The area of a disc of radius `radius` about the origin of its own plane, cut back by the half planes
// `dot(normals[k], p) <= offsets[k]`. It is the face above, and it is exported because the same shape turns
// up away from any cell of any diagram: the pieces of the solvent excluded surface that lie on a plane are
// discs cut by half planes too. Exact for any number of lines, including none.
export double clippedDiscArea(double radius, const std::vector<double2>& normals, const std::vector<double>& offsets);
