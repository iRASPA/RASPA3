module;

export module apollonius_accessibility;

import std;

import double3;
import unit_cell;
import apollonius_network;
import pore_accessibility;

// The accessibility classifier of the Voronoi analyses, built on the Apollonius diagram.
//
// It is the same classifier and the same question of every sample point: which cell does it fall in,
// is it inside an atom, and is the nearest node of that cell in a channel or in a pocket. What the
// diagram changes is the answer to the first part. The cells of the radical diagram are cut by the
// power distance, which is not the distance a probe experiences; the cells here are cut by the
// clearance, which is. A point therefore lands in the cell of the atom that is really nearest to it,
// and the nodes it is then tested against are the true maxima of the clearance around it.
export struct ApolloniusAccessibility
{
  PoreAccessibility accessibility;

  // What the diagram came to. Its network has been handed to `accessibility`, which owns it now, so
  // only the counts and the verification are left here, for the header of an output file.
  ApolloniusPoreNetwork diagram;

  // `radii` are the bare atom radii; the atoms are inflated by `probeRadius` before the diagram is
  // built, exactly as for the radical classifier, so that a node's clearance is the room left for the
  // centre of the probe.
  static ApolloniusAccessibility create(const UnitCell& unitCell,
                                        const std::vector<double3>& fractionalPositions,
                                        const std::vector<double>& radii, double probeRadius);
};
