module;

module apollonius_network;

import std;

import int3;
import double3;
import double3x3;
import unit_cell;
import skapolloniusdiagram;
import voronoi_network;

bool ApolloniusPoreNetwork::networkIsComplete() const
{
  return verification.vertexCount > 0 && verification.coincidentVertices == 0 && verification.unpairedTriples == 0 &&
         verification.overpairedTriples == 0 && verification.unsampledArcs == 0 &&
         verification.verticesOfFullValence + verification.truncatedTriples >= verification.vertexCount;
}

void ApolloniusPoreNetwork::writeHeader(std::ostream& stream) const
{
  std::print(stream, "# Number of Apollonius vertices: {}\n", numberOfVertices);
  std::print(stream, "# Number of Apollonius arcs: {}\n", numberOfArcs);
  if (numberOfRings > 0)
  {
    std::print(stream, "# Ring passages left out of the network (closed arcs without a vertex): {}\n", numberOfRings);
  }

  if (!networkIsComplete())
  {
    std::print(stream, "# WARNING: the network did not verify as consistent, so peaks or passages may be\n");
    std::print(stream, "#          missing and what is read from it below is a lower bound\n");
    std::print(stream, "#          unpaired triples: {}, overpaired: {}, coincident vertices: {}\n",
               verification.unpairedTriples, verification.overpairedTriples, verification.coincidentVertices);
    std::print(stream, "#          vertices: {}, of full valence: {}, arcs clipped by the free region: {}\n",
               verification.vertexCount, verification.verticesOfFullValence, verification.truncatedTriples);
    if (verification.unsampledArcs > 0)
    {
      std::print(stream, "#          arcs whose narrowest point could not be sampled: {}\n",
                 verification.unsampledArcs);
    }
  }
  else if (!verification.isComplete())
  {
    // The nodes and arcs are sound; what is left is degeneracy, which a symmetric framework has plenty
    // of and which the analyses never look at. Reported all the same, since it is a property of the
    // structure worth knowing.
    std::print(stream, "# Note: the input is degenerate: {} vertices touch more than four sites, {} branch\n",
               verification.degenerateVertices, verification.ambiguousBranches);
    std::print(stream, "#       directions came out a tie and {} bisector patches do not close into a cycle.\n",
               verification.unclosedFaces);
    std::print(stream, "#       The nodes and arcs read below are consistent.\n");
  }
}

ApolloniusPoreNetwork ApolloniusPoreNetwork::create(const UnitCell& unitCell,
                                                    const std::vector<double3>& fractionalPositions,
                                                    const std::vector<double>& radii)
{
  ApolloniusPoreNetwork result;
  VoronoiNetwork& network = result.network;

  const double3x3 cell = unitCell.cell;
  const double3x3 inverseCell = unitCell.inverseCell;
  const std::size_t numberOfAtoms = fractionalPositions.size();

  network.unitCell = unitCell;
  network.atomRadii = radii;
  network.atomPositionsFractional.reserve(numberOfAtoms);
  for (const double3& fractional : fractionalPositions)
  {
    network.atomPositionsFractional.push_back(double3::fract(fractional));
  }
  network.atomNodeVectors.assign(numberOfAtoms, {});

  std::vector<double3> atomPositions(numberOfAtoms);
  for (std::size_t i = 0; i < numberOfAtoms; ++i) atomPositions[i] = cell * network.atomPositionsFractional[i];

  // The free-space diagram is the one pore geometry wants: it stops at the surface of the union of the
  // atoms and never explores their interiors, where a probe cannot go anyway.
  SKApolloniusDiagram diagram =
      SKApolloniusDiagram::create(cell, network.atomPositionsFractional, radii, 1, SKApolloniusRegion::FreeSpace);

  result.numberOfVertices = diagram.vertices.size();
  result.verification = diagram.verification;

  network.nodes.reserve(diagram.vertices.size());
  for (const SKApolloniusVertex& vertex : diagram.vertices)
  {
    std::size_t nodeIndex = network.nodes.size();

    VoronoiNode node;
    node.position = vertex.position;
    node.fractional = double3::fract(inverseCell * vertex.position);
    node.radius = vertex.radius;
    // The vertex is the local maximum of the clearance, so it is its own peak and needs no ascent.
    node.maximalRadius = vertex.radius;
    node.maximalPosition = vertex.position;

    for (std::size_t s = 0; s < vertex.siteIndices.size(); ++s)
    {
      std::size_t atomIndex = vertex.siteIndices[s];
      const int3& image = vertex.siteImages[s];

      if (std::find(node.atomIndices.begin(), node.atomIndices.end(), atomIndex) == node.atomIndices.end())
      {
        node.atomIndices.push_back(atomIndex);
      }

      // The tangency is between the atom at `image` and this vertex; seen from the atom in the home
      // cell it is the matching image of the vertex that is tangent, which is what the analyses that
      // walk from an atom to the vertices of its own cell expect.
      double3 shift = cell * double3(image.x, image.y, image.z);
      network.atomNodeVectors[atomIndex].push_back({nodeIndex, vertex.position - shift - atomPositions[atomIndex]});
    }

    network.nodes.push_back(std::move(node));
  }

  network.edges.reserve(2 * diagram.edges.size());
  for (const SKApolloniusEdge& edge : diagram.edges)
  {
    if (edge.isLoop)
    {
      ++result.numberOfRings;
      continue;
    }

    ++result.numberOfArcs;

    // The arc's bottleneck is given in the frame in which its `from` vertex sits in the home cell, so
    // the reverse edge, whose `from` is the other end, sees it one image away.
    double3 reverseShift = cell * double3(edge.toImage.x, edge.toImage.y, edge.toImage.z);

    network.edges.push_back(VoronoiEdge{edge.from, edge.to, edge.toImage, edge.bottleneckRadius, edge.length,
                                        edge.bottleneckPosition, edge.bottleneckDirection, true});
    network.edges.push_back(VoronoiEdge{edge.to, edge.from, -edge.toImage, edge.bottleneckRadius, edge.length,
                                        edge.bottleneckPosition - reverseShift, -edge.bottleneckDirection, true});
  }

  return result;
}
