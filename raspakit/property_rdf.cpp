module;

module property_rdf;

import <cstddef>;
import <string>;
import <iostream>;
import <sstream>;
import <tuple>;
import <vector>;
import <algorithm>;

import double3;
import atom;
import simulationbox;

void PropertyRadialDistributionFunction::sample(const SimulationBox &simulationBox, std::span<Atom> frameworkAtoms, std::span<Atom> moleculeAtoms, size_t block)
{
  double3 dr, posA, posB, f;
  double rr, r;

  if(moleculeAtoms.empty()) return;

  for (std::span<const Atom>::iterator it1 = frameworkAtoms.begin(); it1 != frameworkAtoms.end(); ++it1)
  {
    posA = it1->position;
    size_t typeA = static_cast<size_t>(it1->type);
    for (std::span<const Atom>::iterator it2 = moleculeAtoms.begin(); it2 != moleculeAtoms.end(); ++it2)
    {
      posB = it2->position;
      size_t typeB = static_cast<size_t>(it2->type);

      dr = posA - posB;
      dr = simulationBox.applyPeriodicBoundaryConditions(dr);
      rr = double3::dot(dr, dr);
      r = std::sqrt(rr);

      if (typeA < typeB) std::swap(typeA, typeB);
      size_t index_pseudo_atoms = typeA * (typeA + 1) / 2 + typeB;
      size_t bin = static_cast<size_t>(r / deltaR);
      size_t index = bin + numberOfBins * (index_pseudo_atoms + block * numberOfPseudoAtomsSymmetricMatrix);
      sumProperty.at(index) += 1.0;
    }
  }

  for (std::span<const Atom>::iterator it1 = moleculeAtoms.begin(); it1 != moleculeAtoms.end() - 1; ++it1)
  {
    posA = it1->position;
    size_t molA = static_cast<size_t>(it1->moleculeId);
    size_t compA = static_cast<size_t>(it1->componentId);
    size_t typeA = static_cast<size_t>(it1->type);
    for (std::span<const Atom>::iterator it2 = it1 + 1; it2 != moleculeAtoms.end(); ++it2)
    {
      size_t molB = static_cast<size_t>(it2->moleculeId);
      size_t compB = static_cast<size_t>(it2->componentId);

      // skip interactions within the same molecule
      if (!((compA == compB) && (molA == molB)))
      {
        posB = it2->position;
        size_t typeB = static_cast<size_t>(it2->type);

        dr = posA - posB;
        dr = simulationBox.applyPeriodicBoundaryConditions(dr);
        rr = double3::dot(dr, dr);
        r = std::sqrt(rr);

        if (typeA < typeB) std::swap(typeA, typeB);
        size_t index_pseudo_atoms = typeA * (typeA + 1) / 2 + typeB;
        size_t bin = static_cast<size_t>(r / deltaR);
        size_t index = bin + numberOfBins * (index_pseudo_atoms + block * numberOfPseudoAtomsSymmetricMatrix);
        sumProperty.at(index) += 1.0;
      }
    }
  }

  numberOfCounts++;
}

void PropertyRadialDistributionFunction::writeOutput([[maybe_unused]] size_t systemId)
{
}
