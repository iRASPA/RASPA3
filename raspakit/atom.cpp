module;

module atom;

import archive;
import double3;
import stringutils;

import <istream>;
import <ostream>;
import <sstream>;
import <fstream>;
import <print>;
import <format>;
import <exception>;
import <source_location>;
import <complex>;

//std::string Atom::repr() const
//{
//  std::ostringstream stream;
//
//  std::print(stream, "({}, {}, {}, [{}, {}, {}, {}])\n", position.x, position.y, position.z, 
//      moleculeId, type, componentId, groupId);
//
//  return stream.str();
//}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const Atom &atom)
{
  archive << atom.position;
  archive << atom.velocity;
  archive << atom.gradient;
  archive << atom.charge;
  archive << atom.scalingVDW;
  archive << atom.scalingCoulomb;
  archive << atom.moleculeId;
  archive << atom.type;
  archive << atom.componentId;
  archive << atom.groupId;

  return archive;
};

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, Atom &atom)
{
  archive >> atom.position;
  archive >> atom.velocity;
  archive >> atom.gradient;
  archive >> atom.charge;
  archive >> atom.scalingVDW;
  archive >> atom.scalingCoulomb;
  archive >> atom.moleculeId;
  archive >> atom.type;
  archive >> atom.componentId;
  archive >> atom.groupId;

  return archive;
}


