module;

module blocking_spheres;

import std;

void writeBlockFile(const std::string& frameworkName, const std::vector<BlockingSphere>& spheres)
{
  std::ofstream blockFile;
  blockFile.open(frameworkName + ".block");
  std::print(blockFile, "{}\n", spheres.size());
  for (const BlockingSphere& sphere : spheres)
  {
    std::print(blockFile, "{} {} {} {}\n", sphere.centerFractional.x, sphere.centerFractional.y,
               sphere.centerFractional.z, sphere.radius);
  }
  blockFile.close();
}
