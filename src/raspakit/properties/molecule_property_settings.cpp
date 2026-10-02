module;

module molecule_property_settings;

import std;

import archive;

namespace
{
void checkVersion(std::uint64_t read, std::uint64_t current, const char* name)
{
  if (read != current)
  {
    const std::source_location& location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading '{}' at line {} in file {}\n", name, location.line(),
                                         location.file_name()));
  }
}
}  // namespace

Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const MoleculePropertiesSettings& s)
{
  archive << s.versionNumber;
  archive << s.sampleEvery;
  archive << s.writeEvery;
  archive << s.numberOfBins;
  archive << s.bondRange;
  archive << s.endToEndRange;
  archive << s.endToEndBinWidth;
  return archive;
}

Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, MoleculePropertiesSettings& s)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  checkVersion(versionNumber, s.versionNumber, "MoleculePropertiesSettings");
  archive >> s.sampleEvery;
  archive >> s.writeEvery;
  archive >> s.numberOfBins;
  archive >> s.bondRange;
  archive >> s.endToEndRange;
  archive >> s.endToEndBinWidth;
  return archive;
}

Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const MoleculeShapeSettings& s)
{
  archive << s.versionNumber;
  archive << s.sampleEvery;
  archive << s.writeEvery;
  archive << s.numberOfBins;
  archive << s.massWeighted;
  archive << s.radiusOfGyrationRange;
  archive << s.radiusOfGyrationBinWidth;
  return archive;
}

Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, MoleculeShapeSettings& s)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  checkVersion(versionNumber, s.versionNumber, "MoleculeShapeSettings");
  archive >> s.sampleEvery;
  archive >> s.writeEvery;
  archive >> s.numberOfBins;
  archive >> s.massWeighted;
  archive >> s.radiusOfGyrationRange;
  archive >> s.radiusOfGyrationBinWidth;
  return archive;
}

Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const MoleculeBackboneSettings& s)
{
  archive << s.versionNumber;
  archive << s.sampleEvery;
  archive << s.writeEvery;
  archive << s.numberOfWaveVectors;
  archive << s.waveVectorLowerLimit;
  archive << s.waveVectorUpperLimit;
  return archive;
}

Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, MoleculeBackboneSettings& s)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  checkVersion(versionNumber, s.versionNumber, "MoleculeBackboneSettings");
  archive >> s.sampleEvery;
  archive >> s.writeEvery;
  archive >> s.numberOfWaveVectors;
  archive >> s.waveVectorLowerLimit;
  archive >> s.waveVectorUpperLimit;
  return archive;
}

Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const EndToEndACFSettings& s)
{
  archive << s.versionNumber;
  archive << s.sampleEvery;
  archive << s.writeEvery;
  archive << s.numberOfBlockElements;
  return archive;
}

Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, EndToEndACFSettings& s)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  checkVersion(versionNumber, s.versionNumber, "EndToEndACFSettings");
  archive >> s.sampleEvery;
  archive >> s.writeEvery;
  archive >> s.numberOfBlockElements;
  return archive;
}
