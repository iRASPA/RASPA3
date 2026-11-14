module;

module skposcarparser;

import std;

import skelement;
import skatom;
import double3;
import double3x3;
import scanner;
import characterset;
import skparser;
import skstructure;
import skcell;

std::vector<std::string> splitTokens(const std::string& line)
{
  std::vector<std::string> tokens{};
  std::size_t position = 0;

  while (position < line.size())
  {
    while (position < line.size() && std::isspace(static_cast<unsigned char>(line[position])))
    {
      ++position;
    }

    if (position >= line.size())
    {
      break;
    }

    const std::size_t start = position;
    while (position < line.size() && !std::isspace(static_cast<unsigned char>(line[position])))
    {
      ++position;
    }

    tokens.push_back(line.substr(start, position - start));
  }

  return tokens;
}

// '#' and '!' both start a comment (VASP input uses either)
void stripInlineComment(std::string& line)
{
  const std::size_t commentPosition = line.find_first_of("#!");
  if (commentPosition != std::string::npos)
  {
    line.erase(commentPosition);
  }
}

std::string trim(const std::string& string)
{
  std::size_t start = 0;
  while (start < string.size() && std::isspace(static_cast<unsigned char>(string[start])))
  {
    ++start;
  }

  std::size_t end = string.size();
  while (end > start && std::isspace(static_cast<unsigned char>(string[end - 1])))
  {
    --end;
  }

  return string.substr(start, end - start);
}

bool equalsIgnoreCase(const std::string& left, const std::string& right)
{
  return left.size() == right.size() &&
         std::equal(left.begin(), left.end(), right.begin(),
                    [](char a, char b) { return std::tolower(static_cast<unsigned char>(a)) ==
                                                std::tolower(static_cast<unsigned char>(b)); });
}

bool startsWithIgnoreCase(const std::string& line, std::string_view prefix)
{
  if (line.size() < prefix.size())
  {
    return false;
  }

  return std::equal(prefix.begin(), prefix.end(), line.begin(),
                    [](char a, char b) { return std::tolower(static_cast<unsigned char>(a)) ==
                                                std::tolower(static_cast<unsigned char>(b)); });
}

bool isIntegerToken(const std::string& token)
{
  if (token.empty())
  {
    return false;
  }

  std::size_t index = 0;
  if (token[0] == '+' || token[0] == '-')
  {
    if (token.size() == 1)
    {
      return false;
    }
    index = 1;
  }

  return std::all_of(token.begin() + static_cast<std::ptrdiff_t>(index), token.end(),
                     [](char character) { return std::isdigit(static_cast<unsigned char>(character)); });
}

bool isAllIntegerTokens(const std::vector<std::string>& tokens)
{
  return !tokens.empty() &&
         std::all_of(tokens.begin(), tokens.end(),
                     [](const std::string& token) { return isIntegerToken(token); });
}

// The next line with content; blank lines are skipped by the scanner, and a line holding only a comment is skipped
// here (it used to come back empty, which the callers take for the end of the file).
std::string readLine(Scanner& scanner)
{
  std::string line;
  while (scanner.scanUpToCharacters(CharacterSet::newlineCharacterSet(), line))
  {
    stripInlineComment(line);
    line = trim(line);
    if (!line.empty()) return line;
  }
  return {};
}

std::string readRequiredLine(Scanner& scanner, const std::string& context)
{
  const std::string line = readLine(scanner);
  if (line.empty())
  {
    throw std::runtime_error(std::format("POSCAR parse error: {}", context));
  }

  return line;
}

double readScaleFactor(Scanner& scanner)
{
  const std::string scaleLine = readRequiredLine(scanner, "missing universal scaling factor");
  const std::vector<std::string> tokens = splitTokens(scaleLine);
  if (tokens.empty())
  {
    throw std::runtime_error("POSCAR parse error: missing universal scaling factor");
  }

  return std::stod(tokens.front());
}

double3x3 readLatticeVectors(Scanner& scanner, double scale)
{
  double3x3 unitCell{};

  for (int vectorIndex = 0; vectorIndex < 3; ++vectorIndex)
  {
    const std::string line =
        readRequiredLine(scanner, std::format("missing lattice vector {}", vectorIndex + 1));
    const std::vector<std::string> tokens = splitTokens(line);
    if (tokens.size() < 3)
    {
      throw std::runtime_error(std::format("POSCAR parse error: incomplete lattice vector {}", vectorIndex + 1));
    }

    double* components[3][3] = {{&unitCell.ax, &unitCell.ay, &unitCell.az},
                                {&unitCell.bx, &unitCell.by, &unitCell.bz},
                                {&unitCell.cx, &unitCell.cy, &unitCell.cz}};
    for (std::size_t componentIndex = 0; componentIndex < 3; ++componentIndex)
    {
      *components[vectorIndex][componentIndex] = scale * std::stod(tokens[componentIndex]);
    }
  }

  return unitCell;
}

std::vector<int> parseSpeciesCounts(const std::vector<std::string>& tokens, const std::string& context)
{
  if (tokens.empty())
  {
    throw std::runtime_error(std::format("POSCAR parse error: {}", context));
  }

  std::vector<int> counts{};
  counts.reserve(tokens.size());
  for (const std::string& token : tokens)
  {
    if (!isIntegerToken(token))
    {
      throw std::runtime_error(std::format("POSCAR parse error: expected species count in {}", context));
    }
    const int count = std::stoi(token);
    // a negative count would be summed into the number of atoms and cast to an unsigned size downstream
    if (count < 0)
    {
      throw std::runtime_error(std::format("POSCAR parse error: negative species count {} in {}", count, context));
    }
    counts.push_back(count);
  }

  return counts;
}

std::pair<std::vector<std::string>, std::vector<int>> readSpeciesLines(Scanner& scanner)
{
  const std::string firstSpeciesLine =
      readRequiredLine(scanner, "missing element symbols or species counts");
  const std::vector<std::string> firstTokens = splitTokens(firstSpeciesLine);

  if (isAllIntegerTokens(firstTokens))
  {
    return {{}, parseSpeciesCounts(firstTokens, "species counts line")};
  }

  const std::string countsLine = readRequiredLine(scanner, "missing species counts after element symbols");
  const std::vector<std::string> countTokens = splitTokens(countsLine);
  std::vector<int> counts = parseSpeciesCounts(countTokens, "species counts line");

  if (firstTokens.size() != counts.size())
  {
    throw std::runtime_error(std::format(
        "POSCAR parse error: number of element symbols ({}) does not match number of species counts ({})",
        firstTokens.size(), counts.size()));
  }

  return {firstTokens, counts};
}

bool isSelectiveDynamicsLine(const std::string& line)
{
  return startsWithIgnoreCase(line, "selective");
}

void readCoordinateModeLine(Scanner& scanner, bool& cartesian, bool& hasSelectiveDynamics)
{
  hasSelectiveDynamics = false;

  std::string line = readRequiredLine(scanner, "missing coordinate mode or selective dynamics flag");
  if (isSelectiveDynamicsLine(line))
  {
    hasSelectiveDynamics = true;
    line = readRequiredLine(scanner, "missing coordinate mode after selective dynamics flag");
  }

  const std::vector<std::string> tokens = splitTokens(line);
  if (tokens.empty())
  {
    throw std::runtime_error("POSCAR parse error: missing Direct/Cartesian coordinate mode");
  }

  const std::string& mode = tokens.front();
  if (equalsIgnoreCase(mode, "direct") || equalsIgnoreCase(mode, "d"))
  {
    cartesian = false;
    return;
  }

  if (equalsIgnoreCase(mode, "cartesian") || equalsIgnoreCase(mode, "c") || equalsIgnoreCase(mode, "k"))
  {
    cartesian = true;
    return;
  }

  throw std::runtime_error(
      std::format("POSCAR parse error: expected Direct or Cartesian coordinate mode, got '{}'", mode));
}

double3 parseAtomPosition(const std::vector<std::string>& tokens, bool hasSelectiveDynamics)
{
  const std::size_t minimumTokenCount = hasSelectiveDynamics ? 6 : 3;
  if (tokens.size() < minimumTokenCount)
  {
    throw std::runtime_error("POSCAR parse error: incomplete atom coordinate line");
  }

  return double3(std::stod(tokens[0]), std::stod(tokens[1]), std::stod(tokens[2]));
}

std::string readFileContents(const std::filesystem::path& path)
{
  std::ifstream input(path);
  if (!input)
  {
    throw std::runtime_error(std::format("POSCAR parse error: could not open '{}'", path.string()));
  }

  const std::string content((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());
  if (content.empty())
  {
    throw std::runtime_error(std::format("POSCAR parse error: file '{}' is empty", path.string()));
  }

  return content;
}

void assignAtomElement(std::shared_ptr<SKAsymmetricAtom>& atom, const std::string& chemicalElementString,
                     std::size_t speciesIndex)
{
  atom->setDisplayName(chemicalElementString);
  // a VASP species label carries its POTCAR variant ('Fe_pv', 'O_s', 'Fe_pv/5ae3c0' in VASP 5.4 and later): the
  // element is the part before it (the whole label missed the table and took the species index as atomic number,
  // hydrogen for the first species)
  const std::string element = chemicalElementString.substr(0, chemicalElementString.find_first_of("_/"));
  if (const auto elementIterator = PredefinedElements::atomicNumberData.find(element);
      elementIterator != PredefinedElements::atomicNumberData.end())
  {
    atom->setElementIdentifier(elementIterator->second);
  }
  else
  {
    atom->setElementIdentifier(speciesIndex + 1);
  }
}


SKPOSCARParser::SKPOSCARParser(const std::string& content, bool proteinOnlyAsymmetricUnitCell, bool asMolecule,
                               CharacterSet charactersToBeSkipped)
    : SKParser(),
      _scanner(content, charactersToBeSkipped),
      _proteinOnlyAsymmetricUnitCell(proteinOnlyAsymmetricUnitCell),
      _asMolecule(asMolecule),
      _frame(std::make_shared<SKStructure>()),
      _spaceGroupHallNumber(1)
{
  _frame->kind = SKStructure::Kind::crystal;
}

SKPOSCARParser SKPOSCARParser::fromFile(const std::filesystem::path& path, bool onlyAsymmetricUnitCell,
                                        bool asMolecule)
{
  return fromContent(readFileContents(path), onlyAsymmetricUnitCell, asMolecule);
}

SKPOSCARParser SKPOSCARParser::fromContent(const std::string& content, bool onlyAsymmetricUnitCell, bool asMolecule)
{
  return SKPOSCARParser(content, onlyAsymmetricUnitCell, asMolecule);
}

void SKPOSCARParser::startParsing() noexcept(false)
{
  // the first line is the comment, which may be empty (reading it as a token would skip a blank line and take the
  // scaling factor for the comment)
  std::string comment;
  _scanner.scanLine(comment);
  const double scaleFactor = readScaleFactor(_scanner);
  // a negative scaling factor is the volume of the cell: the lattice vectors are scaled to it
  const double3x3 rawUnitCell = readLatticeVectors(_scanner, 1.0);
  // a zero scaling factor or a cell without volume gave infinite or NaN coordinates instead of an error
  if (scaleFactor == 0.0)
  {
    throw std::runtime_error("POSCAR parse error: the universal scaling factor is zero");
  }
  const double rawVolume = std::abs(rawUnitCell.determinant());
  if (!(rawVolume > 0.0))
  {
    throw std::runtime_error("POSCAR parse error: the lattice vectors do not span a cell (zero volume)");
  }
  const double scale = scaleFactor < 0.0 ? std::cbrt(-scaleFactor / rawVolume) : scaleFactor;
  const double3x3 unitCell = scale * rawUnitCell;
  const auto [elementLabels, speciesCounts] = readSpeciesLines(_scanner);

  bool cartesian = false;
  bool hasSelectiveDynamics = false;
  readCoordinateModeLine(_scanner, cartesian, hasSelectiveDynamics);

  _frame->cell = std::make_shared<SKCell>(unitCell);
  const double3x3 inverseUnitCell = unitCell.inverse();

  for (std::size_t speciesIndex = 0; speciesIndex < speciesCounts.size(); ++speciesIndex)
  {
    const int numberOfAtoms = speciesCounts[speciesIndex];
    const std::string chemicalElementString =
        speciesIndex < elementLabels.size() ? trim(elementLabels[speciesIndex])
                                            : std::to_string(speciesIndex + 1);

    for (int atomIndex = 0; atomIndex < numberOfAtoms; ++atomIndex)
    {
      const std::string atomLine =
          readRequiredLine(_scanner, "unexpected end of file while reading atom positions");
      const std::vector<std::string> tokens = splitTokens(atomLine);
      const double3 position = parseAtomPosition(tokens, hasSelectiveDynamics);

      std::shared_ptr<SKAsymmetricAtom> atom = std::make_shared<SKAsymmetricAtom>();
      assignAtomElement(atom, chemicalElementString, speciesIndex);
      // Cartesian positions are scaled by the scaling factor like the lattice vectors
      atom->setPosition(cartesian ? inverseUnitCell * (scale * position) : position);
      _frame->atoms.push_back(atom);
    }
  }

  _movies.push_back({_frame});
}
