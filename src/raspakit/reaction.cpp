module;

module reaction;

import std;

import archive;
import stringutils;
import mc_moves_move_types;

std::string Reaction::printStatus() const
{
  std::ostringstream stream;

  std::string r, p;
  for (const std::size_t &i : reactantStoichiometry) r += (r.empty() ? "" : ",") + std::to_string(i);
  for (const std::size_t &i : productStoichiometry) p += (p.empty() ? "" : ",") + std::to_string(i);
  std::print(stream, "    reaction [{}]: Stoichiometry reactants {} -> products {}\n", id, r, p);

  return stream.str();
}

nlohmann::json Reaction::jsonStatus() const
{
  nlohmann::json status;
  status["id"] = id;

  std::string r, p;
  for (const std::size_t &i : reactantStoichiometry) r += (r.empty() ? "" : ",") + std::to_string(i);
  for (const std::size_t &i : productStoichiometry) p += (p.empty() ? "" : ",") + std::to_string(i);
  status["reactants"] = r;
  status["products"] = p;

  return status;
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const Reaction &r)
{
  archive << r.versionNumber;

  archive << r.id;
  archive << r.reactantStoichiometry;
  archive << r.productStoichiometry;
  archive << r.lambda;
  archive << r.lambdaProductSide;
  archive << r.currentLambda;
  archive << r.maximumLambdaChange;
  archive << r.maximumLambdaChangeProducts;
  archive << r.lambdaSwitchPoint;
  archive << static_cast<std::uint64_t>(r.reactionMove);
  archive << r.fractionalSideIsReactants;
  archive << r.reactantFractionalMoleculeIds;
  archive << r.productFractionalMoleculeIds;

#if DEBUG_ARCHIVE
  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, Reaction &r)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > r.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'Reaction' at line {} in file {}\n", location.line(),
                                         location.file_name()));
  }

  archive >> r.id;
  archive >> r.reactantStoichiometry;
  archive >> r.productStoichiometry;
  archive >> r.lambda;
  if (versionNumber >= 3)
  {
    archive >> r.lambdaProductSide;
    archive >> r.currentLambda;
    archive >> r.maximumLambdaChange;
    archive >> r.maximumLambdaChangeProducts;
    archive >> r.lambdaSwitchPoint;
    if (versionNumber >= 4)
    {
      std::uint64_t reactionMove;
      archive >> reactionMove;
      r.reactionMove = static_cast<Move::Types>(reactionMove);
    }
    else
    {
      // version 3 stored a per-reaction serial/parallel flag; map it onto the driving move
      bool legacySerialRxCFC;
      archive >> legacySerialRxCFC;
      r.reactionMove =
          legacySerialRxCFC ? Move::Types::ReactionCFCMC : Move::Types::ReactionConventionalCFCMC;
    }
    archive >> r.fractionalSideIsReactants;
    archive >> r.reactantFractionalMoleculeIds;
    archive >> r.productFractionalMoleculeIds;
  }
  else if (versionNumber >= 2)
  {
    archive >> r.currentLambda;
    archive >> r.maximumLambdaChange;
    archive >> r.reactantFractionalMoleculeIds;
    archive >> r.productFractionalMoleculeIds;
  }

#if DEBUG_ARCHIVE
  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("Reaction: Error in binary restart\n"));
  }
#endif

  return archive;
}
