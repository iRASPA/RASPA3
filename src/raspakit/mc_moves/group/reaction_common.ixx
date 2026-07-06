module;

export module mc_moves_reaction_common;

import std;

import atom;
import randomnumbers;
import running_energy;
import cbmc_chain_data;
import reaction;
import system;
import mc_moves_move_types;

export namespace MC_Moves::ReactionCommon
{
enum class ReactionMoveKind : std::uint8_t
{
  LambdaChange = 0,
  ForwardInsert = 1,
  BackwardDelete = 2
};

struct MoleculeGroupGrowData
{
  std::vector<ChainGrowData> molecules;
  double RosenbluthWeight{1.0};
  RunningEnergy energies{};
};

struct MoleculeGroupRetraceData
{
  std::vector<ChainRetraceData> molecules;
  double RosenbluthWeight{1.0};
  RunningEnergy energies{};
};

void appendAllReactionFractionalMoleculeExclusions(const System& system,
                                                   std::vector<std::pair<std::size_t, std::size_t>>& exclude) noexcept;

[[nodiscard]] std::size_t totalStoichiometry(std::span<const std::size_t> stoichiometry) noexcept;

[[nodiscard]] bool selectRandomIntegerMolecules(RandomNumber& random, System& system,
                                                std::span<const std::size_t> stoichiometry,
                                                std::vector<std::pair<std::size_t, std::size_t>>& selected) noexcept;

[[nodiscard]] std::optional<MoleculeGroupGrowData> growMoleculeGroupInsertion(
    RandomNumber& random, System& system, std::span<const std::size_t> stoichiometry,
    std::span<const std::pair<std::size_t, std::size_t>> excludeMolecules, double scaling = 1.0,
    bool isFractional = false, std::uint8_t dUdlambdaGroupId = 0) noexcept;

[[nodiscard]] std::optional<MoleculeGroupRetraceData> retraceMoleculeGroupDeletion(
    RandomNumber& random, System& system,
    std::span<const std::pair<std::size_t, std::size_t>> selectedMolecules) noexcept;

[[nodiscard]] double idealGasRosenbluthWeightProduct(const System& system,
                                                     std::span<const std::size_t> stoichiometry) noexcept;

void setReactionFractionalScaling(System& system, Reaction& reaction, double lambda) noexcept;

[[nodiscard]] double computeReactionEquilibriumLogTerm(const System& system, const Reaction& reaction,
                                                       ReactionMoveKind moveKind) noexcept;

void deleteSelectedMolecules(System& system,
                             std::span<const std::pair<std::size_t, std::size_t>> selectedMolecules) noexcept;

void insertGrownMolecules(System& system, std::span<const ChainGrowData> growData,
                          std::span<const std::size_t> productStoichiometry) noexcept;

[[nodiscard]] std::optional<RunningEnergy> parallelReactionMove(RandomNumber& random, System& system,
                                                                Move::Types move) noexcept;

enum class SerialMoveKind : std::uint8_t
{
  LambdaChange = 0,
  FractionalReaction = 1,
  WholeMoleculeReaction = 2
};

void setSerialReactionFractionalScaling(System& system, Reaction& reaction, double lambda) noexcept;

// Monolithic serial Rx/CFC move: contains the lambda-change, fractional-reaction, and
// whole-molecule-reaction sub-moves. The sub-move is picked randomly using the lambda switch point,
// unless 'forcedKind' selects one deterministically (used by the tests).
[[nodiscard]] std::optional<RunningEnergy> serialReactionMove(
    RandomNumber& random, System& system, Reaction& reaction, Move::Types move, bool useCBMC,
    std::optional<SerialMoveKind> forcedKind = std::nullopt) noexcept;

}  // namespace MC_Moves::ReactionCommon
