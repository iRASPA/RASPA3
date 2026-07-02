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

void applyLinearReactionScaling(std::span<Atom> atoms, bool isReactant, double lambda) noexcept;

void appendAllReactionFractionalMoleculeExclusions(
    const System& system, std::vector<std::pair<std::size_t, std::size_t>>& exclude) noexcept;

[[nodiscard]] std::size_t totalStoichiometry(std::span<const std::size_t> stoichiometry) noexcept;

[[nodiscard]] bool selectRandomIntegerMolecules(RandomNumber& random, System& system,
                                                std::span<const std::size_t> stoichiometry,
                                                std::vector<std::pair<std::size_t, std::size_t>>& selected) noexcept;

[[nodiscard]] std::optional<MoleculeGroupGrowData> growMoleculeGroupInsertion(
    RandomNumber& random, System& system, std::span<const std::size_t> stoichiometry,
    std::span<const std::pair<std::size_t, std::size_t>> excludeMolecules, double scaling = 1.0,
    bool isFractional = false) noexcept;

[[nodiscard]] std::optional<MoleculeGroupRetraceData> retraceMoleculeGroupDeletion(
    RandomNumber& random, System& system,
    std::span<const std::pair<std::size_t, std::size_t>> selectedMolecules) noexcept;

[[nodiscard]] std::vector<Atom> collectReactionFractionalAtoms(System& system, Reaction& reaction) noexcept;

[[nodiscard]] std::vector<Atom> collectReactionFractionalAtomsAtLambda(System& system, Reaction& reaction,
                                                                       double lambda) noexcept;

void setReactionFractionalScaling(System& system, Reaction& reaction, double lambda) noexcept;

[[nodiscard]] std::optional<RunningEnergy> computeReactionFractionalScalingEnergyDifference(
    System& system, Reaction& reaction, double lambdaOld, double lambdaNew,
    bool includeEwaldCorrections = true) noexcept;

[[nodiscard]] double computeReactionEquilibriumLogTerm(const System& system, const Reaction& reaction,
                                                       ReactionMoveKind moveKind) noexcept;

[[nodiscard]] std::optional<RunningEnergy> computeGroupSwapEnergyDifference(
    System& system, std::span<const Atom> newAtoms, std::span<const Atom> oldAtoms,
    bool includeTailCorrections = true, bool includeEwaldCorrections = true) noexcept;

[[nodiscard]] RunningEnergy computeGroupSwapTailEnergyDifference(System& system, std::span<const Atom> newAtoms,
                                                                 std::span<const Atom> oldAtoms) noexcept;

void deleteSelectedMolecules(System& system,
                             std::span<const std::pair<std::size_t, std::size_t>> selectedMolecules) noexcept;

void insertGrownMolecules(System& system, std::span<const ChainGrowData> growData,
                          std::span<const std::size_t> productStoichiometry) noexcept;

void acceptReactionForwardInsert(System& system, Reaction& reaction, double lambdaNew,
                                 std::span<const std::pair<std::size_t, std::size_t>> selectedMolecules,
                                 std::span<const ChainGrowData> growData) noexcept;

void acceptReactionBackwardDelete(System& system, Reaction& reaction, double lambdaNew,
                                  std::span<const std::pair<std::size_t, std::size_t>> selectedMolecules,
                                  std::span<const ChainGrowData> growData) noexcept;

[[nodiscard]] std::size_t numberOfReactionMoleculesForComponent(const Reaction& reaction, std::size_t componentId,
                                                                ReactionMoveKind moveKind) noexcept;

enum class SerialMoveKind : std::uint8_t
{
  LambdaChange = 0,
  FractionalReaction = 1,
  WholeMoleculeReaction = 2
};

void applySerialFractionalScaling(std::span<Atom> atoms, double lambda) noexcept;

void setSerialReactionFractionalScaling(System& system, Reaction& reaction, double lambda) noexcept;

[[nodiscard]] std::optional<RunningEnergy> serialLambdaChangeMove(RandomNumber& random, System& system,
                                                                   Reaction& reaction, Move::Types move) noexcept;

[[nodiscard]] std::optional<RunningEnergy> serialFractionalReactionMove(RandomNumber& random, System& system,
                                                                        Reaction& reaction, Move::Types move,
                                                                        bool useCBMC) noexcept;

[[nodiscard]] std::optional<RunningEnergy> serialWholeMoleculeReactionMove(RandomNumber& random, System& system,
                                                                           Reaction& reaction, Move::Types move,
                                                                           bool useCBMC) noexcept;

}  // namespace MC_Moves::ReactionCommon
