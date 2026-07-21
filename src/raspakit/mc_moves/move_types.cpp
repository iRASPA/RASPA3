module;

module mc_moves_move_types;


std::array<std::string, std::to_underlying(Move::Types::Count)> Move::moveNames = 
{
  "Translation",
  "Random translation",
  "Rotation",
  "Random rotation",
  "Translation smart MC",
  "Translation smart MC (all)",
  "Rotation smart MC",
  "Rotation smart MC (all)",
  "Trans/rot smart MC",
  "Volume change",
  "Anisotropic volume change",
  "Reinsertion (CBMC)",
  "Partial reinsertion (CBMC)",
  "Identity change (CBMC)",
  "Swap",
  "Swap (CBMC)",
  "Swap (CFCMC)",
  "Swap (CB/CFCMC)",
  "Gibbs volume",
  "Gibbs swap (CBMC)",
  "Gibbs conventional CFCMC",
  "Gibbs conventional CFCMC/CBMC",
  "Gibbs swap (CFCMC)",
  "Gibbs swap (CB/CFCMC)",
  "Gibbs identity change (CBMC)",
  "Widom",
  "Widom (CFCMC)",
  "Widom (CB/CFCMC)",
  "Parallel tempering",
  "Hybrid MC",
  "Reaction (CBMC)",
  "Reaction (conventional CFCMC)",
  "Reaction (conventional CFCMC/CBMC)",
  "Reaction (CFCMC)",
  "Reaction (CFCMC/CBMC)",
  "Pair swap",
  "Pair swap (CBMC)",
  "Pair swap (CFCMC)",
  "Pair swap (CB/CFCMC)",
  "Group swap",
  "Group swap (CBMC)",
  "Group swap (CFCMC)",
  "Group swap (CB/CFCMC)"
};

std::array<std::string, std::to_underlying(Move::Timing::Count)> Move::timingNames =
{
  "Total",
  "ExternalField",
  "Framework",
  "Molecule",
  "NonEwald",
  "Ewald",
  "Tail",
  "ExternalField-Molecule",
  "Framework-Molecule",
  "Molecule-Molecule",
  "Insertion-Total",
  "Insertion-ExternalField",
  "Insertion-Framework",
  "Insertion-Molecule",
  "Insertion-NonEwald",
  "Insertion-Ewald",
  "Insertion-Tail",
  "Deletion-Total",
  "Deletion-ExternalField",
  "Deletion-Framework",
  "Deletion-Molecule",
  "Deletion-NonEwald",
  "Deletion-Ewald",
  "Deletion-Tail",
  "Lambda-ExternalField",
  "Lambda-Framework",
  "Lambda-Molecule",
  "Lambda-NonEwald",
  "Lambda-Ewald",
  "Lambda-Tail",
  "LambdaInterchange-NonEwald",
  "LambdaInterchange-Ewald",
  "LambdaInterchange-Tail",
  "LambdaChange-NonEwald",
  "LambdaChange-Ewald",
  "LambdaChange-Tail",
  "LambdaShuffle-NonEwald",
  "LambdaShuffle-Ewald",
  "LambdaShuffle-Tail",
  "Energy",
  "Fugacity",
  "Integration"
};

std::array<std::vector<Move::Timing>, std::to_underlying(Move::Types::Count)> Move::timingKeys =
{
  // Translation
  std::vector<Move::Timing>{Move::Timing::ExternalFieldMolecule, Move::Timing::FrameworkMolecule,
                            Move::Timing::MoleculeMolecule, Move::Timing::Ewald},
  // RandomTranslation
  std::vector<Move::Timing>{Move::Timing::ExternalFieldMolecule, Move::Timing::FrameworkMolecule,
                            Move::Timing::MoleculeMolecule, Move::Timing::Ewald},
  // Rotation
  std::vector<Move::Timing>{Move::Timing::ExternalFieldMolecule, Move::Timing::FrameworkMolecule,
                            Move::Timing::MoleculeMolecule, Move::Timing::Ewald},
  // RandomRotation
  std::vector<Move::Timing>{Move::Timing::ExternalFieldMolecule, Move::Timing::FrameworkMolecule,
                            Move::Timing::MoleculeMolecule, Move::Timing::Ewald},
  // TranslationSmartMC
  std::vector<Move::Timing>{Move::Timing::Integration},
  // TranslationSmartMCAll
  std::vector<Move::Timing>{Move::Timing::Integration},
  // RotationSmartMC
  std::vector<Move::Timing>{Move::Timing::Integration},
  // RotationSmartMCAll
  std::vector<Move::Timing>{Move::Timing::Integration},
  // TranslationRotationSmartMC
  std::vector<Move::Timing>{Move::Timing::Integration},
  // VolumeChange
  std::vector<Move::Timing>{Move::Timing::NonEwald, Move::Timing::Tail, Move::Timing::Ewald},
  // AnisotropicVolumeChange
  std::vector<Move::Timing>{Move::Timing::NonEwald, Move::Timing::Tail, Move::Timing::Ewald},
  // ReinsertionCBMC
  std::vector<Move::Timing>{Move::Timing::NonEwald, Move::Timing::Ewald},
  // PartialReinsertionCBMC
  std::vector<Move::Timing>{Move::Timing::NonEwald, Move::Timing::Ewald},
  // IdentityChangeCBMC
  std::vector<Move::Timing>{},
  // Swap
  std::vector<Move::Timing>{Move::Timing::InsertionTotal, Move::Timing::DeletionTotal, Move::Timing::NonEwald,
                            Move::Timing::Tail, Move::Timing::Ewald},
  // SwapCBMC
  std::vector<Move::Timing>{Move::Timing::InsertionTotal, Move::Timing::DeletionTotal, Move::Timing::NonEwald,
                            Move::Timing::Tail, Move::Timing::Ewald},
  // SwapCFCMC
  std::vector<Move::Timing>{Move::Timing::InsertionExternalField, Move::Timing::InsertionFramework,
                            Move::Timing::InsertionMolecule, Move::Timing::InsertionEwald, Move::Timing::InsertionTail,
                            Move::Timing::DeletionExternalField, Move::Timing::DeletionFramework,
                            Move::Timing::DeletionMolecule, Move::Timing::DeletionEwald, Move::Timing::DeletionTail,
                            Move::Timing::LambdaExternalField, Move::Timing::LambdaFramework,
                            Move::Timing::LambdaMolecule, Move::Timing::LambdaEwald, Move::Timing::LambdaTail},
  // SwapCBCFCMC
  std::vector<Move::Timing>{Move::Timing::InsertionExternalField, Move::Timing::InsertionFramework,
                            Move::Timing::InsertionMolecule, Move::Timing::InsertionNonEwald,
                            Move::Timing::InsertionEwald, Move::Timing::InsertionTail,
                            Move::Timing::DeletionExternalField, Move::Timing::DeletionFramework,
                            Move::Timing::DeletionMolecule, Move::Timing::DeletionNonEwald,
                            Move::Timing::DeletionEwald, Move::Timing::DeletionTail,
                            Move::Timing::LambdaExternalField, Move::Timing::LambdaFramework,
                            Move::Timing::LambdaMolecule, Move::Timing::LambdaNonEwald, Move::Timing::LambdaEwald,
                            Move::Timing::LambdaTail},
  // GibbsVolume
  std::vector<Move::Timing>{Move::Timing::NonEwald, Move::Timing::Tail, Move::Timing::Ewald},
  // GibbsSwapCBMC
  std::vector<Move::Timing>{Move::Timing::NonEwald, Move::Timing::Tail, Move::Timing::Ewald},
  // GibbsConventionalCFCMC
  std::vector<Move::Timing>{},
  // GibbsConventionalCBCFCMC
  std::vector<Move::Timing>{},
  // GibbsSwapCFCMC
  std::vector<Move::Timing>{Move::Timing::LambdaInterchangeNonEwald, Move::Timing::LambdaInterchangeEwald,
                            Move::Timing::LambdaInterchangeTail, Move::Timing::LambdaChangeNonEwald,
                            Move::Timing::LambdaChangeEwald, Move::Timing::LambdaChangeTail,
                            Move::Timing::LambdaShuffleNonEwald, Move::Timing::LambdaShuffleEwald,
                            Move::Timing::LambdaShuffleTail},
  // GibbsSwapCBCFCMC
  std::vector<Move::Timing>{Move::Timing::LambdaInterchangeNonEwald, Move::Timing::LambdaInterchangeEwald,
                            Move::Timing::LambdaInterchangeTail, Move::Timing::LambdaChangeNonEwald,
                            Move::Timing::LambdaChangeEwald, Move::Timing::LambdaChangeTail,
                            Move::Timing::LambdaShuffleNonEwald, Move::Timing::LambdaShuffleEwald,
                            Move::Timing::LambdaShuffleTail},
  // GibbsIdentityChangeCBMC
  std::vector<Move::Timing>{},
  // Widom
  std::vector<Move::Timing>{Move::Timing::NonEwald, Move::Timing::Tail, Move::Timing::Ewald},
  // WidomCFCMC
  std::vector<Move::Timing>{Move::Timing::ExternalField, Move::Timing::Molecule, Move::Timing::Framework,
                            Move::Timing::Ewald, Move::Timing::Tail},
  // WidomCBCFCMC
  std::vector<Move::Timing>{Move::Timing::ExternalField, Move::Timing::Molecule, Move::Timing::Framework,
                            Move::Timing::Ewald, Move::Timing::NonEwald, Move::Timing::Tail},
  // ParallelTempering
  std::vector<Move::Timing>{Move::Timing::Energy, Move::Timing::Fugacity},
  // HybridMC
  std::vector<Move::Timing>{Move::Timing::Integration},
  // ReactionCBMC
  std::vector<Move::Timing>{},
  // ReactionConventionalCFCMC
  std::vector<Move::Timing>{},
  // ReactionConventionalCBCFCMC
  std::vector<Move::Timing>{},
  // ReactionCFCMC
  std::vector<Move::Timing>{},
  // ReactionCBCFCMC
  std::vector<Move::Timing>{},
  // PairSwap
  std::vector<Move::Timing>{Move::Timing::InsertionTotal, Move::Timing::DeletionTotal, Move::Timing::NonEwald,
                            Move::Timing::Tail, Move::Timing::Ewald},
  // PairSwapCBMC
  std::vector<Move::Timing>{Move::Timing::InsertionTotal, Move::Timing::DeletionTotal, Move::Timing::NonEwald,
                            Move::Timing::Tail, Move::Timing::Ewald},
  // PairSwapCFCMC
  std::vector<Move::Timing>{Move::Timing::InsertionNonEwald, Move::Timing::InsertionEwald, Move::Timing::InsertionTail,
                            Move::Timing::DeletionNonEwald, Move::Timing::DeletionEwald, Move::Timing::DeletionTail,
                            Move::Timing::LambdaNonEwald, Move::Timing::LambdaEwald, Move::Timing::LambdaTail},
  // PairSwapCBCFCMC
  std::vector<Move::Timing>{Move::Timing::InsertionNonEwald, Move::Timing::InsertionEwald, Move::Timing::InsertionTail,
                            Move::Timing::DeletionNonEwald, Move::Timing::DeletionEwald, Move::Timing::DeletionTail,
                            Move::Timing::LambdaNonEwald, Move::Timing::LambdaEwald, Move::Timing::LambdaTail},
  // GroupSwap
  std::vector<Move::Timing>{Move::Timing::InsertionTotal, Move::Timing::DeletionTotal, Move::Timing::NonEwald,
                            Move::Timing::Tail, Move::Timing::Ewald},
  // GroupSwapCBMC
  std::vector<Move::Timing>{Move::Timing::InsertionTotal, Move::Timing::DeletionTotal, Move::Timing::NonEwald,
                            Move::Timing::Tail, Move::Timing::Ewald},
  // GroupSwapCFCMC
  std::vector<Move::Timing>{Move::Timing::InsertionNonEwald, Move::Timing::InsertionEwald, Move::Timing::InsertionTail,
                            Move::Timing::DeletionNonEwald, Move::Timing::DeletionEwald, Move::Timing::DeletionTail,
                            Move::Timing::LambdaNonEwald, Move::Timing::LambdaEwald, Move::Timing::LambdaTail},
  // GroupSwapCBCFCMC
  std::vector<Move::Timing>{Move::Timing::InsertionNonEwald, Move::Timing::InsertionEwald, Move::Timing::InsertionTail,
                            Move::Timing::DeletionNonEwald, Move::Timing::DeletionEwald, Move::Timing::DeletionTail,
                            Move::Timing::LambdaNonEwald, Move::Timing::LambdaEwald, Move::Timing::LambdaTail}
};

