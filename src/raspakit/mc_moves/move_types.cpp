module;

module mc_moves_move_types;


std::array<std::string, std::to_underlying(Move::Types::Count)> Move::moveNames = 
{
  "Translation",
  "Random translation",
  "Rotation",
  "Random rotation",
  "Volume change",
  "Reinsertion (CBMC)",
  "Partial reinsertion (CBMC)",
  "Identity change (CBMC)",
  "Swap",
  "Swap (CBMC)",
  "Swap (CFCMC)",
  "Swap (CB/CFCMC)",
  "Gibbs volume",
  "Gibbs swap (CBMC)",
  "Gibbs swap (CFCMC)",
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
  "Gibbs swap (CB/CFCMC)",
  "Gibbs conventional CFCMC",
  "Gibbs conventional CFCMC/CBMC",
  "Pair swap",
  "Pair swap (CBMC)",
  "Anisotropic volume change",
  "Pair swap (CFCMC)",
  "Pair swap (CB/CFCMC)"
};

