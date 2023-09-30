export module move_statistics;

import <algorithm>;

export template<typename T>
struct MoveStatistics
{
  T counts{};
  T constructed{};
  T accepted{};
  T totalCounts{};
  T totalConstructed{};
  T totalAccepted{};
  T maxChange{};
  T targetAcceptance{ 0.5};

  void clear()
  {
    counts = T();
    constructed = T();
    accepted = T();
  }

  void optimizeAcceptance(T lowerLimit = T(0.0), T upperLimit = T(1.0))
  {
    T ratio = accepted / (counts + T(1.0));
    if constexpr (std::is_same_v<double, T>) 
    {
      T scaling = std::clamp( ratio / targetAcceptance, T(0.5), T(1.5) );
      maxChange = std::clamp( maxChange * scaling, lowerLimit, upperLimit );
    }
    else
    {
      T scaling = clamp(ratio / targetAcceptance, T(0.5), T(1.5));
      maxChange = clamp(maxChange * scaling, lowerLimit, upperLimit);
    }
    counts = T();
    constructed = T();
    accepted = T();
  }
};
