export module move_statistics;

import <algorithm>;

export template<typename T>
struct MoveStatistics
{
    T counts{};
    T constructed{};
    T accepted{};
    T maxChange{};
    T targetAcceptance{ 0.5};

    void clear()
    {
        counts = T();
        constructed = T();
        accepted = T();
    }

    void optimizeAcceptance()
    {
      T ratio = accepted / (counts + T(1.0));
      if constexpr (std::is_same_v<double, T>) 
      {
        T scaling = std::clamp( ratio / targetAcceptance, T(0.5), T(1.5) );
        maxChange = std::clamp( maxChange * scaling, T(0.01), T(1.5) );
      }
      else
      {
        T scaling = clamp(ratio / targetAcceptance, T(0.5), T(1.5));
        maxChange = clamp(maxChange * scaling, T(0.01), T(1.5));
      }
      counts = T();
      constructed = T();
      accepted = T();
    }

};
