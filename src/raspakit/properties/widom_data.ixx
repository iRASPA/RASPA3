module;

export module widom_data;

import std;

import archive;

export struct WidomData
{
  WidomData():
    total(0.0),
    excess(0.0),
    idealGas(0.0)
  {
  };

  WidomData(double total, double excess, double idealGas):
    total(total),
    excess(excess),
    idealGas(idealGas)
  {
  };

  inline WidomData& operator+=(const WidomData& b)
  {
    total += b.total;
    excess += b.excess;
    idealGas += b.idealGas;

    return *this;
  }

  std::uint64_t versionNumber{1};

  double total{};
  double excess{};
  double idealGas{};

  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const WidomData& l);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, WidomData& l);
};

export inline WidomData operator+(const WidomData& a, const WidomData& b)
{
  WidomData m{}; 

  m.total = a.total + b.total;
  m.excess = a.excess + b.excess;
  m.idealGas = a.idealGas + b.idealGas;

  return m;
}

export inline WidomData operator-(const WidomData& a, const WidomData& b)
{
  WidomData m{}; 

  m.total = a.total - b.total;
  m.excess = a.excess - b.excess;
  m.idealGas = a.idealGas - b.idealGas;

  return m;
}

export inline WidomData operator*(const WidomData& a, const WidomData& b)
{
  WidomData m{}; 

  m.total = a.total * b.total;
  m.excess = a.excess * b.excess;
  m.idealGas = a.idealGas * b.idealGas;

  return m;
}

export inline WidomData operator*(const double& a, const WidomData& b)
{
  WidomData m{}; 

  m.total = a * b.total;
  m.excess = a * b.excess;
  m.idealGas = a * b.idealGas;

  return m;
}

export inline WidomData operator/(const WidomData& a, const double& b)
{
  WidomData m{}; 

  double inv_b = 1.0 / b;
  m.total = inv_b * a.total;
  m.excess = inv_b * a.excess;
  m.idealGas = inv_b * a.idealGas;

  return m;
}

export inline WidomData operator/(const double& a, const WidomData& b)
{
  WidomData m{}; 

  m.total = a / b.total;
  m.excess = a / b.excess;
  m.idealGas = a / b.idealGas;

  return m;
}

export inline WidomData sqrt(const WidomData& a)
{
  WidomData m{};

  m.total = std::sqrt(a.total);
  m.excess = std::sqrt(a.excess);
  m.idealGas = std::sqrt(a.idealGas);

  return m;
}

export inline WidomData log(const WidomData& a)
{
  WidomData m{};

  m.total = std::log(a.total);
  m.excess = std::log(a.excess);
  m.idealGas = std::log(a.idealGas);

  return m;
}
