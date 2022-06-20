export module probability_histogram;

import <vector>;
import <iostream>;
import <cmath>;
import <string>;
import <sstream>;
import <format>;
import <algorithm>;
import <numeric>;
import <numbers>;

export struct ProbabilityHistogram
{
  ProbabilityHistogram(size_t numberOfBins) :
      numberOfBins(numberOfBins),
      histogram(numberOfBins)
  {
  }

  size_t numberOfBins;

  std::vector<double> histogram;

  inline ProbabilityHistogram& operator+=(const ProbabilityHistogram& b)
  {
    for (size_t i = 0; i != this->histogram.size(); ++i)
    {
      histogram[i] += b.histogram[i];
    }
    return *this;
  }

  inline ProbabilityHistogram& operator-=(const ProbabilityHistogram& b)
  {
    for (size_t i = 0; i != this->histogram.size(); ++i)
    {
      histogram[i] -= b.histogram[i];
    }
    return *this;
  }

  inline ProbabilityHistogram operator-() const
  {
    ProbabilityHistogram v(histogram.size());
    for (size_t i = 0; i < this->histogram.size(); ++i)
    {
      v.histogram[i] = -histogram[i];
    }

    return v;
  }

  inline ProbabilityHistogram compositeProperty [[nodiscard]] () const
  {
    ProbabilityHistogram v(histogram.size());

    for (size_t i = 0; i < this->histogram.size(); ++i)
    {
      v.histogram[i] = -std::log(histogram[i]);
    }
    return v;
  }

};

export inline ProbabilityHistogram operator+(const ProbabilityHistogram& a, const ProbabilityHistogram& b)
{
  ProbabilityHistogram m(a.numberOfBins);
  for (size_t i = 0; i < a.histogram.size(); ++i)
  {
    m.histogram[i] = a.histogram[i] + b.histogram[i];
  }

  return m;
}

export inline ProbabilityHistogram operator-(const ProbabilityHistogram& a, const ProbabilityHistogram& b)
{
  ProbabilityHistogram m(a.numberOfBins);
  for (size_t i = 0; i < a.histogram.size(); ++i)
  {
    m.histogram[i] = a.histogram[i] - b.histogram[i];
  }

  return m;
}

export inline ProbabilityHistogram operator*(const ProbabilityHistogram& a, const ProbabilityHistogram& b)
{
  ProbabilityHistogram m(a.numberOfBins);
  for (size_t i = 0; i < a.histogram.size(); ++i)
  {
    m.histogram[i] = a.histogram[i] * b.histogram[i];
  }

  return m;
}

export inline ProbabilityHistogram operator*(const double& a, const ProbabilityHistogram& b)
{
  ProbabilityHistogram m(b.numberOfBins);
  for (size_t i = 0; i < b.histogram.size(); ++i)
  {
    m.histogram[i] = a * b.histogram[i];
  }

  return m;
}

export inline ProbabilityHistogram operator/(const ProbabilityHistogram& a, const double& b)
{
  ProbabilityHistogram m(a.numberOfBins);
  for (size_t i = 0; i < a.histogram.size(); ++i)
  {
    m.histogram[i] = a.histogram[i] / b;
  }

  return m;
}

export inline ProbabilityHistogram sqrt(const ProbabilityHistogram& a)
{
  ProbabilityHistogram m(a.numberOfBins);
  for (size_t i = 0; i < a.histogram.size(); ++i)
  {
    m.histogram[i] = sqrt(a.histogram[i]);
  }

  return m;
}
