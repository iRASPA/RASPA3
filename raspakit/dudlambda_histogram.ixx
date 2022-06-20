export module dudlambda_histogram;

import <vector>;
import <iostream>;
import <cmath>;
import <string>;
import <sstream>;
import <format>;
import <algorithm>;
import <numeric>;
import <numbers>;

export struct DUdlambdaHistogram
{
  DUdlambdaHistogram(size_t numberOfBins) :
      numberOfBins(numberOfBins),
      histogram(numberOfBins)
  {
  }

  size_t numberOfBins;

  std::vector<double> histogram;

  inline DUdlambdaHistogram& operator+=(const DUdlambdaHistogram& b)
  {
    for (size_t i = 0; i != this->histogram.size(); ++i)
    {
      histogram[i] += b.histogram[i];
    }
    return *this;
  }

  inline DUdlambdaHistogram& operator-=(const DUdlambdaHistogram& b)
  {
    for (size_t i = 0; i != this->histogram.size(); ++i)
    {
      histogram[i] -= b.histogram[i];
    }
    return *this;
  }

  inline DUdlambdaHistogram operator-() const
  {
    DUdlambdaHistogram v(histogram.size());
    for (size_t i = 0; i < this->histogram.size(); ++i)
    {
      v.histogram[i] = -histogram[i];
    }

    return v;
  }

  inline DUdlambdaHistogram compositeProperty [[nodiscard]] () const
  {
    DUdlambdaHistogram v(histogram.size());

    double delta = 1.0 / static_cast<double>(numberOfBins - 1);
    for (size_t i = 1; i < this->histogram.size(); ++i)
    {
      v.histogram[i] += delta * histogram[i -1];
    }
    return v;
  }

};

export inline DUdlambdaHistogram operator+(const DUdlambdaHistogram& a, const DUdlambdaHistogram& b)
{
  DUdlambdaHistogram m(a.numberOfBins);
  for (size_t i = 0; i < a.histogram.size(); ++i)
  {
    m.histogram[i] = a.histogram[i] + b.histogram[i];
  }

  return m;
}

export inline DUdlambdaHistogram operator-(const DUdlambdaHistogram& a, const DUdlambdaHistogram& b)
{
  DUdlambdaHistogram m(a.numberOfBins);
  for (size_t i = 0; i < a.histogram.size(); ++i)
  {
    m.histogram[i] = a.histogram[i] - b.histogram[i];
  }

  return m;
}

export inline DUdlambdaHistogram operator*(const DUdlambdaHistogram& a, const DUdlambdaHistogram& b)
{
  DUdlambdaHistogram m(a.numberOfBins);
  for (size_t i = 0; i < a.histogram.size(); ++i)
  {
    m.histogram[i] = a.histogram[i] * b.histogram[i];
  }

  return m;
}

export inline DUdlambdaHistogram operator*(const double& a, const DUdlambdaHistogram& b)
{
  DUdlambdaHistogram m(b.numberOfBins);
  for (size_t i = 0; i < b.histogram.size(); ++i)
  {
    m.histogram[i] = a * b.histogram[i];
  }

  return m;
}

export inline DUdlambdaHistogram operator/(const DUdlambdaHistogram& a, const double& b)
{
  DUdlambdaHistogram m(a.numberOfBins);
  for (size_t i = 0; i < a.histogram.size(); ++i)
  {
    m.histogram[i] = a.histogram[i] / b;
  }

  return m;
}

export inline DUdlambdaHistogram sqrt(const DUdlambdaHistogram& a)
{
  DUdlambdaHistogram m(a.numberOfBins);
  for (size_t i = 0; i < a.histogram.size(); ++i)
  {
    m.histogram[i] = sqrt(a.histogram[i]);
  }

  return m;
}
