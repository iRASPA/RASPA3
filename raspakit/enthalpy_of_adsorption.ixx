export module enthalpy_of_adsorption;

import energy_status;
import averages;
import matrix;
import units;

import <vector>;
import <numeric>;
import <fstream>;
import <utility>;
import <string>;
import <cmath>;
import <iostream>;
#if defined(_WIN32)
  import <cassert>;
#else
  #include <assert.h>
#endif

export struct EnthalpyOfAdsorption
{
	EnthalpyOfAdsorption(size_t size) :
		size(size),
		values(size)
	{

	}

	EnthalpyOfAdsorption(std::vector<double> values) :
		size(values.size()),
		values(values)
	{

	}

	size_t size;
	std::vector<double> values;

	inline EnthalpyOfAdsorption& operator+=(const EnthalpyOfAdsorption& b)
	{
		for (size_t i = 0; i < size; ++i)
		{
			values[i] += b.values[i];
		}
		return *this;
	}
};


export inline EnthalpyOfAdsorption operator+(const EnthalpyOfAdsorption& a, const EnthalpyOfAdsorption& b)
{
	EnthalpyOfAdsorption m(a.size);

	for (size_t i = 0; i < m.size; ++i)
	{
		m.values[i] = a.values[i] + b.values[i];
	}

	return m;
}

export inline EnthalpyOfAdsorption operator-(const EnthalpyOfAdsorption& a, const EnthalpyOfAdsorption& b)
{
	EnthalpyOfAdsorption m(a.size);

	for (size_t i = 0; i < m.size; ++i)
	{
		m.values[i] = a.values[i] - b.values[i];
	}

	return m;
}

export inline EnthalpyOfAdsorption operator*(const EnthalpyOfAdsorption& a, const EnthalpyOfAdsorption& b)
{
	EnthalpyOfAdsorption m(a.size);

	for (size_t i = 0; i < m.size; ++i)
	{
		m.values[i] = a.values[i] * b.values[i];
	}

	return m;
}

export inline EnthalpyOfAdsorption operator*(const double& a, const EnthalpyOfAdsorption& b)
{
	EnthalpyOfAdsorption m(b.size);

	for (size_t i = 0; i < m.size; ++i)
	{
		m.values[i] = a * b.values[i];
	}

	return m;
}

export inline EnthalpyOfAdsorption sqrt(const EnthalpyOfAdsorption& a)
{
	EnthalpyOfAdsorption m(a.size);

	for (size_t i = 0; i < m.size; ++i)
	{
		m.values[i] = std::sqrt(a.values[i]);
	}

	return m;
}


export struct EnthalpyOfAdsorptionTerms
{
	EnthalpyOfAdsorptionTerms(size_t size) :
		size(size),
		swapableComponents(size),
		totalEnergyTimesNumberOfMolecules(size),
		numberOfMoleculesSquared(size, std::vector<double>(size)),
		numberOfMolecules(size),
                temperature(0.0),
		totalEnergy(0.0)
	{
	
	}

	EnthalpyOfAdsorptionTerms(const EnthalpyOfAdsorptionTerms& a) noexcept = default;
	EnthalpyOfAdsorptionTerms& operator=(const EnthalpyOfAdsorptionTerms& a) noexcept = default;


	EnthalpyOfAdsorptionTerms(const std::vector<size_t> swapableComponents,const std::vector<size_t> numberOfIntegerMolecules, double totalEnergy, double temperature) :
		size(swapableComponents.size()),
		swapableComponents(swapableComponents),
		totalEnergyTimesNumberOfMolecules(swapableComponents.size()),
		numberOfMoleculesSquared(swapableComponents.size(), std::vector<double>(swapableComponents.size())),
		numberOfMolecules(swapableComponents.size()),
		temperature(temperature * Units::KelvinToEnergy),
		totalEnergy(totalEnergy)
	{
		for (size_t i = 0; i < swapableComponents.size(); ++i)
		{
			size_t index_i = swapableComponents[i];
			numberOfMolecules[i] = static_cast<double>(numberOfIntegerMolecules[index_i]);
			
			totalEnergyTimesNumberOfMolecules[i] = totalEnergy * static_cast<double>(numberOfIntegerMolecules[index_i]);
			for (size_t j = 0; j < swapableComponents.size(); ++j)
			{
				size_t index_j = swapableComponents[j];
				numberOfMoleculesSquared[i][j] = static_cast<double>(numberOfIntegerMolecules[index_i]) * static_cast<double>(numberOfIntegerMolecules[index_j]);
			}
		}
	}

	EnthalpyOfAdsorptionTerms() = default;

	size_t size;
	std::vector<size_t> swapableComponents;
	std::vector<double> totalEnergyTimesNumberOfMolecules;
	std::vector<std::vector<double>> numberOfMoleculesSquared;
	std::vector<double> numberOfMolecules;
	double temperature;
	double totalEnergy;

	
	inline EnthalpyOfAdsorptionTerms& operator+=(const EnthalpyOfAdsorptionTerms& b)
	{
		totalEnergy += b.totalEnergy;
		temperature += b.temperature;
		for (size_t i = 0; i < size; ++i)
		{
			numberOfMolecules[i] += b.numberOfMolecules[i];
			totalEnergyTimesNumberOfMolecules[i] += b.totalEnergyTimesNumberOfMolecules[i];
			for (size_t j = 0; j < size; ++j)
			{
				numberOfMoleculesSquared[i][j] += b.numberOfMoleculesSquared[i][j];
			}
		}
		return *this;
	}

	inline EnthalpyOfAdsorption compositeProperty [[nodiscard]] () const
	{
		EnthalpyOfAdsorption v(size);
		if (size > 0)
		{
			// Symmetric matrix
			Matrix m(size, size, 0.0);
			for (size_t i = 0; i < size; ++i)
			{
				for (size_t j = 0; j < size; ++j)
				{
					m(i, j) = numberOfMoleculesSquared[i][j] - numberOfMolecules[i] * numberOfMolecules[j];
				}
			}

			m.inverse();

			for (size_t i = 0; i < size; ++i)
			{
				v.values[i] = 0.0;
				for (size_t j = 0; j < size; ++j)
				{
					v.values[i] += (totalEnergyTimesNumberOfMolecules[j] - totalEnergy * numberOfMolecules[j]) * m(j, i) - temperature;
				}
			}
		}
		return v;
	}
};

export inline EnthalpyOfAdsorptionTerms operator+(const EnthalpyOfAdsorptionTerms& a, const EnthalpyOfAdsorptionTerms& b)
{
	EnthalpyOfAdsorptionTerms m(a.size);

	m.totalEnergy = a.totalEnergy + b.totalEnergy;
	m.temperature = a.temperature + b.temperature;
	for (size_t i = 0; i < m.size; ++i)
	{
		m.numberOfMolecules[i] = a.numberOfMolecules[i] + b.numberOfMolecules[i];
		m.totalEnergyTimesNumberOfMolecules[i] = a.totalEnergyTimesNumberOfMolecules[i] + b.totalEnergyTimesNumberOfMolecules[i];
		for (size_t j = 0; j < m.size; ++j)
		{
			m.numberOfMoleculesSquared[i][j] = a.numberOfMoleculesSquared[i][j] + b.numberOfMoleculesSquared[i][j];
		}
	}
	
	return m;
}

export inline EnthalpyOfAdsorptionTerms operator*(const double& a, const EnthalpyOfAdsorptionTerms& b)
{
	EnthalpyOfAdsorptionTerms m(b.size);

	m.totalEnergy = a * b.totalEnergy;
	m.temperature = a * b.temperature;
	for (size_t i = 0; i < m.size; ++i)
	{
		m.numberOfMolecules[i] = a * b.numberOfMolecules[i];
		m.totalEnergyTimesNumberOfMolecules[i] = a * b.totalEnergyTimesNumberOfMolecules[i];
		for (size_t j = 0; j < m.size; ++j)
		{
			m.numberOfMoleculesSquared[i][j] = a * b.numberOfMoleculesSquared[i][j];
		}
	}

	return m;
}

export inline EnthalpyOfAdsorptionTerms operator/(const EnthalpyOfAdsorptionTerms& a, const double& b)
{
	EnthalpyOfAdsorptionTerms m(a.size);

	m.totalEnergy = a.totalEnergy / b;
	m.temperature = a.temperature / b;
	for (size_t i = 0; i < m.size; ++i)
	{
		m.numberOfMolecules[i] = a.numberOfMolecules[i] / b;
		m.totalEnergyTimesNumberOfMolecules[i] = a.totalEnergyTimesNumberOfMolecules[i] / b;
		for (size_t j = 0; j < m.size; ++j)
		{
			m.numberOfMoleculesSquared[i][j] = a.numberOfMoleculesSquared[i][j] / b;
		}
	}

	return m;
}
