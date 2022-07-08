export module simulationbox;

import double3;
import int3;
import double3x3;
import units;

import <numbers>;
import <string>;
import <iostream>;
import <cmath>;
import <algorithm>;

#if defined(__GNUC__)
#define ALWAYS_INLINE __attribute__((__always_inline__)) 
#elif defined(_MSC_VER)
#define ALWAYS_INLINE __forceinline
#endif


export struct SimulationBox
{
    enum class Type : int
    {
        Rectangular = 0,
        Triclinic = 1
    };

    SimulationBox() :
        lengthA(0.0), lengthB(0.0), lengthC(0.0), 
        angleAlpha(0.0), angleBeta(0.0), angleGamma(0.0),
        unitCell(double3x3(double3(0.0, 0.0, 0.0), double3(0.0, 0.0, 0.0), double3(0.0, 0.0, 0.0))),
        inverseUnitCell(double3x3(double3(0.0, 0.0, 0.0), double3(0.0, 0.0, 0.0), double3(0.0, 0.0, 0.0))),
        volume(0.0)
    {
    };

    SimulationBox(double a, double b, double c, double alpha, double beta, double gamma);
        

    ALWAYS_INLINE inline double3 applyPeriodicBoundaryConditions(const double3& dr) const
    {
        switch (type)
        {
        case SimulationBox::Type::Rectangular:
        {
            double3 s;
            s.x = dr.x - static_cast<int>(dr.x * inverseUnitCell.ax + ((dr.x >= 0.0) ? 0.5 : -0.5)) * unitCell.ax;
            s.y = dr.y - static_cast<int>(dr.y * inverseUnitCell.by + ((dr.y >= 0.0) ? 0.5 : -0.5)) * unitCell.by;
            s.z = dr.z - static_cast<int>(dr.z * inverseUnitCell.cz + ((dr.z >= 0.0) ? 0.5 : -0.5)) * unitCell.cz;
            return s;
        }
        default:
        {
            double3 s = inverseUnitCell * dr;
            s.x -= static_cast<int>(s.x + ((s.x >= 0.0) ? 0.5 : -0.5));
            s.y -= static_cast<int>(s.y + ((s.y >= 0.0) ? 0.5 : -0.5));
            s.z -= static_cast<int>(s.z + ((s.z >= 0.0) ? 0.5 : -0.5));
            return unitCell * s;
        }
        }
    }

    void setBoxLengths(double3 lengths);
    void setBoxAngles(double3 angles);
    double3 lengths();
    double3 angles();
    double3 randomPosition() const;

    std::string printParameters();
    std::string printStatus();
    std::string printStatus(const SimulationBox& average, const SimulationBox& error);

    double lengthA;
    double lengthB;
    double lengthC;
    double angleAlpha;
    double angleBeta;
    double angleGamma;
    double3x3 unitCell;
    double3x3 inverseUnitCell;
    double volume;
    double temperature = 300.0;
    double pressure = 1e4;
    double Beta = 1.0 / (Units::KB * 300.0);
    double alpha{ 0.265058 };
    int3 kmax{ 8, 8, 8 };
    Type type = Type::Rectangular;

    inline SimulationBox& operator+=(const SimulationBox& b)
    {
        lengthA += b.lengthA;
        lengthB += b.lengthB;
        lengthC += b.lengthC;
        angleAlpha += b.angleAlpha;
        angleBeta += b.angleBeta;
        angleGamma += b.angleGamma;
        unitCell += b.unitCell;
        volume += volume;
        temperature += b.temperature;
        pressure += b.pressure;

        return *this;
    }
};

export inline SimulationBox operator+(const SimulationBox& a, const SimulationBox& b)
{
    SimulationBox m;
    
    m.lengthA = a.lengthA + b.lengthA;
    m.lengthB = a.lengthB + b.lengthB;
    m.lengthC = a.lengthC + b.lengthC;
    m.angleAlpha = a.angleAlpha + b.angleAlpha;
    m.angleBeta = a.angleBeta + b.angleBeta;
    m.angleGamma = a.angleGamma + b.angleGamma;
    m.unitCell = a.unitCell + b.unitCell;
    m.volume = a.volume + b.volume;
    m.temperature = a.temperature + b.temperature;
    m.pressure = a.pressure + b.pressure;

    return m;
}

export inline SimulationBox operator-(const SimulationBox& a, const SimulationBox& b)
{
    SimulationBox m;

    m.lengthA = a.lengthA - b.lengthA;
    m.lengthB = a.lengthB - b.lengthB;
    m.lengthC = a.lengthC - b.lengthC;
    m.angleAlpha = a.angleAlpha - b.angleAlpha;
    m.angleBeta = a.angleBeta - b.angleBeta;
    m.angleGamma = a.angleGamma - b.angleGamma;
    m.unitCell = a.unitCell - b.unitCell;
    m.volume = a.volume - b.volume;
    m.temperature = a.temperature - b.temperature;
    m.pressure = a.pressure - b.pressure;

    return m;
}

export inline SimulationBox operator*(const SimulationBox& a, const SimulationBox& b)
{
    SimulationBox m;

    m.lengthA = a.lengthA * b.lengthA;
    m.lengthB = a.lengthB * b.lengthB;
    m.lengthC = a.lengthC * b.lengthC;
    m.angleAlpha = a.angleAlpha * b.angleAlpha;
    m.angleBeta = a.angleBeta * b.angleBeta;
    m.angleGamma = a.angleGamma * b.angleGamma;
    m.unitCell = a.unitCell * b.unitCell;
    m.volume = a.volume * b.volume;
    m.temperature = a.temperature * b.temperature;
    m.pressure = a.pressure * b.pressure;
    return m;
}

export inline SimulationBox operator*(const double& a, const SimulationBox& b)
{
    SimulationBox m;

    m.lengthA = a * b.lengthA;
    m.lengthB = a * b.lengthB;
    m.lengthC = a * b.lengthC;
    m.angleAlpha = a * b.angleAlpha;
    m.angleBeta = a * b.angleBeta;
    m.angleGamma = a * b.angleGamma;
    m.unitCell = a * b.unitCell;
    m.volume = a * b.volume;
    m.temperature = a * b.temperature;
    m.pressure = a * b.pressure;
    return m;
}


export inline SimulationBox operator/(const SimulationBox& a, const double& b)
{
    SimulationBox m;

    double temp = 1.0 / b;
    m.lengthA = a.lengthA * temp;
    m.lengthB = a.lengthB * temp;
    m.lengthC = a.lengthC * temp;
    m.angleAlpha = a.angleAlpha * temp;
    m.angleBeta = a.angleBeta * temp;
    m.angleGamma = a.angleGamma * temp;
    m.unitCell = a.unitCell * temp;
    m.volume = a.volume * temp;
    m.temperature = a.temperature * temp;
    m.pressure = a.pressure * temp;
    return m;
}


export inline SimulationBox sqrt(const SimulationBox& a)
{
    SimulationBox m;

    m.lengthA = sqrt(a.lengthA);
    m.lengthB = sqrt(a.lengthB);
    m.lengthC = sqrt(a.lengthC);
    m.angleAlpha = sqrt(a.angleAlpha);
    m.angleBeta = sqrt(a.angleBeta);
    m.angleGamma = sqrt(a.angleGamma);
    m.unitCell = sqrt(a.unitCell);
    m.volume = sqrt(a.volume);
    m.temperature = sqrt(a.temperature);
    m.pressure = sqrt(a.pressure);

    return m;
}

export inline SimulationBox max(const SimulationBox& a, const SimulationBox& b)
{
    return SimulationBox(std::max(a.lengthA, b.lengthA),
                         std::max(a.lengthB, b.lengthB),
                         std::max(a.lengthC, b.lengthC),
                         std::max(a.angleAlpha, b.angleAlpha),
                         std::max(a.angleBeta, b.angleBeta),
                         std::max(a.angleGamma, b.angleGamma));
}
