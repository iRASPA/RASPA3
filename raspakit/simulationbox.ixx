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
import <ostream>;

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

    explicit SimulationBox(double a, double b, double c);
    explicit SimulationBox(double a, double b, double c, double alpha, double beta, double gamma);
        

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

    void printParameters(std::ostream &stream) const;
    void printStatus(std::ostream &stream) const;
    std::string printStatus(const SimulationBox& average, const SimulationBox& error) const;

    double lengthA;
    double lengthB;
    double lengthC;
    double angleAlpha;
    double angleBeta;
    double angleGamma;
    double3x3 unitCell;
    double3x3 inverseUnitCell;
    double volume;
    double temperature{ 300.0 };
    double pressure{ 1e4 };
    double input_pressure{ 1e4 };
    double Beta{ 1.0 / (Units::KB * 300.0)};
    double alpha{ 0.265058 };
    int3 kmax{ 11, 11, 7 };
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

    void scale(double scale)
    {
      lengthA *= scale;
      lengthB *= scale;
      lengthC *= scale;

      double temp = (cos(angleAlpha) - cos(angleGamma) * cos(angleBeta)) / sin(angleGamma);
      double3 v1 = double3(lengthA, 0.0, 0.0);
      double3 v2 = double3(lengthB * cos(angleGamma), lengthB * sin(angleGamma), 0.0);
      double3 v3 = double3(lengthC * cos(angleBeta), lengthC * temp, lengthC * sqrt(1.0 - cos(angleBeta) * cos(angleBeta) - temp * temp));
      unitCell = double3x3(v1, v2, v3);
      inverseUnitCell = unitCell.inverse();
      volume = unitCell.determinant();
    }

    void scale(int3 scale)
    {
        lengthA *= static_cast<double>(scale.x);
        lengthB *= static_cast<double>(scale.y);
        lengthC *= static_cast<double>(scale.z);

        double temp = (cos(angleAlpha) - cos(angleGamma) * cos(angleBeta)) / sin(angleGamma);
        double3 v1 = double3(lengthA, 0.0, 0.0);
        double3 v2 = double3(lengthB * cos(angleGamma), lengthB * sin(angleGamma), 0.0);
        double3 v3 = double3(lengthC * cos(angleBeta), lengthC * temp, lengthC * sqrt(1.0 - cos(angleBeta) * cos(angleBeta) - temp * temp));
        unitCell = double3x3(v1, v2, v3);
        inverseUnitCell = unitCell.inverse();
        volume = unitCell.determinant();
    }

    SimulationBox scaled(double scale)
    {
      SimulationBox v;

      v.lengthA = scale * lengthA;
      v.lengthB = scale * lengthB;
      v.lengthC = scale * lengthC;
      v.angleAlpha = angleAlpha;
      v.angleBeta = angleBeta;
      v.angleGamma = angleGamma;

      double temp = (cos(v.angleAlpha) - cos(v.angleGamma) * cos(v.angleBeta)) / sin(v.angleGamma);
      double3 v1 = double3(v.lengthA, 0.0, 0.0);
      double3 v2 = double3(v.lengthB * cos(v.angleGamma), v.lengthB * sin(v.angleGamma), 0.0);
      double3 v3 = double3(v.lengthC * cos(v.angleBeta), v.lengthC * temp, v.lengthC * sqrt(1.0 - cos(v.angleBeta) * cos(v.angleBeta) - temp * temp));
      v.unitCell = double3x3(v1, v2, v3);
      v.inverseUnitCell = v.unitCell.inverse();
      v.volume = v.unitCell.determinant();

      v.temperature = temperature;
      v.pressure = pressure;
      v.Beta = Beta;
      v.alpha = alpha;
      v.kmax = kmax;
      v.type = type;

      return v;
    }

    SimulationBox scaled(int3 scale)
    {
        SimulationBox v;

        v.lengthA = static_cast<double>(scale.x) * lengthA;
        v.lengthB = static_cast<double>(scale.y) * lengthB;
        v.lengthC = static_cast<double>(scale.z) * lengthC;
        v.angleAlpha = angleAlpha;
        v.angleBeta = angleBeta;
        v.angleGamma = angleGamma;

        double temp = (cos(v.angleAlpha) - cos(v.angleGamma) * cos(v.angleBeta)) / sin(v.angleGamma);
        double3 v1 = double3(v.lengthA, 0.0, 0.0);
        double3 v2 = double3(v.lengthB * cos(v.angleGamma), v.lengthB * sin(v.angleGamma), 0.0);
        double3 v3 = double3(v.lengthC * cos(v.angleBeta), v.lengthC * temp, v.lengthC * sqrt(1.0 - cos(v.angleBeta) * cos(v.angleBeta) - temp * temp));
        v.unitCell = double3x3(v1, v2, v3);
        v.inverseUnitCell = v.unitCell.inverse();
        v.volume = v.unitCell.determinant();

        v.temperature = temperature;
        v.pressure = pressure;
        v.Beta = Beta;
        v.alpha = alpha;
        v.kmax = kmax;
        v.type = type;

        return v;
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

    m.lengthA = std::sqrt(a.lengthA);
    m.lengthB = std::sqrt(a.lengthB);
    m.lengthC = std::sqrt(a.lengthC);
    m.angleAlpha = std::sqrt(a.angleAlpha);
    m.angleBeta = std::sqrt(a.angleBeta);
    m.angleGamma = std::sqrt(a.angleGamma);
    m.unitCell = sqrt(a.unitCell);
    m.volume = std::sqrt(a.volume);
    m.temperature = std::sqrt(a.temperature);
    m.pressure = std::sqrt(a.pressure);

    return m;
}

export inline SimulationBox max(const SimulationBox& a, const SimulationBox& b)
{
    SimulationBox c = SimulationBox(std::max(a.lengthA, b.lengthA),
                         std::max(a.lengthB, b.lengthB),
                         std::max(a.lengthC, b.lengthC),
                         std::max(a.angleAlpha, b.angleAlpha),
                         std::max(a.angleBeta, b.angleBeta),
                         std::max(a.angleGamma, b.angleGamma));
    c.temperature = a.temperature;
    c.Beta = a.Beta;
    c.pressure = a.pressure;
    c.alpha = a.alpha;
    c.kmax = a.kmax;

    return c;
}
