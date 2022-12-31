module;

module simulationbox;

import <cmath>;
import <numbers>;
import <string>;
import <iostream>;
import <ostream>;
import <sstream>;

import randomnumbers;
import double3x3;
import double3;
import units;
import print;

SimulationBox::SimulationBox(double a, double b, double c, double alpha, double beta, double gamma):
    lengthA(a), lengthB(b), lengthC(c), angleAlpha(alpha), angleBeta(beta), angleGamma(gamma)
{
    double temp = (cos(alpha) - cos(gamma) * cos(beta)) / sin(gamma);

    double3 v1 = double3(a, 0.0, 0.0);
    double3 v2 = double3(b * cos(gamma), b * sin(gamma), 0.0);
    double3 v3 = double3(c * cos(beta), c * temp, c * sqrt(1.0 - cos(beta) * cos(beta) - temp * temp));
    unitCell = double3x3(v1, v2, v3);
    if(a != 0.0 && b != 0.0 && c != 0.0)
    {
      inverseUnitCell = unitCell.inverse();
      volume = unitCell.determinant();
    }
    else
    {
      volume = 0.0;
    }
}

double3 SimulationBox::randomPosition() const
{
    return double3(unitCell.ax * RandomNumber::Uniform(),
                   unitCell.by * RandomNumber::Uniform(),
                   unitCell.cz * RandomNumber::Uniform());
}

void SimulationBox::setBoxLengths(double3 lengths)
{
    lengthA = lengths.x;
    lengthB = lengths.y;
    lengthC = lengths.z;
    double temp = (cos(angleAlpha) - cos(angleGamma) * cos(angleBeta)) / sin(angleGamma);
    double3 v1 = double3(lengths.x, 0.0, 0.0);
    double3 v2 = double3(lengths.y * cos(angleGamma), lengths.y * sin(angleGamma), 0.0);
    double3 v3 = double3(lengths.z * cos(angleBeta), lengths.z * temp, lengths.z * sqrt(1.0 - cos(angleBeta) * cos(angleBeta) - temp * temp));
    unitCell = double3x3(v1, v2, v3);
    inverseUnitCell = unitCell.inverse();
    volume = unitCell.determinant();
}

void SimulationBox::setBoxAngles(double3 angles)
{
    angleAlpha = angles.x;
    angleBeta = angles.y;
    angleGamma = angles.z;
    double3 lengths = this->lengths();
    lengthA = lengths.x;
    lengthB = lengths.y;
    lengthC = lengths.z;
    double temp = (cos(angles.x) - cos(angles.z) * cos(angles.y)) / sin(angles.z);
    double3 v1 = double3(lengths.x, 0.0, 0.0);
    double3 v2 = double3(lengths.y * cos(angles.z), lengths.y * sin(angles.z), 0.0);
    double3 v3 = double3(lengths.z * cos(angles.y), lengths.z * temp, lengths.z * sqrt(1.0 - cos(angles.y) * cos(angles.y) - temp * temp));
    unitCell = double3x3(v1, v2, v3);
    inverseUnitCell = unitCell.inverse();
    volume = unitCell.determinant();
}


double3 SimulationBox::lengths()
{
    return double3(unitCell[0].length(), unitCell[1].length(), unitCell[2].length());
}

double3 SimulationBox::angles()
{
    double3 column1 = unitCell[0];
    double3 column2 = unitCell[1];
    double3 column3 = unitCell[2];
    double length1 = column1.length();
    double length2 = column2.length();
    double length3 = column3.length();

    return double3(acos(double3::dot(column2, column3) / (length2 * length3)),
                   acos(double3::dot(column1, column3) / (length1 * length3)),
                   acos(double3::dot(column1, column2) / (length1 * length2)));
}

void SimulationBox::printParameters(std::ostream &stream) const
{
  std::print(stream, "Simulation parameters:\n");
  std::print(stream, "===============================================================================\n\n");
  std::print(stream, "Temperature:     {:9.5f} \n", temperature);
  std::print(stream, "Beta:            {:9.5f} \n", Beta);
  std::print(stream, "Pressure:        {:9.5f} \n", Units::PressureConversionFactor * pressure);
  std::print(stream, "\n\n");
}

void SimulationBox::printStatus(std::ostream &stream) const
{
  std::print(stream, "Box:     {:9.5f} {:9.5f} {:9.5f}\n", unitCell.ax, unitCell.bx, unitCell.cx);
  std::print(stream, "         {:9.5f} {:9.5f} {:9.5f}\n", unitCell.ay, unitCell.by, unitCell.cy);
  std::print(stream, "         {:9.5f} {:9.5f} {:9.5f}\n", unitCell.az, unitCell.bz, unitCell.cz);
  std::print(stream, "Lengths: {:9.5f} {:9.5f} {:9.5f}\n", lengthA, lengthB, lengthC);
  double conv = 180.0 / std::numbers::pi;
  std::print(stream, "Angles:  {:9.5f} {:9.5f} {:9.5f}\n", conv * angleAlpha, conv * angleBeta, conv * angleGamma);
}

std::string SimulationBox::printStatus(const SimulationBox& average, [[maybe_unused]] const SimulationBox& error) const
{
    std::ostringstream stream;
    std::print(stream, "Box:     {:9.5f} {:9.5f} {:9.5f}  Average: {:9.5f} {:9.5f} {:9.5f}\n", unitCell.ax, unitCell.bx, unitCell.cx,
        average.unitCell.ax, average.unitCell.bx, average.unitCell.cx);
    std::print(stream, "         {:9.5f} {:9.5f} {:9.5f}           {:9.5f} {:9.5f} {:9.5f}\n", unitCell.ay, unitCell.by, unitCell.cy,
        average.unitCell.ay, average.unitCell.by, average.unitCell.cy);
    std::print(stream, "         {:9.5f} {:9.5f} {:9.5f}           {:9.5f} {:9.5f} {:9.5f}\n", unitCell.az, unitCell.bz, unitCell.cz,
        average.unitCell.az, average.unitCell.bz, average.unitCell.cz);
    std::print(stream, "Lengths: {:9.5f} {:9.5f} {:9.5f}  Average: {:9.5f} {:9.5f} {:9.5f}\n", lengthA, lengthB, lengthC,
        average.lengthA, average.lengthB, average.lengthC);
    double conv = 180.0 / std::numbers::pi;
    std::print(stream, "Angles:  {:9.5f} {:9.5f} {:9.5f}  Average: {:9.5f} {:9.5f} {:9.5f}\n", 
        conv * angleAlpha, conv * angleBeta, conv * angleGamma,
         conv * average.angleAlpha, conv * average.angleBeta, conv * average.angleGamma);

    return stream.str();
}
