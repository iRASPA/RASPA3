module;

#include <cmath>

module double3;

import <iostream>;

double3 double3::normalise()
{
    double magnitude = sqrt((x * x) + (y * y) + (z * z));

    if (magnitude != 0)
    {
        x /= magnitude;
        y /= magnitude;
        z /= magnitude;
    }
    return *this;
}

double3 double3::fract() const
{
    double3 s = double3(x, y, z);
    s.x -= std::rint(x);
    s.y -= std::rint(y);
    s.z -= std::rint(z);

    if (s.x < 0.0)
    {
        s.x += 1.0;
    }
    if (s.x > 1.0)
    {
        s.x -= 1.0;
    }

    if (s.y < 0.0)
    {
        s.y += 1.0;
    }
    if (s.y > 1.0)
    {
        s.y -= 1.0;
    }

    if (s.z < 0.0)
    {
        s.z += 1.0;
    }
    if (s.z > 1.0)
    {
        s.z -= 1.0;
    }
    return s;
}

double3 double3::randomVectorOnUnitSphere()
{
    double ran1, ran2, ranh, ransq = 0.0;

    do
    {
        ran1 = 2.0 * (double(rand()) / RAND_MAX) - 1.0;
        ran2 = 2.0 * (double(rand()) / RAND_MAX) - 1.0;
        ransq = ran1 * ran1 + ran2 * ran2;
    } while (ransq >= 1.0);

    ranh = 2.0 * sqrt(1.0 - ransq);
    return double3(ran1 * ranh, ran2 * ranh, 1.0 - 2.0 * ransq);
}

std::ostream& operator<<(std::ostream& out, const double3& vec) // output
{
    out << vec.x;
    out << ",";
    out << vec.y;
    out << ",";
    out << vec.z;
    return out;
}
