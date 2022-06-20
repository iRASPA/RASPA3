module;

#include <cmath>

module double4;

void double4::normalise()
{
    double magnitude = sqrt((x * x) + (y * y) + (z * z) * (w * w));

    if (magnitude != 0)
    {
        x /= magnitude;
        y /= magnitude;
        z /= magnitude;
        w /= magnitude;
    }
}