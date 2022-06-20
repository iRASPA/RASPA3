module;

#include <cmath>

module float4;

void float4::normalise()
{
    float magnitude = sqrt((x * x) + (y * y) + (z * z) * (w * w));

    if (magnitude != 0)
    {
        x /= magnitude;
        y /= magnitude;
        z /= magnitude;
        w /= magnitude;
    }
}