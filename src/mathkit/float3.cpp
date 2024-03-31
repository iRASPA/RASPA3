module;

module float3;

import <cmath>;
import <iostream>;

float3 float3::normalise()
{
    float magnitude = sqrt((x * x) + (y * y) + (z * z));

    if (magnitude != 0.0f)
    {
        x /= magnitude;
        y /= magnitude;
        z /= magnitude;
    }
    return *this;
}

float3 float3::fract()
{
    float3 s = float3(x, y, z);
    s.x -= rint(x);
    s.y -= rint(y);
    s.z -= rint(z);

    if (s.x < 0.0f)
    {
        s.x += 1.0f;
    }
    if (s.x > 1.0f)
    {
        s.x -= 1.0f;
    }

    if (s.y < 0.0f)
    {
        s.y += 1.0f;
    }
    if (s.y > 1.0f)
    {
        s.y -= 1.0f;
    }

    if (s.z < 0.0f)
    {
        s.z += 1.0f;
    }
    if (s.z > 1.0f)
    {
        s.z -= 1.0f;
    }
    return s;
}

std::ostream& operator<<(std::ostream& out, const float3& vec) // output
{
    out << vec.x;
    out << vec.y;
    out << vec.z;
    return out;
}

