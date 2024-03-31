export module float4;

import double4;
import double3;
import float3;

export union float4
{
    float v[4];
    struct { float x, y, z, w; };

    float4(float x = 0, float y = 0, float z = 0, float w = 0) :x(x), y(y), z(z), w(w) {}
    float4(double4 a) :x(float(a.x)), y(float(a.y)), z(float(a.z)), w(float(a.w)) {}
    float4(double3 a, double b) :x(float(a.x)), y(float(a.y)), z(float(a.z)), w(float(b)) {}
    //float4(QColor a, double b) :x(float(a.redF())), y(float(a.greenF())), z(float(a.blueF())), w(float(b)) {}
    float4(float3 a, float b) :x(a.x), y(a.y), z(a.z), w(b) {}

    inline float& operator [] (int i) { return v[i]; }
    inline const float& operator [] (int i) const { return v[i]; }

    void normalise();
    inline static float dot(const float4& a, const float4& b) { return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w; }

    // http://www.gamedev.net/topic/456301-cross-product-vector-4d/
    inline static float4 cross(const float4& a, const float4& b) { return float4(a.y * b.z - a.z * b.y, -a.x * b.z + a.z * b.x, a.x * b.y - a.y * b.x, 0.0); }
    inline static float4 cross(const float4& a, const float4& b, const float4& c)
    {
        return float4(a.y * (b.z * c.w - c.z * b.w) - a.z * (b.y * c.w - c.y * b.w) + a.w * (b.y * c.z - c.y * b.z),
            -a.x * (b.z * c.w - c.z * b.w) + a.z * (b.x * c.w - c.x * b.w) - a.w * (b.x * c.z - c.x * b.z),
            a.x * (b.y * c.w - c.y * b.w) - a.y * (b.x * c.w - c.x * b.w) + a.w * (b.x * c.y - c.x * b.y),
            -a.x * (b.y * c.z - c.y * b.z) + a.y * (b.x * c.z - c.x * b.z) - a.z * (b.x * c.y - c.x * b.y));
    }

    float4 operator-() const { return float4(-this->x, -this->y, -this->z, -this->w); }
    float4& operator+=(const float4& b) { this->x += b.x, this->y += b.y, this->z += b.z; this->w += b.w; return *this; }
    float4& operator-=(const float4& b) { this->x -= b.x, this->y -= b.y, this->z -= b.z; this->w += b.w; return *this; }
};

export inline float4 operator+(const float4& a, const float4& b)
{
    return float4(a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w);
}

inline float4 operator-(const float4& a, const float4& b)
{
    return float4(a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w);
}

export inline float4 operator*(const float4& a, const float4& b)
{
    return float4(a.x * b.x, a.y * b.y, a.z * b.z, a.w * b.w);
}

export inline float4 operator*(const float4& a, const double& b)
{
    return float4(a.x * static_cast<float>(b), a.y * static_cast<float>(b), a.z * static_cast<float>(b), a.w * static_cast<float>(b));
}

inline float4 operator*(const double& a, float4 const& b)
{
    return float4(static_cast<float>(a) * b.x, static_cast<float>(a) * b.y, static_cast<float>(a) * b.z, static_cast<float>(a) * b.w);
}

export inline float4 operator*(const float4& a, const float& b)
{
    return float4(a.x * b, a.y * b, a.z * b, a.w * b);
}

export inline float4 operator*(const float& a, float4 const& b)
{
    return float4(a * b.x, a * b.y, a * b.z, a * b.w);
}
