export module float2;

export union float2
{
  float v[2];
  struct
  {
    float x, y;
  };

  inline float2(float x = 0, float y = 0) : x(x), y(y) {};

  inline float& operator[](int i) { return v[i]; }
  inline const float& operator[](int i) const { return v[i]; }

  inline float2& operator+=(const float2& rhs)
  {
    x += rhs.x;
    y += rhs.y;
    return *this;
  }
  inline float2& operator-=(const float2& rhs)
  {
    x -= rhs.x;
    y -= rhs.y;
    return *this;
  }

  inline static float dot(const float2& v1, const float2& v2) { return v1.x * v2.x + v1.y * v2.y; }
};

export inline float2 operator+(const float2& a, const float2& b) { return float2(a.x + b.x, a.y + b.y); }

export inline float2 operator-(const float2& a, const float2& b) { return float2(a.x - b.x, a.y - b.y); }

export inline float2 operator*(const float2& a, const float2& b) { return float2(a.x * b.x, a.y * b.y); }
