export module double2;

export union double2
{
  double v[2];
  struct
  {
    double x, y;
  };

  inline double2(double x = 0, double y = 0) : x(x), y(y) {};
  inline double& operator[](int i) { return v[i]; }
  inline const double& operator[](int i) const { return v[i]; }

  inline double2& operator+=(const double2& rhs)
  {
    x += rhs.x;
    y += rhs.y;
    return *this;
  }
  inline double2& operator-=(const double2& rhs)
  {
    x -= rhs.x;
    y -= rhs.y;
    return *this;
  }

  inline static double dot(const double2& v1, const double2& v2) { return v1.x * v2.x + v1.y * v2.y; }
};

export inline double2 operator+(const double2& a, const double2& b) { return double2(a.x + b.x, a.y + b.y); }

export inline double2 operator-(const double2& a, const double2& b) { return double2(a.x - b.x, a.y - b.y); }