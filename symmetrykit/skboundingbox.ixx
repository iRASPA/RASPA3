export module skboundingbox;

import double3;
import double4x4;
import <array>;
import <cstdlib>;

export class SKBoundingBox
{
public:
    SKBoundingBox();
    SKBoundingBox(double3 minimum, double3 maximum);
    SKBoundingBox(const double3 center, const double3 width, const double scale);
    double3 widths() const;
    std::array<double3, 8> const corners() const;
    std::array<std::pair<double3, double3>, 12> const sides() const;
    double3 center();
    double volume();
    double shortestEdge();
    double longestEdge();
    double3 aspectRatio();
    double boundingSphereRadius();
    double3 maximum() const { return _maximum; }
    double3 minimum() const { return _minimum; }
    SKBoundingBox adjustForTransformation(double4x4 transformation);
    friend SKBoundingBox operator+(const SKBoundingBox left, const SKBoundingBox right);
    friend SKBoundingBox operator+(const SKBoundingBox left, double3 right);
    friend SKBoundingBox operator-(const SKBoundingBox left, double3 right);
private:
    int64_t _versionNumber{ 1 };
    double3 _minimum = double3(0.0, 0.0, 0.0);
    double3 _maximum = double3(0.0, 0.0, 0.0);
};