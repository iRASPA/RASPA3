module;

#ifdef USE_LEGACY_HEADERS
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>
#endif

export module skseitzmatrix;

#ifndef USE_LEGACY_HEADERS
import <string>;
import <vector>;
import <cmath>;
#endif

import int3;
import int3x3;
import double3;
import double3x3;
import skrotationmatrix;
import skonethirdseitzmatrix;

export struct SKSeitzMatrix
{
  SKRotationMatrix rotation;
  double3 translation;

  SKSeitzMatrix();
  SKSeitzMatrix(SKRotationMatrix rotation, double3 translation);

  inline bool operator==(const SKSeitzMatrix& b)
  {
    double3 dr = (this->translation - b.translation);
    dr.x -= rint(dr.x);
    dr.y -= rint(dr.y);
    dr.z -= rint(dr.z);

    return (this->rotation == b.rotation) && (dr.length_squared() < 1e-5);
  }
};
