export module skseitzmatrix;

import <string>;
import <vector>;
import int3;
import int3x3;
import double3;
import double3x3;
import skrotationmatrix;
import skonethirdseitzmatrix;

import <cmath>;

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

        return (this->rotation == b.rotation) &&
            (dr.length_squared() < 1e-5);
    }

};

