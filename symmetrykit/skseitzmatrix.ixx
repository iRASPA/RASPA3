export module skseitzmatrix;

import <string>;
import <vector>;
import mathkit;
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

        return (this->rotation == b.rotation) &&
            (dr.length_squared() < 1e-5);
    }

};

