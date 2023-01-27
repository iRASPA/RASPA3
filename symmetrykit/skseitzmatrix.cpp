module;

module skseitzmatrix;

import int3;
import int3x3;
import double3;
import double3x3;

import skrotationmatrix;

SKSeitzMatrix::SKSeitzMatrix()
{

}

SKSeitzMatrix::SKSeitzMatrix(SKRotationMatrix rotation, double3 translation)
{
    this->rotation = rotation;
    this->translation = translation;
}
