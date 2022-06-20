module;

module skseitzmatrix;

SKSeitzMatrix::SKSeitzMatrix()
{

}

SKSeitzMatrix::SKSeitzMatrix(SKRotationMatrix rotation, double3 translation)
{
	this->rotation = rotation;
	this->translation = translation;
}