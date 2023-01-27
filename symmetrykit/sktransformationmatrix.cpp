module;

module sktransformationmatrix;

import int3;
import int3x3;
import double3;
import double3x3;

SKTransformationMatrix::SKTransformationMatrix() : transformation(), translation(int3(0, 0, 0))
{

}

SKTransformationMatrix::SKTransformationMatrix(int3x3 m) : transformation(m), translation(int3(0, 0, 0))
{
}

SKTransformationMatrix::SKTransformationMatrix(int3 v1, int3 v2, int3 v3)
{
	this->transformation.m11 = v1.x; this->transformation.m21 = v1.y; this->transformation.m31 = v1.z;
	this->transformation.m12 = v2.x; this->transformation.m22 = v2.y; this->transformation.m32 = v2.z;
	this->transformation.m13 = v3.x; this->transformation.m23 = v3.y; this->transformation.m33 = v3.z;
}

SKTransformationMatrix SKTransformationMatrix::zero = SKTransformationMatrix(int3(0, 0, 0), int3(0, 0, 0), int3(0, 0, 0));
SKTransformationMatrix SKTransformationMatrix::identity = SKTransformationMatrix(int3(1, 0, 0), int3(0, 1, 0), int3(0, 0, 1));
SKTransformationMatrix SKTransformationMatrix::inversionIdentity = SKTransformationMatrix(int3(-1, 0, 0), int3(0, -1, 0), int3(0, 0, -1));

// based on the centering, convert conventional cell to primitive using conventionally used transformation matrices
// Taken from: Table 2.C.1, page 141, Fundamentals of Crystallography, 2nd edition, C. Giacovazzo et al. 2002
// Tranformation matrices M, conventionally used to generate centered from primitive lattices, and vice versa, accoording to: A' = M A

SKTransformationMatrix SKTransformationMatrix::primitiveToPrimitive = SKTransformationMatrix(int3(1, 0, 0), int3(0, 1, 0), int3(0, 0, 1));           // P -> P
SKTransformationMatrix SKTransformationMatrix::primitiveToBodyCentered = SKTransformationMatrix(int3(0, 1, 1), int3(1, 0, 1), int3(1, 1, 0));           // P -> I
SKTransformationMatrix SKTransformationMatrix::primitiveToFaceCentered = SKTransformationMatrix(int3(-1, 1, 1), int3(1, -1, 1), int3(1, 1, -1));        // P -> F
SKTransformationMatrix SKTransformationMatrix::primitiveToACentered = SKTransformationMatrix(int3(-1, 0, 0), int3(0, -1, 1), int3(0, 1, 1));            // P -> A
SKTransformationMatrix SKTransformationMatrix::primitiveToBCentered = SKTransformationMatrix(int3(-1, 0, 1), int3(0, -1, 0), int3(1, 0, 1));            // P -> B
SKTransformationMatrix SKTransformationMatrix::primitiveToCCentered = SKTransformationMatrix(int3(1, 1, 0), int3(1, -1, 0), int3(0, 0, -1));            // P -> C
SKTransformationMatrix SKTransformationMatrix::primitiveToRhombohedral = SKTransformationMatrix(int3(1, -1, 0), int3(0, 1, -1), int3(1, 1, 1));       // P -> R
SKTransformationMatrix SKTransformationMatrix::primitiveToHexagonal = SKTransformationMatrix(int3(1, -1, 0), int3(1, 2, 0), int3(0, 0, 3));          // P -> H
SKTransformationMatrix SKTransformationMatrix::rhombohedralObverseHexagonal = SKTransformationMatrix(int3(1, 0, 1), int3(-1, 1, 1), int3(0, -1, 1));    // Robv -> Rh
SKTransformationMatrix SKTransformationMatrix::rhombohedralHexagonalToReverse = SKTransformationMatrix(int3(1, 1, -2), int3(-1, 0, 1), int3(1, 1, -1));   // Rh -> Rrev
SKTransformationMatrix SKTransformationMatrix::rhombohedralReverseToHexagonal = SKTransformationMatrix(int3(-1, -1, 1), int3(0, 1, 1), int3(-1, 0, 1)); // Rrev -> Rh

SKTransformationMatrix SKTransformationMatrix::monoclinicAtoC = SKTransformationMatrix(int3(0, 0, 1), int3(0, -1, 0), int3(1, 0, 0));
SKTransformationMatrix SKTransformationMatrix::AtoC = SKTransformationMatrix(int3(0, 1, 0), int3(0, 0, 1), int3(1, 0, 0));
SKTransformationMatrix SKTransformationMatrix::monoclinicItoC = SKTransformationMatrix(int3(1, 0, 1), int3(0, 1, 0), int3(-1, 0, 0));
SKTransformationMatrix SKTransformationMatrix::BtoC = SKTransformationMatrix(int3(0, 0, 1), int3(1, 0, 0), int3(0, 1, 0));
SKTransformationMatrix SKTransformationMatrix::primitiveRhombohedralToTripleHexagonalCell_R2 = SKTransformationMatrix(int3(0, -1, 1), int3(1, 0, -1), int3(1, 1, 1));
SKTransformationMatrix SKTransformationMatrix::primitiveRhombohedralToTripleHexagonalCell_R1_Obverse = SKTransformationMatrix(int3(1, -1, 0), int3(0, 1, -1), int3(1, 1, 1));

double3x3 SKTransformationMatrix::bodyCenteredToPrimitive = double3x3(double3(-0.5, 0.5, 0.5), double3(0.5, -0.5, 0.5), double3(0.5, 0.5, -0.5));  // I -> P
double3x3 SKTransformationMatrix::faceCenteredToPrimitive = double3x3(double3(0, 0.5, 0.5), double3(0.5, 0, 0.5), double3(0.5, 0.5, 0));           // F -> P
double3x3 SKTransformationMatrix::ACenteredToPrimitive = double3x3(double3(-1.0, 0, 0), double3(0, -0.5, 0.5), double3(0, 0.5, 0.5));              // A -> P
double3x3 SKTransformationMatrix::BCenteredToPrimitive = double3x3(double3(-0.5, 0, 0.5), double3(0, -1.0, 0), double3(0.5, 0, 0.5));              // B -> P
double3x3 SKTransformationMatrix::CCenteredToPrimitive = double3x3(double3(0.5, 0.5, 0), double3(0.5, -0.5, 0), double3(0, 0, -1.0));              // C -> P
double3x3 SKTransformationMatrix::rhombohedralToPrimitive = double3x3(double3(2.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0), double3(-1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0), double3(-1.0 / 3.0, -2.0 / 3.0, 1.0 / 3.0));  // R -> P

// CHECK
double3x3 SKTransformationMatrix::rhombohedralReverseToPrimitive = double3x3(double3(1.0 / 3.0, 2.0 / 3.0, 1.0 / 3.0), double3(-2.0 / 3.0, -1.0 / 3.0, 1.0 / 3.0), double3(1.0 / 3.0, -1.0 / 3.0, 1.0 / 3.0));  // R -> P

double3x3 SKTransformationMatrix::hexagonalToPrimitive = double3x3(double3(2.0 / 3.0, 1.0 / 3.0, 0), double3(-1.0 / 3.0, 1.0 / 3.0, 0), double3(0, 0, 1.0 / 3.0));  // H -> P
double3x3 SKTransformationMatrix::rhombohedralHexagonalToObverse = double3x3(double3(2.0 / 3.0, -1.0 / 3.0, -1.0 / 3.0), double3(1.0 / 3.0, 1.0 / 3.0, -2.0 / 3.0), double3(1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0));   // Rh -> Robv

SKTransformationMatrix SKTransformationMatrix::monoclinicB1toA1 = SKTransformationMatrix(int3(0, 0, 1), int3(1, 0, 0), int3(0, 1, 0));
SKTransformationMatrix SKTransformationMatrix::monoclinicB1toA2 = SKTransformationMatrix(int3(0, 1, 0), int3(1, 0, 0), int3(0, -1, -1));
SKTransformationMatrix SKTransformationMatrix::monoclinicB1toA3 = SKTransformationMatrix(int3(0, -1, -1), int3(1, 0, 0), int3(0, 0, 1));
SKTransformationMatrix SKTransformationMatrix::monoclinicB1toB2 = SKTransformationMatrix(int3(0, 0, 1), int3(0, 1, 0), int3(-1, 0, -1));
SKTransformationMatrix SKTransformationMatrix::monoclinicB1toB3 = SKTransformationMatrix(int3(-1, 0, -1), int3(0, 1, 0), int3(1, 0, 0));
SKTransformationMatrix SKTransformationMatrix::monoclinicB1toC1 = SKTransformationMatrix(int3(0, 1, 0), int3(0, 0, 1), int3(1, 0, 0));
SKTransformationMatrix SKTransformationMatrix::monoclinicB1toC2 = SKTransformationMatrix(int3(1, 0, 0), int3(0, 0, 1), int3(-1, -1, 0));
SKTransformationMatrix SKTransformationMatrix::monoclinicB1toC3 = SKTransformationMatrix(int3(-1, -1, 0), int3(0, 0, 1), int3(0, 1, 0));

SKTransformationMatrix SKTransformationMatrix::orthorhombicCABtoABC = SKTransformationMatrix(int3(0, 1, 0), int3(0, 0, 1), int3(1, 0, 0));
SKTransformationMatrix SKTransformationMatrix::orthorhombicBCAtoABC = SKTransformationMatrix(int3(0, 0, 1), int3(1, 0, 0), int3(0, 1, 0));
SKTransformationMatrix SKTransformationMatrix::orthorhombicBAmCtoABC = SKTransformationMatrix(int3(0, 1, 0), int3(1, 0, 0), int3(0, 0, -1));
SKTransformationMatrix SKTransformationMatrix::orthorhombicAmCBtoABC = SKTransformationMatrix(int3(1, 0, 0), int3(0, 0, 1), int3(0, -1, 0));
SKTransformationMatrix SKTransformationMatrix::orthorhombicmCBAtoABC = SKTransformationMatrix(int3(0, 0, 1), int3(0, 1, 0), int3(-1, 0, 0));
