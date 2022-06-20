export module skintegersymmetryoperationset;

import skseitzintegermatrix;
import <set>;
import <unordered_set>;
import <vector>;
import <tuple>;
import mathkit;
import skdefinitions;


export struct SKIntegerSymmetryOperationSet
{
	std::unordered_set<SKSeitzIntegerMatrix, SKSeitzIntegerMatrix::hashFunction> operations;
	Centring centring;

	SKIntegerSymmetryOperationSet();
	//SKIntegerSymmetryOperationSet(std::unordered_set<SKSeitzIntegerMatrix> &operations);
	SKIntegerSymmetryOperationSet(std::vector<SKSeitzIntegerMatrix> operations);

	inline size_t size() { return operations.size(); }
	SKIntegerSymmetryOperationSet fullSeitzMatrices();

	std::vector<std::tuple<double3, int, double>> symmetrize(double3x3 lattice, std::vector<std::tuple<double3, int, double>> atoms, double symmetryPrecision);
	std::vector<std::tuple<double3, int, double>> asymmetricAtoms(int HallNumber, std::vector<std::tuple<double3, int, double>>& atoms, double3x3 lattice, bool allowPartialOccupancies, double symmetryPrecision);
};