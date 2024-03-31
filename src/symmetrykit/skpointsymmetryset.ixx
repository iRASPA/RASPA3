export module skpointsymmetryset;

import <vector>;
import skrotationmatrix;

export class SKPointSymmetrySet
{
public:
    SKPointSymmetrySet(std::vector<SKRotationMatrix> rotations);
    const std::vector<SKRotationMatrix>& rotations() { return _rotations; }
private:
    std::vector<SKRotationMatrix> _rotations;
};