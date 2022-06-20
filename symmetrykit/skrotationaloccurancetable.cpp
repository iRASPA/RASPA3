module;

module skrotationaloccurancetable;

import skrotationmatrix;

SKRotationalOccuranceTable::SKRotationalOccuranceTable(int axis_6m, int axis_4m, int axis_3m,
    int axis_2m, int axis_1m, int axis_1,
    int axis_2, int axis_3, int axis_4, int axis_6)
{
    occurance[static_cast<typename std::underlying_type<SKRotationMatrix::RotationType>::type>(SKRotationMatrix::RotationType::axis_6m)] = axis_6m;
    occurance[static_cast<typename std::underlying_type<SKRotationMatrix::RotationType>::type>(SKRotationMatrix::RotationType::axis_4m)] = axis_4m;
    occurance[static_cast<typename std::underlying_type<SKRotationMatrix::RotationType>::type>(SKRotationMatrix::RotationType::axis_3m)] = axis_3m;
    occurance[static_cast<typename std::underlying_type<SKRotationMatrix::RotationType>::type>(SKRotationMatrix::RotationType::axis_2m)] = axis_2m;
    occurance[static_cast<typename std::underlying_type<SKRotationMatrix::RotationType>::type>(SKRotationMatrix::RotationType::axis_1m)] = axis_1m;
    occurance[static_cast<typename std::underlying_type<SKRotationMatrix::RotationType>::type>(SKRotationMatrix::RotationType::axis_1)] = axis_1;
    occurance[static_cast<typename std::underlying_type<SKRotationMatrix::RotationType>::type>(SKRotationMatrix::RotationType::axis_2)] = axis_2;
    occurance[static_cast<typename std::underlying_type<SKRotationMatrix::RotationType>::type>(SKRotationMatrix::RotationType::axis_3)] = axis_3;
    occurance[static_cast<typename std::underlying_type<SKRotationMatrix::RotationType>::type>(SKRotationMatrix::RotationType::axis_4)] = axis_4;
    occurance[static_cast<typename std::underlying_type<SKRotationMatrix::RotationType>::type>(SKRotationMatrix::RotationType::axis_6)] = axis_6;
}