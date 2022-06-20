export module skrotationaloccurancetable;

import <map>;


export struct SKRotationalOccuranceTable
{
    std::map<int, int> occurance;

    SKRotationalOccuranceTable(int axis_6m, int axis_4m, int axis_3m, int axis_2m, int axis_1m, int axis_1, int axis_2, int axis_3, int axis_4, int axis_6);

    inline bool operator==(const SKRotationalOccuranceTable& b) const
    {
        return (this->occurance == b.occurance);
    }
};