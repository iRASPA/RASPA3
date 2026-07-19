#pragma once

#include <string_view>

// Embedded molecule / force-field JSON used by raspa3 unit tests.
namespace molecule_fixtures
{
inline constexpr std::string_view kPropaneJson =
R"({
  "CriticalTemperature" : 369.825,
  "CriticalPressure" : 4247660.0,
  "AcentricFactor" : 0.1524,
  "Type" : "flexible",
  "pseudoAtoms" :
    [
      ["CH3", [0.0, 0.0, 0.0]],
      ["CH2", [0.0, 0.0, 0.0]],
      ["CH3", [0.0, 0.0, 0.0]]
    ],
  "Connectivity" : [
    [0, 1],
    [1, 2]
  ],
  "Bonds" : [
    [["CH3", "CH2"], "FIXED", [1.54]]
  ],
  "Bends" : [
    [["CH3", "CH2", "CH3"], "HARMONIC", [62500.0, 114]]
  ],
  "Partial-reinsertion" : [
    [0, 1],
    [1, 2],
    [0],
    [2]
  ]
}
)";

inline constexpr std::string_view kButaneJson =
R"({
  "CriticalTemperature" : 425.125,
  "CriticalPressure" : 3796000.0,
  "AcentricFactor" : 0.201,
  "Type" : "flexible",
  "pseudoAtoms" :
    [
      ["CH3", [0.0, 0.0, 0.0]],
      ["CH2", [0.0, 0.0, 0.0]],
      ["CH2", [0.0, 0.0, 0.0]],
      ["CH3", [0.0, 0.0, 0.0]]
    ],
  "Connectivity" : [
    [0, 1],
    [1, 2],
    [2, 3]
  ],
  "Bonds" : [
    [["CH3", "CH2"], "FIXED", [1.54]],
    [["CH2", "CH2"], "FIXED", [1.54]]
  ],
  "Bends" : [
    [["CH3", "CH2", "CH2"], "HARMONIC", [62500.0, 114]]
  ],
  "Torsions" : [
    [["CH3", "CH2", "CH2", "CH3"], "TRAPPE", [0.0, 355.03, -68.19, 791.32]]
  ],
  "Partial-reinsertion" : [
    [0, 1],
    [2, 3],
    [0],
    [3]
  ]
}

)";

inline constexpr std::string_view kEthaneJson =
R"({
  "CriticalTemperature" : 305.32,
  "CriticalPressure" : 4872000.0,
  "AcentricFactor" : 0.099,
  "Type" : "flexible",
  "pseudoAtoms" :
    [
      ["CH3", [0.0, 0.0, 0.0]],
      ["CH3", [0.0, 0.0, 0.0]]
    ],
  "Connectivity" : [
    [0, 1]
  ],
  "Bonds" : [
    [["CH3", "CH3"], "FIXED", [1.54]]
  ],
  "Partial-reinsertion" : [
    [0],
    [1]
  ]
}
)";

inline constexpr std::string_view kCO2Json =
R"({
  "CriticalTemperature": 304.1282,
  "CriticalPressure": 7377300.0,
  "AcentricFactor": 0.22394,
  "Type": "rigid",
  "pseudoAtoms": [
    ["O_co2", [0.0, 0.0, 1.149]],
    ["C_co2", [0.0, 0.0, 0.0]],
    ["O_co2", [0.0, 0.0, -1.149]]
  ]
}   
)";

inline constexpr std::string_view kN2Json =
R"({
  "CriticalTemperature": 126.192,
  "CriticalPressure": 3395800.0,
  "AcentricFactor": 0.0372,
  "Type": "rigid",
  "pseudoAtoms": [
    ["N_n2", [0.0, 0.0, 0.55]],
    ["N_com", [0.0, 0.0, 0.0]],
    ["N_n2", [0.0, 0.0, -0.55]]
  ]
}   
)";

inline constexpr std::string_view kCo2N2ForceFieldJson =
R"({
  "PseudoAtoms" :
  [
    {
      "name" : "C_co2",
      "framework" : false,
      "print_to_output" : true,
      "element" : "C",
      "print_as" : "C",
      "mass" : 12.0,
      "charge" :  0.6512,
      "source" : "A. Garcia-Sanchez et al., J. Phys. Chem. C 2009, 113, 8814-8820"
    },
    {
      "name" : "O_co2",
      "framework" : false,
      "print_to_output" : true,
      "element" : "O",
      "print_as" : "O",
      "mass" : 15.9994,
      "charge" : -0.3256,
      "source" : "A. Garcia-Sanchez et al., J. Phys. Chem. C 2009, 113, 8814-8820"
    },
    {
      "name" : "N_n2",
      "framework" : false,
      "print_to_output" : true,
      "element" : "N",
      "print_as" : "N",
      "mass" : 14.00674,
      "charge" : -0.405,
      "source" : "A. Martin-Calvo et al. , Phys. Chem. Chem. Phys. 2011, 13, 11165-11174"
    },
    {
      "name" : "N_com",
      "framework" : false,
      "print_to_output" : false,
      "element" : "N",
      "print_as" : "-",
      "mass" : 0.0,
      "charge" :  0.810,
      "source" : "A. Martin-Calvo et al. , Phys. Chem. Chem. Phys. 2011, 13, 11165-11174"
    }
  ],
  "SelfInteractions" : 
  [
    {
      "name" : "O_co2",
      "type" : "lennard-jones",    
      "parameters" : [85.671, 3.017],
      "source" : "A. Garcia-Sanchez et al., J. Phys. Chem. C 2009, 113, 8814-8820"
    },
    {
      "name" : "C_co2",
      "type" : "lennard-jones",    
      "parameters" : [29.933, 2.745],
      "source" : "A. Garcia-Sanchez et al., J. Phys. Chem. C 2009, 113, 8814-8820"
    },
    {
      "name" : "N_n2",
      "type" : "lennard-jones",    
      "parameters" : [38.298, 3.306],
      "source" : "A. Martin-Calvo et al. , Phys. Chem. Chem. Phys. 2011, 13, 11165-11174"
    },
    {
      "name" : "N_com",
      "type" : "none",    
      "parameters" : [0.0, 1.0],
      "source" : "A. Martin-Calvo et al. , Phys. Chem. Chem. Phys. 2011, 13, 11165-11174"
    }
  ],
  "MixingRule" : "Lorentz-Berthelot",
  "TruncationMethod" : "shifted",
  "TailCorrections" : false
}
)";

inline constexpr std::string_view kCyclohexaneJson =
R"({
  "CriticalTemperature" : 553.6,
  "CriticalPressure" : 4073000.0,
  "AcentricFactor" : 0.211,
  "Type" : "flexible",
  "pseudoAtoms" :
    [
      ["CH2_c", [1.465493, 0.000000, 0.236606]],
      ["CH2_c", [0.732747, 1.269154, -0.236606]],
      ["CH2_c", [-0.732747, 1.269154, 0.236606]],
      ["CH2_c", [-1.465493, 0.000000, -0.236606]],
      ["CH2_c", [-0.732747, -1.269154, 0.236606]],
      ["CH2_c", [0.732747, -1.269154, -0.236606]]
    ],
  "Connectivity" : [
    [0, 1],
    [1, 2],
    [2, 3],
    [3, 4],
    [4, 5],
    [5, 0]
  ],
  "Groups" : [
    { "Type" : "Cycle", "Atoms" : [0, 1, 2, 3, 4, 5] }
  ],
  "Bonds" : [
    [["CH2_c", "CH2_c"], "HARMONIC", [96500.0, 1.54]]
  ],
  "Bends" : [
    [["CH2_c", "CH2_c", "CH2_c"], "HARMONIC", [62500.0, 114.0]]
  ],
  "Torsions" : [
    [["CH2_c", "CH2_c", "CH2_c", "CH2_c"], "TRAPPE", [0.0, 355.03, -68.19, 791.32]]
  ],
  "VanDerWaals" : "auto"
}
)";

inline constexpr std::string_view kSemiFlexiblePentaneJson =
R"({
  "CriticalTemperature" : 469.7,
  "CriticalPressure" : 3370000.0,
  "AcentricFactor" : 0.251,
  "Type" : "flexible",
  "pseudoAtoms" :
    [
      ["CH3", [0.0, 0.0, 0.0]],
      ["CH2", [0.0, 0.0, 0.0]],
      ["CH2", [1.54, 0.0, 0.0]],
      ["CH2", [2.16637443033673, 1.40686000476961, 0.0]],
      ["CH3", [0.0, 0.0, 0.0]]
    ],
  "Connectivity" : [
    [0, 1],
    [1, 2],
    [2, 3],
    [3, 4]
  ],
  "Groups" : [
    { "Type" : "Flexible", "Atoms" : [0] },
    { "Type" : "Rigid",    "Atoms" : [1, 2, 3] },
    { "Type" : "Flexible", "Atoms" : [4] }
  ],
  "Bonds" : [
    [["CH3", "CH2"], "HARMONIC", [96500.0, 1.54]],
    [["CH2", "CH2"], "FIXED", [1.54]]
  ],
  "Bends" : [
    [["CH3", "CH2", "CH2"], "HARMONIC", [62500.0, 114]],
    [["CH2", "CH2", "CH2"], "FIXED", [114]]
  ],
  "Torsions" : [
    [["CH3", "CH2", "CH2", "CH2"], "TRAPPE", [0.0, 355.03, -68.19, 791.32]],
    [["CH2", "CH2", "CH2", "CH2"], "TRAPPE", [0.0, 355.03, -68.19, 791.32]]
  ],
  "VanDerWaals" : "auto"
}
)";

inline constexpr std::string_view kDiethylBiphenylJson =
R"({
  "CriticalTemperature" : 800.0,
  "CriticalPressure" : 3000000.0,
  "AcentricFactor" : 0.4,
  "Type" : "flexible",
  "pseudoAtoms" :
    [
      ["CH3",   [0.0, 0.0, 0.0]],
      ["CH2",   [0.0, 0.0, 0.0]],
      ["C_ar",  [-1.40, 0.0, 0.0]],
      ["CH_ar", [-0.70, 1.21243556529821, 0.0]],
      ["CH_ar", [0.70, 1.21243556529821, 0.0]],
      ["C_ar",  [1.40, 0.0, 0.0]],
      ["CH_ar", [0.70, -1.21243556529821, 0.0]],
      ["CH_ar", [-0.70, -1.21243556529821, 0.0]],
      ["C_ar",  [-1.40, 0.0, 0.0]],
      ["CH_ar", [-0.70, 1.21243556529821, 0.0]],
      ["CH_ar", [0.70, 1.21243556529821, 0.0]],
      ["C_ar",  [1.40, 0.0, 0.0]],
      ["CH_ar", [0.70, -1.21243556529821, 0.0]],
      ["CH_ar", [-0.70, -1.21243556529821, 0.0]],
      ["CH2",   [0.0, 0.0, 0.0]],
      ["CH3",   [0.0, 0.0, 0.0]]
    ],
  "Connectivity" : [
    [0, 1],
    [1, 2],
    [2, 3],
    [3, 4],
    [4, 5],
    [5, 6],
    [6, 7],
    [7, 2],
    [5, 8],
    [8, 9],
    [9, 10],
    [10, 11],
    [11, 12],
    [12, 13],
    [13, 8],
    [11, 14],
    [14, 15]
  ],
  "Groups" : [
    { "Type" : "Flexible", "Atoms" : [0, 1] },
    { "Type" : "Rigid",    "Atoms" : [2, 3, 4, 5, 6, 7] },
    { "Type" : "Rigid",    "Atoms" : [8, 9, 10, 11, 12, 13] },
    { "Type" : "Flexible", "Atoms" : [14, 15] }
  ],
  "Bonds" : [
    [["CH3", "CH2"],  "HARMONIC", [96500.0, 1.54]],
    [["CH2", "C_ar"], "HARMONIC", [96500.0, 1.51]],
    [["C_ar", "C_ar"], "HARMONIC", [96500.0, 1.48]]
  ],
  "Bends" : [
    [["CH3", "CH2", "C_ar"],   "HARMONIC", [62500.0, 114]],
    [["CH2", "C_ar", "CH_ar"], "HARMONIC", [62500.0, 120]],
    [["CH_ar", "C_ar", "C_ar"], "HARMONIC", [62500.0, 120]]
  ],
  "Torsions" : [
    [["CH3", "CH2", "C_ar", "CH_ar"],   "TRAPPE", [0.0, 0.0, 0.0, 0.0]],
    [["CH2", "C_ar", "CH_ar", "CH_ar"], "TRAPPE", [0.0, 0.0, 0.0, 0.0]],
    [["CH_ar", "CH_ar", "C_ar", "C_ar"], "TRAPPE", [0.0, 0.0, 0.0, 0.0]],
    [["CH_ar", "C_ar", "C_ar", "CH_ar"], "TRAPPE", [0.0, 0.0, 200.0, 0.0]]
  ],
  "VanDerWaals" : "auto",
  "Partial-reinsertion" : [
    [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13],
    [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15],
    [0, 1, 2, 3, 4, 5, 6, 7],
    [8, 9, 10, 11, 12, 13, 14, 15]
  ]
}
)";

inline constexpr std::string_view kBadFusedRingJson =
R"({
  "CriticalTemperature" : 553.6,
  "CriticalPressure" : 4073000.0,
  "AcentricFactor" : 0.211,
  "Type" : "flexible",
  "pseudoAtoms" :
    [
      ["CH2_c", [0.0, 0.0, 0.0]],
      ["CH2_c", [1.54, 0.0, 0.0]],
      ["CH2_c", [0.77, 1.33, 0.0]],
      ["CH2_c", [-0.77, 1.33, 0.0]]
    ],
  "Connectivity" : [
    [0, 1],
    [1, 2],
    [2, 0],
    [2, 3],
    [3, 0]
  ],
  "Groups" : [
    { "Type" : "Cycle", "Atoms" : [0, 1, 2, 3] }
  ],
  "Bonds" : [
    [["CH2_c", "CH2_c"], "HARMONIC", [96500.0, 1.54]]
  ],
  "VanDerWaals" : "auto"
}
)";

inline constexpr std::string_view kMethaneJson =
R"({
  "CriticalTemperature" : 190.564,
  "CriticalPressure" : 4599200.0,
  "AcentricFactor" : 0.01142,
  "Type" : "rigid",
  "pseudoAtoms" : 
    [
      ["CH4",[0.0, 0.0, 1.0]]
    ]
} 
)";

inline constexpr std::string_view alkaneJson(std::string_view name)
{
  if (name == "propane") return kPropaneJson;
  if (name == "butane") return kButaneJson;
  if (name == "ethane") return kEthaneJson;
  return {};
}

}  // namespace molecule_fixtures
