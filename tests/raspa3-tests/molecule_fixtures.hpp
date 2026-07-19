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
  "RigidBodies" : [
    [1, 2, 3]
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
  "RigidBodies" : [
    [2, 3, 4, 5, 6, 7],
    [8, 9, 10, 11, 12, 13]
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

// United-atom trans-decalin: two fused six-membered rings sharing the 0-1 bridgehead bond. The
// reference geometry lies on a diamond lattice (all bonds 1.54 Angstrom, all bends 109.47 degrees),
// so every ring starts in a perfect chair. One cyclic cluster of all ten beads with two closure
// bonds; rings are inferred from the connectivity, no declaration needed.
inline constexpr std::string_view kDecalinJson =
R"({
  "CriticalTemperature" : 687.0,
  "CriticalPressure" : 3200000.0,
  "AcentricFactor" : 0.286,
  "pseudoAtoms" :
    [
      ["CH2_c", [0.978020, -0.088910, -0.444560]],
      ["CH2_c", [0.088920, -0.978010, 0.444540]],
      ["CH2_c", [-0.800180, -1.867110, -0.444560]],
      ["CH2_c", [-1.689380, -0.978010, -1.333660]],
      ["CH2_c", [-0.800180, -0.088910, -2.222760]],
      ["CH2_c", [0.088920, 0.800190, -1.333660]],
      ["CH2_c", [-0.800180, -0.088910, 1.333640]],
      ["CH2_c", [0.088920, 0.800190, 2.222840]],
      ["CH2_c", [0.978020, 1.689290, 1.333640]],
      ["CH2_c", [1.867120, 0.800190, 0.444540]]
    ],
  "Connectivity" : [
    [0, 1],
    [0, 5],
    [0, 9],
    [1, 2],
    [1, 6],
    [2, 3],
    [3, 4],
    [4, 5],
    [6, 7],
    [7, 8],
    [8, 9]
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

// A fused pair of three-membered rings (bicyclobutane-like connectivity): two cycles sharing the
// 0-2 edge. Under the fragment-graph model this is a valid molecule (one cyclic cluster with two
// closure bonds) instead of a parse error.
inline constexpr std::string_view kFusedRingJson =
R"({
  "CriticalTemperature" : 553.6,
  "CriticalPressure" : 4073000.0,
  "AcentricFactor" : 0.211,
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
  "Bonds" : [
    [["CH2_c", "CH2_c"], "HARMONIC", [96500.0, 1.54]]
  ],
  "Bends" : [
    [["CH2_c", "CH2_c", "CH2_c"], "HARMONIC", [62500.0, 60.0]]
  ],
  "Torsions" : [
    [["CH2_c", "CH2_c", "CH2_c", "CH2_c"], "TRAPPE", [0.0, 355.03, -68.19, 791.32]]
  ],
  "VanDerWaals" : "auto"
}
)";

// United-atom norbornane (bicyclo[2.2.1]heptane): two bridgehead beads (0 and 3) joined by three
// bridges of two, two, and one bead. A bridged bicyclic is one cyclic cluster of all seven beads
// with two closure bonds; like fused rings it needs no declaration.
inline constexpr std::string_view kNorbornaneJson =
R"({
  "CriticalTemperature" : 638.0,
  "CriticalPressure" : 3900000.0,
  "AcentricFactor" : 0.195,
  "pseudoAtoms" :
    [
      ["CH2_c", [1.116, 0.000, 0.325]],
      ["CH2_c", [0.749, 1.250, -0.470]],
      ["CH2_c", [-0.749, 1.250, -0.470]],
      ["CH2_c", [-1.116, 0.000, 0.325]],
      ["CH2_c", [-0.749, -1.250, -0.470]],
      ["CH2_c", [0.749, -1.250, -0.470]],
      ["CH2_c", [0.000, 0.000, 1.253]]
    ],
  "Connectivity" : [
    [0, 1],
    [1, 2],
    [2, 3],
    [3, 4],
    [4, 5],
    [5, 0],
    [0, 6],
    [3, 6]
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

// United-atom methylcyclohexane: a six-membered flexible ring with one CH3 tail on bead 0. Grown as
// a ring-closure step for the ring followed by an ordinary flexible attachment for the tail.
inline constexpr std::string_view kMethylcyclohexaneJson =
R"({
  "CriticalTemperature" : 572.2,
  "CriticalPressure" : 3470000.0,
  "AcentricFactor" : 0.235,
  "pseudoAtoms" :
    [
      ["CH_c",  [1.465493, 0.000000, 0.236606]],
      ["CH2_c", [0.732747, 1.269154, -0.236606]],
      ["CH2_c", [-0.732747, 1.269154, 0.236606]],
      ["CH2_c", [-1.465493, 0.000000, -0.236606]],
      ["CH2_c", [-0.732747, -1.269154, 0.236606]],
      ["CH2_c", [0.732747, -1.269154, -0.236606]],
      ["CH3",   [2.905493, 0.000000, 0.786606]]
    ],
  "Connectivity" : [
    [0, 1],
    [1, 2],
    [2, 3],
    [3, 4],
    [4, 5],
    [5, 0],
    [0, 6]
  ],
  "Bonds" : [
    [["CH2_c", "CH2_c"], "HARMONIC", [96500.0, 1.54]],
    [["CH_c", "CH2_c"],  "HARMONIC", [96500.0, 1.54]],
    [["CH_c", "CH3"],    "HARMONIC", [96500.0, 1.54]]
  ],
  "Bends" : [
    [["CH2_c", "CH2_c", "CH2_c"], "HARMONIC", [62500.0, 114.0]],
    [["CH2_c", "CH_c", "CH2_c"],  "HARMONIC", [62500.0, 114.0]],
    [["CH_c", "CH2_c", "CH2_c"],  "HARMONIC", [62500.0, 114.0]],
    [["CH3", "CH_c", "CH2_c"],    "HARMONIC", [62500.0, 114.0]]
  ],
  "Torsions" : [
    [["CH2_c", "CH2_c", "CH2_c", "CH2_c"], "TRAPPE", [0.0, 355.03, -68.19, 791.32]],
    [["CH_c", "CH2_c", "CH2_c", "CH2_c"],  "TRAPPE", [0.0, 355.03, -68.19, 791.32]],
    [["CH2_c", "CH_c", "CH2_c", "CH2_c"],  "TRAPPE", [0.0, 355.03, -68.19, 791.32]],
    [["CH3", "CH_c", "CH2_c", "CH2_c"],    "TRAPPE", [0.0, 355.03, -68.19, 791.32]]
  ],
  "VanDerWaals" : "auto"
}
)";

inline constexpr std::string_view kMethaneJson =
R"({
  "CriticalTemperature" : 190.564,
  "CriticalPressure" : 4599200.0,
  "AcentricFactor" : 0.01142,
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
