#pragma once

#include <string_view>

// Embedded IRMOF fixture data for raspakit unit tests.
namespace irmof_fixtures
{
inline constexpr std::string_view kIrmof1Cif =
R"CIF(data_IRMOF-1

_audit_creation_method RASPA-1.0
_audit_creation_date 2011-3-3
_audit_author_name 'David Dubbeldam'

_citation_author_name        'M. Eddaoudi, J. Kim, N. Rosi, D. Vodak, J. Wachter, M. O. Keeffe, and O.M. Yaghi'
_citation_title              'Systematic design of pore size and functionality in isoreticular MOFs and their application in methane storage'
_citation_journal_abbrev     'Science'
_citation_journal_volume     295
_citation_journal_number     5554
_citation_page_first         469
_citation_page_last          472
_citation_year               2002

_cell_length_a    25.832
_cell_length_b    25.832
_cell_length_c    25.832
_cell_angle_alpha 90
_cell_angle_beta  90
_cell_angle_gamma 90
_cell_volume      17237.5

_symmetry_cell_setting          cubic
_symmetry_space_group_name_Hall '-F 4 2 3'
_symmetry_space_group_name_H-M  'F m -3 m'
_symmetry_Int_Tables_number     225

loop_
_symmetry_equiv_pos_as_xyz
 'x,y,z'
 '-x,-y,z'
 '-x,y,-z'
 'x,-y,-z'
 'z,x,y'
 'z,-x,-y'
 '-z,-x,y'
 '-z,x,-y'
 'y,z,x'
 '-y,z,-x'
 'y,-z,-x'
 '-y,-z,x'
 'y,x,-z'
 '-y,-x,-z'
 'y,-x,z'
 '-y,x,z'
 'x,z,-y'
 '-x,z,y'
 '-x,-z,-y'
 'x,-z,y'
 'z,y,-x'
 'z,-y,x'
 '-z,y,x'
 '-z,-y,-x'
 '-x,-y,-z'
 'x,y,-z'
 'x,-y,z'
 '-x,y,z'
 '-z,-x,-y'
 '-z,x,y'
 'z,x,-y'
 'z,-x,y'
 '-y,-z,-x'
 'y,-z,x'
 '-y,z,x'
 'y,z,-x'
 '-y,-x,z'
 'y,x,z'
 '-y,x,-z'
 'y,-x,-z'
 '-x,-z,y'
 'x,-z,-y'
 'x,z,y'
 '-x,z,-y'
 '-z,-y,x'
 '-z,y,-x'
 'z,-y,-x'
 'z,y,x'
 'x,y+1/2,z+1/2'
 '-x,-y+1/2,z+1/2'
 '-x,y+1/2,-z+1/2'
 'x,-y+1/2,-z+1/2'
 'z,x+1/2,y+1/2'
 'z,-x+1/2,-y+1/2'
 '-z,-x+1/2,y+1/2'
 '-z,x+1/2,-y+1/2'
 'y,z+1/2,x+1/2'
 '-y,z+1/2,-x+1/2'
 'y,-z+1/2,-x+1/2'
 '-y,-z+1/2,x+1/2'
 'y,x+1/2,-z+1/2'
 '-y,-x+1/2,-z+1/2'
 'y,-x+1/2,z+1/2'
 '-y,x+1/2,z+1/2'
 'x,z+1/2,-y+1/2'
 '-x,z+1/2,y+1/2'
 '-x,-z+1/2,-y+1/2'
 'x,-z+1/2,y+1/2'
 'z,y+1/2,-x+1/2'
 'z,-y+1/2,x+1/2'
 '-z,y+1/2,x+1/2'
 '-z,-y+1/2,-x+1/2'
 '-x,-y+1/2,-z+1/2'
 'x,y+1/2,-z+1/2'
 'x,-y+1/2,z+1/2'
 '-x,y+1/2,z+1/2'
 '-z,-x+1/2,-y+1/2'
 '-z,x+1/2,y+1/2'
 'z,x+1/2,-y+1/2'
 'z,-x+1/2,y+1/2'
 '-y,-z+1/2,-x+1/2'
 'y,-z+1/2,x+1/2'
 '-y,z+1/2,x+1/2'
 'y,z+1/2,-x+1/2'
 '-y,-x+1/2,z+1/2'
 'y,x+1/2,z+1/2'
 '-y,x+1/2,-z+1/2'
 'y,-x+1/2,-z+1/2'
 '-x,-z+1/2,y+1/2'
 'x,-z+1/2,-y+1/2'
 'x,z+1/2,y+1/2'
 '-x,z+1/2,-y+1/2'
 '-z,-y+1/2,x+1/2'
 '-z,y+1/2,-x+1/2'
 'z,-y+1/2,-x+1/2'
 'z,y+1/2,x+1/2'
 'x+1/2,y,z+1/2'
 '-x+1/2,-y,z+1/2'
 '-x+1/2,y,-z+1/2'
 'x+1/2,-y,-z+1/2'
 'z+1/2,x,y+1/2'
 'z+1/2,-x,-y+1/2'
 '-z+1/2,-x,y+1/2'
 '-z+1/2,x,-y+1/2'
 'y+1/2,z,x+1/2'
 '-y+1/2,z,-x+1/2'
 'y+1/2,-z,-x+1/2'
 '-y+1/2,-z,x+1/2'
 'y+1/2,x,-z+1/2'
 '-y+1/2,-x,-z+1/2'
 'y+1/2,-x,z+1/2'
 '-y+1/2,x,z+1/2'
 'x+1/2,z,-y+1/2'
 '-x+1/2,z,y+1/2'
 '-x+1/2,-z,-y+1/2'
 'x+1/2,-z,y+1/2'
 'z+1/2,y,-x+1/2'
 'z+1/2,-y,x+1/2'
 '-z+1/2,y,x+1/2'
 '-z+1/2,-y,-x+1/2'
 '-x+1/2,-y,-z+1/2'
 'x+1/2,y,-z+1/2'
 'x+1/2,-y,z+1/2'
 '-x+1/2,y,z+1/2'
 '-z+1/2,-x,-y+1/2'
 '-z+1/2,x,y+1/2'
 'z+1/2,x,-y+1/2'
 'z+1/2,-x,y+1/2'
 '-y+1/2,-z,-x+1/2'
 'y+1/2,-z,x+1/2'
 '-y+1/2,z,x+1/2'
 'y+1/2,z,-x+1/2'
 '-y+1/2,-x,z+1/2'
 'y+1/2,x,z+1/2'
 '-y+1/2,x,-z+1/2'
 'y+1/2,-x,-z+1/2'
 '-x+1/2,-z,y+1/2'
 'x+1/2,-z,-y+1/2'
 'x+1/2,z,y+1/2'
 '-x+1/2,z,-y+1/2'
 '-z+1/2,-y,x+1/2'
 '-z+1/2,y,-x+1/2'
 'z+1/2,-y,-x+1/2'
 'z+1/2,y,x+1/2'
 'x+1/2,y+1/2,z'
 '-x+1/2,-y+1/2,z'
 '-x+1/2,y+1/2,-z'
 'x+1/2,-y+1/2,-z'
 'z+1/2,x+1/2,y'
 'z+1/2,-x+1/2,-y'
 '-z+1/2,-x+1/2,y'
 '-z+1/2,x+1/2,-y'
 'y+1/2,z+1/2,x'
 '-y+1/2,z+1/2,-x'
 'y+1/2,-z+1/2,-x'
 '-y+1/2,-z+1/2,x'
 'y+1/2,x+1/2,-z'
 '-y+1/2,-x+1/2,-z'
 'y+1/2,-x+1/2,z'
 '-y+1/2,x+1/2,z'
 'x+1/2,z+1/2,-y'
 '-x+1/2,z+1/2,y'
 '-x+1/2,-z+1/2,-y'
 'x+1/2,-z+1/2,y'
 'z+1/2,y+1/2,-x'
 'z+1/2,-y+1/2,x'
 '-z+1/2,y+1/2,x'
 '-z+1/2,-y+1/2,-x'
 '-x+1/2,-y+1/2,-z'
 'x+1/2,y+1/2,-z'
 'x+1/2,-y+1/2,z'
 '-x+1/2,y+1/2,z'
 '-z+1/2,-x+1/2,-y'
 '-z+1/2,x+1/2,y'
 'z+1/2,x+1/2,-y'
 'z+1/2,-x+1/2,y'
 '-y+1/2,-z+1/2,-x'
 'y+1/2,-z+1/2,x'
 '-y+1/2,z+1/2,x'
 'y+1/2,z+1/2,-x'
 '-y+1/2,-x+1/2,z'
 'y+1/2,x+1/2,z'
 '-y+1/2,x+1/2,-z'
 'y+1/2,-x+1/2,-z'
 '-x+1/2,-z+1/2,y'
 'x+1/2,-z+1/2,-y'
 'x+1/2,z+1/2,y'
 '-x+1/2,z+1/2,-y'
 '-z+1/2,-y+1/2,x'
 '-z+1/2,y+1/2,-x'
 'z+1/2,-y+1/2,-x'
 'z+1/2,y+1/2,x'

loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_charge
Zn1      Zn     0.2934     0.2066     0.2066     0 
O1       O      0.25       0.25       0.25       0 
O2       O      0.2819     0.2181     0.134      0 
C1       C      0.25       0.25       0.1113     0 
C2       C      0.25       0.25       0.0538     0 
C3       C      0.2829     0.2171     0.0269     0 
H1       H      0.3049     0.1951     0.0448     0 


)CIF";

inline constexpr std::string_view kFlexibleIrmofForceFieldJson =
R"CIF({
  "MixingRule": "Lorentz-Berthelot",
  "TruncationMethod": "shifted",
  "TailCorrections": false,
  "CutOffVDW": 12.0,
  "EwaldPrecision" : 1e-10,
  "PseudoAtoms": [
    {
      "name": "Zn1",
      "framework": true,
      "print_to_output": true,
      "element": "Zn",
      "print_as": "Zn",
      "mass": 65.37,
      "charge": 1.275
    },
    {
      "name": "O1",
      "framework": true,
      "print_to_output": true,
      "element": "O",
      "print_as": "O",
      "mass": 15.9994,
      "charge": -1.5
    },
    {
      "name": "O2",
      "framework": true,
      "print_to_output": true,
      "element": "O",
      "print_as": "O",
      "mass": 15.9994,
      "charge": -0.6
    },
    {
      "name": "C1",
      "framework": true,
      "print_to_output": true,
      "element": "C",
      "print_as": "C",
      "mass": 12.0107,
      "charge": 0.475
    },
    {
      "name": "C2",
      "framework": true,
      "print_to_output": true,
      "element": "C",
      "print_as": "C",
      "mass": 12.0107,
      "charge": 0.125
    },
    {
      "name": "C3",
      "framework": true,
      "print_to_output": true,
      "element": "C",
      "print_as": "C",
      "mass": 12.0107,
      "charge": -0.15
    },
    {
      "name": "H1",
      "framework": true,
      "print_to_output": true,
      "element": "H",
      "print_as": "H",
      "mass": 1.00794,
      "charge": 0.15
    },
    {
      "name": "CH4",
      "framework": false,
      "print_to_output": true,
      "element": "C",
      "print_as": "C",
      "mass": 16.04246,
      "charge": 0.0
    }
  ],
  "SelfInteractions": [
    {
      "name": "Zn1",
      "type": "lennard-jones",
      "parameters": [0.42, 2.7],
      "source": "D. Dubbeldam, K.S. Walton, D.E. Ellis, R.Q. Snurr, Angew. Chem. Int. Ed. 2007, 46, 4496-4499."
    },
    {
      "name": "O1",
      "type": "lennard-jones",
      "parameters": [700.0, 2.98],
      "source": "D. Dubbeldam et al., Angew. Chem. Int. Ed. 2007, 46, 4496-4499."
    },
    {
      "name": "O2",
      "type": "lennard-jones",
      "parameters": [70.5, 3.11],
      "source": "D. Dubbeldam et al., Angew. Chem. Int. Ed. 2007, 46, 4496-4499."
    },
    {
      "name": "C1",
      "type": "lennard-jones",
      "parameters": [47.0, 3.74],
      "source": "D. Dubbeldam et al., Angew. Chem. Int. Ed. 2007, 46, 4496-4499."
    },
    {
      "name": "C2",
      "type": "lennard-jones",
      "parameters": [47.86, 3.47],
      "source": "D. Dubbeldam et al., Angew. Chem. Int. Ed. 2007, 46, 4496-4499."
    },
    {
      "name": "C3",
      "type": "lennard-jones",
      "parameters": [47.86, 3.47],
      "source": "D. Dubbeldam et al., Angew. Chem. Int. Ed. 2007, 46, 4496-4499."
    },
    {
      "name": "H1",
      "type": "lennard-jones",
      "parameters": [7.65, 2.85],
      "source": "D. Dubbeldam et al., Angew. Chem. Int. Ed. 2007, 46, 4496-4499."
    },
    {
      "name": "CH4",
      "type": "lennard-jones",
      "parameters": [158.5, 3.72],
      "source": "M. G. Martin et al., J. Chem. Phys. 2001, 114, 7174-7181."
    }
  ]
}
)CIF";

inline constexpr std::string_view kFlexibleIrmofFrameworkJson =
R"CIF({
  "Type": "Flexible",
  "ExcludeIntra12Interactions": false,
  "ExcludeIntra13Interactions": false,
  "ExcludeIntraBondInteractions": true,
  "ExcludeIntraBendInteractions": true,
  "Intra14VanDerWaalsScalingValue": 1.0,
  "Intra14ChargeChargeScalingValue": 1.0,
  "Bonds": [
    [["C3", "H1"], "HARMONIC", [366001.13136396, 0.95]],
    [["C3", "C3"], "HARMONIC", [483413.91047488, 1.36]],
    [["C2", "C3"], "HARMONIC", [483413.91047488, 1.36]],
    [["C1", "C2"], "HARMONIC", [353750.919316375, 1.42]],
    [["O2", "C1"], "HARMONIC", [543840.64928424, 1.25]]
  ],
  "Bends": [
    [["C1", "C2", "C3"], "HARMONIC", [34926.5543205787, 120.0]],
    [["C2", "C3", "H1"], "HARMONIC", [37263.15559911, 120.0]],
    [["C3", "C3", "H1"], "HARMONIC", [37263.15559911, 120.0]],
    [["C3", "C2", "C3"], "HARMONIC", [90640.10821404, 120.0]],
    [["C3", "C3", "C2"], "HARMONIC", [90640.10821404, 120.0]],
    [["O2", "C1", "O2"], "HARMONIC", [135960.162321060, 130.0]],
    [["O2", "C1", "C2"], "HARMONIC", [54882.4848123699, 115.0]]
  ],
  "Torsions": [
    [["O2", "C1", "C2", "C3"], "TRAPPE", [0.0, 0.0, 1258.890391861, 0.0]],
    [["C1", "C2", "C3", "H1"], "TRAPPE", [0.0, 0.0, 1510.668470234, 0.0]],
    [["C1", "C2", "C3", "C3"], "TRAPPE", [0.0, 0.0, 1510.668470234, 0.0]],
    [["H1", "C3", "C3", "H1"], "TRAPPE", [0.0, 0.0, 1510.668470234, 0.0]],
    [["C2", "C3", "C3", "H1"], "TRAPPE", [0.0, 0.0, 1510.668470234, 0.0]],
    [["C2", "C3", "C3", "C2"], "TRAPPE", [0.0, 0.0, 1510.668470234, 0.0]],
    [["H1", "C3", "C2", "C3"], "TRAPPE", [0.0, 0.0, 1510.668470234, 0.0]],
    [["C3", "C2", "C3", "C3"], "TRAPPE", [0.0, 0.0, 1510.668470234, 0.0]]
  ],
  "ImproperTorsions": [
    [["C2", "C1", "O2", "O2"], "TRAPPE", [0.0, 0.0, 5035.561567446, 0.0]],
    [["C3", "C2", "C3", "C1"], "TRAPPE", [0.0, 0.0, 5035.561567446, 0.0]],
    [["C2", "C3", "C3", "H1"], "TRAPPE", [0.0, 0.0, 186.3157779955, 0.0]]
  ]
}
)CIF";

inline constexpr std::string_view kMinimizationForceFieldJson =
R"CIF({
  "PseudoAtoms": [
    {
      "name": "Zn1", "framework": true, "print_to_output": true,
      "element": "Zn", "print_as": "Zn", "mass": 65.37, "charge": 1.275
    },
    {
      "name": "O1", "framework": true, "print_to_output": true,
      "element": "O", "print_as": "O", "mass": 15.9994, "charge": -1.5
    },
    {
      "name": "O2", "framework": true, "print_to_output": true,
      "element": "O", "print_as": "O", "mass": 15.9994, "charge": -0.6
    },
    {
      "name": "C1", "framework": true, "print_to_output": true,
      "element": "C", "print_as": "C", "mass": 12.0107, "charge": 0.475
    },
    {
      "name": "C2", "framework": true, "print_to_output": true,
      "element": "C", "print_as": "C", "mass": 12.0107, "charge": 0.125
    },
    {
      "name": "C3", "framework": true, "print_to_output": true,
      "element": "C", "print_as": "C", "mass": 12.0107, "charge": -0.15
    },
    {
      "name": "H1", "framework": true, "print_to_output": true,
      "element": "H", "print_as": "H", "mass": 1.00794, "charge": 0.15
    },
    {
      "name": "C_co2", "framework": false, "print_to_output": true,
      "element": "C", "print_as": "C", "mass": 12.0, "charge": 0.6512
    },
    {
      "name": "O_co2", "framework": false, "print_to_output": true,
      "element": "O", "print_as": "O", "mass": 15.9994, "charge": -0.3256
    }
  ],
  "SelfInteractions": [
    {
      "name": "Zn1", "type": "lennard-jones", "parameters": [0.42, 2.7],
      "source": "D. Dubbeldam et al., Angew. Chem. Int. Ed. 2007, 46, 4496-4499."
    },
    {
      "name": "O1", "type": "lennard-jones", "parameters": [700.0, 2.98],
      "source": "D. Dubbeldam et al., Angew. Chem. Int. Ed. 2007, 46, 4496-4499."
    },
    {
      "name": "O2", "type": "lennard-jones", "parameters": [70.5, 3.11],
      "source": "D. Dubbeldam et al., Angew. Chem. Int. Ed. 2007, 46, 4496-4499."
    },
    {
      "name": "C1", "type": "lennard-jones", "parameters": [47.0, 3.74],
      "source": "D. Dubbeldam et al., Angew. Chem. Int. Ed. 2007, 46, 4496-4499."
    },
    {
      "name": "C2", "type": "lennard-jones", "parameters": [47.86, 3.47],
      "source": "D. Dubbeldam et al., Angew. Chem. Int. Ed. 2007, 46, 4496-4499."
    },
    {
      "name": "C3", "type": "lennard-jones", "parameters": [47.86, 3.47],
      "source": "D. Dubbeldam et al., Angew. Chem. Int. Ed. 2007, 46, 4496-4499."
    },
    {
      "name": "H1", "type": "lennard-jones", "parameters": [7.65, 2.85],
      "source": "D. Dubbeldam et al., Angew. Chem. Int. Ed. 2007, 46, 4496-4499."
    },
    {
      "name": "O_co2", "type": "lennard-jones", "parameters": [85.671, 3.017],
      "source": "A. Garcia-Sanchez et al., J. Phys. Chem. C 2009, 113, 8814-8820."
    },
    {
      "name": "C_co2", "type": "lennard-jones", "parameters": [29.933, 2.745],
      "source": "A. Garcia-Sanchez et al., J. Phys. Chem. C 2009, 113, 8814-8820."
    }
  ],
  "MixingRule": "Lorentz-Berthelot",
  "TruncationMethod": "shifted",
  "TailCorrections": false,
  "CutOffVDW": 11.8
}
)CIF";

inline constexpr std::string_view kCO2Json =
R"CIF({
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
)CIF";

}  // namespace irmof_fixtures
