module;

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Weverything"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/chrono.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>

#pragma clang diagnostic pop

export module pybindings;

import <span>;
import <vector>;

import int3;
import double3;
import double3x3;

import atom;
import simulationbox;
import forcefield;

import running_energy;

import component;
import system;

PYBIND11_MODULE(RaspaKit, m)
{
  pybind11::class_<int3>(m, "int3")
        .def(pybind11::init<int32_t, int32_t, int32_t>());

  pybind11::class_<double3>(m, "double3")
        .def(pybind11::init<double, double, double>())
        .def_readwrite("x", &double3::x)
        .def_readwrite("y", &double3::y)
        .def_readwrite("z", &double3::z);

  pybind11::class_<RunningEnergy>(m, "RunningEnergy")
        .def(pybind11::init<>())
        .def_readwrite("moleculeMoleculeVDW", &RunningEnergy::moleculeMoleculeVDW)
        .def_readwrite("frameworkMoleculeVDW", &RunningEnergy::frameworkMoleculeVDW);

  pybind11::class_<Atom>(m, "Atom")
        .def(pybind11::init<double3, double, double, uint16_t, uint8_t, uint32_t>())
        .def_readwrite("position", &Atom::position)
        .def("__repr__", &Atom::repr);

  pybind11::class_<SimulationBox> simulationBox(m, "SimulationBox");
  pybind11::enum_<SimulationBox::Type>(simulationBox, "Type")
    .value("Rectangular", SimulationBox::Type::Rectangular)
    .value("Triclinic", SimulationBox::Type::Triclinic)
    .export_values();
  simulationBox.def(pybind11::init<double, double, double, SimulationBox::Type>(), 
                    pybind11::arg("a"), pybind11::arg("b"), pybind11::arg("c"), 
                    pybind11::arg("type") = SimulationBox::Type::Rectangular)
               .def(pybind11::init<double, double, double, double, double, double, SimulationBox::Type>(),
                    pybind11::arg("a"), pybind11::arg("b"), pybind11::arg("c"), 
                    pybind11::arg("alpha"), pybind11::arg("beta"), pybind11::arg("gamma"), 
                    pybind11::arg("type") = SimulationBox::Type::Rectangular)
               .def(pybind11::init<double3x3,  SimulationBox::Type>())
               .def_readwrite("type", &SimulationBox::type)
               .def_readonly("lengthA", &SimulationBox::lengthA)
               .def_readonly("lengthB", &SimulationBox::lengthB)
               .def_readonly("lengthC", &SimulationBox::lengthC);


  pybind11::class_<VDWParameters>(m, "VDWParameters")
        .def(pybind11::init<double, double>(), pybind11::arg("epsilon"), pybind11::arg("sigma"));

  pybind11::class_<PseudoAtom>(m, "PseudoAtom")
        .def(pybind11::init<std::string, double, double, size_t, bool>(), 
             pybind11::arg("name"), pybind11::arg("mass"), pybind11::arg("charge"), 
             pybind11::arg("atomicNumber"), pybind11::arg("printToPDB"));

  pybind11::class_<ForceField> forceField(m, "ForceField");
  pybind11::enum_<ForceField::MixingRule>(forceField, "MixingRule")
    .value("Lorentz_Berthelot", ForceField::MixingRule::Lorentz_Berthelot)
    .export_values();
  forceField.def(pybind11::init<std::vector<PseudoAtom>, std::vector<VDWParameters>, 
                                ForceField::MixingRule, double, bool, bool>(), 
                 pybind11::arg("pseudoAtoms"), pybind11::arg("parameters"), pybind11::arg("mixingRule"), 
                 pybind11::arg("cutOff"), pybind11::arg("shifted"), pybind11::arg("tailCorrecions"));
  forceField.def("__repr__", &ForceField::repr);


  pybind11::class_<Component>(m, "Component")
        .def(pybind11::init<size_t, std::string, double, SimulationBox, 
                            double, double, double, std::vector<Atom>, size_t, size_t>())
        .def(pybind11::init<size_t, std::string, double, SimulationBox, size_t, 
                            std::vector<Atom>, int3, size_t, size_t>())
        .def_readonly("name", &Component::name)
        .def("__repr__", &Component::repr);

  pybind11::class_<System>(m, "System")
        .def(pybind11::init<size_t, double, double, ForceField, std::vector<Component>, std::vector<size_t>, size_t>())
        .def("precomputeTotalRigidEnergy", &System::precomputeTotalRigidEnergy)
        .def("computeTotalEnergies", &System::computeTotalEnergies)
        .def_readwrite("atomPositions", &System::atomPositions);

}

