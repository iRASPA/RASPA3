module;

module energy_status;

import units;
import component;
import energy_status_intra;
import energy_status_inter;

import <string>;
import <iostream>;
import <sstream>;
import <vector>;
import print;

std::string EnergyStatus::printEnergyStatus(const std::vector<Component>& components, const std::string &label)
{
	std::ostringstream stream;

    double conv = Units::EnergyToKelvin;
    std::print(stream, "Energy status {}\n", label);
    std::print(stream, "===============================================================================\n\n");
    std::print(stream, "Total potential energy : {: .6e}\n", conv *( intraEnergy.total() + interEnergy.total()));
    std::print(stream, "    Van der Waals:        {: .6e}\n", conv * interEnergy.VanDerWaals);
    std::print(stream, "    Van der Waals (Tail): {: .6e}\n", conv * interEnergy.VanDerWaalsTailCorrection);
    std::print(stream, "    Coulombic Real:       {: .6e}\n", conv * interEnergy.CoulombicReal);
    std::print(stream, "    Coulombic Fourier:    {: .6e}\n", conv * interEnergy.CoulombicFourier);
    std::print(stream, "    dU/dlambda:           {: .6e}\n\n", conv * dUdlambda);
	
	for (size_t i = 0; i < components.size(); ++i)
	{
        std::print(stream, "    Component: {} [{}]\n", i, components[i].name);
        std::print(stream, "    ---------------------------------------------------------------------------\n\n");
        std::print(stream, "    Molecule bond:             {: .6e}\n", conv * (*this)(i).bond);
        std::print(stream, "    Molecule bend:             {: .6e}\n", conv * (*this)(i).bend);
        std::print(stream, "    Molecule inversionBend     {: .6e}\n", conv * (*this)(i).inversionBend);
        std::print(stream, "    Molecule ureyBradley       {: .6e}\n", conv * (*this)(i).ureyBradley);
        std::print(stream, "    Molecule torsion           {: .6e}\n", conv * (*this)(i).torsion);
        std::print(stream, "    Molecule improperTorsion   {: .6e}\n", conv * (*this)(i).improperTorsion);
        std::print(stream, "    Molecule bondBond          {: .6e}\n", conv * (*this)(i).bondBond);
        std::print(stream, "    Molecule bondBend          {: .6e}\n", conv * (*this)(i).bondBend);
        std::print(stream, "    Molecule bondTorsion       {: .6e}\n", conv * (*this)(i).bondTorsion);
        std::print(stream, "    Molecule bendBend          {: .6e}\n", conv * (*this)(i).bendBend);
        std::print(stream, "    Molecule bendTorsion       {: .6e}\n", conv * (*this)(i).bendTorsion);
        std::print(stream, "    Molecule intraVDW          {: .6e}\n", conv * (*this)(i).intraVDW);
        std::print(stream, "    Molecule intraChargeCharge {: .6e}\n\n", conv * (*this)(i).intraChargeCharge);
		for (size_t j = 0; j < components.size(); ++j)
		{
			std::print(stream, "    Component: {}-{} [{}-{}]\n", i, j, components[i].name, components[j].name);
            std::print(stream, "        Van der Waals:        {: .6e}\n", conv * (*this)(i,j).VanDerWaals);
            std::print(stream, "        Van der Waals (Tail): {: .6e}\n", conv * (*this)(i,j).VanDerWaalsTailCorrection);
            std::print(stream, "        Coulombic Real:       {: .6e}\n", conv * (*this)(i,j).CoulombicReal);
            std::print(stream, "        Coulombic Fourier:    {: .6e}\n", conv * (*this)(i,j).CoulombicFourier);
            std::print(stream, "        -----------------------------------------------------------------------\n");
            std::print(stream, "        Sum                   {: .6e}\n\n", conv * (*this)(i,j).totalInter);
		}
	}
    std::print(stream, "\n");
	return stream.str();
}

