
import <exception>;
import <iostream>;
import monte_carlo;
import input_reader;

int main()
{
	try
	{
		InputReader inputReader("pseudo_atoms.def", "force_field_mixing_rules.def",
			                    "force_field.def", "simulation.input");

		switch (inputReader.simulationType)
		{
		case InputReader::SimulationType::MonteCarlo:
		{
			MonteCarlo mc(inputReader);
    		mc.run();
			break;
		}
		default:
			break;
		}
	}
	catch (std::exception const& e)
	{
		std::cerr << e.what();
		exit(-1);
	}
}
