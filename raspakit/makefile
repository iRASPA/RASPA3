CPP    := $(wildcard *.cpp)
PCM    := $(CPP:.cpp=.pcm)
OBJ    := $(CPP:.cpp=.o)
OBJPCM := $(CPP:.cpp=.pcm.o)

.PHONY: default
default: libraspakit;

.SUFFIXES:
.SUFFIXES: .pcm .ixx *.cpp

.SECONDARY:

$(OBJ): %.o: %.cpp %.pcm
	$(CXX) $(CXXFLAGS) $< -c -o $@


%.pcm: %.ixx
	$(CXX) $(CXXFLAGS) -x c++-module --precompile -o $@ $<
	$(CXX) $(CXXFLAGS) -x c++-module $< -c -o $@.o

mc_moves_gibbs_swap_cfcmc.pcm mc_moves_gibbs_swap_cbmc.pcm mc_moves_insertion.pcm mc_moves_gibbs_volume.pcm mc_moves_reaction.pcm mc_moves_reaction_cfcmc_cbmc.pcm mc_moves_swap_cfcmc.pcm mc_moves_translation.pcm mc_moves_random_rotation.pcm mc_moves_swap_cfcmc_cbmc.pcm mc_moves_reaction_cfcmc.pcm mc_moves_deletion.pcm mc_moves_widom.pcm mc_moves_identity_change.pcm mc_moves_random_translation.pcm mc_moves_volume.pcm mc_moves_reinsertion.pcm mc_moves_rotation.pcm:
	: do nothing for $@
system_cbmc_rigid_deletion.pcm system_cbmc_rigid_insertion.pcm system_cbmc_flexible_insertion.pcm system_cbmc.pcm system_hardware_info.pcm system_compute_fluid_properties.pcm system_interactions_ewald.pcm  system_interactions_framework_molecule_pressure.pcm system_integration.pcm system_interactions_framework_molecule.pcm system_interactions_ewald_pressure.pcm system_md_loop.pcm system_interactions_tailcorrection.pcm system_interactions_intermolecular_pressure.pcm system_interactions.pcm system_cbmc_rigid_reinsertion.pcm system_reservoir.pcm system_interactions_intermolecular.pcm system_print_movestatistics.pcm:
	: do nothing for $@

# dependencies
reaction.pcm: property_lambda_probability_histogram.pcm
reactions.pcm: reaction.pcm
atom.pcm: scaling.pcm
running_energy.pcm: atom.pcm
multi_site_isotherm.pcm: isotherm.pcm
simulationbox.pcm: units.pcm averages.pcm
forcefield.pcm: units.pcm
property_lambda_probability_histogram.pcm: simulationbox.pcm averages.pcm
property_widom.pcm: averages.pcm
mc_moves_probabilities_particles.pcm: move_statistics.pcm
mc_moves_probabilities_system.pcm: move_statistics.pcm
component.pcm: property_lambda_probability_histogram.pcm mc_moves_cputime.pcm mc_moves_count.pcm mc_moves_probabilities_particles.pcm move_statistics.pcm cif_reader.pcm bond_potential.pcm multi_site_isotherm.pcm atom.pcm averages.pcm property_widom.pcm forcefield.pcm
energy_status_inter.pcm: energy_factor.pcm
energy_status_intra.pcm: energy_factor.pcm
energy_status.pcm: component.pcm energy_status_inter.pcm energy_status_intra.pcm
property_energy.pcm: energy_status.pcm energy_status_inter.pcm energy_status_intra.pcm averages.pcm
loadings.pcm: simulationbox.pcm component.pcm
property_loading.pcm: loadings.pcm averages.pcm
enthalpy_of_adsorption.pcm: units.pcm averages.pcm energy_status.pcm
property_enthalpy.pcm: enthalpy_of_adsorption.pcm averages.pcm loadings.pcm
property_pressure.pcm: averages.pcm
cif_reader.pcm: forcefield.pcm simulationbox.pcm atom.pcm
potential_energy_vdw.pcm: energy_factor.pcm forcefield.pcm
potential_gradient_vdw.pcm: force_factor.pcm forcefield.pcm
potential_correction_vdw.pcm: forcefield.pcm
potential_energy_coulomb.pcm: energy_factor.pcm forcefield.pcm
potential_gradient_coulomb.pcm: force_factor.pcm forcefield.pcm
cbmc_growing_status.pcm: atom.pcm component.pcm running_energy.pcm bond_potential.pcm
cbmc.pcm: atom.pcm running_energy.pcm energy_status.pcm energy_status_intra.pcm energy_status_inter.pcm
sample_movies.pcm: units.pcm atom.pcm simulationbox.pcm forcefield.pcm
property_simulationbox.pcm: averages.pcm simulationbox.pcm
system.pcm: transition_matrix.pcm reaction.pcm reactions.pcm mc_moves_cputime.pcm mc_moves_count.pcm mc_moves_probabilities_system.pcm move_statistics.pcm atom.pcm component.pcm simulationbox.pcm forcefield.pcm units.pcm loadings.pcm averages.pcm sample_movies.pcm enthalpy_of_adsorption.pcm energy_status.pcm energy_status_inter.pcm energy_status_intra.pcm cbmc.pcm property_rdf.pcm property_simulationbox.pcm property_energy.pcm property_pressure.pcm property_loading.pcm property_enthalpy.pcm property_widom.pcm energy_factor.pcm force_factor.pcm running_energy.pcm isotherm.pcm multi_site_isotherm.pcm pressure_range.pcm bond_potential.pcm
system_hardware_info.cpp: system.pcm
system_compute_fluid_properties.cpp: system.pcm
system_print_movestatistics.cpp: system.pcm 
system_cbmc.cpp: system.pcm 
system_cbmc_rigid_deletion.cpp: system.pcm 
system_cbmc_rigid_insertion.cpp: system.pcm 
system_cbmc_rigid_reinsertion.cpp: system.pcm 
system_cbmc_flexible_insertion.cpp: system.pcm 
system_interactions.cpp: system.pcm
system_interactions_framework_molecule.cpp: system.pcm
system_interactions_framework_molecule_pressure.cpp: system.pcm
system_interactions_intermolecular.cpp: system.pcm
system_interactions_intermolecular_pressure.cpp: system.pcm
system_interactions_tailcorrection.cpp: system.pcm
system_interactions_ewald.cpp: system.pcm
system_interactions_ewald_pressure.cpp: system.pcm
system_md_loop.cpp: system.pcm
system_integration.cpp: system.pcm
system_reservoir.cpp: system.pcm
mc_moves.pcm: transition_matrix.pcm system.pcm cbmc.pcm component.pcm atom.pcm simulationbox.pcm energy_status.pcm running_energy.pcm
mc_moves_translation.cpp: mc_moves.pcm
mc_moves_random_translation.cpp: mc_moves.pcm
mc_moves_rotation.cpp: mc_moves.pcm
mc_moves_random_rotation.cpp: mc_moves.pcm
mc_moves_reinsertion.cpp: mc_moves.pcm
mc_moves_insertion.cpp: mc_moves.pcm
mc_moves_deletion.cpp: mc_moves.pcm
mc_moves_volume.cpp: mc_moves.pcm
mc_moves_gibbs_volume.cpp: mc_moves.pcm
mc_moves_identity_change.cpp: mc_moves.pcm
mc_moves_swap_cfcmc.cpp: mc_moves.pcm
mc_moves_swap_cfcmc_cbmc.cpp: mc_moves.pcm
mc_moves_gibbs_swap_cbmc.cpp: mc_moves.pcm
mc_moves_gibbs_swap_cfcmc.cpp: mc_moves.pcm
mc_moves_reaction.cpp: mc_moves.pcm
mc_moves_reaction_cfcmc.cpp: mc_moves.pcm
mc_moves_reaction_cfcmc_cbmc.cpp: mc_moves.pcm
mc_moves_widom.cpp: mc_moves.pcm
monte_carlo.pcm: averages.pcm system.pcm property_enthalpy.pcm mc_moves.pcm input_reader.pcm energy_status.pcm
monte_carlo_transition_matrix.pcm: averages.pcm system.pcm property_enthalpy.pcm mc_moves.pcm input_reader.pcm energy_status.pcm
input_reader.pcm: transition_matrix.pcm system.pcm atom.pcm component.pcm simulationbox.pcm forcefield.pcm loadings.pcm enthalpy_of_adsorption.pcm energy_status.pcm averages.pcm
isotherm_fitting.pcm: multi_site_isotherm.pcm component.pcm input_reader.pcm system.pcm
isotherm_fitting_simulation.pcm: isotherm_fitting.pcm input_reader.pcm component.pcm system.pcm multi_site_isotherm.pcm
mixture_prediction.pcm: atom.pcm isotherm.pcm multi_site_isotherm.pcm component.pcm bond_potential.pcm system.pcm
mixture_prediction_simulation.pcm: input_reader.pcm component.pcm system.pcm multi_site_isotherm.pcm mixture_prediction.pcm
breakthrough.pcm: component.pcm input_reader.pcm mixture_prediction.pcm system.pcm
breakthrough_simulation.pcm: breakthrough.pcm input_reader.pcm mixture_prediction.pcm system.pcm simulationbox.pcm component.pcm
pybindings.pcm: system.pcm component.pcm atom.pcm running_energy.pcm

libraspakit: $(OBJ)
	ar -rs libraspakit.a *.o
	ranlib -c libraspakit.a

clean:
	rm -f *.pcm *.o *.a
