/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2018, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "custom.h"
#include "../BioFVM/BioFVM.h"  
using namespace BioFVM;

// declare cell definitions here 

std::vector<bool> nodes;

void create_cell_types( void )
{
	// set the random seed 
	if (parameters.ints.find_index("random_seed") != -1)
	{
		SeedRandom(parameters.ints("random_seed"));
	}
	
	/* 
	   Put any modifications to default cell definition here if you 
	   want to have "inherited" by other cell types. 
	   
	   This is a good place to set default functions. 
	*/ 

	initialize_default_cell_definition(); 
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 

	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;
	cell_defaults.functions.update_phenotype = NULL; 
	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.pre_update_intracellular = pre_update_intracellular; 
	cell_defaults.functions.post_update_intracellular = post_update_intracellular; 
	cell_defaults.functions.custom_cell_rule = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 
	
	cell_defaults.custom_data.add_variable(parameters.strings("node_to_visualize"), "dimensionless", 0.0 ); //for paraview visualization

	/*
	   This parses the cell definitions in the XML config file. 
	*/
	
	initialize_cell_definitions_from_pugixml(); 
	
	/* 
	   Put any modifications to individual cell definitions here. 
	   
	   This is a good place to set custom functions. 
	*/ 
	
	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/

	Cell_Definition* pFounderCD2 = find_cell_definition("Synthetic_SMC");
	if (pFounderCD2 != NULL) { pFounderCD2->functions.update_velocity = founder_barrier_velocity; }

	Cell_Definition* pSMC_CD = find_cell_definition("Smooth_Muscle_Cell");
	if (pSMC_CD != NULL) { pSMC_CD->functions.update_velocity = contraction_velocity; }

	build_cell_definitions_maps(); 

	/*
	   This intializes cell signal and response dictionaries 
	*/

	setup_signal_behavior_dictionaries();

	/*
	   This summarizes the setup. 
	*/
	
	display_cell_definitions( std::cout ); 


	return; 
}

void setup_microenvironment( void )
{
	// set domain parameters 
	
	// put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here. 
	
	// initialize BioFVM 
	
	initialize_microenvironment(); 	
	
	return; 
}

void setup_tissue( void )
{
	// load cells from your CSV file
	load_cells_from_pugixml(); 	
}

void apply_contraction_pull( Cell* pCell, double dt )
{
	// Reversible functional vasoconstriction: while the Boolean node
	// CONTR_1 (myosin phosphorylation / MLCK-MLCP balance, requires
	// Myocd -- i.e. only contractile-phenotype SMCs can do this) is ON,
	// add a small inward radial velocity pulling the agent toward the
	// vessel centerline. This is a VELOCITY term, not a position clamp:
	// when CONTR_1 turns back OFF (e.g. amlodipine blocking Ltcc), no
	// inward term is added and the existing repulsion/pressure
	// mechanics relax the ring back outward on their own -- genuine
	// vasorelaxation, distinct from the irreversible structural
	// narrowing driven by Synthetic_SMC accumulation.
	if ( pCell->phenotype.intracellular == NULL ) { return; }
	static double contraction_speed = parameters.doubles("contraction_speed");
	if ( !pCell->phenotype.intracellular->get_boolean_variable_value("CONTR_1") ) { return; }
	double r = sqrt( pCell->position[0]*pCell->position[0] + pCell->position[1]*pCell->position[1] );
	if (r > 1e-6) {
		double nx = pCell->position[0]/r, ny = pCell->position[1]/r;
		pCell->velocity[0] -= contraction_speed*nx;
		pCell->velocity[1] -= contraction_speed*ny;
	}
}

void apply_lumen_barrier( Cell* pCell, double dt )
{
	// Basement membrane / internal elastic lamina as an impenetrable
	// barrier: no SMC-lineage agent (contractile OR synthetic/founder)
	// may ever occupy the flowing-blood lumen, whether pushed there by
	// neointimal migration or pulled there by its own CONTR_1-driven
	// contraction -- crossing an intact endothelium/IEL in vivo needs
	// active proteolytic degradation, not simple crowding or myosin
	// contraction. The Endothelial_Cell ring itself is NOT subject to
	// this clamp -- it is free to be mechanically compressed inward
	// (below the nominal lumen_radius) by pressure from the SMC-lineage
	// mass outside it, which is exactly how structural/functional
	// lumen narrowing is read out (5th-percentile EC radial position).
	// Implemented as a VELOCITY correction, never a position teleport --
	// an earlier version called Cell::assign_position() every mechanics
	// step to clamp position directly, which calls register_agent()
	// (just a list push_back) without removing the cell from its
	// previous voxel, silently duplicating its neighbor-list
	// registration on every call and causing a runaway force/position
	// blowup within 1-2 simulated days. Cancelling the inward velocity
	// component here instead lets position/voxel bookkeeping go through
	// PhysiCell's own (correct) update_position() pathway, same as
	// every other cell.
	static double lumen_r = parameters.doubles("lumen_radius");
	double r = sqrt( pCell->position[0]*pCell->position[0] + pCell->position[1]*pCell->position[1] );
	if (r < lumen_r && r > 1e-6) {
		double nx = pCell->position[0]/r, ny = pCell->position[1]/r;
		double radial_v = pCell->velocity[0]*nx + pCell->velocity[1]*ny;
		if (radial_v < 0) { pCell->velocity[0] -= radial_v*nx; pCell->velocity[1] -= radial_v*ny; }
		double penetration = lumen_r - r;
		double push = penetration * 0.5;
		pCell->velocity[0] += push*nx; pCell->velocity[1] += push*ny;
	}
}

void contraction_velocity( Cell* pCell, Phenotype& phenotype, double dt )
{
	// Wired to Smooth_Muscle_Cell. Same three-stage composition as
	// founder_barrier_velocity below: standard mechanics, then the
	// CONTR_1-driven inward pull, then the impenetrable-lumen clamp
	// (which always has the final say -- contraction cannot push a
	// contractile SMC through the endothelium any more than migration
	// can push a founder cell through it).
	standard_update_cell_velocity( pCell, phenotype, dt );
	apply_contraction_pull( pCell, dt );
	apply_lumen_barrier( pCell, dt );
}

void founder_barrier_velocity( Cell* pCell, Phenotype& phenotype, double dt )
{
	// Normal mechanics first (adhesion/repulsion/EC contact, etc.) --
	// this is the SAME function cell_defaults.functions.update_velocity
	// already uses, just called explicitly so we can add a hard
	// constraint on top for this one cell type.
	standard_update_cell_velocity( pCell, phenotype, dt );

	// Reversible functional vasoconstriction (CONTR_1), same as
	// Smooth_Muscle_Cell's contraction_velocity -- applied BEFORE the
	// lumen clamp below, so the clamp still has the final say if
	// contraction would otherwise pull this cell past the boundary.
	apply_contraction_pull( pCell, dt );
	apply_lumen_barrier( pCell, dt );
}

void update_amlodipine_pk( double current_time, double dt )
{
	static bool periodic_dosing = parameters.bools("amlo_periodic_dosing");
	if ( !periodic_dosing ) { return; }

	static double ka = parameters.doubles("amlo_pk_ka");
	static double ke = parameters.doubles("amlo_pk_ke");
	static double dose_interval = parameters.doubles("amlo_dose_interval_min");
	static double dose_amount = parameters.doubles("amlo_pk_dose_amount");
	static double ss_peak = parameters.doubles("amlo_pk_ss_peak");
	static double base_rate = parameters.doubles("amlo_pk_base_secretion_rate");

	static double A_gut = 0.0;
	static double C_level = 0.0;
	static double next_dose_time = 0.0;
	static std::vector<Cell*> ec_cells;
	static bool ec_cells_found = false;
	static int amlo_index = microenvironment.find_density_index("AMLO");

	if ( !ec_cells_found )
	{
		for ( auto pC : *all_cells ) { if ( pC->type_name == "Endothelial_Cell" ) { ec_cells.push_back(pC); } }
		ec_cells_found = true;
	}

	if ( current_time >= next_dose_time - 1e-6 )
	{
		A_gut += dose_amount;
		next_dose_time += dose_interval;
	}

	double dA = -ka * A_gut;
	double dC = ka * A_gut - ke * C_level;
	A_gut += dA * dt;
	C_level += dC * dt;
	if ( A_gut < 0 ) { A_gut = 0; }
	if ( C_level < 0 ) { C_level = 0; }

	double secretion_rate = base_rate * (C_level / ss_peak);
	for ( auto pC : ec_cells ) { pC->phenotype.secretion.secretion_rates[amlo_index] = secretion_rate; }
}

void pre_update_intracellular( Cell* pCell, Phenotype& phenotype, double dt )
{
	// no-op: sample's $time_scale desync logic doesn't apply to this model
}

void post_update_intracellular( Cell* pCell, Phenotype& phenotype, double dt )
{
	color_node(pCell);
	if ( pCell->phenotype.intracellular == NULL ) { return; }
	static double p_half = parameters.doubles("pressure_half_max");
	static double p_hill = parameters.doubles("pressure_hill_power");
	double pressure = pCell->state.simple_pressure;
	pCell->custom_data["pressure_debug"] = pressure;
	double ratio_n = pow( pressure / p_half, p_hill );
	double gate = 1.0 / (1.0 + ratio_n); // Hill-type suppression: 1 at zero pressure, 0.5 at p_half, ->0 when crowded
	static double n_half = parameters.doubles("population_half_max");
	static double n_hill = parameters.doubles("population_hill_power");
	double n_total = (double) (*all_cells).size();
	double n_ratio_n = pow( n_total / n_half, n_hill );
	double n_gate = 1.0 / (1.0 + n_ratio_n);
	double current_rate = PhysiCell::get_single_behavior(pCell, "cycle entry");
	PhysiCell::set_single_behavior(pCell, "cycle entry", current_rate * gate * n_gate);
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{
	std::vector< std::string > output( 4 , "rgb(0,0,0)" );
	
	if ( pCell->phenotype.intracellular == NULL ) { return output; }
	if ( !pCell->phenotype.intracellular->get_boolean_variable_value( parameters.strings("node_to_visualize") ) )
	{
		output[0] = "rgb(255,0,0)";
		output[2] = "rgb(125,0,0)";
		
	}
	else{
		output[0] = "rgb(0, 255,0)";
		output[2] = "rgb(0, 125,0)";
	}
	
	return output;
}

void color_node(Cell* pCell){
	if ( pCell->phenotype.intracellular == NULL ) { return; }
	std::string node_name = parameters.strings("node_to_visualize");
	pCell->custom_data[node_name] = pCell->phenotype.intracellular->get_boolean_variable_value(node_name);
}