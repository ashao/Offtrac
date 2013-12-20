/*
 * ********+*********+*********+*********+*********+*********+*********+*******
 * Offtrac Version 2.0
 * By Andrew Shao
 * 9 April 2013
 * MODEL DESCRIPTION:
 * 	This code performs advection/diffusion of passive tracers based on the
 * 	3D flow fields, temperature, salinity, etc. fields derived from the
 * 	Geophysical Fluid Dynamics Laboratory's Hallberg Isopycnal Model
 * 	(deprecated as of 2012) and the General Ocean Layer Dynamics Model.
 *
 * VERSION DESCRIPTION:
 * 	This version represents a major rewrite of the Offtrac code as developed by
 * 	LuAnne Thompson, David Darr, Curtis Deutsch, Holger Brix, David Trossman,
 * 	and Andrew Shao.
 *
 * 	Most of the arrays used now are statically allocated at
 * 	runtime. This was done to more efficiently allocate memory at the expense
 * 	of flexibility. This was justified because the code compiles quickly and
 * 	it simplifies the code.
 *
 * 	Also, structures were used as much as possible to improve readability and
 * 	organization.
 *
 * 	IO is primarily done via netCDF.
 *
 * 	By default, Offtrac runs using 12-month climatologies though the capability
 * 	for using hindcasts and/or longer climatologies (not recommended due to
 * 	issues with wrap-around) is included.
 *
 * ********+*********+*********+*********+*********+*********+*********+*******
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vartypes.h"
#include "initialize_offtrac.h"
#include "user_options.h"

// Declare all the datatypes in vartypes.h
ocean_geometry_t ocean_geometry;
ocean_3D_t ocean_3D;
ocean_2D_t ocean_2D;
io_t io;
iterinfo_t iterinfo;
filenames_t filenames;

int main ( void )
{

	strcpy(filenames.metrics,"ocean_geometry.nc");
	strcpy(filenames.namefile,"input.nml");
	// Read in and parse a namelist file with runtime configurable options
	read_parse_runtime;

	// Initialize arrays and variables used in Offtrac
	initialize_ocean_geometry;

	// Set initial condition for tracer(s)

	/*
	set_initial_condition;

	for(iterno=0;iterno<itermax;iterno++)
	{

		// Update time information
		update_time_information;

		// Read in physical fields from GOLD/HIM
		read_physical_fields;

		// Perform advection/diffusion for one timestep
		tracer_advect_diffuse;

		// Apply tracer sources/sinks
		tracer_sources_sinks;

		// Update arrays used for averaging
		update_tracer_averages;

		// Write output array
		write_tracer_output;
	}

	// Write restart file
	write_restart_file;
*/
	return 0;
}
