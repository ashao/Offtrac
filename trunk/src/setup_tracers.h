/*
 * setup_tracers.h
 *
 *  Created on: Oct 16, 2014
 *      Author: ashao
 */
extern struct vardesc vars[NOVARS];
extern int flags[NOVARS];
extern int rflags[NOVARS];

void initialize_vardesc();
void set_output_flags();
void set_restart_flags();
void calculate_tracer_averages();
void write_variables_to_netcdf();
void write_restart_netcdf();
