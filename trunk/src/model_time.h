/*
 * model_time.h
 *
 *  Created on: Oct 16, 2014
 *      Author: ashao
 */

struct time_model_structure {
	// All times are in days since 1900-1-1 00:00:00

	int begin_time;
	int end_time;

	double dt;
	int step_begin_time; // Time at beginning of the step
	int step_end_time;  // Time at the end of the monthly step

	int forcing_midx; // Month index for monthly averaged fields
	int calendar_month; // Month of the current time step (0-January, 1-February)



} model_time;
