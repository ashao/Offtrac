/*
 * tracer_utilities.h
 *
 *  Created on: Apr 15, 2014
 *      Author: ashao
 */


void z_depth(double h[NZ][NXMEM][NYMEM], double depth[NZ][NXMEM][NYMEM]);
void allocate_ts( );
void read_temp_and_salt( int imon, char *fieldtype);
double linear_interpolation(const double xin[], const double yin[], double xi, int numin, double *yout);

