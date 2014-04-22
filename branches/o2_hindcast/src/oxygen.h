/*
 * oxygen.h
 *
 *  Created on: Apr 15, 2014
 *      Author: ashao
 */

#include "init.h"
extern int mOXYGEN;
// Output arrays
extern double ***mn_oxygen;
extern double ***mn_o2sat;
extern double ***mn_jo2;
extern double mn_oxyflux[NXMEM][NYMEM];
// Working arrays
extern double ***oxy_init;
extern double o2_sat[NZ][NXMEM][NYMEM];
extern double oxyflux[NXMEM][NYMEM];
extern double jo2[NZ][NXMEM][NYMEM];

/* SUBROUTINE PROTOTYPES */
void allocate_oxygen(  );
void initialize_oxygen( int imon );
void oxygen_saturation(double ***T, double ***S,
		double o2_sat[NZ][NXMEM][NYMEM]);

void surface_oxygen();
void apply_oxygen_jterms();


# ifdef SPONGE
extern double sp_oxygen[NZ][NXMEM][NYMEM];
# endif


