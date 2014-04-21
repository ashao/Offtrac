/*
 * oxygen.h
 *
 *  Created on: Apr 15, 2014
 *      Author: ashao
 */

#include "init.h"

/* SUBROUTINE PROTOTYPES */
void allocate_oxygen(  );
void initialize_oxygen( int imon );
void oxygen_saturation(double T[NZ][NXMEM][NYMEM], double S[NZ][NXMEM][NYMEM],
		double o2_sat[NZ][NXMEM][NYMEM]);

void surface_oxygen();
void apply_oxygen_jterms();


# ifdef SPONGE
extern double sp_oxygen[NZ][NXMEM][NYMEM];
# endif


