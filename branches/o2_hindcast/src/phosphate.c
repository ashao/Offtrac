/*
 * phosphate.c
 *
 *  Created on: Apr 15, 2014
 *      Author: ashao
 */
#include "init.h"
#include "netcdf.h"
#include "phosphate.h"
#include "alloc.h"
#include <stdio.h>
/* PHOSPHATE VARIABLE DECLARATIONS */
// Auxiliary variables
int mPHOSPHATE;
int mDOP;
// Output arrays
double ***mn_phos;
double ***mn_dop;
double ***mn_pobs;
double ***mn_jpo4;
double ***mn_jremin;
double ***mn_jremdop;
double mn_pflux[NXMEM][NYMEM];
double ***mn_jprod;

// Working arrays
double ***phosphate_init;
double po4[NZ][NXMEM][NYMEM];
double jpo4[NZ][NXMEM][NYMEM];
double ***po4_star_lay;
double jprod[NZ][NXMEM][NYMEM];
double jremin[NZ][NXMEM][NYMEM];
double jremdop[NZ][NXMEM][NYMEM];
double jdop[NZ][NXMEM][NYMEM];
double flux_pop[NXMEM][NYMEM];
double ***dop_init;
void initialize_phosphate( int imon )
{
	
	extern char restart_filename[200];
	// Base the indices of the main tracer array off of the oxygen index
	extern int mOXYGEN;
	extern double hstart[NZ][NXMEM][NYMEM];
	mDOP = mOXYGEN++;
	mPHOSPHATE = mDOP++;

	mn_phos = alloc3d(NZ,NXMEM,NYMEM);
	mn_dop = alloc3d(NZ,NXMEM,NYMEM);
	mn_pobs = alloc3d(NZ,NXMEM,NYMEM);
	mn_jpo4 = alloc3d(NZ,NXMEM,NYMEM);
	mn_jremin = alloc3d(NZ,NXMEM,NYMEM);
	mn_jprod = alloc3d(NZ,NXMEM,NYMEM);
	phosphate_init = alloc3d(NZ,NXMEM,NYMEM);
	dop_init = alloc3d(NZ,NXMEM,NYMEM);
	po4_star_lay = alloc3d(NZ,NXMEM,NYMEM);
#ifdef WOA_PHOS
	printf("Initializing phosphate from WOA09\n");
	read_woa_file(imon, hstart, phosphate_init, "woalevpo4.nc", "WOAPO4");
	printf("Initializing DOP to zero\n");
	set_darray3d_zero(dop_init, NZ, NXMEM, NYMEM);
#endif

#ifdef RESTART
	printf("Initializing phosphate from restart: %s\n",restart_filename);
	read_var3d( restart_filename, "mn_phos", imon, double phosphate_init);
	printf("Initializing dop from restart: %s\n",restart_filename);
	read_var3d( restart_filename, "mn_dop", imon, double dop_init);
#endif
}

void apply_phosphate_jterms( ) {
	int i,j,k;
	extern double dt;
	extern double ****tr;
	extern double D[NXMEM][NYMEM];
	// j terms here are calculated from biotic_sms routine in biotic.c
	for (i = 0; i <= NXMEM - 1; i++) {
		for (j = 0; j <= NYMEM - 1; j++) {
			//BX - reinstated by HF
			if (D[i][j]>MINIMUM_DEPTH) {
				for (k = 0; k < NZ; k++) {
					tr[mPHOSPHATE][k][i][j] += dt * jpo4[k][i][j];
					tr[mDOP][k][i][j] += dt * jdop[k][i][j];
				}
			} else {
				tr[mPHOSPHATE][k][i][j] = 0.0;
				tr[mDOP][k][i][j] = 0.0;
				jpo4[k][i][j] = 0.0;
			}
		}
	}


}