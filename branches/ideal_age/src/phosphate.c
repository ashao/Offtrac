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
// double po4[NZ][NXMEM][NYMEM];
double jpo4[NZ][NXMEM][NYMEM];
double ***po4_star_lay;
double jprod[NZ][NXMEM][NYMEM];
double jremin[NZ][NXMEM][NYMEM];
double jremdop[NZ][NXMEM][NYMEM];
double jdop[NZ][NXMEM][NYMEM];
double flux_pop[NXMEM][NYMEM];
double ***dop_init;
void allocate_phosphate(  )
{
	int i,j,k;
	// Base the indices of the main tracer array off of the oxygen index
	extern int mOXYGEN;
	mDOP = mOXYGEN+1;
	mPHOSPHATE = mDOP+1;
	printf("mDOP: %d mPHOSPHATE: %d\n",mDOP,mPHOSPHATE);
	mn_phos = alloc3d(NZ,NXMEM,NYMEM);
	mn_dop = alloc3d(NZ,NXMEM,NYMEM);
	mn_pobs = alloc3d(NZ,NXMEM,NYMEM);
	mn_jpo4 = alloc3d(NZ,NXMEM,NYMEM);
	mn_jremin = alloc3d(NZ,NXMEM,NYMEM);
	mn_jprod = alloc3d(NZ,NXMEM,NYMEM);
	phosphate_init = alloc3d(NZ,NXMEM,NYMEM);
	dop_init = alloc3d(NZ,NXMEM,NYMEM);
	po4_star_lay = alloc3d(NZ,NXMEM,NYMEM);
}
void initialize_phosphate( int imon ) {
	int i,j,k;
	extern double hend[NZ][NXMEM][NYMEM];
	extern double ****tr;
	extern char restart_filename[200];
#ifdef WOA_PHOS
	printf("Initializing phosphate from WOA09\n");
	read_woa_file(imon, hend, phosphate_init, "woa09.phos.nc", "p_an",1e-3);
	printf("Initializing DOP to zero\n");
	set_darray3d_zero(dop_init, NZ, NXMEM, NYMEM);
#endif

#ifdef RESTART
	printf("Initializing phosphate from restart: %s\n",restart_filename);
	read_var3d( restart_filename, "mn_phos", 0, phosphate_init);
	printf("Initializing dop from restart: %s\n",restart_filename);
	read_var3d( restart_filename, "mn_dop", 0, dop_init);
#endif

	for (k=0;k<NZ;k++)
		for (i=0;i<NXMEM;i++)
			for(j=0;j<NYMEM;j++) {
			if (dop_init[k][i][j] > 0.0) tr[mDOP][k][i][j] = dop_init[k][i][j];
			else tr[mDOP][k][i][j] = 0.0;
			if (phosphate_init[k][i][j]>0.0) tr[mPHOSPHATE][k][i][j] = phosphate_init[k][i][j];
			else tr[mPHOSPHATE][k][i][j] = 0.0;

			}
	free3d(dop_init,NZ);
	free3d(phosphate_init,NZ);
	printf("Phosphate Initialization:%e\n",tr[mPHOSPHATE][20][100][100]);
}

void apply_phosphate_jterms( ) {
	int i,j,k;
	extern double dt;
	extern double ****tr;
	extern int oceanmask[NXMEM][NYMEM];
	extern double hend[NZ][NXMEM][NYMEM];
	extern double misval;
	// j terms here are calculated from biotic_sms routine in biotic.c
	printf("Applying j terms for phosphate\n");
	printf("Example jpo4/po4: %e/%f\n",jpo4[10][127][127],tr[mPHOSPHATE][10][127][127]);
	for (i = 0; i <= NXMEM - 1; i++) {
		for (j = 0; j <= NYMEM - 1; j++) {
			//BX - reinstated by HF
			if (oceanmask[i][j]) {
				for (k = 0; k < NZ; k++) {
					if (hend[k][i][j]>EPSILON) tr[mPHOSPHATE][k][i][j] += dt * jpo4[k][i][j];
					if (hend[k][i][j]>EPSILON) tr[mDOP][k][i][j] += dt * jdop[k][i][j];
	//				if (tr[mPHOSPHATE][k][i][j] > 0.1) tr[mPHOSPHATE][k][i][j] = 0.1;
	//				if (tr[mPHOSPHATE][k][i][j] < 0.0) tr[mPHOSPHATE][k][i][j] = 0.0;
	//				if (tr[mDOP][k][i][j] > 0.1) tr[mDOP][k][i][j] = 0.1;
	//				if (tr[mDOP][k][i][j] < 0.0) tr[mDOP][k][i][j] = 0.0;
					}
				
			} else {

				for (k = 0; k < NZ; k++) {
				tr[mPHOSPHATE][k][i][j] = misval;
				tr[mDOP][k][i][j] = misval;
				}
			}
		}
	}


}
