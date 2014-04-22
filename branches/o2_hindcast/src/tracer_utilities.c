/*
 * tracer_utilities.c
 *
 *  Created on: Apr 15, 2014
 *      Author: ashao
 */
#include "init.h"
#include "tracer_utilities.h"
#include "alloc.h"

extern double D[NXMEM][NYMEM];
double ***Temptm,***Salttm;

void allocate_ts( ) {

	Temptm = alloc3d(NZ,NXMEM,NYMEM);
	Salttm = alloc3d(NZ,NXMEM,NYMEM);

}

void read_temp_and_salt( int imon, char *fieldtype) {
	extern char directory[75];
	char inpath[200];
	strcpy(directory,inpath);
	strcat(inpath,sprintf("ts-%s.nc",fieldtype));

	read_var3d( inpath, "temp", imon, Temptm);
	read_var3d( inpath, "salt", imon, Salttm);

}

void z_depth(double h[NZ][NXMEM][NYMEM], double depth[NZ][NXMEM][NYMEM]) {
	/* compute depth in meters for use in tracer depth dependent functions
	 * Ivan Lima - Nov 2002 */
	int i, j, k;
	double hsum;
	/* initialize variables with zero */

	for (k = 0; k < NZ; k++) {
		for (i = X1; i <= nx; i++) {
			for (j = Y1; j <= ny; j++) {
				depth[k][i][j] = 0.0;
			}
		}
	}
	/* compute depth at layer k as the the depth of the point in the
	 * middle of that layer */
	for (i = X1; i <NXMEM; i++) {
		for (j = Y1; j <NYMEM; j++) {
			//BX - reinstated by HF
			if (D[i][j]>MINIMUM_DEPTH) {
				hsum = h[0][i][j];
				depth[0][i][j] = h[0][i][j] * 0.5;
				for (k = 1; k < NZ; k++) {
					depth[k][i][j] = hsum + h[k][i][j] * 0.5;
					hsum += h[k][i][j];
				}
				//BX - reinstated by HF
			} else {
				for (k = 0; k < NZ; k++) {
					depth[k][i][j] = 0.0;
				}
			}
		}
	}


}
