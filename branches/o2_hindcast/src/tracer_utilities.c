/*
 * tracer_utilities.c
 *
 *  Created on: Apr 15, 2014
 *      Author: ashao
 */
#include "init.h"
#include "tracer_utilities.h"
extern double D[NXMEM][NYMEM];
void conc_obs_layer(double h[NZ][NXMEM][NYMEM],
		double conc_lev[NZWOA][NXMEM][NYMEM],
		double conc_lay[NZ][NXMEM][NYMEM]) {
	int i, j, k;
	int nzlevitus = NZWOA;
	double depth[NZ][NXMEM][NYMEM];
	double concobsprof[NZWOA];
	double levitus_depths[NZWOA] = { 0, 10, 20, 30, 50, 75, 100, 120, 150,
			200, 250, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200,
			1300, 1400, 1500, 1750, 2000, 2500, 3000, 3500, 4000, 4500, 5000,
			5500 };

	/*---------------------------------------------------------------------
	 *     calculate vertically interpolated concentration:
	 *            levitus levels -> model layers
	 *             conc_lev  -> conc_lay
	 *---------------------------------------------------------------------*/
	z_depth(h, depth);

	//printf("Initialize 2a for month \n");
	//HF for (i=X1;i<=nx;i++) {
	//HF for (j=Y1;j<=ny;j++) {
	for (i = 0; i <= NXMEM - 1; i++) {
		for (j = 0; j <= NYMEM - 1; j++) {
			//BX - reinstated by HF
			if (D[i][j]>MINIMUM_DEPTH) {
				for (k = 0; k < nzlevitus; k++) {
					concobsprof[k] = conc_lev[k][i][j];
				}
				for (k = 0; k < NZ; k++) {
					conc_lay[k][i][j] = lin_interp(depth[k][i][j], concobsprof,
							levitus_depths, 0, nzlevitus);
					if (conc_lay[k][i][j] < 0.e0)
						conc_lay[k][i][j] = 0.;
				}
				//BX - reinstated by HF
			} else {
				for (k = 0; k < NZ; k++) {
					conc_lay[k][i][j] = 0.0;
				}
			}
		}
	}
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

	for (k=0;k<NZ;k++)
		printf("depth[%d]=%f\n",depth[k][127][127]);

}
