/*
 * tracer_utilities.h
 *
 *  Created on: Apr 15, 2014
 *      Author: ashao
 */

void conc_obs_layer(double h[NZ][NXMEM][NYMEM],double conc_lev[NZWOA][NXMEM][NYMEM],
		double conc_lay[NZ][NXMEM][NYMEM]);
void z_depth(double h[NZ][NXMEM][NYMEM], double depth[NZ][NXMEM][NYMEM]);

