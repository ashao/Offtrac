/*
 * tracadv.new.c
 *
 *  Created on: Jun 20, 2014
 *      Author: ashao
 */


void advect_tracer( ) {

	double ***hprev; 		// The cell volume at the end of
							// the previous tracer change, in m3
	double ***uhr;			// The remaining zonal thickness flux, in m3.
	double ***vhr;			// The remaining meridional thickness fluxes, in m3.

	double ***uh_neglect, ***vh_neglect;
							// uh_neglect and vh_neglect are the magnitude of
							// remaining transports that can be simply discarded
							// in m3 or kg.

	double landvolfill;		// An arbitrary? nonzero cell volume, m3.
	double Idt; 			// 1/dt in s-1.

	int **domore_u;			// domore__ indicate whether there is more
	int **domore_v;			// advection to be done in the corresponding
	int domore_k[NZ];		// row or column

	int x_first;			// If true, advect in the x-direction first.
	int max_iter;			// The maximum number of iterations in each layer.

	int stensil;			// The stensil of the advection scheme.
	int nsten_halo;			// The number of stensils that fit in the halos.
	int i, j, k, m, is, ie, js, je, isd, ied, jsd, jed, nz, itt, ntr, do_any;
	int isv, iev, jsv, jev; // The valid range of the indices.

	is = 0; ie = NXMEM; js = 0; je = NYMEM; nz = NZ;
	isd = 0; ied = NXMEM; jsd = 0; jsd = NYMEM;

	landvolfill = 1.0e-20;
	stensil = 2;

	x_first = 1; // Do zonal transport first

	ntr = NTR;
	Idt = 1.0/dt;

	max_iter = NUM_ADV_ITER;

	hprev = alloc3d(NZ,NXMEM,NYMEM);
	uhr = alloc3d(NZ,NXMEM,NYMEM);
	vhr = alloc3d(NZ,NXMEM,NYMEM);

	domore_u = alloc2d_int(NXMEM,NZ);
	domore_v = alloc2d_int(NYMEM,NZ);

	// This initializes the halos of uhr and vhr because pass_vector might do
	// calculations on them, even though they are never used.
	for (k=0;k<nz;k++)
		for (i=is;i<ie;i++)
			for (j=js;j<je;j++) {
				uhr[k][i][j] = 0.0;
				vhr[k][i][j] = 0.0;
				hprev[k][i][j] = landvolfill;
			}

	for (k=0;k<nz;k++) {
		domore_k[k] = 1;
		//  Put the remaining (total) thickness fluxes into uhr and vhr.
		for (j=js;j<je;j++) {



		}


	}
}
