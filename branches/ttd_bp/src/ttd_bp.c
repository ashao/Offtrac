/*
 * ttd_bp.c
 *
 *  Created on: May 1, 2014
 *      Author: ashao
 */
#include "init.h"
#include "read.h"
#include "alloc.h"
#include <stdlib.h>
int ttdidx;
int mTTD;
double *mn_ttdinv, *ttdinv;
double ***mn_ttd;
double ***ttd_init;
extern double ****tr;
extern double misval;
void allocate_ttd( ){


	int i, j, k;
	// Set index in tracer array
	mTTD = 0;
	mn_ttd = alloc3d(NZ,NXMEM,NYMEM);
	mn_ttdinv = malloc(sizeof(double));
	ttdinv = malloc(sizeof(double));
	ttd_init = alloc3d(NZ,NXMEM,NYMEM);

}

void initialize_ttd(  ) {

	int i,j,k;
	extern char restart_filename[200];
	extern int oceanmask[NXMEM][NYMEM];

	for (i=0;i<NXMEM;i++)
		for (j=0;j<NYMEM;j++)
			for (k=0;k<NZ;k++)
				if (oceanmask[i][j] && k < 2) {
					ttd_init[k][i][j] = 1.0;
				} else
					ttd_init[k][i][j] = 0.0;


	#ifdef RESTART
	printf("Initializing TTD from restart %s\n",restart_filename);
	read_var3d( restart_filename, "mn_ttd", 0, ttd_init);
	#endif

	for (i=0;i<NXMEM;i++)
		for (j=0;j<NYMEM;j++)
			for (k=0;k<NZ;k++)
				tr[mTTD][k][i][j] = ttd_init[k][i][j];


}

void surface_ttd( ) {

	int i,j,k;
	if (ttdidx < NMONTHSTTD){
		ttdidx++;
		printf("Setting surface to 1\n");
		for (i=0;i<NXMEM;i++)
			for (j=0;j<NYMEM;j++)
				for (k=0;k<2;k++)
					tr[mTTD][k][i][j] = ttd_init[k][i][j];
	}
	else {
	/*
		for (i=0;i<NXMEM;i++)
			for (j=0;j<NYMEM;j++)
				for (k=0;k<2;k++)
					tr[mTTD][k][i][j] = 0.0;
	*/	
	}

	*ttdinv = tracer_inventory(mTTD);
	printf("TTD Inventory: %e\n",&ttdinv);
}
