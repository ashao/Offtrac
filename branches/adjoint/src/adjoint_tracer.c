/*
 * adjoint_tracer.c
 *
 *  Created on: May 1, 2014
 *      Author: ashao
 */
#include "init.h"
#include "read.h"
#include "alloc.h"
#include <stdlib.h>
#include "tracer_utilities.h"
#include "util.h"
int adjointidx;
int madjoint;
double *mn_adjointinv, *adjointinv;
double **mn_adjoint, **adjoint_tracer;
double ***adjoint_init, ***adjoint_source;
extern double ****tr;
extern double misval;
void allocate_adjoint( ){


	int i, j, k;
	// Set index in tracer array
	madjoint = 0;
	mn_adjoint = alloc2d(NXMEM,NYMEM);
	adjoint_tracer = alloc2d(NXMEM,NYMEM);
	adjoint_source = alloc3d(NZ,NXMEM,NYMEM);
	mn_adjointinv = malloc(sizeof(double));
	adjointinv = malloc(sizeof(double));
	adjoint_init = alloc3d(NZ,NXMEM,NYMEM);
	
}

void initialize_adjoint(  ) {

	int i,j,k;
	extern char restart_filename[200];
	extern int oceanmask[NXMEM][NYMEM];

	// TEMPORARY MANUAL SOURC EFOR ADJOINT TTD, REPLACE WITH NETCDF READ
	adjoint_source[18][29][56] = 1.0;
	
	for (i=0;i<NXMEM;i++)
		for (j=0;j<NYMEM;j++)
			for (k=0;k<NZ;k++) {
				if (adjoint_source[k][i][j]>0.) adjoint_init[k][i][j]=1.0;
				else adjoint_init[k][i][j] = 0.0;
			}

	#ifdef RESTART
	printf("Initializing adjoint from restart %s\n",restart_filename);
	read_var3d( restart_filename, "mn_adjoint", 0, adjoint_init);
	#endif

	for (i=0;i<NXMEM;i++)
		for (j=0;j<NYMEM;j++)
			for (k=0;k<NZ;k++)
				tr[madjoint][k][i][j] = adjoint_init[k][i][j];

}

void step_adjoint( ) {

	int i,j,k;
	if (adjointidx<NMONTHSADJOINT) {
		printf("Adding adjoint tracer to source waters\n");
		for (i=0;i<NXMEM;i++)
			for (j=0;j<NYMEM;j++)
				for (k=0;k<NZ;k++) {
					if (adjoint_source[k][i][j]>0)
						tr[madjoint][k][i][j] = adjoint_source[k][i][j];
				}
		adjointidx++;
	}

	printf("Copying adjoint tracer from main tracer array\n");
	for (i=0;i<NXMEM;i++)
		for (j=0;j<NYMEM;j++)
				adjoint_tracer[i][j] = tr[madjoint][0][i][j];

	printf("Setting surface to 0\n");
		for (i=0;i<NXMEM;i++)
			for (j=0;j<NYMEM;j++)
				for (k=0;k<2;k++)
					tr[madjoint][k][i][j] = 0.0;

//	*adjointinv = tracer_inventory(madjoint);
//	printf("Adjoint Inventory: %e\n",&adjointinv);
}
