#include <stdio.h>
#include "init.h"
#include "netcdf.h"
#include "alloc.h"
#include "init.h"
#include <math.h>
/* OXYGEN VARIABLE DECLARATIONS */
// Auxiliary variables
int mAGE;
// Output arrays
double ***mn_age;
// Working arrays
double ***age_init;
extern double ****tr;
extern int oceanmask[NXMEM][NYMEM];

void allocate_age (  ) {
	printf("Allocating age arrays\n");

	int i, j, k;

	// Set index in tracer array
	mAGE = 0;

	// Allocate working and output arrays
	mn_age = alloc3d(NZ,NXMEM,NYMEM);
	age_init = alloc3d(NZ,NXMEM,NYMEM);
}

void initialize_age( ) {
	int i,j,k;
	extern char restart_filename[200];
	extern double misval;

#ifdef RESTART
	printf("Initializing oxygen from restart %s\n",restart_filename);
	read_var3d( restart_filename, "mn_age", 0, age_init);
#endif
	printf("Age init example: %f\n",age_init[10][100][100]);
	// Copy the initialized tracer value over to main trace array
	for (i=0;i<NXMEM;i++) {
		for (j=0;j<NYMEM;j++) {
				for (k=0;k<NZ;k++) {
					tr[mAGE][k][i][j] = age_init[k][i][j];
			}
		}
	}

	free3d(age_init,NZ);

}

void step_age( double dt ){

	int i,j,k;
	const double numsecsinyear = 365.0*86400.0;

	for(i=0;i<NXMEM;i++)
		for(j=0;j<NYMEM;j++)
			if(oceanmask[i][j]) {
				// Set the mixed layer to zero
				tr[mAGE][0][i][j] = 0.0;
				tr[mAGE][1][i][j] = 0.0;
				for(k=2;k<NZ;k++)
					tr[mAGE][k][i][j] += dt/numsecsinyear;
			}


	for(i=0;i<NXMEM;i++)
		for(j=0;j<NYMEM;j++)
			for(k=0;k<NZ;k++)
				mn_age[k][i][j]+=tr[mAGE][k][i][j];
}
