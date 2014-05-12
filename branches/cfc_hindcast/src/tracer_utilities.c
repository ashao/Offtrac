/*
 * tracer_utilities.c
 *
 *  Created on: Apr 15, 2014
 *      Author: ashao
 */
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "init.h"
#include "tracer_utilities.h"
#include "alloc.h"
#include "read.h"
#include <math.h>
extern int oceanmask[NXMEM][NYMEM];
extern double D[NXMEM][NYMEM];
double ***Temptm,***Salttm;

void allocate_ts( ) {

	Temptm = alloc3d(NZ,NXMEM,NYMEM);
	Salttm = alloc3d(NZ,NXMEM,NYMEM);

}

void read_temp_and_salt( int imon, char *fieldtype) {
	extern char directory[75];
	char filename[20];
	char saltpath[300];
	char temppath[300];

	if ( imon % NMONTHS == 0 )
		imon = 0;

	strcpy(saltpath,directory);
	sprintf(filename,"salt.%s.nc",fieldtype);
	strcat(saltpath,filename);
	strcpy(temppath,directory);
	sprintf(filename,"temp.%s.nc",fieldtype);
	strcat(temppath,filename);
	printf("Reading temperature and salinity from month %d\n",imon);
	read_var3d( temppath, "temp", imon, Temptm);
	read_var3d( saltpath, "salt", imon, Salttm);

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
			if (oceanmask[i][j]) {
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


double linear_interpolation(const double xin[], const double yin[], double xi, int numin, double *yout) {

	int i,j,flipidx;
	int intpidx1,intpidx2;
	int decreasing;
	double deltax,dist;
	double *x,*y;
	double y0,y1,x0,x1,yi;

	x = (double *) malloc(numin * sizeof(double));
	y = (double *) malloc(numin * sizeof(double));

	// first check to see if xin increases or decreases
	if (xin[numin-1]>xin[0]) decreasing = 0;
	if (xin[numin-1]<xin[0]) decreasing = 1;
	
	// flip the input vectors if it is decreasing
	if (decreasing) {
		printf("Flipping vector\n");
		for (flipidx = 0,i=numin-1; i>=0; i--,flipidx++) {
			x[flipidx] = *(xin+i);
			y[flipidx] = *(yin+i);
		}
	}
	if (!decreasing) {
		for (i=0;i<numin;i++) {
			x[i] = xin[i];
			y[i] = yin[i];
		}	
	}
	
	// Figure out which elements of the vector to use
	deltax = fabs(xi-x[0]);
	intpidx1 = 0;
	intpidx2 = 1;
	for (i=0;i<numin;i++) {
		dist = fabs(xi-x[i]); // Calculate how far away the current x value is from the desired number
		if (dist<deltax) {
			deltax = dist;
			intpidx1 = i; // We now know at least one bound of the interpolation interval
			if ( xi<x[i] ){
				intpidx2 = i-1; 
			}
			else{
				intpidx2 = i + 1;
			}
		}
	}
	y0 = y[ intpidx1 ];
	y1 = y[ intpidx2 ];
	x0 = x[ intpidx1 ];
	x1 = x[ intpidx2 ];
	yi = y0 + (y1-y0)*(xi-x0)/(x1-x0);
	yout = &yi+1;
	printf("y0: %f y1: %f yi: %f x0: %f x1: %f xi: %f\n",y0,y1,yi,x0,x1,xi);

	free(x);
	free(y);
	return(yi);

}
