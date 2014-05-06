/*
 * cfc12.c
 *
 *  Created on: May 2, 2014
 *      Author: ashao
 */
#include "init.h"
#include "alloc.h"
#include "read.h"
#include "math.h"
#include "cfcs_sf6.h"
double ***mn_cfc12;
double ***cfc12_init;
double **cfc12_sat;
double **mn_cfc12sat;
double **cfc12_atmconc;

int mCFC12;
extern double ****tr;

void allocate_cfc12 ( ) {

	mn_cfc12 = alloc3d(NZ,NXMEM,NYMEM);
	cfc12_init = alloc3d(NZ,NXMEM,NYMEM);
	cfc12_sat = alloc2d(NXMEM,NYMEM);
	mn_cfc12sat = alloc2d(NXMEM,NYMEM);
	cfc12_atmconc = alloc2d(NXMEM,NYMEM);
}

void initialize_cfc12 ( ) {
	int i, j, k;
	extern char restart_filename[200];
	mCFC12 = mCFC11+1;
	printf("mCFC12: %d\n",mCFC12);
#ifdef RESTART
	printf("Initializing CFC-11 from restart %s\n",restart_filename);
	read_var3d( restart_filename, "mn_cfc12", 0, cfc12_init);
#else
	set_darray3d_zero(cfc12_init,NZ,NXMEM,NYMEM);
#endif

	for(i=0;i<NXMEM;i++)
		for(j=0;j<NYMEM;j++)
			for(k=0;k<NZ;k++)
				tr[mCFC12][k][i][j] = cfc12_init[k][i][j];

	free3d(cfc12_init,NZ);
}

void cfc12_saturation( double **sat ) {

	const double solcoeffs[7] = {-218.0971,298.9702,113.8049,-1.39165,-0.143566,0.091015,-0.0153924};
	int i, j;
	double work;
	double TempK;
	extern double ***Temptm, ***Salttm;

	for (i=0;i<NXMEM;i++)
		for (j=0;j<NYMEM;j++)	{
			TempK = Temptm[0][i][j]+273.15;
			work = solcoeffs[0] + solcoeffs[1]*(100/TempK) +solcoeffs[2]*log( TempK/100 ) +
				solcoeffs[3]*pow(TempK/100,2) + Salttm[0][i][j]*(solcoeffs[4]+
				solcoeffs[5]*(TempK/100)+solcoeffs[6]*pow(TempK/100,2));
				sat[i][j] = exp(work)*cfc12_atmconc[i][j];
			}



}

void cfc12_find_atmconc(  ) {

	int i,j;
	extern double **geolat;
	extern double **atmpres;
	const double equatorbound[2] = {10,-10}; // Set the latitudes where to start interpolating atmospheric concentrations
	extern double currtime;
	double hemisphere_concentrations[2];

	// Interpolate in time to find the atmospheric concentration
	hemisphere_concentrations[0] = linear_interpolation(
				atmconc[mCFC12].time, atmconc[mCFC12].nval, currtime,atmconc[mCFC12].ntime);
		hemisphere_concentrations[1] = linear_interpolation(
				atmconc[mCFC12].time, atmconc[mCFC12].sval, currtime,atmconc[mCFC12].ntime);

	for (i=0;i<NXMEM;i++)
		for (j=0;j<NYMEM;j++) {

			if (geolat[i][j] < equatorbound[0] && geolat[i][j] > equatorbound[1]) {
				cfc12_atmconc[i][j] = linear_interpolation(
						equatorbound,hemisphere_concentrations,geolat[i][j],2);
			}
			if (geolat[i][j]>equatorbound[0] ) {
				cfc12_atmconc[i][j] = hemisphere_concentrations[0];
			}
			if (geolat[i][j]<-equatorbound[1] ) {
				cfc12_atmconc[i][j] = hemisphere_concentrations[1];
			}
		}
	printf("cfc12_atmconc[100][100] = %f\n",cfc12_atmconc[100][100]);
}

void surface_cfc12( ) {

	int i,j,k;

	printf("Setting CFC-12 surface condition\n");
	// Set oxygen values to saturation at the mixed layer to mimic equilibrium with the atmosphere
	cfc12_find_atmconc( );
	cfc12_saturation( cfc12_sat);
	for (k=0;k<NML;k++)
		for (i=0;i<NXMEM;i++)
			for (j=0;j<NYMEM;j++)
					tr[mCFC12][k][i][j]=cfc12_sat[i][j];

}

