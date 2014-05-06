/*
 * sf6.c
 *
 *  Created on: May 2, 2014
 *      Author: ashao
 */
#include "init.h"
#include "alloc.h"
#include "read.h"
#include "math.h"
#include "cfcs_sf6.h"
double ***mn_sf6;
double ***sf6_init;
double **sf6_sat;
double **mn_sf6sat;
double **sf6_atmconc;
int mSF6;
extern double ****tr;

void allocate_sf6 ( ) {

	mn_sf6 = alloc3d(NZ,NXMEM,NYMEM);
	sf6_init = alloc3d(NZ,NXMEM,NYMEM);
	sf6_sat = alloc2d(NXMEM,NYMEM);
	mn_sf6sat = alloc2d(NXMEM,NYMEM);
	sf6_atmconc = alloc2d(NXMEM,NYMEM);
}

void initialize_sf6 ( ) {
	int i, j, k;
	
	mSF6 = mCFC12+1;
	printf("mSF6: %d\n",mSF6);
	extern char restart_filename[200];

#ifdef RESTART
	printf("Initializing CFC-11 from restart %s\n",restart_filename);
	read_var3d( restart_filename, "mn_sf6", 0, sf6_init);
#else
	set_darray3d_zero(sf6_init,NZ,NXMEM,NYMEM);
#endif

	for(i=0;i<NXMEM;i++)
		for(j=0;j<NYMEM;j++)
			for(k=0;k<NZ;k++)
				tr[mSF6][k][i][j] = sf6_init[k][i][j];

	free3d(sf6_init,NZ);
}

void sf6_saturation( double **sat ) {

	const double solcoeffs[6] = {-80.0343,117.232,29.5817,0.0335183,-0.0373942,0.00774862};
	int i, j;
	double work;
	double TempK;
	extern double ***Temptm, ***Salttm;

	for (i=0;i<NXMEM;i++)
		for (j=0;j<NYMEM;j++)	{
					TempK = Temptm[0][i][j]+273.15;
					work = solcoeffs[0] + solcoeffs[1]*(100/TempK) +
						solcoeffs[2]*log( TempK/100 ) +
						Salttm[0][i][j]*(solcoeffs[3]+solcoeffs[4]*(TempK/100)+solcoeffs[5]*pow(TempK/100,2));
					sat[i][j] = exp(work)*sf6_atmconc[i][j];
			}



}

void sf6_find_atmconc(  ) {

	int i,j;
	extern double **geolat;
	extern double **atmpres;
	const double equatorbound[2] = {10,-10}; // Set the latitudes where to start interpolating atmospheric concentrations
	extern double currtime;
	double hemisphere_concentrations[2];

	// Interpolate in time to find the atmospheric concentration
	hemisphere_concentrations[0] = linear_interpolation(
				atmconc[mSF6].time, atmconc[mSF6].nval, currtime,atmconc[mSF6].ntime);
		hemisphere_concentrations[1] = linear_interpolation(
				atmconc[mSF6].time, atmconc[mSF6].sval, currtime,atmconc[mSF6].ntime);


	for (i=0;i<NXMEM;i++)
		for (j=0;j<NYMEM;j++) {

			if (geolat[i][j] < equatorbound[0] && geolat[i][j] > equatorbound[1]) {
				sf6_atmconc[i][j] = lin_interpp(
						equatorbound,hemisphere_concentrations,geolat[i][j],2);
			}
			if (geolat[i][j]>equatorbound[0] ) {
				sf6_atmconc[i][j] = hemisphere_concentrations[0];
			}
			if (geolat[i][j]<-equatorbound[1] ) {
				sf6_atmconc[i][j] = hemisphere_concentrations[1];
			}
		}

}

void surface_sf6( ) {

	int i,j,k;

	printf("Setting SF6 surface condition\n");
	// Set oxygen values to saturation at the mixed layer to mimic equilibrium with the atmosphere
	sf6_find_atmconc( );
	sf6_saturation( sf6_sat);
	for (k=0;k<NML;k++)
		for (i=0;i<NXMEM;i++)
			for (j=0;j<NYMEM;j++)
					tr[mSF6][k][i][j]=sf6_sat[i][j];

}

