/*
 * ttd_bp.c
 * Contains all subroutines for needed to estimate the TTD using the boundary propagator method
 *  Created on: Jan 9, 2014
 *      Author: ashao
 */
#define DAYSINYEAR 365.00

void initialize_boundary_propagator( int imon )
{

	int i,j,k;
	double volume, mldepth,totvolume;
	extern double areagr[NXMEM][NYMEM];
	extern double hstart[NZ][NXMEM][NYMEM];
	extern int wetmask[NXMEM][NYMEM];
	extern int mlen[12];
	extern int mTTDBP;
	extern double ****tr;


	totvolume = 0;



	for ( i=0; i<NXMEM; i++ ) {
		for ( j=0; j<NYMEM; j++ ) {

			mldepth = hstart[0][i][j]+hstart[1][i][j];
			volume = mldepth * areagr[i][j] * (double) wetmask[i][j];
			totvolume += volume;
			tr[mTTDBP][0][i][j] = volume;
			tr[mTTDBP][1][i][j] = volume;

		}
	}

	for ( i=0; i<NXMEM; i++ ) {
		for ( j=0; j<NYMEM; j++ ) {

			tr[mTTDBP][0][i][j] = tr[mTTDBP][0][i][j] * mlen[imon] / DAYSINYEAR / totvolume;
			tr[mTTDBP][1][i][j] = tr[mTTDBP][1][i][j] * mlen[imon] / DAYSINYEAR / totvolume;

		}
	}
}
