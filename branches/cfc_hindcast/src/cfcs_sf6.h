/*
 * cfcs_sf6.h
 *
 *  Created on: May 3, 2014
 *      Author: ashao
 */

extern double ***mn_cfc11;
extern double ***cfc11_init;
extern double **cfc11_sat;
extern double **mn_cfc11sat;
extern int mCFC11;

extern double ***mn_cfc12;
extern double ***cfc12_init;
extern double **cfc12_sat;
extern double **mn_cfc12sat;
extern int mCFC12;

extern double ***mn_sf6;
extern double ***sf6_init;
extern double **sf6_sat;
extern double **mn_sf6sat;
extern int mSF6;
#define NUMATMVALS 102
struct tracer_boundary {

	int ntime;
	double *time;
	double *nval;
	double *sval;

};

#define NUMTRANSIENT 3
extern struct tracer_boundary atmconc[NUMTRANSIENT];


