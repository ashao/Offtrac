/********+*********+*********+*********+*********+*********+*********+*
 *      							      *
 *                source/sink terms for offtrac tracers               * 
 *                                                                    *
 ********+*********+*********+*********+*********+*********+*********+*/

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "io.h"
#include "init.h"
#include "metrics.h"
#include "alloc.h"
#include "biotic.h"
#include "step.h"
#include "tracadv.h"
#include "util.h"

/*---------------------------------------------------------------------
 *     define variables and subroutines
 *---------------------------------------------------------------------*/

void oxygen_saturation(double T[NZ][NXMEM][NYMEM], double S[NZ][NXMEM][NYMEM],
		double o2_sat[NZ][NXMEM][NYMEM]);
void co2_saturation(double T[NZ][NXMEM][NYMEM], double S[NZ][NXMEM][NYMEM],
		double co2_sat[NZ][NXMEM][NYMEM]);
void argon_saturation(double T[NZ][NXMEM][NYMEM], double S[NZ][NXMEM][NYMEM],
		double ar_sat[NZ][NXMEM][NYMEM]);
//void cfc_solubility(double T[NZ][NXMEM][NYMEM],double S[NZ][NXMEM][NYMEM],
//                    int method, double cfc_sol[NZ][NXMEM][NYMEM], int icfc);
void
		co2calc(double t, double s, double dic_in, double ta_in, double pt_in,
				double sit_in, double phlo, double phhi, double xco2_in,
				double atmpres);

void z_depth(double h[NZ][NXMEM][NYMEM], double depth[NZ][NXMEM][NYMEM]);
void calc_ksub(double h[NZ][NXMEM][NYMEM]);
void conc_obs_layer(double h[NZ][NXMEM][NYMEM],
		double conc_lev[NZPHOS][NXMEM][NYMEM],
		double conc_lay[NZ][NXMEM][NYMEM]);
void tracer_integral(int trnum, double ***hvol);
void tracer_profile(double tr3d[NZ][NXMEM][NYMEM]);
void print_to_screen(int pquant[10], int pstage);
void isotope_R_del(int trnum1, int trnum2, double Rstandard,
		double R[NZ][NXMEM][NYMEM], double del[NZ][NXMEM][NYMEM]);
#ifdef MERGED_ML
void merge_ml_tr();
void merge_ml_j();
#endif 
//double cfc_atmospheric(double now, int icfc);

extern double hend[NZ][NXMEM][NYMEM];
extern double h[NZ][NXMEM][NYMEM];
extern double h_clim[NZ][NXMEM][NYMEM];
extern double depth[NZ][NXMEM][NYMEM];
extern double rml[2][NXMEM][NYMEM];
extern double D[NXMEM][NYMEM];
extern double ***uhtm, ***vhtm;
extern double Temptm[NZ][NXMEM][NYMEM];
extern double Salttm[NZ][NXMEM][NYMEM];
extern double ***ea, ***eb;
extern double windspeed[NXMEM][NYMEM];
extern double xkw[NXMEM][NYMEM];
extern double atmpres[NXMEM][NYMEM];
extern double fice[NXMEM][NYMEM];
extern double salt_woa[NXMEM][NYMEM];
extern double ****tr;
extern double dt;
extern double day, daystart, day1948;
extern double trwarn[NTR][2];
extern double trintegral[NTR];
extern double trprofile[NZ];
extern int ksub[NXMEM][NYMEM];
extern int ksub_clim[NXMEM][NYMEM];
extern double currtime; //ashao: get the current time
/*   Begin added DT    */
extern int beginyear;
extern int theyear;
extern int lastsave;
extern int estTTD;
/*   End added DT      */
int month;
int ii;
int pquant[10], pstage, doprint;

#ifdef AGE
extern double ***age;
extern int mAGE;
# ifdef SPONGE
extern double sp_age[NZ][NXMEM][NYMEM];
# endif
#endif



# ifdef SPONGE
extern double sp_age[NZ][NXMEM][NYMEM];
# endif

// begin ashao

double dt_gasex;
double h_mixedlyr;

// end ashao




#ifdef SPONGE
extern double Idamp[NXMEM][NYMEM];
extern double sp_e[NZ][NXMEM][NYMEM];
#endif

extern double misval;

/*---------------------------------------------------------------------
 *
 *     BEGIN MAIN EXECUTABLE
 *        Each tracer is stepped forward in 3 consecutive steps:
 *         1) physical transport
 *         2) biological/chemical interior source-sink
 *         3) fluxes across air-sea interface
 *         
 *---------------------------------------------------------------------*/

void step_fields(int iyear, int itts, int imon, int iterno) {
	int i, j, k, l, m;
	int mdt;
	double inv_mdt;
	int ibiodt;

	//DT
	double Tmptr[NXMEM][NYMEM];
	int n;
	double secperyr = 365.25 * 60 * 60 * 24;
	theyear = iyear;
	//DT

	// begin ashao
	double thetime;
	extern double lath[NXTOT];
	// end ashao


	/*-----------------------------------------
	 *
	 *     PARAMETER VALUES
	 *
	 *-----------------------------------------*/

	/*-----------------------------------------
	 *
	 *     SCREEN PRINT CONTROL
	 *
	 *-----------------------------------------*/

	doprint = 0;
	pquant[1] = 0; /* print mass conservation checks? (1 = yes) */
	pquant[2] = 1; /* print single grid point at each step? (1 = yes) */
	pquant[3] = 1; /* print negative concentrations? */
	/*   pq3=1 => print warning                */
	/*   pq3=0 => set concentration to trwarn   */
	/*   pq3=2 => do both                      */
	/*   pq3=3 => do nothing                   */

	/*-----------------------------------------
	 *
	 *     CALCULATE PRELIMINARY QUANTITIES
	 *
	 *-----------------------------------------*/

	z_depth(h, depth);


	/*-----------------------------------------
	 *
	 *     PRINT TO SCREEN #1
	 *
	 *-----------------------------------------*/

	pstage = 1;
	print_to_screen(pquant, pstage);

	/*-----------------------------------------
	 *
	 *     STEP 1: PHYSICAL TRANSPORT
	 *
	 *-----------------------------------------*/



	printf("Calculate tracer transport. \n");


	tracer(itts); /* perform transport time step */

#ifdef MERGED_ML
	merge_ml_tr();
#endif
	/*-----------------------------------------
	 *
	 *     PRINT TO SCREEN #2
	 *
	 *-----------------------------------------*/


	pstage = 2;
	print_to_screen(pquant, pstage);

	/*-----------------------------------------
	 *
	 *     STEP 2: BIOTIC SOURCE/SINK TERMS
	 *
	 *-----------------------------------------*/

	/*-----------------------------------------
	 *
	 *  CALCULATE BIOTIC SOURCE/SINK TERMS
	 *    based on OCMIP scheme in which
	 *    biological pump is linked to PO4
	 *
	 *-----------------------------------------*/

	ibiodt = 1; // number of biotic time steps per transport time step (dt)

	printf(
			"Calculate biotic source and sink terms with %i time step(s) per dt. \n",
			ibiodt);



# ifdef MERGED_ML
	merge_ml_j();
# endif

	/*-----------------------------------------
	 *
	 *  APPLY BIOTIC SOURCE/SINK TERMS
	 *   step tr array forward using Euler method
	 *   calculate global integral for mass conservation checking
	 *
	 *-----------------------------------------*/

	for (i = 0; i <= NXMEM - 1; i++) {
		for (j = 0; j <= NYMEM - 1; j++) {
			//BX - reinstated by HF
			if (D[i][j] > MINIMUM_DEPTH) {
				for (k = 0; k < NZ; k++) {
				}
				//BX - reinstated by HF
			} else {
				for (k = 0; k < NZ; k++) {
#ifdef AGE
# if defined AGEEXP || defined AGE1 || defined AGE2
					//BX   do nothing here for tracer advection routine test.
# else
					tr[mAGE][k][i][j] = misval;
# endif
#endif
				}
			}
		}
	}


#ifdef AGE
	for (i = 0; i <= NXMEM - 1; i++) {
		for (j = 0; j <= NYMEM - 1; j++) {
# if defined AGEEXP || defined AGE1 || defined AGE2
			//BX   do nothing here for tracer advection routine test.
# else
			if (D[i][j] > MINIMUM_DEPTH) {
				tr[mAGE][0][i][j] = 0.0;
				for (k = 1; k < NZ; k++) {
					tr[mAGE][k][i][j] += dt;
				}
			}
# endif
		}
	}
#endif


#ifdef MERGED_ML
	merge_ml_tr();
#endif




	/*-----------------------------------------
	 *
	 *     OVERWRITE RESULTS IN SPECIAL CASES
	 *
	 *-----------------------------------------*/




#ifdef AGEEXP
	for (i=0;i<=NXMEM-1;i++) {
		for (j=0;j<=NYMEM-1;j++) {
			//BX - reinstated by HF
			if (D[i][j]>MINIMUM_DEPTH) {
				tr[mAGE][5][i][j] = 1.;
				age[5][i][j] = 1.;
				//BX - reinstated by HF
			} else {
				tr[mAGE][5][i][j] = misval;
				age[5][i][j] = misval;
			}
		}
	}
#endif

	/*-----------------------------------------
	 *
	 *    Sponges and Re-Entrant Boundaries
	 *
	 *-----------------------------------------*/

#ifdef SPONGE

	apply_sponge(h, dt, ea, eb);

#endif

	/* undo effect of oxygen sponge on delo18 */
	/* reset tr[mo18] to maintain old delo18  */


	/* zonal, meridional re-entrance    */
	for (m = 0; m < NTR; m++) {
		for (k = 0; k < NZ; k++) {
			for (j = 0; j <= NYMEM - 1; j++) {
				tr[m][k][0][j] = tr[m][k][nx - 1][j];
				tr[m][k][1][j] = tr[m][k][nx][j];
				tr[m][k][nx + 1][j] = tr[m][k][2][j];
				tr[m][k][nx + 2][j] = tr[m][k][3][j];
			}
		}
		for (i = 2; i <= nx; i++) {
			ii = 363 - i;
			for (k = 0; k < NZ; k++) {
				tr[m][k][ii][ny + 1] = tr[m][k][i][ny];
				tr[m][k][ii][ny + 2] = tr[m][k][i][ny - 1];
			}
		}
	}

	/*-----------------------------------------
	 *
	 *     PRINT TO SCREEN #3
	 *
	 *-----------------------------------------*/


	pstage = 3;
	print_to_screen(pquant, pstage);

	/*-----------------------------------------
	 *
	 *   store new tr array into named arrays for output
	 *
	 *-----------------------------------------*/

	for (k = 0; k < NZ; k++) {
		for (i = 0; i <= NXMEM - 1; i++) {
			for (j = 0; j <= NYMEM - 1; j++) {
#ifdef AGE
# if defined AGEEXP || defined AGE1 || defined AGE2
				age[k][i][j] = tr[mAGE][k][i][j];
# else
				age[k][i][j] = tr[mAGE][k][i][j] / (365.0 * 86400.0); /*years*/
# endif
#endif
			}
		}
	}

} /* end time step */

/*---------------------------------------------------------------------
 *
 *     BEGIN SUBROUTINES
 *
 *---------------------------------------------------------------------*/

void print_to_screen(int pquant[10], int pstage) {
#ifdef UNUSED
	int i,j,k;
#endif
	int m;

	/*----------------------------------
	 *     check points within domain
	 *----------------------------------*/

#ifdef UNUSED
	//HF bad loop order, revise if used
	for (i=0;i<=NXMEM-1;i++) {
		for (j=0;j<=NYMEM-1;j++) {
			for (k=0;k<NZ;k++) {
				if (D[i][j]>MINIMUM_DEPTH && h[k][i][j]>1.e-9) {
					if (tr[mOXYGEN][k][i][j] < trwarn[mOXYGEN][1]) {
						if (pquant[3] == 9) {
							printf("WARNING: RAISING LOW [O2] TO %g \n", trwarn[mOXYGEN][1]);
							printf(" i,j,k,h,depth,tr[o2] = %i,%i,%i,%g,%g,%g \n",
									i,j,k,h[k][i][j],depth[k][i][j],tr[mOXYGEN][k][i][j]);
						}
						tr[mOXYGEN][k][i][j] = trwarn[mOXYGEN][1];
						tr[mo18][k][i][j] = trwarn[mOXYGEN][1]*Rsmow*(1.e0+24/1.e3);
						Ro18o16[k][i][j] = tr[mo18][k][i][j] / tr[mOXYGEN][k][i][j];
						delo18[k][i][j] = (Ro18o16[k][i][j]/Rsmow - 1.e0) * 1.e3;
					}
					if (delo18[k][i][j] > trwarn[mo18][2] || delo18[k][i][j] < trwarn[mo18][1]) {
						if (pquant[3] == 9) {
							printf("step.c WARNING: EXTREME d18O VALUES \n");
							printf("i,j,k,D,o2,o18,Ro18o16,delo18 = %i,%i,%i,%g,%g,%g,%g,%g \n",
									i,j,k,D[i][j],tr[mOXYGEN][k][i][j],tr[mo18][k][i][j],
									Ro18o16[k][i][j],delo18[k][i][j]);
						}
					}
				}
			}
		}
	}

	/*----------------------------------
	 *     check points outside domain
	 *----------------------------------*/

	/*
	 for (i=0;i<=NXMEM-1;i++) {
	 for (j=0;j<=NYMEM-1;j++) {
	 for (k=0;k<NZ;k++) {
	 for (m=0;m<NTR;m++) {
	 if (h[k][i][j]<1.e-9 && tr[m][k][i][j]<1e-15) {
	 tr[m][k][i][j]=1e-15;
	 }
	 }
	 }
	 }
	 }
	 */


#endif /* UNUSED */

	for (m = 0; m < NTR; m++) {
		//BX	tracer_integral(m);
	}

}

void isotope_R_del(int trnum1, int trnum2, double Rstandard,
		double R[NZ][NXMEM][NYMEM], double del[NZ][NXMEM][NYMEM]) {

	int i, j, k;

	for (i = 0; i <= NXMEM - 1; i++) {
		for (j = 0; j <= NYMEM - 1; j++) {
			for (k = 0; k < NZ; k++) {
				//BX - reinstated by HF
				if (D[i][j] > MINIMUM_DEPTH) {
					R[k][i][j] = tr[trnum2][k][i][j] / tr[trnum1][k][i][j];
					del[k][i][j] = (R[k][i][j] / Rstandard - 1.e0) * 1.e3;
					//BX - reinstated by HF
				} else {
					R[k][i][j] = misval;
					del[k][i][j] = misval;
				}
			}
		}
	}

	/*  printf("isotope_R_del: o2,o18,Rsmow,del = %g,%g,%g,%g \n",
	 tr[trnum1][0][100][40],tr[trnum2][0][100][40],Rstandard,
	 del[0][100][40]);
	 */
}

void tracer_integral(int trnum, double ***hvol) {

	int i, j, k;

	trintegral[trnum] = 0.e0;
	for (k = 0; k < NZ; k++) {
		for (i = 2; i <= NXMEM - 1; i++) {
			for (j = 2; j <= NYMEM - 1; j++) {
				if (hvol[k][i][j] > 1e-6) {
					trintegral[trnum] += tr[trnum][k][i][j] * hvol[k][i][j];
				}
			}
		}
	}

}

void tracer_profile(double tr3d[NZ][NXMEM][NYMEM]) {

	int i, j, k;
	double layerarea;

	for (k = 0; k < NZ; k++) {
		trprofile[k] = 0.0;
		layerarea = 0.0;
		for (i = 2; i <= NXMEM - 1; i++) {
			for (j = 2; j <= NYMEM - 1; j++) {
				if (h[k][i][j] > 1e-6) {
					trprofile[k] += tr3d[k][i][j] * dxdyh[i][j];
					layerarea += dxdyh[i][j];
				}
			}
		}
		trprofile[k] = trprofile[k] / layerarea;
	}

}

void calc_ksub(double htemp[NZ][NXMEM][NYMEM]) {

	int i, j, kjunk;
	double dcrit;

	dcrit = 1.e0;

	for (i = X1; i <= nx; i++) {
		for (j = Y1; j <= ny; j++) {
			kjunk = 2;
			ksub[i][j] = kjunk;
			while (htemp[kjunk][i][j] < dcrit && kjunk < NZ) {
				ksub[i][j] = kjunk + 1;
				kjunk = kjunk + 1;
			}
		}
	}

}

void oxygen_saturation(double T[NZ][NXMEM][NYMEM], double S[NZ][NXMEM][NYMEM],
		double o2_sat[NZ][NXMEM][NYMEM]) {
	/* compute oxygen saturation concentrations in the mixed layer in ml/l
	 * using Garcia and Gordon 1992 (L&O,1307) and convert to mol/m^3.
	 * Based on Matlab code by Allison Shaw and Steve Emerson.
	 * Ivan Lima (Nov 2002) */

	int i, j, k;
	double TS[NZ][NXMEM][NYMEM], o2s_ml[NZ][NXMEM][NYMEM];
	double work[NZ][NXMEM][NYMEM]; /* temporary work array */
	double molvolo2 = 22.3916; /* molar volumes (liters/mole) */
	double A0, A1, A2, A3, A4, A5, B0, B1, B2, B3, C0; /* constants */

	A0 = 2.00907;
	A1 = 3.22014;
	A2 = 4.05010;
	A3 = 4.94457;
	A4 = -.256847;
	A5 = 3.88767;
	B0 = -.00624523;
	B1 = -.00737614;
	B2 = -.0103410;
	B3 = -.00817083;
	C0 = -4.88682E-07;

	for (k = 0; k < NZ; k++) {
		for (i = X1; i <= nx; i++) {
			for (j = Y1; j <= ny; j++) {
				// HF this was set to 298.15, but T is in C, not K here!
				// HF the Garcia and Gordon paper uses t=40C as upper bound
				if (T[k][i][j] > 40.0) {
					printf(
							"WARNING: TEMPERATURE OUT OF RANGE in oxygen saturation: %g \n",
							T[k][i][j]);
					printf(
							"  setting Temp to 40 deg C at array position %i,%i,%i \n",
							i, j, k);
					T[k][i][j] = 40.0;
				}
				TS[k][i][j]
						= log((298.15 - T[k][i][j]) / (273.15 + T[k][i][j]));
				work[k][i][j] = A0 + A1 * TS[k][i][j] + A2 * pow(TS[k][i][j],
						2.) + A3 * pow(TS[k][i][j], 3.) + A4 * pow(TS[k][i][j],
						4.) + A5 * pow(TS[k][i][j], 5.) + S[k][i][j] * (B0 + B1
						* TS[k][i][j] + B2 * pow(TS[k][i][j], 2.) + B3 * pow(
						TS[k][i][j], 3.)) + C0 * pow(S[k][i][j], 2.);
				o2s_ml[k][i][j] = exp(work[k][i][j]);
				o2_sat[k][i][j] = o2s_ml[k][i][j] / molvolo2;
			}
		}
	}
}

void argon_saturation(double T[NZ][NXMEM][NYMEM], double S[NZ][NXMEM][NYMEM],
		double ar_sat[NZ][NXMEM][NYMEM]) {
	/* compute argon saturation concentrations in the mixed layer in umol/kg
	 * using Roberta Hamme's thesis
	 * Based on Matlab code by Allison Shaw and Steve Emerson.
	 * Curtis Deutsch, 01/2004 */

	int i, j, k;
	double TS[NZ][NXMEM][NYMEM];
	double work[NZ][NXMEM][NYMEM]; /* temporary work array */
	double A0, A1, A2, A3, A4, A5, B0, B1, B2, B3, C0; /* constants */

	A0 = 2.79169;
	A1 = 3.17715;
	A2 = 4.13663;
	A3 = 4.86532;
	A4 = 0.e0;
	A5 = 0.e0;
	B0 = -6.96315e-3;
	B1 = -7.68391e-3;
	B2 = -1.19113e-2;
	B3 = 0.e0;
	C0 = 0.e0;

	for (k = 0; k < NZ; k++) {
		for (i = X1; i <= nx; i++) {
			for (j = Y1; j <= ny; j++) {
				TS[k][i][j]
						= log((298.15 - T[k][i][j]) / (273.15 + T[k][i][j]));
				work[k][i][j] = A0 + A1 * TS[k][i][j] + A2 * pow(TS[k][i][j],
						2.) + A3 * pow(TS[k][i][j], 3.) + A4 * pow(TS[k][i][j],
						4.) + A5 * pow(TS[k][i][j], 5.) + S[k][i][j] * (B0 + B1
						* TS[k][i][j] + B2 * pow(TS[k][i][j], 2.) + B3 * pow(
						TS[k][i][j], 3.)) + C0 * pow(S[k][i][j], 2.);
				ar_sat[k][i][j] = exp(work[k][i][j]);
			}
		}
	}

}

void co2_saturation(double T[NZ][NXMEM][NYMEM], double S[NZ][NXMEM][NYMEM],
		double co2_sat[NZ][NXMEM][NYMEM]) {

	/* compute co2 saturation concentrations in the mixed layer in umol/kg */

	int i, j, k;
	double TS[NZ][NXMEM][NYMEM];
	double work[NZ][NXMEM][NYMEM]; /* temporary work array */
	double A1, A2, A3, B1, B2, B3; /* constants */

	A1 = -58.0931e0;
	A2 = 90.5069e0;
	A3 = 22.2940e0;
	B1 = 0.027766e0;
	B2 = -0.025888e0;
	B3 = 0.0050578e0;

	for (k = 0; k < NZ; k++) {
		for (i = X1; i <= nx; i++) {
			for (j = Y1; j <= ny; j++) {
				TS[k][i][j] = (273.15 + T[k][i][j]) / 100.e0;
				work[k][i][j] = A1 + A2 / TS[k][i][j] + A3 * log(TS[k][i][j])
						+ S[k][i][j] * (B1 + B2 * TS[k][i][j] + B3
								* TS[k][i][j] * TS[k][i][j]);
				co2_sat[k][i][j] = exp(work[k][i][j]);
			}
		}
	}

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
	for (i = X1; i <= nx; i++) {
		for (j = Y1; j <= ny; j++) {
			//BX - reinstated by HF
			if (D[i][j] > MINIMUM_DEPTH) {
				hsum = h[0][i][j];
				depth[0][i][j] = h[0][i][j] * 0.5;
				for (k = 1; k < NZ; k++) {
					depth[k][i][j] = hsum + h[k][i][j] * 0.5;
					hsum += h[k][i][j];
				}
				//BX - reinstated by HF
			} else {
				for (k = 0; k < NZ; k++) {
					depth[k][i][j] = misval;
				}
			}
		}
	}

}

#ifdef USE_CALC_H
void z_sum(double h[NZ][NXMEM][NYMEM],double tot_depth[NXMEM][NYMEM])
{
	/* compute total depth in meters; derived from z_depth
	 * H. Frenzel - Oct 2009 */
	int i, j, k;
	double hsum;

	for (i=X1;i<=nx;i++) {
		for (j=Y1;j<=ny;j++) {
			tot_depth[i][j] = h[0][i][j];
			for (k=1;k<NZ;k++)
			tot_depth[i][j] += h[k][i][j];
		}
	}
}
#endif /* USE_CALC_H */

double lin_interp(double pleth, const double x[], const double z[], int istart,
		int npts) {
	/* finds first place where the value--PLETH--is crossed in the Z-array
	 * and returns X-array value corresponding to that level, KOUNT (=0 for
	 * no value found) and K (= array element #).
	 * Ivan Lima productions - Oct 2002 */
	double val1, val2;
	double xout; // NEW = x[npts];
	int kount;
	int k;

	kount = 0;
	for (k = istart + 1; k < npts; ++k) {
		val1 = pleth - z[k - 1];
		val2 = pleth - z[k];
		if (val1 < 0.0 && val2 < 0.0)
			goto end;
		if (val1 > 0.0 && val2 > 0.0)
			goto end;
		if (val1 == val2)
			goto end;
		kount = 1;
		xout = x[k - 1] + (x[k] - x[k - 1]) * val1 / (z[k] - z[k - 1]);
		break;
		end: ;
	}
	if (kount == 0) {
		xout = x[npts - 1]; //BX  x[npts] is not defined
		//BX    xout = x[npts];
		//BX    k = npts;
	}

	return xout;
}

#ifdef SMOOTH_REST
double lin_interpp(double pleth, const double x[], const double z[],
		int istart, int npts)
{
	/* finds first place where the value--PLETH--is crossed in the Z-array
	 * and returns X-array value corresponding to that level, KOUNT (=0 for
	 * no value found) and K (= array element #).
	 * Ivan Lima productions - Oct 2002 */
	double val1, val2;
	double xout; // NEW = x[npts];
	int kount;
	int k;

	// HF
	if (pleth < z[istart])
	return x[0];
	else if (pleth > z[npts-1])
	return x[npts-1];

	kount=0;
	for(k=istart+1;k<npts;++k) {
		val1=pleth-z[k-1];
		val2=pleth-z[k];
		//HF if(val1<0.0 && val2<0.0) {kount=2;break;}
		if(val1>0.0 && val2>0.0) continue;
		//HF if(val1 == val2) continue;
		kount=1;
		xout=x[k-1]+(x[k]-x[k-1])*val1/(z[k]-z[k-1]);
		break;
	}
	if (kount==0) {
		xout = x[npts-1]; //BX  x[npts] is not defined
		//BX    xout = x[npts];
		//BX    k = npts;
	}
	if (kount==2) {
		xout = x[0]; //BX  if shallower than first value, use first value
		//BX    xout = x[npts];
		//BX    k = npts;
	}

	return xout;
}
#endif

/* This routine maps observed concentrations onto isopycnal layers */
void conc_obs_layer(double h[NZ][NXMEM][NYMEM],
		double conc_lev[NZPHOS][NXMEM][NYMEM],
		double conc_lay[NZ][NXMEM][NYMEM]) {
	int i, j, k;
	int nzlevitus = NZPHOS;
	double depth[NZ][NXMEM][NYMEM];
	double concobsprof[NZPHOS];
	double levitus_depths[NZPHOS] = { 0, 10, 20, 30, 50, 75, 100, 120, 150,
			200, 250, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200,
			1300, 1400, 1500, 1750, 2000, 2500, 3000, 3500, 4000, 4500, 5000,
			5500 };

	/*---------------------------------------------------------------------
	 *     calculate vertically interpolated concentration:
	 *            levitus levels -> model layers
	 *             conc_lev  -> conc_lay
	 *---------------------------------------------------------------------*/
	z_depth(h, depth);

	//printf("Initialize 2a for month \n");
	//HF for (i=X1;i<=nx;i++) {
	//HF for (j=Y1;j<=ny;j++) {
	for (i = 0; i <= NXMEM - 1; i++) {
		for (j = 0; j <= NYMEM - 1; j++) {
			//BX - reinstated by HF
			if (D[i][j] > MINIMUM_DEPTH) {
				for (k = 0; k < nzlevitus; k++) {
					concobsprof[k] = conc_lev[k][i][j];
				}
				for (k = 0; k < NZ; k++) {
					conc_lay[k][i][j] = lin_interp(depth[k][i][j], concobsprof,
							levitus_depths, 0, nzlevitus);
					if (conc_lay[k][i][j] < 0.e0)
						conc_lay[k][i][j] = 0.;
				}
#ifdef MERGED_ML
				//BX overwrite the top layers
				for (k = 0; k <= 2; k = k + 2) {
					conc_lay[k][i][j] = lin_interp(((depth[k][i][j] + depth[k
							+ 1][i][j]) / 2.), concobsprof, levitus_depths, 0,
							nzlevitus);
					if (conc_lay[k][i][j] < 0.e0)
						conc_lay[k][i][j] = 0.;
					conc_lay[k + 1][i][j] = conc_lay[k][i][j];
				}
#endif

				//BX - reinstated by HF
			} else {
				for (k = 0; k < NZ; k++) {
					conc_lay[k][i][j] = misval;
				}
			}
		}
	}
	/*-----------------------------------------------------------------------
	 *     end of subroutine conc_obs_layer
	 *-----------------------------------------------------------------------*/
}

double drtsafe(void(*funcd)(double, double *, double *), double x1, double x2,
		double xacc) {

	//  void nrerror(char error_text[]);

	int j;
	double df, dx, dxold, f, fh, fl;
	double temp, xh, xl, rts;
	int MAXIT;

	MAXIT = 100;

	(*funcd)(x1, &fl, &df);
	(*funcd)(x2, &fh, &df);

	//  if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0))
	//  nrerror("Root must be bracketed in rtsafe");
	//if (fl == 0.0) return x1;
	//if (fh == 0.0) return x2;

	if (fl < 0.0) {
		xl = x1;
		xh = x2;
	} else {
		xh = x1;
		xl = x2;
	}
	rts = 0.5 * (x1 + x2);
	dxold = fabs(x2 - x1);
	dx = dxold;
	(*funcd)(rts, &f, &df);
	for (j = 1; j <= MAXIT; j++) {
		if ((((rts - xh) * df - f) * ((rts - xl) * df - f) > 0.0) || (fabs(2.0
				* f) > fabs(dxold * df))) {
			dxold = dx;
			dx = 0.5 * (xh - xl);
			rts = xl + dx;
			if (xl == rts)
				return rts;
		} else {
			dxold = dx;
			dx = f / df;
			temp = rts;
			rts -= dx;
			if (temp == rts)
				return rts;
		}
		if (fabs(dx) < xacc)
			return rts;
		(*funcd)(rts, &f, &df);
		if (f < 0.0)
			xl = rts;
		else
			xh = rts;
	}
	//  nrerror("Maximum number of iterations exceeded in rtsafe");
	return 0.0;
}

void ta_iter_1(double x, double *fn, double *df)

{

	//
	// This routine expresses TA as a function of DIC, htotal and constants.
	// It also calculates the derivative of this function with respect to
	// htotal. It is used in the iterative solution for htotal. In the call
	// "x" is the input value for htotal, "fn" is the calculated value for TA
	// and "df" is the value for dTA/dhtotal
	//

	double x2, x3, k12, k12p, k123p, c, a, a2, da, b, b2, db;
	extern double k1, k2, kw, kb, ks, kf, k1p, k2p, k3p, ksi;
	extern double bt, st, ft, sit, pt, dic1, ta;

	double x_inv, x2_inv;
	double a_inv, a2_inv;
	double b_inv, b2_inv;
	double c_inv;
	double kb_inv;
	double ksi_inv;

	x2 = x * x;
	x3 = x2 * x;
	k12 = k1 * k2;
	k12p = k1p * k2p;
	k123p = k12p * k3p;
	c = 1.0 + st / ks;
	a = x3 + k1p * x2 + k12p * x + k123p;
	a2 = a * a;
	da = 3.0 * x2 + 2.0 * k1p * x + k12p;
	b = x2 + k1 * x + k12;
	b2 = b * b;
	db = 2.0 * x + k1;

	//HF
	x_inv = 1.0 / x;
	x2_inv = x_inv * x_inv;
	a_inv = 1.0 / a;
	a2_inv = a_inv * a_inv;
	b_inv = 1.0 / b;
	b2_inv = b_inv * b_inv;
	c_inv = 1.0 / c;
	kb_inv = 1.0 / kb;
	ksi_inv = 1.0 / ksi;

	//
	// fn = hco3+co3+borate+oh+hpo4+2*po4+silicate+hfree+hso4+hf+h3po4-ta
	//

	*fn = k1 * x * dic1 * b_inv + 2.0 * dic1 * k12 * b_inv + bt / (1.0 + x
			* kb_inv) + kw * x_inv + pt * k12p * x * a_inv + 2.0 * pt * k123p
			* a_inv + sit / (1.0 + x * ksi_inv) - x * c_inv - st / (1.0 + ks
			* x_inv * c_inv) - ft / (1.0 + kf * x_inv) - pt * x3 * a_inv - ta;

	//
	//df = dfn/dx
	//

	*df = ((k1 * dic1 * b) - k1 * x * dic1 * db) * b2_inv - 2.0 * dic1 * k12
			* db * b2_inv - bt * kb_inv / pow(1.0 + x * kb_inv, 2) - kw
			* x2_inv + (pt * k12p * (a - x * da)) * a2_inv - 2.0 * pt * k123p
			* da * a2_inv - sit * ksi_inv / pow(1.0 + x * ksi_inv, 2) - 1.0
			* c_inv + st * pow(1.0 + ks * x_inv * c_inv, -2) * (ks * c_inv
			* x2_inv) + ft * pow(1.0 + kf * x_inv, -2) * kf * x2_inv - pt * x2
			* (3.0 * a - x * da) * a2_inv;

}

void merge_ml_tr() {
	//*************************************************************************
	// Merge the first and second and the third and fourth layer for all tracer
	// fields
	//*************************************************************************
	int i, j, k, m;
	for (j = Y1; j <= ny; j++) {
		for (i = X0; i <= nx + 1; i++) {
			for (m = 0; m < NTR; m++) {
				for (k = 0; k <= 2; k = k + 2) {
					//HF use z-weighted average
					tr[m][k][i][j] = (tr[m][k][i][j] * h[k][i][j]
							+ tr[m][k + 1][i][j] * h[k + 1][i][j])
							/ (h[k][i][j] + h[k + 1][i][j]);
					tr[m][k + 1][i][j] = tr[m][k][i][j];
				}
			}
		}
	}
}

void merge_ml_j() {
	//*************************************************************************
	// Merge the first and second and the third and fouth layer for all dic
	// related fields
	//*************************************************************************
	int i, j, k, m;
	for (i = X1; i <= nx; i++) {
		for (j = Y1; j <= ny; j++) {
			for (k = 0; k <= 2; k = k + 2) {
			}
		}
	}
}
