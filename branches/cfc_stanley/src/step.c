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
#include "gas_tracer_type.h"
#include "gas_exchange.h"
/*---------------------------------------------------------------------
 *     define variables and subroutines
 *---------------------------------------------------------------------*/

void z_depth(double h[NZ][NXMEM][NYMEM], double depth[NZ][NXMEM][NYMEM]);
void calc_ksub(double h[NZ][NXMEM][NYMEM]);
void conc_obs_layer(double h[NZ][NXMEM][NYMEM],
		double conc_lev[NZPHOS][NXMEM][NYMEM],
		double conc_lay[NZ][NXMEM][NYMEM]);
void tracer_integral(int trnum, double ***hvol);
// void print_to_screen(int pquant[10], int pstage);
#ifdef MERGED_ML
void merge_ml_tr();
void merge_ml_j();
#endif 

extern double h[NZ][NXMEM][NYMEM];
extern double hstart[NZ][NXMEM][NYMEM];
extern double hend[NZ][NXMEM][NYMEM];
extern double depth[NZ][NXMEM][NYMEM];
extern double D[NXMEM][NYMEM];
extern double ***ea, ***eb;
extern double ****tr;
extern double dt;
extern double trintegral[NTR];
extern int ksub[NXMEM][NYMEM];
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

#ifdef CFC
extern Gas_transient_t cfc11_props,cfc12_props,sf6_props;
extern double ***mn_cfc11;
extern double ***mn_cfc12;
extern double ***mn_sf6;
extern double ***cfc11;
extern double ***cfc12;
extern double ***sf6;
extern double mn_cfc11_sat[NXMEM][NYMEM];
extern double mn_cfc12_sat[NXMEM][NYMEM];
extern double mn_sf6_sat[NXMEM][NYMEM];
#endif

extern double atmpres[NXMEM][NYMEM];
extern double windspeed[NXMEM][NYMEM];
extern double fice[NXMEM][NYMEM];
extern double Temptm[NZ][NXMEM][NYMEM];
extern double Salttm[NZ][NXMEM][NYMEM];


# ifdef SPONGE
extern double sp_age[NZ][NXMEM][NYMEM];
# endif




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

	double mldepth_start, mldepth_end;
	//DT
	double Tmptr[NXMEM][NYMEM];
	int n;
	double secperyr = 365.25 * 60 * 60 * 24;
	int ndt_gasex = 100;
	double dt_gasex;
	extern double lath[NYMEM];
	//DT
#ifdef CFC
	double cfc11_nval,cfc12_nval,sf6_nval;
	double cfc11_sval,cfc12_sval,sf6_sval;
	double cfc11_atmconc,cfc12_atmconc,sf6_atmconc;
	double snval[2];
	double eqlat[2];
	double mon_time;
	extern double currtime;
	double Csurf, Csat;	
#endif


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
//	print_to_screen(pquant, pstage);


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
//	print_to_screen(pquant, pstage);

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
					tr[mAGE][k][i][j] = misval;
#endif
				}
			}
		}
	}


#ifdef AGE
	for (i = 0; i <= NXMEM - 1; i++) {
		for (j = 0; j <= NYMEM - 1; j++) {
			if (D[i][j] > MINIMUM_DEPTH) {
				tr[mAGE][0][i][j] = 0.0;
				for (k = 1; k < NZ; k++) {
					tr[mAGE][k][i][j] += dt;
				}
			}
		}
	}
#endif

#ifdef CFC 
	mon_time = currtime - dt/(2*secperyr);
	printf("Performing gas exchange for time %f\n",mon_time);
	cfc11_nval = lin_interp(cfc11_props.timestamp,cfc11_props.nval,mon_time,cfc11_props.ntime);
	cfc12_nval = lin_interp(cfc12_props.timestamp,cfc12_props.nval,mon_time,cfc12_props.ntime);
	sf6_nval = lin_interp(sf6_props.timestamp,sf6_props.nval,mon_time,sf6_props.ntime);

	cfc11_sval = lin_interp(cfc11_props.timestamp,cfc11_props.sval,mon_time,cfc11_props.ntime);
	cfc12_sval = lin_interp(cfc12_props.timestamp,cfc12_props.sval,mon_time,cfc12_props.ntime);
	sf6_sval = lin_interp(sf6_props.timestamp,sf6_props.sval,mon_time,sf6_props.ntime);

	dt_gasex = dt/ndt_gasex;
	copy_from_main_tracer_array(cfc11,cfc11_props.tridx,tr);
	copy_from_main_tracer_array(cfc12,cfc12_props.tridx,tr);
	copy_from_main_tracer_array(sf6,sf6_props.tridx,tr);

	for (i = 0; i < NXMEM; i++)
		for (j = 0; j < NYMEM; j++) {
			if (D[i][j] > MINIMUM_DEPTH) {// Perform gas exchange only on ocean points

				mldepth_start = hstart[0][i][j]+hstart[1][i][j];
				mldepth_end = hend[0][i][j]+hend[1][i][j];

				Csurf = cfc11[0][i][j];
				cfc11_atmconc = choose_atmconc_val(cfc11_sval,cfc11_nval,lath[j]);
				if(cfc11_atmconc > 0) {
				gas_exchange( i, j, cfc11_props, 
						cfc11_atmconc, mldepth_start, mldepth_end, dt_gasex,
						ndt_gasex, &cfc11[0][i][j], &Csat );
				mn_cfc11_sat[i][j] = Csat;
				}
				else cfc11[0][i][j]=0.0;

				cfc12_atmconc = choose_atmconc_val(cfc12_sval,cfc12_nval,lath[j]);
				if(cfc12_atmconc > 0) {
				gas_exchange( i, j, cfc12_props,
						cfc12_atmconc, mldepth_start, mldepth_end, dt_gasex,
						ndt_gasex, &cfc12[0][i][j], &Csat );
				mn_cfc12_sat[i][j] = Csat;
				}
				else cfc12[0][i][j]=0.0;

				Csurf = sf6[0][i][j];
				sf6_atmconc = choose_atmconc_val(sf6_sval,sf6_nval,lath[j]);
				if (sf6_atmconc > 0) {
				gas_exchange( i, j, sf6_props,
						sf6_atmconc, mldepth_start, mldepth_end, dt_gasex,
						ndt_gasex, &sf6[0][i][j],&Csat );
				sf6[0][i][j] = Csurf;
				mn_sf6_sat[i][j] = Csat;
				}
				else sf6[0][i][j]=0.0;
			}
			else {
				cfc11[0][i][j]=0.0;
				cfc12[0][i][j]=0.0;
				sf6[0][i][j]=0.0;
			}
			cfc11[1][i][j]=cfc11[0][i][j];
			cfc12[1][i][j]=cfc12[0][i][j];
			sf6[1][i][j]=sf6[0][i][j];
		}


	printf("CFC11: %f CFC12: %f SF6: %f\n",cfc11[0][100][100],cfc12[0][100][100],sf6[0][100][100]);
	printf("Northern Hemisphere: CFC11: %f CFC12: %f SF6: %f\n",cfc11_nval,cfc12_nval,sf6_nval);
	printf("Southern Hemisphere: CFC11: %f CFC12: %f SF6: %f\n",cfc11_sval,cfc12_sval,sf6_sval);
	printf("Temp: %f Salt: %f MLD: %f\n",Temptm[0][100][100],Salttm[0][100][100],hstart[0][100][100]);
	printf("wind: %f atmpres: %f fice: %f\n",windspeed[100][100],atmpres[100][100],fice[100][100]);

	copy_to_main_tracer_array(cfc11,cfc11_props.tridx,tr);
	copy_to_main_tracer_array(cfc12,cfc12_props.tridx,tr);
	copy_to_main_tracer_array(sf6,sf6_props.tridx,tr);

#endif

#ifdef MERGED_ML
	merge_ml_tr();
#endif




	/*-----------------------------------------
	 *
	 *     OVERWRITE RESULTS IN SPECIAL CASES
	 *
	 *-----------------------------------------*/

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
//	print_to_screen(pquant, pstage);

	/*-----------------------------------------
	 *
	 *   store new tr array into named arrays for output
	 *
	 *-----------------------------------------*/

	for (k = 0; k < NZ; k++) {
		for (i = 0; i <= NXMEM - 1; i++) {
			for (j = 0; j <= NYMEM - 1; j++) {
#ifdef AGE
				age[k][i][j] = tr[mAGE][k][i][j] / (365.0 * 86400.0); /*years*/
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
