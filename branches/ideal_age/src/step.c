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
#include "step.h"
#include "tracadv.h"
#include "util.h"

#if defined(OXYGEN) && defined(PHOSPHATE)
#include "tracer_utilities.h"
#include "oxygen.h"
#include "phosphate.h"
#endif

#ifdef AGE
#include "ideal_age.h"
#endif
/*---------------------------------------------------------------------
 *     define variables and subroutines
 *---------------------------------------------------------------------*/

void calc_ksub(double h[NZ][NXMEM][NYMEM]);

// void print_to_screen(int pquant[10], int pstage);
#ifdef MERGED_ML
void merge_ml_tr();
void merge_ml_j();
#endif 
//double cfc_atmospheric(double now, int icfc);

extern double h[NZ][NXMEM][NYMEM];
extern double hstart[NZ][NXMEM][NYMEM];
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
extern int mAGE;
# ifdef SPONGE
extern double sp_age[NZ][NXMEM][NYMEM];
# endif
#endif



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

	//DT
	double Tmptr[NXMEM][NYMEM];
	int n;
	double secperyr = 365.25 * 60 * 60 * 24;
	//DT



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

	// z_depth(h, depth);


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
#if defined(PHOSPHATE) && defined(OXYGEN)
	ibiodt = 1; // number of biotic time steps per transport time step (dt)

	// Update phosphate concentration for surface value restoration
	read_woa_file(imon, hstart, po4_star_lay, "woa09.phos.nc", "p_an",1e-3);
/*
	set_fix_darray3d_zero(jo2,NZ);
	set_fix_darray3d_zero(jdop,NZ);
	set_fix_darray3d_zero(jpo4,NZ);
	set_fix_darray3d_zero(jprod,NZ);
	set_fix_darray3d_zero(jremin,NZ);
	set_fix_darray3d_zero(jremdop,NZ);
	set_fix_darray2d_zero(flux_pop);
*/	
	printf("Calculating biotic sources/sinks\n");
	biotic_sms(ibiodt);

	merge_ml_tr();
	merge_ml_j();
	apply_oxygen_jterms();
	apply_phosphate_jterms();
	surface_oxygen();


#endif




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
				}
			}
		}
	}


#ifdef AGE
	step_age(dt);
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
				//for (k = 0; k < 1; k = k + 2) { // ashao: Buffer layers do not need to have uniform concentration
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
			for (k = 0; k <= 2; k = k + 2) { //ashao: buffer layer do not need to be uniform
//			for (k = 0; k < 1; k = k + 2) { //ashao: buffer layer do not need to be uniform

# ifdef PHOSPHATE 
                                jpo4[k][i][j] = (jpo4[k][i][j]*h[k][i][j]+
                                                jpo4[k+1][i][j]*h[k+1][i][j])/
                                (h[k][i][j] + h[k+1][i][j]);
                                jpo4[k+1][i][j] = jpo4[k][i][j];

                                jdop[k][i][j] = (jdop[k][i][j]*h[k][i][j]+
                                                jdop[k+1][i][j]*h[k+1][i][j])/
                                (h[k][i][j] + h[k+1][i][j]);
                                jdop[k+1][i][j] = jdop[k][i][j];

                                jremdop[k][i][j] = (jremdop[k][i][j]*h[k][i][j]+
                                                jremdop[k+1][i][j]*h[k+1][i][j])/
                                (h[k][i][j] + h[k+1][i][j]);
                                jremdop[k+1][i][j]= jremdop[k][i][j];

                                jprod[k][i][j] = (jprod[k][i][j]*h[k][i][j]+
                                                jprod[k+1][i][j]*h[k+1][i][j])/
                                (h[k][i][j] + h[k+1][i][j]);
                                jprod[k+1][i][j] = jprod[k][i][j];

                                jremin[k][i][j] = (jremin[k][i][j]*h[k][i][j]+
                                                jremin[k+1][i][j]*h[k+1][i][j])/
                                (h[k][i][j] + h[k+1][i][j]);
                                jremin[k+1][i][j] = jremin[k][i][j];
# endif // PHOSPHATE
# ifdef OXYGEN
                                jo2[k][i][j] = (jo2[k][i][j]*h[k][i][j]+
                                                jo2[k+1][i][j]*h[k+1][i][j])/
                                (h[k][i][j] + h[k+1][i][j]);
                                jo2[k+1][i][j] = jo2[k][i][j];
# endif
			}
		}
	}
}
