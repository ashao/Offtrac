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

double calc_tracer_inventory( int idx );

extern double hend[NZ][NXMEM][NYMEM];
extern double h[NZ][NXMEM][NYMEM];
extern double h_clim[NZ][NXMEM][NYMEM];
extern double depth[NZ][NXMEM][NYMEM];
extern double rml[2][NXMEM][NYMEM];
extern double D[NXMEM][NYMEM];
extern int wetmask[NXMEM][NYMEM];
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

#ifdef INERT
extern double inert1[NZ][NXMEM][NYMEM];
extern double inert2[NZ][NXMEM][NYMEM];
extern double ar_sat[NZ][NXMEM][NYMEM];
extern int mINERT1;
extern int mINERT2;
# ifdef SPONGE
extern double sp_inert1[NZ][NXMEM][NYMEM];
extern double sp_inert2[NZ][NXMEM][NYMEM];
# endif
#endif
#ifdef AGE
extern double ***age;
extern int mAGE;
# ifdef SPONGE
extern double sp_age[NZ][NXMEM][NYMEM];
# endif
#endif

#ifdef OXYGEN
extern double ***oxy_init;
extern double o2_sat[NZ][NXMEM][NYMEM];
extern double oxyflux[NXMEM][NYMEM];
extern double jo2[NZ][NXMEM][NYMEM];
extern double jo2_spec[NZ][NXMEM][NYMEM];
extern double aou_spec[NZ][NXMEM][NYMEM];
extern int mOXYGEN;
double jo2_ocmip[NZ][NXMEM][NYMEM];
double Sc_o2[NXMEM][NYMEM], Sk_o2[NXMEM][NYMEM];
double ratio_sw,alpha_resp,ratio_o2_sw,Rsmow;
# ifdef SPONGE
extern double sp_oxygen[NZ][NXMEM][NYMEM];
# endif
#endif

#ifdef OXY18
extern double o18[NZ][NXMEM][NYMEM];
extern double o18flux[NXMEM][NYMEM];
extern double jo18[NZ][NXMEM][NYMEM];
extern double delo18[NZ][NXMEM][NYMEM];
extern double Ro18o16[NZ][NXMEM][NYMEM];
extern int mo18;
# ifdef SPONGE
extern double sp_oxygen[NZ][NXMEM][NYMEM];
# endif
#endif

#ifdef CFC
extern double cfc11[NZ][NXMEM][NYMEM];
extern double cfc12[NZ][NXMEM][NYMEM];
extern double c11flux[NXMEM][NYMEM];
extern double c12flux[NXMEM][NYMEM];
extern double c11monflux[NXMEM][NYMEM];
extern double c12monflux[NXMEM][NYMEM];
extern double c11res[NZ][NXMEM][NYMEM];
extern double c12res[NZ][NXMEM][NYMEM];
extern int mCFC11, mCFC12;
/*    Begin added DT     */
extern double Sc_cfc11[NXMEM][NYMEM];
extern double Sc_cfc12[NXMEM][NYMEM];
extern double Sk_cfc11[NXMEM][NYMEM];
extern double Sk_cfc12[NXMEM][NYMEM];
extern double cfc11_sol[NZ][NXMEM][NYMEM];
extern double cfc12_sol[NZ][NXMEM][NYMEM];
extern double cfc11_sat[NXMEM][NYMEM];
extern double cfc12_sat[NXMEM][NYMEM];
extern double cfcatm11[NXMEM][NYMEM];
extern double cfcatm12[NXMEM][NYMEM];

// begin ashao
extern int cfc11nrecs;
extern int cfc11year[MAXRECS];
extern double cfc11N[MAXRECS], cfc11S[MAXRECS];
extern double cfc11sol_coeffs[7], cfc11sc_coeffs[4];
extern int cfc12nrecs;
extern int cfc12year[MAXRECS];
extern double cfc12N[MAXRECS], cfc12S[MAXRECS];
extern double cfc12sol_coeffs[7], cfc12sc_coeffs[4];
extern double mn_cfc11_relsat[NZ][NXMEM][NYMEM], mn_cfc12_relsat[NZ][NXMEM][NYMEM];
extern int mCFC11_sat, mCFC12_sat;
double cfc11_relsat, cfc12_relsat;
// end ashao
/*    End added DT       */

#endif

# ifdef SPONGE
extern double sp_age[NZ][NXMEM][NYMEM];
# endif

// begin ashao
#ifdef SF6
extern double sf6[NZ][NXMEM][NYMEM];
extern double sf6flux[NXMEM][NYMEM];
extern double sf6monflux[NXMEM][NYMEM];
extern double sf6res[NZ][NXMEM][NYMEM];
extern int mSF6;
extern double Sc_sf6[NXMEM][NYMEM];
extern double Sk_sf6[NXMEM][NYMEM];
extern double sf6_sol[NZ][NXMEM][NYMEM];
extern double sf6_sat[NXMEM][NYMEM];
extern double sf6atm[NXMEM][NYMEM];
extern int sf6nrecs;
extern int sf6year[MAXRECS];
extern double sf6N[MAXRECS], sf6S[MAXRECS];
extern double sf6sol_coeffs[7], sf6sc_coeffs[4];
extern double mn_sf6_relsat[NZ][NXMEM][NYMEM];
extern int mSF6_sat;
double sf6_relsat;
#endif

double dt_gasex;
double h_mixedlyr;
double TempK;
#if defined(SF6) || defined(CFC)

double cfc11tmp, cfc12tmp;
double yearfrac_gasex;
double sf6tmp;
#endif
#ifdef RADIOACTIVE
extern double i131[NZ][NXMEM][NYMEM];
extern double cs134[NZ][NXMEM][NYMEM];
extern double cs136[NZ][NXMEM][NYMEM];
extern double cs137[NZ][NXMEM][NYMEM];
extern double ba140[NZ][NXMEM][NYMEM];
extern double la140[NZ][NXMEM][NYMEM];

extern int mi131;
extern int mcs134;
extern int mcs136;
extern int mcs137;
extern int mba140;
extern int mla140;

#endif
#ifdef TTD_BP
extern double ttd[NZ][NXMEM][NYMEM];
extern double ttd_mask[NXMEM][NYMEM];
extern int mttd;
#endif
// end ashao

#ifdef PHOSPHATE
extern int mPHOSPHATE;
extern double ***phosphate_init;
extern double po4_star_lev[NZPHOS][NXMEM][NYMEM];
extern double po4_star_lay[NZ][NXMEM][NYMEM];
extern double jpo4[NZ][NXMEM][NYMEM];
extern int mDOP;
extern double jdop[NZ][NXMEM][NYMEM];
extern double jremdop[NZ][NXMEM][NYMEM];
extern double jprod[NZ][NXMEM][NYMEM];
extern double jremin[NZ][NXMEM][NYMEM];
extern double flux_pop[NXMEM][NYMEM];
# ifdef PROGNOSTIC
extern double ***lightlim;
#  ifdef NITRATE
extern double ***nitrlim;
#  endif
extern double ***ironlim;
# endif /* PROGNOSTIC */
# ifdef SPONGE
extern double sp_phosphate[NZ][NXMEM][NYMEM];
extern double sp_dop[NZ][NXMEM][NYMEM];
# endif
#endif /* PHOSPHATE */

#ifdef NITRATE
extern int mNITRATE;
extern double ***nitrate_init;
extern double no3_lev[NZPHOS][NXMEM][NYMEM];
extern double no3_lay[NZ][NXMEM][NYMEM];
extern double jno3[NZ][NXMEM][NYMEM];
extern double jn2fix[NZ][NXMEM][NYMEM];
extern double jdenits[NZ][NXMEM][NYMEM];
extern double jdenitw[NZ][NXMEM][NYMEM];
extern int mDON;
extern double jdon[NZ][NXMEM][NYMEM];
extern double jremdon[NZ][NXMEM][NYMEM];
# ifdef SPONGE
extern double sp_nitrate[NZ][NXMEM][NYMEM];
extern double sp_don[NZ][NXMEM][NYMEM];
# endif
# ifdef N15_CYCLE
extern int mNITRATE_15n;
extern double ***nitrate_init_15n;
extern double jno3_15n[NZ][NXMEM][NYMEM];
extern int mDON_15n;
extern double jdon_15n[NZ][NXMEM][NYMEM];
extern double jremdon_15n[NZ][NXMEM][NYMEM];
#  ifdef SPONGE
extern double sp_nitrate_15n[NZ][NXMEM][NYMEM];
extern double sp_don_15n[NZ][NXMEM][NYMEM];
#  endif
# endif /* N15_CYCLE */
#endif /* NITRATE */

#ifdef DIC
extern double dic[NZ][NXMEM][NYMEM];
extern double alk[NZ][NXMEM][NYMEM];
extern double deltaco2[NXMEM][NYMEM];
extern int mDIC;
extern int mALK;
extern double co2_sat[NZ][NXMEM][NYMEM];
extern double co2flux[NXMEM][NYMEM];
extern double ***jdic;
extern double ***jalk;
extern double k0,k1,k2,kb,kw,ks,kf,k1p,k2p,k3p,ksi;
extern double ff,htotal,htotal2,xco2,xacc,x1,x2;
extern double bt,st,ft,sit,pt,dic1,ta,co2starair;
extern double dco2star;
extern double pco2surf;
extern double pco2s[NXMEM][NYMEM];
extern double Sk_co2[NXMEM][NYMEM];
double Sc_co2[NXMEM][NYMEM];
double phlo,phhi;
double tmpdic[NXMEM][NYMEM];
double tmpalk[NXMEM][NYMEM];
double tmpco2flux;
# ifdef SPONGE
extern double sp_dic[NZ][NXMEM][NYMEM];
extern double sp_alk[NZ][NXMEM][NYMEM];
# endif
# ifdef REST_ARCTIC
extern double dic_lay[NZ][NXMEM][NYMEM];
extern double alk_lay[NZ][NXMEM][NYMEM];
# endif
#endif

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
#ifdef REST_ARCTIC
	extern double qlat[NYMEM];
	double wfunct[NYMEM];
#endif
	//DT
	double Tmptr[NXMEM][NYMEM];
	int n;
	double secperyr = 365.25 * 60 * 60 * 24;
	theyear = iyear;
	//DT
	double thetime;
#if defined(SF6) | defined(CFC)
	// begin ashao
	extern double lath[NXTOT];
	// end ashao
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

#ifdef OXYGEN
	trwarn[mOXYGEN][0] = 1.0e-15;
	trwarn[mOXYGEN][1] = 1.e0;
#endif
#ifdef OXY18
	trwarn[mo18][0] = 1.e0;
	trwarn[mo18][1] = 1.e2;
#endif
#ifdef PHOSPHATE
	trwarn[mPHOSPHATE][0] = 0.e0;
	trwarn[mPHOSPHATE][1] = 5.e0;
	trwarn[mDOP][0] = 0.e0;
	trwarn[mDOP][1] = 1.e0;
#endif
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

#ifdef OXY18
	isotope_R_del(mOXYGEN,mo18,Rsmow,Ro18o16,delo18);
#endif

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

#ifdef obsPsupply

	conc_obs_layer(h,po4_star_lev,po4_star_lay);
	conc_obs_layer(h,no3_lev,no3_lay);

	for (i=0;i<=NXMEM-1;i++) {
		for (j=0;j<=NYMEM-1;j++) {
			if (D[i][j] > MINIMUM_DEPTH) {
				for (k=0; k<NZ; k++)
				if (depth[k][i][j] > 200.) { // HF 09282009
					tr[mPHOSPHATE][k][i][j] = po4_star_lay[k][i][j];
					tr[mNITRATE][k][i][j] = no3_lay[k][i][j];
				}
			}
		}
	}

#endif

#ifdef TTD_BP

		if ( iterno < 12) {
			printf("Setting Mixed Layer to 1 in month %d\n",iterno);
	                for (k=0;k<=1;k++)
	                    for (i=0;i<NXMEM;i++)
	                	for (j=0;j<NYMEM;j++) {

	                	    if (ttd_mask[i][j]){
		                        ttd[k][i][j]=ttd_mask[i][j]/12;
		                        tr[mttd][k][i][j]=1.0/12;
	                	    }
	                	    else {
		                        ttd[k][i][j]=0.0;
		                        tr[mttd][k][i][j]=0.0;
	                	    }

	                }
		}

#endif

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

#ifdef OXY18
	isotope_R_del(mOXYGEN,mo18,Rsmow,Ro18o16,delo18);
#endif

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

#ifdef PHOSPHATE 

	//  set input terms of biotic_sms to original value, output terms to 0.
	set_fix_darray3d_zero(jpo4, NZ);
	set_fix_darray3d_zero(jdop, NZ);
	set_fix_darray3d_zero(jremdop, NZ);
	set_fix_darray3d_zero(jprod, NZ);
	set_fix_darray3d_zero(jremin, NZ);
# ifdef PROGNOSTIC
	set_darray3d_zero(lightlim, NZ, NXMEM, NYMEM);
#  ifdef NITRATE
	set_darray3d_zero(nitrlim, NZ, NXMEM, NYMEM);
#  endif
	set_darray3d_zero(ironlim, NZ, NXMEM, NYMEM);
# endif /* PROGNOSTIC */
# ifdef NITRATE
	set_fix_darray3d_zero(jno3, NZ);
	set_fix_darray3d_zero(jn2fix, NZ);
	set_fix_darray3d_zero(jdenits, NZ);
	set_fix_darray3d_zero(jdenitw, NZ);
	set_fix_darray3d_zero(jdon, NZ);
	set_fix_darray3d_zero(jremdon, NZ);
#  ifdef N15_CYCLE
	set_fix_darray3d_zero(jno3_15n, NZ);
	set_fix_darray3d_zero(jdon_15n, NZ);
	set_fix_darray3d_zero(jremdon_15n, NZ);
#  endif
# endif /* NITRATE */
# ifdef OXYGEN
	set_fix_darray3d_zero(jo2, NZ);
# endif
# ifdef DIC
	set_darray3d_zero(jdic, NZ, NXMEM, NYMEM);
	set_darray3d_zero(jalk, NZ, NXMEM, NYMEM);
# endif
	set_fix_darray2d_zero(flux_pop);
# ifndef NO_BIOTIC
	biotic_sms(ibiodt);
# endif
#endif  /*  PHOSPHATE  */

# ifdef  JO2_SPEC
	for (k=0;k<NZ;k++)
	for (i=X1;i<=nx;i++)
	for (j=Y1;j<=ny;j++)
	jo2[k][i][j] = jo2_spec[k][i][j];
# endif

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
			if (D[i][j]>MINIMUM_DEPTH) {
				for (k = 0; k < NZ; k++) {
#ifdef PHOSPHATE
					tr[mPHOSPHATE][k][i][j] += dt * jpo4[k][i][j];
					tr[mDOP][k][i][j] += dt * jdop[k][i][j];
//					if ( tr[mPHOSPHATE][k][i][j] < 0.0 ) tr[mPHOSPHATE][k][i][j] = 0.0;
//					if ( tr[mDOP][k][i][j] < 0.0 ) tr[mDOP][k][i][j] = 0.0;
#endif
#ifdef NITRATE
					tr[mNITRATE][k][i][j] += dt * jno3[k][i][j];
					tr[mDON][k][i][j] += dt * jdon[k][i][j];
# ifdef N15_CYCLE
					tr[mNITRATE_15n][k][i][j] += dt * jno3_15n[k][i][j];
					tr[mDON_15n][k][i][j] += dt * jdon_15n[k][i][j];
# endif
#endif /* NITRATE */
#ifdef OXYGEN
					tr[mOXYGEN][k][i][j] += dt * jo2[k][i][j];
//					if (tr[mOXYGEN][k][i][j] < 0.0) tr[mOXYGEN][k][i][j] = 0.0;
#endif
#ifdef OXY18
					tr[mo18][k][i][j] += dt * jo18[k][i][j];
#endif
				}
				//BX - reinstated by HF
			} else {
				for (k = 0; k < NZ; k++) {
#ifdef PHOSPHATE
					tr[mPHOSPHATE][k][i][j] = misval;
					tr[mDOP][k][i][j] = misval;
					jpo4[k][i][j] = misval;
#endif
#ifdef OXYGEN
					tr[mOXYGEN][k][i][j] = misval;
					jo2[k][i][j] = misval;
#endif
#ifdef OXY18
					tr[mo18][k][i][j] = misval;
#endif
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

#ifdef DIC
	for (i=0;i<=NXMEM-1;i++) {
		for (j=0;j<=NYMEM-1;j++) {
			if (D[i][j]>MINIMUM_DEPTH) {
				for (k=1;k<NZ;k++) {
					tr[mDIC][k][i][j] += dt * jdic[k][i][j];
					tr[mALK][k][i][j] += dt * jalk[k][i][j];
				}
			} else {
				for (k=0;k<NZ;k++) {
					tr[mDIC][k][i][j] = misval;
					tr[mALK][k][i][j] = misval;
				}
			}
		}
	}
#endif

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

#ifdef INERT

	/*    argon_saturation(Temptm, Salttm, ar_sat); */
	argon_saturation(inert2, Salttm, ar_sat);

	for (i=0;i<=NXMEM-1;i++) {
		for (j=0;j<=NYMEM-1;j++) {
			ar_sat[0][i][j]= ar_sat[0][i][j] * rml[0][i][j] * 1e-3;
			ar_sat[1][i][j]= ar_sat[1][i][j] * rml[1][i][j] * 1e-3;
		}
	}

	for (i=0;i<=NXMEM-1;i++) {
		for (j=0;j<=NYMEM-1;j++) {
			//BX - reinstated by HF
			if (D[i][j]>MINIMUM_DEPTH) {
				tr[mINERT1][0][i][j] = ar_sat[0][i][j]; // set boundary condition only
				tr[mINERT2][0][i][j] = Temptm[0][i][j];
				//BX - reinstated by HF
			} else {
				for (k=0;k<NZ;k++) {
					tr[mINERT2][k][i][j] = misval;
				}
			}
		}
	}

#endif

#ifdef MERGED_ML
	merge_ml_tr();
#endif

	// #define GASEX ashao: moved to init.h

#ifdef GASEX
	/*-----------------------------------------
	 *
	 *     STEP 3: AIR-SEA FLUXES
	 *
	 *-----------------------------------------*/

	/*---------------------------------------------------------------------
	 *
	 *  GAS EXCHANGE:
	 *
	 *  The gas exchange time step applies only to the DIC and Alk terms.
	 *
	 * a) Compute saturation concentration in mixed layer
	 * b) Compute the Schmidt number (Sc: unitless)
	 * c) Compute the piston velocity (Sk: cm/hr), convert to m/s
	 * d) Compute tracer flux into mixed layer (m/s)*(mol/m^3) -> mol/m2/s
	 * e) Compute the tracer concentration
	 *
	 *---------------------------------------------------------------------*/

	mdt = 30; // number of gasex time steps per transport time step (dt)
	dt_gasex = dt / ((double) mdt);
	inv_mdt = 1.0 / (double) mdt;

	printf("Calculate gas exchange terms with %i time step(s) per dt. \n\n",
			mdt);
	printf("Current Time %f\n",currtime);

	//begin ashao


#if defined(SF6) || defined(CFC)

#if !defined(RESTCONC) && !defined(ZEROENTRAIN)
	printf("Using concentrations after full 3D advection\n");
#endif
#ifdef RESTCONC
	printf("Restoring concentration in mixed layer to previous timestep\n");
#endif
#ifdef ZEROENTRAIN
	printf("Assuming entrainment of zero concentration water\n");
#endif
	printf("CFC-12 used in gas exchange: %e\n",tr[mCFC12][0][87][149]);

	// Calculate variables universal to all tracers:

	yearfrac_gasex=dt_gasex/31536000;

//        printf("Gas exchange h=%g\n",h[1][100][100]);

	for (i = X1; i < nx; i++) {
		for (j = Y1; j < ny; j++) {
			// Calculations for 2D arrays (LAT x LON)
			// Property values at grid location
			h_mixedlyr = h[0][i][j] + h[1][i][j]; // Thickness of mixed layer (no buffer layers)
			TempK = Temptm[0][i][j] + 273.15; // Temperature at current grid location
			// Calculate tracer solubility, Schmidt number, and initialize variables for gas exchange
#ifdef CFC
			cfc11_sol[0][i][j] = trac_calcsol(TempK, Salttm[0][i][j],
					cfc11sol_coeffs);
			Sc_cfc11[i][j] = trac_schmidt(Temptm[0][i][j], cfc11sc_coeffs);

			cfc12_sol[0][i][j] = trac_calcsol(TempK, Salttm[0][i][j],
					cfc12sol_coeffs);
			Sc_cfc12[i][j] = trac_schmidt(Temptm[0][i][j], cfc12sc_coeffs);

			c11flux[i][j] = 0.0;
			c12flux[i][j] = 0.0;
			c11monflux[i][j]=0.0;
			c12monflux[i][j]=0.0;

#if !defined(RESTCONC) && !defined(ZEROENTRAIN)
			cfc11tmp = tr[mCFC11][0][i][j];
			cfc12tmp = tr[mCFC12][0][i][j];
#endif //REGULAR
#if defined(RESTCONC) || defined(ZEROENTRAIN)
			cfc11tmp = cfc11[0][i][j];
			cfc12tmp = cfc12[0][i][j];
#endif //RESTCONC or ZEROENTRAIN

#endif // CFC
#ifdef SF6
			sf6_sol[0][i][j] = trac_calcsol(TempK, Salttm[0][i][j],
					sf6sol_coeffs);
			Sc_sf6[i][j] = trac_schmidt(Temptm[0][i][j], sf6sc_coeffs);
			sf6flux[i][j] = 0.0;
			sf6monflux[i][j]=0.0;
#if !defined(RESTCONC) && !defined(ZEROENTRAIN)
			sf6tmp = tr[mSF6][0][i][j];
#endif //REGULAR

#if defined(RESTCONC) || defined(ZEROENTRAIN)
			sf6tmp = sf6[0][i][j];
#endif //RESTCONCi or ZEROENTRAIN
#endif // SF6
			// Gas exchange time loop
			if (D[i][j] > MINIMUM_DEPTH) { // Calculate gas exchange for ocean tile
				// Store previous monthly value for use as tracer concentration

				for (m = 0; m < mdt; m++) {
					thetime=currtime+yearfrac_gasex*(m+1);
					// Calculate atmospheric concentration for tracers
#ifdef SF6
					sf6atm[i][j] = trac_atmospheric(thetime,
							sf6nrecs, sf6year, sf6N, sf6S, j);
#endif // SF6
#ifdef CFC
					cfcatm11[i][j] = trac_atmospheric(thetime,
							cfc11nrecs, cfc11year, cfc11N, cfc11S, j);
					cfcatm12[i][j] = trac_atmospheric(thetime,
							cfc12nrecs, cfc12year, cfc12N, cfc12S, j);
#endif

					// Calculate tracer flux and store to temporary tracer concentration
#ifdef SF6
					sf6flux[i][j] = trac_calcflux(Sc_sf6[i][j], sf6atm[i][j],
							sf6_sol[0][i][j], fice[i][j], xkw[i][j],
							atmpres[i][j], sf6tmp, dt_gasex);
					sf6tmp += sf6flux[i][j];
					sf6monflux[i][j]+=sf6flux[i][j];
#endif
#ifdef CFC
					c11flux[i][j] = trac_calcflux(Sc_cfc11[i][j],
							cfcatm11[i][j], cfc11_sol[0][i][j], fice[i][j],
							xkw[i][j], atmpres[i][j], cfc11tmp, dt_gasex);
					cfc11tmp += c11flux[i][j];
					c11monflux[i][j]+=c11flux[i][j];

					c12flux[i][j] = trac_calcflux(Sc_cfc12[i][j],
							cfcatm12[i][j], cfc12_sol[0][i][j], fice[i][j],
							xkw[i][j], atmpres[i][j], cfc12tmp, dt_gasex);
					cfc12tmp += c12flux[i][j];
					c12monflux[i][j]+=c12flux[i][j];

#endif

				} // Time loop
#ifdef SF6
				sf6_sat[i][j]=atmpres[i][j]*sf6_sol[0][i][j]*sf6atm[i][j];
#endif
#ifdef CFC
				cfc11_sat[i][j]=atmpres[i][j]*cfc11_sol[0][i][j]*cfcatm11[i][j];
				cfc12_sat[i][j]=atmpres[i][j]*cfc12_sol[0][i][j]*cfcatm12[i][j];
#endif
			} else { // No gas exchange for land tile
#ifdef SF6
				sf6flux[i][j] = 0.0;
				sf6tmp = 0.0;
				sf6_sat[i][j]=0.0;
#endif
#ifdef CFC
				c11flux[i][j] = 0.0;
				c12flux[i][j] = 0.0;

				cfc11tmp = 0.0;
				cfc12tmp = 0.0;
				cfc11_sat[i][j]=0.0;
				cfc12_sat[i][j]=0.0;
#endif
			}
			
			// Save result of gas exchange to main tracer array
#ifdef SF6
			tr[mSF6][0][i][j] = sf6tmp;
			tr[mSF6][1][i][j] = sf6tmp;

			sf6[0][i][j]=sf6tmp;
			sf6[0][i][j]=sf6tmp;


#ifdef SAT_TRACER
			    sf6_relsat=sf6tmp/sf6_sat[i][j];
			    if (isnan(sf6_relsat) || isinf(sf6_relsat)) sf6_relsat=1.0;

			    tr[mSF6_sat][0][i][j]=sf6_relsat;
			    tr[mSF6_sat][1][i][j]=sf6_relsat;
			    mn_sf6_relsat[0][i][j]=sf6_relsat;
			    mn_sf6_relsat[1][i][j]=sf6_relsat;
#endif
			

#endif
#ifdef CFC
			tr[mCFC11][0][i][j] = cfc11tmp;
			tr[mCFC11][1][i][j] = cfc11tmp;
			cfc11[0][i][j] = cfc11tmp;
			cfc11[1][i][j] = cfc11tmp;

			tr[mCFC12][0][i][j] = cfc12tmp;
			tr[mCFC12][1][i][j] = cfc12tmp;
			cfc12[0][i][j]=cfc12tmp;
			cfc12[1][i][j]=cfc12tmp;

#ifdef SAT_TRACER
			cfc11_relsat=cfc11tmp/cfc11_sat[i][j];
			if (isnan(cfc11_relsat) || isinf(cfc11_relsat)) cfc11_relsat=1.0;

			mn_cfc11_relsat[0][i][j]=cfc11_relsat;
			mn_cfc11_relsat[1][i][j]=cfc11_relsat;
			tr[mCFC11_sat][0][i][j]=cfc11_relsat;
			tr[mCFC11_sat][1][i][j]=cfc11_relsat;


			cfc12_relsat=cfc12tmp/cfc12_sat[i][j];
			if (isnan(cfc12_relsat) || isinf(cfc12_relsat)) cfc12_relsat=1.0;

			mn_cfc12_relsat[0][i][j]=cfc12_relsat;
			mn_cfc12_relsat[1][i][j]=cfc12_relsat;
			tr[mCFC12_sat][0][i][j]=cfc12_relsat;
			tr[mCFC12_sat][1][i][j]=cfc12_relsat;
#endif

#endif

		} // Y-dim loop
	} // X-dim loop

	#define TESTLON 101
	#define TESTLAT 81
	extern double lonh[NXMEM];
	printf("Lat/Lon: %g %g\n",lath[TESTLAT],lonh[TESTLON]);
	printf("Atmospheric Concentration for year %f: %e\n",thetime,cfcatm12[TESTLON][TESTLAT]);
	printf("T=%g S=%g P=%g\n",Temptm[0][TESTLON][TESTLAT],Salttm[0][TESTLON][TESTLAT],atmpres[TESTLON][TESTLAT]);
	printf("Solubility %g\n",cfc12_sol[0][TESTLON][TESTLAT]);
	printf("Mixed Layer Concentration: %g\n",tr[mCFC12][0][TESTLON][TESTLAT]);

#endif // SF6 or CFC
	// end ashao
#ifdef RADIOACTIVE
#define fukulatidx 142
#define fukulonidx 62
printf("MONTH %d\n",iterno);
	if ( iterno < 12)
	{
                for (k=0;k<=1;k++) {
                        i131[k][fukulonidx][fukulatidx]=1.;
                        cs134[k][fukulonidx][fukulatidx]=1.;
                        cs136[k][fukulonidx][fukulatidx]=1.;
                        cs137[k][fukulonidx][fukulatidx]=1.;
                        ba140[k][fukulonidx][fukulatidx]=1.;
                        la140[k][fukulonidx][fukulatidx]=1.;

                        tr[mi131][k][fukulonidx][fukulatidx]=1.;
                        tr[mcs134][k][fukulonidx][fukulatidx]=1.;
                        tr[mcs136][k][fukulonidx][fukulatidx]=1.;
                        tr[mcs137][k][fukulonidx][fukulatidx]=1.;
                        tr[mba140][k][fukulonidx][fukulatidx]=1.;
                        tr[mla140][k][fukulonidx][fukulatidx]=1.;
                }

	}
#endif

#ifdef OXYGEN
	oxygen_saturation(Temptm, Salttm, o2_sat);
	// ashao: Added time step loop to simulate gas exchange
	for (m=0;m<mdt;m++) {
	    for (i=X1;i<=nx;i++) {
		for (j=Y1;j<=ny;j++) {
		    if (D[i][j] > MINIMUM_DEPTH) {

			//ashao: Removed, Unknown source of coefficients
			/* Sc_o2[i][j] = 1638.0 + Temptm[0][i][j] *
				(-81.83 + Temptm[0][i][j] * (1.483 + Temptm[0][i][j] * (-0.008004))); */
			// ashao: Wanninkhof 1992 for Schmidt number coefficients
			
			h_mixedlyr = h[0][i][j] + h[1][i][j]; // Thickness of mixed layer (no buffer layers)
			Sc_o2[i][j] = 1953.4 - 128*Temptm[0][i][j] + 3.9918*pow(Temptm[0][i][j],2) + 0.050091*pow(Temptm[0][i][j],3);
			Sk_o2[i][j] = (1-fice[i][j])*xkw[i][j]*pow(Sc_o2[i][j]/660.0,-0.5);
			oxyflux[i][j] = Sk_o2[i][j] * (o2_sat[0][i][j] - tr[mOXYGEN][0][i][j])*dt_gasex/h_mixedlyr;
			tr[mOXYGEN][0][i][j]+=oxyflux[i][j];
			tr[mOXYGEN][1][i][j]+=oxyflux[i][j];
		    } else {
			oxyflux[i][j] = 0.e0;
		    }
		}
	    }
	}

	printf("OXYGEN SAT: %f\n",tr[mOXYGEN][0][100][100]/o2_sat[0][100][100]);
#endif

#ifdef OXY18
	for (i=X1;i<=nx;i++) {
		for (j=Y1;j<=ny;j++) {
			if (D[i][j] > MINIMUM_DEPTH) {
				o18flux[i][j] = Sk_o2[i][j] * (o2_sat[0][i][j] * ratio_o2_sw
						- tr[mo18][0][i][j]);
			} else {
				o18flux[i][j] = 0.e0;
			}
			jo18[0][i][j] = jo18[0][i][j] + o18flux[i][j] / h[0][i][j];
		}
	}
#endif

#ifdef DIC

	/* ---------------------------------------------------
	 *  calculate quantities not dependent on dic
	 * -------------------------------------------------*/

	for (i=X1;i<=nx;i++) {
		for (j=Y1;j<=ny;j++) {
			if (D[i][j] > MINIMUM_DEPTH) {
				Sc_co2[i][j] = 2073.1 + Temptm[0][i][j] *
				(-125.62 + Temptm[0][i][j] * (3.6276 + Temptm[0][i][j] * (-0.043219))); // Schmidt number
				Sk_co2[i][j] = (1.0-fice[i][j])*xkw[i][j]*pow(Sc_co2[i][j]/660.0,-0.5); // piston velocity
				tmpdic[i][j] = dic[0][i][j]; // create temporary dic value
				co2flux[i][j] = 0.e0; // set co2flux for this loop to zero
				tmpalk[i][j] = alk[0][i][j]; // create temporary alk value
			}
		}
	}

	/* ---------------------------------------------------
	 *  loop over gas exchange time steps
	 * -------------------------------------------------*/
	for (i=X1;i<=nx;i++) {
		for (j=Y1;j<=ny;j++) {
			if (D[i][j] > MINIMUM_DEPTH) {
				for (l=1;l<=mdt;l++) { // break gas exchange flux into mdt time steps
					co2calc(Temptm[0][i][j],Salttm[0][i][j],tmpdic[i][j],tmpalk[i][j],
							0.0,40.e-3,6.0,10.0,278.0,atmpres[i][j]); // calculate dco2star
					tmpco2flux = Sk_co2[i][j] * dco2star; // calculate air-sea co2 flux
					// save pco2surf for output (debugging)
					pco2s[i][j] = pco2surf;
					// add up temporary co2fluxes for net gas exchange
					co2flux[i][j] += tmpco2flux * inv_mdt;
					tmpdic[i][j] += (tmpco2flux/ h[0][i][j] + jdic[0][i][j]) * dt * inv_mdt;
					tmpalk[i][j] += jalk[0][i][j] * dt * inv_mdt;
#ifdef MERGED_ML
					tmpdic[i][j] = (tmpdic[i][j]*h[0][i][j] + tr[mDIC][1][i][j]*h[1][i][j]) /
					(h[0][i][j] + h[1][i][j]);
					tr[mDIC][1][i][j] = tmpdic[i][j];
					tmpalk[i][j] = (tmpalk[i][j]*h[0][i][j] + tr[mALK][1][i][j]*h[1][i][j]) /
					(h[0][i][j] + h[1][i][j]);
					tr[mALK][1][i][j] = tmpalk[i][j];
#endif
				} // end of loop over mdt time steps
				tr[mDIC][0][i][j] = tmpdic[i][j]; // save final "temporary" dic into tracer array
				tr[mALK][0][i][j] = tmpalk[i][j]; // save final "temporary" alk into tracer array
			}
		}
	}

	/* ---------------------------------------------------
	 *  write zero fluxes to land points
	 * -------------------------------------------------*/
	for (i=X1;i<=nx;i++) {
		for (j=Y1;j<=ny;j++) {
			if (D[i][j] < MINIMUM_DEPTH) {
				co2flux[i][j] = 0.e0;
				pco2s[i][j] = misval;
			}
		}
	}

# ifdef MERGED_ML
	merge_ml_tr();
# endif
#endif
#else /* GASEX */
	// begin ashao
	printf("No gas exchange: setting mixed layer to saturation\n");
	thetime = currtime+dt/31536000;
	//yearfrac = currtime-floor(currtime)+dt/31536000; // calculate what fraction of year the current timestep is

#ifdef OXYGEN
	oxygen_saturation(Temptm, Salttm, o2_sat);
	// ashao: Added time step loop to simulate gas exchange
	    for (i=X1;i<=nx;i++) {
		for (j=Y1;j<=ny;j++) {
		    if (D[i][j]>MINIMUM_DEPTH) {
			
			tr[mOXYGEN][0][i][j]=o2_sat[0][i][j];
			tr[mOXYGEN][1][i][j]=o2_sat[0][i][j];
		    } else {
			oxyflux[i][j] = 0.e0;
		    }
		}
	    }

	printf("OXYGEN SAT: %f\n",tr[mOXYGEN][0][100][100]);
#endif
	for (i = 0; i < NXMEM; i++) {
		for (j = 0; j < NYMEM; j++) {
			// Calculations for 2D arrays (LAT x LON)
			// Property values at grid location
			TempK = Temptm[0][i][j] + 273.15; // Temperature at current grid location
			// Calculate tracer solubility and saturation value
#ifdef CFC
			cfc11_sol[0][i][j] = trac_calcsol(TempK, Salttm[0][i][j],
					cfc11sol_coeffs);
			cfc12_sol[0][i][j] = trac_calcsol(TempK, Salttm[0][i][j],
					cfc12sol_coeffs);
#endif // CFC
#ifdef SF6
			sf6_sol[0][i][j] = trac_calcsol(TempK, Salttm[0][i][j],
					sf6sol_coeffs);
#endif //SF6
			// Calculate atmospheric concentration
#ifdef CFC
			cfcatm11[i][j] = trac_atmospheric(thetime, cfc11nrecs,
					cfc11year, cfc11N, cfc11S, j);
			cfcatm12[i][j] = trac_atmospheric(thetime,cfc12nrecs,
					cfc12year, cfc12N, cfc12S, j);
#endif
#ifdef SF6
			sf6atm[i][j] = trac_atmospheric(thetime, sf6nrecs,
					sf6year, sf6N, sf6S, j);
#endif // SF6

			if (D[i][j]>MINIMUM_DEPTH) {
				 // Calculate saturation value for ocean tile
#ifdef SF6
				sf6tmp = atmpres[i][j] * sf6atm[i][j] * sf6_sol[0][i][j];
#endif

#ifdef CFC
				cfc11tmp = atmpres[i][j] * cfcatm11[i][j] * cfc11_sol[0][i][j];
				cfc12tmp = atmpres[i][j] * cfcatm12[i][j] * cfc12_sol[0][i][j];
#endif

			} else { // No gas exchange for land tile
#ifdef SF6
				sf6flux[i][j] = 0.0;
				sf6tmp = 0.0;
#endif
#ifdef CFC
				c11flux[i][j] = 0.0;
				c12flux[i][j] = 0.0;

				cfc11tmp = 0.0;
				cfc12tmp = 0.0;
#endif
			}
			// Save saturation value to tracer array
#ifdef SF6
			tr[mSF6][0][i][j] = sf6tmp;
			tr[mSF6][1][i][j] = sf6tmp;
#endif
#ifdef CFC
			tr[mCFC11][0][i][j] = cfc11tmp;
			tr[mCFC11][1][i][j] = cfc11tmp;

			tr[mCFC12][0][i][j] = cfc12tmp;
			tr[mCFC12][1][i][j] = cfc12tmp;
#endif

		} // Y-dim loop
	} // X-dim loop
/*
	#define TESTLON 101
	#define TESTLAT 81
	extern double lonh[NXMEM];
	printf("Lat/Lon: %g %g\n",lath[TESTLAT],lonh[TESTLON]);
	printf("Atmospheric Concentration for year %f: %e\n",thetime,cfcatm12[TESTLON][TESTLAT]);
	printf("T=%g S=%g P=%g\n",Temptm[0][TESTLON][TESTLAT],Salttm[0][TESTLON][TESTLAT],atmpres[TESTLON][TESTLAT]);
	printf("Solubility %g\n",cfc12_sol[0][TESTLON][TESTLAT]);
	printf("Mixed Layer Concentration: %g\n",tr[mCFC12][0][TESTLON][TESTLAT]);
	// end ashao
*/

#endif /* GASEX */

// begin ashao
#ifdef RADIOACTIVE
// Calculate decay constants in inverse seconds;

double lambda_i131=log(2.)/(8.04*60*60*24);
double lambda_cs134=log(2.)/(2.0652*60*60*24*365);
double lambda_cs136=log(2.)/(13.1*60*60*24);
double lambda_cs137=log(2.)/(30*60*60*24*365);
double lambda_ba140=log(2.)/(12.74*60*60*24);
//double lambda_la140=log(2.)/(40.272*60*60);
double lambda_la140=0;

printf("lambda_i131=%e\n",lambda_i131);
printf("cs137 idx=%d",mcs137);
printf("Starting Conc dt=%e\n",dt);
printf("cs137=%e after %e\n",cs137[0][65][150],cs137[0][62][145]*exp(-lambda_cs137*dt));
for (k=0;k<NZ;k++){
	for (i = 0; i < NXMEM; i++) {
		for (j = 0; j < NYMEM; j++) {


			tr[mi131][k][i][j]*=exp(-lambda_i131*dt);
			tr[mcs134][k][i][j]*=exp(-lambda_cs134*dt);
			tr[mcs136][k][i][j]*=exp(-lambda_cs136*dt);
			tr[mcs137][k][i][j]*=exp(-lambda_cs137*dt);
			tr[mba140][k][i][j]*=exp(-lambda_ba140*dt);
			tr[mla140][k][i][j]*=exp(-lambda_la140*dt);
			
			i131[k][i][j]=tr[mi131][k][i][j];
			cs134[k][i][j]=tr[mcs134][k][i][j];
			cs136[k][i][j]=tr[mcs136][k][i][j];
			cs137[k][i][j]=tr[mcs137][k][i][j];
			ba140[k][i][j]=tr[mba140][k][i][j];
			la140[k][i][j]=tr[mla140][k][i][j];
		}
	}
}

printf("Concentrations at Fukushima\n");
printf("131-I: %e\n",tr[mi131][0][62][142]);
printf("134-Cs: %e\n",tr[mcs134][0][62][142]);
printf("136-Cs: %e\n",tr[mcs136][0][62][142]);
printf("137-Cs: %e\n",tr[mcs137][0][62][142]);
printf("140-Ba: %e\n",tr[mba140][0][62][142]);
printf("140-La: %e\n",tr[mla140][0][62][142]);



#endif

#ifdef TTD_BP
for (k=0;k<NZ;k++)
	for (i = 0; i < NXMEM; i++)
		for (j = 0; j < NYMEM; j++)
			ttd[k][i][j]=tr[mttd][k][i][j];
#endif

// end ashao
	/*-----------------------------------------
	 *
	 *     OVERWRITE RESULTS IN SPECIAL CASES
	 *
	 *-----------------------------------------*/

#ifdef INERT
	for (i=0;i<=NXMEM-1;i++) {
		for (j=0;j<=NYMEM-1;j++) {
			if (D[i][j]>MINIMUM_DEPTH) {
				tr[mINERT2][0][i][j] += (1.e0 / h[0][i][j]) * dt;
				for (k=0;k<NZ;k++) {
					tr[mINERT2][k][i][j] -= tr[mINERT2][k][i][j]/(53.e0 * 86400.e0) * dt;
				}
			} else {
				for (k=0;k<NZ;k++)
				tr[mINERT2][k][i][j] = misval;
			}
		}
	}
#endif

/* DEPRECATED ashao
#ifdef OXYGEN
	for (i=0;i<=NXMEM-1;i++) {
		for (j=0;j<=NYMEM-1;j++) {
			//BX - reinstated by HF
			if (D[i][j]>MINIMUM_DEPTH) {
				tr[mOXYGEN][0][i][j] = o2_sat[0][i][j];
				//BX - reinstated by HF
			} else {
				tr[mOXYGEN][0][i][j] = misval;
			}
		}
	}
#endif
*/

#ifdef OXY18
	for (i=0;i<=NXMEM-1;i++) {
		for (j=0;j<=NYMEM-1;j++) {
			//BX - reinstated by HF
			if (D[i][j]>MINIMUM_DEPTH) {
				tr[mo18][0][i][j] = tr[mOXYGEN][0][i][j] * Rsmow * (1.0+24e-3);
				//BX - reinstated by HF
			} else {
				tr[mo18][0][i][j] = misval;
			}
		}
	}
#endif

#ifdef UNIVENT
	// Set AOU at isopycnal "outcrops" to climatological value
	// This removes the effect of variable ventilation
	calc_ksub(h_clim);
	for (i=0;i<=NXMEM-1;i++) {
		for (j=0;j<=NYMEM-1;j++) {
			ksub_clim[i][j] = ksub[i][j];
		}
	}
	calc_ksub(h);
	for (k=0;k<NZ-1;k++) {/* ksub=NZ on land => exclude bottom layer in loop */
		for (i=0;i<=NXMEM-1;i++) {
			for (j=0;j<=NYMEM-1;j++) {
				if (k == ksub_clim[i][j] || k == ksub[i][j]) {
					tr[mOXYGEN][k][i][j] = o2_sat[k][i][j]-aou_spec[k][i][j];
				}
			}
		}
	}
#endif

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

#ifdef OXY18
	for (k=0;k<NZ;k++)
	for (i=0;i<=NXMEM-1;i++)
	for (j=0;j<=NYMEM-1;j++)
	tr[mo18][k][i][j] = Ro18o16[k][i][j]*tr[mOXYGEN][k][i][j];
#endif

#ifdef REST_ARCTIC
	/*-----------------------------------------
	 *
	 *    Restore Arctic to climatologies
	 *    No time step is included here:
	 *    Restoring with one time step is assumed
	 *
	 *-----------------------------------------*/
	//BX     filter with value of 1 north of 72 degrees and
	//BX     a value of 0 at and south of starting latitude
	//BX        (slat in following formula)
	//BX     wfunct = sin[(PI* (lat - slat)) / (2 * (72 - slat))]
	//BX        slat = 66 --> North of Bering Strait and Iceland
	for (i=0;i<=NXMEM-1;i++) {
		for (j=0;j<=NYMEM-1;j++) {
			if (qlat[j] >= 66. && qlat[j] <= 72.) {wfunct[j] = sin((M_PI*(qlat[j]-66.))/12.);}
			else if (qlat[j] > 72.) {wfunct[j] = 1.;}
			else {wfunct[j] = 0.;}
			//if (i==45) printf("Restoring with wfunct %f at j= %i, qlat= %f \n",wfunct[j],j,qlat[j]);
			for(k=0;k<NZ;k++) {
# ifdef DIC
				tr[mDIC][k][i][j] += wfunct[j] * (dic_lay[k][i][j]- tr[mDIC][k][i][j]);
				tr[mALK][k][i][j] += wfunct[j] * (alk_lay[k][i][j]- tr[mALK][k][i][j]);
# endif
# ifdef PHOSPHATE
				tr[mPHOSPHATE][k][i][j] += wfunct[j] * (po4_star_lay[k][i][j]- tr[mPHOSPHATE][k][i][j]);
# endif
			}
		}
	}
# ifdef DIC
	//printf("Restoring Arctic for DIC and Alk. \n");
# endif
# ifdef PHOSPHATE
	//printf("Restoring Arctic for PO4. \n");
# endif

#endif

	// Calculate inventories of tracers
#ifdef OXYGEN
	printf("Oxygen inventory: %e\n", calc_inventory(mOXYGEN));
#endif
#ifdef PHOSPHATE
	printf("Phosphate inventory: %e\n", calc_inventory(mPHOSPHATE));
	printf("DOP inventory: %e\n", calc_inventory(mDOP));
#endif

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

#ifdef OXY18
	isotope_R_del(mOXYGEN,mo18,Rsmow,Ro18o16,delo18);
#endif

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
#ifdef OXY18
				o18[k][i][j] = tr[mo18][k][i][j];
#endif
#ifdef DIC
				dic[k][i][j] = tr[mDIC][k][i][j];
				alk[k][i][j] = tr[mALK][k][i][j];
#endif
#ifdef AGE
# if defined AGEEXP || defined AGE1 || defined AGE2
				age[k][i][j] = tr[mAGE][k][i][j];
# else
				age[k][i][j] = tr[mAGE][k][i][j] / (365.0 * 86400.0); /*years*/
# endif
#endif
#ifdef CFC
				cfc11[k][i][j] = tr[mCFC11][k][i][j];

#ifdef SAT_TRACER
				if ( isnan(tr[mCFC11_sat][k][i][j]) ) tr[mCFC11_sat][k][i][j]=1.0;
				mn_cfc11_relsat[k][i][j] = tr[mCFC11_sat][k][i][j];
#endif
				if (k == 0 && i == (int) (3 * NXMEM / 4) && j == (int) (3
						*NYMEM / 4)) {
					printf("Example CFC-11 = %e %e %e %e\n",
							cfc11[0][(int) (NXMEM / 4)][(int) (NYMEM / 4)],
							cfc11[0][(int) (NXMEM / 4)][(int) (2 * NYMEM / 4)],
							cfc11[0][(int) (3 * NXMEM / 4)][(int) (NYMEM / 4)],
							cfc11[0][(int) (3 * NXMEM / 4)][(int) (2 * NYMEM
									/ 4)]);
				}
#ifdef SAT_TRACER
				if ( isnan(tr[mCFC12_sat][k][i][j]) ) tr[mCFC12_sat][k][i][j]=1.0;
				mn_cfc12_relsat[k][i][j] = tr[mCFC12_sat][k][i][j];
#endif
				cfc12[k][i][j] = tr[mCFC12][k][i][j];
				if (k == 0 && i == (int) (3 * NXMEM / 4) && j == (int) (3
						*NYMEM / 4)) {
					printf("Example CFC-12 = %e %e %e %e\n",
							//cfc12[0][(int) (NXMEM / 4)][(int) (NYMEM / 4)],
							cfc12[0][100][100],
							cfc12[0][(int) (NXMEM / 4)][(int) (2 * NYMEM / 4)],
							cfc12[0][(int) (3 * NXMEM / 4)][(int) (NYMEM / 4)],
							cfc12[0][(int) (3 * NXMEM / 4)][(int) (2 * NYMEM
									/ 4)]);
				}
#endif
				// begin ashao
#ifdef SF6
#ifdef SAT_TRACER
				if ( isnan(tr[mSF6_sat][k][i][j])) tr[mSF6_sat][k][i][j]=1.0;
				mn_sf6_relsat[k][i][j] = tr[mSF6_sat][k][i][j];
#endif
				sf6[k][i][j] = tr[mSF6][k][i][j];
				if (k == 0 && i == (int) (3 * NXMEM / 4) && j == (int) (3
						*NYMEM / 4)) {
					printf(
							"Example SF6 = %e %e %e %e\n",
							sf6[0][(int) (NXMEM / 4)][(int) (NYMEM / 4)],
							sf6[0][(int) (NXMEM / 4)][(int) (2 * NYMEM / 4)],
							sf6[0][(int) (3 * NXMEM / 4)][(int) (NYMEM / 4)],
							sf6[0][(int) (3 * NXMEM / 4)][(int) (2 * NYMEM / 4)]);
				}

#endif
				// end ashao
#ifdef INERT
				inert1[k][i][j] = tr[mINERT1][k][i][j];
				inert2[k][i][j] = tr[mINERT2][k][i][j];
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
# ifdef PHOSPHATE
					if (tr[mPHOSPHATE][k][i][j] < trwarn[mPHOSPHATE][1]) {
						if (pquant[3] == 9) {
							printf("WARNING: RAISING LOW [PO4] TO %g \n", trwarn[mPHOSPHATE][1]);
							printf(" i,j,k,h,depth,tr[po4] = %i,%i,%i,%g,%g,%g \n",
									i,j,k,h[k][i][j],depth[k][i][j],tr[mPHOSPHATE][k][i][j]);
						}
						tr[mPHOSPHATE][k][i][j]=trwarn[mPHOSPHATE][1];
					}
# endif
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

# ifdef PHOSPHATE
	printf("SINGLE GRID POINT: po4,dop,o2,o18,del18o = %g,%g,%g,%g,%g \n",
			tr[mPHOSPHATE][10][100][40],tr[mDOP][10][100][40],tr[mOXYGEN][10][100][40],
			tr[mo18][10][100][40],delo18[10][100][40]);
# endif

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
				if (D[i][j]>MINIMUM_DEPTH) {
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
			if (D[i][j]>MINIMUM_DEPTH) {
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
			if (D[i][j]>MINIMUM_DEPTH) {
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

#ifdef DIC
void co2calc(double t,double s,double dic_in,double ta_in,
		double pt_in,double sit_in,
		double phlo,double phhi,double xco2_in,double atmpres)

{

	/********+*********+*********+*********+*********+*********+*********+*
	 *
	 * SUBROUTINE CO2CALC
	 *
	 * PURPOSE
	 *	Calculate delta co2* from total alkalinity and total CO2 at
	 * temperature (t), salinity (s) and "atmpres" atmosphere total pressure.
	 *
	 * USAGE
	 *       call co2calc(t,s,dic_in,ta_in,pt_in,sit_in
	 *    &                  ,phlo,phhi,ph,xco2_in,atmpres
	 *    &                  ,co2star,dco2star,pCO2surf,dpco2)
	 *
	 * INPUT
	 *	dic_in = total inorganic carbon (mol/m^3)
	 *                where 1 T = 1 metric ton = 1000 kg
	 *	ta_in  = total alkalinity (eq/m^3)
	 *	pt_in  = inorganic phosphate (mol/m^3)
	 *	sit_in = inorganic silicate (mol/m^3)
	 *	t      = temperature (degrees C)
	 *	s      = salinity (PSU)
	 *	phlo   = lower limit of pH range
	 *	phhi   = upper limit of pH range
	 *	xco2_in=atmospheric mole fraction CO2 in dry air (ppmv)
	 *	atmpres= atmospheric pressure in atmospheres (1 atm==1013.25mbar)
	 *
	 *       Note: arguments dic_in, ta_in, pt_in, sit_in, and xco2_in are
	 *             used to initialize variables dic, ta, pt, sit, and xco2.
	 *             * Variables dic, ta, pt, and sit are in the common block
	 *               "species".
	 *             * Variable xco2 is a local variable.
	 *             * Variables with "_in" suffix have different units
	 *               than those without.
	 *
	 * OUTPUT
	 *	co2star  = CO2*water (mol/m^3)
	 *	dco2star = delta CO2 (mol/m^3)
	 *       pco2surf = oceanic pCO2 (ppmv)
	 *       dpco2    = Delta pCO2, i.e, pCO2ocn - pCO2atm (ppmv)
	 *
	 * IMPORTANT: Some words about units - (J.C.Orr, 4/4/1999)
	 *     - Models carry tracers in mol/m^3 (on a per volume basis)
	 *     - Conversely, this routine, which was written by observationalists
	 *       (C. Sabine and R. Key), passes input arguments in umol/kg
	 *       (i.e., on a per mass basis)
	 *     - I have changed things slightly so that input arguments are in mol/m^3,
	 *     - Thus, all input concentrations (dic_in, ta_in, pt_in, and st_in)
	 *       should be given in mol/m^3; output arguments "co2star" and "dco2star"
	 *       are likewise in mol/m^3.
	 *
	 * FILES and PROGRAMS NEEDED
	 *	drtsafe
	 *	ta_iter_1
	 *
	 * $Source: /home/ashao/uw-apl/CVS/offsf6/src/step.c,v $
	 *   converted to C code, cdeutsch 2006
	 *
	 *************************************************************************/

	extern double k0,k1,k2,kb,kw,ks,kf,k1p,k2p,k3p,ksi;
	extern double ff,htotal,htotal2,xco2,xacc,x1,x2;
	extern double bt,st,ft,sit,pt,dic1,ta,co2starair;

	double tk,tk100,tk1002,invtk,dlogtk;
	double is,is2,s2,sqrts,sqrtis,s15,scl;
	double permil,permeg;
	// these varaiables are now just local (removed from offtrac.c)
	double dpco2,ph,co2star;

	void ta_iter_1(double xacc, double *fn, double *df);
	double drtsafe(void (*ta_iter_1)(double, double *, double *),
			double x1, double x2, double xacc);

	//*************************************************************************
	//       Change units from the input of mol/m^3 -> mol/kg:
	//       (1 mol/m^3)  x (1 m^3/1024.5 kg)
	//       where the ocean's mean surface density is 1024.5 kg/m^3
	//       Note: mol/kg are actually what the body of this routine uses
	//       for calculations.
	//*************************************************************************

	permil = 1.0 / 1024.5;
	pt=pt_in*permil;
	sit=sit_in*permil;
	ta=ta_in*permil;
	dic1=dic_in*permil;

	//*************************************************************************
	//       Change units from uatm to atm. That is, atm is what the body of
	//       this routine uses for calculations.
	//*************************************************************************

	permeg=1.e-6;
	xco2=xco2_in*permeg;

	//*************************************************************************
	// Calculate all constants needed to convert between various measured
	// carbon species. References for each equation are noted in the code.
	// Once calculated, the constants are
	// stored and passed in the common block "const". The original version of this
	// code was based on the code by Dickson in Version 2 of "Handbook of Methods
	// for the Analysis of the Various Parameters of the Carbon Dioxide System
	// in Seawater", DOE, 1994 (SOP No. 3, p25-26).
	//*************************************************************************

	// Quantities used more than once

	tk = 273.15 + t;
	tk100 = tk/100.0;
	tk1002=tk100*tk100;
	invtk=1.0/tk;
	dlogtk=log(tk);
	is=19.924*s/(1000.-1.005*s);
	is2=is*is;
	sqrtis=sqrt(is);
	s2=s*s;
	sqrts=sqrt(s);
	s15=pow(s,1.5);
	scl=s/1.80655;

	// f = k0(1-pH2O)*correction term for non-ideality
	// Weiss & Price (1980, Mar. Chem., 8, 347-359; Eq 13 with table 6 values)

	ff = exp(-162.8301 + 218.2968/tk100 +
			90.9241*log(tk100) - 1.47696*tk1002 +
			s * (.025695 - .025225*tk100 +
					0.0049867*tk1002));

	// K0 from Weiss 1974

	k0 = exp(93.4517/tk100 - 60.2409 + 23.3585 * log(tk100) +
			s * (.023517 - 0.023656 * tk100 + 0.0047036 * tk1002));

	// k1 = [H][HCO3]/[H2CO3]
	// k2 = [H][CO3]/[HCO3]
	// Millero p.664 (1995) using Mehrbach et al. data on seawater scale

	k1=pow(10,-1*(3670.7*invtk - 62.008 + 9.7944*dlogtk -
					0.0118 * s + 0.000116*s2));

	k2=pow(10,-1*(1394.7*invtk + 4.777 -
					0.0184*s + 0.000118*s2));

	// kb = [H][BO2]/[HBO2]
	// Millero p.669 (1995) using data from Dickson (1990)

	kb=exp((-8966.90 - 2890.53*sqrts - 77.942*s +
					1.728*s15 - 0.0996*s2)*invtk +
			(148.0248 + 137.1942*sqrts + 1.62142*s) +
			(-24.4344 - 25.085*sqrts - 0.2474*s) *
			dlogtk + 0.053105*sqrts*tk);

	// k1p = [H][H2PO4]/[H3PO4]
	// DOE(1994) eq 7.2.20 with footnote using data from Millero (1974)

	k1p = exp(-4576.752*invtk + 115.525 - 18.453 * dlogtk +
			(-106.736*invtk + 0.69171) * sqrts +
			(-0.65643*invtk - 0.01844) * s);

	// k2p = [H][HPO4]/[H2PO4]
	// DOE(1994) eq 7.2.23 with footnote using data from Millero (1974))

	k2p = exp(-8814.715*invtk + 172.0883 - 27.927 * dlogtk +
			(-160.340*invtk + 1.3566) * sqrts +
			(0.37335*invtk - 0.05778) * s);

	// k3p = [H][PO4]/[HPO4]
	// DOE(1994) eq 7.2.26 with footnote using data from Millero (1974)

	k3p = exp(-3070.75*invtk - 18.141 +
			(17.27039*invtk + 2.81197) *
			sqrts + (-44.99486*invtk - 0.09984) * s);

	// ksi = [H][SiO(OH)3]/[Si(OH)4]
	// Millero p.671 (1995) using data from Yao and Millero (1995)

	ksi = exp(-8904.2*invtk + 117.385 - 19.334 * dlogtk +
			(-458.79*invtk + 3.5913) * sqrtis +
			(188.74*invtk - 1.5998) * is +
			(-12.1652*invtk + 0.07871) * is2 +
			log(1.0-0.001005*s));

	// kw = [H][OH]
	// Millero p.670 (1995) using composite data

	kw = exp(-13847.26*invtk + 148.9652 - 23.6521 * dlogtk +
			(118.67*invtk - 5.977 + 1.0495 * dlogtk) *
			sqrts - 0.01615 * s);

	// ks = [H][SO4]/[HSO4]
	// Dickson (1990, J. chem. Thermodynamics 22, 113)

	ks=exp(-4276.1*invtk + 141.328 - 23.093*dlogtk +
			(-13856*invtk + 324.57 - 47.986*dlogtk) * sqrtis +
			(35474*invtk - 771.54 + 114.723*dlogtk) * is -
			2698*invtk*pow(is,1.5) + 1776*invtk*is2 +
			log(1.0 - 0.001005*s));

	// kf = [H][F]/[HF]
	// Dickson and Riley (1979) -- change pH scale to total

	kf=exp(1590.2*invtk - 12.641 + 1.525*sqrtis +
			log(1.0 - 0.001005*s) +
			log(1.0 + (0.1400/96.062)*(scl)/ks));

	// Calculate concentrations for borate, sulfate, and fluoride
	// Uppstrom (1974)
	bt = 0.000232 * scl/10.811;
	// Morris & Riley (1966)
	st = 0.14 * scl/96.062;
	// Riley (1965)
	ft = 0.000067 * scl/18.9984;

	/*************************************************************************
	 *
	 * Calculate [H+] total when DIC and TA are known at T, S and 1 atm.
	 * The solution converges to err of xacc. The solution must be within
	 * the range x1 to x2.
	 *
	 * If DIC and TA are known then either a root finding or iterative method
	 * must be used to calculate htotal. In this case we use the Newton-Raphson
	 * "safe" method taken from "Numerical Recipes" (function "rtsafe.f" with
	 * error trapping removed).
	 *
	 * As currently set, this procedure iterates about 12 times. The x1 and x2
	 * values set below will accomodate ANY oceanographic values. If an initial
	 * guess of the pH is known, then the number of iterations can be reduced to
	 * about 5 by narrowing the gap between x1 and x2. It is recommended that
	 * the first few time steps be run with x1 and x2 set as below. After that,
	 * set x1 and x2 to the previous value of the pH +/- ~0.5. The current
	 * setting of xacc will result in co2star accurate to 3 significant figures
	 * (xx.y). Making xacc bigger will result in faster convergence also, but this
	 * is not recommended (xacc of 10**-9 drops precision to 2 significant figures).
	 *
	 * Parentheses added around negative exponents (Keith Lindsay)
	 *
	 *************************************************************************/

	x1 = pow(10.0,-phhi);
	x2 = pow(10.0,-phlo);
	xacc = 1.e-10;
	htotal = (*drtsafe)(ta_iter_1,x1,x2,xacc);

	// Calculate [CO2*] as defined in DOE Methods Handbook 1994 Ver.2,
	// ORNL/CDIAC-74, Dickson and Goyet, eds. (Ch 2 p 10, Eq A.49)

	htotal2=htotal*htotal;
	co2star=dic1*htotal2/(htotal2 + k1*htotal + k1*k2);
	co2starair=xco2*ff*atmpres;
	dco2star=co2starair-co2star;
	ph=-log10(htotal);

	//     ---------------------------------------------------------------
	//      Add two output arguments for storing pCO2surf
	//      Should we be using K0 or ff for the solubility here?
	//     ---------------------------------------------------------------

	pco2surf = co2star / ff;
	dpco2 = pco2surf - xco2*atmpres;

	//  Convert units of output arguments
	//      Note: co2star and dco2star are calculated in mol/kg within this routine
	//      Thus Convert now from mol/kg -> mol/m^3

	co2star = co2star / permil;
	dco2star = dco2star / permil;

	//      Note: pCO2surf and dpCO2 are calculated in atm above.
	//      Thus convert now to uatm

	pco2surf = pco2surf / permeg;
	dpco2 = dpco2 / permeg;

}
#endif   // DIC for subroutine co2calc
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
				/* #  ifdef PROGNOSTIC
				 nitrlim
				 ironlim[k][i][j] = (ironlim[k][i][j]*h[k][i][j]+
				 ironlim[k+1][i][j]*h[k+1][i][j])/
				 (h[k][i][j] + h[k+1][i][j]);
				 ironlim[k+1][i][j] = ironlim[k][i][j];
				 #  endif */
# endif /* PHOSPHATE */
# ifdef NITRATE 
				jno3[k][i][j] = (jno3[k][i][j]*h[k][i][j]+
						jno3[k+1][i][j]*h[k+1][i][j])/
				(h[k][i][j] + h[k+1][i][j]);
				jno3[k+1][i][j] = jno3[k][i][j];

				jn2fix[k][i][j] = (jn2fix[k][i][j]*h[k][i][j]+
						jn2fix[k+1][i][j]*h[k+1][i][j])/
				(h[k][i][j] + h[k+1][i][j]);
				jn2fix[k+1][i][j] = jn2fix[k][i][j];

				jdenits[k][i][j] = (jdenits[k][i][j]*h[k][i][j]+
						jdenits[k+1][i][j]*h[k+1][i][j])/
				(h[k][i][j] + h[k+1][i][j]);
				jdenits[k+1][i][j] = jdenits[k][i][j];

				jdenitw[k][i][j] = (jdenitw[k][i][j]*h[k][i][j]+
						jdenitw[k+1][i][j]*h[k+1][i][j])/
				(h[k][i][j] + h[k+1][i][j]);
				jdenitw[k+1][i][j] = jdenitw[k][i][j];

				jdon[k][i][j] = (jdon[k][i][j]*h[k][i][j]+
						jdon[k+1][i][j]*h[k+1][i][j])/
				(h[k][i][j] + h[k+1][i][j]);
				jdon[k+1][i][j] = jdon[k][i][j];
				jremdon[k][i][j] = (jremdon[k][i][j]*h[k][i][j]+
						jremdon[k+1][i][j]*h[k+1][i][j])/
				(h[k][i][j] + h[k+1][i][j]);
				jremdon[k+1][i][j]= jremdon[k][i][j];

# ifdef N15_CYCLE
				jno3_15n[k][i][j] = (jno3_15n[k][i][j]*h[k][i][j]+
						jno3_15n[k+1][i][j]*h[k+1][i][j])/
				(h[k][i][j] + h[k+1][i][j]);
				jno3_15n[k+1][i][j] = jno3_15n[k][i][j];

				jdon_15n[k][i][j] = (jdon_15n[k][i][j]*h[k][i][j]+
						jdon_15n[k+1][i][j]*h[k+1][i][j])/
				(h[k][i][j] + h[k+1][i][j]);
				jdon_15n[k+1][i][j] = jdon_15n[k][i][j];
				jremdon_15n[k][i][j] = (jremdon_15n[k][i][j]*h[k][i][j]+
						jremdon_15n[k+1][i][j]*h[k+1][i][j])/
				(h[k][i][j] + h[k+1][i][j]);
				jremdon_15n[k+1][i][j]= jremdon_15n[k][i][j];

# endif
#endif /* NITRATE */
# ifdef DIC
				jdic[k][i][j] = (jdic[k][i][j]*h[k][i][j]+
						jdic[k+1][i][j]*h[k+1][i][j])/
				(h[k][i][j] + h[k+1][i][j]);
				jdic[k+1][i][j] = jdic[k][i][j];

				jalk[k][i][j] = (jalk[k][i][j]*h[k][i][j]+
						jalk[k+1][i][j]*h[k+1][i][j])/
				(h[k][i][j] + h[k+1][i][j]);
				jalk[k+1][i][j] = jalk[k][i][j];
# endif
# ifdef OXYGEN
				jo2[k][i][j] = (jo2[k][i][j]*h[k][i][j]+
						jo2[k+1][i][j]*h[k+1][i][j])/
				(h[k][i][j] + h[k+1][i][j]);
				jo2[k+1][i][j] = jo2[k][i][j];
#  ifdef OXY18
				jo18[k][i][j] = (jo18[k][i][j]*h[k][i][j]+
						jo18[k+1][i][j]*h[k+1][i][j])/
				(h[k][i][j] + h[k+1][i][j]);
				jo18[k+1][i][j] = jo18[k][i][j];
#  endif
# endif
			}
		}
	}
}

// begin ashao
# if defined(SF6) || defined(CFC)// begin ashao
double trac_calcsol(double TempK, double S, double sol_coeffs[7]) {

	double trac_sol;

	// Calculate solubility
	trac_sol = sol_coeffs[0] + sol_coeffs[1] * (100 / TempK) + sol_coeffs[2]
			* log(TempK / 100) + sol_coeffs[3] * (TempK / 100) * (TempK / 100)
			+ S * (sol_coeffs[4] + sol_coeffs[5] * (TempK / 100)
					+ sol_coeffs[6] * (TempK / 100) * (TempK / 100));
	trac_sol = exp(trac_sol);

	return trac_sol;
}
double trac_atmospheric(double thetime, int tracrecs,
		int tracyear[MAXRECS], double tracN[MAXRECS], double tracS[MAXRECS],
		int latidx) {
	int intpidx;
	int year, yidx; // Actual year, corresponding index
	double Nval, Sval, Nlat=10.2, Slat=-10.2; // North/south concentrations/latitudes for the given year
	double m; // Slope of interpolating line
	double tracatm;
	double currentlat;
	double yearfrac;
	extern double lath[NXTOT];

	yearfrac=thetime-floor(thetime);
	currentlat=lath[latidx];
	year = floor(thetime);
	// Find index of the current year
	yidx = year - tracyear[0];
	if (yidx < 0) {
		yidx = -1;
	}
	if (yidx > tracrecs-1) {
		yidx = tracrecs - 1;
	}

	// Determine whether to use the next year or previous year in the linear interpolation, assuming midyear values
	if (yidx > 0 && yidx < (tracrecs - 1)) {

		if (thetime<=(tracyear[tracrecs-1]+0.5)) {
	
			if (yearfrac > 0.5) intpidx = yidx + 1;
			if (yearfrac < 0.5) intpidx = yidx - 1;

			Nval = linear_interp(tracyear[yidx] + 0.5, tracN[yidx], tracyear[intpidx]
					+ 0.5, tracN[intpidx], thetime);
			Sval = linear_interp(tracyear[yidx] + 0.5, tracS[yidx], tracyear[intpidx]
					+ 0.5, tracS[intpidx], thetime);
		}
		else {
			Nval=tracN[tracrecs-1];
			Sval=tracS[tracrecs-1];
		}
		if (currentlat < Slat) {
			tracatm = Sval;
		}
		if (currentlat > Nlat) {
			tracatm = Nval;
		}
		if (currentlat >= Slat && currentlat <= Nlat) {
			tracatm = linear_interp(Nlat, Nval, Slat, Sval, currentlat); //Linearly interpolate across equator
			//printf("I'M INTERPOLATING! %f %f %f\n",tracatm, Nval, Sval);
		}
	}
	else if (yidx==-1) {
		tracatm=0.;
	}
	else  printf("ERROR!!!\n YEAR=%f",thetime);


	return tracatm;

}

double trac_schmidt(double T, double sc_coeffs[4]) {
	double trac_sc;
	trac_sc = sc_coeffs[0] - sc_coeffs[1] * T + sc_coeffs[2] * T * T
			- sc_coeffs[3] * T * T * T;
	return trac_sc;
}
double trac_calcflux(double trac_Sc, double trac_atm, double trac_sol,
		double icefrac, double pistonvel, double pres, double trac_con,
		double dt) {
	double trac_Sk; // Calculated Piston Velocity
	double trac_sat; // Saturation value
	double trac_flux; // Flux into mixed layer

	trac_Sk = (1 - icefrac) * pistonvel * pow(trac_Sc / 660.0, -0.5);
	trac_sat = pres * trac_atm * trac_sol;
	trac_flux = trac_Sk * (trac_sat - trac_con) / h_mixedlyr * dt;

	return trac_flux;
}
# endif // SF6 || CFC

double calc_inventory( int idx ) {

	int i,j,k;
	double inventory;

	inventory = 0;
	for (i=0;i<NXMEM;i++)
		for(j=0;j<NYMEM;j++)
			for(k=0;k<NZ;k++)
				inventory += wetmask[i][j]*areagr[i][j]*hend[k][i][j]*tr[idx][k][i][j];

	return inventory;

}
