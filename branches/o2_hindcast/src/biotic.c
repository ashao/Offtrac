/* compute biological source minus sink terms according to OCMIP protocol */

#include <stdio.h>
#include <math.h>
#include "init.h"
#include "step.h"
#include "phosphate.h"

extern const double misval;
extern int wetmask[NXMEM][NYMEM];

# ifdef N15_CYCLE
extern const double parm_n15_std_fraction;
extern const double parm_alpha_n15_prod;
extern const double parm_alpha_n15_pon_remin;
extern const double parm_alpha_n15_don_remin;
extern const double parm_alpha_n15_denitw;
extern const double parm_alpha_n15_denits;
extern const double parm_alpha_n15_n2fix;
# endif

#ifdef PHOSPHATE 

extern double ****tr;
extern double D[NXMEM][NYMEM];
extern double h[NZ][NXMEM][NYMEM];
// extern double po4_star_lev[NZPHOS][NXMEM][NYMEM];
// extern double po4_star_lay[NZ][NXMEM][NYMEM];
extern double Temptm[NZ][NXMEM][NYMEM];
# ifdef PROGNOSTIC
extern double fe_lev[NZPHOS][NXMEM][NYMEM];
extern double fe_lay[NZ][NXMEM][NYMEM];
extern double par_lay[NZ][NXMEM][NYMEM];
extern double sfc_swr[NXMEM][NYMEM];
extern double ***lightlim;
#  ifdef NITRATE
extern double ***nitrlim;
#  endif
extern double ***ironlim;
# endif /* PROGNOSTIC */
# ifdef NITRATE
extern double no3_lev[NZPHOS][NXMEM][NYMEM];
extern double no3_lay[NZ][NXMEM][NYMEM];
# endif
extern double dt;
/*
extern double jpo4[NZ][NXMEM][NYMEM];
extern double jdop[NZ][NXMEM][NYMEM];
extern double jremdop[NZ][NXMEM][NYMEM];
extern double jprod[NZ][NXMEM][NYMEM];
extern double jremin[NZ][NXMEM][NYMEM];
extern double flux_pop[NXMEM][NYMEM];
*/
#ifdef OXYGEN
extern double jo2[NZ][NXMEM][NYMEM];
#endif
#ifdef OXY18
extern double jo18[NZ][NXMEM][NYMEM];
extern double Ro18o16[NZ][NXMEM][NYMEM];
#endif
#ifdef NITRATE
extern double jno3[NZ][NXMEM][NYMEM];
extern double jdon[NZ][NXMEM][NYMEM];
extern double jremdon[NZ][NXMEM][NYMEM];
extern double jprod_n[NZ][NXMEM][NYMEM];
extern double jn2fix[NZ][NXMEM][NYMEM];
extern double jremin_n[NZ][NXMEM][NYMEM];
extern double jdenitw[NZ][NXMEM][NYMEM];
extern double jdenits[NZ][NXMEM][NYMEM];
extern double flux_pon[NXMEM][NYMEM];
# ifdef N15_CYCLE
extern double jno3_15n[NZ][NXMEM][NYMEM];
extern double jdon_15n[NZ][NXMEM][NYMEM];
extern double jremdon_15n[NZ][NXMEM][NYMEM];
extern double jprod_15n[NZ][NXMEM][NYMEM];
extern double jremin_15n[NZ][NXMEM][NYMEM];
extern double flux_pon_15n[NXMEM][NYMEM];
# endif
#endif /* NITRATE */
#ifdef DIC
extern double ***jdic;
extern double ***jalk;
#endif

extern int mPHOSPHATE;
extern int mDOP;

#ifdef OXYGEN
extern int mOXYGEN;
#endif

#ifdef OXY18
extern int mo18;
#endif

#ifdef NITRATE
extern int mNITRATE;
extern int mDON;
# ifdef N15_CYCLE
extern int mNITRATE_15n;
extern int mDON_15n;
# endif
#endif /* NITRATE */

void biotic_sms(int ibiodt)
{
    int i, j, k, l, kjunk;
    int kcomp, kmax; 
    double po4obsprof[NZPHOS], junk;
# ifdef PROGNOSTIC
    double feobsprof[NZPHOS];
# endif
    double compensation_depth, martin_coeff, fracprod;
    double phi_p, kappa_p, phi_ca;
    extern const double r_bio_tau;
    double c_2_p, caco3_2_c, n_2_p;
    double depth[NZ];
    double dtop[NZ], dbot[NZ];
    double fp[NZ], zremin[NZ];
    double fca[NZ], zreminca[NZ];
    double jdop_ij[NZ], jpo4_ij[NZ], jca[NZ];
    double jdic_ij[NZ], jalk_ij[NZ], jprod_ij[NZ];
    double jremin_ij[NZ], jremdop_ij[NZ];
    double flux_caco3;
    double imbal, flux_rem;
    double D_ij, dt_bio, frac_dt_bio;
    double po4[NZ], dop[NZ]; 
    double flux_pop_ij, flux_sed_ij, fclog, fdlog;
# ifdef NITRATE
    const double r_N_P = 16.0;
    const double r_N_P_diaz = r_N_P;
    const double r_N_P_denit = 104.0;
    double phi_n, kappa_n;
    double jno3_ij[NZ], jdon_ij[NZ], jprodN_ij[NZ];
    double jn2fix_ij[NZ], jdenitw_ij[NZ], jdenits_ij[NZ], jreminN_ij[NZ], jremdon_ij[NZ];
    double fp_n[NZ];
    double no3[NZ], don[NZ];
    double flux_pon_ij;
#  ifdef N15_CYCLE
    double jno3_15n_ij[NZ], jdon_15n_ij[NZ], jprodN_15n_ij[NZ], jreminN_15n_ij[NZ], jremdon_15n_ij[NZ];
    double jdenitw_15n_ij[NZ], jdenits_15n_ij[NZ], jn2fix_15n_ij[NZ], Rn15n14[NZ];
    double no3_15n[NZ], don_15n[NZ];
    double flux_pon_15n_ij, fp_15n[NZ];
#  endif
# endif /* NITRATE */
#ifdef OXYGEN
    double o2[NZ], jo2_ij[NZ]; 
    const double o2_min = 5.e-3; // -1000.;
    const double o_2_p = -150.;
#endif
#ifdef OXY18
    double o18[NZ], jo18_ij[NZ], Ro18o16[NZ];
    const double Rsmow = 0.0020052; 
    const double alpha_resp = 0.982; 
    //    ratio_sw = 0.97704;
    //ratio_o2_sw = 1.00075;
#endif

    int nzlevitus = NZPHOS;
    double levitus_depths[NZPHOS] = {0, 10, 20, 30, 50, 75, 100,
				    120, 150, 200, 250, 300, 400, 500, 600,
				    700, 800, 900, 1000, 1100, 1200, 1300, 
				    1400, 1500, 1750, 2000, 2500, 3000, 
				    3500, 4000, 4500, 5000, 5500};
# ifdef PROGNOSTIC
    double light_lim[NZ];
    double t_lim[NZ];
    double fe_lim[NZ];
    double po4_lim[NZ];
#  ifdef NITRATE 
    double no3_lim[NZ];
#  endif
    extern double lambda0; // mol/(m3 s)
    double lambda[NZ];
    const double par_frac = 0.45; // as in BEC model
    const double k_sw = 0.04; // 1/m
    extern double k_I; // W/m2
    extern double par_b; // 1/dgC
    extern double k_PO4; // mol/m3
#  if defined NITRATE && defined CALC_R_N_P
    const double k_NO3 = 0.5e-3; // mol/m3
#  endif
#  ifdef NITRATE
    const double k_NO3_diaz = 5.e-3; // mol/m3
    const double gamma_diaz = 0.1;
#  endif
    extern double k_Fe; // mol/m3
# endif /* PROGNOSTIC */
//BX # ifdef NITRATE
    const double eps_lim = 1.e-9; // since (po4/dt)*dt is NOT exactly == po4!!
//BX # endif

    phi_p = 0.67;
    kappa_p = 1.0 / (1.0 * 360.0 * 86400.0);
    phi_ca = phi_p;
# ifdef NITRATE
    phi_n = phi_p; // 0.67;
    kappa_n = kappa_p; //1.0 / (1.0 * 360.0 * 86400.0);
# endif
    martin_coeff = 0.9;
    //r_bio_tau    = 1.0 / (30.0 * 86400.0);    // restoring time step
    compensation_depth = 75.0;
    n_2_p     = 16.0;
    caco3_2_c = 0.07; 
    c_2_p     = 117.0;

    dt_bio = dt / (double) ibiodt;
    frac_dt_bio = 1.0 / (double) ibiodt;

//    printf("conc_obs_layer(h,po4_star_lev,po4_star_lay)\n");
//    conc_obs_layer(h,po4_star_lev,po4_star_lay);
# ifdef PROGNOSTIC
#  ifdef NITRATE
    // nitrate
    printf("conc_obs_layer(h,no3_lev,no3_lay)\n");
    conc_obs_layer(h,no3_lev,no3_lay);
#  endif

    // iron
    printf("conc_obs_layer(h,fe_lev,fe_lay)\n");
    conc_obs_layer(h,fe_lev,fe_lay);
# endif /* PROGNOSTIC */

// the outer loops
    for (i=X1;i<=nx;i++) {
	for (j=Y1;j<=ny;j++) {
	    double hsum;

/*-----------------------------------------------------------------------
 *  Compute depth at middle, top and bottom of layers (depth, dtop, dbot)
 *  Find index of compensation depth (kcomp) 
 *  Find index of deepest layer with h>5m (kmax)
 *-----------------------------------------------------------------------*/

	    for (k=0;k<NZ;k++) {
		depth[k] = 0.0;
	    }
// compute depth at layer k as the depth of the point in the 
//  middle of that layer
	    D_ij = D[i][j]; // create local variable
	    if (D_ij>MINIMUM_DEPTH) {
		hsum = h[0][i][j];
		depth[0] = hsum * 0.5;
		for (k=1;k<NZ;k++) {
		    depth[k] = hsum + 0.5 * h[k][i][j]; 
		    hsum += h[k][i][j];
		}
	    } else { 
		for (k=0;k<NZ;k++) {
		    depth[k] = misval;
		}
	    }
	    for (k=0;k<NZ;k++) {
		dtop[k] = depth[k] - 0.5 * h[k][i][j];
		dbot[k] = depth[k] + 0.5 * h[k][i][j];
	    }

	    kmax = NZ - 1;
	    kjunk = NZ - 1;
	    while (h[kjunk][i][j] < 5.0 && kjunk > 0) {
		kmax = --kjunk;
	    }
	  
	    kcomp = 0;
	    kjunk = 0;
	    while (dtop[kjunk] < compensation_depth && 
		   kjunk <= kmax) {
		kcomp = kjunk++;
	    }
	  
/*-----------------------------------------------------------------------
 *     calculate remin curves for OM and CaCO3
 *-----------------------------------------------------------------------*/

	    for (k=0;k<NZ;k++) {
		if (k < kcomp) {
		    zremin[k] = 0.0;
		    zreminca[k] = 0.0;
		} else if (k >= kcomp && k < kmax) {
		    junk = (dbot[k] / compensation_depth); 
		    zremin[k] = pow(junk,-martin_coeff);
		    zreminca[k] = exp(-(dbot[k]-compensation_depth)/3500.0);
		} else if (k >= kmax) {
		    zremin[k] = 0.0;
		    zreminca[k] = 0.0;
		}
	    }

/*---------------------------------------------------------------------
 *     copy tracers to local arrays
 *---------------------------------------------------------------------*/

	    for (k=0;k<NZ;k++) {
		po4[k] = tr[mPHOSPHATE][k][i][j];
		dop[k] = tr[mDOP][k][i][j];
# ifdef NITRATE
		no3[k] = tr[mNITRATE][k][i][j];
		don[k] = tr[mDON][k][i][j];
#  ifdef N15_CYCLE
		no3_15n[k] = tr[mNITRATE_15n][k][i][j];
		don_15n[k] = tr[mDON_15n][k][i][j];
		if (no3[k] > 0.0)
		  Rn15n14[k] = no3_15n[k] / no3[k];
		else
		  Rn15n14[k] = 0.0;
#  endif
# endif /* NITRATE */
# ifdef OXYGEN
		o2[k] = tr[mOXYGEN][k][i][j];
# endif
#  ifdef OXY18
		o18[k] = tr[mO18][k][i][j];
		Ro18o16[k] = o18[k] / o2[k];
#  endif

	    }

/* ---------------------------------------------------
 *  loop over ibiodt time steps
 * -------------------------------------------------*/

	    for (l=1; l<=ibiodt; l++) { 

/*-----------------------------------------------------------------------
 *     Production
 *-----------------------------------------------------------------------*/

# ifdef PROGNOSTIC
		if (D_ij > compensation_depth) {
		  // light
		  for (k=0; k <= kcomp && k < NZ; k++) {
		    par_lay[k][i][j] = par_frac * sfc_swr[i][j] /
		      ((dbot[k] - dtop[k]) * k_sw) *
		      (exp(-k_sw * dtop[k]) - exp(-k_sw * dbot[k]));
		    
		    light_lim[k] = par_lay[k][i][j] / 
		      (par_lay[k][i][j] + k_I);
		    
		    t_lim[k] = exp(par_b * Temptm[k][i][j]);

		    if (po4[k] < 0.0) {
		      po4_lim[k] = 0.0;
		    } else
		      po4_lim[k] = po4[k] / (po4[k] + k_PO4);


		    fe_lim[k] = fe_lay[k][i][j] /
		      (fe_lay[k][i][j] + k_Fe);

		    lambda[k] = lambda0 * t_lim[k] * light_lim[k] * fe_lim[k];
		    jprod_ij[k] = lambda[k]  * po4_lim[k];

		    // limit production to available PO4
		    if (jprod_ij[k] > 0.0 && jprod_ij[k] > ((po4[k]-eps_lim)/dt) ) {
		      if ( (po4[k]-eps_lim)/dt > 0.0)
			jprod_ij[k] = (po4[k]-eps_lim)/dt;
		      else // prevent jprod from being negative
			jprod_ij[k] = 0.0;
		    }
		  } // end of loop over k
		} else { // D_ij <= compensation_depth
		  for (k=0; k<NZ; k++) {
		    jprod_ij[k] = 0.0;
		    light_lim[k] = 0.0;
		    po4_lim[k] = 0.0;
		    fe_lim[k] = 0.0;
		  }		  
		}

# else /* not PROGNOSTIC, i.e., diagnostic model */

		for (k = 0; k <= kcomp && k < NZ; k++) {
		    if (po4[k] > po4_star_lay[k][i][j] && 
			po4_star_lay[k][i][j] >= 0.0 &&
			D_ij > compensation_depth)   {
			jprod_ij[k] = (po4[k] - po4_star_lay[k][i][j]) * r_bio_tau;
		    } else
			jprod_ij[k] = 0.0;
		}

# endif /* PROGNOSTIC */

		for (k = kcomp+1; k < NZ; k++) {
		    jprod_ij[k] = 0.0;
# ifdef PROGNOSTIC
		    light_lim[k] = 0.0;
		    po4_lim[k] = 0.0;		    
		    fe_lim[k] = 0.0;		    
# endif /* PROGNOSTIC */
		}

# ifdef NITRATE
#  ifdef PROGNOSTIC
		if (D_ij > compensation_depth) {
		  for (k = 0; k <= kcomp && k < NZ; k++) {
#   ifdef SPEC_R_N_P
		    no3_lim[k] = 1.0; // set mostly for output
		    
		    jprodN_ij[k] = jprod_ij[k] * r_N_P;
		    if (jprodN_ij[k] > 0.0 && jprodN_ij[k] > ((no3[k]-eps_lim)/dt) ) {
		      if ( (no3[k]-eps_lim)/dt > 0.0)
			jprodN_ij[k] = (no3[k]-eps_lim)/dt;
		      else // prevent jprodN from being negative
			jprodN_ij[k] = 0.0;
		      // adjust PO4 drawdown accordingly
		      jprod_ij[k] = jprodN_ij[k] / r_N_P;
		    }
#   endif /* SPEC_R_N_P */
#   ifdef CALC_R_N_P
		    if (no3[k] > 0.0)
		      no3_lim[k] = no3[k] / (no3[k] + k_NO3);
		    else
		      no3_lim[k] = 0.0;
		    jprodN_ij[k] = lambda[k] * no3_lim[k];
		    // limit production to available NO3
		    if (jprodN_ij[k] > 0.0 && jprodN_ij[k] > ((no3[k]-eps_lim)/dt) ) {
		      if ( (no3[k]-eps_lim)/dt > 0.0)
			jprodN_ij[k] = (no3[k]-eps_lim)/dt;
		      else // prevent jprodN from being negative
			jprodN_ij[k] = 0.0;
		    }
#   endif /* CALC_R_N_P */
		    // N2 fixation
		    jn2fix_ij[k] = lambda[k] * po4_lim[k] * 
		      (1.0 - no3[k]/(no3[k] + k_NO3_diaz)) *
		      gamma_diaz * r_N_P_diaz;
		  } // end of loop over k
		} else { // D_ij <= compensation_depth
		  for (k=0; k<NZ; k++) {
		    jprodN_ij[k] = 0.0;
		    no3_lim[k] = 0.0;
		    jn2fix_ij[k] = 0.0;
#   ifdef N15_CYCLE
		    jprodN_15n_ij[k] = 0.0;
#   endif
		  }		  
		}
#  else /* diagnostic model */
		for (k = 0; k <= kcomp && k < NZ; k++) {
		  double prod_diff;
		  if (no3[k] > no3_lay[k][i][j] && no3_lay[k][i][j] >= 0.0 &&
			D_ij > compensation_depth) {
		      jprodN_ij[k] =  (no3[k] - no3_lay[k][i][j]) * r_bio_tau;
		    } else
		      jprodN_ij[k] = 0.0;

		  // N2 fixation
		  prod_diff = jprod_ij[k] - jprodN_ij[k] / r_N_P;
		  //		  jn2fix_ij[k] = prod_diff > 0.0 ? prod_diff * r_N_P_diaz : 0.0;
		  jn2fix_ij[k] =  Temptm[k][i][j] > 25.0 ? 100e12/14/3.15e7/1.4e16 : 0.0;
		}
#  endif /* PROGNOSTIC */
#  ifdef N15_CYCLE
		for (k = 0; k < NZ; k++) {
		  if (no3[k] > 0.0 && no3_15n[k] > 0.0)
		    jprodN_15n_ij[k] = jprodN_ij[k] * Rn15n14[k] * parm_alpha_n15_prod;
		  else
		    jprodN_15n_ij[k] = 0.0;
		}
#  endif
		for (k = kcomp+1; k < NZ; k++) {
		    jprodN_ij[k] = 0.0;
		    jn2fix_ij[k] = 0.0;
#  ifdef N15_CYCLE
		    jprodN_15n_ij[k] = 0.0;
#  endif
		}		    
# endif /* NITRATE */

		if (kcomp >= 0) {
		    fracprod = (compensation_depth - dtop[kcomp]) /
			h[kcomp][i][j];
		    jprod_ij[kcomp] *= fracprod;
# ifdef NITRATE		    
		    jprodN_ij[kcomp] *= fracprod;
		    jn2fix_ij[k] *= fracprod;
#  ifdef N15_CYCLE
		    jprodN_15n_ij[kcomp] *= fracprod;
#  endif
# endif /* NITRATE */
		}

/*-----------------------------------------------------------------------
 *     Sinking Particle Flux
 *-----------------------------------------------------------------------*/

/*       calculate the sinking flux at the compensation depth */
		flux_pop_ij = (1.0 - phi_p) * jprod_ij[0] * h[0][i][j];
		for (k=1; k<NZ; k++) {
		    flux_pop_ij += (1.0 - phi_p) * jprod_ij[k] * h[k][i][j];
		}
/*       distribute sinking flux with depth */
		for (k=0;k<NZ;k++)
		    fp[k] = flux_pop_ij * zremin[k];
# ifdef NITRATE
		flux_pon_ij = (1.0 - phi_n) * jprodN_ij[0] * h[0][i][j];
		for (k=1; k<NZ; k++) {
		    flux_pon_ij += (1.0 - phi_n) * jprodN_ij[k] * h[k][i][j];
		}
		for (k=0;k<NZ;k++)
		    fp_n[k] = flux_pon_ij * zremin[k];
#  ifdef N15_CYCLE
		flux_pon_15n_ij = (1.0 - phi_n) * jprodN_15n_ij[0] * h[0][i][j];
		for (k=1; k<NZ; k++) {
		    flux_pon_15n_ij += (1.0 - phi_n) * jprodN_15n_ij[k] * h[k][i][j];
		}
		for (k=0;k<NZ;k++)
		    fp_15n[k] = flux_pon_15n_ij * zremin[k];
#  endif /* N15_CYCLE */
# endif /* NITRATE */
	    
/*-----------------------------------------------------------------------
 *     Particle Remineralization
 *-----------------------------------------------------------------------*/

		for (k=0;k<NZ;k++) {
		    if (k < kcomp) 
			jremin_ij[k] = 0.0;
		    else if (k == kcomp) 
			jremin_ij[k] = (flux_pop_ij - fp[k]) / h[k][i][j];
		    else // no need to check that k > kcomp
			jremin_ij[k] = (fp[k-1] - fp[k]) / h[k][i][j];

		    //CAD - Are these checks really needed?
		    /* prevent large concentrations in layers of negligible thickness */
		    if (h[k][i][j] < 1.0e-9)
			jremin_ij[k] = 0.0;
		}
	    
/*       bottom boundary condition for shallow grid points (D < 75m) */

		if (D_ij < compensation_depth && D[i][j]>MINIMUM_DEPTH) {
		    jremin_ij[kmax] = flux_pop_ij / h[kmax][i][j];
		    // printf("D < comp_depth at grid point %i,%i \n",i,j);
		}
# ifdef NITRATE
		for (k = 0; k < NZ; k++) {
		    if (k < kcomp) 
			jreminN_ij[k] = 0.0;
		    else if (k == kcomp) 
			jreminN_ij[k] = (flux_pon_ij - fp_n[k]) / h[k][i][j];
		    else // no need to check that k > kcomp
			jreminN_ij[k] = (fp_n[k-1] - fp_n[k]) / h[k][i][j];

		    /* prevent large concentrations in layers of negligible thickness */
		    if (h[k][i][j] < 1.0e-9)
			jreminN_ij[k] = 0.0;
		}
	    
/*       bottom boundary condition for shallow grid points (D < 75m) */

		if (D_ij < compensation_depth && D_ij > MINIMUM_DEPTH)
		    jreminN_ij[kmax] = flux_pon_ij / h[kmax][i][j];

		/* water column denitrification */
		for (k=0;k<NZ;k++) {
		  if (o2[k] < o2_min) {
		    jdenitw_ij[k] = jremin_ij[k] * r_N_P_denit;
		    if (jdenitw_ij[k] > ((no3[k]-eps_lim)/dt)-jprodN_ij[k] ) {
		      if ( (no3[k]-eps_lim)/dt > 0.0)
			jdenitw_ij[k] = (no3[k]-eps_lim)/dt - jprodN_ij[k];
		      else
			jdenitw_ij[k] = 0.0;
		    }
		  } else
		    jdenitw_ij[k] = 0.0;
		}

		/* sediment denitrification */
		for (k=0;k<NZ;k++) {
		  jdenits_ij[k] = 0.0;
		}		    
		if (D_ij > compensation_depth && flux_pop_ij > 0.0) {
		  junk = dbot[kmax] / compensation_depth; 
		  flux_sed_ij = flux_pop_ij * pow(junk,-martin_coeff);
		  fclog = log(flux_sed_ij * c_2_p * 1.e2 * 86400.0); // log of Corg flux in umol/cm2/d
		  fdlog = -0.9543 + 0.7662 * fclog - 0.2350 * pow(fclog,2.0); // log of sed dn flux Middelburg
		  jdenits_ij[kmax] = exp(fdlog) / (1.e2 * 86400.0 * h[kmax][i][j]);  // denit rate uM/sec
		  if (jdenits_ij[kmax] > 0.0 && 
		      jdenits_ij[kmax] > 
		      ((no3[kmax]-eps_lim)/dt - jdenitw_ij[kmax] - jprodN_ij[kmax]) ) {
		    if ( (no3[kmax]-eps_lim)/dt > 0.0)
			jdenits_ij[kmax] = (no3[kmax]-eps_lim)/dt - jdenitw_ij[kmax] - jprodN_ij[kmax];
		      else
			jdenits_ij[k] = 0.0;
		  }
		}
		

#  ifdef N15_CYCLE
		for (k=0;k<NZ;k++) {
		    if (k < kcomp) 
			jreminN_15n_ij[k] = 0.0;
		    else if (k == kcomp) 
			jreminN_15n_ij[k] = (flux_pon_15n_ij - fp_15n[k]) / h[k][i][j] *
			  parm_alpha_n15_pon_remin;
		    else // no need to check that k > kcomp
			jreminN_15n_ij[k] = (fp_15n[k-1] - fp_15n[k]) / h[k][i][j] *
			  parm_alpha_n15_pon_remin;

		    /* prevent large concentrations in layers of negligible thickness */
		    if (h[k][i][j] < 1.0e-9)
			jreminN_15n_ij[k] = 0.0;
		}
	    
/*       bottom boundary condition for shallow grid points (D < 75m) */

		if (D_ij < compensation_depth && D_ij > MINIMUM_DEPTH) {
		    jreminN_15n_ij[kmax] = flux_pon_15n_ij / h[kmax][i][j] *
		      parm_alpha_n15_pon_remin;  // CAD - Will this conserve 15N mass?
		}

		/* water column denitrification */
		for (k=0;k<NZ;k++) {
		    jdenitw_15n_ij[k] = jdenitw_ij[k] * Rn15n14[k] * parm_alpha_n15_denitw;
		}

		/* sediment denitrification */
		for (k=0;k<NZ;k++) {
		  jdenits_15n_ij[k] = 0.0;
		}
		jdenits_15n_ij[kmax] = jdenits_ij[kmax] * Rn15n14[k] * parm_alpha_n15_denits;

		/* N2 fixation */
		for (k=0;k<NZ;k++) {
		    jn2fix_15n_ij[k] = jn2fix_ij[k] * parm_n15_std_fraction *
		      parm_alpha_n15_n2fix;
		}

#  endif /* N15_CYCLE */
# endif /* NITRATE */

/*-----------------------------------------------------------------------
 *     PO4 and DOP
 *-----------------------------------------------------------------------*/
		for (k=0;k<NZ;k++) {
		  jremdop_ij[k] = kappa_p * dop[k];
		  jpo4_ij[k] = -jprod_ij[k] + jremin_ij[k] + jremdop_ij[k];
		  jdop_ij[k] = phi_p * jprod_ij[k] - jremdop_ij[k];
		}

/*-----------------------------------------------------------------------
 *     NO3 and DON
 *-----------------------------------------------------------------------*/
# ifdef NITRATE

		for (k=0;k<NZ;k++) {
		  jremdon_ij[k] = kappa_n * don[k];
		  jno3_ij[k] = -jprodN_ij[k] + 
		    jreminN_ij[k] + jremdon_ij[k] + 
		    jn2fix_ij[k] - jdenitw_ij[k] - jdenits_ij[k];
		  jdon_ij[k] = phi_n * jprodN_ij[k] - jremdon_ij[k];
		}

/*-----------------------------------------------------------------------
 *     NO3_15N and DON_15N
 *-----------------------------------------------------------------------*/

# ifdef N15_CYCLE

		for (k=0;k<NZ;k++) {
		  if (don[k] > 0.0 && don_15n[k] > 0.0)
		    jremdon_15n_ij[k] = jremdon_ij[k] * don_15n[k] / don[k] *
		      parm_alpha_n15_don_remin;
		  else
		    jremdon_15n_ij[k] = 0.0;

		  jdon_15n_ij[k] = phi_n * jprodN_15n_ij[k] - jremdon_15n_ij[k];
		  jno3_15n_ij[k] = -jprodN_15n_ij[k] + 
		    jreminN_15n_ij[k] + jremdon_15n_ij[k] +
		    jn2fix_15n_ij[k] - jdenitw_15n_ij[k] - jdenits_15n_ij[k];
		}

# endif /* N15_CYCLE */

#endif /* NITRATE */

/*-----------------------------------------------------------------------
 *     O2 / O18
 *-----------------------------------------------------------------------*/

# ifdef OXYGEN
		for (k=0;k<NZ;k++) {
		  if (o2[k] < o2_min && jpo4_ij[k] > 0.0) {
		    jo2_ij[k] = 0.0;
		  } else {
		    jo2_ij[k] = o_2_p * jpo4_ij[k];
		  }
		}

# endif /* OXYGEN */


#  ifdef OXY18
		for (k=0;k<NZ;k++) {
		  if (o2[k] < o2_min && jpo4_ij[k] > 0.0) {
		    jo18_ij[k] = 0.0;
		  } else {
		    jo18_ij[k] = -o_2_p * (jprod_ij[k] * Rsmow -
					(jremin_ij[k] + jremdop_ij[k]) *
					Ro18o16_ij[k] * alpha_resp);
		  }
		}

#  endif

/*-----------------------------------------------------------------------
 *     CaCO3
 *-----------------------------------------------------------------------*/

/*  calculate the flux of CaCO3 from the compensation depth to bottom*/

		flux_caco3 = caco3_2_c * c_2_p * flux_pop_ij;
		for (k=0; k<NZ; k++) {
		    fca[k] = flux_caco3 * zreminca[k];
		}

/*  calculate formation/dissolution of CaCO3 above/below compensation depth */
		for (k=0; k<NZ; k++) {
		    if (k < kcomp) 
			jca[k] = -caco3_2_c * c_2_p * (1.0 - phi_ca) * jprod_ij[k];
		    else if (k == kcomp) 
			jca[k] = -caco3_2_c * c_2_p * (1.0 - phi_ca) * jprod_ij[k] + 
			    (flux_caco3 - fca[k]) / h[k][i][j];
		    else // no need to check that k > kcomp
			jca[k] = (fca[k-1] - fca[k]) / h[k][i][j];
		    /* prevent large concentrations in layers of negligible thickness */
		    if (h[k][i][j] < 1.0e-9)
			jca[k] = 0.0;

		}

/*       bottom boundary condition for shallow grid points (D < 75m) */
		if (D_ij < compensation_depth && D[i][j]>MINIMUM_DEPTH) {
		    jca[kmax] = flux_caco3 / h[kmax][i][j];
		    // printf("D < comp_depth at grid point %i,%i \n",i,j);
		}

/*-----------------------------------------------------------------------
 *     DIC
 *-----------------------------------------------------------------------*/
		for (k=0; k<NZ; k++)
		    jdic_ij[k] = c_2_p * jpo4_ij[k] + jca[k];

/*-----------------------------------------------------------------------
 *     ALK
 *-----------------------------------------------------------------------*/
		for (k=0; k<NZ; k++)
		    jalk_ij[k] = -n_2_p * jpo4_ij[k] + 2.0 * jca[k];


/*-----------------------------------------------------------------------
 *       check for mass conservation in sinking/remin fluxes
 *-----------------------------------------------------------------------*/


		flux_rem = jremin_ij[0] * h[0][i][j];
		for (k=1; k<NZ; k++) {
		    flux_rem += jremin_ij[k] * h[k][i][j];
		}

		imbal = flux_rem - flux_pop_ij;
		if (flux_pop_ij < 0.0)
		  printf("flux_pop[%d][%d] = %g\n", i, j, k, flux_pop_ij);
		if (fabs(imbal) > flux_pop_ij*1e-4) { /* error tolerance = 0.01% */
		    printf("mass imbalance (P) at grid point %i,%i \n",i,j);
		    printf("flux_pop,flux_rem,imbal,D,kcomp,kmax = %g,%g,%g,%g,%i,%i \n",
			   flux_pop_ij,flux_rem,imbal,D_ij,
			   kcomp,kmax);
		}

#ifdef NITRATE

		flux_rem = jreminN_ij[0] * h[0][i][j];
		for (k=1; k<NZ; k++) {
		    flux_rem += jreminN_ij[k] * h[k][i][j];
		}

		imbal = flux_rem - flux_pon_ij;
		if (flux_pon_ij < 0.0)
		  printf("flux_pon[%d][%d] = %g\n", i, j, flux_pon_ij);
		if (fabs(imbal) > flux_pon_ij*1e-4) { /* error tolerance = 0.01% */
		    printf("mass imbalance (N) at grid point %i,%i \n",i,j);
		    printf("flux_pon,flux_rem,imbal,D,kcomp,kmax = %g,%g,%g,%g,%i,%i \n",
			   flux_pon_ij,flux_rem,imbal,D_ij,
			   kcomp,kmax);
		}
# ifdef N15_CYCLE
		flux_rem = jreminN_15n_ij[0] * h[0][i][j];
		for (k=1; k<NZ; k++) {
		    flux_rem += jreminN_15n_ij[k] * h[k][i][j];
		}

		imbal = flux_rem - flux_pon_15n_ij;
		if (flux_pon_15n_ij < 0.0)
		  printf("flux_pon_15n[%d][%d] = %g\n", i, j, flux_pon_15n_ij);
		if (fabs(imbal) > flux_pon_15n_ij*1e-4) { /* error tolerance = 0.01% */
		    printf("mass imbalance (15N) at grid point %i,%i \n",i,j);
		    printf("flux_pon_15n,flux_rem,imbal,D,kcomp,kmax = %g,%g,%g,%g,%i,%i \n",
			   flux_pon_15n_ij,flux_rem,imbal,D_ij,
			   kcomp,kmax);
		}
# endif /* N15_CYCLE */
# endif /* NITRATE */

/*-----------------------------------------------------------------------
 *     save fluxes to 3d fields
 *-----------------------------------------------------------------------*/

		for (k=0; k<NZ; k++) {
		  po4[k]  += dt_bio * jpo4_ij[k];  //CAD - why are these here? 
		    dop[k]  += dt_bio * jdop_ij[k];

		    // copy fluxes from local 1D to global 3D variables
		    jpo4[k][i][j]    += (jpo4_ij[k] * frac_dt_bio);
		    jdop[k][i][j]    += jdop_ij[k] * frac_dt_bio;
		    jremdop[k][i][j] += jremdop_ij[k] * frac_dt_bio;
		    jprod[k][i][j]   += jprod_ij[k] * frac_dt_bio;
		    jremin[k][i][j]  += jremin_ij[k] * frac_dt_bio;
# ifdef PROGNOSTIC
		    lightlim[k][i][j]  += light_lim[k] * frac_dt_bio;
#  ifdef NITRATE
		    nitrlim[k][i][j]  += no3_lim[k] * frac_dt_bio;
#  endif
	 	    ironlim[k][i][j]  += fe_lim[k] * frac_dt_bio;
# endif /* PROGNOSTIC */
# ifdef NITRATE
		    jno3[k][i][j]    += jno3_ij[k] * frac_dt_bio;
		    jn2fix[k][i][j]  += jn2fix_ij[k] * frac_dt_bio;
		    jdenitw[k][i][j] += jdenitw_ij[k] * frac_dt_bio;
		    jdenits[k][i][j] += jdenits_ij[k] * frac_dt_bio;
		    jdon[k][i][j]    += jdon_ij[k] * frac_dt_bio;
		    jremdon[k][i][j] += jremdon_ij[k] * frac_dt_bio;
#  ifdef N15_CYCLE
		    jno3_15n[k][i][j]    += jno3_15n_ij[k] * frac_dt_bio;
		    jdon_15n[k][i][j]    += jdon_15n_ij[k] * frac_dt_bio;
		    jremdon_15n[k][i][j] += jremdon_15n_ij[k] * frac_dt_bio;
#  endif
# endif /* NITRATE */
# ifdef DIC
		    jdic[k][i][j]    += jdic_ij[k] * frac_dt_bio;
		    jalk[k][i][j]    += jalk_ij[k] * frac_dt_bio;
# endif
# ifdef OXYGEN
		    jo2[k][i][j]    += (jo2_ij[k] * frac_dt_bio);
# endif
# ifdef OXY18
		    jo18[k][i][j]    += jo18_ij[k] * frac_dt_bio;
# endif

		}
//BX added if loop - only positive statement before
//BX   not used any more if mn_fields are set to misval in offtrac.c
//BX		if (D_ij > MINIMUM_DEPTH) {
		flux_pop[i][j] += flux_pop_ij * frac_dt_bio;
//BX		} else {
//BX		    flux_pop[i][j] = misval;
//BX		}
//BX-a
/*       bottom boundary condition for shallow grid points (D < 75m) */
//BX		if (D_ij < compensation_depth && D_ij > MINIMUM_DEPTH) {
//BX		  for (k=0; k<NZ; k++) {jalk[k][i][j] = 0;}
//BX		    //printf("D < comp_depth at grid point %i,%i\n",i,j);
//BX		}
//BX-e

# ifdef NITRATE
		flux_pon[i][j] += flux_pon_ij * frac_dt_bio;
#  ifdef N15_CYCLE
		flux_pon_15n[i][j] += flux_pon_15n_ij * frac_dt_bio;
#  endif
# endif /* NITRATE */



	    } // end of loop over l
	} // end of loop over j
    } // end of loop over i
} // end of subroutine biotic_sms

# endif /* PHOSPHATE */
