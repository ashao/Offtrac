double lin_interp(double pleth, const double x[], const double z[], int istart,
		int npts);

double lin_interpp(double pleth, const double x[], const double z[],
		int istart, int npts);

void conc_obs_layer(double h[NZ][NXMEM][NYMEM],
		double conc_lev[NZPHOS][NXMEM][NYMEM],
		double conc_lay[NZ][NXMEM][NYMEM]);

void z_depth(double h[NZ][NXMEM][NYMEM], double depth[NZ][NXMEM][NYMEM]);

// begin ashao
/*     Begin added DT     */
/* REMOVED FOR MORE GENERAL ROUTINES (SEE BELOW)
 void schmidt_cfc(double t[NXMEM][NYMEM], int kn, double Sc_cfc[NXMEM][NYMEM]);
 void cfc_saturation(double T[NZ][NXMEM][NYMEM],double S[NZ][NXMEM][NYMEM],
 int method, double cfc_sat[NZ][NXMEM][NYMEM], int icfc);
 void cfc_atmospheric(int iyear, int icfc, double cfcatm[NXMEM][NYMEM]);
 /*      End added DT      */
// end ashao

#ifdef USE_CALC_H
void z_sum(double h[NZ][NXMEM][NYMEM],double tot_depth[NXMEM][NYMEM]);
#endif

//begin ashao
# if defined(SF6) || defined(CFC)
double trac_calcsol(double TempK, double S, double sol_coeffs[7]);
double trac_atmospheric( double thetime, int tracrecs,
		int tracyear[MAXRECS], double tracN[MAXRECS], double tracS[MAXRECS],
		int latidx);
double trac_schmidt(double T, double sc_coeffs[4]);
double trac_calcflux(double trac_Sc, double trac_atm, double trac_sol,
		double icefrac, double pistonvel, double pres, double trac_con, double dt);

# endif
//end ashao

