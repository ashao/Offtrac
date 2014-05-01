double lin_interp(double pleth, const double x[], const double z[], int istart,
		int npts);

double lin_interpp(double pleth, const double x[], const double z[],
		int istart, int npts);



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

