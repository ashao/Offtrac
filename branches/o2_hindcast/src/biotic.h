extern double ****tr;
extern double D[NXMEM][NYMEM];
extern double h[NZ][NXMEM][NYMEM];
extern double po4_star_lev[NZWOA][NXMEM][NYMEM];
extern double po4_star_lay[NZ][NXMEM][NYMEM];
extern double Temptm[NZ][NXMEM][NYMEM];
extern double dt;
extern double jpo4[NZ][NXMEM][NYMEM];
extern double jdop[NZ][NXMEM][NYMEM];
extern double jremdop[NZ][NXMEM][NYMEM];
extern double jprod[NZ][NXMEM][NYMEM];
extern double jremin[NZ][NXMEM][NYMEM];
extern double flux_pop[NXMEM][NYMEM];
#ifdef OXYGEN
extern double jo2[NZ][NXMEM][NYMEM];
#endif
extern int mPHOSPHATE;
extern int mDOP;

#ifdef OXYGEN
extern int mOXYGEN;
#endif


void biotic_sms(int ibiodt);
