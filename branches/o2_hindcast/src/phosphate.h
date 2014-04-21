/*
 * phosphate.h
 *
 *  Created on: Apr 15, 2014
 *      Author: ashao
 */

extern int mPHOSPHATE;
extern int mDOP;
// Output arrays
extern double ***mn_phos;
extern double ***mn_dop;
extern double ***mn_pobs;
extern double ***mn_jpo4;
extern double ***mn_jremin;
extern double ***mn_jremdop;
extern double mn_pflux[NXMEM][NYMEM];
extern double ***mn_jprod;
//
// // Working arrays
extern double ***phosphate_init;
extern double po4[NZ][NXMEM][NYMEM];
extern double jpo4[NZ][NXMEM][NYMEM];
extern double ***po4_star_lay;
extern double jprod[NZ][NXMEM][NYMEM];
extern double jremin[NZ][NXMEM][NYMEM];
extern double jremdop[NZ][NXMEM][NYMEM];
extern double jdop[NZ][NXMEM][NYMEM];
extern double flux_pop[NXMEM][NYMEM];
//

/* SUBROUTINE prototypes */
void initialize_phosphate( );
void apply_phosphate_jterms();

