#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include "init.h"
#include "io.h"
#include "alloc.h"
#include "par_IO.h"


int alloc_arrays()
{
    extern double ***u;
    extern double ***v;
    extern double ***mn_u;
    extern double ***mn_v;
    extern double ***mn_h;
    extern double ***uhtm;
    extern double ***vhtm;
    extern double ***ea;
    extern double ***eb;
    extern double ***mn_rml;
    extern double ***mn_uhtm;
    extern double ***mn_vhtm;
    extern double ***mn_ea;
    extern double ***mn_eb;
#ifdef INERT
    extern double ***mn_inert1;
    extern double ***mn_inert2;
#endif
#ifdef AGE
    extern double ***age;
    extern double ***mn_age;
# if defined (AGE2) || defined (AGE3)
    extern double ***hnew;
    extern double ***mn_hnew;
# endif
#endif
    /*#ifdef SMOOTH_REST
    extern double ***h_init;
    #endif*/
#ifdef OXYGEN
    extern double ***oxy_init;
    extern double ***mn_oxygen;
    extern double ***mn_o2sat;
    extern double ***mn_jo2;
#endif
#ifdef OXY18
    extern double ***mn_o18;
    extern double ***mn_jo18;
#endif

#ifdef CFC
    extern double ***mn_cfc11;
//    extern double ***mn_cfc11_relsat;
    extern double ***mn_cfc12;
//    extern double ***mn_cfc12_relsat;
#endif
// begin ashao
#ifdef SF6
	extern double ***mn_sf6;
	extern double ***mn_sf6_relsat;
#endif
#ifdef RADIOACTIVE
	extern double ***mn_i131;
	extern double ***mn_cs134;
	extern double ***mn_cs136;
	extern double ***mn_cs137;
	extern double ***mn_ba140;
	extern double ***mn_la140;
#endif
#ifdef TTD_BP
	extern double ***mn_ttd;
#endif
// end ashao

#ifdef PHOSPHATE
    extern double ***phosphate_init;
    extern double ***mn_phos;
    extern double ***mn_dop;
    extern double ***mn_pobs;
    extern double ***mn_jpo4;
# ifdef PROGNOSTIC
    extern double ***lightlim;
    extern double ***ironlim;
    extern double ***mn_lightlim;
    extern double ***mn_ironlim;
#  ifdef NITRATE
    extern double ***nitrlim;
    extern double ***mn_nitrlim;
#  endif
# endif /* PROGNOSTIC */
#endif /* PHOSPHATE */
#ifdef NITRATE
    extern double ***nitrate_init;
    extern double ***mn_nitr;
    extern double ***mn_don;
    extern double ***mn_nobs;
    extern double ***mn_jno3;
    extern double ***mn_n2fix;
    extern double ***mn_denitw;
    extern double ***mn_denits;
# ifdef N15_CYCLE
    extern double ***nitrate_init_15n;
    extern double ***mn_nitr_15n;
    extern double ***mn_don_15n;
    extern double ***mn_jno3_15n;
# endif /* N15_CYCLE */
#endif /* NITRATE */
#ifdef DIC
    extern double ***jdic;
    extern double ***mn_dic;
    extern double ***jalk;
    extern double ***mn_alk;
    extern double ***mn_jdic;
#endif
    if (! (u = alloc3d(NZ,NXMEM,NYMEM))) handle_error("u");
    if (! (v = alloc3d(NZ,NXMEM,NYMEM))) handle_error("v");
    if (! (mn_u = alloc3d(NZ,NXMEM,NYMEM))) handle_error("mn_u");
    if (! (mn_v = alloc3d(NZ,NXMEM,NYMEM))) handle_error("mn_v");
    if (! (mn_h = alloc3d(NZ,NXMEM,NYMEM))) handle_error("mn_h");
    uhtm = alloc3d(NZ,NXMEM,NYMEM);
    if(uhtm == NULL) {
	fprintf(stderr,"not enough memory for uhtm!\n");
	return 1;
    }
    vhtm = alloc3d(NZ,NXMEM,NYMEM);
    if(vhtm == NULL) {
	fprintf(stderr,"not enough memory for vhtm!\n");
	return 1;
    }
    ea = alloc3d(NZ,NXMEM,NYMEM);
    if(ea == NULL) {
	fprintf(stderr,"not enough memory for ea!\n");
	return 1;
    }    eb = alloc3d(NZ,NXMEM,NYMEM);
    if(eb == NULL) {
	fprintf(stderr,"not enough memory for eb!\n");
	return 1;
    }
    mn_rml = alloc3d(2,NXMEM,NYMEM);
    if(mn_rml == NULL) {
	fprintf(stderr,"not enough memory for mn_rml!\n");
	return 1;
    }
    mn_uhtm = alloc3d(NZ,NXMEM,NYMEM);
    if(mn_uhtm == NULL) {
	fprintf(stderr,"not enough memory for mn_uhtm!\n");
	return 1;
    }
    mn_vhtm = alloc3d(NZ,NXMEM,NYMEM);
    if(mn_vhtm == NULL) {
	fprintf(stderr,"not enough memory for mn_vhtm!\n");
	return 1;
    }
    mn_ea = alloc3d(NZ,NXMEM,NYMEM);
    if(mn_ea == NULL) {
	fprintf(stderr,"not enough memory for mn_ea!\n");
	return 1;
    }
    mn_eb = alloc3d(NZ,NXMEM,NYMEM);
    if(mn_eb == NULL) {
	fprintf(stderr,"not enough memory for mn_eb!\n");
	return 1;
    }
#ifdef INERT
    mn_inert1 = alloc3d(NZ,NXMEM,NYMEM);
    if(mn_inert1 == NULL) {
	fprintf(stderr,"not enough memory for mn_inert1!\n");
	return 1;
    }
    mn_inert2 = alloc3d(NZ,NXMEM,NYMEM);
    if(mn_inert2 == NULL) {
	fprintf(stderr,"not enough memory for mn_inert2!\n");
	return 1;
    }
#endif
    
#ifdef AGE
    age = alloc3d(NZ,NXMEM,NYMEM);
    if(age == NULL) {
	fprintf(stderr,"not enough memory for age!\n");
	return 1;
    }
    mn_age = alloc3d(NZ,NXMEM,NYMEM);
    if(mn_age == NULL) {
	fprintf(stderr,"not enough memory for mn_age!\n");
	return 1;
    }
# if defined (AGE2) || defined (AGE3)
    hnew = alloc3d(NZ,NXMEM,NYMEM);
    if(hnew == NULL) {
	fprintf(stderr,"not enough memory for hnew!\n");
	return 1;
    }
    mn_hnew = alloc3d(NZ,NXMEM,NYMEM);
    if(mn_hnew == NULL) {
	fprintf(stderr,"not enough memory for mn_hnew!\n");
	return 1;
    }
# endif
#endif /* AGE */

    /*#ifdef SMOOTH_REST
    h_init = alloc3d(NZ,NXMEM,NYMEM);
    if(h_init == NULL) {
	fprintf(stderr,"not enough memory for h_init!\n");
	return 1;
    }
    #endif*/

#ifdef OXYGEN
    oxy_init = alloc3d(NZ,NXMEM,NYMEM);
    if(oxy_init == NULL) {
	fprintf(stderr,"not enough memory for oxy_init!\n");
	return 1;
    }
    mn_oxygen = alloc3d(NZ,NXMEM,NYMEM);
    if(mn_oxygen == NULL) {
	fprintf(stderr,"not enough memory for mn_oxygen!\n");
	return 1;
    }
    mn_o2sat = alloc3d(NZ,NXMEM,NYMEM);
    if(mn_o2sat == NULL) {
	fprintf(stderr,"not enough memory for mn_o2sat!\n");
	return 1;
    }
    mn_jo2 = alloc3d(NZ,NXMEM,NYMEM);
    if(mn_jo2 == NULL) {
	fprintf(stderr,"not enough memory for mn_jo2!\n");
	return 1;
    }
#endif

#ifdef OXY18
    mn_o18 = alloc3d(NZ,NXMEM,NYMEM);
    if(mn_o18 == NULL) {
	fprintf(stderr,"not enough memory for mn_o18!\n");
	return 1;
    }
    mn_jo18 = alloc3d(NZ,NXMEM,NYMEM);
    if(mn_jo18 == NULL) {
	fprintf(stderr,"not enough memory for mn_jo18!\n");
	return 1;
    }
#endif
#ifdef CFC
    mn_cfc11 = alloc3d(NZ,NXMEM,NYMEM);
    if(mn_cfc11 == NULL) {
	fprintf(stderr,"not enough memory for mn_cfc11!\n");
	return 1;
    }
    mn_cfc12 = alloc3d(NZ,NXMEM,NYMEM);
    if(mn_cfc12 == NULL) {
	fprintf(stderr,"not enough memory for mn_cfc12!\n");
	return 1;
    }
/*
    mn_cfc11_relsat = alloc3d(NZ,NXMEM,NYMEM);
    if(mn_cfc11_relsat == NULL) {
	fprintf(stderr,"not enough memory for mn_cfc11!\n");
	return 1;
    }
    mn_cfc12_relsat = alloc3d(NZ,NXMEM,NYMEM);
    if(mn_cfc12_relsat == NULL) {
	fprintf(stderr,"not enough memory for mn_cfc12!\n");
	return 1;
    }
*/
#endif

// begin ashao
#ifdef SF6
    mn_sf6 = alloc3d(NZ,NXMEM,NYMEM);
    if(mn_sf6 == NULL) {
	fprintf(stderr,"not enough memory for mn_cfc11!\n");
	return 1;
    }
/*
    mn_sf6_relsat = alloc3d(NZ,NXMEM,NYMEM);
    if(mn_sf6_relsat == NULL) {
	fprintf(stderr,"not enough memory for mn_cfc11!\n");
	return 1;
    }
*/

#endif

#ifdef RADIOACTIVE

    mn_i131 = alloc3d(NZ,NXMEM,NYMEM);
    if(mn_i131 == NULL) {
	fprintf(stderr,"not enough memory for mn_cfc11!\n");
	return 1;
    }
    mn_cs134 = alloc3d(NZ,NXMEM,NYMEM);
    if(mn_cs134 == NULL) {
	fprintf(stderr,"not enough memory for mn_cfc11!\n");
	return 1;
    }
    mn_cs136 = alloc3d(NZ,NXMEM,NYMEM);
    if(mn_cs136 == NULL) {
	fprintf(stderr,"not enough memory for mn_cfc11!\n");
	return 1;
    }
    mn_cs137 = alloc3d(NZ,NXMEM,NYMEM);
    if(mn_cs137 == NULL) {
	fprintf(stderr,"not enough memory for mn_cfc11!\n");
	return 1;
    }
    mn_ba140 = alloc3d(NZ,NXMEM,NYMEM);
    if(mn_ba140 == NULL) {
	fprintf(stderr,"not enough memory for mn_cfc11!\n");
	return 1;
    }
    mn_la140 = alloc3d(NZ,NXMEM,NYMEM);
    if(mn_la140 == NULL) {
	fprintf(stderr,"not enough memory for mn_cfc11!\n");
	return 1;
    }
#endif
#ifdef TTD_BP
    mn_ttd = alloc3d(NZ,NXMEM,NYMEM);
    if(mn_ttd == NULL) {
	fprintf(stderr,"not enough memory for mn_ttd!\n");
	return 1;
    }
#endif
// end ashao
#ifdef PHOSPHATE
    phosphate_init = alloc3d(NZ,NXMEM,NYMEM);
    if(phosphate_init == NULL) {
	fprintf(stderr,"not enough memory for phosphate_init!\n");
	return 1;
    }
    mn_phos = alloc3d(NZ,NXMEM,NYMEM);
    if(mn_phos == NULL) {
	fprintf(stderr,"not enough memory for mn_phos!\n");
	return 1;
    }
    mn_dop = alloc3d(NZ,NXMEM,NYMEM);
    if(mn_dop == NULL) {
	fprintf(stderr,"not enough memory for mn_dop!\n");
	return 1;
    }
    mn_pobs = alloc3d(NZ,NXMEM,NYMEM);
    if(mn_pobs == NULL) {
	fprintf(stderr,"not enough memory for mn_pobs!\n");
	return 1;
    }
    mn_jpo4 = alloc3d(NZ,NXMEM,NYMEM);
    if(mn_jpo4 == NULL) {
	fprintf(stderr,"not enough memory for mn_jpo4!\n");
	return 1;
    }
# ifdef PROGNOSTIC
    lightlim = alloc3d(NZ,NXMEM,NYMEM);
    if(lightlim == NULL) {
	fprintf(stderr,"not enough memory for lightlim!\n");
	return 1;
    }
    ironlim = alloc3d(NZ,NXMEM,NYMEM);
    if(ironlim == NULL) {
	fprintf(stderr,"not enough memory for ironlim!\n");
	return 1;
    }
    mn_lightlim = alloc3d(NZ,NXMEM,NYMEM);
    if(mn_lightlim == NULL) {
	fprintf(stderr,"not enough memory for mn_lightlim!\n");
	return 1;
    }
    mn_ironlim = alloc3d(NZ,NXMEM,NYMEM);
    if(mn_ironlim == NULL) {
	fprintf(stderr,"not enough memory for mn_ironlim!\n");
	return 1;
    }
#  ifdef NITRATE
    nitrlim = alloc3d(NZ,NXMEM,NYMEM);
    if(nitrlim == NULL) {
	fprintf(stderr,"not enough memory for nitrlim!\n");
	return 1;
    }
    mn_nitrlim = alloc3d(NZ,NXMEM,NYMEM);
    if(mn_nitrlim == NULL) {
	fprintf(stderr,"not enough memory for mn_nitrlim!\n");
	return 1;
    }
#  endif  /* NITRATE */
# endif /* PROGNOSTIC */
#endif /* PHOSPHATE */
#ifdef NITRATE
    nitrate_init = alloc3d(NZ,NXMEM,NYMEM);
    if(nitrate_init == NULL) {
	fprintf(stderr,"not enough memory for nitrate_init!\n");
	return 1;
    }
    mn_nitr = alloc3d(NZ,NXMEM,NYMEM);
    if(mn_nitr == NULL) {
	fprintf(stderr,"not enough memory for mn_nitr!\n");
	return 1;
    }
    mn_don = alloc3d(NZ,NXMEM,NYMEM);
    if(mn_don == NULL) {
	fprintf(stderr,"not enough memory for mn_don!\n");
	return 1;
    }
    mn_nobs = alloc3d(NZ,NXMEM,NYMEM);
    if(mn_nobs == NULL) {
	fprintf(stderr,"not enough memory for mn_nobs!\n");
	return 1;
    }
    mn_jno3 = alloc3d(NZ,NXMEM,NYMEM);
    if(mn_jno3 == NULL) {
	fprintf(stderr,"not enough memory for mn_jno3!\n");
	return 1;
    }
    mn_n2fix = alloc3d(NZ,NXMEM,NYMEM);
    if(mn_n2fix == NULL) {
	fprintf(stderr,"not enough memory for mn_n2fix!\n");
	return 1;
    }
    mn_denitw = alloc3d(NZ,NXMEM,NYMEM);
    if(mn_denitw == NULL) {
	fprintf(stderr,"not enough memory for mn_denitw!\n");
	return 1;
    }
    mn_denits = alloc3d(NZ,NXMEM,NYMEM);
    if(mn_denits == NULL) {
	fprintf(stderr,"not enough memory for mn_denits!\n");
	return 1;
    }
# ifdef N15_CYCLE
    nitrate_init_15n = alloc3d(NZ,NXMEM,NYMEM);
    if(nitrate_init_15n == NULL) {
	fprintf(stderr,"not enough memory for nitrate_init_15n!\n");
	return 1;
    }
    mn_nitr_15n= alloc3d(NZ,NXMEM,NYMEM);
    if(mn_nitr_15n == NULL) {
	fprintf(stderr,"not enough memory for mn_nitr_15n!\n");
	return 1;
    }
    mn_don_15n = alloc3d(NZ,NXMEM,NYMEM);
    if(mn_don_15n == NULL) {
	fprintf(stderr,"not enough memory for mn_don_15n!\n");
	return 1;
    }
    mn_jno3_15n = alloc3d(NZ,NXMEM,NYMEM);
    if(mn_jno3_15n == NULL) {
	fprintf(stderr,"not enough memory for mn_jno3_15n!\n");
	return 1;
    }
# endif /* N15_CYCLE */
#endif /* NITRATE */
#ifdef DIC
    jdic = alloc3d(NZ,NXMEM,NYMEM);
    if(jdic == NULL) {
	fprintf(stderr,"not enough memory for jdic!\n");
	return 1;
    }
    mn_dic = alloc3d(NZ,NXMEM,NYMEM);
    if(mn_dic == NULL) {
	fprintf(stderr,"not enough memory for mn_dic!\n");
	return 1;
    }
    jalk = alloc3d(NZ,NXMEM,NYMEM);
    if(jalk == NULL) {
	fprintf(stderr,"not enough memory for jalk!\n");
	return 1;
    }
    mn_alk = alloc3d(NZ,NXMEM,NYMEM);
    if(mn_alk == NULL) {
	fprintf(stderr,"not enough memory for mn_alk!\n");
	return 1;
    }
    mn_jdic = alloc3d(NZ,NXMEM,NYMEM);
    if(mn_jdic == NULL) {
	fprintf(stderr,"not enough memory for mn_jdic!\n");
	return 1;
    }
#endif /* DIC */

  return 0;
}


int alloc_tracadv(int NZED, int NY, int NX, double ***vartrac)
{

    vartrac = alloc3d(NZED,NY,NY);
    if(vartrac == NULL) {
	fprintf(stderr,"not enough memory for tracer variables!\n");
	return 1;
    }

  return 0;
}

int alloc_error(char* arr_name)
{
  fprintf(stderr,"not enough memory for %s!\n", arr_name);
  quit(2);
  return 2;
}
