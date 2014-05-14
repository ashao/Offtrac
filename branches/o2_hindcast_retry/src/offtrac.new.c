/********+*********+*********+*********+*********+*********+*********+*
 *                                                                    *
 *                OFFTRAC - off-line tracer advection code            *
 *                                                                    *
 *                    David Darr - April/May 2002                     *
 *                                                                    *
 *                    OCMIP biogeochemistry added                     *
 *                     cdeutsch 2004                                  *
 ********+*********+*********+*********+*********+*********+*********+*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <errno.h>
#include <netcdf.h>

#include "init.h"
#define MAINFILE
#include "metrics.h"
#include "io.h"
#include "util.h"
#include "alloc_arrays.h"
#include "alloc_trac.h"
#include "initialize.h"
#include "iocdf.h"
#include "read.h"

/*-------------------------------------------------------------------*
 *                                                                   *
 *     define variables and subroutines
 *                                                                   *
 *-------------------------------------------------------------------*/

/*  To add a variable, increase NOVARS (in init.h)  */
/*  then add information to struct vardesc (below) */
//BX  31AUG07  added missing value - last fields
//BX  31OCT07  Changed some variables to double precision (needed for restart)
//BX           All variables were float ('f') before
//BX           This is a bug fix as memory restrictions do not allow to create
//BX           a separate vardesc field for restarts.

//BX-a
//BX  missing values are set here and to the same value in "vars" declaration
//BX  they will be written on their respective mn_ fields with mult_.._mv calls
//HF the parameter and the macro need to be the same!
const double misval = -1.e+6;
#define MISVAL -1.e+6

const double r_bio_tau = 1.0 / (30.0 * 86400.0); // restoring time step
#ifdef PROGNOSTIC
int use_default_progn;
double lambda0;
double k_I;
double par_b;
double k_PO4;
double k_Fe;
#endif 

# ifdef N15_CYCLE
const double parm_n15_std_fraction = 1.0;
const double parm_alpha_n15_prod = 0.995;
const double parm_alpha_n15_pon_remin = 1.0;
const double parm_alpha_n15_don_remin = 1.0;
const double parm_alpha_n15_denitw = 0.985;
const double parm_alpha_n15_denits = 1.0;
const double parm_alpha_n15_n2fix = 1.0;
# endif

//BX-e
struct vardesc vars[NOVARS] =
	{
		{
			"D", "Basin Depth", 'h', '1', '1', "meter", 'f', 0
		}, // 0
		{
			"mn_u",
			"Zonal Velocity",
			'u',
			'L',
			's',
			"meter second-1",
			'f',
			0
		}, // 1
		{
			"mn_v",
			"Meridional Velocity",
			'v',
			'L',
			's',
			"meter second-1",
			'f',
			0
		}, // 2
		{
			"mn_h", "Layer Thickness", 'h', 'L', 's', "meter", 'd', 9.e-10
		}, // 3
		{
			"mn_uhtm",
			"zonal vel - cumul",
			'u',
			'L',
			's',
			"meter second-1",
			'f',
			0
		}, // 4
		{
			"mn_vhtm",
			"merid vel - cumul",
			'v',
			'L',
			's',
			"meter second-1",
			'f',
			0
		}, // 5
		{
			"mn_ea",
			"downward entrainment",
			'h',
			'L',
			's',
			"meters",
			'd',
			0
		}, // 6
		{
			"mn_eb", "upward entrainment", 'h', 'L', 's', "meters", 'd', 0
		}, // 7
		{
			"mn_eaml",
			"downward ML detrain",
			'h',
			'1',
			's',
			"meters",
			'd',
			0
		}, // 8
		{
			"mn_age",
			"age tracer",
			'h',
			'L',
			's',
			"years",
			'd',
			-3.17097930758233E-14
		},// 9
		{
			"mn_oxygen",
			"dissolved oxygen tracer",
			'h',
			'L',
			's',
			"mole m-3",
			'd',
			MISVAL
		},// 10
		{
			"mn_o2sat",
			"oxygen saturation",
			'h',
			'L',
			's',
			"mole m-3",
			'f',
			MISVAL
		}, // 11
		{
			"mn_jo2",
			"oxygen source-sink",
			'h',
			'L',
			's',
			"mole m-3 s-1",
			'f',
			5.465535149299944E-15
		}, // 12
		{
			"mn_oxyflux",
			"oxygen gasex flux into ML",
			'h',
			'1',
			's',
			"(mole m-3)*(m s-1)",
			'f',
			0
		}, // 13
		{
			"mn_cfc11",
			"dissolved cfc11 tracer",
			'h',
			'L',
			's',
			"mole liter-1",
			'd',
			MISVAL
		}, // 14
		{
			"mn_cfc12",
			"dissolved cfc12 tracer",
			'h',
			'L',
			's',
			"mole liter-1",
			'd',
			MISVAL
		}, // 15
		{
			"mn_c11flux",
			"cfc fluxed into ML",
			'h',
			'1',
			's',
			"(mole l-1)*(m s-1)",
			'f',
			MISVAL
		}, // 16
		{
			"mn_c12flux",
			"cfc fluxed into ML",
			'h',
			'1',
			's',
			"(mole l-1)*(m s-1)",
			'f',
			MISVAL
		}, // 17
		{
			"mn_rml",
			"ML potential density",
			'h',
			'2',
			's',
			"kg m-3",
			'f',
			MISVAL
		}, // 18
		{
			"mn_phos",
			"phosphate tracer",
			'h',
			'L',
			's',
			"mol/m3",
			'd',
			MISVAL
		}, // 19
		{
			"mn_nitr",
			"nitrate tracer",
			'h',
			'L',
			's',
			"mol/m3",
			'd',
			MISVAL
		}, // 20
		{
			"mn_dic",
			"dissolved inorganic carbon tracer",
			'h',
			'L',
			's',
			"mole m-3",
			'd',
			MISVAL
		}, // 21
		{
			"mn_dop", "dop tracer", 'h', 'L', 's', "mole m-3", 'd', MISVAL
		}, // 22
		{
			"mn_don", "don tracer", 'h', 'L', 's', "mole m-3", 'd', MISVAL
		}, // 23
		{
			"mn_pflux",
			"pop flux",
			'h',
			'1',
			's',
			"mmol/m2/sec",
			'f',
			MISVAL
		}, // 24
		{
			"mn_o18",
			"delta-18O of O2",
			'h',
			'L',
			's',
			"permil",
			'd',
			MISVAL
		}, // 25
		{
			"mn_jo18",
			"oxygen 18 source-sink",
			'h',
			'L',
			's',
			"mole m-3 s-1",
			'f',
			MISVAL
		},// 26
		{
			"mn_pobs",
			"observed P on sigma-layers",
			'h',
			'L',
			's',
			"mole m-3",
			'f',
			MISVAL
		},// 27
		{
			"mn_inert1",
			"inert tracer 1",
			'h',
			'L',
			's',
			"mole m-3",
			'd',
			MISVAL
		}, // 28
		{
			"mn_inert2",
			"inert tracer 2",
			'h',
			'L',
			's',
			"mole m-3",
			'd',
			MISVAL
		}, // 29
		//#if defined AGE2 || defined AGE3
		// HF WARN: this doesn't make any sense right now!!
		{
			"mn_hnew",
			"Layer Thickness",
			'h',
			'L',
			's',
			"meter",
			'f',
			9.e-10
		}, // 30
		//#endif
		{
			"mn_alk", "Alkalinity", 'h', 'L', 's', "mole m-3", 'd', MISVAL
		}, // 31
		{
			"mn_co2flux",
			"air-sea co2 flux",
			'h',
			'1',
			's',
			"mol/m2/sec",
			'f',
			0.
		}, // 32
		{
			"mn_jpo4",
			"phosphate source-sink",
			'h',
			'L',
			's',
			"mole m-3 s-1",
			'f',
			MISVAL
		},// 33
		{
			"mn_temp", "temperature", 'h', 'L', 's', "dgC", 'f', MISVAL
		}, // 34
		{
			"mn_lightlim",
			"light limitation",
			'h',
			'L',
			's',
			"-",
			'f',
			MISVAL
		}, // 35
		{
			"mn_nitrlim",
			"nitrate limitation",
			'h',
			'L',
			's',
			"-",
			'f',
			MISVAL
		}, // 36
		{
			"mn_ironlim",
			"iron limitation",
			'h',
			'L',
			's',
			"-",
			'f',
			MISVAL
		}, // 37
		{
			"mn_jno3",
			"nitrate source-sink",
			'h',
			'L',
			's',
			"mole m-3 s-1",
			'f',
			MISVAL
		},// 38
		{
			"mn_nitr_15n",
			"nitrate tracer 15N",
			'h',
			'L',
			's',
			"mol/m3",
			'd',
			MISVAL
		}, // 39
		{
			"mn_don_15n",
			"don tracer 15N",
			'h',
			'L',
			's',
			"mole m-3",
			'd',
			MISVAL
		}, // 40
		{
			"mn_jno3_15n",
			"nitrate 15N source-sink",
			'h',
			'L',
			's',
			"mole m-3 s-1",
			'f',
			MISVAL
		},// 41
		{
			"mn_nflux",
			"pon flux",
			'h',
			'1',
			's',
			"mmol/m2/sec",
			'f',
			MISVAL
		}, // 42
		{
			"mn_nobs",
			"observed N on sigma-layers",
			'h',
			'L',
			's',
			"mole m-3",
			'f',
			MISVAL
		},// 43
		{
			"mn_nflux_15n",
			"pon_15N flux",
			'h',
			'1',
			's',
			"mmol/m2/sec",
			'f',
			MISVAL
		},// 44
		{
			"mn_n2fix",
			"nitrate fixation",
			'h',
			'L',
			's',
			"mole m-3 s-1",
			'f',
			MISVAL
		},// 45
		{
			"mn_denitw",
			"denitrification",
			'h',
			'L',
			's',
			"mole m-3 s-1",
			'f',
			MISVAL
		}, // 46
		{
			"mn_denits",
			"sediment denitrification",
			'h',
			'L',
			's',
			"mole m-3 s-1",
			'f',
			MISVAL
		}, // 47
		{
			"mn_jdic",
			"dic source-sink",
			'h',
			'L',
			's',
			"mole m-3 s-1",
			'f',
			MISVAL
		}, // 48
		{
			"mn_pco2s",
			"surface pco2",
			'h',
			'1',
			's',
			"mole m-2 s-1",
			'f',
			MISVAL
		}, // 49
		// begin ashao
		{
			"mn_sf6",
			"dissolved sf6 tracer",
			'h',
			'L',
			's',
			"mole liter-1",
			'd',
			MISVAL
		}, // 50
		{
			"mn_sf6flux",
			"dissolved sf6 tracer flux",
			'h',
			'1',
			's',
			"mole m-2",
			'd',
			MISVAL
		}, // 51

		{
			"mn_i131",
			"Concentration of 131-I",
			'h',
			'L',
			's',
			"mole liter-1",
			'd',
			MISVAL
		}, //52
		{
			"mn_cs134",
			"Concentration of 134-Cs",
			'h',
			'L',
			's',
			"mole liter-1",
			'd',
			MISVAL
		}, //53
		{
			"mn_cs136",
			"Concentration of 136-Cs",
			'h',
			'L',
			's',
			"mole liter-1",
			'd',
			MISVAL
		}, //54
		{
			"mn_cs137",
			"Concentration of 137-Cs",
			'h',
			'L',
			's',
			"mole liter-1",
			'd',
			MISVAL
		}, //55
		{
			"mn_ba140",
			"Concentration of 131-Ba",
			'h',
			'L',
			's',
			"mole liter-1",
			'd',
			MISVAL
		}, //56
		{
			"mn_inert",
			"Concentration of a non-radioactive inert tracert",
			'h',
			'L',
			's',
			"mole liter-1",
			'd',
			MISVAL
		}, //57
		{
			"cfc11_sat",
			"Saturation Concentration of CFC-11",
			'h',
			'1',
			's',
			"mole liter-1",
			'd',
			MISVAL
		}, //58
		{
			"cfc12_sat",
			"Saturation Concentration of CFC-12",
			'h',
			'1',
			's',
			"mole liter-1",
			'd',
			MISVAL
		}, //59
		{
			"sf6_sat",
			"Saturation Concentration of SF6",
			'h',
			'1',
			's',
			"mole liter-1",
			'd',
			MISVAL
		}, //60
		{
			"TTD",
			"TTD Inferred from Boundary Impulse",
			'h',
			'L',
			's',
			"mole liter-1",
			'd',
			MISVAL
		} //61

		// end ashao
	};

void alloc_fields(void);
void set_metrics(void);

/* Begin edit DT */
void step_fields(int iyear, int itts, int imon, int iterno); // modified ashao
/*  End edit DT  */
void initialize_sponge(int from_file, char filename[]);

// mn_ refers to mean value for output interval, WRINT (see init.h)

double ***u;
double ***v;
double h[NZ][NXMEM][NYMEM];
double hend[NZ][NXMEM][NYMEM];
#ifdef HTEST
double htest[NZ][NXMEM][NYMEM];
#endif
double hstart[NZ][NXMEM][NYMEM];
#if defined AGE2 || defined AGE3
//  for debugging
//double hnew[NZ][NXMEM][NYMEM];
double ***hnew;
#endif
double depth[NZ][NXMEM][NYMEM];
double rml[2][NXMEM][NYMEM];
double ***mn_u;
double ***mn_v;
double ***mn_h;
double ***mn_rml;

double ***uhtm;
double ***vhtm;
double ***ea;
double ***eb;
double eaml[NXMEM][NYMEM];
double ***mn_uhtm;
double ***mn_vhtm;
double ***mn_ea;
double ***mn_eb;
double mn_eaml[NXMEM][NYMEM];
double wd[NZ + 1][NXMEM][NYMEM];

double Temptm[NZ][NXMEM][NYMEM];
double Salttm[NZ][NXMEM][NYMEM];
double salt_woa[NXMEM][NYMEM];
double windspeed[NXMEM][NYMEM];
double xkw[NXMEM][NYMEM];
double fice[NXMEM][NYMEM];
double atmpres[NXMEM][NYMEM];

double areagr[NXMEM][NYMEM];
double D[NXMEM][NYMEM];
double ****tr;
double trintegral[NTR];
double trprofile[NZ];
double trwarn[NTR][2];

double umask[NXMEM][NYMEM];
double vmask[NXMEM][NYMEM];

int ksub[NXMEM][NYMEM];
int ksub_clim[NXMEM][NYMEM];

int month;
int mlen[12];
int icnt;

double qlat[NYMEM];
double hlat[NYMEM];

#ifdef ENTRAIN
double ea_init[NZ][NXMEM][NYMEM];
double eaml_init[NXMEM][NYMEM];
#endif

/* Begin added DT */
int beginyear;
int theyear;
int lastsave;
/*  End added DT  */

#ifdef INERT
double inert1[NZ][NXMEM][NYMEM];
double inert2[NZ][NXMEM][NYMEM];
double inert1_init[NZ][NXMEM][NYMEM];
double inert2_init[NZ][NXMEM][NYMEM];
double ***mn_inert1;
double ***mn_inert2;
double ar_sat[NZ][NXMEM][NYMEM];
int mINERT1;
int mINERT2;
#endif

#ifdef AGE
//HF  double age[NZ][NXMEM][NYMEM];
double ***age;
double age_init[NZ][NXMEM][NYMEM];
double ***mn_age;
int mAGE;
# if defined (AGE2) || defined (AGE3)
double ***mn_hnew;
# endif
#endif

#ifdef SMOOTH_REST
//double ***h_init;
double h_init[NZ][NXMEM][NYMEM];
#endif

#ifdef OXYGEN
double ***oxy_init;
double o2_sat[NZ][NXMEM][NYMEM];
double oxyflux[NXMEM][NYMEM];
double jo2[NZ][NXMEM][NYMEM];
double jo2_spec[NZ][NXMEM][NYMEM];
double aou_spec[NZ][NXMEM][NYMEM];
int mOXYGEN;
double ***mn_oxygen;
double ***mn_o2sat;
double mn_oxyflux[NXMEM][NYMEM];
double ***mn_jo2;
# ifdef REST_ARCTIC
double oxy_lay[NZ][NXMEM][NYMEM];
# endif
#endif

#ifdef OXY18
double o18_init[NZ][NXMEM][NYMEM];
double o18[NZ][NXMEM][NYMEM];
double o18flux[NXMEM][NYMEM];
double jo18[NZ][NXMEM][NYMEM];
double delo18[NZ][NXMEM][NYMEM];
double Ro18o16[NZ][NXMEM][NYMEM];
int mo18;
double ***mn_o18;
double ***mn_jo18;
#endif

#ifdef CFC
double cfc11_init[NZ][NXMEM][NYMEM];
double cfc12_init[NZ][NXMEM][NYMEM];
double cfc11[NZ][NXMEM][NYMEM];
double cfc12[NZ][NXMEM][NYMEM];
double c11flux[NXMEM][NYMEM];
double c12flux[NXMEM][NYMEM];
double c11monflux[NXMEM][NYMEM];
double c12monflux[NXMEM][NYMEM];
double c11res[NZ][NXMEM][NYMEM];
double c12res[NZ][NXMEM][NYMEM];
int mCFC11, mCFC12;
double ***mn_cfc11;
double ***mn_cfc12;
double mn_c11flux[NXMEM][NYMEM];
double mn_c12flux[NXMEM][NYMEM];
/*    Begin added DT     */
double Sc_cfc11[NXMEM][NYMEM];
double Sc_cfc12[NXMEM][NYMEM];
double Sk_cfc11[NXMEM][NYMEM];
double Sk_cfc12[NXMEM][NYMEM];
double cfc11_sol[NZ][NXMEM][NYMEM];
double cfc12_sol[NZ][NXMEM][NYMEM];
double cfc11_sat[NXMEM][NYMEM];
double cfc12_sat[NXMEM][NYMEM];
double cfcatm11[NXMEM][NYMEM];
double cfcatm12[NXMEM][NYMEM];
// begin ashao

double cfc11sc_coeffs[4], cfc11sol_coeffs[7];
int cfc11year[MAXRECS];
double cfc11N[MAXRECS];
double cfc11S[MAXRECS];

double cfc12sc_coeffs[4], cfc12sol_coeffs[7];
int cfc12year[MAXRECS];
double cfc12N[MAXRECS];
double cfc12S[MAXRECS];
/*    End added DT       */

#ifdef SAT_TRACER
int mCFC11_sat,mCFC12_sat;
#endif

#endif
// begin ashao
#ifdef SF6
double sf6_init[NZ][NXMEM][NYMEM];
double sf6[NZ][NXMEM][NYMEM];
double sf6flux[NXMEM][NYMEM];
double sf6monflux[NXMEM][NYMEM];
double sf6res[NZ][NXMEM][NYMEM];
int mSF6;
double ***mn_sf6;
double mn_sf6flux[NXMEM][NYMEM];
double Sc_sf6[NXMEM][NYMEM];
double Sk_sf6[NXMEM][NYMEM];
double sf6_sol[NZ][NXMEM][NYMEM];
double sf6_sat[NXMEM][NYMEM];
double sf6atm[NXMEM][NYMEM];
double sf6sc_coeffs[4], sf6sol_coeffs[7];
int sf6year[MAXRECS];
double sf6N[MAXRECS];
double sf6S[MAXRECS];

#ifdef SAT_TRACER
int mSF6_sat;
#endif

#endif
double currtime;
int iterno;

#ifdef RADIOACTIVE
double i131[NZ][NXMEM][NYMEM];
double ***mn_i131;
double cs134[NZ][NXMEM][NYMEM];
double ***mn_cs134;
double cs136[NZ][NXMEM][NYMEM];
double ***mn_cs136;
double cs137[NZ][NXMEM][NYMEM];
double ***mn_cs137;
double ba140[NZ][NXMEM][NYMEM];
double ***mn_ba140;
double la140[NZ][NXMEM][NYMEM];
double ***mn_la140;

int mi131;
int mcs134;
int mcs136;
int mcs137;
int mba140;
int mla140;
#endif

#ifdef TTD_BP
double ttd[NZ][NXMEM][NYMEM];
double ttd_mask[NXMEM][NYMEM];
double ***mn_ttd;
int mttd;
#endif
// end ashao

/* Begin edit DT */
int estTTD;
/*  End edit DT  */

#ifdef PHOSPHATE
int mPHOSPHATE;
double ***phosphate_init;
double po4[NZ][NXMEM][NYMEM];
double jpo4[NZ][NXMEM][NYMEM];
double po4_star_lev[NZPHOS][NXMEM][NYMEM];
double po4_star_lay[NZ][NXMEM][NYMEM];
# ifdef PROGNOSTIC
double no3_lev[NZPHOS][NXMEM][NYMEM];
double no3_lay[NZ][NXMEM][NYMEM];
double fe_lev[NZPHOS][NXMEM][NYMEM];
double fe_lay[NZ][NXMEM][NYMEM];
double par_lay[NZ][NXMEM][NYMEM];
double sfc_swr[NXMEM][NYMEM];
# else
#  ifdef NITRATE
double no3_lev[NZPHOS][NXMEM][NYMEM];
double no3_lay[NZ][NXMEM][NYMEM];
#  endif
# endif
double jprod[NZ][NXMEM][NYMEM];
double jremin[NZ][NXMEM][NYMEM];
double jremdop[NZ][NXMEM][NYMEM];
int mDOP;
double jdop[NZ][NXMEM][NYMEM];
double flux_pop[NXMEM][NYMEM];
double ***mn_phos;
double ***mn_dop;
double ***mn_pobs;
double ***mn_jpo4;
double mn_pflux[NXMEM][NYMEM];
# ifdef PROGNOSTIC
double ***lightlim;
double ***mn_lightlim;
#  ifdef NITRATE
double ***nitrlim;
double ***mn_nitrlim;
#  endif
double ***ironlim;
double ***mn_ironlim;
# endif /* PROGNOSTIC */
#endif /* PHOSPHATE */

#ifdef NITRATE
int mNITRATE;
double ***nitrate_init;
double no3[NZ][NXMEM][NYMEM];
double jno3[NZ][NXMEM][NYMEM];
double jremdon[NZ][NXMEM][NYMEM];
double jprod_n[NZ][NXMEM][NYMEM];
double jn2fix[NZ][NXMEM][NYMEM];
double jremin_n[NZ][NXMEM][NYMEM];
double jdenitw[NZ][NXMEM][NYMEM];
double jdenits[NZ][NXMEM][NYMEM];
double flux_pon[NXMEM][NYMEM];

int mDON;
double jdon[NZ][NXMEM][NYMEM];

double ***mn_nitr;
double ***mn_don;
double ***mn_nobs;
double ***mn_jno3;
double ***mn_n2fix;
double ***mn_denitw;
double ***mn_denits;

double mn_nflux[NXMEM][NYMEM];
# ifdef N15_CYCLE
int mNITRATE_15n;
double ***nitrate_init_15n;
double no3_15n[NZ][NXMEM][NYMEM];
double jno3_15n[NZ][NXMEM][NYMEM];
double jremdon_15n[NZ][NXMEM][NYMEM];
double jprod_15n[NZ][NXMEM][NYMEM];
double jremin_15n[NZ][NXMEM][NYMEM];
double flux_pon_15n[NXMEM][NYMEM];

int mDON_15n;
double jdon_15n[NZ][NXMEM][NYMEM];

double ***mn_nitr_15n;
double ***mn_don_15n;

double ***mn_jno3_15n;

double mn_nflux_15n[NXMEM][NYMEM];
# endif /* N15_CYCLE */
#endif /* NITRATE */

#ifdef DIC
int mDIC;
double co2flux[NXMEM][NYMEM];
double co2_sat[NZ][NXMEM][NYMEM];
double dic[NZ][NXMEM][NYMEM];
double ***jdic;
double ***mn_dic;
double Sk_co2[NXMEM][NYMEM];

double deltaco2[NXMEM][NYMEM];
double dco2star;
double pco2surf;
double pco2s[NXMEM][NYMEM];
double mn_pco2s[NXMEM][NYMEM];

double mn_co2flux[NXMEM][NYMEM];
int mALK;
double alk[NZ][NXMEM][NYMEM];
double ***jalk;
double ***mn_alk;
double ***mn_jdic;
# ifdef REST_ARCTIC
double dic_lay[NZ][NXMEM][NYMEM];
double alk_lay[NZ][NXMEM][NYMEM];
# endif
#endif
double k0, k1, k2, kb, kw, ks, kf, k1p, k2p, k3p, ksi;
double ff, htotal, htotal2, xco2, xacc, x1, x2;
double bt, st, ft, sit, pt, dic1, ta, co2starair;

#ifdef SPONGE
double Idamp[NXMEM][NYMEM];
double num_fields;
double sp_e[NZ][NXMEM][NYMEM];
# if defined(AGE) || defined(CFC) || defined(SF6) // ashao
double sp_age[NZ][NXMEM][NYMEM];
# endif
# ifdef OXYGEN
double sp_oxygen[NZ][NXMEM][NYMEM];
double sp_o18[NZ][NXMEM][NYMEM];
# endif
# ifdef INERT
double sp_inert1[NZ][NXMEM][NYMEM];
double sp_inert2[NZ][NXMEM][NYMEM];
# endif
# ifdef PHOSPHATE
double sp_phosphate[NZ][NXMEM][NYMEM];
double sp_dop[NZ][NXMEM][NYMEM];
# endif
# ifdef NITRATE
double sp_nitrate[NZ][NXMEM][NYMEM];
double sp_don[NZ][NXMEM][NYMEM];
#  ifdef N15_CYCLE
double sp_nitrate_15n[NZ][NXMEM][NYMEM];
double sp_don_15n[NZ][NXMEM][NYMEM];
#  endif
# endif
# ifdef DIC
double sp_dic[NZ][NXMEM][NYMEM];
double sp_alk[NZ][NXMEM][NYMEM];
# endif
#endif

double dt;

double junk[(NZ + 1) * (NXMEM) * (NYMEM)];

double *var[NOVARS];
long varsize[NOVARS];
int flags[NOVARS];
int rflags[NOVARS];
char directory[75];
char fname[75];
char inittrac[200];
char kw_varname[100]; //ashao: Specify the name of the kw variable in the gasx input file
#ifdef WRTTS
struct varcdfinfo varinfo[NOVARS];
int varmap[NOVARS];
int itts; /* tracer time step counter */
FILE *fn;
char output_filename[200];
double wrts;
int nvar = 0, cdfid, timeid[2];
size_t nrec = 0;
#endif

int cfc11nrecs, cfc12nrecs, sf6nrecs;

/*-------------------------------------------------------------------*
 *
 *     begin main executable
 *
 *-------------------------------------------------------------------*/

int main(void)
    {
    double inmon, tmon;
    double mon, nxt, lst;
    double *iyr;
    double *nyr;
    double dyr, day;
    double ndyr;
    double *dy;
    double dmon[12];
    double dbmon;
    int imon, inxt, ilst;
# ifndef WRTTS
    int itts; /* tracer time step counter */
# endif
    int nmnfirst;
#ifdef SEPFILES
    double smon, snxt;
    int ismon, isnxt;
#endif
    /* Begin edit DT */
# ifdef TTDEST
    estTTD=1;
# endif
# ifndef TTDEST
    estTTD = 0;
# endif
    /*  End edit DT  */

#ifndef WRTTS
    size_t nrec = 0;
#endif
    int err, i, j, k;
    int cmon;
    int nmn;
    double frac;
    static int m;
#ifndef WRTTS
    int varmap[NOVARS];

    FILE *fn;
    char output_filename[200];
#endif
    char run_name[200];
    char restart_filename[200];
    struct vardesc var_out[NOVARS];
#ifndef WRTTS
    struct varcdfinfo varinfo[NOVARS];
    int nvar = 0, cdfid, timeid[2];
#endif

    extern int flags[NOVARS];
    extern int rflags[NOVARS];

    //BX-a  for testing only
    int status;
    char message[144];
    //BX-e

    //BX  allocate tracer fields
    err = alloc_arrays();
    if (err)
	printf("Error allocating arrays.\n");

    err = alloc_trac();
    if (err)
	printf("Error allocating tracer field.\n");

    iyr = malloc(sizeof(double));
    nyr = malloc(sizeof(double));
    dy = malloc(sizeof(double));

    mlen[0] = 31; /* January      */
    mlen[1] = 28; /* February     */
    mlen[2] = 31; /* March        */
    mlen[3] = 30; /* April        */
    mlen[4] = 31; /* May          */
    mlen[5] = 30; /* June         */
    mlen[6] = 31; /* July         */
    mlen[7] = 31; /* August       */
    mlen[8] = 30; /* September    */
    mlen[9] = 31; /* October      */
    mlen[10] = 30; /* November     */
    mlen[11] = 31; /* December     */

    dmon[0] = 0.0;

    for (i = 1; i <= 11; i++)
	{
	dmon[i] = dmon[i - 1] + mlen[i - 1];
	}

    /*----------------------------------*
     *
     *     get user input
     *
     *----------------------------------*/

    {
    char line[100];
    int scan_count, done = 1;

    printf("Enter directory to use to read.\n");
    fgets(directory, sizeof(directory), stdin);

    directory[strlen(directory) - 1] = '\0';
    k = strlen(directory);
    if (k > 0)
	if (directory[k - 1] != '/')
	    {
	    directory[k] = '/';
	    directory[k + 1] = '\0';
	    printf("Using directory %s first.\n", directory);
	    }

    strcat(directory, fname);
    strcpy(fname, directory);
    printf("file path = %s\n", fname);

    while (done)
	{

	printf(
		"\nEnter the starting month and the total months to integrate.\n");

	fgets(line, sizeof(line), stdin);
	line[strlen(line) - 1] = '\0';
	scan_count = sscanf(line, "%lg, %lg,", &inmon, &tmon);
	if (scan_count == 2)
	    {
	    if (inmon < 0 || tmon < 0)
		printf("Negative values not allowed\n");
	    else
		done = 0;
	    }
	else
	    printf("Incorrect number of values, %d, read.\n", scan_count);
	}

    printf("\ninitial month = %g \n", inmon);
    printf("final month = %g \n", inmon + tmon - 1);
    printf("total months = %g \n\n", tmon);

    /*-----------------------------------
     *
     *     set output print flags to 0
     *
     *     added restart flags as a restart bug fix until
     *     memory restriction problem is solved 31OCT07 BX
     *
     *----------------------------------*/

    for (i = 0; i <= NOVARS - 1; i++)
	flags[i] = 0;
    for (i = 0; i <= NOVARS - 1; i++)
	rflags[i] = 0;

    flags[1] = 0;
    flags[2] = 0; /* u,v */
    rflags[1] = 0;
    rflags[2] = 0; /* u,v */
    flags[0] = 0;
    flags[3] = 0; /* D, h */
    rflags[0] = 0;
    rflags[3] = 0; /* D, h */
    flags[4] = 0;
    flags[5] = 0; /* uhtm, vhtm */
    rflags[4] = 0;
    rflags[5] = 0; /* uhtm, vhtm */
    flags[6] = 0;
    flags[7] = 0;
    flags[8] = 0; /* ea, eb, eaml */
#ifdef ENTRAIN
    rflags[6]=1; rflags[7]=1; rflags[8]=1; /* ea, eb, eaml */
#endif
    flags[18] = 0; /* ML potential density */
    rflags[18] = 0; /* ML potential density */
#ifdef INERT
    flags[28]=1; rflags[28]=1; /* inert1 */
    flags[29]=1; rflags[29]=1; /* inert2 */
#endif
#ifdef AGE
    flags[9] = 0;
    rflags[9] = 1; /* ideal age tracer*/
# if defined (AGE2) || defined (AGE3)
    flags[30]=1; rflags[30]=0; /* alk tracer used for hnew field*/
# endif
#endif

#ifdef OXYGEN
    flags[10]=1; rflags[10]=1; /* oxygen tracer */
    flags[11]=1; rflags[11]=0; /* oxygen saturation */
    flags[12]=1; rflags[12]=0; /* oxygen source-sink */
    flags[13]=1; rflags[13]=0; /* oxygen fluxed into ML */
#endif
#ifdef OXY18
    flags[25]=1; rflags[25]=1; /* oxygen 18 tracer */
    flags[26]=0; rflags[26]=0; /* oxygen 18 source-sink */
#endif

#ifdef CFC
    flags[14] = 1;
    rflags[14] = 1; /* cfc11 tracer */
    flags[15] = 1;
    rflags[15] = 1; /* cfc12 tracer */
    flags[16] = 1;
    rflags[16] = 1; /* cfc11 fluxed into ML */
    flags[17] = 1;
    rflags[17] = 1; /* cfc12 fluxed into ML */
    flags[58] = 1;
    rflags[58] = 1; /* cfc11 ML sat*/
    flags[59] = 1;
    rflags[59] = 1; /* cfc12 ML sat */

#ifdef SAT_TRACER
    flags[62]=1;
    rflags[62]=1; /* cfc11 percent saturation */
    flags[63]=1;
    rflags[63]=1; /* cfc11 percent saturation */
#endif
#endif

    // begin ashao
#ifdef SF6
    flags[50] = 1;
    rflags[50] = 1; /* sf6 tracer */
    flags[51] = 1;
    rflags[51] = 1; /* sf6 fluxed into ML */
    flags[60] = 1;
    rflags[60] = 1; /* sf6 ML sat */
#ifdef SAT_TRACER
    flags[64] = 1;
    rflags[64] = 1; /* sf6 percent saturation */
#endif

#endif
#ifdef RADIOACTIVE
    flags[52]=1; rflags[52]=1; /* 131-I tracer */
    flags[53]=1; rflags[53]=1; /* 134-Cs tracer */
    flags[54]=1; rflags[54]=1; /* 136-Cs tracer */
    flags[55]=1; rflags[55]=1; /* 137-Cs tracer */
    flags[56]=1; rflags[56]=1; /* 140-Ba tracer */
    flags[57]=1; rflags[57]=1; /* 140-La tracer */
#endif
#ifdef TTD_BP
    flags[61] = 1; rflags[61]=1; /* TTD Boundary Propagator */
#endif

    // end ashao

#ifdef PHOSPHATE
    flags[19]=1; rflags[19]=1; /* phosphate tracer */
    flags[22]=1; rflags[22]=1; /* dop tracer */
    flags[24]=1; rflags[24]=0; /* pop flux */
    flags[27]=1; rflags[27]=0; /* po4 obs */
#endif

#ifdef NITRATE
    flags[20]=1; rflags[20]=1; /* nitrate tracer */
    flags[23]=1; rflags[23]=1; /* don tracer */
# ifdef N15_CYCLE
    flags[39]=1; rflags[39]=1; /* nitrate tracer 15N*/
    flags[40]=1; rflags[40]=1; /* don tracer 15N */
# endif
#endif /* NITRATE */
#ifdef DIC
    flags[21]=1; rflags[21]=1; /* dic tracer */
    flags[31]=1; rflags[31]=1; /* alk tracer */
    flags[32]=1; rflags[32]=0; /* air-sea co2 flux */
    flags[48]=1; rflags[48]=0; /* dic flux */
    flags[49]=1; rflags[49]=0; /* surface pco2 */
#endif

#ifdef PHOSPHATE
    flags[33]=1; rflags[33]=0; /* po4 flux */
# ifdef PROGNOSTIC
    flags[34]=1; rflags[34]=0; /* temperature */
    flags[35]=1; rflags[35]=0; /* light limitation of po4 prod. */
#  ifdef NITRATE
    flags[36]=1; rflags[36]=0; /* nitrate limitation of po4 prod. */
#  endif
    flags[37]=1; rflags[37]=0; /* iron limitation of po4 prod. */
# endif /* PROGNOSTIC */
#endif /* PHOSPHATE */
#ifdef NITRATE
    flags[38]=1; rflags[38]=0; /* no3 flux */
    flags[42]=1; rflags[42]=0; /* pon flux */
    flags[43]=1; rflags[43]=0; /* no3 obs */
    flags[45]=1; rflags[45]=0; /* n2 fixation */
    flags[46]=1; rflags[46]=0; /* denitrification */
    flags[47]=1; rflags[47]=0; /* sed. denitrification */
# ifdef N15_CYCLE
    flags[41]=1; rflags[41]=0; /* no3_15n flux */
    flags[44]=1; rflags[44]=0; /* pon_15n flux */
# endif
#endif /* NITRATE */

    printf("Enter base name for output\n");

    fgets(run_name, sizeof(run_name), stdin);
    run_name[strlen(run_name) - 1] = '\0';
    sprintf(output_filename, "%s.%04g.nc", run_name, inmon + tmon - 1);
    printf("Create NetCDF output file '%s'.\n", output_filename);

#ifdef GASEX
    printf("Enter kw variable name (e.g. xkw_wanninkhof)");
    gets(kw_varname);
#endif

#ifdef TTD_BP

    read_patch( ttd_mask );

#endif

# ifdef PROGNOSTIC
    //HF
    do
	{
	printf("Do you want to use the default values for the prognostic model?\n");
	printf("Enter \"0\" for 'NO' or a positive number for 'YES'.\n");
	fgets(line,sizeof(line),stdin);
	printf("read:%s<--\n", line);
	scan_count = sscanf(line, "%d", &use_default_progn);
	if (scan_count < 1)
	    printf("Please enter an integer number!\n");
	else
	    done = 1;
	}while (! done);

    if (! use_default_progn)
	{
	// get value for lambda0
	done = 0;
	do
	    {
	    printf("Please enter the value for lambda0 [mol/(m3 y)] (default=3e-3):\n");
	    fgets(line,sizeof(line),stdin);
	    scan_count = sscanf(line, "%lg", &lambda0);
	    if (scan_count < 1 || lambda0 <= 0.0)
		{
#ifdef INTERACTIVE
		printf("Please enter a positive number!\n");
#else
		lambda0 = 0;
		done = 1;
#endif
		}
	    else
		{
		lambda0 /= (365.0 * 86400.0); // convert to mol/(m3 s)
		done = 1;
		}
	    }while (! done);
	// get value for k_I
	done = 0;
	do
	    {
	    printf("Please enter the value for k_I [W/m2] (default=50):\n");
	    fgets(line,sizeof(line),stdin);
	    scan_count = sscanf(line, "%lg", &k_I);
	    if (scan_count < 1 || k_I <= 0.0)
		{
#ifdef INTERACTIVE
		printf("Please enter a positive number!\n");
#else
		k_I = 0;
		done = 1;
#endif
		}
	    else
		done = 1;
	    }while (! done);
	// get value for par_b
	done = 0;
	do
	    {
	    printf("Please enter the value for par_b [1/dgC] (default=0.04):\n");
	    fgets(line,sizeof(line),stdin);
	    scan_count = sscanf(line, "%lg", &par_b);
	    if (scan_count < 1 || par_b <= 0.0)
		{
#ifdef INTERACTIVE
		printf("Please enter a positive number!\n");
#else
		par_b = 0;
		done = 1;
#endif
		}
	    else
		done = 1;
	    }while (! done);
	// get value for k_PO4
	done = 0;
	do
	    {
	    printf("Please enter the value for k_PO4 [mol/m3] (default=0.2e-3):\n");
	    fgets(line,sizeof(line),stdin);
	    scan_count = sscanf(line, "%lg", &k_PO4);
	    if (scan_count < 1 || k_PO4 <= 0.0)
		{
#ifdef INTERACTIVE
		printf("Please enter a positive number!\n");
#else
		k_PO4 = 0;
		done = 1;
#endif
		}
	    else
		done = 1;
	    }while (! done);
	// get value for k_Fe
	done = 0;
	do
	    {
	    printf("Please enter the value for k_Fe [mol/m3] (default=0.2e-6):\n");
	    fgets(line,sizeof(line),stdin);
	    scan_count = sscanf(line, "%lg", &k_Fe);
	    if (scan_count < 1 || k_Fe <= 0.0)
		{
#ifdef INTERACTIVE
		printf("Please enter a positive number!\n");
#else
		k_Fe = 0;
		done = 1;
#endif
		}
	    else
		done = 1;
	    }while (! done);
	}
    else
	{
	// set default values
	lambda0 = 3e-3 / (365.0 * 86400.0); // mol/(m3 s)
	k_I = 50.0; // W/m2
	par_b = 0.04; // 1/dgC
	k_PO4 = 0.1e-3; // mol/m3
	k_Fe = 0.2e-6; // mol/m3
	}
    printf("----- Parameter for prognostic model:\n\n");
    printf("lambda0=%g\n",lambda0);
    printf("k_I=%g\n",k_I);
    printf("par_b=%g\n",par_b);
    printf("k_PO4=%g\n",k_PO4);
    printf("\n-------------------------------------\n\n");
    //HF-e
# endif

    } // end of block "get user input"

    //DT
    lastsave = -1;
    //DT

    /*-----------------------------------
     *
     *     allocate and initialize fields
     *
     *----------------------------------*/

    read_grid();
    printf("Done reading grid or metric file.\n");

    set_metrics();
    printf("Done setting metrics.\n");
    // begin ashao
# ifdef SF6

    sf6nrecs = read_trac(sf6sc_coeffs, sf6sol_coeffs, sf6year, sf6N, sf6S);
    printf("Read %d records for SF6\n\n", sf6nrecs);
# endif

# ifdef CFC
    cfc11nrecs = read_trac(cfc11sc_coeffs, cfc11sol_coeffs, cfc11year, cfc11N,
	    cfc11S);
    printf("Read %d records for CFC-11\n", cfc11nrecs);
    cfc12nrecs = read_trac(cfc12sc_coeffs, cfc12sol_coeffs, cfc12year, cfc12N,
	    cfc12S);
    printf("Read %d records for CFC-12\n\n", cfc12nrecs);

# endif
    /* Copy the variable descriptions to a list of the actual output variables. */
    for (i = 0; i < NOVARS; i++)
	if (flags[i] > 0)
	    {
	    var_out[nvar] = vars[i];
	    varmap[i] = nvar;
	    nvar++;
	    }
    // force float precision output with last argument
    printf("Making NETCDF %s file\n", output_filename);
    create_file(output_filename, NETCDF_FILE, var_out, nvar, &fn, &cdfid,
	    timeid, varinfo, 1);
    // don't force
    // create_file(output_filename, NETCDF_FILE, var_out, nvar, &fn, &cdfid, timeid, varinfo, 0);
    printf("Closing file \n");
    close_file(&cdfid, &fn);

    /* Allocate the memory for the fields to be calculated.		*/
    alloc_fields();

    /* initialize tracer pointers					*/

    for (m = 0; m < NOVARS; m++)
	{
	if (flags[m])
	    for (k = 0; k < varsize[m]; k++)
		var[m][k] = 0.0;
	}

    /********** set means to zero                                   */
    /*  2d variables */
    if (flags[8])
	set_fix_darray2d_zero(mn_eaml);
#ifdef CFC
    /*    Begin altered DT     */
    //  if (flags[16]) set_fix_darray2d_zero(mn_c11flux);
    //  if (flags[17]) set_fix_darray2d_zero(mn_c12flux);
    if (flags[16])
	{
	for (i = 0; i < NXMEM; i++)
	    for (j = 0; j < NYMEM; j++)
		mn_c11flux[i][j] = 0.0;
	}
    if (flags[17])
	{
	for (i = 0; i < NXMEM; i++)
	    for (j = 0; j < NYMEM; j++)
		mn_c12flux[i][j] = 0.0;
	}
    /*    End altered DT       */
#endif
    // begin ashao
#ifdef SF6
    if (flags[51])
	{
	for (i = 0; i < NXMEM; i++)
	    for (j = 0; j < NYMEM; j++)
		mn_sf6flux[i][j] = 0.0;
	}
#endif
    // end ashao


#ifdef PHOSPHATE
    if (flags[24]) set_fix_darray2d_zero(mn_pflux);
    if (flags[33]) set_darray3d_zero(mn_jpo4, NZ, NXMEM, NYMEM);
# ifdef PROGNOSTIC
    set_darray3d_zero(lightlim, NZ, NXMEM, NYMEM);
    set_darray3d_zero(ironlim, NZ, NXMEM, NYMEM);
    if (flags[35]) set_darray3d_zero(mn_lightlim, NZ, NXMEM, NYMEM);
    if (flags[37]) set_darray3d_zero(mn_ironlim, NZ, NXMEM, NYMEM);
#  ifdef NITRATE
    set_darray3d_zero(nitrlim, NZ, NXMEM, NYMEM);
    if (flags[36]) set_darray3d_zero(mn_nitrlim, NZ, NXMEM, NYMEM);
#  endif
# endif /* PROGNOSTIC */
#endif /* PHOSPHATE */
#ifdef NITRATE
    if (flags[42]) set_fix_darray2d_zero(mn_nflux);
    if (flags[38]) set_darray3d_zero(mn_jno3, NZ, NXMEM, NYMEM);
    if (flags[45]) set_darray3d_zero(mn_n2fix, NZ, NXMEM, NYMEM);
    if (flags[46]) set_darray3d_zero(mn_denitw, NZ, NXMEM, NYMEM);
    if (flags[47]) set_darray3d_zero(mn_denits, NZ, NXMEM, NYMEM);
# ifdef N15_CYCLE
    if (flags[44]) set_fix_darray2d_zero(mn_nflux_15n);
    if (flags[41]) set_darray3d_zero(mn_jno3_15n, NZ, NXMEM, NYMEM);
# endif
#endif
#ifdef OXYGEN
    if (flags[13]) set_fix_darray2d_zero(mn_oxyflux);
#endif
#ifdef DIC
    if (flags[32]) set_fix_darray2d_zero(mn_co2flux);
    if (flags[49]) set_fix_darray2d_zero(mn_pco2s);
#endif
    /*  3d variables */
    if (flags[1])
	set_darray3d_zero(mn_u, NZ, NXMEM, NYMEM);
    if (flags[2])
	set_darray3d_zero(mn_v, NZ, NXMEM, NYMEM);
    if (flags[3])
	set_darray3d_zero(mn_h, NZ, NXMEM, NYMEM);
    if (flags[4])
	set_darray3d_zero(mn_uhtm, NZ, NXMEM, NYMEM);
    if (flags[5])
	set_darray3d_zero(mn_vhtm, NZ, NXMEM, NYMEM);
    if (flags[6])
	set_darray3d_zero(mn_ea, NZ, NXMEM, NYMEM);
    if (flags[7])
	set_darray3d_zero(mn_eb, NZ, NXMEM, NYMEM);
#ifdef AGE
    if (flags[9])
	set_darray3d_zero(mn_age, NZ, NXMEM, NYMEM);
# if defined (AGE2) || defined (AGE3)
    if (flags[30]) set_darray3d_zero(mn_hnew, NZ, NXMEM, NYMEM);
# endif
#endif
#ifdef OXYGEN
    if (flags[10]) set_darray3d_zero(mn_oxygen, NZ, NXMEM, NYMEM);
    if (flags[11]) set_darray3d_zero(mn_o2sat, NZ, NXMEM, NYMEM);
    if (flags[12]) set_darray3d_zero(mn_jo2, NZ, NXMEM, NYMEM);
#endif
#ifdef CFC
    if (flags[14])
	set_darray3d_zero(mn_cfc11, NZ, NXMEM, NYMEM);
    if (flags[15])
	set_darray3d_zero(mn_cfc11, NZ, NXMEM, NYMEM);
#endif
    // begin ashao
# ifdef SF6
    if (flags[50])
	set_darray3d_zero(mn_sf6, NZ, NXMEM, NYMEM);
# endif
#ifdef RADIOACTIVE
    if (flags[52]) set_darray3d_zero(mn_i131, NZ, NXMEM, NYMEM);
    if (flags[53]) set_darray3d_zero(mn_cs134, NZ, NXMEM, NYMEM);
    if (flags[54]) set_darray3d_zero(mn_cs136, NZ, NXMEM, NYMEM);
    if (flags[55]) set_darray3d_zero(mn_cs137, NZ, NXMEM, NYMEM);
    if (flags[56]) set_darray3d_zero(mn_ba140, NZ, NXMEM, NYMEM);
    if (flags[57]) set_darray3d_zero(mn_la140, NZ, NXMEM, NYMEM);
#endif
#ifdef TTD_BP
    if (flags[61]) set_darray3d_zero(mn_ttd,NZ,NXMEM,NYMEM);
#endif
    // end ashao
#ifdef PHOSPHATE
    if (flags[19]) set_darray3d_zero(mn_phos, NZ, NXMEM, NYMEM);
    if (flags[22]) set_darray3d_zero(mn_dop, NZ, NXMEM, NYMEM);
    if (flags[27]) set_darray3d_zero(mn_pobs, NZ, NXMEM, NYMEM);
#endif
#ifdef NITRATE
    if (flags[20]) set_darray3d_zero(mn_nitr, NZ, NXMEM, NYMEM);
    if (flags[23]) set_darray3d_zero(mn_don, NZ, NXMEM, NYMEM);
    if (flags[43]) set_darray3d_zero(mn_nobs, NZ, NXMEM, NYMEM);
# ifdef N15_CYCLE
    if (flags[39]) set_darray3d_zero(mn_nitr_15n, NZ, NXMEM, NYMEM);
    if (flags[40]) set_darray3d_zero(mn_don_15n, NZ, NXMEM, NYMEM);
# endif
#endif /* NITRATE */
#ifdef DIC
    if (flags[21]) set_darray3d_zero(mn_dic, NZ, NXMEM, NYMEM);
    if (flags[31]) set_darray3d_zero(mn_alk, NZ, NXMEM, NYMEM);
    if (flags[48]) set_darray3d_zero(mn_jdic, NZ, NXMEM, NYMEM);
#endif
#ifdef OXY18
    if (flags[25]) set_darray3d_zero(mn_o18, NZ, NXMEM, NYMEM);
    if (flags[26]) set_darray3d_zero(mn_jo18, NZ, NXMEM, NYMEM);
#endif
#ifdef INERT
    if (flags[28]) set_darray3d_zero(mn_inert1, NZ, NXMEM, NYMEM);
    if (flags[29]) set_darray3d_zero(mn_inert2, NZ, NXMEM, NYMEM);
#endif

    printf("Reading bathymetry, D.\n");

    // initialize D to be all ocean first

    for (i = 0; i <= NXMEM - 1; i++)
	{
	for (j = 0; j <= NYMEM - 1; j++)
	    {
	    D[i][j] = MINIMUM_DEPTH;
	    }
	}
#ifndef USE_CALC_H
    printf("calling read_D\n");
    read_D();

    for (i = 0; i <= NXMEM - 1; i++)
	for (j = 0; j <= NYMEM - 1; j++)
	    if (D[i][j] < 10.0)
		D[i][j] = MINIMUM_DEPTH;
#endif

    read_grid();
    set_metrics();

    dyr = inmon / 12;
#ifdef SEPFILES
    smon = (double) ((int) inmon % NMONTHS);
#endif
    mon = 12 * modf(dyr, iyr);

#ifdef SEPFILES
    printf("\n initial mon = %g - %g - %g \n\n", inmon, smon, mon);
#else
    printf("\n initial mon = %g - %g \n\n",inmon,mon);
#endif

    imon = (int) (mon + 0.00001);
    /* Begin edit DT */
    int iyear = floor((inmon + imon) / 12);
    theyear = iyear;
#ifdef RESTART
    theyear++;
    iyear++;
#endif
    /*  End edit DT  */
    inxt = imon + 1;
#ifdef SEPFILES
    ismon = (int) (smon + 0.00001);
    if ((ismon + 1) == NMONTHS)
	{
	isnxt = 0;
	}
    else
	{
	isnxt = ismon + 1;
	}
#endif

    dbmon = (double) inmon;
    lst = 12 * modf((dbmon - 1 + 12) / 12, iyr);
    ilst = (int) (lst + 0.00001);

#ifndef VARIAB_FORC
#  ifdef SEPFILES
    //BX  files are not in regular order (uvw are midmonth and h starts with last month)
    read_uvw(isnxt, 1);
    read_h(ismon, isnxt);
    // for files in regular order (h before uvw) use code here
    //BX
    // read_uvw(imon,1);
    //BX
    // read_h(imon,inxt);
#  else
    read_clim(imon,inxt,ilst);
#  endif
#endif
#ifdef USE_CALC_H
    z_sum(h, D);
#endif
    printf("\nSetting up and initializing\n");

#ifdef RESTART
    initialize(inmon,run_name);
#else
    initialize(imon);
#endif
    nmn = 0;

    //HF the next line should be in init.h (and be optional!)
#undef OUTPUT_IC
#ifdef OUTPUT_IC
    /*-----------------------------------
     *
     *     write tracer initial conditions
     *
     *----------------------------------*/

    printf("Writing initial condition variables out to netCDF\n\n",cmon);

    //  open netcdf file for each writing
    status = nc_open(output_filename, NC_WRITE, &cdfid);
    if (status != NC_NOERR)
	{
	strcpy(message,"Opening file"); strcat(message,output_filename);
	handle_error(message,status);
	}

    err = write_time(cdfid,fn,timeid[0],nrec, dy);
    if (err == -1) printf("Error writing day.\n");

    if (flags[3])
	{
	for (k=0;k<NZ;k++)
	    for (i=0;i<NXMEM;i++)
		for (j=0;j<NYMEM;j++)
		    mn_h[k][i][j] += h[k][i][j];
	}

#ifdef AGE
    if (flags[9]) add_darray3d(mn_age, age, NZ, NXMEM, NYMEM);
    /*{
     for (k=0;k<NZ;k++)
     for (i=0;i<NXMEM;i++)
     for (j=0;j<NYMEM;j++)
     mn_age[k][i][j] += age[k][i][j];
     }*/
#endif /* AGE */
#ifdef OXYGEN
    if (flags[10])
	{
	for (k=0;k<NZ;k++)
	    for (i=0;i<NXMEM;i++)
		for (j=0;j<NYMEM;j++)
		    mn_oxygen[k][i][j] += tr[mOXYGEN][k][i][j];
	}
#endif /* OXYGEN */
#ifdef PHOSPHATE
    if (flags[19])
	{
	for (k=0;k<NZ;k++)
	    for (i=0;i<NXMEM;i++)
		for (j=0;j<NYMEM;j++)
		    mn_phos[k][i][j] += tr[mPHOSPHATE][k][i][j];
	}
    if (flags[22])
	{
	for (k=0;k<NZ;k++)
	    for (i=0;i<NXMEM;i++)
		for (j=0;j<NYMEM;j++)
		    mn_dop[k][i][j] += tr[mDOP][k][i][j];
	}
    if (flags[27])
	{
	for (k=0;k<NZ;k++)
	    for (i=0;i<NXMEM;i++)
		for (j=0;j<NYMEM;j++)
		    mn_pobs[k][i][j] += po4_star_lay[k][i][j];
	}
#endif /* PHOSPHATE */
#ifdef NITRATE
    if (flags[20])
	{
	for (k=0;k<NZ;k++)
	    for (i=0;i<NXMEM;i++)
		for (j=0;j<NYMEM;j++)
		    mn_nitr[k][i][j] += tr[mNITRATE][k][i][j];
	}
    if (flags[23])
	{
	for (k=0;k<NZ;k++)
	    for (i=0;i<NXMEM;i++)
		for (j=0;j<NYMEM;j++)
		    mn_don[k][i][j] += tr[mDON][k][i][j];
	}
    if (flags[43])
	{
	for (k=0;k<NZ;k++)
	    for (i=0;i<NXMEM;i++)
		for (j=0;j<NYMEM;j++)
		    mn_nobs[k][i][j] += no3_lay[k][i][j];
	}
# ifdef N15_CYCLE
    if (flags[39])
	{
	for (k=0;k<NZ;k++)
	    for (i=0;i<NXMEM;i++)
		for (j=0;j<NYMEM;j++)
		    mn_nitr_15n[k][i][j] += tr[mNITRATE_15n][k][i][j];
	}
    if (flags[40])
	{
	for (k=0;k<NZ;k++)
	    for (i=0;i<NXMEM;i++)
		for (j=0;j<NYMEM;j++)
		    mn_don_15n[k][i][j] += tr[mDON_15n][k][i][j];
	}
# endif /* N15_CYCLE */
#endif /* NITRATE */

#ifdef DIC
    if (flags[21])
	{
	for (k=0;k<NZ;k++)
	    for (i=0;i<NXMEM;i++)
		for (j=0;j<NYMEM;j++)
		    mn_dic[k][i][j] += dic[k][i][j];
	}
    if (flags[31])
	{
	for (k=0;k<NZ;k++)
	    for (i=0;i<NXMEM;i++)
		for (j=0;j<NYMEM;j++)
		    mn_alk[k][i][j] += alk[k][i][j];
	}
#endif /* DIC */
#ifdef OXY18
    if (flags[25])
	{
	for (k=0;k<NZ;k++)
	    for (i=0;i<NXMEM;i++)
		for (j=0;j<NYMEM;j++)
		    mn_o18[k][i][j] += delo18[k][i][j];
	}
#endif /* OXY18 */

    for (m=0;m<NOVARS;m++) if (flags[m])
	{
	err = write_field(cdfid, fn, vars[m], varinfo[varmap[m]],
		nrec,var[m]);
	if (err == -1) printf("Error writing %s.\n",vars[m].name);
	}
    //  close file after each writing
    close_file(&cdfid, &fn);
    printf("netcdf record = %d\n",nrec++);
#endif /* OUTPUT_IC */

    // begin ashao
    /*#ifdef RADIOACTIVE
     #define fukulonidx 62
     #define fukulatidx 145
     for (k=0;k<4;k++) {
     i131[k][fukulonidx][fukulatidx]=1.0e-6*1.3e2;
     cs134[k][fukulonidx][fukulatidx]=1.0e-6*3.1e1;
     cs136[k][fukulonidx][fukulatidx]=1.0e-6*2.8e0;
     cs137[k][fukulonidx][fukulatidx]=1.0e-6*3.2e1;
     ba140[k][fukulonidx][fukulatidx]=1.0e-6*5.0e0;
     la140[k][fukulonidx][fukulatidx]=1.0e-6*2.5e0;

     tr[mi131][k][fukulonidx][fukulatidx]=1.0e-6*1.3e2;
     tr[mcs134][k][fukulonidx][fukulatidx]=1.0e-6*3.1e1;
     tr[mcs136][k][fukulonidx][fukulatidx]=1.0e-6*2.8e0;
     tr[mcs137][k][fukulonidx][fukulatidx]=1.0e-6*3.2e1;
     tr[mba140][k][fukulonidx][fukulatidx]=1.0e-6*5.0e0;
     tr[mla140][k][fukulonidx][fukulatidx]=1.0e-6*2.5e0;
     }
     #endif
     */

    /*-------------------------------------------------------------------*
     *
     *     begin integration loop
     *
     *-------------------------------------------------------------------*/

    mon = 0.0; /* reiniti */
    nmnfirst = 1;

    currtime = BEGHIND - SPINUP; //ashao: keep track of the current time
    for (cmon = inmon; cmon < inmon + tmon; cmon++)
	{
	nmn++;

	iterno = cmon - inmon;
	printf("Iteration %d/%d\n", iterno, (int) (inmon + tmon - cmon));

	dyr = cmon / 12.0;
	ndyr = (cmon + 1) / 12.0;

	dbmon = (double) cmon;
	lst = 12 * modf((dbmon - 1 + 12) / 12, iyr);
	ilst = (int) (lst + 0.00001);
	printf("double-mon=%g,lastmon=%g,ilst=%i\n", dbmon, lst, ilst);

#ifdef SEPFILES
	smon = (double) ((int) cmon % NMONTHS);
	snxt = (double) ((int) (cmon + 1) % NMONTHS);
#endif   // files in regular order (h before uvw)
	mon = 12.0 * modf(dyr, iyr);
	nxt = 12.0 * modf(ndyr, nyr);
#ifdef SEPFILES
	printf("the current month is %i-%g-%g \n", cmon, smon, mon);
#else   // files in regular order (h before uvw)
	printf("the current month is %i-%g \n",cmon,mon);
#endif
	printf("the current year is %g \n", *iyr);

	imon = (int) (mon + 0.00001);
	inxt = (int) (nxt + 0.00001);
#ifdef SEPFILES
	ismon = (int) (smon + 0.00001);
	isnxt = (int) (snxt + 0.00001);
#endif
	day = (*iyr) * 365.0 + dmon[imon];
	*dy = day;
	printf("the current day is -%g- \n", *dy);
#ifdef SEPFILES
	printf("the current ismon/smon is -%i-%g- \n", ismon, smon);
	printf("the next smonth/mean index: -%i- \n", isnxt);
#else
	printf("the current imon/mon is -%i-%g- \n", imon,mon);
	printf("the next month/mean index: -%i-%i- \n",inxt,nmn);
#endif

	/*   Begin added DT     */
	if (mon == 0)
	    {
	    beginyear = 1;
	    }
	else
	    {
	    beginyear = 0;
	    }
	if (beginyear == 1 && lastsave != -1)
	    theyear++;
	lastsave++;
	iyear = theyear;
	/*   End added DT       */

	dt = ((double) mlen[imon]);
	dt = dt * 86400 / (double) NTSTEP;
# ifdef SEPFILES
	read_h(ismon, isnxt);
# else
	read_h(imon,inxt);
# endif
	for (itts = 1; itts <= NTSTEP; itts++)
	    {
	    printf("\nSub time step number= %i \n", itts);
	    /*-----------------------------------
	     *
	     *     get physical fields and bc's
	     *
	     *----------------------------------*/
#ifndef VARIAB_FORC
# ifdef SEPFILES
	    //BX  files are not in regular order (uvw are midmonth and h starts with last month)
	    read_uvw(isnxt, itts);
	    // for files in regular order (h before uvw) use code below
	    //BX read_uvw(imon,itts);
# else
	    read_clim(imon,inxt,ilst);
# endif
	    read_ts(isnxt, itts);
	    read_biotic_bc(imon, itts);
	    read_fields(inxt, itts);
#endif
# ifdef REST_ARCTIC
	    read_restoring(imon,itts);
# endif

#ifdef VARIAB_FORC
	read_hind( (int) cmon);
	    read_biotic_bc(imon, itts);
	    read_fields(imon, itts);
#endif

	    /*-----------------------------------
	     *
	     *     integrate 1 time step
	     *
	     *----------------------------------*/

	    printf("step fields - day = %g\n\n", day);

	    step_fields(iyear, itts, imon, iterno); // modified ashao
	    currtime += dt / 31536000; //ashao: advance currenttime

	    /*-------------------------------------------------------------------*
	     *
	     *     end integration loop
	     *
	     *-------------------------------------------------------------------*/

	    /*-----------------------------------
	     *
	     *     calculate tracer averages
	     *
	     *----------------------------------*/

	    printf("calculate tracer averages\n");

	    if (flags[8])
		add_fix_darray2d(mn_eaml, eaml);
#ifdef PHOSPHATE
	    if (flags[24]) add_fix_darray2d(mn_pflux, flux_pop);
#endif
#ifdef CFC
	    /*    Begin altered DT       */
	    //    if (flags[16]) add_fix_darray2d(mn_c11flux, c11flux);
	    //    if (flags[17]) add_fix_darray2d(mn_c12flux, c12flux);
	    if (flags[16])
		{
		for (i = 0; i < NXMEM; i++)
		    for (j = 0; j < NYMEM; j++)
			mn_c11flux[i][j] = c11monflux[i][j];
		}
	    if (flags[17])
		{
		for (i = 0; i < NXMEM; i++)
		    for (j = 0; j < NYMEM; j++)
			mn_c12flux[i][j] = c12monflux[i][j];
		}
	    /*       End altered DT       */
#endif
	    // begin ashao
#ifdef SF6
	    //    if (flags[16]) add_fix_darray2d(mn_c11flux, c11flux);
	    //    if (flags[17]) add_fix_darray2d(mn_c12flux, c12flux);
	    if (flags[51])
		{
		for (i = 0; i < NXMEM; i++)
		    for (j = 0; j < NYMEM; j++)
			mn_sf6flux[i][j] = sf6monflux[i][j];
		}
#endif
	    // end ashao
#ifdef DIC
	    if (flags[32]) add_fix_darray2d(mn_co2flux, co2flux);
	    if (flags[49]) add_fix_darray2d(mn_pco2s, pco2s);
#endif
#ifdef NITRATE
	    if (flags[42]) add_fix_darray2d(mn_nflux, flux_pon);
# ifdef N15_CYCLE
	    if (flags[44]) add_fix_darray2d(mn_nflux_15n, flux_pon_15n);
# endif
#endif
	    if (flags[1])
		add_darray3d(mn_u, u, NZ, NXMEM, NYMEM);
	    if (flags[2])
		add_darray3d(mn_v, v, NZ, NXMEM, NYMEM);
	    if (flags[3])
		{
		for (k = 0; k < NZ; k++)
		    for (i = 0; i < NXMEM; i++)
			for (j = 0; j < NYMEM; j++)
			    {


#ifdef HTEST
			    mn_h[k][i][j] += htest[k][i][j];
#else
			    mn_h[k][i][j] += h[k][i][j];
#endif
			    }
		}
	    if (flags[4])
		add_darray3d(mn_uhtm, uhtm, NZ, NXMEM, NYMEM);
	    if (flags[5])
		add_darray3d(mn_vhtm, vhtm, NZ, NXMEM, NYMEM);
	    if (flags[6])
		add_darray3d(mn_ea, ea, NZ, NXMEM, NYMEM);
	    if (flags[7])
		add_darray3d(mn_eb, eb, NZ, NXMEM, NYMEM);
#ifdef AGE
	    //DT
	    //    if (flags[9]) add_darray3d(mn_age, age, NZ, NXMEM, NYMEM);
	    if (flags[9])
		{
		for (k = 0; k < NZ; k++)
		    {
		    for (i = 0; i < NXMEM; i++)
			{
			for (j = 0; j < NYMEM; j++)
			    {
			    if (age[k][i][j] > 0.0)
				{
				if (estTTD == 0)
				    {
				    if (theyear == 0 && beginyear == 1)
					{
					mn_age[k][i][j] = age[k][i][j];
					}
				    else
					{
					mn_age[k][i][j] += age[k][i][j];
					}
				    }
				else
				    {
				    if (estTTD == 1)
					{
					if (theyear == BPTIME && beginyear == 1)
					    {
					    mn_age[k][i][j] = age[k][i][j];
					    }
					else
					    {
					    mn_age[k][i][j] += age[k][i][j];
					    }
					}
				    }
				}
			    }
			}
		    }
		}
	    //DT
# if defined (AGE2) || defined (AGE3)
	    if (flags[30]) add_darray3d(mn_hnew, hnew, NZ, NXMEM, NYMEM);
# endif
#endif
	    if (flags[18])
		{
		for (k = 0; k < 2; k++)
		    for (i = 0; i < NXMEM; i++)
			for (j = 0; j < NYMEM; j++)
			    mn_rml[k][i][j] += rml[k][i][j];
		}
#ifdef OXYGEN
	    if (flags[10])
		{
		for (k=0;k<NZ;k++)
		    for (i=0;i<NXMEM;i++)
			for (j=0;j<NYMEM;j++)
			    mn_oxygen[k][i][j] += tr[mOXYGEN][k][i][j];
		}
	    if (flags[11])
		{
		for (k=0;k<NZ;k++)
		    for (i=0;i<NXMEM;i++)
			for (j=0;j<NYMEM;j++)
			    mn_o2sat[k][i][j] += o2_sat[k][i][j];
		}
	    if (flags[12])
		{
		for (k=0;k<NZ;k++)
		    for (i=0;i<NXMEM;i++)
			for (j=0;j<NYMEM;j++)
			    mn_jo2[k][i][j] += jo2[k][i][j];
		}
	    if (flags[13])
		{
		for (i=0;i<NXMEM;i++)
		    for (j=0;j<NYMEM;j++)
			mn_oxyflux[i][j] += oxyflux[i][j];
		}
#endif /* OXYGEN */
#ifdef CFC
	    if (estTTD == 1)
		{
		if (flags[14])
		    {
		    for (k = 0; k < NZ; k++)
			{
			for (i = 0; i < NXMEM; i++)
			    {
			    for (j = 0; j < NYMEM; j++)
				{
				mn_cfc11[k][i][j] += cfc11[k][i][j];
				}
			    }
			}
		    }
		if (flags[15])
		    {
		    for (k = 0; k < NZ; k++)
			{
			for (i = 0; i < NXMEM; i++)
			    {
			    for (j = 0; j < NYMEM; j++)
				{
				if (theyear >= BPTIME)
				    {
				    if (theyear == BPTIME && beginyear == 1)
					{
					mn_cfc12[k][i][j] = cfc12[k][i][j];
					}
				    else
					{
					mn_cfc12[k][i][j] += cfc12[k][i][j];
					}
				    }
				else
				    {
				    mn_cfc12[k][i][j] += cfc12[k][i][j];
				    }
				}
			    }
			}
		    }
		}
	    else
		{
		if (flags[14])
		    {
		    for (k = 0; k < NZ; k++)
			{
			for (i = 0; i < NXMEM; i++)
			    {
			    for (j = 0; j < NYMEM; j++)
				{
				mn_cfc11[k][i][j] += cfc11[k][i][j];
				tr[mCFC11][k][i][j] = cfc11[k][i][j];
				}
			    }
			}
		    }
		if (flags[15])
		    {
		    for (k = 0; k < NZ; k++)
			{
			for (i = 0; i < NXMEM; i++)
			    {
			    for (j = 0; j < NYMEM; j++)
				{
				mn_cfc12[k][i][j] += cfc12[k][i][j];
				tr[mCFC12][k][i][j] = cfc12[k][i][j];
				}
			    }
			}
		    }
		}
#endif
	    // begin ashao
#ifdef SF6
	    if (estTTD == 1)
		{
		if (flags[50])
		    {
		    for (k = 0; k < NZ; k++)
			{
			for (i = 0; i < NXMEM; i++)
			    {
			    for (j = 0; j < NYMEM; j++)
				{
				if (theyear >= BPTIME)
				    {
				    if (theyear == BPTIME && beginyear == 1)
					{
					mn_sf6[k][i][j] = sf6[k][i][j];
					}
				    else
					{
					mn_sf6[k][i][j] += sf6[k][i][j];
					}
				    }
				else
				    {
				    mn_sf6[k][i][j] += sf6[k][i][j];
				    }
				}
			    }
			}
		    }
		}
	    else
		{
		if (flags[50])
		    {
		    for (k = 0; k < NZ; k++)
			{
			for (i = 0; i < NXMEM; i++)
			    {
			    for (j = 0; j < NYMEM; j++)
				{
				mn_sf6[k][i][j] += sf6[k][i][j];
				tr[mSF6][k][i][j] = sf6[k][i][j];
				}
			    }
			}
		    }
		}
#endif
#ifdef RADIOACTIVE 
	    for (k=0;k<NZ;k++)
		{
		for (i=0;i<NXMEM;i++)
		    {
		    for (j=0;j<NYMEM;j++)
			{
			mn_i131[k][i][j] += i131[k][i][j];
			tr[mi131][k][i][j] = i131[k][i][j];

			mn_cs134[k][i][j] += cs134[k][i][j];
			tr[mcs134][k][i][j] = cs134[k][i][j];

			mn_cs136[k][i][j] += cs136[k][i][j];
			tr[mcs136][k][i][j] = cs136[k][i][j];

			mn_cs137[k][i][j] += cs137[k][i][j];
			tr[mcs137][k][i][j] = cs137[k][i][j];

			mn_ba140[k][i][j] += ba140[k][i][j];
			tr[mba140][k][i][j] = ba140[k][i][j];

			mn_la140[k][i][j] += la140[k][i][j];
			tr[mla140][k][i][j] = la140[k][i][j];
			}
		    }
		}

#endif
#ifdef TTD_BP
	    for (k=0;k<NZ;k++)
		for (i=0;i<NXMEM;i++)
		    for (j=0;j<NYMEM;j++)
			{
			mn_ttd[k][i][j] = tr[mttd][k][i][j];
			// tr[mttd][k][i][j] = ttd[k][i][j];

			}
#endif

	    // end ashao
#ifdef PHOSPHATE
	    if (flags[19])
		{
		for (k=0;k<NZ;k++)
		    for (i=0;i<NXMEM;i++)
			for (j=0;j<NYMEM;j++)
			    mn_phos[k][i][j] += tr[mPHOSPHATE][k][i][j];
		}
	    if (flags[22])
		{
		for (k=0;k<NZ;k++)
		    for (i=0;i<NXMEM;i++)
			for (j=0;j<NYMEM;j++)
			    mn_dop[k][i][j] += tr[mDOP][k][i][j];
		}
	    if (flags[27])
		{
		for (k=0;k<NZ;k++)
		    for (i=0;i<NXMEM;i++)
			for (j=0;j<NYMEM;j++)
			    mn_pobs[k][i][j] += po4_star_lay[k][i][j];
		}
	    if (flags[33])
		{
		for (k=0;k<NZ;k++)
		    for (i=0;i<NXMEM;i++)
			for (j=0;j<NYMEM;j++)
			    mn_jpo4[k][i][j] += jpo4[k][i][j];
		}
# ifdef PROGNOSTIC
	    if (flags[35]) add_darray3d(mn_lightlim, lightlim, NZ, NXMEM, NYMEM);
#  ifdef NITRATE
	    if (flags[36])
		{
		for (k=0;k<NZ;k++)
		    for (i=0;i<NXMEM;i++)
			for (j=0;j<NYMEM;j++)
			    mn_nitrlim[k][i][j] += nitrlim[k][i][j];
		}
#  endif
	    if (flags[37]) add_darray3d(mn_ironlim, ironlim, NZ, NXMEM, NYMEM);
# endif /* PROGNOSTIC */
#endif /* PHOSPHATE */

#ifdef NITRATE
	    if (flags[20])
		{
		for (k=0;k<NZ;k++)
		    for (i=0;i<NXMEM;i++)
			for (j=0;j<NYMEM;j++)
			    mn_nitr[k][i][j] += tr[mNITRATE][k][i][j];
		}
	    if (flags[23])
		{
		for (k=0;k<NZ;k++)
		    for (i=0;i<NXMEM;i++)
			for (j=0;j<NYMEM;j++)
			    mn_don[k][i][j] += tr[mDON][k][i][j];
		}
	    if (flags[38])
		{
		for (k=0;k<NZ;k++)
		    for (i=0;i<NXMEM;i++)
			for (j=0;j<NYMEM;j++)
			    mn_jno3[k][i][j] += jno3[k][i][j];
		}
	    if (flags[43])
		{
		for (k=0;k<NZ;k++)
		    for (i=0;i<NXMEM;i++)
			for (j=0;j<NYMEM;j++)
			    mn_nobs[k][i][j] += no3_lay[k][i][j];
		}
	    if (flags[45])
		{
		for (k=0;k<NZ;k++)
		    for (i=0;i<NXMEM;i++)
			for (j=0;j<NYMEM;j++)
			    mn_n2fix[k][i][j] += jn2fix[k][i][j];
		}
	    if (flags[46])
		{
		for (k=0;k<NZ;k++)
		    for (i=0;i<NXMEM;i++)
			for (j=0;j<NYMEM;j++)
			    mn_denitw[k][i][j] += jdenitw[k][i][j];
		}
	    if (flags[47])
		{
		for (k=0;k<NZ;k++)
		    for (i=0;i<NXMEM;i++)
			for (j=0;j<NYMEM;j++)
			    mn_denits[k][i][j] += jdenits[k][i][j];
		}
# ifdef N15_CYCLE
	    if (flags[39])
		{
		for (k=0;k<NZ;k++)
		    for (i=0;i<NXMEM;i++)
			for (j=0;j<NYMEM;j++)
			    mn_nitr_15n[k][i][j] += tr[mNITRATE_15n][k][i][j];
		}
	    if (flags[40])
		{
		for (k=0;k<NZ;k++)
		    for (i=0;i<NXMEM;i++)
			for (j=0;j<NYMEM;j++)
			    mn_don_15n[k][i][j] += tr[mDON_15n][k][i][j];
		}
	    if (flags[41])
		{
		for (k=0;k<NZ;k++)
		    for (i=0;i<NXMEM;i++)
			for (j=0;j<NYMEM;j++)
			    mn_jno3_15n[k][i][j] += jno3_15n[k][i][j];
		}
# endif
#endif /* NITRATE */

#ifdef DIC
	    if (flags[21])
		{
		for (k=0;k<NZ;k++)
		    for (i=0;i<NXMEM;i++)
			for (j=0;j<NYMEM;j++)
			    mn_dic[k][i][j] += dic[k][i][j];
		}
	    if (flags[31])
		{
		for (k=0;k<NZ;k++)
		    for (i=0;i<NXMEM;i++)
			for (j=0;j<NYMEM;j++)
			    mn_alk[k][i][j] += alk[k][i][j];
		}
	    if (flags[48]) add_darray3d(mn_jdic, jdic, NZ, NXMEM, NYMEM);
#endif
#ifdef OXY18
	    if (flags[25])
		{
		for (k=0;k<NZ;k++)
		    for (i=0;i<NXMEM;i++)
			for (j=0;j<NYMEM;j++)
			    mn_o18[k][i][j] += delo18[k][i][j];
		}
	    if (flags[26])
		{
		for (k=0;k<NZ;k++)
		    for (i=0;i<NXMEM;i++)
			for (j=0;j<NYMEM;j++)
			    mn_jo18[k][i][j] += jo18[k][i][j];
		}
#endif
#ifdef INERT
	    if (flags[28])
		{
		for (k=0;k<NZ;k++)
		    for (i=0;i<NXMEM;i++)
			for (j=0;j<NYMEM;j++)
			    mn_inert1[k][i][j] += inert1[k][i][j];
		}
	    if (flags[29])
		{
		for (k=0;k<NZ;k++)
		    for (i=0;i<NXMEM;i++)
			for (j=0;j<NYMEM;j++)
			    mn_inert2[k][i][j] += inert2[k][i][j];
		}
#endif /* INERT */

	    /* calculate the mean */

	    if (nmn == WRINT && itts == 1)
		{
		printf("***nmn= %i, itts= %i, nmnfirst= %i\n", nmn, itts,
			nmnfirst);
		frac = 1.0 / ((double) nmn * (double) NTSTEP);

		if (flags[8])
		    mult_fix_darray2d(mn_eaml, frac);
#ifdef CFC
		if (flags[16])
		    mult_fix_darray2d_mv(mn_c11flux, frac, D, misval);
		if (flags[17])
		    mult_fix_darray2d_mv(mn_c12flux, frac, D, misval);
#endif
#ifdef SF6
		if (flags[51])
		    mult_fix_darray2d_mv(mn_sf6flux, frac, D, misval);
#endif
#ifdef PHOSPHATE
		if (flags[24]) mult_fix_darray2d_mv(mn_pflux, frac, D, misval);
#endif
#ifdef DIC
		if (flags[32]) mult_fix_darray2d(mn_co2flux, frac);
		if (flags[49]) mult_fix_darray2d(mn_pco2s, frac);
#endif
#ifdef NITRATE
		if (flags[42]) mult_fix_darray2d_mv(mn_nflux, frac, D, misval);
# ifdef N15_CYCLE
		if (flags[44]) mult_fix_darray2d_mv(mn_nflux_15n, frac, D, misval);
# endif
#endif
		if (flags[1])
		    mult_darray3d(mn_u, NZ, NXMEM, NYMEM, frac);
		if (flags[2])
		    mult_darray3d(mn_v, NZ, NXMEM, NYMEM, frac);
		if (flags[3])
		    mult_darray3d(mn_h, NZ, NXMEM, NYMEM, frac);
		if (flags[4])
		    mult_darray3d(mn_uhtm, NZ, NXMEM, NYMEM, frac);
		if (flags[5])
		    mult_darray3d(mn_vhtm, NZ, NXMEM, NYMEM, frac);
		if (flags[6])
		    mult_darray3d(mn_ea, NZ, NXMEM, NYMEM, frac);
		if (flags[7])
		    mult_darray3d(mn_eb, NZ, NXMEM, NYMEM, frac);
#ifdef AGE
		if (flags[9])
		    mult_darray3d(mn_age, NZ, NXMEM, NYMEM, frac);
# if defined (AGE2) || defined (AGE3)
		if (flags[30]) mult_darray3d(mn_hnew, NZ, NXMEM, NYMEM, frac);
# endif
#endif
#ifdef OXYGEN
		if (flags[10]) mult_darray3d_mv(mn_oxygen, NZ, NXMEM, NYMEM, frac, D, misval);
		if (flags[11]) mult_darray3d_mv(mn_o2sat, NZ, NXMEM, NYMEM, frac, D, misval);
		if (flags[12]) mult_darray3d(mn_jo2, NZ, NXMEM, NYMEM, frac);
		if (flags[13]) mult_fix_darray2d(mn_oxyflux, frac);
#endif
#ifdef CFC
		if (flags[14])
		    mult_darray3d_mv(mn_cfc11, NZ, NXMEM, NYMEM, frac, D,
			    misval);
		if (flags[15])
		    mult_darray3d_mv(mn_cfc12, NZ, NXMEM, NYMEM, frac, D,
			    misval);
		if (flags[16])
		    mult_fix_darray2d_mv(mn_c11flux, frac, D, misval);
		if (flags[17])
		    mult_fix_darray2d_mv(mn_c12flux, frac, D, misval);
#endif
		// begin ashao
#ifdef SF6
		if (flags[50])
		    mult_darray3d_mv(mn_sf6, NZ, NXMEM, NYMEM, frac, D, misval);
		if (flags[51])
		    mult_fix_darray2d_mv(mn_sf6flux, frac, D, misval);
#endif
#ifdef RADIOACTIVE
		if (flags[52]) mult_darray3d_mv(mn_i131, NZ, NXMEM, NYMEM, frac, D, misval);
		if (flags[53]) mult_darray3d_mv(mn_cs134, NZ, NXMEM, NYMEM, frac, D, misval);
		if (flags[54]) mult_darray3d_mv(mn_cs136, NZ, NXMEM, NYMEM, frac, D, misval);
		if (flags[55]) mult_darray3d_mv(mn_cs137, NZ, NXMEM, NYMEM, frac, D, misval);
		if (flags[56]) mult_darray3d_mv(mn_ba140, NZ, NXMEM, NYMEM, frac, D, misval);
		if (flags[57]) mult_darray3d_mv(mn_la140, NZ, NXMEM, NYMEM, frac, D, misval);
#endif
#ifdef TTD_BP
		if (flags[61]) mult_darray3d_mv(mn_ttd, NZ, NXMEM, NYMEM, frac, D, misval);
#endif
		// end ashao
		if (flags[18])
		    mult_darray3d_mv(mn_rml, 2, NXMEM, NYMEM, frac, D, misval);
#ifdef PHOSPHATE
		if (flags[19]) mult_darray3d_mv(mn_phos, NZ, NXMEM, NYMEM, frac, D, misval);
		if (flags[22]) mult_darray3d_mv(mn_dop, NZ, NXMEM, NYMEM, frac, D, misval);
		if (flags[27]) mult_darray3d_mv(mn_pobs, NZ, NXMEM, NYMEM, frac, D, misval);
		if (flags[33]) mult_darray3d_mv(mn_jpo4, NZ, NXMEM, NYMEM, frac, D, misval);
# ifdef PROGNOSTIC
		if (flags[35]) mult_darray3d_mv(mn_lightlim, NZ, NXMEM, NYMEM, frac, D, misval);
#  ifdef NITRATE
		if (flags[36]) mult_darray3d_mv(mn_nitrlim, NZ, NXMEM, NYMEM, frac, D, misval);
#  endif
		if (flags[37]) mult_darray3d_mv(mn_ironlim, NZ, NXMEM, NYMEM, frac, D, misval);
# endif /* PROGNOSTIC */
#endif /* PHOSPHATE */
#ifdef NITRATE
		if (flags[20]) mult_darray3d_mv(mn_nitr, NZ, NXMEM, NYMEM, frac, D, misval);
		if (flags[23]) mult_darray3d_mv(mn_don, NZ, NXMEM, NYMEM, frac, D, misval);
		if (flags[38]) mult_darray3d_mv(mn_jno3, NZ, NXMEM, NYMEM, frac, D, misval);
		if (flags[43]) mult_darray3d_mv(mn_nobs, NZ, NXMEM, NYMEM, frac, D, misval);
		if (flags[45]) mult_darray3d_mv(mn_n2fix, NZ, NXMEM, NYMEM, frac, D, misval);
		if (flags[46]) mult_darray3d_mv(mn_denitw, NZ, NXMEM, NYMEM, frac, D, misval);
		if (flags[47]) mult_darray3d_mv(mn_denits, NZ, NXMEM, NYMEM, frac, D, misval);
# ifdef N15_CYCLE
		if (flags[39]) mult_darray3d_mv(mn_nitr_15n, NZ, NXMEM, NYMEM, frac, D, misval);
		if (flags[40]) mult_darray3d_mv(mn_don_15n, NZ, NXMEM, NYMEM, frac, D, misval);
		if (flags[41]) mult_darray3d_mv(mn_jno3_15n, NZ, NXMEM, NYMEM, frac, D, misval);
# endif
#endif /* NITRATE */
#ifdef DIC
		if (flags[21]) mult_darray3d_mv(mn_dic, NZ, NXMEM, NYMEM, frac, D, misval);
		if (flags[31]) mult_darray3d_mv(mn_alk, NZ, NXMEM, NYMEM, frac, D, misval);
		if (flags[48]) mult_darray3d_mv(mn_jdic, NZ, NXMEM, NYMEM, frac, D, misval);
#endif
#ifdef OXY18
		if (flags[25]) mult_darray3d_mv(mn_o18, NZ, NXMEM, NYMEM, frac, D, misval);
		if (flags[26]) mult_darray3d_mv(mn_jo18, NZ, NXMEM, NYMEM, frac, D, misval);
#endif
#ifdef INERT
		if (flags[28]) mult_darray3d_mv(mn_inert1, NZ, NXMEM, NYMEM, frac, D, misval);
		if (flags[29]) mult_darray3d_mv(mn_inert2, NZ, NXMEM, NYMEM, frac, D, misval);
#endif

		/*-----------------------------------
		 *
		 *     write tracer averages and reset to 0
		 *
		 *----------------------------------*/
		printf("Writing month %i variables out to netCDF\n\n", cmon);

		status = nc_open(output_filename, NC_WRITE, &cdfid);
		if (status != NC_NOERR)
		    {
		    strcpy(message, "Opening file");
		    strcat(message, output_filename);
		    handle_error(message, status);
		    }

		err = write_time(cdfid, fn, timeid[0], nrec, dy);
		if (err == -1)
		    printf("Error writing day.\n");

		for (m = 0; m < NOVARS; m++)
		    if (flags[m]==1)
			{
			printf("m = %d\n",m);
			err = write_field(cdfid, fn, vars[m],
				varinfo[varmap[m]], nrec, var[m]);
			if (err == -1)
			    printf("Error writing %s.\n", vars[m].name);
			}
		close_file(&cdfid, &fn);

		/*    reset all means to zero */
		nmn = 0;
		if (flags[8])
		    set_fix_darray2d_zero(mn_eaml);
#ifdef CFC
		/*  Begin altered DT    */
		//      if (flags[16]) set_darray2d_zero(mn_c11flux, NXMEM, NYMEM);
		//      if (flags[17]) set_darray2d_zero(mn_c12flux, NXMEM, NYMEM);
		if (flags[16])
		    {
		    for (i = 0; i < NXMEM; i++)
			for (j = 0; j < NYMEM; j++)
			    mn_c11flux[i][j] = 0.0;
		    }
		if (flags[17])
		    {
		    for (i = 0; i < NXMEM; i++)
			for (j = 0; j < NYMEM; j++)
			    mn_c12flux[i][j] = 0.0;
		    }
		/*  End altered DT    */
#endif
		// begin ashao
# ifdef SF6
		//      if (flags[16]) set_darray2d_zero(mn_c11flux, NXMEM, NYMEM);
		//      if (flags[17]) set_darray2d_zero(mn_c12flux, NXMEM, NYMEM);
		if (flags[16])
		    {
		    for (i = 0; i < NXMEM; i++)
			for (j = 0; j < NYMEM; j++)
			    mn_sf6flux[i][j] = 0.0;
		    }
# endif
		// end ashao
#ifdef PHOSPHATE
		if (flags[24]) set_fix_darray2d_zero(mn_pflux);
#endif
#ifdef DIC
		if (flags[32]) set_fix_darray2d_zero(mn_co2flux);
		if (flags[49]) set_fix_darray2d_zero(mn_pco2s);
#endif
#ifdef NITRATE
		if (flags[42]) set_fix_darray2d_zero(mn_nflux);
# ifdef N15_CYCLE
		if (flags[44]) set_fix_darray2d_zero(mn_nflux_15n);
# endif
#endif
		if (flags[1])
		    set_darray3d_zero(mn_u, NZ, NXMEM, NYMEM);
		if (flags[2])
		    set_darray3d_zero(mn_v, NZ, NXMEM, NYMEM);
		if (flags[3])
		    set_darray3d_zero(mn_h, NZ, NXMEM, NYMEM);
		if (flags[4])
		    set_darray3d_zero(mn_uhtm, NZ, NXMEM, NYMEM);
		if (flags[5])
		    set_darray3d_zero(mn_vhtm, NZ, NXMEM, NYMEM);
		if (flags[6])
		    set_darray3d_zero(mn_ea, NZ, NXMEM, NYMEM);
		if (flags[7])
		    set_darray3d_zero(mn_eb, NZ, NXMEM, NYMEM);
		if (flags[18])
		    set_darray3d_zero(mn_rml, 2, NXMEM, NYMEM);
#ifdef INERT
		if (flags[28]) set_darray3d_zero(mn_inert1, NZ, NXMEM, NYMEM);
		if (flags[29]) set_darray3d_zero(mn_inert2, NZ, NXMEM, NYMEM);
#endif
#ifdef AGE
		if (flags[9])
		    set_darray3d_zero(mn_age, NZ, NXMEM, NYMEM);
# if defined (AGE2) || defined (AGE3)
		if (flags[30]) set_darray3d_zero(mn_hnew, NZ, NXMEM, NYMEM);
# endif
#endif
#ifdef OXYGEN
		if (flags[10]) set_darray3d_zero(mn_oxygen, NZ, NXMEM, NYMEM);
		if (flags[11]) set_darray3d_zero(mn_o2sat, NZ, NXMEM, NYMEM);
		if (flags[12]) set_darray3d_zero(mn_jo2, NZ, NXMEM, NYMEM);
		if (flags[13]) set_fix_darray2d_zero(mn_oxyflux);
#endif
#ifdef OXY18
		if (flags[25]) set_darray3d_zero(mn_o18, NZ, NXMEM, NYMEM);
		if (flags[26]) set_darray3d_zero(mn_jo18, NZ, NXMEM, NYMEM);
#endif
#ifdef CFC
		if (flags[14])
		    set_darray3d_zero(mn_cfc11, NZ, NXMEM, NYMEM);
		if (flags[15])
		    set_darray3d_zero(mn_cfc12, NZ, NXMEM, NYMEM);
#endif
		// begin ashao
#ifdef SF6
		if (flags[50])
		    set_darray3d_zero(mn_sf6, NZ, NXMEM, NYMEM);
#endif
#ifdef RADIOACTIVE
		if (flags[52]) set_darray3d_zero(mn_i131, NZ, NXMEM, NYMEM);
		if (flags[53]) set_darray3d_zero(mn_cs134, NZ, NXMEM, NYMEM);
		if (flags[54]) set_darray3d_zero(mn_cs136, NZ, NXMEM, NYMEM);
		if (flags[55]) set_darray3d_zero(mn_cs137, NZ, NXMEM, NYMEM);
		if (flags[56]) set_darray3d_zero(mn_ba140, NZ, NXMEM, NYMEM);
		if (flags[57]) set_darray3d_zero(mn_la140, NZ, NXMEM, NYMEM);
#endif
#ifdef TTD_BP
		if (flags[61]) set_darray3d_zero(mn_ttd, NZ, NXMEM, NYMEM);
#endif
		// end ashao
#ifdef PHOSPHATE
		if (flags[19]) set_darray3d_zero(mn_phos, NZ, NXMEM, NYMEM);
		if (flags[22]) set_darray3d_zero(mn_dop, NZ, NXMEM, NYMEM);
		if (flags[27]) set_darray3d_zero(mn_pobs, NZ, NXMEM, NYMEM);
		if (flags[33]) set_darray3d_zero(mn_jpo4, NZ, NXMEM, NYMEM);
# ifdef PROGNOSTIC
		set_darray3d_zero(lightlim, NZ, NXMEM, NYMEM);
#  ifdef NITRATE
		set_darray3d_zero(nitrlim, NZ, NXMEM, NYMEM);
		if (flags[36]) set_darray3d_zero(mn_nitrlim, NZ, NXMEM, NYMEM);
#  endif
		set_darray3d_zero(ironlim, NZ, NXMEM, NYMEM);
		if (flags[35]) set_darray3d_zero(mn_lightlim, NZ, NXMEM, NYMEM);
		if (flags[37]) set_darray3d_zero(mn_ironlim, NZ, NXMEM, NYMEM);
# endif /* PROGNOSTIC */
#endif /* PHOSPHATE */
#ifdef NITRATE
		if (flags[20]) set_darray3d_zero(mn_nitr, NZ, NXMEM, NYMEM);
		if (flags[23]) set_darray3d_zero(mn_don, NZ, NXMEM, NYMEM);
		if (flags[38]) set_darray3d_zero(mn_jno3, NZ, NXMEM, NYMEM);
		if (flags[43]) set_darray3d_zero(mn_nobs, NZ, NXMEM, NYMEM);
		if (flags[45]) set_darray3d_zero(mn_n2fix, NZ, NXMEM, NYMEM);
		if (flags[46]) set_darray3d_zero(mn_denitw, NZ, NXMEM, NYMEM);
		if (flags[47]) set_darray3d_zero(mn_denits, NZ, NXMEM, NYMEM);
# ifdef N15_CYCLE
		if (flags[39]) set_darray3d_zero(mn_nitr_15n, NZ, NXMEM, NYMEM);
		if (flags[40]) set_darray3d_zero(mn_don_15n, NZ, NXMEM, NYMEM);
		if (flags[41]) set_darray3d_zero(mn_jno3_15n, NZ, NXMEM, NYMEM);
# endif
#endif /* NITRATE */
#ifdef DIC
		if (flags[21]) set_darray3d_zero(mn_dic, NZ, NXMEM, NYMEM);
		if (flags[31]) set_darray3d_zero(mn_alk, NZ, NXMEM, NYMEM);
		if (flags[48]) set_darray3d_zero(mn_jdic, NZ, NXMEM, NYMEM);
#endif

		printf("netcdf record = %d\n", nrec + 1);
		nrec++;

		} /*  end if nmn==WRITEINT */

	    /*-------------------------------------------------------------------*
	     *
	     *     end integration loop
	     *
	     *-------------------------------------------------------------------*/

	    } /* end itts loop over NTSTEP */
	} /* end while */

    //BX-a
    /*-----------------------------------
     *
     *     write restart file
     *
     *----------------------------------*/

    //   First write the up-to-date tracer values to the mean fields.
    //   In order to save memory the instantaneous values for the
    //   restart will be written on mean fields in the restart file.

    printf("start of array copying\n");

#ifdef ENTRAIN
    copy_fix_darray2d(mn_eaml, eaml);
#endif
#ifdef ENTRAIN
    copy_darray3d(mn_ea, ea, NZ, NXMEM, NYMEM);
    //    copy_darray3d(mn_eb, Salttm, NZ, NXMEM, NYMEM);
#endif
    //HF #ifdef SMOOTH_REST
    copy_fix_darray3d(mn_h, h, NZ, NXMEM, NYMEM);
    //HF #endif
#ifdef OXYGEN
    copy_darray3d(mn_oxygen, tr[mOXYGEN], NZ, NXMEM, NYMEM);
#endif
#ifdef PHOSPHATE
    copy_darray3d(mn_phos, tr[mPHOSPHATE], NZ, NXMEM, NYMEM);
    copy_darray3d(mn_dop, tr[mDOP], NZ, NXMEM, NYMEM);
#endif
#ifdef NITRATE
    copy_darray3d(mn_nitr, tr[mNITRATE], NZ, NXMEM, NYMEM);
    copy_darray3d(mn_don, tr[mDON], NZ, NXMEM, NYMEM);
# ifdef N15_CYCLE
    copy_darray3d(mn_nitr_15n, tr[mNITRATE_15n], NZ, NXMEM, NYMEM);
    copy_darray3d(mn_don_15n, tr[mDON_15n], NZ, NXMEM, NYMEM);
# endif
#endif /* NITRATE */
#ifdef DIC
    for (k=0;k<NZ;k++)
	for (i=0;i<NXMEM;i++)
	    for (j=0;j<NYMEM;j++)
		{
		mn_dic[k][i][j] = dic[k][i][j];
		mn_alk[k][i][j] = alk[k][i][j];
		}
#endif
#ifdef AGE
    copy_darray3d(mn_age, age, NZ, NXMEM, NYMEM);
#endif
    for (k = 0; k < NZ; k++)
	{
	for (i = 0; i < NXMEM; i++)
	    {
	    for (j = 0; j < NYMEM; j++)
		{
#ifdef INERT
		mn_inert1[k][i][j] = inert1[k][i][j];
		mn_inert2[k][i][j] = inert2[k][i][j];
#endif
#ifdef OXY18
		mn_o18[k][i][j] = o18[k][i][j];
#endif
#ifdef CFC
		mn_cfc11[k][i][j] = cfc11[k][i][j];
		mn_cfc12[k][i][j] = cfc12[k][i][j];
#endif
		// begin ashao
#ifdef SF6
		mn_sf6[k][i][j] = sf6[k][i][j];
#endif
#ifdef RADIOACTIVE
		mn_i131[k][i][j] = i131[k][i][j];
		mn_cs134[k][i][j] = cs134[k][i][j];
		mn_cs136[k][i][j] = cs136[k][i][j];
		mn_cs137[k][i][j] = cs137[k][i][j];
		mn_ba140[k][i][j] = ba140[k][i][j];
		mn_la140[k][i][j] = la140[k][i][j];
#endif
#ifdef TTD_BP
		mn_ttd[k][i][j] = ttd[k][i][j];
#endif
		// end ashao
		} /* for j loop */
	    } /* for i loop */
	} /* for k loop */

    //  Second, create restart file name and open file

    sprintf(restart_filename, "restart.%s.%04d.nc", run_name, cmon);
    printf("Writing NetCDF restart file '%s'.\n\n", restart_filename);

    /* Copy the variable descriptions to a list of the actual restart variables. */
    nvar = 0;
    for (i = 0; i < NOVARS; i++)
	if (rflags[i] > 0)
	    {
	    var_out[nvar] = vars[i];
	    varmap[i] = nvar;
	    nvar++;
	    }

    // do NOT force float precision output with last argument
    create_file(restart_filename, NETCDF_FILE, var_out, nvar, &fn, &cdfid,
	    timeid, varinfo, 0);

    for (m = 0; m < NOVARS; m++)
	if (rflags[m])
	    {
	    err = write_field(cdfid, &fn, vars[m], varinfo[varmap[m]], 0,
		    var[m]);
	    if (err == -1)
		printf("Error writing %s.\n", vars[m].name);
	    }

    close_file(&cdfid, &fn);

    printf("\n programme termine normallement. \n");
    return (0);
    }

/*-------------------------------------------------------------------*
 *
 *     begin subroutines
 *                                                                   *
 *-------------------------------------------------------------------*/

#ifdef WRTTS
void write_ts(double wrts)
    {
    int i,j,k,m;
    int status;
    char message[144];
    int err;
    double *wrdy;
#ifdef WRTTS
    wrdy = malloc(sizeof(double));
#endif
    /********** set means to zero                                   */
    /*  2d variables */
    if (flags[8]) set_fix_darray2d_zero(mn_eaml);
#ifdef CFC
    if (flags[16]) set_fix_darray2d_zero(mn_c11flux);
    if (flags[17]) set_fix_darray2d_zero(mn_c12flux);
#endif
    // begin ashao
#ifdef SF6
    if (flags[51]) set_fix_darray2d_zero(mn_sf6flux);
#endif
    // end ashao
#ifdef PHOSPHATE
    if (flags[24]) set_fix_darray2d_zero(mn_pflux);
#endif
#ifdef OXYGEN
    if (flags[13]) set_fix_darray2d_zero(mn_oxyflux);
#endif
#ifdef DIC
    if (flags[32]) set_fix_darray2d_zero(mn_co2flux);
    if (flags[49]) set_fix_darray2d_zero(mn_pco2s);
#endif
#ifdef NITRATE
    if (flags[42]) set_fix_darray2d_zero(mn_nflux);
# ifdef N15_CYCLE
    if (flags[44]) set_fix_darray2d_zero(mn_nflux_15n);
# endif
#endif
    /*---------------------------------------
     *
     *     write tracer values on output field
     *
     *--------------------------------------*/

    printf("write output field\n");

    if (flags[8]) add_fix_darray2d(mn_eaml, eaml);
#ifdef PHOSPHATE
    if (flags[24]) add_fix_darray2d(mn_pflux, flux_pop);
#endif
#ifdef CFC
    if (flags[16]) add_fix_darray2d(mn_c11flux, c11flux);
    if (flags[17]) add_fix_darray2d(mn_c12flux, c12flux);
#endif
    // begin ashao
#ifdef SF6
    if (flags[51]) add_fix_darray2d(mn_sf6flux, sf6flux);
#endif
    // end ashao
#ifdef DIC
    if (flags[32]) add_fix_darray2d(mn_co2flux, co2flux);
    if (flags[49]) add_fix_darray2d(mn_pco2s, pco2s);
#endif
#ifdef NITRATE
    if (flags[42]) add_fix_darray2d(mn_nflux, flux_pon);
# ifdef N15_CYCLE
    if (flags[44]) add_fix_darray2d(mn_nflux_15n, flux_pon_15n);
# endif
#endif
    //HF outer if, inner for loops, instead of vice versa!
    //HF think about writing appropriate subroutine(s)!!!
    if (flags[1]) copy_darray3d(mn_u, u, NZ, NXMEM, NYMEM);
    /*{
     for (k=0;k<NZ;k++)
     for (i=0;i<NXMEM;i++)
     for (j=0;j<NYMEM;j++)
     mn_u[k][i][j] = u[k][i][j];
     } */
    if (flags[2]) copy_darray3d(mn_v, v, NZ, NXMEM, NYMEM);
    /* {
     for (k=0;k<NZ;k++)
     for (i=0;i<NXMEM;i++)
     for (j=0;j<NYMEM;j++)
     mn_v[k][i][j] = v[k][i][j];
     }*/
    if (flags[3])
	{
	for (k=0;k<NZ;k++)
	    for (i=0;i<NXMEM;i++)
		for (j=0;j<NYMEM;j++)
		    mn_h[k][i][j] = h[k][i][j];
	}
    if (flags[4]) copy_darray3d(mn_uhtm, uhtm, NZ, NXMEM, NYMEM);
    /*{
     for (k=0;k<NZ;k++)
     for (i=0;i<NXMEM;i++)
     for (j=0;j<NYMEM;j++)
     mn_uhtm[k][i][j] = uhtm[k][i][j];
     }*/
    if (flags[5]) copy_darray3d(mn_vhtm, vhtm, NZ, NXMEM, NYMEM);
    /*{
     for (k=0;k<NZ;k++)
     for (i=0;i<NXMEM;i++)
     for (j=0;j<NYMEM;j++)
     mn_vhtm[k][i][j] = vhtm[k][i][j];
     }*/
    if (flags[6])
	{
	for (k=0;k<NZ;k++)
	    for (i=0;i<NXMEM;i++)
		for (j=0;j<NYMEM;j++)
		    mn_ea[k][i][j] = ea[k][i][j];
	}
    if (flags[7])
	{
	for (k=0;k<NZ;k++)
	    for (i=0;i<NXMEM;i++)
		for (j=0;j<NYMEM;j++)
		    mn_eb[k][i][j] = Salttm[k][i][j];
	}
#ifdef AGE
    if (flags[9]) copy_darray3d(mn_age, age, NZ, NXMEM, NYMEM);
# if defined (AGE2) || defined (AGE3)
    if (flags[30]) copy_darray3d(mn_hnew, hnew, NZ, NXMEM, NYMEM);
# endif
#endif
    if (flags[18])
	{
	for (k=0;k<2;k++)
	    for (i=0;i<NXMEM;i++)
		for (j=0;j<NYMEM;j++)
		    mn_rml[k][i][j] = rml[k][i][j];
	}
#ifdef OXYGEN
    if (flags[10])
	{
	for (k=0;k<NZ;k++)
	    for (i=0;i<NXMEM;i++)
		for (j=0;j<NYMEM;j++)
		    mn_oxygen[k][i][j] = tr[mOXYGEN][k][i][j];
	}
    if (flags[11])
	{
	for (k=0;k<NZ;k++)
	    for (i=0;i<NXMEM;i++)
		for (j=0;j<NYMEM;j++)
		    mn_o2sat[k][i][j] = o2_sat[k][i][j];
	}
    if (flags[12])
	{
	for (k=0;k<NZ;k++)
	    for (i=0;i<NXMEM;i++)
		for (j=0;j<NYMEM;j++)
		    mn_jo2[k][i][j] = jo2[k][i][j];
	}
    if (flags[13])
	{
	for (i=0;i<NXMEM;i++)
	    for (j=0;j<NYMEM;j++)
		mn_oxyflux[i][j] = oxyflux[i][j];
	}
#endif /* OXYGEN */
#ifdef CFC
    if (flags[14])
	{
	for (k=0;k<NZ;k++)
	    for (i=0;i<NXMEM;i++)
		for (j=0;j<NYMEM;j++)
		    mn_cfc11[k][i][j] = cfc11[i][j];
	}
    if (flags[15])
	{
	for (k=0;k<NZ;k++)
	    for (i=0;i<NXMEM;i++)
		for (j=0;j<NYMEM;j++)
		    mn_cfc12[k][i][j] = cfc12[i][j];
	}
#endif
    // begin ashao
#ifdef SF6
    if (flags[50])
	{
	for (k=0;k<NZ;k++)
	    for (i=0;i<NXMEM;i++)
		for (j=0;j<NYMEM;j++)
		    mn_sf6[k][i][j] = sf6[i][j];
	}
#endif
#ifdef RADIOACTIVE
    for (k=0;k<NZ;k++)
	for (i=0;i<NXMEM;i++)
	    for (j=0;j<NYMEM;j++)
		{
		mn_i131[k][i][j] = i131[i][j];
		mn_cs134[k][i][j] = cs134[i][j];
		mn_cs136[k][i][j] = cs136[i][j];
		mn_cs137[k][i][j] = cs137[i][j];
		mn_ba140[k][i][j] = ba140[i][j];
		mn_la140[k][i][j] = la140[i][j];
		}
#endif
    // end ashao
#ifdef PHOSPHATE
    if (flags[19])
	{
	for (k=0;k<NZ;k++)
	    for (i=0;i<NXMEM;i++)
		for (j=0;j<NYMEM;j++)
		    mn_phos[k][i][j] = tr[mPHOSPHATE][k][i][j];
	}
    if (flags[22])
	{
	for (k=0;k<NZ;k++)
	    for (i=0;i<NXMEM;i++)
		for (j=0;j<NYMEM;j++)
		    mn_dop[k][i][j] = tr[mDOP][k][i][j];
	}
    if (flags[27])
	{
	for (k=0;k<NZ;k++)
	    for (i=0;i<NXMEM;i++)
		for (j=0;j<NYMEM;j++)
		    mn_pobs[k][i][j] = po4_star_lay[k][i][j];
	if (flags[33])
	    {
	    for (k=0;k<NZ;k++)
		for (i=0;i<NXMEM;i++)
		    for (j=0;j<NYMEM;j++)
			mn_jpo4[k][i][j] = jpo4[k][i][j];
	    }
# ifdef PROGNOSTIC
	if (flags[35]) copy_darray3d(mn_lightlim, lightlim, NZ, NXMEM, NYMEM);
#  ifdef NITRATE
	if (flags[36])
	    {
	    for (k=0;k<NZ;k++)
		for (i=0;i<NXMEM;i++)
		    for (j=0;j<NYMEM;j++)
			mn_nitrlim[k][i][j] = nitrlim[k][i][j];
	    }
#  endif
	if (flags[37]) copy_darray3d(mn_ironlim, ironlim, NZ, NXMEM, NYMEM);
# endif /* PROGNOSTIC */
#endif /* PHOSPHATE */
#ifdef NITRATE
	if (flags[20])
	    {
	    for (k=0;k<NZ;k++)
		for (i=0;i<NXMEM;i++)
		    for (j=0;j<NYMEM;j++)
			mn_nitr[k][i][j] = tr[mNITRATE][k][i][j];
	    }
	if (flags[23])
	    {
	    for (k=0;k<NZ;k++)
		for (i=0;i<NXMEM;i++)
		    for (j=0;j<NYMEM;j++)
			mn_don[k][i][j] = tr[mDON][k][i][j];
	    }
	if (flags[38])
	    {
	    for (k=0;k<NZ;k++)
		for (i=0;i<NXMEM;i++)
		    for (j=0;j<NYMEM;j++)
			mn_jno3[k][i][j] = jno3[k][i][j];
	    }
	if (flags[43])
	    {
	    for (k=0;k<NZ;k++)
		for (i=0;i<NXMEM;i++)
		    for (j=0;j<NYMEM;j++)
			mn_nobs[k][i][j] = no3_lay[k][i][j];
	    if (flags[45])
		{
		for (k=0;k<NZ;k++)
		    for (i=0;i<NXMEM;i++)
			for (j=0;j<NYMEM;j++)
			    mn_n2fix[k][i][j] = jn2fix[k][i][j];
		}
	    if (flags[46])
		{
		for (k=0;k<NZ;k++)
		    for (i=0;i<NXMEM;i++)
			for (j=0;j<NYMEM;j++)
			    mn_denitw[k][i][j] = jdenitw[k][i][j];
		}
	    if (flags[47])
		{
		for (k=0;k<NZ;k++)
		    for (i=0;i<NXMEM;i++)
			for (j=0;j<NYMEM;j++)
			    mn_denits[k][i][j] = jdenits[k][i][j];
		}
# ifdef N15_CYCLE
	    if (flags[39])
		{
		for (k=0;k<NZ;k++)
		    for (i=0;i<NXMEM;i++)
			for (j=0;j<NYMEM;j++)
			    mn_nitr_15n[k][i][j] = tr[mNITRATE_15n][k][i][j];
		}
	    if (flags[40])
		{
		for (k=0;k<NZ;k++)
		    for (i=0;i<NXMEM;i++)
			for (j=0;j<NYMEM;j++)
			    mn_don_15n[k][i][j] = tr[mDON_15n][k][i][j];
		}
	    if (flags[41])
		{
		for (k=0;k<NZ;k++)
		    for (i=0;i<NXMEM;i++)
			for (j=0;j<NYMEM;j++)
			    mn_jno3_15n[k][i][j] = jno3_15n[k][i][j];
		}
# endif
#endif /* NITRATE */
#ifdef DIC
	    if (flags[21])
		{
		for (k=0;k<NZ;k++)
		    for (i=0;i<NXMEM;i++)
			for (j=0;j<NYMEM;j++)
			    mn_dic[k][i][j] = dic[k][i][j];
		}
	    if (flags[31])
		{
		for (k=0;k<NZ;k++)
		    for (i=0;i<NXMEM;i++)
			for (j=0;j<NYMEM;j++)
			    mn_alk[k][i][j] = alk[k][i][j];
		}
	    if (flags[48]) copy_darray3d(mn_jdic, jdic, NZ, NXMEM, NYMEM);
#endif
#ifdef OXY18
	    if (flags[25])
		{
		for (k=0;k<NZ;k++)
		    for (i=0;i<NXMEM;i++)
			for (j=0;j<NYMEM;j++)
			    mn_o18[k][i][j] = delo18[k][i][j];
		}
	    if (flags[26])
		{
		for (k=0;k<NZ;k++)
		    for (i=0;i<NXMEM;i++)
			for (j=0;j<NYMEM;j++)
			    mn_jo18[k][i][j] = jo18[k][i][j];
		}
#endif
#ifdef INERT
	    if (flags[28])
		{
		for (k=0;k<NZ;k++)
		    for (i=0;i<NXMEM;i++)
			for (j=0;j<NYMEM;j++)
			    mn_inert1[k][i][j] = inert1[k][i][j];
		}
	    if (flags[29])
		{
		for (k=0;k<NZ;k++)
		    for (i=0;i<NXMEM;i++)
			for (j=0;j<NYMEM;j++)
			    mn_inert2[k][i][j] = inert2[k][i][j];
		}
#endif /* INERT */
	    printf("Writing variables for sub timestep %i out to netCDF\n\n",itts);
	    //  open netcdf file for each writing
	    status = nc_open(output_filename, NC_WRITE, &cdfid);
	    if (status != NC_NOERR)
		{
		strcpy(message,"Opening file"); strcat(message,output_filename);
		handle_error(message,status);
		}
	    *wrdy = wrts;
	    err = write_time(cdfid,fn,timeid[0],nrec, wrdy);
	    if (err == -1) printf("Error writing day.\n");

	    for (m=0;m<NOVARS;m++) if (flags[m])
		{
		err = write_field(cdfid, fn, vars[m], varinfo[varmap[m]],
			nrec,var[m]);
		if (err == -1) printf("Error writing %s.\n",vars[m].name);
		}
	    close_file(&cdfid, &fn);

	    printf("netcdf record = %d\n",nrec++);
	    }
#endif //WRTTS
	void alloc_fields(void)
	    {

	    int m;

	    extern double *var[NOVARS];
	    extern double junk[(NZ + 1) * (NXMEM) * (NYMEM)];
	    extern long varsize[NOVARS];
	    extern int flags[NOVARS];

	    for (m = 0; m < NOVARS; m++)
		{

		//HF if ( m>3 )    {
		//HF: added mn_h
		if (m >= 3)
		    {
		    switch (vars[m].z_grid)
			{
		    case 'L':
			varsize[m] = NZ * (NXMEM) * (NYMEM);
			break;
		    case 'i':
			varsize[m] = (NZ + 1) * (NXMEM) * (NYMEM);
			break;
		    case '2':
			varsize[m] = 2 * (NXMEM) * (NYMEM);
			break;
		    case '1':
			varsize[m] = (NXMEM) * (NYMEM);
			break;
		    default:
			printf("Unknown layer axis[%d] %c.\n", m, vars[m].z_grid);
			}

		    var[m] = calloc((size_t) varsize[m], sizeof(double));
		    if (var[m] == NULL)
			{
			printf("Unable to allocate memory for var[%d].\n", m);
			exit(-21);
			}

		    //	temporarily disabling this output
		    //      printf("Allocated memory for var[%d],%s.\n\n",m,vars[m].name);

		    }
		else
		    {
		    var[m] = junk;
		    varsize[m] = 0;
		    }
		}

	    var[0] = &D[0][0];
	    var[1] = &mn_u[0][0][0];
	    var[2] = &mn_v[0][0][0];
	    var[3] = &mn_h[0][0][0];
	    var[4] = &mn_uhtm[0][0][0];
	    var[5] = &mn_vhtm[0][0][0];
	    var[6] = &mn_ea[0][0][0];
	    var[7] = &mn_eb[0][0][0];
	    var[8] = &mn_eaml[0][0];
	    var[18] = &mn_rml[0][0][0];

#ifdef INERT
	    var[28] = &mn_inert1[0][0][0];
	    var[29] = &mn_inert2[0][0][0];
#endif
#ifdef AGE
	    var[9] = &mn_age[0][0][0];
# if defined AGE2 || defined AGE3
	    var[30] = &mn_hnew[0][0][0];
# endif
#endif

#ifdef OXYGEN
	    var[10] = &mn_oxygen[0][0][0];
	    var[11] = &mn_o2sat[0][0][0];
	    var[12] = &mn_jo2[0][0][0];
	    var[13] = &mn_oxyflux[0][0];
#endif
#ifdef OXY18
	    var[25] = &mn_o18[0][0][0];
	    var[26] = &mn_jo18[0][0][0];
#endif

#ifdef CFC
	    var[14] = &mn_cfc11[0][0][0];
	    var[15] = &mn_cfc12[0][0][0];
	    var[16] = &mn_c11flux[0][0];
	    var[17] = &mn_c12flux[0][0];
	    var[58] = &cfc11_sat[0][0];
	    var[59] = &cfc12_sat[0][0];
#endif

	    // begin ashao
#ifdef SF6
	    var[50] = &mn_sf6[0][0][0];
	    var[51] = &mn_sf6flux[0][0];
	    var[60] = &sf6_sat[0][0];
#endif
#ifdef RADIOACTIVE
	    var[52] = &mn_i131[0][0][0];
	    var[53] = &mn_cs134[0][0][0];
	    var[54] = &mn_cs136[0][0][0];
	    var[55] = &mn_cs137[0][0][0];
	    var[56] = &mn_ba140[0][0][0];
	    var[57] = &mn_la140[0][0][0];
#endif
#ifdef TTD_BP
	    var[61] = &mn_ttd[0][0][0];
#endif

	    // end ashao
#ifdef PHOSPHATE
	    var[19] = &mn_phos[0][0][0];
	    var[22] = &mn_dop[0][0][0];
	    var[24] = &mn_pflux[0][0];
	    var[27] = &mn_pobs[0][0][0];
	    var[33] = &mn_jpo4[0][0][0];
# ifdef PROGNOSTIC
	    var[34] = &Temptm[0][0][0];
	    var[35] = &mn_lightlim[0][0][0];
#  ifdef NITRATE
	    var[36] = &mn_nitrlim[0][0][0];
#  endif
	    var[37] = &mn_ironlim[0][0][0];
# endif /* PROGNOSTIC */
#endif /* PHOSPHATE */
#ifdef NITRATE
	    var[20] = &mn_nitr[0][0][0];
	    var[23] = &mn_don[0][0][0];
	    var[38] = &mn_jno3[0][0][0];
	    var[42] = &mn_nflux[0][0];
	    var[43] = &mn_nobs[0][0][0];
	    var[45] = &mn_n2fix[0][0][0];
	    var[46] = &mn_denitw[0][0][0];
	    var[47] = &mn_denits[0][0][0];
# ifdef N15_CYCLE
	    var[39] = &mn_nitr_15n[0][0][0];
	    var[40] = &mn_don_15n[0][0][0];
	    var[41] = &mn_jno3_15n[0][0][0];
	    var[44] = &mn_nflux_15n[0][0];
# endif
#endif /* NITRATE */
#ifdef DIC
	    var[21] = &mn_dic[0][0][0];
	    var[31] = &mn_alk[0][0][0];
	    var[32] = &mn_co2flux[0][0];
	    var[48] = &mn_jdic[0][0][0];
	    var[49] = &mn_pco2s[0][0];
#endif

	    }
