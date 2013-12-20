/********+*********+*********+*********+*********+*********+*********+*
 *         Initialize                                                 *
 *                                                                    *
 ********+*********+*********+*********+*********+*********+*********+*/

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "init.h"
#include "offtrac.h"
#include "alloc.h"
#include "read.h"
#ifdef SMOOTH_REST
#include "step.h"
#endif

extern const double misval;

# ifdef N15_CYCLE
extern const double parm_n15_std_fraction;
# endif

extern void initializemasks(void);

extern double ****tr;
extern double umask[NXMEM][NYMEM];   /* _mask are 1 over ocean and 0  */
extern double vmask[NXMEM][NYMEM];   /* over land on the u & v grids. */

extern double D[NXMEM][NYMEM];
extern double lonh[NXMEM],lath[NYMEM];
extern double dxh[NYMEM],dy;

extern double ***u;
extern double ***v;
extern double h[NZ][NXMEM][NYMEM];

extern double ***uhtm;
extern double ***vhtm;

#ifdef INERT
extern double inert1_init[NZ][NXMEM][NYMEM];
extern double inert1[NZ][NXMEM][NYMEM];
extern int mINERT1;
extern double inert2_init[NZ][NXMEM][NYMEM];
extern double inert2[NZ][NXMEM][NYMEM];
extern int mINERT2;
#endif

#ifdef AGE
extern double age_init[NZ][NXMEM][NYMEM];
//extern double age[NZ][NXMEM][NYMEM];
extern double ***age;
extern int mAGE;
#endif

#ifdef SMOOTH_REST
//BXextern double ***h_init;
extern double h_init[NZ][NXMEM][NYMEM];
extern double hstart[NZ][NXMEM][NYMEM];
double tr_in[NZ];
double z_in[NZ];
double depth[NZ][NXMEM][NYMEM];
double hdepth[NZ][NXMEM][NYMEM];
double tmp_in;
#endif

#ifdef OXYGEN
extern double ***oxy_init;
extern double oxyflux[NXMEM][NYMEM];
extern double our[NZ][NXMEM][NYMEM];
extern int mOXYGEN;
extern double jo2_spec[14][NXMEM][NYMEM];
#endif

#ifdef OXY18
extern double o18_init[NZ][NXMEM][NYMEM];
extern double o18[NZ][NXMEM][NYMEM];
extern int mo18;
#endif

/*   Begin added DT    */
extern int beginyear;
extern int theyear;
extern int lastsave;
extern int estTTD;
/*   End added DT      */

#ifdef CFC
extern double cfc11_init[NZ][NXMEM][NYMEM];
extern double cfc12_init[NZ][NXMEM][NYMEM];
extern double cfc11[NZ][NXMEM][NYMEM];
extern double cfc12[NZ][NXMEM][NYMEM];
extern double c11flux[NXMEM][NYMEM];
extern double c12flux[NXMEM][NYMEM];
extern double c11monflux[NXMEM][NYMEM];
extern double c12monflux[NXMEM][NYMEM];

extern int mCFC11,mCFC12;
#ifdef SAT_TRACER
extern int mCFC11_sat, mCFC12_sat;
#endif
#endif

// begin ashao
#ifdef SF6
extern double sf6_init[NZ][NXMEM][NYMEM];
extern double sf6[NZ][NXMEM][NYMEM];
extern double sf6flux[NXMEM][NYMEM];
extern double sf6monflux[NXMEM][NYMEM];
extern int mSF6;
#ifdef SAT_TRACER
extern int mSF6_sat;
#endif
#endif

#ifdef TTD_BP
extern double ttd[NZ][NXMEM][NYMEM];
extern double ttd_mask[NXMEM][NYMEM];
extern int mttd;
#endif
// end ashao

#ifdef PHOSPHATE
extern double ***phosphate_init;
extern int mPHOSPHATE;
double ***dop_init;
extern int mDOP;
#endif

#ifdef NITRATE
extern double ***nitrate_init;
extern int mNITRATE;
double ***don_init;
extern int mDON;
# ifdef N15_CYCLE
extern double ***nitrate_init_15n;
extern int mNITRATE_15n;
double ***don_init_15n;
extern int mDON_15n;
# endif
#endif /* NITRATE */

#ifdef DIC
double ***dic_init;
extern double dic[NZ][NXMEM][NYMEM];
extern int mDIC;
double ***alk_init;
extern double alk[NZ][NXMEM][NYMEM];
extern int mALK;
#endif

//BX-a for restart reading
#ifdef ENTRAIN
extern double ea_init[NZ][NXMEM][NYMEM];
extern double eaml_init[NXMEM][NYMEM];
#endif
//BX-e

#ifdef SPONGE
extern double Idamp[NXMEM][NYMEM];
extern int num_fields;
extern double sp_e[NZ][NXMEM][NYMEM];
# if defined (AGE) || defined(CFC) || defined(SF6) // ashao
extern double sp_age[NZ][NXMEM][NYMEM];
# endif
# ifdef OXYGEN
extern double h_clim[NZ][NXMEM][NYMEM];
extern double sp_oxygen[NZ][NXMEM][NYMEM];
# endif
# ifdef OXY18
extern double sp_o18[NZ][NXMEM][NYMEM];
# endif
# ifdef INERT
extern double sp_inert1[NZ][NXMEM][NYMEM];
extern double sp_inert2[NZ][NXMEM][NYMEM];
# endif
# ifdef PHOSPHATE
extern double sp_phosphate[NZ][NXMEM][NYMEM];
extern double sp_dop[NZ][NXMEM][NYMEM];
# endif
# ifdef NITRATE
extern double sp_nitrate[NZ][NXMEM][NYMEM];
extern double sp_don[NZ][NXMEM][NYMEM];
#  ifdef N15_CYCLE
extern double sp_nitrate_15n[NZ][NXMEM][NYMEM];
extern double sp_don_15n[NZ][NXMEM][NYMEM];
#  endif
# endif /* NITRATE */
# ifdef DIC
extern double sp_dic[NZ][NXMEM][NYMEM];
# endif
#endif

#ifdef RESTART
void initialize(int imon, char *run_name)
#else
void initialize(int imon)
#endif
{
/* Arguments: mon - Time of the start of the run in months.             */
/* Arguments: inmon - Time of the restart of the run in months.         */

    int i, j, k, m;

#ifdef SPONGE

    num_fields=1;  /* that e thing */
# ifdef AGE
    num_fields++;
# endif
# ifdef OXYGEN
    num_fields++;
    num_fields++;
# endif
# ifdef CFC
    num_fields++;
    num_fields++;
# endif
// begin ashao
# ifdef SF6
	num_fields++;
# endif

// end ashao
# ifdef PHOSPHATE
    num_fields++;
    num_fields++;
# endif
# ifdef NITRATE
    num_fields++;
    num_fields++;
#  ifdef N15_CYCLE
    num_fields++;
    num_fields++;
#  endif
# endif /* NITRATE */
# ifdef DIC
    num_fields++;
# endif

    printf("num_fields = %d\n\n",num_fields);
    read_sponge();

    initialize_sponge(Idamp,num_fields);
    set_up_sponge_field(sp_e,NULL,NZ);
 
#endif  /* SPONGE */

    initializemasks();   /* create land/ocean masks based on topography array (D) */
#if defined JO2_SPEC || defined UNIVENT
    read_misc_fields();  /* reading miscellaneous fields */
#endif

    m=0; 

    //  printf("Initialize 173 for month %i.\n",imon);
    //  printf("Initialize 174 for run %s.\n",run_name);

#ifdef PHOSPHATE
    dop_init = alloc3d(NZ,NXMEM,NYMEM);
#endif
#ifdef NITRATE
    don_init = alloc3d(NZ,NXMEM,NYMEM);
# ifdef N15_CYCLE
    don_init_15n = alloc3d(NZ,NXMEM,NYMEM);
# endif
#endif /* NITRATE */
#ifdef DIC
    dic_init = alloc3d(NZ,NXMEM,NYMEM);
    alk_init = alloc3d(NZ,NXMEM,NYMEM);
#endif

#ifdef RESTART
    read_tracer_init(imon,run_name);  /* reading all tracer initial values from restart file*/
#else
    /*  for (i=0;i<=NXMEM-1;i++)
    for (j=0;j<=NYMEM-1;j++)
    printf("l202 ini D[%d][%d]=%g\n",i,j,D[i][j]); */

    read_tracer_init(imon);  /* reading all tracer initial values */
#endif
    // printf("Initialize 173 for month %i.\n",imon);

#ifdef SMOOTH_REST
    //BX interpolate restart files to new depth field associated 
    //BX with time step wrap around at the end of a climatological forcing file
    printf("Restart using smoothed transition of depth fields. \n");
    z_depth(h_init,depth);
    z_depth(hstart,hdepth);
    //printf("Initialize 199,   h_init= %g\n", h_init[20][90][146]);
    //printf("Initialize 199,   hstart= %g\n", hstart[20][90][146]);
    //printf("Initialize 199,    depth= %g\n", depth[20][90][146]);
    //printf("Initialize 199,   hdepth= %g\n", hdepth[20][90][146]);
    for (i=X1;i<=nx;i++) {
      for (j=Y1;j<=ny;j++) {
	  if (D[i][j]>MINIMUM_DEPTH) {
	      for (k=0;k<NZ;k++)
		  z_in[k]  = depth[k][i][j];
# ifdef ENTRAIN
	      for (k=0;k<NZ;k++) tr_in[k] = ea_init[k][i][j];
	      for (k=0;k<NZ;k++) {
		  tmp_in = lin_interpp(hdepth[k][i][j], tr_in, 
		  				z_in, 0, NZ);
		  if (tmp_in < 0.e0) tmp_in = 0.;
		  ea_init[k][i][j] = tmp_in;
	      }
	      //BX not action neccessary for eaml_init - is only 2D
# endif
# ifdef OXYGEN
	      for (k=0;k<NZ;k++) tr_in[k] = oxy_init[k][i][j];
	      for (k=0;k<NZ;k++) {
		  tmp_in = lin_interpp(hdepth[k][i][j], tr_in, 
		  				z_in, 0, NZ);
		  if (tmp_in < 0.e0) tmp_in = 0.;
		  oxy_init[k][i][j] = tmp_in;
	      }
#  ifdef OXY18
	      for (k=0;k<NZ;k++) tr_in[k] = o18_init[k][i][j];
	      for (k=0;k<NZ;k++) {
		  tmp_in = lin_interpp(hdepth[k][i][j], tr_in, 
		  				z_in, 0, NZ);
		  if (tmp_in < 0.e0) tmp_in = 0.;
		  o18_init[k][i][j] = tmp_in;
	      }
#  endif
# endif
# ifdef PHOSPHATE
	      for (k=0;k<NZ;k++) tr_in[k] = phosphate_init[k][i][j];
	      for (k=0;k<NZ;k++) {
		  tmp_in = lin_interpp(hdepth[k][i][j], tr_in, 
		  				z_in, 0, NZ);
		  if (tmp_in < 0.e0) tmp_in = 0.;
		  phosphate_init[k][i][j] = tmp_in;
	      }
	      for (k=0;k<NZ;k++) tr_in[k] = dop_init[k][i][j];
	      for (k=0;k<NZ;k++) {
		  tmp_in = lin_interpp(hdepth[k][i][j], tr_in, 
		  				z_in, 0, NZ);
		  if (tmp_in < 0.e0) tmp_in = 0.;
		  dop_init[k][i][j] = tmp_in;
	      }
# endif
# ifdef NITRATE
	      for (k=0;k<NZ;k++) tr_in[k] = nitrate_init[k][i][j];
	      for (k=0;k<NZ;k++) {
		  tmp_in = lin_interpp(hdepth[k][i][j], tr_in, 
		  				z_in, 0, NZ);
		  if (tmp_in < 0.e0) tmp_in = 0.;
		  nitrate_init[k][i][j] = tmp_in;
	      }
	      for (k=0;k<NZ;k++) tr_in[k] = don_init[k][i][j];
	      for (k=0;k<NZ;k++) {
		  tmp_in = lin_interpp(hdepth[k][i][j], tr_in, 
		  				z_in, 0, NZ);
		  if (tmp_in < 0.e0) tmp_in = 0.;
		  don_init[k][i][j] = tmp_in;
	      }
# endif /* NITRATE */
# ifdef INERT
	      for (k=0;k<NZ;k++) tr_in[k] = inert1_init[k][i][j];
	      for (k=0;k<NZ;k++) {
		  tmp_in = lin_interpp(hdepth[k][i][j], tr_in, 
		  				z_in, 0, NZ);
		  if (tmp_in < 0.e0) tmp_in = 0.;
		  inert1_init[k][i][j] = tmp_in;
	      }
	      for (k=0;k<NZ;k++) tr_in[k] = inert2_init[k][i][j];
	      for (k=0;k<NZ;k++) {
		  tmp_in = lin_interpp(hdepth[k][i][j], tr_in, 
		  				z_in, 0, NZ);
		  if (tmp_in < 0.e0) tmp_in = 0.;
		  inert2_init[k][i][j] = tmp_in;
	      }
# endif
# ifdef AGE
	      for (k=0;k<NZ;k++) tr_in[k] = age_init[k][i][j];
	      for (k=0;k<NZ;k++) {
		  tmp_in = lin_interpp(hdepth[k][i][j], tr_in, 
		  				z_in, 0, NZ);
		  if (tmp_in < 0.e0) tmp_in = 0.;
		  age_init[k][i][j] = tmp_in;
	      }
# endif
# ifdef DIC
	      for (k=0;k<NZ;k++) tr_in[k] = dic_init[k][i][j];
	      for (k=0;k<NZ;k++) {
		  tmp_in = lin_interpp(hdepth[k][i][j], tr_in, 
		  				z_in, 0, NZ);
		  //if(i==90 && j==146) {
		  //	    printf("Initialize 220, k= %i,tmp_in= %g\n",k, tmp_in);
		  //	}
		  if (tmp_in < 0.e0) tmp_in = 0.;
		  dic_init[k][i][j] = tmp_in;
	      }
	      for (k=0;k<NZ;k++) tr_in[k] = alk_init[k][i][j];
	      for (k=0;k<NZ;k++) {
		  tmp_in = lin_interpp(hdepth[k][i][j], tr_in, 
		  				z_in, 0, NZ);
		  if (tmp_in < 0.e0) tmp_in = 0.;
		  alk_init[k][i][j] = tmp_in;
	      }
# endif
# ifdef CFC
	      for (k=0;k<NZ;k++) tr_in[k] = cfc11_init[k][i][j];
	      for (k=0;k<NZ;k++) {
		  tmp_in = lin_interpp(hdepth[k][i][j], tr_in, 
		  				z_in, 0, NZ);
		  if (tmp_in < 0.e0) tmp_in = 0.;
		  cfc11_init[k][i][j] = tmp_in;
	      }
	      for (k=0;k<NZ;k++) tr_in[k] = cfc12_init[k][i][j];
	      for (k=0;k<NZ;k++) {
		  tmp_in = lin_interpp(hdepth[k][i][j], tr_in, 
		  				z_in, 0, NZ);
		  if (tmp_in < 0.e0) tmp_in = 0.;
		  cfc12_init[k][i][j] = tmp_in;
	      }
# endif
// begin ashao
# ifdef SF6
	      for (k=0;k<NZ;k++) tr_in[k] = sf6_init[k][i][j];
	      for (k=0;k<NZ;k++) {
		  tmp_in = lin_interpp(hdepth[k][i][j], tr_in, 
		  				z_in, 0, NZ);
		  if (tmp_in < 0.e0) tmp_in = 0.;
		  sf6_init[k][i][j] = tmp_in;
	      }
# endif


//end ashao
	  } else {   //BX if(D>MINIMUM_DEPTH)
#ifdef WE_DONT_USE_THIS
	      for (k=0;k<NZ;k++) {
# ifdef ENTRAIN
		  ea_init[k][i][j] = misval;
# endif
# ifdef OXYGEN
		  oxy_init[k][i][j] = misval;
#  ifdef OXY18
		  o18_init[k][i][j] = misval;
#  endif
# endif
# ifdef PHOSPHATE
		  phosphate_init[k][i][j] = misval;
		  dop_init[k][i][j] = misval;
# endif
# ifdef NITRATE
		  nitrate_init[k][i][j] = misval;
		  don_init[k][i][j] = misval;
#  ifdef N15_CYCLE
		  nitrate_init_15n[k][i][j] = misval;
		  don_init_15n[k][i][j] = misval;
#  endif
# endif /* NITRATE */
# ifdef INERT
		  inert1_init[k][i][j] = misval;
		  inert2_init[k][i][j] = misval;
# endif
# ifdef AGE
		  age_init[k][i][j] = misval;
# endif
# ifdef DIC
		  dic_init[k][i][j] = misval;
		  alk_init[k][i][j] = misval;
# endif
# ifdef CFC
		  cfc11_init[k][i][j] = misval;
		  cfc12_init[k][i][j] = misval;
# endif
// begin ashao
# ifdef SF6
	sf6_init[k][i][j] = misval;
# endif
// end ashao
	      }
#endif /* NOT USED */
	  }
      }
  }
#endif /* SMOOTH_REST */
#ifdef N15_CYCLE
		for (i=0; i<NXMEM; i++)
		  for (j=0; j<NYMEM; j++)
		    if (fabs(nitrate_init[0][i][j] - misval) > 1.e-3)
		      for (k=0; k<NZ; k++) {
			nitrate_init_15n[k][i][j] = nitrate_init[k][i][j] *
			  parm_n15_std_fraction * 1.005;
			//printf("no3_init_15n(%d,%d,%d)=%g,no3_init=%g\n",
			//     k,i,j,nitrate_init_15n[k][i][j],
			//     nitrate_init[k][i][j]);
		      }
		    else
		      for (k=0; k<NZ; k++)
			nitrate_init_15n[k][i][j] = misval;

		for (i=0; i<NXMEM; i++)
		  for (j=0; j<NYMEM; j++)
		    if (fabs(don_init[0][i][j] - misval) > 1.e-3)
		      for (k=0; k<NZ; k++)
			don_init_15n[k][i][j] = don_init[k][i][j] *
			  parm_n15_std_fraction;
		    else
		      for (k=0; k<NZ; k++)
			don_init_15n[k][i][j] = misval;
#endif /* N15_CYCLE */

    //BX printf("Initialize 239, dic_init= %g\n", dic_init[0][90][146]);

#ifdef INERT
    for (i=0;i<=NXMEM-1;i++) {
	for (j=0;j<=NYMEM-1;j++) {
	    for (k=0;k<NZ;k++) {
		if (D[i][j]>MINIMUM_DEPTH) {
		    inert1[k][i][j] = inert1_init[k][i][j];
		    tr[m][k][i][j] = inert1_init[k][i][j];
		    inert2[k][i][j] = inert2_init[k][i][j];
		    tr[m+1][k][i][j] = inert2_init[k][i][j];
		} else {
		    inert1[k][i][j] = misval;
		    tr[m][k][i][j] = misval;
		    inert2[k][i][j] = misval;
		    tr[m+1][k][i][j] = misval;
		}
	    }                            
	}
    }

# ifdef SPONGE
    set_up_sponge_field(sp_inert1,inert1,NZ);
    set_up_sponge_field(sp_inert2,inert2,NZ);
# endif

    mINERT1=m;
    mINERT2=m+1;
    printf("index of INERT1 tracer m = %d\n\n",mINERT1);
    printf("index of INERT2 tracer m = %d\n\n",mINERT2);
    m++;m++;
#endif /* INERT */

#ifdef AGE

# ifdef TRNCINIT

#  ifdef AGEEXP
  for (i=0;i<=NXMEM-1;i++) {
    for (j=0;j<=NYMEM-1;j++) {
      for (k=0;k<NZ;k++) {
        if (D[i][j]>MINIMUM_DEPTH) {
	    if (k==39 && lath[j] < -50.) {
		age[k][i][j] = 1.;
		tr[m][k][i][j] = 1.;
	    } else {
		age[k][i][j] = 0.;
		tr[m][k][i][j] = 0.;
	    }
        } else {
          age[k][i][j] = misval;
          tr[m][k][i][j] = misval;
        }
      }
    }
  }
#  else  // AGEEXP
#   ifdef AGE2
  for (i=0;i<=NXMEM-1;i++) {
    for (j=0;j<=NYMEM-1;j++) {
      for (k=0;k<NZ;k++) {
        if (D[i][j]>MINIMUM_DEPTH) {
	    age[k][i][j] = 1.;
	    tr[m][k][i][j] = 1.;
        } else {
	    age[k][i][j] = misval;
	    tr[m][k][i][j] = misval;
        }
      }
    }
  }
#   else  // AGE1
#    ifdef AGE1
  for (i=0;i<=NXMEM-1;i++) {
    for (j=0;j<=NYMEM-1;j++) {
      for (k=0;k<NZ;k++) {
        if (D[i][j]>MINIMUM_DEPTH) {
	    //BX for the Eq.Pac.:
	    //BX if (i>=128 && i<=132 && j>=91 && j<=101 && k>=8) {
	    //BX if (i>=128 && i<=132 && j>=91 && j<=101 && k==0) {
	    //BX for the Central NPac.:
	    //BX	    if (i>=128 && i<=132 && j>=151 && j<=155 && k>=8) {
	    //BX for south of Greenland: if (i>=258 && i<=262 && j>=161 && j<=165 && k==0) {

//	    if (j>=0 && j<=NYMEM-1 && i>=90 && i<=100 && k==0) {
	    if (i>=0 && i<=NXMEM-1 && j>=80 && j<=90 && k==0) {
		age[k][i][j] = 1.;
		tr[m][k][i][j] = 1.;
	    } else {
		age[k][i][j] = 0.;
		tr[m][k][i][j] = 0.;
	    }
        } else {
	    age[k][i][j] = misval;
	    tr[m][k][i][j] = misval;
        }

      }
    }
  }
#    else  // AGE2
    printf("initialize tracer from netcdf file\n\n");

    //  read_age_tracer_init();
    for (i=0;i<=NXMEM-1;i++) {
	for (j=0;j<=NYMEM-1;j++) {
	    for (k=0;k<NZ;k++) {
		if (D[i][j]>MINIMUM_DEPTH) {
		    age[k][i][j] = age_init[k][i][j];  /* in years */
		    tr[m][k][i][j] = age_init[k][i][j]*365.0*86400.0; /* in seconds */
		} else {
		    age[k][i][j] = misval;
		    tr[m][k][i][j] = misval;
		}
	    }
	}
    }
#    endif  // AGE2
#   endif  // AGE1
#  endif  // AGEEXP

# else  // TRNCINIT
  
    printf("initialize tracer internally\n\n");

    for (i=0;i<=NXMEM-1;i++) {
	for (j=0;j<=NYMEM-1;j++) {
	    for (k=0;k<NZ;k++) {
		age_init[k][i][j] = 10.0;
		uhtm[k][i][j] = 0.0;
		vhtm[k][i][j] = 0.0;
		u[k][i][j] = 0.0;
		v[k][i][j] = 0.0;
	    }
	}
    }


    for (i=0;i<=NXMEM-1;i++) {
      for (j=0;j<=NYMEM-1;j++) {
	if (D[i][j] > MINIMUM_DEPTH) {
	  for (k=0; k<NZ; k++) {
	    age[k][i][j] = age_init[k][i][j];  /* in years */
	    tr[m][k][i][j] = age_init[k][i][j]*365.0*86400.0; /* in seconds */
	  }
	} else {
	  for (k=0; k<NZ; k++) {
	    age[k][i][j] = misval;
	    tr[m][k][i][j] = misval;
	  }                            
	}
      }
    }
    
# endif  // TRNCINIT

# ifdef SPONGE
    set_up_sponge_field(sp_age,age,NZ);
# endif

    mAGE=m;
    printf("index of AGE tracer m = %d\n\n",mAGE);
    m++;
#endif   // AGE

#ifdef OXYGEN
    for (i=0;i<=NXMEM-1;i++) {
	for (j=0;j<=NYMEM-1;j++) {
	    if (D[i][j]>MINIMUM_DEPTH) {
		for (k=0;k<=NZ-1;k++) {
		    tr[m][k][i][j] = oxy_init[k][i][j];
		}
	    } else {
		for (k=0;k<=NZ-1;k++)
		    tr[m][k][i][j] = misval;
	    }
	}                  
    }
# ifdef SPONGE
    set_up_sponge_field(sp_oxygen,oxygen,NZ);
# endif
    mOXYGEN=m;
    printf("index of OXYGEN tracer m = %d\n\n",mOXYGEN);
    m++;
#endif   // OXYGEN

#ifdef OXY18
    for (i=0;i<=NXMEM-1;i++) {
	for (j=0;j<=NYMEM-1;j++) {
	    for (k=0;k<=NZ-1;k++) {
		if (D[i][j]>MINIMUM_DEPTH) {
		    tr[m][k][i][j] = o18_init[k][i][j];
		    o18[k][i][j] = o18_init[k][i][j];
		} else {
		    tr[m][k][i][j] = misval;
		    o18[k][i][j] = misval;
		}
	    }                            
	}
    }  
# ifdef SPONGE
    set_up_sponge_field(sp_o18,o18,NZ);
# endif
    mo18=m;
    printf("index of O18 tracer m = %d\n\n",mo18);
    m++;
#endif   // OXY18

#ifdef CFC
    for (i=0;i<=NXMEM-1;i++) {
	for (j=0;j<=NYMEM-1;j++) {
	    for (k=0;k<=NZ-1;k++) {
		if (D[i][j]>MINIMUM_DEPTH) {
		    tr[m][k][i][j]   = cfc11_init[k][i][j];
		    cfc11[k][i][j]   = cfc11_init[k][i][j];
		    tr[m+1][k][i][j] = cfc12_init[k][i][j];
		    cfc12[k][i][j]   = cfc12_init[k][i][j];
		} else {
		    tr[m][k][i][j]   = misval;
		    cfc11[k][i][j]   = misval;
		    tr[m+1][k][i][j] = misval;
		    cfc12[k][i][j]   = misval;
		}
	    }                            
	}
    }
# ifdef SPONGE
    set_up_sponge_field(sp_age,cfc11,NZ);
    set_up_sponge_field(sp_age,cfc12,NZ);
# endif

    mCFC11=m;
    mCFC12=m+1;

    printf("index of CFC11 tracer m = %d\n",mCFC11);
    printf("index of CFC12 tracer m = %d\n\n",mCFC12);
    m++;m++;
#ifdef SAT_TRACER
    mCFC11_sat=m;
    mCFC12_sat=m+1;
    printf("index of CFC11 saturation pseudo-tracer m = %d\n",mCFC11_sat);
    printf("index of CFC12 saturation pseudo-tracer m = %d\n\n",mCFC12_sat);
    m++;m++;
#endif
#endif   // CFC

// begin ashao
#ifdef SF6
    for (i=0;i<=NXMEM-1;i++) {
	for (j=0;j<=NYMEM-1;j++) {
	    for (k=0;k<=NZ-1;k++) {
		if (D[i][j]>MINIMUM_DEPTH) {
		    tr[m][k][i][j]   = sf6_init[k][i][j];
		    sf6[k][i][j]   = sf6_init[k][i][j];
		} else {
		    tr[m][k][i][j]   = misval;
		    sf6[k][i][j]   = misval;
		}
	    }                            
	}
    }
# ifdef SPONGE
    set_up_sponge_field(sp_age,sf6,NZ);
# endif

    mSF6=m;
#ifdef SAT_TRACER
    mSF6_sat=m+1;
    m++;
#endif

    printf("Index of SF6 tracer m = %d\n",mSF6);
    m++;
#endif   // SF6
#ifdef TTD_BP
	mttd=m;
	m++;
	printf("Index of TTD m = %d\n",mttd);
	for (k=0;k<=1;k++)
	    for (i=0;i<NXMEM;i++)
		for (j=0;j<NYMEM;j++) {
		    if (ttd_mask[i][j]) {
		    ttd[k][i][j]=1.0/12;
		    tr[mttd][k][i][j]=1.0/12;
		    }
		    else {
			    ttd[k][i][j]=0.0;
			    tr[mttd][k][i][j]=0.0;
		    }
	}
#endif
// end ashao
#ifdef PHOSPHATE
    for (i=0;i<=NXMEM-1;i++) {
	for (j=0;j<=NYMEM-1;j++) {
	    if (D[i][j]>MINIMUM_DEPTH) {
		for (k=0;k<=NZ-1;k++) {
		    if (phosphate_init[k][i][j] < 0.0)
		      printf("po4_init[%d][%d][%d]=%g\n",k,i,j,phosphate_init[k][i][j]);
		    tr[m][k][i][j]     = phosphate_init[k][i][j];
		    if (dop_init[k][i][j] < 0.0)
		      printf("dop_init[%d][%d][%d]=%g\n",k,i,j,dop_init[k][i][j]);
		    tr[m+1][k][i][j]     = dop_init[k][i][j];
		}
	    } else {
		for (k=0;k<=NZ-1;k++) {
		  tr[m][k][i][j]   = 0; // HF misval;
		  tr[m+1][k][i][j] = 0; // HF misval;
		}
	    }
	}
    }
    free3d(dop_init, NZ);

# ifdef SPONGE
    set_up_sponge_field(sp_phosphate,phosphate,NZ);
    set_up_sponge_field(sp_dop,dop,NZ);
# endif

    mPHOSPHATE=m;
    mDOP=m+1;
    printf("index of PHOSPHATE tracer m = %d\n",mPHOSPHATE);
    printf("index of DOP tracer m = %d\n",mDOP);
    m++;m++;
#endif   // PHOSPHATE

#ifdef NITRATE
    for (i=0;i<=NXMEM-1;i++) {
	for (j=0;j<=NYMEM-1;j++) {
	    if (D[i][j]>MINIMUM_DEPTH) {
		for (k=0;k<=NZ-1;k++) {
		  if (nitrate_init[k][i][j] < 0.0)
		      printf("no3_init[%d][%d][%d]=%g\n",k,i,j,nitrate_init[k][i][j]);
		    tr[m][k][i][j]     = nitrate_init[k][i][j];
		    if (don_init[k][i][j] < 0.0)
		      printf("don_init[%d][%d][%d]=%g\n",k,i,j,don_init[k][i][j]);
		    tr[m+1][k][i][j]     = don_init[k][i][j];
		}
	    } else {
		for (k=0;k<=NZ-1;k++) {
		  tr[m][k][i][j]   = 0; // HF misval;
		  tr[m+1][k][i][j] = 0; // HF misval;
		}
	    }
	}
    }
    free3d(don_init, NZ);
# ifdef SPONGE
    set_up_sponge_field(sp_nitrate,nitrate,NZ);
    set_up_sponge_field(sp_don,don,NZ);
# endif

    mNITRATE=m;
    mDON=m+1;
    printf("index of NITRATE tracer m = %d\n",mNITRATE);
    printf("index of DON tracer m = %d\n",mDON);
    m++;m++;

# ifdef N15_CYCLE
    for (i=0;i<=NXMEM-1;i++) {
	for (j=0;j<=NYMEM-1;j++) {
	    if (D[i][j]>MINIMUM_DEPTH) {
		for (k=0;k<=NZ-1;k++) {
		  if (nitrate_init_15n[k][i][j] < 0.0)
		    printf("no3_init_15n[%d][%d][%d]=%g\n",k,i,j,nitrate_init_15n[k][i][j]);
		  tr[m][k][i][j]     = nitrate_init_15n[k][i][j];
		  if (don_init_15n[k][i][j] < 0.0)
		    printf("don_init_15n[%d][%d][%d]=%g\n",k,i,j,don_init_15n[k][i][j]);
		  tr[m+1][k][i][j]     = don_init_15n[k][i][j];
		}
	    } else {
		for (k=0;k<=NZ-1;k++) {
		  tr[m][k][i][j]   = 0; // HF misval;
		  tr[m+1][k][i][j] = 0; // HF misval;
		}
	    }
	}
    }
    free3d(don_init_15n, NZ);
# ifdef SPONGE
    set_up_sponge_field(sp_nitrate_15n,nitrate_15n,NZ);
    set_up_sponge_field(sp_don_15n,don_15n,NZ);
# endif

    mNITRATE_15n=m;
    mDON_15n=m+1;
    printf("index of NITRATE_15N tracer m = %d\n",mNITRATE_15n);
    printf("index of DON_15N tracer m = %d\n",mDON_15n);
    m++;m++;
# endif   /* N15_CYCLE */
#endif   /* NITRATE */

#ifdef DIC
    for (i=0;i<=NXMEM-1;i++) {
	for (j=0;j<=NYMEM-1;j++) {
	    if (D[i][j]>MINIMUM_DEPTH) {
		for (k=0;k<=NZ-1;k++) {
		    tr[m][k][i][j]     = dic_init[k][i][j];
		    dic[k][i][j] = dic_init[k][i][j];
		    tr[m+1][k][i][j]     = alk_init[k][i][j];
		    alk[k][i][j] = alk_init[k][i][j];
		}
	    } else {
		for (k=0;k<=NZ-1;k++) {
		    tr[m][k][i][j]   = misval;
		    dic[k][i][j]     = misval;
		    tr[m+1][k][i][j] = misval;
		    alk[k][i][j]     = misval;
		}
	    }                            
	}
    }
    free3d(dic_init, NZ);
    free3d(alk_init, NZ);
  
# ifdef SPONGE
    set_up_sponge_field(sp_dic,dic,NZ);
    set_up_sponge_field(sp_alk,alk,NZ);
# endif

    mDIC=m;
    mALK=m+1;
    printf("index of DIC tracer m = %d\n",mDIC);
    printf("index of ALK tracer m = %d\n",mALK);
    m++; m++;
#endif   // DIC

#ifdef MERGED_ML
    //BX merge 1st-2nd and 3rd-4th layer for auxiliary fields
    for (i=0;i<=NXMEM-1;i++) {
	for (j=0;j<=NYMEM-1;j++) {
	    for (k=0;k<=2;k=k+2) {
# ifdef OXY18
	      o18[k][i][j]       = (o18[k][i][j]*h[k][i][j]+
				    o18[k+1][i][j]*h[k+1][i][j])/
		(h[k][i][j] + h[k+1][i][j]);
	    o18[k+1][i][j]     = o18[k][i][j];
# endif
# ifdef DIC
	    dic[k][i][j]       = (dic[k][i][j]*h[k][i][j]+
				  dic[k+1][i][j]*h[k+1][i][j])/
	      (h[k][i][j] + h[k+1][i][j]);
	    dic[k+1][i][j]     = dic[k][i][j];

	    alk[k][i][j]       = (alk[k][i][j]*h[k][i][j]+
				  alk[k+1][i][j]*h[k+1][i][j])/
	      (h[k][i][j] + h[k+1][i][j]);
	    alk[k+1][i][j]     = alk[k][i][j];
# endif
# ifdef AGE
	    age[k][i][j]       = (age[k][i][j]*h[k][i][j]+
				  age[k+1][i][j]*h[k+1][i][j])/
	      (h[k][i][j] + h[k+1][i][j]);
	    age[k+1][i][j]     = age[k][i][j];
# endif
# ifdef CFC
	    cfc11[k][i][j]       = (cfc11[k][i][j]*h[k][i][j]+
				    cfc11[k+1][i][j]*h[k+1][i][j])/
	      (h[k][i][j] + h[k+1][i][j]);
	    cfc11[k+1][i][j]     = cfc11[k][i][j];

	    cfc12[k][i][j]       = (cfc12[k][i][j]*h[k][i][j]+
				    cfc12[k+1][i][j]*h[k+1][i][j])/
	      (h[k][i][j] + h[k+1][i][j]);
	    cfc12[k+1][i][j]     = cfc12[k][i][j];
# endif
// begin ashao
# ifdef SF6
	    sf6[k][i][j]       = (sf6[k][i][j]*h[k][i][j]+
				    sf6[k+1][i][j]*h[k+1][i][j])/
	      (h[k][i][j] + h[k+1][i][j]);
	    sf6[k+1][i][j]     = sf6[k][i][j];

# endif
// end ashao

# ifdef INERT
	    inert1[k][i][j]       = (inert1[k][i][j]*h[k][i][j]+
				     inert1[k+1][i][j]*h[k+1][i][j])/
	      (h[k][i][j] + h[k+1][i][j]);
	    inert1[k+1][i][j]     = inert1[k][i][j];
	    inert2[k][i][j]       = (inert2[k][i][j]*h[k][i][j]+
				     inert2[k+1][i][j]*h[k+1][i][j])/
	      (h[k][i][j] + h[k+1][i][j]);
	    inert2[k+1][i][j]     = inert2[k][i][j];
# endif
	    }
	}
    }
# ifdef SMOOTH_REST
    //BX Needed here as fields have been altered
    merge_ml_tr();
# endif
#endif
    /* zonal, meridional re-entrance    */
    for (m=0;m<NTR;m++) {
	for (k=0;k<NZ;k++) {
	    for (j=0;j<=NYMEM-1;j++) {
		tr[m][k][0][j] = tr[m][k][nx-1][j];
		tr[m][k][1][j] = tr[m][k][nx][j];	
		tr[m][k][nx+1][j] = tr[m][k][2][j];
		tr[m][k][nx+2][j] = tr[m][k][3][j];
	    }
	}
	for (k=0;k<NZ;k++) {
	    for (i=2;i<=nx;i++) {
		tr[m][k][363-i][ny+1] = tr[m][k][i][ny];
		tr[m][k][363-i][ny+2] = tr[m][k][i][ny-1];
	    }
	}
    }
}
