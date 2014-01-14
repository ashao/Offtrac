#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include "netcdf.h"
#include "init.h"
#include "io.h"
#include "iocdf.h"
#include "alloc.h"
#include "read.h"
#include "step.h"
#include "util.h"
#include "metrics.h"

extern const double misval;

/* HF extern  double dxu[NYTOT][NXTOT];
extern  double dyu[NYTOT][NXTOT];

extern  double dxv[NYTOT][NXTOT];
extern  double dyv[NYTOT][NXTOT];

extern  double dxh[NYTOT][NXTOT];
extern  double dyh[NYTOT][NXTOT];

extern  double dxq[NYTOT][NXTOT];
extern  double dyq[NYTOT][NXTOT]; HF-e */

extern  double areagr[NXMEM][NYMEM];

extern double ***u;
extern double ***v;
extern double h[NZ][NXMEM][NYMEM];
extern double hend[NZ][NXMEM][NYMEM];
extern double hstart[NZ][NXMEM][NYMEM];
extern double rml[2][NXMEM][NYMEM];
extern double po4_star_lev[NZPHOS][NXMEM][NYMEM];
extern double po4_star_lay[NZ][NXMEM][NYMEM];
#ifdef PROGNOSTIC
extern double no3_lev[NZPHOS][NXMEM][NYMEM];
extern double no3_lay[NZ][NXMEM][NYMEM];
extern double fe_lev[NZPHOS][NXMEM][NYMEM];
extern double fe_lay[NZ][NXMEM][NYMEM];
extern double sfc_swr[NXMEM][NYMEM];
#else
# ifdef NITRATE
extern double no3_lev[NZPHOS][NXMEM][NYMEM];
extern double no3_lay[NZ][NXMEM][NYMEM];
# endif
#endif
extern double salt_woa[NXMEM][NYMEM];
#ifdef REST_ARCTIC
# ifdef DIC
extern double dic_lay[NZ][NXMEM][NYMEM];
extern double alk_lay[NZ][NXMEM][NYMEM];
# endif
# ifdef OXYGEN
extern double oxy_lay[NZ][NXMEM][NYMEM];
# endif
#endif
extern double dt;
/*   Begin added DT    */
extern int beginyear;
extern int theyear;
extern int lastsave;
/*   End added DT      */

extern double ***uhtm,***vhtm;
extern double ***ea,***eb;
extern double eaml[NXMEM][NYMEM];
extern double Temptm[NZ][NXMEM][NYMEM];
extern double Salttm[NZ][NXMEM][NYMEM];
extern double atmpres[NXMEM][NYMEM];
extern double xkw[NXMEM][NYMEM];
extern double fice[NXMEM][NYMEM];
extern double wd[NZ+1][NXMEM][NYMEM];

extern double D[NXMEM][NYMEM];
#ifdef RESTART
extern int rflags[NOVARS];
#endif
extern char directory[75];
extern char inittrac[200];
extern char kw_varname[100]; // ashao

extern double qlat[NYMEM],hlat[NYMEM];

#ifdef INERT
extern double inert1_init[NZ][NXMEM][NYMEM];
extern double inert2_init[NZ][NXMEM][NYMEM];
#endif

#ifdef AGE
extern double age_init[NZ][NXMEM][NYMEM];
#endif

#ifdef OXYGEN
extern double ***oxy_init;
#endif
#ifdef OXY18
extern double o18_init[NZ][NXMEM][NYMEM];
#endif

#ifdef CFC
extern double cfc11_init[NZ][NXMEM][NYMEM];
extern double cfc12_init[NZ][NXMEM][NYMEM];
#endif

// begin ashao
#ifdef SF6
extern double sf6_init[NZ][NXMEM][NYMEM];
#endif
// end ashao

#ifdef PHOSPHATE
extern double ***phosphate_init;
extern double ***dop_init;

#endif

#ifdef NITRATE
extern double ***nitrate_init;
extern double ***don_init;
# ifdef N15_CYCLE
extern double ***nitrate_init_15n;
extern double ***don_init_15n;
# endif
#endif /* NITRATE */

#ifdef DIC
extern double ***dic_init;
extern double ***alk_init;
#endif

#ifdef ENTRAIN
extern double ea_init[NZ][NXMEM][NYMEM];
extern double eaml_init[NXMEM][NYMEM];
#endif

#ifdef SMOOTH_REST
//extern double ***h_init;
extern double h_init[NZ][NXMEM][NYMEM];
#endif

#ifdef SPONGE
extern double Idamp[NXMEM][NYMEM];
extern double sp_e[NZ][NXMEM][NYMEM];
# if defined(AGE) || defined(CFC) || defined(SF6) //ashao
extern double sp_age[NZ][NXMEM][NYMEM];
# endif
# ifdef OXYGEN
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
extern double sp_alk[NZ][NXMEM][NYMEM];
# endif
#endif

#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(2);}

void read_fields(int imon,int itts)
    {

    /* Read ocmip gas exchange parameter fields */
    int i,j;
    int err, cdfid, timeid;
    char infile[25], inpath[200];
    FILE *file;
    int status;

    int xkwid, ficeid, atmpid;

    size_t start[MAX_NC_VARS];
    size_t count[MAX_NC_VARS];

    float** tmp2d;
    int inxt,iprv;
    double fact0,fact1,fact2;
    float** tmp2dp;
    float** tmp2dn;

    if (imon==11) {
	inxt = 0;
	iprv = 10;
    } else if (imon==0) {
	inxt = 1;
	iprv = 11;
    } else {
	inxt = imon + 1;
	iprv = imon - 1;
    }
    //  calculate weight factors for input files - fact0 for previous month,
    //   fact1 for imon, fact2 for next month
    fact1 = 0.5 + ((2.0*(double)itts-1.0) / (2.0*(double)NTSTEP));
    if (fact1 <= 1.0) {
	fact0 = 1.0 - fact1;
	fact2 = 0.0;
    } else {
	fact0 = 0.0;
	fact2 = fact1 - 1.0;
	fact1 = 2.0 - fact1;
    }

    // ashao: Modifying so that other things other than OCMIP fields can be used
    printf("Reading gasex fields: \n");
    sprintf(infile,"gasx_himgrid.nc");

    strcpy(inpath, directory);
    strcat(inpath, infile);
    printf("Reading gasex fields: %s\n",inpath);

    err = open_input_file(inpath,&file,&cdfid,&timeid);
    if (err != 0) {
	strcat(inpath, ".cdf");
	err = open_input_file(inpath,&file,&cdfid,&timeid);
	if (err != 0) {
	    printf("Unable to find gasex file.\n");
	    exit(-73);
	}
    }

    printf("read gasex month=%i\n",imon);

    if ((status = nc_inq_varid(cdfid, "OCMIP_FICE", &ficeid)))
	ERR(status);
    if ((status = nc_inq_varid(cdfid, kw_varname, &xkwid)))
	ERR(status);
    if ((status = nc_inq_varid(cdfid, "OCMIP_ATMP", &atmpid)))
	ERR(status);

    bzero(start, MAX_NC_VARS * sizeof(long));

    count[0] = 1;
    count[1] = NYTOT;
    count[2] = NXTOT;

    tmp2d = alloc2d_f(NYTOT, NXTOT);
    tmp2dp = alloc2d_f(NYTOT, NXTOT);
    tmp2dn = alloc2d_f(NYTOT, NXTOT);
    start[0] = imon;
    if ((status = nc_get_vara_float(cdfid,atmpid,start,count,tmp2d[0])))
	ERR(status);
    start[0] = iprv;
    if ((status = nc_get_vara_float(cdfid,atmpid,start,count,tmp2dp[0])))
	ERR(status);
    start[0] = inxt;
    if ((status = nc_get_vara_float(cdfid,atmpid,start,count,tmp2dn[0])))
	ERR(status);
    for (i=0;i<NXTOT;i++)
	for (j=0;j<NYTOT;j++)
	    atmpres[i+2][j+2]= fact1 * tmp2d[j][i] +
	    fact0 * tmp2dp[j][i] + fact2 * tmp2dn[j][i];

    start[0] = imon;
    if ((status = nc_get_vara_float(cdfid,ficeid,start,count,tmp2d[0])))
	ERR(status);
    start[0] = iprv;
    if ((status = nc_get_vara_float(cdfid,ficeid,start,count,tmp2dp[0])))
	ERR(status);
    start[0] = inxt;
    if ((status = nc_get_vara_float(cdfid,ficeid,start,count,tmp2dn[0])))
	ERR(status);
    for (i=0;i<NXTOT;i++)
	for (j=0;j<NYTOT;j++)
	    fice[i+2][j+2]= fact1 * tmp2d[j][i] +
	    fact0 * tmp2dp[j][i] + fact2 * tmp2dn[j][i];

    start[0] = imon;
    if ((status = nc_get_vara_float(cdfid,xkwid,start,count,tmp2d[0])))
	ERR(status);
    start[0] = iprv;
    if ((status = nc_get_vara_float(cdfid,xkwid,start,count,tmp2dp[0])))
	ERR(status);
    start[0] = inxt;
    if ((status = nc_get_vara_float(cdfid,xkwid,start,count,tmp2dn[0])))
	ERR(status);
    for (i=0;i<NXTOT;i++)
	for (j=0;j<NYTOT;j++)
	    xkw[i+2][j+2]= fact1 * tmp2d[j][i] +
	    fact0 * tmp2dp[j][i] + fact2 * tmp2dn[j][i];

    free2d_f(tmp2dp, NYTOT);
    free2d_f(tmp2dn, NYTOT);
    free2d_f(tmp2d, NYTOT);
    //BX  temporary bug fix for northern-most row (j=211)
    for (i=2;i<=NXMEM;i++) {
	atmpres[i][211] = atmpres[i][210];
	fice[i][211]    = fice[i][210];
	xkw[i][211]     = xkw[i][210];
    }
    close_file(&cdfid,&file);
    }

void read_fieldave(int imon)
    {
    }

void read_grid()
    {
    int err, i, j, cdfid, timeid;
    int status, varid; // HF
    char infile[25], inpath[200];
    FILE *file;
    double latq[NYTOT];
    double lath[NYTOT];
    double Ah[NYTOT][NXTOT];
    //HF
    double dxh_in[NYTOT][NXTOT], dxq_in[NYTOT][NXTOT];
    double dxu_in[NYTOT][NXTOT], dxv_in[NYTOT][NXTOT];
    double dyh_in[NYTOT][NXTOT], dyq_in[NYTOT][NXTOT];
    double dyu_in[NYTOT][NXTOT], dyv_in[NYTOT][NXTOT];

    long  start[MAX_NC_VARS];
    long  end[MAX_NC_VARS];

    sprintf(infile,"metrics.nc");

    strcpy(inpath, directory);
    strcat(inpath, infile);
    printf("\nLooking for file '%s'.\n",inpath);

    err = open_input_file(inpath,&file,&cdfid,&timeid);
    if (err != 0) {
	strcat(inpath, ".cdf");
	err = open_input_file(inpath,&file,&cdfid,&timeid);
	if (err != 0) {
	    printf("Unable to find grid file.\n");
	    exit(-73);
	}
    }

    end[0] = NYTOT;
    end[1] = NXTOT;

    //HF
    status = nc_inq_varid(cdfid, "latq", &varid);
    if (status != NC_NOERR) handle_error("inq latq", status);
    status = nc_get_vara_double(cdfid, varid, start, end, latq);
    if (status != NC_NOERR) handle_error("read latq", status);
    //BX-a  ncvarget(cdfid, 2, start, end, latq);

    //CAD should these loop indices be hardwired (to i=0;i<NYTOT) ?
    for (i=0;i<NYTOT;i++) {
	qlat[i+2]=(double)latq[i];
    }

    //HF ncvarget(cdfid, 2, start, end, lath);
    status = nc_inq_varid(cdfid, "lath", &varid);
    if (status != NC_NOERR) handle_error("inq lath", status);
    status = nc_get_vara_double(cdfid, varid, start, end, lath);
    if (status != NC_NOERR) handle_error("read lath", status);

    //CAD should these loop indices be hardwired (to i=0;i<NYTOT) ?
    for (i=0;i<NYTOT;i++) {
	hlat[i+2]=(double)lath[i];
    }
    //BX-e

    //   ncvarget(cdfid,10, start, end, dxv);
    status = nc_inq_varid(cdfid, "dxv", &varid);
    if (status != NC_NOERR) handle_error("inq dxv", status);
    status = nc_get_vara_double(cdfid, varid, start, end, &dxv_in[0][0]);
    if (status != NC_NOERR) handle_error("read dxv", status);

    //   ncvarget(cdfid,12, start, end, dyu);
    status = nc_inq_varid(cdfid, "dyu", &varid);
    if (status != NC_NOERR) handle_error("inq dyu", status);
    status = nc_get_vara_double(cdfid, varid, start, end, &dyu_in[0][0]);
    if (status != NC_NOERR) handle_error("read dyu", status);

    //HF     ncvarget(cdfid,13, start, end, Ah);
    status = nc_inq_varid(cdfid, "Ah", &varid);
    if (status != NC_NOERR) handle_error("inq Ah", status);
    status = nc_get_vara_double(cdfid, varid, start, end, &Ah[0][0]);
    if (status != NC_NOERR) handle_error("read Ah", status);

    //   ncvarget(cdfid,14, start, end, dyv);
    status = nc_inq_varid(cdfid, "dyv", &varid);
    if (status != NC_NOERR) handle_error("inq dyv", status);
    status = nc_get_vara_double(cdfid, varid, start, end, &dyv_in[0][0]);
    if (status != NC_NOERR) handle_error("read dyv", status);

    //   ncvarget(cdfid,15, start, end, dxh);
    status = nc_inq_varid(cdfid, "dxh", &varid);
    if (status != NC_NOERR) handle_error("inq dxh", status);
    status = nc_get_vara_double(cdfid, varid, start, end, &dxh_in[0][0]);
    if (status != NC_NOERR) handle_error("read dxh", status);

    //   ncvarget(cdfid,16, start, end, dyq);
    status = nc_inq_varid(cdfid, "dyq", &varid);
    if (status != NC_NOERR) handle_error("inq dyq", status);
    status = nc_get_vara_double(cdfid, varid, start, end, &dyq_in[0][0]);
    if (status != NC_NOERR) handle_error("read dyq", status);

    //   ncvarget(cdfid,17, start, end, dxq);
    status = nc_inq_varid(cdfid, "dxq", &varid);
    if (status != NC_NOERR) handle_error("inq dxq", status);
    status = nc_get_vara_double(cdfid, varid, start, end, &dxq_in[0][0]);
    if (status != NC_NOERR) handle_error("read dxq", status);

    //   ncvarget(cdfid,18, start, end, dyh);
    status = nc_inq_varid(cdfid, "dyh", &varid);
    if (status != NC_NOERR) handle_error("inq dyh", status);
    status = nc_get_vara_double(cdfid, varid, start, end, &dyh_in[0][0]);
    if (status != NC_NOERR) handle_error("read dyh", status);

    //   ncvarget(cdfid,19, start, end, dxu);
    status = nc_inq_varid(cdfid, "dxu", &varid);
    if (status != NC_NOERR) handle_error("inq dxu", status);
    status = nc_get_vara_double(cdfid, varid, start, end, &dxu_in[0][0]);
    if (status != NC_NOERR) handle_error("read dxu", status);


    for (i=0;i<NXTOT;i++) {
	for (j=0;j<NYTOT;j++) {
	    areagr[i+2][j+2]= Ah[j][i];
	    DXq(i+2,j+2) = dxq_in[j][i];
	    DYq(i+2,j+2) = dyq_in[j][i];
	    DXv(i+2,j+2) = dxv_in[j][i];
	    DYv(i+2,j+2) = dyv_in[j][i];
	    DXh(i+2,j+2) = dxh_in[j][i];
	    DYh(i+2,j+2) = dyh_in[j][i];
	    DXu(i+2,j+2) = dxu_in[j][i];
	    DYu(i+2,j+2) = dyu_in[j][i];
	}
	// HF
	DXv(i+2,1) = DXv(i+2,2);
	DYv(i+2,1) = DYv(i+2,2);
	// HF-e
    }

    close_file(&cdfid,&file);
    printf("\nFinished reading file '%s'.\n",inpath);


    /* zonal re-entrance            */

    for (j=0;j<=NYMEM-1;j++) {
	areagr[nx+1][j] = areagr[2][j];
	areagr[nx+2][j] = areagr[3][j];
	areagr[0][j] = areagr[nx-1][j];
	areagr[1][j] = areagr[nx][j];

	DXq(nx+1,j) = DXq(2,j);
	DXq(nx+2,j) = DXq(3,j);
	DXq(0,j) = DXq(nx-1,j);
	DXq(1,j) = DXq(nx,j);

	DYq(nx+1,j) = DYq(2,j);
	DYq(nx+2,j) = DYq(3,j);
	DYq(0,j) = DYq(nx-1,j);
	DYq(1,j) = DYq(nx,j);

	DXv(nx+1,j) = DXv(2,j);
	DXv(nx+2,j) = DXv(3,j);
	DXv(0,j) = DXv(nx-1,j);
	DXv(1,j) = DXv(nx,j);

	DYv(nx+1,j) = DYv(2,j);
	DYv(nx+2,j) = DYv(3,j);
	DYv(0,j) = DYv(nx-1,j);
	DYv(1,j) = DYv(nx,j);

	DXh(nx+1,j) = DXh(2,j);
	DXh(nx+2,j) = DXh(3,j);
	DXh(0,j) = DXh(nx-1,j);
	DXh(1,j) = DXh(nx,j);

	DYh(nx+1,j) = DYh(2,j);
	DYh(nx+2,j) = DYh(3,j);
	DYh(0,j) = DYh(nx-1,j);
	DYh(1,j) = DYh(nx,j);

	DXu(nx+1,j) = DXu(2,j);
	DXu(nx+2,j) = DXu(3,j);
	DXu(0,j) = DXu(nx-1,j);
	DXu(1,j) = DXu(nx,j);

	DYu(nx+1,j) = DYu(2,j);
	DYu(nx+2,j) = DYu(3,j);
	DYu(0,j) = DYu(nx-1,j);
	DYu(1,j) = DYu(nx,j);
    }

    /* meridional re-entrance            */

    for (i=2;i<=nx;i++) {
	int ii = 363 - i;
	areagr[ii][ny+1] = areagr[i][ny];
	areagr[ii][ny+2] = areagr[i][ny-1];

	DXq(ii,ny+1) = DXq(i,ny);
	DXq(ii,ny+2) = DXq(i,ny-1);

	DYq(ii,ny+1) = DYq(i,ny);
	DYq(ii,ny+2) = DYq(i,ny-1);

	DXv(ii,ny+1) = DXv(i,ny);
	DXv(ii,ny+2) = DXv(i,ny-1);

	DYv(ii,ny+1) = DYv(i,ny);
	DYv(ii,ny+2) = DYv(i,ny-1);

	DXh(ii,ny+1) = DXh(i,ny);
	DXh(ii,ny+2) = DXh(i,ny-1);

	DYh(ii,ny+1) = DYh(i,ny);
	DYh(ii,ny+2) = DYh(i,ny-1);

	DXu(ii,ny+1) = DXu(i,ny);
	DXu(ii,ny+2) = DXu(i,ny-1);

	DYu(ii,ny+1) = DYu(i,ny);
	DYu(ii,ny+2) = DYu(i,ny-1);
    }

    }

void read_D()
    {
    int i,ii,j;
    int err, cdfid, timeid;
    char infile[25], inpath[200];
    FILE *file;

    /*  Read the depth and other data in from the saved binary files.     */

    sprintf(infile,"topo.nc");

    strcpy(inpath, directory);
    strcat(inpath, infile);
    printf("Looking for file '%s'.\n",inpath);

    err = open_input_file(inpath,&file,&cdfid,&timeid);
    if (err != 0) {
	strcat(inpath, ".cdf");
	err = open_input_file(inpath,&file,&cdfid,&timeid);
	if (err != 0) {
	    printf("Unable to find depth file.\n");
	    exit(-73);
	}
    }

    //HF  WARN read_field(cdfid,file,"D", NXTOT,NYTOT,1, 0,0,0, 1,1,0, 1, D[0]);
    read_field(cdfid,file,"D", NXTOT,NYTOT,1, 0,0,0, X1,Y1,0, 1, D[0]);
    close_file(&cdfid,&file);


    //HF10022009	zonal re-entrance
    for (j=0;j<=NYMEM-1;j++) {
	D[nx+1][j] = D[2][j];
	D[nx+2][j] = D[3][j];
	D[0][j] =   D[nx-1][j];
	D[1][j] =   D[nx][j];
    }

    //      meridional re-entrance
    for (i=0;i<=nx+2;i++) {
	ii = 363 - i;
	D[ii][ny+1] = D[i][ny];
	D[ii][ny+2] = D[i][ny-1];
    }
    }


void read_clim(int imon, int inxt, int ilst)
    {

    int i,j,k;
    int ii;
    int err, cdfid, timeid;
    char infile[25], inpath[200];
    FILE *file;
    int status;

    int hclimid,uhclimid,vhclimid,wdclimid,pmeid;

    double PME[NYTOT][NXTOT];

    size_t start[MAX_NC_VARS];
    size_t count[MAX_NC_VARS];

    float*** tmp3d;
    float**  tmp2d;

    sprintf(infile,"oceanclim.nc");

    strcpy(inpath, directory);
    strcat(inpath, infile);

    err = open_input_file(inpath,&file,&cdfid,&timeid);
    if (err != 0) {
	strcat(inpath, ".cdf");
	err = open_input_file(inpath,&file,&cdfid,&timeid);
	if (err != 0) {
	    printf("Unable to find oceanclim file.\n");
	    exit(-73);
	}
    }

    printf("read_clim month=  %i\n",imon);

    if ((status = nc_inq_varid(cdfid, "h", &hclimid)))
	ERR(status);
    if ((status = nc_inq_varid(cdfid, "uh", &uhclimid)))
	ERR(status);
    if ((status = nc_inq_varid(cdfid, "vh", &vhclimid)))
	ERR(status);
    if ((status = nc_inq_varid(cdfid, "wd", &wdclimid)))
	ERR(status);
    if ((status = nc_inq_varid(cdfid, "pme", &pmeid)))
	ERR(status);

    bzero(start, MAX_NC_VARS * sizeof(long));

    count[0] = 1;
    count[1] = NZ;
    count[2] = NYTOT;
    count[3] = NXTOT;

    start[0] = ilst;

    printf("1-start0=%i \n",start[0]);

    tmp3d = alloc3d_f(NZ+1,NYTOT,NXTOT);
#ifndef USE_CALC_H
    if ((status = nc_get_vara_float(cdfid,hclimid,start,count,tmp3d[0][0])))
	ERR(status);
    for (k=0;k<NZ;k++)
	for (i=0;i<NXTOT;i++)
	    for (j=0;j<NYTOT;j++)
		h[k][i+2][j+2]= tmp3d[k][j][i];
#else
    printf("not reading hclim, h: %g\n", h[0][50][50]);
#endif

    start[0] = imon;

    printf("2-start0=%d\n",start[0]);

    if ((status = nc_get_vara_float(cdfid,uhclimid,start,count,tmp3d[0][0])))
	ERR(status);
    for (k=0;k<NZ;k++)
	for (i=0;i<NXTOT;i++)
	    for (j=0;j<NYTOT;j++)
		uhtm[k][i+2][j+2]= tmp3d[k][j][i];
    if ((status = nc_get_vara_float(cdfid,vhclimid,start,count,tmp3d[0][0])))
	ERR(status);
    for (k=0;k<NZ;k++)
	for (i=0;i<NXTOT;i++)
	    for (j=0;j<NYTOT;j++)
		vhtm[k][i+2][j+2]= tmp3d[k][j][i];

    printf("read_clim next month=%i\n",inxt);

    start[0] = inxt;
    if ((status = nc_get_vara_float(cdfid,hclimid,start,count,&tmp3d[0][0][0])))
	ERR(status);

    for (k=0;k<NZ;k++)
	for (i=0;i<NXTOT;i++)
	    for (j=0;j<NYTOT;j++)
		hend[k][i+2][j+2]= tmp3d[k][j][i];


    count[1] = NZ+1;
    start[0] = imon;
    if ((status = nc_get_vara_float(cdfid,wdclimid,start,count,tmp3d[0][0])))
	ERR(status);

    for (k=0;k<NZ+1;k++)
	for (i=0;i<NXTOT;i++)
	    for (j=0;j<NYTOT;j++)
		wd[k][i+2][j+2]= dt * tmp3d[k][j][i];

    tmp2d = alloc2d_f(NYTOT, NXTOT);

    count[0] = 1;
    start[0]= imon;
    count[1] = NYTOT;
    start[1]=0;
    count[2] = NXTOT;
    start[2]=0;

    if ((status = nc_get_vara_float(cdfid,pmeid,start,count,tmp2d[0])))
	ERR(status);

    for (i=0;i<NXTOT;i++)
	for (j=0;j<NYTOT;j++)
	    wd[0][i+2][j+2]= dt *  tmp2d[j][i];

    printf("dt=%g,imon=%i,inxt=%i\n",dt,imon,inxt);

    free3d_f(tmp3d, NZ);
    free2d_f(tmp2d,NYTOT);

    //	zonal re-entrance
    for (k=0;k<NZ;k++) {
	for (j=0;j<=NYMEM-1;j++) {
	    h[k][nx+1][j] = h[k][2][j];
	    h[k][nx+2][j] = h[k][3][j];
	    h[k][0][j] =   h[k][nx-1][j];
	    h[k][1][j] =   h[k][nx][j];
	    hend[k][nx+1][j] = hend[k][2][j];
	    hend[k][nx+2][j] = hend[k][3][j];
	    hend[k][0][j] =   hend[k][nx-1][j];
	    hend[k][1][j] =   hend[k][nx][j];
	}
    }

    //      meridional re-entrance
    for (i=2;i<=nx;i++) {
	ii = 363 - i;
	for (k=0;k<NZ;k++) {
	    h[k][ii][ny+1] = h[k][i][ny];
	    h[k][ii][ny+2]   = h[k][i][ny-1];
	    hend[k][ii][ny+1] = hend[k][i][ny];
	    hend[k][ii][ny+2] = hend[k][i][ny-1];
	}
    }

    close_file(&cdfid,&file);
    }

void read_uvw(int imon,int itts)
    {

    int i,j,k;
    int err, cdfid, timeid;
    char infile[25], inpath[200];
    FILE *file;
    int status;
    int inxt,iprv;
    int nmonths;
    double fact0,fact1,fact2;

    int uhclimid,vhclimid,wdclimid;

    size_t start[MAX_NC_VARS];
    size_t count[MAX_NC_VARS];

    float*** tmp3d;
    float*** tmp3dp;
    float*** tmp3dn;

    nmonths=NMONTHS;

    if (imon == nmonths) imon=0;
# ifdef SEPFILES
    if (imon==nmonths-1) {
	inxt = 0;
	iprv = imon - 1;
    } else if (imon==0) {
	inxt = 1;
	iprv = nmonths-1;
    } else {
	inxt = imon + 1;
	iprv = imon - 1;
    }
# else
    if (imon==11) {
	inxt = 0;
	iprv = 10;
    } else if (imon==0) {
	inxt = 1;
	iprv = 11;
    } else {
	inxt = imon + 1;
	iprv = imon - 1;
    }
# endif
    printf("READ_UVW: iprv=  %i, imon= %i, inxt= %i\n",iprv,imon,inxt);


    //  calculate weight factors for input files - fact0 for previous month,
    //   fact1 for imon, fact2 for next month
    fact1 = 0.5 + ((2.0*(double)itts-1.0) / (2.0*(double)NTSTEP));
    if (fact1 <= 1.0) {
	fact0 = 1.0 - fact1;
	fact2 = 0.0;
    } else {
	fact0 = 0.0;
	fact2 = fact1 - 1.0;
	fact1 = 2.0 - fact1;
    }
    printf("READ_UVW; facts are %g, %g, %g.\n",fact0,fact1,fact2);

# ifdef SEPFILES
    //
    //   read in separate files for U, V, and W
    //
    //    Start with uhtm
#  ifdef VARIAB_FORC
    sprintf(infile,"UH.nc");
#  else
    sprintf(infile,"UH-clim.nc");
#  endif
    strcpy(inpath, directory);
    strcat(inpath, infile);

    //    printf("Looking for file '%s'.\n",inpath);

    err = open_input_file(inpath,&file,&cdfid,&timeid);
    if (err != 0) {
	strcat(inpath, ".cdf");
	err = open_input_file(inpath,&file,&cdfid,&timeid);
	if (err != 0) {
	    printf("Unable to find UH file.\n");
	    exit(-73);
	}
    }

#  ifdef VARIAB_FORC
    if ((status = nc_inq_varid(cdfid, "uh", &uhclimid)))
	ERR(status);
#  else
    if ((status = nc_inq_varid(cdfid, "uh", &uhclimid)))
	ERR(status);
#  endif
    //    if ((status = nc_inq_varid(cdfid, "UHCLIM", &uhclimid)))
    //if ((status = nc_inq_varid(cdfid, "UH", &uhclimid)))
    //   ERR(status);

    bzero(start, MAX_NC_VARS * sizeof(long));

    count[0] = 1; 
    count[1] = NZ;
    count[2] = NYTOT; 
    count[3] = NXTOT;

#  ifdef VARIAB_FORC
    printf("read UH month=  %i\n",imon);
#  else
    printf("read_clim UH month=  %i\n",imon);
#  endif

    // allocate NZ+1 so that it can be used for WDCLIM as well
    tmp3d  = alloc3d_f(NZ+1,NYTOT,NXTOT);
    tmp3dn = alloc3d_f(NZ+1,NYTOT,NXTOT);
    tmp3dp = alloc3d_f(NZ+1,NYTOT,NXTOT);

    start[0] = imon;
    if ((status = nc_get_vara_float(cdfid,uhclimid,start,count,tmp3d[0][0])))
	ERR(status);
    //printf("read_clim UH month=  %i done\n",imon);
    //    start[0] = iprv;
    //    if ((status = nc_get_vara_float(cdfid,uhclimid,start,count,tmp3dp[0][0])))
    //       ERR(status);
    //printf("read_clim UH month=  %i done\n",iprv);
    //    start[0] = inxt;
    //    if ((status = nc_get_vara_float(cdfid,uhclimid,start,count,tmp3dn[0][0])))
    //       ERR(status);
    //printf("read_clim UH month=  %i done\n",inxt);

    for (k=0;k<NZ;k++)
	for (i=0;i<NXTOT;i++)
	    for (j=0;j<NYTOT;j++)
		uhtm[k][i+2][j+2]= tmp3d[k][j][i];
    //uhtm[k][i+2][j+2]= fact1 * tmp3d[k][j][i] +
    //    fact0 * tmp3dp[k][j][i] + fact2 * tmp3dn[k][j][i];

    close_file(&cdfid,&file);
    //printf("READ uhtm=%g;%g,%g,%g\n",uhtm[0][190][26],tmp3d[0][24][188],
    //	   tmp3dp[0][24][188],tmp3dn[0][24][188]);

    //    Next vhtm
#  ifdef VARIAB_FORC
    sprintf(infile,"VH.nc");
#  else
    sprintf(infile,"VH-clim.nc");
#  endif

    strcpy(inpath, directory);
    strcat(inpath, infile);

    //  printf("Looking for file '%s'.\n",inpath);

    err = open_input_file(inpath,&file,&cdfid,&timeid);
    if (err != 0) {
	strcat(inpath, ".cdf");
	err = open_input_file(inpath,&file,&cdfid,&timeid);
	if (err != 0) {
	    printf("Unable to find VH file.\n");
	    exit(-73);
	}
    }

#  ifdef VARIAB_FORC
    if ((status = nc_inq_varid(cdfid, "vh", &vhclimid)))
	ERR(status);
#  else
    if ((status = nc_inq_varid(cdfid, "vh", &vhclimid)))
	ERR(status);
#  endif
    //    if ((status = nc_inq_varid(cdfid, "VHCLIM", &vhclimid)))
    //if ((status = nc_inq_varid(cdfid, "VH", &vhclimid)))
    //   ERR(status);

    bzero(start, MAX_NC_VARS * sizeof(long));

    count[0] = 1; 
    count[1] = NZ;
    count[2] = NYTOT; 
    count[3] = NXTOT;

#  ifdef VARIAB_FORC
    printf("read VH month=  %i\n",imon);
#  else
    printf("read_clim VH month=  %i\n",imon);
#  endif

    start[0] = imon;
    if ((status = nc_get_vara_float(cdfid,vhclimid,start,count,tmp3d[0][0])))
	ERR(status);
    //    start[0] = iprv;
    //    if ((status = nc_get_vara_float(cdfid,vhclimid,start,count,tmp3dp[0][0])))
    //       ERR(status);
    //    start[0] = inxt;
    //    if ((status = nc_get_vara_float(cdfid,vhclimid,start,count,tmp3dn[0][0])))
    //         ERR(status);

    for (k=0;k<NZ;k++)
	for (i=0;i<NXTOT;i++)
	    for (j=0;j<NYTOT;j++)
		vhtm[k][i+2][j+2]= tmp3d[k][j][i];
    //vhtm[k][i+2][j+2]= fact1 * tmp3d[k][j][i] +
    //    fact0 * tmp3dp[k][j][i] + fact2 * tmp3dn[k][j][i];

    close_file(&cdfid,&file);

    // finally wd
#  ifdef VARIAB_FORC
    sprintf(infile,"WD.nc");
#  else
    sprintf(infile,"WD-clim.nc");
#  endif

    strcpy(inpath, directory);
    strcat(inpath, infile);

    //  printf("Looking for file '%s'.\n",inpath);

    err = open_input_file(inpath,&file,&cdfid,&timeid);
    if (err != 0) {
	strcat(inpath, ".cdf");
	err = open_input_file(inpath,&file,&cdfid,&timeid);
	if (err != 0) {
	    printf("Unable to find WD file.\n");
	    exit(-73);
	}
    }

#  ifdef VARIAB_FORC
    if ((status = nc_inq_varid(cdfid, "wd", &wdclimid)))
	ERR(status);
#  else
    if ((status = nc_inq_varid(cdfid, "wd", &wdclimid)))
	ERR(status);
#  endif
    //    if ((status = nc_inq_varid(cdfid, "WDCLIM", &wdclimid)))
    //if ((status = nc_inq_varid(cdfid, "WD", &wdclimid)))
    //   ERR(status);

    bzero(start, MAX_NC_VARS * sizeof(long));

    count[0] = 1; 
    count[1] = NZ;
    count[2] = NYTOT; 
    count[3] = NXTOT;

#  ifdef VARIAB_FORC
    printf("read WD month=  %i\n",imon);
#  else
    printf("read_clim WD month=  %i\n",imon);
#  endif

    count[1] = NZ+1;
    start[0] = imon;
    if ((status = nc_get_vara_float(cdfid,wdclimid,start,count,tmp3d[0][0])))
	ERR(status);
    //    start[0] = iprv;
    //    if ((status = nc_get_vara_float(cdfid,wdclimid,start,count,tmp3dp[0][0])))
    //       ERR(status);
    //    start[0] = inxt;
    //    if ((status = nc_get_vara_float(cdfid,wdclimid,start,count,tmp3dn[0][0])))
    //       ERR(status);

    for (k=0;k<NZ+1;k++)
	for (i=0;i<NXTOT;i++)
	    for (j=0;j<NYTOT;j++)
		wd[k][i+2][j+2]= dt * tmp3d[k][j][i];
    //wd[k][i+2][j+2]= dt * (fact1 * tmp3d[k][j][i] +
    //    fact0 * tmp3dp[k][j][i] + fact2 * tmp3dn[k][j][i]);

    //BX-a
    // take out PME which is read in onto wd[0][i][j]
    //for (i=0;i<NXTOT;i++)
    //	for (j=0;j<NYTOT;j++)
    //	    wd[k][i+2][j+2]= 0;
    //    printf("ATTENTION: EmP is set to 0. for this run.\n");
    //BX-e


    free3d_f(tmp3d, NZ);
    free3d_f(tmp3dp, NZ);
    free3d_f(tmp3dn, NZ);


    close_file(&cdfid,&file);
# else     //  SEPFILES
    sprintf(infile,"oceanclim.nc");

    strcpy(inpath, directory);
    strcat(inpath, infile);

    //  printf("Looking for file '%s'.\n",inpath);

    err = open_input_file(inpath,&file,&cdfid,&timeid);
    if (err != 0) {
	strcat(inpath, ".cdf");
	err = open_input_file(inpath,&file,&cdfid,&timeid);
	if (err != 0) {
	    printf("Unable to find oceanclim file.\n");
	    exit(-73);
	}
    }

    if ((status = nc_inq_varid(cdfid, "uh", &uhclimid)))
	ERR(status);
    if ((status = nc_inq_varid(cdfid, "vh", &vhclimid)))
	ERR(status);
    if ((status = nc_inq_varid(cdfid, "wd", &wdclimid)))
	ERR(status);

    bzero(start, MAX_NC_VARS * sizeof(long));

    count[0] = 1; 
    count[1] = NZ;
    count[2] = NYTOT; 
    count[3] = NXTOT;

    printf("read_clim month=  %i\n",imon);

    // allocate NZ+1 so that it can be used for WDCLIM as well
    tmp3d  = alloc3d_f(NZ+1,NYTOT,NXTOT);
    tmp3dn = alloc3d_f(NZ+1,NYTOT,NXTOT);
    tmp3dp = alloc3d_f(NZ+1,NYTOT,NXTOT);

    //    Start with uhtm
    start[0] = imon;
    if ((status = nc_get_vara_float(cdfid,uhclimid,start,count,tmp3d[0][0])))
	ERR(status);
    start[0] = iprv;
    if ((status = nc_get_vara_float(cdfid,uhclimid,start,count,tmp3dp[0][0])))
	ERR(status);
    start[0] = inxt;
    if ((status = nc_get_vara_float(cdfid,uhclimid,start,count,tmp3dn[0][0])))
	ERR(status);

    for (k=0;k<NZ;k++)
	for (i=0;i<NXTOT;i++)
	    for (j=0;j<NYTOT;j++)
		uhtm[k][i+2][j+2]= fact1 * tmp3d[k][j][i] +
		fact0 * tmp3dp[k][j][i] + fact2 * tmp3dn[k][j][i];

    //    Next vhtm
    start[0] = imon;
    if ((status = nc_get_vara_float(cdfid,vhclimid,start,count,tmp3d[0][0])))
	ERR(status);
    start[0] = iprv;
    if ((status = nc_get_vara_float(cdfid,vhclimid,start,count,tmp3dp[0][0])))
	ERR(status);
    start[0] = inxt;
    if ((status = nc_get_vara_float(cdfid,vhclimid,start,count,tmp3dn[0][0])))
	ERR(status);

    for (k=0;k<NZ;k++)
	for (i=0;i<NXTOT;i++)
	    for (j=0;j<NYTOT;j++)
		vhtm[k][i+2][j+2]= fact1 * tmp3d[k][j][i] +
		fact0 * tmp3dp[k][j][i] + fact2 * tmp3dn[k][j][i];

    // finally wd
    count[1] = NZ+1;
    start[0] = imon;
    if ((status = nc_get_vara_float(cdfid,wdclimid,start,count,tmp3d[0][0])))
	ERR(status);
    start[0] = iprv;
    if ((status = nc_get_vara_float(cdfid,wdclimid,start,count,tmp3dp[0][0])))
	ERR(status);
    start[0] = inxt;
    if ((status = nc_get_vara_float(cdfid,wdclimid,start,count,tmp3dn[0][0])))
	ERR(status);

    for (k=0;k<NZ+1;k++)
	for (i=0;i<NXTOT;i++)
	    for (j=0;j<NYTOT;j++)
		wd[k][i+2][j+2]= dt * (fact1 * tmp3d[k][j][i] +
			fact0 * tmp3dp[k][j][i] + fact2 * tmp3dn[k][j][i]);

    free3d_f(tmp3d, NZ);
    free3d_f(tmp3dp, NZ);
    free3d_f(tmp3dn, NZ);

    close_file(&cdfid,&file);
# endif    //  SEPFILES
    }

void read_h(int imon,int inxt) 
    {

    int i,j,k;
    int ii;
    int err, cdfid, timeid;
    char infile[25], inpath[200];
    FILE *file;
    int status;

    int hclimid;

    size_t start[MAX_NC_VARS];
    size_t count[MAX_NC_VARS];

    float*** tmp3d;

# ifdef SEPFILES
#  ifdef VARIAB_FORC
    sprintf(infile,"H.nc");
#  else
    sprintf(infile,"H-clim.nc");
#  endif

    strcpy(inpath, directory);
    strcat(inpath, infile);

    //  printf("Looking for file '%s'.\n",inpath);

    err = open_input_file(inpath,&file,&cdfid,&timeid);
    if (err != 0) {
	strcat(inpath, ".cdf");
	err = open_input_file(inpath,&file,&cdfid,&timeid);
	if (err != 0) {
	    printf("Unable to find H oceanclim file.\n");
	    exit(-73);
	}
    }

#  ifdef VARIAB_FORC
    printf("read_h H month=  %i\n",imon);
#  else
    printf("read_hclim H month=  %i\n",imon);
#  endif

#  ifdef VARIAB_FORC
    if ((status = nc_inq_varid(cdfid, "h", &hclimid)))
	ERR(status);
#  else
    if ((status = nc_inq_varid(cdfid, "h", &hclimid)))
	ERR(status);
#  endif
    //   if ((status = nc_inq_varid(cdfid, "HCLIM", &hclimid)))
    //if ((status = nc_inq_varid(cdfid, "H", &hclimid)))
    //   ERR(status);

    bzero(start, MAX_NC_VARS * sizeof(long));

    count[0] = 1;
    count[1] = NZ;
    count[2] = NYTOT;
    count[3] = NXTOT;

    start[0] = imon;
    // allocate NZ+1 so that it can be used for WDCLIM as well
    tmp3d = alloc3d_f(NZ+1,NYTOT,NXTOT);
    if ((status = nc_get_vara_float(cdfid,hclimid,start,count,tmp3d[0][0])))
	ERR(status);
    for (k=0;k<NZ;k++) {
	for (i=0;i<NXTOT;i++) {
	    for (j=0;j<NYTOT;j++) {
		h[k][i+2][j+2]= tmp3d[k][j][i];
		hstart[k][i+2][j+2]= tmp3d[k][j][i];
	    }
	}
    }

#  ifdef VARIAB_FORC
    printf("read_h H next month=%i\n",inxt);
#  else
    printf("read_hclim H next month=%i\n",inxt);
#  endif

    start[0] = inxt;
    if ((status = nc_get_vara_float(cdfid,hclimid,start,count,&tmp3d[0][0][0])))
	ERR(status);

    for (k=0;k<NZ;k++)
	for (i=0;i<NXTOT;i++)
	    for (j=0;j<NYTOT;j++)
		hend[k][i+2][j+2]= tmp3d[k][j][i];

    printf("Read in h=%g, hstart=%g, hend=%g\n",h[1][100][100],hstart[1][100][100],hend[1][100][100]);

# else   //   SEPFILES
    sprintf(infile,"oceanclim.nc");

    strcpy(inpath, directory);
    strcat(inpath, infile);

    //  printf("Looking for file '%s'.\n",inpath);

    err = open_input_file(inpath,&file,&cdfid,&timeid);
    if (err != 0) {
	strcat(inpath, ".cdf");
	err = open_input_file(inpath,&file,&cdfid,&timeid);
	if (err != 0) {
	    printf("Unable to find oceanclim file.\n");
	    exit(-73);
	}
    }

    printf("read_hclim month=  %i\n",imon);

    if ((status = nc_inq_varid(cdfid, "HCLIM", &hclimid)))
	ERR(status);

    bzero(start, MAX_NC_VARS * sizeof(long));

    count[0] = 1;
    count[1] = NZ;
    count[2] = NYTOT;
    count[3] = NXTOT;

    start[0] = imon;
    // allocate NZ+1 so that it can be used for WDCLIM as well
    tmp3d = alloc3d_f(NZ+1,NYTOT,NXTOT);
    if ((status = nc_get_vara_float(cdfid,hclimid,start,count,tmp3d[0][0])))
	ERR(status);
    for (k=0;k<NZ;k++) {
	for (i=0;i<NXTOT;i++) {
	    for (j=0;j<NYTOT;j++) {
		h[k][i+2][j+2]= tmp3d[k][j][i];
		hstart[k][i+2][j+2]= tmp3d[k][j][i];
	    }
	}
    }

    printf("read_hclim next month=%i\n",inxt);

    start[0] = inxt;
    if ((status = nc_get_vara_float(cdfid,hclimid,start,count,&tmp3d[0][0][0])))
	ERR(status);

    for (k=0;k<NZ;k++)
	for (i=0;i<NXTOT;i++)
	    for (j=0;j<NYTOT;j++)
		hend[k][i+2][j+2]= tmp3d[k][j][i];

# endif   //   SEPFILES
    free3d_f(tmp3d, NZ);


    //	zonal re-entrance
    for (k=0;k<NZ;k++) {
	for (j=0;j<=NYMEM-1;j++) {
	    h[k][nx+1][j] = h[k][2][j];
	    h[k][nx+2][j] = h[k][3][j];
	    h[k][0][j] =   h[k][nx-1][j];
	    h[k][1][j] =   h[k][nx][j];
	    hstart[k][nx+1][j] = hstart[k][2][j];
	    hstart[k][nx+2][j] = hstart[k][3][j];
	    hstart[k][0][j] =   hstart[k][nx-1][j];
	    hstart[k][1][j] =   hstart[k][nx][j];
	    hend[k][nx+1][j] = hend[k][2][j];
	    hend[k][nx+2][j] = hend[k][3][j];
	    hend[k][0][j] =   hend[k][nx-1][j];
	    hend[k][1][j] =   hend[k][nx][j];
	}
    }

    //      meridional re-entrance
    for (i=2;i<=nx;i++) {
	ii = 363 - i;
	for (k=0;k<NZ;k++) {
	    h[k][ii][ny+1] = h[k][i][ny];
	    h[k][ii][ny+2]   = h[k][i][ny-1];
	    hstart[k][ii][ny+1] = hstart[k][i][ny];
	    hstart[k][ii][ny+2] = hstart[k][i][ny-1];
	    hend[k][ii][ny+1] = hend[k][i][ny];
	    hend[k][ii][ny+2] = hend[k][i][ny-1];
	}
    }

    close_file(&cdfid,&file);
    }

void read_biotic_bc(int imon,int itts)
    {
#ifdef PHOSPHATE
    /* Read po4_star_lev from climatology. */
    //BX  For changes in the po4 climatology field the number of depth levels
    //BX   in the original field (NZPHOS in init.h) and the values for the depth
    //BX   levels in step.c (two locations) must be changed.

    int i,j,k,ii;
    int err, cdfid, timeid;
    char infile[25], inpath[200];
    FILE *file;
    int status;

    int woapo4id;
# if defined PROGNOSTIC || defined NITRATE
    int woano3id, feid, lightid;
# endif

    size_t start[MAX_NC_VARS];
    size_t count[MAX_NC_VARS];

    float*** tmp3d;
    int inxt,iprv;
    double fact0,fact1,fact2;
    float*** tmp3dp;
    float*** tmp3dn;

    if (imon==11) {
	inxt = 0;
	iprv = 10;
    } else if (imon==0) {
	inxt = 1;
	iprv = 11;
    } else {
	inxt = imon + 1;
	iprv = imon - 1;
    }
    //  calculate weight factors for input files - fact0 for previous month,
    //   fact1 for imon, fact2 for next month
    fact1 = 0.5 + ((2.0*(double)itts-1.0) / (2.0*(double)NTSTEP));
    if (fact1 <= 1.0) {
	fact0 = 1.0 - fact1;
	fact2 = 0.0;
    } else {
	fact0 = 0.0;
	fact2 = fact1 - 1.0;
	fact1 = 2.0 - fact1;
    }

    printf("Reading observed PO4 fields: \n");

    sprintf(infile,"woalevpo4.nc"); // ashao: changing to default woa09 file name
    // sprintf(infile,"p_an.himgrid.monthly.nc");
    strcpy(inpath, directory);
    strcat(inpath, infile);

    printf("Looking for file '%s'.\n",inpath);

    err = open_input_file(inpath,&file,&cdfid,&timeid);
    if (err != 0) {
	strcat(inpath, ".cdf");
	err = open_input_file(inpath,&file,&cdfid,&timeid);
	if (err != 0) {
	    printf("Unable to find world ocean atlas PO4 file.\n");
	    exit(-73);
	}
    }

    printf("read po4 month=%i\n",imon);

    // if ((status = nc_inq_varid(cdfid, "p_an", &woapo4id)))
    if ((status = nc_inq_varid(cdfid, "WOAPO4", &woapo4id)))
	ERR(status);

    bzero(start, MAX_NC_VARS * sizeof(long));

    count[0] = 1;
    count[1] = NZPHOS;
    count[2] = NYTOT;
    count[3] = NXTOT;

    tmp3d = alloc3d_f(NZPHOS,NYTOT,NXTOT);
    start[0] = imon;
    if ((status = nc_get_vara_float(cdfid,woapo4id,start,count,tmp3d[0][0])))
	ERR(status);
    tmp3dp = alloc3d_f(NZPHOS,NYTOT,NXTOT);
    tmp3dn = alloc3d_f(NZPHOS,NYTOT,NXTOT);
    start[0] = iprv;
    if ((status = nc_get_vara_float(cdfid,woapo4id,start,count,tmp3dp[0][0])))
	ERR(status);
    start[0] = inxt;
    if ((status = nc_get_vara_float(cdfid,woapo4id,start,count,tmp3dn[0][0])))
	ERR(status);
    close_file(&cdfid,&file);

    for (k=0;k<NZPHOS;k++)
	for (i=0;i<NXTOT;i++)
	    for (j=0;j<NYTOT;j++)
		po4_star_lev[k][i+2][j+2]= 1e-3 * (fact1 * tmp3d[k][j][i] +
			fact0 * tmp3dp[k][j][i] + fact2 * tmp3dn[k][j][i]);
    //	zonal re-entrance
    // for (k=0;k<NZ;k++) {
    for (k=0;k<NZPHOS;k++) {
	for (j=0;j<=NYMEM-1;j++) {
	    po4_star_lev[k][nx+1][j] = po4_star_lev[k][2][j];
	    po4_star_lev[k][nx+2][j] = po4_star_lev[k][3][j];
	    po4_star_lev[k][0][j] =   po4_star_lev[k][nx-1][j];
	    po4_star_lev[k][1][j] =   po4_star_lev[k][nx][j];
	}
    }
    //      meridional re-entrance
    for (i=2;i<=nx;i++) {
	ii = 363 - i;
	//      for (k=0;k<NZ;k++) {
	for (k=0;k<NZPHOS;k++) {
	    po4_star_lev[k][ii][ny+1] = po4_star_lev[k][i][ny];
	    po4_star_lev[k][ii][ny+2] = po4_star_lev[k][i][ny-1];
	}
    }
# ifdef NITRATE
    printf("Reading observed NO3 fields: \n");
    sprintf(infile,"woalevno3.nc");
    strcpy(inpath, directory);
    strcat(inpath, infile);

    printf("Looking for file '%s'.\n",inpath);

    err = open_input_file(inpath,&file,&cdfid,&timeid);
    if (err != 0) {
	strcat(inpath, ".cdf");
	err = open_input_file(inpath,&file,&cdfid,&timeid);
	if (err != 0) {
	    printf("Unable to find world ocean atlas NO3 file.\n");
	    exit(-73);
	}
    }

    printf("read no3 month=%i\n",imon);

    if ((status = nc_inq_varid(cdfid, "nitrate", &woano3id)))
	ERR(status);

    bzero(start, MAX_NC_VARS * sizeof(long));

    count[0] = 1;
    count[1] = NZPHOS;
    count[2] = NYTOT;
    count[3] = NXTOT;

    //WARN  tmp3d = alloc3d_f(NZPHOS,NYTOT,NXTOT);
    start[0] = imon;
    if ((status = nc_get_vara_float(cdfid,woano3id,start,count,tmp3d[0][0])))
	ERR(status);
    //WARN tmp3dp = alloc3d_f(NZPHOS,NYTOT,NXTOT);
    //WARN tmp3dn = alloc3d_f(NZPHOS,NYTOT,NXTOT);
    start[0] = iprv;
    if ((status = nc_get_vara_float(cdfid,woano3id,start,count,tmp3dp[0][0])))
	ERR(status);
    start[0] = inxt;
    if ((status = nc_get_vara_float(cdfid,woano3id,start,count,tmp3dn[0][0])))
	ERR(status);
    close_file(&cdfid,&file);

    for (k=0;k<NZPHOS;k++)
	for (i=0;i<NXTOT;i++)
	    for (j=0;j<NYTOT;j++)
		no3_lev[k][i+2][j+2]= 1e-3 * (fact1 * tmp3d[k][j][i] +
			fact0 * tmp3dp[k][j][i] + fact2 * tmp3dn[k][j][i]);
    //	zonal re-entrance
    // for (k=0;k<NZ;k++) {
    for (k=0;k<NZPHOS;k++) {
	for (j=0;j<=NYMEM-1;j++) {
	    no3_lev[k][nx+1][j] = no3_lev[k][2][j];
	    no3_lev[k][nx+2][j] = no3_lev[k][3][j];
	    no3_lev[k][0][j] =   no3_lev[k][nx-1][j];
	    no3_lev[k][1][j] =   no3_lev[k][nx][j];
	}
    }
    //      meridional re-entrance
    for (i=2;i<=nx;i++) {
	ii = 363 - i;
	//      for (k=0;k<NZ;k++) {
	for (k=0;k<NZPHOS;k++) {
	    no3_lev[k][ii][ny+1] = no3_lev[k][i][ny];
	    no3_lev[k][ii][ny+2] = no3_lev[k][i][ny-1];
	}
    }
# endif /* NITRATE */

# ifdef PROGNOSTIC
    printf("Reading modeled Fe fields: \n");

    sprintf(infile,"whoi_fe.nc");
    strcpy(inpath, directory);
    strcat(inpath, infile);

    printf("Looking for file '%s'.\n",inpath);

    err = open_input_file(inpath,&file,&cdfid,&timeid);
    if (err != 0) {
	strcat(inpath, ".cdf");
	err = open_input_file(inpath,&file,&cdfid,&timeid);
	if (err != 0) {
	    printf("Unable to find model Fe file.\n");
	    exit(-73);
	}
    }

    printf("read Fe month=%i\n",imon);

    if ((status = nc_inq_varid(cdfid, "Fe", &feid)))
	ERR(status);

    bzero(start, MAX_NC_VARS * sizeof(long));

    count[0] = 1;
    count[1] = NZPHOS;
    count[2] = NYTOT;
    count[3] = NXTOT;

    // tmp3d = alloc3d_f(NZPHOS,NYTOT,NXTOT);
    start[0] = imon;
    if ((status = nc_get_vara_float(cdfid,feid,start,count,tmp3d[0][0])))
	ERR(status);
    // tmp3dp = alloc3d_f(NZPHOS,NYTOT,NXTOT);
    // tmp3dn = alloc3d_f(NZPHOS,NYTOT,NXTOT);
    start[0] = iprv;
    if ((status = nc_get_vara_float(cdfid,feid,start,count,tmp3dp[0][0])))
	ERR(status);
    start[0] = inxt;
    if ((status = nc_get_vara_float(cdfid,feid,start,count,tmp3dn[0][0])))
	ERR(status);
    close_file(&cdfid,&file);

    for (k=0;k<NZPHOS;k++)
	for (i=0;i<NXTOT;i++)
	    for (j=0;j<NYTOT;j++)
		fe_lev[k][i+2][j+2]= 1e-3 * (fact1 * tmp3d[k][j][i] +
			fact0 * tmp3dp[k][j][i] + fact2 * tmp3dn[k][j][i]);
    //	zonal re-entrance
    // for (k=0;k<NZ;k++) {
    for (k=0;k<NZPHOS;k++) {
	for (j=0;j<=NYMEM-1;j++) {
	    fe_lev[k][nx+1][j] = fe_lev[k][2][j];
	    fe_lev[k][nx+2][j] = fe_lev[k][3][j];
	    fe_lev[k][0][j] =   fe_lev[k][nx-1][j];
	    fe_lev[k][1][j] =   fe_lev[k][nx][j];
	}
    }
    //      meridional re-entrance
    for (i=2;i<=nx;i++) {
	ii = 363 - i;
	//      for (k=0;k<NZ;k++) {
	for (k=0;k<NZPHOS;k++) {
	    fe_lev[k][ii][ny+1] = fe_lev[k][i][ny];
	    fe_lev[k][ii][ny+2] = fe_lev[k][i][ny-1];
	}
    }

    printf("Reading light fields: \n");

    sprintf(infile,"isccp_swrad.nc");
    strcpy(inpath, directory);
    strcat(inpath, infile);

    printf("Looking for file '%s'.\n",inpath);

    err = open_input_file(inpath,&file,&cdfid,&timeid);
    if (err != 0) {
	strcat(inpath, ".cdf");
	err = open_input_file(inpath,&file,&cdfid,&timeid);
	if (err != 0) {
	    printf("Unable to find model sw radiation file.\n");
	    exit(-73);
	}
    }

    printf("read light month=%i\n",imon);

    if ((status = nc_inq_varid(cdfid, "isccp_DWN_SRF_clim", &lightid)))
	ERR(status);
    printf("reading light month, status=%i\n",status);

    bzero(start, MAX_NC_VARS * sizeof(long));

    // the (surface) light field does not have a z dimension
    count[0] = 1;
    count[1] = NYTOT;
    count[2] = NXTOT;

    // tmp3d = alloc3d_f(NZPHOS,NYTOT,NXTOT);
    start[0] = imon;
    if ((status = nc_get_vara_float(cdfid,lightid,start,count,tmp3d[0][0])))
	ERR(status);

    // tmp3dp = alloc3d_f(NZPHOS,NYTOT,NXTOT);
    // tmp3dn = alloc3d_f(NZPHOS,NYTOT,NXTOT);
    start[0] = iprv;
    if ((status = nc_get_vara_float(cdfid,lightid,start,count,tmp3dp[0][0])))
	ERR(status);
    start[0] = inxt;
    if ((status = nc_get_vara_float(cdfid,lightid,start,count,tmp3dn[0][0])))
	ERR(status);
    close_file(&cdfid,&file);

    for (i=0;i<NXTOT;i++)
	for (j=0;j<NYTOT;j++)
	    sfc_swr[i+2][j+2]= fact1 * tmp3d[0][j][i] +
	    fact0 * tmp3dp[0][j][i] + fact2 * tmp3dn[0][j][i];
    //	zonal re-entrance
    // for (k=0;k<NZ;k++) {
    for (j=0;j<=NYMEM-1;j++) {
	sfc_swr[nx+1][j] = sfc_swr[2][j];
	sfc_swr[nx+2][j] = sfc_swr[3][j];
	sfc_swr[0][j] =   sfc_swr[nx-1][j];
	sfc_swr[1][j] =   sfc_swr[nx][j];
    }
    //      meridional re-entrance
    for (i=2;i<=nx;i++) {
	ii = 363 - i;
	sfc_swr[ii][ny+1] = sfc_swr[i][ny];
	sfc_swr[ii][ny+2] = sfc_swr[i][ny-1];
    }

# endif /* PROGNOSTIC */

    free3d_f(tmp3d, NZPHOS);
    free3d_f(tmp3dp, NZPHOS);
    free3d_f(tmp3dn, NZPHOS);
    close_file(&cdfid,&file);
    printf("end read biotic bc for month=%i\n",imon);
#endif /* PHOSPHATE */
    }

#ifdef LEV_OXY
void read_oxy_ic()
    {
    int i,j,k;
    int err, cdfid, timeid;
    char infile[25], inpath[200];
    FILE *file;
    int status;

    int levo2id;

    size_t start[MAX_NC_VARS];
    size_t count[MAX_NC_VARS];

    float*** tmp3d;
    double*** oxytmp;

    int nzlevitus = NZPHOS;
    double depth[NZ][NXMEM][NYMEM];
    double po4obsprof[NZPHOS];
    double levitus_depths[NZPHOS] = {0, 10, 20, 30, 50, 75, 100,
	    120, 150, 200, 250, 300, 400, 500, 600,
	    700, 800, 900, 1000, 1100, 1200, 1300,
	    1400, 1500, 1750, 2000, 2500, 3000,
	    3500, 4000, 4500, 5000, 5500};

    printf("Reading Levitus O2 climatology: \n");

    sprintf(infile,"lev94_o2.nc");
    strcpy(inpath, directory);
    strcat(inpath, infile);

    printf("Looking for file '%s'.\n",inpath);

    err = open_input_file(inpath,&file,&cdfid,&timeid);
    if (err != 0) {
	strcat(inpath, ".cdf");
	err = open_input_file(inpath,&file,&cdfid,&timeid);
	if (err != 0) {
	    printf("Unable to find Levitus O2 file.\n");
	    exit(-73);
	}
    }

    if ((status = nc_inq_varid(cdfid, "LEVO2", &levo2id)))
	ERR(status);

    bzero(start, MAX_NC_VARS * sizeof(long));

    count[0] = 1;
    count[1] = NZPHOS;
    count[2] = NYTOT;
    count[3] = NXTOT;

    tmp3d = alloc3d_f(NZPHOS,NYTOT,NXTOT);
    oxytmp = alloc3d(NZPHOS,NXMEM,NYMEM);

    start[0] = 0;

    if ((status = nc_get_vara_float(cdfid,levo2id,start,count,tmp3d[0][0])))
	ERR(status);

    for (k=0;k<NZPHOS;k++) {
	for (i=0;i<NXTOT;i++) {
	    for (j=0;j<NYTOT;j++) {
		oxytmp[k][i+2][j+2]= tmp3d[k][j][i]*1.e-3;
	    }
	}
    }
    //  temporary bug fix for northern-most row (j=211)
    for (i=2;i<NXMEM;i++) {
	for (k=0;k<NZPHOS;k++) {
	    oxytmp[k][i][211] = oxytmp[k][i][210];
	}
    }

    z_depth(h,depth);

    for (i=X1;i<=nx;i++) {
	for (j=Y1;j<=ny;j++) {
	    if (D[i][j]>MINIMUM_DEPTH) {
		for (k=0;k<nzlevitus;k++)
		    po4obsprof[k] = oxytmp[k][i][j];
		for (k=0;k<NZ;k++) {
		    oxy_init[k][i][j] = lin_interp(depth[k][i][j], po4obsprof,
			    levitus_depths, 0, nzlevitus);
		    if (oxy_init[k][i][j] < 0.e0) oxy_init[k][i][j] = 0.;
		}
	    } else {
		for (k=0;k<NZ;k++) {
		    oxy_init[k][i][j] = misval;
		}
	    }
	}
    }

    free3d_f(tmp3d, NZPHOS);
    free3d(oxytmp, NZPHOS);
    close_file(&cdfid,&file);
    }
#endif

#ifdef GLODAP_IC
void read_glodap_ic()
    {
    int i,j,k;
    int err, cdfid, timeid;
    char infile[25], inpath[200];
    FILE *file;
    int status;

    int glodapid;

    size_t start[MAX_NC_VARS];
    size_t count[MAX_NC_VARS];

    float*** tmp3d;
    double*** dictmp;

    //BX-a need this part for extrapolation on layers
    //BX (po4_obs_layer does not work due to field definition
    //BX of oxytmp and oxy_init)
    int nzlevitus = NZPHOS;
    double depth[NZ][NXMEM][NYMEM];
    double po4obsprof[NZPHOS];
    double levitus_depths[NZPHOS] = {0, 10, 20, 30, 50, 75, 100,
	    120, 150, 200, 250, 300, 400, 500, 600,
	    700, 800, 900, 1000, 1100, 1200, 1300,
	    1400, 1500, 1750, 2000, 2500, 3000,
	    3500, 4000, 4500, 5000, 5500};

    printf("\nReading GLODAP DIC and Alk climatology: \n");

    sprintf(infile,"glodap_ic.nc");
    strcpy(inpath, directory);
    strcat(inpath, infile);

    printf("Looking for file '%s'.\n",inpath);

    err = open_input_file(inpath,&file,&cdfid,&timeid);
    if (err != 0) {
	strcat(inpath, ".cdf");
	err = open_input_file(inpath,&file,&cdfid,&timeid);
	if (err != 0) {
	    printf("Unable to find GLODAP file.\n");
	    exit(-73);
	}
    }

    // start with DIC

    if ((status = nc_inq_varid(cdfid, "DIC", &glodapid)))
	ERR(status);

    bzero(start, MAX_NC_VARS * sizeof(long));

    //  count[0] = 1; // no time dimension!
    count[0] = NZPHOS;
    count[1] = NYTOT;
    count[2] = NXTOT;

    tmp3d = alloc3d_f(NZPHOS,NYTOT,NXTOT);
    dictmp = alloc3d(NZPHOS,NXMEM,NYMEM);

    start[0] = 0;

    if ((status = nc_get_vara_float(cdfid,glodapid,start,count,tmp3d[0][0])))
	ERR(status);

    for (k=0;k<NZPHOS;k++) {
	for (i=0;i<NXTOT;i++) {
	    for (j=0;j<NYTOT;j++) {
		dictmp[k][i+2][j+2]= tmp3d[k][j][i]*1.e-3;
	    }
	}
    }
    //BX  temporary bug fix for northern-most row (j=211)
    for (i=2;i<NXMEM;i++) {
	for (k=0;k<NZPHOS;k++) {
	    dictmp[k][i][211] = dictmp[k][i][210];
	}
    }

    //BX  interpolation of GLODAP DIC fields with po4_obs_layer
    //BX   does not work - so we need the entire code here
    //BX po4_obs_layer(h,oxytmp,oxy_init);

    z_depth(h,depth);

#ifdef DIC
    for (i=X1;i<=nx;i++) {
	for (j=Y1;j<=ny;j++) {
	    if (D[i][j]>MINIMUM_DEPTH) {
		for (k=0;k<nzlevitus;k++)
		    po4obsprof[k] = dictmp[k][i][j];
		for (k=0;k<NZ;k++) {
		    dic_init[k][i][j] = lin_interp(depth[k][i][j], po4obsprof,
			    levitus_depths, 0, nzlevitus);
		    if (dic_init[k][i][j] < 0.e0) dic_init[k][i][j] = 0.;
		}
#ifdef MERGED_ML
		// overwrite the top layers
		k=0;
		dic_init[k][i][j] = lin_interp(((depth[k][i][j]+depth[k+1][i][j])/2.), 
			po4obsprof, levitus_depths, 0, nzlevitus);
		if (dic_init[k][i][j] < 0.e0) dic_init[k][i][j] = 0.;
		dic_init[k+1][i][j]=dic_init[k][i][j];
		k=2;
		dic_init[k][i][j] = lin_interp(((depth[k][i][j]+depth[k+1][i][j])/2.), 
			po4obsprof, levitus_depths, 0, nzlevitus);
		if (dic_init[k][i][j] < 0.e0) dic_init[k][i][j] = 0.;
		dic_init[k+1][i][j]=dic_init[k][i][j];
#endif
	    } else {
		for (k=0;k<NZ;k++) {
		    dic_init[k][i][j] = misval;
		}
	    }
	}
    }
#endif /* DIC */


    // Alkalinity is next

    if ((status = nc_inq_varid(cdfid, "TALK", &glodapid)))
	ERR(status);

    bzero(start, MAX_NC_VARS * sizeof(long));

    count[0] = NZPHOS;
    count[1] = NYTOT;
    count[2] = NXTOT;

    start[0] = 0;

    if ((status = nc_get_vara_float(cdfid,glodapid,start,count,tmp3d[0][0])))
	ERR(status);

    // use dictmp here even if it is for alk

    for (k=0;k<NZPHOS;k++) {
	for (i=0;i<NXTOT;i++) {
	    for (j=0;j<NYTOT;j++) {
		dictmp[k][i+2][j+2]= tmp3d[k][j][i]*1.e-3;
	    }
	}
    }
    //BX  temporary bug fix for northern-most row (j=211)
    for (i=2;i<NXMEM;i++) {
	for (k=0;k<NZPHOS;k++) {
	    dictmp[k][i][211] = dictmp[k][i][210];
	}
    }

    //BX  interpolation of GLODAP TAlk fields with po4_obs_layer
    //BX   does not work - so we need the entire code here
    //BX po4_obs_layer(h,oxytmp,oxy_init);

    z_depth(h,depth);

#ifdef DIC
    for (i=X1;i<=nx;i++) {
	for (j=Y1;j<=ny;j++) {
	    if (D[i][j]>MINIMUM_DEPTH) {
		for (k=0;k<nzlevitus;k++)
		    po4obsprof[k] = dictmp[k][i][j];
		for (k=0;k<NZ;k++) {
		    alk_init[k][i][j] = lin_interp(depth[k][i][j], po4obsprof,
			    levitus_depths, 0, nzlevitus);
		    if (alk_init[k][i][j] < 0.e0) alk_init[k][i][j] = 0.;
		}
#ifdef MERGED_ML
		// overwrite the top layers
		k=0;
		alk_init[k][i][j] = lin_interp(((depth[k][i][j]+depth[k+1][i][j])/2.), 
			po4obsprof, levitus_depths, 0, nzlevitus);
		if (alk_init[k][i][j] < 0.e0) alk_init[k][i][j] = 0.;
		alk_init[k+1][i][j]=alk_init[k][i][j];
		k=2;
		alk_init[k][i][j] = lin_interp(((depth[k][i][j]+depth[k+1][i][j])/2.), 
			po4obsprof, levitus_depths, 0, nzlevitus);
		if (alk_init[k][i][j] < 0.e0) alk_init[k][i][j] = 0.;
		alk_init[k+1][i][j]=alk_init[k][i][j];
#endif
	    } else {
		for (k=0;k<NZ;k++) {
		    alk_init[k][i][j] = misval;
		}
	    }
	}
    }
#endif /* DIC */

    free3d_f(tmp3d, NZPHOS);
    free3d(dictmp, NZPHOS);
    close_file(&cdfid,&file);
    }
#endif

void read_ts(int imon,int itts)
    {

#ifdef TS_VARIAB
    /* Read temperature and salinity from time varying records. */
#else
    /* Read temperature and salinity from climatology. */
#endif
    int i,j,k;
    int err, cdfid, timeid;
    char infile[25], inpath[200];
    FILE *file;
    int status;

    int tempclimid,saltclimid,woasaltid;

    size_t start[MAX_NC_VARS];
    size_t count[MAX_NC_VARS];

    float*** tmp3d;
    int inxt,iprv;
    double fact0,fact1,fact2;
    float*** tmp3dp;
    float*** tmp3dn;

# ifdef SEPFILES
    if (imon==NMONTHS) imon=0;
    if (imon==NMONTHS-1) {
	inxt = 0;
	iprv = imon - 1;
    } else if (imon==0) {
	inxt = 1;
	iprv = NMONTHS-1;
    } else {
	inxt = imon + 1;
	iprv = imon - 1;
    }
# else
    if (imon==11) {
	inxt = 0;
	iprv = 10;
    } else if (imon==0) {
	inxt = 1;
	iprv = 11;
    } else {
	inxt = imon + 1;
	iprv = imon - 1;
    }
# endif
    fact1 = 0.5 + ((2.0*(double)itts-1.0) / (2.0*(double)NTSTEP));
    if (fact1 <= 1.0) {
	fact0 = 1.0 - fact1;
	fact2 = 0.0;
    } else {
	fact0 = 0.0;
	fact2 = fact1 - 1.0;
	fact1 = 2.0 - fact1;
    }

    printf("Reading temperature and salinity fields: \n");
#ifdef TS_VARIAB

    sprintf(infile,"ts.nc");
#ifdef TSMEAN
    sprintf(infile,"ts.annual.mean.nc");
#endif

    strcpy(inpath, directory);
    strcat(inpath, infile);

    err = open_input_file(inpath,&file,&cdfid,&timeid);
    if (err != 0) {
	strcat(inpath, ".cdf");
	err = open_input_file(inpath,&file,&cdfid,&timeid);
	if (err != 0) {
	    printf("Unable to find T file.\n");
	    exit(-73);
	}
    }
    bzero(start, MAX_NC_VARS * sizeof(long));

    count[0] = 1; 
    count[1] = NZ;
    count[2] = NYTOT; 
    count[3] = NXTOT;

    printf("read T month=%i\n",imon);
    if ((status = nc_inq_varid(cdfid, "temp", &tempclimid)))
	ERR(status);
    tmp3d  = alloc3d_f(NZ,NYTOT,NXTOT);
    tmp3dn = alloc3d_f(NZ,NYTOT,NXTOT);
    tmp3dp = alloc3d_f(NZ,NYTOT,NXTOT);

    start[0] = imon;
    if ((status = nc_get_vara_float(cdfid,tempclimid,start,count,tmp3d[0][0])))
	ERR(status);
    //printf("read_clim UH month=  %i done\n",imon);
    start[0] = iprv;
    if ((status = nc_get_vara_float(cdfid,tempclimid,start,count,tmp3dp[0][0])))
	ERR(status);
    //printf("read_clim UH month=  %i done\n",iprv);
    start[0] = inxt;
    if ((status = nc_get_vara_float(cdfid,tempclimid,start,count,tmp3dn[0][0])))
	ERR(status);
    //printf("read_clim UH month=  %i done\n",inxt);

    for (k=0;k<NZ;k++)
	for (i=0;i<NXTOT;i++)
	    for (j=0;j<NYTOT;j++)
		Temptm[k][i+2][j+2]= fact1 * tmp3d[k][j][i] +
		fact0 * tmp3dp[k][j][i] + fact2 * tmp3dn[k][j][i];

    close_file(&cdfid,&file);  // for temperature

    //NOW SALT

    sprintf(infile,"ts.nc");
#ifdef TSMEAN
    sprintf(infile,"ts.annual.mean.nc");
#endif
    strcpy(inpath, directory);
    strcat(inpath, infile);

    err = open_input_file(inpath,&file,&cdfid,&timeid);
    if (err != 0) {
	strcat(inpath, ".cdf");
	err = open_input_file(inpath,&file,&cdfid,&timeid);
	if (err != 0) {
	    printf("Unable to find S file.\n");
	    exit(-73);
	}
    }
    bzero(start, MAX_NC_VARS * sizeof(long));

    count[0] = 1; 
    count[1] = NZ;
    count[2] = NYTOT; 
    count[3] = NXTOT;
    printf("read S month=%i\n",imon);
    if ((status = nc_inq_varid(cdfid, "salt", &saltclimid)))
	ERR(status);

    start[0] = imon;
    if ((status = nc_get_vara_float(cdfid,saltclimid,start,count,tmp3d[0][0])))
	ERR(status);
    //printf("read_clim UH month=  %i done\n",imon);
    start[0] = iprv;
    if ((status = nc_get_vara_float(cdfid,saltclimid,start,count,tmp3dp[0][0])))
	ERR(status);
    //printf("read_clim UH month=  %i done\n",iprv);
    start[0] = inxt;
    if ((status = nc_get_vara_float(cdfid,saltclimid,start,count,tmp3dn[0][0])))
	ERR(status);
    //printf("read_clim UH month=  %i done\n",inxt);

    for (k=0;k<NZ;k++)
	for (i=0;i<NXTOT;i++)
	    for (j=0;j<NYTOT;j++)
		Salttm[k][i+2][j+2]= fact1 * tmp3d[k][j][i] +
		fact0 * tmp3dp[k][j][i] + fact2 * tmp3dn[k][j][i];

    /*    // do not use salt_woa from file--use SSS from variability file
	for (i=0;i<NXTOT;i++) {
	  for (j=0;j<NYTOT;j++) {
	    salt_woa[i+2][j+2]= Salttm[0][j][i];
	  }
	  } */

    free3d_f(tmp3d, NZ);
    free3d_f(tmp3dp, NZ);
    free3d_f(tmp3dn, NZ);

    close_file(&cdfid,&file);

#else // TS_VARIAB
    sprintf(infile,"ts.nc");
#ifdef TSMEAN
    sprintf(infile,"ts.annual.mean.nc");
#endif

    strcpy(inpath, directory);
    strcat(inpath, infile);

    err = open_input_file(inpath,&file,&cdfid,&timeid);
    if (err != 0) {
	strcat(inpath, ".cdf");
	err = open_input_file(inpath,&file,&cdfid,&timeid);
	if (err != 0) {
	    printf("Unable to find ts file.\n");
	    exit(-73);
	}
    }

    printf("read T, S, month=%i\n",imon);

    if ((status = nc_inq_varid(cdfid, "salt", &saltclimid)))
	ERR(status);
    if ((status = nc_inq_varid(cdfid, "temp", &tempclimid)))
	ERR(status);

    bzero(start, MAX_NC_VARS * sizeof(long));

    count[0] = 1;
    count[1] = NZ;
    count[2] = NYTOT;
    count[3] = NXTOT;
    //HF
    if (NZ < 24) {
	tmp3d = alloc3d_f(24,NYTOT,NXTOT);
	tmp3dp = alloc3d_f(24,NYTOT,NXTOT);
	tmp3dn = alloc3d_f(24,NYTOT,NXTOT);
    } else {
	tmp3d = alloc3d_f(NZ,NYTOT,NXTOT);
	tmp3dp = alloc3d_f(NZ,NYTOT,NXTOT);
	tmp3dn = alloc3d_f(NZ,NYTOT,NXTOT);
    }

    start[0] = imon;
    if ((status = nc_get_vara_float(cdfid,saltclimid,start,count,tmp3d[0][0])))
	ERR(status);
    start[0] = iprv;
    if ((status = nc_get_vara_float(cdfid,saltclimid,start,count,tmp3dp[0][0])))
	ERR(status);
    start[0] = inxt;
    if ((status = nc_get_vara_float(cdfid,saltclimid,start,count,tmp3dn[0][0])))
	ERR(status);

    for (k=0;k<NZ;k++)
	for (i=0;i<NXTOT;i++)
	    for (j=0;j<NYTOT;j++)
		Salttm[k][i+2][j+2]= fact1 * tmp3d[k][j][i] +
		fact0 * tmp3dp[k][j][i] + fact2 * tmp3dn[k][j][i];

    start[0] = imon;
    if ((status = nc_get_vara_float(cdfid,tempclimid,start,count,tmp3d[0][0])))
	ERR(status);
    start[0] = iprv;
    if ((status = nc_get_vara_float(cdfid,tempclimid,start,count,tmp3dp[0][0])))
	ERR(status);
    start[0] = inxt;
    if ((status = nc_get_vara_float(cdfid,tempclimid,start,count,tmp3dn[0][0])))
	ERR(status);

    for (k=0;k<NZ;k++)
	for (i=0;i<NXTOT;i++)
	    for (j=0;j<NYTOT;j++)
		Temptm[k][i+2][j+2]= fact1 * tmp3d[k][j][i] +
		fact0 * tmp3dp[k][j][i] + fact2 * tmp3dn[k][j][i];

    close_file(&cdfid,&file);

/* ashao: Seems unused
    printf("Reading observed SST, SSS: \n");
    sprintf(infile,"woa01_salt.nc");

    strcpy(inpath, directory);
    strcat(inpath, infile);

    err = open_input_file(inpath,&file,&cdfid,&timeid);
    if (err != 0) {
	strcat(inpath, ".cdf");
	err = open_input_file(inpath,&file,&cdfid,&timeid);
	if (err != 0) {
	    printf("Unable to find world ocean atlas salt file.\n");
	    exit(-73);
	}
    }

    printf("read salt month=%i\n",imon);

    if ((status = nc_inq_varid(cdfid, "WOASALT", &woasaltid)))
	ERR(status);

    bzero(start, MAX_NC_VARS * sizeof(long));

    count[0] = 1;
    count[1] = 24;
    count[2] = NYTOT;
    count[3] = NXTOT;

    start[0] = imon;
    if ((status = nc_get_vara_float(cdfid,woasaltid,start,count,tmp3d[0][0])))
	ERR(status);
    start[0] = iprv;
    if ((status = nc_get_vara_float(cdfid,woasaltid,start,count,tmp3dp[0][0])))
	ERR(status);
    start[0] = inxt;
    if ((status = nc_get_vara_float(cdfid,woasaltid,start,count,tmp3dn[0][0])))
	ERR(status);

    for (i=0;i<NXTOT;i++)
	for (j=0;j<NYTOT;j++)
	    salt_woa[i+2][j+2]= fact1 * tmp3d[0][j][i] +
	    fact0 * tmp3dp[0][j][i] + fact2 * tmp3dn[0][j][i];
    //BX  temporary bug fix for northern-most row (j=211)
    for (i=2;i<=NXMEM;i++) {
	salt_woa[i][211] = salt_woa[i][210];
    }
*/

    if (NZ < 24) {
	free3d_f(tmp3d, 24);
	free3d_f(tmp3dp, 24);
	free3d_f(tmp3dn, 24);
    } else {
	free3d_f(tmp3d, NZ);
	free3d_f(tmp3dp, NZ);
	free3d_f(tmp3dn, NZ);
    }

    close_file(&cdfid,&file);
#endif // TS_VARIAB
    }

#if defined JO2_SPEC || defined UNIVENT
void read_misc_fields(void)
    {
    double dayin, readday;
    char infile[150], inpath[300];
    int cdfid, timeid, err;
    int i,j,k,l;
    FILE *file;
    double work2[NZ][NXMEM][NYMEM];
    char input_dir[75],ocmip_input_dir[75];

    //CAD This code will need to be changed for global domain
    //CAD  if vent/circ/bio components of variability are wanted

#ifdef JO2_SPEC
    /* Read jo2_spec_mon from the 12 binary climatological files. */
    printf("Reading O2 source-sink fields: \n");
    strcpy(ocmip_input_dir,"/bird3/cdeutsch/climdata/");
    for (l=0;l<12;l++) {
	readday = l * 30;
	find_input_file("biosms", "", "", readday, ocmip_input_dir, &file,
		&cdfid,&timeid);
	err = read_time(cdfid, file, timeid, 0, &dayin);
	read_field(cdfid,file,"jo2spec", NXTOT,NYTOT,NZ,
		0,0,0, X1,Y1,0, 0, &work2[0][0][0]);
	//HF	       0,0,0, 1,1,0, 0, &work2[0][0][0]);
	for (j=0;j<=NYMEM-1;j++) {
	    for (i=0;i<=NXMEM-1;i++) {
		for (k=0;k<NZ;k++)
		    jo2_spec_mon[l][k][i][j]=work2[k][i][j];
	    }
	}
    }
#endif

#ifdef UNIVENT
    /* Read aou_spec_mon from the 12 binary climatological files. */
    printf("Reading climatological AOU fields: \n");
    strcpy(ocmip_input_dir,"/bird3/cdeutsch/climdata/");
    for (l=0;l<12;l++) {
	readday = l * 30;
	find_input_file("biosms", "", "", readday, ocmip_input_dir, &file,
		&cdfid,&timeid);
	err = read_time(cdfid, file, timeid, 0, &dayin);
	read_field(cdfid,file,"aouspec", NXTOT,NYTOT,NZ,
		0,0,0, X1,Y1,0, 0, &work2[0][0][0]);
	//HF	       0,0,0, 1,1,0, 0, &work2[0][0][0]);
	for (j=0;j<=NYMEM-1;j++) {
	    for (i=0;i<=NXMEM-1;i++) {
		for (k=0;k<NZ;k++)
		    aou_spec_mon[l][k][i][j]=work2[k][i][j];
	    }
	}
    }
#endif
    }
#endif


#ifdef RESTART
void read_tracer_init(int imon, char *run_name)
#else
void read_tracer_init(int imon)
#endif
    {
    int i,j,k;
#if defined RESTART || defined INERT
    char infile[150], inpath[300];
    int cdfid, timeid, err;
    FILE *file;
#endif

    /* open file containing tracer initial conditions */
    /* see prep_restart_fields.m to make fields for a restart */
    /* for (i=0;i<=NXMEM-1;i++)
    for (j=0;j<=NYMEM-1;j++)
    printf("begin read_t_init D[%d][%d]=%g\n",i,j,D[i][j]); */

#ifdef RESTART

    /*   OPEN RESTART FILE    */

    sprintf(infile,"restart.%s.%04i.nc",run_name,imon);
    printf("Reading NetCDF restart file '%s'.\n\n",infile);

    err = open_input_file(infile,&file,&cdfid,&timeid);
    if (err != 0) {
	strcat(infile, ".cdf");
	err = open_input_file(infile,&file,&cdfid,&timeid);
	if (err != 0) {
	    printf("Unable to find restart file.\n");
	    exit(-73);
	}
    }

    /* Read tracer initial conditions as appropriate. */

# ifdef ENTRAIN
    if (rflags[6]) read_field(cdfid,file,"mn_ea", NXTOT,NYTOT,NZ,
	    0,0,0, X1,Y1,0, 0, &ea_init[0][0][0]);
    //HF             0,0,0, 1,1,0, 0, &ea_init[0][0][0]);
    if (rflags[8]) read_field(cdfid,file,"mn_eaml", NXTOT,NYTOT,NZ,
	    0,0,0, X1,Y1,0, 0, &eaml_init[0][0]);
    //HF             0,0,0, 1,1,0, 0, &eaml_init[0][0]);
# endif

# ifdef SMOOTH_REST
    if (rflags[3]) read_field(cdfid,file,"mn_h", NXTOT,NYTOT,NZ,
	    0,0,0, X1,Y1,0, 0, &h_init[0][0][0]);
    //HF             0,0,0, 1,1,0, 0, &h_init[0][0][0]);
# endif

# ifdef OXYGEN
    if (rflags[10]) read_field(cdfid,file,"mn_oxygen", NXTOT,NYTOT,NZ,
	    0,0,0, X1,Y1,0, 0, &oxy_init[0][0][0]);
    //HF             0,0,0, 1,1,0, 0, &oxy_init[0][0][0]);
#  ifdef OXY18
    if (rflags[25]) read_field(cdfid,file,"mn_o18", NXTOT,NYTOT,NZ,
	    0,0,0, X1,Y1,0, 0, &o18_init[0][0][0]);
    //HF             0,0,0, 1,1,0, 0, &o18_init[0][0][0]);
#  endif 
# endif

# ifdef PHOSPHATE
    if (rflags[19]) read_field(cdfid,file,"mn_phos", NXTOT,NYTOT,NZ,
	    0,0,0, X1,Y1,0, 0, &phosphate_init[0][0][0]);
    //HF	     0,0,0, 1,1,0, 0, &phosphate_init[0][0][0]);
    if (rflags[22]) read_field(cdfid,file,"mn_dop", NXTOT,NYTOT,NZ,
	    0,0,0, X1,Y1,0, 0, &dop_init[0][0][0]);
    //HF	     0,0,0, 1,1,0, 0, &dop_init[0][0][0]);
# endif

# ifdef NITRATE
    if (rflags[20]) read_field(cdfid,file,"mn_nitr", NXTOT,NYTOT,NZ,
	    0,0,0, X1,Y1,0, 0, &nitrate_init[0][0][0]);
    //HF	     0,0,0, 1,1,0, 0, &nitrate_init[0][0][0]);
    if (rflags[23]) read_field(cdfid,file,"mn_don", NXTOT,NYTOT,NZ,
	    0,0,0, X1,Y1,0, 0, &don_init[0][0][0]);
    //HF	     0,0,0, 1,1,0, 0, &don_init[0][0][0]);
#  ifdef N15_CYCLE
    if (rflags[39]) read_field(cdfid,file,"mn_nitr_15n", NXTOT,NYTOT,NZ,
	    0,0,0, X1,Y1,0, 0, &nitrate_init_15n[0][0][0]);
    if (rflags[40]) read_field(cdfid,file,"mn_don_15n", NXTOT,NYTOT,NZ,
	    0,0,0, X1,Y1,0, 0, &don_init_15n[0][0][0]);
#  endif
# endif /* NITRATE */

# ifdef INERT
    if (rflags[28]) read_field(cdfid,file,"mn_inert1", NXTOT,NYTOT,NZ,
	    0,0,0, X1,Y1,0, 0, &inert1_init[0][0][0]);
    //HF 0,0,0, 1,1,0, 0, &inert1_init[0][0][0]);
    if (rflags[29]) read_field(cdfid,file,"mn_inert2", NXTOT,NYTOT,NZ,
	    0,0,0, X1,Y1,0, 0, &inert2_init[0][0][0]);
    //HF             0,0,0, 1,1,0, 0, &inert2_init[0][0][0]);
# endif

# ifdef AGE  /* all zeros */
    if (rflags[9]) read_field(cdfid,file,"mn_age", NXTOT,NYTOT,NZ,
	    0,0,0, X1,Y1,0, 0, &age_init[0][0][0]);
    //HF             0,0,0, 1,1,0, 0, &age_init[0][0][0]);
# endif

# ifdef DIC
    if (rflags[21]) read_field(cdfid,file,"mn_dic", NXTOT,NYTOT,NZ,
	    0,0,0, X1,Y1,0, 0, &dic_init[0][0][0]);
    //HF             0,0,0, 1,1,0, 0, &dic_init[0][0][0]);
    if (rflags[31]) read_field(cdfid,file,"mn_alk", NXTOT,NYTOT,NZ,
	    0,0,0, X1,Y1,0, 0, &alk_init[0][0][0]);
    //HF             0,0,0, 1,1,0, 0, &alk_init[0][0][0]);
# endif

# ifdef CFC  /* all zeros */
    if (rflags[14]) read_field(cdfid,file,"mn_cfc11", NXTOT,NYTOT,NZ,
	    0,0,0, X1,Y1,0, 0, &cfc11_init[0][0][0]);
    //HF             0,0,0, 1,1,0, 0, &alk_init[0][0][0]);
    if (rflags[15]) read_field(cdfid,file,"mn_cfc12", NXTOT,NYTOT,NZ,
	    0,0,0, X1,Y1,0, 0, &cfc12_init[0][0][0]);
    //HF             0,0,0, 1,1,0, 0, &cfc12_init[0][0][0]);
# endif
    // begin ashao
# ifdef SF6  /* all zeros */
    if (rflags[50]) read_field(cdfid,file,"mn_cfc11", NXTOT,NYTOT,NZ,
	    0,0,0, X1,Y1,0, 0, &sf6_init[0][0][0]);
    //HF             0,0,0, 1,1,0, 0, &alk_init[0][0][0]);
# endif
    // end ashao

#else  // RESTART

    /* Read tracer initial conditions as appropriate. */

# ifdef OXYGEN

#  ifdef LEV_OXY
    read_oxy_ic();
#  endif
    for (j=0;j<=NYMEM-1;j++) {
	for (i=0;i<=NXMEM-1;i++) {
	    for (k=0;k<=NZ-1;k++) {
#  ifndef LEV_OXY
		oxy_init[k][i][j] = 300.e-3;
#  endif
#  ifdef OXY18
		o18_init[k][i][j] = oxy_init[k][i][j] * 2.005e-3 * (1.e0 + 24.e-3);
#  endif
	    }
	}
    }
# endif

# ifdef PHOSPHATE
    // Try initializing with observed po4
    // start for time step "1"
    read_biotic_bc(imon,1);

    printf("Initialize 2 for month %i.\n",imon);

    conc_obs_layer(h,po4_star_lev,po4_star_lay);

    printf("Initialize 3 for month %i.\n",imon);
    for (i=0;i<=NXMEM-1;i++) {
	for (j=0;j<=NYMEM-1;j++) {
	    if (D[i][j]>MINIMUM_DEPTH) {
		for (k=0;k<=NZ-1;k++)
		    phosphate_init[k][i][j]   = po4_star_lay[k][i][j];
	    } else {
		for (k=0;k<=NZ-1;k++)
		    phosphate_init[k][i][j]   = misval;
	    }
	}
    }


    set_darray3d_zero(dop_init, NZ, NXMEM, NYMEM);
# endif /* PHOSPHATE */

# ifdef NITRATE
    // Try initializing with observed po4

    printf("Initialize NO3: Step 2 for month %i.\n",imon);

    conc_obs_layer(h,no3_lev,no3_lay);

    printf("Initialize NO3: Step 3 for month %i.\n",imon);
    for (i=0;i<=NXMEM-1;i++) {
	for (j=0;j<=NYMEM-1;j++) {
	    if (D[i][j]>MINIMUM_DEPTH) {
		for (k=0;k<=NZ-1;k++)
		    nitrate_init[k][i][j]   = no3_lay[k][i][j];
	    } else {
		for (k=0;k<=NZ-1;k++)
		    nitrate_init[k][i][j]   = misval;
	    }
	}
    }
    printf("nitrate_init=%g\n",nitrate_init[0][55][25]);
    set_darray3d_zero(don_init, NZ, NXMEM, NYMEM);
# endif /* NITRATE */

# ifdef INERT  
#  ifdef INERT_UNUSED
    /* use sponge as i.c. for inert1 */
    //HF strcpy(inpath,"/bird3/cdeutsch/input/sponge_inert.cdf");
    err = open_input_file(inpath,&file,&cdfid,&timeid);
    read_field(cdfid,file,"init_argon",NXTOT,NYTOT,NZ,
	    0,0,0, X1,Y1,0, 0, &inert1_init[0][0][0]);
    //HF             0,0,0, 1,1,0, 0, &inert1_init[0][0][0]);
    read_field(cdfid,file,"init_temp",NXTOT,NYTOT,NZ,
	    0,0,0, X1,Y1,0, 0, &inert2_init[0][0][0]);
    //HF             0,0,0, 1,1,0, 0, &inert2_init[0][0][0]);
    printf("initial inert1 = %g \n", inert1_init[12][40][100]);
    printf("initial inert2 = %g \n", inert2_init[12][40][100]);

#  endif
    set_fix_darray3d_zero(inert1_init,NZ);
    set_fix_darray3d_zero(inert2_init,NZ);
#endif

# ifdef AGE  /* all zeros */
    for (j=0;j<=NYMEM-1;j++) {
	for (i=0;i<=NXMEM-1;i++) {
	    for (k=0;k<=NZ-1;k++)
		age_init[k][i][j] = 0.0;
	}
    }
# endif


# ifdef DIC
#  ifdef GLODAP_IC
    read_glodap_ic();
#  else
    //strcpy(inpath,"/home/cdeutsch/OFF/INPUT/init_dic_1deg.cdf");
    //err = open_input_file(inpath,&file,&cdfid,&timeid);
    //read_field(cdfid,file,"dic", NXTOT,NYTOT,NZ,
    //           0,0,0, 1,1,0, 0, &dic_init[0][0][0]);
    //strcpy(inpath,"/home/cdeutsch/OFF/INPUT/init_dic_1deg.cdf");
    //err = open_input_file(inpath,&file,&cdfid,&timeid);
    //read_field(cdfid,file,"dic", NXTOT,NYTOT,NZ,
    //           0,0,0, 1,1,0, 0, &alk_init[0][0][0]);
    for (k=0;k<=NZ-1;k++) {
	for (i=0;i<=NXMEM-1;i++) {
	    for (j=0;j<=NYMEM-1;j++) {
		dic_init[k][i][j] = 2250.e-3;   //  model calculates in mol/m^3 not mmol/m^3
		alk_init[k][i][j] = dic_init[k][i][j] + 100.e-3;   //
	    }
	}
    }
#  endif
# endif

# ifdef CFC  /* all zeros */
    for (j=0;j<=NYMEM-1;j++) {
	for (i=0;i<=NXMEM-1;i++) {
	    for (k=0;k<=NZ-1;k++) {
		cfc11_init[k][i][j] = 0.0;
		cfc12_init[k][i][j] = 0.0;
	    }
	}
    }
# endif

    // begin ashao
# ifdef SF6  /* all zeros */
    for (j=0;j<=NYMEM-1;j++) {
	for (i=0;i<=NXMEM-1;i++) {
	    for (k=0;k<=NZ-1;k++) {
		sf6_init[k][i][j] = 0.0;
	    }
	}
    }
# endif
    // end ashao
#endif  //  RESTART

    }

void read_restoring(int imon,int itts)
    {
#ifdef REST_ARCTIC
# if defined DIC || defined OXYGEN
    //BX  Read in CCSM derived climatology for DIC and ALK used for restoring
    //BX  It is assumed that the infile has NZPHOS depth levels
    //BX  28NOV07 BX
    int i,j,k;
    int err, cdfid, timeid;
    char infile[25], inpath[200];
    FILE *file;
    int status;

    int dicid,alkid,oxyid;

    size_t start[MAX_NC_VARS];
    size_t count[MAX_NC_VARS];

    double var_lev[NZPHOS][NXMEM][NYMEM];

    float*** tmp3d;
    int inxt,iprv;
    double fact0,fact1,fact2;
    float*** tmp3dp;
    float*** tmp3dn;

    if (imon==11) {
	inxt = 0;
	iprv = 10;
    } else if (imon==0) {
	inxt = 1;
	iprv = 11;
    } else {
	inxt = imon + 1;
	iprv = imon - 1;
    }
    fact1 = 0.5 + ((2.0*(double)itts-1.0) / (2.0*(double)NTSTEP));
    if (fact1 <= 1.0) {
	fact0 = 1.0 - fact1;
	fact2 = 0.0;
    } else {
	fact0 = 0.0;
	fact2 = fact1 - 1.0;
	fact1 = 2.0 - fact1;
    }

#  ifdef DIC
    printf("Reading CCSM restoring fields for DIC and ALK: \n");
#  endif
#  ifdef OXYGEN
    printf("Reading CCSM restoring fields for OXYGEN: \n");
#  endif
    //BX  sprintf(infile,"CCSM_PI_CLIM_DIC_ALK.nc");
    sprintf(infile,"CCSM_PI_CLIM_DIC_ALK_O2.nc");

    strcpy(inpath, directory);
    strcat(inpath, infile);

    err = open_input_file(inpath,&file,&cdfid,&timeid);
    if (err != 0) {
	strcat(inpath, ".cdf");
	err = open_input_file(inpath,&file,&cdfid,&timeid);
	if (err != 0) {
	    printf("Unable to find CCSM restoring file.\n");
	    exit(-73);
	}
    }

    printf("read restoring file month=%i\n",imon);

    /*   Read DIC first - use WOAPO4 field to save memory  */

#  ifdef DIC
    if ((status = nc_inq_varid(cdfid, "RESTDIC", &dicid)))
	ERR(status);
    if ((status = nc_inq_varid(cdfid, "RESTALK", &alkid)))
	ERR(status);
#  endif
#  ifdef OXYGEN
    if ((status = nc_inq_varid(cdfid, "RESTO2", &oxyid)))
	ERR(status);
#  endif

    bzero(start, MAX_NC_VARS * sizeof(long));

    count[0] = 1;
    count[1] = NZPHOS;
    count[2] = NYTOT;
    count[3] = NXTOT;

    tmp3d = alloc3d_f(NZPHOS,NYTOT,NXTOT);
    tmp3dp = alloc3d_f(NZPHOS,NYTOT,NXTOT);
    tmp3dn = alloc3d_f(NZPHOS,NYTOT,NXTOT);
#  ifdef DIC
    start[0] = imon;
    if ((status = nc_get_vara_float(cdfid,dicid,start,count,tmp3d[0][0])))
	ERR(status);
    start[0] = iprv;
    if ((status = nc_get_vara_float(cdfid,dicid,start,count,tmp3dp[0][0])))
	ERR(status);
    start[0] = inxt;
    if ((status = nc_get_vara_float(cdfid,dicid,start,count,tmp3dn[0][0])))
	ERR(status);

    for (k=0;k<NZPHOS;k++)
	for (i=0;i<NXTOT;i++)
	    for (j=0;j<NYTOT;j++)
		var_lev[k][i+2][j+2]= fact1 * tmp3d[k][j][i] +
		fact0 * tmp3dp[k][j][i] + fact2 * tmp3dn[k][j][i];

    //BX  temporary bug fix for northern-most row (j=211)
    for (k=0;k<NZPHOS;k++)
	for (i=2;i<=NXMEM;i++)
	    var_lev[k][i][211] = var_lev[k][i][210];

    /*  Regrid restoring fields to model grid
     (conc_obs_layer does regridding not just for po4_star) */
    conc_obs_layer(h,var_lev,dic_lay);

    start[0] = imon;
    if ((status = nc_get_vara_float(cdfid,alkid,start,count,tmp3d[0][0])))
	ERR(status);
    start[0] = iprv;
    if ((status = nc_get_vara_float(cdfid,alkid,start,count,tmp3dp[0][0])))
	ERR(status);
    start[0] = inxt;
    if ((status = nc_get_vara_float(cdfid,alkid,start,count,tmp3dn[0][0])))
	ERR(status);

    for (k=0;k<NZPHOS;k++)
	for (i=0;i<NXTOT;i++)
	    for (j=0;j<NYTOT;j++)
		var_lev[k][i+2][j+2]= fact1 * tmp3d[k][j][i] +
		fact0 * tmp3dp[k][j][i] + fact2 * tmp3dn[k][j][i];

    //BX  temporary bug fix for northern-most row (j=211)
    for (i=2;i<=NXMEM;i++) {
	for (k=0;k<NZPHOS;k++) {
	    var_lev[k][i][211] = var_lev[k][i][210];
	}
    }

    /*  Regrid restoring fields to model grid
     (conc_obs_layer does regridding not just for po4_star) */

    conc_obs_layer(h,var_lev,alk_lay);
#  endif    // DIC (and Alk)

    /*  Now Oxygen  */
#  ifdef OXYGEN
    start[0] = imon;
    if ((status = nc_get_vara_float(cdfid,oxyid,start,count,tmp3d[0][0])))
	ERR(status);
    start[0] = iprv;
    if ((status = nc_get_vara_float(cdfid,oxyid,start,count,tmp3dp[0][0])))
	ERR(status);
    start[0] = inxt;
    if ((status = nc_get_vara_float(cdfid,oxyid,start,count,tmp3dn[0][0])))
	ERR(status);

    for (k=0;k<NZPHOS;k++)
	for (i=0;i<NXTOT;i++)
	    for (j=0;j<NYTOT;j++)
		var_lev[k][i+2][j+2]= fact1 * tmp3d[k][j][i] +
		fact0 * tmp3dp[k][j][i] + fact2 * tmp3dn[k][j][i];

    //BX  temporary bug fix for northern-most row (j=211)
    for (i=2;i<=NXMEM;i++) {
	for (k=0;k<NZPHOS;k++) {
	    var_lev[k][i][211] = var_lev[k][i][210];
	}
    }

    /*  Regrid restoring fields to model grid
     (conc_obs_layer does regridding not just for po4_star) */

    conc_obs_layer(h,var_lev,oxy_lay);
#  endif    // OXYGEN

    close_file(&cdfid,&file);

    free3d_f(tmp3d, NZPHOS);
    free3d_f(tmp3dp, NZPHOS);
    free3d_f(tmp3dn, NZPHOS);
# endif
#endif
    }

void read_buoy(int imon)
    {

    }

#ifdef SPONGE
void read_sponge(void)
    {
    char infile[150], inpath[300];
    int cdfid, timeid, err;
    int i,j,k;
    FILE *file;

    /* ******************************************************** */
    /*                                                          */
    /* IF THIS SUBROUTINE IS TO BE USED MAKE SURE FILES WILL BE */
    /* CLOSED!!!!!  bx 22OCT07                                  */
    /*                                                          */
    /* ******************************************************** */

    strcpy(inpath,"/dsk1/suzanne/OFF/input/sponge_1deg.cdf");
    err = open_input_file(inpath,&file,&cdfid,&timeid);

    WARNING NOT YET FIXED
    read_field(cdfid,file,"IDAMP",NXTOT,NYTOT,1,
	    0,0,0, 1,1,0, 0, &Idamp[0][0]);
    read_field(cdfid,file,"E",NXTOT,NYTOT,NZ,
	    0,0,0, 1,1,0, 0, &sp_e[0][0][0]);

#if defined(AGE) || defined(CFC) || defined(SF6) // ashao
    read_field(cdfid,file,"SP_AGE",NXTOT,NYTOT,NZ,
	    0,0,0, 1,1,0, 0, &sp_age[0][0][0]);
#endif
#ifdef OXYGEN
    read_field(cdfid,file,"SP_OXYGEN",NXTOT,NYTOT,NZ,
	    0,0,0, 1,1,0, 0, &sp_oxygen[0][0][0]);
    for (j=0;j<=NYMEM-1;j++) {
	for (i=0;i<=NXMEM-1;i++) {
	    for (k=0;k<=NZ-1;k++)
		sp_oxygen[k][i][j] = ox_init[k][i][j];
#ifdef OXY18
	    sp_o18[k][i][j] = o18_init[k][i][j];
#endif
	}
    }
#endif
#ifdef INERT_UNUSED
    strcpy(inpath,"/bird3/cdeutsch/OFF/input/sponge_argon.cdf");
    err = open_input_file(inpath,&file,&cdfid,&timeid);
    read_field(cdfid,file,"init_argon",NXTOT,NYTOT,NZ,
	    0,0,0, 1,1,0, 0, &sp_inert1[0][0][0]);
    read_field(cdfid,file,"init_temp",NXTOT,NYTOT,NZ,
	    0,0,0, 1,1,0, 0, &sp_inert2[0][0][0]);
#endif

#ifdef PHOSPHATE
    strcpy(inpath,"/dsk1/suzanne/OFF/input/init_po4_1deg.cdf");
    err = open_input_file(inpath,&file,&cdfid,&timeid);
    read_field(cdfid,file,"sp_po4",NXTOT,NYTOT,NZ,
	    0,0,0, 1,1,0, 0, &sp_phosphate[0][0][0]);
    for (j=0;j<=NYMEM-1;j++) {
	for (i=0;i<=NXMEM-1;i++) {
	    for (k=0;k<=8;k++)
		sp_dop[k][i][j] = 0.2e0;
	    for (k=9;k<=NZ-1;k++)
		sp_dop[k][i][j] = 0.0e0;
	}
    }
#endif

#ifdef DIC
    strcpy(inpath,"/dsk1/suzanne/OFF/input/init_dic_1deg.cdf");
    err = open_input_file(inpath,&file,&cdfid,&timeid);
    read_field(cdfid,file,"sp_dic",NXTOT,NYTOT,NZ,
	    0,0,0, 1,1,0, 0, &sp_dic[0][0][0]);
    strcpy(inpath,"/dsk1/suzanne/OFF/input/init_alk_1deg.cdf");
    err = open_input_file(inpath,&file,&cdfid,&timeid);
    read_field(cdfid,file,"sp_alk",NXTOT,NYTOT,NZ,
	    0,0,0, 1,1,0, 0, &sp_alk[0][0][0]);
#endif
    }
#endif /* SPONGE */


// begin ashao
# if defined(SF6) || defined(CFC)
// Read in tracer year, concentrations, and coefficients for schmidt and solubility from csv file
int read_trac(double sc_coeffs[4], double sol_coeffs[7], int tracyear[MAXRECS], double Ncon[MAXRECS], double Scon[MAXRECS])
    {

    int     i,tracrecs;
    FILE    *infile;
    char    line[200],fname[500];
    double  inval;

    printf("Name of tracer data file: \n");
    scanf("%s",fname);
    printf("%s\n",fname);
    infile=fopen(fname,"r");

    fgets(line,sizeof(line),infile);
    sscanf(line,"%lf",&inval);
    tracrecs=(int) floor(inval);


    printf("TRACRECS=%d\n",tracrecs);

    fgets(line,sizeof(line),infile);
    sscanf(line,"%lf%lf%lf\t%lf",&sc_coeffs[0],&sc_coeffs[1],&sc_coeffs[2],&sc_coeffs[3]);
    //	sscanf(line,"%lf",sc_coeffs[0]);

    fgets(line,sizeof(line),infile);
    sscanf(line,"%lf%lf%lf%lf%lf%lf%lf",&sol_coeffs[0],&sol_coeffs[1],&sol_coeffs[2],&sol_coeffs[3],&sol_coeffs[4],&sol_coeffs[5],&sol_coeffs[6]);

    for (i=0;i<tracrecs;i++)
	{
	fgets(line,sizeof(line),infile);
	sscanf(line,"%lf%lf%lf",&inval,&Ncon[i],&Scon[i]);
	tracyear[i]=(int) floor(inval);
	}

    printf("%d\n",tracrecs);
    return tracrecs;

    }
# endif
// end ashao

/* void read_patch
 * Read in which grid cells to patch for estimation of TTD
 */
void read_patch( double ttd_mask[NXMEM][NYMEM] )
    {

    /* Read ocmip gas exchange parameter fields */
    int err, cdfid, timeid;
    char infile[25], inpath[200];
    FILE *file;
    int status;

    int patchid;


    printf("Enter name of file with TTD patchmask\n");
    gets(infile);

    strcpy(inpath, directory);
    strcat(inpath, infile);
    printf("Reading TTD Patch: %s\n",inpath);

    err = open_input_file(inpath,&file,&cdfid,&timeid);
    if (err != 0) {
	strcat(inpath, ".cdf");
	err = open_input_file(inpath,&file,&cdfid,&timeid);
	if (err != 0) {
	    printf("Unable to find TTD Patch file.\n");
	    exit(-73);
	}
    }

    read_field(cdfid,file,"mask", NXTOT,NYTOT,1, 0,0,0, 0,0,0, 1, ttd_mask[0]);
    close_file(&cdfid,&file);

    return;

    }

/* ashao: Read data routines (UH, VH, WD, T, S) for hindcast runs */
void read_hind( int cmon )
    {
    extern double currtime;
    char uhfile[100], vhfile[100], wdfile[100], tsfile[100], hfile[100];
    int readclim; // If 0: read hindcast. If 1: read climatology
    int imon, inxt, iprv;

    int i,j,k, ii;
    int err, cdfid, timeid;
    FILE *file;
    int status;

    int uhclimid,vhclimid,wdclimid;

    strcpy(uhfile,directory);
    strcpy(vhfile,directory);
    strcpy(wdfile,directory);
    strcpy(hfile,directory);
    strcpy(tsfile,directory);
    // Read different files for hindcast/climatology
    // Calculate different monthly index
    if ((currtime>BEGHIND) & (currtime<(BEGHIND+MAXYEARS+1))) { // Hindcast
	readclim=0;
	strcat(uhfile,"UH-hind.nc");
	strcat(vhfile,"VH-hind.nc");
	strcat(wdfile,"WD-hind.nc");
	strcat(tsfile,"ts-hind.nc");
	strcat(hfile,"H-hind.nc");
	imon=cmon-SPINUP*12; // Set index to months since start of spinup
	inxt=imon+1;
	if (imon==0) iprv=MAXYEARS*12-1;
	else iprv=imon-1;
    }
    else { // Climatology
	strcat(uhfile,"UH-clim.nc");
	strcat(vhfile,"VH-clim.nc");
	strcat(wdfile,"WD-clim.nc");
	strcat(tsfile,"ts-clim.nc");
	strcat(hfile,"H-clim.nc");

	imon=cmon%NMONTHS;
	inxt=(cmon+1)%NMONTHS;
	if (imon==0) iprv=NMONTHS-1;
	else iprv=imon-1;

    }

    read_var3d(uhfile,"uh",inxt,uhtm);
    read_var3d(vhfile,"vh",inxt,vhtm);
    read_var3d(wdfile,"wd",inxt,wd);
    read_var3d(hfile,"h",imon,h);
    read_var3d(hfile,"h",imon,hstart);
    read_var3d(hfile,"h",inxt,hend);
    read_var3d(tsfile,"temp",inxt,Temptm);
    read_var3d(tsfile,"salt",inxt,Salttm);

    // zonal re-entrance
    for (k=0;k<NZ;k++) {
        for (j=0;j<=NYMEM-1;j++) {
            h[k][nx+1][j] = h[k][2][j];
            h[k][nx+2][j] = h[k][3][j];
            h[k][0][j] =   h[k][nx-1][j];
            h[k][1][j] =   h[k][nx][j];
            hstart[k][nx+1][j] = hstart[k][2][j];
            hstart[k][nx+2][j] = hstart[k][3][j];
            hstart[k][0][j] =   hstart[k][nx-1][j];
            hstart[k][1][j] =   hstart[k][nx][j];
            hend[k][nx+1][j] = hend[k][2][j];
            hend[k][nx+2][j] = hend[k][3][j];
            hend[k][0][j] =   hend[k][nx-1][j];
            hend[k][1][j] =   hend[k][nx][j];
        }
    }

    //      meridional re-entrance
    for (i=2;i<=nx;i++) {
        ii = 363 - i;
        for (k=0;k<NZ;k++) {
            h[k][ii][ny+1] = h[k][i][ny];
            h[k][ii][ny+2]   = h[k][i][ny-1];
            hstart[k][ii][ny+1] = hstart[k][i][ny];
            hstart[k][ii][ny+2] = hstart[k][i][ny-1];
            hend[k][ii][ny+1] = hend[k][i][ny];
            hend[k][ii][ny+2] = hend[k][i][ny-1];
        }
    }




    }

void read_var3d( char inpath[200], char varname[200], int imon, double ***data)
    {

    int i,j,k;
    int err, cdfid, timeid, varid;
    char infile[25];
    FILE *file;
    int status;
    int inxt,iprv;
    size_t start[MAX_NC_VARS];
    size_t count[MAX_NC_VARS];
    float ***tmp3d;

    start[0] = imon;
    start[1] = 0;
    start[2] = 0;
    start[3] = 0;

    err = open_input_file(inpath,&file,&cdfid,&timeid);
    if ((status = nc_inq_varid(cdfid, varname, &varid))) 
    bzero(start, MAX_NC_VARS * sizeof(long));

    count[0] = 1;
    if (strcmp(varname,"wd")) count[1]=NZ;
    else count[1] = NZ;
    count[2] = NYTOT;
    count[3] = NXTOT;

    tmp3d  = alloc3d_f(NZ+1,NYTOT,NXTOT);
    if ((status = nc_get_vara_float(cdfid,varid,start,count,tmp3d[0][0])))
	ERR(status);
    for (k=0;k<count[1];k++)
	for (i=0;i<NXTOT;i++)
	    for (j=0;j<NYTOT;j++)
		data[k][i+2][j+2]= tmp3d[k][j][i];

    }
