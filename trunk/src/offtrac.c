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
#include "setup_tracers.h"
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


//BX-e

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
double hstart[NZ][NXMEM][NYMEM];
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
int wetmask[NXMEM][NYMEM];

int ksub[NXMEM][NYMEM];
int ksub_clim[NXMEM][NYMEM];

int month;
int mlen[12];
int icnt;

double qlat[NYMEM];
double hlat[NYMEM];

/* Begin added DT */
int beginyear;
int theyear;
int lastsave;
/*  End added DT  */


#ifdef AGE
//HF  double age[NZ][NXMEM][NYMEM];
double ***age;
double age_init[NZ][NXMEM][NYMEM];
double ***mn_age;
int mAGE;
#endif

#ifdef SMOOTH_REST
//double ***h_init;
double h_init[NZ][NXMEM][NYMEM];
#endif



// begin ashao
double currtime;
int iterno;
#ifdef HINDCAST
int usehindcast = 1;
#else
int usehindcast = 0;
#endif

int hindindex = 0;
// Calculate what iteration to start/stop reading hindcast
int starthindindex = (BEGHIND-BEGYEAR)*12; //cmon added later
int numhindmonths = (ENDHIND-BEGHIND+1)*12;
// end ashao

double k0, k1, k2, kb, kw, ks, kf, k1p, k2p, k3p, ksi;
double ff, htotal, htotal2, xco2, xacc, x1, x2;
double bt, st, ft, sit, pt, dic1, ta, co2starair;

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
	double yearday;
	int imon, inxt, ilst;
# ifndef WRTTS
int itts; /* tracer time step counter */
# endif
int nmnfirst;
#ifdef SEPFILES
double smon, snxt;
int ismon, isnxt, ihnxt;
#endif

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

	printf("Enter base name for output\n");

	fgets(run_name, sizeof(run_name), stdin);
	run_name[strlen(run_name) - 1] = '\0';
	sprintf(output_filename, "%s.%04g.nc", run_name, inmon + tmon - 1);
	printf("Create NetCDF output file '%s'.\n", output_filename);



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
smon = (double) ((int) inmon % NMONTHS);
mon = 12 * modf(dyr, iyr);

printf("\n initial mon = %g - %g - %g \n\n", inmon, smon, mon);


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

ismon = (int) (smon + 0.00001) % NMONTHS;

isnxt = (ismon+1) % NMONTHS;
ihnxt = (ismon+2) % NMONTHS;

dbmon = (double) inmon;
lst = 12 * modf((dbmon - 1 + 12) / 12, iyr);
ilst = (int) (lst + 0.00001);

// ashao: Read in next month's isopycnal thickness fields
// (will be copied at the beginning of the integration)
// Done this way so that if hindcasts are used the physical fields switch smoothly
// to/from climatological fields
//BX  files are not in regular order (uvw are midmonth and h starts with last month)

currtime = BEGYEAR;
if (usehindcast) {
	// Check to see if simulation started within the hindcast years
	if ( (currtime >= BEGHIND) || (currtime < (ENDHIND+1) ) ) {
		hindindex=inmon;
		read_h(ismon,hend,"hind");
	}
}
else {
	read_h(ismon, hend,"clim");
}
// for files in regular order (h before uvw) use code here
//BX
// read_uvw(imon,1);
//BX
// read_h(imon,inxt);



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




/*-------------------------------------------------------------------*
 *
 *     begin integration loop
 *
 *-------------------------------------------------------------------*/

mon = 0.0; /* reiniti */
nmnfirst = 1;


for (cmon = inmon; cmon < inmon + tmon; cmon++)
{
	nmn++;

	dyr = cmon / 12.0;
	ndyr = (cmon + 1) / 12.0;

	dbmon = (double) cmon;
	lst = 12 * modf((dbmon - 1 + 12) / 12, iyr);
	ilst = (int) (lst + 0.00001);
	printf("double-mon=%g,lastmon=%g,ilst=%i\n", dbmon, lst, ilst);


	smon = (double) ((int) cmon % NMONTHS);
	snxt = (double) ((int) (cmon + 1) % NMONTHS);

	mon = 12.0 * modf(dyr, iyr);
	nxt = 12.0 * modf(ndyr, nyr);

	printf("the current month is %i-%g-%g \n", cmon, smon, mon);
	printf("the current year is %g \n", *iyr);

	imon = (int) (mon + 0.00001);
	inxt = (int) (nxt + 0.00001);

	ismon = (int) (smon + 0.00001);
	isnxt = (int) (snxt + 0.00001);

	yearday = 0;
	for (i=0;i<=imon;i++) yearday += dmon[imon];
	currtime = *iyr + BEGYEAR + yearday/365.0;

	day = (*iyr) * 365.0 + dmon[imon];
	*dy = currtime;

	printf("the current day is -%g- \n", *dy);
	printf("the current ismon/smon is -%i-%g- \n", ismon, smon);
	printf("the next smonth/mean index: -%i- \n", isnxt);

	dt = ((double) mlen[imon]);
	dt = dt * 86400 / (double) NTSTEP;


	for (itts = 1; itts <= NTSTEP; itts++)
	{
		printf("\nSub time step number= %i \n", itts);
		/*-----------------------------------
		 *
		 *     get physical fields and bc's
		 *
		 *----------------------------------*/

		printf("Month %d timestamp Start: %4.2f End: %4.2f\n",cmon,currtime,currtime+dt/31536000);
		currtime += dt / 31536000; //ashao: advance currenttime
		copy_2fix_darray3d(hstart,hend,NZ,NXMEM,NYMEM);
		copy_2fix_darray3d(h,hstart,NZ,NXMEM,NYMEM);
		if (usehindcast) {
			if ( ( cmon >= starthindindex) && ( hindindex <= (numhindmonths-1) )) {

				printf("Reading in UVW from hindcast\n");
				read_uvw(hindindex,"hind");
				read_h(hindindex,hend,"hind");
				hindindex++;
			}
			else {

				printf("Reading in UVW from climatology\n");
				read_uvw(isnxt,"clim");
				read_h(isnxt, hend,"clim");
			}
		}
		else {
			printf("Reading in UVW from climatology\n");
			read_uvw(isnxt,"clim");
			read_h(isnxt, hend,"clim");
		}
		printf("Month- hstart:%d hend:%d\n",ismon,isnxt);

		/*-----------------------------------
		 *
		 *     integrate 1 time step
		 *
		 *----------------------------------*/

		printf("step fields - day = %g\n\n", day);
		step_fields(iyear, itts, imon, iterno); // modified ashao


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
		calculate_tracer_average();


			/*-----------------------------------
			 *
			 *     write tracer averages and reset to 0
			 *
			 *----------------------------------*/
			printf("Writing month %i variables out to netCDF\n\n", cmon);
			write_variables_to_netcdf();
		/*-------------------------------------------------------------------*
		 *
		 *     end integration loop
		 *
		 *-------------------------------------------------------------------*/

	} /* end itts loop over NTSTEP */
} /* end while */

//BX-a

write_restart_netcdf();

printf("\n programme termine normallement. \n");
return (0);
}
