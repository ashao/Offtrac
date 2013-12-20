/*
 * ar_module.c
 *
 *  Created on: Apr 17, 2013
 *      Author: ashao
 *      Contains all subroutines related to the argon tracer
 */
#include <stdlib.h>
#include "init.h"
#include "gas_tracer_type.h"
#include "cfc12_module.h"
#include "alloc.h"
#include "read.h"
#include <string.h>
#include <stdio.h>
#include <netcdf.h>
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(2);}
extern Gas_transient_t cfc12_props;

extern double ***cfc12;

void initialize_cfc12_properties( )
{

    const int num_Sc_coeffs = 4;
    const int num_F_sol_coeffs = 7;
    const int num_equil_coeffs = 6;

    printf("Setting all the chemical/physical properties for CFC-12\n");
    // Store the number of coefficients in the tracer type
    cfc12_props.num_Sc_coeffs = num_Sc_coeffs;
    cfc12_props.num_equil_coeffs = num_equil_coeffs;
    cfc12_props.num_F_sol_coeffs = num_F_sol_coeffs;

    // Allocate memory for coefficients
    cfc12_props.Sc_coeffs = (double *)malloc(num_Sc_coeffs * sizeof(double));
    cfc12_props.F_sol_coeffs = (double *)malloc(num_F_sol_coeffs * sizeof(double));
    cfc12_props.equil_coeffs = (double *)malloc(num_equil_coeffs * sizeof(double));

    // Check that they've all been alocated
    alloc_check_1d(cfc12_props.Sc_coeffs,"N2 Sc_coeffs");
    alloc_check_1d(cfc12_props.F_sol_coeffs,"N2 F_sol_coeffs");
    alloc_check_1d(cfc12_props.equil_coeffs,"N2 equil_coeffs");

    // Set the actual gas properties

    // Zheng et al.
    cfc12_props.Sc_coeffs[0] = 3845.4;
    cfc12_props.Sc_coeffs[1] = 228.95;
    cfc12_props.Sc_coeffs[2] = 6.1908;
    cfc12_props.Sc_coeffs[3] = 0.067430;

    // Warner and Weiss 1985
    cfc12_props.equil_coeffs[0] = -122.3246;
    cfc12_props.equil_coeffs[1] = 182.5306;
    cfc12_props.equil_coeffs[2] = 50.5898;
    cfc12_props.equil_coeffs[3] = -0.145633;
    cfc12_props.equil_coeffs[4] = 0.092509;
    cfc12_props.equil_coeffs[5] = -0.0156627;

    cfc12_props.F_sol_coeffs[0] = -218.0971;
    cfc12_props.F_sol_coeffs[1] = 298.9702;
    cfc12_props.F_sol_coeffs[2] = 113.8049;
    cfc12_props.F_sol_coeffs[3] = -1.39165;
    cfc12_props.F_sol_coeffs[4] = -0.143566;
    cfc12_props.F_sol_coeffs[5] = 0.091015;
    cfc12_props.F_sol_coeffs[6] = -0.0153924;

    cfc12_props.Di_coefs[0]=0.036; // Zheng et al.
    cfc12_props.Di_coefs[1]=20.1;
}

void read_cfc12_bc()
{


	extern char directory[75];
	int err, cdfid, timeid;
	char inpath[200];
	FILE *file;
	int status;
	int tempid,varid;
	char varname_time[20] = "Year\0";
	char varname_nval[20]="CFC12NH\0";
	char varname_sval[20]="CFC12SH\0";
	size_t ntime;
	size_t start[2]={0, 0};
	size_t count[2];

	printf("Reading CFC-12 boundary condition\n");

	// Read in the timevarying boundary of atmospheric CFC/SF6 concentration
	strcpy(inpath,directory);
	strcat(inpath,"cfc_sf6_bc.nc");

	err = open_input_file(inpath,&file,&cdfid,&timeid);
	if (err != 0) {
		strcat(inpath, ".cdf");
		err = open_input_file(inpath,&file,&cdfid,&timeid);
		if (err != 0) {
			printf("Unable to find cfc_sf6_bc file.\n");
			exit(-73);
		}
	}
	// Get length of the time dimension
	status = nc_inq_dimid(cdfid, varname_time, &timeid);
	if (status != NC_NOERR) handle_error(status);
	status = nc_inq_dimlen(cdfid, timeid, &ntime);
	if (status != NC_NOERR) handle_error(status);
	    cfc12_props.ntime = ntime;
	count[0] = ntime;
	count[1] = ntime;

	// Allocate the arrays to hold atmospheric values and timestamps
    cfc12_props.timestamp = (double *)malloc(ntime * sizeof(double));
    cfc12_props.nval = (double *)malloc(ntime * sizeof(double));
    cfc12_props.sval = (double *)malloc(ntime * sizeof(double));
    // Check that they've all been alocated

    alloc_check_1d(cfc12_props.timestamp,"CFC-12 timestamp");
    alloc_check_1d(cfc12_props.nval,"CFC-12 NVAL");
    alloc_check_1d(cfc12_props.sval,"CFC-12 SVAL");

	// Read in timestamp of atmospheric values

	status = nc_inq_varid(cdfid, varname_time, &varid);
	status = nc_get_vara_double(cdfid,varid,start,count,cfc12_props.timestamp);
	// Read in northern hemisphere values

	status = nc_inq_varid(cdfid, varname_nval, &varid);
	status = nc_get_vara_double(cdfid,varid,start,count,cfc12_props.nval);
	// Read in southern hemisphere values

	status = nc_inq_varid(cdfid, varname_sval, &varid);
	status = nc_get_vara_double(cdfid,varid,start,count,cfc12_props.sval);
	close_file(&cdfid,&file);
}

