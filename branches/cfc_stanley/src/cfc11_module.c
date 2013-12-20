/*
 * cfc11_module.c
 *
 *  Created on: Apr 17, 2013
 *      Author: ashao
 *      Contains all subroutines related to the cfc11 tracer
 */
#include <stdlib.h>
#include "init.h"
#include "gas_tracer_type.h"
#include "cfc11_module.h"
#include "alloc.h"
#include "read.h"
#include <string.h>
#include <netcdf.h>
#include <stdio.h>

extern Gas_transient_t cfc11_props;

extern double ***cfc11;

void initialize_cfc11_properties( )
{

    const int num_Sc_coeffs = 4;
    const int num_F_sol_coeffs = 7;
    const int num_equil_coeffs = 6;

    printf("Setting all the chemical/physical properties for CFC-11\n");
    // Store the number of coefficients in the tracer type
    cfc11_props.num_Sc_coeffs = num_Sc_coeffs;
    cfc11_props.num_equil_coeffs = num_equil_coeffs;
    cfc11_props.num_F_sol_coeffs = num_F_sol_coeffs;

    // Allocate memory for coefficients
    cfc11_props.Sc_coeffs = (double *)malloc(num_Sc_coeffs * sizeof(double));
    cfc11_props.F_sol_coeffs = (double *)malloc(num_F_sol_coeffs * sizeof(double));
    cfc11_props.equil_coeffs = (double *)malloc(num_equil_coeffs * sizeof(double));

    // Check that they've all been alocated
    alloc_check_1d(cfc11_props.Sc_coeffs,"N2 Sc_coeffs");
    alloc_check_1d(cfc11_props.F_sol_coeffs,"N2 F_sol_coeffs");
    alloc_check_1d(cfc11_props.equil_coeffs,"N2 equil_coeffs");

    // Set the actual gas properties

    // Zheng et al.
    cfc11_props.Sc_coeffs[0] = 3501.8;
    cfc11_props.Sc_coeffs[1] = 210.31;
    cfc11_props.Sc_coeffs[2] = 6.1851;
    cfc11_props.Sc_coeffs[3] = 0.07513;

    // Warner and Weiss 1985
    cfc11_props.equil_coeffs[0] = -134.1536;
    cfc11_props.equil_coeffs[1] = 203.2156;
    cfc11_props.equil_coeffs[2] = 56.2320;
    cfc11_props.equil_coeffs[3] = -0.144449;
    cfc11_props.equil_coeffs[4] = 0.092952;
    cfc11_props.equil_coeffs[5] = -0.0159977;


    cfc11_props.F_sol_coeffs[0] = -229.9261;
    cfc11_props.F_sol_coeffs[1] = 319.6552;
    cfc11_props.F_sol_coeffs[2] = 119.4471;
    cfc11_props.F_sol_coeffs[3] = -1.39165;
    cfc11_props.F_sol_coeffs[4] = -0.142382;
    cfc11_props.F_sol_coeffs[5] = 0.091459;
    cfc11_props.F_sol_coeffs[6] = -0.0157274;

    cfc11_props.Di_coefs[0]=0.015; // Zheng et al.
    cfc11_props.Di_coefs[1]=18.1;


}

void read_cfc11_bc()
{

	extern char directory[75];
	int err, cdfid, timeid;
	char inpath[200];
	FILE *file;
	int status;
	int tempid,varid;
	char varname_time[20] = "Year\0";
	char varname_nval[20]="CFC11NH\0";
	char varname_sval[20]="CFC11SH\0";
	size_t ntime;
	size_t start[2]={0, 0};
	size_t count[2];

	printf("Reading CFC-11 boundary condition\n");

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
	cfc11_props.ntime = ntime;
	count[0] = ntime;
	count[1] = ntime;

	// Allocate the arrays to hold atmospheric values and timestamps
    cfc11_props.timestamp = (double *)malloc(ntime * sizeof(double));
    cfc11_props.nval = (double *)malloc(ntime * sizeof(double));
    cfc11_props.sval = (double *)malloc(ntime * sizeof(double));
    // Check that they've all been alocated
	
    alloc_check_1d(cfc11_props.timestamp,"CFC-11 timestamp");
    alloc_check_1d(cfc11_props.nval,"CFC-11 NVAL");
    alloc_check_1d(cfc11_props.sval,"CFC-11 SVAL");
	
	// Read in timestamp of atmospheric values

	status = nc_inq_varid(cdfid, varname_time, &varid);
	status = nc_get_vara_double(cdfid,varid,start,count,cfc11_props.timestamp);
	// Read in northern hemisphere values

	status = nc_inq_varid(cdfid, varname_nval, &varid);
	status = nc_get_vara_double(cdfid,varid,start,count,cfc11_props.nval);
	// Read in southern hemisphere values

	status = nc_inq_varid(cdfid, varname_sval, &varid);
	status = nc_get_vara_double(cdfid,varid,start,count,cfc11_props.sval);
	close_file(&cdfid,&file);

}



