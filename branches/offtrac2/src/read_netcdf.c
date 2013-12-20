/*
 * ********+*********+*********+*********+*********+*********+*********+*******
 * Offtrac Version 2.0
 * By Andrew Shao
 * 9 April 2013
 * FILE DESCRIPTION:
 * Contains routines to read fields from netCDF fields (2D and 3D)
 * NOTE: ALL THESE ROUTINES TAKE INTO ACCOUNT MEMORY HALO
 * ********+*********+*********+*********+*********+*********+*********+*******
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "netcdf.h"
#include "read_netcdf.h"
#include "user_options.h"

#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(2);}

void read_var1d_from1d( char inpath[MAXLEN], char varname[MAXLEN], int varlen,
		double data[varlen])
{
	// This routine is for reading in 2D fields from a 2D netCDF field
	// For example, geolon, geolat, Ah.

	int i,j,k;
	int err, cdfid, timeid, varid;
	char infile[MAXLEN];
	FILE *file;
	int status;
	int inxt,iprv;
	size_t start[MAX_NC_VARS];
	size_t count[MAX_NC_VARS];
	float tmp1d[varlen];

	start[0] = 0;

	err = open_input_file(inpath,&file,&cdfid,&timeid);
	if ((status = nc_inq_varid(cdfid, varname, &varid)))
		bzero(start, MAX_NC_VARS * sizeof(long));

	count[0] = varlen;

	if ((status = nc_get_vara_float(cdfid,varid,start,count,tmp1d)))
		ERR(status);
		for (i=0;i<varlen;i++)
				data[i+2]= tmp1d[i];

	free(tmp1d);

}

void read_var2d_from2d( char inpath[MAXLEN], char varname[MAXLEN],
		double data[NXMEM][NYMEM])
{
	// This routine is for reading in 2D fields from a 2D netCDF field
	// For example, geolon, geolat, Ah.

	int i,j,k;
	int err, cdfid, timeid, varid;
	char infile[25];
	FILE *file;
	int status;
	int inxt,iprv;
	size_t start[MAX_NC_VARS];
	size_t count[MAX_NC_VARS];
	float tmp2d[NYTOT][NXTOT];

	start[0] = 0;
	start[1] = 0;

	err = open_input_file(inpath,&file,&cdfid,&timeid);
	if ((status = nc_inq_varid(cdfid, varname, &varid)))
		bzero(start, MAX_NC_VARS * sizeof(long));

	count[0] = NXTOT;
	count[1] = NYTOT;

	if ((status = nc_get_vara_float(cdfid,varid,start,count,tmp2d[0])))
		ERR(status);
		for (i=0;i<NXTOT;i++)
			for (j=0;j<NYTOT;j++)
				data[i+2][j+2]= tmp2d[j][i];

	free(tmp2d);

}

void read_var2d_from3d( char inpath[MAXLEN], char varname[MAXLEN], int tidx,
		double data[NXMEM][NYMEM])
{
	// This routine is for reading in 2D fields from a 3D netCDF field
	// For example, January (tidx = 0) SST from 12-month climatology

	int i,j,k;
	int err, cdfid, timeid, varid;
	char infile[25];
	FILE *file;
	int status;
	int inxt,iprv;
	size_t start[MAX_NC_VARS];
	size_t count[MAX_NC_VARS];
	float tmp2d[NYTOT][NXTOT];

	start[0] = tidx;
	start[1] = 0;
	start[2] = 0;

	err = open_input_file(inpath,&file,&cdfid,&timeid);
	if ((status = nc_inq_varid(cdfid, varname, &varid)))
		bzero(start, MAX_NC_VARS * sizeof(long));

	count[0] = 1;
	count[1] = NYTOT;
	count[2] = NXTOT;

	if ((status = nc_get_vara_float(cdfid,varid,start,count,tmp2d[0])))
		ERR(status);

	for (i=0;i<NXTOT;i++)
		for (j=0;j<NYTOT;j++)
			data[i+2][j+2]= tmp2d[j][i];

	free(tmp2d);

}

void read_var3d_from4d( char inpath[MAXLEN], char varname[MAXLEN], int monidx,
		double data[NZ+1][NXMEM][NYMEM])
{
	// This routine is for reading in 3D fields from a 4D netCDF field
	// For example, January UH from 12-month climatology
	int i,j,k;
	int err, cdfid, timeid, varid;
	char infile[25];
	FILE *file;
	int status;
	int inxt,iprv;
	size_t start[MAX_NC_VARS];
	size_t count[MAX_NC_VARS];
	float tmp3d[NZ+1][NYTOT][NXTOT];

	start[0] = monidx;
	start[1] = 0;
	start[2] = 0;
	start[3] = 0;

	err = open_input_file(inpath,&file,&cdfid,&timeid);
	if ((status = nc_inq_varid(cdfid, varname, &varid)))
		bzero(start, MAX_NC_VARS * sizeof(long));

	count[0] = 1;

	// Put in an exception for the WD field since it wants number of interfaces
	if (strcmp(varname,"wd")) count[1]=NZ;
	else count[1] = NZ;
	count[2] = NYTOT;
	count[3] = NXTOT;

	if ((status = nc_get_vara_float(cdfid,varid,start,count,tmp3d[0][0])))
		ERR(status);
	for (k=0;k<count[1];k++)
		for (i=0;i<NXTOT;i++)
			for (j=0;j<NYTOT;j++)
				data[k][i+2][j+2]= tmp3d[k][j][i];

}
