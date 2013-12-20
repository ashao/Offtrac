/*
 * ********+*********+*********+*********+*********+*********+*********+*******
 * Offtrac Version 2.0
 * By Andrew Shao
 * 9 April 2013
 * FILE DESCRIPTION:
 * Contains the routines necessary to initialize Offtrac including parsing the
 * runtime-option namefile and reading in ocean grid information
 * ********+*********+*********+*********+*********+*********+*********+*******
*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "initialize_offtrac.h"
#include "vartypes.h"
#include "read_netcdf.h"
#include "user_options.h"

extern ocean_geometry_t ocean_geometry;
extern ocean_3D_t ocean_3D;
extern ocean_2D_t ocean_2D;
extern io_t io;
extern iterinfo_t iterinfo;
extern filenames_t filenames;

// USER CHANGEABLE
#define NUMOPTIONS 8 // Number of entries in the namelist

void read_parse_runtime( )
{
	int nlines = 0;
	int lineidx;
	int optidx;
	char namefile[]="input.nml";
	char line[MAXLEN];
	char parsename[MAXLEN];
	char parseval[MAXLEN];

	FILE *namefile_fid;
	char runtime_options[NUMOPTIONS][MAXLEN];

	printf("Parsing namefile\n");

	strcpy(runtime_options[0],"inpath");
	strcpy(runtime_options[1],"outpath");
	strcpy(runtime_options[2],"outfile");
	strcpy(runtime_options[3],"startmonth");
	strcpy(runtime_options[4],"startyear");
	strcpy(runtime_options[5],"numiter");
	strcpy(runtime_options[6],"restartflag");
	strcpy(runtime_options[7],"hindcastflag");
	// USER: Add additional runtime options here

	namefile_fid = fopen(namefile,"rt");

	// Parse the namefile. Probably a more efficient way, but it shouldn't add
	// much to the runtime

	// Calculate total number of lines
	while(fgets(line, MAXLEN , namefile_fid) != NULL) nlines++;
	if (nlines != NUMOPTIONS)
		printf("WARNING: Number of lines in namefile does not equal number"
				" of options\n");
	rewind(namefile_fid);

	for (optidx=0;optidx<NUMOPTIONS;optidx++) {
		for (lineidx=0;lineidx<nlines;lineidx++) {

			fgets(line, MAXLEN , namefile_fid);
			sscanf(line,"%s,%s",parsename,parseval);

			if (strcasecmp(parsename,runtime_options[optidx])) {
				printf("\tRuntime option: %s = %s\n",parsename,parseval);
				switch (optidx) {
				case 0:
					strcpy(io.inpath,parseval);
					break;
				case 1:
					strcpy(io.outpath,parseval);
					break;
				case 2:
					strcpy(io.outfile,parseval);
					break;
				case 3:
					iterinfo.startmonth = strtol(parseval,NULL,10);
					break;
				case 4:
					iterinfo.startyear = strtol(parseval,NULL,10);
					break;
				case 5:
					iterinfo.numiter = strtol(parseval,NULL,10);
					break;
				case 6:
					io.restart_flag=parseval[0];
					break;
				case 7:
					io.hindcast_flag=parseval[0];
					break;
				default:
					printf("WARNING: %s not recognized\n",parsename);
				}

				rewind(namefile_fid); // Reset file for reading
				break; // Start parsing other options
			}

		}
	}

	// Update relevant iteration information variables
	iterinfo.currentmonth_idx = iterinfo.startmonth;
	iterinfo.currentmonth_total = iterinfo.startmonth;
	iterinfo.currentyear = iterinfo.startyear;


	fclose(namefile_fid);
}

void initialize_ocean_geometry() {

	char fullpath[MAXLEN];
	strcpy(fullpath,strcat(io.inpath,filenames.metrics));

	// 1D Grid-info
	read_var1d_from1d(fullpath,"lath",NYMEM,ocean_geometry.lath[0]);
	read_var1d_from1d(fullpath,"lonh",NXMEM,ocean_geometry.lonh[0]);

	read_var1d_from1d(fullpath,"latq",NYMEM,ocean_geometry.latq[0]);
	read_var1d_from1d(fullpath,"lonq",NXMEM,ocean_geometry.lonq[0]);

	// 2D-Grid information
	read_var2d_from2d(fullpath,"geolat",ocean_geometry.geolat);
	read_var2d_from2d(fullpath,"geolon",ocean_geometry.geolon);
	read_var2d_from2d(fullpath,"wet", ocean_geometry.wet);
	read_var2d_from2d(fullpath,"Ah",ocean_geometry.Ah);

	// Grid spacing q-points
	read_var2d_from2d(fullpath,"dxq",ocean_geometry.dxq);
	read_var2d_from2d(fullpath,"dyq",ocean_geometry.dyq);

	// Grid spacing h-points
	read_var2d_from2d(fullpath,"dxh",ocean_geometry.dxh);
	read_var2d_from2d(fullpath,"dyh",ocean_geometry.dyh);

	// Grid spacing u-points
	read_var2d_from2d(fullpath,"dxu",ocean_geometry.dxu);
	read_var2d_from2d(fullpath,"dyu",ocean_geometry.dyu);

	// Grid spacing u-points
	read_var2d_from2d(fullpath,"dxv",ocean_geometry.dxv);
	read_var2d_from2d(fullpath,"dyv",ocean_geometry.dyv);


}
