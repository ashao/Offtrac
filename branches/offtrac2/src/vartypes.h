/*
 * ********+*********+*********+*********+*********+*********+*********+*******
 * By Andrew Shao
 * 9 April 2013
 * Description: This file defines various structures used to store the physical
 * fields derived from HIM/GOLD, ocean geometry, and other model data
 * ********+*********+*********+*********+*********+*********+*********+*******
*/
#include <stdio.h>
#include "user_options.h"



// Ocean Geometry structure
typedef struct {

	// Wet Masks
	double wet[NXMEM][NYMEM];
	int hmask[NXMEM][NYMEM];
	int umask[NXMEM][NYMEM];
	int vmask[NXMEM][NYMEM];

	// Topography
	double D[NXMEM][NYMEM];

	// Grid information
	double geolon[NXMEM][NYMEM];
	double geolat[NXMEM][NYMEM];
	double lonh[NXMEM][NYMEM];
	double lath[NXMEM][NYMEM];
	double lonq[NXMEM][NYMEM];
	double latq[NXMEM][NYMEM];
	double Ah[NXMEM][NYMEM];

	// Grid Spacing q-points
	double dxq[NXMEM][NYMEM], dyq[NXMEM][NYMEM];
	// Grid Spacing h-points
	double dxh[NXMEM][NYMEM], dyh[NXMEM][NYMEM];
	// Grid Spacing u-points
	double dxu[NXMEM][NYMEM], dyu[NXMEM][NYMEM];
	// Grid Spacing v-points
	double dxv[NXMEM][NYMEM], dyv[NXMEM][NYMEM];

} ocean_geometry_t;


// 3-dimensional physical fields
typedef struct {

	// Isopycnal thickness
	double h[NZ][NXMEM][NYMEM]; // Monthly averaged thickness
	double hstart[NZ][NXMEM][NYMEM];// Instantaneous isopycnal thickness,
									// at start of month.
	double hend[NZ][NXMEM][NYMEM];// Instantaneous isopycnal thickness,
										// at start of month.

	// Mass transports UH + UHGM + UHML
	double uh[NZ][NXMEM][NYMEM];
	double vh[NZ][NXMEM][NYMEM];
	double wd[NZ+1][NXMEM][NYMEM]; //

	// Temperature and salinity
	double temp[NZ][NXMEM][NYMEM];
	double salt[NZ][NXMEM][NYMEM];
} ocean_3D_t;


// 2-dimensional physical fields
typedef struct {

	double icefrac[NXMEM][NYMEM];
	double xkw[NXMEM][NYMEM];

} ocean_2D_t;

// I/O
typedef struct {
	char inpath[300];	// path where Offtrac input files are stored
	char outpath[300];	// where to save the output
	char outfile[300];	// Output filename
	char restart_flag; 	// Specify whether the run is a restart ('y') or
						// cold start ('n')
	char hindcast_flag; // Specify whether to use hindcast ('y') or not ('n')

} io_t;

// Iteration information
typedef struct {

	// Set at runtime
	int numiter;		// Total number of iterations
	int numsubsteps;	// Number of substeps to take in between iterations
	int startyear;		// Used primarily for hindcast, but also for
						// climatology when using time-dependent tracers
	int startmonth;		// Used to start Offtrac on a different month

	// Set internally and updated each timestep
	int iteration;			// Current iteration number
	int currentmonth_total;	// Current month total (e.g. Month 24)
	int currentmonth_idx;	// Month index (e.g. 0 = Jan, ranges from 0-11)
	int currentyear;		// Current year
	int modelyearday;		// Timestamp of model output in yeardays
	double modelyearfrac; 	// Timestamp as in yearfraction (e.g. 1950.5)


	int lastmonth;
	int nextmonth;

} iterinfo_t;


// Default filenames
typedef struct {
	char namefile[MAXLEN];
	char metrics[MAXLEN];
} filenames_t;
