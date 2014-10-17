/*
 * setup_tracers.c
 *
 *  Created on: Oct 16, 2014
 *      Author: ashao
 */
#include "setup_tracers.h"
#include "init.h"

int flags[NOVARS];
int rflags[NOVARS];

double misval = MISVAL;

struct vardesc vars[NOVARS] =
{
		{
				"mn_h", "Layer Thickness", 'h', 'L', 's', "meter", 'd', 9.e-10
		},
		{
				"mn_uh",
				"Zonal Mass Transport",
				'u',
				'L',
				's',
				"meter second^-3",
				'f',
				misval
		}, // 1
		{
				"mn_vh",
				"Meridional Mass Transport",
				'v',
				'L',
				's',
				"meter second^-3",
				'f',
				misval
		}, // 2

		{
				"mn_wd",
				"Diapycnal velocity",
				'u',
				'i',
				's',
				"meter second-1",
				'f',
				misval
		}, // 3
		// Add new tracers below here
		{
				"mn_age",
				"age tracer",
				'h',
				'L',
				's',
				"years",
				'd',
				misval
		}// 4

};

void set_output_flags() {
	// flags[VARINDEX] is an array controlling whether the variable specified
	// by VARINDEX. Defaults are
	// [0]: h
	// [1]: UH
	// [2]: VH
	// [3]: WD
	// [4]: Mean Age
	int i;
	for (i = 0; i <= NOVARS - 1; i++)
		flags[i] = 0;

#ifdef AGE
	flags[4] = 1;
#endif

}

void set_restart_flags() {

	int i;
	 // Default is not to read in for restart
	for (i = 0; i <= NOVARS - 1; i++)  rflags[i] = 0;

	// Do the exceptions
#ifdef AGE
	rflags[4] = 1;
#endif

}

void alloc_fields(void)
{

	int m;

	extern double *var[NOVARS];
	extern double junk[(NZ + 1) * (NXMEM) * (NYMEM)];
	extern long varsize[NOVARS];

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

	var[0] = &mn_h[0][0][0];
	var[1] = &mn_uh[0][0][0];
	var[2] = &mn_vh[0][0][0];
	var[3] = &mn_wd[0][0][0];
#ifdef AGE
	var[4] = &mn_age[0][0][0];
#endif

// end ashao

}

void calculate_tracer_averages() {
	if (flags[0]) add_darray3d(mn_h, h, NZ, NXMEM, NYMEM);
	if (flags[1]) add_darray3d(mn_uh, u, NZ, NXMEM, NYMEM);
	if (flags[2]) add_darray3d(mn_vh, v, NZ, NXMEM, NYMEM);
	if (flags[3]) add_darray3d(mn_wd, wd, NZ+1, NXMEM, NYMEM);
	if (flags[4]) add_darray3d(mn_age, age, NZ, NXMEM, NYMEM);

	/* calculate the mean */
	if (nmn == WRINT && itts == 1)
	{
		printf("***nmn= %i, itts= %i, nmnfirst= %i\n", nmn, itts,
				nmnfirst);
		frac = 1.0 / ((double) nmn * (double) NTSTEP);

		if (flags[0]) mult_darray3d(mn_h, NZ, NXMEM, NYMEM, frac);
		if (flags[1]) mult_darray3d(mn_uh, NZ, NXMEM, NYMEM, frac);
		if (flags[2]) mult_darray3d(mn_vh, NZ, NXMEM, NYMEM, frac);
		if (flags[3]) mult_darray3d(mn_wd, NZ+1, NXMEM, NYMEM, frac);
		if (flags[4]) mult_darray3d(mn_age, NZ, NXMEM, NYMEM, frac);
	}
}

void write_variables_to_netcdf(){
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
#ifdef AGE
	if (flags[9])
		set_darray3d_zero(mn_age, NZ, NXMEM, NYMEM);

#endif
	// begin ashao
	// end ashao

	printf("netcdf record = %d\n", nrec + 1);
	nrec++;
} /*  end if nmn==WRITEINT */
