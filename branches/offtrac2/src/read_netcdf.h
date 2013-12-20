#include "user_options.h"
void read_var1d_from1d( char inpath[MAXLEN], char varname[MAXLEN], int varlen,
		double data[varlen]);
void read_var2d_from2d( char inpath[MAXLEN], char varname[MAXLEN],
		double data[NXMEM][NYMEM]);
void read_var2d_from3d( char inpath[MAXLEN], char varname[MAXLEN], int tidx,
		double data[NXMEM][NYMEM]);
void read_var3d_from4d( char inpath[MAXLEN], char varname[MAXLEN], int monidx,
		double data[NZ][NXMEM][NYMEM]);
