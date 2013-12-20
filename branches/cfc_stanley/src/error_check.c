/*
 * error_check.c
 * Contains subroutines which will throw errors and cause the program to quit
 * if certain operations fail
 *  Created on: Aug 14, 2013
 *      Author: ashao
 */
#include<stdio.h>
#include<stdlib.h>

void alloc_check_1d( double *ptr, char *varname )
{

	if (ptr == NULL) {
		printf("FATAL: Error allocating %s\n",varname);
		exit(-1);
	}

}

void alloc_check_2d( double **ptr, char *varname )
{

	if (ptr == NULL) {
		printf("FATAL: Error allocating %s\n",varname);
		exit(-1);
	}

}

void alloc_check_3d( double ***ptr, char *varname )
{

	if (ptr == NULL) {
		printf("FATAL: Error allocating %s\n",varname);
		exit(-1);
	}

}
