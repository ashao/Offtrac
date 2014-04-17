/*
 * indian_ocean.c
 *  Does a TTD via BIR but only in the Indian Ocean wher mixed
 *  layer depths exceed 200m
 *  Created on: Apr 11, 2014
 *      Author: ashao
 */
#include "init.h"
#include "alloc_arrays.h"

double ***mn_vent, ***vent;
extern double ****tr;

void initialize_ventilation_array( )
{

    mn_vent = alloc3d(NZ,NXMEM,NYMEM);
    vent = alloc3d(NZ,NXMEM,NYMEM);

}

void vent_tracer()
{



}
