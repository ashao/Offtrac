/*
 * oxygen.h
 *
 *  Created on: Apr 15, 2014
 *      Author: ashao
 */

#include "init.h"
extern int mAGE;
// Output arrays
extern double ***mn_age;
// Working arrays
extern double ***age_init;

/* SUBROUTINE PROTOTYPES */
void allocate_age(  );
void initialize_age(  );
void step_age( double dt );



