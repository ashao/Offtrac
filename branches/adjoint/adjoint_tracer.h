/*
 * ttd_bp.h
 *
 *  Created on: May 1, 2014
 *      Author: ashao
 */
extern int madjoint;
extern double ***mn_adjoint;
extern double *mn_adjointinv;
extern double *adjointinv;
void allocate_adjoint( );
void initialize_adjoint( );
void apply_adjoint_source( );
void surface_adjoint( );
