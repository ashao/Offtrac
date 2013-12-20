/*
 * ar_module.c
 *
 *  Created on: Apr 17, 2013
 *      Author: ashao
 *      Contains all subroutines related to the argon tracer
 */
#include <stdlib.h>
#include "init.h"
#include "gas_tracer_type.h"
#include "ar_module.h"
#include "alloc.h"
#include "gas_exchange.h"

extern Gas_const_t ar_props;
extern double ***ar;

void initialize_ar_properties ( )
{

    const int num_Sc_coeffs = 4;
    const int num_F_sol_coeffs = 7;
    const int num_bunsen_coeffs = 6;

    // Store the number of coefficients in the tracer type
    ar_props.num_Sc_coeffs = num_Sc_coeffs;
    ar_props.num_bunsen_coeffs = num_bunsen_coeffs;
    ar_props.num_F_sol_coeffs = num_F_sol_coeffs;

    // Allocate memory for coefficients
    ar_props.Sc_coeffs = (double *)malloc(num_Sc_coeffs * sizeof(double));
    ar_props.F_sol_coeffs = (double *)malloc(num_F_sol_coeffs * sizeof(double));
    ar_props.bunsen_coeffs = (double *)malloc(num_bunsen_coeffs * sizeof(double));

    // Check that they've all been alocated
    alloc_check_1d(ar_props.Sc_coeffs,"N2 Sc_coeffs");
    alloc_check_1d(ar_props.F_sol_coeffs,"N2 F_sol_coeffs");
    alloc_check_1d(ar_props.bunsen_coeffs,"N2 bunsen_coeffs");

//    ar_props.average = alloc3d(NZ,NXMEM,NYMEM);
//    ar_props.current = alloc3d(NZ,NXMEM,NYMEM);
    // Set the actual gas properties

    ar_props.atmconc = 0.00934;

    // Wanninkhof 1992
    ar_props.Sc_coeffs[0] = 1909.1;
    ar_props.Sc_coeffs[1] = 125.09;
    ar_props.Sc_coeffs[2] = 3.9012;
    ar_props.Sc_coeffs[3] = 0.048953;

    // Weiss 1970
    ar_props.bunsen_coeffs[0] = -55.6578;
    ar_props.bunsen_coeffs[1] = 82.0262;
    ar_props.bunsen_coeffs[2] = 22.5929;
    ar_props.bunsen_coeffs[3] = -0.036267;
    ar_props.bunsen_coeffs[4] = 0.016241;
    ar_props.bunsen_coeffs[5] = -0.0020114;


    ar_props.F_sol_coeffs[0] = -173.5146;
    ar_props.F_sol_coeffs[1] = 245.4510;
    ar_props.F_sol_coeffs[2] = 141.8222;
    ar_props.F_sol_coeffs[3] = -21.8020;
    ar_props.F_sol_coeffs[4] = -0.034474;
    ar_props.F_sol_coeffs[5] = 0.014934;
    ar_props.F_sol_coeffs[6] = -0.0017729;

    // Based on volume (given) van der waal's radius times avogadro's nubmer
    ar_props.mol_vol = 18.00;
}

void initialize_ar( double ***array ) {
	// For now this function initializes the argon tracer to saturation everywhere
	int i, j, k;
	extern double Temptm[NZ][NXMEM][NYMEM];
	extern double Salttm[NZ][NXMEM][NYMEM];
	double F;

	for (k=0;k<NZ;k++)
		for (i=0;i<NXMEM;i++)
			for (j=0;j<NYMEM;j++) {
				F = calc_F_sol( ar_props.num_F_sol_coeffs, ar_props.F_sol_coeffs,
						Temptm[k][i][j], Salttm[k][i][j]);
				array[k][i][j]=F*ar_props.atmconc;
			}


}



