/*
 * n2_module.c
 *
 *  Created on: Apr 17, 2013
 *      Author: ashao
 *      Contains all subroutines related to the nitrogen tracer
 */
#include <stdlib.h>
#include "gas_tracer_type.h"
#include "n2_module.h"
#include "alloc.h"
#include "init.h"
#include "gas_exchange.h"
extern Gas_const_t n2_props;

void initialize_n2_properties ( )
{

    const int num_Sc_coeffs = 4;
    const int num_F_sol_coeffs = 7;
    const int num_bunsen_coeffs = 6;

    // Store the number of coefficients in the tracer type
    n2_props.num_Sc_coeffs = num_Sc_coeffs;
    n2_props.num_bunsen_coeffs = num_bunsen_coeffs;
    n2_props.num_F_sol_coeffs = num_F_sol_coeffs;

    // Allocate memory for coefficients
    n2_props.Sc_coeffs = (double *)malloc(num_Sc_coeffs * sizeof(double));
    n2_props.F_sol_coeffs = (double *)malloc(num_F_sol_coeffs * sizeof(double));
    n2_props.bunsen_coeffs = (double *)malloc(num_bunsen_coeffs * sizeof(double));

    // Check that they've all been alocated
    alloc_check_1d(n2_props.Sc_coeffs,"N2 Sc_coeffs");
    alloc_check_1d(n2_props.F_sol_coeffs,"N2 F_sol_coeffs");
    alloc_check_1d(n2_props.bunsen_coeffs,"N2 bunsen_coeffs");
//    n2_props.average = alloc3d(NZ,NXMEM,NYMEM);
 //   n2_props.current = alloc3d(NZ,NXMEM,NYMEM);
    // Set the actual gas properties

    n2_props.atmconc = 0.780840;

    // Wanninkhof 1992
    n2_props.Sc_coeffs[0] = 2206.1;
    n2_props.Sc_coeffs[1] = 144.86;
    n2_props.Sc_coeffs[2] = 4.5413;
    n2_props.Sc_coeffs[3] = 0.056988;

    // Weiss 1970
    n2_props.F_sol_coeffs[0] = -172.4965;
    n2_props.F_sol_coeffs[1] = 248.4262;
    n2_props.F_sol_coeffs[2] = 143.0738;
    n2_props.F_sol_coeffs[3] = -21.7120;
    n2_props.F_sol_coeffs[4] = -0.049781;
    n2_props.F_sol_coeffs[5] = 0.025018;
    n2_props.F_sol_coeffs[6] = -0.0034861;

    n2_props.bunsen_coeffs[0] = -59.6274;
    n2_props.bunsen_coeffs[1] = 85.7661;
    n2_props.bunsen_coeffs[2] = 24.3696;
    n2_props.bunsen_coeffs[3] = -0.051580;
    n2_props.bunsen_coeffs[4] = -0.026329;
    n2_props.bunsen_coeffs[5] = -0.0037252;

    // Wilke and Chang 1955
    n2_props.mol_vol = 31.2;
}

void initialize_n2( double ***array ) {
	// For now this function initializes the n2 tracer to saturation everywhere
	int i, j, k;
	extern double Temptm[NZ][NXMEM][NYMEM];
	extern double Salttm[NZ][NXMEM][NYMEM];
	double F;

	for (k=0;k<NZ;k++)
		for (i=0;i<NXMEM;i++)
			for (j=0;j<NYMEM;j++) {
				F = calc_F_sol( n2_props.num_F_sol_coeffs, n2_props.F_sol_coeffs,
						Temptm[k][i][j], Salttm[k][i][j]);
				array[k][i][j]=F*n2_props.atmconc;
			}


}
