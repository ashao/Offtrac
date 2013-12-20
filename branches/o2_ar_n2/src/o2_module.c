/*
 * o2_module.c
 *
 *  Created on: Apr 17, 2013
 *      Author: ashao
 *      Contains all subroutines related to the oxygen tracer
 */
#include <stdlib.h>
#include "gas_tracer_type.h"
#include "o2_module.h"
#include "alloc.h"
#include "init.h"
#include "gas_exchange.h"
extern Gas_const_t o2_props;

void initialize_o2_properties (  )
{

    int i;
    const int num_Sc_coeffs = 4;
    const int num_F_sol_coeffs = 7;
    const int num_bunsen_coeffs = 6;

    // Store the number of coefficients in the tracer type
    o2_props.num_Sc_coeffs = num_Sc_coeffs;
    o2_props.num_bunsen_coeffs = num_bunsen_coeffs;
    o2_props.num_F_sol_coeffs = num_F_sol_coeffs;

    // Allocate memory for coefficients
    o2_props.Sc_coeffs = (double *)malloc(num_Sc_coeffs * sizeof(double));
    o2_props.F_sol_coeffs = (double *)malloc(num_F_sol_coeffs * sizeof(double));
    o2_props.bunsen_coeffs = (double *)malloc(num_bunsen_coeffs * sizeof(double));


    // Check that they've all been alocated
    alloc_check_1d(o2_props.Sc_coeffs,"N2 Sc_coeffs");
    alloc_check_1d(o2_props.F_sol_coeffs,"N2 F_sol_coeffs");
    alloc_check_1d(o2_props.bunsen_coeffs,"N2 bunsen_coeffs");

//    o2_props.average = alloc3d(NZ,NXMEM,NYMEM);
//    o2_props.current= alloc3d(NZ,NXMEM,NYMEM);

    // Set the actual gas properties

    o2_props.atmconc = 0.209476;

    // Wanninkhof 1992
    o2_props.Sc_coeffs[0] = 1953.4;
    o2_props.Sc_coeffs[1] = 128.00;
    o2_props.Sc_coeffs[2] = 3.9918;
    o2_props.Sc_coeffs[3] = 0.050091;

    // Weiss 1970
    o2_props.bunsen_coeffs[0] = -58.3877;
    o2_props.bunsen_coeffs[1] = 85.8079;
    o2_props.bunsen_coeffs[2] = 23.8439;
    o2_props.bunsen_coeffs[3] = -0.034892;
    o2_props.bunsen_coeffs[4] = 0.015568;
    o2_props.bunsen_coeffs[5] = -0.0019387;


    o2_props.F_sol_coeffs[0] = -173.4292;
    o2_props.F_sol_coeffs[1] = 249.6339;
    o2_props.F_sol_coeffs[2] = 143.3483;
    o2_props.F_sol_coeffs[3] = -21.8492;
    o2_props.F_sol_coeffs[4] = -0.033096;
    o2_props.F_sol_coeffs[5] = 0.014259;
    o2_props.F_sol_coeffs[6] = -0.0017000;

    // Wilke and Chang 1955
    o2_props.mol_vol = 25.6;

    return;
}

void initialize_o2( double ***array ) {
	// For now this function initializes the oxygen tracer to saturation everywhere
	int i, j, k;
	extern double Temptm[NZ][NXMEM][NYMEM];
	extern double Salttm[NZ][NXMEM][NYMEM];
	double F;

	for (k=0;k<NZ;k++)
		for (i=0;i<NXMEM;i++)
			for (j=0;j<NYMEM;j++) {
				F = calc_F_sol( o2_props.num_F_sol_coeffs, o2_props.F_sol_coeffs,
						Temptm[k][i][j], Salttm[k][i][j]);
				array[k][i][j]=F*o2_props.atmconc;
			}

}
