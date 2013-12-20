/*
 * gas_tracer_type.h
 *
 *  Created on: Apr 17, 2013
 *      Author: ashao
 *
 *  Defines a structure to be used for gas tracers with allocatable values for
 *  transient tracers
 */

// A variable type for gases with constant atmospheric concentration (e.g. O2,
// N2, Ar)
typedef struct {

        int num_Sc_coeffs, num_bunsen_coeffs, num_F_sol_coeffs;
        							 // Number of coeffiicients used in the
        							 // Schmidt number, bunsen solubility, and solubility function
        							 // coefficient equations
        double atmconc; // Northern hemisphere value
        double *Sc_coeffs; // Coefficients for Schmidt number
        double *F_sol_coeffs; // Coefficients for the solubility function of a gas (volumetric)
        double *bunsen_coeffs; // Bunsen solubility coefficients
        double mol_vol; // Molal volume (L/kg)

        /*
        double ***average; // Array to do multi-timestep averaging
        double ***current; // Array to hold the current, working array
		*/
        int tridx; //Index in main tracer array
} Gas_const_t;

// A variable type for gases with transient atmospheric concentration (e.g. CFC-11, CFC-12)
typedef struct {

        int num_Sc_coeffs, num_equil_coeffs, num_F_sol_coeffs;
        							 // Number of coeffiicients used in the
        							 // Schmidt number, bunsen solubility, and solubility function
        							 // coefficient equations

        int ntime; // Number of timestamps
        double *timestamp; // Timestamp for atmospheric values
        double *nval; // Northern hemisphere value
        double *sval; // Southern hemisphere values
        double *Sc_coeffs; // Coefficients for Schmidt number
        double *F_sol_coeffs; // Coefficients for the solubility function of a gas (volumetric)
        double *equil_coeffs; // Partial pressure equilbrium coefficients
        double Di_coeffs[2]; // Coefficients for exponential form

        /*
        double ***average; // Array to do multi-timestep averaging
        double ***current; // Array to hold the current, working array
		*/
        int tridx; //Index in main tracer array
} Gas_transient_t;
