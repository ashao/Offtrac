/*
 * gas_exchange.c
 * 17 April 2013
 *
 * MAIN FUNCTION: gas_exchange
 * Calculates the mixed layer concentration of a gas via 1D air-sea gas exchange
 * using an RK2 method with bubble injection all according to
 * Stanley et al. 2009
 *
 * Inputs:
 * 	Csurf: Tracer concentration AFTER transport
 * 	h_start: Beginning of month mixed layer thickness
 * 	h_end: End of month mixed layer thickness
 * 	u10_start: Beginning of month 10m windspped
 * 	u10_end: End of month windspeed
 * 	T_start: Beginning of month SST
 * 	T_end: End of month SST
 * 	S_start: Beginning of month SSS
 * 	S_end: End of month SSS
 * 	dt: Timestep to use for gas exchange integration
 * 	nt: number of timesteps
 *
 *  REFERENCES:
 *
*/
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include "gas_tracer_type.h"
#include "init.h"
#include "seawater_properties.h"
#include "gas_exchange.h"
#define R  8.31 
#define g  9.80665
#define GAMMA_G 0.97*8.6e-7 
//#define GAMMA_G 9.176472738247816e-07 
#define A_C 9.1e-11
#define A_p 2.3e-3

	extern double Temptm[NZ][NXMEM][NYMEM];
	extern double Salttm[NZ][NXMEM][NYMEM];
	extern double windspeed[NXMEM][NYMEM];
	extern double  atmpres[NXMEM][NYMEM];
	extern double  fice[NXMEM][NYMEM];
	void gas_exchange( int i, int j, Gas_transient_t gas, double atmconc,
			double h_start, double h_end, double dt, int nt, double *Csurf, double *Csat) {

		int iter;
		double k1, k2; // Runge Kutta factors
		double avgfact; // Factor calcualted to approximate daily values
		double h, rho; // "Daily" value interpolated from start to end of month
		double Sc, Di, bunsen_sol, F_sol; // Thermodynamic/chemical quantities
		double F_ge, F_c, F_p, F_tot; // Fluxes from each part of the bubble parameterization

		double u10,T,S,P,icefrac;


		T=Temptm[0][i][j];
		S=Salttm[0][i][j];
		icefrac = fice[i][j];

		P=atmpres[i][j];
		u10 = windspeed[i][j];

	//	P=101325;
	//	u10 = 5;

		for(iter=0;iter<nt;iter++) {


			avgfact = (iter)/nt;

			// Calculate the "daily" mixed layer thickness
			h =  avgfact * h_end + (1-avgfact)*h_start;

			// Calculate all other necessary quantities
			Sc = calc_schmidt_number( gas.Sc_coeffs, T);
			Di = calc_Di(T,S,gas.Di_coeffs);
			rho = gsw_rho(S,T,0);

			bunsen_sol = calc_bunsen_sol(
					gas.num_equil_coeffs, gas.equil_coeffs, T, S);
			F_sol = calc_F_sol(
					gas.num_F_sol_coeffs, gas.F_sol_coeffs, T, S);
			*Csat = F_sol*P*atmconc/101325; //Need to convert from Pascal to atmospheres

			F_ge = pow(Sc/660,-0.5)*pow(u10,2)*(*Csat-*Csurf);

		//	F_ge = calc_F_ge(u10, Sc, Csat, Csurf);
			if (u10>2.27) {
				F_c = calc_F_c(u10, P*atmconc, T);
				F_p = calc_F_p(u10, bunsen_sol, Di, atmconc, T, P, rho, *Csurf, *Csat);
			}
			else {
				F_c = 0.0;
				F_p = 0.0;
			}

	//		F_c = 0.0;
	//		F_p = 0.0;

			F_tot = (1-icefrac)*(GAMMA_G*F_ge + A_C*F_c + A_p*F_p);
	/*		if (i==100 && j==100){
			printf("F_tot: %f F_ge: %e F_c: %e F_p: %e k1:%f k2: %f Csurf: %f Csat: %f\n",F_tot/h*dt,GAMMA_G*F_ge,A_C*F_c,A_p*F_p,k1,k2,Csurf,Csat);
			printf("h_start: %f h_end: %f u10: %f T: %f S: %f P: %f fice: %f\n",h_start,h_end,u10,T,S,P,fice);
	}*/
			*Csurf = *Csurf + F_tot/h*dt; // For now just do simple Euler integration later do RK2
	/*
	 * 		// FIRST STEP OF RUNGE-KUTTA
			k1 = 0.5*dt/h * F_tot;

			// SECOND STEP OF RUNGE-KUTTA
			avgfact = (iter+0.5)/nt;

			// Calculate the "daily" varying forcing fields
			h =  avgfact * h_end + (1-avgfact)*h_start;
			T = avgfact * T_end + (1-avgfact)*T_start;
			S = avgfact * S_end + (1-avgfact)*S_start;
			u10 = avgfact * u10_end + (1-avgfact)*u10_start;
			P = avgfact * P_end + (1-avgfact)*P_start;
			P = P * 101325;

			// Calculate all other necessary quantities
			Sc = calc_schmidt_number( gas.Sc_coeffs, T);
			Di = calc_Di(T,S,gas.mol_vol);
			bunsen_sol = calc_bunsen_sol(
					gas.num_bunsen_coeffs, gas.bunsen_coeffs, T, S);
			F_sol = calc_F_sol(
					gas.num_F_sol_coeffs, gas.F_sol_coeffs, T, S);
			Csat = F_sol*P*gas.atmconc;
			rho = gsw_rho(S,T,0);

			F_ge = calc_F_ge(u10, Sc, Csat, Csurf+k1);
			F_c = calc_F_c(u10, P, T);
			F_p = calc_F_p(u10, bunsen_sol, Di, gas.atmconc, T, P*gas.atmconc, rho, Csurf+k1, Csat);
			F_tot = GAMMA_G*F_ge + A_C*F_c + A_p*F_p;
			k2 = dt/h * F_tot;
			printf("Di: %f F_ge: %f F_c: %f F_p: %f k1:%f k2: %f Csat: %f\n",Di,GAMMA_G*F_ge,A_C*F_c,A_p*F_p,k1,k2,Csat);

			Csurf = Csurf + k2;
	*/
		}

	//	printf("End: %f\n",Csurf);
		if (Csurf<0)
			Csurf=0;
	
	}

double calc_F_ge( double u10, double Sc, double *Csat, double *Csurf )
{
// Equation 2 of Stanley et al. 2009

	double F_ge;
	F_ge = pow(Sc/660,-0.5)*pow(u10,2)*(*Csat-*Csurf);
	return F_ge;

}

double calc_F_c(double u10, double P, double T)
{
// Equation 3 of Stanley t al. 2009
	double F_C;
	// Convert T to Kelvin
	T = T + 273.15;

	F_C = pow(u10-2.27,3)*P/(R*T);
	return F_C;

}

double calc_F_p(double u10, double bunsen_sol, double D, double X, double T, double P,
		double rho, double Csurf, double Csat)
{
	double F_p, zbub, TempK;
	TempK = T + 273.15;

	zbub = (0.15*u10-0.55); // Graham et al. 2004
	if (zbub<0)
		zbub=0;
	F_p = pow(u10-2.27,3)*bunsen_sol*pow(D,2/3)*X*P/ (R*TempK)*(1.0 + rho*g*zbub/P - Csurf/Csat);

	if (F_p < 0)
		F_p=0;
	return F_p;
}
double calc_schmidt_number ( double *Sc_coeffs, double T )
{

	double Sc;
	Sc = Sc_coeffs[0]-Sc_coeffs[1]*T+Sc_coeffs[2]*pow(T,2)-Sc_coeffs[3]*pow(T,3);

	return Sc;

}

double calc_Di( double T, double S, double Di_coeffs[2])
{

	/*
	double Di; // Diffusion coefficient
	double mu_sw; // Dynamic viscosity in seawater

	double TempK;
	const double sqrt_q_Mb = sqrt(2.6*18.01528);

	Tempk = T+273.15;
	mu_sw = sw_dyn_visc(T,S); // this is in kg m^-1 s^-1
	mu_sw = mu_sw*10.0; // Converted to g cm^-1 s^-1

	Di = (7.4e-8*sqrt_q_Mb*TempK)/(mu_sw*pow(mol_vol,0.6)); // cm^2 s^-1
	Di = Di*1.0e-4; // m^2 s^-1
	return Di;
	*/

	// Use a different form for diffusion constant

	double Di;
	Di = Di_coeffs[0]*exp(-Di_coeffs[1]/(R*1e-3*(T+273.15)));
	return Di;

}

double calc_bunsen_sol( int num_sol_coeffs, double sol_coeffs[],
		double T, double S)
{

	double bunsen_sol, K_equil; // Bunsen coefficient, partial pressure equilibrium constnat
	const double molar_volume_ideal_gas = 22.71098; // Litres

	double TempK;

	TempK=T+273.15; // Convert to Kelvin
	switch(num_sol_coeffs)
	{
	case 6:
		K_equil = sol_coeffs[0] + sol_coeffs[1]*(100/TempK) +
			sol_coeffs[2]*log(TempK/100) +
			S*(sol_coeffs[3] + sol_coeffs[4]*(TempK/100) +
			sol_coeffs[5]*pow(TempK/100,2));
		K_equil = exp(K_equil); // mol L^-1 atm^-1
		bunsen_sol = K_equil*molar_volume_ideal_gas;
		break;
	default:
		printf("ERROR: Unknown solubility function\n");
		exit(-1);
		break;
	}

	return bunsen_sol; // Should be dimensionless

}

double calc_F_sol( int num_sol_coeffs, double sol_coeffs[],
		double T, double S)
{

	double F_sol;
	T=T+273.15;
	switch(num_sol_coeffs)
	{
	case 7:
		F_sol = sol_coeffs[0] + sol_coeffs[1]*(100/T) + sol_coeffs[2]*log(T/100) +
			sol_coeffs[3]*pow(T/100,2) + S*(sol_coeffs[4] + sol_coeffs[5]*(T/100) +
			sol_coeffs[6]*pow(T/100,2));
		F_sol = exp(F_sol);
		break;
	default:
		printf("ERROR: Unknown solubility function\n");
		exit(-1);
		break;
	}

	return F_sol;

}

void read_gas_exchange_fields(char inpath[200], int monidx, char type_force[4]) {

	int i,j,k;

	char filepath[400];

	strcpy(filepath,inpath);
	strcat(filepath,"auxiliary_input.");
	strcat(filepath,type_force);
	strcat(filepath,".nc");

	printf("Gas exchange fields: %s\n",type_force);
	// Readin variables

	strcpy(filepath,inpath);
	strcat(filepath,"temp.");
	strcat(filepath,type_force);
	strcat(filepath,".nc");
	read_var3d(filepath,"temp",monidx,Temptm, NZ);

	strcpy(filepath,inpath);
	strcat(filepath,"salt.");
	strcat(filepath,type_force);
	strcat(filepath,".nc");
	read_var3d(filepath,"salt",monidx,Salttm, NZ);

	strcpy(filepath,inpath);
	strcat(filepath,"CN.");
	strcat(filepath,type_force);
	strcat(filepath,".nc");
	read_var2d(filepath,"CN",monidx,fice);

	strcpy(filepath,inpath);
	strcat(filepath,"wind.");
	strcat(filepath,type_force);
	strcat(filepath,".nc");
	read_var2d(filepath,"wind",monidx, windspeed);

	strcpy(filepath,inpath);
	strcat(filepath,"p_surf.");
	strcat(filepath,type_force);
	strcat(filepath,".nc");
	read_var2d(filepath,"p_surf",monidx,atmpres);

	// Do the meridional/zonal reentrance
	reentrance_3d(Temptm,NZ);
	reentrance_3d(Salttm,NZ);
	reentrance_2d(atmpres);
	reentrance_2d(windspeed);
	reentrance_2d(fice);

}

double choose_atmconc_val( double sval, double nval, double lat ) {

	double atmconc;

	if (lat>=10.0)
		atmconc = nval;
	if (lat<=-10)
		atmconc = sval;
	// Interpolate across the equator
	if (lat<10.0 && lat>-10.0) 
		atmconc = (nval-sval)/(10+10)*(lat+10)+sval;

	return atmconc;

}
