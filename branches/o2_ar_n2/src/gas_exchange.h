
double gas_exchange_const( Gas_const_t gas, int i, int j, double Csurf, double h_start, double h_end,
		double dt, int nt );

double calc_F_ge( double u10, double Sc, double Csat, double Csurf );

double calc_F_c(double u10, double P, double T);

double calc_F_p(double u10, double bunsen_sol, double D, double X, double T, double P,
                double rho, double Csurf, double Csat);

double calc_schmidt_number ( double *Sc_coeffs, double T );

double calc_Di( double T, double S, double mol_vol);

double calc_bunsen_sol( int num_sol_coeffs, double sol_coeffs[],
                double T, double S);

double calc_F_sol( int num_sol_coeffs, double sol_coeffs[],
                double T, double S);

void read_gas_exchange_fields(char inpath[200], int monidx);
