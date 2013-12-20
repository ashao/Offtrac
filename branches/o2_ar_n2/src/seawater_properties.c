/*
 * seawater_properties.c
 * Contains routines to calculate various properties of seawater
 *  Created on: Aug 14, 2013
 *      Author: ashao
 */
#include<math.h>
#include<stdio.h>
// Dynamic viscosity following equation 22 and 23 of Sharqawy et al. 2010
double sw_dyn_visc(double T, double S)
{

	double mu_w, mu_sw; // Dynamic viscosity of water/seawater
	double A, B; // Temporary values

	A = 1.541 + 1.998e-2*T - 9.52e-5*T*T;
	B = 7.974 - 7.561e-2*T + 4.724e-4*T*T;
	S = S/1000; // (Convert PSU to kg/kg)

	mu_w = 4.2844e-5 + 1/(0.157*pow(T+64.993,2)-91.296);
	mu_sw = mu_w * (1.0 + A*S + B*S*S);

	return mu_sw;

}

/*
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
! density and enthalpy, based on the 48-term expression for density
!--------------------------------------------------------------------------

!==========================================================================
function gsw_rho(sa,ct,p)  !bind(C, name='gsw_rho')
!==========================================================================

!  Calculates in-situ density from Absolute Salinity and Conservative
!  Temperature, using the computationally-efficient 48-term expression for
!  density in terms of SA, CT and p (McDougall et al., 2011).
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature                        [deg C]
! p      : sea pressure                                    [dbar]
!
! gsw_rho  : in-situ density (48 term equation)
*/
double gsw_rho(double sa, double ct, double p)
{
	double	v01 =  9.998420897506056e+2, v02 =  2.839940833161907e0,
		v03 = -3.147759265588511e-2, v04 =  1.181805545074306e-3,
		v05 = -6.698001071123802e0,  v06 = -2.986498947203215e-2,
		v07 =  2.327859407479162e-4, v08 = -3.988822378968490e-2,
		v09 =  5.095422573880500e-4, v10 = -1.426984671633621e-5,
		v11 =  1.645039373682922e-7, v12 = -2.233269627352527e-2,
		v13 = -3.436090079851880e-4, v14 =  3.726050720345733e-6,
		v15 = -1.806789763745328e-4, v16 =  6.876837219536232e-7,
		v17 = -3.087032500374211e-7, v18 = -1.988366587925593e-8,
		v19 = -1.061519070296458e-11,v20 =  1.550932729220080e-10,
		v21 =  1.0e0,
		v22 =  2.775927747785646e-3, v23 = -2.349607444135925e-5,
		v24 =  1.119513357486743e-6, v25 =  6.743689325042773e-10,
		v26 = -7.521448093615448e-3, v27 = -2.764306979894411e-5,
		v28 =  1.262937315098546e-7, v29 =  9.527875081696435e-10,
		v30 = -1.811147201949891e-11, v31 = -3.303308871386421e-5,
		v32 =  3.801564588876298e-7, v33 = -7.672876869259043e-9,
		v34 = -4.634182341116144e-11, v35 =  2.681097235569143e-12,
		v36 =  5.419326551148740e-6, v37 = -2.742185394906099e-5,
		v38 = -3.212746477974189e-7, v39 =  3.191413910561627e-9,
		v40 = -1.931012931541776e-12, v41 = -1.105097577149576e-7,
		v42 =  6.211426728363857e-10, v43 = -1.119011592875110e-10,
		v44 = -1.941660213148725e-11, v45 = -1.864826425365600e-14,
		v46 =  1.119522344879478e-14, v47 = -1.200507748551599e-15,
		v48 =  6.057902487546866e-17;
	double	sqrtsa, v_hat_denominator, v_hat_numerator;

	sqrtsa			= sqrt(sa);

	v_hat_denominator	=
			v01 + ct*(v02 + ct*(v03 + v04*ct))
			+ sa*(v05 + ct*(v06 + v07*ct)
			+ sqrtsa*(v08 + ct*(v09 + ct*(v10 + v11*ct))))
			+ p*(v12 + ct*(v13 + v14*ct) + sa*(v15 + v16*ct)
			+ p*(v17 + ct*(v18 + v19*ct) + v20*sa));

	v_hat_numerator		=
			v21 + ct*(v22 + ct*(v23 + ct*(v24 + v25*ct)))
			+ sa*(v26 + ct*(v27 + ct*(v28 + ct*(v29 + v30*ct)))
			+ v36*sa
			+ sqrtsa*(v31 + ct*(v32 + ct*(v33 + ct*(v34+v35*ct)))))
			+ p*(v37 + ct*(v38 + ct*(v39 + v40*ct))
			+ sa*(v41 + v42*ct)
			+ p*(v43 + ct*(v44 + v45*ct + v46*sa)
			+ p*(v47 + v48*ct)));

	return (v_hat_denominator/v_hat_numerator);
}
