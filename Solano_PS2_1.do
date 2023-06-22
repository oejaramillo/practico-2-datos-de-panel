// RGS - UTDT

// https://www.statalist.org/forums/forum/general-stata-discussion/general/1453333-how-to-calculate-the-inverse-mill-ratio-imr-as-the-stata-does
// https://www.stata.com/statalist/archive/2005-04/msg00109.html
// https://www.statalist.org/forums/forum/general-stata-discussion/general/1395655-right-method-for-unbalanced-panel-data
// chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://www.statisticslab.org/resources/Bendig_Hoke_2022_Statistitcs_Lab_Muenster_Presentation_Heckman_two_stage_estimation.pdf
// https://www.statalist.org/forums/forum/general-stata-discussion/general/1574497-heckman-model-for-selection-bias-with-panel-data
// PS4 ex 2
// how to handle FE/FE in second step?

cls
clear all
set more off

// Seed number (last digits of DNI)
set seed 5026

// Sizes of entity and time dimensions
global N = 20
global T = 2
global NT = $N * $T

// Model selector (A = 1, B = 2, C = 3)
global model = 1

// Coefficients in structural model
global beta = 1
global gamma1 = 1
global gamma2 = 1

// Number of simulations
global simulations = 1

// Matrices to store coefficients
matrix beta_hat_stata = J($simulations , 1, .)
matrix gamma1_hat_stata = J($simulations , 1, .)
matrix gamma2_hat_stata = J($simulations , 1, .)

// Bootstrapping program
program onebootrep, rclass
// 	args s_resid y_resid
// 	drop _all
	preserve
	bsample
	mkmat s_resid, matrix(s_resid)
	mkmat y_resid, matrix(y_resid)
	restore
	svmat s_resid
	svmat y_resid

	list s_resid1 y_resid1, sepby(j)
	
	gen s2_jt = gamma1_hat * w_jt + gamma2_hat * z_jt + alpha_j + s_resid1
	gen s2 = 1 if s2_jt > 0
	replace s2 = 0 if missing(s2)		

	if $model == 1 {
		probit s2 w_jt z_jt
	}
	if $model == 2 | $model == 3 {
		probit s2 w_jt z_jt z_j
	}
	return scalar gamma1_hat2 = e(b)[1, 1]
	return scalar gamma2_hat2 = e(b)[1, 2]	

	predict s2_hat, xb
	gen imr2 = normalden(s2_hat)/normal(s2_hat) if s == 1

	gen y2_jt = beta_hat * x_jt + c_j + y_resid1
	gen y2 = y2_jt * s2

	if $model == 1 {
		reg y2 x_jt i.t#c.imr2 if s2 == 1
	}
	if $model == 2 | $model == 3 {
		reg y2 x_jt x_j i.t#c.imr2 if s2 == 1
	}	
	return scalar beta_hat2 = e(b)[1, 1]

	drop s_resid1 y_resid1 s2_jt s2 s2_hat imr2 y2_jt y2
end

// For loop
// quietly
forvalues i = 1 (1) $simulations {
//	To know where we are in the loop
	if `i' == 1 | floor((`i')/50) == (`i')/50 {
		noisily display "Loop `i' out of $simulations ($S_TIME)" 
		}
	
//	Sets arrays size
	set obs $NT
	
//	Generates indexes, one for entity one for time
	egen j = seq(), f(1) t($N) b($T)
	egen t = seq(), f(1) t($T)
	
//	Sets panel indexes
	xtset j t
	
//	Generates random variable 
	gen psi1 = rnormal(0, sqrt(1))
	
//	Generates variables and error for the selection equation
	gen w_jt = rnormal(0, sqrt(1))
	gen z_jt = rnormal(0, sqrt(1))
	
//	Generates regressor for the equation of interest
	gen x_jt = w_jt	
	
//	Generates error for the selection equation	
	gen epsilon_jt = rnormal(0, sqrt(1))	
		
//	Generates error for the equation of interest
	gen u_jt = 0.6 * epsilon_jt + 0.8 * psi1
	
//	Generates averages of x and z variables
	by id, sort: egen z_j = mean(z_jt)
	by id, sort: egen x_j = mean(x_jt)
	
//	If statement that creates fixed effects according to the chosen model
	gen alpha_j = .
	gen c_j = .
	
	forvalues k = 1/$N {
	    scalar psi2 = rnormal(0, sqrt(1))
		scalar psi3 = rnormal(0, sqrt(1))
		scalar psi4 = rnormal(0, sqrt(1))
		if $model == 1 {
			replace alpha_j = psi2 + psi4 if j == `k'
			replace c_j = psi3 + psi4 if j == `k'
		}
		if $model == 2 {
		    replace alpha_j = psi2 + z_j if j == `k'
			replace c_j = psi3 + x_j if j == `k'
		}
		if $model == 3 {
		    replace alpha_j = psi2 + z_j + psi4 if j == `k'
			replace c_j = psi3 + x_j + psi4 if j == `k'
		}
	}	
		
//	Generates selection variable
	gen s_jt = $gamma1 * w_jt + $gamma2 * z_jt + alpha_j + epsilon_jt

//	Generates selector
	gen s = 1 if s_jt > 0	
	replace s = 0 if missing(s)

//	Generates time-dummies
	forvalues l = 1/$T {
		gen d`l' = 1 if t == `l'
		replace d`l' = 0 if missing(d`l')
	}
	
//	Generates dependant variable for all sample
	gen y_jt = $beta * x_jt + c_j + u_jt
					
//	Generates dependant variable after the incidental truncation
	gen y = y_jt * s
	
// //	Mine
// 	if $model == 1 {
// 		forvalues m = 1/$T {
// 			probit s w_jt z_jt if t == `m'
// 			predict s_hat, xb
// 			gen imr`m' = normalden(s_hat)/normal(s_hat) if t == `m' & s == 1
// 			replace imr`m' = 0 if missing(imr`m')
// 			drop s_hat
// 		}
// 	}
// 	if $model == 2 | $model == 3 {
// 		forvalues m = 1/$T {
// 			probit s w_jt z_jt z_j if t == `m'
// 			predict s_hat, xb
// 			gen imr`m' = normalden(s_hat)/normal(s_hat) if t == `m' & s == 1
// 			replace imr`m' = 0 if missing(imr`m')
// 			drop s_hat
// 		}
// 	}
//
// 	if $model == 1 {
// 		reg y x_jt imr* if s == 1
// 	}
// 	if $model == 2 | $model == 3 {
// 		reg y x_jt x_j imr* if s == 1		
// 	}
//
//	
// // Anastasia
// 	gen lambda = .
//
// 	local i = 1
// 	while `i' <= $T {
// 		probit s w_jt z_jt if t == `i'
// 		predict xb, xb
// 		replace lambda = normalden(xb)/normal(xb) if s == 1 & t ==`i'
// 		drop xb
// 		local i = `i' + 1
// 	}
//	
// 	local i = 2
// 	while `i' <= $T {
// 		gen lam`i'= d`i'*lambda
// 		local i = `i' + 1
// 	}
//	
// 	reg y x_jt lam* if s == 1
	
// // Romi PUNTO 1
// 	if $model == 1 {
// 		probit s w_jt z_jt
// 	}
// 	if $model == 2 | $model == 3 {
// 		probit s w_jt z_jt z_j
// 	}
//	
// 	matrix gamma1_hat_stata[`i', 1] = _b[w_jt]
// 	matrix gamma2_hat_stata[`i', 1] = _b[z_jt]	
//	
// 	predict s_hat, xb
// 	gen imr = normalden(s_hat)/normal(s_hat) if s == 1
// 	replace imr = 0 if missing(imr)
//		
// 	if $model == 1 {
// 		reg y x_jt i.t#c.imr if s == 1
// 	}
// 	if $model == 2 | $model == 3 {
// 		reg y x_jt x_j i.t#c.imr if s == 1
// 	}
//	
// 	matrix beta_hat_stata[`i', 1] = _b[x_jt]
	
	
//	PUNTO 2 (simulations = 1)
	if $model == 1 {
		probit s w_jt z_jt
	}
	if $model == 2 | $model == 3 {
		probit s w_jt z_jt z_j
	}
	scalar gamma1_hat = e(b)[1, 1]
	scalar gamma2_hat = e(b)[1, 2]	
	predict s_hat, xb
	gen s_resid = s - s_hat
	gen imr = normalden(s_hat)/normal(s_hat) if s == 1
	replace imr = 0 if missing(imr)
	if $model == 1 {
		reg y x_jt i.t#c.imr if s == 1
	}
	if $model == 2 | $model == 3 {
		reg y x_jt x_j i.t#c.imr if s == 1
	}	
	scalar beta_hat = e(b)[1, 1]
	predict y_hat, xb
	gen y_resid = y - y_hat

	
** Ahora hacemos 50 replicaciones de bootstrap
// 	preserve

	global reps = 10
	simulate gamma1_hat2 = r(gamma1_hat2) gamma2_hat2 = r(gamma2_hat2) beta_hat2 = r(beta_hat2), reps($reps) nodots: onebootrep
	
// ** Computo de los errores estandar de bootstrap 
	
	list beta_hat2
	
	
	
//	BOOTSTRAP CONFIDENCE INTERVALS
//	SIGNIFICATIVITY (?) INTEVARLS FROM ASYMPTOTIC VAR/COVAR MATRIX
	
	
	
// 	Stata things
// 	drop j t psi1 w_jt z_jt x_jt epsilon_jt u_jt z_j x_j alpha_j c_j s_jt s d* y_jt y imr*
// 	drop j t psi1 w_jt z_jt x_jt epsilon_jt u_jt z_j x_j alpha_j c_j s_jt s d* y_jt y lam*
// 	drop j t psi1 w_jt z_jt x_jt epsilon_jt u_jt z_j x_j alpha_j c_j s_jt s d* y_jt y s_hat imr
}

// matlist beta_hat_stata
// matlist gamma1_hat_stata
// matlist gamma2_hat_stata
// matlist fitted_vals

// // Python magic (finally, Home)
// python
// # Imports libraries
// import numpy as np
// import pandas as pd
// from sfi import Matrix
// from sklearn.metrics import mean_squared_error
//
// # Get parameters from Stata, which are stored in global variables
// beta = $beta
// gamma1 = $gamma1
// gamma2 = $gamma2
//
// # Gets data from Stata, which is stored in matrixes
// beta_hat_data = Matrix.get('beta_hat_stata')
// gamma1_hat_data = Matrix.get('gamma1_hat_stata')
// gamma2_hat_data = Matrix.get('gamma2_hat_stata') 
//
// # Mean values
// print(np.mean(beta_hat_data))
// print(np.mean(gamma1_hat_data))
// print(np.mean(gamma2_hat_data))
//
// end
