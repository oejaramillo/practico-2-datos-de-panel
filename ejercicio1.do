clear all
set more off
set seed 1263


*Generamos panel 
local rep = 1000
global T = 2				// cantidad de periodos 
global N = 20				// cantidad de individuos 
global NT = $N * $T	

*Parámetros
scalar beta = 1 
scalar gamma1 = 1
scalar gamma2 = 1

* Creamos variables para guardar las estimaciones
set obs `rep' 

*Modelo A
gen modelo_a_b = .
gen modelo_a_g1 = .
gen modelo_a_g2 = .
gen modelo_a_rmse_sia = .
gen modelo_a_rmse_ya = .

*Modelo B
gen modelo_b_b = .
gen modelo_b_g1 = .
gen modelo_b_g2 = .
gen modelo_a_rmse_sib = .
gen modelo_a_rmse_yb = .

*Modelo C
gen modelo_c_b = .
gen modelo_c_g1 = .
gen modelo_c_g2 = .
gen modelo_a_rmse_sic = .
gen modelo_a_rmse_yc = .


preserve

clear

* Simulaciones
forvalues i = 1(1)`rep' {
	
	* Generamos el panel
	set obs $NT
	
	egen id = seq(), f(1) t($N) b($T)
	egen tt = seq(), f(1) t($T)
	xtset id tt

	*Variables necesarias
	gen xjt = rnormal(0,1)
	gen zjt = rnormal(0,1)
	gen psi1 = rnormal(0,1)

	gen psi2 = .
	gen psi3 = .
	gen psi4 = .
	foreach var of varlist psi2 psi3 psi4 {
		forvalues r = 1/$N {
			gen psi = rnormal(0,1)
			replace `var' = psi if id == `r'
			drop psi
		}
	}
	gen epsilonjt = rnormal(0,1)
	gen ujt = 0.6*epsilonjt + 0.8*psi1

	bysort id: egen zj = mean(zjt)
	bysort id: egen xj = mean(xjt)

	* Vamos a generar variables diferentes para cada modelo
	*Modelo A
	gen alfaja = psi2 + psi4
	gen cja = psi3 + psi4

	*Modelo B
	gen alfajb = psi2 + zj
	gen cjb = psi3 + xj

	*Modelo C
	gen alfajc = psi2 + zj + psi4
	gen cjc = psi3 + xj + psi4

	* Ecuaciones de selección para cada modelo
	gen sjta = gamma1*xjt + gamma2*zjt + alfaja + epsilonjt
	gen sjtb = gamma1*xjt + gamma2*zjt + alfajb + epsilonjt
	gen sjtc = gamma1*xjt + gamma2*zjt + alfajc + epsilonjt

	* Indicadoras para cada modelo
	gen sia = (sjta > 0)
	gen sib = (sjtb > 0)
	gen sic = (sjtc > 0)

	* Ecuaciones de interés para cada modelo
	gen yjta = beta*xjt + cja + ujt
	gen yjtb = beta*xjt + cjb + ujt
	gen yjtc = beta*xjt + cjc + ujt
	 
	* Variables perdidas 
	gen ya = yjta*sia
	gen yb = yjtb*sib
	gen yc = yjtc*sic

	* Estimación de Wooldridge
	* Probit Modelo A
	probit sia xjt zjt
		local gamma1a_est = _b[xjt]
		local gamma2a_est = _b[zjt]
		predict sia_est, xb
		gen sia_resid = sia - sia_est
	gen lambdait_a = normalden(sia_est)/normal(sia_est) if sia == 1
	replace lambdait_a = 0 if missing(lambdait_a)

	* Probit Modelo B
	probit sib xjt zjt zj
		local gamma1b_est = _b[xjt]
		local gamma2b_est = _b[zjt]
		predict sib_est, xb
		gen sib_resid = sib - sib_est
	gen lambdait_b = normalden(sib_est)/normal(sib_est) if sib == 1
	replace lambdait_b = 0 if missing(lambdait_b)

	* Probit Modelo C
	probit sic xjt zjt zj
		local gamma1c_est = _b[xjt]
		local gamma2c_est = _b[zjt]
		predict sic_est, xb
		gen sic_resid = sic - sic_est
	gen lambdait_c = normalden(sic_est)/normal(sic_est) if sic == 1
	replace lambdait_c = 0 if missing(lambdait_c)

	* Estimación de Woolridge 
	* Modelo A
	reg ya xjt i.tt#c.lambdait_a if sia == 1
		local betaa_est = _b[xjt]
		predict ya_resid, resid

	* Modelo B
	reg yb xjt xj i.tt#c.lambdait_b if sib == 1
		local betab_est = _b[xjt]
		predict yb_resid, resid
	* Modelo C
	reg yc xjt xj i.tt#c.lambdait_c if sic == 1
		local betac_est = _b[xjt]
		predict yc_resid, resid
		
	* Calculo de los residuos al cuadrado para el RMSE
	gen sia2_resid = sia_resid*sia_resid if sia == 1
		sum sia2_resid
		local rmse_a = sqrt(r(mean))
	gen sib2_resid = sib_resid*sib_resid if sib == 1
		sum sib2_resid
		local rmse_b = sqrt(r(mean))
	gen sic2_resid = sic_resid*sic_resid if sic == 1
		sum sic2_resid
		local rmse_c = sqrt(r(mean))

	gen ya2_resid = ya_resid*ya_resid if sia == 1
		sum ya2_resid
		local rmse_ya = sqrt(r(mean))
	gen yb2_resid = yb_resid*yb_resid if sib == 1
		sum yb2_resid
		local rmse_yb = sqrt(r(mean))
	gen yc2_resid = yc_resid*yc_resid if sic == 1
		sum yc2_resid
		local rmse_yc = sqrt(r(mean))

	* Decidimos las salidas
	restore
	
	*Modelo A
	replace modelo_a_b = `betaa_est' in `i'
	replace modelo_a_g1 = `gamma1a_est' in `i'
	replace modelo_a_g2 = `gamma2a_est' in `i'
	replace modelo_a_rmse_sia = `rmse_a' in `i'
	replace modelo_a_rmse_ya = `rmse_ya'

	*Modelo B
	replace modelo_b_b = `betab_est' in `i'
	replace modelo_b_g1 = `gamma1b_est' in `i'
	replace modelo_b_g2 = `gamma2b_est' in `i'
	replace modelo_a_rmse_sib = `rmse_b' in `i'
	replace modelo_a_rmse_yb = `rmse_yb' in `i'

	*Modelo C
	replace modelo_c_b = `betac_est' in `i'
	replace modelo_c_g1 = `gamma1c_est' in `i'
	replace modelo_c_g2 = `gamma2c_est' in `i'
	replace modelo_a_rmse_sic = `rmse_c' in `i'
	replace modelo_a_rmse_yc = `rmse_yc' in `i'

	preserve
	clear
}
