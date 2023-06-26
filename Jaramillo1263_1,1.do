*-----------------------------------------------------
* Econometría de datos de panel - UTDT 
* Trabajo práctico 2 2023
* Maestría en economía
* Oscar Jaramillo  
*-----------------------------------------------------

clear all
set more off
cap log close
set seed 1263

*-----------------------------------------------------------------------------
* Ejercicio 1. Propiedades de muestra finita en paneles no balanceados
*-----------------------------------------------------------------------------
*-----------------------------------------------------------------------------
* Instrucciones
*---------------
* El dofile está hecho para ser corrido de forma completa, se realizarán
* todas las simulaciones y al final hay una tabla resumen, no olvidar definir
* la cantidad de simulaciones y el tamaño de las muestras, en el caso que 
* ocurra algún problema, se recomienda correr el do hasta el final del loop y
* luego generar el resumen del experimento en la segunda parte del código.
*-----------------------------------------------------------------------------

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
gen modelo_a_rmse_si = .
gen modelo_a_rmse_y = .

*Modelo B
gen modelo_b_b = .
gen modelo_b_g1 = .
gen modelo_b_g2 = .
gen modelo_b_rmse_si = .
gen modelo_b_rmse_y = .

*Modelo C
gen modelo_c_b = .
gen modelo_c_g1 = .
gen modelo_c_g2 = .
gen modelo_c_rmse_si = .
gen modelo_c_rmse_y = .


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
	replace modelo_a_rmse_si = `rmse_a' in `i'
	replace modelo_a_rmse_y = `rmse_ya'

	*Modelo B
	replace modelo_b_b = `betab_est' in `i'
	replace modelo_b_g1 = `gamma1b_est' in `i'
	replace modelo_b_g2 = `gamma2b_est' in `i'
	replace modelo_b_rmse_si = `rmse_b' in `i'
	replace modelo_b_rmse_y = `rmse_yb' in `i'

	*Modelo C
	replace modelo_c_b = `betac_est' in `i'
	replace modelo_c_g1 = `gamma1c_est' in `i'
	replace modelo_c_g2 = `gamma2c_est' in `i'
	replace modelo_c_rmse_si = `rmse_c' in `i'
	replace modelo_c_rmse_y = `rmse_yc' in `i'

	preserve
	clear
}

*----------------------------------------------------------------------------
* Generamos un resumen del experimento
*----------------------------------------------------------------------------

*Desviación media absoluta
*Modelo A
egen mad_a_beta = mad(modelo_a_b)
egen mad_a_gamma1 = mad(modelo_a_g1)
egen mad_a_gamma2 = mad(modelo_a_g2)

*Modelo B
egen mad_b_beta = mad(modelo_b_b)
egen mad_b_gamma1 = mad(modelo_b_g1)
egen mad_b_gamma2 = mad(modelo_b_g2)

*Modelo C
egen mad_c_beta = mad(modelo_c_b)
egen mad_c_gamma1 = mad(modelo_c_g1)
egen mad_c_gamma2 = mad(modelo_c_g2)

*Sesgo de los estimadores
*Modelo A
gen sesgo_a_b = modelo_a_b - beta
gen sesgo_a_g1 = modelo_a_g1 - gamma1
gen sesgo_a_g2 = modelo_a_g2 - gamma2

*Modelo B
gen sesgo_b_b = modelo_b_b - beta
gen sesgo_b_g1 = modelo_b_g1 - gamma1
gen sesgo_b_g2 = modelo_b_g2 - gamma2

*Modelo C
gen sesgo_c_b = modelo_c_b - beta
gen sesgo_c_g1 = modelo_c_g1 - gamma1
gen sesgo_c_g2 = modelo_c_g2 - gamma2


*Resumen del sesgo de los estimadores
*beta
tabstat sesgo_a_b sesgo_b_b sesgo_c_b, statistics(mean, median, sd)
*gamma1
tabstat sesgo_a_g1 sesgo_b_g1 sesgo_c_g1, statistics(mean, median, sd)
*gamma2
tabstat sesgo_a_g2 sesgo_b_g2 sesgo_c_g2 , statistics(mean, median, sd)

*RMSE
*Modelo A
tabstat modelo_a_rmse_si modelo_a_rmse_y , statistics(mean)
*Modelo B
tabstat modelo_b_rmse_si modelo_b_rmse_y , statistics(mean)
*Modelo C
tabstat modelo_c_rmse_si modelo_c_rmse_y , statistics(mean)

*DMA
*beta
tabstat mad_a_beta mad_b_beta mad_c_beta
*gamma1
tabstat mad_a_gamma1 mad_b_gamma1 mad_c_gamma1
*gamma2
tabstat mad_a_gamma2 mad_b_gamma2 mad_c_gamma2
