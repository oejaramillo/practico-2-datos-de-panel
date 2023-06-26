*-----------------------------------------------------
* Econometría de datos de panel - UTDT 
* Trabajo práctico 2 2023
* Maestría en economía
* Oscar Jaramillo  
*-----------------------------------------------------

clear all
set more off
cap log close
set matsize 1000
set seed 1263

*-----------------------------------------------------------------------------
* Ejercicio 1. Propiedades de muestra finita en paneles no balanceados
* Caso 2. Bootstrap
*-----------------------------------------------------------------------------
*-----------------------------------------------------------------------------
* Instrucciones
*---------------
* El dofile está hecho para ser corrido en partes, cada literal de forma 
* separada, calcula los literales a y b y en la parte c realiza el remuestreo 
* de bootstrap a través de una simulación.
* Al final hay una tabla resumen de los intervalos de confianza, no olvidar 
* definir la cantidad de simulaciones en la parte c y el tamaño de las muestras
* en la parte inicial.
*-----------------------------------------------------------------------------

*Generamos panel 
global T = 10				// cantidad de periodos 
global N = 100				// cantidad de individuos 
global NT = $N * $T	
 
*Parámetros
scalar beta = 1 
scalar gamma1 = 1
scalar gamma2 = 1

*------------------------------------
* a)
*------------------------------------
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

*------------------------------------
* b)
*------------------------------------
* Estimación de Wooldridge de la ecuación de selección
* Probit Modelo A
probit sia xjt zjt
* Matrices para guardar los intervalos de confianza
matrix probita = r(table)
matrix ci_gamma1_a = J(2,2,.)
matrix ci_gamma1_a[1,1] = (probita[5, 1], probita[6, 1])
matrix ci_gamma2_a = J(2,2,.)
matrix ci_gamma2_a[1,1] = (probita[5, 2], probita[6, 2])
*--------------------------------------------------------	
	scalar gamma1a_est = _b[xjt]
	scalar gamma2a_est = _b[zjt]
	predict sia_est, xb
	gen sia_resid = sia - sia_est
gen lambdait_a = normalden(sia_est)/normal(sia_est) if sia == 1
replace lambdait_a = 0 if missing(lambdait_a)

* Probit Modelo B
probit sib xjt zjt zj
* Matrices para guardar los intervalos de confianza
matrix probitb = r(table)
matrix ci_gamma1_b = J(2,2,.)
matrix ci_gamma1_b[1,1] = (probitb[5, 1], probitb[6, 1])
matrix ci_gamma2_b = J(2,2,.)
matrix ci_gamma2_b[1,1] = (probitb[5, 2], probitb[6, 2])
*--------------------------------------------------------	
	scalar gamma1b_est = _b[xjt]
	scalar gamma2b_est = _b[zjt]
	predict sib_est, xb
	gen sib_resid = sib - sib_est
gen lambdait_b = normalden(sib_est)/normal(sib_est) if sib == 1
replace lambdait_b = 0 if missing(lambdait_b)

* Probit Modelo C
probit sic xjt zjt zj
* Matrices para guardar los intervalos de confianza
matrix probitc = r(table)
matrix ci_gamma1_c = J(2,2,.)
matrix ci_gamma1_c[1,1] = (probitc[5, 1], probitc[6, 1])
matrix ci_gamma2_c = J(2,2,.)
matrix ci_gamma2_c[1,1] = (probitc[5, 2], probitc[6, 2])
*--------------------------------------------------------	
	scalar gamma1c_est = _b[xjt]
	scalar gamma2c_est = _b[zjt]
	predict sic_est, xb
	gen sic_resid = sic - sic_est
gen lambdait_c = normalden(sic_est)/normal(sic_est) if sic == 1
replace lambdait_c = 0 if missing(lambdait_c)

* Estimación de Woolridge de la ecuación de interés
* Modelo A
reg ya xjt i.tt#c.lambdait_a if sia == 1
* Matrices para guardar los intervalos de confianza
matrix woolda = r(table)
matrix ci_beta_a = J(2,2,.)
matrix ci_beta_a[1,1] = (woolda[5, 1], woolda[6, 1])
*--------------------------------------------------------	
	scalar betaa_est = _b[xjt]
	predict ya_resid, resid

* Modelo B
reg yb xjt xj i.tt#c.lambdait_b if sib == 1
* Matrices para guardar los intervalos de confianza
matrix wooldb = r(table)
matrix ci_beta_b = J(2,2,.)
matrix ci_beta_b[1,1] = (wooldb[5, 1], wooldb[6, 1])
*--------------------------------------------------------	
	scalar betab_est = _b[xjt]
	predict yb_resid, resid

* Modelo C
reg yc xjt xj i.tt#c.lambdait_c if sic == 1
* Matrices para guardar los intervalos de confianza
matrix wooldc = r(table)
matrix ci_beta_c = J(2,2,.)
matrix ci_beta_c[1,1] = (wooldc[5, 1], wooldc[6, 1])
*--------------------------------------------------------	
	scalar betac_est = _b[xjt]
	predict yc_resid, resid
	

*------------------------------------
* C)
*------------------------------------
* Opciones de la simulación
global sim = 1000					// cantidad de simulaciones
global modelo = 3				// Elegimos el modelo a simular (A=1, B=2 o C=3)

* Programa de bootstrap
program drop _all

* Programa diferente para cada modelo, es necesario correr todos
*Modelo A
{
program onebootrep_a, rclass
	*(i)
	set more off
	
	preserve
	
	bsample if id
	
	mkmat sia_resid, matrix(sia_resid)
	mkmat ya_resid, matrix(ya_resid)
	
	matrix b_sia_resid = sia_resid
	matrix b_ya_resid = ya_resid
	
	restore
	
	svmat b_sia_resid
	svmat b_ya_resid

	list b_sia_resid b_ya_resid, sepby(id)
	
	* Ecuación de selección
	gen sjt2 = gamma1a_est * xjt + gamma2a_est * zjt + alfaja + b_sia_resid
	gen si2 = (sjt2 > 0)	

	*(ii)
	probit si2 xjt zjt   // modelo A
	
	return scalar gamma1_est2 = _b[xjt]
	return scalar gamma2_est2 = _b[zjt]	

	predict si2_est, xb
	gen lambdait_2 = normalden(si2_est)/normal(si2_est) if si2 == 1
	
	* Ecuación de interés
	gen yjt2 = betaa_est * xjt + cja + b_ya_resid
	
	*(iii)
	gen y2 = yjt2 * si2
	
	reg y2 xjt i.t#c.lambdait_2 if si2 == 1
		
	return scalar beta_est2 = _b[xjt]

	drop b_sia_resid b_ya_resid sjt2 si2 si2_est lambdait_2 yjt2 y2
end
}

*--------------------------------------------------------------
*Modelo B
{
program onebootrep_b, rclass
	*(i)
	set more off
	
	preserve
	
	bsample if id
	
	mkmat sib_resid, matrix(sib_resid)
	mkmat yb_resid, matrix(yb_resid)
	
	matrix b_sib_resid = sib_resid
	matrix b_yb_resid = yb_resid
	
	restore
	
	svmat b_sib_resid
	svmat b_yb_resid

	list b_sib_resid b_yb_resid, sepby(id)
	
	* Ecuación de selección
	gen sjt2 = gamma1b_est * xjt + gamma2b_est * zjt + alfajb + b_sib_resid
	gen si2 = (sjt2 > 0)	

	*(2)
	probit si2 xjt zjt zj   // modelo B
	
	return scalar gamma1_est2 = _b[xjt]
	return scalar gamma2_est2 = _b[zjt]	

	predict si2_est, xb
	gen lambdait_2 = normalden(si2_est)/normal(si2_est) if si2 == 1
	
	* Ecuación de interés
	gen yjt2 = betab_est * xjt + cjb + b_yb_resid
	
	*(3)
	gen y2 = yjt2 * si2
	
	reg y2 xjt xj i.t#c.lambdait_2 if si2 == 1
		
	return scalar beta_est2 = _b[xjt]

	drop b_sib_resid b_yb_resid sjt2 si2 si2_est lambdait_2 yjt2 y2
end
}

*--------------------------------------------------------------
*Modelo C
{
program onebootrep_c, rclass
	*(i)
	set more off
	
	preserve
	
	bsample if id
	
	mkmat sic_resid, matrix(sic_resid)
	mkmat yc_resid, matrix(yc_resid)
	
	matrix b_sic_resid = sic_resid
	matrix b_yc_resid = yc_resid
	
	restore
	
	svmat b_sic_resid
	svmat b_yc_resid

	list b_sic_resid b_yc_resid, sepby(id)
	
	* Ecuación de selección
	gen sjt2 = gamma1c_est * xjt + gamma2c_est * zjt + alfajb + b_sic_resid
	gen si2 = (sjt2 > 0)	

	*(2)
	probit si2 xjt zjt zj   // modelo C
	
	return scalar gamma1_est2 = _b[xjt]
	return scalar gamma2_est2 = _b[zjt]	

	predict si2_est, xb
	gen lambdait_2 = normalden(si2_est)/normal(si2_est) if si2 == 1
	
	* Ecuación de interés
	gen yjt2 = betac_est * xjt + cjc + b_yc_resid
	
	*(3)
	gen y2 = yjt2 * si2
	
	reg y2 xjt xj i.t#c.lambdait_2 if si2 == 1
		
	return scalar beta_est2 = _b[xjt]

	drop b_sic_resid b_yc_resid sjt2 si2 si2_est lambdait_2 yjt2 y2
end
}

*---------------------------------------------------------------
* Corremos la simulación según el modelo elegido
if $modelo == 1 {
	simulate a_gamma1 = r(gamma1_est2) a_gamma2 = r(gamma2_est2) a_beta = r(beta_est2), reps($sim) nodots: onebootrep_a
}
if $modelo == 2 {
	simulate b_gamma1 = r(gamma1_est2) b_gamma2 = r(gamma2_est2) b_beta = r(beta_est2), reps($sim) nodots: onebootrep_b
}
if $modelo == 3 {
	simulate c_gamma1 = r(gamma1_est2) c_gamma2 = r(gamma2_est2) c_beta = r(beta_est2), reps($sim) nodots: onebootrep_c
}
else {
	display "*************************** Elegir un modelo correcto A=1, B=2 o C=3 ***************************"
}
*------------------------------------------------------------------------
*------------------------------------------------------------------------
* Buscamos el intervalo de confianza  para bootstrap
if $modelo == 1 {
	sort a_gamma1
	centile a_gamma1, centile(2.5 97.5)
	matrix ci_gamma1_a[2,1] = r(ub_1)
	matrix ci_gamma1_a[2,2] = r(ub_2)
	
	sort a_gamma2
	centile a_gamma2, centile(2.5 97.5)
	matrix ci_gamma2_a[2,1] = r(ub_1)
	matrix ci_gamma2_a[2,2] = r(ub_2)
	
	sort a_beta
	centile a_beta, centile(2.5 97.5)
	matrix ci_beta_a[2,1] = r(ub_1)
	matrix ci_beta_a[2,2] = r(ub_2)
}
if $modelo == 2 {
	sort b_gamma1
	centile b_gamma1, centile(2.5 97.5)
	matrix ci_gamma1_b[2,1] = r(ub_1)
	matrix ci_gamma1_b[2,2] = r(ub_2)
	
	sort b_gamma2
	centile b_gamma2, centile(2.5 97.5)
	matrix ci_gamma2_b[2,1] = r(ub_1)
	matrix ci_gamma2_b[2,2] = r(ub_2)
	
	sort b_beta
	centile b_beta, centile(2.5 97.5)
	matrix ci_beta_b[2,1] = r(ub_1)
	matrix ci_beta_b[2,2] = r(ub_2)
}
if $modelo == 3 {
	sort c_gamma1
	centile c_gamma1, centile(2.5 97.5)
	matrix ci_gamma1_c[2,1] = r(ub_1)
	matrix ci_gamma1_c[2,2] = r(ub_2)
	
	sort c_gamma2
	centile c_gamma2, centile(2.5 97.5)
	matrix ci_gamma2_c[2,1] = r(ub_1)
	matrix ci_gamma2_c[2,2] = r(ub_2)
	
	sort c_beta
	centile c_beta, centile(2.5 97.5)
	matrix ci_beta_c[2,1] = r(ub_1)
	matrix ci_beta_c[2,2] = r(ub_2)
}
*Comparamos los intervalos de confianza para cada estimador y el modelo elegido
if $modelo == 1 {
	matrix rownames ci_gamma1_a = "Real" "Boootstrap"
	matrix colnames ci_gamma1_a = "Inferior" "Superior"
	matlist ci_gamma1_a
	
	matrix rownames ci_gamma2_a = "Real" "Boootstrap"
	matrix colnames ci_gamma2_a = "Inferior" "Superior"
	matlist ci_gamma2_a
	
	matrix rownames ci_beta_a = "Real" "Boootstrap"
	matrix colnames ci_beta_a = "Inferior" "Superior"
	matlist ci_beta_a
}
if $modelo == 2 {
	matrix rownames ci_gamma1_b = "Real" "Boootstrap"
	matrix colnames ci_gamma1_b = "Inferior" "Superior"
	matlist ci_gamma1_b
	
	matrix rownames ci_gamma2_b = "Real" "Boootstrap"
	matrix colnames ci_gamma2_b = "Inferior" "Superior"
	matlist ci_gamma2_b
	
	matrix rownames ci_beta_b = "Real" "Boootstrap"
	matrix colnames ci_beta_b = "Inferior" "Superior"
	matlist ci_beta_b
}
if $modelo == 3 {
	matrix rownames ci_gamma1_c = "Real" "Boootstrap"
	matrix colnames ci_gamma1_c = "Inferior" "Superior"
	matlist ci_gamma1_c
	
	matrix rownames ci_gamma2_c = "Real" "Boootstrap"
	matrix colnames ci_gamma2_c = "Inferior" "Superior"
	matlist ci_gamma2_c
	
	matrix rownames ci_beta_c = "Real" "Boootstrap"
	matrix colnames ci_beta_c = "Inferior" "Superior"
	matlist ci_beta_c
}
