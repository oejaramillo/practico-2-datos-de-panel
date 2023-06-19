clear all
set more off
set seed 1263
*Generamos panel 
global T = 2				// cantidad de periodos 
global N = 20				// cantidad de individuos 
global NT = $N * $T	
set obs $NT

egen id = seq(), f(1) t($N) b($T)
egen tt = seq(), f(1) t($T)
xtset id tt
		  
*Parámetros
scalar beta = 1 
scalar gamma1 = 1
scalar gamma2 = 1

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

* Estimación de Wooldridge bajo el enfoque de Chamberlain
* Probit Modelo A
probit sia xjt zjt
predict sia_est, xb
gen lambdait_a = normalden(sia_est)/normal(sia_est) if sia == 1
replace lambdait_a = 0 if missing(lambdait_a)

* Probit Modelo B
probit sib xjt zjt zj
predict sib_est, xb
gen lambdait_b = normalden(sib_est)/normal(sib_est) if sib == 1
replace lambdait_b = 0 if missing(lambdait_b)

* Probit Modelo C
probit sic xjt zjt zj
predict sic_est, xb
gen lambdait_c = normalden(sic_est)/normal(sic_est) if sic == 1
replace lambdait_c = 0 if missing(lambdait_c)

* Estimación de Woolridge bajo el enfoque de Mundlak
* Modelo A
reg ya xjt i.tt#c.lambdait_a if sia == 1
* Modelo B
reg yb xjt xj i.tt#c.lambdait_b if sib == 1
* Modelo C
reg yc xjt xj i.tt#c.lambdait_c if sic == 1





 


gen xjt = .

forvalues r = 1/$N {
	gen c = rnormal(0,1)
	replace ci = c if id == `r'
	drop c
}
		 
gen uit = rnormal(0, 1)
		
gen vit = rnormal(0,0.9)
		
gen xit = .
replace xit = 0 if tt == 1
replace xit = 0.8 * L1.xit + vit if tt != 1

gen yit = .
replace yit = 0 if tt == 1
replace yit = alfa * L1.yit + beta * xit + ci + uit if tt! = 1
	
drop if tt <= 10
		
	* Realizamos las estimaciones
	*** LSDV
	xtreg yit L1.yit xit, fe
		
	* guardamos los betas de LSDV
	local betalsdv_yjt = _b[L1.yit]
	local betalsdv_xjt = _b[xit]
		
	test L1.yit = 0.5
		
	* Guardamos el test
	local lsdvtest = r(p)  
		
	*** AB-GMM1 one-step GMM 
	xtabond2 yit L1.yit xit, gmm(L1.yit) iv(xit) nolevel robust
		
	* guardamos los betas de AB-GMM1 one-step GMM 
	local betaab1_yjt = _b[L1.yit]
	local betaab1_xjt = _b[xit]
		
	test L1.yit = 0.5
		
	* Guardamos el test
	local ab1test = r(p)

	*** AB-GMM2 two-step GMM  
	xtabond2 yit L1.yit xit, gmm(L1.yit) iv(xit) nolevel twostep
		
	* guardamos los betas de AB-GMM2 one-step GMM 
	local betaab2_yjt = _b[L1.yit]
	local betaab2_xjt = _b[xit]
		
	test L1.yit = 0.5
		
	* Guardamos el test
	local ab2test = r(p)

	*** BB-GMM1 one-step GMM 
	xtabond2 yit L1.yit xit, gmm(L1.yit) iv(xit) twostep
		
	* guardamos los betas de AB-GMM2 one-step GMM 
	local betabb_yjt = _b[L1.yit]
	local betabb_xjt = _b[xit]
		
	test L1.yit = 0.5
		
	* Guardamos el test
	local bbtest = r(p)

	*** AH
	ivregress 2sls D.yit (D.L1.yit = L2.yit) D.(xit), nocons
		
	* guardamos los betas de AH
	local betaah_yjt = _b[D.L1.yit]
	local betaah_xjt = _b[D.xit]
		
	test D.L1.yit = 0.5
		
	* Guardamos el test
	local ahtest = r(p)	
