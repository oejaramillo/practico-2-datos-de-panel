clear all
*Generamos panel 
global T = 2				// cantidad de periodos 
global N = 20				// cantidad de individuos 
global NT = $N * $T	
set obs $NT

egen id   = seq(), f(1) t($N) b($T)
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
gen psi2 = rnormal(0,1)
gen psi3 = rnormal(0,1)
gen psi4 = rnormal(0,1)

gen epsilonjt = rnormal(0,1)
gen ujt = 0.6*epsilonjt + 0.8*psi1

*Modelo A
gen alfaj = psi2 + psi4
gen cj = psi3 + psi4

*Modelo B
gen alfaj = psi2 + (L.zjt + L2.zjt)/2
gen cj = psi3 + (L.zjt + L2.zjt)/2

*Modelo C
gen alfaj = psi2 + (L.zjt + L2.zjt)/2 + psi4
gen cj = psi3 + (L.zjt + L2.zjt)/2 + psi4




 


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
