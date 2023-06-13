*-----------------------------------------------------
* Econometría de datos de panel - UTDT 
* Trabajo práctico 2 2023
* Maestría en economía
* Oscar Jaramillo  
*-----------------------------------------------------

clear all
set more off 
cap log close
set matsize 500
set seed 1263

*-----------------------------------------------------------------------------
* Ejercicio 1. Propiedades de muestra finita en paneles no balanceados
*-----------------------------------------------------------------------------
*-----------------------------------------------------------------------------
* Instrucciones
*---------------
* El dofile está hecho para ser corrido de forma completa, se realizarán
* todas las simulaciones y al final hay una tabla resumen, no olvidar definir
* la cantidad de simulaciones y el tamaño de las muestras
* Lastimosamente el estimador Kiviet requiere mucha cpacidad de computo y puede
* ser necesario estimarlo de forma separada, el código para simular Kiviet 
* está separado 
*-----------------------------------------------------------------------------

* Instalamos los paquetes en caso de no tenerlos
*ssc install xtabond2, all replace
*ssc install xtlsdvc, all replace 

* Empezamos por definir lo que vamos a necesitar
local rep = 1000
global T = 2				// cantidad de periodos 
global N = 20				// cantidad de individuos 
global NT = $N * $T			
 
scalar beta = 1 
scalar gamma1 = 1
scalar gamma2 = 1

* Creamos variables para guardar las estimaciones
set obs `rep' 

*Modelo A
gen modelo_a_b = .
gen modelo_a_g1 = .
gen modelo_a_g2 = .

*Modelo B
gen modelo_b_b = .
gen modelo_b_g1 = .
gen modelo_b_g2 = .

*Modelo C
gen modelo_c_b = .
gen modelo_c_g1 = .
gen modelo_c_g2 = .


preserve

clear




*Simulaciones
forvalues i = 1(1)`rep' {
		
	preserve
		
	clear
		
	*Generamos panel 
	set obs $NT

	egen id   = seq(), f(1) t($N) b($T)
	egen tt = seq(), f(1) t($T)
	xtset id tt
		  
	*Generamos variables
	gen ci = .
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
		
	* Decidimos las salidas
	restore
		
	*LSDV
	replace betalsdv_yjt = `betalsdv_yjt' in `i'
	replace betalsdv_xjt = `betalsdv_xjt' in `i'
	replace lsdvtest = `lsdvtest' in `i'

	* AB-GMM1
	replace betaab1_yjt = `betaab1_yjt' in `i'
	replace betaab1_xjt = `betaab1_xjt' in `i'
	replace ab1test = `ab1test' in `i'
		
	* AB-GMM2
	replace betaab2_yjt = `betaab2_yjt' in `i'
	replace betaab2_xjt = `betaab2_xjt' in `i'
	replace ab2test = `ab2test' in `i'


	* BB-GMM1
	replace betabb_yjt = `betabb_yjt' in `i'
	replace betabb_xjt = `betabb_xjt' in `i'
	replace bbtest = `bbtest' in `i'

	* AH
	replace betaah_yjt = `betaah_yjt' in `i'
	replace betaah_xjt = `betaah_xjt' in `i'
	replace ahtest = `ahtest' in `i'
}
 ** Simulaciones estimador Kiviet
 
 forvalues i = 1(1)`rep' {

	preserve
		
	clear
		
	*Generamos panel 
	set obs $NT

	egen id   = seq(), f(1) t($N) b($T)
	egen tt = seq(), f(1) t($T)
	xtset id tt
		  
	*Generamos variables
	gen ci = .
	forvalues r = 1/$N {
		gen c = rnormal(0, 1)
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

	*** Kiviet
	** LSDVC (Kiviet, 1995)
	xtlsdvc yit , initial(ah) bias(2) vcov(50)
	
	* guardamos los betas de Kiviet
	matrix define kiv = e(b)
	scalar ki = kiv[1,1]
	local betaki_yjt = ki

	test L1.yit = 1
		
	* Guardamos el test
	local kitest = r(p)
	
	*decidimos las salidas
	restore
	
	* Kiviet
	replace betaki_yjt = `betaki_yjt' in `i'
	replace kitest = `kitest' in `i'
}
*Tabla media, mediana y desviaciones para yjt-1
tabstat betalsdv_yjt betaab1_yjt betaab2_yjt betabb_yjt betaah_yjt betaki_yjt , statistics( mean median sd ) columns(statistics)

*Tabla media, mediana y desviaciones para xjt
tabstat betalsdv_xjt betaab1_xjt betaab2_xjt betabb_xjt betaah_xjt betaki_xjt , statistics( mean median sd ) columns(statistics)

*Tabla tamaño del test alfa=0.5
tabstat lsdvtest ab1test ab2test bbtest ahtest kitest , statistics( mean median sd ) columns(statistics)
