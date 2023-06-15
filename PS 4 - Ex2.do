/* Datos de Panel  - UTDT - PS 4 - Ejercicio 2

Profesor: Martín González-Rozada
Ayudante: Fiona Franco Churruarín

*/

clear all
set more off

*** Abrimos la base de datos 'keane'
cd "C:\Users\fionafch\Dropbox\Teaching\UTDT MAECO Datos de Panel\2022\Clases Panel\Problem Set 4"
u ./data/keane
des

*** Declaramos el panel
xtset id year

*** La variable wage no se observa para todos los individuos en todos los periodos
count if lwage==.

*** FE
xtreg lwage exper educ, fe
outreg, store(t1) ctitle("", "FE") keep(exper educ) noautosumm se


***************************************************************************
*** Correccion del Sesgo de Seleccion Muestral: Truncamiento Incidental ***
***************************************************************************

* Inciso (a)

*** Wooldridge (1995). Enfoque Chamberlain (1980) ***
sum id
global g = r(max)

forvalues t = 81/83{
	gen exper`t'=.
	gen educ`t'=.
	
	forvalues i=1/$g{
	qui sum exper if id==`i'&year==`t'
	qui replace exper`t'=r(mean) if id==`i'
	qui sum educ if id==`i'&year==`t'
	qui replace educ`t'=r(mean) if id==`i'	
	}
}


**** 1. Pooled Probit (Chamberlain) y obtener inv. del coef de Mills
probit obswage exper educ exper81-educ83						// pooled probit
predict sit, xb
gen lambdait=normalden(sit)/normal(sit)							// obtenemos la inversa del cociente de Mills estimado

**** 2. Pooled OLS de la ec. de interes + X(chamberlain) + inv. Coef Mills  
reg lwage exper educ exper81-educ83 i.year#c.lambdait			// pooled ols
outreg, merge(t1) ctitle("", "W(95)-C(80)") keep(exper educ) noautosumm se

drop sit lambdait



* Inciso (b)

*** Wooldridge (1995). Enfoque Mundlak (1978) ***

egen mean_experi = mean(exper), by(id)
egen mean_educi = mean(educ), by(id)

**** 1. Pooled Probit (Mundlak) y obtener inv. del coef de Mills
probit obswage exper educ mean_experi mean_educi 						// pooled probit
predict sit, xb
gen lambdait=normalden(sit)/normal(sit)									// obtenemos la inversa del cociente de Mills estimado

**** 2. Pooled OLS de la ec. de interes + X(Munlak) + inv. Coef Mills  
reg lwage exper educ mean_experi mean_educi i.year#c.lambdait			// pooled ols
outreg, merge(t1) ctitle("", "W(95)-M(80)") keep(exper educ) noautosumm se

drop sit lambdait



* Inciso (d)


*** Computamos los SE via Bootstrapping ***

cap drop ebsel xb
** Paso 1. Computo coeficientes con la muestra original 

probit obswage exper educ mean_experi mean_educi

* Guardo los coeficientes de la estimacion de la ec. de seleccion
scalar bsel1=e(b)[1,1]
scalar bsel2=e(b)[1,2]
scalar bsel3=e(b)[1,3]
scalar bsel4=e(b)[1,4]
scalar bsel5=e(b)[1,5]

* Genero la serie de residuos de la ecuacion de seleccion
predict xb, xb			
gen ebsel = obswage - xb 

* Guardo los coeficientes de la estimacion de la ec. de interes 
gen lambdait=normalden(xb)/normal(xb)			
reg lwage exper educ mean_experi mean_educi i.year#c.lambdait
mat BETA = e(b)[1,1..2]

scalar beqw1=e(b)[1,1]
scalar beqw2=e(b)[1,2]
scalar beqw3=e(b)[1,3]
scalar beqw4=e(b)[1,4]
scalar beqw5=e(b)[1,12]

* Genero la serie de residuos de la ecuacion de interes
predict resid_eqw, resid

drop xb lambdait


** Paso 2. Programa para hacer una replicacion de bootstrap 
program drop _all

program onebootrep, rclass

* Utilizo un procedimiento de muestreo aleatorio simple con reemplazo para obtener una nueva muestra de errores de la ec. de seleccion y de la ec. de interes
preserve
bsample
mkmat ebsel, matrix(ebsel)
mkmat resid_eqw, matrix(resid_eqw)
restore
svmat ebsel
svmat resid_eqw

* Construyo una nueva variable dependiente de la ec. de seleccion
gen bs_obswage = bsel5+bsel1*exper+bsel2*educ+bsel3*mean_experi+bsel4*mean_educi+ebsel1

replace bs_obswage = 1 if bs_obswage>1
replace bs_obswage = 0 if bs_obswage<0
replace bs_obswage = round(bs_obswage)

* Estimacion de la ecuacion de seleccion
probit bs_obswage exper mean_experi educ mean_educi
return scalar bsel_exper = e(b)[1,1]
return scalar bsel_educ = e(b)[1,3]

* Construyo la inversa del cociente de Mills estimada
predict sit, xb
gen lambdait=normalden(sit)/normal(sit)									

* Construyo una nueva variable dependiente de la ec. de salarios
gen bs_lwage = .
replace bs_lwage = beqw5+beqw1*exper+beqw2*educ+beqw3*mean_experi+beqw4*mean_educi+resid_eqw1 if bs_obswage==1

*Estimacion de la ecuacion de salarios
reg bs_lwage exper educ mean_experi mean_educi i.year#c.lambdait
return scalar beqw_exper = e(b)[1,1]
return scalar beqw_educ = e(b)[1,2]

drop ebsel1 resid_eqw1 bs_obswage sit lambdait bs_lwage
end

** Ahora hacemos 50 replicaciones de bootstrap
preserve

global b = 50
simulate bsel_exper=r(bsel_exper) bsel_educ=r(bsel_educ) beqw_exper=r(beqw_exper) beqw_educ=r(beqw_educ), seed(2021) reps($b) saving(bootdata, replace): onebootrep

** Computo de los errores estandar de bootstrap 
for any beqw_exper beqw_educ: egen meanX=mean(X)

gen exper= beqw_exper- meanbeqw_exper
gen educ = beqw_educ- meanbeqw_educ

mat accum Vaux= exper educ , nocons

mat V =  Vaux/($b-1)

mat list V

mat b = BETA

ereturn post b V

ereturn display

outreg, merge(t1) ctitle("", "W(95)-M(80)-BS") keep(exper educ) noautosumm se

restore



