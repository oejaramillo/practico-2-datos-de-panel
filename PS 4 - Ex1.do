/* Datos de Panel  - UTDT - PS 4 - Ejercicio 1

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
outreg, store(t1) keep(exper educ) noautosumm se



*********************************************************
*** Test de Wooldridge (1995). Enfoque Mundlak (1978) ***
*********************************************************

egen mean_experi = mean(exper), by(id)
egen mean_educi = mean(educ), by(id)

*1) Pooled probit
probit obswage exper mean_experi educ mean_educi 	// 'obswage' representa la variable 'sit' de las slides

*2) Generamos la estimacion de la inversa del cociente de Mills
predict xb, xb
gen lambdait=normalden(xb)/normal(xb)			

*3) Contraste t para la estimacion de la inversa del cociente de Mills
xtreg lwage exper educ lambdait, fe vce(cluster id)	// hay evidencia en favor de rechazar la h0 

drop xb lambdait



*************************************************************
*** Test de Wooldridge (1995). Enfoque Chamberlain (1980) ***
*************************************************************

sum id
global g=r(max)			// Numero de individuos en la muestra

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

*1) Pooled probit
probit obswage exper educ exper81 educ81 exper82 educ82 exper83 educ83

*2) Generamos la estimacion de la inversa del cociente de Mills
predict xb, xb
gen lambdait=normalden(xb)/normal(xb)			

*3) Contraste t para la estimacion de la inversa del cociente de Mills
xtreg lwage exper educ lambdait, fe vce(cluster id)	// hay evidencia en favor de rechazar la h0 

