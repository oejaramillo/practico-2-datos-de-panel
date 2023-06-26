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
* Ejercicio 2. Modelos de Panel no Lineales
*-----------------------------------------------------------------------------
*-----------------------------------------------------------------------------
* Datos
*---------------
cd "C:\Oscar Local\di tella\practico-2-datos-de-panel"
use ".\keane.dta"

* Declaramos el panel
*----------------------
order id year, first
xtset id year

xtdescribe

*-----------------------------------------------------------------------------
* a) Use Pooled Probit para estimar el modelo...

probit employ l.employ if black==1, vce(cluster id)

* b) Estime P(employjt = 1| employj,t−1 = 1) y P(employjt = 1| employj,t−1 = 0). Explique cómo obtendría errores estándar para estas estimaciones.

disp in red "La probabilidad elegir trabajar como white collar, blue collar o actividades de servicios dado que estaba en dichos sectores el periodo anterior es " normal(_b[_cons]+_b[l.employ])
disp in red "La probabilidad elegir trabajar como white collar, blue collar o actividades de servicios dado que no estaba en dichos sectores el periodo anterior es " normal(_b[_cons])



/* c) Adicione al modelo estimado en 1. el conjunto completo de variables binarias temporales y estime las probabilidades pedidas en (b). para 1987. Hay diferencias importantes con las estimaciones de (b).?*/

* Como ya tenemos las variables creadas (omitimos el año 81 por el lag y el 82 como bechmarck)

probit employ l.employ y83-y87 if black==1, vce(cluster id)

forvalues t = 83/87{

disp in blue "Año " `t'
disp in red "La probabilidad de elegir trabajar como white collar, blue collar o actividades de servicios en el año "19`t' " dado que estaba en dichos sectores el periodo anterior es " normal(_b[_cons]+_b[l.employ]+_b[y`t'])
disp in red "La probabilidad de elegir trabajar como white collar, blue collar o actividades de servicios en el año "19`t' " dado que estaba en dichos sectores el periodo anterior es " normal(_b[_cons]+_b[y`t'])

}


/* (d) Ahora estime un modelo de efectos no observables dinámico. En particular, adicione employj,81 como una variable explicativa adicional y use el modelo Probit de efectos aleatorios. Use el conjunto completo de variables binarias temporales.*/

gen employ81 = employ if year==81
replace employ81 = employ[_n-1] if year==82
replace employ81 = employ[_n-2] if year==83
replace employ81 = employ[_n-3] if year==84
replace employ81 = employ[_n-4] if year==85
replace employ81 = employ[_n-5] if year==86
replace employ81 = employ[_n-6] if year==87

xtprobit employ l.employ employ81 y83-y87 if black==1, re

/* (f) Promedie las probabilidades estimadas a lo largo de employj,81 para obtener el efecto parcial promedio para 1987. Compare estas estimaciones con aquellas en (c).
*/


gen prob = normal((_b[_cons]+_b[l.employ]+_b[employ81]*employ81+_b[y87])/sqrt(1+e(sigma_u)^2)) if black==1 & year==87
sum prob


gen prob = normal((_b[_cons]+_b[employ81]*employ81+_b[y87])/sqrt(1+e(sigma_u)^2)) if black==1 & year==87
sum prob
