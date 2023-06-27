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
cd "C:\edu\maestria\Torcuatto di Tella\cuarto trimestre\datos de panel\segunda cursada\practico-2-datos-de-panel"
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

matrix b_probit = J(2,1,.)
matrix rownames b_probit = "employj,t−1=1" "employj,t−1=0"
matrix colnames b_probit = "Probabilidad"

matrix b_probit[1,1] = normal(_b[_cons]+_b[l.employ])
matrix b_probit[2,1] = normal(_b[_cons])

matlist b_probit

* c) Adicione al modelo estimado en 1. el conjunto completo de variables binarias temporales...

probit employ l.employ y83-y87 if black==1, vce(cluster id)

matrix c_probit = J(5, 2, .)
matrix rownames c_probit = "Año 83" "Año 84" "Año 85" "Año 86" "Año 87"
matrix colnames c_probit = "employj,t−1=1" "employj,t−1=0"


forvalues i = 83/87{
	scalar p1`i' = normal(_b[_cons]+_b[l.employ]+_b[y`i'])
	scalar p2`i' = normal(_b[_cons]+_b[y`i'])
}
matrix c_probit[1,1] = p183
matrix c_probit[2,1] = p184
matrix c_probit[3,1] = p185
matrix c_probit[4,1] = p186
matrix c_probit[5,1] = p187

matrix c_probit[1,2] = p283
matrix c_probit[2,2] = p284
matrix c_probit[3,2] = p285
matrix c_probit[4,2] = p286
matrix c_probit[5,2] = p287

matrix list c_probit

* d) Ahora estime un modelo de efectos no observables dinámico...

gen employ81 = employ if year==81

replace employ81 = employ[_n-1] if year==82
replace employ81 = employ[_n-2] if year==83
replace employ81 = employ[_n-3] if year==84
replace employ81 = employ[_n-4] if year==85
replace employ81 = employ[_n-5] if year==86
replace employ81 = employ[_n-6] if year==87

xtprobit employ l.employ employ81 y83-y87 if black==1, re

* f) Promedie las probabilidades estimadas a lo largo de employj,81...

gen prob1 = normal((_b[_cons]+_b[l.employ]+_b[employ81]*employ81+_b[y87])/sqrt(1+e(sigma_u)^2)) if black==1 & year==87
sum prob1

gen prob2 = normal((_b[_cons]+_b[employ81]*employ81+_b[y87])/sqrt(1+e(sigma_u)^2)) if black==1 & year==87
sum prob2
