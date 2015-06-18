#Script for the analytical calculations

import os, sys, math
##         COEFFICIENT OF HEAD LOSSES: LAMBDA
##
## 	"LAMBDA" coeffcient as a power low of the Reynolds number (Handbook of heat transfert)
##
##	==> The power law is of the form LAMBDA(Re) = a + b.Re^(c)
##
##	Examples (available notably for Re < 10000): 
##		* Blasius: a = 0 ; b = 0.0791 ; c = -0.25 
##		* Drew and al.: a = 0.00128 ; b = 0.1143 ; c = -0.311
##		* Drew and al. (bis): a = 0.0014 ; b = 0.125 ; c = -0.332
##		* Bhatti and Shah: a = 0.00128 ; b = 0.1143 ; c = -0.311
##		* Laminar: a = 0 ; b = 64 ; c = -1
##
##
##	The head losses rho*g*Delta_H is expressed as : Delta_P = rho * g * LAMBDA * PDC * U**2
##	with: PDC = L/(4g*D_H), and U is the velocity
##
######################
def properties(PG):
	# ouverture des fichiers
	fic = open(PG,'r')
	fichier = fic.readlines()
	l = fichier [0]
	tLi = l.split()
	mu = float(tLi[0])
	D = float(tLi[5])
	P = float(tLi[6])
	L = float(tLi[7])
	Cp = float(tLi[2])
	Pth = float(tLi[4])
	Tinlet = float(tLi[8])
	Pvol = float(tLi[9])
	Uinlet = float(tLi[10])
	gg = float(tLi[11])
	Q = float(tLi[12])
	fic.close()
	return mu,D,P,L,Cp,Pth,Tinlet,Pvol,Uinlet,gg,Q



################
def analytic1(mu,D,P,L,Cp,Pth,Tinlet,Pvol,Uinlet,gg,Q):
#parameters:
	R = 287.06
	S = D*P
## Hydraulic diamter and quantities associated:
	D_H = 2*D*P/(D+P)
	PDC = L/(2*gg*D_H)
	Re = 2*Q/(mu*(D+P))
##
# yy --> array for y coordinate
# TT --> array for temperature
# UU --> array for velocity
# DP --> array for pressure without head losses
# DP_lam --> 	array for pressure with head losses
# 400 points have been choosen for the curve
	yy = [0]*400
	TT = [0]*400
	Rho = [0]*400
	UU = [0]*400
	DP = [0]*400
	DP_lam = [0]*400
# coefficient of the power law
	a = 0.00128
	b = 0.1143
	c = -0.311
	LAMBDA = (a + b*Re**(c))
#
	f1  = open('analytic.dat', 'w')	
	f2 = open('Tbulk.dat', 'w')	
	for p in range(400):
		yy[p] = p*0.8768/400
		TT[p] = Tinlet + (Pvol*D*P)/(Q*Cp) * yy[p]
		Rho[p] = Pth / (R*TT[p])
		UU[p] = Q/(S*Rho[p])
		DP[p] = (Rho[p]*UU[p]*UU[p] - Rho[0]*UU[0]*UU[0]) 
		DP_lam[p] = DP[p] + 0.5*Rho[p]*gg*(LAMBDA*PDC*UU[p]*UU[p])
#
		f1.write('%18.3f %18.3f %18.6f %18.2f\n' % (yy[p], TT[p], Rho[p], UU[p]))
	f2.write('%18.3f\n' % (TT[p]))
#
	return DP_lam, yy, DP
##
def analytic2(DP_lam, yy, DP):
	f2  = open('P_drop.dat', 'w')
	for p in range(400):
		q =400 - p-1
		XX = DP_lam[q]
		YY = yy[p]
		ZZ = DP[q]
		f2.write('%18.2f %18.3f %18.3f\n' % (yy[p], ZZ, XX))
##
##
if __name__ == '__main__':
	args = sys.argv
	PG = args[1]
	mu,D,P,L,Cp,Pth,Tinlet,Pvol,Uinlet,gg,Q = properties(PG)
	DP_lam, yy, DP = analytic1(mu,D,P,L,Cp,Pth,Tinlet,Pvol,Uinlet,gg,Q)
	analytic2(DP_lam, yy, DP)

