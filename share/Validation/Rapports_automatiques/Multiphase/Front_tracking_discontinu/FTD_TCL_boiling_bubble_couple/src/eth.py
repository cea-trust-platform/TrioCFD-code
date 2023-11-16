import math
mu = 2.8e-4
rho = 958.37
lda = 0.679
cp = 4.21e3
beta_th = 0.0007504815417629351

nul=mu/rho
al=lda/(rho*cp)
g=9.81
DT=8.5
delta_th=7.14*math.pow(nul*al/(g*beta_th*DT), 1./3.)
print("delta_th=",delta_th)
