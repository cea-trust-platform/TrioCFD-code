#!/usr/bin/python
# -*- coding: utf-8 -*-

import matplotlib as mpl
mpl.use('Agg')
from numpy import *
from pylab import *
import os, glob, sys
import DNSTools as dns
import numpy as np

coord=5.62500000e-04-0.0026875

schemes = ["EULER", "RK3"]
#niv=np.array([100, 50, 20, 10, 5, 1])
#Nx = np.array([5,10,25,50,100,500])
niv=np.array([100, 50, 20, 10, 5])
Nx = np.array([5,10,25,50,100])
errs = ones((len(niv),len(schemes)))

f = figure()
ss = ['r-', 'b-']
for j, sch in enumerate(schemes):
   for i, raff in enumerate(niv):
      fic = os.path.join(sch,"Niv%d"%raff,"calcul_PRESSION_S.son")
      positions, nsondes = dns.getSondesCoords(fic)
      mat = loadtxt(fic)
      time = mat[:,0]
      for isonde in [3]:
         xp0, yp0, zp0 = positions[1,:]	
         xp, yp, zp = positions[isonde,:]
         print 'sonde en', zp, 'premier point en ', zp0
         S=1e10
         sol_ana = [S*t*t*t*(zp-zp0) for t in time]
	 sol_reel= abs(mat[:,1+isonde]-mat[:,2])
         err_rel = abs((sol_reel - sol_ana)/sol_ana)
	 err_abs = abs(sol_reel - sol_ana)
      errs[i,j] = err_rel[-1]
      pass
   loglog(Nx, errs[:,j], "x", label=r"%s"%(sch))
   pass

grid('on')
legend(loc=0)
m = errs.max()
loglog(Nx, 5.*m/Nx, 'k-.', label="Ordre -1")
loglog(Nx, 25.*m/(Nx*Nx), 'k--',  label="Ordre -2")
loglog(Nx, 125.*m/(Nx*Nx*Nx), 'k-', label="Ordre -3")
ylabel('Erreur relative sur le gradient de pression')
xlabel('Nt (s)')
title('Erreur relative en fonction du nombre de pas de temps')
print("avant savefig") 
savefig("cvg.png")

