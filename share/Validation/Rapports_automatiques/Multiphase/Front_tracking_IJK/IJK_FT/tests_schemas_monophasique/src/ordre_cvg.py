#!/usr/bin/python
# -*- coding: utf-8 -*-

import matplotlib as mpl
mpl.use('Agg')
from numpy import *
from pylab import *
import os, glob, sys
from DNSTools import *
from math import *

# Coordonnees de la sonde : 
xp = 5.86025000e-03  # ix = 5? 
yp = 2.930125e-3 # jy = 2.5
zp = 0.0014650625/2.
Lx = 0.093764
A = 2*pi / Lx

cd("ORDRE_T")
schemes = ["EULER", "RK3"]
Nx = array([1, 2, 4, 8, 16])

errs = ones((len(Nx),len(schemes)))
cwd = pwd()

f = figure()
ss = ['r-', 'g-', 'b-.', 'k-.', 'm-.', 'c-.']
ss+=ss
for j, sch in enumerate(schemes):
   for i, raff in enumerate(Nx):
      cd(os.path.join(sch,"DT%.2d"%raff))
      head=getJddName()
      fic = head+"_P.son"
      positions, nsondes = getSondesCoords(fic)
      mat = loadtxt(fic)
      time = mat[:,0]
      for isonde in range(nsondes):
         xp, yp, zp = positions[isonde,:]
         sol_ana = [(2+2*sin((xp - yp)*A-2*A*t)) for t in time]
         err_abs = (mat[:,1+isonde] - sol_ana) #/ sol_ana
      subplot(221)
      plot(time, err_abs, ss[isonde+nsondes*(2*i+j)], label="Scheme %s -- DT = %d -- Sonde %d"%(sch,raff, isonde))
      subplot(223)
      if i+j == 0:
         plot(time, sol_ana, 'k-o', ms=2, label="Solution analytique")
      plot(time, mat[:,1+isonde], label="%s -- DT = %d -- Sonde %d"%(sch,raff, isonde))
      #
      # Store error at final time for last sonde : 
      errs[i,j] = err_abs[-1]
      cd(cwd)
      pass
   subplot(222)
   loglog(Nx, errs[:,j], ss[j]+"x", label="%s"%(sch))
   pass

grid('on')
legend(loc=0)
m = errs.max()
loglog(Nx, m/Nx, 'k-.', label="Ordre -1")
loglog(Nx, m/(Nx*Nx), 'k--',  label="Ordre -2")
loglog(Nx, m/(Nx*Nx*Nx), 'k-', label="Ordre -3")
ylabel('Erreur relative sur U')
xlabel('t (s)')

subplot(221)
title('erreur absolue')
grid('on')
legend(loc=0)
#
subplot(223)
title('Comparaison solution analytique a solution numerique')
grid('on')
legend(loc=0)

show()
f.set_size_inches([8,4])
savefig("cvg.png")

