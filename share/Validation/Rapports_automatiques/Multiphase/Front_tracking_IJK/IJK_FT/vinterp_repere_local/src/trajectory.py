#!/usr/bin/python
# -*- coding: utf-8 -*-
import os, glob, sys
import matplotlib as mpl
mpl.use('Agg') 
import matplotlib.pyplot as plt
from commons.DNSTools import *

f = plt.figure()
for j, fold in enumerate(["EULER", "RK3"]): #, "RK3_fs1.20"]):
   coords, t, nb = getBarys(os.path.join(fold,"vinterp_repere_local_avec_remaillage"))
   #
   ss = ['r-', 'g-', 'b-.', 'k-.', 'm-.', 'c-.']
   #for ib in range(nb):
   for ib in [1]: # seulement la bulle 1. 
      # plot(xb, yb, ss[j], label="Bulle %d -- Scheme %s"%(ib,fold))
      xb = coords[0,:,ib]
      yb = coords[1,:,ib]
      plt.subplot(121)
      plt.plot(xb, yb, ss[j], label="Scheme %s"%(fold))
      plt.grid('on')
      plt.subplot(122)
      plt.grid('on')
      plt.plot(xb, yb, ss[j], label="Scheme %s"%(fold))
      plt.xlim([-0.01, 0.11])
      plt.ylim([-0.11, -0.01])
      pass
   pass

plt.legend(loc=0)
#show()
f.set_size_inches([8,4])
plt.savefig("trajectory.png")
