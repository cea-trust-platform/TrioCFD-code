#!/usr/bin/python
# -*- coding: utf-8 -*-

import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.use('Agg')
import numpy as np
import commons.DNSTools as dns

head = dns.getJddName()
jdd = head+".data"

rhol = dns.getValue("rho_liquide", jdd)
mul = dns.getValue("mu_liquide", jdd)
Lz = dns.getValue("uniform_domain_size_k", jdd) 
h =  Lz / 2.

# Recuperation des donnees moyennes : 
tm, z, dvar, moys, ltintegration  = dns.getStatsMoy()
ti, _, _, insts, _  = dns.getStatsInstant()

# En instantane, il y a un dt de plus, le premier... 
# ti[1:] == tm

nz  = len(z)
ntm = len(tm)
dz = Lz/nz
iv = dvar["UI"]
tmax = tm[-1]

Um = np.zeros(ntm)
# figure()
for it, t in enumerate(tm[::5]):
   Umoy = moys[it, :,iv]
   # plot(z, Umoy, "%s"%(1.-t/tmax), label="moy : t = %g" %t)
   # Vitesse debitante : 
   Um[it] = Umoy.mean() # car dz = cste.
   pass

plt.figure()
# Ne trace qu'un sur 10
for it, t in enumerate(ti[::5]):
   Uins = insts[it, :,iv]
   plt.plot(z, Uins, "%s"%(1.-t/tmax), label="inst : t = %g" %t)
   pass

plt.title("Vitesse instantannee")
plt.savefig(head+"_Uinst.png")
plt.close()

Ub = insts[:, :,iv].mean(axis=1)
Um = moys[:, :,iv].mean(axis=1)

#Reynolds canal : 
Reb_inst = rhol*Ub*2*Lz/mul
Reb_moy = rhol*Um*2*Lz/mul

plt.figure()
plt.subplot(221)
plt.title("Vitesse debitante")
plt.plot(tm, Um, label="moy")
plt.plot(ti, Ub, label="inst")
plt.legend(loc = 0)
plt.grid('on')

plt.subplot(222)
plt.title("Re canal")
plt.plot(tm, Reb_moy, label="moy")
plt.plot(ti, Reb_inst, label="inst")
plt.legend(loc = 0)
plt.grid('on')

plt.subplot(223)
plt.title("Re tau")
Ret_moy, _, _ = dns.evaluateRetau(moys[:, :,iv], z, jdd)
Ret_inst, _, _ = dns.evaluateRetau(insts[:, :,iv],  z, jdd)

plt.plot(tm, Ret_moy, label="moy")
plt.plot(ti, Ret_inst, label="inst")
plt.legend(loc = 0)
plt.grid('on')

plt.subplot(224)
upup_m = moys[:, :,dvar['UUI']] - moys[:, :,dvar['UI']]*moys[:, :,dvar['UI']]
vpvp_m = moys[:, :,dvar['VVI']] - moys[:, :,dvar['VI']]*moys[:, :,dvar['VI']]
wpwp_m = moys[:, :,dvar['WWI']] - moys[:, :,dvar['WI']]*moys[:, :,dvar['WI']]

plt.plot(z,upup_m[-1,:], label="uu")
plt.plot(z,vpvp_m[-1,:], label="vv")
plt.plot(z,wpwp_m[-1,:], "k--", label="ww")
plt.legend(loc = 0)
plt.grid('on')

plt.savefig(head+"_multi.png")
plt.close("all")

print("Udebitant : ", Um.mean(), " ou ", Ub.mean())
print("Reb       : ", Reb_moy.mean(), " ou ", Reb_inst.mean())
print("Retau     : ", Ret_moy.mean(), " ou ", Ret_inst.mean())
