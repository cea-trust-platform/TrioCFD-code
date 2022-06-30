# -*- coding: utf-8 -*-

import matplotlib as mpl
mpl.use('Agg')
from DNSTools import *

head = getJddName()
jdd = head+".data"

rhol = getValue("rho_liquide", jdd)
mul = getValue("mu_liquide", jdd)
Lz = getValue("uniform_domain_size_k", jdd) 
h =  Lz / 2.

# Recuperation des donnees moyennes : 
tm, z, dvar, moys, ltintegration  = getStatsMoy()
ti, _, _, insts, _  = getStatsInstant()

# En instantane, il y a un dt de plus, le premier... 
# ti[1:] == tm

nz  = len(z)
ntm = len(tm)
dz = Lz/nz
iv = dvar["UI"]
tmax = tm[-1]

Um = zeros(ntm)
# figure()
for it, t in enumerate(tm[::5]):
   Umoy = moys[it, :,iv]
   # plot(z, Umoy, "%s"%(1.-t/tmax), label="moy : t = %g" %t)
   # Vitesse debitante : 
   Um[it] = Umoy.mean() # car dz = cste.
   pass

figure()
# Ne trace qu'un sur 10
for it, t in enumerate(ti[::5]):
   Uins = insts[it, :,iv]
   plot(z, Uins, "%s"%(1.-t/tmax), label="inst : t = %g" %t)
   pass

title("Vitesse instantannee")
savefig(head+"_Uinst.png")
close()
clf()

Ub = insts[:, :,iv].mean(axis=1)
Um = moys[:, :,iv].mean(axis=1)

#Reynolds canal : 
Reb_inst = rhol*Ub*2*Lz/mul
Reb_moy = rhol*Um*2*Lz/mul

figure()
subplot(221)
title("Vitesse debitante")
plot(tm, Um, label="moy")
plot(ti, Ub, label="inst")
legend(loc = 0)
grid('on')


subplot(222)
title("Re canal")
plot(tm, Reb_moy, label="moy")
plot(ti, Reb_inst, label="inst")
legend(loc = 0)
grid('on')


subplot(223)
title("Re tau")
Ret_moy, _, _ = evaluateRetau(moys[:, :,iv], z, jdd)
Ret_inst, _, _ = evaluateRetau(insts[:, :,iv],  z, jdd)

plot(tm, Ret_moy, label="moy")
plot(ti, Ret_inst, label="inst")
legend(loc = 0)
grid('on')

subplot(224)
upup_m = moys[:, :,dvar['UUI']] - moys[:, :,dvar['UI']]*moys[:, :,dvar['UI']]
vpvp_m = moys[:, :,dvar['VVI']] - moys[:, :,dvar['VI']]*moys[:, :,dvar['VI']]
wpwp_m = moys[:, :,dvar['WWI']] - moys[:, :,dvar['WI']]*moys[:, :,dvar['WI']]

plot(z,upup_m[-1,:], label="uu")
plot(z,vpvp_m[-1,:], label="vv")
plot(z,wpwp_m[-1,:], "k--", label="ww")
legend(loc = 0)
grid('on')

savefig(head+"_multi.png")
close("all")

print("Udebitant : ", Um.mean(), " ou ", Ub.mean())
print("Reb       : ", Reb_moy.mean(), " ou ", Reb_inst.mean())
print("Retau     : ", Ret_moy.mean(), " ou ", Ret_inst.mean())
