# -*- coding: utf-8 -*-
#!/bin/python

import os, sys, glob, math
import DNSTools3 as dtool
import numpy as np
from matplotlib.pyplot import *

import readline, rlcompleter
readline.parse_and_bind('tab:complete')

import subprocess

# regle un pb de mpl pour enregistrer des figures trop lourdes
matplotlib.rcParams['agg.path.chunksize'] = 100000

print("Provide a name for headers... as in defo180_Eo2.34_Db0213 or defo250_Eo2.34_Db03")
tdeb=0.
tmax=5.e6
jdd=glob.glob('*.data')[0]
print "jdd selected : ", jdd

if (len(sys.argv) >=2):
    if (len(sys.argv) >=3):
        tdeb=float(sys.argv[2])
        if (len(sys.argv) == 4):
            repr_file=sys.argv[3]    
            if (not os.path.isfile(repr_file)):
                raise Exception("usage: Provided file %s does not exist"%repr_file)
            else:
                print("repr_file= %s [given]"%repr_file)
        else:
            repr_file=dtool.getParam(jdd, 'nom_reprise', string=True)
            print("repr_file= %s [read]"%repr_file)
else:
    raise Exception("usage: Missing header string to name the plots")

fig=sys.argv[1]
print("provided header for name_fig is %s."%fig)
head=jdd.replace('.data',"")
print("Selected jdd : %s."%jdd)

fic="%s_acceleration.out"%head
# Ici on récupère tous les chemins des ..._acceleration.out dans des dossiers parents
# fics=glob.glob(fic)
# p=subprocess.Popen("find . -maxdepth 15 -name "+fic, stdout=subprocess.PIPE, shell=True)
# lis1=p.stdout.read().split()
#HACK: ne retient que le premier fichier en fait 
lis1=[fic]
print("File list to read: ", lis1)

rhol, rhov, mul, muv, alv, _, sigma, Lz, rhom, Eo, g, Svx = \
        dtool.get_prop(jdd=jdd, repr_file=repr_file, Ret=0.)

S0 = dtool.getParam(jdd, "terme_force_init")
vb = dtool.getParam(jdd, "vol_bulle_monodisperse")
db = 2.*pow(3./(4*math.pi)*vb,1./3.)
try:
    facv=1./(alv*(rhol-rhov))
    facl=1./((1.-alv)*(rhol-rhov))
except:
    facv=0.
    facl=1.

print("rhom = ", rhom)
# print("tstep\ttime\tVx\trhoVx\ttauw\tda/dt\tNewT\tacceleration")
print(lis1)

for fic in lis1:
    mat=np.genfromtxt(fic, usecols=(0,1,2,3,4,5,6,7))
    name=os.path.split(fic)[0].strip('./')
    tdeb=max(tdeb,mat[0,1])
    tmaxx=min(tmax,mat[-1,1])
    itdeb=np.argmin(abs(mat[:,1]-tdeb))
    itmax=np.argmin(abs(mat[:,1]-tmaxx))
    print("iterations selected: ", itdeb, " : " , itmax)
    x=[mat[itdeb,1],mat[itmax,1]]
    t=mat[itdeb:itmax,1]
    nlast=min(1000,len(t))
    #
    f0 = figure(0)
    plot(t, rhol*np.sqrt(abs(mat[itdeb:itmax,4])/rhol)*Lz/2./mul, label=fic)
    ##########
    f2 = figure(2)
    # plot(x, [0.104070213285, 0.104070213285], '--', label="Ret=250")
    # plot(x, [1.00828349999999997e-01,1.00828349999999997e-01], '--', label="Ret=180")
    plot([t[0],t[-1]], [S0,S0], "--", label="Sjdd")
    S=mat[itdeb:itmax,7]
    plot(t, S, label="S") #fic)
    Sm = S[-nlast:].mean()
    errS = Sm-S0
    print("Mean S last %d step : %g"%(nlast,S[:-nlast].mean()))
    print("Error on S (<S>-S_theo) : %g"%(errS))
    legend(loc=0)
    ##########
    f3 = figure(3)
    ul=facl*(mat[:,3]-rhov*mat[:,2])
    uv=facv*(rhol*mat[:,2]-mat[:,3])
    plot(mat[itdeb:,1], ul[itdeb:], '-', label='ul')
    plot(mat[itdeb:,1], uv[itdeb:], '-', label='uv')
    uvm = uv[-nlast:].mean()
    ulm = ul[-nlast:].mean()
    urm = uvm -ulm
    plot([t[-nlast],t[-1]], [uvm, uvm], '--', label='uvm[-%d:] = %f'%(nlast,uvm))
    print("uv-last %d: %f"%(nlast, uvm))
    print("ur= %g"%(urm))
    print("Reb= %g"%(rhol*urm*db/mul))
    print("rhom*uv[-%d:] = %f"%(nlast,rhom*uvm))
    legend(loc=0)
    ##########
    f4 = figure(4)
    debl = (1.-alv)*rhol*ul
    debv = alv*rhov*uv
    debt= debl+debv
    plot(t, debl[itdeb:itmax], '-', label=r'$\alpha_l\,\rho_l\,u_l$')
    plot(t, debv[itdeb:itmax], '-', label=r'$\alpha_v\,\rho_v\,u_v$')
    print("debl [min:max]=",  debl[itdeb:itmax].min(),  debl[itdeb:itmax].max())
    print("debv [min:max]=",  debv[itdeb:itmax].min(),  debv[itdeb:itmax].max())
    print("debt [min:max]=",  debt[itdeb:itmax].min(),  debt[itdeb:itmax].max())
    plot(t, debt[itdeb:itmax], '-', label=r'$\rho\,u$')
    xlabel(r'$t$ [s]')
    ylabel(r'$\rho\,u\:\mathrm{[kg.m}^{-2}\mathrm{.s}^{-1}]$')
    legend(loc=0)
    ##########
    f5 = figure(5)
    plot(t, mat[itdeb:itmax,4], '-', label='tauw')
    legend(loc=0)
    #import pdb; pdb.set_trace()
    ##########
    f6 = figure(6)
    plot(mat[itdeb:,1], uv[itdeb:]-ul[itdeb:], '-', label='ur')
    legend(loc=0)
    pass

f0.gca().grid(True)
f0.savefig(fig+'_Retau.png')

f2.gca().grid(True)
f2.savefig(fig+'_S.png')

f3.gca().grid(True)
f3.savefig(fig+'_ul_uv.png')

f4.gca().grid(True)
f4.savefig(fig+'_debits.png')

f5.gca().grid(True)
f5.savefig(fig+'_tauw.png')

f6.gca().grid(True)
f6.savefig(fig+'_ur.png')

fic="%s_bulles_external_force_every_0.out"%head
if os.path.isfile(fic):
    mat=np.genfromtxt(fic)
    t=mat[:,0]
    Fmean=mat[:,1:].mean(axis=1)
    f7 = figure(7)
    plot(t,Fmean, '-', label='Fmean')
    legend(loc=0)
    f7.gca().grid(True)
    f7.savefig(fig+'_force.png')
    close(f7)
    n=min(len(t),1000)
    print("Mean force last %d step : %g"%(n,Fmean[:-n].mean()))
    pass

fic="%s_bulles_centre_x.out"%head
if os.path.isfile(fic):
  for w in ["x", "y", "z"]:
    mat=np.genfromtxt(fic.replace("_x.out","_%s.out"%w))
    # mat=np.genfromtxt(fic)
    t=mat[:,0]
    xmean=mat[:,1:].mean(axis=1)
    f7 = figure(7)
    plot(t,xmean, '-', label='%smean'%w)
    legend(loc=0)
    f7.gca().grid(True)
    f7.savefig(fig+'_%sb.png'%w)
    close(f7)
    n=min(len(t),1000)
    print("Mean %s position bubble last %d step : %g"%(w,n,xmean[:-n].mean()))
    pass
  pass

exit()

