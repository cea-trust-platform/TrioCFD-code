# -*- coding: utf-8 -*-
#!/bin/python
##############################################################################
# mode d'emploi : python accplot.py figname tdeb sauvfile datafile 
#                            0        1       2        3     4
# /!\ 3 : sauvfile -> Bien donner le FILE, pas le DIRECTORY
##############################################################################
import os, sys, glob, math
import DNSTools3 as dtool
import numpy as np
import matplotlib.pyplot as plt

import readline, rlcompleter
readline.parse_and_bind('tab:complete')

import subprocess

# regle un pb de mpl pour enregistrer des figures trop lourdes
plt.rcParams['agg.path.chunksize'] = 100000

print("Provide a name for headers... as in defo180_Eo2.34_Db0213 or defo250_Eo2.34_Db03")
tdeb=0.
tmax=5.e6
jdd=glob.glob('*.data')[-1]

########################################################################
## Recuperation des méta arguments et fichiers 
print(len(sys.argv))
if (len(sys.argv) >=2):
    if (len(sys.argv) >=3):
        # Temps initial
        tdeb=float(sys.argv[2])
        if (len(sys.argv) >= 4):
            # Fichier .sauv
            repr_file=sys.argv[3]    
            if (len(sys.argv) == 5):
                # Fichier .data
                print("#######",sys.argv[4],"#########")
                jdd=sys.argv[4]
        else:
            repr_file=dtool.getParam(jdd, 'nom_reprise', string=True)
            print("repr_file : %s [read]"%repr_file)
            

else:
    raise Exception("usage: Missing header string to name the plots")

print("Fichier .sauv : "+repr_file)


# Nom des figures
fig=sys.argv[1]
print("provided header for name_fig is %s."%fig)
head=jdd.replace('.data',"")

print("Selected jdd : %s."%jdd)

fic="OUT/%s_acceleration.out"%head
lis1=[fic]
print("File list to read: ", lis1)
########################################################################
rhol, rhov, mul, muv, alv, _, sigma, Lz, rhom, Eo, g, Svx = \
        dtool.get_prop(jdd=jdd, repr_file=repr_file, Ret=0.)
#
try:
    S0 = dtool.getParam(jdd, "terme_force_init")
except:
    S0 = 1.
try:
    vb = dtool.getParam(jdd, "vol_bulle_monodisperse")
except:
    vb = 1.
db = 2.*pow(3./(4*math.pi)*vb,1./3.)
#
try:
    facv=1./(alv*(rhol-rhov))
    facl=1./((1.-alv)*(rhol-rhov))
except:
    facv=0.
    facl=1.
#
print("rhom = ", rhom)
print(lis1)
########################################################################
for fic in lis1:
    # Lecture du fichier out
    matrice=np.genfromtxt(fic, usecols=(0,1,2,3,4,5,6,7))
    time_step = np.genfromtxt(fic, usecols=(0))
    time = np.genfromtxt(fic, usecols=(1))
    v_moy = np.genfromtxt(fic, usecols=(2))
    rhov_moy = np.genfromtxt(fic, usecols=(3))
    tau_w = np.genfromtxt(fic, usecols=(4))
    acc_v = np.genfromtxt(fic, usecols=(5))
    terme_source = np.genfromtxt(fic, usecols=(7))
    ####################################################################
    
    name = os.path.split(fic)[0].strip('./')
    # Temps initial et final 
    tdeb = max(tdeb,time[0])
    tmaxx = min(tmax,time[-1])
    itdeb = np.argmin(abs(time-tdeb))
    itmax = np.argmin(abs(time-tmaxx))
    print("iterations selected: ", itdeb, " : " , itmax)
    x = [time[itdeb],time[itmax]]
    nlast = min(1000,len(time)-1)
    ####################################################################
    
    # Ajustement des tableaux
    time = time[itdeb:itmax]
    v_moy = v_moy[itdeb:itmax]
    rhov_moy = rhov_moy[itdeb:itmax]
    tau_w = tau_w[itdeb:itmax]
    acc_v = acc_v[itdeb:itmax]
    terme_source = terme_source[itdeb:itmax]
    ####################################################################
    
    # Tracer Re_tau
    fig1, ax1 = plt.subplots(1)
    Re_tau = rhol*np.sqrt(abs(tau_w)/rhol)*Lz/2./mul 
    ax1.plot(time, Re_tau, label=fic)
    ##########
    
    # Tracer terme_source
    fig2, ax2 = plt.subplots(1)
    ax2.plot([time[0],time[-1]], [S0,S0], "--", label="Sjdd")
    ax2.plot(time, terme_source, label="S")
    moyenne_source = terme_source[-nlast:].mean()
    erreur_source = moyenne_source - S0
    print("Mean S last %d step : %g"%(nlast,terme_source[:-nlast].mean()))
    print("Error on S (<S>-S_theo) : %g"%(erreur_source))
    ax2.legend(loc=0)
    ##########
    
    # Tracer vitesses
    fig3, ax3 = plt.subplots(1)
    ul=facl*(rhov_moy-rhov*v_moy) # =0 si on est dans vap., =1 si on est dans liq.
    uv=facv*(rhol*v_moy-rhov_moy) # =1 si on est dans vap., =0 si on est dans liq.
    ax3.plot(time, ul, '-', label='ul')
    ax3.plot(time, uv, '-', label='uv')
    uv_moy = uv[-nlast:].mean()
    ul_moy = ul[-nlast:].mean()
    ur_moy = uv_moy -ul_moy
    ax3.plot([time[-nlast],time[-1]], [uv_moy, uv_moy], '--', label='uvm[-%d:] = %f'%(nlast,uv_moy))
    print("uv-last %d: %f"%(nlast, uv_moy))
    print("ur = %g"%(ur_moy))
    print("Reb = %g"%(rhol*ur_moy*db/mul))
    print("rhom*uv[-%d:] = %f"%(nlast,rhom*uv_moy))
    ax3.legend(loc=0)
    ##########
    
    # Tracer rho*vitesse
    fig4, ax4 = plt.subplots(1)
    debit_l = (1.-alv)*rhol*ul
    debit_v = alv*rhov*uv
    debit_t = debit_l+debit_v
    ax4.plot(time, debit_l, '-', label=r'$\alpha_l\,\rho_l\,u_l$')
    ax4.plot(time, debit_v, '-', label=r'$\alpha_v\,\rho_v\,u_v$')
    print("debl [min:max] = ",  debit_l.min(),  debit_l.max())
    print("debv [min:max] = ",  debit_v.min(),  debit_v.max())
    print("debt [min:max] = ",  debit_t.min(),  debit_t.max())
    ax4.plot(time, debit_t, '-', label=r'$\rho\,u$')
    ax4.set_xlabel(r'$t$ [s]')
    ax4.set_ylabel(r'$\rho\,u\:\mathrm{[kg.m}^{-2}\mathrm{.s}^{-1}]$')
    ax4.legend(loc=0)
    ##########
    
    # Tracer tau_w
    fig5, ax5 = plt.subplots(1)
    ax5.plot(time, tau_w, '-', label='tauw')
    ax5.legend(loc=0)
    #import pdb; pdb.set_trace()
    ##########
    
    # vitesse relative
    fig6, ax6 = plt.subplots(1)
    ax6.plot(time, uv-ul, '-', label='ur')
    ax6.plot(time[0],(uv-ul)[0],'or')
    ax6.text(time[itdeb],(uv-ul)[itdeb],'t = %s'%(time[itdeb]),fontsize=22)
    ax6.legend(loc=0)
    ##########
    pass

print("fig retau")
ax1.grid(True)
fig1.savefig(fig+'_Retau.png')
plt.close()
print("fig source")
ax2.grid(True)
fig2.savefig(fig+'_S.png')
plt.close()
print("fig vitesse")
ax3.grid(True)
fig3.savefig(fig+'_ul_uv.png')
plt.close()
print("debit")
ax4.grid(True)
fig4.savefig(fig+'_debits.png')
plt.close()
print("tauw")
ax5.grid(True)
fig5.savefig(fig+'_tauw.png')
plt.close()
print("relative")
ax6.grid(True)
fig6.savefig(fig+'_ur.png')
plt.close()

fic="OUT/%s_bulles_external_force_every_0.out"%head
if os.path.isfile(fic):
    mat=np.genfromtxt(fic)
    tdeb=float(sys.argv[2])
    tdeb=max(tdeb,mat[0,1])
    tmaxx=min(tmax,mat[-1,0])
    itdeb=np.argmin(abs(mat[:,0]-tdeb)) 
    t=mat[itdeb:,0]
    Fmean=mat[:,0:].mean(axis=1)
    fig7, ax7 = plt.subplots(1)
    ax7.plot(t,Fmean[itdeb:], '-', label='Fmean')
    ax7.legend(loc=0)
    ax7.grid(True)
    fig7.savefig(fig+'_force.png')
    plt.close()
    n=min(len(t),1000)
    print("Mean force last %d step : %g"%(n,Fmean[:-n].mean()))
    pass

fic="OUT/%s_bulles_centre_x.out"%head
if os.path.isfile(fic):
  for w in ["x", "y", "z"]:
    mat=np.genfromtxt(fic.replace("_x.out","_%s.out"%w))
    #mat=np.genfromtxt(fic)
    tdeb=float(sys.argv[2])
    tdeb=max(tdeb,mat[0,0])
    tmaxx=min(tmax,mat[-1,0])
    itdeb=np.argmin(abs(mat [:,0]-tdeb)) 
    t=mat[itdeb:,0]
    xmean=mat[:,1:].mean(axis=1)
    fig7, ax7 = plt.subplots(1)
    ax7.plot(t,xmean[itdeb:], '-', label='%smean'%w)
    ax7.legend(loc=0)
    ax7.grid(True)
    fig7.savefig(fig+'_%sb.png'%w)
    plt.close()
    n=min(len(t),1000)
    print("Mean %s position bubble last %d step : %g"%(w,n,xmean[:-(n-1)].mean()))
    pass
  pass

exit()

