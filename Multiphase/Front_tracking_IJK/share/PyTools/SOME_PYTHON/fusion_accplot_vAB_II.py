# -*- coding: utf-8 -*-
#!/bin/python

""" Trace l'evolution temporelle des vitesses, debits, Re_tau et source pour ue liste de RUN terminés normalement"""
##############################################################################
# mode d'emploi : python accplot.py figname tdeb header liste_RUNS lite_symbols
#                            0        1       2    3         4        5
# /!\ 3 : liste_RUNS : donner chemin jusqu'a RUN0X INCLUS
# /!\ 4 : header : nom du jdd, qui est aussi le mm que celui du sauv
#                  peut donner "auto",le script prend le premier .data 
#                  et le premier .sauv
# 5 : ce qui apparaitra dans la legende
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

########################################################################
## Recuperation des méta arguments et fichiers 
fig=sys.argv[1]
tdeb=float(sys.argv[2])
head=sys.argv[3]
liste_RUNS=sys.argv[4]
print(liste_RUNS)
if not(".txt" in sys.argv[4]):
    # donne la listes dans les arguments, sous forme [RUN01,RUN02,...,RUN0X]
    liste_RUNS=liste_RUNS.rstrip("]").lstrip("[").split(",")
else:
    # donne la liste dans un fichier, fait avec ls ...*... > nom_fichier
    f=open(liste_RUNS)
    liste=f.readlines()
    f.close()
    liste_RUNS=liste
    for i,l in enumerate(liste_RUNS):
        liste_RUNS[i]=liste_RUNS[i].rstrip("\n")
print(liste_RUNS)

liste_symb=sys.argv[5]
print(liste_symb)
if not(".txt" in sys.argv[5]):
    # donne la listes dans les arguments, sous forme [legend1,legend2,...,legendeX]
    liste_symb=liste_symb.rstrip("]").lstrip("[").split(";")
else:
    # donne la liste dans un fichier, fait avec ls ...*... > nom_fichier
    f=open(liste_symb)
    liste=f.readlines()
    f.close()
    liste_symb=liste
    for i,l in enumerate(liste_symb):
        liste_symb[i]=liste_symb[i].rstrip("\n")
print(liste_symb)

# Preparation des graphiques
print("provided header for name_fig is %s."%fig)
# Tracer Re_tau
fig1, ax1 = plt.subplots(1,figsize=(12,6))
ax1.set_title(r"$Re_{\tau} =  $")
# Tracer terme_source
fig2, ax2 = plt.subplots(1,figsize=(12,6))
ax2.set_title("Source")
# Tracer vitesses
fig3, ax3 = plt.subplots(1,figsize=(12,6))
ax3.set_title(r'$u_l = \frac{\langle \rho u \rangle - \rho_v \langle u \rangle}{\alpha_l (\rho_l - \rho_v)}$'+", "+
              r'$u_v = \frac{\rho_l \langle u \rangle - \langle \rho u \rangle}{\alpha_v (\rho_l - \rho_v)}$')
# Tracer rho*vitesse
fig4, ax4 = plt.subplots(1,figsize=(12,6))
ax4.set_title(r'$\alpha_l\,\rho_l\,u_l$'+", "+r'$\alpha_v\,\rho_v\,u_v$'+", "+r'$\rho\,u$')
# Tracer Re_b
fig5, ax5 = plt.subplots(1,figsize=(12,6))
ax5.set_title(r"$Re_{b} = \frac{\rho_l * u_r * d_b}{\nu_l} $")
# vitesse relative
fig6, ax6 = plt.subplots(1,figsize=(12,6))
ax6.set_title(r'$u_r = u_v - u_l$')
########################################################################


for n_run,run in enumerate(liste_RUNS):
    run.rstrip("\n")
    symbol=liste_symb[n_run]+'\n'
    print(" ")
    print("#############################################")
    print("###### "+run+" ###########")
    
    repr_file = glob.glob(run+"/"+"NEXT/"+"*sauv")[0]
    jdd = glob.glob(run+"/"+"*data")[0]
    head=jdd.replace('.data',"")
    
    print("Selected jdd : %s."%jdd)
    print(dtool.getParam(jdd, "terme_force_init"))
    try: # expression_derivee_acceleration_
        # Le fichier acceleration n'est ecrit que si la derivee de la force n'est pas nulle...
        out_file=glob.glob(run+"/"+"OUT/"+"*acceleration.out")[0]
        matrice=np.genfromtxt(out_file, usecols=(0,1,2,3,4,5,6,7))
    except:
        print("Le fichier acceleration.out n'est pas trouve.")
        print("Le fichier check_etapes_et_termes existe ? (1:oui, 0:non)", dtool.getFlag(jdd, "test_etapes_et_bilan"))
    if dtool.getFlag(jdd, "test_etapes_et_bilan")!=0:
        # Si on n'a pas le fichier acceleration, on va devoir utiliser le fichier check etapes et temres
        print("#########################################")
        print("Le fichier acceleration out n'existe pas (d_t S = 0). Mais le fichier check_etapes_et_termes existe.")
        print("use_check_etapes_et_termes.py")
        print("#########################################")
        arguments = "a completer"
        print("Execution de : python AJOUTER LES COMMANDES POUR EXECUTER use_check_etapes_et_termes.py")
        # os.system('pyhton '+arguments)
    lis1=[out_file]
    print("File list to read: ", lis1)
    ########################################################################
    rhol, rhov, mul, muv, alv, nb_bbl, sigma, Lz, rhom, Eo, g, Svx = \
            dtool.get_prop(jdd=jdd, repr_file=repr_file, Ret=0.)
    
    S0 = dtool.getParam(jdd, "terme_force_init")
    
    # En monophasique, le volume de bulle et le diametre de bulle n'ont pas de sens
    print("disable_diphasique : ",dtool.getFlag(jdd, "disable_diphasique"))
    if not dtool.getFlag(jdd, "disable_diphasique"):
        vb = dtool.getParam(jdd, "vol_bulle_monodisperse")
        db = 2.*pow(3./(4*math.pi)*vb,1./3.)
        print("Now, vb = ",vb, db)
    else:
        vb = 0.
        db = 0.
    try:
        facv=1./(alv*(rhol-rhov))
        facl=1./((1.-alv)*(rhol-rhov))
    except:
        facv=0.
        facl=1.
    
    print("rhom = ", rhom)
    print(lis1)
    ########################################################################
    for out_file in lis1:
        # couleur pour tous les plots
        color3 = next(ax3._get_lines.prop_cycler)['color']
        
        # Lecture du fichier out
        matrice=np.genfromtxt(out_file, usecols=(0,1,2,3,4,5,6,7))
        time_step = np.genfromtxt(out_file, usecols=(0))
        time = np.genfromtxt(out_file, usecols=(1))
        v_moy = np.genfromtxt(out_file, usecols=(2))
        rhov_moy = np.genfromtxt(out_file, usecols=(3))    # rho*u. /!\ rhov = \rho_v
        tau_w = np.genfromtxt(out_file, usecols=(4))
        acc_v = np.genfromtxt(out_file, usecols=(5))
        terme_source = np.genfromtxt(out_file, usecols=(7))
        ####################################################################
        
        name = os.path.split(out_file)[0].strip('./')
        # Temps initial et final 
        tdeb = max(tdeb,time[0])
        tmaxx = min(tmax,time[-1])
        itdeb = np.argmin(abs(time-tdeb))
        itmax = np.argmin(abs(time-tmaxx))
        print("iterations selected: ", itdeb, " : " , itmax)
        x = [time[itdeb],time[itmax]]
        nlast = int(len(time)/2.)# min(1000,len(time)-1)
        ####################################################################
        
        # Ajustement des tableaux
        t0 = time[itdeb]          # on garde la valeur du 1er instant, pour ladonner qq part
        time = time[itdeb:itmax]-time[itdeb]
        t0 = np.round(t0,2)             
        v_moy = v_moy[itdeb:itmax]
        rhov_moy = rhov_moy[itdeb:itmax]
        tau_w = tau_w[itdeb:itmax]
        acc_v = acc_v[itdeb:itmax]
        terme_source = terme_source[itdeb:itmax]
        ####################################################################
        
        # Adoucissement eventuel
        N_smooth = 80
        filtre = np.ones(N_smooth)/float(N_smooth)
        ####################################################################
        
        # Adimensionnement
        echelle_vit = 1.
        echelle_lng = 1.
        echelle_tps = 1.
        if "adim_d" in fig:
            """ Adimentionnement avec diametre de bulle """
            echelle_vit = math.sqrt(abs(g*db))
            echelle_lng = db
            echelle_tps = echelle_lng / echelle_vit
        elif "adim_e" in fig:
            """ 
            Adimentionnement avec la distance entre 2 bulles
                -> a affiner, dist min, moy, max ? 
                -> dist. selon x, selon (y,z) ?
            """
            echelle_lng = (Lz / nb_bbl)**(1./3.)
            echelle_vit = math.sqrt(abs(g*echelle_lng))
            echelle_tps = echelle_lng / echelle_vit 
        print("echelle_lng : ",echelle_lng)
        print("echelle_vit : ",echelle_vit)
        print("echelle_tps : ",echelle_tps)
        ####################################################################
        
        # Construction des grandeurs
        # Re_tau
        Re_tau = rhol*np.sqrt(abs(tau_w)/rhol)*Lz/2./mul 
        # source
        
        # Vitesses
        ul = facl*(rhov_moy-rhov*v_moy) # =0 si on est dans vap., =1 si on est dans liq.
        uv = facv*(rhol*v_moy-rhov_moy) # =1 si on est dans vap., =0 si on est dans liq.
        ur = uv-ul
        # Re_b
        print("rho_l, db, mu_l",rhol,db,mul)
        Re_b = rhol*ur*db / mul
        print("Reb[-1] : ",Re_b[-1])
        # debits
        debit_l = (1.-alv)*rhol*ul
        debit_v = alv*rhov*uv     
        debit_t = debit_l+debit_v 
        # adim
        if "adim" in fig:
            time /= echelle_tps
            ul /= echelle_vit
            uv /= echelle_vit
            ur /= echelle_vit
            print(Re_b[-1])
        ####################################################################
        
        # la longeur du plus petit vecteur dicte combien de points de trace on a 
        n_max = min(len(time),len(np.convolve( ul,filtre,mode='same')))
        
        # Trace
        # Re_tau
        ax1.plot(time, Re_tau, label=out_file)
        # Re_b
        ax5.plot(time, Re_b, label=symbol)
        # terme source
        ax2.plot([time[0],time[-1]], [S0,S0], "--")
        ax2.plot(time, terme_source, label=symbol+r' $t_0=$'+str(t0))
        moyenne_source = terme_source[-nlast:].mean()
        erreur_source = moyenne_source - S0
        print("Mean S last %d step : %g"%(nlast,terme_source[:-nlast].mean()))
        print("Error on S (<S>-S_theo) : %g"%(erreur_source))
        
        # vitesses
        if (n_max-N_smooth > 0):
            ax3.plot(time[N_smooth:n_max-N_smooth], np.convolve(ul,filtre,mode='same')[N_smooth:n_max-N_smooth], '-',                color=color3, label="liq. , "+symbol+r' $t_0=$'+str(t0))
            ax3.plot(time[N_smooth:n_max-N_smooth], np.convolve(uv,filtre,mode='same')[N_smooth:n_max-N_smooth], linestyle='dotted', color=color3)
        else:
            ax3.plot(time[:], ul[:], '-',                color=color3, label=symbol+r' $t_0=$'+str(t0))
            ax3.plot(time[:], uv[:], linestyle='dotted', color=color3)
        uv_moy = uv[-nlast:].mean()
        ul_moy = ul[-nlast:].mean()
        ur_moy = uv_moy -ul_moy
        print("nlast : ",nlast)
        # ax3.plot([time[-nlast],time[-1]], [uv_moy, uv_moy], '--', label='uvm[-%d:] = %f'%(nlast,uv_moy))
        print("uv-last %d: %f"%(nlast, uv_moy))
        print("uv = %g"%(uv_moy))
        print("ur = %g"%(ur_moy))
        print("Reb = %g"%(rhol*ur_moy*db/mul))
        print("rhom*uv[-%d:] = %f"%(nlast,rhom*uv_moy))
        
        # debits
        ax4.plot(time, debit_l, '-', label=r'liq., '+r'%s'%(symbol)+r' $t_0=$'+str(t0))
        ax4.plot(time, debit_v, linestyle='dotted')
        ax4.plot(time, debit_t, '.-')
        print("debl [min:max] = ",  debit_l.min(),  debit_l.max())
        print("debv [min:max] = ",  debit_v.min(),  debit_v.max())
        print("debt [min:max] = ",  debit_t.min(),  debit_t.max())
        
        # Vitesse relative
        ax6.plot(time, ur, '-', color=color3, label=r'$u_r$, '+symbol)
        ax6.plot(time[0],ur[0],'or')
        
        # Reb a partir de la vitesse relative
        # ax6.text(time[itdeb],(uv-ul)[itdeb],'t = %s'%(time[itdeb]),fontsize=22)
        
        
        

# terme_source
ax2.legend(loc=0)
##########
        
# Vitesse
ax3.legend(loc='center left', bbox_to_anchor=(1.04, 0.5),fontsize=18,framealpha=0)
for label in (ax3.get_xticklabels() + ax3.get_yticklabels()):
    label.set_fontsize(18)
##########
        

# Debit = u*rho
ax4.set_xlabel(r'$t$ [s]')
ax4.set_ylabel(r'$\rho\,u\:\mathrm{[kg.m}^{-2}\mathrm{.s}^{-1}]$')
ax4.legend(loc='center left', bbox_to_anchor=(1.04, 0.5),fontsize=18,framealpha=0)
for label in (ax4.get_xticklabels() + ax4.get_yticklabels()):
    label.set_fontsize(18)  
##########
        
# Re_b = rho_l * u_r * db / mu_l
ax5.legend(loc='center left', bbox_to_anchor=(1.04, 0.5),fontsize=18,framealpha=0)
for label in (ax5.get_xticklabels() + ax6.get_yticklabels()):
    label.set_fontsize(18)
##########
    
# vitesse relative, ur = ul-uv
ax6.legend(loc='center left', bbox_to_anchor=(1.04, 0.5),fontsize=18,framealpha=0)
for label in (ax6.get_xticklabels() + ax6.get_yticklabels()):
    label.set_fontsize(18)
##########
pass
    
print("fig retau")
ax1.grid(True)
fig1.tight_layout()
fig1.savefig(fig+'_Retau.eps')
plt.close()
print("fig source")
ax2.grid(True)
fig2.tight_layout()
fig2.savefig(fig+'_S.eps')
plt.close()
print("fig vitesse")
ax3.grid(True)
fig3.tight_layout()
fig3.savefig(fig+'_ul_uv.eps')
plt.close()
print("debit")
ax4.grid(True)
fig4.savefig(fig+'_debits.eps')
fig4.tight_layout()
plt.close()
print("Re_b")
ax5.grid(True,"both")
fig5.tight_layout()
fig5.savefig(fig+'_Reb.eps')
plt.close()
print("relative")
ax6.grid(True)
fig6.tight_layout()
fig6.savefig(fig+'_ur.eps')
plt.close()
    
    
exit()

