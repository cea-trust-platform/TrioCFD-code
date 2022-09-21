# -*- coding: utf-8 -*-
import numpy as np
import math as mt
import matplotlib.pyplot as plt
import glob
import sys
import DNSTools3 as dtool

"""
 Trace l'evolution temporelle de le moyenne de phase de l'énergie cinétique et de la dissipation
 NB : La moyenne de phase est calculee selon l'equation (1.22) du manuscript de these d'Antoine duCluzeau
            moy_ph(q)^k = moy_ens(chi_k q) / moy_ens(chi_k)
      Dans le cas précis de ce script, la moyenne d'ensemble est effectuée sur l'ensemble des tranches de cote Z.
      On aurai pu décider que la moyenne d'ensemble s'operait sur l'ensemble des tranches Z ET des instants de releves t.
      On aurai aussi pu décider que la moyenne d'ensemble s'operait sur l'ensemble des instants de releves t uniquement,
      une fois un certain régime établi atteint.
      Ce dernier point de vue est adopte dans le script vitesses_phases_Dodd (au 19 avril 2021). Il n'est pas adopte ici pour des
      raisons evidente (évolution avec le temps attendue).
"""
############################################################## Mode d'emploi #######################################
# mode d'emploi : python energy.py liq_vap_mel nom_des_figures fichier.data chemin_vers_diphasique_moyenne_spatiale**.txt 
#                          0              1           2                 3                                   4             
#                        fichier_out_ou_log_pour_Reb...   N_smooth  [t_deb,t_fin]
#                            5                                6           7
# La on recommence from scratch en s'inspirant de ce qui a déjà été fait
############################################################## fin Mode d'emploi #######################################
message = "\
###################################################################\n\
# Script energy.py \n\
# -> Evolution temporelle de la moyenne de phase de l'énergie cinétique et de la dissipation\n\
####\n"
print(message)

## Traitement des entrees utilisateur
phase = sys.argv[1]
print('phase : ', phase)
if phase=='liq' : print('On n\'a pas IvdUdx ni dissip vapeur (chiv*dissip) \n grep -r \"AJOUT*.DISSIP\" dans IJK')

figname = sys.argv[2]
print('nom des figures : ', figname)

chemin_txt = sys.argv[4]
print('fichiers txt : ', chemin_txt)

jdd = sys.argv[3]
print('jeu de données',jdd)

accplot = sys.argv[5]
print('fichier de sortie du premier post-trt ', accplot)

print("smoothing : ",sys.argv[6])
N_smooth = max(0,int(sys.argv[6]))
filtre = np.ones(N_smooth)/float(N_smooth)
print('temps',sys.argv[7])

t_deb, t_fin = float(sys.argv[7].rstrip("]").lstrip("[").split(",")[-2]),float(sys.argv[7].rstrip("]").lstrip("[").split(",")[-1])
if (t_deb <= 0.):
    t_deb = 0
if (t_fin <= 0.):
    t_fin = 1e6

## Lecture des parametres physique de la simulaiton
mu = dtool.getParam(jdd, "mu_liquide")
rho = dtool.getParam(jdd, "rho_liquide")
rho_v = dtool.getParam(jdd, "rho_vapeur")
nu = mu/rho

## Preparation nom generique des fichiers utilises
fstat = glob.glob(chemin_txt+"/diphasique_moyenne_spatiale_*.txt")
if fstat == []:
    fstat = glob.glob(chemin_txt+"/monophasique_moyenne_spatiale_*.txt")
    print("Simulation monophasique")
else :
    print("Simulation diphasique")

fstat.sort()

################ On récupère quelques donnees du fichier de sortie du premier post-trt ##########
lines = open(accplot,"r").readlines()
Reb = float(lines[[i for i, x in enumerate(lines) if (("Reb= " in x) or ("Reb = " in x))][-1]].rstrip("\n").split("= ")[-1])
ur = float(lines[[i for i, x in enumerate(lines) if (("ur= " in x) or ("ur = " in x))][-1]].rstrip("\n").split("= ")[-1])
av = float(lines[[i for i, x in enumerate(lines) if (("alv= " in x) or ("alv = " in x))][-1]].rstrip("\n").split("= ")[-1])
g = float(lines[[i for i, x in enumerate(lines) if (("g= " in x) or ("g = " in x))][-1]].rstrip("\n").split("= ")[-1])
delta_rho = rho-rho_v
estimation_Prd = av*delta_rho*g*ur # car av est en %
print("Taux de vide : "+str(av))
print("delta_rho : "+str(delta_rho))
print("g : "+str(g))
print("ur : "+str(ur))
print("Taux de production du aux bulles (estimation très grossière) : "+str(estimation_Prd))
################# Pour récuperer les numéros de colonne qu'il nous faut #########
if phase=="liq":
    n_col_I = (2-1)
    n_col_uI = (3-1,4-1,5-1)
    n_col_uIv = (7-1,8-1,9-1)
    n_col_uuI = (11-1,12-1,13-1,14-1,15-1,16-1)
    n_col_dissip = (228-1)
    n_col_diujI = (72-1,73-1,74-1,75-1,76-1,77-1,78-1,79-1,80-1)
if phase=="vap":
    n_col_I = (2-1)
    n_col_uI = (7-1,8-1,9-1)
    n_col_uuI = (229-1,230-1,231-1,232-1,233-1,234-1)
    # n_col_uuI = (235-1,236-1,237-1,232-1,233-1,234-1)
    # n_col_dissip = (228-1) ON N'A PAS LA chi_v*dissip !!!
    # n_col_diujI = (72-1,73-1,74-1,75-1,76-1,77-1,78-1,79-1,80-1)
if phase=="mel":
    n_col_I = (2-1)
    n_col_uI = [(3-1,4-1,5-1) , (7-1,8-1,9-1)]
    n_col_uuI = [(11-1,12-1,13-1,14-1,15-1,16-1) , (229-1,230-1,231-1,232-1,233-1,234-1)]
    # n_col_uuI = [(11-1,12-1,13-1,14-1,15-1,16-1),]
    # n_col_dissip = [(228-1),]
    # n_col_diujI = [(72-1,73-1,74-1,75-1,76-1,77-1,78-1,79-1,80-1),]
#################################################################################
# temps
t = []
# energie cinetique
KT,KTB,KM,KF = [],[],[],[]
# dissipation visqueuse
ET,EM,EFB,EF = [],[],[],[]
# production / injection de qdm
Pr = []

# print((t[-1],t_deb) , (t[-1],t_fin))
if phase!="mel":
    for f in fstat:
        t.append(float("."+f.split(".")[-2]))
        sec = float(f.split(".")[-3].split("_")[-1])
        t[-1]+=sec
        #
        if ((t[-1]>=t_deb) and (t[-1]<=t_fin)): 
            if phase=="vap":
                I = 1-np.loadtxt(f,usecols=(n_col_I)).T
            elif phase=="liq":
                I = np.loadtxt(f,usecols=(n_col_I)).T
            u_l = np.loadtxt(f,usecols=(n_col_uI)).T
            uu_l = np.loadtxt(f,usecols=(n_col_uuI)).T
            if phase=='liq':
                    dissip_l = np.loadtxt(f,usecols=(n_col_dissip)).T
                    diuj_l = np.loadtxt(f,usecols=(n_col_diujI)).T
            # 
            m_I = np.mean(I,axis=-1)
            m_u_l = np.mean(u_l,axis=-1)
            m_uu_l = np.mean(uu_l,axis=-1)
            if phase=='liq':
                    m_dissip_l = np.mean(dissip_l,axis=-1)
                    m_diuj_l = np.mean(diuj_l,axis=-1)
            # 
            m_u_ll = m_u_l/m_I
            m_uu_ll = m_uu_l/m_I
            if phase=='liq':
                    # Pour construire ur, pour avoir Production
                    m_Iv = 1-m_I
                    u_v = np.loadtxt(f,usecols=(n_col_uIv)).T
                    m_u_v = np.mean(u_v,axis=-1)
                    m_u_vv = m_u_v/m_Iv
                    m_dissip_ll = m_dissip_l/m_I
                    m_diuj_ll =  m_diuj_l/m_I
                # 
                    # print(r"$\overline{s_{ij}=$",np.sum(m_diuj_ll))
            # 
            u_l_u_l = np.einsum("...i,...j->...ij",m_u_ll,m_u_ll)
            u_l_u_l = np.array((u_l_u_l[0][0],u_l_u_l[1][1],u_l_u_l[2][2],u_l_u_l[0][1],u_l_u_l[1][2],u_l_u_l[0][2]))
            m_Rij_ll = m_uu_ll - u_l_u_l
            # 
            k_m = 0.5*(u_l_u_l[0] + u_l_u_l[1] + u_l_u_l[2])
            k_t = 0.5*(m_uu_ll[0] + m_uu_ll[1] + m_uu_ll[2])
            k_f = 0.5*(m_Rij_ll[0] + m_Rij_ll[1] + m_Rij_ll[2])
            k_t_bis = k_f + k_m
            #
            if phase=='liq': 
                    e_t = -2*nu*(m_dissip_ll)
                    e_m = -2*nu*(np.sum(m_diuj_ll))
                    e_f = e_t-e_m
                    u_r = np.sqrt((m_u_vv[0]-m_u_ll[0])**2+(m_u_vv[1]-m_u_ll[1])**2+(m_u_vv[2]-m_u_ll[2])**2)
                    prd = (1./rho)*av*delta_rho*g*u_r
            # 
            KT.append(k_t);KM.append(k_m);KF.append(k_f);KTB.append(k_t_bis);
            if phase=='liq':
                    ET.append(e_t);EM.append(e_m);EF.append(e_f);#EFB.append(e_f_bis)
                    Pr.append(prd)
        else:
            t.pop()
else:
    for f in fstat:
        t.append(float("."+f.split(".")[-2]))
        sec = float(f.split(".")[-3].split("_")[-1])
        t[-1]+=sec
        #
        if ((t[-1]>=t_deb) and (t[-1]<=t_fin)): 
            Iv = 1-np.loadtxt(f,usecols=(n_col_I)).T
            I = np.loadtxt(f,usecols=(n_col_I)).T

            u_l = np.loadtxt(f,usecols=(n_col_uI[0])).T
            u_v = np.loadtxt(f,usecols=(n_col_uI[1])).T
            u_t = u_l+u_v
            uu_l = np.loadtxt(f,usecols=(n_col_uuI[0])).T
            uu_v = np.loadtxt(f,usecols=(n_col_uuI[1])).T
            uu_t = uu_l+uu_v
            # 
            m_I = np.mean(I,axis=-1)
            m_Iv = np.mean(Iv,axis=-1)
            m_uu_t = np.mean(uu_t,axis=-1)
            m_u_t = np.mean(u_t,axis=-1)
            # 
            # 
            u_t_u_t = np.einsum("...i,...j->...ij",m_u_t,m_u_t)
            u_t_u_t = np.array((u_t_u_t[0][0],u_t_u_t[1][1],u_t_u_t[2][2],u_t_u_t[0][1],u_t_u_t[1][2],u_t_u_t[0][2]))
            m_Rij_t = m_uu_t - u_t_u_t
            # 
            k_m = 0.5*(u_t_u_t[0] + u_t_u_t[1] + u_t_u_t[2])
            k_t = 0.5*(m_uu_t[0] + m_uu_t[1] + m_uu_t[2])
            k_f = 0.5*(m_Rij_t[0] + m_Rij_t[1] + m_Rij_t[2])
            k_t_bis = k_f + k_m
            #
            # 
            KT.append(k_t);KM.append(k_m);KF.append(k_f);KTB.append(k_t_bis);
        else:
            t.pop()    

if (phase=='liq') : del(diuj_l,dissip_l)
del(fstat,uu_l,u_l,I)

##################### Pour les  Plots #########################
#Objectif : gérer l'affichage des ylabel en fonction des valeurs de KF
KF_coeff = np.floor(mt.log10(max(KF)))
KF_ylabs = np.linspace(min(KF),max(KF),7)
KF_xlabs = np.floor(100*(np.linspace(min(t),max(t),7)))*10**(-2)

coeff = np.floor(mt.log10(max(KT)))
ylabs = np.linspace(min(KF),max(KT),7)
xlabs = KF_xlabs

mKF = np.mean(np.array(KF[-int(len(KF)/2):]))
mKF_coeff = np.floor(mt.log10(mKF))


if phase=='liq':
    print("maximum dissipation totale : ",max(ET))
    print("maximum dissipation fluctu : ",max(EF))
    E_coeff = np.floor(mt.log10(max([abs(max(EM)),abs(max(EF))])))
    E_ylabs = np.linspace(min([min(EF),min(EM)]),max([max(EF),max(EM)]),7)
    E_xlabs = KF_xlabs
    mEF = np.mean(np.array(EF[-int(len(EF)/2):]))
    mPr = np.mean(np.array(Pr[-int(len(Pr)/2):]))

# print("min(KF),max(KT),KF_coeff,min(t),max(t)")
# print(min(KF),max(KT),KF_coeff,min(t),max(t))
############################### Plots #########################
print("KT =\t",KT)
print("KM =\t",KM)
print("KF =\t",KF)
print("KTB= \t",KTB)

if phase=="liq":
    print("ET =\t ",ET)
    print("EM =\t ",EM)
    print("EF =\t ",EF)
    # print("EFB",EFB)

print("t =\t ",t)
fig, ax = plt.subplots(1,figsize=(9,7))
ax.set_title('Energies cinetiques - moyenne de phase '+phase+r' : $\overline{k}^k = \frac{\overline{\chi_k u \cdot u}}{\overline{\chi_k}}$',fontsize=18)

# print(n_max,len(t),len(np.convolve( KT,filtre,mode='same')))
if N_smooth > 0:
    # la longeur du plus petit vecteur dicte combien de points de trace on a 
    n_max = min(len(t),len(np.convolve( KT,filtre,mode='same')))
    # On sait que les energies et dissipation sont très bruitees : on va les adoucir
    ax.plot(t[:n_max],np.convolve( KT,filtre,mode='same')[:n_max],'k',label=r'total : ($U \cdot U$)')
    ax.plot(t[:n_max],np.convolve(KTB,filtre,mode='same')[:n_max],'+r',label=r"m+f : ($\overline{u} \cdot \overline{u} + u' \cdot u' $)")
    ax.plot(t[:n_max],np.convolve( KM,filtre,mode='same')[:n_max],'--',label='moyen : ($\overline{u} \cdot \overline{u}$)')
    ax.plot(t[:n_max],np.convolve( KF,filtre,mode='same')[:n_max],'-.',color='orange',label=r"fluct : ($ u' \cdot u'$); " + r'$\langle{K_F} \rangle = $ '+str(np.round(mKF*10**(-mKF_coeff),5))+r'$\cdot 10^{%s}$'%(mKF_coeff))
else :
    # Les energies et dissipation ne sont pas trop bruitees, pas besoin de les adoucir
    ax.plot(t[:], KT[:],'k',label=r'total : ($U \cdot U$)')
    ax.plot(t[:],KTB[:],'+r',label=r"m+f : ($\overline{u} \cdot \overline{u} + u' \cdot u' $)")
    ax.plot(t[:], KM[:],'--',label='moyen : ($\overline{u} \cdot \overline{u}$)')
    ax.plot(t[:], KF[:],'-.',color='orange',label=r"fluct : ($ u' \cdot u'$); " + r'$\langle{K_F} \rangle = $ '+str(np.round(mKF*10**(-mKF_coeff),5))+r'$\cdot 10^{%s}$'%(mKF_coeff))

# ax.text(0.99*(t[-1]+t[0])/2, 1.8*mKF, r'$\langle{K_F} \rangle = $ '+str(np.round(mKF*10**(-mKF_coeff),5))+r'$\cdot 10^{%s}$'%(mKF_coeff), fontsize=18)
ax.set_xlabel('temps (s)',fontsize=18)
ax.set_xticks(xlabs)
ax.set_xticklabels(np.round(xlabs,2),fontsize=18)
ax.set_ylabel(r'$\overline{k}^k \times 10^{%s}$'%(str(coeff)),fontsize=18)
ax.set_yticks(ylabs)
ax.set_yticklabels(np.round(ylabs*10**(-coeff),3),fontsize=18)

ax.tick_params(labelsize=18)
ax.legend(fontsize=18)
plt.tight_layout()
fig.savefig(figname+"_K.png")
fig.savefig(figname+"_K.pdf")

if (np.mean(np.array(KM[-int(len(KM)/2):]))/np.mean(np.array(KF[-int(len(KF)/2):])) > 6.):
    fig, ax = plt.subplots(1,figsize=(9,7))
    ax.set_title(r'Energie cinetique des fluctuations : $\overline{k}^k = \frac{\overline{\chi_k u \cdot u}}{\overline{\chi_k}}$',fontsize=18)
    


    # On adoucit le signal
    if N_smooth > 0 :
        # la longeur du plus petit vecteur dicte combien de points de trace on a 
        n_max = min(len(t),len(np.convolve( KT,filtre,mode='same')))
        ax.plot(t[:n_max],np.convolve( KF, filtre, mode='same')[:,n_max],'-.',color='orange',label='fluct')
    # On n'adoucit pas le signal
    else :
        ax.plot(t[:], KF[:],'-.',color='orange',label='fluct')
    
    ax.set_xlabel('temps (s)',fontsize=18)
    ax.set_ylabel(r'$\overline{k}^k \times 10^{%s}$'%(str(KF_coeff)),fontsize=18)
    ax.set_xticks(KF_xlabs)
    ax.set_xticklabels(np.round(KF_xlabs,2),fontsize=18)
    ax.set_yticks(KF_ylabs)
    ax.set_yticklabels(np.round(KF_ylabs*10**(-KF_coeff),3),fontsize=18)
    ax.legend(fontsize=18)
    plt.tight_layout()
    fig.savefig(figname+"_KF.png")
    fig.savefig(figname+"_KF.pdf")

if phase=='liq':
    fig, ax = plt.subplots(1,figsize=(9,7))
    ax.set_title(r'Dissipation turbulente - moyenne de phase '+phase+r' : $\overline{\epsilon}^k = \frac{\overline{\chi_k 2 \nu s_{ij} s_{ij}}}{\overline{\chi_k}}$ ',fontsize=18)
    

    if N_smooth > 0:
        # la longeur du plus petit vecteur dicte combien de points de trace on a 
        n_max = min(len(t),len(np.convolve( ET,filtre,mode='same')))    
        # Adoucissons le signal
        ax.plot(t[:n_max],np.convolve(ET, filtre, mode='same')[:n_max],'k',label="totale")
        ax.plot(t[:n_max],np.convolve(EM, filtre, mode='same')[:n_max],'--',label='moyen : ($\overline{s_{ij}} \cdot \overline{s_{ij}}$)')
        ax.plot(t[:n_max],np.convolve(EF, filtre, mode='same')[:n_max],'-.',label='fluctuant : ($\overline{s\' _{ij} \cdot s\' _{ij} }$)')
        ax.plot(t[:n_max],np.convolve(Pr, filtre, mode='same')[:n_max],'-.',label=r'$\mathcal{P} : \alpha \Delta \rho g u_r \sim $'+str(np.round(mPr*10**(-E_coeff),5))+r'$\cdot 10^{%s}$'%(E_coeff))
        # ax.text(0.99*(t[-1]+t[0])/2, 0.3*mEF, r'$\mathcal{P} \sim \alpha \Delta \rho g u_r =$ '+str(np.round(mPr*10**(-E_coeff),5))+r'$\cdot 10^{%s}$'%(E_coeff), fontsize=18)
    else:
        # N'adoucissons pas le signal
        ax.plot(t[:],ET[:],'k',label="totale")
        ax.plot(t[:],EM[:],'--',label='moyen : ($\overline{s_{ij}} \cdot \overline{s_{ij}}$)')
        ax.plot(t[:],EF[:],'-.',label='fluctuant : ($\overline{s\' _{ij} \cdot s\' _{ij} }$)')
        ax.plot(t[:],Pr[:],'-.',label=r'$\mathcal{P} : \alpha \Delta \rho g u_r \sim $'+str(np.round(mPr*10**(-E_coeff),5))+r'$\cdot 10^{%s}$'%(E_coeff))
        # ax.text(0.99*(t[-1]+t[0])/2, 0.3*mEF, r'$\mathcal{P} \sim \alpha \Delta \rho g u_r =$ '+str(np.round(mPr*10**(-E_coeff),5))+r'$\cdot 10^{%s}$'%(E_coeff), fontsize=18)    
    
    ax.set_xlabel("temps (s)",fontsize=18)
    ax.set_ylabel(r'$\overline{\epsilon}^k \times 10^{%s}$'%(str(E_coeff)),fontsize=18)
    ax.set_xticks(xlabs)
    ax.set_xticklabels(np.round(xlabs,2),fontsize=18)
    ax.set_yticks(E_ylabs)
    ax.set_yticklabels(np.round(E_ylabs*10**(-E_coeff),3),fontsize=18)
    ax.legend(fontsize=18)
    plt.tight_layout()
    fig.savefig(figname+"_Dissip.png")
    fig.savefig(figname+"_Dissip.pdf")
    
# Grandeurs affichees sur les plots 
print("####")
print("Moyenne de l'énergie cinétique des fluctuations :")
print(r'$\langle{K_F} \rangle = $ '+str(np.round(mKF*10**(-mKF_coeff),5))+r'$\cdot 10^{%s}$'%(mKF_coeff))    
print(" ==> <K_F> = ",(np.round(mKF*10**(-mKF_coeff),5))*10**(mKF_coeff))
print("Estimation de la puissance injectee par les bulles :")
print(r'$\mathcal{P} : \alpha \Delta \rho g u_r \sim $'+str(np.round(mPr*10**(-E_coeff),5))+r'$\cdot 10^{%s}$'%(E_coeff))
print(" ==> P_inj = ",(np.round(mPr*10**(-E_coeff),5))*10**(E_coeff))
message = "\
####\n\
# Fin du script energy.py \n\
######################################################################\n"
print(message)
