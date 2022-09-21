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

################################################################################################
# How to correctly sort a string with a number inside? [duplicate] ;
# https://stackoverflow.com/questions/5967500/how-to-correctly-sort-a-string-with-a-number-inside
import re

def atof(text):
    try:
        retval = float(text)
    except ValueError:
        retval = text
    return retval

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    float regex comes from https://stackoverflow.com/a/12643073/190597
    '''
    return [ atof(c) for c in re.split(r'[+-]?([0-9]+(?:[.][0-9]*)?|[.][0-9]+)', text) ]
################################################################################################

################################################################################################
def write_me_in_raw(chemin_sauvegarde,array_to_write,array_name):
            the_shapes_of_array_to_write = np.array(array_to_write.shape)
            np.save(chemin_sauvegarde+"the_shapes_of_"+array_name,the_shapes_of_array_to_write)
            array_to_write.tofile(chemin_sauvegarde+array_name+".raw")
            del(array_to_write)
################################################################################################
            
############################################################## Mode d'emploi #######################################
# mode d'emploi : python energy.py liq_vap_mel nom_des_figures fichier.data chemin_vers_diphasique_moyenne_spatiale**.txt 
#                          0              1           2                 3                                   4             
#                        fichier_out_ou_log_pour_Reb...   N_smooth  mode   [t_deb,t_fin]
#                            5                                6      7         8
# /!\ 2 : Si on met adim ou ADIM dans le nom, alors on trace les evolution de K/<K>, E/<E> en fonction de t/<T>
#         avec <T> = <K>/<E>
# 7 : draw -> dessine les graphiques directement
#     write -> ecrits les fichiers raw
#     both -> desisne les graphiques et ecrit les fichiers raw
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

adim = ("adim".upper() in figname.upper()) & (phase=='liq')
if adim :
    print("################ ADIM ####################")
    print('Traces en adimentionnel : K/<K>; E/<E>; t/(<K>/<E>)')
else :
    print("############## NOT ADIM #################")
    
chemin_txt = sys.argv[4]
print('fichiers txt : ', chemin_txt)
chemin_sauvegarde = chemin_txt.replace("TXT/","")

jdd = sys.argv[3]
print('jeu de données',jdd)

accplot = sys.argv[5]
print('fichier de sortie du premier post-trt ', accplot)

print("smoothing : ",sys.argv[6])
N_smooth = max(0,int(sys.argv[6]))
filtre = np.ones(N_smooth)/float(N_smooth)


# mode = draw / raw /both
mode = sys.argv[7]

print('temps',sys.argv[8])
t_deb, t_fin = float(sys.argv[8].rstrip("]").lstrip("[").split(",")[-2]),float(sys.argv[8].rstrip("]").lstrip("[").split(",")[-1])
if (t_deb <= 0.):
    t_deb = 0
if (t_fin <= 0.):
    t_fin = 1e6
print("-> t_deb : %f, t_fin : %f",t_deb,t_fin)

## Lecture des parametres physique de la simulaiton
mu = dtool.getParam(jdd, "mu_liquide")
rho = dtool.getParam(jdd, "rho_liquide")
diphasique = (not (dtool.getFlag(jdd,"disable_diphasique")))
if (diphasique):
    rho_v = dtool.getParam(jdd, "rho_vapeur")
nu = mu/rho

## Preparation nom generique des fichiers utilises
fstat = glob.glob(chemin_txt+"/diphasique_moyenne_spatiale_*.txt")
if fstat == []:
    fstat = glob.glob(chemin_txt+"/monophasique_moyenne_spatiale_*.txt")
    diphasique=False
    print("Simulation monophasique")
else :
    diphasique=True
    print("Simulation diphasique")
print("## DIPHASIQUE : ", diphasique)
# Tri des fichiers fstat par ordre croissant des seccondes. 
# fstat.sort() --> Mauvais car place la seconde 10 juste apres la fin de la seconde 2, et pas 9
fstat.sort(key=natural_keys)
################ On récupère quelques donnees du fichier de sortie du premier post-trt ##########
if diphasique:
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
else:
    print("simulation monophasique")
################# Pour récuperer les numéros de colonne qu'il nous faut #########
if phase=="liq":
    n_col_I = (2-1)
    n_col_uI = (3-1,4-1,5-1)
    n_col_uIv = (7-1,8-1,9-1)
    # """interpolation 1 : uu = [ 0.5(u_i + u_{i+1/2}) ]^2"""
    n_col_uuI = (11-1,12-1,13-1,14-1,15-1,16-1)
    # """interpolation 2 : uu = 0.5(u_i^2 + u_{i+1/2}^2)"""
    # n_col_uuI = (238-1,239-1,240-1,14-1,15-1,16-1)
    n_col_dissip = (228-1)
    n_col_dissip_vap = (796-1)
    try:
        n_col_true_dissip = (797-1)
        n_col_true_dissip_vap = (798-1)
    except:
        print("Old version of IJK : True_dissip not post-processed directly.")
    n_col_diujI = (72-1,73-1,74-1,75-1,76-1,77-1,78-1,79-1,80-1)  # ordre : dudx,dvdx,dwdx,dudy,dvdy,dwdy,dudz,dvdz,dwdz
if phase=="vap":
    n_col_I = (2-1)
    n_col_uI = (7-1,8-1,9-1)
    n_col_uuI = (229-1,230-1,231-1,232-1,233-1,234-1)
    # n_col_uuI = (235-1,236-1,237-1,232-1,233-1,234-1)
    # n_col_dissip_vap = (796-1) 
    # n_col_true_dissip_vap = (798-1) 
    # n_col_diujI = (72-1,73-1,74-1,75-1,76-1,77-1,78-1,79-1,80-1)
if phase=="mel":
    n_col_I = (2-1)
    n_col_uI = [(3-1,4-1,5-1) , (7-1,8-1,9-1)]
    n_col_uuI = [(11-1,12-1,13-1,14-1,15-1,16-1) , (229-1,230-1,231-1,232-1,233-1,234-1)]
    # n_col_uuI = [(11-1,12-1,13-1,14-1,15-1,16-1),]
    # n_col_dissip = [(228-1),(796-1) ]
    # n_col_true_dissip = [(797-1),(798-1) ]
    # n_col_diujI = [(72-1,73-1,74-1,75-1,76-1,77-1,78-1,79-1,80-1),]
#################################################################################
# temps
t = []
# energie cinetique
KT,KTB,KM,KF = [],[],[],[]
# dissipation visqueuse
ET,EM,EFB,EF = [],[],[],[]
T_ET,T_EM,T_EFB,T_EF = [],[],[],[]
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
                diuj_l = np.loadtxt(f,usecols=(n_col_diujI)).T    # ordre     :  dudx,dvdx,dwdx,dudy,dvdy,dwdy,dudz,dvdz,dwdz
                sij_l = diuj_l.reshape((3,3,diuj_l.shape[-1]))    # [i,j,...] -> diuj[...]
                sij_l = 0.5*(sij_l+sij_l.transpose(1,0,2))        # sij = djui + diuj
                try:
                    true_dissip_l = np.loadtxt(f,usecols=(n_col_true_dissip)).T
                except:
                    true_dissip_l = -2*np.einsum("ij...,ij...->...",sij_l,sij_l)
            # 
            m_I = np.mean(I,axis=-1)
            m_u_l = np.mean(u_l,axis=-1)
            m_uu_l = np.mean(uu_l,axis=-1)
            if phase=='liq':
                m_dissip_l = np.mean(dissip_l,axis=-1)
                m_true_dissip_l = np.mean(true_dissip_l,axis=-1)
                m_diuj_l = np.mean(diuj_l,axis=-1)
                m_sij_l = np.mean(sij_l,axis=-1)
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
                    m_true_dissip_ll = m_true_dissip_l/m_I
                    m_diuj_ll =  m_diuj_l/m_I
                    m_sij_ll = m_sij_l/m_I
                # 
                    # print(r"$\overline{s_{ij}=$",np.sum(m_diuj_ll))
            # 
            ############################################################
            # Il y a un grain dans la puree ici 
            ############################################################
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
                    # pseudo-e = nu*( djui djui )
                    # true-e   = 2*nu*( sij sij) = pseudo-e + nu*( diuj djui)
                    e_t = -1*nu*(m_dissip_ll)
                    true_e_t = -1*nu*(m_true_dissip_ll) 
                    e_m = -1*nu*(np.einsum("ji...,ji->...",m_diuj_ll.reshape(3,3),m_diuj_ll.reshape(3,3)))
                    true_e_m = -1*nu*(np.einsum("ij...,ij...->...",m_sij_ll,m_sij_ll))  # double-produit contracte de sij moyen.
                    e_f = e_t-e_m
                    true_e_f = true_e_t-true_e_m
                    u_r = np.sqrt((m_u_vv[0]-m_u_ll[0])**2+(m_u_vv[1]-m_u_ll[1])**2+(m_u_vv[2]-m_u_ll[2])**2)
                    if diphasique:
                        prd = (1./rho)*av*delta_rho*g*u_r
            # 
            KT.append(k_t);KM.append(k_m);KF.append(k_f);KTB.append(k_t_bis);
            if phase=='liq':
                    ET.append(e_t);EM.append(e_m);EF.append(e_f);#EFB.append(e_f_bis)
                    T_ET.append(true_e_t);T_EM.append(true_e_m);T_EF.append(true_e_f);#EFB.append(e_f_bis)
                    if diphasique:
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
print("m_I : ",m_I)
##################### Pour les  Plots #########################
#Objectif : gérer l'affichage des ylabel en fonction des valeurs de KF
KF_coeff = np.floor(mt.log10(max(KF)))
KF_ylabs = np.linspace(min(KF),max(KF),7)
KF_xlabs = np.floor(100*(np.linspace(min(t),max(t),7)))*10**(-2)

coeff = np.floor(mt.log10(max(KT)))
ylabs = np.linspace(min(KF),max(KT),7)
xlabs = KF_xlabs

mKF = np.abs(np.mean(np.array(KF[-int(len(KF)/2):]))) # moyenne de l'Ec de fluctuations : <u'.u'>
mKF_coeff = np.floor(mt.log10(mKF))


if phase=='liq':
    print("maximum dissipation totale (pseudo et true) : ",max(ET),max(T_ET))
    print("maximum dissipation fluctu (pseudo et true) : ",max(EF),max(T_EF))
    E_coeff = np.floor(mt.log10(max([abs(max(EM)),abs(max(EF))])))
    T_E_coeff = np.floor(mt.log10(max([abs(max(T_EM)),abs(max(T_EF))])))
    E_ylabs = np.linspace(min([min(EF),min(EM)]),max([max(EF),max(EM)]),7)
    T_E_ylabs = np.linspace(min([min(T_EF),min(T_EM)]),max([max(T_EF),max(T_EM)]),7)
    E_xlabs = KF_xlabs
    T_E_xlabs = KF_xlabs
    mEF = np.mean(np.array(EF[-int(len(EF)/2):])) # moyenne de la dissipations due aux fluctuations : < \ol{s'_ij. s'_ij} >
    mT_EF = np.mean(np.array(T_EF[-int(len(T_EF)/2):])) # moyenne de la true_dissipations due aux fluctuations : < \ol{s'_ij. s'_ij} >
    if diphasique:
        mPr = np.mean(np.array(Pr[-int(len(Pr)/2):])) # estimation de la dissipation due a la flotabilite : \a \D\r g u_r
    
    print("mEF,mEF")
    print(mEF,mT_EF)
# print("min(KF),max(KT),KF_coeff,min(t),max(t)")
# print(min(KF),max(KT),KF_coeff,min(t),max(t))
########################################################

### REPENSER L'ECRITURE POUR QUE CE SOIT FACILEMENT TRACABLE AVEC GNUPLOT,
### OU AVEC PYHTON EN LECTURE DIRECTE.
print("KT =\t",KT)
print("KM =\t",KM)
print("KF =\t",KF)
print("KTB= \t",KTB)

if phase=="liq":
    print("### pseudo dissipations")
    print("ET =\t ",ET)
    print("EM =\t ",EM)
    print("EF =\t ",EF)
    print("### true dissipations")
    print("T_ET =\t ",T_ET)
    print("T_EM =\t ",T_EM)
    print("T_EF =\t ",T_EF)
    # print("EFB",EFB)

print("t =\t ",t)

t_t = t
if diphasique:
    T_Pr = Pr
if adim :
    T_L = mKF/np.abs(mEF)
    T_T_L = mKF/np.abs(mT_EF)
    t = t/T_L
    t_t = t_t/T_T_L
    KT,KM,KF,KTB = KT/mKF,KM/mKF,KF/mKF,KTB/mKF
    ET,EM,EF = ET/mEF,EM/mEF,EF/mEF
    T_ET,T_EM,T_EF = T_ET/mT_EF,T_EM/mT_EF,T_EF/mT_EF
    if diphasique:
        Pr = Pr/mEF
        T_Pr = T_Pr/mT_EF

#################################################################################
# TRACE DES FIGURES OU ECRITURE DES DONNEES - ENERGIES CINETIQUES
#################################################################################
draw_the_graphs = (("draw" in mode) or ("both" in mode))
write_the_raws =  (("raw" in mode) or ("both" in mode) or ("write" in mode))
print("mode : ",mode)
#################################################################################

if draw_the_graphs:
    fig, ax = plt.subplots(1,figsize=(9,7))
    ax.set_title('Energies cinetiques - moyenne de phase '+phase+r' : $\overline{k}^k = \frac{\overline{\chi_k u \cdot u}}{\overline{\chi_k}}$',fontsize=18)
    
    if N_smooth > 0:
        # la longeur du plus petit vecteur dicte combien de points de trace on a 
        n_max = min(len(t),len(np.convolve( KT,filtre,mode='same'))-N_smooth)
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
    
    if adim:
        ax.set_xlabel(r't / $ T_L ; T_L = \langle{K_F} \rangle /\langle{\epsilon _f} \rangle = $'+str(np.round(T_L,2)),fontsize=18)
        ax.set_ylabel(r'$\overline{k}^k / \langle{K_F} \rangle$',fontsize=18)
    else:
        ax.set_xlabel('temps (s)',fontsize=18)
        ax.set_xticks(xlabs)
        ax.set_xticklabels(np.round(xlabs,2),fontsize=18)
        ax.set_ylabel(r'$\overline{k}^k \times 10^{%s}$'%(str(coeff)),fontsize=18)
        ax.set_yticks(ylabs)
        ax.set_yticklabels(np.round(ylabs*10**(-coeff),3),fontsize=18)
    
    ax.tick_params(labelsize=22)
    for label in (ax.get_xticklabels() + ax.get_yticklabels()):
            label.set_fontsize(22)
    ax.legend(fontsize=22)
    plt.tight_layout()
    fig.savefig(chemin_sauvegarde+figname+"_K.png")
    fig.savefig(chemin_sauvegarde+figname+"_K.pdf")
    
    # ENERGIE CINETIQUE DES FLUCTUATIONS (si KM et KF sont tres differents uniquement)
    if (np.mean(np.array(KM[-int(len(KM)/2):]))/np.mean(np.array(KF[-int(len(KF)/2):])) > 6.):
        fig, ax = plt.subplots(1,figsize=(9,7))
        ax.set_title(r'Energie cinetique des fluctuations : $\overline{k}^k = \frac{\overline{\chi_k u \cdot u}}{\overline{\chi_k}}$',fontsize=18)
        
        # On adoucit le signal
        if N_smooth > 0 :
            # la longeur du plus petit vecteur dicte combien de points de trace on a 
            n_max = min(len(t),len(np.convolve( KT,filtre,mode='same'))-N_smooth)
            print("KM et KF tres differents")
            if adim :
                print(t.shape, KF.shape)
            print(np.convolve( KF, filtre, mode='same'))
            ax.plot(t[:n_max],np.convolve( KF, filtre, mode='same')[:n_max],'-.',color='orange',label='fluct')
        # On n'adoucit pas le signal
        else :
            ax.plot(t[:], KF[:],'-.',color='orange',label='fluct')
    
    
        if adim:
            ax.set_xlabel(r't / $ T_L ; T_L = \langle{K_F} \rangle /\langle{\epsilon _f} \rangle = $'+str(np.round(T_L,2)),fontsize=18)
            ax.set_ylabel(r'$\overline{k}^k / \langle{K_F} \rangle$',fontsize=18)
        else:    
            ax.set_xlabel('temps (s)',fontsize=18)
            ax.set_ylabel(r'$\overline{k}^k \times 10^{%s}$'%(str(KF_coeff)),fontsize=18)
            ax.set_xticks(KF_xlabs)
            ax.set_xticklabels(np.round(KF_xlabs,2),fontsize=22)
            ax.set_yticks(KF_ylabs)
            ax.set_yticklabels(np.round(KF_ylabs*10**(-KF_coeff),3),fontsize=22)
        for label in (ax.get_xticklabels() + ax.get_yticklabels()):
            label.set_fontsize(22)
        ax.legend(fontsize=22)
        plt.tight_layout()
        fig.savefig(chemin_sauvegarde+figname+"_KF.png")
        fig.savefig(chemin_sauvegarde+figname+"_KF.pdf")
        
        
        
    
    # PSEUDO DISSIPATION TURBULENTE \tilde{\epsilon}
    if phase=='liq':
        fig, ax = plt.subplots(1,figsize=(9,7))
        ax.set_title(r'Pseudo - Dissipation turbulente - moyenne de phase '+phase+r' : $\overline{\epsilon}^k = \frac{\overline{\chi_k 2 \nu \partial_j u_i \partial_j u_i }}{\overline{\chi_k}}$ ',fontsize=18)
        
    
        if N_smooth > 0:
            # la longeur du plus petit vecteur dicte combien de points de trace on a 
            n_max = min(len(t),len(np.convolve( ET,filtre,mode='same'))-N_smooth)    
            # Adoucissons le signal
            ax.plot(t[:n_max],np.convolve(ET, filtre, mode='same')[:n_max],'k',label="totale")
            ax.plot(t[:n_max],np.convolve(EM, filtre, mode='same')[:n_max],'--',label='moyen : ($\overline{s_{ij}} \cdot \overline{s_{ij}}$)')
            ax.plot(t[:n_max],np.convolve(EF, filtre, mode='same')[:n_max],'-.',label='fluct. : ($\overline{s\' _{ij} \cdot s\' _{ij} }$) $\sim$' +str(np.round(mEF*10**(-E_coeff),5))+r'$\cdot 10^{%s}$'%(E_coeff))
            if diphasique:
                ax.plot(t[:n_max],np.convolve(Pr, filtre, mode='same')[:n_max],'-.',label=r'$\mathcal{P} : \alpha \Delta \rho g u_r \sim $'+str(np.round(mPr*10**(-E_coeff),5))+r'$\cdot 10^{%s}$'%(E_coeff))
            # ax.text(0.99*(t[-1]+t[0])/2, 0.3*mEF, r'$\mathcal{P} \sim \alpha \Delta \rho g u_r =$ '+str(np.round(mPr*10**(-E_coeff),5))+r'$\cdot 10^{%s}$'%(E_coeff), fontsize=18)
        else:
            # N'adoucissons pas le signal
            ax.plot(t[:],ET[:],'k',label="totale")
            ax.plot(t[:],EM[:],'--',label='moyen : ($\overline{s_{ij}} \cdot \overline{s_{ij}}$)')
            ax.plot(t[:],EF[:],'-.',label='fluct. : ($\overline{s\' _{ij} \cdot s\' _{ij} }$) $\sim$' +str(np.round(mEF*10**(-E_coeff),5))+r'$\cdot 10^{%s}$'%(E_coeff))
            if diphasique:
                ax.plot(t[:],Pr[:],'-.',label=r'$\mathcal{P} : \alpha \Delta \rho g u_r \sim $'+str(np.round(mPr*10**(-E_coeff),5))+r'$\cdot 10^{%s}$'%(E_coeff))
            # ax.text(0.99*(t[-1]+t[0])/2, 0.3*mEF, r'$\mathcal{P} \sim \alpha \Delta \rho g u_r =$ '+str(np.round(mPr*10**(-E_coeff),5))+r'$\cdot 10^{%s}$'%(E_coeff), fontsize=18)    
    
        if adim:
            ax.set_xlabel(r't / $ T_L ; T_L = \langle{K_F} \rangle /\langle{\epsilon _f} \rangle = $'+str(np.round(T_L,2)),fontsize=18)
            ax.set_ylabel(r'$\overline{\epsilon}^k / \langle{\epsilon} \rangle, \mathcal{P} / \langle{\epsilon} \rangle$',fontsize=18)
        else:
            ax.set_xlabel("temps (s)",fontsize=18)
            ax.set_ylabel(r'$\overline{\epsilon}^k \times 10^{%s}$'%(str(E_coeff)),fontsize=18)
            ax.set_xticks(xlabs)
            ax.set_xticklabels(np.round(xlabs,2),fontsize=22)
            ax.set_yticks(E_ylabs)
            ax.set_yticklabels(np.round(E_ylabs*10**(-E_coeff),3),fontsize=22)
        for label in (ax.get_xticklabels() + ax.get_yticklabels()):
            label.set_fontsize(22)
        ax.legend(fontsize=22)
        plt.tight_layout()
        fig.savefig(chemin_sauvegarde+figname+"_Dissip.png")
        fig.savefig(chemin_sauvegarde+figname+"_Dissip.pdf")
    
    if adim:
        fig, ax = plt.subplots(1,figsize=(9,7))
        ax.set_title(r'Moyenne de phase '+phase,fontsize=22)
    
        ax.plot(t[:],EF[:],'-',color='red',label=' $\epsilon = (\overline{s\' _{ij} \cdot s\' _{ij} }$);'+r'$\langle \epsilon \rangle \sim$' +str(np.round(mEF*10**(-E_coeff),5))+r'$\cdot 10^{%s}$'%(E_coeff))
        ax.plot(t[:], KF[:],'-',color='black',label=r"$K_F =$  ($ u' \cdot u'$); " + r'$\langle {K_F} \rangle = $ '+str(np.round(mKF*10**(-mKF_coeff),5))+r'$\cdot 10^{%s}$'%(mKF_coeff))
    
        ax.set_xlabel(r't / $ T_L ; T_L = \frac{\langle{K_F} \rangle }{\langle{\epsilon _f} \rangle} = $'+str(np.round(T_L,2)),fontsize=18)
        ax.set_ylabel(
                r'$\overline{k}^k / \langle{K_F} \rangle$'
                +', '+
                r'$\overline{\epsilon}^k / \langle{\epsilon} \rangle $'
                ,fontsize=22)
        for label in (ax.get_xticklabels() + ax.get_yticklabels()):
            label.set_fontsize(22)
        ax.legend(fontsize=22)
        plt.tight_layout()
        fig.savefig(chemin_sauvegarde+figname+"_K_D.png")
        fig.savefig(chemin_sauvegarde+figname+"_K_D.pdf")
        
        
    # VRAIE DISSIPATION TURBULENTE \epsilon 
    if phase=='liq':   
        fig, ax = plt.subplots(1,figsize=(9,7))
        ax.set_title(r'Vraie dissipation turbulente - moyenne '+phase+
                    r' : $\overline{{\epsilon}}^k = \frac{\overline{\chi_k 2 \nu s_{ij} s_{ij}}}{\overline{\chi_k}} $'+
                    r'$= \frac{\overline{\chi_k \tilde{\epsilon} + \nu \partial_j u_i \partial_i u_j}}{\overline{\chi_k}} $',fontsize=18)
            
        if N_smooth > 0:
            # la longeur du plus petit vecteur dicte combien de points de trace on a 
            n_max = min(len(t),len(np.convolve( T_ET,filtre,mode='same'))-N_smooth)    
            # Adoucissons le signal
            ax.plot(t_t[:n_max],np.convolve(T_ET, filtre, mode='same')[:n_max],'k',label="totale")
            ax.plot(t_t[:n_max],np.convolve(T_EM, filtre, mode='same')[:n_max],'--',label=r'moyen : ($\nu \overline{s_{ij}} \cdot \overline{s_{ij}}$)')
            ax.plot(t_t[:n_max],np.convolve(T_EF, filtre, mode='same')[:n_max],'-.',label=r"fluct. : ($\nu \overline{s' _{ij} \cdot s' _{ij} }$) $\sim$" +str(np.round(mT_EF*10**(-T_E_coeff),5))+r'$\cdot 10^{%s}$'%(T_E_coeff))
            if diphasique:
                ax.plot(t_t[:n_max],np.convolve(T_Pr, filtre, mode='same')[:n_max],'-.',label=r'$\mathcal{P} : \alpha \Delta \rho g u_r \sim $'+str(np.round(mPr*10**(-T_E_coeff),5))+r'$\cdot 10^{%s}$'%(T_E_coeff))
        else:
            # N'adoucissons pas le signal
            ax.plot(t_t[:],T_ET[:],'k',label="totale")
            ax.plot(t_t[:],T_EM[:],'--',label=r'moyen : ($\nu \overline{s_{ij}} \cdot \overline{s_{ij}}$)')
            ax.plot(t_t[:],T_EF[:],'-.',label=r"fluct. : ($\nu \overline{s' _{ij} \cdot s' _{ij} }$) $\sim$" +str(np.round(mT_EF*10**(-T_E_coeff),5))+r'$\cdot 10^{%s}$'%(T_E_coeff))
            if diphasique:
                ax.plot(t_t[:],T_Pr[:],'-.',label=r'$\mathcal{P} : \alpha \Delta \rho g u_r \sim $'+str(np.round(mPr*10**(-T_E_coeff),5))+r'$\cdot 10^{%s}$'%(T_E_coeff))
        
        if adim:
            ax.set_xlabel(r't / $ T_L ; T_L = \langle{K_F} \rangle /\langle{\tilde{\epsilon _f}} \rangle = $'+str(np.round(T_T_L,2)),fontsize=18)
            ax.set_ylabel(r'$\overline{\tilde{\epsilon}}^k / \langle{{\epsilon}} \rangle $',fontsize=18)
        else:
            ax.set_xlabel("temps (s)",fontsize=18)
            ax.set_ylabel(r'$\overline{\tilde{\epsilon}}^k \times 10^{%s}$'%(str(T_E_coeff)),fontsize=18)
            ax.set_xticks(xlabs)
            ax.set_xticklabels(np.round(xlabs,2),fontsize=22)
            ax.set_yticks(T_E_ylabs)
            ax.set_yticklabels(np.round(T_E_ylabs*10**(-T_E_coeff),3),fontsize=22)
        for label in (ax.get_xticklabels() + ax.get_yticklabels()):
            label.set_fontsize(22)
        ax.legend(fontsize=22)
        plt.tight_layout()
        fig.savefig(chemin_sauvegarde+figname+"_TrDissip.png")
        fig.savefig(chemin_sauvegarde+figname+"_TrDissip.pdf")
        
        
    if adim:
        fig, ax = plt.subplots(1,figsize=(9,7))
        ax.set_title(r'Moyenne '+phase,fontsize=18)
    
        ax.plot(t[:],T_EF[:],'-',color='red',label=r" $\epsilon = (\overline{\nu s' _{ij} \cdot s' _{ij} }$);"+r'$\langle \epsilon \rangle \sim$' +str(np.round(mT_EF*10**(-T_E_coeff),5))+r'$\cdot 10^{%s}$'%(T_E_coeff))
        ax.plot(t[:], KF[:],'-',color='black',label=r"$K_F =$  ($ u' \cdot u'$); " + r'$\langle {K_F} \rangle = $ '+str(np.round(mKF*10**(-mKF_coeff),5))+r'$\cdot 10^{%s}$'%(mKF_coeff))
    
        ax.set_xlabel(r't / $ T_L ; T_L = \frac{\langle{K_F} \rangle }{\langle{\epsilon _f} \rangle} = $'+str(np.round(T_L,2)),fontsize=18)
        ax.set_ylabel(
                r'$\overline{k}^k / \langle{K_F} \rangle$'
                +', '+
                r'$\overline{\epsilon}^k / \langle{\epsilon} \rangle $'
                ,fontsize=18)
        for label in (ax.get_xticklabels() + ax.get_yticklabels()):
            label.set_fontsize(22)
        ax.legend(fontsize=22)
        plt.tight_layout()
        fig.savefig(chemin_sauvegarde+figname+"_K_TrD.png")
        fig.savefig(chemin_sauvegarde+figname+"_K_TrD.pdf")

if write_the_raws:
    
    # Ecriture des information utiles pour les legendes
    if adim:
        TL__TTL = np.array((T_L, T_T_L))
        np.save(chemin_sauvegarde+"TL__TTL",TL__TTL)        
    E_coeff__T_E_coeff__mEF__mT_EF = np.array((E_coeff, T_E_coeff, mEF, mT_EF))
    np.save(chemin_sauvegarde+"Eco__TEco__mEF__mTEF",E_coeff__T_E_coeff__mEF__mT_EF)
    # E_coeff__mEF = np.array((E_coeff, mEF))
    # np.save(chemin_sauvegarde+"Eco__mEF",E_coeff__mEF)
    KF_coeff__coeff__mKF__mKF_coeff = np.array((KF_coeff, coeff, mKF, mKF_coeff)) 
    np.save(chemin_sauvegarde+"KFco__co__mKF__mKFco",KF_coeff__coeff__mKF__mKF_coeff)
    # if diphasique:
        # mPr = np.m
    
    # Ecriture des champs interessants
    if N_smooth > 0:
        # la longeur du plus petit vecteur dicte combien de points de trace on a 
        n_max = min(len(t),len(np.convolve( KT,filtre,mode='same'))-N_smooth)
        # On sait que les energies et dissipation sont très bruitees : on va les adoucir
        temps_KT_KTB_KM_KF = np.array((t[:n_max],
                                      np.convolve( KT,filtre,mode='same')[:n_max],
                                      np.convolve(KTB,filtre,mode='same')[:n_max],
                                      np.convolve( KM,filtre,mode='same')[:n_max],
                                      np.convolve( KF,filtre,mode='same')[:n_max]))
                                      
    else :
        # Les energies et dissipation ne sont pas trop bruitees, pas besoin de les adoucir
        temps_KT_KTB_KM_KF = np.array((t[:],KT[:],KTB[:],KM[:],KF[:]))
    
    # ecriture en raw. efface temps_KT_KTB_KM_KF de la memoire ensuite
    write_me_in_raw(chemin_sauvegarde,temps_KT_KTB_KM_KF,"temps_KT_KTB_KM_KF")
    
    # ENERGIE CINETIQUE DES FLUCTUATIONS (si KM et KF sont tres differents uniquement)
    if (np.mean(np.array(KM[-int(len(KM)/2):]))/np.mean(np.array(KF[-int(len(KF)/2):])) > 6.):
        print("KF et KM sont tres differents, ca vaut le coup de les plotter separement")
    
    # PSEUDO DISSIPATION TURBULENTE \tilde{\epsilon}
    if phase=='liq':
    
        if N_smooth > 0:
            # la longeur du plus petit vecteur dicte combien de points de trace on a 
            n_max = min(len(t),len(np.convolve( ET,filtre,mode='same'))-N_smooth)    
            # Adoucissons le signal
            temps_ET_EM_EF = np.array((t[:n_max],
                                      np.convolve(ET, filtre, mode='same')[:n_max],
                                      np.convolve(EM, filtre, mode='same')[:n_max],
                                      np.convolve(EF, filtre, mode='same')[:n_max]))
            write_me_in_raw(chemin_sauvegarde,temps_ET_EM_EF,"temps_ET_EM_EF")


            if diphasique:
                ax.plot(t[:n_max],np.convolve(Pr, filtre, mode='same')[:n_max],'-.',label=r'$\mathcal{P} : \alpha \Delta \rho g u_r \sim $'+str(np.round(mPr*10**(-E_coeff),5))+r'$\cdot 10^{%s}$'%(E_coeff))
            # ax.text(0.99*(t[-1]+t[0])/2, 0.3*mEF, r'$\mathcal{P} \sim \alpha \Delta \rho g u_r =$ '+str(np.round(mPr*10**(-E_coeff),5))+r'$\cdot 10^{%s}$'%(E_coeff), fontsize=18)
        else:
            # N'adoucissons pas le signal
            temps_ET_EM_EF = np.array((t[:],ET[:],EM[:],EF[:]))
            write_me_in_raw(chemin_sauvegarde,temps_ET_EM_EF,"temps_ET_EM_EF")
            
            if diphasique:
                ax.plot(t[:],Pr[:],'-.',label=r'$\mathcal{P} : \alpha \Delta \rho g u_r \sim $'+str(np.round(mPr*10**(-E_coeff),5))+r'$\cdot 10^{%s}$'%(E_coeff))
            # ax.text(0.99*(t[-1]+t[0])/2, 0.3*mEF, r'$\mathcal{P} \sim \alpha \Delta \rho g u_r =$ '+str(np.round(mPr*10**(-E_coeff),5))+r'$\cdot 10^{%s}$'%(E_coeff), fontsize=18)    
        
    if adim:
        temps_EF_KF = np.array((t[:],EF[:],KF[:]))
        write_me_in_raw(chemin_sauvegarde,temps_EF_KF,"temps_EF_KF")
                                   
        
    # VRAIE DISSIPATION TURBULENTE \epsilon 
    if phase=='liq':   
            
        if N_smooth > 0:
            # la longeur du plus petit vecteur dicte combien de points de trace on a 
            n_max = min(len(t),len(np.convolve( T_ET,filtre,mode='same'))-N_smooth)    
            # Adoucissons le signal
            temps_T_ET_T_EM_T_EF = np.array((t_t[:n_max],
                                            np.convolve(T_ET, filtre, mode='same')[:n_max],
                                            np.convolve(T_EM, filtre, mode='same')[:n_max],
                                            np.convolve(T_EF, filtre, mode='same')[:n_max]))
                                            
            if diphasique:
                ax.plot(t_t[:n_max],np.convolve(T_Pr, filtre, mode='same')[:n_max],'-.',label=r'$\mathcal{P} : \alpha \Delta \rho g u_r \sim $'+str(np.round(mPr*10**(-T_E_coeff),5))+r'$\cdot 10^{%s}$'%(T_E_coeff))
        else:
            # N'adoucissons pas le signal
            temps_T_ET_T_EM_T_EF = np.array((t_t[:],T_ET[:],T_EM[:],T_EF[:]))
        write_me_in_raw(chemin_sauvegarde,temps_T_ET_T_EM_T_EF,"temps_T_ET_T_EM_T_EF")    
        
        if diphasique:
            ax.plot(t_t[:],T_Pr[:],'-.',label=r'$\mathcal{P} : \alpha \Delta \rho g u_r \sim $'+str(np.round(mPr*10**(-T_E_coeff),5))+r'$\cdot 10^{%s}$'%(T_E_coeff))
        
        
    if adim:
        temps_T_EF_KF = np.array((t[:],T_EF[:],KF[:]))
        write_me_in_raw(chemin_sauvegarde,temps_T_EF_KF,"temps_T_EF_KF")

# Grandeurs affichees sur les plots 
print("####")
print("Moyenne de l'énergie cinétique des fluctuations :")
print(r'$\langle{K_F} \rangle = $ '+str(np.round(mKF*10**(-mKF_coeff),5))+r'$\cdot 10^{%s}$'%(mKF_coeff))    
print(" ==> <K_F> = ",(np.round(mKF*10**(-mKF_coeff),5))*10**(mKF_coeff))
print("Estimation de la puissance injectee par les bulles :")
if diphasique:
    print(r'$\mathcal{P} : \alpha \Delta \rho g u_r \sim $'+str(np.round(mPr*10**(-E_coeff),5))+r'$\cdot 10^{%s}$'%(E_coeff))
    print(" ==> P_inj = ",(np.round(mPr*10**(-E_coeff),5))*10**(E_coeff))
print("Estimation de la dissipation turbulente : ")
print(r'$ \epsilon = $',(np.round(E_ylabs*10**(-E_coeff),3)))
message = "\
####\n\
# Fin du script energy.py \n\
######################################################################\n"
print(message)
