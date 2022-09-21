# -*- coding: utf-8 -*-
import numpy as np
import math as mt
import matplotlib.pyplot as plt
import glob
import sys
import DNSTools3 as dtool
import datetime

""" REMARQUE : voir sur vap si le code s'execute jusqu'au bout """

"""
    Sur chaque ligne des fichiers moyenne_spatiale.txt on a la moyenne sur le plan (x,y)
    de la grandeur de la colonne j. Cette grandeur est a lire en en-tete du fichier .txt
    Chaque ligne compte 789 colonnes.
    => Verifier si c'est bien la moyenne sur le plan (x,y) ou si c'est une moyenne sur le 
    plan (y,z)...
"""

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
print("Execute le :"+str(datetime.date.today().strftime("%d %b. %Y")))
print("Execute a : "+str(datetime.datetime.now().strftime("%H:%M:%S")))
print("Commandes :"+str(sys.argv[:]))

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
def write_me_in_raw(chemin_sauvegarde,array_to_write,adim,array_name):
            the_shapes_of_array_to_write = np.array(array_to_write.shape)
            np.save(chemin_sauvegarde+"the_shapes_of_"+array_name,the_shapes_of_array_to_write)
            if adim: array_to_write.tofile(chemin_sauvegarde+"_adim_"+array_name+".raw")
            else: array_to_write.tofile(chemin_sauvegarde+"_no-adim_"+array_name+".raw")
            del(array_to_write)
            
def read_me_from_raw(array_name):
    the_shapes_of_array_to_read = np.load(chemin_sauvegarde+"the_shapes_of_"+array_name+".npy")
    print("the_shapes_of_array_to_read, ",the_shapes_of_array_to_read)
    if adim : 
        print(glob.glob(chemin_sauvegarde+"_adim_"+array_name+".raw"))
        array_to_read = np.fromfile(glob.glob(chemin_sauvegarde+"_adim_"+array_name+".raw")[0])
    else : 
        print(glob.glob(chemin_sauvegarde+"_no-adim_"+array_name+".raw"))
        print((chemin_sauvegarde+"_no-adim_"+array_name+".raw"))
        array_to_read = np.fromfile(glob.glob(chemin_sauvegarde+"_no-adim_"+array_name+".raw")[0])
    print(array_to_read.shape)
    print(the_shapes_of_array_to_read[0])
    print(the_shapes_of_array_to_read[1])
    array_to_read = array_to_read.reshape(the_shapes_of_array_to_read[0],the_shapes_of_array_to_read[1])
    list_of_arrays = []
    for i in range(the_shapes_of_array_to_read[0]):
        list_of_arrays.append(array_to_read[i][:])
    print(list_of_arrays)
    return(list_of_arrays)
################################################################################################
            
############################################################## Mode d'emploi #######################################
# mode d'emploi : python energy.py liq_vap_mel nom_des_figures fichier.data chemin_vers_diphasique_moyenne_spatiale**.txt 
#                          0              1           2                 3                                   4             
#                        fichier_out_ou_log_pour_Reb...   N_smooth  mode   [t_deb,t_fin]
#                            5                                6      7         8
# /!\ 2 : Si on met adim ou ADIM dans le nom, alors on trace les evolution de K/<K>, E/<E> en fonction de t/<T>
#         avec <T> = <K>/<E>
# 7 : # mode = draw(pas pret) / draw_from_txt / draw_from_raw / 
#              write_raw_from_txt / both(=draw_from_txt_and_write_raw_from_txt)
#              adim_draw_frow_dim_raw
# 7 : draw -> dessine les graphiques directement
#     write -> ecrits les fichiers raw
#     both -> desisne les graphiques et ecrit les fichiers raw
# /!\ liq-vap fig_ene_adim donnera nptk pou rles evoluitons de la phase vapeur. DEFINIR PLUSIEURS GRANDEURS ADIMENTIONELLES DU COUP QUAND ON MET PLUS D'UNE PHASE
#     ET QU'ON MET QUAND MEME ADIM.
############################################################## fin Mode d'emploi #######################################
message = "\
###################################################################\n\
# Script energy.py \n\
# -> Evolution temporelle de la moyenne de phase de l'énergie cinétique et de la dissipation\n\
####\n"
print(message)

## Traitement des entrees utilisateur ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
phase = sys.argv[1]
print('phase : ', phase)
if phase=='liq' : print('On n\'a pas IvdUdx ni dissip vapeur (chiv*dissip) \n grep -r \"AJOUT*.DISSIP\" dans IJK')

figname = sys.argv[2]
print('nom des figures : ', figname)

adim = ("adim".upper() in figname.upper())
if adim :
    print("################ ADIM ####################")
    print('Traces en adimentionnel : K/<K>; E/<E>; t/(<K>/<E>)')
else :
    print("############## NOT ADIM #################")

# jdd name
jdd = sys.argv[3]
print('jeu de données',jdd)

chemin_txt = sys.argv[4]
print('fichiers txt : ', chemin_txt)
chemin_sauvegarde = chemin_txt.replace("TXT/","")

# accplot
accplot = sys.argv[5]
print('fichier de sortie du premier post-trt ', accplot)

# smoothing
print("smoothing : ",sys.argv[6])
N_smooth = max(0,int(sys.argv[6]))
filtre = np.ones(N_smooth)/float(N_smooth)


# mode = draw(pas pret) / draw_from_txt / draw_from_raw / write_raw_from_txt / both(=draw_from_txt_and_write_raw_from_txt)
mode = sys.argv[7]

if "both" in mode: mode = "draw_from_txt_and_write_raw_from_txt"
draw_the_graphs = (("draw" in mode) or ("both" in mode))
write_the_raws =  (("both" in mode) or ("write" in mode))
from_txt = (("from_txt" in mode) or ("both" in mode))
from_raw = (("from_raw" in mode) or ("from_dim_raw" in mode))
adim_from_dim = ("adim_draw_from_dim_raw" in mode)
mark_repr = (("mark_repr" in mode) or ("marking_restarts" in mode))
print("#### mode ###")
print("draw_the_graphs"       ,draw_the_graphs)
print("write_the_raws"        ,write_the_raws)
print("from_txt"              ,from_txt)
print("from_raw"              ,from_raw) 
print("marking_restarts (mark_repr)"      ,mark_repr)
print(" ####### ")

# temps
print('temps',sys.argv[8])
t_deb, t_fin = float(sys.argv[8].rstrip("]").lstrip("[").split(",")[-2]),float(sys.argv[8].rstrip("]").lstrip("[").split(",")[-1])
if (t_deb <= 0.):
    t_deb = 0
if (t_fin <= 0.):
    t_fin = 1e6
print("-> t_deb : %f, t_fin : %f",t_deb,t_fin)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Lecture des parametres physiques de la simulaiton
mu = dtool.getParam(jdd, "mu_liquide")
rho = dtool.getParam(jdd, "rho_liquide")
diphasique = (not (dtool.getFlag(jdd,"disable_diphasique")))
if (diphasique):
    rho_v = dtool.getParam(jdd, "rho_vapeur")
nu = mu/rho

## Preparation nom generique des fichiers utilises
if from_txt:
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
elif from_raw:
    if (("liq" in phase) and ("vap" in phase)) : 
        fraw = glob.glob("./*Pr*"+liq+".raw")
    else : 
        fraw = glob.glob("./*Pr*"+phase+".raw")
    if fraw == [] : diphasique = False
    else : diphasique = True
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
    Reb = 0
    ur = 0
    av = 0
    g = 0
    delta_rho = 0
    estimation_Prd = 0
    print("simulation monophasique")
##### Pour récuperer les numéros de colonne qu'il faut lire dans les .txt ######
def get_n_cols():
    n_col_I = (2-1)
    n_col_uI = [(3-1,4-1,5-1) , (7-1,8-1,9-1)]
    # """interpolation 1 : uu = [ 0.5(u_i + u_{i+1/2}) ]^2"""
    n_col_uuI = [(11-1,12-1,13-1,14-1,15-1,16-1) , (229-1,230-1,231-1,232-1,233-1,234-1)]
    # """interpolation 2 : uu = 0.5(u_i^2 + u_{i+1/2}^2)"""
    # n_col_uuI = [(11-1,12-1,13-1,14-1,15-1,16-1),]
    n_col_diujI = [(72-1,73-1,74-1,75-1,76-1,77-1,78-1,79-1,80-1),
                   (706-1,707-1,708-1,709-1,710-1,711-1,712-1,713-1,714-1)]  # ordre : dudx,dvdx,dwdx,dudy,dvdy,dwdy,dudz,dvdz,dwdz
    n_col_dissip = [(228-1),(796-1)]
    try:
        n_col_true_dissip = [(797-1),(798-1)]
    except:
        print("Old version of IJK : True_dissip not post-processed directly.")
    
    return(n_col_I,n_col_uI,n_col_uuI,n_col_diujI,n_col_dissip,n_col_true_dissip)
#################################################################################
# temps
t = []
# energie cinetique
KT,KTB,KM,KF = [],[],[],[]
KT_l,KTB_l,KM_l,KF_l = [],[],[],[]
KT_v,KTB_v,KM_v,KF_v = [],[],[],[]
# dissipation visqueuse
ET,EM,EFB,EF = [],[],[],[]
ET_l,EM_l,EFB_l,EF_l = [],[],[],[]
ET_v,EM_v,EFB_v,EF_v = [],[],[],[]
T_ET,T_EM,T_EFB,T_EF = [],[],[],[]
T_ET_l,T_EM_l,T_EFB_l,T_EF_l = [],[],[],[]
T_ET_v,T_EM_v,T_EFB_v,T_EF_v = [],[],[],[]
# production / injection de qdm
Pr = [] 

##### Pour pour lire les donnees #######################################        
def recupere_les_grandeurs_de_melange_depuis_liste_de_txt(phase="mel",fstat=[],diphasique=False,rho=1000,av=0.03,delta_rho=100,g=9.81):
    """
        A partir de la liste des fichiers moyenne_spatiale.txt, retourne 
        les listes de grandeurs qui nous interessent. Les grandeurs qui 
        nous interessent sont :
         - t : la liste du temps
         - KT(t) : la liste de l'Ec => 1/2 * (U_i + u'_i)((U_i + u'_i))
         - KM(t) : la liste de l'Ec moyenne => 1/2 * (U_i)(U_i)
         - KF(t) : la liste de l'Ec fluctuante => 1/2 * (u'_i)(u'_i)
         - ET(t) : la liste de la pseudo dissipation totale => - nu * dj (U_i + u'i)
         - EM(t) : la liste de la pseudo dissipation de l'ecl moyen => - nu * dj (U_i)
         - EF(t) : la liste de la pseudo dissipation de l'ecl fluctuant => - nu * dj (u'_i)
         - T_ET(t);T_EM(t);T_EF(t) : les vraies dissipations => - 2 nu * (s_ij:s_ij) 
        if diphasique:
         - Pr(t) :  la production d'énergie due aux bulles => alpha*Delta_rho*g
    """
    # temps
    t = []
    # energie cinetique
    KT,KTB,KM,KF = [],[],[],[]
    # dissipation visqueuse
    ET,EM,EFB,EF = [],[],[],[]
    # production / injection de qdm
    Pr = []   
     
    n_col_I,n_col_uI,n_col_uuI,n_col_diujI,n_col_dissip,n_col_true_dissip = get_n_cols()
    
    if phase == "mel":
        for f in fstat:
            t.append(float("."+f.split(".")[-2]))
            sec = float(f.split(".")[-3].split("_")[-1])
            t[-1]+=sec
            #
            if ((t[-1]>=t_deb) and (t[-1]<=t_fin)): 
                # = Indicatrices
                Iv = 1-np.loadtxt(f,usecols=(n_col_I)).T
                I = np.loadtxt(f,usecols=(n_col_I)).T
                # = Vitesses et Vitesse.X.Vitesse, gradients, somme de gradients et produits de gradients
                # VItesses
                u_l = np.loadtxt(f,usecols=(n_col_uI[0])).T
                u_v = np.loadtxt(f,usecols=(n_col_uI[1])).T
                u_t = u_l+u_v
                # Vitesse.X.Vitesse
                uu_l = np.loadtxt(f,usecols=(n_col_uuI[0])).T
                uu_v = np.loadtxt(f,usecols=(n_col_uuI[1])).T
                uu_t = uu_l+uu_v
                # gradients
                diuj_l = np.loadtxt(f,usecols=(n_col_diujI[0])).T    # ordre     :  dudx,dvdx,dwdx,dudy,dvdy,dwdy,dudz,dvdz,dwdz
                diuj_v = np.loadtxt(f,usecols=(n_col_diujI[1])).T    # ordre     :  dudx,dvdx,dwdx,dudy,dvdy,dwdy,dudz,dvdz,dwdz
                diuj_t = diuj_l + diuj_v
                # Somme de gradients
                sij_l = diuj_l.reshape((3,3,diuj_l.shape[-1]))    # [i,j,...] -> diuj[...]
                sij_v = diuj_v.reshape((3,3,diuj_v.shape[-1]))    # [i,j,...] -> diuj[...]
                sij_l = 0.5*(sij_l+sij_l.transpose(1,0,2))        # sij = djui + diuj
                sij_v = 0.5*(sij_v+sij_v.transpose(1,0,2))        # sij = djui + diuj
                sij_t = sij_l + sij_v
                # Somme de produits de gradients
                dissip_l = np.loadtxt(f,usecols=(n_col_dissip[0])).T
                dissip_v = np.loadtxt(f,usecols=(n_col_dissip[1])).T
                dissip_t = dissip_l + dissip_v
                try:
                    true_dissip_l = np.loadtxt(f,usecols=(n_col_true_dissip[0])).T
                    true_dissip_v = np.loadtxt(f,usecols=(n_col_true_dissip[1])).T
                except:
                    print("Old IJK : true_dissip is not written...")
                    true_dissip_l = -2*np.einsum("ij...,ij...->...",sij_l,sij_l)
                    print("Constructed for liq. phase.")
                    true_dissip_v = -2*np.einsum("ij...,ij...->...",sij_v,sij_v)
                    print("Constructed for vap. phase.")
                true_dissip_t = true_dissip_l + true_dissip_v 
                # = Moyennes (par rapport aux tranches)
                m_I = np.mean(I,axis=-1)
                m_Iv = np.mean(Iv,axis=-1)
                m_u_l = np.mean(u_l,axis=-1) # essentiel pour la vitesse relative
                m_u_v = np.mean(u_v,axis=-1) # essentiel pour la vitesse relative
                m_uu_t = np.mean(uu_t,axis=-1)
                m_u_t = np.mean(u_t,axis=-1)
                m_diuj_t = np.mean(diuj_t,axis=-1)
                m_sij_t = np.mean(sij_t,axis=-1)
                m_dissip_t = np.mean(dissip_t,axis=-1)
                m_true_dissip_t = np.mean(true_dissip_t,axis=-1)
                # = Moyennes de phase (pour la vitesse relative)
                m_u_ll = m_u_l/m_I
                m_u_vv = m_u_v/m_Iv
                # = /!\ Construction tension de Reynolds
                u_t_u_t = np.einsum("...i,...j->...ij",m_u_t,m_u_t)
                u_t_u_t = np.array((u_t_u_t[0][0],u_t_u_t[1][1],u_t_u_t[2][2],u_t_u_t[0][1],u_t_u_t[1][2],u_t_u_t[0][2]))
                m_Rij_t = m_uu_t - u_t_u_t
                # = Energies cinetique
                k_m = 0.5*(u_t_u_t[0] + u_t_u_t[1] + u_t_u_t[2])
                k_t = 0.5*(m_uu_t[0] + m_uu_t[1] + m_uu_t[2])
                k_f = 0.5*(m_Rij_t[0] + m_Rij_t[1] + m_Rij_t[2])
                k_t_bis = k_f + k_m
                # = Dissipations
                e_t = -1*nu*(m_dissip_t)
                true_e_t = -1*nu*(m_true_dissip_t) 
                e_m = -1*nu*(np.einsum("ji...,ji->...",m_diuj_t.reshape(3,3),m_diuj_t.reshape(3,3)))
                true_e_m = -1*nu*(np.einsum("ij...,ij...->...",m_sij_t,m_sij_t))  # double-produit contracte de sij moyen.
                e_f = e_t-e_m
                true_e_f = true_e_t-true_e_m
                # = Production d'Energie cinetique due aux bulles
                u_r = np.sqrt((m_u_vv[0]-m_u_ll[0])**2+(m_u_vv[1]-m_u_ll[1])**2+(m_u_vv[2]-m_u_ll[2])**2)
                prd = (1./rho)*av*delta_rho*g*u_r
                # = Appending to the lists
                KT.append(k_t);KM.append(k_m);KF.append(k_f);KTB.append(k_t_bis);
                ET.append(e_t);EM.append(e_m);EF.append(e_f);
                T_ET.append(true_e_t);T_EM.append(true_e_m);T_EF.append(true_e_f);
                if diphasique:
                    prd = (1./rho)*av*delta_rho*g*u_r                
                    Pr.append(prd)        
            else:
                t.pop()  
    # print("m_I : ",m_I)
    return(t,KT,KTB,KM,KF,ET,EM,EF,T_ET,T_EM,T_EF,Pr)

def recupere_les_grandeurs_de_phase_depuis_liste_de_txt(phase="vap",fstat=[],diphasique=False,rho=1000,av=0.03,delta_rho=100,g=9.81):
    """
        A partir de la liste des fichiers moyenne_spatiale.txt, retourne 
        les listes de grandeurs qui nous interessent. Les grandeurs qui 
        nous interessent sont :
         - t : la liste du temps
         - KT(t) : la liste de l'Ec => 1/2 * (U_i + u'_i)((U_i + u'_i))
         - KM(t) : la liste de l'Ec moyenne => 1/2 * (U_i)(U_i)
         - KF(t) : la liste de l'Ec fluctuante => 1/2 * (u'_i)(u'_i)
         - ET(t) : la liste de la pseudo dissipation totale => - nu * dj (U_i + u'i)
         - EM(t) : la liste de la pseudo dissipation de l'ecl moyen => - nu * dj (U_i)
         - EF(t) : la liste de la pseudo dissipation de l'ecl fluctuant => - nu * dj (u'_i)
         - T_ET(t);T_EM(t);T_EF(t) : les vraies dissipations => - 2 nu * (s_ij:s_ij) 
        if diphasique:
         - Pr(t) :  la production d'énergie due aux bulles => alpha*Delta_rho*g
    """
    
    # temps
    t = []
    # energie cinetique
    KT_l,KTB_l,KM_l,KF_l = [],[],[],[]
    KT_v,KTB_v,KM_v,KF_v = [],[],[],[]
    # dissipation visqueuse
    ET_l,EM_l,EFB_l,EF_l = [],[],[],[]
    ET_v,EM_v,EFB_v,EF_v = [],[],[],[]
    T_ET_l,T_EM_l,T_EFB_l,T_EF_l = [],[],[],[]
    T_ET_v,T_EM_v,T_EFB_v,T_EF_v = [],[],[],[]
    # production / injection de qdm
    Pr = [] 
    
    n_col_I,n_col_uI,n_col_uuI,n_col_diujI,n_col_dissip,n_col_true_dissip = get_n_cols()

    if ("liq" in phase) or ("vap" in phase):
        for f in fstat:
            t.append(float("."+f.split(".")[-2]))
            sec = float(f.split(".")[-3].split("_")[-1])
            t[-1]+=sec
            #
            if ((t[-1]>=t_deb) and (t[-1]<=t_fin)): 
                # ~~~ Vitesse relative ~~~~~~~~~~~~~~~~~~~~~
                I = np.loadtxt(f,usecols=(n_col_I)).T
                Iv = 1-I
                u_l, u_v = np.loadtxt(f,usecols=(n_col_uI[0])).T,np.loadtxt(f,usecols=(n_col_uI[1])).T
                m_I, m_Iv = np.mean(I,axis=-1),np.mean(Iv,axis=-1)
                m_u_l, m_u_v = np.mean(u_l,axis=-1), np.mean(u_v,axis=-1)            
                m_u_ll, m_u_vv = m_u_l/m_I, m_u_v/m_Iv
                u_r = np.sqrt((m_u_vv[0]-m_u_ll[0])**2+(m_u_vv[1]-m_u_ll[1])**2+(m_u_vv[2]-m_u_ll[2])**2)
                # ~~~ Liquide ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                if "liq" in phase:
                    # = Indicatrices : deja fait
                    # = Vitesses, Vitesse.X.Vitesse, tenseur de deformations, produits de gradients, somme de produits de somme de gradients
                    # Vitesse : deja fait
                    # VItesse X Vitesse
                    uu_l = np.loadtxt(f,usecols=(n_col_uuI[0])).T
                    # gradients
                    diuj_l = np.loadtxt(f,usecols=(n_col_diujI[0])).T    # ordre     :  dudx,dvdx,dwdx,dudy,dvdy,dwdy,dudz,dvdz,dwdz
                    # tenseur des deformations
                    sij_l = diuj_l.reshape((3,3,diuj_l.shape[-1]))    # [i,j,...] -> diuj[...]
                    sij_l = 0.5*(sij_l+sij_l.transpose(1,0,2))        # sij = djui + diuj
                    # somme de produits de gradients
                    dissip_l = np.loadtxt(f,usecols=(n_col_dissip[0])).T
                    try:
                        true_dissip_l = np.loadtxt(f,usecols=(n_col_true_dissip[0])).T
                    except:
                        print("Old IJK : true_dissip is not written...")
                        true_dissip_l = -2*np.einsum("ij...,ij...->...",sij_l,sij_l)
                        print("Constructed for liq. phase.")
                    # = Moyennes
                    #"m_I = np.mean(I,axis=-1)
                    m_uu_l = np.mean(uu_l,axis=-1)
                    #"m_u_l = np.mean(u_l,axis=-1)
                    m_diuj_l = np.mean(diuj_l,axis=-1)
                    m_sij_l = np.mean(sij_l,axis=-1)
                    m_dissip_l = np.mean(dissip_l,axis=-1)            
                    m_true_dissip_l = np.mean(true_dissip_l,axis=-1)            
                    # = Moyennes de phase (inutil en melange)
                    #"m_u_ll = m_u_l/m_I
                    m_uu_ll = m_uu_l/m_I
                    m_dissip_ll = m_dissip_l/m_I
                    m_diuj_ll = m_diuj_l/m_I
                    m_sij_ll = m_sij_l/m_I
                    m_true_dissip_ll = m_true_dissip_l/m_I
                    # = /!\ Construction des tensions de Reynolds
                    u_l_u_l = np.einsum("...i,...j->...ij",m_u_ll,m_u_ll)
                    u_l_u_l = np.array((u_l_u_l[0][0],u_l_u_l[1][1],u_l_u_l[2][2],u_l_u_l[0][1],u_l_u_l[1][2],u_l_u_l[0][2]))
                    m_Rij_l = m_uu_ll - u_l_u_l
                    # = Energies cinetique
                    k_m_l = 0.5*(u_l_u_l[0] + u_l_u_l[1] + u_l_u_l[2])
                    k_t_l = 0.5*(m_uu_ll[0] + m_uu_ll[1] + m_uu_ll[2])
                    k_f_l = 0.5*(m_Rij_l[0] + m_Rij_l[1] + m_Rij_l[2])
                    k_t_bis_l = k_f_l + k_m_l
                    # = Dissipations
                    e_t_l = -1*nu*(m_dissip_ll)            # (djui.djui)
                    true_e_t_l = -1*nu*(m_true_dissip_ll)  # (djui.djui) + (diuj.djui) = 2 sij:sij
                    e_m_l = -1*nu*(np.einsum("ji...,ji->...",m_diuj_ll.reshape(3,3),m_diuj_ll.reshape(3,3)))
                    true_e_m_l = -2*nu*(np.einsum("ij...,ij...->...",m_sij_ll,m_sij_ll))  # double-produit contracte de sij moyen.
                    e_f_l = e_t_l - e_m_l
                    true_e_f_l = true_e_t_l - true_e_m_l
                    # = Appending to the lists
                    KT_l.append(k_t_l);KM_l.append(k_m_l);KF_l.append(k_f_l);KTB_l.append(k_t_bis_l);
                    ET_l.append(e_t_l);EM_l.append(e_m_l);EF_l.append(e_f_l);
                    T_ET_l.append(true_e_t_l);T_EM_l.append(true_e_m_l);T_EF_l.append(true_e_f_l);
                    
                # ~~~ Vapeur ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                elif "vap" in phase:
                    print("Impossible for vap. phase : I_v du_idx_j missing in IJK")
                    # = Indicatrices : deja fait
                    # = Vitesses, Vitesse.X.Vitesse, tenseur de deformations, produits de gradients, somme de produits de somme de gradients
                    # Vitesse : deja fait
                    # VItesse X Vitesse
                    uu_v = np.loadtxt(f,usecols=(n_col_uuI[1])).T
                    # gradients
                    diuj_v = np.loadtxt(f,usecols=(n_col_diujI[1])).T    # ordre     :  dudx,dvdx,dwdx,dudy,dvdy,dwdy,dudz,dvdz,dwdz
                    # tenseur des deformations
                    sij_v = diuj_v.reshape((3,3,diuj_v.shape[-1]))    # [i,j,...] -> diuj[...]
                    sij_v = 0.5*(sij_v+sij_v.transpose(1,0,2))        # sij = djui + diuj
                    # somme de produits de gradients
                    dissip_v = np.loadtxt(f,usecols=(n_col_dissip[0])).T
                    try:
                        true_dissip_v = np.loadtxt(f,usecols=(n_col_true_dissip[0])).T
                    except:
                        print("Old IJK : true_dissip is not written...")
                        true_dissip_v = -2*np.einsum("ij...,ij...->...",sij_v,sij_v)
                        print("Constructed for vap. phase.")
                    # = Moyennes
                    #"m_Iv = np.mean(Iv,axis=-1)
                    m_uu_v = np.mean(uu_v,axis=-1)
                    m_diuj_v = np.mean(diuj_v,axis=-1)
                    m_dissip_v = np.mean(dissip_v,axis=-1)
                    m_true_dissip_v = np.mean(true_dissip_v,axis=-1)
                    m_sij_v = np.mean(sij_v,axis=-1)
                    #"m_u_v = np.mean(u_v,axis=-1)            
                    # = Moyennes de phase (inutil en melange)
                    #"m_u_vv = m_u_v/m_Iv
                    m_uu_vv = m_uu_v/m_Iv
                    m_dissip_vv = m_dissip_v/m_Iv
                    m_diuj_vv = m_diuj_v/m_Iv
                    m_sij_vv = m_sij_v/m_Iv
                    m_true_dissip_vv = m_true_dissip_v/m_Iv
                    # = /!\ Construction des tensions de Reynolds
                    u_v_u_v = np.einsum("...i,...j->...ij",m_u_vv,m_u_vv)
                    u_v_u_v = np.array((u_v_u_v[0][0],u_v_u_v[1][1],u_v_u_v[2][2],u_v_u_v[0][1],u_v_u_v[1][2],u_v_u_v[0][2]))
                    m_Rij_v = m_uu_vv - u_v_u_v
                    # = Energies cinetique
                    k_m_v = 0.5*(u_v_u_v[0] + u_v_u_v[1] + u_v_u_v[2])
                    k_t_v = 0.5*(m_uu_vv[0] + m_uu_vv[1] + m_uu_vv[2])
                    k_f_v = 0.5*(m_Rij_v[0] + m_Rij_v[1] + m_Rij_v[2])
                    k_t_bis_v = k_f_v + k_m_v
                    # = Dissipations
                    e_t_v = -1*nu*(m_dissip_vv)            # (djui.djui)
                    true_e_t_v = -1*nu*(m_true_dissip_vv)  # (djui.djui) + (diuj.djui) = 2 sij:sij
                    e_m_v = -1*nu*(np.einsum("ji...,ji->...",m_diuj_vv.reshape(3,3),m_diuj_vv.reshape(3,3)))
                    true_e_m_v = -2*nu*(np.einsum("ij...,ij...->...",m_sij_vv,m_sij_vv))  # double-produit contracte de sij moyen.
                    e_f_v = e_t_v - e_m_v
                    true_e_f_v = true_e_t_v - true_e_m_v
                    # = Appending to the lists
                    KT_v.append(k_t_v);KM_v.append(k_m_v);KF_v.append(k_f_v);KTB_v.append(k_t_bis_v);
                    ET_v.append(e_t_v);EM_v.append(e_m_v);EF_v.append(e_f_v);
                    T_ET_v.append(true_e_t_v);T_EM_v.append(true_e_m_v);T_EF_v.append(true_e_f_v);
                
                # = Production d'Energie cinetique due aux bulles
                if diphasique:
                    prd = (1./rho)*av*delta_rho*g*u_r                
                    Pr.append(prd)
            else:
                t.pop()
    return( t,
            KT_l,KTB_l,KM_l,KF_l,ET_l,EM_l,EF_l,T_ET_l,T_EM_l,T_EF_l,
            KT_v,KTB_v,KM_v,KF_v,ET_v,EM_v,EF_v,T_ET_v,T_EM_v,T_EF_v,
            Pr)

def recupere_les_grandeurs_depuis_raw_et_npy(phase="mel"):
    # temps
    t = []
    # energie cinetique
    KT,KTB,KM,KF = [],[],[],[]
    KT_l,KTB_l,KM_l,KF_l = [],[],[],[]
    KT_v,KTB_v,KM_v,KF_v = [],[],[],[]
    # dissipation visqueuse
    ET,EM,EFB,EF = [],[],[],[]
    ET_l,EM_l,EFB_l,EF_l = [],[],[],[]
    ET_v,EM_v,EFB_v,EF_v = [],[],[],[]
    T_ET,T_EM,T_EFB,T_EF = [],[],[],[]
    T_ET_l,T_EM_l,T_EFB_l,T_EF_l = [],[],[],[]
    T_ET_v,T_EM_v,T_EFB_v,T_EF_v = [],[],[],[]
    # production / injection de qdm
    Pr, Pr_v, Pr_l = [], [], [] 
    
    if "mel" in phase:
        t,KT,KTB,KM,KF   = read_me_from_raw("temps_KT_KTB_KM_KF_"+"mel")
        _,ET,EM,EF       = read_me_from_raw("temps_ET_EM_EF_"+"mel")
        _,T_ET,T_EM,T_EF = read_me_from_raw("temps_T_ET_T_EM_T_EF_"+"mel")
        if diphasique:
            if adim:
                _, Pr = read_me_from_raw("temps_Pr_"+"mel")
                _, T_Pr = read_me_from_raw("temps_T_Pr_"+"mel")
            else:
                _, Pr = read_me_from_raw("temps_Pr_"+"mel")
    
    if "liq" in phase:
        t,  KT_l, KTB_l,  KM_l, KF_l = read_me_from_raw("temps_KT_KTB_KM_KF_"+"liq")
        _,  ET_l,  EM_l,  EF_l       = read_me_from_raw("temps_ET_EM_EF_"+"liq")
        _,T_ET_l,T_EM_l,T_EF_l       = read_me_from_raw("temps_T_ET_T_EM_T_EF_"+"liq")
        if diphasique:
            if adim :
                _, Pr_l = read_me_from_raw("temps_Pr_"+"liq")            
                _, T_Pr_l = read_me_from_raw("temps_T_Pr_"+"liq")            
            else :
                _, Pr = read_me_from_raw("temps_Pr_"+"liq")
    if "vap" in phase:
        t,  KT_v, KTB_v,  KM_v, KF_v = read_me_from_raw("temps_KT_KTB_KM_KF_"+"vap")
        _,  ET_v,  EM_v,  EF_v       = read_me_from_raw("temps_ET_EM_EF_"+"vap")
        _,T_ET_v,T_EM_v,T_EF_v       = read_me_from_raw("temps_T_ET_T_EM_T_EF_"+"vap")
        if diphasique:
            if adim :
                _, Pr_v = read_me_from_raw("temps_Pr_"+"vap")
                _, T_Pr_v = read_me_from_raw("temps_T_Pr_"+"vap")
            else:
                _, Pr = read_me_from_raw("temps_Pr_"+"vap")
    return (t,
            KT,KTB,KM,KF, 
            KT_l,KTB_l,KM_l,KF_l,
            KT_v,KTB_v,KM_v,KF_v,
            ET,EM,EFB,EF,
            ET_l,EM_l,EFB_l,EF_l,
            ET_v,EM_v,EFB_v,EF_v,
            T_ET,T_EM,T_EFB,T_EF,
            T_ET_l,T_EM_l,T_EFB_l,T_EF_l,
            T_ET_v,T_EM_v,T_EFB_v,T_EF_v,
            Pr, Pr_v, Pr_l )

def iterations_et_instants_de_reprise(case):
    """
        Renvoi, a prtir du dt_ev :
         - les iterations correspondant a une reprise
         - les instants correspondants a une reprise.
    """
    print(jdd+"*dt_ev")
    print(glob.glob(jdd+"*dt_ev"))
    it = np.loadtxt(glob.glob(case+"*dt_ev")[0],usecols=([0, 1]))
    i,t = it[:,0], it[:,1]
    i = np.argwhere(i==0)
    t = t[i]
    return(i,t)
########################################################################
        
        

# LECTURE DES FICHIERS ET CONSTRUCITON DES ENERGIE, DISSIPATIONS #######
if from_txt:
    if "mel" in phase:
        (   t,
            KT,KTB,KM,KF,ET,EM,EF,T_ET,T_EM,T_EF,
            Pr
        ) = recupere_les_grandeurs_de_melange_depuis_liste_de_txt(phase=phase,fstat=fstat,diphasique=diphasique,rho=rho,av=av,delta_rho=delta_rho,g=g)
    else:
        (   t,
            KT_l,KTB_l,KM_l,KF_l,ET_l,EM_l,EF_l,T_ET_l,T_EM_l,T_EF_l,
            KT_v,KTB_v,KM_v,KF_v,ET_v,EM_v,EF_v,T_ET_v,T_EM_v,T_EF_v,
            Pr
        ) = recupere_les_grandeurs_de_phase_depuis_liste_de_txt(phase=phase,fstat=fstat,diphasique=diphasique,rho=rho,av=av,delta_rho=delta_rho,g=g)
elif from_raw:
    (t,
            KT,KTB,KM,KF, 
            KT_l,KTB_l,KM_l,KF_l,
            KT_v,KTB_v,KM_v,KF_v,
            ET,EM,EFB,EF,
            ET_l,EM_l,EFB_l,EF_l,
            ET_v,EM_v,EFB_v,EF_v,
            T_ET,T_EM,T_EFB,T_EF,
            T_ET_l,T_EM_l,T_EFB_l,T_EF_l,
            T_ET_v,T_EM_v,T_EFB_v,T_EF_v,
            Pr, Pr_v, Pr_l 
    ) = recupere_les_grandeurs_depuis_raw_et_npy(phase)
    if adim and diphasique and not("mel" in phase) :
        if "liq" in phase : Pr = Pr_l
        elif "vap" in phase : Pr = Pr_v
print("diphasique : ",diphasique)
print("Pr : ",Pr)

if mark_repr:
    i_repr,t_repr = iterations_et_instants_de_reprise(jdd.replace(".data",""))
########################################################################

##################### Pour les  Plots #########################
#Objectif : gérer l'affichage des ylabel en fonction des valeurs de KF
xlabs = np.floor(100*(np.linspace(min(t),max(t),7)))*10**(-2)


# ~~~~~~~ Energie cinetique ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if "mel" in phase:
    KF_coeff = np.floor(mt.log10(max(KF)))
    KF_ylabs = np.linspace(min(KF),max(KF),7)
    KF_xlabs = xlabs
    
    K_coeff = np.floor(mt.log10(max(KT)))
    K_ylabs = np.linspace(min(min(KF),min(KT),min(KM)),max(max(KT),max(KF),max(KM)),7)
    
    mKF = np.abs(np.mean(np.array(KF[-int(len(KF)/2):]))) # moyenne de l'Ec de fluctuations : <u'.u'>
    mKF_coeff = np.floor(mt.log10(mKF))
if "liq" in phase:
    KF_coeff_l = np.floor(mt.log10(max(KF_l)))
    KF_ylabs_l = np.linspace(min(KF_l),max(KF_l),7)
    
    K_coeff_l = np.floor(mt.log10(max(KT_l)))
    K_ylabs_l = np.linspace(min(min(KF_l),min(KT_l),min(KM_l)),max(max(KT_l),max(KF_l),max(KM_l)),7)
    
    mKF_l = np.abs(np.mean(np.array(KF_l[-int(len(KF_l)/2):]))) # moyenne de l'Ec de fluctuations : <u'.u'>
    mKF_coeff_l = np.floor(mt.log10(mKF_l))
if "vap" in phase:    
    KF_coeff_v = np.floor(mt.log10(max(KF_v)))
    KF_ylabs_v = np.linspace(min(KF_v),max(KF_v),7)
    
    K_coeff_v = np.floor(mt.log10(max(KT_v)))
    K_ylabs_v = np.linspace(min(min(KF_v),min(KT_v),min(KM_v)),max(max(KT_v),max(KF_v),max(KM_v)),7)
    
    mKF_v = np.abs(np.mean(np.array(KF_v[-int(len(KF_v)/2):]))) # moyenne de l'Ec de fluctuations : <u'.u'>
    mKF_coeff_v = np.floor(mt.log10(mKF_v))

# ~~~~~~~ Dissipation et affichage dans le terminal ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if diphasique :
    print("Pr =\t",Pr)
if "mel" in phase:
    E_coeff = np.floor(mt.log10(max([abs(max(EM)),abs(max(EF))])))
    T_E_coeff = np.floor(mt.log10(max([abs(max(T_EM)),abs(max(T_EF))])))
    E_ylabs = np.linspace(min([min(EF),min(EM)]),max([max(EF),max(EM)]),7)
    if diphasique:
        T_E_ylabs = np.linspace(min([min(T_EF),min(T_EM),min(Pr)]),max([max(T_EF),max(T_EM),max(Pr)]),7)
    else:
        T_E_ylabs = np.linspace(min([min(T_EF),min(T_EM)]),max([max(T_EF),max(T_EM)]),7)
    E_xlabs = KF_xlabs
    T_E_xlabs = KF_xlabs
    mEF = np.mean(np.array(EF[-int(len(EF)/2):])) # moyenne de la dissipations due aux fluctuations : < \ol{s'_ij. s'_ij} >
    mT_EF = np.mean(np.array(T_EF[-int(len(T_EF)/2):])) # moyenne de la true_dissipations due aux fluctuations : < \ol{s'_ij. s'_ij} >
    print("maximum dissipation totale (pseudo et true) : ",max(ET),max(T_ET))
    print("maximum dissipation fluctu (pseudo et true) : ",max(EF),max(T_EF))
    print("mEF,m_True_EF")

    print(mEF,mT_EF)

    print("KT =\t",KT)
    print("KM =\t",KM)
    print("KF =\t",KF)
    print("KTB= \t",KTB)
    
    print("### pseudo dissipations")
    print("ET =\t ",ET)
    print("EM =\t ",EM)
    print("EF =\t ",EF)
    print("### true dissipations")
    print("T_ET =\t ",T_ET)
    print("T_EM =\t ",T_EM)
    print("T_EF =\t ",T_EF)        
if "liq" in phase:
    E_coeff_l = np.floor(mt.log10(max([abs(max(EM_l)),abs(max(EF_l))])))
    T_E_coeff_l = np.floor(mt.log10(max([abs(max(T_EM_l)),abs(max(T_EF_l))])))
    if diphasique:
        E_ylabs_l =   np.linspace(min([min(EF_l),  min(EM_l), min(Pr)]),max([max(EF_l),max(EM_l),max(Pr)]),7)
        T_E_ylabs_l = np.linspace(min([min(T_EF_l),min(T_EM_l),min(Pr)]),max([max(T_EF_l),max(T_EM_l),max(Pr)]),7)
    else:
        E_ylabs_l =   np.linspace(min([min(EF_l),  min(EM_l)]),max([max(EF_l),max(EM_l)]),7)
        T_E_ylabs_l = np.linspace(min([min(T_EF_l),min(T_EM_l)]),max([max(T_EF_l),max(T_EM_l)]),7)
    mEF_l = np.mean(np.array(EF_l[-int(len(EF_l)/2):])) # moyenne de la dissipations due aux fluctuations : < \ol{s'_ij. s'_ij} >
    mT_EF_l = np.mean(np.array(T_EF_l[-int(len(T_EF_l)/2):])) # moyenne de la true_dissipations due aux fluctuations : < \ol{s'_ij. s'_ij} >
    print("maximum dissipation totale (pseudo et true) : ",max(ET_l),max(T_ET_l))
    print("maximum dissipation fluctu (pseudo et true) : ",max(EF_l),max(T_EF_l))
    print("mEF,m_True_EF")

    print(mEF_l,mT_EF_l)

    print("KT =\t",KT_l)
    print("KM =\t",KM_l)
    print("KF =\t",KF_l)
    print("KTB= \t",KTB_l)
    
    print("### pseudo dissipations")
    print("ET =\t ",ET_l)
    print("EM =\t ",EM_l)
    print("EF =\t ",EF_l)
    print("### true dissipations")
    print("T_ET =\t ",T_ET_l)
    print("T_EM =\t ",T_EM_l)
    print("T_EF =\t ",T_EF_l)
if "vap" in phase:
    E_coeff_v = np.floor(mt.log10(max([abs(max(EM_v)),abs(max(EF_v))])))
    print("T_EM_v,T_EF_v")
    print(T_EM_v,T_EF_v)
    T_E_coeff_v = np.floor(mt.log10(max([abs(max(T_EM_v)),abs(max(T_EF_v))])))
    E_ylabs_v = np.linspace(min([min(EF_v),min(EM_v)]),max([max(EF_v),max(EM_v)]),7)
    if diphasique:
        T_E_ylabs_v = np.linspace(min([min(T_EF_v),min(T_EM_v),min(Pr)]),max([max(T_EF_v),max(T_EM_v),max(Pr)]),7)
    else:
        T_E_ylabs_v = np.linspace(min([min(T_EF_v),min(T_EM_v)]),max([max(T_EF_v),max(T_EM_v)]),7)
    mEF_v = np.mean(np.array(EF_v[-int(len(EF_v)/2):])) # moyenne de la dissipations due aux fluctuations : < \ol{s'_ij. s'_ij} >
    mT_EF_v = np.mean(np.array(T_EF_v[-int(len(T_EF_v)/2):])) # moyenne de la true_dissipations due aux fluctuations : < \ol{s'_ij. s'_ij} >
    print("maximum dissipation totale (pseudo et true) : ",max(ET_v),max(T_ET_v))
    print("maximum dissipation fluctu (pseudo et true) : ",max(EF_v),max(T_EF_v))
    print("mEF,m_True_EF")

    print(mEF_v,mT_EF_v)

    print("KT =\t",KT_v)
    print("KM =\t",KM_v)
    print("KF =\t",KF_v)
    print("KTB= \t",KTB_v)
    
    print("### pseudo dissipations")
    print("ET =\t ",ET_v)
    print("EM =\t ",EM_v)
    print("EF =\t ",EF_v)
    print("### true dissipations")
    print("T_ET =\t ",T_ET_v)
    print("T_EM =\t ",T_EM_v)
    print("T_EF =\t ",T_EF_v)
if diphasique:
    mPr = np.mean(np.array(Pr[-int(len(Pr)/2):])) # estimation de la dissipation due a la flotabilite : \a \D\r g u_r

if (from_raw and adim) :
    """
    Si on adimensionne, alors on converge vers E_adim = K_adim = 1
    Ce n'est pas ce qu'on attends.
    Il faut donc aller lire les coeffs dans les .npy ecrit auparavent.
    Les labels correspondent a ce qui est attendu en revanche
    """
    if "mel" in phase:
        K_coeff, K_coeff, mKF, mKF_coeff = np.load(chemin_sauvegarde+"KFco__Kco__mKF__mKFco"+"mel"+".npy")
        E_coeff, T_E_coeff, mEF, mT_EF = np.load(chemin_sauvegarde+"Eco__TEco__mEF__mTEF"+"mel"+".npy")
    if "liq" in phase:
        K_coeff_l, K_coeff_l, mKF_l, mKF_coeff_l = np.load(chemin_sauvegarde+"KFco__Kco__mKF__mKFco"+"liq"+".npy")
        E_coeff_l, T_E_coeff_l, mEF_l, mT_EF_l = np.load(chemin_sauvegarde+"Eco__TEco__mEF__mTEF"+"liq"+".npy")
    if "vap" in phase:
        K_coeff_v, K_coeff_v, mKF_v, mKF_coeff_v = np.load(chemin_sauvegarde+"KFco__Kco__mKF__mKFco"+"vap"+".npy")
        E_coeff_v, T_E_coeff_v, mEF_v, mT_EF_v = np.load(chemin_sauvegarde+"Eco__TEco__mEF__mTEF"+"vap"+".npy")
    

# print("min(KF),max(KT),KF_coeff,min(t),max(t)")
# print(min(KF),max(KT),KF_coeff,min(t),max(t))
########################################################

### REPENSER L'ECRITURE POUR QUE CE SOIT FACILEMENT TRACABLE AVEC GNUPLOT,
### OU AVEC PYHTON EN LECTURE DIRECTE.

# print("EFB",EFB)

if from_txt or adim_from_dim:
    """
    Travail d'adimensionnement
    """
    print("t =\t ",t)
    
    t_t = t
    if diphasique:
        T_Pr = Pr
        
    if adim:
        if ("liq" in phase):
            T_L = mKF_l/np.abs(mEF_l)
            T_T_L = mKF_l/np.abs(mT_EF_l)
            t = t/T_L
            t_t = t_t/T_T_L
            KT_l,KM_l,KF_l,KTB_l = KT_l/mKF_l,KM_l/mKF_l,KF_l/mKF_l,KTB_l/mKF_l
            ET_l,EM_l,EF_l = ET_l/mEF_l,EM_l/mEF_l,EF_l/mEF_l
            T_ET_l,T_EM_l,T_EF_l = T_ET_l/mT_EF_l,T_EM_l/mT_EF_l,T_EF_l/mT_EF_l
            if diphasique:
                Pr_l = Pr/mEF_l
                T_Pr_l = T_Pr/mT_EF_l
        elif ("vap" in phase):
            T_L = mKF_v/np.abs(mEF_v)
            T_T_L = mKF_v/np.abs(mT_EF_v)
            t = t/T_L
            t_t = t_t/T_T_L
            KT_v,KM_v,KF_v,KTB_v = KT_v/mKF_v,KM_v/mKF_v,KF_v/mKF_v,KTB_v/mKF_v
            ET_v,EM_v,EF_v = ET_v/mEF_v,EM_v/mEF_v,EF_v/mEF_v
            T_ET_v,T_EM_v,T_EF_v = T_ET_v/mT_EF_v,T_EM_v/mT_EF_v,T_EF_v/mT_EF_v
            if diphasique:
                Pr_v = Pr/mEF_v
                T_Pr_v = T_Pr/mT_EF_v
        elif ("mel" in phase):        
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
if from_raw and adim:
    """ 
    Si les raw sont ecrits en non adim et qu'on on plot en non adim
    ou
    Si les raw sont ecrits en adim et qu'on plot en adim
    Alors, il n'y a pas besoin de repasser dnas les adimensionnements ci-dessus.
    """
    T_L, T_T_L = np.load(chemin_sauvegarde+"TL__TTL.npy")


#################################################################################
# TRACE DES FIGURES OU ECRITURE DES DONNEES - ENERGIES CINETIQUES
#################################################################################
print("mode : "+str(mode))
#################################################################################
def get_adim_instants_reprise(instants_reprises, adim):
    if adim : instants_reprises = instants_reprises/T_L
    return instants_reprises

if mark_repr: t_repr = get_adim_instants_reprise(t_repr,adim)

def add_repr_vlines(ax,instant_reprises):
    for t in instant_reprises:
        ax.axvline(x=t,linewidth=0.3,color='black',linestyle='dotted')
#################################################################################
    
def get_K_label(phase):
    if "mel" in phase:
        mF,mF_coeff = mKF,mKF_coeff
    elif "vap" in phase:
        mF,mF_coeff = mKF_v,mKF_coeff_v
    elif "liq" in phase:        
        mF,mF_coeff = mKF_l,mKF_coeff_l
    label_KT  = r'total : $ \frac{1}{2} (U_i \cdot U_i)$'
    label_KTB = r"m+f : $\frac{1}{2} (\overline{u} \cdot \overline{u} + u' \cdot u' )$"
    label_KM  = r'moyen : $ \frac{1}{2} (\overline{u} \cdot \overline{u})$'
    label_KF  = r"fluct : $ \frac{1}{2} (u' \cdot u'$); " + r'$\langle{K_F} \rangle = $ '+str(np.round(mF*10**(-mF_coeff),2))+r'$\cdot 10^{%s}$'%(mF_coeff)
    print("Energie fluctuations = "+str(mF)+", (same) K_F = "+str(mF))

    return(label_KT,label_KTB,label_KM,label_KF)
# ~~~~~~ ENERGIE CINETIQUE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def do_the_K_plot(phase):
    
    label_T,label_TB,label_M,label_F = get_K_label(phase)
    if "mel" in phase:
        TOT,TOT_bis,MOY,FLC,YLABS,COEFF = KT,KTB,KM,KF,K_ylabs,K_coeff
    elif "liq" in phase:         
        TOT,TOT_bis,MOY,FLC,YLABS,COEFF = KT_l,KTB_l,KM_l,KF_l,K_ylabs_l,K_coeff_l
    elif "vap" in phase:         
        TOT,TOT_bis,MOY,FLC,YLABS,COEFF = KT_v,KTB_v,KM_v,KF_v,K_ylabs_v,K_coeff_v
    
    fig, ax = plt.subplots(1,figsize=(9,9))
    ax.set_title('Energies cinetiques - moyenne de phase '+phase+r' : $\overline{k}^k = \frac{\overline{\chi_k u \cdot u}}{\overline{\chi_k}}$',fontsize=18)
    if N_smooth > 0:
        # la longeur du plus petit vecteur dicte combien de points de trace on a 
        n_max = min(len(t),len(np.convolve( TOT,filtre,mode='same'))-N_smooth)
        # On sait que les energies et dissipation sont très bruitees : on va les adoucir
        ax.plot(t[:n_max],np.convolve(TOT    ,filtre,mode='same')[:n_max],'k', label=label_T)
        if not(mark_repr): ax.plot(t[:n_max],np.convolve(TOT_bis,filtre,mode='same')[:n_max],'+r',label=label_TB)
        if mark_repr : add_repr_vlines(ax,t_repr)
        ax.plot(t[:n_max],np.convolve(MOY    ,filtre,mode='same')[:n_max],'--',label=label_M)
        ax.plot(t[:n_max],np.convolve(FLC    ,filtre,mode='same')[:n_max],'-.',color='orange',label=label_F)
    else :
        # Les energies et dissipations ne sont pas trop bruitees, pas besoin de les adoucir
        ax.plot(t[:],TOT    [:],'k', label=label_T)
        if not (mark_repr) : ax.plot(t[:],TOT_bis[:],'+r',label=label_TB)
        if mark_repr : add_repr_vlines(ax,t_repr)
        ax.plot(t[:],MOY    [:],'--',label=label_M)
        ax.plot(t[:],FLC    [:],'-.',color='orange',label=label_F)
    
    
    if adim:
        ax.set_xlabel(r't / $ T_L ; T_L = \langle{K_F} \rangle /\langle{\epsilon _f} \rangle = $'+str(np.round(T_L,2)),fontsize=18)
        ax.set_ylabel(r'$\overline{k}^k / \langle{K_F} \rangle$',fontsize=18)
    else:
        ax.set_xlabel('temps (s)',fontsize=18)
        ax.set_xticks(xlabs)
        ax.set_xticklabels(np.round(xlabs,2),fontsize=18)
        ax.set_ylabel(r'$\overline{k}^k \times 10^{%s}$'%(str(COEFF)),fontsize=18)
        ax.set_yticks(YLABS)
        ax.set_yticklabels(np.round(YLABS*10**(-COEFF),3),fontsize=18)
    
    ax.tick_params(labelsize=22)
    for label in (ax.get_xticklabels() + ax.get_yticklabels()):
            label.set_fontsize(22)
    ax.legend(fontsize=22,loc='upper center', ncol=1,bbox_to_anchor=(0.5, -0.5))
    plt.tight_layout()
    fig.savefig(chemin_sauvegarde+figname+"_K_"+phase+".png")
    fig.savefig(chemin_sauvegarde+figname+"_K_"+phase+".pdf")

if draw_the_graphs:
    if "mel" in phase:
        do_the_K_plot("mel")
    if "liq" in phase:
        do_the_K_plot("liq")
    if "vap" in phase:
        do_the_K_plot("vap")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
# ENERGIE CINETIQUE DES FLUCTUATIONS SEULE (si KM et KF sont tres differents uniquement)
def do_only_KF_plot(phase):
    
    _,_,_,label_F = get_K_label(phase)
    if "mel" in phase: TOT,_,MOY,FLC,YLABS,COEFF = KT,KTB,KM,KF,KF_ylabs,KF_coeff
    elif "liq" in phase: TOT,_,MOY,FLC,YLABS,COEFF = KT_l,KTB_l,KM_l,KF_l,KF_ylabs_l,KF_coeff_l
    elif "vap" in phase: TOT,_,MOY,FLC,YLABS,COEFF = KT_v,KTB_v,KM_v,KF_v,KF_ylabs_v,KF_coeff_v
    
    if (np.mean(np.array(MOY[-int(len(MOY)/2):]))/np.mean(np.array(FLC[-int(len(FLC)/2):])) > 6.):
        fig, ax = plt.subplots(1,figsize=(9,9))
        ax.set_title(r'Energie cinetique des fluctuations : $\overline{k}^k = \frac{\overline{\chi_k u \cdot u}}{\overline{\chi_k}}$',fontsize=18)
        
        # On adoucit le signal
        if N_smooth > 0 :
            # la longeur du plus petit vecteur dicte combien de points de trace on a 
            n_max = min(len(t),len(np.convolve( TOT,filtre,mode='same'))-N_smooth)
            print("KM et KF tres differents")
            if adim :
                print(t.shape,FLC.shape)
            print(np.convolve(FLC, filtre, mode='same'))
            ax.plot(t[:n_max],np.convolve(FLC, filtre, mode='same')[:n_max],'-.',color='orange',label=label_F)
        # On n'adoucit pas le signal
        else :
            ax.plot(t[:],FLC[:],'-.',color='orange',label=label_F)
        
        # Lignes marquant les reprises
        if mark_repr : add_repr_vlines(ax,t_repr)

        if adim:
            ax.set_xlabel(r't / $ T_L ; T_L = \langle{K_F} \rangle /\langle{\epsilon _f} \rangle = $'+str(np.round(T_L,2)),fontsize=18)
            ax.set_ylabel(r'$\overline{k}^k / \langle{K_F} \rangle$',fontsize=18)
        else:    
            ax.set_xlabel('temps (s)',fontsize=18)
            ax.set_ylabel(r'$\overline{k}^k \times 10^{%s}$'%(str(COEFF)),fontsize=18)
            ax.set_xticks(xlabs)
            ax.set_xticklabels(np.round(xlabs,2),fontsize=22)
            ax.set_yticks(YLABS)
            ax.set_yticklabels(np.round(YLABS*10**(-COEFF),3),fontsize=22)
        for label in (ax.get_xticklabels() + ax.get_yticklabels()):
            label.set_fontsize(22)
        ax.legend(fontsize=22, loc='upper center', ncol=1,bbox_to_anchor=(0.5, -0.5))
        plt.tight_layout()
        fig.savefig(chemin_sauvegarde+figname+"_KF_"+phase+".png")
        fig.savefig(chemin_sauvegarde+figname+"_KF_"+phase+".pdf")

if draw_the_graphs:
    if "mel" in phase:
        do_only_KF_plot(phase)
    if "vap" in phase:
        do_only_KF_plot("vap")
    if "liq" in phase:
        do_only_KF_plot("liq")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def get_E_label(phase):
    if "mel" in phase:
        mF, coeff = mEF, E_coeff
    elif "liq" in phase:
        mF, coeff = mEF_l, E_coeff_l
    elif "vap" in phase:        
        mF, coeff = mEF_v, E_coeff_v
    
    label_ET = r'totale : $ \nu {\partial_i u_j} {\partial_j u_i}$'
    label_EM = r'moyen : $ \nu \overline{\partial_i u_j} \overline{\partial_j u_i}$'
    label_EF = (r'fluct. : $ \nu \overline{\partial_i u^p _j \partial_j u^p _i} $) $\sim$'+
                str(np.round(mF*10**(-coeff),2))+r'$\cdot 10^{%s}$'%(coeff))
    print("Dissipation Fluctuations = "+str(mF)+", (same) EF = "+str(mF))
    if diphasique:
        label_Pr = (r'$\mathcal{P} : \alpha \Delta \rho g u_r \sim $'
                    +str(np.round(mPr*10**(-coeff),2))+r'$\cdot 10^{%s}$'%(coeff))
        print("Production bulles = "+str(mPr)+", (same) Pr = "+str(mPr))

    else: label_Pr=""

    return (label_ET,label_EM,label_EF,label_Pr)
    
def get_T_E_label(phase):
    # si on l'appelle apres la declaration de TOT, MOY, ... ON PEUT DONNER EXACTMEENT LES COEFFICIENTS APPROPRIES
    if "mel" in phase:
        mF, coeff = mT_EF, T_E_coeff
    elif "liq" in phase:
        mF, coeff = mT_EF_l, T_E_coeff_l
    elif "vap" in phase:        
        mF, coeff = mT_EF_v, T_E_coeff_v

    label_ET,label_EM,_,label_Pr = get_E_label(phase)
    label_T_ET = label_ET + r'$ + \nu \partial_i u_j \partial_j u_i$'
    label_T_EM = label_EM + r'$ + \nu \overline{\partial_j u_i} \overline{\partial_j u_i}$'
    label_T_EF = (r'fluct. : $ \nu (\overline{\partial_i u^p_j \partial_j u^p_i + \partial_i u^p_j \partial_j u^p_i} ) \sim$' 
                +str(np.round(mF*10**(-coeff),2))+r'$\cdot 10^{%s}$'%(coeff))
    print("True Dissipation fluctuaitons = "+str(mF)+", (same) T_EF = "+str(mF))
    
    return(label_T_ET,label_T_EM,label_T_EF,label_Pr)

# ~~~~~~ PSEUDO DISSIPATION ET VRAIE DISSIPATION TURBULENTE \tilde{\epsilon} ~~~~~~~~~~~~
def do_the_E_plot(phase,which_one="true"):
    
    # Preparation des donnees, coefficient et labels ~~~~~~~~~~~~~~~~~~~
    if "pseudo" in which_one:
        """ pseudo-dissipaton = - nu diuj.diuj"""
        label_T,label_M,label_F,label_P = get_E_label(phase)
        if "mel" in phase: TOT,MOY,FLC,YLABS,COEFF = ET,EM,EF,E_ylabs,E_coeff
        elif "liq" in phase: TOT,MOY,FLC,YLABS,COEFF = ET_l,EM_l,EF_l,E_ylabs_l,E_coeff_l
        elif "vap" in phase: TOT,MOY,FLC,YLABS,COEFF = ET_v,EM_v,EF_v,E_ylabs_v,E_coeff_v
        if diphasique:
            if adim:
                if "mel" in phase: PRD = Pr
                elif "liq" in phase: PRD = Pr_l
                elif "vap" in phase: PRD = Pr_v
            else: 
                PRD = Pr
    else:
        """ dissipation = -2nu sij:sij = - nu diuj.diuj - nu diuj.djui"""
        label_T,label_M,label_F,label_P = get_T_E_label(phase)
        if "mel" in phase: TOT,MOY,FLC,YLABS,COEFF = T_ET,T_EM,T_EF,T_E_ylabs,T_E_coeff
        elif "liq" in phase: TOT,MOY,FLC,YLABS,COEFF = T_ET_l,T_EM_l,T_EF_l,T_E_ylabs_l,T_E_coeff_l
        elif "vap" in phase: TOT,MOY,FLC,YLABS,COEFF = T_ET_v,T_EM_v,T_EF_v,T_E_ylabs_v,T_E_coeff_v
        
        if diphasique:
            if adim:
                if "mel" in phase: PRD = Pr
                elif "liq" in phase: PRD = Pr_l
                elif "vap" in phase: PRD = Pr_v
            else: 
                PRD = Pr
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    # Creation du graphique ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    fig, ax = plt.subplots(1,figsize=(9,9))
    if "pseudo" in which_one : 
        ax.set_title(r'Pseudo - Dissipation turbulente - moyenne '+phase+r' : $\overline{\epsilon}^k = \frac{\overline{\chi_k 2 \nu \partial_j u_i \partial_j u_i }}{\overline{\chi_k}}$ ',fontsize=18)
    else : 
        ax.set_title(r'Dissipation turbulente - moyenne '+phase+
                r' : $\overline{{\epsilon}}^k = \frac{\overline{\chi_k 2 \nu s_{ij} s_{ij}}}{\overline{\chi_k}} $'+
                r'$= \frac{\overline{\chi_k \tilde{\epsilon} + \nu \partial_j u_i \partial_i u_j}}{\overline{\chi_k}} $',fontsize=18)
        
    if N_smooth > 0:
        # la longeur du plus petit vecteur dicte combien de points de trace on a 
        n_max = min(len(t),len(np.convolve( TOT,filtre,mode='same'))-N_smooth)    
        # Adoucissons le signal
        ax.plot(t[:n_max],np.convolve(TOT, filtre, mode='same')[:n_max],'k', label=label_T)
        ax.plot(t[:n_max],np.convolve(MOY, filtre, mode='same')[:n_max],'--',label=label_M)
        ax.plot(t[:n_max],np.convolve(FLC, filtre, mode='same')[:n_max],'-.',label=label_F)
        if diphasique:
            ax.plot(t[:n_max],np.convolve(PRD, filtre, mode='same')[:n_max],'-.',label=label_P)
    else:
        # N'adoucissons pas le signal
        ax.plot(t[:],TOT[:],'k', label=label_T)
        ax.plot(t[:],MOY[:],'--',label=label_M)
        ax.plot(t[:],FLC[:],'-.',label=label_F)
        if diphasique:
            ax.plot(t[:],PRD[:],'-.',label=label_P)
    # Lignes marquant les reprises
        if mark_repr : add_repr_vlines(ax,t_repr)

    if adim:
        ax.set_xlabel(r't / $ T_L ; T_L = \langle{K_F} \rangle /\langle{\epsilon _f} \rangle = $'+str(np.round(T_L,2)),fontsize=18)
        if "pseudo" in which_one:
            ax.set_ylabel(r'$\overline{\epsilon}^k / \langle{\epsilon} \rangle, \mathcal{P} / \langle{\epsilon} \rangle$',fontsize=18)
        else:
            ax.set_ylabel(r'$\overline{\tilde{\epsilon}}^k / \langle{{\epsilon}} \rangle $',fontsize=18)
    else:
        ax.set_xlabel("temps (s)",fontsize=18)
        if "pseudo" in which_one:
            ax.set_ylabel(r'$\overline{\epsilon}^k \times 10^{%s}$'%(str(COEFF)),fontsize=18)
        else:
             ax.set_ylabel(r'$\overline{\tilde{\epsilon}}^k \times 10^{%s}$'%(str(COEFF)),fontsize=18)
        ax.set_xticks(xlabs)
        ax.set_xticklabels(np.round(xlabs,2),fontsize=22)
        ax.set_yticks(YLABS)
        ax.set_yticklabels(np.round(YLABS*10**(-COEFF),3),fontsize=22)
    for label in (ax.get_xticklabels() + ax.get_yticklabels()):
        label.set_fontsize(22)
    ax.legend(fontsize=22, loc='upper center', ncol=1,bbox_to_anchor=(0.5, -0.5))
    plt.tight_layout()
    if "pseudo" in which_one:
        fig.savefig(chemin_sauvegarde+figname+"_Dissip_"+phase+".png")
        fig.savefig(chemin_sauvegarde+figname+"_Dissip_"+phase+".pdf")
    else:
        fig.savefig(chemin_sauvegarde+figname+"_TrDissip_"+phase+".png")
        fig.savefig(chemin_sauvegarde+figname+"_TrDissip_"+phase+".pdf")
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
if draw_the_graphs:
    if "mel" in phase:
        do_the_E_plot(phase,"pseudo turbu")
        do_the_E_plot(phase,"true turbu")
    if "vap" in phase:     
        do_the_E_plot("vap","pseudo turbu")
        do_the_E_plot("vap","true turbu")
    if "liq" in phase:     
        do_the_E_plot("liq","pseudo turbu")
        do_the_E_plot("liq","true turbu")
print("after do the E plot : ",t)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~ EC ET DISSIPATION SUR LE MEME GRAPHIQUE ~~~~~~~~~~~~~~~~~~~~~~
def do_EK_plot(phase,which_one):
    if adim:
        
        # Preparation des donnees, coefficient et labels ~~~~~~~~~~~~~~~~~~~
        # Ec
        _,_,_,label_KF = get_K_label(phase)
        if "mel" in phase: K_FLC = KF
        elif "liq" in phase: K_FLC = KF_l
        elif "vap" in phase: K_FLC = KF_v
        # Dissipation
        if "pseudo" in which_one:
            """ pseudo-dissipaton = - nu diuj.diuj"""
            _,_,label_EF,_ = get_E_label(phase)
            if "mel" in phase: E_FLC,YLABS,COEFF = EF,E_ylabs,E_coeff
            elif "liq" in phase: E_FLC,YLABS,COEFF = EF_l,E_ylabs_l,E_coeff_l
            elif "vap" in phase: E_FLC,YLABS,COEFF = EF_v,E_ylabs_v,E_coeff_v
        else:
            """ dissipation = -2nu sij:sij = - nu diuj.diuj - nu diuj.djui"""
            _,_,label_EF,_ = get_T_E_label(phase)
            if "mel" in phase: E_FLC,YLABS,COEFF = T_EF,T_E_ylabs,T_E_coeff
            elif "liq" in phase: E_FLC,YLABS,COEFF = T_EF_l,T_E_ylabs_l,T_E_coeff_l
            elif "vap" in phase: E_FLC,YLABS,COEFF = T_EF_v,T_E_ylabs_v,T_E_coeff_v
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        
            
        fig, ax = plt.subplots(1,figsize=(9,9))
        ax.set_title(r'Moyenne de phase '+phase,fontsize=22)
    
        ax.plot(t[:],E_FLC[:],'-',color='red',label=label_EF)
        ax.plot(t[:], K_FLC[:],'-',color='black',label=label_KF)
        
        # Lignes marquant les reprises
        if mark_repr : add_repr_vlines(ax,t_repr)

        ax.set_xlabel(r't / $ T_L ; T_L = \frac{\langle{K_F} \rangle }{\langle{\epsilon _f} \rangle} = $'+str(np.round(T_L,2)),fontsize=18)
        if "pseudo" in which_one:
            ax.set_ylabel(
                r'$\overline{k}^k / \langle{K_F} \rangle$'
                +', '+
                r'$\overline{\epsilon}^k / \langle{\epsilon} \rangle $'
                ,fontsize=22)
        else:
            ax.set_ylabel(
                    r'$\overline{k}^k / \langle{K_F} \rangle$'
                    +', '+
                    r'$\overline{\tilde{\epsilon}}^k / \langle{\tilde{\epsilon}} \rangle $'
                    ,fontsize=18)                
        for label in (ax.get_xticklabels() + ax.get_yticklabels()):
            label.set_fontsize(22)
        ax.legend(fontsize=22, loc='upper center', ncol=1,bbox_to_anchor=(0.5, -0.5))
        plt.tight_layout()
        if "pseudo" in which_one:
            fig.savefig(chemin_sauvegarde+figname+"_K_D_"+phase+".png")
            fig.savefig(chemin_sauvegarde+figname+"_K_D_"+phase+".pdf")
        else:
            fig.savefig(chemin_sauvegarde+figname+"_K_TrD_"+phase+".png")
            fig.savefig(chemin_sauvegarde+figname+"_K_TrD_"+phase+".pdf")  
            
if draw_the_graphs:
    if "mel" in phase:
        do_EK_plot(phase,"pseudo turbu")
        do_EK_plot(phase,"true turbu")
    if "vap" in phase:     
        do_EK_plot("vap","pseudo turbu")
        do_EK_plot("vap","true turbu")
    if "liq" in phase:     
        print("after do the E plot : ",t)
        do_EK_plot("liq","pseudo turbu")
        do_EK_plot("liq","true turbu")                          
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
def do_the_write_the_raws(phase):
    if "mel" in phase:
        my_E_coeff, my_T_E_coeff, my_mEF, my_mT_EF  = E_coeff, T_E_coeff, mEF, mT_EF
        my_KF_coeff, my_K_coeff, my_mKF, my_mKF_coeff = KF_coeff, K_coeff, mKF, mKF_coeff
        my_KT,my_KTB,my_KM,my_KF                    = KT,KTB,KM,KF
        my_ET,my_EM,my_EF                           = ET,EM,EF
        my_T_ET,my_T_EM,my_T_EF                     = T_ET,T_EM,T_EF
        my_phase                                    = "mel"
    if "liq" in phase:
        my_E_coeff, my_T_E_coeff, my_mEF, my_mT_EF  = E_coeff_l, T_E_coeff_l, mEF_l, mT_EF_l
        my_KF_coeff, my_K_coeff, my_mKF, my_mKF_coeff = KF_coeff_l, K_coeff_l, mKF_l, mKF_coeff_l
        my_KT,my_KTB,my_KM,my_KF                    = KT_l,KTB_l,KM_l,KF_l
        my_ET,my_EM,my_EF                           = ET_l,EM_l,EF_l
        my_T_ET,my_T_EM,my_T_EF                     = T_ET_l,T_EM_l,T_EF_l
        my_phase                                    = "liq"
    if "vap" in phase:
        my_E_coeff, my_T_E_coeff, my_mEF, my_mT_EF  = E_coeff_v, T_E_coeff_v, mEF_v, mT_EF_v
        my_KF_coeff, my_K_coeff, my_mKF, my_mKF_coeff = KF_coeff_v, K_coeff_v, mKF_v, mKF_coeff_v
        my_KT,my_KTB,my_KM,my_KF                    = KT_v,KTB_v,KM_v,KF_v
        my_ET,my_EM,my_EF                           = ET_v,EM_v,EF_v        
        my_T_ET,my_T_EM,my_T_EF                     = T_ET_v,T_EM_v,T_EF_v
        my_phase                                    = "vap"
    
    # Ecriture des information utiles pour les legendes
    if adim:
        TL__TTL = np.array((T_L, T_T_L))
        np.save(chemin_sauvegarde+"TL__TTL",TL__TTL)        
    E_coeff__T_E_coeff__mEF__mT_EF = np.array((my_E_coeff, my_T_E_coeff, my_mEF, my_mT_EF))
    np.save(chemin_sauvegarde+"Eco__TEco__mEF__mTEF"+my_phase,E_coeff__T_E_coeff__mEF__mT_EF)
    KF_coeff__K_coeff__mKF__mKF_coeff = np.array((my_KF_coeff, my_K_coeff, my_mKF, my_mKF_coeff)) 
    np.save(chemin_sauvegarde+"KFco__Kco__mKF__mKFco"+my_phase,KF_coeff__K_coeff__mKF__mKF_coeff)
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    # Ecriture des champs interessants
    # ENERGIE CINETIQUE
    temps_KT_KTB_KM_KF = np.array((t[:],my_KT[:],my_KTB[:],my_KM[:],my_KF[:]))
    write_me_in_raw(chemin_sauvegarde,temps_KT_KTB_KM_KF,adim,"temps_KT_KTB_KM_KF_"+phase)    
    # PSEUDO DISSIPATION TURBULENTE \tilde{\epsilon}
    temps_ET_EM_EF = np.array((t[:],my_ET[:],my_EM[:],my_EF[:]))
    write_me_in_raw(chemin_sauvegarde,temps_ET_EM_EF,adim,"temps_ET_EM_EF_"+phase)
    # PRODUCITON PAR LES BULLES
    if diphasique:
        temps_Pr = np.array((t[:],Pr[:]))
        write_me_in_raw(chemin_sauvegarde,temps_Pr,adim,"temps_Pr_"+phase)
    # VRAIE DISSIPATION TURBULENTE \epsilon 
    temps_T_ET_T_EM_T_EF = np.array((t_t[:],my_T_ET[:],my_T_EM[:],my_T_EF[:]))
    write_me_in_raw(chemin_sauvegarde,temps_T_ET_T_EM_T_EF,adim,"temps_T_ET_T_EM_T_EF_"+phase)    
    # "VRAIES" PRODUCTION PAR LES BULLES (dans le cas adim
    if diphasique and adim:
        temps_T_Pr = np.array((t_t[:],T_Pr[:]))
        write_me_in_raw(chemin_sauvegarde,temps_T_Pr,adim,"temps_T_Pr_"+phase)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
if write_the_raws: 
    if "mel" in phase:
        do_the_write_the_raws("mel")
    if "liq" in phase:
        do_the_write_the_raws("liq")
    if "vap" in phase:
        do_the_write_the_raws("vap")

# Grandeurs affichees sur les plots 
def do_the_final_print(phase):
    
    if "mel" in phase:
        my_E_coeff, my_T_E_coeff, my_mEF, my_mT_EF  = E_coeff, T_E_coeff, mEF, mT_EF
        my_KF_coeff, my_K_coeff, my_mKF, my_mKF_coeff = KF_coeff, K_coeff, mKF, mKF_coeff
        my_KT,my_KTB,my_KM,my_KF                    = KT,KTB,KM,KF
        my_ET,my_EM,my_EF                           = ET,EM,EF
        my_T_ET,my_T_EM,my_T_EF                     = T_ET,T_EM,T_EF
    if "liq" in phase:
        my_E_coeff, my_T_E_coeff, my_mEF, my_mT_EF  = E_coeff_l, T_E_coeff_l, mEF_l, mT_EF_l
        my_KF_coeff, my_K_coeff, my_mKF, my_mKF_coeff = KF_coeff_l, K_coeff_l, mKF_l, mKF_coeff_l
        my_KT,my_KTB,my_KM,my_KF                    = KT_l,KTB_l,KM_l,KF_l
        my_ET,my_EM,my_EF                           = ET_l,EM_l,EF_l
        my_T_ET,my_T_EM,my_T_EF                     = T_ET_l,T_EM_l,T_EF_l
    if "vap" in phase:
        my_E_coeff, my_T_E_coeff, my_mEF, my_mT_EF  = E_coeff_v, T_E_coeff_v, mEF_v, mT_EF_v
        my_KF_coeff, my_K_coeff, my_mKF, my_mKF_coeff = KF_coeff_v, K_coeff_v, mKF_v, mKF_coeff_v
        my_KT,my_KTB,my_KM,my_KF                    = KT_v,KTB_v,KM_v,KF_v
        my_ET,my_EM,my_EF                           = ET_v,EM_v,EF_v        
        my_T_ET,my_T_EM,my_T_EF                     = T_ET_v,T_EM_v,T_EF_v
        
    print("####"+phase+", adim : "+str(adim)+"####")
    print("Les moyennes sont obtenues sur : t = ["+str(t[-int(len(t)/2)])+"; "+str(t[-1])+"]")
    print("Moyenne de l'énergie cinétique des fluctuations :")
    print(r'$\langle{K_F} \rangle = $ '+str(np.round(my_mKF*10**(-my_mKF_coeff),5))+r'$\cdot 10^{%s}$'%(my_mKF_coeff))    
    print(" ==> <K_F> = ",(np.round(my_mKF*10**(-my_mKF_coeff),5))*10**(my_mKF_coeff))
    print("Estimation de la puissance injectee par les bulles :")
    if diphasique:
        print(r'$\mathcal{P} : \alpha \Delta \rho g u_r \sim $'+str(np.round(mPr*10**(-my_E_coeff),2))+r'$\cdot 10^{%s}$'%(my_E_coeff))
        print(" ==> P_inj = ",(np.round(mPr*10**(-my_E_coeff),5))*10**(my_E_coeff))
    print("Moyenne de la pseudo-dissipation, calculée sur "+str(t[-int(len(t)/2)])+"< t < "+str(t[-1]))
    print(r'$\langle{E_F} \rangle = $ '+str(np.round(my_mEF*10**(-my_E_coeff),5))+r'$\cdot 10^{%s}$'%(my_E_coeff))  
    print(" ==> <E_F> = ",(np.round(my_mEF*10**(-my_E_coeff),5))*10**(my_E_coeff))  
    print("Moyenne de la vraie dissipation, calculée sur "+str(t[-int(len(t)/2)])+"< t < "+str(t[-1]))
    print(r'$\langle{K_F} \rangle = $ '+str(np.round(my_mT_EF*10**(-my_T_E_coeff),5))+r'$\cdot 10^{%s}$'%(my_T_E_coeff))   
    print(" ==> <TE_F> = ",(np.round(my_mEF*10**(-my_T_E_coeff),5))*10**(my_T_E_coeff)) 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if "mel" in phase:
    do_the_final_print("mel")
if "liq" in phase:
    do_the_final_print("liq")
if "vap" in phase:
    do_the_final_print("vap")
    
message = "\
####\n\
# Fin du script energy_for_fusion.py \n\
######################################################################\n"
print(message)
