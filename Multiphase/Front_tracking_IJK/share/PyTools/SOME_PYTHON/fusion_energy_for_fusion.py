# -*- coding: utf-8 -*-
import numpy as np
import math as mt
import matplotlib.pyplot as plt
import glob
import sys
import os
import DNSTools3 as dtool


"""
    Voir les histoires de fontsize de label : je donne 2 fois les commandes pour avoir des grosses polices et ca me donne des petites polices whla je craque    
"""
############################################################## Mode d'emploi #######################################
# mode d'emploi : python fusion_energy_for_fusion.py liq_vap_mel nom_des_figures liste_fichier.data liste_chemin_vers_diphasique_moyenne_spatiale**.txt 
#                          0                              1           2                 3                                   4             
#                        fichier_out_ou_log_pour_Reb...   N_smooth  mode   [t_deb,t_fin]  liste_legendes
#                            5                                6      7         8              9
#
# mode d'emploi2 : python fusion_energy_for_fusion.py liq_vap_mel  nom_des_figures  liste_chemins liste_legendes
#                                                           1*               2*            3*             4*
# /!\ 2 : Si on met adim ou ADIM dans le nom, alors on trace les evolution de K/<K>, E/<E> en fonction de t/<T>
#         avec <T> = <K>/<E>
# 3 et 4 : listes de la forme [..., ..., ...]
# /!\ 4 : chaque chemin se temrine par un slash "/"
# 5 : si la chaine contient "do_the_accplot", le script accplot_vAB.py sera execute avec les arguments standards
#     /!\ pour une simulation monophasique, accplot_vAB.py execute use_check_etapes_et_termes, enfin cette 
#         execution en chaine est a finir de coder
#     /!\ pour une simulation monphasique, le fichier de sortie de use_check_etapes_et_termes.py n'est pas utilise
# 7 : draw -> dessine les graphiques directement
#     write -> ecrits les fichiers raw
#     both -> desisne les graphiques et ecrit les fichiers raw
# 9 : Ce qui va apparaitre dans la legend des plots
# 3*: liste des chemins ou sont ecrits les .raw
# 4*: liste des legendes qui vont apparaitre sur les plots
############################################################## fin Mode d'emploi #######################################

####################################################################
# 0 - Preparation des meta-donnees / donnees utilisateur, chemins, ...
####################################################################
# 0.0 - Donnees utilisateur
chemin_python = '$project_directory/share/PyTools/SOME_PYTHON/'
print(sys.argv)

if len(sys.argv)==10:
    # Ecriture des energies et dissipations au format raw
    phase =sys.argv[1]
    nom_figures = sys.argv[2]
    liste_jdd = sys.argv[3].replace("[","").replace("]","").split(",")
    liste_chemin_txt = sys.argv[4].replace("[","").replace("]","").split(",")
    fichier_acc = sys.argv[5]
    N_smooth = int(sys.argv[6])
    mode = sys.argv[7]
    temps = sys.argv[8]
    liste_chemins = []
    for chemin_txt in liste_chemin_txt:
        liste_chemins.append(chemin_txt.replace("TXT/",""))
    liste_legendes = sys.argv[9].replace("[","").replace("]","").split(",")
    
    for index, chemin_txt in enumerate(liste_chemin_txt):
        if "do_the_accplot" in fichier_acc:
            arguments_accplot = liste_chemins[index]+'fig_acc 0. '+liste_chemins[index]+'NEXT/'+liste_jdd[index].replace("data","sauv")+" "+liste_chemins[index]+liste_jdd[index]
            os.system('python '+chemin_python+'accplot_vAB.py '+arguments_accplot)
        arguments = phase+" "+liste_chemins[index]+nom_figures+" "+liste_chemins[index]+liste_jdd[index]+" "+liste_chemin_txt[index]+" "+fichier_acc+" "+str(N_smooth)+" "+mode+" "+temps
        print("arguments" + arguments)
        os.system('python '+chemin_python+'energy_for_fusion.py '+arguments)
elif len(sys.argv)==5:
    # L'ecriture des energies et dissipation en raw est deja faite
    phase = sys.argv[1]
    nom_figures = sys.argv[2]
    liste_chemins = sys.argv[3].replace("[","").replace("]","").split(",")
    liste_legendes = sys.argv[4].replace("[","").replace("]","").split(",")
adim = ('adim' in nom_figures)
# legendes associees a la vraie dissipation
liste_T_legendes = liste_legendes.copy()
# 0.1 Donnees pour les graphiques
# a remplacer certainement par un 'liste_linewidth' pour les daltoniens
liste_couleurs = ['r','g','b','c','m','y','k']
# ENERGIES CINETIQUES
fig_k, ax_k = plt.subplots(1,figsize=(9,7))
ax_k.set_title('k - moyenne de phase '+phase+r' : $\overline{k}^k = \frac{\overline{\chi_k \frac{1}{2} u \cdot u}}{\overline{\chi_k}}$',fontsize=18)
ax_k0, ax_kN = ax_k.twinx(), ax_k.twinx()
# ENERGIE CINETIQUE FLUCTUANTE UNIQUEMENT
fig_kf, ax_kf = plt.subplots(1,figsize=(9,7))
ax_kf.set_title(r"k_f : $\overline{k}^k = \frac{\overline{\chi_k \frac{1}{2} u' \cdot u'}}{\overline{\chi_k}}$",fontsize=18)
ax_kf0, ax_kfN = ax_kf.twinx(), ax_kf.twinx()
# PSEUDO DISSIPATION
fig_pe, ax_pe = plt.subplots(1,figsize=(9,7))
ax_pe.set_title(r'$\tilde{\epsilon}$ - moyenne de phase '+phase+r' : $\frac{\overline{\chi_k \nu_k \partial_j u_i \partial_j u_i }}{\overline{\chi_k}}$ ',fontsize=18)
ax_pe0, ax_peN = ax_pe.twinx(), ax_pe.twinx()
# ENERGIE CINETIQUE ET PSEUDO DISSIPATION
if adim :
    fig_kpe, ax_kpe = plt.subplots(1,figsize=(9,7))
    ax_kpe.set_title(r'Moyenne de phase '+phase,fontsize=22)
    ax_kpe0, ax_kpeN = ax_kpe.twinx(), ax_kpe.twinx()
# TRUE DISSIPATION
fig_te, ax_te = plt.subplots(1,figsize=(9,7))
ax_te.set_title(r'$\epsilon$ - moyenne '+phase+
                 r' : $\frac{\overline{\chi_k 2 \nu_k s_{ij} s_{ij}}}{\overline{\chi_k}} $'+
                 r'$= \frac{\overline{\chi_k \tilde{\epsilon} + \nu_k \partial_j u_i \partial_i u_j}}{\overline{\chi_k}} $',fontsize=18)
ax_te0, ax_teN = ax_te.twinx(), ax_te.twinx()
# ENERGIE CINETIQUE ET TRUE DISSIPATION
if adim :
    fig_kte, ax_kte = plt.subplots(1,figsize=(9,7))
    ax_kte.set_title(r'Moyenne de phase '+phase,fontsize=22)
    ax_kte0, ax_kteN = ax_kte.twinx(), ax_kte.twinx()
# GRANDEURS MOYENNES ET INTEGRALES
fig_mi, ax_mi = plt.subplots(1,figsize=(9,7))
ax_mi.set_title(r'$\langle K_F \rangle$, $\langle \epsilon \rangle$ et $T_L$ '+phase,fontsize=22)
ax_mi0, ax_miN = ax_mi.twinx(), ax_mi.twinx()

for nc,chemin in enumerate(liste_chemins):
    diphasique=False
    Pr = np.array(2)
    ####################################################################
    # 1 - Preparation des donnees
    ####################################################################
    # 1.0 - Chargement des dimensions des champs
    the_shapes_of_temps_KT_KTB_KM_KF   = np.load(chemin+"the_shapes_of_temps_KT_KTB_KM_KF.npy")
    the_shapes_of_temps_ET_EM_EF       = np.load(chemin+"the_shapes_of_temps_ET_EM_EF.npy")  
    if adim:
        the_shapes_of_temps_EF_KF          = np.load(chemin+"the_shapes_of_temps_EF_KF.npy")
    the_shapes_of_temps_T_ET_T_EM_T_EF = np.load(chemin+"the_shapes_of_temps_T_ET_T_EM_T_EF.npy")
    if adim:
        the_shapes_of_temps_T_EF_KF        = np.load(chemin+"the_shapes_of_temps_T_EF_KF.npy")
    # 1.1 - Chargement et mise en forme des champs
    temps_KT_KTB_KM_KF    = np.fromfile(chemin+"temps_KT_KTB_KM_KF.raw").reshape(the_shapes_of_temps_KT_KTB_KM_KF)
    temps_ET_EM_EF        = np.fromfile(chemin+"temps_ET_EM_EF.raw").reshape(the_shapes_of_temps_ET_EM_EF)
    if adim:
        temps_EF_KF           = np.fromfile(chemin+"temps_EF_KF.raw").reshape(the_shapes_of_temps_EF_KF)
    temps_T_ET_T_EM_T_EF  = np.fromfile(chemin+"temps_T_ET_T_EM_T_EF.raw").reshape(the_shapes_of_temps_T_ET_T_EM_T_EF)
    if adim:
        temps_T_EF_KF         = np.fromfile(chemin+"temps_T_EF_KF.raw").reshape(the_shapes_of_temps_T_EF_KF)
    # 1.2 - Calcul des coefficients, moyennes, etc.
    if adim:
        T_L, T_T_L = np.load(chemin+"TL__TTL.npy")
    E_coeff,T_E_coeff,mEF,mT_EF = np.load(chemin+"Eco__TEco__mEF__mTEF.npy")
    KF_coeff, coeff, mKF, mKF_coeff = np.load(chemin+"KFco__co__mKF__mKFco.npy")
    
    
    ####################################################################
    # 2 - Traces des evolutions
    ####################################################################
    # 2.1 - ENERGIES CINETIQUES
    # On trace le grapique d evoluiton des donnees.
    ax_k.plot(temps_KT_KTB_KM_KF[0,:],temps_KT_KTB_KM_KF[1,:],color=liste_couleurs[nc],linestyle='-')
    ax_k.plot(temps_KT_KTB_KM_KF[0,:],temps_KT_KTB_KM_KF[3,:],color=liste_couleurs[nc],linestyle='--')
    ax_k.plot(temps_KT_KTB_KM_KF[0,:],temps_KT_KTB_KM_KF[4,:],color=liste_couleurs[nc],linestyle=':')    
    if nc==0:
        # On ecrit la legende [total,moyen,fluctuant] uniquement pour premier chemin
        ax_k0.plot(np.NaN, np.NaN,'k',label=r'$\frac{1}{2} (U \cdot U)$')
        ax_k0.plot(np.NaN, np.NaN,'--k',label=r'$ \frac{1}{2} (\overline{u} \cdot \overline{u}) $')
        ax_k0.plot(np.NaN, np.NaN,':k',color='black',label=r"$\frac{1}{2} ( u' \cdot u')$")
        ax_k0.get_yaxis().set_visible(False)
        # Personalise les titres des axes
        if adim:
            ax_k.set_xlabel(r't / $ T_L ; T_L = \langle{K_F} \rangle /\langle{\epsilon _f} \rangle $',fontsize=18)
            ax_k.set_ylabel(r'$\overline{k}^k / \langle{K_F} \rangle$',fontsize=18)
        else:
            ax_k.set_xlabel('temps (s)',fontsize=18)
            ax_k.set_ylabel(r'$\overline{k}^k$',fontsize=18)    
    # On ecrit la legende [numero_du_run] pour tous les chemins
    if adim:
        liste_legendes[nc] += r' : $T_L = $'+str(np.round(T_L,2))+r'; $\langle{K_F} \rangle = $ '+str(np.round(mKF*10**(-mKF_coeff),3))+r'$\cdot 10^{%s}$'%(mKF_coeff)
        liste_T_legendes[nc] += r' : $T_L = $'+str(np.round(T_T_L,2))+r'; $\langle{K_F} \rangle = $ '+str(np.round(mKF*10**(-mKF_coeff),3))+r'$\cdot 10^{%s}$'%(mKF_coeff)
    else:
        liste_legendes[nc] += r' : $\langle{K_F} \rangle = $ '+str(np.round(mKF*10**(-mKF_coeff),3))+r'$\cdot 10^{%s}$'%(mKF_coeff)
        liste_T_legendes[nc] += r' : $\langle{K_F} \rangle = $ '+str(np.round(mKF*10**(-mKF_coeff),3))+r'$\cdot 10^{%s}$'%(mKF_coeff)
    ax_kN.plot(np.NaN, np.NaN,color=liste_couleurs[nc],label=liste_legendes[nc])
    ax_kN.get_yaxis().set_visible(False)  
    # Pour que les tites des axes soient avec la taille de police choisie      
    ax_k.tick_params(labelsize=22)
    for label in (ax_k.get_xticklabels() + ax_k.get_yticklabels()):
            label.set_fontsize(22)

    # 2.2 - ENERGIE CINETIQUE DES FLUCTUATIONS    
    if (np.mean(np.array(temps_KT_KTB_KM_KF[3,:][-int(len(temps_KT_KTB_KM_KF[3,:])/2):]))/np.mean(np.array(temps_KT_KTB_KM_KF[4,:][-int(len(temps_KT_KTB_KM_KF[4,:])/2):])) > 6.):
        # On trace le grapique d evoluiton des donnees.
        ax_kf.plot(temps_EF_KF[0,:],temps_EF_KF[1,:],color=liste_couleurs[nc],linestyle=':')   
        if nc==0:
            # On ecrit la legende [total,moyen,fluctuant] uniquement pour premier chemin
            ax_kf0.plot(np.NaN, np.NaN,':k',color='black',label=r"$\frac{1}{2} ( u' \cdot u')$")
            ax_kf0.get_yaxis().set_visible(False)
            # Personalise les titres des axes
            if adim:
                ax_kf.set_xlabel(r'$t /  T_L ; T_L = \langle{K_F} \rangle /\langle{\epsilon _f} \rangle $',fontsize=18)
                ax_kf.set_ylabel(r'$\overline{k}^k / \langle{K_F} \rangle$',fontsize=18)
            else:
                ax_kf.set_xlabel('temps (s)',fontsize=18)
                ax_kf.set_ylabel(r'$\overline{k}^k$',fontsize=18)    
        # On ecrit la legende [numero_du_run] pour tous les chemins, uniquement pour ce plot. Les autres re-piocherons
        liste_legendes[nc] += r' : $T_L = $'+str(np.round(T_L,2))+r'; $\langle{K_F} \rangle = $ '+str(np.round(mKF*10**(-mKF_coeff),3))+r'$\cdot 10^{%s}$'%(mKF_coeff)
        ax_kfN.plot(np.NaN, np.NaN,color=liste_couleurs[nc],label=liste_legendes[nc])
        ax_kfN.get_yaxis().set_visible(False)  
        # Pour que les tites des axes soient avec la taille de police choisie      
        ax_kf.tick_params(labelsize=22)
        for label in (ax_kf.get_xticklabels() + ax_kf.get_yticklabels()):
            label.set_fontsize(22)
            
    # 2.3 PSEUDO DISSIPATION TURBULENTE \tilde{\epsilon}
    # On trace le grapique d evoluiton des donnees.
    ax_pe.plot(temps_ET_EM_EF[0,:],temps_ET_EM_EF[1,:],color=liste_couleurs[nc],linewidth=0.5,linestyle='-')
    ax_pe.plot(temps_ET_EM_EF[0,:],temps_ET_EM_EF[2,:],color=liste_couleurs[nc],linewidth=0.5,linestyle='--')
    ax_pe.plot(temps_ET_EM_EF[0,:],temps_ET_EM_EF[3,:],color=liste_couleurs[nc],linewidth=0.5,linestyle=':')    
    if nc==0:
        # On ecrit la legende [total,moyen,fluctuant] uniquement pour premier chemin
        ax_pe0.plot(np.NaN, np.NaN,'k',linewidth=0.5,label=r'$\nu (\partial_j U_i \cdot \partial_j U_i )$')
        ax_pe0.plot(np.NaN, np.NaN,'--k',linewidth=0.5,label=r'$\nu (\overline{\partial_j U_i} \cdot \overline{\partial_j U_i})$')
        ax_pe0.plot(np.NaN, np.NaN,':k',linewidth=0.5,color='black',label=r"$\nu (\partial_j U_i ' \cdot \partial_j U_i ')$")
        ax_pe0.get_yaxis().set_visible(False)
        # Personalise les titres des axes
        if adim:
            ax_pe.set_xlabel(r't / $ T_L ; T_L = \langle{K_F} \rangle /\langle{\epsilon _f} \rangle $',fontsize=18)
            ax_pe.set_ylabel(r'$\overline{\epsilon}^k / \langle{\epsilon} \rangle$',fontsize=18)
        else:
            ax_pe.set_xlabel('temps (s)',fontsize=18)
            ax_pe.set_ylabel(r'$\overline{\epsilon}^k$',fontsize=18)    
    # On ecrit la legende [numero_du_run] pour tous les chemins
    ax_peN.plot(np.NaN, np.NaN,color=liste_couleurs[nc],linewidth=0.5,label=liste_legendes[nc])
    ax_peN.get_yaxis().set_visible(False)  
    # Pour que les tites des axes soient avec la taille de police choisie      
    ax_pe.tick_params(labelsize=22)
    for label in (ax_pe.get_xticklabels() + ax_pe.get_yticklabels()):
            label.set_fontsize(22)
            
    # 2.4 ENERGIE CINETIQUE ET PSEUDO DISSIPATION TURBULENTE \tilde{\epsilon}
    # On trace le grapique d evolution des donnees.
    if adim:
        ax_kpe.plot(temps_EF_KF[0,:],temps_EF_KF[2,:],color=liste_couleurs[nc],linestyle='-')
        ax_kpe.plot(temps_EF_KF[0,:],temps_EF_KF[1,:],color=liste_couleurs[nc],linestyle='--',linewidth=0.5)
        if nc==0:
            # On ecrit la legende [total,moyen,fluctuant] uniquement pour premier chemin
            ax_kpe0.plot(np.NaN, np.NaN,'k',label=r"$\frac{1}{2} ( u' _i \cdot u' _i)$")
            ax_kpe0.plot(np.NaN, np.NaN,'k',linewidth=0.5,label=r"$\nu (\partial_j u' _i \cdot \partial_j u' _i )$")
            ax_kpe0.get_yaxis().set_visible(False)
            # Personalise les titres des axes
            if adim:
                ax_kpe.set_xlabel(r't / $ T_L ; T_L = \langle{K_F} \rangle /\langle{\epsilon _f} \rangle $',fontsize=18)
                ax_kpe.set_ylabel(r'$\overline{\epsilon}^k / \langle{\epsilon} \rangle, \overline{k}^k / \langle{K_F} \rangle$',fontsize=18)
            else:
                ax_kpe.set_xlabel('temps (s)',fontsize=18)
                ax_kpe.set_ylabel(r'$\overline{\epsilon}^k, overline{k}^k$',fontsize=18)    
        # On ecrit la legende [numero_du_run] pour tous les chemins
        ax_kpeN.plot(np.NaN, np.NaN,color=liste_couleurs[nc],label=liste_legendes[nc])
        ax_kpeN.get_yaxis().set_visible(False)  
        # Pour que les tites des axes soient avec la taille de police choisie      
        ax_kpe.tick_params(labelsize=22)
        for label in (ax_kpe.get_xticklabels() + ax_kpe.get_yticklabels()):
                label.set_fontsize(22)
           
    # 2.5 VRAIE DISSIPATION TURBULENTE {\epsilon}
    # On trace le grapique d evoluiton des donnees.
    ax_te.plot(temps_T_ET_T_EM_T_EF[0,:],temps_T_ET_T_EM_T_EF[1,:],color=liste_couleurs[nc],linewidth=0.5,linestyle='-')
    ax_te.plot(temps_T_ET_T_EM_T_EF[0,:],temps_T_ET_T_EM_T_EF[2,:],color=liste_couleurs[nc],linewidth=0.5,linestyle='--')
    ax_te.plot(temps_T_ET_T_EM_T_EF[0,:],temps_T_ET_T_EM_T_EF[3,:],color=liste_couleurs[nc],linewidth=0.5,linestyle=':')    
    if nc==0:
        # On ecrit la legende [total,moyen,fluctuant] uniquement pour premier chemin
        ax_te0.plot(np.NaN, np.NaN,'k',linewidth=0.5,label=r'$2 \nu ( S_{ij} \cdot S_{ij})$')
        ax_te0.plot(np.NaN, np.NaN,'--k',linewidth=0.5,label=r'$2\nu (\overline{s_{ij}} \cdot \overline{s_{ij}})$')
        ax_te0.plot(np.NaN, np.NaN,':k',linewidth=0.5,color='black',label=r"$2 \nu ( s' _{ij} \cdot s' _{ij})$")
        ax_te0.get_yaxis().set_visible(False)
        # Personalise les titres des axes
        if adim:
            ax_te.set_xlabel(r't / $ T_L ; T_L = \langle{K_F} \rangle /\langle{\epsilon _f} \rangle $',fontsize=18)
            ax_te.set_ylabel(r'$\overline{\epsilon}^k / \langle{\epsilon} \rangle$',fontsize=18)
        else:
            ax_te.set_xlabel('temps (s)',fontsize=18)
            ax_te.set_ylabel(r'$\overline{\epsilon}^k$',fontsize=18)    
    # On ecrit la legende [numero_du_run] pour tous les chemins
    ax_teN.plot(np.NaN, np.NaN,color=liste_couleurs[nc],linewidth=0.5,label=liste_T_legendes[nc])
    ax_teN.get_yaxis().set_visible(False)  
    # Pour que les tites des axes soient avec la taille de police choisie      
    ax_te.tick_params(labelsize=22)
    for label in (ax_te.get_xticklabels() + ax_te.get_yticklabels()):
            label.set_fontsize(22)            
            
    # 2.6 ENERGIE CINETIQUE ET VRAIE DISSIPATION TURBULENTE {\epsilon}
    if adim:
        # On trace le grapique d evoluiton des donnees.
        ax_kte.plot(temps_T_EF_KF[0,:],temps_T_EF_KF[2,:],color=liste_couleurs[nc],linestyle='-')
        ax_kte.plot(temps_T_EF_KF[0,:],temps_T_EF_KF[1,:],color=liste_couleurs[nc],linestyle='--',linewidth=0.5)
        if nc==0:
            # On ecrit la legende [total,moyen,fluctuant] uniquement pour premier chemin
            ax_kte0.plot(np.NaN, np.NaN,'k',label=r"$\frac{1}{2} (u' _i \cdot u' _i U)$")
            ax_kte0.plot(np.NaN, np.NaN,'--k',linewidth=0.5,label=r"$ 2 \nu (s' _{ij} \cdot s' _{ij})$")
            ax_kte0.get_yaxis().set_visible(False)
            # Personalise les titres des axes
            if adim:
                ax_kte.set_xlabel(r't / $ T_L ; T_L = \langle{K_F} \rangle /\langle{\epsilon _f} \rangle $',fontsize=18)
                ax_kte.set_ylabel(
                        r'$\overline{k}^k / \langle{K_F} \rangle$'
                        +', '+
                        r'$\overline{\epsilon}^k / \langle{\epsilon} \rangle $'
                        ,fontsize=18)
            else:
                ax_kte.set_xlabel(r'temps $(s)$',fontsize=18)
                ax_kte.set_ylabel(
                        r'$\overline{k}^k $'
                        +', '+
                        r'$\overline{\epsilon}^k $'
                        ,fontsize=18)   
        # On ecrit la legende [numero_du_run] pour tous les chemins
        # liste_legendes[nc] += r' : $T_L = $'+str(np.round(T_L,2))+r'; $\langle{K_F} \rangle = $ '+str(np.round(mKF*10**(-mKF_coeff),3))+r'$\cdot 10^{%s}$'%(mKF_coeff)
        ax_kteN.plot(np.NaN, np.NaN,color=liste_couleurs[nc],label=liste_T_legendes[nc])
        ax_kteN.get_yaxis().set_visible(False)  
        # Pour que les tites des axes soient avec la taille de police choisie      
        ax_kte.tick_params(labelsize=22)
        for label in (ax_kte.get_xticklabels() + ax_kte.get_yticklabels()):
                label.set_fontsize(22)                                        
        
    # 2.7 EVOLUTION DE <K_F>, <\epsilon>, et T_L EN FONCTION DES CAS
    # On trace le grapique d evoluiton des donnees.
    ax_mi.plot(nc,mKF,color=liste_couleurs[nc],linestyle='-')
    ax_mi.plot(nc,mT_EF,color=liste_couleurs[nc],linestyle='--',linewidth=0.5)
    if nc==0:
        # On ecrit la legende [total,moyen,fluctuant] uniquement pour premier chemin
        ax_mi0.plot(np.NaN, np.NaN,'k',label=r"$\frac{1}{2} \langle (u' _i \cdot u' _i U) \rangle $")
        ax_mi0.plot(np.NaN, np.NaN,'--k',linewidth=0.5,label=r"$ 2 \nu \langle (s' _{ij} \cdot s' _{ij}) \rangle $")
        ax_mi0.get_yaxis().set_visible(False)
        # Personalise les titres des axes
        if adim:
            ax_mi.set_xlabel(r'Simulation',fontsize=18)
            ax_mi.set_ylabel(
                    r'$\langle{K_F} \rangle$'
                    +', '+
                    r'$ \langle{\epsilon} \rangle $'
                    ,fontsize=18)
        else:
            ax_mi.set_xlabel(r'simulation',fontsize=18)
            ax_mi.set_ylabel(
                    r'$\langle{K_F} \rangle$'
                    +', '+
                    r'$\langle{K_F} \rangle$'
                    ,fontsize=18)   
    # On ecrit la legende [numero_du_run] pour tous les chemins
    # liste_legendes[nc] += r' : $T_L = $'+str(np.round(T_L,2))+r'; $\langle{K_F} \rangle = $ '+str(np.round(mKF*10**(-mKF_coeff),3))+r'$\cdot 10^{%s}$'%(mKF_coeff)
    ax_miN.plot(np.NaN, np.NaN,color=liste_couleurs[nc],label=liste_T_legendes[nc])
    ax_miN.get_yaxis().set_visible(False)  
    # Pour que les tites des axes soient avec la taille de police choisie      
    ax_mi.tick_params(labelsize=22)
    for label in (ax_mi.get_xticklabels() + ax_mi.get_yticklabels()):
            label.set_fontsize(22)                                        
                

ax_kN.legend(loc='lower right',fontsize=22)
ax_k0.legend(fontsize=22)
    ax_k0.tick_params(labelsize=22)
    for label in (ax_k0.get_xticklabels() + ax_k0.get_yticklabels()):
            label.set_fontsize(22)
fig_k.tight_layout()
fig_k.savefig(nom_figures+"_K.png")
fig_k.savefig(nom_figures+"_K.pdf")

if (np.mean(np.array(temps_KT_KTB_KM_KF[3,:][-int(len(temps_KT_KTB_KM_KF[3,:])/2):]))/np.mean(np.array(temps_KT_KTB_KM_KF[4,:][-int(len(temps_KT_KTB_KM_KF[4,:])/2):])) > 6.):
    ax_kfN.legend(loc='lower right',fontsize=22)
    ax_kf0.legend(fontsize=22)
    ax_kfN.tick_params(labelsize=22)
    for label in (ax_kfN.get_xticklabels() + ax_kfN.get_yticklabels()):
            label.set_fontsize(22)
    fig_kf.tight_layout()
    fig_kf.savefig(nom_figures+"_KF.png")
    fig_kf.savefig(nom_figures+"_KF.pdf")

ax_peN.legend(loc='lower right',fontsize=22)
ax_pe0.legend(fontsize=22)
    ax_peN.tick_params(labelsize=22)
    for label in (ax_peN.get_xticklabels() + ax_peN.get_yticklabels()):
            label.set_fontsize(22)
fig_pe.tight_layout()
fig_pe.savefig(nom_figures+"_Dissip.png")
fig_pe.savefig(nom_figures+"_Dissip.pdf")

if adim:
    ax_kpeN.legend(loc='lower right',fontsize=22)
    ax_kpe0.legend(fontsize=22)
    ax_kpeN.tick_params(labelsize=22)
    for label in (ax_kpeN.get_xticklabels() + ax_kpeN.get_yticklabels()):
            label.set_fontsize(22)
    fig_kpe.tight_layout()
    fig_kpe.savefig(nom_figures+"_K_D.png")
    fig_kpe.savefig(nom_figures+"_K_D.pdf")

ax_teN.legend(loc='lower right',fontsize=22)
ax_te0.legend(fontsize=22)
    ax_teN.tick_params(labelsize=22)
    for label in (ax_teN.get_xticklabels() + ax_teN.get_yticklabels()):
            label.set_fontsize(22)
fig_te.tight_layout()
fig_te.savefig(nom_figures+"_TrDissip.png")
fig_te.savefig(nom_figures+"_TrDissip.pdf")

if adim:
    ax_kteN.legend(loc='lower right',fontsize=22)
    ax_kte0.legend(fontsize=22)
    ax_kteN.tick_params(labelsize=22)
    for label in (ax_kteN.get_xticklabels() + ax_kteN.get_yticklabels()):
            label.set_fontsize(22)    
    fig_kte.tight_layout()
    fig_kte.savefig(nom_figures+"_K_TrD.png")
    fig_kte.savefig(nom_figures+"_K_TrD.pdf")
        
ax_miN.legend(loc='lower right',fontsize=22)
ax_mi0.legend(fontsize=22)
ax_miN.tick_params(labelsize=22)
for label in (ax_miN.get_xticklabels() + ax_miN.get_yticklabels()):
        label.set_fontsize(22)
fig_mi.tight_layout()
fig_mi.savefig(nom_figures+"_K_e_T_L.png")
fig_mi.savefig(nom_figures+"_K_e_T_L.pdf")
        
plt.plot()          
