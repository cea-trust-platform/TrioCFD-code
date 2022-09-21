import sys
import datetime
import matplotlib.pyplot as plt
import numpy as np
import DNSTools3 as dtool

"""
L'objectif de ce script est de :
 V- Afficher l'evolution temporelle de la correction orthogonale a g (le patch de qdm)   : jdd_acceleration.out
 V- Afficher l'evolution temporelle de la correction colineaire a g (la source de qdm)   : jdd_acceleration.out
 X- Afficher la PDF des deux corrections precedentes, donnant la moyenne, et la variance : jdd_acceleration.out
 V- Afficher l'evolution temporelle des termes du bilan de qdm                           : jdd_check_etapes_et_termes.out
 X- Afficher la PDF  des termes du bilan de qdm                                          : jdd_check_etapes_et_termes.out
 X- Comparer les corrections avec les termes pertinents.
"""

"""
temps     :      n+1      n       /                 1           1                    1                                       \ n    /       n               n   \
grandeur  :    u      = u   + dt | -u   d   u    - ___ d   P - ___ <r>g I         + ___ d   mu s   + g I        + k s I       |  + | Sqdm I        + Pqdm I      |
direction :      i        i       \   i   j   i     r    i      r         i,||g      r    j     ij       i,||g          intf./      \       i,||g           i,+g/
"""

print("Execute le :"+str(datetime.date.today().strftime("%d %b. %Y")))
print("Execute a : "+str(datetime.datetime.now().strftime("%H:%M:%S")))
print("Commandes :"+str(sys.argv[:]))
######################################################################################
# emploi python python_for_qdm_patch_and_source.py chemin_fichier_out combien_de_runs from_scratch [temps0,temps1]
#                                                      1                      2           3              4
# 1 : DOIT contenir RUN0X et le nom du fichier acceleration.out, explicitement.
# 2 : Combien de RUNXX sont ecrits dans ce fichiers out. 
#     Si combien_de_run<0 l'info est cherchee dans chemin_fichier
# 3 : Est-ce que on a l'en-tete du ficher ? (si oui alors le premier ss pdt est a rk=2)
#     scratch -> on a l'en-tete
# /!\ 4: [temps0,temps1] -> releve les donnees des outs entre les 2 temps
# /!\ 4: [temps0,temps0] -> releve la donnee des outs a l'instant le plus proche de temps0
#         [-1,-1] -> releve tous les temps disponibles
# ATTENTION : pour utiliser ce script il faut etre au niveau des dossiers RUN0X car on se sert du nom du RUN pour obtenir des informations
######################################################################################

########################################################################
# PIERRES DE ROSETTE pour savoir quelle colonne correspond a quelle grandeur
#              car dans le cas d'une reprise il se peut que l'on perde cette information
#              meme dans le cas d'un calcul initial, il a de fortes chances pour que ce
#              ne soit pas "bien" fait.
# - rosette_acceleration : correspondance colonne-grandeur pour le fichier jdd_acceleration.out
# - rosette_check        : correspondance colonne-grandeur pour le fichier jdd_check_etapes_et_termes.out
# Le seul moyen de verifier ces infrmations est d'ouvrir IJK, car l'en tete de check_et... etait mal renseigne !

rosette_acceleration = { "it":0,
        "t":1,
        "qdm_source"        : {"x":9,"y":None, "z":None},  # qdm_source_           UNIQUEMENT DANS LA DIRECTION DE LA GRAVITE   
        "u_v"               : {"x":10,"y":None, "z":None}, # vap_velocity_tmoy_    UNIQUEMENT DANS LA DIRECTION DE LA GRAVITE
        "u_l"               : {"x":11,"y":None, "z":None}, # liq_velocity_tmoy_    UNIQUEMENT DANS LA DIRECTION DE LA GRAVITE
        "qdm_patch"         : {"x":None,"y":13, "z":14}    # qdm_patch_correction_ UNIQUEMENT DANS LES AUTRES DIRECTIONS QUE CELLE DE LA GRAVITE
}

rosette_check       = { "it":0,
        "t":1,
        "rk":2,
        "ru_av_pred"         : {"x":3,  "y":4,  "z":5},   # rho_u_euler_av_pred
        "rdu_ap_pred"        : {"x":6,  "y":7,  "z":8},   # rho_du_euler_ap_pred
        "ru_ap_proj"         : {"x":9,  "y":10, "z":11},  # rho_u_euler_av_prediction
        "rdu_ap_proj"        : {"x":12, "y":13, "z":14},  # rho_du_euler_ap_prediction
        "ru_av_rmi"          : {"x":15, "y":16, "z":17},  # rho_u_euler_av_rho_mu_ind
        "ru_ap_rmi"          : {"x":18, "y":19, "z":20},  # rho_u_euler_ap_rho_mu_ind
        "u_ap_rmi"           : {"x":21, "y":22, "z":23},  # u_euler_ap_rho_mu_ind
        "t_intf"             : {"x":24, "y":25, "z":26},  # terme_interfaces               /!\/!\ il est nul celui-la, on le remplit que sous certaines conditions... ne le regardons pas
        "t_conv"             : {"x":27, "y":28, "z":29},  # terme_convection : defaut -> Moyenne_spatiale{ \rho S [u_i u_j] }   @\rho du/dt || /!\ au 24.01.22 : @ \r du/dt / vol.cell.
        "t_diff"             : {"x":30, "y":31, "z":32},  # terme_diffusion : defaut -> Moyenne_spatiale{ S [\mu s_{ij}] } *    @\rho du/dt || /!\ au 24.01.22 : @ \r du/dt / vol.cell.
        "t_pr_2"             : {"x":33, "y":34, "z":35},  # terme_pression_bis : Moyenne_spatiale{ grad(p*) }                   @\rho du/dt || /!\ au 24.01.22 : @ \r du/dt / vol.cell.
        "t_pr_3"             : {"x":36, "y":37, "z":38},  # terme_pression_ter : Moyenne_spatiale{ 1/rho * grad(p*) }           @     du/dt || /!\ au 24.01.22 : @    du/dt 

        "t_intf_bf_ms_bis"   : {"x":39, "y":40, "z":41},  # terme_interfaces_bf_mass_solveur_bis : Moyenne_spatiale{ dv-tc-td }                 @\rho du/dt
        "t_intf_bf_ms"       : {"x":42, "y":43, "z":44},  # terme_interfaces_bf_mass_solver : Moyenne_spatiale{ terme_source_interfaces_ns_ }   @\rho du/dt
        "t_intf_af_ms"       : {"x":45, "y":46, "z":47},  # terme_interfaces_af_mass_solver : Moyenne_spatiale{ terme_source_interfaces_ns_ }   @      du/dt (le mass_solver a ete applique)
        "pres_ap_proj"       : {"x":48, "y":48, "z":48},  # pression_ap_proj, ce n'est qu'UN double
        "t_conv_mass_sol"    : {"x":52, "y":53, "z":54},  # terme_moyen_convection_mass_solver_ : Moyenne_spatiale{ 1/rho * C }  @      du/dt || /!\ au 24.01.22 : @ \r du/dt / vol.cell.
        "t_diff_mass_sol"    : {"x":55, "y":56, "z":57}   # terme_moyen_diffusion_mass_solver_  : Moyenne_spatiale{ 1/rho * D }  @      du/dt || /!\ au 24.01.22 : @ \r du/dt / vol.cell.
}
######################################################################################

# = DONNEES UTILISATEUR
taux_vide = 0.03
g = 9.81
chemin_acc=sys.argv[1]
combien_de_runs = int(sys.argv[2])
from_scratch=(bool(sys.argv[3]=="scratch") or bool(sys.argv[3]=="yes"))
t_deb, t_fin = float(sys.argv[4].rstrip("]").lstrip("[").split(",")[-2]),float(sys.argv[4].rstrip("]").lstrip("[").split(",")[-1])
print("# t_deb = %s \n# t_fin = %s"%(str(t_deb),str(t_fin)))

# = CONSTRUCTION DE DONNEES
chemin_check=chemin_acc.replace("acceleration","check_etapes_et_termes")
jdd = chemin_acc.split('/')[-1].split('_acceleration')[0]
jdd = chemin_acc.split('RUN')[0]+'RUN'+chemin_acc.split('RUN')[-1][:2]+'/'+jdd+'.data'
schema_temps = dtool.getParam(jdd, "time_scheme", string=True)

volume_cellule = 1.
for letter in ["i","j","k"]:
    volume_cellule *= dtool.getParam(jdd, "uniform_domain_size_"+letter, string=True)
for letter in ["i","j","k"]:
    volume_cellule /= dtool.getParam(jdd, "nbelem_"+letter, string=True)
print("volume d'une cellule : "+str(volume_cellule))

At =    dtool.getParam(jdd, "rho_liquide", string=True) - dtool.getParam(jdd, "rho_vapeur", string=True)
At /= ( taux_vide*dtool.getParam(jdd, "rho_liquide", string=True) + (1-taux_vide)*dtool.getParam(jdd, "rho_vapeur", string=True) )

print("coeff immo :", end="", flush=True)
try :
    coeff_immobilisation = dtool.getParam(jdd, "coeff_imobilisation")
except:
    coeff_immobilisation = 0
if coeff_immobilisation != 0:
    # Si on a un coefficient d'immobilisation non nul, il est ecrit dans acceleration.out AVANT notre patch
    rosette_acceleration["qdm_patch"]["y"] += 3
    rosette_acceleration["qdm_patch"]["z"] += 3

# = CHARGEMENT DONNEES : temps et iterations par ACCELERATION.OUT
i = np.loadtxt(chemin_acc)[:,0]

def supprime_les_sous_iterations_du_schema_RK3(chemin_fichier,i):
    """
    # Si simu RK3, on a les valeurs aux trois sous-pas de temps (dont on se moque)
    # -> On va donc supprimer les valeurs du premier et deuxieme sous pas de temps
    #    pour ne garder plus que les valeurs correspondant au rk3 = 2
    """
    # if schema_temps=="RK3":
        # Recupere le numero de run. Les noms des RUN sont : INIT, RUN00, RUN01, ...
    if combien_de_runs < 0:
        if "INIT" in chemin_fichier:
            nrun=0
        else:
            nrun=int(chemin_fichier.split("RUN")[-1][:2])+2
    else:
        nrun=combien_de_runs
    if from_scratch:
        # On a les trois sous pas de temps pour toutes les iterations sauf l'it 0, on ajoute donc deux fois un 0 au debut du vecteur d'iterations
        i=np.insert(i,0,0)
        i=np.insert(i,0,0)
    # Chaque ligne correspond aux it d'un RUN. On reshape pour vectoriser le selecteur
    if nrun!=0:
        i=i.reshape(nrun,int((len(i))/nrun))
        # Pour chaque RUN enregistre, on ne garde qu'une occurence de l'iteration
        uniq_i,indices = np.unique(i,return_index=True,axis=1)
    else:
        # Pour chaque RUN enregistre, on ne garde qu'une occurence de l'iteration
        uniq_i,indices = np.unique(i,return_index=True)
    
    # Remise a plat de la liste des iterations epurees
    i=uniq_i.flatten()
    # La liste indices ne permet de recuperer que pour le premier RUN -le INIT- . 
    # -> Il faut concatener pour obtenir la liste COMPLETE des indices a travers les RUN.
    indices_tot=indices

    for r in range(nrun-1):
        indices_tot=np.concatenate((indices_tot,indices_tot[-len(indices):]+indices[-1]+3))
    
    return(indices_tot,i)

def traitement_t_deb_t_fin(temps,t_deb,t_fin):
    """
    ########################################################################
    # Travail sur les iterations ###########################################
    # /!\ dans le cas de RK3, it_deb/fin sont plutot ligne_deb/fin
    ########################################################################
    """
    ## Recuperation des iterations a considerer dans les fichiers out
    # Si les instants renseignés ne sont pas inclus dans la plage de releves
    #   -> l'instant initial est le premier instant de releve. L'instant final est le dernier instant de releve
    if t_deb<=temps[0] : t_deb = temps[0]; print("Debut des releves : "+str(t_deb))
    if (t_fin<=0 or t_fin>temps[-1]) : t_fin = temps[-1]; print("Fin des releves : "+str(t_fin))
    
    # On associe les instants donnes a l'iteration correspondante
    it_deb = np.argwhere(np.abs(temps-t_deb)==np.min(np.abs(temps-t_deb)))
    #print("1 : ", it_deb)
    """
    ON N'A PAS BESOIN DE CE DETOUR PUISQU'ON APPLIQUE supprime_les_sous_iterations_du_schema_RK3
    AU VECTEUR DES ITERATIONS AVANT DE PASSER ICI
    if "RK3_enftnn" in schema_temps:
        # On souhaite que les tableaux commencent par un rk_step = 2. Un point c'est out
        it=it_deb[0][0]
        vrai_it_deb = np.argwhere( np.abs(temps_2-temps[it]) == np.min(np.abs(temps_2-temps[it])))[0][0]
        #print("vrai_it_deb : ",vrai_it_deb)
        vrai_t_deb = temps_2[vrai_it_deb]
        #print("vrai_t_deb : ",vrai_t_deb)
        it_deb = np.argwhere(np.abs(temps-vrai_t_deb)==np.min(np.abs(temps-vrai_t_deb)))
    #/!\ il faut sans doute faire la meme chose pour it_fin si le fichier out ne se termine pas par un rk = 2
    """
    it_fin = np.argwhere(np.abs(temps-t_fin)==np.min(np.abs(temps-t_fin)))
    
    # Si les instants donnes sont exactement a mi-temsp entre deux instants de releve :
    #    -> on choisi la plus petite valeur de debut et la plus grande valeur de fin
    if len(it_deb)>1:
        it_deb = it_deb[0]
    if len(it_fin)>1:
        it_fin = it_fin[-1]
    
    # Si on a donne le meme instant de debut et de fin :
    #   -> c'est comme si on avait donne un interval de longeur 1
    # if it_fin==it_deb : it_fin +=1
    it_deb, it_fin = int(it_deb[0][0]), int(it_fin[0][0])
    ########################################################################
    return(it_deb,it_fin)

# /!\ SI L'ECRITURE DES DONNES N'EST PLUS FAITE A CHAQUE RK3 PDT,
# IL NE FAUT PAS RENTRER DANS LA CONDITION CI-DESSOUS.
if "RK3" in schema_temps :
    indices_tot,i = supprime_les_sous_iterations_du_schema_RK3(chemin_acc,i)
else :
    indices_tot = i
# /!\ SI L'ECRITURE DES DONNES N'EST PLUS FAITE A CHAQUE RK3 PDT,
# IL NE FAUT PAS RENTRER DANS LA CONDITION CI-DESSUS.
print("Iterations epurees recuperees")

print("Passons au chargement des donnees")
# = CHARGEMENT DONNEES - TEMPS par ACCELERATION.OUT
t = np.loadtxt(chemin_acc)[indices_tot,rosette_acceleration["t"]:rosette_acceleration["t"]+1]
dt = t[1:]-t[0:-1]
to = t
t = t[1:]
# = SELECTION DONNEES - SELECTION DES ITERATIONS ET DES TEMPS DEMANDES PAR UTILISATEUR
it_deb,it_fin = traitement_t_deb_t_fin(t,t_deb,t_fin)
ist = i[it_deb:it_fin]
ts = t[it_deb:it_fin]
tos = to[it_deb:it_fin]

# = CHARGEMENT DONNEES - SOURCE ET PATCH, par ACCELERATION.OUT
print("source, correction colineaire a g ...      ", end="", flush=True)
source = np.loadtxt(chemin_acc)[indices_tot[:],rosette_acceleration["qdm_source"]["x"]:rosette_acceleration["qdm_source"]["x"]+1]
dsource = source[1:]-source[0:-1]
print("OK. size : ",source.shape)
print("patch, correciton orthogonale a g ...      ", end="", flush=True)
patch = np.loadtxt(chemin_acc)[indices_tot[:],rosette_acceleration["qdm_patch"]["y"]:rosette_acceleration["qdm_patch"]["z"]+1]
dpatch = patch[1:,:] - patch[0:-1,:]
print("OK. size : ",patch.shape)

def integration_de_la_grandeur_sur_un_pas_de_temps(quel_fichier=chemin_check,quelle_rosette=rosette_check,quelle_grandeur="t_conv_mass_sol"):
    """
    # Les grandeurs en "terme_bidule" sont evaluees a chaque sous-pas de temps du schema (pour RK3).
    # Or les corrections de QdM se font a la fin de chaque pas de temps GLOBAL, pas de chaque sous-pas de temps.
    # De plus, ces "termes_bidule" ont une signification physique intelligible uniquement s'il sont integres sur tout un
    # pas de temps.
    # Le but de cette fonction est donc de recuperer la valeur de "terme_bidule" a chaque sous-pas de temps du
    # schema temporel, et d'integrer ce "terme_bidule" pour obtenir sa valeur sur le pas de temps COMPLET
    # REMARQUE : D'un terme          @ du/dt 
    #            on passe a un terme @ du
    """
    # Coefficients specifiques au RK3 de IJK, pour ecrire : voir /BILAN_II/build/bilan_II.pdf
    #     v = dt * (dv_0 * a_0 + dv_1 * a_1 + dv_2 * a_2)
    alpha_0 = 1./3. - 25./48. + (153*8.)/(9.*128*3)  # Erreur jusqu'au 24.01.2022 : 128 etait a 28 a cause d'une coquille dans BILAN_II/build/bilan_II.pdf, eq. 59...
    alpha_1 = 45/48. - (153*8.)/(128.*15)         
    alpha_2 = 8/15.                               
    
    ### Les lignes suivantes sont commentees car demandent de creer 3 tableaux inutilement
    # grandeur2 = np.loadtxt(quel_fichier)[indices_tot[1:]-0,quelle_rosette[quelle_grandeur]["x"]:quelle_rosette[quelle_grandeur]["z"]+1]
    # grandeur0 = np.loadtxt(quel_fichier)[indices_tot[1:]-1,quelle_rosette[quelle_grandeur]["x"]:quelle_rosette[quelle_grandeur]["z"]+1]
    # grandeur1 = np.loadtxt(quel_fichier)[indices_tot[1:]-2,quelle_rosette[quelle_grandeur]["x"]:quelle_rosette[quelle_grandeur]["z"]+1]
    
    # grandeur = grandeur0*alpha_0[:,np.newaxis,np.newaxis] + grandeur1*alpha_1[:,np.newaxis,np.newaxis] + grandeur2*alpha_2[:,np.newaxis,np.newaxis]
    ### Les lignes precedentes sont commentees car demandent de creer 3 tableaux inutilement
    
    # Integration
    # je crois qu'il y avait une ligne ou 2 en plus, pour gerer les cas ou on travail avec {la vitesse (ou QdM) directement}
    grandeur =  alpha_2 * np.loadtxt(quel_fichier)[indices_tot[1:]-0,quelle_rosette[quelle_grandeur]["x"]:quelle_rosette[quelle_grandeur]["z"]+1] # Correspond à rk=2
    grandeur += alpha_0 * np.loadtxt(quel_fichier)[indices_tot[1:]-1,quelle_rosette[quelle_grandeur]["x"]:quelle_rosette[quelle_grandeur]["z"]+1] # Correspond à rk=0
    grandeur += alpha_1 * np.loadtxt(quel_fichier)[indices_tot[1:]-2,quelle_rosette[quelle_grandeur]["x"]:quelle_rosette[quelle_grandeur]["z"]+1] # Correspond à rk=1
    grandeur = grandeur*dt[:]
    
    return(grandeur)
   
   
# = QUELLES SONT LES ACTIONS QUE L'ON VA MENER
trace_evolutions_temporelles = True
densites_de_probabilites =  False
trace_densites_de_probabilites =  False
###############################################   
   
   
    
# = CHARGEMENT DES DONNES - TERMES DE CONVECTION {@ du/ (volume cellule)}, DIFFUSION {@ du/ (volume cellule)},
#   PRESSION {@ du}, INTERFACES {@ du/ (volume cellule)}, GRAVITE par CHECK_ETAPES_ET_TEMRES.OUT
print("terme de convection dt * u.grad(u) @ du/ (volume cellule) ...      ", end="", flush=True)
# conv = np.loadtxt(chemin_check)[indices_tot,rosette_check["t_conv_mass_sol"]["x"]:rosette_check["t_conv_mass_sol"]["z"]+1]
conv = integration_de_la_grandeur_sur_un_pas_de_temps(quel_fichier=chemin_check,quelle_rosette=rosette_check,quelle_grandeur="t_conv_mass_sol")
print("OK. size : ",conv.shape)
print("terme de diffusion dt * nu.lapplacien(u) @ du/ (volume cellule)  ...      ", end="", flush=True)
# diff = np.loadtxt(chemin_check)[indices_tot,rosette_check["t_diff_mass_sol"]["x"]:rosette_check["t_diff_mass_sol"]["z"]+1]
diff = integration_de_la_grandeur_sur_un_pas_de_temps(quel_fichier=chemin_check,quelle_rosette=rosette_check,quelle_grandeur="t_diff_mass_sol")
print("OK. size : ",diff.shape)
print("terme des interfaces dt * int_surf {kappa.sigma}@ du/ (volume cellule) ...      ", end="", flush=True)
# intf =  np.loadtxt(chemin_check)[indices_tot,rosette_check["t_intf_af_ms"]["x"]:rosette_check["t_intf_af_ms"]["z"]+1]
intf = integration_de_la_grandeur_sur_un_pas_de_temps(quel_fichier=chemin_check,quelle_rosette=rosette_check,quelle_grandeur="t_intf_af_ms")
print("OK. size : ",intf.shape)
print("terme de pression dt * 1/rho grad(p) @ du  ...      ", end="", flush=True)
# press = np.loadtxt(chemin_check)[indices_tot,rosette_check["t_pr_3"]["x"]:rosette_check["t_pr_3"]["z"]+1]
press = integration_de_la_grandeur_sur_un_pas_de_temps(quel_fichier=chemin_check,quelle_rosette=rosette_check,quelle_grandeur="t_pr_3")
print("OK. size : ",press.shape)
print("la vitesse u @ u...      ", end="", flush=True)
vit = np.loadtxt(chemin_check)[indices_tot[1:],rosette_check["u_ap_rmi"]["x"]:rosette_check["u_ap_rmi"]["z"]+1]
dvit = vit; dvit[1:,:] -= vit[:-1,:]
print("OK. size : ",vit.shape)
#print(conv.shape,source.shape,vit.shape)

# = FONCTIONS POUR UNIFORMISER LE PLOTTING TEMPOREL
def beautiful_two_plots_gr(insert=True, serie1=source, qui1=0, label1="source", couleur1='k', pastille1='oy', serie2=source, qui2=1, label2="source", couleur2='r', pastille2='og'):
    if len(serie1.shape) == 1:
        serie1 = serie1[:,np.newaxis]
    if len(serie2.shape) == 1:
        serie2 = serie2[:,np.newaxis]
        
    serie1s = serie1[it_deb:it_fin,:]
    serie2s = serie2[it_deb:it_fin,:]
    
    fig, (ax1,ax2) = plt.subplots(2,1)
    ins1, ins2  = ax1.inset_axes([0.3, 0.2, 0.6, 0.6]), ax2.inset_axes([0.3, 0.2, 0.6, 0.6])
    
    ax1.plot(t,serie1[:,qui1],couleur1,label=label1);                                      # evolution temporelle de tte la serie
    ax1.plot(t[np.argwhere(i==0)][:,0],serie1[np.argwhere(i==0),qui1],pastille1);         # on marque les reprises d'une pastille
    ins1.plot(ts,serie1s[:,qui1],couleur1)                                                # dans le plot insere on ne plot que la fin de la serie
    ins1.plot(ts[np.argwhere(ist==0)][:,0],serie1s[np.argwhere(ist==0),qui1],pastille1);  # dans le plot insere on marque les reprises d'une pastille
    ins1.tick_params(axis="y",direction="in")#,pad=-0)
    ins1.tick_params(axis="x",direction="in")#,pad=-15)
    ax1.legend();
    
    ax2.plot(t,serie2[:,qui2],couleur2,label=label2);                                       # evolution temporelle de tte la serie
    ax2.plot(t[np.argwhere(i==0)][:,0],serie2[np.argwhere(i==0),qui2],pastille2);          # on marque les reprises d'une pastille
    ins2.plot(ts,serie2s[:,qui2],couleur2)                                                 # dans le plot insere on ne plot qu'une partie de la serie
    ins2.plot(ts[np.argwhere(ist==0)][:,0],serie2s[np.argwhere(ist==0),qui2],pastille2);   # dans le plot insere on marque les reprises d'une pastille
    ins2.tick_params(axis="y",direction="in")#,pad=-0)
    ins2.tick_params(axis="x",direction="in")#,pad=-15)
    ax2.set_xlabel("temps (s)")
    ax2.legend();
    
    for ax in (ax1,ax2):
        ax.tick_params(labelsize=22)
        for label in (ax.get_xticklabels() + ax.get_yticklabels()):
                label.set_fontsize(22)
    
    fig.tight_layout();
    
    name_label1 = ''.join(e for e in label1 if (e.isalnum() or e=="_"))
    name_label2 = ''.join(e for e in label2 if (e.isalnum() or e=="_"))
    fig.savefig(name_label1+"_"+name_label2+".png")
    
    del(fig,ax1,ax2)

def beautiful_two_plots_no_insert_gr(insert=True, serie1=source, qui1=0, label1="source", couleur1='k', pastille1='oy', serie2=source, qui2=1, label2="source", couleur2='r', pastille2='og'):
    # ~~~~~~~~~ Plot toute la serie ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    fig, (ax1,ax2) = plt.subplots(2,1)
    
    ax1.plot(t,serie1[:,qui1],couleur1,label=label1);                                      # evolution temporelle de tte la serie
    ax1.plot(t[np.argwhere(i==0)][:,0],serie1[np.argwhere(i==0),qui1],pastille1);         # on marque les reprises d'une pastille
    ax1.legend();
    
    ax2.plot(t,serie2[:,qui2],couleur2,label=label2);                                       # evolution temporelle de tte la serie
    ax2.plot(t[np.argwhere(i==0)][:,0],serie2[np.argwhere(i==0),qui2],pastille2);          # on marque les reprises d'une pastille
    ax2.set_xlabel("temps (s)")
    ax2.legend();
    
    for ax in (ax1,ax2):
        ax.tick_params(labelsize=22)
        for label in (ax.get_xticklabels() + ax.get_yticklabels()):
                label.set_fontsize(22)
    
    fig.tight_layout();
    
    name_label1 = ''.join(e for e in label1 if (e.isalnum() or e=="_"))
    name_label2 = ''.join(e for e in label2 if (e.isalnum() or e=="_"))
    fig.savefig(name_label1+"_"+name_label2+".png")
    
    del(fig,ax1,ax2)
    # ~~~~~~~~~ Plot uniquement le zoom ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if len(serie1.shape) == 1:
        serie1 = serie1[:,np.newaxis]
    if len(serie2.shape) == 1:
        serie2 = serie2[:,np.newaxis]
        
    serie1s = serie1[it_deb:it_fin,:]
    serie2s = serie2[it_deb:it_fin,:]
    
    fig, (ax1,ax2) = plt.subplots(2,1)
    
    ax1.plot(ts,serie1s[:,qui1],couleur1,label=label1)                                   # dans le plot insere on ne plot que la fin de la serie
    ax1.plot(ts[np.argwhere(ist==0)][:,0],serie1s[np.argwhere(ist==0),qui1],pastille1);  # dans le plot insere on marque les reprises d'une pastille
    ax1.legend();
    
    ax2.plot(ts,serie2s[:,qui2],couleur2,label=label2)                                    # dans le plot insere on ne plot qu'une partie de la serie
    ax2.plot(ts[np.argwhere(ist==0)][:,0],serie2s[np.argwhere(ist==0),qui2],pastille2);   # dans le plot insere on marque les reprises d'une pastille
    ax2.set_xlabel("temps (s)")
    ax2.legend();
    
    for ax in (ax1,ax2):
        ax.tick_params(labelsize=22)
        for label in (ax.get_xticklabels() + ax.get_yticklabels()):
                label.set_fontsize(22)
    
    fig.tight_layout();
    
    name_label1 = ''.join(e for e in label1 if (e.isalnum() or e=="_"))
    name_label2 = ''.join(e for e in label2 if (e.isalnum() or e=="_"))
    fig.savefig(name_label1+"_"+name_label2+"_zoom.png")
    
    del(fig,ax1,ax2)

def beautiful_to_two_plots_gr(insert=True, serie1=source, qui1=0, label1="source", couleur1='k', pastille1='oy', serie2=source, qui2=1, label2="source", couleur2='r', pastille2='og'):
    if len(serie1.shape) == 1:
        serie1 = serie1[:,np.newaxis]
    if len(serie2.shape) == 1:
        serie2 = serie2[:,np.newaxis]
        
    serie1s = serie1[it_deb:it_fin,:]
    serie2s = serie2[it_deb:it_fin,:]
    
    fig, (ax1,ax2) = plt.subplots(2,1)
    ins1, ins2  = ax1.inset_axes([0.3, 0.2, 0.6, 0.6]), ax2.inset_axes([0.3, 0.2, 0.6, 0.6])
    
    ax1.plot(to,serie1[:,qui1],couleur1,label=label1);                                      # evolution temporelle de tte la serie
    ax1.plot(to[np.argwhere(i==0)][:,0],serie1[np.argwhere(i==0),qui1],pastille1);         # on marque les reprises d'une pastille
    ins1.plot(tos,serie1s[:,qui1],couleur1)                                                # dans le plot insere on ne plot que la fin de la serie
    ins1.plot(tos[np.argwhere(ist==0)][:,0],serie1s[np.argwhere(ist==0),qui1],pastille1);  # dans le plot insere on marque les reprises d'une pastille
    ins1.tick_params(axis="y",direction="in")#,pad=-0)
    ins1.tick_params(axis="x",direction="in")#,pad=-15)
    ax1.legend();
    
    ax2.plot(to,serie2[:,qui2],couleur2,label=label2);                                       # evolution temporelle de tte la serie
    ax2.plot(to[np.argwhere(i==0)][:,0],serie2[np.argwhere(i==0),qui2],pastille2);          # on marque les reprises d'une pastille
    ins2.plot(tos,serie2s[:,qui2],couleur2)                                                 # dans le plot insere on ne plot qu'une partie de la serie
    ins2.plot(tos[np.argwhere(ist==0)][:,0],serie2s[np.argwhere(ist==0),qui2],pastille2);   # dans le plot insere on marque les reprises d'une pastille
    ins2.tick_params(axis="y",direction="in")#,pad=-0)
    ins2.tick_params(axis="x",direction="in")#,pad=-15)
    ax2.set_xlabel("temps (s)")
    ax2.legend();
    
    for ax in (ax1,ax2):
        ax.tick_params(labelsize=22)
        for label in (ax.get_xticklabels() + ax.get_yticklabels()):
                label.set_fontsize(22)
    
    fig.tight_layout();
    
    name_label1 = ''.join(e for e in label1 if (e.isalnum() or e=="_"))
    name_label2 = ''.join(e for e in label2 if (e.isalnum() or e=="_"))
    fig.savefig(name_label1+"_"+name_label2+".png")
    
    del(fig,ax1,ax2)

def beautiful_to_two_plots_no_insert_gr(insert=True, serie1=source, qui1=0, label1="source", couleur1='k', pastille1='oy', serie2=source, qui2=1, label2="source", couleur2='r', pastille2='og'):
    # ~~~~ Plot de toute la serie ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    fig, (ax1,ax2) = plt.subplots(2,1)
    
    ax1.plot(to,serie1[:,qui1],couleur1,label=label1);                                      # evolution temporelle de tte la serie
    ax1.plot(to[np.argwhere(i==0)][:,0],serie1[np.argwhere(i==0),qui1],pastille1);         # on marque les reprises d'une pastille
    ax1.legend();
    
    ax2.plot(to,serie2[:,qui2],couleur2,label=label2);                                       # evolution temporelle de tte la serie
    ax2.plot(to[np.argwhere(i==0)][:,0],serie2[np.argwhere(i==0),qui2],pastille2);          # on marque les reprises d'une pastille
    ax2.set_xlabel("temps (s)")
    ax2.legend();
    
    for ax in (ax1,ax2):
        ax.tick_params(labelsize=22)
        for label in (ax.get_xticklabels() + ax.get_yticklabels()):
                label.set_fontsize(22)
    
    fig.tight_layout();
    
    name_label1 = ''.join(e for e in label1 if (e.isalnum() or e=="_"))
    name_label2 = ''.join(e for e in label2 if (e.isalnum() or e=="_"))
    fig.savefig(name_label1+"_"+name_label2+".png")
    del(fig,ax1,ax2)
    
    # ~~~~ Plot du zoom uniquement ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if len(serie1.shape) == 1:
        serie1 = serie1[:,np.newaxis]
    if len(serie2.shape) == 1:
        serie2 = serie2[:,np.newaxis]
        
    serie1s = serie1[it_deb:it_fin,:]
    serie2s = serie2[it_deb:it_fin,:]
    
    fig, (ax1,ax2) = plt.subplots(2,1)
    
    ax1.plot(tos,serie1s[:,qui1],couleur1,label=label1)                                                # dans le plot insere on ne plot que la fin de la serie
    ax1.plot(tos[np.argwhere(ist==0)][:,0],serie1s[np.argwhere(ist==0),qui1],pastille1);  # dans le plot insere on marque les reprises d'une pastille
    ax1.legend();
    
    ax2.plot(tos,serie2s[:,qui2],couleur2,label=label2)                                                 # dans le plot insere on ne plot qu'une partie de la serie
    ax2.plot(tos[np.argwhere(ist==0)][:,0],serie2s[np.argwhere(ist==0),qui2],pastille2);   # dans le plot insere on marque les reprises d'une pastille
    ax2.set_xlabel("temps (s)")
    ax2.legend();
    
    for ax in (ax1,ax2):
        ax.tick_params(labelsize=22)
        for label in (ax.get_xticklabels() + ax.get_yticklabels()):
                label.set_fontsize(22)
    
    fig.tight_layout();
    
    name_label1 = ''.join(e for e in label1 if (e.isalnum() or e=="_"))
    name_label2 = ''.join(e for e in label2 if (e.isalnum() or e=="_"))
    fig.savefig(name_label1+"_"+name_label2+"zoom.png")
    
    del(fig,ax1,ax2)
    
def beautiful_one_plot_gr(insert=True, serie1=source, qui1=0, label1="oui", couleur1='k', pastille1='oy'):
    if len(serie1.shape) == 1:
        serie1 = serie1[:,np.newaxis]
        
    serie1s = serie1[it_deb:it_fin,:]
    
    fig, ax1 = plt.subplots(1,1)
    ins1 = ax1.inset_axes([0.3, 0.2, 0.6, 0.6])
    
    ax1.plot(t,serie1[:,qui1],couleur1,label=label1);                                      # evolution temporelle de tte la serie
    ax1.plot(t[np.argwhere(i==0)][:,0],serie1[np.argwhere(i==0),qui1],pastille1);         # on marque les reprises d'une pastille
    ins1.plot(ts,serie1s[:,qui1],couleur1)                                                # dans le plot insere on ne plot que la fin de la serie
    ins1.plot(ts[np.argwhere(ist==0)][:,0],serie1s[np.argwhere(ist==0),qui1],pastille1);  # dans le plot insere on marque les reprises d'une pastille
    ins1.tick_params(axis="y",direction="in")#,pad=-0)
    ins1.tick_params(axis="x",direction="in")#,pad=-15)
    ax1.legend();
    
    ax1.tick_params(labelsize=22)
    for label in (ax1.get_xticklabels() + ax1.get_yticklabels()):
            label.set_fontsize(22)
    
    fig.tight_layout();
    
    name_label1 = ''.join(e for e in label1 if (e.isalnum() or e=="_"))
    fig.savefig(name_label1+"_"+".png")
    
    del(fig,ax1)

def beautiful_one_plot_no_insert_gr(insert=True, serie1=source, qui1=0, label1="oui", couleur1='k', pastille1='oy'):
    # ~~~ Plot toute la serie ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if len(serie1.shape) == 1:
        serie1 = serie1[:,np.newaxis]
        
    serie1s = serie1[it_deb:it_fin,:]
    
    fig, ax1 = plt.subplots(1,1)
    
    ax1.plot(t,serie1[:,qui1],couleur1,label=label1);                                      # evolution temporelle de tte la serie
    ax1.plot(t[np.argwhere(i==0)][:,0],serie1[np.argwhere(i==0),qui1],pastille1);         # on marque les reprises d'une pastille
    ax1.legend();
    
    ax1.tick_params(labelsize=22)
    for label in (ax1.get_xticklabels() + ax1.get_yticklabels()):
            label.set_fontsize(22)
    
    fig.tight_layout();
    
    name_label1 = ''.join(e for e in label1 if (e.isalnum() or e=="_"))
    fig.savefig(name_label1+"_"+".png")
    del(fig,ax1)
    
    # ~~~ Plot insert uniquement ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if len(serie1.shape) == 1:
        serie1 = serie1[:,np.newaxis]
        
    serie1s = serie1[it_deb:it_fin,:]
    
    fig, ax1 = plt.subplots(1,1)
    
    ax1.plot(ts,serie1s[:,qui1],couleur1,label=label1)                                                # dans le plot insere on ne plot que la fin de la serie
    ax1.plot(ts[np.argwhere(ist==0)][:,0],serie1s[np.argwhere(ist==0),qui1],pastille1);  # dans le plot insere on marque les reprises d'une pastille
    ax1.legend();
    
    ax1.tick_params(labelsize=22)
    for label in (ax1.get_xticklabels() + ax1.get_yticklabels()):
            label.set_fontsize(22)
    
    fig.tight_layout();
    
    name_label1 = ''.join(e for e in label1 if (e.isalnum() or e=="_"))
    fig.savefig(name_label1+"_"+"zoom.png")
    del(fig,ax1)    

def beautiful_to_one_plot_gr(insert=True, serie1=source, qui1=0, label1="oui", couleur1='k', pastille1='oy'):
    if len(serie1.shape) == 1:
        serie1 = serie1[:,np.newaxis]
        
    serie1s = serie1[it_deb:it_fin,:]
    
    fig, ax1 = plt.subplots(1,1)
    ins1 = ax1.inset_axes([0.3, 0.2, 0.6, 0.6])
    
    ax1.plot(to,serie1[:,qui1],couleur1,label=label1);                                      # evolution temporelle de tte la serie
    ax1.plot(to[np.argwhere(i==0)][:,0],serie1[np.argwhere(i==0),qui1],pastille1);         # on marque les reprises d'une pastille
    ins1.plot(tos,serie1s[:,qui1],couleur1)                                                # dans le plot insere on ne plot que la fin de la serie
    ins1.plot(tos[np.argwhere(ist==0)][:,0],serie1s[np.argwhere(ist==0),qui1],pastille1);  # dans le plot insere on marque les reprises d'une pastille
    ins1.tick_params(axis="y",direction="in")#,pad=-0)
    ins1.tick_params(axis="x",direction="in")#,pad=-15)
    ax1.legend();
    
    ax1.tick_params(labelsize=22)
    for label in (ax1.get_xticklabels() + ax1.get_yticklabels()):
            label.set_fontsize(22)
    
    fig.tight_layout();
    
    name_label1 = ''.join(e for e in label1 if (e.isalnum() or e=="_"))
    fig.savefig(name_label1+"_"+".png")
    
    del(fig,ax1)

def beautiful_to_one_plot_no_insert_gr(insert=True, serie1=source, qui1=0, label1="oui", couleur1='k', pastille1='oy'):    
    # ~~~~ Plot de toute la sere ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    fig, ax1 = plt.subplots(1,1)
    
    ax1.plot(to,serie1[:,qui1],couleur1,label=label1);                                      # evolution temporelle de tte la serie
    ax1.plot(to[np.argwhere(i==0)][:,0],serie1[np.argwhere(i==0),qui1],pastille1);         # on marque les reprises d'une pastille
    ax1.legend();
    
    ax1.tick_params(labelsize=22)
    for label in (ax1.get_xticklabels() + ax1.get_yticklabels()):
            label.set_fontsize(22)
    
    fig.tight_layout();
    
    name_label1 = ''.join(e for e in label1 if (e.isalnum() or e=="_"))
    fig.savefig(name_label1+"_"+".png")
    del(fig,ax1)
    # ~~~~ Plot du zoom uniquement ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if len(serie1.shape) == 1:
        serie1 = serie1[:,np.newaxis]
        
    serie1s = serie1[it_deb:it_fin,:]
    
    fig, ax1 = plt.subplots(1,1)
    
    ax1.plot(tos,serie1s[:,qui1],couleur1,label=label1)                                                # dans le plot insere on ne plot que la fin de la serie
    ax1.plot(tos[np.argwhere(ist==0)][:,0],serie1s[np.argwhere(ist==0),qui1],pastille1);  # dans le plot insere on marque les reprises d'une pastille
    ax1.legend();
    
    ax1.tick_params(labelsize=22)
    for label in (ax1.get_xticklabels() + ax1.get_yticklabels()):
            label.set_fontsize(22)
    
    fig.tight_layout();
    
    name_label1 = ''.join(e for e in label1 if (e.isalnum() or e=="_"))
    fig.savefig(name_label1+"_"+"zoom.png")
    
    del(fig,ax1)
    

# = TRACONS LES EVOLUTIONS TEMPORELLES
if trace_evolutions_temporelles:
    # GRAPHIQUE I : qdm_source; selon X
    beautiful_one_plot_no_insert_gr(insert=True, serie1=dsource, qui1=0, label1=r"$\partial_t S_{qdm} \times dt$", couleur1='r', pastille1='og')
    beautiful_to_one_plot_no_insert_gr(insert=True, serie1=source, qui1=0, label1=r"$S_{qdm}$", couleur1='r', pastille1='og')
    # GRAPHIQUE II : qdm_patch; selon Y et Z
    beautiful_two_plots_no_insert_gr(insert=True, serie1=dpatch, qui1=0, label1=r"$\partial_t P_{qdm;y} \times dt$", couleur1='g', pastille1='ob', serie2=dpatch, qui2=1, label2=r"$\partial_t P_{qdm;z} \times dt$", couleur2='b', pastille2='oy')
    beautiful_to_two_plots_no_insert_gr(insert=True, serie1=patch, qui1=0, label1=r"$P_{qdm;y}$", couleur1='g', pastille1='ob', serie2=patch, qui2=1, label2=r"$P_{qdm;z}$", couleur2='b', pastille2='oy')
    # GRAPHIQUE III : conv - diff; selon X
    beautiful_two_plots_no_insert_gr(insert=True, serie1=conv*volume_cellule, qui1=0, label1=r"$\mathcal{V}_c. \times dt \times t_C \vert_X$", couleur1='r', pastille1='og', 
                                                  serie2=diff*volume_cellule, qui2=0, label2=r"$\mathcal{V}_c. \times dt \times t_D \vert_X$", couleur2='r', pastille2='og')
    # GRAPHIQUE IV : conv - diff; selon Y
    beautiful_two_plots_no_insert_gr(insert=True, serie1=conv*volume_cellule, qui1=1, label1=r"$\mathcal{V}_c. \times dt \times t_C \vert_Y$", couleur1='g', pastille1='ob',
                                                  serie2=diff*volume_cellule, qui2=1, label2=r"$\mathcal{V}_c. \times dt \times t_D \vert_Y$", couleur2='g', pastille2='ob')
    # GRAPHIQUE V : intf - press; selon X
    beautiful_two_plots_no_insert_gr(insert=True, serie1=intf*volume_cellule, qui1=0, label1=r"$\mathcal{V}_c. \times dt \times t_I \vert_X$", couleur1='r', pastille1='og', 
                                                  serie2=press,               qui2=0, label2=r"$dt \times t_P \vert_X$", couleur2='r', pastille2='og')
    # GRAPHIQUE VI : intf - press; selon Y
    beautiful_two_plots_no_insert_gr(insert=True, serie1=intf*volume_cellule, qui1=1, label1=r"$\mathcal{V}_c. dt \times t_I \vert_Y$", couleur1='g', pastille1='ob', 
                                                  serie2=press,               qui2=1, label2=r"$dt \times t_P \vert_Y$", couleur2='g', pastille2='ob')
    # ~~~ Bis repetita sans multiplier par volume cellule - Debut ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # GRAPHIQUE III : conv - diff; selon X
    beautiful_two_plots_no_insert_gr(insert=True, serie1=conv, qui1=0, label1=r"$dt \times t_C \vert_X$", couleur1='r', pastille1='og', 
                                                  serie2=diff, qui2=0, label2=r"$dt \times t_D \vert_X$", couleur2='r', pastille2='og')
    # GRAPHIQUE IV : conv - diff; selon Y
    beautiful_two_plots_no_insert_gr(insert=True, serie1=conv, qui1=1, label1=r"$dt \times t_C \vert_Y$", couleur1='g', pastille1='ob',
                                                  serie2=diff, qui2=1, label2=r"$dt \times t_D \vert_Y$", couleur2='g', pastille2='ob')
    # GRAPHIQUE V : intf - press; selon X
    beautiful_two_plots_no_insert_gr(insert=True, serie1=intf, qui1=0, label1=r"$dt \times t_I \vert_X$", couleur1='r', pastille1='og', 
                                                  serie2=press,               qui2=0, label2=r"$dt \times t_P \vert_X$", couleur2='r', pastille2='og')
    # GRAPHIQUE VI : intf - press; selon Y
    beautiful_two_plots_no_insert_gr(insert=True, serie1=intf, qui1=1, label1=r"$dt \times t_I \vert_Y$", couleur1='g', pastille1='ob', 
                                                  serie2=press,               qui2=1, label2=r"$dt \times t_P \vert_Y$", couleur2='g', pastille2='ob')
    # ~~~ Bis repetita sans multiplier par volume cellule - Fin ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # GRAPHIQUE VII : vitesse; selon X et Y
    beautiful_two_plots_no_insert_gr(insert=True, serie1=vit, qui1=0, label1=r"$u_x$", couleur1='r', pastille1='og',serie2=vit, qui2=1, label2=r"$u_y$", couleur2='g', pastille2='ob')
    beautiful_two_plots_no_insert_gr(insert=True, serie1=dvit, qui1=0, label1=r"$du_x$", couleur1='r', pastille1='og',serie2=dvit, qui2=1, label2=r"$du_y$", couleur2='g', pastille2='ob')
    beautiful_two_plots_no_insert_gr(insert=True, serie1=vit*At, qui1=0, label1=r"$At \times u_x$", couleur1='r', pastille1='og',serie2=vit*At, qui2=1, label2=r"$At \times u_y$", couleur2='g', pastille2='ob')
    beautiful_one_plot_no_insert_gr(insert=True, serie1=g*dt, qui1=0, label1=r"$ g \times dt$", couleur1='r', pastille1='og')
    beautiful_one_plot_no_insert_gr(insert=True, serie1=At*g*dt, qui1=0, label1=r"$At \times g \times dt$", couleur1='r', pastille1='og')
    print("Graphiques d'evolution temporelle termines")
    
########################################################################
########################################################################
########################################################################
# = TRAVAIL POUR OBTENIR LES DENSITES DE PROBABILITES
print("Traitement des donnees pour obtenir des densites de probabilites")

def construit_l_histogramme(reduced_sorted_series):
    """
    ## Construction des histogrammes.
    """
    n_bins = 20
    n_series = reduced_sorted_series.shape[-1]
    histos, groups = np.zeros((n_bins,n_series)), np.zeros((n_bins,n_series))
    for i_serie in range(n_series):
        reduced_sorted_serie = reduced_sorted_series[:,i_serie]
        histo, group = np.histogram(reduced_sorted_serie,bins=n_bins)
        somme_histo = np.float64(np.sum(histo))
        histo, group = histo/somme_histo, (group[1:]+group[:-1])/2.
        
        histos[:,i_serie] = histo
        groups[:,i_serie] = group
        
    return(histos,groups)
    
def statistiques_de_la_serie(reduced_sorted_series,histos,groups):
    """
    ### Statistiques des series entres
    """    
    ## Moments 1 et 2
    means = np.mean(reduced_sorted_series,axis=0)
    stds  = np.std(reduced_sorted_series,axis=0)
    ## Indice de la moyenne : Premier indice de sorted_vitesse ou sorted_vitesse vaut sa moyenne 
    ids_means = np.argwhere(
               np.abs(reduced_sorted_series-means[np.newaxis,:])
               ==
               np.abs(reduced_sorted_series-means[np.newaxis,:]).min(axis=0)
                         )
    ## Indice de la mediane : Premier indice de l'histogramme tq l'histogramme vaut son maximum
    ids_medianes = np.argwhere(histos==histos.max(axis=0))
    medianes = groups[ids_medianes[:,0],ids_medianes[:,1]]
    
    return(means,ids_means,stds,medianes,ids_medianes)

# = Tri dans l'ordre croissant (selon l'axe du temps ofc, pas selon l'axe de quelle est la grandeur) pour la construction de l'histogramme
# On ôte de la serie le 1000 eme le plus élevé et le 1000 eme le plus bas
dsource = np.sort(dsource[it_deb:it_fin,:],axis=0)#[int(0.001*it_fin):int((1-0.001)*it_fin)]
h_dsource, g_dsource = construit_l_histogramme(dsource)
dsource_means, dsource_id_means, dsource_stds, dsource_medianes, dsource_id_medianes = statistiques_de_la_serie(dsource,h_dsource,g_dsource)

if densites_de_probabilites :
    dpatch  = np.sort(dpatch[it_deb:it_fin,:] ,axis=0)#[int(0.001*it_fin):int((1-0.001)*it_fin)]
    conv   = np.sort(conv[it_deb:it_fin,:]  ,axis=0)[int(0.001*it_fin):int((1-0.001)*it_fin)]
    diff   = np.sort(diff[it_deb:it_fin,:]  ,axis=0)[int(0.001*it_fin):int((1-0.001)*it_fin)]
    intf   = np.sort(intf[it_deb:it_fin,:]  ,axis=0)[int(0.001*it_fin):int((1-0.001)*it_fin)]
    press  = np.sort(press[it_deb:it_fin,:] ,axis=0)[int(0.001*it_fin):int((1-0.001)*it_fin)]
    grav   = np.sort(grav[it_deb:it_fin,:]  ,axis=0)[int(0.001*it_fin):int((1-0.001)*it_fin)]
    
    h_dpatch, g_dpatch = construit_l_histogramme(dpatch)
    h_conv, g_conv = construit_l_histogramme(conv)
    h_diff, g_diff = construit_l_histogramme(diff)
    h_intf, g_intf = construit_l_histogramme(intf)
    h_press, g_press = construit_l_histogramme(press)
    h_grav, g_grav = construit_l_histogramme(grav)
    
    dpatch_means, dpatch_id_means, dpatch_stds, dpatch_medianes, dpatch_id_medianes = statistiques_de_la_serie(dpatch,h_dpatch,g_dpatch)
    conv_means, conv_id_means, conv_stds, conv_medianes, conv_id_medianes = statistiques_de_la_serie(conv,h_conv,g_conv)
    diff_means, diff_id_means, diff_stds, diff_medianes, diff_id_medianes = statistiques_de_la_serie(diff,h_diff,g_diff)
    intf_means, intf_id_means, intf_stds, intf_medianes, intf_id_medianes = statistiques_de_la_serie(intf,h_intf,g_intf)
    press_means, press_id_means, press_stds, press_medianes, press_id_medianes = statistiques_de_la_serie(press,h_press,g_press)
    grav_means, grav_id_means, grav_stds, grav_medianes, grav_id_medianes = statistiques_de_la_serie(grav,h_grav,g_grav)

    ### Remplacement des valeurs extrêmes par la valeur moyenne de la serie
    # vitesse_x[np.where(vitesse_x<=sorted_vitesse_x[mini])]=moy_x
    # vitesse_y[np.where(vitesse_y<=sorted_vitesse_y[mini])]=moy_y
    # vitesse_z[np.where(vitesse_z<=sorted_vitesse_z[mini])]=moy_z
    
def plot_comme_Dodd(serie1=dsource,group1=g_dsource,histo1=h_dsource,qui1=0,label1="dsourceX",couleur1='r',serie2=dsource,group2=g_dsource,histo2=h_dsource,qui2=0,label2='dsourceX',couleur2='g'):
    if len(serie1.shape) == 1:
        serie1 = serie1[:,np.newaxis]
    if len(serie2.shape) == 1:
        serie2 = serie2[:,np.newaxis]

    t_deb, t_fin = float(t[it_deb]), float(t[it_fin])

    gauche, bas = 0.15, 0.1
    largeur, hauteur = 0.5, 0.5
    espace = 0.005
    rect = [gauche, bas, largeur ,hauteur]
    rect_x = [gauche, bas+hauteur+espace, largeur, 0.3*hauteur]
    rect_y = [gauche+largeur+espace, bas, 0.3*largeur, hauteur]

    fig = plt.figure(figsize=(8,6))
    ax = fig.add_axes(rect)
    ax_hx = fig.add_axes(rect_x, sharex=ax)
    ax_hy = fig.add_axes(rect_y, sharey=ax)
    # nom des axes et titre
    ax.set_xlabel(label1,fontsize=26)
    ax_hx.set_ylabel(r'$\mathcal{PDF}($'+label1+r'$)$',fontsize=18)
    ax.set_ylabel(label2,fontsize=18)
    ax_hy.set_xlabel(r'$\mathcal{PDF}($'+label2+r'$)$',fontsize=18)
    ax_hx.set_title("Repartition de "+label1+" et "+label2)
    # trace des valeurs
    ax_hx.bar(group1[:,qui1],histo1[:,qui1],group1[1,qui1]-group1[0,qui1],color=couleur1,label=label1)
    ax_hx.tick_params(axis="x", labelbottom=False)
    ax_hy.barh(group2[:,qui2],histo1[:,qui2],group2[1,qui2]-group2[0,qui2],color=couleur2,label=label2)
    ax_hy.tick_params(axis="y", labelleft=False)
    ax.scatter(serie1/0.001, serie2/0.001,s=0.02, label="t0={0}, t1={1}, nt={2}".format(str(np.round(t_deb,2)),str(np.round(t_fin,2)),str(it_fin-it_deb)))
    # ajout de lines pour se reperer
    # ax.axvline(x=1, color='k',linewidth=0.5)
    # ax_hx.axvline(x=1, color='k',linewidth=0.5)
    # ax.axhline(y=1, color='k',linewidth=0.5)
    # ax_hy.axhline(y=1, color='k',linewidth=0.5)
    # ajout d'une grille pour se reperer plutot
    """ A FAIRE """
    # un peu de paillettes et cest dans la boite
    ax.legend()
    ax_hx.legend()
    ax_hy.legend()
    plt.legend()
    
    fig.tight_layout();
    
    label1 = ''.join(e for e in label1 if (e.isalnum() or e=="_"))
    label2 = ''.join(e for e in label2 if (e.isalnum() or e=="_"))
    plt.savefig('PDF_'+label1+"_"+label2+'.png')
    

# = TRACONS LES PDF COMME DODD
if trace_densites_de_probabilites :
    plot_comme_Dodd(serie1=dsource,group1=g_dsource,histo1=h_dsource,qui1=0,label1="dsourceX",couleur1='r',serie2=dsource,group2=g_dsource,histo2=h_dsource,qui2=0,label2='sourceX',couleur2='g')
    plot_comme_Dodd(serie1=dpatch,group1=g_dpatch,histo1=h_dpatch,qui1=0,label1="dpatchY",couleur1='r',serie2=dpatch,group2=g_dpatch,histo2=h_dpatch,qui2=1,label2='patchZ',couleur2='g')



"""
### Calcul des Gaussiennes associées aux histogrammes precedents
# 1. Gaussienne centrée sur la moyenne de la série, avec la meme dispersion que la serie
Gauss1_x = Gauss(group_x,moy_x,std_x)#scipy.stats.norm.pdf(group_x,group_x[id_moy_x][-1],std_x) * (group_x[1:]-group_x[:-1])[-1]
Gauss1_y = Gauss(group_y,moy_y,std_y)#scipy.stats.norm.pdf(group_y,group_y[id_moy_y][-1],std_y) * (group_y[1:]-group_y[:-1])[-1]
Gauss1_z = Gauss(group_z,moy_z,std_z)#scipy.stats.norm.pdf(group_z,group_z[id_moy_z][-1],std_z) * (group_z[1:]-group_z[:-1])[-1]
# 2. Gaussiene centrée sur la valeur médiane, avec la meme dispersion que la serie
Gauss2_x = Gauss(group_x,med_x,std_x)#scipy.stats.norm.pdf(group_x,group_x[id_max_x],np.std(histo_x)) * (group_x[1:]-group_x[:-1])[-1]
Gauss2_y = Gauss(group_y,med_y,std_y)#scipy.stats.norm.pdf(group_y,group_y[id_max_y],np.std(histo_y)) * (group_y[1:]-group_y[:-1])[-1]
Gauss2_z = Gauss(group_z,med_z,std_z)#scipy.stats.norm.pdf(group_z,group_z[id_max_z],np.std(histo_z)) * (group_z[1:]-group_z[:-1])[-1]
# 3. Gaussiene centrée sur la valeur maximale de la serie, avec la meme valeur max
G3_std_x = group_x[id_max_histo_x][-1] #1./ (np.sqrt(2*np.pi)*group_x[id_max_histo_x][-1])
G3_std_y = group_y[id_max_histo_y][-1] #1./ (np.sqrt(2*np.pi)*group_y[id_max_histo_y][-1])
G3_std_z = group_z[id_max_histo_z][-1] #1./ (np.sqrt(2*np.pi)*group_z[id_max_histo_z][-1])
Gauss3_x = Gauss(group_x,med_x,G3_std_x) # scipy.stats.norm.pdf(group_x,group_x[id_max_x],G3_std_x) * (group_x[1:]-group_x[:-1])[-1]
Gauss3_y = Gauss(group_y,med_y,G3_std_y) # scipy.stats.norm.pdf(group_y,group_y[id_max_y],G3_std_y) * (group_y[1:]-group_y[:-1])[-1]
Gauss3_z = Gauss(group_z,med_z,G3_std_z) # scipy.stats.norm.pdf(group_z,group_z[id_max_z],G3_std_z) * (group_z[1:]-group_z[:-1])[-1]
"""

# = TRACES DE DENSITES DE PROBABILITES (comme Dodd, a finir de piocher du python dans /volatile/RANDOM_TRIPERIO/SOME_PYTHON/vitesses_bulles_Dodd.py
# plt.figure()
# plt.bar(g_source[:,0],h_source[:,0],g_source[1,0]-g_source[0,0],label="source")
# plt.legend()
# plt.figure()
# plt.bar(g_patch[:,0],h_patch[:,0],g_patch[1,0]-g_patch[0,0],label="patch_y")
# plt.legend()
# plt.figure()
# plt.barh(g_patch[:,1],h_patch[:,1],g_patch[1,1]-g_patch[0,1],label="patch_z")
# plt.legend()


# plt.show()


