import matplotlib.pyplot as plt
import numpy as np
import sys
import DNSTools3 as dtool

"""
L'objectif de ce script est de :
 - Afficher l'evolution temporelle de la correction orthogonale a g (le patch de qdm)
 - Afficher l'evolution temporelle de la correction colineaire a g (la source de qdm)
 - Afficher la PDF des deux corrections precedentes, donnant la moyenne, et la variance
 - Comparer ces grandeurs de facon pertinente avec une autre grandeur du probleme.
"""
######################################################################################
# emploi python python_for_qdm_patch_and_source.py chemin_fichier_out nb_iteration_a_ignorer 
#                                                      1                    2                
# 1 : DOIT contenir RUN0X et le nom du fichier, explicitement.
# 2 : un entier
# ATTENTION : pour utiliser ce script il faut etre au niveau des dossiers RUN0X car on se sert du nom du RUN pour obtenir des informations
######################################################################################
######################################################################################
# Bonne note : 
# x[0] -> qdm_source, x[1] -> uv, x[2] -> ul, x[3] -> qdm_patch[0], x[4] -> qdm_patch[1], x[5] -> qdm_patch[2] 
######################################################################################

# DONNEES UTILISATEUR
chemin_acc=sys.argv[1]
start=int(sys.argv[2])
#schema_temps=sys.argv[3]

# CONSTRUCTION DE DONNEES
chemin_check=chemin_acc.replace("acceleration","check_etapes_et_termes")
jdd = chemin_acc.split('/')[-1].split('_acceleration')[0]
jdd = chemin_acc.split('RUN')[0]+'RUN'+chemin_acc.split('RUN')[-1][:2]+'/'+jdd+'.data'
schema_temps = dtool.getParam(jdd, "time_scheme", string=True)

# CHARGEMENT DONNEES - ACCELERATION.OUT
x,t,i = np.loadtxt(chemin_acc)[:,-6:], np.loadtxt(chemin_acc)[:,1], np.loadtxt(chemin_acc)[:,0]
# CHARGEMENT DONNEES - CHECK_ETAPES_ET_TERMES.OUT
conv = np.loadtxt(chemin_check)[:,-4:]

def supprime_les_sous_iterations_du_schema_RK3(chemin_fichier):
    # Si simu RK3, on a les valeurs aux trois sous-pas de temps (dont on se moque)
    # -> On va donc supprimer les valeurs du premier et deuxieme sous pas de temps
    #
    # if schema_temps=="RK3":
        # Recupere le numero de run. Les noms des RUN sont : INIT, RUN00, RUN01, ...
        if "INIT" in chemin_fichier:
            nrun=0
        else:
            nrun=int(chemin_fichier.split("RUN")[-1][:2])+2
        # On a les trois sous pas de temps pour toutes les iterations sauf l'it 0, on ajoute donc une fois un 0 au debut du vecteur d'iterations
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
        
        return(indices_tot)

# ATTENTION SI L'ECRITURE DES DONNES N'EST PLUS FAITE A CHAQUE RK3 PDT,
# IL NE FAUT PLUS RENTRER DANS LA CONDITION DE SES MORTS CI-DESSOUS.
if schema_temps=="RK3":
    indices_tot = supprime_les_sous_iterations_du_schema_RK3(chemin_acc)
    # Selection des donnees du dernier sous pas de temps.
    x,t=x[indices_tot,:],t[indices_tot]
    dcp=dcp[indices_tot,:]
# ATTENTION SI L'ECRITURE DES DONNES N'EST PLUS FAITE A CHAQUE RK3 PDT, 
# IL NE FAUT PLUS RENTRER DANS LA CONDITION DE SES MORTS CI-DESSUS.

# Selection partielle
xs = x[start:,:]
udcps = udcp[start,:]
d = x[:,1] - x[:,2]
ds = x[start:,1] - x[start:,2]
e = x[:,0] - x[:,2]
es = x[start:,0] - x[start:,2]
ts,ist = t[start:],i[start:]
######################################################################################
def beautiful_plot_gr(insert=True, serie1=x, qui1=0, label1="oui", couleur1='k', pastille1='oy', serie2=x, qui2=1, label2="non", couleur2='r', pastille2='og'):
    if len(serie1.shape) == 1:
        serie1 = serie1[:,np.newaxis]
    if len(serie2.shape) == 1:
        serie2 = serie2[:,np.newaxis]
        
    serie1s = serie1[start:,:]
    serie2s = serie2[start:,:]
    
    fig, (ax1,ax2) = plt.subplots(2,1)
    ins1, ins2  = ax1.inset_axes([0.3, 0.2, 0.6, 0.6]), ax2.inset_axes([0.3, 0.2, 0.6, 0.6])
    
    ax1.plot(t,serie1[:,qui1],couleur1,label=label1);                                      # evolution temporelle de tte la serie
    ax1.plot(t[np.argwhere(i==0)][:,0],serie1[np.argwhere(i==0),qui1],pastille1);         # on marque les reprises d'une pastille
    ins1.plot(ts,serie1s[:,qui1],couleur1)                                                # dans le plot insere on ne plot que la fin de la serie
    ins1.plot(ts[np.argwhere(ist==0)][:,0],serie1s[np.argwhere(ist==0),qui1],pastille1);  # dans le plot insere on marque les reprises d'une pastille
    ins1.tick_params(axis="y",direction="in")#,pad=-0)
    ins1.tick_params(axis="x",direction="in")#,pad=-15)
    ax1.legend();
    
    ax2.plot(t,x[:,qui2],couleur2,label=label2);                                       # evolution temporelle de tte la serie
    ax2.plot(t[np.argwhere(i==0)][:,0],x[np.argwhere(i==0),qui2],pastille2);          # on marque les reprises d'une pastille
    ins2.plot(ts,serie2s[:,qui2],couleur2)                                                 # dans le plot insere on ne plot que la fin de la serie
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
    
    label1 = ''.join(e for e in label1 if (e.isalnum() or e=="_"))
    label2 = ''.join(e for e in label2 if (e.isalnum() or e=="_"))
    fig.savefig(label1+"_"+label2+".png")
    
    del(fig,ax1,ax2)
    
def beautiful_plot_gr_bis(insert=True, serie1=d, label1="oui", couleur1='g', pastille1='ob', serie2=e, label2="non", couleur2='y', pastille2='ok'):
    """Comme beautiful_plot_gr mais pour un tableau a une colonne au lieu d'etre pour un tableau a pls colonnes. """
    
    serie1s = serie1[start:]
    serie2s = serie2[start:]
    
    fig2, (ax12,ax22) = plt.subplots(2,1)
    ins12, ins22  = ax12.inset_axes([0.3, 0.2, 0.6, 0.6]), ax22.inset_axes([0.3, 0.2, 0.6, 0.6])
    
    ax12.plot(t,serie1,couleur1,label='u_r');
    ax12.plot(t[np.argwhere(i==0)][:,0],serie1[np.argwhere(i==0)][:,0],pastille1);
    ins12.plot(ts,serie1s,couleur1)
    ins12.plot(ts[np.argwhere(ist==0)][:,0],serie1[np.argwhere(ist==0)][:,0],pastille1);
    ins12.tick_params(axis="y",direction="in")#,pad=-0)
    ins12.tick_params(axis="x",direction="in")#,pad=-15)
    ax12.legend();
    
    ax22.plot(t,serie2,couleur2,label='qdm-u_l');
    ax22.plot(t[np.argwhere(i==0)][:,0],serie2[np.argwhere(i==0)][:,0],pastille2);
    ins22.plot(ts,serie2s,couleur2)
    ins22.plot(ts[np.argwhere(ist==0)][:,0],serie2s[np.argwhere(ist==0)][:,0],pastille2);
    ins22.tick_params(axis="y",direction="in")#,pad=-0)
    ins22.tick_params(axis="x",direction="in")#,pad=-15)
    ax22.set_xlabel("temps (s)")
    ax22.legend();
    
    del(fig,ax1,ax2)

######################################################################################
# Bonne note : 
# x[0] -> qdm_source, x[1] -> uv, x[2] -> ul, x[3] -> qdm_patch[0], x[4] -> qdm_patch[1], x[5] -> qdm_patch[2] 
######################################################################################
#####################################################################################
# GRAPHIQUE I : uv - qdm_source
beautiful_plot_gr(insert=True,qui1=1,label1="u_v",qui2=0,label2="qdm_source")

########################################################################################
# GRAPHIQUE II : uv - ul
beautiful_plot_gr(True,x,1,"u_v",'k','oy',x,2,"u_l",'b','om')

########################################################################################
# GRAPHIQUE III : ur=uv-ul - qdm_source-ul
beautiful_plot_gr(True,d,0,"u_r",'g','ob',e,0,"u_l-qdm",'y','ok')
# beautiful_plot_gr_bis(True,d,1,"u_v",'k','oy',x,2,"u_l",'b','om')

########################################################################################
# GRAPHIQUE IV : qdm_qource adimentionne : qdm_qource/(terme_force_init) = qdm_qource/(rho_m g) et qdm_qource/(u_v')
terme_force_init = float(dtool.getParam(jdd, "terme_force_init"))
if terme_force_init != 0:
    qdm_adim0 = x[:,0] / terme_force_init 
    qdm_adim0s = xs[:,0] / terme_force_init 
# qdm_adim0 = x[:,0] / dtool.getParam("terme_force_init",jdd) 

beautiful_plot_gr(True,qdm_adim0,0,r'qdm/$\overline{\rho} g$','g','ob',e,0,"u_l-qdm",'y','ok')

########################################################################################
# GRAPHIQUE V : qdm - terme_convection
beautiful_plot_gr(insert=True, serie1=x, qui1=0, label1="q_source", couleur1='r', pastille1='og', serie2=udcp, qui2=2, label2='conv.', couleur2='b', pastille2='og')

########################################################################################
# GRAPHIQUE VI : qdm - terme_diffusion
beautiful_plot_gr(insert=True,serie1=x,qui1=0, label1="q_source", couleur1='r', pastille1='og', serie2=udcp, qui2=1, label2='diff.', couleur2='b', pastille2='og')

########################################################################################
# GRAPHIQUE VII : qdm - terme_pression
beautiful_plot_gr(insert=True,serie1=x,qui1=0, label1="q_source", couleur1='r', pastille1='og', serie2=udcp, qui2=3, label2='intf.', couleur2='b', pastille2='og')

########################################################################################
# GRAPHIQUE VIII : qdm - terme_gravite
# terme_gravite=

# MONTRE TOUS
plt.show()

