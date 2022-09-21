import numpy as np
import math as mt
import matplotlib.pyplot as plt
import DNSTools3 as dtool
import glob
import sys
import datetime

"""
    PDF de la distribution des distance entre les bulles
"""
# emploi : python paire_dist_inter_bulles.py chemin_fichier_out mode [bulle_min,bulle_max] [tdeb,tfin]:pas_it

print("Execute le :"+str(datetime.date.today().strftime("%d %b. %Y")))
print("Execute a : "+str(datetime.datetime.now().strftime("%H:%M:%S")))
print("Commandes :"+str(sys.argv[:]))

# ######################################################################
# ######################################################################
# DONNEES UTILISATEUR
chemin = sys.argv[1]
mode = sys.argv[2]
from_out = "out" in mode
from_raw = "raw" in mode
id_min,id_max = int(sys.argv[3].rstrip("]").lstrip("[").split(",")[-2]),int(sys.argv[3].rstrip("]").lstrip("[").split(",")[-1])
t_deb, t_fin, pas_t = float(sys.argv[4].split(':')[0].rstrip("]").lstrip("[").split(",")[-2]),float(sys.argv[4].split(':')[0].rstrip("]").lstrip("[").split(",")[-1]), int(sys.argv[4].split(':')[1])
print("# t_deb = %s \n# t_fin = %s"%(str(t_deb),str(t_fin)))


def traitement_t_deb_t_fin(temps,t_deb,t_fin):
    """
    ########################################################################
    # Travail sur les iterations ###########################################
    # /!\ dans le cas de RK3, it_deb/fin sont plutot ligne_deb/fin
    ########################################################################
    """
    ## Recuperation des iterations a considerer dans les fichiers out
    # Si les instants renseignes ne sont pas inclus dans la plage de releves
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

# ######################################################################
# ######################################################################
# LECTURE DES DONNEES
test=False
# Longueurs : dmaine, diametre bulle, force repulsion, ...
jdd = glob.glob("*data")[0]
fichier = "*_centre_dir.out"
L_domaine = dtool.getParam(jdd,"uniform_domain_size_i")
try : repulsion_range = dtool.getParam(jdd,"portee_force_repulsion")
except : repulsion_range = 0
d_bulle = (6*dtool.getParam(jdd,"vol_bulle_monodisperse")/np.pi)**(1./3.)

if from_out:
    print("FROM OUT")
    # Temps
    t = np.loadtxt(glob.glob(chemin+fichier.replace("dir.out","x.out"))[0],usecols=(0))[::pas_t]
    it_deb,it_fin = traitement_t_deb_t_fin(t,t_deb,t_fin)
    t = t[it_deb:it_fin]
    print("Lignes chargees : ",it_fin-it_deb)
    n_instants = t.shape[0]
    # Nombre de bulles
    with open(glob.glob(chemin+fichier.replace("dir.out","x.out"))[0], 'r') as f:
        line = f.readline() # read 1 line
        n_col = len(line.split(' '))
    if id_min<=0 : id_min=1
    if id_min>n_col-1 : id_min=n_col-2
    if id_max<=0 : id_max=n_col
    if id_max==1 : id_max=2
    if id_max>n_col-1 : id_max=n_col-1
    print("Colonnes dans le fichier : ",n_col)
    print("Colonnes chargees : ",id_max-id_min)
    colonnes = np.arange(id_min,id_max)
    n_bulles = id_max-id_min
    n_paire_bulles = int( ((n_bulles-1)*(n_bulles)) / 2)
    print("Nombre de bulles : ",n_bulles," = ",id_max-id_min)
    print("Nombre de paires de bulles : ",n_paire_bulles)
    # Positions
    # X
    X = np.loadtxt(glob.glob(chemin+fichier.replace("dir.out","x.out"))[0],usecols=colonnes)[it_deb:it_fin,:]
    # Y
    Y = np.loadtxt(glob.glob(chemin+fichier.replace("dir.out","y.out"))[0],usecols=colonnes)[it_deb:it_fin,:]
    # Z
    Z = np.loadtxt(glob.glob(chemin+fichier.replace("dir.out","z.out"))[0],usecols=colonnes)[it_deb:it_fin,:]
    # AJUSTEMENT 
    # Si en plus d'avoir reduit avec les colonnes on veut encore reduire le nb de bulles 
    id_min_bulles = 0
    id_max_bulles = n_bulles-1

# ######################################################################
# ######################################################################
# TRAITEMENT DES DONNEES
def generate_sample(nb_sample):
    """
    tirer nb_sample position (de bulles) pour test
    en sortir les distances, plotter les distributions
    """
    return(sample_pos)
# = CONSTRUCTION DES PAIRES
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def get_paire_distances_sup(pos_X,pos_Y,pos_Z):
    """ NON UTILISEE AU 24.03.22
     Tableau des distance inter-bulles
      -> En visualisant le tableau a double entree {bulleXbulle}, 
      cette focntion sort la partie superieure strict du tableau.
      Entree : Trois vecteurs position de bulle [nombre_instants,nombre_bulles]
      Sortie : Tableau des distances [nombre_instant,nombre_de_paire_de_bulles, type_de_distance_enregistree]
               ->[:,:,0] : distances dnas le repere spherique
               ->[:,:,1] : distances dans le repere cylindrique, normal a X
               ->[:,:,2:]: distances selon l'axe {X,Y,Z}
    """
    store_distance = np.zeros((pos_X.shape[0],n_paire_bulles,4))
    cpt2=int( (2*id_max_bulles-id_min_bulles)*(id_min_bulles)/2)
    for i1 in range(id_min_bulles,id_max_bulles):
        # print("i1 : ",i1)
        for i2 in range(i1+1,id_max_bulles+1):
            # Ecart dimension par dimension, entre -L/2 et L/2
            dX = (pos_X[:,i2] - pos_X[:,i1])%(L_domaine / 2)
            dY = (pos_Y[:,i2] - pos_Y[:,i1])%(L_domaine / 2)
            dZ = (pos_Z[:,i2] - pos_Z[:,i1])%(L_domaine / 2)
            # distance : Norme du vecteur b1b2
            distance = np.sqrt(dX**2 + dY**2 + dZ**2)
            # Construction indice
            # Rouge -> Idince pour le parcourt d'une matrice caree
            rouge = n_bulles*i1 + i2
            # Gris -> ... La diagonale
            gris = i1+1
            # Vert -> ... Le triangle inferieur
            vert = i1*(i1+1)/2
            # idx -> Indice pour le parcourt du triangle superieur de la matrice
            idx0 = int( rouge - gris - vert) 
            store_distance[:,idx0,0] = distance
            store_distance[:,idx0,1] = dX
            store_distance[:,idx0,2] = dY
            store_distance[:,idx0,3] = dZ
            if test:
                print("i1,i2,idx0",i1,i2,idx0)
                print("rouge : ",rouge)
                print("gris : ",gris)
                print("vert : ",vert)
    return(store_distance)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def get_paire_distances_inf(pos_X,pos_Y,pos_Z):
    """
     Tableau des distance inter-bulles
      -> En visualisant le tableau a double entree {bulleXbulle}, 
      cette focntion sort la partie inferieure strict du tableau.
      Entree : Trois vecteurs position de bulle [nombre_instants,nombre_bulles]
      Sortie : Tableau des distances [nombre_instant,nombre_de_paire_de_bulles, type_de_distance_enregistree]
               ->[:,:,0] : distances dnas le repere spherique
               ->[:,:,1] : distances dans le repere cylindrique, normal a X
               ->[:,:,2:]: distances selon l'axe {X,Y,Z}
    """
    print("Computing pairs : (%)   ", end="",flush=True)
    # Tableau des distance inter-bulles
    P,old_P=0.,0.
    store_distance = np.zeros((pos_X.shape[0],n_paire_bulles,5))
    cpt2=int( (2*id_max_bulles-id_min_bulles)*(id_min_bulles)/2)
    for i1 in range(id_min_bulles+1,id_max_bulles+1):
        for i2 in range(id_min_bulles,i1):
            P+= (100./(n_paire_bulles+1))
            if P >= old_P+10 : 
                print(str(np.round(P,1))+" ... ", end="",flush=True)
                old_P = P
            cpt2+=1
            # Indice dans le triangle inferieur d'une matrice carree
            idx0 = int( i1*(i1-1)/2) + i2
            # Ecart dimension par dimension, entre -L/2 et L/2
            dX = np.abs(pos_X[:,i2] - pos_X[:,i1])
            dY = np.abs(pos_Y[:,i2] - pos_Y[:,i1])
            dZ = np.abs(pos_Z[:,i2] - pos_Z[:,i1])
            # distance : Norme du vecteur b1b2
            distance = np.sqrt(dX**2 + dY**2 + dZ**2)
            distanceP = np.sqrt(dY**2 + dZ**2)
            """ ajouter la distance comme ce que fait Remi, les bandeaux sur sphere"""
            store_distance[:,idx0,0] = distance  
            store_distance[:,idx0,1] = distanceP 
            store_distance[:,idx0,2] = dX
            store_distance[:,idx0,3] = dY
            store_distance[:,idx0,4] = dZ
            if test:
                print("i1,i2,idx0",i1,i2,idx0)
                print("vert : ",i1*(i1-1)/2)
    print("Done")
    return(store_distance)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if from_out:
    distances = get_paire_distances_inf(X,Y,Z)
    shapes_distances_DXYZ = distances.shape
    np.save(chemin+"shapes_distances_DXYZ",shapes_distances_DXYZ)
    distances.tofile(chemin+"/distances_DXYZ"+".raw")
elif from_raw:
    print("FROM RAW")
    n_instants,n_paire_bulles,n_distances = np.load(chemin+"shapes_distances_DXYZ"+".npy")
    print("Nombre d'instants : ",n_instants)
    print("Nombre de paires de bulles : ",n_paire_bulles)
    print("Nombre de types de distances : ",n_distances)
    print(n_instants,n_paire_bulles,n_distances)
    distances = np.zeros((n_instants,n_paire_bulles,n_distances))
    distances = np.fromfile(glob.glob(chemin+"/distances_DXYZ"+".raw")[0]).reshape(n_instants,n_paire_bulles,n_distances)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
def modifie_distances_avant_reshape(store_distances):
    """
        Les distance relevees sont :
              (-L < )distance < L
        Or le domaine est (tri-)periodique. Les distances reelles sont : 
              (-L/2 < )distance < L/2
        Cette fonction ramene les distances entre 0 et L/2,
        pour une entree : [n_releves, nombre_de_paires, nb_distanes]
    """
    for idx0,_ in enumerate(store_distances[0,:,0]):
        # Distance X
        store_distances[:,idx0,2] = np.min((L_domaine-store_distances[:,idx0,2]),store_distances[:,idx0,2],axis=0)
        # Distance Y                                                                                     
        store_distances[:,idx0,3] = np.min((L_domaine-store_distances[:,idx0,3]),store_distances[:,idx0,3],axis=0)
        # Distance Z                                                                                     
        store_distances[:,idx0,4] = np.min((L_domaine-store_distances[:,idx0,4]),store_distances[:,idx0,4],axis=0)
        # Distance dans lespace
        store_distances[:,idx0,0] = np.sqrt(store_distances[:,idx0,2]**2+store_distances[:,idx0,3]**2+store_distances[:,idx0,4]**2)
        # Distance dans le plan (Y,Z) 
        store_distances[:,idx0,1] = np.sqrt(store_distances[:,idx0,3]**2+store_distances[:,idx0,4]**2)
    return(store_distances)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~            
def modifie_distances_apres_reshape(store_distances):
    """
        Les distance relevees sont :
              (-L < )distance < L
        Or le domaine est (tri-)periodique. Les distances reelles sont : 
              (-L/2 < )distance < L/2
        Cette fonction ramene les distances entre 0 et L/2,
        pour une entree : [n_releves*nombre_de_paires, nb_distanes]
    """
    print(store_distances.shape)
    # Distance X
    store_distances[:,2] = np.min(((L_domaine-store_distances[:,2]),store_distances[:,2]),axis=0)
    # Distance Y                                                                         
    store_distances[:,3] = np.min(((L_domaine-store_distances[:,3]),store_distances[:,3]),axis=0)
    # Distance Z                                                                         
    store_distances[:,4] = np.min(((L_domaine-store_distances[:,4]),store_distances[:,4]),axis=0)
    # Distance dans l'espace
    store_distances[:,0] = np.sqrt(store_distances[:,2]**2+store_distances[:,3]**2+store_distances[:,4]**2)
    # Distance dans le plan (Y,Z) 
    store_distances[:,1] = np.sqrt(store_distances[:,3]**2+store_distances[:,4]**2)
    return(store_distances)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~            
if test:
    # AFFICHAGE POUR TESTS
    y_min,y_max, instant = 0.2, 0.8, 0
    bulle_min = X[instant,id_min_bulles]
    
    fig_pos, (ax_pos_brut,ax_pos_relatif) = plt.subplots(2,1,figsize=(20,10))
    fig_dis, (ax_dis_brut,ax_dis_relatif) = plt.subplots(2,1,figsize=(20,10))
    # Titres
    ax_pos_brut.set_title("brut");    ax_dis_brut.set_title("brut")
    ax_pos_relatif.set_title("relatif");    ax_dis_relatif.set_title("relatif")
    # Positions    
    for b in range(id_min_bulles,id_max_bulles):
        bulle = X[instant,b]
        print("X[%s] = "%b,bulle)
        # Position des bulles
        ax_pos_brut.axvline(bulle,                   color='blue')
        ax_pos_brut.text(bulle,y_min,r'$b_{%s}$'%(b),color="k",fontsize=20)
        # Position relativement a la premiere bulle
        ax_pos_relatif.axvline(bulle-bulle_min,                      color='blue')
        ax_pos_relatif.text   (bulle-bulle_min,y_min,r'$b_{%s}$'%(b),color="k",fontsize=20)
    ax_pos_brut.axvline(bulle_min,color='k')
    ax_pos_relatif.axvline(0.,color='k')
    # Limites du domaine
    for ax in [ax_pos_brut,ax_pos_relatif]:
        for limite in [-L_domaine/2,L_domaine/2]:
            ax.axvline(limite,color='k',linestyle='dotted')
    fig_pos.savefig("pos.eps")
    
    # Distances
    store_db=[]
    for b in range(id_min_bulles+1,id_max_bulles):
        db = int( i1*(i1-1)/2) + i2
        distance_bulle = distances[instant,db,1]
        store_db.append(distance_bulle)
        print("dX[%s] = "%db,distance_bulle)
        # distances entre la premiere bulle et les autres bulles (b_autre-id_min_bulles...)
        ax_dis_brut.axvline(distance_bulle,color='r')    
        ax_dis_brut.text(distance_bulle,y_max,r'$\Delta_{{%s,%s}}$'%(id_min_bulles,b),color="r",fontsize=20)
        # distance, ramenee a la position de la premiere bulle (id_min_bulles + b_autre-id_min_bulles) : retombe sur posiiton brute
        ax_dis_relatif.axvline(bulle_min+distance_bulle,color='r')    
        ax_dis_relatif.text(bulle_min+distance_bulle,y_max,r'$\Delta_{{%s,%s}}$'%(id_min_bulles,b),color="r",fontsize=20)
    ax_dis_relatif.axvline(bulle_min+0,color='r')    
    ax_dis_relatif.text(bulle_min+0,y_max,r'$\Delta_{{%s,%s}}$'%(id_min_bulles,id_min_bulles),color="r",fontsize=20)
    # Limites du domaine
    for ax in [ax_dis_brut,ax_dis_relatif]:
        for limite in [-L_domaine/2,L_domaine/2]:
            ax.axvline(limite,color='k',linestyle='dotted')    
    fig_dis.savefig("dist.eps")

    print("###################################")
    nt,nb=3,1
    d=distances[nt,nb,2]
    print("distances[%s,%s,:] : "%(nt,nb),distances[nt,nb,:])

# = MANIPULAIOTN DES DONNEES : perte de l'information "iteration" et "numero paire bulle"
print("Dimension matrices des distances avant reshape : ",distances.shape)
distances = distances.reshape(-1,distances.shape[-1])
print("Dimension matrices des distances apres reshape : ",distances.shape)

if test:
    print("distances[%s:%s,:] : "%(n_paire_bulles*nt,n_paire_bulles*nt+n_paire_bulles))
    print(distances[n_paire_bulles*nt:n_paire_bulles*nt+n_paire_bulles,:])

# = CONSTRUCTION DENSITES DE PROBABILITES
print("Traitement des donnees pour obtenir des densites de probabilites")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def construit_l_histogramme(reduced_sorted_series):
    """
    ## Construction des histogrammes.
    Entree : Tableau de n_series a n_releves
             [n_releves,n_series]
    Sortie : histograme et groupes
             [n_bins,n_series] et [n_bins,n_series]
    """
    n_bins = 100
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
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
def construit_l_histogramme_compense(reduced_sorted_series):
    """
    ## Construction des histogrammes : repartition des distances des bulles
    Les distances sdnas le repere spheriques sont normalisees par la surface
    de la sphere de rayon correspondant.
    Les distances sdnas le repere cylindrique sont normalisees par le perimetre
    du cylindre de rayon correspondant.
    Entree : Tableau de n_series a n_releves
             [n_releves,n_series]
    Sortie : histograme et groupes
             [n_bins,n_series] et [n_bins,n_series]
    ERREUR : TypeError: ufunc 'true_divide' output (typecode 'd') could not be coerced to provided output parameter (typecode 'l') according to the casting rule ''same_kind''
    ABANDON -> fonciton modifie_histogramme fait la meme chose
    """
    n_bins = 200
    n_series = reduced_sorted_series.shape[-1]
    histos, groups = np.zeros((n_bins,n_series)), np.zeros((n_bins,n_series))
    for i_serie in range(n_series):
        reduced_sorted_serie = reduced_sorted_series[:,i_serie]
        histo, group = np.histogram(reduced_sorted_serie,bins=n_bins)
        somme_histo = np.float64(np.sum(histo))
        print(histo)
        print(group)
        # Corrections
        if i_serie == 0 : histo /= (4*np.pi*(group[1:])**2)
        if i_serie == 1 : histo /= (2*np.pi*group[1:])
        histo, group = histo/somme_histo, (group[1:]+group[:-1])/2.
        
        histos[:,i_serie] = histo
        groups[:,i_serie] = group
        
    return(histos,groups)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
def ecrit_l_histogramme(histos,groups,nom_fichier):
    """
    ## Ecriture des histogrammes des series statistiques relevees
    # grandeur pdf(grandeur)
    """
    n_series = histos.shape[-1]
    longueur_serie = histos.shape[0]
    data = np.zeros((longueur_serie,2*n_series))
    head_info = "1.r 2.pdf(|d|) 3.r_(y;z) 4.pdf(|d_(y;z)|) 5.x 6.pdf(dx) 7.y 8.pdf(dy) 9.z 10.pdf(dz)"
    for i_serie in range(n_series):
        data[:,2*i_serie]   = groups[:,i_serie]
        data[:,2*i_serie+1] = histos[:,i_serie]
    np.savetxt(nom_fichier,data,fmt='%1.3e',delimiter=' ',header=head_info)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def ecrit_les_statistiques(means,ids_means,stds,medianes,ids_medianes,nom_fichier):
    """ NULLE
    ## Ecriture des statistiques des histogrammes des series relevees
    # moyenne indice_moyenne mediane indice_mediande ecart-type
    """
    n_series = 5
    n_stats = 5
    data = np.zeros((n_series,7))
    
    # data[:,0] = ["sphe","cyl","x","y","z"]
    data[:,1] = means
    data[:,2] = ids_means[:][0]
    data[:,3] = ids_means[:][1]
    data[:,4] = medianes
    data[:,5] = ids_medianes[:][0]
    data[:,6] = ids_medianes[:][1]
    data[:,7] = stds
    head_info = "1.mean 2.id_mean 3.med 4.id_med 5.sq_dev"

    np.savetxt(nom_fichier,data,fmt='%1.3e',delimiter=' ',header=head_info)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def plot_histo(histo,group,cest_quoi,nom):
    plt.plot(group,histo,color='g',label=r'PDF(d('+cest_quoi+r'))')
    plt.bar(group,histo,group[1]-group[0],color='g',alpha=0.5)
    plt.axvline(x=L_domaine/2              ,color='k',linestyle='dotted',linewidth=0.9)
    plt.axvline(x=d_bulle/2                ,color='b',linestyle='solid',linewidth=0.5)
    plt.axvline(x=d_bulle/2+repulsion_range,color='r',linestyle='dotted',linewidth=0.5)
    plt.xlim(-0.001,1.1*L_domaine/2)
    plt.xlabel(cest_quoi)
    plt.ylabel(r'PDF(d('+cest_quoi+r'))')
    
    plt.savefig("hist"+nom+".pdf")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
def plot_cloud_cartesien_2D(x_data,y_data,x_cest_quoi,y_cest_quoi,nom):
    plt.scatter(x_data, y_data,s=0.02, label="t0={0}, t1={1}, nt={2}".format(str(np.round(t_deb,2)),str(np.round(t_fin,2)),str(n_instants)))
    plt.xlabel(x_cest_quoi)
    plt.ylabel(y_cest_quoi)
    
    plt.savefig("cloud"+nom+".pdf")
    
from mpl_toolkits.mplot3d import Axes3D
def plot_cloud_cartesien_3D(x_data,y_data,z_data,x_cest_quoi,y_cest_quoi,z_cest_quoi,nom):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(x_data, y_data, z_data ,s=0.09, label="t0={0}, t1={1}, nt={2}".format(str(np.round(t_deb,2)),str(np.round(t_fin,2)),str(n_instants)))
    ax.set_xlabel(x_cest_quoi)
    ax.set_ylabel(y_cest_quoi)
    plt.savefig("cloud"+nom+".pdf")

def plot_3D_and_histo(xyz_data,xyz_histo,xyz_group,nom):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_title("PDF distance")
    ax.scatter(xyz_data[:,0], xyz_data[:,1], xyz_data[:,2] ,s=0.09, label="t0={0}, t1={1}, nt={2}".format(str(np.round(t_deb,2)),str(np.round(t_fin,2)),str(n_instants)))
    ax.plot(xyz_group[:,0],xyz_histo[:,0],0*xyz_group[:,0],color="r",label="PDF(X)")
    ax.plot(0*xyz_group[:,1],xyz_group[:,1],xyz_histo[:,1],color="g",label="PDF(Y)")
    ax.plot(xyz_histo[:,2],0*xyz_group[:,2],xyz_group[:,2],color="k",label="PDF(Z)")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.legend()
    plt.show()
   
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
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
qui0=0;qui1=1;qui2=2;qui3=3;qui4=4
noms = ["Sphe","Cyl","X","Y","Z"]
quoi =      [r'$|r|=\sqrt{x^2 + y^2 + z^2}$',r'$|r_{y,z}|=\sqrt{y^2 + z^2}$',r'$x$',         r'$y$',         r'$z$']
quoi_mod1 = [r'$|r|=\sqrt{x^2 + y^2 + z^2}$',r'$|r_{y,z}|=\sqrt{y^2 + z^2}$',r'$min(L-x;x)$',r'$min(L-y,y)$',r'$min(L-z,z)$']
quoi_mod2 = [r'$|r|=\sqrt{x^2 + y^2 + z^2}$',r'$|r_{y,z}|=\sqrt{y^2 + z^2}$',r'$min(L-x;x)$',r'$min(L-y,y)$',r'$min(L-z,z)$']

# HISTOGRAMME DES DISTANCES BRUTES
h_distance, g_distance = construit_l_histogramme(distances)
#triche pour normer
h_distance = h_distance/(np.sum(h_distance*1*(g_distance<L_domaine/2.),axis=0))
ecrit_l_histogramme(h_distance, g_distance,"crude_histo.txt")
for qui in range(0,4):
    plot_histo(h_distance[:,qui],g_distance[:,qui],quoi[qui],noms[qui])
    plt.close()

# HISTOGRAMME DES DISTANCES CORRIGEES NON COMPENSEES
distances = modifie_distances_apres_reshape(distances)
h_distance, g_distance = construit_l_histogramme(distances)
# triche pour normer
h_distance = h_distance/(np.sum(h_distance*1*(g_distance<L_domaine/2.),axis=0))
for qui in range(0,4):
    plot_histo(h_distance[:,qui],g_distance[:,qui],quoi[qui],noms[qui]+"_mod1")
    plt.close()
    
# HISTOGRAMME DES DISTANCES CORRIGEES COMPENSEES
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def modifie_histogramme(histos,groups):
    # Histograme de la distance dans l'espace : division par surface sphere
    histos[:,0] /= (4.*np.pi*groups[:,0]**2)
    # Histogramme de la distance dans le plan Y,Z : divison par la longueur perimetre
    histos[:,1] /= (2.*np.pi*groups[:,1])
    return(histos,groups)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
h_distance, g_distance = modifie_histogramme(h_distance, g_distance)
# triche pour normer
h_distance = h_distance/(np.sum(h_distance*1*(g_distance<L_domaine/2.),axis=0))
means,ids_means,stds,medianes,ids_medianes = statistiques_de_la_serie(distances,h_distance,g_distance)
ecrit_l_histogramme(h_distance, g_distance,"tuned_histo.txt")
# ecrit_les_statistiques(means,ids_means,stds,medianes,ids_medianes,"tuned_histo_s_stats.txt")
for qui in range(0,4):
    plot_histo(h_distance[:,qui],g_distance[:,qui],quoi[qui],noms[qui]+"_mod2")
    plt.close()
    
    
pas_de_bulles = int(n_paire_bulles/300)
plot_cloud_cartesien_2D(distances[::pas_de_bulles,3],distances[::pas_de_bulles,4],"Y","Z","YZ")
plt.close()
# visu 3D
plot_cloud_cartesien_3D(distances[::pas_de_bulles,2],distances[::pas_de_bulles,3],distances[::pas_de_bulles,4],"X","Y","Z","YZ")
plt.close()
plot_3D_and_histo(distances[::pas_de_bulles,2:5],h_distance[:,2:5],g_distance[:,2:5],"XYZ")
    
# FAIRE UN PLOT MASHALLA ++ comme Dodd
def plot_comme_Dodd(serie1=distances,group1=g_distance,histo1=h_distance,qui1=0,label1="distanceX",couleur1='r',serie2=distances,group2=g_distance,histo2=h_distance,qui2=0,label2='distanceX',couleur2='g'):
    if len(serie1.shape) == 1:
        serie1 = serie1[:,np.newaxis]
    if len(serie2.shape) == 1:
        serie2 = serie2[:,np.newaxis]

    # t_deb, t_fin = float(t[it_deb]), float(t[it_fin])

    gauche, bas = 0.15, 0.1
    largeur, hauteur = 0.5, 0.5
    espace = 0.005
    rect = [gauche, bas, largeur ,hauteur]
    rect_x = [gauche, bas+hauteur+espace, largeur, 0.3*hauteur]
    rect_y = [gauche+largeur+espace, bas, 0.3*largeur, hauteur]

    fig = plt.figure(figsize=(8,8))
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
    ax.scatter(serie1[:,qui1], serie2[:,qui2],s=0.02, label="t0={0}, t1={1}, nt={2}".format(str(np.round(t_deb,2)),str(np.round(t_fin,2)),str(n_instants)))
    # un peu de paillettes et cest dans la boite
    ax.legend()
    ax_hx.legend()
    ax_hy.legend()
    plt.legend()
    
    # fig.tight_layout();
    
    label1 = ''.join(e for e in label1 if (e.isalnum() or e=="_"))
    label2 = ''.join(e for e in label2 if (e.isalnum() or e=="_"))
    plt.savefig('PDF_'+label1+"_"+label2+'.pdf')
    
# = TRACONS LES PDF COMME DODD
trace_densites_de_probabilites=True
if trace_densites_de_probabilites :
    plot_comme_Dodd(serie1=distances[::pas_de_bulles,:],group1=g_distance,histo1=h_distance,qui1=3,label1="dY",couleur1='r',serie2=distances[::pas_de_bulles,:],group2=g_distance,histo2=h_distance,qui2=4,label2='dZ',couleur2='g')
    plot_comme_Dodd(serie1=distances[::pas_de_bulles,:],group1=g_distance,histo1=h_distance,qui1=2,label1="dX",couleur1='r',serie2=distances[::pas_de_bulles,:],group2=g_distance,histo2=h_distance,qui2=3,label2='dY',couleur2='g')
    plot_comme_Dodd(serie1=distances[::pas_de_bulles,:],group1=g_distance,histo1=h_distance,qui1=1,label1=r'$\sqrt{dy^2+dz^2}$',couleur1='r',serie2=distances[::pas_de_bulles,:],group2=g_distance,histo2=h_distance,qui2=2,label2='dX',couleur2='g')
    
