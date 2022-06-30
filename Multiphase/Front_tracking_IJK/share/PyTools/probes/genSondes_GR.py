# -*- coding: utf-8 -*-
import os
import readline
import sys
import glob
import numpy as np
from argparse import ArgumentParser, RawDescriptionHelpFormatter



# ######################################################################
# Pour bien s'assurer de ce qui est fait : aller voir 
# IJK_FT : chercher "sondes" --> sondes_
# IJK_FT_Post : --> les_sondes_
# Sonde_IJK : dans Entrees&; case 6, on a verifie que la variable les_positions_
#             contient bien les points de relevés en mode segments
# AAAAARG gen sonde il est cassé quoi !!!!!!!!!
# ==> Le C++ releve aux bon endroits les valeurs des champs a priori
# ==> Pourquoi les x,y,z relevé dans les fichiers *.son n'ont aucun sens ???
# ######################################################################
# mode d'emploi : genSondes.py data_file.data -aX -pX -nX -fF -mM
#     -> nombre de sondes-segments alignées selon X
#     -> nombre de sondes-segments alignées selon Y
#     -> nombre de sondes-segments alignées selon Z
#     -> F : nom du fichier ou ecrire le texte pour les sondes
#     -> M : ".", "v" ou "+" pour indiquer quelle disposition de sondes on veut
# # ######################################################################
readline.parse_and_bind('tab:complete')
proj = os.getenv("project_directory")
dir_module = os.path.join(proj, "share", "PyTools")
sys.path.append(dir_module)
dir_module = os.path.join(dir_module, "commons")
sys.path.append(dir_module)
import DNSTools as dtool


def process_cmd_line():
    text = """Génération d'un bloc de texte de sonde pour trioIJK."""
    parser = ArgumentParser(description=text,
                            formatter_class=RawDescriptionHelpFormatter)

    # lJDD = glob.glob("*.data")
    # if len(tJDD) == 1:
    #     defJDD = lJDD[0]
    # else:
    #     defJDD = ""
    ###################################################################
    # Dans les help : "nX" pour plan normal à X, "nY" pour plan normal à Y, "nZ" pour plan normal à Z.
    ###################################################################
    # Nom du JDD
    parser.add_argument(dest="jdd", default="DNS.data",
                        help="Jeu de données trio par une virgule.",
                        type=str)

    # pas de temps
    parser.add_argument("-t", "--dt", dest="dt", default=1e-6,
                        help="Numéro de fichiers à ouvrir si plusieurs résultats.") # GAB comprend plutot qu cest la frequence décriture sur la sonde

    # Nombre de sondes axiales
    # Nombre de sondes-segments alignees selon X
    parser.add_argument("-a", "--axial", dest="axial", default=1,
                        help="Racine caree du nombre de segment selon Ox : nombre de plan nY et nZ contenant un segment Ox",
                        type=int)

    # Nombre de plans yOz où mettre des sondes Oy et Oz
    parser.add_argument("-p", "--plan", dest="plan", default=1,
                        help="Nombre de plans nX contenant des segments-Oy et Oz (ils son tdistribues sur les meme plan nX).",type=int)

    # Nombre de sondes par plan
    parser.add_argument("-n", "--splan", dest="splan", default=1,
                        help="Nombre de plan nZ contenant un segment-Oy.\nAussi : Nombre de plan nY contenant un segment-Oz.",
                        type=int)

    # Nombre de sondes axiales
    parser.add_argument("-f", "--file", dest="nomfile", default="/tmp/sondes.txt",
                        help="Nom du fichier de destination de sondes.")
                        
    # Mode : v, +, .
    parser.add_argument("-m", "--mode", dest="mode", default=".",
                        help="Mode de disposition des sondes : coin 'v', point '.', croix '+'.")                        

    options = parser.parse_args()
    return options
# Fin process_cmd_line


def sondePoint(nom, var, dt, coord):
    '''Definit une sonde point pour un fichier de donnees Trio'''
    sep = " "
    text = nom + sep + var + sep + str(dt) + sep + "points" + \
           sep + 1 + str(coord[0]) + sep + str(coord[1]) + sep + str(coord[2])
    return text
# Fin sonde point


def sondeSegment(nom, var, dt, npoints, coordDebut, coordFin):
    '''Definit une sonde segment pour un fichier de donnees Trio'''
    sep = " "
    text = nom + sep + var + sep + "periode" + sep + str(dt) + sep + "Segment" + \
           sep + str(npoints) + sep + str(coordDebut[0]) + sep + \
           str(coordDebut[1]) + sep + str(coordDebut[2]) + sep + \
           str(coordFin[0]) + sep + str(coordFin[1]) + sep + str(coordFin[2])
    return text
# Fin sonde segment


if __name__ == '__main__':
    options = process_cmd_line()
    jdd = options.jdd
    dt = options.dt
    axial = options.axial
    plan = options.plan
    splan = options.splan
    nomfic = options.nomfile
    mode = options.mode
    # if len(sys.argv) != 2:
    #     raise Exception("usage: python %s jdd.data\n Datafile is required to read the domain.\nOutput is in /tmp/sondes.txt" % sys.argv[0])
    
    # Distionnaire pour associer une valeur de decalage à un message "o", "p", "m" qui sera dans le nom de la sonde
    rosette_mode = {0:"o", 1:"p", -1:"m"}
    # jdd = sys.argv[1]
    print("Reading domain information from: ", jdd)
    Lx = dtool.getParam(jdd, 'uniform_domain_size_i')
    Ly = dtool.getParam(jdd, 'uniform_domain_size_j')
    Lz = dtool.getParam(jdd, 'uniform_domain_size_k')
    Nx = dtool.getParam(jdd, 'nbelem_i'); Nx = int(Nx)
    Ny = dtool.getParam(jdd, 'nbelem_j'); Ny = int(Ny)
    Nz = dtool.getParam(jdd, 'nbelem_k'); Nz = int(Nz)
    try :
        Ox = dtool.getParam(jdd, 'origin_i')
    except :
        Ox = 0.0
    try :
        Oy = dtool.getParam(jdd, 'origin_j')
    except :
        Oy = 0.0
    try :
        Oz = dtool.getParam(jdd, 'origin_k')
    except :
        Oz = 0.0
    #
    deltaEspace = Lx / Nx
    # Sondes dans des plans yOz. a NsondeX[0] hauteures differentes, on
    # place NsondeY[1] sondes selon Oy et NsondeZ[2] selon Oz
    # NsondeX = [9, 4, 4]
    NsondeX = [plan, splan, splan]  # -p -n -n
    # Sondes selon l'axe Ox. On place des segments sur NsondeY points en
    # Oy et NsondeZ en Oz
    NsondeY = axial # correspond au -a
    NsondeZ = axial # correspond au -a
    
    # Decalages pour construire les "faisceaux e sonde plutot que juste les segments seuls. 
    # 1 pour se decalder de +1 maille; -1 por se décaler de -1 maille
    if mode=="v" : decalages = np.array([[1, 0], [0, 0], [0,1]])
    elif mode=="+" : decalages = np.array([[1, 0], [0, 0], [-1, 0], [0, 1], [0, -1]])
    else : decalages = np.array([[0, 0]])

    print("###########################################")
    print("mode : %s"%(mode))
    print("###########################################")
    print(" You want %iX%i Ox-probes in the Y-Z section"%(axial,axial))
    print(" You want %i plans of Oy-probes and Oz-probes"%(plan))
    print(" For each plan, you want %i Oy-probes or %i Oz-probes"%(splan,splan))
    print(" You want %iX%i Oy-probes in the X-Z section"%(plan,splan))
    print(" You want %iX%i Oz-probes in the X-Y section"%(plan,splan))
    print("###########################################")

    ####################
    # Script
    dx = [Lx / (plan + 1), Ly / (splan + 1), Lz / (splan + 1)]
    dy = Ly / (axial + 1)
    dz = Lz / (axial + 1)
    # Petit decalage pour eviter d'avoir de tomber exactement sur les faces
    epsilon_x = Lx/float(Nx) * 1/3.
    epsilon_y = Ly/float(Ny) * 1/3.
    epsilon_z = Lz/float(Nz) * 1/3.
    # En vrai la il faut faire attn si on a mis le centre au milieu de la boite hein...
    xmin = Ox + epsilon_x
    xmax = Lx + Ox
    ymin = Oy + epsilon_y
    ymax = Ly + Oy
    zmin = Oz + epsilon_z
    zmax = Lz + Oz


    print("###########################################")
    print(" The domain goes from (%f,%f,%f) \n\t\t to (%f,%f,%f)"%(xmin,ymin,zmin,xmax,ymax,zmax))
    print(" The space step is (%f, %f, %f)"%(Lx/Nx, Ly/Ny, Lz/Nz))
    print("###########################################")

    NL = "\n  "
    texte = "  Sondes" + NL
    texte += "{" + NL
    
    for nd, d in enumerate(decalages):
        marqueur_decalage0 = rosette_mode[d[0]]
        marqueur_decalage1 = rosette_mode[d[1]]
        texte += "# -> Sondes %s-%s #"%(marqueur_decalage0, marqueur_decalage1) + NL
        
        #### Sondes segment-Ox
        texte += "# Sondes le long de l'axe Ox #" + NL
        for i in range(1, NsondeY + 1):
            for j in range(1, NsondeZ + 1):
                posY = ymin + (i)*dy + d[0]*Ly/float(Ny)
                posZ = zmin + (j)*dz + d[1]*Lz/float(Nz)
                nomSonde = "  S_X_" + str(i) +str(j)+marqueur_decalage0+marqueur_decalage1+ "_VX"
                # sondeSegment(nom, var, dt, npoints, coordDebut, coordFin)
                texte += sondeSegment(nomSonde, "velocity_x", dt, Nx,
                                    [xmin, posY, posZ], [xmax, posY, posZ])
                texte += NL
                nomSonde = "  S_X_" + str(i) +str(j)+marqueur_decalage0+marqueur_decalage1+ "_VY"
                texte += sondeSegment(nomSonde, "velocity_y", dt, Nx,
                                    [xmin, posY, posZ], [xmax, posY, posZ])
                texte += NL
                nomSonde = "  S_X_" + str(i) +str(j)+marqueur_decalage0+marqueur_decalage1+ "_VZ"
                texte += sondeSegment(nomSonde, "velocity_z", dt, Nx,
                                    [xmin, posY, posZ], [xmax, posY, posZ])
                texte += NL
                nomSonde = "  S_X_" + str(i) +str(j)+marqueur_decalage0+marqueur_decalage1+ "_CHI"
                texte += sondeSegment(nomSonde, "indicatrice", dt, Nx,
                                    [xmin, posY, posZ], [xmax, posY, posZ])
                texte += NL
                
        #### Sondes segment-Oy
        texte += "# Sondes le long de l'axe Oy #" + NL
        for i in range(1, NsondeX[0] + 1):
            for j in range(1, NsondeX[1] + 1):
                posZ = zmin + (j)*dx[2] + d[0]*Lz/float(Nz)
                posX = xmin + (i)*dx[0] + d[1]*Lx/float(Nx)
                nomSonde = "  S_Y_" + str(i) +str(j)+marqueur_decalage0+marqueur_decalage1+ "_VX"
                # sondeSegment(nom, var, dt, npoints, coordDebut, coordFin)
                texte += sondeSegment(nomSonde, "velocity_x", dt, Ny,
                                    [posX, ymin, posZ], [posX, ymax, posZ])
                texte += NL
                nomSonde = "  S_Y_" + str(i) +str(j)+marqueur_decalage0+marqueur_decalage1+ "_VY"
                texte += sondeSegment(nomSonde, "velocity_y", dt, Ny,
                                    [posX, ymin, posZ], [posX, ymax, posZ])
                texte += NL
                nomSonde = "  S_Y_" + str(i) +str(j)+marqueur_decalage0+marqueur_decalage1+ "_VZ"
                texte += sondeSegment(nomSonde, "velocity_z", dt, Ny,
                                    [posX, ymin, posZ], [posX, ymax, posZ])
                texte += NL
                nomSonde = "  S_Y_" + str(i) +str(j)+marqueur_decalage0+marqueur_decalage1+ "_CHI"
                texte += sondeSegment(nomSonde, "indicatrice", dt, Ny,
                                    [posX, ymin, posZ], [posX, ymax, posZ])
                texte += NL
        #### Sondes segment-Oz
        texte += "# Sondes le long de l'axe Oz #" + NL
        for i in range(1, NsondeX[0] + 1):
            for j in range(1, NsondeX[2] + 1):
                posX = xmin + (i)*dx[0] + d[0]*Lx/float(Nx)
                posY = ymin + (j)*dx[1] + d[1]*Lx/float(Nx)
                nomSonde = "  S_Z_" + str(i) +str(j)+marqueur_decalage0+marqueur_decalage1+ "_VX"
                texte += sondeSegment(nomSonde, "velocity_x", dt, Nz,
                                    [posX, posY, zmin], [posX, posY, zmax])
                texte += NL
                nomSonde = "  S_Z_" + str(i) +str(j)+marqueur_decalage0+marqueur_decalage1+ "_VY"
                texte += sondeSegment(nomSonde, "velocity_y", dt, Nz,
                                    [posX, posY, zmin], [posX, posY, zmax])
                texte += NL
                nomSonde = "  S_Z_" + str(i) +str(j)+marqueur_decalage0+marqueur_decalage1+ "_VZ"
                texte += sondeSegment(nomSonde, "velocity_z", dt, Nz,
                                    [posX, posY, zmin], [posX, posY, zmax])
                texte += NL
                nomSonde = "  S_Z_" + str(i) +str(j)+marqueur_decalage0+marqueur_decalage1+ "_CHI"
                texte += sondeSegment(nomSonde, "indicatrice", dt, Nz,
                                    [posX, posY, zmin], [posX, posY, zmax])
                texte += NL
    #
    texte += "# Sondes axiales au centre avec decalage pour Vy et Vz #" + NL
    epsilon = deltaEspace / 2
    nomSonde = "  S_Xc_CHI"
    texte += sondeSegment(nomSonde, "indicatrice", dt, Nx,
                          [xmin, Ly / 2, Lz / 2], [xmax, Ly / 2, Lz / 2])
    texte += NL
    nomSonde = "  S_Xc_VX"
    texte += sondeSegment(nomSonde, "velocity_x", dt, Nx,
                          [xmin, Ly / 2, Lz / 2], [xmax, Ly / 2, Lz / 2])
    texte += NL
    nomSonde = "  S_Xc_VYm"
    texte += sondeSegment(nomSonde, "velocity_y", dt, Nx,
                          [xmin, Ly / 2 - epsilon, Lz / 2], [xmax, Ly / 2 - epsilon, Lz / 2])
    texte += NL
    nomSonde = "  S_Xc_VYp"
    texte += sondeSegment(nomSonde, "velocity_y", dt, Nx,
                          [xmin, Ly / 2 + epsilon, Lz / 2], [xmax, Ly / 2 + epsilon, Lz / 2])
    texte += NL
    nomSonde = "  S_Xc_VZm"
    texte += sondeSegment(nomSonde, "velocity_z", dt, Nx,
                          [xmin, Ly / 2, Lz / 2 - epsilon], [xmax, Ly / 2, Lz / 2 - epsilon])
    texte += NL
    nomSonde = "  S_Xc_VZp"
    texte += sondeSegment(nomSonde, "velocity_z", dt, Nx,
                          [xmin, Ly / 2, Lz / 2 + epsilon], [xmax, Ly / 2, Lz / 2 + epsilon])
    texte += NL
    texte += "}"
    texte += NL

    with open(nomfic, "w") as fic:
        fic.write(texte)

    print("See output in : ", nomfic)
