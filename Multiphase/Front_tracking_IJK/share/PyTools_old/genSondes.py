# -*- coding : utf8
import readline, rlcompleter
readline.parse_and_bind('tab:complete')
import DNSTools as dtool

def sciform(x):
    return '{:.6E}'.format(x)

def sondePoint(nom, var, dt, coord):
    '''Definit une sonde point pour un fichier de donnees Trio'''

    sep = " "
    text = nom + sep + var + sep + str(dt) + sep + "points" +\
        sep + 1 + str(sciform(coord[0])) + sep + str(sciform(coord[1])) +\
        sep + str(sciform(coord[2]))
    return text


def sondeSegment(nom, var, dt, npoints, coordDebut, coordFin):
    '''Definit une sonde segment pour un fichier de donnees Trio'''
    sep = " "
    text = nom + sep + var + sep + "periode" + sep + str(dt) +\
        sep + "Segment" + sep + str(npoints) + sep +\
        str(sciform(coordDebut[0])) + sep + str(sciform(coordDebut[1])) + sep +\
        str(sciform(coordDebut[2])) + sep + str(sciform(coordFin[0])) + sep +\
        str(sciform(coordFin[1])) + sep + str(sciform(coordFin[2]))
    return text

if (len(sys.argv) != 2):
   raise Exception("usage: python %s jdd.data\n Datafile is required to read the domain.\nOutput is in /tmp/sondes.txt"%sys.argv[0])

jdd=sys.argv[1]
print "Reading domain information from: ", jdd
# Donnees utilisateur
Lx=dtool.getParam(jdd, 'uniform_domain_size_i')
Ly=dtool.getParam(jdd, 'uniform_domain_size_j')
Lz=dtool.getParam(jdd, 'uniform_domain_size_k')
Nx = (int)(dtool.getParam(jdd, 'nbelem_i'))
Ny = (int)(dtool.getParam(jdd, 'nbelem_j'))
Nz = (int)(dtool.getParam(jdd, 'nbelem_k'))
###
# Lx = 0.02
# Ly = 0.005
# Lz = 0.005
# Nx = 512
# Ny = 128
# Nz = 128
deltaEspace = Lx/Nx
dt = 1e-6  # pas de temps d'ecriture de la sonde
# Sondes dans des plans yOz. a NsondeX[0] hauteures differentes, on
# place NsondeY[1] sondes selon Oy et NsondeZ[2] selon Oz
NsondeX = [9, 4, 4]
# Sondes selon l'axe Ox. On place des segments sur NsondeY points en
# Oy et NsondeZ en Oz
NsondeY = 7
NsondeZ = 7

####################
# Script
dx = [Lx/(NsondeX[0]+1), Ly/(NsondeX[1]+1), Lz/(NsondeX[2]+1)]
dy = Ly/(NsondeY+1)
dz = Lz/(NsondeZ+1)
xmin = 0.0
xmax = Lx
ymin = 0.0
ymax = Ly
zmin = 0
zmax = Lz

NL = "\n  "
texte = "  Sondes"+NL
texte += "{"+NL
texte += "# Sondes axiales #"+NL
#
for i in range(1, NsondeY+1):
    for j in range(1, NsondeZ+1):
        nomSonde = "  S_X"+str(i)+str(j)+"_VX"
        texte += sondeSegment(nomSonde, "velocity_x", dt, Nx,
                              [xmin, i*dy, j*dz], [xmax, i*dy, j*dz])
        texte += NL
        nomSonde = "  S_X"+str(i)+str(j)+"_VY"
        texte += sondeSegment(nomSonde, "velocity_y", dt, Nx,
                              [xmin, i*dy, j*dz], [xmax, i*dy, j*dz])
        texte += NL
        nomSonde = "  S_X"+str(i)+str(j)+"_VZ"
        texte += sondeSegment(nomSonde, "velocity_z", dt, Nx,
                              [xmin, i*dy, j*dz], [xmax, i*dy, j*dz])
        texte += NL
        nomSonde = "  S_X"+str(i)+str(j)+"_CHI"
        texte += sondeSegment(nomSonde, "indicatrice", dt, Nx,
                              [xmin, i*dy, j*dz], [xmax, i*dy, j*dz])
        texte += NL
#
texte += "# Sondes transverses, axe Oy #"+NL
for i in range(1, NsondeX[0]+1):
    for j in range(1, NsondeX[1]+1):
        nomSonde = "  S_Y"+str(i)+str(j)+"_VX"
        texte += sondeSegment(nomSonde, "velocity_x", dt, Ny,
                              [i*dx[0], ymin, j*dx[2]], [i*dx[0], ymax, j*dx[2]])
        texte += NL
        nomSonde = "  S_Y"+str(i)+str(j)+"_VY"
        texte += sondeSegment(nomSonde, "velocity_y", dt, Ny,
                              [i*dx[0], ymin, j*dx[2]], [i*dx[0], ymax, j*dx[2]])
        texte += NL
        nomSonde = "  S_Y"+str(i)+str(j)+"_VZ"
        texte += sondeSegment(nomSonde, "velocity_z", dt, Ny,
                              [i*dx[0], ymin, j*dx[2]], [i*dx[0], ymax, j*dx[2]])
        texte += NL
        nomSonde = "  S_Y"+str(i)+str(j)+"_CHI"
        texte += sondeSegment(nomSonde, "indicatrice", dt, Ny,
                              [i*dx[0], ymin, j*dx[2]], [i*dx[0], ymax, j*dx[2]])
        texte += NL
#
texte += "# Sondes transverses, axe Oz #"+NL
for i in range(1, NsondeX[0]+1):
    for j in range(1, NsondeX[2]+1):
        nomSonde = "  S_Z"+str(i)+str(j)+"_VX"
        texte += sondeSegment(nomSonde, "velocity_x", dt, Nz,
                              [i*dx[0], j*dx[1], zmin], [i*dx[0], j*dx[1], zmax])
        texte += NL
        nomSonde = "  S_Z"+str(i)+str(j)+"_VY"
        texte += sondeSegment(nomSonde, "velocity_y", dt, Nz,
                              [i*dx[0], j*dx[1], zmin], [i*dx[0], j*dx[1], zmax])
        texte += NL
        nomSonde = "  S_Z"+str(i)+str(j)+"_VZ"
        texte += sondeSegment(nomSonde, "velocity_z", dt, Nz,
                              [i*dx[0], j*dx[1], zmin], [i*dx[0], j*dx[1], zmax])
        texte += NL
        nomSonde = "  S_Z"+str(i)+str(j)+"_CHI"
        texte += sondeSegment(nomSonde, "indicatrice", dt, Nz,
                              [i*dx[0], j*dx[1], zmin], [i*dx[0], j*dx[1], zmax])
        texte += NL
#
texte += "# Sondes axiales au centre avec decalage pour Vy et Vz #"+NL
epsilon = deltaEspace/2
nomSonde = "  S_Xc_CHI"
texte += sondeSegment(nomSonde, "indicatrice", dt, Nx,
                      [xmin, Ly/2, Lz/2], [xmax, Ly/2, Lz/2])
texte += NL
nomSonde = "  S_Xc_VX"
texte += sondeSegment(nomSonde, "velocity_x", dt, Nx,
                      [xmin, Ly/2, Lz/2], [xmax, Ly/2, Lz/2])
texte += NL
nomSonde = "  S_Xc_VYm"
texte += sondeSegment(nomSonde, "velocity_y", dt, Nx,
                      [xmin, Ly/2-epsilon, Lz/2], [xmax, Ly/2-epsilon, Lz/2])
texte += NL
nomSonde = "  S_Xc_VYp"
texte += sondeSegment(nomSonde, "velocity_y", dt, Nx,
                      [xmin, Ly/2+epsilon, Lz/2], [xmax, Ly/2+epsilon, Lz/2])
texte += NL
nomSonde = "  S_Xc_VZm"
texte += sondeSegment(nomSonde, "velocity_z", dt, Nx,
                      [xmin, Ly/2, Lz/2-epsilon], [xmax, Ly/2, Lz/2-epsilon])
texte += NL
nomSonde = "  S_Xc_VZp"
texte += sondeSegment(nomSonde, "velocity_z", dt, Nx,
                      [xmin, Ly/2, Lz/2+epsilon], [xmax, Ly/2, Lz/2+epsilon])
texte += NL
texte += "}"
texte += NL

with open("/tmp/sondes.txt", "w") as fic:
    fic.write(texte)

print "See output in : ", "/tmp/sondes.txt"
