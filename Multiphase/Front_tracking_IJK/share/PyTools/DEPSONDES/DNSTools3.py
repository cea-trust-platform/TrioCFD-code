# -*- coding: utf-8 -*-
#
# Boite a outils pour les post-traitement du IJK...
# 
#
import numpy as np

from math import *
import glob, os, re

alwidth = 0.75
mewidth = 0.75
errlw=0.25
lwidth = 1.5
pltmewidth = 0
msize = 7.5
#### Utilisation de LaTeX

# try:
#     rc('text', usetex = True)
#     #rcParams['text.latex.preamble']=[r"\usepackage{gensymb}",r"\usepackage{txfonts}",r"\usepackage{nicefrac}",r"\usepackage{amssymb}",r"\usepackage{sistyle}"]
#     rcParams['text.latex.preamble']=[r"\usepackage{txfonts}",r"\usepackage{amssymb}"]
#     rc('legend', fontsize='medium',numpoints=2)
# except:

def cd(fold):
    os.chdir(fold)
    return
    
def pwd():
    return os.getcwd()
    
def getJddName():
    l = glob.glob("*data")
    if len(l)>1:
        l  = [s for s in l if "repr" not in s]
    jdd = l[0][:-5]
    return jdd


# Recuperation des cdg des bulles : 
# t : Le temps
# nb : Le nombre de bulles
# coords[:,it,ib] : Les coords du cdg de la bulle ib a l'instant it
def getBarys(head=None):
    if head == None:
        head = getJddName()
    fics = glob.glob(head.replace(".data","")+"_bulles_centre_[xyz].out")
    if (len(fics)!=3):
        raise Exception("Error in getBarys: Missing files :"+fics)
    for fic in fics:
        ax = fic[-5]
        exec("%s = loadtxt(fic)" %ax)
    #
    t = x[:,0]
    x = x[:,1:]
    y = y[:,1:]
    z = z[:,1:]
    print(len(x), shape(x)[0])
    lmin=min(min(shape(x)[0],shape(y)[0]),shape(z)[0])#-1
    print(lmin)
    coords = np.array([x[:lmin,:], y[:lmin,:], z[:lmin,:]])
    _, nb = shape(x)
    #print(coords)
    return coords, t[:lmin], nb

def getDimensions(jdd=None):
    if jdd == None:
        jdd = getJddName()+".data"
    Lx = getValue("uniform_domain_size_i", jdd) 
    Ly = getValue("uniform_domain_size_j", jdd) 
    Lz = getValue("uniform_domain_size_k", jdd) 
    return np.array([Lx, Ly, Lz])

def getVelocity(coords, time, nb, jdd=None):
    if jdd == None:
        jdd = getJddName()+".data"
    Lx = getValue("uniform_domain_size_i", jdd) 
    Ly = getValue("uniform_domain_size_j", jdd) 
    Lz = getValue("uniform_domain_size_k", jdd) 
    L = np.array([Lx, Ly, Lz])
    dt = time[1:]-time[:-1]
    tvel = (time[1:]+time[:-1])/2.
    ndim, nt, nb = shape(coords)
    nt-=1
    velocity = np.zeros((ndim, nt, nb))
    for direction in [0, 1, 2]:
        LL = L[direction]
        dx = coords[direction,1:,:]-coords[direction,:-1,:]
        # Correction du dx si traversee de frontiere perio : 
        nt, nb = shape(dx)
        for i in range(nt):
            for j in range(nb):
                x = dx[i,j] 
                if (abs(x)>LL*0.5):
                    print("Crossing periodic boundary ", direction, " at iter=", i, " at t=", time[i], )
                    print("Old ", x, )
                    # Un deplacement superieur a la moitie de la taille du domaine, 
                    # C'est forcement un passage par une frontiere perio :
                    while x - LL*0.5> 0.:
                        x-= LL
                    while x + LL*0.5< 0.:
                        x+= LL
                    print(" New ", x)
                dx[i,j] = x
        # dx = [ sign(x)*(x % L[direction]) for x in dx ]
        velocity[direction, :,:] = dx/dt[:,newaxis]
    return tvel, velocity
# 
# Retourne : temps, valeurs, nombre de bulles
# val[it, ib] : valeur de var au temps it pour la bulle ib
# time, val, nbulles = getOnBulles("centre_x")
def getOnBulles(var, head=None):
    if head == None:
        head = getJddName()
    mat = loadtxt(head+"_bulles_"+var+".out")
    time = mat[:,0]
    val = mat[:,1:]
    _, nbulles = shape(val)
    return time, val, nbulles

def getTempsIntegrationFile(fic):
    f = open(fic,"r")
    lines = f.readlines()
    tintegration = lines[0].split()[2]
    f.close()
    return float(tintegration)

def getTempsIntegrationList(fics):
    lt = []
    for fic in fics:
        lt.append(getTempsIntegrationFile(fic))
    return np.array(lt)

def buildDicoColonnes(fic):
    d = {}
    f = open(fic,"r")
    lines = f.readlines()
    for val,key in [(int(s.split()[2])-1,s.split()[4])  for s in lines if "# colonne" in s]:
        d[key] = val
    f.close()
    return d

def getStatsInstant():
    return getStats("ins")

def getStatsMoy():
    return getStats()


# caract = moy or ins : instantané ou moyen. 
# Retourne une matrice avec les stats :
# resu[i,j,k] : i --> temps
#                    j --> Position en z
#                    k --> var.
#
# ltimes : liste des temps.
# lz      : liste des coords en z.
# lvar    : liste des variables stockees.
# val[it, iz, ivar] : La matrice 3D contenant tout ca. 
#
# tintegration : le temps d'integration des resulats si caract == moy.
def getStats(caract = "moy"):
    if caract == "moy":
        fic_stats = glob.glob("statistiques_*.txt")
        if (fic_stats == []): 
            fic_stats = glob.glob("*phasique_statistiques_*.txt")
        fic_stats.sort(key=lambda item: float(item.strip("monodiphasique_statistiques_.txt")))
        ltimes = np.array([float(f.strip("monodiphasique_statistiques_.txt")) for f in fic_stats])
    else:
        fic_stats = glob.glob("moyenne_spatiale_*.txt")
        if (fic_stats == []): 
            fic_stats = glob.glob("*phasique_moyenne_spatiale_*.txt")
        fic_stats.sort(key=lambda item: float(item.strip("monodiphasique_moyenne_spatiale_.txt")))
        ltimes = np.array([float(f.strip("monodiphasique_moyenne_spatiale_.txt")) for f in fic_stats])
    if (fic_stats == []): 
        raise Exception("La liste fic_stats est restee vide... ")
    lvar = buildDicoColonnes(fic_stats[0])
    nt = len(ltimes)
    if caract == "moy":
        tintegration = getTempsIntegrationList(fic_stats)
    else:
        tintegration = [0]
    for i, fic in enumerate(fic_stats):
        mat = loadtxt(fic)
        t = ltimes[i]
        lz = mat[:,0]
        if i == 0 :
            # Initialise la taille du resu : 
            nz, nvar = shape(mat)
            resu = np.zeros((nt, nz, nvar))
        resu[i, :,:] = mat[:,:]
    return ltimes, lz, lvar, resu, tintegration

def getValue(key, fic, pre="^", default=None):
    f = open(fic, 'r')
    lines = f.readlines()
    f.close()
    nb=len(lines)
    rc=re.compile(pre+"\s*"+key+"\s*(?P<value>[\-]?[\d]*.?[\d]*[eE]?[+\-]?[\d]*)")
    for i, st in enumerate(lines):
        m=rc.match(st)
        if m:
            val = float(m.group("value"))
            return val
            break
    if (default != None): return default
    print("On a rien trouve dans ", fic)
    raise Exception("Etonnant, non?")
    return -1


def getValues(key, fic, pre=""):
    f = open(fic, 'r')
    lines = f.readlines()
    f.close()
    nb=len(lines)
    rc=re.compile(pre+key+"\s*(?P<value>[\-]?[\d]*.?[\d]*[eE]?[+\-]?[\d]*)")
    for i, st in enumerate(lines):
        it = re.finditer(rc,st)
        val = []
        for match in it:
            val.append(float(match.group("value")))
        if len(val):
            return np.array(val)
    print("On a rien trouve dans ", fic)
    raise Exception("Etonnant, non?")
    return np.array([])

def getSondesCoords(fic):
    x=getValues("x=", fic)
    y=getValues("y=", fic)
    z=getValues("z=", fic)
    return np.array([x,y,z]).T, len(x)

# En z+ et En z- : 
def evaluateRetau(U, z, jdd=None):
    if jdd == None:
        jdd = getJddName()+".data"
    rhol = getValue("rho_liquide", jdd)
    mul = getValue("mu_liquide", jdd)
    Lz = getValue("uniform_domain_size_k", jdd) 
    h =  Lz / 2.
    tauwp = mul*U[:,0]/z[0]
    tauwm = mul*U[:,-1]/(Lz-z[-1])
    #
    tauw = (tauwp+tauwm)/2.
    utau = np.array([sqrt(abs(x)/rhol) for x in tauw])
    Retau = rhol*utau*h/mul
    return Retau, tauwp, tauwm

def getFloatingAverage(tintegration, Um, n):
    if (n>=len(tintegration)):
        print("Too long averaging ", n, " required!")
        return tintegration.mean(), Um.mean() # Whatever, in the good range
    tf = tintegration[n:]  # temps de fin
    ti = tintegration[:-n] # temps de debut
    #print("getFloatingAverage:: Attention : il faut prendre ltintegration et pas tm")
    Uf = Um[n:]
    Ui = Um[:-n]
    dt = tf-ti
    t = (tf+ti)/2.
    var = (Uf*tf-Ui*ti)/dt
    return t, var

def getFloatBtwIters(tab, ltintegration, it,it2):
    dt = ltintegration[it2]-ltintegration[it]
    res = (ltintegration[it2]*tab[it2,:] - ltintegration[it]*tab[it,:])/dt
    return res




# v : NumPy np.array
# dz : pas du maillage. On le suppose constant et on suppose que les bords sont a dz/2 du premier et dernier noeuds
#                              (ce qui correspond a un cell_center regulier). 
# bc_type = 0 : Pour la vitesse, suppose une valeur nulle sur les bords. neuman_homogene
#              2 : Suppose un gradient normal nul a la paroi, par exemple pour la pression
#              1 : Ne suppose rien, decentre interieur... (ordre 1)
def cell_to_cell_gradient(v,dz,bc_type):
     # if ( (z[1:]-z[:-1]).max() -  (z[1:]-z[:-1]).min() ) > 1e-11): raise Exception("Fonction pour pas uniforme!")
     # Formule centree (ordre 2) pour pas cste dans le domaine : 
     grad = np.zeros(len(v))
     grad[1:-1] = (v[2:]-v[:-2])/(2.*dz)
     # Aux bords :
     # Au bord en z- :
     Uc = v[0]
     Up = v[1]
     if (bc_type == 0):
          Um = 0.    
          grad[0]  = (- 4*Um + 3*Uc +Up  ) / (3*dz)
     elif (bc_type == 1):
          grad[0]  = (Up - Uc) / dz
     else:
          Um = Uc+(Uc-Up)/8.
          grad[0]  = (- 4*Um + 3*Uc +Up  ) / (3*dz)
     #
     # Au bord en z+ :
     Uc = v[-1]
     Um = v[-2]
     if (bc_type == 0):
          Up = 0.
          grad[-1] = (- Um - 3*Uc +4*Up ) / (3*dz)
     elif (bc_type == 1):
          grad[-1]  = (Uc - Um) / dz
     else:
          Up = Uc+(Uc-Um)/8.
          grad[-1] = (- Um - 3*Uc +4*Up ) / (3*dz)
     #
     return grad     
     
def deriv(v,z, Lz):
     raise Exception("Tentative pour maillage a pas variable incorrecte sur les bords!")
     xf = (z[1:] + z[:-1])/2. # Centre des faces
     xf = np.append(np.append(0.,xf),Lz)
     dx = xf[1:] - xf[:-1] # dx_i = x_i+1 - x_i
     # L'evaluation du dx entre 2 cellules :
     h = (dx[:-1] + dx[1:]) / 2.# contient (Nsommets - 2) valeurs = Ncellules - 1
     h1 = h[:-1]
     h2 = h[1:]
     u_pl = v[2:]
     u_m = v[:-2]
     u_c = v[1:-1]
     # Formule centree (ordre 2) pour pas variable dans le domaine : 
     grad = np.zeros(len(v))
     grad[1:-1] = (h1/h2*u_pl - h2/h1*u_m + (h2**2-h1**2)/(h1*h2)*u_c) / (h1+h2)
     # Formule decentree (ordre 2) pour pas variable sur le bord gauche :
     h1 = h[0] 
     h2 = h[1]
     u_g    = v[0]
     u_pl  = v[1]
     u_pl2 = v[2]
     grad[0] = (1/(h1*h2*(h1+h2))) * (-(2*h1*h2+h2**2)*u_g + ((h1+h2)**2)*u_pl - (h1**2)*u_pl2)
     # Formule decentree (ordre 2) pour pas variable sur le bord droit : 
     h1 = h[-2] 
     h2 = h[-1]
     u_d  = v[-1]
     u_m  = v[-2]
     u_m2 = v[-3]
     grad[-1] = 1/(h1*h2*(h1+h2)) * ((h2**2)*u_m2 - ((h1+h2)**2)*u_m - ((h2**2)-(h1+h2)**2)*u_d)
     return grad


# v : NumPy np.array
# dz : pas du maillage. On le suppose constant et on suppose que les bords sont a dz/2 du premier et dernier noeuds
#                              (ce qui correspond a un cell_center regulier). 
# bc_type = 0 : Pour la vitesse, suppose une valeur nulle sur les bords. neuman_homogene
#              2 : Suppose un gradient normal nul a la paroi, par exemple pour la pression
#              1 : Ne suppose rien, decentre interieur... (ordre 1)
def cell_to_cell_second_gradient(v,dz,bc_type):
     # Formule centree (ordre 2) pour pas cste dans le domaine :
     dz2 = dz*dz
     grad2 = np.zeros(len(v))
     # L'erreur est d'ordre 2
     grad2[1:-1] = (v[:-2] - 2*v[1:-1] + v[2:]) / dz2
     # Aux bords :
     # Au bord en z- :
     Uc = v[0]
     Up = v[1]
     if (bc_type == 0):
          Um = 0.    
          # L'erreur est d'ordre 1
          grad2[0]  = 4./3. * (2*Um - 3*Uc + Up) / dz2
     elif (bc_type == 1):
          Upp = v[2]
          # L'erreur est d'ordre 1
          grad2[0]  = (-2*Up + Uc + Upp) / (2*dz2)
     else:
          Um = Uc+(Uc-Up)/8.
          # L'erreur est d'ordre 1, mais plus grande qu'avec bc_type=0
          grad2[0]  = 4./3. * (2*Um - 3*Uc + Up) / dz2
     #
     # Au bord en z+ :
     Uc = v[-1]
     Um = v[-2]
     if (bc_type == 0):
          Up = 0.
          grad2[-1]  = 4./3. * (Um - 3*Uc + 2*Up) / dz2
     elif (bc_type == 1):
          Umm = v[-3]
          grad2[-1]  = (-2*Um + Uc + Umm) / (2*dz2)
     else:
          Up = Uc+(Uc-Um)/8.
          grad2[-1]  = 4./3. * (Um - 3*Uc + 2*Up) / dz2
     #
     return grad2 

     
def test_cell_to_cell_gradient(z,Lz):
     # Test d'une fonction nulle aux bords :
     u = np.array([sin(2.*pi*zz/Lz) for zz in z])
     nz = len(u) 
     dz = Lz/nz
     dudz = cell_to_cell_gradient(u,dz,0) # bc_type = champ nul au bord
     dudz_ana = np.array([2.*pi/Lz*cos(2.*pi*zz/Lz) for zz in z])
     err_max = np.fabs(dudz - dudz_ana).max()
     tol = 2*pi*pi/(nz*nz*dz)*2*pi/nz
     if (err_max > tol):
          raise Exception("Failing function u test_deriv error = %g > %g = tol"%(err_max,tol))
     print("Error max u : ", err_max, "(tol=", tol, ")")
     #
     # Test d'une fonction de derivee normale nulle aux bords :
     w = np.array([cos(2.*pi*zz/Lz)+10 for zz in z])
     nz = len(w) 
     dz = Lz/nz
     dwdz = cell_to_cell_gradient(w,dz,2) # bc_type = derivee normale du champ nulle au bord
     dwdz_ana = np.array([-2.*pi/Lz*sin(2.*pi*zz/Lz) for zz in z])
     err_max = np.fabs(dwdz - dwdz_ana).max()
     tol = 2*pi*pi/(nz*nz*dz)*2*pi/nz
     if (err_max > tol):
          raise Exception("Failing function w test_deriv error = %g > %g = tol"%(err_max,tol))
     print("Error max w : ", err_max, "(tol=", tol, ")")
     #
     # Test d'une fonction sans propriete aux bords : 
     v = np.array([sin(2.*pi*zz/Lz)+10 for zz in z])
     nz = len(v) 
     dz = Lz/nz
     dvdz = cell_to_cell_gradient(v,dz,1) # bc_type = quelconque
     dvdz_ana = np.array([2.*pi/Lz*cos(2.*pi*zz/Lz) for zz in z])
     err_max = np.fabs(dvdz - dvdz_ana).max()
     tol = 2*pi*pi/(nz*nz*dz)
     if (err_max > tol):
          raise Exception("Failing function v test_deriv error = %g > %g = tol"%(err_max,tol))
     print("Error max v : ", err_max, "(tol=", tol, ")")
     #
     return True
     
def test_cell_to_cell_second_gradient(z,Lz):
     # Test d'une fonction nulle aux bords :
     u = np.array([sin(2.*pi*zz/Lz)+zz*(Lz-zz)*1000. for zz in z])
     nz = len(u) 
     dz = Lz/nz
     dduddz = cell_to_cell_second_gradient(u,dz,0) # bc_type = champ nul au bord
     dduddz_ana = np.array([-(2.*pi/Lz)**2*sin(2.*pi*zz/Lz)-2000. for zz in z])
     err_max = np.fabs(dduddz - dduddz_ana).max()
     tol = 1./6.*dz*(2*pi/Lz)**3
     if (err_max > tol):
          raise Exception("Failing function u test_second_gradient error = %g > %g = tol"%(err_max,tol))
     print("Error max u : ", err_max, "(tol=", tol, ")")
     #
     # Test d'une fonction de derivee normale nulle aux bords :
     w = np.array([cos(2.*pi*zz/Lz)+1000 for zz in z])
     nz = len(w) 
     dz = Lz/nz
     ddwddz = cell_to_cell_second_gradient(w,dz,2) # bc_type = derivee normale du champ nulle au bord
     ddwddz_ana = np.array([-(2.*pi/Lz)**2*cos(2.*pi*zz/Lz) for zz in z])
     err_max = np.fabs(ddwddz - ddwddz_ana).max()
     tol = (1./6.-1./8.)*dz*(2*pi/Lz)**3
     if (err_max > tol):
          raise Exception("Failing function w test_second_gradient error = %g > %g = tol"%(err_max,tol))
     print("Error max w : ", err_max, "(tol=", tol, ")")
     #
     # Test d'une fonction sans propriete aux bords : 
     v = np.array([sin(2.*pi*zz/Lz)+10+1000.*zz for zz in z])
     nz = len(v) 
     dz = Lz/nz
     ddvddz = cell_to_cell_second_gradient(v,dz,1) # bc_type = quelconque
     ddvddz_ana = np.array([-(2.*pi/Lz)**2*sin(2.*pi*zz/Lz) for zz in z])
     err_max = np.fabs(ddvddz - ddvddz_ana).max()
     tol = 1./2.*dz*(2*pi/Lz)**3
     if (err_max > tol):
          raise Exception("Failing function v test_second_gradient error = %g > %g = tol"%(err_max,tol))
     print("Error max v : ", err_max, "(tol=", tol, ")")
     #
     return True


# v : NumPy np.array
# Derivee seconde (Laplacien)
# On suppose un pas uniforme
def deriv2(v,z, Lz):
     raise Exception("Tentative de derivee seconde pour maillage a pas variable non validee!")
     xf = (z[1:] + z[:-1])/2. # Centre des faces
     xf = np.append(np.append(0.,xf),Lz)
     dx = xf[1:] - xf[:-1] # dx_i = x_i+1 - x_i
     # L'evaluation du dx entre 2 cellules :
     h = (dx[:-1] + dx[1:]) / 2.# contient (Nsommets - 2) valeurs = Ncellules - 1
     if (h.max()-h.min()> 0.01*h.min()) : raise Exception("Bad use of deriv2 : varying step size!")
     #
     # On prend le h moyen 
     h1 = h.mean()
     u_pl = v[2:]
     u_m = v[:-2]
     u_c = v[1:-1]
     # Formule centree (ordre 2) pour pas cste dans le domaine : 
     laplacien = np.zeros(len(v))
     laplacien[1:-1] = (u_pl - 2*u_c + u_m) / (h1*h1)
     # Formule decentree (ordre 1) pour pas cste sur le bord gauche :
     u_g    = v[0]
     u_pl  = v[1]
     u_pl2 = v[2]
     laplacien[0] = (u_g - 2*u_pl + u_pl2) / (h1*h1)
     # Formule decentree (ordre 1) pour pas cste sur le bord droit : 
     u_d  = v[-1]
     u_m  = v[-2]
     u_m2 = v[-3]
     laplacien[-1] =  (u_d - 2*u_m + u_m2) / (h1*h1)
     return laplacien

def integrate(vec,dz):
    nz = len(vec)
    integral = np.zeros(nz)
    old = 0
    for i,val in enumerate(vec):
        old += val*dz
        integral[i] = old
    return integral

def getParam(fic, name, vec=False, compo=0, string=False):
    f=open(fic)
    lines = f.readlines()
    f.close()
    for i, line in enumerate(lines):
        if line.find("#")>=0:
            lp = line.find("#")
            rp = line.rfind("#")
            if rp>lp:
                s = line[:lp-1]+line[rp+1:]
                # print(lines[i], "replaced by ", s)
                lines[i]= s
    goodline=[line for line in lines if name in line.split()]
    #print(goodline)
    # Suppression des lignes contenant un commentaire : 
    # goodline=[line for line in goodline if "#" not in line]
    # print(goodline)
    if len(goodline)==1:
        if vec:
            val=goodline[0].split()[2+compo]
        else:
            val=goodline[0].split()[1]
        try:
            return float(val)
        except:
            if (not string): raise Exception('Trying to return a string? %s'%val)
            return val
    else:
        raise Exception('ohoh, looking for name=%s commented line: %s'%(name,goodline))
    return 0/0

def get_prop(jdd='DNS.data', repr_file=None, Ret=0.):
    Lx=getParam(jdd, 'uniform_domain_size_i')
    Ly=getParam(jdd, 'uniform_domain_size_j')
    Lz=getParam(jdd, 'uniform_domain_size_k')
    mul=getParam(jdd, 'mu_liquide')
    rhol=getParam(jdd, 'rho_liquide')
    try:
        muv=getParam(jdd, 'mu_vapeur')
        rhov=getParam(jdd, 'rho_vapeur')
    except:
        muv=mul
        rhov=rhol
    #
    print("######## Dnstool3")
    print(jdd)
    print(repr_file)
    print(getParam(jdd, 'vol_bulle_monodisperse'))
    print(getParam(repr_file, 'bubble_groups'))
    print(getParam(jdd, 'sigma'))
    try:
        vb=getParam(jdd, 'vol_bulle_monodisperse')
        nb=getParam(repr_file, 'bubble_groups')
        sigma=getParam(jdd, 'sigma')
    except:
        print("Bubble assumed volume 0")
        vb=0.;nb=0.; sigma=1.
    #
    g=getParam(jdd, 'gravite', vec=True,compo=0)
    #
    Db=pow(vb*6./pi, 1./3.)
    alv=nb*vb/(Lx*Ly*Lz)
    h=Lz/2.
    rhom=alv*rhov+(1-alv)*rhol
    utau=Ret*mul/(rhol*h) 
    tauw=rhol*utau*utau 
    Eo=rhol*abs(g)*Db**2/sigma
    beta=tauw/h-rhom*g  
    try:
        str_Svx=getParam(jdd, 'expression_variable_source_x')
        print("Variable source read as: ", str_Svx)
        Svx=float(str_Svx)
    except:
        # The default value is 0. 
        Svx=0.
    #
    verb=1
    if verb:
        print("rhol = ", rhol)
        print("rhov = ", rhov)
        print("mul = ", mul)
        print("muv = ", muv)
        print("vb = ", vb)
        print("nb = ", nb)
        print("sigma = ", sigma)
        print("g = ", g)
        print("alv = ", alv)
        print("Ret = ", Ret)
        print("rhom = ", rhom)
        print("beta = ", beta)
        print("Svx = ", Svx)
        print("tauw = ", tauw)
        print("Db = ", Db)
        print("Eo = ", Eo)
    #
    return rhol, rhov, mul, muv, alv, beta, sigma, Lz, rhom, Eo, g, Svx



class Simu:
    def __init__(self, jdd='DNS.data', repr_file=None, Ret=0.):
        rhol, rhov, mul, muv, alv, beta, sigma, Lz, rhom, Eo, g, Svx = get_prop(jdd, repr_file, Ret)
        self.rhol = rhol
        self.rhov = rhov
        self.rhom = rhom
        self.mul = mul
        self.muv = muv
        self.alv = alv
        self.beta = beta
        self.sigma = sigma
        self.Lz = Lz
        self.Eo = Eo
        self.g = g
        self.Svx = Svx
    




