# -*- coding: utf-8 -*-
from optparse import OptionParser
import MEDCouplingRemapper # Pour l'interpolateur
import numpy as np
import time # Pour la création d'un nouveau maillage irrégulier
import MEDLoader as ml
import methodes_2020 # comme ca on est certain que si on prend juste une fonction
                     # qui en appelait d'autres de ce mm fichier, elle ocntinuera de ofnctionner.
import outils_physique
import matplotlib.pyplot as plt


# ######################################################################
# ################### README ###########################################
# Scénario d'utilisation:
# 1. On dispose d'un fichier .med, on y applique MED2NumPy pour transformer
#     le champ voulu en NumPy array
# 2. On doit pouvoir appliquer à ce NumPy TOUTES les fonctions définies
#    ci-après. (en cours)
# 3. Si on veut ressortir un fichier .med à partir du Numpy, on se servira
#    EXCLUSIVEMENT de la méthode NumPy2MED.
#
# ATTENTION : Sur les conventions
#     - Convention NumPy : Nx, Ny, Nz, Ncompo1, Ncompo2, Ntemps
#          --> Ncompo2 existe uniquement si on considère le gradient du
#              vecteur vitesse, ou si on considère les tensions de Reynolds
#     - Convetion MED    : Nz, Ny, Nx, Ncompo
#          --> Si Ncompo2 existe, alors Ncompo = Ncompo1*Ncompo2.
#              Sinon,                   Ncompo = Ncompo1
#          --> On accede au temps différement
# ######################################################################
# ######################################################################

def en_tete(MEDchamp):
    #################################################################################
    # Passe d'un MEDchamp sur UMesh à MEDchamp sur CMesh
    #################################################################################
    if not MEDchamp.getMesh().isStructured():
     mC, idx = fromUtoC(MEDchamp.getMesh())
     f2 = UMeshToCMesh(MEDchamp, mC, idx)
    elif MEDchamp.getMesh().isStructured():
     mC = MEDchamp.getMesh()
     f2 = MEDchamp.deepCopy()
     pass
    return (f2,mC)

def MEDFileToMEDField_MultiTS(MEDFileName, MEDFieldName):
    #################################################################################
    # Passe d'un fichier MED à un champ MEDFileField, défini sur plsr pas de temps
    # On donne aussi le mesh en MEDFileMesh...
    # --> Pour récuperer un MEDField : Ffield.getFieldOnMeshAtLevel(...,Fmesh,...)
    # --> Pour avoir les iterations disponibles : Ffield.getIterations(...)
    # --> Temps associé à tel pas de temps : ml.GetTimeAttachedOnFieldIteration(...)
    # --> Pour récuperer un MEDMesh  : Fmesh.getMeshAtLevel(...)
    #################################################################################
    Ffield = ml.MEDFileFieldMultiTS.New(MEDFileName,MEDFieldName,True)

    # it = field.getIterations()

    # On veut récupérer la discrétisation spatiale du champ
    # TypeOfField : 0-->ml.ON_CELLS, 1-->ml.ON_NODES, ... 4-->ml.ON_NODES_KR
    if "DOM_dual" in MEDFieldName:
        meshName = "DOM_dual"
    elif (("DOM" in MEDFieldName) and not("_dual" in MEDFieldName)):
        meshName = "DOM"
    elif "Field" in MEDFieldName:
        meshName = "Mesh"
    elif "velocity" in MEDFieldName:
        meshName = "Mesh"
        """ a completer avec les autres maillages possibles (DOM_Ext), ... """
        pass
        
    Fmesh = ml.MEDFileMesh.New(MEDFileName,meshName)
    TypeOfField = ml.GetTypesOfField(MEDFileName,meshName,MEDFieldName)[0]
    print('Detected meshName : ',meshName)
    # mesh = Fmesh.getMeshAtLevel(MeshLevel)
    # field = Ffield.blabla
    return (Ffield,Fmesh,TypeOfField)


def MEDFileToMEDField_1TS(MEDFileName, MEDFieldName):
    #################################################################################
    # Passe d'un fichier MED à un champ MED, défini sur 1 seul pas de temps
    # pour le moment on ne sert pas de cette fonction --> elle ne fonctionne
    # sans doute pas bien du coup
    #################################################################################

    mesh = ml.MEDFileMesh.New(MEDFileName)
    field = ml.MEDFileField1TS.New(MEDFileName,MEDFieldName,True)

    # On veut récupérer la discrétisation spatiale du champ
    # TypeOfField : 0-->ml.ON_CELLS, 1-->ml.ON_NODES, ... 4-->ml.ON_NODES_KR
    if "DOM_dual" in MEDFieldName:
        meshName = "DOM_dual"
    elif (("DOM" in MEDFieldName) and not("_dual" in MEDFieldName)):
        meshName = "DOM"
    elif "Field" in FieldName:
        meshName = "Mesh"
        """ a completer avec les autres maillages possible (DOM_Ext), ... """
        pass
    TypeOfField = ml.GetTypesOfField(MEDFileName,meshName,MEDFieldName)[0]
    
    return (field,mesh,TypeOfField)

def MEDFieldToNumpy(CMEDField,mC):
    #################################################################################
    # Passe d'un champ MED - AUQUEL ON A APPLIQUÉ en_tete - à un Numpy.
    # On renverra aussi les coordonnees du maillage, car c pratique quand meme.
    #################################################################################
    # print("11111111111",CMEDField.getNumberOfTuples())
    f = CMEDField.deepCopy()
    # print("2",CMEDField.getNumberOfTuples())
    # print("3",f.getNumberOfTuples())
    # print("4",f.getArray())
    # print("5",f.getArray().getNumberOfTuples())

    # Récupère les coordonnées des noeuds d'un maillage cartésien par direction : 
    xn = mC.getCoordsAt(0)
    yn = mC.getCoordsAt(1)
    zn = mC.getCoordsAt(2)
    # Nombre de mailles/elements : 
    Nx = len(xn) - 1
    Ny = len(yn) - 1
    Nz = len(zn) - 1
    Ncompo = f.getNumberOfComponents()
    print("Ncompo",Ncompo)
    # Internal storage of a Cfield's order : (1st spatial component varies the quickest):
    champ = np.array(f.getArray().getValues())

    print("MEDFieldToNumpy : ",champ.shape)
    print("f.getArray()")
    print(f.getNumberOfTuples())
    print("champ")
    print(champ)
    # SI ON A DES ERREURS, voir s'il faudrai pas plutot ces lignes la :
    # v = f.getArray()
    # vals = v.deepCopy()
    # vals.rearrange(1)  # rearrange() erase compo names, etc .. hence the deep copy above.
    # champ = vals.toNumPyArray()
    
    champ = champ.reshape(Nz, Ny, Nx, Ncompo)     # space coord on the last axis
    champ = np.transpose(champ, (2,1,0,3)) # Tourne compo pour avoir l'ordre : (Nx, Ny, Nz, Ncompo)    
    return (xn,yn,zn,champ)

def MEDFieldToNumpy_SansChamp(mC):
    #################################################################################
    # On renvoie les coordonnees du maillage, car c pratique quand meme,
    # a partir du champ MED auquel on a appliqué en_tete.
    #################################################################################

    # Récupère les coordonnées des noeuds d'un maillage cartésien par direction : 
    xn = mC.getCoordsAt(0)
    yn = mC.getCoordsAt(1)
    zn = mC.getCoordsAt(2)
   
    return (xn,yn,zn)

def NumpyToMEDField(champ,CMEDmesh,temps=[0,0,0],NomChamp="Mon_champ"):
    #################################################################################
    # Passage d'un champ Numpy avec son Cmesh, son temps et son nom à un
    # MEDField
    # Remarque : o pourrai se passer des 2 lignes à fillFromAnalytic
    #################################################################################

    sh = champ.shape
    Ndim = CMEDmesh.getSpaceDimension()

    MEDchamp = ml.MEDCouplingFieldDouble(ml.ON_CELLS) # a changer ca
    MEDchamp.setNature(ml.IntensiveMaximum)
    MEDchamp.setName(NomChamp)
    MEDchamp.setMesh(CMEDmesh)
    MEDchamp.setTime(temps[0],temps[1],temps[2])
    if len(sh) == Ndim+1:
        MEDNcompo = sh[-1]
        MEDchamp.fillFromAnalytic(sh[-1], "0")
        # sh[-1]=3 pour un champ vectoriel; sh[-1]=1 pour un champ scalaire
        arr = np.transpose(champ, (2,1,0,3))
    elif len(sh) == Ndim+2 :
        MEDNcompo = sh[-1]*sh[-2]
        MEDchamp.fillFromAnalytic(sh[-1]*sh[-2], "0")
        # sh[-1]*sh[-2] = 9 (3D) ou 4 (2D). On est ici pour Rij ou pour jacob
        arr = np.transpose(champ, (2,1,0,3,4))   # Nz, Ny, Nx, Ncomposante, NtermeParComposante
    else:
        MEDNcompo = 1
        print("NumpyToMEDField : Numpy has too many dimensions")
    arr = arr.flatten()
    arr = ml.DataArrayDouble(arr)
    arr.rearrange(MEDNcompo)
    MEDchamp.setArray(arr)
    return (MEDchamp)


def NumpyToMEDField_VRij(champ,MEDchamp,CMEDmesh,temps=[0,0,0]):
    #################################################################################
    # Pour un MEDchamp déja instancié (ml.MEDCouplingFieldDouble.New) a qui on a
    # déjà associé un maillage, cette fonction sert à associer temps, nature et array
    # (Utile quand on boucle sur le temps, pour ne pas recharger le mesh a chq fois)
    # Normalement on ne doit pas avoir de len(sh)==Ndim+2 dans Rij,
    # Mais on le garde si on souhaite utiliser la fct ailleurs...
    # La difference avec le precedent c'est qu'on ne set ni Mesh, ni Name.
    # On ne set que Nature (utile ?) et Time.
    #################################################################################

    sh = champ.shape
    Ndim = CMEDmesh.getSpaceDimension()

    print("dans NumpyToMEDField, shape :",sh)
    # MEDchamp = ml.MEDCouplingFieldDouble(ml.ON_CELLS) # a changer ca
    MEDchamp.setNature(ml.IntensiveMaximum)
    MEDchamp.setTime(temps[0],temps[1],temps[2])
    if len(sh) == Ndim+1:
        MEDNcompo = sh[-1]
        # MEDchamp.fillFromAnalytic(sh[-1], "0")
        # sh[-1]=3 pour un champ vectoriel; sh[-1]=1 pour un champ scalaire
        arr = np.transpose(champ, (2,1,0,3))
    elif len(sh) == Ndim+2 :
        MEDNcompo = sh[-1]*sh[-2]
        # MEDchamp.fillFromAnalytic(sh[-1]*sh[-2], "0")
        # sh[-1]*sh[-2] = 9 (3D) ou 4 (2D). On est ici pour Rij ou pour jacob
        arr = np.transpose(champ, (2,1,0,3,4))   # Nz, Ny, Nx, Ncomposante, NtermeParComposante
    else:
        MEDNcompo = 1
        print("NumpyToMEDField : Numpy has too many dimensions")
    arr = arr.flatten()
    arr = ml.DataArrayDouble(arr)
    arr.rearrange(MEDNcompo)
    MEDchamp.setArray(arr)
    return (MEDchamp)


def MEDFieldToMEDFile(MEDchamp,NomMEDFile="test.med",New=True,StructuredForced=False):
    #################################################################################
    # ATTN : il faut avoir fait un setName du champ hein ...
    #        temps = [temps (en s), iteration(int), order (-1 c carré)]
    #        MEDchamp est un DataArrayDouble bien entendu ...
    #        Le mieux c'est d'avoir créé le MEDField avec NumpyToMEDField
    #################################################################################
    if StructuredForced :
        field,mesh = en_tete(MEDchamp)
    else:
        field,mesh = MEDchamp,MEDchamp.getMesh()
        
    if New :
        if mesh.isStructured():
            ml.WriteMesh(NomMEDFile,mesh,True)
        else:
            ml.WriteUMesh(NomMEDFile,mesh,True)
    else:
        pass

    ml.WriteFieldUsingAlreadyWrittenMesh(NomMEDFile,MEDchamp)


def get_toldim(mesh,tolratio=0.01):
    ####################################################################
    # Get the tolerance parameter to give in getDifferentValues().
    #    Since the mesh is regular, the tolerance will be defined as
    #      t[i] = L[i]/N[i]*0.1
    #    t[i] : tolerance in direction i
    #    L[i] : length of th ebox in direction i
    #    N[i] : number of cells in direction i 
    #     --> Find the number of cells in one direction seems impossible 
    #     to get in an unstructured mesh : the cells are pilled up 
    #     together without predefined order I understand.
    ####################################################################
    nb_dim = mesh.getSpaceDimension()
    coord = mesh.getCoords()
    min_max = coord.getMinMaxPerComponent() # finally useless
    res = nb_dim*[None]
    for dim in xrange(nb_dim):
        coord_i = np.array(coord[:,dim].getValues())
        
        delta = np.abs(coord_i[0:-1]-coord_i[1:])
        delta = [d for d in delta if d!=0]
        
        res[dim] = tolratio*min(delta)     # 0.1 for security
    return (res)
    
    
def fromUtoC(umesh):
    ####################################################################
    # Passage d'un maillage Unstructured  un maillage Cartesien structured
    #    """ Build a MC Cartesian mesh from an unstructured one.
    #    This assumes the UMesh has all its cells properly defined, like a Cartesian mesh.
    #    @param umesh a MEDCouplingUMesh
    #    @return a MEDCouplingCMesh
    #    @return an array idx. Given 'i' a cell index in the new Cartesian
    # mesh, idx[i] is the corresponding cell index in the old UMesh.
    ####################################################################
    
    # Create a Cartesian mesh
    ret = ml.MEDCouplingCMesh.New()
    ret.setName(umesh.getName())
    axis=umesh.getSpaceDimension()*[None]
    toldim = methodes_2020.get_toldim(umesh)
    # Extract coords
    for i in xrange(umesh.getSpaceDimension()):
        toli = toldim[i]
        axis[i]=umesh.getCoords()[:,i].getDifferentValues(toli)
        axis[i].sort()
        pass
    ret.setCoords(*tuple(axis))
    # Cells renumbering - use barycenter comparison
    nc = ret.getNumberOfCells()
    # assert that we have the proper number of cells in the UMesh:
    # ##################################################################
    # #### si l'assertion est fausse c'est que toutes les coordonnees de
    # #### umesh ne osnt pas distinctes. (ce qui semble non physique a priori)
    # ##################################################################
    assert (nc == umesh.getNumberOfCells())
    baryC = ret.computeCellCenterOfMass()
    baryU = umesh.computeCellCenterOfMass()
    bary=ml.DataArrayDouble.Aggregate(baryC,baryU)
    # On fixe la tolérance pour l'association des barycentres.
    # On la fixe inférieure à la moitie de la taille de la plus petite des mailles
    # en nous assurant que les barycentres sont bien associés un à un.
    tolbary = 0.5*min(methodes_2020.get_toldim(umesh,1.))
    c,cI=bary.findCommonTuples(tolbary)
    # c  : numero des cellules qui ont les mm coordonnees à tolbary pres.
    # cI : de c[cI[i]] a c[cI[i+1]-1] on parcourt les c qui ont les mm coordonnees.
    
    # Serie de tests pour s'assurer de la bonne association des barycentres correspondants
    assert ( cI.deltaShiftIndex().isUniform(2) )
    assert ( cI.getNumberOfTuples() == (c.getNumberOfTuples() / 2) + 1 ) ; del cI
    c_copy = c.deepCopy()
    c_copy.rearrange(2)
    c_copyC = c_copy[:,0]
    c_copyU = c_copy[:,1]
    for i in range(c.getNumberOfTuples()/2):
     if c_copyC[i] >= (bary.getNumberOfTuples() / 2):
       print("Tri des barycentres C non correct")
       raise Exception("Tri des barycentres C non correct")
       pass
     elif c_copyU[i] < (bary.getNumberOfTuples() / 2):
       print("Tri des barycentres U non correct")
       raise Exception("Tri des barycentres non correct")
       pass
     pass
    c.rearrange(2)       # c'est la meme chose que c_copy a ce stade ...
    c0=c[:,0]
    try: 
        for i in range(c0.getNumberOfTuples()):
            if (c0[i] != i ): raise Exception("test failed")
        pass
    except: 
       assert ( c0.isIdentity() ) ; del c0
       pass
    # Fin des tests
    # On récupère le tableau index donnant la correspondance cartésien vers unstructured
    # idx[a] = b : a est la maille cartésien et b la maille non structurée
    idx=c[:,1]-nc
    # Autres fonction permettant de construire l'index
    # Méthode beaucoup plus couteuse en temps de calcul, et moins sûre
    # car on ne peut pas faire de test sur le résultat obtenu
    #idx = baryU.findClosestTupleId(baryC)
    return ret, idx


def UMeshToCMesh(f, mC, idx):
    ####################################################################
    # Pour verifier qu'on est bien passe d'un U mesh a un C mesh
    # MODIFICATIONS POUR CREEPING_FLOW :
    # f2.getArray().renumberInPlaceR(idx) commentee car
    # MEDLoader.InterpKernelException: DataArrayDouble::renumberInPlaceR 
    ####################################################################
    f2 = f.deepCopy()
    f2.getArray().renumberInPlaceR(idx)
    f2.setMesh(mC)   
    # Test proper assignment of the Cartesian field:
    p = methodes_2020.getCentreCoords(mC)
    print ("Did we transfer the field correctly ? ",
    f2.getValueOn(p) == f.getValueOn(p))   # on regarde si la valeur au centre est bonne (mais on regarde ni l'orientation ni rien d'autre ?)
    pass
    if not f2.getValueOn(p) == f.getValueOn(p):
        raise Exception("Field not correctly transfered")
        pass
    pass
    return f2

def getCentreCoords(mesh):
    """ ATTENTION : NE SURTOUT PAS RENVOYER UN ARRAY !!!"""
    res = tuple([ (x[0]+x[1])/2. for x in mesh.getBoundingBox()])
    return res



def MeshPlat(CMEDmesh,direction):
    ############################################################################
    # Utilisé uniquement dans Applatit (moyenne par plan). Sert à sortit un
    # Maillage élevé sur une seule maille selon la normale du plan de moyenne
    # pour supporter ce champ qui à été ramené à une seule maille dans le plan
    # de moyenne. 
    # direction : "X" ou "Y" ou "Z"
    ############################################################################
    CoordX = CMEDmesh.getCoordsAt(0)
    CoordY = CMEDmesh.getCoordsAt(1)
    CoordZ = CMEDmesh.getCoordsAt(2)
    if direction == "X":
        Xmin,Xmax = min(CoordX.getValues()),max(CoordX.getValues())
        CoordX = ml.DataArrayDouble([Xmin,Xmax])
    if direction == "Y":
        Ymin,Ymax = min(CoordY.getValues()),max(CoordY.getValues())
        CoordY = ml.DataArrayDouble([Ymin,Ymax])
    if direction == "Z":
        Zmin,Zmax = min(CoordZ.getValues()),max(CoordZ.getValues())
        CoordZ = ml.DataArrayDouble([Zmin,Zmax])
        
    CMEDmeshPlat = ml.MEDCouplingCMesh.New()
    CMEDmeshPlat.setCoords(CoordX,CoordY,CoordZ)
    CMEDmeshPlat.setName(CMEDmesh.getName()+"Plat")

    return(CMEDmeshPlat)


def get_mini_maxi(bary,delta):
    """ En fait faudrai faire ICI MEME une condition pour la périodicité, un ruc qui reboucle non ?"""

    N = len(bary)
    ind_mini_i = np.zeros(N)
    ind_maxi_i = np.zeros(N)
    
    for (i,pos) in enumerate(bary):
        shift_imin = (pos-delta)
        shift_imax = (pos+delta)

        # On le jette pck il y a un break et c mla les break # Calcul des mins : 
        # for (id_i, val) in enumerate(bary-shift_imin):
            # if val >=0.:
                # ind_mini_i[i] = id_i
                # break       ### ATTN : break : des qu'on met une valeur, piouf on ressort de la boucle if !

        # Calcul des mins :
        for (id_i, val) in enumerate(bary-shift_imin):
            if val <0.:
                ind_mini_i[i] = int(id_i+1)
        # Calcul des maxs : 
        for (id_i, val) in enumerate(bary-shift_imax):
            if val <0.:
                ind_maxi_i[i] = int(id_i)
                pass
            pass
    # 
    return ind_mini_i, ind_maxi_i



def filtre_k(champ, vol, mini_x, maxi_x, mini_y, maxi_y, mini_z, maxi_z):
    ########################################################################################
    # Moyenne du champ scalaire ou  vecteur (pas tenseur) f sur le volume vol
    # @champ : Le champ à filtrer (peut etre un champ sur maillage structuré ou non), ou un array,
    #          rempli commme il faut (par une des fonctions grad ou jacobi par exemple).
    # @mini_x : sortie de get_mini_maxi : indices des cellules ocncernées
    # @vol :
    # @fk : champ moyenné, type array. C'est plus pratique avec l'utilisation qu'on en fait dans filtre
    ########################################################################################
    
    if type(champ)!=type(np.array(351)):
        print("\n\n CE N'EST PAS UN ARRAY !!!")
        f2 = champ.deepCopy()
        if not f2.getMesh().isStructured():
            print("Unstructured in grad vect")
            mC, idx = methodes_2020.fromUtoC(f2.getMesh())
            # test_fromUtoC(f2.getMesh(), mC, idx)
            f2 = methodes_2020.UMeshToCMesh(f2, mC, idx)
        elif f2.getMesh().isStructured():
            print("Structured in grad vect")
            mC = f2.getMesh()
            pass
        champ = f2.deepCopy()

        xn = mC.getCoordsAt(0)
        yn = mC.getCoordsAt(1)
        zn = mC.getCoordsAt(2)
        # Nombre de mailles/elements : 
        nx = len(xn) - 1
        ny = len(yn) - 1
        nz = len(zn) - 1

        # Extract field in a NumPy array - we know that the internal storage
        # of a MCartesian field follows the natural (C) order (first spatial component
        # varies the quickest):
        na = np.array(champ.getArray().getValues())    
        # Re-shape to restore field component structure:
        # Test - Le gradient doit travailler sur un scalaire :
        Nco = champ.getNumberOfComponents()
        na = na.reshape(-1, Nco)
        na = na.reshape(nz, ny, nx, -1)     # le MED stocke cm ca
        champ = np.transpose(na, (2,1,0,3)) # le Numpy stocke cm ca

    else :
        pass


    Nx = len(mini_x)/10  # nombre de cellules (sommets -1 )
    Ny = len(mini_y)/10
    Nz = len(mini_z)/10
    
    sh = np.shape(champ)
    # sh_fk = np.array(sh)-1; sh_fk[-1] = sh[-1]

    """la on a quand meme une triple boucle imbriquée qui peut
        parcourir du 10 000 * 10 000 * 10 000 non ?"""
    fk = np.zeros(sh)
    for i in range(Nx):
        imin = int(mini_x[i])
        imax = int(maxi_x[i]+1) # Attention... Il ne prends pas tout sinon... c'est que le delta est plus petit que la taille d'une maille ?
        for j in range(Ny):
            jmin = int(mini_y[j])
            jmax = int(maxi_y[j]+1) # Attention... Il ne prends pas tout sinon...
            for k in range(Nz):
                if k==Nz-1 and j==Ny-1 and i==Nx-1 : print('in the loop',imin,imax,i,j,k)
                kmin = int(mini_z[k])
                kmax = int(maxi_z[k]+1) # Attention... Il ne prends pas tout sinon...
                v = vol[imin:imax,jmin:jmax,kmin:kmax]
                som_vol =  v.sum()
                # print np.shape(champ[imin:imax,jmin:jmax,kmin:kmax,:]), " et " , np.shape(v)
                ch = champ[imin:imax,jmin:jmax,kmin:kmax,:] * v[:,:,:,np.newaxis]
                som_ch =  ch.sum(axis=0).sum(axis=0).sum(axis=0) #.sum(axis=(0,1,2)) # Ne somme pas les composantes
                # print som_ch
                fk[i,j,k, :] = som_ch/som_vol
                pass
            pass
        pass

    # if 'perioX' in BC:
        # fkPerio = np.zeros(sh)
        # fkPerio[:-2,:-2,:-2,:] = fk
        
        
    #
    return fk


""" Utilise dans outils_physique.filtre_gab uniquement"""
def filtre_un_shift(Shift,champ,shiftDirection,vol):
    #########################################################################
    # Shift : [shiftX, shiftY,shiftZ] ou shiftI dicte de combien de cellules on se deplace
    # champ : champ NumPy
    # shiftDirection : "X" ou "Y" ou "Z" dicte dans quelle(s) directions on se deplace
    # vol : volume d'une cellule du maillage (regulier)
    # --> Somme le champ initial avec le champ translate dans la direction et
    #     de la distance demandee : champSum
    # --> Compte le nombre de nouveaux champs cree, dans le but de faire une
    #     moyenne : Nchamp
    #########################################################################
    sh = np.shape(champ)
    shiftX,shiftY,shiftZ = Shift
    
    # Creation des champs translatés pour faire la moyenne apres
    champShiftX,  champShiftY,  champShiftZ   = np.zeros(sh),np.zeros(sh),np.zeros(sh)
    champShiftXY, champShiftXZ, champShiftYZ  = np.zeros(sh),np.zeros(sh),np.zeros(sh)
    champShiftXYZ = np.zeros(sh)
    champSum = 0   # si on ne met pas .copy(), lorsqu'on modifie l'un, on modifie l'autre !!!! (pour contrer ce maléfice, on peut aussi faire champM = 1*champ, et les deux ne seront plus unis par ce lien mystique.)
    Nchamp = 0   #

    DeplacementX =  [shiftX, -shiftX]
    DeplacementY =  [shiftY, -shiftY]
    DeplacementZ =  [shiftZ, -shiftZ]

    # ##################################################################
    # ### X,Y,Z
    # ##################################################################
    if "X" in shiftDirection :
        for sX in DeplacementX:
            # translation selon [x] (il y en aura 2 de créées en tout)
            champShiftX[:sX,:,:,:] = champ[-sX:,:,:,:].copy()
            champShiftX[sX:,:,:,:] = champ[:-sX,:,:,:].copy()
            champSum += champShiftX * vol
            Nchamp+=1
            
    if "Y" in shiftDirection:
        for sY in DeplacementY:
            # translation selon [y] (il y en aura 2 de créées en tout)
            champShiftY[:,:sY,:,:] = champ[:,-sY:,:,:].copy()
            champShiftY[:,sY:,:,:] = champ[:,:-sY,:,:].copy()
            champSum += champShiftY * vol
            Nchamp+=1
            
    if "Z" in shiftDirection:
        for sZ in DeplacementZ:
            # translation selon [z] (il y en aura 2 de créées en tout)
            champShiftZ[:,:,:sZ,:] = champ[:,:,-sZ:,:].copy()
            champShiftZ[:,:,sZ:,:] = champ[:,:,:-sZ,:].copy()
            champSum += champShiftZ * vol
            Nchamp+=1

    # ##################################################################
    # ### XY,YZ,XZ
    # ##################################################################            
    if ("X" in shiftDirection and "Y" in shiftDirection) :
        for sX in DeplacementX:
            for sY in DeplacementY:
                # translation selon [x+y] (il y en aura 4 de créées en tout)
                champShiftXY[-sX:,:-sY,:,:] = champ[:sX,sY:,:,:].copy()
                champShiftXY[:-sX,-sY:,:,:] = champ[sX:,:sY,:,:].copy()
                champShiftXY[-sX:,-sY:,:,:] = champ[:sX,:sY,:,:].copy()
                champShiftXY[:-sX,:-sY,:,:] = champ[sX:,sY:,:,:].copy()
                champSum += champShiftXY * vol
                Nchamp+=1

    if ("Y" in shiftDirection and "Z" in shiftDirection) :
        for sY in DeplacementY:
            for sZ in DeplacementZ: 
                # translation selon [y+z] (il y en aura 4 de créées en tout)
                champShiftYZ[:,-sY:,:-sZ,:] = champ[:,:sY,  sZ:,:].copy()
                champShiftYZ[:,:-sY,-sZ:,:] = champ[:, sY:,:sZ,:].copy()
                champShiftYZ[:,-sY:,-sZ:,:] = champ[:,:sY, :sZ,:].copy()
                champShiftYZ[:,:-sY,:-sZ,:] = champ[:, sY:, sZ:,:].copy() 
                champSum += champShiftYZ * vol
                Nchamp+=1
    
    if ("X" in shiftDirection and "Z" in shiftDirection) :
        for sX in DeplacementX:
            for sZ in DeplacementZ:    
                # translation selon [x+z] (il y en aura 4 de créées en tout)
                champShiftXZ[-sX:,:,:-sZ,:] = champ[:sX,:,sZ:,:].copy()
                champShiftXZ[:-sX,:,-sZ:,:] = champ[sX:,:,:sZ,:].copy()
                champShiftXZ[-sX:,:,-sZ:,:] = champ[:sX,:,:sZ,:].copy()
                champShiftXZ[:-sX,:,:-sZ,:] = champ[sX:,:,sZ:,:].copy()
                champSum += champShiftXZ * vol
                Nchamp+=1
    
    # ##################################################################
    # ### XYZ
    # ##################################################################
    if ("X" in shiftDirection and "Y" in shiftDirection and "Z" in shiftDirection) :
        for sX in DeplacementX:
            for sY in DeplacementY:
                for sZ in DeplacementZ:
                    # translation selon [x+y+z] (il y en aura 8 de créées en tout)
                    champShiftXYZ[:sX,:sY,:sZ,:] = champ[-sX:,-sY:,-sZ:,:].copy()
                    champShiftXYZ[:sX,:sY,sZ:,:] = champ[-sX:,-sY:,:-sZ,:].copy()
                    champShiftXYZ[:sX,sY:,:sZ,:] = champ[-sX:,:-sY,-sZ:,:].copy()
                    champShiftXYZ[:sX,sY:,sZ:,:] = champ[-sX:,:-sY,:-sZ,:].copy()
                    champShiftXYZ[sX:,:sY,:sZ,:] = champ[:-sX,-sY:,-sZ:,:].copy()
                    champShiftXYZ[sX:,:sY,sZ:,:] = champ[:-sX,-sY:,:-sZ,:].copy()
                    champShiftXYZ[sX:,sY:,:sZ,:] = champ[:-sX,:-sY,-sZ:,:].copy()
                    champShiftXYZ[sX:,sY:,sZ:,:] = champ[:-sX,:-sY,:-sZ,:].copy()
                    champSum += champShiftXYZ * vol
                    Nchamp+=1

    return (champSum,Nchamp)
