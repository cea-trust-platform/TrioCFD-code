# -*- coding: utf-8 -*-
import numpy as np
import time # Pour la création d'un nouveau maillage irrégulier
import MEDLoader as ml
import methodes_2020
import methodes_mai_2020
import tests_2020
import outils_physique
import matplotlib.pyplot as plt

"""
On met ici en musique methodes_2020 et tests_2020
"""

# ################# PREPARATION ########################################
"""
On cree ici les maillages et les champs dont on va se servir pour la suite
"""

coord1 = np.linspace(0,10,21); coord2 = np.linspace(3,19,45)
coordX = ml.DataArrayDouble(coord1)
coordY = ml.DataArrayDouble(coord2)
coordZ = ml.DataArrayDouble(0.2*coord1)
# ~ xarr = ml.DataArrayDouble.New(11,1)
# ~ xarr.iota(0.)                       # Generate s, s+1, s+2, ... with a given start value s
# ~ yarr = ml.DataArrayDouble.New(21,1)
# ~ yarr.iota(0.)
# ~ zarr = ml.DataArrayDouble.New(31,1)
# ~ zarr.iota(0.)
cmesh = ml.MEDCouplingCMesh.New()
cmesh.setCoords(coordX,coordY,coordZ)
mesh = cmesh.buildUnstructured()
my_dudx = ml.MEDCouplingFieldDouble(ml.ON_CELLS, ml.ONE_TIME)
my_dudx.setMesh(mesh)
my_dudx.setName("my_dudx")
my_dudx.fillFromAnalytic(3,"(cos(x)+sin(x))*IVec+(cos(y)+sin(y))*JVec+(cos(z)-sin(z))*KVec")    # 1 means that the field should have one component

"""
On importe aussi des med déjà générés (notamment pour tester les from U to C)
et ce genre de chose.
"""
N_mesh = ml.MEDFileMesh.New('./result_N20_0000.med')
N_DUDX = ml.MEDFileField1TS('./result_N20_0000.med','DUD_X_ELEM_DOM')
N_DUDX_on_mesh = N_DUDX.getFieldOnMeshAtLevel(ml.ON_CELLS,0,N_mesh)

"""
On ecrits les champs en .vtu pour regarder sous paraview
"""
my_dudx.writeVTK('my_dudx.vtu')
N_DUDX_on_mesh.writeVTK('N_DUDX_on_mesh.vtu')

# ################### FIN PREPARATION ###################################


if True:
        T = 1*np.pi
        X = np.linspace(-2,-2+T,10);        Nx = len(X)
        Y = np.linspace(-1,-1+T,11);        Ny = len(Y)
        Z = np.linspace(-0,-0+T,12);        Nz = len(Z)
        coordX,coordY,coordZ = ml.DataArrayDouble(X),ml.DataArrayDouble(Y),ml.DataArrayDouble(Z)
        cmesh = ml.MEDCouplingCMesh.New()
        cmesh.setCoords(coordX,coordY,coordZ)
        mesh = cmesh.buildUnstructured()
        mesh.setName("Mesh")
    
        MEDchamp = ml.MEDCouplingFieldDouble.New(ml.ON_CELLS,ml.ONE_TIME)
        MEDchamp.setName("Field");          MEDchamp.setMesh(mesh)
        w = 2*np.pi/10
        valeur = np.zeros((Nx-1,Ny-1,Nz-1,3))
        for i,x in enumerate(X[:-1]):
            for j,y in enumerate(Y[:-1]):
                for k,z in enumerate(Z[:-1]):
                    valeur[i,j,k,0] = 10*np.cos(x)
                    valeur[i,j,k,1] = np.cos(z)
                    valeur[i,j,k,2] = np.cos(10*y)
        mod_valeur = valeur.transpose(2,1,0,3);
        mod_valeur = mod_valeur.flatten()
        champ = ml.DataArrayDouble(mod_valeur)
        champ.rearrange(3)
        MEDchamp.setArray(champ)
        MEDchamp.setTime(0,0,-1)

        testMEDchamp = methodes_2020.NumpyToMEDField_VRij(valeur,MEDchamp,cmesh,[0,0,-1])

        dArray = MEDchamp.getArray()-testMEDchamp.getArray()

        dAx = np.array(dArray[:,0].getValues());
        dAy = np.array(dArray[:,1].getValues());
        dAz = np.array(dArray[:,2].getValues())
        lx = np.linalg.norm(dAx); ly = np.linalg.norm(dAy); lz = np.linalg.norm(dAz);
        succes = (lx+ly+lz)==0
        

        message = ("\nMauvaise converion numpy->med. 1 pour OK. 0 sinon : \n"+
            "Composante I \t {0} \nComposante J \t {1} \nComposante K \t {2}".format(
            int(lx==0),int(ly==0),int(lz==0))
            )
        self.assertTrue(succes,message)






if False :
    ####################################################################
    ########################                    ########################
    ############               TEST RIJ                  ###############
    ########################                    ########################
    ####################################################################

    T = 2*np.pi
    X, Y, Z = np.linspace(-2,-2+T,20), np.linspace(-100,-100+1*T,20), np.linspace(80,80+1*T,20)
    Nx,Ny,Nz = len(X),len(Y),len(Z)
    coordX = ml.DataArrayDouble(X)
    coordY = ml.DataArrayDouble(Y)
    coordZ = ml.DataArrayDouble(Z)#X+Y); Z = X+Y
    cmesh = ml.MEDCouplingCMesh.New()
    cmesh.setCoords(coordX,coordY,coordZ)
    mesh = cmesh.buildUnstructured()
    mesh.setName("Mesh")

    ml.WriteUMesh("Rij_test.med",mesh,True)
    iterations = [(i,-1) for i in range(61)] # pour avoir jusqu'à 30 inclus

    MEDchamp = ml.MEDCouplingFieldDouble.New(ml.ON_CELLS,ml.ONE_TIME)
    MEDchamp.setName("velocity")
    MEDchamp.setMesh(mesh)

    for it in iterations:
        w = 2*np.pi/60
        valeur = np.zeros((Nx-1,Ny-1,Nz-1,3))
        goodRij = np.zeros((Nx-1,Ny-1,Nz-1,3,3))
        for i,x in enumerate(X[:-1]):
            for j,y in enumerate(Y[:-1]):
                for k,z in enumerate(Z[:-1]):
                    valeur[i,j,k,0] = np.cos(x) + np.cos(1.0*w*it[0]) 
                    valeur[i,j,k,1] = np.cos(y) + np.cos(10.*w*it[0])
                    valeur[i,j,k,2] = np.cos(z) + np.cos(100*w*it[0])
    
                    goodRij[i,j,k,0,0] = np.cos(1.0*w*it[0])*np.cos(1.0*w*it[0])#np.cos(x)*np.cos(x)
                    goodRij[i,j,k,0,1] = np.cos(1.0*w*it[0])*np.cos(10.*w*it[0])#np.cos(x)*np.cos(y)
                    goodRij[i,j,k,0,2] = np.cos(1.0*w*it[0])*np.cos(100*w*it[0])#np.cos(x)*np.cos(z)
                    goodRij[i,j,k,1,0] = np.cos(10.*w*it[0])*np.cos(1.0*w*it[0])#np.cos(y)*np.cos(x)
                    goodRij[i,j,k,1,1] = np.cos(10.*w*it[0])*np.cos(10.*w*it[0])#np.cos(y)*np.cos(y)
                    goodRij[i,j,k,1,2] = np.cos(10.*w*it[0])*np.cos(100*w*it[0])#np.cos(y)*np.cos(z)
                    goodRij[i,j,k,2,0] = np.cos(100*w*it[0])*np.cos(1.0*w*it[0])#np.cos(z)*np.cos(x)
                    goodRij[i,j,k,2,1] = np.cos(100*w*it[0])*np.cos(10.*w*it[0])#np.cos(z)*np.cos(y)
                    goodRij[i,j,k,2,2] = np.cos(100*w*it[0])*np.cos(100*w*it[0])#np.cos(z)*np.cos(z)
    
        valeur = valeur.transpose(2,1,0,3)
        valeur = valeur.flatten()
        champ = ml.DataArrayDouble(valeur)
        champ.rearrange(3)
        MEDchamp.setArray(champ)
        MEDchamp.setTime(it[0],it[0],it[1])
        
        ml.WriteFieldUsingAlreadyWrittenMesh("Rij_test.med",MEDchamp)

                            
    NPRij,MEDRij = outils_physique.Rij("Rij_test.med","velocity")
    print(np.shape(NPRij))

    """Norme L2 poue evaluer de la justesse comme d'hab"""
    err = goodRij - NPRij
    L2errTT = np.linalg.norm(err)           ;L2TT = np.linalg.norm(goodRij)           
    L2errXX = np.linalg.norm(err[:,:,:,0,0]);L2XX = np.linalg.norm(goodRij[:,:,:,0,0])
    L2errXY = np.linalg.norm(err[:,:,:,0,1]);L2XY = np.linalg.norm(goodRij[:,:,:,0,1])
    L2errXZ = np.linalg.norm(err[:,:,:,0,2]);L2XZ = np.linalg.norm(goodRij[:,:,:,0,2])
    L2errYX = np.linalg.norm(err[:,:,:,1,0]);L2YX = np.linalg.norm(goodRij[:,:,:,1,0])
    L2errYY = np.linalg.norm(err[:,:,:,1,1]);L2YY = np.linalg.norm(goodRij[:,:,:,1,1])
    L2errYZ = np.linalg.norm(err[:,:,:,1,2]);L2YZ = np.linalg.norm(goodRij[:,:,:,1,2])
    L2errZX = np.linalg.norm(err[:,:,:,2,0]);L2ZX = np.linalg.norm(goodRij[:,:,:,2,0])
    L2errZY = np.linalg.norm(err[:,:,:,2,1]);L2ZY = np.linalg.norm(goodRij[:,:,:,2,1])
    L2errZZ = np.linalg.norm(err[:,:,:,2,2]);L2ZZ = np.linalg.norm(goodRij[:,:,:,2,2])

    pcTT = round(100*L2errTT / L2TT , 2); succesTT = pcTT<10
    pcXX = round(100*L2errXX / L2XX , 2); succesXX = pcXX<10
    pcXY = round(100*L2errXY / L2XY , 2); succesXY = pcXY<10
    pcXZ = round(100*L2errXZ / L2XZ , 2); succesXZ = pcXZ<10
    pcYX = round(100*L2errYX / L2YX , 2); succesYX = pcYX<10
    pcYY = round(100*L2errYY / L2YY , 2); succesYY = pcYY<10
    pcYZ = round(100*L2errYZ / L2YZ , 2); succesYZ = pcYZ<10
    pcZX = round(100*L2errZX / L2ZX , 2); succesZX = pcZX<10
    pcZY = round(100*L2errZY / L2ZY , 2); succesZY = pcZY<10
    pcZZ = round(100*L2errZZ / L2ZZ , 2); succesZZ = pcZZ<10

    message = ('ERREUR DANS LES TENSIONS DE REYNOLDS \t\t (%)'+'\n \n'+
    'Totale : \t\t {0} \n'.format(pcTT)+
    'Selon X : \t\t {0}, {1}, {2} \n'.format(pcXX,pcXY,pcXZ)+
    'Selon Y : \t\t {0}, {1}, {2} \n'.format(pcYX,pcYY,pcYZ)+
    'Selon Z : \t\t {0},  {1}, {2}\n'.format(pcZX,pcZY,pcZZ))
    print(succesTT,
          succesXX,
          succesXY,
          succesXZ,
          succesYX,
          succesYY,
          succesYZ,
          succesZX,
          succesZY,
          succesZZ
    )
    print(message)


if False:
    """ test pour appltit (moyenne par plan) """
    T = 2*np.pi
    X, Y, Z = np.linspace(-2,-2+T,60), np.linspace(-100,-100+1*T,60), np.linspace(80,80+1*T,60)
    Nx,Ny,Nz = len(X),len(Y),len(Z)
    coordX = ml.DataArrayDouble(X)
    coordY = ml.DataArrayDouble(Y)
    coordZ = ml.DataArrayDouble(Z)#X+Y); Z = X+Y
    cmesh = ml.MEDCouplingCMesh.New()
    cmesh.setCoords(coordX,coordY,coordZ)
    mesh = cmesh.buildUnstructured()
    mesh.setName("Mesh")

    temps = [0,0,-1]
    MEDchamp = ml.MEDCouplingFieldDouble.New(ml.ON_CELLS,ml.ONE_TIME)
    MEDchamp.setName("Field")
    MEDchamp.setMesh(mesh)
    w = 2*np.pi/30
    valeur = np.zeros((Nx-1,Ny-1,Nz-1,3))
    for i,x in enumerate(X[:-1]):
        for j,y in enumerate(Y[:-1]):
            for k,z in enumerate(Z[:-1]):
                valeur[i,j,k,0] = np.cos(x)         
                valeur[i,j,k,1] = np.cos(y) 
                valeur[i,j,k,2] = np.cos(z)     
    valeur = valeur.transpose(2,1,0,3)
    valeur = valeur.flatten()
    champ = ml.DataArrayDouble(valeur)
    champ.rearrange(3)
    MEDchamp.setArray(champ)
    MEDchamp.setTime(temps[0],temps[0],temps[1])

    NPchampPlatX,MEDXchampPlatX = outils_physique.Applatit(MEDchamp,temps,"X")
    NPchampPlatY,MEDXchampPlatY = outils_physique.Applatit(MEDchamp,temps,"Y")
    NPchampPlatZ,MEDXchampPlatZ = outils_physique.Applatit(MEDchamp,temps,"Z")
    
    ##########################################################################
    # Remarque : on pourrait définir un "goodchamPlat" et fair ela norme L2 de
    # l'erreur. Ici ce "goodchamp" doit re remplit de 0, et soustraire 0 à qqch
    # c'est comme ne rien faire. On se passera donc de ce "goodchamp".
    ##########################################################################
    l2A = 100*round(np.linalg.norm(NPchampPlatX[0,:,:,0]),3)
    l2B = 100*round(np.linalg.norm(NPchampPlatY[:,0,:,1]),3)
    l2C = 100*round(np.linalg.norm(NPchampPlatZ[:,:,0,2]),3)

    succesA = l2A<5;    succesB = l2B<5;    succesC = l2C<5
    succes = (succesA and succesB and succesC)

    message = ('Problème dans la moyene par plan, exigence à 5% près :\n'+
                '--> 1 pour OK, 0 sinon\n'+
                'Moyenne selon X (résultat dans le plan Y,Z)\t{0}\n'.format(int(succesA))+
                'Moyenne selon Y (résultat dans le plan X,Z)\t{0}\n'.format(int(succesB))+
                'Moyenne selon Z (résultat dans le plan X,Y)\t{0}\n'.format(int(succesC))
                )

    print(succes)
    print(message)
    if False:
        plt.figure(2)
        plt.plot(NPchampPlatZ[:,10,0,0])
        plt.plot(NPchampPlatZ[:,10,0,1])
        plt.plot(NPchampPlatZ[:,10,0,2])
        plt.legend(fontsize=10)
        plt.xticks(fontsize='15')
        plt.yticks(fontsize='15')
        
        plt.figure(3)
        plt.plot(NPchampPlatZ[10,:,0,0])
        plt.plot(NPchampPlatZ[10,:,0,1])
        plt.plot(NPchampPlatZ[10,:,0,2])
        plt.legend(fontsize=10)
        plt.xticks(fontsize='15')
        plt.yticks(fontsize='15')

    if False:
        xx,yy = np.meshgrid(Y[:-1],X[:-1])

        plt.figure(1)
        plt.subplot(1,3,1)
        a = plt.contour(xx,yy,NPchampPlatZ[:,:,0,0],label='i')
        plt.clabel(a,fontsize=10, inline=1)
        plt.legend(fontsize=10)
        plt.xticks(fontsize='15')
        plt.yticks(fontsize='15')
        
        plt.subplot(1,3,2)
        b = plt.contour(xx,yy,NPchampPlatZ[:,:,0,1],label='j')
        plt.clabel(b,fontsize=10, inline=1)
        plt.legend(fontsize=10)
        plt.xticks(fontsize='15')
        plt.yticks(fontsize='15')
        
        plt.subplot(1,3,3)
        c = plt.contour(xx,yy,NPchampPlatZ[:,:,0,2],label='k')
        plt.clabel(c,fontsize=10, inline=1)
        plt.legend(fontsize=10)
        plt.xticks(fontsize='15')
        plt.yticks(fontsize='15')
    plt.show()
    
if False:
    MEDfileField,MEDmesh,TypeOfField = methodes_mai_2020.MEDFileToMEDField_MultiTS("ptit_test.med","Field")
    MEDfileField22 = MEDfileField.getTimeStep(22,-1)
    MEDfield22 = MEDfileField.getFieldOnMeshAtLevel(TypeOfField,22,-1,0,MEDmesh,0)
    print()
    
   
    
if False:
    """ creation d'un champ sur plusieurs TS"""
    T = 2*np.pi
    X, Y, Z = np.linspace(-2,-2+T,120), np.linspace(-100,-100+1*T,60), np.linspace(80,80+1*T,60)
    Nx,Ny,Nz = len(X),len(Y),len(Z)
    coordX = ml.DataArrayDouble(X)
    coordY = ml.DataArrayDouble(Y)
    coordZ = ml.DataArrayDouble(Z)#X+Y); Z = X+Y
    cmesh = ml.MEDCouplingCMesh.New()
    cmesh.setCoords(coordX,coordY,coordZ)
    mesh = cmesh.buildUnstructured()
    mesh.setName("Mesh")

    ml.WriteUMesh("3004_test.med",mesh,True)

    iterations = [(i,-1) for i in range(61)] # pour avoir jusqu'à 60 inclus

    MEDchamp = ml.MEDCouplingFieldDouble.New(ml.ON_CELLS,ml.ONE_TIME)
    MEDchamp.setName("Field")
    MEDchamp.setMesh(mesh)
    for it in iterations:
        w = 2*np.pi/60
        valeur = np.zeros((Nx-1,Ny-1,Nz-1,3))
        for i,x in enumerate(X[:-1]):
            for j,y in enumerate(Y[:-1]):
                for k,z in enumerate(Z[:-1]):
                    valeur[i,j,k,0] = 3+np.cos(w*it[0])
                    valeur[i,j,k,1] = np.cos(x)+np.cos(w*it[0])
                    valeur[i,j,k,2] = it[0]+np.cos(w*it[0])
        valeur = valeur.transpose(2,1,0,3)
        valeur = valeur.flatten()
        champ = ml.DataArrayDouble(valeur)
        champ.rearrange(3)
        MEDchamp.setArray(champ)
        MEDchamp.setTime(it[0],it[0],it[1])
        ml.WriteFieldUsingAlreadyWrittenMesh("3004_test.med",MEDchamp)


if False:
    T = 2*np.pi
    X, Y, Z = np.linspace(-2,-2+T,120), np.linspace(-100,-100+1*T,60), np.linspace(80,80+1*T,60)
    Nx,Ny,Nz = len(X),len(Y),len(Z)
    
    champMean,MEDMean = methodes_2020.TimeAvg('3004_test.med','Field')
    ml.WriteField('3004_mean.med',MEDMean,True)

    goodMean = np.zeros((Nx-1,Ny-1,Nz-1,3))
    goodMean[...,0] = 3; goodMean[...,2] = 30
    for i,x in enumerate(X[:-1]):
        goodMean[i,:,:,1] = np.cos(x)
        
    err = goodMean-champMean
    print(err.shape)
    
    L2errT = np.linalg.norm(err);       L2T = np.linalg.norm(champMean)
    L2errX = np.linalg.norm(err[...,0]);L2X = np.linalg.norm(champMean[...,0])
    L2errY = np.linalg.norm(err[...,1]);L2Y = np.linalg.norm(champMean[...,1])
    L2errZ = np.linalg.norm(err[...,2]);L2Z = np.linalg.norm(champMean[...,2])
    
    pcT = round(100*L2errT / L2T , 2)
    pcX = round(100*L2errX / L2X , 2)
    pcY = round(100*L2errY / L2Y , 2)
    pcZ = round(100*L2errZ / L2Z , 2)

    message = ('ERREUR DANS LA MOYENNE TEMPORELLE \t\t (%)'+'\n \n'+
    'Totale : \t\t {0} \n'.format(pcT)+ 'Selon X : \t\t {0} \n'.format(pcX)+
    'Selon Y : \t\t {0} \n'.format(pcY)+'Selon Z : \t\t {0} \n'.format(pcZ))

    print(message)
    
    plt.figure(1)
    plt.plot(goodMean[:,2,2,0],'k',label='X:3')
    plt.plot(goodMean[:,2,2,1],'r',label='Y:cos')
    plt.plot(goodMean[:,2,2,2],'b',label='Z:3')
    plt.plot(champMean[:,2,2,0],'--k',label='MX')
    plt.plot(champMean[2,2,:,1],'--m',label='MY')
    plt.plot(champMean[:,2,2,1],'--r',label='MY')
    plt.plot(champMean[:,2,2,2],'--b',label='MZ')
    plt.legend(fontsize=10)
    plt.xticks(fontsize='15')
    plt.yticks(fontsize='15')
    plt.show()


    


"""
    # ECRITURE DU CHAMP
    champ = ml.MEDCouplingFieldDouble(ml.ON_CELLS,ml.LINEAR_TIME)
    champ.setNature(ml.IntensiveMaximum)
    champ.setMesh(mesh)
    t = 12
    champ.setTime(t,0,0)   #val,iteration,order
    champ.fillFromAnalytic(1,"pl")
"""

if False:
    MEDFileName =   './MEDFILES/res_monte_2804.med'
    FieldName =     "VELOCITY_ELEM_DOM_dual"
    champTAvg, MEDcTAvg = methode_2020.TimeAvg(MEDFileName,FieldName)
    MEDcTAvg.writeVTK("TEEEEST",MEDcTAvg,True)


if False:

    monteMesh = ml.MEDFileMesh.New('./MEDFILES/res_monte_2804.med')
    print(monteMesh.getIteration())
    monteField = ml.MEDFileFieldMultiTS.New("./MEDFILES/res_monte_2804.med","VELOCITY_ELEM_DOM_dual",True)
    it = monteField.getIterations()

    champSomme = 0
    for i in it:
        champSomme += monteField.getFieldOnMeshAtLevel(ml.ON_CELLS,i[0],i[0],0,monteMesh,0)
    champMean = champSomme/len(it)

    
    champ1 = monteField.getFieldOnMeshAtLevel(ml.ON_CELLS,102,-1,0,monteMesh,0)
    champ1 = champ1.getArray
    print(champ1)
    sh = champ1.shape
    print(sh)

if False:
    """ moyenne temporelle """
    T = 2*np.pi
    X, Y, Z = np.linspace(-2,-2+T,60), np.linspace(-100,-100+1*T,60), np.linspace(80,80+1*T,60)
    Nx,Ny,Nz = len(X),len(Y),len(Z)
    coordX = ml.DataArrayDouble(X)
    coordY = ml.DataArrayDouble(Y)
    coordZ = ml.DataArrayDouble(Z)#X+Y); Z = X+Y
    cmesh = ml.MEDCouplingCMesh.New()
    cmesh.setCoords(coordX,coordY,coordZ)
    mesh = cmesh.buildUnstructured()

    # ECRITURE DU CHAMP
    champ = ml.MEDCouplingFieldDouble(ml.ON_CELLS,ml.LINEAR_TIME)
    champ.setNature(ml.IntensiveMaximum)
    champ.setMesh(mesh)
    t = 12
    champ.setTime(t,0,0)   #val,iteration,order
    champ.fillFromAnalytic(1,"pl")

    # LECTURE DU CHAMP
    champNumpy = np.array(champ.getArray().getValues())
    champNumpy = champNumpy.reshape(Nz-1,Ny-1,Nx-1,1)
    champNumpy = champNumpy.transpose(2,1,0,3)
    print(champNumpy[:3,:3,:3,:])
    
if False:
    """filtre"""
    
    T = 2*np.pi
    X, Y, Z = np.linspace(-2,-2+T,630), np.linspace(-100,-100+1*T,60), np.linspace(80,80+1*T,60)
    Nx,Ny,Nz = len(X),len(Y),len(Z)
    coordX = ml.DataArrayDouble(X)
    coordY = ml.DataArrayDouble(Y)
    coordZ = ml.DataArrayDouble(Z)#X+Y); Z = X+Y
    cmesh = ml.MEDCouplingCMesh.New()
    cmesh.setCoords(coordX,coordY,coordZ)
    mesh = cmesh.buildUnstructured()

    MainchampScalNmpy = np.zeros((Nx,Ny,Nz,1))
    PertchampScalNmpy = np.zeros((Nx,Ny,Nz,1))
    TOTchampScalNmpy = np.zeros((Nx,Ny,Nz,1))

    Lambda = 0.1
    for i,x in enumerate(X):
        for j,y in enumerate(Y):
            for k,z in enumerate(Z):
                Mc =      np.cos(2*np.pi       *x)#+2*np.cos(y)+3*np.cos(z)
                Pc = 0.1*(np.cos(2*np.pi/Lambda*x))#+np.cos(d*y)+np.cos(d*z))
                Tc = Mc+Pc
                MainchampScalNmpy[i,j,k] = Mc
                PertchampScalNmpy[i,j,k] = Pc
                TOTchampScalNmpy[i,j,k] = Tc


    delta = Lambda
    BC='perioX,perioY,perioZ'
    Lx = X[-1]-X[0]
    Ly = Y[-1]-Y[0]
    Lz = Z[-1]-Z[0]
    tailles = np.array([Lx, Ly, Lz])
    # filtre1 = methodes_2020.filtre_gab(delta,TOTchampScalNmpy,tailles,BC)
    filtredemi = methodes_2020.filtre_gab(delta/2,TOTchampScalNmpy,tailles,BC)


    import matplotlib.pyplot as plt
    plt.figure(1)
    plt.plot(X,TOTchampScalNmpy[:,2,2,0],'.',label='pasfi')
    plt.plot(X,MainchampScalNmpy[:,2,2,0],'-',label='filtre1')
    plt.plot(X,filtredemi[:,2,2,0],'-',label='filtredemi')
    plt.legend(fontsize=10)
    plt.xticks(fontsize='15')
    plt.yticks(fontsize='15')

    plt.show()
    

    # print("champU_filtre", np.shape(champU_filtre))
    # print("MainchampScalNmpy",np.shape(MainchampScalNmpy))

    """ on va encore prendre la norme L2 de l'erreur ppour voir si c bon ou pas."""

    L2main = np.linalg.norm(np.abs(MainchampScalNmpy))

    print(L2main)
    errTOT = filtredemi-MainchampScalNmpy; L2errTOT = np.linalg.norm(errTOT)
    print(L2errTOT)
    errX = errTOT[:,1,1,:]; L2errX   = np.linalg.norm(errX)
    errY = errTOT[1,:,1,:]; L2errY   = np.linalg.norm(errY)
    errZ = errTOT[1,1,:,:]; L2errZ   = np.linalg.norm(errZ)

    pcT = round(100*L2errTOT / L2main , 2)
    pcX = round(100*L2errX / L2main , 2)
    pcY = round(100*L2errY / L2main , 2)
    pcZ = round(100*L2errZ / L2main , 2)

    message = ('ERREUR DANS LE FILTRE \t\t (%)'+'\n \n'+
    'Totale : \t\t {0} \n'.format(pcT)+ 'Selon X : \t\t {0} \n'.format(pcX)+
    'Selon Y : \t\t {0} \n'.format(pcY)+'Selon Z : \t\t {0} \n'.format(pcZ))
    print(message)

if False :
        fig, axs = plt.subplots(1, sharex=True, sharey=True)
        ###   ###   ###
        plt.figure(1)
        plt.suptitle('Derivee par rapport a X',fontsize='15')
        plt.subplot(1,3,1)
        plt.title('along X',fontsize=15)
        plt.xlabel('X',fontsize=15)
        plt.plot(MainchampScalNmpy[:,2,2],'.',label='main')
        plt.plot(PertchampScalNmpy[:,2,2],label='pert')
        plt.plot(TOTchampScalNmpy[:,2,2],label='tot')
        plt.plot(champU_filtre[:,2,2],label='filtre')
        plt.xticks(fontsize='15')
        plt.yticks(fontsize='15')
        plt.legend(fontsize=15)
    
        plt.subplot(1,3,2)
        plt.title('along Y',fontsize=15)
        plt.xlabel('Y',fontsize=15)
        plt.plot(MainchampScalNmpy[2,:,2],'.',label='main')
        plt.plot(PertchampScalNmpy[2,:,2],label='pert')
        plt.plot(TOTchampScalNmpy[2,:,2],label='tot')
        plt.plot(champU_filtre[2,:,2],label='filtre')
        plt.xticks(fontsize='15')
        plt.yticks(fontsize='15')
        plt.legend(fontsize=15)
    
        plt.subplot(1,3,3)
        plt.title('along Z',fontsize=15)
        plt.xlabel('Z',fontsize=15)
        plt.plot(MainchampScalNmpy[2,2,:],'.',label='main')
        plt.plot(PertchampScalNmpy[2,2,:],label='pert')
        plt.plot(TOTchampScalNmpy[2,2,:],label='tot')
        plt.plot(champU_filtre[2,2,:],label='filtre')
        plt.xticks(fontsize='15')
        plt.yticks(fontsize='15')
        plt.legend(fontsize=15)
        plt.show()


if False:
    """get_mini_maxi"""
    delta = 2
    
    # Récupération des corrdonnées des noeuds du maillage
    noeudx = cmesh.getCoordsAt(0).toNumPyArray()
    noeudy = cmesh.getCoordsAt(1).toNumPyArray()
    noeudz = cmesh.getCoordsAt(2).toNumPyArray()
    #
    Nx = len(noeudx)-1
    Ny = len(noeudy)-1
    Nz = len(noeudz)-1
    # Récupération des coordonnées des centre des cellules
    baryx = (noeudx[1:] + noeudx[:-1])/2
    baryy = (noeudy[1:] + noeudy[:-1])/2
    baryz = (noeudz[1:] + noeudz[:-1])/2
    #
    # Calcul de la longeur des cotés des cellules
    dx = (noeudx[1:] - noeudx[:-1])
    dxx = dx.reshape(Nx,1,1)
    dy = (noeudy[1:] - noeudy[:-1])
    dyy = dy.reshape(1,Ny,1)
    dz = (noeudz[1:] - noeudz[:-1])
    dzz = dz.reshape(1,1,Nz)
    #
    # Calcul du volume des cellules
    vol = dxx*dyy*dzz
    #
    # Création des tableaux mins et maxs
    mini_x, maxi_x = methodes_2020.get_mini_maxi(baryx,delta)
    mini_y, maxi_y = methodes_2020.get_mini_maxi(baryy,delta)
    mini_z, maxi_z = methodes_2020.get_mini_maxi(baryz,delta)

    
    # ##################################################################
    # ## test get_mini_maxi ###
    good_mini_z = [0, 0, 0, 0, 0 ]+[i for i in range(25)]
    good_maxi_z = [i for i in range(4,29)]+[29, 29, 29, 29, 29, 29]

    succes_mini = True
    for (i,mz) in enumerate(mini_z):
        succes_mini = (succes_mini and mz==good_mini_z[i])

    succes_maxi = True
    for (i,mz) in enumerate(maxi_z):
        succes_maxi = (succes_maxi and mz==good_maxi_z[i])

    # ## fin test get_mini_maxi ###
    # ##################################################################


       


if False :
    """ jacob"""


    """  POINT SUR L ECRITURE DE JAOB :

                           | dxf    dyf    dzf |     | j[0,0]  j[0,1]  j[0,2] |
        jacob((f,g,h)^T) = | dxg    dyg    dzg | "=" | j[1,0]  j[1,1]  j[1,2] |
                           | dxh    dyh    dzh |     | j[2,0]  j[2,1]  j[2,2] |

        --> respecte les conventions usuelles d'écriture.
    """
    def df(Coord):
        x, y, z = Coord[0], Coord[1], Coord[2]
        """
          f(x,y,z) = cos(x)+cos(y)+cos(z)
        """
        dxf = -np.sin(x)
        dyf = -np.sin(y/2)/2
        dzf = -np.sin(z/3)/3

        return(dxf,dyf,dzf)

    T = 2*np.pi
    X, Y, Z = np.linspace(-2,-2+T,20), np.linspace(-1,-1+2*T,21), np.linspace(0,0+3*T,22)
    coordX = ml.DataArrayDouble(X)
    coordY = ml.DataArrayDouble(Y)
    coordZ = ml.DataArrayDouble(Z)#X+Y); Z = X+Y
    cmesh = ml.MEDCouplingCMesh.New()
    cmesh.setCoords(coordX,coordY,coordZ)
    mesh = cmesh.buildUnstructured()
    
    champU = ml.MEDCouplingFieldDouble(ml.ON_CELLS)
    champU.setNature(ml.IntensiveMaximum)
    champU.setMesh(mesh)
    champU.fillFromAnalytic(3,  "1*(cos(x)+cos(y/2)+cos(z/3))*IVec +"+
                                "2*(cos(x)+cos(y/2)+cos(z/3))*JVec +"+
                                "3*(cos(x)+cos(y/2)+cos(z/3))*KVec")


    BC = 'perioX, perioY, perioZ'
    jacob, MEDjacob = methodes_2020.jacob(champU,BC)
    print(np.shape(jacob))
    
    good_jacob = np.zeros(np.shape(jacob))
    xf = (X[0]+ X[-1])/2; X = (X[:-1]+X[1:])/2; X = np.concatenate((X,xf),axis = None)
    yf = (Y[0]+ Y[-1])/2; Y = (Y[:-1]+Y[1:])/2; Y = np.concatenate((Y,yf),axis = None)
    zf = (Z[0]+ Z[-1])/2; Z = (Z[:-1]+Z[1:])/2; Z = np.concatenate((Z,zf),axis = None)
    for i in range(len(X)-1):
        x = X[i]
        for j in range(len(Y)-1):
            y = Y[j]
            for k in range(len(Z)-1):
                z = Z[k]
                good_jacob[i,j,k,0,:] = 1*df([x,y,z])
                good_jacob[i,j,k,1,:] = 2*np.array(df([x,y,z]))
                good_jacob[i,j,k,2,:] = 3*np.array(df([x,y,z]))

    """on recupere les normes pour évaluer si on est bon ou pas : """
    L2_good_jacob = np.linalg.norm(good_jacob)
    L2_err_jacob =  np.linalg.norm(jacob - good_jacob)

    L2_good_jacob_cI =     np.linalg.norm(good_jacob[:,:,:,0,:])
    L2_err_jacob_cI = np.linalg.norm(jacob[:,:,:,0,:] - good_jacob[:,:,:,0,:])
    L2_good_jacob_cJ =     np.linalg.norm(good_jacob[:,:,:,1,:])
    L2_err_jacob_cJ = np.linalg.norm(jacob[:,:,:,1,:] - good_jacob[:,:,:,1,:])    
    L2_good_jacob_cK =     np.linalg.norm(good_jacob[:,:,:,2,:])
    L2_err_jacob_cK = np.linalg.norm(jacob[:,:,:,2,:] - good_jacob[:,:,:,2,:])

    L2_good_jacob_dX =     np.linalg.norm(good_jacob[:,:,:,:,0])
    L2_err_jacob_dX = np.linalg.norm(jacob[:,:,:,:,0] - good_jacob[:,:,:,:,0])
    L2_good_jacob_dY =     np.linalg.norm(good_jacob[:,:,:,:,1])
    L2_err_jacob_dY = np.linalg.norm(jacob[:,:,:,:,1] - good_jacob[:,:,:,:,1])
    L2_good_jacob_dZ =     np.linalg.norm(good_jacob[:,:,:,:,2])
    L2_err_jacob_dZ = np.linalg.norm(jacob[:,:,:,:,2] - good_jacob[:,:,:,:,2])

    pcT = round(100*L2_err_jacob/L2_good_jacob      ,2)
    pcI = round(100*L2_err_jacob_cI/L2_good_jacob_cI,2) 
    pcJ = round(100*L2_err_jacob_cJ/L2_good_jacob_cJ,2) 
    pcK = round(100*L2_err_jacob_cK/L2_good_jacob_cK,2) 
    pcX = round(100*L2_err_jacob_dX/L2_good_jacob_dX,2) 
    pcY = round(100*L2_err_jacob_dY/L2_good_jacob_dY,2) 
    pcZ = round(100*L2_err_jacob_dZ/L2_good_jacob_dZ,2)


    message = ('ERREURS : \t\t   (%)'+'\n \n'+
    'Totale : \t\t {0} \n'.format(pcT)+'Composante I : \t\t {0} \n'.format(pcI)+
    'Composante J : \t\t {0} \n'.format(pcJ)+'Composante K : \t\t {0} \n'.format(pcK)+
    'Dérivée selon X : \t {0} \n'.format(pcX)+'Dérivée selon Y : \t {0} \n'.format(pcY)+
    'Dérivée selon Z : \t {0} \n'.format(pcZ))
    print(message)

    
    import matplotlib.pyplot as plt

    fig, axs = plt.subplots(1, sharex=True, sharey=True)
    ###   ###   ###
    plt.figure(1)
    plt.suptitle('Derivee par rapport a X',fontsize='15')

    plt.subplot(3,1,1)
    plt.plot(X[:-1],jacob[:,4,4,0,0],'.r',label='grad meth')   # dxf = -sin(x)
    plt.plot(X[:-1],good_jacob[:,4,4,0,0],'.k',label='good grad')
    plt.xticks(fontsize='15')
    plt.yticks(fontsize='15')
    plt.xlabel('X',fontsize='10')
    plt.legend(fontsize='15')

    plt.subplot(3,1,2)
    plt.plot(Y[:-1],jacob[4,:,4,0,1],'.r',label='grad meth')
    plt.plot(Y[:-1],good_jacob[4,:,4,0,1],'k',label='good grad')
    plt.xticks(fontsize='15')
    plt.yticks(fontsize='15')
    plt.xlabel('Y',fontsize='10')

    plt.subplot(3,1,3)
    plt.plot(Z[:-1],jacob[4,4,:,0,2],'.r',label='grad meth')
    plt.plot(Z[:-1],good_jacob[4,4,:,0,2],'k',label='good grad')
    plt.xticks(fontsize='15')
    plt.yticks(fontsize='15')
    plt.xlabel('Z',fontsize='10')

    
    fig, axs = plt.subplots(2, sharex=True, sharey=True)
    ###   ###   ###
    plt.figure(2)
    plt.suptitle('Derivee par rapport a X',fontsize='15')

    plt.subplot(3,1,1)
    plt.plot(X[:-1],jacob[:,4,4,1,0],'.r',label='grad meth')   # dxf = -sin(x)
    plt.plot(X[:-1],good_jacob[:,4,4,1,0],'.k',label='good grad')
    plt.xticks(fontsize='15')
    plt.yticks(fontsize='15')
    plt.xlabel('X',fontsize='10')
    plt.legend(fontsize='15')

    plt.subplot(3,1,2)
    plt.plot(Y[:-1],jacob[4,:,4,1,1],'.r',label='grad meth')
    plt.plot(Y[:-1],good_jacob[4,:,4,1,1],'k',label='good grad')
    plt.xticks(fontsize='15')
    plt.yticks(fontsize='15')
    plt.xlabel('Y',fontsize='10')

    plt.subplot(3,1,3)
    plt.plot(Z[:-1],jacob[4,4,:,1,2],'.r',label='grad meth')
    plt.plot(Z[:-1],good_jacob[4,4,:,1,2],'k',label='good grad')
    plt.xticks(fontsize='15')
    plt.yticks(fontsize='15')
    plt.xlabel('Z',fontsize='10')

    
    fig, axs = plt.subplots(3, sharex=True, sharey=True)
    ###   ###   ###
    plt.figure(3)
    plt.suptitle('Derivee par rapport a X',fontsize='15')

    plt.subplot(3,1,1)
    plt.plot(X[:-1],jacob[:,4,4,2,0],'.r',label='grad meth')   # dxf = -sin(x)
    plt.plot(X[:-1],good_jacob[:,4,4,2,0],'.k',label='good grad')
    plt.xticks(fontsize='15')
    plt.yticks(fontsize='15')
    plt.xlabel('X',fontsize='10')
    plt.legend(fontsize='15')

    plt.subplot(3,1,2)
    plt.plot(Y[:-1],jacob[4,:,4,2,1],'.r',label='grad meth')
    plt.plot(Y[:-1],good_jacob[4,:,4,2,1],'k',label='good grad')
    plt.xticks(fontsize='15')
    plt.yticks(fontsize='15')
    plt.xlabel('Y',fontsize='10')

    plt.subplot(3,1,3)
    plt.plot(Z[:-1],jacob[4,4,:,2,2],'.r',label='grad meth')
    plt.plot(Z[:-1],good_jacob[4,4,:,2,2],'k',label='good grad')
    plt.xticks(fontsize='15')
    plt.yticks(fontsize='15')
    plt.xlabel('Z',fontsize='10')

    

    plt.show()












if False :
    """ grad_vect"""
    def df(Coord):
        x, y, z = Coord[0], Coord[1], Coord[2]
        """
          f(x,y,z) = cos(x)+cos(y)+cos(z)
        """
        dxf = -np.sin(x)
        dyf = -np.sin(y/2)/2
        dzf = -np.sin(z/3)/3

        return(dxf,dyf,dzf)

    T = 2*np.pi
    X, Y, Z = np.linspace(-2,-2+T,20), np.linspace(-1,-1+2*T,21), np.linspace(0,0+3*T,22)
    coordX = ml.DataArrayDouble(X)
    coordY = ml.DataArrayDouble(Y)
    coordZ = ml.DataArrayDouble(Z)#X+Y); Z = X+Y
    cmesh = ml.MEDCouplingCMesh.New()
    cmesh.setCoords(coordX,coordY,coordZ)
    mesh = cmesh.buildUnstructured()
    
    champU = ml.MEDCouplingFieldDouble(ml.ON_CELLS)
    champU.setNature(ml.IntensiveMaximum)
    champU.setMesh(mesh)
    champU.fillFromAnalytic(3,  "1*(cos(x)+cos(y/2)+cos(z/3))*IVec +"+
                                "2*(cos(x)+cos(y/2)+cos(z/3))*JVec +"+
                                "3*(cos(x)+cos(y/2)+cos(z/3))*KVec")


    BC = 'perioX, perioY, perioZ'
    gradX = methodes_2020.grad_vect(champU,0,BC)
    gradY = methodes_2020.grad_vect(champU,1,BC)
    gradZ = methodes_2020.grad_vect(champU,2,BC)
    print(np.shape(gradX))
    
    good_grad = np.zeros(np.shape(gradX))
    xf = (X[0]+ X[-1])/2; X = (X[:-1]+X[1:])/2; X = np.concatenate((X,xf),axis = None)
    yf = (Y[0]+ Y[-1])/2; Y = (Y[:-1]+Y[1:])/2; Y = np.concatenate((Y,yf),axis = None)
    zf = (Z[0]+ Z[-1])/2; Z = (Z[:-1]+Z[1:])/2; Z = np.concatenate((Z,zf),axis = None)
    for i in range(len(X)-1):
        x = X[i]
        for j in range(len(Y)-1):
            y = Y[j]
            for k in range(len(Z)-1):
                z = Z[k]
                good_grad[i,j,k,:] = df([x,y,z])

    good_gradX =   good_grad
    good_gradY = 2*good_grad
    good_gradZ = 3*good_grad

    
    """ Pour dire si oui ou non norte gradient est juste on va s'appuyer
    sur la norme L2 de l'erreur """

    """ Si on s'interesse aux normes de certaines composantes, on en a
        pletor ici (grad_compo correspondent a ce qu'on doit avoir avec grad
        scal).
    """
    

    if False :
        grad_compoI = np.zeros(np.shape(gradX))
        grad_compoI[:,:,:,0] = gradX[:,:,:,0]
        grad_compoI[:,:,:,1] = gradY[:,:,:,0]
        grad_compoI[:,:,:,2] = gradZ[:,:,:,0]
    
        grad_compoJ = np.zeros(np.shape(gradX))
        grad_compoJ[:,:,:,0] = gradX[:,:,:,1]
        grad_compoJ[:,:,:,1] = gradY[:,:,:,1]
        grad_compoJ[:,:,:,2] = gradZ[:,:,:,1]
    
        grad_compoK = np.zeros(np.shape(gradX))
        grad_compoK[:,:,:,0] = gradX[:,:,:,2]
        grad_compoK[:,:,:,1] = gradY[:,:,:,2]
        grad_compoK[:,:,:,2] = gradZ[:,:,:,2]
    
    
        """idem mais pour les derivees analytiques"""
        L2_good_derX = np.linalg.norm(good_gradX )
        L2_good_derY = np.linalg.norm(good_gradY )
        L2_good_derZ = np.linalg.norm(good_gradZ )
    
        good_grad_compoI = np.zeros(np.shape(good_gradX))
        good_grad_compoI[:,:,:,0] = good_gradX[:,:,:,0]
        good_grad_compoI[:,:,:,1] = good_gradY[:,:,:,0]
        good_grad_compoI[:,:,:,2] = good_gradZ[:,:,:,0]
        L2_good_grad_compoI = np.linalg.norm(good_grad_compoI)
    
        good_grad_compoJ = np.zeros(np.shape(good_gradX))
        good_grad_compoJ[:,:,:,0] = good_gradX[:,:,:,1]
        good_grad_compoJ[:,:,:,1] = good_gradY[:,:,:,1]
        good_grad_compoJ[:,:,:,2] = good_gradZ[:,:,:,1]
        L2_good_grad_compoJ = np.linalg.norm(good_grad_compoJ)
    
        good_grad_compoK = np.zeros(np.shape(good_gradX))
        good_grad_compoK[:,:,:,0] = good_gradX[:,:,:,2]
        good_grad_compoK[:,:,:,1] = good_gradY[:,:,:,2]
        good_grad_compoK[:,:,:,2] = good_gradZ[:,:,:,2]
        L2_good_grad_compoK = np.linalg.norm(good_grad_compoK)
    
        L2_good_tot = np.linalg.norm(good_gradX+good_gradY+good_gradZ)
    
        """Et maintenant on fait les normes des erreurs"""
    
        L2_err_tot = np.linalg.norm(gradX+gradY+gradZ-(good_gradX+good_gradY+good_gradZ))
    
        L2_err_derX = np.linalg.norm(gradX-good_gradX)
        L2_err_derY = np.linalg.norm(gradY-good_gradY)
        L2_err_derZ = np.linalg.norm(gradZ-good_gradZ)
    
        L2_err_grad_compoI = np.linalg.norm(grad_compoI-good_grad_compoI)
        L2_err_grad_compoJ = np.linalg.norm(grad_compoJ-good_grad_compoJ)
        L2_err_grad_compoK = np.linalg.norm(grad_compoK-good_grad_compoK)
    
        """puis les pourcentages d'erreur."""
    
        pcTot = round((100*L2_err_tot/L2_good_tot                ),2)
        pcI =   round((100*L2_err_grad_compoI/L2_good_grad_compoI),2)
        pcJ =   round((100*L2_err_grad_compoJ/L2_good_grad_compoJ),2)
        pcK =   round((100*L2_err_grad_compoK/L2_good_grad_compoK),2)
        pcX =   round((100*L2_err_derX/L2_good_derX              ),2)
        pcY =   round((100*L2_err_derY/L2_good_derY              ),2)
        pcZ =   round((100*L2_err_derZ/L2_good_derZ              ),2)
    
        print('ERREURS : \t\t   (%)'+'\n \n'+
        'Totale : \t\t {0} \n'.format(pcTot)+'Composante I : \t\t {0} \n'.format(pcI)+
        'Composante J : \t\t {0} \n'.format(pcJ)+'Composante K : \t\t {0} \n'.format(pcK)+
        'Dérivée selon X : \t {0} \n'.format(pcX)+'Dérivée selon Y : \t {0} \n'.format(pcY)+
        'Dérivée selon Z : \t {0} \n'.format(pcZ))


    import matplotlib.pyplot as plt

    fig, axs = plt.subplots(3, sharex=True, sharey=True)
    ###   ###   ###
    plt.figure(1)
    plt.suptitle('Derivee par rapport a X',fontsize='15')

    plt.subplot(3,1,1)
    plt.plot(X[:-1],gradX[:,4,4,0],'.r',label='grad meth')   # dxf = -sin(x)
    plt.plot(X[:-1],good_gradX[:,4,4,0],'.k',label='good grad')
    plt.xticks(fontsize='15')
    plt.yticks(fontsize='15')
    plt.xlabel('X',fontsize='10')
    plt.legend(fontsize='15')

    plt.subplot(3,1,2)
    plt.plot(Y[:-1],gradX[4,:,4,1],'.r',label='grad meth')
    plt.plot(Y[:-1],good_gradX[4,:,4,1],'k',label='good grad')
    plt.xticks(fontsize='15')
    plt.yticks(fontsize='15')
    plt.xlabel('Y',fontsize='10')

    plt.subplot(3,1,3)
    plt.plot(Z[:-1],gradX[4,4,:,2],'.r',label='grad meth')
    plt.plot(Z[:-1],good_gradX[4,4,:,2],'k',label='good grad')
    plt.xticks(fontsize='15')
    plt.yticks(fontsize='15')
    plt.xlabel('Z',fontsize='10')

    fig, axs = plt.subplots(3, sharex=True, sharey=True)
    ###   ###   ###
    plt.figure(2,sharey=True)
    plt.suptitle('Derivee par rapport a Y',fontsize='15')

    plt.subplot(3,1,1)
    plt.plot(X[:-1],gradY[:,4,4,0],'.r',label='grad meth')
    plt.plot(X[:-1],good_gradY[:,4,4,0],'k',label='good grad')
    plt.xticks(fontsize='15')
    plt.yticks(fontsize='15')
    plt.xlabel('X',fontsize='10')
    plt.legend(fontsize='15')

    plt.subplot(3,1,2)
    plt.plot(Y[:-1],gradY[4,:,4,1],'.r',label='grad meth')
    plt.plot(Y[:-1],good_gradY[4,:,4,1],'k',label='good grad')
    plt.xticks(fontsize='15')
    plt.yticks(fontsize='15')
    plt.xlabel('Y',fontsize='10')

    plt.subplot(3,1,3)
    plt.plot(Z[:-1],gradY[4,4,:,2],'.r',label='grad meth')
    plt.plot(Z[:-1],good_gradY[4,4,:,2],'k',label='good grad')
    plt.xticks(fontsize='15')
    plt.yticks(fontsize='15')
    plt.xlabel('Z',fontsize='10')


    fig, axs = plt.subplots(3, sharex=True, sharey=True)
    ###   ###   ###
    plt.figure(3,sharey=True)
    plt.suptitle('Derivee par rapport a Z',fontsize='15')

    plt.subplot(3,1,1)
    plt.plot(X[:-1],gradZ[:,4,4,0],'.r',label='grad meth')
    plt.plot(X[:-1],good_gradZ[:,4,4,0],'k',label='good grad')
    plt.xticks(fontsize='15')
    plt.yticks(fontsize='15')
    plt.xlabel('X',fontsize='10')
    plt.legend(fontsize='15')

    plt.subplot(3,1,2)
    plt.plot(Y[:-1],gradZ[4,:,4,1],'.r',label='grad meth')
    plt.plot(Y[:-1],good_gradZ[4,:,4,1],'k',label='good grad')
    plt.xticks(fontsize='15')
    plt.yticks(fontsize='15')
    plt.xlabel('Y',fontsize='10')

    plt.subplot(3,1,3)
    plt.plot(Z[:-1],gradZ[4,4,:,2],'.r',label='grad meth')
    plt.plot(Z[:-1],good_gradZ[4,4,:,2],'k',label='good grad')
    plt.xticks(fontsize='15')
    plt.yticks(fontsize='15')
    plt.xlabel('Z',fontsize='10')



    
    plt.show()




if False :
    """ grad_scal"""

    def df(Coord):
        x, y, z = Coord[0], Coord[1], Coord[2]
        """
          f(x,y,z) = cos(x)+cos(y)+cos(z)
        """
        dxf = -np.sin(x)
        dyf = -np.sin(y)
        dzf = -np.sin(z)

        return(dxf,dyf,dzf)

    T = 2*np.pi
    X, Y, Z = np.linspace(-2,-2+T,20), np.linspace(-1,-1+T,21), np.linspace(0,0+T,22)
    coordX = ml.DataArrayDouble(X)
    coordY = ml.DataArrayDouble(Y)
    coordZ = ml.DataArrayDouble(Z)#X+Y); Z = X+Y
    cmesh = ml.MEDCouplingCMesh.New()
    cmesh.setCoords(coordX,coordY,coordZ)
    mesh = cmesh.buildUnstructured()
    
    champU = ml.MEDCouplingFieldDouble(ml.ON_CELLS)
    champU.setNature(ml.IntensiveMaximum)
    champU.setMesh(mesh)
    champU.fillFromAnalytic(1, "cos(x)+cos(y)+cos(z)")


    BC = 'perioX, perioY, perioZ'
    grad, MEDgrad = methodes_2020.grad_scal(champU,BC)
    
    good_grad = np.zeros(np.shape(grad))
    xf = (X[0]+ X[-1])/2; X = (X[:-1]+X[1:])/2; X = np.concatenate((X,xf),axis = None)
    yf = (Y[0]+ Y[-1])/2; Y = (Y[:-1]+Y[1:])/2; Y = np.concatenate((Y,yf),axis = None)
    zf = (Z[0]+ Z[-1])/2; Z = (Z[:-1]+Z[1:])/2; Z = np.concatenate((Z,zf),axis = None)
    for i in range(len(X)-1):
        x = X[i]
        for j in range(len(Y)-1):
            y = Y[j]
            for k in range(len(Z)-1):
                z = Z[k]
                good_grad[i,j,k,:] = df([x,y,z])

    """ Pour dire si oui ou non norte gradient est juste on va s'appuyer
    sur la norme L2 de l'erreur """
    L2_good = np.linalg.norm(np.abs(good_grad));
    L2_grad = np.linalg.norm(np.abs(grad)); 
    L2_err  = np.linalg.norm(np.abs(grad-good_grad));

    import matplotlib.pyplot as plt

    plt.figure(1)
    plt.title('selon l\'axe Z',fontsize='15')
    plt.plot(Z[:-1],grad[4,4,:,2],'r',label='grad meth')
    plt.plot(Z[:-1],good_grad[0,0,:,2],'k',label='good grad')
    plt.xlabel('z',fontsize='25')
    plt.xticks(fontsize='15')
    plt.yticks(fontsize='15')
    plt.legend(fontsize='15')


    plt.figure(2)
    plt.title('selon l\'axe X',fontsize='15')
    plt.plot(X[:-1],grad[:,4,4,0],'r',label='grad meth')
    plt.plot(X[:-1],good_grad[:,4,4,0],'k',label='good grad')
    plt.xlabel('x',fontsize='25')
    plt.xticks(fontsize='15')
    plt.yticks(fontsize='15')
    plt.legend(fontsize='15')

    plt.figure(3)
    plt.title('selon l\'axe Y',fontsize='15')
    plt.plot(Y[:-1],grad[4,:,4,1],'r',label='grad meth')
    plt.plot(Y[:-1],good_grad[4,:,4,1],'k',label='good grad')
    plt.xlabel('y',fontsize='25')
    plt.xticks(fontsize='15')
    plt.yticks(fontsize='15')
    plt.legend(fontsize='15')
    
    plt.show()
        
        
    



if False :
    """ UMeshToCMesh"""
    coord1 = np.logspace(-3,1,5); coord2 = np.linspace(-2,1,5)
    coordX = ml.DataArrayDouble(coord1); print('coordX'); print(coordX)
    coordY = ml.DataArrayDouble(coord2); print('coordY'); print(coordY)
    coordZ = ml.DataArrayDouble(coord1+coord2); print('coordZ'); print(coordZ)
    cmesh = ml.MEDCouplingCMesh.New()
    cmesh.setCoords(coordX,coordY,coordZ)
    mesh = cmesh.buildUnstructured()
    
    champU = ml.MEDCouplingFieldDouble(ml.ON_CELLS)
    champU.setNature(ml.IntensiveMaximum)
    champU.setMesh(mesh)
    champU.fillFromAnalytic(2, "2*x*y*IVec+6*y*z*JVec")
    
    
    mC, idx = methodes_2020.fromUtoC(champU.getMesh())
    champC = methodes_2020.UMeshToCMesh(champU,mC,idx)
    
    point1 = [1e-2, 0.7, -0.3]
    point2 = [0.8, -1, 8]
    point3 = [1, -0.1, 3.4]
    
    eq1 = (champC.getValueOn(point1) == champU.getValueOn(point1))
    eq2 = (champC.getValueOn(point2) == champU.getValueOn(point2))
    eq3 = (champC.getValueOn(point3) == champU.getValueOn(point3))
    
    print(eq1,eq2,eq3,eq1 and eq2 and eq3)

if False:
    """ toldim et from U to C """
    res = methodes_2020.get_toldim(N_mesh)
    mesh = methodes_2020.fromUtoC(mesh)
    
    coord1 = np.logspace(-2,1,4); coord2 = np.linspace(-2,1,4)
    coordX = ml.DataArrayDouble(coord1)
    coordY = ml.DataArrayDouble(coord2)
    coordZ = ml.DataArrayDouble(coord1*coord2)
    cmesh = ml.MEDCouplingCMesh.New()
    cmesh.setCoords(coordX,coordY,coordZ)
    mesh = cmesh.buildUnstructured()
    
    """ toldim et grd_scal et getCentreCoords"""
    toldim = methodes_2020.get_toldim(mesh,1)
    methodes_2020.grad_scal(my_dudx)
    centre = methodes_2020.getCentreCoords(mesh)


print("end")


