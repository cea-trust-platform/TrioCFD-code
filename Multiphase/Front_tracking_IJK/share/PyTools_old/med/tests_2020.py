# -*- coding: utf-8 -*-
"""
https://docs.python.org/fr/3/library/unittest.html
"""
import rlcompleter, readline
readline.parse_and_bind('tab:complete')
#############################""
import unittest
# import Tools
import numpy as np
import MEDLoader as ml
import methodes_2020
import outils_physique


class GuillaumeTest(unittest.TestCase):  
    def setUp(self):
        pass  

    def tearDown(self):
        pass

    def testXYZ(self):
        a = 1
        b = 2
        f = True
        self.assertEqual(3, a+b, 'equal')
        self.assertTrue(f, 'true')
        #self.assertRaises(ValueError, self.toto, 1, 2)
        return
    
    def testABC(self):
        return
      
    def test_en_tete(self):

        succes = 1
        message = ''
        self.assertTrue(succes,message)

    def test_MEDFileToMEDField_MultiTS(self):

        T = 1*np.pi
        X = np.linspace(-2,-2+T,10);        Nx = len(X)
        coordX,coordY,coordZ = ml.DataArrayDouble(X),ml.DataArrayDouble(X),ml.DataArrayDouble(X)
        cmesh = ml.MEDCouplingCMesh.New()
        cmesh.setCoords(coordX,coordY,coordZ)
        mesh = cmesh.buildUnstructured()
        mesh.setName("Mesh")
    
        ml.WriteUMesh("test.med",mesh,True)
    
        iterations = [(i,-1) for i in range(11)] # pour avoir jusqu'à 10 inclus
    
        MEDchamp = ml.MEDCouplingFieldDouble.New(ml.ON_CELLS,ml.ONE_TIME)
        MEDchamp.setName("Field");          MEDchamp.setMesh(mesh)
        for it in iterations:
            w = 2*np.pi/10
            valeur = np.zeros((Nx-1,Nx-1,Nx-1,3))
            for i,x in enumerate(X[:-1]):
                for j,y in enumerate(X[:-1]):
                    for k,z in enumerate(X[:-1]):
                        valeur[i,j,k,0] = 3+np.cos(w*it[0])
                        valeur[i,j,k,1] = np.cos(x)+np.cos(w*it[0])
                        valeur[i,j,k,2] = it[0]+np.cos(w*it[0])
            valeur = valeur.transpose(2,1,0,3);
            valeur = valeur.flatten()
            champ = ml.DataArrayDouble(valeur)
            champ.rearrange(3)
            MEDchamp.setArray(champ)
            MEDchamp.setTime(it[0],it[0],it[1])
            ml.WriteFieldUsingAlreadyWrittenMesh("test.med",MEDchamp)

        Ffield,Fmesh,TypeOfField = methodes_2020.MEDFileToMEDField_MultiTS("test.med","Field")

        """condition de succes : les valeurs de  dernier champ du .med correspondent
        à celles du dernier champ généré. C'est léger, mais ça le fait disons """
        testchamp = Ffield.getFieldOnMeshAtLevel(
                        TypeOfField,
                        iterations[-1][0],iterations[-1][1],
                        0,Fmesh,0).getArray()
        succes = (testchamp.getValues()==champ.getValues())
        message = 'mauvaise récupération lecture du fichier .med'
        self.assertTrue(succes,message)


    def test_MEDFileToMEDField_1TS(self):
        """On ne se sert jamais de cette fonction pour le moment,
           on la testera quand elle sera utile...                """
        succes = 1
        message = ''
        self.assertTrue(succes,message)


    def test_MEDFieldToNumpy(self):

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

        f,mC = methodes_2020.en_tete(MEDchamp)
        xn,yn,zn,NPchamp = methodes_2020.MEDFieldToNumpy(f,mC)

        dx = X-xn.getValues();  dy = Y-yn.getValues(); dz = Z-zn.getValues()
        lx = np.linalg.norm(dx); ly = np.linalg.norm(dy) ; lz = np.linalg.norm(dz)
        dchamp = valeur-NPchamp
        lc = np.linalg.norm(dchamp)
        succes = lx+ly+lz+lc 
        succes = succes==0

        message = ("\nMauvaise converion med->numpy. 1 pour OK. 0 sinon : \n"+
           ' X \t {0} \n Y \t {1} \n Z \t {2} \n'.format(int(lx==0),int(ly==0),int(lz==0))+
           'Champ \t {0} \n (ATTN à la transposition : x,y,z --> z,y,x) '.format(int(lc==0)))

        self.assertTrue(succes,message)

    def test_MEDFieldToNumpy_SansChamp(self):

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
        ####################################
        # Pas besoin de remplir le champ...
        ####################################
        mod_valeur = valeur.transpose(2,1,0,3);
        mod_valeur = mod_valeur.flatten()
        champ = ml.DataArrayDouble(mod_valeur)
        champ.rearrange(3)
        MEDchamp.setArray(champ)
        MEDchamp.setTime(0,0,-1)

        f,mC = methodes_2020.en_tete(MEDchamp)
        xn,yn,zn = methodes_2020.MEDFieldToNumpy_SansChamp(mC)

        dx = X-xn.getValues();  dy = Y-yn.getValues(); dz = Z-zn.getValues()
        lx = np.linalg.norm(dx); ly = np.linalg.norm(dy) ; lz = np.linalg.norm(dz)
        succes = lx+ly+lz
        print(succes)
        succes = succes==0

        message = ("\nMauvaise converion med->numpy. 1 pour OK. 0 sinon : \n"+
           ' X \t {0} \n Y \t {1} \n Z \t {2} \n'.format(int(lx==0),int(ly==0),int(lz==0))
           )

        self.assertTrue(succes,message)


    def test_NumpyToMEDField(self):

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

        testMEDchamp = methodes_2020.NumpyToMEDField(valeur,cmesh,[0,0,-1],"Field")

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


    def test_NumpyToMEDField_VRij(self):

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

        testMEDchamp = methodes_2020.NumpyToMEDField_VRij(valeur,MEDchamp,cmesh,[0,0,-1])
        
        MEDchamp.setArray(champ)
        MEDchamp.setTime(0,0,-1)

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



    def test_MEDFieldToMEDFile(self):

        T = 1*np.pi
        X = np.linspace(-2,-2+T,10);        Nx = len(X)
        Y = np.linspace(-1,-1+T,11);        Ny = len(Y)
        Z = np.linspace(-0,-0+T,12);        Nz = len(Z)
        coordX,coordY,coordZ = ml.DataArrayDouble(X),ml.DataArrayDouble(Y),ml.DataArrayDouble(Z)
        cmesh = ml.MEDCouplingCMesh.New()
        cmesh.setCoords(coordX,coordY,coordZ)
        mesh = cmesh.buildUnstructured()
        mesh.setName("Mesh")

        ml.WriteUMesh("good.med",mesh,True)
        iterations = [(i,-1) for i in range(5)] # pour avoir jusqu'à 4 inclu
    
        MEDchamp = ml.MEDCouplingFieldDouble.New(ml.ON_CELLS,ml.ONE_TIME)
        MEDchamp.setName("Field");
        MEDchamp.setMesh(mesh)
        for t,it in enumerate(iterations):
            w = 2*np.pi/4
            valeur = np.zeros((Nx-1,Ny-1,Nz-1,3))
            for i,x in enumerate(X[:-1]):
                for j,y in enumerate(Y[:-1]):
                    for k,z in enumerate(Z[:-1]):
                        valeur[i,j,k,0] = 10*np.cos(x)  +np.cos(w*it[0])
                        valeur[i,j,k,1] = np.cos(z)     +np.cos(w*it[0])
                        valeur[i,j,k,2] = np.cos(10*y)  +np.cos(w*it[0])
            mod_valeur = valeur.transpose(2,1,0,3);
            mod_valeur = mod_valeur.flatten()
            champ = ml.DataArrayDouble(mod_valeur)
            champ.rearrange(3)
            MEDchamp.setArray(champ)
            MEDchamp.setTime(it[0],it[0],it[1])
            if t==0:
                methodes_2020.MEDFieldToMEDFile(MEDchamp,"test.med",True,False)  
            ml.WriteFieldUsingAlreadyWrittenMesh("good.med",MEDchamp)
            methodes_2020.MEDFieldToMEDFile(MEDchamp,"test.med",False,False)

        """
            En soi, les 2 binaires (.med) sont différents. (commande diff du terminal)
            Ce qui nous intéresse nous, c'est leur CONTENU. Pour y acceder, on va se
            servir des fonctions testées plus haut. C'est pas propre mais ça fait le
            travail on espère.

            Nom des variables :
                 - g : good         - t : test
                 - f : file         - m : mesh
                 - F : File
                 - TOF:  TypeOfFile
                 - V : Value        - C : Coords        - N : Name
                 - Its : Iterations
                 - s : succes
                 
        """
        gFf,gFm,gTOF = methodes_2020.MEDFileToMEDField_MultiTS("good.med","Field")
        tFf,tFm,tTOF = methodes_2020.MEDFileToMEDField_MultiTS("test.med","Field")


        # Iterations
        gIts = gFf.getIterations();        tIts = tFf.getIterations()
        # maillage
        gm   = gFm.getMeshAtLevel(0);      tm   = tFm.getMeshAtLevel(0)
        gmC  = gm.getCoords().getValues(); tmC  = tm.getCoords().getValues()    
        gmN  = gm.getName();               tmN  = tm.getName()
        # champ de valeurs
        gf   = gFf.getFieldOnMeshAtLevel(gTOF,gIts[1][0],gIts[1][1],0,gFm,0)
        tf   = tFf.getFieldOnMeshAtLevel(tTOF,tIts[1][0],tIts[1][1],0,tFm,0)
        gfV  = gf.getArray().getValues(); tfV  = tf.getArray().getValues()
        gfN  = gf.getName();              tfN  = tf.getName()
        # Succes
        sTOF = (gTOF==tTOF)
        sIts = (gIts==tIts)
        smN  = (gmN==tmN);             smC = (gmC==tmC)
        sfN  = (gfN==tfN);             sfV = (gfV==tfV)

        succes = (sTOF and sIts and smN and smC and sfN and sfV )
        message = ('\nMauvaise ecriture du champ dans le fichier (field--> file)\n'+
                   '--> 1 pour OK, 0 sinon\n'+ 
                   '  Maillage :\n\tNoms \t {0} \n\tCoords \t {1}\n'.format(int(smN),int(smC))+
                   '  Champ :\n\tNoms \t {0} \n\tVals \t {1}'.format(int(sfN),int(sfV))
                   )
        self.assertTrue(succes,message)

      
    def test_get_toldim(self):
        coord1 = np.logspace(-2,1,4); coord2 = np.linspace(-2,1,4)
        coordX = ml.DataArrayDouble(coord1*2)
        coordY = ml.DataArrayDouble(coord2/2)
        coordZ = ml.DataArrayDouble(coord1*coord2)
        cmesh = ml.MEDCouplingCMesh.New()
        cmesh.setCoords(coordX,coordY,coordZ)
        mesh = cmesh.buildUnstructured()
      
        toldim = np.array(methodes_2020.get_toldim(mesh,1))
        goodtoldim = np.array([0.18000000000000002, 0.5, 0.08])
        
        diff = np.sum(np.abs(goodtoldim-toldim))
        self.assertEqual(diff, 0, 'get_toldim calcule mal.')
        return
      
      
      
    def test_fromUtoC(self):
        # Le paramètre tol est reglé pour le cas du maillage CFD à 3088800 de mailles
        # Test de la fonction fromUtoC
        """
        @param mU a MEDCouplingUMesh
        """
        coord1 = np.linspace(-2,1,4)
        coordX = ml.DataArrayDouble(coord1*2)
        coordY = ml.DataArrayDouble(coord1/2)
        coordZ = ml.DataArrayDouble(coord1/4)
        cmesh = ml.MEDCouplingCMesh.New()
        cmesh.setCoords(coordX,coordY,coordZ)
        
        mU = cmesh.buildUnstructured()
        mC, idx = methodes_2020.fromUtoC(mU)
        
        # Test cell renumbering
        bU = mU.computeCellCenterOfMass()
        bC = mC.computeCellCenterOfMass()
        diff = bU[idx] - bC
        
        Ndim = mU.getMeshDimension()
        diff.rearrange(Ndim)


        ret = ml.MEDCouplingCMesh.New()
        ret.setName(mU.getName())
        axis=mU.getSpaceDimension()*[None]
        tol = methodes_2020.get_toldim(mU)
        
        matching = True
        for i in range(Ndim):
            toli = tol[i]
            axis[i]=mU.getCoords()[:,i].getDifferentValues(toli)
            axis[i].sort()
            
            matching = matching and diff[:,0].isUniform(0.0, tol[0])
            
            pass
        
        ret.setCoords(*tuple(axis)) 
        nc = ret.getNumberOfCells() 
        
        baryC = ret.computeCellCenterOfMass()
        baryU = mU.computeCellCenterOfMass()
        bary=ml.DataArrayDouble.Aggregate(baryC,baryU)
        tolbary = 0.5*min(methodes_2020.get_toldim(mU,1.))
        c,cI=bary.findCommonTuples(tolbary)
        # c  : numero des cellules qui ont les mm coordonnees à tolbary pres.
        # cI : de c[cI[i]] a c[cI[i+1]-1] on parcours les c qui ont les mm coordonnees.
        
        
        # Serie de tests pour s'assurer de la bonne association des barycentres correspondants
        assert ( c.getNumberOfTuples() == nc * 2 )
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
        # Fin des tests
    
        
        # Have we the proper number of cells ?
        assert (nc == mU.getNumberOfCells()) 
        # Are the barycenter matching ?
        print "Are the barycenters matching ? ", matching
        if not matching:
            print "Max : ", diff.getMaxValue()
            print "Min : ", diff.getMinValue()
            raise Exception("Barycenters not matching")
            pass
      
    def test_getCentreCoords(self):
        coord1 = np.logspace(-2,1,4); coord2 = np.linspace(-2,1,4)
        coordX = ml.DataArrayDouble(coord1)
        coordY = ml.DataArrayDouble(coord2)
        coordZ = ml.DataArrayDouble(coord1*coord2)
        cmesh = ml.MEDCouplingCMesh.New()
        cmesh.setCoords(coordX,coordY,coordZ)
        mesh = cmesh.buildUnstructured()

        CC = methodes_2020.getCentreCoords(mesh)
        goodCC = tuple((5.005, -0.5, 4.95))

        match = True
        for i in range(3):
            match = match and (CC[i] == goodCC[i])

        self.assertTrue(match,
        'MESH centers not corresponding')
        pass

    def test_UMeshToCMesh(self):

        coord1 = np.logspace(-3,1,5); coord2 = np.linspace(-2,1,5)
        coordX = ml.DataArrayDouble(coord1);
        coordY = ml.DataArrayDouble(coord2);
        coordZ = ml.DataArrayDouble(coord1+coord2);
        cmesh = ml.MEDCouplingCMesh.New()
        cmesh.setCoords(coordX,coordY,coordZ)
        mesh = cmesh.buildUnstructured()
        
        champU = ml.MEDCouplingFieldDouble(ml.ON_CELLS)
        champU.setNature(ml.IntensiveMaximum)
        champU.setMesh(mesh)
        champU.fillFromAnalytic(2, "2*x*y*IVec+6*y*z*JVec")
        
        mC, idx = methodes_2020.fromUtoC(champU.getMesh())
        champC = methodes_2020.UMeshToCMesh(champU,mC,idx)

        """consider 3 random laocation in the mesh, get the value
        at those locations from champC and champU. The values should
        be strictly equal."""

        point1 = [1e-2, 0.7, -0.3]
        point2 = [0.8, -1, 8]
        point3 = [1, -0.1, 3.4]
        
        eq1 = (champC.getValueOn(point1) == champU.getValueOn(point1))
        eq2 = (champC.getValueOn(point2) == champU.getValueOn(point2))
        eq3 = (champC.getValueOn(point3) == champU.getValueOn(point3))

        EQ = (eq1 and eq2 and eq3)
        self.assertTrue(EQ,
        '\n The values fetched on the Cartesian do not match the Unstructured ones.')
        pass

    def test_grad_scal(self):
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

        """
            1. on va tester le gradient en conditions périodiques
            ATTN : on  prend garde à ce que les champs et "les bords"
            soient bien periodoqies
        """
        champU.fillFromAnalytic(1, "cos(x)+cos(y)+cos(z)")
        BC = 'perioX, perioY, perioZ'
        grad, MEDgrad = outils_physique.grad_scal(champU,BC)

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
                    good_grad[i,j,k,:] = -np.sin([x, y, z])

        """
            Pour dire si oui ou non notre gradient est ok, on va s'ap-
            puyer sur la norme L2 de l'erreur.
        """
        L2good = np.linalg.norm(np.abs(good_grad))
        L2err  = np.linalg.norm(np.abs(grad-good_grad))
        succes = L2err/L2good <= 0.1
        self.assertTrue(succes,'PERIODIQUE : La norme de l\'erreur est de plus de 10 % de la norme du gradient ({0} %)'.format(str(L2err/L2good*100)))

        """
            2. on va tester le gradient en conditions non-periodiques
        """
            # AFFAIRE A SUIVRE

    def test_grad_vect(self):
    
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
        champU.fillFromAnalytic(3,  "1*(cos(x)+cos(y)+cos(z))*IVec +"+
                                    "2*(cos(x)+cos(y)+cos(z))*JVec +"+
                                    "3*(cos(x)+cos(y)+cos(z))*KVec")
    
    
        BC = 'perioX, perioY, perioZ'
        gradX = outils_physique.grad_vect(champU,0,BC)
        gradY = outils_physique.grad_vect(champU,1,BC)
        gradZ = outils_physique.grad_vect(champU,2,BC)
        
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
    
        grad_compoI = np.zeros(np.shape(gradX))
        grad_compoI[:,:,:,0] = gradX[:,:,:,0];
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

        message = ('ERREURS : \t\t   (%)'+'\n \n'+
        'Totale : \t\t {0} \n'.format(pcTot)+'Composante I : \t\t {0} \n'.format(pcI)+
        'Composante J : \t\t {0} \n'.format(pcJ)+'Composante K : \t\t {0} \n'.format(pcK)+
        'Dérivée selon X : \t {0} \n'.format(pcX)+'Dérivée selon Y : \t {0} \n'.format(pcY)+
        'Dérivée selon Z : \t {0} \n'.format(pcZ)+
        '\n IL FAUT QUE LES ERREURS SOIENT DE MOINS DE 10% POUR ETRE VALIDE')

        succes = (pcTot<=10 and pcI<=10 and pcJ<=10 and pcK<=10
                    and pcX<=10 and pcY<=10 and pcZ<=10)
                    
        self.assertTrue(succes,'\n'+message)

    def test_jacob(self):

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
        jacob, MEDjacob = outils_physique.jacob(champU,BC)
        # gradX = outils_physique.grad_vect(champU,0,BC)
        # gradY = outils_physique.grad_vect(champU,1,BC)
        # gradZ = outils_physique.grad_vect(champU,2,BC)
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
        'Dérivée selon Z : \t {0} \n'.format(pcZ)+
        '\n IL FAUT QUE LES ERREURS SOIENT DE MOINS DE 10% POUR ETRE VALIDE')
        succes = (pcT<=10 and pcI<=10 and pcJ<=10 and pcK<=10 and
                pcX<=10 and pcY<=10 and pcZ<=10)
                
        self.assertTrue(succes,message)

    def test_get_mini_maxi(self):
        T = 2*np.pi
        X, Y, Z = np.linspace(-2,-2+T,20), np.linspace(-1,-1+2*T,21), np.linspace(0,30,31)
        coordX = ml.DataArrayDouble(X)
        coordY = ml.DataArrayDouble(Y)
        coordZ = ml.DataArrayDouble(Z)#X+Y); Z = X+Y
        cmesh = ml.MEDCouplingCMesh.New()
        cmesh.setCoords(coordX,coordY,coordZ)
        mesh = cmesh.buildUnstructured()
    
        delta = 5
        # Récupération des corrdonnées des noeuds du maillage
        noeudz = cmesh.getCoordsAt(2).toNumPyArray()
        # Récupération des coordonnées des centre des cellules
        baryz = (noeudz[1:] + noeudz[:-1])/2

        mini_z, maxi_z = methodes_2020.get_mini_maxi(baryz,delta)

        good_mini_z = [0, 0, 0, 0, 0 ]+[i for i in range(25)]
        good_maxi_z = [i for i in range(4,29)]+[29, 29, 29, 29, 29, 29]

        message_mini = ("probleme dans la construction de mini_z : \n nous : {0},\n bon : {1}"
        .format(mini_z,good_mini_z))
        succes_mini = True
        for (i,miz) in enumerate(mini_z):
            succes_mini = (succes_mini and miz==good_mini_z[i])


        message_maxi = ("probleme dans la construction de maxi_z : \n nous : {0},\n bon : {1}"
        .format(maxi_z,good_maxi_z))    
        succes_maxi = True
        for (i,maz) in enumerate(maxi_z):
            succes_maxi = (succes_maxi and maz==good_maxi_z[i])
    
        self.assertTrue(succes_mini,message_mini)
        self.assertTrue(succes_maxi,message_maxi)
        pass


    def test_filtre_gab(self):

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
        filtredemi = outils_physique.filtre_gab(delta/2,TOTchampScalNmpy,tailles,BC)
    
        L2main = np.linalg.norm(np.abs(MainchampScalNmpy))
        errTOT = filtredemi-MainchampScalNmpy; L2errTOT = np.linalg.norm(errTOT)
    
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
    
        succes = (pcT<=10 and pcX<=10 and pcY<=10 and pcZ<=10)
        self.assertTrue(succes,message)
        pass



    def test_TimeAvg(self):
        """ creation d'un champ sur plusieurs TS"""
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
    
        ml.WriteUMesh("test_TA.med",mesh,True)
    
        iterations = [(i,-1) for i in range(61)]
    
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
            ml.WriteFieldUsingAlreadyWrittenMesh("test_TA.med",MEDchamp)
            pass

        champMean,MEDMean = outils_physique.TimeAvg('test_TA.med','Field')

        """On evalue les normes L2 des erreurs"""
        goodMean = np.zeros((Nx-1,Ny-1,Nz-1,3))
        goodMean[...,0] = 3                     # <3     +cos(wt)>_[0;60]
        for i,x in enumerate(X[:-1]):
            goodMean[i,:,:,1] = np.cos(x)       # <cos(x)+cos(wt)>_[0;60]
        goodMean[...,2] = 30                    # <t     +cos(wt)>_[0;60]
        
        err = goodMean-champMean
        
        L2errT = np.linalg.norm(err);       L2T = np.linalg.norm(goodMean)
        L2errX = np.linalg.norm(err[...,0]);L2X = np.linalg.norm(goodMean[...,0])
        L2errY = np.linalg.norm(err[...,1]);L2Y = np.linalg.norm(goodMean[...,1])
        L2errZ = np.linalg.norm(err[...,2]);L2Z = np.linalg.norm(goodMean[...,2])
        
        pcT = round(100*L2errT / L2T , 5)
        pcX = round(100*L2errX / L2X , 5)
        pcY = round(100*L2errY / L2Y , 5)
        pcZ = round(100*L2errZ / L2Z , 5)
        
        message = ('ERREUR DANS LA MOYENNE TEMPORELLE \t\t (%)'+'\n \n'+
        'Totale : \t\t {0} \n'.format(pcT)+ 'Selon X : \t\t {0} \n'.format(pcX)+
        'Selon Y : \t\t {0} \n'.format(pcY)+'Selon Z : \t\t {0} \n'.format(pcZ))
    
        succes = (pcT<=10 and pcX<=10 and pcY<=10 and pcZ<=10)
        self.assertTrue(succes,message)



    def test_Applatit(self):
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
        # Remarque : on pourrait définir un "goodchamPlat" et faire la norme L2 de
        # l'erreur. Ici ce "goodchamp" doit etre remplit de 0, et soustraire 0 à qqch
        # c'est comme ne rien faire. On se passera donc de ce "goodchamp".
        ##########################################################################
        l2A = round(100*np.linalg.norm(NPchampPlatX[0,:,:,0]),3)
        l2B = round(100*np.linalg.norm(NPchampPlatY[:,0,:,1]),3)
        l2C = round(100*np.linalg.norm(NPchampPlatZ[:,:,0,2]),3)
    
        succesA = l2A<5;    succesB = l2B<5;    succesC = l2C<5
        succes = (succesA and succesB and succesC)
    
        message = ('Problème dans la moyene par plan, exigence à 5% près :\n'+
                    '--> 1 pour OK, 0 sinon\n'+
                    'Moyenne selon X (résultat dans le plan Y,Z)\t{0}\n'.format(int(succesA))+
                    'Moyenne selon Y (résultat dans le plan X,Z)\t{0}\n'.format(int(succesB))+
                    'Moyenne selon Z (résultat dans le plan X,Y)\t{0}\n'.format(int(succesC))
                    )
        self.assertTrue(succes,message)   
        pass



    def test_Rij(self):
    
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

        succes = (succesTT and succesXX and succesXY and succesXZ and
                               succesYX and succesYY and succesYZ and
                               succesZX and succesZY and succesZZ)
        message = ('ERREUR DANS LES TENSIONS DE REYNOLDS \t\t '+
        '\nErreurs relatives en % (erreur max acceptée : 10%) :\n'+
        'Totale : \t\t {0} \n'.format(pcTT)+
        'Selon X : \t\t {0}, {1}, {2} \n'.format(pcXX,pcXY,pcXZ)+
        'Selon Y : \t\t {0}, {1}, {2} \n'.format(pcYX,pcYY,pcYZ)+
        'Selon Z : \t\t {0},  {1}, {2}\n'.format(pcZX,pcZY,pcZZ))

        self.assertTrue(succes,message)
        
    # def assertNumpyEqual(self, a, b):
        # diff = np.sum(np.abs(a-b)) 
        # self.assertEqual(diff, 0, 'hohoho')
        # return diff < 1.0e-6
  
    # def testScalarFieldFromFile(self):
        # x = np.arange(0,1,0.1)
        # y = x*x
        # np.savetxt("test.dt_ev", np.c_[x,y], header="x y")
        # f = Tools.Field.LoadFromFile("test.dt_ev", ["y"])
        # self.assertTrue(self.assertNumpyEqual(y, f._npa))
        # return


def suite():
    suite = unittest.TestSuite()
    # suite.addTest(GuillaumeTest('testXYZ'))
    # suite.addTest(GuillaumeTest('testABC'))
    # suite.addTest(GuillaumeTest('test_en_tete'))
    # suite.addTest(GuillaumeTest('test_MEDFileToMEDField_MultiTS'))
    # suite.addTest(GuillaumeTest('test_MEDFileToMEDField_1TS'))
    # suite.addTest(GuillaumeTest('test_MEDFieldToNumpy'))
    # suite.addTest(GuillaumeTest('test_MEDFieldToNumpy_SansChamp'))
    # suite.addTest(GuillaumeTest('test_NumpyToMEDField'))
    # suite.addTest(GuillaumeTest('test_NumpyToMEDField_VRij'))
    # suite.addTest(GuillaumeTest('test_MEDFieldToMEDFile'))
    # suite.addTest(GuillaumeTest('test_get_toldim'))
    # suite.addTest(GuillaumeTest('test_fromUtoC'))
    # suite.addTest(GuillaumeTest('test_getCentreCoords'))
    # suite.addTest(GuillaumeTest('test_UMeshToCMesh'))
    # suite.addTest(GuillaumeTest('test_grad_scal'))
    # suite.addTest(GuillaumeTest('test_grad_vect'))
    # suite.addTest(GuillaumeTest('test_jacob'))
    # suite.addTest(GuillaumeTest('test_get_mini_maxi'))
    # suite.addTest(GuillaumeTest('test_filtre_gab'))
    # suite.addTest(GuillaumeTest('test_TimeAvg'))
    # suite.addTest(GuillaumeTest('test_Applatit'))
    suite.addTest(GuillaumeTest('test_Rij'))


    
    # ~ suite.addTest(GuillaumeTest('testScalarFieldFromFile'))
    return suite

if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    runner.run(suite())
