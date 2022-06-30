#!visit -nowin -cli
#import numpy as np
import os, glob, sys, re

def getValue(key, fic, pre="^", default=None):
    f = open(fic, 'r')
    lines = f.readlines()
    f.close()
    nb = len(lines)
    rc = re.compile(pre + "\s*" + key + "\s*(?P<value>[\-]?[\d]*.?[\d]*[eE]?[+\-]?[\d]*)")
    for i, st in enumerate(lines):
        m = rc.match(st)
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
    nb = len(lines)
    rc = re.compile(pre + key + "\s*(?P<value>[\-]?[\d]*.?[\d]*[eE]?[+\-]?[\d]*)")
    for i, st in enumerate(lines):
        it = re.finditer(rc, st)
        val = []
        for match in it:
            val.append(float(match.group("value")))
        if len(val):
            return np.array(val)
    print("On a rien trouve dans ", fic)
    raise Exception("Etonnant, non?")
    return np.array([])

# v=sys.argv[-1]
# if (v not in ["P", "U", "V", "W"]): raise Exception("Option should be given in P, U, V or W")

def myClip2(z,nz, Lx=0.018):
   AddOperator("Clip",1)
   ClipAtts2 = ClipAttributes()
   ClipAtts2.quality = ClipAtts2.Accurate  # Fast, Accurate
   ClipAtts2.funcType = ClipAtts2.Plane  # Plane, Sphere
   ClipAtts2.plane1Status = 0
   ClipAtts2.plane2Status = 0
   #
   ClipAtts2.plane2Status = 1
   ClipAtts2.plane2Origin = (Lx, 0, 0)
   ClipAtts2.plane2Normal = (-1, 0, 0)
   #
   ClipAtts2.plane3Status = 1
   ClipAtts2.plane3Origin = (0, 0, z)
   ClipAtts2.plane3Normal = (0, 0, -nz)
   ClipAtts2.planeInverse = 1
   ClipAtts2.planeToolControlledClipPlane = ClipAtts2.Plane3  # None, Plane1, Plane2, Plane3
   SetOperatorOptions(ClipAtts2, 1, 1)
   return ClipAtts2


def myClip(z,nz, Lx=0.018):
   AddOperator("Clip", 0)
   ClipAtts = ClipAttributes()
   ClipAtts.quality = ClipAtts.Accurate  # Fast, Accurate
   ClipAtts.funcType = ClipAtts.Plane  # Plane, Sphere
   ClipAtts.plane1Status = 0
   ClipAtts.plane2Status = 0
   #
   ClipAtts.plane2Status = 1
   ClipAtts.plane2Origin = (0., 0, 0)
   ClipAtts.plane2Normal = (+1, 0, 0)
   #
   ClipAtts.plane3Status = 1
   ClipAtts.plane3Origin = (0, 0, z)
   ClipAtts.plane3Normal = (0, 0, -nz)
   ClipAtts.planeInverse = 1
   ClipAtts.planeToolControlledClipPlane = ClipAtts.Plane3  # None, Plane1, Plane2, Plane3
   SetOperatorOptions(ClipAtts, 0, 1)
   return ClipAtts


def zeros(n):
   l = [0. for i in range(n)]
   return l

# Collecting databases:
dbs = glob.glob("GEOM_*sphere/TRANS_*/ijkft_stat_diph_AI.lata")
dbs.extend(glob.glob("GEOM_*sphere/TRANS_*/PAR8/ijkft_stat_diph_AI_par8.lata"))
# dbs = ["GEOM_hemisphere/TRANS_+0.0000/ijkft_stat_diph_AI.lata"]

ndbs = len(dbs)
# Numerical sort : 
# dbs.sort(key=lambda item: int(item.split("/")[0].strip("N")))

jdd="LINKS/sphere_0._SEQ/prepare.data"
Lx=getValue("uniform_domain_size_i", jdd)
Ly=getValue("uniform_domain_size_j", jdd)
Lz=getValue("uniform_domain_size_k", jdd)
nx=int(getValue("nbelem_i", jdd))
ny=int(getValue("nbelem_j", jdd))
nz=int(getValue("nbelem_k", jdd))
dz = Lz/nz
dx = Lx/nx
dy = Ly/ny
vol_cell = dx*dy*dz

# Coordinates of faces and cell centres:
zf=[n*dz for n in range(nz+1)]
zc=[(n+0.5)*dz for n in range(nz)]

first = dbs[0]
OpenDatabase("localhost:"+first, 0)
AddPlot("Pseudocolor", "AIRE_INTERF_ELEM_DOM_EXT", 1, 1)
# AddOperator("Box", 1)
SetActivePlots(0)
# BoxAtts = BoxAttributes()
# BoxAtts.amount = BoxAtts.Some  # Some, All
# BoxAtts.minx = 0
# BoxAtts.maxx = 0.018
# BoxAtts.miny = 0
# BoxAtts.maxy = 0.015
# BoxAtts.minz = zf[0]
# BoxAtts.maxz = zf[1]
# BoxAtts.inverse = 0
# SetOperatorOptions(BoxAtts, 1)
SetQueryFloatFormat("%g")
DrawPlots()

# Le query n'est plus vrai car on est sur DOM_EXT
# Query("Volume"); vol_tranche = GetQueryOutputValue()
vol_tranche = Lx*Ly*dz

DefineScalarExpression("Nx_", "<normals/NORMALE_ELEM_INTERFACES>[0]")
DefineScalarExpression("Ny_", "<normals/NORMALE_ELEM_INTERFACES>[1]")
DefineScalarExpression("Nz_", "<normals/NORMALE_ELEM_INTERFACES>[2]")
DefineScalarExpression("N_", "sqrt(Nx_*Nx_+Ny_*Ny_+Nz_*Nz_)")
DefineScalarExpression("Nx", "Nx_/N_")
DefineScalarExpression("Ny", "Ny_/N_")
DefineScalarExpression("Nz", "Nz_/N_")
DefineScalarExpression("N", "Nx*Nx+Ny*Ny+Nz*Nz")
DefineScalarExpression("ai", "area(INTERFACES)")
DefineScalarExpression("unit", "COURBURE_SOM_INTERFACES*0.+1.")
DefineScalarExpression("aiNx", "ai*Nx")
DefineScalarExpression("aiNy", "ai*Ny")
DefineScalarExpression("aiNz", "ai*Nz")
DefineScalarExpression("kNx", "Nx*COURBURE_SOM_INTERFACES")
DefineScalarExpression("kNy", "Ny*COURBURE_SOM_INTERFACES")
DefineScalarExpression("kNz", "Nz*COURBURE_SOM_INTERFACES")
DefineScalarExpression("kaiNx", "ai*Nx*COURBURE_SOM_INTERFACES")
DefineScalarExpression("kaiNy", "ai*Ny*COURBURE_SOM_INTERFACES")
DefineScalarExpression("kaiNz", "ai*Nz*COURBURE_SOM_INTERFACES")
# DefineScalarExpression("virtuelle", "if(eq(FACETTE_PE_LOCAL_ELEM_INTERFACES-FACETTE_PE_OWNER_ELEM_INTERFACES, 0.*FACETTE_PE_LOCAL_ELEM_INTERFACES), 0., 1.)")
DefineScalarExpression("virtuelle", "if(eq(FACETTE_PE_LOCAL_ELEM_INTERFACES-FACETTE_PE_OWNER_ELEM_INTERFACES, 0.), 0., 1.)")


sth = "# zc AI aiNx aiNy aiNz kaiNx kaiNy kaiNz   ai2 ai3 ?kaiNx2 ?kaiNz2\n"

# ai = ai2 = kai = aiNx = aiNy = aiNz = kaiNx = kaiNy = kaiNz = zeros(nz)
# ai = ai2 = ai3 = aiNx = aiNx2 = kaiNx = zeros(nz)
ai = zeros(nz)
ai2 = zeros(nz)
ai3 = zeros(nz)
aiNx = zeros(nz)
kaiNx = zeros(nz)
kaiNx2 = zeros(nz)
aiNy = zeros(nz)
kaiNy = zeros(nz)
aiNz = zeros(nz)
kaiNz = zeros(nz)
kaiNz2 = zeros(nz)
for i, db in enumerate(dbs):
   OpenDatabase("localhost:"+db, 0)
   ReplaceDatabase("localhost:"+db, 0)
   fold = os.path.dirname(db)
   cas = db.split('/')[0][5:]
   print fold
   print cas
   f = open("%s/%s.dat"%(fold,cas),'w')
   f.write(sth)
   st = "# %s\n"%db
   f.write(st)
   #
   TimeSliderPreviousState()
   if i>0:
      CloseDatabase("localhost:"+dbs[i-1])
      pass
   #
   for j, z in enumerate(zc):
      zf_g = zf[j]   # Face gauche
      zf_d = zf[j+1] # Face droite
      # BoxAtts.minz = zf_g
      # BoxAtts.maxz = zf_d
      # BoxAtts.inverse = 0
      # SetOperatorOptions(BoxAtts, 1)
      RemoveAllOperators()
      #if j==6: import pdb; pdb.set_trace()
      ClipAtts = myClip(zf_g,-1)
      ClipAtts2 = myClip2(zf_d,1)
      DrawPlots()
      SetActivePlots(1)
      # Aire interfaciale : 
      _ = ChangeActivePlotsVar("AIRE_INTERF_ELEM_DOM_EXT")
      Query("Volume"); vol_ext = GetQueryOutputValue()
      _ = Query("Weighted Variable Sum"); ai[j] = GetQueryOutputValue()/vol_tranche # For 3D, it will weight by volume
      print "zf_g zf_d j aij = ",zf_g, zf_d,  j, ai[j], vol_ext, vol_tranche
      # Il ne faut pas faire de variable sum car ca fait comme si la maille etait globalement incluse dans le domaine retenu apres le clip.
      # C'est donc tres sensible selon la position du plan, ca retient ou jete toute la cellule.
      # FAUX = Query("Variable Sum"); ai[j] = GetQueryOutputValue() / (nx*ny)
      #
      # On passe aux grandeurs surfaciques : 
      #
      if (ai[j] >1.):
	 print "ai[j]", ai[j]
	 # _ = ChangeActivePlotsVar("ai")
         # FAUX = Query("Variable Sum"); ai3[j] = GetQueryOutputValue() / vol_tranche
	 if (cas=='hemisphere'):
	    _ = ChangeActivePlotsVar("COURBURE_SOM_INTERFACES"); _ = Query("MinMax", use_actual_data=1); kmin, kmax= GetQueryOutputValue() 
	    pass
	 _ = ChangeActivePlotsVar("unit")
	 #
	 # Seuillage pour supprimer les mailles virtuelles : 
	 #_ = ChangeActivePlotsVar("virtuelle")
	 if 1:
            #import pdb; pdb.set_trace()
	    AddOperator("Threshold", 2)
            ThresholdAtts = ThresholdAttributes()
            ThresholdAtts.outputMeshType = 0
            ThresholdAtts.listedVarNames = ("virtuelle")
            ThresholdAtts.zonePortions = (1)
            ThresholdAtts.lowerBounds = (-1e+37)
            ThresholdAtts.upperBounds = (0.5)
            ThresholdAtts.defaultVarName = "virtuelle"
            ThresholdAtts.defaultVarIsScalar = 1
            SetOperatorOptions(ThresholdAtts, 2, 1)
	    # import pdb; pdb.set_trace()
            DrawPlots()
	    pass
	 #_ = ChangeActivePlotsVar("unit")
	 # Fin du Seuillage pour supprimer les mailles virtuelles.
	 #
         _ = Query("3D surface area"); ai2[j] = GetQueryOutputValue() / vol_tranche
         _ = Query("Weighted Variable Sum"); ai3[j] = GetQueryOutputValue() / vol_tranche
         # Ai * Normale : 
         _ = ChangeActivePlotsVar("Nx")
         _ = Query("Weighted Variable Sum"); aiNx[j] = GetQueryOutputValue() / vol_tranche
         # _ = ChangeActivePlotsVar("aiNx")
         # FAUX _ = Query("Variable Sum"); aiNx2[j] = GetQueryOutputValue() / vol_tranche
         # kappa * Ai * Normale : 
         #_ = ChangeActivePlotsVar("kaiNx")
         # FAUX _ = Query("Variable Sum"); kaiNx[j] = GetQueryOutputValue() / vol_tranche
         _ = ChangeActivePlotsVar("kNx")
         _ = Query("Weighted Variable Sum"); kaiNx[j] = GetQueryOutputValue() / vol_tranche
         SetActivePlots((1, 4))
         # Nouveau bloc pour supprimer la partie courbe:
	 courbure_lower_bound = -600.# -800 -1000
	 courbure_upper_bound = -400.# -200 0
	 if ((cas=='hemisphere') and (kmax>courbure_lower_bound)):
            print "Thresholding on curvature!!!!!!!"
            # AddOperator("Threshold", 2)
            ThAtts = ThresholdAttributes()
            ThAtts.outputMeshType = 0
            ThAtts.listedVarNames = ("COURBURE_SOM_INTERFACES", "virtuelle")
            ThAtts.zonePortions = (1, 1)
            ThAtts.lowerBounds = (courbure_lower_bound, -1e+37)
            ThAtts.upperBounds = (courbure_upper_bound, 0.5)
            ThAtts.defaultVarName = "COURBURE_SOM_INTERFACES"
            ThAtts.defaultVarIsScalar = 1
            SetOperatorOptions(ThAtts, 2, 1)
            DrawPlots()
	    #  GetQueryParameters("MinMax") >>> {'use_actual_data': 0}
	    #  GetQueryParameters("Weighted Variable Sum") >>> "Nothing" meaning no input parameters. I suppose it only uses actual data then.
            _ = Query("Weighted Variable Sum"); kaiNx2[j] = GetQueryOutputValue() / vol_tranche
	    # 
            _ = ChangeActivePlotsVar("kNz")
            _ = Query("Weighted Variable Sum"); kaiNz2[j] = GetQueryOutputValue() / vol_tranche
            # RemoveLastOperator(1, 0)
	    # Reset threshold to former options:
            SetOperatorOptions(ThresholdAtts, 2, 1)
            DrawPlots()
            pass
         #
         _ = ChangeActivePlotsVar("Ny")
         _ = Query("Weighted Variable Sum"); aiNy[j] = GetQueryOutputValue() / vol_tranche
         _ = ChangeActivePlotsVar("kNy")
         _ = Query("Weighted Variable Sum"); kaiNy[j] = GetQueryOutputValue() / vol_tranche
         #
         _ = ChangeActivePlotsVar("Nz")
         _ = Query("Weighted Variable Sum"); aiNz[j] = GetQueryOutputValue() / vol_tranche
         _ = ChangeActivePlotsVar("kNz")
         _ = Query("Weighted Variable Sum"); kaiNz[j] = GetQueryOutputValue() / vol_tranche
      #except:
      else:
         print "Skipping the %dth slice because there is no interface in it."%j
	 pass
      st = "%0.7f\t%10.5g\t%10.5g\t%10.5g\t%10.5g\t%10.5g\t%10.5g\t%10.5g\t%10.5g\t%10.5g\t%10.5g\n"%(z, ai[j], aiNx[j], aiNy[j], aiNz[j], kaiNx[j], kaiNy[j], kaiNz[j], ai2[j], ai3[j], kaiNx2[j])
      f.write(st)
      RemoveAllOperators()
      del ClipAtts; del ClipAtts2
      pass
   f.close()
   pass

exit(0)


