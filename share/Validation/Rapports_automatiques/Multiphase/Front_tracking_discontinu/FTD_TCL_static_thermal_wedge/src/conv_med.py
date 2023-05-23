
import medcoupling as mc
import numpy as np
from math import *

fName = "lata/post_0000.med"

euler = mc.ReadUMeshFromFile(fName, "dom")
interf = mc.ReadUMeshFromFile(fName, "INTERFACES")

eps = 1.0e-12

c = interf.getNodalConnectivity()
cI =interf.getNodalConnectivityIndex()

# Orienting the front to be going forward:
#print(c.getValues())
a = c[1::3]
b = c[2::3]
c[1::3] = b
c[2::3] = a
#print(c.getValues())

interf.setConnectivity(c, cI)
# f = interf.getMeasureField(True)
# print(f.getArray().getValues())

# Creation of interf2 which is the cut version of interf, 
#    cut by elements from dom (euler)
euler2, interf2, mp1, mp2 = euler.Intersect2DMeshWith1DLine(euler,interf, eps)

# To clean interf2 of unnecessary points:
interf2.zipCoords()
cdg_interf2 = interf2.computeCellCenterOfMass()

mp = 'MPOINT_THERMIQUE_ELEM_dom'
ai = 'INTERFACIAL_AREA_ELEM_dom'
if mp in mc.GetAllFieldNames(fName):
   lst_dt = mc.GetAllFieldIterations(fName,mp) # take last time-step
   fC = mc.ReadFieldCell(fName,"dom",0,mp,lst_dt[-1][0],lst_dt[-1][1])
   fai = mc.ReadFieldCell(fName,"dom",0,ai,lst_dt[-1][0],lst_dt[-1][1])
   mpa = fC.getArray()
   aia = fai.getArray()
   cells,cellsIndex = euler.getCellsContainingPoints(cdg_interf2, len(cdg_interf2),eps )
   #print(cells)
   mpa_part = mpa[cells]
   aia_part = aia[cells]
   mat = np.c_[cdg_interf2.toNumPyArray(),mpa_part.toNumPyArray(),aia_part.toNumPyArray()]
   np.savetxt("mpoint.txt", mat)
   pass

coord2 = interf2.getCoords()
xminmax, yminmax = interf2.getBoundingBox()
xmin, _ = xminmax
ymin, _ = yminmax

xysm = np.zeros((0,5))
xysm_cells = np.zeros((0,5))
s=0
c1pre, _ = interf2.getNodeIdsOfCell(0) # Initialize
xpre, ypre = interf2.getCoordinatesOfNode(c1pre)
#print("Do we start at lower left corner of the bounding box?")
#print(xmin, xpre)
#print(ymin, ypre)

for idx in range(interf2.getNumberOfCells()): 
   c0, c1 = interf2.getNodeIdsOfCell(idx)
   m = mpa_part[idx] # The mpoint of the element
   ai = aia_part[idx] # Corresponding interface area of interf2 in the element
   if c1pre == c0:
      # Fine, its contiguous
      c1pre = c1 # Increment it for next passage
   else:
      print("Element: ", idx)
      print("Previous element Edge #1", c1pre)
      print("New element Edge #0", c0)
      raise Exception("Interfacial mesh is not contiguous") 
   x0 , y0 = interf2.getCoordinatesOfNode(c0)
   x1 , y1 = interf2.getCoordinatesOfNode(c1)
   xi , yi  = cdg_interf2[idx].getValues() # of the bary in the element idx
   s0 = s
   edge_length = sqrt((x0-x1)**2+(y0-y1)**2)
   s += edge_length
   s1 = s
   # First segment point:
   tmp1 = np.array((x0,y0,s0,m, ai))
   # second segment point:
   tmp2 = np.array((x1,y1,s1,m, ai))
   # At the center:
   tmpi = np.array((xi,yi,s0+0.5*edge_length,m, ai))
   # Dumps : 
   xysm = np.vstack([xysm, tmp1, tmp2])
   xysm_cells = np.vstack([xysm_cells,tmpi])
   pass

np.savetxt("xysm.dat", xysm, header="x y s m ai")
np.savetxt("xysm_cells.dat", xysm_cells, header="xi yi si m ai")
# Reloading : 
xysm = np.loadtxt("xysm.dat")
xysm_cells = np.loadtxt("xysm_cells.dat")

vx,vy,vs,vm,va = xysm.T
vxi,vyi,vsi,vmi,vai = xysm_cells.T

def readFileInfo():
    data = {}
    with open("info.txt") as f:
        for line in f.readlines():
            key, value = line.rstrip("\n").split("=")
            data[key] = float(value)
    return data

data = readFileInfo()
lda = data["lda"]
dT = data["DT"]
Lvap = data["Lvap"]
theta_app = data["theta"]*pi/180
print("CASE: dT=", dT, "  -- theta_app=",data["theta"],)
# Analytical solution from Vadim:
smin = 2e-7
vs_ana = np.linspace(smin,vs.max(),201)
qi_ana = lda*dT/(vs_ana*theta_app)
mp_ana = qi_ana/Lvap
np.savetxt("sm_ana.dat", np.vstack([vs_ana,mp_ana]).T, header="s m_ana")

import matplotlib.pyplot as plt
plt.figure()
plt.title("mp")
plt.xlabel(r"$s$ [m]")
plt.ylabel(r"$\dot{m}$ [kg/m²/s]")
plt.plot(vs, vm,"k-", label="Simu")
plt.plot(vs_ana, mp_ana,"r-", label="Analytical")
#plt.plot(vsi, vmi,"ko", label="at centre")
plt.legend(loc=0)
plt.savefig("mp.png")

plt.figure()
plt.title("ai")
plt.xlabel(r"$s$ [m]")
plt.ylabel(r"$a_i$ [m²]")
plt.plot(vs, va,"k-", label="over element")
plt.plot(vsi, vai,"ko", label="at centre")
plt.legend(loc=0)
plt.savefig("ai.png")

plt.figure()
plt.title("mp*ai")
plt.xlabel(r"$s$ [m]")
plt.ylabel(r"$\dot{m} a_i$ [kg/s]")
plt.plot(vs, vm*va,"k-", label="over element")
plt.plot(vsi, vmi*vai,"ko", label="at centre")
plt.legend(loc=0)
plt.savefig("mpai.png")
#plt.show()

#print("THE END")

