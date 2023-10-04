import medcoupling as mc
from math import radians, sin, cos, ceil
import numpy as np
def build_mesh(alpha_d, folder, med_name):

    # rayon maillage
    r = 0.0096
    rs = 0.0106
    # hauteur maillage
    h = 5.
    # quelque chose
    n =  20
    ns = 10
    nz = 200

    theta_mesh_size = 2. # en degres
    nt = 1 #ceil(alpha_d / theta_mesh_size) # equivalent axisymetrique (une maille 3D d'epaisseur)
    alpha = radians(alpha_d)

    # maillage de la section de 1 degre
    meshes = []

    x = [i / nt for i in range(nt + 1)]
    y = [i / (n+ns) for i in range(n + ns + 1)]
    arrx = mc.DataArrayDouble(x)
    arry = mc.DataArrayDouble(y)

    # maillage cartesien
    mccm = mc.MEDCouplingCMesh.New("section")
    mccm.setCoords(arrx, arry)
    # passage en maillage non structure
    mesh_e = mccm.buildUnstructured()
    # conversion des elements en polyhedres pour PolyMAC
    mesh_e.convertAllToPoly()
    # transformation des coordonnees des noeuds
    coo = mesh_e.getCoords()
    for i in range(coo.getNumberOfTuples()):
#        r_ = (coo[i, 1] * n/(n+1) + 1/(n+1) ) * r
        if coo[i, 1] <= n / (n+ns) + 1.e-6 : r_ = coo[i, 1] * (n+ns) / n * r
        else : 
            print(coo[i, 1])
            print(coo[i, 1]-n/(n+ns))
            r_ = (coo[i, 1]-n/(n+ns)) * (n+ns) / ns * ( rs - r ) + r
        x = r_ * cos(alpha * coo[i, 0])
        y = r_ * sin(alpha * coo[i, 0])
        coo[i, 0] = x
        coo[i, 1] = y
    meshes.append(mesh_e)

    mesh = mc.MEDCouplingUMesh.MergeMeshes(meshes)
    mesh.mergeNodes(1e-6)
    mesh.zipCoords()
    mesh.convertDegeneratedCells()

    # maillage d'extrusion
    z = [0]
    z.extend([i * h / nz for i in range(1, nz + 1)])

    mesh1d = mc.MEDCouplingUMesh("ex", 1)
    mesh1d.allocateCells(len(z) - 1)
    coo = [0., 0., z[0]]
    for i, z_ in enumerate(z[1:]):
        mesh1d.insertNextCell(mc.NORM_SEG2, 2, [i, i + 1])
        coo.extend([0., 0., z_])
    coo_arr = mc.DataArrayDouble(coo, len(z), 3)
    mesh1d.setCoords(coo_arr)

    # extrusion de mesh par mesh1d
    mesh.changeSpaceDimension(3, z[0])
    mesh = mesh.buildExtrudedMesh(mesh1d, 0)
    mesh.setName("mesh")
    mesh.convertAllToPoly()
    mesh.orientCorrectlyPolyhedrons()

    # groupes de cellules 3D
    bc = mesh.computeCellCenterOfMass()
    g_solide = []
    g_liquide = []
    for i, b in enumerate(bc):
        if (b[0]*b[0]+b[1]*b[1]) > r * r: g_solide.append(i)
        else: g_liquide.append(i)
        
    # maillage des faces
    mf, desc, descIndx, F_E, F_Ei = mesh.buildDescendingConnectivity()
    mf.setName("mesh")
    # groupes de cellules 2D
    bf = mf.computeCellCenterOfMass()
    z_b, z_t = 0, h
    g_t, g_b, g_s, g_o = [], [], [], [] #top, bottom, symmetry of liquid, outside of solid
    for i, b in enumerate(bf):
        if b[1]*b[1]+b[0]*b[0] < r*r : #if liquid    	
            if   abs(b[2] - z_t) < 1e-5: g_t.append(i)
            elif abs(b[2] - z_b) < 1e-5: g_b.append(i)
            elif ( b[1]<1.e-8 ) : g_s.append(i)
            elif ( b[0]>1.e-8 ) and (abs(alpha - np.arctan(b[1]/b[0]))) < 1e-5: g_s.append(i) # use of pi/2-angle because we know b[1]>0
        elif   abs(b[2] - z_t) < 1e-5: g_o.append(i)
        elif abs(b[2] - z_b) < 1e-5: g_o.append(i)
        elif (F_Ei[i + 1] == F_Ei[i] + 1) : g_o.append(i) #if solid
        
    # structure pour stocker les maillages avec leurs groupes
    mm = mc.MEDFileUMesh.New()
    mm.setMeshAtLevel(0, mesh)
    for (na, g) in [("fluide", g_liquide), ("solide", g_solide)]:
        grp = mc.DataArrayInt.New(g)
        grp.setName(na)
        mm.addGroup(0, grp)
        
    mm.setMeshAtLevel(-1, mf)
    # ajout des groupes
    for (na, g) in [("top", g_t), ("bottom", g_b), ("symetrie", g_s), ("outside", g_o)]:
        grp = mc.DataArrayInt.New(g)
        grp.setName(na)
        mm.addGroup(-1, grp)

    mm.write(f"{folder}/{med_name}", 2)

if __name__ == "__main__":
    build_mesh(2.0, ".", "Debora_20_mailles_solide.med")
