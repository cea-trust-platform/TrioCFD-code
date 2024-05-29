# -*- coding: utf-8 -*-
# Suppor team: Christophe Bourcier, Adrien Bruneton, Aymeric Solonet and Guillaume Surat.

from salome.geom import geomBuilder
from salome.smesh import smeshBuilder
from math import radians, sqrt, pi, ceil, cos, sin
import sys, os, pdb
#sys.path.insert(0, os.environ["ModeleCoeur_project_directory"] + "/share/geometrie_maillage")
from manage_groups import make_geometry
from mesh_generation import make_and_export_mesh_in_med_from_geom
import medcoupling as mc

geompy = geomBuilder.New()
smesh = smeshBuilder.New()

p0 = geompy.MakeVertex(0, 0, 0)
x_axis = geompy.MakeVector(p0, geompy.MakeVertex(1, 0, 0))
y_axis = geompy.MakeVector(p0, geompy.MakeVertex(0, 1, 0))
z_axis = geompy.MakeVector(p0, geompy.MakeVertex(0, 0, 1))
axes = [x_axis, y_axis, z_axis]

def build_geom(step_name, x_bottom_tet, x_top_tet, eps, radius):
    # the geometry is imported
    grille = geompy.ImportSTEP(step_name, False, True)
    _, _, ymin, ymax, zmin, zmax = geompy.BoundingBox(grille, True)
    fluid_domain = geompy.MakeBox(x_bottom_tet - eps, ymin, zmin, x_top_tet + eps, ymax, zmax)
    # bearing parameters
    s = 1.65 * radius               # should be enough
    bearing_thickness = 1e-3 + 1e-5   # thickness of the bearing
    pitch = 2 * radius + 0.0028     # distance between rodes
    # create one bearing
    bearing_single = geompy.MakeBox(0.43, 0, -bearing_thickness * 0.5, 0.44, s, bearing_thickness * 0.5)
    bearing_single = geompy.MakeRotation(bearing_single, x_axis, pi * 0.75)
    # create liste of bearing and add them
    bearing = [geompy.MakeTranslation(bearing_single, 0, dy, dz) for dy in [-pitch, 0., pitch] for dz in [-pitch, 0., pitch]]
    geompy.addToStudy(geompy.MakeCompound(bearing), "bearing")
    # create fluid domain
    fluid_domain = geompy.MakeCutList(fluid_domain, [grille] + bearing)
    fluid_domain = geompy.UnionFaces(fluid_domain)

    return fluid_domain, []

def extrude_top_and_bottom(mesh, d, dxh, dx):
    
    mesh3d = mc.ReadUMeshFromFile(mesh)

    ####################
    # The mesh contains only the volumes and vertex. It is necessary to build faces.
    # mf: face mesh ; F_E: face to elements ; F_Ei: indice related to F_E
    # F_E is needed to identified boundaries
    mf, _, _, F_E, F_Ei = mesh3d.buildDescendingConnectivity()
    (xmin, xmax), (_, _), (_, _) = mf.getBoundingBox()

    groups = {"bottom" : [], "top" : []}
    # bf : barycenter of cells
    bf = mf.computeIsoBarycenterOfNodesPerCell()
    for i, b in enumerate(bf):
        # check is the face belongs to a boundary
        if F_Ei[i + 1] == F_Ei[i] + 1:  
            if abs(b[0] - xmin) < 1e-5:
                groups["bottom"].append(i)
            elif abs(b[0] - xmax) < 1e-5:
                groups["top"].append(i)

    meshes = []
    
    
    ####################
    # build a speudo 1D meshes (in 3D) needed for the extrusion. This 1D mesh consider the two zone with different reffinement criteria
    # and build mesh from the 1D meshes created
    for g, x, sgn in [("bottom", xmin, -1.0), ("top", xmax, 1.0)]:
        mesh_bord = mf.buildPartOfMySelf(sorted(groups[g]), False)
        nz1, nz2 = int(d[g][0] / dxh + 1e-7), int(d[g][1] / dx + 1e-7)
        dx1, dx2 = d[g][0] / nz1, d[g][1] / nz2
        print(f"dxh redefined to {dx1}, dx redefined to {dx2}")
        mesh1d = mc.MEDCouplingUMesh("ex", 1)
        mesh1d.allocateCells(nz1 + nz2)
        coo = [x, 0., 0.]
        x_ = x
        i_ = 0
        for dx_, nz in [(dx1, nz1), (dx2, nz2)]:
            for i in range(nz):
                x_ += sgn * dx_
                print(x_)
                mesh1d.insertNextCell(mc.NORM_SEG2, 2, [i_, i_ + 1])
                i_ += 1
                coo.extend([x_, 0., 0.])
        # 1D table of DataArrayDouble needed for medcoupling. 3 if for 3D.
        coo_arr = mc.DataArrayDouble(coo, nz1 + nz2 + 1, 3)
        mesh1d.setCoords(coo_arr)

        mesh = mesh_bord.buildExtrudedMesh(mesh1d, 0)
        mesh.setName("mesh")

        meshes.append(mesh)

    # merge meshes to the 3D original mesh
    mesh = mc.MEDCouplingUMesh.MergeUMeshes([mesh3d] + meshes)
    mesh.mergeNodes(1e-8)
    # merge coincident vertex of the different meshes: mandatory
    mesh.zipCoords()
    mesh.convertAllToPoly()

    return mesh

def create_groups(mesh, r, x1, x2):
    mf, _, _, F_E, F_Ei = mesh.buildDescendingConnectivity()
    (xmin, xmax), (_, _), (_, _) = mf.getBoundingBox()

    groups = {"inlet" : [], "outlet" : [], "wall_hot": [], "wall" : []}
    bf = mf.computeIsoBarycenterOfNodesPerCell()
    for i, b in enumerate(bf):
        if F_Ei[i + 1] == F_Ei[i] + 1:  # bord
            if abs(b[0] - xmin) < 1e-5:
                groups["inlet"].append(i)
            elif abs(b[0] - xmax) < 1e-5:
                groups["outlet"].append(i)
            else:
                # 1e-6 tolerance to the barycentre detection
                if (b[0] < x2 and b[0] > x1 and abs(b[1] * b[1] + b[2] * b[2] - r * r) < 1e-6):
                    groups["wall_hot"].append(i)
                else:
                    groups["wall"].append(i)

    mm = mc.MEDFileUMesh.New()
    # 0 for volume
    mm.setMeshAtLevel(0, mesh)
    # -1 for face
    mm.setMeshAtLevel(-1, mf)
    # add boundary conditions
    for (n, g) in [(k, sorted(v)) for k, v in groups.items()]:
        grp = mc.DataArrayInt.New(g)
        grp.setName(n)
        mm.addGroup(-1, grp)

    return mm

if __name__ == "__main__":

    ####################
    # geometry
    x_bottom_tet, x_top_tet = 0.42, 0.45      # bottom on top axis of the tetrahedrons zone
    x_bottom_dom, x_top_dom = -220e-3, 110e-3 # bottom on top axis of the tetrahedrons zone
    eps = 1e-2                                # distance added upstream and downstream of the 3*3 rod bundle mesh with tetrahedrons, needed for the extrusion.
    r = 5e-3                                  # radius of the rod
    #----------
    dx = 8e-4 # 4e-4                          # mesh refinement of the tetrahedrons zone
    step_name = "../../CAD/070721_full_novane.stp"

    ####################
    # build the geometry for the tetrahedrons zone 
    master, compos = build_geom(step_name, x_bottom_tet, x_top_tet, eps, r)
    d = {"FACE" :
         ["wall_hot", [geompy.MakeVertex(0.421, 0.8 * r, 0)]]
         }
    
    mesh_name="mesh_tetrahedrons.med"
    geometry, groups = make_geometry(master, compos, [], groups_near=d, faces_of_vol_groups=["fluide"])
    make_and_export_mesh_in_med_from_geom(
        geometry,
        groups,
        dx=dx,
        angle_mesh=8,
        # doms = ["fluide", "solide"],
        # elem2d=[("faces_of_fluide", 0.025)]
        # volmesh_ijk=[("ihx2", nr_ihx, nt_ihx), ("periph", nr_periph, nt_periph)]
        mesh_name=mesh_name
    )
    maillage = "mesh/"+mesh_name

    # à ce stade on a un maillage avec des CLs. Mais les CLs vont disparaître à l'extrusion. On pourrait donc simplifier 
    # On pourrait donc simplifier la partie consernant "group" précédement défini.

    ####################
    # extrude the tetrahedrons zone in two directions
    # two zones are defined in both directions
    dx_zone1 = 5*dx
    dx_zone2 = 25* dx
    # extrude the bottom and top face and merge meshes to the tetrahedron mesh.
    mesh = extrude_top_and_bottom(maillage, {"bottom" : [220e-3 - 3e-2 - eps, 0.2], "top" : [110e-3 - eps, 0.2]}, dx_zone1, dx_zone2)
    mm = create_groups(mesh, r, x_top_tet + x_bottom_dom, x_top_tet + x_top_dom)
    mm.write("mesh_novane.med", 2)
