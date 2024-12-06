from __future__ import division

import os
import pdb as pdb
import GEOM
import SMESH
import MEDLoader as ml
from salome.geom import geomBuilder
from salome.smesh import smeshBuilder
from salome.HYBRIDPlugin import HYBRIDPluginBuilder
import shutil

geompy = geomBuilder.New()
smesh = smeshBuilder.New()

z_axis = geompy.MakeVector(geompy.MakeVertex(0, 0, 0), geompy.MakeVertex(0, 0, 1))


def split_groups_doms(shape, doms):
    # faces touchant chaque domaine
    dface = {}
    for g in geompy.GetGroups(shape):
        if g.GetName() in doms:
            dface[g.GetName()] = set([geompy.GetSubShapeID(shape, subshape) for subshape in geompy.SubShapeAll(g, geompy.ShapeType["FACE"])])

    # est-ce que les groupes de face touchent fluide et/ou solide?
    dgr = dict([(dom, set()) for dom in doms])
    for g in geompy.GetGroups(shape):
        if g.GetMaxShapeType() != GEOM.SOLID:  # groupe de face
            set_of_faces = set([geompy.GetSubShapeID(shape, subshape) for subshape in geompy.SubShapeAll(g, geompy.ShapeType["FACE"])])
            for dom, dfaces in dface.items():
                if bool(dfaces.intersection(set_of_faces)):
                    gn = g.GetName()
                    if gn.startswith("DUPLIC_"):
                        gn = gn[7:]
                    if not gn.startswith("INTERN_"):
                        dgr[dom].add(gn)

    return dgr


def make_and_export_mesh_in_med_from_xao(
    xao_file,
    nb_boundary_layers=0,
    elem2d=[],
    def_elem2d=0,
    elem3d="tets",
    dx=0.3,
    inner_b=True,
    angle_mesh=13,
    scaling=[1, 1, 1],
    trans=[0, 0, 0],
    rename_dup=True,
    mesh_coeur=None,
    volmesh_tets=[],
    volmesh_hybrid=[],
    volmesh_prism=[],
    volmesh_ijk=[],
    doms=["fluide", "solide"],
):
    (imported, geometry, [], groups, []) = geompy.ImportXAO(xao_file)
    make_and_export_mesh_in_med_from_geom(geometry, groups, nb_boundary_layers, elem2d, def_elem2d, elem3d, dx, inner_b, angle_mesh, scaling, trans, rename_dup, mesh_coeur, volmesh_tets, volmesh_hybrid, volmesh_prism, volmesh_ijk, doms)


def make_and_export_mesh_in_med_from_geom(
    geometry,
    groups,
    nb_boundary_layers=0,
    elem2d=[],
    def_elem2d=0,
    elem3d="tets",
    dx=0.3,
    inner_b=True,
    angle_mesh=13,
    scaling=[1, 1, 1],
    trans=[0, 0, 0],
    rename_dup=True,
    mesh_coeur=None,
    volmesh_tets=[],
    volmesh_hybrid=[],
    volmesh_prism=[],
    volmesh_ijk=[],
    doms=["fluide", "solide"],
    mesh_name="primary.med"
):
    primary = mesh_name #"primary.med"

    # Path of the directory to check, remove, and create
    directory = "mesh"

    # Check if the directory exists
    if os.path.exists(directory):
        # Check if it is a directory
        if os.path.isdir(directory):
            # Remove the directory and its contents
            shutil.rmtree(directory)
            print(f"The directory '{directory}' has been removed.")
        else:
            print(f"'{directory}' exists but is not a directory.")
    # Create the directory
    os.makedirs(directory)
    print(f"The directory '{directory}' has been created.")

    
    #os.system("rm -f mesh/{}".format(primary))

    geompy.addToStudy(geometry, "geometrie")
    with open("volumes.csv", "w") as f:
        for g in groups:
            props = geompy.BasicProperties(g)
            print("surface_{:s} : {}".format(g.GetName(), props[1]), file=f)
            if props[2] > 1e-6:
                print("volume_{:s} : {}".format(g.GetName(), props[2]), file=f)
            geompy.addToStudyInFather(geometry, g, g.GetName())
    os.system("cat volumes.csv")
    dgr = split_groups_doms(geometry, doms)
    # write_groups(geometry)

    # 2 - maillage
    mcol = smesh.Mesh(geometry, "mesh")

    # 2.1 - maillage surfacique
    msurf = mcol.Triangle(algo=smeshBuilder.MG_CADSurf)
    msurf.Parameters().SetPreCADProcess3DTopology(True)
    msurf.Parameters().SetGeometricMesh(1)
    msurf.Parameters().SetPhySize(dx)
    msurf.Parameters().SetAngleMesh(angle_mesh)
    msurf.Parameters().SetElementType(def_elem2d)
    msurf.Parameters().AddOption("split_overconstrained_surface_elements", "yes")

    # 2.2 - eventuels raffinements
    for gn, raff in elem2d:
        for g in groups:
            if g.GetName() == gn:
                print("Reglages custom sur groupe {}".format(gn))
                msurf.Parameters().SetSizeMap(g, "def f(u,v): return %f" % (raff))
                msurf.Parameters().SetGradation(1.1)

    # 2.3 - maillage volumique
    if elem3d == "tets":
        mvol = mcol.Tetrahedron(algo=smeshBuilder.MG_Tetra)
        mvol.Parameters().SetOptimizationLevel(4)
        mvol.Parameters().SetGradation(1.001)
        # mvol.Parameters().SetAdvancedOption("--target_quality 0.39")  # --recover_sharp_angles yes --ridge_angle 360')
    elif elem3d == "hybrid":
        mvol = mcol.Tetrahedron(algo=smeshBuilder.HYBRID)
        # mvol.Parameters().SetElementGeneration(HYBRIDPluginBuilder.Generation_Hexa_Dominant)
        mvol.Parameters().SetElementGeneration(HYBRIDPluginBuilder.Generation_Cartesian_Core)
        # mvol.Parameters().SetAdvancedOption('--boundary_layer_max_element_angle 120')
        # mvol.Parameters().SetAdvancedOption('--boundary_layer_height_relative_to_local_surface_size yes')
        # mvol.SetHeightFirstLayer(1)
        mvol.Parameters().SetNbOfBoundaryLayers(0)
        mvol.Parameters().SetCoreSize(dx)
    elif elem3d == "hexa":
        mvol = mcol.Hexahedron(algo=smeshBuilder.MG_Hexa)
        mvol.SetMinMaxHexes(5, 16).SetHexoticSdMode(3)
        mvol.Parameters().SetAdvancedOption("--max_memory 20000 --allow_invalid_elements yes --ridge_angle 30")  # --recover_sharp_angles yes --ridge_angle 360')
        # mvol.Parameters().SetAdvancedOption('--ridge_angle 60')
        # mvol = mcol.Tetrahedron(algo=smeshBuilder.HYBRID)
        # mvol.Parameters().SetElementGeneration(HYBRIDPluginBuilder.Generation_Cartesian_Core)
        # mvol.Parameters().SetAdvancedOption('--boundary_layer_max_element_angle 120')
        # mvol.Parameters().SetAdvancedOption('--boundary_layer_height_relative_to_local_surface_size yes')
        # mvol.SetHeightFirstLayer(1)
        # mvol.Parameters().SetNbOfBoundaryLayers(0)
        # mvol.Parameters().SetCoreSize(dx)

    for g in groups:
        mcol.Group(g)
        if g.GetName() in volmesh_tets:
            mvol2 = mcol.Tetrahedron(algo=smeshBuilder.MG_Tetra, geom=g)
            mvol2.Parameters().SetOptimizationLevel(3)
            mvol2.Parameters().SetGradation(1.001)
        elif g.GetName() in volmesh_hybrid:
            mvol2 = mcol.Tetrahedron(algo=smeshBuilder.HYBRID, geom=g)
            mvol2.Parameters().SetElementGeneration(HYBRIDPluginBuilder.Generation_Cartesian_Core)
            mvol2.Parameters().SetAdvancedOption("--boundary_layer_max_element_angle 120")
            mvol2.Parameters().SetAdvancedOption("--boundary_layer_height_relative_to_local_surface_size yes")
            mvol2.SetHeightFirstLayer(1)
            mvol2.Parameters().SetNbOfBoundaryLayers(0)
            mvol2.Parameters().SetCoreSize(dx)
        elif g.GetName() in [vn[0] for vn in volmesh_prism]:
            _vn, _sn0, _sn1 = volmesh_prism[[vn[0] for vn in volmesh_prism].index(g.GetName())]
            gs0 = groups[[gs.GetName() for gs in groups].index(_sn0)]
            gs1 = groups[[gs.GetName() for gs in groups].index(_sn1)]
            mvol2 = mcol.Prism(geom=g)
            # define projection between the inner and outer surfaces
            msurf2 = mcol.Triangle(algo=smeshBuilder.MG_CADSurf, geom=gs0)
            msurf2.Parameters().SetElementType(1)  # quadrangles
            msurf2.Parameters().SetPhySize(dx)
            msurf2.Parameters().SetAngleMesh(angle_mesh)
            mcol.Projection1D2D(gs1).SourceFace(gs0)  # projection faces[0] -> faces[1]
        elif g.GetName() in [vn[0] for vn in volmesh_ijk]:
            _vn, nr, nt = volmesh_ijk[[vn[0] for vn in volmesh_ijk].index(g.GetName())]
            for wire in geompy.ExtractShapes(g, geompy.ShapeType["EDGE"]):
                kos = geompy.KindOfShape(wire)
                mseg = mcol.Segment(geom=wire)
                if kos[0] == GEOM.GEOM_IKindOfShape.ARC_CIRCLE:
                    mseg.NumberOfSegments(nt)
                else:
                    if nr < 0:
                        mseg.NumberOfSegments(nt)
                    else:
                        if kos[0] != GEOM.GEOM_IKindOfShape.SEGMENT:
                            kos = geompy.KindOfShape(geompy.MakeWire([wire]))
                        # on detecte les segments positionnes sur un diametre
                        is_hori = abs(kos[3] - kos[-1]) < 1e-5
                        if is_hori:
                            if isinstance(nr, list):
                                mseg.FixedPoints1D(nr, [1], [])
                            else:
                                if kos[0] != GEOM.GEOM_IKindOfShape.SEGMENT:
                                    mseg.NumberOfSegments(nt)
                                else:
                                    mseg.NumberOfSegments(nr)
                        else:
                            mseg.LocalLength(dx)

            mcol.Quadrangle(algo=smeshBuilder.QUADRANGLE, geom=g)
            mcol.Hexahedron(algo=smeshBuilder.Hexa, geom=g)

    if mesh_coeur is not None:
        print("prise en compte du maillage MC sur :")
        faces_g = []
        for g in groups:
            if g.GetName().startswith("core_boundaries"):
                print(g.GetName())
                faces_g += geompy.SubShapeAll(g, geompy.ShapeType["FACE"])
        assert faces_g is not None, "Probleme : pas de faces coeur dans la geometrie collecteur"
        # faces_g = geompy.MakeCompound(faces_g)
        gr_faces_g = geompy.CreateGroup(geometry, geompy.ShapeType["FACE"])
        geompy.UnionList(gr_faces_g, faces_g)
        geompy.addToStudyInFather(geometry, gr_faces_g, "core_boundaries")
        gface = []
        for g in mesh_coeur.GetGroups():
            if g.GetName().startswith("faces"):
                gface.append(g)
        assert gface, "Probleme : pas de groupe de faces dans le maillage MC"
        mcol.UseExisting2DElements(geom=gr_faces_g).SourceFaces(gface)

    if not mcol.Compute():
        raise Exception("le maillage n'est pas correctement genere!")

    # 3 - eventuellement : translation et conversion [mm] -> [m]
    mcol.TranslateObject(mcol, trans, 0)
    mcol.Scale(mcol, SMESH.PointStruct(0, 0, 0), scaling, 0)
    mcol.ExportMED("mesh/tmp.med")

    # 4 - traitement Ksing, frontieres internes
    g2dup = []
    for g in mcol.GetGroups():
        if g.GetName().startswith("INTERN_"):
            g.SetName(g.GetName()[7:])
            mcol.ExportMED("mesh/{}".format(primary), overwrite=0, meshPart=g, autoDimension=False)
            print("{} --> {}".format(g.GetName(), primary))
            mcol.RemoveGroup(g)
        elif g.GetName().startswith("DUPLIC_"):
            g.SetName(g.GetName()[7:])
            g2dup.append(g.GetName())
        elif g.GetName().startswith("faces_of_"):
            mcol.RemoveGroup(g)

    # Get Information About Mesh by GetMeshInfo
    print("\nInformation about mesh by GetMeshInfo:")
    info = smesh.GetMeshInfo(mcol)
    keys = list(info.keys())
    keys.sort()
    for i in keys:
        print("  %s   :  %d" % (i, info[i]))

    mcol.ExportMED("mesh/tmp2.med")
    mfum = ml.MEDFileUMesh.New("mesh/tmp2.med", "mesh")

    if inner_b:
        # mailles des groupes de volumes
        vols = [(g.GetName(), set(mfum.getGroupArr(0, g.GetName()).getValues())) for g in mcol.GetGroups(elemType=SMESH.VOLUME) if g.GetName() not in doms]
        for n in g2dup:
            print("{} --> duplication".format(n))
            _, cells, cells_dup = mfum.buildInnerBoundaryAlongM1Group(n)  # elements touchant le groupe n, le groupe n_dup
            faces, faces_dup = mfum.getGroupArr(-1, n), mfum.getGroupArr(-1, n + "_dup")
            if len(faces) != len(faces_dup):
                print("{} --> duplication ratee ({} vs {})".format(n, faces.getNumberOfTuples(), faces_dup.getNumberOfTuples()))
            for l in [l for l in dgr.values() if n in l]:
                l.add(n + "_dup")
            cells = set(cells.getValues())
            cells_dup = set(cells_dup.getValues())
            names = set([volname for (volname, volcells) in vols if cells.intersection(volcells)])
            names_dup = set([volname for (volname, volcells) in vols if cells_dup.intersection(volcells)])
            if names != names_dup:  # on peut changer le nom des deux groupes
                names, names_dup = names.difference(names_dup), names_dup.difference(names)
                for old in [n, n + "_dup"]:
                    mfum.removeGroup(old)
                for old, flist, noms in [
                    (n, faces, names),
                    (n + "_dup", faces_dup, names_dup),
                ]:
                    new = "_".join([n] + sorted(list(noms))[:1])
                    print("{} --> {}".format(old, new))
                    for l in [l for l in dgr.values() if n in l]:
                        l.add(new)
                    flist.setName(new)
                    mfum.addGroup(-1, flist)

    for dom, bords in dgr.items():
        with open("mesh/bords_" + dom + ".txt", "w") as f:
            for n in sorted(list(bords)):
                f.write("{} ".format(n))
            f.write("\n")

    mfum.write("mesh/{}".format(primary), 0)
    md = ml.MEDFileData("mesh/{}".format(primary))
    os.system("rm -f mesh/{}".format(primary))
    md.write33("mesh/{}".format(primary), 0)
