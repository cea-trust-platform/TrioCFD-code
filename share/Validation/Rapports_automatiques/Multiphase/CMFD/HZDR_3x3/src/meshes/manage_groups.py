# -*- coding: utf-8 -*-
from __future__ import division

from math import cos, pi, radians, sin, sqrt, atan2
import GEOM
from salome.geom import geomBuilder
from salome.smesh import smeshBuilder
from time import time

geompy = geomBuilder.New()
smesh = smeshBuilder.New()

p0 = geompy.MakeVertex(0, 0, 0)
z_axis = geompy.MakeVector(p0, geompy.MakeVertex(0, 0, 1))
x_axis = geompy.MakeVector(p0, geompy.MakeVertex(1, 0, 0))


class Point:
    def __init__(self, x, y, group=None, poly=False, fillet=None):
        self.x = x
        self.z = y
        self.group = group
        self.poly = poly
        self.f = fillet


def get_xy(M, N, pitch):
    x = pitch * ((M - 30) + (N - 30) / 2.0)
    y = pitch * sqrt(3) / 2.0 * (N - 30)
    return x, y


def import_mc_mesh(boundary="mesh/boundary.med"):
    mc_boundary, m_coarse = None, None
    (mc_meshes, status) = smesh.CreateMeshesFromMED(boundary)
    for m in mc_meshes:
        if m.GetName() == "fine":
            mc_boundary = m
        if m.GetName() == "coarse":
            m_coarse = m

    assert mc_boundary is not None, "on n'a pas trouve le maillage de frontiere mc"
    assert m_coarse is not None, "on n'a pas trouve le maillage grossier de la frontiere mc"

    # import du .stl dans geom, et transformation en volume
    m_coarse.ExportSTL("mesh/boundary_coarse.stl", 1)
    core_boundaries = geompy.UnionFaces(geompy.ImportSTL("mesh/boundary_coarse.stl"))
    core_volume = geompy.UnionFaces(geompy.MakeSolid([core_boundaries]))
    core_boundaries = geompy.MakeCompound(geompy.ExtractShapes(core_volume, geompy.ShapeType["FACE"]))

    core_boundaries.SetName("core_boundaries")
    core_volume.SetName("cut_coeur")
    return core_boundaries, core_volume, mc_boundary


def create_group_from(name, mother_shape, list_elem, type="FACE"):
    new = geompy.CreateGroup(mother_shape, geompy.ShapeType[type])
    geompy.UnionList(new, list_elem)
    new.SetName(name)
    geompy.addToStudyInFather(mother_shape, new, name)
    return new


def create_group_near(name, mother_shape, pos, list_grp, type, add=True):
    shapes = [geompy.GetShapesNearPoint(mother_shape, P, geompy.ShapeType[type]) for P in pos]
    new = create_group_from(name, mother_shape, shapes, type)
    if add:
        list_grp.append(new)
    return new, list_grp


def create_defaut(name, mother_shape):
    all_faces_id, ids_groups = set(geompy.GetFreeFacesIDs(mother_shape)), set()
    for grp in geompy.GetGroups(mother_shape):
        ids_groups |= set(grp.GetSubShapeIndices())

    if all_faces_id - ids_groups:
        faces_g = [geompy.GetSubShape(mother_shape, [x]) for x in all_faces_id - ids_groups]
        defaut = create_group_from(name, mother_shape, faces_g, "FACE")
        return defaut


def create_group_sharing(shape, list_of_names, name):
    # liste de set d'IDs de face pour chaque groupe
    list_ids = []
    for s in geompy.GetGroups(shape):
        if s.GetName() in list_of_names:  # list_ids.append(set(s.GetSubShapeIndices()))
            list_ids.append(set([geompy.GetSubShapeID(shape, subshape) for subshape in geompy.SubShapeAll(s, geompy.ShapeType["FACE"])]))

    # on cherche le volume ayant un nombre > 0 de face en commun avec chaque groupe
    for subshape in geompy.SubShapeAll(shape, geompy.ShapeType["SOLID"]):  # set_of_faces = set(subshape.GetSubShapeIndices())
        set_of_faces = set([geompy.GetSubShapeID(shape, face) for face in geompy.SubShapeAll(subshape, geompy.ShapeType["FACE"])])
        if all([len(setg.intersection(set_of_faces)) for setg in list_ids]):
            create_group_from(name, shape, [subshape], "SOLID")


def create_group_sharing2(shape, list_of_names, name):
    # liste de set d'IDs de face pour chaque groupe
    candidates = []
    list_ids, list_not = [], []
    for s in geompy.GetGroups(shape):
        set_of_faces = set([geompy.GetSubShapeID(shape, subshape) for subshape in geompy.SubShapeAll(s, geompy.ShapeType["FACE"])])
        if s.GetName() in list_of_names:  # list_ids.append(set(s.GetSubShapeIndices()))
            list_ids.append(set_of_faces)
        else:
            list_not.append(set_of_faces)

    # on cherche le volume ayant un nombre > 0 de face en commun avec chaque groupe
    for subshape in geompy.SubShapeAll(shape, geompy.ShapeType["SOLID"]):  # set_of_faces = set(subshape.GetSubShapeIndices())
        set_of_faces = set([geompy.GetSubShapeID(shape, face) for face in geompy.SubShapeAll(subshape, geompy.ShapeType["FACE"])])
        if all([len(setg.intersection(set_of_faces)) for setg in list_ids]):
            candidates.append([subshape, sum([len(setg.intersection(set_of_faces)) for setg in list_not])])
    assert candidates, "Aucun volume en contact avec {} trouve".format(list_of_names)
    return create_group_from(name, shape, [sorted(candidates, key=lambda l: l[1])[0][0]], "SOLID")


def create_group_sharing3(shape, list_of_names, name, doms):
    # 1/ verification que les noms dans list_of_names sont bien des groupes de la geometrie
    gn = [gg for gg in geompy.GetGroups(shape) if (gg.GetName().split("__")[0] not in doms)]
    for n in list_of_names:
        if n not in [g.GetName() for g in gn]:
            raise Exception("vgroups : le groupe {} n'est pas dans la geometrie!!".format(n))

    # 2/ stockage liste de set d'IDs de face pour chaque groupe
    list_ids, candidates = [], []
    for s in geompy.GetGroups(shape):
        set_of_faces = set([geompy.GetSubShapeID(shape, subshape) for subshape in geompy.SubShapeAll(s, geompy.ShapeType["FACE"])])
        if s.GetName() in list_of_names:  # list_ids.append(set(s.GetSubShapeIndices()))
            list_ids.append(set_of_faces)

    # 3/ on cherche les groupes de volume qui touchent les groupes de face
    for subshape in geompy.SubShapeAll(shape, geompy.ShapeType["SOLID"]):  # set_of_faces = set(subshape.GetSubShapeIndices())
        set_of_faces = set([geompy.GetSubShapeID(shape, face) for face in geompy.SubShapeAll(subshape, geompy.ShapeType["FACE"])])
        if all([len(setg.intersection(set_of_faces)) for setg in list_ids]):
            candidates.append(subshape)

    # 4/ si on en a trouve un seul : facile
    if len(candidates) >= 1:
        return create_group_from(name, shape, candidates, "SOLID")
    # si on en a trouve un aucun : bizarre...
    elif len(candidates) == 0:
        raise Exception("aucun candidat pour creer le groupe de volume {} a partir de {} :".format(name, list_of_names))
    # si on en a trouve plusieurs, on aide l'utilisateur a choisir!!
    elif len(candidates) > 1:
        print("il y a {} candidats pour creer le groupe {} a partir de {}".format(len(candidates), name, list_of_names))
        dicg = {}
        for g in gn:
            dicg[g.GetName()] = set([geompy.GetSubShapeID(shape, g) for g in geompy.SubShapeAll(g, geompy.ShapeType["FACE"])])
        for i, s in enumerate(candidates):
            set_of_faces = set([geompy.GetSubShapeID(shape, face) for face in geompy.SubShapeAll(s, geompy.ShapeType["FACE"])])
            for k, v in dicg.items():
                if len(v.intersection(set_of_faces)):
                    print("le volume {} / {} a en contact : {}".format(i + 1, len(candidates), k))
        raise


def create_bounding_groups(shape, excl):
    bb = geompy.BoundingBox(shape, True)
    eps = 1e-6
    prefix = ("left", "right", "front", "back", "down", "top")
    vol_grps = [g for g in geompy.GetGroups(shape) if g.GetMaxShapeType() == GEOM.SOLID and g.GetName()]
    vol_grps_ids = {}
    for g in vol_grps:
        vol_grps_ids[g.GetName()] = set([geompy.GetSubShapeID(shape, subshape) for subshape in geompy.SubShapeAll(g, geompy.ShapeType["FACE"])])

    for face in [geompy.GetSubShape(shape, [x]) for x in geompy.GetFreeFacesIDs(shape)]:
        for i, x in enumerate(bb):
            try:
                au_bord = max([abs(x_f - x) for x_f in geompy.BoundingBox(face, True)[2 * (i // 2) : 2 * (i // 2) + 2]]) < eps
            except:
                au_bord = max([abs(x_f - x) for x_f in geompy.BoundingBox(face)[2 * (i // 2) : 2 * (i // 2) + 2]]) < eps

            if au_bord:
                found = False
                face_id = geompy.GetSubShapeID(shape, face)
                # TODO: check pas dans un  groupe
                for name, faces_vol in vol_grps_ids.items():
                    if not found and face_id in faces_vol and name + "_" + prefix[i] not in excl:
                        print("Bounding group : {}_{}".format(name, prefix[i]))
                        found = True
                        create_group_from("{}_{}".format(name, prefix[i]), shape, [face], "FACE")


def redefine_surface_group(shape, name, list_v):
    f_v = set()
    for s in geompy.GetGroups(shape):
        if s.GetName() == name:  # set_ids_surf = s.GetSubShapeIndices()
            set_ids_surf = set([geompy.GetSubShapeID(shape, subshape) for subshape in geompy.SubShapeAll(s, geompy.ShapeType["FACE"])])
            for g in geompy.GetGroups(shape):
                if g.GetMaxShapeType() == GEOM.SOLID and g.GetName() in list_v:
                    faces_v = set([geompy.GetSubShapeID(shape, subshape) for subshape in geompy.SubShapeAll(g, geompy.ShapeType["FACE"])])
                    if f_v:
                        f_v.intersection_update(faces_v)
                    else:
                        f_v.update(faces_v)
            for id in list(set_ids_surf.difference(f_v)):
                geompy.RemoveObject(s, id)


def create_group_between_volumes(shape, name, list_v):
    f_v = []
    for g in geompy.GetGroups(shape):
        if g.GetMaxShapeType() == GEOM.SOLID and g.GetName() in list_v:
            f_v.append(set([geompy.GetSubShapeID(shape, subshape) for subshape in geompy.SubShapeAll(g, geompy.ShapeType["FACE"])]))
    faces = set.intersection(*f_v)
    create_group_from(name, shape, list(faces), "FACE")


def find_first_group_containing(face, vol_grps_ids, other_name):
    for name, ids in vol_grps_ids.items():
        if name != other_name and face in ids:
            return name
    return None


def create_group_between_all_volumes(shape, groups_of_interest):
    vol_grps = [g for g in geompy.GetGroups(shape) if g.GetMaxShapeType() == GEOM.SOLID and g.GetName()]
    vol_grps_ids = {}

    faces_done = set()
    faces_groups = {}
    for g in vol_grps:
        vol_grps_ids[g.GetName()] = set([geompy.GetSubShapeID(shape, subshape) for subshape in geompy.SubShapeAll(g, geompy.ShapeType["FACE"])])

    for name, ids in vol_grps_ids.items():
        for face in ids:
            if face not in faces_done:
                other_grp = find_first_group_containing(face, vol_grps_ids, name)
                if other_grp is not None:
                    faces_done.add(face)
                    grp_name = f"face_{name}_{other_grp}"
                    if grp_name in faces_groups:
                        faces_groups[grp_name].add(face)
                    else:
                        faces_groups[grp_name] = {face}

    for name, faces_ids in faces_groups.items():
        if groups_of_interest is not None and name in groups_of_interest:
            faces = [geompy.GetSubShape(shape, [f]) for f in faces_ids]
            create_group_from(name, shape, faces, "FACE")


def create_faces_of_vol_groups(shape, vol_groups, groups):
    vol_grps = [g for g in geompy.GetGroups(shape) if g.GetMaxShapeType() == GEOM.SOLID and g.GetName() in vol_groups]

    for g in vol_grps:
        groups.append(create_group_from(f"faces_of_{g.GetName()}", shape, geompy.SubShapeAll(g, geompy.ShapeType["FACE"]), "FACE"))


def revol(profil_pt, name, r=0, theta=0, vol=False, cut=False, split=False):
    # 1/ creation du profil
    poly, straight, pts_fillet = [], [], []
    for i, P in enumerate(profil_pt):
        if P.poly:
            poly.append(positionner(geompy.MakeVertex(P.x, 0, P.z), r, 0, radians(theta)))
        else:
            straight.append(positionner(geompy.MakeVertex(P.x, 0, P.z), r, 0, radians(theta)))
        if P.f:
            pts_fillet.append([positionner(geompy.MakeVertex(0, P.x, P.z), r, 0, theta), P.f])

    # 2/ creation de la revolution
    profil = []
    if poly:
        profil.append(geompy.MakeInterpol(poly))
    profil.extend([geompy.MakeEdge(straight[i], straight[i + 1]) for i in range(len(straight) - 1)])
    profil = geompy.MakeWire(profil)
    if vol:
        profil = geompy.MakeFace(profil, 1)
    shape = geompy.MakeRevolution(profil, geompy.MakeRotation(geompy.MakeTranslation(z_axis, r, 0, 0), z_axis, radians(theta)), 2.0 * pi)

    # 3/ creation des fillets
    for [point, f] in pts_fillet:
        e_id = geompy.GetSubShapeID(shape, geompy.GetEdgeNearPoint(shape, point))
        shape = geompy.MakeFillet(shape, f, geompy.ShapeType["EDGE"], [e_id])

    if split:
        face_1 = geompy.MakeRotation(geompy.MakeTranslation(geompy.MakeFaceHW(100, 100, 3), r, 0, 0), z_axis, radians(theta))
        face_2 = geompy.MakeRotation(geompy.MakeTranslation(geompy.MakeFaceHW(100, 100, 2), r, 0, 0), z_axis, radians(theta))
        shape = geompy.MakePartition([shape], [face_1, face_2])
        if split == 2:
            face_1 = geompy.MakeRotation(geompy.MakeTranslation(geompy.MakeRotation(geompy.MakeFaceHW(100, 100, 3), z_axis, radians(45.0)), r, 0, 0), z_axis, radians(theta))
            face_2 = geompy.MakeRotation(geompy.MakeTranslation(geompy.MakeRotation(geompy.MakeFaceHW(100, 100, 2), z_axis, radians(45.0)), r, 0, 0), z_axis, radians(theta))
            cyl_faces = geompy.MakeCompound([s for s in geompy.ExtractShapes(shape, geompy.ShapeType["FACE"]) if geompy.KindOfShape(s)[0] == GEOM.GEOM_IKindOfShape.CYLINDER2D])
            cyl_faces = geompy.MakePartition([cyl_faces], [face_1, face_2])
            shape = geompy.MakePartition([shape], [cyl_faces])

    shape.SetName("cut_{}".format(name) if cut else name)
    geompy.addToStudy(shape, shape.GetName())
    return shape

def revol2(profil, vol = False, name = None):
    pt0 = profil.pop(0) #premier point
    ptd = pt0           #point de debut de la courbe courante
    edges=[]
    while len(profil):
        curve = profil.pop(0) #type de courbe
        ptf = profil.pop(0) if len(profil) else pt0 #point de fin : le prochain dans la liste si il existe, pt0 sinon
        coin = abs(ptf[0] - ptd[0]) < 1e-8 and abs(ptf[1] - ptd[1]) < 1e-8
        edge = None
        if coin : pass #rien a faire
        elif curve == 'S':
            edges += [geompy.MakeEdge(geompy.MakeVertex(ptd[0], 0, ptd[1]), geompy.MakeVertex(ptf[0], 0, ptf[1]))]
        elif curve[0] == 'T': #courbes avec tangentes de debut / fin
            _, td, tf, pts = curve
            vpts = []
            for pt in pts:
                vpts.append(geompy.MakeVertex(pt[0], 0, pt[1]))
            edges += [geompy.MakeInterpolWithTangents([geompy.MakeVertex(ptd[0], 0, ptd[1])] + vpts + [geompy.MakeVertex(ptf[0], 0, ptf[1])],
                                                   geompy.MakeVectorDXDYDZ(td[0], 0, td[1]), geompy.MakeVectorDXDYDZ(tf[0], 0, tf[1]))]
        elif curve[0] == 'E': #bout d'ellipse avec centre
            edges += [geompy.MakeArcOfEllipse(geompy.MakeVertex(curve[1][0], 0, curve[1][1]),
                                             geompy.MakeVertex(ptd[0], 0, ptd[1]), geompy.MakeVertex(ptf[0], 0, ptf[1]))]
        elif curve[0] == 'I': #bout d'ellipse d'iso-surface en coordonnees cylindriques
            r0, z0 = curve[1][0], curve[1][1]
            #coordonnees polaires par rapport a (r0, z0)
            rthi = [[sqrt((pt[0] - r0)**2 + (pt[1] - z0)**2), atan2(pt[1] - z0, pt[0] - r0)] for pt in [ptd, ptf]]
            if len(curve) >= 3: #sens de parcours impose
                rthi[1][1] += 2 * pi * ((curve[2] == '+' and rthi[1][1] < rthi[0][1]) - (curve[2] == '-' and rthi[1][1] > rthi[0][1]))
            else: #on parcourt l'arc de cercle le plus court
                rthi[1][1] -= 2 * pi * round((rthi[1][1] - rthi[0][1]) / (2 * pi))
            #limites interne / externe / bas / haut
            lim = { "i" : 0, "e" : 1e8, "b" : -1e8, "h" : 1e8 }
            for s, d, m in [("i", 0, 0), ("e", 0, 1), ("b", 1, 0), ("h", 1, 1)]:
                if len(curve) == 4 and s in curve[3]:
                    lim[s] = min(lim[s], max(ptd[d], ptf[d])) if m else max(lim[s], min(ptd[d], ptf[d]))
            #invariant preservant la surface
            ci = [r0 * r + r**2 * cos(th) / 2 for (r, th) in rthi]
            #liste de points de la courbe
            nth = int(50 * abs(rthi[1][1] - rthi[0][1]) / pi) #nombre
            lth = [rthi[0][1] + i * (rthi[1][1] - rthi[0][1]) / nth for i in range(nth + 1)] #angles
            lc = [ci[0] + i * (ci[1] - ci[0]) / nth for i, th in zip(range(nth + 1), lth)] #invariants
            lr = [(sqrt(max(r0**2 + 2 * c * cos(th), 0.)) - r0) / cos(th) if abs(cos(th)) > 1e-8 else c / r0 for (c, th) in zip(lc, lth)] #rayon local
            lp = [(r0 + r * cos(th), z0 + r * sin(th)) for (r, th) in zip(lr, lth)] #points
            isc = [p[0] < lim["i"] or p[0] > lim["e"] or p[1] < lim["b"] or p[1] > lim["h"] for p in lp] #le point est-il clampe?
            lp = [(min(max(p[0], lim["i"]), lim["e"]), min(max(p[1], lim["b"]), lim["h"])) for p in lp]
            llp = [geompy.MakeVertex(lp[0][0], 0, lp[0][1])]
            for i in range(1, len(lp)):
                llp += [geompy.MakeVertex(lp[i][0], 0, lp[i][1])]
                if i > 1 and i + 2 < len(lp) and isc[i] and isc[i - 1] != isc[i + 1]:
                    edges += [geompy.MakeInterpol(llp)]
                    llp = [llp[-1]]
            edges += [geompy.MakeInterpol(llp)] #polyline si on touche r = 0
        else:
            assert 0, "revol2: cannot handle " + str(curve)
        ptd = ptf #point de debut pour la prochaine courbe
    shape = geompy.MakeWire(edges)
    if vol:
        shape = geompy.MakeFace(shape, 1)
    shape = geompy.MakeRevolution(shape, z_axis, 2.0 * pi)
    if name:
        shape.SetName(name)
    return shape


def positionner(shape, r, z, theta):
    return geompy.MakeRotation(geompy.MakeTranslation(shape, r, 0, z), z_axis, theta)


def traversee(shape, r, rr, theta):
    name = shape.GetName()
    if isinstance(r, float) or isinstance(r, int):
        shape_cut = positionner(geompy.MakeCylinderRH(r, 1e4), 0, -5e3, 0)
    else:
        shape_cut = r
    newshape = geompy.MakeCut(shape, positionner(shape_cut, rr, 0, radians(theta)))
    newshape.SetName(name)
    return newshape


def make_and_export_geometry(
    master,
    compos,
    vgroups,
    sgroups=[],
    groups_near=None,
    sec=None,
    z_niveau_libre=1e6,
    rotation_angle=0,
    core_boundaries=None,
    folder=".",
    excl=[],
    doms=["fluide", "solide"],
    create_groups_between_all_vol=False,
    groups_of_interest=None,
    faces_of_vol_groups=[],
):
    geometry, groups = make_geometry(master, compos, vgroups, sgroups, groups_near, sec, z_niveau_libre, rotation_angle, core_boundaries, excl, doms, create_groups_between_all_vol, groups_of_interest, faces_of_vol_groups)

    geompy.ExportXAO(geometry, groups, [], "Francois Pecquery", "{}/{}.xao".format(folder, master.GetName()))
    for g in groups:
        if g.GetMaxShapeType() == GEOM.SOLID:
            geompy.ExportSTL(g, f"{folder}/{g.GetName()}.stl", GEOM.LU_METER)
    geompy.ExportSTEP(geometry, "{}/{}.step".format(folder, master.GetName()), GEOM.LU_METER)


def make_geometry(
    master,
    compos,
    vgroups,
    sgroups=[],
    groups_near=None,
    sec=None,
    z_niveau_libre=1e6,
    rotation_angle=0,
    core_boundaries=None,
    excl=[],
    doms=["fluide", "solide"],
    create_groups_between_all_vol=False,
    groups_of_interest=None,
    faces_of_vol_groups=[],
):
    cut, par = [], []
    dim_max = max([abs(x) for x in geompy.BoundingBox(master, True)])
    geometry = geompy.MakeCut(master, geompy.MakeTranslation(geompy.MakeCylinderRH(2 * dim_max, 2 * dim_max), 0, 0, z_niveau_libre))
    if sec:
        pts = [geompy.MakeVertex(2 * dim_max * cos(radians(s)), 2 * dim_max * sin(radians(s)), -1.5 * dim_max) for s in sec]
        tokeep = geompy.MakeWire(
            [
                geompy.MakeEdge(pts[0], geompy.MakeVertex(0, 0, -1.5 * dim_max)),
                geompy.MakeEdge(geompy.MakeVertex(0, 0, -1.5 * dim_max), pts[1]),
                geompy.MakeArc(
                    pts[1],
                    geompy.MakeVertex(2 * dim_max * cos(radians((sec[0] + sec[1]) / 2)), 2 * dim_max * sin(radians((sec[0] + sec[1]) / 2)), -1.5 * dim_max),
                    pts[0],
                ),
            ]
        )
        tokeep = geompy.MakePrismDXDYDZ(geompy.MakeFace(tokeep, 1), 0, 0, 3 * dim_max)
        geompy.addToStudy(tokeep, "tokeep")
        geometry = geompy.MakeCommon(geometry, tokeep)
    for c in compos:
        if c.GetName().startswith("cut_"):
            cut.append(c)
        else:
            par.append(c)

    if cut:
        geometry = geompy.MakeCutList(geometry, cut)
    if par:
        geometry = geompy.MakePartition([geometry], par)

    # groupes
    t0 = time()
    glob = dict([(dom, []) for dom in doms])
    geompy.addToStudy(geometry, master.GetName())
    for c in par:
        if c.GetName():
            if c.GetMaxShapeType() != GEOM.SOLID:  # l'objet a partitionner est 2D
                print("creating face group {}".format(c.GetName()))
                try:
                    create_group_from(c.GetName(), geometry, [geompy.GetInPlace(geometry, c, 1)])
                except:
                    try:
                        create_group_from(c.GetName(), geometry, [geompy.GetInPlaceByHistory(geometry, c)])
                    except:
                        raise Exception("The group {} can't be created...".format(c.GetName()))

            else:  # l'objet a partitionner est un volume
                assert c.GetName().split("__")[0] in doms, "un groupe de volume doit commencer par".join([dom + "__ ou " for dom in doms])
                print("creating volume group {}".format(c.GetName()))
                globn = c.GetName().split("__")[0]
                grpn = c.GetName().split("__")[1]
                # groupe de faces bordant le solide
                if grpn.endswith("_abg"):  # add boundary group
                    grpn = grpn[:-4]
                    create_group_from(
                        "bord_" + grpn,
                        geometry,
                        [geompy.GetInPlace(geometry, geompy.MakeCompound(geompy.ExtractShapes(c, geompy.ShapeType["FACE"])), 1)],
                    )
                # group de volume
                g = create_group_from(grpn, geometry, [geompy.GetInPlace(geometry, c, 1)], type="SOLID")
                glob[globn].append(g)

    t1 = time()
    print("Groupes : {}s".format(t1 - t0))
    # groupes de volumes
    for sg, vn in vgroups:
        assert vn.split("__")[0] in doms, "{}: un groupe de volume doit commencer par ".join([dom + "__ ou " for dom in doms]).format(vn)
        glob[vn.split("__")[0]].append(create_group_sharing3(geometry, sg, vn.split("__")[1], doms))

    t2 = time()
    print("Groupes volume : {}s".format(t2 - t1))

    # groupes au bord de la bounding_box
    if excl != "*":
        create_bounding_groups(geometry, excl)

    t3 = time()
    print("Groupes bord : {}s".format(t3 - t2))

    # rustine en attendant mieux
    if groups_near:
        for key, gr in groups_near.items():
            name = gr[0]
            pts = gr[1]
            if key == "SOLID":
                assert name.split("__")[0] in doms, "un groupe de volume doit commencer par ".join([dom + "__ ou " for dom in doms])
                globn = name.split("__")[0]
                grpn = name.split("__")[1]
                g, tmp_list = create_group_near(grpn, geometry, pts, [], key)
                glob[globn].append(g)
            else:
                create_group_near(name, geometry, pts, [], key)
    t4 = time()
    print("Groupes groups_near : {}s".format(t4 - t3))

    # groupes de faces entre volumes
    for n, l in sgroups:
        create_group_between_volumes(geometry, n, l)

    if create_groups_between_all_vol:
        create_group_between_all_volumes(geometry, groups_of_interest)

    t5 = time()
    print("Groupes create_group_between_volumes : {}s".format(t5 - t4))

    # regroupement des groupes de meme nom
    groups = []
    for group in geompy.GetGroups(geometry):
        if group.GetName():
            if group.GetName() in [gr.GetName() for gr in groups]:
                idx = [gr.GetName() for gr in groups].index(group.GetName())
                groups[idx] = geompy.UnionGroups(group, groups[idx])
                groups[idx].SetName(group.GetName())
            else:
                groups.append(group)
    for gn, g in glob.items():
        if g:
            gfs = geompy.UnionListOfGroups(g)
            gfs.SetName(gn)
            groups.append(gfs)

    t6 = time()
    print("Groupes meme nom : {}s".format(t6 - t5))

    if core_boundaries is not None:
        for g in groups:
            if g.GetName() in doms:
                print("{}_{}".format(core_boundaries.GetName(), g.GetName()), end=" -> ")
                common_part = geompy.GetInPlace(geometry, geompy.MakeCommon(core_boundaries, g), 1)
                if common_part:
                    print("ajout")
                    grp = create_group_from("{}_{}".format(core_boundaries.GetName(), g.GetName()), geometry, [common_part])
                    groups.append(grp)
                else:
                    print("rien")

    t7 = time()
    print("Groupes core_boundaries : {}s".format(t7 - t6))

    # faces restantes
    groups.append(create_defaut("defaut", geometry))

    create_faces_of_vol_groups(geometry, faces_of_vol_groups, groups)

    t8 = time()
    print("Groupe defaut : {}s".format(t8 - t7))

    geompy.Rotate(geometry, z_axis, radians(rotation_angle))
    groups = [g for g in groups if g is not None]
    return geometry, groups
