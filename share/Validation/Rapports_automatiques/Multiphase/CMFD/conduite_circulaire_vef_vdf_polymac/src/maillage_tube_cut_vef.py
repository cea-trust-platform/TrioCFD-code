#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.10.0 with dump python functionality
###

import sys
import salome
import os

wd = os.getcwd()
print(wd)

salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()
sys.path.insert(0, wd)

###
### GEOM component
###

import GEOM
from salome.geom import geomBuilder
import math
import SALOMEDS

#############################################################################
############################# Mesh parameters ###############################
#############################################################################

Ny=[11,13] # nombre de noeud sur le rayon du tube (direction radiale = normal à la paroi)
R = 2 # [m] rayon cylindre
L = 100 # [m] hauteur cylindre
cell_size = [R/(N-1) for N in Ny] # [m]

#############################################################################
#############################################################################
#############################################################################

def main(R,L,name,cell_size):
  geompy = geomBuilder.New()

  O = geompy.MakeVertex(0, 0, 0)
  OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
  OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
  OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)
  Cylinder_1 = geompy.MakeCylinderRH(R, L)
  outlet = geompy.CreateGroup(Cylinder_1, geompy.ShapeType["FACE"])
  geompy.UnionIDs(outlet, [10])
  wall = geompy.CreateGroup(Cylinder_1, geompy.ShapeType["FACE"])
  geompy.UnionIDs(wall, [3])
  inlet = geompy.CreateGroup(Cylinder_1, geompy.ShapeType["FACE"])
  geompy.UnionIDs(inlet, [12])
  [outlet, wall, inlet] = geompy.GetExistingSubObjects(Cylinder_1, False)
  geompy.addToStudy( O, 'O' )
  geompy.addToStudy( OX, 'OX' )
  geompy.addToStudy( OY, 'OY' )
  geompy.addToStudy( OZ, 'OZ' )
  geompy.addToStudy( Cylinder_1, 'Cylinder_1' )
  geompy.addToStudyInFather( Cylinder_1, outlet, 'outlet' )
  geompy.addToStudyInFather( Cylinder_1, wall, 'wall' )
  geompy.addToStudyInFather( Cylinder_1, inlet, 'inlet' )

  ###
  ### SMESH component
  ###

  import  SMESH, SALOMEDS
  from salome.smesh import smeshBuilder

  smesh = smeshBuilder.New()
  #smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                  # multiples meshes built in parallel, complex and numerous mesh edition (performance)

  mesh = smesh.Mesh(Cylinder_1,name)
  MG_CADSurf = mesh.Triangle(algo=smeshBuilder.MG_CADSurf)
  MG_CADSurf_Parameters_1 = MG_CADSurf.Parameters()
  MG_CADSurf_Parameters_1.SetPhySize( cell_size )
  MG_CADSurf_Parameters_1.SetChordalError( 5 )
  MG_Tetra = mesh.Tetrahedron(algo=smeshBuilder.MG_Tetra)
  MG_Tetra_Parameters_1 = MG_Tetra.Parameters()
  MG_Tetra_Parameters_1.SetOptimizationLevel( 2 )
  MG_Tetra_Parameters_1.SetMinSize( cell_size/2 )
  MG_Tetra_Parameters_1.SetMaxSize( cell_size*4 )
  MG_Tetra_Parameters_1.SetToMeshHoles( 0 )
  MG_Tetra_Parameters_1.SetToMakeGroupsOfDomains( 1 )
  MG_Tetra_Parameters_1.SetMaximumMemory( -1 )
  MG_Tetra_Parameters_1.SetInitialMemory( -1 )
  MG_Tetra_Parameters_1.SetKeepFiles( 0 )
  MG_Tetra_Parameters_1.SetWorkingDirectory( '/tmp/' )
  MG_Tetra_Parameters_1.SetVerboseLevel( 10 )
  MG_Tetra_Parameters_1.SetStandardOutputLog( 0 )
  MG_Tetra_Parameters_1.SetRemoveLogOnSuccess( 0 )
  outlet_1 = mesh.GroupOnGeom(outlet,'outlet',SMESH.FACE)
  wall_1 = mesh.GroupOnGeom(wall,'wall',SMESH.FACE)
  inlet_1 = mesh.GroupOnGeom(inlet,'inlet',SMESH.FACE)
  isDone = mesh.Compute()
  Domain_1 = mesh.GetGroups()[ 3 ]
  [ outlet_1, wall_1, inlet_1, Domain_1 ] = mesh.GetGroups()
  smesh.SetName(mesh, name)
  med_file = wd+"/"+name+".med"
  try:
    mesh.ExportMED( med_file, 0, 41, 1, mesh, 1, [], '',-1, 1 )
    pass
  except:
    print('ExportPartToMED() failed. Invalid file name?')


  ## Set names of Mesh objects
  smesh.SetName(MG_CADSurf.GetAlgorithm(), 'MG-CADSurf')
  smesh.SetName(MG_Tetra.GetAlgorithm(), 'MG-Tetra')
  smesh.SetName(MG_CADSurf_Parameters_1, 'MG-CADSurf Parameters_1')
  smesh.SetName(MG_Tetra_Parameters_1, 'MG-Tetra Parameters_1')
  smesh.SetName(outlet_1, 'outlet')
  smesh.SetName(wall_1, 'wall')
  smesh.SetName(inlet_1, 'inlet')
  smesh.SetName(mesh.GetMesh(), name)
  smesh.SetName(Domain_1, 'Domain_1')


  if salome.sg.hasDesktop():
    salome.sg.updateObjBrowser()

for i in range(len(cell_size)):
  name = f"vef_mesh_n{Ny[i]-1}"
  main(R,L,name,cell_size[i])
