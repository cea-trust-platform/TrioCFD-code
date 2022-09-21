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

# MISE EN DONNEES ######################################################
MEDFileName = "./result_sR_0000.med"
FieldName = "VELOCITY_ELEM_DOM_dual"

MEDFileField1, MEDFileMesh1, TOF1 = methodes_2020.MEDFileToMEDField_1TS(MEDFileName, FieldName)
order = MEDFileField1.getOrder()

MEDFileField, MEDFileMesh, TOF = methodes_2020.MEDFileToMEDField_MultiTS(MEDFileName, FieldName)
it = MEDFileField.getIterations(); tps = it[3]
MEDchamp = MEDFileField.getFieldOnMeshAtLevel(TOF,tps[0],tps[1],0,MEDFileMesh,0)
MEDmesh  = MEDFileMesh.getMeshAtLevel(0)
temps = [tps[0], tps[1], order]
# ######################################################################

# TEST DES FONCTIONS SUR LE CREEPING FLOW ##############################
if False :
    NPTV,MEDTV = outils_physique.TimeAvg(MEDFileName,FieldName)
    methodes_2020.MEDFieldToMEDFile(MEDTV,"./result_sR_AVG.med")
    # Puis on regarde sur Paraview : ca a bien une tete de moyenne temporelle
    # Faudrai demander a GB, ou c'est déjà que j'ai la moyenne avec des chiffres quoi .

    
if False :
    # ATTENTION : Applatit ne boucle pas tout seul sur le temps
    #             au cas ou on fasse une moyenne temporelle d'abord
    #             puis une moyenne par plan. La on a pris le 3 instant
    NPCPX,MEDCPX = outils_physique.Applatit(MEDchamp,temps,"X")
    NPCPY,MEDCPY = outils_physique.Applatit(MEDchamp,temps,"Y")
    NPCPZ,MEDCPZ = outils_physique.Applatit(MEDchamp,temps,"Z")
    methodes_2020.MEDFieldToMEDFile(MEDCPX,"./result_sR_PLX.med")
    methodes_2020.MEDFieldToMEDFile(MEDCPY,"./result_sR_PLY.med")
    methodes_2020.MEDFieldToMEDFile(MEDCPZ,"./result_sR_PLZ.med")

if False :
    # ATTENTION : Rij, calcule les tensions de Reynolds à tous les
    #             instants de la simulation.
    NPRij, MEDRij = outils_physique.Rij(MEDFileName, FieldName)
    methodes_2020.MEDFieldToMEDFile(MEDRij,"./result_sR_Rij.med")

if True :
    #  ATTENTION : jacob prend directemet en entrée un MEDField : ce
    #              n'est pas forcément pratique si on souhaite boucler
    #              sur le temps comme on le voit ici.
    #              -->(Ameliorer jacob ?...)
    
    ####################################################################
    # Preparation pour le MEDFile
    NewMEDFileName = MEDFileName.split('/')[-1]     # get only file name if path given tambien
    NewMEDFileName = "jac_"+NewMEDFileName
    ml.WriteUMesh(NewMEDFileName,MEDmesh,True)  #SI ICI ON N'A PAS UN UMESH, REMPLACER PAR WriteMesh
    
    for i,t in it:
        MEDc = MEDFileField.getFieldOnMeshAtLevel(TOF,i,t,0,MEDFileMesh,0)
        NPjac, MEDjac = outils_physique.jacob(MEDc)
        ml.WriteFieldUsingAlreadyWrittenMesh(NewMEDFileName,MEDjac     )
        
    

