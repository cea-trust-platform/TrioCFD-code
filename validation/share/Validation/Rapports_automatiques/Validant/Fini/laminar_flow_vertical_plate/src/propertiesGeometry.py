#Script de recuperation des propretes physiques et de la geometrie

import optparse
import math
import os

def properties(nomFic):
    # ouverture des fichiers
    fic = open(nomFic,'r')

    chaines = ["Longueurs",
               "mu Champ_Uniforme",
               "Temperature Champ_Uniforme",
               "Gauche  paroi_echange_externe_impose T_ext Champ_Front_Uniforme",
               "Gauche  paroi_temperature_imposee Champ_Front_Uniforme",
               "Nombre_de_Noeuds",
               "vitesse Champ_Uniforme",
               "Cp  ",
               "Prandtl ",
               "rho Champ_Uniforme    ",
               "vitesse Champ_Uniforme "] # Texte a rechercher
    Pr=1
    rho=1
    GrL=1000000000

    for ligne in fic:
        for chaine in chaines:
            if chaine in ligne:
                tLigne = ligne.split()
                if chaine=="Longueurs":
                    W=float(tLigne[1])
                    L=float(tLigne[2])
                if chaine=="mu Champ_Uniforme":
                    mu=float(tLigne[3])
                if chaine=="Temperature Champ_Uniforme":
                    Temp=float(tLigne[-2])
                if chaine=="Gauche  paroi_echange_externe_impose T_ext Champ_Front_Uniforme":
                    TempLim=float(tLigne[5])
                if chaine=="Gauche  paroi_temperature_imposee Champ_Front_Uniforme":
                    TempLim=float(tLigne[4])
                if chaine=="Nombre_de_Noeuds":
                    nbMailles=float(tLigne[-1])
                if chaine=="Cp  ":
                    Cp=float(tLigne[-1])
                if chaine=="Prandtl ":
                    Pr=float(tLigne[-1])
                if chaine=="rho Champ_Uniforme    ":
                    rho=float(tLigne[-1])
                if chaine=="vitesse Champ_Uniforme ":
                    U=float(tLigne[-1])
    fic.close()
    return mu,GrL,W,L,Temp,nbMailles,Cp,Pr,rho,U,TempLim


def ecritureFichier(mu,GrL,W,L,Temp,nbMailles,Cp,Pr,rho,U,TempLim):
    #ecriture du fichier
    nomFic = 'propertiesGeometry.dat'
    fichier = open(nomFic, 'w')
    fichier.write('%18.8f %s %18.8f %18.8f %18.8f %18.8f %18.8f %18.8f %18.8f %18.8f %18.8f\n' % (mu,GrL,W,L,Temp,nbMailles,Cp,Pr,rho,U,TempLim))
    fichier.close()

if __name__ == '__main__':

    #recuperation du fichier data
    #derniere ligne du ls
    ficLS = os.popen('ls *.data')
    lignes = ficLS.readlines()
    derLigne = lignes[-1]
    #suppression du \n en fin de nom
    nomFic = derLigne[:len(derLigne)-1]
    mu,GrL,W,L,Temp,nbMailles,Cp,Pr,rho,U,TempLim = properties(nomFic)

    #ecriture du fichier
    ecritureFichier(mu,GrL,W,L,Temp,nbMailles,Cp,Pr,rho,U,TempLim)
