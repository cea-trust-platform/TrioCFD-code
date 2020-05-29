#Script de post-traitement des donnees.

import optparse
import math
import os


def properties():
    # ouverture des fichiers
    nomFic = 'test_VEF.data'
    nomFic_prop='properties.dat'

    fic = open(nomFic,'r')
    fic_prop = open(nomFic_prop,'w')

    chaines = ["# hauteur_marche",
               "BOTTO      Paroi_flux_impose",
               "rho Champ_Uniforme",
               "Cp Champ_Uniforme",
               "# U0",
               "IN         Frontiere_ouverte_temperature_imposee"] # Texte a rechercher

    for ligne in fic:
        for chaine in chaines:
            if chaine in ligne:
                tLigne = ligne.split()
                if chaine=="# hauteur_marche":
                    h=float(tLigne[2])
                if chaine=="BOTTO      Paroi_flux_impose":
                    Qw=float(tLigne[4])
                if chaine=="rho Champ_Uniforme":
                    rho=float(tLigne[3])
                if chaine=="Cp Champ_Uniforme":
                    Cp=float(tLigne[3])
                if chaine=="# U0":
                    U0=float(tLigne[2])
                if chaine=="IN         Frontiere_ouverte_temperature_imposee":
                    T0=float(tLigne[4])


    fic_prop.write('%18.4f %18.4f %18.4f %18.4f %18.4f %18.4f' % (h,Qw,rho,Cp,U0,T0))

    fic.close()
    fic_prop.close()

    return h,Qw,rho,Cp,U0,T0

def ecrire_Stanton(sonde,h,Qw,rho,Cp,U0,T0):

    fichier_read=open(sonde,'r')
    fichier_write=open('Stanton.dat','w')

    ligne = fichier_read.readlines()
    n=len(ligne)

    for i in range (1,n):

        tligne=ligne[i].split()
        fichier_write.write('%18.8f %18.8f\n' % (float(tligne[0])/h,Qw/(rho*Cp*U0*(float(tligne[1])-T0))))

    fichier_read.close()
    fichier_write.close()

def normalise_Stanton(fichier):

    fichier_read=open(fichier,'r')
    fichier_write=open('Stanton_normalise.dat','w')

    ligne=fichier_read.readlines()
    n=len(ligne)

    tligne0=ligne[1].split()
    max=float(tligne0[1])

    for i in range (2,n-1):
        tligne=ligne[i].split()
        if (float(tligne[1])>max):
            max=float(tligne[1])

    for i in range (1,n):
        tligne=ligne[i].split()
        fichier_write.write('%18.8f %18.8f\n' %(float(tligne[0]),float(tligne[1])/max))

    fichier_read.close()
    fichier_write.close()


def normalise_profil_vitesse(fichier,h,U0):

    fichier_read=open(fichier,'r')
    fichier_write=open(fichier+'_normalise.dat','w')

    ligne = fichier_read.readlines()
    n=len(ligne)

    for i in range (1,n):

        tligne=ligne[i].split()
        fichier_write.write('%18.8f %18.8f\n' % (float(tligne[0])/h,float(tligne[1])/U0))

    fichier_read.close()
    fichier_write.close()

def normalise_profil_temperature(fichier,h,T0):

    fichier_read=open(fichier,'r')
    fichier_write=open(fichier+'_normalise.dat','w')

    ligne = fichier_read.readlines()
    n=len(ligne)

    for i in range (1,n):

        tligne=ligne[i].split()
        fichier_write.write('%18.8f %18.8f\n' % (float(tligne[0])/h,float(tligne[1])-T0))

    fichier_read.close()
    fichier_write.close()

def extract_utau(h):
# extrait les valeurs de x et utau sur la face BOTTOM du domaine, au dernier temps. Ecrit fichier avec utau=f(x/h)
# ATTENTION: PROCEDURE ADAPTEE A LA TAILLE DU MAILLAGE UTILISE

# ouverture des fichiers
    nomFicUstar = 'test_VEF_pb_Ustar.face'

    ficUstar = open(nomFicUstar,'r')
    fichier_write=open('Ustar_bottom.dat','w')

    # lecture de ligne -> entetes
    fichier = ficUstar.readlines()

    #tant que la derniere ligne est vide
    while fichier[-1]=="" or fichier[-1]=="\n":
        del fichier [-1]

    #on enleve les 55 dernieres lignes (jusqu'a la face bottom) avant les donnees sur la face qui nous interesse: attention a mettre le bon nb de lignes
    i=0
    while i < 55:

        del fichier [-1]
        i=i+1

    commentaire=""
    i=len(fichier)
    cumul=0
    u=0.


    #on va prendre les valeurs de x et utau au dernier temps et sur la face "bottom"
    while commentaire!="----------------":
        ligne=fichier[i-1]
        if ligne != "\n" and ligne != "":
            tLigne = ligne.split("|")
            commentaire = tLigne[0]
            if commentaire!="----------------":
                #on eclate la ligne en un tableau de valeurs
                tLigne=ligne.split("|")
                #on ecrit x et utau:
                x=float(tLigne[0])
                utau=float(tLigne[4])


                fichier_write.write('%18.8f %18.8f\n' %(x/h,utau))
        i=i-1

    fichier_write.close()
    ficUstar.close()

def ecrire_xr(fichier,h):
# fonction qui calcule la longueur de reattachement

    xr_expe=0.2508
    xr_expe_adim=6.6
    fichier_read=open(fichier,'r')
    fichier_write=open('longueur_reattachement.dat','w')

    ligne = fichier_read.readlines()
    n=len(ligne)

    i=n-1
    tligne2=ligne[i].split()
    tligne1=ligne[i-1].split()
    V2=float(tligne2[1])
    V1=float(tligne1[1])

    while (V1 > 0 and i>2):
        i=i-1
        tligne2=ligne[i].split()
        tligne1=ligne[i-1].split()
        V2=float(tligne2[1])
        V1=float(tligne1[1])

    x2=float(tligne2[0])
    x1=float(tligne1[0])

    if i==2:
        xr=0
        print('attention: pas de point de reattachement trouve, il semble ne pas y avoir de zone de recirculation')
    else:
        xr = x1 - (V1*((x2-x1)/(V2-V1)))

    error=abs(((xr-xr_expe)/xr_expe)*100)
    fichier_write.write('%18.2f %18.2f %18.2f %18.2f %18.2f\n' %(xr,xr/h,error,xr_expe,xr_expe_adim))

    fichier_read.close()
    fichier_write.close()

if __name__ == '__main__':

    parser = optparse.OptionParser()
    (options, args) = parser.parse_args()

    h,Qw,rho,Cp,U0,T0=properties()

    sonde_Twall='test_VEF_SONDE_TWALL.coupe'

    ecrire_Stanton(sonde_Twall,h,Qw,rho,Cp,U0,T0)
    normalise_Stanton('Stanton.dat')

    extract_utau(h)

    #profils de vitesse
    normalise_profil_vitesse('test_VEF_SONDE_V1_LONGI.coupe',h,U0)
    normalise_profil_vitesse('test_VEF_SONDE_V2_LONGI.coupe',h,U0)
    normalise_profil_vitesse('test_VEF_SONDE_V3_LONGI.coupe',h,U0)
    normalise_profil_vitesse('test_VEF_SONDE_V4_LONGI.coupe',h,U0)
    normalise_profil_vitesse('test_VEF_SONDE_V5_LONGI.coupe',h,U0)
    normalise_profil_vitesse('test_VEF_SONDE_V6_LONGI.coupe',h,U0)

    #profils de temperature
    normalise_profil_temperature('test_VEF_SONDE_T1.coupe',h,T0)
    normalise_profil_temperature('test_VEF_SONDE_T2.coupe',h,T0)
    normalise_profil_temperature('test_VEF_SONDE_T3.coupe',h,T0)
    normalise_profil_temperature('test_VEF_SONDE_T4.coupe',h,T0)
    normalise_profil_temperature('test_VEF_SONDE_T5.coupe',h,T0)


    #calcul de la longueur de reattachement
    ecrire_xr('test_VEF_SONDE_BOTTOM_GRAV.coupe',h)
