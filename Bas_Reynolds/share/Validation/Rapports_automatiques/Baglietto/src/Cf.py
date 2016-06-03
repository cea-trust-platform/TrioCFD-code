import optparse
import math
import os

parser = optparse.OptionParser()
(options, args) = parser.parse_args()
# ouverture des fichiers
#nomFic = 'marche_VEF_efstab_alpha02_amont_pb_Ustar.face'
nomFic = args[0]+'_pb_Ustar.face'
#nomFic_bag = 'marche_VEF_efstab_alpha02_amont_Baglietto_pb_Ustar.face'
fic = open(nomFic,'r')
#fic_bag = open(nomFic_bag,'r')
chaine="Bas"
basligne=0
i=0
X=[]
Cf=[]
#Cf_bag=[]
# lecture fichier
for i in range(4):
    while 1:
        ligne=fic.readline()
        if chaine in ligne:
            break

#for i in range(4):
#    while 1:
#        ligne=fic_bag.readline()
#        if chaine in ligne:
#            break

for i in range(5):
    ligne=fic.readline()
while 1:
    ligne=fic.readline()
    if ligne=="\n" :
        break
    liste=ligne.split()
    X.append(float(liste[0]))
    Cf.append(float(liste[10])*2)
fic.close()

#for i in range(5):
#    ligne=fic_bag.readline()
#while 1:
#    ligne=fic_bag.readline()
#    if ligne=="\n" :
#        break
#    liste=ligne.split()
#    Cf_bag.append(float(liste[10])*2)
#fic_bag.close()

# ecriture fichier
nomFic = 'Cf.dat'
fichier = open(nomFic, 'w')
for i in range(len(X)):
    fichier.write('%18.8f %18.8f\n' % (X[i],Cf[i]))
fichier.close()


#nomFic = 'Cf_bag.dat'
#fichier = open(nomFic, 'w')
#for i in range(len(X)):
#    fichier.write('%18.8f %18.8f\n' % (X[i],Cf_bag[i]))
#fichier.close()
