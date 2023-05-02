# -*- coding: utf-8 -*-

from math import *
from numpy import *

Db = 1e-3
delta = 5e-3
Lx = delta*2*pi
Ly = delta*2*pi/3.
Lz = delta*2

nby = 6# Nombre de bulle par direction (Lx est plus long, donc 3*..)
nbz = nby
nbx = 3*nby

# Dans la direction normale aux parois : 
ez = float(Lz  - nbz * Db)/(nbz+1)

# Dans les autres directions : 
ex = Lx/nbx - Db
ey = Ly/nby - Db

# Les origines des premieres bulles : 
ox = ex/2. + Db/2.
oy = ey/2. + Db/2.
oz = ez + Db/2.

# Les espaces inter_bulles : 
mvx = ex + Db
mvy = ey + Db
mvz = ez + Db

vecx = array([ox + i*mvx for i in range(nbx)])
vecy = array([oy + i*mvy for i in range(nby)])
vecz = array([oz + i*mvz for i in range(nbz)])

coordinates = [(x, y, z) for x in vecx for y in vecy for z in vecz]
# Option manquante dans certaines versions : 
#savetxt("deplacement.txt", coordinates, newline=' bulle1mm.msh\n')
savetxt("deplacements.txt", coordinates, fmt='%g %g %g bulle1mm.msh')

