import math
import sys
if len(sys.argv) !=3 :
   raise Exception(f"usage: python {sys.argv[0]} Nx Ny")

per_regu = 0.333   # percentage of regular mesh
Nx_regu = int(sys.argv[1])    # Number of grid = Num of node - 1
fac = 1.05       # factor progressive > 1

Nx_prog = math.log(1. - (1.- per_regu)/per_regu * (1. -fac)*Nx_regu, fac)
print ("[DIR X] number of NODES in progeresive zone :", int(Nx_prog) + 1)

per_regu = 0.333 # percentage of regular mesh
Ny_regu = int(sys.argv[2])    # Number of grid = Num of node - 1
fac = 1.05       # factor progressive > 1

Ny_prog = math.log(1. - (1.- per_regu)/per_regu * (1. -fac)*Ny_regu, fac)
print ("[DIR Y] number of NODES in progeresive zone :", int(Ny_prog) + 1)

