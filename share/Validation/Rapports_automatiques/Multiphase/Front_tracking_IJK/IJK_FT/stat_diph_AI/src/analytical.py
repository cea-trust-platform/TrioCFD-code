#!/usr/bin/python3
# -*- coding: utf-8 -*-

import numpy as np
import os, glob, sys
import commons.DNSTools as dns

pi = np.pi
R = 0.004
kappa=-2/R

cz = 0.005 # Coord z du Centre de la bulle

svar = "zc AI aiNx aiNy aiNz kaiNx kaiNy kaiNz"
lvar=svar.split()

Lx=dns.getValue("uniform_domain_size_i", "model.data")
Ly=dns.getValue("uniform_domain_size_j", "model.data")
Lz=dns.getValue("uniform_domain_size_k", "model.data")
nz=int(dns.getValue("nbelem_k", "model.data"))
dz = Lz/nz
# Coordinates of faces and cell centres:
zf=np.linspace(0.,Lz, nz+1, endpoint=True)
zc=(zf[1:]+zf[:-1])/2

# The theta angle at faces:
# thetaf=np.arccos(zf/R)
xx=np.maximum(-np.ones(nz+1),np.minimum(np.ones(nz+1),(zf-cz)/R))
thetaf=np.arccos(xx)

inside=abs((zc-cz)/R)<=np.ones(nz)
inside.astype(float)

#########################################
# SPHERE                                #
#########################################
print("Computing analytical values for the sphere...")
C=2*pi/(Lx*Ly)
AI=np.ones(nz)*2*pi*R/(Lx*Ly)*inside
aiNx=np.zeros(nz)
aiNy=np.zeros(nz)
aiNz=2*pi/(Lx*Ly)*(zc-cz)*inside # zc-cz: recentre sur la sphere...
kaiNx=np.zeros(nz)
kaiNy=np.zeros(nz)
kaiNz=aiNz*kappa

mat = np.zeros((len(lvar),nz))
for i,var in enumerate(lvar):
   st="mat[%d]=%s"%(i,var)
   print("\tRunning "+st)
   exec(st)
   pass

try:
   np.savetxt("sphere.ana", mat.T, header=svar)
except:
   np.savetxt("sphere.ana", mat.T)
   pass

#########################################
# HEMISPHERE                            #
#########################################
print("Computing analytical values for the hemisphere...")
C=2*pi/(Lx*Ly)
# L'aire du segment circulaire (cf wiki) sur la face:
Af = 2*thetaf-np.sin(2*thetaf)
# aip : L'aire du troncon de plan.
# ais : L'aire du troncon de la partie spherique. (la moitie du cas precedent)
aip=1/(Lx*Ly*dz)*R*R/2*(Af[:-1]-Af[1:])*inside
ais=pi*R/(Lx*Ly)*inside
AI=aip+ais

# 
aiNx_s = aip # La projection de l'aire sur le plan de coupe de la sphere... 

aiNx= np.zeros(nz) # Bizarre mais surement a cause de la projection par Nx.
aiNy=np.zeros(nz) # par symetrie
aiNz=pi/(Lx*Ly)*(zc-cz)*inside # (la moitie du cas precedent pour la parie spherique et 0 pour la partie plane car Nz=0 sur le plan)

kaiNx= kappa*aiNx_s
kaiNy=np.zeros(nz) # par symetrie
kaiNz=aiNz*kappa # La moitie du cas precedent... 

lvar=svar.split()
mat = np.zeros((len(lvar),nz))
for i,var in enumerate(lvar):
   st="mat[%d]=%s"%(i,var)
   print("\tRunning "+st)
   exec(st)
   pass

try:
   np.savetxt("hemisphere.ana", mat.T, header=svar)
except:
   np.savetxt("hemisphere.ana", mat.T)
   pass

