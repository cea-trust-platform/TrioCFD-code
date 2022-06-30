#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@author: GB
"""
#%%
import matplotlib.pyplot as plt
import numpy as np
import sys, os,glob
cwd = os.getcwd()

phiref=1. # reference value to the flux
Ti_ref=0.0015 # reference value to the interfacial temperature
# Hard coded ratios list : 
ratios=np.array([1., 1.1, 1.5, 2., 2., 5., 10., 100.])
nb_ratios=len(ratios)

lis=glob.glob("ascii.lata.INTERFACES_INTERFACE_*_ELEM.2")
lis.sort()
for idx,fic in enumerate(lis):
  name=fic.replace("ascii.lata.INTERFACES_INTERFACE_","").replace("_ELEM.2", "")
  name=name.replace("_","\t")
  mat = np.loadtxt(fic)
  l=mat[0] # The list length
  if ("PHIN" in name):
    values = mat[1:-1]-phiref
  else:
    values = mat[1:-1]-Ti_ref
    pass
  print("%15s\t%8.8g\t%8.8g\t%8.8g\t%8.8g"%(name,ratios[idx%nb_ratios], (np.abs(values)).max(), values.mean(), values.std()))
  pass


