# -*- coding: utf-8 -*-
#import rlcompleter
from Tools import *
import readline
import unittest
import pickle
import os


from matplotlib.pyplot import clf
#from scipy.sparse.linalg.isolve.minres import Ainfo

readline.parse_and_bind('tab:complete')

clf
global compteur
compteur=0
os.system('rm -f *.png')
fstat="R_exeAnna/Stats.dt_ev"
fstat2="R_exebrMarcellin/Stats.dt_ev"

mu = 1.39e-3 
rho = 986.51


# Using set() 
def Diff(li1, li2): 
    return (list(set(li1) - set(li2))) 

############################################################## Loading Field #######################################
dvar = Field.getEntries(fstat)
dvar2 = Field.getEntries(fstat2)

for fname in dvar:
   if fname in dvar2:
      f1=Field.LoadFromFile(fstat,["coordonnee_K"] , [fname],fname, 'z', r'%s'%fname)
      f2=Field.LoadFromFile(fstat2,["coordonnee_K"] , [fname],fname+' 2', 'z', r'%s exeN'%fname)
      err=(f1-f2).abs()._npa.max()
      den=(f1).abs()._npa.max()
      err_rel=err/den
      if (err>1e-15):
         print " ************** %s **************"%fname
         print f1.name, "absolute/relative error max are : ", err, err_rel
         tracer([f1, f2], 'diff')
         import subprocess; subprocess.call(["display diff_out.png"], shell=True)
         #input("Press Enter to continue...")
         pass
   else:
      print "%s is missing in dvar2"%fname
   pass


print "list of names only in dvar:",  Diff(dvar, dvar2)
print "list of names only in dvar2:",  Diff(dvar2, dvar)

only_in_dvar_2 =Diff(dvar2, dvar)
aa=[x for x in only_in_dvar_2 if 'T' not in x]
print "but without T reduces to:", aa
