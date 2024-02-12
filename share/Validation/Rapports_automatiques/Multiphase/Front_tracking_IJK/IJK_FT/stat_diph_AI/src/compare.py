#!/usr/bin/python3
# -*- coding: utf-8 -*-

import numpy as np
import os, glob, sys
import commons.DNSTools as dns

cwd = dns.pwd()
Lz=dns.getValue("uniform_domain_size_k", "model.data")
nz=int(dns.getValue("nbelem_k", "model.data"))

# Collecting databases:
fictxt = "diphasique_moyenne_spatiale_0.000001.txt"
svar = "coordonnee_K AI aiNx aiNy aiNz kaiNx kaiNy kaiNz"
lvar=svar.split()
msg = "#   cas       translate  run variable      err.max         err.std      err_rel.max(%)     err_rel.std(%)\n"
for cas in ["sphere", "hemisphere"]:
   zc, AI, aiNx, aiNy, aiNz, kaiNx, kaiNy, kaiNz = np.loadtxt("%s.ana"%cas,unpack=True)
   for translate in ["+0.0000", "-0.0002", "-0.0041"]:
      for subfold in [".", "PAR8"]:
         fold="GEOM_%s/TRANS_%s/%s"%(cas,translate,subfold)
         f = os.path.join(fold,fictxt)
         if os.path.isfile(f):
            dvar =dns.buildDicoColonnes(f)
            fp = os.path.join(fold,"PAR8", fictxt)
            m=np.loadtxt(f)
            mat = np.zeros((len(lvar),nz)) # for output results
            for i,key in enumerate(lvar):
               icol=dvar[key]
               st="%sx_%s = m[:,%d]" %(cas[0],key,icol)
               print(("\tRunning "+st))
               exec(st)
               st="mat[%d]=%sx_%s"%(i,cas[0],key)
               print(("\tRunning "+st))
               exec(st)
               if icol!=0:
                  exec("err = np.abs((%sx_%s - %s))" %(cas[0],key,key))
                  exec("adim=np.abs(%s).max()" %(key))
                  msg+="%10s\t%s\t%4s\t%5s\t%10.4f\t%10.4f"%(cas,translate,subfold.replace(".","SEQ"),key,err.max(),err.std())
                  if (adim!=0):
                     err_rel = 100.*err/adim
                     msg +="\t%11.1f\t%11.1f\n"%(err_rel.max(),err_rel.std())
                  else:
                     msg +="            undef           undef\n"
                     pass
                  exec("err_rel = 100.*err/np.abs((1e-8+%s))" %(key))
                  #msg+="%10s\t%s\t%4s\t%5s\t%10.4f\t%10.4f\t%11.1f\t%11.1f\n"%(cas,translate,subfold.replace(".","SEQ"),key,err.max(),err.std(),min(err_rel.max(),123456),min(err_rel.std(),123456))
                  pass
               try:
                  np.savetxt(os.path.join(fold,"%s.num"%cas), mat.T, header=svar.replace("coordonnee_K", "zc"))
               except:
                  np.savetxt(os.path.join(fold,"%s.num"%cas), mat.T)
                  pass
               pass
            pass
         pass
      pass
   pass


print(("*"*120))
#msg.replace("123456.0", "     NaN")
print(msg)
print(("*"*120))
f=open('resu.txt', 'w')
f.write(msg)
f.close()
