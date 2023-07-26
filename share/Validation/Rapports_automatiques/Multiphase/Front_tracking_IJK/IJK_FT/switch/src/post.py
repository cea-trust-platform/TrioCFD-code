#!python
import matplotlib.pylab as plt
import numpy as np
import glob, os

import commons.DNSTools as dtool

def getSondesCoords(fic):
    x=getValues("x=", fic)
    y=getValues("y=", fic)
    z=getValues("z=", fic)
    return np.array([x,y,z]).T, len(x)

l = glob.glob("post_perio_*.son")
for i, d in enumerate(["X", "Y", "Z"]): 
   for case in ["perio", "wall"]:
      for subc in ["", "_par8"]:
         if ((case=="wall") and(d!="Z")):
            print("skiping direction %s for case wall"%(d))
            continue;
         case=case+subc # ajout du suffixe _par8 si //
         lfic=glob.glob("post_%s_S%s_T[0-9].son"%(case,d))
         lfic.sort()
         for k, fic in enumerate(lfic):
            id_T=(int)(fic[-5])
            coords,_ = dtool.getSondesCoords(fic)
            x=coords[:,i] # or Y or z... 
            T  = np.loadtxt(fic).T
            T = T[1:]
            Ta = np.loadtxt("post_%s_S%s_TANA%d.son"%(case,d,id_T)).T
            Ta = Ta[1:]
            E  = np.loadtxt("post_%s_S%s_E%d.son"%(case,d,id_T)).T
            E = E[1:]
            mat=np.c_[x, T,Ta,E] # line 0 is time...
            out="%s_S%s_T%d"%(case,d,id_T)
            np.savetxt(out+".txt", mat)
            print("file %s.txt created from info in %s"%(out,fic))
            plt.figure(out)
            plt.plot(x,T, label="interp")
            plt.plot(x,Ta, label="ana")
            plt.plot(x,E, label="ecart")
            plt.grid(True)
            plt.legend(loc=0)
            plt.savefig(out+".png")
            pass
         # Bloc for velocity 
         fic="post_%s_S%s_V%s.son"%(case,d,d)
         if os.path.isfile(fic):
            coords,_ = dtool.getSondesCoords(fic)
            x=coords[:,i] # or Y or z... 
            V  = np.loadtxt(fic).T
            V = V[1:]
            Va = np.loadtxt(fic.replace(".son","ANA.son")).T
            Va = Va[1:]
            Ev = np.loadtxt(fic.replace("V%s.son"%d,"E%s.son"%d)).T
            Ev = Ev[1:]
            mat=np.c_[x, V,Va,Ev] # line 0 was time...
            out="%s_S%s_V%s"%(case,d,d)
            np.savetxt(out+".txt", mat)
            print("file %s.txt created from info in %s"%(out,fic))
            plt.figure(out)
            plt.plot(x,V, label="interp")
            plt.plot(x,Va, label="ana")
            plt.plot(x,Ev, label="ecart")
            plt.grid(True)
            plt.legend(loc=0)
            plt.savefig(out+".png")
            plt.close('all')
            pass
         pass
      pass
   pass

#plt.show()
#l = ["perio_par8_SX_T1.txt", \
#     "perio_par8_SY_T2.txt", \
#     "perio_par8_SZ_T3.txt", \
#     "wall_par8_SZ_T0.txt", \
#     "wall_par8_SZ_T1.txt"]
out="compa-seq-par.txt"
f = open(out, "w")
l = glob.glob("*_par8_S*txt")
print("Comparison of cases writen in %s : "%out,l)
for fpar in l: 
   fseq=fpar.replace("par8_","")
   x, seq=np.loadtxt(fseq, usecols=(0,1)).T
   par=np.loadtxt(fpar, usecols=(1)).T
   f.write("%s\t%g\n"%(fseq,abs(seq-par).max()))
   pass

f.close()
