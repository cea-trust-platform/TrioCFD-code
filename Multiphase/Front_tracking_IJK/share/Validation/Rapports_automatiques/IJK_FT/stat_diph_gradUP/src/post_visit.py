#!visit -nowin -cli -s 

v=sys.argv[-1]
if (v not in ["P", "U", "V", "W"]): raise Exception("Option should be given in P, U, V or W")
import glob

def myQuery(var):
   ret = ChangeActivePlotsVar(var)
   ret2 = Query("Max")
   maxi = GetQueryOutputValue()
   ret3 = Query("Weighted Variable Sum") # For 3D, it will weight by volume
   norm_l2 = GetQueryOutputValue()
#   if (not (ret & ret2 & ret3)): raise Exception("Failed to change field for %s maximum!"%var)
   return maxi, norm_l2

if (v == "P"):
   # Cas tests en pression : 
   #dbs = glob.glob("N*/P111_U00/ijkft_stat_diph_gradUP[_par8].lata")
   dbs = glob.glob("N*/P111_U11/ijkft_stat_diph_gradUP.lata")
   dbs.extend(glob.glob("N*/P111_U11/ijkft_stat_diph_gradUP_par8.lata"))
else:   
   # Cas tests en vitesse : 
   #dbs = glob.glob("N*/P000_U11/ijkft_stat_diph_gradUP[_par8].lata")
   dbs = glob.glob("N*/P111_U11/ijkft_stat_diph_gradUP.lata")
   dbs.extend(glob.glob("N*/P111_U11/ijkft_stat_diph_gradUP_par8.lata"))
   pass

# Numerical sort : 
dbs.sort(key=lambda item: int(item.split("/")[0].strip("N")))

first = dbs[0]
OpenDatabase("localhost:"+first, 0)
if (v == "P"):
   # Cas tests en pression : not everything is cell-based.
   # Face-based vector : 
   DefineVectorExpression("err_d%sd_"%v, "abs(ANA_d%sd_FACES_DOM_dual-d%sd_FACES_DOM_dual)"%(v,v))
   DefineScalarExpression("err_d%sdx"%v, "err_d%sd_[0]"%v)
   DefineScalarExpression("err_d%sdy"%v, "err_d%sd_[1]"%v)
   DefineScalarExpression("err_d%sdz"%v, "err_d%sd_[2]"%v)
   errs = ["err_d%sdx"%v, "err_d%sdy"%v, "err_d%sdz"%v]
   # Cell-based vectors and tensors : 
   l = ["dd%sdd_X"%v, "dd%sdd_Y"%v, "dd%sdd_Z"%v, "dd%sdxdy"%v, "dd%sdxdz"%v, "dd%sdydz"%v] 
   suite = []
else: 
   # Cas tests en vitesse : all gradients are cell-based (even the first gradients!
   # Cell-based vectors and tensors : 
   l = ["d%sd_X"%v, "d%sd_Y"%v, "d%sd_Z"%v, "dd%sdd_X"%v, "dd%sdd_Y"%v, "dd%sdd_Z"%v, "dd%sdxdy"%v, "dd%sdxdz"%v, "dd%sdydz"%v]  
   errs = []
   #
   #
   DefineVectorExpression("coords_DOM_dual", "coord(DOM_dual)")
   DefineScalarExpression("x", "coord(DOM_dual)[0]")
   DefineScalarExpression("y", "coord(DOM_dual)[1]")
   DefineScalarExpression("z", "coord(DOM_dual)[2]")
   DefineScalarExpression("ANA_Ux", "10.*(1+1*cos(2.*3.141592653589793*x/2.+0.))*(1+1*cos(2.*3.141592653589793*y/3.+0.))*sin(3.141592653589793*z/5.)+(-1)")
   DefineScalarExpression("ANA_Uy", "7.*(1+1*cos(2.*3.141592653589793*x/2.+0.))*(1+1*cos(2.*3.141592653589793*y/3.+0.))*sin(3.141592653589793*z/5.)+(-1)")
   DefineScalarExpression("ANA_Uz", "3.*(1+1*cos(2.*3.141592653589793*x/2.+0.))*(1+1*cos(2.*3.141592653589793*y/3.+0.))*sin(3.141592653589793*z/5.)+(-1)")
   DefineScalarExpression("err_Ux", "abs(ANA_Ux-VELOCITY_X_FACES_DOM_dual)")
   DefineScalarExpression("err_Uy", "abs(ANA_Uy-VELOCITY_Y_FACES_DOM_dual)")
   DefineScalarExpression("err_Uz", "abs(ANA_Uz-VELOCITY_Z_FACES_DOM_dual)")
   suite = ["err_Ux", "err_Uy", "err_Uz"]
   #
   pass

for txt in l:
   ret = DefineScalarExpression("err_%s"%txt, "abs(ANA_%s_ELEM_DOM-%s_ELEM_DOM)"%(txt,txt))
   print ret
   errs.append("err_%s"%txt)
   pass

errs.extend(suite)
AddPlot("Pseudocolor", errs[0], 1, 1)
DrawPlots()

N = []
Linf = [] # Will become a 2D list...
L2 = [] # Will become a 2D list...
for i, db in enumerate(dbs):
   Linf.append([])
   L2.append([])
   N.append(int(db.split("/")[0].strip("N")))
   OpenDatabase("localhost:"+db, 0)
   ReplaceDatabase("localhost:"+db, 0)
   if i>0:
      CloseDatabase("localhost:"+dbs[i-1])
      pass
   for j, var in enumerate(errs):
      maxi, norm_l2 = myQuery(var)
      Linf[i].append(maxi)
      L2[i].append(norm_l2)
      pass
   pass

for nor in ["max","L2"]:
   f = open("cvg_%s_%s.txt"%(v,nor),'w')
   st = "# Nx\t"
   f.write(st)
   for txt in errs:
      f.write("%s\t "%txt)
      pass
   f.write("\n")
   for i, n in enumerate(N):
      f.write("%d\t"%n)
      for val in (Linf[i]):
         st = "%g\t"%(val)
         f.write(st)
         pass
      f.write("\n")
      pass
   f.close()
   pass


print Linf
print L2

exit(0)
