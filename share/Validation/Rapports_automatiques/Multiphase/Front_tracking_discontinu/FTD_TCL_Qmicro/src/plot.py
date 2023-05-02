import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import os, math, sys
def is_non_zero_file(fpath):  
    return os.path.isfile(fpath) and os.path.getsize(fpath) > 0

def get_vi():
    res = np.ones((0,4))
    with open('vifiles.txt') as f:
        for fvit in f.readlines():
           fvit=fvit.strip()
           t=float(fvit.replace("post-ascii.lata.VITESSE.SOM.INTERFACES.",""))
           mat=np.genfromtxt(os.path.join("lata",fvit),skip_header=1, skip_footer=1).T
           if (np.shape(mat) == (0,)):
               res=np.r_[res,np.array([[t,0,0,0]])]
               # Empty no interface:
           else:
               #import pdb; pdb.set_trace()
               v=np.sqrt(mat[0,:]*mat[0,:]+mat[1,:]*mat[1,:])
               res=np.r_[res,np.array([[t,v.min(),v.mean(),v.max()]])]
               pass
           pass
    return res[:,0], res[:,1], res[:,2], res[:,3]

t, v = np.loadtxt('vol_vap.txt').T
_, ai = np.loadtxt('ai.txt').T
ti, rmin, rmax, xmax, ymax, ymin, rcl = np.loadtxt('interface-position.txt').T
tlata, min_vi, mean_vi, max_vi = get_vi()
if is_non_zero_file('sum_dIdt_or_dV.txt'):
  tt, dIdt, dV, dIdt_brm, dV_brm, Vlagrange = np.loadtxt('sum_dIdt_or_dV.txt').T
else:
  tt = dIdt = dV = dIdt_brm = dV_brm = Vlagrange = np.array([]) 

if is_non_zero_file("sum_dIdt_after-PCH.txt"):
  _, dIdt2 = np.loadtxt('sum_dIdt_after-PCH.txt').T
else:
  dIdt2 = np.array([])

Qcl = 50. # J/(s*m)
Lvap = -2.256e6 # J/kg
rhov=800.  # kg/m3

# Analytical solution : 
ti = np.r_[0.,ti]
Gcl = Qcl* 2.*math.pi *rcl / Lvap # kg/(s)
dVbdt_ana = Gcl/rhov # m3/(s)
v_ana  = v[0]-np.r_[0.,integrate.cumtrapz(dVbdt_ana,ti[:-1])] # the np.r_ is to put the integral at the end of each
ti  = ti[:-1] # to delete the n+1

np.savetxt("v_ana.txt", np.c_[ti,v_ana/v[0]])
np.savetxt("v_num.txt", np.c_[t,v/v[0]])

cwd = os.getcwd()
name=os.path.basename(cwd)
kk=cwd.split('/')[-2]
print("python dealing with ", kk, "/", name)
plt.figure(figsize=(10,5))
plt.suptitle("Comparison to analytical solutions. Case %s / %s"%(kk,name))
plt.subplot(131)
plt.plot(t,v, 'r-', label='Simu.')
plt.plot(ti,v_ana, 'ko', markevery=1, label='ana')
plt.xlabel('time [s]')
plt.title("vapor volume")
plt.legend(loc=0)

plt.subplot(132)
plt.title("Interfacial velocity")
plt.plot(tlata, min_vi, ':', label='min')
plt.plot(tlata, mean_vi, '-',label='resu (mean)')
plt.plot(tlata, max_vi, ':', label='max')
#plt.plot(tlata,vi_ana, 'ko', markevery=1, label='ana')
plt.legend(loc=0)

plt.subplot(133)
plt.title("contact line")
plt.plot(ti,rcl, label='rcl')
plt.legend(loc=0)
plt.savefig("plot.png")
plt.close()

if os.path.isfile("results.dat"):
   _, time, xcl, mpai, Vb, dVbdt_from_mpai, dVbdt_ana = np.loadtxt("results.dat").T
   tv = 0.5*(time[1:]+time[:-1])
   dVbdt = (Vb[1:]-Vb[:-1])/(time[1:]-time[:-1])
   plt.figure(figsize=(8,5))
   plt.suptitle("Comparison to analytical solutions. Case %s / %s"%(kk,name))
   plt.plot(tv,dVbdt, 'r-', label='Simu.')
   plt.plot(time,-dVbdt_from_mpai, 'bx--', markevery=1, label=r'Simu. (from $\dot{m}$)')
   plt.plot(time,-dVbdt_ana, 'ko', markevery=1, label='ana')
   plt.xlabel('time [s]')
   plt.ylabel(r"$\frac{d(V_b)}{dt}$ [m3/s]")
   plt.legend(loc=0)
   plt.savefig("plot-dVdt.png")
   pass
   
#plt.show()
