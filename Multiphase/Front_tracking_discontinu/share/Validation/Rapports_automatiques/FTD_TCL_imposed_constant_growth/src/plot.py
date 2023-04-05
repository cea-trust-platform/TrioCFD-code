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
tlata, min_vi, mean_vi, max_vi = get_vi()
if is_non_zero_file('sum_dIdt_or_dV.txt'):
  tt, dIdt, dV, dIdt_brm, dV_brm, Vlagrange = np.loadtxt('sum_dIdt_or_dV.txt').T
else:
  tt = dIdt = dV = dIdt_brm = dV_brm = Vlagrange = np.array([]) 

if is_non_zero_file("sum_dIdt_after-PCH.txt"):
  _, dIdt2 = np.loadtxt('sum_dIdt_after-PCH.txt').T
else:
  dIdt2 = np.array([])

show=False
if (len(sys.argv)==4):
   show=True
   pass

t = np.r_[0.,t]

mp=-10.
rhol=1000.
rhov=800.
r0=0.00012

# Analytical solution : 
r_ana = r0-mp/rhov*t
tlata = np.r_[0.,tlata]
vi_ana =-mp/rhov*np.ones(len(tlata)) # dr/dt

pwd=os.getcwd().split('/')
cas=pwd[-5]
if (cas=="2D_axi"):
   v0=4./3.*math.pi*pow(r0,3)
   v_ana=4./3.*math.pi*pow(r_ana,3)
   a0=4.*math.pi*pow(r0,2)
   a_ana=4.*math.pi*pow(r_ana,2)
elif (cas=="2D"):
   v0=math.pi*pow(r0,2)/2
   v_ana=math.pi*pow(r_ana,2)/2
   a0=math.pi*r0
   a_ana=math.pi*r_ana
   pass

v       = np.r_[v0,v]
ai      = np.r_[a0,ai]
min_vi  = np.r_[vi_ana[0],min_vi]
mean_vi = np.r_[vi_ana[0],mean_vi]
max_vi  = np.r_[vi_ana[0],max_vi]
v_theo  = v[0]-1./rhov*np.r_[0.,integrate.cumtrapz(mp*ai,t)] # the np.r_ is to put the integral at the end of each timestep.
if (len(tt)):
   # I add '0.' before the list of dV so that no dVolume is applied to the value v[0]; otherwise, vbis[0] would recieve the first dV 
   # integrated over the first timestep. 
   vbis = v[0]+np.r_[0.,integrate.cumtrapz(dIdt,tt)]
   vbis_brm = v[0]+np.r_[0.,integrate.cumtrapz(dIdt_brm,tt)]
   vlagrange = v[0]+np.r_[0.,integrate.cumtrapz(dIdt_brm,tt)]
   v_dIdt2 = v[0]-np.r_[0.,integrate.cumtrapz(dIdt2,tt)] # because it already has mp*ai/rhov = int_vol (mp*deltai/rhov) dv
else:
   vbis = vbis_brm = v_dIdt2 = vlagrange = np.array([])
   pass


remesh=pwd[-1].replace("Remesh","").replace("Nothing", "None").replace("Best", "Opt.")
interp=pwd[-2].replace("Interp","").replace("AiBased", "ai").replace("Standard", "Std.")
mass_source=pwd[-3].replace("MassSource","").replace("Historical", "Hist.")
vof=pwd[-4]
name="%s/%s - MS:%s - Int:%s - Rem:%s"%(cas,vof,mass_source,interp,remesh)
print("python dealing with ", name)
plt.figure(figsize=(10,5))
plt.suptitle("Comparisons to analytical solutions. Case %s"%(name))
plt.subplot(131)
plt.plot(t,v,label='resu')
plt.plot(tt,v_dIdt2,'-.',label='theory VoF (dIdt_after-PCH)')
#plt.plot(t,v_theo, 'k--', label='theo rebuilt from actual ai')
#plt.plot(tt,vbis_brm, '--', label='before ramasse-miette')
#plt.plot(tt,vbis, label='after ramasse-miette')
plt.plot(tt,vlagrange, '--', label='lagrangien')
plt.plot(t,v_ana, 'ko', markevery=10, label='ana')
plt.title("vapor volume")
plt.legend(loc=0)

plt.subplot(132)
plt.title("Interfacial velocity")
plt.plot(tlata, min_vi, ':', label='min')
plt.plot(tlata, mean_vi, '-',label='resu (mean)')
plt.plot(tlata, max_vi, ':', label='max')
plt.plot(tlata,vi_ana, 'ko', markevery=1, label='ana')
plt.legend(loc=0)

plt.subplot(133)
plt.title("ai")
plt.plot(t,ai, label='ai')
plt.plot(t,a_ana, 'ko', markevery=10,label='ana')
plt.legend(loc=0)
plt.savefig("plot.png")
if (show):
   plt.show()
