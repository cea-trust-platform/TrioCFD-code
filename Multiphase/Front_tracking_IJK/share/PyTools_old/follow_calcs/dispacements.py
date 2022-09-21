import numpy as np 
import pylab as plt
import sys, os
Lx=Ly=0.005
Lz=0.02
print("usage: python %s <jdd>"%sys.argv[0])
if (len(sys.argv)==2):
   head=sys.argv[1].replace(".data","")
else:
   head="DNS"
   pass

f='%s_dmin.out'%head
if (not os.path.isfile(f)):
   raise Exception("Check usage. Missing file %s."%f)
time,dmin = np.loadtxt(f, usecols=(1,2)).T
xx = np.loadtxt('%s_bulles_centre_x.out'%head)
yy = np.loadtxt('%s_bulles_centre_y.out'%head)
zz = np.loadtxt('%s_bulles_centre_z.out'%head)

tt = xx[:,0]
# Compute bubble displacement : 
dx = xx[:,1:]-xx[0,1:]
dy = yy[:,1:]-yy[0,1:]
dz = zz[:,1:]-zz[0,1:]

# Translate de L si besoin:
dx[dx>Lx/2.]  -=Lx
dx[dx<-Lx/2.] +=Lx
dy[dy>Ly/2.]  -=Ly
dy[dy<-Ly/2.] +=Ly
dz[dz>Lz/2.]  -=Lz
dz[dz<-Lz/2.] +=Lz

deplacements = np.sqrt(dx*dx+dy*dy+dz*dz)
ddmin = deplacements.min(axis=1)
ddmax = deplacements.max(axis=1)
print ddmin

plt.plot(time, dmin, label="dmin")
plt.plot(tt,ddmin, label="displacement min")
plt.plot(tt,ddmax, label="displacement max")

ms = 0.005/128.
for i in range(2,8): 
   plt.plot([tt[0],tt[-1]], [i*ms, i*ms], 'k--')
plt.plot([tt[0],tt[-1]], [ms, ms], 'k--', label='dx')
plt.plot([0.012, 0.012], [3*ms, 5*ms], 'ko', label='sauv')
plt.grid(True)
plt.legend(loc=0)

plt.savefig("displacement.png")
plt.show()
