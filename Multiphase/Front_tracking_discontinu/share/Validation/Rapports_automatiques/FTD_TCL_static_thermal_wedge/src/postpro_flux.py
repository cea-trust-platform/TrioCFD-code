import numpy as np
from numpy import *
from pylab import *
f = open("calc_pb_Diffusion_chaleur.face")
lines = f.readlines()
nb_lines = len(lines)

cnt =0
for ibl, line in enumerate(lines):
   if 'temps' in line:
      cnt+=1
      if cnt==2:
         break;

nx=ibl
nt=int(nb_lines/ibl)

mat = zeros((nt,nx,5))
it = -1
ix = 0
t = []
for line in lines:
   if 'temps' in line:
       tt = float(line.split('temps')[1].replace(":", '').strip())
       t.append(tt)
       it +=1
       ix = 0
       #print "Reading time t=", tt
   else:
       li = line.split(' ')
       x = li[4]
       y= li[6]
       s= li[8]
       phi=li[10] # flux_par_surface(W/m2)
       pui=li[12] # flux(W)
       mat[it,ix,0] = x
       mat[it,ix,1] = y
       mat[it,ix,2] = s
       mat[it,ix,3] = phi
       mat[it,ix,4] = pui
       ix +=1
       pass
   pass

print("Number of timesteps was correct: ", nt == len(t), " (nt=%d)"%nt)

xx = mat[0,:,0]
savetxt('xx.txt', xx)
savetxt('mphit.txt', mat[:,:,3])
savetxt('xphi.txt', np.c_[xx,mat[-1,:,3]])

cum=mat*0
cumsum(mat, axis=1, out=cum)
#plot(xx,cum[-1,:,4])

itdemi=int(nt/2)
print("saving medium step : t_medium=",t[itdemi])
savetxt('xphi_half.txt', np.c_[xx,mat[itdemi,:,3],cum[itdemi,:,4]])
print("saving last step : t_end=",t[nt-1])
savetxt('xphi_last.txt', np.c_[xx,mat[-1,:,3],cum[-1,:,4]])

f = figure()

for it in range(nt):
   plot(mat[it,:,0],mat[it,:,3],'k-', label=r'$t=%g'%t[it])
   pass

legend(loc=0)
# f.savefig("fig.png")
close(f)
show()

