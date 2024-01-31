import sys
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

if (len(sys.argv) != 4):
   raise Exception("One argument is compulsory to set offset, Nx, and angle")

print(f"Offset set to {sys.argv[1]} µm")
x_shift = float(sys.argv[1])*1.e-6

print(f"Number of node in x-dir set to {sys.argv[2]}")
Nx = int(sys.argv[2])

print(f"Contact angle set to {sys.argv[3]}")
angle = int(sys.argv[3])


# for offset 0 and -0.5, the Front will deplace into the left cell due to remaillage
# The start position will be smaller
if (x_shift == 0.0 or x_shift == -0.5e-6) and (Nx == 48) :
   # for some contact angle, the remaillage will glisse the Front...
   if (angle not in [5, 80, 90]) :
      x_start = 5e-6 - 0.05 * 1e-6
   else:
      x_start = 5e-6
else:
    x_start = 5e-6


print(x_start)
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
       shif_ref = float(s) if (5.e-6-x_start) > 0 else 0.
       mat[it,ix,0] = float(x)+shif_ref
       mat[it,ix,1] = y
       mat[it,ix,2] = s
       mat[it,ix,3] = phi if ( x_start+x_shift < float(x)+float(s)/2.) else 0.
       mat[it,ix,4] = pui if ( x_start+x_shift < float(x)+float(s)/2.) else 0.
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

s, phi, phi_ana = np.loadtxt("s_qi_qiana.txt").T


f = figure()
xlabel("s [m]")
ylabel("Flux [W/m²]")
plot(s,phi,'ro', label=r'qi')
plot(s,phi_ana,'k-', label=r'ana')
savefig("flux.png")

close(f)
show()

