# -*- coding : utf8
import readline, rlcompleter, math
readline.parse_and_bind('tab:complete')
import DNSTools as dtool

def sciform(x):
    return '{:.6E}'.format(x)

output="/tmp/deplacements.txt"
bubble_name="bulle_r0.0005.msh"
if (len(sys.argv) != 3):
   raise Exception("usage: python %s jdd.data alpha_wanted\n Datafile is required to read the domain and bubble volume.\nOutput is in %s"%(sys.argv[0],output))

jdd=sys.argv[1]
al=(float)(sys.argv[2])
print "Reading domain information from: ", jdd
if (not os.path.isfile(jdd)): raise Exception("jdd %s is missing!"%jdd)
# Donnees utilisateur
Lx=dtool.getParam(jdd, 'uniform_domain_size_i')
Ly=dtool.getParam(jdd, 'uniform_domain_size_j')
Lz=dtool.getParam(jdd, 'uniform_domain_size_k')
Nx = (int)(dtool.getParam(jdd, 'nbelem_i'))
Ny = (int)(dtool.getParam(jdd, 'nbelem_j'))
Nz = (int)(dtool.getParam(jdd, 'nbelem_k'))
vb = (float)(dtool.getParam(jdd, 'vol_bulle_monodisperse'))
nb =(int)(math.floor(al*Lx*Ly*Lz/vb))
periok=dtool.getParam(jdd, 'perio_k',compo=-1) # compo=-1 for a flag.
print "k-periodic? ", periok
#nb=2; al=nb*vb/(Lx*Ly*Lz); print(al); exit(0)
rb=pow(vb/(4./3.*math.pi), 1./3.)
print("The script will generate %d bubbles of volume %g[m^3] (rb=%g) to reach a void fraction alpha=%g"%(nb,vb,rb,al))

eta=math.pi/(3.*math.sqrt(2.))
if (al>eta): raise Exception("Void fraction %g above maximum achievable Bubble packing (%g)"%(al,eta))
# Trying to do FCC packing (able to reach up to 74%)
# https://mathworld.wolfram.com/CubicClosePacking.html

text = ""
count = 0

eps=0.01
delta_min=(1.+eps)*rb*2.*math.sqrt(2.)

vbmax=eta/al*vb
rmax=pow(vbmax/(4./3.*math.pi), 1./3.)
delta_max=(1.-eps)*rmax*2.*math.sqrt(2.)
fac=1.0 # between 0 and 1. 
delta=delta_min+(fac)*(delta_max-delta_min) 
delta*=1.0 # 1.12805 # Etrange que l'on puisse depasser delta_max! 
print("Bubble spacing: %g in range[%g:%g]"%(delta, delta_min, delta_max))
nx=(int)(math.floor(Lx/delta))
ny=(int)(math.floor(Ly/delta))
nz=(int)(math.floor(Lz/delta))
print("Setting bubbles on FCC network : %dx%dx%d -> %d bubbles (2 bubbles per FCC)"%(nx,ny,nz,2*nx*ny*nz)) 

lim=math.ceil(4.*math.sqrt(2)) # Pourquoi 2sqrt(2) pas suffisant?
a=-math.ceil(lim*max(ny,nz))
b=math.ceil(lim*max(ny,nz))
irange = range((int)(a),(int)(b+1))
print("loop i %d %d"%(a,b+1))
c=-math.ceil(lim*max(nx,nz))
d=math.ceil(lim*max(nx,nz))
print("loop j %d %d"%(c,d+1))
cc=-math.ceil(lim*max(nx,ny))
dd=math.ceil(lim*max(nx,ny))
print("loop k %d %d"%(cc,dd+1))
for k in range((int)(cc),(int)(dd+1)):
   for j in range((int)(c),(int)(d+1)):
      for i in irange:
         x=delta*i
         y=delta*j
         z=delta*k
         xb=(y+z)/2.
         yb=(x+z)/2.
         zb=(x+y)/2.
         # if center in box:
	 dwall=0.
	 if (not periok):
	    dwall = (1.+eps)*rb
	    pass
	 #
         if (xb>=0 and yb >=0 and zb >= dwall and \
             xb<Lx and yb <Ly and zb <Lz-dwall):
            st="%g %g %g %s"%(xb,yb,zb,bubble_name)
            #print(st)
            text+=st+"\n"
            count+=1
            if (count>=nb): 
               print("total number of bubbles %d reached"%nb)
               break
            pass
         pass
      if (count>=nb): break
      pass
   if (count>=nb): break
   pass

if (count<nb):
   print("last bubble %s"%text)
   raise Exception("All bubbles (%d) have been positionned and still nb=%d is not reached (end of loops have not been checked)"%(count, nb))

if (count==((nx+1)*(ny+1)*(nz))):
   raise Exception("The last bubble (%d) has been positionned on an existing bubble via periodic bounds, right? (it has not been checked)"%(count))

with open(output, "w") as fic:
   fic.write(text)

print("see file %s"%(output))
print("The end")
exit(0)

