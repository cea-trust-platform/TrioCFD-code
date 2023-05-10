import numpy as np
import os, re, shutil, glob, math, subprocess

refNx = 50 
refNy = 24
refdegliq = 50
refoffset = 0 # in microns
refDT = 7.
loffsets = ['-0.5', '-0.4', '-0.3', '-0.2', '-0.1', '0.0', '0.1', '0.2', '0.3', '0.4', '0.5']
lraf = [0.25, 0.5, 1, 2, 4, 8]
langle = [5, 10, 15, 20, 30, 50, 80, 90]
ldt = [3,5,7,9,15]
tstep=1.e-6
lts = [0.1e-6, 0.2e-6, 0.3e-6]

def makeCase(dest="TEST", Nx=refNx, Ny=refNy,degliq=refdegliq,offset=refoffset, dt=refDT, tstep=tstep):
   if os.path.isdir(dest):
      shutil.rmtree(dest)
      pass
   os.mkdir(dest)
   # Read in the file and make the replacements:
   with open('template.data') as f:
      data = f.read()
   data = re.sub(r'@Nx@',str(int(Nx+1)),data)
   data = re.sub(r'@Ny@',str(int(Ny+1)),data)
   data = re.sub(r'@degliq@',str(degliq),data)
   data = re.sub(r'@offset@',str(offset),data)
   data = re.sub(r'@DT@',str(dt),data)
   data = re.sub(r'@tstep@',str(tstep),data)
   #data = re.sub(r'(\d+):',lambda m: repl[m.group(1)]+':',data)
   #
   # Write it back out:
   with open(os.path.join(dest,'calc.data'),'w') as f:
      f.write(data)
   #
   # Write it back out:
   with open(os.path.join(dest,'info.txt'),'w') as f:
      f.write("lda=0.67897\n")
      f.write("DT=%s\n"%str(dt))
      f.write("theta=%s\n"%str(degliq))
      f.write("Lvap=%s\n"%str(-2.2574e6))
   #
   n=1
   print("\tCasTest %s calc %d"%(dest,1))
   shutil.copy("post_run", dest)
   # shutil.copy("data_extract.sh", dest)
   pass

# Offset serie : 
for o in loffsets : 
   makeCase(dest="CAS_OFFSET%s"%o, Nx=refNx, Ny=refNy,degliq=refdegliq,offset=o)
   pass 

# Instabilities with offset : 
for ts in lts : 
   for insta in ["-0.2", "0.2"]:
      makeCase(dest="CAS_OFFSET%s_TIMESTEP%s"%(insta,ts), Nx=refNx, Ny=refNy,degliq=refdegliq,offset=o)
      pass
   pass 

if not os.path.exists("REF"): subprocess.call(["ln","-sf", "CAS_OFFSET0.0", "REF"])

# refinement serie : 
for r in lraf : 
   if (r == 1):
       if not os.path.exists("CAS_REFINE%s"%r): subprocess.call(["ln","-sf", "REF", "CAS_REFINE%s"%r])
   else:
      makeCase(dest="CAS_REFINE%s"%r, Nx=int(refNx*r), Ny=(refNy*r),degliq=refdegliq,offset=refoffset)
      pass
   pass 

# Angle serie : 
for angle in langle : 
   if (angle == refdegliq):
       if not os.path.exists("CAS_ANGLE%s"%angle): subprocess.call(["ln","-sf", "REF", "CAS_ANGLE%s"%angle])
   else:
      makeCase(dest="CAS_ANGLE%s"%angle, Nx=refNx, Ny=refNy,degliq=angle,offset=refoffset)
      pass
   pass 

# Overheat series : 
for dt in ldt: 
   if (dt == refDT):
       if not os.path.exists("CAS_DT%s"%dt): subprocess.call(["ln","-sf", "REF", "CAS_DT%s"%dt])
   else:
      makeCase(dest="CAS_DT%s"%dt, Nx=refNx, Ny=refNy,degliq=refdegliq,offset=refoffset, dt=dt)
      pass
   pass 

exit()
