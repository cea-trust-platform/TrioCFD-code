import numpy as np
import os, re, shutil, glob, math, subprocess

refNx = 48 
refNy = 24
refdegliq = 50
refoffset = 0 # in microns
refDT = 8.5
loffsets = ['-0.5', '-0.4', '-0.3', '-0.2', '-0.1', '0.0', '0.1', '0.2', '0.3', '0.4', '0.5']
lraf = [0.25, 0.5, 1, 2, 4, 8]
langle = [5, 10, 15, 20, 30, 50, 80, 90]
ldT = [3,5,8.5,10,15]
tstep=1.e-6  # Reference timestep
lts = [1e-8, 1e-7, 2e-7, 3e-7] # List of Timesteps tested
lda=0.67897  # conductivity liq 

def makeCase(dest="TEST", Nx=refNx, Ny=refNy,degliq=refdegliq,offset=refoffset, dT=refDT, tstep=tstep):
   if os.path.isdir(dest):
      # shutil.rmtree(dest)
      print("Directory %s not updated"%dest)
      return
   else:
      os.mkdir(dest)
      pass
   # Read in the file and make the replacements:
   with open('template.data') as f:
      data = f.read()
   ym = 3e-7
   Ri=6.38209e-08
   theta_app=degliq*math.pi/180
   sm = ym/math.tan(theta_app)
   # Qi=lda*dT/theta_app*math.log(sm*theta_app/(Ri*lda)+1.)
   Qi=lda*dT/theta_app*math.log(sm*theta_app/(Ri*lda)+1.)
   data = re.sub(r'@Nx@',str(int(Nx+1)),data)
   data = re.sub(r'@Ny@',str(int(Ny+1)),data)
   data = re.sub(r'@degliq@',str(degliq),data)
   data = re.sub(r'@offset@',str(offset),data)
   data = re.sub(r'@DT@',str(dT),data)
   data = re.sub(r'@tstep@',str(tstep),data)
   data = re.sub(r'@sm@',str(sm),data)
   data = re.sub(r'@Qi@',str(Qi),data)
   #data = re.sub(r'(\d+):',lambda m: repl[m.group(1)]+':',data)
   #
   # Write it back out:
   with open(os.path.join(dest,'calc.data'),'w') as f:
      f.write(data)

   with open('postpro_flux.py') as f:
      fpost = f.read()
   fpost = re.sub(r'@offset@',str(offset),fpost)
   with open(os.path.join(dest,'postpro_flux.py'),'w') as f:
      f.write(fpost)
   
   #
   # Write it back out:
   file_info = os.path.join(dest,'info.txt')
   with open(file_info,'w') as f:
      f.write("lda=%g\n"%lda)
      f.write("DT=%s\n"%str(dT))
      f.write("theta=%s\n"%str(degliq))
      f.write("Lvap=%s\n"%str(-2.2574e6))
      f.write("offset=%s\n"%str(offset))
      f.write("Nx=%s\n"%str(Nx))
      f.write("sm=%s\n"%str(sm))
   #
   n=1
   print("\tCasTest %s calc %d"%(dest,1))
   shutil.copy("post_run", dest)
   # shutil.copy("data_extract.sh", dest)
   for nmeso in [2,4,6,8]:
      data2 = re.sub(r'deactivate',str(""),data) # To activate TCL
      data2 = re.sub(r'n_extend_meso 4',"n_extend_meso %d"%nmeso,data2) # To activate TCL
      dest2 = os.path.join(dest,"TCL%d"%nmeso)
      os.mkdir(dest2)
      with open(os.path.join(dest2,'calc.data'),'w') as f:
      	f.write(data2)
      with open(os.path.join(dest2,'postpro_flux.py'),'w') as f:
      	f.write(fpost)
      shutil.copy(file_info, dest2)
      shutil.copy("post_run", dest2)
      pass
   pass

# Offset serie : 
for o in loffsets : 
   makeCase(dest="CAS_OFFSET%s"%o, Nx=refNx, Ny=refNy,degliq=refdegliq,offset=o)
   pass 

# Instabilities with offset : 
for ts in lts : 
   for insta in ["-0.2", "0.2", "0.0"]:
      makeCase(dest="CAS_OFFSET%s_TIMESTEP%s"%(insta,ts), Nx=refNx, Ny=refNy,degliq=refdegliq,offset=insta)
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
for dT in ldT: 
   if (dT == refDT):
       if not os.path.exists("CAS_DT%s"%dT): subprocess.call(["ln","-sf", "REF", "CAS_DT%s"%dT])
   else:
      makeCase(dest="CAS_DT%s"%dT, Nx=refNx, Ny=refNy,degliq=refdegliq,offset=refoffset, dT=dT)
      pass
   pass 

exit()
