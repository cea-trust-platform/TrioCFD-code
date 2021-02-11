#*********************** Force and Coeff. Num.*********************
#******************************************************************************

# Import library
import numpy as np

test_case = "3dOscillatingBeam_pb_"

rho = 1000. #fluid density
nu = 1.007e-6 #fluid viscosity

# Definition of displacement
Lmax  = 0.7 ; #length of the beam
Ri=0.02
u= 0.5*1.e-4; #module of displacement 
w = 1000.9   #pulsation of displacement
f=w/(2.*np.pi) #frequency of displacement
T=1./f #period of displacement
nb_p=5.; # number of periods


# Calculation of the fluid force
t= np.loadtxt(test_case+'Force_pression.out', unpack=True, usecols=[0])

Fpy = np.loadtxt(test_case+'Force_pression.out', unpack=True, usecols=[8])
Fvy = np.loadtxt(test_case+'Contrainte_visqueuse.out', unpack=True, usecols=[8])
Fy = Fvy + Fpy 

DataOut = np.column_stack((t,Fy))
np.savetxt('Numerical_force.txt', DataOut)




# Compute the fluid added mass and fluid added damping coeff.



CMNumy=2./nb_p/T*np.trapz(np.sin(w*t)*Fy,t)
CVNumy =- 2./nb_p/T*np.trapz(np.cos(w*t)*Fy,t)


CmNumy = CMNumy/(Lmax*rho*Ri**2.*np.pi*u*w**2.)
CvNumy = CVNumy/(Lmax*rho*Ri**2.*np.pi*u*w**2.)

# print('Added mass coefficient Cmy =', CmNumy)
# print('Added dumping coefficient Cvy =', CvNumy)

DataOut1 = np.column_stack((CmNumy,CvNumy))
np.savetxt('Numerical_coefficients.txt', DataOut1)


# Save DataOut in .txt file

Cm=10.16
Cv=0.
F = -rho*np.pi*Lmax*Ri**2.*u*w**2.*(-Cm*np.sin(w*t) + Cv*np.cos(w*t)) #Force vector
DataOut = np.column_stack((t,F))
np.savetxt('Theoretical_force.txt', DataOut) #this command overwrite the previus .txt file


