# This is a python for the validation test case TwoCylinders.data

import numpy as np

rho = 1000. #fluid density
nu = 1.007e-6 #fluid viscosity

# Definition of diameters
D1 = 1. # inner diameter D_1
D2 = 2. # outer diameter D_2

# Definition of displacement
u = 1.0e-2 # dimensionless module of displacement  (with respect to the D_1)
Sk = 10000. #Stokes number
f = nu*Sk/D1**2. #frequency of displacement
T = 1./f #period of displacement
w = 2.*np.pi*f   #angular frequency of displacement

t, Fpx = np.loadtxt('TwoCylinders_pb_Force_pression.out', unpack=True, usecols=[0,3])
Fvx = np.loadtxt('TwoCylinders_pb_Contrainte_visqueuse.out', unpack=True, usecols=[3])
    
Fx = Fvx + Fpx #force per unit length of cylinder
m_self = 2./5./T*np.trapz(np.sin(w*t)*Fx,t)/(rho*np.pi*(D1/2.)**2.*u*D1*w**2.)
c_self = -2./5./T*np.trapz(np.cos(w*t)*Fx,t)/(rho*np.pi*(D1/2.)**2.*u*D1*w**2.)

DataOut = np.column_stack((t,Fx))
np.savetxt('Numerical_force.txt', DataOut)

DataOut1 = np.column_stack((m_self,c_self))
np.savetxt('Numerical_coefficients.txt', DataOut1)
