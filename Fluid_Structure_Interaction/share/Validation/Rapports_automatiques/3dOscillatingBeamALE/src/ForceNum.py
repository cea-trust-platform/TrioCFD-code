#*********************** force and coeff. num.*********************
#******************************************************************************

# import library
import numpy as np

test_case = "3doscillatingbeamALE_pb_"

rho = 1000. #fluid density
nu = 1.007e-6 #fluid viscosity

# definition of displacement
lmax  = 0.7 ; #length of the beam
ri=0.02
u= 0.5*1.e-4; #module of displacement 
w = 1000.9   #pulsation of displacement
f=w/(2.*np.pi) #frequency of displacement
t=1./f #period of displacement
nb_p=5.; # number of periods


# calculation of the fluid force
t= np.loadtxt(test_case+'force_pression.out', unpack=true, usecols=[0])

fpy = np.loadtxt(test_case+'force_pression.out', unpack=true, usecols=[8])
fvy = np.loadtxt(test_case+'contrainte_visqueuse.out', unpack=true, usecols=[8])
fy = fvy + fpy 

dataout = np.column_stack((t,fy))
np.savetxt('numerical_force.txt', dataout)




# compute the fluid added mass and fluid added damping coeff.



cmnumy=2./nb_p/t*np.trapz(np.sin(w*t)*fy,t)
cvnumy =- 2./nb_p/t*np.trapz(np.cos(w*t)*fy,t)


cmnumy = cmnumy/(lmax*rho*ri**2.*np.pi*u*w**2.)
cvnumy = cvnumy/(lmax*rho*ri**2.*np.pi*u*w**2.)

# print('added mass coefficient cmy =', cmnumy)
# print('added dumping coefficient cvy =', cvnumy)

dataout1 = np.column_stack((cmnumy,cvnumy))
np.savetxt('numerical_coefficients.txt', dataout1)


# save dataout in .txt file

cm=10.16
cv=0.
f = -rho*np.pi*lmax*ri**2.*u*w**2.*(-cm*np.sin(w*t) + cv*np.cos(w*t)) #force vector
dataout = np.column_stack((t,f))
np.savetxt('theoretical_force.txt', dataout) #this command overwrite the previus .txt file


