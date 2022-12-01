# This is a python script for the validation test case DIVA.data

# Import library
import numpy as np

rho = 1000. #fluid density
nu = 1.007e-6 #fluid viscosity

# Definition of diameter
d = 0.03

# Definition of displacement
u = 0.1 #module of dimensionless displacement (KC = u = U/d) 
f = 20. #frequency of displacement
T = 1./f #period of displacement
w = 2.*np.pi*f   #angular frequency of displacement

t = np.loadtxt('DIVA_pb_Force_pression.out', unpack=True, usecols=[0]) 
Fpx_centre, Fpy_centre, Fpx_nord, Fpy_nord, Fpx_sud, Fpy_sud, Fpx_est, Fpy_est, Fpx_ouest, Fpy_ouest = np.loadtxt('DIVA_pb_Force_pression.out', unpack=True, usecols=[11, 12, 13, 14, 15, 16, 17, 18, 19, 20])
Fvx_centre, Fvy_centre, Fvx_nord, Fvy_nord, Fvx_sud, Fvy_sud, Fvx_est, Fvy_est, Fvx_ouest, Fvy_ouest = np.loadtxt('DIVA_pb_Contrainte_visqueuse.out', unpack=True, usecols=[11, 12, 13, 14, 15, 16, 17, 18, 19, 20])

Fx_centre = Fpx_centre + Fvx_centre  
Fy_centre = Fpy_centre + Fvy_centre 

Fx_nord = Fpx_nord + Fvx_nord  
Fy_nord = Fpy_nord + Fvy_nord

Fx_sud = Fpx_sud + Fvx_sud 
Fy_sud = Fpy_sud + Fvy_sud

Fx_est = Fpx_est + Fvx_est
Fy_est = Fpy_est + Fvy_est

Fx_ouest = Fpx_ouest + Fvx_ouest
Fy_ouest = Fpy_ouest + Fvy_ouest

mx_self_num =  2./5./T*np.trapz(np.sin(w*t)*Fx_centre,t)/(rho*np.pi*(d/2.)**2.*u*d*w**2.)
cx_self_num = -2./5./T*np.trapz(np.cos(w*t)*Fx_centre,t)/(rho*np.pi*(d/2.)**2.*u*d*w**2.)
my_self_num =  2./5./T*np.trapz(np.sin(w*t)*Fy_centre,t)/(rho*np.pi*(d/2.)**2.*u*d*w**2.)
cy_self_num = -2./5./T*np.trapz(np.cos(w*t)*Fy_centre,t)/(rho*np.pi*(d/2.)**2.*u*d*w**2.)

mx_nord_num =  2./5./T*np.trapz(np.sin(w*t)*Fx_nord,t)/(rho*np.pi*(d/2.)**2.*u*d*w**2.)
cx_nord_num = -2./5./T*np.trapz(np.cos(w*t)*Fx_nord,t)/(rho*np.pi*(d/2.)**2.*u*d*w**2.)
my_nord_num =  2./5./T*np.trapz(np.sin(w*t)*Fy_nord,t)/(rho*np.pi*(d/2.)**2.*u*d*w**2.)
cy_nord_num = -2./5./T*np.trapz(np.cos(w*t)*Fy_nord,t)/(rho*np.pi*(d/2.)**2.*u*d*w**2.)

mx_sud_num =  2./5./T*np.trapz(np.sin(w*t)*Fx_sud,t)/(rho*np.pi*(d/2.)**2.*u*d*w**2.)
cx_sud_num = -2./5./T*np.trapz(np.cos(w*t)*Fx_sud,t)/(rho*np.pi*(d/2.)**2.*u*d*w**2.)
my_sud_num =  2./5./T*np.trapz(np.sin(w*t)*Fy_sud,t)/(rho*np.pi*(d/2.)**2.*u*d*w**2.)
cy_sud_num = -2./5./T*np.trapz(np.cos(w*t)*Fy_sud,t)/(rho*np.pi*(d/2.)**2.*u*d*w**2.)

mx_est_num =  2./5./T*np.trapz(np.sin(w*t)*Fx_est,t)/(rho*np.pi*(d/2.)**2.*u*d*w**2.)
cx_est_num = -2./5./T*np.trapz(np.cos(w*t)*Fx_est,t)/(rho*np.pi*(d/2.)**2.*u*d*w**2.)
my_est_num =  2./5./T*np.trapz(np.sin(w*t)*Fy_est,t)/(rho*np.pi*(d/2.)**2.*u*d*w**2.)
cy_est_num = -2./5./T*np.trapz(np.cos(w*t)*Fy_est,t)/(rho*np.pi*(d/2.)**2.*u*d*w**2.)

mx_ouest_num =  2./5./T*np.trapz(np.sin(w*t)*Fx_ouest,t)/(rho*np.pi*(d/2.)**2.*u*d*w**2.)
cx_ouest_num = -2./5./T*np.trapz(np.cos(w*t)*Fx_ouest,t)/(rho*np.pi*(d/2.)**2.*u*d*w**2.)
my_ouest_num =  2./5./T*np.trapz(np.sin(w*t)*Fy_ouest,t)/(rho*np.pi*(d/2.)**2.*u*d*w**2.)
cy_ouest_num = -2./5./T*np.trapz(np.cos(w*t)*Fy_ouest,t)/(rho*np.pi*(d/2.)**2.*u*d*w**2.)

DataOut = np.column_stack((t, Fx_centre))
np.savetxt('Numerical_force.txt', DataOut)

DataOut1 = np.column_stack((mx_self_num, cx_self_num, my_self_num, cy_self_num)) 
np.savetxt('Numerical_self_coefficients.txt', DataOut1)

DataOut2 = np.column_stack((mx_nord_num, cx_nord_num, my_nord_num, cy_nord_num)) 
np.savetxt('Numerical_nord_coefficients.txt', DataOut2)

DataOut3 = np.column_stack((mx_sud_num, cx_sud_num, my_sud_num, cy_sud_num)) 
np.savetxt('Numerical_sud_coefficients.txt', DataOut3)

DataOut4 = np.column_stack((mx_est_num, cx_est_num, my_est_num, cy_est_num)) 
np.savetxt('Numerical_est_coefficients.txt', DataOut4)

DataOut5 = np.column_stack((mx_ouest_num, cx_ouest_num, my_ouest_num, cy_ouest_num)) 
np.savetxt('Numerical_west_coefficients.txt', DataOut5)
