#***********************Case test : Force and Coeff. Num.*********************
#******************************************************************************

# Tested on Python 2.7.16

# Import library
import numpy as np

rho = 1000.0 # fluid density
nu = 1.007e-6 # fluid viscosity

# Definition of diameters
R1 = 1.0 # Rayon left static cylinder
R2 = 1.0 # Rayon right moving cylinder

# Definition of displacement
KC = 1.0e-2 # Keulegan-Carpenter number KC= u/R2
u = KC * R2 # module of displacement with u < 1e-2
Sk = 900.0 # Stokes number
f = nu*Sk/(2*np.pi) # frequency of displacement (s^-1)
T = 1.0/f # period of displacement
Omega = Sk * nu / (R2**2.)   # angular frequency of displacement (rad/s)

print('Stokes number Sk =', Sk)
print('Frequency  F[1/s] =', f)
print('Pulsation w[rad/s] =', Omega)
print('Period T[s] =', T)



t = np.loadtxt('TwoOscillatingCylinders_pb_Force_pression.out', unpack=True, usecols=[0]) # time

t_adim = Omega*t # dimensionless time

# Numerical

Fpx_left = np.loadtxt('TwoOscillatingCylinders_pb_Force_pression.out', unpack=True, usecols=[3]) # pressure x-axis for left cylinder
Fvx_left = np.loadtxt('TwoOscillatingCylinders_pb_Contrainte_visqueuse.out', unpack=True, usecols=[3]) # viscous force x-axis for left cylinder

Fpx_right = np.loadtxt('TwoOscillatingCylinders_pb_Force_pression.out', unpack=True, usecols=[5]) # pressure x-axis for right cylinder
Fvx_right = np.loadtxt('TwoOscillatingCylinders_pb_Contrainte_visqueuse.out', unpack=True, usecols=[5]) # viscous force x-axis for right cylinder

Fx_left = Fvx_left + Fpx_left # total fluid force x-axis for left cylinder
Fx_right = Fvx_right + Fpx_right # total fluid force x-axis for right cylinder


Fx_left_adim = Fx_left / (rho*R1**2.*u*Omega**2.) # total dimensionless fluid force x-axis for left cylinder
Fx_right_adim = Fx_right / (rho*R1**2.*u*Omega**2.) # total dimensionless fluid force x-axis for right cylinder

DataOut_left = np.column_stack((t_adim,Fx_left_adim))
np.savetxt('Numerical_force_left.txt', DataOut_left) #self

DataOut_right = np.column_stack((t_adim,Fx_right_adim))
np.savetxt('Numerical_force_right.txt', DataOut_right) # cross

# Determining m_self, m_cross, c_cross and c_self (numerically)
# Fourier's inner product : <f,g> = (2/5T) * integral from 0 to 5T f(t)*g(t) dt
# m_self_x = (2/5T) * integral from 0 to 5T sin(Omega*t) Fx_right(t) dt/rho*pi*R1^2*u*Omega^2

m_self_x = 2./5./T*np.trapz(np.sin(Omega*t)*Fx_right,t)/(rho*np.pi*R1**2.*u*Omega**2.) 
c_self_x = -2./5./T*np.trapz(np.cos(Omega*t)*Fx_right,t)/(rho*np.pi*R1**2.*u*Omega**2.) 
m_cross_x = 2./5./T*np.trapz(np.sin(Omega*t)*Fx_left,t)/(rho*np.pi*R1**2.*u*Omega**2.) 
c_cross_x = -2./5./T*np.trapz(np.cos(Omega*t)*Fx_left,t)/(rho*np.pi*R1**2.*u*Omega**2.) 

DataOut1 = np.column_stack((m_self_x,c_self_x, m_cross_x, c_cross_x))
np.savetxt('Numerical_coefficients.txt', DataOut1)

# Theory 


JCP_m_self = 1.10
LS_m_self = 1.11
COL_m_self = 1.11

JCP_c_self = 0.117
LS_c_self = 0.105
COL_c_self = 0.106

JCP_m_cross = -0.116
LS_m_cross = -0.138
COL_m_cross = -0.138

JCP_c_cross = -0.0135
LS_c_cross = -0.0132
COL_c_cross = -0.0136

DataOut2 = np.column_stack((COL_m_self,COL_c_self, COL_m_cross, COL_c_cross))
np.savetxt('COL_coefficients.txt',  DataOut2)
DataOut3 = np.column_stack((LS_m_self,LS_c_self, LS_m_cross, LS_c_cross))
np.savetxt('LS_coefficients.txt',  DataOut3)

Fx_left_adim_LS = (LS_m_cross*np.sin(t_adim)-LS_c_cross*np.cos(t_adim))*np.pi
Fx_right_adim_LS = (LS_m_self*np.sin(t_adim)-LS_c_self*np.cos(t_adim))*np.pi

DataOut_left_LS = np.column_stack((t_adim,Fx_left_adim_LS))
np.savetxt('LS_force_left.txt', DataOut_left_LS)
DataOut_right_LS = np.column_stack((t_adim,Fx_right_adim_LS))
np.savetxt('LS_force_right.txt', DataOut_right_LS)

Fx_left_adim_COL = (COL_m_cross*np.sin(t_adim)-COL_c_cross*np.cos(t_adim))*np.pi
Fx_right_adim_COL = (COL_m_self*np.sin(t_adim)-COL_c_self*np.cos(t_adim))*np.pi

DataOut_left_COL = np.column_stack((t_adim,Fx_left_adim_COL))
np.savetxt('COL_force_left.txt', DataOut_left_COL)
DataOut_right_COL = np.column_stack((t_adim,Fx_right_adim_COL))
np.savetxt('COL_force_right.txt', DataOut_right_COL)
