# This is a python for the validation test case Beam_Free_Vibration.data

import numpy as np



t, Fpx = np.loadtxt('Beam_Free_Vibration_pb_Force_pression.out', unpack=True, usecols=[0,1])
Fvx = np.loadtxt('Beam_Free_Vibration_pb_Contrainte_visqueuse.out', unpack=True, usecols=[1])
    
Fx = Fvx + Fpx #force per unit length of cylinder

Fpy = np.loadtxt('Beam_Free_Vibration_pb_Force_pression.out', unpack=True, usecols=[2])
Fvy = np.loadtxt('Beam_Free_Vibration_pb_Contrainte_visqueuse.out', unpack=True, usecols=[2])

Fy = Fvy + Fpy #force per unit length of cylinder

Fpz = np.loadtxt('Beam_Free_Vibration_pb_Force_pression.out', unpack=True, usecols=[3])
Fvz = np.loadtxt('Beam_Free_Vibration_pb_Contrainte_visqueuse.out', unpack=True, usecols=[3])

Fz = Fvz + Fpz #force per unit length of cylinder

DataOut = np.column_stack((t,Fx, Fy, Fz))
np.savetxt('Numerical_force.txt', DataOut)


