import numpy as np
from decimal import *

f_air = np.loadtxt('CF_Air/Freq.txt', unpack=True, usecols=[0]) 
f_water=np.loadtxt('CF_Water/Freq.txt', unpack=True, usecols=[0]) 


DataOut = np.column_stack((f_air, f_water))
np.savetxt('Freq_CF.txt', DataOut)

f_air = np.loadtxt('PP_Air/Freq.txt', unpack=True, usecols=[0]) 
f_water=np.loadtxt('PP_Water/Freq.txt', unpack=True, usecols=[0]) 

DataOut = np.column_stack((f_air, f_water))
np.savetxt('Freq_PP.txt', DataOut)

