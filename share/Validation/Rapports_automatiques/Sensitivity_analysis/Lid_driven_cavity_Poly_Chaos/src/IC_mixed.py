# Sensitivity
test_case_v = "Velocity/Lid_driven_cavity_Poly_Chaos_"
test_case_mu = "Mu/Lid_driven_cavity_Poly_Chaos_"
variable = "U"
direction = "y" #used only for vectorial variables (i.e. velocity). x, y or z
dir=2
cross_sec = "H"
confidence = 0.95 #level of confidence for the interval
std_dev_v = 0.1 #std deviation of the input
std_dev_mu = 0.001 #std deviation of the input


import numpy as np
import matplotlib.pyplot as plt

########## Sensitivity section  Velocity ##########

dx= np.loadtxt(test_case_v+variable+'_'+cross_sec+'.coupe', unpack=True, usecols=[0])
state_v= np.loadtxt(test_case_v+variable+'_'+cross_sec+'.coupe', unpack=True, usecols=[dir])
sens_v= np.loadtxt(test_case_v+variable+'_'+cross_sec+'_SENS.coupe', unpack=True, usecols=[dir])


state_mu= np.loadtxt(test_case_mu+variable+'_'+cross_sec+'.coupe', unpack=True, usecols=[dir])
sens_mu= np.loadtxt(test_case_mu+variable+'_'+cross_sec+'_SENS.coupe', unpack=True, usecols=[dir])



variance =np.abs(sens_v*std_dev_v + sens_mu*std_dev_mu)



upper_bound = state_v+1/np.sqrt(1-confidence)*variance
lower_bound = state_v-1/np.sqrt(1-confidence)*variance

DataOut = np.column_stack((dx,state_v, upper_bound, lower_bound))
np.savetxt('IC_mixed_'+cross_sec+'_'+direction+'.txt', DataOut)



variable = "U"
direction = "x" #used only for vectorial variables (i.e. velocity). x, y or z
dir=1
cross_sec = "H"
confidence = 0.95 #level of confidence for the interval
std_dev_v = 0.1 #std deviation of the input
std_dev_mu = 0.001 #std deviation of the input



########## Sensitivity section  Velocity ##########

dx= np.loadtxt(test_case_v+variable+'_'+cross_sec+'.coupe', unpack=True, usecols=[0])
state_v= np.loadtxt(test_case_v+variable+'_'+cross_sec+'.coupe', unpack=True, usecols=[dir])
sens_v= np.loadtxt(test_case_v+variable+'_'+cross_sec+'_SENS.coupe', unpack=True, usecols=[dir])


state_mu= np.loadtxt(test_case_mu+variable+'_'+cross_sec+'.coupe', unpack=True, usecols=[dir])
sens_mu= np.loadtxt(test_case_mu+variable+'_'+cross_sec+'_SENS.coupe', unpack=True, usecols=[dir])



variance =np.abs(sens_v*std_dev_v + sens_mu*std_dev_mu)



upper_bound = state_v+1/np.sqrt(1-confidence)*variance
lower_bound = state_v-1/np.sqrt(1-confidence)*variance

DataOut = np.column_stack((dx,state_v, upper_bound, lower_bound))
np.savetxt('IC_mixed_'+cross_sec+'_'+direction+'.txt', DataOut)

variable = "U"
direction = "y" #used only for vectorial variables (i.e. velocity). x, y or z
dir=2
cross_sec = "V"
confidence = 0.95 #level of confidence for the interval
std_dev_v = 0.1 #std deviation of the input
std_dev_mu = 0.001 #std deviation of the input



########## Sensitivity section  Velocity ##########

dx= np.loadtxt(test_case_v+variable+'_'+cross_sec+'.coupe', unpack=True, usecols=[0])
state_v= np.loadtxt(test_case_v+variable+'_'+cross_sec+'.coupe', unpack=True, usecols=[dir])
sens_v= np.loadtxt(test_case_v+variable+'_'+cross_sec+'_SENS.coupe', unpack=True, usecols=[dir])


state_mu= np.loadtxt(test_case_mu+variable+'_'+cross_sec+'.coupe', unpack=True, usecols=[dir])
sens_mu= np.loadtxt(test_case_mu+variable+'_'+cross_sec+'_SENS.coupe', unpack=True, usecols=[dir])



variance =np.abs(sens_v*std_dev_v + sens_mu*std_dev_mu)



upper_bound = state_v+1/np.sqrt(1-confidence)*variance
lower_bound = state_v-1/np.sqrt(1-confidence)*variance

DataOut = np.column_stack((dx,state_v, upper_bound, lower_bound))
np.savetxt('IC_mixed_'+cross_sec+'_'+direction+'.txt', DataOut)

variable = "U"
direction = "x" #used only for vectorial variables (i.e. velocity). x, y or z
dir=1
cross_sec = "V"
confidence = 0.95 #level of confidence for the interval
std_dev_v = 0.1 #std deviation of the input
std_dev_mu = 0.001 #std deviation of the input



########## Sensitivity section  Velocity ##########

dx= np.loadtxt(test_case_v+variable+'_'+cross_sec+'.coupe', unpack=True, usecols=[0])
state_v= np.loadtxt(test_case_v+variable+'_'+cross_sec+'.coupe', unpack=True, usecols=[dir])
sens_v= np.loadtxt(test_case_v+variable+'_'+cross_sec+'_SENS.coupe', unpack=True, usecols=[dir])


state_mu= np.loadtxt(test_case_mu+variable+'_'+cross_sec+'.coupe', unpack=True, usecols=[dir])
sens_mu= np.loadtxt(test_case_mu+variable+'_'+cross_sec+'_SENS.coupe', unpack=True, usecols=[dir])



variance =np.abs(sens_v*std_dev_v + sens_mu*std_dev_mu)



upper_bound = state_v+1/np.sqrt(1-confidence)*variance
lower_bound = state_v-1/np.sqrt(1-confidence)*variance

DataOut = np.column_stack((dx,state_v, upper_bound, lower_bound))
np.savetxt('IC_mixed_'+cross_sec+'_'+direction+'.txt', DataOut)
