# Sensitivity
test_case_v = "Taylor/Velocity/Lid_v_"
test_case_mu = "Taylor/Mu/Lid_mu_"
variable = "U"
direction = "x" #used only for vectorial variables (i.e. velocity). x, y or z
cross_sec = "H"
confidence = 0.95 #level of confidence for the interval
std_dev_v = 0.1 #std deviation of the input
std_dev_mu = 0.001 #std deviation of the input


import numpy as np
import matplotlib.pyplot as plt

########## Sensitivity section  Velocity ##########

abscissa= np.loadtxt(test_case_v+variable+'_'+cross_sec+'.coupe', unpack=True, usecols=[0])
state_v= np.loadtxt(test_case_v+variable+'_'+cross_sec+'.coupe', unpack=True, usecols=[1])
sens_v= np.loadtxt(test_case_v+variable+'_'+cross_sec+'_SENS.coupe', unpack=True, usecols=[1])


state_mu= np.loadtxt(test_case_mu+variable+'_'+cross_sec+'.coupe', unpack=True, usecols=[1])
sens_mu= np.loadtxt(test_case_mu+variable+'_'+cross_sec+'_SENS.coupe', unpack=True, usecols=[1])



variance =np.abs(sens_v*std_dev_v + sens_mu*std_dev_mu)


test_case_v = "Velocity/Lid_driven_cavity_Poly_Chaos_"
test_case_mu = "Mu/Lid_driven_cavity_Poly_Chaos_"


state_v_pcm= np.loadtxt(test_case_v+variable+'_'+cross_sec+'.coupe', unpack=True, usecols=[1])
sens_v_pcm= np.loadtxt(test_case_v+variable+'_'+cross_sec+'_SENS.coupe', unpack=True, usecols=[1])


state_mu_pcm= np.loadtxt(test_case_mu+variable+'_'+cross_sec+'.coupe', unpack=True, usecols=[1])
sens_mu_pcm= np.loadtxt(test_case_mu+variable+'_'+cross_sec+'_SENS.coupe', unpack=True, usecols=[1])



variance_pcm =np.abs(sens_v_pcm*std_dev_v + sens_mu_pcm*std_dev_mu)

ref= np.loadtxt('Ref_sens/mixed.txt', unpack=True, usecols=[1])
t_ref= np.loadtxt('Ref_sens/mixed.txt', unpack=True, usecols=[0])

fig, ax = plt.subplots()

ax.plot(t_ref, ref, abscissa, variance_pcm, abscissa, variance,'r8', linewidth=2)
ax.legend([ 'El-Beltagy and Wafa', 'Taylor', 'PCM'])
plt.xlabel("x")
plt.ylabel("Horizontal velocity standard deviation")
plt.savefig('Taylor_PCM_Variance_'+variable+'_'+cross_sec+direction+'mixed_v_plus_mu.png')
#plt.show()



######
# Sensitivity
test_case_v = "Velocity/Lid_driven_cavity_Poly_Chaos_"
test_case_mu = "Mu/Lid_driven_cavity_Poly_Chaos_"
variable = "U"
direction = "x" #used only for vectorial variables (i.e. velocity). x, y or z
cross_sec = "H"
confidence = 0.95 #level of confidence for the interval
std_dev_v = 0.1 #std deviation of the input
std_dev_mu = 0.001 #std deviation of the input


import numpy as np
import matplotlib.pyplot as plt

########## Sensitivity section  Velocity ##########

abscissa= np.loadtxt(test_case_v+variable+'_'+cross_sec+'.coupe', unpack=True, usecols=[0])
state_v= np.loadtxt(test_case_v+variable+'_'+cross_sec+'.coupe', unpack=True, usecols=[1])
sens_v= np.loadtxt(test_case_v+variable+'_'+cross_sec+'_SENS.coupe', unpack=True, usecols=[1])


state_mu= np.loadtxt(test_case_mu+variable+'_'+cross_sec+'.coupe', unpack=True, usecols=[1])
sens_mu= np.loadtxt(test_case_mu+variable+'_'+cross_sec+'_SENS.coupe', unpack=True, usecols=[1])



variance =np.abs(sens_v*std_dev_v + sens_mu*std_dev_mu)


ref= np.loadtxt('Ref_sens/mixed.txt', unpack=True, usecols=[1])
t_ref= np.loadtxt('Ref_sens/mixed.txt', unpack=True, usecols=[0])


ref_v= np.loadtxt('Ref_sens/velocity_H_x.txt', unpack=True, usecols=[1])
t_ref_v= np.loadtxt('Ref_sens/velocity_H_x.txt', unpack=True, usecols=[0])

ref_mu= np.loadtxt('Ref_sens/ux_H.txt', unpack=True, usecols=[1])
t_ref_mu= np.loadtxt('Ref_sens/ux_H.txt', unpack=True, usecols=[0])





variance_v =np.abs(sens_v*std_dev_v )
variance_mu =np.abs(sens_mu*std_dev_mu)

fig, ax = plt.subplots()
ax.plot(abscissa, variance_v, 'ro', t_ref_v, ref_v, 'r', abscissa, variance_mu,'bo', t_ref_mu, ref_mu, 'b', abscissa, variance, 'go', t_ref, ref, 'g', linewidth=2, fillstyle='none')
ax.legend([ 'Velocity 10%', 'El-Beltagy and Wafa', 'Viscosity 10 %', 'El-Beltagy and Wafa', 'Mixed', 'El-Beltagy and Wafa'])
#plt.ylim([0, 0.05])
plt.xlabel("x")
plt.ylabel("Horizontal velocity standard deviation")
plt.savefig('Variance_'+variable+'_'+cross_sec+'_'+direction+'.png')
#plt.show()




