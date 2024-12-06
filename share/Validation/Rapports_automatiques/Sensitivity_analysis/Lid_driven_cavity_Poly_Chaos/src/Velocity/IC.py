# Sensitivity
import numpy as np
import matplotlib.pyplot as plt


variable = "U"

confidence = 0.95 #level of confidence for the interval
std_dev = 0.1 #std deviation of the input





cross_sec = "H"

test_case = "../Taylor/Velocity/Lid_v_"

dx= np.loadtxt(test_case+variable+'_'+cross_sec+'.coupe', unpack=True, usecols=[0])
state_x= np.loadtxt(test_case+variable+'_'+cross_sec+'.coupe', unpack=True, usecols=[1])
sens_x= np.loadtxt(test_case+variable+'_'+cross_sec+'_SENS.coupe', unpack=True, usecols=[1])
state_y= np.loadtxt(test_case+variable+'_'+cross_sec+'.coupe', unpack=True, usecols=[2])
sens_y= np.loadtxt(test_case+variable+'_'+cross_sec+'_SENS.coupe', unpack=True, usecols=[2])


variance_x = np.abs(sens_x)*std_dev
variance_y = np.abs(sens_y)*std_dev




test_case = "Lid_driven_cavity_Poly_Chaos_"

state_x_pcm= np.loadtxt(test_case+variable+'_'+cross_sec+'.coupe', unpack=True, usecols=[1])
sens_x_pcm= np.loadtxt(test_case+variable+'_'+cross_sec+'_SENS.coupe', unpack=True, usecols=[1])
state_y_pcm= np.loadtxt(test_case+variable+'_'+cross_sec+'.coupe', unpack=True, usecols=[2])
sens_y_pcm= np.loadtxt(test_case+variable+'_'+cross_sec+'_SENS.coupe', unpack=True, usecols=[2])


variance_x_pcm = np.abs(sens_x_pcm)*std_dev
variance_y_pcm = np.abs(sens_y_pcm)*std_dev
state_x_pcm = state_x_pcm + sens_x_pcm*std_dev
state_y_pcm = state_y_pcm + sens_y_pcm*std_dev


#plot standard deviation x

ref= np.loadtxt('../Ref_sens/velocity_H_x.txt', unpack=True, usecols=[1])
t_ref= np.loadtxt('../Ref_sens/velocity_H_x.txt', unpack=True, usecols=[0])



fig, ax = plt.subplots()
ax.plot(t_ref, ref, dx, variance_x, dx, variance_x_pcm, 'r8', linewidth=2)
ax.legend(['El-Beltagy and Wafa', 'Taylor', 'PCM'])
plt.xlabel("x")
plt.ylabel("Horizontal velocity component standard deviation")
plt.savefig('Taylor_PCM_Variance_Velocity_H_x.png')
#plt.show()



#plot average

ref= np.loadtxt('../Ref_etat/Ref5_y_velocity_horizontal.dat', unpack=True, usecols=[1])
t_ref= np.loadtxt('../Ref_etat/Ref5_y_velocity_horizontal.dat', unpack=True, usecols=[0])



fig, ax = plt.subplots()
ax.plot(t_ref, ref, dx, state_y, dx, state_y_pcm, 'r8', linewidth=2)
ax.legend(['Marchi et al.', 'Taylor', 'PCM'])
plt.xlabel("x")
plt.ylabel("Horizontal velocity component average")
plt.savefig('Taylor_PCM_Mean_Velocity_H_y.png')
#plt.show()



cross_sec = "V"

test_case = "../Taylor/Velocity/Lid_v_"


dx= np.loadtxt(test_case+variable+'_'+cross_sec+'.coupe', unpack=True, usecols=[0])
state_x= np.loadtxt(test_case+variable+'_'+cross_sec+'.coupe', unpack=True, usecols=[1])
sens_x= np.loadtxt(test_case+variable+'_'+cross_sec+'_SENS.coupe', unpack=True, usecols=[1])
state_y= np.loadtxt(test_case+variable+'_'+cross_sec+'.coupe', unpack=True, usecols=[2])
sens_y= np.loadtxt(test_case+variable+'_'+cross_sec+'_SENS.coupe', unpack=True, usecols=[2])


variance_x = np.abs(sens_x)*std_dev
variance_y = np.abs(sens_y)*std_dev


test_case = "Lid_driven_cavity_Poly_Chaos_"

state_x_pcm= np.loadtxt(test_case+variable+'_'+cross_sec+'.coupe', unpack=True, usecols=[1])
sens_x_pcm= np.loadtxt(test_case+variable+'_'+cross_sec+'_SENS.coupe', unpack=True, usecols=[1])
state_y_pcm= np.loadtxt(test_case+variable+'_'+cross_sec+'.coupe', unpack=True, usecols=[2])
sens_y_pcm= np.loadtxt(test_case+variable+'_'+cross_sec+'_SENS.coupe', unpack=True, usecols=[2])

variance_x_pcm = np.abs(sens_x_pcm)*std_dev
variance_y_pcm = np.abs(sens_y_pcm)*std_dev
state_x_pcm = state_x_pcm + sens_x_pcm*std_dev
state_y_pcm = state_y_pcm + sens_y_pcm*std_dev
#plot average

ref= np.loadtxt('../Ref_etat/Ref5_x_velocity_vertical.dat', unpack=True, usecols=[1])
t_ref= np.loadtxt('../Ref_etat/Ref5_x_velocity_vertical.dat', unpack=True, usecols=[0])



fig, ax = plt.subplots()
ax.plot(t_ref, ref, dx, state_x, dx, state_x_pcm, 'r8', linewidth=2)
ax.legend([ 'Marchi et al.', 'Taylor', 'PCM'])
plt.xlabel("y")
plt.ylabel("Vertical velocity component average")
plt.savefig('Taylor_PCM_Mean_Velocity_V_x.png')
#plt.show()

