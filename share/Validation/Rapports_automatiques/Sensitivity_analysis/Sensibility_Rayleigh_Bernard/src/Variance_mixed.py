# Sensitivity
import numpy as np
import matplotlib.pyplot as plt

params = ["Tbottom", "Ttop", "Beta", "Mu", "Lambda_prime"]
mean_params = [363, 313, 3e-3, 1.54e-5, 0.21847e-4]
sigma = mean_params
for i in range(len(params)):
	sigma[i] = mean_params[i]/100.

variable_1 = "VITESSE"
variable_2 = "TEMPERATURE"
variable_3 = "PRESSION"

confidence = 0.95 #level of confidence for the interval

var_vx  = 0
var_vy = 0
var_t = 0 
var_p = 0
pcm_var_vx=0
pcm_var_vy=0
pcm_var_t=0
pcm_var_p=0


cross_sec = "X_C"

for i in range(len(params)):

		test_case = "Sensibility_Rayleigh_Bernard_Taylor/"+params[i]+"/Sensibility_Rayleigh_Bernard_"

		dx= np.loadtxt(test_case+variable_1+'_SENSIBILITE''_'+cross_sec+'.coupe', unpack=True, usecols=[0])
		sens_vx= np.loadtxt(test_case+variable_1+'_SENSIBILITE''_'+cross_sec+'.coupe', unpack=True, usecols=[1])
		sens_vy= np.loadtxt(test_case+variable_1+'_SENSIBILITE''_'+cross_sec+'.coupe', unpack=True, usecols=[2])
		sens_t= np.loadtxt(test_case+variable_2+'_SENSIBILITE''_'+cross_sec+'.coupe', unpack=True, usecols=[1])
		sens_p= np.loadtxt(test_case+variable_3+'_SENSIBILITE''_'+cross_sec+'.coupe', unpack=True, usecols=[1])

		var_vx += abs(sens_vx)*sigma[i]
		var_vy += abs(sens_vy)*sigma[i]
		var_t += abs(sens_t)*sigma[i]
		var_p += abs(sens_p)*sigma[i]

		test_case = "Sensibility_Rayleigh_Bernard_PCM/"+params[i]+"/Sensibility_Rayleigh_Bernard_"

		dx= np.loadtxt(test_case+variable_1+'_SENSIBILITE''_'+cross_sec+'.coupe', unpack=True, usecols=[0])
		pcm_sens_vx= np.loadtxt(test_case+variable_1+'_SENSIBILITE''_'+cross_sec+'.coupe', unpack=True, usecols=[1])
		pcm_sens_vy= np.loadtxt(test_case+variable_1+'_SENSIBILITE''_'+cross_sec+'.coupe', unpack=True, usecols=[2])
		pcm_sens_t= np.loadtxt(test_case+variable_2+'_SENSIBILITE''_'+cross_sec+'.coupe', unpack=True, usecols=[1])
		pcm_sens_p= np.loadtxt(test_case+variable_3+'_SENSIBILITE''_'+cross_sec+'.coupe', unpack=True, usecols=[1])

		pcm_var_vx += abs(pcm_sens_vx)*sigma[i]
		pcm_var_vy += abs(pcm_sens_vy)*sigma[i]
		pcm_var_t += abs(pcm_sens_t)*sigma[i]
		pcm_var_p += abs(pcm_sens_p)*sigma[i]


test_case = "Sensibility_Rayleigh_Bernard_PCM/"+params[0]+"/Sensibility_Rayleigh_Bernard_"
state_vx= np.loadtxt(test_case+variable_1+'_ETAT'+'_'+cross_sec+'.coupe', unpack=True, usecols=[1])
upper_bound = state_vx+1/np.sqrt(1-confidence)*pcm_var_vx
lower_bound = state_vx-1/np.sqrt(1-confidence)*pcm_var_vx
DataOut = np.column_stack((dx,state_vx, upper_bound, lower_bound))
np.savetxt('IC_mixed_'+cross_sec+'_'+'Vx'+'.txt', DataOut)

state_vy= np.loadtxt(test_case+variable_1+'_ETAT'+'_'+cross_sec+'.coupe', unpack=True, usecols=[2])
upper_bound = state_vy+1/np.sqrt(1-confidence)*pcm_var_vy
lower_bound = state_vy-1/np.sqrt(1-confidence)*pcm_var_vy
DataOut = np.column_stack((dx,state_vy, upper_bound, lower_bound))
np.savetxt('IC_mixed_'+cross_sec+'_'+'Vy'+'.txt', DataOut)            
                

state_p= np.loadtxt(test_case+variable_2+'_ETAT'+'_'+cross_sec+'.coupe', unpack=True, usecols=[1])
upper_bound = state_p+1/np.sqrt(1-confidence)*pcm_var_p
lower_bound = state_p-1/np.sqrt(1-confidence)*pcm_var_p
DataOut = np.column_stack((dx,state_p, upper_bound, lower_bound))
np.savetxt('IC_mixed_'+cross_sec+'_'+'P'+'.txt', DataOut) 

state_t= np.loadtxt(test_case+variable_3+'_ETAT'+'_'+cross_sec+'.coupe', unpack=True, usecols=[1])
upper_bound = state_t+1/np.sqrt(1-confidence)*pcm_var_t
lower_bound = state_t-1/np.sqrt(1-confidence)*pcm_var_t
DataOut = np.column_stack((dx,state_t, upper_bound, lower_bound))
np.savetxt('IC_mixed_'+cross_sec+'_'+'T'+'.txt', DataOut) 
		

for i in range(len(params)):

		#plot sensibility
		test_case = "Sensibility_Rayleigh_Bernard_Taylor/"+params[i]+"/Sensibility_Rayleigh_Bernard_"
		fig, ax = plt.subplots()
		ax.plot(dx, var_vx, dx, pcm_var_vx, 'r8', linewidth=2)
		ax.legend(['Taylor', 'PCM'])
		plt.title('Horizontal velocity component standard deviation')
		plt.xlabel("x")
		plt.ylabel("Horizontal velocity component standard deviation")
		plt.savefig('Taylor_PCM_Var_Vx_'+cross_sec+'.png')
		#plt.show()
		plt.close()

		fig, ax = plt.subplots()
		ax.plot(dx, var_vy, dx, pcm_var_vy, 'r8', linewidth=2)
		ax.legend(['Taylor', 'PCM'])
		plt.title('Vertical velocity component standard deviation')
		plt.xlabel("x")
		plt.ylabel("Vertical velocity component standard deviation")
		plt.savefig('Taylor_PCM_Var_Vy_'+cross_sec+'.png')
		#plt.show()
		plt.close()

		fig, ax = plt.subplots()
		ax.plot(dx, var_t, dx, pcm_var_t, 'r8', linewidth=2)
		ax.legend(['Taylor', 'PCM'])
		plt.title('Temperature standard deviation')
		plt.xlabel("x")
		plt.ylabel("Temperature standard deviation")
		plt.savefig('Taylor_PCM_Var_T_'+cross_sec+'.png')
		#plt.show()
		plt.close()

		fig, ax = plt.subplots()
		ax.plot(dx, var_p, dx, pcm_var_p, 'r8', linewidth=2)
		ax.legend(['Taylor', 'PCM'])
		plt.title('Pressure standard deviation')
		plt.xlabel("x")
		plt.ylabel("Pressure standard deviation")
		plt.savefig('Taylor_PCM_Var_P_'+cross_sec+'.png')
		#plt.show()
		plt.close()



var_vx  = 0
var_vy = 0
var_t = 0 
var_p = 0
pcm_var_vx=0
pcm_var_vy=0
pcm_var_t=0
pcm_var_p=0
		
cross_sec = "Y_C"
for i in range(len(params)):

		test_case = "Sensibility_Rayleigh_Bernard_Taylor/"+params[i]+"/Sensibility_Rayleigh_Bernard_"

		dx= np.loadtxt(test_case+variable_1+'_SENSIBILITE''_'+cross_sec+'.coupe', unpack=True, usecols=[0])
		sens_vx= np.loadtxt(test_case+variable_1+'_SENSIBILITE''_'+cross_sec+'.coupe', unpack=True, usecols=[1])
		sens_vy= np.loadtxt(test_case+variable_1+'_SENSIBILITE''_'+cross_sec+'.coupe', unpack=True, usecols=[2])
		sens_t= np.loadtxt(test_case+variable_2+'_SENSIBILITE''_'+cross_sec+'.coupe', unpack=True, usecols=[1])
		sens_p= np.loadtxt(test_case+variable_3+'_SENSIBILITE''_'+cross_sec+'.coupe', unpack=True, usecols=[1])

		var_vx += abs(sens_vx)*sigma[i]
		var_vy += abs(sens_vy)*sigma[i]
		var_t += abs(sens_t)*sigma[i]
		var_p += abs(sens_p)*sigma[i]

		test_case = "Sensibility_Rayleigh_Bernard_PCM/"+params[i]+"/Sensibility_Rayleigh_Bernard_"

		dx= np.loadtxt(test_case+variable_1+'_SENSIBILITE''_'+cross_sec+'.coupe', unpack=True, usecols=[0])
		pcm_sens_vx= np.loadtxt(test_case+variable_1+'_SENSIBILITE''_'+cross_sec+'.coupe', unpack=True, usecols=[1])
		pcm_sens_vy= np.loadtxt(test_case+variable_1+'_SENSIBILITE''_'+cross_sec+'.coupe', unpack=True, usecols=[2])
		pcm_sens_t= np.loadtxt(test_case+variable_2+'_SENSIBILITE''_'+cross_sec+'.coupe', unpack=True, usecols=[1])
		pcm_sens_p= np.loadtxt(test_case+variable_3+'_SENSIBILITE''_'+cross_sec+'.coupe', unpack=True, usecols=[1])

		pcm_var_vx += abs(pcm_sens_vx)*sigma[i]
		pcm_var_vy += abs(pcm_sens_vy)*sigma[i]
		pcm_var_t += abs(pcm_sens_t)*sigma[i]
		pcm_var_p += abs(pcm_sens_p)*sigma[i]




test_case = "Sensibility_Rayleigh_Bernard_PCM/"+params[0]+"/Sensibility_Rayleigh_Bernard_"
state_vx= np.loadtxt(test_case+variable_1+'_ETAT'+'_'+cross_sec+'.coupe', unpack=True, usecols=[1])
upper_bound = state_vx+1/np.sqrt(1-confidence)*pcm_var_vx
lower_bound = state_vx-1/np.sqrt(1-confidence)*pcm_var_vx
DataOut = np.column_stack((dx,state_vx, upper_bound, lower_bound))
np.savetxt('IC_mixed_'+cross_sec+'_'+'Vx'+'.txt', DataOut)

state_vy= np.loadtxt(test_case+variable_1+'_ETAT'+'_'+cross_sec+'.coupe', unpack=True, usecols=[2])
upper_bound = state_vy+1/np.sqrt(1-confidence)*pcm_var_vy
lower_bound = state_vy-1/np.sqrt(1-confidence)*pcm_var_vy
DataOut = np.column_stack((dx,state_vy, upper_bound, lower_bound))
np.savetxt('IC_mixed_'+cross_sec+'_'+'Vy'+'.txt', DataOut)            
                

state_p= np.loadtxt(test_case+variable_2+'_ETAT'+'_'+cross_sec+'.coupe', unpack=True, usecols=[1])
upper_bound = state_p+1/np.sqrt(1-confidence)*pcm_var_p
lower_bound = state_p-1/np.sqrt(1-confidence)*pcm_var_p
DataOut = np.column_stack((dx,state_p, upper_bound, lower_bound))
np.savetxt('IC_mixed_'+cross_sec+'_'+'P'+'.txt', DataOut) 

state_t= np.loadtxt(test_case+variable_3+'_ETAT'+'_'+cross_sec+'.coupe', unpack=True, usecols=[1])
upper_bound = state_t+1/np.sqrt(1-confidence)*pcm_var_t
lower_bound = state_t-1/np.sqrt(1-confidence)*pcm_var_t
DataOut = np.column_stack((dx,state_t, upper_bound, lower_bound))
np.savetxt('IC_mixed_'+cross_sec+'_'+'T'+'.txt', DataOut) 

for i in range(len(params)):

		#plot sensibility
		test_case = "Sensibility_Rayleigh_Bernard_Taylor/"+params[i]+"/Sensibility_Rayleigh_Bernard_"
		fig, ax = plt.subplots()
		ax.plot(dx, var_vx, dx, pcm_var_vx, 'r8', linewidth=2)
		ax.legend(['Taylor', 'PCM'])
		plt.title('Horizontal velocity component standard deviation')
		plt.xlabel("y")
		plt.ylabel("Horizontal velocity component standard deviation")
		plt.savefig('Taylor_PCM_Var_Vx_'+cross_sec+'.png')
		#plt.show()
		plt.close()

		fig, ax = plt.subplots()
		ax.plot(dx, var_vy, dx, pcm_var_vy, 'r8', linewidth=2)
		ax.legend(['Taylor', 'PCM'])
		plt.title('Vertical velocity component standard deviation')
		plt.xlabel("y")
		plt.ylabel("Vertical velocity component standard deviation")
		plt.savefig('Taylor_PCM_Var_Vy_'+cross_sec+'.png')
		#plt.show()
		plt.close()

		fig, ax = plt.subplots()
		ax.plot(dx, var_t, dx, pcm_var_t, 'r8', linewidth=2)
		ax.legend(['Taylor', 'PCM'])
		plt.title('Temperature standard deviation')
		plt.xlabel("y")
		plt.ylabel("Temperature standard deviation")
		plt.savefig('Taylor_PCM_Var_T_'+cross_sec+'.png')
		#plt.show()
		plt.close()

		fig, ax = plt.subplots()
		ax.plot(dx, var_p, dx, pcm_var_p, 'r8', linewidth=2)
		ax.legend(['Taylor', 'PCM'])
		plt.title('Pressure standard deviation')
		plt.xlabel("y")
		plt.ylabel("Pressure standard deviation")
		plt.savefig('Taylor_PCM_Var_P_'+cross_sec+'.png')
		#plt.show()
		plt.close()		
