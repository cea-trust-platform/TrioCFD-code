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
cross_sec = [ "X_C",  "Y_C"]


for j in range(len(cross_sec)):
	for i in range(len(params)):

		test_case = "Sensibility_Rayleigh_Bernard_Taylor/"+params[i]+"/Sensibility_Rayleigh_Bernard_"

		dx= np.loadtxt(test_case+variable_1+'_SENSIBILITE''_'+cross_sec[j]+'.coupe', unpack=True, usecols=[0])
		sens_vx= np.loadtxt(test_case+variable_1+'_SENSIBILITE''_'+cross_sec[j]+'.coupe', unpack=True, usecols=[1])
		sens_vy= np.loadtxt(test_case+variable_1+'_SENSIBILITE''_'+cross_sec[j]+'.coupe', unpack=True, usecols=[2])
		sens_t= np.loadtxt(test_case+variable_2+'_SENSIBILITE''_'+cross_sec[j]+'.coupe', unpack=True, usecols=[1])
		sens_p= np.loadtxt(test_case+variable_3+'_SENSIBILITE''_'+cross_sec[j]+'.coupe', unpack=True, usecols=[1])

		var_vx = abs(sens_vx)*sigma[i]
		var_vy = abs(sens_vy)*sigma[i]
		var_t = abs(sens_t)*sigma[i]
		var_p = abs(sens_p)*sigma[i]

		test_case = "Sensibility_Rayleigh_Bernard_PCM/"+params[i]+"/Sensibility_Rayleigh_Bernard_"

		dx= np.loadtxt(test_case+variable_1+'_SENSIBILITE''_'+cross_sec[j]+'.coupe', unpack=True, usecols=[0])
		pcm_sens_vx= np.loadtxt(test_case+variable_1+'_SENSIBILITE''_'+cross_sec[j]+'.coupe', unpack=True, usecols=[1])
		pcm_sens_vy= np.loadtxt(test_case+variable_1+'_SENSIBILITE''_'+cross_sec[j]+'.coupe', unpack=True, usecols=[2])
		pcm_sens_t= np.loadtxt(test_case+variable_2+'_SENSIBILITE''_'+cross_sec[j]+'.coupe', unpack=True, usecols=[1])
		pcm_sens_p= np.loadtxt(test_case+variable_3+'_SENSIBILITE''_'+cross_sec[j]+'.coupe', unpack=True, usecols=[1])

		pcm_var_vx = abs(pcm_sens_vx)*sigma[i]
		pcm_var_vy = abs(pcm_sens_vy)*sigma[i]
		pcm_var_t = abs(pcm_sens_t)*sigma[i]
		pcm_var_p = abs(pcm_sens_p)*sigma[i]

		#plot sensibility




		fig, ax = plt.subplots()
		ax.plot(dx, var_vx, dx, pcm_var_vx, 'r8', linewidth=2)
		ax.legend(['Taylor', 'PCM'])
		plt.title('Sonde '+cross_sec[j]+', param : '+params[i])
		plt.xlabel("x")
		plt.ylabel("Horizontal velocity component")
		plt.savefig(params[i]+'_Taylor_PCM_Var_Vx_'+cross_sec[j]+'.png')
		#plt.show()
		plt.close()

		fig, ax = plt.subplots()
		ax.plot(dx, var_vy, dx, pcm_var_vy, 'r8', linewidth=2)
		ax.legend(['Taylor', 'PCM'])
		plt.title('Sonde '+cross_sec[j]+', param : '+params[i])
		plt.xlabel("x")
		plt.ylabel("Vertical velocity component")
		plt.savefig(params[i]+'_Taylor_PCM_Var_Vy_'+cross_sec[j]+'.png')
		#plt.show()
		plt.close()

		fig, ax = plt.subplots()
		ax.plot(dx, var_t, dx, pcm_var_t, 'r8', linewidth=2)
		ax.legend(['Taylor', 'PCM'])
		plt.title('Sonde '+cross_sec[j]+', param : '+params[i])
		plt.xlabel("x")
		plt.ylabel("Temperature")
		plt.savefig(params[i]+'_Taylor_PCM_Var_T_'+cross_sec[j]+'.png')
		#plt.show()
		plt.close()

		fig, ax = plt.subplots()
		ax.plot(dx, var_p, dx, pcm_var_p, 'r8', linewidth=2)
		ax.legend(['Taylor', 'PCM'])
		plt.title('Sonde '+cross_sec[j]+', param : '+params[i])
		plt.xlabel("x")
		plt.ylabel("Pressure")
		plt.savefig(params[i]+'_Taylor_PCM_Var_P_'+cross_sec[j]+'.png')
		#plt.show()
		plt.close()
