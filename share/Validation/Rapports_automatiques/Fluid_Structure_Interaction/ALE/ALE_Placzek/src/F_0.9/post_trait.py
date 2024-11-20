import numpy as np
import matplotlib.pyplot as plt
from scipy import signal


##Import des données

F_p=np.loadtxt('ALE_Placzek_pb_Force_pression.out') #Force de pression

Fv_p = np.loadtxt('ALE_Placzek_pb_Contrainte_visqueuse.out') #Contrainte visqueuse

placzek = np.loadtxt('../ReferenceSolution/Placzek_Cl_alpha_0.9.dat')

##Paramètres de l'étude

Reynolds=100       #Reynolds number 
D=0.001             #cylinder diameter
nu=1.e-6            #cinematic viscosity
rho=1000            #volumic mass

V=nu*Reynolds/D     #inlet velocity
f_v = 16.87        #Strouhal frequency

## Paramètres à changer
KC=0.25            #Kalegan Karpenter number (A/D)
f0=f_v*0.9           #imposed frequency on the cylinder
nb_p_i=25        #number of period to ignore

omega=2*np.pi*f0      #pulsation
A=KC*D             #displacment amplitude
T=1/f0             #oscillation period 
ti=T*nb_p_i        #start time to compute coefficient 


# How many period are available ?
nb_p=1
for i in range(len(F_p[:,0])):
	t=F_p[i,0]
	if t > nb_p*T:
		nb_p=nb_p+1
		l=i            #the indice of the time corresponding to the period-ending is kept
tf=nb_p*T          #final time to compute coefficient 

#data selection in the time interval asked
t=0
i=0
while t < ti:
	t=F_p[i,0]
	i=i+1




##sélection des données
t= F_p[i:l,0]       #time
Fx=F_p[i:l,1]      #pressure force x direction
Fv_x=Fv_p[i:l,1]    #viscous force x direction
Fy=F_p[i:l,2]       #pressure force y direction
Fv_y=Fv_p[i:l,2]    #viscous force y direction

alpha=KC*np.sin(omega*t) #+np.full_like(t,np.pi)# #vertical displacement
v_ale=KC*omega*np.cos(omega*t) #vertical speed

Fx_t=Fx+Fv_x                #total force x direction
Fy_t=Fy+Fv_y                #total force y direction



##-------Drag------- 
Drag_mean=np.mean(Fx_t) 
Drag_std=np.std(Fx_t) 

Cd=Fx_t/(0.5*rho*V*V*D) 
Cd_mean=np.mean(Cd) 
Cd_std=np.std(Cd) 

# print(f'Cd moyen = {Cd_mean}')

 ##-------Lift------- 
Lift_mean=np.mean(Fy_t) 
Lift_std=np.std(Fy_t) 

Cl_brut=(Fy+Fv_y)/(0.5*rho*V*V*D) 

# Fenêtre de lissage (nombre de points à utiliser pour la moyenne mobile)
window = 12

# Calcul de la moyenne mobile
Cl = np.convolve(Cl_brut, np.ones(window)/window, mode='valid')

Cl_mean=np.mean(Cl) 
Cl_std=np.std(Cl) 
Cl_peak=np.max(Cl)

#print(f'Cl max = {Cl_peak}')

DataOut = np.column_stack((round(Cd_mean,3),round(Cl_peak,3)))
np.savetxt('Cd_Cl.txt', DataOut)


# Color settings (CEA colors)

# Primary colours

cea_red = (229/255,0,25/255)
cea_black = (0,0,0)
cea_white = (1,1,1)
cea_darkblue = (62/255,74/255,131/255)
cea_lightblue = (126/255,156/255,187/255)
cea_darkgrey = (38/255,38/255,38/255)
cea_yellow = (1,205/255,49/255)

# Additional ones
macaron = (218/255,131/255,123/255)
archipel = (0,147/255,157/255)
glycine = (167/255,37/255,135/255)
opera = (189/255,152/255,122/255)

# Plot Settings
A = 5.5 # Want figures to be A6
plt.rc('figure', figsize=[46.82 * .5**(.5 * A), 33.11 * .5**(.5 * A)])
plt.rc('text', usetex=True)

plt.rc('font', family='serif', size=18) #font of axislabels and title, write r before label string


#Cl en fonction de t
plt.plot(t[window-1:]/T, Cl,color =cea_darkblue)
plt.xlabel(r'$t/T_0$')
plt.ylabel(r'$C_L$')
plt.ylim(-1.,1.)
plt.xlim(nb_p_i,19)
plt.grid(True)
plt.tight_layout()
plt.savefig('Cl.png', bbox_inches='tight')
plt.close()

# Calculer la transformée de Fourier discrète (DFT) du signal
fft_result = np.fft.fft(Cl)

# Calculer la fréquence associée à chaque composante de la DFT
frequencies = np.fft.fftfreq(len(Cl), np.mean(np.diff(t)))

# Calculer la densité de puissance spectrale (PSD)
psd = np.abs(fft_result)**2 / len(Cl)


plt.plot(frequencies/f0, psd/max(psd),color = cea_darkblue)
plt.xlabel(r'$F^*=f/f_0$')
plt.xticks(np.arange(0, 5, 0.5))
plt.ylim(0,1.1)
plt.xlim(0,5)
plt.ylabel(r'$PSD_n$')
plt.grid(True)
plt.tight_layout()
plt.savefig('PSD.png', bbox_inches='tight')
plt.close()

#Cl en fonction de alpha

plt.plot(alpha[window-1:] ,Cl, label=r'Étude actuelle', color=cea_darkblue, ls='-', zorder=1)
plt.scatter(placzek[:,0], placzek[:,1], label = r'Article Placzek',marker='+', color = cea_red, s=30, zorder=2)

plt.xlabel(r'$\alpha$')
plt.ylabel(r'$C_L$')
#plt.ylim(-2,2)
plt.grid(True)
plt.xticks(np.arange(-0.25, 0.3, 0.1))
plt.yticks(np.arange(-0.15, 0.2, 0.05))
plt.ylim(-0.16,0.16)
plt.xlim(-0.3,0.3)
plt.legend(loc='upper center', bbox_to_anchor=(0.45, -0.15), ncol=2, fontsize = 18)
plt.tight_layout()
plt.savefig("Cl_alpha.png", bbox_inches='tight')
plt.close()
