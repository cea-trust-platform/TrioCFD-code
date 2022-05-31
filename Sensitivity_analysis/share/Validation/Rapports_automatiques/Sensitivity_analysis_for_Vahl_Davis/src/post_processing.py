test_case = "VAHL_DAVIS"
params = ["Lambda", "Tg", "Td", "Beta", "Cp"]
norm_params = [0.0262, 283.15, 273.15, 3.6-3, 1006]
state = "State"
variable = "TEMPERATURE"
direction = 2 #only for vectorial variables (i.e. velocity). 1 for x, 2 for y
cross_sec = "X_C"

import numpy as np
import math
import matplotlib.pyplot as plt
x = []
sens = []
for i in range(len(params)):
    sens.append([])
    with open(params[i]+'/'+test_case+'_'+variable+'_SENSIBILITE_'+cross_sec+'.coupe', 'r') as f:
        for line in f:
            l = line.split(' ')
            if l[0] != '\n':
                if(i == 0):
                    x.append(float(l[0]))
                if variable != "VITESSE":
                    sens[i].append((float(l[1])*norm_params[i])**2)
                else:
                    sens[i].append((float(l[direction])*norm_params[i])**2)

sens_tot = np.zeros(np.shape(x))
sens = np.array(sens)
for i in range(len(params)):
    sens_tot = sens_tot + sens[i]
    
for i in range(len(params)):
   plt.plot(x, sens[i]/sens_tot, label = params[i])
plt.legend()
plt.savefig(variable+'_'+cross_sec+'.png')
plt.close()
cross_sec = "X_H"
x = []
sens = []
for i in range(len(params)):
    sens.append([])
    with open(params[i]+'/'+test_case+'_'+variable+'_SENSIBILITE_'+cross_sec+'.coupe', 'r') as f:
        for line in f:
            l = line.split(' ')
            if l[0] != '\n':
                if(i == 0):
                    x.append(float(l[0]))
                if variable != "VITESSE":
                    sens[i].append((float(l[1])*norm_params[i])**2)
                else:
                    sens[i].append((float(l[direction])*norm_params[i])**2)

sens_tot = np.zeros(np.shape(x))
sens = np.array(sens)
for i in range(len(params)):
    sens_tot = sens_tot + sens[i]
    
for i in range(len(params)):
   plt.plot(x, sens[i]/sens_tot, label = params[i])
plt.legend()
plt.savefig(variable+'_'+cross_sec+'.png')
plt.close()

variable = "VITESSE"
direction = 2 #only for vectorial variables (i.e. velocity). 1 for x, 2 for y
cross_sec = "X_C"
x = []
sens = []
for i in range(len(params)):
    sens.append([])
    with open(params[i]+'/'+test_case+'_'+variable+'_SENSIBILITE_'+cross_sec+'.coupe', 'r') as f:
        for line in f:
            l = line.split(' ')
            if l[0] != '\n':
                if(i == 0):
                    x.append(float(l[0]))
                if variable != "VITESSE":
                    sens[i].append((float(l[1])*norm_params[i])**2)
                else:
                    sens[i].append((float(l[direction])*norm_params[i])**2)

sens_tot = np.zeros(np.shape(x))
sens = np.array(sens)
for i in range(len(params)):
    sens_tot = sens_tot + sens[i]
    
for i in range(len(params)):
   plt.plot(x, sens[i]/sens_tot, label = params[i])
plt.legend()
plt.savefig(variable+'_'+cross_sec+'.png')
plt.close()

variable = "VITESSE"
direction = 1 #only for vectorial variables (i.e. velocity). 1 for x, 2 for y
cross_sec = "X_H"
x = []
sens = []
for i in range(len(params)):
    sens.append([])
    with open(params[i]+'/'+test_case+'_'+variable+'_SENSIBILITE_'+cross_sec+'.coupe', 'r') as f:
        for line in f:
            l = line.split(' ')
            if l[0] != '\n':
                if(i == 0):
                    x.append(float(l[0]))
                if variable != "VITESSE":
                    sens[i].append((float(l[1])*norm_params[i])**2)
                else:
                    sens[i].append((float(l[direction])*norm_params[i])**2)

sens_tot = np.zeros(np.shape(x))
sens = np.array(sens)
for i in range(len(params)):
    sens_tot = sens_tot + sens[i]
    
for i in range(len(params)):
   plt.plot(x, sens[i]/sens_tot, label = params[i])
plt.legend()
plt.savefig(variable+'_'+cross_sec+'.png')
plt.close()

