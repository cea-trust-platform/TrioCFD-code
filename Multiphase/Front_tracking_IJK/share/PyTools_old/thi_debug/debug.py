import glob
import numpy as np
import matplotlib.pyplot as plt

# emploi : python degug.py
# /!\ Se placer dans un repertoire de RUN

length    = [] 
d_time    = []
time__    = []
velocity_ = []
runs      = []

for i in range(9):
    run="RUN0"+str(i)
    try:
        print(glob.glob("../"+run+"/debug_4.txt"))
        length.append(np.loadtxt("../"+run+"/debug_4.txt",usecols=4))
        d_time.append(np.loadtxt("../"+run+"/debug_3.txt",usecols=4))
        time__.append(np.loadtxt("../"+run+"/debug_1.txt",usecols=3))
        velocity_.append(length[-1]/time__[-1])
        runs.append(run)
    except:
        print(run+" not available for debugs.")
        
tf = 0

plt.figure(1)
# ~ Preparation ~~~~~~~~~~~~
plt.xlabel("temps physique")
plt.ylabel("vitesse d'advection du champ THI")
plt.grid(True,'both')
plt.tight_layout()
# ~ Traces ~~~~~~~~~~~~~~~~~
for i in range(len(runs)):
    plt.plot(time__[i],velocity_[i]   ,label=run)
    tf += time__[i][-1]              
plt.legend(framealpha=0)
plt.savefig("vitesse_advection_thi.pdf")

plt.figure(2)
# ~ Preparation ~~~~~~~~~~~~~~~~
plt.xlabel("temps physique")
plt.ylabel("deplacement articiciel du champ THI")
plt.grid(True,'both')
plt.tight_layout()
# ~ Traces ~~~~~~~~~~~~~~~~~~~~~
for i in range(len(runs)):
    plt.plot(time__[i],length[i]   ,label=run)
    tf += time__[i][-1]              
plt.legend(framealpha=0)
plt.legend(framealpha=0)
plt.savefig("longueur_advection_thi.pdf")

plt.show()
