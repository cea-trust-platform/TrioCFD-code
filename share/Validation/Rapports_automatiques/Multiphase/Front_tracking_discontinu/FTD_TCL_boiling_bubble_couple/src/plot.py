import numpy as np
import matplotlib.pyplot as plt
mat= np.loadtxt("SEQ_Steady_State_T_SEG.son").T
T=mat[1:,-1]
matopen= np.loadtxt("OPEN/SEQ_Steady_State_Open_T_SEG.son").T
Topen=matopen[1:,-1]
y=np.linspace(1.00000000e-05,1.01000000e-03,51)
plt.xlabel('y [mm]')
plt.ylabel('T [°C]')

plt.title("Temperature Probe at x=1.79mm")
plt.plot(y[:-1]*1000,T[:-1],label="Sym")
plt.plot(y[:-1]*1000,Topen[:-1],label="Open")
plt.legend(loc=0)
plt.savefig("plot.png")
plt.close()
