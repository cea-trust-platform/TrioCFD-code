import glob
import sys
import datetime
import numpy as np
import math as mt
import matplotlib.pyplot as plt
import DNSTools3 as dtool

print("Execute le :"+str(datetime.date.today().strftime("%d %b. %Y")))
print("Execute a : "+str(datetime.datetime.now().strftime("%H:%M:%S")))
print("Commandes :"+str(sys.argv[:]))

""" plotter la distance minimale bulle-bulle au cours du temps"""
# emploi : plot_min_dist.py jdd_radical nb_skiprows mesh_name

jdd=sys.argv[1]
sr=int(sys.argv[2])
try: mesh_name=sys.argv[3]
except: mesh_name=""
fichier="_dmin.out"

try:    data = np.loadtxt(jdd+fichier,usecols=(0,1,2),skiprows=sr)
except:    data = np.loadtxt("OUT/"+jdd+fichier,usecols=(0,1,2),skiprows=sr)
try:    l_fr = dtool.getParam(jdd+".data", "portee_force_repulsion")
except:    pass
    
print(data.shape)
if (data.shape[0])/3.!=int((data.shape[0])/3):
    print(data.shape[0]/3.)
    print(int(data.shape[0]/3.))
    data = np.insert(data,0,0,axis=0)
    data = np.insert(data,0,0,axis=0)

print(data.shape)
print(data[2,:])
i = data[:,0]
t = data[:,1]
d = data[:,2]

# nb_points = int((len(i)+2)/3)
nb_points = int((len(i))/3)
print(nb_points)
good_i = 3*np.linspace(0,i[-1],nb_points)
# t = t.reshape(nb_points,3);t = t[:,2]
# d = d.reshape(nb_points,3);d = d[:,2]
t = t.reshape(-1,3);t = t[:,2]
d = d.reshape(-1,3);d = d[:,2]

L=dtool.getParam(jdd+".data", "uniform_domain_size_i")
n=dtool.getParam(jdd+".data", "nbelem_i")
dx=L/n

fig, (ax1) = plt.subplots(1,1)
ax1.semilogy(t,d,color='r',label=r'$d_{min}$ bulles')                                      
ax1.axhline(y=dx/1,color='k',linestyle='dashed',label=mesh_name+" mesh")                       
# ax1.axhline(y=dx/2,color='k',linestyle='dotted',label="med mesh")                          
# ax1.axhline(y=dx/4,color='k',linestyle='dashdot',label="fine mesh")                        
try:    ax1.axhline(y=l_fr,color='k',linestyle='dashdot',linewidth=0.6,label="repulsion force range")      
except: pass  
ax1.set_ylabel('distance')
ax1.set_xlabel('temps')
ax1.legend(framealpha=0.);

fig.tight_layout();
plt.savefig("d_min_semi_log_y.eps")
plt.close()

fig, (ax1) = plt.subplots(1,1)
ax1.plot(t,d,color='r',label=r'$d_{min}$ bulles')                                                      
ax1.axhline(y=dx/1,color='k',linestyle='dotted',label="coarse mesh")                                   
ax1.axhline(y=dx/2,color='k',linestyle='dashed',label="med mesh")                                      
ax1.axhline(y=dx/4,color='k',linestyle='dashdot',label="fine mesh")                                    
try:    ax1.axhline(y=l_fr,color='k',linestyle='dashdot',linewidth=0.6,label="repulsion force range")      
except: pass  
ax1.set_ylabel('distance')
ax1.set_xlabel('temps')
ax1.legend(framealpha=0.);

fig.tight_layout();
plt.savefig("d_min_linlin.eps")
