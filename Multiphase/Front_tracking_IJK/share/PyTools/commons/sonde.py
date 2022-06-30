#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
# from scipy import *
from scipy import signal
from pylab import *
import matplotlib.pyplot as plt

direct='DNS_PP_VX.son'
col=[1,2,3,4]

fichier=open(direct, "r")
text=fichier.read() 

N=len(open(direct, "rU").readlines())
N=N-4

i=0
cx=[]
cy=[]
cz=[]
c=[]
a=[]

caractere = "\n";
caractere2=" ";
for x in text.split(caractere):
	i+=1
	if i==1:
		j=0
		for y in x.split(caractere2):
			try:
				str(y)
				cx.append(y)						
	else:
		pass




vec= map(float,a)
coord_x=map(float,cx)
coord_y=map(float,cy)
coord_z=map(float,cz)


N_sonde=len(coord_x)
lon=N
lar=N_sonde+1
s = (lon,lar)
mat=np.zeros(s)


k=0
for x in vec:
	i=k/lar
	j=k%lar
	mat[i,j]=x
	k+=1
	

for x in col:
	
	figure(1)
	plot(mat[:,0],mat[:,x], label='En x= %s, \n y=%s, z=%s' %(str(coord_x[x-1]),str(coord_y[x-1]),str(coord_z[x-1])))
	plt.title('%s' %direct)
	plt.legend(prop={'size':10})
	plt.xlabel('Temps')
	plt.ylabel('vitesse')
		
savefig('sondes.png')



def moyenne(m):
	return sum(m)/len(m)

def spectre(col):
	global coord_x
	global coord_y
	global coord_z
	global mat
	global direct

	for x in col:
	
		SamplingFreq=100000      
		Fluct=[]
	
		m=moyenne(mat[:,x])
		Fluct[:] = mat[:,x]-m
		f, Pxx = signal.periodogram(Fluct[:], SamplingFreq)
	
		figure(2)
		loglog(f,Pxx, label='En x= %s, \n y=%s, z=%s' %(str(coord_x[x-1]),str(coord_y[x-1]),str(coord_z[x-1])))
		plt.title('%s' %direct)
		plt.legend(prop={'size':9})
		plt.axis([10, 100000,0.000000000001,0.01])
		plt.xlabel('Frequence (Hz)')
		plt.ylabel('Energie')
			
	y=np.linspace(1, 10000, num=100)		
	loglog(y,y**(-3), label='Theoretical energy cascade slope (-3)')	
	plt.legend(prop={'size':9})
	savefig('spectre.png')

	return 
	
spectre(col)







