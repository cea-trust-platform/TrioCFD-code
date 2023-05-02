#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@author: mathis
"""
#%%
import matplotlib.pyplot as plt
import numpy as np
import sys
import os
cwd = os.getcwd()

#%%Fluid and liquid properties
case=sys.argv[1]
case_type=sys.argv[2]
post_process_all=bool(int(sys.argv[3]))
if case_type == 'CONVERGENCE':
    convergence = True #0 for seq/par study and 1 for convergence study
case_dir=sys.argv[4]
ldomain = float(sys.argv[5]) 
nbelem = float(sys.argv[6])
rho_0 = float(sys.argv[7]) 
Cp_0 = float(sys.argv[8]) 
lambda_0 = float(sys.argv[9])      
rho_1 = float(sys.argv[10])
Cp_1 = float(sys.argv[11]) 
lambda_1 = float(sys.argv[12])    
tf_scope=float(sys.argv[13])
t_study_tmp=sys.argv[14]
t_study=np.array(t_study_tmp.split(',')).astype(np.float) 

#Diffusion coefficient
alpha_0 = lambda_0/(rho_0*Cp_0)
alpha_1 = lambda_1/(rho_1*Cp_1)


#%% Data from Panda et al.
#Grid propeties
r_sphere = 1.0e-3 #same as Panda et al.
d_sphere = 2.0*r_sphere
grid_per_cells = np.array([[20.0],[40.0],[80.0]])
dr = d_sphere/grid_per_cells
domain = (ldomain/2)/r_sphere
#Choose the number of cells in a diameter index_dr = 0 for debug, index_dr = 2 for accuracy
index_dr = 2

r1 = np.arange(0,r_sphere,dr[index_dr])
r2 = np.arange(r_sphere+dr[index_dr],domain*r_sphere+dr[index_dr],dr[index_dr])
r = np.concatenate((r1,np.reshape(r_sphere,1,1),r2))
r_inf = r[-1]
r_norm = r/r_sphere

# theta = (T-T0)/(T1-T0) = (T-T_vap)/(T_liq-T_vap)
T0 = 373.0 #T_liq
T1 = 293.0 #T_vap

#Temporal discretisation
dt = 2e-4
nb_step = round(tf_scope/dt)+1
t = np.arange(0.0,nb_step*dt,dt)
index_t_study = np.zeros((len(t_study)),dtype = int)
for i in range (0,len(t_study)):
    index_t_study[i] = np.argmin(abs(t-t_study[i]))

#Initial condition
dT_dr = 0
T_inf = T1

#Initialisation champ T et matrice pour resolution
T = np.zeros((len(r),len(t)))
T[0:len(r1),0] = T0
T[len(r1)+1::,0] = T1

Beta = np.zeros((len(r),len(t)))
Beta[len(r)-1,::] = T1

A = np.zeros((len(r),len(r),len(t)))

alpha_eq_0 = -alpha_0/np.power(dr[index_dr],2)
alpha_eq_1 = -alpha_1/np.power(dr[index_dr],2)
coeff_T_star = 1/(8*(lambda_0+lambda_1))

#%% Construction de la matrice A pour chaque iteration (aucun terme ne depend du temps)
# 1ere iteration 1st order spatial (idee perso), apres second order spatial (central difference cf Das et al.)
# 1ere iteration 1st order temporal (idee perso), apres second order temporal (central Baker Oliphant cf Das et al.)
print(case_type + ' : 1D Finite difference START')
for i in range (0,len(t)):
    A[0,0:3,i] = dr[index_dr]*np.array([-3.0/2.0,2.0,-1.0/2.0]) #dT/dr (r = 0) = 0
    A[-1,-1,i] = 1 # T(r = r_inf) = T_inf = T1
    
    # Dans la sphere domaine vapeur (0)
    for j in range (1,len(r1)):
        A_N = alpha_eq_0*np.power(0.5*(r[j]+r[j-1]),2)/np.power(r[j],2)
        A_P = alpha_eq_0*np.power(0.5*(r[j]+r[j+1]),2)/np.power(r[j],2)
       
        if i == 0:
            A_C = 1/dt-A_N-A_P
        else:
            A_C = 1.5/dt-A_N-A_P
            
        A[j,j-1,i] = A_N
        A[j,j,i] = A_C
        A[j,j+1,i] = A_P
      
    # A l'interface : heat flux et temperature continus T_star-(..)/(..) = 0 cf Das et al.
    j = len(r1)
#2nd order with interface at dr/2 from n1
#A[j,j-2,i] = -coeff_T_star*lambda_0
#A[j,j-1,i] = coeff_T_star*9*lambda_0
#A[j,j,i] = -1
#A[j,j+1,i] = coeff_T_star*9*lambda_1
#A[j,j+2,i] = -coeff_T_star*lambda_1
#2nd order
#A[j,j-2,i] = 0.5*lambda_0
#A[j,j-1,i] = -2.0*lambda_0
#A[j,j,i] = 1.5*(lambda_0+lambda_1)
#A[j,j+1,i] = -2.0*lambda_1
#A[j,j+2,i] = 0.5*lambda_1
    
    #3rd order
    A[j,j-3,i] = -2.0*lambda_0
    A[j,j-2,i] = 9.0*lambda_0
    A[j,j-1,i] = -18.0*lambda_0
    A[j,j,i] = 11.0*(lambda_0+lambda_1)
    A[j,j+1,i] = -18.0*lambda_1
    A[j,j+2,i] = 9.0*lambda_1
    A[j,j+3,i] = -2.0*lambda_1

    
    # Dans le domaine liquide (1)
    for j in range (len(r1)+1,len(r)-1):
        A_N = alpha_eq_1*np.power(0.5*(r[j]+r[j-1]),2)/np.power(r[j],2)
        A_P = alpha_eq_1*np.power(0.5*(r[j]+r[j+1]),2)/np.power(r[j],2)
       
        if i == 0:
            A_C = 1/dt-A_N-A_P
        else:
            A_C = 1.5/dt-A_N-A_P            

        A[j,j-1,i] = A_N
        A[j,j,i] = A_C
        A[j,j+1,i] = A_P    
        
#%% Resolution au cours du temps
for i in range (0,len(t)-1):
#for i in range (0,1):        
    if i == 0:
        Beta[1:len(r1),i] = T[1:len(r1),i]/dt
        Beta[len(r1)+1:len(r)-1,i] = T[len(r1)+1:len(r)-1,i]/dt

    else:
        Beta[1:len(r1),i] = (4*T[1:len(r1),i]-T[1:len(r1),i-1])/(2*dt)
        Beta[len(r1)+1:len(r)-1,i] = (4*T[len(r1)+1:len(r)-1,i]-T[len(r1)+1:len(r)-1,i-1])/(2*dt)
        
    T[:,i+1] = np.matmul(np.linalg.inv(A[:,:,i]),Beta[:,i])    

print(case_type + ' : 1D Finite difference END')

#%% Affichage au cours du temps t0 = 0 t1 = 8e-4 t2 = 3.2e-3 t3 = 6.4e-3 t4 = 9.6e-3
T_norm = (T-T1)/(T0-T1)
    
#%%
if convergence:
    if post_process_all:
        process_type = np.array(['COARSE','MEDIUM','FINE'])
        T_interp_save_convergence=[]
        T_interp_time_save_convergence=[]
        r_probes_save_convergence=[]
        T_probes_save_convergence=[]
        rmse_save_convergence=[]
        r_squared_correction_save_convergence=[]
        elem=[]
    else:
        process_type = [case]
else:
    if post_process_all:
        process_type = np.array([i + j for i, j in zip([case_type + '_']*4,['seq','par8','repr','repr_par8'])])
    else:
        process_type = [case]
#Trace independant pour chaque cas

for k in range (0,len(process_type)):
    #%% Extracting data from .son
    if convergence:
        name_type_tmp = process_type[k] + ''.join(case_dir.partition("_")[1:])
        process_type_tmp = process_type[k]
    else:
        name_type_tmp = case_dir
        process_type_tmp = process_type[k]
    
    probes_names_tmp = process_type_tmp + '_SONDE_TEMP'
    fig1_name_tmp = name_type_tmp + '/' + process_type_tmp +'_TEMPERATURE' + '.png'
    fig2_name_tmp = name_type_tmp + '/' + process_type_tmp +'_ERRORS' + '.png'
    
    print(probes_names_tmp + '  START')
    #%% Open .son file
    direct = cwd + '/' + name_type_tmp + '/' + probes_names_tmp + '.son'    
    fichier = open(direct, 'r')
    text = fichier.read()
    
    #%% Split the text file
    N = 4 #4 first rows are not results
    caractere = '\n'
    caractere2 = ' '
    text_sp = text.split(caractere)
    if text_sp[-1]  ==  '':
        text_sp = text_sp[:len(text_sp)-1]
        
    #%% Probes positions
    text_position = text_sp[1].split(caractere2)
    m = 0
    x = []
    y = []
    z = []
    position = np.empty((3,1))
    while text_position[m] !='x=':
        m = m+1
    for i in range (m+1,len(text_position),6):
        x.append(float(text_position[i]))
        y.append(float(text_position[i+2]))
        z.append(float(text_position[i+4]))   
    x = np.array(x)
    y = np.array(y)
    z = np.array(z)
    
    #%% Gets the temperature over iterations
    t_simu = []
    T_simu = np.empty((len(x),len(t_study)))
    #Get the index at each time step
    for i in range (N,len(text_sp)):
        t_simu.append(float(text_sp[i].split(caractere2)[0]))
    t_simu = np.array(t_simu)
    #Find the time steps matching with 1D calculations
    index_t_simu = np.zeros((len(t_study)),dtype = int)
    for i in range (0,len(t_study)):
        index_t_simu[i] = np.argmin(abs(t_simu-t_study[i]))
    t_simu_display = t_simu[index_t_simu]
    # Extracing Temperature value at this index
    for i in range (0,len(index_t_simu)):
        T_simu[:,i] = np.array(text_sp[index_t_simu[i]+N].split(caractere2)[1::]).astype(np.float)   
                  
    #%% Normalisation des resultats 3D
    r_probes = (x-r[-1])/r_sphere
    T_probes = (T_simu-T1)/(T0-T1)   
    
    T_probes = T_probes[r_probes>= 0,:]
    r_probes = r_probes[r_probes>= 0]
    
    #%% Interpolation en temps car la simu TrioIJK ne tombe pas parfaitement sur un multiple de dt
    
    T_interp_time = np.empty((len(r),len(t_study)))
    for i in range(0,len(T_interp_time[:,0])):
        T_interp_time[i,:] = np.interp(t_simu_display, t, T_norm[i,:], left = None, right = None, period = None)
          
    #%% Interpolation des resultats 1D sur les coordonnes des sondes 3D
    T_interp = np.empty((len(r_probes),len(t_study)))
    for i in range(0,len(T_interp[0,:])):
        #T_interp[:,i] = np.interp(r_probes, r_norm, T_norm[:,i], left = None, right = None, period = None)
        T_interp[:,i] = np.interp(r_probes, r_norm, T_interp_time[:,i], left = None, right = None, period = None)
    
    #%% Calcul ecart type sur tout le domaine
    ecart_abs = abs(T_interp-T_probes)
    ssres = np.sum(ecart_abs**2,axis = 0)
    probes_mean = np.mean(T_interp,axis = 0)
    probes_mean = probes_mean.reshape(1,len(probes_mean))
    probes_mean = np.repeat(probes_mean,len(T_probes[:,0]),axis = 0)
    sstot = np.sum((T_probes-probes_mean)**2,axis = 0)    
    
    # MSE, RMSE, R^2
    mse = np.mean(ecart_abs**2,axis = 0)
    rmse = np.sqrt(mse)   
    r_squared = 1-ssres/sstot
    
    #%% Correction en ne prenant en compte que le champ proche de la bulle
    r_squared_correction = np.empty(len(t_simu_display))
    ssres_correction = np.empty(len(t_simu_display))
    sstot_correction = np.empty(len(t_simu_display))
    #for i in range(0,1):
    for i in range(0,len(T_interp[0,:])):
        T_interp_stats = T_interp[T_probes[:,i]>1e-2,i]
        T_probes_stats = T_probes[T_probes[:,i]>1e-2,i]
        ssres_correction[i] = np.sum((T_interp_stats-T_probes_stats)**2,axis = 0)
        sstot_correction[i] = np.sum((T_probes_stats-np.mean(T_probes_stats))**2,axis = 0)  
        r_squared_correction[i] =  1-ssres_correction[i]/sstot_correction[i]
    
    #%% Trace des courbes a differents instants (1D et 3D) 
    fig1=plt.figure(figsize=(5,5))
    
    #Trace de la condition initiale
    T_ini=np.concatenate((T0*np.ones((len(r1)+1)),T1*np.ones((len(r)-(len(r1))))))
    T_ini_norm=(T_ini-T1)/(T0-T1)
    #plt.plot(np.concatenate((r_norm[:len(r1)+1],r_norm[len(r1):])),T_ini_norm , label='$t=0$', marker='', linestyle='--', color  = 'black',linewidth=1)
    plt.axvspan(0, 1, alpha=0.25, color='grey')
    color_str = ['red','green','blue','black']
    plt.plot([],[],label = '1D-FD', marker = '', linestyle = '-',color = 'grey')
    plt.plot([],[],label = 'TrioIJK', marker = '+', linestyle = '',color = 'grey')
    for i in range (0,len(t_simu_display)):
        legend_str = ('t = ' + str('{:.1e}'.format(t_simu_display[i])))
        plt.plot(r_norm,T_interp_time[:,i] , label = legend_str, marker = '', linestyle = '-',color = color_str[i],linewidth=0.5)  # Plot some data on the (implicit) axes.
        plt.plot(r_probes,T_probes[:,i] , label = '_nolegend_', marker = '+', linestyle = '',color = color_str[i]) 
    plt.xlim(0, domain) 
    plt.ylim(-0.2,1.2)
    plt.axes().set_aspect('auto')
    plt.xlabel(' Non-dimensionalised radius $\widetilde{r}$ (-) ')
    plt.ylabel(' Non-dimensionalised temperature $\Theta$ (-)')
    plt.title('Transient non-dimensionalised radial temperature profiles: \n sphere kept in a pool of continuous phase \n(' + str(int(nbelem))  + ' elements)')
    plt.legend(loc = 'upper right', frameon = True, shadow = True, borderpad = 0.5)
    plt.grid(b = True, which = 'major', color = '#666666', linestyle = '-', alpha = 1)
    plt.minorticks_on()
    plt.grid(b = True, which = 'minor', color = '#666666', linestyle = '-', alpha = 0.2)
    plt.savefig(fig1_name_tmp,dpi = 300,bbox_inches = 'tight')  
    #plt.show()
    plt.close(fig1)

    #%% Trace de l'ecart type et valeur

    fig2, ax1 = plt.subplots(figsize = (5,5))
    ax1.tick_params(axis = 'y', labelcolor = 'red')
    ax1.loglog(t_simu_display, rmse , label = '_nolegend_', marker = '+', linestyle = '-',color = 'red',linewidth=0.5)
    ax1.plot([],[],label = 'RMSE', marker = '+', linestyle = '-',color = 'red') 
    ax1.plot([],[],label = '$R^2$', marker = '+', linestyle = '--',color = 'blue')

    ax1.set_xlim(t_simu_display[0]*0.0001, t_simu_display[-1]*1.1)
    #ax1.set_ylim(0,1)        
    ax1.set_xlabel('Time t (s)')
    ax1.set_ylabel('Root mean square error RMSE')        
    
    plt.grid(b = True, which = 'major', color = '#666666', linestyle = '-', alpha = 1)
    plt.minorticks_on()
    plt.grid(b = True, which = 'minor', color = '#666666', linestyle = '-', alpha = 0.2)
    
    ax2 = ax1.twinx()
    ax2.set_ylabel('Determination coefficient $R^2$')
    ax2.tick_params(axis = 'y', labelcolor = 'blue')   
    #ax2.set_ylim(0,0.1)
    ax2.semilogx(t_simu_display, r_squared_correction , label = '_nolegend_', marker = '+', linestyle = '--',color = 'blue',linewidth=0.5)
        
    ax1.set_title('Temperature residuals between 1D-FD and TrioIJK simulations: \n sphere kept in a pool of continuous phase \n(' + str(int(nbelem)) + ' elements)')
    ax1.legend(loc = 'lower left', frameon = True, shadow = True, borderpad = 0.5)
    
    #plt.ticklabel_format(style='sci', axis='x', scilimits=(-4,-3))
    plt.xlim([t_simu_display[0]/10, t_simu_display[-1]*10])	
    plt.savefig(fig2_name_tmp,dpi = 300,bbox_inches = 'tight')  
    #plt.show()
    plt.close(fig2)

    print(probes_names_tmp + '  END')

    #%% Enregistrement tmp de MEDIUM et COARSE pour figures commune 
    if convergence and post_process_all:
        T_probes_save_convergence.append(T_probes)
        T_interp_save_convergence.append(T_interp)
        T_interp_time_save_convergence.append(T_interp_time)
        r_probes_save_convergence.append(r_probes)
        rmse_save_convergence.append(rmse)
        r_squared_correction_save_convergence.append(r_squared_correction)

#%% Trace figure temperature commune
if convergence and post_process_all:   
    fig3_name_tmp = case_type + '_' +'TEMPERATURE' + '.png'
    fig3=plt.figure(figsize=(5,5))
    plt.axvspan(0, 1, alpha=0.25, color='grey')
    color_str = ['red','green','blue','black']
    marker_tab = ['+','1','x']
    plt.plot([],[],label = '1D-FD', marker = '', linestyle = '-',color = 'grey')
    for j in range (0,len(T_interp_save_convergence)):
        plt.plot([],[],label = process_type[j], marker = marker_tab[j], linestyle = '',color = 'grey')
    for j in range (0,len(T_interp_save_convergence)):
        for i in range (0,len(t_simu_display)):
            if j==0:
                legend_str = ('t = ' + str('{:.1e}'.format(t_simu_display[i])))
                plt.plot(r_norm,T_interp_time_save_convergence[j][:,i] , label = legend_str, marker = '', linestyle = '-',color = color_str[i],linewidth=0.5)  # Plot some data on the (implicit) axes.
            plt.plot(r_probes_save_convergence[j],T_probes_save_convergence[j][:,i] , label = '_nolegend_', marker = marker_tab[j], linestyle = '',color = color_str[i])      
    plt.xlim(0, domain) 
    plt.ylim(-0.2,1.2)
    plt.axes().set_aspect('auto')
    plt.xlabel(' Non-dimensionalised radius $\widetilde{r}$ (-) ')
    plt.ylabel(' Non-dimensionalised temperature $\Theta$ (-)')
    plt.title('Transient non-dimensionalised radial temperature profiles: \n sphere kept in a pool of continuous phase \n(Convergence study)')
    plt.legend(loc = 'upper right', frameon = True, shadow = True, borderpad = 0.5)
    plt.grid(b = True, which = 'major', color = '#666666', linestyle = '-', alpha = 1)
    plt.minorticks_on()
    plt.grid(b = True, which = 'minor', color = '#666666', linestyle = '-', alpha = 0.2)
    plt.savefig(fig3_name_tmp,dpi = 300,bbox_inches = 'tight')  
    #plt.show()
    plt.close(fig3)
    
#%% Trace figure erreur commune    
    rmse_mean=np.empty((len(T_interp_save_convergence),))
    r_squared_mean=np.empty((len(T_interp_save_convergence),))
    for j in range (0,len(T_interp_save_convergence)):
        elem.append(len(r_probes_save_convergence[j]))
        rmse_mean[j]=np.mean(rmse_save_convergence[j])
        r_squared_mean[j]=np.mean(r_squared_correction_save_convergence[j])
    
    fig4_name_tmp = case_type + '_' +'ERRORS' + '.png'
    fig4, ax1 = plt.subplots(figsize = (5,5))
    ax1.tick_params(axis = 'y', labelcolor = 'red')
    for j in range (0,len(r_squared_mean)):
        ax1.loglog(elem, rmse_mean , label = '_nolegend_', marker = '+', linestyle = '-',color = 'red',linewidth=0.5)
        
    ax1.plot([],[],label = 'RMSE', marker = '', linestyle = '-',color = 'red')
    ax1.plot([],[],label = '$R^2$', marker = '', linestyle = '--',color = 'blue')
    #ax1.set_ylim(0,1)        
    ax1.set_xlabel('Number of elements (-)')
    ax1.set_ylabel('Root mean square error RMSE')
      
    
    plt.grid(b = True, which = 'major', color = '#666666', linestyle = '-', alpha = 1)
    plt.minorticks_on()
    plt.grid(b = True, which = 'minor', color = '#666666', linestyle = '-', alpha = 0.2)
    
    ax2 = ax1.twinx()
    ax2.set_ylabel('Determination coefficient $R^2$') 
    ax2.tick_params(axis = 'y', labelcolor = 'blue')   
    #ax2.set_ylim(0,0.1)
    for j in range (0,len(r_squared_mean)):
        ax2.semilogx(elem, r_squared_mean, label = '_nolegend_', marker = '+', linestyle = '--',color = 'blue',linewidth=0.5)
        
    ax1.set_title('Temperature residuals between 1D-FD and TrioIJK simulations: \n sphere kept in a pool of continuous phase\n(Convergence study)')
    ax1.legend(loc = 'lower left', frameon = True, shadow = True, borderpad = 0.5)
    plt.xlim([min(elem)/10,max(elem)*10])	
    #plt.ticklabel_format(style='sci', axis='x', scilimits=(-4,-3))
    plt.savefig(fig4_name_tmp,dpi = 300,bbox_inches = 'tight')  
    #plt.show()
    plt.close(fig4)
