'''
This script converts the .txt files into .dt_ev using BuildStat and checks the values over time that should remain unchanged in this steady state case. Finally, the last snapshot is compared to analytical values of the spatial average (x,y directions) at each z plane coordinates
'''
##################
# Import libraries
##################

import os
import sys
import glob
import sympy as sym
from sympy.parsing.sympy_parser import parse_expr
from sympy.parsing.sympy_parser import standard_transformations,implicit_multiplication_application
transformations = (standard_transformations + (implicit_multiplication_application,))

##########################
# Load Tools from PyTools3 
##########################

pjd = os.getenv("project_directory")
commons = pjd + 'share/PyTools3'

from commons.DNSTools import *
from commons.BuildStat import *

cur_dir=os.getcwd() + '/'
name_out=['Unsteady','Final_snapshot','Stats']

#####################################################
# diphasique_statistiques_*.txt into the dt_ev format
#####################################################

#No init or end file specified leads to calculate the stats between the first and last files
#diph = StatReader(temperature=True)
print('     -> Transform the .txt files into .dt_ev and sort them into folders')
#Look at the list of stats .txt files in the working dir
fic_stats = glob.glob("*phasique_statistiques_[0-9]*.txt")
# Sort the file by time 
fic_stats.sort(key=lambda item: float(item.strip("monodiphasique_statistiques_.txt")))
# Get the time of each file
time = [float(fic_stats[i].strip("monodiphasique_statistique_.txt")) for i in range(0,len(fic_stats))]

for i in range(1,len(fic_stats)): 
    # Option to specify an initial or final file
    options = process_cmd_line(['-i', fic_stats[0],'-e',fic_stats[i]])
    diph = StatReader(finstatfile=options.endfile, inistatfile=options.inifile, temperature=True)
    #Rename all dt_ev
    diph.saveStats()
    for name in name_out:
        name_dtev = glob.glob(name+'*.dt_ev')
        for name_save in name_dtev:
            #print(name_save)
            name_dir=name_save.replace('.dt_ev','')
            if not os.path.exists(name_dir):      
                os.mkdir(name_dir)
            os.rename(cur_dir + name_save, cur_dir+ name_dir + '/'+ name_dir+'_'+str(time[0])+'_to_'+str(time[i])+'.dt_ev') 
 
############################################################################
# Declare few functions for error calculation from theroretical field values 
############################################################################

def getDomain(jdd=None):
    if jdd == None:
        jdd = getJddName() + ".data"
    Lx = getValue("uniform_domain_size_i", jdd)
    Ly = getValue("uniform_domain_size_j", jdd)
    Lz = getValue("uniform_domain_size_k", jdd)
    nx = getValue('nbelem_i',jdd)
    ny = getValue('nbelem_j',jdd)
    nz = getValue('nbelem_k',jdd)
    return np.array([Lx, Ly, Lz, nx, ny, nz])

def getVelTempInit(jdd=None):
    if jdd == None:
        jdd = getJddName() + ".data"
    vx = getParam(fic=jdd,name='expression_vx_init',string=True)
    vy = getParam(fic=jdd,name='expression_vy_init',string=True)
    vz = getParam(fic=jdd,name='expression_vz_init',string=True)
    T_xyz = getParam(fic=jdd,name='expression_T_init',string=True)
    return np.array([vx, vy, vz, T_xyz])

def intStats(list_snapshot_header,list_snapshot_temp):
    x = sym.symbols('x')
    y = sym.symbols('y')
    z = sym.symbols('z')
    [vx, vy, vz, T_xyz] = getVelTempInit()
    [Lx, Ly, Lz, nx, ny, nz] = getDomain()
    # Transform the str into equations
    T = parse_expr(T_xyz, transformations=transformations)
    U = parse_expr(vx, transformations=transformations) 
    V = parse_expr(vy, transformations=transformations)
    W = parse_expr(vz, transformations=transformations)
    I = sym.Function('1')
    # Check the expressions
    print('Vx_init = ' + str(U))
    print('Vy_init = ' + str(V))
    print('Vz_init = ' + str(W))
    print('T_init = ' + str(T))
    # Expression involving T,u,v,w combined
    f=[]
    for i in range(1,12):
        f_str=list_snapshot_header[i].replace('T',str(T)).replace('U',str(U)).replace('V',str(V)).replace('W',str(W)).replace('I',str(I))
        f.append(parse_expr(f_str,transformations=transformations))
    print('Functions calculated : ' + str(f))
    # Averaging on x,y directions at the last time step, at every z coordinates 
    stats_from_init=np.empty((int(nz),len(f)))
    for i in range(0,int(nz)): 
        for j in range(0,len(f)):
            stats_from_init[i,j]=float(sym.integrate(f[j].subs(z,list_snapshot_temp[i,0]),(x,0,Lx),(y,0,Ly))/(Lx*Ly))
    # Calculate error for each variable and at each coordinates z
    err=np.abs(list_snapshot_temp[:,1:12]-stats_from_init)
    err_mean=np.mean(err,axis=0)
    err_rel_mean=np.mean(err/stats_from_init*100.0,axis=0)
    threshold_abs=1e-10
    threshold_rel=5
    err_bool = (err_mean<threshold_abs) + (err_rel_mean<threshold_rel)
    if 'z'not in str(U)+str(V)+str(W)+str(T):
        err_header=np.array(['Parameters','Functions','Numerical','Theoretical','AbsErrors','RelErrors(percent)','Test'])
        err_out=np.transpose(np.array([list_snapshot_header[1:12],np.array(f).astype(str),list_snapshot_temp[0,1:12].astype(str),stats_from_init[0,:].astype(str),err_mean.astype(str),err_rel_mean.astype(str),err_bool.astype(str)]))
    else:
        err_header=np.array(['Parameters','Functions','AbsErrors','RelErrors(percent)','Test'])
        err_out=np.transpose(np.array([list_snapshot_header[1:12],np.array(f).astype(str),err_mean.astype(str),err_rel_mean.astype(str),err_bool.astype(str)]))
    err_out=np.r_[[err_header],err_out]
    print(err_out)
    return(err_out)

###################################################################
# .txt to tex table in .prm 
###################################################################

def txtToTexTable(name,adjust):
    file_err=open(name,'r')
    file_txt=file_err.read()
    n_column=len(file_txt.rstrip().lstrip().split('\n')[0].split(' '))
    c_column='|c'*n_column
    file_tex='\\begin{tabular}{'+str(c_column)+'|} \hline $'
    file_tex=file_tex + file_txt.rstrip().lstrip().replace(' ','$ & $').replace('\n','$ \\\ \hline $')
    file_tex=file_tex + '$ \\\ \hline \\end{tabular}'
    if adjust:
        file_tex='\\begin{center}\\resizebox{\\textwidth}{!}{'+file_tex+'}\\end{center}'
    else:
        file_tex='\\begin{center}'+file_tex+'\\end{center}'
    prm = getJddName() + ".prm"
    reading_prm = open(prm, "r")
    modify_prm = reading_prm.read().replace(name,file_tex)
    writing_prm = open(prm,"w")
    writing_prm.write(modify_prm)

####################################################################
# Get data from .dt_ev, calculate the value analytically and compare 
####################################################################
print('     ------------------------------------------')
print('     ------------------------------------------')
print('\n   -> Check the statistics of temperature')
for name in name_out:
    name_dir = glob.glob(name+'*')
    for name_dir_tmp in name_dir:
        # Retrieve the .dt_ev
        name_files=glob.glob(name_dir_tmp+'/*.dt_ev')
        ini=[]
        end=[]
        for i in range(0,len(name_files)):
            name_files[i]=name_files[i].replace('.dt_ev','')
            name_split_tmp=name_files[i].split('_')
            ini.append(float(name_split_tmp[-3]))
            end.append(float(name_split_tmp[-1]))        
        
        #Then sort depending on tfin
        files_end=list(zip(end,name_files))
        name_files = [name_files for end, name_files in files_end]
        end.sort()
        if 'temperature' in name_dir_tmp:
            #Check statistics are constants
            if name == 'Stats':
                print('\n       -> Statistics temperature')
                list_array=[]
                for name_stats in name_files:
                    # Get values from dt_ev
                    list_stats_temp_file=open(name_stats+'.dt_ev','r')
                    # Split string into a tab
                    list_stats_temp=list_stats_temp_file.read().replace('#','').rstrip().lstrip().split('\n')
                    # Isolate the header
                    list_stats_header=np.array(list_stats_temp[0].rstrip().lstrip().split(' '))
                    # Put the values into an array
                    list_stats_temp=list_stats_temp[1:]
                    for i in range (0,len(list_stats_temp)):
                        list_stats_temp[i]=list_stats_temp[i].rstrip().lstrip().split(' ')                     
                    list_array.append(np.array(list_stats_temp).astype(float))                            
                list_array=np.array(list_array)
                # Calulate the differences to the initial values, should tend towards zero (constant fields)
                # One tab for each z coordinates
                diff_k=np.sum(np.abs(list_array[1:]-list_array[0]),axis=0)
                diff=np.sum(diff_k,axis=0)
                #print(diff)
                threshold_diff=1e-5
                diff_bool = (diff<threshold_diff)
                diff_header=np.array(['Parameters','FirstRow','SumErrorK','Test'])
                #print(list_array[0][0])
                diff_out=np.transpose(np.array([list_stats_header,list_array[0][0].astype(str),diff.astype(str),diff_bool.astype(str)]))
                diff_out=np.r_[[diff_header],diff_out]
                print(diff_out)
                np.savetxt(name_dir_tmp+ '/' + name_dir_tmp + '_err.out',diff_out,fmt='%s')
                txtToTexTable(name_dir_tmp+ '/' + name_dir_tmp + '_err.out',False)

           # Check average in the plane
            elif name == 'Final_snapshot':
                print('\n       -> Final snapshot temperature')
                # Get the last file 
                list_snapshot_temp_file=open(name_files[-1]+'.dt_ev','r')
                # Split string into a tab
                list_snapshot_temp=list_snapshot_temp_file.read().replace('#','').rstrip().lstrip().split('\n')
                # Isolate the header
                list_snapshot_header=np.array(list_snapshot_temp[0].rstrip().lstrip().split(' '))
                # Put the values into an array
                list_snapshot_temp=list_snapshot_temp[1:]
                for i in range (0,len(list_snapshot_temp)):
                    list_snapshot_temp[i]=list_snapshot_temp[i].rstrip().lstrip().split(' ') 
                list_snapshot_temp=np.array(list_snapshot_temp).astype(float)
                # Calculate the errors
                err_out = intStats(list_snapshot_header,list_snapshot_temp)
                np.savetxt(name_dir_tmp+ '/' + name_dir_tmp + '_err.out',err_out,fmt='%s')
                txtToTexTable(name_dir_tmp+ '/' + name_dir_tmp + '_err.out',True)

            # Unsteady not used 
            else:
                print('\n       -> Unsteady temperature not used')

####################################################################
# Modify the prm data
####################################################################   

[vx, vy, vz, T_xyz] = getVelTempInit()
prm = getJddName() + ".prm"
reading_prm = open(prm, "r")
modify_prm = reading_prm.read().replace('vx_ana',vx).replace('vy_ana',vy).replace('vz_ana',vz).replace('T_ana',T_xyz)
writing_prm = open(prm,"w")
writing_prm.write(modify_prm)


