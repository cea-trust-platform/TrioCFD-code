#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 23 11:03:38 2023

@author: Moncef El Moatamid
"""
from string import Template
import numpy as np
from math import *
from scipy.optimize import newton, brentq

# # paramètres fixés pour tous les calculs
tmax = 2 # [s]

def find_zero(x_sonde,velocity_sonde,x_min, x_max):
    f = lambda x : np.interp(x,x_sonde,velocity_sonde)
    try:
        zero = brentq(f,x_min, x_max)
        #zero = newton(f,0.2)
    except:
        zero = 0
    return zero


###############################################################################
######################### Substitution in data file ###########################
###############################################################################

# #
# # substituion 
# #

# Mot clé de l'équation à utiliser pour la dissipation (tau ou omega)
equation = ["echelle_temporelle_turbulente"]+["taux_dissipation_turbulent"]
# mot clé : terme de diffusion dans l'eq de qdm
diffusion = [ "k_tau"]+["k_omega"]
# Mot clé terme de dissipation turbulente
diss =         ["tau"]+['omega']

diss_ext = ["tau_ext"]+["omega_ext"]
# Condition de convergence du calcul sur le terme de dissipation turbulente
diss_conv= [" tau 1.e-5 "]+[" omega 1. "]
# Conditions initiales pour tau et omega
IC_diss =      ["0.00548"]+["182.399"]
# Condition limite à la paroi pour l'eq de transport de la dissipation turbulente (tau et omega)
CL_diss = ["scalaire_impose_paroi champ_front_uniforme 1 0 ", "Cond_lim_omega_demi { }  "]
# Condition limite à la paroi pour l'eq de transport k
CL_k = ["Cond_lim_k_simple_flux_nul ", "Cond_lim_k_simple_flux_nul  "]
# terme source supp dans l'équation de transport de tau (dissipation turbu)
diffusion_sup = [" "] *2+[" , Diffusion_supplementaire_echelle_temp_turb "]+[" "]+[" , Diffusion_supplementaire_echelle_temp_turb "]
# Autres param de simu
facsec     = 1
nb_pas_dt_max = "500"
options_vdf = ["option_vdf { all_options }", "option_vdf { all_options }", "", "", ""]

params_pb_multi = [equation, diffusion, diss, diss_conv, IC_diss, CL_k, diffusion_sup, facsec, nb_pas_dt_max]


def get_value_at_x(var,x,x_tar):
    # interpolation linéaire
    
    return np.interp(np.array(x_tar), x, var)

def method_name(config):
    # identifier le mot clé de la discretisation utilisée à partir du nom de la config entrée par l'utilisateur
    names = []
    for method in config:
        if "k-epsilon" in method:
            name = "VDF"
            if "VEF" in method:
                name = "VEFPreP1b"
            names.append(name)
        else:
            name = "VDF"
            if "polymac" in method:
                name = "PolyMAC_P0"
            names.append(name)
    return names

def substitution(config,mesh_param,fluid_param):
    Nx, Ny, growth_rate_x_inlet, growth_rate_y_inlet, growth_rate_x_step, growth_rate_y_step = mesh_param
    mu, rho = fluid_param
    methods = method_name(config)

    dict_list = [[0 for x in range(len(Nx))] for y in range(len(config))]
    for i in range(len(config)):
        for j in range(len(Nx)):
            options_vdf = ""
            diffusion_sup = ""
            if config[i] == "VEF_k-epsilon": # si VEF alors deux fois moins de maille dans chaque direction car "trianguler_h" coupe les mailles en 4
                dict_list[i][j] = {"Nx_inlet" : int(2*Nx[j]/3/2)+1,
                    "Ny_inlet" : int(Ny[j]/2)+1,
                    "Nx_step" : int(Nx[j]/3/2)+1,
                    "Ny_step" : int(Ny[j]/4/2)+1,
                    "fact_x_inlet" : growth_rate_x_inlet[j],
                    "fact_y_inlet" : growth_rate_y_inlet[j],
                    "fact_x_step" : growth_rate_x_step[j],
                    "fact_y_step" : growth_rate_y_step[j],

                    "tmax" : tmax,
                    "nb_pas_dt_max" : nb_pas_dt_max,
                    "mu" : mu,
                    "rho" : rho}
                continue
            diffusion_i, diss_i, diss_ext_i, diss_conv_i, IC_diss_i, CL_diss_i, equation_i = diffusion[0], diss[0], diss_ext[0], diss_conv[0], IC_diss[0], CL_diss[0], equation[0]
            solveur_temp = 'sets'
            if "k-omega" in config[i]:
                diffusion_i, diss_i, diss_ext_i, diss_conv_i, IC_diss_i, CL_diss_i, equation_i = diffusion[1], diss[1], diss_ext[1], diss_conv[1], IC_diss[1], CL_diss[1], equation[1]
            if "vdf" in config[i]:
                options_vdf = "option_vdf { all_options }"
            if "k-tau" in config[i]:
                diffusion_sup = " , Diffusion_supplementaire_echelle_temp_turb "
                solveur_temp = 'ice'
            dict = {"Nx_inlet" : int(2*Nx[j]/3)+1, # deux tiers des mailles avant la marche
                    "Ny_inlet" : Ny[j]+1,
                    "Nx_step" : int(Nx[j]/3)+1, # un tiers des mailles après la marche
                    "Ny_step" : int(Ny[j]/4)+1, # un quart des mailles selon y dans la marche (choix arbitraire)
                    "fact_x_inlet" : growth_rate_x_inlet[j],
                    "fact_y_inlet" : growth_rate_y_inlet[j],
                    "fact_x_step" : growth_rate_x_step[j],
                    "fact_y_step" : growth_rate_y_step[j],

                    "method" : methods[i],
                    "options_vdf" : options_vdf,
                    "tmax" : tmax,
                    "nb_pas_dt_max" : nb_pas_dt_max,
                    "solveur_temp" : solveur_temp,
                    "diss_conv" : diss_conv_i,

                    "mu" : mu,
                    "rho" : rho,

                    "diffusion" : diffusion_i,
                    "equation" : equation_i,
                    "diss" : diss_i,
                    "IC_diss" : IC_diss_i,
                    "CL_diss" : CL_diss_i,
                    "diss_ext" : diss_ext_i,
                    "diffusion_sup" : diffusion_sup
                    }
            dict_list[i][j] = dict
    return dict_list


def GenerateInputFile(dir,build,substitutions_dict,name, new_name):
    with open(f"{build}/{name}.data", "r") as file: 
        filedata = Template(file.read())
    result = filedata.substitute(substitutions_dict)
    with open(dir+f"/marche_{new_name}.data", "w") as file:
        file.write(result)
        
###############################################################################
################################# Read files ##################################
###############################################################################

def read_facefile(file):
    """ read data in .face file at différent time values and stores the data 
    for each value of t as a matrix in a dictionnary with the corresponding
    time value as the dictionnary key"""
    data = {}
    with open(file) as f:
       for l in f:
           if l.startswith("Temps"):
               key = float(l.strip("\n").split(" : ")[-1])
               data[key] = []
           line = l.strip().split("\t|")
           try:
               data[key].append([float(s) for s in line])
           except:
               pass
    for k in data:
        data[k] = np.array(data[k])
    return data

def get_wall_param(array):
    # seperate the data in .face files into different arrays
    array = np.transpose(array)
    x = array[0]
    y = array[1]
    u_plus = array[2]
    d_plus = array[3]
    u_star = array[4]
    tau = array[5] # m²/s²
    tau_x = array[6] # m²/s²
    tau_y = array[7] # m²/s²
    sort_index = np.argsort(x)
    x, y, u_plus, d_plus, u_star, tau, tau_x, tau_y = x[sort_index], y[sort_index], u_plus[sort_index], d_plus[sort_index], u_star[sort_index], tau[sort_index], tau_x[sort_index], tau_y[sort_index]
    return x, y, u_plus, d_plus, u_star, tau, tau_x, tau_y

def get_values_at_x(array, x_val):
    # get the values of variables in .face file at specific x values using linear interpolation
    x, y, u_plus, d_plus, u_star, tau, tau_x, tau_y = get_wall_param(array)
    u_plus = np.interp(np.array(x_val), x, u_plus)
    d_plus = np.interp(np.array(x_val), x, d_plus)
    u_star = np.interp(np.array(x_val), x, u_star)
    tau = np.interp(np.array(x_val), x, tau) # m²/s²
    tau_x = np.interp(np.array(x_val), x, tau_x) # m²/s²
    tau_y = np.interp(np.array(x_val), x, tau_y) # m²/s²
    return u_plus, d_plus, u_star, tau, tau_x, tau_y

def get_residu_name(filename, indice):
    # trouve le nom du residu le plus grand à la dernière iteration
    # afin de l'afficher dans les courbes des résidus maximaux
    f = open(filename,"r")
    lines = f.readlines()
    line = lines[0].strip("#").strip("\n").strip(" ").split("\t")
    line = [x for x in line if x!='']
    return line[4+indice]

def read_perf(filename,old_name, new_name):
    f = open(filename, "r")
    lines = f.readlines()
    perf = lines[0]
    perf = perf.replace(old_name, new_name)
    return perf

def get_coord_sonde(filename):
    x,y = [], []
    f = open(filename, "r")
    lines = f.readlines()
    coord = lines[1].split(" ")
    for i in range(len(coord)):
        if coord[i] == "x=":
            x.append(float(coord[i+1]))
        elif coord[i] == "y=":
            y.append(float(coord[i+1]))
    return np.array(x), np.array(y)

###############################################################################
############################# First calculations ##############################
###############################################################################
    
""" We want to compute based on the simulation setup the following :
    - Reynolds number, Re_H and Re_Dh (diametre hydraulique)
    - prediction of cf (and therefore tau, u_tau ...)
    - make sure the mesh can be built
    - y value at wall (y at which the wall loi is used)
    - y⁺ for the different numerical schemes and meshes
"""
def Reynolds(u, rho, H, mu):
    return rho*u*H/mu

cf_guess = 0.1

def cf_prediction(kappa, B, Re_H, cf_guess):
    """We assume the log law to be true in the canal
    we can solve for cf use the log law
    1/sqrt(cf/2) = (1/kappa)ln(ReHsqrt(cf/2))+B-(1/kappa)
    we use an iterative solver (non linear equation)"""
    
    tol = 1e-3
    e = 1
    while e>tol:
        loglaw = (1/kappa)*log(Re_H*sqrt(cf_guess/2))+B-(1/kappa)
        cf = 2*(1/loglaw)**2
        e = abs((cf-cf_guess)/cf_guess)
        cf_guess = cf
    return cf

def tau_and_u_tau(cf,u,rho):
    # get tau_wall and u_tau from cf (friction coef) values
    tau_w = cf*0.5*rho*u**2
    u_tau = sqrt(tau_w/rho)
    return tau_w, u_tau


def check_mesh_param(Ny,taux1,taux2):
    # check if the mesh parameters entered by the user can build a mesh
    if Ny<2:
        print("Attention maillage à refaire!:\n\nLes valeurs de Ny doivent etre supérieures à 2 !!\n")
    elif 2*int((Ny-1)/(taux1))+1<2:
        print("Attention maillage à refaire!:\n\nLes valeurs de taux1 sont trop grandes, il faut un Ny plus grand!!\n")
    elif 2*int((3*(Ny-1))/(taux1*taux2))+1<2:
        print("Attention maillage à refaire!:\n\nLes valeurs de taux2 sont trop grandes, il faut un Ny plus grand!!\n")
    else:
        print("       Ce maillage peut être construit")

def y_wall(y1,Ny):
    # get y value at which the wall velocity is computed
    # this helps us know the y+ value at which the wall law is used
    
    y = y1/(2*(Ny-1))
    y_VDF = y # frontière sup de la 1ere cellule à la paroi
    y_PolyMAC = y # frontière sup de la 1ere cellule à la paroi
    y_VEF = y/2 # hauteur du centre de la face des éléments à la paroi
    return y_VDF, y_VEF, y_PolyMAC

def y_plus(y,u_tau,mu,rho):
    return y*u_tau*rho/mu

def y_plus_prediction(y1,Ny,u_tau,mu,rho):
    y_VDF, y_VEF, y_PolyMAC = y_wall(y1,Ny)
    y_VDF_plus, y_VEF_plus, y_PolyMAC_plus = y_plus(y_VDF,u_tau,mu,rho), y_plus(y_VEF,u_tau,mu,rho), y_plus(y_PolyMAC,u_tau,mu,rho)
    return y_VDF_plus, y_VEF_plus, y_PolyMAC_plus

###############################################################################
###############################################################################
###############################################################################
