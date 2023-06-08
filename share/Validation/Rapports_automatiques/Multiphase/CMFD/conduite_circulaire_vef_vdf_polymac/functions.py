#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 23 11:03:38 2023

@author: me274984
"""
from string import Template
import numpy as np
from math import *

description_fiche = """On simule un écoulement monophasique turbulent en RANS 
    dans une conduite circulaire (tube) avec les différentes configurations:\n
    \t- VDF, VEF et PolyMAC (schémas numériques)\n
    \t- k-epsilon, k-omega et k-tau (modèles de turbulence)\n
    La simulation est faite sur un tube de longueur L=100m et rayon R=2m\n
    Les maillages utilisées avec VEF et Polymac sont générés avec salomé (voir script python dans src/)\n
    Pour VDF, on utilise un maillages 2D avec une condition limite d'axisymétrie en x=0\n
    
                                            paroi
                _______________________________________________________________
                |                                                             |
    Inlet:      |                                                             | Outlet:
      velocity->|                                                             |-> pressure
      k       ->|                                                             |-> k
      epsilon ->|                                                             |-> epsilon
                |                                                             |
                |_____________________________________________________________|
                                         axisymétrie
    """


# # paramètres fixés pour tous les calculs
    
kappa = 0.41
B = 5.2
R = 2 # rayon conduit [m]
AR = 10 # aspect ratio
AR_polymac = 1 # trainguler avec polymac
L = 100 # [m] Longueur du domaine (canal)


def method_name(config):
    # identifier le mot clé de la discretisation utilisée à partir du nom de la config
    method_Pb_hydr = []
    method_Pb_multi = []
    for method in config:
        if "k-epsilon" in method:
            name = "VDF"
            if "VEF" in method:
                name = "VEFPreP1b"
            method_Pb_hydr.append(name)
        else:
            name = "VDF"
            if "PolyMAC" in method:
                name = "PolyMAC_P0"
            method_Pb_multi.append(name)
    return method_Pb_hydr, method_Pb_multi

def file_name(config):
    names = []
    for method in config:
        name = ""
        if "VDF" in method:
            name = "vdf_axi_pb_multi"
            if "k-epsilon" in method:
                name = "vdf_axi_pb_hydr"
        elif "PolyMAC" in method:
            name = "polymac_pb_multi"
        else:
            name = "vef_pb_hydr"
        names.append(name)
    return names

###############################################################################
######################### Substitution in data file ###########################
###############################################################################

def sonde_firstpoint(R,Ny):
    # N : nombre de noeuds selon y 
    # on veut un 1er point à 1/2 de la maille
    first_point = []
    for N in Ny:
        dy = R/(N-1) # taille maille
        first_point.append(dy/2) 
    return np.array(first_point)

# #
# # substituion for pb hydraulique
# #

def substitution_pb_hydr(params):
    # récupérer les paramètres entrés par l'utilisateur
    config, Ny, N_vef, method_Pb_hydr, method_Pb_multi, inlet_velocity,outlet_pressure,inlet_k,inlet_epsilon,x_prof,mu,rho,y_min_prof,y_max_prof,nb_points_prof,tmax = params
    nb_mesh = len(Ny) # nombre de maillages différents du domaines (raffinements choisis)
    nb_method_Pb_hydr = len(method_Pb_hydr) # nombre de configuration avec Pb_hydrau
    dict_list = [[0 for x in range(nb_mesh)] for y in range(nb_method_Pb_hydr)] # matrice qui sera remplie avec les dictionnaires
    for i in range(nb_method_Pb_hydr):
        theta = 5 # deg
        if "VDF" in config[i]:
            theta = 90 # deg
        for j in range(nb_mesh):
            AR=20 # aspect ratio
            N = Ny[j]  
            x_min_yplus, x_max_yplus, nb_points_yplus, y_yplus = sonde_yplus(N,AR,L,R)
            dict = {"Ny" : N,
                   "Nx" : int(L/R/AR*(N-1))+1,
                   "file_mesh" : f"vef_mesh_t45_n{N_vef[j]-1}.med",
                    "name_mesh" : f"vef_mesh_t45_n{N_vef[j]-1}",

                    "method" : method_Pb_hydr[i],
                    "inlet_velocity" : inlet_velocity,
                    "outlet_pressure" : outlet_pressure,
                    "inlet_k" : inlet_k,
                    "inlet_epsilon" : inlet_epsilon,

                   "mu" : mu,
                   "rho" : rho,
                   "tmax" : tmax,
                   
                   "z_prof" : x_prof,                 
                   "x_min_prof" : (R-y_min_prof[j])*np.cos(theta*np.pi/180),
                   "y_min_prof" : (R-y_min_prof[j])*np.sin(theta*np.pi/180),
                   "x_max_prof" : (R-y_max_prof[j])*np.cos(theta*np.pi/180),
                   "y_max_prof" : (R-y_max_prof[j])*np.sin(theta*np.pi/180),
                   "nb_points_prof" : nb_points_prof[j],

                   "z_min_yplus" : x_min_yplus,
                   "z_max_yplus" : x_max_yplus,
                   "x_yplus" : (R-y_yplus)*np.cos(theta*np.pi/180),
                   "y_yplus" : (R-y_yplus)*np.sin(theta*np.pi/180),
                   "nb_points_yplus" : nb_points_yplus}
            dict_list[i][j] = dict
    return dict_list


# #
# # substituion for pb multiphase
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
IC_diss =      ["2"]+[str(1./2)]
# Condition limite à la paroi pour l'eq de transport de la dissipation turbulente (tau et omega)
CL_diss = ["scalaire_impose_paroi champ_front_uniforme 1 0 ", "Cond_lim_omega_demi { }  "]
# Condition limite à la paroi pour l'eq de transport k
CL_k = ["Cond_lim_k_simple_flux_nul ", "Cond_lim_k_simple_flux_nul  "]
# terme source supp dans l'équation de transport de tau (dissipation turbu)
diffusion_sup = [" "] *2+[" , Diffusion_supplementaire_echelle_temp_turb "]+[" "]+[" , Diffusion_supplementaire_echelle_temp_turb "]
# Autres param de simu
facsec     = 1
nb_pas_dt_max = "1000000"
options_vdf = ["option_vdf { all_options }", "option_vdf { all_options }", "", "", ""]

params_pb_multi = [equation, diffusion, diss, diss_conv, IC_diss, CL_k, diffusion_sup, facsec, nb_pas_dt_max]

def sonde_yplus(N,AR,L,H):
    # On veut avoir y_plus au centre de masse des éléments de maillage
    dx = L/(AR*(N-1)) # taille d'une maille
    x_min = dx/2
    x_max = L-x_min
    nb_points = int(L/R/AR*(N-1)) # nb points sonde
    dy = H/(N-1) # position selon x du centre de maille
    return x_min, x_max, nb_points, dy/2

def abscisses_sonde_yplus(config,Ny):
    N = Ny
    AR = 20
    if "triangle" in config:
        AR = AR_polymac
        N = Ny
    x_min, x_max, nb_points, y_demi = sonde_yplus(N,AR,L,R)
    return np.linspace(x_min, x_max, nb_points), y_demi

def get_tau_from_yplus(y_plus, y_demi, rho, mu):
    u_tau = y_plus*mu/(rho*y_demi)
    tau = rho*u_tau**2
    return u_tau, tau

def get_value_at_x(var,x,x_tar):
    return np.interp(np.array(x_tar), x, var)

def substitution_pb_multi(params):
    config, Ny, N_vef, method_Pb_hydr, method_Pb_multi, inlet_velocity,outlet_pressure,inlet_k,inlet_epsilon,x_prof,mu,rho,y_min_prof,y_max_prof,nb_points_prof,tmax = params
    nb_mesh = len(Ny)
    nb_method_Pb_multi = len(method_Pb_multi)
    nb_method_Pb_hydr = len(method_Pb_hydr)
    dict_list = [[0 for x in range(nb_mesh)] for y in range(nb_method_Pb_multi)]
    for i in range(nb_method_Pb_multi):
        theta = 1 # deg
        if "VDF" in config[i+nb_method_Pb_hydr]:
            theta = 90 # deg
        for j in range(nb_mesh):
            AR=20
            N=Ny[j]
            x_min_yplus, x_max_yplus, nb_points_yplus, y_yplus = sonde_yplus(N,AR,L,R)
            dict = {"name_mesh" : f"polymac_n{N-1}_nz{int(L/R/AR*(N-1))}.med",
                    "Ny" : N,
                   "Nx" : int(L/R/AR*(N-1))+1,
                    "method" : method_Pb_multi[i],

                    "tmax" : tmax,
                    "facsec" : str(facsec),
                    "nb_pas_dt_max" : nb_pas_dt_max,

                    "mu" : mu,
                    "rho" : rho,

                    "inlet_velocity" : inlet_velocity,
                    "outlet_pressure" : outlet_pressure,

                    "diffusion_sup" : diffusion_sup[i],
                   "diffusion" : diffusion[i%2] ,
                   "diss": diss[i%2],
                   "diss_ext": diss_ext[i%2],
                   "diss_conv" : diss_conv[i%2],                                  
                   "IC_diss" : IC_diss[i%2] ,
                   "CL_diss" : CL_diss[i%2] ,
                   "CL_k" : CL_k[i%2] ,
                   "equation" : equation[i%2],
                   
                    "z_prof" : x_prof,                 
                   "x_min_prof" : (R-y_min_prof[j])*np.cos(theta*np.pi/180),
                   "y_min_prof" : (R-y_min_prof[j])*np.sin(theta*np.pi/180),
                   "x_max_prof" : (R-y_max_prof[j])*np.cos(theta*np.pi/180),
                   "y_max_prof" : (R-y_max_prof[j])*np.sin(theta*np.pi/180),
                   "nb_points_prof" : nb_points_prof[j],

                   "z_min_yplus" : x_min_yplus,
                   "z_max_yplus" : x_max_yplus,
                   "x_yplus" : (R-y_yplus)*np.cos(theta*np.pi/180),
                   "y_yplus" : (R-y_yplus)*np.sin(theta*np.pi/180),
                   "nb_points_yplus" : nb_points_yplus,
                   }
            dict_list[i][j] = dict
    return dict_list


def GenerateInputFile(dir,build,substitutions_dict, name):
    with open(f"{build}/{name}.data", "r") as file: 
        filedata = Template(file.read())
    result = filedata.substitute(substitutions_dict)
    with open(dir+f"/{name}.data", "w") as file:
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
    f = open(filename,"r")
    lines = f.readlines()
    line = lines[0].strip("#").strip("\n").strip(" ").split("\t")
    line = [x for x in line if x!='']
    return line[4+indice]

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
    tau_w = cf*0.5*rho*u**2
    u_tau = sqrt(tau_w/rho)
    return tau_w, u_tau



def y_wall(R,Ny):
    y = R/(Ny-1)
    y_VDF = y
    y_PolyMAC = y
    y_VEF = y/2
    return y_VDF, y_VEF, y_PolyMAC

def y_plus(y,u_tau,mu,rho):
    return y*u_tau*rho/mu

def y_plus_prediction(y1,Ny,u_tau,mu,rho):
    y_VDF, y_VEF, y_PolyMAC = y_wall(R,Ny)
    y_VDF_plus, y_VEF_plus, y_PolyMAC_plus = y_plus(y_VDF,u_tau,mu,rho), y_plus(y_VEF,u_tau,mu,rho), y_plus(y_PolyMAC,u_tau,mu,rho)
    return y_VDF_plus, y_VEF_plus, y_PolyMAC_plus

###############################################################################
############################# First calculations ##############################
###############################################################################
