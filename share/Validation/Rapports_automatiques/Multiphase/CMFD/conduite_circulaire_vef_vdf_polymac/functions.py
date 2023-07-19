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

###############################################################################
######################### Substitution in data file ###########################
###############################################################################

# # A mettre à jour lors de l'ajout d'une configuration (nouveau modèle ou nouveau problème ou les deux) !!!

available_config = ["VEF_k-epsilon", "VDF_k-epsilon", "VDF_k-tau", "VDF_k-omega", "PolyMAC_k-tau", "PolyMAC_k-omega"]
available_config_pb_hydr = ["VEF_k-epsilon", "VDF_k-epsilon"]
available_config_pb_multi = ["VDF_k-tau", "VDF_k-omega", "PolyMAC_k-tau", "PolyMAC_k-omega"]

def file_name(config):
    # identifier le nom du jeu de donnée à utiliser (dans src/) à partir de la config
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

"""
On utilise un dictionnaire qui définit pour chaque configuration (ex : "VEF_k-epsilon") les valeurs
des variables qu'il faut remplacer dans les jdd
"""
param_config={} 

# # initialisation :
param_config["VEF_k-epsilon"] = {}
param_config["VDF_k-epsilon"] = {}
param_config["VDF_k-omega"] = {}
param_config["VDF_k-tau"] = {}
param_config["PolyMAC_k-omega"] = {}
param_config["PolyMAC_k-tau"] = {}

# # attribution des valeurs
# method de discretisation
for config in available_config:
    if "VEF" in config:
        param_config[config]["method"] = "VEFPreP1b"
    elif "VDF" in config:
        param_config[config]["method"] = "VDF"
    else:
        param_config[config]["method"] = "PolyMAC_P0"
# Option VDF 
for config in available_config:
    if "VDF" in config:
        param_config[config]["options_vdf"] = "option_vdf { all_options }"
    else:
        param_config[config]["options_vdf"] = ""
# facsec
for config in available_config:
    param_config[config]["facsec"] = "1"
# nb_pas_dt_max 
for config in available_config:
    param_config[config]["nb_pas_dt_max"] = "1000000"
# solveur temps
for config in available_config:
    if "k-tau" in config:
        param_config[config]["solveur_temps"] = "sets"
    else:
        param_config[config]["solveur_temps"] = "ice"
# diss_conv : Condition de convergence du calcul sur le terme de dissipation turbulente
for config in available_config:
    if "k-tau" in config:
        param_config[config]["diss_conv"] = " tau 1.e-5 "
    elif "k-omega" in config:
        param_config[config]["diss_conv"] = " omega 1. "
    else:
        param_config[config]["diss_conv"] = ""
# difffusion : terme de diffusion dans l'eq de qdm
for config in available_config:
    if "k-tau" in config:
        param_config[config]["diffusion"] = "k_tau"
    elif "k-omega" in config:
        param_config[config]["diffusion"] = "k_omega"
    else:
        param_config[config]["diffusion"] = ""
# equation : équation à utiliser pour la dissipation (tau ou omega)
for config in available_config:
    if "k-tau" in config:
        param_config[config]["equation"] = "echelle_temporelle_turbulente"
    elif "k-omega" in config:
        param_config[config]["equation"] = "taux_dissipation_turbulent"
    else:
        param_config[config]["equation"] = ""
# IC_diss : Conditions initiales pour tau et omega
for config in available_config:
    if "k-tau" in config:
        param_config[config]["IC_diss"] = "2"
    elif "k-omega" in config:
        param_config[config]["IC_diss"] = str(1./2)
    else:
        param_config[config]["IC_diss"] = ""
# diss_conv : Condition de convergence du calcul sur le terme de dissipation turbulente
for config in available_config:
    if "k-tau" in config:
        param_config[config]["diss_ext"] = "tau_ext"
    elif "k-omega" in config:
        param_config[config]["diss_ext"] = "omega_ext"
    else:
        param_config[config]["diss_ext"] = ""
# diss : Mot clé terme de dissipation turbulente
for config in available_config:
    if "k-tau" in config:
        param_config[config]["diss"] = "tau"
    elif "k-omega" in config:
        param_config[config]["diss"] = "omega"
    else:
        param_config[config]["diss"] = ""
# CL_diss : Condition limite à la paroi pour l'eq de transport de la dissipation turbulente (tau et omega)
for config in available_config:
    if "k-tau" in config:
        param_config[config]["CL_diss"] = "scalaire_impose_paroi champ_front_uniforme 1 0 "
    elif "k-omega" in config:
        param_config[config]["CL_diss"] = "Cond_lim_omega_demi { }  "
    else:
        param_config[config]["CL_diss"] = ""
# diffusion_sup : terme source supp dans l'équation de transport de tau
for config in available_config:
    if "k-tau" in config:
        param_config[config]["diffusion_sup"] = " , Diffusion_supplementaire_echelle_temp_turb "
    elif "k-omega" in config:
        param_config[config]["diffusion_sup"] = " "
    else:
        param_config[config]["diffusion_sup"] = ""

def substitution(params):
    config, Ny, N_vef, inlet_velocity,outlet_pressure,inlet_k,inlet_epsilon,x_prof,mu,rho,y_min_prof,y_max_prof,y_min_profvef,y_max_profvef,nb_points_prof,tmax = params
    nb_mesh = len(Ny)
    nb_method = len(config)

    dict_list = [[0 for x in range(nb_mesh)] for y in range(nb_method)]
    for i in range(nb_method):
        theta = 1 # deg : coordonnée selon theta de la sonde
        if "VDF" in config[i]:
            theta = 90 # deg : car le maillage est 2D 
        for j in range(nb_mesh):
            AR=20
            N=Ny[j]
            x_min_yplus, x_max_yplus, nb_points_yplus, y_yplus = sonde_yplus(N,AR,L,R)
            y_min = y_min_prof
            y_max = y_max_prof
            nb_points = nb_points_prof
            file_mesh = ""
            name_mesh = f"polymac_n{N-1}_nz{int(L/R/AR*(N-1))}.med"
            if "VEF" in config[i]:
                y_min = y_min_profvef
                y_max = y_max_profvef
                nb_points = N_vef
                file_mesh = f"vef_mesh_n{N_vef[j]-1}.med"
                name_mesh = f"vef_mesh_n{N_vef[j]-1}"
            dict = {"name_mesh" : name_mesh,
                    "file_mesh" : file_mesh,
                    "Ny" : N,
                   "Nx" : int(L/R/AR*(N-1))+1,
                    "method" : param_config[config[i]]["method"],

                    "tmax" : tmax,
                    "facsec" : param_config[config[i]]["facsec"],
                    "nb_pas_dt_max" : param_config[config[i]]["nb_pas_dt_max"],

                    "mu" : mu,
                    "rho" : rho,

                    "inlet_velocity" : inlet_velocity,
                    "outlet_pressure" : outlet_pressure,
                    "inlet_k" : inlet_k,
                    "inlet_epsilon" : inlet_epsilon,

                    "diffusion_sup" : param_config[config[i]]["diffusion_sup"],
                   "diffusion" : param_config[config[i]]["diffusion"] ,
                   "diss": param_config[config[i]]["diss"],
                   "diss_ext": param_config[config[i]]["diss_ext"],
                   "diss_conv" : param_config[config[i]]["diss_conv"],                                  
                   "IC_diss" : param_config[config[i]]["IC_diss"] ,
                   "CL_diss" : param_config[config[i]]["CL_diss"] ,
                   "CL_k" : "Cond_lim_k_simple_flux_nul " ,
                   "equation" : param_config[config[i]]["equation"],
                   
                   # sonde profis (dans la section du conduit)
                   "z_prof" : x_prof,                 
                   "x_min_prof" : (R-y_min[j])*np.cos(theta*np.pi/180),
                   "y_min_prof" : (R-y_min[j])*np.sin(theta*np.pi/180),
                   "x_max_prof" : (R-y_max[j])*np.cos(theta*np.pi/180),
                   "y_max_prof" : (R-y_max[j])*np.sin(theta*np.pi/180),
                   "nb_points_prof" : nb_points[j],
                    # sonde y+ : le long de la paroi (1 point par cellule de mailllage)
                   "z_min_yplus" : x_min_yplus,
                   "z_max_yplus" : x_max_yplus,
                   "x_yplus" : (R-y_yplus)*np.cos(theta*np.pi/180),
                   "y_yplus" : (R-y_yplus)*np.sin(theta*np.pi/180),
                   "nb_points_yplus" : nb_points_yplus,
                   # sonde de pression : calcul de la perte de charge
                   "x_p" : (R/4)*np.cos(theta*np.pi/180),
                   "y_p" : (R/4)*np.sin(theta*np.pi/180),
                   "z_min_p" : 5,
                   "z_max_p" : L
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
######################### fonctions pré-traitement ############################
###############################################################################

def sonde_firstpoint(R,Ny):
    # coordonnée du 1er point de la sonde
    # N : nombre de noeuds selon y 
    # on veut un 1er point à 1/2 de la maille
    first_point = []
    for N in Ny:
        dy = R/(N-1) # taille maille
        first_point.append(dy/2) 
    return np.array(first_point)

def sonde_yplus(N,AR,L,H):
    # paramètre de la sonde de y+ le lond de la paroi
    # On veut avoir y_plus au centre de masse des éléments de maillage
    dx = L/(AR*(N-1)) # taille d'une maille
    x_min = dx/2
    x_max = L-x_min
    nb_points = int(L/R/AR*(N-1)) # nb points sonde
    dy = H/(N-1) # position selon x du centre de maille
    return x_min, x_max, nb_points, dy/2

def abscisses_sonde_yplus(config,Ny):
    # cette fonction donne une liste des coordonnées le long de la paroi de la sonde de y+
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
    # interpolation linéaire
    return np.interp(np.array(x_tar), x, var)
        
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

def get_wall_param2D(array):
    # récupérer les données du fichier Ustar.face
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

def get_values_at_x2D(array, x_val):
    # récupérer les données du fichier Ustar.face en x_val (interpolation linéaire)
    x, y, u_plus, d_plus, u_star, tau, tau_x, tau_y = get_wall_param2D(array)
    u_plus = np.interp(np.array(x_val), x, u_plus)
    d_plus = np.interp(np.array(x_val), x, d_plus)
    u_star = np.interp(np.array(x_val), x, u_star)
    tau = np.interp(np.array(x_val), x, tau) # m²/s²
    tau_x = np.interp(np.array(x_val), x, tau_x) # m²/s²
    tau_y = np.interp(np.array(x_val), x, tau_y) # m²/s²
    return u_plus, d_plus, u_star, tau, tau_x, tau_y

def get_wall_param3D(array):
    # récupérer les données du fichier Ustar.face
    array = np.transpose(array)
    x = array[0]
    y = array[1]
    z = array[2]
    u_plus = array[3]
    d_plus = array[4]
    u_star = array[5]
    tau = array[6] # m²/s²
    tau_x = array[7] # m²/s²
    tau_y = array[8] # m²/s²
    tau_z = array[9] # m²/s²
    sort_index = np.argsort(z)
    x, y, z, u_plus, d_plus, u_star, tau, tau_x, tau_y, tau_z = x[sort_index], y[sort_index], z[sort_index], u_plus[sort_index], d_plus[sort_index], u_star[sort_index], tau[sort_index], tau_x[sort_index], tau_y[sort_index], tau_z[sort_index]
    return x, y, z, u_plus, d_plus, u_star, tau, tau_x, tau_y, tau_z

def get_values_at_x3D(array, z_val):
    # récupérer les données du fichier Ustar.face en z_val (interpolation linéaire)
    x, y, z, u_plus, d_plus, u_star, tau, tau_x, tau_y, tau_z = get_wall_param3D(array)
    u_plus = np.interp(np.array(z_val), z, u_plus)
    d_plus = np.interp(np.array(z_val), z, d_plus)
    u_star = np.interp(np.array(z_val), z, u_star)
    tau = np.interp(np.array(z_val), z, tau) # m²/s²
    tau_x = np.interp(np.array(z_val), z, tau_x) # m²/s²
    tau_y = np.interp(np.array(z_val), z, tau_y) # m²/s²
    tau_z = np.interp(np.array(z_val), z, tau_z) # m²/s²
    return u_plus, d_plus, u_star, tau, tau_x, tau_y, tau_z

def get_residu_name(filename, indice):
    # trouve le nom du residu le plus grand à la dernière iteration
    f = open(filename,"r")
    lines = f.readlines()
    line = lines[0].strip("#").strip("\n").strip(" ").split("\t")
    line = [x for x in line if x!='']
    return line[4+indice]

def read_perf(filename,old_name, new_name):
    # lire les performance de calcul
    f = open(filename, "r")
    lines = f.readlines()
    perf = lines[0]
    perf = perf.replace(old_name, new_name)
    return perf

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
############################## Last calculations ##############################
###############################################################################

def facteur_frottement(Re,f_guess):
    # Cette fonction calcule le facteur de frottement f en fonction du Re pour une conduite lisse
    # Au lieu d'utiliser le diagramme de Moody (compliqué pour un calcul numérique)
    # On utilise la relation de Blasius pour 4e3<Re<1e5
    # et la formule de Von Karman et Nikuradse pour Re>1e5 (résolution itérative car relation non linéaire)
    if Re>1e5:
        tol = 1e-3 # tolerance sur l'erreur
        e = 1 # erreur
        while e>tol:
            f = 1/(2*log10(sqrt(f_guess)*Re)-0.8)**2 # formule de Von Karman et Nikuradse
            e = abs((f-f_guess)/f_guess)
            f_guess = f
    elif Re>=4e3 and Re<=1e5:
        f = 0.316*Re**(-1/4) # formule de Blasius
    else: # laminaire
        f = 64/Re
    return f
