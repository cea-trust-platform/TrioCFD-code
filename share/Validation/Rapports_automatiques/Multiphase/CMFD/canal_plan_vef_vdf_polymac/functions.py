#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 23 11:03:38 2023

@author: Moncef El Moatamid
"""
from string import Template
import numpy as np
from math import *

# # paramètres fixés pour tous les calculs
    
kappa = 0.41 # loi log
B = 5.2 # loi log
y1 = 0.2 # [m] height of mesh 1st block
y2 = 0.4 # [m] height of mesh 1st block
y3 = 1 # [m] height of mesh 1st block
AR = 20 # aspect ratio
AR_polymac = 1 # trainguler avec polymac
L = 100 # [m] Longueur du domaine (canal)
H = 1 # [m] demi hauteur du canal (CL : symétrie)

description_fiche = f"""On simule un écoulement monophasique turbulent en RANS 
    dans un canal plan avec les différentes configurations:\n
    \t- VDF, VEF et PolyMAC (schémas numériques)\n
    \t- k-epsilon, k-omega et k-tau (modèles de turbulence)\n
    
    La simulation est faite sur un canal de longueur L=100m et hauteur 2H=2m
                                         symétrie
                _______________________________________________________________
                |                                                             |
    Inlet:      |             Bloc maillage n°3 : hauteur=0.6m                | Outlet:
      velocity->|-------------------------------------------------------------|-> pressure
      k       ->|             Bloc maillage n°2 : hauteur=0.2m                |-> k
      epsilon ->|-------------------------------------------------------------|-> epsilon
                |             Bloc maillage n°1 : hauteur=0.2m                |
                |_____________________________________________________________|
                                         Paroi
    
    
    il est possible de lancer le calcul sur plusieurs maillages : \n
    L'utilisateur choisi le nombre de noeuds selon y dans le bloc n°1 : \n
            liste du nombre de noeuds par maillage Ny\n
    Comme les rectangles sont coupés en 4 triangles pour VEF, \n
    il y aura 2 fois plus de mailles que la valeur entrée\n
    Attention : le rapport d'aspect AR={AR} est défini dans function.py.\n
                pour VEF, si Ny=2 le rapport d'aspect est diviser par 2\n
    il est possible de déraffiner le maillage dans les bloc n°2 et 3 avec \n
    des facteurs multiplictaifs des tailles de mailles : taux1 et taux2\n
    """


def method_name(config):
    # identifier le mot clé de la discretisation utilisée à partir du nom de la config entrée par l'utilisateur
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

###############################################################################
######################### Substitution in data file ###########################
###############################################################################

def sonde_firstpoint(y1,Ny):
    # y1 : taille 1ere bloc de maillage
    # N : nombre de noeuds selon y dans le 1er bloc
    # on veut un 1er point à 1/2 de la maille
    first_point = []
    for N in Ny:
        N=2*(N-1)+1
        dy = y1/(N-1) # taille maille
        first_point.append(dy/2) 
    return np.array(first_point)

def sonde_nbpoints(y1,Ny):
    # y1 : taille 1ere bloc de maillage
    # N : nombre de noeuds selon y dans le 1er bloc
    # on veut un point par maille dans la sonde
    nb_points = []
    for N in Ny:
        N=2*(N-1)+1
        n1 = N-1 # nombre de mailles bloc 1
        n_tot = n1*(ceil(H/y1)) # nombre de mailles tout le domaine selon y (hyp:pas d'inflation)
        nb_points.append(n_tot)
    return nb_points

# #
# # substituion for pb hydraulique
# #

def substitution_pb_hydr(params):
    # Cette fonction crée une liste de dictionnaire pour substituer les variables entrée par l'utilisateur dans le jeu de données
    # Cette fonction ne concerne que les jeux de données utilisant pb_hydraulique_turbulent
    
    
    # récupérer les paramètres entrés par l'utilisateur
    config, Ny, taux1, taux2, method_Pb_hydr, method_Pb_multi, inlet_velocity,outlet_pressure,inlet_k,inlet_epsilon,x_prof,mu,rho,y_min_prof,y_max_prof,nb_points_prof,tmax = params
    nb_mesh = len(Ny) # nombre de maillages différents du domaines (raffinements choisis)
    nb_method_Pb_hydr = len(method_Pb_hydr) # nombre de configuration avec Pb_hydrau
    dict_list = [[0 for x in range(nb_mesh)] for y in range(nb_method_Pb_hydr)] # matrice qui sera remplie avec les dictionnaires
    for i in range(nb_method_Pb_hydr):
        for j in range(nb_mesh):
            AR=20 # aspect ratio
            trianguler = ""
            N = 2*int(Ny[j]-1)+1  # nombre de noeuds bloc 1 : en vef trianguler_h coupe les cellules en 4, donc on fait x2 en y
            if method_Pb_hydr[i] == "VEFPreP1b":
                trianguler = "trianguler_H dom"
                N = Ny[j] # en vef trianguler_h coupe les cellule en 4
                if N==2: # si on a une seule maille au bloc 1 vef ne converge pas, il faut baisser l'aspect ratio
                    AR = int(AR/2)
            dict = {"Ny1" : N,
                    "Ny2" : int((N-1)/(taux1[j]))+1,
                    "Ny3" : int((3*(N-1))/(taux1[j]*taux2[j]))+1,
                   "Nx" : int(L/y1/AR*(N-1))+1,
                    "trianglemesh" : trianguler,

                    "method" : method_Pb_hydr[i],
                    "inlet_velocity" : inlet_velocity,
                    "outlet_pressure" : outlet_pressure,
                    "inlet_k" : inlet_k,
                    "inlet_epsilon" : inlet_epsilon,

                   "x_prof" : x_prof,
                   "mu" : mu,
                   "rho" : rho,
                   "y_min_prof" : y_min_prof[j],
                   "y_max_prof" : y_max_prof[j],
                   "nb_points_prof" : nb_points_prof[j],
                   "tmax" : tmax}
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

def sonde_yplus(N,AR,L,y1):
    # On veut avoir y_plus au centre de masse des éléments de maillage le long de l'axe x pour récupérer tau_wall
    # cette fonction donne les paramètres de la sonde en fonction du maillage choisi
    
    dx = L/(AR*(N-1)) # taille d'une maille
    x_min = dx/2
    x_max = L-x_min
    nb_points = AR*(N-1) # nb points sonde
    dy = y1/(N-1) # position selon y du centre de maille
    return x_min, x_max, nb_points, dy/2

def abscisses_sonde_yplus(config,Ny):
    # Cette fonction donne la liste des abscisses de la sonde de y+ (pour faire les tracé) et la valeur de y_demi
    # qui est utile pour récupérer tau_wall car y+=y_demi*u_tau/nu
    
    N = 2*int(Ny-1)+1
    AR = 20
    if "triangle" in config:
        AR = AR_polymac
        N = Ny
    x_min, x_max, nb_points, y_demi = sonde_yplus(N,AR,L,y1)
    return np.linspace(x_min, x_max, nb_points), y_demi

def get_tau_from_yplus(y_plus, y_demi, rho, mu):
    # cette fonction calcule tau_wall et u_tau à partir de y+
    
    u_tau = y_plus*mu/(rho*y_demi)
    tau = rho*u_tau**2
    return u_tau, tau

def get_value_at_x(var,x,x_tar):
    # interpolation linéaire
    
    return np.interp(np.array(x_tar), x, var)

def substitution_pb_multi(params):
    # Cette fonction crée une liste de dictionnaire pour substituer les variables entrée par l'utilisateur dans le jeu de données
    # Cette fonction ne concerne que les jeux de données utilisant pb_hydraulique_turbulent
    
    config, Ny, taux1, taux2, method_Pb_hydr, method_Pb_multi, inlet_velocity,outlet_pressure,inlet_k,inlet_epsilon,x_prof,mu,rho,y_min_prof,y_max_prof,nb_points_prof,tmax = params
    nb_mesh = len(Ny)
    nb_method_Pb_multi = len(method_Pb_multi)
    dict_list = [[0 for x in range(nb_mesh)] for y in range(nb_method_Pb_multi)]
    for i in range(nb_method_Pb_multi):
        trianguler = ""
        for j in range(nb_mesh):
            AR=20
            N = 2*int(Ny[j]-1)+1  # en vef trianguler_h coupe les cellules en 4
            comment_yplus_sonde = ""
            if "triangle" in config[i+2]: # si on veut un maillage triangulaire avec polymac
                trianguler = "trianguler_H dom"
                AR = AR_polymac # Il faut réduit le rapport d'aspect pour que ça marche
                N=Ny[j]
                comment_yplus_sonde = " # "
            x_min_yplus, x_max_yplus, nb_points_yplus, y_yplus = sonde_yplus(N,AR,L,y1)
            dict = {"Ny1" : N,
                    "Ny2" : int((N-1)/(taux1[j]))+1,
                    "Ny3" : int((3*(N-1))/(taux1[j]*taux2[j]))+1,
                   "Nx" : int(L/y1/AR*(N-1))+1,
                    "method" : method_Pb_multi[i],
                    "inlet_velocity" : inlet_velocity,
                    "outlet_pressure" : outlet_pressure,
                   "x_prof" : x_prof,
                   "mu" : mu,
                   "rho" : rho,
                   "y_min_prof" : y_min_prof[j],
                   "y_max_prof" : y_max_prof[j],
                   "nb_points_prof" : nb_points_prof[j],
                   "x_min_yplus" : x_min_yplus,
                   "x_max_yplus" : x_max_yplus,
                   "y_yplus" : y_yplus,
                   "nb_points_yplus" : nb_points_yplus,
                   "tmax" : tmax,
                   "facsec" : str(facsec),
                   "options_vdf" : options_vdf[i],
                   "diffusion_sup" : diffusion_sup[i],
                   "diffusion" : diffusion[i%2] ,
                   "diss": diss[i%2],
                   "diss_ext": diss_ext[i%2],
                   "diss_conv" : diss_conv[i%2],                                  
                   "IC_diss" : IC_diss[i%2] ,
                   "CL_diss" : CL_diss[i%2] ,
                   "CL_k" : CL_k[i%2] ,
                   "equation" : equation[i%2],
                   "nb_pas_dt_max" : nb_pas_dt_max,
                   "trianglemesh" : trianguler,
                   "comment" : comment_yplus_sonde}
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
