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


description_fiche = f"""On simule un écoulement diphasique turbulent en RANS 
    dans une conduite verticale de section circulaire avec les différentes configurations:\n
    \t- VDF, VEF et PolyMAC (schémas numériques)\n
    \t- k-epsilon, k-omega et k-tau (modèles de turbulence)\n
    """


###############################################################################
######################### Substitution in data file ###########################
###############################################################################

# #
# # substituion 
# #
exp = "Hibiki" # experience

def substitution_mesh_sch(nb_mesh, param_sch):
    # this function returns a list of dictionnaries used to substitute
    # parameters of the mesh and time scheme "sch" in the data files
    filename = [f"mesh_{exp}_{i}.med" for i in range(1, nb_mesh+1)]
    tmax, facsec, nb_pas_dt_max, seuil_statio, solveur = param_sch
    dict_list = [0 for x in range(nb_mesh)]
    for i in range(nb_mesh):
            dict_list[i] = {
                "name_mesh" : filename[i],

                "tmax" : tmax,
                "facsec" : facsec[i],
                "nb_pas_dt_max" : nb_pas_dt_max,
                "seuil_statio" : seuil_statio,
                "solveur" : solveur,
            }
    return dict_list # as much dictionaries as there are meshes (requested by the user of the notebook)


#Physical quantities
interface = "interface_eau_air interface_sigma_constant  { tension_superficielle 0.0728 }"
air_properties =   " gaz_air    Fluide_Incompressible { mu champ_uniforme 1 1.85e-5  rho champ_uniforme 1   1.2  lambda Champ_Uniforme  1 0.023 Cp Champ_Uniforme 1 1006   beta_th Champ_Uniforme 1 0 } "
water_properties = " liquide_eau Fluide_Incompressible { mu champ_uniforme 1 1.002e-3 rho champ_uniforme 1 998.30 lambda Champ_Uniforme  1 0.604 Cp Champ_Uniforme 1 75.366 beta_th Champ_Uniforme 1 0 } "

#Turbulent quantities
tau_ou_omega = [ "omega" ]
diff_u     = ["k_omega"]
substitution_Wlu = ['paroi_frottante_loi { }']
echelle_ou_taux = ['taux_dissipation_turbulent']
substitution_diffusion_diss = ["turbulente"]
substitution_diffusion_supplementaire_tau = [" "]
substitution_Wldiss = ["Cond_lim_omega_demi { } "]
substitution_Wlk = [' Cond_lim_k_simple_flux_nul ']
        
#Two-phase quantities
frottement_interfacial = [" Tomiyama { contamination 2 } " ]
masse_ajoutee = [" coef_constant { } " ]
portance_interfaciale = [" Sugrue { } " ]
dispersion_bulles = [" turbulente_burns { } " ]
beta_portance = ["1"]
beta_disp     = ["1"]
beta_wall_disp= ["1"]
beta_wall_lift= ["1"]

def substitution_pb(carac_essai):
    # this function returns a list of dictionnaries used to substitute
    # parameters of the problem in trust (pb_multiphase)
    nb_tests = len(carac_essai["SetNumber"]) # numbers of flow configurations in the experiment
    dict_list = [0 for x in range(nb_tests)]
    for index, row in carac_essai.iterrows():
        dict_list[index] = {
            "carrying_phase" :  water_properties , 
            "dispersed_phase" : air_properties,
            "interface" : interface,
                                          
            "diametre_bulles" : str(row["Dbubble"]),
            "u_0"  : row["u_0"] ,
            "grav" : row["g"] ,
            "alpha_l0" : row["alpha_l0"] ,
            "alpha_v0" : row["alpha_g0"] ,
            "frottement_interfacial" : frottement_interfacial[0],
            "masse_ajoutee" : masse_ajoutee[0],
            "portance_interfaciale" : portance_interfaciale[0],
            "dispersion_bulles" : dispersion_bulles[0],
            "beta_portance" : beta_portance[0],
            "beta_disp" : beta_disp[0],
            "beta_wall_disp" : beta_wall_disp[0],
            "beta_wall_lift" : beta_wall_lift[0],
            
            "tau_ou_omega" : tau_ou_omega[0] ,
            "tau_ext_ou_omega_ext" : f'{tau_ou_omega[0]}_ext' ,
            "diff_u" : diff_u[0] ,
            "WLu"  : substitution_Wlu[0],
            "echelle_ou_taux" : echelle_ou_taux[0] ,
            "diffusion_diss"  : substitution_diffusion_diss[0] ,
            "diffusion_supplementaire_tau" : substitution_diffusion_supplementaire_tau[0] ,
            "WLdiss" : substitution_Wldiss[0] ,
            "WLk" : substitution_Wlk[0],
            "CI_diss": row["CL_om"] ,
            "CI_k": row["CL_k"] ,

            "h_sonde" : str(row["h"]*0.947986) ,
            "x_sonde" : str(row["D_h"]/2*np.cos(5.6979*np.pi/180)) ,
            "y_sonde" : str(row["D_h"]/2*np.sin(5.6979*np.pi/180)) ,
            "x_sonde_diag" : str(row["D_h"]/2*np.cos(25.6979*np.pi/180)) ,
            "y_sonde_diag" : str(row["D_h"]/2*np.sin(25.6979*np.pi/180)) ,

        }
    return dict_list # as much dictionnaries as there are tests in the experimental study used

def GenerateInputFile(dir,build,substitutions_dict,name, new_name):
    with open(f"{build}/{name}.data", "r") as file: 
        filedata = Template(file.read())
    result = filedata.substitute(substitutions_dict)
    with open(dir+f"/{new_name}.data", "w") as file:
        file.write(result)


def poly(x, coeffs) :
    rep = x-x
    for i in range(len(coeffs)):
        rep += coeffs[i] * x**(len(coeffs)-1-i)
    return rep

nfitd = 4

def str_r_loc_d_loc(rl, dl):
    r_l_l = []
    d_l_l = []
    for i in range(len(rl)):
        if (isnan(dl[i])==False) :
            r_l_l += [rl[i]]
            d_l_l += [dl[i]]
            
    tab_rmin_max = [0,0]
    tab_dmin_max = [0,0]
            
    tab_polyfitd = np.polyfit(r_l_l, d_l_l, nfitd)
    tab_rmin_max[0] = (r_l_l[0] +r_l_l[1] )/2.
    tab_rmin_max[1] = (r_l_l[-1]+r_l_l[-2])/2.
    tab_dmin_max[0] = (d_l_l[0] +d_l_l[1] )/2.
    tab_dmin_max[1] = (d_l_l[-1]+d_l_l[-2])/2.

    str_diam = "0."

    #First we take care of what happens above and under the highest values
    str_diam+= f"+((X*X+Y*Y)<({tab_rmin_max[0]}*{tab_rmin_max[0]}))"
    str_diam+= f"*{tab_dmin_max[0]}"
    str_diam+= f"+((X*X+Y*Y)]({tab_rmin_max[1]}*{tab_rmin_max[1]}))"
    str_diam+= f"*{tab_dmin_max[1]}"

    #Then we take care of the middle
    str_loc = "0."
    for i in range(nfitd+1) :
        str_loc += f"+({tab_polyfitd[i]})"
        for j in range(nfitd-i) :
            str_loc+= "*sqrt(x*x+y*y)"
    str_diam += f"+((X*X+Y*Y)]({tab_rmin_max[0]}*{tab_rmin_max[0]}))" # put to zero outside the right interval
    str_diam += f"*((X*X+Y*Y)<({tab_rmin_max[1]}*{tab_rmin_max[1]}))" 
    str_diam += f"*({str_loc})"          
    str_diam = f"1.e-3*({str_diam})" # mm => m
    return str_diam

def loadText(data, index_column=0, nb_column=-1, transpose=True, dtype="float", skiprows=0):
        
    if nb_column == -1:
        nb = None
    else:
        nb = index_column + nb_column

    try:
        if transpose:
            matrix = np.loadtxt(f'build/{data}', dtype=dtype, skiprows=skiprows).T[index_column:nb]
        else:
            matrix = np.loadtxt(data, dtype=dtype, skiprows=skiprows)[index_column:nb]
    except:
        matrix = np.loadtxt(data, dtype=dtype, skiprows=skiprows)        
        
    return matrix


########!!!!!!!!!!!!le reste n'est Pas utile pour cette fiche!!!!!!!!!!!!!!!!!
        
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
