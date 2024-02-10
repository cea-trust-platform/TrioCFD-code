# pylint: disable=line-too-long
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

kappa = 0.41  # loi log
B = 5.2  # loi log
y1 = 0.2  # [m] height of mesh 1st block
y2 = 0.4  # [m] height of mesh 1st block
y3 = 1  # [m] height of mesh 1st block
AR = 20  # aspect ratio (à modifier dans quelques fonctions aussi :( )
AR_polymac = 1  # trainguler avec polymac
L = 100  # [m] Longueur du domaine (canal)
H = 1  # [m] demi hauteur du canal (CL : symétrie)

###############################################################################
######################### Substitution in data file ###########################
###############################################################################

# # A mettre à jour lors de l'ajout d'une configuration (nouveau modèle ou nouveau problème ou les deux) !!!
vef, vdf, plm = "VEFPreP1b", "VDF", "PolyMAC_P0"
keps, komg, ktau = "k-epsilon", "k-omega", "k-tau"
pbt, pbh, pbm = "turb", "hydr", "multi"

available_config = [
    (vef, keps, pbt),
    (vdf, keps, pbt),
    (vdf, ktau, pbh),
    (vdf, ktau, pbm),
    (vdf, komg, pbh),
    (vdf, komg, pbm),
    (plm, ktau, pbh),
    (plm, ktau, pbm),
    (plm, komg, pbh),
    (plm, komg, pbm),
]
available_config_turb = [
    (vef, keps, pbt),
    (vdf, keps, pbt),
]

available_config_pb_hydr = [
    (vdf, ktau, pbh),
    (vdf, komg, pbh),
    (plm, ktau, pbh),
    (plm, komg, pbh),
]

available_config_pb_multi = [
    (vdf, ktau, pbm),
    (vdf, komg, pbm),
    (plm, ktau, pbm),
    (plm, komg, pbm),
]

"""
On utilise un dictionnaire qui définit pour chaque configuration (ex : "VEF_k-epsilon") les valeurs
des variables qu'il faut remplacer dans les jdd
"""

param_config = {
    (vef, keps, pbt): {},
    (vdf, keps, pbt): {},
    (vdf, ktau, pbh): {},
    (vdf, ktau, pbm): {},
    (vdf, komg, pbh): {},
    (vdf, komg, pbm): {},
    (plm, ktau, pbh): {},
    (plm, ktau, pbm): {},
    (plm, komg, pbh): {},
    (plm, komg, pbm): {},
}


# # attribution des valeurs
# triangle_mesh
for dis, mod, pb in available_config:
    param_config[(dis, mod, pb)]["trianglemesh"] = ""
    if dis == vef:
        param_config[(dis, mod, pb)]["trianglemesh"] = "trianguler_H dom"

# Option VDF
for dis, mod, pb in available_config:
    if dis == vdf:
        param_config[(dis, mod, pb)]["options_vdf"] = "option_vdf { all_options }"
    else:
        param_config[(dis, mod, pb)]["options_vdf"] = ""
# facsec
for dis, mod, pb in available_config:
    param_config[(dis, mod, pb)]["facsec"] = "1"
# nb_pas_dt_max
for dis, mod, pb in available_config:
    param_config[(dis, mod, pb)]["nb_pas_dt_max"] = "10000"
# solveur temps
for dis, mod, pb in available_config:
    if mod == ktau:
        param_config[(dis, mod, pb)]["solveur_temps"] = "sets"
    else:
        param_config[(dis, mod, pb)]["solveur_temps"] = "ice"
# diss_conv : Condition de convergence du calcul sur le terme de dissipation turbulente
for dis, mod, pb in available_config:
    if mod == ktau:
        param_config[(dis, mod, pb)]["diss_conv"] = " tau 1.e-5 "
    elif mod == komg:
        param_config[(dis, mod, pb)]["diss_conv"] = " omega 1. "
    else:
        param_config[(dis, mod, pb)]["diss_conv"] = ""
# difffusion : terme de diffusion dans l'eq de qdm
for dis, mod, pb in available_config:
    if mod == ktau:
        param_config[(dis, mod, pb)]["diffusion"] = "k_tau"
    elif mod == komg:
        param_config[(dis, mod, pb)]["diffusion"] = "k_omega"
    else:
        param_config[(dis, mod, pb)]["diffusion"] = ""
# equation : équation à utiliser pour la dissipation (tau ou omega)
for dis, mod, pb in available_config:
    if mod == ktau:
        param_config[(dis, mod, pb)]["equation"] = "echelle_temporelle_turbulente"
    elif mod == komg:
        param_config[(dis, mod, pb)]["equation"] = "taux_dissipation_turbulent"
    else:
        param_config[(dis, mod, pb)]["equation"] = ""
# IC_diss : Conditions initiales pour tau et omega
for dis, mod, pb in available_config:
    if mod == ktau:
        param_config[(dis, mod, pb)]["IC_diss"] = "2"
    elif mod == komg:
        param_config[(dis, mod, pb)]["IC_diss"] = str(1.0 / 2)
    else:
        param_config[(dis, mod, pb)]["IC_diss"] = ""
# diss_conv : Condition de convergence du calcul sur le terme de dissipation turbulente
for dis, mod, pb in available_config:
    if mod == ktau:
        param_config[(dis, mod, pb)]["diss_ext"] = "tau_ext"
    elif mod == komg:
        param_config[(dis, mod, pb)]["diss_ext"] = "omega_ext"
    else:
        param_config[(dis, mod, pb)]["diss_ext"] = ""
# diss : Mot clé terme de dissipation turbulente
for dis, mod, pb in available_config:
    if mod == ktau:
        param_config[(dis, mod, pb)]["diss"] = "tau"
    elif mod == komg:
        param_config[(dis, mod, pb)]["diss"] = "omega"
    else:
        param_config[(dis, mod, pb)]["diss"] = "eps"
# CL_diss : Condition limite à la paroi pour l'eq de transport de la dissipation turbulente (tau et omega)
for dis, mod, pb in available_config:
    if mod == ktau:
        param_config[(dis, mod, pb)]["CL_diss"] = "scalaire_impose_paroi champ_front_uniforme 1 0 "
    elif mod == komg:
        param_config[(dis, mod, pb)]["CL_diss"] = "Cond_lim_omega_demi { }  "
    else:
        param_config[(dis, mod, pb)]["CL_diss"] = ""
# diffusion_sup : terme source supp dans l'équation de transport de tau
for dis, mod, pb in available_config:
    if mod == ktau:
        param_config[(dis, mod, pb)]["diffusion_sup"] = " , Diffusion_supplementaire_echelle_temp_turb "
    elif mod == komg:
        param_config[(dis, mod, pb)]["diffusion_sup"] = " "
    else:
        param_config[(dis, mod, pb)]["diffusion_sup"] = ""
## !! Only for pb hydraulique !!
for dis, mod, pb in available_config:
    if mod == keps:
        param_config[(dis, mod, pb)]["modele_turb"] = "k_epsilon"
        param_config[(dis, mod, pb)]["equation_hydr"] = "transport_k_epsilon"
        param_config[(dis, mod, pb)]["turb"] = "k_eps"
        param_config[(dis, mod, pb)]["CL_turb"] = "frontiere_ouverte_K_eps_impose"
        param_config[(dis, mod, pb)]["turb_ext"] = "k_eps_ext"
        param_config[(dis, mod, pb)]["sources"] = "sources			{ source_transport_k_eps { C1_eps 1.44 C2_eps 1.92 } }"
    elif mod == komg:
        param_config[(dis, mod, pb)]["modele_turb"] = "k_omega"
        param_config[(dis, mod, pb)]["equation_hydr"] = "transport_k_omega"
        param_config[(dis, mod, pb)]["turb"] = "k_omega"
        param_config[(dis, mod, pb)]["CL_turb"] = "frontiere_ouverte_K_omega_impose"
        param_config[(dis, mod, pb)]["turb_ext"] = "k_omega_ext"
        param_config[(dis, mod, pb)]["sources"] = ""
    else:
        param_config[(dis, mod, pb)]["modele_turb"] = ""
        param_config[(dis, mod, pb)]["equation_hydr"] = ""
        param_config[(dis, mod, pb)]["turb"] = ""
        param_config[(dis, mod, pb)]["CL_turb"] = ""
        param_config[(dis, mod, pb)]["turb_ext"] = ""
        param_config[(dis, mod, pb)]["sources"] = ""


# #
# # substituion for pb hydraulique and pb multiphase
# #


def substitution(params):
    # Cette fonction crée une liste de dictionnaire pour substituer les variables entrée par l'utilisateur dans le jeu de données

    # récupérer les paramètres entrés par l'utilisateur
    (
        configs,
        Ny,
        taux1,
        taux2,
        inlet_velocity,
        outlet_pressure,
        inlet_k,
        inlet_epsilon,
        x_prof,
        mu,
        rho,
        y_min_prof,
        y_max_prof,
        nb_points_prof,
        tmax,
    ) = params
    nb_mesh = len(Ny)  # nombre de maillages différents du domaines (raffinements choisis)
    nb_method = len(configs)  # nombre de configuration avec Pb_hydrau
    dict_list = [[0 for x in range(nb_mesh)] for y in range(nb_method)]  # matrice qui sera remplie avec les dictionnaires
    for i, (dis, mod, pb) in enumerate(configs):
        for j in range(nb_mesh):
            AR = 20  # aspect ratio
            N = 2 * int(Ny[j] - 1) + 1  # nombre de noeuds bloc 1 : en vef trianguler_h coupe les cellules en 4, donc on fait x2 en y
            if dis == vef:
                N = Ny[j]  # en vef trianguler_h coupe les cellule en 4
                if N == 2:  # si on a une seule maille au bloc 1 vef ne converge pas, il faut baisser l'aspect ratio
                    AR = int(AR / 2)
            x_min_yplus, x_max_yplus, nb_points_yplus, y_yplus = sonde_yplus(N, AR, L, y1)

            dict = {
                "Ny1": N,
                "Ny2": int((N - 1) / (taux1[j])) + 1,
                "Ny3": int((3 * (N - 1)) / (taux1[j] * taux2[j])) + 1,
                "Nx": int(L / y1 / AR * (N - 1)) + 1,
                "trianglemesh": param_config[(dis, mod, pb)]["trianglemesh"],
                "method": dis,
                "inlet_velocity": inlet_velocity,
                "outlet_pressure": outlet_pressure,
                "inlet_k": inlet_k,
                "inlet_diss": CL_diss((dis, mod, pb), inlet_k, inlet_epsilon),
                "x_prof": x_prof,
                "mu": mu,
                "rho": rho,
                "y_min_prof": y_min_prof[j],
                "y_max_prof": y_max_prof[j],
                "nb_points_prof": nb_points_prof[j],
                "tmax": tmax,
                "x_min_yplus": x_min_yplus,
                "x_max_yplus": x_max_yplus,
                "y_yplus": y_yplus,
                "nb_points_yplus": nb_points_yplus,
                "facsec": param_config[(dis, mod, pb)]["facsec"],
                "solveur_temps": param_config[(dis, mod, pb)]["solveur_temps"],
                "options_vdf": param_config[(dis, mod, pb)]["options_vdf"],
                "diffusion_sup": param_config[(dis, mod, pb)]["diffusion_sup"],
                "diffusion": param_config[(dis, mod, pb)]["diffusion"],
                "diss": param_config[(dis, mod, pb)]["diss"],
                "diss_ext": param_config[(dis, mod, pb)]["diss_ext"],
                "diss_conv": param_config[(dis, mod, pb)]["diss_conv"],
                "IC_diss": param_config[(dis, mod, pb)]["IC_diss"],
                "CL_diss": param_config[(dis, mod, pb)]["CL_diss"],
                "CL_k": "Cond_lim_k_simple_flux_nul ",
                "equation": param_config[(dis, mod, pb)]["equation"],
                "nb_pas_dt_max": param_config[(dis, mod, pb)]["nb_pas_dt_max"],
                # bloc turbulence pb hydr
                "modele_turb": param_config[(dis, mod, pb)]["modele_turb"],
                "equation_hydr": param_config[(dis, mod, pb)]["equation_hydr"],
                "turb": param_config[(dis, mod, pb)]["turb"],
                "CL_turb": param_config[(dis, mod, pb)]["CL_turb"],
                "turb_ext": param_config[(dis, mod, pb)]["turb_ext"],
                "sources": param_config[(dis, mod, pb)]["sources"],
            }
            dict_list[i][j] = dict
    return dict_list


###############################################################################
######################### fonctions pré-traitement ############################
###############################################################################


def CL_diss(config_name, k, eps):
    # donne la condition limite en omega et tau à partir de k et epsilon
    diss = eps
    if "omega" in config_name:
        diss = eps / 0.09 / k
    elif "tau" in config_name:
        diss = k * 0.09 / eps
    return diss


def sonde_firstpoint(y1, Ny):
    # y1 : taille 1ere bloc de maillage
    # N : nombre de noeuds selon y dans le 1er bloc
    # on veut un 1er point à 1/2 de la maille
    first_point = []
    for N in Ny:
        N = 2 * (N - 1) + 1
        dy = y1 / (N - 1)  # taille maille
        first_point.append(dy / 2)
    return np.array(first_point)


def sonde_nbpoints(y1, Ny):
    # y1 : taille 1ere bloc de maillage
    # N : nombre de noeuds selon y dans le 1er bloc
    # on veut un point par maille dans la sonde
    nb_points = []
    for N in Ny:
        N = 2 * (N - 1) + 1
        n1 = N - 1  # nombre de mailles bloc 1
        n_tot = n1 * (ceil(H / y1))  # nombre de mailles tout le domaine selon y (hyp:pas d'inflation)
        nb_points.append(n_tot)
    return nb_points


def sonde_yplus(N, AR, L, y1):
    # On veut avoir y_plus au centre de masse des éléments de maillage le long de l'axe x pour récupérer tau_wall
    # cette fonction donne les paramètres de la sonde en fonction du maillage choisi

    dx = L / (AR * (N - 1))  # taille d'une maille
    x_min = dx / 2
    x_max = L - x_min
    nb_points = AR * (N - 1)  # nb points sonde
    dy = y1 / (N - 1)  # position selon y du centre de maille
    return x_min, x_max, nb_points, dy / 2


def abscisses_sonde_yplus(config, Ny):
    # Cette fonction donne la liste des abscisses de la sonde de y+ (pour faire les tracé) et la valeur de y_demi
    # qui est utile pour récupérer tau_wall car y+=y_demi*u_tau/nu

    N = 2 * int(Ny - 1) + 1
    AR = 20
    if "triangle" in config:
        AR = AR_polymac
        N = Ny
    x_min, x_max, nb_points, y_demi = sonde_yplus(N, AR, L, y1)
    return np.linspace(x_min, x_max, nb_points), y_demi


def get_tau_from_yplus(y_plus, y_demi, rho, mu):
    # cette fonction calcule tau_wall et u_tau à partir de y+

    u_tau = y_plus * mu / (rho * y_demi)
    tau = rho * u_tau**2
    return u_tau, tau


def get_value_at_x(var, x, x_tar):
    # interpolation linéaire

    return np.interp(np.array(x_tar), x, var)


###############################################################################
################################# Read files ##################################
###############################################################################


def read_facefile(file):
    """read data in .face file at différent time values and stores the data
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
    tau = array[5]  # m²/s²
    tau_x = array[6]  # m²/s²
    tau_y = array[7]  # m²/s²
    sort_index = np.argsort(x)
    x, y, u_plus, d_plus, u_star, tau, tau_x, tau_y = (
        x[sort_index],
        y[sort_index],
        u_plus[sort_index],
        d_plus[sort_index],
        u_star[sort_index],
        tau[sort_index],
        tau_x[sort_index],
        tau_y[sort_index],
    )
    return x, y, u_plus, d_plus, u_star, tau, tau_x, tau_y


def get_values_at_x(array, x_val):
    # get the values of variables in .face file at specific x values using linear interpolation
    x, y, u_plus, d_plus, u_star, tau, tau_x, tau_y = get_wall_param(array)
    u_plus = np.interp(np.array(x_val), x, u_plus)
    d_plus = np.interp(np.array(x_val), x, d_plus)
    u_star = np.interp(np.array(x_val), x, u_star)
    tau = np.interp(np.array(x_val), x, tau)  # m²/s²
    tau_x = np.interp(np.array(x_val), x, tau_x)  # m²/s²
    tau_y = np.interp(np.array(x_val), x, tau_y)  # m²/s²
    return u_plus, d_plus, u_star, tau, tau_x, tau_y


def get_residu_name(filename, indice):
    # trouve le nom du residu le plus grand à la dernière iteration
    # afin de l'afficher dans les courbes des résidus maximaux
    f = open(filename, "r")
    lines = f.readlines()
    line = lines[0].strip("#").strip("\n").strip(" ").split("\t")
    line = [x for x in line if x != ""]
    return line[4 + indice]


def read_perf(filename, old_name, new_name):
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
    return rho * u * H / mu


cf_guess = 0.1


def cf_prediction(kappa, B, Re_H, cf_guess):
    """We assume the log law to be true in the canal
    we can solve for cf use the log law
    1/sqrt(cf/2) = (1/kappa)ln(ReHsqrt(cf/2))+B-(1/kappa)
    we use an iterative solver (non linear equation)"""

    tol = 1e-3
    e = 1
    while e > tol:
        loglaw = (1 / kappa) * log(Re_H * sqrt(cf_guess / 2)) + B - (1 / kappa)
        cf = 2 * (1 / loglaw) ** 2
        e = abs((cf - cf_guess) / cf_guess)
        cf_guess = cf
    return cf


def tau_and_u_tau(cf, u, rho):
    # get tau_wall and u_tau from cf (friction coef) values
    tau_w = cf * 0.5 * rho * u**2
    u_tau = sqrt(tau_w / rho)
    return tau_w, u_tau


def check_mesh_param(Ny, taux1, taux2):
    # check if the mesh parameters entered by the user can build a mesh
    if Ny < 2:
        print("Attention maillage à refaire!:\n\nLes valeurs de Ny doivent etre supérieures à 2 !!\n")
    elif 2 * int((Ny - 1) / (taux1)) + 1 < 2:
        print("Attention maillage à refaire!:\n\nLes valeurs de taux1 sont trop grandes, il faut un Ny plus grand!!\n")
    elif 2 * int((3 * (Ny - 1)) / (taux1 * taux2)) + 1 < 2:
        print("Attention maillage à refaire!:\n\nLes valeurs de taux2 sont trop grandes, il faut un Ny plus grand!!\n")
    else:
        print("       Ce maillage peut être construit")


def y_wall(y1, Ny):
    # get y value at which the wall velocity is computed
    # this helps us know the y+ value at which the wall law is used

    y = y1 / (2 * (Ny - 1))
    y_VDF = y  # frontière sup de la 1ere cellule à la paroi
    y_PolyMAC = y  # frontière sup de la 1ere cellule à la paroi
    y_VEF = y / 2  # hauteur du centre de la face des éléments à la paroi
    return y_VDF, y_VEF, y_PolyMAC


def y_plus(y, u_tau, mu, rho):
    return y * u_tau * rho / mu


def y_plus_prediction(y1, Ny, u_tau, mu, rho):
    y_VDF, y_VEF, y_PolyMAC = y_wall(y1, Ny)
    y_VDF_plus, y_VEF_plus, y_PolyMAC_plus = y_plus(y_VDF, u_tau, mu, rho), y_plus(y_VEF, u_tau, mu, rho), y_plus(y_PolyMAC, u_tau, mu, rho)
    return y_VDF_plus, y_VEF_plus, y_PolyMAC_plus


###############################################################################
############################## Last calculations ##############################
###############################################################################


def facteur_frottement(Re, f_guess):
    # Cette fonction calcule le facteur de frottement f en fonction du Re pour une conduite lisse
    # Au lieu d'utiliser le diagramme de Moody (compliqué pour un calcul numérique)
    # On utilise la relation de Blasius pour 4e3<Re<1e5
    # et la formule de Von Karman et Nikuradse pour Re>1e5 (résolution itérative car relation non linéaire)
    if Re > 1e5:
        tol = 1e-3  # tolerance sur l'erreur
        e = 1  # erreur
        while e > tol:
            f = 1 / (2 * log10(sqrt(f_guess) * Re) - 0.8) ** 2  # formule de Von Karman et Nikuradse
            e = abs((f - f_guess) / f_guess)
            f_guess = f
    elif Re >= 4e3 and Re <= 1e5:
        f = 0.316 * Re ** (-1 / 4)  # formule de Blasius
    else:  # laminaire
        f = 64 / Re
    return f


###############################################################################
###############################################################################
###############################################################################
