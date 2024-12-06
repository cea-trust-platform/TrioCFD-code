#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 23 11:03:38 2023

@author: Moncef El Moatamid
"""
from string import Template
import numpy as np
from math import *
import pandas as pd


description_fiche = f"""On simule un écoulement diphasique turbulent en RANS 
    dans une conduite verticale de section circulaire avec les différentes configurations:\n
    
    """


###############################################################################
######################### Substitution in data file ###########################
###############################################################################

# #
# # substituion 
# #
exp = "Liu" # experience
g = 9.81 # gravite
rho_l = 998.30
rho_g = 1.2
carac_essai = pd.read_csv(f"donnees/LiuBankoff.csv")

def substitution_mesh_sch(nb_mesh, param_sch):
    # this function returns a list of dictionnaries used to substitute
    # parameters of the mesh and time scheme "sch" in the data files
    filename = [f"mesh_{exp}_{i}.med" for i in range(1, nb_mesh+1)]
    tmax, facsec, max_facsec, nb_pas_dt_max, seuil_statio, solveur = param_sch
    dict_list = [0 for x in range(nb_mesh)]
    for i in range(nb_mesh):
            dict_list[i] = {
                "name_mesh" : filename[i],

                "tmax" : tmax,
                "facsec" : facsec,
                "max_facsec" : max_facsec,
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
#diff_u     = ["k_omega"]
substitution_Wlu = ['paroi_frottante_loi { }']
echelle_ou_taux = ['taux_dissipation_turbulent']
substitution_diffusion_diss = ["turbulente"]
substitution_diffusion_supplementaire_tau = [" "]
substitution_Wldiss = ["Cond_lim_omega_demi { } "]
substitution_Wlk = [' Cond_lim_k_simple_flux_nul ']

C_k = 0.002                   # HZDR model coefficient for k equation
C_epsilon = 1.0               # HZDR model coefficient for omega equation
        
#Two-phase quantities
frottement_interfacial = [" Tomiyama { contamination 2 } " ]
masse_ajoutee = [" coef_constant { } " ]
portance_interfaciale = [" Tomiyama { } " ]
dispersion_bulles = [" turbulente_burns { } " ]
beta_portance = ["1"]
beta_disp     = ["1"]
beta_wall_disp= ["1"]
beta_wall_lift= ["1"]

def substitution_turbulence(turbulence_model):
    # this function returns a dictionnary used to substitute
    # parameters of the turbulence model of two-phase flows
    source_terme = ""
    diffusion_qdm = ""
    equation_kWIT = "" 
    source_BIF = ""
    champ_BIF = ""
    sonde_BIF = ""
    sources_WIT = ""
    sonde_WIT = ""
    champ_WIT = ""
    source_k_HZDR = ""
    source_omega_HZDR = ""
    champ_k_HZDR = "" # ajouter le champ du terme source dans le post-traitement
    champ_omega_HZDR = ""  # ajouter le champ du terme source dans le post-traitement
    sonde_k_HZDR = ""
    sonde_omega_HZDR = ""
    champ_post_k_HZDR = "" 
    champ_post_omega_HZDR = ""

    # # viscosité turbulente : mots clés
    models = turbulence_model.split("_")
    if len(models)>1: # on utilise viscosit_multiple s'il y a plusieurs modeles
        diffusion_qdm += "multiple { "
    for model in models:
        if model == "k-omega" or model == "k-omega-HZDR":
            model = "k_omega"
        diffusion_qdm += f"{model} {model} {{ }} "
    if len(models)>1: # on utilise viscosit_multiple s'il y a plusieurs modeles
        diffusion_qdm += " }"
    else:
        diffusion_qdm = f"{model} {{ }} "

    # # qdm source term
    if "sato" in turbulence_model:
        source_terme = ",\n		source_BIF"
    elif "WIF" in turbulence_model or "WIT" in turbulence_model:
        source_terme = ",\n		source_BIF"

    # # k and omega source terms
    if "HZDR" in turbulence_model:
        source_k_HZDR = f",\n	Production_HZDR {{ constante_gravitation 9.81 C_k {C_k} }}"
        source_omega_HZDR = f",\n	Source_Dissipation_HZDR {{ constante_gravitation 9.81 C_k {C_k} C_epsilon {C_epsilon} }}"

    # # terme source à post traiter
    if "sato" in turbulence_model or "WIF" in turbulence_model or "WIT" in turbulence_model:
        source_BIF = """BIF 	operateur_eqn	{
									numero_source 5
									sources { refChamp { pb_champ pb vitesse } }
									}"""
        champ_BIF = "BIF elem"
        sonde_BIF = "BIF   	BIF  	periode 1e8 position_like k"
    
    if "WIT" in turbulence_model:
        sources_WIT = """prod_k_WIT 	operateur_eqn	{
                numero_source 0
                sources { refChamp { pb_champ pb k_WIT } }
                  }
          diss_k_WIT 	operateur_eqn	{
              numero_source 1
              sources { refChamp { pb_champ pb k_WIT } }
                  }
          diff_k_WIT	       operateur_eqn	{
              numero_op 0
              sources { refChamp { pb_champ pb k_WIT } }
                  }
          conv_k_WIT 	operateur_eqn	{
              numero_op 1
              sources { refChamp { pb_champ pb k_WIT } }
                  } """
        sonde_WIT = """k_WIT   	k_WIT  	periode 1e8 position_like k
        prod_k_WIT	prod_k_WIT		periode 1e8 position_like k
        diss_k_WIT	diss_k_WIT		periode 1e8 position_like k"""
        champ_WIT = """k_WIT elem
          prod_k_WIT elem
          diss_k_WIT elem"""
          
    if "HZDR" in turbulence_model: 
        champ_k_HZDR = """prod_k_hzdr 	operateur_eqn	{
              numero_source 2
              sources { refChamp { pb_champ pb k } }
                  }""" 
        champ_omega_HZDR = """source_omega_hzdr  operateur_eqn     {
              numero_source 3
              sources { refChamp { pb_champ pb omega } }
                  }"""  
                  
        sonde_k_HZDR = "prod_k_hzdr   	prod_k_hzdr  	periode 1e8 position_like k"
        sonde_omega_HZDR = "source_omega_hzdr   	source_omega_hzdr  	periode 1e8 position_like k"
        
        champ_post_k_HZDR = "prod_k_hzdr elem"
        champ_post_omega_HZDR = "source_omega_hzdr elem"

    nb_tests = len(carac_essai["SetNumber"]) # numbers of flow configurations in the experiment
    dict_list = [0 for x in range(nb_tests)]
    for index, row in carac_essai.iterrows():
        alpha = carac_essai["JG"][index]/(carac_essai["JG"][index] + carac_essai["JF"][index])
        d_b = carac_essai["Dbubble"][index]
        u_r = sqrt(alpha * d_b * (rho_l - rho_g) * g / rho_l)
        CL_kWIT = 0.1* u_r**2 # A verifier pour voir si c'est pas bete pcq cette estimation sort de ma tete !!!!!!!!!!!!!!!!!!!!
        # # equation k_WIT
        if "WIT" in turbulence_model:
            equation_kWIT = f"""energie_cinetique_turbulente_WIT
        {{
            diffusion {{ turbulente SGDH_WIT {{  }} }}
            convection {{ amont }}
            initial_conditions {{ k_WIT champ_fonc_xyz dom 1 {CL_kWIT} }}
            boundary_conditions
            {{
                wall  paroi 
                bottom	frontiere_ouverte k_ext Champ_Front_Uniforme 1 {CL_kWIT}
                top	frontiere_ouverte k_ext Champ_Front_Uniforme 1 {CL_kWIT}
                    symetrie neumann_paroi champ_front_uniforme 1 0
            }}
            sources
            {{
            Production_WIT {{ g {g} }}  ,
            Dissipation_WIT {{ constante_gravitation {g} }}
            }}
        }}"""
        dict_list[index] = {
            "diffusion_qdm" : diffusion_qdm,
            "source_BIT" : source_terme,
            "equation_kWIT" : equation_kWIT,
            "source_BIF" : source_BIF,
            "champ_BIF" : champ_BIF,
            "sonde_BIF" : sonde_BIF,
            "sources_WIT" : sources_WIT,
            "sonde_WIT" : sonde_WIT,
            "champ_WIT" : champ_WIT,
            "source_k_HZDR" : source_k_HZDR,
            "source_omega_HZDR" : source_omega_HZDR,
            "champ_k_HZDR" : champ_k_HZDR,
            "champ_omega_HZDR" : champ_omega_HZDR,
            "sonde_k_HZDR" : sonde_k_HZDR,
            "sonde_omega_HZDR" : sonde_omega_HZDR,
            "champ_post_k_HZDR" : champ_post_k_HZDR, 
            "champ_post_omega_HZDR" : champ_post_omega_HZDR
        }
    return dict_list


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
            "WLu"  : substitution_Wlu[0],
            "echelle_ou_taux" : echelle_ou_taux[0] ,
            "diffusion_diss"  : substitution_diffusion_diss[0] ,
            "diffusion_supplementaire_tau" : substitution_diffusion_supplementaire_tau[0] ,
            "WLdiss" : substitution_Wldiss[0] ,
            "WLk" : substitution_Wlk[0],
            "CI_diss": row["CL_om"] ,
            "CI_k": row["CL_k"] ,
            "u_r" : "0.1",

            "h_sonde" : str(row["h"]*0.947986) ,
            "x_sonde" : str(row["D_h"]/2*np.cos(2.5*np.pi/180)) ,
            "y_sonde" : str(row["D_h"]/2*np.sin(2.5*np.pi/180)) ,
            #"x_sonde_diag" : str(row["D_h"]/2*np.cos(25.6979*np.pi/180)) ,
            #"y_sonde_diag" : str(row["D_h"]/2*np.sin(25.6979*np.pi/180)) ,

        }
    return dict_list # as much dictionnaries as there are tests in the experimental study used

def GenerateInputFile(dir,build,substitutions_dict,name, new_name):
    with open(f"{build}/{name}.data", "r") as file: 
        filedata = Template(file.read())
    result = filedata.substitute(substitutions_dict)
    with open(dir+f"/{new_name}.data", "w") as file:
        file.write(result)

def GenerateInputFile_multidict(dir,build,dict_list,name, new_name):
    with open(f"{build}/{name}.data", "r") as file: 
        filedata = Template(file.read())
    for substitutions_dict in dict_list:
        result = filedata.substitute(substitutions_dict)
        filedata = result
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

def read_perf(filename,old_name, new_name):
    f = open(filename, "r")
    lines = f.readlines()
    perf = lines[0]
    perf = perf.replace(old_name, new_name)
    return perf
