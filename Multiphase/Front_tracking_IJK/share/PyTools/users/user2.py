# -*- coding: utf-8 -*-
import sys
sys.path.insert(0, "/export/home/ab248677/CALCULS_IJK/RUN_GUILLAUME")

from Tools import *
from user import *
import readline
import unittest
import os
import pickle
from matplotlib.pyplot import clf
readline.parse_and_bind('tab:complete')
clf
global compteur
compteur=0
os.system('rm -f *.png')
print ''


##################################################################################################################################
########################################################## CALCULS ############################################################
##################################################################################################################################
if 0:
    fstat = "Stats.dt_ev_mono_nursaf"   
    rho=594.38
    mu=6.8327e-5
    Source=2.0358956532803925
     
    Mono=RunMonophasicCase(fstat, rho, mu, Source)
     
    fstat1 = "Stats_vreman1.dt_ev"
    fstat2 = "Stats_vreman2.dt_ev"
    rho=1.0
    mu=0.001694915
    Source=1
    Vreman=RunVremanCase(fstat1, fstat2, rho, mu, Source)
    
    gravity=-0.1
    rho_l=1.0
    rho_v=0.1
    alpha_v=0.0304
    mu=0.00016666
    Source=0.0017    
    sigma=0.002  #sigma=0.00512    
    fstat = "Stats_deform.dt_ev"   #"Stats_deform/sphere.dt_ev/Stats.dt_ev_gui"
    DeformableRun=RunDiphasicCase(gravity, rho_l, rho_v, alpha_v, mu, Source, sigma, fstat)
    save(DeformableRun, 'data_DeformableRun')
     
    sigma=0.02
    fstat = "Stats_sphere.dt_ev"   #"Stats_deform/sphere.dt_ev/Stats.dt_ev_gui"
    SphericalRun=RunDiphasicCase(gravity, rho_l, rho_v, alpha_v, mu, Source, sigma, fstat)
      
    fstat = "Stats.dt_ev_gui"
    gravity=-9.81
    rho_l=594.38
    rho_v=101.93
    alpha_v=0.1
    mu=6.8327e-5
    Source=2.03589565328039246
    sigma=0.0046695
    GuillaumeRun=RunDiphasicCase(gravity, rho_l, rho_v, alpha_v, mu, Source, sigma, fstat)
     

    save(GuillaumeRun, 'data_GuillaumeRun')
    save(SphericalRun, 'data_SphericalRun')
    save(Vreman, 'data_Vreman')
    save(Mono, 'data_Mono')


if 0:

    gravity=-0.1
    rho_l=1.0
    rho_v=0.1
    alpha_v=0.06
    mu=0.00016666
    Source=0.098128365  
    sigma=0.002320  #sigma=0.00512    
    fstat0 = "STAT_GROUP0.dt_ev" 
    fstat1 = "STAT_GROUP1.dt_ev" 
    fstat2 = "STAT_GROUP2.dt_ev" 
    fstat012 = "STAT_GROUP012.dt_ev" 
    fstatref = "VALID_BIDI_REF.dt_ev" 
    DeformableRun=RunDiphasicCase(gravity, rho_l, rho_v, alpha_v, mu, Source, sigma, fstat0)
    save(DeformableRun, 'STAT_GROUP0')
    DeformableRun=RunDiphasicCase(gravity, rho_l, rho_v, alpha_v, mu, Source, sigma, fstat1)
    save(DeformableRun, 'STAT_GROUP1')
    DeformableRun=RunDiphasicCase(gravity, rho_l, rho_v, alpha_v, mu, Source, sigma, fstat2)
    save(DeformableRun, 'STAT_GROUP2')
    DeformableRun=RunDiphasicCase(gravity, rho_l, rho_v, alpha_v, mu, Source, sigma, fstat012)
    save(DeformableRun, 'STAT_GROUP012')
    DeformableRun=RunDiphasicCase(gravity, rho_l, rho_v, alpha_v, mu, Source, sigma, fstatref)
    save(DeformableRun, 'STAT_REF_BIDI')


##################################################################################################################################
########################################################## TRAITEMENT ############################################################
##################################################################################################################################

#Vreman=load('data_Vreman')
#Mono=load('data_Mono')
GuillaumeRun=load('data_GuillaumeRun')
SphericalRun=load('data_SphericalRun')
DeformableRun=load('data_DeformableRun')


#Mk(GuillaumeRun)

tracerQdmDiff(SphericalRun)


if 0:
    ref=load('STAT_REF_BIDI')
    gr0=load('STAT_GROUP0')
    gr1=load('STAT_GROUP1')
    gr2=load('STAT_GROUP2')
    gr012=load('STAT_GROUP012')
    tracer([ref['variable']['I'],gr0['variable']['I'], gr1['variable']['I'], gr2['variable']['I'],  gr012['variable']['I']], 'Igroup')
    tracer([ref['variable']['aiN'],gr0['variable']['aiN'], gr1['variable']['aiN'], gr2['variable']['aiN'],  gr012['variable']['aiN']], 'aiNxgroup', [0])
    tracer([ref['variable']['aiN'],gr0['variable']['aiN'], gr1['variable']['aiN'], gr2['variable']['aiN'],  gr012['variable']['aiN']], 'aiNygroup', [1])
    tracer([ref['variable']['aiN'],gr0['variable']['aiN'], gr1['variable']['aiN'], gr2['variable']['aiN'],  gr012['variable']['aiN']], 'aiNzgroup', [2])
    tracer([ref['variable']['aii'],gr0['variable']['aii'], gr1['variable']['aii'], gr2['variable']['aii'],  gr012['variable']['aii']], 'aiixgroup', [0])
    tracer([ref['variable']['aii'],gr0['variable']['aii'], gr1['variable']['aii'], gr2['variable']['aii'],  gr012['variable']['aii']], 'aiiygroup', [1])
    tracer([ref['variable']['aii'],gr0['variable']['aii'], gr1['variable']['aii'], gr2['variable']['aii'],  gr012['variable']['aii']], 'aiizgroup', [2])
    tracer([ref['qdm_diph']['source'], gr0['qdm_diph']['source'], gr1['qdm_diph']['source'], gr2['qdm_diph']['source'], gr012['qdm_diph']['source']], 'sourceqdmgroup', [0])
    tracer([ref['qdm_diph']['pression'], gr0['qdm_diph']['pression'], gr1['qdm_diph']['pression'], gr2['qdm_diph']['pression'], gr012['qdm_diph']['pression']], 'pressionxqdmgroup', [0])
    tracer([ref['qdm_diph']['pression'], gr0['qdm_diph']['pression'], gr1['qdm_diph']['pression'], gr2['qdm_diph']['pression'], gr012['qdm_diph']['pression']], 'pressionyqdmgroup', [1])
    tracer([ref['qdm_diph']['pression'], gr0['qdm_diph']['pression'], gr1['qdm_diph']['pression'], gr2['qdm_diph']['pression'], gr012['qdm_diph']['pression']], 'pressionzqdmgroup', [2])
    tracer([ref['qdm_diph']['inertie'], gr0['qdm_diph']['inertie'], gr1['qdm_diph']['inertie'], gr2['qdm_diph']['inertie'], gr012['qdm_diph']['inertie']], 'inertiexqdmgroup', [0])
    tracer([ref['qdm_diph']['inertie'], gr0['qdm_diph']['inertie'], gr1['qdm_diph']['inertie'], gr2['qdm_diph']['inertie'], gr012['qdm_diph']['inertie']], 'inertieyqdmgroup', [1])
    tracer([ref['qdm_diph']['inertie'], gr0['qdm_diph']['inertie'], gr1['qdm_diph']['inertie'], gr2['qdm_diph']['inertie'], gr012['qdm_diph']['inertie']], 'inertiezqdmgroup', [2])
    tracer([ref['qdm_diph']['interface'], gr0['qdm_diph']['interface'], gr1['qdm_diph']['interface'], gr2['qdm_diph']['interface'], gr012['qdm_diph']['interface']], 'interfacexqdmgroup', [0])
    tracer([ref['qdm_diph']['interface'], gr0['qdm_diph']['interface'], gr1['qdm_diph']['interface'], gr2['qdm_diph']['interface'], gr012['qdm_diph']['interface']], 'interfaceyqdmgroup', [1])
    tracer([ref['qdm_diph']['interface'], gr0['qdm_diph']['interface'], gr1['qdm_diph']['interface'], gr2['qdm_diph']['interface'], gr012['qdm_diph']['interface']], 'interfacezqdmgroup', [2])
    tracer([ref['qdm_diph']['turbulence'], gr0['qdm_diph']['turbulence'], gr1['qdm_diph']['turbulence'], gr2['qdm_diph']['turbulence'], gr012['qdm_diph']['turbulence']], 'turbulencexqdmgroup', [0])
    tracer([ref['qdm_diph']['turbulence'], gr0['qdm_diph']['turbulence'], gr1['qdm_diph']['turbulence'], gr2['qdm_diph']['turbulence'], gr012['qdm_diph']['turbulence']], 'turbulenceyqdmgroup', [1])
    tracer([ref['qdm_diph']['turbulence'], gr0['qdm_diph']['turbulence'], gr1['qdm_diph']['turbulence'], gr2['qdm_diph']['turbulence'], gr012['qdm_diph']['turbulence']], 'turbulencezqdmgroup', [2])
    tracer([ref['qdm_diph']['viscous'], gr0['qdm_diph']['viscous'], gr1['qdm_diph']['viscous'], gr2['qdm_diph']['viscous'], gr012['qdm_diph']['viscous']], 'viscousxqdmgroup', [0])
    tracer([ref['qdm_diph']['viscous'], gr0['qdm_diph']['viscous'], gr1['qdm_diph']['viscous'], gr2['qdm_diph']['viscous'], gr012['qdm_diph']['viscous']], 'viscousyqdmgroup', [1])
    tracer([ref['qdm_diph']['viscous'], gr0['qdm_diph']['viscous'], gr1['qdm_diph']['viscous'], gr2['qdm_diph']['viscous'], gr012['qdm_diph']['viscous']], 'viscouszqdmgroup', [2])
    tracer([ref['qdm_diph']['rhog'], gr0['qdm_diph']['rhog'], gr1['qdm_diph']['rhog'], gr2['qdm_diph']['rhog'], gr012['qdm_diph']['rhog']], 'rhogxqdmgroup', [0])
    tracer([ref['qdm_diph']['rhog'], gr0['qdm_diph']['rhog'], gr1['qdm_diph']['rhog'], gr2['qdm_diph']['rhog'], gr012['qdm_diph']['rhog']], 'rhogyqdmgroup', [1])
    tracer([ref['qdm_diph']['rhog'], gr0['qdm_diph']['rhog'], gr1['qdm_diph']['rhog'], gr2['qdm_diph']['rhog'], gr012['qdm_diph']['rhog']], 'rhogzqdmgroup', [2])











