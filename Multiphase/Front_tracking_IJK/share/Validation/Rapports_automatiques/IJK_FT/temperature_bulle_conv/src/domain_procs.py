# -*- coding: utf-8 -*-

#import sympy as sym
import numpy as np
import matplotlib.pyplot as plt
import sympy as sym
from sympy import solveset, Interval, FiniteSet
import sys

#%%
a=float(sys.argv[1])/2 # Radius
mu_liq=float(sys.argv[2]) # Liquid viscosity
mu_vap=float(sys.argv[3]) # Vapour viscosity
vel_percent=float(sys.argv[4]) # Theoretical percentage of Terminal Velocity 
nb_elem_per_radius=float(sys.argv[5])

tolerance=1-(vel_percent/100.0)

yz_width=5*a
yz_origin=yz_width/2

yz_elements=round(yz_width/a*nb_elem_per_radius)
if yz_elements%2 == 1:
    yz_elements+=1
yz_procs=round(yz_elements*(1/32+1/64)/2)
if yz_procs%2 == 1:
    yz_procs+=1
#%%
mu_tilde=mu_vap/mu_liq

print('Percentage of Hadamard terminal velocity: ' + str(vel_percent) +'\%')
print('Viscosity ratio: '+ str(mu_tilde))

#%%
z=sym.symbols('z', real=True)
rel_velocity=sym.Abs(-a**3/(2*z**3)*mu_tilde/(1+mu_tilde) + a/(2*z)*(2+3*mu_tilde)/(1+mu_tilde))
channel_height=float(solveset(rel_velocity-tolerance,z,Interval(a,1000*a)).args[0])

#%%
print('Channel length to reach ' + str(vel_percent) + ' percent of terminal velocity: ' + str(channel_height/a) + 'a')

#%%
channel_origin=channel_height/2
channel_elements=round(channel_height/a*nb_elem_per_radius)
if channel_elements%2 == 1:
    channel_elements+=1
channel_procs=round(channel_elements*(1/32+1/64)/2)
if channel_procs%2 == 1:
    channel_procs+=1

channel_elements=channel_elements+(channel_procs-channel_elements%channel_procs)
yz_elements=yz_elements+(yz_procs-yz_elements%yz_procs)

n_elem_over_n_procs=(yz_elements**2*channel_elements)/(yz_procs**2*channel_procs)
print('Check elem/processors ratio: ' + str(n_elem_over_n_procs))
print('Procs for the channel height: ' + str(channel_procs) + ', for the transverse direction: ' + str(yz_procs))

#%%
parameters='( "' + str(channel_height) + '" "' + str(channel_origin) + '" "' + str(channel_elements) + '" "' + str(yz_width) + '" "' + str(yz_origin) + '" "' + str(yz_elements) + '" "' + str(channel_procs) + '" "' + str(yz_procs) + '" "' + str(n_elem_over_n_procs) + '" )'
sys.exit(parameters)

