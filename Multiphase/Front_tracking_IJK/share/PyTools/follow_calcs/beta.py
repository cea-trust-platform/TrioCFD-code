#!usr/bin/python
import os, sys, readline, rlcompleter
readline.parse_and_bind('tab:complete')
import DNSTools as dtool
from math import pi

jdd='DNS.data'
try:
   repr_file=dtool.getParam(jdd, 'nom_reprise', string=True)
except:
   repr_file="diph.sauv"
   pass

print "repr_file= ", repr_file


vb=dtool.getParam(jdd, 'vol_bulle_monodisperse')
nb=dtool.getParam(repr_file, 'bubble_groups')
sigma=dtool.getParam(jdd, 'sigma')


if (len(sys.argv) !=2): 
   raise Exception('the script %s expects exactly one parameter giving the value of Re_tau.'%sys.argv[0])

Ret=float(sys.argv[1])
rhol, rhov, mul, muv, alv, beta, sigma, _ = dtool.get_prop(jdd='DNS.data', repr_file=repr_file, Ret=Ret)
simu = dtool.Simu(jdd='DNS.data', repr_file=repr_file, Ret=Ret)
