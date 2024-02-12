#!/usr/bin/python3

# ------------------------------------------------------------------------------
# Boite a outils pour les post-traitement du IJK...
# ------------------------------------------------------------------------------

import numpy as np
import math
import glob
import os
import re
import json


def pwd():
    return os.getcwd()


def getJddName():
    l = glob.glob("*data")
    if len(l) > 1:
        l = [s for s in l if "repr" not in s]
    jdd = l[0][:-5]
    return jdd


def getTempsIntegrationFile(fic):
    f = open(fic, "r")
    lines = f.readlines()
    tintegration = lines[0].split()[2]
    f.close()
    return float(tintegration)


def getTempsIntegrationList(fics):
    lt = []
    for fic in fics:
        lt.append(getTempsIntegrationFile(fic))
    return np.array(lt)


def buildDicoColonnes(fic):
    d = {}
    f = open(fic, "r")
    lines = f.readlines()
    for val, key in [
        (int(s.split()[2]) - 1, s.split()[4]) for s in lines if "# colonne" in s
    ]:
        d[key] = val
    f.close()
    return d


def getStatsInstant():
    return getStats("ins")


def getStatsMoy():
    return getStats()


# caract = moy or ins : instantant ou moyen.
# Retourne une matrice avec les stats :
# resu[i,j,k] : i --> temps
#                    j --> Position en z
#                    k --> var.
#
# ltimes : liste des temps.
# lz      : liste des coords en z.
# lvar    : liste des variables stockees.
# val[it, iz, ivar] : La matrice 3D contenant tout ca.
#
# tintegration : le temps d'integration des resulats si caract == moy.
def getStats(caract="moy"):
    if caract == "moy":
        fic_stats = glob.glob("statistiques_*.txt")
        if fic_stats == []:
            fic_stats = glob.glob("*phasique_statistiques_*.txt")
        fic_stats.sort(
            key=lambda item: float(item.strip("monodiphasique_statistiques_.txt"))
        )
        ltimes = np.array(
            [float(f.strip("monodiphasique_statistiques_.txt")) for f in fic_stats]
        )
    else:
        fic_stats = glob.glob("moyenne_spatiale_*.txt")
        if fic_stats == []:
            fic_stats = glob.glob("*phasique_moyenne_spatiale_*.txt")
        fic_stats.sort(
            key=lambda item: float(item.strip("monodiphasique_moyenne_spatiale_.txt"))
        )
        ltimes = np.array(
            [float(f.strip("monodiphasique_moyenne_spatiale_.txt")) for f in fic_stats]
        )
    if fic_stats == []:
        raise Exception("La liste fic_stats est restee vide... ")
    lvar = buildDicoColonnes(fic_stats[0])
    nt = len(ltimes)
    if caract == "moy":
        tintegration = getTempsIntegrationList(fic_stats)
    else:
        tintegration = [0]
    for i, fic in enumerate(fic_stats):
        mat = np.loadtxt(fic)
        t = ltimes[i]
        lz = mat[:, 0]
        if i == 0:
            # Initialise la taille du resu :
            nz, nvar = np.shape(mat)
            resu = np.zeros((nt, nz, nvar))
        resu[i, :, :] = mat[:, :]
    return ltimes, lz, lvar, resu, tintegration


def getValue(key, fic, pre="^", default=None):
    f = open(fic, "r")
    lines = f.readlines()
    f.close()
    nb = len(lines)
    rc = re.compile(
        pre + "\s*" + key + "\s*(?P<value>[\-]?[\d]*.?[\d]*[eE]?[+\-]?[\d]*)"
    )
    for i, st in enumerate(lines):
        m = rc.match(st)
        if m:
            val = float(m.group("value"))
            return val
            break
    if default != None:
        return default
    print("On n'a rien trouve dans ", fic)
    raise Exception("Etonnant, non?")


# En z+ et En z- :
def evaluateRetau(U, z, jdd=None):
    if jdd == None:
        jdd = getJddName() + ".data"
    rhol = getValue("rho_liquide", jdd)
    mul = getValue("mu_liquide", jdd)
    Lz = getValue("uniform_domain_size_k", jdd)
    h = Lz / 2.0
    tauwp = mul * U[:, 0] / z[0]
    tauwm = mul * U[:, -1] / (Lz - z[-1])
    #
    tauw = (tauwp + tauwm) / 2.0
    utau = np.array([np.sqrt(abs(x) / rhol) for x in tauw])
    Retau = rhol * utau * h / mul
    return Retau, tauwp, tauwm

