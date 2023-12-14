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

