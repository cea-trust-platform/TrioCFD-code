from sys import *
import os

"""
Ecrit dans un fichier .curve tout la donnee du champ demande par l'utilisateur, a l'instant demande par l'utilisateur.
"""

##########################################################################################
# mode d'emploi : visit -s vs_extract_all_volume.py
##########################################################################################
import sys

pwd = os.getcwd()

OpenDatabase("./cube.lata")

AddPlot("Spreadsheet", "INDICATRICE_ELEM_DOM", 1, 1)
DrawPlots()
ExportDBAtts = ExportDBAttributes()
ExportDBAtts.allTimes = 1
ExportDBAtts.dirname = "./"
ExportDBAtts.filename = "INDICATRICE_ELEM_DOM"
ExportDBAtts.timeStateFormat = "_%04d"
ExportDBAtts.db_type = "Curve2D"
ExportDBAtts.db_type_fullname = "Curve2D_1.0"
ExportDBAtts.variables = ("INDICATRICE_ELEM_DOM")
ExportDBAtts.writeUsingGroups = 0
ExportDBAtts.groupSize = 48
ExportDBAtts.opts.types = (5)
ExportDBAtts.opts.help = ""
ExportDatabase(ExportDBAtts)

exit()
