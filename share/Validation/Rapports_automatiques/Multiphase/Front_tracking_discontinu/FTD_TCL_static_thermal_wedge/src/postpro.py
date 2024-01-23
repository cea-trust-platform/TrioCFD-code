#!visit -nowin -cli -s 
OpenDatabase("./lata/post.lata", 0)
DefineScalarExpression("dTdni", "recenter(pos_cmfe(<[0]id:TEMPERATURE_GRAD_THERMIQUE_ELEM_dom>, INTERFACES, -1e8), \"nodal\")")
AddPlot("Pseudocolor", "dTdni", 1, 1)
DrawPlots()
#DefineScalarExpression("xi", "coord(INTERFACES)[0]")
#DefineScalarExpression("yi", "coord(INTERFACES)[1]")
nb = TimeSliderGetNStates()
SetTimeSliderState(nb-1)
ExportDBAtts = ExportDBAttributes()
ExportDBAtts.allTimes = 0
ExportDBAtts.db_type = "XYZ"
ExportDBAtts.db_type_fullname = "XYZ_1.0"
ExportDBAtts.filename = "dTdni"
ExportDBAtts.dirname = "."
ExportDBAtts.variables = ("dTdni")
ExportDBAtts.writeUsingGroups = 0
ExportDBAtts.groupSize = 48
ExportDBAtts.opts.types = ()
ExportDBAtts.opts.help = ""
ExportDatabase(ExportDBAtts)

import math
import numpy as np
import sys

if (len(sys.argv) != 2):
   raise Exception("One argument is compulsory to set sm")

# print(f"sm set to {sys.argv[1]}")
sm = float(sys.argv[1])
print(sm*1.e6)

x, y, z, dTdn = np.loadtxt("dTdni.xyz", skiprows=2, usecols=(1,2,3,4)).T
def readFileInfo():
    data = {}
    with open("info.txt") as f:
        for line in f.readlines():
            key, value = line.rstrip("\n").split("=")
            data[key] = float(value)
    return data

data = readFileInfo()
# Curvilinear absissa : 
s = np.sqrt(x*x+y*y+z*z)
s -= s[0]

lda = data["lda"]
dT = data["DT"]
theta_app = data["theta"]*math.pi/180
qi = lda*dTdn
# Analytical solution from Vadim:
# qi_ana = lda*dT/(s[1:]*theta_app)
Ri = 6.38209e-08

qi_ana = []
for element in s[1:]:
    if element <= sm:
        qi_ana.append(lda * dT / (element * theta_app))
    else:
        qi_ana.append(lda * dT / (element * theta_app + Ri * lda))

# If s[1:] is a NumPy array or a list, you can convert qi_ana back to the same type
qi_ana = np.array(qi_ana)  # For NumPy array

np.savetxt("s_qi_qiana.txt", np.c_[s[1:],qi[1:],qi_ana], header='s qi qi_ana')

exit()
