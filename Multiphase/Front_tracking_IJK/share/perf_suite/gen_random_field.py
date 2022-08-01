from math import *
from random import *

lx=0.093764*.03777953*8
ly=0.029846
lz=0.027*.03777953*16
s=""
for ni in range(1,4):
    for nj in range(1,4):
        for nk in range(1,2):
            if (len(s)>0):
                s=s+"+"
                pass
            a=random()*2.*pi
            b=random()*10.0/ni
            s=s+("sin(x*%.2f+%.2f)*%.2f" % (ni*ni*2*pi/lx,a,b))
            a=random()*2.*pi
            b=random()*10.0/nj
            s=s+("*sin(y*%.2f+%.2f)*%.2f" % (nj*nj*2*pi/ly,a,b))
            a=random()*2.*pi
            b=random()*10.0/nk
            s=s+("*sin(z*%.2f+%.2f)*%.2f" % (nk*nk*2*pi/lz,a,b))
            pass
        pass
    pass
print "z=0"
print "splot", s
print "pause -1"
