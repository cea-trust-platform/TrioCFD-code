# -*- coding: utf-8 -*-
import rlcompleter, readline
from matplotlib.pyplot import clf
from cmath import sin
from scipy import misc
from matplotlib.axes._subplots import Subplot


readline.parse_and_bind('tab:complete')
# import Tools
from Tools import Field
import numpy as np
import matplotlib.pyplot as plt
#############################""


####### Produit simple scalaire/scalaire scalaire/vecteur ###########################
def product(aField, dField, sym='*'):
    
    product=Field()
    product.__init__(aField.name + sym + dField.name ,aField.label + sym + dField.label , aField._ax + sym + dField._ax)
    aa = len(aField._npa[0,:])    
    ab = len(aField._npa[:,0])
    da = len(dField._npa[0,:])    
    db = len(dField._npa[:,0]) 
  
    product._npa=np.zeros((max(ab,db),max(aa,da)))
               
    if ab==1 and aa==1:
        for i in range(0, len(dField._npa[:,0]) ):
            for j in range(0, len(dField._npa[0,:]) ):
                if sym=='*': 
                    product._npa[i,j]=aField._npa[0,0]*dField._npa[i,j]
                if sym=='/': 
                    product._npa[i,j]=dField._npa[i,j]/aField._npa[0,0]                   
    
    elif db%ab==0 and aa==da:
        dim=db/ab
        for i in range(0, len(dField._npa[:,0]) ):
            a=i/dim
            for j in range(0, len(dField._npa[0,:]) ):
                if sym=='*':
                    product._npa[i,j]=aField._npa[a,j]*dField._npa[i,j]
                if sym=='/':  
                    product._npa[i,j]=dField._npa[i,j]/aField._npa[a,j]
    else :
        print 'Ce produit est impossible a effectuer'
                
    return product



################## tracer des courbes resultats ########################
def tracer(afield, bfield, title, subplot=0, ylim=None, xlim=None, grid=True, label=True):
    global compteur
    
    if subplot==0:
        compteur=compteur+1
        plt.figure(compteur)
    i=0
    while i<len(bfield):
        
        j=0
        while j<len(bfield[i]._npa):
            if label:
                plt.plot(afield._npa[0,:],bfield[i]._npa[j,:], label=bfield[i].label+str(j))
                j+=1
            else:
                plt.plot(afield._npa[0,:],bfield[i]._npa[j,:])
                j+=1
        i+=1        
        
    plt.legend(loc = 0)
    if grid:
        plt.grid('on')
    if ylim != None:
        plt.ylim(ylim)
    if xlim != None:
        plt.xlim(xlim)     
    #plt.title(title)   
    plt.xlabel(afield._ax)
    plt.ylabel(bfield[0]._ax)
    sauv=title+'.png'
    if subplot==0:
        plt.savefig(sauv)
    return

#################### classement des derivees ##################

def derivees(aField, dField, n='1', ordre='1', bc_type=0):
    
    grad=Field()
    grad.__init__('derivee de '+aField.name , 'grad_'+aField.label, 'D'+aField._ax+'_D'+dField._ax)
    
    a = len(aField._npa[0,:])
    b = len(aField._npa[:,0])
    grad._npa=np.zeros((b,a))
        
    i=0
    while i<=len(aField._npa[:,0])-1:
        if n=='1':
            if ordre=='1':
                grad._npa[i,:]= cell_to_cell_gradient_scalar(aField._npa[i,:], dField._npa[0,:], bc_type)
                i+=1   
    return grad   

def cell_to_cell_gradient_scalar(Field, dField, bc_type):
    
    gradi = np.zeros(len(Field))
  
    gradi[1:-1] = (Field[2:]-Field[:-2])/(dField[2:]-dField[:-2])
    
    if bc_type==0:
        gradi[0]=(-0.166666666*Field[6]+1.2*Field[5]-3.75*Field[4]+6.66666666*Field[3]-7.5*Field[2]+6*Field[1]-2.45*Field[0])/(dField[1]-dField[0])
        gradi[-1]=(0.166666666*Field[-7]-1.2*Field[-6]+3.75*Field[-5]-6.66666666*Field[-4]+7.5*Field[-3]-6*Field[-2]+2.45*Field[-1])/(dField[-1]-dField[-2])
       
    return gradi   

def reshapeScalarField(field):
    a = len(field._npa[:])
   # newField=Field()
    field._npa=field._npa.reshape(1,a)
    return field


def somme(aField, dField, sym='+'):
    som=Field()
    som.__init__(aField.name + sym + dField.name ,aField.label + sym + dField.label , aField._ax + sym + dField._ax)
    aa = len(aField._npa[0,:])
    ad = len(dField._npa[0,:])
    ba = len(aField._npa[:,0])
    bd = len(dField._npa[:,0])
    som._npa=np.zeros((ad,aa))
    if aa!=ad or ba!=bd:
        print 'ces champs ne peuvent pas être additionés'
        return
    for i in range(0, len(aField._npa[:,0]) ):
        for j in range(0, len(aField._npa[0,:]) ):
            if sym=='+':
                som._npa[i,j]=aField._npa[i,j]+dField._npa[i,j]
            if sym=='-':
                som._npa[i,j]=aField._npa[i,j]-dField._npa[i,j]    
    return som

def integral(aField, dField, j_min=0, j_max=0):
    
    j_min=0
    j_max=len(aField._npa[0,:])
    integral=Field()
    dim = len(aField._npa[:,0])
    integral._npa=np.zeros((dim,1))
       
    for i in range(0, len(aField._npa[:,0]) ):
            somme=0
            for j in range(0, len(aField._npa[i,:])-1 ):  
                if  j_min < j < j_max :
                    somme=somme+aField._npa[i,j]*(dField._npa[0,j+1]-dField._npa[0,j])
            integral._npa[i,0]=somme
    return integral
#######################################################################################
###################################### DEBUT DU CODE #####################################"
#########################################################################################


############# Initialisation des objets lu, lecture et reshape des scalaires ###############
global compteur

compteur=0

fstat = "Stats.dt_ev"
dvar = Field.getEntries(fstat)

k = Field.LoadFromFile(fstat, ["coordonnee_K"])
I = Field.LoadFromFile(fstat, ["I"])
U = Field.LoadFromFile(fstat, ["UI", "VI", "WI"]) 
Mat = Field.LoadFromFile(fstat, ["UUI", "UVI", "UWI", "UVI", "VVI", "VWI", "UWI", "VWI", "WWI"]) 
k=reshapeScalarField(k)
I=reshapeScalarField(I)
############################# Calcul des grandeurs #############################
# 
# 
# gradI = derivees(I,k)
# grad2I= derivees(gradI,k)
# grad = derivees(U,k)
# productU = product(I,U)
# product=product(I,I)
# 
# ############################### OUTPUT #####################################
# 
# tracer(k,[grad], 'Gradient_U')
# tracer(k,[gradI], 'Gradient_I')
# tracer(k,[grad2I], 'Gradient2_I')
# tracer(k,[product], 'Product_II')
# tracer(k,[productU], 'Product_IU')

#############################################################################
################################ VALIDATION #################################
#############################################################################


##################initialisation des champs scalaire/vecteur/tenseur ##########
x=Field()
y=Field()
V=Field()
M=Field()
T=Field()

x.__init__('x', 'x', 'x')
y.__init__('champ scalaire', 'y', 'y')
V.__init__('champ vectoriel', 'V', 'V')
M.__init__('champ vectoriel', 'M', 'M')
T.__init__('champ tensoriel', 'T', 'T')


rho=Field()
rho._npa=np.array([4], float)
x._npa=np.linspace(-5,5,100)
y._npa=np.sin(x._npa)
x=reshapeScalarField(x)
y=reshapeScalarField(y)
rho=reshapeScalarField(rho)
a = len(x._npa[0,:])
M._npa=np.zeros((3,a))
M._npa[0,:]=1
M._npa[1,:]=np.cos(x._npa)
M._npa[2,:]=np.sqrt(np.abs(x._npa))
a = len(x._npa[0,:])
V._npa=np.zeros((3,a))
V._npa[0,:]=x._npa[:]
V._npa[1,:]=0.5*x._npa*x._npa
V._npa[2,:]=0.166666*x._npa*x._npa*x._npa

T._npa=np.zeros((9,a))
uu=np.zeros((1,a))
vv=np.zeros((1,a))
ww=np.zeros((1,a))


uu=np.cos(2*x._npa)
vv=1.3
ww=np.sqrt(np.sqrt(np.abs(x._npa)))
# T._npa[0,:]=uu*uu
# T._npa[1,:]=uu*vv
# T._npa[2,:]=uu*ww
# T._npa[3,:]=vv*uu
# T._npa[4,:]=vv*vv
# T._npa[5,:]=vv*ww
# T._npa[6,:]=ww*uu
# T._npa[7,:]=ww*vv
# T._npa[8,:]=ww*ww


T._npa[0,:]=1
T._npa[1,:]=2
T._npa[2,:]=3
T._npa[3,:]=4
T._npa[4,:]=5
T._npa[5,:]=6
T._npa[6,:]=7
T._npa[7,:]=8
T._npa[8,:]=9
###### validation sclaire ###########
 
# grady = derivees(y,x)
# grad2y = derivees(grady,x)
# grad3y = derivees(grad2y,x)
# grad4y = derivees(grad3y,x)
# product1= product(y,y)
# product2= product(y,y,'/')
# somme1= somme(y,y)
# somme2=somme(y,y,'-')
#  
# plt.figure(111)
# Nx=2
# Ny=3
# plt.subplot(Nx,Ny,1)
# tracer(x, [y],'Field', 1 )
# plt.subplot(Nx,Ny,2)
# tracer(x, [product1],'y*y', 1 )
# plt.subplot(Nx,Ny,3)
# tracer(x, [product2],'y/y', 1 )
# plt.subplot(Nx,Ny,4)
# tracer(x, [somme1],'y+y', 1 )
# plt.subplot(Nx,Ny,5)
# tracer(x, [somme2],'y-y', 1 )
# plt.subplot(Nx,Ny,6)
# tracer(x, [y, grady, grad2y, grad3y, grad4y], 'derivee_sca',1)
# plt.savefig('scalaire.png')



########### validation vecteur ##########################

# gradV = derivees(V,x)
# grad2V = derivees(gradV,x)
# grad3V = derivees(grad2V,x)
# som1=somme(V,V)
# diff=somme(V,V,'-')
# product2 = product(V,V)
# product3 = product(y,V)
# product4 = product(V,V,'/')
# product5 = product(y,V, '/')
# # integral=integral(V, x)
# # print integral._npa
# fonction1 = somme(V,product(y,derivees(V,x)))
#  
# plt.figure(112)
# Nx=4
# Ny=3
# plt.subplot(Nx,Ny,1)
# tracer(x, [V],'Field', 1 , [-5,5])
# plt.subplot(Nx,Ny,2)
# tracer(x, [gradV],'dv_dx', 1 , [-5,5])
# plt.subplot(Nx,Ny,3)
# tracer(x, [grad2V],'d2v_dx2', 1 , [-5,5])
# plt.subplot(Nx,Ny,4)
# tracer(x, [grad3V],'d3v_dx3', 1 , [-5,5])
# plt.subplot(Nx,Ny,5)
# tracer(x, [som1],'vec+vec', 1 , [-5,5])
# plt.subplot(Nx,Ny,6)
# tracer(x, [diff],'vec-vec', 1 , [-5,5])
# plt.subplot(Nx,Ny,7)
# tracer(x, [product2],'v*v', 1, [0,5] )
# plt.subplot(Nx,Ny,8)
# tracer(x, [product3],'y*V', 1 , [-5,5])
# plt.subplot(Nx,Ny,9)
# tracer(x, [product4],'V/V', 1 , [-5,5])
# plt.subplot(Nx,Ny,10)
# tracer(x, [product5],'V/y', 1 , [-20,20])
# plt.subplot(Nx,Ny,11)
# tracer(x, [y],'y', 1)
# plt.subplot(Nx,Ny,12)
# tracer(x, [fonction1],'y',[-5,5])
# plt.savefig('vecteur.png')


###################### validation tenseur ####################################


gradT = derivees(T,x)


grad2T = derivees(gradT,x)
grad3T = derivees(grad2T,x)
sommeT=somme(T,T)
diffT=somme(T,T,'-')
productT = product(T,T)
productT3 = product(y,T)
productT6 = product(M,T)
productT4 = product(T,T,'/')
productT5 = product(y,T, '/')
productT7 = product(M,T, '/')
integral=integral(T, x)
print integral._npa


tracer(x, [T],'T')
tracer(x, [gradT],'dT_dx')
tracer(x, [sommeT],'T_p_T')
tracer(x, [diffT],'T_m_T')
tracer(x, [productT],'TT')
tracer(x, [productT3],'ST')
tracer(x, [productT4],'T_d_T')
tracer(x, [productT5],'T_d_S', [-20,20])
tracer(x, [y],'S')
tracer(x, [M],'V')
tracer(x, [productT6],'VT')
tracer(x, [productT7],'T_d_V', [-20,20])
















#UUl = Field.LoadFromFile(fstat, ["UI", "VI", "WI","UI", "VI", "WI","UI", "VI", "WI"]) # --> TensorField (2nd order, general, does not assume symetry...)

#IUl = I*Ul

#rho = ScalarField()
#rho.label = 
#u = VectorField()

#new_term = rho*u

#derivate(Rij) + div() 

