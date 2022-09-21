# -*- coding: utf-8 -*-
#import rlcompleter
from Tools import *
import readline
import unittest
import pickle
import os


from matplotlib.pyplot import clf


readline.parse_and_bind('tab:complete')


clf
global compteur
compteur=0
os.system('rm -f *.png')




def Init_case() : 
    """ 
    Créer un dictionnaire de sauvegarde
    """        

    res = {}
    return res

def Add_section(dic, name) :

    """ 
    Créee une section dans un dictionnaire
    @param un dictionnaire
    @param un nom
    """  

    section = {}
    dic[name] = section
    return 
################################################################################
########### PLOT DES PROFIL DE VITESSE ET INDICATRICE ##########################
################################################################################
def UsualVariableuluv() :

    """ 
    Comparaison de différents cas (Ss127, D127, SP127, B127) codés en dur. Ne pas appeler cette fonction sur de nouveaux cas, elles ne fonctionnera pas.
    """ 
    ################################################################################
    ########### PLOT DES PROFIL DE VITESSE ET INDICATRICE ##########################
    ################################################################################


    ################## vitesse plot ###############
    fstat = "Stats_deform.dt_ev"
    II=Field.LoadFromFile(fstat,["coordonnee_K"] , ["I"],'I', r'$z$', r'I')
    Iv=(II-1.0)*(-1.0)
    Iinv=II.inv(1e-8)
    Ivinv=Iv.inv(1e-8)
    uD_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UI", "VI", "WI"],r'$U_l$', r'$z$', r'$u_l$')*Iinv
    uD_v=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UIv", "VIv", "WIv"],'U_v', r'$z$', r'$u_v$')*Ivinv
    Vr_D=uD_v-uD_l

    fstat = "Stats_sphere.dt_ev"
    II=Field.LoadFromFile(fstat,["coordonnee_K"] , ["I"],'I', r'$z$', r'I')
    Iv=(II-1.0)*(-1.0)
    Iinv=II.inv(1e-8)
    Ivinv=Iv.inv(1e-8)
    uSb_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UI", "VI", "WI"],r'$U_l$', r'$z$', r'$u_l$')*Iinv
    uSb_v=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UIv", "VIv", "WIv"],'U_v', r'$z$', r'$u_v$')*Ivinv
    Vr_Sb=uSb_v-uSb_l

    fstat = "Stats_bidi127.dt_ev"
    II=Field.LoadFromFile(fstat,["coordonnee_K"] , ["I"],'I', r'$z$', r'I')
    Iv=(II-1.0)*(-1.0)
    Iinv=II.inv(1e-8)
    Ivinv=Iv.inv(1e-8)
    uB_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UI", "VI", "WI"],r'$U_l$', r'$z$', r'$u_l$')*Iinv
    uB_v=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UIv", "VIv", "WIv"],'U_v', r'$z$', r'$u_v$')*Ivinv
    Vr_B=uB_v-uB_l

    fstat = "Stats_petitesphere.dt_ev"
    II=Field.LoadFromFile(fstat,["coordonnee_K"] , ["I"],'I', r'$z$', r'I')
    Iv=(II-1.0)*(-1.0)
    Iinv=II.inv(1e-8)
    Ivinv=Iv.inv(1e-8)
    uSs_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UI", "VI", "WI"],r'$U_l$', r'$z$', r'$u_l$')*Iinv
    uSs_v=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UIv", "VIv", "WIv"],'U_v', r'$z$', r'$u_v$')*Ivinv
    Vr_Ss=uSs_v-uSs_l

    fstat = "Stats_singlephase.dt_ev"
    uSP=Field.LoadFromFile(fstat,["coordonnee_K"] , ["U", "V", "W"],r'$U_l$', r'$z$', r'$u_l$')

    uSb_l.settex(r'$\mathbf{Sb}$')
    uD_l.settex(r'$\mathbf{D}$')
    uSs_l .settex(r'$\mathbf{Ss}$')
    uB_l .settex(r'$\mathbf{B}$')
    uSP.settex(r'$\mathbf{SP}$')
    tracer([uSb_l, uD_l, uSs_l, uB_l, uSP], 'ul', [0] ,couleur=[5,2,1], markevry=[0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03], xlim=[0,1])
    tracer([uSb_l, uD_l, uSs_l, uB_l], 'ul_JFM', [0] ,couleur=[5,2,1], markevry=[0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03], xlim=[0,1])
    Vr_Sb.settex(r'$\mathbf{Sb}$')
    Vr_D.settex(r'$\mathbf{D}$')
    Vr_Ss.settex(r'$\mathbf{Ss}$')
    Vr_B.settex(r'$\mathbf{B}$')
    tracer([Vr_Sb, Vr_D, Vr_Ss, Vr_B], 'ur', [0] ,couleur=[4,2,1], markevry=[0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03], xlim=[0,1])

   ####################################################

    #sym_tens=[1,-1,1,-1,1,-1,1,-1,1]
    sym_tens=[1,1,-1,1,1,-1,-1,-1,1]
    xlims=[0,1]
    fstat = "Stats_sphere.dt_ev"
    I=Field.LoadFromFile(fstat,["coordonnee_K"] , ["I"],'I', r'$z$', r'I')
    Iv=(I-1.0)*(-1.0)
    Iinv=I.inv(1e-8)
    Ivinv=Iv.inv(1e-8)
    u_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UI", "VI", "WI"],r'$U_l$', r'$z$', r'$u_l$')*Iinv
    u_v=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UIv", "VIv", "WIv"],'U_v', r'$z$', r'$u_v$')*Ivinv
    uu_l=Iinv*Field.LoadFromFile(fstat,["coordonnee_K"] , ["UUI", "UVI", "UWI","UVI", "VVI", "VWI","UWI", "VWI", "WWI"],'UU_l', r'$z$', r'$\langle(u\times u)\rangle)_l$')
    uu_v=Ivinv*Field.LoadFromFile(fstat,["coordonnee_K"] , ["UUIv", "UVIv", "UWIv","UVIv", "VVIv", "VWIv","UWIv", "VWIv", "WWIv"],'UU_v', r'$z$', r'\langle(u\times u)\rangle)_v')

    ufuv=u_v.MatProduct(u_v)
    Rijv=uu_v-ufuv
    uful=u_l.MatProduct(u_l)
    Rijl=uu_l-uful
    Rij=Iv*Rijv+I*Rijl
    ## def de Jiacan
    #u_mono=u_l.mean()
    #u_mono=u_l*1.0
    u_mono=u_l*I+u_v*Iv

    uu_mono=I*uu_l+Iv*uu_v
    ufu_mono=u_mono.MatProduct(u_mono)
    Rij_mono=uu_mono-ufu_mono

    fstat = "Stats_deform.dt_ev"
    II=Field.LoadFromFile(fstat,["coordonnee_K"] , ["I"],'I', r'$z$', r'I')
    I2v=(II-1.0)*(-1.0)
    I2inv=II.inv(1e-8)
    I2vinv=I2v.inv(1e-8)
    u2_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UI", "VI", "WI"],r'$U_l$', r'$z$', r'$u_l$')*I2inv
    u2_v=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UIv", "VIv", "WIv"],'U_v', r'$z$', r'$u_v$')*I2vinv
    uu2_l=I2inv*Field.LoadFromFile(fstat,["coordonnee_K"] , ["UUI", "UVI", "UWI","UVI", "VVI", "VWI","UWI", "VWI", "WWI"],'UU_l', r'$z$', r'$\langle(u\times u)\rangle)_l$')
    uu2_v=I2vinv*Field.LoadFromFile(fstat,["coordonnee_K"] , ["UUIv", "UVIv", "UWIv","UVIv", "VVIv", "VWIv","UWIv", "VWIv", "WWIv"],'UU_v', r'$z$', r'\langle(u\times u)\rangle)_v')

    ufuv=u2_v.MatProduct(u2_v)
    R2ijv=uu2_v-ufuv
    uful=u2_l.MatProduct(u2_l)
    R2ijl=uu2_l-uful
    R2ij=I2v*R2ijv+II*R2ijl

    u_mono=u2_l*II+u2_v*I2v
    uu_mono=II*uu2_l+I2v*uu2_v
    ufu_mono=u_mono.MatProduct(u_mono)
    R2ij_mono=uu_mono-ufu_mono



    fstat = "Stats_singlephase.dt_ev"
    u=Field.LoadFromFile(fstat,["coordonnee_K"] , ["U", "V", "W"],r'$U_l$', r'$z$', r'$u_l$')
    uu=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UU", "UV", "UW","UV", "VV", "VW","UW", "VW", "WW"],'UU_l', r'$z$', r'$\langle(u\times u)\rangle)_l$')

    ufu=u.MatProduct(u)
    R5ijl=uu-ufu

    fstat = "Stats_bidi127.dt_ev"
    I3=Field.LoadFromFile(fstat,["coordonnee_K"] , ["I"],'I', r'$z$', r'I')
    I3v=(I3-1.0)*(-1.0)
    I3inv=I3.inv(1e-8)
    I3vinv=I3v.inv(1e-8)
    u3_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UI", "VI", "WI"],r'$U_l$', r'$z$', r'$u_l$')*I3inv
    u3_v=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UIv", "VIv", "WIv"],'U_v', r'$z$', r'$u_v$')*I3vinv
    uu3_l=I3inv*Field.LoadFromFile(fstat,["coordonnee_K"] , ["UUI", "UVI", "UWI","UVI", "VVI", "VWI","UWI", "VWI", "WWI"],'UU_l', r'$z$', r'$\langle(u\times u)\rangle)_l$')
    uu3_v=I3vinv*Field.LoadFromFile(fstat,["coordonnee_K"] , ["UUIv", "UVIv", "UWIv","UVIv", "VVIv", "VWIv","UWIv", "VWIv", "WWIv"],'UU_v', r'$z$', r'\langle(u\times u)\rangle)_v')
    
    ufuv=u3_v.MatProduct(u3_v)
    R3ijv=uu3_v-ufuv
    uful=u3_l.MatProduct(u3_l)
    R3ijl=uu3_l-uful
    R3ij=I3v*R3ijv+I3*R3ijl
    R3ij.settex(r'bidisperse (TrioCFD)')
    u_mono=u3_l*I3+u3_v*I3v
    uu_mono=I3*uu3_l+I3v*uu3_v
    ufu_mono=u_mono.MatProduct(u_mono)
    R3ij_mono=uu_mono-ufu_mono


    fstat = "Stats_petitesphere.dt_ev"
    I4=Field.LoadFromFile(fstat,["coordonnee_K"] , ["I"],'I', r'$z$', r'I')
    I4v=(I4-1.0)*(-1.0)
    I4inv=I4.inv(1e-8)
    I4vinv=I4v.inv(1e-8)
    u4_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UI", "VI", "WI"],r'$U_l$', r'$z$', r'$u_l$')*I4inv
    u4_v=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UIv", "VIv", "WIv"],'U_v', r'$z$', r'$u_v$')*I4vinv
    uu4_l=I4inv*Field.LoadFromFile(fstat,["coordonnee_K"] , ["UUI", "UVI", "UWI","UVI", "VVI", "VWI","UWI", "VWI", "WWI"],'UU_l', r'$z$', r'$\langle(u\times u)\rangle)_l$')
    uu4_v=I4vinv*Field.LoadFromFile(fstat,["coordonnee_K"] , ["UUIv", "UVIv", "UWIv","UVIv", "VVIv", "VWIv","UWIv", "VWIv", "WWIv"],'UU_v', r'$z$', r'\langle(u\times u)\rangle)_v')
    
    ufuv=u4_v.MatProduct(u4_v)
    R4ijv=uu4_v-ufuv
    uful=u4_l.MatProduct(u4_l)
    R4ijl=uu4_l-uful
    R4ij=I4v*R4ijv+I4*R4ijl
    R4ij.settex(r'petite sphérique (TrioCFD)')
    u_mono=u4_l*I4+u4_v*I4v
    uu_mono=I4*uu4_l+I4v*uu4_v
    ufu_mono=u_mono.MatProduct(u_mono)
    R4ij_mono=uu_mono-ufu_mono


    # Set names : 
    u_l.settex(r'$u_{liquide}$ spherical (TrioCFD)')
    u_v.settex(r'$u_{vapeur}$ spherical (TrioCFD)')
    u2_l.settex(r'$u_{liquide}$ deformable (TrioCFD)')
    u2_v.settex(r'$u_{vapeur}$ deformable (TrioCFD)')
    u3_l.settex(r'$u_{liquide}$ deformable (TrioCFD)')
    u3_l.settex(r'$u_{vapeur}$ deformable (TrioCFD)')
    I=I.add(1.0, '--')
    II=II.add(1.0, '--')
    I3=I3.add(1.0, '--')
    I4=I4.add(1.0, '--')
    I_tryg=Field.LoadFromFile('I_sphere_TRYG.dt_ev',["coordonnee_K"] , ["I"],'I', r'$z$', r'I')
    I_tryg.settex(r'Sb (LT2008)')
    II_tryg=Field.LoadFromFile('I_defo_TRYG.dt_ev',["coordonnee_K"] , ["I"],'I', r'$z$', r'I')
    II_tryg.settex(r'D (LT2008)')
    I.settex(r'Sb (TrioCFD)')
    I4.settex(r'Ss (TrioCFD)')
    II.settex(r'D (TrioCFD)')
    I3.settex(r'B (TrioCFD)')
    I.symetriser()
    II.symetriser()
    I3.symetriser()
    I4.symetriser()

    IssD=(I4).ChangerMaillage(II)+II
    tracer([I, I_tryg,II, II_tryg,I3],'Void_fraction',textitle=r'Void fraction',couleur=[2,2,1], xlim=xlims, markevry=[1])
    tracer([I, I_tryg,II, II_tryg],'Void_fraction_compa',textitle=r'Void fraction',couleur=[2,2,1], xlim=xlims, markevry=[1])
    tracer([I,I4, II_tryg,I3],'Void_fraction_bidi',textitle=r'Void fraction',couleur=[1,2,1], xlim=xlims, markevry=[1])
    #tracerJFM([I,II, I4,I3, I_tryg, II_tryg, IssD])
    #tracer([I,II, I4,I3, I_tryg, II_tryg],'Void_fraction_all',textitle=r'Void fraction',couleur=[4,2,1], markevry=[0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03], xlim=[0,1])  
    


    ############################################################################################
    ############################## PLOT DES U_RMS ##############################################
    ############################################################################################
    
    # Adimensionnement par utau

    Rij=Rij*(47.24/2.)*(47.24/2.)
    R2ij=R2ij*(47.24/2.)*(47.24/2.)
    R3ij=R3ij*(47.24/2.)*(47.24/2.)
    R4ij=R4ij*(47.24/2.)*(47.24/2.)
    Rij_mono=Rij_mono*(47.24/2.)*(47.24/2.)
    R2ij_mono=R2ij_mono*(47.24/2.)*(47.24/2.)
    R3ij_mono=R3ij_mono*(47.24/2.)*(47.24/2.)
    R4ij_mono=R4ij_mono*(47.24/2.)*(47.24/2.)
    Rijv=Rijv*(47.24/2.)*(47.24/2.)
    R2ijv=R2ijv*(47.24/2.)*(47.24/2.)
    R3ijv=R3ijv*(47.24/2.)*(47.24/2.)
    R4ijv=R4ijv*(47.24/2.)*(47.24/2.)
    Rijl=Rijl*(47.24/2.)*(47.24/2.)
    R2ijl=R2ijl*(47.24/2.)*(47.24/2.)
    R3ijl=R3ijl*(47.24/2.)*(47.24/2.)
    R4ijl=R4ijl*(47.24/2.)*(47.24/2.)
    R5ijl=R5ijl*(47.24/2.)*(47.24/2.)


    aa=Rij*1.0
    Rij.symetriser(sym_tens)
    R2ij.symetriser(sym_tens)
    R3ij.symetriser(sym_tens)
    R4ij.symetriser(sym_tens)
    Rij_mono.symetriser(sym_tens)
    R2ij_mono.symetriser(sym_tens)
    R3ij_mono.symetriser(sym_tens)
    R4ij_mono.symetriser(sym_tens)

    Rijl.symetriser(sym_tens)
    R2ijl.symetriser(sym_tens)
    R3ijl.symetriser(sym_tens)
    R4ijl.symetriser(sym_tens)
    R5ijl.symetriser(sym_tens)
    Rijv.symetriser(sym_tens)
    R2ijv.symetriser(sym_tens)
    R3ijv.symetriser(sym_tens)
    R4ijv.symetriser(sym_tens)

    aa.settex(r'aa (TrioCFD)')
    Rij.settex(r'Sb (TrioCFD)')
    R2ij.settex(r'D (TrioCFD)')
    R3ij.settex(r'B')
    R4ij.settex(r'Ss')
    Rij_mono.settex(r'Sb (TrioCFD)')
    R2ij_mono.settex(r'D (TrioCFD)')
    R3ij_mono.settex(r'B')
    R4ij_mono.settex(r'Ss')
    Rijv.settex(r'Sb (TrioCFD)')
    R2ijv.settex(r'D (TrioCFD)')
    R3ijv.settex(r'B')
    R4ijv.settex(r'Ss')
    Rijl.settex(r'Sb (TrioCFD)')
    R2ijl.settex(r'D (TrioCFD)')
    R3ijl.settex(r'B')
    R4ijl.settex(r'Ss')
    R5ijl.settex(r'SP')



# Comparaison a la référence
    Rij_tryg00=Field.LoadFromFile('UU_sphere_TRYG.dt_ev',["coordonnee_K"] , ["I"],'I', r'$z$', r'I')
    R2ij_tryg00=Field.LoadFromFile('UU_defo_TRYG.dt_ev',["coordonnee_K"] , ["I"],'I', r'$z$', r'I')
    Rij_tryg00=Rij_tryg00.power(2)
    R2ij_tryg00=R2ij_tryg00.power(2)
    Rij_tryg00.settex(r'Sb (LT2008)')
    R2ij_tryg00.settex(r'D (LT2008)')

    Rij_tryg22=Field.LoadFromFile('VV_sphere_TRYG.dt_ev',["coordonnee_K"] , ["I"],'I', r'$z$', r'I')
    R2ij_tryg22=Field.LoadFromFile('VV_defo_TRYG.dt_ev',["coordonnee_K"] , ["I"],'I', r'$z$', r'I')
    Rij_tryg22=Rij_tryg22.power(2)
    R2ij_tryg22=R2ij_tryg22.power(2)
    Rij_tryg22.settex(r'Sb (LT2008)')
    R2ij_tryg22.settex(r'D (LT2008)')

    Rij_tryg20=Field.LoadFromFile('UV_sphere_TRYG.dt_ev',["coordonnee_K"] , ["I"],'I', r'$z$', r'I')
    R2ij_tryg20=Field.LoadFromFile('UV_defo_TRYG.dt_ev',["coordonnee_K"] , ["I"],'I', 'z', r'I')
    Rij_tryg20.settex(r'Sb (LT2008)')
    R2ij_tryg20.settex(r'D (LT2008)')
    

    tracer([R2ijl.compo(0,0), R2ij.compo(0,0), R2ij_mono.compo(0,0), R2ij_tryg00, R3ijl.compo(0,0), R3ij.compo(0,0), R3ij_mono.compo(0,0), Rijl.compo(0,0), Rij.compo(0,0), Rij_mono.compo(0,0), Rij_tryg00, R4ijl.compo(0,0), R4ij.compo(0,0), R4ij_mono.compo(0,0)], 'Rij_v_l_vl_11', couleur=[4,3,4,3,3], xlim=[0,1], doubleLegend=True, casename=["$\overline{\chi_l\mathbf{u_l^,u_l^,}}$", "$\overline{\chi_l\mathbf{u_l^,u_l^,}}+\overline{\chi_v\mathbf{u_v^,u_v^,}}$", "$\overline{\mathbf{u^,u^,}}$ ", 'ref'], markevry=[0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03])

    #tracerJFM([Rij_mono.compo(0,0), R2ijl.compo(0,0), R4ij_mono.compo(0,0),  R3ij_mono.compo(0,0), R5ijl.compo(0,0), Rij_tryg00, R2ij_tryg00])
    tracer([Rij_mono.compo(0,0), R2ijl.compo(0,0), R4ij_mono.compo(0,0),  R3ij_mono.compo(0,0), R5ijl.compo(0,0), Rij_tryg00, R2ij_tryg00], 'R11_all', couleur=[5,2], xlim=[0,1], doubleLegend=False, markevry=[0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03])


    if 1 :

	color2="#8dd3c7"
	color4="#ffffb3"
	color1="#bebada"
	color3="#fb8072"

        ###### figure 1
	scalepng=1.5
	f, (ax, ax1, ax2) = plt.subplots(figsize=(4, 4), ncols=3)


        x=R3ij_mono.compo(0,0)._ax[0,:]
	y1=R3ij_mono.compo(0,0)._npa[0,:]
	y2=R4ij_mono.compo(0,0)._npa[0,:]
        zero=y2*0.


        ax.set_xlim([0.,1.0], emit=False)
        ax.set_ylim([0.,3.5], emit=False)
	ax.set_xlabel(r"${y/h}$", fontsize=8)
	ax.set_ylabel(r"$\overline{uu}/u_\tau^2$", fontsize=8)
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.xaxis.set_ticks_position('bottom')
	ax.yaxis.set_ticks_position('left') 
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.tight_layout()
	dpii=400
	
        ax.plot(x, y1, linestyle='-', linewidth=0.4, marker=u'.', label='B', color='black', ms=4., markevery=0.03, markerfacecolor='white', markeredgecolor = 'black', markeredgewidth=1.0)
        ax.plot(x, y2, linestyle='-', linewidth=0.4, marker=u'.', label='Ss', color='black', ms=scalepng*1.5, markevery=0.03, markerfacecolor='black', markeredgecolor = 'black',markeredgewidth=1.0)
        ax.plot(x, zero, color='black')
        ax.fill_between(x, y1, y2, where=y1 >= y2, facecolor=color1, edgecolor=color1, interpolate=True)
        ax.fill_between(x, y2, zero, where=y2 >= zero, facecolor=color2, edgecolor=color2, interpolate=True)


	#f.set_size_inches(scalepng,scalepng)
	#box = ax.get_position()
	#ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
	#ax.legend(loc='center left', bbox_to_anchor=(1.1, 0.5), frameon=True)
	#inSizeLegend=int(scalepng*2.0)
	#plt.legend(loc=0,prop={'size':inSizeLegend}, frameon=True)
	#plt.savefig("blabla", bbox_inches='tight', dpi=dpii)

        ###### figure 2
        R3ij_mono_11=(R3ij_mono.compo(0,0)).ChangerMaillage(R2ij_tryg00)
        R4ij_mono_11=(R4ij_mono.compo(0,0)).ChangerMaillage(R2ij_tryg00)
        R5ijl_11=(R5ijl.compo(0,0)).ChangerMaillage(R2ij_tryg00)
        #f=plt.figure(1)
	#ax = plt.subplot(111)
	y3=R3ij_mono_11._npa[0,:]-R4ij_mono_11._npa[0,:]
        x=R2ij_tryg00._ax[0,:]
	y4=R2ij_tryg00._npa[0,:]
	zero=y3*0.
	

        ax1.set_xlim([0.,1.0], emit=False)
        ax1.set_ylim([0.,8.0], emit=False)
	ax1.set_xlabel(r"${y/h}$", fontsize=8)
	#ax1.set_ylabel(r"$\overline{u'u'}/u_\tau^2$", fontsize=8)
	ax1.spines['top'].set_visible(False)
	ax1.spines['right'].set_visible(False)
	ax1.xaxis.set_ticks_position('bottom')
	ax1.yaxis.set_ticks_position('left') 
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.tight_layout()
	dpii=400
	
        ax1.plot(x, y3, linestyle='-', linewidth=0.4, marker=u'+', label='B-Ss', color='black', ms=scalepng*1.5, markevery=0.03, markerfacecolor='white', markeredgecolor = 'black', markeredgewidth=1.0)
        ax1.plot(x, y4, linestyle='-', linewidth=0.4, marker=u'd', label='D', color='black', ms=scalepng*1.5, markevery=0.03, markerfacecolor='black', markeredgecolor = 'black',markeredgewidth=1.0)
        ax1.plot(x, zero, color='black')
        ax1.fill_between(x, y3, y4, where=y4 >= y3, facecolor=color3, edgecolor=color3, interpolate=True)
        ax1.fill_between(x, y3, zero, where=y3 >= zero, facecolor=color1, edgecolor=color1, interpolate=True)


	#f.set_size_inches(scalepng,scalepng)
	#box = ax.get_position()
	#ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
	#ax.legend(loc='center left', bbox_to_anchor=(1.1, 0.5), frameon=True)
	#inSizeLegend=int(scalepng*2.0)
	#plt.legend(loc=0,prop={'size':inSizeLegend}, frameon=True)
	#plt.savefig("zaroufl", bbox_inches='tight', dpi=dpii)


        ###### figure 3


        #f=plt.figure(2)
	#ax = plt.subplot(111)
	

	y5=y4-y3
        y6=R5ijl_11._npa[0,:]
	zero=y5*0.
	




        ax2.set_xlim([0.,1.0], emit=False)
        ax2.set_ylim([0.,8.0], emit=False)
	ax2.set_xlabel(r"${y/h}$", fontsize=8)
	#ax2.set_ylabel(r"$\overline{u'u'}/u_\tau^2$", fontsize=8)
	ax2.spines['top'].set_visible(False)
	ax2.spines['right'].set_visible(False)
	ax2.xaxis.set_ticks_position('bottom')
	ax2.yaxis.set_ticks_position('left') 
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.tight_layout()
	dpii=400
	
        ax2.plot(x, y5, linestyle='-', linewidth=0.4, marker=u'x', label='D-(B-Ss)', color='black', ms=scalepng*1.5, markevery=0.03, markerfacecolor='white', markeredgecolor = 'black', markeredgewidth=1.0)
        ax2.plot(x, y6, linestyle='-', linewidth=0.4, marker=None, label='SP', color='black', ms=scalepng*1.5, markevery=0.03, markerfacecolor='black', markeredgecolor = 'black',markeredgewidth=1.0)
        ax2.plot(x, zero, color='black')
        ax2.fill_between(x, y5, y6, where=y6 >= y5, facecolor=color4, edgecolor=color4, interpolate=True)
        ax2.fill_between(x, y5, zero, where=y5 >= zero, facecolor=color3, edgecolor=color3, interpolate=True)

	plt.text(-2.,2.5, r"$\mathbf{R_{ij,B}^{BIF}}$", fontsize=8)
	plt.text(-0.8,0.7, r"$\mathbf{R_{ij,D}^{BIF}}$", fontsize=8)
	plt.text(-2.38,0.2, r"$\mathbf{R_{ij,Ss}^{SPT}+R_{ij,Ss}^{bif}}$", fontsize=6)
	plt.text(-1.1,3.0, r"$\mathbf{R_{ij,D}^{SPT}}$", fontsize=8)
	plt.text(0.1,1., r"$\mathbf{R_{ij,D}^{SPT}}$", fontsize=8)
	plt.text(0.7,1., "interaction", fontsize=5)
	plt.text(0.8,0.5, "SPT/BIF", fontsize=5)


	f.set_size_inches(scalepng*2.0,scalepng)
	box = ax.get_position()
	#ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
	ax.legend(loc="upper right")
	#ax1.set_position([box.x0, box.y0, box.width * 0.8, box.height])
	ax1.legend(loc="upper right")
	#ax2.set_position([box.x0, box.y0, box.width * 0.8, box.height])
	ax2.legend(loc="upper right")
	inSizeLegend=int(scalepng*2.0)
	#plt.legend(loc=0,prop={'size':inSizeLegend}, frameon=True)
	plt.savefig("schema_interaction", bbox_inches='tight', dpi=dpii)





	#plt.plot(Rij_mono.compo(0,0)._ax[0,:], Rij_mono.compo(0,0)._npa[0,:], linestyle='-', linewidth=0.4, marker=None, label='SP', color='k', ms=scalepng*1.5, markevery=0.03, markerfacecolor='k', markeredgecolor = 'k', markeredgewidth=1.0)


	#plt.semilogy(pdf_mono_bulk._ax[0,155:],pdf_mono_bulk._npa[0,155:], linestyle='-', linewidth=0.4, marker=None, label='SP', color='k', ms=scalepng*1.5, markevery=0.03, markerfacecolor='k', markeredgecolor = 'k', markeredgewidth=1.0)
	#plt.semilogy(pdf_bidi_bulk._ax[0,30:],pdf_bidi_bulk._npa[0,30:], linewidth=0.4, marker='.', label='B', color='b', ms=4., markevery=0.05, markerfacecolor='white', markeredgecolor = 'b', markeredgewidth=1.0)
	#plt.semilogy(pdf_defo_bulk2._ax[0,118:],pdf_defo_bulk2._npa[0,118:], linestyle='-', linewidth=0.4, marker='d', label='D', color='r', ms=scalepng, markevery=0.05, markerfacecolor='k', markeredgecolor = 'r', markeredgewidth=1.0)





    #tracerJFM([Rij_mono.compo(2,2), R2ijl.compo(2,2), R4ij_mono.compo(2,2),  R3ij_mono.compo(2,2), R5ijl.compo(2,2), Rij_tryg22, R2ij_tryg22])
    tracer([Rij_mono.compo(2,2), R2ijl.compo(2,2), R4ij_mono.compo(2,2),  R3ij_mono.compo(2,2), R5ijl.compo(2,2), Rij_tryg22, R2ij_tryg22], 'R22_all', couleur=[5,2], xlim=[0,1], doubleLegend=False, markevry=[0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03])


    R3ij_mono=R3ij_mono.ChangerMaillage(R2ijl)
    R4ij_mono=R4ij_mono.ChangerMaillage(R2ijl)
    Rij_BIF_D=R2ijl-(R3ij_mono-R4ij_mono)
    Rij_BIF_B=R5ijl
    Rij_BIF_D.settex(r'$\mathbf{D-B+Ss}$')
    Rij_BIF_B.settex(r'$\mathbf{SP}$')
    tracer([Rij_BIF_B, Rij_BIF_D], 'R11_BIT', [0], couleur=[2,2], xlim=[0,1], doubleLegend=False, markevry=[0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03])
    tracer([Rij_BIF_B, Rij_BIF_D], 'R22_BIT', [4], couleur=[2,2], xlim=[0,1], doubleLegend=False, markevry=[0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03])
    tracer([Rij_BIF_B, Rij_BIF_D], 'R33_BIT', [8], couleur=[2,2], xlim=[0,1], doubleLegend=False, markevry=[0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03])
    #Rij_BIF_B=Rij_BIF_B*(-1.0)
    #Rij_BIF_D=Rij_BIF_D*(-1.0)
    #Rij_BIF_D.settex(r'$\mathbf{D-B+Ss}$')
    #Rij_BIF_B.settex(r'$\mathbf{SP}$')
    tracer([Rij_BIF_B, Rij_BIF_D], 'R13_BIT', [2], couleur=[2,2], xlim=[0,1], doubleLegend=False, markevry=[0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03])




    R3ij_mono_11=(R3ij_mono.compo(0,0)).ChangerMaillage(R2ij_tryg00)
    R4ij_mono_11=(R4ij_mono.compo(0,0)).ChangerMaillage(R2ij_tryg00)

    Rij_BIF_D=R2ij_tryg00-(R3ij_mono_11-R4ij_mono_11)
    Rij_BIF_D.settex(r'$\mathbf{D-B+Ss}$')
    tracer([Rij_BIF_B.compo(0,0), Rij_BIF_D], 'R11_BIT_other', couleur=[2,2], xlim=[0,1], doubleLegend=False, markevry=      [0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03])



    Rij_mono=Rij_mono*(-1.0)
    R2ijl=R2ijl*(-1.0)
    R4ij_mono=R4ij_mono*(-1.0)
    R3ij_mono=R3ij_mono*(-1.0)
    R5ijl=R5ijl*(-1.0)
    #tracerJFM([Rij_mono.compo(0,2), R2ijl.compo(0,2), R4ij_mono.compo(0,2),  R3ij_mono.compo(0,2), R5ijl.compo(0,2), Rij_tryg20, R2ij_tryg20])
    tracer([Rij_mono.compo(0,2), R2ijl.compo(0,2), R4ij_mono.compo(0,2),  R3ij_mono.compo(0,2), R5ijl.compo(0,2), Rij_tryg20, R2ij_tryg20], 'R20_all', couleur=[5,2], xlim=[0,1], doubleLegend=False, markevry=[0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03])

    Rij_BIF_D=R2ijl-R5ijl
    Rij_BIF_B=R3ij_mono-R4ij_mono
    Rij_BIF_D.settex(r'$\mathbf{D-SP}$')
    Rij_BIF_B.settex(r'$\mathbf{B-Ss}$')
    tracer([Rij_BIF_D, Rij_BIF_B], 'R13_BIT_2', [2], couleur=[2,2,2,2,2], xlim=[0,1], doubleLegend=False, markevry=[0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03])



    ### pour tryggvason
   # tracer([R2ij.compo(0,0), R2ij_tryg00, Rij.compo(0,0), Rij_tryg00], 'R11_bifluid', couleur=[2,2,1], xlim=[0,1], markevry=[0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03])
   # tracer([R2ij_mono.compo(0,0), R2ij_tryg00, Rij_mono.compo(0,0), Rij_tryg00], 'R11_mono', couleur=[2,2,1], xlim=[0,1], markevry=[0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03])
   # tracer([R2ijl.compo(0,0), R2ij_tryg00, Rijl.compo(0,0), Rij_tryg00], 'R11_liquide', couleur=[2,2,1], xlim=[0,1], markevry=[0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03])

   # tracer([R2ij.compo(2,2), R2ij_tryg22, Rij.compo(2,2), Rij_tryg22], 'R22_bifluid', couleur=[2,2,1], xlim=[0,1], markevry=[0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03])
   # tracer([R2ij_mono.compo(2,2), R2ij_tryg22, Rij_mono.compo(2,2), Rij_tryg22], 'R22_mono', couleur=[2,2,1], xlim=[0,1], markevry=[0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03])
   # tracer([R2ijl.compo(2,2), R2ij_tryg22, Rijl.compo(2,2), Rij_tryg22], 'R22_liquide', couleur=[2,2,1], xlim=[0,1], markevry=[0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03])
	

    tracer([R2ijl.compo(2,2), R2ij.compo(2,2), R2ij_mono.compo(2,2), R2ij_tryg22, R3ijl.compo(2,2), R3ij.compo(2,2), R3ij_mono.compo(2,2), Rijl.compo(2,2), Rij.compo(2,2), Rij_mono.compo(2,2), Rij_tryg22, R4ijl.compo(2,2), R4ij.compo(2,2), R4ij_mono.compo(2,2)], 'Rij_v_l_vl_22', couleur=[4,3,4,3,3], xlim=[0,1], doubleLegend=True, casename=["$\overline{\chi_l\mathbf{u_l^,u_l^,}}$", "$\overline{\chi_l\mathbf{u_l^,u_l^,}}+\overline{\chi_v\mathbf{u_v^,u_v^,}}$", "$\overline{\mathbf{u^,u^,}}$ ", 'ref'], markevry=[0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03])

    tracer([R2ijl.compo(2,0), R2ij.compo(2,0), R2ij_mono.compo(2,0), R2ij_tryg20, R3ijl.compo(2,0), R3ij.compo(2,0), R3ij_mono.compo(2,0), Rijl.compo(2,0), Rij.compo(2,0), Rij_mono.compo(2,0), Rij_tryg20, R4ijl.compo(2,0), R4ij.compo(2,0), R4ij_mono.compo(2,0)], 'Rij_v_l_vl_20', couleur=[4,3,4,3,3], xlim=[0,1], doubleLegend=True, casename=["$\overline{\chi_l\mathbf{u_l^,u_l^,}}$", "$\overline{\chi_l\mathbf{u_l^,u_l^,}}+\overline{\chi_v\mathbf{u_v^,u_v^,}}$", "$\overline{\mathbf{u^,u^,}}$ ", 'ref'], markevry=[0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03])

    tracer([R2ijl.compo(1,1), R2ij.compo(1,1), R3ijl.compo(1,1), R3ij.compo(1,1), Rijl.compo(1,1), Rij.compo(1,1), R4ijl.compo(1,1), R4ij.compo(1,1)], 'Rij_v_l_vl_33', couleur=[2,2,2,2,2], xlim=[0,1], doubleLegend=True, casename=["$\overline{\mathbf{u_l^,u_l^,}}^l$", "$\overline{\mathbf{u^,u^,}}$", 'ref'], markevry=[0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03])



    #for i in range(2):Rij, Rijv, Rijl,
    #    for j in range(2):
    #        Rij.compo(i,j).settex(r'spherical (TrioCFD)')
    #        R2ij.compo(i,j).settex(r'deformable (TrioCFD)')
    #        aa.compo(i,j).settex(r'deformable (TrioCFD)')
    #        pass
    #    pass
    


    



    return


def UsualVariable() :

    """ 
    Comparaison de différents cas (Ss127, D127) avec la littérature (lu & tryggvason) codée en dur. Ne pas appeler cette fonction sur de nouveaux cas, elle ne fonctionnera pas.
    """ 

    ################################################################################
    ########### PLOT DES PROFIL DE VITESSE ET INDICATRICE ##########################
    ################################################################################
    #sym_tens=[1,-1,1,-1,1,-1,1,-1,1]
    sym_tens=[1,1,-1,1,1,-1,-1,-1,1]
    xlims=[0,1]
    fstat = "Stats_sphere.dt_ev"
    I=Field.LoadFromFile(fstat,["coordonnee_K"] , ["I"],'I', r'$z$', r'I')
    u_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UI", "VI", "WI"],r'$U_l$', r'$z$', r'$u_l$')
    u_v=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UIv", "VIv", "WIv"],'U_v', r'$z$', r'$u_v$')
    u=u_l+u_v
    uu_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UUI", "UVI", "UWI","UVI", "VVI", "VWI","UWI", "VWI", "WWI"],'UU_l', r'$z$', r'$\langle(u\times u)\rangle)_l$')
    uu_v=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UUIv", "UVIv", "UWIv","UVIv", "VVIv", "VWIv","UWIv", "VWIv", "WWIv"],'UU_v', r'$z$', r'\langle(u\times u)\rangle)_v')
    uu=uu_v+uu_l
    ufu=u.MatProduct(u)
    Rij=uu-ufu
    fstat = "Stats_deform.dt_ev"
    II=Field.LoadFromFile(fstat,["coordonnee_K"] , ["I"],'I', r'$z$', r'I')
    u2_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UI", "VI", "WI"],r'$U_l$', r'$z$', r'$u_l$')
    u2_v=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UIv", "VIv", "WIv"],'U_v', r'$z$', r'$u_v$')
    u2=u2_l+u2_v
    uu2_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UUI", "UVI", "UWI","UVI", "VVI", "VWI","UWI", "VWI", "WWI"],'UU_l', r'$z$', r'$\langle(u\times u)\rangle)_l$')
    uu2_v=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UUIv", "UVIv", "UWIv","UVIv", "VVIv", "VWIv","UWIv", "VWIv", "WWIv"],'UU_v', r'$z$', r'\langle(u\times u)\rangle)_v')
    uu2=uu2_v+uu2_l
    ufu2=u2.MatProduct(u2)
    R2ij=uu2-ufu2
    fstat = "Stats_bidi127.dt_ev"
    I3=Field.LoadFromFile(fstat,["coordonnee_K"] , ["I"],'I', r'$z$', r'I')
    u3_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UI", "VI", "WI"],r'$U_l$', r'$z$', r'$u_l$')
    u3_v=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UIv", "VIv", "WIv"],'U_v', r'$z$', r'$u_v$')
    u3=u3_l+u3_v
    uu3_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UUI", "UVI", "UWI","UVI", "VVI", "VWI","UWI", "VWI", "WWI"],'UU_l', r'$z$', r'$\langle(u\times u)\rangle)_l$')
    uu3_v=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UUIv", "UVIv", "UWIv","UVIv", "VVIv", "VWIv","UWIv", "VWIv", "WWIv"],'UU_v', r'$z$', r'\langle(u\times u)\rangle)_v')
    uu3=uu3_v+uu3_l
    ufu3=u3.MatProduct(u3)
    R3ij=uu3-ufu3
    R3ij.settex(r'bidisperse (TrioCFD)')

    fstat = "Stats_petitesphere.dt_ev"
    I4=Field.LoadFromFile(fstat,["coordonnee_K"] , ["I"],'I', r'$z$', r'I')
    u4_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UI", "VI", "WI"],r'$U_l$', r'$z$', r'$u_l$')
    u4_v=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UIv", "VIv", "WIv"],'U_v', r'$z$', r'$u_v$')
    u4=u4_l+u4_v
    uu4_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UUI", "UVI", "UWI","UVI", "VVI", "VWI","UWI", "VWI", "WWI"],'UU_l', r'$z$', r'$\langle(u\times u)\rangle)_l$')
    uu4_v=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UUIv", "UVIv", "UWIv","UVIv", "VVIv", "VWIv","UWIv", "VWIv", "WWIv"],'UU_v', r'$z$', r'\langle(u\times u)\rangle)_v')
    uu4=uu4_v+uu4_l
    ufu4=u4.MatProduct(u4)
    R4ij=uu4-ufu4
    R4ij.settex(r'petite sphérique (TrioCFD)')

    # Set names : 
    u_l.settex(r'$u_{liquide}$ spherical (TrioCFD)')
    u_v.settex(r'$u_{vapeur}$ spherical (TrioCFD)')
    u2_l.settex(r'$u_{liquide}$ deformable (TrioCFD)')
    u2_v.settex(r'$u_{vapeur}$ deformable (TrioCFD)')
    u3_l.settex(r'$u_{liquide}$ deformable (TrioCFD)')
    u3_l.settex(r'$u_{vapeur}$ deformable (TrioCFD)')
    I=I.add(1.0, '--')
    II=II.add(1.0, '--')
    I3=I3.add(1.0, '--')
    I4=I4.add(1.0, '--')
    I_tryg=Field.LoadFromFile('I_sphere_TRYG.dt_ev',["coordonnee_K"] , ["I"],'I', r'$z$', r'I')
    I_tryg.settex(r'Sb $\mathbf{d_b=0.3}$ (LT2008)')
    II_tryg=Field.LoadFromFile('I_defo_TRYG.dt_ev',["coordonnee_K"] , ["I"],'I', r'$z$', r'I')
    II_tryg.settex(r'D (LT2008)')
    I.settex(r'Sb $\mathbf{d_b=0.3}$ (TrioCFD)')
    I4.settex(r'Ss $\mathbf{d_b=0.165}$ (TrioCFD)')
    II.settex(r'D (TrioCFD)')
    I3.settex(r'B (TrioCFD)')
    I.symetriser()
    II.symetriser()
    I3.symetriser()
    I4.symetriser()
    tracer([I, I_tryg,II, II_tryg,I3],'Void_fraction',textitle=r'Void fraction',couleur=[2,2,1], xlim=xlims, markevry=[1])
    tracer([I, I_tryg,II, II_tryg],'Void_fraction_compa',textitle=r'Void fraction',couleur=[2,2,1], xlim=xlims, markevry=[1])
    tracer([I,I4, II_tryg,I3],'Void_fraction_bidi',textitle=r'Void fraction',couleur=[1,2,1], xlim=xlims, markevry=[1])
    tracer([I, I_tryg, I4,II, II_tryg,I3],'Void_fraction_all',textitle=r'Void fraction',couleur=[3,2,1], xlim=xlims, markevry=[1,1,2,1,1,1,1,1,1,1,1])    
    ############################################################################################
    ############################## PLOT DES U_RMS ##############################################
    ############################################################################################
    
    # Adimensionnement par utau

    Rij=Rij*(47.24/2.)*(47.24/2.)
    R2ij=R2ij*(47.24/2.)*(47.24/2.)
    R3ij=R3ij*(47.24/2.)*(47.24/2.)
    R4ij=R4ij*(47.24/2.)*(47.24/2.)
    aa=Rij*1.0
    #Rij.symetriser(sym_tens)
    #R2ij.symetriser(sym_tens)
    #R3ij.symetriser(sym_tens)
    aa.settex(r'aa (TrioCFD)')
    Rij.settex(r'Sb : $\mathbf{d_b=0.3}$ (TrioCFD)')
    R2ij.settex(r'D : (TrioCFD)')
    R3ij.settex(r'B : (TrioCFD)')
    R4ij.settex(r'Ss : $\mathbf{d_b=0.165}$ (TrioCFD)')
    #for i in range(2):
    #    for j in range(2):
    #        Rij.compo(i,j).settex(r'spherical (TrioCFD)')
    #        R2ij.compo(i,j).settex(r'deformable (TrioCFD)')
    #        aa.compo(i,j).settex(r'deformable (TrioCFD)')
    #        pass
    #    pass
    


    # Comparaison a la référence
    Rij_tryg00=Field.LoadFromFile('UU_sphere_TRYG.dt_ev',["coordonnee_K"] , ["I"],'I', r'$z$', r'I')
    R2ij_tryg00=Field.LoadFromFile('UU_defo_TRYG.dt_ev',["coordonnee_K"] , ["I"],'I', r'$z$', r'I')
    Rij_tryg00=Rij_tryg00.power(2)
    R2ij_tryg00=R2ij_tryg00.power(2)
    Rij_tryg00.settex(r'Sb : (LT2008)')
    R2ij_tryg00.settex(r'D : (LT2008)')
    tracer([Rij.compo(0,0), Rij_tryg00 , R2ij.compo(0,0), R2ij_tryg00, R3ij.compo(0,0)], 'uu',textitle=r"Axial $\overline{u'u'}$", couleur=[2,2,1], xlim=xlims, sym=[1], markevry=[1])
    tracer([Rij.compo(0,0), Rij_tryg00 , R2ij.compo(0,0), R2ij_tryg00], 'uu_compa',textitle=r"Axial $\overline{u'u'}$", couleur=[2,2,1], xlim=xlims, sym=[1], markevry=[1])
    tracer([R4ij.compo(0,0), R2ij.compo(0,0), R2ij_tryg00, R3ij.compo(0,0)], 'uu_bidi',textitle=r"Axial $\overline{u'u'}$", couleur=[1,2,1], xlim=xlims, sym=[1], markevry=[1])
    tracer([Rij.compo(0,0), Rij_tryg00,R4ij.compo(0,0) , R2ij.compo(0,0), R2ij_tryg00, R3ij.compo(0,0)], 'uu_all',textitle=r"Axial $\overline{u'u'}$", couleur=[3,2,1], xlim=xlims, sym=[1], markevry=[1,1,10,1,1,1,1,1,1,1,1])


    Rij_tryg22=Field.LoadFromFile('VV_sphere_TRYG.dt_ev',["coordonnee_K"] , ["I"],'I', r'$z$', r'I')
    R2ij_tryg22=Field.LoadFromFile('VV_defo_TRYG.dt_ev',["coordonnee_K"] , ["I"],'I', r'$z$', r'I')
    Rij_tryg22=Rij_tryg22.power(2)
    R2ij_tryg22=R2ij_tryg22.power(2)
    Rij_tryg22.settex(r'Sb : $\mathbf{d_b=0.3}$ (LT2008)')
    R2ij_tryg22.settex(r'D : (LT2008)')
    tracer([Rij.compo(2,2), Rij_tryg22, R2ij.compo(2,2), R2ij_tryg22, R3ij.compo(2,2)], 'vv',textitle=r"Wall-normal $\overline{v'v'}$", couleur=[2,2,1], xlim=xlims, markevry=[1], sym=[1])
    tracer([Rij.compo(2,2), Rij_tryg22, R2ij.compo(2,2), R2ij_tryg22], 'vv_compa',textitle=r"Wall-normal $\overline{v'v'}$", couleur=[2,2,1], xlim=xlims, markevry=[1], sym=[1])
    tracer([R4ij.compo(2,2), R2ij.compo(2,2), R2ij_tryg22, R3ij.compo(2,2)], 'vv_bidi',textitle=r"Wall-normal $\overline{v'v'}$", couleur=[1,2,1], xlim=xlims, markevry=[1], sym=[1])
    tracer([Rij.compo(2,2), Rij_tryg22, R4ij.compo(2,2), R2ij.compo(2,2), R2ij_tryg22, R3ij.compo(2,2)], 'vv_all',textitle=r"Wall-normal $\overline{v'v'}$", couleur=[3,2,1], xlim=xlims, markevry=[1,1,10,1,1,1,1,1,1,1,1], sym=[1])

    Rij_tryg20=Field.LoadFromFile('UV_sphere_TRYG.dt_ev',["coordonnee_K"] , ["I"],'I', r'$z$', r'I')
    R2ij_tryg20=Field.LoadFromFile('UV_defo_TRYG.dt_ev',["coordonnee_K"] , ["I"],'I', 'z', r'I')
    Rij_tryg20.settex(r'spherical (LT2008)')
    R2ij_tryg20.settex(r'deformable (LT2008)')
    Rij20=Rij.compo(2,0)*(-1.)
    Rij220=R2ij.compo(2,0)*(-1.)
    Rij320=R3ij.compo(2,0)*(-1.)
    Rij420=R4ij.compo(2,0)*(-1.)
    Rij20.settex(r'spherical $\mathbf{d_b=0.3}$ (TrioCFD)')
    Rij220.settex(r'deformable (TrioCFD)')
    Rij320.settex(r'bidisperse (TrioCFD)')
    Rij420.settex(r'spherical $\mathbf{d_b=0.094}$ (TrioCFD)')
    tracer([Rij20, Rij_tryg20 , Rij220, R2ij_tryg20, Rij320], 'uv',textitle=r"Cross-correlation $\overline{u'v'}$", couleur=[2,2,1], xlim=xlims, markevry=[1])
    tracer([Rij20, Rij_tryg20 , Rij220, R2ij_tryg20], 'uv_compa',textitle=r"Cross-correlation $\overline{u'v'}$", couleur=[2,2,1], xlim=xlims, markevry=[1])
    tracer([Rij420 , Rij220, R2ij_tryg20, Rij320], 'uv_bidi',textitle=r"Cross-correlation $\overline{u'v'}$", couleur=[1,2,1], xlim=xlims, markevry=[1])
    tracer([Rij20, Rij_tryg20, Rij420 , Rij220, R2ij_tryg20, Rij320], 'uv_all',textitle=r"Cross-correlation $\overline{u'v'}$", couleur=[3,2,1], xlim=xlims,markevry=[1,1,10,1,1,1,1,1,1,1,1])




    return

############################################################################################
############################## CHARGEMENT DES CHAMPS #######################################
############################################################################################



def RunDiphasicCase(gravity, rho_l, rho_v, alpha_v, mu, Source, sigma, fstat):
    """ 
    Utilise les profils statistiques d'un cas diphasique pour calculer chaque terme des bilans de qdm (one fluid, two fluid formulation), ainsi que chaque terme de l'équation de transport des tensions de Reynolds. Retourne un dictionnaire dans lequel les résultats sont compilés pour pouvoir être utilisés dans une autre routine. 
    @param constante de gravitation
    @param masse volumique du liquide
    @param masse volumique de la vapeur
    @param taux de vide moyen du cas étudié
    @param viscosité
    @param terme source de qdm (gradient de pression moyen)
    @param tension de surface
    @param nom du fichier statistique a utilisé
    """ 

    dvar = Field.getEntries(fstat)
    I=Field.LoadFromFile(fstat,["coordonnee_K"] , ["I"],'I', 'z', r'I')
    pression_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["PI"],'P_l', 'z', r'$P_l$')
    pression_v=Field.LoadFromFile(fstat,["coordonnee_K"] , ["PIv"],'P_v', 'z', r'$P_v$')
    uu_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UUI", "UVI", "UWI","UVI", "VVI", "VWI","UWI", "VWI", "WWI"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    uu_v=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UUIv", "UVIv", "UWIv","UVIv", "VVIv", "VWIv","UWIv", "VWIv", "WWIv"],'UU_v', 'z', r'\langle(u\times u)\rangle)_v')
    uuu_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UUUI", "UUVI", "UUWI","UUVI", "UVVI", "UVWI","UUWI", "UVWI", "UWWI","UUVI", "UVVI", "UVWI","UVVI", "VVVI", "VVWI","UVWI", "VVWI", "VWWI","UUWI", "UVWI", "UWWI","UVWI", "VVWI", "VWWI","UWWI", "VWWI", "WWWI" ],'UUU', 'z', r'\langle(u\times u\times u)\rangle)_v')
    u_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UI", "VI", "WI"],r'$U_l$', r'$z$', r'$u_l$')
    u_v=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UIv", "VIv", "WIv"],'U_v', 'z', r'$u_v$')
    aii=Field.LoadFromFile(fstat,["coordonnee_K"] , ["kaiNx", "kaiNy", "kaiNz"],'f_interf', 'z', r'$\sigma n\kappa\delta^i$')
    pu_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UPI", "VPI", "WPI"],r'$Up_l$', r'$z$', r'$up_l$')
    dudxdudx_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["IdUdxdUdx", "IdUdxdVdx", "IdUdxdWdx","IdUdxdVdx", "IdVdxdVdx", "IdVdxdWdx","IdUdxdWdx", "IdVdxdWdx", "IdWdxdWdx"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    dudydudy_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["IdUdydUdy", "IdUdydVdy", "IdUdydWdy","IdUdydVdy", "IdVdydVdy", "IdVdydWdy","IdUdydWdy", "IdVdydWdy", "IdWdydWdy"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    dudzdudz_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["IdUdzdUdz", "IdUdzdVdz", "IdUdzdWdz","IdUdzdVdz", "IdVdzdVdz", "IdVdzdWdz","IdUdzdWdz", "IdVdzdWdz", "IdWdzdWdz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    PaiN=Field.LoadFromFile(fstat,["coordonnee_K"] , ["PaiNx", "PaiNy", "PaiNz"],'f_interf', 'z', r'$\sigma n\kappa\delta^i$')
    aiN=Field.LoadFromFile(fstat,["coordonnee_K"] , ["aiNx", "aiNy", "aiNz"],'f_interf', 'z', r'$\sigma n\kappa\delta^i$')
    vdvdxaiN=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UdUdxaiNx", "UdVdxaiNx", "UdWdxaiNx","UdVdxaiNx", "VdVdxaiNx", "VdWdxaiNx","UdWdxaiNx", "VdWdxaiNx", "WdWdxaiNx"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    vdvdyaiN=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UdUdyaiNy", "UdVdyaiNy", "UdWdyaiNy","UdVdyaiNy", "VdVdyaiNy", "VdWdyaiNy","UdWdyaiNy", "VdWdyaiNy", "WdWdyaiNy"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    vdvdzaiN=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UdUdzaiNz", "UdVdzaiNz", "UdWdzaiNz","UdVdzaiNz", "VdVdzaiNz", "VdWdzaiNz","UdWdzaiNz", "VdWdzaiNz", "WdWdzaiNz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    dvdxaiN=Field.LoadFromFile(fstat,["coordonnee_K"] , ["dUdxaiNx", "dVdxaiNx", "dWdxaiNx"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    dvdyaiN=Field.LoadFromFile(fstat,["coordonnee_K"] , ["dUdyaiNy", "dVdyaiNy", "dWdyaiNy"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    dvdzaiN=Field.LoadFromFile(fstat,["coordonnee_K"] , ["dUdzaiNz", "dVdzaiNz", "dWdzaiNz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    pdu_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["IPdUdx","IPdUdy","IPdUdz","IPdVdx", "IPdVdy","IPdVdz","IPdWdx","IPdWdy","IPdWdz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    vaiN=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UaiNx","UaiNy","UaiNz","VaiNx","VaiNy","VaiNz","WaiNx","WaiNy","WaiNz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    vPaiN=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UPaiNx", "UPaiNy", "UPaiNz", "VPaiNx","VPaiNy", "VPaiNz", "WPaiNx", "WPaiNy", "WPaiNz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    vvaiN=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UUaiNx","UUaiNy","UUaiNz","UVaiNx","UVaiNy","UVaiNz","UWaiNx","UWaiNy","UWaiNz", "UVaiNx","UVaiNy","UVaiNz","VVaiNx","VVaiNy","VVaiNz","VWaiNx","VWaiNy","VWaiNz","UWaiNx","UWaiNy","UWaiNz","VWaiNx","VWaiNy","VWaiNz","WWaiNx","WWaiNy","WWaiNz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    variable={'I':I,'aiN':aiN,'aii':aii, 'PaiN':PaiN, 'g':gravity, 'rho_l':rho_l, 'rho_v':rho_v, 'alpha':alpha_v,'mu': mu,'S':Source,'sigma':sigma, 'p_l':pression_l, 'p_v':pression_v, 'u_l':u_l, 'u_v':u_v}
    ############################################################################################
    ############################## PARAMETRES  #################################################
    ############################################################################################

    
    #tracer([((aii.integ()*sigma)*(1/mu)).sqrt()], 'retau_erreur_sigma_conservatisme')
    tracer([I], 'I')
    g=Field.initgravity([gravity,0,0], 'g', I)
    Iv=(I-1.0)*(-1.0)
    alpha_l=1-alpha_v
    rho_ll=rho_l*alpha_l
    rho_vv=rho_v*alpha_v
    rho_av=rho_ll+rho_vv
    S=Field.initsource([Source,0,0], 'beta', I)


    ############################################################################################
    ############################## BILAN DE QDM ################################################
    ############################################################################################
    
    rho1=I*rho_l
    rho2=Iv*rho_v
    rho=rho1+rho2
    u=u_v+u_l
    tracer([u],'u0', [0])
    tracer([u],'u1', [1])
    tracer([u],'u2', [2])
	
    ull=u_l*(I.inv(0.00001))
    uvv=u_v*(Iv.inv(0.00001))

    tracer([ull,uvv],'u_phase0', [0])
    tracer([ull,uvv],'u_phase1', [1])
    tracer([ull,uvv],'u_phase2', [2])
       
       
    uu=uu_v+uu_l
    pression=pression_l+pression_v         
    rhog_av=g*rho_av
    rhog_tmp=rho*g
    rhog=rhog_tmp-rhog_av
    ufu=u.MatProduct(u)
    Rij=uu-ufu
    ufu_l=u_l.MatProduct(u_l)
    Rij_l=uu_l-ufu_l
    ufu_v=u_v.MatProduct(u_v)
    Rij_v=uu_v-ufu_v
    variable['Rij']=Rij
    variable['Rij_l']=Rij_l
    variable['Rij_v']=Rij_v
    tracer([pression_l], 'pression')

    qdm_diph=qdm_diphase('diph',u, mu, pression, Rij, S, rho, aii, rho_av, g, sigma, I)
    qdm_liq=qdm_diphase('liquide',u_l, mu, pression_l, Rij_l, S, rho_l, aii, rho_av, g, sigma, I)
    qdm_gaz=qdm_diphase('gaz',u_v, mu, pression_v, Rij_v, S, rho_v, aii, rho_av, g, sigma, Iv)


    ############################################################################################
    ############################## BILAN DE RIJ ################################################
    ############################################################################################
    
#     uu=u_l.tensorProduct(u_l)
#     gradu=u_l.grad()
#     Rij_l=uu_l-uu
#     Rijk_l=uuu_l-uu_l.tensorProduct(u_l)+u_l.tensorProduct(u_l).tensorProduct(u_l)*2.0
#     pu_ll=pu_l-pression_l*u_l
#     pdu_ll=pdu_l-pression_l*u_l.grad()
#     dudu0=dudxdudx_l+dudydudy_l+dudzdudz_l
#     duduijkl=gradu.tensorProduct(gradu)
#     dudu_l=dudu0-duduijkl.contract('ij', 'ikjk')
    vpainn=vPaiN-pression_l*vaiN-u_l.tensorProduct(PaiN)+pression_l*u_l.tensorProduct(aiN)
    uduain=vdvdxaiN+vdvdyaiN+vdvdzaiN
    duain=dvdxaiN+dvdyaiN+dvdzaiN
    
    ## calcul des fluctuations
    Rij_l=fluctuij(u_l, u_l, uu_l)
    pdu_ll=fluctuij(pression_l, u_l.grad(), pdu_l)
    pu_ll=fluctuij(pression_l, u_l, pu_l)
    dudu0=dudxdudx_l+dudydudy_l+dudzdudz_l
    dudu_l=fluctuij(u_l.grad(), u_l.grad(), dudu0)
    Rijk_l=fluctuijk(uuu_l, u_l, uu_l)         
    
    variable2={'Rij':Rij_l, 'Rijk':Rijk_l}
    Rij_diph=Rij_diphase(u_l, Rij_l, Rijk_l, pu_ll, pdu_ll, rho_l, dudu_l, mu, vpainn, uduain, duain, vaiN, aiN, vvaiN)
    res={'qdm_diph':qdm_diph,'qdm_liq':qdm_liq,'qdm_gaz':qdm_gaz, 'Rij':Rij_diph, 'variable':variable, 'variable2':variable2}
    return res
############################################################################################
############################## COMPARAISON MONOPHASIQUE ####################################
############################################################################################

def RunVremanCase(Run, fstat1, fstat2, rho, mu, Source):
    """ 
    Utilise les profils statistiques d'un cas de la litterature (vreman) pour calculer chaque terme des bilans de qdm (one fluid, two fluid formulation), ainsi que chaque terme de l'équation de transport des tensions de Reynolds. Retourne un dictionnaire dans lequel les résultats sont compilés pour pouvoir être utilisés dans une autre routine. 
    """
    dvar1 = Field.getEntries(fstat1)
    dvar2 = Field.getEntries(fstat2)
    u=Field.LoadFromFile(fstat1,["y"] , ["u1", "u3", "u2"],r'$U_l$', r'$z$', r'$u_l$')
    Rij=Field.LoadFromFile(fstat2,["y"] , ['u1u1', 'u1u3', 'u1u2', 'u1u3', 'u3u3', 'u2u3', 'u1u2', 'u2u3', 'u2u2', ],'P_l', 'z', r'$P_l$')
    Rijk=Field.LoadFromFile(fstat2,["y"] , ['0','0','u1u1u2','0','0','u1u3u2','u1u1u2','u1u3u2','u1u2u2','0','0','u1u3u2','0','0','u3u3u2','u1u3u2','u3u3u2','u3u2u2','u1u1u2','u1u3u2','u1u2u2','u1u3u2','u3u3u2','u3u2u2','u1u2u2','u3u2u2','u2u2u2'],'P_l', 'z', r'$P_l$')
    pression=Field.LoadFromFile(fstat1,["y"] , ["p"],'P_l', 'z', r'$P_l$')
    pu=Field.LoadFromFile(fstat2,["y"] , ["u1p", "u3p", "u2p"],r'$Up_l$', r'$z$', r'$up_l$')
    dudxdudx=Field.LoadFromFile(fstat2,["y"] , ["u11u11", "u11u31", "u11u21","u11u31", "u31u31", "u21u31","u11u21", "u21u31", "u21u21"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    dudydudy=Field.LoadFromFile(fstat2,["y"] , ["u13u13", "u13u33", "u13u23","u13u33", "u33u33", "u23u33","u13u23", "u23u33", "u23u23"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    dudzdudz=Field.LoadFromFile(fstat2,["y"] , ["u12u12", "u12u32", "u12u22","u12u32", "u32u32", "u22u32","u12u22", "u22u32", "u22u22"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    pdu=Field.LoadFromFile(fstat2,["y"] , ["u11p","u13p","u12p","u31p","u33p","u32p","u21p","u23p","u22p"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')

    S=Field.initsource([Source,0,0], 'beta', u)
    dudu0=dudxdudx+dudydudy+dudzdudz
    u._ax[:,:]=u._ax[:,:]*mu
    Rij._ax[:,:]=Rij._ax[:,:]*mu
    Rijk._ax[:,:]=Rijk._ax[:,:]*mu
    pression._ax[:,:]=pression._ax[:,:]*mu
    pu._ax[:,:]=pu._ax[:,:]*mu
    dudxdudx._ax[:,:]=dudxdudx._ax[:,:]*mu
    dudydudy._ax[:,:]=dudydudy._ax[:,:]*mu
    dudzdudz._ax[:,:]=dudzdudz._ax[:,:]*mu
    pdu._ax[:,:]=pdu._ax[:,:]*mu
    qdm_vreman=qdm_monophase(Run,'mono', u, mu, pression, Rij, S, rho)
    Rij_vreman=Rij_monophase(Run, u, Rij, Rijk, pu, pdu, rho, dudu0, mu)
    res={'qdm':qdm_vreman, 'Rij':Rij_vreman, 'Urms':Rij}
    
    return res

    
def RunMonophasicCase(Run, fstat, rho, mu, Source):
    
    """ 
    Utilise les profils statistiques d'un cas monophasique pour calculer chaque terme du bilan de qdm, ainsi que chaque terme de l'équation de transport des tensions de Reynolds. Retourne un dictionnaire dans lequel les résultats sont compilés pour pouvoir être utilisés dans une autre routine. 
    """
    dvar = Field.getEntries(fstat)
    I=Field.LoadFromFile(fstat,["coordonnee_K"] , ["I"],'I', 'z', r'I')
    pression=Field.LoadFromFile(fstat,["coordonnee_K"] , ["PI"],'P_l', 'z', r'$P_l$')
    uu=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UUI", "UVI", "UWI","UVI", "VVI", "VWI","UWI", "VWI", "WWI"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    uuu=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UUUI", "UUVI", "UUWI","UUVI", "UVVI", "UVWI","UUWI", "UVWI", "UWWI","UUVI", "UVVI", "UVWI","UVVI", "VVVI", "VVWI","UVWI", "VVWI", "VWWI","UUWI", "UVWI", "UWWI","UVWI", "VVWI", "VWWI","UWWI", "VWWI", "WWWI" ],'UUU', 'z', r'\langle(u\times u\times u)\rangle)_v')
    u=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UI", "VI", "WI"],r'$U_l$', r'$z$', r'$u_l$')
    aii=Field.LoadFromFile(fstat,["coordonnee_K"] , ["kaiNx", "kaiNy", "kaiNz"],'f_interf', 'z', r'$\sigma n\kappa\delta^i$')
    pu=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UPI", "VPI", "WPI"],r'$Up_l$', r'$z$', r'$up_l$')
    dudxdudx=Field.LoadFromFile(fstat,["coordonnee_K"] , ["IdUdxdUdx", "IdUdxdVdx", "IdUdxdWdx","IdUdxdVdx", "IdVdxdVdx", "IdVdxdWdx","IdUdxdWdx", "IdVdxdWdx", "IdWdxdWdx"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    dudydudy=Field.LoadFromFile(fstat,["coordonnee_K"] , ["IdUdydUdy", "IdUdydVdy", "IdUdydWdy","IdUdydVdy", "IdVdydVdy", "IdVdydWdy","IdUdydWdy", "IdVdydWdy", "IdWdydWdy"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    dudzdudz=Field.LoadFromFile(fstat,["coordonnee_K"] , ["IdUdzdUdz", "IdUdzdVdz", "IdUdzdWdz","IdUdzdVdz", "IdVdzdVdz", "IdVdzdWdz","IdUdzdWdz", "IdVdzdWdz", "IdWdzdWdz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    pdu=Field.LoadFromFile(fstat,["coordonnee_K"] , ["IPdUdx","IPdUdy","IPdUdz","IPdVdx", "IPdVdy","IPdVdz","IPdWdx","IPdWdy","IPdWdz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    S=Field.initsource([Source,0,0], 'beta', u)
    ## calcul des fluctuations
    Rij=fluctuij(u, u, uu)
    pdu=fluctuij(pression, u.grad(), pdu)
    pu=fluctuij(pression, u, pu)
    dudu0=dudxdudx+dudydudy+dudzdudz
    dudu0=fluctuij(u.grad(), u.grad(), dudu0)
    Rijk=fluctuijk(uuu, u, uu)
    
    qdm_mono=qdm_monophase(Run, 'mono', u, mu, pression, Rij, S, rho)
    Rij_mono=Rij_monophase(Run, u, Rij, Rijk, pu, pdu, rho, dudu0, mu)
    
    res={'qdm':qdm_mono, 'Rij':Rij_mono, 'Urms':Rij}
    return res  




def RunMonophasicCase127(Run, fstat, rho, mu, Source):
    
    """ 
    Utilise les profils statistiques d'un cas monophasique pour calculer chaque terme du bilan de qdm, ainsi que chaque terme de l'équation de transport des tensions de Reynolds. Retourne un dictionnaire dans lequel les résultats sont compilés pour pouvoir être utilisés dans une autre routine. 
    """
    dvar = Field.getEntries(fstat)
    pression=Field.LoadFromFile(fstat,["coordonnee_K"] , ["P"],'P_l', 'z', r'$P_l$')
    uu=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UU", "UV", "UW","UV", "VV", "VW","UW", "VW", "WW"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    uuu=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UUU", "UUV", "UUW","UUV", "UVV", "UVW","UUW", "UVW", "UWW","UUV", "UVV", "UVW","UVV", "VVV", "VVW","UVW", "VVW", "VWW","UUW", "UVW", "UWW","UVW", "VVW", "VWW","UWW", "VWW", "WWW" ],'UUU', 'z', r'\langle(u\times u\times u)\rangle)_v')
    u=Field.LoadFromFile(fstat,["coordonnee_K"] , ["U", "V", "W"],r'$U_l$', r'$z$', r'$u_l$')
    pu=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UP", "VP", "WP"],r'$Up_l$', r'$z$', r'$up_l$')
    dudxdudx=Field.LoadFromFile(fstat,["coordonnee_K"] , ["dUdxdUdx", "dUdxdVdx", "dUdxdWdx","dUdxdVdx", "dVdxdVdx", "dVdxdWdx","dUdxdWdx", "dVdxdWdx", "dWdxdWdx"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    dudydudy=Field.LoadFromFile(fstat,["coordonnee_K"] , ["dUdydUdy", "dUdydVdy", "dUdydWdy","dUdydVdy", "dVdydVdy", "dVdydWdy","dUdydWdy", "dVdydWdy", "dWdydWdy"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    dudzdudz=Field.LoadFromFile(fstat,["coordonnee_K"] , ["dUdzdUdz", "dUdzdVdz", "dUdzdWdz","dUdzdVdz", "dVdzdVdz", "dVdzdWdz","dUdzdWdz", "dVdzdWdz", "dWdzdWdz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    pdu=Field.LoadFromFile(fstat,["coordonnee_K"] , ["PdUdx","PdUdy","PdUdz","PdVdx", "PdVdy","PdVdz","PdWdx","PdWdy","PdWdz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    S=Field.initsource([Source,0,0], 'beta', u)
    ## calcul des fluctuations
    Rij=fluctuij(u, u, uu)
    pdu=fluctuij(pression, u.grad(), pdu)
    
    debit=(u*rho).mean_compo(0)
    reb=debit*(4/mu)
    vitesse=u.mean_compo(0)
    rebb=vitesse*(4/mu)
    print 'Reb=', reb
    print 'Rebb=', rebb, 'pour un ul=', vitesse
    
    pu=fluctuij(pression, u, pu)
    dudu0=dudxdudx+dudydudy+dudzdudz
    dudu0=fluctuij(u.grad(), u.grad(), dudu0)
    Rijk=fluctuijk(uuu, u, uu)
    
    qdm_monophase(Run, 'mono', u, mu, pression, Rij, S, rho)
    Rij_monophase(Run, u, Rij, Rijk, pu, pdu, rho, dudu0, mu)
    
    return  



#########################################################################################
###################################### MAIN #############################################
#########################################################################################



def Mk(DeformableRun):
    """ 
    Trace les différents termes de transfert de qdm entre les phases à partir du dictionnaire de sauvegarde d'un cas. Necessite d'avoir fait tourné un Run[...]case avant pour obtenir le      dictionnaire. 
    """
    
    #### chargement des donnees
    I=DeformableRun['variable']['I']
    M_l=DeformableRun['qdm_liq']['interface']
    M_v=DeformableRun['qdm_gaz']['interface']
    M_sigma_code=DeformableRun['variable']['aii']*DeformableRun['variable']['sigma']*(-1.0)
    M_sigma_code.settex('interface monofluide code')    
    p_v=DeformableRun['variable']['p_v']
    p_l=DeformableRun['variable']['p_l']
    aiN=DeformableRun['variable']['aiN']
    M_sigma_mono=DeformableRun['qdm_diph']['interface']
    
    DP=p_v*I.grad()*(-1.0)+p_l*I.grad()
    Mi_v=M_v+p_v*I.grad()*(-1.0)#### ATTENTION Ain change de signe en fonction de la phase
    Mi_l=M_l+p_l*I.grad()
    ###### Calculs
    
    M_sigma_bifluide=M_l+M_v
    Mi_sigma=Mi_l+Mi_v
    
    MlNet=NetContribution(M_l)
    MvNet=NetContribution(M_v)
    MsigNet=NetContribution(M_sigma_bifluide)
    ###### PRINT
    MlNet.settex(r'$\langle M_{l}\rangle$')
    MvNet.settex(r'$\langle M_{v}\rangle$')
    MsigNet.settex(r'$\langle M_{\sigma}\rangle$')
    M_l.settex(r'$M_l=M_{il}+M_{\sigma l}$')
    M_v.settex(r'$M_v=M_{iv}+M_{\sigma v}$')
    M_sigma_bifluide.settex(r'$M_{\sigma}=M_l+M_v$')
    M_sigma_code.settex(r'$M_{\sigma}=\overline{\sigma\kappa\mathbf{n}}^i$')
    M_l.labelRef='z'
    
    tracer([M_l, MlNet, M_v, MvNet, M_sigma_bifluide, M_sigma_code, MsigNet], 'Mx', [0], couleur=[2,2,3], markevry=[4])
    tracer([M_l, MlNet, M_v, MvNet, M_sigma_bifluide, M_sigma_code, MsigNet], 'My', [1], couleur=[2,2,3], markevry=[4])  
    tracer([M_l, MlNet, M_v, MvNet, M_sigma_bifluide, M_sigma_code, MsigNet], 'Mz', [2], couleur=[2,2,3], markevry=[4])   
      
    return

def tracerQdmDiff(Mono):
    msig=Mono['variable']['aii']*Mono['variable']['sigma']*(-1.0)
    pl=Mono['variable']['p_l']
    pv=Mono['variable']['p_v']
    pl.settex(r'$p_l$')
    pv.settex(r'$p_v$')
    msig.settex(r'interface $\sigma\kappa\mathbf{n}\delta^i$')
    
    tracer([pl, pv], 'p')
    tracer([msig, Mono['qdm_diph']['interface']], 'valid_qdm_x', [0])
    tracer([msig, Mono['qdm_diph']['interface']], 'valid_qdm_y', [1])
    tracer([msig, Mono['qdm_diph']['interface']], 'valid_qdm_z', [2])        
    
    tracer([msig, Mono['qdm_diph']['interface'],Mono['qdm_diph']['viscous'],Mono['qdm_diph']['turbulence'],Mono['qdm_diph']['inertie'],Mono['qdm_diph']['pression'],Mono['qdm_diph']['source'],Mono['qdm_diph']['rhog']], 'qdm_i', [0], couleur=[2,1,1,1,1,1,1,1,1,1,1,1,1])
    tracer([msig, Mono['qdm_diph']['interface'],Mono['qdm_diph']['viscous'],Mono['qdm_diph']['turbulence'],Mono['qdm_diph']['inertie'],Mono['qdm_diph']['pression'],Mono['qdm_diph']['source'],Mono['qdm_diph']['rhog']], 'qdm_j', [1], couleur=[2,1,1,1,1,1,1,1,1,1,1,1,1])
    tracer([msig, Mono['qdm_diph']['interface'],Mono['qdm_diph']['viscous'],Mono['qdm_diph']['turbulence'],Mono['qdm_diph']['inertie'],Mono['qdm_diph']['pression'],Mono['qdm_diph']['source'],Mono['qdm_diph']['rhog']], 'qdm_k', [2], couleur=[2,1,1,1,1,1,1,1,1,1,1,1,1])
    tracer([msig, Mono['qdm_gaz']['interface'],Mono['qdm_gaz']['viscous'],Mono['qdm_gaz']['turbulence'],Mono['qdm_gaz']['inertie'],Mono['qdm_gaz']['pression'],Mono['qdm_gaz']['source'],Mono['qdm_gaz']['rhog']], 'qdmgaz_i', [0], couleur=[2,1,1,1,1,1,1,1,1,1,1,1,1])
    tracer([msig, Mono['qdm_gaz']['interface'],Mono['qdm_gaz']['viscous'],Mono['qdm_gaz']['turbulence'],Mono['qdm_gaz']['inertie'],Mono['qdm_gaz']['pression'],Mono['qdm_gaz']['source'],Mono['qdm_gaz']['rhog']], 'qdmgaz_j', [1], couleur=[2,1,1,1,1,1,1,1,1,1,1,1,1])
    tracer([msig, Mono['qdm_gaz']['interface'],Mono['qdm_gaz']['viscous'],Mono['qdm_gaz']['turbulence'],Mono['qdm_gaz']['inertie'],Mono['qdm_gaz']['pression'],Mono['qdm_gaz']['source'],Mono['qdm_gaz']['rhog']], 'qdmgaz_k', [2], couleur=[2,1,1,1,1,1,1,1,1,1,1,1,1])
    tracer([msig, Mono['qdm_liq']['interface'],Mono['qdm_liq']['viscous'],Mono['qdm_liq']['turbulence'],Mono['qdm_liq']['inertie'],Mono['qdm_liq']['pression'],Mono['qdm_liq']['source'],Mono['qdm_liq']['rhog']], 'qdmliq_i', [0], couleur=[2,1,1,1,1,1,1,1,1,1,1,1,1])
    tracer([msig, Mono['qdm_liq']['interface'],Mono['qdm_liq']['viscous'],Mono['qdm_liq']['turbulence'],Mono['qdm_liq']['inertie'],Mono['qdm_liq']['pression'],Mono['qdm_liq']['source'],Mono['qdm_liq']['rhog']], 'qdmliq_j', [1], couleur=[2,1,1,1,1,1,1,1,1,1,1,1,1])
    tracer([msig, Mono['qdm_liq']['interface'],Mono['qdm_liq']['viscous'],Mono['qdm_liq']['turbulence'],Mono['qdm_liq']['inertie'],Mono['qdm_liq']['pression'],Mono['qdm_liq']['source'],Mono['qdm_liq']['rhog']], 'qdmliq_k', [2], couleur=[2,1,1,1,1,1,1,1,1,1,1,1,1])
        
    
    return

   
 

def Rij_monophase(Run, u, Rij, Rijk, pu, pdu, rho, dudu0, mu, vpain=None, uduain=None, duain=None, vaiN=None, aiN=None, vvaiN=None):
    """ 
    Calcule, trace et sauvegarde les différents termes de l'équation de transport des tensions de Reynolds pour un cas monophasique
    """
    
    a=u.tensorProduct(Rij)
    b=a.grad()
    nonlin1=b.contract('ij','kijk')
    gradRij=Rij.grad()
    gradu=u.grad()
    nonlin2=Rij.tensorProduct(gradu)
    nonlin2=nonlin2.contract('ij','kijk')
    nonlin3=nonlin2.transpo()      
    nonlin=nonlin1+nonlin2+nonlin3
    gradRijk=Rijk.grad()
    nonlin4=gradRijk.contract('ij','kijk')

    ####### terme de pression ###############
    
    pijk=pu.grad()
    p1ijk=pijk.transpo()+pijk

    p2ijk=pdu.transpo()+pdu
    p=p2ijk-p1ijk
    
    ####### viscosite #########
    
    
    ddRij0=gradRij.grad()
    ddRij=ddRij0.contract('ij', 'ijkk')
    
    
    #################### rho, mu et nu #############################
    nu=mu/rho
    advection=nonlin1
    production=nonlin2+nonlin3
    redistribution=p2ijk*(1/rho)
    dissipation=dudu0*nu*2.0
    diffusion=nonlin4
    transpression=p1ijk*(1/rho)
    diffusion_mol=ddRij*nu
    
    #################### label latex #############################
    
    advection.settex(r'advection $\frac{DR_{i,j}}{Dt}$')
    production.settex(r'production $\overline{v^,_lv^,_j}\frac{\partial\overline{v_i}}{\partial x_l}+\overline{v^,_lv^,_i}\frac{\partial\overline{v_j}}{\partial x_l}$')
    redistribution.settex(r'redistribution $\overline{p^,\left(\frac{\partial v^,_j}{\partial x_i}+\frac{\partial v^,_i}{\partial x_j}\right)}$')
    dissipation.settex(r'dissipation $2\nu\overline{\frac{\partial v^,_i}{\partial x_l} \frac{\partial v^,_j}{\partial x_l}}$')
    diffusion.settex(r'turbulent diffusion $\frac{\partial\overline{v^,_lv^,_iv^,_j}}{\partial x_l}$')
    transpression.settex(r'pressure diffusion $\frac{\partial}{\partial x_l}(\overline{p^,v^,_j}\delta_{il}+\overline{p^,v^,_i}\delta_{jl})$')
    diffusion_mol.settex(r'molecular diffusion  $\nu\frac{\partial^2\overline{v^,_iv^,_j}}{\partial^2 x_l}$')
    
    residu=diffusion_mol-diffusion-advection-production-dissipation+redistribution-transpression
    #residu=visc-nonlin-nonlin4+p
    residu.settex(r'residu')
    tracer([production, redistribution, dissipation, diffusion, transpression,  diffusion_mol,residu], "Rij_mono_11", [0])
    tracer([production, redistribution, dissipation, diffusion, transpression,  diffusion_mol,residu], "Rij_mono_12", [1])
    tracer([production, redistribution, dissipation, diffusion, transpression,  diffusion_mol,residu], "Rij_mono_13", [2])
    tracer([production, redistribution, dissipation, diffusion, transpression,  diffusion_mol,residu], "Rij_mono_21", [3])
    tracer([production, redistribution, dissipation, diffusion, transpression,  diffusion_mol,residu], "Rij_mono_22", [4])
    tracer([production, redistribution, dissipation, diffusion, transpression,  diffusion_mol,residu], "Rij_mono_23", [5])
    tracer([production, redistribution, dissipation, diffusion, transpression,  diffusion_mol,residu], "Rij_mono_31", [6])
    tracer([production, redistribution, dissipation, diffusion, transpression,  diffusion_mol,residu], "Rij_mono_32", [7])
    tracer([production, redistribution, dissipation, diffusion, transpression,  diffusion_mol,residu], "Rij_mono_33", [8])
    
    Add_section(Run, 'Rij')
    Run['Rij']['production']=production
    Run['Rij']['redistribution']=redistribution
    Run['Rij']['dissipation']=dissipation
    Run['Rij']['diffusion']=diffusion
    Run['Rij']['transpression']=transpression
    Run['Rij']['diffusion_mol']=diffusion_mol
    Run['Rij']['residu']=residu

    ###################### plot ##################################
    if vpain==None and uduain==None and duain==None and vaiN==None and aiN==None and vvaiN==None and 0: 
        
        
            ###### Elliptic blending reynolds stress model ######
        #### donnees d'entree 
        g1=3.4
        g1_star=1.8
        g3=0.8
        g3_star=1.3
        g4=1.25
        g5=0.4
        C_eta=80.0
        C_l=0.133
        Prod=production.contract(0, 'kk')
        epsilon=dissipation.contract(0, 'kk') ### divise par 2 ou pas ?  
        k_cine=Rij.contract(0, 'kk')*0.5 #### negatif ou positif ?
        eps=epsilon
        eps_inv=eps.inv()
        k=k_cine
        k_inv=k.inv()
        P=Prod
        
        #eps=epsilon.mean()  ### *0.5 ?
        #k=k_cine.mean()
        #P=Prod.mean()
        ### matrice identite
        Ident=Rij.transpo()
        Ident._npa[:,:,:]=0.0
        Ident._npa[0,0,:]=1.0
        Ident._npa[1,1,:]=1.0
        Ident._npa[2,2,:]=1.0 
        tracer([Ident], 'k')
            
        #### calcul de phij_h (modele SSG)
        bij=Rij*(k_inv*0.5)-Ident*(0.33333333)
    
        
        gradu=u.grad()
        Sij=gradu+gradu.transpo()
        Sij=Sij*0.5
        Wij=gradu-gradu.transpo()
        Wij=Wij*0.5
        
        phij_h_11=bij*g1
        phij_h_12=bij*P*eps_inv*g1_star
        phij_h_1=phij_h_11+phij_h_12
        phij_h_1=phij_h_1*eps
        
        tracer([bij], 'bij')
        tracer([Rij], 'Rij')
        tracer([Sij], 'Sij')
        tracer([Wij], 'Wij')
            
        etape1=bij.tensorProduct(bij)
        etape2=etape1.contract(0,'klkl')
        etape4=etape2.sqrt()
        phij_h_21=Sij*k*g3
        phij_h_22=Sij*k*etape4*g3_star
        phij_h_2=phij_h_21-phij_h_22
        
    
        a=bij.tensorProduct(Sij)
        a=a.contract('ij','ikjk')
        b=bij.tensorProduct(Sij)
        b=b.contract(0,'klkl')
        b=b.mean()
        phij_h_3=a+a.transpo()-Ident*b*0.66666666
        phij_h_3=phij_h_3*k*g4
        
        c=bij.tensorProduct(Wij)
        c=c.contract('ij','ikjk')
        phij_h_4=c+c.transpo()
        phij_h_4=phij_h_4*k*g5
        
        phij_h=phij_h_2+phij_h_3+phij_h_4-phij_h_1
        phij_h.settex('SSG model')
        
    ##### calcul de phij_w (modele aux parois)
    ########## longueur caracteristique correlation vitesse
        L_1=k.power(1.5)*eps_inv*C_l
        eps2=eps.power(0.25)
        eps2=eps2.inv()
        L_2=eps2*mu**0.75*C_eta*C_l
        L=L_1+L_2
        for  i in range(len(L_1._npa[0,:])):
            L._npa[0,i]=max(L_1._npa[0,i],L_2._npa[0,i])
#        L=L.symetriser()
        L2=L*L
        tracer([L2], 'L')
        
    #####################################################
    ###### resolution de alpha par une methode explicite
    ####################################################
#         alpha=L_1+L_2
#         alpha._npa[0,:]=0
#         dx=alpha._ax[0,10]-alpha._ax[0,9]
#         for i in range(len(L_1._npa[0,:])-1):
#             #L22=L2._npa[0,i]
#             L22=L2.mean()           
#             print dx, L22, dx*dx/L22, alpha._npa[0,i+1], alpha._npa[0,i], alpha._npa[0,i-1]
#             print ''
#             if i==0:
#                 alpha._npa[0,i+1]=-dx*dx/L22          
#             else :
#                 alpha._npa[0,i+1]=-dx*dx/L22+alpha._npa[0,i]*(dx*dx/L22+2)-alpha._npa[0,i-1]
#         tracer([alpha], 'alphaaa')
    ####################################################
    ###### resolution alpha par une methode implicite
    ####################################################
#         dx=L2._ax[0,10]-L2._ax[0,9]
#         n=len(L_1._npa[0,:])
#         n=n/2
#         alpha=np.zeros(n*n).reshape(n,n)
#         for i in range(n):
#             for j in range(n):
#                 if i==j :
#                     alpha[i,j]=-2.0-dx*dx/L2._npa[0,i]                      
#                 if i==j+1:
#                     alpha[i,j]=1.0
#                 if i==j-1:
#                     alpha[i,j]=1.0  
#                      
#         ### CL au bord et au centre
#         alpha[0,0]=1
#         alpha[0,1]=0
#         alpha[n-1,n-1]=1
#         alpha[n-1,n-2]=0
#         
#         ### inversion matrice
#         inver=np.linalg.inv(alpha)                
#          
#         #### remplissage vecteur b
#         b=np.zeros(n)  
#         for i in range(n):
#             b[i]=dx*dx/L2._npa[0,i] 
#         ### CL vecteur b
#         b[0]=0
#         b[n-1]=1
#         print alpha
#         print ''
#         print ''
#         print inver[:,:]
#         print ''
#         print ''    
#         print b[:]
#          
#         ### calcul de la solution 
#         solu=np.zeros(n)
#         for i in range(n):
#             for j in range(n):
#                 solu[i]=solu[i]+b[i]*inver[i,j]
#         print solu
#          
# 
#         ## met le resultat dans un objet a tracer
#         alpha=L2+L2
#         alpha._npa[0,0:n]=solu[:]
#         tracer([alpha], 'alpha')   
#     
#         tata    
    ########################################################
    ########################################################
    
    
    ########### solution analytique pour alpha avec coeff constant
    
        L=0.025
        #L=L2.mean()
         
        B=math.exp(1/L)/(math.exp(-1/L)-math.exp(1/L))
        A=-(1+B)
        alpha=ScalarField('alpha', u._ax, u.labelRef)
        aa = len(u._npa[0,:])    
        alpha._npa=np.zeros((1,aa))
        alpha._npa[0,:]=0
        for i in range(0, len(alpha._npa[0,:])/2+1 ):
            alpha._npa[0,i]=A*math.exp(alpha._ax[0,i]/L)+B*math.exp(-alpha._ax[0,i]/L)+1
        ##### ON A TRICHE ICI. SYMETRISER LA SOLUTION POUR SATISFAIRE LES CONDITIONS
            alpha._npa[0,-i-1]=A*math.exp(alpha._ax[0,i]/L)+B*math.exp(-alpha._ax[0,i]/L)+1 
             
    ###########################################################################################
    ###########################################################################################
        n=alpha.grad()
        norme=max(n._npa[2,:])
        n=n*(1./norme)
        alpha3=alpha
        alpha3._npa[0,:]=alpha._npa[0,:]*alpha._npa[0,:]*alpha._npa[0,:]    
        tracer([n, alpha3], "alpha")  
         
        phij_w_1=Rij.tensorProduct(n).tensorProduct(n)
        phij_w_1=phij_w_1.contract('ij','ikjk')
        phij_w_1=phij_w_1+phij_w_1.transpo()
        c=Rij.tensorProduct(n).tensorProduct(n)
        c=c.contract(0,'klkl')
        
        nn=n.tensorProduct(n)
        phij_w_2=c*nn+c*Ident
        phij_w_2=phij_w_2*(-0.5)
        phij_w=phij_w_1+phij_w_2
        phij_w=phij_w*eps*k_inv*(5.0)
    #    phij_w=phij_w*eps*k_inv*(-5.0)   
        #phij_w=phij_w*eps.mean()*k_inv.mean()*(-5.0)
        phij_w.settex('model paroi')
        phij_w_1.settex('1') 
        phij_w_2.settex('2')  
    #     ###### phi_tot 
         
        phij=phij_h+phij_w
        phij._npa[:,:,:]=0.0
        for i in range(len(phij._npa[0,0,:])):
            a=alpha3._npa[0,i]
            phij._npa[:,:,i]=phij_w._npa[:,:,i]*(1.0-a)+phij_h._npa[:,:,i]*a
        phij.settex('elliptic blending model')
        
        ########## 
     
        phi_Rij=Rij*eps*k_inv*5.0
        phi_Rij._npa[0,1,:]=0.0
        phi_Rij._npa[1,0,:]=0.0
    ##### nouvelle repartition de la trace propose
        phi_Rij._npa[0,0,:]=phi_Rij._npa[2,2,:]*(-2.0)
        phi_Rij._npa[1,1,:]=phi_Rij._npa[2,2,:]*(1.0)
        #### trace nulle classique
        phi_Rij_MANC=phi_Rij*(1.0)
        phi_Rij_MANC._npa[0,0,:]=phi_Rij._npa[2,2,:]*(-0.5)
        phi_Rij_MANC._npa[1,1,:]=phi_Rij._npa[2,2,:]*(-0.5) 
    ####### expression rigoureuse
        phi_Rij_RIG=phi_Rij*(1.0)
        phi_Rij_RIG._npa[0,0,:]=0.0
        phi_Rij_RIG._npa[1,1,:]=0.0  
        
        phi_Rij.settex(r'$\phi_{ij}^W$ avec $\phi_{11}^W=-2\phi_{33}^W$ et $\phi_{22}^W=\phi_{33}^W$')
        phi_Rij_RIG.settex(r'$\phi_{ij}^W$ avec $\phi_{11}^W=0$ et $\phi_{22}^W=0$')     
        phi_Rij_MANC.settex(r'$\phi_{ij}^W$ avec $\phi_{11}^W=-1/2\phi_{33}^W$ et $\phi_{22}^W=-1/2\phi_{33}^W$')
        
        #########################################
        tracer([redistribution, phij_h, phij], "model_11_mono", [0], log=True, xlim=[0.01,1])
        tracer([redistribution, phij_h, phij, phi_Rij], "model_12_mono", [1], log=True, xlim=[0.01,1])
        tracer([redistribution, phij_h, phij, phi_Rij], "model_13_mono", [2], log=True, xlim=[0.01,1])
        tracer([redistribution, phij_h, phij, phi_Rij], "model_21_mono", [3], log=True, xlim=[0.01,1])
        tracer([redistribution, phij_h, phij, phi_Rij, phi_Rij_RIG, phi_Rij_MANC], "model_22_mono", [4], log=True, xlim=[0.01,1])
        tracer([redistribution, phij_h, phij, phi_Rij], "model_23_mono", [5], log=True, xlim=[0.01,1])
        tracer([redistribution, phij_h, phij, phi_Rij], "model_31_mono", [6], log=True, xlim=[0.01,1])
        tracer([redistribution, phij_h, phij, phi_Rij], "model_32_mono", [7], log=True, xlim=[0.01,1])
        tracer([redistribution, phij_h, phij, phi_Rij], "model_33_mono", [8], log=True, xlim=[0.01,1])   
        

        
        production_R11=Field.LoadFromFile('Stats_R11.dt_ev',["y"] , ["production"],'P_l', 'z', r'$P_l$')
        pressurestrain_R11=Field.LoadFromFile('Stats_R11.dt_ev',["y"] , ["pressurestrain"],'P_l', 'z', r'$P_l$')
        dissipation_R11=Field.LoadFromFile('Stats_R11.dt_ev',["y"] , ["dissipation"],'P_l', 'z', r'$P_l$')
        turbulenttransport_R11=Field.LoadFromFile('Stats_R11.dt_ev',["y"] , ["turbulenttransport"],'P_l', 'z', r'$P_l$')
        pressuretransport_R11=Field.LoadFromFile('Stats_R11.dt_ev',["y"] , ["pressuretransport"],'P_l', 'z', r'$P_l$')
        viscoustransport_R11=Field.LoadFromFile('Stats_R11.dt_ev',["y"] , ["viscoustransport"],'P_l', 'z', r'$P_l$')
        total_R11=Field.LoadFromFile('Stats_R11.dt_ev',["y"] , ["total"],'P_l', 'z', r'$P_l$')
        production_R11.settex(r'production')
        pressurestrain_R11.settex(r'pressure strain')
        dissipation_R11.settex(r'dissipation')
        turbulenttransport_R11.settex(r'turbulent transport$')
        pressuretransport_R11.settex(r'pressure transport')
        viscoustransport_R11.settex(r'viscous transport')
        total_R11.settex(r'total')
          
        production_R12=Field.LoadFromFile('Stats_R12.dt_ev',["y"] , ["production"],'P_l', 'z', r'$P_l$')
        pressurestrain_R12=Field.LoadFromFile('Stats_R12.dt_ev',["y"] , ["pressurestrain"],'P_l', 'z', r'$P_l$')
        dissipation_R12=Field.LoadFromFile('Stats_R12.dt_ev',["y"] , ["dissipation"],'P_l', 'z', r'$P_l$')
        turbulenttransport_R12=Field.LoadFromFile('Stats_R12.dt_ev',["y"] , ["turbulenttransport"],'P_l', 'z', r'$P_l$')
        pressuretransport_R12=Field.LoadFromFile('Stats_R12.dt_ev',["y"] , ["pressuretransport"],'P_l', 'z', r'$P_l$')
        viscoustransport_R12=Field.LoadFromFile('Stats_R12.dt_ev',["y"] , ["viscoustransport"],'P_l', 'z', r'$P_l$')
        total_R12=Field.LoadFromFile('Stats_R12.dt_ev',["y"] , ["total"],'P_l', 'z', r'$P_l$')
        production_R12.settex(r'production')
        pressurestrain_R12.settex(r'pressure strain')
        dissipation_R12.settex(r'dissipation')
        turbulenttransport_R12.settex(r'turbulent transport$')
        pressuretransport_R12.settex(r'pressure transport')
        viscoustransport_R12.settex(r'viscous transport')
        total_R12.settex(r'total')       
        
        production_R22=Field.LoadFromFile('Stats_R22.dt_ev',["y"] , ["production"],'P_l', 'z', r'$P_l$')
        pressurestrain_R22=Field.LoadFromFile('Stats_R22.dt_ev',["y"] , ["pressurestrain"],'P_l', 'z', r'$P_l$')
        dissipation_R22=Field.LoadFromFile('Stats_R22.dt_ev',["y"] , ["dissipation"],'P_l', 'z', r'$P_l$')
        turbulenttransport_R22=Field.LoadFromFile('Stats_R22.dt_ev',["y"] , ["turbulenttransport"],'P_l', 'z', r'$P_l$')
        pressuretransport_R22=Field.LoadFromFile('Stats_R22.dt_ev',["y"] , ["pressuretransport"],'P_l', 'z', r'$P_l$')
        viscoustransport_R22=Field.LoadFromFile('Stats_R22.dt_ev',["y"] , ["viscoustransport"],'P_l', 'z', r'$P_l$')
        total_R22=Field.LoadFromFile('Stats_R22.dt_ev',["y"] , ["total"],'P_l', 'z', r'$P_l$')
        production_R22.settex(r'production')
        pressurestrain_R22.settex(r'pressure strain')
        dissipation_R22.settex(r'dissipation')
        turbulenttransport_R22.settex(r'turbulent transport$')
        pressuretransport_R22.settex(r'pressure transport')
        viscoustransport_R22.settex(r'viscous transport')
        total_R22.settex(r'total')    
        
        production_R33=Field.LoadFromFile('Stats_R33.dt_ev',["y"] , ["production"],'P_l', 'z', r'$P_l$')
        pressurestrain_R33=Field.LoadFromFile('Stats_R33.dt_ev',["y"] , ["pressurestrain"],'P_l', 'z', r'$P_l$')
        dissipation_R33=Field.LoadFromFile('Stats_R33.dt_ev',["y"] , ["dissipation"],'P_l', 'z', r'$P_l$')
        turbulenttransport_R33=Field.LoadFromFile('Stats_R33.dt_ev',["y"] , ["turbulenttransport"],'P_l', 'z', r'$P_l$')
        pressuretransport_R33=Field.LoadFromFile('Stats_R33.dt_ev',["y"] , ["pressuretransport"],'P_l', 'z', r'$P_l$')
        viscoustransport_R33=Field.LoadFromFile('Stats_R33.dt_ev',["y"] , ["viscoustransport"],'P_l', 'z', r'$P_l$')
        total_R33=Field.LoadFromFile('Stats_R33.dt_ev',["y"] , ["total"],'P_l', 'z', r'$P_l$')
        production_R33.settex(r'production')
        pressurestrain_R33.settex(r'pressure strain')
        dissipation_R33.settex(r'dissipation')
        turbulenttransport_R33.settex(r'turbulent transport$')
        pressuretransport_R33.settex(r'pressure transport')
        viscoustransport_R33.settex(r'viscous transport')
        total_R33.settex(r'total')         
        
        
        tracer([production_R11,pressurestrain_R11, dissipation_R11,turbulenttransport_R11,pressuretransport_R11,viscoustransport_R11, total_R11], "Rij_vrenam_11")    
        tracer([production_R12,pressurestrain_R12, dissipation_R12,turbulenttransport_R12,pressuretransport_R12,viscoustransport_R12, total_R12], "Rij_vrenam_13")    
        tracer([production_R22,pressurestrain_R22, dissipation_R22,turbulenttransport_R22,pressuretransport_R22,viscoustransport_R22, total_R22], "Rij_vrenam_33")    
        tracer([production_R33,pressurestrain_R33, dissipation_R33,turbulenttransport_R33,pressuretransport_R33,viscoustransport_R33, total_R33], "Rij_vrenam_22")    
       
        
        ############################## k-equation #################################
        
        advection=advection.contract(0, 'kk')
        production=production.contract(0, 'kk')
        redistribution=redistribution.contract(0, 'kk')
        dissipation=dissipation.contract(0, 'kk')
        diffusion=diffusion.contract(0, 'kk')
        diffusion_mol=diffusion_mol.contract(0, 'kk')
        residu_e=residu.contract(0, 'kk')
        
        advection.settex(r'advection $\frac{D\epsilon}{Dt}$')
        production.settex(r'production $2\overline{v^,_lv^,_k}\frac{\partial\overline{v_k}}{\partial x_l}$')
        redistribution.settex(r'redistribution $2\overline{p^,\frac{\partial v^,_k}{\partial x_k}}$')
        dissipation.settex(r'dissipation $2\nu\overline{\frac{\partial v^,_k}{\partial x_l} \frac{\partial v^,_k}{\partial x_l}}$')
        diffusion.settex(r'diffusion $\frac{\partial\overline{v^,_lv^,_kv^,_k}}{\partial x_l}+2\frac{\partial\overline{p^,v^,_k}}{\partial x_k}$')
        diffusion_mol.settex(r'diffusion moleculaire $\nu\frac{\partial^2\epsilon}{\partial^2 x_l}$')
        residu_e.settex(r'residu')
        
        tracer([advection, production, redistribution, dissipation, diffusion, diffusion_mol,residu_e], "epsilon_mono")

        
    return

def Rij_diphase(u, Rij, Rijk, pu, pdu, rho, dudu0, mu, vpain=None, uduain=None, duain=None, vaiN=None, aiN=None, vvaiN=None):
    """ 
    Calcule, trace et sauvegarde les différents termes de l'équation de transport des tensions de Reynolds pour un cas diphasique
    """
    a=u.tensorProduct(Rij)
    b=a.grad()
    nonlin1=b.contract('ij','kijk')
    gradRij=Rij.grad()
    gradu=u.grad()
    nonlin2=Rij.tensorProduct(gradu)
    nonlin2=nonlin2.contract('ij','kijk')
    nonlin3=nonlin2.transpo()
    nonlin=nonlin1+nonlin2+nonlin3
    gradRijk=Rijk.grad()
    nonlin4=gradRijk.contract('ij','kijk')
    
    ####### terme de pression ###############
    pijk=pu.grad()
    p1ijk=pijk.transpo()+pijk
    p2ijk=pdu.transpo()+pdu
    p=p2ijk-p1ijk
    ####### viscosite #########
    
    
    ddRij0=gradRij.grad()
    ddRij=ddRij0.contract('ij', 'ijkk')
    visc=ddRij*mu-dudu0*mu*2.0

    
    #################### label latex #############################
    
    advection=nonlin1*(-1.0)
    production=(nonlin2+nonlin3)*(-1.0)
    redistribution=p2ijk
    dissipation=dudu0*mu*2.0*(-1.0)
    diffusion=nonlin4*(-1.0)
    transpression=p1ijk*(-1.0)
    diffusion_mol=ddRij*mu
    advection.settex(r'advection $\frac{DR_{i,j}}{Dt}$')
    production.settex(r'production $\overline{v^,_lv^,_j}\frac{\partial\overline{v_i}}{\partial x_l}+\overline{v^,_lv^,_i}\frac{\partial\overline{v_j}}{\partial x_l}$')
    redistribution.settex(r'redistribution $\overline{p^,\left(\frac{\partial v^,_j}{\partial x_i}+\frac{\partial v^,_i}{\partial x_j}\right)}$')
    dissipation.settex(r'dissipation $2\nu\overline{\frac{\partial v^,_i}{\partial x_l} \frac{\partial v^,_j}{\partial x_l}}$')
    diffusion.settex(r'turbulent diffusion $\frac{\partial\overline{v^,_lv^,_iv^,_j}}{\partial x_l}$')
    transpression.settex(r'pressure diffusion $\frac{\partial}{\partial x_l}(\overline{p^,v^,_j}\delta_{il}+\overline{p^,v^,_i}\delta_{jl})$')
    diffusion_mol.settex(r'molecular diffusion  $\nu\frac{\partial^2\overline{v^,_iv^,_j}}{\partial^2 x_l}$')
    
    
    ###################### plot ##################################

    pmoyijk=vpain+vpain.transpo()
    
    visc1=uduain+uduain.transpo()
    b=u.tensorProduct(duain)
    c=gradu.tensorProduct(vaiN)
    c=c.contract('ij', 'ikjk')
    d=b+c
    visc2=d+d.transpo()
    e=u.tensorProduct(u).grad()
    f=e.tensorProduct(aiN)
    visc3=f.contract('ij', 'ijkk')
    viscinter1=visc1-visc2+visc3
     
    a=u.tensorProduct(vaiN)
    b=u.tensorProduct(u).tensorProduct(aiN)
    term1=vvaiN.grad()
    term2=a.grad()
    term3=b.grad()
    term1=term1.contract('ij', 'ijkk')
    term2=term2.contract('ij', 'ijkk')
    term2=term2+term2.transpo()
    term3=term3.contract('ij', 'ijkk')
    viscinter2=term1-term2+term3
    viscinter=viscinter2*mu+viscinter1*mu
    
    for i in [viscinter, p, visc, redistribution, dissipation, diffusion_mol, pmoyijk, transpression]:  
        try:
            i=i*(1/rho)
        except:
            i=i*rho.inv()
    
    residu=(visc-nonlin-nonlin4+p)*(-1.0)
    InterPlusRedi=residu+redistribution
    InterPlusRedi11=InterPlusRedi.compo(0,0)
    InterPlusRedi22=InterPlusRedi.compo(1,1)
    InterPlusRedi33=InterPlusRedi.compo(2,2)
    tracer([InterPlusRedi11,InterPlusRedi22,InterPlusRedi33], 'interredi')


    #residu=visc-nonlin-nonlin4+p-pmoyijk+viscinter
    residu.settex(r'interface $\langle\frac{\partial v^,_iv^,_jn_l\delta^I}{\partial x_l}\rangle+\langle\frac{\partial v^,_iv^,_j}{\partial x_l}n_l\delta^I\rangle-\langle p^,(v^,_jn_i\delta^I+v^,_in_j\delta^I)\rangle$ ')
    pmoyijk.settex(r'$\langle p^,(v^,_jn_i\delta^I+v^,_in_j\delta^I)\rangle$')
    viscinter.settex(r'$\langle\frac{\partial v^,_iv^,_jn_l\delta^I}{\partial x_l}\rangle+\langle\frac{\partial v^,_iv^,_j}{\partial x_l}n_l\delta^I\rangle$')


#     tracer([production, redistribution, dissipation, diffusion, transpression, diffusion_mol, residu], "Rij_diph_11", [0], log=True, xlim=[0.03,1], ylim=[-0.004,0.004])
#     tracer([production, redistribution, dissipation, diffusion, transpression, diffusion_mol, residu], "Rij_diph_12", [1], log=True, xlim=[0.03,1])
#     tracer([production, redistribution, dissipation, diffusion, transpression, diffusion_mol, residu], "Rij_diph_13", [2], log=True, xlim=[0.03,1], ylim=[-0.01,0.01])
#     tracer([production, redistribution, dissipation, diffusion, transpression, diffusion_mol, residu], "Rij_diph_21", [3], log=True, xlim=[0.03,1])
#     tracer([production, redistribution, dissipation, diffusion, transpression, diffusion_mol, residu], "Rij_diph_22", [4], log=True, xlim=[0.03,1], ylim=[-0.001,0.003])
#     tracer([production, redistribution, dissipation, diffusion, transpression, diffusion_mol, residu], "Rij_diph_23", [5], log=True, xlim=[0.03,1])
#     tracer([production, redistribution, dissipation, diffusion, transpression, diffusion_mol, residu], "Rij_diph_31", [6], log=True, xlim=[0.03,1])
#     tracer([production, redistribution, dissipation, diffusion, transpression, diffusion_mol, residu], "Rij_diph_32", [7], log=True, xlim=[0.03,1])
#     tracer([production, redistribution, dissipation, diffusion, transpression, diffusion_mol, residu], "Rij_diph_33", [8], log=True, xlim=[0.03,1], ylim=[-0.005,0.005])
    sym=[1,-1,-1,-1,1,-1,-1,-1,1]
    if 0:
	    production.symetriser(sym)
	    redistribution.symetriser(sym)
	    dissipation.symetriser(sym)
	    diffusion.symetriser(sym)
	    transpression.symetriser(sym)
	    diffusion_mol.symetriser(sym)
	    residu.symetriser(sym)
    
    listing=[production, redistribution, dissipation, diffusion, transpression, diffusion_mol, residu]
    #listing=[production, dissipation, diffusion, diffusion_mol]

    tracer(listing, "Rij_diph_11", [0], xlim=[0,0.35], ylim=[-0.01,0.01])
    tracer([residu*(-1.0), visc ,nonlin ,nonlin4, p], "Rij_diph_11_bonus", [0])
    tracer(listing, "Rij_diph_12", [1])
    tracer(listing, "Rij_diph_13", [2])
    tracer(listing, "Rij_diph_21", [3])
    tracer(listing, "Rij_diph_22", [4])
    tracer(listing, "Rij_diph_23", [5])
    tracer(listing, "Rij_diph_31", [6])
    tracer(listing, "Rij_diph_32", [7])
    tracer(listing, "Rij_diph_33", [8])

    res={'production':production, 'redistribution':redistribution, 'dissipation':dissipation, 'diffusion':diffusion, 'transpression':transpression,  'diffusion_mol':diffusion_mol,'residu':residu}

#     production=production.symetriser()
#     redistribution=redistribution.symetriser()
#     dissipation=dissipation.symetriser()
#     diffusion=diffusion.symetriser()
#     transpression=transpression.symetriser()
#     diffusion_mol=diffusion_mol.symetriser()
#     residu=residu.symetriser()
    

    if 0 : ## if 1 : pour mettre le blendin
        ###### Elliptic blending reynolds stress model ######
        #### donnees d'entree 
        g1=3.4
        g1_star=1.8
        g3=0.8
        g3_star=1.3
        g4=1.25
        g5=0.4
        C_eta=80.0
        C_l=0.133
        Prod=production.contract(0, 'kk')
        epsilon=dissipation.contract(0, 'kk') ### divise par 2 ou pas ?  
        k_cine=Rij.contract(0, 'kk')*0.5 #### negatif ou positif ?
        eps=epsilon
        eps_inv=eps.inv()
        k=k_cine
        k_inv=k.inv()
        P=Prod
        
        #eps=epsilon.mean()  ### *0.5 ?
        #k=k_cine.mean()
        #P=Prod.mean()
        ### matrice identite
        Ident=Rij.transpo()
        Ident._npa[:,:,:]=0.0
        Ident._npa[0,0,:]=1.0
        Ident._npa[1,1,:]=1.0
        Ident._npa[2,2,:]=1.0 
        tracer([Ident], 'k')
            
        #### calcul de phij_h (modele SSG)
        bij=Rij*(k_inv*0.5)-Ident*(0.33333333)
    
        
        gradu=u.grad()
        Sij=gradu+gradu.transpo()
        Sij=Sij*0.5
        Wij=gradu-gradu.transpo()
        Wij=Wij*0.5
        
        phij_h_11=bij*g1
        phij_h_12=bij*P*eps_inv*g1_star
        phij_h_1=phij_h_11+phij_h_12
        phij_h_1=phij_h_1*eps
        
        tracer([bij], 'bij')
        tracer([Rij], 'Rij')
        tracer([Sij], 'Sij')
        tracer([Wij], 'Wij')
            
        etape1=bij.tensorProduct(bij)
        etape2=etape1.contract(0,'klkl')
        etape4=etape2.sqrt()
        phij_h_21=Sij*k*g3
        phij_h_22=Sij*k*etape4*g3_star
        phij_h_2=phij_h_21-phij_h_22
        
    
        a=bij.tensorProduct(Sij)
        a=a.contract('ij','ikjk')
        b=bij.tensorProduct(Sij)
        b=b.contract(0,'klkl')
        b=b.mean()
        phij_h_3=a+a.transpo()-Ident*b*0.66666666
        phij_h_3=phij_h_3*k*g4
        
        c=bij.tensorProduct(Wij)
        c=c.contract('ij','ikjk')
        phij_h_4=c+c.transpo()
        phij_h_4=phij_h_4*k*g5
        
        phij_h=phij_h_2+phij_h_3+phij_h_4-phij_h_1
        phij_h.settex('SSG model')
        
    ##### calcul de phij_w (modele aux parois)
    ########## longueur caracteristique correlation vitesse
        L_1=k.power(1.5)*eps_inv*C_l
        eps2=eps.power(0.25)
        eps2=eps2.inv()
        L_2=eps2*mu**0.75*C_eta*C_l
        L=L_1+L_2
        for  i in range(len(L_1._npa[0,:])):
            L._npa[0,i]=max(L_1._npa[0,i],L_2._npa[0,i])
    #     L=L.symetriser()
        L2=L*L
        tracer([L2], 'L2', log=True)
        
    #####################################################
    ###### resolution de alpha par une methode explicite
    ####################################################
    #     alpha=L_1+L_2
    #     alpha._npa[0,:]=0
    #     dx=alpha._ax[0,10]-alpha._ax[0,9]
    #     for i in range(len(L_1._npa[0,:])):
    #         #L2=L2._npa[0,i-1]
    #         L2=L2.mean()
    #         if i==0:
    #             alpha._npa[0,i]=0
    #         elif i==1:
    #             alpha._npa[0,i]=-dx*dx/L2+alpha._npa[0,i-1]*dx*dx/L2+2*alpha._npa[0,i-1]            
    #         else :
    #             alpha._npa[0,i]=-dx*dx/L2+alpha._npa[0,i-1]*dx*dx/L2+2*alpha._npa[0,i-1]-alpha._npa[0,i-2]
    
    ####################################################
    ###### resolution alpha par une methode implicite
    ####################################################
    #     dx=L2._ax[0,10]-L2._ax[0,9]
    #     n=len(L_1._npa[0,:])
    #     n=n/2
    #     alpha=np.zeros(n*n).reshape(n,n)
    #     for i in range(n):
    #         for j in range(n):
    #             if i==j :
    #                 alpha[i,j]=-2.0-L2._npa[0,i] 
    #             if i==j+1:
    #                 alpha[i,j]=1.0
    #             if i==j-1:
    #                 alpha[i,j]=1.0  
    #                 
    #     ### CL au bord et au centre
    #     alpha[0,0]=1
    #     alpha[0,1]=0
    #     alpha[n-1,n-1]=1
    #     alpha[n-1,n-2]=0
    #    
    #     ### inversion matrice
    #     inver=np.linalg.inv(alpha)                
    #     
    #     #### remplissage vecteur b
    #     b=np.zeros(n)  
    #     for i in range(n):
    #         b[i]=dx*dx/L2._npa[0,i] 
    #     ### CL vecteur b
    #     b[0]=0
    #     b[n-1]=1
    #     print alpha
    #     print ''
    #     print ''
    #     print inver[:,:]
    #     print ''
    #     print ''    
    #     print b[:]
    #     
    #     ### calcul de la solution 
    #     solu=np.zeros(n)
    #     for i in range(n):
    #         for j in range(n):
    #             solu[i]=solu[i]+b[i]*inver[i,j]
    #     print solu
    #     
    #     ## met le resultat dans un objet a tracer
    #     alpha=L2+L2
    #     alpha._npa[0,0:n]=solu[:]
    #     tracer([alpha], 'alpha')   
    
    ########################################################
    ########################################################
    
    
    ########### solution analytique pour alpha avec coeff constant
    
    #    L=200.0
        L=L2.mean()
    #    L=1.0    
        B=math.exp(1/L)/(math.exp(-1/L)-math.exp(1/L))
        A=-(1+B)
        alpha=ScalarField('alpha', u._ax, u.labelRef)
        aa = len(u._npa[0,:])    
        alpha._npa=np.zeros((1,aa))
        alpha._npa[0,:]=0
        for i in range(0, len(alpha._npa[0,:])/2+1 ):
            alpha._npa[0,i]=A*math.exp(alpha._ax[0,i]/L)+B*math.exp(-alpha._ax[0,i]/L)+1
        ##### ON A TRICHE ICI. SYMETRISER LA SOLUTION POUR SATISFAIRE LES CONDITIONS
            alpha._npa[0,-i-1]=A*math.exp(alpha._ax[0,i]/L)+B*math.exp(-alpha._ax[0,i]/L)+1 
        tracer([alpha], 'alpha', log=True, xlim=[0,1])           
        n=alpha.grad()
        norme=max(n._npa[2,:])
        n=n*(1./norme)
        alpha3=alpha
        alpha3._npa[0,:]=alpha._npa[0,:]*alpha._npa[0,:]*alpha._npa[0,:]    
    
         
        phij_w_1=Rij.tensorProduct(n).tensorProduct(n)
        phij_w_1=phij_w_1.contract('ij','ikjk')
        phij_w_1=phij_w_1+phij_w_1.transpo()
        c=Rij.tensorProduct(n).tensorProduct(n)
        c=c.contract(0,'klkl')
        
        nn=n.tensorProduct(n)
        phij_w_2=c*nn+c*Ident
        phij_w_2=phij_w_2*(-0.5)
        phij_w=phij_w_1+phij_w_2
        phij_w=phij_w*eps*k_inv*(5.0)
    #    phij_w=phij_w*eps*k_inv*(-5.0)   
        #phij_w=phij_w*eps.mean()*k_inv.mean()*(-5.0)
        phij_w.settex('model paroi')
        phij_w_1.settex('1') 
        phij_w_2.settex('2')  
    #     ###### phi_tot 
         
        phij=phij_h+phij_w
        phij._npa[:,:,:]=0.0
        for i in range(len(phij._npa[0,0,:])):
            a=alpha3._npa[0,i]
            phij._npa[:,:,i]=phij_w._npa[:,:,i]*(1.0-a)+phij_h._npa[:,:,i]*a
        phij.settex('elliptic blending model')
        
        ########## 
     
        phi_Rij=Rij*eps*k_inv*5.0
        phi_Rij._npa[0,1,:]=0.0
        phi_Rij._npa[1,0,:]=0.0
    ##### nouvelle repartition de la trace propose
        phi_Rij._npa[0,0,:]=phi_Rij._npa[2,2,:]*(-2.0)
        phi_Rij._npa[1,1,:]=phi_Rij._npa[2,2,:]*(1.0)
        #### trace nulle classique
        phi_Rij_MANC=phi_Rij*(1.0)
        phi_Rij_MANC._npa[0,0,:]=phi_Rij._npa[2,2,:]*(-0.5)
        phi_Rij_MANC._npa[1,1,:]=phi_Rij._npa[2,2,:]*(-0.5) 
    ####### expression rigoureuse
        phi_Rij_RIG=phi_Rij*(1.0)
        phi_Rij_RIG._npa[0,0,:]=0.0
        phi_Rij_RIG._npa[1,1,:]=0.0  
        
        phi_Rij.settex(r'$\phi_{ij}^W$ avec $\phi_{11}^W=-2\phi_{33}^W$ et $\phi_{22}^W=\phi_{33}^W$')
        phi_Rij_RIG.settex(r'$\phi_{ij}^W$ avec $\phi_{11}^W=0$ et $\phi_{22}^W=0$')     
        phi_Rij_MANC.settex(r'$\phi_{ij}^W$ avec $\phi_{11}^W=-1/2\phi_{33}^W$ et $\phi_{22}^W=-1/2\phi_{33}^W$')
        
        #########################################
        tracer([redistribution, phij_h, phij, phi_Rij, phi_Rij_RIG, phi_Rij_MANC], "model_diph_11", [0], [-0.003,0.003], log=True, xlim=[0.0,1])
        tracer([redistribution, phij_h, phij, phi_Rij], "model_diph_12", [1], log=True, xlim=[0.0,1])
        tracer([redistribution, phij_h, phij, phi_Rij], "model_diph_13", [2], log=True, xlim=[0.0,1])
        tracer([redistribution, phij_h, phij, phi_Rij], "model_diph_21", [3], log=True, xlim=[0.0,1])
        tracer([redistribution, phij_h, phij, phi_Rij, phi_Rij_RIG, phi_Rij_MANC], "model_diph_22", [4], [-0.003,0.003], log=True, xlim=[0.0,1])
        tracer([redistribution, phij_h, phij, phi_Rij], "model_diph_23", [5], log=True, xlim=[0.0,1])
        tracer([redistribution, phij_h, phij, phi_Rij], "model_diph_31", [6], log=True, xlim=[0.0,1])
        tracer([redistribution, phij_h, phij, phi_Rij], "model_diph_32", [7], log=True, xlim=[0.0,1])
        tracer([redistribution, phij_h, phij, phi_Rij], "model_diph_33", [8], [-0.003,0.003], log=True, xlim=[0.0,1])   
        

        ############################## k-equation #################################
    
    advection=advection.contract(0, 'kk')
    production=production.contract(0, 'kk')
    redistribution=redistribution.contract(0, 'kk')
    dissipation=dissipation.contract(0, 'kk')
    diffusion=diffusion.contract(0, 'kk')
    diffusion_mol=diffusion_mol.contract(0, 'kk')
    residu_e=residu.contract(0, 'kk')
    transpression=transpression.contract(0, 'kk')
#    pmoyijk=pmoyijk.contract(0, 'kk')
#    viscinter=viscinter.contract(0, 'kk')
    
    advection.settex(r'advection $\frac{D\epsilon}{Dt}$')
    production.settex(r'production $2\overline{v^,_lv^,_k}\frac{\partial\overline{v_k}}{\partial x_l}$')
    redistribution.settex(r'redistribution $2\overline{p^,\frac{\partial v^,_k}{\partial x_k}}$')
    dissipation.settex(r'dissipation $2\nu\overline{\frac{\partial v^,_k}{\partial x_l} \frac{\partial v^,_k}{\partial x_l}}$')
    diffusion.settex(r'diffusion $\frac{\partial\overline{v^,_lv^,_kv^,_k}}{\partial x_l}+2\frac{\partial\overline{p^,v^,_k}}{\partial x_k}$')
    diffusion_mol.settex(r'diffusion moleculaire $\nu\frac{\partial^2\epsilon}{\partial^2 x_l}$')
    residu_e.settex(r'interface $\langle\frac{\partial v^,_kv^,_kn_l\delta^I}{\partial x_l}\rangle+\langle\frac{\partial v^,_kv^,_k}{\partial x_l}n_l\delta^I\rangle-\langle p^,(v^,_kn_k\delta^I)\rangle$ ')
    transpression.settex(r'pressure diffusion $2\frac{\partial\overline{p^,v^,_k}}{\partial x_k}$')
#    pmoyijk.settex(r'2$\langle p^,(v^,_kn_k\delta^I)\rangle$')
#    viscinter.settex(r'$\langle\frac{\partial v^,_kv^,_kn_l\delta^I}{\partial x_l}\rangle+\langle\frac{\partial v^,_kv^,_k}{\partial x_l}n_l\delta^I\rangle$')

    tracer([advection, production, redistribution, dissipation, diffusion, transpression, diffusion_mol,residu_e], "epsilon_diph")
    
    return res

def qdm_monophase(Run, phase, u, mu, pression, Rij, S, rho, aii=None, rho_av=None, g=None, sigma=None):
    
    """ 
    return a dictionary composed by all the qdm profil (turbulence, inertia; buoyancy etc...)
    @param the phase. Could be anything
    @param the velocity field (vector field)
    @param viscosity (scalar)
    @param pressure field (scalar field)
    @param Reynolds stress tensor (tensor field)
    @param source (scalar field) 
    @param volumic mass (scalar)
    """
    # bon avec rho \= 1 et mu \=1
    nonlin=u.tensorProduct(u)
    tau1=u.grad()
    tau2=tau1+tau1.transpo()
    tau3=tau2.grad()
    tau4=tau3.contract('i', 'ikk')
    gradp=pression.grad()
    divuu=nonlin.grad()
    divuu=divuu.contract('i', 'ikk')
    divRij=Rij.grad()
    divRij=divRij.contract('i', 'ikk')
    
    ### rho et mu ponderation
    
    divuu=divuu*rho
    divRij=divRij*rho
    divT=tau4*mu
    
    
    qdm=divT-divRij-gradp+S-divuu
    divuu.settex(r'inertia $\nabla .\left(\rho\overline{\mathbf{v}}\otimes\overline{\mathbf{v}}\right)$') 
    divRij.settex(r'turbulence $\nabla .\left(\overline{\rho\mathbf{vv}\otimes\mathbf{vv}}\right)$')     
    divT.settex(r'viscous stress $\nabla .\left(\nu\left(\nabla\overline{\mathbf{v}}+\nabla\overline{\mathbf{v}}^{t}\right)\right)$')  
    gradp.settex(r'pressure $\nabla\overline{P}$')
    S.settex(r'source $\beta$')
    qdm.settex('residue')
    tracer([divT, divuu, divRij, gradp,S, qdm], "mono_qdm_i", [0])
    tracer([divT, divuu, divRij, gradp,S, qdm], "mono_qdm_j", [1])
    tracer([divT, divuu, divRij, gradp,S, qdm], "mono_qdm_k", [2])
    Add_section(Run, 'qdm')
    Run['qdm']['inertie']=divuu
    Run['qdm']['turbulence']=divRij
    Run['qdm']['viscous']=divT
    Run['qdm']['pression']=gradp
    Run['qdm']['source']= S 
    Run['qdm']['residu']= qdm
    return 




def TestDebit(fstat1,fstat2, DT):
    """ 
    Vérifie le terme instationnaire de qdm entre deux fichiers stats pris sur des intervals de temps différents. 
    """
   
    ## rhou init
    rho_l=1.0
    rho_v=0.1
    fstat=fstat1
    I=Field.LoadFromFile(fstat,["coordonnee_K"] , ["I"],'I', 'z', r'I')
    Iv=(I-1.0)*(-1.0)
    u_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UI", "VI", "WI"],r'$U_l$', r'$z$', r'$u_l$')
    u_v=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UIv", "VIv", "WIv"],'U_v', 'z', r'$u_v$')
    rho1=I*rho_l
    rho2=Iv*rho_v
    rho=rho1+rho2
    u=u_v+u_l
    rhovinit=u*rho

    ## rhou fin
    rho_l=1.0
    rho_v=0.1
    fstat=fstat2
    I=Field.LoadFromFile(fstat,["coordonnee_K"] , ["I"],'I', 'z', r'I')
    Iv=(I-1.0)*(-1.0)
    u_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UI", "VI", "WI"],r'$U_l$', r'$z$', r'$u_l$')
    u_v=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UIv", "VIv", "WIv"],'U_v', 'z', r'$u_v$')
    rho1=I*rho_l
    rho2=Iv*rho_v
    rho=rho1+rho2
    u=u_v+u_l
    rhovfin=u*rho

    result=((rhovfin-rhovinit)*(-1.0/DT)).integ()

    return result

def qdm_diphase(phase, u, mu, pression, Rij, S, rho, aii=None, rho_av=None, g=None, sigma=None, I=None):
    """ 
    return a dictionary composed by all the qdm profil (turbulence, inertia; buoyancy etc...)
    @param the phase. Could be 'diph', 'liq' or 'gaz'
    @param the velocity field (vector field)
    @param viscosity scalar
    @param pressure field (scalar field)
    @param Reynolds stress tensor (tensor field)
    @param source (scalar field) 
    @param volumic mass (scalar)
    """
    
    
    
    

    
#     fstat='Stats.dt_ev'
#     ddUdxx=Field.LoadFromFile(fstat,["coordonnee_K"] , ["IddUdxx", "IddVdxx","IddWdxx"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
#     ddUdyy=Field.LoadFromFile(fstat,["coordonnee_K"] , ["IddUdyy", "IddVdyy","IddWdyy"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
#     ddUdzz=Field.LoadFromFile(fstat,["coordonnee_K"] , ["IddUdzz", "IddVdzz","IddWdzz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
#     tau4=ddUdxx+ddUdyy+ddUdzz
    tau1=u.grad()
    tau2=tau1+tau1.transpo()
    tau3=tau2.grad()
    tau4=tau3.contract('i', 'ikk')


    nonlin=u.tensorProduct(u)
    gradp=pression.grad()
    divuu=nonlin.grad()
    divuu=divuu.contract('i', 'ikk')
    divRij=Rij.grad()
    divRij=divRij.contract('i', 'ikk')
  
    
        ### rho et mu ponderation
    if phase=='diph':    
        divuu=divuu*rho
        divRij=rho*divRij
        divT=tau4*mu
        mutau=tau2*mu
    else:
        divuu=divuu*rho
        divRij=divRij*rho
        divT=tau4*mu
        mutau=tau2*mu       
  
    
    
    rhog_av=g*rho_av
    if phase=='diph':
        rhog_tmp=rho*g ## diphasique rho est un profil
    else:
        rhog_tmp=g*rho 
        
    if phase=='diph':
        rhog=rhog_tmp-rhog_av
        ai=aii*sigma
    else:
        rhog=rhog_tmp-rhog_av
        rhog=I*rhog
        S=I*S
        ai=divT+divRij*(-1.0)+gradp*(-1.0)+S+divuu*(-1.0)+rhog

    if 0:
	    divT.symetriser([1,1,1])
	    divRij.symetriser([1,1,1])
	    gradp.symetriser([1,1,1])
	    rhog.symetriser([1,1,1])
	    ai.symetriser([1,1,1])
    gradp=gradp*(-1.0)
    divRij=divRij*(-1.0)
    divuu=divuu*(-1.0)

    if phase=='diph':
    	insta=divT+divRij+gradp+S+divuu+rhog+ai
    else :
        insta=divT*0.0
    
    divuu.settex(r'inertia $-\nabla .\left(\rho\overline{\mathbf{v}}\otimes\overline{\mathbf{v}}\right)$') 
    divRij.settex(r'turbulence $-\nabla .\left(\overline{\rho\mathbf{v^,}\otimes\mathbf{v^,}}\right)$')     
    divT.settex(r'viscous stress $\nabla .\left(\nu\left(\nabla\overline{\mathbf{v}}+\nabla\overline{\mathbf{v}}^{t}\right)\right)$')  
    gradp.settex(r'pressure $-\nabla\overline{P}$')
    S.settex(r'source $-\beta\mathbf{e_x}$')
    rhog.settex(r'buoyancy $(\rho-\langle\rho\rangle)\mathbf{g}$')
    insta.settex(r'residu')
         
    if phase=='diph':
        ai.settex(r'residue $\sigma\kappa\mathbf{n}\delta^i$')
    else:
        ai.settex(r'residue $\langle\frac{p}{\rho}\nabla\chi-\nu\nabla\mathbf{v}\nabla\chi\rangle$')

    if 1 :
        tracer([divT, divRij, gradp, S,rhog, ai, insta], phase+'_qdm_i', [0])
        tracer([divT, divRij, gradp, S,rhog, ai, insta], phase+'_qdm_j', [1])
        tracer([divT, divRij, gradp, S,rhog, ai, insta], phase+'_qdm_k', [2])
    
    res={'inertie':divuu, 'turbulence':divRij, 'viscous':divT, 'pression':gradp, 'source':S, 'interface':ai, 'rhog':rhog }
    
    stress_divRij=divRij.integ()
    stress_divuu=divuu.integ()
    stress_divT=divT.integ()
    stress_gradp=gradp.integ()
    stress_rhog=rhog.integ()
    stress_ai=ai.integ()
    stress_S=S.integ()
    stress_S=stress_S*(-1.0)
    stress_insta=insta.integ()
    
    # regler constante integration S, divT, ai.
    #stress_divT=stress_divT.integ_const()
    #stress_divT=stress_divT*(-1.0)
    #stress_ai=stress_ai.integ_const()
    #stress_S=stress_S.integ_const()
    #stress_divRij=stress_divRij.integ_const()
    #stress_gradp=stress_gradp.integ_const()
    #stress_rhog=stress_rhog.integ_const()
    #stress_rhog=stress_rhog*(-1.0)


    stress_divuu.settex(r'inertia $\int^y_0\nabla .\left(\rho\overline{\mathbf{v}}\otimes\overline{\mathbf{v}}\right) dy$') 
    stress_divRij.settex(r'turbulence $\int^y_0\nabla .\left(\overline{\rho\mathbf{v^,}\otimes\mathbf{v^,}}\right) dy$')     
    stress_divT.settex(r'viscous stress $\int^y_0\nabla .\left(\nu\left(\nabla\overline{\mathbf{v}}+\nabla\overline{\mathbf{v}}^{t}\right)\right) dy$')  
    stress_gradp.settex(r'pressure $\int^y_0\nabla\overline{P}dy$')
    stress_S.settex(r'source $\int^y_0\beta dy$')
    stress_rhog.settex(r'buoyancy $\int^y_0(\rho-\langle\rho\rangle)\mathbf{g}dy$')     
    stress_ai.settex(r'interface $\int^y_0\sigma\kappa\mathbf{n}\delta^idy$')
    stress_insta.settex(r'residu')

    if phase=='diph':
        stress_ai.settex(r'interface (residue) $\int^y_0\sigma\kappa\mathbf{n}\delta^idy$')
    else:
        stress_ai.settex(r'interface (residue) $\int^y_0\langle\frac{p}{\rho}\nabla\chi-\nu\nabla\mathbf{v}\nabla\chi\rangle dy$')
    

    listing=[stress_divT, stress_divRij, stress_gradp, stress_S,stress_rhog, stress_ai, stress_insta]

    if 1 :   
        tracer(listing, phase+'_stress_i', [0], xlim=[0,stress_divT._ax[0,-1]], zerocentre=True)
        tracer(listing, phase+'_stress_j', [1], xlim=[0,stress_divT._ax[0,-1]])   
        tracer(listing, phase+'_stress_k', [2], xlim=[0,stress_divT._ax[0,-1]])  
        
         
    #### Bilan des pression
    db=0.3
    p_iner=u.compo(0)*u.compo(0)*rho*(1./3.)*(-1.0)
    p_turb=(Rij.compo(0,0)+Rij.compo(1,1)+Rij.compo(2,2))*rho*(1./3.)*(-1.0)
    p_pi=rhog.compo(0)
    p_sig=(ai.compo(0)+ai.compo(1)+ai.compo(2))*(1./3.)
    p_tau=(mutau.compo(0,0)+mutau.compo(1,1)+mutau.compo(2,2))*(1./3.)
    p=p_turb+p_sig+p_tau+p_pi
    preal=pression
    p.settex(r'$p_{tot}$')
    preal.settex(r'$p_{DNS}$')
    p_iner.settex(r'$p_{inertie}$')
    p_turb.settex(r'$p_{turbulence}$') 
    p_pi.settex(r'$p_{\pi}$') 
    p_sig.settex(r'$p_{\sigma}$') 
    p_tau.settex(r'$p_{\tau}$')
    tracer([preal, p_turb, p_tau], phase+'pres.png')
    tracer([preal, p, p_turb, p_sig, p_tau, p_pi], phase+'pression.png')
    
    
    return res
   
    
    
    
def MkMedNonLocal(Run, tol=5):
    """ 
    Routine probablement obsolete
    """

    I=Run['variable']['I']
    
    gradI=I.grad()
    listx=[]
    listym=[]
    listyp=[]
    listIndx=[]
    y_i1=0
    y_i2=0 
    y_ip=0
    for i in range(len(gradI._npa[2,:])):
        if i>2: 
            y_i2=gradI._npa[2,i-2]
        if i>1:
            y_i1=gradI._npa[2,i-1]
        if i<len(gradI._npa[2,:])-1:
            y_ip=gradI._npa[2,i+1]
        x_i0=gradI._ax[0,i]
        y_i0=gradI._npa[2,i]
        if i==0:
            listx.append(x_i0)
            listIndx.append(i)
            listym.append(y_i1)
            listyp.append(y_ip)
        if y_i0>0 and y_i1<0 and i>1:
            listx.append(x_i0)
            listIndx.append(i)
            listym.append(y_i1)
            listyp.append(y_ip)
        elif y_i0<0 and y_i1>0 and i>1:
            listx.append(x_i0)
            listIndx.append(i)
            listym.append(y_i1)
            listyp.append(y_ip)
        elif (abs(y_i0-y_i1))>tol*(abs(y_i1-y_i2)) and i>2:
            listx.append(x_i0)
            listIndx.append(i)
            listym.append(y_i1)
            listyp.append(y_ip)
        elif i==len(gradI._npa[2,:])-1:
            listx.append(x_i0)
            listIndx.append(i)
            listym.append(y_i1)
            listyp.append(y_ip)            
             

    
    alpha=Run['variable']['alpha']
    Mv=Run['qdm_gaz']['interface']
    Ml=Run['qdm_liq']['interface']
    sigma=Run['variable']['sigma']
    epsilon=Run['variable']['dnx']
    rc=Run['variable']['rc']
    db=Run['variable']['db']  

    MlNet=NetContribution(Ml, db, L=2.0)
    MvNet=NetContribution(Mv, db, L=2.0) 
    Msig=Ml+Mv    
    Msigapriori=Msig*0.0  
    gradInet=NetContribution(gradI, db, L=2.0)
    Inet=NetContribution(I, db, L=2.0)
    
    
### relève les demie-bulle cohérente 
    passe=False
    demiBubbles=[]
    for i in range(len(listx[:])-1):
        if passe:
            pass

        bulle=FictiveBubble()
        bulle.Start[0]=listIndx[i]
        bulle.Start[1]=I._npa[0,listIndx[i]]
        bulle.Start[2]=gradI._npa[0,listIndx[i]]
        bulle.Start[3]=I._ax[0,listIndx[i]]
        bulle.fin[0]=listIndx[i+1]
        bulle.fin[1]=I._npa[0,listIndx[i+1]]
        bulle.fin[2]=gradI._npa[0,listIndx[i+1]]
        bulle.fin[3]=I._ax[0,listIndx[i+1]]
        bulle._npa=I._npa[0,bulle.Start[0]:bulle.fin[0]]
        bulle._ax=I._ax[0,bulle.Start[0]:bulle.fin[0]]
        bulle.findPic(I)
        
        #### verifie la coherence du saut de taux de vide 
        if abs(max(abs(bulle._npa))-min(abs(bulle._npa)))<alpha:
            passe=True
            continue   
        else:
            bulle.amplitude=bulle.Start[1]-bulle.fin[1]
            print 'demie bulle fictive d amplitude ', bulle.amplitude, 'de',bulle.Start[0], 'a', bulle.fin[0]
            demiBubbles.append(bulle)     
            passe=False
            
    print 'Il y a en tout', len(demiBubbles), 'demie-bulles fictives soit' ,  len(demiBubbles)/2, 'bulles fictives'
    
### former les paires de bulles
    
    Bubbles=[]
    for i in demiBubbles:
        trouver=False
        for j in demiBubbles:
            if trouver:
                pass
            if i!=j and abs(i.amplitude+j.amplitude)<0.1*abs(i.amplitude):
                print 'deux demies bulles fictives ont formees une bulle fictive'
                ### demie bulle de gauche
                bulle=FictiveBubble()
                if i.amplitude>0.0:
                    bulle.Start[:]=i.Start[:]
                    bulle.fin[:]=j.fin[:]
                    bulle.flatStart[:]=i.fin[:]
                    bulle.flatEnd[:]=j.Start[:]                     
                ### demie bulle de gauche
                if j.amplitude>0.0:
                    bulle.Start[:]=j.Start[:]
                    bulle.fin[:]=i.fin[:]
                    bulle.flatStart[:]=j.fin[:]
                    bulle.flatEnd[:]=i.Start[:]              
                bulle._npa=I._npa[0,bulle.Start[0]:bulle.fin[0]]
                bulle._ax=I._ax[0,bulle.Start[0]:bulle.fin[0]]
                bulle.findPic(I) 
                Bubbles.append(bulle) 
                demiBubbles.remove(i)
                demiBubbles.remove(j)
                trouver=True                   
        

    for i in Bubbles[:]:
 
        db_real=Run['variable']['db']
        db=abs(i.Start[3]-i.fin[3])
        ponde=NetContribution(absoluteValue(I.grad()), db=db)
        ponderation=ponde._npa[2,int(abs(i.Start[0]+i.fin[0])/2)]*2.0
        
        db=abs(i.Start[3]-i.flatStart[3])+abs(i.fin[3]-i.flatEnd[3])
        
        if i.Start[3]<db/10.0:
            rc=(i.flatStart[3]-i.Start[3])-(1./7.)*db
            decal=(1./7.)*db
        else:
            rc=(i.flatStart[3]-i.Start[3])
            decal=0.0
            
        e=2.0*2.0*sigma*epsilon/db
        a=-(4*e)/(db**2-2*db*rc+rc**2)
        b=(4*e*db+4*rc*e)/(db**2-2*db*rc+rc**2)
        c=a*rc*db
        d=-(-2*a*db**3-3*b*db**2-6*c*db+2*a*rc**3+3*b*rc**2+6*c*rc)/(rc**3)
        f=-d*rc

    
        for j in range(len(Msigapriori._npa[:,0])):
            for k in range(len(Msigapriori._npa[j,:])):
                if j==0:
                    if MlNet._ax[0,k]<i.Start[3] or MlNet._ax[0,k]>i.fin[3] :
                        pass             
                    elif (MlNet._ax[0,k]-i.Start[3])<rc:
                        Msigapriori._npa[j,k]=(d*(MlNet._ax[0,k]-i.Start[3])**2+f*(MlNet._ax[0,k]-i.Start[3]))*ponderation
                        #Msigapriori._npa[j,-k]=-(d*(MlNet._ax[0,k]-i.Start[3])**2+f*(MlNet._ax[0,k]-i.Start[3]))*ponderation
                    elif MlNet._ax[0,k]>(i.flatEnd[3]-decal):   
                        #print   MlNet._ax[0,k],  i.flatEnd[3]   
                        Msigapriori._npa[j,k]=(a*(MlNet._ax[0,k]-(i.flatEnd[3]-rc-decal))**2+b*(MlNet._ax[0,k]-(i.flatEnd[3]-rc-decal))+c)*ponderation
                        #Msigapriori._npa[j,-k]=-((a*(MlNet._ax[0,k]-(i.flatEnd[3]-rc-decal))**2+b*(MlNet._ax[0,k]-(i.flatEnd[3]-rc-decal))+c)*ponderation)
                elif j==2:
                    Msigapriori._npa[j,k]=gradI._npa[j,k]*(4*sigma/db_real)             
                

    A=NetContribution(I.derivate().derivate(), 0.05)*0.002*0.1
    Msig.settex(r'$\mathbf{M_\sigma}$ DNS')
    A.settex(r'$\mathbf{M_\sigma}$ a priori')
    tracer([Msig, A], Run['variable']['name']+'MsigxApriori', [0], couleur=[0,0,2,1,1], markevry=Run['variable']['markevery'], xlim=[0,1.0])           


    Msig_num=Msigapriori
    f=MlNet*1.0
    Mil_apriori=MlNet*0.0
    Miv_apriori=MvNet*0.0
    Mlapriori=Mil_apriori*0.0
    Mvapriori=Mil_apriori*0.0    
    
    for k in range(len(listx[:])-2):   
        
        db=listx[k+2]-listx[k]
        L_cara=db/2  ### L_cara fonction de blending
    
        for j in range(len(Mil_apriori._npa[:,0])):
            for i in range(len(Mil_apriori._npa[j,:])):
                f._npa[j,i]=math.exp(-(MvNet._ax[0,i]/L_cara)**2)
                if MlNet._ax[0,i]<listx[k] or MlNet._ax[0,i]>listx[k+2] or Mil_apriori._npa[j,i]!=-0.0:
                    pass
                else:
                    Mil_apriori._npa[j,i]=((6.0*MlNet._npa[j,i]*(MlNet._ax[0,i]-listx[k]))/db)*(1.0-((MlNet._ax[0,i]-listx[k])/db))
                    Miv_apriori._npa[j,i]=((6.0*MvNet._npa[j,i]*(MvNet._ax[0,i]-listx[k]))/db)*(1.0-((MvNet._ax[0,i]-listx[k])/db))
                    
            
            

        for j in range(len(Mil_apriori._npa[:,0])):
            for i in range(len(Mil_apriori._npa[j,:])):
                if MlNet._ax[0,i]<listx[k] or MlNet._ax[0,i]>listx[k+2] or Mlapriori._npa[j,i]!=0.0:
                    pass
                else:              
                    Mlapriori._npa[j,i]=Mil_apriori._npa[j,i]+Msig_num._npa[j,i]*0.5*(1.0-f._npa[j,i])+f._npa[j,i]*Msig_num._npa[j,i]*1.0
                    Mvapriori._npa[j,i]=Miv_apriori._npa[j,i]+Msig_num._npa[j,i]*0.5*(1.0-f._npa[j,i])+f._npa[j,i]*Msig_num._npa[j,i]*0.0    
                        
                
            
    f.settex(r'Blending fonction')        
    Mlapriori.settex(r'$M_{l}$ a priori')
    Mvapriori.settex(r'$M_{v}$ a priori')    
    Mil_apriori.settex(r'$M_{il}$ a priori')   
    Miv_apriori.settex(r'$M_{iv}$ a priori')
    MlNet.settex(r'$F_{lD}$ DNS')
    MvNet.settex(r'$F_{vD}$ DNS')    
    Mv.settex(r'$M_{v}$ DNS')
    Ml.settex(r'$M_{l}$ DNS')        
    Msig_num.settex(r'$M_{\sigma}$ DNS')
    Run['variable']['Milapriori']=Mil_apriori
    Run['variable']['Mivapriori']=Miv_apriori    
    
    ### composante x
    Msig_num.settex(r'$M_{\sigma}$ a priori')
    Msig.settex(r'$M_{\sigma}$ DNS')
    tracer([Msig, Msig_num], 'Mx_apriori_sphe_a', [0], couleur=[0,0,3,1,1], markevry=Run['variable']['markevery'], xlim=[0,1.0])
    #tracer([Ml, Mlapriori, Mv, Mvapriori, Msig, Msig_num], Run['variable']['name']+'MsigxAprioriSphe', [0], couleur=[2,2,1,1,1], markevry=Run['variable']['markevery'], xlim=[0,1], ylim=[-0.8,0.8])
    tracer([Ml, Mlapriori,Mil_apriori, Mv, Mvapriori, Miv_apriori,Msig, Msig_num], Run['variable']['name']+'MsigxAprioriSphe', [0], couleur=[3,3,1,1,1], markevry=Run['variable']['markevery'], xlim=[0,1])
   
    
    
    MLiftApriori(Run)    
    MDragApriori(Run)

    return
    
    
def MkMed(Run, tol=5): 
    MLiftApriori(Run)    
    MDragApriori(Run)     
    return  
    
    
    
def MlVsForcesApriori(Run):
    """ 
    Comparaison du M_k avec les fermetures a priori de lift, drag et masse ajoutée. Routine probablement obsolète
    """

    
    MLiftApriori(Run, model='semiAnalytique')
    MDragApriori(Run)
    MAddedMassTurbApriori(Run)
    Iv=Run['variable']['I']
    pl=Run['variable']['p_l']
    Ml=Run['qdm_liq']['interface']-pl*Iv.grad()
    
    Ml.settex(r'$\mathbf{M_v^{Std}}=-\mathbf{M_l}-(p_l-p_w)\nabla\alpha_v$')
    ### force stationnaire
    Fl=Run['qdm_liq']['M_Lift']
    Fd=Run['qdm_liq']['M_drag_apriori']
    FMA=Run['qdm_liq']['M_AM_apriori']    
    
    ### Dispersion turbulente
    Fdt_iso=Run['qdm_liq']['MTD_theo_iso']
    Fdt_aniso=Run['qdm_liq']['MTD_theo_aniso']
    FMA_turb=Run['qdm_liq']['turb_added_mass_apriori']
    Ftotx=Fd+FMA
    Ftotz=Fl+Fdt_iso
    Ftotx.settex(r'$\mathbf{M^D}+\mathbf{M^{AM}}$')
    Ftotz.settex(r'$\mathbf{M^L}+\mathbf{M^{TD}}$')

    tracer([Ftotx, Ml.compo(0), Fd, FMA], Run['variable']['name']+'SommeForceVsMl_x', couleur=[2,1,1,1,1,1,1,1,1,1], xlim=[0,1.0])
    tracer([Ftotz, Ml.compo(2), Fl, Fdt_iso, FMA_turb.compo(2)], Run['variable']['name']+'SommeForceVsMl_z', couleur=[2,1,1,1,1,1,1,1,1,1], xlim=[0,1.0])
   
   
def MlVsForcesAposteriori(Run, fstat):
    """ 
    Comparaison du M_k avec les résultats a posteriori de neptune. fstat est ici un fichier de stats neptune.
    """

    I=Field.LoadFromFile(fstat,["y"] , ["I"],'dP', 'z', r'$\nabla P$')
    fd=Field.LoadFromFile(fstat,["y"] , ["fdx","fdz","fdy"],'dP', 'z', r'$\nabla P$')
    fl=Field.LoadFromFile(fstat,["y"] , ["flx","flz","fly"],'dP', 'z', r'$\nabla P$')
    fdt=Field.LoadFromFile(fstat,["y"] , ["fdtx","fdtz","fdty"],'dP', 'z', r'$\nabla P$')
    ftot=fl+fd+fdt
    Iv=(Run['variable']['I']-1.0)*(-1.0)
    pl=Run['variable']['p_l']
    Ml=Run['qdm_liq']['interface']-pl*Iv.grad()
    Ml.settex(r'$\mathbf{M_v^{Std}}=-\mathbf{M_l}-(p_l-p_w)\nabla\alpha_v$')
    fd.settex(r'$\mathbf{M_v^{D}}$ a posteriori')
    fl.settex(r'$\mathbf{M_v^{L}}$ a posteriori')
    fdt.settex(r'$\mathbf{M_v^{TD}}$ a posteriori')   
    ftot.settex(r'$\mathbf{M_v^{TD}}+\mathbf{M_v^{L}}+\mathbf{M_v^{D}}$')
    tracer([Ml, ftot, fd, fl, fdt], Run['variable']['name']+'ForcePostVsMl_x', [0], couleur=[2,1,1,1,1,1,1,1,1,1], xlim=[0,1.0])
    tracer([Ml, ftot, fd, fl, fdt], Run['variable']['name']+'ForcePostVsMl_z', [2], couleur=[2,1,1,1,1,1,1,1,1,1], xlim=[0,1.0])
    ### comparaison des coef de dispersion turbulente
    Iv.symetriser([1])
    Il_inv=Run['variable']['I'].inv()
    Il_nept_inv=((I-1.0)*(-1.0)).inv()
    MLiftApriori(Run)
    Rij=Run['Rij']['Rij']
    Klinv=(Rij.compo(0,0)*0.5+Rij.compo(1,1)*0.5+Rij.compo(2,2)*0.5).inv()
    R22inv=(Rij.compo(2,2)).inv()
    RnnSurRij=(Rij.compo(0,0)+Rij.compo(1,1)+Rij.compo(2,2))*(1./3.)*R22inv
    Klinv_proj=Klinv.ProjectionChampsDifferentesTaille(I)
    RnnSurRij_proj=RnnSurRij.ProjectionChampsDifferentesTaille(I)
    rhol=Run['variable']['rho_l']
    
    coef_TD_REF=Ml*((Iv.grad()).compo(2)).inv()*Il_inv*Klinv*(1/rhol)*RnnSurRij
    coef_TD_post=fdt*((I.grad()).compo(2)).inv()*Il_nept_inv*Klinv_proj*(1/rhol)*RnnSurRij_proj
    coef_TD_REF=((coef_TD_REF).compo(2)).abs() 
    coef_TD_post=((coef_TD_post).compo(2)).abs() 
    #tracer([Klinv,Klinv_proj], 'bla')
    #tracer([fdt*((I.grad()).compo(2)).inv()*(1/rhol)*Klinv_proj*(1/rhol) ], 'blablaa')
    #tracer([Ml*((Iv.grad()).compo(2)).inv()*(1/rhol)*Klinv*(1/rhol) ], 'blablaaaa', xlim=[0,1.0], ylim=[-10,10])   
    coef_TD_REF.settex(r'$C^{TD}_{REF}$')
    coef_TD_post.settex(r'$C^{TD}_{POST}$')
    tracer([coef_TD_REF, coef_TD_post], Run['variable']['name']+'CTD', couleur=[2,1,1,1,1,1], xlim=[0.0,1.0], ylim=[0,10]) 

def MLiftApriori(Run, model='semiAnalytique', ffstat=None):
    """ 
    Calcul du lift (et autres) a priori. Routine probablement obsolete. 
    """

    alpha=Run['variable']['alpha']
    fstat=Run['variable']['fstat']
    db=Run['variable']['db']    
    uu_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UUI", "UVI", "UWI","UVI", "VVI", "VWI","UWI", "VWI", "WWI"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    uu_v=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UUIv", "UVIv", "UWIv","UVIv", "VVIv", "VWIv","UWIv", "VWIv", "WWIv"],'UU_v', 'z', r'\langle(u\times u)\rangle)_v')
    u_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UI", "VI", "WI"],r'$U_l$', r'$z$', r'$u_l$')
    u_v=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UIv", "VIv", "WIv"],'U_v', 'z', r'$u_v$')
    F_lturb=NetContribution(Run['qdm_liq']['turbulence'].compo(2), db, 2.0)  
    F_pl=NetContribution(Run['qdm_liq']['pression'].compo(2), db, 2.0)
    F_pv=NetContribution(Run['qdm_gaz']['pression'].compo(2), db, 2.0) 
    I=Run['variable']['I']
    gradI=I.grad().compo(2)
    Ml=Run['qdm_liq']['interface'].compo(2)*(-1.0) 
    Mv=Run['qdm_gaz']['interface'].compo(2)*(-1.0)
    Msig_num=Run['qdm_diph']['interface'].compo(2)
    ul=Run['variable']['u_l'].compo(0) 
    uv=Run['variable']['u_v'].compo(0) 
    rho_l=Run['variable']['rho_l'] 
    rho_v=Run['variable']['rho_v']
    g=abs(Run['variable']['g'])
    sigma=Run['variable']['sigma']
    mu=Run['variable']['mu']
    Iv=(I-1.0)*(-1.0)
    uv=uv*Iv.inv(tol=0.0000001)
    ul=ul*I.inv()

    Eo=(abs(g)*(rho_l-rho_v)*db*db)/sigma
    dH=db*((1+0.163*Eo**(0.757))**(0.33333))
    Eod=(abs(g)*(rho_l-rho_v)*dH*dH)/sigma 
    Msig=Msig_num
    rotul=ul.grad()
    rotul=rotul.compo(2)
    vr=uv-ul
    ClEo=0.00105*Eo**3-0.0159*Eo**2-0.0204*Eo+0.474
    Cl=0.00105*Eod**3-0.0159*Eod**2-0.0204*Eod+0.474
#     print 'le diamètre modifié est', dH, 'au lieu de', db
#     print 'le Cl modifié est', Cl, 'au lieu de', ClEo
#     print 'le Eo modifié est', Eod, 'au lieu de', Eo
    Fll=NetContribution(Ml, db, 2.0)
    Flv=NetContribution(Mv, db, 2.0)

    
    Flcorr=Iv*vr*rotul*rho_l*Cl#*(-1.0)
    FlcorrEo=Iv*vr*rotul*rho_l*ClEo*(-1.0)
    Flcorr.settex(r'$\mathbf{M_v^L}=\alpha_v C_l\rho_l(\mathbf{u_v}-\mathbf{u_l})\frac{\partial u_{lx}}{\partial z}$')
    Run['qdm_liq']['M_Lift']=Flcorr
    Flref=Iv*vr*rotul*rho_l*0.27*(-1.0)
    F_lturb.settex(r'$\mathbf{F_{turb}}$') 
    F_pl.settex(r'$\mathbf{F_{pl}}$')
    F_pv.settex(r'$\mathbf{F_{pv}}$')

    Fll.settex(r'$F_{Ll}$ DNS')
    Flv.settex(r'$F_{Lv}$ DNS')    
    rotul.settex(r'$rotul$')
    vr.settex(r'$vr$')
    
    ######## turbulent dispersion force

    u=u_v+u_l
    uu=uu_v+uu_l 
    ufu=u.MatProduct(u)
    Rij=uu-ufu
    Kl=Rij.compo(0,0)*0.5+Rij.compo(1,1)*0.5+Rij.compo(2,2)*0.5
    rnn=Rij.compo(2,2)
    kk=(Rij.compo(0,0)+Rij.compo(1,1)+Rij.compo(2,2))*(1./3.)
    Run['Rij']['Rij']=Rij
    racine=math.sqrt(sigma/(abs(g)*(rho_l-rho_v)))
    cts=(2./3.)*db/racine
    f=(1-alpha)**(1.5)
    num=1+17.67*f**(6./7.)
    denum=18.67*f
    Cd=cts*(num/denum)**2
    print 'le coefficient de drag vaut', Cd
    gradII=I.grad().compo(2)
    F_tdf_old=(Kl*gradII)*(rho_l*Cd)*(-1.0)    
    
    #### coeff GTD calcul ####
    eij=Run['Rij']['dissipation']
    rhol=Run['variable']['rho_l']
    rhov=Run['variable']['rho_v']   
    eps=eij.compo(0,0)+eij.compo(1,1)+eij.compo(2,2)
    CAM=0.5
    Cmu=0.09
    beta=2.7

    b=(rhol+rhol*CAM)/(rhov+rhol*CAM)
    I=Run['variable']['I']
    Iinv=I.inv(1e-8)
    Iv=(Run['variable']['I']-1.0)*(-1.0)
    Ivinv=Iv.inv(1e-8)
    v_l=Iinv*Field.LoadFromFile(fstat,["coordonnee_K"] , ["UI"],'P_l', 'z', r'$P_l$')
    v_v=Ivinv*Field.LoadFromFile(fstat,["coordonnee_K"] , ["UIv"],'P_v', 'z', r'$P_v$')
    vr=v_v-v_l
    Fd=vr*Cd 
    tau12f=(Fd.inv())*((rhov/rhol)+CAM)
    tau12t=Kl*eps.inv()*(3./2.)*Cmu * ((vr*vr*Kl.inv()*beta+1.0).power(-0.5))
    eta=tau12t*(tau12f.inv())
    GTD=(Fd*tau12t-1.0)*(eta+b)*((eta+1.0).inv())+(eta+b*b)*((eta+1.0).inv())*CAM
    cdtur=GTD*(3./2.)
    id12=Kl*cdtur*(2./3.)*rho_l

    coef_turb=id12
    F_tdf_theo=gradII*coef_turb    
    Run['qdm_liq']['MTD_theo_iso']=F_tdf_theo*(1.0)
    Run['qdm_liq']['MTD_theo_iso'].settex(r'$\mathbf{M_v^{TD}}=-GTD\rho_lK_l\nabla\alpha_v$')
    
    coef_turb=id12*rnn*kk.inv()
    F_tdf_theo_aniso=gradII*coef_turb 
    Run['qdm_liq']['MTD_theo_aniso']=F_tdf_theo_aniso
    #### coeff GTD comme dans Neptune ####
    db=0.3
    coefic12=(Iv*2.0+1.0)*Iinv
    ic12=Iinv*coefic12*0.5*rhol
    cma=I*ic12*(1./rho_l)
    cpt=cma*2.0
    bb=(cpt+cma)*((cma+(rhov/rhol)).inv())
    if12=vr*Iinv*(3./4.)*rhol*Cd*(1/db)
    fd=I*if12*(1/rhol)
    t12f=(cma+(rhov/rhol))*fd.inv()
    ksue=Kl*eps.inv()
    vr2=vr*vr
    t1=Kl*(2./3.)
    zetar2=vr2*t1.inv()
    CBT=1.8
    t12t=ksue*((zetar2*CBT+1.0).sqrt()).inv()*Cmu
    etar=tau12t*(tau12f.inv())
    c12=ic12
    cdtur1=(if12*t12t*(1/rhol))*((bb+etar)*(etar+1.0).inv()+Iv*Iinv)
    cdtur2=((bb+etar)*(etar+1.0).inv())*I*(-1.0)-Iv*Iinv*Iinv
    cdtur3=c12*(1.0/rhol)*(bb.power(2.0)-bb)*(etar+1.0).inv()
    cdtur4=(bb.power(2.0)+etar)*(etar+1.0).inv()*rhov*(1.0/rhol)
    cdtur=cdtur1+cdtur2+cdtur3+cdtur4
    id12=cdtur*rhol*Kl*(2./3.)
    #coef_turb=id12*rnn*kk.inv()*(-1.0)
    coef_turb=id12
    F_tdf_nept_iso=gradII*coef_turb  
    coef_turb=id12*rnn*kk.inv()
    F_tdf_nept_aniso=gradII*coef_turb  
    F_tdf_theo.settex(r'Theo GTD isotrope a priori')
    F_tdf_theo_aniso.settex(r'Theo GTD anisotrope a priori')    
    F_tdf_nept_iso.settex(r'NCFD GTD isotrope a priori')
    F_tdf_nept_aniso.settex(r'NCFD GTD anisotrope a priori')
    F_tdf=F_tdf_nept_aniso
      
    Ml=Ml*(-1.0)
    Mv=Mv*(-1.0)  
    Ml.settex(r'Ref DNS')
    
    
    tracer([Ml , F_tdf_theo, F_tdf_theo_aniso, F_tdf_nept_iso, F_tdf_nept_aniso], Run['variable']['name']+'TurbDispApriori', couleur=[5,3,2], markevry=[0.03], xlim=[0,1.0])


#   
#   
#     if fstat!=None:
#         Ml=Ml.integ()
#         F_tdf_theo=F_tdf_theo.integ()
#         F_tdf_theo_aniso=F_tdf_theo_aniso.integ()
#         F_tdf_nept_iso=F_tdf_nept_iso.integ()
#         F_tdf_nept_aniso=F_tdf_nept_aniso.integ()        
#         fdt=Field.LoadFromFile(ffstat,["y"] , ["fdty"],'dP', 'z', r'$\nabla P$')
#         fdt_iso=Field.LoadFromFile('profil_QdM_iso.dat',["y"] , ["fdty"],'dP', 'z', r'$\nabla P$')
#         fdt_aniso=Field.LoadFromFile('profil_QdM_aniso.dat',["y"] , ["fdty"],'dP', 'z', r'$\nabla P$')       
#         fdt.settex(r'NCFD')
#         fdt_iso=fdt_iso.integ()
#         fdt_aniso=fdt_aniso.integ()   
#         fdt_aniso.settex(r'NCFD anisotrope posteriori')
#         fdt_iso.settex(r'NCFD isotrope posteriori')
#         tracer([Ml , F_tdf_theo, F_tdf_nept_iso, fdt_iso, F_tdf_theo_aniso, F_tdf_nept_aniso, fdt_aniso], Run['variable']['name']+'TurbDispApriori', couleur=[1,3,3,2,2], markevry=0.03, xlim=[0,1.0])
# 
#     
    
    F_tdf.settex(r'$F_L^{turb}=-GTD\rho_l K_l\mathbf{\nabla.\chi_l}$ a priori')
    F_tot=F_tdf+Flcorr
    F_tot_v=F_tot*(-1.0)
    F_tot.settex(r'$F_{Ll}=F_L^{turb}+F_L^{lam}$ a priori')
    F_tot_v.settex(r'$F_{Lv}=-F_{Ll}$ a priori')

#     tracer([Fll , F_tot, Flcorr, F_tdf, Flv, F_tot_v], Run['variable']['name']+'Lift', couleur=[4,3,2], markevry=0.03, xlim=[0,1.0])


    if model=='old':   
        Fl=NetContribution(Run['qdm_diph']['interface'], db, 2.0).compo(2)       
        anum=gradI*(0.09*sigma*4.0/db)+Fl*0.5
        adenom=gradI*(0.9*sigma*4.0/db)+Fl*0.5   
        
        a=anum*adenom.inv()
        coefMl=a*((a+1.0).inv())
        coefMv=(a+1.0).inv()   
        
        Mlapriori=Msig_num*coefMl
        Mvapriori=Msig_num*coefMv


    elif model=='semiAnalytique':
        dp=-10.0
        sigkap=Msig_num*0.0+2*sigma/db
        numA=sigkap*(dp/(dp-1.0))
        numB=Msig_num*(1.0/(2*db))
        denA=sigkap*(1.0/(dp-1.0))
        denB=Msig_num*(1.0/(2*db))
        num=numA-numB
        den=denA+denB
        numA.settex(r'numA')
        numB.settex(r'numB')
        denA.settex(r'denA')
        denB.settex(r'denB')       
        #tracer([numA, numB, denA, denB], 'blabla')
        #tracer([num, den], 'num', xlim=[0,1.0])
        #tracer([Msig_num], 'msig')
        rap=den.inv(1e-5)*num*(-1.0)
        coefMv=rap*((rap+1.0).inv())
        coefMl=(rap+1.0).inv()   
        
        #tracer([coefMv, coefMl], 'coef', xlim=[0,1.0])
        Mlapriori=Msig_num*coefMl
        Mvapriori=Msig_num*coefMv    
    
    elif model=='Analytique':
        Mlapriori=Msig_num*0.0
        Mvapriori=Msig_num*0.0    
        m1=Msig_num+sigkap
        m2=(Msig_num+sigkap).power(2)+Msig_num*sigkap*((8)/(dp-1.0))
        m3=m2.power(0.5)
        m4=m1+(Msig_num*m1.inv())*sigkap*((4)/(dp-1.0))
        Mvapriori=m1*(1./4.)+m4*(1./4.) #-sigkap*0.5
        Mlapriori=Msig_num-Mvapriori
 
    
    
    Mll=NetContribution(Ml, db, 2.0)
    Mvv=NetContribution(Mv, db, 2.0)
    Mlapriori=NetContribution(Mlapriori, db, 2.0) 
    Mvapriori=NetContribution(Mvapriori, db, 2.0)
    Mlapriori.settex(r'$M_{l}$ a priori')
    Mvapriori.settex(r'$M_{v}$ a priori')     
    Mll.settex(r'$F_{Ll}$ DNS')
    Mlapriori.settex(r'$F_{Ll}$ a priori')
    Mvapriori.settex(r'$F_{Lv}$ a priori')
    Mvv.settex(r'$F_{Lv}$ DNS')
#     tracer([Mll, Mlapriori, Mvv, Mvapriori, Msig, Msig_num], Run['variable']['name']+'liftModel', couleur=[2,2,3], markevry=Run['variable']['markevery'], xlim=[0,2.0])
    
    
    
    ### check decomposition pression 
    
        
    I=Run['variable']['I']
    Iinv=I.inv(1e-8)
    Iv=Run['variable']['I']-1.0
    Iv=Iv*(-1.0)
    Ivinv=Iv.inv(1e-8)
    gradI=I.grad()
    gradIv=Iv.grad()
    p_v=Ivinv*Run['variable']['p_v']
    p_l=Iinv*Run['variable']['p_l']
    p_v.settex(r'$p_v$')
    p_l.settex(r'$p_l$')
    kappa_sig=p_l*0.0+sigma*(4./db)
    kappa_sig.settex(r'$\kappa\sigma$')
    pvl=p_v-p_l
    pvl.settex(r'$p_v-p_l$')
    tracer([p_v, p_l, pvl, kappa_sig], 'pression', xlim=[0.5,1.5])
    dpl=p_l.grad()
    dI=I.grad()
    grad_ap_l=Run['variable']['p_l'].grad()
    a_gradp_l=dpl*I
    p_l_grada=p_l*dI
    a_gradp_v=p_v.grad()*Iv
    grad_ap_v=Run['variable']['p_v'].grad()
    p_v_grada=p_v*gradIv
    som_v=a_gradp_v+p_v_grada
    som_l=a_gradp_l+p_l_grada
    a_gradp_v.settex(r'$\alpha_v\nabla p_v$')
    a_gradp_l.settex(r'$\alpha_l\nabla p_l$')
    grad_ap_v.settex(r'$\nabla \alpha_v p_v$')
    grad_ap_l.settex(r'$\nabla \alpha_l p_l$')
    p_l_grada.settex(r'$p_l\nabla\alpha_l $')
    p_v_grada.settex(r'$p_v\nabla\alpha_v $')
    som_v.settex(r'$p_v\nabla\alpha_v +\alpha_v\nabla p_v $')
    som_l.settex(r'$p_l\nabla\alpha_l +\alpha_l\nabla p_l $')
    tab=[a_gradp_v, grad_ap_v, p_v_grada, som_v, a_gradp_l, grad_ap_l, p_l_grada, som_l]
    bb=[4,4]
#     tracer(tab, 'dpx', [0], xlim=[0.2,1.8], couleur=bb, markevry=0.05)
#     tracer(tab, 'dpy', [1], xlim=[0.2,1.8], couleur=bb, markevry=0.05)
#     tracer(tab, 'dpz', [2], xlim=[0.2,1.8], couleur=bb, markevry=0.05) 

    ### comparaison Neptune DNS
  

    M_nept=F_tdf
    #Mv_nept=M_nept+grad_ap_v.compo(2)-(Iv*dpl).compo(2)
    Mv_nept=M_nept+(p_v-p_l)*dI.compo(2)
    Mv_nept_sigkappa=M_nept+dI.compo(2)*(sigma*(4./db))
    Ml_nept=M_nept*(-1.0)#+p_l_grada.compo(2)
    Msig_nept=Ml_nept+Mv_nept_sigkappa
    Ml_new=Mlapriori
    Mv_new=Mvapriori
    Msig_new=Ml_new+Mv_new
    Ml=Run['qdm_liq']['interface'].compo(2)
    Mv=Run['qdm_gaz']['interface'].compo(2)
    Msig_num=Run['qdm_diph']['interface'].compo(2)
    Mv_nept.settex(r'$-\mathbf{M^{DT}}-\left(p_v-p_l\right)\nabla\alpha_v$') 
    Mv_nept_sigkappa.settex(r'$-\mathbf{M^{DT}}-\sigma\kappa\nabla\alpha_v$')
    Ml_nept.settex(r'$\mathbf{M^{DT}}$') 
    Msig_nept.settex(r'NEPTUNE $\mathbf{M_\sigma}$')
    Ml_new.settex(r'MODEL $\mathbf{M_l}$')
    Mv_new.settex(r'MODEL $\mathbf{M_v}$')
    Msig_new.settex(r'MODEL $\mathbf{M_\sigma}$')
    Ml.settex(r'$\mathbf{M_l}$ DNS') 
    Mv.settex(r'$\mathbf{M_v}$ DNS')
    Msig_num.settex(r'DNS $\mathbf{M_\sigma}$')
#     tracer([Ml, Ml_nept, Ml_new, Mv, Mv_nept, Mv_new, Mv_nept_sigkappa, Msig_num, Msig_nept, Msig_new], Run['variable']['name']+'CompaModelTransverse', couleur=[3,4,4], markevry=Run['variable']['markevery'], xlim=[0,1.0])
#     tracer([Ml, Ml_nept, Mv, Mv_nept_sigkappa, Msig_num, Msig_nept], Run['variable']['name']+'CompaModelTransverseSelect', couleur=[2,2,2], markevry=Run['variable']['markevery'], xlim=[0,1.0])
#     tracer([Ml, Ml_nept, Ml_new], Run['variable']['name']+'CompaModelTransverseMl', couleur=[3,4,4], markevry=Run['variable']['markevery'], xlim=[0,1.0])

    a=Ml.integ()
    aa=Ml_nept.integ()
    b=Mv.integ()
    bb=Mv_nept.integ()
    c=Msig_num.integ()
    aaa=Ml_new.integ()
    bbb=Mv_new.integ()
    a.settex(r'$\int_z\mathbf{M_l}$ DNS') 
    aa.settex(r'$\int_z\left(\mathbf{M^{DT}}+p_l\nabla\alpha_l\right)$') 
    b.settex(r'$\int_z\mathbf{M_v}$ DNS') 
    bb.settex(r'$\int_z\left(-\mathbf{M^{TD}}+p_v\nabla\alpha_v\right)$') 
    c.settex(r'$\int_z\mathbf{M_\sigma}=\int\mathbf{M_v}+\mathbf{M_l}$ DNS') 
    aaa.settex(r'modelisation a priori') 
    bbb.settex(r'modelisation a priori') 
       
            
            
            
#     tracer([a,aa,b,bb,c], Run['variable']['name']+'IntegralMkNept', sym=[1], couleur=[2,2,2], markevry=Run['variable']['markevery'])
#     tracer([a,aa,aaa,b,bb,bbb], Run['variable']['name']+'IntegralMkCompa', sym=[1], couleur=[3,3,1], markevry=Run['variable']['markevery'])
#     
    
    return
    
    

    
def MDragApriori(Run):   
    """ 
    Calcul du drag a priori
    """ 
    

    I=Run['variable']['I']
    alpha=Run['variable']['alpha']
    uv=Run['variable']['u_v'].compo(0)
    ul=Run['variable']['u_l'].compo(0)
    mu=Run['variable']['mu']
    db=Run['variable']['db']
    sigma=Run['variable']['sigma']
    rho_l=Run['variable']['rho_l']
    rho_v=Run['variable']['rho_v']
    g=abs(Run['variable']['g'])
    drho=Run['variable']['rho_l']-Run['variable']['rho_v']
    F1=Run['qdm_gaz']['inertie'].compo(0)
    F2=Run['qdm_liq']['inertie'].compo(0)     
    Eo=(g*rho_l*db*db)/sigma

    ##### BONUS #####
    if 1:
        Mv=Run['qdm_gaz']['interface']*(-1.0)
        Ml=Run['qdm_liq']['interface']*(-1.0)  
    ################### force de mass ajoutée ########
    Ca=-0.5*(1+2*alpha)/(1-alpha)*alpha*rho_l
    F0=F1-F2
    M_ma=F0*Ca    

    ################## force de drag #################
        
    E=1/(1.0+0.163*Eo**0.757)
    rd=(db/2.0)/(E**(1./3.))        
    Iv=(I-1.0)*(-1.0)
    uv=uv*Iv.inv(tol=0.0000001)
    ul=ul*I.inv()
    ur=uv-ul
    tracer([uv, ul, ur], 'u') 
    #vr=ur*1.0


    if 1:
   	 correlation=((I-1.0)*(-1.0)).grad()*ur*ur*rho_l*0.25
         correlation.symetriser([1,1,-1])
         Ml.symetriser([1,1,-1])
         Mv.symetriser([1,1,-1])
         Ml.settex(r'liquid')
         Mv.settex(r'vapour') 
 	 correlation.settex(r'Stuhmiller 1977')     
    	 tracer([Mv, Ml, correlation], 'pauchon_corr', [2], xlim=[0.0,1.0])
    ############## drift velocity vr #########
    uv=Run['variable']['u_v'].compo(0)
    ul=Run['variable']['u_l'].compo(0)  
    uv=uv*Iv.inv(tol=0.0000001)
    ul=ul*I.inv()
    V2j=(uv-ul)*((drho)/(alpha*rho_v+(1-alpha)*rho_l))
    jz=uv+ul+V2j
    chijz=(Iv*jz).mean()
    chifjz=(Iv.mean())*(jz.mean())
    Co=chijz/chifjz
    
    
    av=Iv.mean()
    vr=((1.0-Co*av)/(1.0-av))*(uv.mean())-Co*(ul.mean())
    #vr=uv*((1.0-Co*av)/(1.0-av))-ul*Co
    print 'la vitesse relative est', vr

    
    #vr=ur.mean()
    #vr=ur._npa[0,I.IndxInX(1.0)]/1.2
#     Nre=vr*((2*rho_l*rd)/mu)
#     Nmu=mu/(rho_l*sigma*((sigma/(g*drho))**0.5)**0.5)
    #Cd_choix=Nre*Nmu*(2./9.)**0.5*Cd_dpr 
    Cd_choix=uv*0.0+(4./3.)*rd*math.sqrt(g*drho/sigma)
    Ab=(4./3.)*3.1415*((db/2.0)**3)
    Ad=3.1415*rd**2
    Fd=Iv*Cd_choix*vr*vr*rho_l*(-Ad/(2.0*Ab))
    Fd.settex(r'$\mathbf{M_v^D}=-\frac{3}{8r_b}\alpha_vC_D\rho_l|\mathbf{v_r}|\mathbf{v_r}$')
    Run['qdm_liq']['M_drag_apriori']=Fd

                
    Miv_apriori=Run['qdm_gaz']['interface'].compo(0)*(-1.0)
    Miv_apriori.settex(r'Drag DNS')
    M_ma.settex(r'$\mathbf{M_v^{AM}}=C_A\rho_l\alpha_v(\frac{D\mathbf{u_v}}{dt}-\frac{D\mathbf{u_l}}{dt})$')
    Run['qdm_liq']['M_AM_apriori']=M_ma

#     tracer([Miv_apriori, Fd, M_ma], Run['variable']['name']+'DragModel', couleur=[4,1], markevry=Run['variable']['markevery'], xlim=[0,2.0])

    return


def PressureInterfacialApriori(Run): 
    sigma=Run['variable']['sigma']
    db=Run['variable']['db']
    fstat=Run['variable']['fstat']
    I=Run['variable']['I']
    rhol=Run['variable']['rho_l']
    rhov=Run['variable']['rho_v']
    MlDNSapriori=Run['qdm_liq']['interface']*(-1.0)  
    turb_l=Run['qdm_liq']['turbulence']    
    turb_v=Run['qdm_gaz']['turbulence']      
    kappatoile=calculerPression(Run)
    kappa=Run['variable']['kappa']
    Iinv=I.inv(1e-8)
    Iv=(Run['variable']['I']-1.0)*(-1.0)
    Ivinv=Iv.inv(1e-8)    
    pl=Iinv*Field.LoadFromFile(fstat,["coordonnee_K"] , ["PI"],'P_l', 'z', r'$P_l$')
    pv=Ivinv*Field.LoadFromFile(fstat,["coordonnee_K"] , ["PIv"],'P_v', 'z', r'$P_v$') 
    Msig=Run['qdm_diph']['interface']*(-1.0)
    g=Field.initgravity([Run['variable']['g'],0,0], 'g', I)
    viscousl=Run['qdm_liq']['viscous']
    viscousv=Run['qdm_gaz']['viscous']
    DIvDP=((pv-pl)*Iv).grad()
    MlDNSapriori.symetriser([1,1,-1])
    MlRANSapriori=MlDNSapriori-Msig+DIvDP
    MlPRANSapriori=DIvDP*Iv
    MlPDNSapriori=DIvDP*I*(-1.0)
    MlsigRANSapriori=Msig*Iv*(-1.0)
    MlsigDNSapriori=Msig*I

###### fermeture Stuhmiller
    uv=Run['variable']['u_v'].compo(0)
    ul=Run['variable']['u_l'].compo(0)
    uv=uv*Ivinv
    ul=ul*I.inv()
    ur=uv-ul
    MlRANSstuhmiller=(Iv.grad())*ur*ur*0.25*rhol*(-1.0)*Iv
    MlRANSstuhmiller.settex('pression stu')


###### fermeture antoine interfacial pressure

    rb=db/2.0
    gg=Run['variable']['g']
    rho_l=Run['variable']['rho_l']
    rho_v=Run['variable']['rho_v']
    Eo=(abs(gg)*rho_l*db*db)/sigma
    rbb=rb*((1+0.163*Eo**(0.757))**(0.33333))
    a=rbb
    b=(rb*rb*rb)/(a*a)
    mincurv=(b/(a*a))
    maxcurv=(a/(b*b))
    pi=3.1415
    drho=rho_l-rho_v
    Cd_choix=I*0.0+(4./3.)*a*math.sqrt(gg*(-1.0)*drho/sigma)
    e=math.sqrt(a**2-b**2)/a
    Ab=2.0*pi*a**2+pi*b**2/e*math.log((1+e)/(1-e))
    Ad=pi*a**2
    ketoi0=I*0.0+2.0*mincurv*sigma
    ketoimax=I*0.0+2.0*maxcurv*sigma
    C1=0.75
    model=False
    if model :
    	Md=Ivinv*Cd_choix*ur*ur*rho_l*Ad*pi*db*db*db*(1./6.)
    else:
	Md=Ivinv*MlDNSapriori.compo(0)*pi*db*db*db*(1./6.)

    plestimC1=Iv*(Md*(2.0/Ab)*(2*C1-1.0)+ketoi0-ketoimax)
    pvestimC1=Md*(2.0/Ab)*(2*C1-1.0)+ketoi0
    MlPantoine=((pvestimC1-plestimC1)*Iv).grad()*Iv
    MlPDNSantoine=((pvestimC1-plestimC1)*Iv).grad()*I*(-1.0)
    MlPantoine.settex('p antoine')

###### fermeture antoine surface tension

    kappaeff=4/rb-2*mincurv
    Mlsigantoine=(Iv.grad())*Iv*sigma*kappaeff*(-1.0)
    MlsigDNSantoine=(Iv.grad())*I*sigma*kappaeff
    Mlsigantoine.settex('$\sigma$ antoine')

###### fermeture totale antoine
    tot_antoine=Mlsigantoine+MlPantoine
    totDNS_antoine=MlsigDNSantoine+MlPDNSantoine
    tot_antoine.settex('tot antoine')

#### autre terme du bilan et residu
    flotta=g*I*Iv*(rhol-rhov)
    turb=turb_v*I-turb_l*Iv
    viscous=viscousv*I-viscousl*Iv   
    residuRANS=MlPRANSapriori+flotta+turb+MlsigRANSapriori+MlRANSapriori
    residuDNS=MlPDNSapriori+flotta+turb+MlsigDNSapriori-MlDNSapriori


    residuRANS.settex(r'residu')    
    residuDNS.settex(r'residu')   
    MlsigRANSapriori.settex(r"$\sigma\overline{\kappa\nabla\chi_k}$")
    MlPRANSapriori.settex(r'$\mathbf{M_l^P}$') 
    flotta.settex(r'$\alpha_v\alpha_l(\rho_l-\rho_v)\mathbf{g}$')
    turb.settex(r'$\mathbf{M_l^{\tau}}+\mathbf{M_l^{\Re}}$')
    MlRANSapriori.settex(r'$-\mathbf{M_l^{RANS}}$')
    MAddedMassTurbApriori(Run)
    fmaturb=Run['qdm_liq']['turb_added_mass_apriori']
    MlRANSapriori.symetriser([1,1,-1])


    ## pressure interfacial closure relation classique 

    tracer([MlPRANSapriori, MlPantoine, MlsigRANSapriori, Mlsigantoine, MlRANSapriori, tot_antoine, MlRANSstuhmiller, residuRANS], Run['variable']['name']+'equilibreRANS_antx', [0], couleur=[2,2,3,1,1,1,1], xlim=[0,1])
    tracer([MlPRANSapriori, MlPantoine, MlsigRANSapriori, Mlsigantoine, MlRANSapriori, tot_antoine, MlRANSstuhmiller, residuRANS] , Run['variable']['name']+'equilibreRANS_antz', [2], couleur=[2,2,3,1,1,1], xlim=[0,1], ylim=[-0.0005, 0.0005], sym=[-1.0])


    tracer([MlPDNSapriori, MlPDNSantoine, MlsigDNSapriori, MlsigDNSantoine, MlDNSapriori, totDNS_antoine, residuDNS], Run['variable']['name']+'equilibreDNS_antx', [0], couleur=[2,2,2,1,1,1,1], xlim=[0,1])
    tracer([MlPDNSapriori, MlPDNSantoine, MlsigDNSapriori, MlsigDNSantoine, MlDNSapriori, totDNS_antoine, residuDNS] , Run['variable']['name']+'equilibreDNS_antz', [2], couleur=[2,2,2,1,1,1], xlim=[0,1], ylim=[-0.005, 0.005], sym=[-1.0])


   

    return

def MAddedMassTurbApriori(Run):
    """ 
    Calcul du la force de masse ajoutée a priori
    """
    
    fstat=Run['variable']['fstat']
    uu_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UUI", "UVI", "UWI","UVI", "VVI", "VWI","UWI", "VWI", "WWI"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    uu_v=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UUIv", "UVIv", "UWIv","UVIv", "VVIv", "VWIv","UWIv", "VWIv", "WWIv"],'UU_v', 'z', r'\langle(u\times u)\rangle)_v')
    u_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UI", "VI", "WI"],r'$U_l$', r'$z$', r'$u_l$')
    u_v=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UIv", "VIv", "WIv"],'U_v', 'z', r'$u_v$')
    ufu_l=u_l.MatProduct(u_l)
    Rij_l=uu_l-ufu_l
    ufu_v=u_v.MatProduct(u_v)
    Rij_v=uu_v-ufu_v
    I=Run['variable']['I']
    rhol=Run['variable']['rho_l']
    Iv=(I-1.0)*(-1.0)
    ca=0.5
    Rij=Run['variable2']['Rij']
    Run['qdm_liq']['turb_added_mass_apriori']=((Iv*(Rij_v-Rij_l)).grad()).contract('i', 'ikk')*ca*rhol#*(-1.0)
    Run['qdm_liq']['turb_added_mass_apriori'].settex(r'$\mathbf{M^{AMT}}=-C_A\rho_l\nabla.(\alpha_v(R_{ij}^v-R_{ij}^l))$')
    



def SelectionStatOut(Run_imput):
    Il=Run_imput['variable']['I']
    Iv=(Run_imput['variable']['I']-1.0)*(-1.0)

    variables ={'vapour volume fraction': Iv}
    variables ={'liquid volume fraction': Il}
    variables['liquid velocity']=Run_imput['variable']['u_l']*(Il.inv(0.00000001))
    variables['vapour velocity']=Run_imput['variable']['u_v']*(Iv.inv(0.00000001))
    variables['liquid pressure']=Run_imput['variable']['p_l']*(Il.inv(0.00000001))
    variables['vapour pressure']=Run_imput['variable']['p_v']*(Iv.inv(0.00000001))
    variables['name']=Run_imput['variable']['name']
    variables['Rij total']=Run_imput['variable']['Rij']
    variables['Rij liquid']=Run_imput['variable']['Rij_l']
    variables['Rij vapour']=Run_imput['variable']['Rij_v']


    qdm_monofluid = {'viscous operator': Run_imput['qdm_diph']['viscous']}
    qdm_monofluid['buoyancy']=Run_imput['qdm_diph']['rhog']
    qdm_monofluid['source term']=Run_imput['qdm_diph']['source']
    qdm_monofluid['pressure gradient']=Run_imput['qdm_diph']['pression']
    qdm_monofluid['advection operator']=Run_imput['qdm_diph']['inertie']
    qdm_monofluid['Reynolds stresses']=Run_imput['qdm_diph']['turbulence']
    qdm_monofluid['surface tension source term']=Run_imput['qdm_diph']['interface']


    qdm_vapour = {'viscous operator': Run_imput['qdm_gaz']['viscous']}
    qdm_vapour['buoyancy']=Run_imput['qdm_gaz']['rhog']
    qdm_vapour['source term']=Run_imput['qdm_gaz']['source']
    qdm_vapour['pressure gradient']=Run_imput['qdm_gaz']['pression']
    qdm_vapour['advection operator']=Run_imput['qdm_gaz']['inertie']
    qdm_vapour['Reynolds stresses']=Run_imput['qdm_gaz']['turbulence']
    qdm_vapour['Interfacial force']=Run_imput['qdm_gaz']['interface']


    qdm_liquid = {'viscous operator': Run_imput['qdm_liq']['viscous']}
    qdm_liquid['buoyancy']=Run_imput['qdm_liq']['rhog']
    qdm_liquid['source term']=Run_imput['qdm_liq']['source']
    qdm_liquid['pressure gradient']=Run_imput['qdm_liq']['pression']
    qdm_liquid['advection operator']=Run_imput['qdm_liq']['inertie']
    qdm_liquid['Reynolds stresses']=Run_imput['qdm_liq']['turbulence']
    qdm_liquid['Interfacial force']=Run_imput['qdm_liq']['interface']

    Rij_transport_equation = {'redistribution tensor': Run_imput['Rij']['redistribution']} 
    Rij_transport_equation['molecular diffusion'] = Run_imput['Rij']['diffusion_mol']
    Rij_transport_equation['pressure diffusion'] = Run_imput['Rij']['transpression']
    Rij_transport_equation['production'] = Run_imput['Rij']['production']
    Rij_transport_equation['interfacial production'] = Run_imput['Rij']['residu']
    Rij_transport_equation['dissipation'] = Run_imput['Rij']['dissipation']


    Run={'Mean_variables' : variables, 'QdM_monofluid' : qdm_monofluid, 'QdM_vapour' : qdm_vapour, 'QdM_liquid' : qdm_liquid, 'Rij_transport_equation' : Rij_transport_equation}
    return Run

def outDataForStat(Run_imput):
    """ 
    print all the profils in a .txt files  
    @param a dictionnary  
    """   

    Run=SelectionStatOut(Run_imput)  
    L=len(Run['Mean_variables']['liquid volume fraction']._ax[0,:])

    for cle in Run.keys():
    	fichier = open(Run['Mean_variables']['name']+"_"+cle+".txt", 'w')
    	i=1
    	fichier.write('# '+'['+str(i)+']'+"coord_z"+" "+'\n')
	
        for cle2 in Run[cle].keys():
            if isinstance(Run[cle][cle2], ScalarField):
                i=i+1
                fichier.write('# '+'['+str(i)+']'+cle2+" "+'\n')
            if isinstance(Run[cle][cle2], VectorField):
                i=i+1
                fichier.write('# '+'['+str(i)+']'+cle2+'_x'+" "+'\n') 
                i=i+1  
                fichier.write('# '+'['+str(i)+']'+cle2+'_y'+" "+'\n')
                i=i+1
                fichier.write('# '+'['+str(i)+']'+cle2+'_z'+" "+'\n')
            if isinstance(Run[cle][cle2], Tensor2Field):
                i=i+1
                fichier.write('# '+'['+str(i)+']'+cle2+'_xx'+" "+'\n') 
                i=i+1  
                fichier.write('# '+'['+str(i)+']'+cle2+'_xy'+" "+'\n')
                i=i+1
                fichier.write('# '+'['+str(i)+']'+cle2+'_xz'+" "+'\n')   
                i=i+1
                fichier.write('# '+'['+str(i)+']'+cle2+'_yx'+" "+'\n')
                i=i+1   
                fichier.write('# '+'['+str(i)+']'+cle2+'_yy'+" "+'\n')
                i=i+1
                fichier.write('# '+'['+str(i)+']'+cle2+'_yz'+" "+'\n') 
                i=i+1
                fichier.write('# '+'['+str(i)+']'+cle2+'_zx'+" "+'\n')  
                i=i+1 
                fichier.write('# '+'['+str(i)+']'+cle2+'_zy'+" "+'\n')
                i=i+1
                fichier.write('# '+'['+str(i)+']'+cle2+'_zz'+" "+'\n')                                             
    	fichier.write('\n')
    	for j in range(L):
		fichier.write('\n')
		fichier.write(str(Run['Mean_variables']['liquid volume fraction']._ax[0,j])+'  ')

	    	for cle2 in Run[cle].keys():
	        	if isinstance(Run[cle][cle2], ScalarField):
	            	    fichier.write(str(Run[cle][cle2]._npa[0,j])+'  ')      
	        	if isinstance(Run[cle][cle2], VectorField):
			    fichier.write(str(Run[cle][cle2]._npa[0,j])+'  ')
			    fichier.write(str(Run[cle][cle2]._npa[1,j])+'  ')        
			    fichier.write(str(Run[cle][cle2]._npa[2,j])+'  ')
	        	if isinstance(Run[cle][cle2], Tensor2Field):
			    fichier.write(str(Run[cle][cle2]._npa[0,0,j])+'  ')
			    fichier.write(str(Run[cle][cle2]._npa[0,1,j])+'  ')        
			    fichier.write(str(Run[cle][cle2]._npa[0,2,j])+'  ')
			    fichier.write(str(Run[cle][cle2]._npa[1,0,j])+'  ')
			    fichier.write(str(Run[cle][cle2]._npa[1,1,j])+'  ')        
			    fichier.write(str(Run[cle][cle2]._npa[1,2,j])+'  ')
			    fichier.write(str(Run[cle][cle2]._npa[2,0,j])+'  ')
			    fichier.write(str(Run[cle][cle2]._npa[2,1,j])+'  ')        
			    fichier.write(str(Run[cle][cle2]._npa[2,2,j])+'  ')                                        
                                       
    
    return

def outDataForAutovnv(Run_imput):
    """ 
    print all the profils in a .txt files  
    @param a dictionnary  
    """   
    sym_tens=[1,1,-1,1,1,-1,-1,-1,1]
    Run=Run_imput
    for cle in Run.keys():
        for cle2 in Run[cle].keys():
            if isinstance(Run[cle][cle2], ScalarField):
                Run[cle][cle2].symetriser()
            elif isinstance(Run[cle][cle2], VectorField):
                Run[cle][cle2].symetriser()    
            elif isinstance(Run[cle][cle2], Tensor2Field):
                Run[cle][cle2].symetriser(sym_tens)    
    L=len(Run['variable']['I']._ax[0,:])/2
    fichier = open(Run['variable']['name']+"_out.txt", 'w')
    i=1
    fichier.write('# '+'['+str(i)+']'+"coord_z"+" "+'\n')
    for cle in Run.keys():
        for cle2 in Run[cle].keys():
            if isinstance(Run[cle][cle2], ScalarField):
                i=i+1
                fichier.write('# '+'['+str(i)+']'+cle+'_'+cle2+" "+'\n')
            if isinstance(Run[cle][cle2], VectorField):
                i=i+1
                fichier.write('# '+'['+str(i)+']'+cle+'_'+cle2+'_x'+" "+'\n') 
                i=i+1  
                fichier.write('# '+'['+str(i)+']'+cle+'_'+cle2+'_y'+" "+'\n')
                i=i+1
                fichier.write('# '+'['+str(i)+']'+cle+'_'+cle2+'_z'+" "+'\n')
            if isinstance(Run[cle][cle2], Tensor2Field):
                i=i+1
                fichier.write('# '+'['+str(i)+']'+cle+'_'+cle2+'_xx'+" "+'\n') 
                i=i+1  
                fichier.write('# '+'['+str(i)+']'+cle+'_'+cle2+'_xy'+" "+'\n')
                i=i+1
                fichier.write('# '+'['+str(i)+']'+cle+'_'+cle2+'_xz'+" "+'\n')   
                i=i+1
                fichier.write('# '+'['+str(i)+']'+cle+'_'+cle2+'_yx'+" "+'\n')
                i=i+1   
                fichier.write('# '+'['+str(i)+']'+cle+'_'+cle2+'_yy'+" "+'\n')
                i=i+1
                fichier.write('# '+'['+str(i)+']'+cle+'_'+cle2+'_yz'+" "+'\n') 
                i=i+1
                fichier.write('# '+'['+str(i)+']'+cle+'_'+cle2+'_zx'+" "+'\n')  
                i=i+1 
                fichier.write('# '+'['+str(i)+']'+cle+'_'+cle2+'_zy'+" "+'\n')
                i=i+1
                fichier.write('# '+'['+str(i)+']'+cle+'_'+cle2+'_zz'+" "+'\n')                                             
    fichier.write('\n')
    for j in range(L):
        fichier.write('\n')
        fichier.write(str(Run['variable']['I']._ax[0,j])+'  ')
        for cle in Run.keys():
            for cle2 in Run[cle].keys():
                if isinstance(Run[cle][cle2], ScalarField):
                    fichier.write(str(Run[cle][cle2]._npa[0,j])+'  ')      
                if isinstance(Run[cle][cle2], VectorField):
                    fichier.write(str(Run[cle][cle2]._npa[0,j])+'  ')
                    fichier.write(str(Run[cle][cle2]._npa[1,j])+'  ')        
                    fichier.write(str(Run[cle][cle2]._npa[2,j])+'  ')
                if isinstance(Run[cle][cle2], Tensor2Field):
                    fichier.write(str(Run[cle][cle2]._npa[0,0,j])+'  ')
                    fichier.write(str(Run[cle][cle2]._npa[0,1,j])+'  ')        
                    fichier.write(str(Run[cle][cle2]._npa[0,2,j])+'  ')
                    fichier.write(str(Run[cle][cle2]._npa[1,0,j])+'  ')
                    fichier.write(str(Run[cle][cle2]._npa[1,1,j])+'  ')        
                    fichier.write(str(Run[cle][cle2]._npa[1,2,j])+'  ')
                    fichier.write(str(Run[cle][cle2]._npa[2,0,j])+'  ')
                    fichier.write(str(Run[cle][cle2]._npa[2,1,j])+'  ')        
                    fichier.write(str(Run[cle][cle2]._npa[2,2,j])+'  ')                                        
                                       
    
    return


def compaTrustNept(Run, fstat):
    Run=calculerPression(Run)
    tab=["", "v"]
    tab2=["liq", "gaz"]
    for i in range(len(tab)):
        I=Field.LoadFromFile(fstat,["y"] , ["I"],'dP', 'z', r'$\nabla P$')
        P=Field.LoadFromFile(fstat,["y"] , ["P"],'P', 'z', r'$\nabla P$')
        dP=Field.LoadFromFile(fstat,["y"] , ["dPdx","dPdz","dPdy"],'dP', 'z', r'$\nabla P$')
        conv=Field.LoadFromFile(fstat,["y"] , ["convx"+tab[i],"convz"+tab[i],"convy"+tab[i]],'dP', 'z', r'$\nabla P$')
        drij=Field.LoadFromFile(fstat,["y"] , ["drijx"+tab[i],"drijz"+tab[i],"drijy"+tab[i]],'dP', 'z', r'$\nabla P$')
        dtau=Field.LoadFromFile(fstat,["y"] , ["dtaux"+tab[i],"dtauz"+tab[i],"dtauy"+tab[i]],'dP', 'z', r'$\nabla P$')
        fd=Field.LoadFromFile(fstat,["y"] , ["fdx","fdz","fdy"],'dP', 'z', r'$\nabla P$')
        fl=Field.LoadFromFile(fstat,["y"] , ["flx","flz","fly"],'dP', 'z', r'$\nabla P$')
        fdt=Field.LoadFromFile(fstat,["y"] , ["fdtx","fdtz","fdty"],'dP', 'z', r'$\nabla P$')
        fl=fl
        fdt=fdt
        
        #s=Field.LoadFromFile(fstat,["y"] , ["s","dPdx","dPdx"],'dP', 'z', r'$\nabla P$')
        #rhog=Field.LoadFromFile(fstat,["y"] , ["pi","dPdx","dPdx"],'dP', 'z', r'$\nabla P$')
        f=fd+fl+fdt
        dP=dP*(-1.0)
        if i==0:
            fd=fd*(-1.0)
            fdt=fdt*(-1.0)
            fl=fl*(-1.0)
            f=f*(-1.0)

        av_dpv= Run['variable']['av_dpv']*(-1.0)
        pv_dav= Run['variable']['pv_dav']*(-1.0)
        al_dpl= Run['variable']['al_dpl']*(-1.0)
        pl_dal= Run['variable']['pl_dal']*(-1.0)
        
        if i==0:   
            dPref=al_dpl
            dPref2=pl_dal
        elif i==1:
            dPref=av_dpv
            dPref2=pv_dav      
        dPPP=Run['qdm_'+tab2[i]]['pression']
        drijref=Run['qdm_'+tab2[i]]['turbulence']
        dtauref=Run['qdm_'+tab2[i]]['viscous']
        fref=Run['qdm_'+tab2[i]]['interface']
        sref=Run['qdm_'+tab2[i]]['source']
        s=dtau*0.0+Run['variable']['S']
        rhogref=Run['qdm_'+tab2[i]]['rhog']
        rho_l=Run['variable']['rho_l']
        rho_v=Run['variable']['rho_v']
        alpha=Run['variable']['alpha']
        g=Run['variable']['g']
        rhospa=alpha*rho_v+(1.0-alpha)*rho_l
        Il=(I-1.0)*(-1.0)
        Ilinv=Il.inv(1e-5)
        dPP=dP*Ilinv
        gradI=I.grad()
        gradIl=Il.grad()
        if i==0:
            I=(I-1.0)*(-1.0)
            dP=dP*Il
            drho=rho_l-rhospa
            dPP=gradIl*Run['variable']['sigma']*(4.0/0.3)
        elif i==1: 
            drho=rho_v-rhospa
            coef=I*Ilinv
            #dP=dP*coef
            dP=dP*I
            dPP=gradI*Run['variable']['sigma']*(4.0/0.3)
        s=I*s
        s._npa[2,:]=0.0
        s._npa[1,:]=0.0
        #dPP=dPP*0.0
        rhog=dP*0.0
        rhog._npa[0,:]=I._npa[0,:]*drho*g  
        fdt.settex(r'RANS dispersion')
        fl.settex(r'RANS lift') 
        fd.settex(r'RANS drag')
        dP.settex(r'RANS $\alpha\nabla p$')
        dPref.settex(r'DNS $\alpha\nabla p$')
        dPref2.settex(r'DNS $p \nabla\alpha$')        
        drij.settex(r'RANS turbulence')
        drijref.settex(r'DNS turbulence')
        dtau.settex(r'RANS viscous')
        dtauref.settex(r'DNS viscous')
        f.settex(r'RANS force')
        fref.settex(r'DNS force')
        s.settex(r'RANS source')
        sref.settex(r'DNS source')
        rhog.settex(r'RANS buoyancy')
        rhogref.settex (r'DNS buoyancy')
        residu=dP+dtau+drij+rhog-f+s-conv

        
        residuref=dPref+dPref2+dtauref+drijref+rhogref+sref-fref
        residuref.settex (r'DNS residu')
        residu.settex (r'RANS residu')
        dPPP.settex (r'DNS $\nabla\alpha p$')
        tracer([dPref, dP, dPref2, drijref,drij, dtauref, dtau,  fref, f, fd,  sref, s, rhogref, rhog, residuref, residu, conv], 'CompaQdMDNSNeptune_x_'+tab2[i], [0] ,couleur=[2,1,2,2,3,2,2,2,2], xlim=[0.0,1.0])
        #tracer([dPref, dP, dPref2, drijref,drij, dtauref, dtau,  fref, f, Run['qdm_diph']['interface'], sref, s, rhogref, rhog, residuref, residu], 'CompaQdMDNSNeptune_y_'+tab2[i], [1] ,couleur=[2,1,2,2,3,2,2,2,2,2], xlim=[0.0,1.0])
        tracer([dPref, dP, dPref2, drijref,drij, dtauref, dtau,  fref, f, sref, s, rhogref, rhog, residuref, residu, conv ], 'CompaQdMDNSNeptune_z_'+tab2[i], [2] ,couleur=[2,1,2,2,2,2,2,2,2,2], xlim=[0.0,1.0])      
        tracer([f, fdt, fd, fl], 'forcesNept_x_'+tab2[i], [0])
        tracer([f, fdt, fd, fl], 'forcesNept_y_'+tab2[i], [1])
        tracer([f, fdt, fd, fl], 'forcesNept_z_'+tab2[i], [2])
        tracer([fdt], 'dturb_z_'+tab2[i], [2])      
        
        dP=dP.integ()
        dPref=dPref.integ()
        dPref2=dPref2.integ()
        drij=drij.integ()
        drijref=drijref.integ()
        dtau=dtau.integ()
        dtauref=dtauref.integ()
        f=f.integ()
        fref=fref.integ()
        s=s.integ()
        sref=sref.integ()
        rhog=rhog.integ()
        rhogref=rhogref.integ()
        fdt=fdt.integ()
        fd=fd.integ()
        residuref=residuref.integ()
        residu=residu.integ()  
        conv=conv.integ()            
        tracer([dPref,dP,  dPref2, drijref,drij, dtauref, dtau,  fref, f, fd, sref, s, rhogref, rhog, residuref, residu, conv], 'CompaStressDNSNeptune_x_'+tab2[i], [0] ,couleur=[2,1,2,2,3,2,2,2,2], xlim=[0.0,1.0])
        #tracer([dPref,dP,  dPref2, drijref,drij, dtauref, dtau,  fref, f, sref, s, rhogref, rhog, residuref, residu], 'CompaStressDNSNeptune_y_'+tab2[i], [1] ,couleur=[2,1,2,2,2,2,2,2,2], xlim=[0.0,1.0])
        tracer([dPref,dP,  dPref2, drijref,drij,  fref, f, fdt, residuref, residu, conv], 'CompaStressDNSNeptune_z_'+tab2[i], [2] ,couleur=[2,1,2,3,2,2,2,2], xlim=[0.0,1.0])
        MLiftApriori(Run, model='semiAnalytique', ffstat=fstat)


def EquationVrDNS(Run):
    sigma=Run['variable']['sigma']
    db=Run['variable']['db']
    fstat=Run['variable']['fstat']
    I=Run['variable']['I']
    rhol=Run['variable']['rho_l']
    rhov=Run['variable']['rho_v']
    Ml=Run['qdm_liq']['interface']*(-1.0)  
    turb_l=Run['qdm_liq']['turbulence']    
    turb_v=Run['qdm_gaz']['turbulence']      
    Iinv=I.inv(1e-8)
    Iv=(Run['variable']['I']-1.0)*(-1.0)
    Ivinv=Iv.inv(1e-8)    
    pl=Iinv*Field.LoadFromFile(fstat,["coordonnee_K"] , ["PI"],'P_l', 'z', r'$P_l$')
    pv=Ivinv*Field.LoadFromFile(fstat,["coordonnee_K"] , ["PIv"],'P_v', 'z', r'$P_v$') 
    fstat="Stats_kai.dt_ev"
    kai=(Field.LoadFromFile(fstat,["coordonnee_K"] , ["kai"],'P_v', 'z', r'$P_v$'))
    ai=(Field.LoadFromFile(fstat,["coordonnee_K"] , ["AI"],'P_v', 'z', r'$P_v$'))
    fstat=Run['variable']['fstat']


    g=Field.initgravity([Run['variable']['g'],0,0], 'g', I)
    pressionl=Run['qdm_liq']['pression']
    pressionv=Run['qdm_gaz']['pression']
    viscousl=Run['qdm_liq']['viscous']
    viscousv=Run['qdm_gaz']['viscous']  
    Msig=Run['qdm_diph']['interface']*(-1.0)
    sigkappa=sigma*(4./db)
    
    pression=((pv-pl)*I*Iv).grad()+pl*Iv.grad()+Iv*(pv-pl)*Iv.grad()
    pression2=(I*Iv).grad()*sigkappa+pl*Iv.grad()+Iv*sigkappa*Iv.grad()
    pression3=pl*Iv.grad()+I*Msig
    
    
#### autre terme du bilan et residu
    
    #interface=Msig*I*(-1.0)
    interface=Msig*I
    flotta=g*I*Iv*(rhol-rhov)
    turb=turb_v*I-turb_l*Iv
    viscous=viscousv*I-viscousl*Iv   
     
    interface.settex(r"$\sigma\overline{\kappa\nabla\chi_k}$")
    pression.settex(r'$\mathbf{M_l^P}$')
    pression2.settex(r'Pression hypothese $p_v-p_l=\sigma\kappa=Cte$')
    pression3.settex(r'Pression hypothese particulaire')    
    flotta.settex(r'$\alpha_v\alpha_l(\rho_l-\rho_v)\mathbf{g}$')
    turb.settex(r'$\mathbf{M_l^{\tau}}+\mathbf{M_l^{\Re}}$')
    Ml.settex(r'$\mathbf{M_l^{RANS}}$')
    Ml.symetriser([1,1,-1])


    ## pressure interfacial closure relation classique 


    tracer([I], Run['variable']['name']+'I')
    tracer([ Ml, pression, interface, flotta, turb], Run['variable']['name']+'equilibreDNSx', [0], couleur=[1,1,1,1,1,1,1], xlim=[0,1])
    tracer([ Ml, pression, interface, flotta, turb] , Run['variable']['name']+'equilibreDNSz', [2], couleur=[1,1,1,1,1,1,1], xlim=[0,1], sym=[-1.0])

    #tracer([flotta, Ml, turb, dpkappaterm, kappatermdns, residu4], Run['variable']['name']+'equilibreRANSxNew', [0], couleur=[1,1,1,1,1,1,3,1], xlim=[0,1])
    #tracer([flotta, Ml, turb, kappaterm1, kappaterm2+kappaterm3, residu4], Run['variable']['name']+'equilibreRANSzNew', [2], couleur=[1,1,1,1,1,1,1], xlim=[0,1], ylim=[-0.005, 0.005], sym=[-1.0])

    return





def EquationVrDNS_old(Run):
    sigma=Run['variable']['sigma']
    db=Run['variable']['db']
    fstat=Run['variable']['fstat']
    I=Run['variable']['I']
    #tracer([(((I-1.)*(-1.)).grad()).compo(2)], 'gradI')
    rhol=Run['variable']['rho_l']
    rhov=Run['variable']['rho_v']
    Ml=Run['qdm_liq']['interface']*(-1.0)  
    turb_l=Run['qdm_liq']['turbulence']*(-1.)   
    turb_v=Run['qdm_gaz']['turbulence']*(-1.)
    kappatoile=calculerPression(Run)      
    kappa=Run['variable']['kappa']
    Iinv=I.inv(1e-8)
    Iv=(Run['variable']['I']-1.0)*(-1.0)
    Ivinv=Iv.inv(1e-8)    
    pl=Iinv*Field.LoadFromFile(fstat,["coordonnee_K"] , ["PI"],'P_l', 'z', r'$P_l$')
    pv=Ivinv*Field.LoadFromFile(fstat,["coordonnee_K"] , ["PIv"],'P_v', 'z', r'$P_v$') 
    try :
	    kai=(Field.LoadFromFile(fstat,["coordonnee_K"] , ["kai"],'P_v', 'z', r'$P_v$'))
	    ai=(Field.LoadFromFile(fstat,["coordonnee_K"] , ["AI"],'P_v', 'z', r'$P_v$'))
    except :
	    print "k ai non post-traite"


    fstat=Run['variable']['fstat']

    try :
	    kdns=kai*ai.inv(0.000000000000001)
	    #tracer([kai, kdns], 'kappa')
    except :
	    print "k ai non post-traite"

    IDpl=Run['variable']['al_dpl']
    plDI=Run['variable']['pl_dal']
    g=Field.initgravity([Run['variable']['g'],0,0], 'g', I)
    pressionl=Run['qdm_liq']['pression']
    pressionv=Run['qdm_gaz']['pression']
    #viscousl=Run['qdm_liq']['viscous']
    #viscousv=Run['qdm_gaz']['viscous']  

    ul=Run['variable']['u_l']
    uv=Run['variable']['u_v']
    mu=Run['variable']['mu']

    tau1_l=ul.grad()
    tau2_l=tau1_l+tau1_l.transpo()
    tau3_l=tau2_l.grad()
    tau4_l=tau3_l.contract('i', 'ikk')
    tau1_v=uv.grad()
    tau2_v=tau1_v+tau1_v.transpo()
    tau3_v=tau2_v.grad()
    tau4_v=tau3_v.contract('i', 'ikk')
    viscousl=tau4_l*mu
    viscousv=tau4_v*mu


    Msig=Run['qdm_diph']['interface']*(-1.0)
    sigkappa=sigma*(4./db)
    
    #pression=((pv-pl)*I*Iv).grad()+pl*Iv.grad()+Iv*(pv-pl)*Iv.grad()
    pression=((pv-pl)*Iv).grad()*I-pl*(Iv.grad())
    pression_v_1=((pv-pl)*Iv).grad()*Iv*(-1.)
    pression_v_2=pl*(Iv.grad())*(-1.)

    pression_gradi=((pv-pl)*I+pl)*(Iv.grad())*(-1.)
    pression_nogradi=I*Iv*((pv-pl).grad())*(-1.)

    pression2=(I*Iv).grad()*sigkappa+pl*Iv.grad()+Iv*sigkappa*Iv.grad()
    pression3=pl*Iv.grad()+I*Msig
    
    
    #### construction d'un gradient de curvature
    db=Run['variable']['db']

    #### Triche 
    #db=0.35 

    gg=Run['variable']['g']    
    rho_l=Run['variable']['rho_l'] 
    rho_v=Run['variable']['rho_v']     
    rb=db/2.0
    ### Pour le demi grand axe, on peut utiliser la correlation de wellek (attention en réalité il ya le pitching etc...)
    wellek=True
    if wellek:
        Eo=(abs(gg)*rho_l*db*db)/sigma
        rbb=rb*((1+0.163*Eo**(0.757))**(0.33333))
        a=rbb
        print 'Eo=', Eo, 'demi grand axe=', a        
    else:
        a=0.18 #demi grand axe visible (à la main)
    bubble_wall=0.05## limite en deca de laquelle le taux de vide est negligeable (à la main)
    b=(rb*rb*rb)/(a*a)
    pi=3.1415
    ttt=I*0.0
    cost=I*0.0    
    sint=I*0.0   
    for i in range(len(ttt._npa[0,:])):
        if ttt._ax[0,i]>=bubble_wall and ttt._ax[0,i]<=a*2+bubble_wall:
            ttt._npa[0,i]=(pi/(a*2))*ttt._ax[0,i]-(pi*(bubble_wall/(2*a)))
            cost._npa[0,i]=math.cos(ttt._npa[0,i])
            sint._npa[0,i]=math.sin(ttt._npa[0,i])      
            
    curvature=((cost*cost*b*b+sint*sint*a*a).power(1.5)).inv(0.0000001)*a*b
    mincurv=((b/(a*a))+(4.0/db))/2.0
    print 'maxcurv=', (b/(a*a)), 'mincurv=', (a/(b*b))
    k_av=curvature.MoyenneGlissante(db=a*2.0, maxcurv=(a/(b*b)), wall=bubble_wall)


    try :
	    kdns.symetriser([1.0])
	    k_db=kdns*0.0+(4./db)
    except :
	    print "k ai non post-traite"

    #k_av=NetContribution(curvature, a)
    dkappa=k_av.derivate()
    kappaterm=Msig*(0.0) 
    kappaterm._npa[2,:]=I._npa[0,:]*Iv._npa[0,:]*sigma*dkappa._npa[0,:]
    kappaterm_tronc=kappaterm.tronc(val=0.0)
    kappatermdns=Msig*(0.0)
    kappaetoile=(pv-pl)*(1/sigma)
    kappatermdns._npa[2,:]=I._npa[0,:]*Iv._npa[0,:]*sigma*(kappaetoile.derivate())._npa[0,:] # le term en grad(kappa)

    
    #kappamoy=CalculerKappaEtoile(Run)
    kappamoy=4./db ### le fit est k=
    kappaprime=kappatoile[1]
    kappaterm1=I*Iv.grad()*sigma*kappatoile[0]
    kappaterm2=I*Iv.grad()*sigma*kappamoy*(-1.0)
    kappaterm3=I*Iv.grad()*sigma*kappaprime*(-1.0)
    kappaterm1.settex(r'$kappaetoile$')
    kappaterm2.settex(r'$kappamoy$')
    kappatermtot=(kappaterm1+kappaterm2)*(-1.0)
    dpkappaterm=I*Iv.grad()*sigma*(mincurv-(4./db))
    #dpkappaterm=I*Iv.grad()*sigma*(2.0/a-kappamoy)#*2.3 ### kappamoy n est pas bien determiner par un model ellipse, ca fit avec 2.3 fois l'interval d'ou le fois 2.3
    #dpkappaterm=I*Iv.grad()*sigma*(2.0/a)+I*Run['qdm_diph']['interface']

    

    
    
    
    
    
 


    ## pressure interfacial closure relation classique 


    gradI_turb=viscousv*0.
    gradI_turb_v=viscousv*0. 
    gradI_turb_l=viscousv*0. 
    rholRl=Run['variable']['Rij_l']*rho_l
    rhovRv=Run['variable']['Rij_v']*rho_v
 
    gradI_turb._npa[0,:]=(rholRl.compo(0,0)._npa[0,:]*(Iv.grad().compo(0)._npa[0,:])
    +rholRl.compo(0,1)._npa[0,:]*(Iv.grad().compo(1)._npa[0,:])
    +rholRl.compo(0,2)._npa[0,:]*(Iv.grad().compo(2)._npa[0,:])
    -rhovRv.compo(0,0)._npa[0,:]*(I.grad().compo(0)._npa[0,:])
    -rhovRv.compo(0,1)._npa[0,:]*(I.grad().compo(1)._npa[0,:])
    -rhovRv.compo(0,2)._npa[0,:]*(I.grad().compo(2)._npa[0,:]))*Iv._npa[0,:]

    gradI_turb_v._npa[0,:]=(rhovRv.compo(0,0)._npa[0,:]*(I.grad().compo(0)._npa[0,:])*(-1.)
    -rhovRv.compo(0,1)._npa[0,:]*(I.grad().compo(1)._npa[0,:])
    -rhovRv.compo(0,2)._npa[0,:]*(I.grad().compo(2)._npa[0,:]))*Iv._npa[0,:]

    gradI_turb_l._npa[0,:]=(rholRl.compo(0,0)._npa[0,:]*(Iv.grad().compo(0)._npa[0,:])
    +rholRl.compo(0,1)._npa[0,:]*(Iv.grad().compo(1)._npa[0,:])
    +rholRl.compo(0,2)._npa[0,:]*(Iv.grad().compo(2)._npa[0,:]))*Iv._npa[0,:]



    ###



    gradI_turb._npa[1,:]=(rholRl.compo(1,0)._npa[0,:]*(Iv.grad().compo(0)._npa[0,:])
    +rholRl.compo(1,1)._npa[0,:]*(Iv.grad().compo(1)._npa[0,:])
    +rholRl.compo(1,2)._npa[0,:]*(Iv.grad().compo(2)._npa[0,:])
    -rhovRv.compo(1,0)._npa[0,:]*(I.grad().compo(0)._npa[0,:])
    -rhovRv.compo(1,1)._npa[0,:]*(I.grad().compo(1)._npa[0,:])
    -rhovRv.compo(1,2)._npa[0,:]*(I.grad().compo(2)._npa[0,:]))*Iv._npa[0,:]

    gradI_turb_v._npa[1,:]=(rhovRv.compo(1,0)._npa[0,:]*(I.grad().compo(0)._npa[0,:])*(-1.)
    -rhovRv.compo(1,1)._npa[0,:]*(I.grad().compo(1)._npa[0,:])
    -rhovRv.compo(1,2)._npa[0,:]*(I.grad().compo(2)._npa[0,:]))*Iv._npa[0,:]

    gradI_turb_l._npa[1,:]=(rholRl.compo(1,0)._npa[0,:]*(Iv.grad().compo(0)._npa[0,:])
    +rholRl.compo(1,1)._npa[0,:]*(Iv.grad().compo(1)._npa[0,:])
    +rholRl.compo(1,2)._npa[0,:]*(Iv.grad().compo(2)._npa[0,:]))*Iv._npa[0,:]

    ###

    gradI_turb._npa[2,:]=(rholRl.compo(2,0)._npa[0,:]*(Iv.grad().compo(0)._npa[0,:])
    +rholRl.compo(2,1)._npa[0,:]*(Iv.grad().compo(1)._npa[0,:])
    +rholRl.compo(2,2)._npa[0,:]*(Iv.grad().compo(2)._npa[0,:])
    -rhovRv.compo(2,0)._npa[0,:]*(I.grad().compo(0)._npa[0,:])
    -rhovRv.compo(2,1)._npa[0,:]*(I.grad().compo(1)._npa[0,:])
    -rhovRv.compo(2,2)._npa[0,:]*(I.grad().compo(2)._npa[0,:]))*Iv._npa[0,:]

    gradI_turb_v._npa[2,:]=(rhovRv.compo(2,0)._npa[0,:]*(I.grad().compo(0)._npa[0,:])*(-1.)
    -rhovRv.compo(2,1)._npa[0,:]*(I.grad().compo(1)._npa[0,:])
    -rhovRv.compo(2,2)._npa[0,:]*(I.grad().compo(2)._npa[0,:]))*Iv._npa[0,:]

    gradI_turb_l._npa[2,:]=(rholRl.compo(2,0)._npa[0,:]*(Iv.grad().compo(0)._npa[0,:])
    +rholRl.compo(2,1)._npa[0,:]*(Iv.grad().compo(1)._npa[0,:])
    +rholRl.compo(2,2)._npa[0,:]*(Iv.grad().compo(2)._npa[0,:]))*Iv._npa[0,:]
   

    reste_turb_v=turb_v*I*(-1.)-gradI_turb_v
    reste_turb_l=turb_l*Iv-gradI_turb_l


    gradI_turb=gradI_turb*(-1.)
    viscous=viscousv*I-viscousl*Iv
    turb=turb_l*Iv-turb_v*I
    reste_turb=turb-gradI_turb
    turb_lift=((turb_v+turb_l)*Iv*(-1.)-gradI_turb*(-1.))*(-1.)
    turb_non_force=reste_turb-turb_lift

    #flotta=g*I*Iv*(rhol-rhov)*(-1.)
    flotta=turb*0.
    flotta._npa[0,:]=I._npa[0,:]*Iv._npa[0,:]*(rhol-rhov)*gg*(-1.)
    #flotta=flotta+Run['qdm_gaz']['source']*(Iv-I)
    pression=((pv-pl)*Iv).grad()*I*(-1.)#+pl*(Iv.grad())
    interface=Msig*I*(-1.)
    Ml=(pression+interface+flotta+turb_non_force+turb_lift+gradI_turb+viscous)
    tracerJFM([viscousl.compo(2)*(-1.), turb_l.compo(2), ((pl*I).grad()).compo(2), Ml.compo(2), viscousl.compo(2)+turb_l.compo(2)+((pl*I).grad()).compo(2)], 'JFMequilibreDNSz_qdml' )

    tracerJFM([Run['qdm_gaz']['viscous'].compo(0), Run['qdm_gaz']['turbulence'].compo(0), Run['qdm_gaz']['pression'].compo(0), Run['qdm_gaz']['interface'].compo(0), Run['qdm_gaz']['source'].compo(0), Run['qdm_gaz']['rhog'].compo(0)], 'JFMequilibreDNSz_qdmv' )

    dragref=Iv*((uv-ul).compo(0))*((uv-ul).compo(0))
    scale=max(abs((flotta).compo(0).symetriser([1.0]) ._npa[0,:])) 
    scale_ul=max(abs(dragref._npa[0,:]))


    Ml.symetriser([1,1,-1])
    Ml_l_DNS=Run['qdm_liq']['interface']*(-1.)+pl*(Iv.grad())
    #tracer([ Ml, pression, interface, flotta, turb], Run['variable']['name']+'equilibreDNSx', [0], couleur=[1,1,1,1,1,1,1], xlim=[0,1])
    #tracer([ Ml, pression, interface, flotta, turb] , Run['variable']['name']+'equilibreDNSz', [2], couleur=[3,1,1,1,1,1], xlim=[0,1], ylim=[-0.005, 0.005], sym=[-1.0])
    #tracer([ Ml, pression, interface, flotta, turb] , Run['variable']['name']+'equilibreDNSz', [2], couleur=[1,1,1,1,1,1], xlim=[0,1], sym=[-1.0])
    tracerJFM([ Ml.compo(0).symetriser([1.0]), interface.compo(0).symetriser([1.0])+pression.compo(0).symetriser([1.0]), flotta.compo(0).symetriser([1.0])+turb_lift.compo(0).symetriser([1.0]), gradI_turb.compo(0).symetriser([1.0]), viscous.compo(0).symetriser([1.0]),turb_non_force.compo(0).symetriser([1.0]), dragref*(scale/scale_ul)], "JFMequilibreDNSx")
    tracerJFM([ Ml.compo(2).symetriser([-1.0]), interface.compo(2).symetriser([-1.0])+pression.compo(2).symetriser([-1.0]), flotta.compo(2).symetriser([-1.0])+turb_lift.compo(2).symetriser([-1.0]), gradI_turb.compo(2).symetriser([-1.0]), viscous.compo(2).symetriser([-1.0]),turb_non_force.compo(2).symetriser([-1.0]), (((pl.grad())*I).compo(2)*(-1.)).symetriser([-1.0])], "JFMequilibreDNSz")





    gradI=(I.grad()).compo(0)
    flotta=turb*0.
    flotta._npa[0,:]=I._npa[0,:]*Iv._npa[0,:]*(rhol-rhov)*gg
    pression_v_1=((pv-pl)*Iv).grad()*Iv*(-1.)
    interface_v=interface*(Iv*I.inv())
    M_autre=(pl.grad())*Iv*(-1.)
    M_autre2=turb_v

    #beta=Run['qdm_diph']['source']*Iv

    turb_lift=turb_lift*(-1.)
    gradI_turb=gradI_turb*(-1.)
    turb_non_force=turb_non_force*(-1.)
    Ml_v_sum=(pression_v_1+interface_v+flotta+turb_lift+gradI_turb+viscous+turb_non_force)
    Ml_v_DNS=Run['qdm_gaz']['interface']-((pv-pl)*Iv).grad()-pl*(Iv.grad())
 

    gradux=(ul.grad()).compo(0,2).symetriser([-1.0])
    liftref=(uv-ul).compo(0)*Iv*gradux
    scale=max(abs((reste_turb_l).compo(2).symetriser([-1.0]) ._npa[0,:])) 
    scale_ul=max(abs(liftref._npa[0,:]))
    dragref=Iv*((uv-ul).compo(0))*((uv-ul).compo(0))
    scale=max(abs((flotta).compo(0).symetriser([1.0]) ._npa[0,:])) 
    scale_ul=max(abs(dragref._npa[0,:]))


    tracerJFM([Ml_v_sum.compo(2).symetriser([-1.0]),(pression_v_1.compo(2).symetriser([-1.0])+interface_v.compo(2).symetriser([-1.0])), flotta.compo(2).symetriser([-1.0])+turb_lift.compo(2).symetriser([-1.0]), gradI_turb.compo(2).symetriser([-1.0]), viscous.compo(2).symetriser([-1.0]), turb_non_force.compo(2).symetriser([-1.0])], "JFMequilibreDNSz_vapeur")

    tracerJFM([Ml_v_sum.compo(0).symetriser([1.0]), (pression_v_1.compo(0).symetriser([1.0])+interface_v.compo(0).symetriser([1.0])), flotta.compo(0).symetriser([1.0])+turb_lift.compo(0).symetriser([-1.0]), gradI_turb.compo(0).symetriser([1.0]), viscous.compo(0).symetriser([1.0]), turb_non_force.compo(0).symetriser([1.0]), dragref*(scale/scale_ul)*(-1.)], "JFMequilibreDNSx_vapeur")



   ####### equilibre des forces #############
    alpha=Run['variable']['alpha']
    rho_m=Iv*rho_v+I*rho_l
    rho_av=alpha*rho_v+(1.-alpha)*rho_l
    buoyancy=Iv*g*(rho_v-rho_av)
    drag_part=Iv*I*g*(rho_l-rho_v)

    #buoyancy=g*0.
    #buoyancy._npa[0,:]=Iv._npa[0,:]*Run['qdm_gaz']['source']

    #buoyancy=(Iv*rho_v+I*rho_l)*g*Iv*Iv*(-1.)
    viscous=(Run['qdm_liq']['viscous']+Run['qdm_gaz']['viscous'])*Iv
    reste_p=((pv-pl)).grad()*Iv*Iv*(-1.)+(pl.grad())*Iv*(-1.)+Run['qdm_gaz']['source']
    gradI_p=(pv-pl)*(Iv.grad())*Iv*(-1.)+Run['qdm_diph']['interface']*Iv
    gradI_turb=gradI_p*0.
    rholRl=Run['variable']['Rij_l']*rho_l
    rhovRv=Run['variable']['Rij_v']*rho_v
 
    gradI_turb._npa[0,:]=(rholRl.compo(0,0)._npa[0,:]*(Iv.grad().compo(0)._npa[0,:])
    +rholRl.compo(0,1)._npa[0,:]*(Iv.grad().compo(1)._npa[0,:])
    +rholRl.compo(0,2)._npa[0,:]*(Iv.grad().compo(2)._npa[0,:])
    -rhovRv.compo(0,0)._npa[0,:]*(Iv.grad().compo(0)._npa[0,:])
    -rhovRv.compo(0,1)._npa[0,:]*(Iv.grad().compo(1)._npa[0,:])
    -rhovRv.compo(0,2)._npa[0,:]*(Iv.grad().compo(2)._npa[0,:]))*Iv._npa[0,:]

    gradI_turb._npa[1,:]=(rholRl.compo(1,0)._npa[0,:]*(Iv.grad().compo(0)._npa[0,:])
    +rholRl.compo(1,1)._npa[0,:]*(Iv.grad().compo(1)._npa[0,:])
    +rholRl.compo(1,2)._npa[0,:]*(Iv.grad().compo(2)._npa[0,:])
    -rhovRv.compo(1,0)._npa[0,:]*(Iv.grad().compo(0)._npa[0,:])
    -rhovRv.compo(1,1)._npa[0,:]*(Iv.grad().compo(1)._npa[0,:])
    -rhovRv.compo(1,2)._npa[0,:]*(Iv.grad().compo(2)._npa[0,:]))*Iv._npa[0,:]

    gradI_turb._npa[2,:]=(rholRl.compo(2,0)._npa[0,:]*(Iv.grad().compo(0)._npa[0,:])
    +rholRl.compo(2,1)._npa[0,:]*(Iv.grad().compo(1)._npa[0,:])
    +rholRl.compo(2,2)._npa[0,:]*(Iv.grad().compo(2)._npa[0,:])
    -rhovRv.compo(2,0)._npa[0,:]*(Iv.grad().compo(0)._npa[0,:])
    -rhovRv.compo(2,1)._npa[0,:]*(Iv.grad().compo(1)._npa[0,:])
    -rhovRv.compo(2,2)._npa[0,:]*(Iv.grad().compo(2)._npa[0,:]))*Iv._npa[0,:]

    reste_turb=(turb_v+turb_l)*Iv*(-1.)-gradI_turb+drag_part

    residu=reste_p+gradI_p+gradI_turb+reste_turb+buoyancy+viscous

    a=gradI_p.compo(0).symetriser([1.0])
    b=gradI_turb.compo(0).symetriser([1.0])
    c=(reste_p).compo(0).symetriser([1.0])
    d=(reste_turb).compo(0).symetriser([1.0])
    e=(viscous).compo(0).symetriser([1.0])
    f=(buoyancy).compo(0).symetriser([1.0])
    g=residu.compo(0).symetriser([1.0])

    rb=Run['variable']['db']*0.5
    dragref=Iv*(uv-ul).compo(0)*(uv-ul).compo(0)*rho_l*(1./8.)*(3./rb)
    scale=max(abs((reste_turb).compo(0).symetriser([1.0]) ._npa[0,:])) 
    scale_ul=max(abs(dragref._npa[0,:]))

    print "Cd", (scale/scale_ul)

    tracerJFM([a,b,e,f,g,c,d,dragref*(scale/scale_ul)*(-1.)], "JFMequilibre_force_x_vapeur")
    #tracerJFM([a,b,c,d,e,f,g], "JFMequilibre_force_x_vapeur")

    a=gradI_p.compo(2).symetriser([-1.0])
    b=gradI_turb.compo(2).symetriser([-1.0])
    c=(reste_p).compo(2).symetriser([-1.0])
    d=(reste_turb).compo(2).symetriser([-1.0])
    e=(viscous).compo(2).symetriser([-1.0])
    f=(buoyancy).compo(2).symetriser([-1.0])
    g=residu.compo(2).symetriser([-1.0])


    gradux=(ul.grad()).compo(0,2).symetriser([-1.0])
    liftref=(uv-ul).compo(0)*Iv*gradux*rho_l
    scale=max(abs((reste_turb).compo(2).symetriser([-1.0]) ._npa[0,:])) 
    scale_P=max(abs((reste_p).compo(2).symetriser([-1.0]) ._npa[0,:])) 
    scale_ul=max(abs(liftref._npa[0,:]))

    print "CL*=", scale/scale_ul
    tracerJFM([a,b,e,f,g,c,d, liftref*(scale/scale_ul), liftref*(scale/scale_ul)*(-1.)], "JFMequilibre_force_z_vapeur")
    tracerJFM([a,b,c,d,g, liftref*(scale/scale_ul), liftref*(scale/scale_ul)*(-1.)], "JFMequilibre_force_z_vapeur_GBois")
    def moyenne_sur_bulle(field, Iv):
	moy=0.
	nb=0.
	for i in range(len(field._ax[0,:])):
		if (Iv._npa[0,i]!=0.):
			nb=nb+1.
			moy=moy+field._npa[0,i]
	return moy/nb
    ## Ivrappelz comparable a g
    #try: 


    alpha=Run['variable']['alpha']
    viscous=(Run['qdm_liq']['viscous']+Run['qdm_gaz']['viscous'])
    reste_p=((pv-pl)).grad()*Iv*(-1.)+(pl.grad())*(-1.)
    gradI_p=(pv-pl)*(Iv.grad())*(-1.)+Run['qdm_diph']['interface']
    gradI_turb=gradI_p*0.
    rholRl=Run['variable']['Rij_l']*rho_l
    rhovRv=Run['variable']['Rij_v']*rho_v
 
    gradI_turb._npa[0,:]=(rholRl.compo(0,0)._npa[0,:]*(Iv.grad().compo(0)._npa[0,:])
    +rholRl.compo(0,1)._npa[0,:]*(Iv.grad().compo(1)._npa[0,:])
    +rholRl.compo(0,2)._npa[0,:]*(Iv.grad().compo(2)._npa[0,:])
    -rhovRv.compo(0,0)._npa[0,:]*(Iv.grad().compo(0)._npa[0,:])
    -rhovRv.compo(0,1)._npa[0,:]*(Iv.grad().compo(1)._npa[0,:])
    -rhovRv.compo(0,2)._npa[0,:]*(Iv.grad().compo(2)._npa[0,:]))

    gradI_turb._npa[1,:]=(rholRl.compo(1,0)._npa[0,:]*(Iv.grad().compo(0)._npa[0,:])
    +rholRl.compo(1,1)._npa[0,:]*(Iv.grad().compo(1)._npa[0,:])
    +rholRl.compo(1,2)._npa[0,:]*(Iv.grad().compo(2)._npa[0,:])
    -rhovRv.compo(1,0)._npa[0,:]*(Iv.grad().compo(0)._npa[0,:])
    -rhovRv.compo(1,1)._npa[0,:]*(Iv.grad().compo(1)._npa[0,:])
    -rhovRv.compo(1,2)._npa[0,:]*(Iv.grad().compo(2)._npa[0,:]))

    gradI_turb._npa[2,:]=(rholRl.compo(2,0)._npa[0,:]*(Iv.grad().compo(0)._npa[0,:])
    +rholRl.compo(2,1)._npa[0,:]*(Iv.grad().compo(1)._npa[0,:])
    +rholRl.compo(2,2)._npa[0,:]*(Iv.grad().compo(2)._npa[0,:])
    -rhovRv.compo(2,0)._npa[0,:]*(Iv.grad().compo(0)._npa[0,:])
    -rhovRv.compo(2,1)._npa[0,:]*(Iv.grad().compo(1)._npa[0,:])
    -rhovRv.compo(2,2)._npa[0,:]*(Iv.grad().compo(2)._npa[0,:]))

    reste_turb=(turb_v+turb_l)*(-1.)-gradI_turb#+drag_part

    try:
	    a=(Run['qdm_diph']['interface']).compo(2)
	    b=((pv-pl)*(Iv.grad())*(-1.)).compo(2)
	    c=(((pv-pl)).grad()*Iv*(-1.)).compo(2)
	    d=gradI_turb.compo(2)
	    e=((pl.grad())*(-1.)).compo(2)
	    f=reste_turb.compo(2)
	    g=viscous.compo(2)
	    rappel=Run['variable']['Ivrappelz']*rho_v
	    residu=a+b+c+d+e+f+g+rappel

	    tracerJFM([a*Iv+b*Iv,c*Iv,d*Iv,e*Iv,f*Iv,residu*Iv,rappel*Iv], "JFMequilibre_force_z_vapeur_2")
	    print "surface tension source", moyenne_sur_bulle(a, Iv)
	    print "interfacial pressure", moyenne_sur_bulle(b, Iv)
	    print "laminar dispersion", moyenne_sur_bulle(a+b, Iv)
	    print "reste pression", moyenne_sur_bulle(c, Iv)
	    print "turbulent disepersion", moyenne_sur_bulle(d, Iv)
	    print "gradP", moyenne_sur_bulle(e, Iv)
	    print "Lift", moyenne_sur_bulle(f, Iv)
	    print "viscous", moyenne_sur_bulle(g, Iv)
	    print "rappel", moyenne_sur_bulle(rappel, Iv)
	    print "residu", moyenne_sur_bulle(residu, Iv)
	    #except:
		    #print "bulle libres"


	    print "Cl_Re theorique=", 6./20.
	    print "Cl_Re=", (scale/scale_ul)
	    beta=9./6.4
	    We=0.15*0.15*rho_l*db/sigma
	    print "We=", We
	    print "Cl_P theorique=", (32*beta)/(1+8*beta*We)
	    print "Cl_P=", (scale_P/scale_ul)

	    tracerJFM([(reste_turb.compo(2)*(liftref).inv()*(-1.))], "JFMequilibre_force_z_lift")
	    #tracerJFM([reste_turb.compo(2), (liftref)*(scale/scale_ul)], "JFMequilibre_force_z_liiift")

    except:
	    print "pas de force de rappel"

    #tracer([flotta, Ml, turb, dpkappaterm, kappatermdns], Run['variable']['name']+'equilibreRANSxNew', [0], couleur=[1,1,1,1,1,1,3,1], xlim=[0,1])
    #tracer([flotta, Ml, turb, kappaterm1, kappaterm2+kappaterm3], Run['variable']['name']+'equilibreRANSzNew', [2], couleur=[1,1,1,1,1,1,1], xlim=[0,1], ylim=[-0.005, 0.005], sym=[-1.0])



    interface2=interface*I#*(-1.0)
    interface2.settex(r"$\alpha_l\sigma\overline{\kappa\nabla\chi_k}$")
    try :
    	msig1=I*kdns*sigma*Iv.grad()
    except :
	msig1=Iv.grad()*0.
    msig1.settex(r'$\alpha_l\sigma\overline{\kappa}^i\nabla\alpha_v$')
    msig2=interface-msig1
    msig2.settex(r"$\alpha_l\sigma\left(\overline{\kappa\nabla\chi_k}-\overline{\kappa}^i\nabla\alpha_v\right)$")
    kappaterm2.settex(r"$\alpha_l\sigma\left(2/r_{0}\right)\nabla\alpha_v$")
    kappaterm3.settex(r"$\alpha_l\sigma\left(1/r_0-r_0^3/a^4\right)\nabla\alpha_v$")
    kappatermfinal=kappaterm2+kappaterm3
    kappatermfinal.settex(r"$\alpha_l\sigma\left(3/r_0-r_0^3/a^4\right)\nabla\alpha_v$")
    #tracer([msig1, kappaterm2, msig2, kappaterm3, interface, kappatermfinal], Run['variable']['name']+'msigcirrelation', [2], couleur=[2,2,2,1,1,1], xlim=[0,1], ylim=[-0.005, 0.0005], sym=[-1.0])
    #tracer([msig1, kappaterm2, msig2, kappaterm3, interface2, kappatermfinal], Run['variable']['name']+'msigcirrelation', [2], couleur=[2,2,2,1,1,1], xlim=[0,1], sym=[-1.0])
    #tracer([flotta, Ml, turb, residu3], Run['variable']['name']+'equilibreRANSxOld', [0], couleur=[1,1,1,1,1,1,3,1], xlim=[0,1])
    #tracer([flotta, Ml, turb, residu3], Run['variable']['name']+'equilibreRANSzOld', [2], couleur=[1,1,1,1,1,1,3,1], xlim=[0,1], ylim=[-0.005, 0.005], sym=[-1.0])
    #tracer([flotta, Ml, turb, pression, pression3, interface, dpkappaterm, residu, residu3, residu4], Run['variable']['name']+'AmouraequilibreDNSx', [0], couleur=[1,1,1,2,1,1,3,1])
    #tracer([flotta, Ml, turb, pression, pression3, interface, dpkappaterm, residu, residu3, residu4], Run['variable']['name']+'AmouraequilibreDNSz', [2], couleur=[1,1,1,2,1,1,3,1])
    #mstd=Ml+pl*Iv.grad()
    #mstd.settex(r'$-\mathbf{M_l^{Std}}$')
    #residu=flotta+mstd+turb
    #residu.settex(r'residu')
    #tracer([flotta, mstd, turb, residu], Run['variable']['name']+'equilibreRANSapriorix', [0], couleur=[1,1,1,1], sym=[1], xlim=[0,1])
    #tracer([flotta, mstd, turb, residu], Run['variable']['name']+'equilibreRANSaprioriz', [2], couleur=[1,1,1,1], sym=[-1], xlim=[0,1])


    #kappatoile=calculerPression(Run, inter=Ml, pltilde=0.0288)
    kappatoile=calculerPression(Run, inter=Ml, pltilde=0.0)

    plestim=Run['variable']['plestimC075']
    pvestim=Run['variable']['pvestimC075']
    #tracer([plestim, pvestim], 'pestim', xlim=[0.1,1.])
    #tracer([I, Iv], 'Ivvv')
    #pression_estim=(((pvestim-plestim)*Iv).grad())*I*(-1.0)+plestim*(Iv.grad())
    pression_estim=I*(pvestim-plestim)*Iv.grad()
    interface_estim=kappatermfinal
    Mlestim=(interface_estim+pression_estim)*(-1.0)
    pression.settex(r"$\mathbf{M_l^P}$ DNS")
    pression_estim.settex(r"$\mathbf{M_l^P}$ model") 
    interface.settex(r"$\mathbf{M_l^\sigma}$ DNS")
    interface_estim.settex(r"$\mathbf{M_l^\sigma}$ model") 
    turb.settex(r"$\mathbf{M_l^{TD}}$ DNS")
    Ml.settex(r"$\mathbf{M_l^{disp}}$ DNS")
    Mlestim.settex(r"$\mathbf{M_l^{disp}}$ model")
    #tracer([pression, pression_estim, interface, interface_estim, turb, Ml, Mlestim] , Run['variable']['name']+'equilibreDNSzapriori', [2], couleur=[2,2,1,2,2,2], xlim=[0,1], ylim=[-0.2, 0.2], sym=[-1.0], markevry=[0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05])

    tracerJFM([pression.compo(2), interface.compo(2), Ml.compo(2), pression_estim.compo(2)*(-1.), interface_estim.compo(2)*(-1.), Mlestim.compo(2)], 'JFMequilibreDNSzapriori')



    return

def CalculerKappaEtoile(Run):

    #### construction d'un gradient de curvature 
    db=Run['variable']['db']
    gg=Run['variable']['g']    
    rho_l=Run['variable']['rho_l'] 
    rho_v=Run['variable']['rho_v']
    sigma=Run['variable']['sigma'] 
    II=Run['variable']['I'] 
    rb=db/2.0
    ### Pour le demi grand axe, on peut utiliser la correlation de wellek (attention en réalité il ya le pitching etc...)
    wellek=True
    if wellek:
        Eo=(abs(gg)*rho_l*db*db)/sigma
        rbb=rb*((1+0.163*Eo**(0.757))**(0.33333))
        a=rbb
        print 'Eo=', Eo, 'demi grand axe=', a        
    else:
        a=0.18 #demi grand axe visible (à la main)
    b=(rb*rb*rb)/(a*a)
    I=II*0.0
    L=I._ax[0,-1]
    pi=3.1415
    for i in range(len(I._ax[0,:])):
	I._ax[0,i]=I._ax[0,i]*pi/L
    cost=I*0.0    
    sint=I*0.0   
    for i in range(len(I._ax[0,:])):
        cost._npa[0,i]=math.cos(I._ax[0,i])
        sint._npa[0,i]=math.sin(I._ax[0,i])
	I._ax[0,i]=I._ax[0,i]*L/pi      
            
    curvature=((cost*cost*b*b+sint*sint*a*a).power(1.5)).inv(0.0000001)*a*b
    integcurv=curvature.integ()*(1.0/pi)
    dd=sint*2.0*b
    dkappa=dd*curvature
    curvature.settex(r'curvature')
    dkappa.settex(r'dkappa')
    dd.settex(r'd')
    kappa=dkappa.integ()*(1/(4.0*b))
    kappa0=kappa*0.0+2.0/db
    kappa1=kappa*0.0+(4.0*b/a)/(4.0*b)
    curvature.settex(r'curvature')
    dkappa.settex(r'dkappa')
    integcurv.settex(r'curvature moyenne')
    dd.settex(r'd')
    dkappa.settex(r'curvature ponderree')
    kappa.settex(r'integ curv ponderree')
    kappa0.settex(r'a/b*b')
    kappa1.settex(r'1/a')
    tracer([curvature, integcurv, dd, dkappa, kappa, kappa0, kappa1], 'kappaetoile')
    kappa_moy=2.0*integcurv._npa[0,-1] ## 2 fois car on est en 3D


    return kappa_moy

def CalculerKappaNprime(Run):

    #### construction d'un gradient de curvature 
    db=Run['variable']['db']
    gg=Run['variable']['g']    
    rho_l=Run['variable']['rho_l'] 
    rho_v=Run['variable']['rho_v']
    sigma=Run['variable']['sigma'] 
    II=Run['variable']['I'] 
    rb=db/2.0
    ### Pour le demi grand axe, on peut utiliser la correlation de wellek (attention en réalité il ya le pitching etc...)
    wellek=True
    if wellek:
        Eo=(abs(gg)*rho_l*db*db)/sigma
        rbb=rb*((1+0.163*Eo**(0.757))**(0.33333))
        a=rbb
        print 'Eo=', Eo, 'demi grand axe=', a        
    else:
        a=0.18 #demi grand axe visible (à la main)
    b=(rb*rb*rb)/(a*a)
    I=II*0.0
    L=I._ax[0,-1]
    pi=3.1415*2.0
    for i in range(len(I._ax[0,:])):
	I._ax[0,i]=I._ax[0,i]*pi/L
    cost=I*0.0    
    sint=I*0.0   
    for i in range(len(I._ax[0,:])):
        cost._npa[0,i]=math.cos(I._ax[0,i])
        sint._npa[0,i]=math.sin(I._ax[0,i])
	I._ax[0,i]=I._ax[0,i]*L/pi      
            
    curvature=((cost*cost*b*b+sint*sint*a*a).power(1.5)).inv(0.0000001)*a*b
    nz=((cost*cost*b*b+sint*sint*a*a).power(0.5)).inv(0.0000001)*b*cost*(-1.0)
    curvnz=curvature*nz
    integcurvnz=curvnz.integ()*(1.0/pi)
    curvature.settex(r'$\kappa$')
    nz.settex(r'$n_z$')
    curvnz.settex(r'$\kappa n_z$')
    integcurvnz .settex(r'$\overline{\kappa n_z}$')
    tracer([curvature, nz, curvnz, integcurvnz], 'kappanzprime')
    kappa_moy=2.0*integcurvnz._npa[0,-1] ## 2 fois car on est en 3D


    return kappa_moy

def EquationVrNeptune(fstat, Run):
    sigma=Run['variable']['sigma']
    db=Run['variable']['db']
    rhol=Run['variable']['rho_l']
    rhov=Run['variable']['rho_v'] 
    Iv=Field.LoadFromFile(fstat,["y"] , ["I"],'dP', 'z', r'$\nabla P$') 
    I=(Iv-1.0)*(-1.0)     
    g=Field.initgravity([Run['variable']['g'],0,0], 'g', I)
    tab=["", "v"]
    tab2=["liq", "gaz"]
    P=Field.LoadFromFile(fstat,["y"] , ["P"],'P', 'z', r'$\nabla P$')
    dP=Field.LoadFromFile(fstat,["y"] , ["dPdx","dPdz","dPdy"],'dP', 'z', r'$\nabla P$')
    conv=Field.LoadFromFile(fstat,["y"] , ["convx"+tab[0],"convz"+tab[0],"convy"+tab[0]],'dP', 'z', r'$\nabla P$')
    drij=Field.LoadFromFile(fstat,["y"] , ["drijx"+tab[0],"drijz"+tab[0],"drijy"+tab[0]],'dP', 'z', r'$\nabla P$')
    drijv=Field.LoadFromFile(fstat,["y"] , ["drijx"+tab[1],"drijz"+tab[1],"drijy"+tab[1]],'dP', 'z', r'$\nabla P$')
    dtau=Field.LoadFromFile(fstat,["y"] , ["dtaux"+tab[0],"dtauz"+tab[0],"dtauy"+tab[0]],'dP', 'z', r'$\nabla P$')
    dtauv=Field.LoadFromFile(fstat,["y"] , ["dtaux"+tab[1],"dtauz"+tab[1],"dtauy"+tab[1]],'dP', 'z', r'$\nabla P$')
    fd=Field.LoadFromFile(fstat,["y"] , ["fdx","fdz","fdy"],'dP', 'z', r'$\nabla P$')
    fl=Field.LoadFromFile(fstat,["y"] , ["flx","flz","fly"],'dP', 'z', r'$\nabla P$')
    fdt=Field.LoadFromFile(fstat,["y"] , ["fdtx","fdtz","fdty"],'dP', 'z', r'$\nabla P$')
    
    
    
    Ml=fl+fd+fdt
    #Ml=fdt
    flotta=g*I*Iv*(rhol-rhov)
    turb=drijv*I-drij*Iv
    viscous=dtauv*I-dtau*Iv   
    residu=flotta-turb+Ml+viscous
#     sigkaDIvI=(Iv*I).grad()*sigma*(4./db)
#     IvpvdI=Iv*pv*I.grad()
#     IpldI=I*pl*I.grad()
#     pression=sigkaDIvI-IvpvdI-IpldI
    residu.settex(r'residu')
    flotta.settex(r'$\alpha_v\alpha_l(\rho_l-\rho_v)g_i$')
    turb.settex(r'$\nabla\mathbf{R_{ij}^{pv-l}}$')
    viscous.settex(r'$\nabla\mathbf{\tau_{ij}^{pv-l}}$')    
    Ml.settex(r'$\mathbf{M_l^{Std}}$')

    

    tracer([flotta, Ml, turb, viscous, residu], Run['variable']['name']+'equilibreRANSx', [0])
    tracer([flotta, Ml, turb, viscous, residu], Run['variable']['name']+'equilibreRANSy', [1])
    tracer([flotta, Ml, turb, viscous, residu], Run['variable']['name']+'equilibreRANSz', [2])


def Compa2A2(Run1,Run2,name1='rk3',name2='euler'):
    
    sym_tens=[1,1,-1,1,1,-1,-1,-1,1]
    for cle in Run1.keys():
        for cle2 in Run1[cle].keys():
            if isinstance(Run1[cle][cle2], ScalarField):
                Run1[cle][cle2].symetriser()
                try:
                    Run2[cle][cle2].symetriser()
                except:    
                    print 'cette clef n existe pas dans run2'
            elif isinstance(Run1[cle][cle2], VectorField):
                Run1[cle][cle2].symetriser()   
                try:
                    Run2[cle][cle2].symetriser() 
                except:    
                    print 'cette clef n existe pas dans run2'                    
            elif isinstance(Run1[cle][cle2], Tensor2Field):
                Run1[cle][cle2].symetriser(sym_tens)
                try: 
                    Run2[cle][cle2].symetriser(sym_tens)  
                except:    
                    print 'cette clef n existe pas dans run2'                
                
    for cle in Run1.keys():
        for cle2 in Run1[cle].keys():
            try :
                Run1[cle][cle2].settex(name1)
                Run2[cle][cle2].settex(name2) 
                if isinstance(Run1[cle][cle2], ScalarField):
                    tracer([Run1[cle][cle2],Run2[cle][cle2]], cle+'_'+cle2+'_compa', couleur=[1,1], plotall=False)
                elif isinstance(Run1[cle][cle2], VectorField):
                    tracer([Run1[cle][cle2],Run2[cle][cle2]], cle+'_'+cle2+'_compa_x',[0], couleur=[1,1], plotall=False)
                    tracer([Run1[cle][cle2],Run2[cle][cle2]], cle+'_'+cle2+'_compa_y',[1], couleur=[1,1], plotall=False)
                    tracer([Run1[cle][cle2],Run2[cle][cle2]], cle+'_'+cle2+'_compa_z',[2], couleur=[1,1], plotall=False)
                elif isinstance(Run1[cle][cle2], Tensor2Field): 
                    tracer([Run1[cle][cle2],Run2[cle][cle2]], cle+'_'+cle2+'_compa_xx',[0], couleur=[1,1], plotall=False)
                    tracer([Run1[cle][cle2],Run2[cle][cle2]], cle+'_'+cle2+'_compa_xy',[1], couleur=[1,1], plotall=False)
                    tracer([Run1[cle][cle2],Run2[cle][cle2]], cle+'_'+cle2+'_compa_xz',[2], couleur=[1,1], plotall=False)
                    tracer([Run1[cle][cle2],Run2[cle][cle2]], cle+'_'+cle2+'_compa_yx',[3], couleur=[1,1], plotall=False)
                    tracer([Run1[cle][cle2],Run2[cle][cle2]], cle+'_'+cle2+'_compa_yy',[4], couleur=[1,1], plotall=False)
                    tracer([Run1[cle][cle2],Run2[cle][cle2]], cle+'_'+cle2+'_compa_yz',[5], couleur=[1,1], plotall=False)                                      
                    tracer([Run1[cle][cle2],Run2[cle][cle2]], cle+'_'+cle2+'_compa_zx',[6], couleur=[1,1], plotall=False)
                    tracer([Run1[cle][cle2],Run2[cle][cle2]], cle+'_'+cle2+'_compa_zy',[7], couleur=[1,1], plotall=False)
                    tracer([Run1[cle][cle2],Run2[cle][cle2]], cle+'_'+cle2+'_compa_zz',[8], couleur=[1,1], plotall=False)                         
            except:
                print 'impossible de tracer', cle, cle2 
    return           




def calculerReBulk(Run):
    h=1.0
    nu=Run['variable']['mu']
    rhol=Run['variable']['rho_l']
    rhov=Run['variable']['rho_v']
    Iul=Run['variable']['u_l']
    Ivuv=Run['variable']['u_v']
    I=Run['variable']['I']
    Iinv=I.inv(1e-8)
    Iv=Run['variable']['I']-1.0
    Iv=Iv*(-1.0)
    Ivinv=Iv.inv(1e-8)
    ul=Iul*Iinv
    uv=Ivuv*Ivinv
    
    debit=(ul*rhol+uv*rhov).mean_compo(0)
    reb=debit*(4*h/nu)
    vitesse=ul.mean_compo(0)
    rebb=vitesse*(4*h/nu)
    print 'Reb=', reb
    print 'Rebb=', rebb, 'pour un ul=', vitesse
    return

def calculerEo(Run):
    h=1.0
    rhol=Run['variable']['rho_l']
    db=Run['variable']['db']
    sigma=Run['variable']['sigma']
    g=Run['variable']['g']
    print 'Eo=', rhol*db*db*g/sigma
    return
    
def calculerPression(DeformableRun, inter=0, pltilde=0.0):
    I=DeformableRun['variable']['I']
    Iinv=I.inv(1e-8)
    Iv=DeformableRun['variable']['I']-1.0
    Iv=Iv*(-1.0)
    Ivinv=Iv.inv(1e-8)

    gradI=I.grad()
    gradIv=Iv.grad()
    fstat=DeformableRun['variable']['fstat']
    print fstat
    #tracer([Field.LoadFromFile(fstat,["coordonnee_K"] , ["PI"],'P_l', 'z', r'$P_l$')], 'blabla')
    p_l=Iinv*Field.LoadFromFile(fstat,["coordonnee_K"] , ["PI"],'P_l', 'z', r'$P_l$')
    #tracer([p_l], "pl")
    p_v=Ivinv*Field.LoadFromFile(fstat,["coordonnee_K"] , ["PIv"],'P_v', 'z', r'$P_v$')  
    p_v.settex(r'DNS')
    p_l.settex(r'DNS')
    sigma=DeformableRun['variable']['sigma']
    fstat=DeformableRun['variable']['fstat']
    ai=Field.LoadFromFile(fstat,["coordonnee_K"] , ["AI"],'P_v', 'z', r'$P_v$')
    #kappa=(NN.grad()).compo(2)
    
    Vmoy=(4./3.)*3.1415*((0.3/2.)**3.0)
    I.symetriser([1])
    Iv.symetriser([1])
    ai=ai.inv()
    DDI=I.derivate().derivate()
    kappa=DDI
    #kappa=((((((I).grad()).compo(2)).grad()).compo(2))*ai*(1./Vmoy)).abs()
    #dkappa=NetContribution(kappa.grad(),0.1)
    dkappa=NetContribution(kappa.derivate()*I*Iv*sigma,0.1)
    DeformableRun['variable']['kappa']=kappa


    db=DeformableRun['variable']['db']
    kappa=kappa*0.0+(4.0/db)
    kappa_sig=kappa*sigma
    #kappa_sig=p_l*0.0+sigma*(4./0.3)
    kappa_sig.settex(r'$\sigma\kappa=\sigma\kappa^{sphe}$')
    pvl=p_v-p_l
    pvl.settex(r'DNS : $\overline{p_v}^v-\overline{p_l}^l$')
    
    #### Triche 
    db=0.3


    ### loi de Wellek
    rb=db/2.0
    gg=DeformableRun['variable']['g']
    rho_l=DeformableRun['variable']['rho_l']
    rho_v=DeformableRun['variable']['rho_v']
    Eo=(abs(gg)*rho_l*db*db)/sigma
    rbb=rb*((1+0.163*Eo**(0.757))**(0.33333))
    a=rbb
    gamma=rbb/rb
    b=(rb*rb*rb)/(a*a)
    mincurv=(b/(a*a))
    maxcurv=(a/(b*b))
    meancurv=1./rb
    kappaprimetype=4./db-(mincurv+2./db)

    ketoi5=DeformableRun['qdm_diph']['interface']*(-1.0)
    ketoi5.symetriser([-1,-1,-1])
    gradIv.symetriser([-1,-1,-1])
    ketoi5*=gradIv.inv(0.0000001)
    ketoi5=ketoi5.compo(2)



    uv=DeformableRun['variable']['u_v'].compo(0)
    ul=DeformableRun['variable']['u_l'].compo(0)
    uv=uv*Ivinv
    ul=ul*I.inv()
    ur=uv-ul
    pi=3.1415
    #tracer([uv, ul, ur], 'blablabla')
    drho=rho_l-rho_v
    Cd_choix=ketoi5*0.0+(4./3.)*a*math.sqrt(gg*(-1.0)*drho/sigma)
    e=math.sqrt(a**2-b**2)/a
    Ab=2.0*pi*a**2+pi*b**2/e*math.log((1+e)/(1-e))
    Ad=pi*a**2
    #ketoi0=kappa*0.0+2.0*mincurv*sigma
    ketoi0=kappa*0.0+(mincurv+mincurv)*sigma

    ketoimax=kappa*0.0+(maxcurv+mincurv)*sigma
    C1=1.0
    C2=0.75
    C3=0.5
    C4=0.25
    C5=0.0
    model=False
    if model :
    	Md=Ivinv*Cd_choix*ur*ur*rho_l*Ad*pi*db*db*db*(1./6.)
    elif (inter==0):
	Md=Ivinv*DeformableRun['qdm_liq']['interface'].compo(0)*(-1.0)*pi*db*db*db*(1./6.)
	print "inter DNS"
    else:
	Md=Ivinv*inter.compo(0)*pi*db*db*db*(1./6.)
        print "inter autre"
    blending=uv*0.0+1.0
    for i in range(len(blending._npa[0,:])):
        if blending._ax[0,i]<db:
                A=1./(math.exp(-db)-1.0)
		blending._npa[0,i]=(1./db)*blending._ax[0,i]
		#blending._npa[0,i]=A*(math.exp(-blending._ax[0,i])-1.0)
	#blending._npa[0,i]=-math.exp(-2*blending._ax[0,i]/db)+1.0


    volume_control_coef=2.0*(2.-gamma)
    #tracer([blending], 'blending')
    ketoi6_1=Md*(1.0/Ab)*(2*C1-1.0)*(Iv*2.0-1.0)*(-1.0)+ketoi0+Iv*(ketoimax-ketoi0)
    ketoi6_2=Md*(1.0/Ab)*(2*C2-1.0)*(Iv*2.0-1.0)*(-1.0)+ketoi0+Iv*(ketoimax-ketoi0)
    ketoi6_3=Md*(1.0/Ab)*(2*C3-1.0)*(Iv*2.0-1.0)*(-1.0)+ketoi0+Iv*(ketoimax-ketoi0)
    plestimC1=Iv*(Md*(4.0/Ab)*(2*C1-1.0)+(ketoi0-ketoimax)*volume_control_coef)
    plestimC2=Iv*(Md*(4.0/Ab)*(2*C2-1.0)+(ketoi0-ketoimax)*volume_control_coef)
    plestimC3=Iv*(Md*(4.0/Ab)*(2*C3-1.0)+(ketoi0-ketoimax)*volume_control_coef)
    plestimC4=Iv*(Md*(4.0/Ab)*(2*C4-1.0)+(ketoi0-ketoimax)*volume_control_coef)
    plestimC5=Iv*(Md*(4.0/Ab)*(2*C5-1.0)+(ketoi0-ketoimax)*volume_control_coef)
    pvestimC1=Md*(2.0/Ab)*(2*C1-1.0)+ketoi0
    pvestimC2=Md*(2.0/Ab)*(2*C2-1.0)+ketoi0
    pvestimC3=Md*(2.0/Ab)*(2*C3-1.0)+ketoi0
    pvestimC4=Md*(2.0/Ab)*(2*C4-1.0)+ketoi0
    pvestimC5=Md*(2.0/Ab)*(2*C5-1.0)+ketoi0
    plestimC1.settex(r'$C=1$')
    plestimC2.settex(r'$C=0.75$')
    plestimC3.settex(r'$C=0.5$')
    plestimC4.settex(r'$C=0.25$')
    plestimC5.settex(r'$C=0.0$')
    pvestimC1.settex(r'$C=1$')
    pvestimC2.settex(r'$C=0.75$')
    pvestimC3.settex(r'$C=0.5$')
    pvestimC4.settex(r'$C=0.25$')
    pvestimC5.settex(r'$C=0.0$')


    ketoi6=ketoi6_1*1.0
    ketoi6.settex(r'$\frac{F^D}{S_b}(2C-1)+\sigma\kappa^{min}$')
    ketoi6_1.settex(r'$\sigma\overline{\kappa^*}$ with $C=1$')
    ketoi6_2.settex(r'$\sigma\overline{\kappa^*}$ with $C=0.75$')
    ketoi6_3.settex(r'$\sigma\overline{\kappa^*}$ with $C=0.5$')



    ketoi2=kappa*0.0+2.0*maxcurv*sigma
    ketoi3=kappa*0.0+2.0*(1.0/rb)*sigma
    ketoi4bis=CalculerKappaEtoile(DeformableRun)
    ketoi4=ketoi3*0.0+ketoi4bis*sigma
    ketoi0.settex(r'$\sigma\kappa^*=\sigma\kappa^{ellipse}_{min}$')
    ketoi2.settex(r'$\sigma\kappa^*=\sigma\kappa^{ellipse}_{max}$')
    ketoi3.settex(r'$2\sigma/r_{0}$')
    ketoi4.settex(r'$\sigma\kappa^*=\sigma\overline{\kappa^{ellipse}}$')
    ketoi5.settex(r'$M_\sigma/\nabla\alpha_v$')
    ketoi=ketoi0*0.0+(2.0/a)*sigma
    ketoi.settex(r'$\sigma\kappa^*=2\sigma/a$')
    #kai=(Field.LoadFromFile('Stats_kai.dt_ev',["coordonnee_K"] , ["kai"],'P_v', 'z', r'$P_v$'))
    #ai=(Field.LoadFromFile('Stats_kai.dt_ev',["coordonnee_K"] , ["AI"],'P_v', 'z', r'$P_v$'))
    #fstat="Stats_kai.dt_ev" ### only for D127

    try :
	    kai=(Field.LoadFromFile(fstat,["coordonnee_K"] , ["kai"],'P_v', 'z', r'$P_v$'))
	    ai=(Field.LoadFromFile(fstat,["coordonnee_K"] , ["AI"],'P_v', 'z', r'$P_v$'))
	    fstat=DeformableRun['variable']['fstat']
	    kdns=kai*ai.inv()*sigma*(-1.0)
	    kdns.symetriser([1])
	    kdns.settex(r'DNS')
	    #tracer([p_v, p_l, pvl, ketoi6, ketoi0, ketoi, ketoi2, ketoi3, ketoi4, kdns, ketoi5], 'LaplaceLaw1', ylim=[-0.002, 0.04], xlim=[0.0,1.0])
	    #tracer([p_v, p_l, pvl, kdns], 'LaplaceLaw2', ylim=[-0.002, 0.04], xlim=[0.0,1.0])
	    #tracer([kdns, ketoi3], 'LaplaceLaw3', couleur=[3,4,2], ylim=[0.02, 0.04], xlim=[0.0,1.0], markevry=[3])
    except :
            print "k ai non post-traite"

    #tracer([p_l, plestimC1, plestimC2, plestimC3, plestimC4, plestimC5], 'LaplaceLaw4', couleur=[1,6,4,2], xlim=[0.0,1.0], markevry=[3])
    #tracer([p_v, pvestimC1, pvestimC2, pvestimC3, pvestimC4, pvestimC5], 'LaplaceLaw5', couleur=[1,6,4,2], xlim=[0.0,1.0], markevry=[3])

    #tracerJFM([p_l, plestimC1, plestimC2, plestimC3, plestimC4, plestimC5], "C_pl")
    #tracerJFM([p_v, pvestimC1, pvestimC2, pvestimC3, pvestimC4, pvestimC5], "C_pv")

    tracerJFM([p_l+pltilde, plestimC1, plestimC2, plestimC3, plestimC4, plestimC5], "Cpl")
    tracerJFM([p_v+pltilde, pvestimC1, pvestimC2, pvestimC3, pvestimC4, pvestimC5], "Cpv")



    tracerJFM([(p_l+pltilde).grad().compo(2)*Iv*(-1.), plestimC2.grad().compo(2)*Iv*(-1.), (Iv)*(((Md*(4.0/Ab)*(2*C2-1.0)+ketoi0-ketoimax)*volume_control_coef).grad().compo(2))*Iv*(-1.), ((Md*(4.0/Ab)*(2*C2-1.0)+ketoi0-ketoimax)*volume_control_coef)*(Iv.grad().compo(2))*Iv*(-1.)], "Cpl_grad")


    qdmdp_l=DeformableRun['qdm_liq']['pression']*(-1.0)
    qdmdp_v=DeformableRun['qdm_gaz']['pression']*(-1.0)
    dpl=p_l.grad()
    dI=I.grad()
    grad_ap_l=DeformableRun['variable']['p_l'].grad()
    a_gradp_l=dpl*I
    p_l_grada=p_l*dI
    a_gradp_v=p_v.grad()*Iv
    grad_ap_v=DeformableRun['variable']['p_v'].grad()
    p_v_grada=p_v*gradIv
    som_v=a_gradp_v+p_v_grada
    som_l=a_gradp_l+p_l_grada
    a_gradp_v.settex(r'$\alpha_v\nabla p_v$')
    a_gradp_l.settex(r'$\alpha_l\nabla p_l$')
    grad_ap_v.settex(r'$\nabla \alpha_v p_v$')
    grad_ap_l.settex(r'$\nabla \alpha_l p_l$')
    p_l_grada.settex(r'$p_l\nabla\alpha_l $')
    p_v_grada.settex(r'$p_v\nabla\alpha_v $')
    som_v.settex(r'$p_v\nabla\alpha_v +\alpha_v\nabla p_v $')
    som_l.settex(r'$p_l\nabla\alpha_l +\alpha_l\nabla p_l $')
    tab=[a_gradp_v, grad_ap_v, p_v_grada, som_v, a_gradp_l, grad_ap_l, p_l_grada, som_l, qdmdp_l, qdmdp_v]
    bb=[4,4,2]
    #tracer(tab, 'dpx', [0], xlim=[0.2,1.8], couleur=bb, markevry=[0.05])
    #tracer(tab, 'dpy', [1], xlim=[0.2,1.8], couleur=bb, markevry=[0.05])
    #tracer(tab, 'dpz', [2], xlim=[0.2,1.8], couleur=bb, markevry=[0.05])
    DeformableRun['variable']['av_dpv']=a_gradp_v
    DeformableRun['variable']['pv_dav']=p_v_grada
    DeformableRun['variable']['al_dpl']=a_gradp_l
    DeformableRun['variable']['pl_dal']=p_l_grada
    DeformableRun['variable']['plestimC075']=plestimC2
    DeformableRun['variable']['pvestimC075']=pvestimC2

    sortie=[ketoi6*(1./sigma), kappaprimetype]
    return sortie


def PositionBulles(Run):
    pi=3.1415
    Lx=2.0*pi
    Ly=pi/2
    Lz=2.0
    Iv=((Run['variable']['I']-1.0)*(-1.0)*Lx*Ly)
    rb=Run['variable']['db']/2.0
    Iv.symetriser([1])

    #Nb=[0.2825, 0.3175, 0.3525, 0.3875, 0.4225, 0.4575, 0.4925, 0.5275, 0.5625, 0.5975, 0.6325, 0.6675, 0.7025, 0.7375, 0.7725, 0.8075, 0.8425, 0.8775, 0.9125, 0.9475, 0.9825]
    #Nb=[0.24, 0.3, 0.3475, 0.3875, 0.4225, 0.4575, 0.4925, 0.5275, 0.5625, 0.5975, 0.6325, 0.6675, 0.7025, 0.7375, 0.7725, 0.8075, 0.8425, 0.8775, 0.9125, 0.9475, 0.9825]
     #   1.76  1.7  1.6525  1.6125  1.5775  1.5425  1.5075  1.4725  1.4375  1.4025  1.3675  1.3325  1.2975  1.2625  1.2275  1.1925  1.1575  1.1225  1.0875  1.0525  1.0175
    Nb=[0.3875, 0.4225, 0.4575, 0.4925, 0.5275, 0.5625, 0.5975, 0.6325, 0.6675, 0.7025, 0.7375, 0.7725, 0.8075, 0.8425, 0.8775, 0.9125, 0.9475, 0.9825]


    Ivnew=Iv*0.0
    tracer([Iv, Ivnew], 'Ivdeb')
    for i in range(len(Nb)):

    	xb=Nb[i]
    	xbm=Lz-Nb[i]
        if xb==1.0:
		xbm=500000
        for j in range(len(Iv._ax[0,:])):
                x=Iv._ax[0,j]
        	if x>(xb-rb) and x<(xb+rb):
		        Ivnew._npa[0,j]=Ivnew._npa[0,j]+(pi*(rb*rb-(x-xb)*(x-xb)))
        	if x>(xbm-rb) and x<(xbm+rb):
		        Ivnew._npa[0,j]=Ivnew._npa[0,j]+(pi*(rb*rb-(x-xbm)*(x-xbm)))

    tracer([Iv, Ivnew], 'Iv')       
    
    return



def outDataForSrcNeptune(Run_imput):
    """ 
    print all the profils in a .txt files  
    @param a dictionnary  
    """   
    sym_tens=[1,1,-1,1,1,-1,-1,-1,1]
    Run=Run_imput
 #   for cle in Run.keys():
 #       for cle2 in Run[cle].keys():
 #           if isinstance(Run[cle][cle2], ScalarField):
 #               Run[cle][cle2].symetriser()
 #           elif isinstance(Run[cle][cle2], VectorField):
 #               Run[cle][cle2].symetriser()    
 #           elif isinstance(Run[cle][cle2], Tensor2Field):
 #               Run[cle][cle2].symetriser(sym_tens)                   
                
                
    L=len(Run['variable']['I']._ax[0,:])/2
    fichier = open(Run['variable']['name']+"_out_ligne.txt", 'w')
    fichier.write('# Coordonnee_z['+str(L)+'] = {')
    for j in range(L):
        if j==L-1:
            fichier.write(str(Run['variable']['I']._ax[0,j])+'};\n')
        else:
            fichier.write(str(Run['variable']['I']._ax[0,j])+', ') 
    
    for cle in Run.keys():
        for cle2 in Run[cle].keys():
            if isinstance(Run[cle][cle2], ScalarField):
                fichier.write('# '+cle+'_'+cle2+'['+str(L)+'] = {')
                for j in range(L):
                    if j==L-1:
                        fichier.write(str(Run[cle][cle2]._npa[0,j])+'};\n')
                    else:
                        fichier.write(str(Run[cle][cle2]._npa[0,j])+', ') 
            if isinstance(Run[cle][cle2], VectorField):
                for k in range(3):
                    fichier.write('# '+cle+'_'+cle2+'_'+str(k)+'['+str(L)+'] = {')
                    for j in range(L):
                        if j==L-1:
                            fichier.write(str(Run[cle][cle2]._npa[k,j])+'};\n')
                        else:
                            fichier.write(str(Run[cle][cle2]._npa[k,j])+', ') 
                    
                  
                    
            if isinstance(Run[cle][cle2], Tensor2Field):
                for k in range(3):
                    for l in range(3):
                        fichier.write('# '+cle+'_'+cle2+'_'+str(k)+str(l)+'['+str(L)+'] = {')
                        for j in range(L):
                            if j==L-1:
                                fichier.write(str(Run[cle][cle2]._npa[k,l,j])+'};\n')
                            else:
                                fichier.write(str(Run[cle][cle2]._npa[k,l,j])+', ') 
                
 
    return
    
    return
