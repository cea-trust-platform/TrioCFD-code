
import matplotlib.pyplot as plt
global compteur
compteur=0
import Tools as tbx

# Methode tracerJFM extraite de Tools.py
def tracerJFM(field, name="blabla"):


	print "tracer JFM pictures"
	color2="#8dd3c7"
	color4="#ffffb3"
	color1="#bebada"
	color3="#fb8072"

        ###### figure 1
	scalepng=1.5
	f, (ax) = plt.subplots()
	markeredgewidths=[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.]
	ax.set_xlim([0.,1.0], emit=False)
        ###### void fraction #######

        #style=['-', '-', '-', '-', '--', '--', '--']
	#dash=[(1000,1),(1000,1),(1000,1),(1000,1),(1,1),(1,1),(1,1)]
	#color=['k', 'k', 'k', 'k', 'r', 'r', 'b']
        #marker=[u'.', u'd', u'.', u'.', u'.', u'd', u'x']
	#label=['Sb (TrioCFD)', 'D (TrioCFD)', 'Ss', 'B', 'Sb (LT2008)', 'D (LT2008)', 'D+Ss (TrioCFD)']
        #markerfacecols = ['k', 'k', 'k', 'white', 'r', 'r', 'b']
        #markersizes=[4., 2., scalepng*1.5, 4., 4., scalepng*1.5, scalepng*1.5]
        #ax.set_xlim([0.,1.0], emit=False)
        #ax.set_ylim([0.,0.22], emit=False)
	#ax.set_xlabel(r"$\mathbf{y/h}$", fontsize=8)
	#ax.set_ylabel(r"$\mathbf{\alpha_v}$", fontsize=8)


        ######### URMS ###################
        #style=['-', '-', '-', '-', '-', '--', '--']
	#dash=[(1000,1),(1000,1),(1000,1),(1000,1), (1000,1),(1,1),(1,1)]
	#color=['k', 'k', 'k', 'k', 'k', 'r', 'r']
        #marker=[u'.', u'd', u'.', u'.', None, u'.', u'd']
	#label=['Sb (TrioCFD)', 'D (TrioCFD)', 'Ss', 'B', 'SP', 'Sb (LT2008)', 'D (LT2008)']
        #markerfacecols = ['k', 'k', 'k', 'white', 'k', 'r', 'r']
        #markersizes=[4., 2., scalepng*1.5, 4., scalepng*1.5, 4., scalepng*1.5]
        #ax.set_xlim([0.,1.0], emit=False)
        #ax.set_ylim([-0.1,1.], emit=False)
	#ax.set_xlabel(r"$\mathbf{y/h}$", fontsize=8)
	#ax.set_ylabel(r"$\mathbf{-{\overline{u'v'}}/{u_\tau^2}}$", fontsize=8)



	######## qdm z single-fluid D et SP ####
        #style=['-', '-', '-', '-', '-', '--']
	#dash=[(1000.,1.),(0.7,0.7),(1.5,1.5),(1000.,1.), (0.7,0.7),(1.5,1.5)]
	#color=['r', 'y', 'm', 'r', 'y', 'm']
        #marker=[u'd', u'd', u'd', None, None, None]
	#label=['r$DR_{ij}$', 'D (TrioCFD)', 'Ss', 'B', 'SP', 'Sb (LT2008)', 'D (LT2008)']
        #markerfacecols = ['r', 'y', 'm', 'r', 'y', 'm']
        #markersizes=[2.,2.,2., scalepng*1.5, scalepng*1.5, scalepng*1.5]
        #ax.set_xlim([0.,1.0], emit=False)
        #ax.set_ylim([-0.002,0.0025], emit=False)
	#ax.set_xlabel(r"$\mathbf{y/h}$", fontsize=8)
 	#ax.plot([], [], linestyle='-', dashes=(1000.,1.), linewidth=0.4, label='$\mathbf{DR_{ij}}$', color='r')
 	#ax.plot([], [], linestyle='-', dashes=(0.7,0.7), linewidth=0.4, label='$\mathbf{DP}$', color='y')
 	#ax.plot([], [], linestyle='-', dashes=(1.5,1.5), linewidth=0.4, label='$\mathbf{M_\sigma}$', color='m')
	#ax.plot([], [], linestyle='-', label='$\mathbf{D}$', color='k', linewidth=0., marker=u'd', ms=2., markerfacecolor='k', markeredgecolor = 'k', markeredgewidth=1.0)
 	#ax.plot([], [], linestyle='-', dashes=(1000.,1.), linewidth=0.4, label='$\mathbf{SP}$', color='k')


	######## qdm z single-fluid Ss et B ####
        #style=['-', '-', '-', '-', '-', '-', '-', '-', '-']
	#dash=[(1000.,1.),(0.7,0.7),(1.5,1.5),(1000.,1.), (0.7,0.7),(1.5,1.5),(1000.,1.), (0.7,0.7),(1.5,1.5)]
	#color=['r', 'y', 'm', 'r', 'y', 'm', 'r', 'y', 'm']
        #marker=[u'.', u'.', u'.',u'.', u'.', u'.', None, None, None]
        #markerfacecols = ['r', 'y', 'm', 'white', 'white', 'white', 'r', 'y', 'm']
        #markersizes=[scalepng*1.5, scalepng*1.5, scalepng*1.5,4.,4.,4., scalepng*1.5, scalepng*1.5, scalepng*1.5]
        #ax.set_xlim([0.,1.0], emit=False)
        #ax.set_ylim([-0.01,0.01], emit=False)
	#ax.set_xlabel(r"$\mathbf{y/h}$", fontsize=8)
 	#ax.plot([], [], linestyle='-', dashes=(1000.,1.), linewidth=0.4, label='$\mathbf{DR_{ij}}$', color='r')
 	#ax.plot([], [], linestyle='-', dashes=(0.7,0.7), linewidth=0.4, label='$\mathbf{DP}$', color='y')
 	#ax.plot([], [], linestyle='-', dashes=(1.5,1.5), linewidth=0.4, label='$\mathbf{M_\sigma}$', color='m')
	#ax.plot([], [], linestyle='-', label='$\mathbf{Ss}$', color='k', linewidth=0., marker=u'.', ms=scalepng*1.5, markerfacecolor='k', markeredgecolor = 'k', markeredgewidth=1.0)
	#ax.plot([], [], linestyle='-', label='$\mathbf{B}$', color='k', linewidth=0., marker=u'.', ms=4., markerfacecolor='white', markeredgecolor = 'k', markeredgewidth=1.0)
 	#ax.plot([], [], linestyle='-', dashes=(1000.,1.), linewidth=0.4, label='$\mathbf{SP}$', color='k')
	

	###### qdm bi fluid D ###########
	#style=['-', '-', '-']
	#dash=[(1000.,1.),(0.7,0.7),(1.5,1.5)]
	#color=['r', 'y', 'm']
        #marker=[u'd', u'd', u'd']
        #markerfacecols = ['r', 'y', 'm']
        #markersizes=[scalepng*1.5, scalepng*1.5, scalepng*1.5]
	#label=['$\mathbf{M_\sigma}$', '$\mathbf{M_l}$', '$\mathbf{M_v}$']
        #ax.set_xlim([0.,1.0], emit=False)
        #ax.set_ylim([-0.0,0.0018], emit=False)
	#ax.set_xlabel(r"$\mathbf{y/h}$", fontsize=8)


	###### qdm bi fluid Ss and B ###########
	#style=['-', '-', '-','-', '-', '-']
	#dash=[(1000.,1.),(0.7,0.7),(1.5,1.5),(1000.,1.),(0.7,0.7),(1.5,1.5)]
	#color=['r', 'y', 'm','r', 'y', 'm']
        #marker=[u'.', u'.', u'.',u'.', u'.', u'.']
        #markerfacecols = ['r', 'y', 'm', 'white', 'white', 'white']
        #markersizes=[scalepng*1.5, scalepng*1.5, scalepng*1.5, 4.,4.,4.]
        #ax.set_xlim([0.,1.0], emit=False)
        #ax.set_ylim([0.0,0.01], emit=False)
	#ax.set_xlabel(r"$\mathbf{y/h}$", fontsize=8)
 	#ax.plot([], [], linestyle='-', dashes=(1000.,1.), linewidth=0.4, label='$\mathbf{M_\sigma}$', color='r')
 	#ax.plot([], [], linestyle='-', dashes=(0.7,0.7), linewidth=0.4, label='$\mathbf{M_l}$', color='y')
 	#ax.plot([], [], linestyle='-', dashes=(1.5,1.5), linewidth=0.4, label='$\mathbf{M_v}$', color='m')
	#ax.plot([], [], linestyle='-', label='$\mathbf{Ss}$', color='k', linewidth=0., marker=u'.', ms=scalepng*1.5, markerfacecolor='k', markeredgecolor = 'k', markeredgewidth=1.0)
	#ax.plot([], [], linestyle='-', label='$\mathbf{B}$', color='k', linewidth=0., marker=u'.', ms=4., markerfacecolor='white', markeredgecolor = 'k', markeredgewidth=1.0)


	### PI avant apres #########
	#style=['-', '-', '-', '-', '-', '-', '-', '-']
	#dash=[(1000.,1.),(1000.,1.),(0.7,0.7),(0.7,0.7),(1.5,1.5),(1.5,1.5),(3.,3.),(3.,3.)]
	#color=['r', 'r', 'm', 'm', 'y', 'y', 'b', 'b']
        #marker=[None, u'x', None, u'x', None, u'x', None, u'x']
	#label=['$\phi_{11}$ direct', '$\phi_{11}$ corrected','$\phi_{22}$ direct', '$\phi_{22}$ corrected','$\phi_{33}$ direct', '$\phi_{33}$ corrected','$\phi_{12}$ direct']
	#label=['$\Pi_{11}$ direct', '$\Pi_{11}$ corrected','$\Pi_{22}$ direct', '$\Pi_{22}$ corrected','$\Pi_{33}$ direct', '$\Pi_{33}$ corrected','$\Pi_{12}$ direct','$\Pi_{12}$ corrected']
        #markerfacecols = ['r', 'r', 'y', 'y', 'm', 'm', 'b', 'b']
        #markersizes=[scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5]
        #markeredgewidths=[1.,1.,1.,1.,0.7,0.7,0.4,0.4]
        #ax.set_xlim([0.,1.0], emit=False)
        #ax.set_ylim([-0.00013,0.0011], emit=False)
	#ax.set_xlabel(r"$\mathbf{y/h}$", fontsize=8)


	### model redistribution ###
	#style=['-', '-', '-', '-', '-', '-', '-']
	#dash=[(1000.,1.),(1000.,1.),(1000.,1.),(1000.,1.),(0.7,0.7),(1.5,1.5), (2.5,2.5)]
	#color=['k', 'k', 'k', 'r', 'y', 'm', 'g']
        #marker=[u'*', u'.', u'x', None, None, None, None]
	#label=['$(\phi_{11}+\Pi_{11})/\Pi_{11}$','$(\phi_{22}+\Pi_{22})/\Pi_{11}$','$(\phi_{33}+\Pi_{33})/\Pi_{11}$', '$3/5$', '$1/2$', '$1/4$', '$1/5$' ]
	#markerfacecols=['k', 'k', 'k', 'r', 'y', 'm', 'g']
        #markersizes=[scalepng*3, scalepng*3, scalepng*1.5, scalepng*3, scalepng*3, scalepng*3, scalepng*3]
        #markeredgewidths=[0.0,0.0,1.0,0.0,0.0,0.0,0.0]
        #ax.set_xlim([0.,1.0], emit=False)
        #ax.set_ylim([0.,1.3], emit=False)
	#ax.set_xlabel(r"$\mathbf{y/h}$", fontsize=8)


	### model interfacial prod
	#style=['-', '-', '-', '-', '-', '-', '-']
	#dash=[(1000.,1.),(1000.,1.),(1000.,1.),(1000.,1.),(0.7,0.7),(1.5,1.5), (2.5,2.5)]
	#color=['k', 'k', 'k', 'r', 'y', 'm', 'g']
        #marker=[u'*', None, u'x', None, None, None, None]
	#label=['$\Pi_{11}$', r"$\alpha_l\alpha_v\Delta\rho gu_r$"]
	#markerfacecols=['k', 'k', 'k', 'r', 'y', 'm', 'g']
        #markersizes=[scalepng*3, scalepng*3, scalepng*1.5, scalepng*3, scalepng*3, scalepng*3, scalepng*3]
        #markeredgewidths=[0.0,0.0,1.0,0.0,0.0,0.0,0.0]
        #ax.set_xlim([0.,1.0], emit=False)
        #ax.set_ylim([-0.0002,0.001], emit=False)
	#ax.set_xlabel(r"$\mathbf{y/h}$", fontsize=8)

	### model epsilon
	#style=['-', '-', '-', '-', '-', '-', '-']
	#dash=[(1000.,1.),(1000.,1.),(1000.,1.),(1000.,1.),(0.7,0.7),(1.5,1.5), (2.5,2.5)]
	#color=['k', 'k', 'k', 'r', 'y', 'm', 'g']
        #marker=[u'*', u'.', u'x', None, None, None, None]
	#label=[r'$\epsilon_{11}/R_{11}$', r'$\epsilon_{22}/R_{22}$', r'$\epsilon_{33}/R_{33}$', r'$u_r/2d_b$', r'$u_r/3d_b$', r'$u_r/4d_b$']
	#markerfacecols=['k', 'k', 'k', 'r', 'y', 'm', 'g']
        #markersizes=[scalepng*3, scalepng*3, scalepng*1.5, scalepng*3, scalepng*3, scalepng*3, scalepng*3]
        #markeredgewidths=[0.0,0.0,1.0,0.0,0.0,0.0,0.0]
        #ax.set_xlim([0.,1.0], emit=False)
        #ax.set_ylim([-0.7,0.], emit=False)
	#ax.set_xlabel(r"$\mathbf{y/h}$", fontsize=8)

	#### compa model risso model antoine
	#style=['-', '-', '-', '-', '-', '-', '-']
	#dash=[(1000.,1.),(0.7,0.7),(1.5,1.5),(0.7,0.7),(1.5,1.5)]
	#color=['k', 'r', 'y', 'r', 'y']
        #marker=[u'*', u'.', u'.', u'x', u'x']
	#label=[r'$\mathbf{B-Ss}$', r'$k_{budget} L_w=4d_b$', r'$k_{budget} L_w=3d_b$', r'$k_{Risso} L_w=4d_b$', r'$k_{Risso} L_w=3d_b$']
	#markerfacecols=['k', 'r', 'y', 'r', 'y']
        #markersizes=[scalepng*3, scalepng*3, scalepng*3, scalepng*1.5, scalepng*1.5]
        #markeredgewidths=[0.0,0.0,0.0,1.0,1.0]
        #ax.set_xlim([0.,1.0], emit=False)
        #ax.set_ylim([0.0,0.005], emit=False)
	#ax.set_xlabel(r"$\mathbf{y/h}$", fontsize=8)

	if("JFMequilibreDNS" in name):
		style=['-', '-', '-', '-', '-', '-', '-', '-', '-', '-']
		dash=[(1000.,1.), (1000.,1.) ,(1000.,1.),(1000.,1.),(1000.,1.),(1000.,1000.),(1.,1.),(1000.,1.),(1.,1.),(2.,2.)]
		color=['k', 'g', 'r', 'b', 'm', 'c', 'k', 'g', 'r', 'b']
		marker=[None, u'd', u'.', u'x', u'd', u'+', None, '_', None, None]
		label=[r'$\mathbf{\dot M_l=\Sigma M_l}$',r'$\mathbf{\dot M_l}$ DNS', r'$\mathbf{M_l^\sigma}+\mathbf{M_l^P}$', r'$\mathbf{M_l^\Pi}$', r'$\mathbf{M_l^\Re}$', r'$\mathbf{M_l^\tau}$', r'$\alpha_vu_r\wedge(\nabla\wedge u_l)$', r'$-\alpha_vu_r\wedge(\nabla\wedge u_l)$']
		if("JFMequilibreDNSx" in name):
			label=[r'$\mathbf{\dot M_l}$', r'$\mathbf{M_l^{LD}}$', r'$\mathbf{M_l^D}$', r'$\mathbf{M_l^{TD}}$', r'$\mathbf{M_l^\tau}$', r'$\mathbf{M_l^{extra}}$', r'$\alpha_v|u_r|u_r$', r'$\alpha_v|u_r|u_r$']

		if("JFMequilibreDNSz" in name):
			label=[r'$\mathbf{M_l^{RANS}-\frac{\alpha_l}{\alpha_v}M^L_{\nabla P}}$', r'$\frac{\alpha_l}{\alpha_v}\mathbf{M^{LD}}$', r'$-\mathbf{M_\Re^L}$', r'$\mathbf{-M^{TD}}$', r'$\mathbf{-M^\tau}$', r'$\mathbf{-M^{extra}}$', r'$-\alpha_l\nabla\overline{p_l^{SP}}-\frac{\alpha_l}{\alpha_v}\mathbf{M^L_{\nabla P}}$', r'$\alpha_v|u_r|u_r$']
			#label=[r'$\mathbf{M_l^{RANS}-\frac{\alpha_l}{\alpha_v}M^L_{\nabla P}}$', r'$\frac{\alpha_l}{\alpha_v}\mathbf{M^{LD}}$', r'$-\mathbf{M_\Re^L}$', r'$\mathbf{-M^{TD}}$', r'$\mathbf{-M^\tau}$', r'$\mathbf{-M^{extra}}$', r'$-\alpha_l\nabla\overline{P_l^{SP}}-\frac{\alpha_l}{\alpha_v}\mathbf{M^L_{\nabla P}}$', r'$\alpha_v|u_r|u_r$']
		
		if("JFMequilibreDNSx_vapeur" in name):
			label=[r'$\mathbf{M_v^{RANS}}-\alpha_v\nabla\left(\overline{p_l^b}-\widetilde{p_l}\right)$', r'$\mathbf{M_v^{LD}}$', r'$\mathbf{M_v^D}$', r'$\mathbf{M_v^{TD}}$', r'$\mathbf{M_v^\tau}$', r'$\mathbf{M_v^{extra}}$', r'$\alpha_v|u_r|u_r$', r'$\alpha_v|u_r|u_r$']

		if("JFMequilibreDNSz_vapeur" in name):
			label=[r'$\mathbf{M_v^{RANS}-M^L_{\nabla P}}$', r'$\mathbf{M^{LD}}$', r'$\mathbf{M^L_\Re}$', r'$\mathbf{M^{TD}}$', r'$\mathbf{M^\tau}$', r'$\mathbf{M^{extra}}$', r'$\alpha_v|u_r|u_r$', r'$\alpha_v|u_r|u_r$']


		markerfacecols=['None','none', 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None']
		markersizes=[scalepng*3, scalepng*2, scalepng*3.5, scalepng*3, scalepng*2, scalepng*3, scalepng*3, scalepng*3, scalepng*3, scalepng*3]
		markeredgewidths=[0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3]



		if("JFMequilibreDNSx" in name):
			#ax.set_ylim([-0.01,0.01], emit=False)
			ax.set_ylim([-0.015,0.04], emit=False) #S180
			#ax.set_ylim([-0.01,0.03], emit=False) #S180g8
			#ax.set_ylim([-0.002,0.01], emit=False) #D180
			#ax.set_ylim([-0.004,0.03], emit=False) #D180g8
			#ax.set_ylim([-0.001,0.004], emit=False) #D127

		if("JFMequilibreDNSz" in name):
			#ax.set_ylim([-0.06,0.06], emit=False) #S180
			#ax.set_ylim([-0.02,0.1], emit=False) #S180g8
			ax.set_ylim([-0.0015,0.008], emit=False) #D127
			#ax.set_ylim([-0.002,0.004], emit=False) #D180
			#ax.set_ylim([-0.005,0.02], emit=False) #D180g8

		if("JFMequilibreDNSx_vapeur" in name):
			#ax.set_ylim([-0.01,0.01], emit=False)
			ax.set_ylim([-0.02,0.015], emit=False) #S180
			#ax.set_ylim([-0.03,0.01], emit=False) #S180g8
			#ax.set_ylim([-0.01,0.002], emit=False) #D180
			#ax.set_ylim([-0.03,0.003], emit=False) #D180g8
			#ax.set_ylim([-0.004,0.001], emit=False) #D127

		if("JFMequilibreDNSz_vapeur" in name):
			#ax.set_ylim([-0.0003,0.0005], emit=False) #D180
			#ax.set_ylim([-0.0006,0.001], emit=False) #D180g8
			#ax.set_ylim([-0.003,0.003], emit=False) #S180
			#ax.set_ylim([-0.002,0.003], emit=False) #S180g8
			ax.set_ylim([-0.0001,0.00006], emit=False) #D127


		ax.set_xlabel(r"$\mathbf{y/h}$", fontsize=8)


	if("JFMequilibre_force" in name):
		style=['-', '-', '-', '-', '-', '-', '-', '-', '-', '-']
		dash=[(1000.,1.), (1000.,1.) ,(1000.,1.),(1000.,1.),(1000.,1.),(1000.,1.),(1000.,1.),(1.,1.), (2.,2.)]
		if("Bois" in name):
			dash=[(1000.,1.),(1000.,1.),(1000.,1.),(1000.,1.),(1000.,1.),(1.,1.), (2.,2.)]
		color=['k', 'g', 'r', 'b', 'k', 'g', 'r', 'b', 'k', 'g']
		if("Bois" in name):
			color=['k', 'g', 'r', 'b', 'k', 'b', 'r']
		marker=[u'd', u'.', u'x', u'*', u'd', u'.', u'*', None, None, None]
		if("Bois" in name):
			marker=[u'd', u'.', u'x', u'*', u'd', None, None, None]



		if("z" in name):
			label=[r'$\mathbf{M^{LD}}$', r'$\mathbf{M^{TD}}$',r'$\mathbf{M^\tau}$',r'$\mathbf{M^\Pi}$', r'residue',r"$-\alpha_v\nabla{\overline{p_l^{SP}}}+\mathbf{M^L_{\nabla P}}$",r'$\mathbf{M_\Re^{L}|_{lam}}+\mathbf{M_\Re^{L}|_{turb}}$', r'$C_L^*\rho_l\alpha_vu_r\frac{\partial\overline{u_l}^l}{\partial y}$', r'$-C_L^*\rho_l\alpha_vu_r\frac{\partial\overline{u_l}^l}{\partial y}$' ]
			#label=[r'$\mathbf{M^{LD}}$', r'$\mathbf{M^{TD}}$',r'$\mathbf{M^\tau}$',r'$\mathbf{M^\Pi}$', r'residue',r"$-\alpha_v\nabla{\overline{P_l^{SP}}}+\mathbf{M^L_{\nabla P}}$",r'$\mathbf{M_\Re^{L}|_{lam}}+\mathbf{M_\Re^{L}|_{turb}}$', r'$C_L^*\rho_l\alpha_vu_r\frac{\partial\overline{U_l}^l}{\partial y}$', r'$-C_L^*\rho_l\alpha_vu_r\frac{\partial\overline{U_l}^l}{\partial y}$' ]
		if("Bois" in name):
			label=[r'$\mathbf{M^{LD}}$', r'$\mathbf{M^{TD}}$',r"$-\alpha_v\nabla{\overline{p_l^{SP}}}+\mathbf{M^L_{\nabla P}}$",r'$\mathbf{M_\Re^{L}|_{lam}}+\mathbf{M_\Re^{L}|_{turb}}$', r'residue', r'$C_L^*\rho_l\alpha_vu_r\frac{\partial\overline{u_l}^l}{\partial y}$', r'$-C_L^*\rho_l\alpha_vu_r\frac{\partial\overline{u_l}^l}{\partial y}$' ]



		if("x" in name):
			label=[r'$\mathbf{M^{LD}}$', r'$\mathbf{M^{TD}}$',r'$\mathbf{M^\tau}$',r'$\mathbf{M^\Pi}$',r'residue',r'$-\alpha_v\nabla{\overline{p_l^{,}}}^l$',r'$\mathbf{M^{D}}$', r'$C_D^*\frac{3}{4d_b}\rho_l\alpha_v|u_r|u_r$']
			#label=[r'$\mathbf{M^{LD}}$', r'$\mathbf{M^{TD}}$',r'$\mathbf{M^\tau}$',r'$\mathbf{M^\Pi}$',r'residue',r'$-\alpha_v\nabla{\overline{P_l^{,}}}^l$',r'$\mathbf{M^{D}}$', r'$C_D^*\frac{3}{4d_b}\rho_l\alpha_v|u_r|u_r$']
		markerfacecols=['None','None', 'None', 'None', 'k', 'g', 'r', 'b', 'k', 'g']
		markersizes=[scalepng*2, scalepng*3.5, scalepng*2.5, scalepng*3.5, scalepng*2, scalepng*3.5, scalepng*3., scalepng*3, scalepng*3, scalepng*3]
		markeredgewidths=[0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3]


		if("x" in name):
			#ax.set_ylim([-0.01,0.01], emit=False)
			#ax.set_ylim([-0.015,0.015], emit=False) #S180
			ax.set_ylim([-0.04,0.04], emit=False) #S180g8
			#ax.set_ylim([-0.006,0.006], emit=False) #D180
			#ax.set_ylim([-0.035,0.035], emit=False) #D180g8
			#ax.set_ylim([-0.005,0.005], emit=False) #D127

		if("z" in name):
			#ax.set_ylim([-0.06,0.06], emit=False) #S180
			#ax.set_ylim([-0.02,0.1], emit=False) #S180g8
			#ax.set_ylim([-0.0005,0.003], emit=False) #D127
			ax.set_ylim([-0.002,0.004], emit=False) #D180
			#ax.set_ylim([-0.005,0.02], emit=False) #D180g8
			if("vapeur" in name):
				#ax.set_ylim([-0.0003,0.0003], emit=False) #D180
				#ax.set_ylim([-0.0006,0.0006], emit=False) #D180g8
				#ax.set_ylim([-0.003,0.003], emit=False) #S180
				ax.set_ylim([-0.002,0.002], emit=False) #S180g8
				#ax.set_ylim([-0.0001,0.0001], emit=False) #D127	



		if ("lift" in name):
			ax.set_ylim([-0.1,0.1], emit=False) #D127		

		ax.set_xlabel(r"$\mathbf{y/h}$", fontsize=8)

	if("JFMequilibre_force_z_vapeur_2" in name):
		style=['-', '-', '-', '-', '-', '-', '-', '-', '-', '-']
		dash=[(1000.,1.), (1000.,1.) ,(1000.,1.),(1000.,1.),(1000.,1.),(1000.,1.),(1000.,1.),(1.,1.), (2.,2.), (3.,1.)]
		color=['k', 'g', 'r', 'b', 'k', 'g', 'r', 'b', 'k', 'g']
		marker=[u'd', u'.', u'x', u'*', u'd', u'.', u'*', None, None, None]

		label=[r'interface+inter pressure',r'reste pressure',r'turb disp',r'gradP',r'lift',r'residue', r'maintien']
		markerfacecols=['None','None', 'None', 'None', 'k', 'g', 'r', 'b', 'k', 'g']
		markersizes=[scalepng*2, scalepng*3.5, scalepng*2.5, scalepng*3.5, scalepng*2, scalepng*3.5, scalepng*3.5, scalepng*3, scalepng*3, scalepng*3]
		markeredgewidths=[0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3]


		ax.set_ylim([-70.,50.], emit=False) #D127	


		if ("lift" in name):
			ax.set_ylim([-0.1,0.1], emit=False) #D127		

		ax.set_xlabel(r"$\mathbf{y/h}$", fontsize=8)

	if("C_p" in name):
		style=['-', '-', '-', '-', '-', '-', '-']
		dash=[(1000.,1.),(2.,2.),(2.,2.),(2.,2.),(2.,2.),(2.,2.)]
		color=['k', 'r', 'b', 'g', 'm', 'c']
		marker=[None, u'.', u'x', u'd', u'*', u'+']
		label=[r'DNS', r'$C=1$', r'$C=0.75$', r'$C=0.5$', r'$C=0.25$', r'$C=0$']
		markerfacecols=['None', 'None', 'None', 'None', 'None', 'None']
		markersizes=[scalepng*3, scalepng*3, scalepng*2, scalepng*2, scalepng*3, scalepng*2]
		markeredgewidths=[0.3,0.3,0.3,0.3,0.3,0.3]


		if("l" in name):
			#ax.set_ylim([-0.03,0.01], emit=False) #S180g8
			#ax.set_ylim([-0.03,0.01], emit=False) #S180
			#ax.set_ylim([-0.035,0.01], emit=False) #D80g8
			ax.set_ylim([-0.003,0.001], emit=False) #D127
			#ax.set_ylim([-0.005,0.0015], emit=False) #D180
		if("v" in name):
			#ax.set_ylim([0.0,3.], emit=False) #S80g8
			#ax.set_ylim([0.0,0.5], emit=False) #S80
			#ax.set_ylim([0.0,0.5], emit=False) #D80g8
			ax.set_ylim([0.0,0.07], emit=False) #D127
			#ax.set_ylim([0.0,0.07], emit=False) #D180


		#ax.set_ylim([-0.001,0.001], emit=False) #S180

		ax.set_xlabel(r"$\mathbf{y/h}$", fontsize=8)


	if("JFMequilibreDNSzapriori" in name):
		style=['-', '-', '-', '-', '-', '-', '-']
		dash=[(1000.,1.),(1000.,1.),(1000.,1.),(2.,2.),(2.,2.),(2.,2.)]
		color=['k', 'b', 'r', 'k', 'b', 'r']
		marker=[None, None, None, u'd', u'd', u'd']
                #label[r"$\mathbf{M_l^P}$", r"$\mathbf{M_l^P}$", r"$\mathbf{M_l^P}$", r"$\mathbf{M_l^P}$", r"$\mathbf{M_l^P}$", r"$\mathbf{M_l^P}$"]
                label=[r"$\mathbf{M^P}$ DNS", r"$\mathbf{M^\sigma}$ DNS",r"$\mathbf{M^{LD}}$ DNS", r"$\mathbf{M^P}$ model", r"$\mathbf{M^\sigma}$ model", r"$\mathbf{M^{LD}}$ model"]
		markerfacecols=['None', 'None', 'None', 'None', 'None', 'None']
		markersizes=[scalepng*2.5, scalepng*2.5, scalepng*2.5, scalepng*2.5, scalepng*2.5, scalepng*2.5]
		markeredgewidths=[0.3,0.3,0.3,0.3,0.3,0.3]
		#ax.set_ylim([-0.005,0.006], emit=False) #D127
		#ax.set_ylim([-0.015,0.015], emit=False) #D180
		#ax.set_ylim([-0.04,0.05], emit=False) #D180g8
		ax.set_ylim([-0.6,0.8], emit=False) #S180g8
		#ax.set_ylim([-0.6,0.8], emit=False) #S180
		ax.set_xlabel(r"$\mathbf{y/h}$", fontsize=8)


	if("JFMI" in name):
		style=['-', '-', '-', '-', '-', '-', '-']
		dash=[(1000.,1.),(1000.,1.),(1000.,1.),(1000.,1.),(1000.,1.),(3.,3.)]
		color=['k', 'b', 'r', 'g', 'm', 'k']
		marker=['None', u'd', u'.', u'x', u'*', None]
                label=[r"D127", r"D180", r"D180g8", r"S180",r"S180g8", ""]
		markerfacecols=['None', 'None', 'None', 'None', 'None', 'None']
		markersizes=[scalepng*3, scalepng*2, scalepng*3, scalepng*2, scalepng*3, scalepng*3]
		markeredgewidths=[0.3,0.3,0.3,0.3,0.3,0.3]
		ax.set_ylim([0.,0.15], emit=False)
		ax.set_xlabel(r"$\mathbf{y/h}$", fontsize=8)
		ax.set_ylabel(r"$\alpha_v$", fontsize=8)


	if("JFM_I_nept" in name):
		style=['-', '-', '-', '-', '-', '-', '-']
		dash=[(1000.,1.),(1.,1.),(2.,2.),(3.,3.),(1.,3.)]
		color=['k', 'b', 'r', 'g', 'm', 'k']
		marker=['None', u'd', u'.', u'x', u'*', None]
                label=[r"DNS", r"$Eo_c=10$ standard", r"$Eo_c=2.5$ standard", r"$Eo_c=10$ model",r"$Eo_c=2.5$ model"]
		markerfacecols=['None', 'None', 'None', 'None', 'None', 'None']
		markersizes=[scalepng*3, scalepng*2, scalepng*3, scalepng*2, scalepng*3, scalepng*3]
		markeredgewidths=[0.3,0.3,0.3,0.3,0.3,0.3]
		ax.set_ylim([0.,0.13], emit=False)
		ax.set_xlabel(r"$\mathbf{y/h}$", fontsize=8)
		ax.set_ylabel(r"$\alpha_v$", fontsize=8)

	if("JFM_I_nept_sphe" in name):
		style=['-', '-', '-', '-', '-', '-', '-']
		dash=[(1000.,1.),(1.,1.),(3.,3.),(1.,3.)]
		color=['k', 'b', 'g', 'm', 'k']
		marker=['None', u'd', u'x', u'*', None]
                label=[r"DNS", r"$Eo_c=10$ standard", r"$Eo_c=10$ model"]
		markerfacecols=['None', 'None', 'None', 'None', 'None', 'None']
		markersizes=[scalepng*3, scalepng*2, scalepng*2, scalepng*3, scalepng*3]
		markeredgewidths=[0.3,0.3,0.3,0.3,0.3,0.3]
		ax.set_ylim([0.,0.08], emit=False)
		ax.set_xlabel(r"$\mathbf{y/h}$", fontsize=8)
		ax.set_ylabel(r"$\alpha_v$", fontsize=8)





	if("JFMWe" in name):
		style=['-', '-', '-', '-', '-', '-', '-']
		dash=[(1000.,1.),(1000.,1.),(1000.,1.),(1000.,1.),(1000.,1.),(3.,3.)]
		color=['k', 'b', 'r', 'g', 'm', 'k']
		marker=['None', u'd', u'.', u'x', u'*', None]
                label=[r"D127", r"D180", r"D180g8", r"S180",r"S180g8", r"$\mathbf{We_c=3}$"]
		markerfacecols=['None', 'None', 'None', 'None', 'None', 'None']
		markersizes=[scalepng*3, scalepng*2, scalepng*3, scalepng*2, scalepng*3, scalepng*3]
		markeredgewidths=[0.3,0.3,0.3,0.3,0.3,0.3]
		ax.set_ylim([0.,13.], emit=False)
		ax.set_xlabel(r"$\mathbf{y/h}$", fontsize=8)
		ax.set_ylabel(r"$\mathbf{We}$", fontsize=8)

	if("JFMvitesse" in name):
		style=['-', '-', '-', '-', '-', '-', '-']
		dash=[(1000.,1.),(1000.,1.),(1000.,1.),(1000.,1.),(1000.,1.),(3.,3.)]
		color=['k', 'b', 'r', 'g', 'm', 'k']
		marker=['None', u'd', u'.', u'x', u'*', None]
                label=[r"D127", r"D180", r"D180g8", r"S180",r"S180g8"]
		markerfacecols=['None', 'None', 'None', 'None', 'None', 'None']
		markersizes=[scalepng*3, scalepng*2, scalepng*3, scalepng*2, scalepng*3, scalepng*3]
		markeredgewidths=[0.3,0.3,0.3,0.3,0.3,0.3]
		ax.set_ylim([-0.1,1.2], emit=False)
		ax.set_xlabel(r"$\mathbf{y/h}$", fontsize=8)
		if("_l" in name):
			ax.set_ylabel(r"$\mathbf{\overline{U_l}^l}$", fontsize=8)
		if("_r" in name):
			ax.set_ylabel(r"$\mathbf{\overline{U_v}^v-\overline{U_l}^l}$", fontsize=8)


	if("Rij" in name):
		style=['-', '-', '-', '-', '-', '-', '-']
		dash=[(1000.,1.),(1000.,1.),(1000.,1.),(1000.,1.),(1000.,1.),(3.,3.)]
		color=['k', 'b', 'r', 'g', 'm', 'k']
		marker=['None', u'd', u'.', u'x', u'*', None]
                label=[r"D127", r"D180", r"D180g8", r"S180",r"S180g8"]
		markerfacecols=['None', 'None', 'None', 'None', 'None', 'None']
		markersizes=[scalepng*3, scalepng*2, scalepng*3, scalepng*2, scalepng*3, scalepng*3]
		markeredgewidths=[0.3,0.3,0.3,0.3,0.3,0.3]
		ax.set_xlabel(r"$\mathbf{y/h}$", fontsize=8)
  		if ("11" in name):
			ax.set_ylabel(r"$\mathbf{R_{11}}$", fontsize=8)
			ax.set_ylim([-0.01,0.1], emit=False)
			ax.set_xlim([0.,2.], emit=False)
  		if ("22" in name):
			ax.set_ylabel(r"$\mathbf{R_{22}}$", fontsize=8)
			ax.set_ylim([0.,0.015], emit=False)
  		if ("33" in name):
			ax.set_ylabel(r"$\mathbf{R_{22}}$", fontsize=8)
			ax.set_ylim([0.,0.015], emit=False)
  		if ("13" in name):
			ax.set_ylabel(r"$\mathbf{-R_{12}}$", fontsize=8)
			ax.set_ylim([0.,0.007], emit=False)


	if("Cp" in name):
		style=['-', '-', '-', '-', '-', '-', '-']
		dash=[(1000.,1.),(2.,2.),(2.,2.),(2.,2.),(2.,2.),(2.,2.)]
		color=['k', 'r', 'b', 'g', 'm', 'c']
		marker=[None, u'.', u'x', u'd', u'*', u'+']
                #label[r"$\mathbf{M_l^P}$", r"$\mathbf{M_l^P}$", r"$\mathbf{M_l^P}$", r"$\mathbf{M_l^P}$", r"$\mathbf{M_l^P}$", r"$\mathbf{M_l^P}$"]
                label=[r'DNS', r'$C=1$', r'$C=0.75$', r'$C=0.5$', r'$C=0.25$', r'$C=0$']
		markerfacecols=['None', 'None', 'None', 'None', 'None', 'None']
		markersizes=[scalepng*3, scalepng*3, scalepng*2, scalepng*2, scalepng*3, scalepng*2]
		markeredgewidths=[0.3,0.3,0.3,0.3,0.3,0.3]
		#ax.set_ylim([-0.005,0.01], emit=False) #D127
		#ax.set_ylim([-0.015,0.03], emit=False) #D180
		#ax.set_ylim([-0.04,0.085], emit=False) #D180g8
		#ax.set_ylim([-0.6,0.9], emit=False) #S180g8
		ax.set_ylim([-0.6,0.9], emit=False) #S180
		ax.set_xlabel(r"$\mathbf{y/h}$", fontsize=8)


		if("l" in name):
			#ax.set_ylim([-0.04,0.0], emit=False) #S180g8
			ax.set_ylim([-0.015,0.01], emit=False) #S180
			#ax.set_ylim([-0.035,0.01], emit=False) #D80g8
			#ax.set_ylim([-0.004,0.0005], emit=False) #D127
			#ax.set_ylim([-0.005,0.0015], emit=False) #D180
		if("v" in name):
			#ax.set_ylim([0.0,3.], emit=False) #S80g8
			#ax.set_ylim([0.0,0.5], emit=False) #S80
			#ax.set_ylim([0.0,0.4], emit=False) #D80g8
			ax.set_ylim([0.0,0.04], emit=False) #D127
			#ax.set_ylim([0.0,0.04], emit=False) #D180

	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.xaxis.set_ticks_position('bottom')
	ax.yaxis.set_ticks_position('left') 
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	#plt.ticklabel_format(style='plain', axis='y', scilimits=(0,0))
	plt.tight_layout()
	dpii=400
	
	for i in range(len(field)):
		print i
		### avec legende
        	ax.plot(field[i]._ax[0,:], field[i]._npa[0,:], linestyle=style[i], dashes=dash[i], linewidth=0.4, marker=marker[i], label=label[i], color=color[i], ms=markersizes[i], markevery=0.06, markerfacecolor=markerfacecols[i], markeredgecolor = color[i], markeredgewidth=markeredgewidths[i])
        	### sans legende
        	#ax.plot(field[i]._ax[0,:], field[i]._npa[0,:], linestyle=style[i], dashes=dash[i], linewidth=0.4, marker=marker[i], color=color[i], ms=markersizes[i], markevery=0.06, markerfacecolor=markerfacecols[i], markeredgecolor = color[i], markeredgewidth=1.0)
	
	f.set_size_inches(scalepng,scalepng)
	#f.set_size_inches(scalepng*1.5,scalepng) # faire des images rectangulaires
	box = ax.get_position()
	ax.set_position([box.x0, box.y0, box.width, box.height])
	ax.legend(loc='center left', bbox_to_anchor=(1.1, 0.5), frameon=True)
	inSizeLegend=int(scalepng*3.5)
	#plt.legend(loc='center right',prop={'size':inSizeLegend}, frameon=False)
	#plt.legend(loc=(0.5,0.15),prop={'size':inSizeLegend}, frameon=True)
	plt.legend(loc=0,prop={'size':inSizeLegend}, frameon=True, ncol=1)
	#plt.legend().set_visible(False) ### deactive la legende
	#plt.legend(loc=0,prop={'size':inSizeLegend}, frameon=True, ncol=2)
	#plt.legend(loc=0,prop={'size':inSizeLegend}, frameon=True)
	plt.savefig(name, bbox_inches='tight', dpi=dpii)
	plt.close(f)
	return




# Methode tracer extraite de Tools.py
def tracer(bfield, title, axe=None, textitle='', 
           ylim=None, xlim=None, 
           grid=False, legend="outside", markerSize=None, 
           log=False, couleur=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
           markevry=[3], yscale=1., xscale=1., scalepng=1.5, applatir=False, plotall=True, sym=[10], zerocentre=False, doubleLegend=False, casename=None):
    """ 
    Creation of a parametrized plotinitialization of a field object
    @param the field's list to plot (can be of any sub-field type)
    @param the name of the png file created
    @param can contain the list of components to be plotted. Default is None for all   
    @param the title of the plot (optional)
    @param limit of y axis. Default is None for all 
    @param limit of x axis. Default is None for all 
    @param Plot the grid background. Default is True
    @param legend ... pas au point
    @param markerSize. Default is an autoscale from the scale of the png file (scalepng).
    @param Plot with logarithmic scale. Default is False
    @param couleur=[3,2,3,1] means the first 3 variables are plotted in the same color, then the following two with an other etc..
    @param markevry
    @param dimenssionless tools for y axis
    @param dimenssionless tools for x axis
    @param scale of the out.png
    @param aspect ratio of the picture = 0.5. Default is False 

    """
    
    if len(markevry)==1 :
	cte=markevry[0]
	tmp=[cte,cte,cte,cte,cte,cte,cte,cte,cte,cte,cte,cte,cte,cte,cte,cte,cte,cte,cte,cte,cte]
	markevry=tmp
    print 'exporting graph '+title+'.png'
    global compteur
    compteur=compteur+1
    f=plt.figure(compteur)
    #f.tightlayout()
    #linestyles = ['_', '-', '--', ':']

    #markers = [u'.', u'x', None, u'.', u'.', u'.', u',', u'x', u'.', u'1', u'p', u'3', u'2', u'4', u'H', u'v', u'8', u'<']
    #markerfacecols = ['full', 'full', 'full', 'full', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r']
    #markersizes=[2., 2., scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5,scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5]

    linestyles = []
    #sens SP / D / Ss / B
    #markers = [None, u'd', u'.', u'.', u'.', None, u'.', u'.', u'.', u',', u'x', u'.', u'1', u'p', u'3', u'2', u'4', u'H', u'v', u'8', u'<']
    #markerfacecols = ['full', 'full', 'full', 'white', 'full', 'full', 'full', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r']
    #markersizes=[scalepng*1.5, 1.5, 2., 4., 2., scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5,scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5]

    # sens D /SP / Ss / B
    markers = [u'd',u'.',  u'.', u'.', u'.', None, u'.', u'.', u'.', u',', u'x', u'.', u'1', u'p', u'3', u'2', u'4', u'H', u'v', u'8', u'<']
    markerfacecols = ['white', 'full', 'full', 'white', 'full', 'full', 'full', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r']
    markersizes=[ 2.5, 4., 2., 4., 2., scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5,scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5]

    #sens Ss / B
    #markers = [u'.', u'.', u'.', None, u'.', u'.', u'.', u',', u'x', u'.', u'1', u'p', u'3', u'2', u'4', u'H', u'v', u'8', u'<']
    #markerfacecols = ['full', 'white', 'full', 'full', 'full', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r']
    #markersizes=[2., 4., 2., scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5,scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5]

    # sens B
    #markers = [u'.', u'.', None, u'.', u'.', u'.', u',', u'x', u'.', u'1', u'p', u'3', u'2', u'4', u'H', u'v', u'8', u'<']
    #markerfacecols = ['white', 'full', 'full', 'full', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r']
    #markersizes=[4., 2., scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5,scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5]

    # sens SP / croix / B 
    #markers = [None, u'x', u'.', u'.', None, u'.', u'.', u'.', u',', u'x', u'.', u'1', u'p', u'3', u'2', u'4', u'H', u'v', u'8', u'<']
    #markerfacecols = ['full', 'full', 'white', 'full', 'full', 'full', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r']
    #markersizes=[scalepng*1.5, scalepng*1.5, 4., 2., scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5,scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5]

    # sens SP / croix / B 
    #markers = [None, u'x', u'x', u'.', None, u'.', u'.', u'.', u',', u'x', u'.', u'1', u'p', u'3', u'2', u'4', u'H', u'v', u'8', u'<']
    #markerfacecols = ['full', 'full', 'full', 'full', 'full', 'full', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r']
    #markersizes=[scalepng*1.5, scalepng*1.5, scalepng*1.5, 2., scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5,scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5]

    # sens croix
    #markers = [u'x', u'x', u'x', u'.', None, u'.', u'.', u'.', u',', u'x', u'.', u'1', u'p', u'3', u'2', u'4', u'H', u'v', u'8', u'<']
    #markerfacecols = ['full', 'full', 'full', 'full', 'full', 'full', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r']
    #markersizes=[scalepng*1.5, scalepng*1.5, scalepng*1.5, 2., scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5,scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5]

    # sens SP / croix / B 
    #markers = [None, u'x', None, None,None, None,None, None,None, None,None, None,None, None,None, None,None, None,None, None,None, None]
    #markerfacecols = ['full', 'full', 'full', 'full', 'full', 'full', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r']
    #markersizes=[scalepng*1.5, scalepng*1.5, scalepng*1.5, 2., scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5,scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5]

    #markers = [None, u'x', u'x', None,None, None,None, None,None, None,None, None,None, None,None, None,None, None,None, None,None, None]
    #markerfacecols = ['full', 'full', 'full', 'full', 'full', 'full', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r']
    #markersizes=[scalepng*1.5, scalepng*1.5, scalepng*1.5, 2., scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5,scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5]


    # sens Sb / D / Ss / B
    #markers = [u'.', u'd', u'.', u'.', u'.', None, u'.', u'.', u'.', u',', u'x', u'.', u'1', u'p', u'3', u'2', u'4', u'H', u'v', u'8', u'<']
    #markerfacecols = ['full', 'full', 'full', 'white', 'full', 'full', 'full', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r']
    #markersizes=[4., 1.5, 2., 4., 2., scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5,scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5, scalepng*1.5]

    markerSize=scalepng*1.5
    
    if couleur=='standard': 
        colors = ['k', 'r', 'g', 'm', 'c', 'y', 'b','0.8','0.5']  
    elif couleur=='grey':
        colors = ['1.0','0.9','0.8','0.7','0.6','0.5','0.4','0.3','0.2','0.1']
    elif couleur=='mychoice1': 
        colors = ['r', 'c', 'm', 'g', 'y', 'k', 'b'] 
    elif couleur=='red':
        colors = ['r', 'r', 'r', 'r', 'r', 'r', 'r']
    else :
        colors=[]
        colori=['k', 'r', 'm', 'g', 'y', 'b','0.8','0.5'] # GB-HERE is default colors
        #colori=['m', 'b', 'c', 'g', 'y', 'r', 'k'] 
        #colori=['k', 'r', 'm', 'g', 'y', 'k', 'b'] 
        #colori=['k', 'm', 'b', 'c', 'g', 'y', 'r', 'b']
        #colori=['0.0', '0.35', '0.57', '0.8', 'r', 'g', 'm', 'c', 'y', 'b','0.25','0.5','0.75'] 
        #colori=['0.0', '0.45', '0.8', '0.65', '0.23', 'g', 'm', 'c', 'y', 'b','0.25','0.5','0.75'] 
        for k in range(len(couleur)):
            try:            
                    for i in range(couleur[k]):
                        colors.append(colori[k])
            except:
                colors.append('k')

    colors.append('k')
    styles = markers           
    axisNum = 0
    LineStyleNum=0
    ax = plt.subplot(111)
    
    if sym==[0]:
        sym_sca=[1]
        sym_vec=[1,1,1]
        sym_tens=[1,1,-1,1,1,-1,-1,-1,1]        
    elif sym==[1]:
        sym_sca=[1]
        sym_vec=[1,1,1]
        sym_tens=[1,1,1,1,1,1,1,1,1]  
    elif sym==[-1]:
        sym_sca=[-1]
        sym_vec=[-1,-1,-1]
        sym_tens=[-1,-1,-1,-1,-1,-1,-1,-1,-1]         
        
#     compt=0      
    for i in range(len(bfield)):

#         if sym!=60:
#             tracer([bfield[-1]], "graph_%d" % (compt), [0], sym=60)  
#             compt=compt+1     
        if sym==[10] or sym==60:
            pass
        else:
            if isinstance(bfield[i], tbx.ScalarField):
                bfield[i].symetriser(sym_sca,axe)
            elif isinstance(bfield[i], tbx.VectorField):
                bfield[i].symetriser(sym_vec,axe)    
            elif isinstance(bfield[i], tbx.Tensor2Field):
                bfield[i].symetriser(sym_tens,axe)  
                
                      
        if isinstance(bfield[i],tbx.Tensor2Field) or isinstance(bfield[i],tbx.Tensor3Field):
            bfield[i].preReshape()
        if zerocentre:
	    ### mettre zero au centre
            #if isinstance(bfield[i],tbx.VectorField):
            #    L=len(bfield[i]._npa[0,:])/2
            #    milieu_x=bfield[i]._npa[0,L]
            #    milieu_y=bfield[i]._npa[1,L]
            #    milieu_z=bfield[i]._npa[2,L]                        
            #    bfield[i]._npa[0,:]-=milieu_x
            #    bfield[i]._npa[1,:]-=milieu_y
            #    bfield[i]._npa[2,:]-=milieu_z
	    # mettre zero en y=0
            if isinstance(bfield[i],tbx.VectorField):
                milieu_x=bfield[i]._npa[0,0]
                milieu_y=bfield[i]._npa[1,0]
                milieu_z=bfield[i]._npa[2,0]                        
                bfield[i]._npa[0,:]-=milieu_x
                bfield[i]._npa[1,:]-=milieu_y
                bfield[i]._npa[2,:]-=milieu_z


        
        # j sont les composantes? 
        if axe != None:
            if log==True:
                    for j in axe:
                        color = colors[axisNum % len(colors)]
			if ((colors[axisNum-1]!=colors[axisNum] or axisNum==0) and doubleLegend) or doubleLegend==False:
				labell=bfield[i].tex
			else :
				labell=''
                        if  LineStyleNum< len(linestyles):
                            ax.semilogx(bfield[i]._ax[0,:]*xscale,bfield[i]._npa[j,:]*yscale, linestyles[LineStyleNum % len(linestyles)], label=labell, color=color, ms=markerSize ,markevery=markevry[i], markerfacecolor= color, markeredgecolor = color)                                    
                        else:
                            style = styles[(LineStyleNum- len(linestyles)) % len(styles)]
			    markerfacecol = markerfacecols[(LineStyleNum- len(linestyles)) % len(styles)]
			    if (markerfacecol=='full'):
				markerfacecol=color
                            markersizess=markersizes[(LineStyleNum- len(linestyles)) % len(styles)]
                            ax.semilogx(bfield[i]._ax[0,:]*xscale,bfield[i]._npa[j,:]*yscale, linestyle='-', linewidth=0.4, marker=style, label=labell, color=color, ms=markersizess, markevery=markevry[i], markerfacecolor=markerfacecol, markeredgecolor = color, markeredgewidth=0.2)                                    
                        
                        if colors[axisNum]!=colors[axisNum+1] and couleur!='standard':
                            axisNum += 1
                            LineStyleNum=0    
                        else:
                            axisNum += 1 
                            LineStyleNum+=1 
            else:
                    for j in axe:
                        color = colors[axisNum % len(colors)]
			if ((colors[axisNum-1]!=colors[axisNum] or axisNum==0) and doubleLegend) or doubleLegend==False:
				labell=bfield[i].tex
			else :
				labell=''
                        if  LineStyleNum< len(linestyles):
			    if (doubleLegend==True):
                            	ax.plot(bfield[i]._ax[0,:]*xscale,bfield[i]._npa[j,:]*yscale, linestyles[LineStyleNum % len(linestyles)], label='', color=color, ms=markerSize ,markevery=markevry[i], markerfacecolor= color, markeredgecolor = color) 
                            	plt.plot([],[],'-', color=color, label=labell)	
   
			    else:
	                        ax.plot(bfield[i]._ax[0,:]*xscale,bfield[i]._npa[j,:]*yscale, linestyles[LineStyleNum % len(linestyles)], label=labell, color=color, ms=markerSize ,markevery=markevry[i], markerfacecolor= color, markeredgecolor = color) 			                                
                        else:
                            style = styles[(LineStyleNum- len(linestyles)) % len(styles)]
			    markerfacecol = markerfacecols[(LineStyleNum- len(linestyles)) % len(styles)]
			    if (markerfacecol=='full'):
				markerfacecol=color
                            markersizess=markersizes[(LineStyleNum- len(linestyles)) % len(styles)]
			    if (doubleLegend==True):
                            	ax.plot(bfield[i]._ax[0,:]*xscale,bfield[i]._npa[j,:]*yscale, linestyle='-', linewidth=0.4, marker=style, label='', color=color, ms=markersizess, markevery=markevry[i], markerfacecolor=markerfacecol, markeredgecolor = color, markeredgewidth=0.2)  
                            	plt.plot([],[],'-', color=color, label=labell)	
                 	    else:                 
                                ax.plot(bfield[i]._ax[0,:]*xscale,bfield[i]._npa[j,:]*yscale, linestyle='-', linewidth=0.4, marker=style, label=labell, color=color, ms=markersizess, markevery=markevry[i], markerfacecolor=markerfacecol, markeredgecolor = color, markeredgewidth=0.2)             
                        if colors[axisNum]!=colors[axisNum+1] and couleur!='standard':
                            axisNum += 1
                            LineStyleNum=0    
                        else:
                            axisNum += 1 
                            LineStyleNum+=1 
                        #ax.plot(bfield[i]._ax[0,:]*xscale,bfield[i]._npa[j,:]*yscale, linestyles[axisNum % len(linestyles)], label=bfield[i].tex, color=color)                                    
        else:
             
            if log==True:      
                for j in range(len(bfield[i]._npa[:,0])):
		    if ((colors[axisNum-1]!=colors[axisNum] or axisNum==0) and doubleLegend) or doubleLegend==False:
			labell=bfield[i].tex
	            else :
		        labell=''
                    color = colors[axisNum % len(colors)]
                    plt.semilogx(bfield[i]._ax[0,:]*xscale,bfield[i]._npa[j,:]*yscale, linestyles[axisNum % len(linestyles)], label=labell, color=color)
                    axisNum += 1
            else:
                for j in range(len(bfield[i]._npa[:,0])):
		    if ((colors[axisNum-1]!=colors[axisNum] or axisNum==0) and doubleLegend) or doubleLegend==False:
			labell=bfield[i].tex
	            else :
		        labell=''
                    color = colors[axisNum % len(colors)]
                    if  LineStyleNum< len(linestyles):
                        plt.plot(bfield[i]._ax[0,:]*xscale,bfield[i]._npa[j,:]*yscale, linestyles[axisNum % len(linestyles)], label=labell, color=color, ms=markerSize, markevery=markevry[i], markerfacecolor='none', markeredgecolor = color)                
                    else:
                        style = styles[(LineStyleNum- len(linestyles)) % len(styles)]
                        markerfacecol = markerfacecols[(LineStyleNum- len(linestyles)) % len(styles)]
			if (markerfacecol=='full'):
				markerfacecol=color
                        markersizess=markersizes[(LineStyleNum- len(linestyles)) % len(styles)]
                        plt.plot(bfield[i]._ax[0,:]*xscale,bfield[i]._npa[j,:]*yscale, linestyle='-', linewidth=0.4, marker=style, label=labell, color=color, ms=markersizess, markevery=markevry[i], markerfacecolor=markerfacecol, markeredgecolor = color, markeredgewidth=0.2)                
                     
                    if colors[axisNum]!=colors[axisNum+1] and couleur!='standard':
                        axisNum += 1
                        LineStyleNum=0    
                    else:
                        axisNum += 1 
                        LineStyleNum+=1 
            
   
   
    if doubleLegend:
	if casename==None:
		print 'il faut remplir casename si option doublelegend est activee'
    	linestyles = ['-']
    	for i in range(len(casename)):
		color = colors[i]
		style = styles[i]
		markerfacecol = markerfacecols[i]
		if (markerfacecol=='full'):
			markerfacecol=color
		markersizess=markersizes[i]
		if i==0:
	    		plt.plot([],[],linestyles[i], label=casename[i], color='k', linewidth=0.4, marker=style, ms=markersizess, markerfacecolor=markerfacecol, markeredgecolor = 'k', markeredgewidth=1.0)			
	    	else:	
			markerlocal=markers[i-1]
			if (markers[i-1]==None):
				markerlocal=u'.'
			plt.plot([],[],markerlocal, label=casename[i], color='k', marker=style, ms=markersizess, markerfacecolor=markerfacecol, markeredgecolor = 'k', markeredgewidth=1.0)




    if grid:
        plt.grid('on')
    if xlim != None:
        ax.set_xlim(xlim, emit=False)
        plt.xlim(xlim)        
    if ylim != None:
        plt.ylim(ylim)
    else:
        try :
            autoscale_y(ax)
        except:
            print 'auto scale fail'
    #plt.xlabel(r'$y/h$')
    ax.set_xlabel(r'$\mathbf{y/h}$', fontsize=8)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left') 
    plt.ticklabel_format(style='plain', axis='y', scilimits=(0,0))
    #plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.tight_layout()
    dpii=400
    #plt.title(textitle)
    if applatir:
        f.set_size_inches(scalepng,scalepng/2.)
    else:
        f.set_size_inches(scalepng,scalepng)      
    inSizeLegend=int(scalepng*3.5)
    sauv=title+'_off.png'
    if plotall:
        plt.savefig(sauv, bbox_inches='tight', dpi=dpii) 
        #plt.legend(loc=(0.2,0.7),prop={'size':inSizeLegend}, frameon=True, ncol=2)
	plt.legend(loc=0,prop={'size':inSizeLegend}, frameon=True, ncol=1)
        sauv=title+'_in.png'
        plt.savefig(sauv, bbox_inches='tight', dpi=dpii)    
      
    plt.legend(loc = 0)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1.1, 0.5), frameon=True)
    sauv=title+'_out.png'
    plt.savefig(sauv, bbox_inches='tight', dpi=dpii)
    
    
    for i in range(len(bfield)):
        if isinstance(bfield[i],tbx.Tensor2Field) or isinstance(bfield[i],tbx.Tensor3Field):
            bfield[i].postReshape()
    plt.close(f)

# Methode importee de Tools_fixedBubbles : 
def tracer_fixedBubbles(Run, bfield, title, axe=None, textitle='', couleur=None):
    """ 
    Creation of a parametrized plotinitialization of a field object
    @param the field's list to plot (can be of any sub-field type)
    @param the name of the png file created
    @param can contain the list of components to be plotted. Default is None for all   
    @param the title of the plot (optional)
    @param limit of y axis. Default is None for all 
    @param limit of x axis. Default is None for all 
    @param Plot the grid background. Default is True
    @param legend ... pas au point
    @param markerSize. Default is an autoscale from the scale of the png file (scalepng).
    @param Plot with logarithmic scale. Default is False
    @param couleur=[3,2,3,1] means the first 3 variables are plotted in the same color, then the following two with an other etc..
    @param markevry
    @param dimenssionless tools for y axis
    @param dimenssionless tools for x axis
    @param scale of the out.png
    @param aspect ratio of the picture = 0.5. Default is False 

    """
    if (couleur==None):
	couleur=Run["plot"]["couleur"]

    xlim=Run["plot"]["xlim"]  
    ylim= Run["plot"]["ylim"]  
    grid=Run["plot"]["grid"]
    legend=Run["plot"]["legend"]
    markerSize=Run["plot"]["markerSize"]
    log=Run["plot"]["log"]
    markevry=Run["plot"]["markevry"]
    yscale=Run["plot"]["yscale"]
    xscale=Run["plot"]["xscale"]
    scalepng=Run["plot"]["scalepng"]
    applatir=Run["plot"]["applatir"]
    plotall=Run["plot"]["plotall"]
    sym=Run["plot"]["sym"]
    zerocentre=Run["plot"]["zerocentre"]
    doubleLegend=Run["plot"]["doubleLegend"]
    casename=Run["plot"]["casename"]

 
    if len(markevry)==1 :
	cte=markevry[0]
	tmp=[cte,cte,cte,cte,cte,cte,cte,cte,cte,cte,cte,cte,cte,cte,cte,cte,cte,cte,cte,cte,cte]
	markevry=tmp

    global compteur
    compteur=compteur+1
    f=plt.figure(compteur)
    #f.tightlayout()
    #linestyles = ['_', '-', '--', ':']
    

    linestyles = []
    markers = [u'.', u'x', u'+', u'_', u'd', u'h', u'+', u'*', u',', u'x', u'.', u'1', u'p', u'3', u'2', u'4', u'H', u'v', u'8', u'<', u'>']
    markerSize=scalepng*1.5
    
    if couleur=='standard': 
        colors = ['k', 'r', 'g', 'm', 'c', 'y', 'b']  
    elif couleur=='grey':
        colors = ['1.0','0.9','0.8','0.7','0.6','0.5','0.4','0.3','0.2','0.1']
    elif couleur=='mychoice1': 
        colors = ['r', 'c', 'm', 'g', 'y', 'k', 'b'] 
    elif couleur=='red':
        colors = ['r', 'r', 'r', 'r', 'r', 'r', 'r']
    else :
        colors=[]
        colori=['k', 'r', 'g', 'g', 'y', 'k', 'b'] 
        #colori=['0.0', '0.35', '0.57', '0.8', 'r', 'g', 'm', 'c', 'y', 'b','0.25','0.5','0.75'] 
        #colori=['0.0', '0.45', '0.8', '0.65', '0.23', 'g', 'm', 'c', 'y', 'b','0.25','0.5','0.75'] 
        for k in range(len(couleur)):
            try:            
                    for i in range(couleur[k]):
                        colors.append(colori[k])
            except:
                colors.append('k')

    colors.append('k')
    styles = markers + [
    r'$\lambda$',
    r'$\bowtie$',
    r'$\circlearrowleft$',
    r'$\clubsuit$',
    r'$\checkmark$']            
    axisNum = 0
    LineStyleNum=0
    ax = plt.subplot(111)
    
    if sym==[0]:
        sym_sca=[1]
        sym_vec=[1,1,1]
        sym_tens=[1,1,-1,1,1,-1,-1,-1,1]        
    elif sym==[1]:
        sym_sca=[1]
        sym_vec=[1,1,1]
        sym_tens=[1,1,1,1,1,1,1,1,1]  
    elif sym==[-1]:
        sym_sca=[-1]
        sym_vec=[-1,-1,-1]
        sym_tens=[-1,-1,-1,-1,-1,-1,-1,-1,-1]         
        
#     compt=0      
    for i in range(len(bfield)):

#         if sym!=60:
#             tracer([bfield[-1]], "graph_%d" % (compt), [0], sym=60)  
#             compt=compt+1     
        if sym==[10] or sym==60:
            pass
        else:
            if isinstance(bfield[i], tbx.ScalarField):
                bfield[i].symetriser(sym_sca,axe)
            elif isinstance(bfield[i], tbx.VectorField):
                bfield[i].symetriser(sym_vec,axe)    
            elif isinstance(bfield[i], tbx.Tensor2Field):
                bfield[i].symetriser(sym_tens,axe)  
                
                      
        if isinstance(bfield[i],tbx.Tensor2Field) or isinstance(bfield[i],tbx.Tensor3Field):
            bfield[i].preReshape()
        if zerocentre:
            if isinstance(bfield[i],tbx.VectorField):
                L=len(bfield[i]._npa[0,:])/2
                milieu_x=bfield[i]._npa[0,L]
                milieu_y=bfield[i]._npa[1,L]
                milieu_z=bfield[i]._npa[2,L]                        
                bfield[i]._npa[0,:]-=milieu_x
                bfield[i]._npa[1,:]-=milieu_y
                bfield[i]._npa[2,:]-=milieu_z
        
        # j sont les composantes? 
        if axe != None:
            if log==True:
                    for j in axe:
                        color = colors[axisNum % len(colors)]
			if ((colors[axisNum-1]!=colors[axisNum] or axisNum==0) and doubleLegend) or doubleLegend==False:
				labell=bfield[i].tex
			else :
				labell=''
                        if  LineStyleNum< len(linestyles):
                            ax.semilogx(bfield[i]._ax[0,:]*xscale,bfield[i]._npa[j,:]*yscale, linestyles[LineStyleNum % len(linestyles)], label=labell, color=color, ms=markerSize ,markevery=markevry[i], markerfacecolor= color, markeredgecolor = color)                                    
                        else:
                            style = styles[(LineStyleNum- len(linestyles)) % len(styles)]
                            ax.semilogx(bfield[i]._ax[0,:]*xscale,bfield[i]._npa[j,:]*yscale, '--', marker=style, label=labell, color=color, ms=markerSize ,markevery=markevry[i], markerfacecolor= color, markeredgecolor = color)                                    
                        
                        if colors[axisNum]!=colors[axisNum+1] and couleur!='standard':
                            axisNum += 1
                            LineStyleNum=0    
                        else:
                            axisNum += 1 
                            LineStyleNum+=1 
            else:
                    for j in axe:
                        color = colors[axisNum % len(colors)]
			if ((colors[axisNum-1]!=colors[axisNum] or axisNum==0) and doubleLegend) or doubleLegend==False:
				labell=bfield[i].tex
			else :
				labell=''
                        if  LineStyleNum< len(linestyles):
                            ax.plot(bfield[i]._ax[0,:]*xscale,bfield[i]._npa[j,:]*yscale, linestyles[LineStyleNum % len(linestyles)], label=labell, color=color, ms=markerSize ,markevery=markevry[i], markerfacecolor= color, markeredgecolor = color)                                    
                        else:
                            style = styles[(LineStyleNum- len(linestyles)) % len(styles)]
                            ax.plot(bfield[i]._ax[0,:]*xscale,bfield[i]._npa[j,:]*yscale, '--', marker=style, label=labell, color=color, ms=markerSize ,markevery=markevry[i], markerfacecolor= color, markeredgecolor = color)                                    
                        
                        if colors[axisNum]!=colors[axisNum+1] and couleur!='standard':
                            axisNum += 1
                            LineStyleNum=0    
                        else:
                            axisNum += 1 
                            LineStyleNum+=1 
                        #ax.plot(bfield[i]._ax[0,:]*xscale,bfield[i]._npa[j,:]*yscale, linestyles[axisNum % len(linestyles)], label=bfield[i].tex, color=color)                                    
        else:
             
            if log==True:      
                for j in range(len(bfield[i]._npa[:,0])):
		    if ((colors[axisNum-1]!=colors[axisNum] or axisNum==0) and doubleLegend) or doubleLegend==False:
			labell=bfield[i].tex
	            else :
		        labell=''
                    color = colors[axisNum % len(colors)]
                    plt.semilogx(bfield[i]._ax[0,:]*xscale,bfield[i]._npa[j,:]*yscale, linestyles[axisNum % len(linestyles)], label=labell, color=color)
                    axisNum += 1
            else:
                #import pdb; pdb.set_trace()
		for j in range(len(bfield[i]._npa[:,0])):
		    if ((colors[axisNum-1]!=colors[axisNum] or axisNum==0) and doubleLegend) or doubleLegend==False:
			labell=bfield[i].tex
	            else :
		        labell=''
                    color = colors[axisNum % len(colors)]
                    if  LineStyleNum< len(linestyles):
                        plt.plot(bfield[i]._ax[0,:]*xscale,bfield[i]._npa[j,:]*yscale, linestyles[axisNum % len(linestyles)], label=labell, color=color, ms=markerSize, markevery=markevry[i], markerfacecolor= color, markeredgecolor = color)                
                    else:
                        style = styles[(LineStyleNum- len(linestyles)) % len(styles)]
                        plt.plot(bfield[i]._ax[0,:]*xscale,bfield[i]._npa[j,:]*yscale, linestyle=':', linewidth=0.2, marker=style, label=labell, color=color, ms=markerSize, markevery=markevry[i], markerfacecolor= color, markeredgecolor = color)                
                     
                    if colors[axisNum]!=colors[axisNum+1] and couleur!='standard':
                        axisNum += 1
                        LineStyleNum=0    
                    else:
                        axisNum += 1 
                        LineStyleNum+=1 
            
   
   
    if doubleLegend:
	if casename==None:
		print 'il faut remplir casename si option doublelegend est activee'
    	linestyles = ['-']
    	markers = [u'.', u'x', u'+', u'_', u'd', u'h', u'+', u'*', u',', u'x', u'.', u'1', u'p', u'3', u'2', u'4', u'H', u'v', u'8', u'<', u'>']
	for i in range(len(casename)):
		if i==0:
	    		plt.plot([],[],linestyles[i], label=casename[i], color='k', ms=markerSize)			
	    	else:	
			plt.plot([],[],markers[i-1], label=casename[i], color='k', ms=markerSize)




    if grid:
        plt.grid('on')
    if xlim != None:
        ax.set_xlim(xlim, emit=False)
        plt.xlim(xlim)        
    if ylim != None:
        plt.ylim(ylim)
    else:
        try :
            autoscale_y(ax)
        except:
            print 'auto scale fail'
    plt.xlabel(r'$y/h$')

    ax.spines['top'].set_visible(True)
    ax.spines['right'].set_visible(True)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.tight_layout()
    ## resolution de l'image
    dpii=400
    #plt.title(textitle)
    if applatir:
        f.set_size_inches(scalepng,scalepng/2.)
    else:
        f.set_size_inches(scalepng,scalepng)      
    inSizeLegend=int(scalepng*2.0)
    sauv=Run["plot"]["name"]+"_"+title+'_off.png'
    if plotall:
        print 'exporting graph '+sauv
        plt.savefig(sauv, bbox_inches='tight', dpi=dpii) 
        plt.legend(loc=0,prop={'size':inSizeLegend}, frameon=True)
        sauv=Run["plot"]["name"]+"_"+title+'_in.png'
        print 'exporting graph '+sauv
        plt.savefig(sauv, bbox_inches='tight', dpi=dpii)    
      
    plt.legend(loc = 0)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1.1, 0.5), frameon=True)
    sauv=Run["plot"]["name"]+"_"+title+'_out.png'
    print 'exporting graph '+sauv
    plt.savefig(sauv, bbox_inches='tight', dpi=dpii)
    
    
    for i in range(len(bfield)):
        if isinstance(bfield[i],tbx.Tensor2Field) or isinstance(bfield[i],tbx.Tensor3Field):
            bfield[i].postReshape()
    plt.close(f)


# Methode importee de Tools_fixedBubbles : 
def tracerJFM_fixedBubbles(field, name="blabla"):

	color2="#8dd3c7"
	color4="#ffffb3"
	color1="#bebada"
	color3="#fb8072"

        ###### figure 1
	scalepng=3.
	f, (ax) = plt.subplots()
	markeredgewidths=[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.]


	### model production
        if ("production_WIT" in name):
		style=['-', '-', '-', '-', '-', '-', '-']
		dash=[(0.001,1000.),(0.001,1000.),(0.001,1000.),(0.001,1000.),(1000.,1.),(1.,1.), (2.5,2.5)]
		color=['g', 'b', 'k', 'r', 'r', 'b', 'g']
		marker=[u'*', u'.', u'x', u'd', None, None, None, None]
		label=[r"$\Pi^{WIF}_{fixed}$",r"$\Pi^{WIT}_{fixed}$",r"$\Pi_{fixed}$",r"$\Pi_{free}$",r"$\alpha_v\Delta\rho gu_r/\rho_l$", r"$\Pi_{fixed}[0.9-e^{-0.006Re_b}]$", r"$\Pi_{fixed}[0.1+e^{-0.006Re_b}]$"]
		markerfacecols=['none','none','none','none','none','none','none']
		markersizes=[scalepng*3, scalepng*3, scalepng*1.5, scalepng*2., scalepng*3, scalepng*3, scalepng*3]
		markeredgewidths=[0.5,0.5,0.5,0.5,0.5,0.5,0.5]
		ax.set_xlim([0.,900], emit=False)
		#ax.set_ylim([0,0.3], emit=False)
		ax.set_ylim([0,0.15], emit=False)
		ax.set_xlabel(r"${Re_b}$", fontsize=15)
		ax.set_ylabel(r"${\Pi_{11}}$", fontsize=15)


        if ("correlation_triple" in name):
		style=['-', '-', '-', '-', '-', '-', '-']
		dash=[(0.001,1000.),(0.001,1000.),(0.001,1000.),(0.001,1000.),(1000.,1.),(1.,1.), (2.5,2.5)]
		color=['g', 'b', 'k', 'r', 'r', 'b', 'g']
		marker=[u'*', u'.', u'x', u'd', None, None, None, None]
		label=[r"$\overline{\chi_lu_l^*u_l^*u_l^*}$",r"$3{\overline{\overline{\chi_lu_l^*}^T\overline{\chi_lu_l^{*'}u_l^{*'}}^T}}$"]
		markerfacecols=['none','none','none','none','none','none','none']
		markersizes=[scalepng*3, scalepng*3, scalepng*1.5, scalepng*2., scalepng*3, scalepng*3, scalepng*3]
		markeredgewidths=[0.5,0.5,0.5,0.5,0.5,0.5,0.5]
		ax.set_xlim([0.,900], emit=False)
		ax.set_ylim([0.,5.e-3], emit=False)
		ax.set_xlabel(r"$Re_b$", fontsize=15)
		ax.set_ylabel(r"$m^3.s^{-3}$", fontsize=15)


	### model epsilon WIT
	if ("lambda_epsilon" in name):
		style=['-', '-', '-', '-', '-', '-', '-']
		dash=[(0.001,1000.),(0.001,1000.),(0.001,1000.),(0.001,1000.),(1000.,1.),(1.,1.), (2.5,2.5)]
		color=['k', 'r', 'r', 'r', 'r', 'b', 'g']
		marker=[u'*', u'.', u'x', u'd', None, None, None, None]
		label=[r"$\Lambda/20$ (Amoura2017)",r"$\sqrt{\nu R_{11}^{WIT}/\epsilon_{11}^{WIT}}/1.7$",r"$\sqrt{\nu R_{22}^{WIT}/\epsilon_{22}^{WIT}}/1.7$",r"$\sqrt{\nu R_{33}^{WIT}/\epsilon_{33}^{WIT}}/1.7$",r"$1/\sqrt{C_dRe_b}$"]
		markerfacecols=['none','none','none','none','none','none','none']
		markersizes=[scalepng*3, scalepng*3, scalepng*1.5, scalepng*2., scalepng*3, scalepng*3, scalepng*3]
		markeredgewidths=[0.5,0.5,0.5,0.5,0.5,0.5,0.5]
		ax.set_xlim([0.,1100.], emit=False)
		ax.set_ylim([0.0,0.15], emit=False)
		ax.set_xlabel(r"${Re_b}$", fontsize=15)
		ax.set_ylabel(r"${\Lambda/d_b}$", fontsize=15)


	### model epsilon WIF
	#style=['-', '-', '-', '-', '-', '-', '-']
	#dash=[(0.001,1000.),(0.001,1000.),(0.001,1000.),(1000.,1.),(1.,1.), (2.5,2.5)]
	#color=['k', 'k', 'k', 'r', 'b', 'g']
        #marker=[u'*', u'.', u'x', None, None, None, None]
	#label=[r"$\sqrt{\nu R_{11}^{WIF}/\epsilon_{11}^{WIF}}$",r"$\sqrt{\nu R_{22}^{WIF}/\epsilon_{22}^{WIF}}$",r"$\sqrt{\nu R_{33}^{WIF}/\epsilon_{33}^{WIF}}$",r"$[0.3+0.2e^{-0.006(Re_b-100)}]$", r"$0.2$"]
	#markerfacecols=['none','none','none','none','none','none','none']
        #markersizes=[scalepng*3, scalepng*3, scalepng*1.5, scalepng*2., scalepng*3, scalepng*3, scalepng*3]
        #markeredgewidths=[0.5,0.5,0.5,0.5,0.5,0.5,0.5]
        #ax.set_xlim([0.,1000], emit=False)
        #ax.set_ylim([0,1.], emit=False)
	#ax.set_xlabel(r"${Re_b}$", fontsize=15)
	#ax.set_ylabel(r"${\Lambda/d_b}$", fontsize=15)

	### redistribution
	#style=['-', '-', '-', '-', '-', '-', '-']
	#dash=[(0.001,1000.),(0.001,1000.),(0.001,1000.),(1000.,1.),(1.,1.), (2.5,2.5)]
	#color=['k', 'k', 'k', 'r', 'b', 'g']
        #marker=[u'*', u'.', u'x', None, None, None, None]
	#label=[r"$(\phi_{11}^{WIF}+\Pi_{11}^{WIF})/\Pi_{11}^{WIF}$",r"$(\phi_{22}^{WIF}+\Pi_{22}^{WIF})/\Pi_{11}^{WIF}$",r"$(\phi_{33}^{WIF}+\Pi_{33}^{WIF})/\Pi_{11}^{WIF}$",r"$1/2$", r"$1/4$"]
	#markerfacecols=['none','none','none','none','none','none','none']
        #markersizes=[scalepng*3, scalepng*3, scalepng*1.5, scalepng*2., scalepng*3, scalepng*3, scalepng*3]
        #markeredgewidths=[0.5,0.5,0.5,0.5,0.5,0.5,0.5]
        #ax.set_xlim([0.,1000], emit=False)
        #ax.set_ylim([0,1.], emit=False)
	#ax.set_xlabel(r"${Re_b}$", fontsize=15)
	#ax.set_ylabel(r"${\Lambda/d_b}$", fontsize=15)

	### model complet

	if ("Fixe_" in name):
		style=['-', '-', '-', '-', '-', '-', '-']
		dash=[(0.001,1000.),(1000.,1.),(0.001,1000.),(1000.,1.),(1.,1.), (2.5,2.5),(5.,5.)]
		color=['k', 'k', 'r', 'r', 'k', 'k', 'k']
		marker=[u'd', None, u'd', None, None, None, None, None]
		label=[r"DNS $R_{11}$",r"model $R_{11}$",r"DNS $R_{22}$",r"model $R_{22}$",r"$R^{WIT}$ : $C_\Lambda=1.6$", r"$R^{WIF}_{11}$ : $C_V=1.8$", r"$R^{WIF}_{22}$"]
		markerfacecols=['none','none','none','none','none','none','none']
		markersizes=[scalepng*2., scalepng*3, scalepng*2., scalepng*3., scalepng*3, scalepng*3, scalepng*3]
		markeredgewidths=[0.5,0.5,0.5,0.5,0.5,0.5,0.5]
		ax.set_xlim([0.,900], emit=False)
		ax.set_ylim([0,0.004], emit=False)
		ax.set_xlabel(r"${Re_b}$", fontsize=15)
		ax.set_ylabel(r"$m^2.s^{-2}$", fontsize=15)


	if ("Libre_" in name):
		style=['-', '-', '-', '-', '-', '-', '-']
		dash=[(0.001,1000.),(1000.,1.),(0.001,1000.),(1000.,1.),(1.,1.), (2.5,2.5),(5.,5.)]
		color=['k', 'k', 'r', 'r', 'k', 'k', 'k']
		marker=[u'd', None, u'd', None, None, None, None, None]
		if ("cd_reel" in name):
			label=[r"$R_{11}$",r"model $C_V=0.36$",r"$R_{22}$",r"model $C_\Lambda=2.7$"]
		else :
			label=[r"$R_{11}$",r"model $C_d^{sphe}$, $C_\Lambda=2.7$, $C_V=0.25\gamma^3$",r"$R_{22}$",r"model $C_d^{sphe}$, $C_\Lambda=2.7$, $C_V=0.25\gamma^3$"]

		markerfacecols=['none','none','none','none','none','none','none']
		markersizes=[scalepng*2., scalepng*3, scalepng*2., scalepng*3., scalepng*3, scalepng*3, scalepng*3]
		markeredgewidths=[0.5,0.5,0.5,0.5,0.5,0.5,0.5]
		ax.set_xlim([0.,900], emit=False)
		ax.set_ylim([0,0.002], emit=False)
		ax.set_xlabel(r"${Re_b}$", fontsize=15)
		ax.set_ylabel(r"$m^2.s^{-2}$", fontsize=15)

	if ("Libre_WIT" in name):
		style=['-', '-', '-', '-', '-', '-', '-']
		dash=[(0.001,1000.),(1000.,1.),(0.001,1000.),(1000.,1.),(1.,1.), (2.5,2.5),(5.,5.)]
		color=['k', 'k', 'r', 'r', 'k', 'k', 'k']
		marker=[u'd', None, u'd', None, None, None, None, None]
		if ("cd_reel" in name):
			label=[r"$R_{22}$",r"model $C_\Lambda=2.7$"]
		else :
			label=[r"$R_{22}$",r"model $C_\Lambda=2.7$"]

		markerfacecols=['none','none','none','none','none','none','none']
		markersizes=[scalepng*2., scalepng*3, scalepng*2., scalepng*3., scalepng*3, scalepng*3, scalepng*3]
		markeredgewidths=[0.5,0.5,0.5,0.5,0.5,0.5,0.5]
		ax.set_xlim([0.,900], emit=False)
		ax.set_ylim([0,0.002], emit=False)
		ax.set_xlabel(r"${Re_b}$", fontsize=15)
		ax.set_ylabel(r"$m^2.s^{-2}$", fontsize=15)


	if ("Libre_WIF" in name):
		style=['-', '-', '-', '-', '-', '-', '-']
		dash=[(0.001,1000.),(1000.,1.),(0.001,1000.),(1000.,1.),(1.,1.), (2.5,2.5),(5.,5.)]
		color=['k', 'k', 'r', 'r', 'k', 'k', 'k']
		marker=[u'd', None, u'd', None, None, None, None, None]
		if ("cd_reel" in name):
			label=[r"$R_{11}$",r"model $C_V=0.25\gamma^3$"]
		else :
			label=[r"$R_{11}$",r"model $C_V=0.25\gamma^3$"]

		markerfacecols=['none','none','none','none','none','none','none']
		markersizes=[scalepng*2., scalepng*3, scalepng*2., scalepng*3., scalepng*3, scalepng*3, scalepng*3]
		markeredgewidths=[0.5,0.5,0.5,0.5,0.5,0.5,0.5]
		ax.set_xlim([0.,900], emit=False)
		ax.set_ylim([0,0.002], emit=False)
		ax.set_xlabel(r"${Re_b}$", fontsize=15)
		ax.set_ylabel(r"$m^2.s^{-2}$", fontsize=15)

	### compa libre fixe
	#style=['-', '-', '-', '-', '-', '-', '-']
	#dash=[(1000.,0.0001),(1.,1.),(1000.,0.0001),(1.,1.),(1000.,0.0001),(1.,1.)]
	#color=['k', 'r', 'k', 'r', 'k', 'r', 'g']
        #marker=[None, None, u'd', u'd', u'*', u'*', None, None]
	#label=[r"$R_{11}$ fixed",r"$R_{11}$ free",r"$R_{22}$ fixed",r"$R_{22}$ free", r"$R_{33}$ fixed",r"$R_{33}$ free"]
	#markerfacecols=['none','none','none','none','none','none','none']
        #markersizes=[scalepng*3, scalepng*3, scalepng*2., scalepng*2., scalepng*2., scalepng*2., scalepng*3]
        #markeredgewidths=[0.5,0.5,0.5,0.5,0.5,0.5,0.5]
        #ax.set_xlim([0.,1000.], emit=False)
        #ax.set_ylim([0.0,0.3], emit=False)
	#ax.set_xlabel(r"${Re_b}$", fontsize=15)
	#ax.set_ylabel(r"${R_{ij}}/u_r^2$", fontsize=15)

	### compa riboux donnees / model
	if ("Riboux" in name):
		style=['-', '-', '-', '-', '-', '-', '-']
		dash=[(0.0001,1000.),(1000.,0.0001),(0.0001,1000.),(1000.,0.0001),(1000.,0.0001),(1.,1.)]
		color=['k', 'k','r', 'r', 'k', 'r', 'g']
		marker=[u'd', None, u'd', None, u'*', u'*', None, None]
		if ("cd_bulle_spherique_riboux" in name):
			if ("0.0025" in name):
				label=[r"$R_{11}$ exp $d_b=0.0025$",r"model $C_d^{sphe}$, $C_\Lambda=2.7$, $C_V=0.12\gamma^3$", r"$R_{22}$ exp $d_b=0.0025$",r"model $C_d^{sphe}$, $C_\Lambda=2.7$, $C_V=0.12\gamma^3$"]
			if ("0.0016" in name):
				label=[r"$R_{11}$ exp $d_b=0.0016$",r"model $C_d^{sphe}$, $C_\Lambda=2.7$, $C_V=0.12\gamma^3$", r"$R_{22}$ exp $d_b=0.0016$",r"model $C_d^{sphe}$, $C_\Lambda=2.7$, $C_V=0.12\gamma^3$"]
			if ("0.0021" in name):
				label=[r"$R_{11}$ exp $d_b=0.0021$",r"model $C_d^{sphe}$, $C_\Lambda=2.7$, $C_V=0.12\gamma^3$", r"$R_{22}$ exp $d_b=0.0021$",r"model $C_d^{sphe}$, $C_\Lambda=2.7$, $C_V=0.12\gamma^3$"]
		elif ("cd_reel" in name):
			if ("0.0025" in name):
				label=[r"$R_{11}$ exp $d_b=0.0025$",r"model $C_d^{real}$, $C_\Lambda=3.8$, $C_V=0.12\gamma^3$", r"$R_{22}$ exp $d_b=0.0025$",r"model $C_d^{real}$, $C_\Lambda=3.8$, $C_V=0.12\gamma^3$"]
			if ("0.0016" in name):
				label=[r"$R_{11}$ exp $d_b=0.0016$",r"model $C_d^{real}$, $C_\Lambda=2.6$, $C_V=0.12\gamma^3$", r"$R_{22}$ exp $d_b=0.0016$",r"model $C_d^{real}$, $C_\Lambda=2.6$, $C_V=0.12\gamma^3$"]
			if ("0.0021" in name):
				label=[r"$R_{11}$ exp $d_b=0.0021$",r"model $C_d^{real}$, $C_\Lambda=3.5$, $C_V=0.12\gamma^3$", r"$R_{22}$ exp $d_b=0.0021$",r"model $C_d^{real}$, $C_\Lambda=3.5$, $C_V=0.12\gamma^3$"]
		else :
			if ("0.0025" in name):
				label=[r"$R_{11}$ exp $d_b=0.0025$",r"model", r"$R_{22}$ exp $d_b=0.0025$",r"model"]
			if ("0.0016" in name):
				label=[r"$R_{11}$ exp $d_b=0.0016$",r"model", r"$R_{22}$ exp $d_b=0.0025$",r"model"]
			if ("0.0021" in name):
				label=[r"$R_{11}$ exp $d_b=0.0021$",r"model", r"$R_{22}$ exp $d_b=0.0025$",r"model"]
	
		markerfacecols=['none','none','none','none','none','none','none']
		markersizes=[scalepng*2, scalepng*2, scalepng*2., scalepng*2., scalepng*2., scalepng*2., scalepng*3]
		markeredgewidths=[0.5,0.5,0.5,0.5,0.5,0.5,0.5]
		ax.set_xlim([0.,0.105], emit=False)
		ax.set_ylim([0.0,0.022], emit=False)
		ax.set_xlabel(r"${\alpha_v}$", fontsize=15)
		ax.set_ylabel(r"$R_{ij}$", fontsize=15)

	if ("Riboux_WIF" in name):
		style=['-', '-', '-', '-', '-', '-', '-']
		dash=[(0.0001,1000.),(1000.,0.0001),(0.0001,1000.),(1000.,0.0001),(1000.,0.0001),(1.,1.)]
		color=['k', 'k','r', 'r', 'k', 'r', 'g']
		marker=[u'd', None, u'd', None, u'*', u'*', None, None]
		if ("cd_bulle_spherique_riboux" in name):
			if ("0.0025" in name):
				label=[r"$R_{11}$",r"model $C_V=0.12\gamma^3$"]
			if ("0.0016" in name):
				label=[r"$R_{11}$",r"model $C_V=0.12\gamma^3$"]
			if ("0.0021" in name):
				label=[r"$R_{11}$",r"model $C_V=0.12\gamma^3$"]
		elif ("cd_reel" in name):
			if ("0.0025" in name):
				label=[r"$R_{11}$",r"model $C_V=0.12\gamma^3$"]
			if ("0.0016" in name):
				label=[r"$R_{11}$",r"model $C_V=0.12\gamma^3$"]
			if ("0.0021" in name):
				label=[r"$R_{11}$",r"model $C_V=0.12\gamma^3$"]
		else :
			if ("0.0025" in name):
				label=[r"$R_{11}$",r"model"]
			if ("0.0016" in name):
				label=[r"$R_{11}$",r"model"]
			if ("0.0021" in name):
				label=[r"$R_{11}$",r"model"]
	
		markerfacecols=['none','none','none','none','none','none','none']
		markersizes=[scalepng*2, scalepng*2, scalepng*2., scalepng*2., scalepng*2., scalepng*2., scalepng*3]
		markeredgewidths=[0.5,0.5,0.5,0.5,0.5,0.5,0.5]
		ax.set_xlim([0.,0.105], emit=False)
		ax.set_ylim([0.0,0.022], emit=False)
		ax.set_xlabel(r"${\alpha_v}$", fontsize=15)
		ax.set_ylabel(r"$m^2.s^{-2}$", fontsize=15)

	if ("Riboux_WIT" in name):
		style=['-', '-', '-', '-', '-', '-', '-']
		dash=[(0.0001,1000.),(1000.,0.0001),(0.0001,1000.),(1000.,0.0001),(1000.,0.0001),(1.,1.)]
		color=['k', 'k','r', 'r', 'k', 'r', 'g']
		marker=[u'd', None, u'd', None, u'*', u'*', None, None]
		if ("cd_bulle_spherique_riboux" in name):
			if ("0.0025" in name):
				label=[r"$R_{22}$",r"model $C_\Lambda=2.7$"]
			if ("0.0016" in name):
				label=[r"$R_{22}$",r"model $C_\Lambda=2.7$"]
			if ("0.0021" in name):
				label=[r"$R_{22}$",r"model $C_\Lambda=2.7$"]
		elif ("cd_reel" in name):
			if ("0.0025" in name):
				label=[r"$R_{22}$",r"model $C_\Lambda=3.8$"]
			if ("0.0016" in name):
				label=[r"$R_{22}$",r"model $C_\Lambda=2.6$"]
			if ("0.0021" in name):
				label=[r"$R_{22}$",r"model $C_\Lambda=3.5$"]
		else :
			if ("0.0025" in name):
				label=[r"$R_{22}$",r"model"]
			if ("0.0016" in name):
				label=[r"$R_{22}$",r"model"]
			if ("0.0021" in name):
				label=[r"$R_{22}$",r"model"]
	
		markerfacecols=['none','none','none','none','none','none','none']
		markersizes=[scalepng*2, scalepng*2, scalepng*2., scalepng*2., scalepng*2., scalepng*2., scalepng*3]
		markeredgewidths=[0.5,0.5,0.5,0.5,0.5,0.5,0.5]
		ax.set_xlim([0.,0.105], emit=False)
		ax.set_ylim([0.0,0.022], emit=False)
		ax.set_xlabel(r"${\alpha_v}$", fontsize=15)
		ax.set_ylabel(r"$m^2.s^{-2}$", fontsize=15)

	if ('cd_compa' in name):
		style=['-', '-', '-', '-', '-', '-', '-', '-', '-']
		dash=[(0.001,1000),(0.001,1000), (1000,0.001),(1000,0.001),(1000,0.001),(3.,3.),(3.,3.),(3.,3.),(3.,3.)]
		color=['k', 'k','r', 'b', 'g', 'r','b','g', 'r', 'y']
		marker=[u'd', u'*', u'x', u'+', u'.', u'+', u'x', u'+']
		label=[r"$C_d$", r"$C_d^*$", r"Moore (1963)",r"Mei (1994)",r"Schiller (1933)",r"Ishii and Zuber (1979)", r"Tomiyama (1998)"]
		markerfacecols=['none','none','none','none','none','none','none','none','none']
		markersizes=[scalepng*2, scalepng*3, scalepng*2., scalepng*2., scalepng*2., scalepng*2., scalepng*2, scalepng*2, scalepng*2]
		markeredgewidths=[0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5]
		ax.set_xlim([0.,1100.], emit=False)
		ax.set_ylim([0.0,2.2], emit=False)
		ax.set_xlabel(r"${Re_b}$", fontsize=15)
		ax.set_ylabel(r"$C_d$", fontsize=15)

	if ('cd_compa_exp' in name):
		style=['-', '-', '-', '-', '-', '-', '-', '-', '-']
		dash=[(0.001,1000),(0.001,1000), (1000,0.001),(1000,0.001),(1000,0.001),(3.,3.),(3.,3.),(3.,3.),(3.,3.)]
		color=['k', 'k','r', 'b', 'g', 'r','b','g', 'r', 'y']
		marker=[u'd', u'*', u'x', u'+', u'.', u'+', u'x', u'+']

		if ("0.0025" in name):
			label=[r"$d_b=2.5mm$ $C_d$ ", r"$d_b=2.5mm$ $C_d^*$", r"Moore (1963)",r"Mei (1994)",r"Schiller (1933)",r"Ishii and Zuber (1979)", r"Tomiyama (1998)"]
		if ("0.0016" in name):
			label=[r"$d_b=1.6mm$ $C_d$ ", r"$d_b=1.6mm$ $C_d^*$", r"Moore (1963)",r"Mei (1994)",r"Schiller (1933)",r"Ishii and Zuber (1979)", r"Tomiyama (1998)"]
		if ("0.0021" in name):
			label=[r"$d_b=2.1mm$ $C_d$ ", r"$d_b=2.1mm$ $C_d^*$", r"Moore (1963)",r"Mei (1994)",r"Schiller (1933)",r"Ishii and Zuber (1979)", r"Tomiyama (1998)"]
		markerfacecols=['none','none','none','none','none','none','none','none','none']
		markersizes=[scalepng*2, scalepng*3, scalepng*2., scalepng*2., scalepng*2., scalepng*2., scalepng*2, scalepng*2, scalepng*2]
		markeredgewidths=[0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5]
		ax.set_xlim([0.,0.105], emit=False)
		ax.set_ylim([0.0,1.2], emit=False)
		ax.set_xlabel(r"${\alpha_v}$", fontsize=15)
		ax.set_ylabel(r"$C_d$", fontsize=15)

	### compa c_lambda / ratio aspect
	if ('lambda_ratio' in name):
		style=['-', '-', '-', '-', '-', '-', '-']
		dash=[(0.0001,1000.), (0.0001,1000.),(1000.,0.0001),(0.0001,1000.),(1000.,0.0001),(1000.,0.0001),(1.,1.)]
		color=['b', 'k', 'r','r', 'r', 'k', 'r', 'g']
		marker=[u'd', u'd', None, u'd', None, u'*', u'*', None, None]
		label=[r"$C_\Lambda$ fit simu", r"$C_\Lambda$ fit exp",r"$3.85$"]
		markerfacecols=['none','none','none','none','none','none','none']
		markersizes=[scalepng*2, scalepng*2, scalepng*2., scalepng*2., scalepng*2., scalepng*2., scalepng*3]
		markeredgewidths=[0.5,0.5,0.5,0.5,0.5,0.5,0.5]
		ax.set_xlim([0.,1.], emit=False)
		ax.set_ylim([2.,4.], emit=False)
		ax.set_xlabel(r"${Eo}$", fontsize=15)
		ax.set_ylabel(r"$C_\Lambda$", fontsize=15)

	### compa c_lambda / ratio aspect
	if ('lambda_eotvos' in name):
		style=['-', '-', '-', '-', '-', '-', '-']
		dash=[(0.0001,1000.), (0.0001,1000.),(0.0001,1000.), (0.0001,1000.),(1000.,0.0001),(0.0001,1000.),(1000.,0.0001),(1000.,0.0001),(1.,1.)]
		color=['b', 'k','b', 'k', 'r','r', 'r', 'k', 'r', 'g']
		marker=[u'.', u'.', u'd', u'd', None, u'd', None, u'*', u'*', None, None]
		label=[r"$C_\Lambda$ DNS $C_d^{real}$", r"$C_\Lambda$ exp $C_d^{real}$",r"$C_\Lambda$ DNS $C_D^{sphe}$", r"$C_\Lambda$ exp $C_D^{sphe}$",r"$2.7$"]
		markerfacecols=['none','none','none','none','none','none','none']
		markersizes=[scalepng*4, scalepng*4, scalepng*2., scalepng*2., scalepng*2., scalepng*2., scalepng*3]
		markeredgewidths=[0.5,0.5,0.5,0.5,0.5,0.5,0.5]
		ax.set_xlim([0.,1.], emit=False)
		ax.set_ylim([2.,4.], emit=False)
		ax.set_xlabel(r"${Eo}$", fontsize=15)
		ax.set_ylabel(r"$C_\Lambda$", fontsize=15)

	if ('_cv' in name):
		style=['-', '-', '-', '-', '-', '-', '-']
		dash=[(100000.,1000000.), (2.,2.),(5,5.), (0.0001,1000.),(1000.,0.0001),(0.0001,1000.),(1000.,0.0001),(1000.,0.0001),(1.,1.)]
		color=['k', 'r','b', 'k', 'r','r', 'r', 'k', 'r', 'g']
		marker=[u'd', u'*', u'.', u'd', None, u'd', None, u'*', u'*', None, None]
		label=[r"$R_{11}$", r"$R_{22}$", r"$R_{33}$",r"$C_\Lambda$ DNS $C_D^{sphe}$", r"$C_\Lambda$ exp $C_D^{sphe}$",r"$2.7$"]
                if ('u' in name):
			label=[r"$u_r$", r"$R_{22}$", r"$R_{33}$",r"$C_\Lambda$ DNS $C_D^{sphe}$", r"$C_\Lambda$ exp $C_D^{sphe}$",r"$2.7$"]
		markerfacecols=['none','none','none','none','none','none','none']
		markersizes=[scalepng*2.5, scalepng*3.5, scalepng*3.5, scalepng*2., scalepng*2., scalepng*2., scalepng*3]
		markeredgewidths=[0.5,0.5,0.5,0.5,0.5,0.5,0.5]
		ax.set_xlim([10.,55.], emit=False)
		ax.set_ylim([0.0005,0.0014], emit=False)
                if ('u' in name):
			ax.set_ylim([0.08,0.11], emit=False)
		ax.set_xlabel(r"$cells/diameter$", fontsize=15)
		ax.set_ylabel(r"$m^2.s^{-2}$", fontsize=15)
                if ('u' in name):
			ax.set_ylabel(r"$m.s^{-1}$", fontsize=15)


	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.xaxis.set_ticks_position('bottom')
	ax.yaxis.set_ticks_position('left') 
	plt.tick_params(axis='both', which='major', labelsize=7)
	plt.rc('font', size=7)
 	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

 	#plt.ticklabel_format(style='plain', axis='y', scilimits=(0,0))
	plt.tight_layout()
	dpii=400
	
	for i in range(len(field)):
		#print i
		### avec legende
        	ax.plot(field[i]._ax[0,:], field[i]._npa[0,:], linestyle=style[i], dashes=dash[i], linewidth=0.4, marker=marker[i], label=label[i], color=color[i], ms=markersizes[i], markevery=0.06, markerfacecolor=markerfacecols[i], markeredgecolor = color[i], markeredgewidth=markeredgewidths[i])
		### avec legende et semilog y
        	#ax.semilogy(field[i]._ax[0,:], field[i]._npa[0,:], linestyle=style[i], dashes=dash[i], linewidth=0.4, marker=marker[i], label=label[i], color=color[i], ms=markersizes[i], markevery=0.06, markerfacecolor=markerfacecols[i], markeredgecolor = color[i], markeredgewidth=markeredgewidths[i])
        	### sans legende
        	#ax.plot(field[i]._ax[0,:], field[i]._npa[0,:], linestyle=style[i], dashes=dash[i], linewidth=0.4, marker=marker[i], color=color[i], ms=markersizes[i], markevery=0.06, markerfacecolor=markerfacecols[i], markeredgecolor = color[i], markeredgewidth=1.0)
	
	f.set_size_inches(scalepng,scalepng)
	box = ax.get_position()
	ax.set_position([box.x0, box.y0, box.width, box.height])
	ax.legend(loc='center left', bbox_to_anchor=(1.1, 0.5), frameon=True)
	inSizeLegend=int(scalepng*3.)
	plt.legend(loc=0,prop={'size':inSizeLegend}, frameon=True, ncol=1)
        #plt.legend().set_visible(False) ### deactive la legende
	print "save ", name+".png"
	plt.savefig(name+'.png', bbox_inches='tight', dpi=dpii)
	plt.close(f)
	return













