import LataLoader
import numpy as np
import rlcompleter, readline
import time 
import os
import MEDLoader as ml
from MEDLoader import MEDLoader
import matplotlib.pyplot as plt
readline.parse_and_bind('tab:complete')
# a=LataLoader.LataLoader("ijkft_bulle_oscillante.lata")
a=LataLoader.LataLoader("ijkft_bulle_oscillante_0000.med")
print a.GetFieldNames()
global compteur
compteur=10000
os.system('rm -f *.png')

#### fonction qui prend en argument la liste des voisins interfaciaux d un element et donne les listes (interfaciaux et phasiques) de tous les voisins en question
def voisinVoisin(surf, desc, descIndx, revDesc, revDescIndx, ai_vi, eul_vi, lien, lien2, phase):
	lien=[]
	lien3=[]
	for k in lien2:
		for i in desc[descIndx[k]:descIndx[k+1]].getValues():                 
			for n in revDesc[revDescIndx[i]:revDescIndx[i+1]].getValues():
				if ai_vi[n]==phase:	
					lien.append(n) 
				else:
					lien3.append(n) 

	return lien, lien3




def listeVoisins(ai_vi,eul_vi, phase):
	surf, desc, descIndx, revDesc, revDescIndx = eul_vi.buildDescendingConnectivity()
	n_maille_3D=eul_vi.getNumberOfCells()  #### recuperer le nombre de mailles 3D
	voisin_lata=LataLoader.DataArrayInt()
	voisin_lata.alloc(eul_vi.getNumberOfCells())
	voisin_lata.fillWithValue(0)
	voisin_complet=[[0]]*eul_vi.getNumberOfCells()
	tol=0.001
	mini=0.0+tol
	maxi=1.0-tol
	print 'recherche des voisins ....'
	print '..........................'
	
	for k in ai_vi.getIdsInRange(mini, maxi).getValues():
		lien=[]
		lien2=[]
		j=0		
		for i in desc[descIndx[k]:descIndx[k+1]].getValues():                   	## pour les element 2D constituant la maille 0 3D
			if lien!=[]:
				break
			for n in revDesc[revDescIndx[i]:revDescIndx[i+1]].getValues():
				if ai_vi[n]==phase:						### si le voisin est dans la bonne phase
					lien.append(n) 						### on l'ajoute a la liste des voisins phasique
				else:								### sinon
					lien2.append(n) 					### on l'ajoute a la liste des voisins interfaciaux
					
												### a ce niveau la, il se peut que certaine liste de voisin soit vide (tous les voisins de la maille sont sur l interface)
												### il faut alors chercher les voisins des voisin
			
				if lien==[]:
					j=j+1						### tant que lien[0] n est pas remplit, on continue l exploration des voisins des voisins jusqu a trouver une valeur phasique
					lien, lien2=voisinVoisin(surf, desc, descIndx, revDesc, revDescIndx, ai_vi, eul_vi, lien, lien2, phase=phase)
				if j>=3:                         ### au dela de trois mailles, le voisin le plus proche ne nous interesse plus
					lien.append(k)					
					break	
		lien=list(set(lien))				### supprimer les doublons	
		voisin_lata[k] = lien[0]  #sum(lien)/len(lien)			### liste du numero de maille du premier voisin phasique dans un format LataLoader.DataArrayInt
		voisin_complet[k] = lien
		
	print 'fin de recherche'
	return voisin_lata, voisin_complet

	


					
def champsInterfacePhase(chi, p, phase):

	voisin_v, voisin_complet=listeVoisins(chi.getArray(),chi.getMesh(), phase)
	n_maille_3D=chi.getMesh().getNumberOfCells() 
	pArray=p.getArray()
	p_v=p*1.0
	tol=0.001
	mini=0.0+tol
	maxi=1.0-tol
	print 'ecriture du champ modifie...'
	#for k in chi.getArray().getIdsInRange(mini, maxi).getValues():	  #### prendre le premier voisin phasiques		     
	#	p_v.getArray()[k]=pArray[voisin_v[k]]
			
	for k in chi.getArray().getIdsInRange(mini, maxi).getValues():    #### faire la moyenne de tous les voisins phasiques
		p_v.getArray()[k]=0.0
		for i in range(len(voisin_complet[k])):
			p_v.getArray()[k]=p_v.getArray()[k]+pArray[voisin_complet[k][i]]
		p_v.getArray()[k]=p_v.getArray()[k]/len(voisin_complet[k])
	print 'fin.........................'
	return p_v






##########################################################################################################

#################################### CALCUL DES M_K ETC... MOYENNE #######################################


def moyStatPlanSurf(p_v, ai):

	coord_z=p_v.getMesh().getCoords()[:,2]
	zLev=coord_z.getDifferentValues(1e-12)
	zLev.sort()
	baryZ=p_v.getMesh().getBarycenterAndOwner()[:,2]
	moy=np.zeros((len(zLev.getValues())-1, 2))
	for i in range(len(zLev.getValues())-1):
		moyenne=0.0
		Surf=0.0
		num=0.0
		for k in baryZ.getIdsInRange(zLev.getValues()[i],zLev.getValues()[i+1]).getValues():
			moyenne=moyenne+p_v.getArray()[k]*ai.getArray()[k]
			#moyenne=moyenne+p_v.getArray()[k]
			Surf=Surf+ai.getArray()[k]
			num=num+1.0
		moy[i,0]=zLev.getValues()[i]
		if Surf==0.0:
			moy[i,1]=0.0
		else:
			moy[i,1]=moyenne/Surf
			#moy[i,1]=moyenne/num
	return moy
	
def moyStatPlan(p_v):

	coord_z=p_v.getMesh().getCoords()[:,2]
	zLev=coord_z.getDifferentValues(1e-12)
	zLev.sort()
	baryZ=p_v.getMesh().getBarycenterAndOwner()[:,2]
	moy=np.zeros((len(zLev.getValues())-1, 2))
	for i in range(len(zLev.getValues())-1):
		moyenne=0.0
		Surf=0.0
		for k in baryZ.getIdsInRange(zLev.getValues()[i],zLev.getValues()[i+1]).getValues():
			moyenne=moyenne+p_v.getArray()[k]
			Surf=Surf+1.0
		moy[i,0]=zLev.getValues()[i]
		moy[i,1]=moyenne/Surf
	
	return moy
	
def tracer(moy, nom):
	global compteur
   	compteur=compteur+1
   	plt.figure(compteur)
	plt.plot(moy[:,0],moy[:,1], label='bla')                          
	plt.legend(loc = 0)
	plt.grid('on')
	plt.savefig(nom, bbox_inches='tight')
	return		



def StatsPressionsInterface(p_l, phase):

	ai=a.GetFieldDouble('AIRE_INTERF_ELEM_DOM_EXT',a.GetNTimesteps()-1 ) # dernier temps
	nx=a.GetFieldDouble('NORMALE_EULER_X_ELEM_DOM_EXT',a.GetNTimesteps()-1 ) # dernier temps
	ny=a.GetFieldDouble('NORMALE_EULER_Y_ELEM_DOM_EXT',a.GetNTimesteps()-1 ) # dernier temps
	nz=a.GetFieldDouble('NORMALE_EULER_Z_ELEM_DOM_EXT',a.GetNTimesteps()-1 ) # dernier temps
	################ calcul (p_l*n)i ##################
	
	if phase=='liq':
		nx=nx*(-1.0)
		ny=ny*(-1.0)
		nz=nz*(-1.0)
	pnx=nx*1.0
	pny=ny*1.0
	pnz=nz*1.0
	pnz.setArray(nz.getArray()*p_l.getArray())
	pnx.setArray(nx.getArray()*p_l.getArray())
	pny.setArray(ny.getArray()*p_l.getArray())
	moy_pnz=moyStatPlanSurf(pnz, ai)
	moy_pny=moyStatPlanSurf(pny, ai)
	moy_pnx=moyStatPlanSurf(pnx, ai)


	############### calcul (p_l)i*(n)i ################



	moy_nx=moyStatPlanSurf(nx, ai)
	moy_ny=moyStatPlanSurf(ny, ai)
	moy_nz=moyStatPlanSurf(nz, ai)
	moy_pl=moyStatPlanSurf(p_l, ai)



	tracer(moyStatPlanSurf(nz, ai), 'nz_ai_'+phase)
	tracer(moyStatPlanSurf(ny, ai), 'ny_ai_'+phase)
	tracer(moyStatPlanSurf(nx, ai), 'nx_ai_'+phase)
	tracer(moyStatPlanSurf(p_l, ai), 'p_ai_'+phase)
	tracer(moyStatPlan(ai), 'ai_'+phase)

	pinxi=np.zeros((len(moy_nx[:,1]),2))
	pinyi=np.zeros((len(moy_nx[:,1]),2))
	pinzi=np.zeros((len(moy_nx[:,1]),2))


	for i in range(len(moy_nx[:,1])):
		pinxi[i,1]=moy_nx[i,1]*moy_pl[i,1]
		pinxi[i,0]=moy_nx[i,0]
		pinyi[i,1]=moy_ny[i,1]*moy_pl[i,1]
		pinyi[i,0]=moy_nx[i,0]
		pinzi[i,1]=moy_nz[i,1]*moy_pl[i,1]
		pinzi[i,0]=moy_nx[i,0]
	


	plt.figure(11)
	plt.plot(pinxi[:,0],pinxi[:,1], label='pinxi') 
	plt.plot(pinyi[:,0],pinyi[:,1], label='pinyi') 
	plt.plot(pinzi[:,0],pinzi[:,1], label='pinzi') 
	plt.plot(moy_pnx[:,0],moy_pnx[:,1], label='pnx')
	plt.plot(moy_pny[:,0],moy_pny[:,1], label='pny')
	plt.plot(moy_pnz[:,0],moy_pnz[:,1], label='pnz')
	plt.legend(loc = 0)
	plt.grid('on')
	plt.savefig('pn_'+phase, bbox_inches='tight')
			

	############## les differents champs de pression ########
	############## saut de pression et tension de surface ###

	moy_v=moyStatPlanSurf(p_v, ai)
	moy_l=moyStatPlanSurf(p_l, ai)
	moy_p=moyStatPlanSurf(p, ai)
	
				
	plt.figure(12)
	plt.plot(moy_v[:,0],moy_v[:,1], label='vapeur') 
	plt.plot(moy_l[:,0],moy_l[:,1], label='liquide') 
	plt.plot(moy_p[:,0],moy_p[:,1], label='interp')                           
	plt.legend(loc = 0)
	plt.grid('on')
	plt.savefig('p_surf', bbox_inches='tight')
	
	
	
	moy_v=moyStatPlan(p_v)
	moy_l=moyStatPlan(p_l)
	moy_p=moyStatPlan(p)
	
				
	plt.figure(15)
	plt.plot(moy_v[:,0],moy_v[:,1], label='vapeur') 
	plt.plot(moy_l[:,0],moy_l[:,1], label='liquide') 
	plt.plot(moy_p[:,0],moy_p[:,1], label='interp')                           
	plt.legend(loc = 0)
	plt.grid('on')
	plt.savefig('p', bbox_inches='tight')
		
	
	######### modele force de drag #######################

	fx=np.zeros((len(moy_nx[:,1]),2))
	fy=np.zeros((len(moy_nx[:,1]),2))
	fz=np.zeros((len(moy_nx[:,1]),2))


	for i in range(len(moy_nx[:,1])):
		fx[i,1]=pinxi[i,1]-moy_pnx[i,1]
		fx[i,0]=pinxi[i,0]
		fy[i,1]=pinyi[i,1]-moy_pny[i,1]
		fy[i,0]=pinxi[i,0]
		fz[i,1]=pinzi[i,1]-moy_pnz[i,1]
		fz[i,0]=pinxi[i,0]

	

			
	plt.figure(13)
	plt.plot(fx[:,0],fx[:,1], label='fx_drag') 
	plt.plot(fy[:,0],fy[:,1], label='fy_drag') 
	plt.plot(fz[:,0],fz[:,1], label='fz_drag')    
	plt.plot(pinxi[:,0],pinxi[:,1], label='fx_sigma') 
	plt.plot(pinyi[:,0],pinyi[:,1], label='fy_sigma') 
	plt.plot(pinzi[:,0],pinzi[:,1], label='fz_sigma') 	
	plt.legend(loc = 0)
	plt.grid('on')
	plt.savefig('force_'+phase, bbox_inches='tight')
	
	return	
	

##########################################################################################################
##################################### MAIN ###############################################################
##########################################################################################################


if 0 :
	p=a.GetFieldDouble('PRESSURE_ELEM_DOM',a.GetNTimesteps()-1 ) # dernier temps
	chi=a.GetFieldDouble('INDICATRICE_ELEM_DOM',a.GetNTimesteps()-1 ) # dernier temps
	p_v=champsInterfacePhase(chi, p, 0.0)
	MEDLoader.WriteField("p_v.med",p_v,True)
	p_l=champsInterfacePhase(chi, p, 1.0)
	MEDLoader.WriteField("p_l.med",p_l,True)
	
	
p=a.GetFieldDouble('PRESSURE_ELEM_DOM',a.GetNTimesteps()-1 ) # dernier temps	
p_v = MEDLoader.ReadFieldCell("p_v.med", 'DOM', 0, 'PRESSURE_ELEM_DOM', 1, -1)
p_l = MEDLoader.ReadFieldCell("p_l.med", 'DOM', 0, 'PRESSURE_ELEM_DOM', 1, -1)



StatsPressionsInterface(p_l, 'liq')
StatsPressionsInterface(p_v, 'gaz')	
	

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
			
					
					
					
					
					
