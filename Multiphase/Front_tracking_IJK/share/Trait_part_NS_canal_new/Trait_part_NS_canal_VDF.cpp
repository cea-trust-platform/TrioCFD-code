//////////////////////////////////////////////////////////////////////////////
//
// File:        Trait_part_NS_canal_VDF.cpp
// Directory:   $TRIO_U_ROOT/VDF
// Version:     /main/22
//
//////////////////////////////////////////////////////////////////////////////

#include <Trait_part_NS_canal_VDF.h>
#include <Zone_VDF.h>
//FA 3/02/11
#include <Zone_VF.h>

#include <Pave.h>

#include <Fluide_Inc.h>
#include <N_S_Turb.h>
#include <VectEsp_Dist.h>
// modified AT 5/06/09
#include <C_D_Turb.h>
#include <C_D_Turb_Chaleur_QC.h>
#include <C_D_Turb_C.h>
#include <C_D_Turb_T.h>
#include <Eqn_base.h>
#include <Pb_base.h>
// fin modif AT 5/06/09
#include <Champ_Don.h>	//modif YB 30/6/09


Implemente_instanciable(Traitement_particulier_NS_canal_VDF,"Traitement_particulier_NS_canal_VDF",Traitement_particulier_NS_canal);


  // Description: 
  //    
  // Precondition: 
  // Parametre: Sortie& is
  //    Signification: un flot de sortie
  //    Valeurs par defaut: 
  //    Contraintes: 
  //    Acces: entree/sortie
  // Retour: Sortie&
  //    Signification: le flot de sortie modifie
  //    Contraintes: 
  // Exception: 
  // Effets de bord: 
  // Postcondition: la methode ne modifie pas l'objet 
Sortie& Traitement_particulier_NS_canal_VDF::printOn(Sortie& is) const
{
  return is;
}


// Description: 
//   
// Precondition: 
// Parametre: Entree& is
//    Signification: un flot d'entree
//    Valeurs par defaut: 
//    Contraintes: 
//    Acces: entree/sortie
// Retour: Entree& 
//    Signification: le flot d'entree modifie
//    Contraintes: 
// Exception: 
// Effets de bord: 
// Postcondition: 
Entree& Traitement_particulier_NS_canal_VDF::readOn(Entree& is)
{
  return is;
}

Entree& Traitement_particulier_NS_canal_VDF::lire(Entree& is)
{
  return Traitement_particulier_NS_canal::lire(is);
}

void Traitement_particulier_NS_canal_VDF::remplir_Y(DoubleVect& Y,  DoubleVect& compt, entier& Ny) const
{
  // On va initialiser les differents parametres membres de la classe 
  // utiles au calcul des differentes moyennes
  // Initialisation de : Y, compt
  
  const Zone_dis_base& zdisbase=mon_equation->inconnue().zone_dis_base();     
  const Zone_VF& zone_VF=ref_cast(Zone_VF, zdisbase); 
  const DoubleTab& xp = zone_VF.xp();
  entier nb_elems = zone_VF.zone().nb_elem();
  entier num_elem,j,indic,trouve;
  double y;
  
  j=0;
  indic = 0;

  Y.resize(1);
  compt.resize(1);

  Y = -100.;
  compt = 0;
  
  //Remplissage du tableau Y
  ////////////////////////////////////////////////////////

  for (num_elem=0;num_elem<nb_elems;num_elem++)
    {
      y = xp(num_elem,1);
      trouve = 0;

      for (j=0;j<indic+1;j++)
	{
	  if(est_egal(y,Y[j]))
	    {
	      compt[j] ++;
	      j=indic+1;
	      trouve = 1;
	    }
	}
      if (trouve==0)
	{	
	  Y[indic]=y;
	  compt[indic] ++;
	  indic++;
	 
	 Y.resize(indic+1);
	 Y(indic)=-100.;
         compt.resize(indic+1);
	}
    }

  Ny = indic;

  Y.resize(Ny);
  compt.resize(Ny);
}

//Ajout F.A 15/02/11 on va faire un gros changement, 
// l'objectif est de réunir des opérations faites et refaite pour profiter pleinement de l'espace mémoire (boucles)
// par la création d'un tableau de grande taille (+/- 7M par proc mais qui ne s'échange pas)
// le tableau aura la structure suivant : La ligne est le numéro de l'élémennt 
// numéro de l'élément au dessus, numéro de l'élément au dessous,position dans le vecteur Y.
// soit un tableau de nelem x 3.
// pour cela après remplir_Y on va appeller la fonction qui fais les différénent calculs,
// en utilisant le tableau comme argument de la fonction.

void Traitement_particulier_NS_canal_VDF::remplir_Tab_recap(DoubleTab& Tab_recap)
{
  const Zone_dis_base& zdisbase=mon_equation->inconnue().zone_dis_base();     
  const Zone_VDF& zone_VDF=ref_cast(Zone_VDF, zdisbase);
  const DoubleTab& xp = zone_VDF.xp();
  const IntTab& elem_faces = zone_VDF.elem_faces();
 
  int face; //récepteur des faces
  int elem_test,elem_test2; // élément de test pour les ficitfs
  int nb_elem_tot = zone_VDF.zone().nb_elem_tot(); // nombre total d'éléments (réel + fict)
  int nb_elems = zone_VDF.zone().nb_elem();
  int dimension=Objet_U::dimension;
  
  IntTab trouve(1);// tableau des éléments déja effectué
  double y=0;
  int i,num_elem; // compteurs
  int q=1; //Curseur pour les tableau haut
  trouve[0]=0;

  
 Tab_recap.resize(nb_elems,3); // On dimenssione le tableau. 
  
    for (num_elem=nb_elems;num_elem<nb_elem_tot;num_elem++) // boucle sur les éléments fictifs
      {
    	face = elem_faces(num_elem,1+dimension);	
	elem_test=zone_VDF.elem_voisin(num_elem,face,0);
	face = elem_faces(num_elem,1);
	elem_test2=zone_VDF.elem_voisin(num_elem,face,1);
	
	if ((elem_test>0) && (elem_test<nb_elems)) // si l'élément en dessus est un élément réel alors
	{
	trouve[q-1]=elem_test;
	q =q +1;
	trouve.resize(q);
	
	Tab_recap(elem_test,0)=num_elem; // on affecte la même valeur aux deux case haut et bas 
	Tab_recap(elem_test,1)=num_elem; //ainsi la fonction qui calcul les valeurs voie un élément normal.
	
	y=xp(elem_test,1);
	for (i=0;i<Ny;i++)
	  if(est_egal(y,Y[i])) break;
	  
	Tab_recap(elem_test,2)=i; // on garde la valeur de i pour ne pas rééfectuer la boucle a chaque pas de temps.
	}
	else if ((elem_test2<nb_elems)&&(elem_test2>0)) //sinon si l'élément en dessous est un élément réel alors
  	{
  	trouve[q-1]=elem_test2;
	q =q +1;
	trouve.resize(q);
	
	Tab_recap(elem_test2,0)=num_elem; // on affecte la même valeur aux deux case haut et bas 
	Tab_recap(elem_test2,1)=num_elem; //ainsi la fonction qui calcul les valeurs voie un élément normal.
	
	y=xp(elem_test2,1);
	for (i=0;i<Ny;i++)
	  if(est_egal(y,Y[i])) break;
	Tab_recap(elem_test2,2)=i; // on garde la valeur de i pour ne pas rééfectuer la boucle a chaque pas de temps.
		
	} 
  	// sinon rien
      }
 
 	Cerr << "Traitement particulier canal : Il y a une amélioration a apporter aux face de bord !! " << finl;
	for (num_elem=0;num_elem<nb_elems;num_elem++)
	{ 
	 q=0;// on utilise le compteur q qui ne nous sert plus pour vérifier si on a trouver un équivalent.
	 for(i=0;i<(trouve.size()-1);i++) // trouve est une case trop grand, mais plutot que de le redimentionner on utilise le critère taille -1
	  if((num_elem==trouve[i]))
	  { q = 0; break; } // on met fixe q qui ne peu répondre au prochain test. // correction on fixe q =0 car c'étais un faux problème.
	  // en réalité nu_t/lambda_smt explose a l'interface.
	 if(q==0) //       
	 {
	face=elem_faces(num_elem,1); //face inférieure
	elem_test=zone_VDF.elem_voisin(num_elem,face,1);
	
	
	if (elem_test+1) {Tab_recap(num_elem,1)=elem_test;} // faux si elem_test=-1 sinon remplit avec l'élément en dessous
	else {Tab_recap(num_elem,1)=zone_VDF.elem_voisin(num_elem,elem_faces(num_elem,1+dimension),0);} // on le traite alors comme un virtuel 
	
	face= elem_faces(num_elem,1+dimension); //face supérieure
	elem_test=zone_VDF.elem_voisin(num_elem,face,0);
	
	if (elem_test+1) {Tab_recap(num_elem,0)=elem_test;} // faux si elem_test=-1 sinon remplit avec l'élément au dessus
	else {Tab_recap(num_elem,0)=zone_VDF.elem_voisin(num_elem,elem_faces(num_elem,1),1);} // on le traite alors comme un virtuel 
	
	 y = xp(num_elem,1);	 
	 for (i=0;i<Ny;i++)
	  if(est_egal(y,Y[i])) break;
	 
	 Tab_recap(num_elem,2)=i;
	 
	 } 
	 
	}
}

void Traitement_particulier_NS_canal_VDF::calculer_moyenne_spatiale_vitesse_rho_mu(DoubleTab& val_moy) const
{
  const Zone_dis_base& zdisbase=mon_equation->inconnue().zone_dis_base();     
  const Zone_VDF& zone_VDF=ref_cast(Zone_VDF, zdisbase);
  const DoubleTab& xp = zone_VDF.xp();
  const IntTab& elem_faces = zone_VDF.elem_faces();
  const DoubleTab& vitesse = mon_equation->inconnue().valeurs();
  double y,u,v,w;
  entier nb_elems = zone_VDF.zone().nb_elem();
  entier num_elem,i;
  entier face_x_0,face_x_1,face_y_0,face_y_1,face_z_0,face_z_1;
  
  int dimension=Objet_U::dimension;

  const Fluide_Incompressible& le_fluide = ref_cast(Fluide_Incompressible,mon_equation->milieu());
  const DoubleTab& visco_dyn = le_fluide.viscosite_dynamique();
  const DoubleTab& tab_rho_elem = le_fluide.masse_volumique();
  int taille_mu=visco_dyn.dimension(0);
  int taille_rho=tab_rho_elem.dimension(0);

  for (num_elem=0;num_elem<nb_elems;num_elem++)
  {
    y=xp(num_elem,1);
    
    face_x_0 = elem_faces(num_elem,0);
    face_x_1 = elem_faces(num_elem,dimension);
    face_y_0 = elem_faces(num_elem,1);
    face_y_1 = elem_faces(num_elem,1+dimension);
    
    // PQ : 12/10 : pour eviter de moyenner localement la vitesse u et w
    //  	    on fait le choix de "deplacer" celles-ci 
    //  	    du centre de la face au centre de l'element
    //  	    en se basant sur le principe que l'ecoulement est homogene
    //		    suivant les plans xz
    //		    Pour v, la moyenne s'impose si l'on veut revenir au centre des elements

//  u = .5*(vitesse[face_x_0]+vitesse[face_x_1]);
    u = vitesse[face_x_0];
    v = .5*(vitesse[face_y_0]+vitesse[face_y_1]);    

    i=Tab_recap(num_elem,2);
      
      val_moy(i,0) += u;	
      val_moy(i,1) += v;    
      val_moy(i,3) += u*u;  
      val_moy(i,4) += v*v;  
      val_moy(i,6) += u*v;  

      if(dimension==2)   val_moy(i,9) += sqrt(u*u);      //vitesse tangentielle pour calcul du frottement
    
      if(dimension==3)
      {
       face_z_0 = elem_faces(num_elem,2);
       face_z_1 = elem_faces(num_elem,2+dimension);
 
//     w = .5*(vitesse[face_z_0]+vitesse[face_z_1]);
       w = vitesse[face_z_0];
     
       val_moy(i,2) += w;   
       val_moy(i,5) += w*w; 
       val_moy(i,7) += u*w; 
       val_moy(i,8) += v*w;
       val_moy(i,9) += sqrt(u*u+w*w);      //vitesse tangentielle pour calcul du frottement
      }

    if (taille_rho==1)  val_moy(i,10) += tab_rho_elem(0,0);
    else		val_moy(i,10) += tab_rho_elem[num_elem];
    

    if (taille_mu==1)   val_moy(i,11) += visco_dyn(0,0);
    else		val_moy(i,11) += visco_dyn[num_elem];
  }
}

void Traitement_particulier_NS_canal_VDF::calculer_moyenne_spatiale_nut(DoubleTab& val_moy) const
{
  const Zone_dis_base& zdisbase=mon_equation->inconnue().zone_dis_base();     
  const Zone_VDF& zone_VDF=ref_cast(Zone_VDF, zdisbase);
  const DoubleTab& xp = zone_VDF.xp();
  const Navier_Stokes_Turbulent& N_S_Turb  = ref_cast(Navier_Stokes_Turbulent,mon_equation.valeur());
  const DoubleTab& nu_t = N_S_Turb.viscosite_turbulente().valeurs();
  entier nb_elems = zone_VDF.zone().nb_elem();
  entier num_elem,i;
  double y;

  for (num_elem=0;num_elem<nb_elems;num_elem++)
  {
   
    y=xp(num_elem,1);
    
   i=Tab_recap(num_elem,2);

    val_moy(i,12) += nu_t[num_elem];	
  }
}

void Traitement_particulier_NS_canal_VDF::calculer_moyenne_spatiale_Temp(DoubleTab& val_moy) const
{
  //modified AT 5/06/09
  const Probleme_base& pb = mon_equation->probleme();
  const Equation_base& eqn_th = pb.equation(1);


 // on identifie le type de problème puis on appel la fonction qui calcule.
    if (same_type(Convection_Diffusion_Temperature_Turbulent,eqn_th)) {
    const Convection_Diffusion_Temperature_Turbulent& eqn_thermo = ref_cast(Convection_Diffusion_Temperature_Turbulent,eqn_th);
    const DoubleTab& diffusivite_turb =eqn_thermo.diffusivite_turbulente().valeurs();
    calculer_Temp(val_moy, diffusivite_turb);
    }
    else if (same_type(Convection_Diffusion_Chaleur_Turbulent_QC,eqn_th)) {
    const Convection_Diffusion_Chaleur_Turbulent_QC& eqn_thermo = ref_cast(Convection_Diffusion_Chaleur_Turbulent_QC,eqn_th);		
    const DoubleTab& diffusivite_turb =eqn_thermo.diffusivite_turbulente().valeurs();
    calculer_Temp(val_moy, diffusivite_turb);
    }
    else if (same_type(Convection_Diffusion_Concentration_Turbulent,eqn_th)) {
    const Convection_Diffusion_Concentration_Turbulent& eqn_thermo = ref_cast(Convection_Diffusion_Concentration_Turbulent,eqn_th);		
    const DoubleTab& diffusivite_turb =eqn_thermo.diffusivite_turbulente().valeurs();
    calculer_Temp(val_moy, diffusivite_turb);
    }
    else {
    Cerr << " on doit forcement avoir un Pb_Thermohydraulique_Turbulent ou un Pb_Thermohydraulique_Turbulent_QC ou un Pb_Thermohydraulique_Concentration_Turbulent" << finl;
    exit(-1);
  	}	
  
}
void Traitement_particulier_NS_canal_VDF::calculer_Temp(DoubleTab& val_moy,const DoubleTab& diffusivite_turb) const
{
const Zone_dis_base& zdisbase=mon_equation->inconnue().zone_dis_base();     
  const Zone_VDF& zone_VDF=ref_cast(Zone_VDF, zdisbase);
  const DoubleTab& xp = zone_VDF.xp();
  const IntTab& elem_faces = zone_VDF.elem_faces();
  const DoubleTab& temperature = Temp.valeur().valeurs();
  const DoubleTab& vitesse = mon_equation->inconnue().valeurs();

  const DoubleTab& vitesse_np1 = mon_equation->inconnue().futur();
  const DoubleTab& vitesse_nm1 = mon_equation->inconnue().passe();

  double y,u,v,v_nm1,v_np1,w,T,lambda_sm,rho,Tb,Th,dyb,dyh; //modified YB 29/06/09
  double lambdah,lambdab,lambdae,lambda_smh,lambda_smb,nuth,nutb; //modified AT 03/09/09
  
  entier nb_elems = zone_VDF.zone().nb_elem();
  entier num_elem,i;
  entier face_x_0,face_x_1,face_y_0,face_y_1,face_z_0,face_z_1;
  
   
  const Fluide_Quasi_Compressible& le_fluide = ref_cast(Fluide_Quasi_Compressible,mon_equation->milieu());//modif AT 03/09/09
  const DoubleTab& tab_rho_elem = le_fluide.masse_volumique();
  //const DoubleTab& tab_rho_face = le_fluide.rho_face_n(); //modif AT 03/09/09
  const DoubleTab& lambda = le_fluide.conductivite(); //modif YB 30/06/09
  int taille_rho=tab_rho_elem.dimension(0);  
  const Navier_Stokes_Turbulent& N_S_Turb  = ref_cast(Navier_Stokes_Turbulent,mon_equation.valeur());	//modif YB 30/06/09
  const DoubleTab& nu_t = N_S_Turb.viscosite_turbulente().valeurs();	//modif YB 30/06/09
  
  w=0.;
   
    for (num_elem=0;num_elem<nb_elems;num_elem++)
      { 
	y=xp(num_elem,1);
	i=Tab_recap(num_elem,2);
	int elem_haut=Tab_recap(num_elem,0);
	int elem_bas=Tab_recap(num_elem,1);
	
	
	if (taille_rho==1)  rho = tab_rho_elem(0,0);
	else		rho = tab_rho_elem[num_elem];
		
	face_x_0 = elem_faces(num_elem,0);
	face_x_1 = elem_faces(num_elem,dimension);
	face_y_0 = elem_faces(num_elem,1);
	face_y_1 = elem_faces(num_elem,1+dimension);
	
	u = .5*(vitesse[face_x_0]+vitesse[face_x_1]);
	v = .5*(vitesse[face_y_0]+vitesse[face_y_1]);
	
	v_nm1 = .5*(vitesse_nm1[face_y_0]+vitesse_nm1[face_y_1]);
	v_np1 = .5*(vitesse_np1[face_y_0]+vitesse_np1[face_y_1]);      
	
	T = temperature[num_elem];
	Th = temperature[elem_haut];
	Tb = temperature[elem_bas];
	dyh = xp(elem_haut,1)-y;
	dyb = y-xp(elem_bas,1);
		
	lambdah = lambda(elem_haut);
	lambdab = lambda(elem_bas);
	lambdae = lambda(num_elem);
	nuth = nu_t[elem_haut];
	nutb = nu_t[elem_bas];
	lambda_smh=diffusivite_turb[elem_haut];
	lambda_sm = diffusivite_turb[num_elem];
	lambda_smb=diffusivite_turb[elem_bas];
	
//      if (i==32) 
//	{
//       Cerr << " oooo " << y << " " << ((lambda_smb-lambda_sm)*100./lambda_sm) << " " << ((lambda_smh-lambda_sm)*100./lambda_sm) << " "  << ((lambdab-lambdae)*100./lambdae) << " " << ((lambdah-lambdae)*100./lambdae) << finl;
//	if(abs(((lambda_smb-lambda_sm)*100./lambda_sm))>50) { Cerr << " IIII "<< num_elem << " Ltb " << ((lambda_smb-lambda_sm)*100./lambda_sm) ;}
//	if(abs(((lambda_smh-lambda_sm)*100./lambda_sm))>50) { Cerr << " IIII "<<  num_elem << " Lth " << ((lambda_smh-lambda_sm)*100./lambda_sm) ;}
//	if(abs(((lambdab-lambdae)*100./lambdae))>50) { Cerr << " III " << num_elem << " Lb " <<((lambdab-lambdae)*100./lambdae);}
//	if(abs(((lambdah-lambdae)*100./lambdae))>50) { Cerr << " III " << num_elem << " Lh " <<((lambdah-lambdae)*100./lambdae);}
//	}
	
	
	val_moy(i,13) += T;	
	val_moy(i,14) += T*T;    
	val_moy(i,15) += u*T;  
	val_moy(i,16) += v*T;  
	val_moy(i,18) += lambda_sm;
	val_moy(i,19) += rho*T;
	val_moy(i,20) += rho*T*v_nm1;//rho*u*T;
	val_moy(i,21) += rho*T*v;
	val_moy(i,23) += rho*T*v_np1;//rho*u;
	val_moy(i,24) += rho*v;
	val_moy(i,26) += rho*u*u;
	val_moy(i,27) += rho*v*v;
	val_moy(i,29) += rho*u*v;
	
	if(Objet_U::dimension==3)
	  {
	    face_z_0 = elem_faces(num_elem,2);
	    face_z_1 = elem_faces(num_elem,2+dimension);
	    
	    w = .5*(vitesse[face_z_0]+vitesse[face_z_1]);
	    
	    val_moy(i,17) += w*T;   
	    val_moy(i,22) += rho*w*T;   
	    val_moy(i,25) += rho*w;
	    val_moy(i,28) += rho*w*w;
	    val_moy(i,30) += rho*u*w;
	    val_moy(i,31) += rho*v*w;
	  }
	val_moy(i,18) += lambda_sm;	
	
	val_moy(i,32) += dyh/(dyb*(dyb+dyh))*(T-Tb)+dyb/(dyh*(dyb+dyh))*(Th-T);		//modif YB 29/06/09
	//modif AT 10/05/2010	
	val_moy(i,33) += -0.25*((lambdae+lambdah)*(Th-T)/dyh+(lambdae+lambdab)*(T-Tb)/dyb);  
	val_moy(i,34) += -0.25*((lambda_sm+lambda_smh)*(Th-T)/dyh+(lambda_sm+lambda_smb)*(T-Tb)/dyb);
      
    val_moy(0,32)=val_moy(1,32);
    val_moy(0,33)=val_moy(1,33);
    val_moy(0,34)=val_moy(1,34);
    val_moy(Ny-1,32)=val_moy(Ny-2,32);
    val_moy(Ny-1,33)=val_moy(Ny-2,33);
    val_moy(Ny-1,34)=val_moy(Ny-2,34);
      
    }
  
}
//fin modif AT 5/06/09

