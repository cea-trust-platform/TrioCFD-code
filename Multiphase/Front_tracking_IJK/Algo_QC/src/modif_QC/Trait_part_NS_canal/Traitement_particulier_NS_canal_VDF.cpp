/****************************************************************************
* Copyright (c) 2015, CEA
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
* 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
* 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
* IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
* OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
*****************************************************************************/
//////////////////////////////////////////////////////////////////////////////
//
// File:        Traitement_particulier_NS_canal_VDF.cpp
// Directory:   $TRIO_U_ROOT/VDF/Turbulence
// Version:     /main/29
//
//////////////////////////////////////////////////////////////////////////////

#include <Traitement_particulier_NS_canal_VDF.h>
#include <Zone_VDF.h>
//FA 3/02/11
#include <Zone_VF.h>

#include <Pave.h>

#include <Fluide_Incompressible.h>
#include <Navier_Stokes_Turbulent.h>
// modified AT 5/06/09
#include <Convection_Diffusion_Turbulent.h>
#include <Convection_Diffusion_Chaleur_Turbulent_QC.h>
#include <Convection_Diffusion_Concentration_Turbulent.h>
#include <Convection_Diffusion_Temperature_Turbulent.h>
#include <Equation_base.h>
#include <Probleme_base.h>
// fin modif AT 5/06/09
#include <Champ_Don.h>	//modif YB 30/6/09


Implemente_instanciable(Traitement_particulier_NS_canal_VDF,"Traitement_particulier_NS_canal_VDF",Traitement_particulier_NS_canal);


/*! @brief
 *
 * @param (Sortie& is) un flot de sortie
 * @return (Sortie&) le flot de sortie modifie
 */
Sortie& Traitement_particulier_NS_canal_VDF::printOn(Sortie& is) const
{
  return is;
}


/*! @brief
 *
 * @param (Entree& is) un flot d'entree
 * @return (Entree&) le flot d'entree modifie
 */
Entree& Traitement_particulier_NS_canal_VDF::readOn(Entree& is)
{
  return is;
}

Entree& Traitement_particulier_NS_canal_VDF::lire(Entree& is)
{
  return Traitement_particulier_NS_canal::lire(is);
}

void Traitement_particulier_NS_canal_VDF::remplir_Y(DoubleVect& Yc,  DoubleVect& compteur, int& N_y) const
{
  // On va initialiser les differents parametres membres de la classe
  // utiles au calcul des differentes moyennes
  // Initialisation de : Y, compt

  const Zone_dis_base& zdisbase=mon_equation->inconnue().zone_dis_base();
  const Zone_VF& zone_VF=ref_cast(Zone_VF, zdisbase);
  const DoubleTab& xp = zone_VF.xp();
  int nb_elems = zone_VF.zone().nb_elem();
  int num_elem,j,indic,trouve;
  double y;

  j=0;
  indic = 0;

  Yc.resize(1);
  compteur.resize(1);

  Yc = -100.;
  compteur = 0;

  //Remplissage du tableau Y
  ////////////////////////////////////////////////////////

  for (num_elem=0; num_elem<nb_elems; num_elem++)
    {
      y = xp(num_elem,1);
      trouve = 0;

      for (j=0; j<indic+1; j++)
        {
          if(est_egal(y,Y[j]))
            {
              compteur[j] ++;
              j=indic+1;
              trouve = 1;
            }
        }
      if (trouve==0)
        {
          Yc[indic]=y;
          compteur[indic] ++;
          indic++;

          Yc.resize(indic+1);
          Yc(indic)=-100.;
          compteur.resize(indic+1);
        }
    }

  N_y = indic;

  Yc.resize(N_y);
  compteur.resize(N_y);
}

//Ajout F.A 15/02/11 on va faire un changement,
// l'objectif est de reunir des operations faites et refaite pour profiter pleinement de l'espace memoire (boucles)
// par la creation d'un tableau de grande taille (+/- 7M par proc mais qui ne s'echange pas)
// le tableau aura la structure suivant : La ligne est le numero de l'elemennt
// numero de l'element au dessus, numero de l'element au dessous,position dans le vecteur Y.
// soit un tableau de nelem x 3.
// pour cela apres remplir_Y on va appeller la fonction qui fais les differenent calculs,
// en utilisant le tableau comme argument de la fonction.

void Traitement_particulier_NS_canal_VDF::remplir_Tab_recap(IntTab& Tab_rec) const
{
  const Zone_dis_base& zdisbase=mon_equation->inconnue().zone_dis_base();
  const Zone_VDF& zone_VDF=ref_cast(Zone_VDF, zdisbase);
  const DoubleTab& xp = zone_VDF.xp();
  const IntTab& elem_faces = zone_VDF.elem_faces();

  int face; //recepteur des faces
  int elem_test,elem_test2; // element de test pour les ficitfs
  int nb_elem_tot = zone_VDF.zone().nb_elem_tot(); // nombre total d'elements (reel + fict)
  int nb_elems = zone_VDF.zone().nb_elem();
  int dimension=Objet_U::dimension;

  IntTab trouve(1);// tableau des elements deja effectue
  double y=0;
  int i,num_elem; // compteurs
  int q=1; //Curseur pour les tableau haut
  trouve[0]=0;


  Tab_rec.resize(nb_elems,3); // On dimenssione le tableau.

  for (num_elem=nb_elems; num_elem<nb_elem_tot; num_elem++) // boucle sur les elements fictifs
    {
      face = elem_faces(num_elem,1+dimension);
      elem_test=zone_VDF.elem_voisin(num_elem,face,0);
      face = elem_faces(num_elem,1);
      elem_test2=zone_VDF.elem_voisin(num_elem,face,1);

      if ((elem_test>0) && (elem_test<nb_elems)) // si l'element en dessus est un element reel alors
        {
          trouve[q-1]=elem_test;
          q =q +1;
          trouve.resize(q);

          Tab_rec(elem_test,0)=num_elem; // on affecte la meme valeur aux deux case haut et bas
          Tab_rec(elem_test,1)=num_elem; //ainsi la fonction qui calcul les valeurs voie un element normal.

          y=xp(elem_test,1);
          for (i=0; i<Ny; i++)
            {
              if(est_egal(y,Y[i]))
                break;
            }

          Tab_rec(elem_test,2)=i; // on garde la valeur de i pour ne pas reefectuer la boucle a chaque pas de temps.
        }
      else if ((elem_test2<nb_elems)&&(elem_test2>0)) //sinon si l'element en dessous est un element reel alors
        {
          trouve[q-1]=elem_test2;
          q =q +1;
          trouve.resize(q);

          Tab_rec(elem_test2,0)=num_elem; // on affecte la meme valeur aux deux case haut et bas
          Tab_rec(elem_test2,1)=num_elem; //ainsi la fonction qui calcul les valeurs voie un element normal.

          y=xp(elem_test2,1);
          for (i=0; i<Ny; i++)
            {
              if(est_egal(y,Y[i]))
                break;
            }
          Tab_rec(elem_test2,2)=i; // on garde la valeur de i pour ne pas reefectuer la boucle a chaque pas de temps.

        }
      // sinon rien
    }

  Cerr << "Traitement particulier canal : Il y a une amelioration a apporter aux face de bord !! " << finl;
  for (num_elem=0; num_elem<nb_elems; num_elem++)
    {
      q=0;// on utilise le compteur q qui ne nous sert plus pour verifier si on a trouver un equivalent.
      for(i=0; i<(trouve.size()-1); i++) // trouve est une case trop grand, mais plutot que de le redimentionner on utilise le critere taille -1
        if((num_elem==trouve[i]))
          {
            q = 0;  // on met fixe q qui ne peu repondre au prochain test. // correction on fixe q =0 car c'etais un faux probleme.
            break;
          }
      // en realite lambda explose a l'interface.
      if(q==0) //
        {
          face=elem_faces(num_elem,1); //face inferieure
          elem_test=zone_VDF.elem_voisin(num_elem,face,1);


          if (elem_test+1)
            {
              Tab_rec(num_elem,1)=elem_test; // faux si elem_test=-1 sinon remplit avec l'element en dessous
            }
          else
            {
              Tab_rec(num_elem,1)=zone_VDF.elem_voisin(num_elem,elem_faces(num_elem,1+dimension),0); // on le traite alors comme un virtuel
            }

          face= elem_faces(num_elem,1+dimension); //face superieure
          elem_test=zone_VDF.elem_voisin(num_elem,face,0);

          if (elem_test+1)
            {
              Tab_rec(num_elem,0)=elem_test; // faux si elem_test=-1 sinon remplit avec l'element au dessus
            }
          else
            {
              Tab_rec(num_elem,0)=zone_VDF.elem_voisin(num_elem,elem_faces(num_elem,1),1); // on le traite alors comme un virtuel
            }

          y = xp(num_elem,1);
          for (i=0; i<Ny; i++)
            if(est_egal(y,Y[i])) break;

          Tab_rec(num_elem,2)=i;

        }

    }
}

void Traitement_particulier_NS_canal_VDF::calculer_moyenne_spatiale_vitesse_rho_mu(DoubleTab& val_moy) const
{
  const Zone_dis_base& zdisbase=mon_equation->inconnue().zone_dis_base();
  const Zone_VDF& zone_VDF=ref_cast(Zone_VDF, zdisbase);
  //  const DoubleTab& xp = zone_VDF.xp();
  const IntTab& elem_faces = zone_VDF.elem_faces();
  const DoubleTab& vitesse = mon_equation->inconnue().valeurs();
  double u,v,wl;
  int nb_elems = zone_VDF.zone().nb_elem();
  int num_elem,i;
  int face_x_0,face_y_0,face_y_1,face_z_0;

  int dimension=Objet_U::dimension;

  const Fluide_Incompressible& le_fluide = ref_cast(Fluide_Incompressible,mon_equation->milieu());
  const DoubleTab& visco_dyn = le_fluide.viscosite_dynamique();
  const DoubleTab& tab_rho_elem = le_fluide.masse_volumique();
  int taille_mu=visco_dyn.dimension(0);
  int taille_rho=tab_rho_elem.dimension(0);

  for (num_elem=0; num_elem<nb_elems; num_elem++)
    {
      //y=xp(num_elem,1);

      face_x_0 = elem_faces(num_elem,0);
      //face_x_1 = elem_faces(num_elem,dimension);
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

      i= Tab_recap(num_elem,2);

      val_moy(i,0) += u;
      val_moy(i,1) += v;
      val_moy(i,3) += u*u;
      val_moy(i,4) += v*v;
      val_moy(i,6) += u*v;

      if(dimension==2)   val_moy(i,9) += sqrt(u*u);      //vitesse tangentielle pour calcul du frottement

      if(dimension==3)
        {
          face_z_0 = elem_faces(num_elem,2);
          //face_z_1 = elem_faces(num_elem,2+dimension);

//     w = .5*(vitesse[face_z_0]+vitesse[face_z_1]);
          wl = vitesse[face_z_0];

          val_moy(i,2) += wl;
          val_moy(i,5) += wl*wl;
          val_moy(i,7) += u*wl;
          val_moy(i,8) += v*wl;
          val_moy(i,9) += sqrt(u*u+wl*wl);      //vitesse tangentielle pour calcul du frottement
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
  //const DoubleTab& xp = zone_VDF.xp();
  const Navier_Stokes_Turbulent& N_S_Turb  = ref_cast(Navier_Stokes_Turbulent,mon_equation.valeur());
  const DoubleTab& nu_t = N_S_Turb.viscosite_turbulente().valeurs();
  int nb_elems = zone_VDF.zone().nb_elem();
  int num_elem,i;
  // double y;

  for (num_elem=0; num_elem<nb_elems; num_elem++)
    {

      //y=xp(num_elem,1);

      i= Tab_recap(num_elem,2);

      val_moy(i,12) += nu_t[num_elem];
    }
}

void Traitement_particulier_NS_canal_VDF::calculer_moyenne_spatiale_Temp(DoubleTab& val_moy) const
{
  //modified AT 5/06/09
  const Probleme_base& pb = mon_equation->probleme();
  const Equation_base& eqn_th = pb.equation(1);
  const DoubleTab& diffusivite_turb =ref_cast(Modele_turbulence_scal_base,eqn_th.get_modele(TURBULENCE).valeur()).diffusivite_turbulente().valeurs();

  calculer_Temp(val_moy, diffusivite_turb);
  /*
    else {
    Cerr << " on doit forcement avoir un Pb_Thermohydraulique_Turbulent ou un Pb_Thermohydraulique_Turbulent_QC ou un Pb_Thermohydraulique_Concentration_Turbulent" << finl;
    exit();
  	}
  */
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

  double y,u,v,v_nm1,v_np1,wl,T,lambda_sm,rho,Tb,Th,dyb,dyh; //modified YB 29/06/09
  double lambdah,lambdab,lambdae,lambda_smh,lambda_smb; //modified AT 03/09/09

  int nb_elems = zone_VDF.zone().nb_elem();
  int num_elem,i;
  int face_x_0,face_x_1,face_y_0,face_y_1,face_z_0,face_z_1;

  const Fluide_Quasi_Compressible& le_fluide = ref_cast(Fluide_Quasi_Compressible,mon_equation->milieu());//modif AT 03/09/09
  const DoubleTab& tab_rho_elem = le_fluide.masse_volumique();
  //const DoubleTab& tab_rho_face = le_fluide.rho_face_n(); //modif AT 03/09/09
  const DoubleVect& lambda = le_fluide.conductivite(); //modif YB 30/06/09
  int taille_rho=tab_rho_elem.dimension(0);
  //  const Navier_Stokes_Turbulent& N_S_Turb  = ref_cast(Navier_Stokes_Turbulent,mon_equation.valeur());	//modif YB 30/06/09
  //const DoubleTab& nu_t = N_S_Turb.viscosite_turbulente().valeurs();	//modif YB 30/06/09

  wl=0.;

  for (num_elem=0; num_elem<nb_elems; num_elem++)
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
      //nuth = nu_t[elem_haut];
      //nutb = nu_t[elem_bas];
      lambda_smh=diffusivite_turb[elem_haut];
      lambda_sm = diffusivite_turb[num_elem];
      lambda_smb=diffusivite_turb[elem_bas];


//      if (num_elem==elem_haut||num_elem==elem_bas)
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

          wl = .5*(vitesse[face_z_0]+vitesse[face_z_1]);

          val_moy(i,17) += wl*T;
          val_moy(i,22) += rho*wl*T;
          val_moy(i,25) += rho*wl;
          val_moy(i,28) += rho*wl*wl;
          val_moy(i,30) += rho*u*wl;
          val_moy(i,31) += rho*v*wl;
        }
      val_moy(i,18) += lambda_sm;

      val_moy(i,32) += 0; //dyh/(dyb*(dyb+dyh))*(T-Tb)+dyb/(dyh*(dyb+dyh))*(Th-T);		//modif YB 29/06/09
      //modif AT 10/05/2010
      val_moy(i,33) += -0.25*((lambdae+lambdah)*(Th-T)/dyh+(lambdae+lambdab)*(T-Tb)/dyb);
      val_moy(i,34) += -0.25*((lambda_sm+lambda_smh)*(Th-T)/dyh+(lambda_sm+lambda_smb)*(T-Tb)/dyb);


    }

}
//fin modif AT 5/06/09

