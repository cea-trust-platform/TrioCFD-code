/****************************************************************************
* Copyright (c) 2015 - 2016, CEA
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
/////////////////////////////////////////////////////////////////////////////
//
// File      : Trait_part_NS_surface_VDF.cpp
// Directory : $NEW_ALGO_QC_ROOT/src/modif_QC/Trait_part_NS_surface
//
/////////////////////////////////////////////////////////////////////////////

#include <Trait_part_NS_surface_VDF.h>
#include <Zone_VDF.h>
#include <Pave.h>

#include <Fluide_Incompressible.h>
#include <Navier_Stokes_Turbulent.h>
#include <Loi_Etat_GP.h>
#include <Equation.h>
#include <DoubleTrav.h>
#include <Equation_base.h>
#include <Fluide_Quasi_Compressible.h>
#include <Convection_Diffusion_Turbulent.h>
#include <Convection_Diffusion_Chaleur_Turbulent_QC.h>
#include <Convection_Diffusion_Concentration_Turbulent.h>
#include <DoubleTab.h>
#include <Convection_Diffusion_Temperature_Turbulent.h>
#include <Navier_Stokes_QC.h>
#include <Navier_Stokes_Turbulent_QC.h>
#include <communications.h>
Implemente_instanciable(Traitement_particulier_NS_surface_VDF,"Traitement_particulier_NS_surface_VDF",Traitement_particulier_NS_surface);

#define CERR(x) \
Cerr << " oooo " << x << finl;
/*! @brief 
 *
 * @param (Sortie& is) un flot de sortie 
 * @return (Sortie&) le flot de sortie modifie 
 */
Sortie& Traitement_particulier_NS_surface_VDF::printOn(Sortie& is) const
{
  return is;
}

/*! @brief 
 *
 * @param (Entree& is) un flot d'entree 
 * @return (Entree&) le flot d'entree modifie 
 */
Entree& Traitement_particulier_NS_surface_VDF::readOn(Entree& is)
{
  return is;
}

Entree& Traitement_particulier_NS_surface_VDF::lire(Entree& is)
{
  return Traitement_particulier_NS_surface::lire(is);
}

void Traitement_particulier_NS_surface_VDF::remplir_XYZ(DoubleVect& Xc,DoubleVect& Yc, DoubleVect& Zc, int& N_x, int& N_y, int& N_z,IntTab& Tab_rec) const
{
  const Zone_dis_base& zdisbase=mon_equation->inconnue().zone_dis_base();
  const Zone_VF& zone_VF=ref_cast(Zone_VF, zdisbase);
  const Zone_VDF& zone_VDF=ref_cast(Zone_VDF, zdisbase);
  const DoubleTab& xp = zone_VF.xp();
  const IntTab& elem_faces = zone_VDF.elem_faces();

  int face; //recepteur des faces
  int nb_elems = zone_VF.zone().nb_elem_tot(); // pour pouvoir faire les derivation le tableau doit comporter les adresses des fictifs.
  int num_elem,i,indicx,indicy,indicz;
  double x,y,z;
  int a=0;

  indicx = 1;
  Xc.resize(2);
  Xc[0]=xp(0,0);
  indicy = 1;
  Yc.resize(2);
  Yc[0]=xp(0,1);
  indicz = 1;
  Zc.resize(2);
  Zc[0]=xp(0,2);

  Tab_rec.resize(nb_elems,6);

  for (num_elem=0; num_elem<nb_elems; num_elem++) //Remplissage du tableau des voisin et des vecteur X et Z
    {
      z = xp(num_elem,2);
      y = xp(num_elem,1);
      x = xp(num_elem,0);

      for(i=0; i<indicx; i++)
        {
          a = est_egal(x,Xc[i]);
          if(a)
            {
              break ;
            }
        }

      if(!a)
        {
          indicx ++;
          Xc[i]=x;
          Xc.resize(Xc.size()+1);
        }

      for (i=0; i<indicy; i++) // boucle sur les surfaces Y
        {
          a = est_egal(y,Yc[i]);
          if(a)
            {
              break;
            }
        }
      if(!a)
        {
          indicy ++;
          Yc[i]=y;
          Yc.resize(Yc.size()+1);
        }

      for (i=0; i<indicz; i++) // boucle sur les surfaces Z
        {
          a = est_egal(z,Zc[i]);
          if(a)
            {
              break;
            }
        }
      if(!a)
        {
          indicz ++;
          Zc[i]=z;
          Zc.resize(Zc.size()+1);
        }
      // on remplit le tableau des voisins
      // elem_voisin(numerode de l'elem, numero de la face, 0 dans la direction oposee du repere ou 1 dans la direction au repere)
      // si le voisin n'existe pas la valeur retourner est -1
      //faces X
      face                   =  elem_faces(num_elem,0);
      Tab_rec(num_elem,0)  =  zone_VDF.elem_voisin(num_elem,face,0);

      face                   =  elem_faces(num_elem,0+dimension);
      Tab_rec(num_elem,1)  =  zone_VDF.elem_voisin(num_elem,face,1);

      //faces Y
      face                   =  elem_faces(num_elem,1); //face inferieure
      Tab_rec(num_elem,2)  =  zone_VDF.elem_voisin(num_elem,face,0);

      face                   =  elem_faces(num_elem,1+dimension); //face superieure
      Tab_rec(num_elem,3)  =  zone_VDF.elem_voisin(num_elem,face,1);

      //face Z
      face                   =  elem_faces(num_elem,2);
      Tab_rec(num_elem,4)  =  zone_VDF.elem_voisin(num_elem,face,0);

      face                   =  elem_faces(num_elem,2+dimension);
      Tab_rec(num_elem,5)  =  zone_VDF.elem_voisin(num_elem,face,1);

    } //fin  for (num_elem=0;num_elem<nb_elems;num_elem++)
  Xc.resize(indicx); // on minimise la taille des tableaux
  Zc.resize(indicz);
  Yc.resize(indicy);
  N_z = indicz;
  N_y = indicy;
  N_x = indicx;
}

// fonction qui cree un tableau de reference permettant dans chaque proc de calculer des valeurs moyennes associe aux elements.
// elle est appelee apres remplir XYZ pour profiter de la construction du vecteur Y

void Traitement_particulier_NS_surface_VDF::recuperation_grandeurs(DoubleTab& val_post) const
{
// cerr.setf(ios::scientific); /// pour le debug a retirer quand tout marche

  const DoubleTab&            		temperature	= Temp.valeur().valeurs();
  const Zone_dis_base&        		zdisbase	= mon_equation->inconnue().zone_dis_base();
  const Zone_VDF&            	 	zone_VDF	= ref_cast(Zone_VDF, zdisbase);
  const IntTab&              	  	elem_faces 	= zone_VDF.elem_faces();
  const Fluide_Incompressible&  	le_fluide 	= ref_cast(Fluide_Incompressible,mon_equation->milieu());
  const DoubleTab& 					vitesse 	= mon_equation->inconnue().valeurs();
  const DoubleTab& 					pression 	= mon_equation->pression().valeurs();
  const DoubleTab& 					visco_dyn 	= le_fluide.viscosite_dynamique();
  const DoubleTab& 					tab_rho_elem= le_fluide.masse_volumique();
  const DoubleTab& 					lambda 		= le_fluide.conductivite();
  const Navier_Stokes_Turbulent& 	N_S_Turb 	= ref_cast(Navier_Stokes_Turbulent,mon_equation.valeur());
  const DoubleTab& 					mu_t		= N_S_Turb.viscosite_turbulente().valeurs();

  int dimension  = Objet_U::dimension;
  int nb_elems   = zone_VDF.zone().nb_elem();
  int num_elem;  // pos,j;
  int face_x_0, face_x_1, face_y_0, face_y_1, face_z_0;

  double U1=0; // vitesses au centre des mailles
  double U2=0; // vitesses au centre des mailles
  double U3=0; // vitesses au centre des mailles
  double T=0;
  double rho=0;
  double P=0;
  double L=0;
  DoubleTab viscosite_effec(visco_dyn); //  mu_eff = mu +mu_sm

  //Se servir seulement des valpost pour envoyer les valeurs necessaires au calcul de la moyenne vers trait_part_surface.cpp
  // la moyenne se fera dans ce dernier fichier a la ligne 495

  for (num_elem=0; num_elem<nb_elems; num_elem++)
    {
      face_x_0 = elem_faces(num_elem,0);
      face_x_1 = elem_faces(num_elem,dimension);
      face_y_0 = elem_faces(num_elem,1);
      face_y_1 = elem_faces(num_elem,1+dimension);
      face_z_0 = elem_faces(num_elem,2);
      // face_z_1 = elem_faces(num_elem,2+dimension);

      // FA 1/09/11 : pour eviter de moyenner localement la vitesse u et w
      // on est obliger de ramener la vitesse au centre des elements.
      U1 = .5*(vitesse[face_x_0]+vitesse[face_x_1]);
      U2 = .5*(vitesse[face_y_0]+vitesse[face_y_1]);
      U3 = vitesse[face_z_0]; // on peu se permettre d'utilise la vitesse d'une seule face grace a la periodicite
      T  = temperature[num_elem];
      rho= tab_rho_elem[num_elem];
      P  = pression[num_elem];
      L  = lambda[num_elem];

      // Grandeurs directes
      val_post(num_elem,0)   =  U1;											// U moyen
      val_post(num_elem,1)   =  U2; 											// V moyen
      val_post(num_elem,2)   =  U3; 											// W moyen
      val_post(num_elem,3)   =  T;  											// T moyen
      val_post(num_elem,4)   =  P; 	  										// P moyen
      val_post(num_elem,5)   =  rho; 											// rho moyen
      val_post(num_elem,6)   =  L; 	  										// lambda moyen
      val_post(num_elem,7)   =  visco_dyn[num_elem];	 						// mu_phys moyen
      val_post(num_elem,8)   =  mu_t[num_elem];  							// mu_turb
      val_post(num_elem,9)   =  visco_dyn[num_elem]/rho;   					// viscosite cinematique mu/rho moyenne

      // Grandeurs au carre pour le calcul des rms
      val_post(num_elem,10)  =  U1*U1;             							// Urms
      val_post(num_elem,11)  =  U2*U2;                      					// Vrms
      val_post(num_elem,12)  =  U3*U3;                       					// Wrms
      val_post(num_elem,13)  =  T*T;  										// Trms
      val_post(num_elem,14)  =  P*P; 											// Prms
      val_post(num_elem,15)  =  rho*rho;										// rhorms
      val_post(num_elem,16)  =  L*L;											// lambdarms
      val_post(num_elem,17)  =  visco_dyn[num_elem]*visco_dyn[num_elem];		// mu_physrms
      val_post(num_elem,18)  =  mu_t[num_elem]*mu_t[num_elem];				// mu_turbrms
      val_post(num_elem,19)  =  val_post(num_elem,16)*val_post(num_elem,16);  // viscosite cinematiquerms

      // calcul des correlation
      val_post(num_elem,20)  =  U1*U2;                          				// correlation U.V
      val_post(num_elem,21)  =  U1*U3;    									// correlation U.W utile ?
      val_post(num_elem,22)  =  U2*U3;    									// correlation V.W utile ?
      val_post(num_elem,23)  =  U1*rho;                   					// correlation U.rho
      val_post(num_elem,24)  =  U2*rho;                   					// correlation V.rho
      val_post(num_elem,25)  =  U1*T;   										// correlation U.T
      val_post(num_elem,26)  =  U2*T;  										// correlation V.T
      val_post(num_elem,27)  =  T*rho;										// correlation rho.T

    }

}
