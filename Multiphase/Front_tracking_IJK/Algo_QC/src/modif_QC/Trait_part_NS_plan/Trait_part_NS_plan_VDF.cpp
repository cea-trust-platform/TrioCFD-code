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
// File      : Trait_part_NS_plan_VDF.cpp
// Directory : $NEW_ALGO_QC_ROOT/src/modif_QC/Trait_part_NS_plan
//
/////////////////////////////////////////////////////////////////////////////

#include <Trait_part_NS_plan_VDF.h>
#include <Domaine_VDF.h>
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
#include <Convection_Diffusion_Temperature_Turbulent.h>
#include <Probleme_base.h>
#include <DoubleTab.h>

#include <Navier_Stokes_QC.h>
#include <Navier_Stokes_Turbulent_QC.h>
#include <communications.h>
#include <IJK_Field.h>

Implemente_instanciable(Traitement_particulier_NS_plan_VDF,"Traitement_particulier_NS_plan_VDF",Traitement_particulier_NS_plan);

#define CERR(x)					\
  Cerr << " oooo " << x << finl;
/*! @brief
 *
 * @param (Sortie& is) un flot de sortie
 * @return (Sortie&) le flot de sortie modifie
 */
Sortie& Traitement_particulier_NS_plan_VDF::printOn(Sortie& is) const
{
  return is;
}


/*! @brief
 *
 * @param (Entree& is) un flot d'entree
 * @return (Entree&) le flot d'entree modifie
 */
Entree& Traitement_particulier_NS_plan_VDF::readOn(Entree& is)
{
  return is;
}

Entree& Traitement_particulier_NS_plan_VDF::lire(Entree& is)
{
  return Traitement_particulier_NS_plan::lire(is);
}

void Traitement_particulier_NS_plan_VDF::calculer_valeur_spatiale_vitesse_rho_mu(DoubleTab& val_post,DoubleTab& Moyennes,double dpth)
{
  cerr.setf(ios::scientific);
  const DoubleTab&            						temperature	= ref_champ_temperature_.valeur().valeurs();
  const Domaine_dis_base&        						zdisbase	= mon_equation->inconnue().domaine_dis_base();
  const Domaine_VDF&            	 					domaine_VDF	= ref_cast(Domaine_VDF, zdisbase);
  const IntTab&              	  					elem_faces 	= domaine_VDF.elem_faces();
  const Fluide_Incompressible&  					le_fluide 	= ref_cast(Fluide_Incompressible,mon_equation->milieu());
  const DoubleTab& 									vitesse 	= mon_equation->inconnue().valeurs();
  const DoubleTab& 									pression 	= mon_equation->pression().valeurs();
  const DoubleTab& 									visco_dyn 	= le_fluide.viscosite_dynamique();
  const DoubleTab& 									tab_rho_elem= le_fluide.masse_volumique();
  const DoubleTab& 									lambda 		= le_fluide.conductivite();
  const Fluide_Quasi_Compressible& 					fluide_QC	= ref_cast(Fluide_Quasi_Compressible,le_fluide);
  const Loi_Etat_GP& 								loi_ 		= ref_cast(Loi_Etat_GP,fluide_QC.loi_etat().valeur());
  const Navier_Stokes_Turbulent& 					N_S_Turb 	= ref_cast(Navier_Stokes_Turbulent,mon_equation.valeur());
  const DoubleTab& 									nu_t		= N_S_Turb.viscosite_turbulente().valeurs();
  const DoubleTab&     								xv        	= domaine_VDF.xv();
  const Probleme_base& 								pb 			= mon_equation->probleme();
  const Equation_base& 								eqn_th 		= pb.equation(1);
  const DoubleTab& cond_turb =ref_cast(Modele_turbulence_scal_base,eqn_th.get_modele(TURBULENCE).valeur()).diffusivite_turbulente().valeurs();
  const DoubleTab& 								    tab_rho_face = fluide_QC.rho_face_n();// rho discretisee aux faces.

  double CP      = loi_.Cp();
  double tampon  = 0, tampon2 =0 ; // variables de travail
  double f_p1 = 0 , f_m1 =0 , f =0;
  double unsurcp = 1/CP;

  int dimension  = Objet_U::dimension;
  int nb_elem= domaine_VDF.domaine().nb_elem(); // nombre total d'elements (reel + fict)
  int nb_elem_tot= domaine_VDF.domaine().nb_elem_tot(); // nombre total d'elements (reel + fict)
  int num_elem = 0,i,j;
  int REFY = 0 ;

  int face_p_0, face_p_1, face_m_0, face_m_1;
  //  int Nxt = X_tot.size() ;
  //  int Nzt = Z_tot.size() ;
  int p = 0 ; // variable qui sera utiliser pour les boucles sur les processeurs.
  /// ici je vais ecrire les numeros d'indices.
  int Tq = 3, Tn = 6 , Ti = 9 ,  Tg = 12 , Tp = 15, TTq = 19 , TTi = 20 , Dq = 21 , Dg = 24, TTd = 27 , K = 28 ;

  // on declare les variables qui vont stocker les derivees.
  double   divu = 0;

  DoubleTrav dnudxi(3) ;
  DoubleTrav dcdxi(3);
  // ici on va stocker les fluctuations de vitesse sur les faces.
  // ca ne coute rien de les recalculees.
  DoubleTrav u1(3);
  DoubleTrav u0(3);
  DoubleTrav v1(3);
  DoubleTrav v0(3);

  DoubleTrav usdelta(3);
  DoubleTrav stock(3);
  DoubleTab distance(6) ;
  distance = -1 ;
  IntTab face(6);
  face = -1 ;

  ////////////////////////
  DoubleTab sources(vitesse);
  sources =0;
  mon_equation.valeur().sources().ajouter(sources);
  DoubleTrav moyennes_par_proc(Moyennes); // tableau qui copie la structure de moyennes_par_proc.
  //

  DoubleTab U(nb_elem_tot,3);  // vitesses au centre des mailles
  DoubleTab V(nb_elem_tot,3);  // vitesses au centre des mailles * rho^(1/2)
  DoubleTab duidxi(nb_elem_tot,3); // derivees de vitesse locales de fluctuation
  DoubleTab dvidxi(nb_elem_tot,3); // derivees de vitesse locales globales
  DoubleTrav DivU(nb_elem_tot,3); // derivee de vitesse locales globales
  DoubleTab b(tab_rho_elem);
  DoubleTab C(temperature);
  DoubleTab viscosite_effec(visco_dyn); // b = rho^(-1/2) C = b^(-1) * T  mu_eff = mu +mu_sm
  DoubleTab conduc_effec(lambda) ;


  /*	// on va tester ici.
  	///////////////////////////
  	if (!ival)
  	{
  		const Domaine_dis_base& zdis2base=mon_equation->inconnue().domaine_dis_base();
  		const Domaine_VF& domaine_VF=ref_cast(Domaine_VF, zdis2base);
  		const DoubleTab& xp = domaine_VF.xp();
  		DoubleTab test_para(vitesse);
  		DoubleTab test_xp(xp);
  	if (me())
  	{
  		envoyer(test_para, Process::me(), 0, 0);
  		envoyer(test_para2, Process::me(), 0, 0);
  		envoyer(test_xp, Process::me(), 0, 0);
  	}
  	if(je_suis_maitre())
      {
          for(int p=0;p<Process::nproc();p++)
          {
  			if (p)
  			{
  				recevoir(test_para,p,0,0);
  				envoyer(test_para2,p, 0, 0);
  				recevoir(test_xp,p,0,0);
  			}
  			for (int i=0; i<nb_elem_vdf; i++)
  				Cerr << "p " << p << " i " <<  test_xp(i,1) << " " << test_para(i)-test_para2(i,0) << endl;
  	    }
      }
  	////////////////////////////////////////////// ici c'est bon ....
  	}
  	*/

  // on dissocie en 2 boucles car il faut que toutes les operations soit finies sur les maille fictives pour post traiter les derivees d'ordre sup.
  for (num_elem=0; num_elem<nb_elem_tot; num_elem++)
    {
      REFY = Ref_Y(num_elem) ;

      for (i = 0 ; i < 3 ; i++)
        {
          face(i)			 = elem_faces(num_elem,i);
          face(i+dimension)= elem_faces(num_elem,i+dimension);
          if ( (!num_elem) || (i == 1) ) // on se sert du fait que la taille de maille ne varie pas sur x et z pour ne faire qu'une fois cette division.
            {
              usdelta(i) = 1. / ( xv(face(i+dimension),i) - xv(face(i),i) ) ;
            }

          if (num_elem >= nb_elem ) // on calcule ce qui n'a pas ete calculer lors de la moyenne.
            {
              Rrhouf[face(i)] 			= sqrt( tab_rho_face[face(i)]			) * vitesse[face(i)] ;
              Rrhouf[face(i+dimension)] 	= sqrt( tab_rho_face[face(i+dimension)] ) * vitesse[face(i+dimension)] ;
              racinederho[num_elem] 		= sqrt( tab_rho_elem[num_elem]			);
            }
          // pour calculer la fluctuation au centre des element ou la derivee; on commence par la calculer aux face puis on interpole au centre
          // on est obliger de ramener la vitesse au centre des elements.
          // mais on ne le fait qu'au dernier moment.
          tampon  			= vitesse[face(i+dimension)	] - Moyennes(REFY,11+i*2);
          tampon2 	 		= vitesse[face(i)			] - Moyennes(REFY,10+i*2);
          U(num_elem,i) 		= ( tampon + tampon2 ) * 0.5 ;
          duidxi(num_elem,i)	= ( tampon - tampon2 ) * usdelta(i);

          tampon  			= Rrhouf[face(i+dimension)	] - Moyennes(REFY,38+i*2);
          tampon2 			= Rrhouf[face(i) 			] - Moyennes(REFY,37+i*2);
          V(num_elem,i) 		= ( tampon + tampon2 ) * 0.5;
          dvidxi(num_elem,i)	= ( tampon - tampon2 ) * usdelta(i);

          DivU(num_elem,i) 	= (vitesse[face(i+dimension)] - vitesse[face(i)] ) * usdelta(i) ;
        }
      if (num_elem >= nb_elem ) // on calcule ce qui n'a pas ete calculer lors de la moyenne.
        racinederho[num_elem] = sqrt(tab_rho_elem[num_elem]);


      C[num_elem] = racinederho[num_elem]*temperature(num_elem)-Moyennes(REFY,8); // idem pour C'
      b[num_elem] = 1./ racinederho[num_elem];
      viscosite_effec[num_elem]+=0*nu_t[num_elem]; // on ne travaille pas en viscosite complete car cela pose trop de probleme.
      conduc_effec[num_elem]   +=0* cond_turb[num_elem];// on ne travaille pas en conductivite complete car cela pose trop de probleme.

    }
  for ( p =0 ; p<distance.size() ; p++) // on calcule et stock les distances.
    distance[p] = distance_elem(0,p); // On ne s'interesse ici qu'au direction x et z.
  // a savoir que les distance impaires sont positives et les paires negatives

  // on calcule toutes les donnees sur les mailles reelles et ont les envois dans Val_post

  // Pour chaque plan de maillage en k, le plan doit-il etre postraite ou non ?
  ArrOfInt plan_a_traiter(ijk_splitting_.get_grid_geometry().get_nb_elem_tot(DIRECTION_K));
  {
    for (int plan = 0; plan < indices_planes_.size_array(); plan++)
      {
        plan_a_traiter[indices_planes_[plan]] = 1; // marque le plan a traiter
      }
  }
  const int nb_elem_local = Ref_Y.size();
  for (num_elem = 0; num_elem < nb_elem_local; num_elem++)
    {

      // Ne pas traiter les elements qui ne sont pas sur un plan a post-traiter:
      if (! plan_a_traiter[Ref_Y[num_elem]])
        continue;

      // for(plan=0;plan< indices_planes_.size_array();plan++)
      // for(cordx=0;cordx<nb_elem_i;cordx++)
      // for(cordz=0;cordz<nb_elem_j;cordz++) 	// boucles sur les element a post traites
      //  {
      //   if( Tab_post(plan,cordx,cordz) +1 ) // faux si on ne post traite pas.
      //	      {
      // num_elem = Tab_post(plan,cordx,cordz) ;
      if ( ( Tab_recap(num_elem,2)!=-1 ) && (Tab_recap(num_elem,3)!=-1) )// // la methode ne calcule pas les valeurs aux bords
        {
          // on calcule les derivees d'ordre 1 (hors div U )
          //derivee de vitesse non locale
          REFY = Ref_Y(num_elem) ;

          for (i=0 ; i < 6 ; i++)
            face(i)	 = elem_faces(num_elem,i);

          for (i=0 ; i < 3 ; i++ ) // calcul des u' v' ux faces
            {
              u0(i) = vitesse[face(i)]		  -Moyennes(REFY,10+2*i) ;
              u1(i) = vitesse[face(i+dimension)]-Moyennes(REFY,11+2*i) ;
              v0(i) = Rrhouf [face(i)]		  -Moyennes(REFY,37+2*i) ;
              v1(i) = Rrhouf [face(i+dimension)]-Moyennes(REFY,38+2*i) ;
            }

          usdelta(1) = 1. /( xv(face(1+dimension),1) - xv(face(1),1) )  ;
          distance[2] = distance_elem(num_elem,2); // on evalue uniquement  les distances sur y
          distance[3] = distance_elem(num_elem,3); // puisqu'on ne change pas sur x et z

          divu = 0 ;
          for(i=0 ; i< 3 ; i++) // on calcule divu des fluctutations. !! DivU != divu
            divu += duidxi(num_elem,i);

          // on prepare le calcul des derivees
          for(i=0 ; i<3 ; i++)
            dnudxi(i)	 = derivee_un_elem(viscosite_effec,num_elem,i) ;

          //L'objectif est de code l'equation suivante :
          //dv'i/dt + u'2d<Vi>/dx2 + <Uj>dvi/dxj + du'iv'j/dxj + 0.5 [v'id<U2>/dx2 + <Vi>du'j/dxj  - (vi'du'j/dxj') ]
          //  = -b dPdyn/dxi + mu b d2 u'i/dxj2 + mu b'd2<Ui>/dxj2 + b dmu /dxj * dUi/dxj + b dmu/dxj * dUj/dxi + 1/3 b mu*d2Uj/dxi*dxj -2b/3 dmu /dxi*dUj/dxj.

          // dv'i/dt et u'2d<Vi>/dx2 ne sont pas code car ils resultes directement des TF de v'i
          for (i = 0 ; i < 3 ; i++)
            {
              val_post(num_elem,i)   =  V(num_elem,i);

              // on va calculer simutanement dv'iu'j/dxj et <Uj>dv'i/dxj pour le deuxieme on calcule il suffit de calculer  j = 1 et j = 2 car <U3> = 0 ;
              for( j = 0 ; j<3 ; j++)
                {
                  if(i==j)
                    {
                      tampon = ( v1(i) * u1(j) - v0(i) *u0(j) ) * usdelta(j) ;  //du'iv'i/dxi
                      if ( j !=2 )
                        {
                          tampon2 = Moyennes(REFY,i)*dvidxi(num_elem,i); //<Uj>dv'i/dxj
                        }
                    }
                  else
                    {
                      f_p1  =  ( V(Tab_recap(num_elem,j*2+1),i)-V(num_elem,i) )  / distance[j*2+1]   ;
                      f_m1  =  ( V(Tab_recap(num_elem,j*2  ),i)-V(num_elem,i) )  / distance[j*2] ;

                      tampon = f_p1 * u1(j) + f_m1 * u0(j) ;
                      tampon *= 0.5 ;
                      tampon += duidxi(num_elem,j) * V(num_elem,i)  ;
                      if ( j!= 2 )
                        {
                          tampon2 = f_p1 * Moyennes(REFY, 11+2*j) + f_m1 * Moyennes(REFY, 10+2*j) ;
                          tampon2 *= 0.5 ;
                        }
                    }
                  if 	(j==1)
                    {
                      val_post(num_elem,Ti+i) += tampon ;
                      val_post(num_elem,Tg+i) += tampon2;
                    }
                  else
                    {
                      val_post(num_elem,Tq+i) += tampon ;
                      if (j !=2 )
                        {
                          val_post(num_elem,Tn+i) += tampon2 ;
                        }
                    }
                }

              // 0.5 v'id<U2>/dx2
              tampon = 0.5 * V(num_elem,i) *( Moyennes(REFY,13) - Moyennes(REFY,12) ) * usdelta(1)  ;
              val_post(num_elem,Tg+i) += tampon ;
              // 0.5 <Vi>du'j/dxj
              val_post(num_elem,Tg+i) += 0.5 * divu * Moyennes(REFY,3+i);
              // -  0.5 (vi'du'j/dxj')
              val_post(num_elem,Tg+i) -= 0.5 * divu * V(num_elem,i);
              /// -b dPdyn/dxi
              val_post(num_elem,Tp+i)  = (-1) * b[num_elem] * derivee_un_elem(pression,num_elem,i) ;
              /// mu b d2u'i/dxj2 ;
              for ( j = 0 ; j < 3 ; j++ )
                {
                  if (i == j )
                    {
                      tampon  = duidxi(Tab_recap(num_elem,2*j+1),i) ;
                      tampon -= duidxi(Tab_recap(num_elem,2*j  ),i) ;
                      tampon /= distance[j*2+1]-distance[j*2]; // donne la bonne distance mais la disparite de maille dans le cas de j = 2 n'est pas geree avec le decentrement
                    }
                  else
                    {
                      tampon  = derivee_deux_vitesse_nonco(U,i,num_elem,j); // methode qui gere le decentrement.
                    }
                  val_post(num_elem,Dq+i) += tampon ;
                }
              val_post(num_elem,Dq+i) *= b[num_elem]*viscosite_effec[num_elem];

              /// mu b d2<Ui>/dxj2 ;
              for ( j = 1 ; j  == 1 ; j++ ) // pas besoin de calculer la derivee si on est pas dans la direction y car sinon = 0 ;
                {
                  if (i == j )
                    {
                      face_m_0 = elem_faces(Tab_recap(num_elem,2*j  ) ,i) ; // dans la direction de la derivee  (j)
                      face_m_1 = elem_faces(Tab_recap(num_elem,2*j  ) ,i+dimension) ; // on va chercher la composante de vitesse i
                      face_p_0 = elem_faces(Tab_recap(num_elem,2*j+1) ,i) ;
                      face_p_1 = elem_faces(Tab_recap(num_elem,2*j+1) ,i+dimension);


                      f_p1  = ( Moyennes(REFY+1,13) - Moyennes(REFY,13) )  / ( xv(face_p_1,j) - xv(face_p_0,j) ) ;// d<U2>/dx2
                      f_m1  = ( Moyennes(REFY,12) - Moyennes(REFY-1,12) )  / ( xv(face_m_1,j) - xv(face_m_0,j) ) ;

                      tampon  = distance(2*j+1) ;
                      tampon2 = fabs(distance(2*j)) ;

                      tampon  = 2*(f_p1 * tampon - f_m1*tampon2 )/( (tampon + tampon2)* (tampon + tampon2) ) ; // on pondere par la distance au centre de la maille.
                    }
                  else
                    {
                      tampon = (Moyennes(REFY+1,i) - Moyennes(REFY,i)) * fabs(distance(3)) + (Moyennes(REFY-1,i) - Moyennes(REFY,i)) * fabs(distance(2));
                      tampon /= fabs(distance(3) * distance(2) ) * ( fabs(distance(3)) + fabs(distance(2)) ) * 0.5 ;
                    }
                  if ( i!=2)
                    {
                      val_post(num_elem,Dg+i) += tampon ;  // <U3> est cense faire 0
                    }
                }
              val_post(num_elem,Dg+i) *= viscosite_effec[num_elem];


              // on va calculer simultanement b dmu /dxj * dUi/dxj et b dmu /dxj * dUj/dxi comme bdmu /dxj ( dUi/dxj + dUj/dxi)
              for ( j = 0 ; j < 3 ; j++ )
                {
                  if (i == j )
                    {
                      tampon  =  2 * DivU(num_elem,i);
                    }
                  else
                    {
                      if ( i - j < 0 ) // si on est dans ce cas on n'a pas encore calcule cette valeur sinon elle est stockee.
                        {
                          face_m_0 = elem_faces(Tab_recap(num_elem,2*j  ) ,i) ; // dans la direction de la derivee  (j)
                          face_m_1 = elem_faces(Tab_recap(num_elem,2*j  ) ,i+dimension) ; // on va chercher la composante de vitesse i
                          face_p_0 = elem_faces(Tab_recap(num_elem,2*j+1) ,i) ;
                          face_p_1 = elem_faces(Tab_recap(num_elem,2*j+1) ,i+dimension);

                          f_p1     =  vitesse[face_p_1] + vitesse[face_p_0];
                          f_m1     =  vitesse[face_m_1] + vitesse[face_m_0];
                          f        =  vitesse[face(i) ] + vitesse[face(i+dimension)] ;  // le calcul de la moyenne n'est pas fait car il s'annule avec la division par deux de la derivee.


                          tampon    = (f_m1-f) /distance(2*j) +  (f_p1 - f) / distance(2*j+1) ;
                          tampon   /= distance(2*j+1)- distance(2*j) ;// on calcule la derivee tout en tenant compte des changement de taille de la maille

                          face_m_0 = elem_faces(Tab_recap(num_elem,2*i  ) ,j) ;
                          face_m_1 = elem_faces(Tab_recap(num_elem,2*i  ) ,j+dimension) ;
                          face_p_0 = elem_faces(Tab_recap(num_elem,2*i+1) ,j) ;
                          face_p_1 = elem_faces(Tab_recap(num_elem,2*i+1) ,j+dimension);

                          f_p1     =  vitesse[face_p_1] + vitesse[face_p_0];
                          f_m1     =  vitesse[face_m_1] + vitesse[face_m_0];
                          f        =  vitesse[face(j) ] + vitesse[face(j+dimension)] ;  // le calcul de la moyenne n'est pas fait car il s'annule avec la division par deux de la derivee.

                          tampon2   = (f_m1-f) /distance(2*i) +  (f_p1 - f) / distance(2*i+1) ;
                          tampon2  /= distance(2*i+1) - distance(2*i);// on calcule la derivee tout en tenant compte des changement de taille de la maille

                          tampon = tampon + tampon2 ;

                          stock(i+j-1) = tampon ; // permet de retrouver la bonne case. le -1 permet de commencer a la case 0.
                        }
                      else
                        {
                          tampon = stock(i+j-1) ;
                        }
                    }
                  val_post(num_elem,Dg+i) += tampon  * dnudxi(j);
                }

              //bmu /3*d2Uj/dxi*dxj
              f =0 ;
              f_p1 =0 ;
              f_m1=0;
              for ( j = 0 ; j < 3 ; j++ ) // calcule le Div U dans les mailles proches (direction i)
                {
                  f_p1     +=  DivU(Tab_recap(num_elem,2*i+1),j) ;
                  f_m1     +=  DivU(Tab_recap(num_elem,2*i  ),j) ;
                  f        +=  DivU(num_elem				   ,j) ; // le calcul de la moyenne n'est pas fait car il s'annule avec la division par deux de la derivee.
                }
              tampon   = (f_m1-f) /distance(2*i) +  (f_p1 - f) / distance(2*i+1) ;
              tampon  /=  distance(2*i+1)- distance(2*i);// on calcule la derivee tout en tenant compte des changement de taille de la maille
              tampon  *= 2 ;

              val_post(num_elem,Dg+i) +=tampon * viscosite_effec[num_elem] * 1/3.;
              // -2b/3 dmu /dxi*dUj/dxj
              tampon = 0 ;
              for ( j = 0 ; j < 3 ; j++ )
                tampon += DivU(num_elem,j);

              val_post(num_elem,Dg+i) += -2./3. * dnudxi(i) * tampon ;


              val_post(num_elem,Dg+i) *= b[num_elem] ; // on multiplie par b ici
            }

          for ( j = 0 ; j < 3 ; j++ )
            dcdxi(j)=derivee_un_elem(C,num_elem,j);

          // dc'/dt + <Uj>dc'/dxj + uj'dc'/dxj + C/2 DivU -b/cp[dpth/dt] =  b/cpd/dxj(ldT/dxj)
          val_post(num_elem,18)  = C[num_elem];  //c'

          for (j = 0 ; j <3 ; j++) //<Uj>dc'/dxj  et u'jdc'/dxj
            {
              f_p1   =  ( C(Tab_recap(num_elem,j*2+1))-C(num_elem) )  / distance[j*2+1]  ;
              f_m1   =  ( C(Tab_recap(num_elem,j*2  ))-C(num_elem) )  / distance[j*2]    ;

              tampon *= 0.5 * ( u1(j) * f_p1 + u0(j) * f_m1 ) ;

              if ( j !=2 )
                {
                  tampon2 = 0.5 * ( Moyennes(REFY, 11+2*j) * f_p1 + Moyennes(REFY, 10+2*j) * f_m1 ); //<Uj>dc'/dxj
                }
              else
                {
                  tampon2 = 0 ;
                }
              if (j == 1 )
                {
                  val_post(num_elem,TTd)  += tampon2 ; 	// <U2>dc'/dx2
                  val_post(num_elem,TTi)  += tampon  ;    // u'2 dc'/dx2
                }
              else
                {
                  val_post(num_elem,TTq)  += tampon + tampon2 ;// Ujdc'/dxj
                }
            }

          // u'2 d<C>/dx2
          tampon  = u1(1) *  ( Moyennes(REFY +1 ,8) - Moyennes(REFY,8) ) * usdelta(1) ;
          tampon -= u0(1) *  ( Moyennes(REFY -1 ,8) - Moyennes(REFY,8) ) * usdelta(1) ; // le signe - viens du sens de derivation.
          tampon *= 0.5 ;
          val_post(num_elem,TTi)  += tampon ;

          // CdUj/dxj
          tampon = 0;
          for (j = 0 ; j <3 ; j++)
            tampon += DivU(num_elem,j);

          val_post(num_elem,TTd) += racinederho[num_elem]*temperature(num_elem) * tampon ; // on ne distingue pas les contribution moyennes des fluctuantes.
          val_post(num_elem,TTd) += (b[num_elem] - Moyennes(REFY,6) ) * unsurcp * dpth ;

          //  b/Cp * d/dxj ( kdT/dxj )
          for ( j = 0 ; j < 3 ; j ++)
            {
              tampon   = (temperature[Tab_recap(num_elem,2*j+1)] - temperature[num_elem] )/ fabs(distance[2*j+1]) ;
              tampon  *=  conduc_effec[num_elem]+conduc_effec[Tab_recap(num_elem,2*j+1)] ; // calcul sur la face 1

              tampon2  = (temperature[num_elem] - temperature[Tab_recap(num_elem,2*j)]  )/ fabs(distance[2*j]) ;
              tampon2 *=  conduc_effec[num_elem]+conduc_effec[Tab_recap(num_elem,2*j)] ; // calcul sur la face 0

              tampon   = ( tampon - tampon2 ) * 0.5 * usdelta(j);// la moyenne des conduc est faite ici // on moyennes les faces ici.

              val_post(num_elem,K) += b[num_elem] * unsurcp * tampon;
            }

        }
      else // pour les bords ont met les valeurs a zero pour ne pas avoir de probleme.
        for (int pos = 0; pos < val_post.dimension(1) ; pos ++)
          val_post(num_elem,pos) = 0 ;
    }

}

double Traitement_particulier_NS_plan_VDF::derivee_un_elem(const DoubleTab& valeur,int element,int direction) const
{

  int    moins   = -1 ;
  int    plus    = -1 ;
  int    centre  =  0 ;

  if (direction == 0)
    {
      moins = 0;
      plus = 1;
    }
  else if(direction == 1)
    {
      moins = 2;
      plus = 3;
    }
  else if(direction == 2)
    {
      moins = 4;
      plus = 5;
    }
  else
    {
      Cerr << "Erreur dans le choix de la direction de derivation" << finl;
      Cerr << " x=0 y=1 z=2" << finl;
      exit();
    }

  double derivee;  // valeur retournee par la fonction

  double f_p1;     // valeur au premier point
  double f_m1;     // valeur au deuxieme point
  double f = 0;    // valeur dans la maille traitee

  double delta_p;  //distance au premier point
  double delta_m;  // distance au deuxieme point
  double alpha,unsalpha; // rapports des distances

  if(Tab_recap(element,plus)==-1) // si l'element plus n'existe pas alors on descentre le schema de cote de moins
    {
      f_p1      =  valeur[Tab_recap(Tab_recap(element,moins),moins)];
      f_m1      =  valeur[Tab_recap(element,moins)];
      f         =  valeur[element];

    }
  else if(Tab_recap(element,moins)==-1)// si l'element moins n'existe pas alors on decentre le schema vers plus
    {
      f_p1      =  valeur[Tab_recap(element,plus)];
      f_m1      =  valeur[Tab_recap(Tab_recap(element,plus),plus)];
      f         =  valeur[element];

    }
  else // si tous va bien
    {
      f_p1      =  valeur[Tab_recap(element,plus)];
      f_m1      =  valeur[Tab_recap(element,moins)];
      f 		  =  valeur[element];
      if( direction-1 ) centre++ ; // on indique le schemasi on est dans l'axe y on est forcement non isotrope
    }

  if(centre)
    {
      delta_p	  = distance_elem(element,plus) ; // distance a l'element proche

      derivee   = f_p1 - f_m1 ; //schema centre
      derivee   /= 2 * delta_p ;
    }
  else
    {
      delta_m   = distance_elem(element,moins); // la methode de calcul gere de decentrement et la periodicite
      delta_p   = distance_elem(element,plus);  // la methode de calcul gere de decentrement et la periodicite
      delta_m *= -1 ; // on prend en valeur absolue ;

      alpha 	  = delta_p/delta_m ; // rapport des distance
      unsalpha  = delta_m/delta_p ; //

      derivee   = (f_p1 -f ) * unsalpha  + (f - f_m1 ) * alpha ;
      derivee  /= delta_m + delta_p;

    }
  return derivee;
}

double Traitement_particulier_NS_plan_VDF::derivee_un_elem(DoubleTab& valeur,int element,int direction) const
{
  int    moins   = -1 ;
  int    plus    = -1 ;
  int    centre  =  0 ;

  if (direction == 0)
    {
      moins = 0;
      plus = 1;
    }
  else if(direction == 1)
    {
      moins = 2;
      plus = 3;
    }
  else if(direction == 2)
    {
      moins = 4;
      plus = 5;
    }
  else
    {
      Cerr << "Erreur dans le choix de la direction de derivation" << finl;
      Cerr << " x=0 y=1 z=2" << finl;
      exit();
    }

  double derivee;  // valeur retournee par la fonction

  double f_p1;     // valeur au premier point
  double f_m1;     // valeur au deuxieme point
  double f = 0;    // valeur dans la maille traitee

  double delta_p;  //distance au premier point
  double delta_m;  // distance au deuxieme point
  double alpha,unsalpha; // rapports des distances

  if(Tab_recap(element,plus)==-1) // si l'element plus n'existe pas alors on descentre le schema de cote de moins
    {
      f_p1      =  valeur[Tab_recap(Tab_recap(element,moins),moins)];
      f_m1      =  valeur[Tab_recap(element,moins)];
      f         =  valeur[element];

    }
  else if(Tab_recap(element,moins)==-1)// si l'element moins n'existe pas alors on decentre le schema vers plus
    {
      f_p1      =  valeur[Tab_recap(element,plus)];
      f_m1      =  valeur[Tab_recap(Tab_recap(element,plus),plus)];
      f         =  valeur[element];

    }
  else // si tous va bien
    {
      f_p1      =  valeur[Tab_recap(element,plus)];
      f_m1      =  valeur[Tab_recap(element,moins)];
      f 		  =  valeur[element];
      if( direction-1 ) centre++ ; // on indique le schemasi on est dans l'axe y on est forcement non isotrope
    }

  if(centre)
    {
      delta_m	  = distance_elem(element,plus) ; // distance a l'element proche

      derivee   = f_p1 - f_m1 ; //schema centre
      derivee  /= 2 * delta_m ;
    }
  else
    {
      delta_m   = distance_elem(element,moins); // la methode de calcul gere de decentrement et la periodicite
      delta_p   = distance_elem(element,plus);  // la methode de calcul gere de decentrement et la periodicite
      delta_m *= -1 ; // on prend en valeur absolue ;

      alpha 	  = delta_p/delta_m ; // rapport des distance
      unsalpha  = delta_m/delta_p ; //

      derivee   = (f_p1 -f ) * unsalpha  + (f - f_m1 )*alpha ;
      derivee  /= delta_m + delta_p;

    }
  return derivee;
}

double Traitement_particulier_NS_plan_VDF::derivee_un_elem(DoubleTab& vitesse,int composante, int element,int direction) const
{
  int    moins   = -1 ;
  int    plus    = -1 ;
  int    centre  =  0 ;

  if (direction == 0)
    {
      moins = 0;
      plus = 1;
    }
  else if(direction == 1)
    {
      moins = 2;
      plus = 3;
    }
  else if(direction == 2)
    {
      moins = 4;
      plus = 5;
    }
  else
    {
      Cerr << "Erreur dans le choix de la direction de derivation" << finl;
      Cerr << " x=0 y=1 z=2" << finl;
      exit();
    }

  double derivee;  // valeur retournee par la fonction

  double f_p1;     // vitesse au premier point
  double f_m1;     // vitesse au deuxieme point
  double f = 0;    // vitesse dans la maille traitee

  double delta_p;  //distance au premier point
  double delta_m;  // distance au deuxieme point
  double alpha,unsalpha; // rapports des distances

  if(Tab_recap(element,plus)==-1) // si l'element plus n'existe pas alors on descentre le schema de cote de moins
    {
      f_p1      =  vitesse(Tab_recap(Tab_recap(element,moins),moins),composante);
      f_m1      =  vitesse(Tab_recap(element,moins),composante);
      f         =  vitesse(element,composante);

    }
  else if(Tab_recap(element,moins)==-1)// si l'element moins n'existe pas alors on decentre le schema vers plus
    {
      f_p1      =  vitesse(Tab_recap(element,plus),composante);
      f_m1      =  vitesse(Tab_recap(Tab_recap(element,plus),plus),composante);
      f         =  vitesse(element,composante);

    }
  else // si tous va bien
    {
      f_p1      =  vitesse(Tab_recap(element,plus),composante);
      f_m1      =  vitesse(Tab_recap(element,moins),composante);
      f 		  =  vitesse(element,composante);
      if( direction-1 ) centre++ ; // on indique le schemasi on est dans l'axe y on est forcement non isotrope
    }

  if(centre)
    {
      delta_m	  = distance_elem(element,plus) ; // distance a l'element proche

      derivee   = f_p1 - f_m1 ; //schema centre
      derivee  /= 2 * delta_m ;
    }
  else
    {
      delta_m   = distance_elem(element,moins); // la methode de calcul gere de decentrement et la periodicite
      delta_p   = distance_elem(element,plus);  // la methode de calcul gere de decentrement et la periodicite
      delta_m *= -1 ; // on prend en valeur absolue ;

      alpha 	  = delta_p/delta_m ; // rapport des distance
      unsalpha  = delta_m/delta_p ; //

      derivee   = (f_p1 -f ) * unsalpha  + (f - f_m1 )*alpha ;
      derivee  /= delta_m + delta_p;

    }
  return derivee;
}

double Traitement_particulier_NS_plan_VDF::derivee_deux_vitesse_coli(const DoubleTab& vitesse, int element,int direction) const
{
  // ici la periodicite est naturellement gerer car les calculs sont fait sur les faces interne a chaque element
  int    moins   = -1 ;
  int    plus    = -1 ;
  int 	centre  =  0 ;

  if (direction == 0)
    {
      moins = 0;
      plus = 1;
    }
  else if(direction == 1)
    {
      moins = 2;
      plus = 3;
      centre ++;
    }
  else if(direction == 2)
    {
      moins = 4;
      plus = 5;
    }
  else
    {
      Cerr << "Erreur dans le choix de la direction de derivation" << finl;
      Cerr << " x=0 y=1 z=2" << finl;
      exit();
    }
  const Domaine_dis_base&  zdisbase          =  mon_equation->inconnue().domaine_dis_base();
  const Domaine_VDF&       domaine_VDF          =  ref_cast(Domaine_VDF, zdisbase);
  const DoubleTab&      xv                =  domaine_VDF.xv();
  const IntTab&         elem_faces        =  domaine_VDF.elem_faces();

  double derivee;  // valeur retournee par la fonction
  double derivee_0; // valeur de la derivee premiere en 0
  double derivee_1; // valeur de la derivee premiere en 1
  double v_p;       // vitesse a la face suivantepoint
  double v_m;       // valeur a la face precedente
  double v_0;       // valeur a la face d'entree
  double v_1;       // vitesse a la face de sortie
  double delta_0m;  // distance entre les faces m et 0
  double delta_10;  // taille de la maille
  double delta_p1;  // distance entre les faces 1 et p
  double alpha; 	   // rapport des distances
  double unsalpha;
  // recupere le numero des faces de l'element et des deux faces voisines dans la direction de derivation
  int  face_0  = elem_faces(element,direction); // face d'entree
  int  face_1  = elem_faces(element,direction+dimension); // face de sortie
  int  face_m  = elem_faces(Tab_recap(element,moins),direction); // face precedente
  int  face_0b = elem_faces(Tab_recap(element,moins),direction+dimension);// face de sortie de la maille precedente
  int  face_p  = elem_faces(Tab_recap(element,plus),direction+dimension);  // face suivante
  int  face_1b  = elem_faces(Tab_recap(element,plus),direction);  // face d'entree de la maille suivante
  // on recupere les vitesses associees aux faces
  // | x | x | x |
  // m   0   1   p
  //
  v_m       =  vitesse[face_m];
  v_0       =  vitesse[face_0];
  v_1       =  vitesse[face_1];
  v_p       =  vitesse[face_p];
  // on calcule les distances
  delta_0m  = xv(face_0b,direction)-xv(face_m,direction);  //
  delta_10  = xv(face_1,direction)-xv(face_0,direction);  //
  delta_p1  = xv(face_p,direction)-xv(face_1b,direction);  //


  alpha 	  = delta_10/delta_0m ; // rapport des distance
  unsalpha  = delta_0m/delta_10 ; //

  // on calcule la derivee premiere en 0 tout en tenant compte des changement de taille de la maille
  derivee_0   = (v_1 - v_0 ) * unsalpha  + (v_0 - v_m ) * alpha ;
  derivee_0  /= delta_0m + delta_10;

  alpha 	  = delta_p1/delta_10 ; // rapport des distance
  unsalpha  = delta_10/delta_p1 ; //

  // on calcule la derivee premiere en 1 tout en tenant compte des changement de taille de la maille
  derivee_1   = (v_p - v_1 ) * unsalpha  + (v_1 - v_0 ) * alpha ;
  derivee_1  /= delta_p1 + delta_10;

  // calcule la derivee d'ordre deux
  derivee     = derivee_1-derivee_0 ;
  derivee    /= delta_10 ;

  return derivee;
}

double Traitement_particulier_NS_plan_VDF::derivee_deux_vitesse_nonco(DoubleTab& vitesse,int composante, int element, int direction) const
{
  // ici la periodicite est ignore car il n'y a pas de periodicite suivant y
  int    moins   = -1 ;
  int    plus    = -1 ;
  int    centre  =  0 ;
  if (direction == 0)
    {
      moins = 0;
      plus = 1;
    }
  else if(direction == 1)
    {
      moins = 2;
      plus = 3;
    }
  else if(direction == 2)
    {
      moins = 4;
      plus = 5;
    }
  else
    {
      Cerr << "Erreur dans le choix de la direction de derivation" << finl;
      Cerr << " x=0 y=1 z=2" << finl;
      exit();
    }

  double derivee;  // valeur retournee par la fonction
  double v_p;       // vitesse a la maille suivante
  double v_m;       // valeur a la maille precedente
  double v_0;       // valeur a la maille de calcul

  double delta_m;  // distance a la maille precedente
  double delta_p;  // distance a la maille suivante

  if(Tab_recap(element,plus)==-1) // si l'element plus n'existe pas alors on descentre le schema de cote de moins
    {
      v_p      =  vitesse(Tab_recap(Tab_recap(element,moins),moins),composante);
      v_m      =  vitesse(Tab_recap(element,moins),composante);
      v_0      =  vitesse(element,composante);

    }
  else if(Tab_recap(element,moins)==-1)// si l'element moins n'existe pas alors on decentre le schema vers plus
    {
      v_p      =  vitesse(Tab_recap(element,plus),composante);
      v_m      =  vitesse(Tab_recap(Tab_recap(element,plus),plus),composante);
      v_0      =  vitesse(element,composante);

    }
  else // si tous va bien
    {
      v_p      =  vitesse(Tab_recap(element,plus),composante) ;
      v_m      =  vitesse(Tab_recap(element,moins),composante) ;
      v_0	     =  vitesse(element,composante) ;
      if( direction-1 ) centre++ ; // on indique le schema si on est dans l'axe y on est forcement non isotrope
    }

  if(centre)
    {
      // on calcule les distances
      delta_m  = distance_elem(element,plus);  // norme de la distance avec la maille superieure

      // on calcule la derivee
      derivee   = v_p + v_m -2 * v_0 ;
      derivee  /= delta_m * delta_m ;
    }
  else
    {
      // on calcule les distances

      delta_p  = distance_elem(element,plus);  // distance avec la maille superieure
      delta_m  = distance_elem(element,moins); //  distance avec la maille inferieure
      if (delta_m < 0) delta_m *= -1;
      // on calcule la derivee tout en tenant compte des changement de taille de la maille
      derivee   = (v_p -v_0 ) / delta_p + ( v_0- v_m ) / delta_m ;
      derivee  /= delta_p +  delta_m ;
      derivee  *= 2 ;
    }
  return derivee;
}

double Traitement_particulier_NS_plan_VDF::distance_elem(int elem,int pos) const
{
  const Domaine_dis_base&  zdisbase  	=  mon_equation->inconnue().domaine_dis_base();
  const Domaine_VDF&       domaine_VDF 	 	=  ref_cast(Domaine_VDF, zdisbase);
  const DoubleTab&      xp       		=  domaine_VDF.xp();
  const DoubleTab&      xv        	=  domaine_VDF.xv();
  const IntTab&         elem_faces 	=  domaine_VDF.elem_faces();
  // on code ici une methode pour calculer la distance entre deux element voisin tenant compte de la periodicite ou des bords.
  // si on cherche la distance avec un element qui n'existe pas, on revois la distance dans la direction opposee avec le deuxieme element.
  // vois schema ci dessous
  //  <---->          //  voulue
  //    || e : x : x :
  //       < ----- >  // renvoyee

  // Rappel des directions dans Tab_recap (ne correspondent pas au direction Trio_U standard (a reguler un jour)
  //  0 : -ox | 1 : ox | 2 : -oy | 3 : oy | 4 : -oz | 5 : oz

  int direction = pos/2;  // la division euclidienne assure la bonne coordonee
  double distance=0; //Cerr << "   " << pos << " dir " << direction << "   " ;

  int voisin;
  if (Tab_recap(elem,pos) == -1) // si l'element n'existe pas
    {
      switch(pos) // on corrige la direction
        {
        case 0 :
          {
            pos = 1;
            break;
          }
        case 1 :
          {
            pos = 0;
            break;
          }
        case 2 :
          {
            pos = 3;
            break;
          }
        case 3 :
          {
            pos = 2;
            break;
          }
        case 4 :
          {
            pos = 5;
            break;
          }
        case 5 :
          {
            pos = 4;
            break;
          }
        default :
          {
            Cerr << pos << " Erreur dans la direction de recherche de la distance " << finl;
            exit();
          }
        }
      // il est impensable qu'une frontiere periodique soit dans le meme plan qu'une CL on calcule donc directement par les position
      voisin    = Tab_recap ( Tab_recap(elem,pos) ,pos );
      distance  = xp(voisin , direction ) - xp(elem,direction) ;
    }
  else
    {
      voisin	  =  Tab_recap(elem,pos);
      // on utilise la technique de calcul de distance pour des elements periodiques

      distance   =  xp(voisin,direction);
      distance  -=  xv(elem_faces(voisin,direction),direction);
      distance  +=  xv(elem_faces(elem,direction+3),direction) ;
      distance  -=  xp(elem,direction);
      if (!(pos%2)) // si on cherche un element dans le sens inverse a l'axe il faut inverser le signe.
        {
          distance *= -1 ;
        }
    }
  return distance;
}
