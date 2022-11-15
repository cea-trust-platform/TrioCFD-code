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
// File      : Trait_part_NS_surface.cpp
// Directory : $NEW_ALGO_QC_ROOT/src/modif_QC/Trait_part_NS_surface
//
/////////////////////////////////////////////////////////////////////////////

#include <Trait_part_NS_surface.h>
#include <LecFicDistribue.h>
#include <EcrFicCollecte.h>
#include <Navier_Stokes_Turbulent.h>
#include <Probleme_base.h>
#include <Discretisation_base.h>
#include <Zone_VF.h>
#include <Fluide_Quasi_Compressible.h>
#include <Zone_Cl_dis_base.h>
#include <Dirichlet_paroi_fixe.h>
#include <Schema_Temps_base.h>
#include <DoubleTrav.h>
#include <communications.h>
#include <Zone_VDF.h>
#include <Loi_Etat_GP.h>
#include <Schema_Temps_base.h>
#include <Schema_Temps.h>
#include <communications.h>

Implemente_base(Traitement_particulier_NS_surface,"Traitement_particulier_NS_surface",Traitement_particulier_NS_base);


/*! @brief 
 *
 * @param (Sortie& is) un flot de sortie 
 * @return (Sortie&) le flot de sortie modifie 
 */
Sortie& Traitement_particulier_NS_surface::printOn(Sortie& is) const
{
  return is;
}


/*! @brief 
 *
 * @param (Entree& is) un flot d'entree 
 * @return (Entree&) le flot d'entree modifie 
 */
Entree& Traitement_particulier_NS_surface::readOn(Entree& is)
{
  return is;
}

Entree& Traitement_particulier_NS_surface::lire(Entree& is)
{
  oui_profil_nu_t = 0;
  oui_profil_Temp = 0;
  oui_repr = 0;
  oui_pulse = 0;

  dt_impr_moy_spat=1e6;
  dt_impr_moy_temp=1e6;
  temps_deb=1e6;
  temps_fin=1e6;
  debut_phase=1e6;
  ind_phase=0;
  Nval=23;
  Nphase=1;
  Nb_ech_phase.resize(1);
  Nb_ech_phase(0)=1;
  Cerr << "Traitement particulier surface actif" << finl;

  if(!sub_type(Navier_Stokes_std,mon_equation.valeur()))
    {
      Cerr << " Traitement_particulier_surface has to be called from a Navier_Stokes problem " << finl;
      Process::exit();
    }

  if (mon_equation.valeur().probleme().nombre_d_equations()>1)
    try
      {
        Temp = mon_equation.valeur().probleme().equation(1).get_champ("temperature_qc"); // a changer en temperature pour ancien algo.
        oui_profil_Temp = 1 ;
        Nval+=5;
      }
    catch (Champs_compris_erreur)
      {
      }

  choix_fichier=0;
  Motcle accouverte = "{" , accfermee = "}" ;
  Motcle motbidon, motlu;
  is >> motbidon ;
  if (motbidon == accouverte)
    {
      Motcles les_mots(7);
      les_mots[0] = "dt_impr";
      les_mots[1] = "debut";
      les_mots[2] = "gnuplot";
      les_mots[3] = "planxy";
      is >> motlu;
      while(motlu != accfermee)
        {
          int rang=les_mots.search(motlu);
          switch(rang)
            {
            case 0 :
              {
                is >> dt_impr_moy_spat;      // intervalle de temps de sorties des moyennes spatiales
                Cerr << "Post-traitement surface imprime tout les : dt_impr_moy_spat = " << dt_impr_moy_spat << finl;
                break;
              }
            case 1 :
              {
                is >> temps_deb;	     // temps de debut des moyennes temporelles
                Cerr << "Post-traitement surface  commence a  : temps_deb = " << temps_deb << finl;
                break;
              }
            case 2 :
              {
                Cerr << "Affichage au format Gnuplot "<< finl;
                choix_fichier=1;
                break;
              }
            case 3 :
              {
                Cerr << "Affichage au format Planxy "<< finl;
                choix_fichier=2;
                break;
              }
            default :
              {
                Cerr << motlu << " is not a keyword for Traitement_particulier_NS_surface use: *dt_impr*, *debut*, *gnuplot* or *planxy*" << finl;
                Process::exit();
              }
            }
          is >> motlu;
        }
      is >> motlu;
      if (motlu != accfermee)
        {
          Cerr << "Error while reading in Traitement_particulier_NS_surface" << finl;
          Cerr << "We were expecting a }" << finl;
          Process::exit();
        }
    }
  else
    {
      Cerr << "Error while reading in Traitement_particulier_NS_surface" << finl;
      Cerr << "We were expecting a {" << finl;
      Process::exit();
    }

  return is;
}

void Traitement_particulier_NS_surface::preparer_calcul_particulier()
{
  const Zone_dis_base& 				zdisbase	= mon_equation->inconnue().zone_dis_base();
  const Zone_VF& 					zone_VF		= ref_cast(Zone_VF, zdisbase);
  const DoubleTab& 					xp 		= zone_VF.xp();
  int 								nb_elems 	= zone_VF.zone().nb_elem();

  int num_elem	= 0; // compteur des elements
  int i,j,l;	  // compteurs des boucles

  remplir_XYZ(X,Y,Z,Nx,Ny,Nz,Tab_recap); // renvoie vers Traitement_particulier_NS_surface_VDF ou Traitement_particulier_NS_surface_VEF

  remplir_reordonne_Z_tot(Z,Z_tot);
  remplir_reordonne_Y_tot(Y,Y_tot);
  remplir_reordonne_X_tot(X,X_tot);
  envoyer_broadcast(X_tot,0);
  envoyer_broadcast(Y_tot,0);
  envoyer_broadcast(Z_tot,0);
  int Nxt = X_tot.size();
  int Nyt = Y_tot.size();

  int maxX = ::mp_max( X.size() ) ;
  int maxY = ::mp_max( Y.size() ) ;

  Tab_post.resize(nb_elems,2);
  IntTab posX(maxX);
  IntTab posY(maxY);
  posX = -1 ;
  posY = -1 ;
  for ( num_elem=0; num_elem<nb_elems; num_elem++) // boucle sur les elements
    {
      for(i=0; i<X.size(); i++) // boucle sur les X
        if(xp(num_elem,0)==X(i))
          {
            break; // on cherche l'element correspondant dans X_tot
          }

      for(j=0; j<Y.size(); j++) // boucles sur les surfaces a post traite
        if(xp(num_elem,1)==Y(j))
          {
            break;
          }

      Tab_post(num_elem,0)=i;
      Tab_post(num_elem,1)=j;

      for(l=0; l<Nxt; l++) // boucle sur les X
        if( X(i) == X_tot(l))
          {
            break; // on cherche l'element correspondant dans X_tot
          }
      posX(i)= l;

      for(l=0; l<Nyt; l++) //Boucle sur les y
        if( Y(j) == Y_tot(l))
          {
            break; // Idem pour Y_tot
          }
      posY(j) = l;
    }

  envoyer(posX ,Process::me(),0,Process::me()); // on dis si on va post traite.
  envoyer(posY ,Process::me(),0,Process::me()); // on dis si on va post traite.
  if(je_suis_maitre())
    {
      PosX.resize(maxX,Process::nproc());
      PosY.resize(maxY,Process::nproc());
      for(int p=0; p<Process::nproc(); p++)
        {
          recevoir(posX,p,0,p);
          recevoir(posY,p,0,p);
          for (i = 0 ; i < posX.size() ; i++ )	PosX(i,p) = (int)posX(i) ;
          for (i = 0 ; i < posY.size() ; i++ )	PosY(i,p) = (int)posY(i) ;
        }
    }
}

void Traitement_particulier_NS_surface::remplir_reordonne_Z_tot(const DoubleVect& Zc, DoubleVect& Z_total) const
{
  DoubleVect Z_p(Zc);
  int j =0;
  envoyer(Z_p, Process::me(), 0, Process::me());

  if(je_suis_maitre())
    {
      Z_total.resize(1);
      Z_total(0) = Z(0);

      for(int p=0; p<Process::nproc(); p++)
        {
          recevoir(Z_p,p,0,p);

          for (int k=0; k<Z_p.size(); k++)
            {
              for ( j=0; j<Z_total.size(); j++)
                {
                  if (est_egal(Z_total(j),Z_p(k))) break;
                }
              if(j==Z_total.size())
                {
                  Z_total.resize(Z_total.size()+1);
                  Z_total(Z_total.size()-1) = Z_p(k);
                }
            }
        }


      if(est_egal(Z_total(0),1.e12))
        {
          Cerr << "Probleme a la construction de Z_tot" << finl;
          Process::exit();
        }

      Z_total.ordonne_array();
    }
}
void Traitement_particulier_NS_surface::remplir_reordonne_Y_tot(const DoubleVect& Yc, DoubleVect& Y_total) const
{
  DoubleVect Y_p(Yc);
  int j =0;
  envoyer(Y_p, Process::me(), 0, Process::me());

  if(je_suis_maitre())
    {
      Y_total.resize(1);
      Y_total(0) = Y(0);

      for(int p=0; p<Process::nproc(); p++)
        {
          recevoir(Y_p,p,0,p);
          for (int k=0; k<Y_p.size(); k++)
            {
              for ( j=0; j<Y_total.size(); j++)
                {
                  if (est_egal(Y_total(j),Y_p(k))) break;
                }
              if(j==Y_total.size())
                {
                  Y_total.resize(Y_total.size()+1);
                  Y_total(Y_total.size()-1) = Y_p(k);
                }
            }
        }

      if(est_egal(Y_total(0),1.e12))
        {
          Cerr << "Probleme a la construction de Y_tot" << finl;
          Process::exit();
        }
      Y_total.ordonne_array();
    }
  envoyer_broadcast(Y_total,0);
}

void Traitement_particulier_NS_surface::remplir_reordonne_X_tot(const DoubleVect& Xc, DoubleVect& X_total) const
{
  DoubleVect X_p(Xc);
  int j=0;

  envoyer(X_p, Process::me(), 0, Process::me());

  if(je_suis_maitre())
    {
      X_total.resize(1);
      X_total(0) = X(0);

      for(int p=0; p<Process::nproc(); p++)
        {
          recevoir(X_p,p,0,p);

          for (int k=0; k<X_p.size(); k++)
            {
              for ( j=0; j<X_total.size(); j++)
                {
                  if (est_egal(X_total(j),X_p(k)))
                    {
                      break;
                    }
                }
              if(j==X_total.size())
                {
                  X_total.resize(X_total.size()+1);
                  X_total(X_total.size()-1) = X_p(k);
                }
            }
        }

      if(est_egal(X_total(0),1.e12))
        {
          Cerr << "Probleme a la construction de X_tot" << finl;
          Process::exit();
        }

      X_total.ordonne_array();
    }
}
void Traitement_particulier_NS_surface::reprendre_stat()
{
  /// Apartient a trait_NS_base donc ne peu etre supprimer

}

void Traitement_particulier_NS_surface::sauver_stat() const
{
  /// Apartient a trait_NS_base donc ne peu etre supprimer
}

void Traitement_particulier_NS_surface::post_traitement_particulier()
{
  const Fluide_Incompressible&      le_fluide 	= ref_cast(Fluide_Incompressible,mon_equation->milieu());
  const Fluide_Quasi_Compressible&  fluide_QC    = ref_cast(Fluide_Quasi_Compressible,le_fluide);
  const Loi_Etat_GP& 		   		loi_		= ref_cast(Loi_Etat_GP,fluide_QC.loi_etat().valeur()); // recupere la loi des GP
  double 			   				R			= loi_.R();
  R=R;
  double tps = mon_equation->inconnue().temps();
  if (tps>temps_deb)
    {
      double stationnaire_atteint = mon_equation->schema_temps().stationnaire_atteint();
      double tps_passe		   = mon_equation->schema_temps().temps_courant();
      double nb_pas_dt_max	   = mon_equation->schema_temps().nb_pas_dt_max();
      double nb_pas_dt		   = mon_equation->schema_temps().nb_pas_dt();
      double temps_max		   = mon_equation->schema_temps().temps_max();
      int impr_inst;

      if ( dt_impr_moy_spat<=(tps-tps_passe) )
        impr_inst=1;
      else
        {
          // Voir Schema_Temps_base::limpr pour information sur epsilon et modf
          double i, j, epsilon = 1.e-8;
          modf(tps/dt_impr_moy_spat + epsilon, &i);
          modf(tps_passe/dt_impr_moy_spat + epsilon, &j);
          impr_inst=(i>j);
        }

      if ( (nb_pas_dt+1<=1) || impr_inst || (temps_max <= tps) || (nb_pas_dt_max <= nb_pas_dt+1) || stationnaire_atteint)
        {
          // Calcul des Moyennes spatiales
          //////////////////////////////////////////////////////////

          const Zone_dis_base& zdisbase=mon_equation->inconnue().zone_dis_base();
          const Zone_VDF& zone_VDF=ref_cast(Zone_VDF, zdisbase);

          int nb_elem_p = zone_VDF.zone().nb_elem();
          int num_elem;

          int i=0,j=0,k =0;
          int val =0 ;
          DoubleTrav val_post(nb_elem_p,Nval);
          recuperation_grandeurs(val_post);
          // ici on va moyenner localement.
          int Nxp = X.size() ;
          int Nyp = Y.size() ;

          DoubleTrav val_plan(Nxp,Nyp,Nval);

          for (num_elem = 0 ; num_elem < nb_elem_p ; num_elem ++ )
            for ( val = 0 ; val < Nval ; val ++ )
              val_plan( Tab_post(num_elem,0) , Tab_post(num_elem,1) , val ) +=  val_post(num_elem,val);

          val_plan /= Z_tot.size() ;

          // Echange des donnees entre processeurs

          envoyer(Nxp,Process::me(),0,Process::me());
          envoyer(Nyp,Process::me(),0,Process::me());
          envoyer(val_plan,Process::me(),0,Process::me());
          int Nxt = X_tot.size();
          int Nyt = Y_tot.size();
          DoubleTrav  Moy_final;

          if(je_suis_maitre())
            {

              Moy_final.resize(Nxt, Nyt, Nval); // Doubletrav initialise le tableau
              Moy_final = 0 ;
              for(int p=0; p<Process::nproc(); p++)
                {
                  recevoir(Nxp,p,0,p);
                  recevoir(Nyp,p,0,p);
                  recevoir(val_plan,p,0,p);
                  for(i=0; i<Nxp; i++)
                    for (j=0; j<Nyp; j++)
                      for ( val = 0 ; val < Nval ; val ++ )
                        if ( ( PosX(i,p) + 1 ) && ( PosY(j,p)+ 1 ) )
                          Moy_final(PosX(i,p),PosY(j,p),val) += val_plan(i,j,val);
                }

              double produit = 0 ;

              for ( k=10; k<28; k++) // on commence a 10 on ne calcul ici que des fluctutations.
                for ( j=0; j<Nyt; j++)
                  for ( i=0; i<Nxt; i++)
                    {
                      if( k < 20 ) // calcul des valeurs RMS en valeurs absolues
                        {
                          produit = Moy_final(i, j, k) - Moy_final(i, j, k-10)*Moy_final(i, j, k-10) ;
                          if ( 0 > produit )
                            Moy_final(i, j, k) = -produit;
                          else
                            Moy_final(i, j, k) =  produit;
                        }
                      else
                        switch(k) // calcul des correlations
                          {
                          case 20 : // <u'.v'> = <u.v> -<u>.<v>
                            {
                              Moy_final(i, j, k) = (Moy_final(i, j, k) - Moy_final(i, j, 0)*Moy_final(i, j, 1));
                              break;
                            }
                          case 21 : // <u'.w'> = <u.w> -<u>.<w>
                            {
                              Moy_final(i, j, k) = (Moy_final(i, j, k) - Moy_final(i, j, 0)*Moy_final(i, j, 2));
                              break;
                            }
                          case 22 : // <v'.w'> = <v.w> -<v>.<w>
                            {
                              Moy_final(i, j, k) = (Moy_final(i, j, k) - Moy_final(i, j, 1)*Moy_final(i, j, 2));
                              break;
                            }
                          case 23 : // <u'.rho'> = <u.rho> -<u>.<rho>
                            {
                              Moy_final(i, j, k) = (Moy_final(i, j, k) - Moy_final(i, j, 0)*Moy_final(i, j, 5));
                              break;
                            }
                          case 24 : // <v'.rho'> = <v.rho> -<v>.<rho>
                            {
                              Moy_final(i, j, k) = (Moy_final(i, j, k) - Moy_final(i, j, 1)*Moy_final(i, j, 5));
                              break;
                            }
                          case 25 : // <u'.T'> = <u.T> -<u>.<T>
                            {
                              Moy_final(i, j, k) = (Moy_final(i, j, k) - Moy_final(i, j, 0)*Moy_final(i, j, 3));
                              break;
                            }
                          case 26 : // <v'.T'> = <v.T> -<v>.<T>
                            {
                              Moy_final(i, j, k) = (Moy_final(i, j, k) - Moy_final(i, j, 1)*Moy_final(i, j, 3));
                              break;
                            }

                          case 27 : // <rho'.T'> = <rho.T> -<rho>.<T>
                            {
                              Moy_final(i, j, k) = (Moy_final(i, j, k) - Moy_final(i, j, 5)*Moy_final(i, j, 3));
                              break;
                            }
                          default :
                            {
                              break;
                            }

                          }

                    }

              // Ecriture des Moyennes spatiales et temporelles
              ///////////////////////////////////////////////////////////////
              i = 0;
              Nom fichier1;
              Nom tampon ;
              Nom temps = Nom(tps);

              for( k=0; k<Nval; k++)
                {
                  switch(k)
                    {
                    case 0 :
                      {
                        tampon  = "U_";
                        tampon += temps;
                        break;
                      }
                    case 1 :
                      {
                        tampon  = "V_";
                        tampon += temps;
                        break;
                      }
                    case 2 :
                      {
                        tampon  = "W_";
                        tampon += temps;
                        break;
                      }
                    case 3 :
                      {
                        tampon  = "TEMP_";
                        tampon += temps;
                        break;
                      }
                    case 4 :
                      {
                        tampon  = "P_dyn_";
                        tampon += temps;
                        break;
                      }
                    case 5 :
                      {
                        tampon  = "RHO_";
                        tampon += temps;
                        break;
                      }
                    case 6 :
                      {
                        tampon  = "LAMBDA_";
                        tampon += temps;
                        break;
                      }
                    case 7 :
                      {
                        tampon  = "MU_";
                        tampon += temps;
                        break;
                      }
                    case 8 :
                      {
                        tampon  = "MU_turb_";
                        tampon += temps;
                        break;
                      }
                    case 9 :
                      {
                        tampon  = "NU_";
                        tampon += temps;
                        break;
                      }
                    case 10 :
                      {
                        tampon  = "u_";
                        tampon += temps;
                        break;
                      }
                    case 11 :
                      {
                        tampon  = "v_";
                        tampon += temps;
                        break;
                      }
                    case 12 :
                      {
                        tampon  = "w_";
                        tampon += temps;
                        break;
                      }
                    case 13 :
                      {
                        tampon  = "temp_";
                        tampon += temps;
                        break;
                      }
                    case 14 :
                      {
                        tampon  = "p_dyn_";
                        tampon += temps;
                        break;
                      }
                    case 15 :
                      {
                        tampon  = "rho_";
                        tampon += temps;
                        break;
                      }
                    case 16 :
                      {
                        tampon  = "lambda_";
                        tampon += temps;
                        break;
                      }
                    case 17 :
                      {
                        tampon  = "mu_";
                        tampon += temps;
                        break;
                      }
                    case 18 :
                      {
                        tampon  = "mu_turb_";
                        tampon += temps;
                        break;
                      }
                    case 19 :
                      {
                        tampon  = "nu_";
                        tampon += temps;
                        break;
                      }
                    case 20 :
                      {
                        tampon  = "u.v_";
                        tampon += temps;
                        break;
                      }
                    case 21 :
                      {
                        tampon  = "u.w_";
                        tampon += temps;
                        break;
                      }
                    case 22 :
                      {
                        tampon  = "v.w_";
                        tampon += temps;
                        break;
                      }
                    case 23 :
                      {
                        tampon  = "u.rho_";
                        tampon += temps;
                        break;
                      }
                    case 24 :
                      {
                        tampon  = "v.rho_";
                        tampon += temps;
                        break;
                      }
                    case 25 :
                      {
                        tampon  = "u.temp_";  //
                        tampon += temps;
                        break;
                      }
                    case 26 :
                      {
                        tampon  = "v.temp_";  //
                        tampon += temps;
                        break;
                      }
                    case 27 :
                      {
                        tampon  = "rho.temp_";  //
                        tampon += temps;
                        break;
                      }
                    default :
                      {
                        Cerr << " pas de valeur definie " << finl ;
                        Process::exit(0);
                        break;
                      }
                    }

                  fichier1 = tampon;
                  fichier1 += "_s_.donnees";

                  if (choix_fichier==2)
                    {
                      ecriture_fichiers_moy_vitesse_rho_mu(Moy_final,fichier1,i,k);
                    }
                  else if (choix_fichier==1)
                    {
                      fichier1 += "gnu";
                      ecriture_fichiers_moy_vitesse_rho_mu_gnuplot(Moy_final,fichier1,i,k);
                    }
                  else
                    {
                      ecriture_fichiers_moy_vitesse_rho_mu(Moy_final,fichier1,i,k);
                      fichier1 += "gnu";
                      ecriture_fichiers_moy_vitesse_rho_mu_gnuplot(Moy_final,fichier1,i,k);
                    }
                } // for(int k=0;k<Nval;k++)
            } // if(je_suis_maitre())
        } //  if ( (nb_pas_dt+1<=1) || impr_inst || (temps_max <= tps) || (nb_pas_dt_max <= nb_pas_dt+1) || stationnaire_atteint)
    } // if (tps>temps_deb)

}

void Traitement_particulier_NS_surface::calcul_reynolds_tau()
{
  /*

   double tps = mon_equation->inconnue().temps();

   // Aspects geometriques du maillage
   /////////////////////////////////////////////////////////

  // !!!!!  Hypotheses : maillage symetrique suivant la demi-hauteur et s'etendant de Z=0 a Z=H

   Equation_base& eqn=mon_equation.valeur();
   const Discretisation_base& discr=eqn.discretisation();
   Nom nom_discr=discr.que_suis_je();
   int k=0;
   if(nom_discr=="VEFPreP1B" || nom_discr=="VEF") k=1;

   int Nzt = Z_tot.size();
   int Nxt = X_tot.size();
   int kmin=k;		    // indice du premier point hors paroi
   int kmax=Nzt-1-k;		    // indice du dernier point hors paroi
   double ymin=Z_tot(kmin);   	    // position du premier point
   double ymax=Z_tot(kmax);   	    // position du dernier point
   double hs2=0.5*(ymin+ymax);       // demi-hauteur

   // definition des grandeurs interveannt dans la suite des calculs
   //////////////////////////////////////////////////////

       SFichier fic1("reynolds_tau.dat",ios::app);
       SFichier fic2("u_tau.dat",ios::app);
       SFichier fic3("tauw.dat",ios::app);

       fic1 << "# (0)  : X " << finl;
       fic1 << "# (1)  : t " << finl;
       fic1 << "# (2)  : Reynolds tau moyen " << finl;
       fic1 << "# (3)  : Reynolds tau bas " << finl;
       fic1 << "# (4)  : Reynolds tau haut " << finl;

       fic2 << "# (0)  : X " << finl;
       fic2 << "# (1)  : t " << finl;
       fic2 << "# (2)  : u tau moyenne  " << finl;
       fic2 << "# (3)  : u tau basse  " << finl;
       fic2 << "# (4)  : u tau haute  " << finl;

       fic3 << "# (0)  : X " << finl;
       fic1 << "# (1)  : t " << finl;
       fic3 << "# (2)  : Constante de cisallement moyenne  " << finl;
       fic3 << "# (3)  : Constante de cisallement basse " << finl;
       fic3 << "# (4)  : Constante de cisallement haute " << finl;


   DoubleVect mu_bas, mu_haut;   	 // viscosite dynamique
   DoubleVect rho_bas, rho_haut;	 // masse volumique
   DoubleVect utang_bas, utang_haut;	 // norme de la vitesse tangente a la paroi
   DoubleVect tauwb, tauwh, tauwm;	 // cisaillement a la paroi
   DoubleVect utaub, utauh, utaum;	 // vitesse de frottement
   DoubleVect retaub, retauh, retaum; // Reynolds de frottement

   mu_bas.resize(Nxt) ;
   mu_haut.resize(Nxt);

   rho_bas.resize(Nxt)  ;
   rho_haut.resize(Nxt) ;

   utang_bas.resize(Nxt)  ;
   utang_haut.resize(Nxt) ;

   utaub.resize(Nxt)  ;
   utauh.resize(Nxt) ;
   utaum.resize(Nxt) ;

   tauwb.resize(Nxt)  ;
   tauwh.resize(Nxt) ;
   tauwm.resize(Nxt) ;

   retaub.resize(Nxt)  ;
   retauh.resize(Nxt) ;
   retaum.resize(Nxt) ;

  for (int l=0;l<Nxt;l++)
  {
   mu_bas[l]  = val_moy_tot(l,kmin,11,0);
   mu_haut[l] = val_moy_tot(l,kmax,11,0);

   rho_bas[l]  = val_moy_tot(l,kmin,10,0);
   rho_haut[l] = val_moy_tot(l,kmax,10,0);

   utang_bas[l]  = val_moy_tot(l,kmin,9,0);
   utang_haut[l] = val_moy_tot(l,kmax,9,0);


   // calcul du cisaillement a la paroi suivant la vitesse tangentielle
   //////////////////////////////////////////////////////

   // approximation lineaire : Tau_w = mu * ||u_t||(y1) / dy1

    tauwb[l]= mu_bas[l] * utang_bas[l]  / ymin;
    tauwh[l]= mu_haut[l]* utang_haut[l] / ymin;


   // calcul et ecritures des differentes grandeurs parietales
   //////////////////////////////////////////////////////

       utaub[l]=sqrt(tauwb[l]/rho_bas[l]);  if(val_moy_tot(l,kmin,0,0)<=0) utaub[l]*=-1.;
       utauh[l]=sqrt(tauwh[l]/rho_haut[l]); if(val_moy_tot(l,kmax,0,0)<=0) utauh[l]*=-1.;
       retaub[l]=rho_bas[l]*utaub[l]*hs2/mu_bas[l];
       retauh[l]=rho_haut[l]*utauh[l]*hs2/mu_haut[l];
       tauwm[l]=(tauwh[l]+tauwb[l])/2.;
       retaum[l]=(retauh[l]+retaub[l])/2.;
       utaum[l]=(utauh[l]+utaub[l])/2.;

       fic1 << X_tot(l) << "   " << tps << " " << retaum[l] << " " << retaub[l] << " " << retauh[l] << finl;
       fic2 << X_tot(l) << "   " << tps << " " << utaum[l] << " "  << utaub[l] << " "  << utauh[l]  << finl;
       fic3 << X_tot(l) << "   " << tps << " " << tauwm[l] << " "  << tauwb[l] << " "  << tauwh[l]  << finl;
       }// fin de la boucle sur les x
       fic1<<flush;
       fic1.close();
       fic2<<flush;
       fic2.close();
       fic3<<flush;
       fic3.close();
  */
}


void Traitement_particulier_NS_surface::ecriture_fichiers_moy_vitesse_rho_mu(const DoubleTab& Moy_final, const Nom& fichier, int& surface, int& k) const
{
  //const scientifique = ; // Affichage au format scientifique (virgule flottante)
  SFichier fic (fichier);
  fic.setf(ios::scientific);
  double tps = mon_equation->inconnue().temps();
  int i,j;

  fic << "# Grandeur calculee:  ";

  switch(k)
    {

    case 0 :
      {
        fic << "<U>";
        break;
      }
    case 1 :
      {
        fic << "<V>";
        break;
      }
    case 2 :
      {
        fic << "<W>";
        break;
      }
    case 3 :
      {
        fic << "<T>";
        break;
      }
    case 4 :
      {
        fic << "<P_dyn>";
        break;
      }
    case 5 :
      {
        fic << "<RHO>";
        break;
      }
    case 6 :
      {
        fic << "<LAMBDA>";
        break;
      }
    case 7 :
      {
        fic << "<MU>";
        break;
      }
    case 8 :
      {
        fic << "<MU_turb>";
        break;
      }
    case 9 :
      {
        fic <<  "<NU>";
        break;
      }
    case 10 :
      {
        fic <<  "abs(<U U>-<U>*<U>)";
        break;
      }
    case 11 :
      {
        fic <<  "abs(<V V>-<V>*<V>)";
        break;
      }
    case 12 :
      {
        fic <<  "abs(<W W>-<W>*<W>)";
        break;
      }
    case 13 :
      {
        fic <<  "abs(<T T>-<T>*<T>)";
        break;
      }
    case 14 :
      {
        fic <<  "abs(<P.P>-<P>*<P>)";
        break;
      }
    case 15 :
      {
        fic <<  "abs(<RHO RHO>-<RHO><RHO>)";
        break;
      }
    case 16 :
      {
        fic <<  "abs(<LAMBDA LAMBDA>-<LAMBDA><LAMBDA>)";
        break;
      }
    case 17 :
      {
        fic <<  "abs(<MU MU>-<MU><MU>)";
        break;
      }
    case 18 :
      {
        fic <<  "abs(<MU_turb MU_turb>-<MU_turb>< MU_turb>";
        break;
      }
    case 19 :
      {
        fic <<  "abs(<NU NU>-<NU><NU>)";
        break;
      }
    case 20 :
      {
        fic <<  "<U V> - <U><V>";    //
        break;
      }
    case 21 :
      {
        fic <<  "<U W> - <U><W>";    //
        break;
      }
    case 22 :
      {
        fic <<  "<V W> - <V><W>";    //
        break;
      }
    case 23 :
      {
        fic <<  "<U RHO> - <U><RHO>";    //
        break;
      }
    case 24 :
      {
        fic <<  "<V RHO> - <V><RHO>";    //
        break;
      }
    case 25 :
      {
        fic <<  "<U T> - <U><T>";
        break;
      }
    case 26 :
      {
        fic <<  "<V T> - <V><T>";
        break;
      }
    case 27 :
      {
        fic <<  "<T RHO> - <T><RHO>";
        break;
      }
    default :
      {
        Cerr <<  "Erreur sur k "<< finl;  // ne peu jamais passer ici.
        break;
      }
    }

  fic << finl;
  fic << "# Temps  "  << "  " << " fichier : " <<finl;
  fic << "# "<< tps << "   "  <<  fichier << finl;
  fic << "#  X =  x1  x2  " << finl;
  fic << "#  Y = y1   " << finl;
  fic << "#  Y = y2   " << finl;
  fic << " #    /          ";
  for (i=0; i<X_tot.size(); i++)
    {
      fic << X_tot[i] << "    ";
    }
  fic << finl;
  for (j=0; j<Y_tot.size(); j++)
    {

      fic << Y_tot(j) << "    " ;
      for (i=0; i<X_tot.size(); i++)
        {

          fic << Moy_final(i, j, k) << "    " ;


        }
      fic << finl;
    }
  fic.flush();
  fic.close();
}

void Traitement_particulier_NS_surface::ecriture_fichiers_moy_vitesse_rho_mu_gnuplot(const DoubleTab& Moy_final, const Nom& fichier, int& surface, int& k) const
{
  //const scientifique = ; // Affichage au format scientifique (virgule flottante)
  SFichier fic (fichier);
  fic.setf(ios::scientific);
  double tps = mon_equation->inconnue().temps();
  int i,j;

  fic << "# Grandeur calculee:  ";

  switch(k)
    {

    case 0 :
      {
        fic << "<U>";
        break;
      }
    case 1 :
      {
        fic << "<V>";
        break;
      }
    case 2 :
      {
        fic << "<W>";
        break;
      }
    case 3 :
      {
        fic << "<T>";
        break;
      }
    case 4 :
      {
        fic << "<P_dyn>";
        break;
      }
    case 5 :
      {
        fic << "<RHO>";
        break;
      }
    case 6 :
      {
        fic << "<LAMBDA>";
        break;
      }
    case 7 :
      {
        fic << "<MU>";
        break;
      }
    case 8 :
      {
        fic << "<MU_turb>";
        break;
      }
    case 9 :
      {
        fic <<  "<NU>";
        break;
      }
    case 10 :
      {
        fic <<  "abs(<U U>-<U>*<U>)";
        break;
      }
    case 11 :
      {
        fic <<  "abs(<V V>-<V>*<V>)";
        break;
      }
    case 12 :
      {
        fic <<  "abs(<W W>-<W>*<W>)";
        break;
      }
    case 13 :
      {
        fic <<  "abs(<T T>-<T>*<T>)";
        break;
      }
    case 14 :
      {
        fic <<  "abs(<P.P>-<P>*<P>)";
        break;
      }
    case 15 :
      {
        fic <<  "abs(<RHO RHO>-<RHO><RHO>)";
        break;
      }
    case 16 :
      {
        fic <<  "abs(<LAMBDA LAMBDA>-<LAMBDA><LAMBDA>)";
        break;
      }
    case 17 :
      {
        fic <<  "abs(<MU MU>-<MU><MU>)";
        break;
      }
    case 18 :
      {
        fic <<  "abs(<MU_turb MU_turb>-<MU_turb>< MU_turb>";
        break;
      }
    case 19 :
      {
        fic <<  "abs(<NU NU>-<NU><NU>)";
        break;
      }
    case 20 :
      {
        fic <<  "<U V> - <U><V>";    //
        break;
      }
    case 21 :
      {
        fic <<  "<U W> - <U><W>";    //
        break;
      }
    case 22 :
      {
        fic <<  "<V W> - <V><W>";    //
        break;
      }
    case 23 :
      {
        fic <<  "<U RHO> - <U><RHO>";    //
        break;
      }
    case 24 :
      {
        fic <<  "<V RHO> - <V><RHO>";    //
        break;
      }
    case 25 :
      {
        fic <<  "<U T> - <U><T>";
        break;
      }
    case 26 :
      {
        fic <<  "<V T> - <V><T>";
        break;
      }
    case 27 :
      {
        fic <<  "<T RHO> - <T><RHO>";
        break;
      }
    default :
      {
        Cerr <<  "Erreur sur k "<< finl;  // ne peu jamais passer ici.
        break;
      }
    }

  fic << finl;
  fic << "# Temps  "  << "  " << " fichier : " <<finl;
  fic << "# "<< tps << "   "  <<  fichier << finl;
  fic <<"#   X             Y             Grandeur " <<finl;

  for (i=0; i<X_tot.size(); i++)
    {
      for (j=0; j<Y_tot.size(); j++)
        {
          fic << X_tot[i]<<"    "<<Y_tot[j] <<"    "<<Moy_final(i, j, k)<<finl;
        }
      fic << "               " << finl;
    }
  fic.flush();
  fic.close();
}

void Traitement_particulier_NS_surface::ecriture_fichiers_moy_nut(const DoubleTab& val_moy, const Nom& fichier, const double& dt, const int& k) const
{
  /*
      SFichier fic (fichier);

      double tps = mon_equation->inconnue().temps();
      int i;

      fic << "# Temps : " << tps << finl;
      fic << " " << finl;

      fic << "# (0)  : X " << finl;
      fic << "# (1)  : Z " << finl;
      fic << "# (2)  : <nu_t> " << finl;

      fic << "#  " << finl;
      fic << "# Nombre minimum de points par segment (X,Z)=cste ayant servi au calcul de la moyenne spatiale : " << min_array(compt_tot) << finl;
      fic << "#  " << finl;
     for (int l=0;l<X_tot.size();l++)
       for (i=0;i<Z_tot.size();i++)
  fic << X_tot(l) << "   " << Z_tot(i) << "    " << val_moy(l,i,12,k)/dt << finl;

      fic.flush();
      fic.close();

      */

}

void Traitement_particulier_NS_surface::ecriture_fichiers_moy_Temp(const DoubleTab& val_moy, const Nom& fichier, const double& dt, const int& k) const
{
  /*
      SFichier fic (fichier);

      double tps = mon_equation->inconnue().temps();
      double u,v,w,T,T2,uT,vT,wT;
      int i;

      fic << "# Temps : " << tps << finl;
      fic << " " << finl;

      fic << "# (0)  : Z " << finl;
      fic << "# (1)  : <T> " << finl;
      fic << "# (2)  : sqrt(<T.T>-<T>*<T>) " << finl;
      fic << "# (3)  : <u.T>-<u>*<T> " << finl;
      fic << "# (4)  : <v.T>-<v>*<T> " << finl;
      fic << "# (5)  : <w.T>-<w>*<T> " << finl;

      fic << "#  " << finl;
      fic << "# Nombre minimum de points par segment (X,Z)=cste ayant servi au calcul de la moyenne spatiale : " << min_array(compt_tot) << finl;
      fic << "#  " << finl;
   for (int l=0;l<X_tot.size();l++)
   {
   fic << " X = " << X_tot(l) << finl;
      for (i=0;i<Z_tot.size();i++)
  {
    u    = val_moy(l,i,0,k)/dt;
    v    = val_moy(l,i,1,k)/dt;
    w    = val_moy(l,i,2,k)/dt;
    T    = val_moy(l,i,13,k)/dt;
    T2   = val_moy(l,i,14,k)/dt;
    uT   = val_moy(l,i,15,k)/dt;
    vT   = val_moy(l,i,16,k)/dt;
    wT   = val_moy(l,i,17,k)/dt;


    fic << Z_tot(i) << "    "  ;
    // Pour eviter NAN, on prend le dmax(,0):
    fic << T << "    "  <<sqrt(dmax(T2,0))<<" ";
    fic << -uT << "    " << -vT  << "    " << -wT  << "    ";
    fic << finl;
  }
    }
      fic.flush();
      fic.close();

    */

}

//Apres correction de l expression de la moyenne temporelle et des methodes d ecriture des moyennes temporelles
//On conserve ces deux methodes d ecriture en version old pour l ecriture des moyennes de phase

void Traitement_particulier_NS_surface::ecriture_fichiers_moy_vitesse_rho_mu_old(const DoubleTab& val_moy, const Nom& fichier, const double& dt, const int& k) const
{
  /*  SFichier fic (fichier);

    double tps = mon_equation->inconnue().temps();
    double u,v,w,u2,v2,w2,uv,uw,vw,rho,mu;
    int i;

    fic << "# Temps : " << tps << finl;
    fic << " " << finl;

    fic << "# (0)  : Z " << finl;
    fic << "# (1)  : <u> " << finl;
    fic << "# (2)  : <v> " << finl;
    fic << "# (3)  : <w> " << finl;
    fic << "# (4)  : sqrt(<u.u>-<u>*<u>) " << finl;
    fic << "# (5)  : sqrt(<v.v>-<v>*<v>) " << finl;
    fic << "# (6)  : sqrt(<w.w>-<w>*<w>) " << finl;
    fic << "# (7)  : <u.v>-<u>*<v> " << finl;
    fic << "# (8)  : <u.w>-<u>*<w> " << finl;
    fic << "# (9)  : <v.w>-<v>*<w> " << finl;
    fic << "# (10) : <rho> " << finl;
    fic << "# (11) : <mu> " << finl;

    fic << "#  " << finl;
    fic << "# Nombre minimum de points par segment (X,Z)=cst ayant servi au calcul de la moyenne spatiale : " << min_array(compt_tot) << finl;
    fic << "#  " << finl;
  for (int l=0;l<X_tot.size();l++)
  {
  fic << " X = " << X_tot(l) << finl;
    for (i=0;i<Z_tot.size();i++)
  {
  u    = val_moy(l,i,0,k)/dt;
  v    = val_moy(l,i,1,k)/dt;
  w    = val_moy(l,i,2,k)/dt;
  u2   = val_moy(l,i,3,k)/dt;
  v2   = val_moy(l,i,4,k)/dt;
  w2   = val_moy(l,i,5,k)/dt;
  uv   = val_moy(l,i,6,k)/dt;
  uw   = val_moy(l,i,7,k)/dt;
  vw   = val_moy(l,i,8,k)/dt;
  rho  = val_moy(l,i,10,k)/dt;
  mu   = val_moy(l,i,11,k)/dt;


  fic << Z_tot(i) << "    "  ;
  fic << u << "    " << v << "    " << w << "    " ;
  // Pour eviter NAN, on prend le dmax(,0):
  fic << sqrt(dmax(0,u2-u*u)) << "    " << sqrt(dmax(v2-v*v,0)) << "    " << sqrt(dmax(w2-w*w,0)) << "    "  ;
  fic << u*v-uv << "    " << u*w-uw << "    " << v*w-vw  << "    " ;


  fic << rho << "    " << mu << "    " ;
  fic << finl;
  }
  }
    fic.flush();
    fic.close();
    */
}

void Traitement_particulier_NS_surface::ecriture_fichiers_moy_Temp_old(const DoubleTab& val_moy, const Nom& fichier, const double& dt, const int& k) const
{
  /*
      SFichier fic (fichier);
      fic.setf(ios::scientific);
      double tps = mon_equation->inconnue().temps();
      double u,v,w,T,T2,uT,vT,wT;
      int i;

      fic << "# Temps : " << tps << finl;
      fic << " " << finl;

      fic << "# (0)  : Z " << finl;
      fic << "# (1)  : <T> " << finl;
      fic << "# (2)  : sqrt(<T.T>-<T>*<T>) " << finl;
      fic << "# (3)  : <u.T>-<u>*<T> " << finl;
      fic << "# (4)  : <v.T>-<v>*<T> " << finl;
      fic << "# (5)  : <w.T>-<w>*<T> " << finl;

      fic << "#  " << finl;
      fic << "# Nombre minimum de points par segment (X,Z)=cst ayant servi au calcul de la moyenne spatiale : " << min_array(compt_tot) << finl;
      fic << "#  " << finl;
  for (int l=0;l<X_tot.size();l++)
   {
   fic << " X = " << X_tot(l) << finl;
      for (i=0;i<Z_tot.size();i++)
  {
    u    = val_moy(l,i,0,k)/dt;
    v    = val_moy(l,i,1,k)/dt;
    w    = val_moy(l,i,2,k)/dt;
    T    = val_moy(l,i,13,k)/dt;
    T2   = val_moy(l,i,14,k)/dt;
    uT   = val_moy(l,i,15,k)/dt;
    vT   = val_moy(l,i,16,k)/dt;
    wT   = val_moy(l,i,17,k)/dt;


    fic << Z_tot(i) << "    "  ;
    // Pour eviter NAN, on prend le dmax(,0):
    fic << T << "    " << sqrt(dmax(T2-T*T,0)) << "    "  ;
    fic << u*T-uT << "    " << v*T-vT << "    " << w*T-wT << "    ";

    fic << finl;
    }
  }
      fic.flush();
      fic.close();

  */

}
