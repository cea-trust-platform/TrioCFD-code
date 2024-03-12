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
//////////////////////////////////////////////////////////////////////////////
//
// File:        Traitement_particulier_NS_canal.cpp
// Directory:   $TRIO_U_ROOT/ThHyd/Turbulence
// Version:     /main/39
//
//////////////////////////////////////////////////////////////////////////////

#include <Traitement_particulier_NS_canal.h>
#include <LecFicDistribue.h>
#include <EcrFicCollecte.h>
#include <Navier_Stokes_Turbulent.h>
#include <Probleme_base.h>
#include <Discretisation_base.h>
#include <Domaine_VF.h>
#include <Fluide_Quasi_Compressible.h>
#include <Domaine_Cl_dis_base.h>
#include <Dirichlet_paroi_fixe.h>
#include <Schema_Temps_base.h>
#include <DoubleTrav.h>
#include <communications.h>

Implemente_base(Traitement_particulier_NS_canal,"Traitement_particulier_NS_canal",Traitement_particulier_NS_base);


/*! @brief
 *
 * @param (Sortie& is) un flot de sortie
 * @return (Sortie&) le flot de sortie modifie
 */
Sortie& Traitement_particulier_NS_canal::printOn(Sortie& is) const
{
  return is;
}


/*! @brief
 *
 * @param (Entree& is) un flot d'entree
 * @return (Entree&) le flot d'entree modifie
 */
Entree& Traitement_particulier_NS_canal::readOn(Entree& is)
{
  return is;
}

Entree& Traitement_particulier_NS_canal::lire(Entree& is)
{
  // Initialisation

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

  Ny=200; //modified AT 5/06/09
  Nval=12;
  Nphase=1;
  Nb_ech_phase.resize(1);
  Nb_ech_phase(0)=1;


  if(!sub_type(Navier_Stokes_std,mon_equation.valeur()))
    {
      Cerr << " Traitement_particulier_canal has to be called from a Navier_Stokes problem " << finl;
      Process::exit();
    }

  if(sub_type(Navier_Stokes_Turbulent,mon_equation.valeur()))
    {
      oui_profil_nu_t = 1 ;
      Nval=13 ;
    }

  if (mon_equation.valeur().probleme().nombre_d_equations()>1)
    try
      {
        Temp = mon_equation.valeur().probleme().equation(1).get_champ("temperature_qc");
        oui_profil_Temp = 1 ;
        Nval=18;
        Nval=35; // Modification par YB (29/06/09)
      }
    catch (Champs_compris_erreur)
      {
      }
  // FIN Initialisation

  Motcle accouverte = "{" , accfermee = "}" ;
  Motcle motbidon, motlu;
  is >> motbidon ;
  if (motbidon == accouverte)
    {
      Motcles les_mots(7);
      les_mots[0] = "dt_impr_moy_spat";
      les_mots[1] = "dt_impr_moy_temp";
      les_mots[2] = "debut_stat";
      les_mots[3] = "fin_stat";
      les_mots[4] = "reprise";
      les_mots[5] = "pulsation_w";
      les_mots[6] = "nb_points_par_phase";

      is >> motlu;
      while(motlu != accfermee)
        {
          int rang=les_mots.search(motlu);
          switch(rang)
            {
            case 0 :
              {
                is >> dt_impr_moy_spat;      // intervalle de temps de sorties des moyennes spatiales
                Cerr << "Spatial averages are printed for : dt_impr_moy_spat = " << dt_impr_moy_spat << finl;
                break;
              }
            case 1 :
              {
                is >> dt_impr_moy_temp;      // intervalle de temps de sorties des moyennes temporelles
                Cerr << "Temporal averages are printed for : dt_impr_moy_temp = " << dt_impr_moy_temp << finl;
                break;
              }
            case 2 :
              {
                is >> temps_deb;	     // temps de debut des moyennes temporelles
                Cerr << "Temporal averages start at : temps_deb = " << temps_deb << finl;
                break;
              }
            case 3 :
              {
                is >> temps_fin;	     // temps de debut des moyennes temporelles
                Cerr << "Temporal averages finish at : temps_fin = " << temps_fin << finl;
                break;
              }
            case 4 :
              {
                oui_repr = 1;
                is  >> fich_repr ; 	    // indication du nom du fichier de reprise des stats
                Cerr << "The time statistics file is : " << fich_repr << finl;
                break;
              }
            case 5 :
              {
                oui_pulse = 1;
                is  >> w ;  		   // pulsation
                freq = w/(2.*3.141592653); // frequence associee
                break;
              }
            case 6 :
              {
                oui_pulse = 1;
                is  >> Nphase ;  	  // nombre de points pour decrire une phase
                break;
              }
            default :
              {
                Cerr << motlu << " is not a keyword for Traitement_particulier_NS_canal " << finl;
                Cerr << "Since the Trio_U v1.5, the syntax of Canal option has changed." << finl;
                Cerr << "Check the 1.5 version (or later) of the user's manual." << finl;
                Process::exit();
              }
            }
          is >> motlu;
        }
      is >> motlu;
      if (motlu != accfermee)
        {
          Cerr << "Error while reading in Traitement_particulier_NS_canal" << finl;;
          Cerr << "We were expecting a }" << finl;
          Process::exit();
        }
    }
  else
    {
      Cerr << "Error while reading in Traitement_particulier_NS_canal" << finl;
      Cerr << "We were expecting a {" << finl;
      Process::exit();
    }

  return is;
}

void Traitement_particulier_NS_canal::remplir_Tab_recap(IntTab& Tab_rec) const
{
// surchargeee dans VDF
  Cerr << "Traitement_particulier_NS_canal::remplir_Tab_recap ne marche pas pour le VEF" << finl;
  Process::exit();
}

void Traitement_particulier_NS_canal::preparer_calcul_particulier()
{
  remplir_Y(Y,compt,Ny); // renvoie vers Traitement_particulier_NS_canal_VDF ou Traitement_particulier_NS_canal_VEF
  remplir_Tab_recap(Tab_recap);
  remplir_reordonne_Y_tot(Y,Y_tot);

  int NN=Y_tot.size();
  compt_tot.resize(NN);
  val_moy_tot.resize(NN,Nval,1);

  if (oui_repr!=1)
    {
      val_moy_temp.resize(Y_tot.size(),Nval,1);
      val_moy_temp=0.;

      if (oui_pulse==1)
        {
          val_moy_phase.resize(Y_tot.size(),Nval,Nphase);
          val_moy_phase=0.;
          Nb_ech_phase.resize(Nphase);
          Nb_ech_phase=0;
        }
    }
}

void Traitement_particulier_NS_canal::remplir_reordonne_Y_tot(const DoubleVect& Yc, DoubleVect& Y_total) const
{
  DoubleVect Y_p(Yc);
  int j_tot=1;

  envoyer(Y_p, Process::me(), 0, Process::me());

  if(je_suis_maitre())
    {
      Y_total.resize(1);
      Y_total(0) = Yc(0);

      for(int p=0; p<Process::nproc(); p++)
        {
          recevoir(Y_p,p,0,p);


          for (int k=0; k<Y_p.size(); k++)
            {
              if(!est_egal(Y_p(k),-100.)) // on recherche si Y_p(k) est deja contenu dans Y_tot
                {
                  int ok_new=1;

                  for (int j=0; j<Y_total.size(); j++)
                    {
                      if (est_egal(Y_total(j),Y_p(k))) ok_new=0;
                    }

                  if(ok_new==1)
                    {
                      Y_total.resize(j_tot+1);
                      Y_total(j_tot) = Y_p(k);
                      j_tot++;
                    }
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
}

void Traitement_particulier_NS_canal::reprendre_stat()
{
  if (je_suis_maitre())
    {
      Cerr << "Traitement_particulier_NS_canal::reprendre_stat!!" << flush << finl;

      if (oui_repr==1)
        {
          reprendre_stat_canal(val_moy_temp,fich_repr);
          Nval=val_moy_temp.dimension(1);
        }

      if (oui_repr==1 && oui_pulse==1)
        {
          Nom fich_repr_phase;
          fich_repr_phase  = fich_repr;
          fich_repr_phase +="_phase";

          reprendre_stat_canal(val_moy_phase,fich_repr_phase);
          Nphase=val_moy_phase.dimension(2);
        }
    }
}

void Traitement_particulier_NS_canal::reprendre_stat_canal(DoubleTab& val, const Nom& fichier)
{
  double tps,tps_deb_moy;
  int NN,Nv,Np;
  Nom temps;

  EFichier fic(fichier);
  fic.setf(ios::scientific);

  if(fic.fail())
    {
      Cerr << "Impossible d'ouvrir le fichier " << fichier << finl;
      Process::exit();
    }

  fic >> tps;
  fic >> tps_deb_moy;
  fic >> debut_phase;
  fic >> ind_phase;
  fic >> NN;
  fic >> Nv;
  fic >> Np;

  val.resize(NN,Nv,Np);
  Nb_ech_phase.resize(Np);

  for (int k=0; k<Np; k++)
    fic >> Nb_ech_phase(k);

  for (int k=0; k<Np; k++)
    for (int j=0; j<Nv; j++)
      for (int i=0; i<NN; i++)
        fic >> val(i,j,k);

  if(!est_egal(tps_deb_moy,temps_deb,0.1))
    {
      Cerr << "ERREUR : le temps de debut des stats figurant dans le jeu de donnees differe de beaucoup" << finl;
      Cerr << "du temps defini dans le fichier de reprise : " << fichier << finl;
      Process::exit();
    }

  temps_deb=tps_deb_moy;
}

void Traitement_particulier_NS_canal::sauver_stat() const
{
  double tps = mon_equation->inconnue().temps();

  if (je_suis_maitre() && (tps>=temps_deb) && (tps<=temps_fin) )
    {
      Cerr << "Traitement_particulier_NS_canal::sauver_stat!!" << flush << finl;

      Nom temps = Nom(tps);
      Nom fich_sauv1 ="val_moy_temp_";
      fich_sauv1+=temps;
      fich_sauv1+=".sauv";
      Nom fich_sauv2 = fich_sauv1;
      fich_sauv2+="_phase";

      sauver_stat_canal(val_moy_temp,fich_sauv1);
      if(oui_pulse==1)  sauver_stat_canal(val_moy_phase,fich_sauv2);
    }
}

void Traitement_particulier_NS_canal::sauver_stat_canal(const DoubleTab& val, const Nom& fichier) const
{
  double tps = mon_equation->inconnue().temps();
  int NN,Nv,Np;
  NN=val.dimension(0);
  Nv=val.dimension(1);
  Np=val.dimension(2);

  SFichier fic(fichier);
  fic.setf(ios::scientific);

  fic << tps << finl;
  fic << temps_deb << finl;
  fic << debut_phase << finl;
  fic << ind_phase << finl;
  fic << NN << finl;
  fic << Nv << finl;
  fic << Np << finl;

  for (int k=0; k<Np; k++)
    fic <<Nb_ech_phase(k) << finl;

  for (int k=0; k<Np; k++)
    for (int j=0; j<Nv; j++)
      for (int i=0; i<NN; i++)
        fic << val(i,j,k) << finl;
}

void Traitement_particulier_NS_canal::post_traitement_particulier()
{

  // Calcul des Moyennes spatiales
  //////////////////////////////////////////////////////////

  DoubleTrav val_moy(Ny,Nval);

  val_moy=0.;

  calculer_moyenne_spatiale_vitesse_rho_mu(val_moy);

  if (oui_profil_nu_t != 0)
    calculer_moyenne_spatiale_nut(val_moy);

  if (oui_profil_Temp != 0)
    calculer_moyenne_spatiale_Temp(val_moy);


  // Echange des donnees entre processeurs
  //////////////////////////////////////////////////////////

  DoubleVect Y_p(Y);
  DoubleVect compt_p(compt);
  DoubleTab val_moy_p(val_moy);

  const int NN = Y_tot.size();

  envoyer(Y_p,Process::me(),0,Process::me());
  envoyer(compt_p,Process::me(),0,Process::me());
  envoyer(val_moy_p,Process::me(),0,Process::me());

  if(je_suis_maitre())
    {
      compt_tot=0;
      val_moy_tot=0.;
      int j;

      for(int p=0; p<Process::nproc(); p++)
        {
          recevoir(Y_p,p,0,p);
          recevoir(compt_p,p,0,p);
          recevoir(val_moy_p,p,0,p);

          int Y_p_size = Y_p.size();
          for (int i=0; i<Y_p_size; i++)
            {

              for (j=0; j<NN; j++)
                if(est_egal(Y_tot[j],Y_p[i])) break;

              compt_tot(j) += compt_p(i);

              for (int k=0; k<Nval; k++)
                val_moy_tot(j,k,0) += val_moy_p(i,k);


              //   if(est_egal(Y_tot[32],Y_p[i])) Cerr  << " proc " << p <<" flux reel par proc " << val_moy_p(i,33) << " flux sm par proc " << val_moy_p(i,34) << finl;

            }
        }

      // Cerr <<" flux reel " << val_moy_tot(32,33,0) << " flux sm " << val_moy_tot(32,34,0) << finl;


      for (int i=0; i<NN; i++)
        for (int k=0; k<Nval; k++)
          val_moy_tot(i,k,0) /= compt_tot[i];



    }


  // Calcul et ecriture des grandeurs parietales
  ///////////////////////////////////////////////////////////////////

  if(je_suis_maitre())
    {
      calcul_reynolds_tau();
    }

  // Calcul des Moyennes temporelles
  //////////////////////////////////////////////////////////////////

  double tps = mon_equation->inconnue().temps();

  if(je_suis_maitre())
    {
      if ((tps>=temps_deb)&&(tps<=temps_fin))
        {
          static int init_stat_temps = 0;

          double dt_v = mon_equation->schema_temps().pas_de_temps();

          if(init_stat_temps==0 && oui_repr!=1) // sinon, les valeurs de val_moy_temp ont ete lues a partir de reprendre_stat()
            {
              temps_deb = tps-dt_v;

              val_moy_temp.resize(NN,Nval,1);
              val_moy_temp=0.;

              init_stat_temps++;
            }

          for (int j=0; j<Nval; j++)
            for (int i=0; i<NN; i++)
              {
                if ((j!=3) && (j!=4) && (j!=5) && (j!=6) && (j!=7) && (j!=8) && (j!=14) && (j!=15) && (j!=16) && (j!=17))
                  val_moy_temp(i,j,0)+= dt_v*val_moy_tot(i,j,0);

                if (j==3)
                  val_moy_temp(i,j,0)+= dt_v*(val_moy_tot(i,3,0)-val_moy_tot(i,0,0)*val_moy_tot(i,0,0));
                if (j==4)
                  val_moy_temp(i,j,0)+= dt_v*(val_moy_tot(i,4,0)-val_moy_tot(i,1,0)*val_moy_tot(i,1,0));
                if (j==5)
                  val_moy_temp(i,j,0)+= dt_v*(val_moy_tot(i,5,0)-val_moy_tot(i,2,0)*val_moy_tot(i,2,0));

                if (j==6)
                  val_moy_temp(i,j,0)+= dt_v*(val_moy_tot(i,6,0)-val_moy_tot(i,0,0)*val_moy_tot(i,1,0));
                if (j==7)
                  val_moy_temp(i,j,0)+= dt_v*(val_moy_tot(i,7,0)-val_moy_tot(i,0,0)*val_moy_tot(i,2,0));
                if (j==8)
                  val_moy_temp(i,j,0)+= dt_v*(val_moy_tot(i,8,0)-val_moy_tot(i,1,0)*val_moy_tot(i,2,0));


                if (j==14)
                  //val_moy_temp(i,j,0)+= dt_v*val_moy_tot(i,13,0)*val_moy_tot(i,13,0);
                  val_moy_temp(i,j,0)+= dt_v*(val_moy_tot(i,14,0)-val_moy_tot(i,13,0)*val_moy_tot(i,13,0));
                if (j==15)
                  val_moy_temp(i,j,0)+= dt_v*(val_moy_tot(i,15,0)-val_moy_tot(i,0,0)*val_moy_tot(i,13,0));
                if (j==16)
                  val_moy_temp(i,j,0)+= dt_v*(val_moy_tot(i,16,0)-val_moy_tot(i,1,0)*val_moy_tot(i,13,0));
                if (j==17)
                  val_moy_temp(i,j,0)+= dt_v*(val_moy_tot(i,17,0)-val_moy_tot(i,2,0)*val_moy_tot(i,13,0));

              }

        }
    }

  // Calcul des Moyennes de phases
  //////////////////////////////////////////////////////////////////

  if(je_suis_maitre() && oui_pulse==1)
    {
      if ((tps>=temps_deb)&&(tps<=temps_fin))
        {
          double dt_v = mon_equation->schema_temps().pas_de_temps();

          if ( (cos(w*tps) > cos(w*(tps+dt_v))) && (cos(w*tps) > cos(w*(tps-dt_v))) ) // debut d'une periode
            {
              debut_phase=tps;
              ind_phase=1;
            }

          if(ind_phase==1)
            {
              for(int k=0; k<Nphase; k++)
                {
                  double tps_k=debut_phase+k/(Nphase*freq);

                  if( (tps > tps_k-0.5*dt_v)  && (tps < tps_k+0.5*dt_v) ) // recherche de la k-phase correspondant
                    {
                      for (int j=0; j<Nval; j++)
                        for (int i=0; i<NN; i++)
                          val_moy_phase(i,j,k) += val_moy_tot(i,j,0);

                      Nb_ech_phase(k)++;

                      if(k==Nphase-1)  ind_phase=2; // marqueur pour ecriture de moyennes de phase (voir plus bas)
                    }
                }
            }
        }
    }

  // Ecriture des Moyennes spatiales et temporelles
  ///////////////////////////////////////////////////////////////

  if(je_suis_maitre())
    {
      double stationnaire_atteint=mon_equation->schema_temps().stationnaire_atteint();
      double tps_passe=mon_equation->schema_temps().temps_courant();
      double nb_pas_dt_max=mon_equation->schema_temps().nb_pas_dt_max();
      double nb_pas_dt=mon_equation->schema_temps().nb_pas_dt();
      double temps_max=mon_equation->schema_temps().temps_max();
      int impr_inst;
      if (dt_impr_moy_spat<=(tps-tps_passe))
        impr_inst=1;
      else
        {
          double epsilon = 1.e-8;
          int i=(int) (tps/dt_impr_moy_spat + epsilon);
          int j=(int) (tps_passe/dt_impr_moy_spat + epsilon);
          impr_inst=(i>j);
        }
      int impr_stat;
      if (dt_impr_moy_temp<=(tps-tps_passe))
        impr_stat=1;
      else
        {
          double epsilon = 1.e-8;
          int i=(int) (tps/dt_impr_moy_temp + epsilon);
          int j=(int) (tps_passe/dt_impr_moy_temp + epsilon);
          impr_stat=(i>j);
        }

      // sauvegarde periodique des moyennes spatiales

      if ((nb_pas_dt+1<=1) || impr_inst || (temps_max <= tps) || (nb_pas_dt_max <= nb_pas_dt+1) || stationnaire_atteint)
        {
          double dt = 1.;

          Nom fichier1 = "Moyennes_spatiales_vitesse_rho_mu_";
          Nom fichier2 = "Moyennes_spatiales_nut_";
          Nom fichier3 = "Moyennes_spatiales_Temp_";
          Nom temps = Nom(tps);

          fichier1 += temps;
          fichier2 += temps;
          fichier3 += temps;

          /*
            			        ecriture_fichiers_moy_vitesse_rho_mu(val_moy_tot,fichier1,dt,0);
                if (oui_profil_nu_t == 1)   ecriture_fichiers_moy_nut(val_moy_tot,fichier2,dt,0);
                if (oui_profil_Temp == 1)   ecriture_fichiers_moy_Temp(val_moy_tot,fichier3,dt,0);
          */

          ecriture_fichiers_moy_vitesse_rho_mu_old(val_moy_tot,fichier1,dt,0);
          if (oui_profil_nu_t == 1)   ecriture_fichiers_moy_nut(val_moy_tot,fichier2,dt,0);
          if (oui_profil_Temp == 1)   ecriture_fichiers_moy_Temp_old(val_moy_tot,fichier3,dt,0);
        }

      // sauvegarde periodique des moyennes temporelles

      if ((tps>temps_deb) && ((nb_pas_dt+1<=1) || impr_stat || (temps_max <= tps) || (nb_pas_dt_max <= nb_pas_dt+1) || stationnaire_atteint))
        {
          double dt = tps-temps_deb;

          Nom fichier1 = "Moyennes_temporelles_vitesse_rho_mu_";
          Nom fichier2 = "Moyennes_temporelles_nut_";
          Nom fichier3 = "Moyennes_temporelles_Temp_"; // modifier AT 5/06/09
          Nom temps = Nom(tps);

          fichier1 += temps;
          fichier2 += temps;
          fichier3 += temps;

          ecriture_fichiers_moy_vitesse_rho_mu(val_moy_temp,fichier1,dt,0);
          if (oui_profil_nu_t == 1)   ecriture_fichiers_moy_nut(val_moy_temp,fichier2,dt,0);
          if (oui_profil_Temp == 1)   ecriture_fichiers_moy_Temp(val_moy_temp,fichier3,dt,0);
        }

      // sauvegarde en continue des moyennes temporelles

      if (tps>temps_deb)
        {
          double dt = tps-temps_deb;

          Nom fichier1 = "Moyennes_temporelles_vitesse_rho_mu";
          Nom fichier2 = "Moyennes_temporelles_nut";
          Nom fichier3 = "Moyennes_temporelles_Temp"; // modifier AT 5/06/09

          ecriture_fichiers_moy_vitesse_rho_mu(val_moy_temp,fichier1,dt,0);
          if (oui_profil_nu_t == 1)   ecriture_fichiers_moy_nut(val_moy_temp,fichier2,dt,0);
          if (oui_profil_Temp == 1)   ecriture_fichiers_moy_Temp(val_moy_temp,fichier3,dt,0);
        }

      // sauvegarde en continue des moyennes de phase (en fin de cycle)

      if ( (ind_phase==2) && (tps>temps_deb) )
        {
          ind_phase=0;

          for(int k=0; k<Nphase; k++)
            {
              double dt = Nb_ech_phase(k);

              Nom fichier1 = "Moyennes_de_phase_vitesse_rho_mu_";
              Nom fichier2 = "Moyennes_de_phase_nut_";
              Nom fichier3 = "Moyennes_de_phase_Temp_";
              Nom phase = Nom(k);

              fichier1 += phase;
              fichier2 += phase;
              fichier3 += phase;


              ecriture_fichiers_moy_vitesse_rho_mu_old(val_moy_phase,fichier1,dt,k);
              if (oui_profil_nu_t == 1)   ecriture_fichiers_moy_nut(val_moy_phase,fichier2,dt,k);
              if (oui_profil_Temp == 1)   ecriture_fichiers_moy_Temp_old(val_moy_phase,fichier3,dt,k);
            }
        }
    }
}

void Traitement_particulier_NS_canal::calcul_reynolds_tau()
{

  double tps = mon_equation->inconnue().temps();

  // Aspects geometriques du maillage
  /////////////////////////////////////////////////////////

// !!!!!  Hypotheses : maillage symetrique suivant la demi-hauteur et s'etendant de Y=0 a Y=H

  Equation_base& eqn=mon_equation.valeur();
  const Discretisation_base& discr=eqn.discretisation();
  Nom nom_discr=discr.que_suis_je();
  int k=0;
  if(nom_discr=="VEFPreP1B" || nom_discr=="VEF") k=1;

  int NN = Y_tot.size();
  int kmin=k;		    // indice du premier point hors paroi
  int kmax=NN-1-k;		    // indice du dernier point hors paroi
  double ymin=Y_tot(kmin);   	    // position du premier point
  double ymax=Y_tot(kmax);   	    // position du dernier point
  double hs2=0.5*(ymin+ymax);       // demi-hauteur


  // definition des grandeurs interveannt dans la suite des calculs
  //////////////////////////////////////////////////////

  double mu_bas, mu_haut;   	 // viscosite dynamique
  double rho_bas, rho_haut;	 // masse volumique
  double utang_bas, utang_haut;	 // norme de la vitesse tangente a la paroi
  double tauwb, tauwh, tauwm;	 // cisaillement a la paroi
  double utaub, utauh, utaum;	 // vitesse de frottement
  double retaub, retauh, retaum; // Reynolds de frottement

  mu_bas  = val_moy_tot(kmin,11,0);
  mu_haut = val_moy_tot(kmax,11,0);

  rho_bas  = val_moy_tot(kmin,10,0);
  rho_haut = val_moy_tot(kmax,10,0);

  utang_bas  = val_moy_tot(kmin,9,0);
  utang_haut = val_moy_tot(kmax,9,0);


  // calcul du cisaillement a la paroi suivant la vitesse tangentielle
  //////////////////////////////////////////////////////

  // approximation lineaire : Tau_w = mu * ||u_t||(y1) / dy1

  tauwb= mu_bas * utang_bas  / ymin;
  tauwh= mu_haut* utang_haut / ymin;


  // calcul du cisaillement a la paroi si loi de paroi
  //////////////////////////////////////////////////////

  if (sub_type(Navier_Stokes_Turbulent,mon_equation.valeur()))
    {
      const Navier_Stokes_Turbulent& eqn_hydr = ref_cast(Navier_Stokes_Turbulent,mon_equation.valeur() ) ;
      const Mod_turb_hyd& mod_turb = eqn_hydr.modele_turbulence();
      const Turbulence_paroi& loipar = mod_turb.loi_paroi();

      Nom type_loi = loipar.valeur().que_suis_je();

      if ( !(type_loi.debute_par("negligeable")) )
        {
          tauwb=0.;
          tauwh=0.;

          // PQ : 13/07/05 : prise en compte des lois de paroi pour le calcul de u_tau

          // Hypotheses : 1ere condition de Dirichlet = paroi basse
          //	 	2eme condition de Dirichlet = paroi haute
          //	 	maillage regulier suivant x
          //	 	calcul sequenciel

          if(Process::nproc()>1) // PQ : 11/06/07
            {
              Cerr << " La prise en compte de la loi de paroi pour le calcul de u_tau paroi haute et u_tau paroi basse " << finl;
              Cerr << " n'est pas parallelise - A FAIRE " << finl;
              Cerr << " autre possibilite : mettre les parois sur le meme proc" << finl;
              Cerr << " et commenter la presente verif dans Traitement_particulier_NS_canal::calcul_reynolds_tau   " << finl;
              Process::exit();
            }

          const Domaine_Cl_dis_base& domaine_Cl_dis_base = ref_cast(Domaine_Cl_dis_base,eqn_hydr.domaine_Cl_dis().valeur());

          const Conds_lim& les_cl = domaine_Cl_dis_base.les_conditions_limites();
          int nb_cl=les_cl.size();
          int num_cl,fac;
          int num_cl_rep=0;

          DoubleTab tau_tan;
          tau_tan.ref(loipar->Cisaillement_paroi());

          for (num_cl=0; num_cl<nb_cl; num_cl++)
            {
              //Boucle sur les bords
              const Cond_lim& la_cl = les_cl[num_cl];
              const Front_VF& la_front_dis = ref_cast(Front_VF,la_cl.frontiere_dis());
              int nbfaces = la_front_dis.nb_faces();
              int ndeb = la_front_dis.num_premiere_face();
              int nfin = ndeb + nbfaces;

              if (sub_type(Dirichlet_paroi_fixe,la_cl.valeur()))
                {
                  if (Objet_U::dimension == 2 )
                    {
                      for (fac=ndeb; fac<nfin ; fac++)
                        {
                          tauwb+=sqrt(tau_tan(fac,0)*tau_tan(fac,0));
                        }
                    }
                  else
                    {
                      for (fac=ndeb; fac<nfin ; fac++)
                        {
                          tauwb+=sqrt(tau_tan(fac,0)*tau_tan(fac,0)+tau_tan(fac,2)*tau_tan(fac,2));
                        }
                    }

                  tauwb/=nbfaces;

                  num_cl_rep=num_cl+1;
                  break;
                }
            } //Boucle sur les bords

          for (num_cl=num_cl_rep; num_cl<nb_cl; num_cl++)
            {
              //Boucle sur les bords
              const Cond_lim& la_cl = les_cl[num_cl];
              const Front_VF& la_front_dis = ref_cast(Front_VF,la_cl.frontiere_dis());
              int nbfaces = la_front_dis.nb_faces();
              int ndeb = la_front_dis.num_premiere_face();
              int nfin = ndeb + nbfaces;

              if (sub_type(Dirichlet_paroi_fixe,la_cl.valeur()))
                {
                  if (Objet_U::dimension == 2 )
                    {
                      for (fac=ndeb; fac<nfin ; fac++)
                        {
                          tauwh+=sqrt(tau_tan(fac,0)*tau_tan(fac,0));
                        }
                    }
                  else
                    {
                      for (fac=ndeb; fac<nfin ; fac++)
                        {
                          tauwh+=sqrt(tau_tan(fac,0)*tau_tan(fac,0)+tau_tan(fac,2)*tau_tan(fac,2));
                        }
                    }

                  tauwh/=nbfaces;

                  break;
                }
            } //Boucle sur les bords


          if(!(mon_equation.valeur().probleme().is_QC()))
            {
              tauwb*=rho_bas;
              tauwh*=rho_haut;
            }

        }// if (!Paroi_negligeable)
    }// if(sub_type(Navier_Stokes_Turbulent,mon_equation.valeur()))


  // calcul et ecritures des differentes grandeurs parietales
  //////////////////////////////////////////////////////

  utaub=sqrt(tauwb/rho_bas);
  if(val_moy_tot(kmin,0,0)<=0) utaub*=-1.;
  utauh=sqrt(tauwh/rho_haut);
  if(val_moy_tot(kmax,0,0)<=0) utauh*=-1.;
  retaub=rho_bas*utaub*hs2/mu_bas;
  retauh=rho_haut*utauh*hs2/mu_haut;
  tauwm=(tauwh+tauwb)/2.;
  retaum=(retauh+retaub)/2.;
  utaum=(utauh+utaub)/2.;

  SFichier fic1("reynolds_tau.dat",ios::app);
  SFichier fic2("u_tau.dat",ios::app);
  SFichier fic3("tauw.dat",ios::app);

  fic1 << tps << " " << retaum << " " << retaub << " " << retauh << finl;
  fic2 << tps << " " << utaum << " "  << utaub << " "  << utauh  << finl;
  fic3 << tps << " " << tauwm << " "  << tauwb << " "  << tauwh  << finl;
  fic1<<flush;
  fic1.close();
  fic2<<flush;
  fic2.close();
  fic3<<flush;
  fic3.close();
}


void Traitement_particulier_NS_canal::ecriture_fichiers_moy_vitesse_rho_mu(const DoubleTab& val_moy, const Nom& fichier, const double& dt, const int& k) const
{
  SFichier fic (fichier);

  double tps = mon_equation->inconnue().temps();
  double u,v,wl,u2,v2,w2,uv,uw,vw,rho,mu;
  int i;

  fic << "# Temps : " << tps << finl;
  fic << " " << finl;

  fic << "# (0)  : Y " << finl;
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
  fic << "# Nombre minimum de points par plan Y=cste ayant servi au calcul de la moyenne spatiale : " << min_array(compt_tot) << finl;
  fic << "#  " << finl;

  for (i=0; i<Y_tot.size(); i++)
    {
      u    = val_moy(i,0,k)/dt;
      v    = val_moy(i,1,k)/dt;
      wl    = val_moy(i,2,k)/dt;
      u2   = val_moy(i,3,k)/dt;
      v2   = val_moy(i,4,k)/dt;
      w2   = val_moy(i,5,k)/dt;
      uv   = val_moy(i,6,k)/dt;
      uw   = val_moy(i,7,k)/dt;
      vw   = val_moy(i,8,k)/dt;
      rho  = val_moy(i,10,k)/dt;
      mu   = val_moy(i,11,k)/dt;


      fic << Y_tot(i) << "    "  ;
      fic << u << "    " << v << "    " << wl << "    " ;
      // Pour eviter NAN, on prend le dmax(,0):
      fic << sqrt(dmax(0,u2)) << "    " << sqrt(dmax(v2,0)) << "    " << sqrt(dmax(w2,0)) << "    "  ;
      fic << -uv << "    " << -uw << "    " << -vw  << "    " ;
      fic << rho << "    " << mu << "    " ;
      fic << finl;
    }
  fic.flush();
  fic.close();
}

void Traitement_particulier_NS_canal::ecriture_fichiers_moy_nut(const DoubleTab& val_moy, const Nom& fichier, const double& dt, const int& k) const
{
  SFichier fic (fichier);

  double tps = mon_equation->inconnue().temps();
  int i;

  fic << "# Temps : " << tps << finl;
  fic << " " << finl;

  fic << "# (0)  : Y " << finl;
  fic << "# (1)  : <nu_t> " << finl;

  fic << "#  " << finl;
  fic << "# Nombre minimum de points par plan Y=cste ayant servi au calcul de la moyenne spatiale : " << min_array(compt_tot) << finl;
  fic << "#  " << finl;

  for (i=0; i<Y_tot.size(); i++)
    fic << Y_tot(i) << "    " << val_moy(i,12,k)/dt << finl;

  fic.flush();
  fic.close();
}

void Traitement_particulier_NS_canal::ecriture_fichiers_moy_Temp(const DoubleTab& val_moy, const Nom& fichier, const double& dt, const int& k) const
{
  SFichier fic (fichier);

  double tps = mon_equation->inconnue().temps();
  double T,T2,uT,vT,wT,kappa,rhoT,rhouT,rhovT,rhowT,rhou,rhov,rhow,rhouu,rhovv,rhoww,rhouv,rhouw,rhovw,gradT,fluxdiff,fluxsm; //modified YB 30/06/09
  int i;

  fic << "# Temps : " << tps << finl;
  fic << " " << finl;

  fic << "# (0)  : Y " << finl;
  fic << "# (1)  : <T> " << finl;
  fic << "# (2)  : sqrt(<T.T>-<T>*<T>) " << finl;
  fic << "# (3)  : <u.T>-<u>*<T> " << finl;
  fic << "# (4)  : <v.T>-<v>*<T> " << finl;
  fic << "# (5)  : <w.T>-<w>*<T> " << finl;
  //modified AT 5/06/09
  fic << "# (6)  : <kappa_sgs> " << finl;
  fic << "# (7)  : <rho.T> " << finl;
  fic << "# (8)  : <rho.u.T> " << finl;
  fic << "# (9)  : <rho.v.T> " << finl;
  fic << "# (10)  : <rho.w.T> " << finl;
  fic << "# (11)  : <rho.u> " << finl;
  fic << "# (12)  : <rho.v> " << finl;
  fic << "# (13)  : <rho.w> " << finl;
  fic << "# (14)  : <rho.u.u> " << finl;
  fic << "# (15)  : <rho.v.v> " << finl;
  fic << "# (16)  : <rho.w.w> " << finl;
  fic << "# (17)  : <rho.u.v> " << finl;
  fic << "# (18)  : <rho.u.w> " << finl;
  fic << "# (19)  : <rho.v.w> " << finl;
  // fin modif AT 5/06/09
  fic << "# (20)  : <dT/dy> " << finl;		//modif YB 29/06/09
  fic << "# (21)  : <-lambda.dT/dy> " << finl;	//modif YB 30/06/09
  fic << "# (22)  : <-lambda_sm.dT/dy> " << finl;	//modif YB 30/06/09

  fic << "#  " << finl;
  fic << "# Nombre minimum de points par plan Y=cste ayant servi au calcul de la moyenne spatiale : " << min_array(compt_tot) << finl;
  fic << "#  " << finl;

  for (i=0; i<Y_tot.size(); i++)
    {
      // u    = val_moy(i,0,k)/dt;
      // v    = val_moy(i,1,k)/dt;
      //w    = val_moy(i,2,k)/dt;
      T    = val_moy(i,13,k)/dt;
      T2   = val_moy(i,14,k)/dt;
      uT   = val_moy(i,15,k)/dt;
      vT   = val_moy(i,16,k)/dt;
      wT   = val_moy(i,17,k)/dt;
      // modified AT 5/06/09
      kappa = val_moy(i,18,k)/dt;
      rhoT	= val_moy(i,19,k)/dt;
      rhouT	= val_moy(i,20,k)/dt;
      rhovT	= val_moy(i,21,k)/dt;
      rhowT	= val_moy(i,22,k)/dt;
      rhou	= val_moy(i,23,k)/dt;
      rhov	= val_moy(i,24,k)/dt;
      rhow	= val_moy(i,25,k)/dt;
      rhouu	= val_moy(i,26,k)/dt;
      rhovv	= val_moy(i,27,k)/dt;
      rhoww	= val_moy(i,28,k)/dt;
      rhouv	= val_moy(i,29,k)/dt;
      rhouw	= val_moy(i,30,k)/dt;
      rhovw	= val_moy(i,31,k)/dt;
      // fin modif AT 5/06/09
      gradT = val_moy(i,32,k)/dt;		//modif YB 29/06/09
      fluxdiff = val_moy(i,33,k)/dt;	//modif YB 30/06/09
      fluxsm = val_moy(i,34,k)/dt;	//modif YB 30/06/09

      fic << Y_tot(i) << "    "  ;
      // Pour eviter NAN, on prend le dmax(,0):
      fic << T << "    "  <<sqrt(dmax(T2,0))<<" ";
      fic << -uT << "    " << -vT  << "    " << -wT  << "    ";
      //modified at 5/06/09
      fic << kappa << "    ";
      fic << rhoT << "    ";
      fic << rhouT << "    ";
      fic << rhovT << "    ";
      fic << rhowT << "    ";
      fic << rhou << "    ";
      fic << rhov << "    ";
      fic << rhow << "    ";
      fic << rhouu << "    ";
      fic << rhovv << "    ";
      fic << rhoww << "    ";
      fic << rhouv << "    ";
      fic << rhouw << "    ";
      fic << rhovw << "    ";
      //fin modif AT 5/06/09
      fic << gradT << "    ";		//modif YB 29/06/09
      fic << fluxdiff << "    ";		//modif YB 30/06/09
      fic << fluxsm << "    ";		//modif YB 30/06/09

      fic << finl;
    }
  fic.flush();
  fic.close();
}

//Apres correction de l expression de la moyenne temporelle et des methodes d ecriture des moyennes temporelles
//On conserve ces deux methodes d ecriture en version old pour l ecriture des moyennes de phase

void Traitement_particulier_NS_canal::ecriture_fichiers_moy_vitesse_rho_mu_old(const DoubleTab& val_moy, const Nom& fichier, const double& dt, const int& k) const
{
  SFichier fic (fichier);

  double tps = mon_equation->inconnue().temps();
  double u,v,wl,u2,v2,w2,uv,uw,vw,rho,mu; //modified AT 5/06/09
  int i;

  fic << "# Temps : " << tps << finl;
  fic << " " << finl;

  fic << "# (0)  : Y " << finl;
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
  fic << "# Nombre minimum de points par plan Y=cste ayant servi au calcul de la moyenne spatiale : " << min_array(compt_tot) << finl;
  fic << "#  " << finl;

  for (i=0; i<Y_tot.size(); i++)
    {
      u    = val_moy(i,0,k)/dt;
      v    = val_moy(i,1,k)/dt;
      wl    = val_moy(i,2,k)/dt;
      u2   = val_moy(i,3,k)/dt;
      v2   = val_moy(i,4,k)/dt;
      w2   = val_moy(i,5,k)/dt;
      uv   = val_moy(i,6,k)/dt;
      uw   = val_moy(i,7,k)/dt;
      vw   = val_moy(i,8,k)/dt;
      rho  = val_moy(i,10,k)/dt;
      mu   = val_moy(i,11,k)/dt;


      fic << Y_tot(i) << "    "  ;
      fic << u << "    " << v << "    " << wl << "    " ;
      // Pour eviter NAN, on prend le dmax(,0):
      fic << sqrt(dmax(0,u2-u*u)) << "    " << sqrt(dmax(v2-v*v,0)) << "    " << sqrt(dmax(w2-wl*wl,0)) << "    "  ;
      fic << u*v-uv << "    " << u*wl-uw << "    " << v*wl-vw  << "    " ;


      fic << rho << "    " << mu << "    " ;
      fic << finl;
    }
  fic.flush();
  fic.close();
}



void Traitement_particulier_NS_canal::ecriture_fichiers_moy_Temp_old(const DoubleTab& val_moy, const Nom& fichier, const double& dt, const int& k) const
{
  SFichier fic (fichier);
  fic.setf(ios::scientific);
  double tps = mon_equation->inconnue().temps();
  double u,v,wl,T,T2,uT,vT,wT,kappa,rhoT,rhouT,rhovT,rhowT,rhou,rhov,rhow,rhouu,rhovv,rhoww,rhouv,rhouw,rhovw,gradT,fluxdiff,fluxsm;
  // modified YB 30/06/09

  fic << "# Temps : " << tps << finl;
  fic << " " << finl;

  fic << "# (0)  : Y " << finl;
  fic << "# (1)  : <T> " << finl;
  fic << "# (2)  : sqrt(<T.T>-<T>*<T>) " << finl;
  fic << "# (3)  : <u.T>-<u>*<T> " << finl;
  fic << "# (4)  : <v.T>-<v>*<T> " << finl;
  fic << "# (5)  : <w.T>-<w>*<T> " << finl;
  //modified AT 5/06/09
  fic << "# (6)  : <kappa_sgs> " << finl;
  fic << "# (7)  : <rho.T> " << finl;
  fic << "# (8)  : <rho.u.T> " << finl;
  fic << "# (9)  : <rho.v.T> " << finl;
  fic << "# (10)  : <rho.w.T> " << finl;
  fic << "# (11)  : <rho.u> " << finl;
  fic << "# (12)  : <rho.v> " << finl;
  fic << "# (13)  : <rho.w> " << finl;
  fic << "# (14)  : <rho.u.u> " << finl;
  fic << "# (15)  : <rho.v.v> " << finl;
  fic << "# (16)  : <rho.w.w> " << finl;
  fic << "# (17)  : <rho.u.v> " << finl;
  fic << "# (18)  : <rho.u.w> " << finl;
  fic << "# (19)  : <rho.v.w> " << finl;
  //fin modif AT 5/06/09
  fic << "# (20)  : <dT/dy> " << finl;		//modif YB 29/06/09
  fic << "# (21)  : <-lambda.dT/dy> " << finl;	//modif YB 30/06/09
  fic << "# (22)  : <-lambda_sm.dT/dy> " << finl;	//modif YB 30/06/09

  fic << "#  " << finl;
  fic << "# Nombre minimum de points par plan Y=cste ayant servi au calcul de la moyenne spatiale : " << min_array(compt_tot) << finl;
  fic << "#  " << finl;

  for (int i=0; i<Y_tot.size(); i++)
    {
      u    = val_moy(i,0,k)/dt;
      v    = val_moy(i,1,k)/dt;
      wl    = val_moy(i,2,k)/dt;
      T    = val_moy(i,13,k)/dt;
      T2   = val_moy(i,14,k)/dt;
      uT   = val_moy(i,15,k)/dt;
      vT   = val_moy(i,16,k)/dt;
      wT   = val_moy(i,17,k)/dt;
      //modified AT 5/06/09
      kappa   = val_moy(i,18,k)/dt;
      rhoT	= val_moy(i,19,k)/dt;
      rhouT	= val_moy(i,20,k)/dt;
      rhovT	= val_moy(i,21,k)/dt;
      rhowT	= val_moy(i,22,k)/dt;
      rhou	= val_moy(i,23,k)/dt;
      rhov	= val_moy(i,24,k)/dt;
      rhow	= val_moy(i,25,k)/dt;
      rhouu	= val_moy(i,26,k)/dt;
      rhovv	= val_moy(i,27,k)/dt;
      rhoww	= val_moy(i,28,k)/dt;
      rhouv	= val_moy(i,29,k)/dt;
      rhouw	= val_moy(i,30,k)/dt;
      rhovw	= val_moy(i,31,k)/dt;
      //fin modif AT 5/06/09
      gradT = val_moy(i,32,k)/dt;		//modif YB 29/06/09
      fluxdiff = val_moy(i,33,k)/dt;	//modif YB 30/06/09
      fluxsm = val_moy(i,34,k)/dt;	//modif YB 30/06/09

      fic << Y_tot(i) << "    "  ;
      // Pour eviter NAN, on prend le dmax(,0):
      fic << T << "    " << sqrt(dmax(T2-T*T,0)) << "    "  ;
      fic << u*T-uT << "    " << v*T-vT << "    " << wl*T-wT << "    ";
      //modified AT 5/06/09
      fic << kappa << "    ";
      fic << rhoT << "    ";
      fic << rhouT << "    ";
      fic << rhovT << "    ";
      fic << rhowT << "    ";
      fic << rhou << "    ";
      fic << rhov << "    ";
      fic << rhow << "    ";
      fic << rhouu << "    ";
      fic << rhovv << "    ";
      fic << rhoww << "    ";
      fic << rhouv << "    ";
      fic << rhouw << "    ";
      fic << rhovw << "    ";
      //fin modif AT 5/06/09
      fic << gradT << "    ";		//modif YB 29/06/09
      fic << fluxdiff << "    ";		//modif YB 30/06/09
      fic << fluxsm << "    ";		//modif YB 30/06/09

      fic << finl;
    }
  fic.flush();
  fic.close();
}



