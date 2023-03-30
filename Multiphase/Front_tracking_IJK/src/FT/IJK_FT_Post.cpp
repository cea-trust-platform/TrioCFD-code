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
// File      : IJK_FT_Post.cpp
// Directory : $IJK_ROOT/src/FT
//
/////////////////////////////////////////////////////////////////////////////

#include <Param.h>
#include <EFichier.h>
#include <SFichier.h>
#include <IJK_FT_Post.h>
#include <IJK_FT.h>
#include <IJK_Lata_writer.h>
#include <IJK_Navier_Stokes_tools.h>
#include <IJK_Splitting.h>
#include <Process.h>    // Process::Journal()
#include <stat_counters.h>

#include <sstream>

//Implemente_liste(IJK_Thermique);
/*
 * Take as main parameter reference to FT to be able to use its members.
 */
IJK_FT_Post::IJK_FT_Post(IJK_FT_double& ijk_ft) :
  statistiques_FT_(ijk_ft), ref_ijk_ft_(ijk_ft), disable_diphasique_(ijk_ft.disable_diphasique_), interfaces_(ijk_ft.interfaces_), pressure_(ijk_ft.pressure_), velocity_(ijk_ft.velocity_),
  d_velocity_(ijk_ft.d_velocity_), splitting_(ijk_ft.splitting_), splitting_ft_(ijk_ft.splitting_ft_), thermique_(ijk_ft.thermique_), energie_(ijk_ft.energie_)
{
  groups_statistiques_FT_.dimensionner(0);
}

void IJK_FT_Post::complete_interpreter(Param& param, Entree& is)
{
  t_debut_statistiques_ = 1.e20;
  check_stats_ = 0;
  expression_pression_analytique_ = "??"; // par defaut, invalide
  fichier_reprise_integrated_velocity_ = "??"; // par defaut, invalide

  fichier_reprise_integrated_pressure_ = "??"; // par defaut, invalide
  fichier_reprise_indicatrice_non_perturbe_ = "??"; // par defaut, invalide
  fichier_reprise_integrated_timescale_ = "??"; // par defaut, invalide
  compteur_post_instantanes_ = 0;
  dt_post_ = 100;
  dt_post_stats_plans_ = 1;
  dt_post_stats_bulles_ = 1;
  //poisson_solver_post_ = xxxx;
  postraiter_sous_pas_de_temps_ = 0;

  param.ajouter_flag("check_stats", &check_stats_);
  param.ajouter("dt_post", &dt_post_);
  param.ajouter("dt_post_stats_plans", &dt_post_stats_plans_);
  param.ajouter("dt_post_stats_bulles", &dt_post_stats_bulles_);
  param.ajouter("champs_a_postraiter", &liste_post_instantanes_);
  param.ajouter_flag("postraiter_sous_pas_de_temps", &postraiter_sous_pas_de_temps_);

  expression_vitesse_analytique_.dimensionner_force(3);
  param.ajouter("expression_vx_ana", &expression_vitesse_analytique_[0]);
  param.ajouter("expression_vy_ana", &expression_vitesse_analytique_[1]);
  param.ajouter("expression_vz_ana", &expression_vitesse_analytique_[2]);

  expression_dvitesse_analytique_.dimensionner_force(3);
  param.ajouter("expression_dvx_ana", &expression_dvitesse_analytique_[0]);
  param.ajouter("expression_dvy_ana", &expression_dvitesse_analytique_[1]);
  param.ajouter("expression_dvz_ana", &expression_dvitesse_analytique_[2]);
  param.ajouter("expression_p_ana", &expression_pression_analytique_);

  expression_gradP_analytique_.dimensionner_force(3);
  param.ajouter("expression_dPdx_ana", &expression_gradP_analytique_[0]);
  param.ajouter("expression_dPdy_ana", &expression_gradP_analytique_[1]);
  param.ajouter("expression_dPdz_ana", &expression_gradP_analytique_[2]);

  expression_gradU_analytique_.dimensionner_force(3);
  param.ajouter("expression_dUdx_ana", &expression_gradU_analytique_[0]);
  param.ajouter("expression_dUdy_ana", &expression_gradU_analytique_[1]);
  param.ajouter("expression_dUdz_ana", &expression_gradU_analytique_[2]);
  expression_gradV_analytique_.dimensionner_force(3);
  param.ajouter("expression_dVdx_ana", &expression_gradV_analytique_[0]);
  param.ajouter("expression_dVdy_ana", &expression_gradV_analytique_[1]);
  param.ajouter("expression_dVdz_ana", &expression_gradV_analytique_[2]);
  expression_gradW_analytique_.dimensionner_force(3);
  param.ajouter("expression_dWdx_ana", &expression_gradW_analytique_[0]);
  param.ajouter("expression_dWdy_ana", &expression_gradW_analytique_[1]);
  param.ajouter("expression_dWdz_ana", &expression_gradW_analytique_[2]);

  // Pour les seconds gradients :
  expression_grad2P_analytique_.dimensionner_force(6);
  param.ajouter("expression_ddPdxdx_ana", &expression_grad2P_analytique_[0]);
  param.ajouter("expression_ddPdydy_ana", &expression_grad2P_analytique_[1]);
  param.ajouter("expression_ddPdzdz_ana", &expression_grad2P_analytique_[2]);
  param.ajouter("expression_ddPdxdy_ana", &expression_grad2P_analytique_[3]);
  param.ajouter("expression_ddPdxdz_ana", &expression_grad2P_analytique_[4]);
  param.ajouter("expression_ddPdydz_ana", &expression_grad2P_analytique_[5]);
  // And for velocities :
  expression_grad2U_analytique_.dimensionner_force(6);
  param.ajouter("expression_ddUdxdx_ana", &expression_grad2U_analytique_[0]);
  param.ajouter("expression_ddUdydy_ana", &expression_grad2U_analytique_[1]);
  param.ajouter("expression_ddUdzdz_ana", &expression_grad2U_analytique_[2]);
  param.ajouter("expression_ddUdxdy_ana", &expression_grad2U_analytique_[3]);
  param.ajouter("expression_ddUdxdz_ana", &expression_grad2U_analytique_[4]);
  param.ajouter("expression_ddUdydz_ana", &expression_grad2U_analytique_[5]);
  //
  expression_grad2V_analytique_.dimensionner_force(6);
  param.ajouter("expression_ddVdxdx_ana", &expression_grad2V_analytique_[0]);
  param.ajouter("expression_ddVdydy_ana", &expression_grad2V_analytique_[1]);
  param.ajouter("expression_ddVdzdz_ana", &expression_grad2V_analytique_[2]);
  param.ajouter("expression_ddVdxdy_ana", &expression_grad2V_analytique_[3]);
  param.ajouter("expression_ddVdxdz_ana", &expression_grad2V_analytique_[4]);
  param.ajouter("expression_ddVdydz_ana", &expression_grad2V_analytique_[5]);
  //
  expression_grad2W_analytique_.dimensionner_force(6);
  param.ajouter("expression_ddWdxdx_ana", &expression_grad2W_analytique_[0]);
  param.ajouter("expression_ddWdydy_ana", &expression_grad2W_analytique_[1]);
  param.ajouter("expression_ddWdzdz_ana", &expression_grad2W_analytique_[2]);
  param.ajouter("expression_ddWdxdy_ana", &expression_grad2W_analytique_[3]);
  param.ajouter("expression_ddWdxdz_ana", &expression_grad2W_analytique_[4]);
  param.ajouter("expression_ddWdydz_ana", &expression_grad2W_analytique_[5]);

  param.ajouter("t_debut_statistiques", &t_debut_statistiques_);

  //param.ajouter("multigrid_solver_post", &poisson_solver_post_);
  // Lecture des sondes :
  param.ajouter("Sondes", &les_sondes_);
}

int IJK_FT_Post::initialise(int reprise)
{
  int nalloc = 0;

  source_spectrale_ = ref_ijk_ft_.forcage_.get_force_ph2();
  //poisson_solver_post_.initialize(splitting_);
  // pour relire les champs de temps integres:
  if (liste_post_instantanes_.contient_("INTEGRATED_TIMESCALE"))
    {
      integrated_timescale_.allocate(splitting_, IJK_Splitting::ELEM, 0);
      nalloc += 1;
      if ((reprise) && (fichier_reprise_integrated_timescale_ != "RESET"))
        {
          if (fichier_reprise_integrated_timescale_ == "??")
            {
              Cerr << "fichier_reprise_integrated_timescale should be specified in the restart file" << finl;
              Process::exit();
            }
          const int timestep_reprise_integrated_timescale_ = 1;
          const Nom& geom_name = integrated_timescale_.get_splitting().get_grid_geometry().le_nom();
          lire_dans_lata(fichier_reprise_integrated_timescale_, timestep_reprise_integrated_timescale_, geom_name, "INTEGRATED_TIMESCALE", integrated_timescale_); // fonction qui lit un champ a partir d'un lata .
        }
      else
        {
          integrated_timescale_.data() = 0.;
          // Question GB pour Antoine : c'est une precaution pour pas qu'il vaille 0 au debut?
          // Mais du coup, on le compte deux fois...
          update_integral_indicatrice(interfaces_.In(), 1. /* Should be the integration timestep */, integrated_timescale_);
        }

    }

  // Pour relire les champs de vitesse et pression integres :
  if (((ref_ijk_ft_.coef_immobilisation_ > 1e-16) && (t_debut_statistiques_ < 1.e10)) || (liste_post_instantanes_.contient_("INTEGRATED_VELOCITY")))
    {
      allocate_velocity(integrated_velocity_, splitting_, 2);
      nalloc += 3;
    }
  if (liste_post_instantanes_.contient_("INTEGRATED_VELOCITY"))
    {
      if ((reprise) && (fichier_reprise_integrated_velocity_ != "RESET"))
        {
          if (fichier_reprise_integrated_velocity_ == "??")
            {
              Cerr << "fichier_reprise_integrated_velocity should be specified in the restart file" << finl;
              Process::exit();
            }
          const int timestep_reprise_integrated_velocity_ = 1;
          Cout << "Lecture vitesse integree initiale dans fichier " << fichier_reprise_integrated_velocity_ << " timestep= " << timestep_reprise_integrated_velocity_ << finl;
          const Nom& geom_name = velocity_[0].get_splitting().get_grid_geometry().le_nom();
          lire_dans_lata(fichier_reprise_integrated_velocity_, timestep_reprise_integrated_velocity_, geom_name, "INTEGRATED_VELOCITY", integrated_velocity_[0], integrated_velocity_[1],
                         integrated_velocity_[2]); // fonction qui lit un champ a partir d'un lata .
        }
      else
        {
          for (int i = 0; i < 3; i++)
            {
              integrated_velocity_[i].data() = 0.;
            }
          velocity_[0].echange_espace_virtuel(velocity_[0].ghost());
          velocity_[1].echange_espace_virtuel(velocity_[1].ghost());
          velocity_[2].echange_espace_virtuel(velocity_[2].ghost());

          update_integral_velocity(velocity_, integrated_velocity_, interfaces_.In(), integrated_timescale_);

        }
    }
  if (((ref_ijk_ft_.coef_immobilisation_ > 1e-16) && (t_debut_statistiques_ < 1.e10)) || (liste_post_instantanes_.contient_("INTEGRATED_PRESSURE")))
    {
      integrated_pressure_.allocate(splitting_, IJK_Splitting::ELEM, 0);
      nalloc += 1;
    }
  if (liste_post_instantanes_.contient_("INTEGRATED_PRESSURE"))
    {
      if ((reprise) && (fichier_reprise_integrated_velocity_ != "RESET"))
        {
          if (fichier_reprise_integrated_pressure_ == "??")
            {
              Cerr << "fichier_reprise_integrated_pressure should be specified in the restart file" << finl;
              Process::exit();
            }
          const int timestep_reprise_integrated_pressure_ = 1;
          Cout << "Lecture pression integree initiale dans fichier " << fichier_reprise_integrated_pressure_ << " timestep= " << timestep_reprise_integrated_pressure_ << finl;
          const Nom& geom_name = pressure_.get_splitting().get_grid_geometry().le_nom();
          lire_dans_lata(fichier_reprise_integrated_pressure_, timestep_reprise_integrated_pressure_, geom_name, "INTEGRATED_PRESSURE", integrated_pressure_); // fonction qui lit un champ a partir d'un lata .

        }
      else
        {
          integrated_pressure_.data() = 0.;

          // Le champ de pression initial ne vaut-il pas forcemment 0?
          update_integral_pressure(pressure_, integrated_pressure_, interfaces_.In(), integrated_timescale_);

        }
    }

  // En reprise, il se peut que le champ ne soit pas dans la liste des posts, mais qu'on l'ait quand meme.
  // Dans ce cas, on choisi de le lire, remplir le field et le re-sauvegarder a la fin (on n'en a rien fait de plus entre temps...)
  if (((ref_ijk_ft_.coef_immobilisation_ > 1e-16) && (t_debut_statistiques_ < 1.e10)) || (liste_post_instantanes_.contient_("INDICATRICE_PERTURBE"))
      || ((reprise) && ((fichier_reprise_indicatrice_non_perturbe_ != "??"))))
    {
      indicatrice_non_perturbe_.allocate(splitting_, IJK_Splitting::ELEM, 0);
      nalloc += 1;
    }

  // Pour le post-traitement de lambda2 :
  if (liste_post_instantanes_.contient_("LAMBDA2"))
    {
      lambda2_.allocate(splitting_, IJK_Splitting::ELEM, 0);
      dudy_.allocate(splitting_, IJK_Splitting::ELEM, 0);
      dvdx_.allocate(splitting_, IJK_Splitting::ELEM, 0);
      dwdy_.allocate(splitting_, IJK_Splitting::ELEM, 0);
      nalloc += 4;
    }

  // Pour le check_stats_ :
  if (check_stats_)
    {
      Cout << "Initialisation champs analytiques (derivee P)" << "\ndPdx = " << expression_gradP_analytique_[0] << "\ndPdy = " << expression_gradP_analytique_[1] << "\ndPdz = "
           << expression_gradP_analytique_[2] << finl;

      Cout << "Initialisation champs analytiques (derivee U)" << "\ndUdx = " << expression_gradU_analytique_[0] << "\ndUdy = " << expression_gradU_analytique_[1] << "\ndUdz = "
           << expression_gradU_analytique_[2] << finl;

      Cout << "Initialisation champs analytiques (derivee V)" << "\ndVdx = " << expression_gradV_analytique_[0] << "\ndVdy = " << expression_gradV_analytique_[1] << "\ndVdz = "
           << expression_gradV_analytique_[2] << finl;

      Cout << "Initialisation champs analytiques (derivee W)" << "\ndWdx = " << expression_gradW_analytique_[0] << "\ndWdy = " << expression_gradW_analytique_[1] << "\ndWdz = "
           << expression_gradW_analytique_[2] << finl;

      Cout << "Initialisation champs analytiques (derivees secondes P) " << "\nddPdxdx = " << expression_grad2P_analytique_[0] << "\nddPdydy = " << expression_grad2P_analytique_[1] << "\nddPdzdz = "
           << expression_grad2P_analytique_[2];
      Cout << "\nddPdxdy = " << expression_grad2P_analytique_[3] << "\nddPdxdz = " << expression_grad2P_analytique_[4] << "\nddPdydz = " << expression_grad2P_analytique_[5] << finl;

      Cout << "Initialisation champs analytiques (derivees secondes U) " << "\nddUdxdx = " << expression_grad2U_analytique_[0] << "\nddUdydy = " << expression_grad2U_analytique_[1] << "\nddUdzdz = "
           << expression_grad2U_analytique_[2];
      Cout << "\nddUdxdy = " << expression_grad2U_analytique_[3] << "\nddUdxdz = " << expression_grad2U_analytique_[4] << "\nddUdydz = " << expression_grad2U_analytique_[5] << finl;

      Cout << "Initialisation champs analytiques (derivees secondes V) " << "\nddVdxdx = " << expression_grad2V_analytique_[0] << "\nddVdydy = " << expression_grad2V_analytique_[1] << "\nddVdzdz = "
           << expression_grad2V_analytique_[2];
      Cout << "\nddVdxdy = " << expression_grad2V_analytique_[3] << "\nddVdxdz = " << expression_grad2V_analytique_[4] << "\nddVdydz = " << expression_grad2V_analytique_[5] << finl;

      Cout << "Initialisation champs analytiques (derivees secondes W) " << "\nddWdxdx = " << expression_grad2W_analytique_[0] << "\nddWdydy = " << expression_grad2W_analytique_[1] << "\nddWdzdz = "
           << expression_grad2W_analytique_[2];
      Cout << "\nddWdxdy = " << expression_grad2W_analytique_[3] << "\nddWdxdz = " << expression_grad2W_analytique_[4] << "\nddWdydz = " << expression_grad2W_analytique_[5] << finl;

      for (int i = 0; i < 3; i++)
        {
          set_field_data(ana_gradP_[i], expression_gradP_analytique_[i]);
          set_field_data(ana_dUd_[i], expression_gradU_analytique_[i]);
          set_field_data(ana_dVd_[i], expression_gradV_analytique_[i]);
          set_field_data(ana_dWd_[i], expression_gradW_analytique_[i]);
          // Pour les deriv secondes :
          set_field_data(ana_grad2Pi_[i], expression_grad2P_analytique_[i]);
          set_field_data(ana_grad2Pc_[i], expression_grad2P_analytique_[3 + i]);
          // And for the 3 components of velocity :
          set_field_data(ana_grad2Ui_[i], expression_grad2U_analytique_[i]);
          set_field_data(ana_grad2Uc_[i], expression_grad2U_analytique_[3 + i]);
          set_field_data(ana_grad2Vi_[i], expression_grad2V_analytique_[i]);
          set_field_data(ana_grad2Vc_[i], expression_grad2V_analytique_[3 + i]);
          set_field_data(ana_grad2Wi_[i], expression_grad2W_analytique_[i]);
          set_field_data(ana_grad2Wc_[i], expression_grad2W_analytique_[3 + i]);

          // Pas necessaire d'echange_espace_virtuel car ghost_ = 0
        }
    }
  return nalloc;
}

void IJK_FT_Post::complete(int reprise)
{
  // Meme if que pour l'allocation.
  // On ne fait le calcul/remplissage du champ que dans un deuxieme temps car on
  // n'avait pas les interfaces avant (lors de l'init)
  if (((ref_ijk_ft_.coef_immobilisation_ > 1e-16) && (t_debut_statistiques_ < 1.e10)) || (liste_post_instantanes_.contient_("INDICATRICE_PERTURBE"))
      || ((reprise) && ((fichier_reprise_indicatrice_non_perturbe_ != "??"))))
    {
      init_indicatrice_non_perturbe();
    }
}

int IJK_FT_Post::initialise_stats(IJK_Splitting& splitting, ArrOfDouble& vol_bulles, const double vol_bulle_monodisperse)
{
  Cout << "Initialisation des statistiques. T_debut_statistiques=" << t_debut_statistiques_ << finl;
  int nalloc = statistiques_FT_.initialize(ref_ijk_ft_, splitting, check_stats_);
  // Si on utilise un seul groupe et qu'on impose un volume unique a toutes les bulles,
  if (vol_bulle_monodisperse >= 0.)
    {
      // on redimensionne le tableau a nb bulles reelles'
      vol_bulles.resize_array(interfaces_.get_nb_bulles_reelles());
      vol_bulles = vol_bulle_monodisperse;
    }
  // S'il n'y a pas qu'un group, on s'occupe des objets stats pour chaque group:
  const int nb_groups = interfaces_.nb_groups();
  if (nb_groups > 1)
    {
      groups_statistiques_FT_.dimensionner(nb_groups);
      for (int igroup = 0; igroup < nb_groups; igroup++)
        {
          groups_statistiques_FT_[igroup].initialize(ref_ijk_ft_, splitting, check_stats_);
        }
    }
  return nalloc;
}

void IJK_FT_Post::init_indicatrice_non_perturbe()
{
  // Est-il deja rempli et stocke?
  // Si on n'est pas en reprise de calcul, le fichier "fichier_reprise_indicatrice_non_perturbe_" est forcement a "??"
  if ((fichier_reprise_indicatrice_non_perturbe_ != "??") && (fichier_reprise_indicatrice_non_perturbe_ != "RESET"))
    {
      const int timestep_reprise_indicatrice_non_perturbe = 1; // 1 ou 0 est le premier? attention au get_db ou latadb...
      Cout << "Lecture indicatrice non perturbee dans fichier " << fichier_reprise_indicatrice_non_perturbe_ << " timestep= " << timestep_reprise_indicatrice_non_perturbe << finl;
      const Nom& geom_name = indicatrice_non_perturbe_.get_splitting().get_grid_geometry().le_nom();
      lire_dans_lata(fichier_reprise_indicatrice_non_perturbe_, timestep_reprise_indicatrice_non_perturbe, geom_name, "INDICATRICE_PERTURBE", indicatrice_non_perturbe_); // fonction qui lit un champ a partir d'un lata .
    }
  else
    {
      // Sinon, on le calcule une fois pour toute (cas bulles fixe = le champ ne varie pas en temps...)
      ArrOfDouble volume_reel;
      DoubleTab position;
      interfaces_.calculer_volume_bulles(volume_reel, position);
      interfaces_.compute_indicatrice_non_perturbe(indicatrice_non_perturbe_, ref_ijk_ft_.itfce().I(), volume_reel, position);
      supprimer_chevauchement(indicatrice_non_perturbe_);
    }
}

// GAB
void IJK_FT_Post::posttraiter_champs_instantanes(const char *lata_name, double current_time, int time_iteration)
{
  statistiques().begin_count(postraitement_counter_);

  const int latastep = compteur_post_instantanes_;
  dumplata_newtime(lata_name, current_time);
  if ((liste_post_instantanes_.contient_("FORCE_PH")) or (liste_post_instantanes_.contient_("CELL_FORCE_PH")))
    {
      source_spectrale_ = ref_ijk_ft_.forcage_.get_force_ph2();
    }
  if (liste_post_instantanes_.contient_("TOUS"))
    {
      liste_post_instantanes_.dimensionner_force(0);
      liste_post_instantanes_.add("VELOCITY");
      liste_post_instantanes_.add("PRESSURE");
      liste_post_instantanes_.add("INDICATRICE_FT");
      liste_post_instantanes_.add("INDICATRICE");
      liste_post_instantanes_.add("RHO");
      liste_post_instantanes_.add("MU");
      liste_post_instantanes_.add("PRESSURE_RHS");
      liste_post_instantanes_.add("VELOCITY_FT");
      liste_post_instantanes_.add("SOURCE_QDM_INTERF");
      liste_post_instantanes_.add("GRAD_INDICATRICE_FT");
      liste_post_instantanes_.add("REBUILT_INDICATRICE_FT");
      liste_post_instantanes_.add("REPULSION_FT"); //("POTENTIEL_FT");
      liste_post_instantanes_.add("AIRE_INTERF");
      liste_post_instantanes_.add("PRESSURE_LIQ");
      liste_post_instantanes_.add("PRESSURE_VAP");
      liste_post_instantanes_.add("GRAD_U");
      liste_post_instantanes_.add("GRAD_V");
      liste_post_instantanes_.add("GRAD_W");
      liste_post_instantanes_.add("SURFACE_VAPEUR_PAR_FACE");
      liste_post_instantanes_.add("BARYCENTRE_VAPEUR_PAR_FACE");
      // GAB, THI sondes
      //liste_post_instantanes_.add("FORCE_PH");

      interfaces_.posttraiter_tous_champs(liste_post_instantanes_);

      {
        int idx_th = 0;
        for (auto&& itr = thermique_.begin(); itr != thermique_.end(); ++itr)
          {
            posttraiter_tous_champs_thermique(liste_post_instantanes_, idx_th);
            ++idx_th;
          }
      }
      {
        int idx_en = 0;
        for (auto&& itr = energie_.begin(); itr != energie_.end(); ++itr)
          {
            posttraiter_tous_champs_energie(liste_post_instantanes_, idx_en);
            ++idx_en;
          }
      }
    }
  int n = liste_post_instantanes_.size();
  if (liste_post_instantanes_.contient_("CURL"))
    {
      // C'est un vecteur mais localise aux elems. On ne peut donc pas le dumper par dumplata_vector :
      n--, dumplata_cellvector(lata_name, "CURL", rot_, latastep);
    }
  if (liste_post_instantanes_.contient_("CRITERE_Q"))
    {
      n--, dumplata_scalar(lata_name, "CRITERE_Q", critere_Q_, latastep);
    }
  if (liste_post_instantanes_.contient_("EXTERNAL_FORCE"))
    {
      n--;
      if (ref_ijk_ft_.coef_immobilisation_ > 1e-16)
        {
          if (!ref_ijk_ft_.interfaces_.get_forcing_method())
            for (int dir = 0; dir < 3; dir++)
              ref_ijk_ft_.redistribute_from_splitting_ft_faces_[dir].redistribute(ref_ijk_ft_.force_rappel_ft_[dir], ref_ijk_ft_.force_rappel_[dir]);
          dumplata_vector(lata_name, "EXTERNAL_FORCE", ref_ijk_ft_.force_rappel_[0], ref_ijk_ft_.force_rappel_[1], ref_ijk_ft_.force_rappel_[2], latastep);
        }
      else
        Cerr << "Posttraitement demande pour EXTERNAL_FORCE but ignored because coef_immobilisation_ <= 1e-16" << finl;
    }
  if (liste_post_instantanes_.contient_("NUM_COMPO"))
    {
      const int ni = num_compo_ft_.ni();
      const int nj = num_compo_ft_.nj();
      const int nk = num_compo_ft_.nk();
      const IntVect& num_compo = interfaces_.get_num_compo();
      for (int k = 0; k < nk; k++)
        for (int j = 0; j < nj; j++)
          for (int i = 0; i < ni; i++)
            {
              const int num_elem = splitting_ft_.convert_ijk_cell_to_packed(i, j, k);
              num_compo_ft_(i, j, k) = num_compo[num_elem];
            }
      n--, dumplata_scalar(lata_name, "NUM_COMPO", num_compo_ft_, latastep);
    }
  if (liste_post_instantanes_.contient_("VELOCITY"))
    n--, dumplata_vector(lata_name, "VELOCITY", velocity_[0], velocity_[1], velocity_[2], latastep);
  // GAB
  if (liste_post_instantanes_.contient_("FORCE_PH"))
    {
      n--, dumplata_vector(lata_name, "FORCE_PH", source_spectrale_[0], source_spectrale_[1], source_spectrale_[2], latastep);
    }
  //
  if (liste_post_instantanes_.contient_("INTEGRATED_VELOCITY"))
    {
      update_integral_velocity(velocity_, integrated_velocity_, interfaces_.In(), integrated_timescale_);
      n--, dumplata_vector(lata_name, "INTEGRATED_VELOCITY", integrated_velocity_[0], integrated_velocity_[1], integrated_velocity_[2], latastep);
    }
  if (liste_post_instantanes_.contient_("INTEGRATED_PRESSURE"))
    {
      //      integrated_pressure_ += timestep_*pressure_;
      update_integral_pressure(pressure_, integrated_pressure_, interfaces_.In(), integrated_timescale_);
      n--, dumplata_scalar(lata_name, "INTEGRATED_PRESSURE", integrated_pressure_, latastep);
    }
  if (liste_post_instantanes_.contient_("INDICATRICE_PERTURBE"))
    {
      //Faut-il faire un update_...( indicatrice_non_perturbe_) avant?
      n--, dumplata_scalar(lata_name, "INDICATRICE_PERTURBE", indicatrice_non_perturbe_, latastep);
    }
  if (liste_post_instantanes_.contient_("INTEGRATED_TIMESCALE"))
    {
      update_integral_indicatrice(interfaces_.In(), 1. /* Should be the integration timestep */, integrated_timescale_);
      n--, dumplata_scalar(lata_name, "INTEGRATED_TIMESCALE", integrated_timescale_, latastep);
    }

  if (liste_post_instantanes_.contient_("COORDS"))
    n--, dumplata_vector(lata_name, "COORDS", coords_[0], coords_[1], coords_[2], latastep);
  if (liste_post_instantanes_.contient_("LAMBDA2"))
    {
      get_update_lambda2();
      n--, dumplata_scalar(lata_name, "LAMBDA2", lambda2_, latastep);
    }
  if (liste_post_instantanes_.contient_("VELOCITY_ANA"))
    {
      for (int i = 0; i < 3; i++)
        set_field_data(velocity_ana_[i], expression_vitesse_analytique_[i], current_time);
      n--, dumplata_vector(lata_name, "VELOCITY_ANA", velocity_ana_[0], velocity_ana_[1], velocity_ana_[2], latastep);
    }
  if (liste_post_instantanes_.contient_("VARIABLE_SOURCE"))
    {
      n--, dumplata_vector(lata_name, "VARIABLE_SOURCE", ref_ijk_ft_.variable_source_[0], ref_ijk_ft_.variable_source_[1], ref_ijk_ft_.variable_source_[2], latastep);
    }

  if (liste_post_instantanes_.contient_("ECART_ANA"))
    {
      // double err[3] = {0., 0., 0.};
      Cerr << "GB: ERROR FIELD " << current_time;
      for (int dir = 0; dir < 3; dir++)
        {
          double err = 0.;
          set_field_data(velocity_ana_[dir], expression_vitesse_analytique_[dir], current_time);
          const int ni = velocity_[dir].ni();
          const int nj = velocity_[dir].nj();
          const int nk = velocity_[dir].nk();
          const int ntot = Process::mp_sum(ni * nj * nk);
          for (int k = 0; k < nk; k++)
            for (int j = 0; j < nj; j++)
              for (int i = 0; i < ni; i++)
                {
                  const double val = velocity_ana_[dir](i, j, k) - velocity_[dir](i, j, k);
                  ecart_ana_[dir](i, j, k) = val;
                  err += val * val;
                }
          err = Process::mp_sum(err);
          err = sqrt(err / ntot);
          Cerr << " " << err;
          if (!Process::je_suis_maitre())
            {
              Process::Journal() << "IJK_FT_Post::posttraiter_champs_instantanes : Champ ECART_ANA sur ce proc (ni,nj,nk,ntot):" << " " << ni << " " << nj << " " << nk << " " << ntot << finl;
            }
        }
      Cerr << finl;
      n--, dumplata_vector(lata_name, "ECART_ANA", ecart_ana_[0], ecart_ana_[1], ecart_ana_[2], latastep);
    }
  if (liste_post_instantanes_.contient_("PRESSURE_ANA"))
    {
      set_field_data(pressure_ana_, expression_pression_analytique_, current_time);
      n--, dumplata_scalar(lata_name, "PRESSURE_ANA", pressure_ana_, latastep);
    }
  if (liste_post_instantanes_.contient_("ECART_P_ANA"))
    {
      double ct = current_time;
      if (ref_ijk_ft_.get_time_scheme() == IJK_FT_double::EULER_EXPLICITE)
        {
          ct -= ref_ijk_ft_.timestep_;
        }
      else if (ref_ijk_ft_.get_time_scheme() == IJK_FT_double::RK3_FT)
        {
          Cerr << "rkstep " << ref_ijk_ft_.rk_step_ << finl;
          int rk_step_before = ref_ijk_ft_.rk_step_;
          if ((rk_step_before == 0) || (rk_step_before == 3))
            rk_step_before = 2;
          else if (rk_step_before == 1)
            rk_step_before = 0;
          else
            /* ici, c'est rk_step_before=2 */
            rk_step_before = 1;
          Cerr << "rkstep_before " << rk_step_before << finl;
          const double intermediate_dt = compute_fractionnal_timestep_rk3(ref_ijk_ft_.timestep_, rk_step_before);
          ct -= intermediate_dt;
        }
      else
        {
          Cerr << "To do for other time scheme" << finl;
        }
      Cerr << "GB: ERROR P FIELD " << ct;
      double err = 0.;
      set_field_data(pressure_ana_, expression_pression_analytique_, ct);
      const int ni = pressure_.ni();
      const int nj = pressure_.nj();
      const int nk = pressure_.nk();
      const int ntot = Process::mp_sum(ni * nj * nk);
      // La pression est definie a une constante pres:
      const double cst_press = pressure_ana_(0, 0, 0) - pressure_(0, 0, 0);
      for (int k = 0; k < nk; k++)
        for (int j = 0; j < nj; j++)
          for (int i = 0; i < ni; i++)
            {
              const double val = pressure_ana_(i, j, k) - pressure_(i, j, k) - cst_press;
              ecart_p_ana_(i, j, k) = val;
              err += val * val;
            }
      err = Process::mp_sum(err);
      err = sqrt(err / ntot);
      Cerr << " " << err;
      if (!Process::je_suis_maitre())
        {
          Process::Journal() << "IJK_FT_Post::posttraiter_champs_instantanes : Champ ECART_P_ANA sur ce proc (ni,nj,nk,ntot):" << " " << ni << " " << nj << " " << nk << " " << ntot << finl;
        }
      ecart_p_ana_.echange_espace_virtuel(ecart_p_ana_.ghost());
      Cerr << finl;
      n--, dumplata_scalar(lata_name, "ECART_P_ANA", ecart_p_ana_, latastep);
    }
  if (liste_post_instantanes_.contient_("D_VELOCITY_ANA"))
    {
      {
        // double err[3] = {0., 0., 0.};
        Cerr << "GB: ERROR DV FIELD " << current_time;
        for (int dir = 0; dir < 3; dir++)
          {
            double err = 0.;
            set_field_data(d_velocity_ana_[dir], expression_dvitesse_analytique_[dir], current_time);
            const int ni = d_velocity_[dir].ni();
            const int nj = d_velocity_[dir].nj();
            const int nk = d_velocity_[dir].nk();
            const int ntot = Process::mp_sum(ni * nj * nk);
            for (int k = 0; k < nk; k++)
              for (int j = 0; j < nj; j++)
                for (int i = 0; i < ni; i++)
                  {
                    const double val = d_velocity_ana_[dir](i, j, k) - d_velocity_[dir](i, j, k);
                    err += val * val;
                  }
            err = Process::mp_sum(err);
            err = sqrt(err / ntot);
            Cerr << " " << err;
            if (!Process::je_suis_maitre())
              {
                Process::Journal() << "IJK_FT_Post::posttraiter_champs_instantanes : Champ ECART_ANA sur ce proc (ni,nj,nk,ntot):" << " " << ni << " " << nj << " " << nk << " " << ntot << finl;
              }
          }
        Cerr << finl;
      }
      n--, dumplata_vector(lata_name, "D_VELOCITY_ANA", d_velocity_ana_[0], d_velocity_ana_[1], d_velocity_ana_[2], latastep);
    }
  if (liste_post_instantanes_.contient_("D_VELOCITY"))
    {
//      d_velocity_[0].echange_espace_virtuel(2 /*, IJK_Field_ST::EXCHANGE_GET_AT_RIGHT_I*/);
//      d_velocity_[1].echange_espace_virtuel(2 /*, IJK_Field_ST::EXCHANGE_GET_AT_RIGHT_J*/);
//      d_velocity_[2].echange_espace_virtuel(2 /*, IJK_Field_ST::EXCHANGE_GET_AT_RIGHT_K*/);
      n--, dumplata_vector(lata_name, "D_VELOCITY", d_velocity_[0], d_velocity_[1], d_velocity_[2], latastep);
    }
  if (liste_post_instantanes_.contient_("OP_CONV"))
    {
      n--, dumplata_vector(lata_name, "DU_DT", op_conv_[0], op_conv_[1], op_conv_[2], latastep);
    }

  if ((liste_post_instantanes_.contient_("GRAD_P")) or (liste_post_instantanes_.contient_("CELL_GRAD_P")))
    n--, dumplata_vector(lata_name, "dPd", grad_P_[0], grad_P_[1], grad_P_[2], latastep);
  // Pour le check_stats_ :
  if (check_stats_)
    {

      //MR
      /*   if (liste_post_instantanes_.contient_("GRAD_T"))
       n--,dumplata_vector(lata_name,"dTd", grad_T_[0], grad_T_[1], grad_T_[2], latastep);*/

      if (liste_post_instantanes_.contient_("ANA_GRAD_P"))
        n--, dumplata_vector(lata_name, "ANA_dPd", ana_gradP_[0], ana_gradP_[1], ana_gradP_[2], latastep);
      if (liste_post_instantanes_.contient_("ANA_GRAD_U"))
        {
          n--, dumplata_cellvector(lata_name, "ANA_dUd", ana_dUd_, latastep);
        }
      if (liste_post_instantanes_.contient_("ANA_GRAD_V"))
        {
          n--, dumplata_cellvector(lata_name, "ANA_dVd", ana_dVd_, latastep);
        }
      if (liste_post_instantanes_.contient_("ANA_GRAD_W"))
        {
          n--, dumplata_cellvector(lata_name, "ANA_dWd", ana_dWd_, latastep);
        }
      if (liste_post_instantanes_.contient_("ANA_GRAD2_P"))
        {
          n--, dumplata_cellvector(lata_name, "ANA_ddPdd", ana_grad2Pi_, latastep);
          dumplata_scalar(lata_name, "ANA_ddPdxdy", ana_grad2Pc_[0], latastep);
          dumplata_scalar(lata_name, "ANA_ddPdxdz", ana_grad2Pc_[1], latastep);
          dumplata_scalar(lata_name, "ANA_ddPdydz", ana_grad2Pc_[2], latastep);
        }
      if (liste_post_instantanes_.contient_("ANA_GRAD2_U"))
        {
          n--, dumplata_cellvector(lata_name, "ANA_ddUdd", ana_grad2Ui_, latastep);
          dumplata_scalar(lata_name, "ANA_ddUdxdy", ana_grad2Uc_[0], latastep);
          dumplata_scalar(lata_name, "ANA_ddUdxdz", ana_grad2Uc_[1], latastep);
          dumplata_scalar(lata_name, "ANA_ddUdydz", ana_grad2Uc_[2], latastep);
        }
      if (liste_post_instantanes_.contient_("ANA_GRAD2_V"))
        {
          n--, dumplata_cellvector(lata_name, "ANA_ddVdd", ana_grad2Vi_, latastep);
          dumplata_scalar(lata_name, "ANA_ddVdxdy", ana_grad2Vc_[0], latastep);
          dumplata_scalar(lata_name, "ANA_ddVdxdz", ana_grad2Vc_[1], latastep);
          dumplata_scalar(lata_name, "ANA_ddVdydz", ana_grad2Vc_[2], latastep);
        }
      if (liste_post_instantanes_.contient_("ANA_GRAD2_W"))
        {
          n--, dumplata_cellvector(lata_name, "ANA_ddWdd", ana_grad2Wi_, latastep);
          dumplata_scalar(lata_name, "ANA_ddWdxdy", ana_grad2Wc_[0], latastep);
          dumplata_scalar(lata_name, "ANA_ddWdxdz", ana_grad2Wc_[1], latastep);
          dumplata_scalar(lata_name, "ANA_ddWdydz", ana_grad2Wc_[2], latastep);
        }
      if (liste_post_instantanes_.contient_("GRAD2_P"))
        {
          const FixedVector<IJK_Field_double, 3>& grad2Pi = statistiques_FT_.get_IJK_vector_field("grad2Pi");
          const FixedVector<IJK_Field_double, 3>& grad2Pc = statistiques_FT_.get_IJK_vector_field("grad2Pc");
          n--, dumplata_cellvector(lata_name, "ddPdd", grad2Pi, latastep);
          dumplata_scalar(lata_name, "ddPdxdy", grad2Pc[0], latastep);
          dumplata_scalar(lata_name, "ddPdxdz", grad2Pc[1], latastep);
          dumplata_scalar(lata_name, "ddPdydz", grad2Pc[2], latastep);
        }
      if (liste_post_instantanes_.contient_("GRAD2_U"))
        {
          const FixedVector<IJK_Field_double, 3>& grad2Ui = statistiques_FT_.get_IJK_vector_field("grad2Ui");
          const FixedVector<IJK_Field_double, 3>& grad2Uc = statistiques_FT_.get_IJK_vector_field("grad2Uc");
          n--, dumplata_cellvector(lata_name, "ddUdd", grad2Ui, latastep);
          dumplata_scalar(lata_name, "ddUdxdy", grad2Uc[0], latastep);
          dumplata_scalar(lata_name, "ddUdxdz", grad2Uc[1], latastep);
          dumplata_scalar(lata_name, "ddUdydz", grad2Uc[2], latastep);
        }
      if (liste_post_instantanes_.contient_("GRAD2_V"))
        {
          const FixedVector<IJK_Field_double, 3>& grad2Vi = statistiques_FT_.get_IJK_vector_field("grad2Vi");
          const FixedVector<IJK_Field_double, 3>& grad2Vc = statistiques_FT_.get_IJK_vector_field("grad2Vc");
          n--, dumplata_cellvector(lata_name, "ddVdd", grad2Vi, latastep);
          dumplata_scalar(lata_name, "ddVdxdy", grad2Vc[0], latastep);
          dumplata_scalar(lata_name, "ddVdxdz", grad2Vc[1], latastep);
          dumplata_scalar(lata_name, "ddVdydz", grad2Vc[2], latastep);
        }
      if (liste_post_instantanes_.contient_("GRAD2_W"))
        {
          const FixedVector<IJK_Field_double, 3>& grad2Wi = statistiques_FT_.get_IJK_vector_field("grad2Wi");
          const FixedVector<IJK_Field_double, 3>& grad2Wc = statistiques_FT_.get_IJK_vector_field("grad2Wc");
          n--, dumplata_cellvector(lata_name, "ddWdd", grad2Wi, latastep);
          dumplata_scalar(lata_name, "ddWdxdy", grad2Wc[0], latastep);
          dumplata_scalar(lata_name, "ddWdxdz", grad2Wc[1], latastep);
          dumplata_scalar(lata_name, "ddWdydz", grad2Wc[2], latastep);
        }
    }
  if (liste_post_instantanes_.contient_("CELL_VELOCITY"))
    {
      const int kmax = cell_velocity_[0].nk();
      const int jmax = cell_velocity_[0].nj();
      const int imax = cell_velocity_[0].ni();
      for (int k = 0; k < kmax; k++)
        for (int j = 0; j < jmax; j++)
          for (int i = 0; i < imax; i++)
            {
              double u = (velocity_[0](i, j, k) + velocity_[0](i + 1, j, k)) * 0.5;
              double v = (velocity_[1](i, j, k) + velocity_[1](i, j + 1, k)) * 0.5;
              double w = (velocity_[2](i, j, k) + velocity_[2](i, j, k + 1)) * 0.5;
              cell_velocity_[0](i, j, k) = u;
              cell_velocity_[1](i, j, k) = v;
              cell_velocity_[2](i, j, k) = w;
            }
      n--, dumplata_cellvector(lata_name, "VELOCITY" /* AT CELL-CENTER */, cell_velocity_, latastep);
    }
  // GAB
  if (liste_post_instantanes_.contient_("CELL_FORCE_PH"))
    {
      const int kmax = cell_source_spectrale_[0].nk();
      const int jmax = cell_source_spectrale_[0].nj();
      const int imax = cell_source_spectrale_[0].ni();
      for (int k = 0; k < kmax; k++)
        for (int j = 0; j < jmax; j++)
          for (int i = 0; i < imax; i++)
            {
              double f;
              double g;
              double h;

              {
                // Interpolation d'ordre 1 des vitesses aux centres des elements
                f = (source_spectrale_[0](i, j, k) + source_spectrale_[0](i + 1, j, k)) * 0.5;
                g = (source_spectrale_[1](i, j, k) + source_spectrale_[1](i, j + 1, k)) * 0.5;
                h = (source_spectrale_[2](i, j, k) + source_spectrale_[2](i, j, k + 1)) * 0.5;

              }

              cell_source_spectrale_[0](i, j, k) = f;
              cell_source_spectrale_[1](i, j, k) = g;
              cell_source_spectrale_[2](i, j, k) = h;
            }
      n--, dumplata_cellvector(lata_name, "FORCE_PH" /* AT CELL-CENTER */, cell_source_spectrale_, latastep);
    }
  //
  // GAB
  if (liste_post_instantanes_.contient_("CELL_GRAD_P"))
    {
      const int kmax = cell_grad_p_[0].nk();
      const int jmax = cell_grad_p_[0].nj();
      const int imax = cell_grad_p_[0].ni();
      for (int k = 0; k < kmax; k++)
        for (int j = 0; j < jmax; j++)
          for (int i = 0; i < imax; i++)
            {
              double p;
              double q;
              double r;

              {
                // Interpolation d'ordre 1 des vitesses aux centres des elements
                p = (grad_P_[0](i, j, k) + grad_P_[0](i + 1, j, k)) * 0.5;
                q = (grad_P_[1](i, j, k) + grad_P_[1](i, j + 1, k)) * 0.5;
                r = (grad_P_[2](i, j, k) + grad_P_[2](i, j, k + 1)) * 0.5;

              }

              cell_grad_p_[0](i, j, k) = p;
              cell_grad_p_[1](i, j, k) = q;
              cell_grad_p_[2](i, j, k) = r;
            }
      n--, dumplata_cellvector(lata_name, "GRAD_P" /* AT CELL-CENTER */, cell_grad_p_, latastep);
    }
  //
  if (liste_post_instantanes_.contient_("GRAD_U"))
    {
      const FixedVector<IJK_Field_double, 3>& gradU = statistiques_FT_.get_IJK_vector_field("gradU");
      n--, dumplata_cellvector(lata_name, "dUd", gradU, latastep);
    }
  if (liste_post_instantanes_.contient_("GRAD_V"))
    {
      const FixedVector<IJK_Field_double, 3>& gradV = statistiques_FT_.get_IJK_vector_field("gradV");
      n--, dumplata_cellvector(lata_name, "dVd", gradV, latastep);
    }
  if (liste_post_instantanes_.contient_("GRAD_W"))
    {
      const FixedVector<IJK_Field_double, 3>& gradW = statistiques_FT_.get_IJK_vector_field("gradW");
      n--, dumplata_cellvector(lata_name, "dWd", gradW, latastep);
    }

  if (liste_post_instantanes_.contient_("PRESSURE"))
    n--, dumplata_scalar(lata_name, "PRESSURE", pressure_, latastep);
  if (liste_post_instantanes_.contient_("D_PRESSURE"))
    n--, dumplata_scalar(lata_name, "D_PRESSURE", ref_ijk_ft_.d_pressure_, latastep);
  if (liste_post_instantanes_.contient_("INDICATRICE"))
    n--, dumplata_scalar(lata_name, "INDICATRICE", interfaces_.In(), latastep);
  if (liste_post_instantanes_.contient_("INDICATRICE_FT"))
    n--, dumplata_scalar(lata_name, "INDICATRICE_FT", interfaces_.In_ft(), latastep);
  if (liste_post_instantanes_.contient_("MU"))
    n--, dumplata_scalar(lata_name, "MU", ref_ijk_ft_.molecular_mu_, latastep);
  if (liste_post_instantanes_.contient_("RHO"))
    n--, dumplata_scalar(lata_name, "RHO", ref_ijk_ft_.rho_field_, latastep);
  if (liste_post_instantanes_.contient_("PRESSURE_RHS"))
    n--, dumplata_scalar(lata_name, "PRESSURE_RHS", ref_ijk_ft_.pressure_rhs_, latastep);
  if (liste_post_instantanes_.contient_("VELOCITY_FT"))
    n--, dumplata_vector(lata_name, "VELOCITY_FT", ref_ijk_ft_.velocity_ft_[0], ref_ijk_ft_.velocity_ft_[1], ref_ijk_ft_.velocity_ft_[2], latastep);
  if (liste_post_instantanes_.contient_("SOURCE_QDM_INTERF"))
    n--, dumplata_vector(lata_name, "SOURCE_QDM_INTERF", ref_ijk_ft_.terme_source_interfaces_ft_[0], ref_ijk_ft_.terme_source_interfaces_ft_[1], ref_ijk_ft_.terme_source_interfaces_ft_[2], latastep);
  if (liste_post_instantanes_.contient_("CELL_SOURCE_QDM_INTERF"))
    n--, dumplata_cellvector(lata_name, "CELL_SOURCE_QDM_INTERF", ref_ijk_ft_.terme_source_interfaces_ns_, latastep);
  if (liste_post_instantanes_.contient_("GRAD_INDICATRICE_FT"))
    n--, dumplata_vector(lata_name, "GRAD_INDICATRICE_FT", grad_I_ft_[0], grad_I_ft_[1], grad_I_ft_[2], latastep);
  if (liste_post_instantanes_.contient_("REBUILT_INDICATRICE_FT"))
    n--, dumplata_scalar(lata_name, "REBUILT_INDICATRICE_FT", rebuilt_indic_, latastep);
  //  if (liste_post_instantanes_.contient_("POTENTIEL_FT"))
  //    n--,dumplata_scalar(lata_name,"POTENTIEL_FT", potentiel_, latastep);
  if (liste_post_instantanes_.contient_("REPULSION_FT"))
    n--, dumplata_scalar(lata_name, "REPULSION_FT", potentiel_, latastep);
  if (liste_post_instantanes_.contient_("AIRE_INTERF"))
    {
      interfaces_.calculer_aire_interfaciale(ai_ft_);
      n--, dumplata_scalar(lata_name, "AIRE_INTERF", ai_ft_, latastep);
    }
  if (liste_post_instantanes_.contient_("COURBURE_AIRE_INTERF"))
    {
      // On suppose implicitement qu'il est bien calcule et mis a jour avant d'arriver ici.
      n--, dumplata_scalar(lata_name, "COURBURE_AIRE_INTERF", kappa_ai_ft_, latastep);
    }
  if (liste_post_instantanes_.contient_("NORMALE_EULER"))
    {
      // On suppose implicitement qu'il est bien calcule et mis a jour avant d'arriver ici.
      n--, dumplata_cellvector(lata_name, "NORMALE_EULER", normale_cell_ft_, latastep);
    }

  if (liste_post_instantanes_.contient_("PRESSURE_LIQ"))
    {
      // On suppose implicitement qu'il est bien calcule et mis a jour avant d'arriver ici.
      n--, dumplata_scalar(lata_name, "PRESSURE_LIQ", extended_pl_, latastep);
    }
  if (liste_post_instantanes_.contient_("PRESSURE_VAP"))
    {
      // On suppose implicitement qu'il est bien calcule et mis a jour avant d'arriver ici.
      n--, dumplata_scalar(lata_name, "PRESSURE_VAP", extended_pv_, latastep);
    }
  if (liste_post_instantanes_.contient_("GROUPS"))
    n--, dumplata_cellvector(lata_name, "GROUPS", interfaces_.groups_indicatrice_n_ns(), latastep);
  if (liste_post_instantanes_.contient_("GROUPS_FT"))
    n--, dumplata_cellvector(lata_name, "GROUPS_FT", interfaces_.groups_indicatrice_n_ft(), latastep);
  if (liste_post_instantanes_.contient_("SURFACE_VAPEUR_PAR_FACE"))
    {
      Cerr << "Tentative de sauvegarder champ surface vapeur par face" << finl;
      n--, dumplata_vector(lata_name, "SURFACE_VAPEUR_PAR_FACE", interfaces_.get_surface_vapeur_par_face()[0], interfaces_.get_surface_vapeur_par_face()[1],
                           interfaces_.get_surface_vapeur_par_face()[2], latastep);
      Cerr << "Reussi" << finl;
    }
  if (liste_post_instantanes_.contient_("BARYCENTRE_VAPEUR_PAR_FACE"))
    {
      n--, dumplata_vector(lata_name, "BARYCENTRE_X_VAPEUR_PAR_FACE", interfaces_.get_barycentre_vapeur_par_face()[0][0], interfaces_.get_barycentre_vapeur_par_face()[0][1],
                           interfaces_.get_barycentre_vapeur_par_face()[0][2], latastep);
      n--, dumplata_vector(lata_name, "BARYCENTRE_Y_VAPEUR_PAR_FACE", interfaces_.get_barycentre_vapeur_par_face()[1][0], interfaces_.get_barycentre_vapeur_par_face()[1][1],
                           interfaces_.get_barycentre_vapeur_par_face()[1][2], latastep);
      n--, dumplata_vector(lata_name, "BARYCENTRE_Z_VAPEUR_PAR_FACE", interfaces_.get_barycentre_vapeur_par_face()[2][0], interfaces_.get_barycentre_vapeur_par_face()[2][1],
                           interfaces_.get_barycentre_vapeur_par_face()[2][2], latastep);
    }

  if (!disable_diphasique_)
    {
      interfaces_.update_surface_normale(); // necessaire avant posttraiter_champs_instantanes_thermique_interfaciaux
      n -= interfaces_.posttraiter_champs_instantanes(liste_post_instantanes_, lata_name, latastep);
    }

  {
    int idx_th = 0;
    for (auto &itr : thermique_)
      {
        int nb = posttraiter_champs_instantanes_thermique(liste_post_instantanes_, lata_name, latastep, current_time, itr, idx_th);
        // Interfacial thermal fields :
        if (!disable_diphasique_)
          nb += posttraiter_champs_instantanes_thermique_interfaciaux(liste_post_instantanes_, lata_name, latastep, current_time, itr, idx_th);

        if (idx_th == 0)
          n -= nb; // On compte comme "un" tous les CHAMPS_N (ou N est la longueur de la liste)
        ++idx_th;
      }
    // TODO: finir post-traitement de l'energie, choisir a quel niveau le faire.
    int idx_en = 0;
    for (auto &itr : energie_)
      {
        int nb = posttraiter_champs_instantanes_energie(liste_post_instantanes_, lata_name, latastep, current_time, itr, idx_en);
        // Interfacial thermal fields :
        if (!disable_diphasique_)
          nb += posttraiter_champs_instantanes_energie_interfaciaux(liste_post_instantanes_, lata_name, latastep, current_time, itr, idx_en);

        if (idx_en == 0)
          n -= nb; // On compte comme "un" tous les CHAMPS_N (ou N est la longueur de la liste)
        ++idx_en;
      }

    Cerr << "les champs postraites sont: " << liste_post_instantanes_ << finl;
  }

  if (n > 0)
    {
      Cerr << "Il y a des noms de champs a postraiter inconnus dans la liste de champs a postraiter" << finl;
      Process::exit();
    }

  // GB : Je crois que ce n'est plus fait ici car les frequences sont differentes...
  //ecrire_statistiques_bulles(compteur_post_instantanes_ == 0);

  compteur_post_instantanes_++;
  statistiques().end_count(postraitement_counter_);
}

// 2020.03.12. CHOIX : Meme en disable_diphasique, on fait appel a la classe fille stats FT
void IJK_FT_Post::posttraiter_statistiques_plans(double current_time)
{
  statistiques().begin_count(postraitement_counter_);

  if (Process::je_suis_maitre())
    {
      Nom n("");
      // post-traitement des stats diphasiques :
      // (si calcul monophasique, on devrait avoir chi=1)
      for (int flag_valeur_instantanee = 0; flag_valeur_instantanee < 2; flag_valeur_instantanee++)
        {
          if (disable_diphasique_)
            n = "monophasique_";
          else
            n = "diphasique_";

          if (flag_valeur_instantanee == 0)
            n += "statistiques_";
          else
            n += "moyenne_spatiale_";

          n += Nom(current_time);
          if ((flag_valeur_instantanee) || (statistiques_FT_.t_integration() > 0.))
            {
              SFichier f(n + Nom(".txt"));
              f.setf(ios::scientific);      // precision pour allez chercher les 4 ordres
              f.precision(15);
              statistiques_FT_.postraiter(f, flag_valeur_instantanee /* flag pour ecrire la moyenne instantanee ou la moyenne */);
              statistiques_FT_.postraiter_thermique(current_time); /* moyenne instantanee et temporelle */

              // S'il n'y a pas qu'un group, on posttraite les objets stats pour chaque group:
              if ((!disable_diphasique_) && (interfaces_.nb_groups() > 1))
                {
                  for (int igroup = 0; igroup < interfaces_.nb_groups(); igroup++)
                    {
                      SFichier figroup(n + Nom("_grp") + Nom(igroup) + Nom(".txt"));
                      figroup.setf(ios::scientific);
                      figroup.precision(15);
                      groups_statistiques_FT_[igroup].postraiter(figroup, flag_valeur_instantanee /* flag pour ecrire la moyenne instantanee ou la moyenne */);
                      if (flag_valeur_instantanee == 1)
                        groups_statistiques_FT_[igroup].postraiter_thermique(current_time); /* moyenne instantanee et temporelle */
                    }
                }
            }
        }
      statistiques_FT_.postraiter_thermique(current_time); /* moyenne instantanee et temporelle */
    }
  statistiques().end_count(postraitement_counter_);

}

// Le nom du fichier est base sur le nom du cas...
// Si reset!=0, on efface le fichier avant d'ecrire, sinon on ajoute...
void IJK_FT_Post::ecrire_statistiques_bulles(int reset, const Nom& nom_cas, const ArrOfDouble& gravite, const double current_time) const
{
  if (disable_diphasique_)
    return;

  statistiques().begin_count(postraitement_counter_);

  ArrOfDouble volume;
  DoubleTab position;
  ArrOfDouble surface;
  const int nbulles = interfaces_.get_nb_bulles_reelles();
  DoubleTab hauteurs_bulles(nbulles, 3);
  DoubleTab bounding_box;
  interfaces_.calculer_bounding_box_bulles(bounding_box);
  for (int ib = 0; ib < nbulles; ib++)
    for (int dir = 0; dir < 3; dir++)
      hauteurs_bulles(ib, dir) = bounding_box(ib, dir, 1) - bounding_box(ib, dir, 0);

  // La methode calcule a present les surfaces meme pour les bulles ghost.
  // Pour les enlever, il suffit simplement de reduire la taille du tableau :
  interfaces_.calculer_surface_bulles(surface);
  surface.resize_array(nbulles);

  // La methode calcule a present les volumes meme pour les bulles ghost.
  // Pour les enlever, il suffit simplement de reduire la taille du tableau :
  interfaces_.calculer_volume_bulles(volume, position);
  volume.resize_array(nbulles);
  position.resize(nbulles, 3);

  DoubleTab poussee;
  interfaces_.calculer_poussee_bulles(gravite, poussee);
  if (1)
    {
      int idx_th = 0;
      for (auto &itr : thermique_)
        {
          ArrOfDouble interfacial_temperature;
          ArrOfDouble interfacial_phin_ai;
          // To transfer the field to FT splitting (because interfaces are there...) !!! NEEDED for compute_interfacial_temperature
          IJK_Field_double& temperature_ft = itr.get_temperature_ft();
          ref_ijk_ft_.redistribute_to_splitting_ft_elem_.redistribute(itr.get_temperature(), temperature_ft);
          temperature_ft.echange_espace_virtuel(temperature_ft.ghost());
          //itr.compute_interfacial_temperature(interfacial_temperature, interfacial_phin_ai, itr.get_storage());
          itr.compute_interfacial_temperature2(interfacial_temperature, interfacial_phin_ai);

          // Compute Bubble mean :
          ArrOfDouble Ti_per_bubble;
          ArrOfDouble phin_per_bubble;
          interfaces_.compute_surface_average_per_bubble(surface, interfacial_phin_ai, phin_per_bubble);
          interfaces_.compute_surface_average_per_bubble(surface, interfacial_temperature, Ti_per_bubble);
          if (Process::je_suis_maitre())
            {
              char s[1000];
              const char *nomcas = nom_cas;
              SFichier fic;
              const int n = Ti_per_bubble.size_array();
              IOS_OPEN_MODE mode = (reset) ? ios::out : ios::app;

#ifndef INT_is_64_
              snprintf(s, 1000, "%s_bulles_Ti_%d.out", nomcas, idx_th);
#else
              snprintf(s, 1000, "%s_bulles_Ti_%ld.out", nomcas, idx_th);
#endif
              // Cerr << "Ecriture des donnees par bulles: fichier " << s << finl;
              fic.ouvrir(s, mode);
              snprintf(s, 1000, "%.16e ", current_time);
              fic << s;
              for (int i = 0; i < n; i++)
                {
                  snprintf(s, 1000, "%.16e ", Ti_per_bubble[i]);
                  fic << s;
                }
              fic << finl;
              fic.close();

              // Cerr << "Ecriture des donnees par bulles: fichier " << s << finl;
#ifndef INT_is_64_
              snprintf(s, 1000, "%s_bulles_phin_%d.out", nomcas, idx_th);
#else
              snprintf(s, 1000, "%s_bulles_phin_%ld.out", nomcas, idx_th);
#endif
              fic.ouvrir(s, mode);
              snprintf(s, 1000, "%.16e ", current_time);
              fic << s;
              for (int i = 0; i < n; i++)
                {
                  snprintf(s, 1000, "%.16e ", phin_per_bubble[i]);
                  fic << s;
                }
              fic << finl;
              fic.close();

              Cerr << "Fin de l'ecriture des stats par bulles pour la temperature " << idx_th << finl;
            }

          ++idx_th;
        }
    }

  if (Process::je_suis_maitre())
    {
      char s[1000];
      const char *nomcas = nom_cas;
      SFichier fic;
      const int n = position.dimension(0);
      IOS_OPEN_MODE mode = (reset) ? ios::out : ios::app;

      snprintf(s, 1000, "%s_bulles_pousseex.out", nomcas);
      // Cerr << "Ecriture des donnees par bulles: fichier " << s << finl;
      fic.ouvrir(s, mode);
      snprintf(s, 1000, "%.16e ", current_time);
      fic << s;
      for (int i = 0; i < n; i++)
        {
          snprintf(s, 1000, "%.16e ", poussee(i, 0));
          fic << s;
        }
      fic << finl;
      fic.close();

      snprintf(s, 1000, "%s_bulles_hx.out", nomcas);
      // Cerr << "Ecriture des donnees par bulles: fichier " << s << finl;
      fic.ouvrir(s, mode);
      snprintf(s, 1000, "%.16e ", current_time);
      fic << s;
      for (int i = 0; i < n; i++)
        {
          snprintf(s, 1000, "%.16e ", hauteurs_bulles(i, 0));
          fic << s;
        }
      fic << finl;
      fic.close();

      snprintf(s, 1000, "%s_bulles_hy.out", nomcas);
      // Cerr << "Ecriture des donnees par bulles: fichier " << s << finl;
      fic.ouvrir(s, mode);
      snprintf(s, 1000, "%.16e ", current_time);
      fic << s;
      for (int i = 0; i < n; i++)
        {
          snprintf(s, 1000, "%.16e ", hauteurs_bulles(i, 1));
          fic << s;
        }
      fic << finl;
      fic.close();

      snprintf(s, 1000, "%s_bulles_hz.out", nomcas);
      // Cerr << "Ecriture des donnees par bulles: fichier " << s << finl;
      fic.ouvrir(s, mode);
      snprintf(s, 1000, "%.16e ", current_time);
      fic << s;
      for (int i = 0; i < n; i++)
        {
          snprintf(s, 1000, "%.16e ", hauteurs_bulles(i, 2));
          fic << s;
        }
      fic << finl;
      fic.close();

      snprintf(s, 1000, "%s_bulles_centre_x.out", nomcas);
      // Cerr << "Ecriture des donnees par bulles: fichier " << s << finl;
      fic.ouvrir(s, mode);
      snprintf(s, 1000, "%.16e ", current_time);
      fic << s;
      for (int i = 0; i < n; i++)
        {
          snprintf(s, 1000, "%.16e ", position(i, 0));
          fic << s;
        }
      fic << finl;
      fic.close();

      snprintf(s, 1000, "%s_bulles_centre_y.out", nomcas);
      fic.ouvrir(s, mode);
      snprintf(s, 1000, "%.16e ", current_time);
      fic << s;
      for (int i = 0; i < n; i++)
        {
          snprintf(s, 1000, "%.16e ", position(i, 1));
          fic << s;
        }
      fic << finl;
      fic.close();

      snprintf(s, 1000, "%s_bulles_centre_z.out", nomcas);
      fic.ouvrir(s, mode);
      snprintf(s, 1000, "%.16e ", current_time);
      fic << s;
      for (int i = 0; i < n; i++)
        {
          snprintf(s, 1000, "%.16e ", position(i, 2));
          fic << s;
        }
      fic << finl;
      fic.close();

      snprintf(s, 1000, "%s_bulles_surface.out", nomcas);
      fic.ouvrir(s, mode);
      snprintf(s, 1000, "%.16e ", current_time);
      fic << s;
      for (int i = 0; i < n; i++)
        {
          snprintf(s, 1000, "%.16e ", surface[i]);
          fic << s;
        }
      fic << finl;
      fic.close();

      snprintf(s, 1000, "%s_bulles_volume.out", nomcas);
      fic.ouvrir(s, mode);
      snprintf(s, 1000, "%.16e ", current_time);
      fic << s;
      for (int i = 0; i < n; i++)
        {
          snprintf(s, 1000, "%.16e ", volume[i]);
          fic << s;
        }
      fic << finl;
      fic.close();

      if (interfaces_.follow_colors())
        {
          const ArrOfInt& colors = interfaces_.get_colors();
          snprintf(s, 1000, "%s_bulles_colors.out", nomcas);
          fic.ouvrir(s, mode);
          snprintf(s, 1000, "%.16e ", current_time);
          fic << s;
          for (int i = 0; i < n; i++)
            {
              snprintf(s, 1000, "%d ", (True_int) colors[i]);
              fic << s;
            }
          fic << finl;
          fic.close();
        }
#if 0
      // Malheureusement, Le tableau individual_forces n'est pas stocke pour l'instant.
      // On ne peut donc pas le post-traiter avec la frequence lue dans le jdd comme les autres.
      for (int idir=0; idir<3; idir++)
        {
          snprintf(s, 1000, "%s_bulles_external_force_%d.out", nomcas, idir);
          fic.ouvrir(s, mode);
          snprintf(s, 1000, "%.16e ", current_time);
          fic << s;
          for (int ib = 0; ib < interfaces_.nb_bulles_reelles_; ib++)
            {
              snprintf(s, 1000,"%.16e ", individual_forces(ib,idir));
              fic << s;
            }
          fic << finl;
          fic.close();
        }
    }
#endif
}
statistiques().end_count(postraitement_counter_);

}

// Methode qui met a jour l'indicatrice, les termes de repulsion
// ainsi que les termes interfaciaux : ai, kappa*ai, n(aux cellules)
//
// Par definition, mettre igroup a -1 pour inclure toutes les bulles
// Dans ce cas, la methode met a jour l'ev de l'indicatrice au lieu de celui de interfaces_.groups_indicatrice_n_ns()[igroup]
//
// Attention: de nombreux tableaux sont modifies par cette methode en sortie.
// Ils peuvent etre des tableaux de travail. Si on veut qu'il soient correctent
// pour la suite, il faut faire l'appel avec les champs globaux (incluant tous
// les groupes a la fin). Sinon, les champs en ai, normale ou grad_I ne contiendront qu'un groupe.
void IJK_FT_Post::update_stat_ft(const double dt)
{
  static Stat_Counter_Id updtstat_counter_ = statistiques().new_counter(2, "update statistiques");
  statistiques().begin_count(updtstat_counter_);
  if (disable_diphasique_)
    {
      // Calcul du champ grad_P_:
      for (int dir = 0; dir < 3; dir++)
        grad_P_[dir].data() = 0.;

      // pressure gradient requires the "left" value in all directions:
      //  pressure_.echange_espace_virtuel(1 /*, IJK_Field_double::EXCHANGE_GET_AT_LEFT_IJK*/);
      add_gradient_times_constant(pressure_, 1. /*constant*/, grad_P_[0], grad_P_[1], grad_P_[2]);
      for (int dir = 0; dir < 3; dir++)
        grad_P_[dir].echange_espace_virtuel(1);

      statistiques_FT_.update_stat(ref_ijk_ft_, dt);
      return;
    }
  int nb_groups = interfaces_.nb_groups();
  // Boucle debute a -1 pour faire l'indicatrice globale.
  // S'il n'y a pas de groupes de bulles (monophasique ou monodisperse), on passe exactement une fois dans la boucle
  if (nb_groups == 1)
    nb_groups = 0; // Quand il n'y a qu'un groupe, on ne posttraite pas les choses pour ce groupe unique puisque c'est identique au cas global
  for (int igroup = -1; igroup < nb_groups; igroup++)
    {
      interfaces_.calculer_normales_et_aires_interfaciales(ai_ft_, kappa_ai_ft_, normale_cell_ft_, igroup);
      // Puis les redistribue sur le ns :
      ref_ijk_ft_.redistribute_from_splitting_ft_elem_.redistribute(ai_ft_, ai_ns_);
      ref_ijk_ft_.redistribute_from_splitting_ft_elem_.redistribute(kappa_ai_ft_, kappa_ai_ns_);
      ref_ijk_ft_.redistribute_from_splitting_ft_elem_.redistribute(normale_cell_ft_, normale_cell_ns_);
      // (pas besoin d'echange EV car ils n'ont pas de ghost).

      // Calcul du gradient de l'indicatrice et de vitesse :
      if (igroup == -1)
        {
          // interfaces_.In().echange_espace_virtuel(1);
          // Calcul des champs grad_P_, grad_I_ns_
          calculer_gradient_indicatrice_et_pression(interfaces_.In());
          ref_ijk_ft_.transfer_ft_to_ns(); // pour remplir : terme_repulsion_interfaces_ft_ et terme_abs_repulsion_interfaces_ft_
          // Calcul des champs grad_P_, grad_I_ns_, terme_repulsion_interfaces_ns_, terme_abs_repulsion_interfaces_ns_
          // a partir de pressure_, interfaces_.In(), et terme_*_ft_
          statistiques_FT_.update_stat(ref_ijk_ft_, dt);
        }
      else
        {
          // interfaces_.groups_indicatrice_n_ns()[igroup].echange_espace_virtuel(1);
          // Calcul des champs grad_P_, grad_I_ns_, terme_repulsion_interfaces_ns_, terme_abs_repulsion_interfaces_ns_
          // a partir de pressure_, interfaces_.In(), et terme_*_ft_
          calculer_gradient_indicatrice_et_pression(interfaces_.groups_indicatrice_n_ns()[igroup]);
          ref_ijk_ft_.transfer_ft_to_ns();
          groups_statistiques_FT_[igroup].update_stat(ref_ijk_ft_, dt);
        }
    }
  statistiques().end_count(updtstat_counter_);
}

// Calcul du lambda2 a partir du gradient.
// A optimiser simplement en mutualisant avec la methode d'update_stats.
// Et en ne faisant le calcul que si besoin, cad si les champs de gradient ne sont pas a jour...
void IJK_FT_Post::get_update_lambda2()
{
  compute_and_store_gradU_cell(velocity_[0], velocity_[1], velocity_[2],
                               /* Et les champs en sortie */
                               dudx_, dvdy_, dwdx_, dudz_, dvdz_, dwdz_, 1 /* yes compute_all */, dudy_, dvdx_, dwdy_, lambda2_);
}

void IJK_FT_Post::get_update_lambda2_and_rot_and_curl()
{
  get_update_lambda2();
  // Nombre local de mailles en K
  const int kmax = lambda2_.nk();
  const int imax = lambda2_.ni();
  const int jmax = lambda2_.nj();
  for (int k = 0; k < kmax; k++)
    for (int j = 0; j < jmax; j++)
      for (int i = 0; i < imax; i++)
        {
          rot_[0](i, j, k) = dwdy_(i, j, k) - dvdz_(i, j, k);
          rot_[1](i, j, k) = dudz_(i, j, k) - dwdx_(i, j, k);
          rot_[2](i, j, k) = dvdx_(i, j, k) - dudy_(i, j, k);
          // Calcul du critere Q selon (Jeong & Hussain 1995)
          critere_Q_(i, j, k) = -0.5
                                * (dudx_(i, j, k) * dudx_(i, j, k) + 2. * dudy_(i, j, k) * dvdx_(i, j, k) + 2. * dudz_(i, j, k) * dwdx_(i, j, k) + dvdy_(i, j, k) * dvdy_(i, j, k) + 2. * dvdz_(i, j, k) * dwdy_(i, j, k)
                                   + dwdz_(i, j, k) * dwdz_(i, j, k));
        }
}

const FixedVector<IJK_Field_double, 3>& IJK_FT_Post::get_IJK_vector_field(const Nom& nom) const
{
  //int idx_th = 0;
  for (const auto &itr : ref_ijk_ft_.thermique_)
    {
      std::ostringstream oss;
      oss << "GRAD_T"; // << idx_th;
      if (nom == Nom(oss.str().c_str()))
        return itr.grad_T_;
      oss.str("");
      // ++idx_th;
    }
  Cerr << "Erreur dans IJK_FT_Post::get_IJK_vector_field : " << "Champ demande : " << nom << " Liste des champs possibles : " << finl;
  Process::exit();
  throw;
}

bool is_number(const std::string& s)
{
  std::string::const_iterator it = s.begin();
  while (it != s.end() && std::isdigit(*it))
    ++it;
  return !s.empty() && it == s.end();
}

int convert_suffix_to_int(const Nom& nom)
{
  std::string nom_str = nom.getString();
  std::size_t index = nom_str.find_last_of("_");
  std::string suffix("");
  if (index != std::string::npos)
    suffix = nom_str.substr(index + 1);

  bool is_suffix_number = is_number(suffix);
  int suffix_int = -1;
  if (is_suffix_number)
    {
      suffix_int = std::stoi(suffix);
      Cerr << suffix_int << finl;
    }
  return suffix_int;
}

const IJK_Field_double& IJK_FT_Post::get_IJK_field(const Nom& nom) const
{
  if (nom == "PRESSURE_ANA")
    return pressure_ana_;
  if (nom == "PRESSURE_LIQ")
    return extended_pl_;
  if (nom == "PRESSURE_VAP")
    return extended_pv_;

  if (nom == "CELL_VELOCITY_X")
    {
      if (!liste_post_instantanes_.contient_("CELL_VELOCITY"))
        {
          Cerr << "A probe is attempting to access a field CELL_VELOCITY while it has not been computed in the post-processed fields" << finl;
          Process::exit();
        }
      return cell_velocity_[0];
    }
  if (nom == "CELL_VELOCITY_Y")
    {
      if (!liste_post_instantanes_.contient_("CELL_VELOCITY"))
        {
          Cerr << "A probe is attempting to access a field CELL_VELOCITY while it has not been computed in the post-processed fields" << finl;
          Process::exit();
        }
      return cell_velocity_[1];
    }
  if (nom == "CELL_VELOCITY_Z")
    {
      if (!liste_post_instantanes_.contient_("CELL_VELOCITY"))
        {
          Cerr << "A probe is attempting to access a field CELL_VELOCITY while it has not been computed in the post-processed fields" << finl;
          Process::exit();
        }
      return cell_velocity_[2];
    }

  // GAB
  if (nom == "CELL_FORCE_PH_X")
    {
      if (!liste_post_instantanes_.contient_("CELL_FORCE_PH"))
        {
          Cerr << "A probe is attempting to access a field CELL_FORCE_PH while it has not been computed in the post-processed fields" << finl;
          Process::exit();
        }
      return cell_source_spectrale_[0];
    }
  if (nom == "CELL_FORCE_PH_Y")
    {
      if (!liste_post_instantanes_.contient_("CELL_FORCE_PH"))
        {
          Cerr << "A probe is attempting to access a field CELL_FORCE_PH while it has not been computed in the post-processed fields" << finl;
          Process::exit();
        }
      return cell_source_spectrale_[1];
    }
  if (nom == "CELL_FORCE_PH_Z")
    {
      if (!liste_post_instantanes_.contient_("CELL_FORCE_PH"))
        {
          Cerr << "A probe is attempting to access a field CELL_FORCE_PH while it has not been computed in the post-processed fields" << finl;
          Process::exit();
        }
      return cell_source_spectrale_[2];
    }

  // GAB
  if (nom == "CELL_GRAD_P_X")
    {
      if (!liste_post_instantanes_.contient_("CELL_GRAD_P"))
        {
          Cerr << "A probe is attempting to access a field CELL_GRAD_P while it has not been computed in the post-processed fields" << finl;
          Process::exit();
        }
      return cell_grad_p_[0];
    }
  if (nom == "CELL_GRAD_P_Y")
    {
      if (!liste_post_instantanes_.contient_("CELL_GRAD_P"))
        {
          Cerr << "A probe is attempting to access a field CELL_GRAD_P while it has not been computed in the post-processed fields" << finl;
          Process::exit();
        }
      return cell_grad_p_[1];
    }
  if (nom == "CELL_GRAD_P_Z")
    {
      if (!liste_post_instantanes_.contient_("CELL_GRAD_P"))
        {
          Cerr << "A probe is attempting to access a field CELL_GRAD_P while it has not been computed in the post-processed fields" << finl;
          Process::exit();
        }
      return cell_grad_p_[2];
    }
  //

  if (nom.debute_par("dU"))
    {
      const FixedVector<IJK_Field_double, 3>& gradU = statistiques_FT_.get_IJK_vector_field("gradU");
      if (nom == "DUDX")
        return gradU[0];
      if (nom == "DUDY")
        return gradU[1];
      if (nom == "DUDZ")
        return gradU[2];
    }
  if (nom.debute_par("dV"))
    {
      const FixedVector<IJK_Field_double, 3>& gradV = statistiques_FT_.get_IJK_vector_field("gradV");
      if (nom == "DVDX")
        return gradV[0];
      if (nom == "DVDY")
        return gradV[1];
      if (nom == "DVDZ")
        return gradV[2];
    }
  if (nom.debute_par("dW"))
    {
      const FixedVector<IJK_Field_double, 3>& gradW = statistiques_FT_.get_IJK_vector_field("gradW");
      if (nom == "DWDX")
        return gradW[0];
      if (nom == "DWDY")
        return gradW[1];
      if (nom == "DWDZ")
        return gradW[2];
    }

  // GAB, sondes THI
  if (nom.debute_par("FORCE_PH"))
    {
      // A priori inutile : tester si sonde ok avec champ_a_postrer sans FORCE_PH
      // Reponse GB : le vrai test a faire c'est si le field force_ph existe, ie s'il y a un forcage
      if (!liste_post_instantanes_.contient_("FORCE_PH"))
        {
          Cerr << "A probe is attempting to access a field FORCE_PH while it has not been computed in the post-processed fields" << finl;
          Process::exit();
        }
//      FixedVector<IJK_Field_double, 3>& source_spectrale = ref_ijk_ft_.forcage_.get_force_ph2();
      if (nom == "FORCE_PH_X")
        return source_spectrale_[0];
      if (nom == "FORCE_PH_Y")
        return source_spectrale_[1];
      if (nom == "FORCE_PH_Z")
        return source_spectrale_[2];
    }
  //
  // if (disable_diphasique_)
  {
    if (nom == "VELOCITY_ANA_X")
      return velocity_ana_[0];
    if (nom == "VELOCITY_ANA_Y")
      return velocity_ana_[1];
    if (nom == "VELOCITY_ANA_Z")
      return velocity_ana_[2];
    if (nom == "ECART_ANA_X")
      return ecart_ana_[0];
    if (nom == "ECART_ANA_Y")
      return ecart_ana_[1];
    if (nom == "ECART_ANA_Z")
      return ecart_ana_[2];
  }

  const int idx_wanted = convert_suffix_to_int(nom);
  const Nom field_name = nom.getPrefix(Nom("_") + Nom(idx_wanted));
  Cerr << "In get_IJK_field by name : " << nom << " read as : (" << field_name << " ; " << idx_wanted << ")" << finl;

// Remplir la liste de tous les possibles :
  Motcles liste_champs_thermiques_possibles;
  posttraiter_tous_champs_thermique(liste_champs_thermiques_possibles, 0);
  int rang = liste_champs_thermiques_possibles.rang(field_name);
  if (rang == -1)
    {
      Cerr << field_name << " not found as possible for field name. Should be in the list: " << liste_champs_thermiques_possibles << finl;
      Process::exit();
    }

  const Motcle& mot = liste_champs_thermiques_possibles[rang];
  if ((mot == field_name) && (idx_wanted >= 0))
    {
      Cerr << "found as planned " << finl;
    }
  else
    {
      Cerr << "Some issue with the name provided for the sonde. Unrecognised." << finl;
      Process::exit();
    }

  int idx_th = 0;
  for (const auto &itr : ref_ijk_ft_.thermique_)
    {
      if (idx_th == idx_wanted)
        {
          if (field_name == Nom("TEMPERATURE"))
            return itr.temperature_;
          if (field_name == Nom("TEMPERATURE_ANA"))
            return itr.temperature_ana_;
          if (field_name == Nom("SOURCE_TEMPERATURE"))
            return itr.source_temperature_;
          if (field_name == Nom("TEMPERATURE_PHYSIQUE_T"))
            return itr.temperature_physique_T_;
          if (field_name == Nom("T_RUST"))
            return itr.T_rust_;
          if (field_name == Nom("div_rho_cp_T"))
            return itr.div_rho_cp_T_;
          if (field_name == Nom("ECART_T_ANA"))
            return itr.ecart_t_ana_;
          if (field_name == Nom("TEMPERATURE_ADIM_BULLES"))
            return itr.temperature_adim_bulles_;
          if (field_name == Nom("GRAD_T0"))
            return itr.grad_T_[0];
          if (field_name == Nom("GRAD_T1"))
            return itr.grad_T_[1];
          if (field_name == Nom("GRAD_T2"))
            return itr.grad_T_[2];
          break;
        }
      else
        {
          ++idx_th;
        }
    }

  Cerr << "Erreur dans IJK_FT_Post::get_IJK_field : " << finl;
  Cerr << "Champ demande : " << nom << finl;
  Cerr << "Index maximal pour la temperature : " << idx_th - 1 << finl;
  Cerr << "Liste des champs possibles pour la thermique : " << liste_champs_thermiques_possibles << finl;
  Process::exit();
  throw;

}

const int& IJK_FT_Post::get_IJK_flag(const Nom& nom) const
{
  const auto& ii = ref_ijk_ft_.thermique_.front();
  if (nom == "WALL_FLUX")
    return ii.wall_flux_;

  Cerr << "Erreur dans IJK_FT_Post::get_IJK_variable : " << "Variable demandee : " << nom << " Liste des variables possibles : " << finl;
  Process::exit();
  throw;
}

void IJK_FT_Post::sauvegarder_post(const Nom& lata_name)
{
  if (liste_post_instantanes_.contient_("INTEGRATED_VELOCITY"))
    dumplata_vector(lata_name, "INTEGRATED_VELOCITY", integrated_velocity_[0], integrated_velocity_[1], integrated_velocity_[2], 0);

  if (liste_post_instantanes_.contient_("INTEGRATED_PRESSURE"))
    dumplata_scalar(lata_name, "INTEGRATED_PRESSURE", integrated_pressure_, 0);

  if (liste_post_instantanes_.contient_("INDICATRICE_PERTURBE"))
    dumplata_scalar(lata_name, "INDICATRICE_PERTURBE", indicatrice_non_perturbe_, 0);

  if (liste_post_instantanes_.contient_("INTEGRATED_TIMESCALE"))
    dumplata_scalar(lata_name, "INTEGRATED_TIMESCALE", integrated_timescale_, 0);
  // Not necessary, but convenient for picturing...
  if (liste_post_instantanes_.contient_("LAMBDA2"))
    {
      get_update_lambda2();
      dumplata_scalar(lata_name, "LAMBDA2", lambda2_, 0);
    }

}

void IJK_FT_Post::sauvegarder_post_maitre(const Nom& lata_name, SFichier& fichier) const
{
  if (liste_post_instantanes_.contient_("INTEGRATED_VELOCITY"))
    fichier << " fichier_reprise_integrated_velocity " << lata_name << "\n";
  if (liste_post_instantanes_.contient_("INTEGRATED_PRESSURE"))
    fichier << " fichier_reprise_integrated_pressure " << lata_name << "\n";
  if (liste_post_instantanes_.contient_("INTEGRATED_TIMESCALE"))
    fichier << " fichier_reprise_integrated_timescale " << lata_name << "\n";
  if (liste_post_instantanes_.contient_("INDICATRICE_PERTURBE"))
    fichier << " fichier_reprise_indicatrice_non_perturbe " << lata_name << "\n";

  if (statistiques_FT_.t_integration() > 0.)
    {
      Cerr << "All bubbles : " << finl;
      fichier << " statistiques_FT " << statistiques_FT_;
      // S'il y a plusieurs groups, on s'occupe des objets stats pour chaque group:
      // (en ecrivant directement le vecteur d'objets)
      if (interfaces_.nb_groups() > 1)
        {
          Cerr << "Group by group :" << finl;
          fichier << " groups_statistiques_FT " << groups_statistiques_FT_;
        }
    }
}

void IJK_FT_Post::reprendre_post(Param& param)
{
  param.ajouter("statistiques_FT", &statistiques_FT_);
  param.ajouter("groups_statistiques_FT", &groups_statistiques_FT_);

  if (ref_ijk_ft_.coef_immobilisation_ > 1e-16)
    {
      param.ajouter("fichier_reprise_integrated_velocity", &fichier_reprise_integrated_velocity_);
      param.ajouter("fichier_reprise_integrated_pressure", &fichier_reprise_integrated_pressure_);
      param.ajouter("fichier_reprise_integrated_timescale", &fichier_reprise_integrated_timescale_);
      param.ajouter("fichier_reprise_indicatrice_non_perturbe", &fichier_reprise_indicatrice_non_perturbe_);
    }
}

void IJK_FT_Post::fill_op_conv()
{
  if (liste_post_instantanes_.contient_("OP_CONV"))
    for (int i = 0; i < 3; i++)
      op_conv_[i].data() = d_velocity_[i].data();
}

// Calcul du gradient de l'indicatrice et de pression :
//   Attention, il faut que la pression et l'indicatrice soient a jour
//   dans leur espaces virtuels avant d'appeler cette methode
// Methode qui calcule des champs grad_P_, grad_I_ns_,
// a partir de pressure_ et indicatrice_ns__
void IJK_FT_Post::calculer_gradient_indicatrice_et_pression(const IJK_Field_double& indic)
{
  // Remise a zero :
  for (int dir = 0; dir < 3; dir++)
    {
      grad_I_ns_[dir].data() = 0.;
      grad_P_[dir].data() = 0.;
    }

  // From IJK_Navier_Stokes_Tools.cpp
  // interfaces_.In().echange_espace_virtuel(1);
  add_gradient_times_constant(indic, 1. /*Constante multiplicative*/, grad_I_ns_[DIRECTION_I], grad_I_ns_[DIRECTION_J], grad_I_ns_[DIRECTION_K]);

  // pressure gradient requires the "left" value in all directions:
  //  pressure_.echange_espace_virtuel(1 /*, IJK_Field_double::EXCHANGE_GET_AT_LEFT_IJK*/);
  add_gradient_times_constant(pressure_, 1. /*constant*/, grad_P_[0], grad_P_[1], grad_P_[2]);

  for (int dir = 0; dir < 3; dir++)
    {
      grad_I_ns_[dir].echange_espace_virtuel(1);
      grad_P_[dir].echange_espace_virtuel(1);
    }
}

int IJK_FT_Post::alloc_fields()
{
  int nalloc = 0;
  rebuilt_indic_.allocate(splitting_ft_, IJK_Splitting::ELEM, 0);
  potentiel_.allocate(splitting_ft_, IJK_Splitting::ELEM, 0);
  if ((!disable_diphasique_) && ((liste_post_instantanes_.contient_("AIRE_INTERF")) || (liste_post_instantanes_.contient_("TOUS")) || ((t_debut_statistiques_ < 1.e10))))
    {
      ai_ft_.allocate(splitting_ft_, IJK_Splitting::ELEM, 0);
      nalloc += 1;
    }
  // Pour les stats, on calcule kappa*ai :
  if ((!disable_diphasique_) && ((t_debut_statistiques_ < 1.e10)))
    {
      kappa_ai_ft_.allocate(splitting_ft_, IJK_Splitting::ELEM, 0);
      kappa_ai_ns_.allocate(splitting_, IJK_Splitting::ELEM, 0);
      ai_ns_.allocate(splitting_, IJK_Splitting::ELEM, 0);
      allocate_cell_vector(normale_cell_ns_, splitting_, 0);
      nalloc += 6;
    }

  // For the pressure field extension:
  if ((!disable_diphasique_)
      && ((liste_post_instantanes_.contient_("PRESSURE_LIQ")) || (liste_post_instantanes_.contient_("PRESSURE_VAP")) || (liste_post_instantanes_.contient_("TOUS")) || (t_debut_statistiques_ < 1.e10)))
    {
      extended_pressure_computed_ = 1;
      pressure_ft_.allocate(splitting_ft_, IJK_Splitting::ELEM, 5);
      extended_pl_.allocate(splitting_, IJK_Splitting::ELEM, 0);
      extended_pv_.allocate(splitting_, IJK_Splitting::ELEM, 0);
      extended_pl_ft_.allocate(splitting_ft_, IJK_Splitting::ELEM, 0);
      extended_pv_ft_.allocate(splitting_ft_, IJK_Splitting::ELEM, 0);
      nalloc += 5;
    }
  else
    extended_pressure_computed_ = 0;

  // Allocation du champ de normale aux cellules :
  if ((!disable_diphasique_) && ((liste_post_instantanes_.contient_("NORMALE_INTERF")) || (liste_post_instantanes_.contient_("PRESSURE_LIQ")) // Je ne suis pas sur que ce soit necessaire. Seulement si on l'utilise dans le calcul de p_ext
                                 || (liste_post_instantanes_.contient_("PRESSURE_VAP")) // Je ne suis pas sur que ce soit necessaire.
                                 || (liste_post_instantanes_.contient_("TOUS")) || ((t_debut_statistiques_ < 1.e10))))
    {
      allocate_cell_vector(normale_cell_ft_, splitting_ft_, 0);
      nalloc += 3;
    }

  // Allocation des champs derivee de vitesse :
  if ((t_debut_statistiques_ < 1.e10) || (liste_post_instantanes_.contient_("LAMBDA2")) || (liste_post_instantanes_.contient_("CRITERE_Q")) || (liste_post_instantanes_.contient_("CURL")))
    {
      dudx_.allocate(splitting_, IJK_Splitting::ELEM, 1);
      dudy_.allocate(splitting_, IJK_Splitting::ELEM, 1);
      dvdx_.allocate(splitting_, IJK_Splitting::ELEM, 1);
      dvdy_.allocate(splitting_, IJK_Splitting::ELEM, 1);
      dwdx_.allocate(splitting_, IJK_Splitting::ELEM, 1);
      dudz_.allocate(splitting_, IJK_Splitting::ELEM, 1);
      dvdz_.allocate(splitting_, IJK_Splitting::ELEM, 1);
      dwdy_.allocate(splitting_, IJK_Splitting::ELEM, 1);
      dwdz_.allocate(splitting_, IJK_Splitting::ELEM, 1);
      nalloc += 9;
      if (liste_post_instantanes_.contient_("LAMBDA2"))
        {
          lambda2_.allocate(splitting_, IJK_Splitting::ELEM, 0);
          nalloc += 1;
        }
      if ((liste_post_instantanes_.contient_("CRITERE_Q")) || (liste_post_instantanes_.contient_("CURL")))
        {
          critere_Q_.allocate(splitting_, IJK_Splitting::ELEM, 0);
          // Le rotationnel, aux elems aussi :
          allocate_cell_vector(rot_, splitting_, 0);
          nalloc += 4;
        }
    }

  if (interfaces_.nb_groups() > 1)
    {
      // On alloue un tableau assez grand pour contenir tous les groupes.
      if (interfaces_.nb_groups() > 3)
        {
          Cerr << "More than 3 groups are planned, but the allocated fields has only 3 components" << finl;
          Process::exit();
        }
      // TODO AYM: allocate fait dans IJK_Interfaces a l init
      // allocate_cell_vector(interfaces_.groups_indicatrice_n_ns(),splitting_, 1); // Besoin d'un ghost pour le calcul du grad
      // allocate_cell_vector(interfaces_.groups_indicatrice_n_ft(),splitting_ft_, 1); // peut-etre qu'un ghost=0 suffit et fonctionne pour le redistribute.
      nalloc += 6;
    }

  // Pour verification des stats :
  if (check_stats_)
    {
      // Le gradient de pression aux faces : (il est deja alloue par defaut)
      //    allocate_velocity(gradP_, splitting_, 1);
      // Le gradient de vitesse aux elems : (il est stocke, seulement si besoin, dans l'objet statistiques_FT_)
      //    allocate_cell_vector(dUd_, splitting_, 1);
      //    allocate_cell_vector(dVd_, splitting_, 1);
      //    allocate_cell_vector(dWd_, splitting_, 1);
      // Et leurs solutions analytiques sans ghost :
      allocate_velocity(ana_gradP_, splitting_, 1);
      allocate_cell_vector(ana_dUd_, splitting_, 0);
      allocate_cell_vector(ana_dVd_, splitting_, 0);
      allocate_cell_vector(ana_dWd_, splitting_, 0);
      // Pour les deriv secondes :
      allocate_cell_vector(ana_grad2Pi_, splitting_, 0);
      allocate_cell_vector(ana_grad2Pc_, splitting_, 0);
      // And for velocities components :
      allocate_cell_vector(ana_grad2Ui_, splitting_, 0);
      allocate_cell_vector(ana_grad2Uc_, splitting_, 0);
      allocate_cell_vector(ana_grad2Vi_, splitting_, 0);
      allocate_cell_vector(ana_grad2Vc_, splitting_, 0);
      allocate_cell_vector(ana_grad2Wi_, splitting_, 0);
      allocate_cell_vector(ana_grad2Wc_, splitting_, 0);
      nalloc += 36;
    }
  if (liste_post_instantanes_.contient_("NUM_COMPO"))
    {
      num_compo_ft_.allocate(splitting_ft_, IJK_Splitting::ELEM, 0);
      nalloc += 1;
    }
  if (liste_post_instantanes_.contient_("CELL_VELOCITY"))
    {
      allocate_cell_vector(cell_velocity_, splitting_, 0);
      nalloc += 3;
    }
  if (liste_post_instantanes_.contient_("CELL_FORCE_PH"))
    {
      allocate_cell_vector(cell_source_spectrale_, splitting_, 0);
      nalloc += 3;
    }
  if (liste_post_instantanes_.contient_("CELL_GRAD_P"))
    {
      allocate_cell_vector(cell_grad_p_, splitting_, 0);
      nalloc += 3;
    }
  return nalloc;
}

int IJK_FT_Post::alloc_velocity_and_co(bool flag_variable_source)
{
  int n = 0;
  // Le mot cle TOUS n'a pas encore ete compris comme tel.
  if ((liste_post_instantanes_.contient_("GRAD_INDICATRICE_FT")) || (liste_post_instantanes_.contient_("TOUS")))
    n += 3, allocate_velocity(grad_I_ft_, splitting_ft_, 2);

  if ((liste_post_instantanes_.contient_("VELOCITY_ANA")) || (liste_post_instantanes_.contient_("ECART_ANA")))
    n += 3, allocate_velocity(velocity_ana_, splitting_, 1);
  if (liste_post_instantanes_.contient_("ECART_ANA"))
    n += 3, allocate_velocity(ecart_ana_, splitting_, 0);
  if (liste_post_instantanes_.contient_("D_VELOCITY_ANA"))
    n += 3, allocate_velocity(d_velocity_ana_, splitting_, 1);
  if ((liste_post_instantanes_.contient_("PRESSURE_ANA")) || (liste_post_instantanes_.contient_("ECART_P_ANA")))
    n++, pressure_ana_.allocate(splitting_, IJK_Splitting::ELEM, 1);
  if (liste_post_instantanes_.contient_("ECART_P_ANA"))
    n++, ecart_p_ana_.allocate(splitting_, IJK_Splitting::ELEM, 1);
  if (liste_post_instantanes_.contient_("OP_CONV"))
    {
      n += 3, allocate_velocity(op_conv_, splitting_, ref_ijk_ft_.d_velocity_[0].ghost()); // Il y a 1 ghost chez d_velocity_
      //                                          On veut qqch d'aligne pour copier les data() l'un dans l'autre
    }
  // Pour le calcul des statistiques diphasiques :
  // (si le t_debut_stat a ete initialise... Sinon, on ne va pas les calculer au cours de ce calcul)
  if ((t_debut_statistiques_ < 1.e10))
    {
      n += 3, allocate_velocity(grad_I_ns_, splitting_, 1);
      n += 3, allocate_velocity(grad_P_, splitting_, 1);
    }
  else if (flag_variable_source)
    {
      n += 3, allocate_velocity(grad_I_ns_, splitting_, 1);
    }
  return n;
}

void IJK_FT_Post::completer_sondes()
{
  les_sondes_.completer_IJK(ref_ijk_ft_);
}

void IJK_FT_Post::postraiter_sondes()
{
  les_sondes_.postraiter();
}

void IJK_FT_Post::improved_initial_pressure_guess(bool imp)
{
  if ((imp) || (liste_post_instantanes_.contient_("COORDS")))
    {
      Noms noms_coords; // on attend trois expressions
      noms_coords.dimensionner_force(3);
      noms_coords[0] = "X";
      noms_coords[1] = "Y";
      noms_coords[2] = "Z";
      allocate_velocity(coords_, splitting_, 1);
      for (int i = 0; i < 3; i++)
        {
          // Cette methode parcours ni(), nj() et nk() et donc pas les ghost...
          set_field_data(coords_[i], noms_coords[i]);
        }
    }

}

void IJK_FT_Post::postraiter_ci(const Nom& lata_name, const double current_time)
{
  dumplata_header(lata_name);
  dumplata_add_geometry(lata_name, velocity_[0]);
  dumplata_add_geometry(lata_name, ref_ijk_ft_.velocity_ft_[0]);

  // Calcul des moyennes spatiales sur la condition initiale:
  if (current_time >= t_debut_statistiques_)
    {
      // FA AT 16/07/2013 pensent que necessaire pour le calcul des derivees dans statistiques_.update_stat_k(...)
      // Je ne sais pas si c'est utile, mais j'assure...
      velocity_[0].echange_espace_virtuel(2 /*, IJK_Field_ST::EXCHANGE_GET_AT_RIGHT_I*/);
      velocity_[1].echange_espace_virtuel(2 /*, IJK_Field_ST::EXCHANGE_GET_AT_RIGHT_J*/);
      velocity_[2].echange_espace_virtuel(2 /*, IJK_Field_ST::EXCHANGE_GET_AT_RIGHT_K*/);
      pressure_.echange_espace_virtuel(1);

      // C'est update_stat_ft qui gere s'il y a plusieurs groupes
      // pour faire la vraie indicatrice + les groupes
      update_stat_ft(0.);
    }
  else if (!(liste_post_instantanes_.contient_("CURL")) && !(liste_post_instantanes_.contient_("CRITERE_Q")) && (liste_post_instantanes_.contient_("LAMBDA2")))
    {
      // On ne calcul pas encore les stats, mais on veut post-traiter Lambda2 seulement...
      get_update_lambda2();
    }
  else if ((liste_post_instantanes_.contient_("CURL")) || (liste_post_instantanes_.contient_("CRITERE_Q")) || (liste_post_instantanes_.contient_("LAMBDA2")))
    {
      // On ne calcul pas encore les stats, mais on veut deja post-traiter le rotationnel ou Lambda2 ou critere_Q...
      get_update_lambda2_and_rot_and_curl();
    }
}

void IJK_FT_Post::postraiter_fin(bool stop, int tstep, double current_time, double timestep, const Nom& lata_name, const ArrOfDouble& gravite, const Nom& nom_cas)
{
  if (tstep % dt_post_ == dt_post_ - 1 || stop)
    {
      Cout << "tstep : " << tstep << finl;
      posttraiter_champs_instantanes(lata_name, current_time, tstep);
    }
  if (tstep % dt_post_stats_bulles_ == dt_post_stats_bulles_ - 1 || stop)
    {
      ecrire_statistiques_bulles(0, nom_cas, gravite, current_time);
    }
  if (tstep % dt_post_stats_plans_ == dt_post_stats_plans_ - 1 || stop)
    {
      if (current_time >= t_debut_statistiques_)
        posttraiter_statistiques_plans(current_time);
    }
  // Pour ne post-traiter qu'a la frequence demandee :
  les_sondes_.mettre_a_jour(current_time, timestep);
  // Pour le faire a chaque pas de temps :
  // les_sondes_.postraiter();
}

// NEW INTERPOLATING FUNCTION TO AVOID EULERIAN POINTS LYING IN THE BUBBLES OR ON THE INTERFACE
static void ijk_interpolate_implementation_bis(const IJK_Field_double& field, const DoubleTab& coordinates, ArrOfDouble& result, int skip_unknown_points, double value_for_bad_points,
                                               const IJK_Field_double& indic)
{
  const int ghost = field.ghost();
  const int ni = field.ni();
  const int nj = field.nj();
  const int nk = field.nk();
  const IJK_Splitting& splitting = field.get_splitting();
  const IJK_Grid_Geometry& geom = splitting.get_grid_geometry();
  const double dx = geom.get_constant_delta(DIRECTION_I);
  const double dy = geom.get_constant_delta(DIRECTION_J);
  const double dz = geom.get_constant_delta(DIRECTION_K);
  const IJK_Splitting::Localisation loc = field.get_localisation();
  double origin_x = geom.get_origin(DIRECTION_I) + ((loc == IJK_Splitting::FACES_J || loc == IJK_Splitting::FACES_K || loc == IJK_Splitting::ELEM) ? (dx * 0.5) : 0.);
  double origin_y = geom.get_origin(DIRECTION_J) + ((loc == IJK_Splitting::FACES_K || loc == IJK_Splitting::FACES_I || loc == IJK_Splitting::ELEM) ? (dy * 0.5) : 0.);
  double origin_z = geom.get_origin(DIRECTION_K) + ((loc == IJK_Splitting::FACES_I || loc == IJK_Splitting::FACES_J || loc == IJK_Splitting::ELEM) ? (dz * 0.5) : 0.);
  const int nb_coords = coordinates.dimension(0);
  const int gh = indic.ghost();
  result.resize_array(nb_coords);
  for (int idx = 0; idx < nb_coords; idx++)
    {
      const double x = coordinates(idx, 0);  // coordinate of the point where the interpolation is needed
      const double y = coordinates(idx, 1);
      const double z = coordinates(idx, 2);
      const double x2 = (x - origin_x) / dx;
      const double y2 = (y - origin_y) / dy;
      const double z2 = (z - origin_z) / dz;
      const int index_i = (int) (floor(x2)) - splitting.get_offset_local(DIRECTION_I);  //index of the closest cell center is this defined only in the single processor?
      const int index_j = (int) (floor(y2)) - splitting.get_offset_local(DIRECTION_J);
      const int index_k = (int) (floor(z2)) - splitting.get_offset_local(DIRECTION_K);
      const int index_kp1 = index_k + 1;
      const double xfact = x2 - floor(x2);  // non-dimension distance inside the cell where the interpolation is requested
      const double yfact = y2 - floor(y2);
      const double zfact = z2 - floor(z2);

      // Compute the phase indicator function in all the eulerian points involved in the interpolation
      double c_1 = 0.;
      double c_2 = 0.;
      double c_3 = 0.;
      double c_4 = 0.;
      double c_5 = 0.;
      double c_6 = 0.;
      double c_7 = 0.;
      double c_8 = 0.;
      //The following loop is used to take into account the non-periodic condition along the z-axis
      // if the index of the point is outside the wall the cell is considered as vapor, therefore putting an invalid value
      if (!geom.get_periodic_flag(DIRECTION_K))
        // In this case the two domains (FT and NS) are equal in the z direction
        {
          const int kmin = splitting.get_offset_local(DIRECTION_K);
          const int nktot = splitting.get_nb_items_global(loc, DIRECTION_K);
          //Serial calculation
          if (nktot == nk)
            {
              if ((kmin == 0) || (kmin + nk == nktot)) //The conditions on walls may be imposed at the same time since there is no misunderstanding on the indices
                {
                  // The phase indicator function is correctly evaluated only inside the domain of the single processor
                  // otherwise its value remains equal to 0, leading to ad invalid value for the interpolation
                  if (index_k >= 0 && index_k < nk && index_kp1 >= 0 && index_kp1 < nk && index_i >= 0 && index_i < ni && index_j >= 0 && index_j < nj)
                    {
                      c_1 = indic(index_i, index_j, index_k);
                      if (index_i != ni - 1)
                        c_2 = indic(index_i + 1, index_j, index_k);
                      if (index_j != nj - 1)
                        c_3 = indic(index_i, index_j + 1, index_k);
                      if (index_i != ni - 1 && index_j != nj - 1)
                        c_4 = indic(index_i + 1, index_j + 1, index_k);

                      c_5 = indic(index_i, index_j, index_kp1);
                      if (index_i != ni - 1)
                        c_6 = indic(index_i + 1, index_j, index_kp1);
                      if (index_j != nj - 1)
                        c_7 = indic(index_i, index_j + 1, index_kp1);
                      if (index_i != ni - 1 && index_j != nj - 1)
                        c_8 = indic(index_i + 1, index_j + 1, index_kp1);
                    }
                }
            }
          //Parallel calculation
          //NB The following procedure does not fit for one processor
          //since the first condition assigns a value to the point outside the right wall which should be invalid
          //THE LOOP SHOULD BE IMPROVED
          //The two walls are addressed separately since they are elaborated by two different processors.
          else
            {
              if (kmin == 0) //on the left wall
                {
                  if (index_k >= 0 && index_kp1 >= 0 && index_i >= -gh && index_i < ni + gh && index_j >= -gh && index_j < nj + gh && index_kp1 < nk && index_k < nk)
                    {
                      c_1 = indic(index_i, index_j, index_k);
                      if (index_i != ni + gh - 1)
                        c_2 = indic(index_i + 1, index_j, index_k);
                      if (index_j != nj + gh - 1)
                        c_3 = indic(index_i, index_j + 1, index_k);
                      if (index_i != ni + gh - 1 && index_j != nj + gh - 1)
                        c_4 = indic(index_i + 1, index_j + 1, index_k);

                      c_5 = indic(index_i, index_j, index_kp1);
                      if (index_i != ni + gh - 1)
                        c_6 = indic(index_i + 1, index_j, index_kp1);
                      if (index_j != nj + gh - 1)
                        c_7 = indic(index_i, index_j + 1, index_kp1);
                      if (index_i != ni + gh - 1 && index_j != nj + gh - 1)
                        c_8 = indic(index_i + 1, index_j + 1, index_kp1);
                    }
                }
              else if (kmin + nk == nktot) // on the right wall
                {
                  if (index_k < nk && index_kp1 < nk && index_i >= -gh && index_i < ni + gh && index_j >= -gh && index_j < nj + gh && index_kp1 >= 0 && index_k >= 0)
                    {
                      //Process::Journal() << "Proc: " << Process::me() << " i,j,k: "
                      //                   << index_i << " " << index_j << " " << index_k << finl;
                      c_1 = indic(index_i, index_j, index_k);
                      if (index_i != ni + gh - 1)
                        c_2 = indic(index_i + 1, index_j, index_k);
                      if (index_j != nj + gh - 1)
                        c_3 = indic(index_i, index_j + 1, index_k);
                      if (index_i != ni + gh - 1 && index_j != nj + gh - 1)
                        c_4 = indic(index_i + 1, index_j + 1, index_k);

                      c_5 = indic(index_i, index_j, index_kp1);
                      if (index_i != ni + gh - 1)
                        c_6 = indic(index_i + 1, index_j, index_kp1);
                      if (index_j != nj + gh - 1)
                        c_7 = indic(index_i, index_j + 1, index_kp1);
                      if (index_i != ni + gh - 1 && index_j != nj + gh - 1)
                        c_8 = indic(index_i + 1, index_j + 1, index_kp1);
                    }
                }
              else  // points in processor not treating walls
                {
                  if (index_k < nk + gh && index_kp1 < nk + gh && index_i >= -gh && index_i < ni + gh && index_j >= -gh && index_j < nj + gh && index_kp1 >= -gh && index_k >= -gh)
                    {
                      c_1 = indic(index_i, index_j, index_k);
                      if (index_i < ni + gh - 1)
                        c_2 = indic(index_i + 1, index_j, index_k);
                      if (index_j < nj + gh - 1)
                        c_3 = indic(index_i, index_j + 1, index_k);
                      if (index_i < ni + gh - 1 && index_j < nj + gh - 1)
                        c_4 = indic(index_i + 1, index_j + 1, index_k);

                      c_5 = indic(index_i, index_j, index_kp1);
                      if (index_i < ni + gh - 1)
                        c_6 = indic(index_i + 1, index_j, index_kp1);
                      if (index_j < nj + gh - 1)
                        c_7 = indic(index_i, index_j + 1, index_kp1);
                      if (index_i < ni + gh - 1 && index_j < nj + gh - 1)
                        c_8 = indic(index_i + 1, index_j + 1, index_kp1);
                    }
                }
            }
        }

      else // ghost cells may be interrogated when the direction is periodic
        {
          if (index_k < nk + gh && index_kp1 < nk + gh && index_i >= -gh && index_i < ni + gh && index_j >= -gh && index_j < nj + gh && index_kp1 >= -gh && index_k >= -gh)
            {
              c_1 = indic(index_i, index_j, index_k);
              if (index_i < ni + gh - 1)
                c_2 = indic(index_i + 1, index_j, index_k);
              if (index_j < nj + gh - 1)
                c_3 = indic(index_i, index_j + 1, index_k);
              if (index_i < ni + gh - 1 && index_j < nj + gh - 1)
                c_4 = indic(index_i + 1, index_j + 1, index_k);

              c_5 = indic(index_i, index_j, index_kp1);
              if (index_i < ni + gh - 1)
                c_6 = indic(index_i + 1, index_j, index_kp1);
              if (index_j < nj + gh - 1)
                c_7 = indic(index_i, index_j + 1, index_kp1);
              if (index_i < ni + gh - 1 && index_j < nj + gh - 1)
                c_8 = indic(index_i + 1, index_j + 1, index_kp1);
            }
        }

      // Where at least of the points used in the interpolation is outside the liquid phase (at the interface or in the vapor phase)
      // the invalid value is assigned by the interpolation function
      bool ok_2 = (c_1 + c_2 + c_3 + c_4 + c_5 + c_6 + c_7 + c_8 >= 7.99);

// is point in the domain ? (ghost cells ok...)
      bool ok = (index_i >= -ghost && index_i < ni + ghost - 1) && (index_j >= -ghost && index_j < nj + ghost - 1) && (index_k >= -ghost && index_k < nk + ghost - 1);
      if (!ok or !ok_2)
        {
          if (skip_unknown_points)
            {
              result[idx] = value_for_bad_points;
              continue; // go to next point
            }
          else
            {
              // Error!
              Cerr << "Error in ijk_interpolate_implementation: request interpolation of point " << x << " " << y << " " << z << " which is outside of the domain on processor " << Process::me()
                   << finl;
              Process::exit();
            }
        }

      double r = (((1. - xfact) * field(index_i, index_j, index_k) + xfact * field(index_i + 1, index_j, index_k)) * (1. - yfact)
                  + ((1. - xfact) * field(index_i, index_j + 1, index_k) + xfact * field(index_i + 1, index_j + 1, index_k)) * (yfact)) * (1. - zfact)
                 + (((1. - xfact) * field(index_i, index_j, index_kp1) + xfact * field(index_i + 1, index_j, index_kp1)) * (1. - yfact)
                    + ((1. - xfact) * field(index_i, index_j + 1, index_kp1) + xfact * field(index_i + 1, index_j + 1, index_kp1)) * (yfact)) * (zfact);
      result[idx] = r;
    }
}
void ijk_interpolate_skip_unknown_points_bis(const IJK_Field_double& field, const DoubleTab& coordinates, ArrOfDouble& result, const double value_for_bad_points, const IJK_Field_double& indic)
{
  ijk_interpolate_implementation_bis(field, coordinates, result, 1 /* yes:skip unknown points */, value_for_bad_points, indic);
}

// Here begins the computation of the extended pressure fields
// ALL the calculations are performed on the extended domain (FT) and the final field is the redistributed on the physical one (Navier-Stokes)
void IJK_FT_Post::compute_extended_pressures(const Maillage_FT_IJK& mesh)
{
  if (!extended_pressure_computed_)
    return; // Leave the function if the extended fields are not necessary...

  statistiques().begin_count(postraitement_counter_);
  // The following calculation is defined on the extended domain ft
  const IJK_Splitting& split_ft = splitting_ft_;
  //const IJK_Grid_Geometry& geom = split_ft.get_grid_geometry();
  const int ni = split_ft.get_nb_elem_local(DIRECTION_I);
  const int nj = split_ft.get_nb_elem_local(DIRECTION_J);
  const int nk = split_ft.get_nb_elem_local(DIRECTION_K);
  const double dx = split_ft.get_grid_geometry().get_constant_delta(DIRECTION_I);
  const double dy = split_ft.get_grid_geometry().get_constant_delta(DIRECTION_J);
  const double dz = split_ft.get_grid_geometry().get_constant_delta(DIRECTION_K);
  //double origin_x = geom.get_origin(DIRECTION_I);
// double origin_y = geom.get_origin(DIRECTION_J);
// double origin_z = geom.get_origin(DIRECTION_K);

  /*
   void IJK_FT_Post::compute_extended_pressures(const Maillage_FT_IJK& mesh,
   const IJK_Field_double& indic,
   const int phase)
   {
   //const IJK_Field_double& pressure = pressure_;
   IJK_Field_double& extended_p = phase == 1? extended_pl_ : extended_pv_;
   // The field we want to fill :
   if (phase < 0 || phase > 1)
   Process::exit();

   Cerr << extended_p(0,0,0) << finl;


   // If we need it, mesh must not be "const"
   //  const DoubleTab& normale_facettes = mesh.get_update_normale_facettes();

   / Do we need a global table for it?
   interfaces_.calculer_normales_et_aires_interfaciales(ai_,
   kappa_ai_,
   normale_cell_ft_,
   -1) ;
   */
  //onst IJK_Splitting& split = ref_ijk_ft_.splitting_;
  //const int ni = split.get_nb_elem_local(DIRECTION_I) ;
  //const int nj = split.get_nb_elem_local(DIRECTION_J) ;
  //const int nk = split.get_nb_elem_local(DIRECTION_K) ;
  //const double dx = split.get_grid_geometry().get_constant_delta(DIRECTION_I);
  //const double dy = split.get_grid_geometry().get_constant_delta(DIRECTION_J);
  //const double dz = split.get_grid_geometry().get_constant_delta(DIRECTION_K);
  //const double delta[3] = {dx,dy,dz};
  int nbsom = 0;
  // ArrOfInt liste_composantes_connexes_dans_element;
  // liste_composantes_connexes_dans_element.set_smart_resize(1);
  DoubleTab positions_liq(2 * nbsom, 3); // Table of coordinates where interpolation needs to be computed
  DoubleTab positions_vap(2 * nbsom, 3);
  IntTab crossed_cells(nbsom, 3); // Table to store i,j,k of cells crossed by the interface.
  positions_liq.set_smart_resize(1);
  positions_vap.set_smart_resize(1);
  crossed_cells.set_smart_resize(1);

  //pressure field has to be extended from ns to ft
  ref_ijk_ft_.redistribute_to_splitting_ft_elem_.redistribute(pressure_, pressure_ft_);
  pressure_ft_.echange_espace_virtuel(pressure_ft_.ghost()); // 5 ghost cells needed to avoid invalid points in the vapor phase
  //dP_ft_.echange_espace_virtuel(gradP_ft_.ghost());

  //initialisation

  // The other points of the domain ( the one outside the crossed cells are still equal to original value of pressure)
  // This generates the discontinuities seen in Visit
  extended_pl_ft_ = pressure_ft_;
  extended_pv_ft_ = pressure_ft_;
  //for (int dir = 0; dir < 3; dir++)
  //  {
  // dP_ft_[dir].data() = 0.;
  //}

  int errcount_pext = 0;
  // i,j,k are the indices if the cells in the extended domain, for each processor
  for (int k = 0; k < nk; k++)
    {
      for (int j = 0; j < nj; j++)
        {
          for (int i = 0; i < ni; i++)
            {
              if ((interfaces_.In_ft()(i, j, k) > 1.e-6) && (1. - interfaces_.In_ft()(i, j, k) > 1.e-6))
                {
                  // Cerr << "Indicatrice[" << i << ", " << j << ", " << k << "] = " << interfaces_.In_ft()(i,j,k) << finl;
                  //The normal may only be computed on the extended domain
                  // A relationship between the indices of the two meshes is defined
                  // Non sono piu necessari perche si lavora solo  una griglia, quella estesa
                  // const int n_ext = ref_ijk_ft_.get_splitting_extension();
                  // const int i_ft = i+ n_ext;
                  // const int j_ft = j+ n_ext;
                  // int k_ft = k;
                  // if (split.get_grid_geometry().get_periodic_flag(DIRECTION_K))
                  //   k_ft+= n_ext;
                  // const int elem = split_ft.convert_ijk_cell_to_packed(i, j, k);
                  // const int nb_compo_traversantes = interfaces_.compute_list_compo_connex_in_element(mesh, elem, liste_composantes_connexes_dans_element); //number of bubbles crossing the cell
                  Vecteur3 bary_facettes_dans_elem;
                  Vecteur3 normale;
                  double norm = 1.;
                  double dist = 0.;
                  // int num_compo;
                  // if ( nb_compo_traversantes == 0)
                  //   {
                  //     normale[0]=1.;
                  //     normale[1]=1.;
                  //     normale[2]=1.;
                  //     Cerr << "Error no compo traversante on proc. " << Process::me() << finl;
                  //     Cerr << "Indicatrice[" << i <<"," << j << "," << k << "] = " << interfaces_.In_ft()(i,j,k) << finl;
                  //     Process::exit();
                  //   }
                  // else if ( nb_compo_traversantes == 1)
                  //   {

                  //     num_compo = liste_composantes_connexes_dans_element[0];
                  //     interfaces_.calculer_normale_et_bary_element_pour_compo(num_compo,
                  //                                                             elem,
                  //                                                             mesh,
                  //                                                             normale,
                  //                                                             bary_facettes_dans_elem);
                  //   }

                  // // If the same cell is crossed by several bubbles the image points coincide with the crossed cells
                  // // the interpolation function in these points will lead to invalid values.
                  // else
                  //   {
                  //     num_compo = liste_composantes_connexes_dans_element[1];
                  //     interfaces_.calculer_normale_et_bary_element_pour_compo(num_compo,
                  //                                                             elem,
                  //                                                             mesh,
                  //                                                             normale,
                  //                                                             bary_facettes_dans_elem);
                  //   }
                  for (int c = 0; c < 3; c++)
                    {
                      normale[c] = interfaces_.get_norm_par_compo_itfc_in_cell_ft()[c](i, j, k);
                      bary_facettes_dans_elem[c] = interfaces_.get_bary_par_compo_itfc_in_cell_ft()[c](i, j, k);
                    }
                  norm = sqrt(normale[0] * normale[0] + normale[1] * normale[1] + normale[2] * normale[2]);
                  //if (norm<0.95)
                  //    Process::Journal() << "[WARNING-NORM] " << "Indicatrice[" << i <<"," << j << "," << k << "] = " << interfaces_.In_ft()(i,j,k)
                  //                       << " norm= " << norm << finl;
                  //  }
                  if (norm < 1.e-8)
                    {
                      // Process::Journal() << " nb_compo_traversantes " << nb_compo_traversantes << finl;
                      Process::Journal() << "Indicatrice[" << i << "," << j << "," << k << "] = " << interfaces_.In_ft()(i, j, k) << finl;
                      Process::Journal() << "[WARNING-Extended-pressure] on Proc. " << Process::me() << "Floating Point Exception is barely avoided (" << " normale " << normale[0] << " " << normale[1]
                                         << " " << normale[2] << " )" << finl;
                      Process::Journal() << " But we have no distance to extrapolate the pressure" << finl;
                      dist = 1.52 * sqrt(dx * dx + dy * dy + dz * dz) / 3.;
                      if (interfaces_.In_ft()(i, j, k) * (1 - interfaces_.In_ft()(i, j, k)) > 1.e-6)
                        {
                          Process::Journal() << "[WARNING-Extended-pressure] " << "Indicatrice[" << i << "," << j << "," << k << "] = " << interfaces_.In_ft()(i, j, k) << finl;
                          if (interfaces_.In_ft()(i, j, k) > 0.99)
                            {
                              Process::Journal() << "[WARNING-Extended-pressure] " << "Pressure_ft_ will be kept as an extension for p_liq pressure_[" << i << "," << j << "," << k << "] = "
                                                 << pressure_ft_(i, j, k) << finl;
                              extended_pv_ft_(i, j, k) = 1.e20;
                            }
                          else
                            {
                              Process::Journal() << "[WARNING-Extended-pressure] " << "Pressure_ft_ will be kept as an extension for p_vap pressure_[" << i << "," << j << "," << k << "] = "
                                                 << pressure_ft_(i, j, k) << finl;
                              extended_pl_ft_(i, j, k) = 1.e20;
                            }
                          continue; // This cell is not added to crossed cells.
                        }
                      else
                        {
                          // We still need a unit normal
                          for (int dir = 0; dir < 3; dir++)
                            if (normale[dir] != 0)
                              normale[dir] = (normale[dir] > 0.) ? 1.0 : -1.0;

                          norm = sqrt(normale[0] * normale[0] + normale[1] * normale[1] + normale[2] * normale[2]);
                          if (std::fabs(norm) < 1.e-10)
                            {
                              Process::Journal() << "[WARNING-Extended-pressure] ||normal|| < 1.e-10" << finl;
                            }
                        }
                    }

                  if (std::fabs(norm) < 1.e-10)
                    {
                      Process::Journal() << "[WARNING-Extended-pressure] Even with precaution, the normal truely is zero in compute_extended_pressures()... " << finl;
                      Cerr << "We ignore the extrapolation and keep the local value... (hopefully rare enough)" << finl;
                      dist = 0.;
                      errcount_pext++;
                    }
                  else
                    {
                      dist = 1.52 * (std::fabs(dx * normale[0]) + std::fabs(dy * normale[1]) + std::fabs(dz * normale[2])) / norm;
                      // 2020.04.15 : GB correction for non-isotropic meshes and closer extrapolation points:
                      // The previous version was looking very far away in the direction where the mesh is fine (dz in channel flows).
                      // Then, a lot of ghost would be required.
                      // This new version still goes to the same value of dist when nx=1 or when nz=1, but is sin between
                      // double eps = 1.e-20;
                      // dist = 1.52*std::min(min(dx/(std::fabs(normale[0])+eps),dy/(std::fabs(normale[1])+eps)),dz/(std::fabs(normale[2])+eps));
                    }

                  nbsom++;
                  crossed_cells.resize(nbsom, 3, Array_base::COPY_INIT);
                  positions_liq.resize(2 * nbsom, 3, Array_base::COPY_INIT);
                  positions_vap.resize(2 * nbsom, 3, Array_base::COPY_INIT);

                  crossed_cells(nbsom - 1, 0) = i;
                  crossed_cells(nbsom - 1, 1) = j;
                  crossed_cells(nbsom - 1, 2) = k;

                  for (int dir = 0; dir < 3; dir++)
                    {
                      // Four image points are calculated, two on each side of the interface
                      //liquid phase
                      positions_liq(2 * nbsom - 2, dir) = bary_facettes_dans_elem[dir] + dist * normale[dir]; // 1st point to be done...
                      positions_liq(2 * nbsom - 1, dir) = bary_facettes_dans_elem[dir] + 2 * dist * normale[dir]; // 2nd point to be done...
                      //vapor_phase
                      positions_vap(2 * nbsom - 2, dir) = bary_facettes_dans_elem[dir] - dist * normale[dir]; // 1st point to be done...
                      positions_vap(2 * nbsom - 1, dir) = bary_facettes_dans_elem[dir] - 2 * dist * normale[dir]; // 2nd point to be done...
                    }
                }
            }
        }
    }

  errcount_pext = Process::mp_sum(errcount_pext);
  if ((Process::je_suis_maitre()) && (errcount_pext))
    {
      Cerr << "[WARNING-Extended-pressure] Error Count = " << errcount_pext << finl;
    }

// Interpolation on the image points
// All the quantities are evaluated on the extended domain, both the pressure field, both the image points coordinates

  ArrOfDouble p_interp_liq(2 * nbsom);
  ArrOfDouble p_interp_vap(2 * nbsom);
  ijk_interpolate_skip_unknown_points(pressure_ft_, positions_vap, p_interp_vap, 1.e5 /*value for unknown points*/);
  ijk_interpolate_skip_unknown_points_bis(pressure_ft_, positions_liq, p_interp_liq, 1.e5 /* value for unknown points */, interfaces_.In_ft());

// Extrapolation in the eulerian cells crossed by the interface
  int inval_pl_count = 0.;
  int inval_pv_count = 0.;
  for (int icell = 0; icell < nbsom; icell++)
    {
      const int i = crossed_cells(icell, 0);
      const int j = crossed_cells(icell, 1);
      const int k = crossed_cells(icell, 2);
      // const int nb_compo_traversantes = interfaces_.compute_list_compo_connex_in_element(mesh, elem, liste_composantes_connexes_dans_element);
      const int nb_compo_traversantes = interfaces_.nb_compo_traversantes(i, j, k);
      if (nb_compo_traversantes != 1)
        {
          extended_pv_ft_(i, j, k) = 1.e20;
          inval_pv_count++;
        }
      else
        {
          extended_pv_ft_(i, j, k) = 2 * p_interp_vap[2 * icell] - 1 * p_interp_vap[2 * icell + 1];
        }

      //for(int dir=0; dir<3; dir++)
      if (p_interp_liq[2 * icell + 1] == 1.e5)
        {
          if (p_interp_liq[2 * icell] == 1.e5)   //May changing the order affect the result??
            {
              extended_pl_ft_(i, j, k) = 1.e20;
              inval_pl_count++;
              //dP_ft_[dir](i,j,k) = 1.e20
            }
          else
            {
              extended_pl_ft_(i, j, k) = p_interp_liq[2 * icell];
              //dP_ft_[dir](i,j,k) = 1.e20
            }
        }
      else
        {
          if (p_interp_liq[2 * icell] != 1.e5)
            {
              extended_pl_ft_(i, j, k) = 2 * p_interp_liq[2 * icell] - 1 * p_interp_liq[2 * icell + 1];
              //dP_ft_[dir](i,j,k) = (p_interp_liq(2*icell) - 1*p_interp_liq(2*icell+1))*normale[dir]/dist;
            }
          else
            {
              extended_pl_ft_(i, j, k) = p_interp_liq[2 * icell + 1];
              //dP_ft_[dir](i,j,k) = 1.e20
            }
        }
    }

  inval_pl_count = Process::mp_sum(inval_pl_count);
  inval_pv_count = Process::mp_sum(inval_pv_count);
  if ((Process::je_suis_maitre()) && (inval_pl_count))
    Cerr << "[WARNING-Extended-pressure] Invalid p_l cells Count = " << inval_pl_count << finl;
  if ((Process::je_suis_maitre()) && (inval_pv_count))
    Cerr << "[WARNING-Extended-pressure] Invalid p_v cells Count = " << inval_pv_count << finl;

  // The previous evaluated extended pressure has to be recomputed on the real NS domain
  ref_ijk_ft_.redistribute_from_splitting_ft_elem_.redistribute(extended_pl_ft_, extended_pl_);
  ref_ijk_ft_.redistribute_from_splitting_ft_elem_.redistribute(extended_pv_ft_, extended_pv_);
  //ref_ijk_ft_.redistribute_from_splitting_ft_elem_.redistribute(dP_ft, dP_);

  extended_pl_.echange_espace_virtuel(extended_pl_.ghost());
  extended_pv_.echange_espace_virtuel(extended_pv_.ghost());
  statistiques().end_count(postraitement_counter_);
}

#if 0
// copy:
extended_pl_ = pressure_;
extended_pv_ = pressure_;

// Change d_velocity + rho_field_ + mu_field_ ...
compute_phase_pressures_based_on_poisson(0);
compute_phase_pressures_based_on_poisson(1);

//Computing the integral of the extended pressure fields on each bubble so to take into account the real jump condition constitued by the surface tension:
//The values obtained can described better the physics of the problem but are not anymore involved in the pressure extension
const Maillage_FT_IJK& maillage= interfaces_.maillage_ft_ijk_;
const DoubleTab& sommets = maillage.sommets() ; // Table of the coordinates of the markers on the front
const IntTab& facettes = maillage.facettes(); // Table of the nodes of the unstructured mesh of the front
int nbsom = sommets.dimension(0); // Number of  points on the front
int n =maillage.nb_facettes(); // Number of elements on the front
const ArrOfDouble& courbure = maillage.get_update_courbure_sommets();
const ArrOfDouble& surface_facettes = maillage.get_update_surface_facettes();
const int nbulles_reelles = interfaces_.get_nb_bulles_reelles();
const int nbulles_ghost =interfaces_.get_nb_bulles_ghost();
const int nbulles_tot = nbulles_reelles + nbulles_ghost;
ArrOfDouble surface_par_bulle;
interfaces_.calculer_surface_bulles(surface_par_bulle);
const ArrOfInt& compo_connex = maillage.compo_connexe_facettes();



// Extended fields interpolation on the front
ArrOfDouble pl_ext_interp(nbsom);
ArrOfDouble pv_ext_interp(nbsom);
ijk_interpolate_skip_unknown_points(extended_pl_, sommets, pl_ext_interp, 1.e5 /* value for unknown points */);
ijk_interpolate_skip_unknown_points(extended_pv_, sommets, pv_ext_interp, 1.e5 /* value for unknown points */);

// Weighted averaged values on each bubble
ArrOfDouble pl(nbulles_tot);
ArrOfDouble pv(nbulles_tot);
ArrOfDouble k(nbulles_tot);
// Weights based on the surface of the elements
for (int fa7 = 0; fa7 < n; fa7++)
  {
    const double sf=surface_facettes[fa7];
    int compo = compo_connex[fa7];
    if (compo < 0)
      {
        const int idx_ghost = interfaces_.get_ghost_number_from_compo(compo);
        compo = nbulles_reelles - 1 -idx_ghost;
      }
    assert(compo >=0);
    assert(compo < nbulles_tot);
    for (int isom = 0; isom< 3; isom++)
      {
        const int num_som = facettes(fa7, isom);
        const double kappa = courbure(num_som);
        const double pl_ext = pl_ext_interp(num_som);
        const double pv_ext = pv_ext_interp(num_som);
        const double fac = sf/3.;
        k(compo) += kappa*fac;
        pl(compo)+= pl_ext*fac;
        pv(compo)+= pv_ext*fac;
      }
  }
// Imposing the correct jump value
ArrOfDouble source(nbulles_tot);
ArrOfDouble diff(nbulles_tot);
ArrOfDouble pv_2(nbulles_tot);

for (int icompo = 0; icompo < nbulles_tot; icompo++)
  {
    k(icompo) /=surface_par_bulle(icompo);
    pl(icompo) /=surface_par_bulle(icompo);
    pv(icompo) /=surface_par_bulle(icompo);
    diff(icompo) = pl(icompo) - pv(icompo);
    source (icompo) = ref_ijk_ft_.sigma_*k(icompo);
    pv_2(icompo) = pv(icompo)-source(icompo)+diff(icompo);
  }
// Printing on screen
Cout<< "courbure: "  << k;
Cout<< "Liquid field: "  << pl;
Cout<< "Vapor field: "  << pv_2;

}
fichier_reprise_vitesse_

void IJK_FT_Post::compute_phase_pressures_based_on_poisson(const int phase)
// Computes a new d_velocity based on a virtual time step
// Takes into account only convection + diffusion (no interfacial source term)
// Takes constant properties per phase
{
  const double fac = 1.e-4;
  const double virtual_timestep = ref_ijk_ft_.timestep_ * fac;

  // Attribution of the value of each phase everywhere
  IJK_Field_double& rho_field_ = ref_ijk_ft_.rho_field_;
  IJK_Field_double& molecular_mu_ = ref_ijk_ft_.molecular_mu_;
  const double mu_phase = (phase==1) ? ref_ijk_ft_.mu_liquide_ : ref_ijk_ft_.mu_vapeur_;
  const double rho_phase = (phase==1) ? ref_ijk_ft_.rho_liquide_ : ref_ijk_ft_.rho_vapeur_;
  molecular_mu_.data() = mu_phase;
  rho_field_.data() = rho_phase;
  const int store_value_disable_diphasique = ref_ijk_ft_.disable_diphasique_;
  IJK_Field_double& extended_p = (phase==1) ? extended_pl_ : extended_pv_;
  ref_ijk_ft_.disable_diphasique_ = true;

  ref_ijk_ft_.calculer_dv(virtual_timestep, ref_ijk_ft_.current_time_, -1 /* rk_step -> not in RK3 */);
  FixedVector<IJK_Field_double, 3>& dvdt = ref_ijk_ft_.d_velocity_; //
  IJK_Field_double& rhs = ref_ijk_ft_.pressure_rhs_;

  // Poisson solver applied to the calculated projection field
  pressure_projection(dvdt[0], dvdt[1],  dvdt[2], extended_p, 1/rho_phase,
                      rhs, ref_ijk_ft_.check_divergence_, poisson_solver_post_);
  ref_ijk_ft_.disable_diphasique_ = store_value_disable_diphasique;
  {
    // Restore the one-fluid variables
    for (int k=0; k < interfaces_.In().nk() ; k++)
      for (int j=0; j< interfaces_.In().nj(); j++)
        for (int i=0; i < interfaces_.In().ni(); i++)
          {
            double chi_l = interfaces_.In()(i,j,k);
            rho_field_(i,j,k)    = ref_ijk_ft_.rho_liquide_ * chi_l + (1.- chi_l) * ref_ijk_ft_.rho_vapeur_;
            molecular_mu_(i,j,k) = ref_ijk_ft_.mu_liquide_  * chi_l + (1.- chi_l) * ref_ijk_ft_.mu_vapeur_ ;
          }
    rho_field_.echange_espace_virtuel(rho_field_.ghost());
    molecular_mu_.echange_espace_virtuel(molecular_mu_.ghost());
  }
}
#endif

// Methode appelee lorsqu'on a mis "TOUS" dans la liste des champs a postraiter.
// Elle ajoute a la liste tous les noms de champs postraitables par IJK_Interfaces
void IJK_FT_Post::posttraiter_tous_champs_thermique(Motcles& liste, const int idx) const
{
  liste.add("TEMPERATURE");
  liste.add("CP");
  liste.add("LAMBDA");
  liste.add("TEMPERATURE_ANA");
  liste.add("ECART_T_ANA");
  liste.add("DIV_LAMBDA_GRAD_T_VOLUME");
  liste.add("SOURCE_TEMPERATURE");
  liste.add("TEMPERATURE_ADIM_BULLES");
  liste.add("TEMPERATURE_PHYSIQUE_T");
  liste.add("TEMPERATURE_ADIMENSIONNELLE_THETA");
  liste.add("SOURCE_TEMPERATURE_ANA");
  liste.add("ECART_SOURCE_TEMPERATURE_ANA");
  liste.add("GRAD_T");
  liste.add("GRAD_T0");
  liste.add("GRAD_T1");
  liste.add("GRAD_T2");
  liste.add("T_RUST");
  liste.add("DIV_RHO_CP_T_V");

}

// Methode appelee lorsqu'on a mis "TOUS" dans la liste des champs a postraiter.
// Elle ajoute a la liste tous les noms de champs postraitables par IJK_Interfaces
void IJK_FT_Post::posttraiter_tous_champs_energie(Motcles& liste, const int idx) const
{
  liste.add("TEMPERATURE");
  liste.add("CP");
  liste.add("LAMBDA");
  liste.add("TEMPERATURE_ANA");
  liste.add("ECART_T_ANA");
  liste.add("DIV_LAMBDA_GRAD_T_VOLUME");
  liste.add("SOURCE_TEMPERATURE");
  liste.add("TEMPERATURE_ADIM_BULLES");
  liste.add("TEMPERATURE_PHYSIQUE_T");
  liste.add("TEMPERATURE_ADIMENSIONNELLE_THETA");
  liste.add("SOURCE_TEMPERATURE_ANA");
  liste.add("ECART_SOURCE_TEMPERATURE_ANA");
  liste.add("GRAD_T");
  liste.add("GRAD_T0");
  liste.add("GRAD_T1");
  liste.add("GRAD_T2");
  liste.add("T_RUST");
  liste.add("DIV_RHO_CP_T_V");
}

// idx is the number of the temperature in the list
int IJK_FT_Post::posttraiter_champs_instantanes_thermique(const Motcles& liste_post_instantanes, const char *lata_name, const int latastep, const double current_time, IJK_Thermique& itr,
                                                          const int idx)
{
  Cerr << liste_post_instantanes << finl;
  int n = 0; // nombre de champs postraites
  std::ostringstream oss;
  oss << "TEMPERATURE_" << idx;
  Nom nom_temp(oss.str().c_str());
  if ((liste_post_instantanes_.contient_("TEMPERATURE")) || (liste_post_instantanes.contient_(nom_temp)))
    {
      n++, dumplata_scalar(lata_name, nom_temp, itr.temperature_, latastep);
    }
  oss.str("");
  if (liste_post_instantanes.contient_("CP"))
    n++, dumplata_scalar(lata_name, "CP", itr.cp_, latastep);
  if (liste_post_instantanes.contient_("LAMBDA"))
    n++, dumplata_scalar(lata_name, "LAMBDA", itr.lambda_, latastep);

  //MR: temperature analytique
  oss << "TEMPERATURE_" << idx << "_ANA";
  Nom nom_ana(oss.str().c_str());
  if ((liste_post_instantanes_.contient_("TEMPERATURE_ANA")) || (liste_post_instantanes.contient_(nom_ana)))
    {
      //set_field_data(itr.temperature_ana_, itr.expression_T_ana_, current_time);
      itr.set_field_T_ana();
      n++, dumplata_scalar(lata_name, nom_ana, itr.temperature_ana_, latastep);
    }
  oss.str("");

  oss << "GRAD_T_" << idx;
  Nom nom1(oss.str().c_str());

  if (liste_post_instantanes_.contient_("GRAD_T") || (liste_post_instantanes.contient_(nom1)))
    {
      const FixedVector<IJK_Field_double, 3>& grad_T = itr.grad_T_;
      n++, dumplata_vector(lata_name, "GRADT", grad_T[0], grad_T[1], grad_T[2], latastep);
      //  IJK_Field_double& grad_T0 = itr.grad_T_[0];
      //  IJK_Field_double& grad_T1 = itr.grad_T_[1];
      // IJK_Field_double& grad_T2 = itr.grad_T_[2];
      //FixedVector<IJK_Field_double, 3>& grad_T = [grad_T0, grad_T1, grad_T2];
      // n++,dumplata_scalar(lata_name,"dTd", grad_T2, latastep);
    }
  oss.str("");

  oss << "TEMPERATURE_" << idx << "_ADIMENSIONNELLE_THETA";
  Nom nom2(oss.str().c_str());
  if ((liste_post_instantanes_.contient_("TEMPERATURE_ADIMENSIONNELLE_THETA")) || (liste_post_instantanes.contient_(nom2)))
    {
      //  set_field_data(temperature_ana_, itr.expression_T_ana_, current_time);
      n++, dumplata_scalar(lata_name, nom2, itr.temperature_adimensionnelle_theta_, latastep);
    }
  oss.str("");

  oss << "SOURCE_TEMPERATURE_" << idx;
  Nom nom3(oss.str().c_str());
  if ((liste_post_instantanes_.contient_("SOURCE_TEMPERATURE")) || liste_post_instantanes.contient_(nom3))
    {
      // set_field_data(source_temperature_ana_, itr.expression_source_temperature_, velocity_[0], current_time);
      n++, dumplata_scalar(lata_name, nom3, itr.source_temperature_, latastep);
    }
  oss.str("");

  oss << "TEMPERATURE_" << idx << "_ADIM_BULLES";
  Nom nom20(oss.str().c_str());
  if ((liste_post_instantanes_.contient_("TEMPERATURE_ADIM_BULLES")) || (liste_post_instantanes.contient_(nom20)))
    {
      n++, dumplata_scalar(lata_name, nom20, itr.temperature_adim_bulles_, latastep);
    }
  oss.str("");

  oss << "TEMPERATURE_PHYSIQUE_T_" << idx;
  Nom nom6(oss.str().c_str());
  if ((liste_post_instantanes_.contient_("TEMPERATURE_PHYSIQUE_T")) || liste_post_instantanes.contient_(nom6))
    {
      // set_field_data(source_temperature_ana_, itr.expression_source_temperature_, velocity_[0], current_time);
      n++, dumplata_scalar(lata_name, nom6, itr.temperature_physique_T_, latastep);
    }
  oss.str("");

  oss << "ECART_T_ANA" << idx;
  Nom nom4(oss.str().c_str());
  if ((liste_post_instantanes_.contient_("ECART_T_ANA") || liste_post_instantanes.contient_(nom4)))
    {
      itr.calculer_ecart_T_ana();
      n++, dumplata_scalar(lata_name, nom4, itr.ecart_t_ana_, latastep);
    }
  oss.str("");

  oss << "DIV_LAMBDA_GRAD_T_VOLUME" << idx;
  Nom nom5(oss.str().c_str());
  if ((liste_post_instantanes_.contient_("DIV_LAMBDA_GRAD_T_VOLUME") || liste_post_instantanes.contient_(nom5)))
    {
      n++, dumplata_scalar(lata_name, nom5, itr.div_lambda_grad_T_volume_, latastep);
    }
  oss.str("");

  oss << "T_" << idx << "_RUST";
  Nom nom7(oss.str().c_str());
  if (((liste_post_instantanes.contient_("T_RUST")) || (liste_post_instantanes.contient_(nom7))) && (itr.conserv_energy_global_))
    {
      //  set_field_data(temperature_ana_, itr.expression_T_ana_, current_time);
      n++, dumplata_scalar(lata_name, nom7, itr.T_rust_, latastep);
    }
  else if ((liste_post_instantanes.contient_("T_RUST")) || (liste_post_instantanes.contient_(nom7)))
    {
      n++;
    }
  oss.str("");

  oss << "DIV_RHO_CP_T_" << idx << "_V";
  Nom nom8(oss.str().c_str());
  if (((liste_post_instantanes.contient_("DIV_RHO_CP_T_V")) || (liste_post_instantanes.contient_(nom8))) && (itr.type_temperature_convection_form_ == 2))
    {
      //  set_field_data(temperature_ana_, itr.expression_T_ana_, current_time);
      n++, dumplata_scalar(lata_name, nom8, itr.div_rho_cp_T_, latastep);
    }
  else if ((liste_post_instantanes.contient_("DIV_RHO_CP_T_V")) || (liste_post_instantanes.contient_(nom8)))
    {
      n++;
    }
  oss.str("");

  return n;
}

int IJK_FT_Post::posttraiter_champs_instantanes_energie(const Motcles& liste_post_instantanes, const char *lata_name, const int latastep, const double current_time, IJK_Energie& itr, const int idx)
{
  Cerr << liste_post_instantanes << finl;
  int n = 0; // nombre de champs postraites
  std::ostringstream oss;
  oss << "TEMPERATURE_EN_" << idx;
  Nom nom_temp(oss.str().c_str());
  if ((liste_post_instantanes_.contient_("TEMPERATURE")) || (liste_post_instantanes.contient_(nom_temp)))
    {
      n++, dumplata_scalar(lata_name, nom_temp, itr.temperature_, latastep);
    }
  oss.str("");

  //MR: temperature analytique
  oss << "TEMPERATURE_EN_" << idx << "_ANA";
  Nom nom_ana(oss.str().c_str());
  if ((liste_post_instantanes_.contient_("TEMPERATURE_ANA")) || (liste_post_instantanes.contient_(nom_ana)))
    {
      //set_field_data(itr.temperature_ana_, itr.expression_T_ana_, current_time);
      itr.set_field_T_ana();
      n++, dumplata_scalar(lata_name, nom_ana, itr.temperature_ana_, latastep);
    }
  oss.str("");

  oss << "GRAD_T_EN_" << idx;
  Nom nom1(oss.str().c_str());

  if (liste_post_instantanes_.contient_("GRAD_T") || (liste_post_instantanes.contient_(nom1)))
    {
      const FixedVector<IJK_Field_double, 3>& grad_T = itr.grad_T_;
      n++, dumplata_vector(lata_name, "GRADT", grad_T[0], grad_T[1], grad_T[2], latastep);
      //  IJK_Field_double& grad_T0 = itr.grad_T_[0];
      //  IJK_Field_double& grad_T1 = itr.grad_T_[1];
      // IJK_Field_double& grad_T2 = itr.grad_T_[2];
      //FixedVector<IJK_Field_double, 3>& grad_T = [grad_T0, grad_T1, grad_T2];
      // n++,dumplata_scalar(lata_name,"dTd", grad_T2, latastep);
    }
  oss.str("");

  oss << "ECART_T_ANA_EN_" << idx;
  Nom nom4(oss.str().c_str());
  if ((liste_post_instantanes_.contient_("ECART_T_ANA") || liste_post_instantanes.contient_(nom4)))
    {
      itr.calculer_ecart_T_ana();
      n++, dumplata_scalar(lata_name, nom4, itr.ecart_t_ana_, latastep);
    }
  oss.str("");

  oss << "DIV_LAMBDA_GRAD_T_VOLUME_EN_" << idx;
  Nom nom5(oss.str().c_str());
  if ((liste_post_instantanes_.contient_("DIV_LAMBDA_GRAD_T_VOLUME") || liste_post_instantanes.contient_(nom5)))
    {
      n++, dumplata_scalar(lata_name, nom5, itr.div_lambda_grad_T_volume_, latastep);
    }
  oss.str("");
  return n;
}

// idx is the number of the temperature in the list
int IJK_FT_Post::posttraiter_champs_instantanes_thermique_interfaciaux(const Motcles& liste_post_instantanes, const char *lata_name, const int latastep, const double current_time, IJK_Thermique& itr,
                                                                       const int idx)
{
  Cerr << liste_post_instantanes << finl;
  int n = 0; // nombre de champs postraites
  std::ostringstream oss;
  oss << "INTERFACE_TEMPERATURE_" << idx;
  Nom nom_temp(oss.str().c_str());
  std::ostringstream oss2;
  oss2 << "INTERFACE_PHIN_" << idx;
  Nom nom_phin(oss2.str().c_str());
  if ((liste_post_instantanes.contient_("INTERFACE_TEMPERATURE")) || (liste_post_instantanes.contient_("INTERFACE_PHIN")) || (liste_post_instantanes.contient_(nom_temp))
      || (liste_post_instantanes.contient_(nom_phin)))
    {
      //  Computing interfacial temperature at fa7 centre :
      const Maillage_FT_IJK& mesh = interfaces_.maillage_ft_ijk();
      const ArrOfDouble& surface_facettes = mesh.get_update_surface_facettes();
      const int nb_facettes = mesh.nb_facettes();
      ArrOfDouble interfacial_temperature;
      ArrOfDouble interfacial_phin;
      // To transfer the field to FT splitting (because interfaces are there...) !!! NEEDED for compute_interfacial_temperature
      IJK_Field_double& temperature_ft = itr.get_temperature_ft();
      ref_ijk_ft_.redistribute_to_splitting_ft_elem_.redistribute(itr.get_temperature(), temperature_ft);
      temperature_ft.echange_espace_virtuel(temperature_ft.ghost());
      // results are prop to the area :
      //itr.compute_interfacial_temperature(interfacial_temperature, interfacial_phin, itr.get_storage());
      itr.compute_interfacial_temperature2(interfacial_temperature, interfacial_phin);
      for (int fa7 = 0; fa7 < nb_facettes; fa7++)
        {
          const double sf = surface_facettes[fa7];
          interfacial_temperature[fa7] /= sf;
          interfacial_phin[fa7] /= sf;
        }
      if ((liste_post_instantanes.contient_("INTERFACE_TEMPERATURE")) || (liste_post_instantanes.contient_(nom_temp)))
        n++, dumplata_ft_field(lata_name, "INTERFACES", nom_temp, "ELEM", interfacial_temperature, latastep);
      if ((liste_post_instantanes.contient_("INTERFACE_PHIN")) || (liste_post_instantanes.contient_(nom_phin)))
        n++, dumplata_ft_field(lata_name, "INTERFACES", nom_phin, "ELEM", interfacial_phin, latastep);
    }
  oss.str("");
  return n;
}

// idx is the number of the temperature in the list
int IJK_FT_Post::posttraiter_champs_instantanes_energie_interfaciaux(const Motcles& liste_post_instantanes, const char *lata_name, const int latastep, const double current_time, IJK_Energie& itr,
                                                                     const int idx)
{
  Cerr << liste_post_instantanes << finl;
  int n = 0; // nombre de champs postraites
  std::ostringstream oss;
  oss << "INTERFACE_TEMPERATURE_EN_" << idx;
  Nom nom_temp(oss.str().c_str());
  std::ostringstream oss2;
  oss2 << "INTERFACE_PHIN_EN_" << idx;
  Nom nom_phin(oss2.str().c_str());
  if ((liste_post_instantanes.contient_("INTERFACE_TEMPERATURE")) || (liste_post_instantanes.contient_("INTERFACE_PHIN")) || (liste_post_instantanes.contient_(nom_temp))
      || (liste_post_instantanes.contient_(nom_phin)))
    {
      //  Computing interfacial temperature at fa7 centre :
      const Maillage_FT_IJK& mesh = interfaces_.maillage_ft_ijk();
      const ArrOfDouble& surface_facettes = mesh.get_update_surface_facettes();
      const int nb_facettes = mesh.nb_facettes();
      ArrOfDouble interfacial_temperature;
      ArrOfDouble interfacial_phin;
      // To transfer the field to FT splitting (because interfaces are there...) !!! NEEDED for compute_interfacial_temperature
      IJK_Field_double& temperature_ft = itr.get_temperature_ft();
      ref_ijk_ft_.redistribute_to_splitting_ft_elem_.redistribute(itr.get_temperature(), temperature_ft);
      temperature_ft.echange_espace_virtuel(temperature_ft.ghost());
      // results are prop to the area :
      //itr.compute_interfacial_temperature(interfacial_temperature, interfacial_phin, itr.get_storage());
      itr.compute_interfacial_temperature2(interfacial_temperature, interfacial_phin);
      for (int fa7 = 0; fa7 < nb_facettes; fa7++)
        {
          const double sf = surface_facettes[fa7];
          interfacial_temperature[fa7] /= sf;
          interfacial_phin[fa7] /= sf;
        }
      if ((liste_post_instantanes.contient_("INTERFACE_TEMPERATURE")) || (liste_post_instantanes.contient_(nom_temp)))
        n++, dumplata_ft_field(lata_name, "INTERFACES", nom_temp, "ELEM", interfacial_temperature, latastep);
      if ((liste_post_instantanes.contient_("INTERFACE_PHIN")) || (liste_post_instantanes.contient_(nom_phin)))
        n++, dumplata_ft_field(lata_name, "INTERFACES", nom_phin, "ELEM", interfacial_phin, latastep);
    }
  oss.str("");
  return n;
}
