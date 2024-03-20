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
// File      : IJK_FT.cpp
// Directory : $IJK_ROOT/src/FT
//
/////////////////////////////////////////////////////////////////////////////

#include <IJK_FT.h>
#include <Param.h>
#include <Interprete_bloc.h>
#include <EFichier.h>
#include <SFichier.h>
#include <stat_counters.h>
#include <IJK_Lata_writer.h>
#include <IJK_Navier_Stokes_tools.h>
#include <communications.h>
#include <LecFicDiffuse_JDD.h>
#include <MaillerParallel.h>
#include <EChaine.h>
#include <SChaine.h>
#include <Probleme_base.h>
#include <Domaine_VF.h>
#include <Sonde_IJK.h>
#include <Ouvrir_fichier.h>
#include <EcritureLectureSpecial.h>
#include <Init_spectral.h>
#include <init_forcage_THI.h>
#include <corrections_qdm.h>
#include <Force_sp.h>
#include <Force_ph.h>
#include <Vecteur3.h>

#define COMPLEMENT_ANTI_DEVIATION_RESIDU
// #define VARIABLE_DZ
// #define PROJECTION_DE_LINCREMENT_DV
// #define SMOOTHING_RHO
// static Stat_Counter_Id cnt_SourceInterf;

#define select(a,x,y,z) ((a==0)?(x):((a==1)?(y):(z)))

//#define SMOOTHING_RHO

Implemente_instanciable_sans_constructeur(IJK_FT_double, "IJK_FT_double", Interprete);

IJK_FT_double::IJK_FT_double():
  post_(IJK_FT_Post(*this)),
  thermals_(IJK_Thermals(*this))
{
  p_seuil_min_ = 0.;
  p_seuil_max_ = 0;

  vitesse_entree_ = 0.;
  vitesse_upstream_ = 0.;
  expression_vitesse_upstream_ = "??";
  upstream_dir_ = 0;
  upstream_stencil_ = 0;
  nb_diam_upstream_ = 0.;

  rho_liquide_ = 0.;
  rho_vapeur_ = 0.;
  rho_moyen_ = 0.;
  mu_liquide_ = 0.;
  mu_vapeur_ = 0.;
  sigma_ = 0.;

  disable_solveur_poisson_ = 0;
  disable_diffusion_qdm_ = 0;
  disable_convection_qdm_ = 0;
  disable_source_interf_ = 0;
  disable_diphasique_ = 0;
  frozen_velocity_ = 0;
  resolution_fluctuations_ = 0;
  velocity_reset_ = 0;

  improved_initial_pressure_guess_ = 0;
  include_pressure_gradient_in_ustar_ = 0;

  use_inv_rho_for_mass_solver_and_calculer_rho_v_ = 0;
  use_inv_rho_in_poisson_solver_ = 0;
  use_inv_rho_ = 0;

  correction_bilan_qdm_ = 0;
  refuse_patch_conservation_QdM_RK3_source_interf_ = 0;
  test_etapes_et_bilan_ = 0;
  add_initial_field_ = 0;

  diffusion_alternative_ = 0;
  suppression_rejetons_ = 0;

  time_scheme_ = 0;
  store_RK3_source_acc_ = 0.;
  store_RK3_fac_sv_ = 0.;
  rk_step_ = 0;

  modified_time_ini_ = 0.;
  current_time_ = 0.;
  current_time_at_rk3_step_ = 0.;
  tstep_ = 0;
  timestep_ = 0.;
  max_simu_time_ = (int) 1e6;
  dt_sauvegarde_ = 0;
  nb_timesteps_ = 0;
  timestep_reprise_vitesse_ = 0;

  timestep_facsec_ = 0.;
  cfl_ = 0.;
  fo_ = 0.;
  oh_ = 0.;
  enable_dt_oh_ideal_length_factor_ = 0;

  ijk_splitting_ft_extension_ = 0;
  thermal_probes_ghost_cells_ = 2;

  pression_ap_proj_ = 0.;
  check_divergence_ = 0;
  coef_force_time_n_ = 0;
  facteur_variable_source_ = 0;
  coef_ammortissement_ = 0;
  coef_rayon_force_rappel_ = 0;
  terme_source_acceleration_ = 0;
  coef_mean_force_ = 0;
  coef_immobilisation_ = 0;

  vol_bulle_monodisperse_ = 0;

  compute_force_init_ = 0;
  direction_gravite_ = 0;

  sauvegarder_xyz_ = 0;
  reprise_ = 0;
}

IJK_FT_double::IJK_FT_double(const IJK_FT_double& x):
  Interprete(x), post_(IJK_FT_Post(*this))
{
  exit();
}

// Fonctions mutualisees avec IJK_switch.
// Pour cela, deplacement vers IJK_Navier_Stokes_tools.cpp.P
// static void extend_array(const IJK_Grid_Geometry &geom1, ...
// void build_extended_splitting(const IJK_Splitting &split1, ...
// Probleme_FT_Disc_gen & creer_domaine_vdf(const IJK_Splitting & splitting, const Nom & nom_domaine)

#ifdef SMOOTHING_RHO
static void smoothing_field(IJK_Field_double& field,
                            const double rhol, const double  rhov,
                            const double ratio_density_max)
{
  const int nx = field.ni();
  const int ny = field.nj();
  const int nz = field.nk();
  const int ng= field.ghost();
  IJK_Field_local_double smooth_rho;
  smooth_rho.allocate(nx, ny, nz, ng);
  double ratio=std::max(rhol/rhov, rhov/rhol);
  int nb_loop_smoothing = 0;
  while (ratio>ratio_density_max)
    {
      nb_loop_smoothing++;
      ratio /=ratio_density_max;
    }
  Cerr << "Density smoothing activated... nb_loops : "
       << nb_loop_smoothing << " -- Ratio density max : "
       << ratio_density_max << finl;


  for (int iloop = 0; iloop <nb_loop_smoothing; iloop++)
    {
      int ncells = 0;
      for (int k=0; k < nz ; k++)
        for (int j=0; j< ny; j++)
          for (int i=0; i < nx; i++)
            {
              double x = field(i,j,k);
              double xmax = x;
              double xmin = x;
              for (int idx=-1; idx < 2; idx+=2)
                {
                  double val;
                  val = field(i+idx,j,k);
                  xmax=std::max(xmax,val);
                  xmin=std::min(xmin,val);
                  val = field(i,j+idx,k);
                  xmax=std::max(xmax,val);
                  xmin=std::min(xmin,val);
                  val = field(i,j,k+idx);
                  xmax=std::max(xmax,val);
                  xmin=std::min(xmin,val);
                }
              if (xmax/x > ratio_density_max || x/xmin > ratio_density_max )
                {
                  ncells++;
                  //	      Cerr << "   i,j,k,xmin,x,xmax " << i << " " << j << " "
                  //		   << k << " " << xmin  << " " << x << " " << xmax << finl;
                  smooth_rho(i,j,k) = sqrt(xmin*xmax);// prendre la moyenne geometrique du min et du max des voisins.
                }
              else
                {
                  smooth_rho(i,j,k) = x;
                }
            }
      // On utilise le smoothed_rho :
      for (int k=0; k < nz ; k++)
        for (int j=0; j< ny; j++)
          for (int i=0; i < nx; i++)
            field(i,j,k) =  smooth_rho(i,j,k);
      field.echange_espace_virtuel(ng);
      Cerr << "Counted cells : " << ncells << " -- Fin loop " << iloop << finl;
    }

  {
    double mmax = 1.;
    for (int k=0; k < nz ; k++)
      for (int j=0; j< ny; j++)
        for (int i=0; i < nx; i++)
          {
            double xmin = field(i,j,k);
            double xmax = xmin;
            double x = field(i,j,k);
            for (int idx=-1; idx < 2; idx+=2)
              {
                double val;
                val = field(i+idx,j,k);
                xmax=std::max(xmax,val);
                xmin=std::min(xmin,val);
                val = field(i,j+idx,k);
                xmax=std::max(xmax,val);
                xmin=std::min(xmin,val);
                val = field(i,j,k+idx);
                xmax=std::max(xmax,val);
                xmin=std::min(xmin,val);
              }
            mmax=std::max(mmax,std::max(xmax/x,x/xmin));
          }
    Cerr << "Ratio max sur le rho filtre : " << mmax << finl;
  }

}
#endif

#if 0
// Codage optimise de copy_field_to_extended_domain


// Attention: algorithme ne supporte pas une extension plus grande
// que la taille d'un sous-domaine. Ne supporte pas non plus un process_grouping
// different de 1,1,1 dans le splitting.
void copy_field_to_extended_domain(const IJK_Field_double& input_field,
                                   IJK_Field_double& output_field,
                                   const int nextension)
{
  // Algorithme vite ecrit mais lent pour redistribuer:
  // Creation d'une copie de input field avec nextension cellules ghost,
  // ensuite, on a localement toutes les donnees dont on a besoin
  IJK_Field_double input_field2;
  input_field2.allocate(input_field.get_splitting(), // sur decoupage initial
                        input_field.get_localisation(), // meme localisation
                        nextension+1); // mais avec nextension ghost cells
  // copie input_field dans input field2:
  int i, j, k;
  const int ni = input_field.ni();
  const int nj = input_field.nj();
  const int nk = input_field.nk();
  for (k = 0; k < nk; k++)
    for (j = 0; j < nj; j++)
      for (i = 0; i < ni; i++)
        {
          input_field2(i,j,k) = input_field(i,j,k);
        }
  // Echange espace virtuel pour avoir les valeurs sur nextension mailles
  // autour du domaine local initial.
  // On espere que le domaine etendu sur ce processeur est a l'interieur...
  // Il nous faut une epaisseur supplementaire car les faces de droite du
  // maillage etendu portent une inconnue qui n'existe pas dans le maillage
  // d'origine, s'il est periodique.
  input_field2.echange_espace_virtuel(nextension+1);
  // On va recopier les donnees de input_field2 dans output_field, en
  // decalant comme il faut.

  // Taille du maillage destination:
  const int ni2 = output_field.ni();
  const int nj2 = output_field.nj();
  const int nk2 = output_field.nk();
  const IJK_Grid_Geometry& geom = input_field.get_splitting().get_grid_geometry();
  const int di = (geom.get_periodic_flag(DIRECTION_I) ? (-nextension) : 0)
                 - input_field.get_splitting().get_offset_local(DIRECTION_I)
                 + output_field.get_splitting().get_offset_local(DIRECTION_I);
  const int dj = (geom.get_periodic_flag(DIRECTION_J) ? (-nextension) : 0)
                 - input_field.get_splitting().get_offset_local(DIRECTION_J)
                 + output_field.get_splitting().get_offset_local(DIRECTION_J);
  const int dk = (geom.get_periodic_flag(DIRECTION_K) ? (-nextension) : 0)
                 - input_field.get_splitting().get_offset_local(DIRECTION_K)
                 + output_field.get_splitting().get_offset_local(DIRECTION_K);

  for (k = 0; k < nk2; k++)
    for (j = 0; j < nj2; j++)
      for (i = 0; i < ni2; i++)
        {
          output_field(i, j, k) = input_field2(i+di, j+dj, k+dk);
        }
  // L'espace virtuel de output_field doit etre rempli.
  // C'est delicat sur les ex-bords periodiques qui ne le sont plus.
  const int ghost = output_field.ghost();
  output_field.echange_espace_virtuel(ghost);
  // L'espace virtuel est a jour. Il est remplit correctement si les
  // flags periodiques sur l'extended splitting sont a faux.
}


void add_field_from_extended_domain(const IJK_Field_double& input_field,
                                    IJK_Field_double& output_field,
                                    const int nextension)
{
  if (nextension > 0 && Process::nproc() > 1)
    {
      Cerr << "Erreur : pas encore code (parallele + nextension>0)" << finl;
      Process::exit();
    }
  for (int k = 0; k < output_field.nk(); k++)
    for (int j = 0; j < output_field.nj(); j++)
      for (int i = 0; i < output_field.ni(); i++)
        output_field(i,j,k) += input_field(i + nextension, j + nextension, k);
}

void copy_field_from_extended_domain(const IJK_Field_double& input_field,
                                     IJK_Field_double& output_field,
                                     const int nextension)
{
  if (nextension > 0 && Process::nproc() > 1)
    {
      Cerr << "Erreur : pas encore code (parallele + nextension>0)" << finl;
      Process::exit();
    }
  for (int k = 0; k < output_field.nk(); k++)
    for (int j = 0; j < output_field.nj(); j++)
      for (int i = 0; i < output_field.ni(); i++)
        output_field(i,j,k) = input_field(i + nextension, j + nextension, k);
}

#endif

IJK_FT_double::TimeScheme IJK_FT_double::get_time_scheme() const
{
  return (TimeScheme) time_scheme_;
}

// XD IJK_FT_double interprete IJK_FT_double 1 not_set
Entree& IJK_FT_double::interpreter(Entree& is)
{
  tstep_ = 0;

  // On force l'attribut dimension a 3 pour ne pas avoir besoin de le mettre dans le jeu de donnees.
  // Cet attribut est utilise dans les routines front-tracking issues de triou
  Objet_U::dimension=3;

  check_divergence_ = 0;
  rk_step_ = -1; // default value

  expression_pression_initiale_ = "??"; // par defaut, invalide
  fichier_reprise_vitesse_ = "??"; // par defaut, invalide
  Param param(que_suis_je());
  Nom ijk_splitting_name;

  dt_sauvegarde_ = 2000000000; // jamais
  current_time_ = 0.;
  nom_sauvegarde_ = nom_du_cas() + ".sauv";
  gravite_.resize_array(3);
  gravite_ = 0.;

  // GAB, rotation
  direction_gravite_ = 0;
  //
  // GAB, qdm
  // terme_diffusion.resize_array(3);
  // terme_convection.resize_array(3);
  // terme_pression.resize_array(3);
  // rho_u_euler_av_prediction.resize_array(3);
  // rho_u_euler_av_projection.resize_array(3);
  // rho_u_euler_ap_prediction.resize_array(3);
  // rho_u_euler_ap_projection.resize_array(3);
  rho_u_euler_av_prediction_ = 0.;
  rho_du_euler_ap_prediction_ = 0.;
  rho_u_euler_ap_projection_ = 0.;
  rho_du_euler_ap_projection_ = 0.;
  rho_u_euler_av_rho_mu_ind_ = 0.;
  rho_u_euler_ap_rho_mu_ind_ = 0.; //7.;
  u_euler_ap_rho_mu_ind_ = 0.;
  terme_diffusion_ = 0.;
  terme_convection_ = 0.;
  terme_pression_ = 0.;
  terme_pression_bis_ = 0.;
  terme_pression_ter_ = 0.;
  terme_interfaces_ = 0.;
  terme_interfaces_bf_mass_solver_ = 0.;
  terme_interfaces_bf_mass_solver_bis_ = 0.;
  terme_interfaces_af_mass_solver_ = 0.;
  terme_interfaces_conv_diff_mass_solver_ = 0.;
  pression_ap_proj_ = 0.;
  terme_moyen_convection_mass_solver_ = 0.;
  terme_moyen_diffusion_mass_solver_ = 0.;
  //
  vitesse_entree_ = -1.1e20;
  vitesse_upstream_ = -1.1e20;
  expression_vitesse_upstream_ = "??";
  nb_diam_upstream_ = 0.;
  upstream_dir_=-1;
  upstream_stencil_=3;
  disable_solveur_poisson_ = 0;
  resolution_fluctuations_ = 0;
  disable_diffusion_qdm_ = 0;
  disable_convection_qdm_ = 0;
  disable_source_interf_ = 0;
  frozen_velocity_ = 0;
  velocity_reset_ = 0;
  disable_diphasique_ = 0;
  improved_initial_pressure_guess_ = 0;
  include_pressure_gradient_in_ustar_ = 0;
  use_inv_rho_ = 0;
  use_inv_rho_for_mass_solver_and_calculer_rho_v_ = 0;
  use_inv_rho_in_poisson_solver_ = 0;
  correction_bilan_qdm_ = 0;
  diffusion_alternative_ = 0;
  suppression_rejetons_ = 0; // By defaults, break-ups are not fixed on restart. (no deletion of smaller fractions)
  refuse_patch_conservation_QdM_RK3_source_interf_ = 0; // Par defaut, on utilise le patch!
  // GAB, qdm
  test_etapes_et_bilan_ = 0;
  //
  time_scheme_ = EULER_EXPLICITE;
  sauvegarder_xyz_ = 0;

  reprise_ = 0; // Indique si on fait une reprise ou pas.

  timestep_facsec_ = 1.;
  cfl_ = 1.;
  fo_ = 1.;
  oh_ = 1.;

  rho_vapeur_ = -1.;
  mu_vapeur_ = -1.;
  sigma_     = 0.;

  //ab-forcage-control-ecoulement-deb
  expression_derivee_acceleration_ = "0"; // par defaut pas de terme d'acceleration
  terme_source_acceleration_ = 0.; // par defaut, zero
  integrated_residu_ = 0.;
  //ab-forcage-control-ecoulement-fin

  expression_potential_phi_ = "??";
  fichier_post_ = "??";

  terme_source_correction_.resize_array(3); // Initialement a zero, puis sera calcule a chaque iter.
  terme_source_correction_ = 0.;
  //facteur_variable_source_.resize_array(3); // Initialement a 1, puis sera calcule a chaque iter si expression est donnee.
  facteur_variable_source_= 1.;
  expression_derivee_facteur_variable_source_ = "0";
  correction_force_.resize_array(3); // Par defaut, les flags d'activations sont a zero (ie inactif).
  correction_force_ = 0;

  gravite_.resize_array(3);
  // GAB, rotation
  direction_gravite_ = 0;
  //
  gravite_ = 0.;
  vol_bulle_monodisperse_ = -1.; // Le volume des bulles n'est pas impose par defaut.
  vol_bulles_.resize_array(0); // Initialement a zero, puis sera calcule a chaque iter.
  vol_bulles_ = 0.;
  store_RK3_source_acc_ = 0.;
  store_RK3_fac_sv_ = 1.;
#ifdef SMOOTHING_RHO
  smooth_density_ = 0;
  ratio_density_max_ = 15;
  param.ajouter_flag("smooth_density", &smooth_density_);
  param.ajouter("ratio_density_max", &ratio_density_max_);
#endif

  // valeurs par default des parametres de bulles fixes
  coef_immobilisation_ = 0.;
  coef_ammortissement_ = 0.;
  coef_mean_force_=0.;
  coef_force_time_n_=0.;
  coef_rayon_force_rappel_ = 0.;
  p_seuil_max_ = 10000000 ;
  p_seuil_min_ = -10000000 ;
  param.ajouter("p_seuil_max", &p_seuil_max_); // XD_ADD_P floattant not_set, default 10000000
  param.ajouter("p_seuil_min", &p_seuil_min_); // XD_ADD_P floattant not_set, default -10000000
  param.ajouter("coef_ammortissement", &coef_ammortissement_); // XD_ADD_P floattant not_set
  param.ajouter("coef_immobilisation", &coef_immobilisation_); // XD_ADD_P floattant not_set
  param.ajouter("coef_mean_force", &coef_mean_force_); // XD_ADD_P floattant not_set
  param.ajouter("coef_force_time_n", &coef_force_time_n_); // XD_ADD_P floattant not_set
  param.ajouter("coef_rayon_force_rappel", &coef_rayon_force_rappel_); // XD_ADD_P floattant not_set
  param.ajouter("tinit", &current_time_); // XD_ADD_P floattant initial time
  param.ajouter("ijk_splitting", &ijk_splitting_name, Param::REQUIRED); // XD_ADD_P chaine(into=["grid_splitting"]) Definition of domain decomposition for parallel computations
  param.ajouter("timestep", &timestep_, Param::REQUIRED); // XD_ADD_P floattant Upper limit of the timestep
  param.ajouter("timestep_facsec", &timestep_facsec_); // XD_ADD_P floattant Security factor on timestep
  param.ajouter("cfl", &cfl_); // XD_ADD_P floattant  To provide a value of the limiting CFL number used for setting the timestep
  param.ajouter("fo", &fo_); // XD_ADD_P floattant not_set
  param.ajouter("oh", &oh_); // XD_ADD_P floattant not_set
  param.ajouter_flag("enable_dt_oh_ideal_length_factor", &enable_dt_oh_ideal_length_factor_);
  param.ajouter("nb_pas_dt_max", &nb_timesteps_, Param::REQUIRED); // XD_ADD_P entier maximum limit for the number of timesteps
  param.ajouter("max_simu_time", &max_simu_time_); // XD_ADD_P entier maximum limit for the number of timesteps
  param.ajouter("multigrid_solver", &poisson_solver_, Param::REQUIRED); // XD_ADD_P multigrid_solver not_set
  param.ajouter_flag("check_divergence", &check_divergence_); // XD_ADD_P rien Flag to compute and print the value of div(u) after each pressure-correction
  param.ajouter("mu_liquide", &mu_liquide_, Param::REQUIRED); // XD_ADD_P floattant liquid viscosity
  param.ajouter("vitesse_entree", &vitesse_entree_); // XD_ADD_P floattant Velocity to prescribe at inlet
  param.ajouter("vitesse_upstream", &vitesse_upstream_); // XD_ADD_P floattant Velocity to prescribe at 'nb_diam_upstream_' before bubble 0.
  param.ajouter("upstream_dir", &upstream_dir_); // XD_ADD_P entier Direction to prescribe the velocity
  param.ajouter("expression_vitesse_upstream", &expression_vitesse_upstream_); // XD_ADD_P chaine Analytical expression to set the upstream velocity
  param.ajouter("upstream_stencil", &upstream_stencil_); // XD_ADD_P int Width on which the velocity is set
  param.ajouter("nb_diam_upstream", &nb_diam_upstream_); // XD_ADD_P floattant Number of bubble diameters upstream of bubble 0 to prescribe the velocity.
  param.ajouter("rho_liquide", &rho_liquide_, Param::REQUIRED); // XD_ADD_P floattant liquid density
  param.ajouter("check_stop_file", &check_stop_file_); // XD_ADD_P chaine stop file to check (if 1 inside this file, stop computation)
  param.ajouter("dt_sauvegarde", &dt_sauvegarde_); // XD_ADD_P entier saving frequency (writing files for computation restart)
  param.ajouter("nom_sauvegarde", &nom_sauvegarde_); // XD_ADD_P chaine Definition of filename to save the calculation
  param.ajouter_flag("sauvegarder_xyz", &sauvegarder_xyz_); // XD_ADD_P rien save in xyz format
  param.ajouter("nom_reprise", &nom_reprise_); // XD_ADD_P chaine Enable restart from filename given
  param.ajouter("gravite", &gravite_); // XD_ADD_P list gravity vector [gx, gy, gz]
  expression_vitesse_initiale_.dimensionner_force(3);
  param.ajouter("expression_vx_init", &expression_vitesse_initiale_[0]); // XD_ADD_P chaine initial field for x-velocity component (parser of x,y,z)
  param.ajouter("expression_vy_init", &expression_vitesse_initiale_[1]); // XD_ADD_P chaine initial field for y-velocity component (parser of x,y,z)
  param.ajouter("expression_vz_init", &expression_vitesse_initiale_[2]); // XD_ADD_P chaine initial field for z-velocity component (parser of x,y,z)
  param.ajouter("expression_derivee_force", &expression_derivee_acceleration_); // XD_ADD_P chaine expression of the time-derivative of the X-component of a source-term (see terme_force_ini for the initial value). terme_force_ini : initial value of the X-component of the source term (see expression_derivee_force  for time evolution)
  param.ajouter_flag("compute_force_init", &compute_force_init_); // XD_ADD_P chaine not_set
  param.ajouter("terme_force_init", &terme_source_acceleration_); // XD_ADD_P chaine not_set
  param.ajouter("correction_force", &correction_force_); // XD_ADD_P chaine not_set
  param.ajouter("vol_bulle_monodisperse", &vol_bulle_monodisperse_); // XD_ADD_P chaine not_set
  param.ajouter("vol_bulles", &vol_bulles_); // XD_ADD_P chaine not_set
  param.ajouter("time_scheme", &time_scheme_); // XD_ADD_P chaine(into=["euler_explicit","RK3_FT"]) Type of time scheme
  param.dictionnaire("euler_explicit", EULER_EXPLICITE);
  param.dictionnaire("RK3_FT", RK3_FT);

  // GAB question : pourquoi expression_variable_source_ est de type nom et pas de type Vecteur3 ??
  expression_variable_source_.dimensionner_force(3);
  param.ajouter("expression_variable_source_x", &expression_variable_source_[0]); // XD_ADD_P chaine not_set
  param.ajouter("expression_variable_source_y", &expression_variable_source_[1]); // XD_ADD_P chaine not_set
  param.ajouter("expression_variable_source_z", &expression_variable_source_[2]); // XD_ADD_P chaine not_set
  param.ajouter("facteur_variable_source_init", &facteur_variable_source_); // XD_ADD_P chaine not_set
  param.ajouter("expression_derivee_facteur_variable_source", &expression_derivee_facteur_variable_source_); // XD_ADD_P chaine not_set

  param.ajouter("expression_p_init", &expression_pression_initiale_); // XD_ADD_P chaine initial pressure field (optional)

  param.ajouter("expression_potential_phi", &expression_potential_phi_); // XD_ADD_P chaine parser to define phi and make a momentum source Nabla phi.

  param.ajouter("velocity_diffusion_op", &velocity_diffusion_op_);
  param.ajouter("velocity_convection_op", &velocity_convection_op_);

  param.ajouter("interfaces", &interfaces_); // XD_ADD_P interfaces not_set
  // GAB, THI
  param.ajouter("forcage", &forcage_);  // XD_ADD_P chaine not_set
  param.ajouter("corrections_qdm", &qdm_corrections_); // XD_ADD_P chaine not_set

  // Read list of thermic equations:
  param.ajouter("thermals", &thermals_);

  /*
   * TODO: Change this block with DERIV CLASS IJK_Thermal
   */
  param.ajouter("thermique", &thermique_); // XD_ADD_P thermique not_set
  param.ajouter("energie", &energie_); // XD_ADD_P chaine not_set


  param.ajouter("ijk_splitting_ft_extension", &ijk_splitting_ft_extension_, Param::REQUIRED); // XD_ADD_P entier Number of element used to extend the computational domain at each side of periodic boundary to accommodate for bubble evolution.

  param.ajouter("fichier_post", &fichier_post_); // XD_ADD_P chaine name of the post-processing file (lata file)
  // ATTENTION les fichiers reprises sont des fichiers .lata ou sauv.lata
  // On peut reprendre uniquement la vitesse ou uniquement rho dans un fichier de post:
  param.ajouter("fichier_reprise_vitesse", &fichier_reprise_vitesse_); // XD_ADD_P chaine not_set
  param.ajouter("timestep_reprise_vitesse", &timestep_reprise_vitesse_); // XD_ADD_P chaine not_set

  param.ajouter("boundary_conditions", &boundary_conditions_, Param::REQUIRED); // XD_ADD_P bloc_lecture BC
  param.ajouter_flag("disable_solveur_poisson", &disable_solveur_poisson_); // XD_ADD_P rien Disable pressure poisson solver
  param.ajouter_flag("resolution_fluctuations", &resolution_fluctuations_); // XD_ADD_P rien Disable pressure poisson solver
  param.ajouter_flag("disable_diffusion_qdm", &disable_diffusion_qdm_); // XD_ADD_P rien Disable diffusion operator in momentum
  param.ajouter_flag("disable_source_interf", &disable_source_interf_); // XD_ADD_P rien Disable computation of the interfacial source term
  param.ajouter_flag("disable_convection_qdm", &disable_convection_qdm_); // XD_ADD_P rien Disable convection operator in momentum
  param.ajouter_flag("disable_diphasique", &disable_diphasique_); // XD_ADD_P rien Disable all calculations related to interfaces (phase properties, interfacial force, ... )
  param.ajouter_flag("frozen_velocity", &frozen_velocity_); // XD_ADD_P chaine not_set
  param.ajouter_flag("velocity_reset", &velocity_reset_); // XD_ADD_P chaine not_set
  param.ajouter_flag("improved_initial_pressure_guess", &improved_initial_pressure_guess_); // XD_ADD_P chaine not_set
  param.ajouter_flag("include_pressure_gradient_in_ustar", &include_pressure_gradient_in_ustar_); // XD_ADD_P chaine not_set
  // param.ajouter_flag("use_inv_rho", &use_inv_rho_);
  param.ajouter_flag("use_inv_rho_for_mass_solver_and_calculer_rho_v", &use_inv_rho_for_mass_solver_and_calculer_rho_v_); // XD_ADD_P chaine not_set
  param.ajouter_flag("use_inv_rho_in_poisson_solver", &use_inv_rho_in_poisson_solver_); // XD_ADD_P chaine not_set
  param.ajouter_flag("diffusion_alternative", &diffusion_alternative_); // XD_ADD_P chaine not_set
  param.ajouter_flag("suppression_rejetons", &suppression_rejetons_); // XD_ADD_P chaine not_set
  // param.ajouter("correction_bilan_qdm", &correction_bilan_qdm_); // XD_ADD_P chaine not_set
  param.ajouter_flag("refuse_patch_conservation_QdM_RK3_source_interf", &refuse_patch_conservation_QdM_RK3_source_interf_); // XD_ADD_P rien experimental Keyword, not for use
  // GAB; qdm
  param.ajouter_flag("test_etapes_et_bilan", &test_etapes_et_bilan_); // XD_ADD_P chaine not_set
  // GAB, champ de reprise + champ initial
  param.ajouter_flag("ajout_init_a_reprise", &add_initial_field_); // XD_ADD_P chaine not_set


  param.ajouter("reprise_vap_velocity_tmoy", &vap_velocity_tmoy_); // XD_ADD_P chaine not_set
  param.ajouter("reprise_liq_velocity_tmoy", &liq_velocity_tmoy_); // XD_ADD_P chaine not_set
  vap_velocity_tmoy_ = reprise_vap_velocity_tmoy_;
  liq_velocity_tmoy_ = reprise_liq_velocity_tmoy_;


  param.ajouter("sigma", &sigma_); // XD_ADD_P floattant surface tension
  param.ajouter("rho_vapeur", &rho_vapeur_); // XD_ADD_P floattant vapour density
  param.ajouter("mu_vapeur", &mu_vapeur_); // XD_ADD_P floattant vapour viscosity

  param.ajouter_flag("first_step_interface_smoothing", &first_step_interface_smoothing_);
  post_.complete_interpreter(param, is);

// XD attr check_stats rien check_stats 1 Flag to compute additional (xy)-plane averaged statistics
// XD attr dt_post entier dt_post 1 Post-processing frequency (for lata output)
// XD attr dt_post_stats_plans entier dt_post_stats_plans 1 Post-processing frequency for averaged statistical files (txt files containing averaged information on (xy) planes for each z-center) both instantaneous, or cumulated time-integration (see file header for variables list)
// XD attr dt_post_stats_bulles entier dt_post_stats_bulles 1 Post-processing frequency for bubble information (for out files as bubble area, centroid position, etc...)
// XD attr champs_a_postraiter listchaine champs_a_postraiter 1 List of variables to post-process in lata files.
// XD attr expression_vx_ana chaine expression_vx_ana 1 Analytical Vx (parser of x,y,z, t) used for post-processing only
// XD attr expression_vy_ana chaine expression_vy_ana 1 Analytical Vy (parser of x,y,z, t) used for post-processing only
// XD attr expression_vz_ana chaine expression_vz_ana 1 Analytical Vz (parser of x,y,z, t) used for post-processing only
// XD attr expression_p_ana chaine expression_p_ana 1 analytical pressure solution (parser of x,y,z, t) used for post-processing only
// XD attr expression_dPdx_ana chaine expression_dPdx_ana 1 analytical expression dP/dx=f(x,y,z,t), for post-processing only
// XD attr expression_dPdy_ana chaine expression_dPdy_ana 1 analytical expression dP/dy=f(x,y,z,t), for post-processing only
// XD attr expression_dPdz_ana chaine expression_dPdz_ana 1 analytical expression dP/dz=f(x,y,z,t), for post-processing only
// XD attr expression_dUdx_ana chaine expression_dUdx_ana 1 analytical expression dU/dx=f(x,y,z,t), for post-processing only
// XD attr expression_dUdy_ana chaine expression_dUdy_ana 1 analytical expression dU/dy=f(x,y,z,t), for post-processing only
// XD attr expression_dUdz_ana chaine expression_dUdz_ana 1 analytical expression dU/dz=f(x,y,z,t), for post-processing only
// XD attr expression_dVdx_ana chaine expression_dVdx_ana 1 analytical expression dV/dx=f(x,y,z,t), for post-processing only
// XD attr expression_dVdy_ana chaine expression_dVdy_ana 1 analytical expression dV/dy=f(x,y,z,t), for post-processing only
// XD attr expression_dVdz_ana chaine expression_dVdz_ana 1 analytical expression dV/dz=f(x,y,z,t), for post-processing only
// XD attr expression_dWdx_ana chaine expression_dWdx_ana 1 analytical expression dW/dx=f(x,y,z,t), for post-processing only
// XD attr expression_dWdy_ana chaine expression_dWdy_ana 1 analytical expression dW/dy=f(x,y,z,t), for post-processing only
// XD attr expression_dWdz_ana chaine expression_dWdz_ana 1 analytical expression dW/dz=f(x,y,z,t), for post-processing only
// XD attr expression_ddPdxdx_ana chaine expression_ddPdxdx_ana 1 analytical expression d2P/dx2=f(x,y,z,t), for post-processing only
// XD attr expression_ddPdydy_ana chaine expression_ddPdydy_ana 1 analytical expression d2P/dy2=f(x,y,z,t), for post-processing only
// XD attr expression_ddPdzdz_ana chaine expression_ddPdzdz_ana 1 analytical expression d2P/dz2=f(x,y,z,t), for post-processing only
// XD attr expression_ddPdxdy_ana chaine expression_ddPdxdy_ana 1 analytical expression d2P/dxdy=f(x,y,z,t), for post-processing only
// XD attr expression_ddPdxdz_ana chaine expression_ddPdxdz_ana 1 analytical expression d2P/dxdz=f(x,y,z,t), for post-processing only
// XD attr expression_ddPdydz_ana chaine expression_ddPdydz_ana 1 analytical expression d2P/dydz=f(x,y,z,t), for post-processing only
// XD attr expression_ddUdxdx_ana chaine expression_ddUdxdx_ana 1 analytical expression d2U/dx2=f(x,y,z,t), for post-processing only
// XD attr expression_ddUdydy_ana chaine expression_ddUdydy_ana 1 analytical expression d2U/dy2=f(x,y,z,t), for post-processing only
// XD attr expression_ddUdzdz_ana chaine expression_ddUdzdz_ana 1 analytical expression d2U/dz2=f(x,y,z,t), for post-processing only
// XD attr expression_ddUdxdy_ana chaine expression_ddUdxdy_ana 1 analytical expression d2U/dxdy=f(x,y,z,t), for post-processing only
// XD attr expression_ddUdxdz_ana chaine expression_ddUdxdz_ana 1 analytical expression d2U/dxdz=f(x,y,z,t), for post-processing only
// XD attr expression_ddUdydz_ana chaine expression_ddUdydz_ana 1 analytical expression d2U/dydz=f(x,y,z,t), for post-processing only
// XD attr expression_ddVdxdx_ana chaine expression_ddVdxdx_ana 1 analytical expression d2V/dx2=f(x,y,z,t), for post-processing only
// XD attr expression_ddVdydy_ana chaine expression_ddVdydy_ana 1 analytical expression d2V/dy2=f(x,y,z,t), for post-processing only
// XD attr expression_ddVdzdz_ana chaine expression_ddVdzdz_ana 1 analytical expression d2V/dz2=f(x,y,z,t), for post-processing only
// XD attr expression_ddVdxdy_ana chaine expression_ddVdxdy_ana 1 analytical expression d2V/dxdy=f(x,y,z,t), for post-processing only
// XD attr expression_ddVdxdz_ana chaine expression_ddVdxdz_ana 1 analytical expression d2V/dxdz=f(x,y,z,t), for post-processing only
// XD attr expression_ddVdydz_ana chaine expression_ddVdydz_ana 1 analytical expression d2V/dydz=f(x,y,z,t), for post-processing only
// XD attr expression_ddWdxdx_ana chaine expression_ddWdxdx_ana 1 analytical expression d2W/dx2=f(x,y,z,t), for post-processing only
// XD attr expression_ddWdydy_ana chaine expression_ddWdydy_ana 1 analytical expression d2W/dy2=f(x,y,z,t), for post-processing only
// XD attr expression_ddWdzdz_ana chaine expression_ddWdzdz_ana 1 analytical expression d2W/dz2=f(x,y,z,t), for post-processing only
// XD attr expression_ddWdxdy_ana chaine expression_ddWdxdy_ana 1 analytical expression d2W/dxdy=f(x,y,z,t), for post-processing only
// XD attr expression_ddWdxdz_ana chaine expression_ddWdxdz_ana 1 analytical expression d2W/dxdz=f(x,y,z,t), for post-processing only
// XD attr expression_ddWdydz_ana chaine expression_ddWdydz_ana 1 analytical expression d2W/dydz=f(x,y,z,t), for post-processing only
// XD attr t_debut_statistiques floattant t_debut_statistiques 1 Initial time for computation, printing and accumulating time-integration
// XD attr sondes bloc_lecture sondes 1 probes

  param.lire_avec_accolades(is);

  // GAB, rotation
  direction_gravite_ = get_direction(gravite_);
  if (frozen_velocity_)
    {
      disable_solveur_poisson_=1; // automatically force the suppression of the poisson solver
      if (!disable_diphasique_)
        {
          interfaces_.freeze();  // Stop the interfacial displacement.
          Cout << "The option frozen_velocity automatically freeze the interface motion "
               <<  "by activating the flag interfaces_.frozen_" << finl;
        }
    }

  if ((expression_potential_phi_ != "??") &&
      ((expression_variable_source_[0] != "??") ||
       (expression_variable_source_[1] != "??")
       || (expression_variable_source_[2] != "??")))
    {
      Cerr << "expression_potential_phi and expression_variable_source are used together"
           << "nabla(phi) will be added to the expression given for the variable source" << finl;
      //Process::exit();
    }

  if ((include_pressure_gradient_in_ustar_) && (expression_pression_initiale_ == "??"))
    {
      Cerr << "When using pressure increment in u^star, expression_p_init " << expression_pression_initiale_
           << "becomes a required parameter. " << finl;
      Process::exit();
    }

//  if ((correction_bilan_qdm_<0) || (correction_bilan_qdm_>4))
//    {
//      Cerr << "Invalid value of correction_bilan_qdm : " << correction_bilan_qdm_  << ". " << finl;
//      Cerr << "Please use 0, 1, or 2 for no_correction, geometric_mean or arithmetic_mean respectively" << finl;
//      Cerr << "Or 3 to anhilate the residual deviation... " << finl;
//      Cerr << "Or 4 to anhilate the residual deviation (except along z)... " << finl;
//      Process::exit();
//    }

  // Si on utilise un seul groupe et qu'on impose un volume unique a toutes les bulles,
  if (vol_bulle_monodisperse_>=0.)
    {
      if (vol_bulles_.size_array() !=0)
        {
          Cerr << "Attention, conflit entre les options : vol_bulle_monodisperse_ et vol_bulles."
               << "Merci de choisir" << finl;
          Process::exit();
        }
    }

  // On utilise inv_rho pour l'un ou l'autre... Il faut donc le calculer :
  use_inv_rho_ = use_inv_rho_for_mass_solver_and_calculer_rho_v_ + use_inv_rho_in_poisson_solver_;

  if (rho_vapeur_ < 0. )
    rho_vapeur_ = rho_liquide_;
  if (mu_vapeur_ < 0. )
    mu_vapeur_ = mu_liquide_;

  // Avec cette option, on travaille avec nu :
  if (diffusion_alternative_)
    {
      Cerr << "Option diffusion_alternative activee : le champ mu contient nu (la viscosite dynamique)" << finl;
      mu_liquide_ /= rho_liquide_;
      mu_vapeur_ /= rho_vapeur_;
    }

  if (gravite_.size_array() != 3)
    {
      Cerr << "Erreur: la gravite doit etre un vecteur de 3 composantes" << finl;
      exit();
    }

  splitting_ = ref_cast(IJK_Splitting, Interprete_bloc::objet_global(ijk_splitting_name));

  Cerr << "Construction du domaine VDF NS pour les sondes..." << finl;
  refprobleme_ns_ = creer_domaine_vdf(splitting_, "DOM_NS_VDF");

  Cerr << "Construction du domaine VDF..." << finl;
  {
    build_extended_splitting(splitting_, splitting_ft_, ijk_splitting_ft_extension_);
    refprobleme_ft_disc_ = creer_domaine_vdf(splitting_ft_, "DOM_VDF");
    for (int dir = 0; dir < 3; dir++)
      {
        VECT(IntTab) map(3);
        IJK_Splitting::Localisation loc = (dir==0) ? IJK_Splitting::FACES_I : (dir==1) ? IJK_Splitting::FACES_J : IJK_Splitting::FACES_K;
        const int n_ext = ijk_splitting_ft_extension_;

        for (int dir2 = 0; dir2 < 3; dir2++)
          {
            const int n = splitting_.get_nb_items_global(loc, dir2);
            if (n_ext == 0 || !splitting_.get_grid_geometry().get_periodic_flag(dir2))
              {
                map[dir2].resize(1,3);
                map[dir2](0,0) = 0;// source index
                map[dir2](0,1) = 0;// dest index
                map[dir2](0,2) = n;// size
              }
            else
              {
                map[dir2].resize(3,3);
                // copy NS field to central domaine of extended field
                map[dir2](0,0) = 0;
                map[dir2](0,1) = n_ext;
                map[dir2](0,2) = n;
                // copy right part of NS field to left part of extended field
                map[dir2](1,0) = n - n_ext;
                map[dir2](1,1) = 0;
                map[dir2](1,2) = n_ext;
                // copy left part of NS field to right of extended field
                map[dir2](2,0) = 0;
                map[dir2](2,1) = n + n_ext;
                map[dir2](2,2) = n_ext;
              }
          }
        redistribute_to_splitting_ft_faces_[dir].initialize(splitting_, splitting_ft_, loc, map);

        for (int dir2 = 0; dir2 < 3; dir2++)
          {
            const int n = splitting_.get_nb_items_global(loc, dir2);
            if (n_ext == 0 || !splitting_.get_grid_geometry().get_periodic_flag(dir2))
              {
                map[dir2].resize(1,3);
                map[dir2](0,0) = 0;
                map[dir2](0,1) = 0;
                map[dir2](0,2) = n;
              }
            else
              {
                map[dir2].resize(1,3);
                // When copying back from extended splitting, ignore extended data, take only central part
                map[dir2](0,0) = n_ext;  // source index
                map[dir2](0,1) = 0; // dest index
                map[dir2](0,2) = n; // size
              }
          }
        redistribute_from_splitting_ft_faces_[dir].initialize(splitting_ft_, splitting_, loc, map);
      }

    // Pour les elements:
    {
      VECT(IntTab) map(3);
      IJK_Splitting::Localisation loc = IJK_Splitting::ELEM;
      const int n_ext = ijk_splitting_ft_extension_;

      for (int dir2 = 0; dir2 < 3; dir2++)
        {
          const int n = splitting_.get_nb_items_global(loc, dir2);
          if (n_ext == 0 || !splitting_.get_grid_geometry().get_periodic_flag(dir2))
            {
              map[dir2].resize(1,3);
              map[dir2](0,0) = 0;// source index
              map[dir2](0,1) = 0;// dest index
              map[dir2](0,2) = n;// size
            }
          else
            {
              map[dir2].resize(3,3);
              // copy NS field to central domaine of extended field
              map[dir2](0,0) = 0;
              map[dir2](0,1) = n_ext;
              map[dir2](0,2) = n;
              // copy right part of NS field to left part of extended field
              map[dir2](1,0) = n - n_ext;
              map[dir2](1,1) = 0;
              map[dir2](1,2) = n_ext;
              // copy left part of NS field to right of extended field
              map[dir2](2,0) = 0;
              map[dir2](2,1) = n + n_ext;
              map[dir2](2,2) = n_ext;
            }
        }
      redistribute_to_splitting_ft_elem_.initialize(splitting_, splitting_ft_, loc, map);


      for (int dir2 = 0; dir2 < 3; dir2++)
        {
          const int n = splitting_.get_nb_items_global(loc, dir2);
          if (n_ext == 0 || !splitting_.get_grid_geometry().get_periodic_flag(dir2))
            {
              map[dir2].resize(1,3);
              map[dir2](0,0) = 0;
              map[dir2](0,1) = 0;
              map[dir2](0,2) = n;
            }
          else
            {
              map[dir2].resize(1,3);
              // When copying back from extended splitting, ignore extended data, take only central part
              map[dir2](0,0) = n_ext;  // source index
              map[dir2](0,1) = 0; // dest index
              map[dir2](0,2) = n; // size
            }
        }
      redistribute_from_splitting_ft_elem_.initialize(splitting_ft_, splitting_, loc, map);

      for (int dir2 = 0; dir2 < 3; dir2++)
        {
          const int ghost_a_redistribute = 2 ;
          const int n = splitting_.get_nb_items_global(loc, dir2);
          if(dir2==2)
            {
              // on ne redistribue les ghost que sur z pour le shear perio
              map[dir2].resize(1,3);
              // envoyer les rangees -2, -1, 0, 1 dans 0, 1, 2, 3
              map[dir2](0,0) = n_ext-ghost_a_redistribute;  // source index
              map[dir2](0,1) = 0; // dest index
              map[dir2](0,2) = ghost_a_redistribute*2; // size
            }
          else
            {
              map[dir2].resize(1,3);
              map[dir2](0,0) = n_ext;  // source index
              map[dir2](0,1) = 0; // dest index
              map[dir2](0,2) = n; // size
            }

        }
      redistribute_from_splitting_ft_elem_ghostz_min_.initialize(splitting_ft_, splitting_, loc, map);

      for (int dir2 = 0; dir2 < 3; dir2++)
        {
          const int ghost_a_redistribute = 2 ;
          const int n = splitting_.get_nb_items_global(loc, dir2);
          const int n_ft = splitting_ft_.get_nb_items_global(loc, dir2);
          if(dir2==2)
            {
              // on ne redistribue les ghost que sur z pour le shear perio
              map[dir2].resize(1,3);
              // envoyer les rangees n-2, n-1, n, n+1 dans 4, 5, 6, 7
              map[dir2](0,0) = n_ft - n_ext - ghost_a_redistribute;  // source index
              map[dir2](0,1) = n - 2*ghost_a_redistribute; // dest index
              map[dir2](0,2) = ghost_a_redistribute*2; // size
            }
          else
            {
              map[dir2].resize(1,3);
              map[dir2](0,0) = n_ext;  // source index
              map[dir2](0,1) = 0; // dest index
              map[dir2](0,2) = n; // size
            }

        }
      redistribute_from_splitting_ft_elem_ghostz_max_.initialize(splitting_ft_, splitting_, loc, map);
    }
  }

  // Preparation de l'expression derivee de l'acceleration
  std::string tmpstring(expression_derivee_acceleration_);
  parser_derivee_acceleration_.setString(tmpstring);
  parser_derivee_acceleration_.setNbVar(9);
  parser_derivee_acceleration_.addVar("rappel_moyen");
  parser_derivee_acceleration_.addVar("force");
  parser_derivee_acceleration_.addVar("v_moyen");
  parser_derivee_acceleration_.addVar("ur");
  parser_derivee_acceleration_.addVar("ul");
  parser_derivee_acceleration_.addVar("uv");
  parser_derivee_acceleration_.addVar("T");
  parser_derivee_acceleration_.addVar("rhov_moyen");
  parser_derivee_acceleration_.addVar("tauw");
  parser_derivee_acceleration_.parseString();

  std::string tmpstring2(expression_derivee_facteur_variable_source_);
  parser_derivee_facteur_variable_source_.setString(tmpstring2);
  parser_derivee_facteur_variable_source_.setNbVar(6);
  parser_derivee_facteur_variable_source_.addVar("rappel_moyen");
  parser_derivee_facteur_variable_source_.addVar("facteur_sv");
  parser_derivee_facteur_variable_source_.addVar("v_moyen");
  parser_derivee_facteur_variable_source_.addVar("T");
  parser_derivee_facteur_variable_source_.addVar("rhov_moyen");
  parser_derivee_facteur_variable_source_.addVar("tauw");
  parser_derivee_facteur_variable_source_.parseString();

  if (nom_reprise_ != "??")
    reprendre_probleme(nom_reprise_);

  interfaces_.associer(*this);

  /*
   * TODO: Change this block with DERIV CLASS IJK_Thermal
   */
  for (auto& itr : thermique_)
    itr.associer(*this);
  for (auto& itr : energie_)
    itr.associer(*this);

  thermals_.associer(*this);
  first_step_interface_smoothing_ = (first_step_interface_smoothing_ && !(*this).reprise_);

  run();
  return is;
}

Sortie& IJK_FT_double::printOn(Sortie& os) const
{
  // Objet_U::printOn(os);
  return os;
}

Entree& IJK_FT_double::readOn(Entree& is)
{
  // Objet_U::readOn(is);
  return is;
}

const IJK_Field_double& IJK_FT_double::get_IJK_field(const Nom& nom) const
{
  /*
  const Motcles & liste_fields = liste_post_instantanes_;
  if (!liste_fields.contient_(nom)) {
    Cerr << "Erreur dans IJK_FT_double::get_IJK_field : "
   << "Champ demande : " << nom
   << "Liste des champs possibles : " << liste_fields << finl;
    Process::exit();
  }
  */
  // FixedVector<IJK_Field_double, 3> velocity_;
  // IJK_Field_double & velocity = select(direction, vx, vy, vz);

  // Dans ce cas, le champ velocity_ft_ n'est pas utilise :
  if (disable_diphasique_)
    {
      if (nom== "VELOCITY_X")
        return velocity_[0];
      if (nom== "VELOCITY_Y")
        return velocity_[1];
      if (nom== "VELOCITY_Z")
        return velocity_[2];
    }

  if (nom== "VELOCITY_X")
    return velocity_ft_[0];
  if (nom== "VELOCITY_Y")
    return velocity_ft_[1];
  if (nom== "VELOCITY_Z")
    return velocity_ft_[2];
  if (nom== "D_VELOCITY_X")
    return d_velocity_[0];
  if (nom== "D_VELOCITY_Y")
    return d_velocity_[1];
  if (nom== "D_VELOCITY_Z")
    return d_velocity_[2];
  if (nom== "INDICATRICE")
    return interfaces_.I_ft();
  if (nom== "PRESSURE")
    return pressure_;

  return post_.get_IJK_field(nom);
}

void IJK_FT_double::force_entry_velocity(IJK_Field_double& vx, IJK_Field_double& vy, IJK_Field_double& vz, double v_imposed)
{
  const IJK_Splitting& splitting = vx.get_splitting();
  const int offset_i = splitting.get_offset_local(DIRECTION_I);
  if (offset_i > 0)
    return;
  {
    double imposed[3] = {0., 0., 0.};
    imposed[0] = v_imposed;
    for (int direction = 0; direction < 3; direction++)
      {
        IJK_Field_double& velocity = select(direction, vx, vy, vz);
        const int imin = 0;
        const int jmin = 0;
        const int kmin = 0;
        const int imax = 3;
        const int jmax = velocity.nj();
        const int kmax = velocity.nk();
        for (int k = kmin; k < kmax; k++)
          for (int j = jmin; j < jmax; j++)
            for (int i = imin; i < imax; i++)
              velocity(i,j,k) = imposed[direction];
      }
  }
}

void IJK_FT_double::force_upstream_velocity(IJK_Field_double& vx, IJK_Field_double& vy, IJK_Field_double& vz,
                                            double v_imposed,
                                            const IJK_Interfaces& interfaces,
                                            double nb_diam,
                                            int upstream_dir,
                                            int gravity_dir,
                                            int upstream_stencil)
{
  int dir = 0;
  if (upstream_dir == -1)
    {
      dir = gravity_dir;
      if (dir == -1)
        dir=0;
    }
  const IJK_Splitting& splitting = vx.get_splitting();
  const IJK_Grid_Geometry& geom = splitting.get_grid_geometry();

  bool perio = geom.get_periodic_flag(dir);

  assert(interfaces.get_nb_bulles_reelles() == 1);
  DoubleTab bounding_box;
  Cerr << "Upstream Velocity - Compute Bounding box" << finl;
  interfaces.calculer_bounding_box_bulles(bounding_box);
  // Calcule la hauteur en x de la permiere bulle et la position de son cdg :
  const double Dbdir = bounding_box(0, dir, 1) - bounding_box(0, dir, 0);
  const double dirb  = ( bounding_box(0, dir, 1) + bounding_box(0, dir, 0) ) / 2.;
  const double ldir = geom.get_domain_length(dir) ;
  if (nb_diam == 0.)
    nb_diam = (ldir/Dbdir) / 2;
  double dirobj = dirb + nb_diam * Dbdir;

  // L'origine est sur un noeud. Donc que la premiere face en I est sur get_origin(DIRECTION_I)
  const double ddir = geom.get_constant_delta(dir);
  const double origin_dir = geom.get_origin(dir) ;
  const int offset_dir = splitting.get_offset_local(dir);

  // FIXME: If nb_diam is too large it will iterate a lot
  if (perio)
    {
      while (dirobj<origin_dir)
        dirobj += ldir;
      while (dirobj>origin_dir+ldir)
        dirobj -= ldir;
    }

  // On devrait avoir xobj dans le domaine, sinon, on a choisi nb_diam trop grand :
  assert( ((dirobj>=origin_dir) && (dirobj <= origin_dir+ldir) ));

  const double x2 = (dirobj - origin_dir)/ ddir;
  int index_dir = (int)(floor(x2)) - offset_dir; // C'est l'index local, donc potentiellement negatif...
  int ndir;
  switch(dir)
    {
    case 0:
      ndir = vx.ni();
      break;
    case 1:
      ndir = vx.nj();
      break;
    case 2:
      ndir = vx.nk();
      break;
    default:
      ndir = vx.ni();
      break;
    }
  // Cerr << "index_dir " << index_dir << finl;
  if ((index_dir >=0) && (index_dir < ndir))
    {
      // On est sur le bon proc...
      if (index_dir+upstream_stencil >= ndir)
        {
          // On ne veut pas s'embeter sur 2 procs...
          index_dir = ndir-upstream_stencil;
        }
    }
  else
    {
      return;
    }
  {
    double imposed[3] = {0., 0., 0.};
    imposed[dir] = v_imposed;
    for (int direction = 0; direction < 3; direction++)
      {
        IJK_Field_double& velocity = select(direction, vx, vy, vz);
        int imin;
        int jmin;
        int kmin;
        int imax;
        int jmax;
        int kmax;
        switch (dir)
          {
          case 0:
            imin = index_dir;
            jmin = 0;
            kmin = 0;
            imax = imin+upstream_stencil;
            jmax = velocity.nj();
            kmax = velocity.nk();
            break;
          case 1:
            imin = 0;
            jmin = index_dir;
            kmin = 0;
            imax = velocity.ni();
            jmax = jmin+upstream_stencil;
            kmax = velocity.nk();
            break;
          case 2:
            imin = 0;
            jmin = 0;
            kmin = index_dir;
            imax = velocity.ni();
            jmax = velocity.nj();
            kmax = kmin+upstream_stencil;
            break;
          default:
            imin = index_dir;
            jmin = 0;
            kmin = 0;
            imax = imin+upstream_stencil;
            jmax = velocity.nj();
            kmax = velocity.nk();
            break;
          }
        for (int k = kmin; k < kmax; k++)
          for (int j = jmin; j < jmax; j++)
            for (int i = imin; i < imax; i++)
              velocity(i,j,k) = imposed[direction];
      }
  }
  Cerr << "Upstream Velocity has been forced" << finl;
}

void IJK_FT_double::ecrire_donnees(const FixedVector<IJK_Field_double, 3>& f3compo, SFichier& le_fichier, const int compo, bool binary) const
{
  const IJK_Field_double& f =  f3compo[compo];

  ArrOfDouble coord_i, coord_j, coord_k;
  build_local_coords(f, coord_i, coord_j, coord_k);

  const int ni = f.ni();
  const int nj = f.nj();
  const int nk = f.nk();

  int cnt = 0;
  for (int k = 0; k < nk; k++)
    for (int j = 0; j < nj; j++)
      for (int i = 0; i < ni; i++)
        {
          // le_fichier << coord_i[i] << Separateur::SPACE << coord_j[j] << Separateur::SPACE << coord_k[k] << f(i,j,k) << Separateur::SPACE;
          le_fichier << coord_i[i] << coord_j[j] << coord_k[k] << f(i,j,k);
          cnt++;
        }

  const IJK_Grid_Geometry& geom = get_geometry();
  const int idx_min = f.get_splitting().get_offset_local(compo);
  if ((idx_min == 0) &&  (geom.get_periodic_flag(compo)))
    {
      double l = geom.get_domain_length(compo) + geom.get_origin(compo);
      if (compo == 0)
        {
          Cerr << "compo 0: " << l << " " << coord_i[0] << finl;
          for (int k = 0; k < nk; k++)
            for (int j = 0; j < nj; j++)
              {
                le_fichier << l << coord_j[j] << coord_k[k] << f(0,j,k);
                cnt++;
              }
        }
      if (compo == 1)
        {
          Cerr << "compo 1: " << l << " " << coord_j[0] << finl;
          for (int k = 0; k < nk; k++)
            for (int i = 0; i < ni; i++)
              {
                le_fichier << coord_i[i] << l << coord_k[k] << f(i,0,k);
                cnt++;
              }
        }
      if (compo == 2)
        {
          for (int j = 0; j < nj; j++)
            for (int i = 0; i < ni; i++)
              {
                le_fichier << coord_i[i] << coord_j[j] << l << f(i,j,0);
                cnt++;
              }
        }
    }
  le_fichier.flush();
  Cerr << "Written " << cnt << " items with (ni, nj, nk)." << ni << " " << nj << " " << nk << finl;
}

// Initialize field with specified string expression (must be understood by Parser class)
void IJK_FT_double::dumpxyz_vector(const FixedVector<IJK_Field_double, 3>& f3compo, const char * filename, bool binary) const
{
  int np = Process::nproc();
  int rank = Process::me();

  Process::barrier();
  int token = 1;
  if (Process::je_suis_maitre())
    {
      // Write and send token to rank+1
      SFichier le_fichier;
      le_fichier.set_bin(1);
      le_fichier.ouvrir(filename);
      for (unsigned compo=0; compo<3; compo++)
        ecrire_donnees(f3compo, le_fichier, compo, binary);
      if (np > 1)
        envoyer(token, rank, rank+1, 2345);
    }
  else
    {
      int rcv;
      recevoir(rcv, rank-1, rank, 2345);

      SFichier le_fichier;
      le_fichier.set_bin(1);
      le_fichier.ouvrir(filename, ios::app);  // in append mode!
      for (unsigned compo=0; compo<3; compo++)
        ecrire_donnees(f3compo, le_fichier, compo, binary);

      if (rank != np-1)
        envoyer(token, rank, rank+1, 2345);
    }

  Process::barrier();
  Cerr << "Fin de l ecriture dans le fichier XYZ: " << filename << finl;
}


void IJK_FT_double::sauvegarder_probleme(const char *fichier_sauvegarde,
                                         const int& stop)//  const
{
  statistiques().begin_count(sauvegarde_counter_);

  Nom lata_name(fichier_sauvegarde);
  lata_name += ".lata";
  dumplata_header(lata_name, velocity_[0] /* on passe un champ pour ecrire la geometrie */);
  dumplata_newtime(lata_name,current_time_);
  dumplata_vector(lata_name,"VELOCITY", velocity_[0], velocity_[1], velocity_[2], 0);

  post_.sauvegarder_post(lata_name);

  if (sauvegarder_xyz_)
    {
      Nom xyz_name(fichier_sauvegarde);
      xyz_name += ".xyz";
      // Nom xyz_name_ascii = xyz_name + "_ascii";
      dumpxyz_vector(velocity_, xyz_name, true);
      // dumpxyz_vector(velocity_, xyz_name_ascii, false);
    }
  if (!disable_diphasique_)
    interfaces_.sauvegarder_interfaces(lata_name);

  // thermique_->sauvegarder_temperature(lata_name);
  int idx =0;
  //TODO: sauvegarde des champs surfaces (vapeur) et barycentre,
  // eventuellement du med pour voir si la conversion marche.
  /*
   * TODO: Change this block with DERIV CLASS IJK_Thermal
   */
  for (auto& itr : thermique_)
    {
      itr.sauvegarder_temperature(lata_name, idx);
      ++idx;
    }
  int idx2 =0;
  for (auto& itr : energie_)
    {
      itr.sauvegarder_temperature(lata_name, idx2);
      ++idx2;
    }

  thermals_.sauvegarder_temperature(lata_name, stop);

  // curseur = thermique_; // RAZ : Remise au depart du curseur. GB -> Anida : Ne marche pas sur une liste vide? Je dois grader le curseur_bis ensuite.
  SFichier fichier;
  if (Process::je_suis_maitre())
    {
      fichier.ouvrir(fichier_sauvegarde);
      Cerr << "T= " << current_time_ << " Checkpointing dans le fichier " << fichier_sauvegarde << finl;
      fichier.precision(17);
      fichier.setf(std::ios_base::scientific);

#ifdef GAUTHIER
      Param param(que_suis_je());
      param.ajouter("tinit", &current_time_);
      param.ajouter("terme_acceleration_init", &terme_source_acceleration_);
      // param.ajouter("force_init", &force_init_);
      param.ajouter("fichier_reprise_vitesse", &fichier_reprise_vitesse_);
      param.ajouter("timestep_reprise_vitesse", &timestep_reprise_vitesse_);
      //  param.ajouter("timestep_reprise_rho", &timestep_reprise_rho_);
      param.ajouter("interfaces", & interfaces_);
      param.ajouter("statistiques", &statistiques_);
      param.ajouter("statistiques_FT", &statistiques_FT_);
      param.ajouter("groups_statistiques_FT", &groups_statistiques_FT_);
      // S'il y a plusieurs groups, on s'occupe des objets stats pour chaque group:
      // (en ecrivant directement le vecteur d'objets)
      param.print(fichier);
#else
      fichier << "{\n"
              << " tinit " << current_time_ << "\n"
              << " terme_acceleration_init " << terme_source_acceleration_ << "\n"
              // GAB : qdm_source. Les valeurs des attributs utiles pour le calcul de source_qdm_gr sont
              //       ecrits dans la reprise. Ils sont ecrits avec des mots-clefs qui n'ont pas vocation a
              //       etre dans un jdd. Ces mots doivent se trouver uniquement dans des fichiers sauv.
              << " reprise_vap_velocity_tmoy " << vap_velocity_tmoy_ << "\n"
              << " reprise_liq_velocity_tmoy " << liq_velocity_tmoy_ << "\n"
              << " fichier_reprise_vitesse " << basename(lata_name) << "\n";
      fichier << " timestep_reprise_vitesse 1\n"
              << " interfaces " << interfaces_  ;
      fichier << " forcage " << forcage_
              << " corrections_qdm " << qdm_corrections_;

      /*
       * TODO: Change this block with DERIV CLASS IJK_Thermal
       */
      /*
       * Temperature
       */
      int flag_list_not_empty = 0;
      if (thermique_.size() > 0)
        {
          fichier << " thermique {\n" ;
          flag_list_not_empty = 1;
        }
      for(auto itr = thermique_.begin(); itr != thermique_.end(); )
        {
          fichier << *itr ;
          ++itr;
          if (itr != thermique_.end())
            fichier << ", \n" ;
          else
            fichier << "\n" ;
        }
      if (flag_list_not_empty)
        fichier << " } \n" ;
      /*
       * Energy
       */
      int flag_list_not_empty_en = 0;
      if (energie_.size() > 0)
        {
          fichier << " energie {\n" ;
          flag_list_not_empty_en = 1;
        }
      for(auto itr = energie_.begin(); itr != energie_.end(); )
        {
          fichier << *itr ;
          ++itr;
          if (itr != energie_.end())
            fichier << ", \n" ;
          else
            fichier << "\n" ;
        }
      if (flag_list_not_empty_en)
        fichier << " } \n" ;
      /*
       * Thermals problems
       */
      thermals_.sauvegarder_thermals(fichier);
      /*
       *
       */
      post_.sauvegarder_post_maitre(lata_name, fichier);
      fichier << "}\n" ;
#endif
      Cerr << "T= " << current_time_ << " Checkpointing dans le fichier l.1168 " << fichier_sauvegarde << finl;
    }
  statistiques().end_count(sauvegarde_counter_);

}

void IJK_FT_double::reprendre_probleme(const char *fichier_reprise)
{
  // Lecture par tous les processeurs, on retire les commentaires etc...
  LecFicDiffuse_JDD fichier(fichier_reprise);
  Param param(que_suis_je());
  param.ajouter("tinit", &current_time_);

  param.ajouter("terme_acceleration_init", &terme_source_acceleration_);
  // param.ajouter("force_init", &force_init_);
  param.ajouter("fichier_reprise_vitesse", &fichier_reprise_vitesse_);
  param.ajouter("timestep_reprise_vitesse", &timestep_reprise_vitesse_);
  param.ajouter("interfaces", & interfaces_);
  // param.ajouter("force_init", &force_init_);

  /*
   * TODO: Change this block with DERIV CLASS IJK_Thermal
   */
  param.ajouter("thermique", &thermique_);
  param.ajouter("energie", &energie_);
  param.ajouter("thermals", &thermals_);
  param.ajouter("forcage", &forcage_);
  param.ajouter("corrections_qdm", &qdm_corrections_);

  // GAB : En chantier, (ce qui suit)
  //  param.ajouter("reprise_qdm_source", &reprise_qdm_source_);
  //  param.ajouter("reprise_qdm_source", &qdm_source_);
  //  param.ajouter("reprise_vap_velocity_tmoy", &reprise_vap_velocity_tmoy_);
  param.ajouter("reprise_vap_velocity_tmoy", &vap_velocity_tmoy_);
  //  param.ajouter("reprise_liq_velocity_tmoy", &reprise_liq_velocity_tmoy_);
  param.ajouter("reprise_liq_velocity_tmoy", &liq_velocity_tmoy_);
  /* GAB : Si on encapsule tout ce qui touche a qdm_source dans une classe,
   * on met ce qui suit dans initialise plutot.
   */
  //  param.ajouter("last_source_qdm_update_time", &last_source_qdm_update_time_);
  //  param.ajouter("offset_list_index_", &offset_list_index_);

  post_.reprendre_post(param);

  param.lire_avec_accolades(fichier);

  // Appeler ensuite initialize() pour lire les fichiers lata etc...
  Cerr << "Reprise des donnees a t=" << current_time_ << "\n" << finl;
  reprise_ = 1;
  interfaces_.set_reprise(1);
  Nom prefix = dirname(fichier_reprise);
  interfaces_.set_fichier_reprise(prefix + interfaces_.get_fichier_reprise());
  if (!thermals_.est_vide())
    thermals_.set_fichier_reprise(prefix + thermals_.get_fichier_reprise());
  fichier_reprise_vitesse_=prefix+fichier_reprise_vitesse_;
}

// Methode de calcul du pas de temps max base sur CFL, Oh et Fo
// Pour les maillages uniformes uniquement.
// On choisi de mettre le facteur 0.5 dans dt_cfl et dt_fo
// pour que le calcul soit stable avec un CFL <=1.0 et Fo <= 1.0.
// Sinon, il faudrait recommander CFL <= 0.5 et Fo <=0.5 ce qui n'est pas la valeur par defaut...
double IJK_FT_double::find_timestep(const double max_timestep,
                                    const double cfl,
                                    const double fo,
                                    const double oh)
{
  statistiques().begin_count(dt_counter_);

  dt_cfl_ = 1.e20;
  double dxmin = 1.e20;
  // On ne connait pas la longueur minimum du maillage lagrangien, mais on en prend une approximation :
  // lg \approx 1.7 delta
  double lg_cube = 1.;
  double lg_cube_raw = 1.;
  for (int dir = 0; dir < 3; dir++)
    {
      const IJK_Field_double& v = velocity_[dir];
      double max_v = 1.e-20; // Pas zero pour la division a la fin...
      const int ni = v.ni();
      const int nj = v.nj();
      const int nk = v.nk();
      for (int k = 0; k < nk; k++)
        for (int j = 0; j < nj; j++)
          for (int i = 0; i < ni; i++)
            max_v = std::max(max_v, fabs(v(i,j,k)));
      max_v = Process::mp_max(max_v);
      const IJK_Grid_Geometry& geom = v.get_splitting().get_grid_geometry();
#ifndef VARIABLE_DZ
      const double delta = geom.get_constant_delta(dir);
#else
      const ArrOfDouble& tab_dz=geom.get_delta(dir);
      const double delta = Process::mp_min(min_array(tab_dz));
#endif
      lg_cube *= 1.7*delta;
      lg_cube_raw *= delta;
      dxmin = std::min(delta, dxmin);
      // QUESTION GAB : pourquoi on ne reprend pas dxmin ?
      if (max_v>0)
        dt_cfl_ = std::min(dt_cfl_, delta / max_v * 0.5);
    }
  dt_cfl_ *= cfl;
  dt_cfl_liq_ = dt_cfl_;
  dt_cfl_vap_ = dt_cfl_;

  // const double nu_max = std::max(mu_liquide_/rho_liquide_, mu_vapeur_/rho_vapeur_);
  // dt_fo_  = dxmin*dxmin/(nu_max + 1.e-20) * fo * 0.125;
  dt_fo_liq_ = dxmin*dxmin/((mu_liquide_/rho_liquide_) + 1.e-20) * fo * 0.125;
  dt_fo_vap_ =dxmin*dxmin/((mu_vapeur_/rho_vapeur_) + 1.e-20) * fo * 0.125;
  dt_fo_ = std::min(dt_fo_liq_, dt_fo_vap_);
  if (disable_diffusion_qdm_)
    dt_fo_ = 1.e20;

  /*
   * Popinet et.al 2018 (review surface tension calculation)
   * Au cas ou sigma = 0, on utilise (sigma + 1e-20)
   */
  double ideal_length_factor = get_remaillage_ft_ijk().get_facteur_longueur_ideale();
  ideal_length_factor = (ideal_length_factor == -1) ? 1.7 : ideal_length_factor;
  if (enable_dt_oh_ideal_length_factor_)
    dt_oh_  = sqrt((rho_liquide_ + rho_vapeur_) / (2 * M_PI) * lg_cube_raw * pow(ideal_length_factor, 3) / (sigma_ + 1e-20)) * oh;
  else
    dt_oh_  = sqrt((rho_liquide_+rho_vapeur_)/2. * lg_cube/(sigma_+1e-20) ) * oh;
  if (disable_convection_qdm_)
    dt_oh_ = 1.e20;

  // Seems underestimated !
  //  const double dt_oh  = sqrt((rho_liquide_+rho_vapeur_)/2. * lg_cube/(sigma_+1e-20) ) * oh;

  const double dt_eq_velocity = 1. / (1./dt_cfl_ + 1./dt_fo_ + 1./dt_oh_);

  /*
   * TODO: Change this block with DERIV CLASS IJK_Thermal
   */
  double dt_thermique = 1.e20;
  for (const auto& itr : thermique_)
    {
      const double dt_th = itr.compute_timestep(dt_thermique, dxmin);
      // We take the most restrictive of all thermal problems and use it for all:
      dt_thermique= std::min(dt_thermique, dt_th);
    }

  double dt_energie = 1.e20;
  for (const auto& itr : energie_)
    {
      const double dt_en = itr.compute_timestep(dt_energie, dxmin);
      // We take the most restrictive of all thermal problems and use it for all:
      dt_energie= std::min(dt_energie, dt_en);
    }
  const double dt_thermique_energie = std::min(dt_energie, dt_thermique);

  double dt_thermals = 1.e20;
  thermals_.compute_timestep(dt_thermals, dxmin);
  dt_thermals = std::min(dt_thermals, dt_thermique_energie);

  const double dt = std::min(max_timestep, timestep_facsec_*std::min(dt_eq_velocity, dt_thermals));

  if (Process::je_suis_maitre())
    {
      int reset = (!reprise_) && (tstep_==0);
      SFichier fic = Ouvrir_fichier(".dt_ev","tstep\ttime\ttimestep\tdt_cfl\tdt_fo\tdt_oh\tdt_diff_th",reset);
      fic << tstep_ << " " << current_time_ << " " << dt;
      fic << " " << dt_cfl_ << " " << dt_fo_ << " " << dt_oh_;
      fic << " " << dt_thermique; // If no thermal equation, value will be large.
      fic << " " << dt_energie; // If no thermal equation, value will be large.
      fic << " " << dt_thermals; // If no thermal equation, value will be large.
      fic << finl;
      fic.close();
    }
  statistiques().end_count(dt_counter_);

  return dt;
}

// Methode appelee par run() une fois la memoire alouee pour les champs.
// Cette fonction remplit les valeurs initiales de vitesse
// Elle debute aussi les compteurs.
int IJK_FT_double::initialise()
{
  Cout << que_suis_je() << "::initialise()" << finl;
  int nalloc = 0;
#if 0
  compute_inital_velocity_spectral(velocity_);
#else

/////////////////////////////////////:
// MODIFICATIONS GAB
  Cout << "forcage_.get_type_forcage() : " <<forcage_.get_type_forcage() << finl;
  if (forcage_.get_type_forcage() > 0)
    {
      const IJK_Splitting& gbz_splitting = velocity_[0].get_splitting();
      const IJK_Grid_Geometry& my_geom = velocity_[0].get_splitting().get_grid_geometry();

      const int my_ni = velocity_[0].ni();
      const int my_nj = velocity_[0].nj();
      const int my_nk = velocity_[0].nk();
      const int nproc_tot = Process::nproc();
      Cout << "BF compute_initial_chouippe" << finl;
      Cout << "ni : "<<my_ni<<" ,nj : "<<my_nj<<" ,nk : "<<my_nk << finl;
      std::cout << "in initialise i_offset : " << gbz_splitting.get_offset_local(DIRECTION_I) << std::endl;
      std::cout << "Process::me()" << Process::me() << std::endl;
      forcage_.compute_initial_chouippe(nproc_tot,my_geom,my_ni,my_nj,my_nk,gbz_splitting,nom_sauvegarde_);
      statistiques().begin_count(m2_counter_);
      Cout << "AF compute_initial_chouippe" << finl;
    }

  if (fichier_reprise_vitesse_ == "??")   // si on ne fait pas une reprise on initialise V
    {
      if (expression_vitesse_initiale_.size() != 3)
        {
          Cerr << "Erreur dans l'initialisation: la vitesse initiale doit etre fournie avec trois expressions" << finl;
          Process::exit();
        }
      else
        {
          Cout << "Initialisation vitesse \nvx = " << expression_vitesse_initiale_[0]
               << "\nvy = " <<  expression_vitesse_initiale_[1]
               << "\nvz = " <<  expression_vitesse_initiale_[2]  << finl;
          for (int i = 0; i < 3; i++)
            {
              // Cette methode parcours ni(), nj() et nk() et donc pas les ghost...
              set_field_data(velocity_[i], expression_vitesse_initiale_[i]);
            }

          velocity_[0].echange_espace_virtuel(2);
          velocity_[1].echange_espace_virtuel(2);
          velocity_[2].echange_espace_virtuel(2);

          // permet de checker les espaces_virtuel
//          for (int dir = 0; dir < 3; dir++)
//          {
//			  std::cout << std::endl;
//			  std::cout << " ########### DIR = " << dir <<"  #############" << std::endl;
//			  std::cout << std::endl;
//        	  for (int j = -2; j < velocity_[dir].nj() + 2; j++)
//        	  {
//        		  std::cout << "j=" << j << std::endl;
//    			  std::cout << std::endl;
//
//        		  for (int i = velocity_[dir].ni() + 1; i > -3; i--)
//        		  {
//        			  std::cout << "i=" << i << " : ";
//        			  for (int k = -2; k < velocity_[dir].nk() + 2; k++)
//        			  {
//        				  std::cout << velocity_[dir](i,j,k) - 1. << " ";
//        			  }
//        			  std::cout << std::endl;
//        		  }
//
//    			  std::cout << std::endl;
//        	  }
//          }
        }
    }
  else
    {
      Cout << "Lecture vitesse initiale dans fichier " << fichier_reprise_vitesse_ << " timestep= " << timestep_reprise_vitesse_ << finl;
      const Nom& geom_name = velocity_[0].get_splitting().get_grid_geometry().le_nom();
      lire_dans_lata(fichier_reprise_vitesse_, timestep_reprise_vitesse_, geom_name, "VELOCITY",
                     velocity_[0], velocity_[1], velocity_[2]); // fonction qui lit un champ a partir d'un lata .

      if (add_initial_field_)
        {
          compose_field_data(velocity_[0], expression_vitesse_initiale_[0]);
          compose_field_data(velocity_[1], expression_vitesse_initiale_[1]);
          compose_field_data(velocity_[2], expression_vitesse_initiale_[2]);
//          velocity_[0].data() += expression_vitesse_initiale_;
//          velocity_[1].data() += expression_vitesse_initiale_;
//          velocity_[2].data() += expression_vitesse_initiale_;
        }

      velocity_[0].echange_espace_virtuel(2, boundary_conditions_.get_dU_perio(boundary_conditions_.get_resolution_u_prime_()));
      velocity_[1].echange_espace_virtuel(2);
      velocity_[2].echange_espace_virtuel(2);
#ifdef CONVERT_AT_READING_FROM_NURESAFE_TO_ADIM_TRYGGVASON_FOR_LIQUID_VELOCITY
      const double coef = 14.353432757182377;
      for (int dir=0; dir< 3; dir++)
        velocity_[dir].data() *= coef;
#endif
    }
#endif

  // Pour le check_stats_ ou pour travailler en increment de pression, il faut connaitre la pression initiale :
  if (expression_pression_initiale_ != "??")
    {
      Cout << "Initialisation pression \nPini = " << expression_pression_initiale_ << finl;
      set_field_data(pressure_, expression_pression_initiale_);
      pressure_.echange_espace_virtuel(pressure_.ghost());
    }

  // On peut recuperer le domainevf:
  const Domaine_dis& domaine_dis = refprobleme_ft_disc_.valeur().domaine_dis();

  // TODO: a valider
  // if (!disable_diphasique_)
  nalloc += interfaces_.initialize(splitting_ft_, splitting_, domaine_dis, thermal_probes_ghost_cells_);

  /*
   * Compute mean rho_g using the indicator function
   */
  if (compute_force_init_)
    {
      terme_source_acceleration_ = 0.;
      double indicator_sum = 0;
      double terme_source_acceleration_increment;
      const int ni = interfaces_.I().ni();
      const int nj = interfaces_.I().nj();
      const int nk = interfaces_.I().nk();
      for (int k = 0; k < nk; k++)
        for (int j = 0; j < nj; j++)
          for (int i = 0; i < ni; i++)
            {
              const double indic = interfaces_.I(i,j,k);
              indicator_sum += indic;
              terme_source_acceleration_increment = rho_liquide_ * indic + rho_vapeur_ * (1-indic);
              terme_source_acceleration_ += terme_source_acceleration_increment;
            }
      indicator_sum = Process::mp_sum(indicator_sum);
      terme_source_acceleration_ = Process::mp_sum(terme_source_acceleration_);
      const double gravite_norm = sqrt(gravite_[0]*gravite_[0] + gravite_[1]*gravite_[1] + gravite_[2]*gravite_[2]);
      terme_source_acceleration_ /= indicator_sum;
      terme_source_acceleration_ *= gravite_norm;
      Cout << "Calculation of force init due to gravity mean(rho_g): " << terme_source_acceleration_ << finl;
    }

  nalloc += post_.initialise(reprise_);

  // statistiques...
  nalloc += post_.initialise_stats(splitting_, vol_bulles_, vol_bulle_monodisperse_);

  if (coef_immobilisation_ > 1e-16)
    {
      nalloc +=3;
      allocate_velocity(force_rappel_, splitting_, 2);
      allocate_velocity(force_rappel_ft_, splitting_ft_, 2);
      // A la reprise, c'est fait par le IJK_Interfaces::readOn
      if (interfaces_.get_flag_positions_reference() == 0) // (!reprise_)
        {
          Cerr << "Saving interfacial positions as references." <<finl;
          interfaces_.set_positions_reference();
        }
    }
  // L'indicatrice non-perturbee est remplie (si besoin, cad si post-traitement) par le post.complete()
  post_.complete(reprise_);

  // On la met a jour 2 fois, une fois next et une fois old
  update_twice_indicator_field();

  if (!disable_diphasique_)
  for (int i=0; i<2; i++)
    {
      interfaces_.switch_indicatrice_next_old();
      interfaces_.calculer_indicatrice_next(
        post_.potentiel(),
        gravite_,
        delta_rho,
        sigma_,
        /*Pour post-traitement : post_.rebuilt_indic()
        */
#ifdef SMOOTHING_RHO
        /* Pour le smoothing : */
        rho_field_ft_,
        rho_vapeur_,
        smooth_density_,
#endif
        current_time_, tstep_
      );
    }

  if (IJK_Splitting::defilement_ == 1)
    {
      redistribute_from_splitting_ft_elem_ghostz_min_.redistribute(interfaces_.I_ft(), I_ns_);
      redistribute_from_splitting_ft_elem_ghostz_max_.redistribute(interfaces_.I_ft(), I_ns_);
      rho_field_.set_indicatrice_ghost_zmin_(I_ns_, 0);
      rho_field_.set_indicatrice_ghost_zmax_(I_ns_, rho_field_.nk()-4);
      if (use_inv_rho_)
        {
          inv_rho_field_.set_indicatrice_ghost_zmin_(I_ns_, 0);
          inv_rho_field_.set_indicatrice_ghost_zmax_(I_ns_, inv_rho_field_.nk()-4);
        }
    }
  maj_indicatrice_rho_mu();

  static Stat_Counter_Id calculer_thermique_prop_counter_= statistiques().new_counter(2, "Calcul des prop thermiques");
  statistiques().begin_count(calculer_thermique_prop_counter_);

  /*
   * TODO: Change this block with DERIV CLASS IJK_Thermal
   */
  int idx =0;
  for (auto& itr : thermique_)
    {
      nalloc += itr.initialize(splitting_, idx);
      if (!disable_diphasique_)
        itr.update_thermal_properties();
      idx++;
    }

  int idx2 =0;
  for (auto& itr : energie_)
    {
      nalloc += itr.initialize(splitting_, idx2);
      if (!disable_diphasique_)
        itr.update_thermal_properties();
      idx2++;
    }

  /*
   * Thermal problems
   */
  interfaces_.initialise_ijk_compo_connex_bubbles_params();

  if (!thermals_.est_vide())
    {
      thermals_.initialize(splitting_, nalloc);
      thermals_.get_rising_velocities_parameters(compute_rising_velocities_,
                                                 fill_rising_velocities_);
    }
  nalloc += interfaces_.associate_rising_velocities_parameters(splitting_,
                                                               compute_rising_velocities_,
                                                               fill_rising_velocities_);

  statistiques().end_count(calculer_thermique_prop_counter_);
  Cout << "End of IJK_FT_double::initialise()" << finl;

  /*
   * TODO: Change this block with DERIV CLASS IJK_Thermal
   */
  //  if ((energie_.size() > 0) or (thermique_.size() >0) or (thermal_subresolution_.size()>0))
  if (energie_.size() > 0)
    {
      interfaces_.set_compute_surfaces_mouillees();
      update_twice_indicator_field();
    }

  if (thermals_.size_thermal_problem(Nom("onefluidenergy")) > 0)
    {
      interfaces_.set_compute_surfaces_mouillees();
      update_twice_indicator_field();
    }

  if (IJK_Splitting::defilement_ == 1)
    {
      redistribute_from_splitting_ft_elem_ghostz_min_.redistribute(interfaces_.I_ft(), I_ns_);
      redistribute_from_splitting_ft_elem_ghostz_max_.redistribute(interfaces_.I_ft(), I_ns_);
      rho_field_.set_indicatrice_ghost_zmin_(I_ns_, 0);
      rho_field_.set_indicatrice_ghost_zmax_(I_ns_, rho_field_.nk()-4);
      if (use_inv_rho_)
        {
          inv_rho_field_.set_indicatrice_ghost_zmin_(I_ns_, 0);
          inv_rho_field_.set_indicatrice_ghost_zmax_(I_ns_, inv_rho_field_.nk()-4);
        }
    }
  return nalloc;
}
/*
static double calculer_force_rappel_moy(const IJK_Field_double& vx, const IJK_Field_double& indic)
{
  const int ni = vx.ni();
  const int nj = vx.nj();
  const int nk = vx.nk();
  double v_moy = 0.;
  double indic_moy=0.;
  for (int k = 0; k < nk; k++)
    {
      for (int j = 0; j < nj; j++)
        {
          for (int i = 0; i < ni; i++)
            {
              if(indic(i,j,k)==0.0)
                {
                  v_moy += vx(i,j,k);
                  indic_moy+=1.;
                }
            }
        }
    }
  v_moy = Process::mp_sum(v_moy);
  indic_moy= Process::mp_sum(indic_moy);
  if (indic_moy!=0.)
    v_moy /= indic_moy;

  return v_moy;
}
*/

#ifdef COMPLEMENT_ANTI_DEVIATION_RESIDU
// left average - right average
static double calculer_wall_difference(const IJK_Field_double& vx)
{
  const int nj = vx.nj();
  const int ni = vx.ni();
  const IJK_Grid_Geometry& geom = vx.get_splitting().get_grid_geometry();
  // Maillage uniforme, il suffit donc de diviser par le nombre total de mailles:
  // cast en double au cas ou on voudrait faire un maillage >2 milliards
  const double n_mailles_plan_xy = ((double) geom.get_nb_elem_tot(0)) * geom.get_nb_elem_tot(1) ;
  const int kmin = vx.get_splitting().get_offset_local(DIRECTION_K);
  const IJK_Splitting::Localisation loc = vx.get_localisation();
  const int nktot = vx.get_splitting().get_nb_items_global(loc, DIRECTION_K);

  double x=0.;
  if (kmin == 0)
    {
      for (int j = 0; j < nj; j++)
        for (int i = 0; i < ni; i++)
          x += vx(i,j,0) ;
    }
  if (kmin + vx.nk() == nktot)
    {
      const int k = vx.nk()-1;
      for (int j = 0; j < nj; j++)
        for (int i = 0; i < ni; i++)
          x -= vx(i,j,k) ;
    }

  // somme sur tous les processeurs.
  x = Process::mp_sum(x);
  // Il faut diviser par n_mailles_plan_xy
  x /= n_mailles_plan_xy;
  return x;
}
#endif

static double calculer_tau_wall(const IJK_Field_double& vx, const double mu_liquide)
{
  const int nj = vx.nj();
  const int ni = vx.ni();
  const IJK_Grid_Geometry& geom = vx.get_splitting().get_grid_geometry();
  // Maillage uniforme, il suffit donc de diviser par le nombre total de mailles:
  // cast en double au cas ou on voudrait faire un maillage >2 milliards
  const double n_mailles_plan_xy = ((double) geom.get_nb_elem_tot(0)) * geom.get_nb_elem_tot(1) ;
  const int kmin = vx.get_splitting().get_offset_local(DIRECTION_K);
  const  IJK_Splitting::Localisation loc = vx.get_localisation();
  const int nktot = vx.get_splitting().get_nb_items_global(loc, DIRECTION_K);
  double tauw=0.;

#ifndef VARIABLE_DZ
  const double dz = geom.get_constant_delta(DIRECTION_K);
#else
  const ArrOfDouble& tab_dz=geom.get_delta(DIRECTION_K);
  //  const int nz = geom.get_nb_elem_tot(DIRECTION_K);
  double dz = -1.; // invalid value
  // On prend le dz sur la paroi basse :
  if (kmin == 0)
    dz = tab_dz[0];
#endif
  if (kmin == 0)
    {
      for (int j = 0; j < nj; j++)
        for (int i = 0; i < ni; i++)
          tauw += vx(i,j,0) ;
    }
  if (kmin + vx.nk() == nktot)
    {
      const int k = vx.nk()-1;
      for (int j = 0; j < nj; j++)
        for (int i = 0; i < ni; i++)
          tauw += vx(i,j,k) ;
    }

  // somme sur tous les processeurs.
  tauw = Process::mp_sum(tauw);
  // Il faut diviser par dz/2. et par 2*n_mailles_plan_xy
#ifdef VARIABLE_DZ
  // Pour definir un dz sur tous les procs, on recupere celui de la paroi basse
  // sur tous les procs, y compris ceux qui n'ont pas de mur, ou pas le mur gauche.
  dz = Process::mp_max(dz);
#endif
  tauw /= (dz*n_mailles_plan_xy);
  tauw *= mu_liquide;
  return tauw;
}


// Inspiree de la methode d'IJK_Navier_Stokes_Tools :
static void runge_kutta3_update_for_float(const double dx, double& store, double& v,
                                          const int step, double dt_tot)
{
  const double coeff_a[3] = { 0., -5./9., -153./128. };
  // Fk[0] = 1; Fk[i+1] = Fk[i] * a[i+1] + 1
  const double coeff_Fk[3] = { 1., 4./9., 15./32. };

  const double facteurF = coeff_a[step];
  const double intermediate_dt = compute_fractionnal_timestep_rk3(dt_tot, step);
  const double delta_t_divided_by_Fk = intermediate_dt / coeff_Fk[step];
  double x;
  switch(step)
    {
    case 0:
      x = dx;
      store = x;
      v += x * delta_t_divided_by_Fk;
      break;
    case 1:
      // general case, read and write F
      x = store * facteurF + dx;
      store = x;
      v += x * delta_t_divided_by_Fk;
      break;
    case 2:
      // do not write F
      x = store * facteurF + dx;
      v += x * delta_t_divided_by_Fk;
      break;
    default:
      Cerr << "Error in runge_kutta_update_for_float: wrong step" << finl;
      Process::exit();
    };
}

void IJK_FT_double::calculer_terme_source_acceleration(IJK_Field_double& vx, const double time, const double timestep,
                                                       const int rk_step)
{
  /*
   * Cette methode calcule la source de qdm a appliquer pour que les bulles soient
   * globalement fixes dans la direction de la gravite
   *  o Si le parametre source_qdm_gr vaut 1, la correction a appliquer est :
   *          vx = vx - moy^{xyz} (rho vx) - moy^{xyz} (rho vx_terminale) / moy^{xyz} (rho)
   *  o Si le parametre source_qdm_gr vaut 0, la source est determinee par :
   *          temre_force_init         --> temre_source_acceleration et par
   *          expression_derivee_force --> expression_derivee_acceleration
   * REMARQUE II : On peut envisager de faire une correction qui n'a pas besoin qu'on lui donne vx_terminale en
   *               entree. C'est ce qui a ete explore, mais qui n'a pas aboutit.
   *  */
  statistiques().begin_count(source_counter_);
  double new_time = time;
  double v_moy = calculer_v_moyen(vx);

  // S'il n'y a pas de derivee, la source est constante donc on peut sortir:
  if (expression_derivee_acceleration_ == Nom("0"))
    {
      //    terme_source_acceleration_ = 0.;
      statistiques().end_count(source_counter_);
      return;
    }

  update_rho_v();
  // GAB, rotation : pas sur de mon coup la
  // double rhov_moy = calculer_v_moyen(rho_v_[DIRECTION_I]);
  double rhov_moy = calculer_v_moyen(rho_v_[direction_gravite_]);

  double moy_rappel = 0.;
  // GAB, rotation
  const IJK_Grid_Geometry& geom = velocity_[direction_gravite_].get_splitting().get_grid_geometry();
  double vol_dom =  geom.get_domain_length(DIRECTION_I)*geom.get_domain_length(DIRECTION_J)*geom.get_domain_length(DIRECTION_K);
  if (coef_immobilisation_ > 1e-16)
    {
      // GAB, rotation, pas sur de mon coup
      // moy_rappel=calculer_v_moyen(force_rappel_[DIRECTION_I])*vol_dom;
      moy_rappel=calculer_v_moyen(force_rappel_[direction_gravite_])*vol_dom;
    }

  double tauw = calculer_tau_wall(vx, mu_liquide_);
  double derivee_acceleration = 0.;
  double derivee_facteur_sv = 0.;

  double alv = 0.;
  if (vol_bulle_monodisperse_>=0.)
    alv = interfaces_.get_nb_bulles_reelles()*vol_bulle_monodisperse_/vol_dom;
  else
    alv = 1.-calculer_v_moyen(interfaces_.I());

  if (Process::je_suis_maitre())
    {
      //
      double drho = rho_liquide_-rho_vapeur_;
      double facv = 0., facl=1.;
      if (std::fabs(alv*drho)>DMINFLOAT)
        {
          facv=1./(alv*drho);
          facl=1./((1.-alv)*drho);
        }
      double ul = facl*(rhov_moy-rho_vapeur_*v_moy);
      double uv = facv*(rho_liquide_*v_moy-rhov_moy);

      // Mise a jour de l'acceleration
      parser_derivee_acceleration_.setVar("rappel_moyen", moy_rappel);
      parser_derivee_acceleration_.setVar("force", terme_source_acceleration_);
      parser_derivee_acceleration_.setVar("v_moyen", v_moy);
      parser_derivee_acceleration_.setVar("ur", uv-ul);
      parser_derivee_acceleration_.setVar("ul", ul);
      parser_derivee_acceleration_.setVar("uv", uv);
      parser_derivee_acceleration_.setVar("T", time);
      parser_derivee_acceleration_.setVar("rhov_moyen", rhov_moy);
      parser_derivee_acceleration_.setVar("tauw", tauw);
      // Pour utiliser rho_v il faudrait deplacer cette mise a jour a un endroit ou rho
      // est a jour en fonction de l'indicatrice
      // parser_derivee_acceleration_.setVar("rho_v_moyen", rho_v_moy);
      derivee_acceleration = parser_derivee_acceleration_.eval();

      // Mise a jour de la source variable
      if (expression_derivee_facteur_variable_source_ != Nom("0"))
        {
          parser_derivee_facteur_variable_source_.setVar("rappel_moyen", moy_rappel);
          parser_derivee_facteur_variable_source_.setVar("facteur_sv", facteur_variable_source_);
          parser_derivee_facteur_variable_source_.setVar("v_moyen", v_moy);
          parser_derivee_facteur_variable_source_.setVar("T", time);
          parser_derivee_facteur_variable_source_.setVar("rhov_moyen", rhov_moy);
          parser_derivee_facteur_variable_source_.setVar("tauw", tauw);
          derivee_facteur_sv = parser_derivee_facteur_variable_source_.eval();
        }//

      if (qdm_corrections_.is_type_gb())
        {
          // calcul de terme_source_acceleration_ et de terme_source_acceleration_

          // ON NE VEUT PAS METTRE A JOUR TERME_SOURCE_ACCELERATION_ AVEC CETTE METHODE
          if ( get_time_scheme() == EULER_EXPLICITE)
            {
              terme_source_acceleration_ += derivee_acceleration * timestep;
              facteur_variable_source_ += derivee_facteur_sv * timestep;
              //facteur_variable_source_ += 0;//derivee_facteur_sv * timestep;
              new_time += timestep;
            }
          else if ( get_time_scheme() == RK3_FT )
            {
              const double intermediate_dt = compute_fractionnal_timestep_rk3( timestep, rk_step);
              runge_kutta3_update_for_float(derivee_acceleration, store_RK3_source_acc_,
                                            terme_source_acceleration_, rk_step, timestep);

              runge_kutta3_update_for_float(derivee_facteur_sv, store_RK3_fac_sv_,
                                            facteur_variable_source_, rk_step, timestep);
              Cout << "facteur_variable_source_" << facteur_variable_source_ << finl;
              //facteur_variable_source_ += 0;
              new_time += intermediate_dt;
            }
        }
      else if ( get_time_scheme() == RK3_FT )
        {
          const double intermediate_dt = compute_fractionnal_timestep_rk3( timestep/*total */, rk_step);
          runge_kutta3_update_for_float(derivee_acceleration, store_RK3_source_acc_,
                                        terme_source_acceleration_, rk_step, timestep/*total */);

          runge_kutta3_update_for_float(derivee_facteur_sv, store_RK3_fac_sv_,
                                        facteur_variable_source_, rk_step, timestep/*total */);
          new_time += intermediate_dt;
        }
    }
  envoyer_broadcast(terme_source_acceleration_, 0);

  // -----------------------------------------------------------
  // Force interface (:"force sigma") seloon x,y,z et u.Force_interface
  double fs0(0),fs1(0),fs2(0),psn(0);
  if (!disable_diphasique_)
    {
      // FORCE INTERFACIALE : on veut un terme homogene a [rho.g]=[N.m^{-3}]
      // terme_source_interfaces_ns_       est homogene a [du/dt]=[m.s^{-2}]
      //   -> in calculer_dv : "~ velocity_ += force_interf * dt ~"
      //   --> force_interf a bien eu un mass_solver_with_rho plus haut
      fs0=calculer_v_moyen(scalar_fields_product(rho_field_,terme_source_interfaces_ns_[0],0));
      fs1=calculer_v_moyen(scalar_fields_product(rho_field_,terme_source_interfaces_ns_[1],1));
      fs2=calculer_v_moyen(scalar_fields_product(rho_field_,terme_source_interfaces_ns_[2],2));
      psn=calculer_v_moyen(scalar_product(velocity_,scalar_times_vector(rho_field_,terme_source_interfaces_ns_)));
    }
  // energie cinetique (monophasique) et diphasique
  double uu(calculer_v_moyen(scalar_product(velocity_,velocity_)));
  double uru(calculer_v_moyen(scalar_product(velocity_,scalar_times_vector(rho_field_,velocity_))));
  // Force exterieur (:"force thi") selon x,y,z et acceleration_thi.acceleration_thi, force_thi.foce_thi, u.Force_THI
  double ft0(0),ft1(0),ft2(0),atat(0),ftft(0),ptn(0);
  if (forcage_.get_type_forcage() > 0)
    {
      // FORCE IMPOSEE : on veut un terme homogene a [rho.g]=[N.m^{-3}]
      // forcage_.get_force_ph2()     est homogene a [du/dt]=[m.s^{-2}]
      //   -> in compute_add_THI_force_sur_d_velocity : "~ d_velocity += forcage_.get_force_ph2() ~"
      ft0 = calculer_v_moyen(scalar_fields_product(rho_field_,forcage_.get_force_ph2()[0],0));
      ft1 = calculer_v_moyen(scalar_fields_product(rho_field_,forcage_.get_force_ph2()[1],1));
      ft2 = calculer_v_moyen(scalar_fields_product(rho_field_,forcage_.get_force_ph2()[2],2));
      atat = calculer_v_moyen(scalar_product(forcage_.get_force_ph2(),forcage_.get_force_ph2()));
      ftft = calculer_v_moyen(scalar_product(scalar_times_vector(rho_field_,forcage_.get_force_ph2()),scalar_times_vector(rho_field_,forcage_.get_force_ph2())));
      ptn = calculer_v_moyen(scalar_product(velocity_,scalar_times_vector(rho_field_,forcage_.get_force_ph2())));
    }
  // -----------------------------------------------------------

  envoyer_broadcast(terme_source_acceleration_z_, 0);

  // Impression dans le fichier _acceleration.out
  if (Process::je_suis_maitre())
    {
      // GR : 07.01.22 : ce serai pas mal de mettre une condition if (tstep % dt_post_stats_acc_ == dt_post_stats_acc_ - 1 || stop)
      //      pour alleger le dossier OUT. Voir avec GB et AB.
      // double ff=0.;
      int reset = (!reprise_) && (tstep_==0);
      SFichier fic=Ouvrir_fichier("_acceleration.out",
                                  "1.tstep\t2.time\t3.Vx\t4.rhoVx\t5.tauw\t6.da/dt\t7.NewT\t8.acceleration\t9.fac_var_source\t10.qdm_source\t11.vap_velocity_tmoy_\t12.liq_velocity_tmoy_\t13.qdm_patch_correction_[0]\t14.qdm_patch_correction_[1]\t15.qdm_patch_correction_[2]\t16.F_sigma_moyen[0]\t17.F_sigma_moyen[1]\t18.F_sigma_moyen[2]\t19.y.F_sigma\t20.u.u\t21.F_THI[0]\t22.F_THI[1]\t23.F_THI[2]\t24.A_THI.A_THI\t25.F_THI.F_THI\t26.u.F_THI\t27.u.rho.u",
                                  reset);
      // la derivee_acceleration n'est connue que sur le maitre
      fic<< tstep_<<" "<< time<<" "<<v_moy<<" "<<rhov_moy <<" "<<tauw ;
      fic  <<" "<<derivee_acceleration <<" "<<new_time <<" "<<terme_source_acceleration_;
      fic <<" "<< facteur_variable_source_;
      fic <<" "<< 0.; //qdm_source_;
      fic <<" "<< vap_velocity_tmoy_;
      fic <<" "<< liq_velocity_tmoy_;

      if (coef_immobilisation_ > 1e-16)
        fic << " " << moy_rappel;

      for (int dir = 0; dir < 3; dir++)
        fic <<" "<< 0.; //qdm_patch_correction_[dir];

      // Force interfaciale et puissance du travail des forces interfaciales
      // rho*terme_source_interfaces_ns_
      fic <<" "<< fs0; // F_sigma_moyen[0]
      fic <<" "<< fs1; // F_sigma_moyen[1]
      fic <<" "<< fs2; // F_sigma_moyen[2]
      // u.rho*terme_source_interfaces_ns_
      fic <<" "<< psn; // velocity.F_sigma

      // Energie cinetique (double)
      // u.u (qui est aussi accessible par les .txt)
      fic <<" "<< uu;

      // Force imposee et puissance du traveil de la force imposee
      // rho*F_THI
      fic <<" "<< ft0; // F_THI[0]
      fic <<" "<< ft1; // F_THI[1]
      fic <<" "<< ft2; // F_THI[2]
      fic <<" "<< atat; // A_THI.A_THI
      fic <<" "<< ftft; // F_THI.F_THI
      // u.rho*F_THI
      fic <<" "<< ptn; // velocity.F_THI
      // Energie cinetique en diphasiqeu (double)
      // u.rho.u (qui est aussi accessible par les .txt)
      fic <<" "<< uru;

      fic<<finl;
      fic.close();
    }
  statistiques().end_count(source_counter_);
}

static int decoder_numero_bulle(const int code)
{
  const int num_bulle = code >>6;
  return num_bulle;
}

void IJK_FT_double::run()
{
  splitting_.get_local_mesh_delta(DIRECTION_K, 2 /* ghost cells */,
                                  delta_z_local_);
  Cerr << "IJK_FT_double::run()" << finl;
  int nalloc = 0;
  thermal_probes_ghost_cells_ = 2;
  thermals_.compute_ghost_cell_numbers_for_subproblems(splitting_, thermal_probes_ghost_cells_);
  thermal_probes_ghost_cells_ = thermals_.get_probes_ghost_cells(thermal_probes_ghost_cells_);
  if (IJK_Splitting::defilement_ == 1)
  {
  allocate_velocity(velocity_, splitting_, 2, boundary_conditions_.get_dU_perio(boundary_conditions_.get_resolution_u_prime_()));
  }
  else
  {
  allocate_velocity(velocity_, splitting_, thermal_probes_ghost_cells_);
  }

  if (IJK_Splitting::defilement_ == 1)
    {
      if (splitting_.get_nb_elem_local(2) < 4 )
        {
          std::cout << "Nb of cells / proc in z-direction must be >=4 for shear periodic run" << std::endl;
          std::cout << "Nb of cells / proc in z-direction must be >=8 for shear periodic run if only one proc on z-direction" << std::endl;
          std::cout << "Or find an other way to stock indic_ghost_zmin and zmax than IJK_Field" << std::endl;
          Process::exit();
        }
      if (splitting_.get_offset_local(0)!=0.)
        {
          std::cout << " Shear_periodic conditions works only without splitting in i-direction " << std::endl;
          std::cout << "if splitting in i-direction --> get_neighbour_processor has to be changed" << std::endl;
          Process::exit();
        }
    }

  allocate_velocity(d_velocity_, splitting_, 1);
  nalloc += 6;
  // GAB, qdm
  if (test_etapes_et_bilan_)
    {
      allocate_velocity(rho_u_euler_av_prediction_champ_, splitting_, 1);
      allocate_velocity(rho_u_euler_av_rho_mu_ind_champ_, splitting_, 1);
      allocate_velocity(rho_du_euler_ap_prediction_champ_, splitting_, 1);
      allocate_velocity(rho_u_euler_ap_projection_champ_, splitting_, 1);
      allocate_velocity(rho_du_euler_ap_projection_champ_, splitting_, 1);
      allocate_velocity(rho_u_euler_ap_rho_mu_ind_champ_, splitting_, 1);
      allocate_velocity(terme_diffusion_local_, splitting_, 1);
      allocate_velocity(terme_pression_local_, splitting_, 1);
      allocate_velocity(terme_pression_in_ustar_local_, splitting_, 1);
      allocate_velocity(d_v_diff_et_conv_, splitting_, 1);
      allocate_velocity(terme_convection_mass_solver_, splitting_, 1);
      allocate_velocity(terme_diffusion_mass_solver_, splitting_, 1);
      nalloc += 36;
    }
  //
  pressure_ghost_cells_.allocate(splitting_, IJK_Splitting::ELEM, thermal_probes_ghost_cells_);
  pressure_ghost_cells_.data() = 0.;
  pressure_ghost_cells_.echange_espace_virtuel(pressure_ghost_cells_.ghost());
  nalloc += 1;

  // if interp_monofluide == 2 --> reconstruction uniquement sur rho, mu. Pas sur P !
  if (!disable_diphasique_ && boundary_conditions_.get_correction_interp_monofluide()==1)
    {
      pressure_.allocate(splitting_, IJK_Splitting::ELEM, 3, 0 ,1, false, 1, rho_vapeur_, rho_liquide_, use_inv_rho_in_poisson_solver_);
    }
  else
    pressure_.allocate(splitting_, IJK_Splitting::ELEM, 3);
  nalloc += 1;

  if (include_pressure_gradient_in_ustar_)
    {
      d_pressure_.allocate(splitting_, IJK_Splitting::ELEM, 1);
      if (get_time_scheme() == RK3_FT)
        {
          RK3_F_pressure_.allocate(splitting_, IJK_Splitting::ELEM, 1);
          nalloc += 1;
        }
      nalloc += 1;
    }

  // On utilise aussi rhov pour le bilan de forces et pour d'autres formes de convection...
  //  if (!(expression_derivee_acceleration_ == Nom("0")))
  allocate_velocity(rho_v_, splitting_, 1);
  nalloc += 3;

  pressure_rhs_.allocate(splitting_, IJK_Splitting::ELEM, 1);
  nalloc += 1;
  I_ns_.allocate(splitting_, IJK_Splitting::ELEM, 2);

  if (!disable_diphasique_ && (boundary_conditions_.get_correction_interp_monofluide()==1 || boundary_conditions_.get_correction_interp_monofluide()==2))
    {
      molecular_mu_.allocate(splitting_, IJK_Splitting::ELEM, 2, 0 ,1, false, 2, mu_vapeur_, mu_liquide_);
      rho_field_.allocate(splitting_, IJK_Splitting::ELEM, 2, 0 ,1, false, 2, rho_vapeur_, rho_liquide_);
      nalloc += 2;
      IJK_Splitting::rho_vap_ref_for_poisson_=rho_vapeur_;
      IJK_Splitting::rho_liq_ref_for_poisson_=rho_liquide_;
      if (use_inv_rho_)
        {
          inv_rho_field_.allocate(splitting_, IJK_Splitting::ELEM, 2, 0 ,1, false, 2, 1./rho_vapeur_, 1./rho_liquide_);
          IJK_Splitting::rho_vap_ref_for_poisson_=1./rho_vapeur_;
          IJK_Splitting::rho_liq_ref_for_poisson_=1./rho_liquide_;
        }
    }
  else
    {

      IJK_Splitting::rho_vap_ref_for_poisson_=rho_vapeur_;
      IJK_Splitting::rho_liq_ref_for_poisson_=rho_liquide_;
      molecular_mu_.allocate(splitting_, IJK_Splitting::ELEM, 2);
      rho_field_.allocate(splitting_, IJK_Splitting::ELEM, 2);
      nalloc += 2;
      if (use_inv_rho_)
        {
          inv_rho_field_.allocate(splitting_, IJK_Splitting::ELEM, 2);
          nalloc += 1;
          IJK_Splitting::rho_vap_ref_for_poisson_=1./rho_vapeur_;
          IJK_Splitting::rho_liq_ref_for_poisson_=1./rho_liquide_;
        }
    }

  //  rho_batard_.allocate(splitting_, IJK_Splitting::ELEM, 2);

  if (first_step_interface_smoothing_)
    {
      allocate_velocity(zero_field_ft_, splitting_ft_, thermal_probes_ghost_cells_);
      for (int dir = 0; dir < 3; dir++)
        zero_field_ft_[dir].data() = 0.;
      zero_field_ft_.echange_espace_virtuel();
      nalloc += 3;
    }

  if (diffusion_alternative_)
    {
      allocate_velocity(laplacien_velocity_, splitting_, 1);
      unit_.allocate(splitting_, IJK_Splitting::ELEM, 2);
      unit_.data() = 1.;
      unit_.echange_espace_virtuel(unit_.ghost());
      nalloc += 4;
    }

  if (velocity_convection_op_.get_convection_op_option_rank() == non_conservative_rhou)
    {
      div_rhou_.allocate(splitting_, IJK_Splitting::ELEM, 1);
      nalloc += 1;
    }
  allocate_velocity(psi_velocity_, splitting_, 2);
  nalloc += 3;

#ifdef SMOOTHING_RHO
  // Pour le smoothing :
  rho_field_ft_.allocate(splitting_ft_, IJK_Splitting::ELEM, 2);
  nalloc += 1;
#endif
  // champs pour post-traitement :
  // post_.alloc_fields();
  nalloc += post_.alloc_fields();

  // Allocation du terme source variable spatialement:
  int flag_variable_source = false;
  if ((expression_variable_source_[0] != "??")
      || (expression_variable_source_[1] != "??")
      || (expression_variable_source_[2] != "??")
      || (expression_potential_phi_ != "??"))
    {
      allocate_velocity(variable_source_, splitting_, 1);
      flag_variable_source = true;
      potential_phi_.allocate(splitting_, IJK_Splitting::ELEM, 1);
      nalloc += 4;
      for (int dir = 0; dir < 3; dir++)
        variable_source_[dir].data() = 0.;
      potential_phi_.data() = 0.;
    }

  // thermals_.initialize(splitting_, nalloc);

  // GB : Je ne sais pas si on a besoin d'un ghost... Je crois que oui. Lequel?
  // Si la a vitesse ft doit transporter les sommets virtuels des facettes reelles,
  // alors il faut un domaine ghost de la taille de la longueur maximale des arretes.
  // allocate_velocity(velocity_ft_, splitting_ft_, 0);
  /*
   * FIXME: Allocate based on the thermal subproblems
   * as the thermal probes necessitates several ghost cells to interpolate velocity !
   * Check the difference between elem and faces ? and for interpolation of the velocity ?
   */
  /*
   * Finally use the ns velocity field for the thermal sub-problems
   */
  // thermals_.compute_ghost_cell_numbers_for_subproblems(splitting_, ft_ghost_cells);
  // ft_ghost_cells = thermals_.get_probes_ghost_cells(ft_ghost_cells);
  int ft_ghost_cells = 4;
  allocate_velocity(velocity_ft_, splitting_ft_, ft_ghost_cells);
  //  allocate_velocity(velocity_ft_, splitting_ft_, 4);
  nalloc += 3;

  if (!disable_diphasique_)
    {
      allocate_velocity(terme_source_interfaces_ft_, splitting_ft_, 2);
      allocate_velocity(backup_terme_source_interfaces_ft_, splitting_, 2);
      // Seulement pour le calcul du bilan de forces :
      allocate_velocity(terme_source_interfaces_ns_, splitting_, 1);
      allocate_velocity(backup_terme_source_interfaces_ns_, splitting_, 1);
      // Seulement pour le calcul des statistiques :
      allocate_velocity(terme_repulsion_interfaces_ns_, splitting_, 1);
      allocate_velocity(terme_repulsion_interfaces_ft_, splitting_ft_, 1);
      allocate_velocity(terme_abs_repulsion_interfaces_ns_, splitting_, 1);
      allocate_velocity(terme_abs_repulsion_interfaces_ft_, splitting_ft_, 1);
      nalloc += 18;
    }

  // FIXME: on a oublie pleins de choses la !
  // int nalloc = 24;
  nalloc += post_.alloc_velocity_and_co(flag_variable_source);
  if (get_time_scheme() == RK3_FT)
    {
      allocate_velocity(RK3_F_velocity_, splitting_, 1);
      nalloc += 3;
      Cout << "Schema temps de type : RK3_FT" << finl;
    }
  else
    Cout << "Schema temps de type : euler_explicite" << finl;


  velocity_diffusion_op_.initialize(splitting_);
  velocity_diffusion_op_.set_bc(boundary_conditions_);
  velocity_convection_op_.initialize(splitting_);

  // Economise la memoire si pas besoin
  if (!disable_solveur_poisson_)
    poisson_solver_.initialize(splitting_);


  // C'est ici aussi qu'on alloue les champs de temperature.
  nalloc += initialise();

//  rho_field_.echange_espace_virtuel(2);
//  recalculer_rho_de_chi(chi_, rho_field_, 2);
  Cerr << " Allocating " << nalloc << " arrays, approx total size= "
       << (double)(molecular_mu_.data().size_array() * (int)sizeof(double) * nalloc)
       * 9.537E-07 << " MB per core" << finl;

// Les champs ont etes alloues.
// On peut completer les sondes car les ijk_field.get_splitting() sont a present remplis.
  post_.completer_sondes();
  post_.improved_initial_pressure_guess(improved_initial_pressure_guess_);

// Cette projection n'est pas utile en reprise.
// Elle sert uniquement a rendre le champ de vitesse initial a divergence nulle
// lorsque son expression est analytique.

  if (!disable_solveur_poisson_)
    {
      if (improved_initial_pressure_guess_)
        {
          Cerr << "Improved initial pressure" << finl;
          maj_indicatrice_rho_mu();
          if (!disable_diphasique_)
            {
              /*
               * TODO: Change this block with DERIV CLASS IJK_Thermal
               */
              for (auto& itr : thermique_)
                itr.update_thermal_properties();

              for (auto& itr : energie_)
                itr.update_thermal_properties();

              thermals_.update_thermal_properties();
            }
          // La pression n'est pas encore initialisee. elle est donc nulle.
          // Avec cette option, on essaye une initialisation basee sur le champ de pression diphasique
          // a l'equilibre, cad sans vitesse, ou a minima pour un champ a div(u)=0.

          if (!disable_diphasique_)
            {
              FixedVector<IJK_Field_double, 3>& coords = post_.coords();

              // Calcul du potentiel.
              for (int dir = 0; dir < 3; dir++)
                {
                  terme_source_interfaces_ft_[dir].data() = 0.;
                  terme_repulsion_interfaces_ft_[dir].data() = 0.;
                  terme_abs_repulsion_interfaces_ft_[dir].data() = 0.;
                }
              const double delta_rho = rho_liquide_ - rho_vapeur_;
              interfaces_.ajouter_terme_source_interfaces(
                terme_source_interfaces_ft_,
                terme_repulsion_interfaces_ft_,
                terme_abs_repulsion_interfaces_ft_
              );

              assert(interfaces_.get_nb_bulles_reelles() == 1);
              DoubleTab bounding_box;
              interfaces_.calculer_bounding_box_bulles(bounding_box);
              // Calcul la hauteur en x de la permiere bulle :
              const double Dbx = bounding_box(0, 0, 1)
                                 - bounding_box(0, 0, 0);
              const double kappa = 2. / (Dbx / 2.);

              const int ni = pressure_.ni();
              const int nj = pressure_.nj();
              const int nk = pressure_.nk();
              for (int k = 0; k < nk; k++)
                for (int j = 0; j < nj; j++)
                  for (int i = 0; i < ni; i++)
                    {
                      double phi = gravite_[0] * coords[0](i, j, k)
                                   + gravite_[1] * coords[1](i, j, k)
                                   + gravite_[2] * coords[2](i, j, k);
                      double potentiel_elem = sigma_ * kappa
                                              - delta_rho * phi;
                      // La pression est hydrostatique, cad : pressure_ = P - rho g z
                      pressure_(i, j, k) = potentiel_elem
                                           * interfaces_.I(i, j, k); // - rho_field_(i,j,k) * phi;
                    }

              // pressure gradient requires the "left" value in all directions:
              pressure_.echange_espace_virtuel(
                1 /*, IJK_Field_double::EXCHANGE_GET_AT_LEFT_IJK*/);

              // Mise a jour du champ de vitesse (avec dv = seulement le terme source)
              d_velocity_[0].data() = 0.;
              d_velocity_[1].data() = 0.;
              d_velocity_[2].data() = 0.;

              for (int dir = 0; dir < 3; dir++)
                redistribute_from_splitting_ft_faces_[dir].redistribute_add(
                  terme_source_interfaces_ft_[dir], d_velocity_[dir]);

              for (int dir = 0; dir < 3; dir++)
                {
                  const int kmax = d_velocity_[dir].nk();
                  for (int k = 0; k < kmax; k++)
                    {
                      euler_explicit_update(d_velocity_[dir], velocity_[dir],
                                            k);
                    }
                }

              if (use_inv_rho_in_poisson_solver_)
                {
                  pressure_projection_with_inv_rho(inv_rho_field_,
                                                   velocity_[0], velocity_[1], velocity_[2], pressure_,
                                                   1., pressure_rhs_, check_divergence_,
                                                   poisson_solver_);

                }
              else
                {



                  pressure_projection_with_rho(rho_field_, velocity_[0],
                                               velocity_[1], velocity_[2], pressure_, 1.,
                                               pressure_rhs_, check_divergence_, poisson_solver_);

                }

            }
          else
            {

              pressure_projection(velocity_[0], velocity_[1], velocity_[2],
                                  pressure_, 1., pressure_rhs_, check_divergence_,
                                  poisson_solver_);
            }
          copy_field_values(pressure_ghost_cells_, pressure_);
        }
    }

  const double max_timestep = timestep_;

// Si calcul monophasique, on initialise correctement rho, mu, I une fois pour toute :
  if (disable_diphasique_)
    {
      rho_field_.data() = rho_liquide_;
      rho_moyen_ = rho_liquide_;
      molecular_mu_.data() = mu_liquide_;

      // C'est deja fait dans l'initialize (aucune raison de ne pas le faire)
      // indicatrice_ns_.data() = 1.;
      // indicatrice_ns_next_.data() = 1.;

      /*
       * TODO: Change this block with DERIV CLASS IJK_Thermal
       */
      for (auto& itr : thermique_)
        {
          // To fill in fields for cp (with cp_liq) and lambda (with lambda_liq)
          itr.update_thermal_properties();
        }
      for (auto& itr : energie_)
        {
          // To fill in fields for cp (with cp_liq) and lambda (with lambda_liq)
          itr.update_thermal_properties();
        }
      thermals_.update_thermal_properties();
    }
  else
    {
      Cerr << "Cas normal diphasique IJK_FT::run()" << finl;

      /*
       * TODO: Change this block with DERIV CLASS IJK_Thermal
       */
      for (auto& itr : thermique_)
        itr.update_thermal_properties();

      for (auto& itr : energie_)
        itr.update_thermal_properties();

      thermals_.update_thermal_properties();

      const double indic_moyen = calculer_v_moyen(interfaces_.I());
      rho_moyen_ = indic_moyen*rho_liquide_ + (1-indic_moyen)*rho_vapeur_;
      if (post_.get_liste_post_instantanes().contient_("EXTERNAL_FORCE"))
        {
          for (int dir=0; dir<3; dir++)
            compute_add_external_forces(dir);
        }
    }

  if ((!disable_diphasique_) && (post_.get_liste_post_instantanes().contient_("VI")
                                 || post_.get_liste_post_instantanes().contient_("TOUS")))
    interfaces_.compute_vinterp();

  // Preparer le fichier de postraitement et postraiter la condition initiale:
  Nom lata_name = nom_du_cas();
  if (fichier_post_ != "??")
    {
      lata_name = fichier_post_;
    }
  lata_name += Nom(".lata");
  post_.postraiter_ci(lata_name, current_time_);

//  if ( 0 && disable_diphasique_
//       && (liste_post_instantanes_.contient_("CURL")))
//    {
//      Cerr << " Dans un calcul monophasique, on ne calcule pas le rotationnel, "
//           << "donc ce n'est pas la peine de le demander dans les posts!" << finl;
//      Process::exit();
//    }

  post_.compute_extended_pressures(interfaces_.maillage_ft_ijk());
//post_.compute_phase_pressures_based_on_poisson(0);
//post_.compute_phase_pressures_based_on_poisson(1);

  modified_time_ini_ = thermals_.get_modified_time();
  if (!reprise_)
    current_time_ = modified_time_ini_;

  if (!first_step_interface_smoothing_)
    {
      Cout << "BF posttraiter_champs_instantanes "
           << current_time_ << " " << tstep_ << finl;
      post_.posttraiter_champs_instantanes(lata_name, current_time_, tstep_);
      thermals_.thermal_subresolution_outputs(); // for thermal counters
      Cout << "AF posttraiter_champs_instantanes" << finl;
    }

// GB 2019.01.01 Why immobilisation? if (!disable_diphasique_ && coef_immobilisation_==0.)
  if ((!disable_diphasique_) && suppression_rejetons_)
    interfaces_.detecter_et_supprimer_rejeton(true);
  if (reprise_)
    {
      // On ecrit a la suite du fichier. Cela suppose qu'il est bien a jour.
      // L'instant initial a deja ete ecrit a la fin du calcul precedent donc on
      // ne le reecrit pas.
    }
  else
    {
      // On creer de nouveaux fichiers :
      Cout << "BF ecrire_statistiques_bulles" << finl;
      post_.ecrire_statistiques_bulles(1 /* reset files */, nom_du_cas(),
                                       gravite_, current_time_);
      Cout << "AF ecrire_statistiques_bulles" << finl;
    }

// Ecrire la valeur initiale dans les sondes :
// Ecriture de la valeur initiale seulement hors reprise
  if (!reprise_)
    post_.postraiter_sondes();

//ab-forcage-control-ecoulement-deb
  update_rho_v(); // Peut-etre pas toujours necessaire selon la formulation pour la convection?
  for (int direction = 0; direction < 3; direction++)
    store_rhov_moy_[direction] = calculer_v_moyen(rho_v_[direction]);
//ab-forcage-control-ecoulement-fin

  statistiques().end_count(initialisation_calcul_counter_);

  if (!disable_TU)
    {
      if(GET_COMM_DETAILS)
        statistiques().print_communciation_tracking_details("Statistiques d'initialisation du calcul", 0);

      statistiques().dump("Statistiques d'initialisation du calcul", 0);
      print_statistics_analyse("Statistiques d'initialisation du calcul", 0);
    }
  statistiques().reset_counters();
  statistiques().begin_count(temps_total_execution_counter_);

  int stop = 0;

// Variation de volume de chaque bulle integree au cours du pas de temps :
  ArrOfDouble var_volume_par_bulle;
  var_volume_par_bulle.set_smart_resize(1);
  for (tstep_ = 0; tstep_ < nb_timesteps_ && stop == 0; tstep_++)
    {
      statistiques().begin_count(timestep_counter_);
      if (!splitting_.get_grid_geometry().get_periodic_flag(DIRECTION_K))
        {
          force_zero_on_walls(velocity_[2]);
        }

      if (timestep_facsec_ > 0.)
        {
          timestep_ = find_timestep(max_timestep, cfl_, fo_, oh_);
        }

      // Tableau permettant de calculer la variation de volume au cours du pas de temps :
      // Si on veut le mettre en optionel, il faut faire attention a faire vivre la taille de ce tableau avec les
      // creations et destructions de ghosts :
      const int nbulles_tot = interfaces_.get_nb_bulles_reelles()
                              + interfaces_.get_nb_bulles_ghost(1/*print=1*/);
      var_volume_par_bulle.resize_array(nbulles_tot);
      var_volume_par_bulle = 0.; // Je ne suis pas sur que ce soit un bon choix. Si on ne le remet pas a zero
      //                          a chaque dt, on corrigera la petite erreur qui pouvait rester d'avant...

      compute_var_volume_par_bulle(var_volume_par_bulle);

      // Choix de l'avancement en temps :
      // euler_explicite ou RK3.
      if (get_time_scheme() == EULER_EXPLICITE)
        {
          // Deplacement des interfaces par le champ de vitesse de l'instant n :
          if (!disable_diphasique_)
            {
              // TODO: aym pour GAB, si tu veux gagner en memoire et virer le doublon n/np1 il faut
              // inserer une methode ici style "mettre_a_jour_valeur_interface_temps_n()"
              do
                {
                  first_step_interface_smoothing_ = first_step_interface_smoothing_ && (counter_first_iter_ == 2);
                  deplacer_interfaces(timestep_,
                                      -1 /* le numero du sous pas de temps est -1 si on n'est pas en rk3 */,
                                      var_volume_par_bulle,
                                      first_step_interface_smoothing_);
                  counter_first_iter_--;
                  if(first_step_interface_smoothing_)
                    {
                      thermals_.set_temperature_ini();
                      Cout << "BF posttraiter_champs_instantanes " << current_time_ << " " << tstep_ << finl;
                      post_.posttraiter_champs_instantanes(lata_name, current_time_, tstep_);
                      Cout << "AF posttraiter_champs_instantanes" << finl;
                      compute_var_volume_par_bulle(var_volume_par_bulle);
                      thermals_.set_post_pro_first_call();
                    }
                }
              while (first_step_interface_smoothing_);
              first_step_interface_smoothing_ = 0;
              parcourir_maillage();
            }
          // Mise a jour de la vitesse (utilise les positions des marqueurs, rho, mu et indic a l'instant n)
          // Retourne une vitesse mise a jour et projetee a div nulle
          euler_time_step(var_volume_par_bulle);

          // Calcul du terme source force acceleration :
          // GAB : question a Guillaume, on fait time + time_step ? 'est pas homogene non ?'
          // GAB, rotation
          // /!\ On  laisse ce calcul active meme pour source_qdm_gr_!=-1 pour toujours avoir un fichier acceleration.out rempli correctement
          calculer_terme_source_acceleration(velocity_[direction_gravite_],
                                             current_time_ + timestep_, timestep_, -1);

          // Deplacement des interfaces par le champ de vitesse :
          // met a jour la position des marqueurs, la vitesse_ft, et gere les duplicatas.
          // Ne met pas a jour rho_mu_indicatrice

          if (!disable_diphasique_) // && !marker_advection_first_)
            {
              // indicatrice (and rho, mu...) are updated from the new interface position.
              // GB 2019.01.01 It is important to keep that calculation, because without it, the interface status would be
              // set to "minimal" where it should be "parcouru".
              // GAB, qdm : rho_n v_n+1
              if (test_etapes_et_bilan_)
                {
                  calculer_rho_v(rho_field_,velocity_,rho_u_euler_av_rho_mu_ind_champ_);
                  for (int dir=0; dir<3; dir++)
                    rho_u_euler_av_rho_mu_ind_[dir] = calculer_v_moyen(rho_u_euler_av_rho_mu_ind_champ_[dir]);
                }
              maj_indicatrice_rho_mu();

              /*
               * TODO: Change this block with DERIV CLASS IJK_Thermal
               */
              for (auto& itr : thermique_)
                {
                  itr.update_thermal_properties();
                  if (itr.conserv_energy_global_)
                    {
                      const double dE = itr.E0_ - itr.compute_global_energy();
                      itr.euler_rustine_step(timestep_, dE);
                    }
                }

              thermals_.euler_rustine_step(timestep_);

              // GAB, qdm rho_n+1 v_n+1 :
              if (test_etapes_et_bilan_)
                {
                  calculer_rho_v(rho_field_,velocity_,rho_u_euler_ap_rho_mu_ind_champ_);
                  for (int dir=0; dir<3; dir++)
                    {
                      rho_u_euler_ap_rho_mu_ind_[dir] = calculer_v_moyen(rho_u_euler_ap_rho_mu_ind_champ_[dir]);
                      u_euler_ap_rho_mu_ind_[dir] = calculer_v_moyen(velocity_[dir]);
                    }
                }
            }
          else
            {
              if (test_etapes_et_bilan_)
                {
                  calculer_rho_v(rho_field_,velocity_,rho_u_euler_ap_rho_mu_ind_champ_);
                  for (int dir=0; dir<3; dir++)
                    {
                      rho_u_euler_ap_rho_mu_ind_[dir] = calculer_v_moyen(rho_u_euler_ap_rho_mu_ind_champ_[dir]);
                      u_euler_ap_rho_mu_ind_[dir] = calculer_v_moyen(velocity_[dir]);
                    }
                }
            }
        }
      else if (get_time_scheme() == RK3_FT)
        {
          double current_time_at_rk3_step = current_time_;
          // GAB, qdm : passe en attribut de classe car utilise au moment de l'ecriture de mon out
          current_time_at_rk3_step_ = current_time_;
          // Evaluation de la variation de volume accumule au cours des sous pas de temps.
          // On la laisse croitre pendant les sous dt 0 et 1 puis on la corrige a la fin du 2eme :

          // Au cas ou on soit dans un cas ou des duplicatas sont necessaires mais n'ont pas ete
          // crees, on les cree :
          if (!interfaces_.get_nb_bulles_ghost() && !disable_diphasique_)
            {
              interfaces_.creer_duplicata_bulles();
            }
          for (rk_step_ = 0; rk_step_ < 3; rk_step_++)
            {
              const double fractionnal_timestep =
                compute_fractionnal_timestep_rk3(timestep_ /* total*/,
                                                 rk_step_);

              // Mise a jour des positions des marqueurs.
              // Deplacement des interfaces par le champ de vitesse au sous pas de temps k :
              if (!disable_diphasique_)
                {
                  deplacer_interfaces_rk3(timestep_ /* total */, rk_step_,
                                          var_volume_par_bulle);
                  parcourir_maillage();
                }
              // Cerr << "RK3 : step " << rk_step << finl;
              // Mise a jour de la temperature et de la vitesse :
              rk3_sub_step(rk_step_, timestep_, fractionnal_timestep,
                           current_time_at_rk3_step);

              // GAB patch qdm : choix 1

              // GAB, qdm : rho_n v_n+1
              if (test_etapes_et_bilan_)
                {
                  calculer_rho_v(rho_field_,velocity_,rho_u_euler_av_rho_mu_ind_champ_);
                  for (int dir=0; dir<3; dir++)
                    rho_u_euler_av_rho_mu_ind_[dir] = calculer_v_moyen(rho_u_euler_av_rho_mu_ind_champ_[dir]);
                }

              // Mise a jour rho, mu et l'indicatrice a partir de la nouvelle position de l'interface :
              // (sauf au dernier sous pas de temps pour lequel c'est fait a la fin du pas de temps)
              // TODO: verifier qu'on doit bien le faire aussi au dernier sous pas de temps : rk_step != 2 &&
              // TODO aym: verifier ce bloc, qui applique les sous pas de temps RK3 de la rustine a la temperature
              if (rk_step_ != 2 && !disable_diphasique_)
                {
                  // Attention, il faut que les duplicatas soient present pour faire maj_indicatrice_rho_mu :
                  maj_indicatrice_rho_mu();
                  for (auto& itr : thermique_)
                    {
                      itr.update_thermal_properties();
                      if (itr.conserv_energy_global_)
                        {
                          const double dE = itr.E0_ - itr.compute_global_energy();
                          itr.rk3_rustine_sub_step(rk_step_, timestep_, fractionnal_timestep,
                                                   current_time_at_rk3_step, dE);
                        }
                    }

                  thermals_.rk3_rustine_sub_step(rk_step_, timestep_, fractionnal_timestep,
                                                 current_time_at_rk3_step);

                }
              // Calcul du terme source force acceleration :
              // GAB, rotation
              calculer_terme_source_acceleration(velocity_[direction_gravite_],
                                                 current_time_at_rk3_step, timestep_ /*total*/, rk_step_);


              current_time_at_rk3_step += fractionnal_timestep;
              // GAB, qdm : passe en attribut de classe car utilise au moment de l'ecriture de mon out
              current_time_at_rk3_step_ += fractionnal_timestep;
              // On ne postraite pas le sous-dt 2 car c'est fait plus bas si on post-traite le pas de temps :
              if (post_.postraiter_sous_pas_de_temps()
                  && (tstep_ % post_.dt_post() == post_.dt_post() - 1)
                  && (rk_step_ != 2))
                {
                  post_.posttraiter_champs_instantanes(lata_name, current_time_at_rk3_step, tstep_);
                }
            }
          if (!disable_diphasique_)
            {
              // Les sous-pas de temps sont termines. Il n'est plus necessaire de gerer le tableau
              // RK3_G_store_vi_. On peut donc transferer les bulles et re-creer les duplicatas :
              interfaces_.supprimer_duplicata_bulles();
              interfaces_.transferer_bulle_perio();
              // On supprime les fragments de bulles.
              //interfaces_.detecter_et_supprimer_rejeton(false);
              interfaces_.creer_duplicata_bulles();

              // Mise a jour rho, mu et l'indicatrice a partir de la nouvelle position de l'interface :
              maj_indicatrice_rho_mu();

              for (auto& itr : thermique_)
                itr.update_thermal_properties();

              for (auto& itr : energie_)
                itr.update_thermal_properties();

              thermals_.update_thermal_properties();
            }
          // GAB, qdm rho_n+1 v_n+1 :
          if (test_etapes_et_bilan_)
            {
              calculer_rho_v(rho_field_,velocity_,rho_u_euler_ap_rho_mu_ind_champ_);
              for (int dir=0; dir<3; dir++)
                {
                  rho_u_euler_ap_rho_mu_ind_[dir] = calculer_v_moyen(rho_u_euler_ap_rho_mu_ind_champ_[dir]);
                  u_euler_ap_rho_mu_ind_[dir] = calculer_v_moyen(velocity_[dir]);
                }
            }
        }
      else
        {
          Cerr << "Erreur dans le run: time_scheme " << time_scheme_
               << " inconnu!" << finl;
          Process::exit();
        }
      // ------------------------------------------------------------------
      // CORRECTION DE QUANTITE DE MOUVEMENT
      // ------------------------------------------------------------------
      // Correction de QdM : permet de controler la QdM globale, dans chaque direction
      if (!(qdm_corrections_.is_type_none()) )
        {
          set_time_for_corrections();
          if (disable_diphasique_)
            compute_and_add_qdm_corrections_monophasic();
          else
            compute_and_add_qdm_corrections();
        }
      else
        {
          Cout << "qdm_corrections_.is_type_none() : " << qdm_corrections_.is_type_none() << finl;
          Cout << "terme_source_acceleration_" << terme_source_acceleration_ << finl;
        }
      if (qdm_corrections_.write_me())
        {write_qdm_corrections_information();}
      else
        {;}
      // ------------------------------------------------------------------

      //ab-forcage-control-ecoulement-deb
      // Quel que soit le schema en temps, on corrige le bilan de qdm par le residu integre :
      // integrated_residu_ est homogene a rho*u.
      // Il faut donc appliquer le solveur masse a integrated_residu_. Pour cela, on a besoin d'un champ.
      //                              On prend psi_velocity_ qui est dispo.
      // Attention, en entree du solveur mass, il faut qqch homogene a rho*u*volume_cell...
      // On rempli donc psi_velocity avec vol * integrated_residu_

      // static Stat_Counter_Id bilanQdM_counter_ = statistiques().new_counter(2, "Bilan QdM & Corrections");
      // statistiques().begin_count(bilanQdM_counter_);
      // statistiques().end_count(bilanQdM_counter_);
      static Stat_Counter_Id bilanQdM_counter_ = statistiques().new_counter(2, "Bilan QdM & Corrections");
      statistiques().begin_count(bilanQdM_counter_);
      if ((correction_bilan_qdm_ == 3) || (correction_bilan_qdm_ == 4))
        {

#ifndef VARIABLE_DZ
          double volume = 1.;
          for (int i = 0; i < 3; i++)
            volume *= splitting_.get_grid_geometry().get_constant_delta(i);
#endif

          for (int dir = 0; dir < 3; dir++)
            {
              if ((dir == 2) && (correction_bilan_qdm_ == 4))
                {
                  // passe, on ne traite pas z...
                }
              else
                {
#ifndef VARIABLE_DZ
                  const double x = volume * integrated_residu_[dir];
                  psi_velocity_[dir].data() = x;
#endif
                  const int kmax = psi_velocity_[dir].nk();
                  for (int k = 0; k < kmax; k++)
                    {
#ifdef VARIABLE_DZ
                      const double volume = get_channel_control_volume(psi_velocity_[dir], k, delta_z_local_);
                      const double x = volume*integrated_residu_[dir];
                      psi_velocity_[dir].data() = x;
#endif
                      if (use_inv_rho_for_mass_solver_and_calculer_rho_v_)
                        {
                          Cerr
                              << "Verifier que inv_rho_field soit valide et a jour ici ... "
                              << finl;
                          Process::exit();
                          mass_solver_with_inv_rho(psi_velocity_[dir],
                                                   inv_rho_field_, delta_z_local_, k);
                        }
                      else
                        {
                          mass_solver_with_rho(psi_velocity_[dir], rho_field_,
                                               delta_z_local_, k);
                        }
                      const int imax = velocity_[dir].ni();
                      const int jmax = velocity_[dir].nj();
                      for (int j = 0; j < jmax; j++)
                        {
                          for (int i = 0; i < imax; i++)
                            {
                              velocity_[dir](i, j, k) -= psi_velocity_[dir](i,
                                                                            j, k);
                            }
                        }
                    }


                  if (boundary_conditions_.get_correction_conserv_qdm()==2)
                    {
                      update_rho_v();
                      rho_field_.echange_espace_virtuel(rho_field_.ghost());
                      update_v_ghost_from_rho_v();
                    }
                  else
                    {
                      velocity_[dir].echange_espace_virtuel(
                        velocity_[dir].ghost());
                    }


                }
            }
          // Ces operations ont modifie le store_rhov_moy_ qu'il faut donc updater :
          for (int dir = 0; dir < 3; dir++)
            {
              store_rhov_moy_[dir] -= integrated_residu_[dir];
            }

          // Remise a zero du residu integre puisqu'il a ete corrige :
          integrated_residu_ = 0.;

        }
      statistiques().end_count(bilanQdM_counter_);

      //ab-forcage-control-ecoulement-fin
      current_time_ += timestep_;
      // stock dans le spliting le decallage periodique total avec condition de shear (current_time_) et celui du pas de temps (timestep_)
      IJK_Splitting::shear_x_time_ = boundary_conditions_.get_dU_perio()*(current_time_ + boundary_conditions_.get_t0_shear());
      //IJK_Splitting::shear_x_DT_ = boundary_conditions_.get_dU_perio()*timestep_;

      if (current_time_ >= post_.t_debut_statistiques())
        {
          if (boundary_conditions_.get_correction_conserv_qdm()==2)
            {
              update_rho_v();
              rho_field_.echange_espace_virtuel(rho_field_.ghost());
              update_v_ghost_from_rho_v();
            }
          else
            {
              // FA AT 16/07/2013 pensent que necessaire pour le calcul des derivees dans statistiques_.update_stat_k(...)
              // Je ne sais pas si c'est utile, mais j'assure...
              velocity_[0].echange_espace_virtuel(
                2 /*, IJK_Field_ST::EXCHANGE_GET_AT_RIGHT_I*/);
              velocity_[1].echange_espace_virtuel(
                2 /*, IJK_Field_ST::EXCHANGE_GET_AT_RIGHT_J*/);
              velocity_[2].echange_espace_virtuel(
                2 /*, IJK_Field_ST::EXCHANGE_GET_AT_RIGHT_K*/);
            }


          pressure_.echange_espace_virtuel(1);

          // C'est update_stat_ft qui gere s'il y a plusieurs groupes
          // pour faire la vraie indicatrice + les groupes
          post_.update_stat_ft(timestep_);
          if (!disable_diphasique_)
            {
              post_.compute_extended_pressures(interfaces_.maillage_ft_ijk());
              //post_.compute_phase_pressures_based_on_poisson(0);
              //post_.compute_phase_pressures_based_on_poisson(1);
            }
        }

      // Calcul du terme source d'acceleration deplacee dans les iterations du rk3 ou dans l'iteration d'euler.

      // verification du fichier stop
      stop = 0;
      if (check_stop_file_ != "??")
        {
          if (je_suis_maitre())
            {
              EFichier f;
              stop = f.ouvrir(check_stop_file_);
              if (stop)
                {
                  // file exists, check if it contains 1:
                  f >> stop;
                }
            }
          envoyer_broadcast(stop, 0);
        }
      if (tstep_ == nb_timesteps_ - 1)
        stop = 1;
      if (current_time_ >= max_simu_time_)
        stop = 1;

      if (tstep_ % dt_sauvegarde_ == dt_sauvegarde_ - 1 || stop)
        {
          // Choix : On supprime les duplicatas pour la sauvegarde.
          // On pourrait tres bien tout garder. ca serait plus leger en CPU, plus lourd en espace disque.
          if (!disable_diphasique_)
            interfaces_.supprimer_duplicata_bulles();

          sauvegarder_probleme(nom_sauvegarde_, stop);
          if (!disable_diphasique_)
            {
              // On les recree :
              interfaces_.creer_duplicata_bulles();

              // Be on the safe side, on met a jour :
              //   A la suppression des duplicatas, on avait fait mesh.supprimer_facettes qui remet le maillage
              //   a l'etat MINIMAL. Pour les post-tt sur l'interface (eg ai_ft_), il faut que le statut du maillage
              //   soit >= PARCOURU. C'est fait au debut de maj_indicatrice_rho_mu dans
              //   IJK_Interfaces::calculer_indicatrice.
              const double delta_rho = rho_liquide_ - rho_vapeur_;
              interfaces_.calculer_indicatrice_next(
                post_.potentiel(),
                gravite_,
                delta_rho,
                sigma_,
                /*Pour post-traitement : post_.rebuilt_indic()
                */
#ifdef SMOOTHING_RHO
                /* Pour le smoothing : */
                rho_field_ft_,
                rho_vapeur_,
                smooth_density_,
#endif
                current_time_, tstep_
              );
            }
        }

      // TODO: on pourrait mutualiser tous les parcourir maillages dans IJK_Interface au moment du transport de l'interface
      // et le supprimer de IJK_FT
      // interfaces_.parcourir_maillage();
      if ((!disable_diphasique_) && (post_.get_liste_post_instantanes().contient_("VI")))
        interfaces_.compute_vinterp();
      post_.postraiter_fin(stop, tstep_, current_time_, timestep_, lata_name,
                           gravite_, nom_du_cas());
      statistiques().end_count(timestep_counter_);

      if(JUMP_3_FIRST_STEPS && tstep_ < 3)
        {
          //demarrage des compteurs CPU
          if(tstep_ == 2)
            statistiques().set_three_first_steps_elapsed(true);
        }
      else
        statistiques().compute_avg_min_max_var_per_step(tstep_);


    }
  if (Process::je_suis_maitre())
    {
      SFichier master_file;
      master_file.ouvrir(lata_name, ios::app);
      master_file << "FIN" << finl;
      master_file.close();
    }
// Pour forcer l'ecriture du dernier pas de temps dans la sonde (peut-etre deja ecrit...)
// Alan 2020/03/02 : effectivement, deja ecrit
// post_.postraiter_sondes();

  if (!disable_TU)
    {
      if(GET_COMM_DETAILS)
        statistiques().print_communciation_tracking_details("Statistiques de resolution du probleme", 1);

      statistiques().dump("Statistiques de resolution du probleme", 1);
      print_statistics_analyse("Statistiques de resolution du probleme", 1);
    }

  statistiques().reset_counters();
  statistiques().begin_count(temps_total_execution_counter_);

}

// force_tot est homogene a une force volumique (comparable a la source S), comme tout
// ce qu'on evalue ici.
// Attention, cette methode ne fait pas que des evaluations, elle calcule aussi
// terme_source_correction_x_ ou _y_
void IJK_FT_double::compute_correction_for_momentum_balance(const int rk_step)
{
  // Toutes les interfaces etant fermees, on a donc la somme des forces
  // donnee par la poussee d'archimede : delta_rho * g * alpha
  double vol_cell = 1., vol_NS = 1., vol_gaz = 0.;
  const IJK_Grid_Geometry& geom_NS = get_geometry();
  for (int direction = 0; direction < 3; direction++)
    {
      vol_NS *= geom_NS.get_domain_length(direction);
      vol_cell *= geom_NS.get_constant_delta(direction);
    }

  ArrOfDouble volume;
  DoubleTab position;
  const int nbulles = interfaces_.get_nb_bulles_reelles();
  // La methode calcule a present les volumes meme pour les bulles ghost.
  // Pour les enlever, il suffit simplement de reduire la taille du tableau :
  interfaces_.calculer_volume_bulles(volume, position);
  volume.resize_array(nbulles);
  position.resize(nbulles,3);
  for (int ib=0; ib < nbulles; ib++)
    vol_gaz += volume[ib];

  update_rho_v();
  Vecteur3 force_tot, force_theo, rhov_moy, tauw;
  Vecteur3 acc, residu; // homogenes a d(rhou)/dt
#ifdef COMPLEMENT_ANTI_DEVIATION_RESIDU
  double moins_delta_Pwall_sur_h = 0.;
#endif
  tauw = 0.;
  // Flottaison opposee a la gravite :
  const double drho_alpha = -(rho_liquide_-rho_vapeur_)*vol_gaz/vol_NS;
  const double un_sur_h = 2./ geom_NS.get_domain_length(DIRECTION_K);
  double fractional_dt;
  if ( get_time_scheme() == RK3_FT )
    {
      fractional_dt = compute_fractionnal_timestep_rk3(timestep_ /* total*/, rk_step);
    }
  else
    {
      // On suppose que c'est euler :
      fractional_dt = timestep_;
    }

  // Convertir en des contributions volumiques (homogene au terme source) :
  for (int direction = 0; direction < 3; direction++)
    {
      if (direction != DIRECTION_K)
        {
          // Calcul du frottement moyen sur un mur (selon X ou Y)
          tauw[direction] = calculer_tau_wall(velocity_[direction], mu_liquide_);
#ifdef COMPLEMENT_ANTI_DEVIATION_RESIDU
        }
      else
        {
          // On n'est pas certain qu'il (le gradient normal de Uz moyenne sur un mur)
          // soit nul en instantane. Pour verifier, on le calcule:
          tauw[direction] = calculer_tau_wall(velocity_[direction], mu_liquide_);
          // La fonction calculer_tau_wall suppose que le champ de vitesse est a dz/2 du mur.
          // Pour la vitesse uz, il est en fait situe a dz. Donc il faut le doubler:
          tauw[direction] *=2. ;
          // Il y a aussi potentiellement une partie due a la pression en instantanne:
          // (dans la direction normale aux parois seulement)
          // +/-?  (-Pw^+ + Pw^- ) / h :
          moins_delta_Pwall_sur_h = -calculer_wall_difference(pressure_) * un_sur_h;
#endif
        }
      force_tot[direction] = calculer_v_moyen(terme_source_interfaces_ns_[direction]) / vol_cell;
      if (interfaces_.is_terme_gravite_rhog())
        {
          force_theo[direction] = 0.;
        }
      else
        {
          force_theo[direction] = drho_alpha*gravite_[direction];
        }
      rhov_moy[direction] = calculer_v_moyen(rho_v_[direction]);
      acc[direction] =  (rhov_moy[direction]-store_rhov_moy_[direction]) /fractional_dt;
      store_rhov_moy_[direction] = rhov_moy[direction];
      // Une partie de la force en x ou en z et la force en y sont des erreurs de discretisation,
      // que l'on peut corriger globalement :
      terme_source_correction_[direction] = force_theo[direction]-force_tot[direction] ;
      // Le frottement moyen doit etre divise par h pour etre rendu "volumique"
      // Il est oriente vers z- donc signe "-"
      tauw[direction] *= -un_sur_h;
    }

  if (tstep_ == 0)
    {
      // L'acceleration n'est pas valide car  store_rhov_moy_^0 = rhov_moy^0 .. donc on ferait (0-0)/0...
      // Donc le residu est invalide.
      // Donc on ne l'ajoute pas a l'integrated_residu_, qui vaut 0 a ce moment la...
      residu = 0.;
      integrated_residu_ = 0.;
    }
  else
    {
      for (int dir = 0; dir < 3; dir++)
        {
          residu[dir] = acc[dir] - (force_tot[dir]+tauw[dir]
                                    +correction_force_[dir]*terme_source_correction_[dir]);
          // if (dir == DIRECTION_I)
          //   {
          //     residu[dir] -= terme_source_acceleration_;
          //   }
          if (dir == direction_gravite_)
            {
              residu[dir] -= terme_source_acceleration_;
            }
          if (interfaces_.is_terme_gravite_rhog())
            {
              residu[dir] = -(rho_liquide_-(rho_liquide_-rho_vapeur_)*vol_gaz/vol_NS)*gravite_[dir];
            }
          integrated_residu_[dir] += residu[dir] * fractional_dt;
        }
    }

  // Impression dans un fichier dedie :
  if (Process::je_suis_maitre())
    {
      int reset = (!reprise_) && (tstep_==0);
      SFichier fic=Ouvrir_fichier("_bilan_qdm.out",
                                  "tstep\ttime\tFs_theo\tFtot\trhov\ttauw\tS\tacceleration\tres\tcumul_res\tCorrFs\nConvection\tDiffusion\tPression\n# Forces have 3 components.",
                                  reset, 20/*prec*/);

      fic << tstep_<<" "<< current_time_<<" "
          << force_theo[0] << " "
          << force_theo[1] << " "
          << force_theo[2] << " ";
      fic << force_tot[0] << " "
          << force_tot[1] << " ";
      fic<< force_tot[2] << " "
         << rhov_moy[0] << " "
         << rhov_moy[1] << " "
         << rhov_moy[2] << " "
         << tauw[0] << " " ;
      fic	<< tauw[1] << " "
          << tauw[2] << " ";
      switch (direction_gravite_)
        {
        case 0:
          fic  << terme_source_acceleration_ << " "
               << "0. "
               << "0. ";
          break;
        case 1:
          fic  << "0. "
               << terme_source_acceleration_ << " "
               << "0. ";
          break;
        case 2:
          fic  << "0. "
               << "0. "
               << terme_source_acceleration_ << " ";
          break;
        default:
          fic  << terme_source_acceleration_ << " "
               << "0. "
               << "0. ";
        }
      fic << acc[0] << " "
          << acc[1] << " "
          << acc[2] << " "
          << residu[0] << " "
          << residu[1] << " " ;
      fic<< residu[2] << " "
         << integrated_residu_[0] << " "
         << integrated_residu_[1] << " "
         << integrated_residu_[2] << " ";
      fic	<< terme_source_correction_[0] << " "
          << terme_source_correction_[1] << " "
          << terme_source_correction_[2] << " ";
      // GAB, qdm
      // fic << terme_interfaces[0] << " "
      // << terme_interfaces[1] << " "
      // << terme_interfaces[2] << " ";
      fic << terme_convection_[0] << " "
          << terme_convection_[1] << " "
          << terme_convection_[2] << " ";
      fic << terme_diffusion_[0] << " "
          << terme_diffusion_[1] << " "
          << terme_diffusion_[2] << " ";
      fic << terme_pression_[0] << " "
          << terme_pression_[1] << " "
          << terme_pression_[2] << " ";
      //
#ifdef COMPLEMENT_ANTI_DEVIATION_RESIDU
      fic	<< moins_delta_Pwall_sur_h << " ";
#endif
      fic<< finl;
      fic.close();
    }
}

// Mettre rk_step = -1 si schema temps different de rk3.
// /!\ rk_step = 0 signifie qu'on est en RK3, mais n'indique pas a quelle ss pas de temps de RK3 on se place !!!
void IJK_FT_double::calculer_dv(const double timestep, const double time, const int rk_step)
{
  // GAB : initialisation. On initialise que pour rk_step<=0, mais on pourrai initialiser a chaque ss pdt
  // Avec le post-traitement que je fais de mes termes, je dois les re-initialiser a chaque ss pdt,
  //   si je ne les re-initialise pas alors je n'aurai pas leur valeur pour chaque avancement (apres je pourrai me passer de ce niveau de detail)
  //   ==> Ce qui est certain c'est que je ne dois pas re initialiser ICI mes rho_u_... puisqu'ils sont
  //       evalues EN DEHORS des boucles de RK3.
  if (rk_step<=0)
    {
      terme_convection_ = 0.;
      terme_pression_bis_ = 0.;
      terme_pression_ter_ = 0.;
      terme_diffusion_ = 0.;
      terme_interfaces_ = 0.;
      terme_interfaces_bf_mass_solver_ = 0.;
//      rho_u_euler_av_rho_mu_ind = 0.;
//      rho_u_euler_ap_rho_mu_ind = 0.;
    }
  // GAB
  // /!\ Valable que pour un domaine uniforme, j'ai vu des choses que je ne comprends pas la ou on defini volume aussi ...
  //     on est dans variable dz et on ne prends pas en compte k dans le calcul du volume...
  double volume_cell_uniforme = 1.;
  for (int i = 0; i < 3; i++)
    volume_cell_uniforme *= splitting_.get_grid_geometry().get_constant_delta(i);
  if (velocity_reset_)
    for (int dir=0; dir<3; dir++)
      velocity_[dir].data() = 0.; //Velocity reset for test

  static Stat_Counter_Id calcul_dv_counter_ = statistiques().new_counter(2, "maj vitesse : calcul derivee vitesse");
  statistiques().begin_count(calcul_dv_counter_);
  // Calcul d_velocity = convection
  if (!disable_convection_qdm_)
    {
      if (velocity_convection_op_.get_convection_op_option_rank() == non_conservative_simple)
        {
          velocity_convection_op_.calculer(velocity_[0], velocity_[1], velocity_[2],
                                           velocity_[0], velocity_[1], velocity_[2],
                                           d_velocity_[0], d_velocity_[1], d_velocity_[2]);
          // Multiplication par rho (on va rediviser a la fin)
          // (a partir de rho aux elements et dv aux faces)
          if (use_inv_rho_for_mass_solver_and_calculer_rho_v_)
            {
              // rho_face = 2*(rho_gauche*rho_droite)/(rho_gauche+rho_droite)
              //          = 1./ (1/2 * (1/rho_g + 1/rho_d))
              // 1/rho_face est donc la moyenne geometrique de inv_rho.
              calculer_rho_harmonic_v(rho_field_, d_velocity_, d_velocity_);
            }
          else
            {
              calculer_rho_v(rho_field_, d_velocity_, d_velocity_);
            }

        }
      else if (velocity_convection_op_.get_convection_op_option_rank() == non_conservative_rhou)
        {
          update_rho_v();
          // Non optimise car la methode calculer_avec_u_div_rhou inexistante.
          // Alors on initialise a 0 puis on fait ajouter :
          d_velocity_[0].data() = 0.;
          d_velocity_[1].data() = 0.;
          d_velocity_[2].data() = 0.;
          if (velocity_convection_op_.get_convection_op() != Nom("Centre4"))
            {
              velocity_convection_op_.ajouter_avec_u_div_rhou(rho_v_[0], rho_v_[1], rho_v_[2], // rhov_
                                                              velocity_[0], velocity_[1], velocity_[2],
                                                              d_velocity_[0], d_velocity_[1], d_velocity_[2],
                                                              div_rhou_);
            }
        }
      else if (velocity_convection_op_.get_convection_op_option_rank() == conservative)
        {
          update_rho_v();
          velocity_convection_op_.calculer(rho_v_[0], rho_v_[1], rho_v_[2],
                                           velocity_[0], velocity_[1], velocity_[2],
                                           d_velocity_[0], d_velocity_[1], d_velocity_[2]);
        }
      else
        {
          Cerr << "Unknown velocity convection type! " << finl;
          Process::exit();
        }
      post_.fill_op_conv();
      // GAB, qdm

      // a priori homogene a int_{volume_cellule} (d rho v / dt) pour le moment
      // on divise le terme_convection par volume_cellule
      // Il ne faudrai pas un mass solveur sur mon terme de convection... avant de le moyenner meme ?

      for (int dir = 0; dir < 3; dir++)
        {
          // if (rk_step==-1 || rk_step==0) // euler ou premier pdt de rk3
          // Pourquoi diviser par volume_cell_uniforme ?
          terme_convection_[dir] = calculer_v_moyen(d_velocity_[dir])/volume_cell_uniforme;

          if (test_etapes_et_bilan_)
            {
              terme_convection_mass_solver_[dir] = d_velocity_[dir];
              if(!disable_diphasique_)
                {
                  for (int k=0; k<d_velocity_[dir].nk(); k++)
                    {
                      mass_solver_with_rho(terme_convection_mass_solver_[dir], rho_field_, delta_z_local_, k);
                    }
                }
              else
                {
                  for (int k=0; k<d_velocity_[dir].nk(); k++)
                    for (int j=0; j<d_velocity_[dir].nj(); j++)
                      for (int i=0; i<d_velocity_[dir].ni(); i++)
                        terme_convection_mass_solver_[dir](i,j,k) = terme_convection_mass_solver_[dir](i,j,k) / (rho_liquide_ * volume_cell_uniforme);
                }

              terme_moyen_convection_mass_solver_[dir] = calculer_v_moyen(terme_convection_mass_solver_[dir]);
            }
          // else
          // terme_convection[dir] += calculer_v_moyen(d_velocity_[dir])/volume_cell_uniforme;
        }          //
    }
  else
    {
      d_velocity_[0].data() = 0.;
      d_velocity_[1].data() = 0.;
      d_velocity_[2].data() = 0.;
    }

  // Calcul diffusion
  if ((!diffusion_alternative_) && (!disable_diffusion_qdm_) )
    {
      velocity_diffusion_op_.set_nu(molecular_mu_);
      velocity_diffusion_op_.ajouter(velocity_[0], velocity_[1], velocity_[2],
                                     d_velocity_[0], d_velocity_[1], d_velocity_[2]);
      // GAB, qdm
      // a priori homogene a int_{volume_cellule} (d rho v / dt) pour le moment
      // mais on le divise par volume_cell_uniforme donc homogene a d rho v / dt maintenant
      // en diffusion "normale", on ajoute la diffusion plus tard
      for (int dir = 0; dir < 3; dir++)
        {
          // if (rk_step==-1 || rk_step==0) // revient a rk_step>0 // euler ou premier pdt de rk3
          terme_diffusion_[dir] = calculer_v_moyen(d_velocity_[dir])/volume_cell_uniforme - terme_convection_[dir];
          terme_diffusion_mass_solver_[dir] = d_velocity_[dir];
          // terme_diffusion_mass_solver_ contient la diffusion et la convection
          if (test_etapes_et_bilan_)
            {
              if (!disable_diphasique_)
                {
                  for (int k=0; k<terme_diffusion_mass_solver_[dir].nk(); k++)
                    {
                      mass_solver_with_rho(terme_diffusion_mass_solver_[dir], rho_field_, delta_z_local_, k);

                      // On retranche la convection
                      for (int j=0; j<terme_diffusion_mass_solver_[dir].nk(); j++ )
                        for (int i=0; i<terme_diffusion_mass_solver_[dir].ni(); i++)
                          terme_diffusion_mass_solver_[dir](i,j,k) = terme_diffusion_mass_solver_[dir](i,j,k) - terme_convection_mass_solver_[dir](i,j,k);
                    }
                }
              else
                {
                  // Dans le cas monophasique, rho_field_ vaut partout rho_liquide
                  for (int k=0; k<terme_diffusion_mass_solver_[dir].nk(); k++)
                    for (int j=0; j<terme_diffusion_mass_solver_[dir].nk(); j++ )
                      for (int i=0; i<terme_diffusion_mass_solver_[dir].ni(); i++)
                        {
                          terme_diffusion_mass_solver_[dir](i,j,k) = terme_diffusion_mass_solver_[dir](i,j,k) / (rho_liquide_*volume_cell_uniforme);
                          terme_diffusion_mass_solver_[dir](i,j,k) = terme_diffusion_mass_solver_[dir](i,j,k) - terme_convection_mass_solver_[dir](i,j,k);
                        }
                }
              // Moyenne sur le volume du domaine de simulation
              terme_moyen_diffusion_mass_solver_[dir] = calculer_v_moyen(terme_diffusion_mass_solver_[dir]);
              // else
              // terme_diffusion[dir] += calculer_v_moyen(d_velocity_[dir])/volume_cell_uniforme - terme_convection[dir];
            }
        }
      //
    }

  // GAB, qdm
  if (test_etapes_et_bilan_)
    for (int i = 0; i < 3; i++)
      {
        // Cout << "BF d_v_diff_et_conv" << finl;
        d_v_diff_et_conv_[i] = d_velocity_[i];
        // Cout << "AF d_v_diff_et_conv" << finl;
      }
  //

  // Calcul et ajout du grad(P^{n})*volume_cell (_entrelace?) :
  if (include_pressure_gradient_in_ustar_)
    {
      // pressure gradient requires the "left" value in all directions:
      pressure_.echange_espace_virtuel(1 /*, IJK_Field_double::EXCHANGE_GET_AT_LEFT_IJK*/);
#ifndef VARIABLE_DZ
      double volume = 1.;
      for (int i = 0; i < 3; i++)
        volume *= splitting_.get_grid_geometry().get_constant_delta(i);
#else
      Cerr << "This methods does not support variable DZ yet... Contact trust support. ";
      Process::exit();
      const double volume = get_channel_control_volume(dv, k, delta_z_local_);
#endif
      add_gradient_times_constant(pressure_, -volume /*constant*/,  d_velocity_[0], d_velocity_[1], d_velocity_[2]);
      Cout << "BF add_gradient_times_constant" << finl;
      // GAB, qdm
      // En utilisation normale de la pression, elle n'est pas ajoutee ici
      if (test_etapes_et_bilan_)
        {
          add_gradient_times_constant(pressure_, -volume,
                                      terme_pression_in_ustar_local_[0],
                                      terme_pression_in_ustar_local_[1],
                                      terme_pression_in_ustar_local_[2]);
          for (int dir_press = 0; dir_press < 3; dir_press++)
            terme_pression_in_ustar_[dir_press] = calculer_v_moyen(terme_pression_in_ustar_local_[dir_press]);
          Cout << "AF add_gradient_times_constant" << finl;
        }
      //
    }

  // Calcul du terme source aux interfaces pour l'ajouter a dv :
  // ATTENZIONE : Questa è una bugia. I termini di interfaccia vengono aggiunti direttamente a velocity_!
  if (!disable_diphasique_)
    {
      for (int dir = 0; dir < 3; dir++)
        {
          terme_source_interfaces_ft_[dir].data() = 0.;
          terme_repulsion_interfaces_ft_[dir].data() = 0.;
          terme_abs_repulsion_interfaces_ft_[dir].data() = 0.;
        }
      interfaces_.ajouter_terme_source_interfaces(terme_source_interfaces_ft_,
                                                  terme_repulsion_interfaces_ft_,
                                                  terme_abs_repulsion_interfaces_ft_);


      // Avant le solveur de masse, il faut un terme homogene a \int_vol {rho v }
      if (!disable_source_interf_)
        {
          for (int dir = 0; dir < 3; dir++)
            {
              if ((rk_step == -1) || (refuse_patch_conservation_QdM_RK3_source_interf_))
                {
                  // Avec le schema Euler, ou si on refuse le patch de conservation de la QdM en RK3 diphasique :
                  redistribute_from_splitting_ft_faces_[dir].redistribute_add(terme_source_interfaces_ft_[dir],
                                                                              d_velocity_[dir]);
                }
              else
                {
                  // On n'ajoute pas le terme source a d_velocity_...
                  // On le prendra en compte directement dans velocity_.
                }
            }

          for (int dir = 0; dir < 3; dir++)
            {
              redistribute_from_splitting_ft_faces_[dir].redistribute(terme_source_interfaces_ft_[dir],
                                                                      terme_source_interfaces_ns_[dir]);
              // GAB, qdm : les forces d'interface sont directement ajoutees a velocity... aucun effet sur d_velocity normalement
              //            Donc terme_interfaces_bf_mass_solver_bis est nul. C'EST A VERIFIER !!
              if (test_etapes_et_bilan_)
                {
                  terme_interfaces_bf_mass_solver_[dir] = calculer_v_moyen(terme_source_interfaces_ns_[dir]);
                  terme_interfaces_bf_mass_solver_bis_[dir] = calculer_v_moyen(d_velocity_[dir])/volume_cell_uniforme - terme_convection_[dir] - terme_diffusion_[dir];
                }
            }
          post_.fill_surface_force(terme_source_interfaces_ns_);
          // Computing force_tot (Attention, il faut le faire avant d'appliquer le solver mass a terme_source_interfaces_ns_) :
          compute_correction_for_momentum_balance(rk_step);
          for (int dir = 0; dir < 3; dir++)
            {
              // On est en RK3 et on utilise le patch de conservation de la QdM (comportement par defaut du RK3)
              // Utilisation directe du terme source interf pour l'ajouter a velocity_.
              // On ne le met plus dans d_velocity_ car ce n'est pas conservatif globalement... (test quand sigm et drho)
              if (post_.get_liste_post_instantanes().contient_("REPULSION_FT") ||
                  post_.get_liste_post_instantanes().contient_("CELL_REPULSION_FT")    )
                {
                  redistribute_from_splitting_ft_faces_[dir].redistribute(
                    terme_repulsion_interfaces_ft_[dir],
                    terme_repulsion_interfaces_ns_[dir]);
                }
              const int kmax = terme_source_interfaces_ns_[dir].nk();
              for (int k = 0; k < kmax; k++)
                {
                  // division par le produit (volume * rho_face)
                  if (use_inv_rho_for_mass_solver_and_calculer_rho_v_)
                    {
                      Cerr << "Je ne sais pas si inv_rho_field_ est a jour ici. A Verifier avant de l'activer." << finl;
                      Process::exit();
                      // mass_solver_with_inv_rho(terme_source_interfaces_ns_[di], inv_rho_field_, delta_z_local_, k);
                    }
                  else
                    {
                      // Division de terme_source_interfaces_ns_ par rho_field et par volume cellule
                      //cout << " code Delta : "<<terme_repulsion_interfaces_ns_[0](6,9,10);
                      //cout << " code Delta : "<<backup_terme_source_interfaces_ns_[0](6,9,10);
                      //cout << " code Delta : "<<post_.get_rho_Ssigma()[2](6,11,11);
                      //cout << endl;
                      mass_solver_with_rho(terme_source_interfaces_ns_[dir], rho_field_, delta_z_local_, k);
                      //cout << " code Nabla : "<<terme_repulsion_interfaces_ns_[0](6,9,10);
                      //cout << " code Nabla : "<<backup_terme_source_interfaces_ns_[0](6,9,10);
                      //cout << " code Nabla : "<<post_.get_rho_Ssigma()[2](6,11,11);
                      //cout << endl;
                      //if (k==10)
                      //  {
                      //    cout << "code 2211 " << rho_field_(5,9,10) << " & " << rho_field_(6,9,10) << " & " << rho_field_(7,9,10)<<endl;
                      //    cout << "code Z,611 " << post_.get_rho_Ssigma()[2](6,11,11);
                      //    cout << ", code Z,611 " << backup_terme_source_interfaces_ns_[2](6,11,11);
                      //    cout << ", code Z,611 " << terme_source_interfaces_ns_[2](6,11,11)<<endl;
                      //    cout << "code X,6910 " << post_.get_rho_Ssigma()[0](6,9,10);
                      //    cout << ", code X,6910 " << backup_terme_source_interfaces_ns_[0](6,9,10);
                      //    cout << ", code X,6910 " << terme_source_interfaces_ns_[0](6,9,10)<<endl;
                      //  }
                      if (post_.get_liste_post_instantanes().contient_("REPULSION_FT") ||
                          post_.get_liste_post_instantanes().contient_("CELL_REPULSION_FT")    )
                        {
                          // Division de terme_repulsion_interfaces_ns_ par rho_field_ et par volume cellule
                          mass_solver_with_rho(terme_repulsion_interfaces_ns_[dir], rho_field_, delta_z_local_, k);
                          //terme_repulsion_interfaces_ns_[0](6,9,10)=50;
                        }
                      // Egalite aux dimensions :
                      // [terme_source_interfaces_ns_]=[terme_repulsion_interfaces_ns_]=[du/dt / Vcell] = m/s^2/m^3
                    }
                  if ((!refuse_patch_conservation_QdM_RK3_source_interf_) && (rk_step>=0) )
                    {
                      // puis
                      // comme euler_explicit_update mais avec un pas de temps partiel :
                      const double delta_t = compute_fractionnal_timestep_rk3(timestep_ /* total*/, rk_step);
                      const int imax = terme_source_interfaces_ns_[dir].ni();
                      const int jmax = terme_source_interfaces_ns_[dir].nj();
                      for (int j = 0; j < jmax; j++)
                        {
                          for (int i = 0; i < imax; i++)
                            {
                              double x = terme_source_interfaces_ns_[dir](i,j,k);
                              velocity_[dir](i,j,k) += x * delta_t;
                            }
                        }
                    }
                  // On est dans une boucle sur les directions la, c ok
                  // GAB, qdm  ATTENTION on ne va ici que si on est en rk3
                  if (test_etapes_et_bilan_)
                    {
                      terme_interfaces_af_mass_solver_[dir] = calculer_v_moyen(terme_source_interfaces_ns_[dir]);
                    }
                }

              //const int kmax = terme_source_interfaces_ns_[0].nk();
              //for (int ii = 0; ii < kmax; ++ii)
              //  {
              //    for (int jj = 0; jj < kmax; ++jj)
              //      {
              //        for (int kk = 0; kk < kmax; ++kk)
              //          {
              //            if (backup_terme_source_interfaces_ns_[dir](ii,jj,kk) != terme_source_interfaces_ns_[dir](ii,jj,kk) )
              //              {
              //                cout << "change after mass solver for : dir,i,j,k = "<<dir<<", "<<ii<<", "<<jj<<", "<<kk<<endl;
              //              }
              //          }
              //      }
              //  }
            }
        }
    }

  // On laisse l'ecriture de ce fichier de sortie A L'INTERIEUR de calculer_dv car on souhaite relever
  // la valeur des differents termes A CHAQUE sous-pas de temps du schema RK3. On peut en revanche,

  fill_variable_source_and_potential_phi(time);


  for (int dir = 0; dir < 3; dir++)
    {
      // GAB question : pour quoi on cree dv, qu'est ce qui empeche de travailler sur d_velocity ???
      //                -> Est ce que c'est uniquement parce qu'on va faire un certain nb d'operations dessus
      //                   et que manipuler dv est plus leger que de manipuler d_velocity ?
      IJK_Field_double& dv = d_velocity_[dir];
      const int kmax = d_velocity_[dir].nk();
      const int ni = dv.ni();
      const int nj = dv.nj();
      // terme source acceleration homogene a une force volumique (gradient de pression uniforme)
      // Si la correction est activee, on oppose la force_moy
      double force_volumique = correction_force_[dir]*terme_source_correction_[dir];
      // if (dir == DIRECTION_I)
      //   {
      //     force_volumique += terme_source_acceleration_;
      //   }
      // GAB, rotation
      if (dir == direction_gravite_)
        {
          force_volumique += terme_source_acceleration_;
        }

#ifndef VARIABLE_DZ
      double volume = 1.;
      for (int i = 0; i < 3; i++)
        volume *= splitting_.get_grid_geometry().get_constant_delta(i);

      // GAB, qdm
      // dans d_velocity_moyen on a la contrib de interfaces, forces ajoutees
      terme_interfaces_conv_diff_mass_solver_[dir] = calculer_v_moyen(d_velocity_[dir]);

      Cerr << "disable_diffusion_qdm_ : "<< disable_diffusion_qdm_ << finl;
      Cerr << "diffusion_alternative_ : "<< diffusion_alternative_ << finl;
      Cerr << "type_velocity_diffusion_form : "<< velocity_diffusion_op_.get_diffusion_op_option() << finl;
      for (int k = 0; k < kmax; k++)
        {
          // #else
          //       for (int k = 0; k < kmax; k++)
          //         {
          //           const double volume = get_channel_control_volume(dv, k, delta_z_local_);
          //         }
          const double f = force_volumique * volume;
          if ((expression_variable_source_[0] != "??") ||
              (expression_variable_source_[1] != "??") ||
              (expression_variable_source_[2] != "??") ||
              (expression_potential_phi_ != "??"))
            {
              for (int j = 0; j < nj; j++)
                for (int i = 0; i < ni; i++)
                  dv(i,j,k) += facteur_variable_source_*variable_source_[dir](i,j,k) * volume + f;
            }
          else
            {
              for (int j = 0; j < nj; j++)
                for (int i = 0; i < ni; i++)
                  dv(i,j,k) += f;
            }

          if (use_inv_rho_for_mass_solver_and_calculer_rho_v_)
            mass_solver_with_inv_rho(d_velocity_[dir], inv_rho_field_, delta_z_local_, k);
          else
            mass_solver_with_rho(d_velocity_[dir], rho_field_, delta_z_local_, k);


          // Terme source gravitaire en version simple "rho_g":
          // GAB : question a Guillaume : l 3235 -> "IJK_Field_double& dv = d_velocity_[dir];" ,  MAIS C'EST APRES
          // LA LIGNE 2865 ou d_velocity subi un calculer_rho_v.
          //       donc dv est homogene a rho*v et donc a rho*g. Or 8 lignes plus bas on a dv(i,j,k) += g; !!!
          // il faut appliquer un mass_solver_with_rho a dv pour que l'operation soit homogene d'apres gr...
          // a moins que le mass_solver soit applique dans une des methodes...
          // mass_solver_with... est applique sur d_velocity_ juste au dessus de ce long commentaire. Donc d_velocity
          // est bien homogene a g, et ainsi (jespere) dv est bien homogene a g.
          if (interfaces_.is_terme_gravite_rhog())
            {
              const double g = gravite_[dir];
              for (int j = 0; j < nj; j++)
                for (int i = 0; i < ni; i++)
                  dv(i,j,k) += g;  // c'est rho*g qu'il faut pour GR, 18.08.2021
            }

          if ((diffusion_alternative_) && (!disable_diffusion_qdm_))
            {
              velocity_diffusion_op_.set_nu(unit_);
              velocity_diffusion_op_.ajouter(velocity_[0], velocity_[1], velocity_[2],
                                             laplacien_velocity_[0], laplacien_velocity_[1], laplacien_velocity_[2]);
              for (int dir2 = 0; dir2 < 3; dir2++)
                {
                  const int kmax2 = d_velocity_[dir2].nk();
                  const int imax = d_velocity_[dir2].ni();
                  const int jmax = d_velocity_[dir2].nj();
                  for (int k2 = 0; k2 < kmax2; k2++)
                    for (int j = 0; j < jmax; j++)
                      for (int i = 0; i < imax; i++)
                        {
                          // double laplacien_u = laplacien_velocity_[dir2](i,j,k2);
                          // GAB, c'est assez dangereux de nommer nu une viscosite dynamique...
                          // const double nu = molecular_mu_(i,j,k2) ; // On stocke nu dans mu dans ce cas.
                          // GAB, qdm
                          if (test_etapes_et_bilan_)
                            {
                              terme_diffusion_local_[dir2](i,j,k2) = molecular_mu_(i,j,k2) * laplacien_velocity_[dir2](i,j,k2);
                              // Cerr << "terme diffusion local" << terme_diffusion_local[dir2](i,j,k2) << finl;
                              d_velocity_[dir2](i,j,k2) += terme_diffusion_local_[dir2](i,j,k2);
                              // d_velocity_[dir2](i,j,k2) += nu * laplacien_u;
                            }
                        }
                  // GAB, qdm
                  d_velocity_[dir2].echange_espace_virtuel(d_velocity_[dir2].ghost());
                  //
                }
            }

        }
      // GAB, TODO bilan qdm
      // d_velocity_moyen_ap_mass_solver[dir] = calculer_v_moyen(d_velocity_[dir]);
      // dans d_velocity_conv_et_diff_moy on a que la contrib de convection et de diffusion
      // d_velocity_conv_et_diff_moy[dir] = calculer_v_moyen(d_velocity_conv_et_diff[dir])

#endif
      // For the addition of the external forces to d_velocity_
      compute_add_external_forces(dir);
    } // end of loop [dir].
  ///////////////////////////////////////////////////////
  // GAB, THI
  // MODIF GAB.. /!\ tstep : le numero d'iteration; timestep : le pas de temps, le dt
  // est-on au moment ou velocity_ contient la vitesse ou est-on au moment ou velocity_ contient la qdm ?
  if (forcage_.get_type_forcage() > 0)// && (rk_step==-1 || rk_step==0))
    {
      if (rk_step==-1)
        compute_add_THI_force_sur_d_velocity(velocity_, tstep_, timestep_, time,
                                             d_velocity_.get_splitting(),
                                             forcage_.get_facteur_forcage());//, rk_step);
      else
        {
          const double intermediate_dt = compute_fractionnal_timestep_rk3(timestep_, rk_step);
          compute_add_THI_force_sur_d_velocity(velocity_, tstep_, intermediate_dt, time,
                                               d_velocity_.get_splitting(),
                                               forcage_.get_facteur_forcage());//, rk_step);
        }
    }
  // verifier si mon terme de thi est bon en integrale
  ///////////////////////////////////////////////////////
  Cout << "G bilan qdm " << finl;
  if (test_etapes_et_bilan_)
    write_check_etapes_et_termes(rk_step);

  //


  // Il est important de s'assurer a la fin que la derivee de la vitesse soit a zero sur les parois:
  if (!splitting_.get_grid_geometry().get_periodic_flag(DIRECTION_K))
    force_zero_on_walls(d_velocity_[2]);

  statistiques().end_count(calcul_dv_counter_);
}

void IJK_FT_double::compute_add_external_forces(const int dir)
{
  if (coef_immobilisation_ > 1e-16)
    {
      // maxValue n'avait pas le mp_max
      double integration_time = max_ijk(post_.integrated_timescale());
      integration_time=std::max(1.,integration_time); // if integrated_timescale is missing, the value of integration will be -1.e30;
      // 											The trick is to set it to 1, as it is in fact numbers ot timesteps stored (see dirty code)
      interfaces_.compute_external_forces_(force_rappel_ft_, force_rappel_, velocity_,interfaces_.I(),interfaces_.I_ft(),
                                           coef_immobilisation_, tstep_, current_time_,
                                           coef_ammortissement_, coef_rayon_force_rappel_,
                                           integration_time, coef_mean_force_, coef_force_time_n_);

      if (interfaces_.get_forcing_method())
        {
          // Si Parser
          const int kmax = d_velocity_[dir].nk();
          const int ni = d_velocity_[dir].ni();
          const int nj = d_velocity_[dir].nj();
          for (int k = 0; k < kmax; k++)
            for (int j = 0; j < nj; j++)
              for (int i = 0; i < ni; i++)
                d_velocity_[dir](i,j,k) += force_rappel_[dir](i,j,k);
        }
      else
        {
          // Si la force de rappel est calculee sans le parser (par la color function), elle est evaluee sur le FT.
          redistribute_from_splitting_ft_faces_[dir].redistribute_add(
            force_rappel_ft_[dir], d_velocity_[dir]);
          // Je pense qu'il faut la redistribuer pour le calcul de la force de rappel dans le terme expression_derivee_acceleration_
          redistribute_from_splitting_ft_faces_[dir].redistribute(
            force_rappel_ft_[dir], force_rappel_[dir]);
        }
    }
  return;
}

// -----------------------------------------------------------------------------------
//  FORCAGE EXTERIEUR, DEFINI DANS L'ESPACE SPECTRAL
// GAB, THI /!\ REMPLACER PAR compute_add_THI_force_sur_d_velocity
void IJK_FT_double::compute_add_THI_force(const FixedVector<IJK_Field_double, 3>& vitesse,
                                          const int time_iteration,
                                          const double dt, //tstep, /!\ ce dt est faux, je ne sais pas pk mais en comparant sa valeur avec celle du dt_ev, je vois que c'est faux
                                          const double current_time,
                                          const IJK_Splitting& my_splitting
                                          // const int rk_step
                                         )
{
  statistiques().begin_count(m2_counter_);
  if (forcage_.get_forced_advection()==-1)
    {
      ArrOfDouble mean_u_liq;
      mean_u_liq.resize_array(3);
      for (int dir=0; dir<3; dir++)
        mean_u_liq[dir] = calculer_moyenne_de_phase_liq(vitesse[dir]);
      forcage_.update_advection_velocity(mean_u_liq);
    }
  if (forcage_.get_forced_advection()!=0)
    {
      Cout << "forced_advection" << forcage_.get_forced_advection() << finl;
      Cout << "BF : update_advection_length" << finl;
      forcage_.update_advection_length(dt);
      Cout << "AF : update_advection_length" << finl;
    }
  forcage_.compute_THI_force(time_iteration,dt,current_time,splitting_);

  statistiques().end_count(m2_counter_);

  statistiques().begin_count(m3_counter_);

  const FixedVector<IJK_Field_double, 3>& force = forcage_.get_force_ph2();
  for(int dir=0; dir<3; dir++)
    {
      // d_velocity_ est deja decoupe sur les differents procs; donc ni, nj, nk =! nb elem.
      const int kmax = d_velocity_[dir].nk();
      const int ni = d_velocity_[dir].ni();
      const int nj = d_velocity_[dir].nj();

      for (int k = 0; k < kmax+1; k++)
        for (int j = 0; j < nj+1; j++)
          for (int i = 0; i < ni+1; i++)
            {
              // GAB error : No continuous phase found ...
              // appeler mass_solver_with_rho si je veux div par rho*V (tant que je garde f localise aux faces)
              // double cell_mass = rho_field_(i,j,k)*my_geom.get_constant_delta(DIRECTION_I)*my_geom.get_constant_delta(DIRECTION_J)*my_geom.get_constant_delta(DIRECTION_K);
              const double inv_cell_mass = 1.; //my_geom.get_constant_delta(DIRECTION_I)*my_geom.get_constant_delta(DIRECTION_J)*my_geom.get_constant_delta(DIRECTION_K);
              //d_velocity_[dir](i,j,k) += force[dir](i,j,k)*inv_cell_mass;
              velocity_[dir](i,j,k) += force[dir](i,j,k)*inv_cell_mass*dt;
            }
    }
  statistiques().end_count(m3_counter_);
}

void IJK_FT_double::compute_add_THI_force_sur_d_velocity(const FixedVector<IJK_Field_double, 3>& vitesse,
                                                         const int time_iteration,
                                                         const double dt, //tstep,  /!\ ce dt est faux, je ne sais pas pk mais en comparant sa valeur avec celle du dt_ev, je vois que c'est faux
                                                         const double current_time,
                                                         const IJK_Splitting& my_splitting,
                                                         const int facteur
                                                         // const int rk_step
                                                        )
{
  statistiques().begin_count(m2_counter_);
  if (forcage_.get_forced_advection()==-1)
    {
      /* Advection du champ de force par mean{u_l}^l */
      ArrOfDouble mean_u_liq;
      mean_u_liq.resize_array(3);
      for (int dir=0; dir<3; dir++)
        {
          // (22.02.22) FAUX : mean_u_liq[dir] contient alors alpha_liq*vitesse_liq[dir]
          // mean_u_liq[dir] = calculer_moyenne_de_phase_liq(vitesse_liq[dir]);
          // (22.02.22) JUSTE : mean_u liq[dir] contient bien vitesse_liq[dir]
          mean_u_liq[dir] = calculer_true_moyenne_de_phase_liq(vitesse[dir]);
        }
      forcage_.update_advection_velocity(mean_u_liq);
    }
  if (forcage_.get_forced_advection()!=0)
    {
      /* Advection du champ de force par advection_velocity_, donnee du jdd */
      forcage_.update_advection_length(dt);
    }
  forcage_.compute_THI_force(time_iteration,dt,current_time,splitting_);
  statistiques().end_count(m2_counter_);

  statistiques().begin_count(m3_counter_);

  const FixedVector<IJK_Field_double, 3>& force = forcage_.get_force_ph2();

  for(int dir=0; dir<3; dir++)
    {
      // d_velocity_ est deja decoupe sur les differents procs; donc ni, nj, nk =! nb elem.
      const int kmax = d_velocity_[dir].nk();
      const int ni = d_velocity_[dir].ni();
      const int nj = d_velocity_[dir].nj();

      // La meme force appliquee partout
      if (facteur==0)
        {
          for (int k = 0; k < kmax; k++)
            for (int j = 0; j < nj; j++)
              for (int i = 0; i < ni; i++)
                {
                  // appeler mass_solver_with_rho si je veux div par rho*V (tant que je garde f localise aux faces)
                  // double cell_mass = rho_field_(i,j,k)*my_geom.get_constant_delta(DIRECTION_I)*my_geom.get_constant_delta(DIRECTION_J)*my_geom.get_constant_delta(DIRECTION_K);
                  //d_velocity_[dir](i,j,k) += force[dir](i,j,k)*inv_cell_mass;
                  d_velocity_[dir](i,j,k) += force[dir](i,j,k);
                }
        }
      // Force active uniquement dans le liquide. Nulle dans le gaz
      else if (facteur==1)
        {
          for (int k = 0; k < kmax; k++)
            for (int j = 0; j < nj; j++)
              for (int i = 0; i < ni; i++)
                {
                  // appeler mass_solver_with_rho si je veux div par rho*V (tant que je garde f localise aux faces)
                  // double cell_mass = rho_field_(i,j,k)*my_geom.get_constant_delta(DIRECTION_I)*my_geom.get_constant_delta(DIRECTION_J)*my_geom.get_constant_delta(DIRECTION_K);
                  // indicatrice == 0 dans le gaz et 1 dans le liquide. ATTN : ~0.5 aux interfaces !!!
                  d_velocity_[dir](i,j,k) += force[dir](i,j,k)*interfaces_.I(i,j,k);
                }
        }
      // Force ponderee par la masse volumique de la phase. Attention, il faut rester homogene a une vitesse
      else if (facteur==2)
        {
          // S'assurer qu'on a bien calcule rho_moyen_
          for (int k = 0; k < kmax; k++)
            for (int j = 0; j < nj; j++)
              for (int i = 0; i < ni; i++)
                {
                  // appeler mass_solver_with_rho si je veux div par rho*V (tant que je garde f localise aux faces)
                  // double cell_mass = rho_field_(i,j,k)*my_geom.get_constant_delta(DIRECTION_I)*my_geom.get_constant_delta(DIRECTION_J)*my_geom.get_constant_delta(DIRECTION_K);
                  d_velocity_[dir](i,j,k) += force[dir](i,j,k)*rho_field_(i,j,k)/rho_moyen_;
                }
        }
      d_velocity_[dir].echange_espace_virtuel(d_velocity_[dir].ghost());
    }
  statistiques().end_count(m3_counter_);
  Cout << "end of from_spect_to_phys_opti2_advection" << finl;
}
// -----------------------------------------------------------------------------------

void IJK_FT_double::euler_time_step(ArrOfDouble& var_volume_par_bulle)
{
  static Stat_Counter_Id euler_rk3_counter_ = statistiques().new_counter(2, "Mise a jour de la vitesse");
  statistiques().begin_count(euler_rk3_counter_);
  if (thermals_.size())
    {
      if (boundary_conditions_.get_correction_conserv_qdm()==2)
        {
          update_rho_v();
          rho_field_.echange_espace_virtuel(rho_field_.ghost());
          update_v_ghost_from_rho_v();
        }
      else
        {
          // Protection to make sure that even without the activation of the flag check_divergence_, the EV of velocity is correctly field.
          // This protection MAY be necessary if convection uses ghost velocity (but I'm not sure it actually does)
          velocity_[0].echange_espace_virtuel(2);
          velocity_[1].echange_espace_virtuel(2);
          velocity_[2].echange_espace_virtuel(2);
        }


    }

  /*
   * TODO: Change this block with DERIV CLASS IJK_Thermal
   */
  for (auto& itr : thermique_)
    itr.euler_time_step(timestep_);

  for (auto& itr : energie_)
    itr.euler_time_step(velocity_);

  thermals_.euler_time_step(timestep_);

  if (!frozen_velocity_)
    {
      if (boundary_conditions_.get_correction_conserv_qdm()==2)
        {
          update_rho_v();
          rho_field_.echange_espace_virtuel(rho_field_.ghost());
          update_v_ghost_from_rho_v();
        }
      else
        {
          velocity_[0].echange_espace_virtuel(2);
          velocity_[1].echange_espace_virtuel(2);
          velocity_[2].echange_espace_virtuel(2);
        }

      // GAB, qdm
      if (test_etapes_et_bilan_)
        {
          calculer_rho_v(rho_field_,velocity_,rho_u_euler_av_prediction_champ_);
          for (int dir = 0; dir<3; dir++)
            rho_u_euler_av_prediction_[dir] = calculer_v_moyen(rho_u_euler_av_prediction_champ_[dir]);
        }
      // GAB, remarque : calculer dv calcule dv, MAIS NE L'APPLIQUE PAS au champ de vitesse !!!
      //                 l'increment de vitesse est ajoute au champ de vitesse avec euler_explicit_update
      calculer_dv(timestep_, current_time_, -1 /*rk_step = -1 pour sch euler... */);
      // GAB, qdm calculer_dv ne fait que l'etape de prediction)
      if (test_etapes_et_bilan_)
        {
          calculer_rho_v(rho_field_,d_velocity_,rho_du_euler_ap_prediction_champ_);
          for (int dir = 0; dir<3; dir++)
            rho_du_euler_ap_prediction_[dir] = calculer_v_moyen(rho_du_euler_ap_prediction_champ_[dir]);
        }
#ifdef PROJECTION_DE_LINCREMENT_DV
      // ajout du gradient de pression a dv
      if (!disable_solveur_poisson_)
        {
          if (!include_pressure_gradient_in_ustar_)
            {
              pressure_projection_with_rho(rho_field_, d_velocity_[0], d_velocity_[1],  d_velocity_[2],
                                           pressure_, 1.,  pressure_rhs_, check_divergence_, poisson_solver_);
              // GAB --> C'est plutot ici que l(on ajoute le terme_pression !!)
            }
          else
            {
              pressure_projection_with_rho(rho_field_, d_velocity_[0], d_velocity_[1],  d_velocity_[2],
                                           d_pressure_, 1.,  pressure_rhs_, check_divergence_, poisson_solver_);

              // Then update the pressure field :
              const int kmax = pressure_.nk();
              for (int k = 0; k < kmax; k++)
                euler_explicit_update(d_pressure_, pressure_, k);
            }
        }

#else
#endif
      // Mise a jour du champ de vitesse (etape de projection et de prediction)
      for (int dir = 0; dir < 3; dir++)
        {
          const int kmax = d_velocity_[dir].nk();
          for (int k = 0; k < kmax; k++)
            {
              // GAB, question : d_velocity est issu de calculer_dv. Il manque pas l'appliquation de l'operateur de divergence avant d'appliquer d_velocity e velocity ?
              //                 il manque au moins la multiplication par les surfaces je pense -> NON, lis bien les etapes de convection et de diffusion.
              euler_explicit_update(d_velocity_[dir], velocity_[dir], k);
            }
        }

      // GAB, qdm : cree le rho_n * v_n+1
      if (test_etapes_et_bilan_)
        {
          calculer_rho_v(rho_field_,d_velocity_,rho_du_euler_ap_projection_champ_);
          for (int dir=0; dir<3; dir++)
            rho_du_euler_ap_projection_[dir] = calculer_v_moyen(rho_du_euler_ap_projection_champ_[dir]);
        }

      if (test_etapes_et_bilan_)
        {
          calculer_rho_v(rho_field_,velocity_,rho_u_euler_ap_projection_champ_);
          for (int dir=0; dir<3; dir++)
            rho_u_euler_ap_projection_[dir] = calculer_v_moyen(rho_u_euler_ap_projection_champ_[dir]);
        }

      // Conditions en entree
      if (vitesse_entree_ > -1e20)
        force_entry_velocity(velocity_[0], velocity_[1], velocity_[2], vitesse_entree_);

      // Forcage de la vitesse en amont de la bulle :
      if (vitesse_upstream_ > -1e20)
        {
          if (expression_vitesse_upstream_ != "??")
            {
              std::string expr(expression_vitesse_upstream_);
              Parser parser;
              parser.setString(expr);
              parser.setNbVar((int) 1);
              parser.addVar("t");
              parser.parseString();
              parser.setVar((int) 0, (*this).current_time_ - modified_time_ini_);
              vitesse_upstream_ = parser.eval();
            }
          Cerr << "Force upstream velocity" << finl;
          force_upstream_velocity(velocity_[0], velocity_[1], velocity_[2],
                                  vitesse_upstream_, interfaces_, nb_diam_upstream_,
                                  upstream_dir_, get_direction_gravite(), upstream_stencil_);
        }
    } // end of if ! frozen_velocity
// static Stat_Counter_Id projection_counter_ = statistiques().new_counter(0, "projection");
#ifdef PROJECTION_DE_LINCREMENT_DV
  if (0)
#else
  if (!disable_solveur_poisson_)
#endif
    {

      //  statistiques().begin_count(projection_counter_);
      if (include_pressure_gradient_in_ustar_)
        {
          Cerr << "Methode incrementale pour le grad(P)" << finl;
          Cerr << "Initialisation du d_pressure_ conservee depuis le pas de temps precedent... " << finl;
          Cerr << "Ce n'est probablement pas optimal. QQ idees dans les sources si divergence . " << finl;
          // Que vaut d_pressure ?
          // Important car c'est l'initialisation du solveur...

          // 1. raz :
          // d_pressure_.data() = 0.; // raz...
          // 2. dp = - timestep_ * (potentiel_elem - delta_rho * phi) * u . grad(I)
          //                       (sigma_ * courbure)                  on a ustar a dispo, par u^n.
          // 3. dp = - timestep_ * u . grad(P^n)
          //                       ici, on a ustar dispo, plus u^n.
          //                       Si on veut tester avec u^n (c mieux je pense), il faut initialiser dp
          //                       juste avant l'euler_explicit_update

          if (use_inv_rho_in_poisson_solver_)
            {

              pressure_projection_with_inv_rho(inv_rho_field_, velocity_[0], velocity_[1],  velocity_[2], d_pressure_, timestep_,
                                               pressure_rhs_, check_divergence_, poisson_solver_);

            }
          else
            {
#ifdef PROJECTION_DE_LINCREMENT_DV
              // On l'a fait avant pour etre sur qu'elle soit bien dans la derivee stockee...
#else

              pressure_projection_with_rho(rho_field_, velocity_[0], velocity_[1],  velocity_[2], d_pressure_, timestep_,
                                           pressure_rhs_, check_divergence_, poisson_solver_);
#endif

            }

          // Mise a jour de la pression :
          for (int dir = 0; dir < 3; dir++)
            {
              const int kmax = pressure_.nk();
              for (int k = 0; k < kmax; k++)
                euler_explicit_update(d_pressure_, pressure_, k);
            }

          Cerr << " Un exit pour voir avec gdb... " << finl;
          Cerr << " Si ca fonctionne, faire le meme en RK3... " << finl;
          // Process::exit();
        }
      else
        {
          if (use_inv_rho_in_poisson_solver_)
            {

              pressure_projection_with_inv_rho(inv_rho_field_, velocity_[0], velocity_[1],  velocity_[2], pressure_, timestep_,
                                               pressure_rhs_, check_divergence_, poisson_solver_);

            }
          else
            {
#ifdef PROJECTION_DE_LINCREMENT_DV
#else

              pressure_projection_with_rho(rho_field_, velocity_[0], velocity_[1],  velocity_[2], pressure_, timestep_,
                                           pressure_rhs_, check_divergence_, poisson_solver_);

#endif
            }
        }
      if (test_etapes_et_bilan_)
        {
          // GAB, qdm : recuperons le temre de pression (1/rho * grad(p)) si on fait le bilan en u (ca a du sens meme?)
          //                                                     grap(p) si on fait le bilan de qdm
          // terme_pression_bis = calculer_inv_rho_grad_p_moyen(rho_field_, pressure_);
          terme_pression_bis_ = calculer_grad_p_moyen(pressure_);
          // GAB, qdm : recuperons le terme de pression (1/rho * grad(p))
          terme_pression_ter_ = calculer_grad_p_over_rho_moyen(pressure_);
          pression_ap_proj_ = calculer_v_moyen(pressure_);
        }
      //statistiques().end_count(projection_counter_);
    }

  Cerr << "Copy pressure on extended field for probes" << finl;
  copy_field_values(pressure_ghost_cells_, pressure_);

  if (Process::je_suis_maitre())
    {
      Cout << "Timings diff=" << statistiques().last_time(diffusion_counter_)
           << " conv=" << statistiques().last_time(convection_counter_);
      Cout << " src=" << statistiques().last_time(source_counter_)
           << finl;
    }
  statistiques().end_count(euler_rk3_counter_);
}

// Perform one sub-step of rk3 for FT algorithm, called 3 times per time step.
// rk_step = 0, 1 or 2
// total_timestep = not the fractionnal timestep !
void IJK_FT_double::rk3_sub_step(const int rk_step, const double total_timestep,
                                 const double fractionnal_timestep, const double time )
{
  assert(rk_step>=0 && rk_step<3);
  static Stat_Counter_Id euler_rk3_counter_ = statistiques().new_counter(2, "Mise a jour de la vitesse");
  statistiques().begin_count(euler_rk3_counter_);

  /*
   * TODO: Change this block with DERIV CLASS IJK_Thermal
   */
  for (auto& itr : thermique_)
    {
      itr.rk3_sub_step(rk_step, total_timestep, time);
    }
  for (auto&& itr = energie_.begin(); itr != energie_.begin(); ++itr)
    {
      // curseur->rk3_sub_step(rk_step, total_timestep, time);
      // ++curseur;
      Cerr << "Le schema RK3 n est pas implemente avec des champs d energie" << finl;
      Process::exit();
    }

  thermals_.rk3_sub_step(rk_step, total_timestep, time);

  if (!frozen_velocity_)
    {
      if (boundary_conditions_.get_correction_conserv_qdm()==2)
        {
          update_rho_v();
          rho_field_.echange_espace_virtuel(rho_field_.ghost());
          update_v_ghost_from_rho_v();
        }
      else
        {
          velocity_[0].echange_espace_virtuel(2);
          velocity_[1].echange_espace_virtuel(2);
          velocity_[2].echange_espace_virtuel(2);
        }
      // GAB TODO : voir dans euler_explicite ce qu'on a dit qu'on ferai pour voir
      // si le calculer_dv s'est bien passe
      Cout << "rk3ss: rk_step " << rk_step << finl;
      calculer_dv(total_timestep, time, rk_step);
      //
#ifdef PROJECTION_DE_LINCREMENT_DV
      // ajout du gradient de pression a dv
      if (!disable_solveur_poisson_)
        {
          if (include_pressure_gradient_in_ustar_)
            {
              pressure_projection_with_rho(rho_field_, d_velocity_[0], d_velocity_[1],  d_velocity_[2],
                                           d_pressure_,1. ,  pressure_rhs_, check_divergence_, poisson_solver_);
            }
          else
            {
              pressure_projection_with_rho(rho_field_, d_velocity_[0], d_velocity_[1],  d_velocity_[2],
                                           pressure_,1. ,  pressure_rhs_, check_divergence_, poisson_solver_);
            }
        }
#else
#endif

      // Mise a jour du champ de vitesse (etape de projection GAB, 28/06/21 : c'est la prediction, pas la projection non ?)
      for (int dir = 0; dir < 3; dir++)
        {
          const int kmax = d_velocity_[dir].nk();
          for (int k = 0; k < kmax; k++)
            {
              runge_kutta3_update(d_velocity_[dir], RK3_F_velocity_[dir], velocity_[dir], rk_step, k, total_timestep);
              // GAB, correction qdm a posteriori
            }
        }

#ifdef PROJECTION_DE_LINCREMENT_DV
      // Mise a jour du champ de pression
      if ((!disable_solveur_poisson_) && (include_pressure_gradient_in_ustar_))
        {
          const int kmax = pressure_.nk();
          for (int k = 0; k < kmax; k++)
            {
              runge_kutta3_update(d_pressure_ /* increment */,
                                  RK3_F_pressure_ /* intermediate storage */,
                                  pressure_ /* variable to update */, rk_step, k, total_timestep);
            }
        }
#else
#endif

      // Conditions en entree
      if (vitesse_entree_ > -1e20)
        force_entry_velocity(velocity_[0], velocity_[1], velocity_[2], vitesse_entree_);


      // Forcage de la vitesse en amont de la bulle :
      if (vitesse_upstream_ > -1e20)
        force_upstream_velocity(velocity_[0], velocity_[1], velocity_[2],
                                vitesse_upstream_, interfaces_, nb_diam_upstream_,
                                upstream_dir_, get_direction_gravite(), upstream_stencil_);

    } // end of if ! frozen_velocity
//static Stat_Counter_Id projection_counter_ = statistiques().new_counter(0, "projection");
#ifdef PROJECTION_DE_LINCREMENT_DV
  if (0)
#else
  if (!disable_solveur_poisson_)
#endif
    {
      //statistiques().begin_count(projection_counter_);
      if (include_pressure_gradient_in_ustar_)
        {
          Cerr << "L'option include_pressure_gradient_in_ustar n'est pas encore implementee en RK3." << finl;

          Process::exit();
        }

      if (include_pressure_gradient_in_ustar_)
        {
          Cerr << "Methode incremental pour le grad(P)" << finl;
          Cerr << " Option codee uniquement pour le sch_euler... Tester et implementer si besoiN. " << finl;
          Process::exit();
        }
      // OPTION A SELECTIONNE DANS LE CAS DUN SHEAR PERIO POUR EVITER LES PBMS D INTERPOLATION DE RHO AU NIVEAU DU BORD PERIO_Z
      if (use_inv_rho_in_poisson_solver_)
        {
          // Cerr << "Methode basee sur inv rho pour le grad(P) en RK3" << finl;
          // Cerr << " Option a tester si besoin. " << finl;

          pressure_projection_with_inv_rho(inv_rho_field_, velocity_[0], velocity_[1],  velocity_[2], pressure_,
                                           fractionnal_timestep,
                                           pressure_rhs_, check_divergence_, poisson_solver_);

        }
      else
        {
#ifdef PROJECTION_DE_LINCREMENT_DV
          // On l'a fait avant pour etre sur qu'elle soit bien dans la derivee stockee...
#else


          pressure_projection_with_rho(rho_field_, velocity_[0], velocity_[1],  velocity_[2], pressure_,
                                       fractionnal_timestep,
                                       pressure_rhs_, check_divergence_, poisson_solver_);

          // GAB TODO : cest a peu pres ici qu'il faudra travailler pour recuperer le
          // terme de pression
#endif
          // GAB TODO : checker si le passage de rho_n a rho_n+1 est bon
          // chercher ca pour l'etape de deplacement de rho : maj_indicatrice_rho_mu
        }

      // GAB, qdm : on recupere ici le terme grad(p),
      terme_pression_bis_ = calculer_grad_p_moyen(pressure_);
      // GAB, qdm : on recupere ici le terme de pression (1/rho * grad(p))
      terme_pression_ter_ = calculer_grad_p_over_rho_moyen(pressure_);
      pression_ap_proj_ += calculer_v_moyen(pressure_);

      //statistiques().end_count(projection_counter_);
    }

  if (Process::je_suis_maitre())
    {
      Cout << "Timings diff=" << statistiques().last_time(diffusion_counter_)
           << " conv=" << statistiques().last_time(convection_counter_);
      Cout << " src=" << statistiques().last_time(source_counter_)
           << finl;
    }
  statistiques().end_count(euler_rk3_counter_);

}

// Hard coded constant pressure gradient in i direction, add contribution in m/s*volume of control volume
void IJK_FT_double::terme_source_gravite(IJK_Field_double& dv, int k_index, int dir) const
{
  const double constant = gravite_[dir];
  const int imax = dv.ni();
  const int jmax = dv.nj();
  for (int j = 0; j < jmax; j++)
    for (int i = 0; i < imax; i++)
      {
        dv(i, j, k_index) += constant;
      }
}

void IJK_FT_double::euler_explicit_update(const IJK_Field_double& dv, IJK_Field_double& v,
                                          const int k_layer) const
{
  const double delta_t = timestep_;
  const int imax = v.ni();
  const int jmax = v.nj();
  for (int j = 0; j < jmax; j++)
    for (int i = 0; i < imax; i++)
      {
        double x = dv(i,j,k_layer);
        v(i,j,k_layer) += x * delta_t;
      }
}

// Deplacement des interfaces par le champ de vitesse :
// 1. Calcule vitesse_ft (etendue) a partir du champ de vitesse.
// 2. Supprime les duplicatas.
// 3. Transporte le maillage avec velocity_ft_.
// 4. Transfere les bulles a travers la periodicite si besoin.
// 5. Cree les nouveaux duplicatas.
//
// ATTENTION : rho_mu_indicatrice ne sont pas mis a jours.
//
// Mettre rk_step = -1 si schema temps different de rk3.
void IJK_FT_double::deplacer_interfaces(const double timestep, const int rk_step,
                                        ArrOfDouble& var_volume_par_bulle,
                                        const int first_step_interface_smoothing)
{
  static Stat_Counter_Id deplacement_interf_counter_ = statistiques().new_counter(1, "Deplacement de l'interface");
  statistiques().begin_count(deplacement_interf_counter_);

  if (disable_diphasique_ || interfaces_.is_frozen())
    return;
  //  Calculer vitesse_ft (etendue) a partir du champ de vitesse.
  {

    // Il faut ajouter le cisaillement moyen au champ de vitesse_ft avant de convecter les marqueurs

    redistribute_to_splitting_ft_faces_[2].redistribute(velocity_[2], velocity_ft_[2]);
    redistribute_to_splitting_ft_faces_[1].redistribute(velocity_[1], velocity_ft_[1]);
    redistribute_to_splitting_ft_faces_[0].redistribute(velocity_[0], velocity_ft_[0]);

    if (IJK_Splitting::defilement_ == 1)
      {
        // after redistribute, velocity in ft domain must be shifted by the shear
        velocity_ft_[0].redistribute_with_shear_domain_ft(velocity_[0], boundary_conditions_.get_dU_perio(boundary_conditions_.get_resolution_u_prime_()), ijk_splitting_ft_extension_);
        velocity_ft_[1].redistribute_with_shear_domain_ft(velocity_[1], 0., ijk_splitting_ft_extension_);
        velocity_ft_[2].redistribute_with_shear_domain_ft(velocity_[2], 0., ijk_splitting_ft_extension_ );
      }

    for (int dir = 0; dir < 3; dir++)
      velocity_ft_[dir].echange_espace_virtuel(velocity_ft_[dir].ghost());

  }

  /*
   * Calculation of intersections on interface at time (n)
   */
  Cerr << "Compute Eulerian distance and curvature fields" << finl;
  thermals_.compute_eulerian_distance_curvature();
  Cerr << "Clean IJK intersections" << finl;
  thermals_.clean_ijk_intersections();
  // thermals_.update_intersections(); // no need as IJK_intersections call interfaces_nI interfaces_xI
  Cerr << "Compute compo_connex from bounding box" << finl;
  interfaces_.compute_compo_connex_from_bounding_box();
  Cerr << "Compute compo_connex from interface compo in mixed cells" << finl;
  interfaces_.compute_compo_connex_from_interface();
  Cerr << "Compute rising velocity from compo connex (barycentre calc)" << finl;
  interfaces_.compute_rising_velocities_from_compo();

  // On supprime les duplicatas avant le transport :
  interfaces_.supprimer_duplicata_bulles();

  /*
   * TODO: Advect with zero velocity at the beggining of the simulation
   * to use the remeshing algo
   */
  // if (counter_first_iter_ && first_step_interface_smoothing_)
  interfaces_.transporter_maillage(timestep/* total meme si RK3*/,
                                   var_volume_par_bulle,
                                   rk_step, current_time_,
                                   first_step_interface_smoothing);

  // Apres le transport, est-ce que certaines bulles reeles sont trop proche du bord
  // du domaine etendu? Si on en trouve, on les transferts :
  interfaces_.transferer_bulle_perio();

  // On supprime les fragments de bulles.
  //interfaces_.detecter_et_supprimer_rejeton(false);

  // On cree les duplicatas pour le prochain pas de temps avant le posttraitement.
  // Comme ca, ils seront visibles a la visu.
  interfaces_.creer_duplicata_bulles();

  // On met a jour l'indicatrice du pas de temps d'apres.
  // On met aussi a jour le surf et bary des faces mouillees,
  // les valeurs moyennes en ijk, les val moy en ijkf, etc

  update_indicator_field();

  if (counter_first_iter_ && first_step_interface_smoothing_)
    thermals_.recompute_temperature_init();
  // mise a jour de l'indicatrice pour les variables monofluides
  const double delta_rho = rho_liquide_ - rho_vapeur_;
  interfaces_.switch_indicatrice_next_old();
  interfaces_.calculer_indicatrice_next(
    post_.potentiel(),
    gravite_,
    delta_rho,
    sigma_,
    /*Pour post-traitement : post_.rebuilt_indic()
    */
#ifdef SMOOTHING_RHO
    /* Pour le smoothing : */
    rho_field_ft_,
    rho_vapeur_,
    smooth_density_,
#endif
    current_time_, tstep_
  );

  if (IJK_Splitting::defilement_ == 1)
    {
      // mise a jour de l'indicatrice pour les variables monofluides
      redistribute_from_splitting_ft_elem_ghostz_min_.redistribute(interfaces_.I_ft(), I_ns_);
      redistribute_from_splitting_ft_elem_ghostz_max_.redistribute(interfaces_.I_ft(), I_ns_);
      rho_field_.set_indicatrice_ghost_zmin_(I_ns_, 0);
      rho_field_.set_indicatrice_ghost_zmax_(I_ns_, rho_field_.nk()-4);
      if (use_inv_rho_)
        {
          inv_rho_field_.set_indicatrice_ghost_zmin_(I_ns_, 0);
          inv_rho_field_.set_indicatrice_ghost_zmax_(I_ns_, inv_rho_field_.nk()-4);
        }
    }

  statistiques().end_count(deplacement_interf_counter_);
}

// Nouvelle version ou le transport se fait avec les ghost...
void IJK_FT_double::deplacer_interfaces_rk3(const double timestep, const int rk_step,
                                            ArrOfDouble& var_volume_par_bulle)
{
  if (disable_diphasique_ || interfaces_.is_frozen())
    return;
  //  Calculer vitesse_ft (etendue) a partir du champ de vitesse.
  static Stat_Counter_Id deplacement_interf_counter_ = statistiques().new_counter(1, "Deplacement de l'interface");
  statistiques().begin_count(deplacement_interf_counter_);

  for (int dir = 0; dir < 3; dir++)
    redistribute_to_splitting_ft_faces_[dir].redistribute(velocity_[dir], velocity_ft_[dir]);

  redistribute_to_splitting_ft_faces_[2].redistribute(velocity_[2], velocity_ft_[2]);
  redistribute_to_splitting_ft_faces_[1].redistribute(velocity_[1], velocity_ft_[1]);
  redistribute_to_splitting_ft_faces_[0].redistribute(velocity_[0], velocity_ft_[0]);

  if (IJK_Splitting::defilement_ == 1)
    {
      // after redistribute, velocity in ft domain must be shifted by the shear
      velocity_ft_[0].redistribute_with_shear_domain_ft(velocity_[0], boundary_conditions_.get_dU_perio(boundary_conditions_.get_resolution_u_prime_()), ijk_splitting_ft_extension_);
      velocity_ft_[1].redistribute_with_shear_domain_ft(velocity_[1], 0., ijk_splitting_ft_extension_);
      velocity_ft_[2].redistribute_with_shear_domain_ft(velocity_[2], 0., ijk_splitting_ft_extension_);
    }

  for (int dir = 0; dir < 3; dir++)
    velocity_ft_[dir].echange_espace_virtuel(velocity_ft_[dir].ghost());


  // On conserve les duplicatas que l'on transporte comme le reste.

  // Normalement, transporter_maillage gere aussi les duplicatas...
  interfaces_.transporter_maillage(timestep/* total meme si RK3*/,
                                   var_volume_par_bulle,
                                   rk_step, current_time_);

  statistiques().end_count(deplacement_interf_counter_);
  // On verra a la fin du pas de temps si certaines bulles reeles sont trop proche du bord
  // du domaine etendu. Pour l'instant, dans les sous dt, on ne les transferts pas.

  // On a conserve les duplicatas donc pas besoin de les re-creer...
  // On calcule l'indicatrice du prochain pas de temps (qui correspond aux interfaces qu'on
  // vient de deplacer.

  update_indicator_field();
  const double delta_rho = rho_liquide_ - rho_vapeur_;
  interfaces_.switch_indicatrice_next_old();
  interfaces_.calculer_indicatrice_next(
    post_.potentiel(),
    gravite_,
    delta_rho,
    sigma_,
    /*Pour post-traitement : post_.rebuilt_indic()
    */
#ifdef SMOOTHING_RHO
    /* Pour le smoothing : */
    rho_field_ft_,
    rho_vapeur_,
    smooth_density_,
#endif
    current_time_, rk_step
  );
  if (IJK_Splitting::defilement_ == 1)
    {
      redistribute_from_splitting_ft_elem_ghostz_min_.redistribute(interfaces_.I_ft(), I_ns_);
      redistribute_from_splitting_ft_elem_ghostz_max_.redistribute(interfaces_.I_ft(), I_ns_);
      rho_field_.set_indicatrice_ghost_zmin_(I_ns_, 0);
      rho_field_.set_indicatrice_ghost_zmax_(I_ns_, rho_field_.nk()-4);
      if (use_inv_rho_)
        {
          inv_rho_field_.set_indicatrice_ghost_zmin_(I_ns_, 0);
          inv_rho_field_.set_indicatrice_ghost_zmax_(I_ns_, inv_rho_field_.nk()-4);
        }
    }

}

//  Parcourir_maillage cree des noeuds et facettes virtuelles.
//  Pour maintenir un tableau a jour entre avant et apres,
//  il suffit de resizer le tableau a la sortie de la methode
//  (forcement plus grand qu'avant) et de faire un echange_espace_virtuel.
void IJK_FT_double::parcourir_maillage()
{
  //const int nbsom_before = interfaces_.maillage_ft_ijk().nb_sommets();
  interfaces_.parcourir_maillage();

  const int nbsom = interfaces_.maillage_ft_ijk().nb_sommets();
  // const int size_store = interfaces_.RK3_G_store_vi().dimension(0);
  // if (!((nbsom >= nbsom_before) &&
  //       ((nbsom_before == size_store) || (0 == size_store) )))
  //   {
  //     Cerr << "Une des tailles de tableau n'est pas bonne... "
  //          << " size_store = " << size_store
  //          << " nbsom_before = " << nbsom_before
  //          << " nbsom = " << nbsom
  //          << finl;
  //     Process::exit();
  //   }
  interfaces_.RK3_G_store_vi_resize(nbsom, 3);
  interfaces_.RK3_G_store_vi_echange_esp_vect();
}

// Maj indicatrice rho mu met indicatrice a indicatrice next
// et maj rho et mu en fonction de la nouvelle indicatrice
void IJK_FT_double::maj_indicatrice_rho_mu(const bool parcourir)
{
  // En monophasique, les champs sont a jours donc on zap :
  if (disable_diphasique_)
    return;

  static Stat_Counter_Id calculer_rho_mu_indicatrice_counter_= statistiques().new_counter(2, "Calcul Rho Mu Indicatrice");
  statistiques().begin_count(calculer_rho_mu_indicatrice_counter_);

  // En diphasique sans bulle (pour cas tests)
  if (interfaces_.get_nb_bulles_reelles() == 0)
    {
      rho_field_.data() = rho_liquide_;
      rho_field_.echange_espace_virtuel(rho_field_.ghost());
      if (use_inv_rho_)
        {
          inv_rho_field_.data() = 1./rho_liquide_;
          inv_rho_field_.echange_espace_virtuel(inv_rho_field_.ghost());
        }
      molecular_mu_.data() = mu_liquide_;
      molecular_mu_.echange_espace_virtuel(molecular_mu_.ghost());

      if (parcourir)
        interfaces_.parcourir_maillage();
      return;
    }

  if (parcourir)
    parcourir_maillage();

  // Nombre de mailles du domaine NS :
  const int nx = interfaces_.I().ni();
  const int ny = interfaces_.I().nj();
  const int nz = interfaces_.I().nk();

  if (use_inv_rho_)
    {
      for (int k=0; k < nz; k++)
        for (int j=0; j < ny; j++)
          for (int i=0; i < nx; i++)
            {
              double chi_l = interfaces_.I(i,j,k);
              rho_field_(i,j,k)    = rho_liquide_ * chi_l + (1.- chi_l) * rho_vapeur_;
              inv_rho_field_(i,j,k) = 1./rho_liquide_ * chi_l + (1.- chi_l) * 1./rho_vapeur_;
              molecular_mu_(i,j,k) = mu_liquide_ * chi_l + (1.- chi_l) * mu_vapeur_ ;
            }
    }
  else
    {
      for (int k=0; k < nz; k++)
        for (int j=0; j < ny; j++)
          for (int i=0; i < nx; i++)
            {
              double chi_l = interfaces_.I(i,j,k);
              rho_field_(i,j,k)    = rho_liquide_ * chi_l + (1.- chi_l) * rho_vapeur_;
              molecular_mu_(i,j,k) = mu_liquide_ * chi_l + (1.- chi_l) * mu_vapeur_ ;
            }
    }

  //Mise a jour des espaces virtuels des champs :
  rho_field_.echange_espace_virtuel(rho_field_.ghost());
  if (use_inv_rho_)
    inv_rho_field_.echange_espace_virtuel(inv_rho_field_.ghost());
  molecular_mu_.echange_espace_virtuel(molecular_mu_.ghost());

#ifdef SMOOTHING_RHO
  if (smooth_density_)
    {
      smoothing_field(rho_field_, rho_liquide_, rho_vapeur_,
                      ratio_density_max_);

      redistribute_to_splitting_ft_elem_.redistribute(rho_field_, rho_field_ft_);
    }
#endif
  statistiques().end_count(calculer_rho_mu_indicatrice_counter_);
}

void IJK_FT_double::update_rho_v()
{
  if (use_inv_rho_for_mass_solver_and_calculer_rho_v_)
    {
      // rho_face = 2*(rho_gauche*rho_droite)/(rho_gauche+rho_droite)
      //          = 1./ (1/2 * (1/rho_g + 1/rho_d))
      // 1/rho_face est donc la moyenne geometrique de inv_rho.
      calculer_rho_harmonic_v(rho_field_, velocity_, rho_v_);
    }
  else
    {
      calculer_rho_v(rho_field_, velocity_, rho_v_);
    }
  rho_v_[0].echange_espace_virtuel(rho_v_[0].ghost());
  rho_v_[1].echange_espace_virtuel(rho_v_[1].ghost());
  rho_v_[2].echange_espace_virtuel(rho_v_[2].ghost());
}

void IJK_FT_double::update_v_ghost_from_rho_v()
{
  for (int dir = 0 ; dir < 3 ; dir++)
    {
      const int imax = velocity_[dir].ni();
      const int jmax = velocity_[dir].nj();
      const int kmax = velocity_[dir].nk();
      const int ghost = velocity_[dir].ghost();
      int last_global_k = splitting_.get_nb_items_global(IJK_Splitting::ELEM, 2);
      for (int j = 0; j < jmax; j++)
        {
          for (int i = 0; i < imax; i++)
            {
              if(splitting_.get_offset_local(2)==0)
                {
                  for (int k = -ghost; k < 0; k++)
                    {
                      double rho = 0.;
                      double DU = 0.;
                      if (dir==0)
                        {
                          rho = 0.5*(rho_field_(i, j, k) + rho_field_(i-1, j, k));
                          DU = boundary_conditions_.get_dU_perio();
                        }
                      else if (dir==1)
                        {
                          rho = 0.5*(rho_field_(i, j, k) + rho_field_(i, j-1, k));
                        }
                      else if (dir==2)
                        {
                          rho= 0.5*(rho_field_(i, j, k) + rho_field_(i, j, k-1));
                        }

                      velocity_[dir](i, j, k) = rho_v_[dir](i, j, k)/rho - DU;
                    }
                }
              if(splitting_.get_offset_local(2)+kmax==last_global_k)
                {
                  for (int k = kmax; k < kmax + ghost; k++)
                    {
                      double rho = 0.;
                      double DU = 0.;
                      if (dir==0)
                        {
                          rho = 0.5*(rho_field_(i, j, k) + rho_field_(i-1, j, k));
                          DU = boundary_conditions_.get_dU_perio();
                        }
                      else if (dir==1)
                        {
                          rho = 0.5*(rho_field_(i, j, k) + rho_field_(i, j-1, k));
                        }
                      else if (dir==2)
                        {
                          rho= 0.5*(rho_field_(i, j, k) + rho_field_(i, j, k-1));
                        }

                      velocity_[dir](i, j, k) = rho_v_[dir](i, j, k)/rho + DU;
                    }
                }
            }
        }
    }
}


// Transfert du maillage ft vers ns de champs aux faces :
void IJK_FT_double::transfer_ft_to_ns()
{
  for (int dir = 0; dir < 3; dir++)
    {
      redistribute_from_splitting_ft_faces_[dir].redistribute(terme_repulsion_interfaces_ft_[dir],
                                                              terme_repulsion_interfaces_ns_[dir]);
      redistribute_from_splitting_ft_faces_[dir].redistribute(terme_abs_repulsion_interfaces_ft_[dir],
                                                              terme_abs_repulsion_interfaces_ns_[dir]);
    }
}

void IJK_FT_double::fill_variable_source_and_potential_phi(const double time)
{
  // Ajout du terme source (force acceleration)
  // Solveur masse pour chaque composante du bilan de QdM
  const FixedVector<IJK_Field_double, 3>& grad_I_ns = post_.get_grad_I_ns();
  for (int dir = 0; dir < 3; dir++)
    {
      // Si on est en presence d'une source analytique variable spatialement:
      if (expression_variable_source_[dir] != "??")
        set_field_data(variable_source_[dir], expression_variable_source_[dir], interfaces_.I(), grad_I_ns[dir], time);
      else if (expression_potential_phi_ != "??")
        {
          // Pour Remettre a zero la source:
          set_field_data(variable_source_[dir], Nom("0."), grad_I_ns[dir], time);
        }
// modification de la valeur du temps (localement) pour s'avancer a la fin du sous-pas de temps :
//          const double end_of_substep_time = time + compute_fractionnal_timestep_rk3(timestep_,rk_step);
//          set_field_data(variable_source_[dir], expression_variable_source_[dir], grad_I_ns_[dir], end_of_substep_time);
    }

  if (expression_potential_phi_ != "??")
    {
      set_field_data(potential_phi_, expression_potential_phi_, time);
      potential_phi_.echange_espace_virtuel(potential_phi_.ghost());
      add_gradient_times_constant(potential_phi_, 1.,variable_source_[0],variable_source_[1],variable_source_[2]);
    }
}

// GAB, rotation
int IJK_FT_double::get_direction(const ArrOfDouble& vecteur)
{
  if (vecteur[0]==0. && vecteur[1]==0. && vecteur[2]!=0)
    return DIRECTION_K;
  else if (vecteur[0]==0. && vecteur[1]!=0. && vecteur[2]==0.)
    return DIRECTION_J;
  else
    return DIRECTION_I;
}

// GAB, qdm : construction du terme de pression pour le bilan de qdm. Je trouve ea plus logique de bouger cette fonction dans
//            IJK_Navier_Stokes_tool, mais etant donne qu'elle se trouve dans IJK_kernel, je ne sais pas trop si c'est propre que j'y touche
//            voir avec Guillaume ce qu'il en dit.
//            inspire de : void IJK_FT_Post::calculer_gradient_indicatrice_et_pression(const IJK_Field_double& indic)
Vecteur3 IJK_FT_double::calculer_inv_rho_grad_p_moyen(const IJK_Field_double& rho,const IJK_Field_double& pression)
{
  FixedVector<IJK_Field_double, 3> champ;
  allocate_velocity(champ, splitting_, 1);
  Vecteur3 resu;

  // Remise a zero :
  for (int dir = 0; dir < 3; dir++)
    {
      champ[dir].data() = 0.;
      // je ne sais pas a quoi ca sert, mais il est appele avant add_gradient_times_constant_over_rho dans IJK_NS_tool
      champ[dir].echange_espace_virtuel(1);
    }
  add_gradient_times_constant_over_rho(pression, rho, 1. /*constant*/,  champ[0], champ[1], champ[2]);
  for (int dir = 0; dir < 3; dir++)
    {
      // Je ne sais pas e quoi cette ligne sert. elle est dans calculer_gradient_indicatrice_et_pression, alors je copie
      champ[dir].echange_espace_virtuel(1);
    }
  // calculer_rho_v(inv_rho, champ, inv_rho_champ);
  for (int dir=0; dir<3; dir++)
    resu[dir] = calculer_v_moyen(champ[dir]);

  return resu;
}


Vecteur3 IJK_FT_double::calculer_grad_p_moyen(const IJK_Field_double& pression)
{
  FixedVector<IJK_Field_double, 3> champ;
  allocate_velocity(champ, splitting_, 1);
  Vecteur3 resu;

  // Remise a zero :
  for (int dir = 0; dir < 3; dir++)
    {
      champ[dir].data() = 0.;
      // je ne sais pas a quoi ca sert, mais il est appele avant add_gradient_times_constant_over_rho dans IJK_NS_tool
      champ[dir].echange_espace_virtuel(1);
    }
  add_gradient_times_constant(pression, 1. /*constant*/,  champ[0], champ[1], champ[2]);
  for (int dir = 0; dir < 3; dir++)
    {
      // Je ne sais pas e quoi cette ligne sert. elle est dans calculer_gradient_indicatrice_et_pression, alors je copie
      champ[dir].echange_espace_virtuel(1);
    }
  for (int dir=0; dir<3; dir++)
    resu[dir] = calculer_v_moyen(champ[dir]);

  return resu;
}

Vecteur3 IJK_FT_double::calculer_grad_p_over_rho_moyen(const IJK_Field_double& pression)
{
  /*
   * Calcule Moyenne_spatiale{ 1/rho * grad(p) }
   * */
  FixedVector<IJK_Field_double, 3> champ;
  allocate_velocity(champ, splitting_, 1);
  Vecteur3 resu;

  // Remise a zero :
  for (int dir = 0; dir < 3; dir++)
    {
      champ[dir].data() = 0.;
      // je ne sais pas a quoi ca sert, mais il est appele avant add_gradient_times_constant_over_rho dans IJK_NS_tool
      champ[dir].echange_espace_virtuel(1);
    }
  add_gradient_times_constant_over_rho(pression, rho_field_, 1. /*constant*/,  champ[0], champ[1], champ[2]);
  for (int dir = 0; dir < 3; dir++)
    {
      // Je ne sais pas e quoi cette ligne sert. elle est dans calculer_gradient_indicatrice_et_pression, alors je copie
      champ[dir].echange_espace_virtuel(1);
    }
  for (int dir=0; dir<3; dir++)
    resu[dir] = calculer_v_moyen(champ[dir]);

  return resu;
}

void IJK_FT_double::write_check_etapes_et_termes(int rk_step)
{
  // GR : 07.01.22 : mettre une condition if (tstep % dt_post_stats_check_ == dt_post_stats_check_ - 1 || stop)
  //                 serai vraiment une bonne chose pour alleger les _check_....out... voir avec AB et GB.
  if (Process::je_suis_maitre())
    {
      if (test_etapes_et_bilan_)
        {
          int reset_test = (!reprise_) && (tstep_==0);
          // GAB, qdm RK3 : accurate_current_time_ = t0(1+ 1/4); t0(1+ 1/4 + 5/12); t0(1+ 1/4 + 5/12 + 1/3)
          double accurate_current_time = 0.0;
          if (get_time_scheme()==EULER_EXPLICITE)
            accurate_current_time = current_time_;
          else if (get_time_scheme()==RK3_FT)
            accurate_current_time = current_time_at_rk3_step_;

          SFichier fic_test=Ouvrir_fichier("_check_etapes_et_termes.out",
                                           "tstep\tcurrent_time\trk_step"
                                           "\trho_u_euler_av_prediction\trho_du_euler_ap_prediction"
                                           "\trho_u_euler_ap_projection\trho_du_euler_ap_projection"
                                           "\trho_u_euler_av_rho_mu_ind\trho_u_euler_ap_rho_mu_ind"
                                           "\tu_euler_ap_rho_mu_ind"
                                           "\ttemre_interfaces\tterme_convection\tterme_diffusion"
                                           "\tterme_pression_bis\tterme_pression_ter"
                                           "\tterme_interfaces_bf_mass_solver_bis\tterme_interfaces_bf_mass_solver\tterme_interfaces_bf_mass_solver"
                                           "\tpression_ap_proj"
                                           "\tdrho_u"
                                           "\tterme_moyen_convection_mass_solver"
                                           "\tterme_moyen_diffusion_mass_solver"
                                           "\n# Forces have 3 components.",
                                           reset_test, 20/*prec*/);


          // GAB, qdm
          Cout << "check etapes et termes" << finl;
          fic_test << tstep_ << " " << accurate_current_time << " ";
          fic_test << rk_step << " ";
          // Inspection a gros grain
          fic_test << rho_u_euler_av_prediction_[0] << " "
                   << rho_u_euler_av_prediction_[1] << " "
                   << rho_u_euler_av_prediction_[2] << " ";  // colone 5
          fic_test << rho_du_euler_ap_prediction_[0] << " "
                   << rho_du_euler_ap_prediction_[1] << " "
                   << rho_du_euler_ap_prediction_[2] << " ";
          fic_test << rho_u_euler_ap_projection_[0] << " "
                   << rho_u_euler_ap_projection_[1] << " "  // colonne 10
                   << rho_u_euler_ap_projection_[2] << " ";
          fic_test << rho_du_euler_ap_projection_[0] << " "
                   << rho_du_euler_ap_projection_[1] << " "
                   << rho_du_euler_ap_projection_[2] << " ";
          fic_test << rho_u_euler_av_rho_mu_ind_[0] << " "  // colonne 15
                   << rho_u_euler_av_rho_mu_ind_[1] << " "
                   << rho_u_euler_av_rho_mu_ind_[2] << " ";
          fic_test << rho_u_euler_ap_rho_mu_ind_[0] << " "
                   << rho_u_euler_ap_rho_mu_ind_[1] << " "
                   << rho_u_euler_ap_rho_mu_ind_[2] << " ";  // colonne 20
          fic_test << u_euler_ap_rho_mu_ind_[0] << " "
                   << u_euler_ap_rho_mu_ind_[1] << " "
                   << u_euler_ap_rho_mu_ind_[2] << " ";          // Dans l'etape de prediction, inspection terme a terme
          fic_test << terme_interfaces_[0] << " "
                   << terme_interfaces_[1] << " "  // colonne 25
                   << terme_interfaces_[2] << " ";
          fic_test << terme_convection_[0] << " "
                   << terme_convection_[1] << " "
                   << terme_convection_[2] << " ";
          fic_test << terme_diffusion_[0] << " "  // colonne 30
                   << terme_diffusion_[1] << " "
                   << terme_diffusion_[2] << " ";
          // Moyenne_spatiale{ grad(p) }
          fic_test << terme_pression_bis_[0] << " "
                   << terme_pression_bis_[1] << " "
                   << terme_pression_bis_[2] << " ";  // colonne 35
          // Moyenne_spatiale{ 1/rho grad(p) }
          fic_test << terme_pression_ter_[0] << " "
                   << terme_pression_ter_[1] << " "
                   << terme_pression_ter_[2] << " ";
          // Inspection de l'effet du mass solver
          fic_test << terme_interfaces_bf_mass_solver_bis_[0] << " "
                   << terme_interfaces_bf_mass_solver_bis_[1] << " "  // colonne 40
                   << terme_interfaces_bf_mass_solver_bis_[2] << " ";
          fic_test << terme_interfaces_bf_mass_solver_[0] << " "
                   << terme_interfaces_bf_mass_solver_[1] << " "
                   << terme_interfaces_bf_mass_solver_[2] << " ";
          fic_test << terme_interfaces_af_mass_solver_[0] << " "  // colonne 45 /!\ au 17.01.22 une coquille a ete corrigee : terme_interfaces_bf_mass_solver-> terme_interfaces_af_mass_solver
                   << terme_interfaces_af_mass_solver_[1] << " "
                   << terme_interfaces_af_mass_solver_[2] << " ";
          fic_test << pression_ap_proj_ << " ";
          // On pourrai se passer de cette sortie et la creer uniquement dans python
          // ==> ( r^{n+1} - r^{n} ) * u^{n+1}
          fic_test << rho_u_euler_av_rho_mu_ind_[0]-rho_u_euler_ap_rho_mu_ind_[0] << " "
                   << rho_u_euler_av_rho_mu_ind_[1]-rho_u_euler_ap_rho_mu_ind_[1] << " "  // colonne 50
                   << rho_u_euler_av_rho_mu_ind_[2]-rho_u_euler_ap_rho_mu_ind_[2] << " ";
          fic_test << terme_moyen_convection_mass_solver_[0] << " "
                   << terme_moyen_convection_mass_solver_[1] << " "
                   << terme_moyen_convection_mass_solver_[2] << " ";
          fic_test << terme_moyen_diffusion_mass_solver_[0] << " "  // colonne 55
                   << terme_moyen_diffusion_mass_solver_[1] << " "
                   << terme_moyen_diffusion_mass_solver_[2] << " ";
          fic_test << finl;
          fic_test.close();
          Cout << "check etapes et termes fini" << finl;
        }
    }
}

// -----------------------------------------------------------------------------------
//  CORRECTION DE QdM
double IJK_FT_double::calculer_true_moyenne_de_phase_liq(const IJK_Field_double& vx)
{
  /* Au 04.11.21 : Renvoi vx_liq */
  double alpha_liq_vx_liq = calculer_moyenne_de_phase_liq(vx);
  double alpha_liq = calculer_v_moyen(interfaces_.I()); //en utilisant calculer_moyenne_de_phase_liq(indicatrice_ns_) on somme des chi^2 ce qui est genant aux mailles diphasiques
  return alpha_liq_vx_liq/alpha_liq;
}

double IJK_FT_double::calculer_moyenne_de_phase_liq(const IJK_Field_double& vx)
{
  /* Au 04.11.21 : Renvoi alpha_liq * vx_liq
   * */
  const IJK_Splitting& splitting = vx.get_splitting();
  const IJK_Grid_Geometry& geom = splitting.get_grid_geometry();
  const int ni = vx.ni();
  const int nj = vx.nj();
  const int nk = vx.nk();
  double v_moy = 0.;

#ifndef VARIABLE_DZ
  for (int k = 0; k < nk; k++)
    for (int j = 0; j < nj; j++)
      for (int i = 0; i < ni; i++)
        v_moy += vx(i,j,k)*interfaces_.I(i,j,k);

  // somme sur tous les processeurs.
  v_moy = Process::mp_sum(v_moy);
  // Maillage uniforme, il suffit donc de diviser par le nombre total de mailles:
  // cast en double au cas ou on voudrait faire un maillage >2 milliards
  const double n_mailles_tot = ((double) geom.get_nb_elem_tot(0)) * geom.get_nb_elem_tot(1) * geom.get_nb_elem_tot(2);
  v_moy /= n_mailles_tot;
#else
  const int offset = splitting.get_offset_local(DIRECTION_K);
  const ArrOfDouble& tab_dz=geom.get_delta(DIRECTION_K);
  for (int k = 0; k < nk; k++)
    {
      const double dz = tab_dz[k+offset];
      for (int j = 0; j < nj; j++)
        for (int i = 0; i < ni; i++)
          v_moy += vx(i,j,k)*dz;
    }
  // somme sur tous les processeurs.
  v_moy = Process::mp_sum(v_moy);
  // Maillage uniforme, il suffit donc de diviser par le nombre total de mailles:
  // cast en double au cas ou on voudrait faire un maillage >2 milliards
  const double n_mailles_xy = ((double) geom.get_nb_elem_tot(0)) * geom.get_nb_elem_tot(1);
  v_moy /= (n_mailles_xy * geom.get_domain_length(DIRECTION_K) );
#endif
  return v_moy;
}

double IJK_FT_double::calculer_true_moyenne_de_phase_vap(const IJK_Field_double& vx)
{
  /* Au 04.11.21 : Renvoi vx_vap */
  double alpha_vap_vx_vap = calculer_moyenne_de_phase_vap(vx);
  double alpha_vap = 1-calculer_v_moyen(interfaces_.I()); //en utilisant 1-calculer_moyenne_de_phase_liq(indicatrice_ns_) on somme des chi^2 ce qui est genant aux mailles diphasiques
  return alpha_vap_vx_vap/alpha_vap;
}

double IJK_FT_double::calculer_moyenne_de_phase_vap(const IJK_Field_double& vx)
{
  /* Au 04.11.21 : Renvoi alpha_vap * vx_vap
   * */
  const IJK_Splitting& splitting = vx.get_splitting();
  const IJK_Grid_Geometry& geom = splitting.get_grid_geometry();
  const int ni = vx.ni();
  const int nj = vx.nj();
  const int nk = vx.nk();
  double v_moy = 0.;

#ifndef VARIABLE_DZ
  for (int k = 0; k < nk; k++)
    for (int j = 0; j < nj; j++)
      for (int i = 0; i < ni; i++)
        v_moy += vx(i,j,k)*(1-interfaces_.I(i,j,k));

  // somme sur tous les processeurs.
  v_moy = Process::mp_sum(v_moy);
  // Maillage uniforme, il suffit donc de diviser par le nombre total de mailles:
  // cast en double au cas ou on voudrait faire un maillage >2 milliards
  const double n_mailles_tot = ((double) geom.get_nb_elem_tot(0)) * geom.get_nb_elem_tot(1) * geom.get_nb_elem_tot(2);
  v_moy /= n_mailles_tot;
#else
  const int offset = splitting.get_offset_local(DIRECTION_K);
  const ArrOfDouble& tab_dz=geom.get_delta(DIRECTION_K);
  for (int k = 0; k < nk; k++)
    {
      const double dz = tab_dz[k+offset];
      for (int j = 0; j < nj; j++)
        for (int i = 0; i < ni; i++)
          v_moy += vx(i,j,k)*dz;
    }
  // somme sur tous les processeurs.
  v_moy = Process::mp_sum(v_moy);
  // Maillage uniforme, il suffit donc de diviser par le nombre total de mailles:
  // cast en double au cas ou on voudrait faire un maillage >2 milliards
  const double n_mailles_xy = ((double) geom.get_nb_elem_tot(0)) * geom.get_nb_elem_tot(1);
  v_moy /= (n_mailles_xy * geom.get_domain_length(DIRECTION_K) );
#endif
  return v_moy;
}

void IJK_FT_double::set_time_for_corrections()
{
  qdm_corrections_.set_time(current_time_,timestep_,tstep_);
}

void IJK_FT_double::compute_and_add_qdm_corrections()
{
  /*
   * Corrections to comply with the momentum budget
   * u_corrected = u - u_correction
   * u_corrected and u_corrrection are computed in the qdm_corrections_ object.
   * */
  double alpha_l = calculer_v_moyen(interfaces_.I());
  double rho_moyen = calculer_v_moyen(rho_field_);
  qdm_corrections_.set_rho_moyen_alpha_l(rho_moyen,alpha_l);
  qdm_corrections_.set_rho_liquide(rho_liquide_);
  for (int dir=0; dir<3; dir++)
    {
      IJK_Field_double& vel = velocity_[dir];
      double rho_vel_moyen = calculer_v_moyen(rho_v_[dir]);
      qdm_corrections_.set_rho_vel_moyen(dir,rho_vel_moyen);
      Cout << "qdm_corrections_.get_need_for_vitesse_relative("<<dir<<")" << qdm_corrections_.get_need_for_vitesse_relative(dir)<< finl;
      if (qdm_corrections_.get_need_for_vitesse_relative(dir))
        {
          double vel_rel = calculer_true_moyenne_de_phase_liq(vel) - calculer_true_moyenne_de_phase_vap(vel);
          qdm_corrections_.set_vitesse_relative(dir, vel_rel);
        }
      // TODO : Demander de l'aide a Guillaume :
      /* Pour moyenne glissante, je fais appel a des ArrOfDouble. En utilisation sequentielle, j'en
       * suis satisfait disons. En utilisation parallele, je ne sais pas vraiment comment sont gerees
       * mes listes. Sont-t-elles dupliquees ? Decoupees ? En tout cas la simu plante pour une liste
       * plus de 10 doubles, sur un maillage a 40^3 mailles, une seule bulle... */
      if (qdm_corrections_.get_need_to_compute_correction_value_one_direction(dir))
        qdm_corrections_.compute_correction_value_one_direction(dir);
      for (int k=0; k<vel.nk(); ++k)
        for (int j=0; j<vel.nj(); ++j)
          for (int i=0; i<vel.ni(); ++i)
            {
              qdm_corrections_.compute_correct_velocity_one_direction(dir, vel(i,j,k));
              velocity_[dir](i,j,k) = qdm_corrections_.get_correct_velocitiy_one_direction(dir);
            }
    }

  Cout << "AF : compute_and_add_qdm_corrections" << finl;
}

void IJK_FT_double::compute_and_add_qdm_corrections_monophasic()
{
  /* For monophasic corrections only */
  /*
   * Corrections to comply with the momentum budget
   * u_corrected = u - u_correction
   * u_corrected and u_corrrection are computed in the qdm_corrections_ object.
   * */
  double alpha_l = 1.; // calculer_v_moyen(interfaces_.I())
  double rho_moyen = rho_liquide_; // calculer_v_moyen(rho_field_)
  qdm_corrections_.set_rho_moyen_alpha_l(rho_moyen,alpha_l);
  qdm_corrections_.set_rho_liquide(rho_liquide_);
  for (int dir=0; dir<3; dir++)
    {
      IJK_Field_double& vel = velocity_[dir];
      double rho_vel_moyen = calculer_v_moyen(rho_v_[dir]);
      qdm_corrections_.set_rho_vel_moyen(dir,rho_vel_moyen);
      Cout << "qdm_corrections_.get_need_for_vitesse_relative("<<dir<<")" << qdm_corrections_.get_need_for_vitesse_relative(dir)<< finl;
      if (qdm_corrections_.get_need_for_vitesse_relative(dir))
        {
          double vel_rel = calculer_true_moyenne_de_phase_liq(vel); // - calculer_true_moyenne_de_phase_vap(vel);
          qdm_corrections_.set_vitesse_relative(dir, vel_rel);
        }
      // TODO : Demander de l'aide a Guillaume :
      /* Pour moyenne glissante, je fais appel a des ArrOfDouble. En utilisation sequentielle, j'en
       * suis satisfait disons. En utilisation parallele, je ne sais pas vraiment comment sont gerees
       * mes listes. Sont-t-elles dupliquees ? Decoupees ? En tout cas la simu plante pour une liste
       * plus de 10 doubles, sur un maillage a 40^3 mailles, une seule bulle... */
      if (qdm_corrections_.get_need_to_compute_correction_value_one_direction(dir))
        qdm_corrections_.compute_correction_value_one_direction(dir);
      for (int k=0; k<vel.nk(); k++)
        for (int j=0; j<vel.nj(); j++)
          for (int i=0; i<vel.ni(); i++)
            {
              qdm_corrections_.compute_correct_velocity_one_direction(dir, vel(i,j,k));
              velocity_[dir](i,j,k) = qdm_corrections_.get_correct_velocitiy_one_direction(dir);
            }
    }
  Cout << "AF : compute_and_add_qdm_corrections_monophasic" << finl;
}

void IJK_FT_double::compute_var_volume_par_bulle(ArrOfDouble& var_volume_par_bulle)
{
  if (vol_bulles_.size_array() > 0.)
    {
      ArrOfDouble volume_reel;
      DoubleTab position;
      interfaces_.calculer_volume_bulles(volume_reel, position);
      const int nb_reelles = interfaces_.get_nb_bulles_reelles();
      for (int ib = 0; ib < nb_reelles; ib++)
        var_volume_par_bulle[ib] = volume_reel[ib] - vol_bulles_[ib];
      // Pour les ghost : on retrouve leur vrai numero pour savoir quel est leur volume...
      for (int i = 0; i < interfaces_.get_nb_bulles_ghost(0 /* no print*/); i++)
        {
          const int ighost = interfaces_.ghost_compo_converter(i);
          const int ibulle_reelle = decoder_numero_bulle(-ighost);
          // Cerr << " aaaa " << i << " " << ighost << " " << ibulle_reelle << finl;
          var_volume_par_bulle[nb_reelles + i] = volume_reel[nb_reelles+ i] - vol_bulles_[ibulle_reelle];
        }
    }
}

void IJK_FT_double::redistribute_to_splitting_ft_elem(const IJK_Field_double& input_field,
                                                      IJK_Field_double& output_field)
{
  redistribute_to_splitting_ft_elem_.redistribute(input_field, output_field);
}

void IJK_FT_double::redistribute_from_splitting_ft_elem(const IJK_Field_double& input_field,
                                                        IJK_Field_double& output_field)
{
  redistribute_from_splitting_ft_elem_.redistribute(input_field, output_field);
}

void IJK_FT_double::copy_field_values(IJK_Field_double& field, const IJK_Field_double& field_to_copy)
{
  const int ni = field.ni();
  const int nj = field.nj();
  const int nk = field.nk();
  const int ni_to_copy = field_to_copy.ni();
  const int nj_to_copy = field_to_copy.nj();
  const int nk_to_copy = field_to_copy.nk();
  int bool_dim = (ni == ni_to_copy && nj == nj_to_copy && nk == nk_to_copy);
  assert(bool_dim);
  if (!bool_dim)
    Process::exit();
  for (int k = 0; k < nk; k++)
    for (int j = 0; j < nj; j++)
      for (int i = 0; i < ni; i++)
        field(i,j,k) = field_to_copy(i,j,k);
  // TODO: M.G. Test with the following operator
  // field.data() = field_to_copy.data();
  field.echange_espace_virtuel(field.ghost());
}

void IJK_FT_double::update_indicator_field()
{
  const double delta_rho = rho_liquide_ - rho_vapeur_;
  interfaces_.switch_indicatrice_next_old();
  interfaces_.calculer_indicatrice_next(post_.potentiel(),
                                        gravite_,
                                        delta_rho,
                                        sigma_,
                                        /*Pour post-traitement : post_.rebuilt_indic()
                                        */
#ifdef SMOOTHING_RHO
                                        /* Pour le smoothing : */
                                        rho_field_ft_,
                                        rho_vapeur_,
                                        smooth_density_,
#endif
                                        current_time_, tstep_
                                       );
}

void IJK_FT_double::update_twice_indicator_field()
{
  for(int i=0; i<2; i++)
    update_indicator_field();
}

void IJK_FT_double::write_qdm_corrections_information()
{
  // Impression dans le fichier qdm_correction.out
  if (get_time_scheme() == RK3_FT) // && (rk3_sub_step!=0)
    Cout << "in write_qdm_corrections_information, rk_step_ = "<<rk_step_<<finl;
  Vecteur3 qdm_cible = qdm_corrections_.get_correction_values();
  Vecteur3 velocity_correction = qdm_corrections_.get_velocity_corrections();
  Vecteur3 rho_vel;
  for (int dir=0; dir<3; ++dir) rho_vel[dir] =  calculer_v_moyen(scalar_fields_product(rho_field_,velocity_[dir],dir));

  if (Process::je_suis_maitre())
    {
      int reset = (!reprise_) && (tstep_==0);
      SFichier fic=Ouvrir_fichier("_qdm_correction.out",
                                  "1.iteration\t2.time\t3.qdm_cible[0]\t4.qdm_cible[1]\t5.qdm_cible[2]\t6.velocity_correction[0]\t7.velocity_correction[1]\t8.velocity_correction[2]\t9.qdm[0]\t10.qdm[1]\t11.qdm[2]",
                                  reset);
      // temps
      fic << tstep_ << " ";
      fic << current_time_ << " ";
      // CIBLE CONSTANE : qdm_cible = al.rl.u_cible
      fic << qdm_cible[0] << " ";
      fic << qdm_cible[1] << " ";
      fic << qdm_cible[2] << " ";
      // velocity_correction = (<r.u> - qdm_cible) / <r>
      fic << velocity_correction[0] << " ";
      fic << velocity_correction[1] << " ";
      fic << velocity_correction[2] << " ";
      // <r.u>
      fic << rho_vel[0] << " ";
      fic << rho_vel[1] << " ";
      fic << rho_vel[2] << " ";
      fic<<finl;
      fic.close();
      //	 << finl;
    }
}


// -----------------------------------------------------------------------------------
//  PRODUITS DE CHAMPS
IJK_Field_double IJK_FT_double::scalar_product(const FixedVector<IJK_Field_double, 3>& V1, const FixedVector<IJK_Field_double, 3>& V2)
{
  /*
   * * ATTENTION : valide pour un maillage cartesien, de maille cubiques uniquement !
   */
  IJK_Field_double resu;
  resu.allocate(splitting_, IJK_Splitting::ELEM, 3);
  int nk = V1[0].nk();
  if (nk != V2[0].nk()) {Cerr << "scalar product of fields with different dimensions (nk)"<< finl;}
  int nj = V1[0].nj();
  if (nj != V2[0].nj()) {Cerr << "scalar product of fields with different dimensions (nj)"<< finl;}
  int ni = V1[0].ni();
  if (ni != V2[0].ni()) {Cerr << "scalar product of fields with different dimensions (ni)"<< finl;}

  for (int k=0; k<nk; ++k)
    for (int j=0; j<nj; ++j)
      for (int i=0; i<ni; ++i)
        {
          resu(i,j,k) = 0.25*(
                          (V1[0](i,j,k)+V1[0](i+1,j,k))*(V2[0](i,j,k)+V2[0](i+1,j,k))
                          +(V1[1](i,j,k)+V1[1](i,j+1,k))*(V2[1](i,j,k)+V2[1](i,j+1,k))
                          +(V1[2](i,j,k)+V1[2](i,j,k+1))*(V2[2](i,j,k)+V2[2](i,j,k+1))
                        );
        }
  // Communication avec tous les process ?
  return resu;
}

FixedVector<IJK_Field_double, 3> IJK_FT_double::scalar_times_vector(const IJK_Field_double& Sca, const FixedVector<IJK_Field_double, 3>& Vec)
{
  /*
   * Produit d'un champ scalaire (Sca) par un champ de vecteur (Vec).
   * Le champ scalaire est aux centre des elements, le champ de vecteur est aux faces
   * Le resultat reste localise au meme endroit que le champ de vecteur passe en entree.
   * ATTENTION : valide pour un maillage cartesien, de maille cubiques uniquement !
   */

  FixedVector<IJK_Field_double, 3> resu;
  allocate_velocity(resu,splitting_,3); // j'ai besoin de mettre des cellules ghost ? non, je ne pense pas
  int nk = Vec[0].nk();
  if (nk != Sca.nk()) {Cerr << "scalar fields has different dimension from vector field  (nk)"<< finl;}
  int nj = Vec[0].nj();
  if (nj != Sca.nj()) {Cerr << "scalar fields has different dimension from vector field  (nj)"<< finl;}
  int ni = Vec[0].ni();
  if (ni != Sca.nk()) {Cerr << "scalar fields has different dimension from vector field  (ni)"<< finl;}

  for (int k=0; k<nk; ++k)
    for (int j=0; j<nj; ++j)
      for (int i=0; i<ni; ++i)
        {
          resu[0](i,j,k) = 0.5*(Sca(i-1,j,k)+Sca(i,j,k))*Vec[0](i,j,k);
          resu[1](i,j,k) = 0.5*(Sca(i,j-1,k)+Sca(i,j,k))*Vec[1](i,j,k);
          resu[2](i,j,k) = 0.5*(Sca(i,j,k-1)+Sca(i,j,k))*Vec[2](i,j,k);
        }
  // Communication avec tous les process ?
  return resu;
}

IJK_Field_double IJK_FT_double::scalar_fields_product(const IJK_Field_double& S1, const IJK_Field_double& S2, int dir)
{
  /*
   * Produit d'un champ scalaire aux centres (S1) par une des composantes d'un champ de vecteur (S2).
   * Le resultat est localise au meme endroit que le champ de vecteur dont est issu S2.
   * ATTENTION : valide pour un maillage cartesien, de maille cubiques uniquement !
   */
  IJK_Field_double resu;
  resu.allocate(splitting_, IJK_Splitting::ELEM, 3);
  int nk = S1.nk();
  if (nk != S2.nk()) {Cerr << "scalar fields have different dimensions for the product (nk)"<< finl;}
  int nj = S1.nj();
  if (nj != S2.nj()) {Cerr << "scalar fields have different dimensions for the product (nj)"<< finl;}
  int ni = S1.ni();
  if (ni != S2.ni()) {Cerr << "scalar fields have different dimensions for the product (ni)"<< finl;}

  for (int k=0; k<nk; ++k)
    for (int j=0; j<nj; ++j)
      for (int i=0; i<ni; ++i)
        {
          if (dir==0) {resu(i,j,k) = 0.5*(S1(i-1,j,k)+S1(i,j,k))*S2(i,j,k);}
          if (dir==1) {resu(i,j,k) = 0.5*(S1(i,j-1,k)+S1(i,j,k))*S2(i,j,k);}
          if (dir==2) {resu(i,j,k) = 0.5*(S1(i,j,k-1)+S1(i,j,k))*S2(i,j,k);}
        }
  // Communication avec tous les process ?
  return resu;
}
