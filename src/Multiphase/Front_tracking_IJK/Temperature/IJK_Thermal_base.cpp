/****************************************************************************
* Copyright (c) 2023, CEA
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
// File      : IJK_Thermal_base.cpp
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/Temperature
//
/////////////////////////////////////////////////////////////////////////////

#include <IJK_Thermal_base.h>
#include <Param.h>
#include <stat_counters.h>
#include <DebogIJK.h>
#include <Corrige_flux_FT.h>
#include <IJK_FT.h>
#include <IJK_FT_Post.h>
#include <IJK_switch_FT.h>
#include <IJK_Ghost_Fluid_tools.h>
#include <IJK_Bubble_tools.h>

Implemente_base_sans_constructeur( IJK_Thermal_base, "IJK_Thermal_base", Objet_U ) ;

/********************************************
 * Methods inherited from Objet_U
 ********************************************/

IJK_Thermal_base::IJK_Thermal_base()
{
  rang_ = 0; // default value, used as an index for the list of thermal sub-problems
  calulate_grad_T_=0;
  calculate_local_energy_=0;

  /*
   * Physical parameters
   */
  dt_fo_=1.e20;
  cp_liquid_=1.;
  cp_vapour_=0.;
  lambda_liquid_=1.;
  lambda_vapour_=0.;
  single_phase_=1.;
  fo_ = 1.; // Fourier number

  /*
   * Initialisation of the problem
   */
  expression_T_init_="??";
  fichier_reprise_temperature_= "??";
  timestep_reprise_temperature_=1;

  conv_temperature_negligible_=0;
  diff_temperature_negligible_=0;

  type_T_source_="??";
  wall_flux_=0;
  lambda_variable_=0; // terme source variable
  // type_thermal_problem_=0;

  /*
   * Ceci est une initialisation des derivees des temperatures moyenne de chaque phase
   * Il n'est peut-etre pas pertinent de les mettre ici
   */
  dTv_ = 0.;
  dTl_ = 1.;
  Tl_ = 0.;
  Tv_ = 1.;
  Tv0_ = 1.;  // Serait-ce plutot Tref (une temperature de reference pour reconstruire le champ dim??)
  kl_ = -100000000000000.;
  kv_ = -200000000000000.;
  T0v_ = 1.;
  T0l_ = 0.;
  expression_T_ana_="??";
  upstream_temperature_=-1.1e20;

  E0_ = 0;
  uniform_lambda_ = 0;
  uniform_alpha_ = 0;
  prandtl_number_=0;
  global_energy_ = 0;
  conserv_energy_global_ = 0;
  vol_ = 0.;
  min_delta_xyz_ = 0.;
  cell_diagonal_ = 0.;
  ghost_cells_ = 2;
  rho_cp_post_ = 0;

  nb_diam_upstream_ = 0;
  side_temperature_ = 0;
  stencil_side_ = 2;

  n_iter_distance_ = 6;

  gfm_recompute_field_ini_ = 1;
  gfm_zero_neighbour_value_mean_ = 0;
  gfm_vapour_mixed_only_ = 1;

  compute_grad_T_interface_ = 0;
  compute_curvature_ = 0;
  compute_distance_= 0;
  ghost_fluid_ = 0;
  compute_grad_T_elem_ = 0;
  compute_hess_T_elem_ = 0;
  compute_hess_diag_T_elem_ = 0;
  compute_hess_cross_T_elem_ = 0;

  compute_eulerian_compo_ = 0;
  compute_rising_velocities_ = 0;
  fill_rising_velocities_ = 0;

  debug_ = 0;
  spherical_approx_ = 1;
  spherical_exact_ = 0;

  eulerian_compo_connex_ft_ = nullptr;
  eulerian_compo_connex_ns_ = nullptr;
  eulerian_compo_connex_ghost_ft_ = nullptr;
  eulerian_compo_connex_ghost_ns_ = nullptr;
  eulerian_compo_connex_from_interface_ns_ = nullptr;
  eulerian_compo_connex_from_interface_ft_ = nullptr;
  eulerian_compo_connex_from_interface_ghost_ft_ = nullptr;
  eulerian_compo_connex_from_interface_ghost_ns_ = nullptr;
  eulerian_compo_connex_from_interface_int_ns_ = nullptr;
  eulerian_compo_connex_from_interface_ghost_int_ns_= nullptr;

  liquid_velocity_ = 0.;
  latastep_reprise_=0;
  latastep_reprise_ini_=0;
}

Sortie& IJK_Thermal_base::printOn( Sortie& os ) const
{
  Nom front_space = "    ";
  Nom end_space = " ";
  Nom escape = "\n";
  Objet_U::printOn( os );
  os << "  {" << escape;

  os << escape;
  os << front_space << "# BASE PARAMS #" << escape;
  os << escape;

  os << "    boundary_conditions {"  << escape;
  /*
   * Boundary conditions (Periodicity or wall)
   */
  Nom bctype_kmin, bctype_kmax, bckmin, bckmax, valeur_kmin, valeur_kmax;
  if( boundary_conditions_.get_bctype_k_max()==boundary_conditions_.Paroi_Temperature_imposee)
    {
      bctype_kmax="Paroi_Temperature_imposee";
      bckmax = "temperature_imposee_kmax";
      valeur_kmax = boundary_conditions_.get_temperature_kmax();
    }
  else if( boundary_conditions_.get_bctype_k_max()==boundary_conditions_.Paroi_Flux_impose)
    {
      bctype_kmax="Paroi_Flux_impose";
      bckmax = "flux_impose_kmax";
      valeur_kmax = boundary_conditions_.get_flux_kmax();
    }
  else if( boundary_conditions_.get_bctype_k_max()==boundary_conditions_.Perio)
    {
      bctype_kmax="Perio";
      bctype_kmin="Perio";
      bckmax = " ";
      bckmin = " ";
      valeur_kmax = " ";
    }
  if( boundary_conditions_.get_bctype_k_min()==boundary_conditions_.Paroi_Temperature_imposee)
    {
      bctype_kmin="Paroi_Temperature_imposee";
      bckmin = "temperature_imposee_kmin";
      valeur_kmin = boundary_conditions_.get_temperature_kmin();
    }
  else if( boundary_conditions_.get_bctype_k_min()==boundary_conditions_.Paroi_Flux_impose)
    {
      bctype_kmin="Paroi_Flux_impose";
      bckmin = "flux_impose_kmin";
      valeur_kmin = boundary_conditions_.get_flux_kmin();
    }
  else if( boundary_conditions_.get_bctype_k_min()==boundary_conditions_.Perio)
    {
      bctype_kmin="Perio";
      bctype_kmin="Perio";
      bckmin = " ";
      bckmin = " ";
      valeur_kmin = " ";
    }
  os<< "      bctype_kmin" << " " << bctype_kmin << " \n";
  os<< "      bctype_kmax" << " " << bctype_kmax << " \n";
  os<< "      " << bckmin << " " << valeur_kmin << " \n";
  os<< "      " << bckmax << " " << valeur_kmax << " \n";
  os<< "    } \n" ;

  /*
   * Physical parameters
   */
  os << front_space << "lambda_liquid" << end_space << lambda_liquid_ << escape;
  os << front_space << "lambda_vapour" << end_space << lambda_vapour_ << escape;
  os << front_space << "cp_liquid" << end_space << cp_liquid_ << escape;
  os << front_space << "cp_vapour" << end_space << cp_vapour_ << escape;

  /*
   * Source term
   */
  if (type_T_source_!="??")
    os<< "    type_T_source " << type_T_source_ << "\n";
  if (type_T_source_=="SWARM")
    {
      os<< "      kl_source " <<  kl_ << "\n";
      os<< "      kv_source " <<  kv_ << "\n";
      os<< "      T0l_source " <<  T0l_ << "\n";
      os<< "      T0v_source " <<  T0v_ << "\n";
    }

  if( wall_flux_)
    os << front_space << "wall_flux" << escape;

  /*
   * Resume calculation
   */
  os << front_space << "fichier_reprise_temperature" << end_space << basename(fichier_reprise_temperature_)  << escape;
  os << front_space << "timestep_reprise_temperature" << end_space << timestep_reprise_temperature_ << escape;
  os << front_space << "latastep_reprise" << end_space << latastep_reprise_ << escape;

  /*
   * Analytical expression of temperature at t_initial
   */
  if ( expression_T_ana_!="??")
    os << front_space << "expression_T_ana" <<  end_space << expression_T_ana_ << escape;


  os << front_space << "upstream_temperature" << end_space << upstream_temperature_ << escape;
  os << front_space << "nb_diam_upstream" << end_space << nb_diam_upstream_ << escape;
  os << front_space << "side_temperature" << end_space << side_temperature_ << escape;
  os << front_space << "stencil_side" << end_space << stencil_side_ << escape;
  os << front_space << "n_iter_distance" << end_space << n_iter_distance_ << escape;


  os << front_space << "temperature_diffusion_op" << end_space << temperature_diffusion_op_ << escape;
  os << front_space << "temperature_convection_op" << end_space << temperature_convection_op_ << escape;

  /*
   * Neglect an operator
   */

  os << escape;
  os << front_space << "# BASE FLAGS #" << escape;
  os << escape;

  if ( conv_temperature_negligible_)
    os << front_space << "conv_temperature_negligible" << escape;
  if ( diff_temperature_negligible_)
    os << front_space << "diff_temp_negligible" << escape;
  if (ghost_fluid_)
    os << front_space << "ghost_fluid" <<  escape;
  if (spherical_exact_)
    os << front_space << "spherical_exact" <<  escape;
  if (debug_)
    os << front_space << "debug" <<  escape;
  if (calculate_local_energy_)
    os << front_space << "calculate_local_energy" <<  escape;

  return os;
}

// XD thermique listobj thermique -1 thermique_bloc 1 to add energy equation resolution if needed
// XD thermique_bloc interprete nul 1 not_set
Entree& IJK_Thermal_base::readOn( Entree& is )
{
  /*
   * Parse the datafile
   */
  Param param(que_suis_je());
  set_param(param);
  param.lire_avec_accolades(is);
  Cout << "IJK_Thermal_base::readOn : Parameters summary. " << finl;
  printOn(Cout);
  return is;
}

void IJK_Thermal_base::set_fichier_reprise(const char *lataname)
{
  fichier_reprise_temperature_ = lataname;
}

void IJK_Thermal_base::set_param(Param& param)
{
  param.ajouter("fo", &fo_); // XD_ADD_P floattant not_set
  param.ajouter("cp_liquid", &cp_liquid_, Param::REQUIRED); // XD_ADD_P floattant Liquid specific heat at constant pressure
  param.ajouter("lambda_liquid", &lambda_liquid_, Param::REQUIRED); // XD_ADD_P floattant Liquid thermal conductivity
  param.ajouter("cp_vapour", &cp_vapour_); // XD_ADD_P floattant Liquid specific heat at constant pressure
  param.ajouter("lambda_vapour", &lambda_vapour_); // XD_ADD_P floattant Liquid thermal conductivity
  param.ajouter("expression_T_init", &expression_T_init_); // XD_ADD_P chaine Expression of initial temperature (parser of x,y,z)
  param.ajouter("boundary_conditions", &boundary_conditions_, Param::REQUIRED); // XD_ADD_P bloc_lecture boundary conditions
  param.ajouter("type_T_source", &type_T_source_); // XD_ADD_P chaine(into=["dabiri","patch_dabiri","unweighted_dabiri"]) source term
  param.ajouter("expression_source_temperature", &expression_source_temperature_); // XD_ADD_P chaine source terms
  param.ajouter_flag("lambda_variable", &lambda_variable_);
  param.ajouter_flag("wall_flux", &wall_flux_); // XD_ADD_P rien not_set
  param.ajouter("kl_source", &kl_);
  param.ajouter("kv_source", &kv_);
  param.ajouter("T0l_source", &T0l_);
  param.ajouter("T0v_source", &T0v_);
  param.ajouter("fichier_reprise_temperature", &fichier_reprise_temperature_);
  param.ajouter("timestep_reprise_temperature", &timestep_reprise_temperature_);
  param.ajouter("latastep_reprise", &latastep_reprise_ini_);
  param.ajouter_flag("conv_temperature_negligible", &conv_temperature_negligible_); // XD_ADD_P rien neglect temperature convection
  param.ajouter_flag("diff_temperature_negligible", &diff_temperature_negligible_); // XD_ADD_P rien neglect temperature diffusion
  param.ajouter("temperature_diffusion_op", &temperature_diffusion_op_);
  param.ajouter("temperature_convection_op", &temperature_convection_op_);
  param.ajouter("expression_T_ana", &expression_T_ana_); // XD_ADD_P chaine Analytical expression T=f(x,y,z,t) for post-processing only
  param.ajouter("calculate_local_energy", &calculate_local_energy_);
  param.ajouter("upstream_temperature", &upstream_temperature_);
  param.ajouter("nb_diam_upstream", &nb_diam_upstream_);
  param.ajouter("side_temperature", &side_temperature_);
  param.ajouter("stencil_side", &stencil_side_);
  param.ajouter("n_iter_distance", &n_iter_distance_);
  param.ajouter_flag("ghost_fluid", &ghost_fluid_);
  param.ajouter_flag("spherical_exact", &spherical_exact_);
  param.ajouter_flag("debug", &debug_);
  //  param.ajouter_flag("gfm_recompute_field_ini", &gfm_recompute_field_ini_);
  //  param.ajouter_flag("gfm_zero_neighbour_value_mean", &gfm_zero_neighbour_value_mean_);
  //  param.ajouter_flag("gfm_vapour_mixed_only", &gfm_vapour_mixed_only_);
}

/********************************************
 * Public methods
 ********************************************/

int IJK_Thermal_base::initialize_switch(const IJK_Splitting& splitting, const int idx)
{
  int nalloc = 0;
  temperature_.allocate(splitting, IJK_Splitting::ELEM, 1);
  nalloc += 1;
  if (fichier_reprise_temperature_ == "??") // si on ne fait pas une reprise on initialise V
    {
      Cerr << "Please provide initial conditions for temperature, either by an expression or a field for restart."
           << "You should consider using either fichier_reprise_temperature or expression_T_init keywords. " << finl;
      Process::exit();
    }
  else
    {
      Cout << "Reading initial temperature field T" << rang_ <<" from file "
           << fichier_reprise_temperature_ << " timestep= " << timestep_reprise_temperature_ << finl;
      const Nom& geom_name = splitting.get_grid_geometry().le_nom();
      lire_dans_lata(fichier_reprise_temperature_, timestep_reprise_temperature_, geom_name, Nom("TEMPERATURE_") + Nom(idx),
                     temperature_); // fonction qui lit un champ a partir d'un lata .
      temperature_.echange_espace_virtuel(temperature_.ghost()); // It is essential to fill the EV because the first call to convection needs them.
    }
  return nalloc;
}

int IJK_Thermal_base::initialize(const IJK_Splitting& splitting, const int idx)
{
  //  Cout << que_suis_je() << "::initialize()" << finl;
  rang_ = idx;
  int nalloc = 0;

  latastep_reprise_ = latastep_reprise_ini_;

  /*
   * Diffusion operator:
   * If temperature_diffusion_op_ is not written in the .data
   * the operator is initialised with uniform_lambda & centre2
   * in Operateur_IJK_elem_diff
   */
  if (single_phase_)
    temperature_diffusion_op_.typer_diffusion_op("uniform");
  else
    temperature_diffusion_op_.typer_diffusion_op("standard");


  /*
   * Convection operator
   * If temperature_convection_op_ is not written in the .data
   * the operator is initialised with quick
   * in Operateur_IJK_elem_conv.h
   */

  /*
   * Initialise the operators
   */
  if (!conv_temperature_negligible_)
    temperature_convection_op_.initialize(splitting);
  if (!diff_temperature_negligible_)
    temperature_diffusion_op_.initialize(splitting);

  /*
   * Corrige Flux FT
   */
  corrige_flux_.typer("Corrige_flux_FT_temperature_conv");

  /*
   * Fields
   */
  temperature_.allocate(splitting, IJK_Splitting::ELEM, ghost_cells_);
  temperature_for_ini_per_bubble_.allocate(splitting, IJK_Splitting::ELEM, 1);
  d_temperature_.allocate(splitting, IJK_Splitting::ELEM, 2);
  nalloc += 3;
  compute_cell_volume();
  compute_min_cell_delta();

  //if (!diff_temperature_negligible_)
  {
    div_coeff_grad_T_volume_.allocate(splitting, IJK_Splitting::ELEM, 0);
    nalloc += 1;
    div_coeff_grad_T_volume_.data() = 0.;
  }
  // if (!conv_temperature_negligible_)
  {
    u_T_convective_volume_.allocate(splitting, IJK_Splitting::ELEM, 0);
    nalloc += 1;
    u_T_convective_volume_.data() = 0.;
  }

  rho_cp_post_ = (liste_post_instantanes_.size() && liste_post_instantanes_.contient_("RHO_CP"));
  if (rho_cp_post_)
    {
      rho_cp_.allocate(splitting, IJK_Splitting::ELEM, 2);
      nalloc += 1;
    }
  if (calculate_local_energy_)
    {
      rho_cp_T_.allocate(splitting, IJK_Splitting::ELEM, 2);
      nalloc += 1;
    }

  /*
   * Storage for temperature gradient post-processing or method
   */
  if ((ref_ijk_ft_.non_nul()) && (!ref_ijk_ft_->get_disable_diphasique()))
    {
      Cout << "Allocating fields temperature_ft_ and storage" << finl;
      allocate_cell_vector(storage_, ref_ijk_ft_->get_splitting_ft(), 1);
      nalloc += 3;
    }

  /*
   * Dimensionless temperature field (thermostat)
   */
  if ((wall_flux_) || liste_post_instantanes_.contient_("SOURCE_TEMPERATURE")
      || liste_post_instantanes_.contient_("TEMPERATURE_PHYSIQUE_T")
      || liste_post_instantanes_.contient_("TEMPERATURE_ADIMENSIONNELLE_THETA")
      || (type_T_source_ != "??"))
    {
      Cout << "Allocating field for the thermal source term & co. " << finl;
      source_temperature_.allocate(splitting, IJK_Splitting::ELEM, 1);
      source_temperature_v_.allocate(splitting, IJK_Splitting::ELEM, 1);
      source_temperature_l_.allocate(splitting, IJK_Splitting::ELEM, 1);
      d_source_Tv_.allocate(splitting, IJK_Splitting::ELEM, 1);
      d_source_Tl_.allocate(splitting, IJK_Splitting::ELEM, 1);
      temperature_physique_T_.allocate(splitting, IJK_Splitting::ELEM, 2);
      temperature_adimensionnelle_theta_.allocate(splitting, IJK_Splitting::ELEM, 2);
      nalloc += 7;
      // par defaut s'il n'y a pas de source renseignee, on utilise la source de Dabiri/Kawamura
      // cela veut dire que dans le cas des SWARMS il faut imperativement renseigner le nom de
      // la source
      if (type_T_source_ == "??")
        {
          Cerr << "Attention on demande des post-traitement sans avoir renseigner type_T_source" << finl;
          throw "Erreur post et type_T_source";
        }
    }
  if (liste_post_instantanes_.contient_("TEMPERATURE_ADIM_BULLES"))
    {
      temperature_adim_bulles_.allocate(splitting, IJK_Splitting::ELEM, 2);
      nalloc += 1;
    }

  /*
   * RK3 sub-steps
   * Check that pointer is not null:
   */
  if (ref_ijk_ft_.non_nul() && ref_ijk_ft_->get_time_scheme()== ref_ijk_ft_->RK3_FT)
    {
      RK3_F_temperature_.allocate(splitting, IJK_Splitting::ELEM, 0);
      nalloc +=1;
    }

  if (liste_post_instantanes_.size() && (liste_post_instantanes_.contient_("TEMPERATURE_ANA")
                                         || liste_post_instantanes_.contient_("ECART_T_ANA") || liste_post_instantanes_.contient_("ECART_T_ANA_REL")))
    {
      temperature_ana_.allocate(splitting, IJK_Splitting::ELEM, 1);
      ecart_t_ana_.allocate(splitting, IJK_Splitting::ELEM, 1);
      ecart_t_ana_rel_.allocate(splitting, IJK_Splitting::ELEM, 1);
      nalloc +=3;
      temperature_ana_.data() = 0.;
      ecart_t_ana_.data() = 0.;
      ecart_t_ana_rel_.data() = 0.;
      temperature_ana_.echange_espace_virtuel(temperature_ana_.ghost());
      ecart_t_ana_.echange_espace_virtuel(ecart_t_ana_.ghost());
      ecart_t_ana_rel_.echange_espace_virtuel(ecart_t_ana_rel_.ghost());
    }

  /*
   * Resume the calculation
   */
  if (fichier_reprise_temperature_ == "??")   // si on ne fait pas une reprise on initialise V
    {
      if (expression_T_init_ != "??")
        {
          compute_temperature_init();
        }
      else
        {
          Cerr << "Please provide initial conditions for temperature, either by an expression or a field for restart. "
               << "You should consider using either fichier_reprise_temperature or expression_T_init keywords. " << finl;
          Process::exit();
        }
    }
  else
    {
      Cout << "Reading initial temperature field T"<< rang_<<" from file " << fichier_reprise_temperature_ << " timestep= " << timestep_reprise_temperature_ << finl;
      const Nom& geom_name = splitting.get_grid_geometry().le_nom();
      lire_dans_lata(fichier_reprise_temperature_, timestep_reprise_temperature_, geom_name, Nom("TEMPERATURE_")+Nom(idx),
                     temperature_); // fonction qui lit un champ a partir d'un lata .
      temperature_.echange_espace_virtuel(temperature_.ghost()); // It is essential to fill the EV because the first call to convection needs them.
    }

  /*
   * List of post-processed data
   */
  Cerr << " Initializing thermal fields dependant on the post-pro list : " ;
  if (liste_post_instantanes_.size())
    Cerr << liste_post_instantanes_;
  else
    Cerr << "empty";
  Cerr << finl;

  if ((liste_post_instantanes_.contient_("SOURCE_TEMPERATURE_ANA")) || (liste_post_instantanes_.contient_("ECART_SOURCE_TEMPERATURE_ANA")) )
    {
      source_temperature_ana_.allocate(splitting, IJK_Splitting::ELEM, 1);
      ecart_source_t_ana_.allocate(splitting, IJK_Splitting::ELEM, 1);
      nalloc += 2;
    }

  // TODO: Check with Aymeric
  calulate_grad_T_ = (liste_post_instantanes_.size() && liste_post_instantanes_.contient_("GRAD_T"))
                     || (ref_ijk_ft_.non_nul() && ref_ijk_ft_->t_debut_statistiques() <  1.e10 );
  if (calulate_grad_T_)
    {
      allocate_velocity(grad_T_, splitting, 1);
      nalloc += 3;
    }

  compute_grad_T_interface_ = ghost_fluid_ || compute_grad_T_interface_ || (liste_post_instantanes_.size() && liste_post_instantanes_.contient_("GRAD_T_INTERFACE"));
  compute_curvature_ = compute_curvature_ || compute_grad_T_interface_ || (liste_post_instantanes_.size() && liste_post_instantanes_.contient_("CURVATURE"));
  compute_distance_ = compute_distance_ || compute_curvature_ || (liste_post_instantanes_.size() && liste_post_instantanes_.contient_("DISTANCE"));

  if (compute_distance_)
    {
      /*
       * TODO: Move to IJK_Interfaces
       */
      // Laplacian(d) necessitates 2 ghost cells like temperature
      eulerian_distance_ft_.allocate(ref_ijk_ft_->get_splitting_ft(), IJK_Splitting::ELEM, 2);
      nalloc += 1;
      // grad(d) necessitates 1 ghost cell ?
      allocate_cell_vector(eulerian_normal_vectors_ft_, ref_ijk_ft_->get_splitting_ft(), 1);
      // allocate_velocity(eulerian_normal_vectors_, ref_ijk_ft_->get_splitting_ft(), 1);
      nalloc += 3;
      allocate_cell_vector(eulerian_facets_barycentre_ft_, ref_ijk_ft_->get_splitting_ft(), 0);
      nalloc += 3;
      eulerian_distance_ft_.echange_espace_virtuel(eulerian_distance_ft_.ghost());
      eulerian_normal_vectors_ft_.echange_espace_virtuel();
      eulerian_facets_barycentre_ft_.echange_espace_virtuel();
      /*
       * TODO: This is already calculated in IJK_Interfaces
       * Keep it for now and clean later
       */
      eulerian_distance_ns_.allocate(splitting, IJK_Splitting::ELEM, 2);
      allocate_cell_vector(eulerian_normal_vectors_ns_, splitting, 1);
      allocate_cell_vector(eulerian_facets_barycentre_ns_, splitting, 0);
      nalloc += 7;
      eulerian_distance_ns_.echange_espace_virtuel(eulerian_distance_ns_.ghost());
      eulerian_normal_vectors_ns_.echange_espace_virtuel();
      eulerian_facets_barycentre_ns_.echange_espace_virtuel();
      allocate_cell_vector(eulerian_normal_vectors_ns_normed_, splitting, 1);
      nalloc += 3;
      eulerian_normal_vectors_ns_normed_.echange_espace_virtuel();
    }
  if (compute_curvature_)
    {
      // Laplacian(d) necessitates 0 ghost cells like div_lambda_grad_T
      // but if calculated using the neighbours maybe 1
      eulerian_curvature_ft_.allocate(ref_ijk_ft_->get_splitting_ft(), IJK_Splitting::ELEM, 1);
      nalloc += 1;
      eulerian_curvature_ft_.echange_espace_virtuel(eulerian_curvature_ft_.ghost());
      // Only calculated in the mixed cells ghost_cells = 0
      eulerian_interfacial_area_ft_.allocate(ref_ijk_ft_->get_splitting_ft(), IJK_Splitting::ELEM, 0);
      nalloc += 1;
      eulerian_interfacial_area_ft_.echange_espace_virtuel(eulerian_interfacial_area_ft_.ghost());
      /*
       * TODO: This is already calculated in IJK_Interfaces
       * Keep it for now and clean later
       */
      eulerian_curvature_ns_.allocate(splitting, IJK_Splitting::ELEM, 1);
      eulerian_curvature_ns_.echange_espace_virtuel(eulerian_curvature_ns_.ghost());
      nalloc += 2;
      eulerian_interfacial_area_ns_.allocate(splitting, IJK_Splitting::ELEM, 0);
      eulerian_interfacial_area_ns_.echange_espace_virtuel(eulerian_interfacial_area_ns_.ghost());
    }
  if (compute_grad_T_interface_)
    {
      // 1 ghost cell for eulerian_grad_T_interface_ and temperature_ft_ to access its neighbour
      eulerian_grad_T_interface_ft_.allocate(ref_ijk_ft_->get_splitting_ft(), IJK_Splitting::ELEM, 1);
      nalloc += 1;
      eulerian_grad_T_interface_ft_.echange_espace_virtuel(eulerian_grad_T_interface_ft_.ghost());
      //
      eulerian_grad_T_interface_ns_.allocate(splitting, IJK_Splitting::ELEM, 1);
      nalloc += 1;
      eulerian_grad_T_interface_ns_.echange_espace_virtuel(eulerian_grad_T_interface_ns_.ghost());
      //
      temperature_ft_.allocate(ref_ijk_ft_->get_splitting_ft(), IJK_Splitting::ELEM, ghost_cells_);
      nalloc += 1;
      temperature_ft_.echange_espace_virtuel(eulerian_grad_T_interface_ft_.ghost());
    }

  compute_eulerian_compo_ = compute_eulerian_compo_ || (liste_post_instantanes_.size() && liste_post_instantanes_.contient_("EULERIAN_COMPO"))
                            || (liste_post_instantanes_.size() && liste_post_instantanes_.contient_("EULERIAN_COMPO_NS"));
  if (compute_eulerian_compo_)
    {
      eulerian_compo_connex_ft_ = &(ref_ijk_ft_->itfce().get_ijk_compo_connex().get_eulerian_compo_connex_ft());
      eulerian_compo_connex_ns_ = &(ref_ijk_ft_->itfce().get_ijk_compo_connex().get_eulerian_compo_connex());
      eulerian_compo_connex_ghost_ft_ = &(ref_ijk_ft_->itfce().get_ijk_compo_connex().get_eulerian_compo_connex_ghost_ft());
      eulerian_compo_connex_ghost_ns_ = &(ref_ijk_ft_->itfce().get_ijk_compo_connex().get_eulerian_compo_connex_ghost());
      eulerian_compo_connex_from_interface_ft_ = &(ref_ijk_ft_->itfce().get_ijk_compo_connex().get_eulerian_compo_connex_from_interface_ft());
      eulerian_compo_connex_from_interface_ns_ = &(ref_ijk_ft_->itfce().get_ijk_compo_connex().get_eulerian_compo_connex_from_interface_ns());
      eulerian_compo_connex_from_interface_ghost_ft_ = &(ref_ijk_ft_->itfce().get_ijk_compo_connex().get_eulerian_compo_connex_from_interface_ghost_ft());
      eulerian_compo_connex_from_interface_ghost_ns_ = &(ref_ijk_ft_->itfce().get_ijk_compo_connex().get_eulerian_compo_connex_from_interface_ghost_ns());
      eulerian_compo_connex_from_interface_int_ns_ = &(ref_ijk_ft_->itfce().get_ijk_compo_connex().get_eulerian_compo_connex_int_from_interface_ns());
      eulerian_compo_connex_from_interface_ghost_int_ns_ = &(ref_ijk_ft_->itfce().get_ijk_compo_connex().get_eulerian_compo_connex_int_from_interface_ghost_ns());
    }

  compute_rising_velocities_ = compute_rising_velocities_ ||
                               (liste_post_instantanes_.size() && liste_post_instantanes_.contient_("RISING_VELOCITIES"));
  fill_rising_velocities_ = compute_rising_velocities_ && (fill_rising_velocities_ ||
                                                           (liste_post_instantanes_.size() && liste_post_instantanes_.contient_("RISING_VELOCITIES")));
  if (fill_rising_velocities_)
    {
      eulerian_rising_velocities_.allocate(splitting, IJK_Splitting::ELEM, 0);
      eulerian_rising_velocities_.data() = 0;
      nalloc += 1;
      eulerian_rising_velocities_.echange_espace_virtuel(eulerian_rising_velocities_.ghost());
    }

  compute_hess_T_elem_ = compute_hess_T_elem_ || liste_post_instantanes_.contient_("HESS_T_ELEM");
  compute_hess_diag_T_elem_ = compute_hess_T_elem_ || compute_hess_diag_T_elem_ || liste_post_instantanes_.contient_("HESS_DIAG_T_ELEM")
                              || liste_post_instantanes_.contient_("HESS_XX_T_ELEM") || liste_post_instantanes_.contient_("HESS_YY_T_ELEM")
                              || liste_post_instantanes_.contient_("HESS_ZZ_T_ELEM");
  compute_hess_cross_T_elem_ = compute_hess_T_elem_ || compute_hess_cross_T_elem_ || liste_post_instantanes_.contient_("HESS_CROSS_T_ELEM")
                               || liste_post_instantanes_.contient_("HESS_XY_T_ELEM") || liste_post_instantanes_.contient_("HESS_XZ_T_ELEM")
                               || liste_post_instantanes_.contient_("HESS_YX_T_ELEM") || liste_post_instantanes_.contient_("HESS_YZ_T_ELEM")
                               || liste_post_instantanes_.contient_("HESS_ZX_T_ELEM") || liste_post_instantanes_.contient_("HESS_ZY_T_ELEM")
                               || liste_post_instantanes_.contient_("HESS_XZ_ZX_T_ELEM") || liste_post_instantanes_.contient_("HESS_ZX_XZ_T_ELEM")
                               || liste_post_instantanes_.contient_("HESS_YZ_ZY_T_ELEM") || liste_post_instantanes_.contient_("HESS_ZY_YZ_T_ELEM")
                               || liste_post_instantanes_.contient_("HESS_XY_YX_T_ELEM") || 	liste_post_instantanes_.contient_("HESS_YX_XY_T_ELEM");

  compute_hess_T_elem_ = compute_hess_diag_T_elem_ && compute_hess_cross_T_elem_;

  compute_grad_T_elem_ = compute_hess_cross_T_elem_ || compute_grad_T_elem_ || liste_post_instantanes_.contient_("GRAD_T_ELEM")
                         || liste_post_instantanes_.contient_("GRAD_T_DIR_X_ELEM") || liste_post_instantanes_.contient_("GRAD_T_DIR_Y_ELEM")
                         || liste_post_instantanes_.contient_("GRAD_T_DIR_Z_ELEM");
  if (compute_grad_T_elem_)
    {
      allocate_cell_vector(grad_T_elem_, splitting, ghost_cells_); // 1 or 0 ?
      nalloc += 3;
      grad_T_elem_.echange_espace_virtuel();
      temperature_grad_op_centre_.initialize(splitting);
    }

  if (compute_hess_diag_T_elem_)
    {
      allocate_cell_vector(hess_diag_T_elem_, splitting, ghost_cells_);  // 1 or 0 ?
      nalloc += 3;
      hess_diag_T_elem_.echange_espace_virtuel();
      temperature_hess_op_centre_.initialize(splitting);
    }

  if (compute_hess_cross_T_elem_)
    {
      allocate_cell_vector(hess_cross_T_elem_, splitting, ghost_cells_);  // 1 or 0 ?
      nalloc += 3;
      hess_cross_T_elem_.echange_espace_virtuel();
      /*
       * TODO: Cross derivatives (adapt the diffusion operator ?)
       * Pb with Finite Volume Op, can not really relate on the fluxes
       * to derive the cross derivatives...
       */
    }
  spherical_approx_ = !spherical_exact_;

  /*
   * FIXME: Temporary need to rewrite IJK_Thermal.cpp posttraitements
   */
  allocate_cell_vector(dummy_int_vect_, splitting, 0); // 1 or 0 ?
  allocate_cell_vector(dummy_double_vect_, splitting, 0); // 1 or 0 ?
  dummy_int_field_.allocate(splitting, IJK_Splitting::ELEM, 0);
  dummy_double_field_.allocate(splitting, IJK_Splitting::ELEM, 0);
  nalloc += 8;
  for (int c=0; c<3; c++)
    {
      dummy_int_vect_[c].data() = 0;
      dummy_int_vect_[c].data() = 0.;
    }
  dummy_int_field_.data() = 0;
  dummy_double_field_.data() = 0;

  // ref_ijk_ft_->redistrib_from_ft_elem().redistribute(eulerian_grad_T_interface_, eulerian_grad_T_interface_);
  // Cout << "End of " << que_suis_je() << "::initialize()" << finl;
  return nalloc;
}

void IJK_Thermal_base::compute_temperature_init()
{
  Cout << "Temperature initialization from expression \nTini = " << expression_T_init_ << finl;
  set_field_data(temperature_, expression_T_init_, ref_ijk_ft_->itfce().I(), 0.);
}

void IJK_Thermal_base::recompute_temperature_init()
{
  Cout << "Temperature initialization from expression \nTini = " << expression_T_init_ << finl;
  set_field_data(temperature_, expression_T_init_, ref_ijk_ft_->itfce().In(), 0.);
}

double IJK_Thermal_base::get_modified_time()
{
  return ref_ijk_ft_->get_current_time();
}

void IJK_Thermal_base::compute_cell_volume()
{
  const IJK_Grid_Geometry& geom = d_temperature_.get_splitting().get_grid_geometry();
  const double dx = geom.get_constant_delta(DIRECTION_I);
  const double dy = geom.get_constant_delta(DIRECTION_J);
  const double dz = geom.get_constant_delta(DIRECTION_K);
  vol_ = dx*dy*dz;
}

void IJK_Thermal_base::compute_min_cell_delta()
{
  const IJK_Grid_Geometry& geom = d_temperature_.get_splitting().get_grid_geometry();
  const double dx = geom.get_constant_delta(DIRECTION_I);
  const double dy = geom.get_constant_delta(DIRECTION_J);
  const double dz = geom.get_constant_delta(DIRECTION_K);
  min_delta_xyz_ = std::min(std::min(dx,dy),dz);
}

void IJK_Thermal_base::compute_cell_diagonal(const IJK_Splitting& splitting)
{
  const IJK_Grid_Geometry& geom = splitting.get_grid_geometry();
  const double dx = geom.get_constant_delta(DIRECTION_I);
  const double dy = geom.get_constant_delta(DIRECTION_J);
  const double dz = geom.get_constant_delta(DIRECTION_K);
  cell_diagonal_ = sqrt(dx*dx + dy*dy + dz*dz);
}

void IJK_Thermal_base::update_thermal_properties()
{
  if (single_phase_)
    cp_vapour_ = 0.;

  const double ene_ini = compute_global_energy();
  const IJK_Field_double& indic = ref_ijk_ft_->itfce().I();

  // Nombre de mailles du domaine NS :
  const int nx = indic.ni();
  const int ny = indic.nj();
  const int nz = indic.nk();
  const double rho_l = ref_ijk_ft_->get_rho_l();
  const double rho_v = ref_ijk_ft_->get_rho_v();
  for (int k=0; k < nz ; k++)
    for (int j=0; j< ny; j++)
      for (int i=0; i < nx; i++)
        {
          double chi_l = indic(i,j,k);
          if (rho_cp_post_)
            rho_cp_(i,j,k) = rho_l*cp_liquid_*chi_l + rho_v*cp_vapour_*(1-chi_l);
          if (calculate_local_energy_)
            rho_cp_T_(i,j,k) = (rho_l*cp_liquid_*chi_l + rho_v*cp_vapour_*(1-chi_l))*temperature_(i,j,k);
        }
  if (rho_cp_post_)
    rho_cp_.echange_espace_virtuel(rho_cp_.ghost());
  if (calculate_local_energy_)
    rho_cp_T_.echange_espace_virtuel(rho_cp_T_.ghost());

  // Semble un endroit approprie pour calculer la variation d'energie due au transport de l'interface:
  const double ene_post = compute_global_energy();
  Cerr << "[Energy-Budget-T"<<rang_<<"-2-TransportIndic] time t=" << ref_ijk_ft_->get_current_time()
       << " " << ene_ini
       << " " << ene_post
       << " delta=" << ene_post-ene_ini << " [W.m-3]." << finl;
}

// Methode de calcul du pas de temps max base sur Fo pour l'equation de thermique.
// CFL value is not computed as it is the same as for the velocity equation.
// The calculation should be stable if Fo <= 1.0 (thanks to the 0.5 in the formula below).
double IJK_Thermal_base::compute_timestep(const double timestep,
                                          const double dxmin)
{
  double alpha_max;
  double rho_l = ref_ijk_ft_->get_rho_l();
  double rho_v= ref_ijk_ft_->get_rho_v();
  if (single_phase_)
    alpha_max = lambda_liquid_ / (rho_l * cp_liquid_);
  else
    alpha_max = std::max(lambda_liquid_ / (rho_l * cp_liquid_), lambda_vapour_ / (rho_v * cp_vapour_));
  dt_fo_ = dxmin * dxmin / (alpha_max + 1.e-20) * fo_ * 0.125; // 1/6 ou 1/8 ?
  if (diff_temperature_negligible_) dt_fo_ = 1.e20;
  return dt_fo_;
}

void IJK_Thermal_base::associer(const IJK_FT_double& ijk_ft)
{
  ref_ijk_ft_ = ijk_ft;
  liste_post_instantanes_ = ijk_ft.get_post().get_liste_post_instantanes();
}

void IJK_Thermal_base::associer_post(const IJK_FT_Post& ijk_ft_post)
{
  ref_ijk_ft_post_ = ijk_ft_post;
}

void IJK_Thermal_base::associer_switch(const Switch_FT_double& ijk_ft_switch)
{
  ref_ijk_ft_switch_ = ijk_ft_switch;
}

void IJK_Thermal_base::associer_interface_intersections(const Intersection_Interface_ijk_cell& intersection_ijk_cell,
                                                        const Intersection_Interface_ijk_face& intersection_ijk_face)
{
  ref_intersection_ijk_cell_ = intersection_ijk_cell;
  ref_intersection_ijk_face_ = intersection_ijk_face;
}

void IJK_Thermal_base::retrieve_ghost_fluid_params(int& compute_distance,
                                                   int& compute_curvature,
                                                   int& n_iter_distance)
{
  compute_distance = compute_distance || compute_distance_;
  compute_curvature = compute_curvature || compute_curvature_;
  n_iter_distance = std::max(n_iter_distance, n_iter_distance_);
}

void IJK_Thermal_base::get_boundary_fluxes(IJK_Field_local_double& boundary_flux_kmin,
                                           IJK_Field_local_double& boundary_flux_kmax)
{
  boundary_flux_kmin = boundary_flux_kmin_;
  boundary_flux_kmax = boundary_flux_kmax_;
}

void IJK_Thermal_base::euler_time_step(const double timestep)
{
  calculer_dT(ref_ijk_ft_->get_velocity());
  // Update the temperature :
  const int kmax = temperature_.nk();
  const double ene_ini = compute_global_energy();
  for (int k = 0; k < kmax; k++)
    {
      ref_ijk_ft_->euler_explicit_update(d_temperature_, temperature_, k);
    }
  /*
   * Erase the temperature increment (second call)
   */
  post_process_after_temperature_increment();

  temperature_.echange_espace_virtuel(temperature_.ghost());
  const double ene_post = compute_global_energy();
  Cerr << "[Energy-Budget-T"<<rang_<<"] time t=" << ref_ijk_ft_->get_current_time()
       << " " << ene_ini
       << " " << ene_post << " [W.m-3]." << finl;
  source_callback();
}

void IJK_Thermal_base::rk3_sub_step(const int rk_step, const double total_timestep,
                                    const double time)
{
  calculer_dT(ref_ijk_ft_->get_velocity());
  // Update the temperature :
  const int kmax = temperature_.nk();
  const double ene_ini = compute_global_energy();
  for (int k = 0; k < kmax; k++)
    {
      runge_kutta3_update(d_temperature_, RK3_F_temperature_, temperature_, rk_step, k, total_timestep);
    }
  temperature_.echange_espace_virtuel(temperature_.ghost());
  const double ene_post = compute_global_energy();
  Cerr << "[Energy-Budget-T"<<rang_<<"] time t=" << ref_ijk_ft_->get_current_time()
       << " " << ene_ini
       << " " << ene_post << " [W.m-3]. [step"<< rk_step << "]" << finl;
  source_callback();
}

void IJK_Thermal_base::sauvegarder_temperature(Nom& lata_name, int idx)
{
  fichier_reprise_temperature_ = lata_name;
  timestep_reprise_temperature_ = 1;
  dumplata_scalar(lata_name, Nom("TEMPERATURE_") + Nom(idx) , temperature_, 0 /*we store a 0 */);
}

/********************************************
 * Protected methods
 ********************************************/


// Mettre rk_step = -1 si schema temps different de rk3.
void IJK_Thermal_base::calculer_dT(const FixedVector<IJK_Field_double, 3>& velocity)
{
  const double current_time = ref_ijk_ft_->get_current_time();
  const double ene_ini = compute_global_energy(d_temperature_);

  /*
   * Clean_subproblems !
   */
  clean_thermal_subproblems();

  // Correct the vapour and mixed cells values
  store_temperature_before_extrapolation();
  correct_temperature_for_eulerian_fluxes();

  /*
   * Correct the temperature field using either the ghost-fluid
   * approach or the laminar sub-resolution approach (and zero values for debug)
   */
  if (debug_)
    Cerr << "Start the Ghost-fluid (GFM) approach" << finl;
  if (debug_)
    Cerr << "Br0 (GFM approach)" << finl;
  compute_eulerian_grad_T_interface();
  if (debug_)
    Cerr << "Br1 (GFM approach)" << finl;
  propagate_eulerian_grad_T_interface();
  if (debug_)
    Cerr << "Br2 (GFM approach)" << finl;
  compute_eulerian_temperature_ghost();
  if (debug_)
    Cerr << "Br3 (GFM approach)" << finl;
  compute_eulerian_bounding_box_fill_compo();
  if (debug_)
    Cerr << "Br4 (GFM approach)" << finl;
  compute_rising_velocities();
  if (debug_)
    Cerr << "End the Ghost-fluid (GFM) approach" << finl;

  /*
   * Compute gradients and hessian of the temperature after the ghost fluid extension
   */
  if (debug_)
    Cerr << "Compute temperature derivatives" << finl;
  compute_temperature_gradient_elem();
  compute_temperature_hessian_diag_elem();
  compute_temperature_hessian_cross_elem();

  /*
   * Compute sub-problems (For Subresolution Child classes only !)
   */
  if (debug_)
    Cerr << "Compute thermal subproblems" << finl;
  compute_thermal_subproblems();

  /*
   * Interpolate a value for the QUICK SCHEME (first call)
   */
  if (debug_)
    Cerr << "Compute temperature mixed cell centres" << finl;
  compute_temperature_cell_centres(0);

  /*
   * Convective and Diffusive fluxes
   */
  if (debug_)
    Cerr << "Compute thermal convective and diffusive fluxes from subproblems" << finl;
  compute_convective_diffusive_fluxes_face_centre();

  if (debug_)
    Cerr << "Prepare ij fluxes" << finl;
  if (!conv_temperature_negligible_ || !diff_temperature_negligible_)
    prepare_ij_fluxes_k_layers();

  /*
   * For post-processing purposes
   */
  enforce_zero_value_eulerian_distance();
  enforce_max_value_eulerian_curvature();
  enforce_max_value_eulerian_field(eulerian_interfacial_area_ft_);

  double nb_diam_upstream_velocity = ref_ijk_ft_->get_nb_diam_upstream();
  if (nb_diam_upstream_ == 0.)
    nb_diam_upstream_ = nb_diam_upstream_velocity;
  if (upstream_temperature_ > -1e20 && ref_ijk_ft_->get_vitesse_upstream() > -1e20)
    force_upstream_temperature(temperature_, upstream_temperature_,
                               ref_ijk_ft_->get_interface(), nb_diam_upstream_,
                               ref_ijk_ft_->get_upstream_dir(), ref_ijk_ft_->get_direction_gravite(),
                               ref_ijk_ft_->get_upstream_stencil());

  compute_temperature_convection(velocity);
  const double ene_postConv = compute_global_energy(d_temperature_);
  add_temperature_diffusion();
  const double ene_postDiffu = compute_global_energy(d_temperature_);
  add_temperature_source();
  const double ene_postSource = compute_global_energy(d_temperature_);

  /*
   * In case of the subresolution or not
   */
  set_zero_temperature_increment();
  // calculer_gradient_temperature(temperature_, grad_T_); Routine Aymeric gradient sur faces

  Cerr << "[Energy-Budget-T"<<rang_<<"-1-TimeResolution] time t=" << current_time
       << " " << ene_ini
       << " " << ene_postConv
       << " " << ene_postDiffu
       << " " << ene_postSource
       << " delta=" << ene_postSource-ene_ini << " [W.m-3]." << finl;

  const IJK_Field_double& T = temperature_;
  double Tmax = -1.e20;
  double Tmin = 1.e20;
  const int ni = T.ni();
  const int nj = T.nj();
  const int nk = T.nk();
  for (int k = 0; k < nk; k++)
    for (int j = 0; j < nj; j++)
      for (int i = 0; i < ni; i++)
        {
          Tmax = std::max(Tmax, T(i,j,k));
          Tmin = std::min(Tmin, T(i,j,k));
        }
  Tmax = Process::mp_max(Tmax);
  Tmin = Process::mp_min(Tmin);
  Cerr <<"[Temperature-MinMax-" << rang_ <<"] t/Tmin/Tmax " << current_time << " "
       << Tmin << " " << Tmax
       << finl;
  return;
}

void IJK_Thermal_base::post_process_after_temperature_increment()
{
  compute_temperature_cell_centres(1);
  enforce_periodic_temperature_boundary_value();
  clip_temperature_values();
  correct_temperature_for_visu();
  set_field_T_ana();
  correct_operators_for_visu();
}

void IJK_Thermal_base::compute_eulerian_distance()
{
  if (compute_distance_)
    {
      // TODO: Do we need to perform an echange_virtuel with interfaces ?
      compute_eulerian_normal_distance_facet_barycentre_field(ref_ijk_ft_->get_interface(),
                                                              eulerian_distance_ft_,
                                                              eulerian_normal_vectors_ft_,
                                                              eulerian_facets_barycentre_ft_,
                                                              n_iter_distance_);
      eulerian_distance_ft_.echange_espace_virtuel(eulerian_distance_ft_.ghost());
      eulerian_distance_ns_.data() = 0.;
      eulerian_distance_ns_.echange_espace_virtuel(eulerian_distance_ns_.ghost());
      ref_ijk_ft_->redistribute_from_splitting_ft_elem(eulerian_distance_ft_, eulerian_distance_ns_);
      eulerian_distance_ns_.echange_espace_virtuel(eulerian_distance_ns_.ghost());
      for(int dir=0; dir<3; dir++)
        {
          ref_ijk_ft_->redistribute_from_splitting_ft_elem(eulerian_normal_vectors_ft_[dir], eulerian_normal_vectors_ns_[dir]);
          ref_ijk_ft_->redistribute_from_splitting_ft_elem(eulerian_facets_barycentre_ft_[dir], eulerian_facets_barycentre_ns_[dir]);
        }
      eulerian_normal_vectors_ns_normed_[0].data() = 0.;
      eulerian_normal_vectors_ns_normed_[1].data() = 0.;
      eulerian_normal_vectors_ns_normed_[2].data() = 0.;
      const int nx = eulerian_normal_vectors_ns_normed_[0].ni();
      const int ny = eulerian_normal_vectors_ns_normed_[0].nj();
      const int nz = eulerian_normal_vectors_ns_normed_[0].nk();
      for (int k=0; k < nz ; k++)
        for (int j=0; j< ny; j++)
          for (int i=0; i < nx; i++)
            {
              double norm_x = eulerian_normal_vectors_ns_[0](i,j,k);
              double norm_y = eulerian_normal_vectors_ns_[1](i,j,k);
              double norm_z = eulerian_normal_vectors_ns_[2](i,j,k);
              norm_x *= norm_x;
              norm_y *= norm_y;
              norm_z *= norm_z;
              const double norm = norm_x + norm_y + norm_z;
              if (norm > 0)
                {
                  eulerian_normal_vectors_ns_normed_[0](i,j,k) = eulerian_normal_vectors_ns_[0](i,j,k) / sqrt(norm);
                  eulerian_normal_vectors_ns_normed_[1](i,j,k) = eulerian_normal_vectors_ns_[1](i,j,k) / sqrt(norm);
                  eulerian_normal_vectors_ns_normed_[2](i,j,k) = eulerian_normal_vectors_ns_[2](i,j,k) / sqrt(norm);
                }
            }
      eulerian_normal_vectors_ns_normed_.echange_espace_virtuel();
    }
  else
    Cerr << "Don't compute the eulerian distance field" << finl;
}

void IJK_Thermal_base::enforce_zero_value_eulerian_distance()
{
  if (compute_distance_)
    {
      enforce_zero_value_eulerian_field(eulerian_distance_ft_);
      enforce_zero_value_eulerian_field(eulerian_distance_ns_);
    }
  else
    Cerr << "Eulerian distance has not been computed" << finl;
}

void IJK_Thermal_base::compute_eulerian_curvature()
{
  if (compute_curvature_)
    {
      /*
       * Laplacian operator may not work properly with FT_field ?
       */
      eulerian_distance_ft_.echange_espace_virtuel(eulerian_distance_ft_.ghost());
      compute_eulerian_curvature_field_from_distance_field(eulerian_distance_ft_,
                                                           eulerian_curvature_ft_,
                                                           boundary_flux_kmin_,
                                                           boundary_flux_kmax_);
      eulerian_curvature_ft_.echange_espace_virtuel(eulerian_curvature_ft_.ghost());
      ref_ijk_ft_->redistribute_from_splitting_ft_elem(eulerian_curvature_ft_, eulerian_curvature_ns_);
    }
  else
    Cerr << "Don't compute the eulerian curvature field" << finl;
}

void IJK_Thermal_base::compute_eulerian_curvature_from_interface()
{
  if (compute_curvature_)
    {
      eulerian_interfacial_area_ft_.echange_espace_virtuel(eulerian_interfacial_area_ft_.ghost());
      eulerian_normal_vectors_ft_.echange_espace_virtuel();
      int nb_groups = ref_ijk_ft_->get_interface().nb_groups();
      // Boucle debute a -1 pour faire l'indicatrice globale.
      // S'il n'y a pas de groupes de bulles (monophasique ou monodisperse), on passe exactement une fois dans la boucle
      if (nb_groups == 1)
        nb_groups = 0; // Quand il n'y a qu'un groupe, on ne posttraite pas les choses pour ce groupe unique puisque c'est identique au cas global
      for (int igroup = -1; igroup < nb_groups; igroup++)
        {
          compute_eulerian_curvature_field_from_interface(eulerian_normal_vectors_ft_,
                                                          ref_ijk_ft_->get_interface(),
                                                          eulerian_interfacial_area_ft_,
                                                          eulerian_curvature_ft_,
                                                          n_iter_distance_,
                                                          igroup);
        }
      eulerian_interfacial_area_ft_.echange_espace_virtuel(eulerian_interfacial_area_ft_.ghost());
      eulerian_curvature_ft_.echange_espace_virtuel(eulerian_curvature_ft_.ghost());
      ref_ijk_ft_->redistribute_from_splitting_ft_elem(eulerian_interfacial_area_ft_, eulerian_interfacial_area_ns_);
      ref_ijk_ft_->redistribute_from_splitting_ft_elem(eulerian_curvature_ft_, eulerian_curvature_ns_);
    }
  else
    Cerr << "Don't compute the eulerian curvature field" << finl;
}

void IJK_Thermal_base::enforce_zero_value_eulerian_curvature()
{
  if (compute_curvature_)
    enforce_zero_value_eulerian_field(eulerian_curvature_ft_);
  else
    Cerr << "Eulerian curvature has not been computed" << finl;
}

void IJK_Thermal_base::enforce_max_value_eulerian_curvature()
{
  if (compute_curvature_)
    {
      enforce_max_value_eulerian_field(eulerian_curvature_ft_);
      enforce_max_value_eulerian_field(eulerian_curvature_ns_);
    }
  else
    Cerr << "Eulerian curvature has not been computed" << finl;
}

void IJK_Thermal_base::compute_eulerian_grad_T_interface()
{
  if (compute_grad_T_interface_)
    {
      temperature_ft_.data() = 0.;
      temperature_ft_.echange_espace_virtuel(temperature_ft_.ghost());
      temperature_.echange_espace_virtuel(temperature_.ghost());
      ref_ijk_ft_->redistribute_to_splitting_ft_elem(temperature_, temperature_ft_);
      temperature_ft_.echange_espace_virtuel(temperature_ft_.ghost());
      eulerian_distance_ft_.echange_espace_virtuel(eulerian_distance_ft_.ghost());
      eulerian_interfacial_area_ft_.echange_espace_virtuel(eulerian_interfacial_area_ft_.ghost());
      eulerian_curvature_ft_.echange_espace_virtuel(eulerian_curvature_ft_.ghost());
      compute_eulerian_normal_temperature_gradient_interface(eulerian_distance_ft_,
                                                             ref_ijk_ft_->itfce().I_ft(),
                                                             eulerian_interfacial_area_ft_,
                                                             eulerian_curvature_ft_,
                                                             temperature_ft_,
                                                             eulerian_grad_T_interface_ft_,
                                                             spherical_approx_);
      /*
       * TODO:
       */
      compute_eulerian_normal_temperature_gradient_interface(eulerian_distance_ns_,
                                                             ref_ijk_ft_->itfce().I(),
                                                             eulerian_interfacial_area_ns_,
                                                             eulerian_curvature_ns_,
                                                             temperature_,
                                                             eulerian_grad_T_interface_ns_,
                                                             spherical_approx_);
    }
  else
    Cerr << "Don't compute the grad_T_interface field" << finl;
}

void IJK_Thermal_base::propagate_eulerian_grad_T_interface()
{
  if (compute_grad_T_interface_)
    {
      eulerian_distance_ft_.echange_espace_virtuel(eulerian_distance_ft_.ghost());
      eulerian_grad_T_interface_ft_.echange_espace_virtuel(eulerian_grad_T_interface_ft_.ghost());
      propagate_eulerian_normal_temperature_gradient_interface(ref_ijk_ft_->itfce(),
                                                               eulerian_distance_ft_,
                                                               eulerian_grad_T_interface_ft_,
                                                               n_iter_distance_,
                                                               gfm_recompute_field_ini_,
                                                               gfm_zero_neighbour_value_mean_,
                                                               gfm_vapour_mixed_only_);
      eulerian_grad_T_interface_ft_.echange_espace_virtuel(eulerian_grad_T_interface_ft_.ghost());

      eulerian_grad_T_interface_ns_.echange_espace_virtuel(eulerian_grad_T_interface_ns_.ghost());
      ref_ijk_ft_->redistribute_from_splitting_ft_elem(eulerian_grad_T_interface_ft_, eulerian_grad_T_interface_ns_);
      /*
       * FIXME: Extrapolate with NS field is not a good idea
       */
      // eulerian_distance_ns_.echange_espace_virtuel(eulerian_distance_ns_.ghost());
      // eulerian_grad_T_interface_ns_.echange_espace_virtuel(eulerian_grad_T_interface_ns_.ghost());
      // propagate_eulerian_normal_temperature_gradient_interface(ref_ijk_ft_->itfce(),
      //                                                          eulerian_distance_ns_,
      //                                                          eulerian_grad_T_interface_ns_,
      //                                                          n_iter_distance_);
    }
  else
    Cerr << "Don't compute the grad_T_interface field" << finl;
}

void IJK_Thermal_base::compute_eulerian_temperature_ghost()
{
  if (ghost_fluid_)
    {
      eulerian_distance_ft_.echange_espace_virtuel(eulerian_distance_ft_.ghost());
      eulerian_curvature_ft_.echange_espace_virtuel(eulerian_curvature_ft_.ghost());
      eulerian_grad_T_interface_ft_.echange_espace_virtuel(eulerian_grad_T_interface_ft_.ghost());
      temperature_ft_.echange_espace_virtuel(temperature_ft_.ghost());
      compute_eulerian_extended_temperature(ref_ijk_ft_->itfce().I_ft(),
                                            eulerian_distance_ft_,
                                            eulerian_curvature_ft_,
                                            eulerian_grad_T_interface_ft_,
                                            temperature_ft_,
                                            spherical_approx_);

      // ref_ijk_ft_->redistribute_from_splitting_ft_elem(temperature_ft_, temperature_);

      eulerian_distance_ns_.echange_espace_virtuel(eulerian_distance_ns_.ghost());
      eulerian_curvature_ns_.echange_espace_virtuel(eulerian_curvature_ns_.ghost());
      eulerian_grad_T_interface_ns_.echange_espace_virtuel(eulerian_grad_T_interface_ns_.ghost());
      temperature_.echange_espace_virtuel(temperature_.ghost());
      compute_eulerian_extended_temperature(ref_ijk_ft_->itfce().I(),
                                            eulerian_distance_ns_,
                                            eulerian_curvature_ns_,
                                            eulerian_grad_T_interface_ns_,
                                            temperature_,
                                            spherical_approx_);
      temperature_.echange_espace_virtuel(temperature_.ghost());
    }
  else
    Cerr << "Don't compute the ghost temperature field" << finl;
}

void IJK_Thermal_base::compute_eulerian_bounding_box_fill_compo()
{
  if (compute_eulerian_compo_)
    {
      bubbles_barycentre_ = ref_ijk_ft_->itfce().get_ijk_compo_connex().get_bubbles_barycentre();
      bubbles_volume_ = ref_ijk_ft_->itfce().get_ijk_compo_connex().get_bubbles_volume();
      bounding_box_= ref_ijk_ft_->itfce().get_ijk_compo_connex().get_bounding_box();
      min_max_larger_box_ = ref_ijk_ft_->itfce().get_ijk_compo_connex().get_min_max_larger_box();
    }
  else
    Cerr << "Don't compute the eulerian bubbles' components (composantes connexes)" << finl;
}

void IJK_Thermal_base::compute_rising_velocities()
{
  if (compute_rising_velocities_)
    {
      int nb_bubbles = ref_ijk_ft_->itfce().get_nb_bulles_reelles();
      rising_velocities_ = ArrOfDouble(nb_bubbles);
      rising_vectors_ = DoubleTab(nb_bubbles, 3);
      compute_rising_velocity(ref_ijk_ft_->get_velocity(), ref_ijk_ft_->itfce(),
                              eulerian_compo_connex_from_interface_int_ns_, ref_ijk_ft_->get_direction_gravite(),
                              rising_velocities_, rising_vectors_,
                              liquid_velocity_);
      // compute_rising_velocity(ref_ijk_ft_->get_velocity(), ref_ijk_ft_->itfce(),
      //                         ref_ijk_ft_->itfce().get_ijk_compo_connex().get_eulerian_compo_connex(),
      //	   										 ref_ijk_ft_->get_direction_gravite(),
      //                         rising_velocities_, rising_vectors_);

      /*
       * FIXME: Use the velocity of the interface instead ?
       */
      if (fill_rising_velocities_)
        {
          eulerian_rising_velocities_.data() = 0.;
          eulerian_rising_velocities_.echange_espace_virtuel(eulerian_rising_velocities_.ghost());
          fill_rising_velocity_int(eulerian_compo_connex_from_interface_int_ns_, rising_velocities_, eulerian_rising_velocities_);
          // fill_rising_velocity_double(eulerian_compo_connex_ns_, rising_velocities_, eulerian_rising_velocities_);
        }
    }
  else
    Cerr << "Don't compute the ghost temperature field" << finl;
}

void IJK_Thermal_base::enforce_zero_value_eulerian_field(IJK_Field_double& eulerian_field)
{
  const int nx = eulerian_field.ni();
  const int ny = eulerian_field.nj();
  const int nz = eulerian_field.nk();
  static const double invalid_distance_value = -1.e30;
  for (int k=0; k < nz ; k++)
    for (int j=0; j< ny; j++)
      for (int i=0; i < nx; i++)
        if (eulerian_field(i,j,k) < invalid_distance_value)
          eulerian_field(i,j,k) = 0.;
}

void IJK_Thermal_base::enforce_max_value_eulerian_field(IJK_Field_double& eulerian_field)
{
  double eulerian_field_max = -1.e20;
  const int nx = eulerian_field.ni();
  const int ny = eulerian_field.nj();
  const int nz = eulerian_field.nk();
  static const double invalid_distance_value = -1.e30;
  for (int k=0; k < nz ; k++)
    for (int j=0; j< ny; j++)
      for (int i=0; i < nx; i++)
        eulerian_field_max = std::max(eulerian_field_max, eulerian_field(i,j,k));
  for (int k=0; k < nz ; k++)
    for (int j=0; j< ny; j++)
      for (int i=0; i < nx; i++)
        if (eulerian_field(i,j,k) < invalid_distance_value)
          eulerian_field(i,j,k) = eulerian_field_max;
}

void IJK_Thermal_base::enforce_min_value_eulerian_field(IJK_Field_double& eulerian_field)
{
  double eulerian_field_min = 1.e20;
  const int nx = eulerian_field.ni();
  const int ny = eulerian_field.nj();
  const int nz = eulerian_field.nk();
  static const double invalid_distance_value = -1.e30;
  for (int k=0; k < nz ; k++)
    for (int j=0; j< ny; j++)
      for (int i=0; i < nx; i++)
        eulerian_field_min = std::min(eulerian_field_min, eulerian_field(i,j,k));
  for (int k=0; k < nz ; k++)
    for (int j=0; j< ny; j++)
      for (int i=0; i < nx; i++)
        if (eulerian_field(i,j,k) < invalid_distance_value)
          eulerian_field(i,j,k) = eulerian_field_min;
}

void IJK_Thermal_base::compute_mixed_cells_number(const IJK_Field_double& indicator)
{
  mixed_cells_number_ = 0.;
  const int nx = indicator.ni();
  const int ny = indicator.nj();
  const int nz = indicator.nk();
  for (int k=0; k < nz ; k++)
    for (int j=0; j< ny; j++)
      for (int i=0; i < nx; i++)
        if (fabs(indicator(i,j,k)) > VAPOUR_INDICATOR_TEST && indicator(i,j,k) < LIQUID_INDICATOR_TEST)
          mixed_cells_number_ += 1;
  Cerr << "There are " << mixed_cells_number_ << "mixed cells." << finl;
}

void IJK_Thermal_base::compute_temperature_gradient_elem()
{
  if (compute_grad_T_elem_)
    {
      temperature_.echange_espace_virtuel(temperature_.ghost());
      for (int dir=0; dir<3; dir++)
        grad_T_elem_[dir].data()=0;
      grad_T_elem_.echange_espace_virtuel();
      temperature_grad_op_centre_.calculer_grad(temperature_, grad_T_elem_);
      grad_T_elem_.echange_espace_virtuel();
    }
  else
    Cerr << "The temperature gradient at the cell centres is not computed" << finl;
}

void IJK_Thermal_base::compute_temperature_hessian_diag_elem()
{
  if (compute_hess_diag_T_elem_)
    {
      temperature_.echange_espace_virtuel(temperature_.ghost());
      for (int dir=0; dir<3; dir++)
        hess_diag_T_elem_[dir].data()=0;
      hess_diag_T_elem_.echange_espace_virtuel();
      temperature_hess_op_centre_.calculer_hess(temperature_, hess_diag_T_elem_,
                                                boundary_flux_kmin_, boundary_flux_kmax_);
      hess_diag_T_elem_.echange_espace_virtuel();
    }
  else
    Cerr << "The temperature gradient at the cell centres is not computed" << finl;
}

void IJK_Thermal_base::compute_temperature_hessian_cross_elem()
{
  if (compute_hess_cross_T_elem_)
    {
      temperature_.echange_espace_virtuel(temperature_.ghost());
      for (int dir=0; dir<3; dir++)
        hess_cross_T_elem_[dir].data()=0;
      hess_cross_T_elem_.echange_espace_virtuel();
      temperature_grad_op_centre_.calculer_grad_z(grad_T_elem_[1], hess_cross_T_elem_[0]);
      temperature_grad_op_centre_.calculer_grad_z(grad_T_elem_[0], hess_cross_T_elem_[1]);
      temperature_grad_op_centre_.calculer_grad_y(grad_T_elem_[0], hess_cross_T_elem_[2]);
      hess_cross_T_elem_.echange_espace_virtuel();
    }
  else
    Cerr << "The temperature gradient at the cell centres is not computed" << finl;
}

// Convect temperature field by velocity.
// The output is stored in d_temperature_ (it is a volume integral over the CV)
void IJK_Thermal_base::compute_temperature_convection(const FixedVector<IJK_Field_double, 3>& velocity)
{
  static Stat_Counter_Id cnt_conv_temp = statistiques().new_counter(1, "FT convection rho");
  statistiques().begin_count(cnt_conv_temp);
  if (conv_temperature_negligible_)
    {
      d_temperature_.data()=0;
      u_T_convective_volume_.data() = 0;
    }
  else
    {
      temperature_convection_op_.calculer(temperature_, velocity[0], velocity[1], velocity[2], d_temperature_);
      const int ni = d_temperature_.ni();
      const int nj = d_temperature_.nj();
      const int nk = d_temperature_.nk();
      for (int k = 0; k < nk; k++)
        for (int j = 0; j < nj; j++)
          for (int i = 0; i < ni; i++)
            {
              d_temperature_(i,j,k) /= vol_ ;
              u_T_convective_volume_(i,j,k) = d_temperature_(i,j,k);
            }
    }
  statistiques().end_count(cnt_conv_temp);
  DebogIJK::verifier("op_conv(rho)", d_temperature_);
  return;
}

void IJK_Thermal_base::add_temperature_diffusion()
{
  if (boundary_conditions_.get_bctype_k_min() == Boundary_Conditions_Thermique::Paroi_Temperature_imposee)
    {
      const double T_paroi_impose_kmin = boundary_conditions_.get_temperature_kmin();
      double lambda_de_t_paroi_kmin = lambda_liquid_;
      // calculer ici div(lambda*grad(T))*volume)
      calculer_flux_thermique_bord(temperature_, lambda_de_t_paroi_kmin,
                                   T_paroi_impose_kmin, boundary_flux_kmin_, 0 /* boundary kmin */);
    }
  else if (boundary_conditions_.get_bctype_k_min() == Boundary_Conditions_Thermique::Paroi_Flux_impose)
    {
      const double flux_paroi_impose_kmin = boundary_conditions_.get_flux_kmin();
      imposer_flux_thermique_bord(temperature_,
                                  flux_paroi_impose_kmin, boundary_flux_kmin_, 0 /* boundary kmin */);
    }
  else
    {
      // Perio... not needed. The flux on the "boundary" (in that case , it's not truely a boundary) will be computed as inside...
    }
  if (boundary_conditions_.get_bctype_k_max() == Boundary_Conditions_Thermique::Paroi_Temperature_imposee)
    {
      const double T_paroi_impose_kmax = boundary_conditions_.get_temperature_kmax();
      double lambda_de_t_paroi_kmax = lambda_liquid_;
      //TODO: calculer ici div(lambda*grad(T))*volume)
      calculer_flux_thermique_bord(temperature_, lambda_de_t_paroi_kmax,
                                   T_paroi_impose_kmax, boundary_flux_kmax_, 1 /* boundary kmax */);
    }
  else if (boundary_conditions_.get_bctype_k_max() == Boundary_Conditions_Thermique::Paroi_Flux_impose)
    {
      const double flux_paroi_impose_kmax = boundary_conditions_.get_flux_kmax();
      imposer_flux_thermique_bord(temperature_,
                                  flux_paroi_impose_kmax, boundary_flux_kmax_, 1 /* boundary kmax */);
      // Cerr << "not coded yet" << finl;
      // Process::exit();
    }
  else
    {
      // Perio... not needed. The flux on the "boundary" (in that case , it's not truely a boundary) will be computed as inside...
    }
  temperature_.echange_espace_virtuel(temperature_.ghost());
  DebogIJK::verifier("temp", temperature_);

  if (!diff_temperature_negligible_)
    {
      // Performance counters:
      static Stat_Counter_Id cnt_diff_temp = statistiques().new_counter(1, "FT diffusion temperature");
      statistiques().begin_count(cnt_diff_temp);
      /*
       * Correct the diffusive fluxes here or in the operator ?
       */
      temperature_diffusion_op_.calculer(temperature_,
                                         div_coeff_grad_T_volume_,
                                         boundary_flux_kmin_,
                                         boundary_flux_kmax_);
      compute_diffusion_increment();
      statistiques().end_count(cnt_diff_temp);
      DebogIJK::verifier("div_coeff_grad_T_volume_", div_coeff_grad_T_volume_);
    }
}

//////////////////////////////////////////
//void IJK_Thermal_base::add_temperature_source() { ; }
double IJK_Thermal_base::compute_rho_cp_u_mean(const IJK_Field_double& vx)
{
  /*
   * By default use only the liquid phase (same for subresolution)
   * Overridden in Onefluid and others
   */
  const double rho_cp = ref_ijk_ft_->get_rho_l() * cp_liquid_;
  return calculer_rho_cp_u_moyen(vx, vx, vx, rho_cp, 2);
}

double IJK_Thermal_base::get_rho_cp_u_ijk(const IJK_Field_double& vx, int i, int j, int k) const
{
  return ref_ijk_ft_->get_rho_l() * cp_liquid_ * vx(i,j,k);
}

void IJK_Thermal_base::add_temperature_source()
{
  static Stat_Counter_Id cnt_source_temp = statistiques().new_counter(1, "FT source temperature");
  statistiques().begin_count(cnt_source_temp);
  // Dans le cas ou les flux entrants et sortants sont identiques :
  // DONE: changer cette condition non adaptee
  if (type_T_source_!="??")
    {
      // int gravity_dir = ref_ijk_ft_->get_direction_gravite();
      const int gravity_dir=0;
      /*
       * Modifications un jour peut-être ?!
       * Adaptation paroi dans une autre direction ?
       */
      const int wall_normal_dir = DIRECTION_K;
      const IJK_Field_double& vx = ref_ijk_ft_->get_velocity()[gravity_dir];
      double rho_cp_u_moy = compute_rho_cp_u_mean(vx);
      const IJK_Splitting& splitting = temperature_.get_splitting();
      const IJK_Grid_Geometry& geom = splitting.get_grid_geometry();
      const double dl = geom.get_constant_delta(gravity_dir);
      const double lwall = geom.get_domain_length(wall_normal_dir) ;
      const double h = lwall / 2.;
      const int nk = d_temperature_.nk();
      const int ni = d_temperature_.ni();
      const int nj = d_temperature_.nj();
      // TODO: faire une methode calculer_rho_cp
      //debut if source = ponderee par la vitesse
      if (type_T_source_=="dabiri")
        {
          const double wall_flux = boundary_conditions_.get_flux_kmax();
          // wall_flux in W.m^{-2}
          const double qw = wall_flux;
          // TODO: Ask Aymeric -> 2 wall with same B.Cs ?
          const double dTm = 2 * qw / rho_cp_u_moy;
          /*
           * Local correction using global dTm
           * TODO: faux, il manque un facteur 1/h pour homogeneite
           * MG: Pas sûr...
           */
          for (int k = 0; k < nk; k++)
            for (int j = 0; j < nj; j++)
              for (int i = 0; i < ni; i++)
                {
//                  const double rho = ref_ijk_ft_->rho_field_(i,j,k);
//                  const double cp = cp_(i,j,k);
//                  const double u = (vx(i,j,k) +vx(i+1,j,k))/2;
//                  const double rho_cp_u = rho*cp*u;
                  /*
                   * TODO: Ask Aymeric : WTF lambda_variable_; at the wall ???
                   */
//                  if(lambda_variable_)
//                    {
//											const double div_lambda = compute_lambda_variations(dl);
//                      const double div_lambda = (lambda_(i+1,j,k)-lambda_(i-1,j,k))/(2*dl);
//                      source_temperature_(i,j,k) = (rho_cp_u - div_lambda) * dTm;
//                    }
//                  else
//                    {
//                      source_temperature_(i,j,k) = rho_cp_u * dTm;
//                    }
                  const double rho_cp_u = get_rho_cp_u_ijk(vx, i, j, k);
                  const double div_lambda = get_div_lambda_ijk(i,j,k) / (2*dl);
                  source_temperature_(i,j,k) = (rho_cp_u - div_lambda) * dTm;
                  // TODO: Each plate generates heat in a volume of half the plates distance ??
                  const double Sc = qw / h ;  // qw / lambda * h;
                  const double source = (source_temperature_(i,j,k) - Sc) * vol_; // -Sc ) * volume;
                  // TODO: faux, que vient faire Sc en plus de source ?
                  d_temperature_(i,j,k) += source;
                }
          //
          calculer_temperature_physique_T(vx, dTm);
          calculer_temperature_adimensionnelle_theta(vx, wall_flux);
          calculer_Nusselt(vx);
          return;
        }
      // TODO: remplacer dabiri par patch_dabiri apres verif
      else if (type_T_source_=="patch_dabiri")
        {

          Cerr << "Type de source : patch_dabiri" << finl;
          const double wall_flux = boundary_conditions_.get_flux_kmax();
          const double qw = wall_flux;
          const double dTm = -2*qw/(2*h*rho_cp_u_moy) ;

          for (int k = 0; k < nk; k++)
            for (int j = 0; j < nj; j++)
              for (int i = 0; i < ni; i++)
                {
                  const double u = (vx(i,j,k) +vx(i+1,j,k))/2;
                  source_temperature_(i,j,k) = (u - get_div_lambda_ijk(i,j,k) / (2*dl)) * dTm;
                  const double source = source_temperature_(i,j,k);
                  d_temperature_(i,j,k) += source;
                }
          //
          calculer_temperature_physique_T(vx, dTm);
          calculer_temperature_adimensionnelle_theta(vx, wall_flux);
          calculer_Nusselt(vx);
          return;
        }
      else if (type_T_source_=="unweighted_dabiri")
        {
          // Sans la ponderation de S par u (S=rhocpu.../<rhocpu>), cela permet d'avoir une source uniforme!
          // Mais ce n'est plus le meme changement de variable, c'est quoi alors?
          // TODO: Que faire de rho_cp en diphasique? Moyen ou local?
          Cerr << "Type de source : unweighted_dabiri" << finl;
          const double wall_flux = boundary_conditions_.get_flux_kmax();
          const double qw = wall_flux;
          const double liquid_fraction = calculer_v_moyen(ref_ijk_ft_->itfce().I());
          const double rho_cp_l = ref_ijk_ft_->get_rho_l() * cp_liquid_;
          const double rho_cp_v = ref_ijk_ft_->get_rho_v() * cp_vapour_;
          const double rhocp_moy =  rho_cp_l*liquid_fraction + rho_cp_v*(1-liquid_fraction);
          const double dTm = -2*qw/(2*h*rhocp_moy) ;
          for (int k = 0; k < nk; k++)
            for (int j = 0; j < nj; j++)
              for (int i = 0; i < ni; i++)
                {
                  if(lambda_variable_)
                    {
                      Cerr << "Veut-on vraiment calculer une partie en lambda spatialement variable??" <<finl;
                      Cerr << "Exit at IJK_Thermal_base::add_temperature_source" << finl;
                      Process::exit();
                    }
                  else
                    {
                      source_temperature_(i,j,k) = dTm;
                    }
                  const double source = source_temperature_(i,j,k);
                  d_temperature_(i,j,k) += source;
                }
          //
          calculer_temperature_physique_T(vx, dTm);
          calculer_temperature_adimensionnelle_theta(vx, wall_flux);
          calculer_Nusselt(vx);
          return;
        }
      // Debut source = SWARM
      else if (type_T_source_=="SWARM")
        {
          // DONE: idem
          double Tv=0;
          double Tl=0;
          double Vv = 0;
          double Vl = 0;

          // calcul de Tv et Tl
          for (int k = 0; k < nk; k++)
            for (int j = 0; j < nj; j++)
              for (int i = 0; i < ni; i++)
                {
                  const double chi = ref_ijk_ft_->itfce().I(i, j, k);
                  const double T = temperature_(i, j, k);
                  Tv += T*(1.-chi);
                  Vv += (1.-chi);
                  Tl += T*chi;
                  Vl += chi;
                }
          Tv = Process::mp_sum(Tv);
          Tl = Process::mp_sum(Tl);
          Vv = Process::mp_sum(Vv);
          Vl = Process::mp_sum(Vl);
          Tv /= Vv;
          Tl /= Vl;
          Cerr << "AY-test_source : " <<  Tv << finl;

          // Calcul de dT a partir de l'expression suivante : dT = k(T - Tm)/(rho*cp)
          const double kv = kv_;
          const double kl = kl_;
          const double T0v = T0v_;
          const double T0l = T0l_;
          const double rho_cp_l = ref_ijk_ft_->get_rho_l() * cp_liquid_;
          const double rho_cp_v = ref_ijk_ft_->get_rho_v() * cp_vapour_;
          for (int k = 0; k < nk; k++)
            for (int j = 0; j < nj; j++)
              for (int i = 0; i < ni; i++)
                {
                  d_source_Tv_(i,j,k) = kv/rho_cp_v * (Tv - T0v);
                  d_source_Tl_(i,j,k) = kl/rho_cp_l * (Tl - T0l);
                }
          Cerr << "AY-test_source1 : " <<  Tv << finl;
          // TODO: Remplacer euler_explicit_update par l'utilisation de timestep_ et utiliser d_source_Tv_ comme une constante
          for (int k = 0; k < nk; k++)
            {
              ref_ijk_ft_->euler_explicit_update(d_source_Tv_, source_temperature_v_, k);
              ref_ijk_ft_->euler_explicit_update(d_source_Tl_, source_temperature_l_, k);
            }
          Cerr << "AY-test_source2 : " <<  Tv << finl;
          for (int k = 0; k < nk; k++)
            for (int j = 0; j < nj; j++)
              for (int i = 0; i < ni; i++)
                {
                  // chi vaut 1. dans le liquide et 0 dans les bulles
                  const double chi = ref_ijk_ft_->itfce().I(i, j, k);
                  source_temperature_(i,j,k) = (1.-chi)*source_temperature_v_(i,j,k) + chi*source_temperature_l_(i,j,k);
                }
          Cerr << "AY-test_source3 : " <<  Tv << finl;
          for (int k = 0; k < nk; k++)
            for (int j = 0; j < nj; j++)
              for (int i = 0; i < ni; i++)
                {
                  d_temperature_(i,j,k) += source_temperature_(i,j,k)*vol_;
                }
          Cerr << "AY-test_source 111 : " <<  source_temperature_(1,1,1) << finl;
          Cerr << "AY-test_source vol : " <<  vol_ << finl;
          Cerr << "source_temp " << " " << d_temperature_(1,1,1) << finl;
          const double current_time = ref_ijk_ft_->get_current_time();

          // GB : Ma comprehension est que ce ne sont pas des champs, mais un scalaire unique
          const double Sl = source_temperature_l_(0,0,0);
          const double Sv = source_temperature_v_(0,0,0);
          const double dSl = d_source_Tl_(0,0,0);
          const double dSv = d_source_Tv_(0,0,0);
          Cerr <<"[ThermalInfo-" << rang_ <<"] t/Tl/Tv/Sl/Sv/dSldt/dSvdt " << current_time << " "
               << Tl << " " << Tv << " "
               << Sl << " " << Sv << " "
               << dSl << " " << dSv
               <<finl;
          /*
           * TODO: M.G c'est quoi ??
           * Si on utilise ca c'est pour remplir dans tous les cas d'utilisation d'une source le champs temperature_physique_T
           * calculer_temperature_physique_T_dummy();
           */
          temperature_physique_T_.data() = 0.;
          return;
        }
    }
  // Dans ce cas la ce ne sont pas des flux thermiques identiques
  else
    {
      Cerr << "no_source_for_temperature" << finl;
      return;
    }
}

void IJK_Thermal_base::source_callback()
{
  if (liste_post_instantanes_.contient_("TEMPERATURE_ADIM_BULLES"))
    calculer_temperature_adim_bulles();
}

/*
 * Aymeric: Renommer pour expliciter qu'il s'agit de la transformation
 * inverse de Kawamura avec le gradient de temperature moyenne
 */
void IJK_Thermal_base::calculer_temperature_physique_T(const IJK_Field_double&  vx, const double dTm)
{
  if (wall_flux_)
    {
      const IJK_Splitting& splitting = temperature_.get_splitting();
      const IJK_Grid_Geometry& geom = splitting.get_grid_geometry();
      double dx =geom.get_constant_delta(DIRECTION_I);
      double origin_x = geom.get_origin(DIRECTION_I) + (dx * 0.5) ;
      const int offset_i = splitting.get_offset_local(DIRECTION_I);

      const int nk = d_temperature_.nk();
      const int ni = d_temperature_.ni();
      const int nj = d_temperature_.nj();

      for (int k = 0; k < nk; k++)
        for (int j = 0; j < nj; j++)
          for (int i = 0; i < ni; i++)
            {
              const double x = (i+ offset_i ) *dx + origin_x;    //MR A GB: ne doit-on pas soustraire par  -dx*0.5???
              temperature_physique_T_(i,j,k) = (x*dTm)-temperature_(i,j,k);
            }
      temperature_physique_T_.echange_espace_virtuel(temperature_physique_T_.ghost());
      DebogIJK::verifier("temperature_physique_T", temperature_physique_T_);
      return;
    }
  else
    {
      Cerr << "No source for the temperature field" << finl;
      return;
    }
}

void IJK_Thermal_base::calculer_temperature_adim_bulles()
{
  const int nk = temperature_.nk();
  const int ni = temperature_.ni();
  const int nj = temperature_.nj();

  // Calcul de moy1 = moy(chi*T+)/moy(chi) et moy2 = moy((1-chi)*T+) / moy(1-chi)
  double Tl = 0.;
  double Tv = 0.;
  double Vl = 0.;
  double Vv = 0.;
  const IJK_Field_double& T = temperature_;
  const IJK_Field_double& chi = ref_ijk_ft_->itfce().I(); // rappel : chi vaut 1. dans le liquide et 0 dans la vapeur

  // assuming uniform mesh : vol_cell=cste.
  for (int k = 0; k < nk; k++)
    for (int j = 0; j < nj; j++)
      for (int i = 0; i < ni; i++)
        {
          Tl += T(i,j,k)*chi(i, j, k);
          Tv += T(i,j,k)*(1.-chi(i, j, k));
          Vl += chi(i, j, k);
          Vv += (1.-chi(i, j, k));
        }
  Tl = Process::mp_sum(Tl);
  Tv = Process::mp_sum(Tv);
  Vv = Process::mp_sum(Vv);
  Vl = Process::mp_sum(Vl);
  Tl /= Vl;
  Tv /= Vv;
  // Calcul de Tl et Tv :
  // const double Tl = Tv0_ / (1 - Tl) * (1 - Tl / (1 + Tv)) / (1 + Tl * Tv / ((1 - Tl) * (1 + Tv)));
  // const double Tv = Tv0_ / (1 + Tv) * (1 + Tv / (1 - Tl)) / (1 + Tv * Tl / ((1 + Tv) * (1 - Tl)));
  const int ntot = temperature_.get_splitting().get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_I)
                   *temperature_.get_splitting().get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_J)
                   *temperature_.get_splitting().get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_K);
  const double time = ref_ijk_ft_->get_current_time();
  Cerr << "Tl_test : time= "<< time << " alpha_l= " << Vl/ntot<< "  TI=" <<Tl*Vl/ntot<< " Tl="<< Tl << finl;
  Cerr << "Tv_test : time= "<< time << " alpha_v= " << Vv/ntot<< " TIv=" <<Tv*Vv/ntot<< " Tv="<< Tv << finl;

  for (int k = 0; k < nk; k++)
    for (int j = 0; j < nj; j++)
      for (int i = 0; i < ni; i++)
        temperature_adim_bulles_(i,j,k) = (temperature_(i,j,k) - Tv0_)/(Tl - Tv);

  //  temperature_adim_bulles_.echange_espace_virtuel(temperature_adim_bulles_.ghost());
  //  TODO: ce qui suit ne devrait pas etre la, mais je le met ici temporairement avant de trouver une meilleure solution

  double E_tot = 0.;
  double E_liq_pure = 0., E_liq = 0;
  double E_vap_pure = 0., E_vap = 0;
  double E_mixt = 0.;
  calculer_energies(E_liq_pure, E_liq,
                    E_vap_pure, E_vap,
                    E_mixt, E_tot);

  /*
   * TODO: voir si on ne doit pas faire mieux, mais a priori les variations de Tl et Tv
   * sont lentes par rapport au reste donc ea devrait aller.
   * DONE: il y a manifestement un pb ici, car on ne peut pas avoir acces e Tv(n+1) encore,
   * donc il faut stocker Tv(n-1)
   */

  // Impression dans le fichier temperature_bulles.out
  if (Process::je_suis_maitre())
    {
      int reset = (!ref_ijk_ft_->get_reprise()) && (ref_ijk_ft_->get_tstep()==0);
      SFichier fic=Ouvrir_fichier(Nom("_source_temperature_")+Nom(rang_)+Nom("_bulles.out"),
                                  "tstep\ttime\tTl\tTv\tEtot\tElpu\tEl\tEvpu\tEv\tEm",
                                  reset);
      // la derivee_acceleration n'est connue que sur le maitre
      fic<< ref_ijk_ft_->get_tstep()<<" "<< ref_ijk_ft_->get_current_time() <<" "<< Tl << " " << Tv << " " << E_tot;
      fic<< " " << E_liq_pure <<  " " << E_liq;
      fic<< " " << E_vap_pure <<  " " << E_vap;
      fic<< " " << E_mixt;
      fic<<finl;
      fic.close();
    }
}

void IJK_Thermal_base::calculer_energies(double& E_liq_pure,
                                         double& E_liq,
                                         double& E_vap_pure,
                                         double& E_vap,
                                         double& E_mixt,
                                         double& E_tot)
{
  const int nk = temperature_.nk();
  const int ni = temperature_.ni();
  const int nj = temperature_.nj();
  const IJK_Field_double& chi = ref_ijk_ft_->itfce().I(); // rappel : chi vaut 1. dans le liquide et 0 dans la vapeur
  const double rhocpl = ref_ijk_ft_->get_rho_l()*cp_liquid_;
  const double rhocpv = ref_ijk_ft_->get_rho_v()*cp_vapour_;
  for (int k = 0; k < nk; k++)
    for (int j = 0; j < nj; j++)
      for (int i = 0; i < ni; i++)
        {
          const double indic=chi(i, j, k);
          const double  rhocpm = indic*rhocpl+(1.-indic)*rhocpv;
          const double T = temperature_(i,j,k);
          E_tot += T * rhocpm;
          E_liq += indic * rhocpl * T;
          E_vap += (1.-indic) * rhocpv * T;
          if (std::fabs(indic)<1.e-8)
            {
              // vap pure
              E_vap_pure += rhocpv * T;
            }
          else if (std::fabs(1.-indic)<1.e-8)
            {
              // liq pure
              E_liq_pure += rhocpl * T;
            }
          else
            {
              // mixte :
              E_mixt += rhocpm * T;
            }
        }

  const int ntot = (temperature_.get_splitting().get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_I)
                    * temperature_.get_splitting().get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_J)
                    * temperature_.get_splitting().get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_K));
  E_vap_pure = Process::mp_sum(E_vap_pure)/ntot;
  E_liq_pure = Process::mp_sum(E_liq_pure)/ntot;
  E_vap = Process::mp_sum(E_vap)/ntot;
  E_liq = Process::mp_sum(E_liq)/ntot;
  E_tot = Process::mp_sum(E_tot)/ntot;
  E_mixt = Process::mp_sum(E_mixt)/ntot;
}

void IJK_Thermal_base::calculer_source_temperature_ana()
{
  if (liste_post_instantanes_.contient_("ECART_SOURCE_TEMPERATURE_ANA"))
    {
      if (!liste_post_instantanes_.contient_("SOURCE_TEMPERATURE_ANA"))
        set_field_data(source_temperature_ana_, expression_source_temperature_, ref_ijk_ft_->get_velocity()[0], ref_ijk_ft_->get_current_time());
      // do some work

      double ct = ref_ijk_ft_->get_current_time();
      Cerr << "MR: ERROR SOURCE T FIELD " << ct;
      double err = 0.;
      //set_field_data(source_temperature_ana_, curseur->expression_source_T_ana_, ct);
      const int ni = source_temperature_.ni();
      const int nj = source_temperature_.nj();
      const int nk = source_temperature_.nk();
      const int ntot=Process::mp_sum(ni*nj*nk);
      // La temperature est definie a une constante pres:
      // const double cst_temp = temperature_ana_(0,0,0) - curseur->temperature_(0,0,0);
      for (int k = 0; k < nk; k++)
        for (int j = 0; j < nj; j++)
          for (int i = 0; i < ni; i++)
            {
              const double val = source_temperature_ana_(i,j,k) - source_temperature_(i,j,k); //- cst_temp;
              ecart_source_t_ana_(i,j,k) = val;
              err += val*val;
            }
      err=Process::mp_sum(err);
      err=sqrt(err/ntot);
      Cerr << " " << err ;
      if (!Process::je_suis_maitre())
        {
          Process::Journal() << "IJK_FT_Post::posttraiter_champs_instantanes : Champ ECART_SOURCE_TEMPERATURE_ANA sur ce proc (ni,nj,nk,ntot):"
                             << " " << ni << " " << nj << " " << nk << " " << ntot << finl;
        }
      // ecart_t_ana_.echange_espace_virtuel(ecart_t_ana_.ghost());
      Cerr << finl ;
      // n++,dumplata_scalar(lata_name,"ECART_SOURCE_TEMPERATURE_ANA", ecart_source_t_ana_, latastep);
    }
}

double IJK_Thermal_base::compute_variable_wall_temperature(const int kmin, const int kmax)
{
  return calculer_variable_wall(temperature_, temperature_, temperature_, cp_liquid_ * ref_ijk_ft_->get_rho_l(), kmin, kmax, 2);
}

void IJK_Thermal_base::calculer_temperature_adimensionnelle_theta(const IJK_Field_double& vx, const double wall_flux)
{
  /*
   * TODO : M.G -> Don't understand anything ask Aymeric
   */
  if(wall_flux_)
    {
      const int wall_normal_dir = DIRECTION_K;
      const IJK_Splitting& splitting = temperature_.get_splitting();
      const IJK_Grid_Geometry& geom = splitting.get_grid_geometry();
      const double lwall = geom.get_domain_length(wall_normal_dir);
      const double h = lwall/2.;
      const double q_w = wall_flux;
      const int kmin = temperature_.get_splitting().get_offset_local(wall_normal_dir);
      const int kmax = splitting.get_grid_geometry().get_nb_elem_tot(wall_normal_dir);
      const int nk = temperature_.nk();
      const int ni = temperature_.ni();
      const int nj = temperature_.nj();
      double T_wall = compute_variable_wall_temperature(kmin, kmax);
      /*   if (Process::je_suis_maitre())
           {
           T_wall = calculer_variable_wall(temperature_, cp_, ref_ijk_ft_->rho_field_, kmin, kmax);
           Cerr << "calcul de T_wall sur maitre" << finl;
           }
           envoyer_broadcast(T_wall, 0); */
      /*
       * if(kmin+ nk == kmax) //|| (kmin ==0)
       {
       rank = Process::me();
       T_wall = calculer_variable_wall(temperature_, cp_, ref_ijk_ft_->rho_field_, kmin, kmax);
       envoyer_broadcast(T_wall, rank);
       }
       */
      for (int k = 0; k < nk; k++)
        {
          //  const double T_mean = compute_spatial_mean(vx, temperature_, cp_, ref_ijk_ft_->rho_field_, kmin, nktot, k);
          for (int j = 0; j < nj; j++)
            for (int i = 0; i < ni; i++)
              {
                double theta = T_wall - temperature_(i,j,k);
                // double theta = T_wall -T_mean;
                // temperature_adimensionnelle_theta_(i,j,k) = theta/theta_tau;
                const double lambda_l = lambda_liquid_;
                temperature_adimensionnelle_theta_(i,j,k) = theta * lambda_l / q_w / h;
                //    temperature_adimensionnelle_theta_(i,j,k) = -T_mean*lambda_l/q_w/h;
              }

        }
      return;
    }
  else
    {
      Cerr << "no_source_for_temperature" << finl;
      return;
    }
}

double IJK_Thermal_base::compute_temperature_dimensionless_theta_mean(const IJK_Field_double& vx)
{
  return calculer_temperature_adimensionnelle_theta_moy(vx, temperature_adimensionnelle_theta_, vx, vx, ref_ijk_ft_->get_rho_l() * cp_liquid_, 2);
}

void IJK_Thermal_base::calculer_Nusselt(const IJK_Field_double& vx)
{
  const double theta_adim_moy = compute_temperature_dimensionless_theta_mean(vx);
  double Nu = 0.;
  if (std::fabs(theta_adim_moy)>1.e-10)
    Nu = 2./theta_adim_moy;
  const double rho_cp_u_moy = compute_rho_cp_u_mean(vx);
  // Impression dans le fichier source_temperature.out
  if (Process::je_suis_maitre())
    {
      int reset = (!ref_ijk_ft_->get_reprise()) && (ref_ijk_ft_->get_tstep() == 0);
      const Nom name = Nom("_temperature_") + Nom(rang_) + Nom(".out");
      SFichier fic = Ouvrir_fichier(name,
                                    "tstep\ttime\ttheta_adim_moy\tNu\trho_cp_u",
                                    reset);
      fic << ref_ijk_ft_->get_tstep() << " " << ref_ijk_ft_->get_current_time() << " " << theta_adim_moy << " " << Nu << " " << rho_cp_u_moy  << finl;
      fic.close();
    }
}

void IJK_Thermal_base::set_field_T_ana()
{
  Cerr << "Setting analytical temperature "<< rang_ <<" field to "<< expression_T_ana_ << finl;
  set_field_data(temperature_ana_, expression_T_ana_, ref_ijk_ft_->get_current_time());
}

void IJK_Thermal_base::calculer_ecart_T_ana()
{
  if (liste_post_instantanes_.contient_("ECART_T_ANA"))
    {
      if (!liste_post_instantanes_.contient_("TEMPERATURE_ANA"))
        {
          set_field_data(temperature_ana_, expression_T_ana_, ref_ijk_ft_->get_current_time());
        }
      // do some work

      double ct = ref_ijk_ft_->get_current_time();
      Cerr << "GB: ERROR T FIELD " << ct;
      double err = 0.;
      set_field_data(temperature_ana_, expression_T_ana_, ct);
      const int ni = temperature_.ni();
      const int nj = temperature_.nj();
      const int nk = temperature_.nk();
      const int ntot=Process::mp_sum(ni*nj*nk);
      // La temperature est definie a une constante pres:
      // const double cst_temp = temperature_ana_(0,0,0) - curseur->temperature_(0,0,0);
      for (int k = 0; k < nk; k++)
        for (int j = 0; j < nj; j++)
          for (int i = 0; i < ni; i++)
            {
              const double val =   temperature_ana_(i,j,k) - temperature_(i,j,k); //- cst_temp;
              ecart_t_ana_(i,j,k) = val;
              err += val*val;
            }
      err=Process::mp_sum(err);
      err=sqrt(err/ntot);
      Cerr << " " << err ;
      if (!Process::je_suis_maitre())
        {
          Process::Journal() << "IJK_FT_Post::posttraiter_champs_instantanes : Champ ECART_T_ANA sur ce proc (ni,nj,nk,ntot):"
                             << " " << ni << " " << nj << " " << nk << " " << ntot << finl;
        }
      ecart_t_ana_.echange_espace_virtuel(ecart_t_ana_.ghost());
      Cerr << finl ;
      //  n++,dumplata_scalar(lata_name,"ECART_T_ANA", ecart_t_ana_, latastep);
    }
}

void IJK_Thermal_base::calculer_gradient_temperature(const IJK_Field_double& temperature, FixedVector<IJK_Field_double, 3>& grad_T)
{
  if ((liste_post_instantanes_.contient_("GRAD_T") || (calulate_grad_T_)))
    {
      /*
       * Re-initialisation of the gradient vector
       */
      for (int dir = 0; dir < 3; dir++)
        grad_T[dir].data() = 0.;

      //  add_gradient_temperature(temperature, 1. /*constant*/,  grad_T[0], grad_T[1], grad_T[2], boundary_conditions_, lambda_);
      for (int dir = 0; dir < 3; dir++)
        grad_T[dir].echange_espace_virtuel(1);
    }
}

// Results are intensive (ie prop to area)
// Method fills storage_ so it changes the class
// Les interfaces connaissent le splitting_ft_ donc la correspondance doit etre appliquee au splitting ft pour convertir :
// convert_packed_to_ijk_cell.
// Donc il faut un champ de T etendu...

double IJK_Thermal_base::compute_global_energy(const IJK_Field_double& temperature)
{
  global_energy_ = 0.;
  const IJK_Field_double& indic = ref_ijk_ft_->itfce().I();
  const double rhocpl = get_rhocp_l();
  const double rhocpv = get_rhocp_v();
  const int nx = temperature.ni();
  const int ny = temperature.nj();
  const int nz = temperature.nk();
  // To be sure we're on a regular mesh
  assert(indic.get_splitting().get_grid_geometry().get_constant_delta(DIRECTION_K) >0);
  for (int k=0; k < nz ; k++)
    for (int j=0; j< ny; j++)
      for (int i=0; i < nx; i++)
        {
          double chi_l = indic(i,j,k);
          global_energy_ += (rhocpl * chi_l + (1.- chi_l) * rhocpv) * temperature(i,j,k);
        }
  const int ntot = temperature.get_splitting().get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_I)
                   *temperature.get_splitting().get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_J)
                   *temperature.get_splitting().get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_K);
  global_energy_ = mp_sum(global_energy_)/(double)(ntot);
  return global_energy_;
}

/*
 * Getters and setters
 */

double IJK_Thermal_base::get_rhocp_l() const
{
  return  ref_ijk_ft_->get_rho_l() * cp_liquid_;
}

double IJK_Thermal_base::get_rhocp_v() const
{
  return  ref_ijk_ft_->get_rho_v() * cp_vapour_;
}

/*
 * Methods that do not belong to the class
 */

// From DNS_QC; Vectorize code later?
int IJK_Thermal_base::calculer_k_pour_bord(const IJK_Field_double& temperature, const bool bord_kmax)
{
  const int kmin = temperature.get_splitting().get_offset_local(DIRECTION_K);
  const int nktot = temperature.get_splitting().get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_K);
  int k;
  // calcul l'indice k de la couche de mailles voisine du bord. Si je n'ai pas de bord, on met k = -1
  if (!bord_kmax)
    {
      // on veut le bord "k_global = 0"
      if (kmin == 0)
        {
          // ce bord est chez moi... et il est en k=0
          k = 0;
        }
      else
        {
          // ce bord n'est pas chez moi
          k = -1;
        }
    }
  else
    {
      // on veut le bord kmax
      if (kmin + temperature.nk() == nktot)
        {
          // ce bord est chez moi... et il est en k= truc...
          k = temperature.nk() - 1;
        }
      else
        {
          k = -1;
        }
    }
  return k;
}

// From DNS_QC; Vectorize code later?
// valeur de retour: indice local du plan de temperature voisin utilise,
//  -1 si on n'a pas le bord sur ce processeur
// Calcule l'integrale sur chaque face du bord demande du flux de chaleur a travers la face
// positif si le flux va vers les k positifs.
int IJK_Thermal_base::calculer_flux_thermique_bord(const IJK_Field_double& temperature,
                                                   const double lambda_de_t_paroi,
                                                   const double T_paroi_impose,
                                                   IJK_Field_local_double& flux_bord,
                                                   const bool bord_kmax)
{
  const int kmin = temperature.get_splitting().get_offset_local(DIRECTION_K);
  int k = calculer_k_pour_bord(temperature, bord_kmax);
  if (k == -1)
    return k;

  // redimensionne flux_bord avec ni * nj:
  const int ni = temperature.ni(); // nombre d'element local sur ce processeur
  const int nj = temperature.nj();
  flux_bord.allocate(ni, nj, 1 /* 1 seule couche de valeurs en k */, 0 /* pas d'elements fantomes */);

  const IJK_Grid_Geometry& geometry = temperature.get_splitting().get_grid_geometry();
  const double delta_k = geometry.get_delta(DIRECTION_K)[k + kmin]; // k+kmin est l'indice global de la maille locale k
  double facteur = 2.0 / delta_k * geometry.get_constant_delta(DIRECTION_I) * geometry.get_constant_delta(DIRECTION_J);
  if (bord_kmax)
    facteur *= -1.; // he he... je vous laisse reflechir a ca :)
  // nan c'est pas simpa: la convention dans l'operateur de diffusion est
  // d/dt = flux(i,j) - flux(i+1,j) + ... + flux(i,j,k) - flux(i,j,k+1)
  // donc si la paroi inferieure (k=0) est plus froide que le fluide, il faut que le flux stocke soit negatif.
  // et si la paroi inferieure (k=kmax) est plus froide que le fluide, il faut que le flux stocke soit positif.
  for (int j = 0; j < nj; j++)
    {
      for (int i = 0; i < ni; i++)
        {
          // Temperature de la maille voisine
          const double t = temperature(i,j,k);
          // le flux est positif s'il va vers les k croissants
          flux_bord(i,j,0) = (T_paroi_impose - t) * lambda_de_t_paroi * facteur;
        }
    }
  return k;
}
int IJK_Thermal_base::imposer_flux_thermique_bord(const IJK_Field_double& temperature,
                                                  const double flux_paroi_impose,
                                                  IJK_Field_local_double& flux_bord,
                                                  const bool bord_kmax)
{
  int k = calculer_k_pour_bord(temperature, bord_kmax);
  if (k == -1)
    return k;

  // redimensionne flux_bord avec ni * nj:
  const int ni = temperature.ni(); // nombre d'element local sur ce processeur
  const int nj = temperature.nj();
  flux_bord.allocate(ni, nj, 1 /* 1 seule couche de valeurs en k */, 0 /* pas d'elements fantomes */);
  // MR je multiplie le flux par la surface dxdy
  const IJK_Grid_Geometry& geometry = temperature.get_splitting().get_grid_geometry();

  double facteur = 1.* geometry.get_constant_delta(DIRECTION_I) * geometry.get_constant_delta(DIRECTION_J);
  if (bord_kmax)
    facteur *= -1.; // he he... je vous laisse reflechir a ca :)
  // nan c'est pas sympa: la convention dans l'operateur de diffusion est
  // d/dt = flux(i,j) - flux(i+1,j) + ... + flux(i,j,k) - flux(i,j,k+1)
  // donc si la paroi inferieure (k=0) est plus froide que le fluide, il faut que le flux stocke soit negatif.
  // et si la paroi superieure (k=kmax) est plus froide que le fluide, il faut que le flux stocke soit positif.
  for (int j = 0; j < nj; j++)
    {
      for (int i = 0; i < ni; i++)
        {
          // le flux est positif s'il va vers les k croissants
          flux_bord(i,j,0) = flux_paroi_impose * facteur;
        }
    }
  return k;
}

/*
 * Patch coming from IJK_Thermique (Aymeric's thesis)
 */

void IJK_Thermal_base::euler_rustine_step(const double timestep, const double dE)
{
  compute_dT_rustine(dE);
  // Update the temperature :
  const int kmax = temperature_.nk();
  const double ene_ini = compute_global_energy();
  for (int k = 0; k < kmax; k++)
    ref_ijk_ft_->euler_explicit_update(d_T_rustine_, temperature_, k);
  temperature_.echange_espace_virtuel(temperature_.ghost());
  const double ene_post = compute_global_energy();
  Cerr << "[Energy-Budget-T"<<rang_<<" euler rustine] time t=" << ref_ijk_ft_->get_current_time()
       << " " << ene_ini
       << " " << ene_post << " [W.m-3]."
       << " dE "<< dE
       << finl;
  source_callback();
}

void IJK_Thermal_base::compute_dT_rustine(const double dE)
{
  const int ni = T_rust_.ni();
  const int nj = T_rust_.nj();
  const int nk = T_rust_.nk();
  const double rho_l = ref_ijk_ft_->get_rho_l();
  const double rho_v = ref_ijk_ft_->get_rho_v();
  const IJK_Field_double indic = ref_ijk_ft_->itfce().I();
  double int_rhocpTrust = 0;
  for (int k = 0; k < nk; k++)
    for (int j = 0; j < nj; j++)
      for (int i = 0; i < ni; i++)
        {
          int_rhocpTrust +=  (rho_l*cp_liquid_*indic(i,j,k) + rho_v*cp_vapour_*(1.-indic(i,j,k)))*T_rust_(i,j,k);
        }
  const int ntot = T_rust_.get_splitting().get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_I)
                   *T_rust_.get_splitting().get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_J)
                   *T_rust_.get_splitting().get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_K);
  int_rhocpTrust = mp_sum(int_rhocpTrust)/(double)(ntot);
  Cerr << "Le coeff de manque d'energie dE/int_rhocpTrust vaut : " << dE/int_rhocpTrust << finl;
  if (int_rhocpTrust)
    {
      for (int k = 0; k < nk; k++)
        for (int j = 0; j < nj; j++)
          for (int i = 0; i < ni; i++)
            {
              d_T_rustine_(i,j,k) = dE/ int_rhocpTrust * T_rust_(i,j,k);
            }
    }
}

void IJK_Thermal_base::rk3_rustine_sub_step(const int rk_step, const double total_timestep,
                                            const double fractionnal_timestep, const double time, const double dE)
{
  compute_dT_rustine(dE);
  // Update the temperature :
  const int kmax = temperature_.nk();
  const double ene_ini = compute_global_energy();
  for (int k = 0; k < kmax; k++)
    {
      runge_kutta3_update(d_T_rustine_, RK3_F_rustine_, temperature_, rk_step, k, total_timestep);
    }
  temperature_.echange_espace_virtuel(temperature_.ghost());
  const double ene_post = compute_global_energy();
  Cerr << "[Energy-Budget-T"<<rang_<<"RK3 rustine step "<<rk_step<<"] time t=" << ref_ijk_ft_->get_current_time()
       << " " << ene_ini
       << " " << ene_post << " [W.m-3]. [step"<< rk_step << "]" << finl;
  source_callback();
}

void IJK_Thermal_base::compute_interfacial_temperature2(ArrOfDouble& interfacial_temperature,
                                                        ArrOfDouble& flux_normal_interp)
{
  const IJK_Grid_Geometry& geom = ref_ijk_ft_->get_geometry();
  const double dist = 1.52 * std::pow(std::pow(geom.get_constant_delta(0), 2.) +
                                      std::pow(geom.get_constant_delta(1), 2.) +
                                      std::pow(geom.get_constant_delta(2), 2.),
                                      0.5);
  const Maillage_FT_IJK& maillage = ref_ijk_ft_->itfce().maillage_ft_ijk();
  const DoubleTab& normale_facettes = maillage.get_update_normale_facettes();
  const ArrOfDouble& surface_facettes = maillage.get_update_surface_facettes();
  const int nb_facettes = maillage.nb_facettes();
  const IntTab& facettes = maillage.facettes();
  const DoubleTab& sommets = maillage.sommets();
  DoubleTab coord_facettes;
  coord_facettes.resize(nb_facettes, 3);
  coord_facettes = 0.;
  for (int fa7 = 0; fa7 < nb_facettes; fa7++)
    for (int dir = 0; dir < 3; dir++)
      for (int som = 0; som < 3; som++)
        coord_facettes(fa7, dir) += sommets(facettes(fa7, som), dir);

  for (int fa7 = 0; fa7 < nb_facettes; fa7++)
    for (int dir = 0; dir < 3; dir++)
      coord_facettes(fa7, dir) /= 3.;

  DoubleTab coo_liqu, coo_vap;
  ArrOfDouble temp_liqu, temp_vap;
  temp_liqu.set_smart_resize(1);
  temp_vap.set_smart_resize(1);
  coo_liqu.set_smart_resize(1);
  coo_vap.set_smart_resize(1);
  corrige_flux_.calcul_temperature_flux_interface(temperature_ft_,
                                                  lambda_liquid_,
                                                  lambda_vapour_,
                                                  dist,
                                                  coord_facettes,
                                                  normale_facettes,
                                                  interfacial_temperature,
                                                  flux_normal_interp,
                                                  temp_liqu,
                                                  temp_vap,
                                                  coo_liqu,
                                                  coo_vap);
  for (int fa7 = 0; fa7 < nb_facettes; fa7++)
    {
      flux_normal_interp(fa7) *= surface_facettes(fa7);
      interfacial_temperature(fa7) *= surface_facettes(fa7);
    }
}

void IJK_Thermal_base::force_upstream_temperature(IJK_Field_double& temperature, double T_imposed,
                                                  const IJK_Interfaces& interfaces, double nb_diam,
                                                  int upstream_dir, int gravity_dir,
                                                  int upstream_stencil)
{
  int dir = 0;
  if (upstream_dir == -1)
    {
      dir = gravity_dir;
      if (dir == -1)
        dir=0;
    }
  const IJK_Splitting& splitting = temperature.get_splitting();
  const IJK_Grid_Geometry& geom = splitting.get_grid_geometry();

  bool perio =  geom.get_periodic_flag(dir);

  assert(interfaces.get_nb_bulles_reelles() == 1);
  DoubleTab bounding_box;
  // interfaces.calculer_bounding_box_bulles(bounding_box);
  bounding_box = interfaces.get_ijk_compo_connex().get_bounding_box();
  // Calcule la hauteur en x de la permiere bulle et la position de son cdg :
  const double Dbdir = bounding_box(0, dir, 1) - bounding_box(0, dir, 0);
  const double dirb  = ( bounding_box(0, dir, 1) + bounding_box(0, dir, 0) ) / 2.;
  const double ldir = geom.get_domain_length(dir) ;
  if (nb_diam == 0.)
    nb_diam = (ldir/Dbdir) / 2;
  double dirobj = dirb + nb_diam*Dbdir;

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

  const double x2 = (dirobj-origin_dir)/ ddir;
  int index_dir = (int)(floor(x2)) - offset_dir; // C'est l'index local, donc potentiellement negatif...
  int ndir;
  switch(dir)
    {
    case 0:
      ndir = temperature.ni();
      break;
    case 1:
      ndir = temperature.nj();
      break;
    case 2:
      ndir = temperature.nk();
      break;
    default:
      ndir = temperature.ni();
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
    imposed[dir] = T_imposed;
    for (int direction = 0; direction < 3; direction++)
      {
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
            jmax = temperature.nj();
            kmax = temperature.nk();
            break;
          case 1:
            imin = 0;
            jmin = index_dir;
            kmin = 0;
            imax = temperature.ni();
            jmax = jmin+upstream_stencil;
            kmax = temperature.nk();
            break;
          case 2:
            imin = 0;
            jmin = 0;
            kmin = index_dir;
            imax = temperature.ni();
            jmax = temperature.nj();
            kmax = kmin+upstream_stencil;
            break;
          default:
            imin = index_dir;
            jmin = 0;
            kmin = 0;
            imax = imin+upstream_stencil;
            jmax = temperature.nj();
            kmax = temperature.nk();
            break;
          }
        for (int k = kmin; k < kmax; k++)
          for (int j = jmin; j < jmax; j++)
            for (int i = imin; i < imax; i++)
              temperature(i,j,k) = imposed[direction];
      }
  }
}
