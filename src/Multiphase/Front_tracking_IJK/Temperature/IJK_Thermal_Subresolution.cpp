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
// File      : IJK_Thermal_Subresolution.cpp
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/Temperature
//
/////////////////////////////////////////////////////////////////////////////

#include <IJK_Thermal_Subresolution.h>
#include <Param.h>
#include <IJK_Navier_Stokes_tools.h>
#include <DebogIJK.h>
#include <stat_counters.h>
#include <IJK_FT.h>
#include <Corrige_flux_FT.h>
#include <OpConvDiscQuickIJKScalar.h>
#include <IJK_Ghost_Fluid_tools.h>

Implemente_instanciable_sans_constructeur( IJK_Thermal_Subresolution, "IJK_Thermal_Subresolution", IJK_Thermal_base ) ;

IJK_Thermal_Subresolution::IJK_Thermal_Subresolution()
{
  disable_mixed_cells_increment_=1;
  enable_mixed_cells_increment_=0;
  allow_temperature_correction_for_visu_=0;
  subproblem_temperature_extension_=0;
  disable_subresolution_=0;
  convective_flux_correction_ = 0;
  diffusive_flux_correction_ = 0;
  impose_fo_flux_correction_ = 1;
  override_vapour_mixed_values_ = 0;
  compute_eulerian_compo_ = 1;

  points_per_thermal_subproblem_ = 32;
  dr_ = 0.;
  probe_length_=0.;

  boundary_condition_interface_ = -1;
  interfacial_boundary_condition_value_ = 0.;
  impose_boundary_condition_interface_from_simulation_=0;
  boundary_condition_end_ = 0;
  end_boundary_condition_value_ = -1;
  impose_user_boundary_condition_end_value_=0;

  source_terms_type_=2;
  source_terms_correction_=0;

  fd_solver_type_ = "Solv_Gen";
  fd_solvers_ = Motcles(4);
  {
    fd_solvers_[0] = "Solv_Gen";
    fd_solvers_[1] = "Solv_GCP";
    fd_solvers_[2] = "Solv_Cholesky";
    fd_solvers_[3] = "Solv_Cholesky";
  }
  fd_solvers_jdd_ = Motcles(4);
  {
    fd_solvers_jdd_[0] = "thermal_fd_solver_1";
    fd_solvers_jdd_[1] = "thermal_fd_solver_2";
    fd_solvers_jdd_[2] = "thermal_fd_solver_3";
    fd_solvers_jdd_[3] = "thermal_fd_solver_4";
  }
  fd_solver_rank_ = 0;
  // one_dimensional_advection_diffusion_thermal_solver_.nommer("finite_difference_solver");
  // one_dimensional_advection_diffusion_thermal_solver_.typer("Solv_GCP");
  discrete_integral_ = 0;
  quadtree_levels_ = 1;
  advected_frame_of_reference_=0;
  neglect_frame_of_reference_radial_advection_=1;
  approximate_temperature_increment_=0;

  is_first_time_step_ = false;
  first_time_step_temporal_ = 0;
  first_time_step_implicit_ = 0;
  first_time_step_explicit_ = 1;
  local_fourier_ = 1.;
  local_cfl_ = 1.;

  distance_cell_faces_from_lrs_=0.;

  pre_initialise_thermal_subproblems_list_ = 0;
  pre_factor_subproblems_number_ = 3.;

  correct_temperature_cell_neighbours_first_iter_ = 0;
  find_temperature_cell_neighbours_ = 0;
  correct_neighbours_using_probe_length_ = 0;
  neighbours_corrected_rank_ = 1;
  neighbours_colinearity_weighting_ = 0.;
  use_temperature_cell_neighbours_ = 0;

  clip_temperature_values_ = 0;
  disable_post_processing_probes_out_files_ = 0;
  enforce_periodic_boundary_value_ = 0;

  store_cell_faces_corrected_ = 0;

  find_cell_neighbours_for_fluxes_spherical_correction_ = 0;
  use_cell_neighbours_for_fluxes_spherical_correction_ = 0;
  find_reachable_fluxes_ = 0;
  use_reachable_fluxes_ = 0;
}

Sortie& IJK_Thermal_Subresolution::printOn( Sortie& os ) const
{
  IJK_Thermal_base::printOn( os );
  os << "  {\n";
  os << "    convective_flux_correction" <<  " " << convective_flux_correction_ << "\n";
  os << "    diffusive_flux_correction" <<  " " << diffusive_flux_correction_ << "\n";
  os << "    override_vapour_mixed_values" <<  " " << override_vapour_mixed_values_ << "\n";
  os << "  \n}";
  return os;
}

Entree& IJK_Thermal_Subresolution::readOn( Entree& is )
{
  IJK_Thermal_base::readOn( is );
  if (ghost_fluid_)
    override_vapour_mixed_values_ = 1;
  if (!ghost_fluid_)
    approximate_temperature_increment_ = 0;
  if (convective_flux_correction_ || diffusive_flux_correction_)
    compute_grad_T_elem_ = 1;
  if (boundary_condition_interface_ == -1 && (!impose_boundary_condition_interface_from_simulation_ && interfacial_boundary_condition_value_!=0.))
    {
      boundary_condition_interface_ = 0; // Dirichlet
      impose_boundary_condition_interface_from_simulation_ = 0;
    }
  if (boundary_condition_end_ == -1 && (!impose_user_boundary_condition_end_value_ && end_boundary_condition_value_!=-1.))
    {
      boundary_condition_end_ = 0; // Dirichlet
      impose_user_boundary_condition_end_value_ = 1;
    }

  return is;
}

void IJK_Thermal_Subresolution::set_param( Param& param )
{
  IJK_Thermal_base::set_param(param);
  param.ajouter_flag("disable_subresolution", &disable_subresolution_);
  param.ajouter_flag("convective_flux_correction", &convective_flux_correction_);
  param.ajouter_flag("diffusive_flux_correction", &diffusive_flux_correction_);
  param.ajouter_flag("disable_fo_flux_correction", &impose_fo_flux_correction_);
  param.ajouter_flag("override_vapour_mixed_values", &override_vapour_mixed_values_);
  param.ajouter("points_per_thermal_subproblem", &points_per_thermal_subproblem_);
  param.ajouter("coeff_distance_diagonal", &coeff_distance_diagonal_);
  param.ajouter("finite_difference_assembler", &finite_difference_assembler_);
  // param.ajouter_flag("disable_mixed_cells_increment", &disable_mixed_cells_increment_);
  param.ajouter_flag("enable_mixed_cells_increment", &enable_mixed_cells_increment_);
  param.ajouter_flag("allow_temperature_correction_for_visu", &allow_temperature_correction_for_visu_);

  param.ajouter("boundary_condition_interface", &boundary_condition_interface_);
  param.dictionnaire("dirichlet", 0);
  param.dictionnaire("neumann", 1);
  param.dictionnaire("flux_jump", 2);
  param.ajouter("interfacial_boundary_condition_value", &interfacial_boundary_condition_value_);
  param.ajouter_flag("impose_boundary_condition_interface_from_simulation", &impose_boundary_condition_interface_from_simulation_);
  param.ajouter("boundary_condition_end", &boundary_condition_end_);
  param.dictionnaire("dirichlet", 0);
  param.dictionnaire("neumann",1);
  param.ajouter("end_boundary_condition_value", &end_boundary_condition_value_);
  param.ajouter_flag("impose_user_boundary_condition_end_value", &impose_user_boundary_condition_end_value_);

  param.ajouter("source_terms_type", &source_terms_type_);
  param.dictionnaire("linear_diffusion", 0);
  param.dictionnaire("spherical_diffusion",1);
  param.dictionnaire("spherical_diffusion_approx",2);
  param.dictionnaire("tangential_conv_2D", 3);
  param.dictionnaire("tangential_conv_3D", 4);
  param.dictionnaire("tangential_conv_2D_tangential_diffusion_3D", 5);
  param.dictionnaire("tangential_conv_3D_tangentual_diffusion_3D", 6);
  param.ajouter_flag("source_terms_correction", &source_terms_correction_);

  // param.ajouter("thermal_fd_solver", &one_dimensional_advection_diffusion_thermal_solver_);
  param.ajouter("thermal_fd_solver", &one_dimensional_advection_diffusion_thermal_solver_);

  param.ajouter_flag("discrete_integral", &discrete_integral_);
  param.ajouter("quadtree_levels", &quadtree_levels_);
  param.ajouter_flag("advected_frame_of_reference", &advected_frame_of_reference_);
  param.ajouter_flag("neglect_frame_of_reference_radial_advection", &neglect_frame_of_reference_radial_advection_);
  param.ajouter_flag("approximate_temperature_increment", &approximate_temperature_increment_);

  param.ajouter_flag("first_time_step_temporal", &first_time_step_temporal_);
  param.ajouter_flag("first_time_step_implicit", &first_time_step_implicit_);
  param.ajouter("local_fourier", &local_fourier_);
  param.ajouter("local_cfl", &local_cfl_);
  param.ajouter("delta_T_subcooled_overheated", &delta_T_subcooled_overheated_);
  param.ajouter_flag("local_diffusion_fourier_priority", &local_diffusion_fourier_priority_);

  param.ajouter_flag("first_time_step_varying_probes", &first_time_step_varying_probes_);
  param.ajouter_flag("probe_variations_priority", &probe_variations_priority_);
  param.ajouter_flag("max_u_radial", &max_u_radial_);

  param.ajouter_flag("distance_cell_faces_from_lrs", &distance_cell_faces_from_lrs_);

  param.ajouter_flag("pre_initialise_thermal_subproblems_list", &pre_initialise_thermal_subproblems_list_);
  param.ajouter("pre_factor_subproblems_number", &pre_factor_subproblems_number_);

  param.ajouter_flag("correct_temperature_cell_neighbours_first_iter", &correct_temperature_cell_neighbours_first_iter_);
  param.ajouter_flag("find_temperature_cell_neighbours", &find_temperature_cell_neighbours_);
  param.ajouter_flag("correct_neighbours_using_probe_length", &correct_neighbours_using_probe_length_);
  param.ajouter("neighbours_corrected_rank", &neighbours_corrected_rank_);

  param.ajouter_flag("neighbours_colinearity_weighting", &neighbours_colinearity_weighting_);
  param.ajouter_flag("use_temperature_cell_neighbours", &use_temperature_cell_neighbours_);

  param.ajouter_flag("clip_temperature_values", &clip_temperature_values_);
  param.ajouter_flag("disable_post_processing_probes_out_files", &disable_post_processing_probes_out_files_);

  param.ajouter_flag("enforce_periodic_boundary_value", &enforce_periodic_boundary_value_);
  param.ajouter_flag("store_cell_faces_corrected", &store_cell_faces_corrected_);

  param.ajouter_flag("find_cell_neighbours_for_fluxes_spherical_correction", &find_cell_neighbours_for_fluxes_spherical_correction_);
  param.ajouter_flag("use_cell_neighbours_for_fluxes_spherical_correction", &use_cell_neighbours_for_fluxes_spherical_correction_);

  param.ajouter_flag("find_reachable_fluxes", &find_reachable_fluxes_);
  param.ajouter_flag("use_reachable_fluxes", &use_reachable_fluxes_);

  //  for (int i=0; i<fd_solvers_jdd_.size(); i++)
  //    param.ajouter_non_std(fd_solvers_jdd_[i], this);
}

int IJK_Thermal_Subresolution::lire_motcle_non_standard(const Motcle& mot, Entree& is)
{
  read_fd_solver(mot, is);
  return 1;
}

void IJK_Thermal_Subresolution::read_fd_solver(const Motcle& mot, Entree& is)
{
  Nom type = "";
  fd_solver_rank_ = fd_solvers_jdd_.search(mot);
  type += fd_solvers_[fd_solver_rank_];
  Cerr << "Finite difference solver: " << type << finl;
  fd_solver_type_ = fd_solvers_[fd_solver_rank_];
  one_dimensional_advection_diffusion_thermal_solver_.typer(type);
  // is >> one_dimensional_advection_diffusion_thermal_solver_;
  // Cerr << "Finite difference solver: " << fd_solver_type_ << finl;
}


int IJK_Thermal_Subresolution::initialize(const IJK_Splitting& splitting, const int idx)
{
  Cout << que_suis_je() << "::initialize()" << finl;

  uniform_lambda_ = lambda_liquid_;
  uniform_alpha_ =	lambda_liquid_ / (ref_ijk_ft_->get_rho_l() * cp_liquid_);
  //  calulate_grad_T_ = 1;
  // TODO: Reused grad T calculation
  calulate_grad_T_=0;
  /*
   * Be careful, it plays a role for allocating the fields in
   * IJK_Thermal_base::initialize
   */
  if (ghost_fluid_)
    {
      compute_grad_T_interface_ = 1;
      compute_curvature_ = 1;
      compute_distance_= 1;
    }

  int nalloc = 0;
  nalloc = IJK_Thermal_base::initialize(splitting, idx);
  corrige_flux_.typer("Corrige_flux_FT_temperature_subresolution");

  temperature_before_extrapolation_.allocate(splitting, IJK_Splitting::ELEM, ghost_cells_);
  nalloc += 1;
  temperature_before_extrapolation_.data() = 0.;

  use_temperature_cell_neighbours_ = (use_temperature_cell_neighbours_ && find_temperature_cell_neighbours_);

  find_cell_neighbours_for_fluxes_spherical_correction_ = find_cell_neighbours_for_fluxes_spherical_correction_ && distance_cell_faces_from_lrs_;
  use_cell_neighbours_for_fluxes_spherical_correction_ = use_cell_neighbours_for_fluxes_spherical_correction_ && distance_cell_faces_from_lrs_;
  if (use_cell_neighbours_for_fluxes_spherical_correction_)
    find_cell_neighbours_for_fluxes_spherical_correction_ = use_cell_neighbours_for_fluxes_spherical_correction_;
  corrige_flux_.set_correction_cell_faces_neighbours(find_cell_neighbours_for_fluxes_spherical_correction_,
                                                     use_cell_neighbours_for_fluxes_spherical_correction_,
                                                     find_reachable_fluxes_);

  temperature_diffusion_op_.set_conductivity_coefficient(uniform_lambda_, temperature_, temperature_, temperature_, temperature_);
  Cerr << "Uniform lambda: " << temperature_diffusion_op_.get_uniform_lambda() << finl;

  impose_fo_flux_correction_ = !impose_fo_flux_correction_;
  convective_flux_correction_ = convective_flux_correction_ && (!conv_temperature_negligible_);
  diffusive_flux_correction_ = diffusive_flux_correction_ && (!diff_temperature_negligible_);

  use_reachable_fluxes_ = (use_reachable_fluxes_ && find_reachable_fluxes_) && (convective_flux_correction_ || diffusive_flux_correction_);
  if (use_reachable_fluxes_)
    {
      find_temperature_cell_neighbours_ = 1;
      use_temperature_cell_neighbours_ = 1;
    }

  disable_mixed_cells_increment_ = (!enable_mixed_cells_increment_);

  corrige_flux_.set_convection_diffusion_correction(convective_flux_correction_, diffusive_flux_correction_);
  corrige_flux_.set_fluxes_feedback_params(discrete_integral_, quadtree_levels_);
  corrige_flux_.set_debug(debug_);
  corrige_flux_.set_distance_cell_faces_from_lrs(distance_cell_faces_from_lrs_);
  corrige_flux_.set_correction_cell_neighbours(find_temperature_cell_neighbours_, neighbours_colinearity_weighting_);
  corrige_flux_.set_eulerian_normal_vectors_ns_normed(eulerian_normal_vectors_ns_normed_);


  if (diffusive_flux_correction_ || use_cell_neighbours_for_fluxes_spherical_correction_)
    {
      // Use an operator that override compute_set() with corrige_flux;
      temperature_diffusion_op_.typer("OpDiffUniformIJKScalarCorrection_double");
      temperature_diffusion_op_.initialize(splitting);
      temperature_diffusion_op_.set_corrige_flux(corrige_flux_);
    }

  if (convective_flux_correction_ || use_cell_neighbours_for_fluxes_spherical_correction_)
    {
      // Use an operator that override compute_set() with corrige_flux;
      temperature_convection_op_.typer("OpConvQuickInterfaceIJKScalar_double");
      temperature_convection_op_.initialize(splitting);
      temperature_convection_op_.set_corrige_flux(corrige_flux_);
    }

  if((convective_flux_correction_ || diffusive_flux_correction_) && store_cell_faces_corrected_)
    {
      allocate_cell_vector(cell_faces_corrected_convective_, splitting, 1);
      allocate_cell_vector(cell_faces_corrected_diffusive_, splitting, 1);
      allocate_cell_vector(cell_faces_corrected_bool_, splitting, 1);
      nalloc += 9;
      for (int c=0; c<3; c++)
        {
          cell_faces_corrected_convective_[c].data() = 0.;
          cell_faces_corrected_diffusive_[c].data() = 0.;
          cell_faces_corrected_bool_[c].data() = 0;
        }
      cell_faces_corrected_convective_.echange_espace_virtuel();
      cell_faces_corrected_diffusive_.echange_espace_virtuel();
      cell_faces_corrected_bool_.echange_espace_virtuel();
    }

  if (approximate_temperature_increment_)
    {
      d_temperature_uncorrected_.allocate(splitting, IJK_Splitting::ELEM, 2);
      div_coeff_grad_T_volume_uncorrected_.allocate(splitting, IJK_Splitting::ELEM, 0);
      nalloc += 2;
      // Centre2
      temperature_diffusion_op_.typer("standard");
      temperature_diffusion_op_.initialize(splitting);
      // Quick
      temperature_convection_op_.typer("standard");
      temperature_convection_op_.initialize(splitting);
    }

  thermal_local_subproblems_.associer(ref_ijk_ft_);
  initialise_thermal_subproblems_list();


  is_first_time_step_ = (!ref_ijk_ft_->get_reprise()) && (ref_ijk_ft_->get_tstep()==0);
  first_time_step_temporal_ = first_time_step_temporal_ && is_first_time_step_;
  first_time_step_explicit_ = !first_time_step_implicit_;
  if (first_time_step_varying_probes_)
    first_time_step_temporal_ = 0;
  probe_variations_enabled_ = first_time_step_varying_probes_;
  compute_overall_probes_parameters();

  if (one_dimensional_advection_diffusion_thermal_solver_.est_nul())
    one_dimensional_advection_diffusion_thermal_solver_.cast_direct_solver_by_default();

  /*
   * Considered constant values
   */
  corrige_flux_.set_physical_parameters(cp_liquid_ * ref_ijk_ft_->get_rho_l(),
                                        cp_vapour_ * ref_ijk_ft_->get_rho_v(),
                                        lambda_liquid_,
                                        lambda_vapour_);
  corrige_flux_.initialize_with_subproblems(
    ref_ijk_ft_->get_splitting_ns(),
    temperature_,
    ref_ijk_ft_->itfce(),
    ref_ijk_ft_,
    ref_ijk_ft_->get_set_interface().get_set_intersection_ijk_face(),
    ref_ijk_ft_->get_set_interface().get_set_intersection_ijk_cell(),
    thermal_local_subproblems_);
  corrige_flux_->set_convection_negligible(!convective_flux_correction_);
  corrige_flux_->set_diffusion_negligible(!diffusive_flux_correction_);

  if (debug_)
    {
      debug_LRS_cells_.allocate(splitting, IJK_Splitting::ELEM, 0);
      nalloc += 1;
      debug_LRS_cells_.data() = -1.;
    }

  if (find_temperature_cell_neighbours_)
    {
      neighbours_temperature_to_correct_.allocate(splitting, IJK_Splitting::ELEM, ghost_cells_);
      temperature_cell_neighbours_.allocate(splitting, IJK_Splitting::ELEM, ghost_cells_);
      nalloc += 2;
      neighbours_temperature_to_correct_.data() = 0;
      temperature_cell_neighbours_.data() = 0.;
      neighbours_temperature_to_correct_.echange_espace_virtuel(neighbours_temperature_to_correct_.ghost());
      temperature_cell_neighbours_.echange_espace_virtuel(temperature_cell_neighbours_.ghost());
      if (neighbours_colinearity_weighting_)
        {
          neighbours_temperature_colinearity_weighting_.allocate(splitting, IJK_Splitting::ELEM, ghost_cells_);
          nalloc += 1;
          neighbours_temperature_colinearity_weighting_.data() = 0.;
          neighbours_temperature_colinearity_weighting_.echange_espace_virtuel(neighbours_temperature_colinearity_weighting_.ghost());
        }
      if (debug_)
        {
          temperature_cell_neighbours_debug_.allocate(splitting, IJK_Splitting::ELEM, ghost_cells_);
          nalloc += 1;
          temperature_cell_neighbours_debug_.data() = 0.;
          temperature_cell_neighbours_debug_.echange_espace_virtuel(temperature_cell_neighbours_debug_.ghost());
        }
    }

  if (find_cell_neighbours_for_fluxes_spherical_correction_)
    {
      allocate_cell_vector(cell_faces_neighbours_corrected_diag_bool_, splitting, ghost_cells_);
      nalloc += 3;
      for (int c=0; c<3; c++)
        cell_faces_neighbours_corrected_diag_bool_[c].data() = 0;
      cell_faces_neighbours_corrected_diag_bool_.echange_espace_virtuel();
      corrige_flux_->set_cell_faces_neighbours_corrected_bool(cell_faces_neighbours_corrected_diag_bool_);
    }


  if (find_reachable_fluxes_)
    {
      allocate_cell_vector(cell_faces_neighbours_corrected_all_bool_, splitting, ghost_cells_);
      nalloc += 3;
      for (int c=0; c<3; c++)
        cell_faces_neighbours_corrected_all_bool_[c].data() = 0;
      cell_faces_neighbours_corrected_all_bool_.echange_espace_virtuel();

      allocate_cell_vector(cell_faces_neighbours_corrected_min_max_bool_, splitting, 1);
      nalloc += 3;
      for (int c=0; c<3; c++)
        cell_faces_neighbours_corrected_min_max_bool_[c].data() = 0;
      cell_faces_neighbours_corrected_min_max_bool_.echange_espace_virtuel();
    }
  if (use_reachable_fluxes_)
    {
      if (convective_flux_correction_)
        {
          allocate_cell_vector(cell_faces_neighbours_corrected_convective_, splitting, 1);
          nalloc += 3;
          for (int c=0; c<3; c++)
            cell_faces_neighbours_corrected_convective_[c].data() = 0.;
          cell_faces_neighbours_corrected_convective_.echange_espace_virtuel();
        }
      if (diffusive_flux_correction_)
        {
          allocate_cell_vector(cell_faces_neighbours_corrected_diffusive_, splitting, 1);
          nalloc += 3;
          for (int c=0; c<3; c++)
            cell_faces_neighbours_corrected_diffusive_[c].data() = 0.;
          cell_faces_neighbours_corrected_diffusive_.echange_espace_virtuel();
        }
      if (neighbours_colinearity_weighting_)
        {
          allocate_cell_vector(neighbours_faces_weighting_colinearity_, splitting, 1);
          nalloc += 3;
          for (int c=0; c<3; c++)
            neighbours_faces_weighting_colinearity_[c].data() = 0.;
          neighbours_faces_weighting_colinearity_.echange_espace_virtuel();
        }
    }

  // double ideal_length_factor = ref_ijk_ft_->get_remaillage_ft_ijk().get_facteur_longueur_ideale();
  if (diffusive_flux_correction_ && fo_ >= 1.)
    fo_ = 0.95 * (sqrt(2) / 2);

  Cout << "End of " << que_suis_je() << "::initialize()" << finl;
  return nalloc;
}

void IJK_Thermal_Subresolution::update_thermal_properties()
{
  IJK_Thermal_base::update_thermal_properties();
}

void IJK_Thermal_Subresolution::compute_diffusion_increment()
{
  /*
   * Update d_temperature
   * d_temperature_ += div_lambda_grad_T_volume_;
   * It depends on the nature of the properties (one-fluid or single-fluid)
   */
  const int ni = div_coeff_grad_T_volume_.ni();
  const int nj = div_coeff_grad_T_volume_.nj();
  const int nk = div_coeff_grad_T_volume_.nk();
  const double rhocp_l = ref_ijk_ft_->get_rho_l() * cp_liquid_;
  for (int k = 0; k < nk; k++)
    for (int j = 0; j < nj; j++)
      for (int i = 0; i < ni; i++)
        {
          const double rhocpVol = rhocp_l * vol_;
          const double ope = div_coeff_grad_T_volume_(i,j,k);
          const double resu = ope / rhocpVol;
          d_temperature_(i,j,k) += resu;
        }
  Cerr << "Uniform lambda: " << temperature_diffusion_op_.get_uniform_lambda() << finl;
}

void IJK_Thermal_Subresolution::correct_temperature_for_eulerian_fluxes()
{
  if (override_vapour_mixed_values_)
    {
      const int ni = temperature_.ni();
      const int nj = temperature_.nj();
      const int nk = temperature_.nk();
      for (int k = 0; k < nk; k++)
        for (int j = 0; j < nj; j++)
          for (int i = 0; i < ni; i++)
            {
              const double indic = ref_ijk_ft_->itfce().I(i,j,k);
              if (fabs(indic) < LIQUID_INDICATOR_TEST) // Mixed cells and pure vapour cells
                temperature_(i,j,k) = 0.;
            }
    }
}

void IJK_Thermal_Subresolution::store_temperature_before_extrapolation()
{
  temperature_before_extrapolation_.data() = 0.;
  temperature_before_extrapolation_.echange_espace_virtuel(temperature_before_extrapolation_.ghost());
  const int ni = temperature_.ni();
  const int nj = temperature_.nj();
  const int nk = temperature_.nk();
  for (int k = 0; k < nk; k++)
    for (int j = 0; j < nj; j++)
      for (int i = 0; i < ni; i++)
        {
          const double indic = ref_ijk_ft_->itfce().I(i,j,k);
          if (fabs(indic) > VAPOUR_INDICATOR_TEST)
            {
              if (ref_ijk_ft_->get_current_time() == 0 && !ref_ijk_ft_->get_reprise())
                temperature_before_extrapolation_(i,j,k) = delta_T_subcooled_overheated_;
              else
                {
                  const double temperature = temperature_(i,j,k);
                  temperature_before_extrapolation_(i,j,k) = temperature;
                }
            }
          else
            temperature_before_extrapolation_(i,j,k) = 0.;
        }
  temperature_before_extrapolation_.echange_espace_virtuel(temperature_before_extrapolation_.ghost());
}

void IJK_Thermal_Subresolution::correct_temperature_increment_for_interface_leaving_cell()
{
  /*
   * Correct only if we have not extended the temperature field across the interface (no fluxes calculation)
   */
  const int ni = d_temperature_.ni();
  const int nj = d_temperature_.nj();
  const int nk = d_temperature_.nk();
  if (disable_mixed_cells_increment_ && disable_subresolution_ && ghost_fluid_)
    {
      for (int k = 0; k < nk; k++)
        for (int j = 0; j < nj; j++)
          for (int i = 0; i < ni; i++)
            {
              const double indic = ref_ijk_ft_->itfce().I(i,j,k);
              if (fabs(indic)<LIQUID_INDICATOR_TEST) // Mixed cells and pure vapour cells
                { d_temperature_(i,j,k) = 0; }
            }
    }
}

void IJK_Thermal_Subresolution::correct_temperature_for_visu()
{
  /*
   * Correct only if the temperature gradient is post-processed
   * If not we may want to reconstruct the gradient field in Python
   * using the ghost temperature !
   */
  if (liste_post_instantanes_.contient_("GRAD_T_ELEM") && allow_temperature_correction_for_visu_)
    {
      const int ni = temperature_.ni();
      const int nj = temperature_.nj();
      const int nk = temperature_.nk();
      for (int k = 0; k < nk; k++)
        for (int j = 0; j < nj; j++)
          for (int i = 0; i < ni; i++)
            {
              // const double temperature = temperature_(i,j,k);
              const double indic = ref_ijk_ft_->itfce().I(i,j,k);
              // if (temperature > 0)
              if (indic < VAPOUR_INDICATOR_TEST)
                temperature_(i,j,k) = 0;
            }
      temperature_.echange_espace_virtuel(temperature_.ghost());
    }
}

void IJK_Thermal_Subresolution::clip_temperature_values()
{
  if (clip_temperature_values_)
    {
      const int ni = temperature_.ni();
      const int nj = temperature_.nj();
      const int nk = temperature_.nk();
      for (int k = 0; k < nk; k++)
        for (int j = 0; j < nj; j++)
          for (int i = 0; i < ni; i++)
            {
              const double temperature = temperature_(i,j,k);
              if (temperature < delta_T_subcooled_overheated_)
                temperature_(i,j,k) = delta_T_subcooled_overheated_;
            }
      temperature_.echange_espace_virtuel(temperature_.ghost());
    }
}

void IJK_Thermal_Subresolution::compute_thermal_subproblems()
{
  is_first_time_step_ = (!ref_ijk_ft_->get_reprise()) && (ref_ijk_ft_->get_tstep()==0);
  if (is_first_time_step_)
    first_time_step_temporal_ = first_time_step_temporal_ && is_first_time_step_;
  /*
   * FIXME: Do the first step only at the first iteration
   */
  if (debug_)
    debug_LRS_cells_.data() = -1.;

  if (debug_)
    Cerr << "Initialise thermal subproblems" << finl;
  initialise_thermal_subproblems();

  pre_initialise_thermal_subproblems_matrix();

  reset_subresolution_distributed_vectors();

  int varying_probes_length = 1;
  const int varying_probes = first_time_step_varying_probes_;
  do
    {
      if (debug_)
        Cerr << "Interpolate the velocities on probes and compute local cfl fourier conditions" << finl;
      interpolate_project_velocities_on_probes();
      if (varying_probes_length)
        reajust_probes_length();
      else
        probe_variations_enabled_ = 0;
      varying_probes_length--;
    }
  while (probe_variations_enabled_);

  if (debug_)
    if(varying_probes && !first_time_step_varying_probes_)
      Cerr << "Probes length is now fixed" << finl;

  if (debug_)
    Cerr << "Compute radial subresolution convection and diffusion operators" << finl;
  compute_radial_subresolution_convection_diffusion_operators();

  if (debug_)
    Cerr << "Compute local substep for the first iter" << finl;
  compute_local_substep();
  prepare_temporal_schemes();

  if (debug_)
    Cerr << "Prepare boundary conditions, compute source terms and impose boundary conditions" << finl;
  temperature_.echange_espace_virtuel(temperature_.ghost());
  compute_source_terms_impose_subresolution_boundary_conditions();

  if (debug_)
    Cerr << "Compute material derivative (modelling)" << finl;
  approximate_temperature_increment_material_derivative();

  if (debug_)
    Cerr << "Solve thermal subproblems" << finl;
  solve_thermal_subproblems();

  if (debug_)
    Cerr << "Prepare thermal flux correction" << finl;
  update_intersections();
  prepare_thermal_flux_correction();

  // apply_thermal_flux_correction();
}

void IJK_Thermal_Subresolution::compute_ghost_cell_numbers_for_subproblems(const IJK_Splitting& splitting, int ghost_init)
{
  int ghost_cells;
  compute_cell_diagonal(splitting);
  const IJK_Grid_Geometry& geom = splitting.get_grid_geometry();
  const double dx = geom.get_constant_delta(DIRECTION_I);
  const double dy = geom.get_constant_delta(DIRECTION_J);
  const double dz = geom.get_constant_delta(DIRECTION_K);
  double maximum_distance = cell_diagonal_ * coeff_distance_diagonal_;
  const int max_dx = (int) trunc(maximum_distance / dx + 1);
  const int max_dy = (int) trunc(maximum_distance / dy + 1);
  const int max_dz = (int) trunc(maximum_distance / dz + 1);
  ghost_cells = std::max(max_dx, std::max(max_dy, max_dz));
  // Add one or more cell ??
  if (ghost_cells < ghost_init)
    ghost_cells = ghost_init;
  ghost_cells_ = ghost_cells;
}

void IJK_Thermal_Subresolution::compute_overall_probes_parameters()
{
  if (!disable_subresolution_)
    {
      probe_length_ = coeff_distance_diagonal_ * cell_diagonal_;
      dr_ = probe_length_ / (points_per_thermal_subproblem_ - 1);
      radial_coordinates_ = DoubleVect(points_per_thermal_subproblem_);
      for (int i=0; i < points_per_thermal_subproblem_; i++)
        radial_coordinates_(i) = i * dr_;

      /*
       * Compute the matrices for Finite-Differences (first iteration only)
       */
      compute_radial_first_second_order_operators(radial_first_order_operator_raw_, radial_second_order_operator_raw_,
                                                  radial_first_order_operator_, radial_second_order_operator_);
      if (first_time_step_temporal_)
        {
          int check_nb_elem;
          check_nb_elem = finite_difference_assembler_.build(identity_matrix_explicit_implicit_, points_per_thermal_subproblem_, -1);
          Cerr << "Check_nb_elem: " << check_nb_elem << finl;
        }

      thermal_subproblems_matrix_assembly_for_solver_.typer("Matrice_Morse");
      thermal_subproblems_matrix_assembly_for_solver_reduced_.typer("Matrice_Morse");
    }
}

void IJK_Thermal_Subresolution::compute_radial_first_second_order_operators(Matrice& radial_first_order_operator_raw,
                                                                            Matrice& radial_second_order_operator_raw,
                                                                            Matrice& radial_first_order_operator,
                                                                            Matrice& radial_second_order_operator)
{
  /*
   * Compute the matrices for Finite-Differences
   */
  compute_first_order_operator_raw(radial_first_order_operator_raw);
  compute_second_order_operator_raw(radial_second_order_operator_raw);
  radial_first_order_operator = Matrice(radial_first_order_operator_raw);
  radial_second_order_operator = Matrice(radial_second_order_operator_raw);
  compute_first_order_operator(radial_first_order_operator, dr_);
  compute_second_order_operator(radial_second_order_operator, dr_);
}

/*
 * FIXME : Should be moved to a class of tools
 */

void IJK_Thermal_Subresolution::compute_first_order_operator_raw(Matrice& radial_first_order_operator)
{
  /*
   * Compute the first-order matrix for Finite-Differences
   */
  // TODO:
  int check_nb_elem;
  check_nb_elem = finite_difference_assembler_.build(radial_first_order_operator, points_per_thermal_subproblem_, 0);
  Cerr << "Check_nb_elem: " << check_nb_elem << finl;

}

void IJK_Thermal_Subresolution::compute_first_order_operator(Matrice& radial_first_order_operator, double dr)
{
  const double dr_inv = 1 / dr;
  radial_first_order_operator *= dr_inv;
}

void IJK_Thermal_Subresolution::compute_second_order_operator_raw(Matrice& radial_second_order_operator)
{
  /*
   * Compute the second-order matrix for Finite-Differences
   */
  // TODO:
  int check_nb_elem;
  check_nb_elem = finite_difference_assembler_.build(radial_second_order_operator, points_per_thermal_subproblem_, 1);
  Cerr << "Check_nb_elem: " << check_nb_elem << finl;
}

void IJK_Thermal_Subresolution::compute_second_order_operator(Matrice& radial_second_order_operator, double dr)
{
  const double dr_squared_inv = 1 / pow(dr, 2);
  radial_second_order_operator *= dr_squared_inv;
}

void IJK_Thermal_Subresolution::initialise_thermal_subproblems_list()
{
  thermal_local_subproblems_.initialise_thermal_subproblems_list_params(pre_initialise_thermal_subproblems_list_,
                                                                        pre_factor_subproblems_number_);
}

void IJK_Thermal_Subresolution::initialise_thermal_subproblems()
{
  // FIXME : Should I use IJK_Field_local_double
  //	IJK_Field_local_double eulerian_normal_vectors_ns;
  //	IJK_Field_local_double eulerian_facets_barycentre_ns;
  //	const int ni = eulerian_normal_vectors_[0].ni;
  //	const int nj = eulerian_normal_vectors_[0].nj;
  //	const int nk = eulerian_normal_vectors_[0].nk;
  //	const int ng= eulerian_normal_vectors_[0].ghost();
  //	eulerian_normal_vectors_ns.allocate(ni, nj, nk, ng);
  //	eulerian_facets_barycentre_ns.allocate(ni, nj, nk, ng);
  // FIXME : Work with FT fields
  //	compute_mixed_cells_number(ref_ijk_ft_->itfce().I());
  //  thermal_local_subproblems_.add_subproblems(mixed_cells_number_);
  //  FixedVector<IJK_Field_double, 3> eulerian_facets_barycentre_ns;
  //  FixedVector<IJK_Field_double, 3> eulerian_normal_vectors_ns;
  //  IJK_Splitting splitting = temperature_.get_splitting();
  //  allocate_cell_vector(eulerian_facets_barycentre_ns, splitting, 0);
  //  allocate_cell_vector(eulerian_normal_vectors_ns, splitting, 0);
  //  eulerian_facets_barycentre_ns.echange_espace_virtuel();
  //  eulerian_normal_vectors_ns.echange_espace_virtuel();
  //  for (int dir=0; dir < 3; dir++)
  //    {
  //      ref_ijk_ft_->redistribute_from_splitting_ft_elem(eulerian_facets_barycentre_[dir], eulerian_facets_barycentre_ns[dir]);
  //      ref_ijk_ft_->redistribute_from_splitting_ft_elem(eulerian_normal_vectors_[dir], eulerian_normal_vectors_ns[dir]);
  //    }
  // FIXME: should I perform this on the extended fields ? or original ?
  // If so, I need to convert, distance, curvature, interfacial_area earlier !
  // Or I need to calculate gradient and hessian fields using the temperature_ft_ mesh !
  if (!disable_subresolution_)
    {
      const IJK_Field_double& indicator = ref_ijk_ft_->itfce().I();
      const int ni = temperature_.ni();
      const int nj = temperature_.nj();
      const int nk = temperature_.nk();
      for (int k = 0; k < nk; k++)
        for (int j = 0; j < nj; j++)
          for (int i = 0; i < ni; i++)
            if (fabs(indicator(i,j,k)) > VAPOUR_INDICATOR_TEST && fabs(indicator(i,j,k)) < LIQUID_INDICATOR_TEST)
              {
                if (debug_)
                  {
                    Cerr << "Liquid Indicator (Old): " << indicator(i,j,k) << finl;
                    Cerr << "Liquid Indicator (Next): " << ref_ijk_ft_->itfce().In()(i,j,k) << finl;
                    for (int i_bulle=0; i_bulle < bounding_box_.dimension(0); i_bulle++)
                      {
                        Cerr << "Bounding box : " << i_bulle << "; " << bounding_box_(i_bulle,0,0) << finl;
                        Cerr << "Bounding box : " << i_bulle << "; " << bounding_box_(i_bulle,0,1) << finl;
                        Cerr << "Bounding box : " << i_bulle << "; " << bounding_box_(i_bulle,1,0) << finl;
                        Cerr << "Bounding box : " << i_bulle << "; " << bounding_box_(i_bulle,1,1) << finl;
                        Cerr << "Bounding box : " << i_bulle << "; " << bounding_box_(i_bulle,2,0) << finl;
                        Cerr << "Bounding box : " << i_bulle << "; " << bounding_box_(i_bulle,2,1) << finl;

                        Cerr << "Min-Max Bounding box : " << i_bulle << "; " << min_max_larger_box_(i_bulle,0,0) << finl;
                        Cerr << "Min-Max Bounding box : " << i_bulle << "; " << min_max_larger_box_(i_bulle,0,1) << finl;
                        Cerr << "Min-Max Bounding box : " << i_bulle << "; " << min_max_larger_box_(i_bulle,1,0) << finl;
                        Cerr << "Min-Max Bounding box : " << i_bulle << "; " << min_max_larger_box_(i_bulle,1,1) << finl;
                        Cerr << "Min-Max Bounding box : " << i_bulle << "; " << min_max_larger_box_(i_bulle,2,0) << finl;
                        Cerr << "Min-Max Bounding box : " << i_bulle << "; " << min_max_larger_box_(i_bulle,2,1) << finl;
                      }
                    debug_LRS_cells_(i,j,k) = indicator(i,j,k);
                  }
                thermal_local_subproblems_.associate_sub_problem_to_inputs(debug_,
                                                                           i, j, k,
                                                                           ref_ijk_ft_->get_timestep(),
                                                                           ref_ijk_ft_->get_current_time(),
                                                                           *eulerian_compo_connex_ns_,
                                                                           eulerian_distance_ns_,
                                                                           eulerian_curvature_ns_,
                                                                           eulerian_interfacial_area_ns_,
                                                                           eulerian_facets_barycentre_ns_,
                                                                           eulerian_normal_vectors_ns_,
                                                                           rising_velocities_,
                                                                           rising_vectors_,
                                                                           bubbles_barycentre_,
                                                                           advected_frame_of_reference_,
                                                                           neglect_frame_of_reference_radial_advection_,
                                                                           points_per_thermal_subproblem_,
                                                                           uniform_alpha_,
                                                                           uniform_lambda_,
                                                                           coeff_distance_diagonal_,
                                                                           cell_diagonal_,
                                                                           dr_,
                                                                           radial_coordinates_,
                                                                           identity_matrix_explicit_implicit_,
                                                                           radial_first_order_operator_raw_,
                                                                           radial_second_order_operator_raw_,
                                                                           radial_first_order_operator_,
                                                                           radial_second_order_operator_,
                                                                           identity_matrix_subproblems_,
                                                                           radial_diffusion_matrix_,
                                                                           radial_convection_matrix_,
                                                                           ref_ijk_ft_->itfce(),
                                                                           indicator(i,j,k),
                                                                           temperature_,
                                                                           temperature_ft_,
                                                                           temperature_before_extrapolation_,
                                                                           ref_ijk_ft_->get_velocity(),
                                                                           ref_ijk_ft_->get_velocity_ft(),
                                                                           ref_ijk_ft_->get_pressure_ghost_cells(),
                                                                           grad_T_elem_,
                                                                           hess_diag_T_elem_,
                                                                           hess_cross_T_elem_,
                                                                           finite_difference_assembler_,
                                                                           thermal_subproblems_matrix_assembly_,
                                                                           thermal_subproblems_rhs_assembly_,
                                                                           thermal_subproblems_temperature_solution_ini_,
                                                                           thermal_subproblems_temperature_solution_,
                                                                           source_terms_type_,
                                                                           source_terms_correction_,
                                                                           is_first_time_step_,
                                                                           first_time_step_temporal_,
                                                                           first_time_step_explicit_,
                                                                           local_fourier_,
                                                                           local_cfl_,
                                                                           min_delta_xyz_,
                                                                           delta_T_subcooled_overheated_,
                                                                           first_time_step_varying_probes_,
                                                                           probe_variations_priority_,
                                                                           disable_interpolation_in_mixed_cells_,
                                                                           max_u_radial_,
                                                                           (convective_flux_correction_ || diffusive_flux_correction_),
                                                                           distance_cell_faces_from_lrs_,
                                                                           pre_initialise_thermal_subproblems_list_,
                                                                           find_temperature_cell_neighbours_,
                                                                           !correct_neighbours_using_probe_length_,
                                                                           neighbours_corrected_rank_,
                                                                           neighbours_colinearity_weighting_,
                                                                           find_reachable_fluxes_,
                                                                           find_cell_neighbours_for_fluxes_spherical_correction_,
                                                                           n_iter_distance_);

              }
    }
}

void IJK_Thermal_Subresolution::pre_initialise_thermal_subproblems_matrix()
{
  if (ref_ijk_ft_->get_tstep()==0 && pre_initialise_thermal_subproblems_list_)
    {
      if (!disable_subresolution_)
        {
          const int nb_subproblems_ini = thermal_local_subproblems_.get_subproblems_counter();
          const int max_subproblems_predicted = (int) ((double) nb_subproblems_ini * pre_factor_subproblems_number_);
          finite_difference_assembler_.pre_initialise_matrix_subproblems(thermal_subproblems_matrix_assembly_,
                                                                         radial_second_order_operator_raw_,
                                                                         max_subproblems_predicted);
          finite_difference_assembler_.pre_initialise_matrix_subproblems(radial_convection_matrix_,
                                                                         radial_first_order_operator_raw_,
                                                                         max_subproblems_predicted);
          finite_difference_assembler_.pre_initialise_matrix_subproblems(radial_diffusion_matrix_,
                                                                         radial_second_order_operator_raw_,
                                                                         max_subproblems_predicted);
          if (first_time_step_temporal_)
            finite_difference_assembler_.pre_initialise_matrix_subproblems(identity_matrix_subproblems_,
                                                                           radial_second_order_operator_raw_,
                                                                           max_subproblems_predicted);
          const int vector_size = points_per_thermal_subproblem_ * max_subproblems_predicted;
          if (first_time_step_temporal_ && first_time_step_explicit_)
            thermal_subproblems_temperature_solution_ini_.resize(vector_size);
          thermal_subproblems_rhs_assembly_.resize(vector_size);
          thermal_subproblems_temperature_solution_.resize(vector_size);

          finite_difference_assembler_.complete_empty_matrices_initialisation(thermal_subproblems_matrix_assembly_,
                                                                              radial_second_order_operator_raw_,
                                                                              nb_subproblems_ini, max_subproblems_predicted);
          finite_difference_assembler_.complete_empty_matrices_initialisation(radial_convection_matrix_,
                                                                              radial_first_order_operator_raw_,
                                                                              nb_subproblems_ini, max_subproblems_predicted);
          finite_difference_assembler_.complete_empty_matrices_initialisation(radial_diffusion_matrix_,
                                                                              radial_second_order_operator_raw_,
                                                                              nb_subproblems_ini, max_subproblems_predicted);
          if (first_time_step_temporal_)
            finite_difference_assembler_.complete_empty_matrices_initialisation(radial_diffusion_matrix_,
                                                                                radial_second_order_operator_raw_,
                                                                                nb_subproblems_ini, max_subproblems_predicted);
        }
    }
}

void IJK_Thermal_Subresolution::reset_subresolution_distributed_vectors()
{
  if (!global_probes_characteristics_ || !pre_initialise_thermal_subproblems_list_)
    {
      thermal_subproblems_rhs_assembly_.reset();
      thermal_subproblems_rhs_assembly_.set_smart_resize(1);
    }

  if (!pre_initialise_thermal_subproblems_list_)
    {
      thermal_subproblems_temperature_solution_.reset();
      thermal_subproblems_temperature_solution_.set_smart_resize(1);
    }
  else
    {
      if (first_time_step_temporal_ && first_time_step_explicit_)
        {
          thermal_subproblems_temperature_solution_ini_for_solver_.reset();
          thermal_subproblems_temperature_solution_ini_for_solver_.set_smart_resize(1);
        }
    }
  thermal_subproblems_rhs_assembly_for_solver_.reset();
  thermal_subproblems_rhs_assembly_for_solver_.set_smart_resize(1);
  thermal_subproblems_temperature_solution_for_solver_.reset();
  thermal_subproblems_temperature_solution_for_solver_.set_smart_resize(1);
}


void IJK_Thermal_Subresolution::interpolate_project_velocities_on_probes()
{
  if (!disable_subresolution_)
    thermal_local_subproblems_.interpolate_project_velocities_on_probes();
}

void IJK_Thermal_Subresolution::reajust_probes_length()
{
  if (!disable_subresolution_)
    if (first_time_step_varying_probes_)
      {
        thermal_local_subproblems_.reajust_probes_length();
        //
        probe_variations_enabled_ = thermal_local_subproblems_.get_probe_variations_enabled(probe_variations_priority_);
        first_time_step_varying_probes_ = probe_variations_enabled_;
        thermal_local_subproblems_.compute_modified_probe_length(probe_variations_enabled_);
      }
}

void IJK_Thermal_Subresolution::compute_radial_subresolution_convection_diffusion_operators()
{
  if (!disable_subresolution_)
    {
      /*
       * Same spatial discretisation on all probes
       */
      if (!pre_initialise_thermal_subproblems_list_)
        {
          if (first_time_step_temporal_)
            initialise_identity_matrices(identity_matrix_explicit_implicit_, identity_matrix_subproblems_);
          initialise_radial_convection_operator(radial_first_order_operator_, radial_convection_matrix_);
          initialise_radial_diffusion_operator(radial_second_order_operator_, radial_diffusion_matrix_);
        }
      thermal_local_subproblems_.compute_radial_convection_diffusion_operators();
    }
}

void IJK_Thermal_Subresolution::initialise_identity_matrices(Matrice& identity_matrix,
                                                             Matrice& identity_matrix_subproblems)
{
  /*
   * Compute the convection matrices for Finite-Differences
   */
  if (debug_)
    Cerr << "Initialise the Identity matrices" << finl;
  const int nb_subproblems = thermal_local_subproblems_.get_subproblems_counter();
  finite_difference_assembler_.initialise_matrix_subproblems(identity_matrix_subproblems,
                                                             identity_matrix,
                                                             nb_subproblems,
                                                             first_time_step_varying_probes_);
  if (debug_)
    Cerr << "Radial convection operator has been initialised." << finl;
}

void IJK_Thermal_Subresolution::initialise_radial_convection_operator(Matrice& radial_first_order_operator,
                                                                      Matrice& radial_convection_matrix)
{
  /*
   * Compute the convection matrices for Finite-Differences
   */
  if (debug_)
    Cerr << "Initialise the Radial convection operator" << finl;
  const int nb_subproblems = thermal_local_subproblems_.get_subproblems_counter();
  finite_difference_assembler_.initialise_matrix_subproblems(radial_convection_matrix, radial_first_order_operator, nb_subproblems, first_time_step_varying_probes_);
  if (debug_)
    Cerr << "Radial convection operator has been initialised." << finl;
}

void IJK_Thermal_Subresolution::initialise_radial_diffusion_operator(Matrice& radial_second_order_operator,
                                                                     Matrice& radial_diffusion_matrix)
{
  /*
   * Compute the diffusion matrices for Finite-Differences
   */
  if (debug_)
    Cerr << "Initialise the Radial diffusion operator" << finl;
  const int nb_subproblems = thermal_local_subproblems_.get_subproblems_counter();
  finite_difference_assembler_.initialise_matrix_subproblems(radial_diffusion_matrix, radial_second_order_operator, nb_subproblems, first_time_step_varying_probes_);
  if (debug_)
    Cerr << "Radial diffusion operator has been initialised." << finl;
}

void IJK_Thermal_Subresolution::compute_local_substep()
{
  if (!disable_subresolution_)
    if (first_time_step_temporal_)
      {
        const double local_time_step = thermal_local_subproblems_.get_min_euler_time_step(nb_iter_explicit_local_);
        Cerr << "Local timestep: " << local_time_step << finl;
        Cerr << "Number of sub-steps: " << nb_iter_explicit_local_ << finl;
        thermal_local_subproblems_.set_local_time_step(local_time_step);
        local_fourier_time_step_probe_length_ = thermal_local_subproblems_.get_local_min_fourier_time_step_probe_length();
        local_cfl_time_step_probe_length_ = thermal_local_subproblems_.get_local_min_cfl_time_step_probe_length();
        local_dt_cfl_ = thermal_local_subproblems_.get_local_dt_cfl();
        local_dt_cfl_min_delta_xyz_ = thermal_local_subproblems_.get_local_dt_cfl_min_delta_xyz();
        local_dt_cfl_counter_ += (ref_ijk_ft_->get_timestep() / local_dt_cfl_min_delta_xyz_);
        local_dt_fourier_counter_ += (ref_ijk_ft_->get_timestep() / local_fourier_time_step_probe_length_);
        if (debug_)
          {
            Cerr << "Compare the thermal time-steps" << finl;
            Cerr << "Current time: "<< ref_ijk_ft_->get_current_time() << finl;
            Cerr << "local_fourier_time_step_probe_length: " << local_fourier_time_step_probe_length_ << finl;
            Cerr << "local_cfl_time_step_probe_length: " << local_cfl_time_step_probe_length_ << finl;
            Cerr << "local_dt_cfl: " << local_dt_cfl_ << finl;
            Cerr << "local_dt_cfl_min_delta_xyz: " << local_dt_cfl_min_delta_xyz_ << finl;
            Cerr << "local_dt_cfl_counter: " << local_dt_cfl_counter_ << finl;
            Cerr << "local_dt_fourier_counter: " << local_dt_fourier_counter_ << finl;
            Cerr << "dt_cfl: " << ref_ijk_ft_->get_dt_cfl_liq() << finl;
          }
        /*
         * Activate Quasi-static when convection is significant
         */
        if (local_dt_cfl_counter_ > 1. && !local_diffusion_fourier_priority_)
          first_time_step_temporal_ = 0;
        else if (local_dt_fourier_counter_ > 1.)
          first_time_step_temporal_ = 0;

        // if (local_fourier_time_step_probe_length_ < local_dt_cfl_ * convection_diffusion_time_scale_factor_)
        // convection_diffusion_time_scale_factor_ = ref_ijk_ft_->get_dt_cfl_liq() / local_cfl_time_step_probe_length_;
//        if (local_fourier_time_step_probe_length_ < (convection_diffusion_time_scale_factor_ * local_cfl_time_step_probe_length_))
//          {
//            //TODO: Check that it works (or ref_ijk_ft_->get_timestep() * ref_ijk_ft_->get_tstep())
//            // double current_time = ref_ijk_ft_->get_current_time();
//            // if ((current_time > local_fourier_time_step_probe_length_) || (local_dt_cfl_counter_ > 1.))
//        		if (local_dt_fourier_counter_ > 1. || (local_dt_cfl_counter_ > 1. && !local_diffusion_fourier_priority_))
//              first_time_step_temporal_ = 0;
//          }
//        else
//        	if (local_dt_cfl_counter_ > 1.)
//        		first_time_step_temporal_ = 0;
      }
}

void IJK_Thermal_Subresolution::prepare_temporal_schemes()
{
  if (!disable_subresolution_)
    {
      if (first_time_step_temporal_)
        thermal_local_subproblems_.prepare_temporal_schemes();
      thermal_subproblems_matrix_assembly_ = radial_convection_matrix_;
      finite_difference_assembler_.sum_matrices_subproblems(thermal_subproblems_matrix_assembly_, radial_diffusion_matrix_);
      if (first_time_step_temporal_)
        finite_difference_assembler_.sum_matrices_subproblems(thermal_subproblems_matrix_assembly_, identity_matrix_subproblems_);
    }
}

void IJK_Thermal_Subresolution::compute_source_terms_impose_subresolution_boundary_conditions()
{
  if (!disable_subresolution_)
    {
      thermal_local_subproblems_.compute_source_terms_impose_boundary_conditions(thermal_subproblems_rhs_assembly_,
                                                                                 thermal_subproblems_temperature_solution_ini_,
                                                                                 boundary_condition_interface_,
                                                                                 interfacial_boundary_condition_value_,
                                                                                 impose_boundary_condition_interface_from_simulation_,
                                                                                 boundary_condition_end_,
                                                                                 end_boundary_condition_value_,
                                                                                 impose_user_boundary_condition_end_value_);
    }
}

void IJK_Thermal_Subresolution::approximate_temperature_increment_material_derivative()
{
  if (!disable_subresolution_)
    {
      thermal_local_subproblems_.approximate_temperature_increment_material_derivative();
    }
}

/*
 * TODO: DEBUG
 */
void IJK_Thermal_Subresolution::solve_thermal_subproblems()
{
  if (!disable_subresolution_)
    {
      int nb_points = copy_local_unknwowns_rhs();
      check_wrong_values_rhs();
      convert_into_sparse_matrix(nb_points);

      if (first_time_step_temporal_ && first_time_step_explicit_)
        {
          Cerr << "Finite-difference Explicit thermal sub-resolution !" << finl;
          DoubleVect thermal_subproblems_temperature_solution_tmp = thermal_subproblems_temperature_solution_ini_;
          for (int i=0; i<nb_iter_explicit_local_; i++)
            {
              thermal_subproblems_temperature_solution_ = thermal_subproblems_matrix_assembly_for_solver_reduced_ * thermal_subproblems_temperature_solution_tmp;
              thermal_subproblems_temperature_solution_ += thermal_subproblems_rhs_assembly_;
              thermal_subproblems_temperature_solution_tmp = thermal_subproblems_temperature_solution_;
            }
          Cerr << "Finite-difference Explicit thermal sub-resolution has finished in "<< nb_iter_explicit_local_ << "iterations!" << finl;
        }
      else
        {

          compute_md_vector();
          if (first_time_step_temporal_ && !first_time_step_explicit_)
            {
              if (one_dimensional_advection_diffusion_thermal_solver_implicit_.est_nul())
                one_dimensional_advection_diffusion_thermal_solver_implicit_.cast_direct_solver_by_default();
              Cerr << "Finite-difference thermal sub-resolution has started (Implicit)!" << finl;
              one_dimensional_advection_diffusion_thermal_solver_implicit_.resoudre_systeme(thermal_subproblems_matrix_assembly_for_solver_reduced_.valeur(),
                                                                                            thermal_subproblems_rhs_assembly_for_solver_,
                                                                                            thermal_subproblems_temperature_solution_for_solver_);
              Cerr << "Finite-difference thermal sub-resolution has finished (Implicit)!" << finl;
            }
          else
            {
              Cerr << "Finite-difference thermal sub-resolution has started !" << finl;
              one_dimensional_advection_diffusion_thermal_solver_.resoudre_systeme(thermal_subproblems_matrix_assembly_for_solver_reduced_.valeur(),
                                                                                   thermal_subproblems_rhs_assembly_for_solver_,
                                                                                   thermal_subproblems_temperature_solution_for_solver_);
              Cerr << "Finite-difference thermal sub-resolution has finished !" << finl;
            }
        }
      retrieve_temperature_solution();
    }
}


int IJK_Thermal_Subresolution::copy_local_unknwowns_rhs()
{
  const int nb_subproblems_ini = thermal_local_subproblems_.get_subproblems_counter();
  const int nb_points = points_per_thermal_subproblem_ * nb_subproblems_ini;
  if (!pre_initialise_thermal_subproblems_list_)
    {
      thermal_subproblems_temperature_solution_.resize(thermal_subproblems_rhs_assembly_.size());
      if (first_time_step_temporal_ && first_time_step_explicit_)
        thermal_subproblems_temperature_solution_ini_for_solver_ = thermal_subproblems_temperature_solution_ini_;
      thermal_subproblems_rhs_assembly_for_solver_ = thermal_subproblems_rhs_assembly_;
      thermal_subproblems_temperature_solution_for_solver_ = thermal_subproblems_temperature_solution_;
    }
  else
    {
      if (first_time_step_temporal_ && first_time_step_explicit_)
        thermal_subproblems_temperature_solution_ini_for_solver_.resize(nb_points);
      thermal_subproblems_rhs_assembly_for_solver_.resize(nb_points);
      thermal_subproblems_temperature_solution_for_solver_.resize(nb_points);
      for (int i=0; i<nb_points; i++)
        {
          if (first_time_step_temporal_ && first_time_step_explicit_)
            thermal_subproblems_temperature_solution_ini_for_solver_(i) = thermal_subproblems_temperature_solution_ini_(i);
          thermal_subproblems_rhs_assembly_for_solver_(i) = thermal_subproblems_rhs_assembly_(i);
          thermal_subproblems_temperature_solution_for_solver_(i) = thermal_subproblems_temperature_solution_(i);
        }
    }
  return nb_points;
}


void IJK_Thermal_Subresolution::convert_into_sparse_matrix(const int& nb_points)
{
  /*
   * Convert into a huge sparse matrix
   */
  Matrice_Morse& sparse_matrix_solver  = ref_cast(Matrice_Morse, thermal_subproblems_matrix_assembly_for_solver_.valeur());
  Matrice_Bloc& bloc_matrix_solver = ref_cast(Matrice_Bloc, thermal_subproblems_matrix_assembly_.valeur());
  bloc_matrix_solver.block_to_morse(sparse_matrix_solver);
  Matrice_Morse& sparse_matrix_solver_reduced  = ref_cast(Matrice_Morse, thermal_subproblems_matrix_assembly_for_solver_reduced_.valeur());
  finite_difference_assembler_.reduce_solver_matrix(sparse_matrix_solver,
                                                    sparse_matrix_solver_reduced,
                                                    nb_points,
                                                    pre_initialise_thermal_subproblems_list_);
}

void IJK_Thermal_Subresolution::compute_md_vector()
{
  one_dimensional_advection_diffusion_thermal_solver_.reinit();
  ArrOfInt pe_voisins_dummy;
  ArrsOfInt items_to_send_dummy;
  ArrsOfInt items_to_recv_dummy;
  ArrsOfInt blocs_to_recv_dummy;
  pe_voisins_dummy.resize(0);
  items_to_send_dummy.resize(0);
  items_to_recv_dummy.resize(0);
  blocs_to_recv_dummy.resize(0);
//          const int last_subproblem_index = thermal_local_subproblems_.get_subproblems_counter();
//          const int subproblems_vector_size = thermal_local_subproblems_.get_end_index_subproblem(last_subproblem_index-1);
  const int subproblems_vector_size = thermal_subproblems_temperature_solution_for_solver_.size();
  MD_Vector_std md_std = MD_Vector_std(subproblems_vector_size, subproblems_vector_size,
                                       pe_voisins_dummy, items_to_send_dummy,
                                       items_to_recv_dummy, blocs_to_recv_dummy);
  md_.copy(md_std);
  MD_Vector_tools::creer_tableau_distribue(md_, thermal_subproblems_rhs_assembly_for_solver_); //, Array_base::NOCOPY_NOINIT);
  MD_Vector_tools::creer_tableau_distribue(md_, thermal_subproblems_temperature_solution_for_solver_); //, Array_base::NOCOPY_NOINIT);
}

void IJK_Thermal_Subresolution::retrieve_temperature_solution()
{
  /*
   * Retrieve the overall vector
   */
  for (int i=0; i<thermal_subproblems_temperature_solution_for_solver_.size(); i++)
    thermal_subproblems_temperature_solution_(i) = thermal_subproblems_temperature_solution_for_solver_(i);
  thermal_local_subproblems_.retrieve_radial_quantities();
}


void IJK_Thermal_Subresolution::check_wrong_values_rhs()
{
  if (debug_)
    {
      Cerr << "Check the modified RHS: INI" << finl;
      const double fd_second_order_magnitude = 1 / pow(dr_,2) * 1e3;
      for (int i=0; i<thermal_subproblems_rhs_assembly_.size(); i++)
        if (fabs(thermal_subproblems_rhs_assembly_[i]) > fd_second_order_magnitude)
          Cerr << "Check the modified RHS values: " << thermal_subproblems_rhs_assembly_[i] << finl;
      Cerr << "Check the modified RHS: END" << finl;
    }
}

/*
 * TODO:
 */

void IJK_Thermal_Subresolution::update_intersections()
{
  if (!disable_subresolution_)
    corrige_flux_->update_intersections();
}

void IJK_Thermal_Subresolution::prepare_thermal_flux_correction()
{
  if (!disable_subresolution_)
    corrige_flux_->update();
}

void IJK_Thermal_Subresolution::compute_min_max_reachable_fluxes()
{
  if (!disable_subresolution_ && find_reachable_fluxes_)
    corrige_flux_->compute_min_max_ijk_reachable_fluxes(cell_faces_neighbours_corrected_all_bool_, cell_faces_neighbours_corrected_min_max_bool_);
}

void IJK_Thermal_Subresolution::compute_convective_fluxes_face_centre()
{
  if (!disable_subresolution_ && convective_flux_correction_)
    corrige_flux_->compute_thermal_convective_fluxes();
}

void IJK_Thermal_Subresolution::compute_diffusive_fluxes_face_centre()
{
  if (!disable_subresolution_ && diffusive_flux_correction_)
    corrige_flux_->compute_thermal_diffusive_fluxes();
}

void IJK_Thermal_Subresolution::compute_temperature_cell_centres(const int first_corr)
{
  if (!disable_subresolution_)
    {
      switch(first_corr)
        {
        case 0:
          compute_temperature_cell_centres_first_correction();
          break;
        case 1:
          compute_temperature_cell_centres_second_correction();
          break;
        default:
          compute_temperature_cell_centres_first_correction();
          break;
        }
    }
}

void IJK_Thermal_Subresolution::compute_temperature_cell_centres_first_correction()
{
  corrige_flux_->compute_temperature_cell_centre(temperature_);
  corrige_flux_->compute_temperature_cell_centre_neighbours(temperature_cell_neighbours_,
                                                            neighbours_temperature_to_correct_,
                                                            neighbours_temperature_colinearity_weighting_);
  corrige_flux_->initialise_cell_neighbours_indices_to_correct();
  corrige_flux_->compute_cell_neighbours_faces_indices_for_spherical_correction(n_iter_distance_);
  corrige_flux_->compute_cell_neighbours_faces_indices_to_correct(cell_faces_neighbours_corrected_all_bool_,
                                                                  cell_faces_neighbours_corrected_convective_,
                                                                  cell_faces_neighbours_corrected_diffusive_,
                                                                  neighbours_faces_weighting_colinearity_,
                                                                  use_reachable_fluxes_);
  compute_min_max_reachable_fluxes();
  if (debug_)
    {
      temperature_cell_neighbours_debug_.data() = 0.;
      corrige_flux_->replace_temperature_cell_centre_neighbours(temperature_cell_neighbours_debug_,
                                                                temperature_cell_neighbours_,
                                                                neighbours_temperature_to_correct_,
                                                                neighbours_temperature_colinearity_weighting_);
    }
  int correct_first_iter = (correct_temperature_cell_neighbours_first_iter_
                            && ref_ijk_ft_->get_tstep() == 0
                            && !ref_ijk_ft_->get_reprise() != 0);
  if (use_temperature_cell_neighbours_)
    if (correct_temperature_cell_neighbours_first_iter_ == correct_first_iter)
      corrige_flux_->replace_temperature_cell_centre_neighbours(temperature_,
                                                                temperature_cell_neighbours_,
                                                                neighbours_temperature_to_correct_,
                                                                neighbours_temperature_colinearity_weighting_);
}

void IJK_Thermal_Subresolution::compute_temperature_cell_centres_second_correction()
{
  if (!disable_subresolution_ && disable_mixed_cells_increment_)
    {
      corrige_flux_->compute_temperature_cell_centre(temperature_);
    }
}

void IJK_Thermal_Subresolution::enforce_periodic_temperature_boundary_value()
{
  if (!disable_subresolution_ && enforce_periodic_boundary_value_)
    {
      const int ni = temperature_.ni();
      const int nj = temperature_.nj();
      const int nk = temperature_.nk();
      std::vector<double> indices_i_to_correct;
      std::vector<double> indices_j_to_correct;
      std::vector<double> indices_k_to_correct;
      int stencil_to_correct = 1;
      for (int i=0; i<=stencil_to_correct; i++)
        {
          indices_i_to_correct.push_back(stencil_to_correct);
          indices_j_to_correct.push_back(stencil_to_correct);
          indices_k_to_correct.push_back(stencil_to_correct);
          indices_i_to_correct.push_back(ni-1-stencil_to_correct);
          indices_j_to_correct.push_back(nj-1-stencil_to_correct);
          indices_k_to_correct.push_back(nk-1-stencil_to_correct);
        }
      int i,j,k;
      for (k = 0; k < nk; k++)
        if (std::find(indices_k_to_correct.begin(), indices_k_to_correct.end(), k) != indices_k_to_correct.end())
          for (j = 0; j < nj; j++)
            for (i = 0; i < ni; i++)
              {
                const double indic = ref_ijk_ft_->itfce().I(i,j,k);
                if (fabs(indic)<LIQUID_INDICATOR_TEST) // Mixed cells and pure vapour cells
                  { temperature_(i,j,k) = delta_T_subcooled_overheated_; }
              }
      for (j = 0; j < nj; j++)
        if (std::find(indices_j_to_correct.begin(), indices_j_to_correct.end(), j) != indices_j_to_correct.end())
          for (k = 0; k < nk; k++)
            for (i = 0; i < ni; i++)
              {
                const double indic = ref_ijk_ft_->itfce().I(i,j,k);
                if (fabs(indic)<LIQUID_INDICATOR_TEST) // Mixed cells and pure vapour cells
                  { temperature_(i,j,k) = delta_T_subcooled_overheated_; }
              }
      for (i = 0; i < ni; i++)
        if (std::find(indices_i_to_correct.begin(), indices_i_to_correct.end(), i) != indices_i_to_correct.end())
          for (k = 0; k < nk; k++)
            for (j = 0; j < nj; j++)
              {
                const double indic = ref_ijk_ft_->itfce().I(i,j,k);
                if (fabs(indic)<LIQUID_INDICATOR_TEST) // Mixed cells and pure vapour cells
                  { temperature_(i,j,k) = delta_T_subcooled_overheated_; }
              }
    }
}

void IJK_Thermal_Subresolution::correct_operators_for_visu()
{
  if (!diff_temperature_negligible_)
    {
      const int ni = div_coeff_grad_T_volume_.ni();
      const int nj = div_coeff_grad_T_volume_.nj();
      const int nk = div_coeff_grad_T_volume_.nk();
      for (int k = 0; k < nk; k++)
        for (int j = 0; j < nj; j++)
          for (int i = 0; i < ni; i++)
            {
              const double indic = ref_ijk_ft_->itfce().I(i,j,k);
              if (fabs(indic)<LIQUID_INDICATOR_TEST) // Mixed cells and pure vapour cells
                { div_coeff_grad_T_volume_(i,j,k) = 0; }
            }
    }
  if (!conv_temperature_negligible_)
    {
      const int ni = u_T_convective_volume_.ni();
      const int nj = u_T_convective_volume_.nj();
      const int nk = u_T_convective_volume_.nk();
      for (int k = 0; k < nk; k++)
        for (int j = 0; j < nj; j++)
          for (int i = 0; i < ni; i++)
            {
              const double indic = ref_ijk_ft_->itfce().I(i,j,k);
              if (fabs(indic)<LIQUID_INDICATOR_TEST) // Mixed cells and pure vapour cells
                { u_T_convective_volume_(i,j,k) = 0; }
            }
    }
}

void IJK_Thermal_Subresolution::prepare_ij_fluxes_k_layers()
{
  if ((!disable_subresolution_) && (convective_flux_correction_ || diffusive_flux_correction_))
    {
      corrige_flux_->compute_ijk_pure_faces_indices();
      corrige_flux_->sort_ijk_intersections_subproblems_indices_by_k_layers();
      if (store_cell_faces_corrected_)
        {
          corrige_flux_->store_cell_faces_corrected(cell_faces_corrected_bool_, cell_faces_corrected_convective_, cell_faces_corrected_diffusive_);
        }
    }
}

void IJK_Thermal_Subresolution::set_zero_temperature_increment()
{
  if (!disable_subresolution_ && disable_mixed_cells_increment_)
    corrige_flux_->set_zero_temperature_increment(d_temperature_);
}

void IJK_Thermal_Subresolution::clean_thermal_subproblems()
{
  if (!disable_subresolution_ && ref_ijk_ft_->get_tstep() > 0)
    {
      thermal_local_subproblems_.clean();
      corrige_flux_->clear_vectors();
    }
}

void IJK_Thermal_Subresolution::clean_ijk_intersections()
{
  if (!disable_subresolution_)
    corrige_flux_->clean();
}

void IJK_Thermal_Subresolution::set_thermal_subresolution_outputs(SFichier& fic)
{
  if (!disable_subresolution_)
    thermal_local_subproblems_.thermal_subresolution_outputs(fic, rang_);
}



