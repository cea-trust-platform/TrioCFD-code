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
#include <Corrige_flux_FT_base.h>
#include <OpConvDiscQuickIJKScalar.h>
#include <IJK_Ghost_Fluid_tools.h>

Implemente_instanciable_sans_constructeur( IJK_Thermal_Subresolution, "IJK_Thermal_Subresolution", IJK_Thermal_base ) ;

IJK_Thermal_Subresolution::IJK_Thermal_Subresolution()
{

  reference_gfm_on_probes_ = 0;
  compute_normal_derivatives_on_reference_probes_ = 0;

  disable_probe_weak_gradient_ = 0;
  disable_probe_weak_gradient_gfm_=0;

  reconstruct_previous_probe_field_=0;
  implicit_solver_from_previous_probe_field_=0;

  enable_probe_collision_detection_=0;
  enable_resize_probe_collision_=0;
  debug_probe_collision_=0;

  disable_spherical_diffusion_start_ = 0;
  probes_end_value_start_ = -1;
  probes_end_value_coeff_ = 0.05;
  single_centred_bubble_ = 1;
  computed_centred_bubble_start_ = 1.;
  single_centred_bubble_radius_ini_ = 1.e-3;
  temperature_ini_type_ = 1;
  time_ini_user_ = 0.;
  nusselt_spherical_diffusion_ = 2.;
  nusselt_spherical_diffusion_liquid_ = 2.;
  heat_flux_spherical_ = 0.;
  mean_liquid_temperature_ = -1;

  disable_mixed_cells_increment_=0;
  enable_mixed_cells_increment_=1;
  allow_temperature_correction_for_visu_=0;
  subproblem_temperature_extension_=0;
  disable_subresolution_=0;
  convective_flux_correction_ = 0;
  diffusive_flux_correction_ = 0;
  impose_fo_flux_correction_ = 1;
  disable_fo_flux_correction_=0;
  override_vapour_mixed_values_ = 0;
  compute_eulerian_compo_ = 1;

  error_temperature_ana_total_ = 0.;
  error_temperature_ana_squared_total_ = 0.;
  error_temperature_ana_rel_total_ = 0.;

  points_per_thermal_subproblem_ = 32;
  dr_ = 0.;
  probe_length_=0.;

  boundary_condition_interface_dict_ = Motcles(3);
  {
    boundary_condition_interface_dict_[0] = "dirichlet";
    boundary_condition_interface_dict_[1] = "neumann";
    boundary_condition_interface_dict_[2] = "flux_jump";
  }
  boundary_condition_interface_ = -1;
  interfacial_boundary_condition_value_ = 0.;
  impose_boundary_condition_interface_from_simulation_=0;
  boundary_condition_end_ = 0;
  boundary_condition_end_dict_ = Motcles(2);
  {
    boundary_condition_end_dict_[0] = "dirichlet";
    boundary_condition_end_dict_[1] = "neumann";
  }
  end_boundary_condition_value_ = -1;
  impose_user_boundary_condition_end_value_=0;

  source_terms_type_=2;
  source_terms_type_dict_ = Motcles(7);
  {
    source_terms_type_dict_[0] = "linear_diffusion";
    source_terms_type_dict_[1] = "spherical_diffusion";
    source_terms_type_dict_[2] = "spherical_diffusion_approx";
    source_terms_type_dict_[3] = "tangential_conv_2D";
    source_terms_type_dict_[4] = "tangential_conv_3D";
    source_terms_type_dict_[5] = "tangential_conv_2D_tangential_diffusion_3D";
    source_terms_type_dict_[6] = "tangential_conv_3D_tangentual_diffusion_3D";
  }
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
  compute_tangential_variables_ = 0;

  discrete_integral_ = 0;
  quadtree_levels_ = 1;
  advected_frame_of_reference_=0;
  neglect_frame_of_reference_radial_advection_=0;
  approximate_temperature_increment_=0;

  is_first_time_step_ = false;
  first_time_step_temporal_ = 0;
  first_time_step_implicit_ = 0;
  first_time_step_explicit_ = 1;
  local_fourier_ = 1.;
  local_cfl_ = 1.;

  distance_cell_faces_from_lrs_=1;
  disable_distance_cell_faces_from_lrs_=0;

  pre_initialise_thermal_subproblems_list_ = 0;
  pre_factor_subproblems_number_ = 3.;
  remove_append_subproblems_ = 0;

  correct_temperature_cell_neighbours_first_iter_ = 0;
  correct_first_iter_deactivate_cell_neighbours_ = 0;
  find_temperature_cell_neighbours_ = 0;
  correct_neighbours_using_probe_length_ = 0;
  neighbours_corrected_rank_ = 1;
  use_temperature_cell_neighbours_ = 0;
  neighbours_weighting_ = 0.;
  neighbours_colinearity_weighting_ = 0.;
  neighbours_distance_weighting_ = 0;
  neighbours_colinearity_distance_weighting_ = 0;
  smooth_temperature_field_=0;
  readjust_probe_length_from_vertices_=0;

  clip_temperature_values_ = 0;
  disable_post_processing_probes_out_files_ = 0;
  enforce_periodic_boundary_value_ = 0;
  stencil_periodic_boundary_value_ = 2;

  store_cell_faces_corrected_ = 0;

  find_cell_neighbours_for_fluxes_spherical_correction_ = 0;
  use_cell_neighbours_for_fluxes_spherical_correction_ = 0;
  find_reachable_fluxes_ = 0;
  use_reachable_fluxes_ = 0;
  keep_first_reachable_fluxes_ = 0;

  modified_time_init_ = 0.;
  spherical_diffusion_ = 1;

  post_process_all_probes_ = 0;
  nb_theta_post_pro_ = 10;
  nb_phi_post_pro_ = 4;
  nb_probes_post_pro_ = 40;

  interp_eulerian_ = 0;
  first_step_thermals_post_=1;
  disable_first_step_thermals_post_=0;

  neighbours_last_faces_weighting_ = 0;
  neighbours_last_faces_colinearity_weighting_ = 0;
  neighbours_last_faces_colinearity_face_weighting_ = 0;
  neighbours_last_faces_distance_weighting_ = 0;
  neighbours_last_faces_distance_colinearity_weighting_ = 0;
  neighbours_last_faces_distance_colinearity_face_weighting_ = 0;

  copy_fluxes_on_every_procs_ = 1;
  copy_temperature_on_every_procs_ = 1;

  initialise_sparse_indices_ = 0;
  use_sparse_matrix_ = 0;

  thermal_subproblems_matrix_assembly_for_solver_ref_ = nullptr;
}

Sortie& IJK_Thermal_Subresolution::printOn( Sortie& os ) const
{
  Nom front_space = "    ";
  Nom end_space = " ";
  Nom escape = "\n";

  IJK_Thermal_base::printOn( os );

  os << escape;
  os << front_space << "# SUBRESOLUTION #" << escape;
  os << escape;
  /*
   * Flags
   */
  os << escape;
  os << front_space << "# SUBRESOLUTION FLAGS #" << escape;
  os << escape;

  if (disable_probe_weak_gradient_)
    os << front_space << "disable_probe_weak_gradient " << escape;
  if (disable_probe_weak_gradient_gfm_)
    os << front_space << "disable_probe_weak_gradient_gfm " << escape;
  if (enable_probe_collision_detection_)
    os << front_space << "enable_probe_collision_detection " << escape;
  if (enable_resize_probe_collision_)
    os << front_space << "enable_resize_probe_collision " << escape;
  if (debug_probe_collision_)
    os << front_space << "debug_probe_collision " << escape;
  if (reconstruct_previous_probe_field_)
    os << front_space << "reconstruct_previous_probe_field " << escape;
  if (implicit_solver_from_previous_probe_field_)
    os << front_space << "implicit_solver_from_previous_probe_field " << escape;
  if (convective_flux_correction_)
    os << front_space << "convective_flux_correction " << escape;
  if (diffusive_flux_correction_)
    os << front_space << "diffusive_flux_correction " << escape;
  if (reference_gfm_on_probes_)
    os << front_space << "reference_gfm_on_probes " << escape;
  if (compute_normal_derivatives_on_reference_probes_)
    os << front_space << "compute_normal_derivatives_on_reference_probes" << escape;
  if (disable_spherical_diffusion_start_)
    os << front_space << "disable_spherical_diffusion_start" << escape;
  if (disable_subresolution_)
    os << front_space << "disable_subresolution" << escape;
  if (disable_fo_flux_correction_)
    os << front_space << "disable_fo_flux_correction" << escape;
  if (override_vapour_mixed_values_)
    os << front_space << "override_vapour_mixed_values" << escape;
  if (enable_mixed_cells_increment_)
    os << front_space << "enable_mixed_cells_increment" << escape;
  if (allow_temperature_correction_for_visu_)
    os << front_space << "allow_temperature_correction_for_visu" << escape;
  if (impose_boundary_condition_interface_from_simulation_)
    os << front_space << "impose_boundary_condition_interface_from_simulation" << escape;
  if (impose_user_boundary_condition_end_value_)
    os << front_space << "impose_user_boundary_condition_end_value" << escape;
  if (discrete_integral_)
    os << front_space << "discrete_integral" << escape;
  if (advected_frame_of_reference_)
    os << front_space << "advected_frame_of_reference" << escape;
  if (neglect_frame_of_reference_radial_advection_)
    os << front_space << "neglect_frame_of_reference_radial_advection" << escape;
  if (approximate_temperature_increment_)
    os << front_space << "approximate_temperature_increment" << escape;
  if (first_time_step_temporal_)
    os << front_space << "first_time_step_temporal" << escape;
  if (first_time_step_implicit_)
    os << front_space << "first_time_step_implicit" << escape;
  if (local_diffusion_fourier_priority_)
    os << front_space << "local_diffusion_fourier_priority" << escape;
  if (first_time_step_varying_probes_)
    os << front_space << "first_time_step_varying_probes" << escape;
  if (probe_variations_priority_)
    os << front_space << "probe_variations_priority" << escape;
  if (max_u_radial_)
    os << front_space << "max_u_radial" << escape;
  if (disable_distance_cell_faces_from_lrs_)
    os << front_space << "disable_distance_cell_faces_from_lrs" << escape;
  if (pre_initialise_thermal_subproblems_list_)
    os << front_space << "pre_initialise_thermal_subproblems_list" << escape;
  if (remove_append_subproblems_)
    os << front_space << "remove_append_subproblems" << escape;
  if (use_sparse_matrix_)
    os << front_space << "use_sparse_matrix" << escape;
  if (correct_temperature_cell_neighbours_first_iter_)
    os << front_space << "correct_temperature_cell_neighbours_first_iter" << escape;
  if (find_temperature_cell_neighbours_)
    os << front_space << "find_temperature_cell_neighbours" << escape;
  if (correct_neighbours_using_probe_length_)
    os << front_space << "correct_neighbours_using_probe_length" << escape;
  if (neighbours_colinearity_weighting_)
    os << front_space << "neighbours_colinearity_weighting" << escape;
  if (neighbours_distance_weighting_)
    os << front_space << "neighbours_distance_weighting" << escape;
  if (neighbours_colinearity_distance_weighting_)
    os << front_space << "neighbours_colinearity_distance_weighting" << escape;
  if (smooth_temperature_field_)
    os << front_space << "smooth_temperature_field" << escape;
  if (readjust_probe_length_from_vertices_)
    os << front_space << "readjust_probe_length_from_vertices" << escape;
  if (use_temperature_cell_neighbours_)
    os << front_space << "use_temperature_cell_neighbours" << escape;
  if (clip_temperature_values_)
    os << front_space << "clip_temperature_values" << escape;
  if (disable_post_processing_probes_out_files_)
    os << front_space << "disable_post_processing_probes_out_files" << escape;
  if (enforce_periodic_boundary_value_)
    os << front_space << "enforce_periodic_boundary_value" << escape;
  if (disable_post_processing_probes_out_files_)
    os << front_space << "disable_post_processing_probes_out_files" << escape;
  if (enforce_periodic_boundary_value_)
    os << front_space << "enforce_periodic_boundary_value" << escape;
  if (store_cell_faces_corrected_)
    os << front_space << "store_cell_faces_corrected" << escape;
  if (find_cell_neighbours_for_fluxes_spherical_correction_)
    os << front_space << "find_cell_neighbours_for_fluxes_spherical_correction" << escape;
  if (use_cell_neighbours_for_fluxes_spherical_correction_)
    os << front_space << "use_cell_neighbours_for_fluxes_spherical_correction" << escape;
  if (find_reachable_fluxes_)
    os << front_space << "find_reachable_fluxes" << escape;
  if (use_reachable_fluxes_)
    os << front_space << "use_reachable_fluxes" << escape;
  if (keep_first_reachable_fluxes_)
    os << front_space << "keep_first_reachable_fluxes" << escape;
  if (post_process_all_probes_)
    os << front_space << "post_process_all_probes" << escape;
  if (interp_eulerian_)
    os << front_space << "interp_eulerian" << escape;
  if (first_step_thermals_post_)
    os << front_space << "first_step_thermals_post" << escape;
  if (neighbours_last_faces_colinearity_weighting_)
    os << front_space << "neighbours_last_faces_colinearity_weighting" << escape;
  if (neighbours_last_faces_colinearity_face_weighting_)
    os << front_space << "neighbours_last_faces_colinearity_face_weighting" << escape;
  if (neighbours_last_faces_distance_weighting_)
    os << front_space << "neighbours_last_faces_distance_weighting" << escape;
  if (neighbours_last_faces_distance_weighting_)
    os << front_space << "neighbours_last_faces_distance_weighting" << escape;
  if (neighbours_last_faces_distance_colinearity_face_weighting_)
    os << front_space << "neighbours_last_faces_distance_colinearity_face_weighting" << escape;

  /*
   * Values
   */
  os << escape;
  os << front_space << "# SUBRESOLUTION PARAMS #" << escape;
  os << escape;

  os << front_space << "points_per_thermal_subproblem" << end_space << points_per_thermal_subproblem_ << escape;
  os << front_space << "coeff_distance_diagonal" << end_space << coeff_distance_diagonal_ << escape;
  os << front_space << "finite_difference_assembler" << end_space << finite_difference_assembler_ << escape;
  if (boundary_condition_interface_ != -1)
    os << front_space << "boundary_condition_interface" << end_space << boundary_condition_interface_dict_[boundary_condition_interface_] << escape;
  os << front_space << "interfacial_boundary_condition_value" << end_space << interfacial_boundary_condition_value_ << escape;
  if (boundary_condition_end_ != -1)
    os << front_space << "boundary_condition_end" << end_space << boundary_condition_end_dict_[boundary_condition_end_] << escape;
  os << front_space << "end_boundary_condition_value" << end_space << end_boundary_condition_value_ << escape;
  if (source_terms_type_!= -1)
    os << front_space << "source_terms_type" << end_space << source_terms_type_dict_[source_terms_type_] << escape;
  os << front_space << "thermal_fd_solver" << end_space << one_dimensional_advection_diffusion_thermal_solver_ << escape;
  os << front_space << "quadtree_levels" << end_space << quadtree_levels_ << escape;
  os << front_space << "local_fourier" << end_space << local_fourier_ << escape;
  os << front_space << "local_cfl" << end_space << local_cfl_ << escape;
  os << front_space << "delta_T_subcooled_overheated" << end_space << delta_T_subcooled_overheated_ << escape;
  os << front_space << "pre_factor_subproblems_number" << end_space << pre_factor_subproblems_number_ << escape;
  os << front_space << "neighbours_corrected_rank" << end_space << neighbours_corrected_rank_ << escape;
  os << front_space << "probes_end_value_coeff" << end_space << probes_end_value_coeff_ << escape;
  os << front_space << "temperature_ini_type" << end_space << temperature_ini_type_ << escape;
  os << front_space << "time_ini_user" << end_space << time_ini_user_ << escape;
  os << front_space << "nb_theta_post_pro" << end_space << nb_theta_post_pro_ << escape;
  os << front_space << "nb_phi_post_pro" << end_space << nb_phi_post_pro_ << escape;
  os << front_space << "nb_probes_post_pro" << end_space << nb_probes_post_pro_ << escape;
  os << front_space << "modified_time_init" << end_space << modified_time_init_ << escape;

  os << "}" << escape;
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

  param.ajouter("modified_time_init", &modified_time_init_);

  param.ajouter_flag("reference_gfm_on_probes", &reference_gfm_on_probes_);
  param.ajouter_flag("compute_normal_derivatives_on_reference_probes", &compute_normal_derivatives_on_reference_probes_);

  param.ajouter_flag("disable_probe_weak_gradient", &disable_probe_weak_gradient_);
  param.ajouter_flag("disable_probe_weak_gradient_gfm", &disable_probe_weak_gradient_gfm_);

  param.ajouter_flag("enable_probe_collision_detection", &enable_probe_collision_detection_);
  param.ajouter_flag("enable_resize_probe_collision", &enable_resize_probe_collision_);
  param.ajouter_flag("debug_probe_collision", &debug_probe_collision_);

  param.ajouter_flag("reconstruct_previous_probe_field", &reconstruct_previous_probe_field_);
  param.ajouter_flag("implicit_solver_from_previous_probe_field", &implicit_solver_from_previous_probe_field_);

  param.ajouter_flag("disable_spherical_diffusion_start", &disable_spherical_diffusion_start_);
  param.ajouter_flag("disable_subresolution", &disable_subresolution_);
  param.ajouter_flag("convective_flux_correction", &convective_flux_correction_);
  param.ajouter_flag("diffusive_flux_correction", &diffusive_flux_correction_);
  param.ajouter_flag("disable_fo_flux_correction", &disable_fo_flux_correction_);
  param.ajouter_flag("override_vapour_mixed_values", &override_vapour_mixed_values_);
  param.ajouter("points_per_thermal_subproblem", &points_per_thermal_subproblem_);
  param.ajouter("coeff_distance_diagonal", &coeff_distance_diagonal_);
  param.ajouter("finite_difference_assembler", &finite_difference_assembler_);
  // param.ajouter_flag("disable_mixed_cells_increment", &disable_mixed_cells_increment_);
  param.ajouter_flag("enable_mixed_cells_increment", &enable_mixed_cells_increment_);
  param.ajouter_flag("disable_mixed_cells_increment", &disable_mixed_cells_increment_);
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

  param.ajouter_flag("disable_distance_cell_faces_from_lrs", &disable_distance_cell_faces_from_lrs_);

  param.ajouter_flag("pre_initialise_thermal_subproblems_list", &pre_initialise_thermal_subproblems_list_);
  param.ajouter("pre_factor_subproblems_number", &pre_factor_subproblems_number_);
  param.ajouter_flag("remove_append_subproblems", &remove_append_subproblems_);
  param.ajouter_flag("use_sparse_matrix", &use_sparse_matrix_);

  param.ajouter_flag("correct_temperature_cell_neighbours_first_iter", &correct_temperature_cell_neighbours_first_iter_);
  param.ajouter_flag("find_temperature_cell_neighbours", &find_temperature_cell_neighbours_);
  param.ajouter_flag("correct_neighbours_using_probe_length", &correct_neighbours_using_probe_length_);
  param.ajouter("neighbours_corrected_rank", &neighbours_corrected_rank_);

  param.ajouter_flag("neighbours_colinearity_weighting", &neighbours_colinearity_weighting_);
  param.ajouter_flag("neighbours_distance_weighting", &neighbours_distance_weighting_);
  param.ajouter_flag("neighbours_colinearity_distance_weighting", &neighbours_colinearity_distance_weighting_);
  param.ajouter_flag("smooth_temperature_field", &smooth_temperature_field_);
  param.ajouter_flag("readjust_probe_length_from_vertices", &readjust_probe_length_from_vertices_);
  param.ajouter_flag("use_temperature_cell_neighbours", &use_temperature_cell_neighbours_);

  param.ajouter_flag("clip_temperature_values", &clip_temperature_values_);
  param.ajouter_flag("disable_post_processing_probes_out_files", &disable_post_processing_probes_out_files_);

  param.ajouter_flag("enforce_periodic_boundary_value", &enforce_periodic_boundary_value_);
  param.ajouter_flag("store_cell_faces_corrected", &store_cell_faces_corrected_);

  param.ajouter_flag("find_cell_neighbours_for_fluxes_spherical_correction", &find_cell_neighbours_for_fluxes_spherical_correction_);
  param.ajouter_flag("use_cell_neighbours_for_fluxes_spherical_correction", &use_cell_neighbours_for_fluxes_spherical_correction_);

  param.ajouter_flag("find_reachable_fluxes", &find_reachable_fluxes_);
  param.ajouter_flag("use_reachable_fluxes", &use_reachable_fluxes_);
  param.ajouter_flag("keep_first_reachable_fluxes", &keep_first_reachable_fluxes_);

  param.ajouter("probes_end_value_coeff", &probes_end_value_coeff_);
  param.ajouter("temperature_ini_type", &temperature_ini_type_);
  param.ajouter("time_ini_user", &time_ini_user_);

  param.ajouter_flag("post_process_all_probes", &post_process_all_probes_);
  param.ajouter("nb_theta_post_pro", &nb_theta_post_pro_);
  param.ajouter("nb_phi_post_pro", &nb_phi_post_pro_);
  param.ajouter("nb_probes_post_pro", &nb_probes_post_pro_);

  param.ajouter_flag("interp_eulerian", &interp_eulerian_);

  param.ajouter_flag("disable_first_step_thermals_post", &disable_first_step_thermals_post_);

  param.ajouter_flag("neighbours_last_faces_colinearity_weighting", &neighbours_last_faces_colinearity_weighting_);
  param.ajouter_flag("neighbours_last_faces_colinearity_face_weighting", &neighbours_last_faces_colinearity_face_weighting_);
  param.ajouter_flag("neighbours_last_faces_distance_weighting", &neighbours_last_faces_distance_weighting_);
  param.ajouter_flag("neighbours_last_faces_distance_colinearity_weighting", &neighbours_last_faces_distance_colinearity_weighting_);
  param.ajouter_flag("neighbours_last_faces_distance_colinearity_face_weighting", &neighbours_last_faces_distance_colinearity_face_weighting_);

//  param.ajouter_flag("copy_fluxes_on_every_procs", &copy_fluxes_on_every_procs_);
//  param.ajouter_flag("copy_temperature_on_every_procs", &copy_temperature_on_every_procs_);

// param.ajouter_flag("enforce_periodic_boundary_value", &enforce_periodic_boundary_value_);
// param.ajouter_non_std("enforce_periodic_boundary_value",(this));
//  for (int i=0; i<fd_solvers_jdd_.size(); i++)
//    param.ajouter_non_std(fd_solvers_jdd_[i], this);
}

int IJK_Thermal_Subresolution::lire_motcle_non_standard(const Motcle& mot, Entree& is)
{
  if (mot == "enforce_periodic_boundary_value")
    {

      Nom accolade_ouverte("{");
      Param param_bloc(que_suis_je());
      Param param(que_suis_je());
      Nom bloc_name_next;
      param_bloc.ajouter("enforce_periodic_boundary_value",&bloc_name_next, Param::REQUIRED);
      param_bloc.lire_sans_accolade(is);
      if (bloc_name_next == accolade_ouverte)
        {
          enforce_periodic_boundary_value_ = 1;
          param.ajouter("stencil_periodic_boundary_value", &stencil_periodic_boundary_value_);
          param.lire_avec_accolades_depuis(is);
        }
      else
        enforce_periodic_boundary_value_ = 1;
    }
  return 1;
}

//int IJK_Thermal_Subresolution::lire_motcle_non_standard(const Motcle& mot, Entree& is)
//{
//  read_fd_solver(mot, is);
//  return 1;
//}
//
//void IJK_Thermal_Subresolution::read_fd_solver(const Motcle& mot, Entree& is)
//{
//  Nom type = "";
//  fd_solver_rank_ = fd_solvers_jdd_.search(mot);
//  type += fd_solvers_[fd_solver_rank_];
//  Cerr << "Finite difference solver: " << type << finl;
//  fd_solver_type_ = fd_solvers_[fd_solver_rank_];
//  one_dimensional_advection_diffusion_thermal_solver_.typer(type);
//  // is >> one_dimensional_advection_diffusion_thermal_solver_;
//  // Cerr << "Finite difference solver: " << fd_solver_type_ << finl;
//}


int IJK_Thermal_Subresolution::initialize(const IJK_Splitting& splitting, const int idx)
{
  Cout << que_suis_je() << "::initialize()" << finl;

  uniform_lambda_ = lambda_liquid_;
  uniform_alpha_ =	lambda_liquid_ / (ref_ijk_ft_->get_rho_l() * cp_liquid_);
  prandtl_number_ = (ref_ijk_ft_->get_mu_liquid() / ref_ijk_ft_->get_rho_l()) / uniform_alpha_;
  //  calulate_grad_T_ = 1;
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
  compute_rising_velocities_ = 1;

  distance_cell_faces_from_lrs_ = !disable_distance_cell_faces_from_lrs_;
  first_step_thermals_post_ = !disable_first_step_thermals_post_;

  if (compute_normal_derivatives_on_reference_probes_)
    reference_gfm_on_probes_ = 1;

  if (use_sparse_matrix_)
    {
      pre_initialise_thermal_subproblems_list_ = 0;
      remove_append_subproblems_ = 1;
    }

  int nalloc = 0;
  nalloc = IJK_Thermal_base::initialize(splitting, idx);
  corrige_flux_.typer("Corrige_flux_FT_temperature_subresolution");

  temperature_before_extrapolation_.allocate(splitting, IJK_Splitting::ELEM, ghost_cells_);
  nalloc += 1;
  temperature_before_extrapolation_.data() = 0.;

  neighbours_weighting_ = (neighbours_colinearity_weighting_ || neighbours_distance_weighting_ || neighbours_colinearity_distance_weighting_);
  neighbours_last_faces_weighting_ = (neighbours_last_faces_colinearity_weighting_ || neighbours_last_faces_colinearity_face_weighting_
                                      || neighbours_last_faces_distance_weighting_ || neighbours_last_faces_distance_colinearity_weighting_
                                      || neighbours_last_faces_distance_colinearity_face_weighting_);
  if (neighbours_weighting_)
    if (!neighbours_last_faces_weighting_)
      {
        neighbours_last_faces_weighting_ = 1;
        neighbours_last_faces_colinearity_weighting_ = 1;
      }
  if(neighbours_last_faces_weighting_)
    if (!neighbours_weighting_)
      {
        neighbours_weighting_ = 1;
        neighbours_colinearity_weighting_ = 1;
      }

  correct_temperature_cell_neighbours_first_iter_ = (correct_temperature_cell_neighbours_first_iter_ && disable_spherical_diffusion_start_);
  if (correct_temperature_cell_neighbours_first_iter_)
    use_temperature_cell_neighbours_ = correct_temperature_cell_neighbours_first_iter_;

  if (use_temperature_cell_neighbours_)
    find_temperature_cell_neighbours_ = use_temperature_cell_neighbours_;

  find_cell_neighbours_for_fluxes_spherical_correction_ = find_cell_neighbours_for_fluxes_spherical_correction_ && distance_cell_faces_from_lrs_;
  use_cell_neighbours_for_fluxes_spherical_correction_ = use_cell_neighbours_for_fluxes_spherical_correction_ && distance_cell_faces_from_lrs_;
  if (use_cell_neighbours_for_fluxes_spherical_correction_)
    find_cell_neighbours_for_fluxes_spherical_correction_ = use_cell_neighbours_for_fluxes_spherical_correction_;
  corrige_flux_->set_correction_cell_faces_neighbours(find_cell_neighbours_for_fluxes_spherical_correction_,
                                                      use_cell_neighbours_for_fluxes_spherical_correction_,
                                                      find_reachable_fluxes_,
                                                      use_reachable_fluxes_,
                                                      keep_first_reachable_fluxes_);

  temperature_diffusion_op_.set_conductivity_coefficient(uniform_lambda_, temperature_, temperature_, temperature_, temperature_);
  if (debug_)
    Cerr << "Uniform lambda: " << temperature_diffusion_op_.get_uniform_lambda() << finl;

  impose_fo_flux_correction_ = !disable_fo_flux_correction_;
  convective_flux_correction_ = convective_flux_correction_ && (!conv_temperature_negligible_);
  diffusive_flux_correction_ = diffusive_flux_correction_ && (!diff_temperature_negligible_);

  use_reachable_fluxes_ = (use_reachable_fluxes_ && find_reachable_fluxes_) && (convective_flux_correction_ || diffusive_flux_correction_);
  if (use_reachable_fluxes_)
    {
      find_temperature_cell_neighbours_ = 1;
      use_temperature_cell_neighbours_ = 1;
    }
  if(find_reachable_fluxes_)
    find_temperature_cell_neighbours_ = 1;

  if (enable_mixed_cells_increment_)
    disable_mixed_cells_increment_ = (!enable_mixed_cells_increment_);
  if (!disable_mixed_cells_increment_)
    enable_mixed_cells_increment_ = (!disable_mixed_cells_increment_);

  corrige_flux_->set_convection_diffusion_correction(convective_flux_correction_, diffusive_flux_correction_);
  corrige_flux_->set_fluxes_feedback_params(discrete_integral_, quadtree_levels_);
  corrige_flux_->set_debug(debug_);
  corrige_flux_->set_distance_cell_faces_from_lrs(distance_cell_faces_from_lrs_);
  corrige_flux_->set_correction_cell_neighbours(find_temperature_cell_neighbours_,
                                                neighbours_weighting_,
                                                smooth_temperature_field_);
  corrige_flux_->set_eulerian_normal_vectors_ns_normed(eulerian_normal_vectors_ns_normed_);
  corrige_flux_->set_temperature_fluxes_periodic_sharing_strategy_on_processors(copy_fluxes_on_every_procs_,
                                                                                copy_temperature_on_every_procs_);


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
  if (readjust_probe_length_from_vertices_)
    probe_variations_enabled_ = 1;

  if (one_dimensional_advection_diffusion_thermal_solver_.est_nul())
    one_dimensional_advection_diffusion_thermal_solver_.cast_direct_solver_by_default();

  /*
   * Considered constant values
   */
  corrige_flux_->set_physical_parameters(cp_liquid_ * ref_ijk_ft_->get_rho_l(),
                                         cp_vapour_ * ref_ijk_ft_->get_rho_v(),
                                         lambda_liquid_,
                                         lambda_vapour_);
  corrige_flux_->initialize_with_subproblems(
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
      neighbours_temperature_to_correct_.allocate(splitting, IJK_Splitting::ELEM, ghost_cells_); // , 1);
      temperature_cell_neighbours_.allocate(splitting, IJK_Splitting::ELEM, ghost_cells_); // , 1);
      nalloc += 2;
      neighbours_temperature_to_correct_.data() = 0;
      temperature_cell_neighbours_.data() = 0.;
      neighbours_temperature_to_correct_.echange_espace_virtuel(neighbours_temperature_to_correct_.ghost());
      temperature_cell_neighbours_.echange_espace_virtuel(temperature_cell_neighbours_.ghost());
      if (neighbours_weighting_)
        {
          neighbours_temperature_colinearity_weighting_.allocate(splitting, IJK_Splitting::ELEM, ghost_cells_); // , 1);
          nalloc += 1;
          neighbours_temperature_colinearity_weighting_.data() = 0.;
          neighbours_temperature_colinearity_weighting_.echange_espace_virtuel(neighbours_temperature_colinearity_weighting_.ghost());
        }
      if (debug_)
        {
          temperature_cell_neighbours_debug_.allocate(splitting, IJK_Splitting::ELEM, ghost_cells_); // , 1);
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
      allocate_cell_vector(cell_faces_neighbours_corrected_all_bool_, splitting, ghost_cells_); // , 1);
      nalloc += 3;
      for (int c=0; c<3; c++)
        cell_faces_neighbours_corrected_all_bool_[c].data() = 0;
      cell_faces_neighbours_corrected_all_bool_.echange_espace_virtuel();

      allocate_cell_vector(cell_faces_neighbours_corrected_min_max_bool_, splitting, 2);
      nalloc += 3;
      for (int c=0; c<3; c++)
        cell_faces_neighbours_corrected_min_max_bool_[c].data() = 0;
      cell_faces_neighbours_corrected_min_max_bool_.echange_espace_virtuel();

      neighbours_temperature_to_correct_trimmed_.allocate(splitting, IJK_Splitting::ELEM, 1);
      nalloc += 1;
      neighbours_temperature_to_correct_trimmed_.data() = 0.;
      neighbours_temperature_to_correct_trimmed_.echange_espace_virtuel(neighbours_temperature_to_correct_trimmed_.ghost());
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
      if (neighbours_weighting_)
        {
          allocate_cell_vector(neighbours_faces_weighting_colinearity_, splitting, 1);
          nalloc += 3;
          for (int c=0; c<3; c++)
            neighbours_faces_weighting_colinearity_[c].data() = 0.;
          neighbours_faces_weighting_colinearity_.echange_espace_virtuel();
        }
    }
  if (debug_probe_collision_)
    {
      probe_collision_debug_field_.allocate(splitting, IJK_Splitting::ELEM, 0); // , 1);
      nalloc += 1;
      probe_collision_debug_field_.data() = 0.;
    }

  if (impose_fo_flux_correction_ || (diffusive_flux_correction_ && fo_ >= 1.)) // By default ?
    fo_ = 0.95 * pow((sqrt(2) / 2), 2); // squared or not?

  Cout << "End of " << que_suis_je() << "::initialize()" << finl;
  return nalloc;
}

double IJK_Thermal_Subresolution::get_modified_time()
{
  if (!disable_spherical_diffusion_start_)
    return modified_time_init_; // ref_ijk_ft_->get_current_time();
  else
    return ref_ijk_ft_->get_current_time();
}

static int decoder_numero_bulle(const int code)
{
  const int num_bulle = code >> 6;
  return num_bulle;
}

void IJK_Thermal_Subresolution::compute_temperature_init()
{
  /*
   * Theses Thiam et Euzenat
   */
  if (!disable_spherical_diffusion_start_)
    {
      if (single_centred_bubble_)
        {
          double time_ini = 0.;
          const double rho_l = ref_ijk_ft_->get_rho_l();
          const double alpha_liq = lambda_liquid_ / (rho_l * cp_liquid_);
          switch(temperature_ini_type_)
            {
            case local_criteria:
              {
                const double R = get_probes_length() +  single_centred_bubble_radius_ini_;
                const double erf_val = 1 - probes_end_value_coeff_ * abs(delta_T_subcooled_overheated_) * R / single_centred_bubble_radius_ini_;
                double erf_inv_val = 1.;
                approx_erf_inverse(erf_val, erf_inv_val);
                const double time_local = pow((R-single_centred_bubble_radius_ini_),2) / (pow(erf_inv_val,2) * 4. * alpha_liq) ;
                time_ini = time_local;
              }
              break;
            case integral_criteria:
              {
                const int max_attempt = 50;
                int attempt_counter = 0;
                double temperature_end_prev = delta_T_subcooled_overheated_;
                double time_integral = 0.;
                double temperature_end_next = 2 * temperature_end_prev;
                while (abs(temperature_end_next - temperature_end_prev) > 0.01 && (attempt_counter < max_attempt))
                  {
                    const double temperature_integral = compute_spherical_steady_dirichlet_left_right_integral(temperature_end_prev);
                    time_integral = find_time_dichotomy_integral(temperature_integral, temperature_end_next);
                    temperature_end_prev = 0.5 * (temperature_end_prev + temperature_end_next);
                    attempt_counter++;
                  }
                time_ini = time_integral;
              }
              break;
            case derivative_criteria:
              {
                double time_derivative=0.;
                const int max_attempt = 100;
                int attempt_counter = 0;
                const double R = get_probes_length() +  single_centred_bubble_radius_ini_;
                double error = 1.;
                double temperature_end_min = 0.;
                get_time_inflection_derivative(temperature_end_min);
                double temperature_limit_left = temperature_end_min;
                double temperature_limit_right = temperature_end_min + 0.25 * abs(temperature_end_min);
                double temperature_middle = 0.5 * (temperature_limit_left + temperature_limit_right);
                while (error > 1e-2 && (attempt_counter < max_attempt))
                  {
                    const double temperature_derivative = compute_spherical_steady_dirichlet_left_right_derivative_value(R, temperature_middle);
                    time_derivative = find_time_dichotomy_derivative(temperature_derivative, temperature_limit_left, temperature_limit_right);
                    error = abs(temperature_limit_right - temperature_limit_left);
                    if (debug_)
                      Cerr << "error: " << error << finl;
                    temperature_middle = 0.5 * (temperature_limit_right + temperature_limit_right);
                    attempt_counter++;
                  }
                time_ini = time_derivative;
              }
              break;
            case time_criteria:
              {
                time_ini = time_ini_user_;
              }
              break;
            default:
              break;
            }
          if (debug_)
            Cerr << "Time ini: " << time_ini << finl;
          const int nb_bubble_tot = ref_ijk_ft_->itfce().get_ijk_compo_connex().get_bubbles_barycentre().dimension(0);
          const int nb_bubbles_real = ref_ijk_ft_->itfce().get_nb_bulles_reelles();
          for (int index_bubble=0; index_bubble<nb_bubble_tot; index_bubble++)
            {
              int index_bubble_real = index_bubble;
              if (index_bubble>=nb_bubbles_real)
                {
                  const int ighost = ref_ijk_ft_->itfce().ghost_compo_converter(index_bubble-nb_bubbles_real);
                  index_bubble_real = decoder_numero_bulle(-ighost);
                }
              Nom expression_T_ini = compute_quasi_static_spherical_diffusion_expression(time_ini, index_bubble, index_bubble_real);
              set_field_data(temperature_for_ini_per_bubble_, expression_T_ini);
              set_field_temperature_per_bubble(index_bubble);
            }
          correct_temperature_for_visu();
          modified_time_init_ = time_ini;
        }
    }
  else
    IJK_Thermal_base::compute_temperature_init();
}

Nom IJK_Thermal_Subresolution::compute_quasi_static_spherical_diffusion_expression(const double& time_scope,
                                                                                   const int index_bubble,
                                                                                   const int index_bubble_real)
{

  if (computed_centred_bubble_start_)
    {
      const DoubleTab& bubbles_centres = ref_ijk_ft_->itfce().get_ijk_compo_connex().get_bubbles_barycentre();
      double x,y,z;
      x = bubbles_centres(index_bubble,0);
      y = bubbles_centres(index_bubble,1);
      z = bubbles_centres(index_bubble,2);
      if (index_bubble != index_bubble_real)
        {
          const double x_real = bubbles_centres(index_bubble_real,0);
          const double y_real = bubbles_centres(index_bubble_real,1);
          const double z_real = bubbles_centres(index_bubble_real,2);
          const double lx = temperature_.get_splitting().get_grid_geometry().get_domain_length(0);
          const double ly = temperature_.get_splitting().get_grid_geometry().get_domain_length(0);
          const double lz = temperature_.get_splitting().get_grid_geometry().get_domain_length(0);
          x = (abs(x - x_real)< (lx / 2.)) ? x_real : ((x_real < (lx / 2.)) ? x_real + lx : x_real - lx);
          y = (abs(y - y_real)< (ly / 2.)) ? y_real : ((y_real < (ly / 2.)) ? y_real + ly : y_real - ly);
          z = (abs(z - z_real)< (lz / 2.)) ? z_real : ((z_real < (lz / 2.)) ? z_real + lz : z_real - lz);
        }
      return generate_expression_temperature_ini(time_scope, x, y, z);
    }
  return generate_expression_temperature_ini(time_scope, 0., 0., 0.);
}

void IJK_Thermal_Subresolution::set_field_temperature_per_bubble(const int index_bubble)
{
  if (!index_bubble)
    temperature_.data() = delta_T_subcooled_overheated_;
  const int sign_delta = signbit(delta_T_subcooled_overheated_);
  const int ni = temperature_.ni();
  const int nj = temperature_.nj();
  const int nk = temperature_.nk();
  for (int k = 0; k < nk; k++)
    for (int j = 0; j < nj; j++)
      for (int i = 0; i < ni; i++)
        {
          const double indic = ref_ijk_ft_->itfce().I()(i,j,k);
          if (indic > VAPOUR_INDICATOR_TEST)
            {
              const double temperature = temperature_(i,j,k);
              if (sign_delta)
                {
                  if (temperature <= (1 - probes_end_value_coeff_) * delta_T_subcooled_overheated_)
                    {
                      const double temperature_per_bubble = temperature_for_ini_per_bubble_(i,j,k);
                      temperature_(i,j,k) = temperature_per_bubble;
                    }
                }
              else
                {
                  if (temperature >= (1 - probes_end_value_coeff_) * delta_T_subcooled_overheated_)
                    {
                      const double temperature_per_bubble = temperature_for_ini_per_bubble_(i,j,k);
                      temperature_(i,j,k) = temperature_per_bubble;
                    }
                }
            }
        }
}

void reinit_streamObj(std::ostringstream& streamObj, const double& param)
{
  streamObj.str("");
  streamObj.clear();
  streamObj << (double) param;
}

Nom IJK_Thermal_Subresolution::generate_expression_temperature_ini(const double& time_scope, const double x, const double y, const double z)
{
  const double rho_l = ref_ijk_ft_->get_rho_l();
  const double alpha_liq = lambda_liquid_ / (rho_l * cp_liquid_);
  std::ostringstream streamObj;
  Nom expression_T = "(";
  streamObj << delta_T_subcooled_overheated_;
  expression_T += streamObj.str().c_str();
  expression_T += ")-(";
  expression_T += streamObj.str().c_str();
  expression_T += ")*(";
  reinit_streamObj(streamObj, single_centred_bubble_radius_ini_);
  expression_T += streamObj.str().c_str();
  Nom expression_tmp = "sqrt((x-(";
  reinit_streamObj(streamObj, (signbit(x) ? x-1e-16 : x+1e-16));
  expression_tmp += streamObj.str().c_str();
  expression_tmp += "))^2+(y-(";
  reinit_streamObj(streamObj, (signbit(y) ? y-1e-16 : y+1e-16));
  expression_tmp += streamObj.str().c_str();
  expression_tmp += "))^2+(z-(";
  reinit_streamObj(streamObj, (signbit(z) ? z-1e-16 : z+1e-16));
  expression_tmp += streamObj.str().c_str();
  expression_tmp += "))^2)";
  expression_T += "/";
  expression_T += expression_tmp;
  expression_T += "*(1-erf((";
  expression_T += expression_tmp;
  expression_T += "-";
  reinit_streamObj(streamObj, single_centred_bubble_radius_ini_);
  expression_T += streamObj.str().c_str();
  expression_T += ")/(2.*sqrt(";
  streamObj.str("");
  streamObj.clear();
  reinit_streamObj(streamObj, ((alpha_liq * time_scope) + 1e-16));
  expression_T += streamObj.str().c_str();
  expression_T += ")))))";
  return expression_T;
}

void IJK_Thermal_Subresolution::approx_erf_inverse(const double& x, double& res)
{
  // https://www.had2know.org/academics/error-function-calculator-erf.html
  const double a = 18.75537;
  const double b = 2.47197;
  const double c = 0.25;
  const double d = 4.33074;
  const double e = 0.5;
  const double w = log(1-pow(x,2));
  res = sqrt(sqrt(a - b * w + c * pow(w,2)) - d - e * w);
  const double sign = signbit(x) ? -1. : 1.;
  res *= sign;
}

double IJK_Thermal_Subresolution::compute_spherical_steady_dirichlet_left_right_value(const double& r)
{
  double temperature_value;
  const double T1 = delta_T_subcooled_overheated_;
  const double R0 = single_centred_bubble_radius_ini_;
  const double R1 = get_probes_length() + R0;
  temperature_value = T1 + T1 * (r * R0 - R1 * R0) / (r * (R1 - R0));
  return temperature_value;
}

double IJK_Thermal_Subresolution::compute_spherical_steady_dirichlet_left_right_derivative_value(const double& r, const double& temperature_end_prev)
{
  double temperature_derivative;
  const double T1 = temperature_end_prev;
  const double R0 = single_centred_bubble_radius_ini_;
  const double R1 = get_probes_length() + R0;
  temperature_derivative = T1 * (R1 * R0) / (pow(r, 2) * (R1 - R0));
  return temperature_derivative;
}

double IJK_Thermal_Subresolution::compute_spherical_steady_dirichlet_left_right_integral(const double& temperature_end_prev)
{
  double temperature_integral;
  const double T1 = temperature_end_prev;
  const double R0 = single_centred_bubble_radius_ini_;
  const double R1 = get_probes_length() + R0;
  const double Delta_R = R1 - R0;
  temperature_integral = (T1 * R1) / (R1 - R0) * (Delta_R - R0 * (log(R1) - log(R0)));
  temperature_integral /= Delta_R;
  return temperature_integral;
}

double IJK_Thermal_Subresolution::find_time_dichotomy_integral(const double& temperature_integral,
                                                               double& temperature_end_prev)
{
  // Arbitrary large time
  const double time_integral = 100.;
  const double T1 = delta_T_subcooled_overheated_;
  const double R0 = single_centred_bubble_radius_ini_;
  const double R1 = get_probes_length() + R0;
  const double rho_l = ref_ijk_ft_->get_rho_l();
  const double Delta_R = R1 - R0;
  const double alpha_liq = lambda_liquid_ / (rho_l * cp_liquid_);
  double time_tmp = time_integral / 2;
  double left_time = 0.;
  double right_time = time_integral;
  double temperature_integral_eval = 1.e20;
  // const int sign_temperature_integral = signbit(temperature_integral);
  auto fflambda = [](const double& r, const double& R, const double& alpha, const double& t)
  { return R / r * (1 - erf( (r - R)/(2 * sqrt(alpha * t)))) ; };
  while(abs(temperature_integral_eval - temperature_integral) > 1e-6)
    {
      temperature_integral_eval = 0.;
      time_tmp = (right_time + left_time) / 2;

      const int max_unknowns = 100;
      const double radial_incr = Delta_R / (max_unknowns - 1);
      for (int i=0; i<max_unknowns-1; i++)
        temperature_integral_eval += (fflambda(R0 + radial_incr * i, R0, alpha_liq, time_tmp)
                                      + fflambda(R0 + radial_incr * (i+1), R0, alpha_liq, time_tmp))
                                     * (radial_incr / 2);
      temperature_integral_eval = T1 - T1 * temperature_integral_eval / Delta_R;
      // temperature_integral_eval is positive
      if (temperature_integral_eval > temperature_integral)
        right_time = time_tmp;
      else
        left_time = time_tmp;
    }
  auto flambda = [](const double& r, const double& R, const double& alpha, const double& t, const double& Tinfty)
  { return Tinfty - Tinfty * R / r * (1- erf( (r - R)/(2 * sqrt(alpha * t)))) ; };
  auto glambda = [](const double& r, const double& R, const double& alpha, const double& t, const double& Tinfty)
  { return Tinfty * R / pow(r,2) * (1- erf( (r - R)/(2 * sqrt(alpha * t)))) ; };
  auto hlambda = [](const double& r, const double& R, const double& alpha, const double& t, const double& Tinfty)
  { return Tinfty * R / r / (2 * sqrt(alpha * t)) * (2 / sqrt(M_PI)) * (exp(-pow((r - R)/(2 * sqrt(alpha * t)),2))) ; };

  const double temperature_end = flambda(R1, R0, alpha_liq, time_tmp, T1);
  temperature_end_prev = temperature_end;
  const double temperature_derivative_end = glambda(R1, R0, alpha_liq, time_tmp, T1)
                                            + hlambda(R1, R0, alpha_liq, time_tmp, T1);
  if (debug_)
    {
      Cerr << "Temperature at the probes end: " << temperature_end_prev << finl;
      Cerr << "Temperature at the probes end: " << temperature_end << finl;
      Cerr << "Temperature derivative at the probes end: " << temperature_derivative_end << finl;
    }
  return time_tmp;
}

void IJK_Thermal_Subresolution::compute_Nusselt_spherical_diffusion()
{
  const double T1 = delta_T_subcooled_overheated_;
  const double R0 = single_centred_bubble_radius_ini_;
  const double rho_l = ref_ijk_ft_->get_rho_l();
  const double alpha_liq = lambda_liquid_ / (rho_l * cp_liquid_);
  auto glambda = [](const double& r, const double& R, const double& alpha, const double& t, const double& Tinfty)
  { return Tinfty * R / pow(r,2) * (1- erf( (r - R)/(2 * sqrt(alpha * (t + 1e-16))))) ; };
  auto hlambda = [](const double& r, const double& R, const double& alpha, const double& t, const double& Tinfty)
  { return Tinfty * R / r / (2 * sqrt(alpha * (t + 1e-16))) * (2 / sqrt(M_PI)) * (exp(-pow((r - R)/(2 * sqrt(alpha * (t + 1e-16))),2))) ; };
  const double temperature_derivative_interface = glambda(R0, R0, alpha_liq, ref_ijk_ft_->get_current_time(), T1)
                                                  + hlambda(R0, R0, alpha_liq, ref_ijk_ft_->get_current_time(), T1);
  heat_flux_spherical_ = lambda_liquid_ * temperature_derivative_interface;
  nusselt_spherical_diffusion_ = abs(temperature_derivative_interface * (2 * single_centred_bubble_radius_ini_) / T1);
  const double temperature_derivative_interface_liquid = glambda(R0, R0, alpha_liq, ref_ijk_ft_->get_current_time(), mean_liquid_temperature_)
                                                         + hlambda(R0, R0, alpha_liq, ref_ijk_ft_->get_current_time(), mean_liquid_temperature_);
  nusselt_spherical_diffusion_liquid_ = abs(temperature_derivative_interface_liquid * (2 * single_centred_bubble_radius_ini_)
                                            / mean_liquid_temperature_);
}


double IJK_Thermal_Subresolution::get_time_inflection_derivative(double& temperature_end_min)
{
  const double T1 = delta_T_subcooled_overheated_;
  const double R0 = single_centred_bubble_radius_ini_;
  const double R1 = get_probes_length() + R0;
  const double rho_l = ref_ijk_ft_->get_rho_l();
  const double alpha_liq = lambda_liquid_ / (rho_l * cp_liquid_);
  // Solve dphidt(dphidr(T)) == 0
  const double inflection_time_temperature_derivative = 0.5 * (pow(R0,2) * R1 - 2 * R0 * pow(R1,2) + pow(R1,3)) / (R0*alpha_liq);
  auto glambda = [](const double& r, const double& R, const double& alpha, const double& t, const double& Tinfty)
  { return - (-Tinfty) * R / pow(r,2) * (1- erf( (r - R)/(2 * sqrt(alpha * t)))) ; };
  auto hlambda = [](const double& r, const double& R, const double& alpha, const double& t, const double& Tinfty)
  { return - (-Tinfty) * R / r / (2 * sqrt(alpha * t)) * (2 / sqrt(M_PI)) * (exp(-pow((r - R)/(2 * sqrt(alpha * t)),2))) ; };
  const double inflection_temperature_derivative = glambda(R1, R0, alpha_liq, inflection_time_temperature_derivative, T1)
                                                   + hlambda(R1, R0, alpha_liq, inflection_time_temperature_derivative, T1);
  temperature_end_min = (inflection_temperature_derivative * R1 * (R1 - R0)) / (R0);
  return inflection_time_temperature_derivative;
}

double IJK_Thermal_Subresolution::find_time_dichotomy_derivative(const double& temperature_derivative,
                                                                 double& temperature_limit_left,
                                                                 double& temperature_limit_right)
{
  const double time_max = 1000.;
  const double T1 = delta_T_subcooled_overheated_;
  const double R0 = single_centred_bubble_radius_ini_;
  const double R1 = get_probes_length() + R0;
  const double rho_l = ref_ijk_ft_->get_rho_l();
  const double alpha_liq = lambda_liquid_ / (rho_l * cp_liquid_);
  double temperature_end_min = 0.;
  // Solve dphidt(dphidr(T)) == 0
  const double inflection_time_temperature_derivative = get_time_inflection_derivative(temperature_end_min);

  double time_tmp = 0.;
  double left_time = inflection_time_temperature_derivative;
  double right_time = time_max;
  double temperature_derivative_eval = 1.e20;
  auto flambda = [](const double& r, const double& R, const double& alpha, const double& t, const double& Tinfty)
  { return Tinfty - Tinfty * R / r * (1- erf( (r - R)/(2 * sqrt(alpha * t)))) ; };
  auto glambda = [](const double& r, const double& R, const double& alpha, const double& t, const double& Tinfty)
  { return - (-Tinfty) * R / pow(r,2) * (1- erf( (r - R)/(2 * sqrt(alpha * t)))) ; };
  auto hlambda = [](const double& r, const double& R, const double& alpha, const double& t, const double& Tinfty)
  { return - (-Tinfty) * R / r / (2 * sqrt(alpha * t)) * (2 / sqrt(M_PI)) * (exp(-pow((r - R)/(2 * sqrt(alpha * t)),2))) ; };
  const double inflection_temperature_derivative = glambda(R1, R0, alpha_liq, inflection_time_temperature_derivative, T1)
                                                   + hlambda(R1, R0, alpha_liq, inflection_time_temperature_derivative, T1);
  const double inflection_temperature = flambda(R1, R0, alpha_liq, inflection_time_temperature_derivative, T1);
  if (debug_)
    {
      Cerr << "Inflection temperature derivative init: " << inflection_temperature_derivative << finl;
      Cerr << "Temperature at the inflection point: " << inflection_temperature << finl;
      Cerr << "Temperature derivative: " << temperature_derivative << finl;
    }
  while(abs(temperature_derivative_eval - temperature_derivative) > 1e-6)
    {
      temperature_derivative_eval = 0.;
      time_tmp = (right_time + left_time) / 2;
      temperature_derivative_eval = glambda(R1, R0, alpha_liq, time_tmp, T1)
                                    + hlambda(R1, R0, alpha_liq, time_tmp, T1);
      if (abs(temperature_derivative_eval) > abs(inflection_temperature_derivative))
        {
          time_tmp = inflection_time_temperature_derivative;
          break;
        }
      if (abs(temperature_derivative_eval - temperature_derivative) > 1e-2 && abs(left_time-right_time) < 1e-2)
        right_time = 2 * right_time;
      // temperature_integral_eval is positive
      if (abs(temperature_derivative_eval) < abs(temperature_derivative))
        right_time = time_tmp;
      else
        left_time = time_tmp;
    }
  const double temperature_end = flambda(R1, R0, alpha_liq, time_tmp, T1);
  double temperature_middle = 0.5 * (temperature_limit_right + temperature_limit_left);
  if (debug_)
    {
      Cerr << "Temperature at the probes end: " << temperature_end << finl;
      Cerr << "Temperature left limit: " << temperature_limit_left << finl;
      Cerr << "Temperature right limit: " << temperature_limit_right << finl;
      Cerr << "Temperature middle: " << temperature_middle << finl;
    }
  if (abs(temperature_end) > abs(temperature_middle))
    temperature_limit_left = temperature_middle;
  else
    temperature_limit_right = temperature_middle;
  const double temperature_derivative_end = glambda(R1, R0, alpha_liq, time_tmp, T1)
                                            + hlambda(R1, R0, alpha_liq, time_tmp, T1);
  if (debug_)
    {
      Cerr << "Temperature left limit: " << temperature_limit_left << finl;
      Cerr << "Temperature right limit: " << temperature_limit_right << finl;
    }

  Cerr << "Temperature derivative at the probes end: " << temperature_derivative_end << finl;
  return time_tmp;
}

void IJK_Thermal_Subresolution::set_field_T_ana()
{
  if (liste_post_instantanes_.contient_("TEMPERATURE_ANA") || liste_post_instantanes_.contient_("ECART_T_ANA"))
    {
      if (spherical_diffusion_)
        {
          Nom expression_T_ana = compute_quasi_static_spherical_diffusion_expression(ref_ijk_ft_->get_current_time(), 0, 0);
          set_field_data(temperature_ana_, expression_T_ana);
          correct_any_temperature_field_for_visu(temperature_ana_);
          if (liste_post_instantanes_.contient_("ECART_T_ANA"))
            compare_temperature_fields(temperature_, temperature_ana_, ecart_t_ana_, ecart_t_ana_rel_);
          correct_any_temperature_fields_for_eulerian_fluxes(ecart_t_ana_);
          correct_any_temperature_fields_for_eulerian_fluxes(ecart_t_ana_rel_);
          evaluate_total_liquid_absolute_parameter(ecart_t_ana_, error_temperature_ana_total_);
          evaluate_total_liquid_parameter_squared(ecart_t_ana_, error_temperature_ana_squared_total_);
          evaluate_total_liquid_absolute_parameter(ecart_t_ana_rel_, error_temperature_ana_rel_total_);
          error_temperature_ana_total_ = abs(error_temperature_ana_total_ / delta_T_subcooled_overheated_);
          error_temperature_ana_squared_total_ = sqrt(error_temperature_ana_squared_total_) / abs(delta_T_subcooled_overheated_);
          error_temperature_ana_rel_total_ = abs(error_temperature_ana_rel_total_);
        }
    }
}

void IJK_Thermal_Subresolution::update_thermal_properties()
{
  IJK_Thermal_base::update_thermal_properties();
}

void IJK_Thermal_Subresolution::post_process_after_temperature_increment()
{
  IJK_Thermal_base::post_process_after_temperature_increment();
  if (debug_)
    Cerr << "Compute mean liquid temperature" << finl;
  compute_mean_liquid_temperature();
  if (debug_)
    Cerr << "Compute Nusselt spherical diffusion" << finl;
  compute_Nusselt_spherical_diffusion();
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
          div_coeff_grad_T_volume_(i,j,k) = ope / rhocp_l;
          d_temperature_(i,j,k) += resu;
          if (liste_post_instantanes_.contient_("DIV_LAMBDA_GRAD_T"))
            div_coeff_grad_T_(i,j,k) = resu;
        }
  if (debug_)
    Cerr << "Uniform lambda: " << temperature_diffusion_op_.get_uniform_lambda() << finl;
}

void IJK_Thermal_Subresolution::correct_any_temperature_fields_for_eulerian_fluxes(IJK_Field_double& temperature)
{
  const int ni = temperature.ni();
  const int nj = temperature.nj();
  const int nk = temperature.nk();
  for (int k = 0; k < nk; k++)
    for (int j = 0; j < nj; j++)
      for (int i = 0; i < ni; i++)
        {
          const double indic = ref_ijk_ft_->itfce().I(i,j,k);
          if (fabs(indic) < LIQUID_INDICATOR_TEST) // Mixed cells and pure vapour cells
            temperature(i,j,k) = 0.;
        }
  temperature.echange_espace_virtuel(temperature.ghost());
}

void IJK_Thermal_Subresolution::correct_temperature_for_eulerian_fluxes()
{
  if (override_vapour_mixed_values_)
    correct_any_temperature_fields_for_eulerian_fluxes(temperature_);
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

void IJK_Thermal_Subresolution::compare_temperature_fields(const IJK_Field_double& temperature,
                                                           const IJK_Field_double& temperature_ana,
                                                           IJK_Field_double& error_temperature_ana,
                                                           IJK_Field_double& error_temperature_ana_rel)
{
  const int ni = temperature.ni();
  const int nj = temperature.nj();
  const int nk = temperature.nk();
  error_temperature_ana.data() = 0.;
  for (int k = 0; k < nk; k++)
    for (int j = 0; j < nj; j++)
      for (int i = 0; i < ni; i++)
        {
          error_temperature_ana(i,j,k) = (temperature(i,j,k) - temperature_ana(i,j,k));
          error_temperature_ana_rel(i,j,k) = error_temperature_ana(i,j,k) / (temperature_ana(i,j,k) + 1.e-16) * 100.;
        }
}

void IJK_Thermal_Subresolution::evaluate_total_liquid_absolute_parameter(const IJK_Field_double& field,
                                                                         double& total_parameter)
{
  const int ni = field.ni();
  const int nj = field.nj();
  const int nk = field.nk();
  total_parameter = 0;
  double liquid_volume = 0.;
  for (int k = 0; k < nk; k++)
    for (int j = 0; j < nj; j++)
      for (int i = 0; i < ni; i++)
        {
          const double indic = ref_ijk_ft_->itfce().I(i,j,k);
          liquid_volume += (vol_ * indic);
          total_parameter += abs(field(i,j,k)) * (vol_ * indic);
        }
  total_parameter = mp_sum(total_parameter);
  liquid_volume = mp_sum(liquid_volume);
  total_parameter = total_parameter / liquid_volume;
}

void IJK_Thermal_Subresolution::evaluate_total_liquid_parameter_squared(const IJK_Field_double& field,
                                                                        double& total_parameter)
{
  const int ni = field.ni();
  const int nj = field.nj();
  const int nk = field.nk();
  total_parameter = 0;
  double liquid_volume = 0.;
  for (int k = 0; k < nk; k++)
    for (int j = 0; j < nj; j++)
      for (int i = 0; i < ni; i++)
        {
          const double indic = ref_ijk_ft_->itfce().I(i,j,k);
          liquid_volume += (vol_ * indic);
          total_parameter += pow(field(i,j,k) * (vol_ * indic), 2);
        }
  total_parameter = mp_sum(total_parameter);
  liquid_volume = mp_sum(liquid_volume);
  total_parameter = total_parameter / pow(liquid_volume, 2);
}

void IJK_Thermal_Subresolution::correct_any_temperature_field_for_visu(IJK_Field_double& temperature)
{
  const int ni = temperature.ni();
  const int nj = temperature.nj();
  const int nk = temperature.nk();
  for (int k = 0; k < nk; k++)
    for (int j = 0; j < nj; j++)
      for (int i = 0; i < ni; i++)
        {
          // const double temperature = temperature_(i,j,k);
          const double indic = ref_ijk_ft_->itfce().I(i,j,k);
          // if (temperature > 0)
          if (indic < VAPOUR_INDICATOR_TEST)
            temperature(i,j,k) = 0;
        }
  temperature.echange_espace_virtuel(temperature.ghost());
}

void IJK_Thermal_Subresolution::correct_temperature_for_visu()
{
  /*
   * Correct only if the temperature gradient is post-processed
   * If not we may want to reconstruct the gradient field in Python
   * using the ghost temperature !
   */
  if (liste_post_instantanes_.contient_("GRAD_T_ELEM") && allow_temperature_correction_for_visu_)
    correct_any_temperature_field_for_visu(temperature_);
  if (liste_post_instantanes_.contient_("U_T_CONVECTIVE") && allow_temperature_correction_for_visu_)
    correct_any_temperature_field_for_visu(u_T_convective_);
  if (liste_post_instantanes_.contient_("U_T_CONVECTIVE_VOLUME") && allow_temperature_correction_for_visu_)
    correct_any_temperature_field_for_visu(u_T_convective_volume_);
  if (liste_post_instantanes_.contient_("DIV_LAMBDA_GRAD_T") && allow_temperature_correction_for_visu_)
    correct_any_temperature_field_for_visu(div_coeff_grad_T_);
  if (liste_post_instantanes_.contient_("DIV_LAMBDA_GRAD_T_VOLUME") && allow_temperature_correction_for_visu_)
    correct_any_temperature_field_for_visu(div_coeff_grad_T_volume_);
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

void IJK_Thermal_Subresolution::clip_max_temperature_values()
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
              const double indic = ref_ijk_ft_->itfce().I(i,j,k);
              if (temperature > 0 && indic > LIQUID_INDICATOR_TEST)
                temperature_(i,j,k) = 0;
            }
      temperature_.echange_espace_virtuel(temperature_.ghost());
    }
}

void IJK_Thermal_Subresolution::compute_mean_liquid_temperature()
{
  double vol_liq = 0.;
  mean_liquid_temperature_ = 0.;
  const int ni = temperature_.ni();
  const int nj = temperature_.nj();
  const int nk = temperature_.nk();
  for (int k = 0; k < nk; k++)
    for (int j = 0; j < nj; j++)
      for (int i = 0; i < ni; i++)
        {
          const double indic = ref_ijk_ft_->itfce().I(i,j,k);
          if (indic > VAPOUR_INDICATOR_TEST)
            {
              vol_liq += (indic * vol_);
              mean_liquid_temperature_ += (temperature_(i,j,k) * indic * vol_);
            }
        }
  vol_liq = Process::mp_sum(vol_liq);
  mean_liquid_temperature_ /= vol_liq;
  mean_liquid_temperature_ = Process::mp_sum(mean_liquid_temperature_);
}

void IJK_Thermal_Subresolution::compute_thermal_subproblems()
{
  is_first_time_step_ = (!ref_ijk_ft_->get_reprise()) && (ref_ijk_ft_->get_tstep()==0);
  first_step_thermals_post_ = (first_step_thermals_post_ && !ref_ijk_ft_->get_tstep());
  if (is_first_time_step_)
    first_time_step_temporal_ = first_time_step_temporal_ && is_first_time_step_;

  if (debug_probe_collision_)
    probe_collision_debug_field_.data() = 0.;

  if (debug_)
    debug_LRS_cells_.data() = -1.;

  if (debug_)
    Cerr << "Initialise thermal subproblems" << finl;
  initialise_thermal_subproblems();

  detect_probe_collision();

  int varying_probes_length = 1;
  const int varying_probes = first_time_step_varying_probes_;
  if (!enable_resize_probe_collision_)
    if (readjust_probe_length_from_vertices_)
      probe_variations_enabled_ = 1;
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

  pre_initialise_thermal_subproblems_any_matrices();

  reset_subresolution_distributed_vectors();

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
    Cerr << "Solve thermal subproblems" << finl;
  solve_thermal_subproblems();

  if (debug_)
    Cerr << "Compute material derivative (modelling)" << finl;
  approximate_temperature_increment_material_derivative();

  if (debug_)
    Cerr << "Store probe temperature" << finl;
  store_previous_temperature_indicator_velocities();

  if (debug_)
    Cerr << "Prepare thermal flux correction" << finl;
  update_intersections();
  prepare_thermal_flux_correction();
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

double IJK_Thermal_Subresolution::get_probes_length()
{
  return coeff_distance_diagonal_ * cell_diagonal_;
}

void IJK_Thermal_Subresolution::compute_overall_probes_parameters()
{
  if (!disable_subresolution_ || reference_gfm_on_probes_)
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
      if (first_time_step_temporal_ || implicit_solver_from_previous_probe_field_)
        {
          int check_nb_elem;
          check_nb_elem = finite_difference_assembler_.build(identity_matrix_explicit_implicit_, points_per_thermal_subproblem_, -1);
          if (debug_)
            Cerr << "Check_nb_elem: " << check_nb_elem << finl;
        }

      if (!use_sparse_matrix_)
        {
          thermal_subproblems_matrix_assembly_for_solver_.typer("Matrice_Morse");
          thermal_subproblems_matrix_assembly_for_solver_reduced_.typer("Matrice_Morse");
        }
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
 * TODO : Should be moved to a class of tools ?
 */
void IJK_Thermal_Subresolution::compute_first_order_operator_raw(Matrice& radial_first_order_operator)
{
  /*
   * Compute the first-order matrix for Finite-Differences
   */
  // TODO: Replace with matrice morse
  int check_nb_elem;
  check_nb_elem = finite_difference_assembler_.build(radial_first_order_operator, points_per_thermal_subproblem_, 0);
  if (debug_)
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
  // TODO: Replace with matrice morse
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
                                                                        pre_factor_subproblems_number_,
                                                                        remove_append_subproblems_,
                                                                        use_sparse_matrix_);
}

void IJK_Thermal_Subresolution::initialise_thermal_subproblems()
{
  if (!disable_subresolution_ || reference_gfm_on_probes_)
    {
      const IJK_Field_double& indicator = ref_ijk_ft_->itfce().I();
      const int ni = temperature_.ni();
      const int nj = temperature_.nj();
      const int nk = temperature_.nk();
      int counter = 0;
      for (int k = 0; k < nk; k++)
        for (int j = 0; j < nj; j++)
          for (int i = 0; i < ni; i++)
            if (fabs(indicator(i,j,k)) > VAPOUR_INDICATOR_TEST && fabs(indicator(i,j,k)) < LIQUID_INDICATOR_TEST)
              {
                if (debug_)
                  {
                    Cerr << "Liquid Indicator (Old): " << indicator(i,j,k) << finl;
                    Cerr << "Liquid Indicator (Next): " << ref_ijk_ft_->itfce().In()(i,j,k) << finl;
                    debug_LRS_cells_(i,j,k) = indicator(i,j,k);
                  }
                thermal_local_subproblems_.associate_sub_problem_to_inputs((*this),
                                                                           i,j,k,
                                                                           indicator(i,j,k),
                                                                           ref_ijk_ft_->get_timestep(),
                                                                           ref_ijk_ft_->get_current_time(),
                                                                           ref_ijk_ft_->itfce(),
                                                                           ref_ijk_ft_->get_velocity(),
                                                                           ref_ijk_ft_->get_velocity_ft(),
                                                                           ref_ijk_ft_->get_pressure_ghost_cells());
                counter++;
              }
      thermal_local_subproblems_.set_effective_subproblems(enable_probe_collision_detection_);
      thermal_local_subproblems_.compute_global_indices();
      if (!counter)
        thermal_local_subproblems_.associate_variables_for_post_processing((*this));
    }
}

void IJK_Thermal_Subresolution::detect_probe_collision()
{
  if (enable_probe_collision_detection_)
    {
      interpolate_indicator_on_probes();
      reajust_probes_length_collisions();
      clear_sort_problems_colliding_bubbles();
    }
}

void IJK_Thermal_Subresolution::pre_initialise_thermal_subproblems_any_matrices()
{
  pre_initialise_thermal_subproblems_matrices();
}

void IJK_Thermal_Subresolution::pre_initialise_thermal_subproblems_matrices()
{
  if (ref_ijk_ft_->get_tstep()==0 && pre_initialise_thermal_subproblems_list_)
    {
      if (!disable_subresolution_)
        {
          // int nb_subproblems_ini = thermal_local_subproblems_.get_subproblems_counter();
          int nb_subproblems_ini = thermal_local_subproblems_.get_effective_subproblems_counter();
          if (!(ref_ijk_ft_->get_disable_convection_qdm() && ref_ijk_ft_->get_disable_diffusion_qdm()))
            nb_subproblems_ini = Process::mp_sum(nb_subproblems_ini);
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
            finite_difference_assembler_.complete_empty_matrices_initialisation(identity_matrix_subproblems_,
                                                                                radial_second_order_operator_raw_,
                                                                                nb_subproblems_ini, max_subproblems_predicted);
        }
    }
}

void IJK_Thermal_Subresolution::reset_subresolution_distributed_vectors()
{
  thermal_subproblems_rhs_assembly_.reset();
  thermal_subproblems_temperature_solution_.detach_vect();
  if (first_time_step_temporal_ && first_time_step_explicit_)
    thermal_subproblems_temperature_solution_ini_.detach_vect();

  // For a constant number of points (we can just change the probe size)
  // const int nb_points = points_per_thermal_subproblem_ * thermal_local_subproblems_.get_subproblems_counter();
  const int nb_points = points_per_thermal_subproblem_ * thermal_local_subproblems_.get_effective_subproblems_counter();
  thermal_subproblems_rhs_assembly_.resize(nb_points);
  thermal_subproblems_temperature_solution_.resize(nb_points);
  if (first_time_step_temporal_ && first_time_step_explicit_)
    thermal_subproblems_temperature_solution_ini_.resize(nb_points);
}

void IJK_Thermal_Subresolution::interpolate_indicator_on_probes()
{
  if (!disable_subresolution_ || reference_gfm_on_probes_)
    thermal_local_subproblems_.interpolate_indicator_on_probes();
}

void IJK_Thermal_Subresolution::clear_sort_problems_colliding_bubbles()
{
  if (!disable_subresolution_) // || reference_gfm_on_probes_)
    thermal_local_subproblems_.clear_sort_problems_colliding_bubbles();
}

void IJK_Thermal_Subresolution::interpolate_project_velocities_on_probes()
{
  if (!disable_subresolution_ || reference_gfm_on_probes_)
    thermal_local_subproblems_.interpolate_project_velocities_on_probes();
}

void IJK_Thermal_Subresolution::reajust_probes_length_collisions()
{
  if (!disable_subresolution_)
    {
      thermal_local_subproblems_.reajust_probes_length(0);
      if (!readjust_probe_length_from_vertices_)
        thermal_local_subproblems_.compute_modified_probe_length(probe_variations_enabled_);
    }
}

void IJK_Thermal_Subresolution::reajust_probes_length()
{
  if (!disable_subresolution_)
    if (readjust_probe_length_from_vertices_ ||
        first_time_step_varying_probes_)
      {
        thermal_local_subproblems_.reajust_probes_length();
        if (first_time_step_varying_probes_)
          {
            probe_variations_enabled_ = thermal_local_subproblems_.get_probe_variations_enabled(probe_variations_priority_);
            first_time_step_varying_probes_ = probe_variations_enabled_;
          }
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
          initialise_sparse_indices_ = 1;
          if (first_time_step_temporal_ || implicit_solver_from_previous_probe_field_)
            {
              if (use_sparse_matrix_)
                initialise_identity_matrices_sparse(identity_matrix_explicit_implicit_, identity_matrix_subproblems_);
              else
                initialise_identity_matrices(identity_matrix_explicit_implicit_, identity_matrix_subproblems_);
            }
          if (use_sparse_matrix_)
            {
              initialise_radial_convection_operator_sparse(radial_first_order_operator_, radial_convection_matrix_);
              initialise_radial_diffusion_operator_sparse(radial_second_order_operator_, radial_diffusion_matrix_);
            }
          else
            {
              initialise_radial_convection_operator(radial_first_order_operator_, radial_convection_matrix_);
              initialise_radial_diffusion_operator(radial_second_order_operator_, radial_diffusion_matrix_);
            }
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
    Cerr << "Initialise the Identity matrix" << finl;
  // const int nb_subproblems = thermal_local_subproblems_.get_subproblems_counter();
  const int nb_subproblems = thermal_local_subproblems_.get_effective_subproblems_counter();
  finite_difference_assembler_.initialise_matrix_subproblems(identity_matrix_subproblems,
                                                             identity_matrix,
                                                             nb_subproblems,
                                                             first_time_step_varying_probes_);
  if (debug_)
    Cerr << "Identity matrix has been initialised." << finl;
}

void IJK_Thermal_Subresolution::initialise_identity_matrices_sparse(Matrice& identity_matrix,
                                                                    Matrice& identity_matrix_subproblems)
{
  /*
   * Compute the convection matrices for Finite-Differences
   */
  const int first_initialisation = (ref_ijk_ft_->get_tstep() == 0);
  if (debug_)
    Cerr << "Initialise the Identity sparse matrix" << finl;
  // const int nb_subproblems = thermal_local_subproblems_.get_subproblems_counter();
  const int nb_subproblems = thermal_local_subproblems_.get_effective_subproblems_counter();
  finite_difference_assembler_.initialise_sparse_matrix_subproblems(identity_matrix_subproblems,
                                                                    identity_matrix,
                                                                    nb_subproblems,
                                                                    first_time_step_varying_probes_,
                                                                    first_indices_sparse_matrix_,
                                                                    first_initialisation,
                                                                    initialise_sparse_indices_);
  if (debug_)
    Cerr << "Identity sparse matrix has been initialised." << finl;
}

void IJK_Thermal_Subresolution::initialise_radial_convection_operator(Matrice& radial_first_order_operator,
                                                                      Matrice& radial_convection_matrix)
{
  /*
   * Compute the convection matrices for Finite-Differences
   */
  if (debug_)
    Cerr << "Initialise the Radial convection operator" << finl;
  // const int nb_subproblems = thermal_local_subproblems_.get_subproblems_counter();
  const int nb_subproblems = thermal_local_subproblems_.get_effective_subproblems_counter();
  finite_difference_assembler_.initialise_matrix_subproblems(radial_convection_matrix,
                                                             radial_first_order_operator,
                                                             nb_subproblems,
                                                             first_time_step_varying_probes_);
  if (debug_)
    Cerr << "Radial convection operator has been initialised." << finl;
}

void IJK_Thermal_Subresolution::initialise_radial_convection_operator_sparse(Matrice& radial_first_order_operator,
                                                                             Matrice& radial_convection_matrix)
{
  /*
   * Compute the convection matrices for Finite-Differences
   */
  const int first_initialisation = (ref_ijk_ft_->get_tstep() == 0);
  if (debug_)
    Cerr << "Initialise the Radial convection sparse operator" << finl;
  // const int nb_subproblems = thermal_local_subproblems_.get_subproblems_counter();
  const int nb_subproblems = thermal_local_subproblems_.get_effective_subproblems_counter();
  finite_difference_assembler_.initialise_sparse_matrix_subproblems(radial_convection_matrix,
                                                                    radial_first_order_operator,
                                                                    nb_subproblems,
                                                                    first_time_step_varying_probes_,
                                                                    first_indices_sparse_matrix_,
                                                                    first_initialisation,
                                                                    initialise_sparse_indices_);
  if (debug_)
    Cerr << "Radial convection sparse operator has been initialised." << finl;
}


void IJK_Thermal_Subresolution::initialise_radial_diffusion_operator(Matrice& radial_second_order_operator,
                                                                     Matrice& radial_diffusion_matrix)
{
  /*
   * Compute the diffusion matrices for Finite-Differences
   */
  if (debug_)
    Cerr << "Initialise the Radial diffusion operator" << finl;
  // const int nb_subproblems = thermal_local_subproblems_.get_subproblems_counter();
  const int nb_subproblems = thermal_local_subproblems_.get_effective_subproblems_counter();
  finite_difference_assembler_.initialise_matrix_subproblems(radial_diffusion_matrix, radial_second_order_operator, nb_subproblems, first_time_step_varying_probes_);
  if (debug_)
    Cerr << "Radial diffusion operator has been initialised." << finl;
}

void IJK_Thermal_Subresolution::initialise_radial_diffusion_operator_sparse(Matrice& radial_second_order_operator,
                                                                            Matrice& radial_diffusion_matrix)
{
  /*
   * Compute the diffusion matrices for Finite-Differences
   */
  const int first_initialisation = (ref_ijk_ft_->get_tstep() == 0);
  if (debug_)
    Cerr << "Initialise the Radial diffusion sparse operator" << finl;
  // const int nb_subproblems = thermal_local_subproblems_.get_subproblems_counter();
  const int nb_subproblems = thermal_local_subproblems_.get_effective_subproblems_counter();
  finite_difference_assembler_.initialise_sparse_matrix_subproblems(radial_diffusion_matrix,
                                                                    radial_second_order_operator,
                                                                    nb_subproblems,
                                                                    first_time_step_varying_probes_,
                                                                    first_indices_sparse_matrix_,
                                                                    first_initialisation,
                                                                    initialise_sparse_indices_);
  initialise_sparse_indices_ = 1;
  if (debug_)
    Cerr << "Radial diffusion sparse operator has been initialised." << finl;
}

void IJK_Thermal_Subresolution::compute_local_substep()
{
  if (!disable_subresolution_)
    if (first_time_step_temporal_)
      {
        const double local_time_step = thermal_local_subproblems_.get_min_euler_time_step(nb_iter_explicit_local_);
        if (debug_)
          {
            Cerr << "Local timestep: " << local_time_step << finl;
            Cerr << "Number of sub-steps: " << nb_iter_explicit_local_ << finl;
          }
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
      }
}

void IJK_Thermal_Subresolution::prepare_temporal_schemes()
{
  if (!disable_subresolution_)
    {
      if (first_time_step_temporal_ || implicit_solver_from_previous_probe_field_)
        thermal_local_subproblems_.prepare_temporal_schemes();
      thermal_subproblems_matrix_assembly_ = radial_convection_matrix_;
      finite_difference_assembler_.sum_any_matrices_subproblems(thermal_subproblems_matrix_assembly_, radial_diffusion_matrix_, use_sparse_matrix_, debug_);
      if (first_time_step_temporal_ || implicit_solver_from_previous_probe_field_)
        finite_difference_assembler_.sum_any_matrices_subproblems(thermal_subproblems_matrix_assembly_, identity_matrix_subproblems_, use_sparse_matrix_, debug_);
    }
}

void IJK_Thermal_Subresolution::compute_source_terms_impose_subresolution_boundary_conditions()
{
  if (!disable_subresolution_)
    {
      thermal_local_subproblems_.compute_source_terms_impose_boundary_conditions(boundary_condition_interface_,
                                                                                 interfacial_boundary_condition_value_,
                                                                                 impose_boundary_condition_interface_from_simulation_,
                                                                                 boundary_condition_end_,
                                                                                 end_boundary_condition_value_,
                                                                                 impose_user_boundary_condition_end_value_);
    }
}

void IJK_Thermal_Subresolution::approximate_temperature_increment_material_derivative()
{
  if (!disable_subresolution_ || reference_gfm_on_probes_)
    thermal_local_subproblems_.approximate_temperature_increment_material_derivative();
}

void IJK_Thermal_Subresolution::solve_thermal_subproblems()
{
  if (!disable_subresolution_)
    {
      check_wrong_values_rhs();
      convert_into_sparse_matrix();

      if (first_time_step_temporal_ && first_time_step_explicit_)
        {
          Cerr << "Finite-difference Explicit thermal sub-resolution !" << finl;
          DoubleVect thermal_subproblems_temperature_solution_tmp = thermal_subproblems_temperature_solution_ini_;
          for (int i=0; i<nb_iter_explicit_local_; i++)
            {
              thermal_subproblems_temperature_solution_ = (*thermal_subproblems_matrix_assembly_for_solver_ref_) * thermal_subproblems_temperature_solution_tmp;
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
              one_dimensional_advection_diffusion_thermal_solver_implicit_.resoudre_systeme((*thermal_subproblems_matrix_assembly_for_solver_ref_).valeur(),
                                                                                            thermal_subproblems_rhs_assembly_,
                                                                                            thermal_subproblems_temperature_solution_);
              Cerr << "Finite-difference thermal sub-resolution has finished (Implicit)!" << finl;
            }
          else
            {
              Cerr << "Finite-difference thermal sub-resolution has started !" << finl;
              one_dimensional_advection_diffusion_thermal_solver_.resoudre_systeme((*thermal_subproblems_matrix_assembly_for_solver_ref_).valeur(),
                                                                                   thermal_subproblems_rhs_assembly_,
                                                                                   thermal_subproblems_temperature_solution_);
              Cerr << "Finite-difference thermal sub-resolution has finished !" << finl;
            }
        }
    }
  retrieve_temperature_solution();
}

void IJK_Thermal_Subresolution::convert_into_sparse_matrix()
{
  /*
   * Convert into a huge sparse matrix
   */
  // const int nb_points = points_per_thermal_subproblem_ * thermal_local_subproblems_.get_subproblems_counter();
  const int nb_points = points_per_thermal_subproblem_ * thermal_local_subproblems_.get_effective_subproblems_counter();
  if (!use_sparse_matrix_)
    {
      Matrice_Morse& sparse_matrix_solver  = ref_cast(Matrice_Morse, thermal_subproblems_matrix_assembly_for_solver_.valeur());
      Matrice_Bloc& bloc_matrix_solver = ref_cast(Matrice_Bloc, thermal_subproblems_matrix_assembly_.valeur());
      bloc_matrix_solver.block_to_morse(sparse_matrix_solver);
      if (pre_initialise_thermal_subproblems_list_)
        {
          Matrice_Morse& sparse_matrix_solver_reduced  = ref_cast(Matrice_Morse, thermal_subproblems_matrix_assembly_for_solver_reduced_.valeur());
          finite_difference_assembler_.reduce_solver_matrix(sparse_matrix_solver,
                                                            sparse_matrix_solver_reduced,
                                                            nb_points,
                                                            pre_initialise_thermal_subproblems_list_);
          thermal_subproblems_matrix_assembly_for_solver_ref_ = &thermal_subproblems_matrix_assembly_for_solver_reduced_;
        }
      else
        thermal_subproblems_matrix_assembly_for_solver_ref_ = &thermal_subproblems_matrix_assembly_for_solver_;

    }
  else
    {
      thermal_subproblems_matrix_assembly_for_solver_ref_ = &thermal_subproblems_matrix_assembly_;
    }
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

  const int subproblems_vector_size = thermal_subproblems_rhs_assembly_.size();
  MD_Vector_std md_std = MD_Vector_std(subproblems_vector_size, subproblems_vector_size,
                                       pe_voisins_dummy, items_to_send_dummy,
                                       items_to_recv_dummy, blocs_to_recv_dummy);
  md_.copy(md_std);
  MD_Vector_tools::creer_tableau_distribue(md_, thermal_subproblems_rhs_assembly_); //, RESIZE_OPTIONS::NOCOPY_NOINIT);
  MD_Vector_tools::creer_tableau_distribue(md_, thermal_subproblems_temperature_solution_); //, RESIZE_OPTIONS::NOCOPY_NOINIT);
}

void IJK_Thermal_Subresolution::retrieve_temperature_solution()
{
  /*
   * Retrieve the overall vector
   */
  if (!disable_subresolution_ || reference_gfm_on_probes_)
    thermal_local_subproblems_.retrieve_radial_quantities();
}

void IJK_Thermal_Subresolution::store_previous_temperature_indicator_velocities()
{
  if (!disable_subresolution_ || reference_gfm_on_probes_)
    if (reconstruct_previous_probe_field_)
      thermal_local_subproblems_.store_previous_temperature_indicator_velocities();
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

void IJK_Thermal_Subresolution::compute_convective_diffusive_fluxes_face_centre()
{

  if (!conv_temperature_negligible_)
    compute_convective_fluxes_face_centre();
  if (!diff_temperature_negligible_)
    compute_diffusive_fluxes_face_centre();

  corrige_flux_->initialise_cell_neighbours_indices_to_correct();
  corrige_flux_->compute_cell_neighbours_faces_indices_for_spherical_correction(n_iter_distance_);

  compute_min_max_reachable_fluxes();
}

void IJK_Thermal_Subresolution::compute_convective_fluxes_face_centre()
{
  if (!disable_subresolution_ && convective_flux_correction_ && !use_reachable_fluxes_)
    corrige_flux_->compute_thermal_convective_fluxes();
}

void IJK_Thermal_Subresolution::compute_diffusive_fluxes_face_centre()
{
  if (!disable_subresolution_ && diffusive_flux_correction_ && !use_reachable_fluxes_)
    corrige_flux_->compute_thermal_diffusive_fluxes();
}

void IJK_Thermal_Subresolution::compute_min_max_reachable_fluxes()
{
  corrige_flux_->compute_cell_neighbours_faces_indices_to_correct(cell_faces_neighbours_corrected_all_bool_,
                                                                  cell_faces_neighbours_corrected_convective_,
                                                                  cell_faces_neighbours_corrected_diffusive_,
                                                                  neighbours_faces_weighting_colinearity_);
  /*
   * TODO: Should we put this in the jdd (just for test)
   */
  const int discontinous_min_max_neighbours_faces = 1;
  const int remove_external_neighbour_values = 1;
  const int check_neighbour_cell_centre = !remove_external_neighbour_values;
  corrige_flux_->compute_min_max_ijk_any_reachable_fluxes(cell_faces_neighbours_corrected_all_bool_,
                                                          neighbours_temperature_to_correct_,
                                                          cell_faces_neighbours_corrected_min_max_bool_,
                                                          discontinous_min_max_neighbours_faces,
                                                          check_neighbour_cell_centre,
                                                          remove_external_neighbour_values,
                                                          neighbours_temperature_to_correct_trimmed_);
  corrige_flux_->replace_cell_neighbours_thermal_convective_diffusive_fluxes_faces(cell_faces_neighbours_corrected_min_max_bool_,
                                                                                   cell_faces_neighbours_corrected_convective_,
                                                                                   0);
  corrige_flux_->replace_cell_neighbours_thermal_convective_diffusive_fluxes_faces(cell_faces_neighbours_corrected_min_max_bool_,
                                                                                   cell_faces_neighbours_corrected_diffusive_,
                                                                                   1);
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
  int correct_first_iter = (correct_temperature_cell_neighbours_first_iter_
                            && ref_ijk_ft_->get_tstep() == 1
                            && ref_ijk_ft_->get_reprise() == 0);
  if (debug_)
    Cerr << "Set correction cell neighbours" << finl;
  if (correct_first_iter_deactivate_cell_neighbours_ && correct_first_iter)
    {

      find_temperature_cell_neighbours_ = 0;
      use_temperature_cell_neighbours_ = 0;
      corrige_flux_->set_correction_cell_neighbours(find_temperature_cell_neighbours_,
                                                    neighbours_weighting_,
                                                    smooth_temperature_field_);
    }


  if (debug_)
    Cerr << "Compute temperature cell centre" << finl;
  corrige_flux_->compute_temperature_cell_centre(temperature_);

  if (debug_)
    Cerr << "Compute temperature cell centre neighbours" << finl;
  corrige_flux_->compute_temperature_cell_centre_neighbours(temperature_cell_neighbours_,
                                                            neighbours_temperature_to_correct_,
                                                            neighbours_temperature_colinearity_weighting_);
  if (debug_)
    {
      temperature_cell_neighbours_debug_.data() = 0.;
      corrige_flux_->replace_temperature_cell_centre_neighbours(temperature_cell_neighbours_debug_,
                                                                temperature_cell_neighbours_,
                                                                neighbours_temperature_to_correct_,
                                                                neighbours_temperature_colinearity_weighting_);
    }

  if (debug_)
    Cerr << "Replace temperature cell centres neighbours" << finl;
  replace_temperature_cell_centres_neighbours(0);
}

void IJK_Thermal_Subresolution::compute_temperature_cell_centres_second_correction()
{
  if (disable_mixed_cells_increment_)
    {
      if (debug_)
        Cerr << "Apply second temperature correction (Case C)" << finl;
      replace_temperature_cell_centres_neighbours(use_reachable_fluxes_ && !keep_first_reachable_fluxes_);
    }
}

void IJK_Thermal_Subresolution::replace_temperature_cell_centres_neighbours(const int& use_neighbours_temperature_to_correct_trimmed)
{
  int correct_first_iter = (correct_temperature_cell_neighbours_first_iter_
                            && ref_ijk_ft_->get_tstep() == 0
                            && ref_ijk_ft_->get_reprise() == 0);
  if (use_temperature_cell_neighbours_)
    {
      if (keep_first_reachable_fluxes_)
        {
          if (correct_first_iter)
            {
              corrige_flux_->replace_temperature_cell_centre_neighbours(temperature_,
                                                                        temperature_cell_neighbours_,
                                                                        neighbours_temperature_to_correct_,
                                                                        neighbours_temperature_colinearity_weighting_);
            }
          else
            corrige_flux_->compute_temperature_cell_centre(temperature_);
        }
      else
        {
          if (correct_temperature_cell_neighbours_first_iter_ == correct_first_iter)
            {
              if (use_neighbours_temperature_to_correct_trimmed)
                corrige_flux_->replace_temperature_cell_centre_neighbours(temperature_,
                                                                          temperature_cell_neighbours_,
                                                                          neighbours_temperature_to_correct_trimmed_,
                                                                          neighbours_temperature_colinearity_weighting_);

              else
                corrige_flux_->replace_temperature_cell_centre_neighbours(temperature_,
                                                                          temperature_cell_neighbours_,
                                                                          neighbours_temperature_to_correct_,
                                                                          neighbours_temperature_colinearity_weighting_);
            }
          else
            corrige_flux_->compute_temperature_cell_centre(temperature_);
        }
    }
  else
    corrige_flux_->compute_temperature_cell_centre(temperature_);
}


void IJK_Thermal_Subresolution::enforce_periodic_temperature_boundary_value()
{
  if (!disable_subresolution_ && enforce_periodic_boundary_value_)
    {
      const int ni = temperature_.ni();
      const int nj = temperature_.nj();
      const int nk = temperature_.nk();
      const int offset_i = temperature_.get_splitting().get_offset_local(0);
      const int offset_j = temperature_.get_splitting().get_offset_local(1);
      const int offset_k = temperature_.get_splitting().get_offset_local(2);
      const int ni_tot = temperature_.get_splitting().get_grid_geometry().get_nb_elem_tot(0);
      const int nj_tot = temperature_.get_splitting().get_grid_geometry().get_nb_elem_tot(1);
      const int nk_tot = temperature_.get_splitting().get_grid_geometry().get_nb_elem_tot(2);
      std::vector<double> indices_i_to_correct;
      std::vector<double> indices_j_to_correct;
      std::vector<double> indices_k_to_correct;
      int stencil_to_correct = stencil_periodic_boundary_value_;
      for (int i=0; i<stencil_to_correct; i++)
        {
          indices_i_to_correct.push_back(i);
          indices_j_to_correct.push_back(i);
          indices_k_to_correct.push_back(i);
          indices_i_to_correct.push_back(ni_tot-1-i);
          indices_j_to_correct.push_back(nj_tot-1-i);
          indices_k_to_correct.push_back(nk_tot-1-i);
        }
      int i,j,k;
      for (k = 0; k < nk; k++)
        if (std::find(indices_k_to_correct.begin(), indices_k_to_correct.end(), k+offset_k) != indices_k_to_correct.end())
          for (j = 0; j < nj; j++)
            for (i = 0; i < ni; i++)
              {
                const double indic = ref_ijk_ft_->itfce().I(i,j,k);
                if (fabs(indic)>LIQUID_INDICATOR_TEST) // Mixed cells and pure vapour cells
                  { temperature_(i,j,k) = delta_T_subcooled_overheated_; }
              }
      for (j = 0; j < nj; j++)
        if (std::find(indices_j_to_correct.begin(), indices_j_to_correct.end(), j+offset_j) != indices_j_to_correct.end())
          for (k = 0; k < nk; k++)
            for (i = 0; i < ni; i++)
              {
                const double indic = ref_ijk_ft_->itfce().I(i,j,k);
                if (fabs(indic)>LIQUID_INDICATOR_TEST) // Mixed cells and pure vapour cells
                  { temperature_(i,j,k) = delta_T_subcooled_overheated_; }
              }
      for (i = 0; i < ni; i++)
        if (std::find(indices_i_to_correct.begin(), indices_i_to_correct.end(), i+offset_i) != indices_i_to_correct.end())
          for (k = 0; k < nk; k++)
            for (j = 0; j < nj; j++)
              {
                const double indic = ref_ijk_ft_->itfce().I(i,j,k);
                if (fabs(indic)>LIQUID_INDICATOR_TEST) // Mixed cells and pure vapour cells
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
  if ((!disable_subresolution_) && (convective_flux_correction_ || diffusive_flux_correction_) && (!use_reachable_fluxes_))
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
  if (!disable_subresolution_)
    {
      if (disable_mixed_cells_increment_)
        corrige_flux_->set_zero_temperature_increment(d_temperature_);
    }
  else
    correct_temperature_increment_for_interface_leaving_cell();
}

void IJK_Thermal_Subresolution::clean_thermal_subproblems()
{
  if (ref_ijk_ft_->get_tstep() > 0)
    {
      if (!disable_subresolution_ || reference_gfm_on_probes_)
        thermal_local_subproblems_.clean();
      if (!disable_subresolution_)
        corrige_flux_->clear_vectors();
    }
}

void IJK_Thermal_Subresolution::clean_ijk_intersections()
{
  if (!disable_subresolution_)
    corrige_flux_->clean();
}

void IJK_Thermal_Subresolution::set_thermal_subresolution_outputs(const Nom& interfacial_quantities_thermal_probes,
                                                                  const Nom& overall_bubbles_quantities,
                                                                  const Nom& local_quantities_thermal_probes_time_index_folder)
{
  if (!disable_subresolution_ || reference_gfm_on_probes_)
    {
      if (debug_)
        Cerr << "Compute bubbles quantities" << finl;
      thermal_local_subproblems_.compute_overall_bubbles_quantities((*this));
      if (debug_)
        Cerr << "Sort spherical coords" << finl;
      thermal_local_subproblems_.sort_limited_probes_spherical_coords_post_processing(post_process_all_probes_,
                                                                                      nb_theta_post_pro_, nb_phi_post_pro_,
                                                                                      1, 1);
      if (debug_)
        Cerr << "Write post-processings" << finl;
      thermal_local_subproblems_.thermal_subresolution_outputs_parallel(rang_,
                                                                        interfacial_quantities_thermal_probes,
                                                                        overall_bubbles_quantities,
                                                                        local_quantities_thermal_probes_time_index_folder);
    }
}



