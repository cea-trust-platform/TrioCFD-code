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
// File      : IJK_One_Dimensional_Subproblem.cpp
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/Temperature
//
/////////////////////////////////////////////////////////////////////////////

#include <IJK_One_Dimensional_Subproblem.h>
#include <IJK_Navier_Stokes_tools.h>
#include <Ouvrir_fichier.h>
#include <IJK_FT.h>
#include <IJK_Thermal_base.h>
#include <IJK_Thermal_Subresolution.h>


Implemente_instanciable_sans_constructeur( IJK_One_Dimensional_Subproblem, "IJK_One_Dimensional_Subproblem", Objet_U ) ;

IJK_One_Dimensional_Subproblem::IJK_One_Dimensional_Subproblem()
{
  init_ = 1;
  debug_ = 0;

  points_per_thermal_subproblem_ = nullptr;
  points_per_thermal_subproblem_base_ = nullptr;
  alpha_ = nullptr;
  lambda_ = nullptr;
  coeff_distance_diagonal_ = nullptr;
  cell_diagonal_ = nullptr;

  interfaces_ = nullptr;
  eulerian_distance_ = nullptr;
  eulerian_curvature_ = nullptr;
  eulerian_interfacial_area_ = nullptr;
  eulerian_normal_vect_ = nullptr;
  eulerian_facets_barycentre_ = nullptr;

  temperature_ = nullptr;
  temperature_ft_ = nullptr;
  temperature_before_extrapolation_ = nullptr;
  velocity_ = nullptr;
  velocity_ft_ = nullptr;
  pressure_ = nullptr;

  grad_T_elem_ = nullptr;
  hess_diag_T_elem_ = nullptr;
  hess_cross_T_elem_ = nullptr;

  radial_coordinates_ = nullptr;
  radial_coordinates_base_ = nullptr;

  radial_first_order_operator_raw_base_ = nullptr;
  radial_second_order_operator_raw_base_ = nullptr;
  radial_first_order_operator_base_ = nullptr;
  radial_second_order_operator_base_ = nullptr;
  radial_diffusion_matrix_base_ = nullptr;
  radial_convection_matrix_base_ = nullptr;
  radial_velocity_convection_matrix_base_=nullptr;

  radial_first_order_operator_=nullptr;
  radial_second_order_operator_=nullptr;
  finite_difference_assembler_=nullptr;
  thermal_subproblems_matrix_assembly_=nullptr;
  thermal_subproblems_rhs_assembly_=nullptr;
  thermal_subproblems_temperature_solution_=nullptr;
  thermal_subproblems_temperature_solution_ini_=nullptr;

  interfacial_boundary_condition_value_ = 0.;
  end_boundary_condition_value_ = 0.;
  start_index_ = 0;
  end_index_ = 0;

  correct_radial_velocity_ = 1;
  correct_tangential_temperature_gradient_ = 0;
  correct_tangential_temperature_hessian_ = 0;

  first_tangential_vector_compo_solver_=nullptr;
  second_tangential_vector_compo_solver_=nullptr;

  first_tangential_velocity_solver_ = nullptr;
  second_tangential_velocity_solver_ = nullptr;

  tangential_temperature_gradient_first_solver_ = nullptr;
  tangential_temperature_gradient_second_solver_ = nullptr;

  max_u_ = 0.9 * INVALID_VELOCITY_CFL;
  identity_matrix_explicit_implicit_base_ = nullptr;
  identity_matrix_explicit_implicit_ = nullptr;
  identity_matrix_subproblems_ = nullptr;

  first_time_step_temporal_ = nullptr;
}

IJK_One_Dimensional_Subproblem::IJK_One_Dimensional_Subproblem(const IJK_FT_double& ijk_ft) : IJK_One_Dimensional_Subproblem()
{
  associer(ijk_ft);
}

Sortie& IJK_One_Dimensional_Subproblem::printOn( Sortie& os ) const
{
  Objet_U::printOn( os );
  return os;
}

Entree& IJK_One_Dimensional_Subproblem::readOn( Entree& is )
{
  Objet_U::readOn( is );
  return is;
}

void IJK_One_Dimensional_Subproblem::get_ijk_indices(int& i, int& j, int& k) const
{
  i = (int) index_i_;
  j = (int) index_j_;
  k = (int) index_k_;
}

void IJK_One_Dimensional_Subproblem::reinit_variable(DoubleVect& vect)
{
  vect.resize(*points_per_thermal_subproblem_);
  vect *= 0.;
}



/*
 * TODO: Remplacer cette methode pour eviter de fournir tous les attributs
 */
void IJK_One_Dimensional_Subproblem::associate_sub_problem_to_inputs(IJK_Thermal_Subresolution& ref_thermal_subresolution,
                                                                     int i, int j, int k,
                                                                     int init,
                                                                     int sub_problem_index,
                                                                     double global_time_step,
                                                                     double current_time,
                                                                     int compo_connex,
                                                                     double distance,
                                                                     double curvature,
                                                                     double interfacial_area,
                                                                     ArrOfDouble facet_barycentre,
                                                                     ArrOfDouble normal_vector,
                                                                     double bubble_rising_velocity,
                                                                     ArrOfDouble bubble_rising_vector,
                                                                     ArrOfDouble bubble_barycentre,
                                                                     const double& indicator,
                                                                     const IJK_Interfaces& interfaces,
                                                                     const FixedVector<IJK_Field_double, 3>& velocity,
                                                                     const FixedVector<IJK_Field_double, 3>& velocity_ft,
                                                                     const IJK_Field_double& pressure)
{
  /*
   * Should not change over iterations
   */
  init_ = init;
  sub_problem_index_ = sub_problem_index;

  if (init_)
    {
      associate_thermal_subproblem_parameters(ref_thermal_subresolution.debug_,
                                              ref_thermal_subresolution.n_iter_distance_,
                                              ref_thermal_subresolution.delta_T_subcooled_overheated_,
                                              ref_thermal_subresolution.pre_initialise_thermal_subproblems_list_);
      associate_eulerian_fields_references(interfaces,
                                           ref_thermal_subresolution.eulerian_distance_ns_,
                                           ref_thermal_subresolution.eulerian_curvature_ns_,
                                           ref_thermal_subresolution.eulerian_interfacial_area_ns_,
                                           ref_thermal_subresolution.eulerian_normal_vectors_ns_,
                                           ref_thermal_subresolution.eulerian_facets_barycentre_ns_,
                                           ref_thermal_subresolution.temperature_,
                                           ref_thermal_subresolution.temperature_ft_,
                                           ref_thermal_subresolution.temperature_before_extrapolation_,
                                           velocity,
                                           velocity_ft,
                                           pressure,
                                           ref_thermal_subresolution.grad_T_elem_,
                                           ref_thermal_subresolution.hess_diag_T_elem_,
                                           ref_thermal_subresolution.hess_cross_T_elem_);
      associate_probe_parameters(ref_thermal_subresolution.points_per_thermal_subproblem_,
                                 ref_thermal_subresolution.uniform_alpha_,
                                 ref_thermal_subresolution.uniform_lambda_,
                                 ref_thermal_subresolution.coeff_distance_diagonal_,
                                 ref_thermal_subresolution.cell_diagonal_,
                                 ref_thermal_subresolution.dr_,
                                 ref_thermal_subresolution.radial_coordinates_);
      associate_finite_difference_operators(ref_thermal_subresolution.radial_first_order_operator_raw_,
                                            ref_thermal_subresolution.radial_second_order_operator_raw_,
                                            ref_thermal_subresolution.radial_first_order_operator_,
                                            ref_thermal_subresolution.radial_second_order_operator_,
                                            ref_thermal_subresolution.identity_matrix_explicit_implicit_,
                                            ref_thermal_subresolution.identity_matrix_subproblems_,
                                            ref_thermal_subresolution.radial_diffusion_matrix_,
                                            ref_thermal_subresolution.radial_convection_matrix_);
      associate_finite_difference_solver_solution(ref_thermal_subresolution.finite_difference_assembler_,
                                                  ref_thermal_subresolution.thermal_subproblems_matrix_assembly_,
                                                  ref_thermal_subresolution.thermal_subproblems_rhs_assembly_,
                                                  ref_thermal_subresolution.thermal_subproblems_temperature_solution_,
                                                  ref_thermal_subresolution.thermal_subproblems_temperature_solution_ini_);
      associate_source_terms_parameters(ref_thermal_subresolution.source_terms_type_,
                                        ref_thermal_subresolution.source_terms_correction_,
                                        ref_thermal_subresolution.source_terms_correction_,
                                        ref_thermal_subresolution.advected_frame_of_reference_,
                                        ref_thermal_subresolution.neglect_frame_of_reference_radial_advection_);
      associate_flux_correction_parameters((ref_thermal_subresolution.convective_flux_correction_
                                            || ref_thermal_subresolution.diffusive_flux_correction_),
                                           ref_thermal_subresolution.distance_cell_faces_from_lrs_,
                                           ref_thermal_subresolution.interp_eulerian_);
      associate_varying_probes_params(ref_thermal_subresolution.first_time_step_varying_probes_,
                                      ref_thermal_subresolution.probe_variations_priority_,
                                      ref_thermal_subresolution.disable_interpolation_in_mixed_cells_);
      associate_flags_neighbours_correction(ref_thermal_subresolution.find_temperature_cell_neighbours_,
                                            !ref_thermal_subresolution.correct_neighbours_using_probe_length_,
                                            ref_thermal_subresolution.neighbours_corrected_rank_,
                                            ref_thermal_subresolution.neighbours_colinearity_weighting_,
                                            ref_thermal_subresolution.neighbours_distance_weighting_,
                                            ref_thermal_subresolution.neighbours_colinearity_distance_weighting_,
                                            ref_thermal_subresolution.neighbours_last_faces_colinearity_weighting_,
                                            ref_thermal_subresolution.neighbours_last_faces_colinearity_face_weighting_,
                                            ref_thermal_subresolution.neighbours_last_faces_distance_weighting_,
                                            ref_thermal_subresolution.neighbours_last_faces_distance_colinearity_weighting_,
                                            ref_thermal_subresolution.neighbours_last_faces_distance_colinearity_face_weighting_,
                                            ref_thermal_subresolution.find_reachable_fluxes_,
                                            ref_thermal_subresolution.find_cell_neighbours_for_fluxes_spherical_correction_);

//      // TODO: If the calculation of the distance is changed in intersection_ijk, it will be useless...
      associate_sub_problem_temporal_params(ref_thermal_subresolution.is_first_time_step_,
                                            ref_thermal_subresolution.first_time_step_temporal_,
                                            ref_thermal_subresolution.first_time_step_explicit_,
                                            ref_thermal_subresolution.local_fourier_,
                                            ref_thermal_subresolution.local_cfl_,
                                            ref_thermal_subresolution.min_delta_xyz_,
                                            ref_thermal_subresolution.max_u_radial_);
    }

  /*
   * Should be reinitialised at each time step.
   */
  clear_vectors();
  reset_counters();
  set_global_index(0);
  reset_post_processing_theta_phi_scope();
  associate_temporal_parameters(global_time_step, current_time);
  associate_cell_ijk(i, j, k);
  associate_eulerian_field_values(compo_connex, indicator);
  associate_interface_related_parameters(distance, curvature, interfacial_area, facet_barycentre, normal_vector);
  associate_rising_velocity(bubble_rising_velocity, bubble_rising_vector, bubble_barycentre);
  initialise_thermal_probe();
  if (!global_probes_characteristics_)
    (*first_time_step_temporal_) = 0;
  recompute_finite_difference_matrices();
}

void IJK_One_Dimensional_Subproblem::clear_vectors()
{
  pure_neighbours_to_correct_.clear();
  pure_neighbours_corrected_distance_.clear();
  pure_neighbours_corrected_colinearity_.clear();
  pure_neighbours_last_faces_to_correct_.clear();
  pure_neighbours_last_faces_corrected_distance_.clear();
  pure_neighbours_last_faces_corrected_colinearity_.clear();
}

void IJK_One_Dimensional_Subproblem::reset_counters()
{
  velocities_calculation_counter_ = 0;
}

void IJK_One_Dimensional_Subproblem::associate_thermal_subproblem_parameters(int debug,
                                                                             const int& n_iter_distance,
                                                                             const double& delta_T_subcooled_overheated,
                                                                             const int& pre_initialise_thermal_subproblems_list)
{
  debug_ = debug;
  n_iter_distance_ = n_iter_distance;
  delta_T_subcooled_overheated_ = delta_T_subcooled_overheated;
  pre_initialise_thermal_subproblems_list_ = pre_initialise_thermal_subproblems_list;
}

void IJK_One_Dimensional_Subproblem::associate_eulerian_fields_references(const IJK_Interfaces& interfaces,
                                                                          const IJK_Field_double& eulerian_distance,
                                                                          const IJK_Field_double& eulerian_curvature,
                                                                          const IJK_Field_double& eulerian_interfacial_area,
                                                                          const FixedVector<IJK_Field_double, 3>& eulerian_normal_vect,
                                                                          const FixedVector<IJK_Field_double, 3>& eulerian_facets_barycentre,
                                                                          const IJK_Field_double& temperature,
                                                                          const IJK_Field_double& temperature_ft,
                                                                          const IJK_Field_double& temperature_before_extrapolation,
                                                                          const FixedVector<IJK_Field_double, 3>& velocity,
                                                                          const FixedVector<IJK_Field_double, 3>& velocity_ft,
                                                                          const IJK_Field_double& pressure,
                                                                          const FixedVector<IJK_Field_double, 3>& grad_T_elem,
                                                                          const FixedVector<IJK_Field_double, 3>& hess_diag_T_elem,
                                                                          const FixedVector<IJK_Field_double, 3>& hess_cross_T_elem)
{
  interfaces_ = &interfaces;
  eulerian_distance_ = &eulerian_distance;
  eulerian_curvature_ = &eulerian_curvature;
  eulerian_interfacial_area_ = &eulerian_interfacial_area;
  eulerian_normal_vect_ = &eulerian_normal_vect;
  eulerian_facets_barycentre_ = &eulerian_facets_barycentre;
  temperature_ = &temperature ;
  temperature_ft_ = &temperature_ft ;
  temperature_before_extrapolation_ = &temperature_before_extrapolation;
  velocity_ = &velocity ;
  velocity_ft_ = &velocity_ft;
  pressure_ = &pressure;
  grad_T_elem_ = &grad_T_elem ;
  hess_diag_T_elem_ = &hess_diag_T_elem ;
  hess_cross_T_elem_ = &hess_cross_T_elem ;
}

void IJK_One_Dimensional_Subproblem::associate_sub_problem_temporal_params(const bool& is_first_time_step,
                                                                           int& first_time_step_temporal,
                                                                           const int& first_time_step_explicit,
                                                                           const double& local_fourier,
                                                                           const double& local_cfl,
                                                                           const double& min_delta_xyz,
                                                                           int max_u_radial)
{
  is_first_time_step_ = is_first_time_step;
  first_time_step_temporal_ = &first_time_step_temporal;
  first_time_step_explicit_ = first_time_step_explicit;
  local_fourier_ = local_fourier;
  local_cfl_ = local_cfl;
  min_delta_xyz_ = min_delta_xyz;
  max_u_radial_ = max_u_radial;
  max_u_cartesian_ = !max_u_radial_;
}

void IJK_One_Dimensional_Subproblem::associate_varying_probes_params(const int& first_time_step_varying_probes,
                                                                     const int& probe_variations_priority,
                                                                     const int& disable_interpolation_in_mixed_cells)
{
  first_time_step_varying_probes_ = first_time_step_varying_probes;
  probe_variations_priority_ = probe_variations_priority;
  disable_interpolation_in_mixed_cells_ = disable_interpolation_in_mixed_cells;
}

void IJK_One_Dimensional_Subproblem::associate_flux_correction_parameters(const int& correct_fluxes,
                                                                          const int& distance_cell_faces_from_lrs,
                                                                          const int& interp_eulerian)
{
  correct_fluxes_ = correct_fluxes;
  distance_cell_faces_from_lrs_ = distance_cell_faces_from_lrs;
  interp_eulerian_ = interp_eulerian;
}

void IJK_One_Dimensional_Subproblem::associate_source_terms_parameters(const int& source_terms_type,
                                                                       const int& correct_tangential_temperature_gradient,
                                                                       const int& correct_tangential_temperature_hessian,
                                                                       int advected_frame_of_reference,
                                                                       int neglect_frame_of_reference_radial_advection)
{
  source_terms_type_ = source_terms_type;
  correct_tangential_temperature_gradient_ = correct_tangential_temperature_gradient;
  correct_tangential_temperature_hessian_ = correct_tangential_temperature_hessian;
  advected_frame_of_reference_ = advected_frame_of_reference;
  neglect_frame_of_reference_radial_advection_ = neglect_frame_of_reference_radial_advection;
}

void 	IJK_One_Dimensional_Subproblem::associate_finite_difference_solver_solution(IJK_Finite_Difference_One_Dimensional_Matrix_Assembler& finite_difference_assembler,
                                                                                  Matrice& thermal_subproblems_matrix_assembly,
                                                                                  DoubleVect& thermal_subproblems_rhs_assembly,
                                                                                  DoubleVect& thermal_subproblems_temperature_solution,
                                                                                  DoubleVect& thermal_subproblems_temperature_solution_ini)
{
  finite_difference_assembler_ = &finite_difference_assembler;
  thermal_subproblems_matrix_assembly_ = &thermal_subproblems_matrix_assembly;
  thermal_subproblems_rhs_assembly_ = &thermal_subproblems_rhs_assembly;
  thermal_subproblems_temperature_solution_ = &thermal_subproblems_temperature_solution;
  thermal_subproblems_temperature_solution_ini_ = &thermal_subproblems_temperature_solution_ini;
}

void IJK_One_Dimensional_Subproblem::associate_temporal_parameters(const double& global_time_step, const double& current_time)
{
  global_time_step_ = global_time_step;
  current_time_ = current_time;
}

void IJK_One_Dimensional_Subproblem::associate_flags_neighbours_correction(const int& correct_temperature_cell_neighbours,
                                                                           const int& correct_neighbours_rank,
                                                                           const int& neighbours_corrected_rank,
                                                                           const int& neighbours_colinearity_weighting,
                                                                           const int& neighbours_distance_weighting,
                                                                           const int& neighbours_colinearity_distance_weighting,
                                                                           const int& neighbours_last_faces_colinearity_weighting,
                                                                           const int& neighbours_last_faces_colinearity_face_weighting,
                                                                           const int& neighbours_last_faces_distance_weighting,
                                                                           const int& neighbours_last_faces_distance_colinearity_weighting,
                                                                           const int& neighbours_last_faces_distance_colinearity_face_weighting,
                                                                           const int& compute_reachable_fluxes,
                                                                           const int& find_cell_neighbours_for_fluxes_spherical_correction)
{
  /*
   * Keep it under IJK_One_Dimensional_Subproblem::associate_flux_correction_parameters()
   */
  correct_temperature_cell_neighbours_ = correct_temperature_cell_neighbours;
  correct_neighbours_rank_ = correct_neighbours_rank;
  neighbours_corrected_rank_ = neighbours_corrected_rank;
  neighbours_colinearity_weighting_ = neighbours_colinearity_weighting;
  neighbours_distance_weighting_ = neighbours_distance_weighting;
  neighbours_colinearity_distance_weighting_ = neighbours_colinearity_distance_weighting;
  neighbours_weighting_ = (neighbours_colinearity_weighting_
                           || neighbours_distance_weighting_
                           || neighbours_colinearity_distance_weighting_);
  neighbours_last_faces_colinearity_weighting_ = neighbours_last_faces_colinearity_weighting;
  neighbours_last_faces_colinearity_face_weighting_ = neighbours_last_faces_colinearity_face_weighting;
  neighbours_last_faces_distance_weighting_ = neighbours_last_faces_distance_weighting;
  neighbours_last_faces_distance_colinearity_weighting_ = neighbours_last_faces_distance_colinearity_weighting;
  neighbours_last_faces_distance_colinearity_face_weighting_ = neighbours_last_faces_distance_colinearity_face_weighting;
  neighbours_last_faces_weighting_ = (neighbours_last_faces_colinearity_weighting_
                                      || neighbours_last_faces_colinearity_face_weighting_
                                      || neighbours_last_faces_distance_weighting_
                                      || neighbours_last_faces_distance_colinearity_weighting_
                                      || neighbours_last_faces_distance_colinearity_face_weighting_);
  compute_reachable_fluxes_ = compute_reachable_fluxes;
  find_cell_neighbours_for_fluxes_spherical_correction_ = find_cell_neighbours_for_fluxes_spherical_correction;
  correct_temperature_cell_neighbours_ = (correct_temperature_cell_neighbours_ && distance_cell_faces_from_lrs_);
}

void IJK_One_Dimensional_Subproblem::associate_probe_parameters(const int& points_per_thermal_subproblem,
                                                                const double& alpha,
                                                                const double& lambda,
                                                                const double& coeff_distance_diagonal,
                                                                const double& cell_diagonal,
                                                                const double& dr_base,
                                                                const DoubleVect& radial_coordinates)
{
  /*
   * If the probe characteristics are the same (don't copy attributes)
   */
  points_per_thermal_subproblem_base_ = &points_per_thermal_subproblem;
  coeff_distance_diagonal_ = &coeff_distance_diagonal;
  alpha_ = &alpha;
  lambda_ = &lambda;
  cell_diagonal_ = &cell_diagonal;
  dr_base_ = &dr_base;
  radial_coordinates_base_ = &radial_coordinates;

  /*
   * If each probe differs (create attributes !)
   */
  if (global_probes_characteristics_)
    points_per_thermal_subproblem_ = points_per_thermal_subproblem_base_;
  else
    points_per_thermal_subproblem_ = increase_number_of_points(); //copy if modified later
}

void IJK_One_Dimensional_Subproblem::associate_finite_difference_operators(const Matrice& radial_first_order_operator_raw,
                                                                           const Matrice& radial_second_order_operator_raw,
                                                                           const Matrice& radial_first_order_operator,
                                                                           const Matrice& radial_second_order_operator,
                                                                           const Matrice& identity_matrix_explicit_implicit_raw,
                                                                           Matrice& identity_matrix_subproblems,
                                                                           Matrice& radial_diffusion_matrix,
                                                                           Matrice& radial_convection_matrix)
{
  radial_first_order_operator_raw_base_ = &radial_first_order_operator_raw;
  radial_second_order_operator_raw_base_ = &radial_second_order_operator_raw;
  radial_first_order_operator_base_ = &radial_first_order_operator;
  radial_second_order_operator_base_ = &radial_second_order_operator;
  identity_matrix_explicit_implicit_base_ = &identity_matrix_explicit_implicit_raw;

  radial_first_order_operator_ = radial_first_order_operator_base_;
  radial_second_order_operator_ = radial_second_order_operator_base_;
  identity_matrix_explicit_implicit_ = identity_matrix_explicit_implicit_base_;

  identity_matrix_subproblems_ = &identity_matrix_subproblems;
  radial_diffusion_matrix_base_ = &radial_diffusion_matrix;
  radial_convection_matrix_base_ = &radial_convection_matrix;
}

void IJK_One_Dimensional_Subproblem::initialise_thermal_probe()
{
  if (debug_)
    Cerr << "Compute interface basis vectors" << finl;
  compute_interface_basis_vectors();

  if (debug_)
    Cerr << "Compute pure spherical basis vectors" << finl;
  compute_pure_spherical_basis_vectors();

  if (debug_)
    Cerr << "Compute probe parameters" << finl;
  /*
   *  Curvature is negative for a convex bubble
   *  but R should be positive in that case
   *  FIXME: What happen with highly deformable bubbles (concave interface portions) ?
   */
  if (fabs(curvature_) > DMINFLOAT)
    osculating_radius_ = fabs(2 / curvature_);

  probe_length_ = (*coeff_distance_diagonal_) * (*cell_diagonal_);

  if (debug_)
    Cerr << "Compute local discretisation" << finl;
  compute_local_discretisation();

  if (distance_cell_faces_from_lrs_)
    {
      if (debug_)
        Cerr << "Compute cell and faces distance to the interface" << finl;
      if (correct_fluxes_ || correct_temperature_cell_neighbours_ || find_cell_neighbours_for_fluxes_spherical_correction_ || compute_reachable_fluxes_)
        {
          compute_distance_cell_centre();
          compute_distance_faces_centres();
          Cerr << "Compute distance cell neighbours" << finl;
          if (correct_temperature_cell_neighbours_ || find_cell_neighbours_for_fluxes_spherical_correction_)
            compute_distance_cell_centres_neighbours();
          Cerr << "Compute distance faces neighbours" << finl;
          if (compute_reachable_fluxes_)
            compute_distance_last_cell_faces_neighbours();
        }
    }

  surface_ = (*eulerian_interfacial_area_)(index_i_, index_j_, index_k_);
  // Resize won't change the values if size is not changing...
  reinit_variable(rhs_assembly_);
//  rhs_assembly_.resize(*points_per_thermal_subproblem_);
//  rhs_assembly_ *= 0.;
}

void IJK_One_Dimensional_Subproblem::compute_interface_basis_vectors()
{
  /*
   * TODO: Associate a basis to each subproblem
   * Use Rodrigues' rotation formula to determine ephi ?
   * Needs an axis of (rotation gravity_dir x relative_vectors)
   * and an angle (gravity_dir dot relative_vectors) / (norm(gravity_dir)*norm(relative_vectors))
   * ephi is determined in the gravity_align rising direction
   * 		 | gravity_dir
   * 		 |
   *   *****
   * ***   ***
   * **     **
   * ***   ***
   *   *****
   *     |
   *     |
   */

  facet_barycentre_relative_ = facet_barycentre_ - bubble_barycentre_;
  if (debug_)
    {
      Cerr << "bubble_barycentre_"<< bubble_barycentre_[0] << " ; " << bubble_barycentre_[1] << " ; " << bubble_barycentre_[2] << finl;
      Cerr << "facet_barycentre_"<< facet_barycentre_[0] << " ; " << facet_barycentre_[1] << " ; " << facet_barycentre_[2] << finl;
      Cerr << "facet_barycentre_relative_"<< facet_barycentre_relative_[0] << " ; " << facet_barycentre_relative_[1] << " ; " << facet_barycentre_relative_[2] << finl;
    }
  Vecteur3 facet_barycentre_relative_normed = facet_barycentre_relative_;
  const double facet_barycentre_relative_norm = facet_barycentre_relative_normed.length();
  facet_barycentre_relative_normed *= (1 / facet_barycentre_relative_norm);
  Vecteur3 normal_contrib;
  const double normal_vector_compo_norm = normal_vector_compo_.length();
  normal_vector_compo_ *= (1 / normal_vector_compo_norm);
  Cerr << "Normal vector norm:" << normal_vector_compo_norm << finl;
  /*
   * First method with tangential direction of maximum velocity variations
   */
  DoubleTab facet_barycentre(1, 3);
  interfacial_velocity_compo_ = 0.;
  for (int dir=0; dir<3; dir++)
    facet_barycentre(0, dir) = facet_barycentre_[dir];
  for (int dir=0; dir<3; dir++)
    {
      DoubleVect interfacial_velocity_component(1);
      ijk_interpolate_skip_unknown_points((*velocity_)[dir], facet_barycentre, interfacial_velocity_component, INVALID_INTERP);
      interfacial_velocity_compo_[dir] = interfacial_velocity_component[0];
    }
  if (interfacial_velocity_compo_.length() < INVALID_VELOCITY)
    {
      normal_contrib = normal_vector_compo_;
      normal_contrib *= Vecteur3::produit_scalaire(facet_barycentre_relative_normed, normal_vector_compo_);
      first_tangential_vector_compo_ = facet_barycentre_relative_normed - normal_contrib;
    }
  else
    {
      // Should I remove the rising velocity ?
      interfacial_velocity_compo_ = interfacial_velocity_compo_ - bubble_rising_velocity_compo_;
      normal_contrib = normal_vector_compo_;
      normal_contrib *= Vecteur3::produit_scalaire(interfacial_velocity_compo_, normal_vector_compo_);
      interfacial_tangential_velocity_compo_ = interfacial_velocity_compo_ - normal_contrib;
      first_tangential_vector_compo_ = interfacial_tangential_velocity_compo_;
    }
  const double norm_first_tangential_vector = first_tangential_vector_compo_.length();
  first_tangential_vector_compo_ *= (1 / norm_first_tangential_vector);
  Vecteur3::produit_vectoriel(normal_vector_compo_, first_tangential_vector_compo_, second_tangential_vector_compo_);
  const double norm_second_tangential_vector = second_tangential_vector_compo_.length();
  second_tangential_vector_compo_ *= (1 / norm_second_tangential_vector);

  /*
   * Second method with rising velocity
   */
  Vecteur3::produit_vectoriel(bubble_rising_vector_, facet_barycentre_relative_, azymuthal_vector_compo_raw_);

  azymuthal_vector_compo_ = azymuthal_vector_compo_raw_;
  const double norm_azymuthal_vector_compo_raw_ = azymuthal_vector_compo_raw_.length();
//  const int sign_vector = signbit(Vecteur3::produit_scalaire(bubble_rising_vector_, normal_vector_compo_));
//  if (sign_vector)
//    azymuthal_vector_compo_ *= -1;
  azymuthal_vector_compo_ *= (1 / norm_azymuthal_vector_compo_raw_);

  normal_contrib = normal_vector_compo_;
  normal_contrib *=	Vecteur3::produit_scalaire(azymuthal_vector_compo_, normal_vector_compo_);
  azymuthal_vector_compo_ = azymuthal_vector_compo_ - normal_contrib;
  Vecteur3::produit_vectoriel(azymuthal_vector_compo_, normal_vector_compo_, first_tangential_vector_compo_from_rising_dir_);
  const double norm_first_tangential_vector_from_rising_dir = first_tangential_vector_compo_from_rising_dir_.length();
  first_tangential_vector_compo_from_rising_dir_ *= (1 / norm_first_tangential_vector_from_rising_dir);


  if (tangential_from_rising_vel_)
    {
      first_tangential_vector_compo_solver_ = &first_tangential_vector_compo_from_rising_dir_;
      second_tangential_vector_compo_solver_ = &azymuthal_vector_compo_;
      first_tangential_velocity_solver_ = &first_tangential_velocity_from_rising_dir_corrected_;
      second_tangential_velocity_solver_ = &azymuthal_velocity_corrected_;
      tangential_temperature_gradient_first_solver_ = &tangential_temperature_gradient_first_from_rising_dir_;
      tangential_temperature_gradient_second_solver_ = &azymuthal_temperature_gradient_;
    }
  else
    {
      // By default
      first_tangential_vector_compo_solver_= &first_tangential_vector_compo_;
      second_tangential_vector_compo_solver_ = &second_tangential_vector_compo_;
      first_tangential_velocity_solver_ = &first_tangential_velocity_corrected_;
      second_tangential_velocity_solver_ = &second_tangential_velocity_corrected_;
      tangential_temperature_gradient_first_solver_ = &tangential_temperature_gradient_first_;
      tangential_temperature_gradient_second_solver_ = &tangential_temperature_gradient_second_;
    }
}

void IJK_One_Dimensional_Subproblem::compute_pure_spherical_basis_vectors()
{
  /*
   * FIXME: It is align with gravity z but it should be modified to be align with the gravity dir ?
   */
  if (debug_)
    Cerr << "r_sph_ calculation"  << finl;
  r_sph_ = sqrt(facet_barycentre_relative_[0] * facet_barycentre_relative_[0]
                + facet_barycentre_relative_[1] * facet_barycentre_relative_[1]
                + facet_barycentre_relative_[2] * facet_barycentre_relative_[2]);
  if (debug_)
    {
      Cerr << "r_sph_ = " << r_sph_ << finl;
      Cerr << "theta_sph_ calculation"  << finl;
    }

//  theta_sph_ = atan(sqrt(facet_barycentre_relative_[0] * facet_barycentre_relative_[0]
//                         + facet_barycentre_relative_[1] * facet_barycentre_relative_[1])/ facet_barycentre_relative_[2]);
  theta_sph_ = atan2(sqrt(facet_barycentre_relative_[0] * facet_barycentre_relative_[0]
                          + facet_barycentre_relative_[1] * facet_barycentre_relative_[1]), facet_barycentre_relative_[2]);
  const double atan_theta_incr_ini = M_PI / 2;
  const double atan_incr_factor = -1;
  theta_sph_ = (theta_sph_ - atan_theta_incr_ini) * atan_incr_factor;

  if (debug_)
    {
      Cerr << "theta_sph_ = " << theta_sph_ << finl;
      Cerr << "phi_sph_ calculation"  << finl;
    }
  phi_sph_ = atan2(facet_barycentre_relative_[1], facet_barycentre_relative_[0]);

  if (debug_)
    {
      Cerr << "phi_sph_ = " << phi_sph_ << finl;
      Cerr << "er_sph_ calculation"  << finl;
    }
  for (int dir=0; dir<3; dir++)
    er_sph_[dir] = facet_barycentre_relative_[dir] / r_sph_;

  const double length = sqrt(facet_barycentre_relative_[0] * facet_barycentre_relative_[0]
                             + facet_barycentre_relative_[1] * facet_barycentre_relative_[1]);

  if (debug_)
    {
      Cerr << "er_sph_ = " << er_sph_[0] << finl;
      Cerr << "etheta_sph_ calculation"  << finl;
    }
  for (int dir=0; dir<2; dir++)
    etheta_sph_[dir] = facet_barycentre_relative_[dir] * facet_barycentre_relative_[2] / (r_sph_ * length);
  etheta_sph_[2] = - facet_barycentre_relative_[2] * length / r_sph_;

  ephi_sph_ = {0., 0., 0.};
  ephi_sph_[0] = - facet_barycentre_relative_[1];
  ephi_sph_[1] = facet_barycentre_relative_[0];
}

void IJK_One_Dimensional_Subproblem::compute_local_discretisation()
{
  int i;
  if (global_probes_characteristics_)
    {
      if (!probe_variations_enabled_)
        {
          if (!velocities_calculation_counter_)
            {
              radial_coordinates_ = radial_coordinates_base_;
              dr_ = *dr_base_;
            }
        }
      else
        {
          dr_ = probe_length_ / (*points_per_thermal_subproblem_ - 1);
          radial_coordinates_modified_.resize(*points_per_thermal_subproblem_);
          for (i=0; i < *points_per_thermal_subproblem_; i++)
            radial_coordinates_modified_(i) = i * dr_;
          radial_coordinates_ = &radial_coordinates_modified_;
        }
    }
  else
    {
      /*
       * coeff_distance_diagonal_ as well as
       * points_per_thermal_subproblem_ could be adapted
       */
      if (!probe_variations_enabled_)
        radial_coordinates_modified_.resize(*points_per_thermal_subproblem_);
      dr_ = probe_length_ / (*points_per_thermal_subproblem_ - 1);
      for (i=0; i < *points_per_thermal_subproblem_; i++)
        radial_coordinates_modified_(i) = i * dr_;
      radial_coordinates_ = &radial_coordinates_modified_;
    }
  /*
   * Following attributes differ anyway !
   */
  if (!velocities_calculation_counter_ || probe_variations_enabled_)
    {
      dr_inv_ = 1 / dr_;
      osculating_radial_coordinates_ = (*radial_coordinates_);
      osculating_radial_coordinates_ += osculating_radius_;
      if (!probe_variations_enabled_)
        {
          radial_coordinates_cartesian_compo_.resize(*points_per_thermal_subproblem_, 3);
          coordinates_cartesian_compo_.resize(*points_per_thermal_subproblem_, 3);
          osculating_radial_coordinates_cartesian_compo_.resize(*points_per_thermal_subproblem_, 3);
          osculating_radial_coordinates_inv_.resize(*points_per_thermal_subproblem_);
        }
      for (i=0; i < *points_per_thermal_subproblem_; i++)
        {
          osculating_radial_coordinates_inv_[i] = 1 / osculating_radial_coordinates_[i];
          for (int dir=0; dir<3; dir++)
            {
              radial_coordinates_cartesian_compo_(i, dir) = (*radial_coordinates_)(i) * normal_vector_compo_[dir];
              osculating_radial_coordinates_cartesian_compo_(i, dir) = osculating_radial_coordinates_(i) * normal_vector_compo_[dir];
              coordinates_cartesian_compo_(i, dir) = radial_coordinates_cartesian_compo_(i, dir) + facet_barycentre_[dir];
            }
        }
    }
}

void IJK_One_Dimensional_Subproblem::compute_local_time_step()
{
  if (*first_time_step_temporal_)
    {
      double max_u_inv = 1.e20;
      if (max_u_ > INVALID_VELOCITY_CFL)
        max_u_inv = 1 / max_u_;
      local_fourier_time_step_probe_length_ = probe_length_ * probe_length_ * local_fourier_ / (*alpha_) * 0.125; // factor 1/8 in 3D ?
      local_cfl_time_step_probe_length_ = probe_length_ * max_u_inv  * local_cfl_ * 0.5; // factor 1/2 in 3D ?
      local_dt_cfl_min_delta_xyz_ = min_delta_xyz_ * max_u_inv  * local_cfl_ * 0.5;
      if (first_time_step_explicit_)
        {
          local_dt_fo_ = dr_ * dr_ * local_fourier_ / (*alpha_);

          local_dt_cfl_ = dr_ * max_u_inv  * local_cfl_;
          local_time_step_ = std::min(local_dt_cfl_, local_dt_fo_);
          local_time_step_ = std::min(local_time_step_, global_time_step_);
          if (is_first_time_step_)
            nb_iter_explicit_ = (int) (global_time_step_ / local_time_step_);
          else
            nb_iter_explicit_ = (int) ((current_time_ + global_time_step_) / local_time_step_);
          nb_iter_explicit_++;
          if (is_first_time_step_)
            local_time_step_round_ = global_time_step_ / (double) nb_iter_explicit_;
          else
            local_time_step_round_ = (current_time_ + global_time_step_) / (double) nb_iter_explicit_;
          local_time_step_round_ /= (double) (*points_per_thermal_subproblem_);
          nb_iter_explicit_ *= (*points_per_thermal_subproblem_);
        }
      else
        {
          local_time_step_ = global_time_step_;
          local_time_step_round_ = local_time_step_;
        }
    }
}

const int * IJK_One_Dimensional_Subproblem::increase_number_of_points()
{
  increased_point_numbers_ = *points_per_thermal_subproblem_base_;
  return &increased_point_numbers_;
}

void IJK_One_Dimensional_Subproblem::compute_identity_matrix_local(Matrice& identity_matrix_explicit_implicit)
{
  int check_nb_elem;
  check_nb_elem = (*finite_difference_assembler_).build(identity_matrix_explicit_implicit, *points_per_thermal_subproblem_, -1);
  Cerr << "Check_nb_elem" << check_nb_elem << finl;
}

void IJK_One_Dimensional_Subproblem::compute_first_order_operator_local(Matrice& radial_first_order_operator)
{
  int check_nb_elem;
  check_nb_elem = (*finite_difference_assembler_).build(radial_first_order_operator, *points_per_thermal_subproblem_, 0);
  Cerr << "Check_nb_elem" << check_nb_elem << finl;
  radial_first_order_operator_local_ *= dr_inv_;
}

void IJK_One_Dimensional_Subproblem::compute_second_order_operator_local(Matrice& radial_second_order_operator)
{
  int check_nb_elem;
  check_nb_elem = (*finite_difference_assembler_).build(radial_second_order_operator, *points_per_thermal_subproblem_, 0);
  Cerr << "Check_nb_elem" << check_nb_elem << finl;
  const double dr_squared_inv = 1 / pow(dr_, 2);
  radial_second_order_operator_local_ *= dr_squared_inv;
}

void IJK_One_Dimensional_Subproblem::recompute_finite_difference_matrices()
{
  if (!global_probes_characteristics_)
    {
      compute_identity_matrix_local(identity_matrix_explicit_implicit_local_);
      compute_first_order_operator_local(radial_first_order_operator_local_);
      compute_second_order_operator_local(radial_second_order_operator_local_);
      identity_matrix_explicit_implicit_local_ = (*identity_matrix_explicit_implicit_base_);
      identity_matrix_explicit_implicit_ = &identity_matrix_explicit_implicit_local_;
      radial_first_order_operator_ = &radial_first_order_operator_local_;
      radial_second_order_operator_ = &radial_second_order_operator_local_;
    }
}

void IJK_One_Dimensional_Subproblem::compute_first_order_operator_local_varying_probe_length(const Matrice * radial_first_order_operator)
{
  radial_first_order_operator_local_ = (*radial_first_order_operator);
  radial_first_order_operator_local_ *= dr_inv_;
}

void IJK_One_Dimensional_Subproblem::compute_second_order_operator_local_varying_probe_length(const Matrice * radial_second_order_operator)
{
  radial_second_order_operator_local_ = (*radial_second_order_operator);
  const double dr_squared_inv = 1 / pow(dr_, 2);
  radial_second_order_operator_local_ *= dr_squared_inv;
}

void IJK_One_Dimensional_Subproblem::recompute_finite_difference_matrices_varying_probe_length()
{
  compute_first_order_operator_local_varying_probe_length(radial_first_order_operator_raw_base_);
  compute_second_order_operator_local_varying_probe_length(radial_second_order_operator_raw_base_);
  radial_second_order_operator_= &radial_second_order_operator_local_;
}

void IJK_One_Dimensional_Subproblem::interpolate_project_velocities_on_probes()
{
  if (!velocities_calculation_counter_ || probe_variations_enabled_)
    {
      radial_velocity_.resize(*points_per_thermal_subproblem_);
      first_tangential_velocity_.resize(*points_per_thermal_subproblem_);
      second_tangential_velocity_.resize(*points_per_thermal_subproblem_);
      azymuthal_velocity_.resize(*points_per_thermal_subproblem_);
      first_tangential_velocity_from_rising_dir_.resize(*points_per_thermal_subproblem_);

      radial_velocity_corrected_.resize(*points_per_thermal_subproblem_);
      first_tangential_velocity_corrected_.resize(*points_per_thermal_subproblem_);
      second_tangential_velocity_corrected_.resize(*points_per_thermal_subproblem_);
      azymuthal_velocity_corrected_.resize(*points_per_thermal_subproblem_);
      first_tangential_velocity_from_rising_dir_corrected_.resize(*points_per_thermal_subproblem_);

      pressure_interp_.resize(*points_per_thermal_subproblem_);

      x_velocity_.resize(*points_per_thermal_subproblem_);
      y_velocity_.resize(*points_per_thermal_subproblem_);
      z_velocity_.resize(*points_per_thermal_subproblem_);

      x_velocity_corrected_.resize(*points_per_thermal_subproblem_);
      y_velocity_corrected_.resize(*points_per_thermal_subproblem_);
      z_velocity_corrected_.resize(*points_per_thermal_subproblem_);

      interpolate_pressure_on_probes();
      interpolate_cartesian_velocities_on_probes();
      compute_velocity_magnitude();
      project_velocities_on_probes();
      velocities_calculation_counter_++;
    }
}

void IJK_One_Dimensional_Subproblem::reajust_probe_length()
{
  if (first_time_step_varying_probes_)
    compute_modified_probe_length_condition();
}

void IJK_One_Dimensional_Subproblem::compute_modified_probe_length_condition()
{
  const double current_time = ref_ijk_ft_->get_current_time();
  cfl_probe_length_ = (max_u_ * current_time) / local_cfl_; // Add 3D constants ?
  fourier_probe_length_ = sqrt(((*alpha_ * (current_time + global_time_step_))) / local_fourier_); // Add 3D constants ?
  max_cfl_fourier_probe_length_ = std::max(cfl_probe_length_, fourier_probe_length_);
  /*
   * TODO: ADD constraint on the temperature because of the displacement of the interface !!!!
   */
  cell_temperature_ = (*temperature_before_extrapolation_)(index_i_, index_j_, index_k_);
  if (max_cfl_fourier_probe_length_ < probe_length_)
    {
      if (debug_)
        Cerr << "Probe length should be modified" << finl;
      if (!correct_fluxes_)
        {
          if (indicator_ < 0.5)
            {
              compute_distance_cell_centre();
              assert(cell_centre_distance_ >= 0.);
              if (max_cfl_fourier_probe_length_ < cell_centre_distance_)
                short_probe_condition_ = 1;
              else
                short_probe_condition_ = 0;
            }
        }
      else
        {
          compute_distance_faces_centres();
          bool has_liquid_neighbours = 1;
          for (int i=0; i<6; i++)
            has_liquid_neighbours = has_liquid_neighbours && pure_liquid_neighbours_[i];
          const double max_distance_pure_face_centre = compute_max_distance_pure_face_centre();
          const double max_distance_pure_vertex_centre = compute_max_distance_pure_face_vertices();
          const double max_distance_face_centre_vertex = std::max(max_distance_pure_face_centre, max_distance_pure_vertex_centre);
          if (max_cfl_fourier_probe_length_ < max_distance_face_centre_vertex || !has_liquid_neighbours)
            {
              short_probe_condition_ = 1;
              if (cell_temperature_ != delta_T_subcooled_overheated_)
                temperature_probe_condition_ = 1;
              else
                temperature_probe_condition_ = 0;
            }
          else
            short_probe_condition_= 0;
        }
      probe_variations_enabled_ = 1;
    }
  else
    probe_variations_enabled_ = 0;
}

/*
 * TODO: Use ijk_intersections_interface instead !
 * and avoid redundancy...
 * Compare Calculation with mean_over_compo()
 * not weighted by the surface
 * Separate geometric attributes and physical attributes -> in the future
 */

void IJK_One_Dimensional_Subproblem::compute_distance_cell_centre()
{
  Vecteur3 centre = ref_ijk_ft_->get_splitting_ns().get_coords_of_dof(index_i_, index_j_, index_k_, IJK_Splitting::ELEM);

  Vecteur3 facet_to_cell_centre = facet_barycentre_;
  facet_to_cell_centre *= -1;
  facet_to_cell_centre += centre;
  cell_centre_distance_ = Vecteur3::produit_scalaire(facet_to_cell_centre, normal_vector_compo_);

  Vecteur3 normal_contrib = normal_vector_compo_;
  normal_contrib *= cell_centre_distance_;
  Vecteur3 tangential_displacement = normal_contrib;
  tangential_displacement *= (-1);
  tangential_displacement += facet_to_cell_centre;
  cell_centre_tangential_distance_ = tangential_displacement.length();
  tangential_distance_vector_ = tangential_displacement;
  if (cell_centre_tangential_distance_ > 1e-16)
    tangential_distance_vector_ *= (1 / cell_centre_tangential_distance_);
}

void IJK_One_Dimensional_Subproblem::compute_distance_faces_centres()
{
  Vecteur3 bary_face {0., 0., .0};
  Vecteur3 vector_relative {0., 0., 0.};
  Vecteur3 bary_vertex {0., 0., 0.};
  int neighbours_i[6] = NEIGHBOURS_I;
  int neighbours_j[6] = NEIGHBOURS_J;
  int neighbours_k[6] = NEIGHBOURS_K;
  int neighbours_faces_i[6] = NEIGHBOURS_FACES_I;
  int neighbours_faces_j[6] = NEIGHBOURS_FACES_J;
  int neighbours_faces_k[6] = NEIGHBOURS_FACES_K;
  int face_dir[6] = FACES_DIR;
  int m;
  for (int l=0; l<6; l++)
    {
      const int ii = neighbours_i[l];
      const int jj = neighbours_j[l];
      const int kk = neighbours_k[l];
      const double indic_neighbour = ref_ijk_ft_->itfce().I()(index_i_+ii, index_j_+jj, index_k_+kk);
      if (fabs(indic_neighbour) > LIQUID_INDICATOR_TEST)
        {
          const int ii_f = neighbours_faces_i[l];
          const int jj_f = neighbours_faces_j[l];
          const int kk_f = neighbours_faces_k[l];
          pure_liquid_neighbours_[l] = 1;
          if (ii)
            bary_face = ref_ijk_ft_->get_splitting_ns().get_coords_of_dof(index_i_+ii_f, index_j_+jj_f, index_k_+kk_f, IJK_Splitting::FACES_I);
          if (jj)
            bary_face = ref_ijk_ft_->get_splitting_ns().get_coords_of_dof(index_i_+ii_f, index_j_+jj_f, index_k_+kk_f, IJK_Splitting::FACES_J);
          if (kk)
            bary_face = ref_ijk_ft_->get_splitting_ns().get_coords_of_dof(index_i_+ii_f, index_j_+jj_f, index_k_+kk_f, IJK_Splitting::FACES_K);
          vector_relative = facet_barycentre_;
          vector_relative *= (-1);
          vector_relative += bary_face;
          const double distance_face_centre = Vecteur3::produit_scalaire(vector_relative, normal_vector_compo_);
          face_centres_distance_[l] = distance_face_centre;
          for (m=0; m<4; m++)
            {
              double distance_vertex_centre = 0.;
              bary_vertex = vector_relative;
              compute_vertex_position(m, face_dir[l], bary_face, distance_vertex_centre, bary_vertex);
              vertices_centres_distance_[l][m] = distance_vertex_centre;
            }
        }
      else
        {
          pure_liquid_neighbours_[l] = 0;
          face_centres_distance_[l] = 0.;
          for (m=0; m<4; m++)
            vertices_centres_distance_[l][m] = 0.;
        }
    }
}

double IJK_One_Dimensional_Subproblem::compute_min_distance_pure_face_centre()
{
  double min_face_centre_distance = 1.e20;
  for (int l=0; l<6; l++)
    if (pure_liquid_neighbours_[l])
      if(face_centres_distance_[l] > 0)
        min_face_centre_distance = std::min(min_face_centre_distance, face_centres_distance_[l]);
  return min_face_centre_distance;
}

double IJK_One_Dimensional_Subproblem::compute_min_distance_pure_face_vertices()
{
  double min_face_vertex_distance = 1.e20;
  int m;
  for (int l=0; l<6; l++)
    if (pure_liquid_neighbours_[l])
      for (m=0; m<4; m++)
        if(vertices_centres_distance_[l][m] > 0)
          min_face_vertex_distance = std::min(min_face_vertex_distance, vertices_centres_distance_[l][m]);
  return min_face_vertex_distance;
}

double IJK_One_Dimensional_Subproblem::compute_max_distance_pure_face_centre()
{
  double max_face_centre_distance = 0.;
  for (int l=0; l<6; l++)
    if (pure_liquid_neighbours_[l])
      if(face_centres_distance_[l] > 0)
        max_face_centre_distance = std::max(max_face_centre_distance, face_centres_distance_[l]);
  return max_face_centre_distance;
}

double IJK_One_Dimensional_Subproblem::compute_max_distance_pure_face_vertices()
{
  double max_face_vertex_distance = 0.;
  int m;
  for (int l=0; l<6; l++)
    if (pure_liquid_neighbours_[l])
      for (m=0; m<4; m++)
        if(vertices_centres_distance_[l][m] > 0)
          max_face_vertex_distance = std::max(max_face_vertex_distance, vertices_centres_distance_[l][m]);
  return max_face_vertex_distance;
}

void IJK_One_Dimensional_Subproblem::compute_vertex_position(const int& vertex_number,
                                                             const int& face_dir,
                                                             const Vecteur3& bary_face,
                                                             double& distance_vertex_centre,
                                                             Vecteur3& bary_vertex)
{
  const IJK_Grid_Geometry& geom = ref_ijk_ft_->get_splitting_ns().get_grid_geometry();
  const double dx = geom.get_constant_delta(DIRECTION_I);
  const double dy = geom.get_constant_delta(DIRECTION_J);
  const double dz = geom.get_constant_delta(DIRECTION_K);
  const double neighbours_first_dir[4] = NEIGHBOURS_FIRST_DIR;
  const double neighbours_second_dir[4] = NEIGHBOURS_SECOND_DIR;
  Vecteur3 point_coords {0., 0., 0.};
  double dl1;
  double dl2;
  switch(face_dir)
    {
    case 0:
      dl1 = dy / 2.;
      dl2 = dz / 2.;
      point_coords[1] = dl1 * neighbours_first_dir[vertex_number];
      point_coords[2] = dl2 * neighbours_second_dir[vertex_number];
      break;
    case 1:
      dl1 = dx / 2.;
      dl2 = dz / 2.;
      point_coords[0] = dl1 * neighbours_first_dir[vertex_number];
      point_coords[2] = dl2 * neighbours_second_dir[vertex_number];
      break;
    case 2:
      dl1 = dx / 2.;
      dl2 = dy / 2.;
      point_coords[0] = dl1 * neighbours_first_dir[vertex_number];
      point_coords[1] = dl2 * neighbours_second_dir[vertex_number];
      break;
    default:
      dl1 = dx / 2.;
      dl2 = dy / 2.;
      point_coords[0] = dl1 * neighbours_first_dir[vertex_number];
      point_coords[1] = dl2 * neighbours_second_dir[vertex_number];
      break;
    }
  bary_vertex += point_coords;
  distance_vertex_centre = Vecteur3::produit_scalaire(bary_vertex, normal_vector_compo_);
}

void IJK_One_Dimensional_Subproblem::compute_distance_cell_centres_neighbours()
{
  const IJK_Grid_Geometry& geom = ref_ijk_ft_->get_splitting_ns().get_grid_geometry();
  const double dx = geom.get_constant_delta(DIRECTION_I);
  const double dy = geom.get_constant_delta(DIRECTION_J);
  const double dz = geom.get_constant_delta(DIRECTION_K);
  int l, m, n;

  int dxyz_increment_max = get_dxyz_increment_max();
  /*
   * 8-1 values for one neighbour in each dir... (Too much, enhance later)
   * Positive OR Negative dir depending on the normal vector
   */
  pure_neighbours_to_correct_.resize(dxyz_increment_max + 1);
  pure_neighbours_corrected_distance_.resize(dxyz_increment_max + 1);
  if (neighbours_colinearity_weighting_)
    pure_neighbours_corrected_colinearity_.resize(dxyz_increment_max + 1);
  for (l=dxyz_increment_max; l>=0; l--)
    {
      pure_neighbours_to_correct_[l].resize(dxyz_increment_max + 1);
      pure_neighbours_corrected_distance_[l].resize(dxyz_increment_max + 1);
      if (neighbours_colinearity_weighting_)
        pure_neighbours_corrected_colinearity_[l].resize(dxyz_increment_max + 1);
      for (m=dxyz_increment_max; m>=0; m--)
        {
          pure_neighbours_to_correct_[l][m].resize(dxyz_increment_max + 1);
          pure_neighbours_corrected_distance_[l][m].resize(dxyz_increment_max + 1);
          if (neighbours_colinearity_weighting_)
            pure_neighbours_corrected_colinearity_[l][m].resize(dxyz_increment_max + 1);
          for (n=dxyz_increment_max; n>=0; n--)
            {
              pure_neighbours_to_correct_[l][m][n] = false;
              pure_neighbours_corrected_distance_[l][m][n] = 0.;
              if (neighbours_colinearity_weighting_)
                pure_neighbours_corrected_colinearity_[l][m][n] = 0.;
            }
        }
    }

  double remaining_distance_diag = probe_length_ - cell_centre_distance_;
  Vecteur3 remaining_distance_diag_projected = normal_vector_compo_;
  remaining_distance_diag_projected *= remaining_distance_diag;
  for (int i=0; i<3; i++)
    pure_neighbours_corrected_sign_[i] = signbit(normal_vector_compo_[i]);
  int dx_increment = (int) abs(remaining_distance_diag_projected[0] / dx);
  int dy_increment = (int) abs(remaining_distance_diag_projected[1] / dy);
  int dz_increment = (int) abs(remaining_distance_diag_projected[2] / dz);
  if (correct_neighbours_rank_)
    {
      dx_increment = std::min(dx_increment, dxyz_increment_max);
      dy_increment = std::min(dy_increment, dxyz_increment_max);
      dz_increment = std::min(dz_increment, dxyz_increment_max);
    }
  dxyz_increment_bool_ = (dx_increment || dy_increment || dz_increment);
  for (l=dx_increment; l>=0; l--)
    for (m=dy_increment; m>=0; m--)
      for (n=dz_increment; n>=0; n--)
        if (l!=0 || m!=0 || n!=0)
          {
            const int l_dir = (pure_neighbours_corrected_sign_[0]) ? l * (-1) : l;
            const int m_dir = (pure_neighbours_corrected_sign_[1]) ? m * (-1) : m;
            const int n_dir = (pure_neighbours_corrected_sign_[2]) ? n * (-1) : n;
            const double indic_neighbour = ref_ijk_ft_->itfce().I()(index_i_ + l_dir, index_j_ + m_dir, index_k_ + n_dir);
            if (indic_neighbour > LIQUID_INDICATOR_TEST)
              {
                pure_neighbours_to_correct_[l][m][n] = true;
                const double dx_contrib = l_dir * dx;
                const double dy_contrib = m_dir * dy;
                const double dz_contrib = n_dir * dz;
                Vecteur3 distance_contrib = {dx_contrib, dy_contrib, dz_contrib};
                pure_neighbours_corrected_distance_[l][m][n] = cell_centre_distance_
                                                               + Vecteur3::produit_scalaire(normal_vector_compo_, distance_contrib);
                if (neighbours_colinearity_weighting_)
                  {
                    const double colinearity = compute_cell_faces_weighting(dx_contrib, dy_contrib, dz_contrib, 0);
                    pure_neighbours_corrected_colinearity_[l][m][n] = colinearity;
                  }
              }
          }
}

double IJK_One_Dimensional_Subproblem::compute_cell_weighting(const double& dx_contrib,
                                                              const double& dy_contrib,
                                                              const double& dz_contrib)
{
  if (neighbours_colinearity_weighting_)
    return compute_colinearity(dx_contrib, dy_contrib, dz_contrib);
  if (neighbours_distance_weighting_)
    return compute_distance_cell_faces(dx_contrib, dy_contrib, dz_contrib);
  if (neighbours_colinearity_distance_weighting_)
    return compute_colinearity(dx_contrib, dy_contrib, dz_contrib) * compute_distance_cell_faces(dx_contrib, dy_contrib, dz_contrib);
  return 1;
}

void IJK_One_Dimensional_Subproblem::compute_distance_last_cell_faces_neighbours()
{
  const IJK_Grid_Geometry& geom = ref_ijk_ft_->get_splitting_ns().get_grid_geometry();
  const double dx = geom.get_constant_delta(DIRECTION_I);
  const double dy = geom.get_constant_delta(DIRECTION_J);
  const double dz = geom.get_constant_delta(DIRECTION_K);
  const double dx_over_two = dx / 2.;
  const double dy_over_two = dy / 2.;
  const double dz_over_two = dz / 2.;
  int l, m, n;
  int l_cell, m_cell, n_cell;

  int dxyz_increment_max = get_dxyz_increment_max();
  int dxyz_over_two_increment_max = get_dxyz_over_two_increment_max();
  const int first_increment[3] = {dxyz_over_two_increment_max + 1, dxyz_increment_max, dxyz_increment_max};
  const int second_increment[3] = {dxyz_increment_max, dxyz_over_two_increment_max + 1, dxyz_increment_max};
  const int third_increment[3] = {dxyz_increment_max, dxyz_increment_max, dxyz_over_two_increment_max + 1};
  //  dxyz_over_two_increment_max *= 2;
  //  if (!dxyz_over_two_increment_max%2)
  //    dxyz_over_two_increment_max -= 1;

  /*
   * 8-1 values for one neighbour in each dir... (Too much, enhance later)
   * Positive OR Negative dir depending on the normal vector
   */
  pure_neighbours_last_faces_to_correct_.resize(3);
  pure_neighbours_last_faces_corrected_distance_.resize(3);
  pure_neighbours_last_faces_corrected_colinearity_.resize(3);
  for (int c=0; c<3; c++)
    {
      const int first_incr = first_increment[c];
      const int second_incr = second_increment[c];
      const int third_incr = third_increment[c];
      pure_neighbours_last_faces_to_correct_[c].resize(first_incr + 1);
      pure_neighbours_last_faces_corrected_distance_[c].resize(first_incr + 1);
      pure_neighbours_last_faces_corrected_colinearity_[c].resize(first_incr + 1);
      for (l=first_incr; l>=0; l--)
        {
          pure_neighbours_last_faces_to_correct_[c][l].resize(second_incr + 1);
          pure_neighbours_last_faces_corrected_distance_[c][l].resize(second_incr + 1);
          pure_neighbours_last_faces_corrected_colinearity_[c][l].resize(second_incr + 1);
          for (m=second_incr; m>=0; m--)
            {
              pure_neighbours_last_faces_to_correct_[c][l][m].resize(third_incr + 1);
              pure_neighbours_last_faces_corrected_distance_[c][l][m].resize(third_incr + 1);
              pure_neighbours_last_faces_corrected_colinearity_[c][l][m].resize(third_incr + 1);
              for (n=third_incr; n>=0; n--)
                {
                  pure_neighbours_last_faces_to_correct_[c][l][m][n] = false;
                  pure_neighbours_last_faces_corrected_distance_[c][l][m][n] = 0.;
                  pure_neighbours_last_faces_corrected_colinearity_[c][l][m][n] = 0.;
                }
            }
        }
    }

  double remaining_distance_diag = probe_length_ - cell_centre_distance_;
  Vecteur3 remaining_distance_diag_projected = normal_vector_compo_;
  remaining_distance_diag_projected *= remaining_distance_diag;
  for (int i=0; i<3; i++)
    pure_neighbours_corrected_sign_[i] = signbit(normal_vector_compo_[i]);
  int dx_over_two_increment = (int) abs(remaining_distance_diag_projected[0] / dx_over_two);
  int dy_over_two_increment = (int) abs(remaining_distance_diag_projected[1] / dy_over_two);
  int dz_over_two_increment = (int) abs(remaining_distance_diag_projected[2] / dz_over_two);
  int dx_increment = (int) (dx_over_two_increment / 2);
  int dy_increment = (int) (dy_over_two_increment / 2);
  int dz_increment = (int) (dz_over_two_increment / 2);
  dx_over_two_increment -= dx_increment;
  dy_over_two_increment -= dy_increment;
  dz_over_two_increment -= dz_increment;
  if (correct_neighbours_rank_)
    {
      dx_over_two_increment = std::min(dx_over_two_increment, dxyz_over_two_increment_max);
      dy_over_two_increment = std::min(dy_over_two_increment, dxyz_over_two_increment_max);
      dz_over_two_increment = std::min(dz_over_two_increment, dxyz_over_two_increment_max);
      dx_increment = std::min(dx_increment, dxyz_increment_max);
      dy_increment = std::min(dy_increment, dxyz_increment_max);
      dz_increment = std::min(dz_increment, dxyz_increment_max);
    }
  dxyz_over_two_increment_bool_ = (dx_over_two_increment >0 || dy_over_two_increment>0 || dz_over_two_increment >0 );
  if (dx_over_two_increment>0)
    dx_over_two_increment--;
  if (dy_over_two_increment>0)
    dy_over_two_increment--;
  if (dz_over_two_increment>0)
    dz_over_two_increment--;

  /*
   * TODO: Should we look to cells faces that are not in the normal vector direction ?
   */
  for (l=dx_over_two_increment + 1; l>=0; l--)
    for (m_cell=dy_increment; m_cell>=0; m_cell--)
      for (n_cell=dz_increment; n_cell>=0; n_cell--)
        {
          const int l_dir = (pure_neighbours_corrected_sign_[0]) ? l * (-1) + 1 : l;
          const int m_dir = (pure_neighbours_corrected_sign_[1]) ? m_cell * (-1) : m_cell;
          const int n_dir = (pure_neighbours_corrected_sign_[2]) ? n_cell * (-1) : n_cell;
          const int l_dir_elem = (pure_neighbours_corrected_sign_[0]) ? (l + 1) * (-1) + 1 : l;
          const double indic_neighbour = ref_ijk_ft_->itfce().I()(index_i_ + l_dir_elem, index_j_ + m_dir, index_k_ + n_dir);
          if (indic_neighbour > LIQUID_INDICATOR_TEST)
            {
              pure_neighbours_last_faces_to_correct_[0][l][m_cell][n_cell] = true;
              const double lmn_zero = (l > 0) ? 1. : 0.;
              const double contrib_factor = (pure_neighbours_corrected_sign_[0]) ? (lmn_zero * (2 * abs(l_dir) + 1) - (1. - lmn_zero)) * (-1):
                                            lmn_zero * (2 * (l_dir - 1) + 1) - (1. - lmn_zero);
              const double dx_contrib = contrib_factor * dx_over_two;
              const double dy_contrib = m_dir * dy;
              const double dz_contrib = n_dir * dz;
              Vecteur3 distance_contrib = {dx_contrib, dy_contrib, dz_contrib};
              pure_neighbours_last_faces_corrected_distance_[0][l][m_cell][n_cell] = cell_centre_distance_
                                                                                     + Vecteur3::produit_scalaire(normal_vector_compo_, distance_contrib);
              if (pure_neighbours_last_faces_corrected_distance_[0][l][m_cell][n_cell] < 0)
                pure_neighbours_last_faces_to_correct_[0][l][m_cell][n_cell] = false;
              if (neighbours_last_faces_weighting_)
                {
                  const double colinearity = compute_cell_faces_weighting(dx_contrib, dy_contrib, dz_contrib, 0);
                  pure_neighbours_last_faces_corrected_colinearity_[0][l][m_cell][n_cell] = colinearity;
                }
            }
        }

  for (l_cell=dx_increment; l_cell>=0; l_cell--)
    for (m=dy_over_two_increment + 1; m>=0; m--)
      for (n_cell=dz_increment; n_cell>=0; n_cell--)
        {
          const int l_dir = (pure_neighbours_corrected_sign_[0]) ? l_cell * (-1) : l_cell;
          const int m_dir = (pure_neighbours_corrected_sign_[1]) ? m * (-1) + 1: m;
          const int n_dir = (pure_neighbours_corrected_sign_[2]) ? n_cell * (-1) : n_cell;
          const int m_dir_elem = (pure_neighbours_corrected_sign_[1]) ? (m + 1) * (-1) + 1 : m;
          const double indic_neighbour = ref_ijk_ft_->itfce().I()(index_i_ + l_dir, index_j_ + m_dir_elem, index_k_ + n_dir);
          if (indic_neighbour > LIQUID_INDICATOR_TEST)
            {
              pure_neighbours_last_faces_to_correct_[1][l_cell][m][n_cell] = true;
              const double lmn_zero = (m > 0) ? 1. : 0.;
              const double contrib_factor = (pure_neighbours_corrected_sign_[1]) ? (lmn_zero * (2 * abs(m_dir) + 1) - (1. - lmn_zero)) * (-1):
                                            lmn_zero * (2 * (m_dir - 1) + 1) - (1. - lmn_zero);
              const double dx_contrib = l_dir * dx;
              const double dy_contrib = contrib_factor * dy_over_two;
              const double dz_contrib = n_dir * dz;
              Vecteur3 distance_contrib = {dx_contrib, dy_contrib, dz_contrib};
              pure_neighbours_last_faces_corrected_distance_[1][l_cell][m][n_cell] = cell_centre_distance_
                                                                                     + Vecteur3::produit_scalaire(normal_vector_compo_, distance_contrib);
              if (pure_neighbours_last_faces_corrected_distance_[1][l_cell][m][n_cell] < 0)
                pure_neighbours_last_faces_to_correct_[1][l_cell][m][n_cell] = false;
              if (neighbours_last_faces_weighting_)
                {
                  const double colinearity = compute_cell_faces_weighting(dx_contrib, dy_contrib, dz_contrib, 1);
                  pure_neighbours_last_faces_corrected_colinearity_[1][l_cell][m][n_cell] = colinearity;
                }
            }
        }

  for (l_cell=dx_increment; l_cell>=0; l_cell--)
    for (m_cell=dy_increment; m_cell>=0; m_cell--)
      for (n=dz_over_two_increment + 1; n>=0; n--)
        {
          const int l_dir = (pure_neighbours_corrected_sign_[0]) ? l_cell * (-1) : l_cell;
          const int m_dir = (pure_neighbours_corrected_sign_[1]) ? m_cell * (-1) : m_cell;
          const int n_dir = (pure_neighbours_corrected_sign_[2]) ? n * (-1) + 1: n;
          const int n_dir_elem = (pure_neighbours_corrected_sign_[2]) ? (n + 1) * (-1) + 1 : n;
          const double indic_neighbour = ref_ijk_ft_->itfce().I()(index_i_ + l_dir, index_j_ + m_dir, index_k_ + n_dir_elem);
          if (indic_neighbour > LIQUID_INDICATOR_TEST)
            {
              pure_neighbours_last_faces_to_correct_[2][l_cell][m_cell][n] = true;
              const double lmn_zero = (n > 0) ? 1. : 0.;
              const double contrib_factor = (pure_neighbours_corrected_sign_[2]) ? (lmn_zero * (2 * abs(n_dir) + 1) - (1. - lmn_zero)) * (-1):
                                            lmn_zero * (2 * (n_dir - 1) + 1) - (1. - lmn_zero);
              const double dx_contrib = l_dir * dx;
              const double dy_contrib = m_dir * dy;
              const double dz_contrib = contrib_factor * dz_over_two;
              Vecteur3 distance_contrib = {dx_contrib, dy_contrib, dz_contrib};
              pure_neighbours_last_faces_corrected_distance_[2][l_cell][m_cell][n] = cell_centre_distance_ +
                                                                                     Vecteur3::produit_scalaire(normal_vector_compo_, distance_contrib);
              if (pure_neighbours_last_faces_corrected_distance_[2][l_cell][m_cell][n] < 0)
                pure_neighbours_last_faces_to_correct_[2][l_cell][m_cell][n] = false;
              if (neighbours_last_faces_weighting_)
                {
                  const double colinearity = compute_cell_faces_weighting(dx_contrib, dy_contrib, dz_contrib, 2);
                  pure_neighbours_last_faces_corrected_colinearity_[2][l_cell][m_cell][n] = colinearity;
                }
            }
        }
//  if (neighbours_last_faces_colinearity_weighting_)
//    {
//      Vecteur3 relative_vector = normal_vector_compo_;
//      relative_vector *= cell_centre_distance_;
//      relative_vector[0] += ((l + 1) * normal_vector_compo_[0] * dx_over_two);
//      relative_vector[1] += ((m + 1) * normal_vector_compo_[1] * dy_over_two);
//      relative_vector[2] += ((n + 1) * normal_vector_compo_[2] * dz_over_two);
//      const double relative_vector_norm = relative_vector.length();
//      relative_vector *= (1 / relative_vector_norm);
//      // pure_neighbours_corrected_colinearity_[l][m][n] = Vecteur3::produit_scalaire(normal_vector_compo_, relative_vector);
//    }
}

double IJK_One_Dimensional_Subproblem::compute_cell_faces_weighting(const double& dx_contrib,
                                                                    const double& dy_contrib,
                                                                    const double& dz_contrib,
                                                                    const int& dir)
{
  if (neighbours_last_faces_colinearity_weighting_)
    return compute_colinearity(dx_contrib, dy_contrib, dz_contrib);
  if (neighbours_last_faces_colinearity_face_weighting_)
    return compute_colinearity_cell_faces(dx_contrib, dy_contrib, dz_contrib, dir);
  if (neighbours_last_faces_distance_weighting_)
    return compute_distance_cell_faces(dx_contrib, dy_contrib, dz_contrib);
  if (neighbours_last_faces_distance_colinearity_weighting_)
    return compute_colinearity(dx_contrib, dy_contrib, dz_contrib) * compute_distance_cell_faces(dx_contrib, dy_contrib, dz_contrib);
  if (neighbours_last_faces_distance_colinearity_face_weighting_)
    return compute_colinearity_cell_faces(dx_contrib, dy_contrib, dz_contrib, dir) * compute_distance_cell_faces(dx_contrib, dy_contrib, dz_contrib);
  return 1;
}

Vecteur3 IJK_One_Dimensional_Subproblem::compute_relative_vector_cell_faces(const double& dx_contrib,
                                                                            const double& dy_contrib,
                                                                            const double& dz_contrib)
{
  Vecteur3 relative_vector = normal_vector_compo_;
  relative_vector *= cell_centre_distance_;
  Vecteur3 tangential_relative_vector = tangential_distance_vector_;
  tangential_relative_vector *= cell_centre_tangential_distance_;
  relative_vector += tangential_relative_vector;
  Vecteur3 dxyz_contrib = {dx_contrib, dy_contrib, dz_contrib};
  relative_vector += dxyz_contrib;
  return relative_vector;
}

double IJK_One_Dimensional_Subproblem::compute_colinearity(const double& dx_contrib, const double& dy_contrib, const double& dz_contrib)
{
  Vecteur3 relative_vector = compute_relative_vector_cell_faces(dx_contrib, dy_contrib, dz_contrib);
  const double relative_vector_norm = relative_vector.length();
  relative_vector *= (1 / relative_vector_norm);
  const double colinearity = Vecteur3::produit_scalaire(normal_vector_compo_, relative_vector);
  return abs(colinearity);
}

double IJK_One_Dimensional_Subproblem::compute_colinearity_cell_faces(const double& dx_contrib,
                                                                      const double& dy_contrib,
                                                                      const double& dz_contrib,
                                                                      const int& dir)
{
  Vecteur3 relative_vector = compute_relative_vector_cell_faces(dx_contrib, dy_contrib, dz_contrib);
  const double relative_vector_norm = relative_vector.length();
  relative_vector *= (1 / relative_vector_norm);
  return abs(relative_vector[dir]);
}

double IJK_One_Dimensional_Subproblem::compute_distance_cell_faces(const double& dx_contrib,
                                                                   const double& dy_contrib,
                                                                   const double& dz_contrib)
{
  Vecteur3 relative_vector = compute_relative_vector_cell_faces(dx_contrib, dy_contrib, dz_contrib);
  const double distance = Vecteur3::produit_scalaire(tangential_distance_vector_, relative_vector);
  return abs(1 / (distance + 1e-16));
}

int IJK_One_Dimensional_Subproblem::get_dxyz_increment_max()
{
  const IJK_Grid_Geometry& geom = ref_ijk_ft_->get_splitting_ns().get_grid_geometry();
  const double dx = geom.get_constant_delta(DIRECTION_I);
  const double dy = geom.get_constant_delta(DIRECTION_J);
  const double dz = geom.get_constant_delta(DIRECTION_K);

  int dxyz_increment_max;
  if (!correct_neighbours_rank_)
    {
      const int dx_increment_max = (int) ((dx + probe_length_) / dx);
      const int dy_increment_max = (int) ((dy + probe_length_) / dy);
      const int dz_increment_max = (int) ((dz + probe_length_) / dz);
      dxyz_increment_max = std::max(std::max(dx_increment_max, dy_increment_max), dz_increment_max);
    }
  else
    {
      dxyz_increment_max = neighbours_corrected_rank_;
    }
  return dxyz_increment_max;
}

int IJK_One_Dimensional_Subproblem::get_dxyz_over_two_increment_max()
{
  const IJK_Grid_Geometry& geom = ref_ijk_ft_->get_splitting_ns().get_grid_geometry();
  const double dx = geom.get_constant_delta(DIRECTION_I);
  const double dy = geom.get_constant_delta(DIRECTION_J);
  const double dz = geom.get_constant_delta(DIRECTION_K);

  int dxyz_over_two_increment_max;
  if (!correct_neighbours_rank_)
    {
      const int dx_increment_max = (int) ((probe_length_ + dx) / (dx / 2.));
      const int dy_increment_max = (int) ((probe_length_ + dy) / (dy / 2.));
      const int dz_increment_max = (int) ((probe_length_ + dz) / (dz / 2.));
      dxyz_over_two_increment_max = std::max(std::max(dx_increment_max, dy_increment_max), dz_increment_max);
    }
  else
    {
      dxyz_over_two_increment_max = neighbours_face_corrected_rank_;
    }
  return dxyz_over_two_increment_max;
}

void IJK_One_Dimensional_Subproblem::compute_modified_probe_length(const int& probe_variations_enabled)
{
  if (probe_variations_enabled && probe_variations_enabled_)
    {
      probe_length_ = max_cfl_fourier_probe_length_;
      compute_local_discretisation();
      recompute_finite_difference_matrices_varying_probe_length();
      // probe_variations_enabled_ = 0;
    }
  else
    {
      probe_variations_enabled_ = 0;
      first_time_step_varying_probes_ = 0;
    }
}

void IJK_One_Dimensional_Subproblem::interpolate_pressure_on_probes()
{
  ijk_interpolate_skip_unknown_points((*pressure_), coordinates_cartesian_compo_, pressure_interp_, INVALID_INTERP);
}

void IJK_One_Dimensional_Subproblem::interpolate_cartesian_velocities_on_probes()
{
  DoubleVect * vel_compo = nullptr;
  for (int dir = 0; dir < 3; dir++)
    {
      switch (dir)
        {
        case 0:
          vel_compo = &x_velocity_;
          break;
        case 1:
          vel_compo = &y_velocity_;
          break;
        case 2:
          vel_compo = &z_velocity_;
          break;
        default:
          vel_compo = &x_velocity_;
          break;
        }
      ijk_interpolate_skip_unknown_points((*velocity_)[dir], coordinates_cartesian_compo_, *vel_compo, INVALID_INTERP);
    }
}

void IJK_One_Dimensional_Subproblem::compute_velocity_magnitude()
{
  velocity_magnitude_ = x_velocity_;
  velocity_magnitude_ *= x_velocity_;
  DoubleVect velocity_magnitude_y = y_velocity_;
  velocity_magnitude_y *= y_velocity_;
  velocity_magnitude_ += velocity_magnitude_y;
  DoubleVect velocity_magnitude_z = z_velocity_;
  velocity_magnitude_z *= z_velocity_;
  velocity_magnitude_ += velocity_magnitude_z;
  for(int i=0; i<(*points_per_thermal_subproblem_); i++)
    {
      const double velocity_magnitude_sqrt = sqrt(velocity_magnitude_[i]);
      velocity_magnitude_[i] = velocity_magnitude_sqrt;
    }
}

void IJK_One_Dimensional_Subproblem::project_velocities_on_probes()
{
  project_cartesian_onto_basis_vector(x_velocity_, y_velocity_, z_velocity_, normal_vector_compo_, radial_velocity_);
  project_cartesian_onto_basis_vector(x_velocity_, y_velocity_, z_velocity_, first_tangential_vector_compo_, first_tangential_velocity_);
  project_cartesian_onto_basis_vector(x_velocity_, y_velocity_, z_velocity_, second_tangential_vector_compo_, second_tangential_velocity_);
  project_cartesian_onto_basis_vector(x_velocity_, y_velocity_, z_velocity_, azymuthal_vector_compo_, azymuthal_velocity_);
  project_cartesian_onto_basis_vector(x_velocity_, y_velocity_, z_velocity_, first_tangential_vector_compo_from_rising_dir_,
                                      first_tangential_velocity_from_rising_dir_);
//  project_cartesian_onto_basis_vector(x_velocity_, y_velocity_, z_velocity_, *first_tangential_vector_compo_solver_, first_tangential_velocity_);
//	project_cartesian_onto_basis_vector(x_velocity_, y_velocity_, z_velocity_, *second_tangential_vector_compo_solver_, second_tangential_velocity_);

  correct_velocities();

  max_u_ = INVALID_VELOCITY_CFL * 0.9;
  for (int i=0; i<radial_velocity_corrected_.size(); i++)
    if (max_u_cartesian_)
      max_u_ = std::max(max_u_, fabs(velocity_magnitude_[i]));
    else
      max_u_ = std::max(max_u_, fabs(radial_velocity_corrected_[i]));
  compute_local_time_step();

  reinit_variable(x_velocity_corrected_);
  reinit_variable(y_velocity_corrected_);
  reinit_variable(z_velocity_corrected_);
  project_basis_vector_onto_cartesian_dir(0, first_tangential_velocity_, second_tangential_velocity_, radial_velocity_corrected_,
                                          *first_tangential_vector_compo_solver_, *second_tangential_vector_compo_solver_, normal_vector_compo_,
                                          x_velocity_corrected_);
  project_basis_vector_onto_cartesian_dir(1, first_tangential_velocity_, second_tangential_velocity_, radial_velocity_corrected_,
                                          *first_tangential_vector_compo_solver_, *second_tangential_vector_compo_solver_, normal_vector_compo_,
                                          y_velocity_corrected_);
  project_basis_vector_onto_cartesian_dir(2, first_tangential_velocity_, second_tangential_velocity_, radial_velocity_corrected_,
                                          *first_tangential_vector_compo_solver_, *second_tangential_vector_compo_solver_, normal_vector_compo_,
                                          z_velocity_corrected_);
}

void IJK_One_Dimensional_Subproblem::correct_velocities()
{
  correct_velocity(radial_velocity_, radial_velocity_advected_frame_);
  correct_velocity(first_tangential_velocity_, first_tangential_velocity_advected_frame_);
  correct_velocity(second_tangential_velocity_, second_tangential_velocity_advected_frame_);
  correct_velocity(first_tangential_velocity_from_rising_dir_, first_tangential_velocity_from_rising_dir_advected_frame_);
  correct_velocity(azymuthal_velocity_, azymuthal_velocity_advected_frame_);

  correct_velocity_rise(radial_velocity_, normal_vector_compo_, radial_velocity_static_frame_);
  correct_velocity_rise(first_tangential_velocity_, first_tangential_vector_compo_, first_tangential_velocity_static_frame_);
  correct_velocity_rise(second_tangential_velocity_, second_tangential_vector_compo_, second_tangential_velocity_static_frame_);
  correct_velocity_rise(first_tangential_velocity_from_rising_dir_, first_tangential_vector_compo_from_rising_dir_, first_tangential_velocity_from_rising_dir_static_frame_);
  correct_velocity_rise(azymuthal_velocity_, azymuthal_vector_compo_, azymuthal_velocity_static_frame_);

  if (advected_frame_of_reference_)
    {
      radial_velocity_corrected_ = radial_velocity_advected_frame_;
      first_tangential_velocity_corrected_ = first_tangential_velocity_advected_frame_;
      second_tangential_velocity_corrected_ = second_tangential_velocity_advected_frame_;
      first_tangential_velocity_from_rising_dir_corrected_ = first_tangential_velocity_from_rising_dir_advected_frame_;
      azymuthal_velocity_corrected_ = azymuthal_velocity_advected_frame_;
    }
  else
    {
      if (neglect_frame_of_reference_radial_advection_)
        radial_velocity_corrected_ = radial_velocity_static_frame_;
      else
        radial_velocity_corrected_ = radial_velocity_advected_frame_;
      first_tangential_velocity_corrected_ = first_tangential_velocity_static_frame_;
      second_tangential_velocity_corrected_ = second_tangential_velocity_static_frame_;
      first_tangential_velocity_from_rising_dir_corrected_ = first_tangential_velocity_from_rising_dir_static_frame_;
      azymuthal_velocity_corrected_ = azymuthal_velocity_static_frame_;
    }
}

void IJK_One_Dimensional_Subproblem::correct_velocity(const DoubleVect& velocity, DoubleVect& velocity_corrected)
{
  velocity_corrected = velocity;
  for (int i=0; i<velocity_corrected.size(); i++)
    velocity_corrected[i] -= velocity[0];
}

void IJK_One_Dimensional_Subproblem::correct_velocity_rise(const DoubleVect& velocity, const Vecteur3& basis, DoubleVect& velocity_corrected)
{
  DoubleVect bubble_rising_velocity_projection(1);
  DoubleVect bubble_rising_velocity_compo_x(1);
  DoubleVect bubble_rising_velocity_compo_y(1);
  DoubleVect bubble_rising_velocity_compo_z(1);
  bubble_rising_velocity_compo_x[0] = bubble_rising_velocity_compo_[0];
  bubble_rising_velocity_compo_y[0] = bubble_rising_velocity_compo_[1];
  bubble_rising_velocity_compo_z[0] = bubble_rising_velocity_compo_[2];
  project_cartesian_onto_basis_vector(bubble_rising_velocity_compo_x,bubble_rising_velocity_compo_y, bubble_rising_velocity_compo_z,
                                      basis, bubble_rising_velocity_projection);
  velocity_corrected = velocity;
  for (int i=0; i<velocity_corrected.size(); i++)
    velocity_corrected[i] -= bubble_rising_velocity_projection[0];
}

void IJK_One_Dimensional_Subproblem::correct_radial_velocity_probe()
{
  correct_velocity(radial_velocity_, radial_velocity_corrected_);
}

void IJK_One_Dimensional_Subproblem::project_cartesian_onto_basis_vector(const DoubleVect& compo_x, const DoubleVect& compo_y, const DoubleVect& compo_z, const Vecteur3& basis, DoubleVect& projection)
{
  const int size_x = compo_x.size();
//  const int size_y = compo_y.size();
//  const int size_z = compo_z.size();
//  const int size_projection = projection.size();
//  assert((size_x == size_y) && (size_x == size_z) && (size_x == size_projection));
//  Cerr << "size_x << size_y << size_z << size_projection" <<  finl;
//  Cerr << size_x << size_y << size_z << size_projection <<  finl;
  for (int i=0; i<size_x; i++)
    projection[i] = compo_x[i] * basis[0] + compo_y[i] * basis[1] + compo_z[i] * basis[2];
}

void IJK_One_Dimensional_Subproblem::project_basis_vector_onto_cartesian_dir(const int& dir, const DoubleVect& compo_u, const DoubleVect& compo_v, const DoubleVect& compo_w,
                                                                             const Vecteur3& basis_u, const Vecteur3& basis_v, const Vecteur3& basis_w,
                                                                             DoubleVect& projection)
{
  const int size_u = compo_u.size();
  const int size_v = compo_v.size();
  const int size_w = compo_w.size();
  for (int i=0; i<projection.size(); i++)
    {
      if (i< size_u)
        projection[i] += (compo_u[i] * basis_u[dir]);
      if (i< size_v)
        projection[i] += (compo_v[i] * basis_v[dir]);
      if (i< size_w)
        projection[i] += (compo_w[i] * basis_w[dir]);
    }
}

void IJK_One_Dimensional_Subproblem::interpolate_temperature_on_probe()
{
  temperature_interp_.resize(*points_per_thermal_subproblem_);
  if (first_time_step_varying_probes_)
    ijk_interpolate_skip_unknown_points(*temperature_before_extrapolation_, coordinates_cartesian_compo_, temperature_interp_, INVALID_INTERP);
  else
    ijk_interpolate_skip_unknown_points(*temperature_, coordinates_cartesian_compo_, temperature_interp_, INVALID_INTERP);
}

void IJK_One_Dimensional_Subproblem::interpolate_temperature_gradient_on_probe()
{
  for (int dir = 0; dir < 3; dir++)
    {
      grad_T_elem_interp_[dir].resize(*points_per_thermal_subproblem_);
      ijk_interpolate_skip_unknown_points((*grad_T_elem_)[dir], coordinates_cartesian_compo_, grad_T_elem_interp_[dir], INVALID_INTERP);
    }
}

void IJK_One_Dimensional_Subproblem::project_temperature_gradient_on_probes()
{
  normal_temperature_gradient_.resize(*points_per_thermal_subproblem_);
  project_cartesian_onto_basis_vector(grad_T_elem_interp_[0], grad_T_elem_interp_[1], grad_T_elem_interp_[2], normal_vector_compo_, normal_temperature_gradient_);

  tangential_temperature_gradient_first_.resize(*points_per_thermal_subproblem_);
  tangential_temperature_gradient_second_.resize(*points_per_thermal_subproblem_);
  project_cartesian_onto_basis_vector(grad_T_elem_interp_[0], grad_T_elem_interp_[1], grad_T_elem_interp_[2],
                                      first_tangential_vector_compo_, tangential_temperature_gradient_first_);
  project_cartesian_onto_basis_vector(grad_T_elem_interp_[0], grad_T_elem_interp_[1], grad_T_elem_interp_[2],
                                      second_tangential_vector_compo_, tangential_temperature_gradient_second_);

  azymuthal_temperature_gradient_.resize(*points_per_thermal_subproblem_);
  tangential_temperature_gradient_first_from_rising_dir_.resize(*points_per_thermal_subproblem_);
  project_cartesian_onto_basis_vector(grad_T_elem_interp_[0], grad_T_elem_interp_[1], grad_T_elem_interp_[2],
                                      azymuthal_vector_compo_, azymuthal_temperature_gradient_);
  project_cartesian_onto_basis_vector(grad_T_elem_interp_[0], grad_T_elem_interp_[1], grad_T_elem_interp_[2],
                                      first_tangential_vector_compo_from_rising_dir_, tangential_temperature_gradient_first_from_rising_dir_);
}

void IJK_One_Dimensional_Subproblem::interpolate_temperature_hessian_on_probe()
{
  for (int dir = 0; dir < 3; dir++)
    {
      hess_diag_T_elem_interp_[dir].resize(*points_per_thermal_subproblem_);
      hess_cross_T_elem_interp_[dir].resize(*points_per_thermal_subproblem_);
      hess_diag_T_elem_spherical_[dir].resize(*points_per_thermal_subproblem_);
      hess_cross_T_elem_spherical_[dir].resize(*points_per_thermal_subproblem_);
      hess_diag_T_elem_spherical_from_rising_[dir].resize(*points_per_thermal_subproblem_);
      hess_cross_T_elem_spherical_from_rising_[dir].resize(*points_per_thermal_subproblem_);
      ijk_interpolate_skip_unknown_points((*hess_diag_T_elem_)[dir], coordinates_cartesian_compo_, hess_diag_T_elem_interp_[dir], INVALID_INTERP);
      ijk_interpolate_skip_unknown_points((*hess_cross_T_elem_)[dir], coordinates_cartesian_compo_, hess_cross_T_elem_interp_[dir], INVALID_INTERP);
    }
  temperature_diffusion_hessian_cartesian_trace_ = hess_diag_T_elem_interp_[0];
  temperature_diffusion_hessian_cartesian_trace_ += hess_diag_T_elem_interp_[1];
  temperature_diffusion_hessian_cartesian_trace_ += hess_diag_T_elem_interp_[2];
}

void IJK_One_Dimensional_Subproblem::project_temperature_hessian_on_probes()
{
  compute_projection_matrix_cartesian_to_local_spherical();
  /*
   * A' = P^tAP
   */
  Matrice33 temperature_hessian_cartesian;
  Matrice33 temperature_hessian_spherical;
  Matrice33 temperature_hessian_spherical_from_rising;
  for (int i=0; i<*points_per_thermal_subproblem_; i++)
    {
      temperature_hessian_cartesian = Matrice33(hess_diag_T_elem_interp_[0](i), hess_cross_T_elem_interp_[2](i), hess_cross_T_elem_interp_[1](i),
                                                hess_cross_T_elem_interp_[2](i), hess_diag_T_elem_interp_[1](i), hess_cross_T_elem_interp_[0](i),
                                                hess_cross_T_elem_interp_[1](i), hess_cross_T_elem_interp_[0](i), hess_diag_T_elem_interp_[2](i));

      project_matrix_on_basis(projection_matrix_, inverse_projection_matrix_, temperature_hessian_cartesian, temperature_hessian_spherical);
      hess_diag_T_elem_spherical_[0](i) = temperature_hessian_spherical(0, 0);
      hess_diag_T_elem_spherical_[1](i) = temperature_hessian_spherical(1, 1);
      hess_diag_T_elem_spherical_[2](i) = temperature_hessian_spherical(2, 2);
      hess_cross_T_elem_spherical_[0](i) = temperature_hessian_spherical(1, 2);
      hess_cross_T_elem_spherical_[1](i) = temperature_hessian_spherical(0, 2);
      hess_cross_T_elem_spherical_[2](i) = temperature_hessian_spherical(0, 1);

      project_matrix_on_basis(projection_matrix_from_rising_, inverse_projection_matrix_from_rising_, temperature_hessian_cartesian, temperature_hessian_spherical_from_rising);
      hess_diag_T_elem_spherical_from_rising_[0](i) = temperature_hessian_spherical_from_rising(0, 0);
      hess_diag_T_elem_spherical_from_rising_[1](i) = temperature_hessian_spherical_from_rising(1, 1);
      hess_diag_T_elem_spherical_from_rising_[2](i) = temperature_hessian_spherical_from_rising(2, 2);
      hess_cross_T_elem_spherical_from_rising_[0](i) = temperature_hessian_spherical_from_rising(1, 2);
      hess_cross_T_elem_spherical_from_rising_[1](i) = temperature_hessian_spherical_from_rising(0, 2);
      hess_cross_T_elem_spherical_from_rising_[2](i) = temperature_hessian_spherical_from_rising(0, 1);
    }
}

void IJK_One_Dimensional_Subproblem::retrieve_temperature_diffusion_spherical_on_probes()
{
  temperature_diffusion_hessian_trace_ = hess_diag_T_elem_spherical_[0];
  temperature_diffusion_hessian_trace_ += hess_diag_T_elem_spherical_[1];
  temperature_diffusion_hessian_trace_ += hess_diag_T_elem_spherical_[2];
  radial_temperature_diffusion_ = normal_temperature_gradient_;
  radial_temperature_diffusion_ *= osculating_radial_coordinates_inv_;
  radial_temperature_diffusion_ *= 2;
  radial_temperature_diffusion_ += hess_diag_T_elem_spherical_[0];
  tangential_temperature_diffusion_ = radial_temperature_diffusion_;
  tangential_temperature_diffusion_ *= -1;
  tangential_temperature_diffusion_ += temperature_diffusion_hessian_trace_;
}

void IJK_One_Dimensional_Subproblem::compute_projection_matrix_cartesian_to_local_spherical()
{
  /*
   * X = PX' <-> X'=P^tX
   * Easier to calculate the transpose/inverse P^t
   * Change of basis matrix send the coordinate to the basis of the LHS
   * X=P(beta->beta')X'
   * A' = P^tAP
   */
  inverse_projection_matrix_ = Matrice33(normal_vector_compo_[0], normal_vector_compo_[1], normal_vector_compo_[2],
                                         first_tangential_vector_compo_[0], first_tangential_vector_compo_[1], first_tangential_vector_compo_[2],
                                         second_tangential_vector_compo_[0], second_tangential_vector_compo_[1], second_tangential_vector_compo_[2]);
  Matrice33::inverse(inverse_projection_matrix_, projection_matrix_, 0);

  inverse_projection_matrix_from_rising_ = Matrice33(normal_vector_compo_[0], normal_vector_compo_[1], normal_vector_compo_[2],
                                                     first_tangential_vector_compo_from_rising_dir_[0], first_tangential_vector_compo_from_rising_dir_[1], first_tangential_vector_compo_from_rising_dir_[2],
                                                     azymuthal_vector_compo_[0], azymuthal_vector_compo_[1], azymuthal_vector_compo_[2]);
  Matrice33::inverse(inverse_projection_matrix_from_rising_, projection_matrix_from_rising_, 0);
}

void IJK_One_Dimensional_Subproblem::project_matrix_on_basis(const Matrice33& projection_matrix, const Matrice33& inverse_projection_matrix, const Matrice33& matrix, Matrice33& projected_matrix)
{
  Matrice33 matrix_ap_tmp = matrix;
  // AP
  Matrice33::produit_matriciel(matrix, projection_matrix, matrix_ap_tmp);
  // P^tAP
  Matrice33::produit_matriciel(inverse_projection_matrix, matrix_ap_tmp, projected_matrix);
}

void IJK_One_Dimensional_Subproblem::initialise_radial_convection_operator_local()
{
  if (!global_probes_characteristics_ || first_time_step_varying_probes_ || pre_initialise_thermal_subproblems_list_)
    (*finite_difference_assembler_).reinitialise_matrix_subproblem(radial_convection_matrix_base_,
                                                                   radial_first_order_operator_,
                                                                   sub_problem_index_);
}

void IJK_One_Dimensional_Subproblem::initialise_radial_diffusion_operator_local()
{
  if (!global_probes_characteristics_ || first_time_step_varying_probes_ || pre_initialise_thermal_subproblems_list_)
    (*finite_difference_assembler_).reinitialise_matrix_subproblem(radial_diffusion_matrix_base_,
                                                                   radial_second_order_operator_,
                                                                   sub_problem_index_);
}

void IJK_One_Dimensional_Subproblem::initialise_identity_operator_local()
{
  if ((!global_probes_characteristics_ || pre_initialise_thermal_subproblems_list_) && (*first_time_step_temporal_))
    (*finite_difference_assembler_).reinitialise_matrix_subproblem(identity_matrix_subproblems_,
                                                                   identity_matrix_explicit_implicit_,
                                                                   sub_problem_index_);
}

void IJK_One_Dimensional_Subproblem::compute_radial_convection_diffusion_operators()
{
  const double alpha_inv = - 1 / *alpha_;
  DoubleVect osculating_radial_coefficient = osculating_radial_coordinates_inv_;
  osculating_radial_coefficient *= 2;
  reinit_variable(radial_convection_prefactor_);
  // radial_convection_prefactor_.resize(*points_per_thermal_subproblem_);
  // interpolate_project_velocities_on_probes();
  if (source_terms_type_ != linear_diffusion
      && source_terms_type_ != spherical_diffusion
      && source_terms_type_ != spherical_diffusion_approx)
    {
      if (correct_radial_velocity_)
        radial_convection_prefactor_ = radial_velocity_corrected_;
      else
        radial_convection_prefactor_ = radial_velocity_;
      radial_convection_prefactor_ *= alpha_inv;
    }
  else
    Cerr << "Diffusion case : don't compute the radial convection pre-factor" << finl;
  if (source_terms_type_ != linear_diffusion)
    {
      if (source_terms_type_ == spherical_diffusion_approx)
        radial_convection_prefactor_ +=	2 / osculating_radius_;
      else
        radial_convection_prefactor_ +=	osculating_radial_coefficient;
    }
  const int boundary_conditions = 0;
  initialise_radial_convection_operator_local();
  initialise_radial_diffusion_operator_local();
  initialise_identity_operator_local();
  (*finite_difference_assembler_).scale_matrix_subproblem_by_vector(radial_convection_matrix_base_,
                                                                    radial_convection_prefactor_,
                                                                    sub_problem_index_,
                                                                    boundary_conditions);
}

void IJK_One_Dimensional_Subproblem::prepare_temporal_schemes()
{
  if (*first_time_step_temporal_)
    {
      if (first_time_step_explicit_)
        (*finite_difference_assembler_).apply_euler_time_step(radial_convection_matrix_base_,
                                                              radial_diffusion_matrix_base_,
                                                              sub_problem_index_,
                                                              local_time_step_overall_,
                                                              (*alpha_));
      else
        (*finite_difference_assembler_).correct_sign_temporal_schemes_subproblems(radial_convection_matrix_base_,
                                                                                  radial_diffusion_matrix_base_,
                                                                                  sub_problem_index_,
                                                                                  local_time_step_overall_,
                                                                                  (*alpha_));
    }
}

void IJK_One_Dimensional_Subproblem::prepare_boundary_conditions(DoubleVect& thermal_subproblems_rhs_assembly,
                                                                 DoubleVect& thermal_subproblems_temperature_solution_ini,
                                                                 const int& boundary_condition_interface,
                                                                 const double& interfacial_boundary_condition_value,
                                                                 const int& impose_boundary_condition_interface_from_simulation,
                                                                 const int& boundary_condition_end,
                                                                 const double& end_boundary_condition_value,
                                                                 const int& impose_user_boundary_condition_end_value)
{
  interpolate_temperature_on_probe();

  for (int i=0; i<temperature_interp_.size(); i++)
    if (fabs(temperature_interp_[i]) > INVALID_INTERP_TEST)
      Cerr << "Error in the temperature_interpolation" << temperature_interp_[i] << finl;

  /*
   * Written for Dirichlet only
   * TODO: Write for other B.Cs or vapour-liquid coupled system
   */
  if (!impose_boundary_condition_interface_from_simulation)
    interfacial_boundary_condition_value_ = interfacial_boundary_condition_value;
  else
    interfacial_boundary_condition_value_ = temperature_interp_[0];
  if (impose_user_boundary_condition_end_value)
    end_boundary_condition_value_ = end_boundary_condition_value;
  else
    end_boundary_condition_value_ = temperature_interp_[temperature_interp_.size() -1];

  const int thermal_subproblems_rhs_size = (*thermal_subproblems_rhs_assembly_).size();
  if (pre_initialise_thermal_subproblems_list_ && global_probes_characteristics_)
    start_index_ = (int) (sub_problem_index_ * (*points_per_thermal_subproblem_));
  else
    {
      start_index_ = thermal_subproblems_rhs_size;
      thermal_subproblems_rhs_assembly.resize(thermal_subproblems_rhs_size + *points_per_thermal_subproblem_);
    }
  end_index_ = start_index_ + (*points_per_thermal_subproblem_);

  if (*first_time_step_temporal_)
    {
      if (first_time_step_explicit_)
        if (!(pre_initialise_thermal_subproblems_list_ && global_probes_characteristics_))
          thermal_subproblems_temperature_solution_ini.resize(thermal_subproblems_rhs_size + *points_per_thermal_subproblem_);
      temperature_ini_temporal_schemes_.resize(*points_per_thermal_subproblem_);
      if (boundary_condition_interface == dirichlet || boundary_condition_interface == default_bc)
        temperature_ini_temporal_schemes_[0] = interfacial_boundary_condition_value;
      else
        temperature_ini_temporal_schemes_[0] = temperature_interp_[0];
      double end_field_value_temporal=0.;
      if (boundary_condition_end == dirichlet || boundary_condition_end == default_bc)
        if (impose_user_boundary_condition_end_value)
          end_field_value_temporal = end_boundary_condition_value;
        else
          end_field_value_temporal = delta_T_subcooled_overheated_;
      else
        end_field_value_temporal = temperature_interp_[(*points_per_thermal_subproblem_) - 1];
      for (int i=1; i<(*points_per_thermal_subproblem_); i++)
        temperature_ini_temporal_schemes_[i] = end_field_value_temporal;
      // temperature_ini_temporal_schemes_[(*points_per_thermal_subproblem_) - 1] = delta_T_subcooled_overheated_;
    }
}

void IJK_One_Dimensional_Subproblem::compute_source_terms_impose_boundary_conditions(DoubleVect& thermal_subproblems_rhs_assembly,
                                                                                     DoubleVect& thermal_subproblems_temperature_solution_ini,
                                                                                     const int& boundary_condition_interface,
                                                                                     const double& interfacial_boundary_condition_value,
                                                                                     const int& impose_boundary_condition_interface_from_simulation,
                                                                                     const int& boundary_condition_end,
                                                                                     const double& end_boundary_condition_value,
                                                                                     const int& impose_user_boundary_condition_end_value)
{

  prepare_boundary_conditions(thermal_subproblems_rhs_assembly,
                              thermal_subproblems_temperature_solution_ini,
                              boundary_condition_interface,
                              interfacial_boundary_condition_value,
                              impose_boundary_condition_interface_from_simulation,
                              boundary_condition_end,
                              end_boundary_condition_value,
                              impose_user_boundary_condition_end_value);

  compute_add_source_terms();

  (*finite_difference_assembler_).impose_boundary_conditions_subproblem(thermal_subproblems_matrix_assembly_,
                                                                        thermal_subproblems_rhs_assembly_,
                                                                        rhs_assembly_,
                                                                        boundary_condition_interface,
                                                                        interfacial_boundary_condition_value_,
                                                                        boundary_condition_end,
                                                                        end_boundary_condition_value_,
                                                                        sub_problem_index_,
                                                                        dr_inv_,
                                                                        (*first_time_step_temporal_),
                                                                        first_time_step_explicit_,
                                                                        temperature_ini_temporal_schemes_,
                                                                        start_index_);

  if ((*first_time_step_temporal_) && first_time_step_explicit_)
    (*finite_difference_assembler_).add_source_terms(thermal_subproblems_temperature_solution_ini_, temperature_ini_temporal_schemes_);
}

void IJK_One_Dimensional_Subproblem::compute_add_source_terms()
{
  if ((source_terms_type_ != linear_diffusion
       && source_terms_type_ != spherical_diffusion
       && source_terms_type_ != spherical_diffusion_approx) || !avoid_post_processing_all_terms_)
    {
      interpolate_temperature_gradient_on_probe();
      project_temperature_gradient_on_probes();

      if (source_terms_type_ == tangential_conv_2D_tangential_diffusion_3D ||
          source_terms_type_ == tangential_conv_3D_tangentual_diffusion_3D ||
          !avoid_post_processing_all_terms_)
        {
          interpolate_temperature_hessian_on_probe();
          project_temperature_hessian_on_probes();
          retrieve_temperature_diffusion_spherical_on_probes();
        }
    }
  if (source_terms_type_ != linear_diffusion
      && source_terms_type_ != spherical_diffusion
      && source_terms_type_ != spherical_diffusion_approx)
    {
      /*
       * Compute at least the tangential convection
       */
      tangential_convection_source_terms_first_ = (*tangential_temperature_gradient_first_solver_);
      if (correct_tangential_temperature_gradient_)
        correct_tangential_temperature_gradient(tangential_convection_source_terms_first_);
      tangential_convection_source_terms_first_ *= (*first_tangential_velocity_solver_);
      const double alpha_inv = 1 / (*alpha_);
      tangential_convection_source_terms_first_ *= alpha_inv;
      switch (source_terms_type_)
        {
        case tangential_conv_2D:
          tangential_convection_source_terms_ = tangential_convection_source_terms_first_;
          source_terms_ = tangential_convection_source_terms_;
          break;
        case tangential_conv_3D:
          tangential_convection_source_terms_second_ = (*tangential_temperature_gradient_second_solver_);
          if (correct_tangential_temperature_gradient_)
            correct_tangential_temperature_gradient(tangential_convection_source_terms_second_);
          tangential_convection_source_terms_second_ *= (*second_tangential_velocity_solver_);
          tangential_convection_source_terms_second_ *= alpha_inv;
          tangential_convection_source_terms_ = tangential_convection_source_terms_first_;
          tangential_convection_source_terms_ += tangential_convection_source_terms_second_;
          source_terms_ = tangential_convection_source_terms_;
          break;
        case tangential_conv_2D_tangential_diffusion_3D:
          tangential_convection_source_terms_ = tangential_convection_source_terms_first_;
          source_terms_ = tangential_convection_source_terms_;
          tangential_diffusion_source_terms_ = tangential_temperature_diffusion_;
          if (correct_tangential_temperature_hessian_)
            correct_tangential_temperature_hessian(tangential_diffusion_source_terms_);
          // Be careful to the sign
          source_terms_ += tangential_diffusion_source_terms_;
          break;
        case tangential_conv_3D_tangentual_diffusion_3D:
          tangential_convection_source_terms_second_ = (*tangential_temperature_gradient_second_solver_);
          if (correct_tangential_temperature_gradient_)
            correct_tangential_temperature_gradient(tangential_convection_source_terms_second_);
          tangential_convection_source_terms_second_ *= second_tangential_velocity_;
          tangential_convection_source_terms_second_ *= alpha_inv;
          tangential_convection_source_terms_ = tangential_convection_source_terms_first_;
          tangential_convection_source_terms_ += tangential_convection_source_terms_second_;
          source_terms_ = tangential_convection_source_terms_;
          // tangential_diffusion_source_terms_;
          // if (correct_tangential_temperature_hessian_)
          // 	correct_tangential_temperature_hessian(tangential_diffusion_source_terms_)
          // Be careful to the sign
          source_terms_ += tangential_diffusion_source_terms_;
          break;
        default:
          tangential_convection_source_terms_ = tangential_convection_source_terms_first_;
          source_terms_ = tangential_convection_source_terms_;
          break;
        }

      if (*first_time_step_temporal_)
        {
          source_terms_ *= (local_time_step_overall_ * (*alpha_));
          rhs_assembly_ = source_terms_;
          if (first_time_step_explicit_)
            rhs_assembly_[0] = 0.;
          else
            rhs_assembly_ += temperature_ini_temporal_schemes_;
        }
      else
        rhs_assembly_ = source_terms_;

      (*finite_difference_assembler_).add_source_terms(thermal_subproblems_rhs_assembly_, rhs_assembly_);
    }
}

void IJK_One_Dimensional_Subproblem::approximate_temperature_increment_material_derivative()
{
  approximate_partial_temperature_time_increment();
  approximate_temperature_material_derivatives();
}

void IJK_One_Dimensional_Subproblem::approximate_partial_temperature_time_increment()
{
  temperature_time_increment_.resize(*points_per_thermal_subproblem_);
  if (!avoid_post_processing_all_terms_)
    {
      // const double dt_iter = 0.;
      const double alpha_liq = *alpha_;
      for (int i=0; i<*points_per_thermal_subproblem_; i++)
        {
          temperature_time_increment_[i] = -(x_velocity_(i) * grad_T_elem_interp_[0](i) +
                                             y_velocity_(i) * grad_T_elem_interp_[1](i) +
                                             z_velocity_(i) * grad_T_elem_interp_[2](i)) +
                                           alpha_liq *
                                           (hess_diag_T_elem_interp_[0](i) +
                                            hess_diag_T_elem_interp_[1](i) +
                                            hess_diag_T_elem_interp_[2](i));
        }
    }
}

void IJK_One_Dimensional_Subproblem::approximate_temperature_material_derivatives()
{
  /*
   * Two frame reference calculations
   * Advected or not by the velocity tangentially
   * -> 4 cases
   */
  if (!avoid_post_processing_all_terms_)
    {
      approximate_temperature_material_derivatives(normal_vector_compo_,
                                                   first_tangential_vector_compo_,
                                                   second_tangential_vector_compo_,
                                                   radial_velocity_advected_frame_,
                                                   first_tangential_velocity_advected_frame_,
                                                   second_tangential_velocity_advected_frame_,
                                                   temperature_time_increment_,
                                                   convective_term_advected_frame_,
                                                   material_derivative_advected_frame_);
      approximate_temperature_material_derivatives(normal_vector_compo_,
                                                   first_tangential_vector_compo_,
                                                   second_tangential_vector_compo_,
                                                   radial_velocity_static_frame_,
                                                   first_tangential_velocity_static_frame_,
                                                   second_tangential_velocity_static_frame_,
                                                   temperature_time_increment_,
                                                   convective_term_static_frame_,
                                                   material_derivative_static_frame_);
      approximate_temperature_material_derivatives(normal_vector_compo_,
                                                   first_tangential_vector_compo_from_rising_dir_,
                                                   azymuthal_vector_compo_,
                                                   radial_velocity_advected_frame_,
                                                   first_tangential_velocity_advected_frame_,
                                                   second_tangential_velocity_advected_frame_,
                                                   temperature_time_increment_,
                                                   convective_term_advected_frame_rising_,
                                                   material_derivative_advected_frame_rising_);
      approximate_temperature_material_derivatives(normal_vector_compo_,
                                                   first_tangential_vector_compo_from_rising_dir_,
                                                   azymuthal_vector_compo_,
                                                   radial_velocity_static_frame_,
                                                   first_tangential_velocity_static_frame_,
                                                   second_tangential_velocity_static_frame_,
                                                   temperature_time_increment_,
                                                   convective_term_static_frame_rising_,
                                                   material_derivative_static_frame_rising_);
    }
}

void IJK_One_Dimensional_Subproblem::approximate_temperature_material_derivatives(const Vecteur3& normal_vector_compo,
                                                                                  const Vecteur3& first_tangential_vector_compo,
                                                                                  const Vecteur3& second_tangential_vector_compo,
                                                                                  const DoubleVect& radial_velocity_frame,
                                                                                  const DoubleVect& first_tangential_velocity_frame,
                                                                                  const DoubleVect& second_tangential_velocity_frame,
                                                                                  const DoubleVect& temperature_time_increment,
                                                                                  DoubleVect& convective_term,
                                                                                  DoubleVect& material_derivative)
{
  material_derivative.resize(*points_per_thermal_subproblem_);
  convective_term.resize(*points_per_thermal_subproblem_);
  for (int i=0; i<*points_per_thermal_subproblem_; i++)
    {
      Vecteur3 cartesian_velocity = {x_velocity_(i), y_velocity_(i), z_velocity_(i)};

      const double radial_velocity = radial_velocity_frame[i];
      const double first_tangential_velocity = first_tangential_velocity_frame[i];
      const double second_tangential_velocity = second_tangential_velocity_frame[i];
      Vecteur3 temperature_gradient = {grad_T_elem_interp_[0](i), grad_T_elem_interp_[1](i), grad_T_elem_interp_[2](i)};
      Vecteur3 radial_convection_frame = normal_vector_compo;
      radial_convection_frame *= -radial_velocity;
      Vecteur3 first_tangent_convection_frame = first_tangential_vector_compo;
      first_tangent_convection_frame *= -first_tangential_velocity;
      Vecteur3 second_tangent_convection_frame = second_tangential_vector_compo;
      second_tangent_convection_frame *= -second_tangential_velocity;
      radial_convection_frame += cartesian_velocity;
      first_tangent_convection_frame += cartesian_velocity;
      second_tangent_convection_frame += cartesian_velocity;
      convective_term[i] = Vecteur3::produit_scalaire(radial_convection_frame, temperature_gradient) +
                           Vecteur3::produit_scalaire(first_tangent_convection_frame, temperature_gradient) +
                           Vecteur3::produit_scalaire(second_tangent_convection_frame, temperature_gradient);
      material_derivative[i] = temperature_time_increment[i] + convective_term[i];
    }
}

void IJK_One_Dimensional_Subproblem::correct_tangential_temperature_gradient(DoubleVect& tangential_convection_source_terms)
{
  DoubleVect tangential_convection_source_terms_tmp = tangential_convection_source_terms;
  for (int i=0; i<tangential_convection_source_terms.size(); i++)
    tangential_convection_source_terms[i] -= tangential_convection_source_terms_tmp[0];
}

void IJK_One_Dimensional_Subproblem::correct_tangential_temperature_hessian(DoubleVect& tangential_diffusion_source_terms)
{
  DoubleVect tangential_diffusion_source_terms_tmp = tangential_diffusion_source_terms_tmp;
  for (int i=0; i<tangential_diffusion_source_terms_tmp.size(); i++)
    tangential_diffusion_source_terms_tmp[i] -= tangential_diffusion_source_terms_tmp[0];
}

void IJK_One_Dimensional_Subproblem::retrieve_radial_quantities()
{
  retrieve_temperature_solution();
  compute_local_temperature_gradient_solution();
  compute_local_velocity_gradient();
  compute_local_pressure_gradient();
  is_updated_ = true;
}

void IJK_One_Dimensional_Subproblem::retrieve_temperature_solution()
{
  temperature_solution_.resize(rhs_assembly_.size());
  for (int i=start_index_; i<end_index_; i++)
    temperature_solution_[i - start_index_] = (*thermal_subproblems_temperature_solution_)[i];
}

void IJK_One_Dimensional_Subproblem::compute_local_temperature_gradient_solution()
{
  //	for (int dir=0; dir<3; dir++)
  //		temperature_gradient_solution_[dir].resize(temperature_solution_.size());
  normal_temperature_gradient_solution_.resize(temperature_solution_.size());
  (*finite_difference_assembler_).compute_operator(radial_first_order_operator_, temperature_solution_, normal_temperature_gradient_solution_);

  normal_temperature_double_derivative_solution_.resize(temperature_solution_.size());
  (*finite_difference_assembler_).compute_operator(radial_second_order_operator_, temperature_solution_, normal_temperature_double_derivative_solution_);

  reinit_variable(temperature_x_gradient_solution_);
  reinit_variable(temperature_y_gradient_solution_);
  reinit_variable(temperature_z_gradient_solution_);
  if (source_terms_type_ == linear_diffusion
      || source_terms_type_ == spherical_diffusion
      || source_terms_type_ == spherical_diffusion_approx)
    {
      DoubleVect dummy_tangential_deriv;
      dummy_tangential_deriv.resize(*points_per_thermal_subproblem_);
      temperature_x_gradient_solution_.resize(normal_temperature_gradient_solution_.size());
      project_basis_vector_onto_cartesian_dir(0, dummy_tangential_deriv, dummy_tangential_deriv, normal_temperature_gradient_solution_,
                                              *first_tangential_vector_compo_solver_, *second_tangential_vector_compo_solver_, normal_vector_compo_,
                                              temperature_x_gradient_solution_);
      temperature_y_gradient_solution_.resize(normal_temperature_gradient_solution_.size());
      project_basis_vector_onto_cartesian_dir(1, dummy_tangential_deriv, dummy_tangential_deriv, normal_temperature_gradient_solution_,
                                              *first_tangential_vector_compo_solver_, *second_tangential_vector_compo_solver_, normal_vector_compo_,
                                              temperature_y_gradient_solution_);
      temperature_z_gradient_solution_.resize(normal_temperature_gradient_solution_.size());
      project_basis_vector_onto_cartesian_dir(2, dummy_tangential_deriv, dummy_tangential_deriv, normal_temperature_gradient_solution_,
                                              *first_tangential_vector_compo_solver_, *second_tangential_vector_compo_solver_, normal_vector_compo_,
                                              temperature_z_gradient_solution_);
    }
  else
    {
      temperature_x_gradient_solution_.resize(normal_temperature_gradient_solution_.size());
      project_basis_vector_onto_cartesian_dir(0, tangential_temperature_gradient_first_, tangential_temperature_gradient_second_, normal_temperature_gradient_solution_,
                                              *first_tangential_vector_compo_solver_, *second_tangential_vector_compo_solver_, normal_vector_compo_,
                                              temperature_x_gradient_solution_);
      temperature_y_gradient_solution_.resize(normal_temperature_gradient_solution_.size());
      project_basis_vector_onto_cartesian_dir(1, tangential_temperature_gradient_first_, tangential_temperature_gradient_second_, normal_temperature_gradient_solution_,
                                              *first_tangential_vector_compo_solver_, *second_tangential_vector_compo_solver_, normal_vector_compo_,
                                              temperature_y_gradient_solution_);
      temperature_z_gradient_solution_.resize(normal_temperature_gradient_solution_.size());
      project_basis_vector_onto_cartesian_dir(2, tangential_temperature_gradient_first_, tangential_temperature_gradient_second_, normal_temperature_gradient_solution_,
                                              *first_tangential_vector_compo_solver_, *second_tangential_vector_compo_solver_, normal_vector_compo_,
                                              temperature_z_gradient_solution_);
    }

  if ((source_terms_type_ == linear_diffusion
       || source_terms_type_ == spherical_diffusion
       || source_terms_type_ == spherical_diffusion_approx) && avoid_post_processing_all_terms_)
    {
      normal_temperature_gradient_.resize(*points_per_thermal_subproblem_);
      tangential_temperature_gradient_first_.resize(*points_per_thermal_subproblem_);
      tangential_temperature_gradient_second_.resize(*points_per_thermal_subproblem_);
      tangential_temperature_gradient_first_from_rising_dir_.resize(*points_per_thermal_subproblem_);
      azymuthal_temperature_gradient_.resize(*points_per_thermal_subproblem_);
      for (int dir = 0; dir < 3; dir++)
        {
          hess_diag_T_elem_interp_[dir].resize(*points_per_thermal_subproblem_);
          hess_cross_T_elem_interp_[dir].resize(*points_per_thermal_subproblem_);
          hess_diag_T_elem_spherical_[dir].resize(*points_per_thermal_subproblem_);
          hess_cross_T_elem_spherical_[dir].resize(*points_per_thermal_subproblem_);
          hess_diag_T_elem_spherical_from_rising_[dir].resize(*points_per_thermal_subproblem_);
          hess_cross_T_elem_spherical_from_rising_[dir].resize(*points_per_thermal_subproblem_);
        }
      temperature_diffusion_hessian_trace_.resize(*points_per_thermal_subproblem_);
      temperature_diffusion_hessian_cartesian_trace_.resize(*points_per_thermal_subproblem_);
      radial_temperature_diffusion_.resize(*points_per_thermal_subproblem_);
      tangential_temperature_diffusion_.resize(*points_per_thermal_subproblem_);
      tangential_convection_source_terms_first_.resize(*points_per_thermal_subproblem_);
      tangential_convection_source_terms_second_.resize(*points_per_thermal_subproblem_);
      tangential_diffusion_source_terms_.resize(*points_per_thermal_subproblem_);
      source_terms_.resize(*points_per_thermal_subproblem_);
    }

  thermal_flux_ = normal_temperature_gradient_solution_;
  thermal_flux_*= ((*lambda_) * surface_);
}

double IJK_One_Dimensional_Subproblem::get_interfacial_gradient_corrected() const
{
  return normal_temperature_gradient_solution_[0];
}

double IJK_One_Dimensional_Subproblem::get_interfacial_double_derivative_corrected() const
{
  return normal_temperature_double_derivative_solution_[0];
}

void IJK_One_Dimensional_Subproblem::compute_local_velocity_gradient()
{
  normal_velocity_normal_gradient_.resize(*points_per_thermal_subproblem_);
  first_tangential_velocity_normal_gradient_.resize(*points_per_thermal_subproblem_);
  second_tangential_velocity_normal_gradient_.resize(*points_per_thermal_subproblem_);
  azymuthal_velocity_normal_gradient_.resize(*points_per_thermal_subproblem_);
  first_tangential_velocity_gradient_from_rising_dir_.resize(*points_per_thermal_subproblem_);

  (*finite_difference_assembler_).compute_operator(radial_first_order_operator_,
                                                   radial_velocity_,
                                                   normal_velocity_normal_gradient_);
  (*finite_difference_assembler_).compute_operator(radial_first_order_operator_,
                                                   first_tangential_velocity_,
                                                   first_tangential_velocity_normal_gradient_);
  (*finite_difference_assembler_).compute_operator(radial_first_order_operator_,
                                                   second_tangential_velocity_,
                                                   second_tangential_velocity_normal_gradient_);
  (*finite_difference_assembler_).compute_operator(radial_first_order_operator_,
                                                   first_tangential_velocity_from_rising_dir_,
                                                   first_tangential_velocity_gradient_from_rising_dir_);
  (*finite_difference_assembler_).compute_operator(radial_first_order_operator_,
                                                   azymuthal_velocity_,
                                                   azymuthal_velocity_normal_gradient_);
  compute_local_shear_stress();
}

void IJK_One_Dimensional_Subproblem::compute_local_shear_stress()
{
  shear_stress_.resize(*points_per_thermal_subproblem_);
  shear_stress_from_rising_dir_.resize(*points_per_thermal_subproblem_);
  for (int i=0; i<(*points_per_thermal_subproblem_); i++)
    {
      Vecteur3 local_shear_stress = first_tangential_vector_compo_;
      Vecteur3 local_shear_stress_second = second_tangential_vector_compo_;
      local_shear_stress *= first_tangential_velocity_normal_gradient_(i);
      local_shear_stress_second *= second_tangential_velocity_normal_gradient_(i);
      local_shear_stress += local_shear_stress_second;
      shear_stress_(i) = local_shear_stress.length();

      local_shear_stress = first_tangential_vector_compo_from_rising_dir_;
      local_shear_stress_second = azymuthal_vector_compo_;
      local_shear_stress *= first_tangential_velocity_gradient_from_rising_dir_(i);
      local_shear_stress_second *= azymuthal_velocity_normal_gradient_(i);
      local_shear_stress += local_shear_stress_second;
      shear_stress_from_rising_dir_(i) = local_shear_stress.length();
    }
  const double mu_liquid = ref_ijk_ft_->get_mu_liquid();
  shear_stress_ *= mu_liquid;
  shear_stress_from_rising_dir_ *= mu_liquid;
  velocity_shear_stress_ = shear_stress_[0];
  velocity_shear_force_ = velocity_shear_stress_ * surface_;
}

void IJK_One_Dimensional_Subproblem::compute_local_pressure_gradient()
{
  pressure_normal_gradient_.resize(*points_per_thermal_subproblem_);
  (*finite_difference_assembler_).compute_operator(radial_first_order_operator_,
                                                   pressure_interp_,
                                                   pressure_normal_gradient_);
  pressure_gradient_ = pressure_normal_gradient_[0];
}

double IJK_One_Dimensional_Subproblem::get_normal_velocity_normal_gradient() const
{
  return normal_velocity_normal_gradient_[0];
}

double IJK_One_Dimensional_Subproblem::get_tangential_velocity_normal_gradient() const
{
  return first_tangential_velocity_normal_gradient_[0];
}

double IJK_One_Dimensional_Subproblem::get_second_tangential_velocity_normal_gradient() const
{
  return second_tangential_velocity_normal_gradient_[0];
}

double IJK_One_Dimensional_Subproblem::get_azymuthal_velocity_normal_gradient() const
{
  return azymuthal_velocity_normal_gradient_[0];
}

double IJK_One_Dimensional_Subproblem::get_field_profile_at_point(const double& dist,
                                                                  const DoubleVect& field,
                                                                  const IJK_Field_double& eulerian_field,
                                                                  const int temp_bool,
                                                                  const int interp_eulerian) const
{
  double field_value = INVALID_TEMPERATURE;
  if (dist >= (*radial_coordinates_)[0] && dist <= (*radial_coordinates_)[*points_per_thermal_subproblem_-1])
    {
      /*
       * Dummy dichotomy and linear interpolation along the probe
       */
      int left_interval = 0;
      int right_interval = *points_per_thermal_subproblem_-1;
      find_interval(dist, left_interval, right_interval);
      const double field_interp = (field[right_interval] - field[left_interval])
                                  / ((*radial_coordinates_)[right_interval] - (*radial_coordinates_)[left_interval]) *
                                  (dist-(*radial_coordinates_)[left_interval]) + field[left_interval];
      field_value = field_interp;
    }
  else if (dist < (*radial_coordinates_)[0])
    {
      if (temp_bool)
        {
          const double interfacial_temperature_gradient_solution = get_interfacial_gradient_corrected();
          const double interfacial_temperature_double_derivative_solution = get_interfacial_double_derivative_corrected();
          switch(order_approx_temperature_ext_)
            {
            case 1:
              field_value = interfacial_temperature_gradient_solution * dist;
              break;
            case 2:
              field_value = (interfacial_temperature_gradient_solution * dist
                             + 0.5 * interfacial_temperature_double_derivative_solution * (dist * dist));
              break;
            default:
              field_value = interfacial_temperature_gradient_solution * dist;
            }
        }
      else
        field_value = field[0];
    }
  else
    {
      if(interp_eulerian)
        {
          DoubleTab coordinates_point;
          DoubleVect field_interp(1);
          coordinates_point.resize(1,3);
          Vecteur3 compo_xyz = normal_vector_compo_;
          compo_xyz *= dist;
          compo_xyz += facet_barycentre_;
          for (int c=0; c<3; c++)
            coordinates_point(0,c) = compo_xyz[0];
          ijk_interpolate_skip_unknown_points(eulerian_field, coordinates_point, field_interp, INVALID_INTERP);
          field_value = field_interp(0);
        }
      else
        field_value = field[*points_per_thermal_subproblem_ - 1];
    }
  return field_value;
}

double IJK_One_Dimensional_Subproblem::get_field_profile_at_point(const double& dist, const DoubleVect& field, const int temp_bool) const
{
  double field_value = INVALID_TEMPERATURE;// temperature_solution_;
  if (debug_ && ((dist < 0 && indicator_>0.5) || (dist > (*radial_coordinates_)[*points_per_thermal_subproblem_-1] && indicator_>0.5)))
    {
      Cerr << "Probe length: " << probe_length_ << finl;
      Cerr << "Distance d: " << dist << finl;
      Cerr << "Indicator I: " << indicator_ << finl;
    }
  if (debug_ && temp_bool)
    {
      Cerr << "Radial coordinate ini: " << (*radial_coordinates_)[0] << finl;
      Cerr << "Radial coordinate end: " << (*radial_coordinates_)[*points_per_thermal_subproblem_-1] << finl;
      Cerr << "Field ini: " << field[0] << finl;
      Cerr << "Field end: " << field[*points_per_thermal_subproblem_-1] << finl;
      Cerr << "Field interp ini: " << temperature_interp_[0] << finl;
      Cerr << "Field interp end: " << temperature_interp_[*points_per_thermal_subproblem_-1] << finl;
      Cerr << "End_boundary_condition_value: " << end_boundary_condition_value_ << finl;
      Cerr << "Curvature: " << curvature_ << finl;
      Cerr << "Osculating radius: " << osculating_radius_ << finl;
    }
  if (dist >= (*radial_coordinates_)[0] && dist <= (*radial_coordinates_)[*points_per_thermal_subproblem_-1])
    {
      /*
       * Dummy dichotomy and linear interpolation along the probe
       */
      int left_interval = 0;
      int right_interval = *points_per_thermal_subproblem_-1;
      find_interval(dist, left_interval, right_interval);
      const double field_interp = (field[right_interval] - field[left_interval])
                                  / ((*radial_coordinates_)[right_interval] - (*radial_coordinates_)[left_interval]) *
                                  (dist-(*radial_coordinates_)[left_interval]) + field[left_interval];
      field_value = field_interp;
    }
  else if (dist < (*radial_coordinates_)[0])
    {
      if (temp_bool)
        {
          if (short_probe_condition_)
            field_value = 0.;
          else
            {
              const double interfacial_temperature_gradient_solution = get_interfacial_gradient_corrected();
              const double interfacial_temperature_double_derivative_solution = get_interfacial_double_derivative_corrected();
              switch(order_approx_temperature_ext_)
                {
                case 1:
                  field_value = interfacial_temperature_gradient_solution * dist;
                  break;
                case 2:
                  field_value = (interfacial_temperature_gradient_solution * dist
                                 + 0.5 * interfacial_temperature_double_derivative_solution * (dist * dist));
                  break;
                default:
                  field_value = interfacial_temperature_gradient_solution * dist;
                }
            }
        }
      else
        field_value = 0.;
    }
  else
    {
      if (temp_bool)
        {
          if (short_probe_condition_)
            {
              if (temperature_probe_condition_)
                field_value = cell_temperature_;
              else
                field_value = delta_T_subcooled_overheated_;
            }
          else
            field_value = field[(*points_per_thermal_subproblem_) - 1];
        }
      else
        field_value = field[(*points_per_thermal_subproblem_) - 1];
    }
  return field_value;
}

double IJK_One_Dimensional_Subproblem::get_temperature_profile_at_point(const double& dist) const
{
  return get_field_profile_at_point(dist, temperature_solution_, *temperature_, 1, interp_eulerian_);
  // return get_field_profile_at_point(dist, temperature_solution_, 1);
}

double IJK_One_Dimensional_Subproblem::get_velocity_component_at_point(const double& dist, const int& dir) const
{
  double velocity = 0;
  switch(dir)
    {
    case 0:
      velocity = get_field_profile_at_point(dist, x_velocity_, (*velocity_)[0] , 0, interp_eulerian_);
      // velocity = get_field_profile_at_point(dist, x_velocity_, 0);
      break;
    case 1:
      velocity = get_field_profile_at_point(dist, y_velocity_, (*velocity_)[1] , 0, interp_eulerian_);
//      velocity = get_field_profile_at_point(dist, y_velocity_, 0);
      break;
    case 2:
      velocity = get_field_profile_at_point(dist, z_velocity_, (*velocity_)[2] , 0, interp_eulerian_);
//      velocity = get_field_profile_at_point(dist, z_velocity_, 0);
      break;
    default:
      velocity = get_field_profile_at_point(dist, x_velocity_, (*velocity_)[0] , 0, interp_eulerian_);
      // velocity = get_field_profile_at_point(dist, x_velocity_, 0);
      break;
    }
  return velocity;
}

double IJK_One_Dimensional_Subproblem::get_temperature_gradient_profile_at_point(const double& dist, const int& dir) const
{
  double temperature_gradient = 0;
  switch(dir)
    {
    case 0:
      temperature_gradient = get_field_profile_at_point(dist, temperature_x_gradient_solution_, (*grad_T_elem_)[0], 0, interp_eulerian_);
//      temperature_gradient = get_field_profile_at_point(dist, temperature_x_gradient_solution_, 0);
      break;
    case 1:
      temperature_gradient = get_field_profile_at_point(dist, temperature_y_gradient_solution_, (*grad_T_elem_)[1], 0, interp_eulerian_);
//      temperature_gradient = get_field_profile_at_point(dist, temperature_y_gradient_solution_, 0);
      break;
    case 2:
      temperature_gradient = get_field_profile_at_point(dist, temperature_z_gradient_solution_, (*grad_T_elem_)[2], 0, interp_eulerian_);
//      temperature_gradient = get_field_profile_at_point(dist, temperature_z_gradient_solution_, 0);
      break;
    default:
      temperature_gradient = get_field_profile_at_point(dist, temperature_x_gradient_solution_, (*grad_T_elem_)[0], 0, interp_eulerian_);
//      temperature_gradient = get_field_profile_at_point(dist, temperature_x_gradient_solution_, 0);
      break;
    }
  return temperature_gradient;
}

double IJK_One_Dimensional_Subproblem::get_temperature_times_velocity_profile_at_point(const double& dist, const int& dir) const
{
//  double temperature_interp = get_field_profile_at_point(dist, temperature_solution_, 1);
  double temperature_interp = get_field_profile_at_point(dist, temperature_solution_, *temperature_, 1, interp_eulerian_);
  double velocity_interp = get_velocity_component_at_point(dist, dir);
  return temperature_interp * velocity_interp;
}

double IJK_One_Dimensional_Subproblem::get_temperature_gradient_times_conductivity_profile_at_point(const double& dist, const int& dir) const
{
  double diffusive_flux = 0;
  diffusive_flux = get_temperature_gradient_profile_at_point(dist, dir);
  diffusive_flux *= (*lambda_);
  return diffusive_flux;
}

DoubleVect IJK_One_Dimensional_Subproblem::get_field_discrete_integral_velocity_weighting_at_point(const double& dist,
                                                                                                   const int& levels,
                                                                                                   const int& dir,
                                                                                                   const DoubleVect& field,
                                                                                                   const IJK_Field_double& eulerian_field,
                                                                                                   const int temp_bool,
                                                                                                   const int vel) const
{
  const int nb_values = (int) pow(4., (double) levels);
  DoubleVect discrete_values(nb_values);
  double surface = get_discrete_surface_at_level(dir, levels);
  double value;
  double velocity;
  int value_counter = 0;
  if (levels==0)
    {
      velocity = get_velocity_weighting(dist, dir, vel);
      value = get_field_profile_at_point(dist, field, eulerian_field, temp_bool, interp_eulerian_);
      discrete_values(0) = value * surface * velocity;
    }
  else
    {
      double dl1_ini = 0.;
      double dl2_ini = 0.;
      Vecteur3 point_coords_ini = {0., 0., 0.};
      get_field_discrete_value_recursive(1, levels + 1, dir, dist, vel, surface, field, dl1_ini, dl2_ini, point_coords_ini, discrete_values, value_counter);
    }
  return discrete_values;
}

void IJK_One_Dimensional_Subproblem::get_field_discrete_value_recursive(const int& ilevel, const int& max_level,
                                                                        const int& dir, const double& dist,
                                                                        const int& vel,
                                                                        const double& surface,
                                                                        const DoubleVect& field,
                                                                        const double dl1_parent,
                                                                        const double dl2_parent,
                                                                        Vecteur3& point_coords_parent,
                                                                        DoubleVect& discrete_values,
                                                                        int& value_counter) const
{
  if (ilevel != max_level)
    {
      const double neighbours_first_dir[4] = NEIGHBOURS_FIRST_DIR;
      const double neighbours_second_dir[4] = NEIGHBOURS_SECOND_DIR;
      for(int i=ilevel; i<max_level; i++)
        {
          if (i==ilevel)
            {
              for(int l=0; l<4; l++)
                {
                  const double first_dir = neighbours_first_dir[l];
                  const double second_dir = neighbours_second_dir[l];
                  double dl1;
                  double dl2;
                  Vecteur3 point_coords = {0., 0., 0.};
                  get_discrete_two_dimensional_spacing(dir, ilevel, first_dir, second_dir, dl1, dl2, point_coords);
                  dl1 += dl1_parent;
                  dl2 += dl2_parent;
                  point_coords += point_coords_parent;
                  get_field_discrete_value_recursive(i+1, max_level, dir, dist, vel, surface, field, dl1, dl2, point_coords, discrete_values, value_counter);
                }
            }
        }
    }
  else
    {
      const double dist_increment = Vecteur3::produit_scalaire(point_coords_parent, normal_vector_compo_);
      const double dist_value = dist + dist_increment;
      const double value = get_field_profile_at_point(dist_value, field, 0);
      const double velocity = get_velocity_weighting(dist, dir, vel);
      discrete_values(value_counter) = value * surface * velocity;
      value_counter++;
    }
}

double IJK_One_Dimensional_Subproblem::get_velocity_weighting(const double& dist, const int& dir, const int vel) const
{
  if (vel)
    return get_velocity_component_at_point(dist, dir);
  else
    return 1.;
}

DoubleVect IJK_One_Dimensional_Subproblem::get_field_discrete_integral_at_point(const double& dist, const int& levels,
                                                                                const int& dir, const DoubleVect& field,
                                                                                const IJK_Field_double& eulerian_field,
                                                                                const int temp_bool) const
{
  return get_field_discrete_integral_velocity_weighting_at_point(dist, levels, dir, field, eulerian_field, temp_bool, 0);
}

DoubleVect IJK_One_Dimensional_Subproblem::get_field_times_velocity_discrete_integral_at_point(const double& dist, const int& levels,
                                                                                               const int& dir, const DoubleVect& field,
                                                                                               const IJK_Field_double& eulerian_field) const
{
  return get_field_discrete_integral_velocity_weighting_at_point(dist, levels, dir, field, eulerian_field, 1, 1);
}

DoubleVect IJK_One_Dimensional_Subproblem::get_temperature_profile_discrete_integral_at_point(const double& dist,
                                                                                              const int& levels,
                                                                                              const int& dir) const
{
  return get_field_discrete_integral_at_point(dist, levels, dir, temperature_solution_, *temperature_, 1);
}

DoubleVect IJK_One_Dimensional_Subproblem::get_temperature_times_velocity_profile_discrete_integral_at_point(const double& dist,
                                                                                                             const int& levels,
                                                                                                             const int& dir) const
{
  return get_field_times_velocity_discrete_integral_at_point(dist, levels, dir, temperature_solution_, *temperature_);
}

DoubleVect IJK_One_Dimensional_Subproblem::get_temperature_gradient_profile_discrete_integral_at_point(const double& dist, const int& levels, const int& dir) const
{
  DoubleVect temperature_gradient;
  switch(dir)
    {
    case 0:
      temperature_gradient = get_field_discrete_integral_at_point(dist, levels, dir, temperature_x_gradient_solution_, (*grad_T_elem_)[0], 0);
      // temperature_gradient = get_field_discrete_integral_at_point(dist, levels, dir, temperature_x_gradient_solution_);
      break;
    case 1:
      temperature_gradient = get_field_discrete_integral_at_point(dist, levels, dir, temperature_y_gradient_solution_, (*grad_T_elem_)[1], 0);
//      temperature_gradient = get_field_discrete_integral_at_point(dist, levels, dir, temperature_y_gradient_solution_);
      break;
    case 2:
      temperature_gradient = get_field_discrete_integral_at_point(dist, levels, dir, temperature_z_gradient_solution_, (*grad_T_elem_)[2], 0);
//      temperature_gradient = get_field_discrete_integral_at_point(dist, levels, dir, temperature_z_gradient_solution_);
      break;
    default:
      temperature_gradient = get_field_discrete_integral_at_point(dist, levels, dir, temperature_x_gradient_solution_, (*grad_T_elem_)[0], 0);
//      temperature_gradient = get_field_discrete_integral_at_point(dist, levels, dir, temperature_x_gradient_solution_);
      break;
    }
  return temperature_gradient;
}

DoubleVect IJK_One_Dimensional_Subproblem::get_temperature_gradient_times_conductivity_profile_discrete_integral_at_point(const double& dist, const int& levels, const int& dir) const
{
  DoubleVect diffusive_flux = get_temperature_gradient_profile_discrete_integral_at_point(dist, levels, dir);
  diffusive_flux *= (*lambda_);
  return diffusive_flux;
}

void IJK_One_Dimensional_Subproblem::find_interval(const double& dist, int& left_interval, int& right_interval) const
{
  int mid_interval = left_interval + (right_interval - left_interval) / 2;
  while ((right_interval - left_interval) != 1)
    {
      if (dist > (*radial_coordinates_)[mid_interval])
        left_interval = mid_interval;
      else
        right_interval = mid_interval;
      mid_interval = left_interval + (right_interval - left_interval) / 2;
    }
}

void IJK_One_Dimensional_Subproblem::get_discrete_two_dimensional_spacing(const int& dir, const int& levels,
                                                                          const double& first_dir, const double& second_dir,
                                                                          double& dl1, double& dl2,
                                                                          Vecteur3& point_coords) const
{
  const IJK_Grid_Geometry& geom = ref_ijk_ft_->get_splitting_ns().get_grid_geometry();
  const double dx = geom.get_constant_delta(DIRECTION_I);
  const double dy = geom.get_constant_delta(DIRECTION_J);
  const double dz = geom.get_constant_delta(DIRECTION_K);
  switch(dir)
    {
    case 0:
      dl1 = dy;
      dl2 = dz;
      point_coords[1] = dl1 * first_dir;
      point_coords[2] = dl2 * second_dir;
      break;
    case 1:
      dl1 = dx;
      dl2 = dz;
      point_coords[0] = dl1 * first_dir;
      point_coords[2] = dl2 * second_dir;
      break;
    case 2:
      dl1 = dx;
      dl2 = dy;
      point_coords[0] = dl1 * first_dir;
      point_coords[1] = dl2 * second_dir;
      break;
    default:
      dl1 = dx;
      dl2 = dy;
      point_coords[0] = dl1 * first_dir;
      point_coords[1] = dl2 * second_dir;
      break;
    }
  dl1 /= pow(2., (double) levels + 1.);
  dl2 /= pow(2., (double) levels + 1.);
  dl1 *= first_dir;
  dl2 *= second_dir;
  point_coords *= (1 / pow(2., (double) levels + 1.));
}

double IJK_One_Dimensional_Subproblem::get_discrete_surface_at_level(const int& dir, const int& level) const
{
  const IJK_Grid_Geometry& geom = ref_ijk_ft_->get_splitting_ns().get_grid_geometry();
  const double dx = geom.get_constant_delta(DIRECTION_I);
  const double dy = geom.get_constant_delta(DIRECTION_J);
  const double dz = geom.get_constant_delta(DIRECTION_K);
  double surface = 0.;
  switch(dir)
    {
    case 0:
      surface = (dy*dz);
      break;
    case 1:
      surface = (dx*dz);
      break;
    case 2:
      surface = (dx*dy);
      break;
    default:
      surface = (dx*dy);
      break;
    }
  surface /= pow(pow(2., (double) level), 2.);
  return surface;
}

void IJK_One_Dimensional_Subproblem::compute_temperature_integral_subproblem_probe()
{
  temperature_integral_ = compute_temperature_integral_subproblem(probe_length_);
}

double IJK_One_Dimensional_Subproblem::compute_temperature_integral_subproblem(const double& distance)
{
  double max_distance = distance;
  if (distance <= 0 || distance > probe_length_)
    max_distance = probe_length_;
  ArrOfDouble discrete_int_eval(*points_per_thermal_subproblem_ - 1);
  double integral_eval = 0.;
  const double radial_incr = max_distance / (*points_per_thermal_subproblem_ - 1);
  for (int i=0; i<(*points_per_thermal_subproblem_) - 1; i++)
    {
      discrete_int_eval(i) = get_temperature_profile_at_point(radial_incr * i);
      discrete_int_eval(i) += get_temperature_profile_at_point(radial_incr * (i + 1));
      discrete_int_eval(i) *= (radial_incr / 2.);
      integral_eval += discrete_int_eval(i);
    }
  integral_eval *= (1 / (radial_incr + 1e-20));
  return integral_eval;
}


void IJK_One_Dimensional_Subproblem::thermal_subresolution_outputs(SFichier& fic, const int rank, const Nom& local_quantities_thermal_probes_time_index_folder)
{
  post_process_interfacial_quantities(fic, rank);
  post_process_radial_quantities(rank, local_quantities_thermal_probes_time_index_folder);
}

void IJK_One_Dimensional_Subproblem::thermal_subresolution_outputs_parallel(const int rank, const Nom& local_quantities_thermal_probes_time_index_folder)
{
  post_process_radial_quantities(rank, local_quantities_thermal_probes_time_index_folder);
}

void IJK_One_Dimensional_Subproblem::retrieve_interfacial_quantities(const int rank,
                                                                     const int& itr,
                                                                     std::vector<std::string> key_results_int,
                                                                     std::vector<std::string> key_results_double,
                                                                     std::map<std::string, ArrOfInt>& results_probes_int,
                                                                     std::map<std::string, ArrOfDouble>& results_probes_double)
{
  const double last_time = ref_ijk_ft_->get_current_time() - ref_ijk_ft_->get_timestep();
  const int last_time_index = ref_ijk_ft_->get_tstep();
  std::vector<int> results_int =
  {
    last_time_index, rank, index_post_processing_, global_subproblem_index_, sub_problem_index_
  };

  std::vector<double> results_double =
  {
    last_time,
    normal_vector_compo_[0], normal_vector_compo_[1], normal_vector_compo_[2],
    first_tangential_vector_compo_[0], first_tangential_vector_compo_[1], first_tangential_vector_compo_[2],
    second_tangential_vector_compo_[0], second_tangential_vector_compo_[1], second_tangential_vector_compo_[2],
    first_tangential_vector_compo_from_rising_dir_[0], first_tangential_vector_compo_from_rising_dir_[1], first_tangential_vector_compo_from_rising_dir_[2],
    azymuthal_vector_compo_[0], azymuthal_vector_compo_[1], azymuthal_vector_compo_[2],
    r_sph_, theta_sph_, phi_sph_,
    temperature_interp_[0], temperature_solution_[0],
    normal_temperature_gradient_[0], normal_temperature_gradient_solution_[0],
    normal_temperature_double_derivative_solution_[0],
    tangential_temperature_gradient_first_[0],
    tangential_temperature_gradient_second_[0],
    tangential_temperature_gradient_first_from_rising_dir_[0],
    azymuthal_temperature_gradient_[0],
    temperature_diffusion_hessian_cartesian_trace_[0],
    temperature_diffusion_hessian_trace_[0],
    radial_temperature_diffusion_[0],
    tangential_temperature_diffusion_[0],
    surface_, thermal_flux_[0], (*lambda_), (*alpha_), Pr_l_,
    velocity_shear_force_, velocity_shear_stress_,
    pressure_interp_[0],
    x_velocity_[0], y_velocity_[0], z_velocity_[0],
    radial_velocity_[0], radial_velocity_corrected_[0],
    radial_velocity_static_frame_[0], radial_velocity_advected_frame_[0],
    first_tangential_velocity_[0], first_tangential_velocity_corrected_[0],
    first_tangential_velocity_static_frame_[0], first_tangential_velocity_advected_frame_[0],
    second_tangential_velocity_[0], second_tangential_velocity_corrected_[0],
    second_tangential_velocity_static_frame_[0], second_tangential_velocity_advected_frame_[0],
    first_tangential_velocity_from_rising_dir_[0], first_tangential_velocity_from_rising_dir_corrected_[0],
    first_tangential_velocity_from_rising_dir_static_frame_[0], first_tangential_velocity_from_rising_dir_advected_frame_[0],
    azymuthal_velocity_[0], azymuthal_velocity_corrected_[0],
    azymuthal_velocity_static_frame_[0], azymuthal_velocity_advected_frame_[0],
    normal_velocity_normal_gradient_[0],
    first_tangential_velocity_normal_gradient_[0],
    second_tangential_velocity_normal_gradient_[0],
    first_tangential_velocity_gradient_from_rising_dir_[0],
    azymuthal_velocity_normal_gradient_[0]
  };
  int i;
  assert(key_results_int.size() == results_int.size());
  int size_int = (int) key_results_int.size();
  for (i=0; i<size_int; i++)
    results_probes_int[key_results_int[i]](itr) = results_int[i];
  assert(key_results_double.size() == results_double.size());
  int size_double = (int) key_results_double.size();
  for (i=0; i<size_double; i++)
    results_probes_double[key_results_double[i]](itr) = results_double[i];

}

void IJK_One_Dimensional_Subproblem::post_process_interfacial_quantities(SFichier& fic, const int rank) //SFichier& fic)
{
  // if (Process::je_suis_maitre())
  {
    if (is_updated_)
      {
        const double last_time = ref_ijk_ft_->get_current_time() - ref_ijk_ft_->get_timestep();
        const int last_time_index = ref_ijk_ft_->get_tstep();
        fic << last_time_index << " ";
        fic << rank << " " << index_post_processing_ << " " << global_subproblem_index_ << " " << sub_problem_index_ << " ";
        fic << last_time << " ";
        fic << normal_vector_compo_[0] << " " << normal_vector_compo_[1] << " " << normal_vector_compo_[2] << " ";
        fic << first_tangential_vector_compo_[0] << " " << first_tangential_vector_compo_[1] << " " << first_tangential_vector_compo_[2] << " ";
        fic << second_tangential_vector_compo_[0] << " " << second_tangential_vector_compo_[1] << " " << second_tangential_vector_compo_[2] << " ";
        fic << first_tangential_vector_compo_from_rising_dir_[0] << " " << first_tangential_vector_compo_from_rising_dir_[1] << " " << first_tangential_vector_compo_from_rising_dir_[2] << " ";
        fic << azymuthal_vector_compo_[0] << " " << azymuthal_vector_compo_[1] << " " << azymuthal_vector_compo_[2] << " ";
        fic << r_sph_ << " " << theta_sph_ << " " << phi_sph_ << " ";
        fic << temperature_interp_[0] << " " << temperature_solution_[0] << " ";
        fic << normal_temperature_gradient_[0] << " " << normal_temperature_gradient_solution_[0] << " ";
        fic << normal_temperature_double_derivative_solution_[0] << " ";
        fic << tangential_temperature_gradient_first_[0] << " ";
        fic << tangential_temperature_gradient_second_[0] << " ";
        fic << tangential_temperature_gradient_first_from_rising_dir_[0] << " ";
        fic << azymuthal_temperature_gradient_[0] << " ";
        fic << temperature_diffusion_hessian_cartesian_trace_[0] << " ";
        fic << temperature_diffusion_hessian_trace_[0] << " ";
        fic << radial_temperature_diffusion_[0] << " ";
        fic << tangential_temperature_diffusion_[0] << " ";
        fic << surface_ << " " << thermal_flux_[0] << " " << *lambda_ << " " << *alpha_ << " " << Pr_l_ << " ";
        fic << velocity_shear_force_ << " " << velocity_shear_stress_ << " ";
        fic << pressure_interp_[0] << " ";
        fic << x_velocity_[0] << " " << y_velocity_[0] << " " << z_velocity_[0] << " ";
        fic << radial_velocity_[0] << " " << radial_velocity_corrected_[0] << " ";
        fic << radial_velocity_static_frame_[0] << " " << radial_velocity_advected_frame_[0] << " ";
        fic << first_tangential_velocity_[0] << " " << first_tangential_velocity_corrected_[0] << " ";
        fic << first_tangential_velocity_static_frame_[0] << " " << first_tangential_velocity_advected_frame_[0] << " ";
        fic << second_tangential_velocity_[0] << " " << second_tangential_velocity_corrected_[0] << " ";
        fic << second_tangential_velocity_static_frame_[0] << " " << second_tangential_velocity_advected_frame_[0] << " ";
        fic << first_tangential_velocity_from_rising_dir_[0] << " " << first_tangential_velocity_from_rising_dir_corrected_[0] << " ";
        fic << first_tangential_velocity_from_rising_dir_static_frame_[0] << " " << first_tangential_velocity_from_rising_dir_advected_frame_[0] << " ";
        fic << azymuthal_velocity_[0] << " " << azymuthal_velocity_corrected_[0] << " ";
        fic << azymuthal_velocity_static_frame_[0] << " " << azymuthal_velocity_advected_frame_[0] << " ";
        fic << normal_velocity_normal_gradient_[0] << " ";
        fic << first_tangential_velocity_normal_gradient_[0] << " ";
        fic << second_tangential_velocity_normal_gradient_[0] << " ";
        fic << first_tangential_velocity_gradient_from_rising_dir_[0] << " ";
        fic << azymuthal_velocity_normal_gradient_[0] << " ";
        fic << finl;
      }
  }
}

void IJK_One_Dimensional_Subproblem::post_process_radial_quantities(const int rank, const Nom& local_quantities_thermal_probes_time_index_folder)
{
  //	Nom probe_header = Nom("tstep\ttime\tthermalproblem\tsubproblem\ttemperature_gradient\ttemperature_double_deriv"
  //												 "ttemperature_gradient_tangential\ttemperature_gradient_tangential2\ttemperature_gradient_azymuthal"
  //												 "\tsurface\tthermal_flux\tlambda\talpha\tprandtl_liq"
  //												 "\tu_x\tu_y\tu_z\tu_r\tu_r_corr\tu_theta\tu_theta2\tu_phi\tdu_r_dr\tdu_theta_dr\tdu_theta2_dr\tdu_phi_dr");
  // if (Process::je_suis_maitre())
  {
    if (is_updated_ && is_post_processed_local_)
      {
        // const int reset = (!ref_ijk_ft_->get_reprise()) && (ref_ijk_ft_->get_tstep()==0);
        const int reset = 1;
        const int max_digit = 8;
        const int last_time_index = ref_ijk_ft_->get_tstep();
        const int nb_digit_index_post_pro = index_post_processing_ < 1 ? 1 : (int) (log10(index_post_processing_) + 1);
        const int nb_digit_index_global = global_subproblem_index_ < 1 ? 1 : (int) (log10(global_subproblem_index_) + 1);
        const int nb_digit_tstep = last_time_index < 1 ? 1 : (int) (log10(last_time_index) + 1);
        const int max_digit_rank = 3;
        const int nb_digit_rank = rank < 1 ? 1 : (int) (log10(rank) + 1);
        Nom probe_name = Nom("_thermal_rank_") + Nom(std::string(max_digit_rank - nb_digit_rank, '0')) + Nom(rank) + Nom("_thermal_subproblem_") + Nom(std::string(max_digit - nb_digit_index_post_pro, '0'))
                         + Nom(index_post_processing_)
                         + Nom("_global_") + Nom(std::string(max_digit - nb_digit_index_global, '0')) + Nom(global_subproblem_index_)
                         + Nom("_radial_quantities_time_index_") +
                         + Nom(std::string(max_digit - nb_digit_tstep, '0')) + Nom(last_time_index) + Nom(".out");
        Nom probe_header = Nom("tstep\tthermalrank\tpostproindex\tglobalsubproblem\tlocalsubproblem\ttime"
                               "\tnx\tny\tnz\tt1x\tt1y\tt2z\tt2x\tt2y\tt2z\ts1x\ts1y\ts1z\ts2x\ts2y\ts2z"
                               "\tradial_coord\ttemperature_interp\ttemperature_solution"
                               "\ttemperature_gradient\ttemperature_gradient_sol"
                               "\ttemperature_double_deriv_sol"
                               "\ttemperature_gradient_tangential\ttemperature_gradient_tangential2"
                               "\ttemperature_gradient_tangential_rise\ttemperature_gradient_azymuthal"
                               "\ttemperature_diffusion_hessian_cartesian_trace"
                               "\ttemperature_diffusion_hessian_trace"
                               "\tradial_temperature_diffusion"
                               "\ttangential_temperature_diffusion"
                               "\tsurface\tthermal_flux\tlambda\talpha\tprandtl_liq"
                               "\tshear\tforce"
                               "\tpressure"
                               "\tu_x\tu_y\tu_z"
                               "\tu_r\tu_r_corr\tu_r_static\tu_r_advected"
                               "\tu_theta\tu_theta_corr\tu_theta_static\tu_theta_advected"
                               "\tu_theta2\tu_theta2_corr\tu_theta2_static\tu_theta2_advected"
                               "\tu_theta_rise\tu_theta_rise_corr\tu_theta_rise_static\tu_theta_rise_advected"
                               "\tu_phi\tu_phi_corr\tu_phi_static\tu_phi_advected"
                               "\tdu_r_dr\tdu_theta_dr\tdu_theta2_dr\tdu_theta_rise_dr\tdu_phi_dr");
        SFichier fic = Open_file_folder(local_quantities_thermal_probes_time_index_folder, probe_name, probe_header, reset);
        const double last_time = ref_ijk_ft_->get_current_time() - ref_ijk_ft_->get_timestep();
        for (int i=0; i<(*points_per_thermal_subproblem_); i++)
          {
            fic << last_time_index << " ";
            fic << rank << " " << index_post_processing_ << " " << global_subproblem_index_ << " " << sub_problem_index_ << " ";
            fic << last_time << " ";
            fic << normal_vector_compo_[0] << " " << normal_vector_compo_[1] << " " << normal_vector_compo_[2] << " ";
            fic << first_tangential_vector_compo_[0] << " " << first_tangential_vector_compo_[1] << " " << first_tangential_vector_compo_[2] << " ";
            fic << second_tangential_vector_compo_[0] << " " << second_tangential_vector_compo_[1] << " " << second_tangential_vector_compo_[2] << " ";
            fic << first_tangential_vector_compo_from_rising_dir_[0] << " " << first_tangential_vector_compo_from_rising_dir_[1] << " " << first_tangential_vector_compo_from_rising_dir_[2] << " ";
            fic << azymuthal_vector_compo_[0] << " " << azymuthal_vector_compo_[1] << " " << azymuthal_vector_compo_[2] << " ";
            fic << r_sph_ << " " << theta_sph_ << " " << phi_sph_ << " ";
            fic << (*radial_coordinates_)[i] << " " << temperature_interp_[i] << " " << temperature_solution_[i] << " ";
            fic << normal_temperature_gradient_[i] << " " << normal_temperature_gradient_solution_[i] << " ";
            fic << normal_temperature_double_derivative_solution_[i] << " ";
            fic << tangential_temperature_gradient_first_[i] << " ";
            fic << tangential_temperature_gradient_second_[i] << " ";
            fic << tangential_temperature_gradient_first_from_rising_dir_[i] << " ";
            fic << azymuthal_temperature_gradient_[i] << " ";
            fic << temperature_diffusion_hessian_cartesian_trace_[i] << " ";
            fic << temperature_diffusion_hessian_trace_[i] << " ";
            fic << radial_temperature_diffusion_[i] << " ";
            fic << tangential_temperature_diffusion_[i] << " ";
            fic << surface_ << " " << thermal_flux_[i] << " " << *lambda_ << " " << *alpha_ << " " << Pr_l_ << " ";
            fic << shear_stress_[i] << " " << (shear_stress_[i] * surface_) << " ";
            fic << pressure_interp_[i] << " ";
            fic << x_velocity_[i] << " " << y_velocity_[i] << " " << z_velocity_[i] << " ";
            fic << radial_velocity_[i] << " " << radial_velocity_corrected_[i] << " ";
            fic << radial_velocity_static_frame_[i] << " " << radial_velocity_advected_frame_[i] << " ";
            fic << first_tangential_velocity_[i] << " " << first_tangential_velocity_corrected_[i] << " ";
            fic << first_tangential_velocity_static_frame_[i] << " " << first_tangential_velocity_advected_frame_[i] << " ";
            fic << second_tangential_velocity_[i] << " " << second_tangential_velocity_corrected_[i] << " ";
            fic << second_tangential_velocity_static_frame_[i] << " " << second_tangential_velocity_advected_frame_[i] << " ";
            fic << first_tangential_velocity_from_rising_dir_[i] << " " << first_tangential_velocity_from_rising_dir_corrected_[i] << " ";
            fic << first_tangential_velocity_from_rising_dir_static_frame_[i] << " " << first_tangential_velocity_from_rising_dir_advected_frame_[i] << " ";
            fic << azymuthal_velocity_[i] << " " << azymuthal_velocity_corrected_[i] << " ";
            fic << azymuthal_velocity_static_frame_[i] << " " << azymuthal_velocity_advected_frame_[i] << " ";
            fic << normal_velocity_normal_gradient_[i] << " ";
            fic << first_tangential_velocity_normal_gradient_[i] << " ";
            fic << second_tangential_velocity_normal_gradient_[i] << " ";
            fic << first_tangential_velocity_gradient_from_rising_dir_[i] << " ";
            fic << azymuthal_velocity_normal_gradient_[i] << " ";
            fic << finl;
          }
        fic.close();
      }
  }
}

double IJK_One_Dimensional_Subproblem::get_min_temperature() const
{
  double min_temperature_value=1e20;
  for (int i=0; i<temperature_solution_.size(); i++)
    min_temperature_value = std::min(min_temperature_value, temperature_solution_[i]);
  return min_temperature_value;
}

double IJK_One_Dimensional_Subproblem::get_max_temperature() const
{
  double max_temperature_value=-1e20;
  for (int i=0; i<temperature_solution_.size(); i++)
    max_temperature_value = std::max(max_temperature_value, temperature_solution_[i]);
  return max_temperature_value;
}

double IJK_One_Dimensional_Subproblem::get_min_temperature_domain_ends() const
{
  double min_temperature_value = std::min(temperature_solution_[0], temperature_solution_[temperature_solution_.size()-1]);
  return min_temperature_value;
}

double IJK_One_Dimensional_Subproblem::get_max_temperature_domain_ends() const
{
  double max_temperature_value = std::max(temperature_solution_[0], temperature_solution_[temperature_solution_.size()-1]);
  return max_temperature_value;
}
