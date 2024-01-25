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
// File      : IJK_Thermal_Subresolution.h
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/Temperature
//
/////////////////////////////////////////////////////////////////////////////

#ifndef IJK_Thermal_Subresolution_included
#define IJK_Thermal_Subresolution_included

#include <IJK_Thermal_base.h>
#include <IJK_Field.h>
#include <Boundary_Conditions_Thermique.h>
#include <IJK_Splitting.h>
#include <Parser.h>
#include <IJK_Lata_writer.h>
#include <OpConvQuickIJKScalar.h>
#include <OpConvCentre2IJKScalar.h>
#include <Ouvrir_fichier.h>
#include <Corrige_flux_FT.h>
#include <TRUST_Ref.h>
#include <Operateur_IJK_elem_diff_base.h>
#include <OpConvAmontIJK.h>
#include <OpConvDiscQuickIJKScalar.h>
#include <OpConvCentre4IJK.h>
#include <IJK_One_Dimensional_Subproblems.h>
#include <IJK_Finite_Difference_One_Dimensional_Matrix_Assembler.h>
#include <IJK_SolveSys_FD_thermal.h>
#include <MD_Vector_tools.h>

/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class IJK_Thermal_Subresolution
//
// <Description of class IJK_Thermal_Subresolution>
//
/////////////////////////////////////////////////////////////////////////////


class IJK_Thermal_Subresolution : public IJK_Thermal_base
{

  Declare_instanciable( IJK_Thermal_Subresolution ) ;
  friend class IJK_One_Dimensional_Subproblems;
  friend class IJK_One_Dimensional_Subproblem;

public :

  int initialize(const IJK_Splitting& splitting, const int idx) override;
  // void sauvegarder_temperature(Nom& lata_name, int idx) override;
  void update_thermal_properties() override;
  void post_process_after_temperature_increment() override;
  void set_param(Param& param) override;
  void compute_ghost_cell_numbers_for_subproblems(const IJK_Splitting& splitting, int ghost_init) override;

  double get_probes_length();
  // Entree& read_fd_solver(Entree& is);
  // void read_fd_solver(const Motcle& mot, Entree& is);
  int lire_motcle_non_standard(const Motcle&, Entree&) override;
  void set_thermal_subresolution_outputs(const Nom& interfacial_quantities_thermal_probes,
                                         const Nom& overall_bubbles_quantities,
                                         const Nom& local_quantities_thermal_probes_time_index_folder) override;

  const IJK_Field_double& get_debug_lrs_cells() const override
  {
    return debug_LRS_cells_;
  }
  const IJK_Field_double& get_temperature_cell_neighbours_debug() const override
  {
    if (find_temperature_cell_neighbours_ && debug_)
      return temperature_cell_neighbours_debug_;
    else
      return dummy_double_field_;
  }
  const IJK_Field_double& get_temperature_cell_neighbours() const override
  {
    if (find_temperature_cell_neighbours_)
      return temperature_cell_neighbours_;
    else
      return dummy_double_field_;
  }
  const IJK_Field_int& get_cell_neighbours_corrected() const override
  {
    if (find_temperature_cell_neighbours_)
      return neighbours_temperature_to_correct_;
    else
      return dummy_int_field_;
  }
  const IJK_Field_double& get_neighbours_temperature_colinearity_weighting() const override
  {
    if (find_temperature_cell_neighbours_)
      return neighbours_temperature_colinearity_weighting_;
    else
      return dummy_double_field_;
  }
  const FixedVector<IJK_Field_double,3>& get_cell_faces_corrected_diffusive() const override
  {
    if((convective_flux_correction_ || diffusive_flux_correction_) && store_cell_faces_corrected_)
      return cell_faces_corrected_diffusive_;
    else
      return dummy_double_vect_;
  }
  const FixedVector<IJK_Field_double,3>& get_cell_faces_corrected_convective() const override
  {
    if((convective_flux_correction_ || diffusive_flux_correction_) && store_cell_faces_corrected_)
      return cell_faces_corrected_convective_;
    else
      return dummy_double_vect_;
  }
  const FixedVector<IJK_Field_int,3>& get_cell_faces_corrected_bool() const override
  {
    if ((convective_flux_correction_ || diffusive_flux_correction_) && store_cell_faces_corrected_)
      return cell_faces_corrected_bool_;
    else
      return dummy_int_vect_;
  }
  const FixedVector<IJK_Field_int,3>& get_cell_faces_neighbours_corrected_diag_bool() const override
  {
    if (find_cell_neighbours_for_fluxes_spherical_correction_)
      return cell_faces_neighbours_corrected_diag_bool_;
    else
      return dummy_int_vect_;
  }
  const FixedVector<IJK_Field_int,3>& get_cell_faces_neighbours_corrected_all_bool() const override
  {
    if (find_reachable_fluxes_)
      return cell_faces_neighbours_corrected_all_bool_;
    else
      return dummy_int_vect_;
  }
  const FixedVector<IJK_Field_int,3>& get_cell_faces_neighbours_corrected_min_max_bool() const override
  {
    if (find_reachable_fluxes_)
      return cell_faces_neighbours_corrected_min_max_bool_;
    else
      return dummy_int_vect_;
  }
  const FixedVector<IJK_Field_double,3>& get_cell_faces_neighbours_corrected_convective() const override
  {
    if (use_reachable_fluxes_ && convective_flux_correction_)
      return cell_faces_neighbours_corrected_convective_;
    else
      return dummy_double_vect_;
  }
  const FixedVector<IJK_Field_double,3>& get_cell_faces_neighbours_corrected_diffusive() const override
  {
    if (use_reachable_fluxes_ && diffusive_flux_correction_)
      return cell_faces_neighbours_corrected_diffusive_;
    else
      return dummy_double_vect_;
  }
  const FixedVector<IJK_Field_double,3>& get_neighbours_faces_weighting_colinearity() const override
  {
    if (use_reachable_fluxes_ && neighbours_colinearity_weighting_)
      return neighbours_faces_weighting_colinearity_;
    else
      return dummy_double_vect_;
  }
  const IJK_Field_int& get_cell_neighbours_corrected_trimmed() const override
  {
    if (find_reachable_fluxes_)
      return neighbours_temperature_to_correct_trimmed_;
    else
      return dummy_int_field_;
  }
  int get_disable_post_processing_probes_out_files() const override
  {
    return disable_post_processing_probes_out_files_;
  }
  int get_first_step_thermals_post() override { return first_step_thermals_post_; }

protected :
  void reset_subresolution_distributed_vectors();
  void compute_thermal_subproblems() override;
  void compute_diffusion_increment() override;
  void correct_temperature_increment_for_interface_leaving_cell() override;
  void correct_any_temperature_fields_for_eulerian_fluxes(IJK_Field_double& temperature);
  void correct_temperature_for_eulerian_fluxes() override;
  void store_temperature_before_extrapolation() override;
  void compare_temperature_fields(const IJK_Field_double& temperature,
                                  const IJK_Field_double& temperature_ana,
                                  IJK_Field_double& error_temperature_ana,
                                  IJK_Field_double& error_temperature_ana_rel);
  void evaluate_total_liquid_absolute_parameter(const IJK_Field_double& field,
                                                double& total_parameter);
  void evaluate_total_liquid_parameter_squared(const IJK_Field_double& field,
                                               double& total_parameter);
  void correct_any_temperature_field_for_visu(IJK_Field_double& temperature);
  void correct_temperature_for_visu() override;
  void clip_temperature_values() override;
  void compute_mean_liquid_temperature();
  void compute_overall_probes_parameters();

  void pre_initialise_thermal_subproblems_any_matrices();
  void pre_initialise_thermal_subproblems_matrices();

  void interpolate_indicator_on_probes();
  void clear_problems_colliding_bubbles();
  void interpolate_project_velocities_on_probes();
  void reajust_probes_length();
  void compute_radial_subresolution_convection_diffusion_operators();
  void compute_local_substep();
  void prepare_temporal_schemes();
  void compute_source_terms_impose_subresolution_boundary_conditions();
  void compute_add_subresolution_source_terms();
  void compute_subresolution_temporal_explicit_implicit_matrices();
  void approximate_temperature_increment_material_derivative();
  void compute_radial_first_second_order_operators(Matrice& radial_first_order_operator_raw,
                                                   Matrice& radial_second_order_operator_raw,
                                                   Matrice& radial_first_order_operator,
                                                   Matrice& radial_second_order_operator);
  void compute_first_order_operator_raw(Matrice& radial_first_order_operator);
  void compute_first_order_operator(Matrice& radial_first_order_operator, double dr);
  void compute_second_order_operator(Matrice& radial_second_order_operator, double dr);
  void compute_second_order_operator_raw(Matrice& radial_second_order_operator);
  void initialise_identity_matrices(Matrice& identity_matrix,
                                    Matrice& identity_matrix_subproblems);
  void initialise_identity_matrices_sparse(Matrice& identity_matrix,
                                           Matrice& identity_matrix_subproblems);
  void initialise_radial_convection_operator(Matrice& radial_first_order_operator,
                                             Matrice& radial_convection_matrix);
  void initialise_radial_convection_operator_sparse(Matrice& radial_first_order_operator,
                                                    Matrice& radial_convection_matrix);
  void initialise_radial_diffusion_operator(Matrice& radial_second_order_operator,
                                            Matrice& radial_diffusion_matrix);
  void initialise_radial_diffusion_operator_sparse(Matrice& radial_second_order_operator,
                                                   Matrice& radial_diffusion_matrix);
  // int copy_local_unknwowns_rhs();
  void convert_into_sparse_matrix();
  void compute_md_vector();
  void retrieve_temperature_solution();
  void check_wrong_values_rhs();
  void initialise_thermal_subproblems_list();
  void initialise_thermal_subproblems();
  void solve_thermal_subproblems();
  void prepare_thermal_flux_correction();
  void compute_min_max_reachable_fluxes();
  void update_intersections() override;
  void compute_convective_diffusive_fluxes_face_centre() override;
  void compute_convective_fluxes_face_centre() override;
  void compute_diffusive_fluxes_face_centre() override;
  void compute_temperature_cell_centres(const int first_corr) override;
  void compute_temperature_cell_centres_first_correction();
  void compute_temperature_cell_centres_second_correction();
  void replace_temperature_cell_centres_neighbours(const int& use_neighbours_temperature_to_correct_trimmed);
  void prepare_ij_fluxes_k_layers() override;
  void set_zero_temperature_increment() override;
  void clean_thermal_subproblems() override;
  void clean_ijk_intersections() override;
  void clean_add_thermal_subproblems();
  void enforce_periodic_temperature_boundary_value() override;
  void correct_operators_for_visu() override;

  double get_modified_time() override;
  void compute_temperature_init() override;
  void set_field_temperature_per_bubble(const int index_bubble);
  Nom compute_quasi_static_spherical_diffusion_expression(const double& time_scope, const int index_bubble, const int index_bubble_real);
  Nom generate_expression_temperature_ini(const double& time_scope, const double x, const double y, const double z);
  void approx_erf_inverse(const double& x, double& res);
  void set_field_T_ana() override;
  void calculer_ecart_T_ana() override { ; };
  double compute_spherical_steady_dirichlet_left_right_value(const double& r);
  double compute_spherical_steady_dirichlet_left_right_derivative_value(const double& r, const double& temperature_prev);
  double compute_spherical_steady_dirichlet_left_right_integral(const double& temperature_end_prev);
  double find_time_dichotomy_integral(const double& temperature_integral, double& temperature_end_prev);
  void compute_Nusselt_spherical_diffusion();
  double get_time_inflection_derivative(double& temperature_end_min);
  double find_time_dichotomy_derivative(const double& temperature_derivative, double& temperature_limit_left, double& temperature_limit_right);

  /* compute_rho_cp_u_mean() May be clearly overridden later */
  double compute_rho_cp_u_mean(const IJK_Field_double& vx) override { return IJK_Thermal_base::compute_rho_cp_u_mean(vx); };

  int reference_gfm_on_probes_;
  int compute_normal_derivatives_on_reference_probes_;

  int disable_spherical_diffusion_start_;
  int single_centred_bubble_;
  int computed_centred_bubble_start_;
  double single_centred_bubble_radius_ini_;
  double probes_end_value_start_;
  double probes_end_value_coeff_;
  int temperature_ini_type_;
  double modified_time_init_;
  int spherical_diffusion_;
  double nusselt_spherical_diffusion_;
  double nusselt_spherical_diffusion_liquid_;
  double heat_flux_spherical_;
  enum temperature_ini_dict { local_criteria, integral_criteria, derivative_criteria };
  double mean_liquid_temperature_;

  int disable_mixed_cells_increment_;
  int enable_mixed_cells_increment_;
  int allow_temperature_correction_for_visu_;
  int disable_subresolution_;
  int diffusive_flux_correction_;
  int convective_flux_correction_;
  int impose_fo_flux_correction_;
  int disable_fo_flux_correction_;
  int subproblem_temperature_extension_; // ghost fluid extension based on the interfacial gradient computed with the subproblem

  int override_vapour_mixed_values_; // For debug purposes

  IJK_One_Dimensional_Subproblems thermal_local_subproblems_;
  int points_per_thermal_subproblem_;
  double coeff_distance_diagonal_ = 2.;
  double probe_length_;
  double dr_;
  DoubleVect radial_coordinates_;
  Matrice identity_matrix_explicit_implicit_;
  Matrice radial_first_order_operator_raw_;
  Matrice radial_second_order_operator_raw_;
  Matrice radial_first_order_operator_;
  Matrice radial_second_order_operator_;
  Matrice identity_matrix_subproblems_;
  Matrice radial_diffusion_matrix_;
  Matrice radial_convection_matrix_;

  Matrice radial_diffusion_matrix_test_;

  FixedVector<ArrOfInt, 6> first_indices_sparse_matrix_;

  int initialise_sparse_indices_;
  /*
   * Thermal subproblems are regrouped in a single linear system AX=b
   * on each processor !
   */
  IJK_Finite_Difference_One_Dimensional_Matrix_Assembler finite_difference_assembler_;
  Matrice thermal_subproblems_matrix_assembly_;
  DoubleVect thermal_subproblems_temperature_solution_ini_;
  DoubleVect thermal_subproblems_rhs_assembly_;
  DoubleVect thermal_subproblems_temperature_solution_;

  /*
   * TODO: Cast the matrice with Matrice Morse directly (not Matrice Bloc)
   */
  Matrice * thermal_subproblems_matrix_assembly_for_solver_ref_;
  Matrice thermal_subproblems_matrix_assembly_for_solver_;
  Matrice thermal_subproblems_matrix_assembly_for_solver_reduced_;

  // SolveurSys one_dimensional_advection_diffusion_thermal_solver_;
  IJK_SolveSys_FD_thermal one_dimensional_advection_diffusion_thermal_solver_;
  IJK_SolveSys_FD_thermal one_dimensional_advection_diffusion_thermal_solver_implicit_;
  MD_Vector md_;
  Motcles fd_solvers_;
  Motcles fd_solvers_jdd_;
  int fd_solver_rank_;
  Nom fd_solver_type_;
  int discrete_integral_;
  int quadtree_levels_;

  int compute_tangential_variables_;

  int boundary_condition_interface_;
  Motcles boundary_condition_interface_dict_;
  double interfacial_boundary_condition_value_;
  int impose_boundary_condition_interface_from_simulation_;
  int boundary_condition_end_;
  Motcles boundary_condition_end_dict_;
  double end_boundary_condition_value_;
  int impose_user_boundary_condition_end_value_;
  int source_terms_type_;
  Motcles source_terms_type_dict_;
  int source_terms_correction_;
  int advected_frame_of_reference_;
  int neglect_frame_of_reference_radial_advection_;
  int approximate_temperature_increment_;

  /*
   * Some tries to do explicit temporal variations at the beginning
   */
  bool is_first_time_step_;
  int first_time_step_temporal_;
  int first_time_step_explicit_;
  int first_time_step_implicit_;
  int local_diffusion_fourier_priority_ = 0;
  int nb_iter_explicit_local_ = 0;
  double local_fourier_;
  double local_cfl_;
  double delta_T_subcooled_overheated_=-1.;
  double convection_diffusion_time_scale_factor_ = 1.;
  double local_fourier_time_step_probe_length_ = 0.;
  double local_cfl_time_step_probe_length_ = 0.;
  double local_dt_cfl_ = 0.;
  double local_dt_cfl_min_delta_xyz_ = 0.;
  double local_dt_cfl_counter_ = 0.;
  double local_dt_fourier_counter_ = 0.;

  double error_temperature_ana_total_;
  double error_temperature_ana_squared_total_;
  double error_temperature_ana_rel_total_;

  /*
   * Some tries to make the probe length varies at the beginning of the simulation
   */
  int first_time_step_varying_probes_ = 0;
  int probe_variations_enabled_ = 0;
  int probe_variations_priority_ = 0;
  int disable_interpolation_in_mixed_cells_ = 0;
  int keep_temperature_extrapolated_from_LRS_ = 0;
  int max_u_radial_=0;

  IJK_Field_double debug_LRS_cells_;
  int distance_cell_faces_from_lrs_;
  int  disable_distance_cell_faces_from_lrs_;

  int pre_initialise_thermal_subproblems_list_;
  double pre_factor_subproblems_number_;
  int remove_append_subproblems_;
  int use_sparse_matrix_;
  int global_probes_characteristics_ = 1;

  int correct_temperature_cell_neighbours_first_iter_;
  int correct_first_iter_deactivate_cell_neighbours_;
  int find_temperature_cell_neighbours_;
  int use_temperature_cell_neighbours_;
  int correct_neighbours_using_probe_length_;
  int neighbours_corrected_rank_;
  int neighbours_weighting_;
  int neighbours_colinearity_weighting_;
  int neighbours_distance_weighting_;
  int neighbours_colinearity_distance_weighting_;
  IJK_Field_double temperature_cell_neighbours_;
  IJK_Field_double temperature_cell_neighbours_debug_;
  IJK_Field_int neighbours_temperature_to_correct_;
  IJK_Field_double neighbours_temperature_colinearity_weighting_;

  int clip_temperature_values_;
  int enforce_periodic_boundary_value_;
  int stencil_periodic_boundary_value_;

  int disable_post_processing_probes_out_files_;

  /*
   * Pure cells corrected for visualisation
   */
  FixedVector<IJK_Field_double,3> cell_faces_corrected_diffusive_;
  FixedVector<IJK_Field_double,3> cell_faces_corrected_convective_;
  FixedVector<IJK_Field_int,3> cell_faces_corrected_bool_;
  int store_cell_faces_corrected_;

  /*
   * Neighbouring faces in the diagonal
   */
  FixedVector<IJK_Field_int,3> cell_faces_neighbours_corrected_diag_bool_;
  int find_cell_neighbours_for_fluxes_spherical_correction_;
  int use_cell_neighbours_for_fluxes_spherical_correction_;

  /*
   * All reachable faces to correct
   */
  int find_reachable_fluxes_;
  int use_reachable_fluxes_;
  int keep_first_reachable_fluxes_;
  FixedVector<IJK_Field_int,3> cell_faces_neighbours_corrected_all_bool_;
  FixedVector<IJK_Field_double,3> cell_faces_neighbours_corrected_convective_;
  FixedVector<IJK_Field_double,3> cell_faces_neighbours_corrected_diffusive_;
  FixedVector<IJK_Field_double, 3> neighbours_faces_weighting_colinearity_;
  FixedVector<IJK_Field_int,3> cell_faces_neighbours_corrected_min_max_bool_;
  IJK_Field_int neighbours_temperature_to_correct_trimmed_;
  int neighbours_last_faces_weighting_;
  int neighbours_last_faces_colinearity_weighting_;
  int neighbours_last_faces_colinearity_face_weighting_;
  int neighbours_last_faces_distance_weighting_;
  int neighbours_last_faces_distance_colinearity_weighting_;
  int neighbours_last_faces_distance_colinearity_face_weighting_;

  int post_process_all_probes_;
  int nb_theta_post_pro_;
  int nb_phi_post_pro_;
  int nb_probes_post_pro_;

  int interp_eulerian_;
  int first_step_thermals_post_;
  int disable_first_step_thermals_post_;

  int copy_fluxes_on_every_procs_;
  int copy_temperature_on_every_procs_;

};

#endif /* IJK_Thermal_Subresolution_included */
