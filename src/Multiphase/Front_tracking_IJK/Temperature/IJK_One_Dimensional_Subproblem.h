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
// File      : IJK_One_Dimensional_Subproblem.h
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/Temperature
//
/////////////////////////////////////////////////////////////////////////////

#ifndef IJK_One_Dimensional_Subproblem_included
#define IJK_One_Dimensional_Subproblem_included

#include <Objet_U.h>
#include <IJK_Field.h>
#include <IJK_Interfaces.h>
#include <Linear_algebra_tools.h>
#include <FixedVector.h>
#include <TRUSTArrays.h>
#include <TRUSTTab.h>
#include <Vecteur3.h>
#include <Matrice33.h>
#include <Matrice.h>
#include <IJK_Finite_Difference_One_Dimensional_Matrix_Assembler.h>

// #include <TRUSTArray.h>
#define INVALID_TEMPERATURE 1e10
#define INVALID_FIELD 1e10
#define INVALID_VELOCITY 1e-12
#define INVALID_INTERP 1.e20
#define INVALID_INTERP_TEST 1.e19
#define INVALID_VELOCITY_CFL 1e-20
#define INVALID_SOURCE_TERM 1e-20
#define NEIGHBOURS_FIRST_DIR {-1., -1., 1., 1.}
#define NEIGHBOURS_SECOND_DIR {-1., 1., -1., 1.}
#define NEIGHBOURS_I {-1, 1, 0, 0, 0, 0}
#define NEIGHBOURS_J {0, 0, -1, 1, 0, 0}
#define NEIGHBOURS_K {0, 0, 0, 0, -1, 1}
#define NEIGHBOURS_FACES_I {0, 1, 0, 0, 0, 0}
#define NEIGHBOURS_FACES_J {0, 0, 0, 1, 0, 0}
#define NEIGHBOURS_FACES_K {0, 0, 0, 0, 0, 1}
#define LIQUID_INDICATOR_TEST 1.-1.e-12
#define VAPOUR_INDICATOR_TEST 1.e-12
#define FACES_DIR {0, 0, 1, 1, 2, 2}
#define FLUX_SIGN {-1, -1, -1, -1, -1, -1}

/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class IJK_One_Dimensional_Subproblem
//
// <Description of class IJK_One_Dimensional_Subproblem>
//
/////////////////////////////////////////////////////////////////////////////
class IJK_FT_double;

class IJK_One_Dimensional_Subproblem : public Objet_U
{

  Declare_instanciable( IJK_One_Dimensional_Subproblem ) ;

public :
  IJK_One_Dimensional_Subproblem(const IJK_FT_double& ijk_ft);
  void associer(const IJK_FT_double& ijk_ft) { ref_ijk_ft_ = ijk_ft; };
  void associate_sub_problem_to_inputs(int debug,
                                       int sub_problem_index,
                                       int i, int j, int k,
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
                                       ArrOfDouble bubbles_barycentre,
                                       int advected_frame_of_reference,
                                       int neglect_frame_of_reference_radial_advection,
                                       const int& points_per_thermal_subproblem,
                                       const double& alpha,
                                       const double& lambda,
                                       const double& coeff_distance_diagonal,
                                       const double& cell_diagonal,
                                       const double& dr_base,
                                       const DoubleVect& radial_coordinates,
                                       const Matrice& identity_matrix_explicit_implicit_raw,
                                       const Matrice& radial_first_order_operator_raw,
                                       const Matrice& radial_second_order_operator_raw,
                                       const Matrice& radial_first_order_operator,
                                       const Matrice& radial_second_order_operator,
                                       Matrice& identity_matrix_subproblems,
                                       Matrice& radial_diffusion_matrix,
                                       Matrice& radial_convection_matrix,
                                       const IJK_Interfaces& interfaces,
                                       const double& indicator,
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
                                       const FixedVector<IJK_Field_double, 3>& grad_T_elem,
                                       const FixedVector<IJK_Field_double, 3>& hess_diag_T_elem,
                                       const FixedVector<IJK_Field_double, 3>& hess_cross_T_elem,
                                       IJK_Finite_Difference_One_Dimensional_Matrix_Assembler& finite_difference_assembler,
                                       Matrice& thermal_subproblems_matrix_assembly,
                                       DoubleVect& thermal_subproblems_rhs_assembly,
                                       DoubleVect& thermal_subproblems_temperature_solution_ini,
                                       DoubleVect& thermal_subproblems_temperature_solution,
                                       const int& source_terms_type,
                                       const int& source_terms_correction,
                                       const bool& is_first_time_step,
                                       int& first_time_step_temporal,
                                       const int& first_time_step_explicit,
                                       const double& local_fourier,
                                       const double& local_cfl,
                                       const double& min_delta_xyz,
                                       const double& delta_T_subcooled_overheated,
                                       const int& first_time_step_varying_probes,
                                       const int& probe_variations_priority,
                                       const int& disable_interpolation_in_mixed_cells,
                                       const int& max_u_radial,
                                       const int& correct_fluxes,
                                       const int& distance_cell_faces_from_lrs,
                                       const int& pre_initialise_thermal_subproblems_list,
                                       const int& correct_temperature_cell_neighbours);
  void interpolate_project_velocities_on_probes();
  void reajust_probe_length();
  void compute_modified_probe_length_condition();
  void compute_distance_cell_centre();
  void compute_distance_faces_centres();
  void compute_distance_cell_centres_neighbours();
  double compute_min_distance_pure_face_centre();
  double compute_min_distance_pure_face_vertices();
  double compute_max_distance_pure_face_centre();
  double compute_max_distance_pure_face_vertices();
  void compute_vertex_position(const int& vertex_number,
                               const int& face_dir,
                               const Vecteur3& bary_face,
                               double& distance_vertex_centre,
                               Vecteur3& bary_vertex);
  void compute_modified_probe_length(const int& probe_variations_enabled);
  void compute_radial_convection_diffusion_operators();
  void prepare_temporal_schemes();
  void prepare_boundary_conditions(DoubleVect& thermal_subproblems_rhs_assembly,
                                   DoubleVect& thermal_subproblems_temperature_solution_ini,
                                   const int& boundary_condition_interface,
                                   const double& interfacial_boundary_condition_value,
                                   const int& impose_boundary_condition_interface_from_simulation,
                                   const int& boundary_condition_end,
                                   const double& end_boundary_condition_value,
                                   const int& impose_user_boundary_condition_end_value);
  void compute_source_terms_impose_boundary_conditions(DoubleVect& thermal_subproblems_rhs_assembly,
                                                       DoubleVect& thermal_subproblems_temperature_solution_ini,
                                                       const int& boundary_condition_interface,
                                                       const double& interfacial_boundary_condition_value,
                                                       const int& impose_boundary_condition_interface_from_simulation,
                                                       const int& boundary_condition_end,
                                                       const double& end_boundary_condition_value,
                                                       const int& impose_user_boundary_condition_end_value);
  void compute_add_source_terms();
  void compute_temporal_explicit_implicit_matrices();
  void approximate_temperature_increment_material_derivative();
  void retrieve_temperature_solution();
  void retrieve_radial_quantities();
  void compute_local_temperature_gradient_solution();
  double get_interfacial_gradient_corrected() const;
  double get_interfacial_double_derivative_corrected() const;

  void compute_local_velocity_gradient();
  double get_normal_velocity_normal_gradient() const;
  double get_tangential_velocity_normal_gradient() const;
  double get_second_tangential_velocity_normal_gradient() const;
  double get_azymuthal_velocity_normal_gradient() const;

  void get_ijk_indices(int& i, int& j, int& k) const;
  double get_field_profile_at_point(const double& dist, const DoubleVect& field, const int temp_bool) const;
  double get_temperature_profile_at_point(const double& dist) const;
  double get_temperature_times_velocity_profile_at_point(const double& dist, const int& dir) const;
  DoubleVect get_field_discrete_integral_velocity_weighting_at_point(const double& dist, const int& levels, const int& dir, const DoubleVect& field, const int vel) const;
  DoubleVect get_field_times_velocity_discrete_integral_at_point(const double& dist, const int& levels, const int& dir, const DoubleVect& field) const;
  DoubleVect get_field_discrete_integral_at_point(const double& dist, const int& levels, const int& dir, const DoubleVect& field) const;
  double get_velocity_weighting(const double& dist, const int& dir, const int vel) const;
  DoubleVect get_temperature_profile_discrete_integral_at_point(const double& dist, const int& levels, const int& dir) const;
  DoubleVect get_temperature_times_velocity_profile_discrete_integral_at_point(const double& dist, const int& levels, const int& dir) const;
  DoubleVect get_temperature_gradient_profile_discrete_integral_at_point(const double& dist, const int& levels, const int& dir) const;
  DoubleVect get_temperature_gradient_times_conductivity_profile_discrete_integral_at_point(const double& dist, const int& levels, const int& dir) const;
  void get_field_discrete_value_recursive(const int& ilevel, const int& max_level,
                                          const int& dir, const double& dist,
                                          const int& vel,
                                          const double& surface,
                                          const DoubleVect& field,
                                          const double dl1_parent,
                                          const double dl2_parent,
                                          Vecteur3& point_coords_parent,
                                          DoubleVect& discrete_values,
                                          int& value_counter) const;
  double get_velocity_component_at_point(const double& dist, const int& dir) const;
  double get_temperature_gradient_profile_at_point(const double& dist, const int& dir) const;
  double get_temperature_gradient_times_conductivity_profile_at_point(const double& dist, const int& dir) const;
  void get_discrete_two_dimensional_spacing(const int& dir, const int& level,
                                            const double& first_dir, const double& second_dir,
                                            double& dl1, double& dl2, Vecteur3& point_coords) const;
  double get_discrete_surface_at_level(const int& dir, const int& level) const;
  void thermal_subresolution_outputs(SFichier& fic, const int rank);

  double get_min_temperature() const;
  double get_max_temperature() const;
  double get_min_temperature_domain_ends() const;
  double get_max_temperature_domain_ends() const;
  const double& get_local_time_step_round() const
  {
    return local_time_step_round_;
  };
  const int& get_nb_iter_explicit() const
  {
    return nb_iter_explicit_;
  };
  void set_local_time_step(const double& local_time_step)
  {
    local_time_step_overall_ = local_time_step;
  };
  const double& get_local_fourier_time_step_probe_length() const
  {
    return local_fourier_time_step_probe_length_;
  };
  const double& get_local_cfl_time_step_probe_length() const
  {
    return local_cfl_time_step_probe_length_;
  };
  const double& get_local_dt_cfl() const
  {
    return local_dt_cfl_;
  };
  const double& get_local_dt_cfl_min_delta_xyz() const
  {
    return local_dt_cfl_min_delta_xyz_;
  };
  const int& get_probe_variations_enabled() const
  {
    return probe_variations_enabled_;
  };
  const double& get_dist_cell() const
  {
    return cell_centre_distance_;
  };
  const FixedVector<double,6>& get_dist_faces() const
  {
    return face_centres_distance_;
  };
  const Vecteur3& get_bary_facet() const
  {
    return facet_barycentre_;
  };
  const int& get_end_index_subproblem() const
  {
    return end_index_;
  }
  const int& get_dxyz_increment_bool() const
  {
    return dxyz_increment_bool_;
  }
  const FixedVector<int,3>& get_pure_neighbours_corrected_sign() const
  {
    return pure_neighbours_corrected_sign_;
  }
  const std::vector<std::vector<std::vector<bool>>>& get_pure_neighbours_to_correct() const
  {
    return pure_neighbours_to_correct_;
  }
  const std::vector<std::vector<std::vector<double>>>& get_pure_neighbours_corrected_distance() const
  {
    return pure_neighbours_corrected_distance_;
  }
protected :
  void associate_cell_ijk(int i, int j, int k) { index_i_ = i; index_j_=j; index_k_=k; };
  void associate_sub_problem_temporal_params(const bool& is_first_time_step,
                                             int& first_time_step_temporal,
                                             const int& first_time_step_explicit,
                                             const double& local_fourier,
                                             const double& local_cfl,
                                             const double& min_delta_xyz);
  void associate_varying_probes_params(const int& first_time_step_varying_probes,
                                       const int& probe_variations_priority,
                                       const int& disable_interpolation_in_mixed_cells);
  void associate_compos(int compo_connex) { compo_connex_ = compo_connex; };
  void associate_compos(int compo_connex, int compo_group) { compo_connex_ = compo_connex; compo_group_ = compo_group; };
  void associate_interface_related_parameters(double distance, double curvature, double interfacial_area, ArrOfDouble facet_barycentre, ArrOfDouble normal_vector)
  {
    distance_ = distance;
    curvature_ = curvature;
    interfacial_area_ = interfacial_area;
    facet_barycentre_ = Vecteur3(facet_barycentre);
    normal_vector_compo_ = Vecteur3(normal_vector);
  };
  void associate_rising_velocity(double bubble_rising_velocity, ArrOfDouble bubble_rising_vector, ArrOfDouble bubble_barycentre)
  {
    bubble_rising_velocity_ = bubble_rising_velocity;
    bubble_rising_vector_ = Vecteur3(bubble_rising_vector);
    bubble_barycentre_ = Vecteur3(bubble_barycentre);
    bubble_rising_velocity_compo_ =  bubble_rising_vector_;
    bubble_rising_velocity_compo_ *= bubble_rising_velocity_;
  };

  void associate_eulerian_fields_references(const IJK_Interfaces& interfaces,
                                            const double& indicator,
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
                                            const FixedVector<IJK_Field_double, 3>& grad_T_elem,
                                            const FixedVector<IJK_Field_double, 3>& hess_diag_T_elem,
                                            const FixedVector<IJK_Field_double, 3>& hess_cross_T_elem);
  void associate_probe_parameters(const int& points_per_thermal_subproblem,
                                  const double& alpha,
                                  const double& lambda,
                                  const double& coeff_distance_diagonal,
                                  const double& cell_diagonal,
                                  const double& dr_base,
                                  const DoubleVect& radial_coordinates);
  void associate_finite_difference_operators(const Matrice& radial_first_order_operator_raw,
                                             const Matrice& radial_second_order_operator_raw,
                                             const Matrice& radial_first_order_operator,
                                             const Matrice& radial_second_order_operator,
                                             const Matrice& identity_matrix_explicit_implicit,
                                             Matrice& identity_matrix_subproblems,
                                             Matrice& radial_diffusion_matrix,
                                             Matrice& radial_convection_matrix);
  void initialise_thermal_probe();
  void compute_interface_basis_vectors();
  void compute_pure_spherical_basis_vectors();
  void compute_local_discretisation();
  void compute_local_time_step();
  const int *  increase_number_of_points();
  void compute_identity_matrix_local(Matrice& identity_matrix_explicit_implicit);
  void compute_first_order_operator_local(Matrice& radial_first_order_operator);
  void compute_second_order_operator_local(Matrice& second_first_order_operator);
  void recompute_finite_difference_matrices();
  void compute_first_order_operator_local_varying_probe_length(const Matrice * radial_first_order_operator);
  void compute_second_order_operator_local_varying_probe_length(const Matrice * radial_second_order_operator);
  void recompute_finite_difference_matrices_varying_probe_length();
  void initialise_radial_convection_operator_local();
  void initialise_radial_diffusion_operator_local();
  void initialise_identity_operator_local();
  void interpolate_cartesian_velocities_on_probes();
  void compute_velocity_magnitude();
  void project_velocities_on_probes();
  void correct_velocities();
  void correct_velocity(const DoubleVect& velocity, DoubleVect& velocity_corrected);
  void correct_velocity_rise(const DoubleVect& velocity, const Vecteur3& basis, DoubleVect& velocity_corrected);
  void correct_radial_velocity();
  void project_cartesian_onto_basis_vector(const DoubleVect& compo_x, const DoubleVect& compo_y, const DoubleVect& compo_z, const Vecteur3& basis, DoubleVect& projection);
  void project_basis_vector_onto_cartesian_dir(const int& dir, const DoubleVect& compo_u, const DoubleVect& compo_v, const DoubleVect& compo_w,
                                               const Vecteur3& basis_u, const Vecteur3& basis_v, const Vecteur3& basis_w,
                                               DoubleVect& projection);
  void interpolate_temperature_on_probe();
  void interpolate_temperature_gradient_on_probe();
  void project_temperature_gradient_on_probes();
  void interpolate_temperature_hessian_on_probe();
  void project_temperature_hessian_on_probes();
  void retrieve_temperature_diffusion_spherical_on_probes();
  void compute_projection_matrix_cartesian_to_local_spherical();
  void project_matrix_on_basis(const Matrice33& projection_matrix, const Matrice33& inverse_projection_matrix, const Matrice33& matrix, Matrice33& projected_matrix);
  void approximate_partial_temperature_time_increment();
  void approximate_temperature_material_derivatives();
  void approximate_temperature_material_derivatives(const Vecteur3& normal_vector_compo,
                                                    const Vecteur3& first_tangential_vector_compo,
                                                    const Vecteur3& second_tangential_vector_compo,
                                                    const DoubleVect& radial_velocity_frame,
                                                    const DoubleVect& first_tangential_velocity_frame,
                                                    const DoubleVect& second_tangential_velocity_frame,
                                                    const DoubleVect& temperature_time_increment,
                                                    DoubleVect& convective_term,
                                                    DoubleVect& material_derivative);
  void correct_tangential_temperature_gradient(DoubleVect& tangential_convection_source_terms);
  void correct_tangential_temperature_hessian(DoubleVect& tangential_diffusion_source_terms);
  void find_interval(const double& dist, int& left_interval, int& right_interval) const;

  void post_process_interfacial_quantities(SFichier& fic, const int rank);
  void post_process_radial_quantities(const int rank);

  enum Boundary_conditions { default_bc=-1, dirichlet, neumann, flux_jump };

  int debug_;
  int advected_frame_of_reference_=1;
  int neglect_frame_of_reference_radial_advection_=1;
  /*
   * FIXME: Should I use only references or just for IJK_Field_double ?
   * Should I use IJK_Field_local_double or IJK_Field_double as pointers ?
   */
  int sub_problem_index_ = 0;
  int index_i_ = 0, index_j_ = 0, index_k_ = 0;
  int compo_connex_ = -1;
  int compo_group_ = -1;
  double distance_ = 0.;
  double curvature_ = 0.;
  double interfacial_area_ = 0.;
  double osculating_radius_ = 0.;
  Vecteur3 facet_barycentre_;
  Vecteur3 normal_vector_compo_;

  double bubble_rising_velocity_ = 0.;
  Vecteur3 bubble_rising_vector_;
  Vecteur3 bubble_rising_velocity_compo_;
  Vecteur3 bubble_barycentre_;
  Vecteur3 facet_barycentre_relative_;

  Vecteur3 osculating_sphere_centre_;

  /*
   * Several ways to calculate the tangential vector !
   * Either by considering a unique tangential vector (pure spherical)
   * Or two tangential vectors (osculating sphere)
   */
  Vecteur3 interfacial_velocity_compo_;
  Vecteur3 interfacial_tangential_velocity_compo_;

  Vecteur3 first_tangential_vector_compo_;
  Vecteur3 second_tangential_vector_compo_;

  Vecteur3 first_tangential_vector_compo_from_rising_dir_;
  Vecteur3 azymuthal_vector_compo_raw_;
  Vecteur3 azymuthal_vector_compo_;

  bool tangential_from_rising_vel_ = false;
  Vecteur3 * first_tangential_vector_compo_solver_;
  Vecteur3 * second_tangential_vector_compo_solver_;

  Vecteur3 interfacial_temperature_gradient_compo_;
  Matrice33 interfacial_temperature_hessian_compo_;
  Matrice33 projection_matrix_;
  Matrice33 inverse_projection_matrix_;
  Matrice33 projection_matrix_from_rising_;
  Matrice33 inverse_projection_matrix_from_rising_;

  double normal_interfacial_gradient_ = 0;
  Vecteur3 normal_interfacial_gradient_compo_;

  // FIXME: Should each probes have their own number of points ?
  bool global_probes_characteristics_ = true;

  const int * points_per_thermal_subproblem_base_;
  const int * points_per_thermal_subproblem_;
  int increased_point_numbers_ = 32;
  // FIXME: Should alpha_liq be constant, or a reference ?
  const double * alpha_;
  const double * lambda_;
  double Pr_l_ = 0.;
  const double * coeff_distance_diagonal_;
  const double * cell_diagonal_;
  double probe_length_ = 0.;
  double surface_ = 0.;

  double r_sph_ = 0.;
  double theta_sph_ = 0.;
  double phi_sph_ = 0.;
  Vecteur3 er_sph_;
  Vecteur3 etheta_sph_;
  Vecteur3 ephi_sph_;

  /*
   * References to IJK_Field_double to avoid copy of large fields
   * Similar to operators !
   */
  const IJK_Interfaces * interfaces_;
  double indicator_=0.5;
  const IJK_Field_double * eulerian_distance_;
  const IJK_Field_double * eulerian_curvature_;
  const IJK_Field_double * eulerian_interfacial_area_;
  const FixedVector<IJK_Field_double, 3> * eulerian_normal_vect_;
  const FixedVector<IJK_Field_double, 3> * eulerian_facets_barycentre_;

  const IJK_Field_double * temperature_;
  const IJK_Field_double * temperature_ft_;
  const IJK_Field_double * temperature_before_extrapolation_;
  const FixedVector<IJK_Field_double, 3> * velocity_;
  const FixedVector<IJK_Field_double, 3> * velocity_ft_;

  const FixedVector<IJK_Field_double, 3> * grad_T_elem_;
  const FixedVector<IJK_Field_double, 3> * hess_diag_T_elem_;
  const FixedVector<IJK_Field_double, 3> * hess_cross_T_elem_;

  const double * dr_base_ = 0;
  // FIXME: Should I use DoubleTab instead ?
  const DoubleVect* radial_coordinates_base_;

  const Matrice *identity_matrix_explicit_implicit_base_;
  const Matrice *radial_first_order_operator_raw_base_;
  const Matrice *radial_second_order_operator_raw_base_;
  const Matrice *radial_first_order_operator_base_;
  const Matrice *radial_second_order_operator_base_;
  const Matrice *radial_first_order_operator_;
  const Matrice *radial_second_order_operator_;
  const Matrice *identity_matrix_explicit_implicit_;
  Matrice identity_matrix_explicit_implicit_local_;
  Matrice radial_first_order_operator_local_;
  Matrice radial_second_order_operator_local_;
  /*
   * Pointers to non-constant matrice
   * FIXME: Should I declare constant pointers ?
   */
  Matrice *identity_matrix_subproblems_;
  Matrice *radial_diffusion_matrix_base_;
  Matrice *radial_convection_matrix_base_;
  const Matrice *radial_velocity_convection_matrix_base_;
//  const Matrice* tangential_velocity_convection_matrix_base_;
//  const Matrice* azymuthal_velocity_convection_matrix_base_;

  double dr_=0.;
  double dr_inv_=0.;
  const DoubleVect * radial_coordinates_;
  DoubleVect radial_coordinates_modified_;
  DoubleVect osculating_radial_coordinates_;
  DoubleVect osculating_radial_coordinates_inv_;
  DoubleTab radial_coordinates_cartesian_compo_;
  DoubleTab osculating_radial_coordinates_cartesian_compo_;
  DoubleTab coordinates_cartesian_compo_;

  DoubleVect x_velocity_;
  DoubleVect y_velocity_;
  DoubleVect z_velocity_;
  DoubleVect velocity_magnitude_;
  DoubleVect x_velocity_corrected_;
  DoubleVect y_velocity_corrected_;
  DoubleVect z_velocity_corrected_;
  DoubleVect radial_velocity_;
  DoubleVect radial_velocity_advected_frame_;
  DoubleVect radial_velocity_static_frame_;
  DoubleVect radial_velocity_corrected_;
  DoubleVect first_tangential_velocity_;
  DoubleVect first_tangential_velocity_advected_frame_;
  DoubleVect first_tangential_velocity_static_frame_;
  DoubleVect first_tangential_velocity_corrected_;
  DoubleVect second_tangential_velocity_;
  DoubleVect second_tangential_velocity_advected_frame_;
  DoubleVect second_tangential_velocity_static_frame_;
  DoubleVect second_tangential_velocity_corrected_;
  DoubleVect first_tangential_velocity_from_rising_dir_;
  DoubleVect first_tangential_velocity_from_rising_dir_advected_frame_;
  DoubleVect first_tangential_velocity_from_rising_dir_static_frame_;
  DoubleVect first_tangential_velocity_from_rising_dir_corrected_;
  DoubleVect azymuthal_velocity_;
  DoubleVect azymuthal_velocity_advected_frame_;
  DoubleVect azymuthal_velocity_static_frame_;
  DoubleVect azymuthal_velocity_corrected_;
  DoubleVect * first_tangential_velocity_solver_;
  DoubleVect * second_tangential_velocity_solver_;
  DoubleVect radial_convection_prefactor_;
  DoubleVect temperature_interp_;
  DoubleVect temperature_time_increment_;
  DoubleVect temperature_time_increment_from_eulerian_;
  DoubleVect material_derivative_advected_frame_;
  DoubleVect material_derivative_static_frame_;
  DoubleVect material_derivative_advected_frame_rising_;
  DoubleVect material_derivative_static_frame_rising_;
  DoubleVect material_derivative_velocity_advected_frame_;
  DoubleVect material_derivative_velocity_static_frame_;
  DoubleVect material_derivative_velocity_advected_frame_rising_;
  DoubleVect material_derivative_velocity_static_frame_rising_;
  DoubleVect convective_term_advected_frame_;
  DoubleVect convective_term_static_frame_;
  DoubleVect convective_term_advected_frame_rising_;
  DoubleVect convective_term_static_frame_rising_;
  FixedVector<DoubleVect, 3> grad_T_elem_interp_;
  FixedVector<DoubleVect, 3> hess_diag_T_elem_interp_;
  FixedVector<DoubleVect, 3> hess_cross_T_elem_interp_;
  FixedVector<DoubleVect, 3> hess_diag_T_elem_spherical_;
  FixedVector<DoubleVect, 3> hess_cross_T_elem_spherical_;
  FixedVector<DoubleVect, 3> hess_diag_T_elem_spherical_from_rising_;
  FixedVector<DoubleVect, 3> hess_cross_T_elem_spherical_from_rising_;
  DoubleVect temperature_diffusion_hessian_cartesian_trace_;
  DoubleVect temperature_diffusion_hessian_trace_;
  DoubleVect radial_temperature_diffusion_;
  DoubleVect tangential_temperature_diffusion_;

  int source_terms_type_=0;
  enum Source_terms { linear_diffusion, spherical_diffusion, spherical_diffusion_approx,
                      tangential_conv_2D, tangential_conv_3D,
                      tangential_conv_2D_tangential_diffusion_3D, tangential_conv_3D_tangentual_diffusion_3D
                    };
  DoubleVect normal_temperature_gradient_;
  DoubleVect tangential_temperature_gradient_first_;
  DoubleVect tangential_temperature_gradient_second_;
  DoubleVect tangential_temperature_gradient_first_from_rising_dir_;
  DoubleVect azymuthal_temperature_gradient_;
  DoubleVect * tangential_temperature_gradient_first_solver_;
  DoubleVect * tangential_temperature_gradient_second_solver_;
  DoubleVect tangential_hessian_contribution_;
  DoubleVect tangential_convection_source_terms_first_;
  DoubleVect tangential_convection_source_terms_second_;
  DoubleVect tangential_convection_source_terms_;
  DoubleVect tangential_diffusion_source_terms_;
  DoubleVect source_terms_;

  int correct_radial_velocity_;
  int correct_tangential_temperature_gradient_;
  int correct_tangential_temperature_hessian_;

  IJK_Finite_Difference_One_Dimensional_Matrix_Assembler * finite_difference_assembler_;
  Matrice * thermal_subproblems_matrix_assembly_;
  DoubleVect * thermal_subproblems_rhs_assembly_;
  DoubleVect * thermal_subproblems_temperature_solution_;
  DoubleVect * thermal_subproblems_temperature_solution_ini_;
  DoubleVect rhs_assembly_;
  double interfacial_boundary_condition_value_;
  double end_boundary_condition_value_;
  int start_index_;
  int end_index_;
  DoubleVect temperature_solution_;
  FixedVector<DoubleVect, 3> temperature_gradient_solution_;
  DoubleVect normal_temperature_gradient_solution_;
  DoubleVect normal_temperature_double_derivative_solution_;
  DoubleVect temperature_x_gradient_solution_;
  DoubleVect temperature_y_gradient_solution_;
  DoubleVect temperature_z_gradient_solution_;
  DoubleVect thermal_flux_;

  DoubleVect normal_velocity_normal_gradient_;
  DoubleVect first_tangential_velocity_normal_gradient_;
  DoubleVect second_tangential_velocity_normal_gradient_;
  DoubleVect azymuthal_velocity_normal_gradient_;
  DoubleVect first_tangential_velocity_gradient_from_rising_dir_;

  double delta_T_subcooled_overheated_ = -1.;

  int order_approx_temperature_ext_=1;
  int avoid_post_processing_all_terms_=0;

  REF(IJK_FT_double) ref_ijk_ft_;
  bool is_updated_ = false;

  /*
   * Some tries to do explicit temporal variations at the beginning
   */
  DoubleVect temperature_ini_temporal_schemes_;
  bool is_first_time_step_ = false;
  int * first_time_step_temporal_;
  int first_time_step_explicit_ = 1;
  double current_time_ = 0.;
  double global_dt_cfl_ = 0.;
  double global_dt_fo_ = 0.;
  double global_time_step_ = 0.;
  double local_dt_cfl_ = 0.;
  double local_dt_fo_ = 0.;
  double local_time_step_ = 0.;
  double local_fourier_ = 1.;
  double local_cfl_ = 1.;
  double min_delta_xyz_=0.;
  double max_u_;
  double local_time_step_round_ = 0.;
  double local_time_step_overall_ = 0.;
  double local_fourier_time_step_probe_length_ = 0.;
  double local_cfl_time_step_probe_length_ = 0.;
  double local_dt_cfl_min_delta_xyz_=0.;
  int nb_iter_explicit_ = 0;
  int max_u_cartesian_ = 1;

  /*
   * Some tries to make the probe length varies at the beginning of the simulation
   */
  int first_time_step_varying_probes_ = 0;
  int probe_variations_enabled_ = 0;
  int probe_variations_priority_ = 0;
  double cfl_probe_length_ = 0.;
  double fourier_probe_length_ = 0.;
  double max_cfl_fourier_probe_length_ = 0.;
  int velocities_calculation_counter_ = 0;
  int disable_interpolation_in_mixed_cells_ = 0;
  int short_probe_condition_ = 0;
  int temperature_probe_condition_ = 0;
  int max_u_radial_=0;
  double cell_centre_distance_ = 0;
  FixedVector<bool,6> pure_liquid_neighbours_;
  FixedVector<double,6> face_centres_distance_;
  FixedVector<FixedVector<double,4>,6> vertices_centres_distance_;
  int correct_fluxes_ = 0;
  double cell_temperature_ = 0.;

  int distance_cell_faces_from_lrs_ = 0;
  int pre_initialise_thermal_subproblems_list_ = 0;

  int correct_temperature_cell_neighbours_ = 0;
  int correct_neighbours_rank_ = 1;
  int neighbours_corrected_rank_ = 1;
  FixedVector<int,3> pure_neighbours_corrected_sign_;
  std::vector<std::vector<std::vector<bool>>> pure_neighbours_to_correct_;
  std::vector<std::vector<std::vector<double>>> pure_neighbours_corrected_distance_;
  int dxyz_increment_bool_=0;
};

#endif /* IJK_One_Dimensional_Subproblem_included */
